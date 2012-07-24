/*
 * tag_plot.h - a tk canvas item to display a tag track
 *              on the Template Display.
 *
 * Andrew Whitwham, June 2012
 *
 */

#include <stdio.h>
#include <tcl.h>
#include <X11/Xlib.h>

#include "active_tags.h"
#include "template_draw.h"
#include "tagdb.h"
#include "gap_range.h"
#include "tag_plot.h"

/*
 * On windows we need to use TkPutImage instead
 */
#ifdef _WIN32
#  include <tkInt.h>
#  define XPutImage(a,b,c,d,e,f,g,h,i,j) TkPutImage(NULL, 0, (a),(b),(c),(d),(e),(f),(g),(h),(i),(j))
#endif


/* item structures and prototypes */

typedef struct TagPlot {
    Tk_Item header; // mandatory entry
    
    GC gc;  	    	    // tag graphic context
    Pixmap pm;	    	    /* pixmap to draw on */
    image_t *image;         /* image drawing goes on to */
    int background;         /* background colour */
    XColor *bg;
    int width, height;      /* pixmap width and height */
    Tk_Anchor anchor;	    /* pixmap anchorpoint */
    double an_x, an_y;      /* postioning points */
    double px, py;          /* current mouse pointer x, y (in world units)*/
    char *tag;              /* tag name under pointer */

    gap_range_t *gr;        /* range info */
    rangec_t *ar;           /* anno range info, use for now */
    tline *tl;              /* line structure for drawing */
    int nr;                 /* number of annotations found */
    
    double contig_start;    /* start value of contig */
    double contig_length;   /* length of contig */
    double wx0, wx1;	    /* world coords */
    double wy0, wy1;        
    double y_start, y_end;  /* world edges in y */
    double yzoom;
    double yz;
    int min_sz;
    double ymax;
 } TagPlot;

/* mandatory prototypes */

static int		tagplot_create(Tcl_Interp *interp,
			    Tk_Canvas canvas, struct Tk_Item *itemPtr,
			    int objc, Tcl_Obj *CONST objv[]);
static int		tagplot_configure(Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr, int objc,
			    Tcl_Obj *CONST objv[], int flags);
static int		tagplot_coords(Tcl_Interp *interp, Tk_Canvas canvas,
    	    	    	    Tk_Item *itemPtr, int objc, Tcl_Obj *CONST objv[]);
static void		tagplot_delete(Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display);
static void		tagplot_display(Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display, Drawable dst,
			    int x, int y, int width, int height);
static double		tagplot_point(Tk_Canvas canvas, Tk_Item *itemPtr,
    	    	    	    double *coordPtr);
static int		tagplot_area(Tk_Canvas canvas, Tk_Item *itemPtr, 
    	    	    	    double *rectPtr);
static int		tagplot_postscript(Tcl_Interp *interp, 
    	    	    	    Tk_Canvas canvas, Tk_Item *itemPtr, int prepass);
static void		tagplot_scale(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double originX, double originY,
			    double scaleX, double scaleY);
static void		tagplot_translate(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double deltaX, double deltaY);
			    
/* config options, some are meant only to be returned with itemcget */
static Tk_ConfigSpec tag_config_specs[] = {
    {TK_CONFIG_ANCHOR, "-anchor", NULL, NULL,
	"nw", Tk_Offset(TagPlot, anchor), TK_CONFIG_DONT_SET_DEFAULT},
    {TK_CONFIG_CUSTOM, "-range", NULL, NULL,
    	NULL, Tk_Offset(TagPlot, gr), TK_CONFIG_NULL_OK, &range_option},
    {TK_CONFIG_DOUBLE, "-wx0", NULL, NULL, "0", Tk_Offset(TagPlot, wx0), 0}, 
    {TK_CONFIG_DOUBLE, "-wx1", NULL, NULL, "0", Tk_Offset(TagPlot, wx1), 0},
    {TK_CONFIG_DOUBLE, "-wy0", NULL, NULL, "0", Tk_Offset(TagPlot, wy0), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-wy1", NULL, NULL, "0", Tk_Offset(TagPlot, wy1), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-yzoom", "yZoom", "YZoom", "100.0", Tk_Offset(TagPlot, yzoom), 0, 0},
    {TK_CONFIG_DOUBLE, "-y_start", NULL, NULL, "0", Tk_Offset(TagPlot, y_start), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-y_end", NULL, NULL, "0", Tk_Offset(TagPlot, y_end), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-px", NULL, NULL, "0", Tk_Offset(TagPlot, px), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-py", NULL, NULL, "0", Tk_Offset(TagPlot, py), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_COLOR, "-background", NULL, NULL, "black", Tk_Offset(TagPlot, bg), TK_CONFIG_NULL_OK},
    {TK_CONFIG_STRING, "-tag", NULL, NULL, NULL, Tk_Offset(TagPlot, tag), TK_CONFIG_DONT_SET_DEFAULT},
    {TK_CONFIG_END, (char *) NULL, (char *) NULL, (char *) NULL, (char *) NULL, 0, 0}
};			    


/* Define the canvas item type */
Tk_ItemType tkTagItem = {
    "tag_plot",
    sizeof(TagPlot),
    tagplot_create,
    tag_config_specs,
    tagplot_configure,
    tagplot_coords,
    tagplot_delete,
    tagplot_display,
    TK_CONFIG_OBJS, /* alwaysRedraw: whether callbacks accept Tcl_Obj */
    tagplot_point,
    tagplot_area,
    tagplot_postscript,
    tagplot_scale,
    tagplot_translate,
    (Tk_ItemIndexProc *)     NULL,
    (Tk_ItemCursorProc *)    NULL,
    (Tk_ItemSelectionProc *) NULL,
    (Tk_ItemInsertProc *)    NULL,
    (Tk_ItemDCharsProc *)    NULL,
    (Tk_ItemType *)          NULL,
};



static void colour_name_to_RGB(Tcl_Interp *interp, Tk_Window tk_win, unsigned short *red,
    	    	unsigned short *green, unsigned short *blue, char *name) {
    XColor *colour = Tk_GetColor(interp, tk_win, name);
    
    if (colour) {
    	*red   = colour->red;
	*green = colour->green;
	*blue  = colour->blue;
    } else {
    	// default to red
	*red   = 255;
	*green = 0;
	*blue  = 0;
    }
    
    return;
}


/* set up the image prior to drawing */

static int initialise_tagplot_image(Tcl_Interp *interp, Tk_Window tk_win, TagPlot *tp) {
    Display *display = Tk_Display(tk_win);
    int i;
    
    if (NULL == (tp->image = initialise_image(display))) {
    	fprintf(stderr, "TagPlot - unable to initialise image.\n");
	return 0;
    }
    
    // set up the colours for the tags
    for (i = 0; i < tag_db_count; i++) {
    	unsigned short red   = 255;
	unsigned short green = 0;
	unsigned short blue  = 0;
	
	
	if (tag_db[i].bg_colour) {
    	    colour_name_to_RGB(interp, tk_win, &red, &green, &blue, tag_db[i].bg_colour);
	}
	
	add_colour(tp->image, red, green, blue);
    }
    
    tp->background = add_colour(tp->image, 0, 0, 0);

    return 1;
}


/* --------------------------------------------------------------------------
 * Canvas interface implementation
 */
static int tagplot_create(Tcl_Interp *interp,
			Tk_Canvas canvas,
			Tk_Item *itemPtr,
			int argc,
			Tcl_Obj *CONST *argv) {
    TagPlot *tp = (TagPlot *)itemPtr;
    int i;
    
    if (argc==1) {
	i = 1;
    } else {
	char *arg = Tcl_GetStringFromObj(argv[1], NULL);
	
	if (((argc>1) && (arg[0] == '-')
		&& (arg[1] >= 'a') && (arg[1] <= 'z'))) {
	    i = 1;
	} else {
	    i = 2;
	}
    }

    if (argc < i) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
		Tk_PathName(Tk_CanvasTkwin(canvas)), " create ",
		itemPtr->typePtr->name, " x y ?options?\"",
		(char *) NULL);
	return TCL_ERROR;
    }
    
    /* initialise item record */
    tp->gc = NULL;
    tp->background = 0;
    tp->bg = NULL;
    tp->pm = None;
    tp->anchor = TK_ANCHOR_NW;
    tp->an_x = 0.0;
    tp->an_y = 0.0;
    tp->wy0 = tp->y_start = 0;
    tp->wy1 = tp->y_end = Tk_Height(Tk_CanvasTkwin(canvas)); /* initial world height */
    tp->width = -1;
    tp->height = -1;
    tp->yz = 100;
    tp->tl = NULL;
    tp->image = NULL;
    tp->ymax = 0;
    tp->tag = NULL;
    
    if (initialise_tagplot_image(interp, Tk_CanvasTkwin(canvas), tp)) {
    	if ((tagplot_coords(interp, canvas, itemPtr, i, argv) == TCL_OK)) {
    	    if (tagplot_configure(interp, canvas, itemPtr, argc - i, argv + i, 0) == TCL_OK) {
	    	return TCL_OK;
	    }
    	}
    }
    
    /* if we get here then something is wrong */
    fprintf(stderr, "TagPlot item creation failed\n");
    tagplot_delete(canvas, itemPtr, Tk_Display(Tk_CanvasTkwin(canvas)));
    
    return TCL_ERROR;
}

/* ---------------------------------------------------------------------------
 * Utility functions used by the main canvas below interface
 */


/* compute the bounding box of the item and return in the Tk_Item header */
static void update_bbox(Tk_Canvas canvas, TagPlot *tp) {
    Tk_Window tkwin;
    int width, height;
    int x, y;
    
    tkwin = Tk_CanvasTkwin(canvas);

    width  = Tk_Width(tkwin);
    height = Tk_Height(tkwin);
    
    x = (int) (tp->an_x + ((tp->an_x >= 0) ? 0.5 : - 0.5));
    y = (int) (tp->an_y + ((tp->an_y >= 0) ? 0.5 : - 0.5));
    
    switch (tp->anchor) {
	case TK_ANCHOR_N:
	    x -= width / 2;
	    break;
	case TK_ANCHOR_NE:
	    x -= width;
	    break;
	case TK_ANCHOR_E:
	    x -= width;
	    y -= height / 2;
	    break;
	case TK_ANCHOR_SE:
	    x -= width;
	    y -= height;
	    break;
	case TK_ANCHOR_S:
	    x -= width / 2;
	    y -= height;
	    break;
	case TK_ANCHOR_SW:
	    y -= height;
	    break;
	case TK_ANCHOR_W:
	    y -= height / 2;
	    break;
	case TK_ANCHOR_NW:
	    break;
    	case TK_ANCHOR_CENTER:
	    x -= width / 2;
	    y -= height / 2;
	    break;
    }
    
    tp->header.x1 = x;
    tp->header.y1 = y;
    tp->header.x2 = x + width;
    tp->header.y2 = y + height;
}


/* -------------------------------------------------------------------------
 * Y Coordinate allocation routines
 */

SPLAY_HEAD(xtag_TREE, xy_pair);
SPLAY_HEAD(ytag_TREE, xy_pair);
SPLAY_PROTOTYPE(xtag_TREE, xy_pair, x_link, x_cmp);
SPLAY_PROTOTYPE(ytag_TREE, xy_pair, y_link, y_cmp);
SPLAY_GENERATE(xtag_TREE, xy_pair, x_link, x_cmp);
SPLAY_GENERATE(ytag_TREE, xy_pair, y_link, y_cmp);

/*
 * This algorithm allocates Y coordinates to the horizontal lines listed
 * in tl[0..ntl-1].
 *
 * To keep track of this we expect sorted data, by X. We then keep track
 * of every display line as an entry in one of two splay trees; one
 * sorted on X and one sorted on Y. (NB, change macros SPLAY_ to RB_ for
 * a red-black tree implementation instead.)
 *
 * The X-sorted tree holds the used portion of active display rows.
 * The Y-sorted tree holds rows that can be considered as inactive, but
 * will be reused if we have to add another row.
 */
 
static int compute_ypos(TagPlot  *tp, int xgap, tline *tl, int ntl) {
    int i;
    struct xy_pair *node, *curr, *next;
    int yn = 0;
    int min_sz = INT_MIN;
    int max_sz = tp->min_sz;
    int yoffset = 0, ymax = 0, nleft;

    /* Create and initialise X and Y trees */
    struct xtag_TREE xtree = SPLAY_INITIALIZER(&xtree);
    struct ytag_TREE ytree = SPLAY_INITIALIZER(&ytree);
    
    tp->ymax = 0;

    nleft = ntl;
    for (i = 0; i < ntl; i++) {
	tl[i].y = -1;
    }

    while (nleft) {
	ymax = 0;
	yn = 0;

	/* Compute Y coords */
	for (i = 0; i < ntl; i++) {
	    if (tl[i].y != -1)
		continue;

	    if (tl[i].x[3] - tl[i].x[0] + 1 < min_sz ||
		tl[i].x[3] - tl[i].x[0] + 1 > max_sz)
		continue;
	    nleft--;

	    if ((node = SPLAY_MIN(xtag_TREE, &xtree)) != NULL && tl[i].x[0] >= node->x) {
		int try_cull = 0;

		/* We found a node, is it the smallest in y? */
		curr = SPLAY_NEXT(xtag_TREE, &xtree, node);
		while (curr && tl[i].x[0] >= curr->x) {
		    if (node->y > curr->y)
			node = curr;
		    curr = SPLAY_NEXT(xtag_TREE, &xtree, curr);
		}
	    
		/* Shift non-smallest y (but < x) to Y-tree */
		curr = SPLAY_MIN(xtag_TREE, &xtree);
		while (curr && tl[i].x[0] >= curr->x) {
		    next = SPLAY_NEXT(xtag_TREE, &xtree, curr);
		    if (curr != node) {
			SPLAY_REMOVE(xtag_TREE, &xtree, curr);
			SPLAY_INSERT(ytag_TREE, &ytree, curr);
			try_cull = 1;
		    }
		    curr = next;
		}
	    
		tl[i].y = node->y + yoffset;
		SPLAY_REMOVE(xtag_TREE, &xtree, node);
		node->x = tl[i].x[3] + xgap;
		SPLAY_INSERT(xtag_TREE, &xtree, node);
	    } else {
		/* Check if we have a free y on ytree */
		if ((node = SPLAY_MIN(ytag_TREE, &ytree)) != NULL) {
		    SPLAY_REMOVE(ytag_TREE, &ytree, node);
		} else {
		    node = (struct xy_pair *)malloc(sizeof(*node));
		    node->y = ++yn;
		    if (ymax < yn)
			ymax = yn;
		}
		tl[i].y = node->y + yoffset;
		node->x = tl[i].x[3] + xgap;
		SPLAY_INSERT(xtag_TREE, &xtree, node);
	    }
	}

	/* Delete trees */
	for (node = SPLAY_MIN(xtag_TREE, &xtree); node; node = next) {
	    next = SPLAY_NEXT(xtag_TREE, &xtree, node);
	    SPLAY_REMOVE(xtag_TREE, &xtree, node);
	}

	for (node = SPLAY_MIN(ytag_TREE, &ytree); node; node = next) {
	    next = SPLAY_NEXT(ytag_TREE, &ytree, node);
	    SPLAY_REMOVE(ytag_TREE, &ytree, node);
	}

	min_sz = max_sz;
	max_sz += tp->min_sz;
	
	if (max_sz < min_sz+10)
	    max_sz = min_sz+10;
	yoffset += ymax;
	
    	if (tp->ymax < ymax) ymax = tp->ymax;
    }
    
   return 0;
}




static int sort_tline_by_x(const void *p1, const void *p2) {
    tline *r1 = (tline *)p1;
    tline *r2 = (tline *)p2;

    return r1->x[0] - r2->x[3];
}


static void tagplot_redraw(TagPlot *tp, Display *display) {
    int i;
    double xgap = 100;
    double ax, bx, ay, by;
    int ymin = INT_MAX;
    int ymax = INT_MIN;
    contig_t *c;
     
    image_remove(tp->image);
    if(!create_image_buffer(tp->image, tp->width, tp->height, tp->background)) return;

    // probably need to put something to not query the range code each time
    c = cache_search(tp->gr->io, GT_Contig, tp->gr->crec);
    cache_incr(tp->gr->io, c);
    tp->ar = contig_anno_in_range(tp->gr->io, &c, tp->wx0, tp->wx1, CSIR_SORT_BY_X, &tp->nr);
    cache_decr(tp->gr->io, c);
    
    if (tp->nr == 0) { // if no tags have a blank track
    	create_image_from_buffer(tp->image);
    	XPutImage(display, (Drawable)tp->pm, tp->gc, tp->image->img, 0, 0, 0, 0, tp->width, tp->height);
	return;
    }
    
    /* world to pixmap conversion values */
    if (tp->wx1 - tp->wx0 == 0) return;
    
    ax = tp->width / (tp->wx1 - tp->wx0);
    bx = tp->wx0;
    
    ay = tp->height / (tp->wy1 - tp->wy0);
    by = tp->wy0;
        

    tp->tl = (tline *)realloc(tp->tl, tp->nr * sizeof(tline));
    
    for (i = 0; i < tp->nr; i++) {
    	char type[5];
    
	tp->tl[i].x[0] = tp->ar[i].start;
	tp->tl[i].x[3] = tp->ar[i].end;
	tp->tl[i].col[0] = idToIndex(type2str(tp->ar[i].mqual, type));
    }
    
    // now make some Y positions
    qsort(tp->tl, tp->nr, sizeof(tline), sort_tline_by_x);
    compute_ypos(tp, xgap, tp->tl, tp->nr);
    

    // plot the lines
    if (tp->ymax) {
    	tp->yz = ((double)tp->height / (double)tp->ymax) * 0.95;
    } else {
    	tp->yz = tp->yzoom / 50.0; // this will do for now
    }
    
    for (i = 0; i < tp->nr; i++) {
    	double y;
	
	y = (tp->tl[i].y + 5) * tp->yz;
	
	if (ymin > y) ymin = y;
	if (ymax < y) ymax = y;
	
	// plot if visible
	if (y >= tp->wy0 && y <= tp->wy1) {
	    int rx1, rx2, ry;
	    
	    rx1 = (tp->tl[i].x[0] - bx) * ax; // world to raster conversion
	    rx2 = ((tp->tl[i].x[3] + 1) - bx) * ax;
	    ry  = (y - by) * ay;
	    
	    
	    if ((rx2 - rx1) <= 1) { // 1 pixel is too small to see easily, make it 2
	    	rx2++;
	    }
	    
	    draw_line(tp->image, rx1, rx2, ry, tp->tl[i].col[0]);
	    draw_line(tp->image, rx1, rx2, ry + 1, tp->tl[i].col[0]);
	}
    }
    
    create_image_from_buffer(tp->image);
    XPutImage(display, (Drawable)tp->pm, tp->gc, tp->image->img, 0, 0, 0, 0, tp->width, tp->height);
    
    /* some last bits of size calculation for scrolling */
    tp->y_start = ymin - 10;
    tp->y_end   = ymax + 10;
    
    if (tp->y_start > 0) tp->y_start = 0;
    if (tp->y_end < tp->height) tp->y_end = tp->height;
}    


static int tagplot_configure(Tcl_Interp *interp,
	       Tk_Canvas canvas,
	       Tk_Item *itemPtr,
	       int argc,
	       Tcl_Obj *CONST *argv,
	       int flags) {
    TagPlot *tp = (TagPlot *) itemPtr;
    XGCValues gcValues;
    Tk_Window tkwin;
    Display *display;
    unsigned long mask = 0;
    int width, height;
        
    tkwin   = Tk_CanvasTkwin(canvas);
    display = Tk_Display(tkwin);
    
    if (Tk_ConfigureWidget(interp, tkwin, tag_config_specs, argc, (char **) argv,
	    (char *) tp, flags|TK_CONFIG_OBJS) != TCL_OK) {

	verror(ERR_WARN, "tagplot_configure", "%s", Tcl_GetStringResult(interp));
	return TCL_ERROR;
    }
    
    width  = Tk_Width(tkwin);
    height = Tk_Height(tkwin);
    
    if (tp->gc == None) {
    	mask = 0;
	tp->gc = Tk_GetGC(tkwin, mask, &gcValues);
    }
	
    if (width != tp->width || height != tp->height) {
    	/* going to need a new pixmap */
 	if (tp->pm) {
	    Tk_FreePixmap(display, tp->pm);
	}

    	tp->width = width;
    	tp->height = height;	
	
	tp->pm = Tk_GetPixmap(display, Tk_WindowId(tkwin), tp->width, tp->height, DefaultDepthOfScreen(Tk_Screen(tkwin)));
    }
    
    tagplot_redraw(tp, display);
    update_bbox(canvas, tp);
    return TCL_OK;
}


static void tagplot_delete(Tk_Canvas canvas, Tk_Item *itemPtr, Display *display) {
    TagPlot *tp = (TagPlot *) itemPtr;
    
    if (tp->pm != None) {
    	Tk_FreePixmap(display, tp->pm);
    }
    
    if (tp->gc != None) {
    	Tk_FreeGC(display, tp->gc);
    }
    
    if (tp->image != NULL) {
    	image_destroy(tp->image);
    }
    
    if (tp->tl != NULL) {
    	free(tp->tl);
    }
}

    
/* invoked by the coords command */

static int tagplot_coords(Tcl_Interp *interp, Tk_Canvas canvas,		
    	    	    	    Tk_Item *itemPtr, int objc, Tcl_Obj *CONST objv[]) {
    TagPlot *tp = (TagPlot *)itemPtr;
    
    if (objc == 0) {
	Tcl_Obj *obj = Tcl_NewObj();

	Tcl_Obj *subobj = Tcl_NewDoubleObj(tp->an_x);
    	Tcl_ListObjAppendElement(interp, obj, subobj);

	subobj = Tcl_NewDoubleObj(tp->an_y);
	Tcl_ListObjAppendElement(interp, obj, subobj);

	Tcl_SetObjResult(interp, obj);
    } else if (objc < 3) {
	if (objc == 1) {
	    if (Tcl_ListObjGetElements(interp, objv[0], &objc,
		    (Tcl_Obj ***) &objv) != TCL_OK) {
		return TCL_ERROR;
	    } else if (objc != 2) {
		char buf[64 + TCL_INTEGER_SPACE];

		sprintf(buf, "wrong # coordinates: expected 2, got %d", objc);
		Tcl_SetResult(interp, buf, TCL_VOLATILE);
		return TCL_ERROR;
	    }
	}
	
	if ((Tk_CanvasGetCoordFromObj(interp, canvas, objv[0],
		&tp->an_x) != TCL_OK)
		|| (Tk_CanvasGetCoordFromObj(interp, canvas, objv[1],
			&tp->an_y) != TCL_OK)) {
	    return TCL_ERROR;
	}
	
	update_bbox(canvas, tp);
    } else {
	char buf[64 + TCL_INTEGER_SPACE];

	sprintf(buf, "wrong # coordinates: expected 0 or 2, got %d", objc);
	Tcl_SetResult(interp, buf, TCL_VOLATILE);
	return TCL_ERROR;
    }
    
    return TCL_OK;
}


static char *tag_under_pointer(TagPlot *tp, double cpy) {
    // lets try with a brute force and ignorance approach
    
    int i;
    double cx;
    double cy, ay, by;
    
    cx = (tp->wx1 - tp->wx0) / tp->width;

    if (cx < 1) cx = 1;
    
    if ((tp->wy1 - tp->wy0) == 0) return NULL;   
    
    ay = tp->height / (tp->wy1 - tp->wy0);
    by = tp->wy0;

    for (i = 0; i < tp->nr; i++) {
    	cy = (tp->tl[i].y + 5) * tp->yz;
	cy  = (cy - by) * ay;
	
    	if ((tp->tl[i].x[0] <= tp->px) && (tp->tl[i].x[3] + cx >= tp->px)) {
	    if (cy == cpy || (cy + 1) == cpy) {
	    	return tag_db[tp->tl[i].col[0]].search_id;
	    }
	}
    }
    
    return NULL;
}
    
    
/* subvert the point function so that it now
   now provides world coords for xhair feedback in depth.tcl
*/

static double tagplot_point(Tk_Canvas canvas, Tk_Item *itemPtr, double *coordPtr) {
    TagPlot *tp = (TagPlot *)itemPtr;
    double ax, bx;

    tp->px = coordPtr[0];
    tp->py = coordPtr[1];

    tp->py = (tp->height - tp->py) / tp->yz;
    
    ax = tp->width / (tp->wx1 - tp->wx0);
    bx = tp->wx0;
    
    tp->px = (tp->px / ax) + bx;
    
    // see if there is a tag under the pointer
    tp->tag = tag_under_pointer(tp, coordPtr[1]);
    

    return 0;
}


/* determines whether an item lies entirely
   inside, entirely outside, or overlapping a given rectangle
*/

static int tagplot_area(Tk_Canvas canvas, Tk_Item *itemPtr, double *rectPtr) {
    TagPlot *tp = (TagPlot *)itemPtr;

    if ((rectPtr[2] <= tp->header.x1)
	    || (rectPtr[0] >= tp->header.x2)
	    || (rectPtr[3] <= tp->header.y1)
	    || (rectPtr[1] >= tp->header.y2)) {
	return -1;
    }
    if ((rectPtr[0] <= tp->header.x1)
	    && (rectPtr[1] <= tp->header.y1)
	    && (rectPtr[2] >= tp->header.x2)
	    && (rectPtr[3] >= tp->header.y2)) {
	return 1;
    }
    
    return 0;
}


/* unimplemented scaling function */

static void tagplot_scale(Tk_Canvas canvas, Tk_Item *itemPtr,
    	    	    	    double originX, double originY,
			    double scaleX, double scaleY) {
    fprintf(stderr, "TagPlot scale - not implemented.\n");
}


/* move item by given amount */

static void tagplot_translate(Tk_Canvas canvas, Tk_Item *itemPtr, 
    	    	    	    	double deltaX, double deltaY) {

    TagPlot *tp = (TagPlot *)itemPtr;

    tp->an_x += deltaX;
    tp->an_y += deltaY;
    
    update_bbox(canvas, tp);
}


/* unimplemented function to generate postscript data for printing */

static int tagplot_postscript(Tcl_Interp *interp, Tk_Canvas canvas, Tk_Item *itemPtr,	
    	    	    	    	int prepass) {
    fprintf(stderr, "TagPlot postscript - not implemented\n");
    
    return TCL_ERROR;
}


static void tagplot_display(Tk_Canvas canvas, Tk_Item *itemPtr,
    	    	    	    Display *display, Drawable drawable,
			    int x, int y, int width, int height) {
    
    TagPlot *tp = (TagPlot *) itemPtr;
    int pix_x, pix_y, pix_w, pix_h;
    short draw_x, draw_y;

    if (tp->pm != None) {
	if (x > tp->header.x1) {
	    pix_x = x - tp->header.x1;
	    pix_w = tp->header.x2 - x;
	} else {
	    pix_x = 0;
	    if ((x + width) < tp->header.x2) {
		pix_w = x + width - tp->header.x1;
	    } else {
		pix_w = tp->header.x2 - tp->header.x1;
	    }
	}
	if (y > tp->header.y1) {
	    pix_y = y - tp->header.y1;
	    pix_h = tp->header.y2 - y;
	} else {
	    pix_y = 0;
	    if ((y + height) < tp->header.y2) {
		pix_h = y + height - tp->header.y1;
	    } else {
		pix_h = tp->header.y2 - tp->header.y1;
	    }
	}
	
	Tk_CanvasDrawableCoords(canvas,
		(double) (tp->header.x1 + pix_x),
		(double) (tp->header.y1 + pix_y),
		&draw_x, &draw_y);

	/*
	 * Must modify the mask origin within the graphics context to line up
	 * with the bitmap's origin (in order to make bitmaps with
	 * "-background {}" work right).
	 */

 	XSetClipOrigin(display, tp->gc, draw_x - pix_x, draw_y - pix_y);
	XCopyArea(display, tp->pm, drawable, tp->gc, pix_x, pix_y,(unsigned int) pix_w,
		(unsigned int) pix_h, draw_x, draw_y); 
	
	XSetClipOrigin(display, tp->gc, 0, 0);
    }
}
    
    
    
 
        
