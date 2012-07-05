/*
 *
 * template_display.c - a tk canvas item to show the template display
 *                      formerly used tkRaster.
 *
 * Andrew Whitwham, July 2010
 *
 */

#include <stdio.h>
#include <math.h>
#include <tcl.h>
#include <X11/Xlib.h>

#include "template_draw.h"
#include "gap_range.h"
#include "template_display.h"

/*
 * On windows we need to use TkPutImage instead
 */
#ifdef _WIN32
#  include <tkInt.h>
#  define XPutImage(a,b,c,d,e,f,g,h,i,j) TkPutImage(NULL, 0, (a),(b),(c),(d),(e),(f),(g),(h),(i),(j))
#endif

/* Define the template display item */

typedef struct TemplateDisplayItem {
    Tk_Item header; 	    /* mandatory entry */
    GC gc;  	    	    /* graphics context */
    Pixmap pm;	    	    /* pixmap to draw on */
    Tk_Anchor anchor;	    /* pixmap anchorpoint */
    double an_x, an_y;      /* postioning points */
    double px, py;          /* current mouse pointer x, y (in world units)*/
    
    gap_range_t *gr;        /* range info */
    double contig_start;    /* start value of contig */
    double contig_length;   /* length of contig */
    double wx0, wx1;	    /* world coords */
    double wy0, wy1;        
    double y_start, y_end;  /* world edges in y */
    
    int logy;
    int cmode;
    int ymode;
    int yoffset;
    int accuracy;
    int spread;
    int reads_only;
    int sep_by_strand;
    int filter;
    int min_qual;
    int max_qual;
    int min_sz;
    double yzoom;
    double yz;
    
    int width, height;      /* pixmap width and height */
    image_t *image;         /* image drawing goes on to */
    int single_col;         /* special colours for image below */
    int span_col;
    int inconsistent_col;
    int fwd_col;
    int rev_col;
    int fwd_col3;
    int rev_col3;
    int background;
    
    int ntl;
             
} TemplateDisplayItem;



/* mandatory prototypes */

static int		template_coords(Tcl_Interp *interp, Tk_Canvas canvas,
    	    	    	    Tk_Item *itemPtr, int argc, Tcl_Obj *CONST argv[]);
static int		template_area(Tk_Canvas canvas, Tk_Item *itemPtr, 
    	    	    	    double *rectPtr);
static double		template_point(Tk_Canvas canvas, Tk_Item *itemPtr,
    	    	    	    double *coordPtr);
static int		template_postscript(Tcl_Interp *interp, 
    	    	    	    Tk_Canvas canvas, Tk_Item *itemPtr, int prepass);
static int		configure_template(Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr, int argc,
			    Tcl_Obj *CONST argv[], int flags);
static int		create_template(Tcl_Interp *interp,
			    Tk_Canvas canvas, struct Tk_Item *itemPtr,
			    int argc, Tcl_Obj *CONST argv[]);
static void		delete_template(Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display);
static void		display_template(Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display, Drawable dst,
			    int x, int y, int width, int height);
static void		scale_template(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double originX, double originY,
			    double scaleX, double scaleY);
static void		translate_template(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double deltaX, double deltaY);
			    
/* none mandatory protoypes */
static void compute_template_bbox(Tk_Canvas canvas, TemplateDisplayItem *tdi);
static int  initialise_template_image(TemplateDisplayItem *tdi, Display *display);
static void redraw_template_image(TemplateDisplayItem *tdi, Display *display);



/* config options, some are meant only to be returned with itemcget */
static Tk_ConfigSpec config_specs[] = {
    {TK_CONFIG_ANCHOR, "-anchor", NULL, NULL,
	"nw", Tk_Offset(TemplateDisplayItem, anchor), TK_CONFIG_DONT_SET_DEFAULT},
    {TK_CONFIG_CUSTOM, "-range", NULL, NULL,
    	NULL, Tk_Offset(TemplateDisplayItem, gr), TK_CONFIG_NULL_OK, &range_option},
    {TK_CONFIG_DOUBLE, "-contig_start", NULL, NULL, "0", Tk_Offset(TemplateDisplayItem, contig_start), 0}, 
    {TK_CONFIG_DOUBLE, "-contig_length", NULL, NULL, "0", Tk_Offset(TemplateDisplayItem, contig_length), 0}, 
    {TK_CONFIG_INT, "-logy", "logy", "LogY", "0", Tk_Offset(TemplateDisplayItem, logy), 0, 0},
    {TK_CONFIG_INT, "-cmode", "colourMode", "ColourMode", "0", Tk_Offset(TemplateDisplayItem, cmode), 0, 0},
    {TK_CONFIG_INT, "-ymode", "YMode", "YMode", "0", Tk_Offset(TemplateDisplayItem, ymode), 0, 0},
    {TK_CONFIG_INT, "-yoffset", "YOffset", "YOffset", "0", Tk_Offset(TemplateDisplayItem, yoffset), 0, 0},
    {TK_CONFIG_INT, "-accuracy", "accuracy", "Accuracy", "0", Tk_Offset(TemplateDisplayItem, accuracy), 0, 0},
    {TK_CONFIG_INT, "-spread", "spread", "Spread", "0", Tk_Offset(TemplateDisplayItem, spread), 0, 0},
    {TK_CONFIG_INT, "-reads_only", "readsOnly", "ReadsOnly", "0", Tk_Offset(TemplateDisplayItem, reads_only), 0, 0},
    {TK_CONFIG_INT, "-by_strand", "byStrand", "ByStrand", "1", Tk_Offset(TemplateDisplayItem, sep_by_strand), 0, 0},
    {TK_CONFIG_INT, "-filter", "filter", "Filter", "0", Tk_Offset(TemplateDisplayItem, filter), 0, 0},
    {TK_CONFIG_DOUBLE, "-yzoom", "yZoom", "YZoom", "10.0", Tk_Offset(TemplateDisplayItem, yzoom), 0, 0},
    {TK_CONFIG_DOUBLE, "-yz", "yZ", "YZ", "0", Tk_Offset(TemplateDisplayItem, yz), 0, 0},
    {TK_CONFIG_INT, "-min_qual", "minQual", "MinQual", "0", Tk_Offset(TemplateDisplayItem, min_qual), 0, 0},
    {TK_CONFIG_INT, "-max_qual", "maxQual", "MaxQual", "255", Tk_Offset(TemplateDisplayItem, max_qual), 0, 0},
    {TK_CONFIG_INT, "-min_y_size", "minYSize", "MinYSize", "512", Tk_Offset(TemplateDisplayItem, min_sz), 0, 0},
    {TK_CONFIG_DOUBLE, "-wx0", NULL, NULL, "0", Tk_Offset(TemplateDisplayItem, wx0), 0}, 
    {TK_CONFIG_DOUBLE, "-wx1", NULL, NULL, "0", Tk_Offset(TemplateDisplayItem, wx1), 0}, 
    {TK_CONFIG_DOUBLE, "-wy0", NULL, NULL, "0", Tk_Offset(TemplateDisplayItem, wy0), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-wy1", NULL, NULL, "0", Tk_Offset(TemplateDisplayItem, wy1), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-y_start", NULL, NULL, "0", Tk_Offset(TemplateDisplayItem, y_start), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-y_end", NULL, NULL, "0", Tk_Offset(TemplateDisplayItem, y_end), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-px", NULL, NULL, "0", Tk_Offset(TemplateDisplayItem, px), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-py", NULL, NULL, "0", Tk_Offset(TemplateDisplayItem, py), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_END, (char *) NULL, (char *) NULL, (char *) NULL, (char *) NULL, 0, 0}
};
    
/* define the TemplateDisplayItem functions */
Tk_ItemType tkTDItem = {
    "template_display",
    sizeof(TemplateDisplayItem),
    create_template,
    config_specs,
    configure_template,
    template_coords,
    delete_template,
    display_template,
    TK_CONFIG_OBJS,          /* needs to be this */
    template_point,
    template_area,
    template_postscript,
    scale_template,
    translate_template,
    (Tk_ItemIndexProc *)     NULL,
    (Tk_ItemCursorProc *)    NULL,
    (Tk_ItemSelectionProc *) NULL,
    (Tk_ItemInsertProc *)    NULL,
    (Tk_ItemDCharsProc *)    NULL,
    (Tk_ItemType *)          NULL,
};


/* this function is called when item is created as part
   of a canvas */

static int create_template(Tcl_Interp *interp,
	    Tk_Canvas canvas,
	    Tk_Item *itemPtr,
	    int argc,
	    Tcl_Obj *CONST *argv) {
    TemplateDisplayItem *tdi = (TemplateDisplayItem *)itemPtr;
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
    
    tdi->gc = None;
    tdi->pm = None;
    tdi->anchor = TK_ANCHOR_CENTER;
    tdi->an_x = 0.0;
    tdi->an_y = 0.0;
    tdi->image = NULL;
    tdi->gr = NULL;
    tdi->wy0 = tdi->y_start = 0;
    tdi->wy1 = tdi->y_end = Tk_Height(Tk_CanvasTkwin(canvas)); /* initial world height */
    tdi->width = -1;
    tdi->height = -1;
   
    if(initialise_template_image(tdi, Tk_Display(Tk_CanvasTkwin(canvas)))) {
	if ((template_coords(interp, canvas, itemPtr, i, argv) == TCL_OK)) {
    	    if (configure_template(interp, canvas, itemPtr, argc - i, argv + i, 0) == TCL_OK) {
		// possibly more initialisation here

	        return TCL_OK;
    	    }
	}
    }
    
    /* if we get here then something is wrong */
    verror(ERR_WARN, "create_template", 
	   "TemplateDisplayItem creation failed\n");
    delete_template(canvas, itemPtr, Tk_Display(Tk_CanvasTkwin(canvas)));
    
    return TCL_ERROR;
}


/* invoked by the coords command */
   	
static int template_coords(Tcl_Interp *interp, Tk_Canvas canvas,		
    	    	    	    Tk_Item *itemPtr, int objc, Tcl_Obj *CONST objv[]) {
    TemplateDisplayItem *tdi = (TemplateDisplayItem *) itemPtr;
    
    if (objc == 0) {
	Tcl_Obj *obj = Tcl_NewObj();

	Tcl_Obj *subobj = Tcl_NewDoubleObj(tdi->an_x);
    	Tcl_ListObjAppendElement(interp, obj, subobj);

	subobj = Tcl_NewDoubleObj(tdi->an_y);
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
		&tdi->an_x) != TCL_OK)
		|| (Tk_CanvasGetCoordFromObj(interp, canvas, objv[1],
			&tdi->an_y) != TCL_OK)) {
	    return TCL_ERROR;
	}
	
	compute_template_bbox(canvas, tdi);
    } else {
	char buf[64 + TCL_INTEGER_SPACE];

	sprintf(buf, "wrong # coordinates: expected 0 or 2, got %d", objc);
	Tcl_SetResult(interp, buf, TCL_VOLATILE);
	return TCL_ERROR;
    }
    
    return TCL_OK;
}


/* configure and redraw the template track */

static int configure_template(Tcl_Interp *interp,
	       Tk_Canvas canvas,
	       Tk_Item *itemPtr,
	       int argc,
	       Tcl_Obj *CONST *argv,
	       int flags) {
    TemplateDisplayItem *tdi = (TemplateDisplayItem *) itemPtr;
    XGCValues gcValues;
    Tk_Window tkwin;
    Display *display;
    unsigned long mask = 0;
    int width, height;
        
    tkwin   = Tk_CanvasTkwin(canvas);
    display = Tk_Display(tkwin);
    
    if (Tk_ConfigureWidget(interp, tkwin, config_specs, argc, (char **) argv,
	    (char *) tdi, flags|TK_CONFIG_OBJS) != TCL_OK) {

	verror(ERR_WARN, "configure_template", "%s",
	       Tcl_GetStringResult(interp));
	return TCL_ERROR;
    }
    
    width  = Tk_Width(tkwin);
    height = Tk_Height(tkwin);
    
    if (tdi->gc == None) {
    	mask = 0;
	tdi->gc = Tk_GetGC(tkwin, mask, &gcValues);
    }
	
    if (width != tdi->width || height != tdi->height) {
    	/* going to need a new pixmap */
 	if (tdi->pm) {
	    Tk_FreePixmap(display, tdi->pm);
	}

    	tdi->width = width;
    	tdi->height = height;	
	
	tdi->pm = Tk_GetPixmap(display, Tk_WindowId(tkwin), tdi->width, tdi->height, DefaultDepthOfScreen(Tk_Screen(tkwin)));
    }
    
    redraw_template_image(tdi, display);
    compute_template_bbox(canvas, tdi);
    return TCL_OK;
}


/* free resources */

static void delete_template(Tk_Canvas canvas, Tk_Item *itemPtr, Display *display) {
    TemplateDisplayItem *tdi = (TemplateDisplayItem *) itemPtr;
    
    if (tdi->pm != None) {
    	Tk_FreePixmap(display, tdi->pm);
    }
    
    if (tdi->gc != None) {
    	Tk_FreeGC(display, tdi->gc);
    }
    
    if (tdi->image != NULL) {
    	image_destroy(tdi->image);
    }
}


/* compute the bounding box of the item */

static void compute_template_bbox(Tk_Canvas canvas, TemplateDisplayItem *tdi) {
    Tk_Window tkwin;
    int width, height;
    int x, y;
    
    tkwin = Tk_CanvasTkwin(canvas);

    width  = Tk_Width(tkwin);
    height = Tk_Height(tkwin);
    
    x = (int) (tdi->an_x + ((tdi->an_x >= 0) ? 0.5 : - 0.5));
    y = (int) (tdi->an_y + ((tdi->an_y >= 0) ? 0.5 : - 0.5));
    
    switch (tdi->anchor) {
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
    
    tdi->header.x1 = x;
    tdi->header.y1 = y;
    tdi->header.x2 = x + width;
    tdi->header.y2 = y + height;
}
    

/* draw the template track onto the screen, copies from the pixmap to
   the canvas */			    
			    
static void display_template(Tk_Canvas canvas, Tk_Item *itemPtr,
    	    	    	    Display *display, Drawable drawable,
			    int x, int y, int width, int height) {
    
    TemplateDisplayItem *tdi = (TemplateDisplayItem *) itemPtr;
    int pix_x, pix_y, pix_w, pix_h;
    short draw_x, draw_y;

    if (tdi->pm != None) {
	if (x > tdi->header.x1) {
	    pix_x = x - tdi->header.x1;
	    pix_w = tdi->header.x2 - x;
	} else {
	    pix_x = 0;
	    if ((x + width) < tdi->header.x2) {
		pix_w = x + width - tdi->header.x1;
	    } else {
		pix_w = tdi->header.x2 - tdi->header.x1;
	    }
	}
	if (y > tdi->header.y1) {
	    pix_y = y - tdi->header.y1;
	    pix_h = tdi->header.y2 - y;
	} else {
	    pix_y = 0;
	    if ((y + height) < tdi->header.y2) {
		pix_h = y + height - tdi->header.y1;
	    } else {
		pix_h = tdi->header.y2 - tdi->header.y1;
	    }
	}
	
	Tk_CanvasDrawableCoords(canvas,
		(double) (tdi->header.x1 + pix_x),
		(double) (tdi->header.y1 + pix_y),
		&draw_x, &draw_y);

	/*
	 * Must modify the mask origin within the graphics context to line up
	 * with the bitmap's origin (in order to make bitmaps with
	 * "-background {}" work right).
	 */

 	XSetClipOrigin(display, tdi->gc, draw_x - pix_x, draw_y - pix_y);
	XCopyArea(display, tdi->pm, drawable, tdi->gc, pix_x, pix_y,(unsigned int) pix_w,
		(unsigned int) pix_h, draw_x, draw_y); 
	
	XSetClipOrigin(display, tdi->gc, 0, 0);
    }
}
    

/* subvert the point function so that it now
   now provides world coords for xhair feedback in depth.tcl
*/

static double template_point(Tk_Canvas canvas, Tk_Item *itemPtr, double *coordPtr) {
    TemplateDisplayItem *tdi = (TemplateDisplayItem *) itemPtr;
    int height2;
    double ax, bx;

    tdi->px = coordPtr[0];
    tdi->py = coordPtr[1];
    
    tdi->py += tdi->wy0;

    height2 = tdi->height / 2;

    if (tdi->sep_by_strand) {
	if (tdi->py > height2)
	    tdi->py = tdi->py - height2;
	else
	    tdi->py = height2 - tdi->py;
    }

    if (tdi->ymode == 1) tdi->py /= 10;

    tdi->py = tdi->py / (tdi->yzoom / 200.0) + tdi->yoffset - 50;
    
    if (tdi->logy && tdi->ymode != 1) {
	tdi->py = exp(tdi->py / 50.0) - 1;
    }
    
    tdi->py++;
    
    ax = tdi->width / (tdi->wx1 - tdi->wx0);
    bx = tdi->wx0;
    
    tdi->px = (tdi->px / ax) + bx;

    return 0;
}

   
/* determines whether an item lies entirely
   inside, entirely outside, or overlapping a given rectangle
*/

static int template_area(Tk_Canvas canvas, Tk_Item *itemPtr, double *rectPtr) {
    TemplateDisplayItem *tdi = (TemplateDisplayItem *) itemPtr;

    if ((rectPtr[2] <= tdi->header.x1)
	    || (rectPtr[0] >= tdi->header.x2)
	    || (rectPtr[3] <= tdi->header.y1)
	    || (rectPtr[1] >= tdi->header.y2)) {
	return -1;
    }
    if ((rectPtr[0] <= tdi->header.x1)
	    && (rectPtr[1] <= tdi->header.y1)
	    && (rectPtr[2] >= tdi->header.x2)
	    && (rectPtr[3] >= tdi->header.y2)) {
	return 1;
    }
    
    return 0;
}


/* unimplemented scaling function */

static void scale_template(Tk_Canvas canvas, Tk_Item *itemPtr,
    	    	    	    double originX, double originY,
			    double scaleX, double scaleY) {

    printf("TemplateDisplayItem scale - not implemented\n");
}


/* move item by given amount */

static void translate_template(Tk_Canvas canvas, Tk_Item *itemPtr, 
    	    	    	    	double deltaX, double deltaY) {

    TemplateDisplayItem *tdi = (TemplateDisplayItem *) itemPtr;

    tdi->an_x += deltaX;
    tdi->an_y += deltaY;
    
    compute_template_bbox(canvas, tdi);
}


/* unimplemented function to generate postscript data for printing */

static int template_postscript(Tcl_Interp *interp, Tk_Canvas canvas, Tk_Item *itemPtr,	
    	    	    	    	int prepass) {
    printf("TemplateDisplayItem postscript - not implemented\n");
    
    /* leave blank for now */
    return TCL_ERROR;
}



/* set up the image to draw the template on */

static int initialise_template_image(TemplateDisplayItem *tdi, Display *display) {
    int i;

    if (NULL == (tdi->image = initialise_image(display))) {
	printf("Unable to initialise image_t\n");
	return 0;
    }

    for (i = 0; i < 32; i++) {
	add_colour(tdi->image, 64+i*5, 64+i*5, 64+i*5);
    }

    tdi->background 	  = add_colour(tdi->image, 0, 0, 0);	     // black
    tdi->span_col         = add_colour(tdi->image, 255, 165, 0);       // orange
    tdi->single_col       = add_colour(tdi->image, 0, 0, 255);         // blue
    tdi->inconsistent_col = add_colour(tdi->image, 255, 0, 0);         // red
    tdi->fwd_col  = tdi->fwd_col3 = add_colour(tdi->image, 0, 139, 0);   // green4
    tdi->rev_col  = tdi->rev_col3 = add_colour(tdi->image, 255, 0, 255); // magenta

    return 1;
}
   
    
/* -------------------------------------------------------------------------
 * Y Coordinate allocation routines
 */

SPLAY_HEAD(XTREE, xy_pair);
SPLAY_HEAD(YTREE, xy_pair);
SPLAY_PROTOTYPE(XTREE, xy_pair, x_link, x_cmp);
SPLAY_PROTOTYPE(YTREE, xy_pair, y_link, y_cmp);
SPLAY_GENERATE(XTREE, xy_pair, x_link, x_cmp);
SPLAY_GENERATE(YTREE, xy_pair, y_link, y_cmp);

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
 
static int compute_ypos(TemplateDisplayItem  *tdi, int xgap, tline *tl, int ntl) {
    int i;
    struct xy_pair *node, *curr, *next;
    int yn = 0;
    int min_sz = INT_MIN;
    int max_sz = tdi->min_sz;
    int yoffset = 0, ymax = 0, nleft;

    /* Create and initialise X and Y trees */
    struct XTREE xtree = SPLAY_INITIALIZER(&xtree);
    struct YTREE ytree = SPLAY_INITIALIZER(&ytree);

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

	    if ((node = SPLAY_MIN(XTREE, &xtree)) != NULL && tl[i].x[0] >= node->x) {
		int try_cull = 0;

		/* We found a node, is it the smallest in y? */
		curr = SPLAY_NEXT(XTREE, &xtree, node);
		while (curr && tl[i].x[0] >= curr->x) {
		    if (node->y > curr->y)
			node = curr;
		    curr = SPLAY_NEXT(XTREE, &xtree, curr);
		}
	    
		/* Shift non-smallest y (but < x) to Y-tree */
		curr = SPLAY_MIN(XTREE, &xtree);
		while (curr && tl[i].x[0] >= curr->x) {
		    next = SPLAY_NEXT(XTREE, &xtree, curr);
		    if (curr != node) {
			SPLAY_REMOVE(XTREE, &xtree, curr);
			SPLAY_INSERT(YTREE, &ytree, curr);
			try_cull = 1;
		    }
		    curr = next;
		}
	    
		tl[i].y = node->y + yoffset;
		SPLAY_REMOVE(XTREE, &xtree, node);
		node->x = tl[i].x[3] + xgap;
		SPLAY_INSERT(XTREE, &xtree, node);
	    } else {
		/* Check if we have a free y on ytree */
		if ((node = SPLAY_MIN(YTREE, &ytree)) != NULL) {
		    SPLAY_REMOVE(YTREE, &ytree, node);
		} else {
		    node = (struct xy_pair *)malloc(sizeof(*node));
		    node->y = ++yn;
		    if (ymax < yn)
			ymax = yn;
		}
		tl[i].y = node->y + yoffset;
		node->x = tl[i].x[3] + xgap;
		SPLAY_INSERT(XTREE, &xtree, node);
	    }
	}

	/* Delete trees */
	for (node = SPLAY_MIN(XTREE, &xtree); node; node = next) {
	    next = SPLAY_NEXT(XTREE, &xtree, node);
	    SPLAY_REMOVE(XTREE, &xtree, node);
	}

	for (node = SPLAY_MIN(YTREE, &ytree); node; node = next) {
	    next = SPLAY_NEXT(YTREE, &ytree, node);
	    SPLAY_REMOVE(YTREE, &ytree, node);
	}

	min_sz = max_sz;
	max_sz += tdi->min_sz;
	
	if (max_sz < min_sz+10)
	    max_sz = min_sz+10;
	yoffset += ymax;
    }
    
    return 0;
}


static int sort_tline_by_x(const void *p1, const void *p2) {
    tline *r1 = (tline *)p1;
    tline *r2 = (tline *)p2;

    return r1->x[0] - r2->x[1];
}

/* do the actual work of drawing the template track, uses the gap_range
   functions for most of the data handling */
   	    
static void redraw_template_image(TemplateDisplayItem *tdi, Display *display) {
    double working_wx0, working_wx1;
    int force_change = 0;
    int mode;
    double ax, bx, ay, by;
    int fwd_col, rev_col;
    int half_height;
    int i;
    int ymin = INT_MAX;
    int ymax = INT_MIN;
    static int last_zoom = 0;
    int tsize = MIN(template_max_size(tdi->gr->io), GR_WINDOW_RANGE);
    
    image_remove(tdi->image);
    if(!create_image_buffer(tdi->image, tdi->width, tdi->height, tdi->background)) return;

    // use some values beyond the window size.
    working_wx0 = tdi->wx0 - tsize;
    working_wx1 = tdi->wx1 + tsize;
 
    mode = tdi->reads_only ? 0 : CSIR_PAIR;
    
    if (tdi->ymode == 1) mode |= CSIR_SORT_BY_Y;
    
    set_filter(tdi->gr, tdi->filter, tdi->min_qual, tdi->max_qual, tdi->cmode, tdi->accuracy);
    
    if (gap_range_recalculate(tdi->gr, tdi->width, working_wx0, working_wx1, mode, force_change)) {
	if (tdi->gr->r == NULL) {
	    // either lack of memory or an empty part of contig, blank to black 
	    create_image_from_buffer(tdi->image);
	    XPutImage(display, (Drawable)tdi->pm, tdi->gc, tdi->image->img, 0, 0, 0, 0, tdi->width, tdi->height);
	    return;
	}
	
	force_change = 1;
    }
    
    fwd_col = tdi->yzoom >= 150 ? tdi->fwd_col3 : tdi->fwd_col;
    rev_col = tdi->yzoom >= 150 ? tdi->rev_col3 : tdi->rev_col;
    
    /* world to pixmap conversion values */
    if (tdi->wx1 - tdi->wx0 == 0) return;
    
    ax = tdi->width / (tdi->wx1 - tdi->wx0);
    bx = tdi->wx0;
    
    ay = tdi->height / (tdi->wy1 - tdi->wy0);
    by = tdi->wy0;    
    
    /* 1) Compute X */
    tdi->ntl = gap_range_x(tdi->gr, ax, bx, fwd_col, rev_col, 
    	    	    	    tdi->single_col, tdi->span_col, tdi->inconsistent_col,
		    	    force_change, tdi->reads_only);
			    
    /* 2) Compute Y coordinates (part 1) */
    if (tdi->ymode == 1) {
    	double xgap = 100;
    
	qsort(tdi->gr->tl, tdi->ntl, sizeof(tline), sort_tline_by_x);
	compute_ypos(tdi, xgap, tdi->gr->tl, tdi->ntl);
    }
    
    /* 3) Plot the lines */
    half_height = tdi->height / 2;
    tdi->yz = tdi->yzoom / 200.0;
    
    for (i = 0; i < tdi->ntl; i++) {
    	int j;
    	tline *tl = &tdi->gr->tl[i];
    
	for (j = 0; j < 3; j++) {
	
	    if (tl->x[j] < tl->x[j+1]) {
		double y;

		/* Compute the Y coordinates (part 2) */
		/* See XHAIR sub-command code too; keep it in sync */
		if (tdi->ymode == 1) {
		    y = tl->y * 10;
		} else {
		    if (tdi->ymode == 0)
			y = tl->x[3] - tl->x[0];
		    else
			y = tl->mq * 4;
			
		    if (tdi->logy) {
			if (y < 0) y = 0;
			
			y = 50 * log(y + 1);
		    }
		}
		
		y = (y + 50 - tdi->yoffset) * tdi->yz;

		if (tdi->spread)
		    y = y - tdi->spread / 2 + ((tl->x[0] + tl->rec) % (tdi->spread));

		if (tdi->sep_by_strand)
		    y = tl->t_strand ? half_height - y : half_height + y;

		if (ymin > y) ymin = y;
		if (ymax < y) ymax = y;

		/* And plot if visible */
		if (y >= tdi->wy0 && y <= tdi->wy1) {
		    int rx1, rx2, ry;
		    
		    rx1 = (tl->x[j] - bx) * ax; // world to raster conversion
		    rx2 = (tl->x[j + 1] - bx) * ax;
		    ry  = (y - by) * ay;
		    
		    draw_line(tdi->image, rx1, rx2, ry, tl->col[j]);
		}
	    }
	}
    }
    
    create_image_from_buffer(tdi->image);
    XPutImage(display, (Drawable)tdi->pm, tdi->gc, tdi->image->img, 0, 0, 0, 0, tdi->width, tdi->height);
    
    /* some last bits of size calculation for scrolling */
    tdi->y_start = ymin - 10;
    tdi->y_end   = ymax + 10;
    
    if (tdi->y_start > 0) tdi->y_start = 0;
    if (tdi->y_end < tdi->height) tdi->y_end = tdi->height;
    
    last_zoom = tdi->yzoom;
}
    
    
    
    
    
    
   

	    

    
    
    
    
    
    
