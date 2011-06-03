/*
 *
 * depth_track.c - a tk canvas item to display depth track data
 *                 formerly used tkRaster to do the same job.
 *
 * Andrew Whitwham, July 2010
 *
 */

#include <stdio.h>
#include <math.h>
#include <tcl.h>
#include <X11/Xlib.h>

#include "gap_range.h"
#include "depth_track.h"

/* define the depth track item */

typedef struct DepthTrackItem {
    Tk_Item header; 	    /* mandatory entry */
    GC template;    	    /* template context */
    GC reads;	    	    /* reads context */
    GC copy;
    XColor *tcolour;
    XColor *scolour;
    XColor *background;
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
    int accuracy;
    int reads_only;
    int filter;
    int min_qual;
    int max_qual;
    int min_sz;
    double yzoom;
    double yz;
    
    int width, height;      /* pixmap width and height */
    int ntl;	    	    /* num lines in viewable area */
} DepthTrackItem;


/* mandatory prototypes */

static int		depth_coords(Tcl_Interp *interp, Tk_Canvas canvas,
    	    	    	    Tk_Item *itemPtr, int argc, Tcl_Obj *CONST argv[]);
static int		depth_area(Tk_Canvas canvas, Tk_Item *itemPtr, 
    	    	    	    double *rectPtr);
static double		depth_point(Tk_Canvas canvas, Tk_Item *itemPtr,
    	    	    	    double *coordPtr);
static int		depth_postscript(Tcl_Interp *interp, 
    	    	    	    Tk_Canvas canvas, Tk_Item *itemPtr, int prepass);
static int		configure_depth(Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr, int argc,
			    Tcl_Obj *CONST argv[], int flags);
static int		create_depth(Tcl_Interp *interp,
			    Tk_Canvas canvas, struct Tk_Item *itemPtr,
			    int argc, Tcl_Obj *CONST argv[]);
static void		delete_depth(Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display);
static void		depth_display(Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display, Drawable dst,
			    int x, int y, int width, int height);
static void		scale_depth(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double originX, double originY,
			    double scaleX, double scaleY);
static void		translate_depth(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double deltaX, double deltaY);

/* none mandatory protoypes */
static void compute_depth_bbox(Tk_Canvas canvas, DepthTrackItem *dti);

static void redraw_depth_image(DepthTrackItem *dti, Display *display);

/* config options, some are meant only to be returned with itemcget */
static Tk_ConfigSpec config_specs[] = {
    {TK_CONFIG_ANCHOR, "-anchor", NULL, NULL,
	"nw", Tk_Offset(DepthTrackItem, anchor), TK_CONFIG_DONT_SET_DEFAULT},
    {TK_CONFIG_CUSTOM, "-range", NULL, NULL,
    	NULL, Tk_Offset(DepthTrackItem, gr), TK_CONFIG_NULL_OK, &range_option},
    {TK_CONFIG_INT, "-logy", "logy", "LogY", "0", Tk_Offset(DepthTrackItem, logy), 0, 0},
    {TK_CONFIG_INT, "-cmode", "colourMode", "ColourMode", "0", Tk_Offset(DepthTrackItem, cmode), 0, 0},
    {TK_CONFIG_INT, "-ymode", "YMode", "YMode", "0", Tk_Offset(DepthTrackItem, ymode), 0, 0},
    {TK_CONFIG_INT, "-accuracy", "accuracy", "Accuracy", "0", Tk_Offset(DepthTrackItem, accuracy), 0, 0},
    {TK_CONFIG_INT, "-reads_only", "readsOnly", "ReadsOnly", "0", Tk_Offset(DepthTrackItem, reads_only), 0, 0},
    {TK_CONFIG_INT, "-filter", "filter", "Filter", "0", Tk_Offset(DepthTrackItem, filter), 0, 0},
    {TK_CONFIG_INT, "-min_qual", "minQual", "MinQual", "0", Tk_Offset(DepthTrackItem, min_qual), 0, 0},
    {TK_CONFIG_INT, "-max_qual", "maxQual", "MaxQual", "255", Tk_Offset(DepthTrackItem, max_qual), 0, 0},
    {TK_CONFIG_DOUBLE, "-yzoom", "yZoom", "YZoom", "100.0", Tk_Offset(DepthTrackItem, yzoom), 0, 0},
    {TK_CONFIG_DOUBLE, "-yz", "yZ", "YZ", "0", Tk_Offset(DepthTrackItem, yz), 0, 0},
    {TK_CONFIG_DOUBLE, "-wx0", NULL, NULL, "0", Tk_Offset(DepthTrackItem, wx0), 0}, 
    {TK_CONFIG_DOUBLE, "-wx1", NULL, NULL, "0", Tk_Offset(DepthTrackItem, wx1), 0},
    {TK_CONFIG_DOUBLE, "-wy0", NULL, NULL, "0", Tk_Offset(DepthTrackItem, wy0), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-wy1", NULL, NULL, "0", Tk_Offset(DepthTrackItem, wy1), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-y_start", NULL, NULL, "0", Tk_Offset(DepthTrackItem, y_start), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-y_end", NULL, NULL, "0", Tk_Offset(DepthTrackItem, y_end), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-px", NULL, NULL, "0", Tk_Offset(DepthTrackItem, px), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-py", NULL, NULL, "0", Tk_Offset(DepthTrackItem, py), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_COLOR, "-tcolour", NULL, NULL,
	"magenta", Tk_Offset(DepthTrackItem, tcolour), TK_CONFIG_NULL_OK},
    {TK_CONFIG_COLOR, "-scolour", NULL, NULL,
	"green4", Tk_Offset(DepthTrackItem, scolour), TK_CONFIG_NULL_OK},
    {TK_CONFIG_COLOR, "-background", NULL, NULL,
        "black", Tk_Offset(DepthTrackItem, background), TK_CONFIG_NULL_OK},
    {TK_CONFIG_END, (char *) NULL, (char *) NULL, (char *) NULL, (char *) NULL, 0, 0}
};			    

/* define the DepthTrackItem functions and name */			    
Tk_ItemType tkDepthItem = {
    "depth_track",
    sizeof(DepthTrackItem),
    create_depth,
    config_specs,
    configure_depth,
    depth_coords,
    delete_depth,
    depth_display,
    TK_CONFIG_OBJS,          /* needs to be this */
    depth_point,
    depth_area,
    depth_postscript,
    scale_depth,
    translate_depth,
    (Tk_ItemIndexProc *)     NULL,
    (Tk_ItemCursorProc *)    NULL,
    (Tk_ItemSelectionProc *) NULL,
    (Tk_ItemInsertProc *)    NULL,
    (Tk_ItemDCharsProc *)    NULL,
    (Tk_ItemType *)          NULL,
};
	    
	    
/* this function is called when item is created as part
   of a canvas */
   			    
static int create_depth(Tcl_Interp *interp,
	    Tk_Canvas canvas,
	    Tk_Item *itemPtr,
	    int argc,
	    Tcl_Obj *CONST *argv) {
    DepthTrackItem *dti = (DepthTrackItem *)itemPtr;
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
    
    dti->template = None;
    dti->reads = None;
    dti->tcolour = NULL;
    dti->scolour = NULL;
    dti->copy = NULL;
    dti->background = NULL;
    dti->pm = None;
    dti->anchor = TK_ANCHOR_NW;
    dti->an_x = 0.0;
    dti->an_y = 0.0;
    dti->gr = NULL;
    dti->wy0 = dti->y_start = 0;
    dti->wy1 = dti->y_end = Tk_Height(Tk_CanvasTkwin(canvas)); /* initial world height */
    dti->width = -1;
    dti->height = -1;
    
    
    if ((depth_coords(interp, canvas, itemPtr, i, argv) == TCL_OK)) {
    	if (configure_depth(interp, canvas, itemPtr, argc - i, argv + i, 0) == TCL_OK) {
	    // possibly more initialisation here

	    return TCL_OK;
	}
    }
    
    /* if we get here then something is wrong */
    printf("DepthTrackItem creation failed\n");
    delete_depth(canvas, itemPtr, Tk_Display(Tk_CanvasTkwin(canvas)));
    
    return TCL_ERROR;
}


/* invoked by the coords command */

static int depth_coords(Tcl_Interp *interp, Tk_Canvas canvas,		
    	    	    	    Tk_Item *itemPtr, int objc, Tcl_Obj *CONST objv[]) {
    DepthTrackItem *dti = (DepthTrackItem *)itemPtr;
    
    if (objc == 0) {
	Tcl_Obj *obj = Tcl_NewObj();

	Tcl_Obj *subobj = Tcl_NewDoubleObj(dti->an_x);
    	Tcl_ListObjAppendElement(interp, obj, subobj);

	subobj = Tcl_NewDoubleObj(dti->an_y);
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
		&dti->an_x) != TCL_OK)
		|| (Tk_CanvasGetCoordFromObj(interp, canvas, objv[1],
			&dti->an_y) != TCL_OK)) {
	    return TCL_ERROR;
	}
	
	compute_depth_bbox(canvas, dti);
    } else {
	char buf[64 + TCL_INTEGER_SPACE];

	sprintf(buf, "wrong # coordinates: expected 0 or 2, got %d", objc);
	Tcl_SetResult(interp, buf, TCL_VOLATILE);
	return TCL_ERROR;
    }
    
    return TCL_OK;
}


/* configure and redraw the depth track */

static int configure_depth(Tcl_Interp *interp,
	       Tk_Canvas canvas,
	       Tk_Item *itemPtr,
	       int argc,
	       Tcl_Obj *CONST *argv,
	       int flags) {
	    
    DepthTrackItem *dti = (DepthTrackItem *)itemPtr;
    Tk_Window tkwin;
    Display *display;
    int width, height;
        
    tkwin   = Tk_CanvasTkwin(canvas);
    display = Tk_Display(tkwin);
    
    if (Tk_ConfigureWidget(interp, tkwin, config_specs, argc, (char **) argv,
	    (char *) dti, flags|TK_CONFIG_OBJS) != TCL_OK) {

	printf("ERROR %s\n", Tcl_GetStringResult(interp));
	return TCL_ERROR;
    }
    
    width  = Tk_Width(tkwin);
    height = Tk_Height(tkwin);
    
    if (dti->template == None) {
    	XGCValues gcValues;
    	gcValues.foreground = dti->tcolour->pixel;
	gcValues.background = dti->background->pixel;
	dti->template = Tk_GetGC(tkwin, GCForeground|GCBackground, &gcValues);
    }

    if (dti->reads == None) {
    	XGCValues gcValues;
    	gcValues.foreground = dti->scolour->pixel;
	gcValues.background = dti->background->pixel;
	dti->reads = Tk_GetGC(tkwin, GCForeground|GCBackground, &gcValues);
    }

    if (dti->copy == None) {
    	XGCValues gcValues;
	gcValues.foreground = dti->background->pixel;
	gcValues.background = dti->background->pixel;
	gcValues.function = GXcopy;
	gcValues.graphics_exposures = False;
	dti->copy = Tk_GetGC(tkwin, GCForeground|GCBackground|GCGraphicsExposures|GCFunction, &gcValues);
    }

    if (width != dti->width || height != dti->height) {
    	/* going to need a new pixmap */
 	if (dti->pm) {
	    Tk_FreePixmap(display, dti->pm);
	}

    	dti->width = width;
    	dti->height = height;	
	
	dti->pm = Tk_GetPixmap(display, Tk_WindowId(tkwin), dti->width, dti->height, DefaultDepthOfScreen(Tk_Screen(tkwin)));
    }
    

    redraw_depth_image(dti, display);
    compute_depth_bbox(canvas, dti);
    return TCL_OK;
}


/* free resources nicely*/

static void delete_depth(Tk_Canvas canvas, Tk_Item *itemPtr, Display *display) {
    DepthTrackItem *dti = (DepthTrackItem *)itemPtr;
    
    if (dti->pm != None) {
    	Tk_FreePixmap(display, dti->pm);
    }
    
    if (dti->template != None) {
    	Tk_FreeGC(display, dti->template);
    }
    
    if (dti->reads != None) {
    	Tk_FreeGC(display, dti->reads);
    }
    
    if (dti->copy != None) {
    	Tk_FreeGC(display, dti->copy);
    }
    
    if (dti->tcolour != NULL) {
	Tk_FreeColor(dti->tcolour);
    }

    if (dti->scolour != NULL) {
	Tk_FreeColor(dti->scolour);
    }

    if (dti->background != NULL) {
	Tk_FreeColor(dti->background);
    }
}


/* compute the bounding box of the item */

static void compute_depth_bbox(Tk_Canvas canvas, DepthTrackItem *dti) {
    Tk_Window tkwin;
    int width, height;
    int x, y;
    
    tkwin = Tk_CanvasTkwin(canvas);

    width  = Tk_Width(tkwin);
    height = Tk_Height(tkwin);
    
    x = (int) (dti->an_x + ((dti->an_x >= 0) ? 0.5 : - 0.5));
    y = (int) (dti->an_y + ((dti->an_y >= 0) ? 0.5 : - 0.5));
    
    switch (dti->anchor) {
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
    
    dti->header.x1 = x;
    dti->header.y1 = y;
    dti->header.x2 = x + width;
    dti->header.y2 = y + height;
}

			    
/* draw the depth track onto the screen, copies from the pixmap to
   the canvas */			    
			    
static void depth_display(Tk_Canvas canvas, Tk_Item *itemPtr,
    	    	    	    Display *display, Drawable drawable,
			    int x, int y, int width, int height) {
    
    DepthTrackItem *dti = (DepthTrackItem *)itemPtr;
    int pix_x, pix_y, pix_w, pix_h;
    short draw_x, draw_y;

    if (dti->pm != None) {
	if (x > dti->header.x1) {
	    pix_x = x - dti->header.x1;
	    pix_w = dti->header.x2 - x;
	} else {
	    pix_x = 0;
	    if ((x + width) < dti->header.x2) {
		pix_w = x + width - dti->header.x1;
	    } else {
		pix_w = dti->header.x2 - dti->header.x1;
	    }
	}
	if (y > dti->header.y1) {
	    pix_y = y - dti->header.y1;
	    pix_h = dti->header.y2 - y;
	} else {
	    pix_y = 0;
	    if ((y + height) < dti->header.y2) {
		pix_h = y + height - dti->header.y1;
	    } else {
		pix_h = dti->header.y2 - dti->header.y1;
	    }
	}
	
	Tk_CanvasDrawableCoords(canvas,
		(double) (dti->header.x1 + pix_x),
		(double) (dti->header.y1 + pix_y),
		&draw_x, &draw_y);

	/*
	 * Must modify the mask origin within the graphics context to line up
	 * with the bitmap's origin (in order to make bitmaps with
	 * "-background {}" work right).
	 */

	XSetClipOrigin(display, dti->copy, draw_x - pix_x, draw_y - pix_y);
	XCopyArea(display, dti->pm, drawable, dti->copy, pix_x, pix_y,(unsigned int) pix_w,
		  (unsigned int) pix_h, draw_x, draw_y); 
	
	XSetClipOrigin(display, dti->copy, 0, 0);
    }
}


/* subvert the point function so that it now
   now provides world coords for xhair feedback in depth.tcl
*/

static double depth_point(Tk_Canvas canvas, Tk_Item *itemPtr, double *coordPtr) {
    DepthTrackItem *dti = (DepthTrackItem *)itemPtr;
    double ax, bx;

    dti->px = coordPtr[0];
    dti->py = coordPtr[1];

    dti->py = (dti->height - dti->py) / dti->yz;
    
    ax = dti->width / (dti->wx1 - dti->wx0);
    bx = dti->wx0;
    
    dti->px = (dti->px / ax) + bx;

    return 0;
}


/* determines whether an item lies entirely
   inside, entirely outside, or overlapping a given rectangle
*/

static int depth_area(Tk_Canvas canvas, Tk_Item *itemPtr, double *rectPtr) {
    DepthTrackItem *dti = (DepthTrackItem *)itemPtr;

    if ((rectPtr[2] <= dti->header.x1)
	    || (rectPtr[0] >= dti->header.x2)
	    || (rectPtr[3] <= dti->header.y1)
	    || (rectPtr[1] >= dti->header.y2)) {
	return -1;
    }
    if ((rectPtr[0] <= dti->header.x1)
	    && (rectPtr[1] <= dti->header.y1)
	    && (rectPtr[2] >= dti->header.x2)
	    && (rectPtr[3] >= dti->header.y2)) {
	return 1;
    }
    
    return 0;
}


/* unimplemented scaling function */

static void scale_depth(Tk_Canvas canvas, Tk_Item *itemPtr,
    	    	    	    double originX, double originY,
			    double scaleX, double scaleY) {
    printf("DepthTrackItem scale - not implemented.\n");
}


/* move item by given amount */

static void translate_depth(Tk_Canvas canvas, Tk_Item *itemPtr, 
    	    	    	    	double deltaX, double deltaY) {

    DepthTrackItem *dti = (DepthTrackItem *)itemPtr;

    dti->an_x += deltaX;
    dti->an_y += deltaY;
    
    compute_depth_bbox(canvas, dti);
}


/* unimplemented function to generate postscript data for printing */

static int depth_postscript(Tcl_Interp *interp, Tk_Canvas canvas, Tk_Item *itemPtr,	
    	    	    	    	int prepass) {
    printf("DepthTrackItem postscript - not implemented\n");
    
    return TCL_ERROR;
}


/* do the actual work of drawing the depth track, uses the gap_range
   functions for most of the data handling */
   	    
static void redraw_depth_image(DepthTrackItem *dti, Display *display) {
    double working_wx0, working_wx1;
    int force_change = 0;
    XPoint *sp = NULL;
    XPoint *tp = NULL;
    int i;
    int ymin = INT_MAX;
    int ymax = INT_MIN;
    double ax, bx;

    /* clear the pixmap */
    XFillRectangle(display, dti->pm, dti->copy, 0, 0, dti->width, dti->height);
    
    working_wx0 = dti->wx0 - GR_WINDOW_RANGE; // use some values beyond the window size.
    working_wx1 = dti->wx1 + GR_WINDOW_RANGE;
    
    if (gap_range_recalculate(dti->gr, dti->width, working_wx0, working_wx1, dti->gr->template_mode, force_change)) {

	if (dti->gr->r == NULL) {
	    return;
	}
	
	force_change = 1;
    }
    
    /* world to pixmap conversion values */
    if (dti->wx1 - dti->wx0 == 0) return;
    
    ax = dti->width / (dti->wx1 - dti->wx0);
    bx = dti->wx0;
    
    /* 1) Compute X */
    dti->ntl = gap_range_x(dti->gr, ax, bx, 0, 0, 0, 0, 0, force_change, dti->reads_only);

    /* 2) Calculate Y scale */
    if (dti->gr->max_height > 0) {    
    	dti->yz = ((double)dti->height / (double)dti->gr->max_height) * 0.95;
    } else {
    	dti->yz = dti->yzoom / 50.0;
    }
    
   /* Plot depth */
    sp = (XPoint *)calloc(dti->width, sizeof(*sp));
    tp = (XPoint *)calloc(dti->width, sizeof(*tp));
    
    // FIXME - need to put in a memory check
    
    for (i = 0; i < dti->width; i++) {
	sp[i].x = i;
	sp[i].y = dti->height - 5 - dti->gr->depth[i].s * dti->yz;
	tp[i].x = i;
	tp[i].y = dti->height - 5 - dti->gr->depth[i].t * dti->yz;
	
	if (ymin > sp[i].y) ymin = sp[i].y;
	if (ymin > tp[i].y) ymin = tp[i].y;
	
	if (ymax < sp[i].y) ymax = sp[i].y;
	if (ymax < tp[i].y) ymax = tp[i].y;
    }
    
    XDrawLines(display, dti->pm, dti->reads, sp, dti->width, CoordModeOrigin);
        
    if (!dti->reads_only) {
	XDrawLines(display, dti->pm, dti->template, tp, dti->width, CoordModeOrigin);
    }
    
    /* Compute scrollbar and bounding box locations */
    dti->y_start = ymin != INT_MAX ? ymin : 0;
    dti->y_end = ymax != INT_MIN ? ymax : 0;
    
    free(sp);
    free(tp);
}
   
