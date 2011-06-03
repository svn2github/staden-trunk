/*
 * A quality plot canvas item type.
 */
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "tg_gio.h"
#include "tg_struct.h"
#include "quality_plot.h"
#include "gap_range.h"
#include "misc.h"
#include "consensus.h"

#define MAX_QUAL 150

/* --------------------------------------------------------------------------
 * Structs, prototypes and configurations
 */

typedef struct QualityPlot {
    /* mandatory entry */
    Tk_Item header; 	    

    /* GC / Colors */
    GC qualGC;
    GC hetGC;
    GC discrepGC;
    GC copyGC;

    XColor *qcolour;
    XColor *hcolour;
    XColor *dcolour;

    /* Other plotting bits */
    int discrep;      /* what to plot */
    int quality;
    int hetero;

    Pixmap pm;	    	    /* pixmap to draw on */
    Tk_Anchor anchor;	    /* pixmap anchorpoint */
    double an_x, an_y;      /* postioning points */
    double px, py;          /* current mouse pointer x, y (in world units)*/
    double contig_start;    /* start value of contig */
    double contig_length;   /* length of contig */
    double wx0, wx1;	    /* world coords */
    double wy0, wy1;        
    double y_start, y_end;  /* world edges in y */
    int width, height;      /* pixmap width and height */
    double yz;              /* y zoom */

    /* Gap5 integration */
    GapIO *io;
    tg_rec contig;
    gap_range_t *gr;        /* range info */
} QualityPlot;


/* Prototypes for mandatory entries of Tk_Item */
static int		qplot_create(Tcl_Interp *interp,
			    Tk_Canvas canvas, struct Tk_Item *itemPtr,
			    int objc, Tcl_Obj *CONST objv[]);
static int		qplot_configure(Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr, int objc,
			    Tcl_Obj *CONST objv[], int flags);
static int		qplot_coords(Tcl_Interp *interp, Tk_Canvas canvas,
    	    	    	    Tk_Item *itemPtr, int objc, Tcl_Obj *CONST objv[]);
static void		qplot_delete(Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display);
static void		qplot_display(Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display, Drawable dst,
			    int x, int y, int width, int height);
static double		qplot_point(Tk_Canvas canvas, Tk_Item *itemPtr,
    	    	    	    double *coordPtr);
static int		qplot_area(Tk_Canvas canvas, Tk_Item *itemPtr, 
    	    	    	    double *rectPtr);
static int		qplot_postscript(Tcl_Interp *interp, 
    	    	    	    Tk_Canvas canvas, Tk_Item *itemPtr, int prepass);
static void		qplot_scale(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double originX, double originY,
			    double scaleX, double scaleY);
static void		qplot_translate(Tk_Canvas canvas,
			    Tk_Item *itemPtr, double deltaX, double deltaY);


/*---- define the non-standard option types */ 

/* parse the tgrec option to get at the data */
static int rec_cmd_parse(ClientData clientData, Tcl_Interp *interp,
			 Tk_Window tkwin, char *value,
			 char *widgRec, int offset) {
    tg_rec rec;
    char *endp;
    
    rec = strtol64(value, &endp, 10);
    *((tg_rec *)(widgRec+offset)) = rec;

    return TCL_OK;
}
    

/* Stringify for the rec option */
static char *rec_cmd_print(ClientData clientData, Tk_Window tkwin,
			   char *widgRec, int offset,
			   Tcl_FreeProc **freeProcPtr) {
    static char buf[1024];
    tg_rec rec = *(tg_rec *)(widgRec + offset);

    sprintf(buf, "%"PRIrec, rec);
    *freeProcPtr = NULL;

    return buf;
}

static Tk_CustomOption rec_option = {
    (Tk_OptionParseProc *)rec_cmd_parse,
    rec_cmd_print, (ClientData)NULL
};


/* parse the tgrec option to get at the data */
static int io_cmd_parse(ClientData clientData, Tcl_Interp *interp,
			Tk_Window tkwin, char *value,
			char *widgRec, int offset) {
    GapIO *io;

    if (0 != strncmp(value, "io=", 3))
	return TCL_ERROR;

    if (1 != sscanf(value+3, "%p", &io))
	return TCL_ERROR;

    *((GapIO **)(widgRec+offset)) = io;

    return TCL_OK;
}
    

/* Stringify for the rec option */
static char *io_cmd_print(ClientData clientData, Tk_Window tkwin,
			  char *widgRec, int offset,
			   Tcl_FreeProc **freeProcPtr) {
    static char buf[1024];
    GapIO *io = *(GapIO **)(widgRec + offset);

    sprintf(buf, "io=%p", io);
    *freeProcPtr = NULL;

    return buf;
}

static Tk_CustomOption io_option = {
    (Tk_OptionParseProc *)io_cmd_parse,
    io_cmd_print, (ClientData)NULL
};

/* Config options */
static Tk_ConfigSpec qplot_config_specs[] = { 
    {TK_CONFIG_CUSTOM, "-io", NULL, NULL,
    	NULL, Tk_Offset(QualityPlot, io), TK_CONFIG_NULL_OK, &io_option},
    {TK_CONFIG_CUSTOM, "-contig", NULL, NULL,
    	NULL, Tk_Offset(QualityPlot, contig), TK_CONFIG_NULL_OK, &rec_option},
    {TK_CONFIG_CUSTOM, "-range", NULL, NULL,
    	NULL, Tk_Offset(QualityPlot, gr), TK_CONFIG_NULL_OK, &range_option},
    {TK_CONFIG_ANCHOR, "-anchor", NULL, NULL,
	"nw", Tk_Offset(QualityPlot, anchor), TK_CONFIG_DONT_SET_DEFAULT},
    {TK_CONFIG_COLOR, "-quality_colour", NULL, NULL,
	"yellow", Tk_Offset(QualityPlot,  qcolour), TK_CONFIG_NULL_OK},
    {TK_CONFIG_COLOR, "-heterozygosity_colour", NULL, NULL,
	"magenta", Tk_Offset(QualityPlot, hcolour), TK_CONFIG_NULL_OK},
    {TK_CONFIG_COLOR, "-discrepancies_colour", NULL, NULL,
	"green4", Tk_Offset(QualityPlot,  dcolour), TK_CONFIG_NULL_OK},
    {TK_CONFIG_INT, "-discrepancies", NULL, NULL,
        "0", Tk_Offset(QualityPlot, discrep), 0}, 
    {TK_CONFIG_INT, "-quality", NULL, NULL,
        "0", Tk_Offset(QualityPlot, quality), 0}, 
    {TK_CONFIG_INT, "-heterozygosity", NULL, NULL,
        "0", Tk_Offset(QualityPlot, hetero), 0}, 
    {TK_CONFIG_DOUBLE, "-wx0", NULL, NULL,
        "0", Tk_Offset(QualityPlot, wx0), 0}, 
    {TK_CONFIG_DOUBLE, "-wx1", NULL, NULL,
        "0", Tk_Offset(QualityPlot, wx1), 0},
    {TK_CONFIG_DOUBLE, "-wy0", NULL, NULL,
        "0", Tk_Offset(QualityPlot, wy0), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-wy1", NULL, NULL,
        "0", Tk_Offset(QualityPlot, wy1), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-px", NULL, NULL,
        "0", Tk_Offset(QualityPlot, px), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-py", NULL, NULL,
        "0", Tk_Offset(QualityPlot, py), TK_CONFIG_DONT_SET_DEFAULT}, 
    {TK_CONFIG_DOUBLE, "-yz", "yZ", "YZ",
        "0", Tk_Offset(QualityPlot, yz), 0, 0},
    {TK_CONFIG_END, (char *) NULL, (char *) NULL,
        (char *) NULL, (char *) NULL, 0, 0}
};


/* Define the canvas item type */
Tk_ItemType consQualityItem = {
    "quality_plot",
    sizeof(QualityPlot),
    qplot_create,
    qplot_config_specs,
    qplot_configure,
    qplot_coords,
    qplot_delete,
    qplot_display,
    TK_CONFIG_OBJS, /* alwaysRedraw: whether callbacks accept Tcl_Obj */
    qplot_point,
    qplot_area,
    qplot_postscript,
    qplot_scale,
    qplot_translate,
    (Tk_ItemIndexProc *)     NULL,
    (Tk_ItemCursorProc *)    NULL,
    (Tk_ItemSelectionProc *) NULL,
    (Tk_ItemInsertProc *)    NULL,
    (Tk_ItemDCharsProc *)    NULL,
    (Tk_ItemType *)          NULL,
};


/* ---------------------------------------------------------------------------
 * Utility functions used by the main canvas below interface
 */

/* compute the bounding box of the item and return in the Tk_Item header */
static void update_bbox(Tk_Canvas canvas, QualityPlot *qp) {
    Tk_Window tkwin;
    int width, height;
    int x, y;
    
    tkwin = Tk_CanvasTkwin(canvas);

    width  = Tk_Width(tkwin);
    height = Tk_Height(tkwin);
    
    x = (int) (qp->an_x + ((qp->an_x >= 0) ? 0.5 : - 0.5));
    y = (int) (qp->an_y + ((qp->an_y >= 0) ? 0.5 : - 0.5));
    
    switch (qp->anchor) {
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
    
    qp->header.x1 = x;
    qp->header.y1 = y;
    qp->header.x2 = x + width;
    qp->header.y2 = y + height;
}


static void qplot_redraw(QualityPlot *qp, Display *display) {
    XPoint *xp = NULL;
    XSegment *xs = NULL;
    consensus_t *cons;
    int i, len;
    float yz, xz;
    int mode = 0;

    if (!qp->pm)
	return;

    /* clear the pixmap */
    XFillRectangle(display, qp->pm, qp->copyGC, 0, 0, qp->width, qp->height);

    if (qp->hetero)
	mode |= CONS_SCORES;
    if (qp->discrep)
	mode |= CONS_DISCREP;
    if (!mode && !qp->quality)
	return;

    /* Fetch consensus */
    len = (int)qp->wx1 - (int)qp->wx0 + 1;
    cons = malloc(len * sizeof(*cons));
    //calculate_consensus_fast(qp->io, qp->contig, qp->wx0, qp->wx1, cons);
    calculate_consensus_bit_het(qp->io, qp->contig, qp->wx0, qp->wx1, mode,
    				qp->gr->r, qp->gr->nr, cons);

    yz = qp->height / (MAX_QUAL+1 + 10.0);
    xz = len / (double)qp->width;
    qp->yz = yz;

    /* Plot it */
    xp = malloc(qp->width * sizeof(*xp));
    xs = malloc(qp->width * sizeof(*xs));

    if (qp->quality) {
	if (xz <= 1) {
	    for (i = 0; i < qp->width; i++) {
		int J = i * (qp->wx1 - qp->wx0) / qp->width;
		int v = cons[J].phred;
		if (v > MAX_QUAL) v = MAX_QUAL;
		if (v < 0)        v = 0;
	
		assert(J >= 0 && J < qp->wx1 - qp->wx0 + 1);
		xp[i].x = i;
		xp[i].y = yz * (MAX_QUAL+1-v + 5);
	    }

	    XDrawLines(display, qp->pm, qp->qualGC, xp, qp->width,
		       CoordModeOrigin);
	} else {
	    int ox = -1;
	    for (i = 0; i < len; i++) {
		int x, v = cons[i].phred;
		if (v > MAX_QUAL) v = MAX_QUAL;
		if (v < 0)        v = 0;
	    
		x = i / xz;
		assert(x >= 0 && x < qp->width);
		if (x == ox) {
		    if (xs[x].y1 < yz * (MAX_QUAL+1-v + 5))
			xs[x].y1 = yz * (MAX_QUAL+1-v + 5);
		    if (xs[x].y2 > yz * (MAX_QUAL+1-v + 5))
			xs[x].y2 = yz * (MAX_QUAL+1-v + 5);
		} else {
		    xs[x].x1 = x;
		    xs[x].x2 = x;
		    xs[x].y1 = yz * (MAX_QUAL+1-v + 5);
		    xs[x].y2 = yz * (MAX_QUAL+1-v + 5);
		    ox = x;
		}
	    }
	    for (i = 1; i < qp->width; i++) {
		if (xs[i].y2 > xs[i-1].y1)
		    xs[i].y2 = xs[i-1].y1;
		else if (xs[i].y1 < xs[i-1].y2)
		    xs[i].y1 = xs[i-1].y2;
	    }

	    XDrawSegments(display, qp->pm, qp->qualGC, xs, qp->width);
	}
    }

    if (qp->hetero) {
	if (xz <= 1) {
	    for (i = 0; i < qp->width; i++) {
		int J = i * (qp->wx1 - qp->wx0) / qp->width;
		int v = cons[J].scores[6];
		//v = cons[J].phred;
		if (v > MAX_QUAL) v = MAX_QUAL;
		if (v < 0)        v = 0;
	
		assert(J >= 0 && J < qp->wx1 - qp->wx0 + 1);
		xp[i].x = i;
		xp[i].y = yz * (MAX_QUAL+1-v + 5);
	    }

	    XDrawLines(display, qp->pm, qp->hetGC, xp, qp->width,
		       CoordModeOrigin);
	} else {
	    int ox = -1;
	    for (i = 0; i < len; i++) {
		int x, v = cons[i].scores[6];
		//v = cons[i].phred;
		if (v > MAX_QUAL) v = MAX_QUAL;
		if (v < 0)        v = 0;
	    
		x = i / xz;
		assert(x >= 0 && x < qp->width);
		if (x == ox) {
		    if (xs[x].y2 > yz * (MAX_QUAL+1-v + 5))
			xs[x].y2 = yz * (MAX_QUAL+1-v + 5);
		} else {
		    xs[x].x1 = x;
		    xs[x].x2 = x;
		    xs[x].y1 = yz * (MAX_QUAL+1 + 5); /* impulses */
		    xs[x].y2 = yz * (MAX_QUAL+1-v + 5);
		    ox = x;
		}
	    }

	    XDrawSegments(display, qp->pm, qp->hetGC, xs, qp->width);
	}
    }

    if (qp->discrep) {
	if (xz <= 1) {
	    for (i = 0; i < qp->width; i++) {
		int J = i * (qp->wx1 - qp->wx0) / qp->width;
		int v = cons[J].discrep * 10;
		if (v > MAX_QUAL) v = MAX_QUAL;
		if (v < 0)        v = 0;
	
		assert(J >= 0 && J < qp->wx1 - qp->wx0 + 1);
		xp[i].x = i;
		xp[i].y = yz * (MAX_QUAL+1-v + 5);
	    }

	    XDrawLines(display, qp->pm, qp->discrepGC, xp, qp->width,
		       CoordModeOrigin);
	} else {
	    int ox = -1;
	    for (i = 0; i < len; i++) {
		int x, v = cons[i].discrep * 10;
		if (v > MAX_QUAL) v = MAX_QUAL;
		if (v < 0)        v = 0;
	    
		x = i / xz;
		assert(x >= 0 && x < qp->width);
		if (x == ox) {
		    if (xs[x].y2 > yz * (MAX_QUAL+1-v + 5))
			xs[x].y2 = yz * (MAX_QUAL+1-v + 5);
		} else {
		    xs[x].x1 = x;
		    xs[x].x2 = x;
		    xs[x].y1 = yz * (MAX_QUAL+1 + 5); /* impulses */
		    xs[x].y2 = yz * (MAX_QUAL+1-v + 5);
		    ox = x;
		}
	    }

	    XDrawSegments(display, qp->pm, qp->discrepGC, xs, qp->width);
	}
    }

    if (xp) free(xp);
    if (xs) free(xs);
    free(cons);
}


/* --------------------------------------------------------------------------
 * Canvas interface implementation
 */
static int qplot_create(Tcl_Interp *interp,
			Tk_Canvas canvas,
			Tk_Item *itemPtr,
			int objc,
			Tcl_Obj *CONST *objv) {
    QualityPlot *qp = (QualityPlot *)itemPtr;
    
    /* Usage is "create depth x y ?-options?". */
    if (objc < 2) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
		Tk_PathName(Tk_CanvasTkwin(canvas)), " create ",
		itemPtr->typePtr->name, " x y ?options?\"",
		(char *) NULL);
	return TCL_ERROR;
    }
    
    /* initialise item record */
    qp->qualGC = NULL;
    qp->hetGC = NULL;
    qp->discrepGC = NULL;
    qp->copyGC = NULL;
    qp->qcolour = NULL;
    qp->hcolour = NULL;
    qp->dcolour = NULL;
    qp->pm = None;
    qp->anchor = TK_ANCHOR_NW;
    qp->an_x = 0.0;
    qp->an_y = 0.0;
    qp->wy0 = qp->y_start = 0;
    qp->wy1 = qp->y_end = Tk_Height(Tk_CanvasTkwin(canvas)); /* initial world height */
    qp->width = -1;
    qp->height = -1;
    qp->yz = 1;
    qp->discrep = 0;
    qp->quality = 0;
    qp->hetero  = 0;
    
    if ((qplot_coords(interp, canvas, itemPtr, 2, objv) != TCL_OK))
	goto error;

    if (qplot_configure(interp, canvas, itemPtr,
			objc - 2, objv + 2, 0) == TCL_OK)
	return TCL_OK;
    
    /* if we get here then something is wrong */
 error:
    qplot_delete(canvas, itemPtr, Tk_Display(Tk_CanvasTkwin(canvas)));
    
    return TCL_ERROR;
}


/* invoked by the coords command */
static int qplot_coords(Tcl_Interp *interp, Tk_Canvas canvas,		
			Tk_Item *itemPtr, int objc, Tcl_Obj *CONST objv[]) {
    QualityPlot *qp = (QualityPlot *)itemPtr;
    
    if (objc == 0) {
	Tcl_Obj *obj = Tcl_NewObj();

	Tcl_Obj *subobj = Tcl_NewDoubleObj(qp->an_x);
    	Tcl_ListObjAppendElement(interp, obj, subobj);

	subobj = Tcl_NewDoubleObj(qp->an_y);
	Tcl_ListObjAppendElement(interp, obj, subobj);

	Tcl_SetObjResult(interp, obj);
    } else if (objc == 2) {
	if ((Tk_CanvasGetCoordFromObj(interp, canvas, objv[0],
		&qp->an_x) != TCL_OK)
	    || (Tk_CanvasGetCoordFromObj(interp, canvas, objv[1],
		&qp->an_y) != TCL_OK)) {
	    return TCL_ERROR;
	}
	
	update_bbox(canvas, qp);
    } else {
	char buf[64 + TCL_INTEGER_SPACE];

	sprintf(buf, "wrong # coordinates: expected 0 or 2, got %d", objc);
	Tcl_SetResult(interp, buf, TCL_VOLATILE);
	return TCL_ERROR;
    }
    
    return TCL_OK;
}


/* canvas itemconfigure callback */
static int qplot_configure(Tcl_Interp *interp,
			   Tk_Canvas canvas,
			   Tk_Item *itemPtr,
			   int objc,
			   Tcl_Obj *CONST *objv,
			   int flags) {
	    
    QualityPlot *qp = (QualityPlot *)itemPtr;
    Tk_Window tkwin;
    Display *display;
    int width, height;

    tkwin   = Tk_CanvasTkwin(canvas);
    display = Tk_Display(tkwin);
    
    if (Tk_ConfigureWidget(interp, tkwin, qplot_config_specs,
			   objc, (char **) objv,
			   (char *)qp, flags|TK_CONFIG_OBJS) != TCL_OK) {

	printf("ERROR %s\n", Tcl_GetStringResult(interp));
	return TCL_ERROR;
    }
    
    width  = Tk_Width(tkwin);
    height = Tk_Height(tkwin);
    
    if (qp->qualGC == None) {
    	XGCValues gcValues;
    	gcValues.foreground = qp->qcolour->pixel;
	qp->qualGC = Tk_GetGC(tkwin, GCForeground, &gcValues);
    }

    if (qp->hetGC == None) {
    	XGCValues gcValues;
    	gcValues.foreground = qp->hcolour->pixel;
	qp->hetGC = Tk_GetGC(tkwin, GCForeground, &gcValues);
    }

    if (qp->discrepGC == None) {
    	XGCValues gcValues;
    	gcValues.foreground = qp->dcolour->pixel;
	qp->discrepGC = Tk_GetGC(tkwin, GCForeground, &gcValues);
    }

    if (qp->copyGC == None) {
    	XGCValues gcValues;
	gcValues.function = GXcopy;
    	gcValues.foreground = 0;
	qp->copyGC = Tk_GetGC(tkwin, GCForeground, &gcValues);
    }

    if (width != qp->width || height != qp->height) {
    	/* going to need a new pixmap */
 	if (qp->pm) {
	    Tk_FreePixmap(display, qp->pm);
	}

    	qp->width = width;
    	qp->height = height;	

	if (Tk_WindowId(tkwin)) {
	    qp->pm = Tk_GetPixmap(display, Tk_WindowId(tkwin),
				  qp->width, qp->height,
				  DefaultDepthOfScreen(Tk_Screen(tkwin)));
	}
    }
    

    qplot_redraw(qp, display);
    update_bbox(canvas, qp);
    return TCL_OK;
}


/* free resources nicely*/

static void qplot_delete(Tk_Canvas canvas, Tk_Item *itemPtr, Display *display) {
    QualityPlot *qp = (QualityPlot *)itemPtr;
    
    if (qp->pm != None) {
    	Tk_FreePixmap(display, qp->pm);
    }
    
    if (qp->qcolour != NULL) {
	Tk_FreeColor(qp->qcolour);
    }

    if (qp->hcolour != NULL) {
	Tk_FreeColor(qp->hcolour);
    }

    if (qp->dcolour != NULL) {
	Tk_FreeColor(qp->dcolour);
    }
}


/* Display method: just copies from the pixmap to the canvas */
static void qplot_display(Tk_Canvas canvas, Tk_Item *itemPtr,
			  Display *display, Drawable drawable,
			  int x, int y, int width, int height) {
    
    QualityPlot *qp = (QualityPlot *)itemPtr;
    int pix_x, pix_y, pix_w, pix_h;
    short draw_x, draw_y;

    if (qp->pm != None) {
	if (x > qp->header.x1) {
	    pix_x = x - qp->header.x1;
	    pix_w = qp->header.x2 - x;
	} else {
	    pix_x = 0;
	    if ((x + width) < qp->header.x2) {
		pix_w = x + width - qp->header.x1;
	    } else {
		pix_w = qp->header.x2 - qp->header.x1;
	    }
	}
	if (y > qp->header.y1) {
	    pix_y = y - qp->header.y1;
	    pix_h = qp->header.y2 - y;
	} else {
	    pix_y = 0;
	    if ((y + height) < qp->header.y2) {
		pix_h = y + height - qp->header.y1;
	    } else {
		pix_h = qp->header.y2 - qp->header.y1;
	    }
	}
	
	Tk_CanvasDrawableCoords(canvas,
				(double) (qp->header.x1 + pix_x),
				(double) (qp->header.y1 + pix_y),
				&draw_x, &draw_y);

	/*
	 * Must modify the mask origin within the graphics context to line up
	 * with the bitmap's origin (in order to make bitmaps with
	 * "-background {}" work right).
	 */
	XSetClipOrigin(display, qp->copyGC, draw_x - pix_x, draw_y - pix_y);
	XCopyArea(display, qp->pm, drawable, qp->copyGC, pix_x, pix_y,
		  (unsigned int) pix_w, (unsigned int) pix_h,
		  draw_x, draw_y); 
	
	XSetClipOrigin(display, qp->copyGC, 0, 0);
    }
}


/* Returns proximity of itemPtr to coord */
/* Here: subverted to give world coords instead */
static double qplot_point(Tk_Canvas canvas, Tk_Item *itemPtr, double *coordPtr) {
    QualityPlot *qp = (QualityPlot *)itemPtr;
    double ax, bx;

    qp->px = coordPtr[0];
    qp->py = coordPtr[1];

    qp->py = (qp->height - qp->py) / qp->yz - 5;
    
    ax = qp->width / (qp->wx1 - qp->wx0);
    bx = qp->wx0;
    
    qp->px = (qp->px / ax) + bx;

    return 0;
}


/*
 * determines whether an item lies entirely
 * inside, entirely outside, or overlapping a given rectangle
 */
static int qplot_area(Tk_Canvas canvas, Tk_Item *itemPtr, double *rectPtr) {
    QualityPlot *qp = (QualityPlot *)itemPtr;

    if ((rectPtr[2] <= qp->header.x1)
	    || (rectPtr[0] >= qp->header.x2)
	    || (rectPtr[3] <= qp->header.y1)
	    || (rectPtr[1] >= qp->header.y2)) {
	return -1;
    }
    if ((rectPtr[0] <= qp->header.x1)
	    && (rectPtr[1] <= qp->header.y1)
	    && (rectPtr[2] >= qp->header.x2)
	    && (rectPtr[3] >= qp->header.y2)) {
	return 1;
    }
    
    return 0;
}


/* unimplemented scaling function */
static void qplot_scale(Tk_Canvas canvas, Tk_Item *itemPtr,
			double originX, double originY,
			double scaleX, double scaleY) {
    printf("QualityPlot scale - not implemented.\n");
}


/* move item by given amount */
static void qplot_translate(Tk_Canvas canvas, Tk_Item *itemPtr, 
			    double deltaX, double deltaY) {

    QualityPlot *qp = (QualityPlot *)itemPtr;

    qp->an_x += deltaX;
    qp->an_y += deltaY;
    
    update_bbox(canvas, qp);
}


/* unimplemented function to generate postscript data for printing */
static int qplot_postscript(Tcl_Interp *interp, Tk_Canvas canvas,
			    Tk_Item *itemPtr, int prepass) {
    printf("QualityPlot postscript - not implemented\n");
    
    return TCL_ERROR;
}

