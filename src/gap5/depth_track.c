#include <string.h>
#include <xalloc.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <float.h>
#include <assert.h>

#include "depth_track.h"
#include "gap_cli_arg.h"
#include "tkRaster.h"

static void dtrack_move_xhair(depth_track_t *t, int x, int y,
			     double *rx, double *ry);

/* From tclIntDecls.h */
extern Tcl_Command Tcl_GetCommandFromObj(Tcl_Interp * interp, Tcl_Obj * objPtr);

/* ------------------------------------------------------------------------ */
/* Tcl_Obj "template_display" type implementation */
/* Not sure if needed for depth_track, trying anyway */

static void dtrack_update_string(Tcl_Obj *obj);
static int  dtrack_from_any(Tcl_Interp *interp, Tcl_Obj *obj);

static Tcl_ObjType dtrack_obj_type = {
    "depth_track",
    (Tcl_FreeInternalRepProc*)NULL,
    (Tcl_DupInternalRepProc*)NULL,
    dtrack_update_string,
    dtrack_from_any
};

static void dtrack_update_string(Tcl_Obj *obj) {
    depth_track_t *t = obj->internalRep.otherValuePtr;
    obj->bytes = ckalloc(30);
    obj->length = sprintf(obj->bytes, "dtrack=%p", t);
}

static int dtrack_from_any(Tcl_Interp *interp, Tcl_Obj *obj) {
    char *bytes;
    int length;
    depth_track_t *t;

    if (NULL == (bytes = Tcl_GetStringFromObj(obj, &length)))
	return TCL_ERROR;

    if (0 != strncmp(bytes, "dtrack=", 3))
	return TCL_ERROR;

    /* Free the old internalRep before setting the new one. */
    if (obj->typePtr && obj->typePtr->freeIntRepProc)
	(*obj->typePtr->freeIntRepProc)(obj);

    /* Convert the hex value to a pointer once more */
    if (1 != sscanf(bytes+3, "%p", &t))
	return TCL_ERROR;

    obj->internalRep.otherValuePtr = t;
    obj->typePtr = &dtrack_obj_type;
    return TCL_OK;
}

/* not all of these needed, prune list later */
static Tk_OptionSpec optionSpecs[] = {
    {TK_CONFIG_INT, "-logy", "logy", "LogY", "0",
     -1, Tk_Offset(depth_track_t, logy), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-cmode", "colourMode", "ColourMode", "0",
     -1, Tk_Offset(depth_track_t, cmode), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-ymode", "YMode", "YMode", "0",
     -1, Tk_Offset(depth_track_t, ymode), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-yoffset", "YOffset", "YOffset", "0",
     -1, Tk_Offset(depth_track_t, yoffset), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-accuracy", "accuracy", "Accuracy", "0",
     -1, Tk_Offset(depth_track_t, accuracy), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-spread", "spread", "Spread", "0",
     -1, Tk_Offset(depth_track_t, spread), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-reads_only", "readsOnly", "ReadsOnly", "0",
     -1, Tk_Offset(depth_track_t, reads_only), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-by_strand", "byStrand", "ByStrand", "1",
     -1, Tk_Offset(depth_track_t, sep_by_strand), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-filter", "filter", "Filter", "0",
     -1, Tk_Offset(depth_track_t, filter), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-plot_depth", "plotDepth", "PlotDepth", "0",
     -1, Tk_Offset(depth_track_t, plot_depth), 0, 0, 0 /* mask */},
    {TK_CONFIG_DOUBLE, "-xzoom", "xZoom", "XZoom", "10.0",
     -1, Tk_Offset(depth_track_t, xzoom), 0, 0, 0 /* mask */},
    {TK_CONFIG_DOUBLE, "-yzoom", "yZoom", "YZoom", "10.0",
     -1, Tk_Offset(depth_track_t, yzoom), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-min_qual", "minQual", "MinQual", "0",
     -1, Tk_Offset(depth_track_t, min_qual), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-max_qual", "maxQual", "MaxQual", "255",
     -1, Tk_Offset(depth_track_t, max_qual), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-min_y_size", "minYSize", "MinYSize", "512",
     -1, Tk_Offset(depth_track_t, min_sz), 0, 0, 0 /* mask */},
    {TK_OPTION_END}
};

static int dtrack_cmd(ClientData clientData, Tcl_Interp *interp,
		      int objc, Tcl_Obj *CONST objv[]) {
    int index, replot = 0;
    depth_track_t *t = (depth_track_t *)clientData;
    Tcl_Obj *res = NULL;

    static char *options[] = {
	"delete",       "cget",    "configure", "replot",
	"ymin",         "ymax",         "yrange",  "xhair", "yz",
	(char *)NULL,
    };

    enum options {
	DELETE,         CGET,      CONFIGURE,   REPLOT,
	YMIN,           YMAX,      YRANGE,      XHAIR,    YZ
    };

    if (objc < 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "option arg ?arg ...?");
	return TCL_ERROR;
    }

    if (Tcl_GetIndexFromObj(interp, objv[1], options, "option", 0,
            &index) != TCL_OK) {
        return TCL_ERROR;
    }

    switch ((enum options)index) {
    case DELETE:
	Tcl_DeleteCommandFromToken(interp,
				   Tcl_GetCommandFromObj(interp, objv[0]));
	break;

    case CGET:
	if (objc != 3) {
	    Tcl_WrongNumArgs(interp, 2, objv, "cget");
	    return TCL_ERROR;
	}
	break;

    case CONFIGURE:
	if (objc == 2) {
	    res = Tk_GetOptionValue(interp, (char *)t, t->optionTable,
				    NULL, t->tkwin);
	    if (!res)
		return TCL_ERROR;
	} else if (objc == 3) {
	    res = Tk_GetOptionValue(interp, (char *)t, t->optionTable,
				    objv[2], t->tkwin);
	    if (!res)
		return TCL_ERROR;
	} else {
	    if (Tk_SetOptions(interp, (char *)t, t->optionTable,
			      objc-2, objv+2, t->tkwin, NULL, NULL)
		!= TCL_OK)
		return TCL_ERROR;
	    replot = 1;
	}
	if (res)
	    Tcl_SetObjResult(interp, res);
	break;

    case REPLOT:
	replot = 1;
	break;

    case YMIN:
	Tcl_SetObjResult(interp, Tcl_NewIntObj(t->ymin));
	break;

    case YMAX:
	Tcl_SetObjResult(interp, Tcl_NewIntObj(t->ymax));
	break;
	
    case YZ:
	Tcl_SetObjResult(interp, Tcl_NewDoubleObj(t->yz));
	break;
    	

    case YRANGE: {
	char buf[1024];
	double top, bottom;

    	/* Line adjusts to raster size, so no scolling needed for the moment.
	   There may be value in having an adustable size, but not yet.
	*/
	
	bottom = 1;
	top = 0;

	sprintf(buf, "%f %f", top, bottom);
	Tcl_SetObjResult(interp, Tcl_NewStringObj(buf, -1));
	break;
    }

    case XHAIR:
	if (objc == 4) {
	    int x, y;
	    double rx, ry;
	    Tcl_Obj *obj[2];
	    int width, height, height2;
	    
	    RasterWinSize(t->raster, &width, &height);
	    height2 = height / 2;

	    Tcl_GetIntFromObj(interp, objv[2], &x);
	    Tcl_GetIntFromObj(interp, objv[3], &y);
	    dtrack_move_xhair(t, x, y, &rx, &ry);

    	    ry  = (height - ry) / t->yz;
	    
	    obj[0] = Tcl_NewDoubleObj(rx);
	    obj[1] = Tcl_NewDoubleObj(ry);
	    Tcl_SetObjResult(interp, Tcl_NewListObj(2, obj));
	} else {
	    dtrack_move_xhair(t, INT_MAX, INT_MAX, NULL, NULL);
	}
	break;
    }

    if (replot)
	depth_track_replot(t);

    return TCL_OK;
}


static void dtrack_cmd_delete(ClientData clientData) {
    depth_track_destroy((depth_track_t *)clientData);
}

static void dt_event_proc(ClientData clientData, XEvent *eventPtr) {
    if (eventPtr->type == DestroyNotify) {
    	depth_track_destroy((depth_track_t *)clientData);
    }
}


typedef struct {
    char *raster_win;
    char *range_obj;
} dtrack_arg;

static int tcl_depth_track(ClientData clientData, Tcl_Interp *interp,
				int objc, Tcl_Obj *CONST objv[]) {
    depth_track_t *dt;
    Tcl_Obj *tobj;
    Tk_Raster *raster;
    Tcl_CmdInfo info;

    dtrack_arg args;
    cli_args a[] = {
	{"-raster", ARG_STR, 1, NULL, offsetof(dtrack_arg, raster_win)},
	{"-range",  ARG_STR, 1, NULL, offsetof(dtrack_arg, range_obj)},
	{NULL,	     0,	      0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;
	
	
    if (!Tcl_GetCommandInfo(interp, args.raster_win, &info))
	return TCL_ERROR;

    raster = (Tk_Raster*)info.clientData;

    if (NULL == (dt = depth_track_new(interp, raster)))
	return TCL_ERROR;

    if (!Tcl_GetCommandInfo(interp, args.range_obj, &info))
	return TCL_ERROR;

    dt->gr = (gap_range_t*)info.objClientData; // N.B. different from rasta
    gap_range_test(dt->gr);
	
    if (NULL == (tobj = Tcl_NewObj())) {
	free(dt);
	return TCL_ERROR;
    }
    
    tobj->internalRep.otherValuePtr = (VOID *)dt;
    tobj->typePtr = &dtrack_obj_type;
    dtrack_update_string(tobj);

    dt->optionTable = Tk_CreateOptionTable(interp, optionSpecs);
    if (Tk_InitOptions(interp, (char *)dt, dt->optionTable, dt->tkwin)
	!= TCL_OK) {
	free(dt);
	return TCL_ERROR;
    }

    /* Register the string form as a new command */
    if (NULL == Tcl_CreateObjCommand(interp, tobj->bytes, dtrack_cmd,
				     (ClientData)dt,
				     (Tcl_CmdDeleteProc *)dtrack_cmd_delete))
	return TCL_ERROR;
	
    Tk_CreateEventHandler(dt->tkwin,
    	    	    	     StructureNotifyMask,
			     dt_event_proc, (ClientData) dt);

    Tcl_SetObjResult(interp, tobj);
    
    return TCL_OK;
}

/* Our only external function */
int DTrack_Init(Tcl_Interp *interp) {
    Tcl_RegisterObjType(&dtrack_obj_type);

    if (NULL == Tcl_CreateObjCommand(interp, "g5::depth_track",
				     tcl_depth_track,
				     (ClientData)NULL,
				     (Tcl_CmdDeleteProc *)NULL))
	return TCL_ERROR;

    return TCL_OK;
}


/* ------------------------------------------------------------------------ */
/* And the C to actually implement it, minus the Tcl interface gubbins above */
depth_track_t *depth_track_new(Tcl_Interp *interp,
			      Tk_Raster *raster) {
    depth_track_t *t = (depth_track_t *)calloc(1, sizeof(*t));
    int i;
    char *opts[10], b[1024];

    if (!t)
	return NULL;

    t->interp = interp;
    t->raster = raster;
    t->tkwin = GetRasterTkWin(raster);

    opts[0] = "-fg";
    opts[1] = b;
    opts[2] = "-linewidth";
    opts[3] = "0";
    opts[4] = "-function";
    opts[5] = "copy";
    opts[6] = NULL;

    for (i = 0; i < 32; i++) {
	sprintf(b, "#%02x%02x%02x", 64+i*5, 64+i*5, 64+i*5);
	t->map_col[i] = CreateDrawEnviron(interp, raster, 6, opts);
	SetDrawEnviron(t->interp, t->raster, t->map_col[i]);
    }
   
    opts[1] = "orange";
    t->span_col = CreateDrawEnviron(interp, raster, 6, opts);

    opts[1] = "blue";
    t->single_col = CreateDrawEnviron(interp, raster, 6, opts);

    opts[1] = "red";
    t->inconsistent_col = CreateDrawEnviron(interp, raster, 6, opts);

    opts[1] = "green4";
    t->fwd_col  = CreateDrawEnviron(interp, raster, 6, opts);
    opts[3] = "3";
    t->fwd_col3 = CreateDrawEnviron(interp, raster, 6, opts);
    opts[3] = "0";

    opts[1] = "magenta";
    t->rev_col  = CreateDrawEnviron(interp, raster, 6, opts);
    opts[3] = "3";
    t->rev_col3 = CreateDrawEnviron(interp, raster, 6, opts);
    opts[3] = "0";
   
    opts[1] = "green";
    opts[5] = "xor";
    t->xhair_col = CreateDrawEnviron(interp, raster, 6, opts);

    t->tdepth = NULL;
    t->sdepth = NULL;
    t->depth_width = 0;
    t->xhair_pos = DBL_MAX;
    t->yhair_pos = DBL_MAX;
    t->wx0 = 0;
    t->wx1 = 0;
    
    return t;
}

void depth_track_destroy(depth_track_t *t) {
    if (!t)
	return;

    if (t->tdepth)
	free(t->tdepth);
    if (t->sdepth)
	free(t->sdepth);
	
    /* Draw Environments automatically freed on raster destroy */

    free(t);
}


int depth_track_replot(depth_track_t *t) {
    double wx0, wy0, wx1, wy1, ax, ay, bx, by;
    int j;
    double tsize = 1000;
    int ymin = INT_MAX;
    int ymax = INT_MIN;
    int width, height, height2;
    Display *rdisp;
    Drawable rdraw;
    GC rgc;
    int force_change = 0;
    XPoint *sp, *tp;
    
    if (!t)
	return -1;

    tk_RasterClear(t->raster);
    RasterWinSize(t->raster, &width, &height);
    GetRasterCoords(t->raster, &wx0, &wy0, &wx1, &wy1);
//    printf("DT Height %f to %f height %d\n", wy0, wy1, height);
 //   printf("DT Width  %f to %f width  %d\n", wx0, wx1, width);
    
    wx0 -= tsize; /* we use some of the reads beyond the window size */
    wx1 += tsize;

    /* if range has changed (or no range) recalculate */
    
    if (gap_range_recalculate(t->gr, width, wx0, wx1, t->gr->template_mode, force_change)) {
    
	if (t->gr->r == NULL) {
	    return 0;
	}
	
	force_change = 1;
	
    }
    

    /*
     * Do this in multiple passes.
     * 1) Compute X (start/end of lines to draw), colours, status.
     * 2) Plot (possibly combined with 2)
     */

    GetWorldToRasterConversion(t->raster, &ax, &ay, &bx, &by);

    t->ntl = gap_range_x(t->gr, ax, bx, t->fwd_col, t->rev_col, 
    	    	    	    t->single_col, t->span_col, t->inconsistent_col,
		    	    force_change, t->reads_only);


    /* 2) Plot the lines */

    /* prepare image */
    height2 = height / 2;
    t->yzoom = 100;
    t->yz = t->yzoom / 50.0;

    if (t->gr->max_height > 0) {    
    	t->yz = ((double)height / (double)t->gr->max_height) * 0.95;
    } 
    
    
    /* Now to draw the show the image */
    rdisp = GetRasterDisplay (t->raster);
    rdraw = GetRasterDrawable(t->raster);
    rgc   = GetRasterGC      (t->raster);


    /* Plot depth */
    sp = (XPoint *)calloc(width, sizeof(*sp));
    tp = (XPoint *)calloc(width, sizeof(*tp));

    for (j = 0; j < width; j++) {
	sp[j].x = j;
	sp[j].y = height - 5 - t->gr->depth[j].s * t->yz;
	tp[j].x = j;
	tp[j].y = height - 5 - t->gr->depth[j].t * t->yz;
	
	if (ymin > sp[j].y) ymin = sp[j].y;
	if (ymin > tp[j].y) ymin = tp[j].y;
	
	if (ymax < sp[j].y) ymax = sp[j].y;
	if (ymax < tp[j].y) ymax = tp[j].y;
    } 

    if (!(t->filter & FILTER_PAIRED)) {
	SetDrawEnviron(t->interp, t->raster, t->fwd_col);
	XDrawLines(rdisp, rdraw, GetRasterGC(t->raster), sp, width, CoordModeOrigin);
    }
    
    SetDrawEnviron(t->interp, t->raster, t->rev_col);
    XDrawLines(rdisp, rdraw, GetRasterGC(t->raster), tp, width, CoordModeOrigin);

    free(sp);
    free(tp);


    /* Compute scrollbar and bounding box locations */
    t->ymin = ymin != INT_MAX ? ymin : 0;
    t->ymax = ymax != INT_MIN ? ymax : 0;
    
    t->xhair_pos = DBL_MAX;
    t->yhair_pos = DBL_MAX;

    return 0;
}

static void dtrack_move_xhair(depth_track_t *t, int x, int y,
			     double *rx, double *ry) {
    double wx0, wx1, wy0, wy1;
    double xh, xh2, yh, yh2;
    int width, height;

    RasterWinSize(t->raster, &width, &height);
    RasterToWorld(t->raster, x, y, &xh, &yh);
    RasterToWorld(t->raster, x+1, y+1, &xh2, &yh2);
    xh = (xh+xh2)/2;
    yh = (yh+yh2)/2;
    RasterToWorld(t->raster, 0, 0, &wx0, &wy0);
    RasterToWorld(t->raster, width, height, &wx1, &wy1);

    SetDrawEnviron(t->interp, t->raster, t->xhair_col);
    if (t->xhair_pos != DBL_MAX) {
	RasterDrawLine(t->raster, t->xhair_pos, wy0, t->xhair_pos, wy1);
    }
    if (t->yhair_pos != DBL_MAX) {
	RasterDrawLine(t->raster, wx0, t->yhair_pos, wx1, t->yhair_pos);
    }

    if (x != INT_MAX) {
	t->xhair_pos = xh;
	RasterDrawLine(t->raster, t->xhair_pos, wy0, t->xhair_pos, wy1);
	if (rx) *rx = t->xhair_pos;
    } else {
	t->xhair_pos = DBL_MAX;
    }

    if (y != INT_MAX) {
	t->yhair_pos = yh;
	RasterDrawLine(t->raster, wx0, t->yhair_pos, wx1, t->yhair_pos);
	if (ry) *ry = t->yhair_pos;
    } else {
	t->yhair_pos = DBL_MAX;
    }

    tk_RasterRefresh(t->raster);

}
