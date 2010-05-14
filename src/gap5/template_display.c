#include <string.h>
#include <xalloc.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <float.h>
#include <assert.h>

#include "template_display.h"
#include "gap_cli_arg.h"
#include "tree.h"
#include "tkRaster.h"

static void tdisp_move_xhair(template_disp_t *t, int x, int y,
			     double *rx, double *ry);

/* From tclIntDecls.h */
extern Tcl_Command Tcl_GetCommandFromObj(Tcl_Interp * interp, Tcl_Obj * objPtr);

/* ------------------------------------------------------------------------ */
/* Tcl_Obj "template_display" type implementation */

static void tdisp_update_string(Tcl_Obj *obj);
static int tdisp_from_any(Tcl_Interp *interp, Tcl_Obj *obj);

static Tcl_ObjType tdisp_obj_type = {
    "template_display",
    (Tcl_FreeInternalRepProc*)NULL,
    (Tcl_DupInternalRepProc*)NULL,
    tdisp_update_string,
    tdisp_from_any
};

static void tdisp_update_string(Tcl_Obj *obj) {
    template_disp_t *t = obj->internalRep.otherValuePtr;
    obj->bytes = ckalloc(30);
    obj->length = sprintf(obj->bytes, "tdisp=%p", t);
}

static int tdisp_from_any(Tcl_Interp *interp, Tcl_Obj *obj) {
    char *bytes;
    int length;
    template_disp_t *t;

    if (NULL == (bytes = Tcl_GetStringFromObj(obj, &length)))
	return TCL_ERROR;

    if (0 != strncmp(bytes, "tdisp=", 3))
	return TCL_ERROR;

    /* Free the old internalRep before setting the new one. */
    if (obj->typePtr && obj->typePtr->freeIntRepProc)
	(*obj->typePtr->freeIntRepProc)(obj);

    /* Convert the hex value to a pointer once more */
    if (1 != sscanf(bytes+3, "%p", &t))
	return TCL_ERROR;

    obj->internalRep.otherValuePtr = t;
    obj->typePtr = &tdisp_obj_type;
    return TCL_OK;
}

static Tk_OptionSpec optionSpecs[] = {
    {TK_CONFIG_INT, "-logy", "logy", "LogY", "0",
     -1, Tk_Offset(template_disp_t, logy), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-cmode", "colourMode", "ColourMode", "0",
     -1, Tk_Offset(template_disp_t, cmode), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-ymode", "YMode", "YMode", "0",
     -1, Tk_Offset(template_disp_t, ymode), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-yoffset", "YOffset", "YOffset", "0",
     -1, Tk_Offset(template_disp_t, yoffset), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-accuracy", "accuracy", "Accuracy", "0",
     -1, Tk_Offset(template_disp_t, accuracy), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-spread", "spread", "Spread", "0",
     -1, Tk_Offset(template_disp_t, spread), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-reads_only", "readsOnly", "ReadsOnly", "0",
     -1, Tk_Offset(template_disp_t, reads_only), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-by_strand", "byStrand", "ByStrand", "1",
     -1, Tk_Offset(template_disp_t, sep_by_strand), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-filter", "filter", "Filter", "0",
     -1, Tk_Offset(template_disp_t, filter), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-plot_depth", "plotDepth", "PlotDepth", "0",
     -1, Tk_Offset(template_disp_t, plot_depth), 0, 0, 0 /* mask */},
    {TK_CONFIG_DOUBLE, "-xzoom", "xZoom", "XZoom", "10.0",
     -1, Tk_Offset(template_disp_t, xzoom), 0, 0, 0 /* mask */},
    {TK_CONFIG_DOUBLE, "-yzoom", "yZoom", "YZoom", "10.0",
     -1, Tk_Offset(template_disp_t, yzoom), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-min_qual", "minQual", "MinQual", "0",
     -1, Tk_Offset(template_disp_t, min_qual), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-max_qual", "maxQual", "MaxQual", "255",
     -1, Tk_Offset(template_disp_t, max_qual), 0, 0, 0 /* mask */},
    {TK_CONFIG_INT, "-min_y_size", "minYSize", "MinYSize", "512",
     -1, Tk_Offset(template_disp_t, min_sz), 0, 0, 0 /* mask */},
    {TK_OPTION_END}
};

static int tdisp_cmd(ClientData clientData, Tcl_Interp *interp,
		      int objc, Tcl_Obj *CONST objv[]) {
    int index, replot = 0;
    template_disp_t *t = (template_disp_t *)clientData;
    Tcl_Obj *res = NULL;

    static char *options[] = {
	"delete",       "io",           "cget",    "configure", "replot",
	"ymin",         "ymax",         "yrange",  "xhair", "yz",
	(char *)NULL,
    };

    enum options {
	DELETE,         IO,     	CGET,      CONFIGURE,   REPLOT,
	YMIN,           YMAX,           YRANGE,    XHAIR,   	YZ
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

    case IO:
	Tcl_SetResult(interp, io_obj_as_string(t->io) , TCL_VOLATILE);
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
	double wx0, wy0, wx1, wy1;
//    	int width, height;
	double top, bottom;
	double wx_start, wx_end, wy_start, wy_end;

	GetRasterCoords(t->raster, &wx0, &wy0, &wx1, &wy1);
    	RasterGetWorldScroll(t->raster, &wx_start, &wy_start, &wx_end, &wy_end); 
//   	RasterWinSize(t->raster, &width, &height);
	
/*
    	printf("TD Raster wx0 %f wx1 %f wy0 %f wy1 %f\n", wx0, wx1, wy0, wy1);
	printf("TD Y min %d Y max %d\n", t->ymin, t->ymax);
    	printf("TD wx_start %f wx_end %f wy_start %f wy_end %f\n", wx_start, wx_end, wy_start, wy_end);
	printf("TD width %d height %d\n", width, height);
*/
    	
	if (t->ymax == 0 || wy1 == 0) {
	    top = 0;
	    bottom = 1;
	} else {
	    top = (wy0 - wy_start) / (t->ymax - wy_start);
	    bottom = (wy1 - wy_start) / (t->ymax - wy_start);
	    
	    if (top < 0) {
	    	top = 0;
	    }
	    
	    if (bottom > 1) {
	    	bottom = 1;
	    }
	}
    	
//	printf("TD top %f bottom %f\n", top, bottom);
	
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
	    tdisp_move_xhair(t, x, y, &rx, &ry);
	    
	    if (t->sep_by_strand) {
		if (ry > height2)
		    ry = ry-height2;
		else
		    ry = height2-ry;
	    }

	    if (t->ymode == 1)
		ry /= 10;

	    ry = ry/(t->yzoom / 200.0) + t->yoffset - 50;
	    if (t->logy && t->ymode != 1) {
		ry = exp(ry/50.0)-1;
	    }
	    ry++;
	    
	    obj[0] = Tcl_NewDoubleObj(rx);
	    obj[1] = Tcl_NewDoubleObj(ry);
	    Tcl_SetObjResult(interp, Tcl_NewListObj(2, obj));
	} else {
	    tdisp_move_xhair(t, INT_MAX, INT_MAX, NULL, NULL);
	}
	break;
    }

    if (replot)
	template_replot(t);

    return TCL_OK;
}


static void tdisp_cmd_delete(ClientData clientData) {
    template_destroy((template_disp_t *)clientData);
}

typedef struct {
    GapIO *io;
    int cnum;
    char *raster_win;
    char *range_obj;
} tdisp_arg;


static void td_event_proc(ClientData clientData, XEvent *eventPtr) {
    template_disp_t *td = (template_disp_t *)clientData;
    
    if (eventPtr->type == DestroyNotify) {
    	template_destroy(td);
    }
}
    
    
static int tcl_template_display(ClientData clientData, Tcl_Interp *interp,
				int objc, Tcl_Obj *CONST objv[]) {
    template_disp_t *td;
    Tcl_Obj *tobj;
    Tk_Raster *raster;
    Tcl_CmdInfo info;

    tdisp_arg args;
    cli_args a[] = {
	{"-io",     ARG_IO,  1, NULL, offsetof(tdisp_arg, io)},
	{"-cnum",   ARG_INT, 1, NULL, offsetof(tdisp_arg, cnum)},
	{"-raster", ARG_STR, 1, NULL, offsetof(tdisp_arg, raster_win)},
	{"-range",  ARG_STR, 1, NULL, offsetof(tdisp_arg, range_obj)},
	{NULL,	     0,	      0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;
	
	
    if (!Tcl_GetCommandInfo(interp, args.raster_win, &info))
	return TCL_ERROR;

    raster = (Tk_Raster*)info.clientData;

    if (NULL == (td = template_new(args.io, args.cnum, interp, raster)))
	return TCL_ERROR;

    if (!Tcl_GetCommandInfo(interp, args.range_obj, &info))
	return TCL_ERROR;

    td->gr = (gap_range_t*)info.objClientData; // N.B. different from rasta
    gap_range_test(td->gr);
	
    if (NULL == (tobj = Tcl_NewObj())) {
	free(td);
	return TCL_ERROR;
    }
    
    tobj->internalRep.otherValuePtr = (VOID *)td;
    tobj->typePtr = &tdisp_obj_type;
    tdisp_update_string(tobj);

    td->optionTable = Tk_CreateOptionTable(interp, optionSpecs);
    if (Tk_InitOptions(interp, (char *)td, td->optionTable, td->tkwin)
	!= TCL_OK) {
	free(td);
	return TCL_ERROR;
    }

    /* Register the string form as a new command */
    if (NULL == Tcl_CreateObjCommand(interp, tobj->bytes, tdisp_cmd,
				     (ClientData)td,
				     (Tcl_CmdDeleteProc *)tdisp_cmd_delete))
	return TCL_ERROR;
	
    Tk_CreateEventHandler(td->tkwin,
    	    	    	     StructureNotifyMask,
			     td_event_proc, (ClientData) td);

    Tcl_SetObjResult(interp, tobj);
    
    return TCL_OK;
}

/* Our only external function */
int TDisp_Init(Tcl_Interp *interp) {
    Tcl_RegisterObjType(&tdisp_obj_type);

    if (NULL == Tcl_CreateObjCommand(interp, "g5::template_display",
				     tcl_template_display,
				     (ClientData)NULL,
				     (Tcl_CmdDeleteProc *)NULL))
	return TCL_ERROR;

    return TCL_OK;
}


/* ------------------------------------------------------------------------ */
/* And the C to actually implement it, minus the Tcl interface gubbins above */
template_disp_t *template_new(GapIO *io, int cnum,
			      Tcl_Interp *interp,
			      Tk_Raster *raster) {
    template_disp_t *t = (template_disp_t *)calloc(1, sizeof(*t));
    int i;
    char *opts[10], b[1024];

    if (!t)
	return NULL;

    t->io = io;
    t->contig = (contig_t *)cache_search(io, GT_Contig, cnum);
    t->crec = t->contig->rec;
    if (!t->contig)
	return NULL;
    cache_incr(io, t->contig);
    
    t->interp = interp;
    t->raster = raster;
    t->tkwin = GetRasterTkWin(raster);

    if (NULL == (t->image = initialise_image(GetRasterDisplay(t->raster)))) return NULL;

    opts[0] = "-fg";
    opts[1] = b;
    opts[2] = "-linewidth";
    opts[3] = "0";
    opts[4] = "-function";
    opts[5] = "copy";
    opts[6] = NULL;

    for (i = 0; i < 32; i++) {
	add_colour(t->image, 64+i*5, 64+i*5, 64+i*5);
    }

    t->span_col         = add_colour(t->image, 255, 165, 0); // orange
    t->single_col       = add_colour(t->image, 0, 0, 255);   // blue
    t->inconsistent_col = add_colour(t->image, 255, 0, 0);   // red
    t->fwd_col  = t->fwd_col3 = add_colour(t->image, 0, 139, 0);   // green4
    t->rev_col  = t->rev_col3 = add_colour(t->image, 255, 0, 255); // magenta
    t->background = add_colour(t->image, 0, 0, 0); // black
    
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

void template_destroy(template_disp_t *t) {
    if (!t)
	return;

    if (t->io && t->contig)
	cache_decr(t->io, t->contig);

    if (t->tdepth)
	free(t->tdepth);
    if (t->sdepth)
	free(t->sdepth);
	
    if (t->image) {	
	image_destroy(t->image);
    }
    
    /* Draw Environments automatically freed on raster destroy */

    free(t);
}



/* -------------------------------------------------------------------------
 * Y Coordinate allocation routines
 */

SPLAY_HEAD(XTREE, xy_pair);
SPLAY_HEAD(YTREE, xy_pair);
SPLAY_PROTOTYPE(XTREE, xy_pair, link, x_cmp);
SPLAY_PROTOTYPE(YTREE, xy_pair, link, y_cmp);
SPLAY_GENERATE(XTREE, xy_pair, link, x_cmp);
SPLAY_GENERATE(YTREE, xy_pair, link, y_cmp);

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
static int compute_ypos(template_disp_t *t, int xgap, tline *tl, int ntl) {
    int i;
    struct xy_pair *node, *curr, *next;
    int yn = 0;
    int min_sz = INT_MIN;
    int max_sz = t->min_sz;
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
#if 0
		/* Cull Y tree if appropriate to remove excess rows */
		if (try_cull) {
		    for (curr = SPLAY_MAX(YTREE, &ytree); curr; curr = next) {
			next = SPLAY_PREV(YTREE, &ytree, curr);
			if (curr->y == yn) {
			    SPLAY_REMOVE(YTREE, &ytree, curr);
			    free(curr);
			    yn--;
			} else {
			    break;
			}
		    }
		}
#endif

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
	max_sz += t->min_sz;
	if (max_sz < min_sz+10)
	    max_sz = min_sz+10;
	yoffset += ymax;
    }
    return 0;
}

int sort_by_mq(void *p1, void *p2) {
    rangec_t *r1 = (rangec_t *)p1;
    rangec_t *r2 = (rangec_t *)p2;

    return (int)((r1->mqual + r1->pair_mqual)/8)
	-  (int)((r2->mqual + r2->pair_mqual)/8);
}

int sort_by_rec(void *p1, void *p2) {
    rangec_t *r1 = (rangec_t *)p1;
    rangec_t *r2 = (rangec_t *)p2;

    return r1->rec - r2->rec;
}

int sort_tline_by_x(const void *p1, const void *p2) {
    tline *r1 = (tline *)p1;
    tline *r2 = (tline *)p2;

    return r1->x[0] - r2->x[1];
}

int template_replot(template_disp_t *t) {
    double wx0, wy0, wx1, wy1, ny0, ny1, ax, ay, bx, by;
    int i, j, mode;
    struct timeval tv1, tv2;
    double t1, t2 = 0, t3, t4;
    double tsize = 1000;
    int ymin = INT_MAX;
    int ymax = INT_MIN;
    int width, height, height2;
    static int last_zoom = 0;
    int fwd_col, rev_col;
    Display *rdisp;
    Drawable rdraw;
    GC rgc;
    int force_change = 0;
    
    if (!t)
	return -1;
    if (!t->io)
	return -1;

    tk_RasterClear(t->raster);
    RasterWinSize(t->raster, &width, &height);
    GetRasterCoords(t->raster, &wx0, &wy0, &wx1, &wy1);

    // printf("templates %f to %f\n", wx0, wx1);
/*
    printf("TD Height %f to %f height %d\n", wy0, wy1, height);
    printf("TD Width  %f to %f width  %d\n", wx0, wx1, width);
*/
    wx0 -= tsize; /* we use some of the reads beyond the window size */
    wx1 += tsize;

    /* Timing test start*/
    gettimeofday(&tv1, NULL);
    
    /* Find sequences on screen */
    mode = t->reads_only ? 0 : CSIR_PAIR;
    
    if (t->ymode == 1)
	mode |= CSIR_SORT_BY_Y;
	
    set_filter(t->gr, t->filter, t->min_qual, t->max_qual, t->cmode, t->accuracy);
    
    if (gap_range_recalculate(t->gr, width, wx0, wx1, mode, force_change)) {
    
	if (t->gr->r == NULL) {
	    return 0;
	}
	
	force_change = 1;
    }
    
    /* more timing */
    gettimeofday(&tv2, NULL);
    t1 = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;

    /*
     * Do this in multiple passes.
     * 1) Compute X (start/end of lines to draw), colours, status.
     * 2) Compute Y
     * 3) Plot (possibly combined with 2)
     */

    fwd_col = t->yzoom >= 150 ? t->fwd_col3 : t->fwd_col;
    rev_col = t->yzoom >= 150 ? t->rev_col3 : t->rev_col;

    /* 1) Compute X */
    GetWorldToRasterConversion(t->raster, &ax, &ay, &bx, &by);

    t->ntl = gap_range_x(t->gr, ax, bx, fwd_col, rev_col, 
    	    	    	    t->single_col, t->span_col, t->inconsistent_col,
		    	    force_change, t->reads_only);
			    
   /* 2) Compute Y coordinates (part 1) */
    if (t->ymode == 1) {
    	double xgap = 100;
    
	puts("Sorting");
	qsort(t->gr->tl, t->ntl, sizeof(tline), sort_tline_by_x);
	puts("computing ypos");
	compute_ypos(t, xgap, t->gr->tl, t->ntl);
	puts("done");
    }



    /* 3) Plot the lines */

    /* prepare image */
    image_remove(t->image);
    if(!create_image_buffer(t->image, width, height, t->background)) return -1;

    height2 = height / 2;
    t->yz = t->yzoom / 200.0;

    for (i = 0; i < t->ntl; i++) {
    	tline *tl = &t->gr->tl[i];
    
	for (j = 0; j < 3; j++) {
	
	    if (tl->x[j] < tl->x[j+1]) {
		double y;

		/* Compute the Y coordinates (part 2) */
		/* See XHAIR sub-command code too; keep it in sync */
		if (t->ymode == 1) {
		    y = tl->y * 10;
		} else {
		    if (t->ymode == 0)
			y = tl->x[3] - tl->x[0];
		    else
			y = tl->mq * 4;
			
		    if (t->logy) {
			if (y < 0) y = 0;
			
			y = 50 * log(y + 1);
		    }
		}
		
		y = (y + 50 - t->yoffset) * t->yz;

		if (t->spread)
		    y = y - t->spread / 2 + ((tl->x[0] + tl->rec) % (t->spread));

		if (t->sep_by_strand)
		    y = tl->t_strand ? height2 - y : height2 + y;

		if (ymin > y) ymin = y;
		if (ymax < y) ymax = y;

		/* And plot if visible */
		if (y >= wy0 && y <= wy1) {
		    int rx1, rx2, ry;
		    
		    rx1 = (tl->x[j] - bx) * ax; // world to raster conversion
		    rx2 = (tl->x[j + 1] - bx) * ax;
		    ry  = (y - by) * ay;
		    
		    draw_line(t->image, rx1, rx2, ry, tl->col[j]);
		}
		
		
	    }
	}
    }
    
    /* Now to draw the show the image */
    create_image_from_buffer(t->image);
    rdisp = GetRasterDisplay (t->raster);
    rdraw = GetRasterDrawable(t->raster);
    rgc   = GetRasterGC      (t->raster);

    XPutImage(rdisp, rdraw, rgc, t->image->img, 0, 0, 0, 0, t->image->width, t->image->height);

    /* Compute scrollbar and bounding box locations */
    t->ymin = ymin != INT_MAX ? ymin : 0;
    t->ymax = ymax != INT_MIN ? ymax : 0;

    tv1 = tv2;
    gettimeofday(&tv2, NULL);
    t4 = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;

//    printf("Query range %d..%d => %d reads, E%5.3fs + X%5.3fs + Y%5.3fs + D%5.3fs\n",
//        	   (int)wx0, (int)wx1, nr, t1, t2, t3, t4);

    ny0 = wy0; ny1 = wy1;
    if (t->yzoom != last_zoom) {
	ny0 = t->ymin - 10;
	ny1 = t->ymax + 10;
    } else {
	if (ny0 > t->ymin - 10)
	    ny0 = t->ymin - 10;
	if (ny1 < t->ymax + 10)
	    ny1 = t->ymax + 10;
    }
    last_zoom = t->yzoom;
    
    if (ny0 > 0) ny0 = 0;
    if (ny1 < 400) ny1 = 400;
    
    wx0 = contig_get_start(&t->contig) - 10;
    wx1 = contig_get_end(&t->contig) + 10;
    RasterSetWorldScroll(t->raster,  wx0,  ny0,  wx1,  ny1);

    t->xhair_pos = DBL_MAX;
    t->yhair_pos = DBL_MAX;

    return 0;
}

static void tdisp_move_xhair(template_disp_t *t, int x, int y,
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
