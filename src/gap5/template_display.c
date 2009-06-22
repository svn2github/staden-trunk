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
	"ymin",         "ymax",         "yrange",  "xhair",
	(char *)NULL,
    };

    enum options {
	DELETE,         IO,     	CGET,      CONFIGURE,   REPLOT,
	YMIN,           YMAX,           YRANGE,    XHAIR
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

    case YRANGE: {
	char buf[1024];
	double wx0, wy0, wx1, wy1;
	GetRasterCoords(t->raster, &wx0, &wy0, &wx1, &wy1);
	sprintf(buf, "%f %f",
		(wy0-t->ymin-10)/(t->ymax-t->ymin+20),
		(wy1-t->ymin-10)/(t->ymax-t->ymin+20));
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
} tdisp_arg;
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
	{NULL,	     0,	      0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    if (!Tcl_GetCommandInfo(interp, args.raster_win, &info))
	return TCL_ERROR;

    raster = (Tk_Raster*)info.clientData;

    if (NULL == (td = template_new(args.io, args.cnum, interp, raster)))
	return TCL_ERROR;

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

    //    opts[3] = "2";
    //    opts[6] = "-capstyle";
    //    opts[7] = "projecting"; // but or round
    //    opts[8] = NULL;

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

    /* Draw Environments automatically freed on raster destroy */

    free(t);
}

typedef struct {
    double y;          /* y coord */
    int col[3];        /* colours */
    int x[4];          /* coordinates */
    int t_strand;      /* template strand */
    int mq;            /* mapping qual */
    int rec;	       /* used for yspread */
} tline;


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
    double wx0, wy0, wx1, wy1, ny0, ny1, y;
    rangec_t *r;
    int nr, i, j, mode;
    struct timeval tv1, tv2;
    double t1, t2 = 0, t3, t4;
    double tsize = 1000;
    double yz;
    int ymin = INT_MAX;
    int ymax = INT_MIN;
    int width, height, height2;
    static int last_zoom = 0;
    tline *tl = NULL;
    int ntl = 0;
    int fwd_col, rev_col;
    double xgap;

    Display *rdisp;
    Drawable rdraw;
    GC rgc;

    rdisp = GetRasterDisplay (t->raster);
    rdraw = GetRasterDrawable(t->raster);
    rgc   = GetRasterGC      (t->raster);

    if (!t)
	return -1;
    if (!t->io)
	return -1;

    yz = t->yzoom / 200.0;

    fwd_col = t->yzoom >= 150 ? t->fwd_col3 : t->fwd_col;
    rev_col = t->yzoom >= 150 ? t->rev_col3 : t->rev_col;

    tk_RasterClear(t->raster);
    RasterWinSize(t->raster, &width, &height);
    height2 = height / 2;

    GetRasterCoords(t->raster, &wx0, &wy0, &wx1, &wy1);
    //RasterToWorld(t->raster, 0, 0, &wx0, &wy0);
    //RasterToWorld(t->raster, width, height, &wx1, &wy1);

    printf("templates %f to %f\n", wx0, wx1);

    //xgap = (int)(10*yz*(wx1-wx0)/width);
    xgap = 100;

    wx0 -= tsize;
    wx1 += tsize;

    if (t->depth_width != width) {
	t->depth_width  = width;
	t->sdepth = (int *)realloc(t->sdepth, width * sizeof(int));
	t->tdepth = (int *)realloc(t->tdepth, width * sizeof(int));
    }
    memset(t->sdepth, 0, t->depth_width * sizeof(int));
    memset(t->tdepth, 0, t->depth_width * sizeof(int));

    /* Find sequences on screen */
    gettimeofday(&tv1, NULL);
    mode = t->reads_only ? 0 : CSIR_PAIR;
    if (t->ymode == 1)
	mode |= CSIR_SORT_BY_Y;

    //mode = CSIR_SORT_BY_X | CSIR_PAIR | CSIR_ALLOCATE_Y_MULTIPLE;

    r = contig_seqs_in_range(t->io, &t->contig, wx0, wx1, mode, &nr);

    if (!r)
	return 0;

    tl = (tline *)malloc(nr * sizeof(*tl));

    gettimeofday(&tv2, NULL);
    t1 = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;

    /* Sort data */
    //qsort(r, nr, sizeof(*r), sort_by_mq);

    /*
     * Do this in multiple passes.
     * 1) Compute X (start/end of lines to draw), colours, status.
     * 2) Compute Y
     * 3) Plot (possibly combined with 2)
     */

    /* 1) Compute X */
    for (i = 0; i < nr; i++) {
	int sta, end;
	int r_sta, r_end, r_y;
	int span = 0;
	int single = 0;
	int col;
	int last_col = -1;
	double mq;

	sta = r[i].start;
	end = r[i].end;
	mq  = r[i].mqual;

	if (t->plot_depth) {
	    WorldToRaster (t->raster, sta, y, &r_sta, &r_y);
	    WorldToRaster (t->raster, end, y, &r_end, &r_y);
	    if (r_sta < 0) r_sta = 0;
	    if (r_end >= width) r_end = width-1;
	    for (j = r_sta; j <= r_end; j++) {
		t->sdepth[j]++;
	    }
	}

	if (t->filter & (FILTER_PAIRED | FILTER_SPANNING) && !r[i].pair_rec)
	    continue;

	if (t->filter & FILTER_SINGLE && r[i].pair_rec)
	    continue;

	tl[ntl].t_strand = ((r[i].flags & GRANGE_FLAG_END_REV) != 0)
	    == ((r[i].flags & GRANGE_FLAG_COMP1) != 0);
	tl[ntl].rec = r[i].rec;

	if (!t->reads_only && r[i].pair_rec) {
	    /* Only draw once */
	    if (!r[i].rec)
		continue;

	    if (r[i].pair_ind != -1) {
		r[r[i].pair_ind].rec = 0;
	    }

	    if (t->accuracy && !(r[i].pair_start || r[i].pair_end)) {
		/* Pair was off-screen, so get the results */
		seq_t *s;
		int crec;
		sequence_get_position(t->io, r[i].pair_rec, &crec,
				      &r[i].pair_start, &r[i].pair_end, NULL);

		s = (seq_t *)cache_search(t->io, GT_Seq, r[i].pair_rec);
		r[i].pair_mqual = s->mapping_qual;
		if (t->crec != crec) {
		    span = 1;
		}
	    }

	    if (!t->accuracy && !(r[i].pair_start || r[i].pair_end)) {
		span = 1;
	    }

	    if (t->filter & FILTER_SPANNING && !span)
		continue;
	    if (t->filter & FILTER_NONSPANNING && span)
		continue;

	    if (!span) {
		sta = r[i].start < r[i].pair_start
		    ? r[i].start : r[i].pair_start;
		end = r[i].end > r[i].pair_end
		    ? r[i].end : r[i].pair_end;
		r[i].flags |= GRANGE_FLAG_CONTIG;

		switch (t->cmode) {
		case 0:
		case 3:
		    mq = (r[i].mqual + r[i].pair_mqual)/2.0;
		    break;
		case 1:
		    mq = r[i].mqual < r[i].pair_mqual
			? r[i].mqual
			: r[i].pair_mqual;
		    break;
		case 2:
		    mq = r[i].mqual < r[i].pair_mqual
			? r[i].pair_mqual
			: r[i].mqual;
		    break;
		}

		if (t->plot_depth) {
		    WorldToRaster (t->raster, sta, y, &r_sta, &r_y);
		    WorldToRaster (t->raster, end, y, &r_end, &r_y);
		    if (r_sta < 0) r_sta = 0;
		    if (r_end >= width) r_end = width-1;
		    for (j = r_sta; j <= r_end; j++) {
			t->sdepth[j]++;
		    }
		}
	    }
	} else {
	    mq = r[i].mqual;
	    if (!t->reads_only)
		single = 1;
	}

	mq *= 3;
	if (mq < 0) mq = 0;
	if (mq > 255) mq = 255;
	tl[ntl].mq = mq;

	col = t->map_col[(int)(mq/8)];
	if (single)
	    col = t->single_col;
	if (span)
	    col = t->span_col;

	/* Check consistency */
	if (r[i].pair_rec && r[i].flags & GRANGE_FLAG_CONTIG) {
	    if (!((r[i].start < r[i].pair_start &&
		   (r[i].flags & GRANGE_FLAG_COMP1) == 0 &&
		   (r[i].flags & GRANGE_FLAG_COMP2) != 0) ||
		  (r[i].pair_start < r[i].start &&
		   (r[i].flags & GRANGE_FLAG_COMP1) != 0 &&
		   (r[i].flags & GRANGE_FLAG_COMP2) == 0)))
		col = t->inconsistent_col;
	}
	
	if (t->filter & FILTER_CONSISTENT && col == t->inconsistent_col)
	    continue;
	if (t->filter & FILTER_INCONSISTENT && col != t->inconsistent_col)
	    continue;

	if (t->plot_depth && !single && !span && col != t->inconsistent_col) {
	    WorldToRaster (t->raster, sta, 0, &r_sta, &r_y);
	    WorldToRaster (t->raster, end, 0, &r_end, &r_y);
	    if (r_sta < 0) r_sta = 0;
	    if (r_end >= width) r_end = width-1;
	    for (j = r_sta; j <= r_end; j++) {
		t->tdepth[j]++;
	    }
	}

	/* Filter */
	if (tl[ntl].mq < t->min_qual || tl[ntl].mq > t->max_qual)
	    continue;

	/* Generate line data */
	tl[ntl].y      = y;
	tl[ntl].col[1] = col;
	tl[ntl].x[0]   = sta;
	/* range in col 0 */
	tl[ntl].x[1]   = sta;
	/* range in col 1 */
	tl[ntl].x[2]   = end; 
	/* range in col 3 */
	tl[ntl].x[3]   = end;

	/* Readings too? */
	if (t->cmode == 3) {
	    //col = (r[i].flags & GRANGE_FLAG_END_MASK
	    //       == GRANGE_FLAG_END_FWD) ? t->fwd_col : t->rev_col;
	    col = (r[i].flags & GRANGE_FLAG_COMP1)
		? fwd_col : rev_col;

	    if (r[i].start == tl[ntl].x[0]) {
		tl[ntl].x[1] = r[i].end;
		tl[ntl].col[0] = col;
	    } else if (r[i].end == tl[ntl].x[3]) {
		tl[ntl].x[2] = r[i].start;
		tl[ntl].col[2] = col;
	    } else {
		printf("error, start/end do not match template pos\n");
	    }

	    if (r[i].pair_rec && (r[i].pair_start || r[i].pair_end)) {
		//col = (r[i].flags & GRANGE_FLAG_PEND_MASK
		//	   == GRANGE_FLAG_PEND_FWD) ? t->fwd_col : t->rev_col;
		col = (r[i].flags & GRANGE_FLAG_COMP2)
		    ? fwd_col : rev_col;

		if (r[i].pair_start == tl[ntl].x[0]) {
		    tl[ntl].x[1] = r[i].pair_end;
		    tl[ntl].col[0] = col;
		} else if (r[i].pair_end == tl[ntl].x[3]) {
		    tl[ntl].x[2] = r[i].pair_start;
		    tl[ntl].col[2] = col;
		} else {
		    printf("error, start/end do not match template pos\n");
		}
	    }
	}

	last_col = col;
	ntl++;
    }
    free(r);

    /* 2) Compute Y coordinates (part 1) */
    if (t->ymode == 1) {
	puts("Sorting");
	qsort(tl, ntl, sizeof(*tl), sort_tline_by_x);
	puts("computing ypos");
    tv1 = tv2;
    gettimeofday(&tv2, NULL);
    t2 = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;
	compute_ypos(t, xgap, tl, ntl);
    tv1 = tv2;
    gettimeofday(&tv2, NULL);
    t3 = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;
	puts("done");
    }


    /* 3) Plot the lines */
    for (i = 0; i < ntl; i++) {
	int last_col = -1;
	for (j = 0; j < 3; j++) {
	    //	    printf("%d: %d..%d in %d\n",
	    //		   i, tl[i].x[j], tl[i].x[j+1], tl[i].col[j]);
	    if (tl[i].x[j] < tl[i].x[j+1]) {
		double y;

		/* Compute the Y coordinates (part 2) */
		/* See XHAIR sub-command code too; keep it in sync */
		if (t->ymode == 1) {
		    y = tl[i].y*10;
		} else {
		    if (t->ymode == 0)
			y = tl[i].x[3] - tl[i].x[0];
		    else
			y = tl[i].mq*4;
		    if (t->logy) {
			if (y < 0) y = 0;
			y = 50*log(y+1);
		    }
		}
		y = (y + 50 - t->yoffset) * yz;

		if (t->spread)
		    y = y-t->spread/2+((tl[i].x[0]+tl[i].rec)%(t->spread));

		if (t->sep_by_strand)
		    y = tl[i].t_strand ? height2 - y : height2 + y;

		if (ymin > y) ymin = y;
		if (ymax < y) ymax = y;

		/* And plot if visible */
		if (y >= wy0 && y <= wy1) {
		    if (last_col != tl[i].col[j]) {
			last_col  = tl[i].col[j];
			SetDrawEnviron(t->interp, t->raster, last_col);
		    }
		    RasterDrawLine(t->raster, tl[i].x[j], y, tl[i].x[j+1], y);
		}
	    }
	}
    }
    free(tl);


    /* Plot depth */
    if (t->plot_depth) {
	XPoint *sp, *tp;
	sp = (XPoint *)calloc(width, sizeof(*sp));
	tp = (XPoint *)calloc(width, sizeof(*tp));
	for (j = 0; j < width; j++) {
	    sp[j].x = j;
	    sp[j].y = height-5 - t->sdepth[j] * yz;
	    tp[j].x = j;
	    tp[j].y = height-5 - t->tdepth[j] * yz;
	} 
	if (!(t->filter & FILTER_PAIRED)) {
	    SetDrawEnviron(t->interp, t->raster, t->fwd_col);
	    XDrawLines(rdisp, rdraw, rgc, sp, width, CoordModeOrigin);
	}
	SetDrawEnviron(t->interp, t->raster, t->rev_col);
	XDrawLines(rdisp, rdraw, rgc, tp, width, CoordModeOrigin);
	free(sp);
	free(tp);
   }


    /* Compute scrollbar and bounding box locations */
    t->ymin = ymin != INT_MAX ? ymin : 0;
    t->ymax = ymax != INT_MIN ? ymax : 0;

    tv1 = tv2;
    gettimeofday(&tv2, NULL);
    t4 = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;

    //    printf("Query range %d..%d => %d reads, %5.3fs + %5.3fs + %5.3fs + %5.3fs\n",
    //    	   (int)wx0, (int)wx1, nr, t1, t2, t3, t4);

    ny0 = wy0; ny1 = wy1;
    if (t->yzoom != last_zoom) {
	ny0 = t->ymin-10;
	ny1 = t->ymax+10;
    } else {
	if (ny0 > t->ymin-10)
	    ny0 = t->ymin-10;
	if (ny1 < t->ymax+10)
	    ny1 = t->ymax+10;
    }
    last_zoom = t->yzoom;
    if (ny0 > 0) ny0 = 0;
    if (ny1 < 400) ny1 = 400;
    wx0 = contig_get_start(&t->contig)-10;
    wx1 = contig_get_end(&t->contig)+10;
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
