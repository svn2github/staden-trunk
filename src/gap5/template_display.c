#include <string.h>
#include <xalloc.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "template_display.h"
#include "gap_cli_arg.h"

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
    {TK_CONFIG_DOUBLE, "-yzoom", "yZoom", "YZoom", "10.0",
     -1, Tk_Offset(template_disp_t, yzoom), 0, 0, 0 /* mask */},
    {TK_OPTION_END}
};

static int tdisp_cmd(ClientData clientData, Tcl_Interp *interp,
		      int objc, Tcl_Obj *CONST objv[]) {
    int index, replot = 0;
    template_disp_t *t = (template_disp_t *)clientData;
    Tcl_Obj *res = NULL;

    static char *options[] = {
	"delete",       "io",           "cget",    "configure", "replot",
	"ymin",         "ymax",         "yrange",
	(char *)NULL,
    };

    enum options {
	DELETE,         IO,     	CGET,      CONFIGURE,   REPLOT,
	YMIN,           YMAX,           YRANGE
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
		(t->ymin-wy0)/(wy1-wy0),
		(t->ymax-wy0)/(wy1-wy0));
	Tcl_SetObjResult(interp, Tcl_NewStringObj(buf, -1));
	break;
    }
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

    opts[1] = "green";
    t->fwd_col = CreateDrawEnviron(interp, raster, 6, opts);

    opts[1] = "magenta";
    t->rev_col = CreateDrawEnviron(interp, raster, 6, opts);

    return t;
}

void template_destroy(template_disp_t *t) {
    if (!t)
	return;

    if (t->io && t->contig)
	cache_decr(t->io, t->contig);

    /* Draw Environments automatically freed on raster destroy */

    free(t);
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

int template_replot(template_disp_t *t) {
    double wx0, wy0, wx1, wy1, y;
    rangec_t *r;
    int nr, i, mode;
    struct timeval tv1, tv2;
    double t1, t2, t3;
    double tsize = 1000;
    double yz;
    int ymin = INT_MAX;
    int ymax = INT_MIN;
    int t_strand;
    int width, height;

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

    yz = t->yzoom / 100.0;

    tk_RasterClear(t->raster);
    GetRasterCoords(t->raster, &wx0, &wy0, &wx1, &wy1);
    RasterWinSize(t->raster, &width, &height);
    height /= 2;

    printf("Coordinates (%f,%f) - (%f,%f) yz %f\n",
	   wx0, wy0, wx1, wy1, yz);

    wx0 -= tsize;
    wx1 += tsize;

    /* Find sequences on screen */
    gettimeofday(&tv1, NULL);
    mode = t->reads_only ? 0 : CSIR_PAIR;
    if (t->ymode == 1) {
	r = contig_seqs_in_range(t->io, &t->contig, wx0, wx1, 
				 mode |
				 CSIR_ALLOCATE_Y_MULTIPLE |
				 CSIR_SORT_BY_Y, &nr);
    } else {
	r = contig_seqs_in_range(t->io, &t->contig, wx0, wx1, mode, &nr);
    }
    if (!r)
	return 0;

    gettimeofday(&tv2, NULL);
    t1 = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;

    /* Sort data */
    //    qsort(r, nr, sizeof(*r), sort_by_mq);

    tv1 = tv2;
    gettimeofday(&tv2, NULL);
    t2 = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;

    /* Plot */
    for (i = 0; i < nr; i++) {
	int sta, end;
	int span = 0;
	int single = 0;
	int col;
	int last_col = -1;
	double mq;

	sta = r[i].start;
	end = r[i].end;
	mq  = r[i].mqual;
	t_strand = ((r[i].flags & GRANGE_FLAG_END_REV) != 0)
	    == ((r[i].flags & GRANGE_FLAG_COMP1) != 0);

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

	    if (!span) {
		sta = r[i].start < r[i].pair_start
		    ? r[i].start : r[i].pair_start;
		end = r[i].end > r[i].pair_end
		    ? r[i].end : r[i].pair_end;
		r[i].flags |= GRANGE_FLAG_CONTIG;

		switch (t->cmode) {
		case 0:
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
		case 3:
		    mq = 99;
		    break;
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

	if (t->ymode == 1) {
	    y = r[i].y*10;
	} else {
	    if (t->ymode == 0)
		y = end - sta;
	    else
		y = mq*4;
	    if (t->logy) {
		if (y < 0) y = 0;
		y = 50*log(y+1);
	    }
	}
	y = (y + 50 - t->yoffset) * yz;

	if (t->spread)
	    y = y-t->spread/2+(sta%(t->spread));

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

	if (ymin > y) ymin = y;
	if (ymax < y) ymax = y;

	if (t->sep_by_strand)
	    y = t_strand ? height - y : height + y;

	if (y >= wy0 && y <= wy1) {
	    if (col != last_col) {
		SetDrawEnviron(t->interp, t->raster, col);
	    }
#if 0
	    {
		int r_sta, r_end, r_y;
		WorldToRaster (t->raster, sta, y, &r_sta, &r_y);
		WorldToRaster (t->raster, end, y, &r_end, &r_y);
		XDrawLine(rdisp, rdraw, rgc, r_sta, r_y, r_end, r_y);
	    }
#endif
	    RasterDrawLine(t->raster, sta, y, end, y);

	    /* Draw readings too */
	    if (t->cmode == 3) {
		//col = (r[i].flags & GRANGE_FLAG_END_MASK
		//       == GRANGE_FLAG_END_FWD) ? t->fwd_col : t->rev_col;
		col = (r[i].flags & GRANGE_FLAG_COMP1)
		    ? t->fwd_col : t->rev_col;
		SetDrawEnviron(t->interp, t->raster, col);
		RasterDrawLine(t->raster, r[i].start, y, r[i].end, y);
		if (r[i].pair_rec && (r[i].pair_start || r[i].pair_end)) {
		    //col = (r[i].flags & GRANGE_FLAG_PEND_MASK
		    //	   == GRANGE_FLAG_PEND_FWD) ? t->fwd_col : t->rev_col;
		    col = (r[i].flags & GRANGE_FLAG_COMP2)
			? t->fwd_col : t->rev_col;
		    SetDrawEnviron(t->interp, t->raster, col);
		    RasterDrawLine(t->raster, r[i].pair_start, y,
				   r[i].pair_end, y);
		}
	    }

	    last_col = col;
	}
    }

    t->ymin = ymin != INT_MAX ? ymin : 0;
    t->ymax = ymax != INT_MIN ? ymax : 0;

    tv1 = tv2;
    gettimeofday(&tv2, NULL);
    t3 = tv2.tv_sec - tv1.tv_sec + (tv2.tv_usec - tv1.tv_usec)/1e6;

    printf("Query range %d..%d => %d reads, %5.3fs + %5.3fs + %5.3fs\n",
	   (int)wx0, (int)wx1, nr, t1, t2, t3);

    free(r);

    return 0;
}
