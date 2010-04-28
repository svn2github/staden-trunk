#include <staden_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <tcl.h>
#include <tk.h>

#include <tg_gio.h>
#include "gap4_compat.h"
#include "io_utils.h"
#include "tk-io-reg.h"
#include "tcl_utils.h"
#include "IO.h"
#include "qual.h"
#include "qualIO.h"
#include "gap_cli_arg.h"
#include "newgap_structs.h"
#include "newgap_cmds.h"
#include "io_handle.h"
#include "tagdb.h"
#include "text_output.h"
#include "active_tags.h"
#include "xalloc.h"
#include "read_matrix.h"
#include "genetic_code.h"
#include "contig_selector.h"
#include "align.h"
#include "list_proc.h"
#include "find_repeats.h"
#include "tkEditor.h"
#include "tkEdNames.h"
#include "read_depth.h"
#include "fij.h"
#include "break_contig.h"
#include "template_display.h"
#include "export_contigs.h"
#include "find_oligo.h"
#include "tg_index_common.h"

#ifdef VALGRIND
#    include <valgrind/memcheck.h>
#endif

int tcl_get_tag_array(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    Tcl_DString tags;
    int i;

    get_tag_types();

    Tcl_DStringInit(&tags);
    for (i = 0; i < tag_db_count; i++) {
	Tcl_DStringStartSublist(&tags);
	Tcl_DStringAppendElement(&tags, tag_db[i].type);
	Tcl_DStringAppendElement(&tags, tag_db[i].search_id);
	Tcl_DStringAppendElement(&tags, tag_db[i].default_text);
	Tcl_DStringEndSublist(&tags);
    }
    Tcl_DStringResult(interp, &tags);

    return TCL_OK;
}


/*****************************************************************************
 *			      db_info
 *****************************************************************************/
/*
 * General syntax:
 * db_info io command ?-options? arguments ?arguments?
 *
 * Commands supported:
 * db_info num_readings io
 * db_info num_contigs io
 * db_info t_contig_length io
 * db_info t_read_length io
 * db_info get_read_num io identifier
 * db_info get_template_num io t_identifier
 * db_info get_contig_num io identifier
 * db_info chain_left io identifer (returns reading number)
 * db_info longest_contig io
 * db_info db_name io
 */
int db_info(ClientData clientData,
	    Tcl_Interp *interp,
	    int objc,
	    Tcl_Obj *CONST objv[])
{
    GapIO *io;
    char *cmd, *a3;

    if (objc < 3) {
	goto db_error;
    }

    cmd = Tcl_GetStringFromObj(objv[1], NULL);
    io = io_from_obj(objv[2]);

    if (strcmp(cmd, "num_readings") == 0) {
	//vTcl_SetResult(interp, "%d", NumReadings(io));
	vTcl_SetResult(interp, "%d", 1);
	return TCL_OK;
    }

    if (strcmp(cmd, "num_contigs") == 0) {
	//vTcl_SetResult(interp, "%d", NumContigs(io));
	vTcl_SetResult(interp, "%d", io->db->Ncontigs);
	return TCL_OK;
    }

    if (strcmp(cmd, "t_contig_length") == 0) {
	//vTcl_SetResult(interp, "%"PRId64, CalcTotalContigLen(io));
	vTcl_SetResult(interp, "%d", 1);

	return TCL_OK;
    }

    if (strcmp(cmd, "t_read_length") == 0) {
	//vTcl_SetResult(interp, "%"PRId64, CalcTotalReadingLen(io, NumReadings(io)));
	vTcl_SetResult(interp, "%d", 1);
	return TCL_OK;
    }

    if (strcmp(cmd, "get_read_num") == 0) {
	if (objc != 4) {
	    goto db_error;
	}
	a3 = Tcl_GetStringFromObj(objv[3], NULL);
	vTcl_SetResult(interp, "%d", get_gel_num(io, a3, GGN_ID));
	return TCL_OK;
    }

    if (strcmp(cmd, "get_template_num") == 0) {
	if (objc != 4) {
	    goto db_error;
	}
	a3 = Tcl_GetStringFromObj(objv[3], NULL);
	vTcl_SetResult(interp, "%d", template_name_to_number(io, a3));
	return TCL_OK;
    }

    if (strcmp(cmd, "get_contig_num") == 0) {
	if (objc != 4) {
	    goto db_error;
	}
	a3 = Tcl_GetStringFromObj(objv[3], NULL);
	vTcl_SetResult(interp, "%d", get_contig_num(io, a3, GGN_ID));
	return TCL_OK;
    }

    if (strcmp(cmd, "get_contig_nums") == 0) {
	int inArgc, outArgc;
	char **inArgv = NULL;
	contig_list_t *outArgv = NULL;
	Tcl_Obj *lobj, *iobj;
	int i;

	/* A list version of the above function */
	if (objc != 4) {
	    goto db_error;
	}

	a3 = Tcl_GetStringFromObj(objv[3], NULL);
	if (Tcl_SplitList(interp, a3, &inArgc, &inArgv) != TCL_OK)
	    return TCL_ERROR;

	if (-1 == lget_contig_num(io, inArgc, inArgv, &outArgc, &outArgv))
	    return TCL_ERROR;

	Tcl_Free((char *)inArgv);
	
	if (NULL == (lobj = Tcl_NewListObj(0, NULL)))
	    return TCL_ERROR;
	Tcl_IncrRefCount(lobj);

	for (i = 0; i < outArgc; i++) {
	    iobj = Tcl_NewIntObj(outArgv[i].contig);
	    Tcl_ListObjAppendElement(interp, lobj, iobj);
	}
	xfree(outArgv);

	Tcl_SetObjResult(interp, lobj);
	Tcl_DecrRefCount(lobj);
	return TCL_OK;
    }

    if (strcmp(cmd, "chain_left") == 0) {
	int i;

	if (objc != 4) {
	    goto db_error;
	}

	a3 = Tcl_GetStringFromObj(objv[3], NULL);
	i = get_gel_num(io, a3, GGN_ID);
	vTcl_SetResult(interp, "%d", i == -1 ? -1 : chain_left(io, i));
	return TCL_OK;
    }

    if (strcmp(cmd, "longest_contig") == 0) {
	//vTcl_SetResult(interp, "%"PRId64, CalcLongContig(io));
	vTcl_SetResult(interp, "%d", arr(GCardinal, io->contig_order, 0));
	return TCL_OK;
    }

    if (strcmp(cmd, "db_name") == 0) {
	vTcl_SetResult(interp, "%s", io->name);
	return TCL_OK;
    }

 db_error:
    Tcl_SetResult(interp,
		  "wrong # args: should be \"db_info command ?args?\"\n",
		  TCL_STATIC);
    return TCL_ERROR;
}

int
ObjGetOps(ClientData clientData,
	  Tcl_Interp *interp,
	  int argc,
	  char *argv[])
{
    int inum, l;
    char *ops;

    if (argc != 3)
	return TCL_ERROR;
    inum = atoi(argv[2]);

    ops = obj_get_ops(inum);
    Tcl_SetVar(interp, argv[1], "", 0);

    if (NULL == ops)
	return TCL_OK;

    while (l = strlen(ops)) {
	Tcl_SetVar(interp, argv[1], ops, TCL_LIST_ELEMENT | TCL_APPEND_VALUE);
	ops += l+1;
    }

    return TCL_OK;
}

int
ObjInvokeOp(ClientData clientData,
	  Tcl_Interp *interp,
	  int argc,
	  char *argv[])
{
    int inum;
    int op;

    if (argc != 3)
	return TCL_ERROR;

    inum = atoi(argv[1]);
    op = atoi(argv[2]);

    obj_invoke_op(inum, op);

    return TCL_OK;
}

int
ObjInvokeNext(ClientData clientData,
	      Tcl_Interp *interp,
	      int argc,
	      char *argv[])
{
    void *ptr;

    if (argc != 2)
	return TCL_ERROR;

    ptr = TclPtr2C(argv[1]);

    obj_invoke_next(ptr);

    return TCL_OK;
}

int
ObjGetNext(ClientData clientData,
	      Tcl_Interp *interp,
	      int argc,
	      char *argv[])
{
    void *ptr;

    if (argc != 2)
	return TCL_ERROR;

    ptr = TclPtr2C(argv[1]);

    vTcl_SetResult(interp, "%d", obj_get_next(ptr));

    return TCL_OK;
}

int
ObjGetBrief(ClientData clientData,
	    Tcl_Interp *interp,
	    int argc,
	    char *argv[])
{
    int inum = atoi(argv[1]);
    char *desc;

    if (NULL != (desc = obj_get_brief(inum)))
	vTcl_SetResult(interp, "%s", desc);

    return TCL_OK;
}


int
tcl_load_alignment_matrix(ClientData clientData,
			  Tcl_Interp *interp,
			  int argc,
			  char *argv[])
{
    static char *nt_order = "ACGTURYMWSKDHVB-*";
    int **nt_matrix;

    if (argc != 2) {
	Tcl_SetResult(interp, "Usage: load_alignment_matrix filename\n",
		      TCL_STATIC);
	return TCL_ERROR;
    }

    if (nt_matrix = create_matrix(argv[1], nt_order)) {
	init_W128(nt_matrix, nt_order, 0);
	free_matrix(nt_matrix, nt_order);
	return TCL_OK;
    }

    vTcl_SetResult(interp, "%s: could not load", argv[1]);
    return TCL_ERROR;
}

int tcl_load_genetic_code(ClientData clientData, Tcl_Interp *interp,
			  int objc, Tcl_Obj *CONST objv[])
{
    read_enz_arg args;
    FILE *fp;

    cli_args a[] = {
	{"-filename",	ARG_STR, 1, NULL, offsetof(read_enz_arg, filename)},
	{NULL,		0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    if (NULL == (fp = fopen(args.filename, "r"))) {
	Tcl_SetResult(interp, "unable to open file\n", TCL_STATIC);
	return TCL_ERROR;
    }

    if (0 == read_global_genetic_code(fp)) {
	verror(ERR_WARN, "load_genetic_code",
	       "Could not read genetic code. Using standard code.");
	init_genetic_code();
	vTcl_SetResult(interp, "%d", -1);
    } else {
	vTcl_SetResult(interp, "%d", 0);
    }
    fclose(fp);


    return TCL_OK;
}


int tcl_contig_order_to_number(ClientData clientData, Tcl_Interp *interp,
			       int objc, Tcl_Obj *CONST objv[]) {
    ord2num_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,	 1, NULL, offsetof(ord2num_arg, io)},
	{"-order",	ARG_INT, 1, NULL, offsetof(ord2num_arg, order)},
	{NULL,		0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;
	
    vTcl_SetResult(interp, "%d",
		   arr(GCardinal, args.io->contig_order, args.order));
    return TCL_OK;
}

int
tcl_save_contig_order(ClientData clientData,
		      Tcl_Interp *interp,
		      int objc,
		      Tcl_Obj *CONST objv[])
{
    list2_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    int i;
    GapIO *io;
    GCardinal *order;
    reg_order ro;
    reg_buffer_start rs;
    reg_buffer_end re;

    /* Parse arguments */
    cli_args a[] = {
	{"-io",	         ARG_IO,  1, NULL, offsetof(list2_arg, io)},
	{"-contigs",     ARG_STR, 1, NULL, offsetof(list2_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    /* Fetch the list of ordered contigs */
    io = args.io;
    active_list_contigs(io, args.inlist, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }

    /* Update the IO structure with the new contig and save it */
    io->contig_order = cache_rw(io, io->contig_order);
    order = ArrayBase(GCardinal, io->contig_order);
    for (i = 0; i < num_contigs; i++) {
	order[i] = contig_array[i].contig;
    }
    xfree(contig_array);

    cache_flush(io);

    /* Notify other displays of the change */
    /* Notify of the start of the flurry of updates */
    rs.job = REG_BUFFER_START;
    for (i = 1; i <= NumContigs(io); i++)
	contig_notify(io, i, (reg_data *)&rs);

    ro.job = REG_ORDER;

    for (i = 1; i <= NumContigs(io); i++) {
	ro.pos = order[i-1];
	contig_notify(io, i, (reg_data *)&ro);
    }

    /* Notify the end of our updates */
    re.job = REG_BUFFER_END;
    for (i = 1; i <= NumContigs(io); i++)
	contig_notify(io, 1, (reg_data *)&re);

    return TCL_OK;
}


int
tcl_flush_contig_order(ClientData clientData,
		      Tcl_Interp *interp,
		      int objc,
		      Tcl_Obj *CONST objv[])
{
    io_arg args;
    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(io_arg, io)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    args.io->contig_order = cache_rw(args.io, args.io->contig_order);
    cache_flush(args.io);

    return TCL_OK;
}


/*
 * contig selector commands
 */
int
DisplayContigSelector(ClientData clientData,
		      Tcl_Interp *interp,
		      int objc,
		      Tcl_Obj *CONST objv[])
{
    display_cs_arg args;
    tick_s *tick;
    tag_s tag;
    cursor_s cursor;

    cli_args a[] = {
	{"-io",	        ARG_IO,	 1, NULL, offsetof(display_cs_arg, io)},
	{"-window",     ARG_STR, 1, NULL, offsetof(display_cs_arg, window)},
	{"-frame",      ARG_STR, 1, NULL, offsetof(display_cs_arg, frame)},
	{"-tick_height",ARG_INT, 1, "-1", offsetof(display_cs_arg, tick_ht)},
	{"-tick_width", ARG_INT, 1, "-1", offsetof(display_cs_arg, tick_wd)},
	{"-tick_fill",  ARG_STR, 1,  "", offsetof(display_cs_arg, tick_fill)},
	{"-tag_width",  ARG_INT, 1, "-1", offsetof(display_cs_arg, tag_wd)},
	{"-tag_offset", ARG_INT, 1, "-1",offsetof(display_cs_arg, tag_offset)},
	{"-cursor_width",ARG_INT,1, "-1", offsetof(display_cs_arg, cursor_wd)},
	{"-cursor_fill",ARG_STR, 1,  "",offsetof(display_cs_arg, cursor_fill)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;


    tag = tag_struct(interp, gap5_defs, "CONTIG_SEL", args.tag_wd,
		     args.tag_offset);
    cursor = cursor_struct(interp, gap5_defs, "CONTIG_SEL", args.cursor_wd,
			   args.cursor_fill);
    tick = tick_struct(interp, gap5_defs, "CONTIG_SEL", args.tick_wd,
		       args.tick_ht, args.tick_fill);

    vTcl_SetResult(interp, "%d",
		   contig_selector_reg(interp, args.io, args.frame,
				       args.window, tag, cursor, tick));

    return TCL_OK;
}

int
UpdateContigOrder(ClientData clientData,
		  Tcl_Interp *interp,
		  int objc,
		  Tcl_Obj *CONST objv[])
{
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    int *contigs;

    update_order_arg args;
    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(update_order_arg, io)},
	{"-id",	     ARG_INT, 1, NULL, offsetof(update_order_arg, id)},
	{"-contigs", ARG_STR, 1, NULL, offsetof(update_order_arg, contigs)},
	{"-x",       ARG_INT, 1, NULL, offsetof(update_order_arg, x)},
	{NULL,	     0,	      0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }

    update_contig_order(interp, args.io, args.id, contig_array, num_contigs,
			args.x);

    xfree(contig_array);
    return TCL_OK;

}

int
ZoomCanvas(ClientData clientData,
	   Tcl_Interp *interp,
	   int objc,
	   Tcl_Obj *CONST objv[])
{
    s_zoom szoom;
    reg_generic gen;

    zoom_arg args;
    cli_args a[] = {
	{"-io",	       ARG_IO,  1, NULL, offsetof(zoom_arg, io)},
	{"-id",        ARG_INT, 1, NULL, offsetof(zoom_arg, id)},
	{"-r_id",      ARG_INT, 1, "-1", offsetof(zoom_arg, r_id)},
	{"-amount",    ARG_FLOAT, 1, "-1", offsetof(zoom_arg, amount)},
	{"-x1",        ARG_INT, 1, "-1", offsetof(zoom_arg, x1)},
	{"-y1",        ARG_INT, 1, "-1", offsetof(zoom_arg, y1)},
	{"-x2",        ARG_INT, 1, "-1", offsetof(zoom_arg, x2)},
	{"-y2",        ARG_INT, 1, "-1", offsetof(zoom_arg, y2)},
	{"-direction", ARG_STR, 1, "b", offsetof(zoom_arg, scroll)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    if (args.amount == -1 && args.x1 == -1 && args.y1 == -1 && args.x2 == -1
	&& args.y2 == -1) {
	gen.job = REG_GENERIC;
	gen.task = TASK_CANVAS_ZOOMBACK;

	result_notify(args.io, args.id, (reg_data *)&gen, 0);
	return TCL_OK;
    }

    gen.job = REG_GENERIC;
    gen.task = TASK_CANVAS_ZOOM;
    gen.data = (void *)&szoom;

    if (NULL == (szoom.zoom = (box *)xmalloc(sizeof(box))))
	return TCL_OK;

    szoom.amount = args.amount;
    szoom.r_id = args.r_id;
    szoom.zoom->x1 = args.x1;
    szoom.zoom->y1 = args.y1;
    szoom.zoom->x2 = args.x2;
    szoom.zoom->y2 = args.y2;
    sscanf(args.scroll, "%c", &szoom.scroll);

    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    xfree(szoom.zoom);
    return TCL_OK;

}


int
DeleteWindow(ClientData clientData,
	     Tcl_Interp *interp,
	     int objc,
	     Tcl_Obj *CONST objv[])
{
    reg_generic gen;
    delete_arg args;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(delete_arg, io)},
	{"-id",     ARG_INT, 1, NULL, offsetof(delete_arg, id)},
	{"-window", ARG_STR, 1, NULL, offsetof(delete_arg, window)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    gen.job = REG_GENERIC;
    gen.task = TASK_WINDOW_DELETE;
    gen.data = (void *)args.window;

    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    return TCL_OK;

}


int
DisplayContigComparator(ClientData clientData,
			Tcl_Interp *interp,
			int objc,
			Tcl_Obj *CONST objv[])
{
    display_cp_arg args;
    obj_cs *cs;

    cli_args a[] = {
	{"-io",	      ARG_IO,	 1, NULL, offsetof(display_cp_arg, io)},
	{"-id",	      ARG_INT,	 1, NULL, offsetof(display_cp_arg, id)},
	{"-window",   ARG_STR,	 1, NULL, offsetof(display_cp_arg, window)},
	{"-win_vertical",ARG_STR,1, NULL, offsetof(display_cp_arg, v_window)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    cs = result_data(args.io, args.id);

    vTcl_SetResult(interp, "%d",
		   contig_comparator_reg(interp, args.io, cs, args.window,
					 args.v_window));

    return TCL_OK;
}


int
DisplayCSDiagonal(ClientData clientData,
		  Tcl_Interp *interp,
		  int objc,
		  Tcl_Obj *CONST objv[])
{
    cs_tags_arg args;
    obj_cs *cs;
    char cmd[1024];
    int t_contig_length;
    cli_args a[] = {
	{"-io",	ARG_IO,  1, NULL, offsetof(cs_tags_arg, io)},
	{"-id", ARG_INT, 1, NULL, offsetof(cs_tags_arg, id)},
	{NULL,	0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    cs = result_data(args.io, args.id);

    t_contig_length = CalcTotalContigLen(args.io);
    sprintf(cmd, "%s create line 1 1 %d %d -tag diagonal",
	    cs->window, t_contig_length, t_contig_length);
    Tcl_Eval(interp, cmd);

    scaleSingleCanvas(interp, cs->world, cs->canvas, cs->window, 'b',
		      "diagonal");
    return TCL_OK;
}


int
FindRepeats(ClientData clientData,
	    Tcl_Interp *interp,
	    int objc,
	    Tcl_Obj *CONST objv[])
{
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    find_repeats_arg args;
    Tcl_DString input_params;
    char *name1;
    char *name2;
    char *name4;
    int mask;

    cli_args a[] = {
	{"-io",	       ARG_IO,  1, NULL, offsetof(find_repeats_arg, io)},
	{"-direction", ARG_INT, 1, "3",  offsetof(find_repeats_arg, idir)},
	{"-min_match", ARG_INT, 1, "25", offsetof(find_repeats_arg, minmat)},
	{"-contigs",   ARG_STR, 1, NULL, offsetof(find_repeats_arg, inlist)},
	{"-outfile",   ARG_STR, 1, "",	 offsetof(find_repeats_arg, outfile)},
	{"-tag_types", ARG_STR,	1, "",   offsetof(find_repeats_arg, tag_list)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("find repeats");
    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }

    mask = *args.tag_list ? 3 : 0;

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "Contigs: %s\n", args.inlist);

    name1 = get_default_string(interp, gap5_defs, "FINDREP.MINREP.NAME");
    name2 = get_default_string(interp, gap5_defs,
			       vw("FINDREP.SELTASK.BUTTON.%d",
				  args.idir));

    if (mask) {
	name4 = get_default_string(interp, gap5_defs,
				   "FINDREP.SELMODE.BUTTON.1");
    } else {
	name4 = get_default_string(interp, gap5_defs,
				   "FINDREP.SELMODE.BUTTON.2");
    }

    vTcl_DStringAppend(&input_params, "%s: %d\n%s\n%s %s\n",
		       name1, args.minmat, name2, name4, args.tag_list);

    if (*args.outfile) {
	vTcl_DStringAppend(&input_params, "Saved tags to file %s\n",
			   args.outfile);
    }
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    if (SetActiveTags(args.tag_list) == -1) {
	return TCL_OK;
    }

    if (find_repeats(args.io, args.idir, args.minmat,
		     mask, consensus_cutoff, num_contigs, contig_array,
		     *args.outfile ? args.outfile : NULL) < 0 ) {
	verror(ERR_WARN, "Find repeats", "Failure in Find Repeats");
	SetActiveTags("");
	return TCL_OK;
    }

    SetActiveTags("");
    if (contig_array)
	xfree(contig_array);
    return TCL_OK;

} /* end FindRepeats */

int tcl_list_confidence(ClientData clientData, Tcl_Interp *interp,
			int objc, Tcl_Obj *CONST objv[])
{
    int rargc;
    contig_list_t *rargv;
    list_conf_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(list_conf_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,  offsetof(list_conf_arg, inlist)},
	{"-summary",    ARG_INT, 1, "1",   offsetof(list_conf_arg, summary)},
	{NULL,	    0,	     0, NULL, 0}
    };
    int i, j;
    int *freqs;
    int freqs_tot[101];
    int length_tot;

    vfuncheader("list confidence");

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    for (j = 0; j <= 100; j++) freqs_tot[j] = 0;
    length_tot = 0;
    for (i = 0; i < rargc; i++) {
	freqs = count_confidence(args.io, rargv[i].contig, rargv[i].start,
				 rargv[i].end);
	if (!freqs) {
	    verror(ERR_WARN, "list_confidence",
		   "Failed in count confidence frequencies");
	    continue;
	}
	for (j = 0; j <= 100; j++) freqs_tot[j] += freqs[j];
	if (!args.summary) {
	    vmessage("---Contig %s---\n",
		     get_contig_name(args.io, rargv[i].contig));
	    list_confidence(freqs, rargv[i].end - rargv[i].start + 1);
	}
	length_tot += rargv[i].end - rargv[i].start + 1;
    }

    if (rargc > 1 || args.summary) {
	vmessage("---Combined totals---\n");
	list_confidence(freqs_tot, length_tot);
    }

    xfree(rargv);
    return TCL_OK;
}

int tcl_list_base_confidence(ClientData clientData, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[])
{
    int rargc;
    contig_list_t *rargv;
    list_conf_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(list_conf_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,  offsetof(list_conf_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };
    int freqmat[256], freqmis[256];
    int i;

    vfuncheader("list base confidence");

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    memset(freqmat, 0, 256 * sizeof(int));
    memset(freqmis, 0, 256 * sizeof(int));

    for (i = 0; i < rargc; i++) {
	if (-1 == get_base_confidences(args.io, rargv[i].contig,
				       freqmat, freqmis)) {
	    verror(ERR_WARN, "list_base_confidence",
		   "Failed to get base confidences");
	    continue;
	}

    }

    vTcl_SetResult(interp, "%f",
		   list_base_confidence(freqmat, freqmis)
		   );

    xfree(rargv);
    return TCL_OK;
}

int tcl_calc_consensus(ClientData clientData, Tcl_Interp *interp,
		       int objc, Tcl_Obj *CONST objv[])
{
    int rargc;
    contig_list_t *rargv;
    list2_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(list2_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,  offsetof(list2_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    if (rargc >= 1) {
	char *buf;

	if (NULL == (buf = (char *)ckalloc(rargv[0].end - rargv[0].start + 2)))
	    return TCL_OK;

	calculate_consensus_simple(args.io, rargv[0].contig, rargv[0].start,
				   rargv[0].end, buf, NULL);
	buf[rargv[0].end - rargv[0].start + 1] = 0;
	Tcl_SetResult(interp, buf, TCL_VOLATILE);
	ckfree(buf);
    }

    xfree(rargv);
    return TCL_OK;
}

int tcl_calc_quality(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[])
{
    int rargc;
    contig_list_t *rargv;
    list2_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(list2_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,  offsetof(list2_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    if (rargc >= 1) {
	char *buf;
	float *flt;
	int len = rargv[0].end - rargv[0].start + 1;
	int i;
	
	if (NULL == (flt = (float *)xmalloc(len * sizeof(float))))
	    return TCL_ERROR;
	if (NULL == (buf = (char *)xmalloc(len)))
	    return TCL_ERROR;

	calculate_consensus_simple(args.io, rargv[0].contig, rargv[0].start,
				   rargv[0].end, NULL, flt);
	for (i = 0; i < len; i++) {
	    int q = rint(flt[i]);
	    if (q < -127) q = -127;
	    if (q > 127) q = 127;
	    buf[i] = q;
	}

	Tcl_SetObjResult(interp, Tcl_NewStringObj(buf, len));

	xfree(flt);
	xfree(buf);
    }

    xfree(rargv);
    return TCL_OK;
}

/*
 * Returns the full consensus information as a Tcl list of lists
 * which each sub-list containing call; probabilities of A, C, G, T, *;
 * and sequence depth.
 */
int tcl_calc_consensus_full(ClientData clientData, Tcl_Interp *interp,
			    int objc, Tcl_Obj *CONST objv[])
{
    int rargc;
    contig_list_t *rargv;
    list2_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(list2_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,  offsetof(list2_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    if (rargc >= 1) {
	Tcl_Obj *cons = Tcl_NewListObj(0, NULL);
	consensus_t *c;
	int len = rargv[0].end - rargv[0].start + 2;
	int i;
	
	if (NULL == (c = (consensus_t *)xcalloc(len, sizeof(*c))))
	    return TCL_ERROR;

	calculate_consensus(args.io, rargv[0].contig, rargv[0].start,
			    rargv[0].end, c);
	for (i = 0; i < len; i++) {
	    Tcl_Obj *items[7], *list;
	    int j;

	    items[0] = Tcl_NewIntObj(c[i].call);
	    for (j = 0; j < 5; j++)
		items[j+1] = Tcl_NewIntObj(rint(c[i].scores[j]));
	    items[6] = Tcl_NewIntObj(c[i].depth);
	    list = Tcl_NewListObj(7, items);
	    Tcl_ListObjAppendElement(interp, cons, list);
	}

	Tcl_SetObjResult(interp, cons);

	xfree(c);
    }

    xfree(rargv);
    return TCL_OK;
}

int tcl_sequence_depth(ClientData clientData, Tcl_Interp *interp,
		       int objc, Tcl_Obj *CONST objv[])
{
    int rargc;
    contig_list_t *rargv;
    list2_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(list2_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,  offsetof(list2_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    if (rargc >= 1) {
	min_max_avg_t *depth;
	Tcl_Obj *l = Tcl_NewListObj(0, NULL);
	int len = rargv[0].end - rargv[0].start + 1;
	int i, start, end, inc;

	depth = sequence_depth(args.io, rargv[0].contig,
			       rargv[0].start, rargv[0].end,
			       &start, &end, &inc);
	if (!depth)
	    return TCL_ERROR;

	Tcl_ListObjAppendElement(interp, l, Tcl_NewIntObj(start));
	Tcl_ListObjAppendElement(interp, l, Tcl_NewIntObj(end));
	Tcl_ListObjAppendElement(interp, l, Tcl_NewIntObj(inc));
	len = (end-start)/inc+1;
	for (i = 0; i < len; i++) {
	    Tcl_ListObjAppendElement(interp, l, Tcl_NewIntObj(depth[i].min));
	    Tcl_ListObjAppendElement(interp, l, Tcl_NewIntObj(depth[i].max));
	    Tcl_ListObjAppendElement(interp, l, Tcl_NewIntObj(depth[i].avg));
	}

	Tcl_SetObjResult(interp, l);

	xfree(depth);
    }

    xfree(rargv);
    return TCL_OK;
}


/* FIXME: this entire function and below needs a complete rewrite */
int
tcl_find_internal_joins(ClientData clientData, Tcl_Interp *interp,
			int objc, Tcl_Obj *CONST objv[])
{
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    fij_arg args;
    Tcl_DString input_params;
    char *name1;
    char *name2;
    int mode = 0, mask = 0;

    cli_args a[] = {
	{"-io",		  ARG_IO,   1, NULL,      offsetof(fij_arg, io)},
	{"-mask",	  ARG_STR,  1, "none",    offsetof(fij_arg, mask)},
	{"-mode",	  ARG_STR,  1, "all_all", offsetof(fij_arg, mode)},
	{"-min_overlap",  ARG_INT,  1, "20",   offsetof(fij_arg, min_overlap)},
	{"-max_pmismatch",ARG_FLOAT,1, "30.0",    offsetof(fij_arg, max_mis)},
	{"-word_length",  ARG_INT,  1, "4",	  offsetof(fij_arg, word_len)},
	{"-max_prob",     ARG_FLOAT,  1, "1.0e-8",offsetof(fij_arg, max_prob)},
	{"-min_match",    ARG_INT,  1, "20",	 offsetof(fij_arg, min_match)},
	{"-band",         ARG_INT,  1, "10",     offsetof(fij_arg, band)},
	{"-win_size",	  ARG_INT,  1, "0",       offsetof(fij_arg, win_size)},
	{"-max_dashes",	  ARG_INT,  1, "0",       offsetof(fij_arg, dash)},
	{"-min_conf",	  ARG_INT,  1, "8",       offsetof(fij_arg, min_conf)},
	{"-tag_types",	  ARG_STR,  1, "",        offsetof(fij_arg, tag_list)},
	{"-contigs",	  ARG_STR,  1, NULL,      offsetof(fij_arg, inlist)},
	{"-use_conf",	  ARG_INT,  1, "1",	  offsetof(fij_arg, use_conf)},
	{"-use_hidden",	  ARG_INT,  1, "1",	  offsetof(fij_arg, use_hidden)},
	{"-max_display",  ARG_INT,  1, "0",       offsetof(fij_arg, max_display)},
	{NULL,		  0,	    0, NULL,      0}
    };

    vfuncheader("find internal joins");

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    /* Parse mode and mask */
    if (strcmp(args.mask, "none") == 0)
	mask = 1;
    else if (strcmp(args.mask, "mark") == 0)
	mask = 2;
    else if (strcmp(args.mask, "mask") == 0)
	mask = 3;
    else {
	Tcl_SetResult(interp, "invalid mask mode", TCL_STATIC);
	return TCL_ERROR;
    }

    if (strcmp(args.mode, "all_all") == 0)
	mode = COMPARE_ALL;
    else if (strcmp(args.mode, "segment") == 0)
	mode = COMPARE_SINGLE;
    else {
	Tcl_SetResult(interp, "invalid fij mode", TCL_STATIC);
	return TCL_ERROR;
    }

    /* create contig name array */
    active_list_contigs(args.io, args.inlist, &num_contigs, &contig_array);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "Contigs: %s\n", args.inlist);

    name1 = get_default_string(interp, gap5_defs,
			       vw("FIJ.SELTASK.BUTTON.%d",
				  mode == COMPARE_ALL ? 1 : 2));
    vTcl_DStringAppend(&input_params, "%s\n", name1);

    name1 = get_default_string(interp, gap5_defs, "FIJ.MINOVERLAP.NAME");
    name2 = get_default_string(interp, gap5_defs, "FIJ.MAXMIS.NAME");
    vTcl_DStringAppend(&input_params, "%s: %d\n%s: %f\n",
		       name1, args.min_overlap,
		       name2, args.max_mis);

#if 0
    /* FIXME: Disabled for now as WINSIZE.NAME no longer exists */
    if ((args.win_size == 0) && (args.dash == 0)) {

	Tcl_DStringAppend(&input_params, "Not using hidden data\n", -1);
    } else {
	name1 = get_default_string(interp, gap5_defs,
				   "FIJ.HIDDEN.WINSIZE.NAME");
	name2 = get_default_string(interp, gap5_defs,
				   "FIJ.HIDDEN.MAXDASH.NAME");
	vTcl_DStringAppend(&input_params, "Hidden data: %s: %d\n%s: %d\n",
			   name1, args.win_size, name2, args.dash);
    }
#endif

    name1 = get_default_string(interp, gap5_defs,
			       vw("FIJ.SELMODE.BUTTON.%d", mask));
    vTcl_DStringAppend(&input_params, "%s %s\n", name1, args.tag_list);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    if (SetActiveTags(args.tag_list) == -1) {
	return TCL_OK;
    }

    if (fij(args.io, mask, mode, args.min_overlap, (double)args.max_mis,
	    args.word_len, (double)args.max_prob, args.min_match, args.band,
	    args.win_size, args.dash, args.min_conf, args.use_conf,
	    args.use_hidden, args.max_display,
	    num_contigs, contig_array) < 0 ) {
	verror(ERR_WARN, "Find internal joins", "Failure in Find Internal Joins");
	SetActiveTags("");
	return TCL_OK;
    }

    SetActiveTags("");
    xfree(contig_array);

    return TCL_OK;

} /* end FIJ */

int tcl_complement_contig(ClientData clientData, Tcl_Interp *interp,
		       int objc, Tcl_Obj *CONST objv[])
{
    int rargc, i;
    contig_list_t *rargv;
    list2_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(list2_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,  offsetof(list2_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("complement contig");

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    /* create contig name array */
    active_list_contigs(args.io, args.inlist, &rargc, &rargv);
    if (rargc == 0) {
        xfree(rargv);
        return TCL_OK;
    }

    for (i = 0; i < rargc; i++) {
	complement_contig(args.io, rargv[i].contig);
    }

    xfree(rargv);
    return TCL_OK;
}

int
tcl_break_contig(ClientData clientData, Tcl_Interp *interp,
		 int objc, Tcl_Obj *CONST objv[])
{
    break_contig_arg args;
    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(break_contig_arg, io)},
	{"-contig", ARG_INT, 1, NULL, offsetof(break_contig_arg, contig)},
	{"-pos",    ARG_INT, 1, NULL, offsetof(break_contig_arg, pos)},
	{NULL,	 0,	  0, NULL, 0}
    };

    vfuncheader("break contig");

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    if (break_contig(args.io, args.contig, args.pos) != 0) {
	Tcl_SetResult(interp, "Failure in Break Contig", TCL_STATIC);
	return TCL_ERROR;
    }

    return TCL_OK;
}

typedef struct {
    int fold;     /* line wrapping, 0 to disable */
    int shift;    /* add 'shift' to all output chars, default 0 */
    int phred;    /* boolean => convert from log-odds to phred */
    Tcl_Obj *str; /* The string to convert */
    int min;      /* Minimum ASCII value to print */
    int max;      /* Maximum ASCII value to print */
} format_sequence_arg;

int
tcl_reformat_sequence(ClientData clientData, Tcl_Interp *interp,
		 int objc, Tcl_Obj *CONST objv[])
{
    int i, j, k, len;
    char *in, *out;

    format_sequence_arg args;
    cli_args a[] = {
	{"-fold",  ARG_INT, 1, "0",   offsetof(format_sequence_arg, fold)},
	{"-shift", ARG_INT, 1, "0",   offsetof(format_sequence_arg, shift)},
	{"-phred", ARG_INT, 0, "0",   offsetof(format_sequence_arg, phred)},
	{"-str",   ARG_OBJ, 1, NULL,  offsetof(format_sequence_arg, str)},
	{"-min",   ARG_INT, 1, "0",   offsetof(format_sequence_arg, min)},
	{"-max",   ARG_INT, 1, "255", offsetof(format_sequence_arg, max)},
	{NULL,	 0,	  0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    in = Tcl_GetStringFromObj(args.str, &len);

    out = (char *)malloc(len + 1 + (args.fold ? len / args.fold + 1 : 0));
    if (!out)
	return TCL_ERROR;
    
    for (i = j = k = 0; i < len; i++, j++) {
	signed int c = in[i] + args.shift;
	if (c < args.min) c = args.min;
	if (c > args.max) c = args.max;
	out[j] = c;
	if (args.fold && ++k == args.fold) {
	    out[++j] = '\n';
	    k = 0;
	}
    }

    Tcl_SetObjResult(interp, Tcl_NewStringObj(out, j));
    free(out);

    return TCL_OK;
}

int
tcl_find_oligo(ClientData clientData,
	       Tcl_Interp *interp,
	       int objc,
	       Tcl_Obj *CONST objv[])
{
    oligo_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    Tcl_DString input_params;
    char *name1;

    cli_args a[] = {
	{"-io",		ARG_IO,	  1, NULL,  offsetof(oligo_arg, io)},
	{"-contigs",	ARG_STR,  1, NULL,  offsetof(oligo_arg, inlist)},
	{"-min_pmatch",	ARG_FLOAT,1, "75.", offsetof(oligo_arg, mis_match)},
	{"-tag_types",	ARG_STR,  1, "",    offsetof(oligo_arg, tag_list)},
	{"-seq",	ARG_STR,  1, "",    offsetof(oligo_arg, seq)},
	{"-consensus_only",
	                ARG_INT,  1, "0",   offsetof(oligo_arg, consensus_only)},
	{"-cutoffs",    ARG_INT,  1, "0",   offsetof(oligo_arg, cutoffs)},
	{"-file",	ARG_STR,  1, "",    offsetof(oligo_arg, file)},
	{NULL,		0,	  0, NULL,  0}
    };

    vfuncheader("sequence search");

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    /* create contig name array */
    active_list_contigs(args.io, args.inlist, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    free(contig_array);
	return TCL_OK;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "Contigs: %s\n", args.inlist);
    name1 = get_default_string(interp, gap5_defs, "FINDOLIGO.MAXMIS.NAME");

    vTcl_DStringAppend(&input_params, "%s: %f\n", name1, args.mis_match);

    if (*args.seq) {
	vTcl_DStringAppend(&input_params, "Sequence: %s\n", args.seq);
    } else if (*args.file) {
	vTcl_DStringAppend(&input_params, "File: %s\n", args.file);
    } else {
	vTcl_DStringAppend(&input_params, "Tags: %s\n", args.tag_list);
    }
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    /* check that if using sequence to probe, then it is DNA */
    /*
     * if (strcmp(args.tag_list, "") == 0) {
     *	 if (1 != get_seq_type (args.seq, strlen(args.seq))) {
     *	    verror(ERR_WARN, "find oligos", "sequence is not DNA");
     *	    return TCL_OK;
     *	  }
     * }
     */

    if (SetActiveTags(args.tag_list) == -1) {
	return TCL_ERROR;
    }

    if (args.file && *args.file) {
	if (-1 == find_oligo_file(args.io, num_contigs, contig_array,
				  args.mis_match, args.file,
				  args.consensus_only, args.cutoffs))
	    verror(ERR_FATAL, "find oligos", "could not search");
    } else {
	if (-1 == find_oligos(args.io, num_contigs, contig_array,
			      args.mis_match, args.seq,
			      args.consensus_only, args.cutoffs))
	    verror(ERR_FATAL, "find oligos", "out of memory");
    }

    SetActiveTags("");
    if (contig_array)
	xfree(contig_array);
    return TCL_OK;
}

typedef struct import_reads_arg {
    GapIO *io;
    char *data_type;
    char *comp_mode;
    char *file;
    char *fmt;
    tg_args a;
    int index_names;
} ir_arg;

int
tcl_import_reads(ClientData clientData,
	       Tcl_Interp *interp,
	       int objc,
	       Tcl_Obj *CONST objv[])
{
    ir_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    Tcl_DString input_params;
    char *name1;
    int fmt;

    /* Parse arguments */
    cli_args a[] = {
	{"-io",		   ARG_IO,  1, NULL,   offsetof(ir_arg, io)},
	{"-file",          ARG_STR, 1, NULL,   offsetof(ir_arg, file)},
	{"-data_type",     ARG_STR, 1, "all",  offsetof(ir_arg, data_type)},
	{"-comp_mode",     ARG_STR, 1, "zlib", offsetof(ir_arg, comp_mode)},
	{"-format",        ARG_STR, 1, "auto", offsetof(ir_arg, fmt)},
	{"-append",        ARG_INT, 1, "1",    offsetof(ir_arg, a.append)},
	{"-index_names",   ARG_INT, 1, "0",    offsetof(ir_arg, index_names)},
	{"-merge_contigs", ARG_INT, 1, "-1",   offsetof(ir_arg, a.merge_contigs)},
	{"-fast_mode",     ARG_INT, 1, "0",    offsetof(ir_arg, a.fast_mode)},
	{"-reserved_seqs", ARG_INT, 1, "0",    offsetof(ir_arg, a.reserved_seqs)},
	{"-repad",         ARG_INT, 1, "0",    offsetof(ir_arg, a.repad)},
	{"-pair_reads",    ARG_INT, 1, "1",    offsetof(ir_arg, a.pair_reads)},
	{NULL,		   0,	    0, NULL,   0}
    };

    vfuncheader("sequence search");

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    args.a.no_tree = args.index_names ? 0 : 1;
    args.a.data_type = parse_data_type(args.data_type);
    if (0 == strcmp(args.comp_mode, "none")) {
	args.a.comp_mode = COMP_MODE_NONE;
    } else if (0 == strcmp(args.comp_mode, "zlib")) {
	args.a.comp_mode = COMP_MODE_ZLIB;
    } else if (0 == strcmp(args.comp_mode, "lzma")) {
	args.a.comp_mode = COMP_MODE_LZMA;
    } else {
	fprintf(stderr, "Unknown compression mode '%s'\n", args.comp_mode);
	return TCL_ERROR;
    }


    /* Initialise io */
    args.io->iface->setopt(args.io->dbh, OPT_COMP_MODE, args.a.comp_mode);
    args.a.tmp = args.a.no_tree ? NULL : bttmp_file_open();


    /* Load data */
    if ((fmt = *args.fmt) == 'a')
	fmt = tg_index_file_type(args.file);

    switch(fmt) {
	case 'm':
	case 'M':
	    parse_maqmap(args.io, args.file, &args.a);
	    break;

	case 'A':
	    parse_ace(args.io, args.file, &args.a);
	    break;

	case 'B':
	    parse_baf(args.io, args.file, &args.a);
	    break;

#ifdef HAVE_SAMTOOLS
	case 'b':
	    parse_bam(args.io, args.file, &args.a);
	    break;
	case 's':	
	    parse_sam(args.io, args.file, &args.a);
	    break;
#endif

	default:
	    fprintf(stderr, "Unknown file type for '%s' - skipping\n",
		    args.file);
	    return TCL_ERROR;
    }


    /* Add to our sequence name B+Tree */
    if (args.a.tmp) {
	char *name;
	int rec;

	puts("Sorting sequence name index");
	bttmp_file_sort(args.a.tmp);

	puts("Building index");
	while (name = bttmp_file_get(args.a.tmp, &rec)) {
	    sequence_index_update(args.io, name, strlen(name), rec);
	}
	
	bttmp_file_close(args.a.tmp);
    }

    return TCL_OK;
}

#ifdef VALGRIND
tcl_leak_check(ClientData clientData,
	       Tcl_Interp *interp,
	       int objc,
	       Tcl_Obj *CONST objv[])
{
    VALGRIND_DO_LEAK_CHECK
}
#endif

/* set up tcl commands which call C procedures */
/*****************************************************************************/
/*				   NewGap_Init				     */
/*****************************************************************************/
int
NewGap_Init(Tcl_Interp *interp) {
    init_globals(interp);

    Tcl_CreateCommand(interp, "get_tag_array", tcl_get_tag_array,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "db_info", db_info,
			 (ClientData) NULL,
			 (Tcl_CmdDeleteProc *) NULL);

    /* Contig registration commands, see tk-io-reg.c */
    Tcl_CreateObjCommand(interp, "result_names", tk_result_names,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "register_id", tk_register_id,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "result_time", tk_result_time,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "result_delete", tk_result_delete,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "result_quit", tk_result_quit,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "contig_register", tk_contig_register,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "contig_deregister", tk_contig_deregister,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "contig_notify", tk_contig_notify,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "result_notify", tk_result_notify,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "reg_get_ops", tk_reg_get_ops,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "reg_invoke_op", tk_reg_invoke_op,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "reg_notify_update", tk_reg_notify_update,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "quit_displays", tcl_quit_displays,
			 (ClientData) NULL,
			 (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "reg_notify_highlight",
			 tk_reg_notify_highlight,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "result_is_2d", tk_result_is_2d,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateCommand(interp, "obj_get_ops", ObjGetOps,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateCommand(interp, "obj_invoke_op", ObjInvokeOp,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateCommand(interp, "obj_invoke_next", ObjInvokeNext,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateCommand(interp, "obj_get_next", ObjGetNext,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateCommand(interp, "obj_get_brief", ObjGetBrief,
			 (ClientData) NULL,
			 NULL);

    Tcl_CreateObjCommand(interp, "complement_contig", tcl_complement_contig,
			 (ClientData) NULL,
			 NULL);
#if 0
    Tcl_CreateObjCommand(interp, "cursor_ref", tk_cursor_ref,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "contig_lock_write", tk_contig_lock_write,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "query_cursor", tk_query_cursor,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "create_cursor", tk_create_cursor,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "delete_cursor", tk_delete_cursor,
			 (ClientData) NULL,
			 NULL);
#endif
    Tcl_CreateCommand(interp, "load_alignment_matrix",
		      tcl_load_alignment_matrix, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "load_genetic_code",
			 tcl_load_genetic_code, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "contig_order_to_number",
			 tcl_contig_order_to_number,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "save_contig_order", tcl_save_contig_order,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "flush_contig_order", tcl_flush_contig_order,
			 (ClientData) NULL,
			 NULL);

    /* Contig selector window */
    Tcl_CreateObjCommand(interp, "display_contig_selector",
			 DisplayContigSelector,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "display_contig_comparator",
			 DisplayContigComparator,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "update_contig_order", UpdateContigOrder,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "zoom_canvas", ZoomCanvas,
			 (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "delete_window", DeleteWindow,
			 (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "clear_cp", tk_clear_cp,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "matchresult_configure",
			 tk_matchresult_configure,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "display_cs_diagonal", DisplayCSDiagonal,
			 (ClientData) NULL,
			 NULL);

    Tcl_CreateObjCommand(interp, "find_repeats", FindRepeats,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "list_confidence",
			 tcl_list_confidence, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "list_base_confidence",
			 tcl_list_base_confidence, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "calc_consensus", tcl_calc_consensus,
			 (ClientData) NULL, NULL);
    Tcl_CreateObjCommand(interp, "calc_quality", tcl_calc_quality,
			 (ClientData) NULL, NULL);
    Tcl_CreateObjCommand(interp, "calc_consensus_full",
			 tcl_calc_consensus_full,
			 (ClientData) NULL, NULL);

    Tcl_CreateObjCommand(interp, "find_internal_joins",
			 tcl_find_internal_joins,
			 (ClientData) NULL, NULL);
    Tcl_CreateObjCommand(interp, "break_contig",
			 tcl_break_contig,
			 (ClientData) NULL, NULL);

    Tcl_CreateObjCommand(interp, "sequence_depth",
			 tcl_sequence_depth,
			 (ClientData) NULL, NULL);

    Tcl_CreateObjCommand(interp, "reformat_sequence",
			 tcl_reformat_sequence,
			 (ClientData) NULL, NULL);

    Tcl_CreateObjCommand(interp, "export_contigs",
			 tcl_export_contigs,
			 (ClientData) NULL, NULL);
    Tcl_CreateObjCommand(interp, "export_tags",
			 tcl_export_tags,
			 (ClientData) NULL, NULL);

    Tcl_CreateObjCommand(interp, "find_oligo",
			 tcl_find_oligo,
			 (ClientData) NULL, NULL);

    Tcl_CreateObjCommand(interp, "import_reads",
			 tcl_import_reads,
			 (ClientData) NULL, NULL);

#ifdef VALGRIND
    Tcl_CreateObjCommand(interp, "leak_check",
			 tcl_leak_check,
			 (ClientData) NULL, NULL);
#endif

    //Ced_Init(interp);
    Editor_Init(interp);
    EdNames_Init(interp);
    TDisp_Init(interp);

    return TCL_OK;

} /* end NewGap_Init */
