/*  Last edited: Jan  7 11:56 2004 (mng) */
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
	vTcl_SetResult(interp, "%d", 1);
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
	//vTcl_SetResult(interp, "%s", io_name(io));
	vTcl_SetResult(interp, "???");
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
#ifdef TO_PORT
    list2_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    int *contigs;
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
    contigs = to_contigs_only(num_contigs, contig_array);

    /* Update the IO structure with the new contig and save it */
    order = ArrayBase(GCardinal, io->contig_order);
    for (i = 0; i < num_contigs; i++) {
	order[i] = contig_array[i].contig;
    }
    xfree(contig_array);

    ArrayDelay(io, io->db.contig_order, io->db.Ncontigs, io->contig_order);
    flush2t(io);

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
#endif
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
#if 0
    Tcl_CreateCommand(interp, "complement_contig", tk_complement_contig,
		      (ClientData) NULL,
		      NULL);
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

    Tcl_CreateObjCommand(interp, "display_contig_selector",
			 DisplayContigSelector,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "display_contig_comparator",
			 DisplayContigComparator,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "zoom_canvas", ZoomCanvas,
			 (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "find_repeats", FindRepeats,
			 (ClientData) NULL,
			 NULL);

    //Ced_Init(interp);
    Editor_Init(interp);
    EdNames_Init(interp);

    return TCL_OK;

} /* end NewGap_Init */
