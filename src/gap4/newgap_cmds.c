/*  Last edited: Jan  7 11:56 2004 (mng) */
#include <stdio.h>
#include <stdlib.h>
#include <tcl.h>
#include <tk.h>

#include "canvas_box.h"
#include "gap_canvas_box.h"
#include "ruler_display.h"
#include "io_utils.h"
#include "io-reg.h"
#include "tcl_utils.h"
#include "IO.h"
#include "qual.h"
#include "qualIO.h"
#include "gap_cli_arg.h"
#include "newgap_structs.h"
#include "tagUtils.h"
#include "newgap_cmds.h"
#include "template.h"
#include "template_display.h"
#include "active_tags.h"
#include "hash.h"
#include "io_handle.h"
#include "gap_array.h"
#include "fij.h"
#include "plot_quality.h"
#include "auto_assemble.h"
#include "dis_readings.h"
#include "find_repeats.h"
#include "break_contig.h"
#include "newgap_cmds.h"
#include "readpair.h"
#include "tk-io-reg.h"
#include "cs-object.h"
#include "complement.h"
#include "reactions.h"
#include "list_proc.h"
#include "misc.h"
#include "text_output.h"
#include "quality_plot.h"
#include "restriction_enzymes.h"
#include "sequence_formats.h"
#include "gap_globals.h"
#include "stop_codon.h"
#include "assemble_direct.h"
#include "locks.h"
#include "check_assembly.h"
#include "probe.h"
#include "find_oligo.h"
#include "show_relationships.h"
#include "extract.h"
#include "preass.h"
#include "copy_db.h"
#include "oligo_sel.h"
#include "consen.h"
#include "contig_selector.h"
#include "contig_order.h"
#include "clip.h"
#include "notes.h"
#include "active_tags.h"
#include "licence.h"
#include "read_matrix.h"
#include "genetic_code.h"
#include "consistency_display.h"
#include "confidence_graph.h"
#include "gap-create.h"
#include "dbcheck.h"
#include "alter_rel.h"
#include "dstrand.h"
#include "tclXkeylist.h"
#include "readpair_coverage.h"
#include "reading_coverage.h"
#include "strand_coverage.h"

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


/*****************************************************************************/
/*				   OpenDB				     */
/*****************************************************************************/
/* open the database of db_name and version version_num
 * initialise global arrays
 * initilaise tag array
 */
int
OpenDB(ClientData clientData,
       Tcl_Interp *interp,
       int argc,
       char *argv[])
{
    int status;
    int access = 0;
    int handle;
    GapIO *io;
    open_db_arg args;
    cli_args a[] = {
	{"-name",    ARG_STR, 1, NULL, offsetof(open_db_arg, db_name)},
	{"-version", ARG_STR, 1, "0", offsetof(open_db_arg, version)},
	{"-create",  ARG_INT, 1, "0", offsetof(open_db_arg, create)},
	{"-access",  ARG_STR, 1, "r", offsetof(open_db_arg, access_str)},
	{NULL,	     0,	      0, NULL, 0}
    };

    vfuncheader("open database");

    if (-1 == gap_parse_args(a, &args, argc, argv)) {
	Tcl_SetResult(interp, "wrong # args:\n", TCL_STATIC);
	return TCL_ERROR;
    }

    if (strcmp(args.access_str, "READONLY") == 0 ||
	strcmp(args.access_str, "r") == 0) {
	access = 1;
    } else if (strcmp(args.access_str, "WRITE") == 0 ||
	       strcmp(args.access_str, "rw") == 0) {
	access = 0;
    }

    io = open_db(args.db_name, args.version, &status, args.create, access);

    /* Licence type may have changed; eg DEMO to VIEWER */
    {
	char c[2];
	c[0] = get_licence_type();
	c[1] = '\0';

	Tcl_SetVar2(interp, "licence","type", c, TCL_GLOBAL_ONLY);
    }

    if (!io) {
	Tcl_ResetResult(interp);
	return TCL_OK; /* caught by calling language */
    }


    if ( (handle = get_free_handle(io)) < 0 ) {
	/* no free handles */
	xfree(io);
	verror(ERR_FATAL, "open_database", "no free io handles");
	return TCL_ERROR;
    }

    /* Set read_only Tcl variable to reflect the status */
    if (access == 1 || status == IO_READ_ONLY) {
	Tcl_SetVar(interp, "read_only", "1", TCL_GLOBAL_ONLY);
    } else {
	Tcl_SetVar(interp, "read_only", "0", TCL_GLOBAL_ONLY);
    }

    vTcl_SetResult(interp, "%d", handle);
    return TCL_OK;

} /* end OpenDB */

/*****************************************************************************/
/*				  CloseDB				     */
/*****************************************************************************/
int
CloseDB(ClientData clientData,
	Tcl_Interp *interp,
	int argc,
	char *argv[])
{
    handle_arg args;
    cli_args a[] = {
	{"-io",	    ARG_INT, 1, NULL, offsetof(handle_arg, handle)},
	{NULL,	    0,	     0, NULL, 0}
    };
    GapIO *io;

    vfuncheader("close database");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (io = io_handle(&args.handle)))
	return TCL_ERROR;

    if (close_db(io) == -1) {
	remove_handle(&args.handle);
	Tcl_SetResult(interp, "Failed to close database", TCL_STATIC);
	return TCL_ERROR;
    }

    remove_handle(&args.handle);

    return TCL_OK;
}

/*****************************************************************************/
/*				  CopyDB				     */
/*****************************************************************************/
int
CopyDB(ClientData clientData,
       Tcl_Interp *interp,
       int argc,
       char *argv[])
{
    char basename[DB_FILELEN], *p;
    copy_db_arg args;
    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(copy_db_arg, io)},
	{"-version",ARG_STR, 1, NULL, offsetof(copy_db_arg, version)},
	{"-collect",ARG_INT, 1, "0",  offsetof(copy_db_arg, collect)},
	{NULL,	    0,	     0, NULL, 0}
    };
    int ret;

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("copy database");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (p = strrchr(io_name(args.io), '.'))) {
	vTcl_SetResult(interp, "Malformed database names", TCL_STATIC);
	return TCL_ERROR;
    }

    strncpy(basename, io_name(args.io), p - io_name(args.io));
    basename[p - io_name(args.io)] = 0;

    if (strcmp(p+1, args.version) == 0) {
	verror(ERR_FATAL, "copy_database", "attempt to copy over ourself!");
	vTcl_SetResult(interp, "-1", TCL_STATIC);
	return TCL_OK;
    }

    if (args.io->client->generic.mode != G_LOCK_RO) {
      /* Save the contig order first */
      ArrayDelay(args.io, args.io->db.contig_order, args.io->db.Ncontigs,
		 args.io->contig_order);
      flush2t(args.io);
    }

    if (args.collect) {
	ret = copy_database_from(args.io, basename, args.version);
    } else {
	ret = cpdb(basename, p+1, args.version);
    }

    if (ret == -1)
	verror(ERR_FATAL, "copy_database", "copy failed");

    vTcl_SetResult(interp, "%d", ret);

    return TCL_OK;
}

int MinimalCoverage(ClientData clientData,
		    Tcl_Interp *interp,
		    int argc,
		    char *argv[])
{
    int rargc;
    contig_list_t *rargv;
    list2_arg args;
    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(list2_arg, io)},
	{"-contigs", ARG_STR, 1, NULL, offsetof(list2_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };
    char *res;

    vfuncheader("minimal coverage");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    res = minimal_coverage(args.io, rargc, rargv);

    xfree(rargv);
    Tcl_SetResult(interp, res, TCL_DYNAMIC);

    return TCL_OK;
}

int UnattachedReadings(ClientData clientData,
		       Tcl_Interp *interp,
		       int argc,
		       char *argv[])
{
    list2_arg args;
    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL,  offsetof(list2_arg, io)},
	{NULL,	    0,	     0, NULL,  0}
    };
    char *res;

    vfuncheader("unattached readings");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    res = unattached_reads(args.io);

    Tcl_SetResult(interp, res, TCL_DYNAMIC);

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
 * db_info get_contig_num io identifier
 * db_info chain_left io identifer (returns reading number)
 * db_info longest_contig io
 * db_info db_name io
 */
int db_info(ClientData clientData,
	    Tcl_Interp *interp,
	    int argc,
	    char *argv[])
{
    GapIO *io;
    int handle;


    if (argc < 3) {
	goto db_error;
    }

    handle = atoi(argv[2]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    if (strcmp(argv[1], "num_readings") == 0) {
	vTcl_SetResult(interp, "%d", NumReadings(io));
	return TCL_OK;
    }

    if (strcmp(argv[1], "num_contigs") == 0) {
	vTcl_SetResult(interp, "%d", NumContigs(io));
	return TCL_OK;
    }

    if (strcmp(argv[1], "t_contig_length") == 0) {
	vTcl_SetResult(interp, "%d", CalcTotalContigLen(io));

	return TCL_OK;
    }

    if (strcmp(argv[1], "t_read_length") == 0) {
	vTcl_SetResult(interp, "%d", CalcTotalReadingLen(io, NumReadings(io)));
	return TCL_OK;
    }

    if (strcmp(argv[1], "get_read_num") == 0) {
	if (argc != 4) {
	    goto db_error;
	}
	vTcl_SetResult(interp, "%d", get_gel_num(io, argv[3], GGN_ID));
	return TCL_OK;
    }

    if (strcmp(argv[1], "get_contig_num") == 0) {
	if (argc != 4) {
	    goto db_error;
	}
	vTcl_SetResult(interp, "%d", get_contig_num(io, argv[3], GGN_ID));
	return TCL_OK;
    }

    if (strcmp(argv[1], "chain_left") == 0) {
	int i;

	if (argc != 4) {
	    goto db_error;
	}

	i = get_gel_num(io, argv[3], GGN_ID);
	vTcl_SetResult(interp, "%d", i == -1 ? -1 : chain_left(io, i));
	return TCL_OK;
    }

    if (strcmp(argv[1], "longest_contig") == 0) {
	vTcl_SetResult(interp, "%d", CalcLongContig(io));
	return TCL_OK;
    }

    if (strcmp(argv[1], "db_name") == 0) {
	vTcl_SetResult(interp, "%s", io_name(io));
	return TCL_OK;
    }

 db_error:
    Tcl_SetResult(interp,
		  "wrong # args: should be \"db_info command ?args?\"\n",
		  TCL_STATIC);
    return TCL_ERROR;
}

/*
 * canvas box commands
 */
int
ScrollCanvas(ClientData clientData,
	     Tcl_Interp *interp,
	     int argc,
	     char *argv[])
{
    reg_generic gen;
    scroll_arg args;

    cli_args a[] = {
	{"-io",	            ARG_IO,  1, NULL, offsetof(scroll_arg, io)},
	{"-id",             ARG_INT, 1, NULL, offsetof(scroll_arg, id)},
	{"-xscrollcommand", ARG_STR, 1, "", offsetof(scroll_arg, xscroll)},
	{"-yscrollcommand", ARG_STR, 1, "", offsetof(scroll_arg, yscroll)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    gen.job = REG_GENERIC;

    /* if (strcmp(args.direction, "x") == 0) { */
    if (strcmp(args.xscroll, "") != 0) {
	gen.data = (void *)args.xscroll;
	gen.task = TASK_CANVAS_SCROLLX;
	result_notify(args.io, args.id, (reg_data *)&gen, 0);
    }
    if (strcmp(args.yscroll, "") != 0) {
	gen.data = (void *)args.yscroll;
	gen.task = TASK_CANVAS_SCROLLY;
	result_notify(args.io, args.id, (reg_data *)&gen, 0);
    }

    return TCL_OK;

}

int
ZoomCanvas(ClientData clientData,
	   Tcl_Interp *interp,
	   int argc,
	   char *argv[])
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

    if (-1 == gap_parse_args(a, &args, argc, argv))
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
ResizeCanvas(ClientData clientData,
	     Tcl_Interp *interp,
	     int argc,
	     char *argv[])
{

    reg_generic gen;
    resize_arg args;

    cli_args a[] = {
	{"-io",	ARG_IO,  1, NULL, offsetof(resize_arg, io)},
	{"-id", ARG_INT, 1, NULL, offsetof(resize_arg, id)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    gen.job = REG_GENERIC;
    gen.task = TASK_CANVAS_RESIZE;

    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    return TCL_OK;
}

int
DeleteCanvasCursor(ClientData clientData,
		   Tcl_Interp *interp,
		   int argc,
		   char *argv[])
{
    resize_arg args;
    reg_generic gen;

    cli_args a[] = {
	{"-io",	ARG_IO,  1, NULL, offsetof(resize_arg, io)},
	{"-id", ARG_INT, 1, NULL, offsetof(resize_arg, id)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    gen.job = REG_GENERIC;
    gen.task = TASK_CANVAS_CURSOR_DELETE;
    result_notify(args.io, args.id, (reg_data *)&gen, 0);
    return TCL_OK;
}

int
DrawCanvasCursor_X(ClientData clientData,
		   Tcl_Interp *interp,
		   int argc,
		   char *argv[])
{
    t_cursor_arg args;
    reg_generic gen;

    cli_args a[] = {
	{"-io",	ARG_IO,  1, NULL, offsetof(t_cursor_arg, io)},
	{"-id", ARG_INT, 1, NULL, offsetof(t_cursor_arg, id)},
	{"-x",  ARG_INT, 1, NULL, offsetof(t_cursor_arg, cx)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    gen.job = REG_GENERIC;
    gen.task = TASK_CANVAS_CURSOR_X;
    gen.data = (void *)&args.cx;
    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    return TCL_OK;
}

int
DrawCanvasCursor_Y(ClientData clientData,
		   Tcl_Interp *interp,
		   int argc,
		   char *argv[])
{
    t_cursor_arg args;
    reg_generic gen;

    cli_args a[] = {
	{"-io",	ARG_IO,  1, NULL, offsetof(t_cursor_arg, io)},
	{"-id", ARG_INT, 1, NULL, offsetof(t_cursor_arg, id)},
	{"-y",  ARG_INT, 1, NULL, offsetof(t_cursor_arg, cy)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    gen.job = REG_GENERIC;
    gen.task = TASK_CANVAS_CURSOR_Y;
    gen.data = (void *)&args.cy;
    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    return TCL_OK;
}


int
DeleteWindow(ClientData clientData,
	     Tcl_Interp *interp,
	     int argc,
	     char *argv[])
{
    reg_generic gen;
    delete_arg args;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(delete_arg, io)},
	{"-id",     ARG_INT, 1, NULL, offsetof(delete_arg, id)},
	{"-window", ARG_STR, 1, NULL, offsetof(delete_arg, window)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    gen.job = REG_GENERIC;
    gen.task = TASK_WINDOW_DELETE;
    gen.data = (void *)args.window;

    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    return TCL_OK;

}

/*
 * template display commands
 */

int DisplayTemplates(ClientData clientData,
		     Tcl_Interp *interp,
		     int argc,
		     char *argv[])
{
    disp_template_arg args;
    contig_list_t *contig_array = NULL;
    int *contigs;
    int num_contigs = 0;
    ruler_s *ruler;
    int width, bold;
    cursor_s cursor;

    cli_args a[] = {
	{"-io",		ARG_IO,	 1, NULL, offsetof(disp_template_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL, offsetof(disp_template_arg, contig)},
	{"-frame",	ARG_STR, 1, NULL, offsetof(disp_template_arg, frame)},
	{"-win_template",ARG_STR, 1, NULL,
	     offsetof(disp_template_arg, win_template)},
	{"-win_ruler", ARG_STR, 1, NULL,
	     offsetof(disp_template_arg, win_ruler)},
	{"-width",  ARG_INT, 1, "-1", offsetof(disp_template_arg, line_width)},
	{"-bold",   ARG_INT, 1, "-1", offsetof(disp_template_arg, line_bold)},
	{"-cursor_width", ARG_INT,1, "-1",
	     offsetof(disp_template_arg, cursor_wd)},
	{"-cursor_fill",  ARG_STR, 1,  "",
	     offsetof(disp_template_arg, cursor_fill)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncgroup(2, "template display");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contig, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }
    contigs = to_contigs_only(num_contigs, contig_array);
    xfree(contig_array);

    cursor = cursor_struct(interp, gap_defs, "TEMPLATE", args.cursor_wd,
			   args.cursor_fill);

    if (NULL == (ruler = (ruler_s *)xmalloc(sizeof(*(ruler)))))
	return -1;

    ruler = ruler_struct(interp, gap_defs, "TEMPLATE", 1);
    if (args.line_width == -1) {
	width = get_default_int(GetInterp(), gap_defs, "TEMPLATE.LINE_WIDTH");
    } else {
	width = args.line_width;
    }

    if (args.line_bold == -1) {
	bold = get_default_int(GetInterp(), gap_defs, "TEMPLATE.LINE_BOLD");
    } else {
	bold = args.line_bold;
    }

    /* 'contigs' is remembered by template_reg, so don't free it */
    vTcl_SetResult(interp, "%d",
		  template_reg(interp, args.io, contigs, num_contigs,
			       args.frame, args.win_template, args.win_ruler,
			       ruler, cursor, width, bold));
    return TCL_OK;
}

int
UpdateTemplateDisplay(ClientData clientData,
		      Tcl_Interp *interp,
		      int argc,
		      char *argv[])
{
    reg_generic gen;
    template_arg args;
    int i;
    obj_template_disp *t;
    int recalc;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(template_arg, io)},
	{"-id",	     ARG_INT, 1, NULL, offsetof(template_arg, id)},
	{"-recalc",  ARG_INT, 1, "0",  offsetof(template_arg, recalc)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /*
     * update template display - MUST do before other notifications because
     * the contig offsets are generated during the template display
     */
    t = result_data(args.io, args.id, 0);
    gen.job = REG_GENERIC;
    gen.task = TASK_CANVAS_REDRAW;
    recalc = args.recalc;
    gen.data = (void *)&recalc;
    result_notify(args.io, t->id, (reg_data *)&gen, 0);

    for (i = 0; i < t->num_wins; i++) {
	if (t->win_list[i]->id != t->id) {
	    result_notify(args.io, t->win_list[i]->id, (reg_data *)&gen, 0);
	}
    }

    return TCL_OK;
}

/*
 * recreate a template plot on a registered display
 */
int
AddTemplatePlot(ClientData clientData,
		Tcl_Interp *interp,
		int argc,
		char *argv[])
{
    reg_generic gen;
    template_arg args;
    int i;
    obj_template_disp *t;
    win winfo;
    int recalc;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(template_arg, io)},
	{"-id",	     ARG_INT, 1, NULL, offsetof(template_arg, id)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    t = result_data(args.io, args.id, 0);
    strcpy(t->window, t->t_win);

    t->canvas->height = get_default_int(GetInterp(), gap_defs,
					"TEMPLATE.PLOT_HEIGHT");
    t->canvas->width  = get_default_int(GetInterp(), gap_defs,
					"TEMPLATE.PLOT_WIDTH");

    /* add to window list */
    winfo.window = t->window;
    winfo.scroll = 'b';
    winfo.id = args.id;

    gen.job = REG_GENERIC;
    gen.task = TASK_WINDOW_ADD;
    gen.data = (void *)&winfo;
    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    gen.job = REG_GENERIC;
    gen.task = TASK_CANVAS_REDRAW;
    recalc = 1;
    gen.data = (void *)&recalc;
    result_notify(args.io, t->id, (reg_data *)&gen, 0);

    for (i = 0; i < t->num_wins; i++) {
	if (t->win_list[i]->id != t->id) {
	    result_notify(args.io, t->win_list[i]->id, (reg_data *)&gen, 0);
	}
    }

    return TCL_OK;
}

/*
 * return 1 (true) if there num_wins is < MAX_NUM_WINS else return 0
 */
int
TemplateWinFree(ClientData clientData,
		Tcl_Interp *interp,
		int argc,
		char *argv[])
{
    template_arg args;
    obj_template_disp *t;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(template_arg, io)},
	{"-id",	     ARG_INT, 1, NULL, offsetof(template_arg, id)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    t = result_data(args.io, args.id, 0);

    if (t->num_wins >= MAX_NUM_WINS) {
	vTcl_SetResult(interp, "%d", 0);
    } else {
	vTcl_SetResult(interp, "%d", 1);
    }
    return TCL_OK;
}

int
DisplayRuler(ClientData clientData,
	     Tcl_Interp *interp,
	     int argc,
	     char *argv[])
{
    ruler_s ruler;
    ruler_arg args;
    reg_generic gen;
    win winfo;

    cli_args a[] = {
	{"-io",	         ARG_IO, 1, NULL,  offsetof(ruler_arg, io)},
	{"-id",          ARG_INT, 1, NULL, offsetof(ruler_arg, id)},
	{"-win_ruler",   ARG_STR, 1, NULL, offsetof(ruler_arg, win_ruler)},
	{"-width",       ARG_INT, 1, "-1", offsetof(ruler_arg, line_width)},
	{"-fill",        ARG_STR, 1, "",   offsetof(ruler_arg, colour)},
	{"-offset",      ARG_INT, 1, "-1", offsetof(ruler_arg, offset)},
	{"-tick_ht",     ARG_INT, 1, "-1", offsetof(ruler_arg, tick_height)},
	{"-tick_wd",     ARG_INT, 1, "-1", offsetof(ruler_arg, tick_width)},
	{"-tick_fill",   ARG_STR, 1, "",   offsetof(ruler_arg, tick_colour)},
	{"-text_offset", ARG_INT, 1, "-1", offsetof(ruler_arg, text_offset)},
	{"-tag_offset",  ARG_INT, 1, "-1", offsetof(ruler_arg, tag_offset)},
	{"-tag_width",   ARG_INT, 1, "-1", offsetof(ruler_arg, tag_width)},
	{NULL,	0,	0, NULL, 0}
    };
    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;


    if (args.line_width == -1) {
	ruler.line_width = get_default_int(interp, gap_defs,
					   "RULER.LINE_WIDTH");
    }

    if (strcmp(args.colour, "") == 0) {
	ruler.colour = get_default_string(interp, gap_defs, "RULER.COLOUR");
    }
    if (args.offset == -1) {
	ruler.offset = get_default_int(interp, gap_defs, "RULER.OFFSET");
    }

    if (args.tick_height == -1) {
	ruler.tick.t.ht = get_default_int(interp, gap_defs,
					  "RULER.TICK_HEIGHT");
    }

    if (args.tick_width == -1) {
	ruler.tick.t.line_width = get_default_int(interp, gap_defs,
						  "RULER.TICK_WIDTH");
    }

    if (strcmp(args.tick_colour, "") == 0) {
	ruler.tick.t.colour = get_default_string(interp, gap_defs,
						 "RULER.TICK_COLOUR");
    }

    if (args.text_offset == -1) {
	ruler.tick.offset = get_default_int(interp, gap_defs,
					    "RULER.TEXT_OFFSET");
    }

    if (args.tag_offset == -1) {
	ruler.tag.offset = get_default_int(interp, gap_defs,
					   "RULER.TAG_OFFSET");
    }

    if (args.tag_width == -1) {
	ruler.tag.width = get_default_int(interp, gap_defs, "RULER.TAG_WIDTH");
    }

    strcpy(ruler.window, args.win_ruler);

    /* add to window list */
    winfo.window = ruler.window;
    winfo.scroll = 'x';
    winfo.id = args.id;

    gen.job = REG_GENERIC;
    gen.task = TASK_WINDOW_ADD;
    gen.data = (void *)&winfo;
    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    gen.job = REG_GENERIC;
    gen.task = TASK_DISPLAY_RULER;
    gen.data = (void *)&ruler;
    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    return TCL_OK;
}

int
DisplayRulerTicks(ClientData clientData,
		   Tcl_Interp *interp,
		   int argc,
		   char *argv[])
{
    r_ticks_arg args;
    reg_generic gen;
    cli_args a[] = {
	{"-io",	   ARG_IO,  1, NULL, offsetof(r_ticks_arg, io)},
	{"-id",    ARG_INT, 1, NULL, offsetof(r_ticks_arg, id)},
	{"-ticks", ARG_INT, 1, "0",  offsetof(r_ticks_arg, ticks)},
	{NULL,	0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;


    gen.job = REG_GENERIC;
    gen.task = TASK_DISPLAY_TICKS;
    gen.data = (void *)&args.ticks;
    result_notify(args.io, args.id, (reg_data *)&gen, 0);
    return TCL_OK;

}

/*****************************************************************************/
/*			      DisplayReadingTags			     */
/*****************************************************************************/
/* plot tags on readings and ruler in the template display.
 * Use create rectangle to get  round the problem of not seeing very small tags
 */
int
DisplayReadingTags(ClientData clientData,
		   Tcl_Interp *interp,
		   int argc,
		   char *argv[])
{
    r_tags_arg args;
    obj_template_disp *t;

    cli_args a[] = {
	{"-io",	ARG_IO,  1, NULL, offsetof(r_tags_arg, io)},
	{"-id", ARG_INT, 1, NULL, offsetof(r_tags_arg, id)},
	{NULL,	0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    t = result_data(args.io, args.id, 0);

    display_reading_tags(interp, args.io, t);

    scaleCanvas(interp, t->win_list, t->num_wins, "tag", t->world->visible,
		t->canvas);

    return TCL_OK;

} /* end DisplayReadingTags */


/*****************************************************************************/
/*			  PrintTemplateReadings				     */
/*****************************************************************************/
int
PrintTemplateReadings(ClientData clientData,
		     Tcl_Interp *interp,
		     int argc,
		     char *argv[])
{
    obj_template_disp *temp;
    template_read_arg args;
    template_c *t;
    gel_cont_t *gc;
    template_d *t_changes = NULL;
    int length;
    GTemplates te;
    Tcl_DString str;
    char tmp[100];
    item_t *it;
    int in_contig_list;

    cli_args a[] = {
	{"-io",	  ARG_IO,  1, NULL, offsetof(template_read_arg, io)},
	{"-id",	  ARG_INT, 1, NULL, offsetof(template_read_arg, id)},
	{"-tnum", ARG_INT, 1, NULL, offsetof(template_read_arg, t_num)},
	{NULL,	      0,       0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    Tcl_DStringInit(&str);

    temp = result_data(args.io, args.id, 0);
    t = temp->tarr[args.t_num];

    if (!t->gel_cont) {
	Tcl_DStringAppend(&str, "Status                  Unknown\n\n", -1);
	Tcl_DStringResult(interp, &str);
	return TCL_OK;
    }

    if (t->flags & TEMP_FLAG_SPANNING) {
	
	for (it = head(t->gel_cont); it; it = it->next) {
	    gel_cont_t *gc = (gel_cont_t *)(it->data);
	    if (inContigList(temp->contig, temp->num_contigs,
			     gc->contig)) {
		in_contig_list = 1;
	    } else {
		in_contig_list = 0;
	    }
	}
    }

    if (t->flags & TEMP_FLAG_SPANNING && in_contig_list) {

	FindTemplatePositions(args.io, temp->contig_offset, temp->contig,
			      temp->num_contigs, temp->tarr, &t_changes);

	length = ABS(t_changes[args.t_num].start-t_changes[args.t_num].end)+1;
	sprintf(tmp, "estimated length        %d\n", length);
	Tcl_DStringAppend(&str, tmp, -1);

	if (t_changes[args.t_num].consist == 0) {
	    Tcl_DStringAppend(&str,
			      "Status                  Inconsistent\n", -1);
	} else {
	    template_read(args.io, args.t_num, te);
	    if ((length < te.insert_length_min) ||
		(length > te.insert_length_max)) {
		Tcl_DStringAppend(&str,
				  "Status                  Inconsistent\n",
				  -1);
	    } else {
		Tcl_DStringAppend(&str, "Status                  Ok\n", -1);
	    }
	}

#ifdef REMOVE
	if (t_changes[args.t_num].consist == 1) {
	    Tcl_DStringAppend(&str, "Status                  Ok\n", -1);
	} else if (t_changes[args.t_num].consist == 2) {
	    template_read(args.io, args.t_num, te);
	    if ((length < te.insert_length_min) ||
		(length > te.insert_length_max)) {
		Tcl_DStringAppend(&str,
				  "Status                  Inconsistent\n",
				  -1);
	    } else {
		Tcl_DStringAppend(&str, "Status                  Ok\n",
				  -1);
	    }
	} else {
	    Tcl_DStringAppend(&str, "Status                  Inconsistent\n",
			      -1);
	}
#endif
	if (t_changes)
	    xfree(t_changes);

    } else {

	if (t->flags & (TEMP_FLAG_GUESSED_START | TEMP_FLAG_GUESSED_END)) {
	    sprintf(tmp, "estimated length        %d\n",
		    ABS(t->end - t->start));
	} else {
	    sprintf(tmp, "observed length         %d\n",
		    ABS(t->end - t->start));
	}
	Tcl_DStringAppend(&str, tmp, -1);

	if (t->consistency) {
	    Tcl_DStringAppend(&str, "Status                  Inconsistent - ",
			      -1);
	    if (t->consistency & TEMP_CONSIST_DIST)
		Tcl_DStringAppend(&str, "Distance ", -1);
	    if (t->consistency & TEMP_CONSIST_PRIMER)
		Tcl_DStringAppend(&str, "Primer ", -1);
	    if (t->consistency & TEMP_CONSIST_STRAND)
		Tcl_DStringAppend(&str, "Strand ", -1);
	    if (t->consistency & TEMP_CONSIST_UNKNOWN)
		Tcl_DStringAppend(&str, "Missing", -1);
	    Tcl_DStringAppend(&str, "\n", -1);
	} else {
	    Tcl_DStringAppend(&str, "Status                  Ok\n", -1);
	}

	if (t->flags & TEMP_FLAG_GUESSED_START)
	    Tcl_DStringAppend(&str, "Start position has been guessed\n", -1);

	if (t->flags & TEMP_FLAG_GUESSED_END)
	    Tcl_DStringAppend(&str, "End position has been guessed\n", -1);
    }
    for (it = head(t->gel_cont); it; it = it->next) {
	char read_name[DB_NAMELEN+1];
	gc = (gel_cont_t *)(it->data);

	strcpy(read_name, get_read_name(args.io, gc->read));

	sprintf(tmp, "Contains reading %s (%d) from contig %s (%d)\n",
		 read_name, gc->read,
		 get_contig_name(args.io, gc->contig),
		 io_clnbr(args.io, gc->contig));
	Tcl_DStringAppend(&str, tmp, -1);
    }

    Tcl_DStringAppend(&str, "\n", -1);
    Tcl_DStringResult(interp, &str);

    return TCL_OK;

} /* end PrintTemplateReadings */

int TemplateContig(ClientData clientData,
		   Tcl_Interp *interp,
		   int argc,
		   char *argv[])
{
    obj_template_disp *t;
    t_cursor_arg args;
    double wx, wy;

    cli_args a[] = {
	{"-io",	ARG_IO,  1, NULL, offsetof(t_cursor_arg, io)},
	{"-id", ARG_INT, 1, NULL, offsetof(t_cursor_arg, id)},
	{"-x",  ARG_INT, 1, NULL, offsetof(t_cursor_arg, cx)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    t = result_data(args.io, args.id, 0);

    /* fill in global position of cursor in label box */
    CanvasToWorld(t->canvas, args.cx, 0, &wx, &wy);
    vTcl_SetResult(interp, "%d",
		   find_cursor_contig(args.io, args.id, t->contig_offset,
				      t->contig, t->num_contigs, wx));

    return TCL_OK;
}

/*****************************************************************************/
/*			 UpdateTemplateContigOrder			     */
/*****************************************************************************/
int UpdateTemplateContigOrder(ClientData clientData,
			   Tcl_Interp *interp,
			   int argc,
			   char *argv[])
{
    t_order_arg args;
    contig_list_t *contig_array = NULL;
    int *contigs;
    int num_contigs = 0;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(t_order_arg, io)},
	{"-id",	     ARG_INT,  1, NULL,offsetof(t_order_arg, id)},
	{"-x",	     ARG_INT,  1, NULL,offsetof(t_order_arg, cx)},
	{"-contigs", ARG_STR, 1, NULL, offsetof(t_order_arg, contig)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contig, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }
    contigs = to_contigs_only(num_contigs, contig_array);
    xfree(contig_array);

    update_template_contig_order(interp, args.io, args.id, args.cx, contigs,
				 num_contigs);

    xfree(contigs);

    /* return the new list */
    {
	int i;
	obj_template_disp *t;
	t = result_data(args.io, args.id, 0);

	for (i = 0; i < t->num_contigs; i++) {
	    Tcl_AppendElement(interp, get_contig_name(args.io,
						      ABS(t->contig[i])));
	}
    }


    return TCL_OK;
}

/*****************************************************************************/
/*			 RefreshContigOrder	  		             */
/*****************************************************************************/
int RefreshContigOrder(ClientData clientData,
		       Tcl_Interp *interp,
		       int argc,
		       char *argv[])
{
    r_tags_arg args;
    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(r_tags_arg, io)},
	{"-id",	     ARG_INT,  1, NULL,offsetof(r_tags_arg, id)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    refresh_contig_order(interp, args.io, args.id);
    return TCL_OK;
}

/*****************************************************************************/
/*			 RemoveContigDuplicates  		             */
/*****************************************************************************/
int RemoveContigDuplicates(ClientData clientData,
			   Tcl_Interp *interp,
			   int argc,
			   char *argv[])
{
    duplicates_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs;
    int i;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(t_order_arg, io)},
	{"-contigs", ARG_STR, 1, NULL, offsetof(duplicates_arg, contig)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contig, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }

    for (i = 0; i < num_contigs; i++) {
	Tcl_AppendElement(interp, get_contig_name(args.io, contig_array[i].contig));
    }
    xfree(contig_array);

    return TCL_OK;
}

int Canvas_To_World(ClientData clientData,
		    Tcl_Interp *interp,
		    int argc,
		    char *argv[])
{
    t_cursor_arg args;
    reg_generic gen;
    double wx;
    task_world_t tw;

    cli_args a[] = {
	{"-io",	  ARG_IO,  1, NULL, offsetof(t_cursor_arg, io)},
	{"-id",   ARG_INT, 1, NULL, offsetof(t_cursor_arg, id)},
	{"-cnum", ARG_INT, 1, "0",  offsetof(t_cursor_arg, cnum)},
	{"-x",    ARG_INT, 1, NULL, offsetof(t_cursor_arg, cx)},
	{NULL,	0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    gen.job = REG_GENERIC;
    gen.task = TASK_CANVAS_WORLD;
    tw.canvasx = args.cx;
    tw.cnum = args.cnum;
    gen.data = (void *)&tw;
    result_notify(args.io, args.id, (reg_data *)&gen, 0);

    wx = ((task_world_t *)gen.data)->basex;
    vTcl_SetResult(interp, "%d", (int)wx);
    return TCL_OK;
}

int DisplayTemplateQuality(ClientData clientData,
			   Tcl_Interp *interp,
			   int argc,
			   char *argv[])
{
    contig_list_t *contig_array = NULL;
    int *contigs;
    int num_contigs = 0;
    quality_arg args;
    int id;

    cli_args a[] = {
	{"-io",	         ARG_IO,  1, NULL, offsetof(quality_arg, io)},
	{"-id",          ARG_INT, 1, NULL, offsetof(quality_arg, id)},
	{"-contigs",     ARG_STR, 1, NULL, offsetof(quality_arg, contig)},
	{"-frame",       ARG_STR, 1, NULL, offsetof(quality_arg, frame)},
	{"-win_quality", ARG_STR, 1, NULL, offsetof(quality_arg, quality)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("calculate quality");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contig, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }
    contigs = to_contigs_only(num_contigs, contig_array);
    xfree(contig_array);

    id = template_quality_reg(args.io, interp, contigs, num_contigs,
			      consensus_cutoff, quality_cutoff,
			      args.frame, args.quality, args.id);

    xfree(contigs);

    vTcl_SetResult(interp, "%d", id);

    return TCL_OK;
}

int DisplayQuality(ClientData clientData,
		   Tcl_Interp *interp,
		   int argc,
		   char *argv[])
{
    disp_quality_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    ruler_s *ruler;
    char *tmp;
    cursor_s cursor;

    cli_args a[] = {
	{"-io",	         ARG_IO,  1, NULL, offsetof(disp_quality_arg, io)},
	{"-contigs",     ARG_STR, 1, NULL, offsetof(disp_quality_arg, contig)},
	{"-frame",       ARG_STR, 1, NULL, offsetof(disp_quality_arg, frame)},
	{"-window",     ARG_STR, 1, NULL, offsetof(disp_quality_arg, quality)},
	{"-cursor_width", ARG_INT,1, "-1",
	     offsetof(disp_quality_arg, cursor_wd)},
	{"-cursor_fill",  ARG_STR, 1,  "",
	     offsetof(disp_quality_arg, cursor_fill)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("display quality");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contig, &num_contigs, &contig_array);

    cursor = cursor_struct(interp, gap_defs, "QUALITY", args.cursor_wd,
			   args.cursor_fill);

    ruler = ruler_struct(interp, gap_defs, "QUALITY", 1);

    ruler->start = contig_array[0].start;
    ruler->end = contig_array[0].end;

    tmp = get_default_string(interp, gap_defs, "QUALITY.RULER.WIN");
    sprintf(ruler->window, "%s%s", args.frame, tmp);

    vTcl_SetResult(interp, "%d",
		   quality_reg(args.io, interp, contig_array[0].contig,
			       contig_array[0].start,
			       contig_array[0].end,
			       consensus_cutoff, quality_cutoff, args.frame,
			       args.quality, ruler, cursor));

    xfree(contig_array);
    return TCL_OK;
}

/*
 * finds and plots cut sites generated by selected enzymes provided in 'inlist'
 * from a file (filename) on a contig consensus sequence given by 'contig_num'
 * Displays results in the template display
 */
int
PlotTemplateREnz(ClientData clientData,
	       Tcl_Interp *interp,
	       int argc,
	       char *argv[])
{
    single_enz_arg args;
    int id;
    contig_list_t *contig_array;
    int *contigs;
    int num_contigs;
    tick_s *tick;

    cli_args a[] = {
	{"-io",	        ARG_IO,	 1, NULL, offsetof(single_enz_arg, io)},
	{"-id",         ARG_INT, 1, NULL, offsetof(single_enz_arg, id)},
	{"-file",       ARG_STR, 1, NULL, offsetof(single_enz_arg, filename)},
	{"-frame",      ARG_STR, 1, "", offsetof(single_enz_arg, frame)},
	{"-window",     ARG_STR, 1, NULL, offsetof(single_enz_arg, plot)},
	{"-enzymes",    ARG_STR, 1, NULL, offsetof(single_enz_arg, inlist)},
	{"-num_enzymes",ARG_INT, 1, NULL, offsetof(single_enz_arg, num_items)},
	{"-contigs",    ARG_STR, 1, NULL, offsetof(single_enz_arg, contig)},
	{"-tick_height",ARG_INT, 1, "-1", offsetof(single_enz_arg, tick_ht)},
	{"-tick_width", ARG_INT, 1, "-1", offsetof(single_enz_arg, tick_wd)},
	{"-tick_fill", ARG_STR,  1,   "", offsetof(single_enz_arg, tick_fill)},
	{"-yoffset",    ARG_INT, 1, "-1", offsetof(single_enz_arg, yoffset)},
	{NULL,	        0,	 0, NULL, 0}
    };

    vfuncgroup(5, "restriction enzymes");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contig, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }
    contigs = to_contigs_only(num_contigs, contig_array);
    xfree(contig_array);

    tick = tick_struct(interp, gap_defs, "R_ENZ", args.tick_wd, args.tick_ht,
		       args.tick_fill);

    id = template_renz_reg(interp, args.io, contigs, num_contigs,
			   args.filename, args.frame,
			   args.plot, args.id, args.inlist, args.num_items,
			   tick, args.yoffset);

    xfree(contigs);

    vTcl_SetResult(interp, "%d", id);

    return TCL_OK;
}

/*
 * finds and plots cut sites generated by selected enzymes provided in 'inlist'
 * from a file (filename) on a contig consensus sequence given by 'contig_num'
 * The names of each enzyme are display in canvas 'win_name' and the cuts in
 * canvas 'plot'
 */
int
PlotREnz(ClientData clientData,
	 Tcl_Interp *interp,
	 int argc,
	 char *argv[])
{
    renz_arg args;
    contig_list_t *contigs;
    int num_contigs;
    ruler_s *ruler;
    int id;
    tick_s *tick;
    cursor_s cursor;

    cli_args a[] = {
	{"-io",		 ARG_IO,  1, NULL, offsetof(renz_arg, io)},
	{"-file",	 ARG_STR, 1, NULL, offsetof(renz_arg, filename)},
	{"-frame",	 ARG_STR, 1, NULL, offsetof(renz_arg, frame)},
	{"-win_names",	 ARG_STR, 1, NULL, offsetof(renz_arg, win_name)},
	{"-window",	 ARG_STR, 1, NULL, offsetof(renz_arg, plot)},
	{"-win_ruler",	 ARG_STR, 1, NULL, offsetof(renz_arg, win_ruler)},
	{"-enzymes",	 ARG_STR, 1, NULL, offsetof(renz_arg, inlist)},
	{"-num_enzymes", ARG_INT, 1, NULL, offsetof(renz_arg, num_items)},
	{"-contigs",	 ARG_STR, 1, NULL, offsetof(renz_arg, contigs)},
	{"-text_offset", ARG_INT, 1, NULL, offsetof(renz_arg, text_offset)},
	{"-text_fill",   ARG_STR, 1, NULL, offsetof(renz_arg, text_fill)},
	{"-tick_height", ARG_INT, 1, "-1", offsetof(renz_arg, tick_ht)},
	{"-tick_width",  ARG_INT, 1, "-1", offsetof(renz_arg, tick_wd)},
	{"-tick_fill",   ARG_STR, 1,   "", offsetof(renz_arg, tick_fill)},
	{"-cursor_width",ARG_INT,1, "-1", offsetof(renz_arg, cursor_wd)},
	{"-cursor_fill", ARG_STR, 1,  "", offsetof(renz_arg, cursor_fill)},
	{"-yoffset",	 ARG_INT, 1, NULL, offsetof(renz_arg, yoffset)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncgroup(5, "restriction enzymes");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &num_contigs, &contigs);

    if (num_contigs != 1) {
	printf("ONLY DEAL WITH SINGLE CONTIG \n");
    }

    cursor = cursor_struct(interp, tk_utils_defs, "R_ENZ", args.cursor_wd,
			   args.cursor_fill);

    tick = tick_struct(interp, tk_utils_defs, "R_ENZ", args.tick_wd,
		       args.tick_ht, args.tick_fill);

    ruler = ruler_struct(interp, tk_utils_defs, "R_ENZ", 0);

    ruler->start = contigs[0].start;
    ruler->end = contigs[0].end;
    strcpy(ruler->window, args.win_ruler);

    id = renz_reg(interp, args.io, args.filename, args.frame,
		  args.win_name, args.plot, args.inlist, args.num_items,
		  contigs[0].contig, contigs[0].start, contigs[0].end,
		  args.text_offset, args.text_fill, tick, args.yoffset,
		  ruler, cursor);
    vTcl_SetResult(interp, "%d", id);
    xfree(contigs);
    return TCL_OK;
}

/*
 * return the name of an enzyme given a file of names and an index into that
 * file
 */
int
GetREnzName(ClientData clientData,
	    Tcl_Interp *interp,
	    int argc,
	    char *argv[])
{
    enz_name_arg args;
    obj_t_renz *r;

    cli_args a[] = {
	{"-io",	    ARG_IO, 1, NULL, offsetof(enz_name_arg, io)},
	{"-id",	    ARG_INT, 1, NULL, offsetof(enz_name_arg, id)},
	{"-enzyme", ARG_INT, 1, NULL, offsetof(enz_name_arg, item)},
	{"-cnum",   ARG_INT, 1, NULL, offsetof(enz_name_arg, contig)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    r = result_data(args.io, args.id, args.contig);
    /* printf("name %s \n", r->r_enzyme[args.item].name); */
    if (r) {
	vTcl_SetResult(interp, "%s", r->r_enzyme[args.item].name);
    } else {
	vTcl_SetResult(interp, "No renz plot for id %d, contig %d\n",
		       args.id, args.contig);
	return TCL_ERROR;
    }

    return TCL_OK;
}

/*
 * return info on an enzyme given a contig_register id and an item index into
 * the r_enzyme array
 */
int
GetREnzInfo(ClientData clientData,
	    Tcl_Interp *interp,
	    int argc,
	    char *argv[])
{
    enz_info_arg args;
    reg_generic gen;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(enz_info_arg, io)},
	{"-id",	    ARG_INT, 1, NULL, offsetof(enz_info_arg, id)},
	{"-enzyme", ARG_INT, 1, NULL, offsetof(enz_info_arg, enzyme)},
	{"-cnum",   ARG_INT, 1, NULL, offsetof(enz_info_arg, contig)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    gen.job = REG_GENERIC;
    gen.task = TASK_RENZ_INFO;
    gen.data = (void *)&args.enzyme;

    vfuncgroup(5, "restriction enzymes");
    result_notify(args.io, args.id, (reg_data *)&gen, args.contig);

    return TCL_OK;
}

int
CreateREnzTags(ClientData clientData,
	       Tcl_Interp *interp,
	       int argc,
	       char *argv[])
{
    contig_list_t *contigs;
    int num_contigs;
    enz_tag_arg args;
    obj_renz *r;
    int result;
    char **enz_ids = NULL;
    int num_ids;

    cli_args a[] = {
	{"-io",	      ARG_IO,  1, NULL, offsetof(enz_tag_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(enz_tag_arg, id)},
	{"-contigs",  ARG_STR, 1, NULL, offsetof(enz_tag_arg, contigs)},
	{"-enum",     ARG_STR, 1, NULL, offsetof(enz_tag_arg, enz_list)},
	{"-tag_types",ARG_STR, 1, NULL, offsetof(enz_tag_arg, id_list)},
	{NULL,	      0,       0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &num_contigs, &contigs);
    if (num_contigs == 0 || !contigs) {
	if (contigs) xfree(contigs);
	return TCL_OK;
    }
    if (num_contigs != 1) {
	printf("ERROR: only supported for single contig. Processing first"
	       " contig only\n");
    }

    r = result_data(args.io, args.id, contigs[0].contig);

    /* create reading name array */
    if (Tcl_SplitList(interp, args.id_list, &num_ids, &enz_ids)
	!= TCL_OK) {
	return TCL_ERROR;
    }

    result = Create_REnz_Tags(args.io, contigs[0].contig, r, args.enz_list,
			      enz_ids, num_ids);
    vTcl_SetResult(interp, "%d", result);
    xfree(contigs);
    Tcl_Free((char *)enz_ids);
    return TCL_OK;
}


/*
 * finds and plots cut sites generated by selected enzymes provided in 'inlist'
 * from a file (filename) on a contig consensus sequence given by 'contig_num'
 * The names of each enzyme are display in canvas 'win_name' and the cuts in
 * canvas 'plot'
 */
int
PlotStopCodons(ClientData clientData,
	      Tcl_Interp *interp,
	      int argc,
	      char *argv[])
{
    stop_codon_arg args;
    int id;
    contig_list_t *contigs;
    int num_contigs;
    ruler_s *ruler;
    char *tmp;
    tick_s *tick;
    cursor_s cursor;

    cli_args a[] = {
	{"-io",	        ARG_IO,	 1, NULL, offsetof(stop_codon_arg, io)},
	{"-frame",      ARG_STR, 1, NULL, offsetof(stop_codon_arg, frame)},
	{"-win_names",  ARG_STR, 1, NULL, offsetof(stop_codon_arg, names)},
	{"-window",     ARG_STR, 1, NULL, offsetof(stop_codon_arg, plot)},
	{"-strand",     ARG_INT, 1, "0", offsetof(stop_codon_arg, strand)},
	{"-contigs",    ARG_STR, 1, NULL, offsetof(stop_codon_arg, contigs)},
	{"-tick_height",ARG_INT, 1, "-1", offsetof(stop_codon_arg, tick_ht)},
	{"-tick_width", ARG_INT, 1, "-1", offsetof(stop_codon_arg, tick_wd)},
	{"-tick_fill",  ARG_STR, 1,   "", offsetof(stop_codon_arg, tick_fill)},
	{"-cursor_width",ARG_INT,1, "-1", offsetof(stop_codon_arg, cursor_wd)},
	{"-cursor_fill", ARG_STR, 1,  "", offsetof(stop_codon_arg, cursor_fill)},
	{"-yoffset",    ARG_INT, 1, NULL, offsetof(stop_codon_arg, yoffset)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("plot stop codons");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &num_contigs, &contigs);
    if (num_contigs == 0 || !contigs) {
	if (contigs) xfree(contigs);
	return TCL_OK;
    }

    if (num_contigs != 1) {
	printf("ONLY DEAL WITH SINGLE CONTIG \n");
    }

    cursor = cursor_struct(interp, gap_defs, "CODON", args.cursor_wd,
			   args.cursor_fill);

    tick = tick_struct(interp, gap_defs, "CODON", args.tick_wd, args.tick_ht,
		       args.tick_fill);

    ruler = ruler_struct(interp, gap_defs, "CODON", 1);

    ruler->start = contigs[0].start;
    ruler->end = contigs[0].end;

    tmp = get_default_string(interp, gap_defs, "CODON.RULER.WIN");
    sprintf(ruler->window, "%s%s", args.frame, tmp);

    id = codon_reg(interp, args.strand, args.io, args.frame,
		   args.names, args.plot,
		   contigs[0].contig, contigs[0].start, contigs[0].end,
		   tick, args.yoffset, ruler, cursor);
    vTcl_SetResult(interp, "%d", id);
    xfree(contigs);
    return TCL_OK;
}

int
GetCodonName(ClientData clientData,
	     Tcl_Interp *interp,
	     int argc,
	     char *argv[])
{
    enz_name_arg args;
    obj_codon *s;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(enz_name_arg, io)},
	{"-codon",  ARG_INT, 1, NULL, offsetof(enz_name_arg, item)},
	{"-id",	    ARG_INT, 1, NULL, offsetof(enz_name_arg, id)},
	{"-cnum",   ARG_INT, 1, NULL, offsetof(enz_name_arg, contig)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    s = result_data(args.io, args.id, args.contig);

    if (args.item > 3)
	args.item -= 3;
    vTcl_SetResult(interp, "%s", s->codon[args.item]);
    return TCL_OK;
}

int
RefreshCodonMap(ClientData clientData,
		Tcl_Interp *interp,
		int argc,
		char *argv[])
{
    refresh_codon_arg args;
    obj_codon *s;
    char *seq = NULL;

    cli_args a[] = {
	{"-io",	    ARG_IO, 1,  NULL, offsetof(refresh_codon_arg, io)},
	{"-id",	    ARG_INT, 1, NULL, offsetof(refresh_codon_arg, id)},
	{"-cnum",   ARG_INT, 1, NULL, offsetof(refresh_codon_arg, contig)},
	{"-strand", ARG_INT, 1, NULL, offsetof(refresh_codon_arg, strand)},
	{"-update", ARG_INT, 1, NULL, offsetof(refresh_codon_arg, update)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("refresh stop codons");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    s = result_data(args.io, args.id, args.contig);

    if (args.update) {
	reg_generic gen;
	task_editor_getcon tc;

	gen.job = REG_GENERIC;
	gen.task = TASK_EDITOR_GETCON;
	gen.data = (void *)&tc;

	tc.lreg = 0;
	tc.rreg = 0;
	tc.con_cut = consensus_cutoff;
	tc.qual_cut = quality_cutoff;

	if (type_contig_notify(args.io, args.contig, REG_TYPE_EDITOR,
			       (reg_data *)&gen, 0) == -1)
	    return TCL_OK;

	seq = tc.con;
    }

    s->strand = args.strand;
    stop_codon_replot(interp, args.io, s, seq);

    if (seq)
	xfree(seq);

    return TCL_OK;
}

int
ShowRelationships(ClientData clientData,
		    Tcl_Interp *interp,
		    int argc,
		    char *argv[])
{
    show_relationships_arg args;
    contig_list_t *contigs;
    int num_contigs;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(show_relationships_arg, io)},
	{"-contigs",ARG_STR, 1, "",   offsetof(show_relationships_arg,
					       contigs)},
	{"-order",  ARG_INT, 1, "1",  offsetof(show_relationships_arg,
					       ordered)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    vfuncheader("show relationships");

    active_list_contigs(args.io, args.contigs, &num_contigs, &contigs);

    show_relationships(args.io, contigs, num_contigs, args.ordered);

    if (contigs)
	xfree(contigs);

    return TCL_OK;
}

int
FindInternalJoins(ClientData clientData,
		  Tcl_Interp *interp,
		  int argc,
		  char *argv[])
{
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    fij_arg args;
    GapIO *io;
    Tcl_DString input_params;
    char *name1;
    char *name2;
    int mode = 0, mask = 0;

    cli_args a[] = {
	{"-io",		  ARG_INT,  1, NULL,      offsetof(fij_arg, handle)},
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

    if (-1 == gap_parse_args(a, &args, argc, argv))
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

    if ( (io = io_handle(&args.handle)) == NULL){
	verror(ERR_FATAL, "find_internal_joins", "invalid io handle");
	return -1;
    }

    /* create contig name array */
    active_list_contigs(io, args.inlist, &num_contigs, &contig_array);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "Contigs: %s\n", args.inlist);

    name1 = get_default_string(interp, gap_defs,
			       vw("FIJ.SELTASK.BUTTON.%d",
				  mode == COMPARE_ALL ? 1 : 2));
    vTcl_DStringAppend(&input_params, "%s\n", name1);

    name1 = get_default_string(interp, gap_defs, "FIJ.MINOVERLAP.NAME");
    name2 = get_default_string(interp, gap_defs, "FIJ.MAXMIS.NAME");
    vTcl_DStringAppend(&input_params, "%s: %d\n%s: %f\n",
		       name1, args.min_overlap,
		       name2, args.max_mis);

#if 0
    /* FIXME: Disabled for now as WINSIZE.NAME no longer exists */
    if ((args.win_size == 0) && (args.dash == 0)) {

	Tcl_DStringAppend(&input_params, "Not using hidden data\n", -1);
    } else {
	name1 = get_default_string(interp, gap_defs,
				   "FIJ.HIDDEN.WINSIZE.NAME");
	name2 = get_default_string(interp, gap_defs,
				   "FIJ.HIDDEN.MAXDASH.NAME");
	vTcl_DStringAppend(&input_params, "Hidden data: %s: %d\n%s: %d\n",
			   name1, args.win_size, name2, args.dash);
    }
#endif

    name1 = get_default_string(interp, gap_defs,
			       vw("FIJ.SELMODE.BUTTON.%d", mask));
    vTcl_DStringAppend(&input_params, "%s %s\n", name1, args.tag_list);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    if (SetActiveTags(args.tag_list) == -1) {
	return TCL_OK;
    }

    if (fij(io, mask, mode, args.min_overlap, (double)args.max_mis,
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

int
FindReadPairs(ClientData clientData,
		  Tcl_Interp *interp,
		  int argc,
		  char *argv[])
{
    contig_list_t *contig_array = NULL;
    int *contigs;
    int num_contigs = 0;
    Tcl_DString input_params;
    readpair_arg args;
    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(readpair_arg, io)},
	{"-contigs", ARG_STR, 1, NULL, offsetof(readpair_arg, inlist)},
	{NULL,	     0,	     0, NULL, 0}
    };

    vfuncheader("find read pairs");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /* create contig name array */
    active_list_contigs(args.io, args.inlist, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	xfree(contig_array);
	return TCL_OK;
    }
    contigs = to_contigs_only(num_contigs, contig_array);
    xfree(contig_array);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "Contigs: %s\n", args.inlist);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    if (find_read_pairs(args.io, num_contigs, contigs) < 0 ) {
	verror(ERR_WARN, "Find read pairs", "Failure in Find Read Pairs");
	return TCL_OK;
    }

    xfree(contigs);
    return TCL_OK;

} /* end FindReadPairs */


int
PlotQuality(ClientData clientData,
	    Tcl_Interp *interp,
	    int argc,
	    char *argv[])
{
    int handle;
    int i;

    vfuncheader("plot quality");

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be \"%.50s option "
		       "[arg arg ...]\"", argv[0]);
	return TCL_ERROR;
    }

    for (i = 1; i < argc; i++){
	if (strcmp(argv[i], "-io") == 0) {
	    handle = atoi(argv[i+1]);
	}
    }

    if (plot_quality(handle, consensus_cutoff) < 0 ) {
	verror(ERR_WARN, "Plot quality", "Failure in plot quality");
	return TCL_OK;
    }

    return TCL_OK;

} /* end PlotQuality */

/*
 * Generic assembly interface - clientData holds the task.
 *
 * 1 = Shotgun (normal and independent modes)
 * 2 = Screen only (set entry to 0)
 * 3 = One contig
 * 4 = New contigs
 * 5 = Single stranded
 * 6 = Independent
 */
int MainAssembly(ClientData clientData, Tcl_Interp *interp,
		 int argc, char *argv[]) {
    aa_arg args;
    int option;
    int entry;
    long mode = (long)clientData;
    Tcl_DString input_params;
    char *name1;
    char *name2;
    char *name3;
    char *name4;
    char *ret;

    cli_args a[] = {
/*1-6*/ {"-io",            ARG_INT,  1, NULL, offsetof(aa_arg, handle)},
/*1-6*/ {"-files",         ARG_STR,  1, NULL, offsetof(aa_arg, inlist)},
/*1256*/{"-output_mode",   ARG_INT,  1, "1",  offsetof(aa_arg, disp_mode)},
/*1256*/{"-min_match",     ARG_INT,  1, "20", offsetof(aa_arg, min_mat)},
/*1256*/{"-min_overlap",   ARG_INT,  1, "0",  offsetof(aa_arg, min_ovr)},
/*1256*/{"-max_pads",      ARG_INT,  1, "25", offsetof(aa_arg, max_pad)},
/*1256*/{"-max_pmismatch", ARG_FLOAT,1, "5.0",offsetof(aa_arg, max_mis)},
/*2*/   {"-save_align",    ARG_INT,  1, "0",  offsetof(aa_arg, align)},
/*156*/ {"-joins",         ARG_INT,  1, "1",  offsetof(aa_arg, joins)},
/*2*/   {"-win_size",      ARG_INT,  1, "0",  offsetof(aa_arg, win_size)},
/*2*/   {"-max_dashes",    ARG_INT,  1, "0",  offsetof(aa_arg, dash)},
/*16*/  {"-enter_failures",ARG_INT,  1, "0",  offsetof(aa_arg, fail_mode)},
/*126*/ {"-tag_types",     ARG_STR,  1, "",   offsetof(aa_arg, tag_list)},
/*16*/  {"-ignore_prev",   ARG_INT,  1, "0",  offsetof(aa_arg, ignore_prev)},
        {NULL,             0,       0,  NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("auto assemble");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (*args.tag_list != 0) {
	if (SetActiveTags(args.tag_list) == -1) {
	    return TCL_ERROR;
	}
    }

    /*
     * Set option and entry parameters
     */
    switch(mode) {
    case 1: /* shotgun */
	option = *args.tag_list ? 2 : 1;
	entry = 1;
	break;

    case 2: /* screen only */
	option = *args.tag_list ? 2 : 1;
	entry = 0;
	break;

    case 3: /* one contig */
    case 4: /* new contig */
    case 5: /* single stranded */
	option = mode;
	entry = 1;
	break;

    case 6: /* independent */
	option = *args.tag_list ? 2 : 1;
	entry = 1;
	args.ignore_prev = 1;
	break;

    default:
	return TCL_OK;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    if (mode == 1 || mode == 2) {
	if (*args.tag_list) {
	    vTcl_DStringAppend(&input_params, "Masking: %s\n", args.tag_list);
	} else {
	    Tcl_DStringAppend(&input_params, "Not using masking\n", -1);
	}
    }

    if (mode != 3 && mode != 4) {
	name1 = get_default_string(interp, gap_defs,
				   "AUTO_ASSEMBLE.MINMATCH.NAME");
	name2 = get_default_string(interp, gap_defs,
				   "AUTO_ASSEMBLE.MAXPADS.NAME");
	name3 = get_default_string(interp, gap_defs,
				   "AUTO_ASSEMBLE.MISMATCH.NAME");
	name4 = get_default_string(interp, gap_defs,
				   vw("AUTO_ASSEMBLE.DISPMODE.BUTTON.%d",
				      args.disp_mode));

	vTcl_DStringAppend(&input_params, "%s\n%s: %d\n%s: %d\n%s: %f\n",
			   name4,
			   name1, args.min_mat,
			   name2, args.max_pad,
			   name3, args.max_mis);
    }

    if (mode == 1) {
	if (args.joins) {
	    Tcl_DStringAppend(&input_params, "Permit joins\n", -1);
	} else {
	    Tcl_DStringAppend(&input_params, "Do not permit joins\n", -1);
	}

	name1 = get_default_string(interp, gap_defs,
				   vw("AUTO_ASSEMBLE.FAILURE.BUTTON.%d",
				      args.fail_mode+1));
	Tcl_DStringAppend(&input_params, name1, -1);
    }

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    if ((ret = auto_assemble(args.handle, args.inlist, option,
                             entry, args.disp_mode, args.min_mat,
                             args.min_ovr, args.max_pad, args.max_mis,
                             1-args.align, args.joins, args.fail_mode+1,
			     args.win_size, args.dash, args.ignore_prev,
                             consensus_cutoff)) == NULL) {
 	verror(ERR_WARN, "Auto assemble", "Failure in Auto Assemble");
	SetActiveTags("");
	return TCL_OK;
    }
    SetActiveTags("");

    Tcl_SetResult(interp, ret, TCL_DYNAMIC);
    return TCL_OK;
}

int
DisReadings(ClientData clientData,
	    Tcl_Interp *interp,
	    int argc,
	    char *argv[])
{
    char *name2 = NULL;
    dis_reading_arg args;
    char **reads = NULL;
    int *rnums, i, j;
    int num_reads;
    cli_args a[] = {
	{"-io",	     ARG_IO ,  1, NULL, offsetof(dis_reading_arg, io)},
	{"-readings",ARG_STR,  1, NULL, offsetof(dis_reading_arg, list)},
	{"-move",    ARG_INT,  1, "1",  offsetof(dis_reading_arg, move)},
	{"-duplicate_tags",
	             ARG_INT,  1, "1",  offsetof(dis_reading_arg,
						 duplicate_tags)},
	{NULL,	     0,	       0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("disassemble readings");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    name2 = get_default_string(interp, gap_defs,
			       vw("DIS_READINGS.SELTASK.BUTTON.%d",
				  args.move+1));

    if (name2) {
        Tcl_DString input_params;
        Tcl_DStringInit(&input_params);
	vTcl_DStringAppend(&input_params, "%s\n%s\n",
			   args.list, name2);
	vfuncparams("%s", Tcl_DStringValue(&input_params));
	Tcl_DStringFree(&input_params);
    }

    if (Tcl_SplitList(interp, args.list, &num_reads, &reads) != TCL_OK)
	return TCL_ERROR;

    if (NULL == (rnums = (int *)xmalloc(num_reads * sizeof(int))))
	return TCL_ERROR;

    for (i = j = 0; i < num_reads; i++) {
	rnums[j] = get_gel_num(args.io, reads[i], GGN_ID);
	if (rnums[j])
	    j++;
    }
    num_reads = j;

    if (disassemble_readings(args.io, rnums, num_reads, args.move,
			     args.duplicate_tags) < 0) {
	verror(ERR_WARN, "Disassemble readings",
	       "Failure in Disassemble Readings");
	return TCL_OK;
    }
    Tcl_Free((char *)reads);
    xfree(rnums);

    db_check(args.io);

    return TCL_OK;
} /* end DisReadings */

int
FindRepeats(ClientData clientData,
	    Tcl_Interp *interp,
	    int argc,
	    char *argv[])
{
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    find_repeats_arg args;
    GapIO *io;
    Tcl_DString input_params;
    char *name1;
    char *name2;
    char *name4;
    int mask;

    cli_args a[] = {
	{"-io",	       ARG_INT, 1, NULL, offsetof(find_repeats_arg, handle)},
	{"-direction", ARG_INT, 1, "3",  offsetof(find_repeats_arg, idir)},
	{"-min_match", ARG_INT, 1, "25", offsetof(find_repeats_arg, minmat)},
	{"-contigs",   ARG_STR, 1, NULL, offsetof(find_repeats_arg, inlist)},
	{"-outfile",   ARG_STR, 1, "",	 offsetof(find_repeats_arg, outfile)},
	{"-tag_types", ARG_STR,	1, "",   offsetof(find_repeats_arg, tag_list)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("find repeats");
    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    if ( (io = io_handle(&args.handle)) == NULL){
	verror(ERR_FATAL, "find_repeats", "invalid io handle");
	return -1;
    }

    active_list_contigs(io, args.inlist, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }

    mask = *args.tag_list ? 1 : 0;

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "Contigs: %s\n", args.inlist);

    name1 = get_default_string(interp, gap_defs, "FINDREP.MINREP.NAME");
    name2 = get_default_string(interp, gap_defs,
			       vw("FINDREP.SELTASK.BUTTON.%d",
				  args.idir));

    if (mask) {
	name4 = get_default_string(interp, gap_defs,
				   "FINDREP.SELMODE.BUTTON.1");
    } else {
	name4 = get_default_string(interp, gap_defs,
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

    if (find_repeats(io, args.handle, args.idir, args.minmat,
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

int
BreakContig(ClientData clientData,
	    Tcl_Interp *interp,
	    int argc,
	    char *argv[])
{
    break_contig_arg args;
    int *readings;
    int num_reads;
    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(break_contig_arg, io)},
	{"-readings",ARG_STR, 1, NULL, offsetof(break_contig_arg, reading)},
	{NULL,	 0,	  0, NULL, 0}
    };
    int i;
    int err = 0;

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("break contig");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    vfuncparams("Readings to be left ends of contigs: %s\n", args.reading);

    active_list_readings(args.io, args.reading, &num_reads, &readings);
    if (num_reads == 0) {
	if (readings) xfree(readings);
	return TCL_OK;
    }

    for (i = 0; i < num_reads; i++) {
	if (break_contig(args.io, readings[i]) != 0) {
	    Tcl_SetResult(interp, "Failure in Break Contig", TCL_STATIC);
	    err = 1;
	}
    }

    xfree(readings);

    return err ? TCL_ERROR : TCL_OK;
} /* end BreakContig */


/*
 * contig selector commands
 */
int
DisplayContigSelector(ClientData clientData,
		      Tcl_Interp *interp,
		      int argc,
		      char *argv[])
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

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;


    tag = tag_struct(interp, gap_defs, "CONTIG_SEL", args.tag_wd,
		     args.tag_offset);
    cursor = cursor_struct(interp, gap_defs, "CONTIG_SEL", args.cursor_wd,
			   args.cursor_fill);
    tick = tick_struct(interp, gap_defs, "CONTIG_SEL", args.tick_wd,
		       args.tick_ht, args.tick_fill);

    vTcl_SetResult(interp, "%d",
		   contig_selector_reg(interp, args.io, args.frame,
				       args.window, tag, cursor, tick));

    return TCL_OK;
}

int
DisplayContigComparator(ClientData clientData,
			Tcl_Interp *interp,
			int argc,
			char *argv[])
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

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    cs = result_data(args.io, args.id, 0);

    vTcl_SetResult(interp, "%d",
		   contig_comparator_reg(interp, args.io, cs, args.window,
					 args.v_window));

    return TCL_OK;
}

int
UpdateContigOrder(ClientData clientData,
		  Tcl_Interp *interp,
		  int argc,
		  char *argv[])
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

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }
    contigs = to_contigs_only(num_contigs, contig_array);
    xfree(contig_array);

    update_contig_order(interp, args.io, args.id, contigs, num_contigs,
			args.x);

    xfree(contigs);
    return TCL_OK;

}

int
FlushContigOrder(ClientData clientData,
		 Tcl_Interp *interp,
		 int argc,
		 char *argv[])
{
    io_arg args;
    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(io_arg, io)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ArrayDelay(args.io, args.io->db.contig_order, args.io->db.Ncontigs,
	       args.io->contig_order);
    flush2t(args.io);
    return TCL_OK;
}

int
DisplayCSTags(ClientData clientData,
	      Tcl_Interp *interp,
	      int argc,
	      char *argv[])
{
    reg_anno ra;
    cs_tags_arg args;
    cli_args a[] = {
	{"-io",	ARG_IO,  1, NULL, offsetof(cs_tags_arg, io)},
	{"-id", ARG_INT, 1, NULL, offsetof(cs_tags_arg, id)},
	{NULL,	0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ra.job = REG_ANNO;
    result_notify(args.io, args.id, (reg_data *)&ra, 0);
    return TCL_OK;
}

int
DisplayCSDiagonal(ClientData clientData,
		  Tcl_Interp *interp,
		  int argc,
		  char *argv[])
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

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    cs = result_data(args.io, args.id, 0);

    t_contig_length = CalcTotalContigLen(args.io);
    sprintf(cmd, "%s create line 1 1 %d %d -tag diagonal",
	    cs->window, t_contig_length, t_contig_length);
    Tcl_Eval(interp, cmd);

    scaleSingleCanvas(interp, cs->world, cs->canvas, cs->window, 'b',
		      "diagonal");
    return TCL_OK;
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
FindLongGels(ClientData clientData,
	     Tcl_Interp *interp,
	     int argc,
	     char *argv[])
{
    Tcl_DString input_params;
    long_gels_arg args;
    char *name2;
    int rargc;
    contig_list_t *rargv;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL,  offsetof(long_gels_arg, io)},
	{"-contigs", ARG_STR, 1, NULL,  offsetof(long_gels_arg, contigs)},
	{"-avg_len", ARG_INT, 1, "500", offsetof(long_gels_arg, avg_len)},
	{NULL,	  0,	   0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("suggest long readings");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &rargc, &rargv);
    if (rargc == 0) {
	xfree(rargv);
	return TCL_OK;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    name2 = get_default_string(interp, gap_defs, "LONGGELS.GLEN.NAME");
    vTcl_DStringAppend(&input_params, "Contigs: %s\n%s: %d\n",
		       args.contigs, name2, args.avg_len);
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    (void)find_long_gels(args.io, rargc, rargv, args.avg_len);

    xfree(rargv);
    return TCL_OK;
}

int
FindTaqTerminator(ClientData clientData,
		   Tcl_Interp *interp,
		   int argc,
		   char *argv[])
{
    Tcl_DString input_params;
    char *name2;
    taq_terms_arg args;
    int rargc;
    contig_list_t *rargv;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL,  offsetof(taq_terms_arg, io)},
	{"-contigs", ARG_STR, 1, NULL,  offsetof(taq_terms_arg, contigs)},
	{"-avg_len", ARG_INT, 1, "350", offsetof(taq_terms_arg, avg_len)},
	{NULL,	  0,	   0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("solve compressions and stops");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &rargc, &rargv);
    if (rargc == 0) {
	xfree(rargv);
	return TCL_OK;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    name2 = get_default_string(interp, gap_defs, "TTERM.TLEN.NAME");

    vTcl_DStringAppend(&input_params, "Contigs: %s\n%s: %d\n",
		       args.contigs, name2, args.avg_len);
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    (void)find_taq_terms(args.io, rargc, rargv, args.avg_len);

    xfree(rargv);
    return TCL_OK;

}

int
CheckDatabase(ClientData clientData,
	      Tcl_Interp *interp,
	      int argc,
	      char *argv[])
{
    io_arg args;

    cli_args a[] = {
	{"-io",	  ARG_IO,  1, NULL, offsetof(io_arg, io)},
	{NULL,	  0,	   0, NULL, 0}
    };

    vfuncheader("check database");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    vTcl_SetResult(interp, "%d", db_check(args.io));
    return TCL_OK;
}

int
tcl_remove_contig_holes(ClientData clientData,
			Tcl_Interp *interp,
			int argc,
			char *argv[])
{
    list2_arg args;
    int rargc, i;
    contig_list_t *rargv;

    cli_args a[] = {
	{"-io",	      ARG_IO,  1, NULL, offsetof(list2_arg, io)},
	{"-contigs",  ARG_STR, 1, NULL, offsetof(list2_arg, inlist)},
	{NULL,	  0,	   0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("remove_contig_holes");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    /*
     * NB this may create, and shuffle contig numbers, so we convert to
     * left-most reading number instead, and swap back to contig number
     * before each call of remove_contig_holes().
     */
    for (i = 0; i < rargc; i++) {
	rargv[i].contig = io_clnbr(args.io, rargv[i].contig);
    }

    for (i = 0; i < rargc; i++) {
	int cnum = rnumtocnum(args.io, rargv[i].contig);
	remove_contig_holes(args.io, cnum);
    }

    xfree(rargv);
    return TCL_OK;
}

int
DeleteContig(ClientData clientData,
	     Tcl_Interp *interp,
	     int argc,
	     char *argv[])
{
    delete_contig_arg args;
    char **contigs;
    int num_contigs;
    int i, contig;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(delete_contig_arg, io)},
	{"-contigs", ARG_STR, 1, NULL, offsetof(delete_contig_arg, contig)},
	{NULL,	  0,	   0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("delete contig");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /* Care must be taken as this changes contig numbers. Delay conversion */
    if (Tcl_SplitList(interp, args.contig, &num_contigs, &contigs) != TCL_OK)
	return TCL_ERROR;

    for (i = 0; i < num_contigs; i++) {
	if (-1 == (contig = get_contig_num(args.io, contigs[i], GGN_ID))) {
	    verror(ERR_WARN, "delete_contig", "unknown contig %s",
		   contigs[i]);
	} else {
	    delete_contig(args.io, contig);
	}
    }

    Tcl_Free((char *)contigs);
    return TCL_OK;
}


int
AnnotationAddress(ClientData clientData,
	     Tcl_Interp *interp,
	     int argc,
	     char *argv[])
{
    ann_addr_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,	 1, NULL, offsetof(ann_addr_arg, io)},
	{"-annotation", ARG_INT, 1, NULL, offsetof(ann_addr_arg, anno)},
	{NULL,	  0,	   0, NULL, 0}
    };
    int type, num, mode = 1, c;
    char buf[100];

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    Tcl_ResetResult(interp);
    while ((c = annotation_address(args.io, mode, args.anno, &type, &num))>0) {
	sprintf(buf, "{%d %d %d} ", type, num, c);
	Tcl_AppendResult(interp, buf, NULL);
	mode = 0;
    }

    if (mode == 1) {
	Tcl_AppendResult(interp, "{}", NULL);
    }

    annotation_address(args.io, 2, 0, NULL, NULL);

    return TCL_OK;
}

int
ResetContigOrder(ClientData clientData,
		 Tcl_Interp *interp,
		 int argc,
		 char *argv[])
{
    io_arg args;

    cli_args a[] = {
	{"-io",	  ARG_IO,  1, NULL, offsetof(io_arg, io)},
	{NULL,	  0,	   0, NULL, 0}
    };

    vfuncheader("reset contig order");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    vTcl_SetResult(interp, "%d", reset_contig_order(args.io));
    return TCL_OK;
}

int
DoubleStrand(ClientData clientData,
	     Tcl_Interp *interp,
	     int argc,
	     char *argv[])
{
    Tcl_DString input_params;
    char *name1;
    char *name2;
    int rargc;
    contig_list_t *rargv;
    dstrand_arg args;
    cli_args a[] = {
	{"-io",	           ARG_IO,  1, NULL, offsetof(dstrand_arg, io)},
	{"-max_nmismatch", ARG_INT, 1, "8",  offsetof(dstrand_arg, maxmis)},
	{"-max_pmismatch", ARG_FLOAT,1, "8.", offsetof(dstrand_arg, maxperc)},
	{"-contigs",       ARG_STR, 1, NULL, offsetof(dstrand_arg, list)},
	{NULL,	  0,	   0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("double strand");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.list, &rargc, &rargv);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "Contigs: %s\n", args.list);

    name1 = get_default_string(interp, gap_defs, "DOUBLE_STRAND.MAXMIS.NAME");
    name2 = get_default_string(interp, gap_defs, "DOUBLE_STRAND.MAXPERC.NAME");

    vTcl_DStringAppend(&input_params, "%s: %d\n%s: %f\n",
		       name1, args.maxmis, name2, args.maxperc);
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    double_strand_list(args.io, rargc, rargv, args.maxmis, args.maxperc);
    xfree(rargv);

    return TCL_OK;
}

int
FindPrimers(ClientData clientData,
	       Tcl_Interp *interp,
	       int argc,
	       char *argv[])
{
    int rargc;
    contig_list_t *rargv;
    primer_arg args;
    cli_args a[] = {
	{"-io",	         ARG_IO,  1, NULL, offsetof(primer_arg, io)},
	{"-contigs",     ARG_STR, 1, NULL, offsetof(primer_arg, inlist)},
	{"-search_from", ARG_INT, 1, "20", offsetof(primer_arg, search_from)},
	{"-search_to",	 ARG_INT, 1, "60", offsetof(primer_arg, search_to)},
	{"-num_primers", ARG_INT, 1, "1",  offsetof(primer_arg, num_primers)},
	{"-primer_start",ARG_INT, 1, "1",  offsetof(primer_arg, primer_start)},
	{"-params",      ARG_STR, 1, "",   offsetof(primer_arg, primer_defs)},
	{NULL,	         0,	  0, NULL, 0}
    };
    char *res, *primer_defs;

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("suggest primers");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (*args.primer_defs == 0) {
	primer_defs = get_default_string(interp, gap_defs, "PRIMER");
    } else {
	primer_defs = args.primer_defs;
    }

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    res = suggest_primers_list(args.io, rargc, rargv,
			       args.search_from, args.search_to,
			       args.num_primers, args.primer_start,
			       primer_defs);
    xfree(rargv);

    Tcl_SetResult(interp, res, TCL_DYNAMIC);

    return TCL_OK;
}

int
ExtractReadings(ClientData clientData,
		Tcl_Interp *interp,
		int argc,
		char *argv[])
{
    extract_arg args;
    char **reading_array = NULL;
    int num_readings;

    cli_args a[] = {
	{"-io",	      ARG_IO,  1, NULL,       offsetof(extract_arg, io)},
	{"-readings", ARG_STR, 1, NULL,       offsetof(extract_arg, list)},
	{"-directory",ARG_STR, 1, "extracts", offsetof(extract_arg, dir)},
	{"-format",   ARG_INT, 1, "0",        offsetof(extract_arg, format)},
	{NULL,	0,	0, NULL, 0}
    };

    vfuncheader("extract readings");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /* create reading name array */
    if (Tcl_SplitList(interp, args.list, &num_readings, &reading_array)
	!= TCL_OK) {
	return TCL_ERROR;
    }

    (void)extract_readings(args.io, num_readings, reading_array, args.dir,
			   args.format);
    Tcl_Free((char *)reading_array);

    return TCL_OK;
}

int
PreAssemble(ClientData clientData,
	    Tcl_Interp *interp,
	    int argc,
	    char *argv[])
{
    pre_ass_arg args;
    char **reading_array = NULL;
    int num_readings;

    cli_args a[] = {
	{"-io",	    ARG_INT, 1, NULL, offsetof(pre_ass_arg, handle)},
	{"-files",  ARG_STR, 1, NULL, offsetof(pre_ass_arg, list)},
	{NULL,	0,	0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("enter preassembled data");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /* create reading name array */
    if (Tcl_SplitList(interp, args.list, &num_readings, &reading_array)
	!= TCL_OK) {
	return TCL_ERROR;
    }

    pre_assemble(args.handle, num_readings, reading_array);
    Tcl_Free((char *)reading_array);
    return TCL_OK;

}

int
Consensus(ClientData clientData,
	  Tcl_Interp *interp,
	  int argc,
	  char *argv[])
{
    consensus_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    int type, mask;

    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(consensus_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,  offsetof(consensus_arg, inlist)},
	{"-type",	ARG_STR, 1, "normal",offsetof(consensus_arg, type)},
	{"-mask",	ARG_STR, 1, "none",offsetof(consensus_arg, mask)},
	{"-win_size",	ARG_INT, 1, "0",   offsetof(consensus_arg, win_size)},
	{"-max_dashes",	ARG_INT, 1, "0",   offsetof(consensus_arg, dash)},
	{"-format",	ARG_INT, 1, "3",   offsetof(consensus_arg, format)},
	{"-annotations",ARG_INT, 1, "0",   offsetof(consensus_arg, gel_anno)},
	{"-notes",      ARG_INT, 1, "0",   offsetof(consensus_arg, gel_notes)},
	{"-truncate",	ARG_INT, 1, "0",   offsetof(consensus_arg, truncate)},
	{"-outfile",	ARG_STR, 1, NULL,  offsetof(consensus_arg, out_file)},
	{"-tag_types",	ARG_STR, 1, "",    offsetof(consensus_arg, tag_list)},
	{"-strip_pads",	ARG_INT, 1, "1",   offsetof(consensus_arg, nopads)},
	{"-min_conf",   ARG_INT, 1, "8",   offsetof(consensus_arg, min_conf)},
	{"-use_conf",	ARG_INT, 1, "1",   offsetof(consensus_arg, use_conf)},
	{"-name_format",ARG_INT, 1, "1",   offsetof(consensus_arg, name_format)},
	{NULL,	        0,	 0, NULL,  0}
    };

    vfuncheader("calculate a consensus");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (strcmp(args.type, "normal") == 0)
	type = 1;
    else if (strcmp(args.type, "extended") == 0)
	type = 2;
    else if (strcmp(args.type, "unfinished") == 0)
	type = 3;
    else if (strcmp(args.type, "quality") == 0)
	type = 4;
    else {
	Tcl_SetResult(interp, "Unknown type", TCL_STATIC);
	return TCL_ERROR;
    }

    /* Parse mode and mask */
    if (strcmp(args.mask, "none") == 0)
	mask = 0;
    else if (strcmp(args.mask, "mark") == 0)
	mask = 2;
    else if (strcmp(args.mask, "mask") == 0)
	mask = 1;
    else {
	Tcl_SetResult(interp, "invalid mask mode", TCL_STATIC);
	return TCL_ERROR;
    }

    active_list_contigs(args.io, args.inlist, &num_contigs, &contig_array);

    /* set up active tag list */
    if (*args.tag_list != 0) {
	if (SetActiveTags(args.tag_list) == -1) {
	    return TCL_ERROR;
	}
    }
    if (-1 == consensus_dialog(args.io, mask, type, args.format,
			       args.gel_anno, args.truncate,
			       args.gel_notes,
			       args.use_conf, args.min_conf,
			       args.win_size, args.dash,
			       args.out_file, num_contigs,
			       contig_array, args.nopads,
			       args.name_format)) {
	verror(ERR_WARN, "consensus_ouput",
	       "failed to calculate or write file");
    }

    if (contig_array)
	free(contig_array);
    SetActiveTags("");

    return TCL_OK;
}

int tcl_calc_consensus(ClientData clientData, Tcl_Interp *interp,
		       int argc, char **argv)
{
    int rargc;
    contig_list_t *rargv;
    list2_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(list2_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL,  offsetof(list2_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    if (rargc >= 1) {
	char *buf;

	if (NULL == (buf = (char *)xmalloc(rargv[0].end - rargv[0].start + 2)))
	    return TCL_OK;

	calc_consensus(rargv[0].contig, rargv[0].start, rargv[0].end,
		       CON_SUM, buf, NULL, NULL, NULL,
		       consensus_cutoff, quality_cutoff,
		       database_info, (void *)args.io);
	buf[rargv[0].end - rargv[0].start + 1] = 0;
	Tcl_SetResult(interp, buf, TCL_DYNAMIC);
    }

    xfree(rargv);
    return TCL_OK;
}

/*int tcl_calc_quality(ClientData clientData, Tcl_Interp *interp,
		       int argc, char **argv)
*/
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
    Tcl_Obj *obj;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    if (rargc >= 1) {
	char *con;
	float *qual;
	char *qualc;
	int i;

	qual = (float *)xmalloc((rargv[0].end - rargv[0].start + 2) *
				sizeof(*qual));
	con = (char *)xmalloc((rargv[0].end - rargv[0].start + 2) *
				sizeof(*con));
	qualc = (char *)xmalloc((rargv[0].end - rargv[0].start + 2) *
				sizeof(*qualc));

	if (!qual || !con || !qualc)
	    return TCL_OK;

	calc_consensus(rargv[0].contig, rargv[0].start, rargv[0].end,
		       CON_SUM, con, NULL, qual, NULL,
		       consensus_cutoff, quality_cutoff,
		       database_info, (void *)args.io);

	for (i = 0; i < rargv[0].end - rargv[0].start + 1; i++) {
	    qualc[i] = (char)(qual[i] + 0.499);
	}

	obj = Tcl_NewStringObj(qualc, rargv[0].end - rargv[0].start + 1);
	Tcl_SetObjResult(interp, obj);

	xfree(qual);
	xfree(con);
	xfree(qualc);
    }

    xfree(rargv);
    return TCL_OK;
}

int tcl_calc_consensus_double(ClientData clientData, Tcl_Interp *interp,
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
    Tcl_Obj *lobj, *obj;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);

    if (rargc >= 1) {
	char *con, *con1, *con2;
	float *qual, *qual1, *qual2;
	int i;
	Tcl_Obj *lobj2;

	if (NULL == (lobj = Tcl_NewListObj(0, NULL)))
	    return TCL_ERROR;
	Tcl_IncrRefCount(lobj);

	qual = (float *)xmalloc((rargv[0].end - rargv[0].start + 2) *
				sizeof(*qual));
	con = (char *)xmalloc((rargv[0].end - rargv[0].start + 2) *
				sizeof(*con));
	qual1 = (float *)xmalloc((rargv[0].end - rargv[0].start + 2) *
				sizeof(*qual1));
	con1 = (char *)xmalloc((rargv[0].end - rargv[0].start + 2) *
				sizeof(*con1));
	qual2 = (float *)xmalloc((rargv[0].end - rargv[0].start + 2) *
				sizeof(*qual2));
	con2 = (char *)xmalloc((rargv[0].end - rargv[0].start + 2) *
				sizeof(*con2));

	if (!qual || !con || !qual1 || !qual2 || !con1 || !con2)
	    return TCL_OK;

	calc_consensus(rargv[0].contig, rargv[0].start, rargv[0].end,
		       CON_SUM, con, NULL, qual, NULL,
		       consensus_cutoff, quality_cutoff,
		       database_info, (void *)args.io);

	calc_consensus(rargv[0].contig, rargv[0].start, rargv[0].end,
		       CON_SUM, con1, con2, qual1, qual2,
		       consensus_cutoff, quality_cutoff,
		       database_info, (void *)args.io);

	for (i = 0; i < rargv[0].end - rargv[0].start + 1; i++) {
	    /* list of con, qual, con+, qual+, con-, con- */
	    lobj2 = Tcl_NewListObj(0, NULL);

	    /* con */
	    obj = Tcl_NewStringObj(&con[i], 1);
	    Tcl_ListObjAppendElement(interp, lobj2, obj);

	    /* qual */
	    obj = Tcl_NewDoubleObj(qual[i]);
	    Tcl_ListObjAppendElement(interp, lobj2, obj);

	    /* con+ */
	    obj = Tcl_NewStringObj(&con1[i], 1);
	    Tcl_ListObjAppendElement(interp, lobj2, obj);

	    /* qual+ */
	    obj = Tcl_NewDoubleObj(qual1[i]);
	    Tcl_ListObjAppendElement(interp, lobj2, obj);

	    /* con- */
	    obj = Tcl_NewStringObj(&con2[i], 1);
	    Tcl_ListObjAppendElement(interp, lobj2, obj);

	    /* qual- */
	    obj = Tcl_NewDoubleObj(qual2[i]);
	    Tcl_ListObjAppendElement(interp, lobj2, obj);

	    Tcl_ListObjAppendElement(interp, lobj, lobj2);
	}

	Tcl_SetObjResult(interp, lobj);
	Tcl_DecrRefCount(lobj);

	xfree(qual);
	xfree(con);
	xfree(qual1);
	xfree(con1);
	xfree(qual2);
	xfree(con2);
    }

    xfree(rargv);
    return TCL_OK;
}

int EnterTags(ClientData clientData,
	      Tcl_Interp *interp,
	      int argc,
	      char *argv[])
{
    f_int ngels, nconts, idbsiz;
    enter_tags_arg args;
    GapIO *io;
    cli_args a[] = {
	{"-io",	      ARG_IO,  1, NULL, offsetof(enter_tags_arg, io)},
	{"-file",     ARG_STR, 1, NULL, offsetof(enter_tags_arg, file)},
	{"-unpadded", ARG_INT, 1, "0",  offsetof(enter_tags_arg, unpadded)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("enter tags");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    io = args.io;

    ngels  = NumReadings(io);
    nconts = NumContigs(io);
    idbsiz = io_dbsize(io);
    tagfil_(&io_relpos(io,1), &io_length(io,1), &io_lnbr(io,1), &io_rnbr(io,1),
	    &ngels, &nconts, &idbsiz, args.file, handle_io(io), NULL,
	    &args.unpadded, strlen(args.file));

    return TCL_OK;
}

int
check_active_lock(ClientData clientData,
		  Tcl_Interp *interp,
		  int argc,
		  char *argv[])
{
    vTcl_SetResult(interp, "%d", semaphoreFree(activeLock));
    return TCL_OK;
}




int
tcl_assemble_direct(ClientData clientData, Tcl_Interp *interp,
		    int argc, char *argv[]) {
    Tcl_DString input_params;
    char *name1, *res;
    ass_direct_arg args;
    cli_args a[] = {
	{"-io",	         ARG_IO,   1, NULL, offsetof(ass_direct_arg, io)},
	{"-files",       ARG_STR,  1, NULL, offsetof(ass_direct_arg, inlist)},
	{"-output_mode", ARG_INT,  1, "0",  offsetof(ass_direct_arg, display)},
	{"-max_pmismatch",ARG_FLOAT,1,"-1.",offsetof(ass_direct_arg, mism)},
	{"-align",       ARG_INT,  1, "1",  offsetof(ass_direct_arg, align)},
	{"-enter_failures",ARG_INT, 1, "0", offsetof(ass_direct_arg,
						     enter_failures)},
	{NULL,	     0,	       0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("directed assembly");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    Tcl_DStringInit(&input_params);
    if (args.display) {
	Tcl_DStringAppend(&input_params, "Display alignments\n", -1);
    } else {
	Tcl_DStringAppend(&input_params, "Do not display alignments\n", -1);
    }

    name1 = get_default_string(interp, gap_defs,
			       "DIRECT_ASSEMBLY.MAXMIS.NAME");

    vTcl_DStringAppend(&input_params, "%s: %f\n", name1, args.mism);
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    res = assemble_direct(args.io, args.display, (double)args.mism,
			  args.inlist, args.align, args.enter_failures);
    vTcl_SetResult(interp, "%s", res ? res : "");
    xfree(res);

    return TCL_OK;
}

int
tcl_check_assembly(ClientData clientData, Tcl_Interp *interp,
		       int argc, char *argv[]) {
    int rargc, *contigs;
    contig_list_t *rargv;
    Tcl_DString input_params;
    char *name1;
    char *name2;
    char *name3;
    check_ass_arg args;

    cli_args a[] = {
	{"-io",	        ARG_IO,	 1, NULL, offsetof(check_ass_arg, io)},
	{"-contigs",    ARG_STR, 1, NULL, offsetof(check_ass_arg, inlist)},
	{"-cutoff",     ARG_INT, 1, "1",  offsetof(check_ass_arg, cutoff)},
	{"-min_len",    ARG_INT, 1, "10", offsetof(check_ass_arg, min_len)},
	{"-win_size",   ARG_INT, 1, "29", offsetof(check_ass_arg, win_size)},
	{"-max_dashes", ARG_INT, 1, "3",  offsetof(check_ass_arg, max_dash)},
	{"-max_pmismatch", ARG_FLOAT, 1, "15.0",
	     offsetof(check_ass_arg, max_mismatch)},
	{NULL,	  0,	   0, NULL, 0}
    };

    vfuncheader("check assembly");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &rargc, &rargv);
    if (rargc == 0) {
	xfree(rargv);
	return TCL_OK;
    }
    contigs = to_contigs_only(rargc, rargv);
    xfree(rargv);

    name1 = get_default_string(interp, gap_defs,
			       "CHECK_ASSEMBLY.MAXPERC.NAME");

    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "Contigs %s\n%s: %f\n",
		       args.inlist, name1, args.max_mismatch);

    if (args.cutoff) {
	name1 = get_default_string(interp, gap_defs,
				   "CHECK_ASSEMBLY.MINLEN.NAME");
	name2 = get_default_string(interp, gap_defs,
				   "CHECK_ASSEMBLY.HIDDEN.WINSIZE.NAME");
	name3 = get_default_string(interp, gap_defs,
				   "CHECK_ASSEMBLY.HIDDEN.MAXDASH.NAME");
	vTcl_DStringAppend(&input_params,
			   "Hidden data: %s: %d\n%s: %d\n%s: %d\n",
			   name1, args.min_len,
			   name2, args.win_size,
			   name3, args.max_dash);
    } else {
	Tcl_DStringAppend(&input_params, "Not using hidden data\n", -1);
    }
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);
    check_assembly(args.io, rargc, contigs, args.cutoff, args.min_len,
		   args.win_size, args.max_dash, args.max_mismatch / 100.0);

    xfree(contigs);

    return TCL_OK;
}


int tcl_anno_list(ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[]) {
    int i;
    anno_list_arg args;
    cli_args a[] = {
	{"-io",	      ARG_IO,	 1, NULL, offsetof(anno_list_arg, io)},
	{"-type",     ARG_STR,	 1, NULL, offsetof(anno_list_arg, type)},
	{NULL,	      0,	 0, NULL, 0}
    };
    Array l;

    vfuncheader("output annotations");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (l = anno_list(args.io, str2type(args.type)))) {
	verror(ERR_FATAL, "tcl_anno_list", "out of memory");
	return TCL_OK;
    }

    Tcl_ResetResult(interp);
    for (i = 0; i < ArrayMax(l); i++) {
	char type[5];
	struct anno_list_t *p;
	char buf[1024];

	p = arrp(struct anno_list_t, l, i);
	sprintf(buf, "%d %s %d %d %d\n",
		p->anno,
		type2str(p->type, type),
		p->position,
		p->length,
		p->strand);

	Tcl_AppendResult(interp, buf, NULL);
    }

    ArrayDestroy(l);

    return TCL_OK;
}


int tcl_delete_anno_list(ClientData clientData, Tcl_Interp *interp,
			 int argc, char *argv[]) {
    delete_anno_list_arg args;
    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(delete_anno_list_arg, io)},
	{"-annos",  ARG_STR, 1, NULL, offsetof(delete_anno_list_arg, annos)},
	{NULL,	    0,	     0, NULL, 0}
    };
    char *cp;
    int anno, pos, count = 0, *anno_av;

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("delete annotations");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /* Count annos */
    cp = args.annos;
    while (1 == sscanf(cp, "%d %*s %*d %*d %*d\n%n", &anno, &pos)) {
	cp += pos;
	count++;
    }

    if (0 == count)
	return TCL_OK;

    /* Alloc and copy anno numbers */
    if (NULL == (anno_av = (int *)xmalloc(count * sizeof(int))))
	return TCL_OK;

    count = 0;
    cp = args.annos;
    while (1 == sscanf(cp, "%d %*s %*d %*d %*d\n%n", &anno, &pos)) {
	cp += pos;
	anno_av[count++] = anno;
    }

    if (-1 == rmanno_list(args.io, count, anno_av))
	verror(ERR_FATAL, "delete_annotations", "out of memory");

    return TCL_OK;
}

int
tcl_find_probes(ClientData clientData, Tcl_Interp *interp,
		    int argc, char *argv[]) {
    Tcl_DString input_params;
    char *name1;
    char *name2;
    char *name3;
    char *name4;
    char *name5;
    char *name6;
    int rargc, *contigs;
    contig_list_t *rargv;
    find_probes_arg args;
    cli_args a[] = {
	{"-io",	       ARG_IO,  1, NULL, offsetof(find_probes_arg, io)},
	{"-contigs",   ARG_STR, 1, NULL, offsetof(find_probes_arg, contigs)},
	{"-min_size",  ARG_INT, 1, "15", offsetof(find_probes_arg, min_size)},
	{"-max_size",  ARG_INT, 1, "19", offsetof(find_probes_arg, max_size)},
	{"-max_pmatch",ARG_FLOAT,1,"90.0",offsetof(find_probes_arg, max_perc)},
	{"-from",      ARG_INT, 1, "10", offsetof(find_probes_arg, from)},
	{"-to",	       ARG_INT, 1, "100",offsetof(find_probes_arg, to)},
	{"-vectors",   ARG_STR, 1, ""   ,offsetof(find_probes_arg, vectors)},
	{"-primer_arg",ARG_STR, 1, ""   ,offsetof(find_probes_arg, primer_arg)},
	{NULL,	      0,       0, NULL, 0}
    };
    Tcl_DString dstr;

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("suggest probes");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (args.vectors && args.vectors[0] == '\0')
	args.vectors = NULL;

    active_list_contigs(args.io, args.contigs, &rargc, &rargv);
    if (rargc == 0) {
	xfree(rargv);
	return TCL_OK;
    }
    contigs = to_contigs_only(rargc, rargv);
    xfree(rargv);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "Contigs: %s\n", args.contigs);

    name1 = get_default_string(interp, gap_defs,
			       "FIND_PROBES.MAXPERC.NAME");
    name2 = get_default_string(interp, gap_defs,
			       "FIND_PROBES.OLIGOSIZE.MIN_NAME");
    name3 = get_default_string(interp, gap_defs,
			       "FIND_PROBES.OLIGOSIZE.MAX_NAME");
    name4 = get_default_string(interp, gap_defs,
			       "FIND_PROBES.OLIGOPOS.MIN_NAME");
    name5 = get_default_string(interp, gap_defs,
			       "FIND_PROBES.OLIGOPOS.MAX_NAME");
    name6 = get_default_string(interp, gap_defs,
			       "FIND_PROBES.VECTOR.NAME");

    vTcl_DStringAppend(&input_params,
		       "%s: %f\n%s: %d\n%s: %d\n%s: %d\n%s: %d\n%s: %s\n",
		       name1, args.max_perc,
		       name2, args.min_size,
		       name3, args.max_size,
		       name4, args.from,
		       name5, args.to,
		       name6, args.vectors ? args.vectors : "<none>");
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);

    Tcl_DStringInit(&dstr);

    /* Find oligos for probing */
    if (-1 == find_probes(args.io, rargc, contigs,
			  args.min_size, args.max_size,
			  args.max_perc / 100.0, args.from, args.to,
			  args.vectors, args.primer_arg, &dstr)) {
	verror(ERR_WARN, "find_probes", "failed");
    }

    Tcl_DStringResult(interp, &dstr);

    xfree(contigs);

    return TCL_OK;
}

int
FindOligo(ClientData clientData,
	  Tcl_Interp *interp,
	  int argc,
	  char *argv[])
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
	{NULL,		0,	  0, NULL,  0}
    };

    vfuncheader("sequence search");

    if (-1 == gap_parse_args(a, &args, argc, argv))
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
    name1 = get_default_string(interp, gap_defs, "FINDOLIGO.MAXMIS.NAME");

    vTcl_DStringAppend(&input_params, "%s: %f\n", name1, args.mis_match);

    if (*args.seq) {
	vTcl_DStringAppend(&input_params, "Sequence: %s\n", args.seq);
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

    if (-1 == find_oligos(args.io, num_contigs, contig_array,
			  args.mis_match, args.seq,
			  args.consensus_only, args.cutoffs))
	verror(ERR_FATAL, "find oligos", "out of memory");

    SetActiveTags("");
    if (contig_array)
	xfree(contig_array);
    return TCL_OK;
}


int tcl_add_tags(ClientData clientData, Tcl_Interp *interp,
		int argc, char **argv) {
    add_tags_arg args;
    int num_tags, i;
    char **tags = NULL;
    int *contigs;
    reg_anno ra;
    reg_buffer_start rs;
    reg_buffer_end re;
    int last;
    int *cache = NULL;
    int cache_len, cache_pos;

    cli_args a[] = {
	{"-io",		ARG_IO,	 1, NULL, offsetof(add_tags_arg, io)},
	{"-tags",	ARG_STR, 1, NULL, offsetof(add_tags_arg, tag_list)},
	{"-unpadded",   ARG_INT, 1, "0",  offsetof(add_tags_arg, unpadded)},
	{NULL,		0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /* create reading name array */
    if (Tcl_SplitList(interp, args.tag_list, &num_tags, &tags)
	!= TCL_OK) {
	return TCL_ERROR;
    }

    /*
     * Allocate our contigs[] array, to store info on which contigs have been
     * modified.
     */
    if (NULL == (contigs = (int *)xcalloc(NumContigs(args.io), sizeof(int))))
	return TCL_ERROR;

    last = 0;
    for (i=0; i<num_tags; i++) {
	int r, n, gellen, cnum;

	sscanf(tags[i], "%d %n", &r, &n);
	if (r >= 0) {
	    GReadings rs;
	    gel_read(args.io, r, rs);
	    gellen = rs.length;
	    cnum = rnumtocnum(args.io, r);
	    if (contigs[cnum-1]&2)
		continue;
	    if (contig_lock_write(args.io, cnum) == -1) {
		verror(ERR_WARN, "add_tags", "Contig is busy");
		contigs[cnum-1]|=2;
		continue;
	    }
	    contigs[cnum-1]|=1;
	} else {
	    gellen = io_clength(args.io, -r);
	    if (contigs[-r-1]&2)
		continue;
	    if (contig_lock_write(args.io, -r) == -1) {
		verror(ERR_WARN, "add_tags", "Contig is busy");
		contigs[-r-1]|=2;
		continue;
	    }
	    contigs[-r-1]|=1;
	}
	if (r != last) {
	    if (cache)
		xfree(cache);
	    cache = (int *)xcalloc(gellen+2, sizeof(int));
	    cache_len = gellen;
	    cache_pos = 0;
	    last = r;
	}
	create_tag_for_gel(args.io, r, gellen, tags[i]+n,
			   cache, cache_len, &cache_pos, args.unpadded);
    }
    if (cache)
	xfree(cache);

    /* Notify of the start of the flurry of updates */
    rs.job = REG_BUFFER_START;
    for (i = 0; i < NumContigs(args.io); i++) {
	if (contigs[i]&1) {
	    contig_notify(args.io, i+1, (reg_data *)&rs);
	}
    }

    /* Now notify all the contigs that we've added tags to */
    ra.job = REG_ANNO;
    for (i = 0; i < NumContigs(args.io); i++) {
	if (contigs[i]&1) {
	    contig_notify(args.io, i+1, (reg_data *)&ra);
	}
    }

    /* Notify of the end of the flurry of updates */
    re.job = REG_BUFFER_END;
    for (i = 0; i < NumContigs(args.io); i++) {
	if (contigs[i]&1) {
	    contig_notify(args.io, i+1, (reg_data *)&re);
	}
    }

    flush2t(args.io);
    xfree(contigs);
    Tcl_Free((char *)tags);

    return TCL_OK;
}


int tcl_contig_order_to_number(ClientData clientData, Tcl_Interp *interp,
			       int argc, char **argv) {
    ord2num_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,	 1, NULL, offsetof(ord2num_arg, io)},
	{"-order",	ARG_INT, 1, NULL, offsetof(ord2num_arg, order)},
	{NULL,		0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (args.order > 0 && args.order <= NumContigs(args.io)) {
	vTcl_SetResult(interp, "%d",
		arr(GCardinal, args.io->contig_order, args.order-1));
	return TCL_OK;
    } else {
	Tcl_SetResult(interp, "Invalid contig number", TCL_STATIC);
	return TCL_ERROR;
    }
}

/*
 * Converts a list of reading identifiers into reading names.
 */
int tcl_get_read_names(ClientData clientData, Tcl_Interp *interp,
		       int argc, char **argv) {
    int i;
    GapIO *io;
    Tcl_DString result;

    if (argc < 3) {
	Tcl_SetResult(interp, "Wrong # args: "
		      "get_read_names -io io_id name ...\n",
		      TCL_STATIC);
	return TCL_ERROR;
    }

    if (strcmp(argv[1], "-io") != 0) {
	Tcl_SetResult(interp, "Usage: get_read_names -io io_id name ...\n",
		      TCL_STATIC);
	return TCL_ERROR;
    }

    i = atoi(argv[2]);
    if (NULL == (io = io_handle(&i))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    Tcl_DStringInit(&result);
    for (i=3; i<argc; i++) {
	int rnum;

	/*
	 * Both convert and validate ids.
	 * It's more efficient to check for '#' and only convert then numbers
	 * to names, but this data is cached in memory these days anyway and
	 * this also allows us to spot invalid reading names.
	 */
	if (-1 == (rnum = get_gel_num(io, argv[i], GGN_ID))) {
	    verror(ERR_WARN, "get_read_names",
		   "reading '%s' not found", argv[i]);
	    continue;
	}

	Tcl_DStringAppendElement(&result, get_read_name(io, rnum));
    }
    Tcl_DStringResult(interp, &result);

    return TCL_OK;
}

int
tcl_order_contigs(ClientData clientData,
	     Tcl_Interp *interp,
	     int argc,
	     char *argv[])
{
    io_arg args;

    cli_args a[] = {
	{"-io",	  ARG_IO,  1, NULL, offsetof(io_arg, io)},
	{NULL,	  0,	   0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("order contigs");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (-1 == find_contig_order(interp, args.io))
	verror(ERR_WARN, "Order Contigs", "Failure in Order Contigs");

    return TCL_OK;
}

int
tcl_quality_clip(ClientData clientData,
		 Tcl_Interp *interp,
		 int argc,
		 char *argv[])
{
    qclip_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(qclip_arg, io)},
	{"-contigs", ARG_STR, 1, NULL, offsetof(qclip_arg, contig)},
	{"-quality", ARG_INT, 1, "20", offsetof(qclip_arg, quality)},
	{NULL,	  0,	   0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("quality clip");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contig, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }

    quality_clip(args.io, num_contigs, contig_array, args.quality);

    xfree(contig_array);

    return TCL_OK;
}

int
tcl_N_clip(ClientData clientData,
	   Tcl_Interp *interp,
	   int argc,
	   char *argv[])
{
    list2_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(list2_arg, io)},
	{"-contigs", ARG_STR, 1, NULL, offsetof(list2_arg, inlist)},
	{NULL,	  0,	   0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("N-base clip");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }

    N_clip(args.io, num_contigs, contig_array);

    xfree(contig_array);

    return TCL_OK;
}

int
tcl_quality_clip_ends(ClientData clientData,
		      Tcl_Interp *interp,
		      int argc,
		      char *argv[])
{
    qclip_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs;
    int i;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(qclip_arg, io)},
	{"-contigs", ARG_STR, 1, NULL, offsetof(qclip_arg, contig)},
	{"-quality", ARG_INT, 1, "15", offsetof(qclip_arg, quality)},
	{NULL,	  0,	   0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("Quality clip ends");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contig, &num_contigs, &contig_array);
    for (i = 0; i < num_contigs; i++) {
	quality_clip_ends(args.io, contig_array[i].contig, args.quality);
    }

    xfree(contig_array);

    return TCL_OK;
}

int
tcl_difference_clip(ClientData clientData,
		    Tcl_Interp *interp,
		    int argc,
		    char *argv[])
{
    dclip_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(dclip_arg, io)},
	{"-contigs", ARG_STR, 1, NULL, offsetof(dclip_arg, inlist)},
	{"-tag",     ARG_INT, 1, "0",  offsetof(dclip_arg, tag)},
	{NULL,	  0,	   0, NULL, 0}
    };

    if (get_licence_type() == LICENCE_VIEWER) return TCL_ERROR;

    vfuncheader("difference clip");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.inlist, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }

    difference_clip(args.io, num_contigs, contig_array, args.tag);

    xfree(contig_array);

    return TCL_OK;
}

int
NumReadingsInContig(ClientData clientData,
		    Tcl_Interp *interp,
		    int argc,
		    char *argv[])
{
    int i, c_num;
    int cnt = 0;
    GapIO *io;
    int handle;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io contig_number\"\n", *argv);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    c_num = atoi(argv[2]);
    for (i = io_clnbr(io, c_num); i; i = io_rnbr(io, i)) {
	cnt++;
    }
    vTcl_SetResult(interp, "%d", cnt);

    return TCL_OK;
}

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



int tcl_find_tags(ClientData clientData, Tcl_Interp *interp,
		  int objc, Tcl_Obj *CONST objv[])
{
    find_tags_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    Array al;

    cli_args a[] = {
	{"-io",		ARG_IO,	 1, NULL, offsetof(find_tags_arg, io)},
	{"-contigs",	ARG_STR, 1, NULL, offsetof(find_tags_arg, contigs)},
	{"-tag_types",	ARG_STR, 1, "",   offsetof(find_tags_arg, tag_list)},
	{NULL,		0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }

    if (SetActiveTags(args.tag_list) == -1) {
	return TCL_ERROR;
    }

    al = find_tags(args.io, contig_array, num_contigs,
		   active_tag_types, number_of_active_tags);

    xfree(contig_array);

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

int
tcl_consistency_display(ClientData clientData,
			Tcl_Interp *interp,
			int argc,
			char *argv[])
{
    consistency_arg args;
    contig_list_t *contig_array = NULL;
    int num_contigs = 0;
    int *contigs;
    ruler_s *ruler;
    cursor_s cursor;
    int start, end;
    int i;

    cli_args a[] = {
	{"-io",	         ARG_IO,  1, NULL, offsetof(consistency_arg, io)},
	{"-contigs",     ARG_STR, 1, NULL, offsetof(consistency_arg, contigs)},
	{"-frame",       ARG_STR, 1, NULL, offsetof(consistency_arg, frame)},
	{"-win_ruler", ARG_STR, 1, NULL,
	     offsetof(consistency_arg, r_win)},
	{"-cursor_width", ARG_INT,1, "-1",
	     offsetof(consistency_arg, cursor_wd)},
	{"-cursor_fill",  ARG_STR, 1,  "",
	     offsetof(consistency_arg, cursor_fill)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    active_list_contigs(args.io, args.contigs, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	if (contig_array)
	    xfree(contig_array);
	return TCL_OK;
    }
    contigs = to_contigs_only(num_contigs, contig_array);

    /* set start and end where end is the total length of all contigs */
    start = contig_array[0].start;
    end = 0;
    for (i = 0; i < num_contigs; i++) {
        end += contig_array[i].end;
    }
    xfree(contig_array);

    cursor = cursor_struct(interp, gap_defs, "CONSISTENCY_DISPLAY", args.cursor_wd,
			   args.cursor_fill);

    ruler = ruler_struct(interp, gap_defs, "CONSISTENCY_DISPLAY", 1);

    ruler->start = start;
    ruler->end = end;
    sprintf(ruler->window, "%s", args.r_win);

    vTcl_SetResult(interp, "%d",
		   consistency_reg(args.io, interp, contigs, num_contigs,
				   start, end, args.frame, ruler, cursor));
    return TCL_OK;
}

int
tcl_confidence_graph(ClientData clientData,
		     Tcl_Interp *interp,
		     int argc,
		     char *argv[])
{
    confidence_arg args;
    ruler_s *ruler;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(confidence_arg, io)},
	{"-id",     ARG_INT, 1, NULL, offsetof(confidence_arg, id)},
	{"-frame",  ARG_STR, 1, NULL, offsetof(confidence_arg, frame)},
	{"-window", ARG_STR, 1, NULL, offsetof(confidence_arg, conf_win)},
	{"-win_ruler", ARG_STR, 1, NULL, offsetof(confidence_arg, r_win)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("confidence graph");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ruler = ruler_struct(interp, gap_defs, "CONFIDENCE_GRAPH", 1);
    sprintf(ruler->window, "%s", args.r_win);

    vTcl_SetResult(interp, "%d",
		   confidence_graph_reg(args.io, interp, args.frame,
					args.conf_win, args.id, ruler,
					CONFIDENCE_GRAPH_QUALITY));
    return TCL_OK;
}

int
tcl_discrepancy_graph(ClientData clientData,
		      Tcl_Interp *interp,
		      int argc,
		      char *argv[])
{
    confidence_arg args;
    ruler_s *ruler;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(confidence_arg, io)},
	{"-id",     ARG_INT, 1, NULL, offsetof(confidence_arg, id)},
	{"-frame",  ARG_STR, 1, NULL, offsetof(confidence_arg, frame)},
	{"-window", ARG_STR, 1, NULL, offsetof(confidence_arg, conf_win)},
	{"-win_ruler", ARG_STR, 1, NULL, offsetof(confidence_arg, r_win)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("discrepancy graph");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ruler = ruler_struct(interp, gap_defs, "CONFIDENCE_GRAPH", 1);
    sprintf(ruler->window, "%s", args.r_win);

    vTcl_SetResult(interp, "%d",
		   confidence_graph_reg(args.io, interp, args.frame,
					args.conf_win, args.id, ruler,
					CONFIDENCE_GRAPH_DISCREP));
    return TCL_OK;
}

int
tcl_reading_coverage(ClientData clientData,
		     Tcl_Interp *interp,
		     int argc,
		     char *argv[])
{
    reading_cov_arg args;
    ruler_s *ruler;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(reading_cov_arg, io)},
	{"-id",     ARG_INT, 1, NULL, offsetof(reading_cov_arg, id)},
	{"-frame",  ARG_STR, 1, NULL, offsetof(reading_cov_arg, frame)},
	{"-window", ARG_STR, 1, NULL, offsetof(reading_cov_arg, conf_win)},
	{"-win_ruler", ARG_STR, 1, NULL, offsetof(reading_cov_arg, r_win)},
	{"-strand",   ARG_INT, 1, "3", offsetof(reading_cov_arg, strand)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("reading coverage");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ruler = ruler_struct(interp, gap_defs, "READING_COVERAGE", 1);
    sprintf(ruler->window, "%s", args.r_win);

    vTcl_SetResult(interp, "%d",
		   reading_coverage_reg(args.io, interp, args.frame,
					args.conf_win, args.id, ruler,
					args.strand));
    return TCL_OK;
}
int
tcl_readpair_coverage(ClientData clientData,
		      Tcl_Interp *interp,
		      int argc,
		      char *argv[])
{
    confidence_arg args;
    ruler_s *ruler;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(confidence_arg, io)},
	{"-id",     ARG_INT, 1, NULL, offsetof(confidence_arg, id)},
	{"-frame",  ARG_STR, 1, NULL, offsetof(confidence_arg, frame)},
	{"-window", ARG_STR, 1, NULL, offsetof(confidence_arg, conf_win)},
	{"-win_ruler", ARG_STR, 1, NULL, offsetof(confidence_arg, r_win)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("readpair coverage");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ruler = ruler_struct(interp, gap_defs, "READPAIR_COVERAGE", 1);
    sprintf(ruler->window, "%s", args.r_win);

    vTcl_SetResult(interp, "%d",
		   readpair_coverage_reg(args.io, interp, args.frame,
					args.conf_win, args.id, ruler));
    return TCL_OK;
}

int
tcl_strand_coverage(ClientData clientData,
		    Tcl_Interp *interp,
		    int argc,
		    char *argv[])
{
    strand_arg args;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(strand_arg, io)},
	{"-id",     ARG_INT, 1, NULL, offsetof(strand_arg, id)},
	{"-frame",  ARG_STR, 1, NULL, offsetof(strand_arg, frame)},
	{"-window", ARG_STR, 1, NULL, offsetof(strand_arg, win)},
	{"-strand",   ARG_INT, 1, "3", offsetof(strand_arg, strand)},
	{"-problems", ARG_INT, 1, "1", offsetof(strand_arg, problems)},
	{NULL,	    0,	     0, NULL, 0}
    };

    vfuncheader("strand coverage");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    vTcl_SetResult(interp, "%d",
		   strand_coverage_reg(args.io, interp, args.frame,
				       args.win, args.id, args.strand,
				       args.problems));
    return TCL_OK;
}

int tcl_consistency_contig(ClientData clientData,
			   Tcl_Interp *interp,
			   int argc,
			   char *argv[])
{
    obj_consistency_disp *c;
    t_cursor_arg args;
    double wx, wy;

    cli_args a[] = {
	{"-io",	ARG_IO,  1, NULL, offsetof(t_cursor_arg, io)},
	{"-id", ARG_INT, 1, NULL, offsetof(t_cursor_arg, id)},
	{"-x",  ARG_INT, 1, NULL, offsetof(t_cursor_arg, cx)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    c = result_data(args.io, args.id, 0);

    /* fill in global position of cursor in label box */
    CanvasToWorld(c->win_list[0]->canvas, args.cx, 0, &wx, &wy);
    vTcl_SetResult(interp, "%d",
		   find_cursor_contig(args.io, args.id, c->contig_offset,
				      c->contigs, c->num_contigs, wx));

    return TCL_OK;
}
int tcl_create_consistency_ruler(ClientData clientData,
				 Tcl_Interp *interp,
				 int argc,
				 char *argv[])
{
    obj_consistency_disp *c;
    cons_world_arg args;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(cons_world_arg, io)},
	{"-id",     ARG_INT, 1, NULL, offsetof(cons_world_arg, id)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    c = result_data(args.io, args.id, 0);

    create_consistency_ruler(args.io, c);
    return TCL_OK;
}

int tcl_delete_consistency_ruler(ClientData clientData,
				 Tcl_Interp *interp,
				 int argc,
				 char *argv[])
{
    obj_consistency_disp *c;
    del_cons_ruler_arg args;
    int win_num;

    cli_args a[] = {
	{"-io",	    ARG_IO,  1, NULL, offsetof(del_cons_ruler_arg, io)},
	{"-id",     ARG_INT, 1, NULL, offsetof(del_cons_ruler_arg, id)},
	{"-window", ARG_STR, 1, NULL, offsetof(del_cons_ruler_arg, window)},
	{NULL,	0,	0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    c = result_data(args.io, args.id, 0);

    win_num = get_consistency_win_num(c, args.id);
    delete_consistency_window(c, win_num);
    deleteWindow(c->win_list, &c->num_wins, args.window);

    if (c->num_wins == 0)
        consistency_shutdown(args.io, c);

    return TCL_OK;
}

int tcl_rightmost_read(ClientData clientData,
		       Tcl_Interp *interp,
		       int argc,
		       char *argv[])
{
    int err;
    GContigs c;
    contig_arg args;
    char *rightmost_read;
    int read, right_read;
    int not_at_end = 1;
    GReadings r;

    cli_args a[] = {
	{"-io",	     ARG_IO,  1, NULL, offsetof(contig_arg, io)},
	{"-contig", ARG_INT, 1, NULL, offsetof(contig_arg, contig)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    err = GT_Read(args.io, arr(GCardinal, args.io->contigs, args.contig-1),
		  &c, sizeof(c), GT_Contigs);
    
    read = c.left;
    /*
    while (read = rr_read(args.io, read, c.length)) {
	right_read = read;
    }
    */
    
    while (not_at_end) {
	gel_read(args.io, read, r);
	if (r.position + r.sequence_length < c.length) {
	    if (r.right != 0) 
		read = r.right;
	} else {
	    not_at_end = 0;
	}
	right_read = read;
    }

    rightmost_read = io_rname(args.io, right_read);

    vTcl_SetResult(interp, "%s", rightmost_read);
    return TCL_OK;
}

int tcl_read_enz_file(ClientData clientData,
		      Tcl_Interp *interp,
		      int argc,
		      char *argv[])
{
    read_enz_arg args;
    int num_enzymes;
    int i;
    char **names;

    cli_args a[] = {
	{"-file",   ARG_STR, 1, NULL, offsetof(read_enz_arg, filename)},
	{NULL,	    0,	     0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /* open_r_enz_file(args.filename, &r_enzyme, &num_enzymes); */
    if (0 == r_enz_file_names(args.filename, &names, &num_enzymes))
	return TCL_OK;

    for (i =0; i < num_enzymes; i++) {
	Tcl_AppendElement(interp, names[i]);
	xfree(names[i]);
    }

    if (num_enzymes)
	xfree(names);

    return TCL_OK;
}


int
tcl_save_contig_order(ClientData clientData,
		      Tcl_Interp *interp,
		      int argc,
		      char *argv[])
{
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

    if (-1 == gap_parse_args(a, &args, argc, argv))
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

    return TCL_OK;
}

#if 0
int tcl_shuffle_pads(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[])
{
    list2_arg args;
    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL,  offsetof(list2_arg, io)},
	{"-contigs",	ARG_STR, 1, "*",  offsetof(list2_arg, inlist)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    shuffle_contigs_io(args.io);

    return TCL_OK;
}
#endif


/* set up tcl commands which call C procedures */
/*****************************************************************************/
/*				   NewGap_Init				     */
/*****************************************************************************/
int
NewGap_Init(Tcl_Interp *interp) {
    init_globals(interp);

    Tcl_CreateCommand(interp, "disassemble_readings", DisReadings,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "break_contig", BreakContig,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "open_db", OpenDB,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "close_db", CloseDB,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "copy_db", CopyDB,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_tag_array", tcl_get_tag_array,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "db_info", db_info,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "display_ruler", DisplayRuler,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "update_template_display", UpdateTemplateDisplay,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "add_template_plot", AddTemplatePlot,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "template_win_free", TemplateWinFree,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "print_template_readings",
		      PrintTemplateReadings,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "display_ruler_ticks", DisplayRulerTicks,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "display_reading_tags", DisplayReadingTags,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "show_relationships", ShowRelationships,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "find_internal_joins", FindInternalJoins,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "find_read_pairs", FindReadPairs,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "plot_quality", PlotQuality,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "find_repeats", FindRepeats,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "update_contig_order", UpdateContigOrder,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "refresh_contig_order", RefreshContigOrder,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "flush_contig_order", FlushContigOrder,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "result_names", tk_result_names,
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
    Tcl_CreateCommand(interp, "complement_contig", tk_complement_contig,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "register_id", tk_register_id,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "result_time", tk_result_time,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "reg_get_ops", tk_reg_get_ops,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "reg_invoke_op", tk_reg_invoke_op,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "reg_notify_update", tk_reg_notify_update,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "reg_notify_highlight", tk_reg_notify_highlight,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "result_is_2d", tk_result_is_2d,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "result_is_consistency", tk_result_is_consistency,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "result_delete", tk_result_delete,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "result_quit", tk_result_quit,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "clear_cp", tk_clear_cp,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "clear_template", tk_clear_template,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "clear_consistency", tk_clear_consistency,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "matchresult_configure",
		      tk_matchresult_configure,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "contig_register", tk_contig_register,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "contig_deregister", tk_contig_deregister,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "contig_notify", tk_contig_notify,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "cursor_ref", tk_cursor_ref,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "query_cursor", tk_query_cursor,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "create_cursor", tk_create_cursor,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "delete_cursor", tk_delete_cursor,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "find_long_gels", FindLongGels,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "find_taq_terminator", FindTaqTerminator,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "check_database", CheckDatabase,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "remove_contig_holes", tcl_remove_contig_holes,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "delete_contig", DeleteContig,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "annotation_address", AnnotationAddress,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "reset_contig_order", ResetContigOrder,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "quit_displays", tcl_quit_displays,
		      (ClientData) NULL,
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "double_strand", DoubleStrand,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "find_primers", FindPrimers,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "extract_readings", ExtractReadings,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "pre_assemble", PreAssemble,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "get_consensus", Consensus,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "calc_consensus", tcl_calc_consensus,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateObjCommand(interp, "calc_quality", tcl_calc_quality,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateObjCommand(interp, "calc_consensus_double",
			 tcl_calc_consensus_double,
			 (ClientData) NULL,
			 NULL);
    Tcl_CreateCommand(interp, "display_quality", DisplayQuality,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "display_template_quality", DisplayTemplateQuality,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "display_templates", DisplayTemplates,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "minimal_coverage", MinimalCoverage,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "unattached_readings", UnattachedReadings,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "enter_tags", EnterTags,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "display_contig_selector", DisplayContigSelector,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "display_contig_comparator",
		      DisplayContigComparator,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "display_cs_tags", DisplayCSTags,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "display_cs_diagonal", DisplayCSDiagonal,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "plot_template_renz", PlotTemplateREnz,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "plot_renz", PlotREnz,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "get_r_enz_name", GetREnzName,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "get_r_enz_info", GetREnzInfo,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "create_renz_tags", CreateREnzTags,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "plot_stop_codons", PlotStopCodons,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "get_codon_name", GetCodonName,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "refresh_codon_map", RefreshCodonMap,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "assemble_direct", tcl_assemble_direct,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "semaphore_free", check_active_lock,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "check_assembly", tcl_check_assembly,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "anno_list", tcl_anno_list,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "delete_anno_list", tcl_delete_anno_list,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "find_probes", tcl_find_probes,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "find_oligo", FindOligo,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "add_tags", tcl_add_tags,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "contig_order_to_number",
		      tcl_contig_order_to_number,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "get_read_names", tcl_get_read_names,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "scroll_canvas", ScrollCanvas,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "zoom_canvas", ZoomCanvas,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "resize_canvas", ResizeCanvas,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "delete_canvas_cursor", DeleteCanvasCursor,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "draw_canvas_cursor_x", DrawCanvasCursor_X,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "draw_canvas_cursor_y", DrawCanvasCursor_Y,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "delete_window", DeleteWindow,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "template_contig", TemplateContig,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "update_template_contig_order",
		      UpdateTemplateContigOrder,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "remove_contig_duplicates",
		      RemoveContigDuplicates,
		      (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "canvas_to_world",
		      Canvas_To_World,
		      (ClientData)NULL, NULL);

    Tcl_CreateCommand(interp, "assemble_shotgun",
		      MainAssembly, (ClientData)1, NULL);
    Tcl_CreateCommand(interp, "assemble_screen",
		      MainAssembly, (ClientData)2, NULL);
    Tcl_CreateCommand(interp, "assemble_one_contig",
		      MainAssembly, (ClientData)3, NULL);
    Tcl_CreateCommand(interp, "assemble_new_contigs",
		      MainAssembly, (ClientData)4, NULL);
    Tcl_CreateCommand(interp, "assemble_single_strand",
		      MainAssembly, (ClientData)5, NULL);
    Tcl_CreateCommand(interp, "assemble_independent",
		      MainAssembly, (ClientData)6, NULL);
    Tcl_CreateCommand(interp, "order_contigs",
		      tcl_order_contigs, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "quality_clip",
		      tcl_quality_clip, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "N_clip",
		      tcl_N_clip, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "quality_clip_ends",
		      tcl_quality_clip_ends, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "difference_clip",
		      tcl_difference_clip, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "num_readings_in_contig",
		      NumReadingsInContig, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "list_confidence",
			 tcl_list_confidence, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "list_base_confidence",
			 tcl_list_base_confidence, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "new_note",
			 tcl_new_note, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "delete_note",
			 tcl_delete_note, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "edit_note",
			 tcl_edit_note, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "find_tags",
			 tcl_find_tags, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "load_alignment_matrix",
		      tcl_load_alignment_matrix, (ClientData)NULL, NULL);
    Tcl_CreateObjCommand(interp, "load_genetic_code",
			 tcl_load_genetic_code, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "confidence_graph",
		      tcl_confidence_graph, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "discrepancy_graph",
		      tcl_discrepancy_graph, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "reading_coverage",
		      tcl_reading_coverage, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "readpair_coverage",
		      tcl_readpair_coverage, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "strand_coverage",
		      tcl_strand_coverage, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "consistency_display",
		      tcl_consistency_display, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "consistency_contig",
		      tcl_consistency_contig, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "create_consistency_ruler",
		      tcl_create_consistency_ruler, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "delete_consistency_ruler",
		      tcl_delete_consistency_ruler, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "rightmost_read",
		      tcl_rightmost_read, (ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "read_enz_file", tcl_read_enz_file,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "save_contig_order", tcl_save_contig_order,
		      (ClientData) NULL,
		      NULL);
    
#if 0
    Tcl_CreateObjCommand(interp, "shuffle_pads", tcl_shuffle_pads,
			 (ClientData) NULL, NULL);
#endif

    return TCL_OK;

} /* end NewGap_Init */
