#include <stdio.h>
#include <tk.h>
#include <string.h>

#include "IO.h"
#include "tk-io-reg.h"
#include "io-reg.h"
#include "gap_cli_arg.h"
#include "newgap_cmds.h" /* GetInterp() */
#include "contig_selector.h" /* csplot_hash */
#include "newgap_structs.h" /* contig_arg */
#include "misc.h"
#include "tcl_utils.h"
#include "template_display.h"
#include "tclXkeylist.h"
#include "consistency_display.h"

static void reg_init_args(char *args);
static char *reg_get_arg(char *name);

typedef struct {
    GapIO *io;
} rn_arg;

/*
 * Create a list of results and associated contig numbers. The format of the
 * list is "{contig regnum id string} ?{contig regnum id string}? ..."
 *
 * The list is stored in the variable specified on the command line.
 */
int tk_result_names(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {
    char *line;
    int contig, reg, id;
    rn_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rn_arg, io)},
        {NULL,        0,       0, NULL, 0}
    };
    char buf[1024];
    Tcl_DString ds;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    contig = -1;
    line = result_names(args.io, &contig, &reg, &id, 1);

    Tcl_DStringInit(&ds);
    while (line) {
	if (*line) {
	    sprintf(buf, "%d %d %d {%s}", contig, reg, id, line);
	    Tcl_DStringAppendElement(&ds, buf);
	}
	contig = -1;
	line = result_names(args.io, &contig, &reg, &id, 0);
    }

    Tcl_DStringResult(interp, &ds);

    return TCL_OK;
}


/*
 * Tk interface to register_id() - all this work for such a miniscule 2 line
 * function!!!
 */
int tk_register_id(ClientData clientData, Tcl_Interp *interp,
		   int argc, char **argv) {
    vTcl_SetResult(interp, "%d", register_id());

    return TCL_OK;
}


typedef struct {
    GapIO *io;
    int contig;
    int id;
} rt_arg;

/*
 * Tk interface to result_time()
 */
int tk_result_time(ClientData clientData, Tcl_Interp *interp,
		   int argc, char **argv) {
    rt_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rt_arg, io)},
	{"-contig",   ARG_INT, 1, NULL, offsetof(rt_arg, contig)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(rt_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    Tcl_SetResult(interp, result_time(args.io, args.contig, args.id),
		  TCL_VOLATILE);

    return TCL_OK;
}


typedef struct {
    GapIO *io;
    int id;
} go_arg;

/*
 * Fetches a list of operations for a particular result (id)
 */
int tk_reg_get_ops(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {
    go_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(go_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(go_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_get_ops ro;
    int l;
    char *ops;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    ro.job = REG_GET_OPS;
    ro.ops = NULL;
    result_notify(args.io, args.id, (reg_data *)&ro, 0);

    if (!ro.ops)
	return TCL_OK; /* blank list */

    Tcl_ResetResult(interp);
    ops = ro.ops;
    while(l = strlen(ops)) {
	Tcl_AppendElement(interp, ops);
	ops += l+1;
    }

    return TCL_OK;
}


typedef struct {
    GapIO *io;
    int id;
    int option;
} iop_arg;

/*
 * Performs an operation for a particular result (id)
 */
int tk_reg_invoke_op(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv) {
    iop_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(iop_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(iop_arg, id)},
	{"-option",   ARG_INT, 1, NULL, offsetof(iop_arg, option)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_invoke_op inv;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    inv.job = REG_INVOKE_OP;
    inv.op = args.option;

    result_notify(args.io, args.id, (reg_data *)&inv, 0);

    return TCL_OK;
}

/*
 * Notifies a contig of a general update. Contig '0' implies all contigs
 */
int tk_reg_notify_update(ClientData clientData, Tcl_Interp *interp,
			 int argc, char **argv) {
    contig_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(contig_arg, io)},
        {"-contig",   ARG_INT, 1, NULL, offsetof(contig_arg, contig)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_length rl;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    rl.job = REG_LENGTH;

    if (args.contig) {
	rl.length = io_clength(args.io, args.contig);
	contig_notify(args.io, args.contig, (reg_data *)&rl);
    } else {
	int i;

	for (i = 0; i <= NumContigs(args.io); i++) {
	    if (i)
		rl.length = io_clength(args.io, i);
	    else
		rl.length = 0;
	    contig_notify(args.io, i, (reg_data *)&rl);
	}
    }

    return TCL_OK;
}


/*
 * Updates the result manager window. For efficiency reasons we only bother
 * updating when things become idle again.
 */
static int result_pending = 0;

void update_results_(ClientData clientData) {
    char buf[100];
    GapIO *io = (GapIO *)clientData;

    sprintf(buf, "result_list_update %d", *handle_io(io));

    Tcl_Eval(GetInterp(), buf);
    result_pending = 0;
}

void update_results(GapIO *io) {
    if (!result_pending) {
	result_pending = 1;

	Tcl_DoWhenIdle(update_results_, (ClientData)io);
    }
}

void busy_dialog(GapIO *io, int contig) {
    char buf[1024];

    sprintf(buf, "tk_messageBox \
			-icon warning \
			-title {Contig is busy} \
			-message {The contig %s is busy, probably due to an "
	    "editor in use for this contig. Changes will not be made for "
	    "this contig.} \
			-type ok",
	    get_contig_name(io, contig));

    Tcl_Eval(GetInterp(), buf);
}

typedef struct {
    char *result;
    char *colour;
    char *csplot;
    int width;
} conf_arg;

/*
 * Configures a result.
 */
int tk_matchresult_configure(ClientData clientData, Tcl_Interp *interp,
			     int argc, char **argv) {
    conf_arg args;
    cli_args a[] = {
	{"-result",  ARG_STR, 1, NULL, offsetof(conf_arg, result)},
	{"-colour",   ARG_STR, 1, "",   offsetof(conf_arg, colour)},
	{"-width",    ARG_INT, 1, "-1", offsetof(conf_arg, width)},
	{"-csplot",   ARG_STR, 1, NULL, offsetof(conf_arg, csplot)},
        {NULL,        0,       0, NULL, 0}
    };
    mobj_repeat *r;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    /*
     * Find the result associated with this tag. 'result' is a C pointer
     * held in an implementation defined manner as a tcl string.
     */
    r = (mobj_repeat *)TclPtr2C(args.result);

    /* Adjust configurations */
    if (*args.colour) {
	strncpy(r->colour, args.colour, COLOUR_LEN-1);
    }

    if (args.width != -1) {
	r->linewidth = args.width;
    }

    return TCL_OK;
}


/*
 * Attempts to shut down all active displays.
 * As used by alter relationships, assembly, etc.
 *
 * Args: quit_display handle function_name
 */
int tcl_quit_displays(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    int ret = 1, handle, i;
    GapIO *io;
    reg_quit rq;
    int cnum;

    if (argc != 3) {
	Tcl_SetResult(interp, "wrong # args:\n", TCL_STATIC);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    io = io_handle(&handle);

    rq.job = REG_QUIT;
    rq.lock = REG_LOCK_WRITE;
    cnum = 0;
    for (i=0; i<=NumContigs(io); i++) {
	contig_notify(io, i, (reg_data *)&rq);
	if (!(rq.lock & REG_LOCK_WRITE)) {
	    ret = 0;
	    cnum = i;
	}
    }

    if (!ret) {
	verror(ERR_WARN, argv[2], "Database busy");
	verror(ERR_WARN, argv[2], "Please shut down editing displays");
	if (cnum)
	    busy_dialog(io, cnum);
    }

    vTcl_SetResult(interp, "%d", ret);

    return TCL_OK;
}

typedef struct {
    GapIO *io;
    int id;
} rd_arg;

/*
 * Send a delete request to a specific result
 */
int tk_result_delete(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv) {
    rd_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rd_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(rd_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_delete rd;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    rd.job = REG_DELETE;
    result_notify(args.io, args.id, (reg_data *)&rd, 0);

    return TCL_OK;
}

/*
 * Send a delete request to a specific result
 */
int tk_result_quit(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv) {
    rd_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rd_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(rd_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_delete rd;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    rd.job = REG_QUIT;
    result_notify(args.io, args.id, (reg_data *)&rd, 0);

    return TCL_OK;
}

/*
 * Delete all contig comparator displays.
 */
int tk_clear_cp(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv) {
    obj_cs *cs;
    rd_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rd_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(rd_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_quit rq;


    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    rq.job = REG_QUIT;
    rq.lock = REG_LOCK_WRITE;

    type_notify(args.io, REG_TYPE_FIJ,      (reg_data *)&rq, 1);
    type_notify(args.io, REG_TYPE_READPAIR, (reg_data *)&rq, 1);
    type_notify(args.io, REG_TYPE_REPEAT,   (reg_data *)&rq, 1);
    type_notify(args.io, REG_TYPE_CHECKASS, (reg_data *)&rq, 1);
    type_notify(args.io, REG_TYPE_OLIGO,    (reg_data *)&rq, 1);

    cs = result_data(args.io, args.id, 0);
    strcpy(cs->window, cs->hori);
    cs->vert[0] = '\0';
    return TCL_OK;
}

/*
 * Delete all template displays.
 */
int tk_clear_template(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    template_arg args;
    int i;
    obj_template_disp *t;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(template_arg, io)},
        {"-id",       ARG_INT, 1, NULL, offsetof(template_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_quit rq;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    rq.job = REG_QUIT;
    rq.lock = REG_LOCK_WRITE;

    t = result_data(args.io, args.id, 0);
    for (i = 0; i < t->num_wins; i++) {
	if (t->win_list[i]->id != t->id) {
	    int num_wins = t->num_wins;
	    result_notify(args.io, t->win_list[i]->id, (reg_data *)&rq, 1);
	    i -= num_wins - t->num_wins;
	}
    }
/*
    type_notify(args.io, REG_TYPE_QUALITY,     (reg_data *)&rq, 1);
    type_notify(args.io, REG_TYPE_RESTRICTION, (reg_data *)&rq, 1);
*/
    return TCL_OK;
}

/*
 * Delete all consistency displays.
 */
int tk_clear_consistency(ClientData clientData, Tcl_Interp *interp,
			 int argc, char **argv) {
    template_arg args;
    obj_consistency_disp *c;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(template_arg, io)},
        {"-id",       ARG_INT, 1, NULL, offsetof(template_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    c = result_data(args.io, args.id, 0);
    clear_consistency(args.io, c);
    return TCL_OK;
}

typedef struct {
    GapIO *io;
    char *reading;
    int value;
} notify_arg;

/*
 * Sends a highlight reading notification
 */
int tk_reg_notify_highlight(ClientData clientData, Tcl_Interp *interp,
			    int argc, char **argv) {
    notify_arg args;
    int r_num;
    int is_name;

    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL, offsetof(notify_arg, io)},
	{"-reading",	ARG_STR, 1, NULL, offsetof(notify_arg, reading)},
	{"-highlight",	ARG_INT, 1, NULL, offsetof(notify_arg, value)},
        {NULL,          0,       0, NULL, 0}
    };
    reg_highlight_read hr;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    if (args.reading[0] == '=' || args.reading[0] == '#') {
	is_name = GGN_ID;
    } else {
	is_name = GGN_NAME;
    }
    r_num = get_gel_num(args.io, args.reading, is_name);
    if (r_num <= 0) {
	verror(ERR_WARN, "reg_notify_hightlight", "Unknown reading '%s'",
	       args.reading);
	return TCL_OK;
    }

    hr.job = REG_HIGHLIGHT_READ;
    hr.seq = r_num;
    hr.val = args.value;

    contig_notify(args.io,
		  rnumtocnum(args.io, chain_left(args.io, hr.seq)),
		  (reg_data *)&hr);

    return TCL_OK;
}

/*
 * Determines whether a result is a component of the 2D contig comparator
 * display. Returns 1 or 0.
 */
int tk_result_is_2d(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {
    rd_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rd_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(rd_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    contig_reg_t **regs;
    int r=0;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    regs = result_to_regs(args.io, args.id);
    if (regs) {
	if ((regs[0]->type == REG_TYPE_FIJ) ||
	    (regs[0]->type == REG_TYPE_READPAIR) ||
	    (regs[0]->type == REG_TYPE_REPEAT) ||
	    (regs[0]->type == REG_TYPE_CHECKASS) ||
	    (regs[0]->type == REG_TYPE_OLIGO)) {
	    r=1;
	}
    }

    xfree(regs);

    vTcl_SetResult(interp, "%d", r);

    return TCL_OK;
}

typedef struct {
    GapIO *io;
    int id;
    int cons_id;
} cons_arg;

/*
 * Determines whether a result is a component of the consistency
 * display. Returns 1 or 0.
 */
int tk_result_is_consistency(ClientData clientData, Tcl_Interp *interp,
			     int argc, char **argv) {
    cons_arg args;
    reg_generic gen;

    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(cons_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(cons_arg, id)},
	{"-cons_id",  ARG_INT, 1, NULL, offsetof(cons_arg, cons_id)},
        {NULL,        0,       0, NULL, 0}
    };
    contig_reg_t **regs;
    int r=0;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    gen.job = REG_GENERIC;
    gen.task = TASK_CONS_ID;

    regs = result_to_regs(args.io, args.id);
    if (regs) {
      if ((regs[0]->type == REG_TYPE_CONFIDENCE) ||
	  (regs[0]->type == REG_TYPE_READING_COVERAGE) ||
	  (regs[0]->type == REG_TYPE_READPAIR_COVERAGE) ||
	  (regs[0]->type == REG_TYPE_STRAND_COVERAGE)) {

	result_notify(args.io, args.id, (reg_data *)&gen, 0);
	if (args.cons_id == (int)gen.data) {
	  r=1;
	}
      }
    }

    xfree(regs);
    vTcl_SetResult(interp, "%d", r);

    return TCL_OK;
}


/*
 *----------------------------------------------------------------------------
 * Register Tcl procedures with contig/notifications
 *----------------------------------------------------------------------------
 */
typedef struct {
    Tcl_Interp *interp;
    char *command;
    int id;
} cr_type;

static void tk_contig_register_cmd(GapIO *io, int contig, void *fdata,
reg_data *jdata);

static int reg_str2flags(Tcl_Interp *interp, char *str) {
    int num_flags;
    int flags;
    char **str_flags;
    int i;

    if (Tcl_SplitList(interp, str, &num_flags, &str_flags) != TCL_OK)
	return 0;

    flags = 0;
    for (i = 0; i < num_flags; i++) {
	if (strcmp(str_flags[i], "GENERIC") == 0) {
	    flags |= REG_GENERIC;
	} else if (strcmp(str_flags[i], "NUMBER_CHANGE") == 0) {
	    flags |= REG_NUMBER_CHANGE;
	} else if (strcmp(str_flags[i], "JOIN_TO") == 0) {
	    flags |= REG_JOIN_TO;
	} else if (strcmp(str_flags[i], "ORDER") == 0) {
	    flags |= REG_ORDER;
	} else if (strcmp(str_flags[i], "LENGTH") == 0) {
	    flags |= REG_LENGTH;
	} else if (strcmp(str_flags[i], "QUERY_NAME") == 0) {
	    flags |= REG_QUERY_NAME;
	} else if (strcmp(str_flags[i], "DELETE") == 0) {
	    flags |= REG_DELETE;
	} else if (strcmp(str_flags[i], "GET_LOCK") == 0) {
	    flags |= REG_GET_LOCK;
	} else if (strcmp(str_flags[i], "SET_LOCK") == 0) {
	    flags |= REG_SET_LOCK;
	} else if (strcmp(str_flags[i], "COMPLEMENT") == 0) {
	    flags |= REG_COMPLEMENT;
	} else if (strcmp(str_flags[i], "PARAMS") == 0) {
	    flags |= REG_PARAMS;
	} else if (strcmp(str_flags[i], "QUIT") == 0) {
	    flags |= REG_QUIT;
	} else if (strcmp(str_flags[i], "CURSOR_NOTIFY") == 0) {
	    flags |= REG_CURSOR_NOTIFY;
	} else if (strcmp(str_flags[i], "GET_OPS") == 0) {
	    flags |= REG_GET_OPS;
	} else if (strcmp(str_flags[i], "INVOKE_OP") == 0) {
	    flags |= REG_INVOKE_OP;
	} else if (strcmp(str_flags[i], "ANNO") == 0) {
	    flags |= REG_ANNO;
	} else if (strcmp(str_flags[i], "REGISTER") == 0) {
	    flags |= REG_REGISTER;
	} else if (strcmp(str_flags[i], "DEREGISTER") == 0) {
	    flags |= REG_DEREGISTER;
	} else if (strcmp(str_flags[i], "HIGHLIGHT_READ") == 0) {
	    flags |= REG_HIGHLIGHT_READ;
	} else if (strcmp(str_flags[i], "BUFFER_START") == 0) {
	    flags |= REG_BUFFER_START;
	} else if (strcmp(str_flags[i], "BUFFER_END") == 0) {
	    flags |= REG_BUFFER_END;
	} else if (strcmp(str_flags[i], "REQUIRED") == 0) {
	    flags |= REG_REQUIRED;
	} else if (strcmp(str_flags[i], "DATA_CHANGE") == 0) {
	    flags |= REG_DATA_CHANGE;
	} else if (strcmp(str_flags[i], "OPS") == 0) {
	    flags |= REG_OPS;
	} else if (strcmp(str_flags[i], "LOCKS") == 0) {
	    flags |= REG_LOCKS;
	} else if (strcmp(str_flags[i], "REGISTERS") == 0) {
	    flags |= REG_REGISTERS;
	} else if (strcmp(str_flags[i], "BUFFER") == 0) {
	    flags |= REG_BUFFER;
	} else if (strcmp(str_flags[i], "NOTE") == 0) {
	    flags |= REG_NOTE;
	} else if (strcmp(str_flags[i], "ALL") == 0) {
	    flags |= REG_ALL;
	}
    }

    Tcl_Free((char *)str_flags);

    return flags;
}

static char *reg_flags2str(int flags) {
    /* str is long enough to hold all registration flags */
    static char str[1024];

    *str = 0;
    if (flags & REG_GENERIC)
	strcat(str, "GENERIC ");
    if (flags & REG_NUMBER_CHANGE)
	strcat(str, "NUMBER_CHANGE ");
    if (flags & REG_JOIN_TO)
	strcat(str, "JOIN_TO ");
    if (flags & REG_ORDER)
	strcat(str, "ORDER ");
    if (flags & REG_LENGTH)
	strcat(str, "LENGTH ");
    if (flags & REG_QUERY_NAME)
	strcat(str, "QUERY_NAME ");
    if (flags & REG_DELETE)
	strcat(str, "DELETE ");
    if (flags & REG_GET_LOCK)
	strcat(str, "GET_LOCK ");
    if (flags & REG_SET_LOCK)
	strcat(str, "SET_LOCK ");
    if (flags & REG_COMPLEMENT)
	strcat(str, "COMPLEMENT ");
    if (flags & REG_PARAMS)
	strcat(str, "PARAMS ");
    if (flags & REG_QUIT)
	strcat(str, "QUIT ");
    if (flags & REG_CURSOR_NOTIFY)
	strcat(str, "CURSOR_NOTIFY ");
    if (flags & REG_GET_OPS)
	strcat(str, "GET_OPS ");
    if (flags & REG_INVOKE_OP)
	strcat(str, "INVOKE_OP ");
    if (flags & REG_ANNO)
	strcat(str, "ANNO ");
    if (flags & REG_REGISTER)
	strcat(str, "REGISTER ");
    if (flags & REG_DEREGISTER)
	strcat(str, "DEREGISTER ");
    if (flags & REG_HIGHLIGHT_READ)
	strcat(str, "HIGHLIGHT_READ ");
    if (flags & REG_BUFFER_START)
	strcat(str, "BUFFER_START ");
    if (flags & REG_BUFFER_END)
	strcat(str, "BUFFER_END ");
    if (flags & REG_NOTE)
	strcat(str, "NOTE ");

    return str;
}

static char *reg_type2str(int type) {
    switch(type) {
    case REG_TYPE_EDITOR:
	return "TYPE_EDITOR";

    case REG_TYPE_FIJ:
	return "TYPE_FIJ";

    case REG_TYPE_READPAIR:
	return "TYPE_READPAIR";

    case REG_TYPE_REPEAT:
	return "TYPE_REPEAT";

    case REG_TYPE_QUALITY:
	return "TYPE_QUALITY";

    case REG_TYPE_TEMPLATE:
	return "TYPE_TEMPLATE";

    case REG_TYPE_RESTRICTION:
	return "TYPE_RESTRICTION";

    case REG_TYPE_STOPCODON:
	return "TYPE_STOPCODON";

    case REG_TYPE_CONTIGSEL:
	return "TYPE_CONTIGSEL";

    case REG_TYPE_CHECKASS:
	return "TYPE_CHECKASS";

    case REG_TYPE_OLIGO:
	return "TYPE_OLIGO";

    default:
	return "TYPE_UNKNOWN";
    }
}

/*
 * Turns a cursor job string into an integer.
 */
static int cjob_str2int(Tcl_Interp *interp, char *str)
{
    int largc;
    char **largv;
    int i;
    int job = 0;

    if (str == NULL)
	return 0;

    if (Tcl_SplitList(interp, str, &largc, &largv) != TCL_OK)
	return 0;

    for (i = 0; i < largc; i++) {
	if (strcmp(largv[i], "INCREMENT") == 0)
	    job |= CURSOR_INCREMENT;
	if (strcmp(largv[i], "DECREMENT") == 0)
	    job |= CURSOR_DECREMENT;
	if (strcmp(largv[i], "MOVE") == 0)
	    job |= CURSOR_MOVE;
	if (strcmp(largv[i], "DELETE") == 0)
	    job |= CURSOR_DELETE;
    }

    Tcl_Free((char *)largv);
    return job;
}

typedef struct {
    GapIO *io;
    int contig;
    int id;
    char *command;
    char *flags;
    char *type;
} cr_arg;

int tk_contig_register(ClientData clientData, Tcl_Interp *interp,
		       int argc, char **argv) {
    cr_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(cr_arg, io)},
	{"-contig",   ARG_INT, 1, NULL, offsetof(cr_arg, contig)},
	{"-command",  ARG_STR, 1, NULL, offsetof(cr_arg, command)},
	{"-flags",    ARG_STR, 1, "",   offsetof(cr_arg, flags)},
	{"-type",     ARG_STR, 1, "TYPE_UNKNOWN", offsetof(cr_arg, type)},
        {NULL,        0,       0, NULL, 0}
    };
    int flags, type;
    cr_type *crt;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    if (NULL == (crt = (cr_type *)xmalloc(sizeof(*crt))))
	return TCL_ERROR;

    crt->interp = interp;
    crt->command = strdup(args.command);
    crt->id = register_id();
    flags = reg_str2flags(interp, args.flags);
    type = REG_TYPE_UNKNOWN; /* FIXME */

    contig_register(args.io, args.contig, tk_contig_register_cmd,
		    (void *)crt, crt->id, flags, type);

    vTcl_SetResult(interp, "%d", crt->id);

    return TCL_OK;
}

typedef struct {
    GapIO *io;
    int id;
} cdr_arg;

/*
 * Deregisters a tcl contig_register call. This requires 'id' returned
 * by the contig_register as tcl doesn't have an easy way of rememebering
 * the 'fdata' paramater.
 */
int tk_contig_deregister(ClientData clientData, Tcl_Interp *interp,
			 int argc, char **argv) {
    cdr_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(cdr_arg, io)},
	{"-id",       ARG_INT, 1, NULL, offsetof(cdr_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    contig_reg_t **regs, *ptr;
    int ret = 0, i, j;
    cr_type *crt;
    int n;
    int *uids;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    regs = result_to_regs(args.io, args.id);
    if (regs) {
	/* Find the number of registered items */
	for (n = 0, ptr = *regs; ptr != NULL; ptr = regs[++n])
	    ;

	/* Copy the unique ids */
	if (NULL == (uids = (int *)xmalloc(n * sizeof(int))))
	    return TCL_OK;
	for (i = 0; i < n; i++)
	    uids[i] = regs[i]->uid;

	/* Deregister each item, searching by uid */
	for (j = 0; j < n; j++) {
	    for (i = 0, ptr = *regs; ptr != NULL; ptr = regs[++i]) {
		if (ptr->uid == uids[j])
		    break;
	    }
	    if (ptr) {
		crt = (cr_type *)ptr->fdata;
		ret |= contig_deregister(args.io, 0, ptr->func, ptr->fdata);
		xfree(crt->command);
		xfree(crt);
	    }
	}

	xfree(uids);
    }
    vTcl_SetResult(interp, "%d", ret);
    xfree(regs);

    return TCL_OK;
}

static void tk_contig_register_cmd(GapIO *io, int contig, void *fdata,
				   reg_data *jdata)
{
    Tcl_DString ds;
    cr_type *crt = (cr_type *)fdata;
    char *type = reg_flags2str(jdata->job);
    char buf[1024];

    Tcl_DStringInit(&ds);
    switch (jdata->job) {
    case REG_NUMBER_CHANGE:
	sprintf(buf, "{number %d}",
		jdata->number.number);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_JOIN_TO:
	sprintf(buf, "{contig %d} {offset %d}",
		jdata->join.contig,
		jdata->join.offset);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_ORDER:
	sprintf(buf, "{pos %d}",
		jdata->order.pos);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_LENGTH:
	sprintf(buf, "{length %d}",
		jdata->length.length);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_QUERY_NAME:
	/* 'line' element is returned */
	break;

    case REG_DELETE:
    case REG_COMPLEMENT:
    case REG_ANNO:
    case REG_BUFFER_START:
    case REG_BUFFER_END:
	/* No extra data */
	break;

    case REG_GET_LOCK:
    case REG_SET_LOCK:
    case REG_QUIT:
	sprintf(buf, "{lock %d}",
		jdata->glock.lock);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_GET_OPS:
	/* 'ops' element is returned */
	break;

    case REG_INVOKE_OP:
	sprintf(buf, "{op %d}",
		jdata->invoke_op.op);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_PARAMS:
	/* 'params' element is returned */
	break;

    case REG_CURSOR_NOTIFY: {
	char job[1024];
	int first = 1;

	*job = 0;
	strcat(job, "{");
	if (jdata->cursor_notify.cursor->job & CURSOR_MOVE) {
	    strcat(job, "MOVE");
	    first = 0;
	}
	if (jdata->cursor_notify.cursor->job & CURSOR_INCREMENT) {
	    strcat(job, first ? "INCREMENT" : " INCREMENT");
	    first = 0;
	}
	if (jdata->cursor_notify.cursor->job & CURSOR_DECREMENT) {
	    strcat(job, first ? "DECREMENT" : " DECREMENT");
	    first = 0;
	}
	if (jdata->cursor_notify.cursor->job & CURSOR_DELETE) {
	    strcat(job, first ? "DELETE" : " DELETE");
	    first = 0;
	}
	strcat(job, "}");
	sprintf(buf, "{id %d} {seq %d} {pos %d} {abspos %d} {job %s}",
		jdata->cursor_notify.cursor->id,
		jdata->cursor_notify.cursor->seq,
		jdata->cursor_notify.cursor->pos,
		jdata->cursor_notify.cursor->abspos,
		job);
	Tcl_DStringAppend(&ds, buf, -1);
	break;
    }
    case REG_REGISTER:
    case REG_DEREGISTER:
	sprintf(buf, "{id %d} {type %s} {contig %d}",
		jdata->c_register.id,
		reg_type2str(jdata->c_register.type),
		jdata->c_register.contig);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_HIGHLIGHT_READ:
	sprintf(buf, "{seq %d} {val %d}",
		jdata->highlight.seq,
		jdata->highlight.val);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_NOTE:
	sprintf(buf, "{note %d} {task %s}",
		jdata->note.note,
		jdata->note.task == REG_NOTE_CREATE ? "CREATE" :
		jdata->note.task == REG_NOTE_DELETE ? "DELETE" : "EDIT");
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    default:
	Tcl_DStringAppend(&ds, "{unknown unknown}", -1);
	break;
    }

    sprintf(buf, "%d", crt->id);
    if (Tcl_VarEval(crt->interp, crt->command, " ", type, " ", buf, " ",
		    Tcl_DStringValue(&ds), NULL) != TCL_OK) {
	verror(ERR_WARN, "registration_callback", "%s",
	       Tcl_GetStringResult(crt->interp));
    }

    /* Now for the return value bit */
    switch(jdata->job) {
	static char buf[1025];
    case REG_QUERY_NAME:
	sprintf(buf, "%.80s", Tcl_GetStringResult(crt->interp));
	jdata->name.line = buf;
	break;

    case REG_GET_OPS:
	{
	    int largc, i, p;
	    char **largv;

	    if (Tcl_SplitList(crt->interp, Tcl_GetStringResult(crt->interp),
			      &largc, &largv) != TCL_OK)
		return;

	    for (p=i=0; i<largc; i++) {
		strcpy(&buf[p], largv[i]);
		p += strlen(largv[i])+1;
	    }
	    buf[p] = buf[p+1] = 0;
	    jdata->get_ops.ops = buf;

	    Tcl_Free((char *)largv);
	    break;
	}

    case REG_PARAMS:
	sprintf(buf, "%.1024s", Tcl_GetStringResult(crt->interp));
	jdata->params.string = buf;
	break;
    }

    Tcl_DStringFree(&ds);
}


typedef struct {
    GapIO *io;
    int cnum;
    char *type;
    char *args;
} cn_arg;

/*
 * An arbitray Tcl interface to contig event notification
 */
int tk_contig_notify(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv)
{
    cn_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(cn_arg, io)},
	{"-cnum",     ARG_INT, 1, NULL, offsetof(cn_arg, cnum)},
	{"-type",     ARG_STR, 1, NULL, offsetof(cn_arg, type)},
	{"-args",     ARG_STR, 1, NULL, offsetof(cn_arg, args)},
        {NULL,        0,       0, NULL, 0}
    };
    int flags;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    flags = reg_str2flags(interp, args.type);

    reg_init_args(args.args);

    switch(flags) {
    case REG_CURSOR_NOTIFY:
	{
	    reg_cursor_notify cn;
	    cursor_t *cp;
	    int cnum = args.cnum;
	    int abspos, job;

	    cp = find_contig_cursor(args.io, &cnum, atoi(reg_get_arg("id")));
	    if (!cp) {
#ifdef DEBUG
		verror(ERR_WARN, "contig_notify", "Unable to find cursor for "
		       "contig %d with id %d\n",
		       cnum, atoi(reg_get_arg("id")));
#endif
		return TCL_OK;
	    }

	    abspos = atoi(reg_get_arg("abspos"));
	    job = cjob_str2int(interp, reg_get_arg("job"));
#if 0
	    /* FIXME */
	    if (abspos == cp->abspos /* && (job & CURSOR_MOVE) */) {
		/* puts("No movement"); */
		break;
	    }
#endif

	    cn.job = REG_CURSOR_NOTIFY;
	    cn.cursor = cp;
	    cp->job = job;
	    cp->seq = atoi(reg_get_arg("seq"));
	    cp->pos = atoi(reg_get_arg("pos"));
	    cp->abspos = abspos;
	    cp->sent_by = atoi(reg_get_arg("sent_by"));

	    contig_notify(args.io, args.cnum, (reg_data *)&cn);
	    break;
	}

    case REG_BUFFER_START:
	{
	    reg_buffer_start rs;
	    rs.job = REG_BUFFER_START;
	    contig_notify(args.io, args.cnum, (reg_data *)&rs);
	    break;
	}

    case REG_BUFFER_END:
	{
	    reg_buffer_start rs;
	    rs.job = REG_BUFFER_END;
	    contig_notify(args.io, args.cnum, (reg_data *)&rs);
	    break;
	}

    case REG_NOTE:
	{
	    reg_note rn;
	    char *task;

	    rn.job = REG_NOTE;
	    rn.note = atoi(reg_get_arg("note"));
	    task = reg_get_arg("task");
	    if (strcmp(task, "CREATE") == 0) {
		rn.task = REG_NOTE_CREATE;
	    } else if (strcmp(task, "DELETE") == 0) {
		rn.task = REG_NOTE_DELETE;
	    } else {
		rn.task = REG_NOTE_EDIT;
	    }

	    contig_notify(args.io, args.cnum, (reg_data *)&rn);
	    break;
	}

    default:
	verror(ERR_WARN, "contig_notify", "Unknown event type '%s'",
	       args.type);
    }

    return TCL_OK;
}

typedef struct {
    GapIO *io;
    int cnum;
    int ref;
    int id;
} c_ref_arg;

static char *arg_names[100];
static char *arg_values[100];
static int arg_count;

/*
 * Process 'args' into arg_names and arg_values pairs.
 * Input format is "{name val} {name val} ..."
 */
static void reg_init_args(char *args)
{
    static char targs_a[8192];
    char *targs = targs_a;
    int bcount;

    strncpy(targs, args, 8191);
    targs[8191] = 0;

    arg_count = 0;
    while (*targs) {
	/* Skip past initial { */
	while (*targs == ' ' || *targs == '{')
	    targs++;

	/* Point to name */
	arg_names[arg_count] = targs;
	while (*++targs != ' ');
	*targs++ = 0;

	/* Point to val */
	arg_values[arg_count] = targs;
	for (bcount = 1; bcount > 0; targs++) {
	    if (*targs == '{')
		bcount++;
	    else if (*targs == '}')
		bcount--;
	}
	*(targs-1)=0;

	/* Skip any prevailing '}' */
	while (*targs && *targs == '}')
	    targs++;

	arg_count++;
    }
}

static char *reg_get_arg(char *name)
{
    int i;

    for (i = 0; i < arg_count; i++) {
	if (strcmp(name, arg_names[i]) == 0)
	    return arg_values[i];
    }

    return NULL;
}

typedef struct {
    GapIO *io;
    int cursorid;
    int cnum;
} qc_arg;

/*
 * An arbitray Tcl interface to contig event notification
 */
int tk_query_cursor(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv)
{
    qc_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(qc_arg, io)},
	{"-cursorid", ARG_INT, 1, NULL, offsetof(qc_arg, cursorid)},
	{"-cnum",     ARG_INT, 1, "0",  offsetof(qc_arg, cnum)},
        {NULL,        0,       0, NULL, 0}
    };
    cursor_t *gc;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    if (NULL == (gc = find_contig_cursor(args.io, &args.cnum, args.cursorid)))
	return TCL_OK;

    vTcl_SetResult(interp,
		   "{id %d} {refs %d} {private %d} {abspos %d} {contig %d}",
		   gc->id, gc->refs, gc->private, gc->abspos, args.cnum);
    return TCL_OK;
}

typedef struct {
    GapIO *io;
    int cnum;
    int private;
    int sent_by;
} ccc_arg;

/*
 * A tcl interface to create_contig_cursor().
 */
int tk_create_cursor(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv)
{
    ccc_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(ccc_arg, io)},
	{"-cnum",     ARG_INT, 1, NULL, offsetof(ccc_arg, cnum)},
	{"-private",  ARG_INT, 1, "0", offsetof(ccc_arg, private)},
	{"-sent_by",  ARG_INT, 1, "0", offsetof(ccc_arg, sent_by)},
        {NULL,        0,       0, NULL, 0}
    };
    cursor_t *cp;

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    cp = create_contig_cursor(args.io, args.cnum, args.private, args.sent_by);
    vTcl_SetResult(interp, "%d", cp->id);
    return TCL_OK;
}


typedef struct {
    GapIO *io;
    int cnum;
    int id;
    int private;
} dc_arg;

/*
 * A tcl interface to delete_contig_cursor().
 */
int tk_delete_cursor(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv)
{
    dc_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(dc_arg, io)},
	{"-cnum",     ARG_INT, 1, "0",  offsetof(dc_arg, cnum)},
	{"-id",       ARG_INT, 1, NULL, offsetof(dc_arg, id)},
	{"-private",  ARG_INT, 1, "0",  offsetof(dc_arg, private)},
        {NULL,        0,       0, NULL, 0}
    };

    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    delete_contig_cursor(args.io, args.cnum, args.id, args.private);
    return TCL_OK;
}


/*
 * Increment or decrement the cursor reference counter
 */
int tk_cursor_ref(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv)
{
    c_ref_arg args;
    reg_cursor_notify cn;
    cursor_t *cursor;

    cli_args a[] = {
        {"-io",   ARG_IO,  1, NULL, offsetof(c_ref_arg, io)},
	{"-cnum", ARG_INT, 1, NULL, offsetof(c_ref_arg, cnum)},
	{"-ref",  ARG_INT, 1, NULL, offsetof(c_ref_arg, ref)},
	{"-id",   ARG_INT, 1, NULL, offsetof(c_ref_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    if (-1 == gap_parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    cursor = find_contig_cursor(args.io, &args.cnum, args.id);
    if (!cursor) {
#ifdef DEBUG
	verror(ERR_WARN, "contig_notify", "Unable to find cursor for "
	       "contig %d with id %d\n",
	       args.cnum, args.id);
#endif
	return TCL_OK;
    }

    cursor->refs += args.ref;

    cn.job = REG_CURSOR_NOTIFY;
    cn.cursor = cursor;
    cn.cursor->job = CURSOR_MOVE;

    contig_notify(args.io, args.cnum, (reg_data *)&cn);
    return TCL_OK;
}
