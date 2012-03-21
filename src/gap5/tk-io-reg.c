#include <stdio.h>
#include <tk.h>
#include <string.h>

#include "gap4_compat.h"
#include "tk-io-reg.h"
#include "gap_cli_arg.h"
#include "newgap_cmds.h" /* GetInterp() */
#include "newgap_structs.h" /* contig_arg */
#include "misc.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"
#include "canvas_box.h"
#include "io_utils.h" /* get_gel_num, rnumtocnum, chain_left */

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
		    int objc, Tcl_Obj *CONST objv[]) {
    result_name_t *res;
    int i, nres;
    rn_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rn_arg, io)},
        {NULL,        0,       0, NULL, 0}
    };
    char buf[1024];
    Tcl_DString ds;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    res = result_names(args.io, &nres);

    Tcl_DStringInit(&ds);
    for (i = 0; i < nres; i++) {
	sprintf(buf, "%"PRIrec" %d {%s}",
		res[i].contig, res[i].id, res[i].name);
	Tcl_DStringAppendElement(&ds, buf);
    }

    Tcl_DStringResult(interp, &ds);
    if (res)
	free(res);

    return TCL_OK;
}


/*
 * Tk interface to register_id() - all this work for such a miniscule 2 line
 * function!!!
 */
int tk_register_id(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[]) {
    vTcl_SetResult(interp, "%d", register_id());

    return TCL_OK;
}


typedef struct {
    GapIO *io;
    int id;
} rt_arg;

/*
 * Tk interface to result_time()
 */
int tk_result_time(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[]) {
    rt_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rt_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(rt_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    Tcl_SetResult(interp, result_time(args.io, args.id),
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
		     int objc, Tcl_Obj *CONST objv[]) {
    go_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(go_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(go_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_get_ops ro;
    int l;
    char *ops;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
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
		     int objc, Tcl_Obj *CONST objv[]) {
    iop_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(iop_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(iop_arg, id)},
	{"-option",   ARG_INT, 1, NULL, offsetof(iop_arg, option)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_invoke_op inv;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
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
		     int objc, Tcl_Obj *CONST objv[]) {
    contig_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(contig_arg, io)},
        {"-contig",   ARG_REC, 1, NULL, offsetof(contig_arg, contig)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_length rl;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    rl.job = REG_LENGTH;

    rl.length = args.contig ? io_clength(args.io, args.contig) : 0;
    contig_notify(args.io, args.contig, (reg_data *)&rl);

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

    sprintf(buf, "result_list_update %s", io_obj_as_string(io));

    Tcl_Eval(GetInterp(), buf);
    result_pending = 0;
}

void update_results(GapIO *io) {
    if (!result_pending) {
	result_pending = 1;

	Tcl_DoWhenIdle(update_results_, (ClientData)io);
    }
}

typedef struct {
    char *result;
    char *colour;
    char *csplot;
    int width;
} conf_arg;


typedef struct {
    GapIO *io;
    char *msg;
} qd_arg;

/*
 * Attempts to shut down all active displays.
 * As used by alter relationships, assembly, etc.
 *
 * Args: quit_display handle function_name
 */
int tcl_quit_displays(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[]) {
    qd_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(qd_arg, io)},
	{"-msg",      ARG_STR, 1, "?",  offsetof(qd_arg, msg)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_quit rq;
    int ret = 1;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    rq.job = REG_QUIT;
    rq.lock = REG_LOCK_WRITE;
    contig_notify(args.io, 0, (reg_data *)&rq);
    if (!(rq.lock & REG_LOCK_WRITE)) {
	ret = 0;
    }

    if (!ret) {
	verror(ERR_WARN, args.msg, "Database busy");
	verror(ERR_WARN, args.msg, "Please shut down editing displays");
	//busy_dialog(args.io, cnum);
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
		     int objc, Tcl_Obj *CONST objv[]) {
    rd_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rd_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(rd_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_delete rd;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    rd.job = REG_DELETE;
    result_notify(args.io, args.id, (reg_data *)&rd, 0);

    return TCL_OK;
}

/*
 * Send a delete request to a specific result
 */
int tk_result_quit(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[]) {
    rd_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rd_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(rd_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_delete rd;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    rd.job = REG_QUIT;
    result_notify(args.io, args.id, (reg_data *)&rd, 0);

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
		     int objc, Tcl_Obj *CONST objv[]) {
    notify_arg args;
    tg_rec r_num;
    int is_name;

    cli_args a[] = {
	{"-io",		ARG_IO,  1, NULL, offsetof(notify_arg, io)},
	{"-reading",	ARG_STR, 1, NULL, offsetof(notify_arg, reading)},
	{"-highlight",	ARG_INT, 1, NULL, offsetof(notify_arg, value)},
        {NULL,          0,       0, NULL, 0}
    };
    reg_highlight_read hr;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    if (args.reading[0] == '=' || args.reading[0] == '#') {
	is_name = GGN_ID;
    } else {
	is_name = GGN_NAME;
    }
    r_num = get_gel_num(args.io, args.reading, is_name);
    if (r_num <= 0) {
	verror(ERR_WARN, "reg_notify_highlight", "Unknown reading '%s'",
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
		     int objc, Tcl_Obj *CONST objv[]) {
    rd_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rd_arg, io)},
	{"-id",	      ARG_INT, 1, NULL, offsetof(rd_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    contig_reg_t *regs;
    int r=0;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    if ((regs = get_reg_by_id(args.io, args.id, NULL))) {
	if ((regs->type == REG_TYPE_FIJ) ||
	    (regs->type == REG_TYPE_READPAIR) ||
	    (regs->type == REG_TYPE_REPEAT) ||
	    (regs->type == REG_TYPE_CHECKASS) ||
	    (regs->type == REG_TYPE_OLIGO)) {
	    r=1;
	}
    }

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
    int ref;
} cr_type;

static void tk_contig_register_cmd(GapIO *io, tg_rec contig, void *fdata,
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
	} else if (strcmp(str_flags[i], "RENAME") == 0) {
	    flags |= REG_RENAME;
	} else if (strcmp(str_flags[i], "CHILD_EDIT") == 0) {
	    flags |= REG_CHILD_EDIT;
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
    if (flags & REG_RENAME)
	strcat(str, "RENAME ");
    if (flags & REG_CHILD_EDIT)
	strcat(str, "CHILD_EDIT ");
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
 * Used as a task parameter to REG_GENERIC. NB putting this here means this
 * code has internal knowledge of the various gap4 windows. For this reason
 * we also allow them to be passed numerically in.
 *
 * For tcl it's often easier to simply use TASK_GENERIC and pass in data
 * as an arbitrary string.
 */
static int reg_str2task(char *str) {
    /* Blank str => task is just TASK_GENERIC */
    if (str[0] == '\0')
	return TASK_GENERIC;

    /* Numerical versions */
    if (str[0] >= '0' && str[0] <= '9')
	return atoi(str);

    /* Or string symbolic ones */
    if (strcmp(str, "TASK_GENERIC") == 0)
	return TASK_GENERIC;
    if (strcmp(str, "TASK_CANVAS_SCROLLX") == 0)
	return TASK_CANVAS_SCROLLX;
    if (strcmp(str, "TASK_CANVAS_SCROLLY") == 0)
	return TASK_CANVAS_SCROLLY;
    if (strcmp(str, "TASK_CANVAS_SCROLLX") == 0)
	return TASK_CANVAS_SCROLLX;
    if (strcmp(str, "TASK_CANVAS_ZOOM") == 0)
	return TASK_CANVAS_ZOOM;
    if (strcmp(str, "TASK_CANVAS_CURSOR_X") == 0)
	return TASK_CANVAS_CURSOR_X;
    if (strcmp(str, "TASK_CANVAS_CURSOR_Y") == 0)
	return TASK_CANVAS_CURSOR_Y;
    if (strcmp(str, "TASK_CANVAS_CURSOR_DELETE") == 0)
	return TASK_CANVAS_CURSOR_DELETE;
    if (strcmp(str, "TASK_CANVAS_RESIZE") == 0)
	return TASK_CANVAS_RESIZE;
    if (strcmp(str, "TASK_CANVAS_WORLD") == 0)
	return TASK_CANVAS_WORLD;
    if (strcmp(str, "TASK_WINDOW_ADD") == 0)
	return TASK_WINDOW_ADD;
    if (strcmp(str, "TASK_WINDOW_DELETE") == 0)
	return TASK_WINDOW_DELETE;
    if (strcmp(str, "TASK_DISPLAY_TICKS") == 0)
	return TASK_DISPLAY_TICKS;
    if (strcmp(str, "TASK_DISPLAY_RULER") == 0)
	return TASK_DISPLAY_RULER;
    if (strcmp(str, "TASK_CONS_WORLD") == 0)
	return TASK_CONS_WORLD;
    if (strcmp(str, "TASK_CONS_JOIN") == 0)
	return TASK_CONS_JOIN;
    if (strcmp(str, "TASK_CONS_CURSOR_DELETE") == 0)
	return TASK_CONS_CURSOR_DELETE;
    if (strcmp(str, "TASK_CONS_ID") == 0)
	return TASK_CONS_ID;
#if 0
    if (strcmp(str, "TASK_EDITOR_SETCON") == 0)
	return TASK_EDITOR_SETCON;
    if (strcmp(str, "TASK_EDITOR_GETCON") == 0)
	return TASK_EDITOR_GETCON;
    if (strcmp(str, "TASK_RENZ_INFO") == 0)
	return TASK_RENZ_INFO;
    if (strcmp(str, "TASK_TEMPLATE_REDRAW") == 0)
	return TASK_TEMPLATE_REDRAW;
#endif

    return -1;
}

static char *reg_task2str(int task) {
    switch (task) {
    case TASK_GENERIC:
	return "TASK_GENERIC";
    case TASK_CANVAS_SCROLLX:
	return "TASK_CANVAS_SCROLLX";
    case TASK_CANVAS_SCROLLY:
	return "TASK_CANVAS_SCROLLY";
    case TASK_CANVAS_ZOOM:
	return "TASK_CANVAS_ZOOM";
    case TASK_CANVAS_CURSOR_X:
	return "TASK_CANVAS_CURSOR_X";
    case TASK_CANVAS_CURSOR_Y:
	return "TASK_CANVAS_CURSOR_Y";
    case TASK_CANVAS_CURSOR_DELETE:
	return "TASK_CANVAS_CURSOR_DELETE";
    case TASK_CANVAS_RESIZE:
	return "TASK_CANVAS_RESIZE";
    case TASK_CANVAS_WORLD:
	return "TASK_CANVAS_WORLD";
    case TASK_WINDOW_ADD:
	return "TASK_WINDOW_ADD";
    case TASK_WINDOW_DELETE:
	return "TASK_WINDOW_DELETE";
    case TASK_DISPLAY_TICKS:
	return "TASK_DISPLAY_TICKS";
    case TASK_DISPLAY_RULER:
	return "TASK_DISPLAY_RULER";
    case TASK_CONS_WORLD:
	return "TASK_CONS_WORLD";
    case TASK_CONS_JOIN:
	return "TASK_CONS_JOIN";
    case TASK_CONS_CURSOR_DELETE:
	return "TASK_CONS_CURSOR_DELETE";
    case TASK_CONS_ID:
	return "TASK_CONS_ID";
    }

    return "unknown";
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
    tg_rec contig;
    int id;
    char *command;
    char *flags;
    char *type;
} cr_arg;

int tk_contig_register(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[]) {
    cr_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(cr_arg, io)},
	{"-contig",   ARG_REC, 1, NULL, offsetof(cr_arg, contig)},
	{"-command",  ARG_STR, 1, NULL, offsetof(cr_arg, command)},
	{"-flags",    ARG_STR, 1, "",   offsetof(cr_arg, flags)},
	{"-type",     ARG_STR, 1, "TYPE_UNKNOWN", offsetof(cr_arg, type)},
        {NULL,        0,       0, NULL, 0}
    };
    int flags, type;
    cr_type *crt;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    if (NULL == (crt = (cr_type *)xmalloc(sizeof(*crt))))
	return TCL_ERROR;

    crt->interp = interp;
    crt->command = strdup(args.command);
    crt->id = register_id();
    crt->ref = 1;
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
		     int objc, Tcl_Obj *CONST objv[]) {
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

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
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
		ret |= contig_deregister(args.io, -args.id,
					 ptr->func, ptr->fdata);
		xfree(crt->command);
		crt->command = NULL;
		if (--crt->ref == 0) {
		    xfree(crt);
		}
	    }
	}

	xfree(uids);
    }
    vTcl_SetResult(interp, "%d", ret);

    if (regs)
      xfree(regs);

    return TCL_OK;
}

static void tk_contig_register_cmd(GapIO *io, tg_rec contig, void *fdata,
				   reg_data *jdata)
{
    Tcl_DString ds;
    cr_type *crt = (cr_type *)fdata;
    char *type = reg_flags2str(jdata->job);
    char buf[1024];

    Tcl_DStringInit(&ds);
    sprintf(buf, "{contig_num %"PRIrec"} ", contig);
    Tcl_DStringAppend(&ds, buf, -1);
    
    switch (jdata->job) {
    case REG_NUMBER_CHANGE:
	sprintf(buf, "{number %d}",
		jdata->number.number);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_JOIN_TO:
	sprintf(buf, "{contig %"PRIrec"} {offset %d}",
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

    case REG_RENAME:
	sprintf(buf, "{name %s}",
		jdata->rename.name);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_QUERY_NAME:
	/* 'line' element is returned */
	break;

    case REG_CHILD_EDIT:
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
	    //first = 0;
	}
	strcat(job, "}");
	sprintf(buf, "{id %d} {seq %"PRIrec"} {pos %d} {abspos %d} {refs %d} "
		"{sent_by %d} {job %s}",
		jdata->cursor_notify.cursor->id,
		jdata->cursor_notify.cursor->seq,
		jdata->cursor_notify.cursor->pos,
		jdata->cursor_notify.cursor->abspos,
		jdata->cursor_notify.cursor->refs,
		jdata->cursor_notify.cursor->sent_by,
		job);
	Tcl_DStringAppend(&ds, buf, -1);
	break;
    }
    case REG_REGISTER:
    case REG_DEREGISTER:
	sprintf(buf, "{id %d} {type %s} {contig %"PRIrec"}",
		jdata->c_register.id,
		reg_type2str(jdata->c_register.type),
		jdata->c_register.contig);
	Tcl_DStringAppend(&ds, buf, -1);
	break;

    case REG_HIGHLIGHT_READ:
	sprintf(buf, "{seq %"PRIrec"} {val %d}",
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

    case REG_GENERIC:
	Tcl_DStringAppendElement(&ds, reg_task2str(jdata->generic.task));
	Tcl_DStringAppendElement(&ds, jdata->generic.data);
	break;

    default:
	Tcl_DStringAppend(&ds, "{unknown unknown}", -1);
	break;
    }

    crt->ref++;
    sprintf(buf, "%d", crt->id);
    if (Tcl_VarEval(crt->interp, crt->command, " ", type, " ", buf, " ",
		    Tcl_DStringValue(&ds), NULL) != TCL_OK) {
	fprintf(stderr, "registration_callback: %s", Tcl_GetStringResult(crt->interp));
	verror(ERR_WARN, "registration_callback", "%s",
	       Tcl_GetStringResult(crt->interp));
    }

    /* Now for the return value bit */
    switch(jdata->job) {
	static char buf[1025];
    case REG_QUERY_NAME:
	sprintf(jdata->name.line, "%.80s", Tcl_GetStringResult(crt->interp));
	break;

    case REG_GET_OPS:
	{
	    int largc, i, p;
	    char **largv;

	    if (Tcl_SplitList(crt->interp, Tcl_GetStringResult(crt->interp),
			      &largc, &largv) == TCL_OK) {
		for (p=i=0; i<largc; i++) {
		    strcpy(&buf[p], largv[i]);
		    p += strlen(largv[i])+1;
		}
		buf[p] = buf[p+1] = 0;
		jdata->get_ops.ops = buf;

		Tcl_Free((char *)largv);
		break;
	    }
	}

    case REG_GET_LOCK:
	jdata->glock.lock = atoi(Tcl_GetStringResult(crt->interp));
	break;

    case REG_PARAMS:
	sprintf(buf, "%.1024s", Tcl_GetStringResult(crt->interp));
	jdata->params.string = buf;
	break;

    case REG_QUIT:
	/* True = ok to quit, false = don't quit */
	if (atoi(Tcl_GetStringResult(crt->interp)) == 0) {
	    jdata->glock.lock &= ~REG_LOCK_WRITE;
	}
	break;
    }

    if (--crt->ref == 0) {
	xfree(crt);
    }

    Tcl_DStringFree(&ds);
}

/*
 * Parses type and str (a tcl list) to fill out a reg_data struct.
 *
 * Returns 0 on successful parse
 *        -1 on error or no notification required
 */
int str2reg_data(Tcl_Interp *interp, GapIO *io, 
		 tg_rec cnum, char *type, char *str, reg_data *rd) {
    int itype = reg_str2flags(interp, type);
    rd->job = itype;

    reg_init_args(str);

    switch (itype) {
    case REG_CHILD_EDIT:
    case REG_BUFFER_START:
    case REG_BUFFER_END:
	break;

    case REG_CURSOR_NOTIFY:
    {
	cursor_t *cp = find_contig_cursor(io, cnum, atoi(reg_get_arg("id")));
	int abspos;
	if (!cp)
	    return -1;

	abspos = atoi(reg_get_arg("abspos"));
	if (abspos == cp->abspos) {
	    /* No movement so skip event */
	    return -1;
	}

	cp->abspos = abspos;
	cp->job = cjob_str2int(interp, reg_get_arg("job"));
	cp->seq = atorec(reg_get_arg("seq"));
	cp->pos = atoi(reg_get_arg("pos"));
	cp->sent_by = atoi(reg_get_arg("sent_by"));

	rd->cursor_notify.cursor = cp;

	break;
    }

    case REG_NOTE:
    {
	char *task;
	rd->note.note = atoi(reg_get_arg("note"));

	task = reg_get_arg("task");
	if (strcmp(task, "CREATE") == 0) {
	    rd->note.task = REG_NOTE_CREATE;
	} else if (strcmp(task, "DELETE") == 0) {
	    rd->note.task = REG_NOTE_DELETE;
	} else {
	    rd->note.task = REG_NOTE_EDIT;
	}

	break;
    }

    case REG_LENGTH:
	rd->length.length = io_clength(io, cnum);
	break;
	
    case REG_RENAME: {
	static char buf[1024];
	contig_t *c = cache_search(io, GT_Contig, cnum);
	if (c)
	    sprintf(buf, "%s", contig_get_name(&c));
	else
	    sprintf(buf, "?");
	rd->rename.name = buf;
	break;
    }
	
    case REG_QUERY_NAME:
    {
	static char buf[81];
	rd->name.line = buf;
	break;
    }

    case REG_GENERIC:
	rd->generic.task = reg_str2task(reg_get_arg("task"));
	switch (rd->generic.task) {
	case TASK_GENERIC:
	    rd->generic.data = reg_get_arg("data");
	    break;

	case TASK_CANVAS_SCROLLX:
	case TASK_CANVAS_SCROLLY:
	    rd->generic.data = reg_get_arg("data");
	    break;
	    
	case TASK_CANVAS_CURSOR_X: {
	    static int x;
	    x = atoi(reg_get_arg("data"));
	    rd->generic.data = (void *)&x;
	    break;
	  }

	case TASK_CANVAS_CURSOR_Y: {
	    static int y;
	    y = atoi(reg_get_arg("data"));
	    rd->generic.data = (void *)&y;
	    break;
	  }
	}
	break;

    case REG_HIGHLIGHT_READ:
	rd->highlight.seq = atorec(reg_get_arg("seq"));
	rd->highlight.val = atoi(reg_get_arg("highlight"));
	break;

    default:
	verror(ERR_WARN, "str2reg_data", "unsupported event type '%s'", type);
	return -1;
    }

    return 0;
}

typedef struct {
    GapIO *io;
    tg_rec cnum;
    char *type;
    char *args;
} cn_arg;

/*
 * An arbitray Tcl interface to contig event notification
 */
int tk_contig_notify(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[])
{
    cn_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(cn_arg, io)},
	{"-cnum",     ARG_REC, 1, NULL, offsetof(cn_arg, cnum)},
	{"-type",     ARG_STR, 1, NULL, offsetof(cn_arg, type)},
	{"-args",     ARG_STR, 1, NULL, offsetof(cn_arg, args)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_data rd;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    if (-1 == str2reg_data(interp, args.io, args.cnum, args.type,
			   args.args, &rd))
	return TCL_OK;

    contig_notify(args.io, args.cnum, &rd);

    switch(rd.job) {
    case REG_QUERY_NAME:
	Tcl_SetResult(interp, rd.name.line, TCL_VOLATILE);
	break;
    }

    return TCL_OK;
}

typedef struct {
    GapIO *io;
    int    id;
    char  *type;
    char  *args;
} rnot_arg;

int tk_result_notify(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[])
{
    rnot_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(rnot_arg, io)},
	{"-id",       ARG_INT, 1, NULL, offsetof(rnot_arg, id)},
	{"-type",     ARG_STR, 1, NULL, offsetof(rnot_arg, type)},
	{"-args",     ARG_STR, 1, NULL, offsetof(rnot_arg, args)},
        {NULL,        0,       0, NULL, 0}
    };
    reg_data rd;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;
    
    if (-1 == str2reg_data(interp, args.io, 0, args.type, args.args, &rd))
	return TCL_OK;

    result_notify(args.io, args.id, &rd, 0);

    return TCL_OK;
}

typedef struct {
    GapIO *io;
    tg_rec cnum;
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
	while (*++targs && *targs != ' ');
	*targs++ = 0;

	/* Point to val */
	while (*targs && *targs == ' ')
	    targs++;
	if (*targs == '{') {
	    arg_values[arg_count] = targs+1;
	    bcount = 1;
	    while (bcount && *++targs) {
		if (*targs == '{')
		    bcount++;
		else if (*targs == '}')
		    bcount--;
	    }
	} else {
	    arg_values[arg_count] = targs;
	    while (*++targs && *targs != ' ' && *targs != '}');
	}
	*targs++ = 0;

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

    return "";
}

typedef struct {
    GapIO *io;
    int cursorid;
    tg_rec cnum;
} qc_arg;

/*
 * An arbitray Tcl interface to contig event notification
 */
int tk_query_cursor(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[])
{
    qc_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(qc_arg, io)},
	{"-cursorid", ARG_INT, 1, NULL, offsetof(qc_arg, cursorid)},
	{"-cnum",     ARG_REC, 1, "0",  offsetof(qc_arg, cnum)},
        {NULL,        0,       0, NULL, 0}
    };
    cursor_t *gc;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    if (NULL == (gc = find_contig_cursor(args.io, args.cnum, args.cursorid)))
	return TCL_OK;

    vTcl_SetResult(interp,
		   "{id %d} {refs %d} {private %d} {abspos %d} {contig %"PRId64"}",
		   gc->id, gc->refs, gc->private, gc->abspos, args.cnum);
    return TCL_OK;
}

typedef struct {
    GapIO *io;
    tg_rec cnum;
    int private;
    int sent_by;
} ccc_arg;

/*
 * A tcl interface to create_contig_cursor().
 */
int tk_create_cursor(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[])
{
    ccc_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(ccc_arg, io)},
	{"-cnum",     ARG_REC, 1, NULL, offsetof(ccc_arg, cnum)},
	{"-private",  ARG_INT, 1, "0", offsetof(ccc_arg, private)},
	{"-sent_by",  ARG_INT, 1, "0", offsetof(ccc_arg, sent_by)},
        {NULL,        0,       0, NULL, 0}
    };
    cursor_t *cp;

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    cp = create_contig_cursor(args.io, args.cnum, args.private, args.sent_by);
    vTcl_SetResult(interp, "%d", cp->id);
    return TCL_OK;
}


typedef struct {
    GapIO *io;
    tg_rec cnum;
    int id;
    int private;
} dc_arg;

/*
 * A tcl interface to delete_contig_cursor().
 */
int tk_delete_cursor(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[])
{
    dc_arg args;
    cli_args a[] = {
        {"-io",       ARG_IO,  1, NULL, offsetof(dc_arg, io)},
	{"-cnum",     ARG_REC, 1, "0",  offsetof(dc_arg, cnum)},
	{"-id",       ARG_INT, 1, NULL, offsetof(dc_arg, id)},
	{"-private",  ARG_INT, 1, "0",  offsetof(dc_arg, private)},
        {NULL,        0,       0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    delete_contig_cursor(args.io, args.cnum, args.id, args.private);
    return TCL_OK;
}


/*
 * Increment or decrement the cursor reference counter
 */
int tk_cursor_ref(ClientData clientData, Tcl_Interp *interp,
		  int objc, Tcl_Obj *CONST objv[])
{
    c_ref_arg args;
    reg_cursor_notify cn;
    cursor_t *cursor;

    cli_args a[] = {
        {"-io",   ARG_IO,  1, NULL, offsetof(c_ref_arg, io)},
	{"-cnum", ARG_REC, 1, NULL, offsetof(c_ref_arg, cnum)},
	{"-ref",  ARG_INT, 1, NULL, offsetof(c_ref_arg, ref)},
	{"-id",   ARG_INT, 1, NULL, offsetof(c_ref_arg, id)},
        {NULL,        0,       0, NULL, 0}
    };
    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    cursor = find_contig_cursor(args.io, args.cnum, args.id);
    if (!cursor) {
	verror(ERR_WARN, "contig_notify", "Unable to find cursor for "
	       "contig %"PRIrec" with id %d\n",
	       args.cnum, args.id);
	return TCL_OK;
    }

    cursor->refs += args.ref;

    cn.job = REG_CURSOR_NOTIFY;
    cn.cursor = cursor;
    cn.cursor->job = CURSOR_MOVE;

    contig_notify(args.io, args.cnum, (reg_data *)&cn);
    return TCL_OK;
}


/*
 * Interface to contig_lock_write()
 */
int tk_contig_lock_write(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[])
{
    cn_arg args;

    cli_args a[] = {
        {"-io",   ARG_IO,  1, NULL, offsetof(cn_arg, io)},
	{"-cnum", ARG_REC, 1, NULL, offsetof(cn_arg, cnum)},
        {NULL,        0,       0, NULL, 0}
    };
    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
        return TCL_ERROR;

    vTcl_SetResult(interp, "%d", contig_lock_write(args.io, args.cnum));

    return TCL_OK;
}
