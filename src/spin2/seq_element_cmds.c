#include <tcl.h>
#include <string.h>
#include <tclXkeylist.h>

#include "text_output.h"
#include "cli_arg.h"
#include "xalloc.h"
#include "seq_element.h"
#include "seq_results.h"
#include "container.h"
#include "spin_cli_arg.h"
#include "tcl_utils.h"

#define setli(s,x) \
do { \
    Tcl_Obj *o = Tcl_NewIntObj((s).x); \
    TclX_KeyedListSet(interp, klist, w(#x), o); \
} while (0)

#define setls(s,x) \
do { \
    Tcl_Obj *o = Tcl_NewStringObj((s).x, -1); \
    TclX_KeyedListSet(interp, klist, w(#x), o); \
} while (0)

typedef struct {
    int container_id;
    int seq_id_h;
    int seq_id_v;
    char *plot_type;
    int frame;
    char *element_type;
} get_element_arg;

int tcl_get_element_name(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{    
    get_element_arg args;
    int plot_type;
    seq_id_dir *seq_ids = NULL;
    int num = 0;
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *klist;
    element_info e_info;
    int new_container = 0;

    cli_args a[] = {
	{"-container_id", ARG_INT, 1, "-1", offsetof(get_element_arg, container_id)},
	{"-seq_id_h", ARG_INT, 1, "-1", offsetof(get_element_arg, seq_id_h)},
	{"-seq_id_v", ARG_INT, 1, "-1", offsetof(get_element_arg, seq_id_v)},
	{"-plot_type", ARG_STR, 1, "", offsetof(get_element_arg, plot_type)},
	{"-frame",  ARG_INT, 1, "0",  offsetof(get_element_arg, frame)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv)) {
	printf("trouble parsing args \n");
	return TCL_ERROR;
    }
    
    if (!list)
	return TCL_ERROR;

    Tcl_IncrRefCount(list);

    if (strcmp(args.plot_type, "") == 0) {
	plot_type = -1;
    } else {
	if (strcmp(args.plot_type, "GENESEARCH") == 0) {
	    plot_type = SEQ_TYPE_GENESEARCH;
	} else if (strcmp(args.plot_type, "STRINGSEARCH") == 0) {
	    plot_type = SEQ_TYPE_STRINGSEARCH;
	} else if (strcmp(args.plot_type, "RESTRICTION") == 0) {
	    plot_type = SEQ_TYPE_RESTRICTION;
	} else if (strcmp(args.plot_type, "BASECOMP") == 0) {
	    plot_type = SEQ_TYPE_BASECOMP;
	} else if (strcmp(args.plot_type, "CODONPREF") == 0) {
	    plot_type = SEQ_TYPE_CODONPREF;
	} else if (strcmp(args.plot_type, "AUTHOR") == 0) {
	    plot_type = SEQ_TYPE_AUTHOR;
	} else if (strcmp(args.plot_type, "BASEBIAS") == 0) {
	    plot_type = SEQ_TYPE_BASEBIAS;
	} else if (strcmp(args.plot_type, "TRNA") == 0) {
	    plot_type = SEQ_TYPE_TRNA;
	} else if (strcmp(args.plot_type, "STOPCODON") == 0) {
	    plot_type = SEQ_TYPE_STOPCODON;
	} else if (strcmp(args.plot_type, "STARTCODON") == 0) {
	    plot_type = SEQ_TYPE_STARTCODON;
	} else if (strcmp(args.plot_type, "SPLICE") == 0) {
	    plot_type = SEQ_TYPE_SPLICE;
	} else if (strcmp(args.plot_type, "WTMATRIXSEARCH") == 0) {
	    plot_type = SEQ_TYPE_WTMATRIXSEARCH;
	} else if (strcmp(args.plot_type, "DOT_PLOT") == 0) {
	    plot_type = SEQ_TYPE_DOT_PLOT;
	} else if (strcmp(args.plot_type, "RULER") == 0) {
	    plot_type = SEQ_TYPE_RULER;
	} else {
	    verror(ERR_WARN, "get_element_name", "Unrecognised type");
	    return TCL_OK;
	}
    }

    if (args.container_id > -1) {
	int element_id;
	char *e_win;
	container *c;

	c = get_container(args.container_id);

	if (-1 == (element_id = new_element(interp, &e_win))) {
	    verror(ERR_WARN, "get_element_name", "Unable to get new element");
	    return TCL_OK;
	}

	e_info.element_id = element_id;
	e_info.element_win = strdup(e_win);
	e_info.container_id = args.container_id;
	e_info.container_win = strdup(c->win);
	e_info.orientation = HORIZONTAL;
	e_info.e_win_to = NULL;
    } else {
	/* force a new container to be created if container_id == -2 */
	if (args.container_id == -2)
	    new_container = 1;

	if (args.seq_id_h != -1)
	    num++;
	if (args.seq_id_v != -1)
	    num++;
	
	if (NULL == (seq_ids = (seq_id_dir *)xmalloc(num * sizeof(seq_id_dir))))
	    return -1;
	
	num = 0;
	if (args.seq_id_h != -1) {
	    seq_ids[num].seq_id = args.seq_id_h;
	    seq_ids[num].direction = HORIZONTAL;
	    num++;
	}
	if (args.seq_id_v != -1) {
	    seq_ids[num].seq_id = args.seq_id_v;
	    seq_ids[num].direction = VERTICAL;
	}
	
	e_info = get_element_info(interp, seq_ids, num, plot_type, 
				  args.frame, new_container); 
    }

    klist = TclX_NewKeyedListObj();
    if (!klist)
	return TCL_ERROR;

    setli(e_info, element_id);
    setli(e_info, container_id);
    setls(e_info, element_win);
    setls(e_info, container_win);
    setli(e_info, orientation);
    setls(e_info, e_win_to);

    /* Set the result and return */
    Tcl_SetObjResult(interp, klist);
    Tcl_DecrRefCount(list);
    xfree(seq_ids);

    return TCL_OK;
}

int tcl_get_element_type(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{    
    element *e; 
    int type;

    if (argc != 2) {
	Tcl_SetResult(interp,
		      "wrong # args: should be \"get_element_type element_id",
		      TCL_STATIC);
	return TCL_ERROR;
    }
    e = get_element(atoi(argv[1]));
    type = get_element_type(e);

    vTcl_SetResult(interp, "%d", type);    
    return TCL_OK;
}
int tcl_get_element_gr_type(ClientData clientData, 
			    Tcl_Interp *interp, 
			    int argc, 
			    char *argv[]) 
{    
    element *e; 
    int type;

    if (argc != 2) {
	Tcl_SetResult(interp,
		      "wrong # args: should be \"get_element_type element_id",
		      TCL_STATIC);
	return TCL_ERROR;
    }
    e = get_element(atoi(argv[1]));
    type = get_element_gr_type(e);

    vTcl_SetResult(interp, "%d", type);    
    return TCL_OK;
}

typedef struct {
    int element_id;
    int seq_id_h;
    int seq_id_v;
} resultid_arg;

int tcl_seq_find_result_id(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc, 
			   char *argv[]) 
{    
    resultid_arg args;
    element *e;
    int result_id = -1;
    int i;
    seq_result *s_result;

    cli_args a[] = {
	{"-element_id", ARG_INT, 1, NULL, offsetof(resultid_arg, element_id)},
	{"-seq_id_h",  ARG_INT, 1, "-1", offsetof(resultid_arg, seq_id_h)},
	{"-seq_id_v",  ARG_INT, 1, "-1", offsetof(resultid_arg, seq_id_v)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv)) {
	printf("trouble parsing args \n");
	return TCL_ERROR;
    }

    e = get_element(args.element_id);
    for (i = 0; i < e->num_results; i++) {
	s_result = seq_id_to_result(e->results[i]->result_id);
	if (s_result->seq_id[HORIZONTAL] == args.seq_id_h &&
	    s_result->seq_id[VERTICAL] == args.seq_id_v)
	    result_id = e->results[i]->result_id;
    }

    vTcl_SetResult(interp, "%d", result_id);
    return TCL_OK;
}

typedef struct {
    int seq_id;
} seq_id_arg;

/* given a seq_id, find the first column associated with it */
int tcl_seq_find_column(ClientData clientData, 
			Tcl_Interp *interp, 
			int argc, 
			char *argv[]) 
{    
    seq_id_arg args;
    int row;

    cli_args a[] = {
	{"-seq_id",  ARG_INT, 1, "-1", offsetof(seq_id_arg, seq_id)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv)) {
	printf("trouble parsing args \n");
	return TCL_ERROR;
    }

    row = seq_find_column(interp, args.seq_id);

    vTcl_SetResult(interp, "%d", row);
    return TCL_OK;
}


typedef struct {
    int element_id;
    int direction;
} edcursor_arg;

int tcl_seq_find_cursor(ClientData clientData, 
			Tcl_Interp *interp, 
			int argc, 
			char *argv[]) 
{    
    edcursor_arg args;
    int seq_id, cursor_id;

    cli_args a[] = {
	{"-element_id", ARG_INT, 1, NULL, offsetof(edcursor_arg, element_id)},
	{"-direction",  ARG_INT, 1, "-1", offsetof(edcursor_arg, direction)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv)) {
	printf("trouble parsing args \n");
	return TCL_ERROR;
    }

    seq_id = -1;
    cursor_id = -1;
    if (args.direction == -1) 
	args.direction = HORIZONTAL;

    seq_find_cursor(args.element_id, args.direction, &seq_id, &cursor_id);

    vTcl_SetResult(interp, "%d %d", cursor_id, seq_id);
    return TCL_OK;
}

typedef struct {
    int element_id;
    int cursor_id;
    int pos;
    int direction;
} move_cursor_arg;

int tcl_canvas_move_cursor(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc, 
			   char *argv[]) 
{    
    move_cursor_arg args;

    cli_args a[] = {
	{"-element_id", ARG_INT, 1, NULL, offsetof(move_cursor_arg, element_id)},
	{"-pos",  ARG_INT, 1, NULL, offsetof(move_cursor_arg, pos)},
	{"-cursor_id",  ARG_INT, 1, NULL, offsetof(move_cursor_arg, cursor_id)},
	{"-direction", ARG_INT, 1, "-1", offsetof(move_cursor_arg, direction)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv)) {
	printf("trouble parsing args \n");
	return TCL_ERROR;
    }

    if (args.direction == -1) 
	args.direction = HORIZONTAL;

    seq_canvas_move_cursor(args.element_id, args.cursor_id, args.pos, args.direction);
    return TCL_OK;
}

typedef struct {
    int new_c_id;
    char *new_c_win;
    int old_c_id;
    int new_e_id;
    char *new_e_win;
    int new_orientation;
    int new_crosshair;
    int old_e_id;
    int result_id;
    char *element_type;
    char *job;
} update_c_arg;

int tcl_update_container(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{    
    update_c_arg args;

    cli_args a[] = {
	{"-new_container_id", ARG_INT, 1, "-1", offsetof(update_c_arg, new_c_id)},
	{"-new_container_win", ARG_STR, 1, "", offsetof(update_c_arg, new_c_win)},
	{"-old_container_id", ARG_INT, 1, "-1", offsetof(update_c_arg, old_c_id)},
	{"-new_element_id", ARG_INT, 1, "-1", offsetof(update_c_arg, new_e_id)},
	{"-new_element_win", ARG_STR, 1, "", offsetof(update_c_arg, new_e_win)},
	{"-new_orientation", ARG_INT, 1, "", offsetof(update_c_arg, new_orientation)},
	{"-new_crosshair", ARG_INT, 1, "", offsetof(update_c_arg, new_crosshair)},
	{"-old_element_id", ARG_INT, 1, "-1", offsetof(update_c_arg, old_e_id)},
	{"-result_id", ARG_INT, 1, "-1", offsetof(update_c_arg, result_id)},	
	{"-element_type", ARG_STR, 1, "CANVAS", offsetof(update_c_arg, element_type)},	
	{"-job", ARG_STR, 1, NULL, offsetof(update_c_arg, job)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv)) {
	printf("update_container trouble parsing args \n");
	return TCL_ERROR;
    }

    update_container(interp, args.new_c_id, args.new_c_win, args.old_c_id, 
		     args.new_e_id, args.new_e_win, args.new_orientation,
		     args.new_crosshair, args.old_e_id,
		     args.result_id, args.element_type, args.job);
    return TCL_OK;
}
typedef struct {
    int element_id;
    int result_id;
    int line_width;
    char *colour;
} result_arg;

int tcl_config_result(ClientData clientData, 
		      Tcl_Interp *interp, 
		      int argc, 
		      char *argv[]) 
{    
    result_arg args;

    cli_args a[] = {
	{"-element_id", ARG_INT, 1, NULL, offsetof(result_arg, element_id)},
	{"-result_id",  ARG_INT, 1, NULL, offsetof(result_arg, result_id)},
	{"-width", ARG_INT, 1, NULL, offsetof(result_arg, line_width)},
	{"-fill",     ARG_STR, 1, NULL, offsetof(result_arg, colour)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    config_result(args.element_id, args.result_id, args.line_width, args.colour);
    return TCL_OK;
}


void
Seq_Element_Init(Tcl_Interp *interp) {

    Tcl_CreateCommand(interp, "get_element_name", tcl_get_element_name,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "get_element_type", tcl_get_element_type,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "get_element_gr_type", tcl_get_element_gr_type,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "seq_find_cursor", tcl_seq_find_cursor, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_find_result_id", tcl_seq_find_result_id, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_find_column", tcl_seq_find_column, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "canvas_move_cursor", tcl_canvas_move_cursor, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "update_container", tcl_update_container, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "config_result", 
		      tcl_config_result, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
}
