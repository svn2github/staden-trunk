#include <tcl.h>
#include <string.h>
#include <tclXkeylist.h>

#include "xalloc.h"
#include "text_output.h"
#include "container.h"
#include "cli_arg.h"
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

#define setld(s,x) \
do { \
    Tcl_Obj *o = Tcl_NewDoubleObj((s).x); \
    TclX_KeyedListSet(interp, klist, w(#x), o); \
} while (0)

typedef struct {
    int container_id;
    int column;
    char *command;
} scroll_x_arg;

int tcl_container_scroll_x(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc, 
			   char *argv[]) 
{    
    scroll_x_arg args;

    cli_args a[] = {
	{"-container_id", ARG_INT, 1, NULL, offsetof(scroll_x_arg, container_id)},
	{"-column", ARG_INT, 1, NULL, offsetof(scroll_x_arg, column)},
	{"-command",ARG_STR, 1, "",   offsetof(scroll_x_arg, command)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    container_scroll_x(interp, args.container_id, args.column, args.command);

    return TCL_OK;
}

typedef struct {
    int container_id;
    int row;
    char *command;
} scroll_y_arg;

int tcl_container_scroll_y(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc, 
			   char *argv[]) 
{    
    scroll_y_arg args;

    cli_args a[] = {
	{"-container_id", ARG_INT, 1, NULL, offsetof(scroll_y_arg, container_id)},
	{"-row", ARG_INT, 1, NULL, offsetof(scroll_y_arg, row)},
	{"-command",ARG_STR, 1, "",   offsetof(scroll_y_arg, command)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    container_scroll_y(interp, args.container_id, args.row, args.command);

    return TCL_OK;
}

typedef struct {
    int container_id;
    float amount;
    int x1;
    int y1;
    int x2;
    int y2;
} zoom_arg;

int tcl_container_zoom(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{    
    zoom_arg args;

    cli_args a[] = {
	{"-container_id", ARG_INT, 1, NULL, offsetof(zoom_arg, container_id)},
	{"-amount",    ARG_FLOAT, 1, "-1.0", offsetof(zoom_arg, amount)},
	{"-x1",        ARG_INT, 1, "-1", offsetof(zoom_arg, x1)},
	{"-y1",        ARG_INT, 1, "-1", offsetof(zoom_arg, y1)},
	{"-x2",        ARG_INT, 1, "-1", offsetof(zoom_arg, x2)},
	{"-y2",        ARG_INT, 1, "-1", offsetof(zoom_arg, y2)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (args.amount == -1 && args.x1 == -1 && args.y1 == -1 && args.x2 == -1 
	&& args.y2 == -1) {
	container_zoomback(interp, args.container_id);
    } else {
	container_zoom(interp, args.container_id, args.amount, args.x1, args.y1,
		       args.x2, args.y2);
    }
    return TCL_OK;
}

typedef struct {
    int container_id;
} cid_arg;

typedef struct {
    int element_id;
} eid_arg;

int tcl_element_resize(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{    
    eid_arg args;

    cli_args a[] = {
	{"-element_id", ARG_INT, 1, NULL, offsetof(eid_arg, element_id)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    element_resize(interp, args.element_id);

    return TCL_OK;
}

typedef struct {
    int element_id;
    int x;
    int y;
} xhair_arg;

int tcl_draw_container_crosshair(ClientData clientData, 
				 Tcl_Interp *interp, 
				 int argc, 
				 char *argv[]) 
{    
    xhair_arg args;

    cli_args a[] = {
	{"-element_id", ARG_INT, 1, NULL, offsetof(xhair_arg, element_id)},
	{"-x", ARG_INT, 1, NULL, offsetof(xhair_arg, x)},
	{"-y", ARG_INT, 1, NULL,   offsetof(xhair_arg, y)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    draw_container_crosshair(interp, args.element_id, args.x, args.y);

    return TCL_OK;
}

typedef struct {
    int element_id;
    int x;
    int y;
} world_arg;

int tcl_pixel_to_world(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{    
    world_arg args;
    double wx, wy;
    element *e;

    cli_args a[] = {
	{"-element_id", ARG_INT, 1, NULL, offsetof(world_arg, element_id)},
	{"-x", ARG_INT, 1, "0", offsetof(world_arg, x)},
	{"-y", ARG_INT, 1, "0", offsetof(world_arg, y)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    e = get_element(args.element_id);
    pixel_to_world(e->pixel, args.x, args.y, &wx, &wy);

    vTcl_SetResult(interp, "%f %f", wx, wy);
    return TCL_OK;
}

int tcl_delete_container_crosshair(ClientData clientData, 
				   Tcl_Interp *interp, 
				   int argc, 
				   char *argv[]) 
{    
    eid_arg args;

    cli_args a[] = {
	{"-element_id", ARG_INT, 1, NULL, offsetof(eid_arg, element_id)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    delete_container_crosshair(interp, args.element_id);

    return TCL_OK;
}

int tcl_container_shutdown(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc, 
			   char *argv[]) 
{    
    cid_arg args;

    cli_args a[] = {
	{"-container_id", ARG_INT, 1, NULL, offsetof(cid_arg, container_id)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    container_start_shutdown(args.container_id);
    return TCL_OK;
}

/*
 * FIXME: to do - only do CANVAS & SEQ_EDITOR for now
 * number of elements of type -type 
 */
int tcl_num_container_elements(ClientData clientData, 
			       Tcl_Interp *interp, 
			       int argc, 
			       char *argv[]) 
{    
    cid_arg args;
    container *c;
    int i, j, num;

    cli_args a[] = {
	{"-container_id", ARG_INT, 1, NULL, offsetof(cid_arg, container_id)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    c = get_container(args.container_id);
    num = 0;
    for (i = 0; i < c->num_rows; i++) {
	for (j = 0; j < c->num_columns; j++) {
	    if (c->matrix[i][j]) {
		if (c->matrix[i][j]->type == CANVAS) {
		    num++;
		} else if (c->matrix[i][j]->type == SEQ_EDITOR) {
		    num++;
		}
	    }
	}
    }
    vTcl_SetResult(interp, "%d", num);
    return TCL_OK;
}

/*
 * FIXME: to do - only do CANVAS & SEQ_EDITOR for now
 * number of elements of type -type 
 */
int tcl_num_container_results(ClientData clientData, 
			      Tcl_Interp *interp, 
			      int argc, 
			      char *argv[]) 
{    
    cid_arg args;
    container *c;
    element *e;
    int i, j, num;

    cli_args a[] = {
	{"-container_id", ARG_INT, 1, NULL, offsetof(cid_arg, container_id)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    c = get_container(args.container_id);
    num = 0;
    for (i = 0; i < c->num_rows; i++) {
	for (j = 0; j < c->num_columns; j++) {
	    e = c->matrix[i][j];
	    if (e) {
		if (c->matrix[i][j]->type == CANVAS) {
		    num += e->num_results;
		} else if (c->matrix[i][j]->type == SEQ_EDITOR) {
		    num += e->num_results;
		}
	    }
	}
    }
    vTcl_SetResult(interp, "%d", num);
    return TCL_OK;
}

int tcl_get_new_container(ClientData clientData, 
			  Tcl_Interp *interp, 
			  int argc, 
			  char *argv[]) 
{    
    Tcl_Obj *klist;
    element_info e_info;
    char *c_win;
    int container_id;
    
    container_id = new_container(interp, &c_win);

    e_info.element_id = -1;
    e_info.element_win = NULL;
    e_info.container_id = container_id;
    e_info.container_win = strdup(c_win);

    klist = TclX_NewKeyedListObj();
    if (!klist)
	return TCL_ERROR;
    Tcl_IncrRefCount(klist);

    setli(e_info, element_id);
    setli(e_info, container_id);
    setls(e_info, element_win);
    setls(e_info, container_win);

    /* Set the result and return */
    Tcl_SetObjResult(interp, klist);
    Tcl_DecrRefCount(klist);

    return TCL_OK;
}

int tcl_get_new_element(ClientData clientData, 
			  Tcl_Interp *interp, 
			  int argc, 
			  char *argv[]) 
{    
    Tcl_Obj *klist;
    element_info e_info;
    char *e_win;
    int element_id;
    
    element_id = new_element(interp, &e_win);

    e_info.container_id = -1;
    e_info.container_win = NULL;
    e_info.element_id = element_id;
    e_info.element_win = strdup(e_win);

    klist = TclX_NewKeyedListObj();
    if (!klist)
	return TCL_ERROR;

    Tcl_IncrRefCount(klist);

    setli(e_info, element_id);
    setli(e_info, container_id);
    setls(e_info, element_win);
    setls(e_info, container_win);

    /* Set the result and return */
    Tcl_SetObjResult(interp, klist);
    Tcl_DecrRefCount(klist);

    xfree(e_win);
    return TCL_OK;
}

typedef struct {
    int element_id;
    double y;
} invert_arg;

int tcl_invert_wy(ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, 
		  char *argv[]) 
{    
    invert_arg args;
    element *e;
    double wy;

    cli_args a[] = {
	{"-element_id", ARG_INT, 1, NULL, offsetof(invert_arg, element_id)},
	{"-y",          ARG_DOUBLE, 1, NULL, offsetof(invert_arg, y)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    e = get_element(args.element_id);

    wy = invert_wy(e, args.y);

    vTcl_SetResult(interp, "%f", invert_wy(e, args.y));
    return TCL_OK;
}


/* return keyed list of element structures */
int tcl_element_struct(ClientData clientData, 
		     Tcl_Interp *interp, 
		     int argc, 
		     char *argv[]) 
{    
    element *e; 
    Tcl_Obj *klist;

    if (argc < 2) {
	goto e_error;
    }

    e = get_element(atoi(argv[1]));

    klist = TclX_NewKeyedListObj();
    if (!klist)
	return TCL_ERROR;

    Tcl_IncrRefCount(klist);

    setli(*e, container_id);
    setli(*e, id);
    setls(*e, win);
    setli(*e, type);
    setli(*e, orientation);
    setli(*e, crosshair);   
    setli(*e, num_results); 
    setld(*e, max_y);
    setld(*e, min_y);
    setli(*e, min_x);
    setli(*e, max_x);
    setli(*e, ruler_id);
    setli(*e, num_seqs);
    setli(*e, cursor_visible);
    setli(*e, status);
    setli(*e, seq_e_id);
    setli(*e, num_children_left);
    setli(*e, num_children_right);
    setli(*e, children_id);
    setli(*e, children_position);
    setli(*e, parent);

    /* Set the result and return */
    Tcl_SetObjResult(interp, klist);
    Tcl_DecrRefCount(klist);

    return TCL_OK;

e_error:
    Tcl_SetResult(interp,
		  "wrong # args: should be \"tcl_element_struct element_id\"\n",
		  TCL_STATIC);
    return TCL_ERROR;

}


int tcl_element_info(ClientData clientData, 
		     Tcl_Interp *interp, 
		     int argc, 
		     char *argv[]) 
{    
    element *e; 
    int i;
    Tcl_Obj *list = Tcl_NewListObj(0, NULL);
    Tcl_Obj *int_obj;

    if (argc < 3) {
	goto e_error;
    }

    e = get_element(atoi(argv[1]));
    if (strcmp(argv[2], "num_results") == 0 ) {
	vTcl_SetResult(interp, "%d", e->num_results);
    } else if (strcmp(argv[2], "seq_id_h") == 0 ) {
	
	if (!list)
	    return TCL_ERROR;
	
	Tcl_IncrRefCount(list);

	for (i = 0; i < e->num_seqs; i++) {
	    int_obj = Tcl_NewIntObj(e->seqs[i].seq_id);
	    Tcl_ListObjAppendElement(interp, list, int_obj);
	}

	Tcl_SetObjResult(interp, list);
	Tcl_DecrRefCount(list);
    }
    return TCL_OK;

e_error:
    Tcl_SetResult(interp,
		  "wrong # args: should be \"tcl_element_info element_id command\"\n",
		  TCL_STATIC);
    return TCL_ERROR;
}


void Container_Init(Tcl_Interp *interp)
{
    
    Tcl_CreateCommand(interp, "container_scroll_x", 
		      tcl_container_scroll_x, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "container_scroll_y", 
		      tcl_container_scroll_y, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);

    Tcl_CreateCommand(interp, "container_zoom", 
		      tcl_container_zoom, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "element_resize", 
		      tcl_element_resize, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "draw_container_crosshair", 
		      tcl_draw_container_crosshair, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "delete_container_crosshair", 
		      tcl_delete_container_crosshair, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);

    Tcl_CreateCommand(interp, "pixel_to_world", 
		      tcl_pixel_to_world, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);

    Tcl_CreateCommand(interp, "container_shutdown", 
		      tcl_container_shutdown, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);    

    Tcl_CreateCommand(interp, "num_container_elements", 
		      tcl_num_container_elements, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "num_container_results", 
		      tcl_num_container_results, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_new_container", 
		      tcl_get_new_container, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_new_element", 
		      tcl_get_new_element, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "invert_wy", 
		      tcl_invert_wy, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "element_info", 
		      tcl_element_info, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "element_struct", 
		      tcl_element_struct, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
}

