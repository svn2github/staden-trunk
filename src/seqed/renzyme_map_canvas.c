#include <tcl.h>
#include <tk.h>
#include <itcl.h>
#include <itk.h>

#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tcl_utils.h"
#include "misc.h"
#include "text_output.h"
#include "canvas_box.h"
#include "editor_reg.h"
#include "seq_results.h"
#include "graphic_editor.h"
#include "renzyme_map_canvas.h"
#include "editor.h"
#include "nip_structs.h"
/******* NEEDED FOR ********
   nip_scroll_arg;
   nip_zoom_arg;
   nip_cursor_arg;
   nip_enz_name_arg;
*/
#include "nip_results.h"
/******* NEEDED FOR *******
  out_canvas
*/


/*
 * ---------------------------------------------------------------------------
 * Cursor manipulation
 * ---------------------------------------------------------------------------
 */


extern void nip_canvas_cursor_delete(Tcl_Interp *interp,
				     cursor_e *cursor,
				     CanvasPtr *canvas, 
				     win **win_list, 
				     int num_wins);

extern int canvas_cursor_show(Tcl_Interp *interp, 
			      CanvasPtr *canvas, 
			      win **win_list, 
			      int num_wins,
			      WorldPtr *world, 
			      int cursorx, 
			      int sent_by, 
			      int reg_id);

/*
 * Moves and existsing cursor - calls canvas_cursor_move in tcl
 *
 * Returns whether we've scrolled the canvas.
 */
int ed_nip_canvas_cursor_move(Tcl_Interp *interp, 
			      int seq_id,
			      cursor_e *cursor, 
			      CanvasPtr *canvas, 
			      win **win_list, 
			      int num_wins, 
			      int reg_id, 
			      WorldPtr *world, 
			      int cursor_show)
{
    int i, apos, ret;
    char cmd[1024];
    double cx, cy;
    int seq_num;
    
    seq_num = GetEdenNum (seq_id);
    ret = 0;
    apos = cursor->abspos;
    if (apos < 1)
	apos = 1;
    if (apos > GetEdenLength (seq_num) + 1)
	apos = GetEdenLength (seq_num) + 1;
    
    for (i = 0; i < num_wins; i++) {
	/* only move cursors in the x direction */
	if (win_list[i]->scroll == 'x' || win_list[i]->scroll == 'b') {

	    WorldToCanvas(canvas, (double)apos, 0, &cx, &cy);
	    /*printf("seq_id=%d cursor->id=%d  reg_id=%d cursor->colour%s cx=%f\n",
	      seq_id, cursor->id,reg_id,cursor->colour, cx);*/
	    sprintf(cmd, "ed_nip_canvas_cursor_move %d %s %d %d %s %f",
		    seq_id, win_list[i]->window, cursor->id, reg_id, 
		    cursor->colour, cx);

	    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
		verror(ERR_WARN, "ed_nip_canvas_cursor_move", "%s\n", 
		       interp->result);
	    }
	}
    }

    /* Make sure the cursor is still visible */
    if (cursor_show) {
	return canvas_cursor_show(interp, canvas, win_list, num_wins,
				  world, apos, cursor->sent_by, reg_id);
    }
    return -1;
}

/*
 * Refreshes the cursor drawing. This involves moving, creating or deleting
 * ed_find_cursor * as the position and ref counts change.
 * Returns whether we've scrolled the canvas and updates the visibility
 * pointer.
 */

int ed_nip_canvas_cursor_refresh(Tcl_Interp *interp, 
				 int seq_num,
				 cursor_e *changed_c, 
				 cursor_e *canvas_c,
				 CanvasPtr *canvas, 
				 win **win_list, 
				 int num_wins,
				 int reg_id, 
				 int *visible, 
				 WorldPtr *world,
				 int cursor_show)
{
    int r;

    /* Check for deletion */
    if (changed_c->job & CURSOR_DELETE) {
	nip_canvas_cursor_delete(interp, changed_c, canvas, win_list, 
				 num_wins);
	*visible = 0;
	return 0;
    }
    
    /* If the cursor is our own with only one reference, don't draw it */
    if (changed_c == canvas_c && changed_c->refs <= 1) {
	if (*visible) {
	    /* In fact, if it's still visible, then remove it */
	    nip_canvas_cursor_delete(interp, changed_c, canvas, win_list, 
				     num_wins);
	    *visible = 0;
	}
	return 0;
    }
   
    /* Move it, creating it in the process if neccesary */
    r = ed_nip_canvas_cursor_move(interp, seq_num, changed_c,
			   canvas, win_list, num_wins, reg_id, world,
			   cursor_show);
    *visible = 1;

    return r;
}

int EdNipScrollCanvas(ClientData clientData, 
		      Tcl_Interp *interp, 
		      int argc, 
		      char *argv[]) 
{
    nip_scroll_arg args;
    seq_result *result;
    seq_reg_info info;
    ed_renz_res *data;

    cli_args a[] = {
	{"-id",             ARG_INT, 1, NULL, offsetof(nip_scroll_arg, id)},
	{"-xscrollcommand", ARG_STR, 1, "", offsetof(nip_scroll_arg, xscroll)},
	{"-yscrollcommand", ARG_STR, 1, "", offsetof(nip_scroll_arg, yscroll)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    editor_result_notify(args.id, (editor_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    data = result->data;

    if (strcmp(args.xscroll, "") != 0) {
	canvasScrollX(interp, data->re_win, data->win_list, 
		      data->num_wins,
		      data->world->visible, data->canvas, args.xscroll);
    }
    if (strcmp(args.yscroll, "") != 0) {
	canvasScrollY(interp, data->re_win, data->win_list, 
		      data->num_wins,
		      data->world->visible, data->canvas, args.yscroll);
    }    
    return TCL_OK;
}

int EdNipZoomCanvas(ClientData clientData, 
		    Tcl_Interp *interp, 
		    int argc, 
		    char *argv[]) 
{
    box *zoom;
    nip_zoom_arg args;
    seq_result *result;

    out_canvas_e *output;
    seq_reg_info info;
    ed_renz_res *data;

    cli_args a[] = {
	{"-seq_id",    ARG_INT, 1, NULL, offsetof(nip_zoom_arg, seq_id)},
	{"-id",        ARG_INT, 1, NULL, offsetof(nip_zoom_arg, id)},
	{"-x1",        ARG_INT, 1, "-1", offsetof(nip_zoom_arg, x1)},
	{"-y1",        ARG_INT, 1, "-1", offsetof(nip_zoom_arg, y1)},
	{"-x2",        ARG_INT, 1, "-1", offsetof(nip_zoom_arg, x2)},
	{"-y2",        ARG_INT, 1, "-1", offsetof(nip_zoom_arg, y2)},
	{"-direction", ARG_STR, 1, "b", offsetof(nip_zoom_arg, scroll)},
	{NULL,	0,	0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    editor_result_notify(args.id, (editor_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    data = result->data;
    output = result->output;

    /* zoom back */
    if (args.x1 == -1 && args.y1 == -1 && args.x2 == -1 && args.y2 == -1) {
	canvasZoomback (interp, data->canvas, data->re_win, data->world,
			data->win_list, data->num_wins, &data->zoom);
    } else {

	/* zoom */
	if (NULL == (zoom = (box *)xmalloc(sizeof(box))))
	    return TCL_OK;
	zoom->x1 = args.x1;
	zoom->y1 = args.y1;
	zoom->x2 = args.x2;
	zoom->y2 = args.y2;

	canvasZoom(interp, data->canvas, data->re_win, data->world,
		   data->win_list, data->num_wins, &data->zoom, zoom, 
		   args.scroll[0]);
	xfree(zoom);
    }

    /* redraw ruler ticks by redrawing ruler */
    
    /* draw_single_ruler(interp, data->ruler, data->canvas, 1); */
    draw_single_ruler(interp, data->ruler, data->canvas, data->ruler->start,
		      data->ruler->end, 1);
    scaleSingleCanvas(interp, data->world, data->canvas, data->ruler->window, 
		      'x', "all");
		      
    /* redraw cursor */
    ed_nip_canvas_cursor_refresh(interp, args.seq_id, output->cursor, 
				 output->cursor, data->canvas, data->win_list, 
				 data->num_wins, result->id, 
				 &output->cursor_visible, data->world, 1);

    return TCL_OK;
}

int EdNipCanvasCursorX(ClientData clientData, 
		     Tcl_Interp *interp, 
		     int argc, 
		     char *argv[]) 
{
    nip_cursor_arg args;
    seq_result *result;
    seq_reg_info info;
    ed_renz_res *data;
    double wx, wy;
    char *label;
    
    cli_args a[] = {
	{"-id", ARG_INT, 1, NULL, offsetof(nip_cursor_arg, id)},
	{"-x",  ARG_INT, 1, NULL, offsetof(nip_cursor_arg, cx)},
	{NULL,	0,	0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    editor_result_notify(args.id, (editor_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    data = result->data;

    CanvasToWorld(data->canvas, args.cx, 0, &wx, &wy);
    label = get_default_astring(interp, tk_utils_defs, w("R_ENZ.CURSOR"));

    canvasCursorX(interp, data->canvas, data->frame, label, 
		  data->cursor.colour, data->cursor.width, args.cx, wx, 
		  data->win_list, data->num_wins);
    
    xfree(label);
    return TCL_OK;
}

int 
EdNipGetREnzName(ClientData clientData,
		 Tcl_Interp *interp, 
		 int argc, 
		 char *argv[])
{
    nip_enz_name_arg args;
    seq_result *result;
    seq_reg_info info;
    ed_renz_res *data;

    cli_args a[] = {
	{"-id",	    ARG_INT, 1, NULL, offsetof(nip_enz_name_arg, id)},
	{"-enzyme", ARG_INT, 1, NULL, offsetof(nip_enz_name_arg, enzyme)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    editor_result_notify(args.id, (editor_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    data = result->data;

    vTcl_SetResult(interp, "%s", data->r_enzyme->renzyme[args.enzyme]->name);

    return TCL_OK;
}

int EdNipGetREnzInfo(ClientData clientData,
	       Tcl_Interp *interp, 
	       int argc, 
	       char *argv[])
{
    nip_enz_name_arg args;
    seq_reg_generic gen;

    cli_args a[] = {
	{"-id",	    ARG_INT, 1, NULL, offsetof(nip_enz_name_arg, id)},
	{"-enzyme", ARG_INT, 1, NULL, offsetof(nip_enz_name_arg, enzyme)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    gen.job = SEQ_GENERIC;
    gen.task = TASK_NIP_RENZ_INFO;
    gen.data = (void *)&args.enzyme;

    vfuncgroup(5, "restriction enzymes");
    editor_result_notify(args.id, (editor_reg_data *)&gen, 0); 

    return TCL_OK;
}

int EdResizeCanvas(ClientData clientData, 
		   Tcl_Interp *interp, 
		   int argc, 
		   char *argv[]) 
{
    nip_resize_arg args;
    seq_result *result;
    seq_reg_info info;
    ed_renz_res *data;
    
    cli_args a[] = {
	{"-id", ARG_INT, 1, NULL, offsetof(nip_resize_arg, id)},
	{NULL,	0,	0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    editor_result_notify(args.id, (editor_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    data = result->data;

    resizeCanvas(interp, data->re_win, data->win_list, data->num_wins,
		 data->world->visible, data->world->total, data->canvas);
    
    return TCL_OK;
}

int EdSLength(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    seq_id_arg args;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, "-1", offsetof(seq_id_arg, seq_id)},
	{NULL,      0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (args.seq_id == -1) {
	/*vTcl_SetResult(interp, "%d", GetSequenceLength(GetActiveSeqNumber(0)));*/
    } else {
	/*printf ("args.seq_id=%d\n", args.seq_id);*/
	vTcl_SetResult(interp, "%d", GetEdenLength(GetEdenNum(args.seq_id)));
    }
    return TCL_OK;
}

int EdCursorNotify(ClientData clientData, 
		   Tcl_Interp *interp, 
		   int argc, 
		   char *argv[]) 
{
    ed_cursor_notify_arg args;
    seq_reg_cursor_notify cn;
    int seq_id, member_id;

    cli_args a[] = {
	{"-seq_num", ARG_INT, 1, NULL, offsetof(ed_cursor_notify_arg, seq_num)},
	{"-id",      ARG_INT, 1, NULL, offsetof(ed_cursor_notify_arg, id)},
	{"-pos",     ARG_INT, 1, NULL, offsetof(ed_cursor_notify_arg, pos)},
	{"-direction", ARG_INT, 1, "0", offsetof(ed_cursor_notify_arg, direction)},
	{NULL,       0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (cn.cursor = ed_find_cursor(&args.seq_num, args.id, args.direction))) {
#ifdef DEBUG
	printf("CursorNotify: unable to find cursor for seq_num %d cursor_id %d\n",
	       args.seq_num, args.id);
#endif
	return TCL_OK;
    }
    seq_id = GetEdenId (args.seq_num);
    member_id = GetMemberIdFromSeqId (seq_id);
    cn.cursor->abspos = args.pos;
    cn.cursor->posy = member_id;
    cn.cursor->job = CURSOR_MOVE;
    cn.cursor->sent_by = -1; /* FIXME: Don't know correct value. */
    cn.job = SEQ_CURSOR_NOTIFY;

    editor_notify(args.seq_num, (editor_reg_data *)&cn);

    return TCL_OK;
}

int EdNipCanvasToWorld(ClientData clientData,
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[])
{
    nip_world_arg args;
    double wx, wy;
    ed_renz_res *r;
    seq_result *result;
    seq_reg_info info;

    cli_args a[] = {
	{"-id", ARG_INT, 1, NULL, offsetof(nip_world_arg, id)},
	{"-x",  ARG_INT, 1, NULL, offsetof(nip_world_arg, cx)},
	{NULL,	0,	 0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    editor_result_notify(args.id, (editor_reg_data *)&info, 0);
    if (!info.result) {
	return TCL_OK;
    }
    result = (seq_result *)info.result;
    r = result->data;
   
    CanvasToWorld(r->canvas, args.cx, 0, &wx, &wy);
    vTcl_SetResult(interp, "%d", (int)wx);
    return TCL_OK;
}

typedef struct {
    int id;
    char *option;
} update_arg;

/*
 * hides or deletes all results on a single raster widget
 */
int
EditorResultUpdate(ClientData clientData, 
		   Tcl_Interp *interp, 
		   int argc, 
		   char *argv[]) 
{
    update_arg args;
    seq_reg_info info;

    cli_args a[] = {
	{"-index", ARG_INT, 1, "-1", offsetof(update_arg, id)},
	{"-job",   ARG_STR, 1, NULL, offsetof(update_arg, option)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (strcmp(args.option, "HIDE") == 0) {
	info.job = SEQ_HIDE;
    } else if (strcmp(args.option, "REVEAL") == 0) {
	info.job = SEQ_REVEAL;
    } else if (strcmp(args.option, "DELETE") == 0) {
	info.job = SEQ_DELETE;
    } else if (strcmp(args.option, "QUIT") == 0) {
	info.job = GRAPHIC_EDITOR_QUIT;
    } else {
	verror(ERR_FATAL, "seq_result_notify_all", "invalid command");
	return TCL_OK;
    }

    editor_result_notify(args.id, (editor_reg_data *)&info, 1);
    /*if (args.id == -1) {
	seq_result_notify_all((editor_reg_data *)&info);
    } else {
	editor_result_notify(args.id, (editor_reg_data *)&info, 1);
	}*/
    return TCL_OK;
}

int EdQueryCursor(ClientData clientData, Tcl_Interp *interp,
		  int argc, char **argv)
{
    eqc_arg args;

    cli_args a[] = {
	{"-cursorid", ARG_INT, 1, NULL, offsetof(eqc_arg, cursorid)},
	{"-seq_num",  ARG_INT, 1, "0",  offsetof(eqc_arg, seq_num)},
        {NULL,        0,       0, NULL, 0}
    };
    cursor_e *gc;

    if (-1 == parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    if (NULL == (gc = ed_find_cursor(&args.seq_num, args.cursorid, -1)))
	return TCL_OK;

    vTcl_SetResult(interp,
		   "{id %d} {refs %d} {private %d} {abspos %d}",
		   gc->id, gc->refs, gc->private, gc->abspos);
    return TCL_OK;
}


int RenzCanvas_Init(Tcl_Interp *interp) {
 
    Tcl_CreateCommand(interp, "ed_nip_scroll_canvas", EdNipScrollCanvas, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "ed_nip_zoom_canvas", EdNipZoomCanvas, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "ed_nip_canvas_cursor_x", EdNipCanvasCursorX, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "ed_nip_get_renz_name", EdNipGetREnzName, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "ed_nip_get_renz_info", EdNipGetREnzInfo, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "ed_nip_resize_canvas", EdResizeCanvas, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "ed_s_length", EdSLength, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL); 
    Tcl_CreateCommand(interp, "ed_cursor_notify", EdCursorNotify, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "ed_nip_canvas_to_world", EdNipCanvasToWorld, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "editor_result_update", EditorResultUpdate, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "ed_query_cursor", EdQueryCursor,
		      (ClientData) NULL,
		      NULL);
    
    return TCL_OK;
}

