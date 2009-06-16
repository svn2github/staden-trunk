#include "seq_reg.h"
#include "seq_results.h"
#include "nip_results.h"
#include "text_output.h"
#include "canvas_box.h"

/*
 * ---------------------------------------------------------------------------
 * Cursor manipulation
 * ---------------------------------------------------------------------------
 */

/*
 * Updates the scroll region to ensure that the cursor is visible.
 * The absolute cursor position (in bases) on the display is 'cursorx'.
 *
 * Returns whether we've scrolled the canvas.
 */
int canvas_cursor_show(Tcl_Interp *interp, 
		       CanvasPtr *canvas, win **win_list, int num_wins,
		       WorldPtr *world, int cursorx, int sent_by, int reg_id)
{
    int x1, x2, dx;
    char cmd[1024];
    double fract;

    /* Is it visible? */
    if (cursorx >= world->visible->x1 && cursorx <= world->visible->x2)
	return 0;

    /* Get new x1 and x2, centred around cursorx */
    dx = world->visible->x2 - world->visible->x1;

    if (reg_id != sent_by) {
	x1 =  cursorx - dx / 2;
	if (x1 < world->total->x1)
	    x1 = world->total->x1;
	if (x1 > world->total->x2 - dx)
	    x1 = world->total->x2 - dx;
	x2 = x1 + dx;
    } else {
	/* scrolling to the right */
	if (cursorx > world->visible->x1) {
	    x2 = cursorx;
	    if (x2 > world->total->x2)
		x2 = world->total->x2;

	    if (x2 < world->total->x1 + dx)
		x2 =  world->total->x1 + dx;
	    
	    x1 = x2 - dx;
	    
	} else {
	    /* scrolling to the left */
	    x1 = cursorx; 
	    
	    if (x1 < world->total->x1)
	    x1 = world->total->x1;
	    
	    if (x1 > world->total->x2 - dx)
		x1 = world->total->x2 - dx;
	    
	    x2 = x1 + dx;
	}
    }

    /* Compute xview fraction and scroll */
    fract = (double)(x1 - world->total->x1) /
	(world->total->x2 - world->total->x1);
    sprintf(cmd, "moveto %f", fract);
    canvasScrollX(interp, win_list[0]->window, win_list, num_wins,
		  world->visible, canvas, cmd);

    return 1;
}

/*
 * Deletes an existing cursor - calls canvas_cursor_delete in tcl
 */
void nip_canvas_cursor_delete(Tcl_Interp *interp,cursor_t *cursor,
			      CanvasPtr *canvas, win **win_list, int num_wins)
{
    int i;
    char cmd[1024];

    for (i = 0; i < num_wins; i++) {
	/* existing cursors will only be in the x direction */
	if (win_list[i]->scroll == 'x' || win_list[i]->scroll == 'b') {

	    sprintf(cmd, "nip_canvas_cursor_delete %s %d",
		    win_list[i]->window, cursor->id);

	    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
		verror(ERR_WARN, "canvas_cursor_delete", "%s\n", 
		       Tcl_GetStringResult(interp));
	    }
	    
	}
    }
}


/*
 * Moves and existsing cursor - calls canvas_cursor_move in tcl
 *
 * Returns whether we've scrolled the canvas.
 */
int nip_canvas_cursor_move(Tcl_Interp *interp, int seq_id,
			   cursor_t *cursor, CanvasPtr *canvas, 
			   win **win_list, int num_wins, int reg_id, 
			   WorldPtr *world, int cursor_show)
{
    int i, apos, ret;
    char cmd[1024];
    double cx, cy;

    ret = 0;
    apos = cursor->abspos;
    if (apos < 1)
	apos = 1;
    if (apos > GetSeqLength(GetSeqNum(seq_id)) + 1)
	apos =  GetSeqLength(GetSeqNum(seq_id)) + 1;
    for (i = 0; i < num_wins; i++) {
	/* only move cursors in the x direction */
	if (win_list[i]->scroll == 'x' || win_list[i]->scroll == 'b') {

	    WorldToCanvas(canvas, (double)apos, 0, &cx, &cy);
	    sprintf(cmd, "nip_canvas_cursor_move %d %s %d %d %s %f",
		    seq_id, win_list[i]->window, cursor->id, reg_id, 
		    cursor->colour, cx);
	    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
		verror(ERR_WARN, "nip_canvas_cursor_move", "%s\n", 
		       Tcl_GetStringResult(interp));
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
 * as the position and ref counts change.
 * Returns whether we've scrolled the canvas and updates the visibility
 * pointer.
 */
int nip_canvas_cursor_refresh(Tcl_Interp *interp, 
			      int seq_num,
			      cursor_t *changed_c, 
			      cursor_t *canvas_c,
			      CanvasPtr *canvas, win **win_list, int num_wins,
			      int reg_id, int *visible, WorldPtr *world,
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
    r = nip_canvas_cursor_move(interp, seq_num, changed_c,
			   canvas, win_list, num_wins, reg_id, world,
			   cursor_show);
    *visible = 1;

    return r;
}
