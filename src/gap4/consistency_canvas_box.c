#include <tcl.h>

#include "io-reg.h"
#include "gap_canvas_box.h"
#include "consistency_display.h"

/*
 * Updates the scroll region to ensure that the cursor is visible.
 * The absolute cursor position (in bases) on the display is 'cursorx'.
 *
 * Returns whether we've scrolled the canvas.
 */
int consistency_cursor_show(Tcl_Interp *interp, GapIO *io,
			    obj_consistency_disp *c,
			    CanvasPtr *canvas, win **win_list, int num_wins,
			    WorldPtr *world, int cursorx, int sent_by, 
			    int reg_id)
{
    int x1, x2, dx;
    char cmd[1024];
    double fract;

#ifdef DEBUG
    printf("consistency_cursor_show %d %f %f\n", cursorx, world->visible->x1,
	   world->visible->x2);
#endif

    /* Is it visible? */
    if (cursorx >= world->visible->x1 && cursorx <= world->visible->x2) {
#ifdef DEBUG
        printf("VISIBLE 0\n");
#endif
	return 0;
    }
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
    sprintf(cmd, "moveto %.20f", fract);
    consistency_canvasScrollX(interp, c, win_list, num_wins, cmd);

#ifdef DEBUG
        printf("VISIBLE 1\n");
#endif
    return 1;
}

/*
 * Moves and existsing cursor - calls canvas_cursor_move in tcl
 *
 * Returns whether we've scrolled the canvas.
 */
int consistency_cursor_move(Tcl_Interp *interp, 
			    GapIO *io, 
			    obj_consistency_disp *c,
 			    int cnum,
			    cursor_t *cursor, 
			    CanvasPtr *canvas, 
			    win **win_list,
			    int num_wins, 
			    int reg_id,
			    int offset, 
			    WorldPtr *world,
			    int cursor_show)
{
    int i, apos, ret;
    double cx, cy;
    char cmd[1024];

#ifdef DEBUG
    printf("canvas_cursor_move %d\n", cursor->abspos);
#endif
    ret = 0;
    apos = cursor->abspos;
    if (apos < 1)
	apos = 1;
    if (apos > io_clength(io, cnum) + 1)
	apos = io_clength(io, cnum) + 1;
    for (i = 0; i < num_wins; i++) {
	/* only move cursors in the x direction */
	if (win_list[i]->scroll == 'x' || win_list[i]->scroll == 'b') {

	    WorldToCanvas(canvas, (double)(apos + offset), 0,
			  &cx, &cy);
	    sprintf(cmd, "canvas_cursor_move %d %d %s %d %d %.20f",
		    *handle_io(io), cnum, win_list[i]->window,
		    cursor->id, reg_id, cx);
#ifdef DEBUG
	    printf("%s\n", cmd);
#endif

	    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
		printf("consistency_cursor_move: %s\n", Tcl_GetStringResult(interp));
	    }
	}
    }

    /* Make sure the cursor is still visible */
    if (cursor_show) 
	return consistency_cursor_show(interp, io, c, canvas, win_list, 
				       num_wins, world, apos+offset,
				       cursor->sent_by, reg_id);

    return 0;
}

/*
 * Refreshes the cursor drawing. This involves moving, creating or deleting
 * as the position and ref counts change.
 * Returns whether we've scrolled the canvas and updates the visibility
 * pointer.
 */
int consistency_cursor_refresh(Tcl_Interp *interp, GapIO *io, 
			       obj_consistency_disp *c, int cnum,
			       cursor_t *changed_c, cursor_t *canvas_c,
			       CanvasPtr *canvas, win **win_list, int num_wins,
			       int reg_id, int offset, int *visible,
			       WorldPtr *world, int cursor_show)
{
    int r;

#ifdef DEBUG
    printf("canvas_cursor_refresh offset %d changed job %d\n", 
	   offset, changed_c->job);
#endif

    /* Check for deletion */
    if (changed_c->job & CURSOR_DELETE) {
	canvas_cursor_delete(interp, io, changed_c,
			     canvas, win_list, num_wins);
	if (canvas_c == changed_c)
	    *visible = 0;
	return 0;
    }
    
    /* If the cursor is our own with only one reference, don't draw it */
    if (changed_c == canvas_c && changed_c->refs <= 1) {
	if (*visible) {
	    /* In fact, if it's still visible, then remove it */
	    canvas_cursor_delete(interp, io, changed_c,
				 canvas, win_list, num_wins);
	    *visible = 0;
	}
	return 0;
    }
    
    /* Move it, creating it in the process if neccesary */
    r = consistency_cursor_move(interp, io, c, cnum, changed_c,
				canvas, win_list, num_wins, reg_id, offset, world,
				cursor_show);
    *visible = 1;

    return r;
}
