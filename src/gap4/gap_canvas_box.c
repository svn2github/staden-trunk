#include <tcl.h>

#include "io-reg.h"
#include "gap_globals.h"
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
int canvas_cursor_show(Tcl_Interp *interp, GapIO *io,
		       CanvasPtr *canvas, win **win_list, int num_wins,
		       WorldPtr *world, int cursorx, int sent_by, int reg_id)
{
    int x1, x2, dx;
    char cmd[1024];
    double fract;

#ifdef DEBUG
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!canvas_cursor_show %d %f %f reg_id %d sent_by %d\n", cursorx, world->visible->x1,
	       world->visible->x2, reg_id, sent_by);
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
    canvasScrollX(interp, win_list[0]->window, win_list, num_wins,
		  world->visible, canvas, cmd);

#ifdef DEBUG
        printf("VISIBLE 1\n");
#endif
    return 1;
}

/*
 * Deletes an existing cursor - calls canvas_cursor_delete in tcl
 */
void canvas_cursor_delete(Tcl_Interp *interp, GapIO *io, cursor_t *cursor,
			  CanvasPtr *canvas, win **win_list, int num_wins)
{
    int i;
    char cmd[1024];

#ifdef DEBUG
    printf("canvas_cursor_delete\n");
#endif
    for (i = 0; i < num_wins; i++) {
	/* existing cursors will only be in the x direction */
	if (win_list[i]->scroll == 'x' || win_list[i]->scroll == 'b') {

	    sprintf(cmd, "canvas_cursor_delete %d %s %d",
		    *handle_io(io), win_list[i]->window, cursor->id);

	    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
		printf("canvas_cursor_delete: %s\n", Tcl_GetStringResult(interp));
	    }

	}
    }
}

/*
 * Moves and existsing cursor - calls canvas_cursor_move in tcl
 *
 * Returns whether we've scrolled the canvas.
 */
int canvas_cursor_move(Tcl_Interp *interp,
		       GapIO *io,
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

	    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
		printf("canvas_cursor_move: %s\n", Tcl_GetStringResult(interp));
	    }
	}
    }

    /* Make sure the cursor is still visible */
    if (cursor_show)
	return canvas_cursor_show(interp, io, canvas, win_list, num_wins,
				  world, apos+offset, cursor->sent_by,
				  reg_id);

    return 0;
}

/*
 * Refreshes the cursor drawing. This involves moving, creating or deleting
 * as the position and ref counts change.
 * Returns whether we've scrolled the canvas and updates the visibility
 * pointer.
 */
int canvas_cursor_refresh(Tcl_Interp *interp, GapIO *io, int cnum,
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
    r = canvas_cursor_move(interp, io, cnum, changed_c,
			   canvas, win_list, num_wins, reg_id, offset, world,
			   cursor_show);
    *visible = 1;

    return r;
}


/*
 * find the contig the cursor is in
 */
int find_cursor_contig(GapIO *io,
		       int id,
		       c_offset *contig_offset,
		       int *contig_array,
		       int num_contigs,
		       double wx)
{
    int i;
    int offset = 0;
    int prev_offset = 0;
    int max_x = 0;
    int last_contig;

    /*
     * a couple of fudges: if num_contigs is 1 then wx is still wx
     * if wx < 0 then must be to the left of the 1st contig and wx is still wx
     */
    if ((num_contigs == 1) || (wx < 0)) {
	return contig_array[0];
    }
    max_x = io_clength(io, contig_array[0]);
    last_contig = contig_array[0];

    for (i = 1; i < num_contigs; i++) {
	prev_offset = offset;

	offset = contig_offset[contig_array[i]].offset;

	if (max_x < offset + io_clength(io, contig_array[i])) {
	    max_x = offset + io_clength(io, contig_array[i]);
	    last_contig = contig_array[i];
	}
	if ((wx > prev_offset) && (wx <= offset)) {
	    return (contig_array[i-1]);
	}
    }

    if (wx < offset + io_clength(io, contig_array[num_contigs-1])) {
	return (contig_array[num_contigs-1]);
    } else {
	return (last_contig);
    }
}
