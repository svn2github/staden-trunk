/*
 * functions specific to the spin registration scheme
 */
#include <string.h>
#include <limits.h>
#include <tcl.h>

#include "tcl_utils.h"
#include "container.h"
#include "text_output.h"
#include "xalloc.h"
#include "seq_reg.h"
#include "seq_results.h"
#include "seq_element.h"
#include "tclCanvGraph.h"
#include "element_canvas.h"
#include "container.h"

void seq_element_callback(int seq_num, void *obj, seq_reg_data *jdata);
int seq_canvas_cursor_refresh(Tcl_Interp *interp, int seq_id,
			      cursor_t *changed_c, cursor_t *canvas_c,
			      CanvasPtr *canvas, element *e, int reg_id,
			      int *visible, WorldPtr *world, int show_all);

void seq_element_start_shutdown(element *e, int job);

int comparison_type(void *v_result,
		    int type)
{
    if (type == SEQ_PLOT_PERM) {
	return 1;
    }
    return 0;
}

int seq_get_seq_type(int seq_id,
		     int type,
		     int frame,
		     element_info *e_info)
{
   int num_elements;
    seq_result **s_result;
    int num_funcs;
    int i;
    char *win;

    num_elements = seq_num_results();

    if (num_elements == 0) {
        return -1;
    }
    if (NULL == (s_result = (seq_result **)xmalloc(num_elements *
						   sizeof(seq_result *))))
        return -1;

     /* find all plotted results */
    search_reg_data(comparison_type, (void **)s_result, &num_funcs);
    for (i = num_funcs-1; i >= 0; i--) {
#ifdef DEBUG
        printf("seqid %d %d type %d %d frame %d %d\n",
               s_result[i]->seq_id[HORIZONTAL], seq_id, s_result[i]->type, type,
               s_result[i]->frame, frame);
#endif

        if (s_result[i]->seq_id[HORIZONTAL] == seq_id && s_result[i]->type & type) {
            if (/*s_result[i]->frame == 0 ||*/ frame == s_result[i]->frame) {
		if (s_result[i]->e) {
		    e_info->element_id = s_result[i]->e->id;
		    if (win = strrchr(s_result[i]->e->win, '.')) {
			e_info->element_win = strdup(win);
		    }
		    e_info->container_id = s_result[i]->e->c->id;
		    e_info->container_win = strdup(s_result[i]->e->c->win);
		    e_info->orientation = s_result[i]->e->orientation;
		    return 0;
		}
	    }
	}
    }
    xfree(s_result);
    return -1;
}

int seq_find_column(Tcl_Interp *interp,
		    int seq_id)
{
    int num_results;
    seq_result **s_result;
    int num_funcs;
    int column = -1;
    int i;

    num_results = seq_num_results();

    if (num_results == 0) {
        return -1;
    }
    if (NULL == (s_result = (seq_result **)xmalloc(num_results *
						   sizeof(seq_result *))))
        return -1;
    /* find all plotted results */
    search_reg_data(comparison_type, (void **)s_result, &num_funcs);

    for (i = 0; i < num_funcs; i++) {
	if (s_result[i]->e && seq_id == s_result[i]->seq_id[HORIZONTAL]) {
	    column = get_element_column(interp, s_result[i]->e->win);
	    return column;
	} else {
	    column = get_default_int(interp, tk_utils_defs, w("CONTAINER.RULER_COL"));
	    return column;

	}
    }
    column = get_default_int(interp, tk_utils_defs, w("CONTAINER.RULER_COL"));
    return column;

}

void seq_update_cursor(Tcl_Interp *interp,
		       cursor_t *cursor,
		       int seq_id,
		       element *e,
		       int direction,
		       int show_all)
{
    cursor_t *e_cursor;
#ifdef DEBUG
    printf("seq_update_cursor\n");
#endif


    e_cursor = find_element_cursor(e, seq_id, direction);

    if (e_cursor)
	seq_canvas_cursor_refresh(interp, seq_id, cursor, e_cursor,
				  e->pixel, e, e->id, &e->cursor_visible,
				  e->world, show_all);
}

/*
 * Updates the scroll region to ensure that the cursor is visible.
 * The absolute cursor position (in bases) on the display is 'cursor_pos'.
 *
 * Returns whether we've scrolled the canvas.
 */
int seq_canvas_cursor_show(Tcl_Interp *interp,
			   CanvasPtr *canvas,
			   element *e,
			   WorldPtr *world,
			   int cursor_pos,
			   int direction,
			   int show_all)
{
    int x1, x2, dx, y1, y2, dy;
    char cmd[1024];
    double fract;
    coord *column = e->c->column[e->column_index];
    coord *row = e->c->row[e->row_index];

#ifdef DEBUG
    printf("SHOW %s %d %f %f %d\n", e->win, cursor_pos, column->visible.min, column->visible.max, show_all);
#endif
    /* when add new element, need to update all elements */
    if (!show_all) {
	/* Is it visible? */
	if (direction == HORIZONTAL) {
	    if (cursor_pos >= column->visible.min && 
		cursor_pos <= column->visible.max) {
		return 0;
	    }
	} else {
	    if (cursor_pos >= row->visible.min && 
		cursor_pos <= row->visible.max) {
		return 0;
	    }
	}
    }
    
    if (direction == HORIZONTAL) {
	/* Get new x1 and x2, centred around cursor_pos */
	dx = column->visible.max - column->visible.min;
	x1 =  cursor_pos - dx / 2;
	if (x1 < column->total.min)
	    x1 = column->total.min;
	if (x1 > column->total.max - dx)
	    x1 = column->total.max - dx;
	x2 = x1 + dx;

	/* Compute xview fraction and scroll */
	fract = (double)(x1 - column->total.min) /
	    (column->total.max - column->total.min);

	sprintf(cmd, "moveto %f", fract);
	container_scroll_x(interp, e->c->id,
			   get_element_column(interp, e->win), cmd);

    } else {
	/* Get new y1 and y2, centred around cursor_pos */
	dy = row->visible.max - row->visible.min;
	y1 =  cursor_pos - dy / 2;
	if (y1 < row->total.min)
	    y1 = row->total.min;
	if (y1 > row->total.max - dy)
	    y1 = row->total.max - dy;
	y2 = y1 + dy;

	/* Compute xview fraction and scroll */
	fract = (double)(y1 - row->total.min) /
	    (row->total.max - row->total.min);

	sprintf(cmd, "moveto %f", fract);

	/* e->scroll_y_func(interp, e, cmd); */
	container_scroll_y(interp, e->c->id,
			   get_element_row(interp, e->win), cmd);


    }
    return 1;
}

/*
 * Deletes an existing cursor - calls canvas_cursor_delete in tcl
 */
void seq_canvas_cursor_delete(Tcl_Interp *interp,
			      cursor_t *cursor,
			      element *e)
{
    int i, j;
    char cmd[1024];
    element *next;

#ifdef DEBUG
    printf("seq_canvas_cursor_delete\n");
#endif

    for (i = 0; i < e->c->num_rows; i++) {
	for (j = 0; j < e->c->num_columns; j++) {
	    next = e->c->matrix[i][j];
	    if (next != NULL) {
		sprintf(cmd, "%s delete cursor_%d", next->win, cursor->id);
		Tcl_Eval(interp, cmd);
	    }
	}
    }
}

/*
 * Moves and existing cursor - calls canvas_cursor_move in tcl
 *
 */
void seq_canvas_cursor_move(Tcl_Interp *interp,
			    int seq_id,
			    cursor_t *cursor,
			    CanvasPtr *pixel,
			    element *e,
			    int reg_id,
			    WorldPtr *world,
			    int show_all)
{
    int cx, cy, apos, ret;
    char cmd[1024];
    coord *column = e->c->column[e->column_index];
    coord *row = e->c->row[e->row_index];

    ret = 0;
    apos = cursor->abspos;

#ifdef DEBUG
    printf("seq_canvas_cursor_move\n");
#endif


    if (cursor->direction == HORIZONTAL) {
	if (!e->scroll_x_func)
	    return;

	world_to_pixel(column->pixel, (double)apos, (double)apos, &cx, &cy);
	sprintf(cmd, "seq_canvas_cursor_move_x %d %s %d %d %s %d %d",
		seq_id, e->win, cursor->id, reg_id,
		cursor->colour, cx, apos);
#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	if (TCL_ERROR == Tcl_Eval(interp, cmd)) {

	    printf("&&&&&&&&&&&&&&&&& %s\n", interp->result);
	    verror(ERR_WARN, "seq_canvas_cursor_move_x", "%s\n",
		   interp->result);
	}
	/* Make sure the cursor is still visible */
	seq_canvas_cursor_show(interp, pixel, e, world, apos,
			       cursor->direction, show_all);
    } else {
	if (!e->scroll_y_func)
	    return;

	world_to_pixel(row->pixel, (double)apos, (double)apos, &cx, &cy);
	cy = invert_cy(e, cy);

	sprintf(cmd, "seq_canvas_cursor_move_y %d %s %d %d %s %d %d",
		seq_id, e->win, cursor->id, reg_id,
		cursor->colour, cy, apos);
#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	    verror(ERR_WARN, "seq_canvas_cursor_move_y", "%s\n",
		   interp->result);
	}
	seq_canvas_cursor_show(interp, pixel, e, world, apos,
			       cursor->direction, show_all);
    }

    /* FIXME dont know why I have these loops, next is never used */
#if 0
    for (i = 0; i < e->c->num_rows; i++) {
	for (j = 0; j < e->c->num_columns; j++) {
	    next = e->c->matrix[i][j];
	    if (cursor->direction == HORIZONTAL) {
		world_to_pixel(column->pixel, (double)apos, (double)apos, &cx, &cy);
		sprintf(cmd, "seq_canvas_cursor_move_x %d %s %d %d %s %d %d",
			seq_id, e->win, cursor->id, reg_id,
			cursor->colour, cx, apos);
#ifdef DEBUG
		printf("%s\n", cmd);
#endif
		if (TCL_ERROR == Tcl_Eval(interp, cmd)) {

		    printf("&&&&&&&&&&&&&&&&& %s\n", interp->result);
		    verror(ERR_WARN, "seq_canvas_cursor_move_x", "%s\n",
			   interp->result);
		}
		/* Make sure the cursor is still visible */
		seq_canvas_cursor_show(interp, pixel, e, world, apos,
				       cursor->direction, show_all);
	    } else {
		world_to_pixel(row->pixel, (double)apos, (double)apos, &cx, &cy);
		cy = invert_cy(e, cy);

		sprintf(cmd, "seq_canvas_cursor_move_y %d %s %d %d %s %d %d",
			seq_id, e->win, cursor->id, reg_id,
			cursor->colour, cy, apos);
#ifdef DEBUG
		printf("%s\n", cmd);
#endif
		if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
		    verror(ERR_WARN, "seq_canvas_cursor_move_y", "%s\n",
			   interp->result);
		}
		seq_canvas_cursor_show(interp, pixel, e, world, apos,
				       cursor->direction, show_all);
	    }
	}
    }
#endif
}

/*
 * Refreshes the cursor drawing. This involves moving, creating or deleting
 * as the position and ref counts change.
 * Returns whether we've scrolled the canvas and updates the visibility
 * pointer.
 */
int seq_canvas_cursor_refresh(Tcl_Interp *interp, 
			      int seq_id,
			      cursor_t *changed_c, 
			      cursor_t *canvas_c,
			      CanvasPtr *canvas,
			      element *e,
			      int reg_id, 
			      int *visible, 
			      WorldPtr *world,
			      int show_all)
{
#ifdef DEBUG
    printf("seq_canvas_cursor_refresh\n");
#endif

    /* Check for deletion */
    if (changed_c->job & CURSOR_DELETE) {
	seq_canvas_cursor_delete(interp, changed_c, e);
	*visible = 0;
	return 0;
    }
    
    /* If the cursor is our own with only one reference, don't draw it */
    if (changed_c == canvas_c && changed_c->refs <= 1) {
	if (*visible) {
	    /* In fact, if it's still visible, then remove it */
	    seq_canvas_cursor_delete(interp, changed_c, e);
	    *visible = 0;
	}
	return 0;
    }
    
    /* Move it, creating it in the process if neccesary */
    seq_canvas_cursor_move(interp, seq_id, changed_c,
			   canvas, e, reg_id, world, show_all);

    *visible = 1;

    return 0;
}

void seq_element_shutdown(seq_element *result,
			  element *e)
{
    int i;
    int num;
    seq_reg_quit info;

#ifdef DEBUG
    printf("seq_element_shutdown %s %d\n", e->win, e->num_results);
    for (i = 0; i < e->num_seqs; i++) {
	printf("SEQS %d\n", e->seqs[i].seq_id);
    }
#endif

    /* remove all plots in element */
    info.job = SEQ_QUIT;
    for (i = 0; i < e->num_results; i++) {
	seq_result_notify(e->results[i]->result_id, (seq_reg_data *)&info, 0);
    }

    /* need to deregister from all seq_num the element is registered with */
    for (i = 0; i < e->num_seqs; i++) {
	num = GetSeqNum(e->seqs[i].seq_id);
	seq_deregister(num, seq_element_callback, result);
	delete_cursor(num, e->cursor[i]->id, 0);
    }
    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
}

void seq_element_callback(int seq_num, void *obj, seq_reg_data *jdata) 
{
    seq_element *result = (seq_element *) obj;
    element *e = result->e;

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "element plot");
	    return;
	}
    case SEQ_QUIT:
	{
	    /* quit element display */
	    /* 
	     * HACK to prevent re-entry into this routine when deleting a 
	     * result in the element which may invoke a QUIT event if the 
	     * result is the last one in the element
	     */
#ifdef DEBUG
	    printf("element SEQ_QUIT %d\n", e->status);
#endif
	    if (e->status) {
		return;
	    } else {
		e->status = 1;
	    }
	    seq_element_shutdown(result, e);
	    return;
	}
    case SEQ_GET_OPS:
	jdata->get_ops.ops = "Remove\0"; 
	return;

    case SEQ_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0: /* remove */
	    /* 
	     * HACK to prevent re-entry into this routine when deleting a 
	     * result in the element which may invoke a QUIT event if the 
	     * result is the last one in the element
	     */
#if 0
	    if (e->status) {
		return;
	    } else {
		e->status = 1;
	    }
#endif
	    seq_element_start_shutdown(e, ALL);
	    
	    return;
	}
	break;
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case WIN_NAME:
	    {
		char *r_win =  e->win;
		jdata->info.result = (void *)r_win;
		break;
	    }
	} 
	break;
    case SEQ_CURSOR_NOTIFY:
	{
	    cursor_t *cursor = (cursor_t *)jdata->cursor_notify.cursor;
	    cursor_t *gc;
	    int show_all = (int)jdata->cursor_notify.show_all;
	    int i;
	    int direction = -1;

	    for (i = 0; i < e->num_seqs; i++) {
		for (gc = e->cursor[i]; gc; gc = gc->next) {
		    if (GetSeqId(seq_num) == e->seqs[i].seq_id &&
			gc == cursor && 
			gc->direction == e->seqs[i].direction) {
			direction = cursor->direction;
			
			if (direction != -1) {
			    seq_update_cursor(e->c->interp, cursor, 
					      GetSeqId(seq_num), 
					      e, direction, show_all);
			} else {
			    printf("ERROR: SEQ_CURSOR_NOTIFY no direction found\n");
			}
		    } 
		}
	    }
	    e->cursor_array[cursor->id].prev_pos = cursor->abspos;
	    break;
	}
    }
}

/* need to register with sequences for the cursor stuff */
int seq_element_reg(Tcl_Interp *interp,
		    element *e)
{
    int i;
    seq_element *seq_e;
    int line_width;
    int *seq_counts_h, *seq_counts_v, nseq;


    if (NULL == (seq_e = (seq_element *)xmalloc(sizeof(seq_element)))) {
	    return -1;
    }

    /* FIXME: remove global variables */
    if (NULL == (e->cursor = (cursor_t **)xmalloc(MAX_NUM_SEQ * 
						  sizeof(cursor_t*))))
	return -1;

    if (NULL == (e->cursor_array = (cursor_info *)xmalloc(NUM_CURSORS * 
							  sizeof(cursor_info))))
	return -1;
    

    seq_e->op_func = seq_element_callback;

    /* initialised raster cursor array */
    for (i = 0; i < NUM_CURSORS; i++) {
	e->cursor_array[i].env = -2;
	e->cursor_array[i].visible[HORIZONTAL] = 0;
	e->cursor_array[i].visible[VERTICAL] = 0;
	e->cursor_array[i].prev_pos = -1;
    }

    line_width = get_default_int(interp, tk_utils_defs, w("CURSOR.LINE_WIDTH"));

    /* need to register with each seq_id in an element */
    nseq = NumSequences();
    if (NULL == (seq_counts_h = (int *)xmalloc(nseq * sizeof(int))))
	return -1;
    if (NULL == (seq_counts_v = (int *)xmalloc(nseq * sizeof(int))))
	return -1;
    for (i = 0; i < nseq; i++) {
	seq_counts_h[i] = 0;
	seq_counts_v[i] = 0;
    }

    for (i = 0; i < e->num_seqs; i++) {
	int snum;

	snum = GetSeqNum(e->seqs[i].seq_id);
	if (e->seqs[i].direction == HORIZONTAL) {
	    e->cursor[i] = 
		create_cursor(snum, 0, NULL, line_width, 
			      ++seq_counts_h[snum],
			      e->seqs[i].direction, 1);
	} else {
	    e->cursor[i] = 
		create_cursor(snum, 0, NULL, line_width, 
			      ++seq_counts_v[snum],
			      e->seqs[i].direction, 1);
	}
	e->cursor_array[e->cursor[i]->id].env = -1;
    }

    e->seq_e_id = get_reg_id();

    seq_e->e = e;

    for (i = 0; i < e->num_seqs; i++) {
	seq_register(GetSeqNum(e->seqs[i].seq_id), 
		     seq_element_callback, 
		     (void *)seq_e, SEQ_ELEMENT, e->seq_e_id);
    }    
    xfree(seq_counts_h);
    xfree(seq_counts_v);
    return 0;
}

/* 
 * return the element_info for either an existing element or create a new one
 */
element_info get_element_info(Tcl_Interp *interp,
			      seq_id_dir *seq_ids,
			      int num_seqs,
			      int plot_type,
			      int frame,
			      int create_new_container)
{
    int container_id, element_id;
    element_info e_info;
    char *e_win, *c_win, *e_win_to;
    int orientation = HORIZONTAL;

    e_win_to = NULL;
    e_info.element_id = -1;
    e_info.element_win = NULL;
    e_info.container_id = -1;
    e_info.container_win = NULL;
    e_info.orientation = HORIZONTAL;
    e_info.e_win_to = NULL;

    /* 
     * if plot_type is set, find an element window containing a plot of that type
     */
#ifdef DEBUG
    printf("plot_type %d\n", plot_type);
#endif

    if (plot_type != -1) {
	if (0 == seq_get_seq_type(seq_ids[0].seq_id, plot_type, frame, &e_info)) {
	    return e_info;
	}
    }

    /* force a new container to be created */
    if (create_new_container) 
	container_id = -1;
    else
	container_id = find_container(seq_ids, num_seqs, &orientation, 
				      &e_win_to, &c_win);
    
    /* create new container window */
    if (container_id == -1) {
	if (-1 == (container_id = new_container(interp, &c_win)))
	    return (e_info);
    }
    
    /* create new element window */
    if (-1 == (element_id = new_element(interp, &e_win)))
	return (e_info);

    e_info.element_id = element_id;
    e_info.element_win = strdup(e_win);
    e_info.container_id = container_id;
    e_info.container_win = strdup(c_win);
    e_info.orientation = orientation;

    if (e_win_to) {
	e_info.e_win_to = strdup(e_win_to);
    }
    xfree(e_win);
#ifdef FIXME
    seq_element_reg(interp, e);
#endif
    return(e_info);
}

/* 
 * return the element_info for either an existing element or create a new one
 * from a container
 */
element_info get_element_container_info(Tcl_Interp *interp,
					int container_id)
{
    int element_id;
    element_info e_info;
    char *e_win;

    e_info.element_id = -1;
    e_info.element_win = NULL;
    e_info.container_id = container_id;
    e_info.container_win = NULL;

    /* create new element window */
    if (-1 == (element_id = new_element(interp, &e_win)))
	return (e_info);

    e_info.element_id = element_id;
    e_info.element_win = strdup(e_win);
    
    return (e_info);
}

void seq_find_cursor (int element_id,
		      int direction,
		      int *seq_id,
		      int *cursor_id)
		      
{
    element *e = get_element(element_id);
    int i;

     for (i = 0; i < e->num_seqs; i++) {
	 if (e->cursor[i]->private && e->seqs[i].direction == direction) {
	     *seq_id = e->seqs[i].seq_id;
	     *cursor_id = e->cursor[i]->id;
	}
    }
   
    /* no editor cursor displayed */
    for (i = 0; i < e->num_seqs; i++) {
	if (e->seqs[i].direction == direction) {
	    *seq_id = e->seqs[i].seq_id;
	    break;
	}
    }

}

/* invoked from tcl move cursor binding command */
int seq_canvas_move_cursor(int element_id,
			   int cursor_id,
			   int pos,
			   int direction)
{
    cursor_t *cursor;
    double wx, wy;
    int seq_num = -1;
    coord *column, *row;
    element *e = get_element(element_id);

    column = e->c->column[e->column_index];
    row = e->c->row[e->row_index];

    if (direction == HORIZONTAL) {
	/* convert pixel coords into world coords */
	pixel_to_world(column->pixel, pos, pos, &wx, &wy);

	if (e->orientation & HORIZONTAL) {
	    if (wx < column->total.min) {
		wx = column->total.min;
	    }
	    if (wx > column->total.max) {
		wx = column->total.max;
	    }
	} else {
	    if (wx < e->world->total->x1) {
		wx = e->world->total->x1;
	    }
	    if (wx > e->world->total->x2) {
		wx = e->world->total->x2;
	    }
	}
	
    } else {
	/* convert pixel coords into world coords */
	pixel_to_world(row->pixel, pos, pos, &wx, &wy);

	/* have to turn world y upside down after pixel_to_world */
	wy = invert_wy(e, wy);

	if (e->orientation & VERTICAL) {
	    if (wy < row->total.min) {
		wy = row->total.min;
	    }
	    if (wy > row->total.max) {
		wy = row->total.max;
	    }
	} else {
	    if (wy < e->world->total->y1) 
		wy = e->world->total->y1;
	    if (wy > e->world->total->y2)
		wy = e->world->total->y2;
	}
    }
    
    cursor = find_cursor(&seq_num, cursor_id, -1);

#ifdef DEBUG
    printf("seq_canvas_move_cursor pos %d wx %f cursorpos %d cursor_id %d %d seq_num %d\n", pos, wx, 
	   cursor->abspos, cursor_id, cursor->id, seq_num);
    printf("seq__move_cursor pos %d wx %f\n", pos, wx);
#endif

    e->cursor_array[cursor->id].prev_pos = cursor->abspos;

    if (direction == HORIZONTAL) {
	cursor->abspos = ROUND(wx);
    } else {
	cursor->abspos = ROUND(wy);
    }
    seq_element_cursor_notify(e->c->interp, seq_num, e, cursor, 0);
    
    return 0;
}

void seq_element_start_shutdown(element *e,
				int job)
{
    seq_reg_info info;
    int c_id = e->c->id;

    info.job = SEQ_QUIT;

#ifdef DEBUG
    printf("START SHUTDOWN %s %d\n", e->win, e->seq_e_id);
#endif
    /* Tcl_VarEval(e->c->interp, "delete_element_menu ", e->win, NULL); */
    seq_result_notify(e->seq_e_id, (seq_reg_data *)&info, 1);
    delete_element(e, job);

    update_container_menu(c_id);
}

void seq_replot_element(element *e)
{
    seq_reg_plot jdata; 
    int i;
    char cmd[1024];

    if (e->type == SEQ_EDITOR)
	return;

    if (e->type == RULER_LEN || e->type == RULER_AMP) {
	sprintf(cmd, "%s delete cursor", e->win);
	Tcl_Eval(e->c->interp, cmd);
    } else {

	sprintf(cmd, "%s delete all", e->win);
	Tcl_Eval(e->c->interp, cmd);
	
	jdata.job = SEQ_PLOT;
	for (i = 0; i < e->num_results; i++) {
	    if (!e->results[i]->hidden) {
		seq_result_notify(e->results[i]->result_id, 
				  (seq_reg_data *)&jdata, 0);
	    }
	}
    }
    
    for (i = 0; i < e->num_seqs; i++) {
	if (e->cursor) {
	    seq_element_cursor_notify(e->c->interp, GetSeqNum(e->seqs[i].seq_id), 
				      e, e->cursor[i], 1);
	}   
    }
}

void set_seq_element_funcs(element *e)
{
    e->shutdown_func = seq_element_start_shutdown;
    e->replot_func = seq_replot_element;
}

void seq_element_cursor_notify(Tcl_Interp *interp,
			       int seq_num,
			       element *e,
			       cursor_t *cursor,
			       int show_all)
{
    seq_cursor_notify cn;

    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = cursor;
    cn.cursor->job = CURSOR_MOVE;
    cn.show_all = show_all; /* hack to ensure all elements are in sync */
    seq_notify(seq_num, (seq_reg_data *)&cn); 

    set_container_column_coords(e->c->interp, e->c->column[e->column_index], e->element_x(e->c->interp, e->win, 0.0));

}

void add_new_cursor(Tcl_Interp *interp,
		    element *e,
		    int seq_id,
		    int direction,
		    double min_x)
{
    cursor_t *cursor;
 
#ifdef DEBUG
    printf("add_new_cursor %s\n", e->win);
#endif
    cursor = find_element_cursor(e, seq_id, direction);
    
    if (cursor->refs == 1 && min_x > cursor->abspos) {
	cursor->abspos = min_x;
    }

    /* 
     * need to ensure done all the necessary plotting before adding the 
     * cursor
     */
    Tcl_VarEval(interp, "update idletasks ", NULL);
    seq_element_cursor_notify(interp, GetSeqNum(seq_id), e, cursor, 1);
}

void update_container(Tcl_Interp *interp,
		      int new_c_id,
		      char *new_c_win,
		      int old_c_id,
		      int new_e_id,
		      char *new_e_win,
		      int new_orientation,
		      int new_crosshair,
		      int old_e_id,
		      int result_id,
		      char *element_type,
		      char *job)
{
    seq_reg_plot jdata;
    element *e_old, *e_new;
    seq_result *s_result;
    plot_data *results;
    Graph *graph;
    double min_x, max_x, min_y, max_y;
    int i;
    Tcl_Obj *graph_obj;
    int e_type;
    int old_num_results;
    

#ifdef DEBUG
    printf("update container %d %d\n", old_e_id, new_e_id);
#endif

    if (result_id != -1) {
	s_result = seq_id_to_result(result_id);
	graph_obj = (Tcl_Obj *) s_result->data;
	graph = Tcl_GetGraphFromObj(graph_obj);
    }

    if (strcmp(element_type, "SEQ_EDITOR") == 0) {
	e_type = SEQ_EDITOR;
    } else if (strcmp(element_type, "CANVAS") == 0) {
	e_type = CANVAS;
    } else if (strcmp(element_type, "RULER_AMP") == 0) {
	e_type = RULER_AMP;
    } else if (strcmp(element_type, "RULER_LEN") == 0) {
	e_type = RULER_LEN;
    } else {
	e_type = -1;
    }

    jdata.job = SEQ_PLOT;
    if (strcmp(job, "NEW") == 0) {
	e_old = get_element(old_e_id);

	/* moving single result to new container */
	if (result_id != -1) {
#ifdef DEBUG
	    printf("MOVING SINGLE RESULT TO NEW CONTAINER\n");
#endif
	    results = find_plot_data(e_old, result_id);

	    /* if (new_orientation == e_old->orientation) { */

	    /* HACK - original data is *always* in horizontal orientation */ 
	    if (new_orientation == HORIZONTAL) {
		min_x = graph->dim.x0;
		max_x = graph->dim.x1;
		min_y = graph->dim.y0;
		max_y = graph->dim.y1;
	    } else {
		min_x = graph->dim.y0;
		max_x = graph->dim.y1;
		min_y = graph->dim.x0;
		max_y = graph->dim.x1;
	    }

	    /* if only scale in one orienation, must swap this if
	     * swapping orientation
	     */
	    if (new_orientation != e_old->orientation) {
		if (results->scale == SCALE_X)
		    results->scale = SCALE_Y;
		else if (results->scale == SCALE_Y)
		    results->scale = SCALE_X;
	    }
	    
	    /*FIXME - do crosshair? */
	    make_new_element(interp, new_c_id, new_c_win, new_e_id, new_e_win, 
			     s_result->seq_id[HORIZONTAL], results,
			     min_x, max_x, min_y, max_y, new_orientation,
			     HORIZONTAL|VERTICAL, e_type, &e_new);
	    set_seq_element_funcs(e_new);

	    /* add element pointer to seq_result structure */
	    s_result->e = e_new;
	    seq_ruler_reg(e_new);
	    seq_element_reg(interp, e_new);

	    /* HACK !!!!!!! */

	    /* if change scrollregion need to move everthing first */
	    canvas_move(interp, e_new->c->matrix[0][0], 0, 0, 72);


	    seq_result_notify(result_id, (seq_reg_data *)&jdata, 0);
	    /*
	    {
		d_box bbox;
		bbox = scale_box(e);
		e_new->scale_func(interp, e, -1, &bbox, e->pixel);
	    }
	    */
	    add_new_cursor(interp, e_new, s_result->seq_id[HORIZONTAL],
			   new_orientation, min_x);
	    
	    old_num_results = e_old->num_results;
	    remove_result_from_element(e_old, result_id);

	    if (old_num_results > 1) {
		seq_replot_element(e_old);
	    }
	} else {
#ifdef DEBUG
	    printf("MOVING ENTIRE ELEMENT TO NEW CONTAINER %d\n", new_e_id);
#endif
	    /* moving entire element */
	    if (e_old->type == SEQ_EDITOR)
		return;

	    for (i = 0; i < e_old->num_results; i++) {
		results = find_plot_data(e_old, e_old->results[i]->result_id);

		/* if only scale in one orienation, must swap this if
		 * swapping orientation
		 */
		if (new_orientation != e_old->orientation) {
		    if (results->scale == SCALE_X)
			results->scale = SCALE_Y;
		    else if (results->scale == SCALE_Y)
			results->scale = SCALE_X;
		}
	    }
	    move_element_to_new_container(interp, old_e_id, new_c_id, 
					  new_c_win, old_c_id, new_e_win,
					  new_e_id, new_orientation); 
	    seq_ruler_reg(e_old);

	    for (i = 0; i < e_old->num_results; i++) {
		Tcl_Obj *graph_obj;
		s_result = seq_id_to_result(e_old->results[i]->result_id);
		graph_obj = (Tcl_Obj *) s_result->data;
		graph = Tcl_GetGraphFromObj(graph_obj);
		/* HACK - original data is *always* in horizontal orientation */
		if (new_orientation == HORIZONTAL) {
		    min_x = graph->dim.x0;
		    max_x = graph->dim.x1;
		    min_y = graph->dim.y0;
		    max_y = graph->dim.y1;
		} else {
		    min_x = graph->dim.y0;
		    max_x = graph->dim.y1;
		    min_y = graph->dim.x0;
		    max_y = graph->dim.x1;
		}
		add_new_cursor(interp, e_old, s_result->seq_id[HORIZONTAL],
			       new_orientation, min_x);

	    }
	}
    } else if (strcmp(job, "SUPERIMPOSE") == 0) {
#ifdef DEBUG
	printf("SUPERIMPOSE \n");
#endif
	/* moving single result to existing element */
	if (result_id != -1) {
	    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
	    graph = Tcl_GetGraphFromObj(graph_obj);
	    e_new = get_element(new_e_id);
	    e_old = get_element(old_e_id);
	    results = find_plot_data(e_old, result_id);
	    
	    min_x = graph->dim.x0;
	    max_x = graph->dim.x1;
	    min_y = graph->dim.y0;
	    max_y = graph->dim.y1;

	    if (e_new->orientation != e_old->orientation) {
		if (results->scale == SCALE_X)
		    results->scale = SCALE_Y;
		else if (results->scale == SCALE_Y)
		    results->scale = SCALE_X;
	    }

	    add_result_to_element(e_new, results,  min_x, min_y, max_x, max_y, 
				  e_new->orientation, e_old->type);

	    /* if superimposing, must maintain same orientation of seqs */
	    for (i = 0; i < e_old->num_seqs; i++) {
		add_seq_id_to_element(e_new, e_old->seqs[i].seq_id, 
				      e_old->seqs[i].direction);		
	    }

	    s_result->e = e_new;

	    old_num_results = e_old->num_results;
	    remove_result_from_element(e_old, result_id);

	    if (old_num_results > 1) 
		e_old->replot_func(e_old);

	    seq_result_notify(result_id, (seq_reg_data *)&jdata, 0);
	} else {
	    /* moving entire element to existing element */
#ifdef DEBUG
	    printf("superimpose element\n");
#endif
	    e_new = get_element(new_e_id);
	    e_old = get_element(old_e_id);
	    
	    for (i = 0; i < e_old->num_results; i++) {
		Tcl_Obj *graph_obj;
		s_result = seq_id_to_result(e_old->results[i]->result_id);
		graph_obj = (Tcl_Obj *) s_result->data;
		graph = Tcl_GetGraphFromObj(graph_obj);
		s_result->e = e_new;
		results = find_plot_data(e_old, e_old->results[i]->result_id);

		if (e_new->orientation == HORIZONTAL) {
		    min_x = graph->dim.x0;
		    max_x = graph->dim.x1;
		    min_y = graph->dim.y0;
		    max_y = graph->dim.y1;
		} else {
		    min_x = graph->dim.y0;
		    max_x = graph->dim.y1;
		    min_y = graph->dim.x0;
		    max_y = graph->dim.x1;
		}

		if (e_new->orientation != e_old->orientation) {
		    if (results->scale == SCALE_X)
			results->scale = SCALE_Y;
		    else if (results->scale == SCALE_Y)
			results->scale = SCALE_X;
		}

		add_result_to_element(e_new, results, min_x, min_y, 
				      max_x, max_y, e_new->orientation, 
				      e_old->type);
	    }
	    for (i = 0; i < e_old->num_seqs; i++) {
		add_seq_id_to_element(e_new, e_old->seqs[i].seq_id,
				      e_old->seqs[i].direction);
	    }
	    /*
	    for (i = 0; i < e_old->num_results; i++) {
		remove_result_from_element(e_old, e_old->results[i]->result_id);
	    }
	    */
	    e_old->num_results = 0;
	    e_old->shutdown_func(e_old, WIN_ONLY); 
	    seq_replot_element(e_new);
	}
    }
}

void seq_plot_graph_func(void *obj, seq_reg_plot *plot)
{
    seq_result *s_result = (seq_result *) obj;
    element *e = s_result->e;
    Tcl_Interp *interp = e->c->interp;
    d_box bbox;
    plot_data *result = find_plot_data(e, s_result->id);
    int i;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
    Graph *graph;
    double m, c, pos;
    int win_height;
    char cmd[1024];

#ifdef DEBUG
    printf("seq_plot_graph_func id %d ewin %s\n", s_result->id, e->win);
#endif

    if (e->orientation & HORIZONTAL) {
	win_height = e->element_height(interp, e->win);
    } else {
	win_height = e->element_width(interp, e->win);
    }

    graph = Tcl_GetGraphFromObj(graph_obj);

    /* if tick_ht is a percentage (0.0 - 1.0), convert to pixels */

    /* FIXME - need to test this */
#if 0
    for (i = 0; i < result->n_configure; i++) {
	if (result->configure[i]->height > 0.0 && result->configure[i]->height <= 1.0) {
	    for (k = 0; k < graph->n_darrays; k++) {
		for (j = 0; j < graph->d_arrays[k].n_dlines; j++) {
		    graph->d_arrays[k].d_array[j].y1 *= result->configure[i]->height * win_height; 
		}
	    }
	    graph->dim.y1 *= win_height;
	}
    }
#endif

    sprintf(cmd, "%s delete %s", e->win, result->tags);
    Tcl_Eval(interp, cmd);

    if (s_result->graph == GRAPH) {
	create_graph(interp, e->win, (Tcl_Obj *) s_result->data, result->line_width, result->colour, result->tags, e->orientation);
    } else {
	create_canv_line(interp, e->win, graph, result, result->line_width, result->colour, result->tags, e->orientation);
    }

    bbox = scale_box(e);

    /* use the data min and max so all graphs are to their own scale */
    if (e->orientation & HORIZONTAL) {
	m = (e->world->visible->y2 - e->world->visible->y1) /
	    (e->world->total->y2 - e->world->total->y1);
	
	c = e->world->visible->y2 - (m * e->world->total->y2);
	
	bbox.y1 = m * graph->dim.y0 + c;
	bbox.y2 = m * graph->dim.y1 + c;
    } 
    if (e->orientation & VERTICAL) {
	m = (e->world->visible->x2 - e->world->visible->x1) /
	    (e->world->total->x2 - e->world->total->x1);
	
	c = e->world->visible->x2 - (m * e->world->total->x2);
	
	/* note that graph->dim never changes orientation */
	if (!(e->orientation & HORIZONTAL)) {
	    /* vertical only */
	    bbox.x1 = m * graph->dim.y0 + c;
	    bbox.x2 = m * graph->dim.y1 + c;
	} else {
	    /* dot plot */
	    bbox.x1 = m * graph->dim.x0 + c;
	    bbox.x2 = m * graph->dim.x1 + c;
	}
    }

    e->scale_func(interp, e, s_result->id, &bbox, e->pixel); 

    e->scrollregion_func(interp, e, e->world->total, 
			 e->c->column[e->column_index]->pixel,
			 e->c->row[e->row_index]->pixel);    
    
    for (i = 0; i < result->n_configure; i++) {
	if (result->configure[i]->position != -1.0) {
	    if (e->orientation & HORIZONTAL) {
		pos = win_height * result->configure[i]->position;
		if (result->configure[i]->y_direction == '+') {
		    if (result->configure[i]->scroll) {
			pos = pos - win_height;
		    } else {
			/* pos = pos - graph->dim.y1; */
			
		    }
		} else {
		    if (s_result->graph == GRAPH) {
			sprintf(cmd, "%s itemconfigure id%d -inverty 0", 
				e->win, s_result->id);
			if (TCL_ERROR == Tcl_Eval(interp, cmd))
			    printf("ERROR in configure %s\n", interp->result);
		    }
		}
#ifdef DEBUG
		printf("POS %f %f %f win_height %d id %d\n", pos, graph->dim.y1, graph->dim.y0, win_height, s_result->id);
#endif
		canvas_move(interp, e, s_result->id, 0, pos);
	    }
	    if (e->orientation & VERTICAL) {
		pos = win_height * result->configure[i]->position;
		if (result->configure[i]->y_direction == '+') {
		    if (s_result->graph == GRAPH) {
			if (result->configure[i]->scroll) {
			    pos = pos - win_height;
			} else {
			    /* pos = pos - graph->dim.y1; */
			}
			sprintf(cmd, "%s itemconfigure id%d -invertx 1", 
				e->win, s_result->id);
			if (TCL_ERROR == Tcl_Eval(interp, cmd))
			    printf("ERROR in configure %s\n", interp->result);
		    }
		}
#ifdef DEBUG
		printf("POSY %f %f %f\n", pos, graph->dim.y1, graph->dim.y0);
#endif
		canvas_move(interp, e, s_result->id, pos, 0);
	    }
	}
    }
}

void config_result(int e_id,
		   int result_id,
		   int line_width,
		   char *colour)
{
    element *e = get_element(e_id);
    plot_data *result = find_plot_data(e, result_id);

    result->line_width = line_width;
    result->colour = strdup(colour);
}

void seq_ruler_reg(element *e) 
{
    container *c = e->c;

    if (e->orientation & HORIZONTAL) {
	if (c->column[e->column_index]->ruler) {
	    if (!c->column[e->column_index]->ruler_reg) {
		c->column[e->column_index]->ruler_reg = 1;
		seq_element_reg(c->interp, c->column[e->column_index]->ruler);
		set_seq_element_funcs(c->column[e->column_index]->ruler);
	    }
	}
    }
    if (e->orientation & VERTICAL) {
	if (c->row[e->row_index]->ruler) {
	    if (!c->row[e->row_index]->ruler_reg) {
		c->row[e->row_index]->ruler_reg = 1;
		seq_element_reg(c->interp, c->row[e->row_index]->ruler);
		set_seq_element_funcs(c->row[e->row_index]->ruler);
	    }
	}
    }

    /* amplitude ruler */
    if (e->ruler_id != -1) {
	element *e_ruler = get_element(e->ruler_id);
	set_seq_element_funcs(e_ruler);
    }
}

/*
 * Superimpose a result onto an existing window
 */
void seq_superimpose(int result_id,
		     element *e_old,
		     element *e_new)
{
    seq_result *s_result = seq_id_to_result(result_id);
    double p2, q2;
    double m, c; 
    double o_wy0, o_wy1, y0, y1;
    double n_wy0, n_wy1;
    Graph *graph;
    plot_data *result = find_plot_data(e_old, result_id);

    if (s_result->graph == GRAPH) {
	Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
	graph = Tcl_GetGraphFromObj(graph_obj);
	
	y0 = graph->dim.y0;
	y1 = graph->dim.y1;
    }    

    o_wy0 = e_old->min_y;
    o_wy1 = e_old->max_y;

    n_wy0 = e_new->min_y;
    n_wy1 = e_new->max_y;

    p2 = (n_wy1 - n_wy0) * (y0 - o_wy0) / (o_wy1 - o_wy0) + n_wy0;
    q2 = (n_wy1 - n_wy0) * (y1 - o_wy0) / (o_wy1 - o_wy0) + n_wy0;
    
    if (y0 - y1 != 0) {
	m = (p2 - q2) / (y0 - y1);
    } else {
	m = 0;
    }
    c = p2 - (m * y0);
 
    result->sf_c = (m * result->sf_c) + c;
    result->sf_m *= m;

}

int init_seq_element(Tcl_Interp *interp,
		     seq_result *s_result,
		     int seq_id_h,
		     int seq_id_v,
		     plot_data *result,
		     int container_id,
		     int element_id,
		     char *e_win,
		     char *c_win,
		     int orientation,
		     int crosshair,
		     Graph *graph,
		     char *element_type)
{
    element *e;
    container *c;
    char tags[10];
    cursor_t *cursor_h, *cursor_v;
    int new_element = 0;
    int e_type;


    if (NULL == (e = (get_element(element_id)))) {
	e = create_element(interp, container_id, element_id, e_win, orientation, crosshair);
	new_element = 1;
    }
    /* add element pointer to seq_result structure */
    s_result->e = e;

    if (strcmp(element_type, "SEQ_EDITOR") == 0) {
	e_type = SEQ_EDITOR;
    } else if (strcmp(element_type, "CANVAS") == 0) {
	e_type = CANVAS;
    } else if (strcmp(element_type, "RULER_AMP") == 0) {
	e_type = RULER_AMP;
    } else if (strcmp(element_type, "RULER_LEN") == 0) {
	e_type = RULER_LEN;
    } else {
	e_type = -1;
    }

    if (-1 == (add_result_to_element(e, result, graph->dim.x0, graph->dim.y0, 
				     graph->dim.x1, graph->dim.y1, orientation,
				     e_type))) {
	return -1;
    }

    if (e->orientation & HORIZONTAL) 
	add_seq_id_to_element(e, seq_id_h, HORIZONTAL);

    if (e->orientation & VERTICAL) 
	add_seq_id_to_element(e, seq_id_v, VERTICAL);

    if (new_element) {
	initCanvas(interp, e->pixel, e->win);
   
	if (e->orientation & HORIZONTAL && e->orientation & VERTICAL) {
	    add_element_to_container(interp, container_id, c_win, e, 
				     graph->dim.x0, graph->dim.x1,  
				     graph->dim.y0, graph->dim.y1);
	} else if (e->orientation & HORIZONTAL) {
	    add_element_to_container(interp, container_id, c_win, e, 
				     graph->dim.x0, graph->dim.x1,
				     INT_MAX, INT_MIN);
	} else if (e->orientation & VERTICAL) {
	    add_element_to_container(interp, container_id, c_win, e, 
				     INT_MAX, INT_MIN,
				     graph->dim.y0, graph->dim.y1);
	} else if (e->orientation & CIRCLE) {
	    add_element_to_container(interp, container_id, c_win, e, 
				     graph->dim.x0, graph->dim.x1, 
				     graph->dim.y0, graph->dim.y1);
	}
	/* set seq specific functions */
	set_seq_element_funcs(e);
    }

    e->world->visible->x1 = e->world->total->x1;
    e->world->visible->x2 = e->world->total->x2;
    e->world->visible->y1 = e->world->total->y1;
    e->world->visible->y2 = e->world->total->y2;
	
    if (e->orientation & HORIZONTAL) {
	e->world->visible->x1 = e->c->column[e->column_index]->total.min;
	e->world->visible->x2 = e->c->column[e->column_index]->total.max;
    } 
    if (e->orientation & VERTICAL) {
	e->world->visible->y1 = e->c->row[e->row_index]->total.min;
	e->world->visible->y2 = e->c->row[e->row_index]->total.max;
    }

    set_pixel_coords(e->world->visible->x1, e->world->visible->y1, 
		     e->world->visible->x2, e->world->visible->y2, e->pixel);
    
    c = get_container(container_id);

    sprintf(tags, "id%d", s_result->id);

    e->replot_func(e);

    if (new_element) {
	seq_ruler_reg(e);
	seq_element_reg(interp, e);
    }

    if (e->orientation & HORIZONTAL) {
	cursor_h = find_element_cursor(e, seq_id_h, HORIZONTAL);
	if (cursor_h->refs <= 2 && graph->dim.x0 > cursor_h->abspos) {
	    cursor_h->abspos = graph->dim.x0;
	}
    }
    if (e->orientation & VERTICAL) {
	cursor_v = find_element_cursor(e, seq_id_v, VERTICAL);
	if (cursor_v->refs <= 2 && graph->dim.x0 > cursor_v->abspos) {
	    cursor_v->abspos = graph->dim.x0;
	}
    }    

    /* 
     * need to ensure done all the necessary plotting before adding the 
     * cursor
     */
    Tcl_VarEval(interp, "update idletasks ", NULL);
    if (e->orientation & HORIZONTAL)
	seq_element_cursor_notify(e->c->interp, GetSeqNum(seq_id_h), e, 
				  cursor_h, 1);
    if (e->orientation & VERTICAL)
	seq_element_cursor_notify(e->c->interp, GetSeqNum(seq_id_v), e, 
				  cursor_v, 1);
    /* remove all current zooming info */
    freeZoom(&e->zoom);

    /* add first zoom */
    pushZoom(&e->zoom, e->world->visible);

    return 0;
}

int get_element_type(element *e)
{
    int i;
    int type = 0;
    seq_result *s_result;

    for (i = 0; i < e->num_results; i++) {
	if (e->results[i]->result_id != -1) {
	    s_result = seq_id_to_result(e->results[i]->result_id);
	    type |= s_result->type;
	}
    }
    return type;
}

int get_element_gr_type(element *e)
{
    int i;
    int type = 0;
    seq_result *s_result;

    if (e->type == SEQ_EDITOR)
	return SEQ_EDITOR;

    for (i = 0; i < e->num_results; i++) {
	if (e->results[i]->result_id != -1) {
	    s_result = seq_id_to_result(e->results[i]->result_id);
	    type |= s_result->gr_type;
	}
    }
#ifdef DEBUG
    printf("final type %d\n", type);
#endif
    return type;
}

int compare_g_pt(const void *p1,
		 const void *p2)
{
    g_pt *i1 = (g_pt *) p1, *i2 = (g_pt *) p2;
    return i1->x - i2->x;
}

int compare_gd_line(const void *p1,
		   const void *p2)
{
    gd_line *i1 = (gd_line *) p1, *i2 = (gd_line *) p2;
    return i1->x0 - i2->x0;
}
