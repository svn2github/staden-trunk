#include <tk.h>
#include <stdlib.h>

#include "canvas_box.h"
#include "misc.h"
#include "io-reg.h"
#include "gap_canvas_box.h"
#include "consistency_canvas_box.h"
#include "consistency_display.h"
#include "tcl_utils.h"
#include "gap_globals.h"
#include "template_display.h"
#include "ruler_tick.h"

static void consistency_callback(GapIO *io, int contig, void *fdata,
			      reg_data *jdata);

int update_consistency_display(Tcl_Interp *interp, GapIO *io,
			       obj_consistency_disp *c);

void consistency_contig_offsets(GapIO *io,
				c_offset *contig_offset,
				int *contigs,
				int num_contigs)
{
    int i;

    contig_offset[contigs[0]].offset = 0;
    contig_offset[contigs[0]].gap = 0;
    for (i = 1; i < num_contigs; i++) {

        contig_offset[contigs[i]].gap = 0;
	contig_offset[contigs[i]].offset = 
	    contig_offset[contigs[i-1]].offset + 
		    ABS(io_clength(io, contigs[i-1]));
    }
}

void print_contig_offset(obj_consistency_disp *c)
{
    int i;

    printf("PRINT_CONTIG_OFFSET\n");
    for (i = 0; i < c->num_contigs; i++) {

      printf("contigs[%d]= %d offset= %d\n", 
	     i, c->contigs[i], c->contig_offset[c->contigs[i]].offset);
    }
}

static void consistency_update_cursor(GapIO *io, 
				      obj_consistency_disp *c,
				      int contig, 
				      cursor_t *cursor,
				      int show_cursor)
{
    int index, i;
    int win_num;

#ifdef DEBUG
    printf("consistency_update_cursor\n");
#endif

    index = -1;
    for (i = 0; i < c->num_contigs; i++) {
	if (c->contigs[i] == contig) {
	    index = i;
	    break;
	}
    }

    if (index == -1)
        return;

    /* doesn't matter which window - only using x coords */
    win_num = 0;

    consistency_cursor_refresh(c->interp, io, c, c->contigs[index],
			       cursor, c->cursor[index],
			       c->win_list[win_num]->canvas, c->win_list, 
			       c->num_wins,
			       c->id, c->contig_offset[c->contigs[index]].offset,
			       &c->cursor_visible[index], 
			       c->win_list[win_num]->world, 
			       show_cursor);
}

void consistency_update_cursors(GapIO *io, 
				obj_consistency_disp *c,
				int show_cursor)
{
    int i;
    for (i = 0; i < c->num_contigs; i++) {
	consistency_update_cursor(io, c, c->contigs[i], c->cursor[i], 
				  show_cursor);
    }
}

static void consistency_join(Tcl_Interp *interp,
			     GapIO *io,
			     obj_consistency_disp *c, 
			     int old_contig,
			     int new_contig)
{
    int i, j;
    int length, old_found = 0;
    int remove_dup = -1;
    static char old_contig_name[DB_NAMELEN+1];
    static char new_contig_name[DB_NAMELEN+1];
    int win_num = 0; /* any win_num will do */

    strcpy(old_contig_name, 
	   get_read_name(io, io_clnbr(io, old_contig)));
    strcpy(new_contig_name, 
	   get_read_name(io, io_clnbr(io, new_contig)));
  
    /* replace the old contig with the new contig */
    for (i = 0; i < c->num_contigs; i++) {
	if (ABS(c->contigs[i]) == old_contig) {
	    old_found = 1;

	    /* check for duplicates and set flag */
	    for (j = 0; j < c->num_contigs; j++) {
		 if (ABS(c->contigs[j]) == new_contig) {
		     remove_dup = i;
		     break;
		 }
	     }
	    /* no duplicates, don't set flag */
	    c->contigs[i] = new_contig;
	    break;
	}
    }

    if (remove_dup != -1) {
	/* remove any duplicates */
	delete_contig_cursor(io, old_contig, c->cursor[remove_dup]->id, 0);
	consistency_cursor_refresh(interp, io, c, new_contig,
				   c->cursor[remove_dup], c->cursor[remove_dup],
				   c->win_list[win_num]->canvas, c->win_list, 
				   c->num_wins,
				   c->id, c->contig_offset[new_contig].offset,
				   &c->cursor_visible[remove_dup], 
				   c->win_list[win_num]->world, 1);
	
	length = c->num_contigs - remove_dup - 1;
	memmove(&(c->contigs)[remove_dup],
		&(c->contigs)[remove_dup+1],
		length * sizeof(int));
	memmove(&(c->cursor)[remove_dup],
		&(c->cursor)[remove_dup+1],
		length * sizeof(int));
	memmove(&(c->cursor_visible)[remove_dup],
		&(c->cursor_visible)[remove_dup+1],
		length * sizeof(int));
	{
	  /* update all displays in consistency plot */
	  reg_generic gen;
	  gen.job = REG_GENERIC;
	  gen.data = (void *)remove_dup;
	  gen.task = TASK_CONS_JOIN;
	  for (i = 0; i < c->num_wins; i++) {
	    result_notify(io, c->win_list[i]->id, (reg_data *)&gen, 0); 
	  }
	}

	c->num_contigs--;
    }
}

void clear_consistency(GapIO *io, obj_consistency_disp *c)
{
    reg_quit rq;
    int i;

    rq.job = REG_QUIT;
    rq.lock = REG_LOCK_WRITE;

    for (i = 0; i < c->num_wins; i++) {
	if (c->win_list[i]->id != c->id) {
	    int num_wins = c->num_wins;
	    result_notify(io, c->win_list[i]->id, (reg_data *)&rq, 1);
	    i -= num_wins - c->num_wins;
	}
    }
}

/*
 * Removes the consistency display (and unplots etc).
 */
void consistency_shutdown(GapIO *io, obj_consistency_disp *c) {
    int i;
    char cmd[1024];

    clear_consistency(io, c);

    for (i = 0; i < c->num_contigs; i++) {
	contig_deregister(io, c->contigs[i], consistency_callback, (void *)c);
	delete_contig_cursor(io, c->contigs[i], c->cursor[i]->id, 0);
    }
    xfree(c->contigs);

    sprintf(cmd, "DeleteConsistencyDisplay %s\n", c->frame);
    if (TCL_ERROR == Tcl_Eval(c->interp, cmd))
	printf("consistency_shutdown:%s\n", Tcl_GetStringResult(c->interp));

    if (c->orig_total)
        xfree(c->orig_total);

    if (c->contig_offset)
	xfree(c->contig_offset);

    if (c->ruler_coord) {
	for (i = 0; i < c->num_contigs; i++) {
	    xfree(c->ruler_coord[i].type);
	}
	xfree(c->ruler_coord);
    }
    
    for (i = 0; i < c->num_wins; i++) {
        delete_consistency_window(c, i);
    } 
    free_win_list(c->win_list, c->num_wins);

    if (c->ruler->colour) free(c->ruler->colour);
    if (c->ruler->tick.t.colour) free(c->ruler->tick.t.colour);
    xfree(c->ruler);

    if (c->xhair.colour) free(c->xhair.colour);

    xfree(c->cursor);
    xfree(c->cursor_visible);
    xfree(c);
}

/* 
 * scale visible world by amount when 10% or 50% buttons are pressed
 */
void consistencyScaleZoom(Tcl_Interp *interp,
			  win *win_list,
			  d_box *zoom,
			  float amount)
{
    int width, height;

    /* 
     * think I prefer the second method since this doesn't have to convert 
     * from world coords
     */
#if 0
    double length;
    double dist;
    double wx1, wy1, wx2, wy2;
    length = win_list->world->visible->x2 - win_list->world->visible->x1 + 1;
    
    dist = length * amount;
    wx1 = win_list->world->visible->x1 + dist;
    wx2 = win_list->world->visible->x2 - dist;
    
    length = win_list->world->visible->y2 - win_list->world->visible->y1 + 1;
    dist = length * amount;
    
    wy1 = win_list->world->visible->y1 + dist;
    wy2 = win_list->world->visible->y2 - dist;
    
    WorldToCanvas(win_list->canvas, wx1, wy1, &zoom->x1, &zoom->y1);
    WorldToCanvas(win_list->canvas, wx2, wy2, &zoom->x2, &zoom->y2);

    zoom->x1 = canvasx(interp, win_list->window, zoom->x1);
    zoom->x2 = canvasx(interp, win_list->window, zoom->x2);
#endif

    width = win_list->canvas->width + 1;
    height = win_list->canvas->height + 1;

    zoom->x1 = canvasx(interp, win_list->window, width * amount);
    zoom->x2 = canvasx(interp, win_list->window, width * (1.0 - amount));
    zoom->y1 = canvasy(interp, win_list->window, height * amount);
    zoom->y2 = canvasy(interp, win_list->window, height * (1.0 - amount));
}

void consistencyZoom(obj_consistency_disp *c,
		     box *zoom,
		     char scroll,
		     float amount,
		     int result_id)
{
    d_box *bbox;
    char cmd[1024];
    double x1, y1, x2, y2;
    int i;
    int win_num;

    if (NULL == (bbox = (d_box *)xmalloc(sizeof(d_box))))
	return;

    for (i = 0; i < c->num_wins; i++) {

      if (amount != -1) {
	  consistencyScaleZoom(c->interp, c->win_list[i], bbox, amount);
      } else {
	  bbox->x1 = (double)zoom->x1;
	  bbox->y1 = (double)zoom->y1;
	  bbox->x2 = (double)zoom->x2;
	  bbox->y2 = (double)zoom->y2;
      }

      x1 = c->win_list[i]->world->visible->x1;
      y1 = c->win_list[i]->world->visible->y1;
      x2 = c->win_list[i]->world->visible->x2;
      y2 = c->win_list[i]->world->visible->y2;

      CanvasToWorld(c->win_list[i]->canvas, (int)bbox->x1, (int)bbox->y1, 
		    &c->win_list[i]->world->visible->x1, 
		    &c->win_list[i]->world->visible->y1);
      CanvasToWorld(c->win_list[i]->canvas, (int)bbox->x2, (int)bbox->y2, 
		    &c->win_list[i]->world->visible->x2, 
		    &c->win_list[i]->world->visible->y2);
    
      /*
       * need to take into account how the originating zooming window scrolled
       * so that eg a window that only scrolls in x only zooms in x
       */  
      if (scroll == 'x' || scroll == 'n') {
	c->win_list[i]->world->visible->y1 = y1;
	c->win_list[i]->world->visible->y2 = y2;
	bbox->y1 = 0.0;
	bbox->y2 = 0.0;
      }
      if (scroll == 'y' || scroll == 'n') {
	c->win_list[i]->world->visible->x1 = x1;
	c->win_list[i]->world->visible->x2 = x2;
	bbox->x1 = 0.0;
	bbox->x2 = 0.0;
      }
      /* 
       * also if drawn a zooming box, only want to zoom in y that window 
       */
      
      if (amount == -1) {
	  win_num = get_consistency_win_num(c, result_id);
	  
	  if (i != win_num) {
	      c->win_list[i]->world->visible->y1 = y1;
	      c->win_list[i]->world->visible->y2 = y2;
	      bbox->y1 = 0.0;
	      bbox->y2 = 0.0;
	  }
      }
      
      SetCanvasCoords(c->interp, c->win_list[i]->world->visible->x1, 
		      c->win_list[i]->world->visible->y1, 
		      c->win_list[i]->world->visible->x2, 
		      c->win_list[i]->world->visible->y2, 
		      c->win_list[i]->canvas);

      scaleCanvas(c->interp, &c->win_list[i], 1, "all", bbox,
		  c->win_list[i]->canvas);
      scrollRegion(c->interp, &c->win_list[i], 1, 
		   c->win_list[i]->world->total, c->win_list[i]->canvas);

      pushZoom(&c->win_list[i]->zoom, c->win_list[i]->world->visible);
    
      sprintf(cmd, "%s canvasx 0\n", c->win_list[i]->window);
      Tcl_Eval(c->interp, cmd);
      c->win_list[i]->canvas->x = atoi(Tcl_GetStringResult(c->interp));
    }
    xfree(bbox);
}

void consistencyZoomback(obj_consistency_disp *c)
{
    int i;
    box *zoom;
    d_box *bbox;
    char cmd[1024];
    
#ifdef DEBUG
    printf("consistencyZoomback\n");
#endif
    if (NULL == (zoom = (box *)xmalloc(sizeof(box))))
	return;
	
    if (NULL == (bbox = (d_box *)xmalloc(sizeof(d_box))))
	return;
  
    for (i = 0; i < c->num_wins; i++) {
      popZoom(&c->win_list[i]->zoom);
      if (examineZoom(c->win_list[i]->zoom) == NULL) {
	return;
      }
      memcpy(c->win_list[i]->world->visible,  
	     examineZoom(c->win_list[i]->zoom), 
	     sizeof(d_box));

      WorldToCanvas(c->win_list[i]->canvas, 
		    c->win_list[i]->world->visible->x1, 
		    c->win_list[i]->world->visible->y1, 
		    &bbox->x1, &bbox->y1);
      WorldToCanvas(c->win_list[i]->canvas, 
		    c->win_list[i]->world->visible->x2, 
		    c->win_list[i]->world->visible->y2, 
		    &bbox->x2, &bbox->y2);
	
      scaleCanvas(c->interp, &c->win_list[i], 1, "all", 
		  bbox, c->win_list[i]->canvas);

      /* set new canvas coords from zoom */
      SetCanvasCoords(c->interp, c->win_list[i]->world->visible->x1, 
		      c->win_list[i]->world->visible->y1, 
		      c->win_list[i]->world->visible->x2, 
		      c->win_list[i]->world->visible->y2, 
		      c->win_list[i]->canvas);
      
      scrollRegion(c->interp, &c->win_list[i], 1, 
		   c->win_list[i]->world->total, c->win_list[i]->canvas);
      
      sprintf(cmd, "%s canvasx 0", c->win_list[i]->window);
      Tcl_Eval(c->interp, cmd);
      c->win_list[i]->canvas->x = atoi(Tcl_GetStringResult(c->interp));
    }

    xfree(zoom);
    xfree(bbox);
}

void consistency_canvasScrollX(Tcl_Interp *interp,
			       obj_consistency_disp *c,
			       win **win_list,
			       int num_wins,
			       char *scroll_args)
{
    char cmd[1024];
    int i;
    double wy;

    if (num_wins == 0) 
	return;

    for (i = 0; i < num_wins; i++) {

	if ((win_list[i]->scroll == 'x') || (win_list[i]->scroll == 'b')) {
	    sprintf(cmd, "eval %s xview %s", win_list[i]->window, scroll_args);
	    Tcl_Eval(interp, cmd);
	}
	/* find new left hand edge of canvas in canvasx coords */
	Tcl_VarEval(interp, win_list[i]->window, " canvasx 0", NULL);
	win_list[i]->canvas->x = atoi(Tcl_GetStringResult(interp));

	/* find new left and right edges of canvas in world coords */
	CanvasToWorld(win_list[i]->canvas, win_list[i]->canvas->x, 0, 
		      &win_list[i]->world->visible->x1, &wy);
	CanvasToWorld(win_list[i]->canvas, 
		      win_list[i]->canvas->x + win_list[i]->canvas->width,
		      0, &win_list[i]->world->visible->x2, &wy);
#ifdef DEBUG
	printf("scrollx %d %f %f %d\n", i, win_list[i]->world->visible->x1,
	       win_list[i]->world->visible->x2, win_list[i]->canvas->x);
#endif
	SetCanvasCoords(interp, win_list[i]->world->visible->x1, 
			win_list[i]->world->visible->y1, 
			win_list[i]->world->visible->x2, 
			win_list[i]->world->visible->y2, 
			win_list[i]->canvas);
    }
}

/* only scroll single window */
void consistencyScrollY(Tcl_Interp *interp,
			char *window,
			int scroll,
			d_box *visible,
			CanvasPtr *canvas,
			char *scroll_args)
{
    char cmd[1024];
    double wx;

    if ((scroll == 'y') || (scroll == 'b')) {
      sprintf(cmd, "eval %s yview %s", window, scroll_args);
      Tcl_Eval(interp, cmd);
    }

    /* find new top edge of canvas in canvasy coords */
    Tcl_VarEval(interp, window, " canvasy 0", NULL);
    canvas->y = atoi(Tcl_GetStringResult(interp));

    /* find new top and bottom edges of canvas in world coords */
    CanvasToWorld(canvas, 0, canvas->y, &wx, &visible->y1);
    CanvasToWorld(canvas, 0, canvas->y + canvas->height, &wx, &visible->y2);
    SetCanvasCoords(interp, visible->x1, visible->y1, 
		    visible->x2, visible->y2, canvas);
}

void consistency_resizeCanvas(Tcl_Interp *interp,
			       obj_consistency_disp *c,
			       win **win_list,
			       int num_wins)
{
    d_box *bbox;
    int height, width;
    int i;
	      
    if (num_wins == 0) 
	return;

    if (NULL == (bbox = (d_box *)xmalloc(sizeof(d_box))))
	return;

    for (i = 0; i < num_wins; i++) {

        bbox->x1 = (double)win_list[i]->canvas->x;
	bbox->y1 = (double)win_list[i]->canvas->y;
	bbox->x2 = (double)(win_list[i]->canvas->x + 
			    win_list[i]->canvas->width);
	bbox->y2 = (double)(win_list[i]->canvas->y + 
			    win_list[i]->canvas->height);

	Tcl_VarEval(interp, "winfo width ", win_list[i]->window, NULL);
	width = atoi(Tcl_GetStringResult(interp)) - 1;

	Tcl_VarEval(interp, "winfo height ", win_list[i]->window, NULL);
	height = atoi(Tcl_GetStringResult(interp)) - 1;

#ifdef DEBUG
	printf("resizeCanvas \n");
	printf("window %s width %d %d height %d %d\n", 
	       win_list[i]->window, width, win_list[i]->canvas->width, 
	       height, win_list[i]->canvas->height);
#endif
	/* if changed height or width */
	if ((height != win_list[i]->canvas->height) || 
	    (width != win_list[i]->canvas->width)) {

	    win_list[i]->canvas->height = height;
	    win_list[i]->canvas->width = width;
	    SetCanvasCoords(interp, win_list[i]->world->visible->x1, 
			    win_list[i]->world->visible->y1, 
			    win_list[i]->world->visible->x2, 
			    win_list[i]->world->visible->y2, 
			    win_list[i]->canvas);

	    scaleCanvas(interp, &win_list[i], 1, "all", bbox, 
		      win_list[i]->canvas);
	    scrollRegion(interp, &win_list[i], 1, win_list[i]->world->total, 
			 win_list[i]->canvas);
	}
    }
    xfree(bbox);
}

void display_consistency_ruler(GapIO *io,
			       Tcl_Interp *interp,
			       obj_consistency_disp *c)
{
    int win_num;
    int i;

    if (c->ruler_coord) {
      for (i = 0; i < c->num_contigs; i++) {
	xfree(c->ruler_coord[i].type);
      }
      xfree(c->ruler_coord);
      c->ruler_coord = NULL;
    }
    if (c->configs[CONS_RULER]) {
        win_num = get_consistency_win_num(c, c->id);
        display_ruler(interp, io, c->win_list[win_num]->canvas, 
		      c->contig_offset, 
		      c->contigs, c->num_contigs, c->configs[CONS_RULER], 
		      c->configs[CONS_TICKS], c->ruler, c->frame, 
		      &c->ruler_coord);

	scaleSingleCanvas(c->interp, c->win_list[win_num]->world,  
			  c->win_list[win_num]->canvas, c->ruler->window, 'x', 
			  "all");
	consistency_update_cursors(io, c, 0);
    }
}

static void consistency_renumber(GapIO *io,
				 obj_consistency_disp *c,
				 int old_contig,
				 int new_contig)
{
    int i;

    for (i = 0; i < c->num_contigs; i++) {
	if (ABS(c->contigs[i]) == old_contig) {
	    c->contigs[i] = c->contigs[i] > 0 ? new_contig : 
		-new_contig;
	    break;
	}
    }
}

/*
 * Callback for the consistency display
 */
static void consistency_callback(GapIO *io, int contig, void *fdata,
			      reg_data *jdata) {
    obj_consistency_disp *c = (obj_consistency_disp *)fdata;
    char cmd[1024];
    int i;

    /* printf("consistency_callback job %d contig %d\n", jdata->job, contig); */
    switch(jdata->job) {
    case REG_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "Consistency display");
	    return;
	}
    case REG_COMPLEMENT:
	{
	    if (!c->do_update)
		update_consistency_display(c->interp, io, c);
	    else
		c->do_update |= REG_LENGTH;
	    return;

	}

    case REG_JOIN_TO:
#ifdef DEBUG
	printf("consistency JOIN_TO contig %d join %d\n", 
	       contig, jdata->join.contig);
#endif
	consistency_join(c->interp, io, c, contig, jdata->join.contig);

	return;

    case REG_LENGTH:
        {
#ifdef DEBUG
	    printf("consistency LENGTH %d\n", contig);
#endif
	    if (!c->do_update)
		update_consistency_display(c->interp, io, c);
	    else
		c->do_update |= REG_LENGTH;
	    return;
	}
    case REG_NUMBER_CHANGE:
#ifdef DEBUG
	printf("consistency NUMBER_CHANGE\n");
#endif
	consistency_renumber(io, c, contig, jdata->number.number);
	return;

     case REG_BUFFER_START:
	{
#ifdef DEBUG
	    printf("consistency REG_BUFFER_START\n");
#endif
	    c->buffer_count++;
	    c->do_update = 1; /* initialise do_update */
	    return;
	}

    case REG_BUFFER_END:
	{
#ifdef DEBUG
	    printf("consistency REG_BUFFER_END count %d \n", c->buffer_count);
#endif
	    c->buffer_count--;

	    if (c->buffer_count <= 0) {
		
		if (c->do_update & REG_LENGTH) {
		    update_consistency_display(c->interp, io, c);
		} else if (c->do_update & TASK_CANVAS_REDRAW) {
		    update_consistency_display(c->interp, io, c);
		} 
		c->buffer_count = 0;
		c->do_update = 0;
	    }
	    return;
	}
    case REG_QUIT:
    case REG_DELETE:
	{
	    /* if deleting a contig shutdown the display */
#ifdef DEBUG
	    printf("consistency QUIT\n");
#endif
	    consistency_shutdown(io, c);
	    return;
	}
    case REG_GET_OPS:
	jdata->get_ops.ops = "Remove\0";
	return;

    case REG_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0:
	    consistency_shutdown(io, c);
	}
	return;

    case REG_PARAMS:
	return;

    case REG_CURSOR_NOTIFY:

#ifdef DEBUG
        printf("consistency REG_CURSOR_NOTIFY\n");
#endif
	consistency_update_cursor(io, c, contig, jdata->cursor_notify.cursor, 1);
	return;
    case REG_GENERIC:
	/* printf("consistency GENERIC %d\n", jdata->generic.task); */
	switch (jdata->generic.task) {
	case TASK_CANVAS_SCROLLX: 
	    {
		char *scroll = (char *)jdata->generic.data;
		consistency_canvasScrollX(c->interp, c, c->win_list, 
					  c->num_wins, scroll);
		break;
	    }
	case TASK_CANVAS_RESIZE:
	  { 
	      char scroll_args[20];
	    /* resize consistency display window */
#ifdef DEBUG
	    printf("consistency TASK_CANVAS_RESIZE\n");
#endif
	    consistency_resizeCanvas(c->interp, c, c->win_list, c->num_wins);
	    sprintf(scroll_args, "scroll 0 units");
	    consistency_canvasScrollX(c->interp, c, c->win_list, 
				      c->num_wins, scroll_args);
	    break;
	  }
	case TASK_CANVAS_WORLD:
	    {
		task_world_t *tw = (task_world_t *)jdata->generic.data;
		int cx = tw->canvasx;
		double wx, wy;
		int win_num;
#ifdef DEBUG
		printf("consistency TASK_CANVAS_WORLD\n");
#endif
		win_num = get_consistency_win_num(c, c->id);
		CanvasToWorld(c->win_list[win_num]->canvas, cx, 0, &wx, &wy);
		tw->basex = wx - c->contig_offset[tw->cnum].offset;
		break;
	    }
	case TASK_CANVAS_ZOOMBACK:
	    {
		reg_generic gen;
#ifdef DEBUG
	        printf("consistency ZOOMBACK\n");
#endif
		consistencyZoomback(c);

		/* redraw ruler ticks */
		Tcl_VarEval(c->interp, c->ruler->window, " delete tick", NULL);
		gen.job = REG_GENERIC;
		gen.task = TASK_DISPLAY_TICKS;
		for (i = 0; i < c->num_wins; i++) {
		  result_notify(io, c->win_list[i]->id, (reg_data *)&gen, 0);
		}

		break;
	    }
	case TASK_CANVAS_ZOOM:
	    {
		s_zoom *szoom = (s_zoom *)jdata->generic.data;
		reg_generic gen;
#ifdef DEBUG
	        printf("consistency ZOOM\n");
#endif
		consistencyZoom(c, szoom->zoom, szoom->scroll, szoom->amount,
				szoom->r_id);

		/* redraw ruler ticks */
		Tcl_VarEval(c->interp, c->ruler->window, " delete tick", NULL);
		gen.job = REG_GENERIC;
		gen.task = TASK_DISPLAY_TICKS;
		for (i = 0; i < c->num_wins; i++) {
		  result_notify(io, c->win_list[i]->id, (reg_data *)&gen, 0);
		}
		break;
	    }
	case TASK_CANVAS_CURSOR_X:
	    {    
		char *label;
		int *cx = (int *)jdata->generic.data;
		double local_pos;
		double wx, wy;
		int win_num;

		/* ruler window */
		win_num = get_consistency_win_num(c, c->id);

		CanvasToWorld(c->win_list[win_num]->canvas, *cx, 0, &wx, &wy);

		label = get_default_string(c->interp, gap_defs,
					   "CONSISTENCY_DISPLAY.CURSOR1");
		canvasCursorX(c->interp, c->win_list[win_num]->canvas, 
			      c->frame, label, 
			      c->xhair.colour, c->xhair.width, *cx, wx, 
			      c->win_list, c->num_wins);

		/* fill in local position of cursor in label box */
		local_pos = TemplateLocalCursor(c->id, c->contig_offset, 
						c->contigs, 
						c->num_contigs, wx);
		
		label = get_default_string(c->interp, gap_defs,
					   "CONSISTENCY_DISPLAY.CURSOR2");
		sprintf(cmd, "%s%s configure -text %d\n", c->frame, label, 
			(int)local_pos);
		Tcl_Eval(c->interp, cmd);
		break;
	    }
	case TASK_CANVAS_CURSOR_DELETE:
	    {    
		int i;
		reg_generic gen;
		gen.job = REG_GENERIC;
		gen.task = TASK_CONS_CURSOR_DELETE;

#ifdef DEBUG
		printf("CONSISTENCY TASK_CANVAS_CURSOR_DELETE\n");
#endif
		for (i = 0; i < c->num_wins; i++) {
		    Tcl_VarEval(c->interp, c->win_list[i]->window, 
				" delete cursor_x", NULL);
		    /* remove y crosshairs if necessary */
		    result_notify(io, c->win_list[i]->id, (reg_data *)&gen, 0);
		}
		break;
	    }
	case TASK_CANVAS_REDRAW:
	    {
		if (!c->do_update) {
		    update_consistency_display(c->interp, io, c);
		} else {
		    c->do_update |= TASK_CANVAS_REDRAW;
		}
		break;
	    }
	case TASK_DISPLAY_RULER:
	    {
#ifdef DEBUG
		printf("TASK_DISPLAY_RULER\n");
#endif
		display_consistency_ruler(io, c->interp, c);
		break;
	    }
	case TASK_DISPLAY_TICKS:
	    {
		/* int *ticks = (int *)jdata->generic.data; */
		char cmd[100];
		int win_num;
		int start, end;
#ifdef DEBUG
		printf("consistency TASK_DISPLAY_TICKS\n");
#endif
		if (!c->configs[CONS_TICKS] || !c->configs[CONS_RULER])
		    return;

		/* get ruler window */
		win_num = get_consistency_win_num(c, c->id);

		for (i = 0; i < c->num_contigs; i++) {
		  if (c->num_contigs > 1) {
		      start = 1;
		      end = ABS(io_clength(io, c->contigs[i]));
		  } else {
		      start = c->start;
		      end = c->end;
		  }
		    display_ruler_ticks(c->interp, 
					c->win_list[win_num]->canvas, 
					c->contig_offset[c->contigs[i]].offset, 
					c->ruler_coord[i].l.y1, 
					c->ruler, start, end);
		}
		scaleSingleCanvas(c->interp, c->win_list[win_num]->world, 
				  c->win_list[win_num]->canvas, 
				  c->ruler->window, 'x', "tick");
		sprintf(cmd, "RulerWindowSize %d %s %s ", 1, c->frame, 
			c->ruler->window);
		Tcl_Eval(c->interp, cmd);
		consistency_update_cursors(io, c, 0);
		break;
	    }
	}
    }
}

/*
 * adding a new consistency window but need to fill in it's zoom list in
 * x direction only
 */
void consistency_update_zoom(obj_consistency_disp *c,
			     int win_num,
			     d_box *orig_zoom)
{
    StackPtr *stackp;

    /* first instance */
    if (win_num == 0) {
        createZoom(&c->win_list[win_num]->zoom);
	pushZoom(&c->win_list[win_num]->zoom, orig_zoom);
	return;
    } 

    /* copy zoom list */
    copyZoom(&c->win_list[win_num]->zoom, c->win_list[0]->zoom);

    /* keep original y values */
    stackp = c->win_list[win_num]->zoom;
    while (stackp != NULL) {
        stackp->data->y1 = orig_zoom->y1;
        stackp->data->y2 = orig_zoom->y2;
        stackp = stackp->next;
    }
}

int add_consistency_window(GapIO *io,
			   obj_consistency_disp *c,
			   char *window,
			   char scroll,
			   int id,
			   double x1, 
			   double y1,
			   double x2,
			   double y2)
{
    int win_num;
    char cmd[1024];
    double fract;
    win **win_list;

    win_num = c->num_wins;
    addWindow(c->win_list, &c->num_wins, window, scroll, id);

    if (NULL == (c->win_list[win_num]->canvas = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    if (NULL == (c->win_list[win_num]->world= (WorldPtr *)xmalloc(sizeof(WorldPtr))))
	return -1;

    if (NULL == (c->win_list[win_num]->world->visible = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    if (NULL == (c->win_list[win_num]->world->total = (d_box *)xmalloc(sizeof(d_box))))
	return -1;
    
    initCanvas(c->interp, c->win_list[win_num]->canvas, 
	       c->win_list[win_num]->window);

    c->win_list[win_num]->world->total->x1 = x1;
    c->win_list[win_num]->world->total->x2 = x2;
    c->win_list[win_num]->world->total->y1 = y1;
    c->win_list[win_num]->world->total->y2 = y2;
    
    /* copy the current zoom for x */
    if (win_num == 0) {
        c->win_list[win_num]->world->visible->x1 = x1;
        c->win_list[win_num]->world->visible->x2 = x2;
    } else {
        c->win_list[win_num]->world->visible->x1 = 
	  c->win_list[0]->world->visible->x1;
        c->win_list[win_num]->world->visible->x2 = 
	  c->win_list[0]->world->visible->x2;
        c->win_list[win_num]->canvas->x =  c->win_list[0]->canvas->x;
    }

    c->win_list[win_num]->world->visible->y1 = y1;
    c->win_list[win_num]->world->visible->y2 = y2;

    SetCanvasCoords(c->interp, 
		    c->win_list[win_num]->world->visible->x1, 
		    c->win_list[win_num]->world->visible->y1, 
		    c->win_list[win_num]->world->visible->x2, 
		    c->win_list[win_num]->world->visible->y2, 
		    c->win_list[win_num]->canvas);

    /* move new window so cursors are in correct place */

    consistency_update_zoom(c, win_num, c->win_list[win_num]->world->visible);
    scaleSingleCanvas(c->interp, c->win_list[win_num]->world, c->win_list[win_num]->canvas, c->win_list[win_num]->window, c->win_list[win_num]->scroll, "all");

    fract = (double)(c->win_list[win_num]->world->visible->x1 - c->win_list[win_num]->world->total->x1) /
	(c->win_list[win_num]->world->total->x2 - c->win_list[win_num]->world->total->x1);
    
    sprintf(cmd, "moveto %f", fract);
    
    if (NULL == (win_list = (win **)xmalloc(sizeof(win*))))
        return -1;

    win_list[0] = c->win_list[win_num];

    consistency_canvasScrollX(c->interp, c, win_list, 1, cmd);
    xfree(win_list);

    Tcl_VarEval(c->interp, c->win_list[win_num]->window, " canvasx 0", NULL);
    c->win_list[win_num]->canvas->x = atoi(Tcl_GetStringResult(c->interp));

    return 0;
}

void delete_consistency_window(obj_consistency_disp *c,
			       int win_num)
{
    int length;
#ifdef DEBUG
    int i;
    for (i = 0; i < c->num_wins; i++) 
      printWorld(c->win_list[i]->world);
#endif

    xfree(c->win_list[win_num]->canvas);
    xfree(c->win_list[win_num]->world->visible);
    xfree(c->win_list[win_num]->world->total);
    xfree(c->win_list[win_num]->world);
    freeZoom(&c->win_list[win_num]->zoom);
    free(c->win_list[win_num]->window);
    xfree(c->win_list[win_num]);

    length = c->num_wins - win_num - 1;
    memmove(&(c->win_list)[win_num], &(c->win_list)[win_num+1], 
	    length * sizeof(win*));

    c->num_wins--;
}

int get_consistency_win_num(obj_consistency_disp *c,
			    int id)
{
    int i;

    for (i = 0; i < c->num_wins; i++) {
        if (c->win_list[i]->id == id)
	    return i;
    }

    return -1;
}


/*
 * draws the consistency display from scratch, ie removes all previous zooming
 * information
 */
int update_consistency_display(Tcl_Interp *interp,
			       GapIO *io,
			       obj_consistency_disp *c)
{
    int ruler_length;
    int i;

/*
    printf("total1 %f %f %f %f\n", c->world->total->x1, 
	   c->world->total->y1, c->world->total->x2, c->world->total->x2);
*/
    consistency_contig_offsets(io, c->contig_offset, c->contigs,
			       c->num_contigs);
    
    ruler_length = c->contig_offset[c->contigs[c->num_contigs-1]].offset + 
	io_clength(io, c->contigs[c->num_contigs-1]);

    /* set start and end */
    c->start = 1;
    c->end = ruler_length;

    /* reset x world for all other win_list in consistency display */
    for (i = 0; i < c->num_wins; i++) {
        c->win_list[i]->world->total->x1 = c->start;
        c->win_list[i]->world->total->x2 = c->end;
    }

#ifdef DEBUG
    printf("RULER length %d %f %f\n", ruler_length, c->cons->world->total->x1,
	   c->cons->world->total->x2);
#endif

    for (i = 0; i < c->num_wins; i++) {
      
        memcpy(c->win_list[i]->world->visible, 
	       c->win_list[i]->world->total, sizeof(d_box));
	SetCanvasCoords(interp, c->win_list[i]->world->visible->x1, 
			c->win_list[i]->world->visible->y1, 
			c->win_list[i]->world->visible->x2, 
			c->win_list[i]->world->visible->y2, 
			c->win_list[i]->canvas); 
	freeZoom(&c->win_list[i]->zoom);
	pushZoom(&c->win_list[i]->zoom, c->win_list[i]->world->visible);
    }

    display_consistency_ruler(io, interp, c);
    consistency_update_cursors(io, c, 0);

    return 0;
}

void consistency_config(Tcl_Interp *interp,
			char *frame,
			int *config_array)
{
    char config[1024];
    int i;

    /* 
     * initialise config_array
     */
    for (i = 0; i < NUM_CONS_CONFIGS; i++) {
	config_array[i] = 0;
    }

    sprintf(config, "config%s.ruler", frame); 
    config_array[CONS_RULER] = atoi(Tcl_GetVar(interp, config, TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[CONS_RULER],
		TCL_LINK_INT);

    sprintf(config, "config%s.ticks", frame); 
    config_array[CONS_TICKS] = atoi(Tcl_GetVar(interp, config, TCL_GLOBAL_ONLY));
    Tcl_LinkVar(interp, config, (char *)&config_array[CONS_TICKS],
		TCL_LINK_INT);
}
 
/*
 * already created a consistency display but redisplaying a ruler
 */
void create_consistency_ruler(GapIO *io,
			      obj_consistency_disp *c)
{
    add_consistency_window(io, c, c->ruler->window, 'x', c->id,
			   c->win_list[0]->world->total->x1, 
			   1, c->win_list[0]->world->total->x2, 1);

    display_consistency_ruler(io, c->interp, c);
}

int consistency_reg(GapIO *io, 
		    Tcl_Interp *interp, 
		    int *contig_array,
		    int num_contigs,
		    int start,
		    int end,
		    char *frame, 
		    ruler_s *ruler, 
		    cursor_s xhair)
{
    obj_consistency_disp *c;
    int id;
    int i;

    if (NULL == (c = (obj_consistency_disp *)xmalloc(sizeof(obj_consistency_disp))))
	return 0;
    
    if (NULL == (c->cursor =
		 (cursor_t **)xmalloc(num_contigs * sizeof(*(c->cursor)))))
	return -1;
    if (NULL == (c->cursor_visible =
		 (int *)xmalloc(num_contigs * sizeof(*(c->cursor_visible)))))
	return -1;

    id = register_id();
    c->num_wins = 0;
      
    strcpy(c->frame, frame);   
    c->contigs = contig_array; 
    c->num_contigs = num_contigs;
    c->id = id;
    c->xhair = xhair;
    c->do_update = 0;
    c->buffer_count = 0;
    c->interp = interp;
    c->ruler = ruler;
    c->start = start;
    c->end = end;
    c->ruler_coord = NULL;

    for (i = 0; i < num_contigs; i++) {
	c->cursor_visible[i] = 0;
	c->cursor[i] = create_contig_cursor(io, contig_array[i], 0, id);
    }

    /* create list of windows in the consistency display */
    if (NULL == (c->win_list = (win **)xmalloc(MAX_NUM_WINS * sizeof(win*))))
	return -1;

    if (NULL == (c->orig_total = (d_box *)xmalloc(sizeof(d_box))))
	return -1;
    if (NULL == (c->contig_offset = (c_offset *)xmalloc((1+io->db.num_contigs)
							* sizeof(c_offset))))
	return -1;
    
    /* initialise */
    for (i = 0; i <= io->db.num_contigs; i++) {
        c->contig_offset[i].gap = 0;
        c->contig_offset[i].offset = 0;
    }
    consistency_config(interp, c->frame, c->configs);

    consistency_contig_offsets(io, c->contig_offset, c->contigs,
			       c->num_contigs);

    /* 
     * initialise the ruler canvas for the basic consistency display - only
     * use the x values
     */
    c->orig_total->x1 = start;
    c->orig_total->x2 = end;
    c->orig_total->y1 = 1;
    c->orig_total->y2 = 1;
    
    /* add ruler window */
    add_consistency_window(io, c, c->ruler->window, 'x', id, start, 1, end, 1);
    
    for (i = 0; i < c->num_wins; i++) {
      
        memcpy(c->win_list[i]->world->total, c->orig_total, sizeof(d_box));
    
        memcpy(c->win_list[i]->world->visible, 
	       c->win_list[i]->world->total, sizeof(d_box));
	SetCanvasCoords(interp, c->win_list[i]->world->visible->x1, 
			c->win_list[i]->world->visible->y1, 
			c->win_list[i]->world->visible->x2, 
			c->win_list[i]->world->visible->y2, 
			c->win_list[i]->canvas); 
	freeZoom(&c->win_list[i]->zoom);
	pushZoom(&c->win_list[i]->zoom, c->win_list[i]->world->visible);
    }
    display_consistency_ruler(io, interp, c);
#ifdef REMOVE
    consistency_update_cursors(io, c, 0);
#endif
    for (i = 0; i < num_contigs; i++) {

	contig_register(io, contig_array[i], consistency_callback, 
			(void *)c, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS | 
			REG_CURSOR_NOTIFY |  REG_NUMBER_CHANGE |
			REG_BUFFER_START | REG_BUFFER_END |
			REG_JOIN_TO | REG_GENERIC, REG_TYPE_CONSISTENCY_DISP);
	consistency_update_cursor(io, c, contig_array[i], c->cursor[i], 1);
    }
    return id;
}
