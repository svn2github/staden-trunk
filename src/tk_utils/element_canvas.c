#include <tcl.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "container.h"
#include "text_output.h"
#include "element_canvas.h"
#include "tclCanvGraph.h"

double canvas_x(Tcl_Interp *interp,
		char *win,
		double val)
{
    Tcl_Obj *obj[3];
    double wx;
    int i;

    /* find new left hand edge of canvas in canvasx coords */
    obj[0] = Tcl_NewStringObj(win, -1);
    obj[1] = Tcl_NewStringObj("canvasx", -1);
    obj[2] = Tcl_NewDoubleObj(val);
    
    for (i = 0; i < 3; i++) {
	Tcl_IncrRefCount(obj[i]);
    }
    if (Tcl_EvalObjv(interp, 3, obj, 0) != TCL_OK) {
	return -1;
    }
    for (i = 0; i < 3; i++) {
	Tcl_DecrRefCount(obj[i]);
    }
   
    Tcl_GetDoubleFromObj(interp, Tcl_GetObjResult(interp), &wx);
    return wx;
}

double canvas_y(Tcl_Interp *interp,
		char *win,
		double val)
{
    Tcl_Obj *obj[3];
    double wy;
    int i;

    /* find new left hand edge of canvas in canvasx coords */
    obj[0] = Tcl_NewStringObj(win, -1);
    obj[1] = Tcl_NewStringObj("canvasy", -1);
    obj[2] = Tcl_NewDoubleObj(val);
    
    for (i = 0; i < 3; i++) {
	Tcl_IncrRefCount(obj[i]);
    }
    if (Tcl_EvalObjv(interp, 3, obj, 0) != TCL_OK) {
	return -1;
    }
    for (i = 0; i < 3; i++) {
	Tcl_DecrRefCount(obj[i]);
    }
    
    Tcl_GetDoubleFromObj(interp, Tcl_GetObjResult(interp), &wy);
    return wy;
}

int canvas_width(Tcl_Interp *interp,
		 char *win)
{
    Tcl_VarEval(interp, "update idletasks", NULL);
    Tcl_VarEval(interp, "winfo width ", win, NULL);
    return(atoi(Tcl_GetStringResult(interp)));
}

int canvas_height(Tcl_Interp *interp,
		 char *win)
{
    Tcl_VarEval(interp, "update idletasks", NULL);
    Tcl_VarEval(interp, "winfo height ", win, NULL);
    return(atoi(Tcl_GetStringResult(interp)));
}

/* scales single result */
void canvas_scale_result(Tcl_Interp *interp,
			 element *e,
			 int result_id,
			 double x_origin,
			 double sf_x,
			 double y_origin,
			 double sf_y)
{
    char cmd[1024];
    plot_data *result = find_plot_data(e, result_id);

    if (!(result->scale & SCALE_X)) {
	x_origin = 0.0;
	sf_x = 1.0;
    }

    if (!(result->scale & SCALE_Y)) {
	y_origin = 0.0;
	sf_y = 1.0;
    }

    if (result_id != -1) {
	/* must scale cursor aswell */
	sprintf(cmd, "%s scale cursor %.20f %.20f %.20f %.20f \n", 
		e->win, x_origin, y_origin, sf_x, sf_y);
	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    verror(ERR_WARN, "canvas_scale_result", "%s\n", Tcl_GetStringResult(interp));

	/* must scale ruler ticks aswell */
	sprintf(cmd, "%s scale tick %.20f %.20f %.20f %.20f \n", 
		e->win, x_origin, y_origin, sf_x, sf_y);
	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    verror(ERR_WARN, "canvas_scale_result", "%s\n", Tcl_GetStringResult(interp));

	sprintf(cmd, "%s scale id%d %.20f %.20f %.20f %.20f \n", 
		e->win, result->result_id, x_origin, y_origin, sf_x, sf_y);
    } else {
	sprintf(cmd, "%s scale all %.20f %.20f %.20f %.20f \n", 
		e->win, x_origin, y_origin, sf_x, sf_y);	
    }
#ifdef DEBUG
    printf("Scale Canvas Result%s\n", cmd);
#endif

    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	verror(ERR_WARN, "canvas_scale_result", "%s\n", Tcl_GetStringResult(interp));
    }
}

void canvas_move(Tcl_Interp *interp,
		 element *e,
		 int result_id,
		 double x_amount,
		 double y_amount)
{
    char cmd[1024];

    if (result_id == -1) {
	sprintf(cmd, "%s move all %.20f %.20f", e->win, x_amount, y_amount);
    } else {
	sprintf(cmd, "%s move id%d %.20f %.20f", e->win, result_id, 
		x_amount, y_amount);
    }

    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "moveCanvas", "%s\n", Tcl_GetStringResult(interp));

}

/*
 * scales each result in an element
 */
void canvas_scale(Tcl_Interp *interp,
		  element *e,
		  int result_id,
		  d_box *current, /* world coords */
		  CanvasPtr *canvas)
{
    double x_origin, y_origin;
    double sf_x, sf_y;
    int p_x1, p_y1, p_x2, p_y2;
    int i;
    double w_x1, w_x2, w_y1, w_y2;

    /* initialise world coords */
    w_x1 = current->x1;
    w_x2 = current->x2;
    w_y1 = current->y1;
    w_y2 = current->y2;

    p_x1 = e->pixel->x;
    p_x2 = e->pixel->x + e->pixel->width;
    p_y1 = (double)e->pixel->y;
    p_y2 = (double)e->pixel->y + e->pixel->height;

    if (e->orientation & HORIZONTAL) {
	p_x1 = e->c->column[e->column_index]->pixel->x;
	p_x2 = e->c->column[e->column_index]->pixel->x + e->c->column[e->column_index]->pixel->width;
    }

    if (e->orientation & VERTICAL) {
	p_y1 = e->c->row[e->row_index]->pixel->y;
	p_y2 = e->c->row[e->row_index]->pixel->y + e->c->row[e->row_index]->pixel->height;
    }
     
    x_origin = calc_zoom_origin(w_x1, p_x1, w_x2, p_x2);
    sf_x = calc_zoom_sf((double)p_x1, w_x1, p_x2, w_x2);
    y_origin = calc_zoom_origin(w_y1, p_y1, w_y2, p_y2);
    sf_y = calc_zoom_sf(p_y1, w_y1, p_y2, w_y2);
    
    if ((check_element_scale(e) & SCALE_X) && x_origin == 0.0 && sf_x == 1.0) {
	/* translation */
#ifdef TODO
	canvas_move(interp, e, result_id, -1.0 * w_x1, 0.0);
	return;
#endif
    }

    if ((check_element_scale(e) & SCALE_Y) && y_origin == 0.0 && sf_y == 1.0) {
#ifdef TODO
	/* translation */
	canvas_move(interp, e, result_id, 0.0, -1.0 * w_y1);
	return;
#endif
    }

    if (result_id != -1) {
	canvas_scale_result(interp, e, result_id, x_origin, sf_x, y_origin, 
	    sf_y);
	return;
    }

    if (!(check_element_scale(e) & SCALE_X)) {
	x_origin = 0.0;
	sf_x = 1.0;
    }

    if (!(check_element_scale(e) & SCALE_Y)) {
	y_origin = 0.0;
	sf_y = 1.0;
    }

    for (i = 0; i < e->num_results; i++) {
#ifdef DEBUG
	printf("canvas_scale_result %d\n", e->results[i]->result_id);
#endif
	canvas_scale_result(interp, e, e->results[i]->result_id, x_origin, sf_x, y_origin, 
			    sf_y);
    }
}

void canvas_scrollregion(Tcl_Interp *interp,
			 element *e,
			 d_box *total,
			 CanvasPtr *pixel_column,
			 CanvasPtr *pixel_row)
{
    char cmd[1024];
    double w_x1, w_y1, w_x2, w_y2;
    int p_x1, p_x2, p_y1, p_y2;
    int dummy;
    int win_height;

    w_x1 = e->world->total->x1;
    w_x2 = e->world->total->x2;
    w_y1 = e->world->total->y1;
    w_y2 = e->world->total->y2;
    world_to_pixel(e->pixel, w_x1, w_y1, &p_x1, &p_y1);
    world_to_pixel(e->pixel, w_x2, w_y2, &p_x2, &p_y2);

    if (e->orientation & HORIZONTAL) {
	w_x1 = e->c->column[e->column_index]->total.min;
	w_x2 = e->c->column[e->column_index]->total.max;
	world_to_pixel(pixel_column, w_x1, w_y1, &p_x1, &dummy);
	world_to_pixel(pixel_column, w_x2, w_y2, &p_x2, &dummy);
    } 

    if (e->orientation & VERTICAL) {
	int ht;

	w_y1 = e->c->row[e->row_index]->total.min;
	w_y2 = e->c->row[e->row_index]->total.max;
	world_to_pixel(pixel_row, w_x1, w_y1, &dummy, &p_y1);
	world_to_pixel(pixel_row, w_x2, w_y2, &dummy, &p_y2);

	/* need to turn scrollregion round */
	win_height = e->element_height(interp, e->win);
	ht = p_y2 - p_y1;
    }

    /* no x scale */
    if (!(check_element_scale(e) & SCALE_X)) {
	p_x1 = 0.0;
	p_x2 = 0.0;
    }

    /* no y scale */
    if (!(check_element_scale(e) & SCALE_Y)) {
	p_y1 = 0.0;
	p_y2 = 0.0;
    }

    sprintf(cmd, "%s configure -scrollregion \"%d %d %d %d\"", 
	    e->win, p_x1, p_y1, p_x2, p_y2);
#ifdef DEBUG
    printf("scroll region %s\n", cmd);
#endif
    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "scrollRegion", "%s\n", Tcl_GetStringResult(interp));
}

void canvas_scroll_x(Tcl_Interp *interp,
		     element *e,
		     char *command)
{
    char cmd[1024];
    double wy;

    sprintf(cmd, "%s xview %s", e->win, command);
#ifdef DEBUG	
    printf("canvas_scroll_x %s\n", cmd);
    fflush(stdout);
#endif
    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "canvas_scroll_x", "%s\n", Tcl_GetStringResult(interp));

    Tcl_VarEval(interp, "update idletasks", NULL);

#ifdef DEBUG
    printf("SCROLL X e %f %f c %f %f\n", e->world->visible->x1,
	   e->world->visible->x2, e->c->column[e->column_index]->visible.min,
	   e->c->column[e->column_index]->visible.max);
#endif

    /* find new left hand edge of canvas in canvasx coords */
    e->pixel->x = e->element_x(interp, e->win, 0.0);

    /* find new left and right edges of canvas in world coords */
    CanvasToWorld(e->pixel, e->pixel->x, 0, &e->world->visible->x1, &wy);
    CanvasToWorld(e->pixel, e->pixel->x + e->pixel->width,
		  0, &e->world->visible->x2, &wy);

    set_pixel_coords(e->world->visible->x1, e->world->visible->y1, 
		     e->world->visible->x2, e->world->visible->y2, e->pixel);
}

void canvas_scroll_y(Tcl_Interp *interp,
		     element *e,
		     char *command)
{
    char cmd[1024];
    double wx;
    int i, j, k;
    double y1;
    Tcl_Obj *get_coords[3];
    Tcl_Obj *set_coords[5];
    double coords[4];
    Tcl_Obj *cobj;
    int argc;
    Tcl_Obj ** argv;

    sprintf(cmd, "%s yview %s", e->win, command);
    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "canvas_scroll_y", "%s\n", Tcl_GetStringResult(interp));

    for (i = 0; i < e->num_results; i++) {
	for (j = 0; j < e->results[i]->n_configure; j++) {
	    if (!e->results[i]->configure[j]->scroll && e->orientation == HORIZONTAL) {
		sprintf(cmd, "id%d", e->results[i]->result_id);

		get_coords[0] = Tcl_NewStringObj(e->win, -1);
		get_coords[1] = Tcl_NewStringObj("coords", -1);
		get_coords[2] = Tcl_NewStringObj(cmd, -1);

		for (k = 0; k < 3; k++) {
		    Tcl_IncrRefCount(get_coords[k]);
		}
		if (Tcl_EvalObjv(interp, 3, get_coords, 0) != TCL_OK) {
		    printf("Failed get_coords\n");
		    return;
		}
		for (k = 0; k < 3; k++) {
		    Tcl_DecrRefCount(get_coords[k]);
		}
		cobj = Tcl_GetObjResult(interp);
		Tcl_IncrRefCount(cobj);

		Tcl_ListObjGetElements(interp, cobj, &argc, (Tcl_Obj ***) &argv);
		for (k = 0; k < 4; k++) {
		    Tcl_GetDoubleFromObj(interp, argv[k], (double *)&coords[k]);
		}
		Tcl_DecrRefCount(cobj);
		y1 = e->element_y(interp, e->win, coords[1]);

		set_coords[0] = Tcl_NewStringObj(e->win, -1);
		set_coords[1] = Tcl_NewStringObj("coords", -1);
		set_coords[2] = Tcl_NewStringObj(cmd, -1);
		set_coords[3] = Tcl_NewDoubleObj(coords[0]);
		set_coords[4] = Tcl_NewDoubleObj(y1);
		for (k = 0; k < 5; k++) {
		    Tcl_IncrRefCount(set_coords[k]);
		}
		if (Tcl_EvalObjv(interp, 5, set_coords, 0) != TCL_OK) {
		    printf("Failed set_coords\n");
		    return;
		}
		for (k = 0; k < 5; k++) {
		    Tcl_DecrRefCount(set_coords[k]);
		}
	    }
	}
    }

    Tcl_VarEval(interp, "update idletasks", NULL);

    /* find new left hand edge of canvas in canvasy coords */
    e->pixel->y = e->element_y(interp, e->win, 0.0);
    
    /* find new top and bottom edges of canvas in world coords */
    CanvasToWorld(e->pixel, 0, e->pixel->y, &wx, &e->world->visible->y1);
    CanvasToWorld(e->pixel, 0, e->pixel->y + e->pixel->height, &wx, 
		  &e->world->visible->y2);
    set_pixel_coords(e->world->visible->x1, e->world->visible->y1, 
		     e->world->visible->x2, e->world->visible->y2, e->pixel);
}

int create_graph(Tcl_Interp *interp,
		 char *e_win,
		 Tcl_Obj *results,
		 int width,
		 char *foreground,
		 char *tags,
		 int orientation)
{
    Tcl_Obj *create_graph[21];
    int i;
    char orient[2];
    /* canvas create graph 0 0 -graph -width -foreground -anchor nw */

    if (orientation & HORIZONTAL) {
	strcpy(orient, "h");
    } else {
	strcpy(orient, "v");
    }

    create_graph[0] = Tcl_NewStringObj(e_win, -1);
    create_graph[1] = Tcl_NewStringObj("create", -1);
    create_graph[2] = Tcl_NewStringObj("graph", -1);
    create_graph[3] = Tcl_NewIntObj(0);
    create_graph[4] = Tcl_NewIntObj(0);
    create_graph[5] = Tcl_NewStringObj("-anchor", -1);
    create_graph[6] = Tcl_NewStringObj("nw", -1);
    create_graph[7] = Tcl_NewStringObj("-graph", -1);
    create_graph[8] = results; 
    create_graph[9] = Tcl_NewStringObj("-width", -1);
    create_graph[10] = Tcl_NewIntObj(width);
    create_graph[11] = Tcl_NewStringObj("-fill", -1);
    create_graph[12] = Tcl_NewStringObj(foreground, -1);
    create_graph[13] = Tcl_NewStringObj("-tags", -1);
    create_graph[14] = Tcl_NewStringObj(tags, -1);
    create_graph[15] = Tcl_NewStringObj("-invertx", -1);
    create_graph[16] = Tcl_NewBooleanObj(0);
    create_graph[17] = Tcl_NewStringObj("-inverty", -1);
    create_graph[18] = Tcl_NewBooleanObj(1);
    create_graph[19] = Tcl_NewStringObj("-orient", -1);
    create_graph[20] = Tcl_NewStringObj(orient, -1);

    
    for (i = 0; i < 21; i++) {
	Tcl_IncrRefCount(create_graph[i]);
    }
    
    if (Tcl_EvalObjv(interp, 21, create_graph, 0) != TCL_OK) {
	printf("Failed create graph\n");
	return -1;
    }
    
    for (i = 0; i < 21; i++) {
	Tcl_DecrRefCount(create_graph[i]);
    }
    
    return 0;
}

void create_canv_item(Tcl_Interp *interp,
		      char *e_win,
		      Graph *graph,
		      plot_data *result,
		      int width,
		      char *foreground,
		      char *tags,
		      int orientation)
{
    int i, j;

    for (i = 0; i < graph->n_darrays; i++) {
	if (graph->d_arrays[i].type == G_LINES) {
	    /* discrete lines */
	    for (j = 0; j < graph->d_arrays[i].n_dlines; j++) {

	    }
	}

    }

}

void create_canv_dot(Tcl_Interp *interp,
		      char *e_win,
		      Graph *graph,
		      plot_data *result,
		      int width,
		      char *foreground,
		      char *tags,
		      int orientation)
{
    int i, j;
    char cmd[1024];

    for (i = 0; i < graph->n_parrays; i++) {
	for (j = 0; j < graph->p_arrays[i].n_pts; j++) {
	    if (result->configure[i]->position != -1.0) {

		if (orientation & HORIZONTAL) {
		    if (result->configure[i]->y_direction == '+') {
			sprintf(cmd,"%s create line %.20f %.20f %.20f %.20f -width %d -fill %s "
				"-tag {%s S d%d.%d}",
				e_win, graph->p_arrays[i].p_array[j].x,
				graph->dim.y1 - graph->p_arrays[i].p_array[j].y + graph->dim.y0, 
				graph->p_arrays[i].p_array[j].x,
				graph->dim.y1 - graph->p_arrays[i].p_array[j].y + graph->dim.y0,
				result->line_width, result->colour, tags, i, j);
		    } else {
			sprintf(cmd,"%s create line %.20f %.20f %.20f %.20f -width %d -fill %s "
				"-tag {%s S d%d.%d}",
				e_win, graph->p_arrays[i].p_array[j].x,
				graph->p_arrays[i].p_array[j].y, 
				graph->p_arrays[i].p_array[j].x,
				graph->p_arrays[i].p_array[j].y,
				result->line_width, result->colour, tags, i, j);

		    }
#ifndef DEBUG
		    printf("canv_dot %s\n", cmd);
#endif
		    Tcl_Eval(interp, cmd);
		}
		if (orientation & VERTICAL) {
		    if (result->configure[i]->x_direction == '+') {
			sprintf(cmd,"%s create line %.20f %.20f %.20f %.20f -width %d -fill %s "
				"-tag {%s S d%d.%d}",
				e_win, 
				graph->p_arrays[i].p_array[j].y,
				graph->dim.x1 - graph->p_arrays[i].p_array[j].x + graph->dim.x0,
				graph->p_arrays[i].p_array[j].y,
				graph->dim.x1 - graph->p_arrays[i].p_array[j].x + graph->dim.x0,
				result->line_width, result->colour, tags, i, j);
		    } else {
			sprintf(cmd,"%s create line %.20f %.20f %.20f %.20f -width %d -fill %s "
				"-tag {%s S d%d.%d}",
				e_win, 
				graph->p_arrays[i].p_array[j].y, 
				graph->p_arrays[i].p_array[j].x,
				graph->p_arrays[i].p_array[j].y,
				graph->p_arrays[i].p_array[j].x,
				result->line_width, result->colour, tags, i, j);

		    }
#ifndef DEBUG
		    printf("canv_dot %s\n", cmd);
#endif
		    Tcl_Eval(interp, cmd);
 
		}
	    }
	}
    }
}

void create_canv_line(Tcl_Interp *interp,
		      char *e_win,
		      Graph *graph,
		      plot_data *result,
		      int width,
		      char *foreground,
		      char *tags,
		      int orientation)
{
    int i, j;
    char cmd[1024];

    for (i = 0; i < graph->n_darrays; i++) {
	for (j = 0; j < graph->d_arrays[i].n_dlines; j++) {
	    if (result->configure[i]->position != -1.0) {

		if (orientation & HORIZONTAL) {
		    if (result->configure[i]->y_direction == '+') {
			sprintf(cmd,"%s create line %.20f %.20f %.20f %.20f -width %d -fill %s "
				"-tag {%s S d%d.%d}",
				e_win, graph->d_arrays[i].d_array[j].x0,
				graph->dim.y1 - graph->d_arrays[i].d_array[j].y0 + graph->dim.y0, 
				graph->d_arrays[i].d_array[j].x1,
				graph->dim.y1 - graph->d_arrays[i].d_array[j].y1 + graph->dim.y0,
				result->line_width, result->colour, tags, i, j);
		    } else {
			sprintf(cmd,"%s create line %.20f %.20f %.20f %.20f -width %d -fill %s "
				"-tag {%s S d%d.%d}",
				e_win, graph->d_arrays[i].d_array[j].x0,
				graph->d_arrays[i].d_array[j].y0, 
				graph->d_arrays[i].d_array[j].x1,
				graph->d_arrays[i].d_array[j].y1,
				result->line_width, result->colour, tags, i, j);

		    }
#ifdef DEBUG
		    printf("canv_line %s\n", cmd);
#endif
		    Tcl_Eval(interp, cmd);
		}
		if (orientation & VERTICAL) {
		    if (result->configure[i]->x_direction == '+') {
			sprintf(cmd,"%s create line %.20f %.20f %.20f %.20f -width %d -fill %s "
				"-tag {%s S d%d.%d}",
				e_win, 
				graph->d_arrays[i].d_array[j].y0,
				graph->dim.x1 - graph->d_arrays[i].d_array[j].x0 + graph->dim.x0,
				graph->d_arrays[i].d_array[j].y1,
				graph->dim.x1 - graph->d_arrays[i].d_array[j].x1 + graph->dim.x0,
				result->line_width, result->colour, tags, i, j);
		    } else {
			sprintf(cmd,"%s create line %.20f %.20f %.20f %.20f -width %d -fill %s "
				"-tag {%s S d%d.%d}",
				e_win, 
				graph->d_arrays[i].d_array[j].y0, 
				graph->d_arrays[i].d_array[j].x0,
				graph->d_arrays[i].d_array[j].y1,
				graph->d_arrays[i].d_array[j].x1,
				result->line_width, result->colour, tags, i, j);

		    }
#ifdef DEBUG
		    printf("canv_line %s\n", cmd);
#endif
		    Tcl_Eval(interp, cmd);
 
		}
	    }
	}
    }
}

void draw_canvas_crosshair(Tcl_Interp *interp,
			   element *e,
			   int pos,
			   int orientation)
{
    char cmd[1024];
    int cx, cy;
    double wx, wy, dummy;
    coord *column;

    if (orientation == HORIZONTAL) {
	sprintf(cmd, "%s canvasx %d", e->win, pos);
	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    verror(ERR_WARN, "draw_canvas_crosshairX", "%s\n", Tcl_GetStringResult(interp));

	cx = atoi(Tcl_GetStringResult(interp));

	column = e->c->column[e->column_index];

	pixel_to_world(column->pixel, cx, 0, &wx, &dummy);

	sprintf(cmd, "draw_canvas_crosshairX %s %s %d %.20f\n", e->c->win, e->win, cx, wx);

	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    verror(ERR_WARN, "draw_canvas_crosshairX", "%s\n", Tcl_GetStringResult(interp));
    }

    if (orientation == VERTICAL) {
	sprintf(cmd, "%s canvasy %d", e->win, pos);
	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    verror(ERR_WARN, "draw_canvas_crosshairY1", "%s\n", Tcl_GetStringResult(interp));
	cy = atoi(Tcl_GetStringResult(interp));
	pixel_to_world(e->pixel, 0, cy, &dummy, &wy);	
	sprintf(cmd, "draw_canvas_crosshairY %s %s %d %.20f\n", e->c->win, e->win, cy, invert_wy(e, wy));
	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    verror(ERR_WARN, "draw_canvas_crossshairY2", "%s\n", Tcl_GetStringResult(interp));
    }
}

void delete_canvas_crosshair(Tcl_Interp *interp,
			     element *e)
{
    Tcl_VarEval(interp, e->win, " delete crosshair_x", NULL);
    Tcl_VarEval(interp, e->win, " delete crosshair_y", NULL);
}

