#include <stdio.h>
#include <limits.h>
#include <string.h>

#include "container.h"
#include "xalloc.h"
#include "tcl_utils.h"
#include "text_output.h"

int add_length_ruler(Tcl_Interp *interp,
		     container *c,
		     int row_index,
		     int column_index,
		     int row_num,
		     int column_num,
		     int orientation)
{
    char cmd[1024];
    int ruler_id;
    char *ruler_win;
    int num = 0;
    char **list;
    element *e_ruler;
    plot_data *output;
    configs *configure;
    d_box bbox;
    seq_id_dir *seq_ids;
    int i, num_seq_id;
    int width, height;

#ifdef DEBUG
    printf("***********add_length_ruler %d\n", orientation);
#endif

    if (orientation == HORIZONTAL) {
	
	row_num = get_default_int(interp, tk_utils_defs, w("CONTAINER.RULER_ROW"));
	height = get_default_int(interp, tk_utils_defs, w("RULER.PLOT_HEIGHT"));
	width = get_default_int(interp, tk_utils_defs, w("ELEMENT.PLOT_WIDTH"));

    } else {
	column_num -= 1;
	width = get_default_int(interp, tk_utils_defs, w("RULER.PLOT_HEIGHT"));
	height = get_default_int(interp, tk_utils_defs, w("ELEMENT.PLOT_HEIGHT"));
    }

    Tcl_ResetResult(interp);
    sprintf(cmd, "create_canvas_ruler %s %d %d %d %d %d %d LENGTH", c->win, c->id, 
	    orientation, row_num, column_num, width, height);

    if (TCL_OK != (Tcl_Eval(interp, cmd))) {
	printf("error create_canvas_ruler: %s\n", Tcl_GetStringResult(interp));
    }
 
    if (Tcl_SplitList(interp, Tcl_GetStringResult(interp), &num, &list) != TCL_OK) {
	return -1;
    }
    
    ruler_id = atoi(list[0]);
    ruler_win = list[1];

    if (NULL == (e_ruler = (get_element(ruler_id))))
	e_ruler = create_element(interp, c->id, ruler_id, ruler_win, orientation, orientation);

    e_ruler->ruler = ruler_struct(interp, tk_utils_defs, "CONTAINER", 0);

   if (NULL == (output = (plot_data *)xmalloc(sizeof(plot_data))))
	return -1;

    if (NULL == (output->configure = (configs **)xmalloc(sizeof(configs*))))
	return -1;

   if (NULL == (configure = (configs *)xmalloc(sizeof(configs))))
	return -1;

    configure->position = 0.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = 1.0;
    configure->zoom = 2;
    configure->scroll = 1;

    output->configure[0] = configure;
    output->n_configure = 1;
    output->sf_m = 1.0;
    output->sf_c = 0.0;
    output->result_id = -1;
    output->hidden = 0;
    output->len_ruler = 0;
    output->amp_ruler = 0;
    output->line_width = 0;
    output->colour = NULL;

    if (orientation == HORIZONTAL) {
	output->scale = SCALE_X;
	get_coord_seq_ids(c, row_index, orientation, &seq_ids, &num_seq_id);
    } else {
	output->scale = SCALE_Y;
	get_coord_seq_ids(c, column_index, orientation, &seq_ids, &num_seq_id);
    }
    
    if (-1 == (add_result_to_element(e_ruler, output, INT_MAX, INT_MAX, 
				     INT_MIN, INT_MIN, 
				     orientation, RULER_LEN))) {
	return -1;
    }

    for (i = 0; i < num_seq_id; i++) {
	add_seq_id_to_element(e_ruler, seq_ids[i].seq_id, orientation);
    }

    initCanvas(interp, e_ruler->pixel, e_ruler->win);
    add_element_to_container(interp, c->id, c->win, e_ruler, 
			     INT_MAX, INT_MIN, INT_MAX, INT_MIN);

    
    if (orientation & HORIZONTAL) {
	c->column[column_index]->ruler = e_ruler;

	e_ruler->world->total->x1 = c->column[column_index]->total.min;
	e_ruler->world->total->x2 = c->column[column_index]->total.max;
	e_ruler->world->total->y1 = 0;
	e_ruler->world->total->y2 = 0;

	e_ruler->ruler->start = e_ruler->world->total->x1;
	e_ruler->ruler->end = e_ruler->world->total->x2;
    } else {
	c->row[row_index]->ruler = e_ruler;

	e_ruler->world->total->x1 = 0;
	e_ruler->world->total->x2 = 0;
	e_ruler->world->total->y1 = c->row[row_index]->total.min;
	e_ruler->world->total->y2 = c->row[row_index]->total.max;

	e_ruler->ruler->start = e_ruler->world->total->y1;
	e_ruler->ruler->end = e_ruler->world->total->y2;
    
    }
    e_ruler->world->visible->x1 = e_ruler->world->total->x1;
    e_ruler->world->visible->x2 = e_ruler->world->total->x2;
    e_ruler->world->visible->y1 = e_ruler->world->total->y1;
    e_ruler->world->visible->y2 = e_ruler->world->total->y2;

    set_pixel_coords(e_ruler->world->visible->x1, e_ruler->world->visible->y1, 
		     e_ruler->world->visible->x2, e_ruler->world->visible->y2,
		     e_ruler->pixel);

    sprintf(e_ruler->ruler->window, "%s", e_ruler->win);

    if (orientation & HORIZONTAL) {
	draw_single_ruler(interp, e_ruler->ruler, e_ruler->pixel, 
			  c->column[column_index]->total.min, 
			  c->column[column_index]->total.max, 1);
    } else {
	draw_single_ruler_vertical(interp, e_ruler->ruler, e_ruler->pixel, 
			  c->row[row_index]->total.min, 
			  c->row[row_index]->total.max, 1);
    }
    bbox.x1 = e_ruler->world->total->x1;
    bbox.x2 = e_ruler->world->total->x2;
    bbox.y1 = e_ruler->world->total->y1;
    bbox.y2 = e_ruler->world->total->y2;

    
    e_ruler->scale_func(interp, e_ruler, -1, &bbox, e_ruler->pixel); 

    e_ruler->scrollregion_func(interp, e_ruler, e_ruler->world->total, 
			       e_ruler->c->column[e_ruler->column_index]->pixel,
			       e_ruler->c->row[e_ruler->row_index]->pixel);
    
    /* remove all current zooming info */
    freeZoom(&e_ruler->zoom);

    /* add first zoom */
    pushZoom(&e_ruler->zoom, e_ruler->world->visible);
    Tcl_Free((char *)list);
    return 0;
}

/* amplitude ruler for element */
int add_element_ruler(Tcl_Interp *interp,
		      container *c,
		      element *e,
		      int orientation)
{
    char cmd[1024];
    plot_data *output;
    configs *configure;
    element *e_ruler;
    int ruler_id;
    int num = 0;
    char **list;
    char *ruler_win;
    d_box bbox;
    int width, height;
    int row_num, column_num;

    row_num = get_element_row(interp, e->win);
    column_num = get_element_column(interp, e->win);

#ifdef DEBUG
    printf("*********** add_element_ruler %s %d r=%d c=%d\n", e->win, e->row_index, row_num, column_num);
#endif

    Tcl_ResetResult(interp);

    if (orientation == HORIZONTAL) {
	height = get_default_int(interp, tk_utils_defs, w("RULER.PLOT_HEIGHT"));
	width = get_default_int(interp, tk_utils_defs, w("ELEMENT.PLOT_HEIGHT"));

	sprintf(cmd, "find_result_position \"\" %s %d", e->win, BOTTOM);
	if (TCL_OK != (Tcl_Eval(interp, cmd))) {
	    verror(ERR_WARN, "add_element_ruler", "create_canvas_ruler: %s\n",
		   Tcl_GetStringResult(interp));
	}

	if (Tcl_SplitList(interp, Tcl_GetStringResult(interp), &num, &list) != TCL_OK) {
	    return -1;
	}
	
	row_num = atoi(list[0]);
	column_num = atoi(list[1]);
	
	sprintf(cmd, "create_canvas_ruler %s %d %d %d %d %d %d AMPLITUDE", 
		c->win, c->id, orientation, row_num, 
		column_num, width, height);
    } else {
	width = get_default_int(interp, tk_utils_defs, w("RULER.PLOT_HEIGHT"));
	height = get_default_int(interp, tk_utils_defs, w("ELEMENT.PLOT_HEIGHT"));
	sprintf(cmd, "create_canvas_ruler %s %d %d %d %d %d %d AMPLITUDE", 
		c->win, c->id, orientation, row_num, column_num - 1, width, 
		height);
    }

    if (TCL_OK != (Tcl_Eval(interp, cmd))) {
	verror(ERR_WARN, "add_element_ruler", "create_canvas_ruler: %s\n",
	       Tcl_GetStringResult(interp));
    }
 
    if (Tcl_SplitList(interp, Tcl_GetStringResult(interp), &num, &list) != TCL_OK) {
	return -1;
    }
    
    ruler_id = atoi(list[0]);
    ruler_win = list[1];

    if (NULL == (e_ruler = (get_element(ruler_id))))
	e_ruler = create_element(interp, c->id, ruler_id, ruler_win, 0, orientation);

    e_ruler->ruler = ruler_struct(interp, tk_utils_defs, "CONTAINER", 0);

   if (NULL == (output = (plot_data *)xmalloc(sizeof(plot_data))))
	return -1;

    if (NULL == (output->configure = (configs **)xmalloc(sizeof(configs*))))
	return -1;

   if (NULL == (configure = (configs *)xmalloc(sizeof(configs))))
	return -1;

    configure->position = 0.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = 1.0;
    configure->zoom = 2;
    configure->scroll = 1;

    output->configure[0] = configure;
    output->n_configure = 1;
    output->sf_m = 1.0;
    output->sf_c = 0.0;
    output->result_id = -1;
    output->hidden = 0;
    output->len_ruler = 0;
    output->amp_ruler = 0;
    output->line_width = 0;
    output->colour = NULL;

    if (orientation == HORIZONTAL) 
	output->scale = SCALE_X;
    else
	output->scale = SCALE_Y;
    
    if (-1 == (add_result_to_element(e_ruler, output, INT_MAX, INT_MAX, 
				     INT_MIN, INT_MIN, 
				     0, RULER_AMP))) {
	return -1;
    }

    initCanvas(interp, e_ruler->pixel, e_ruler->win);
    add_element_to_container(interp, c->id, c->win, e_ruler, 
			     INT_MAX, INT_MIN, INT_MAX, INT_MIN);


    if (orientation == VERTICAL) {
	e_ruler->world->total->x1 = 0;
	e_ruler->world->total->x2 = 0;
	e_ruler->world->total->y1 = e->world->total->y1;
	e_ruler->world->total->y2 = e->world->total->y2;
	e_ruler->world->visible->x1 = 0;
	e_ruler->world->visible->x2 = 0;
	e_ruler->world->visible->y1 = e->world->total->y1;
	e_ruler->world->visible->y2 = e->world->total->y2;
    } else {
	e_ruler->world->total->x1 = e->world->total->y1;
	e_ruler->world->total->x2 = e->world->total->y2;
	e_ruler->world->total->y1 = 0;
	e_ruler->world->total->y2 = 0;
	e_ruler->world->visible->x1 = e->world->total->y1;
	e_ruler->world->visible->x2 = e->world->total->y2;
	e_ruler->world->visible->y1 = 0;
	e_ruler->world->visible->y2 = 0;

    }

    set_pixel_coords(e_ruler->world->visible->x1, e_ruler->world->visible->y1, 
		     e_ruler->world->visible->x2, e_ruler->world->visible->y2,
		     e_ruler->pixel);

    sprintf(e_ruler->ruler->window, "%s", e_ruler->win);

    if (orientation == HORIZONTAL) {
	draw_single_ruler(interp, e_ruler->ruler, e_ruler->pixel, 
			  e->world->total->y1, e->world->total->y2, 1);
    } else {
	draw_single_ruler_vertical(interp, e_ruler->ruler, e_ruler->pixel, 
				   e->world->total->y1, e->world->total->y2, 1);
    }

    bbox.x1 = e_ruler->world->total->x1;
    bbox.x2 = e_ruler->world->total->x2;
    bbox.y1 = e_ruler->world->total->y1;
    bbox.y2 = e_ruler->world->total->y2;

    e_ruler->scale_func(interp, e_ruler, -1, &bbox, e_ruler->pixel); 

    e_ruler->scrollregion_func(interp, e_ruler, e_ruler->world->total, 
			       e_ruler->c->column[e_ruler->column_index]->pixel,
			       e_ruler->c->row[e_ruler->row_index]->pixel);

    /* remove all current zooming info */
    freeZoom(&e_ruler->zoom);

    /* add first zoom */
    pushZoom(&e_ruler->zoom, e_ruler->world->visible);

    e->ruler_id = e_ruler->id;

    e->children_id = e_ruler->id;
    e->children_position = -1; /* left/top of parent */
    e->num_children_left = 1;
    e_ruler->parent = e->id;

    return 0;
}

void update_length_ruler(Tcl_Interp *interp,
			 container *c,
			 coord *coords)
{
    element *e_ruler = NULL;
    d_box bbox;

#ifdef DEBUG
    printf("update_length_ruler\n");
#endif
    e_ruler = coords->ruler;

    if (e_ruler->orientation == HORIZONTAL) {
	draw_single_ruler(interp, e_ruler->ruler, e_ruler->pixel, 
			  coords->total.min, coords->total.max, 1);
    } else {
	draw_single_ruler_vertical(interp, e_ruler->ruler, e_ruler->pixel, 
				   coords->total.min, coords->total.max, 
				   1);
    }

    bbox.x1 = coords->visible.min;
    bbox.x2 = coords->visible.max;
    bbox.y1 = coords->visible.min;
    bbox.y2 = coords->visible.max;
	
    e_ruler->scale_func(interp, e_ruler, -1, &bbox, e_ruler->pixel); 

    e_ruler->scrollregion_func(interp, e_ruler, e_ruler->world->total, 
			       e_ruler->c->column[e_ruler->column_index]->pixel,
			       e_ruler->c->row[e_ruler->row_index]->pixel);
    
    /* remove all current zooming info */
    freeZoom(&e_ruler->zoom);

    /* add first zoom */
    pushZoom(&e_ruler->zoom, e_ruler->world->visible);
}
