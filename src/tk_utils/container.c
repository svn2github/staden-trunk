#include <tcl.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include "container.h"
#include "xalloc.h"
#include "tcl_utils.h"
#include "canvas_box.h"
#include "text_output.h"
#include "element_canvas.h"
#include "container_ruler.h"

static container **container_array; /* array */
static int container_num = 0;  /* current element of array */
static int container_size = 0; /* current size of array */

void delete_element(element *e, int job);
void container_update_scrollregion(Tcl_Interp *interp, container *c);

/*
 * delete entire container array
 */
void shutdown_container(void)
{

    xfree(container_array);
    container_num = 0;
    container_size = 0;
}

int container_id_to_num(int container_id)
{
    int i;

    for (i = 0; i < container_num; i++) {
	if (container_array[i]->id == container_id) {
	    return i;
	}
    }

    return -1;
}

/*
 * shutdown single container
 */
void delete_container(container *c) {
    int c_num;
    char cmd[1024];

#ifdef DEBUG
    printf("DELETE CONTAINER\n");
#endif
    if (-1 == (c_num = container_id_to_num(c->id))) {
	/* already deleted */
	return;
    }

    sprintf(cmd, "destroy %s", c->win);
    Tcl_Eval(c->interp, cmd);

    c->num_rows = 0;
    c->num_columns = 0;

    /* only move if not the last container in array */
    if (c_num < container_num - 1) {

	memmove(&container_array[c_num], &container_array[c_num+1],
		sizeof(container *));
    }
    if (container_num > 0)
	container_num--;
}

void container_start_shutdown(int container_id) 
{
    container *c = get_container(container_id);
    element *e;
    int i, j, cnt;
    int num_rows = c->num_rows;
    int num_columns = c->num_columns;
    int *e_id;

    /* need to set num_rows and num_columns as the shutdown func will alter
     * these
     */
    if (NULL == (e_id = (int *)xmalloc((c->num_rows * c->num_columns) * 
				       sizeof(int))))
	return;

    cnt = 0;
    for (i = 0; i < num_rows; i++) {
	for (j = 0; j < num_columns; j++) {
	    e = c->matrix[i][j];
	    if (e) {
		e_id[cnt++] = e->id;
	    }
	}
    }

    c->status = 1;

    for (i = 0; i < cnt; i++) {
	e = get_element(e_id[i]);
	if (e) {
	    e->shutdown_func(e, ALL);
	}
    }
    xfree(e_id);
    delete_container(c);

}

/* allocates a unique identifier and window name to container */
int new_container(Tcl_Interp *interp,
		  char **c_win)
{
    static int container_id = 0;
    char *win;

    win = get_default_string(interp, tk_utils_defs, w("CONTAINER.WIN"));
    
    if (NULL == (*c_win = (char *)xmalloc((strlen(win) + 10) * 
						 sizeof(char))))
	return -1;
    sprintf(*c_win, "%s%d", win, container_id);
    
    return container_id++;
}

/*
 * creates container_array if necessary
 * adds new container to container_array
 */
container *create_container(Tcl_Interp *interp,
			    char *c_win,
			    int container_id)
{
    int incr_num = 10;
    int i;
    ruler_s *ruler;
    tick_s *tick;


    if (container_size == 0) {
	container_size += incr_num;
	
	if (NULL == (container_array = (container **)xmalloc(container_size * 
							    sizeof(container*)))){
	    shutdown_container();
	    return NULL;
	}
	for (i = container_num; i < container_size; i++) {
	    if (NULL == (container_array[i] = (container *)xmalloc(sizeof(container)))) {
		return NULL;
	    }
	}
    }
    
    /* DEBUG - need to check this! */
    if (container_num >= container_size) {
	container_size += incr_num;
	if (NULL == (container_array = (container **)xrealloc(container_array, 
							     (container_size)  * 
							     sizeof(container*)))) {
	    shutdown_container();
	    return NULL;
	}
	
	for (i = container_num; i < container_size; i++) {
	    if (NULL == (container_array[i] = (container *)xmalloc(sizeof(container)))) {
		return NULL;
	    }
	}
    }
    

    /* intialise container */
    container_array[container_num]->interp = interp;
    container_array[container_num]->win = strdup(c_win);
    container_array[container_num]->id = container_id;
    container_array[container_num]->matrix = NULL;
    container_array[container_num]->row = NULL;
    container_array[container_num]->column = NULL;
    container_array[container_num]->num_rows = 0;
    container_array[container_num]->max_rows = 0;
    container_array[container_num]->num_columns = 0;
    container_array[container_num]->max_columns = 0;
    container_array[container_num]->status = 0;

    ruler = ruler_struct(interp, tk_utils_defs, "CONTAINER", 0);
    tick = tick_struct(interp, tk_utils_defs, "CONTAINER", -1, -1, "");

    container_array[container_num]->ruler = ruler;
    container_array[container_num]->tick = tick;

    return container_array[container_num++];
}

container *get_container(int container_id)
{
    int i;
    
    for (i = 0; i < container_num; i++) {
	if (container_array[i]->id == container_id)
	    return container_array[i];
    }
    return NULL;
}

/*
 * given an array of seqs, find a container which contains one of the 
 * "direction" sequences in the same direction
 */
int find_container(seq_id_dir *seq_ids,
		   int num_seqs,
		   int *direction,
		   char **e_win,
		   char **c_win)
{
    container *c;
    element *e;
    int i, j, k, l, m;

    for (i = 0; i < container_num; i++) {
	c = container_array[i];

	for (j = 0; j < c->num_rows; j++) {
	    for (k = 0; k < c->num_columns; k++) {
		e = c->matrix[j][k];
		if (e) {
		    for (l = 0; l < num_seqs; l++) {
			for (m = 0; m < e->num_seqs; m++) {
#ifdef DEBUG
			    printf("!!!!!!seq_ids %d %d %d %d %d\n", *direction,
				   seq_ids[l].seq_id,
				   seq_ids[l].direction, 
				   e->seqs[m].seq_id,
				   e->seqs[m].direction);
#endif
			    if (seq_ids[l].seq_id == e->seqs[m].seq_id) {
				*c_win = c->win;
				*direction = e->seqs[m].direction;
				*e_win = e->win;
				return container_array[i]->id;
			    }
			}
		    }
		}
	    }
	}
    }   
    return (-1);
}

int init_row(coord *row)
{
    if (NULL == (row->pixel = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    row->pixel->width = 0;
    row->pixel->height = 0;
    
    row->pixel->x = 0;
    row->pixel->y = 0;
    row->pixel->ax = 0;
    row->pixel->ay = 0;
    row->pixel->bx = 0;
    row->pixel->by = 0;

    row->visible.min = INT_MAX;
    row->visible.max = INT_MIN;
    row->total.min = INT_MAX;
    row->total.max = INT_MIN;
    row->ruler = NULL;
    row->ruler_reg = 0;
    createZoom(&row->zoom);
    return 0;
}

int init_column(coord *column)
{
    if (NULL == (column->pixel = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    column->pixel->width = 0;
    column->pixel->height = 0;
    
    column->pixel->x = 0;
    column->pixel->y = 0;
    column->pixel->ax = 0;
    column->pixel->ay = 0;
    column->pixel->bx = 0;
    column->pixel->by = 0;
   
    column->visible.min = INT_MAX;
    column->visible.max = INT_MIN;
    column->total.min = INT_MAX;
    column->total.max = INT_MIN;
    column->ruler = NULL;
    column->ruler_reg = 0;
    createZoom(&column->zoom);
    return 0;
}

int init_container_matrix(container *c,
			  int row_num,
			  int column_num,
			  int *r_index,
			  int *c_index)
{
    int incr = 10;
    int i, j;

    c->max_rows += incr;
    c->max_columns += incr;

    if (NULL == (c->matrix = (element ***)xmalloc(c->max_rows * 
						  sizeof(element **))))
	return -1;

    for (i = 0; i < c->max_rows; i++) {
	if (NULL == (c->matrix[i] = (element **)xmalloc(c->max_columns * 
						       sizeof(element*))))
	    return -1;
    }

    for (i = 0; i < c->max_rows; i++) {
	for (j = 0; j < c->max_columns; j++) {
	    c->matrix[i][j] = NULL;
	}
    }

    if (NULL == (c->row = (coord **)xmalloc(c->max_rows * sizeof(coord*))))
	return -1;

    if (NULL == (c->column = (coord **)xmalloc(c->max_columns * sizeof(coord*))))
	return -1;
    
    for (i = 0; i < c->max_rows; i++) {
	if (NULL == (c->row[i] = (coord *)malloc(sizeof(coord))))
	    return -1;
	init_row(c->row[i]);
    }

    for (i = 0; i < c->max_columns; i++) {
	if (NULL == (c->column[i] = (coord *)malloc(sizeof(coord))))
	    return -1;
	init_column(c->column[i]);
    }

    c->num_rows++;
    c->num_columns++;

    *r_index = 0;
    *c_index = 0;

    /* row number */
    return 0;
}

void init_pixel(Tcl_Interp *interp,
		element *e,
		CanvasPtr *canvas)        /* in, out */
{
    /* assume if one isn't defined, none of them are */
    if (!e->element_width)
	return;

    canvas->width = e->element_width(interp, e->win);
    canvas->height = e->element_height(interp, e->win);
    
    canvas->x = 0;
    canvas->y = 0;
    canvas->ax = 0;
    canvas->ay = 0;
    canvas->bx = 0;
    canvas->by = 0;
}


/* 
 * go through all results in element to accumulate scaling for element
 */
int check_element_scale(element *e)
{
    int i;
    int scale = 0;

    for (i = 0; i < e->num_results; i++) {
#ifdef DEBUG
	printf("RESULTS SCALE %s %d\n", e->win, e->results[i]->scale);
#endif
	if (e->results[i]->scale & SCALE_X) {
	    scale |= SCALE_X;
	}
	if (e->results[i]->scale & SCALE_Y) {
	    scale |= SCALE_Y;
	}
    }
    return scale;
}

/* 
 * go through all results in element to check for length ruler
 */
int check_element_len_ruler(element *e)
{
    int i;

    for (i = 0; i < e->num_results; i++) {
	if (e->results[i]->len_ruler) {
	    return 1;
	}
    }
    return 0;
}

/* 
 * go through all results in element to check for amplitude ruler
 */
int check_element_amp_ruler(element *e)
{
    int i;

    for (i = 0; i < e->num_results; i++) {
	if (e->results[i]->amp_ruler) {
	    return 1;
	}
    }
    return 0;
}

int get_element_row(Tcl_Interp *interp,
		    char *win)
{
    char cmd[1024];

    sprintf(cmd, "get_element_row %s", win);
    Tcl_Eval(interp, cmd);
    return (atoi(Tcl_GetStringResult(interp)));
    
}
int get_element_column(Tcl_Interp *interp,
		    char *win)
{
    char cmd[1024];

    sprintf(cmd, "get_element_column %s", win);
    Tcl_Eval(interp, cmd);
    return (atoi(Tcl_GetStringResult(interp)));
    
}

int alloc_more_rows(container *c)
{
    int incr = 10;
    int prev_max = c->max_rows;
    int i, j;

    
#ifdef DEBUG
    printf("alloc_more_rows num %d max %d\n", c->num_rows, c->max_rows);
#endif

    if (c->num_rows >= c->max_rows) {

	c->max_rows += incr;
	if (NULL == (c->matrix = (element ***)xrealloc(c->matrix, c->max_rows * 
						     sizeof(element **))))
	    return -1;
	
	if (NULL == (c->row = (coord **)xrealloc(c->row, c->max_rows * sizeof(coord*))))
	    return -1;
	/* initialiase rows */
	for (i = prev_max; i < c->max_rows; i++) {
	    if (NULL == (c->matrix[i] = (element **)xmalloc(c->max_columns * sizeof(element*))))
		return -1;

	    if (NULL == (c->row[i] = (coord *)xmalloc(sizeof(coord))))
	    init_row(c->row[i]);
	}

	/* initialise new elements */
	for (i = prev_max; i < c->max_rows; i++) {
	    for (j = 0; j < c->max_columns; j++) {
		c->matrix[i][j] = NULL;
	    }
	}

	if (c->max_columns == 0) {
	    c->max_columns++;
	    c->num_columns++;
	} 
    }
    return 0;
}

/* adds an entire row */
int add_row_to_container(container *c,
			 int row_index,
			 int column_index,
			 int row_num,
			 int min,
			 int max)
{
    int i, j;

#ifdef DEBUG
    printf("add_row_to_container num %d max %d\n", c->num_rows, c->max_rows);
#endif
    
    alloc_more_rows(c);

    /* 
     * need to modify the row_index for all those elements which will be moved
     * by the subsequent memmove
     */
    for (i = row_index; i < c->num_rows; i++) {
	for (j = column_index; j < c->num_columns; j++) {
	    if (c->matrix[i][j])
		c->matrix[i][j]->row_index++;
	}
    }

    /* only move if not the last row in array */
    if (row_index < c->num_rows) {
	memmove(&c->row[row_index + 1],
		&c->row[row_index],
		(c->num_rows - row_index) * sizeof(coord*));
	memmove(&c->matrix[row_index + 1],
		&c->matrix[row_index],
		(c->num_rows - row_index) * sizeof(element*));

    }

    /* allocate new row */
    if (NULL == (c->row[row_index] = (coord *)malloc(sizeof(coord))))
	return -1;
    init_row(c->row[row_index]);
    
    if (NULL == (c->matrix[row_index] = (element **)malloc(c->max_columns * 
							   sizeof(element*))))
	return -1;
    for (i = 0; i < c->max_columns; i++) {
	c->matrix[row_index][i] = NULL;
    }

    c->num_rows++; 

    return 0;
}

int alloc_more_columns(container *c)
{
    int incr = 3;
    int prev_max = c->max_columns;
    int i, j;

    if (c->num_columns >= c->max_columns) {

	c->max_columns += incr;

	if (c->max_rows == 0) {
	    c->max_rows++;
	    c->num_rows++;
	    if (NULL == (c->matrix = (element ***)xrealloc(c->matrix, c->max_rows * sizeof(element **))))
		return -1;
	}

	if (NULL == (c->column = (coord **)xrealloc(c->column, c->max_columns * sizeof(coord*))))
	    return -1;
	
	/* initialise columns */
	for (i = prev_max; i < c->max_columns; i++) {
	    if (NULL == (c->column[i] = (coord *)xmalloc(sizeof(coord))))
		return -1;
	    init_column(c->column[i]);
	}

	for (i = 0; i < c->max_rows; i++) {
	    if (NULL == (c->matrix[i] = (element **)xrealloc(c->matrix[i], 
					 c->max_columns * sizeof(element*))))
		return -1;
	    for (j = prev_max; j < c->max_columns; j++) {
		c->matrix[i][j] = NULL;
	    }
	}
    }
    return 0;
}

int add_column_to_container(container *c,
			    int row_index,
			    int column_index,
			    int column_num,
			    int min,
			    int max)
{
    int i, j;

#ifdef DEBUG
    printf("add_column_to_container\n");
#endif

    alloc_more_columns(c);

    /* 
     * need to modify the row_index for all those elements which will be moved
     * by the subsequent memmove
     */
    for (i = row_index; i < c->num_rows; i++) {
	for (j = column_index; j < c->num_columns; j++) {
	    c->matrix[i][j]->column_index++;
	}
    }

    /* only move if not the last column in array */
    if (column_index < c->num_columns) {
	memmove(&c->column[column_index + 1],
		&c->column[column_index],
		(c->num_columns - column_index) * sizeof(coord*));

	for (i = 0; i < c->num_rows; i++) {
	    memmove(&c->matrix[i][column_index + 1],
		    &c->matrix[i][column_index],
		    (c->num_columns - column_index) * sizeof(element));

	}
    }
    /* allocate new column */

    if (NULL == (c->column[column_index] = (coord *)malloc(sizeof(coord))))
	return -1;
    init_column(c->column[column_index]);
    
    for (i = 0; i < c->num_rows; i++) {
#ifdef REMOVE
	if (NULL == (c->matrix[i][column_index] = (element *)malloc(sizeof(element))))
	    return -1;
#endif
	c->matrix[i][column_index] = NULL;
    }

    c->num_columns++; 
    return 0;
}

/* removes entire row */
void delete_row_from_container(container *c,
			       int row_index,
			       int column_index)
{
    int i, j;
    element *e;

#ifdef DEBUG
    printf("delete row %d %d\n", row_index, c->num_rows);
#endif

    /* 
     * need to modify the row_index for all those elements which will be moved
     * by the subsequent memmove
     */
    for (i = row_index; i < c->num_rows; i++) {
	for (j = column_index; j < c->num_columns; j++) {
	    e = c->matrix[i][j]; 
	    if (e) {
		e->row_index--;
	    }
	}
    }

    xfree(c->row[row_index]->pixel);
    freeZoom(&c->row[row_index]->zoom);
    xfree(c->row[row_index]);

    /* only move if not the last row in array */
    if (row_index < c->num_rows - 1) {
	memmove(&c->row[row_index],
		&c->row[row_index + 1],
		(c->num_rows - row_index - 1) * sizeof(coord *));
    }

    for (i = row_index; i < c->num_rows-1; i++) {
	for (j = 0; j < c->num_columns; j++) {
	    memmove(&c->matrix[i][j],
		    &c->matrix[i + 1][j],
		    sizeof(element *));
	}
    }
     
    for (i = 0; i < c->num_columns; i++) {
	c->matrix[c->num_rows-1][i] = NULL;
    }

    c->num_rows--; 
}

void delete_column_from_container(container *c,
				  int row_index,
				  int column_index)
{
    int i, j;
    element *e;
    int cnt = 0;

#ifdef DEBUG
    printf("delete column %d %d\n", column_index, c->num_columns);
#endif

    if (row_index < 0) {
	c->num_columns--;
	return;
    }

    /*
     * need to modify the column_index for all those elements which will be moved
     * by the subsequent memmove
     */
    for (i = 0; i < c->num_rows; i++) {
	for (j = column_index; j < c->num_columns; j++) {
	    e = c->matrix[i][j];
	    if (e) {
		e->column_index--;
		cnt++;
	    }
	}
    }

    xfree(c->column[column_index]->pixel);
    freeZoom(&c->column[column_index]->zoom);
    xfree(c->column[column_index]);

#ifdef DEBUG
    printf("NUM COLUMNS TO MOVE %d\n", cnt);
#endif

    /* only move if not the last row in array */
    if (column_index < c->num_columns - 1) {
	memmove(&c->column[column_index],
		&c->column[column_index + 1],
		(c->num_columns - column_index - 1) * sizeof(coord *));
	for (i = 0; i < c->num_rows; i++) {
	    memmove(&c->matrix[i][column_index],
		    &c->matrix[i][column_index + 1],
		    (cnt) * sizeof(element *));
	}
    }
    
    for (i = 0; i < c->num_rows; i++) {
	c->matrix[i][c->num_columns] = NULL;
    }
    
    c->num_columns--; 
}

void update_row(element *e)
{
    container *c = e->c;

    if (e->row_index < 0)
	return;

    if (e->orientation & VERTICAL) {
	c->row[e->row_index]->visible.min = e->world->visible->x1;
	c->row[e->row_index]->visible.max = e->world->visible->x2;
    }
    set_pixel_coords(c->column[e->column_index]->visible.min, 
		     c->row[e->row_index]->visible.min, 
		     c->column[e->column_index]->visible.max, 
		     c->row[e->row_index]->visible.max, 
		     c->row[e->row_index]->pixel);

    container_update_scrollregion(e->c->interp, e->c);

}

void update_column(element *e)
{
    container *c = e->c;

    if (e->column_index < 0)
	return;

    if (e->orientation & HORIZONTAL) {
	if (e->world->visible->x1 > c->column[e->column_index]->visible.min)
	    c->column[e->column_index]->visible.min = e->world->visible->x1;
	if (e->world->visible->x2 < c->column[e->column_index]->visible.max)
	    c->column[e->column_index]->visible.max = e->world->visible->x2;
    }
    set_pixel_coords(c->column[e->column_index]->visible.min, 
		     c->row[e->row_index]->visible.min, 
		     c->column[e->column_index]->visible.max, 
		     c->row[e->row_index]->visible.max, 
		     c->column[e->column_index]->pixel);

    container_update_scrollregion(e->c->interp, e->c);
}

void delete_element_from_container(element *e)
{
    int i, j;
    container *c = e->c;
    int row_delete = 1;
    int column_delete = 1;
    char cmd[1024];
    int row_num, column_num;
    int num = 0;

#ifdef DEBUG
    printf("delete_element_from_container %s row_index %d column_index %d num cols %d num rows %d type %d\n", 
	   e->win, e->row_index, e->column_index, c->num_columns, c->num_rows, e->type);

    print_elements_in_container(c);
#endif

    if (c->num_rows == 1 && c->num_columns == 1) {
	delete_container(c);
	return;
    }

    /* if ((e->type == RULER_AMP) || (e->type == RULER_LEN)) {  */
    if (e->type == RULER_LEN) { 
	row_delete = 0;
	column_delete = 0;
    } else if (e->type == RULER_AMP) {

	Tcl_VarEval(c->interp, "is_empty_column ", e->win, NULL);
	if (!atoi(Tcl_GetStringResult(c->interp))) {
	    column_delete = 0;
	}
	Tcl_VarEval(c->interp, "is_empty_row ", e->win, NULL);
	if (!atoi(Tcl_GetStringResult(c->interp))) {
	    row_delete = 0;
	}
#if 0

	for (i = 0; i < c->num_columns; i++) {
	    if (i != e->column_index && c->matrix[e->row_index][i]) {
		row_delete = 0;
	    }
	}    
	for (i = 0; i < c->num_rows; i++) {
	    if (i != e->row_index && c->matrix[i][e->column_index]) {
		column_delete = 0;
	    }
	} 
#endif
    } else {
	/* check if need to delete entire row */
	for (i = 0; i < c->num_columns; i++) {
	    if (i != e->column_index && c->matrix[e->row_index][i] != NULL
		&& (c->matrix[e->row_index][i]->type != RULER_AMP) && (c->matrix[e->row_index][i]->type != RULER_LEN)) {
	    
		row_delete = 0;
		break;
	    }
	}
	/* check if need to delete entire column */
	for (i = 0; i < c->num_rows; i++) {
	    if (i != e->row_index && c->matrix[i][e->column_index] != NULL
		&& (c->matrix[i][e->column_index]->type != RULER_AMP) && (c->matrix[i][e->column_index]->type != RULER_LEN)) {
		column_delete = 0;
		break;
	    }
	}
    }

#ifdef DEBUG
    printf("row_delete %d column_delete %d\n", row_delete, column_delete);
#endif
    
    if (row_delete) {

	/* delete all windows in row */
	row_num = get_element_row(c->interp, e->win);

	sprintf(cmd, "tcl_delete_row %s %d", c->win, row_num);
	if (TCL_OK != (Tcl_Eval(c->interp, cmd)))
	    printf("tcl_delete_row %s\n", Tcl_GetStringResult(c->interp));

	/* length ruler */
	if (c->row[e->row_index]->ruler) {
	    c->row[e->row_index]->ruler->shutdown_func(c->row[e->row_index]->ruler, ALL);
	    c->row[e->row_index]->ruler = NULL;
	}
	/* amplitude ruler */
#ifdef DEBUG
	printf("Delete ruler %s %d\n", e->win, e->ruler_id);
#endif
	if (e->ruler_id > -1) {
	    delete_element(get_element(e->ruler_id), ALL);
	    e->ruler_id = -1;
	}
	delete_row_from_container(e->c, e->row_index, e->column_index);

	/* update_row(e); */
    } 
    if (column_delete) {
	/* delete all windows in column */
	column_num = get_element_column(c->interp, e->win);

	sprintf(cmd, "tcl_delete_column %s %d", c->win, column_num);
	if (TCL_OK != (Tcl_Eval(c->interp, cmd)))
	    printf("tcl_delete_column %s\n", Tcl_GetStringResult(c->interp));

	if (c->column[e->column_index]->ruler && (e != c->column[e->column_index]->ruler)) {

	    sprintf(cmd, "delete_sb_h %s %d", c->win, column_num);
	    Tcl_Eval(c->interp, cmd);

	    c->column[e->column_index]->ruler->shutdown_func(c->column[e->column_index]->ruler, ALL);
	    c->column[e->column_index]->ruler = NULL;
	}
	/* amplitude ruler */
	if (e->ruler_id > -1) {
#if 0
	    sprintf(cmd, "delete_sb_h %s %d", c->win, c->column[e->column_index]->num);
	    Tcl_Eval(c->interp, cmd);
#endif
	    delete_element(get_element(e->ruler_id), ALL);
	    e->ruler_id = -1;
	}

	delete_column_from_container(e->c, e->row_index, e->column_index);
	/* update_column(e); */
    }
    
    if (!row_delete && !column_delete) {
	c->matrix[e->row_index][e->column_index] = NULL;
	Tcl_VarEval(c->interp, "tcl_delete_element ", e->win, NULL);
    }

    if (e->type == RULER_LEN) {
	if (e->orientation == HORIZONTAL)
	    e->c->column[e->column_index]->ruler = NULL;
	else
	    e->c->row[e->row_index]->ruler = NULL;
    }

    /* if only RULER_LEN elements are left, then delete container */
    for (i = 0; i < c->num_rows; i++) {
	for (j = 0; j < c->num_columns; j++) {
	    if (c->matrix[i][j]) {
		if (c->matrix[i][j]->type != RULER_LEN) {
		    num++;
		}
	    }
	}
    }

    if (num == 0) {
	/* need to delete ruler */
	for (i = 0; i < c->num_rows; i++) {
	    for (j = 0; j < c->num_columns; j++) {
		if (c->matrix[i][j]) {
		    e->shutdown_func(c->matrix[i][j], ALL); 
		}
	    }
	}

	delete_container(c);
    }
}

int find_row_index(container *c,
		   int row_num,
		   int *new_row)
{
    int i, j;
    element *e;

    for (i = 0; i < c->num_rows; i++) {
	for (j = 0; j < c->num_columns; j++) {
	    if (e = c->matrix[i][j]) {
		if (row_num == get_element_row(c->interp, e->win)) {
		    *new_row = 0;
		    return e->row_index;
		}
	    }
	}
    }    
    *new_row = 1;
    return c->num_rows;
}

int find_column_index(container *c,
		      int column_num,
		      int *new_column)
{
    int i, j;
    element *e;

    for (i = 0; i < c->num_rows; i++) {
	for (j = 0; j < c->num_columns; j++) {
	    if (e = c->matrix[i][j]) {
		if (column_num == get_element_column(c->interp, e->win)) {
		    *new_column = 0;
		    return e->column_index;
		}
	    }
	}
    }    
    *new_column = 1;
    return c->num_columns;
}

void print_elements_in_container(container *c)
{
    int row, column;

    for (row = 0; row < c->num_rows; row++) {
	for (column = 0; column < c->num_columns; column++) {
	    printf("%p %d %d\n", c->matrix[row], row, column);
	    if (c->matrix[row][column])
		print_element(c->matrix[row][column]);
	}
    }
}

void push_row_zoom(StackPtr **zoom, mm visible)
{
    d_box v;

    v.x1 = 0;
    v.y1 = visible.min;
    v.x2 = 0;
    v.y2 = visible.max;

    pushZoom(zoom, &v);
}

void push_column_zoom(StackPtr **zoom, mm visible)
{
    d_box v;

    v.x1 = visible.min;
    v.y1 = 0;
    v.x2 = visible.max;
    v.y2 = 0;

    pushZoom(zoom, &v);
}

void container_update_scrollregion(Tcl_Interp *interp,
				   container *c)
{
    int i, j;
    element *e;

    for (i = 0; i < c->num_rows; i++) {
	for (j = 0; j < c->num_columns; j++) {
	    e = c->matrix[i][j];
	    if (e){
		/*
		if (e->scale_func) {
		    bbox = scale_box(e);

		    printf("BBOX %s %f %f %f %f\n", e->win, bbox.x1, bbox.x2,
			   bbox.y1, bbox.y2);
		    e->scale_func(interp, e, -1, &bbox, e->pixel);
		}
		*/
		if (e->scrollregion_func)
		    e->scrollregion_func(interp, e, e->world->total, 
					 c->column[e->column_index]->pixel,
					 c->row[e->row_index]->pixel);
	    }
	}
    }
}

int add_element_to_container(Tcl_Interp *interp,
			     int container_id,
			     char *c_win,
			     element *e,
			     int min_h,
			     int max_h,
			     int min_v,
			     int max_v)
{
    int row_num, column_num;
    container *c;
    int new_row = 0;
    int new_column = 0;
    int update_scrollregion = 0;
    int row_index, column_index;
    int draw_row_ruler = 0;
    int draw_column_ruler = 0;

    if (NULL == (c = get_container(container_id))) {
	c = create_container(interp, c_win, container_id);
    }

    /* get row and column numbers from tcl grid */
    row_num = get_element_row(interp, e->win);
    column_num = get_element_column(interp, e->win);

#ifdef DEBUG
    printf("ROW NUM %d COL %d\n", row_num, column_num);
#endif
        
    row_index = find_row_index(c, row_num, &new_row);
    column_index = find_column_index(c, column_num, &new_column);
    
#ifdef DEBUG
    printf("INDEX %d %d\n", row_index, column_index);
    printf("NEW_ROW %d NEW_COLUMN %d\n", new_row, new_column);
#endif

    /* if element already exists, delete it */
    
    if (row_index > 0 && column_index > 0 && c->matrix[row_index][column_index]) {
	delete_element(c->matrix[row_index][column_index], ALL);
    }

    e->c = c;

    if (e->orientation & HORIZONTAL && check_element_len_ruler(e))
	draw_column_ruler = 1;

    if (e->orientation & VERTICAL && check_element_len_ruler(e))
	draw_row_ruler = 1;

    if (c->num_rows == 0 && c->num_columns == 0) {
	init_container_matrix(c, row_num, column_num, &row_index, &column_index);
	new_row = 1;
	new_column = 1;
    } else {

	if (new_row) {
	    add_row_to_container(c, row_index, column_index, row_num, min_v,
				 max_v);
	}
	if (new_column) {
	    add_column_to_container(c, row_index, column_index, column_num, 
				    min_h, max_h); 
	}
    }

#ifdef DEBUG
    printf("INDEX %d %d\n", row_index, column_index);
#endif

    c->matrix[row_index][column_index] = e;
    e->row_index = row_index;
    e->column_index = column_index;
	
#ifdef DEBUG
    printf("add_element_to_container %d %d\n", row_index, column_index);
    print_elements_in_container(c);
    printf("min_h %d max_h %d min %f max %f\n", min_h, max_h,
	   c->column[column_index]->total.min, c->column[column_index]->total.max);
#endif

    if (min_h < c->column[column_index]->total.min) {
	c->column[column_index]->total.min = min_h;
	update_scrollregion = 1;
    }
    if (max_h > c->column[column_index]->total.max) {
	c->column[column_index]->total.max = max_h;    
	update_scrollregion = 1;
    }
    if (min_v < c->row[row_index]->total.min) {
	c->row[row_index]->total.min = min_v;
	update_scrollregion = 1;
    }
    if (max_v > c->row[row_index]->total.max) {
	c->row[row_index]->total.max = max_v;
	update_scrollregion = 1;
    }

    /* initialise visible */
    if (c->row[row_index]->visible.min == INT_MAX)
	c->row[row_index]->visible.min = c->row[row_index]->total.min;
    if (c->row[row_index]->visible.max == INT_MIN)
	c->row[row_index]->visible.max = c->row[row_index]->total.max;
    if (c->column[column_index]->visible.min == INT_MAX)
	c->column[column_index]->visible.min = c->column[column_index]->total.min;
    if (c->column[column_index]->visible.max == INT_MIN)
	c->column[column_index]->visible.max = c->column[column_index]->total.max;

    if (new_row) {
	init_pixel(interp, e, c->row[row_index]->pixel);
	set_pixel_coords(c->column[column_index]->visible.min, 
			 c->row[row_index]->visible.min, 
			 c->column[column_index]->visible.max, 
			 c->row[row_index]->visible.max, 
			 c->row[row_index]->pixel);

	push_row_zoom(&c->row[row_index]->zoom, c->row[row_index]->visible);

	if (draw_row_ruler && !c->row[row_index]->ruler) 
	    add_length_ruler(c->interp, c, row_index, column_index, row_num, column_num, VERTICAL);
    }

    if (new_column) {
	init_pixel(interp, e, c->column[column_index]->pixel);
	set_pixel_coords(c->column[column_index]->visible.min, 0,
			 c->column[column_index]->visible.max, 0,
			 c->column[column_index]->pixel);

	push_column_zoom(&c->column[column_index]->zoom, c->column[column_index]->visible);

	if (draw_column_ruler && !c->column[column_index]->ruler) 
	    add_length_ruler(c->interp, c, row_index, column_index, row_num, column_num, HORIZONTAL);

    }

    /* don't do this first time round */
    if (update_scrollregion && !(new_row && new_column)) {
	container_update_scrollregion(interp, c);
    }
    
    if (check_element_amp_ruler(e) && !e->ruler && e->orientation & HORIZONTAL) 
	add_element_ruler(interp, c, e, VERTICAL);
    if (check_element_amp_ruler(e) && !e->ruler && e->orientation & VERTICAL) 
	add_element_ruler(interp, c, e, HORIZONTAL);


    if (check_element_len_ruler(e) && !(new_row && new_column)) {
	if (e->orientation & HORIZONTAL) 
	    update_length_ruler(interp, c, c->column[column_index]);
	if (e->orientation & VERTICAL) 
	    update_length_ruler(interp, c, c->row[row_index]);
    }

    return 0;
}

int new_element(Tcl_Interp *interp,
		char **e_win)
{
    static int element_id = 0;
    char *win;

    win = get_default_string(interp, tk_utils_defs, w("ELEMENT.WIN"));

    if (NULL == (*e_win = (char *)xmalloc((strlen(win)+10) *  sizeof(char))))
	return -1;

    sprintf(*e_win, "%s%d", win, element_id);

    return element_id++;
}

element *create_element(Tcl_Interp *interp,
			int container_id,
			int element_id,
			char *e_win,
			int orientation,
			int crosshair)

{
    element *e;
#ifdef DEBUG
    printf("create_element %d %s\n", element_id, e_win);
#endif
    if (NULL == (e = (element *)xmalloc(sizeof(element))))
	return NULL;

    e->id = element_id;
    e->win = strdup(e_win);
    e->num_seqs = 0;
    e->max_seqs = 0;
    e->seqs = NULL;
    e->num_results = 0;
    e->max_results = 0;
    e->results = NULL;
    e->status = 0;
    e->container_id = container_id;
    e->c = NULL;
    e->crosshair = crosshair;

    if (NULL == (e->pixel = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return NULL;

    if (NULL == (e->world= (WorldPtr *)xmalloc(sizeof(WorldPtr))))
	return NULL;

    if (NULL == (e->world->visible = (d_box *)xmalloc(sizeof(d_box))))
	return NULL;
    
    if (NULL == (e->world->total = (d_box *)xmalloc(sizeof(d_box))))
	return NULL;

    createZoom(&e->zoom);

    e->world->total->x1 = DBL_MAX;
    e->world->total->x2 = -DBL_MAX;
    e->world->total->y1 = DBL_MAX;
    e->world->total->y2 = -DBL_MAX;

    e->world->visible->x1 = DBL_MAX;
    e->world->visible->x2 = -DBL_MAX;
    e->world->visible->y1 = DBL_MAX;
    e->world->visible->y2 = -DBL_MAX;

    e->type = -1;
    e->orientation = orientation;
    e->min_y = DBL_MAX;
    e->max_y = -DBL_MAX;
    e->min_x = INT_MAX;
    e->max_x = INT_MIN;
    e->cursor = NULL;
    e->cursor_array = NULL;
    e->cursor_visible = 0;
    e->ruler = NULL;
    e->ruler_id = -1;
    e->parent = -1;
    e->num_children_left = 0;
    e->num_children_right = 0;
    e->children_id = -1;
    e->children_position = 0;
    e->seq_e_id = -1;
    e->replot_func = NULL;

    return e;
}

void delete_element(element *e,
		    int job)
{
    int i;
    
    if (!e) {
	return;
    }
#ifdef DEBUG
    printf("***** delete_element %s\n", e->win);
#endif

    delete_element_from_container(e);

    if (e->ruler) {
	xfree(e->ruler->window);
	xfree(e->ruler);
    }
    xfree(e->pixel);
    xfree(e->world->visible);
    xfree(e->world->total);
    xfree(e->world);

    freeZoom(&e->zoom);

    if (job == ALL) {
	for (i = 0; i < e->num_results; i++) {
	    if (e->results[i]->n_configure > 0) {
		xfree(e->results[i]->configure[0]);
		if (e->results[i]->n_configure == 2) 
		    xfree(e->results[i]->configure[1]);
		xfree(e->results[i]->configure);
	    }
	    if (e->results[i]->colour)
		free(e->results[i]->colour);
	    xfree(e->results[i]);
	}
    }
    xfree(e->results);
    xfree(e->cursor);
    xfree(e->cursor_array);
    xfree(e->seqs);

    free(e->win);
    xfree(e);
}

void print_element(element *e) {

    printf("PRINT ELEMENT\n");

    printf("win %s ", e->win);
    printf("element_id %d ", e->id);
    printf("container_id %d ", e->container_id);
    printf("row_index %d ", e->row_index);
    printf("column_index %d ", e->column_index);
    printf("type %d ", e->type);
    printf("orientation %d\n", e->orientation);

    printf("\n\n");
}

element *get_element(int element_id)
{
    int i, r, c;
    
    for (i = 0; i < container_num; i++) {
#ifdef DEBUG
	printf("&&&&&&&&&&&&&&&&& get_element r= %d c= %d\n",
	       container_array[i]->num_rows, container_array[i]->num_columns);
#endif
	for (r = 0; r < container_array[i]->num_rows; r++) {
	    for (c = 0; c < container_array[i]->num_columns; c++) {
		if (container_array[i]->matrix[r][c]) {
#ifdef DEBUG
		    printf("matrix[%d][%d]=%s %d\n", r, c, container_array[i]->matrix[r][c]->win,
			   container_array[i]->matrix[r][c]->id);
#endif
		    if (container_array[i]->matrix[r][c]->id == element_id) {
			return (container_array[i]->matrix[r][c]);
		    }
		}
	    }
	}
    }
    return NULL;
}

char *element_window(char *c_win,
		     char *e_win)
{
    char *window;

    if (NULL == (window = (char *)xmalloc(strlen(c_win)+strlen(e_win) * 

					  sizeof(char))))
	return NULL;

    sprintf(window, "%s%s", c_win, e_win);
    return (window);

}

/*
 * find if seq_id of orientation exists and return container_id & element_id
 */
int find_seq_id(int seq_id,
		int orientation,
		int *container_id,
		int *element_id)
{
    int i, j, r, c;

    for (i = 0; i < container_num; i++) {
	for (r = 0; r < container_array[i]->num_rows; r++) {
	    for (c = 0; c < container_array[i]->num_columns; c++) {
		for (j = 0; j < container_array[i]->matrix[r][c]->num_seqs; j++) {
		    if ((seq_id == container_array[i]->matrix[r][c]->seqs[j].seq_id) && (orientation == container_array[i]->matrix[r][c]->seqs[j].direction)) {
			*container_id = container_array[i]->matrix[r][c]->container_id;
			*element_id = container_array[i]->matrix[r][c]->id;
			return 0;
		    }
		}
	    }
	}
    }
    return -1;
}

int set_element_type(element *e,
		     int type)
{
    e->type = type;

    if (e->type == SEQ_EDITOR) {
	e->scroll_x_func = NULL;
	e->scroll_y_func = NULL;
	e->scale_func = NULL;
	e->scrollregion_func = NULL;
	e->draw_crosshair_func = NULL;
	e->delete_crosshair_func = NULL;
	e->element_width = NULL;
	e->element_height = NULL;
	e->element_x = NULL;
	e->element_y = NULL;
    } else if (e->type == CANVAS) {
	e->scroll_x_func = canvas_scroll_x;
	e->scroll_y_func = canvas_scroll_y;
	e->scale_func = canvas_scale;
	e->scrollregion_func = canvas_scrollregion;
	e->draw_crosshair_func = draw_canvas_crosshair;
	e->delete_crosshair_func = delete_canvas_crosshair;
	e->element_width = canvas_width;
	e->element_height = canvas_height;
	e->element_x = canvas_x;
	e->element_y = canvas_y;
    } else if (e->type == RULER_AMP) {
	e->scroll_x_func = canvas_scroll_x;
	e->scroll_y_func = canvas_scroll_y;
	e->scale_func = canvas_scale;
	e->scrollregion_func = canvas_scrollregion;
	e->draw_crosshair_func = draw_canvas_crosshair;
	e->delete_crosshair_func = delete_canvas_crosshair;
	e->element_width = canvas_width;
	e->element_height = canvas_height;
	e->element_x = canvas_x;
	e->element_y = canvas_y;
    } else if (e->type == RULER_LEN){
	e->scroll_x_func = canvas_scroll_x;
	e->scroll_y_func = canvas_scroll_y;
	e->scale_func = canvas_scale;
	e->scrollregion_func = canvas_scrollregion;
	e->draw_crosshair_func = draw_canvas_crosshair;
	e->delete_crosshair_func = delete_canvas_crosshair;
	e->element_width = canvas_width;
	e->element_height = canvas_height;
	e->element_x = canvas_x;
	e->element_y = canvas_y;
    } else {
	verror(ERR_WARN, "set_element_type", "Invalid element type");
	return -1;
    }
    return 0;

}

int add_seq_id_to_element(element *e,
			  int seq_id,
			  int orientation)
{
    int i;
    int incr = 10;

#ifdef DEBUG
    printf("add_seq_id_to_element %d %d\n", seq_id, orientation);
#endif
    /* if no seq_id exists, need to add it */
    for (i = 0; i < e->num_seqs; i++) {
	if ((e->seqs[i].seq_id == seq_id) && (e->seqs[i].direction & orientation)) {
	    return 0;
	}
    }

    e->num_seqs++;
    if (e->num_seqs > e->max_seqs) {
	e->max_seqs += incr;
	if (NULL == (e->seqs = (seq_id_dir *)realloc(e->seqs, 
						     e->max_seqs * 
						    sizeof(seq_id_dir))))
	    return -1;
    }
    e->seqs[e->num_seqs-1].seq_id = seq_id;
    e->seqs[e->num_seqs-1].direction = orientation;

    return 0;
}

int add_result_to_element(element *e,
			  plot_data *result,
			  double min_x,
			  double min_y,
			  double max_x,
			  double max_y,
			  int orientation,
			  int element_type)
{
    int incr = 10;
    double m, c; 

    if (-1 == set_element_type(e, element_type)) 
	return -1;

    e->num_results++;
    if (e->num_results > e->max_results) {
	e->max_results += incr;
	if (NULL == (e->results = (plot_data **)realloc(e->results, 
							e->max_results * 
							sizeof(plot_data *))))
	    return -1;
    }
    e->results[e->num_results-1] = result;
    e->orientation |= orientation;

    if (e->num_results > 1) {
	m = (e->max_y - e->min_y) / (max_y - min_y);
	c = e->max_y - (m * max_y);
    } else {
	m = 1.0;
	c = 0.0;
    }
    /*
    result->sf_c = (m * result->sf_c) + c;
    result->sf_m *= m;
    */
    result->sf_c = c;
    result->sf_m = m;
#ifdef DEBUG
    printf("min_y %f %f max_y %f %f m %f c %f sf_c %f sf_m %f\n", 
	   min_y, e->min_y, max_y, e->max_y, m, c, result->sf_c, result->sf_m); 
#endif
    min_y = result->sf_m * min_y + result->sf_c;
    max_y = result->sf_m * max_y + result->sf_c;

#ifdef DEBUG
    printf("after %f %f\n", min_y, max_y);
    printf("y1 %f y2 %f\n", e->world->total->y1, e->world->total->y2);
#endif

    if (min_y < e->min_y)
	e->min_y = min_y;

    if (max_y > e->max_y)
	e->max_y = max_y;

    if (min_x < e->world->total->x1)
	e->world->total->x1 = min_x;

    if (max_x > e->world->total->x2)
	e->world->total->x2 = max_x;
  
    if (min_y < e->world->total->y1)
	e->world->total->y1 = min_y;

    if (max_y > e->world->total->y2)
	e->world->total->y2 = max_y;
    return 0;
}

void remove_result_from_element(element *e,
				int result_id)
{
    int e_num = -1;
    int i;

    for (i = 0; i < e->num_results; i++) {
	if (result_id == e->results[i]->result_id) {
	    e_num = i;
	    break;
	}
    }

    if (e_num == -1) {
	/* Error */
	return;
    }

    if (e_num < e->num_results-1) {
	memmove(&e->results[e_num], &e->results[e_num+1], 
		(e->num_results - e_num - 1) * sizeof(plot_data *));
    }

    e->num_results--;
    if (e->num_results == 0) {
#if 0
	if (e->ruler_id != -1) {
	    /* delete amplitude ruler */
	    /* e->shutdown_func(get_element(e->ruler_id), ALL); */
	    delete_element(get_element(e->ruler_id), ALL);
	    e->ruler_id = -1;
	}
#endif
	e->shutdown_func(e, WIN_ONLY); 
    }
}

/*
 * convert between element (pixel) coords and world (base) coords
 * element = world * m + c 
 * world = (element - c) / m
 */
void pixel_to_world(CanvasPtr *pixel,
		    int px, int py,
		    double *wx, double *wy)
{
   *wx = (px - pixel->bx) / pixel->ax;
   *wy = (py - pixel->by) / pixel->ay; 

#ifdef DEBUG
   printf("CanvasToWorld wx %f wy %f px %d py %d\n", *wx, *wy, px, py); 
#endif

}

/*
 * convert between world (base) coords and element (pixel) coords
 * pixel = world * m + c 
 */
void world_to_pixel(CanvasPtr *pixel,
		    double wx, double wy, 
		    int *px, int *py)
{
    *px = (int)ROUND((wx * pixel->ax) + pixel->bx);
    *py = (int)ROUND((wy * pixel->ay) + pixel->by);
/*
    printf("WorldToCanvas wx %f px %d ax %f bx %f\n", 
	   wx, *px, pixel->ax, pixel->bx);
*/ 
}

void set_pixel_coords(double x1, 
		      double y1, 
		      double x2, 
		      double y2,
		      CanvasPtr *pixel)     /* out */
{
    if (x2 - x1 == 0) {
	pixel->ax = 1;
    } else {
	pixel->ax = pixel->width / (x2-x1);
    }
    if (y2 - y1 == 0) {
	pixel->ay = 1;
    } else {
	pixel->ay = pixel->height / (y2-y1);
    }
    pixel->bx = pixel->x - (pixel->ax * x1);
    pixel->by = pixel->y - (pixel->ay * y1);

#ifdef DEBUG
    printf("********************SET PIXEL COORDS \n");
    printf("x1 %f x2 %f y1 %f y2 %f \n", x1, x2, y1, y2);
    printf("x %d ax %f ay %f bx %f by %f\n", pixel->x, pixel->ax, 
	   pixel->ay, pixel->bx, pixel->by);
#endif
}


void container_scroll_x (Tcl_Interp *interp,
			 int container_id,
			 int column_num,
			 char *command)
{
    container *c;
    element *e;
    int i;
    int column_index;
    double dummy;    
    int new_column;

    c = get_container(container_id);
    if (!c)
	return;

    column_index = find_column_index(c, column_num, &new_column);

    for (i = 0; i < c->num_rows; i++) {
	if (e = c->matrix[i][column_index])
	    if (e->scroll_x_func)
		e->scroll_x_func(interp, e, command);
    }

    if (e = c->matrix[0][column_index]) {
	
	c->column[e->column_index]->pixel->x = e->element_x(interp, e->win, 0.0);
	
	pixel_to_world(c->column[e->column_index]->pixel, 
		       c->column[e->column_index]->pixel->x, 0,
		       &c->column[e->column_index]->visible.min, &dummy);
	pixel_to_world(c->column[e->column_index]->pixel, 
		       c->column[e->column_index]->pixel->x + 
		       c->column[e->column_index]->pixel->width, 0,
		       &c->column[e->column_index]->visible.max, &dummy);
	
	set_pixel_coords(c->column[e->column_index]->visible.min, 0,
			 c->column[e->column_index]->visible.max, 0,
			 c->column[e->column_index]->pixel);
    }
}

void container_scroll_y (Tcl_Interp *interp,
			 int container_id,
			 int row_num,
			 char *command)
{
    container *c;
    element *e;
    int i;
    int row_index;
    double dummy;
    int new_row;

#ifdef DEBUG
    printf("container_scroll_y\n");
#endif

    c = get_container(container_id);
    if (!c)
	return;

    row_index = find_row_index(c, row_num, &new_row);

    for (i = 0; i < c->num_columns; i++) {
	if (e = c->matrix[row_index][i])
	    if (e->scroll_y_func)
		e->scroll_y_func(interp, e, command);
    }

#if 0
    /* kfs 6/1/03 don't know why this is here! */
    for (i = 0; i < c->num_rows; i++) {

	if (e = c->matrix[i][0]) {
    
	    c->row[i]->pixel->y = e->element_y(interp, e->win, 0.0);
	    
	    pixel_to_world(c->row[i]->pixel, 
			  0, c->row[i]->pixel->y,
			  &dummy, &c->row[i]->visible.min);
	    pixel_to_world(c->row[i]->pixel, 
			  0, c->row[i]->pixel->y + 
			  c->row[i]->pixel->height,
			  &dummy, &c->row[i]->visible.max);
	    
	    set_pixel_coords(0, c->row[i]->visible.min,
			     0, c->row[i]->visible.max,
			     c->row[i]->pixel);
	}
    }
#endif
    e = c->matrix[row_index][0];
    c->row[e->row_index]->pixel->y = e->element_y(interp, e->win, 0.0);
    
   pixel_to_world(c->row[e->row_index]->pixel, 
		  0, c->row[e->row_index]->pixel->y,
		  &dummy, &c->row[e->row_index]->visible.min);
    pixel_to_world(c->row[e->row_index]->pixel, 
		   0, c->row[e->row_index]->pixel->y + 
		   c->row[e->row_index]->pixel->height,
		   &dummy, &c->row[e->row_index]->visible.max);
    
    set_pixel_coords(0, c->row[e->row_index]->visible.min,
		     0, c->row[e->row_index]->visible.max,
		     c->row[e->row_index]->pixel);
   
}


/*
 * convert a zoom amount (from button) into a bounding box
 */
void container_convert_zoom(element *e,
			    float amount,
			    box *bbox)
{
    double length, dist;
    double wx1, wx2, wy1, wy2;

    length = e->world->visible->x2 - e->world->visible->x1 + 1;
    dist = length * amount;
    wx1 = (int)(e->world->visible->x1 + dist);
    wx2 = (int)(e->world->visible->x2 - dist);

    length = e->world->visible->y2 - e->world->visible->y1 + 1;
    dist = length * amount;
    wy1 = e->world->visible->y1 + dist;
    wy2 = e->world->visible->y2 - dist;

    world_to_pixel(e->pixel, wx1, wy1, &bbox->x1, &bbox->y1);
    world_to_pixel(e->pixel, wx2, wy2, &bbox->x2, &bbox->y2);
}

void element_zoom(Tcl_Interp *interp,
		  element *e,
		  float amount,
		  int x1,
		  int y1,
		  int x2,
		  int y2)
{
    d_box bbox;
    box zoom;
    double w_x1, w_y1, w_x2, w_y2;
    container *c = e->c;

    if (!e->scale_func)
	return;

    if (amount != -1) {
	container_convert_zoom(e, amount, &zoom);
    } else {
	zoom.x1 = x1;
	zoom.y1 = y1;
	zoom.x2 = x2;
	zoom.y2 = y2;
    }

#ifdef DEBUG
    printf("ELEMENT ZOOM %d %d %d %d\n", zoom.x1, zoom.x2, zoom.y1, zoom.y2);
#endif

    if (e->world->visible->x1 == DBL_MAX || e->world->visible->x2 == -DBL_MAX
	|| e->world->visible->y1 == DBL_MAX || e->world->visible->y2 == -DBL_MAX) {
	return;
    }
    
    w_x1 = e->world->visible->x1;
    w_y1 = e->world->visible->y1;
    w_x2 = e->world->visible->x2;
    w_y2 = e->world->visible->y2;
    
#ifdef DEBUG
     printf("ZOOM before e %f %f c %f %f\n", e->world->visible->x1,
	   e->world->visible->x2, c->column[e->column_index]->visible.min,
	   c->column[e->column_index]->visible.max);
#endif

    pixel_to_world(e->pixel, zoom.x1, zoom.y1, 
		   &e->world->visible->x1, 
		   &e->world->visible->y1);
    pixel_to_world(e->pixel, zoom.x2, zoom.y2, 
		   &e->world->visible->x2, 
		   &e->world->visible->y2);
    
#ifdef DEBUG
    printf("ZOOM after e %f %f c %f %f\n", e->world->visible->x1,
	   e->world->visible->x2, c->column[e->column_index]->visible.min,
	   c->column[e->column_index]->visible.max);
#endif

    bbox.x1 = (double)zoom.x1;
    bbox.y1 = (double)zoom.y1;
    bbox.x2 = (double)zoom.x2;
    bbox.y2 = (double)zoom.y2;
	    
    /* 
     * also if drawn a zooming box, only want to zoom in y that window 
     */
    if (amount == -1) {
#ifdef FIXME
	win_num = get_consistency_win_num(c, result_id);
	
	if (i != win_num) {
	    c->win_list[i]->world->visible->y1 = w_y1;
	    c->win_list[i]->world->visible->y2 = w_y2;
	    bbox.y1 = 0;
	    bbox.y2 = 0;
	}
#endif
    }
    
    set_pixel_coords(e->world->visible->x1, e->world->visible->y1, 
		     e->world->visible->x2, e->world->visible->y2, 
		     e->pixel);

    e->scale_func(interp, e, -1, &bbox, e->pixel);

    e->scrollregion_func(interp, e, e->world->total, 
			 c->column[e->column_index]->pixel,
			 c->row[e->row_index]->pixel);
    
    e->pixel->x = e->element_x(interp, e->win, 0.0);
    e->pixel->y = e->element_y(interp, e->win, 0.0);
    
    pushZoom(&e->zoom, e->world->visible);
    
}

void element_zoomback(Tcl_Interp *interp,
		      element *e)
{
    d_box bbox;
    int dummy;
    coord *column = e->c->column[e->column_index];
    coord *row = e->c->row[e->row_index];
    CanvasPtr *pixel_column = NULL;
    CanvasPtr *pixel_row = NULL;
    box zoom;
    
#ifdef DEBUG
    printf("element_zoomback %s\n", e->win);
#endif
    if (!e->scale_func)
	return;

    popZoom(&e->zoom);
    if (examineZoom(e->zoom) == NULL) {
	return;
    }
    memcpy(e->world->visible,  examineZoom(e->zoom), sizeof(d_box));

    world_to_pixel(e->pixel, e->world->visible->x1, e->world->visible->y1, 
		  &zoom.x1, &zoom.y1);
    world_to_pixel(e->pixel, e->world->visible->x2, e->world->visible->y2, 
		   &zoom.x2, &zoom.y2);

    if (e->orientation & HORIZONTAL) {
	world_to_pixel(column->pixel, column->visible.min, 0,
		      &zoom.x1, &dummy);
	world_to_pixel(column->pixel, column->visible.max, 0,
		      &zoom.x2, &dummy);
    }
    if (e->orientation & VERTICAL) {
	world_to_pixel(row->pixel, 0, row->visible.min,
		       &dummy, &zoom.y1);
	world_to_pixel(row->pixel, 0, row->visible.max,
		       &dummy, &zoom.y2);
    }
    
    bbox.x1 = (double)zoom.x1;
    bbox.y1 = (double)zoom.y1;
    bbox.x2 = (double)zoom.x2;
    bbox.y2 = (double)zoom.y2;
    
    e->scale_func(interp, e, -1, &bbox, e->pixel);
    
    set_pixel_coords(e->world->visible->x1, e->world->visible->y1, 
		     e->world->visible->x2, e->world->visible->y2, 
		     e->pixel);

    if (e->orientation & HORIZONTAL) {
	if (NULL == (pixel_column = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	    return;
	
	memcpy(pixel_column, e->c->column[e->column_index]->pixel, sizeof(CanvasPtr));
	set_pixel_coords(column->visible.min, 0, column->visible.max, 0,
			 pixel_column);
    }

    if (e->orientation & VERTICAL) {
	if (NULL == (pixel_row = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	    return;
	
	memcpy(pixel_row, e->c->row[e->row_index]->pixel, sizeof(CanvasPtr));
	
	set_pixel_coords(0, row->visible.min, 0, row->visible.max,
			 pixel_row);
    }
    
    e->scrollregion_func(interp, e, e->world->total, pixel_column, pixel_row);

    e->pixel->x = e->element_x(interp, e->win, 0.0);
    e->pixel->y = e->element_y(interp, e->win, 0.0);

    if (pixel_column)
	xfree(pixel_column);
    if (pixel_row)
	xfree(pixel_row);
}

void container_zoom(Tcl_Interp *interp,
		    int container_id,
		    float amount,
		    int x1,
		    int y1,
		    int x2,
		    int y2)
{
    int i, j;
    container *c = get_container(container_id);
    element *e;
    box **zoom;
    double length, dist;
    /* int dummy;*/
    double ddummy;

    if (NULL == (zoom = (box **)xmalloc(c->num_rows * sizeof(box *))))
	return;

    for (i = 0; i < c->num_rows; i++) {
	if (NULL == (zoom[i] = (box *)xmalloc(c->num_columns * sizeof(box))))
	    return;
    }

    for (i = 0; i < c->num_rows; i++) {
	if (c->row[i]->ruler) {
	    if (amount != -1.0) {
		length = c->row[i]->visible.max - c->row[i]->visible.min + 1;
		dist = length * amount;
		c->row[i]->visible.min += dist;
		c->row[i]->visible.max -= dist;
	    } else {
		pixel_to_world(c->row[i]->pixel, x1, y1, 
			       &ddummy, &c->row[i]->visible.min);

		pixel_to_world(c->row[i]->pixel, x2, y2, 
			       &ddummy, &c->row[i]->visible.max);
	    }
#ifdef REMOVE
	    /* kfs 11/04/2003 
	     * can't see a reason for setting zoom here, it is overwritten
	     * in element_zoom anyway
	     */
	    world_to_pixel(c->row[i]->pixel, 0, c->row[i]->visible.min, 
			   &dummy, &zoom[i][0].y1);
	    world_to_pixel(c->row[i]->pixel, 0, c->row[i]->visible.max, 
		       &dummy, &zoom[i][0].y2);
#endif
	    set_pixel_coords(0, c->row[i]->visible.min,
			     0, c->row[i]->visible.max,
			     c->row[i]->pixel);
	}
    }

    for (i = 0; i < c->num_columns; i++) {
	if (c->column[i]->ruler) {
	    if (amount != -1.0) {
		length = c->column[i]->visible.max - c->column[i]->visible.min + 1;
		dist = length * amount;
		c->column[i]->visible.min += dist;
		c->column[i]->visible.max -= dist;
	    } else {
		pixel_to_world(c->column[i]->pixel, x1, y1, 
			       &c->column[i]->visible.min, &ddummy);

		pixel_to_world(c->column[i]->pixel, x2, y2, 
			       &c->column[i]->visible.max, &ddummy);
	    }
#ifdef REMOVE
	    /* kfs 11/04/2003 
	     * can't see a reason for setting zoom here, it is overwritten
	     * in element_zoom anyway
	     */
	    world_to_pixel(c->column[i]->pixel, c->column[i]->visible.min, 
			   0, &zoom[0][i].x1, &dummy);
	    world_to_pixel(c->column[i]->pixel, c->column[i]->visible.max, 
			   0, &zoom[0][i].x2, &dummy);
#endif
	    
	    set_pixel_coords(c->column[i]->visible.min, 0,
			     c->column[i]->visible.max, 0,
			     c->column[i]->pixel);
	}
    }
    
    for (i = 0; i < c->num_rows; i++) {
	for (j = 0; j < c->num_columns; j++) {
	    if (e = c->matrix[i][j])
		element_zoom(interp, e, amount, x1, y1, x2, y2);
	}
    }
    
    for (i = 0; i < c->num_rows; i++) {
	/* need a window, use the first! */
	e = c->matrix[i][0];
	if (e != NULL && e->element_y) {
	    c->row[i]->pixel->y = e->element_y(interp, e->win, 0.0);
	    push_row_zoom(&c->row[i]->zoom, c->row[i]->visible);
	}
    }

    for (i = 0; i < c->num_columns; i++) {
	/* need a window, use the first! */
	e = c->matrix[0][i];
	if (e != NULL && e->element_x) {
	    c->column[i]->pixel->x = e->element_x(interp, e->win, 0.0);
	    push_column_zoom(&c->column[i]->zoom, c->column[i]->visible);
	}
    }
    for (i = 0; i < c->num_rows; i++)
	xfree(zoom[i]);
    xfree(zoom);
}

/* zoomback all elements in container */
void container_zoomback(Tcl_Interp *interp,
			int container_id)
{
    int i, j;
    element *e;
    container *c = get_container(container_id);
    d_box *visible;
   
    /* set visible to be previous zoom level for each row */
    for (i = 0; i < c->num_rows; i++) {
	if (c->row[i]->ruler) {
	    popZoom(&c->row[i]->zoom);
	    if (examineZoom(c->row[i]->zoom) == NULL) {
		return;
	    }
	    visible = examineZoom(c->row[i]->zoom);
	    
	    c->row[i]->visible.min = visible->y1;
	    c->row[i]->visible.max = visible->y2;
	}
    }

    /* set visible to be previous zoom level for each column */
    for (i = 0; i < c->num_columns; i++) {
	if (c->column[i]->ruler) {
	    popZoom(&c->column[i]->zoom);
	    if (examineZoom(c->column[i]->zoom) == NULL) {
		return;
	    }
	    visible = examineZoom(c->column[i]->zoom);
	    
	    c->column[i]->visible.min = visible->x1;
	    c->column[i]->visible.max = visible->x2;
	}
    }

    for (i = 0; i < c->num_rows; i++) {
	for (j = 0; j < c->num_columns; j++) {
	    e = c->matrix[i][j];
	    if (e != NULL) 
		element_zoomback(interp, e);
	}
    }   

    for (i = 0; i < c->num_rows; i++) {
	e = c->matrix[i][0];
	if (e != NULL && e->element_y) {
	    c->row[i]->pixel->y = e->element_y(interp, e->win, 0.0);
	    set_pixel_coords(0, c->row[i]->visible.min,
			     0, c->row[i]->visible.max,
			     c->row[i]->pixel);
	}
    }

    for (i = 0; i < c->num_columns; i++) {
	e = c->matrix[0][i];
	if (e != NULL && e->element_x) {
	    c->column[i]->pixel->x = e->element_x(interp, e->win, 0.0);
	    set_pixel_coords(c->column[i]->visible.min, 0, 
			     c->column[i]->visible.max, 0,
			     c->column[i]->pixel);
	}
    }

}

d_box scale_box(element *e)
{
    d_box bbox;
    
    bbox.x1 = e->world->visible->x1;
    bbox.x2 = e->world->visible->x2;
    bbox.y1 = e->world->visible->y1;
    bbox.y2 = e->world->visible->y2;

    /* horizontal length ruler */
    if (e->orientation & HORIZONTAL) {
	bbox.x1 = e->c->column[e->column_index]->visible.min;
	bbox.x2 = e->c->column[e->column_index]->visible.max;
    }
    
    /* vertical length ruler */
    if (e->orientation & VERTICAL) {
	bbox.y1 = e->c->row[e->row_index]->visible.min;
	bbox.y2 = e->c->row[e->row_index]->visible.max;
    }

#ifdef DEBUG
    printf("scale_box %f %f %f %f\n", bbox.x1, bbox.y1, bbox.x2, bbox.y2);
#endif
    return bbox;
}

void draw_container_crosshair(Tcl_Interp *interp,
			      int element_id,
			      int x,
			      int y)
{
    int i;
    element *e = get_element(element_id);
    element *next;

    /* draw crosshair in all elements in same column and row */
    if (e->crosshair & HORIZONTAL) {
	for (i = 0; i < e->c->num_rows; i++) {
	    next = e->c->matrix[i][e->column_index];
	    if (next) 
		e->draw_crosshair_func(interp, next, x, HORIZONTAL);
	}
    } 
    if (e->crosshair & VERTICAL) {
	for (i = 0; i < e->c->num_columns; i++) {
	    next = e->c->matrix[e->row_index][i];
	    if (next)
		e->draw_crosshair_func(interp, next, y, VERTICAL);
	}
    }
}

void delete_container_crosshair(Tcl_Interp *interp,
				int element_id)
{
    int i, j;
    element *e = get_element(element_id);
    element *next;

    if (!e)
	return;

    for (i = 0; i < e->c->num_rows; i++) {
	for (j = 0; j < e->c->num_columns; j++) {
	    next = e->c->matrix[i][j];
	    if (next != NULL)
		e->delete_crosshair_func(interp, next);
	}
    }
	
}

double invert_wy (element *e,
		 double y) 
{

    return(e->world->total->y2 - y + e->world->total->y1);
}

int invert_cy (element *e,
	       int y) 
{

    return(e->pixel->height - y);
}

cursor_t *find_element_cursor(element *e,
			      int seq_id,
			      int direction)
{
    int i;

    for (i = 0; i < e->num_seqs; i++) {
	if (e->seqs[i].seq_id == seq_id && 
	    e->seqs[i].direction == direction) {
	    return e->cursor[i];
	}
    }
    return NULL;
}

void set_container_column_coords(Tcl_Interp *interp,
				 coord *column,
				 int x)
{
    double dummy;

    column->pixel->x = x;
    
#ifdef DEBUG
    printf("set_container_column_coords %d\n", column->pixel->x);
#endif
    pixel_to_world(column->pixel, 
		  column->pixel->x, 0,
		  &column->visible.min, &dummy);
    
    pixel_to_world(column->pixel, 
		  column->pixel->x + 
		  column->pixel->width, 0,
		  &column->visible.max, &dummy);
    
    set_pixel_coords(column->visible.min, 0, column->visible.max, 0,
		     column->pixel);
    
}

int make_new_element(Tcl_Interp *interp,
		     int c_id,
		     char *c_win,
		     int e_id,
		     char *e_win, 
		     int seq_id,
		     plot_data *result,
		     double min_x,
		     double max_x,
		     double min_y,
		     double max_y,
		     int orientation,
		     int crosshair,
		     int element_type,
		     element **e_new)
{
    element *e;
    double x0, x1, y0, y1;

    if (NULL == (e = (get_element(e_id))))
	e = create_element(interp, c_id, e_id, e_win, orientation, crosshair);
    
    if (-1 == (add_result_to_element(e, result, min_x, min_y, max_x, max_y, 
				     orientation, element_type))) {
	return -1;
    }
    add_seq_id_to_element(e, seq_id, orientation);

    x0 = min_x;
    x1 = max_x;
    y0 = min_y;
    y1 = max_y;

    if (!(orientation & HORIZONTAL)) {
	x0 = INT_MAX;
	x1 = INT_MIN;
    }
    if (!(orientation & VERTICAL)) {
	y0 = INT_MAX;
	y1 = INT_MIN;
    }

    e->world->total->x1 = min_x;
    e->world->total->x2 = max_x;
    e->world->total->y1 = min_y;
    e->world->total->y2 = max_y;

    init_pixel(interp, e, e->pixel);

    add_element_to_container(interp, c_id, c_win, e, x0, x1, y0, y1);

#ifdef DEBUG    
    printf("COLUMN TOTAL %f %f VISIBLE %f %f\n", 
	   e->c->column[e->column_index]->total.min,
	   e->c->column[e->column_index]->total.max,
	   e->c->column[e->column_index]->visible.min,
	   e->c->column[e->column_index]->visible.max);
#endif

    if (orientation == HORIZONTAL) {
	e->world->visible->x1 = e->c->column[e->column_index]->total.min;
	e->world->visible->x2 = e->c->column[e->column_index]->total.max;
	e->world->visible->y1 = e->world->total->y1;
	e->world->visible->y2 = e->world->total->y2;
    } else if (orientation == VERTICAL) {
	e->world->visible->x1 = e->world->total->x1;
	e->world->visible->x2 = e->world->total->x2;
	e->world->visible->y1 = e->c->row[e->row_index]->total.min;
	e->world->visible->y2 = e->c->row[e->row_index]->total.max;
    } else {
	e->world->visible->x1 = e->c->column[e->column_index]->total.min;
	e->world->visible->x2 = e->c->column[e->column_index]->total.max;
	e->world->visible->y1 = e->c->row[e->row_index]->total.min;
	e->world->visible->y2 = e->c->row[e->row_index]->total.max;
    }

    set_pixel_coords(e->world->visible->x1, e->world->visible->y1, 
		     e->world->visible->x2, e->world->visible->y2, e->pixel);

    /* remove all current zooming info */
    freeZoom(&e->zoom);

    /* add first zoom */
    pushZoom(&e->zoom, e->world->visible);
     
    *e_new = e;
    return 0;
}

plot_data *find_plot_data(element *e,
			  int result_id)
{
    int i;

    for (i = 0; i < e->num_results; i++) {
	if (result_id == e->results[i]->result_id)
	    return e->results[i];
    }
    return NULL;
}

void move_element_to_new_container(Tcl_Interp *interp,
				   int old_e_id,
				   int new_c_id,
				   char *new_c_win,
				   int old_c_id,
				   char *new_e_win,
				   int new_e_id,
				   int new_orientation)
{
    element *e_new;
    int i;

    if (NULL == (e_new = get_element(old_e_id))) {
	printf("ERROR in move_element_to_new_container\n");
	return;
    }

    delete_element_from_container(e_new);

    /* add element to new container */

    /* change orientation of all seqs in element */
    if (new_orientation != e_new->orientation) {
	for (i = 0; i < e_new->num_seqs; i++) {
	    e_new->seqs[i].direction = new_orientation;
	}
    }

    e_new->win = strdup(new_e_win);
    e_new->id = new_e_id;
    e_new->orientation = new_orientation;


    add_element_to_container(interp, new_c_id, new_c_win, e_new, 
			     e_new->world->total->x1, e_new->world->total->x2, 
			     e_new->world->total->y1, e_new->world->total->y2);
    e_new->replot_func(e_new);
    
}

int get_coord_seq_ids(container *c,
		      int index,
		      int orientation,
		      seq_id_dir **seq_ids,
		      int *num_seqs) {
    int i, j;
    element *e;
    int cnt = 0;
    
    /* iterate down a column */
    if (orientation == VERTICAL) {

	/* first time round to see how much room to allocate */
	for (i = 0; i < c->num_rows; i++) {
	    if (!(e = c->matrix[i][index]))
		break;
	    for (j = 0; j < e->num_seqs; j++) {
		if (e->seqs[j].direction == orientation) 
		    cnt++;
	    }
	}

	if (NULL == (*seq_ids = (seq_id_dir *)xmalloc(cnt * sizeof(seq_id_dir))))
	    return -1;

	cnt = 0;
	/* fill in seq_id structure */
 	for (i = 0; i < c->num_rows; i++) {
	    if (!(e = c->matrix[i][index]))
		break;
	    for (j = 0; j < e->num_seqs; j++) {
		if (e->seqs[j].direction == orientation) {
		    (*seq_ids)[cnt].seq_id =  e->seqs[j].seq_id;
		    (*seq_ids)[cnt++].direction =  e->seqs[j].direction;
		}
	    }
	}
    } else {

	/* first time round to see how much room to allocate */
	for (i = 0; i < c->num_columns; i++) {
	    if (!(e = c->matrix[index][i]))
		break;
	    for (j = 0; j < e->num_seqs; j++) {
		if (e->seqs[j].direction == orientation) 
		    cnt++;
	    }
	}

	if (NULL == (*seq_ids = (seq_id_dir *)xmalloc(cnt * sizeof(seq_id_dir))))
	    return -1;

	cnt = 0;
	/* fill in seq_id structure */
 	for (i = 0; i < c->num_columns; i++) {
	    if (!(e = c->matrix[index][i]))
		break;
	    for (j = 0; j < e->num_seqs; j++) {
		if (e->seqs[j].direction == orientation) {
		    (*seq_ids)[cnt].seq_id =  e->seqs[j].seq_id;
		    (*seq_ids)[cnt++].direction =  e->seqs[j].direction;
		}
	    }
	}
    }

    *num_seqs = cnt;
    return 0;
}

void element_resize(Tcl_Interp *interp,
		      int e_id)
{
    element *e = get_element(e_id);
    int x1, x2, y1, y2;
    int width, height;

#ifdef DEBUG
    printf("element_resize %d\n", e_id);
#endif

    if (!e)
	return;

    x1 = e->pixel->x;
    y1 = e->pixel->y;
    x2 = e->pixel->x + e->pixel->width;
    y2 = e->pixel->y + e->pixel->height;
   
    width = e->element_width(interp, e->win);
    height = e->element_height(interp, e->win);

    if (width !=  e->pixel->width || height != e->pixel->height) {
	e->pixel->width = width;
	e->pixel->height = height;

	if (e->orientation & HORIZONTAL) {
	    e->c->column[e->column_index]->pixel->width = width;
	    e->c->column[e->column_index]->pixel->height = height;	
	}
	if (e->orientation & VERTICAL) {
	    e->c->row[e->row_index]->pixel->width = width;
	    e->c->row[e->row_index]->pixel->height = height;
	}

	element_zoom(interp, e, -1, x1, y1, x2, y2);

	if (e->orientation & HORIZONTAL) {
	    set_pixel_coords(e->c->column[e->column_index]->visible.min, 0,
			     e->c->column[e->column_index]->visible.max, 0,
			     e->c->column[e->column_index]->pixel);
	}
	if (e->orientation & VERTICAL) {
	    set_pixel_coords(0, e->c->row[e->row_index]->visible.min,
			     0, e->c->row[e->row_index]->visible.max,
			     e->c->row[e->row_index]->pixel);
	}
	if (e->replot_func) 
	    e->replot_func(e);
    }
}

void rotate_element(element *e,
		    int seq_id, 
		    int result_id)
{
    int i;
    container *c = e->c;
    int row_index = -1;
    char cmd[1024];
    int row_num;

    if (e->orientation == HORIZONTAL) {
	/* can only rotate an element if there is a vertical ruler in the same
	   row */
	for (i = 0; i < c->num_rows; i++) {
	    if (c->row[i]->ruler && c->row[i]->ruler->orientation == VERTICAL) {
		row_index = i;
		break;
	    }
	}
	if (row_index == -1)
	    return;
	
	row_num = get_element_row(c->interp, e->win);
	sprintf(cmd, "rotate_element %s %s %d %d %d %d", e->win, c->row[i]->ruler->win, seq_id, result_id, VERTICAL, row_num);
	
	if (TCL_OK != (Tcl_Eval(c->interp, cmd)))
	    printf("rotate_element!!! %s\n", Tcl_GetStringResult(c->interp));
    }
}

void update_container_menu(int container_id) {
    int i, j;
    container *c;
    char cmd[1024];
    
    if (!(c = get_container(container_id)))
	return;

    /* in the process of shutting down container */
    if (c->status)
	return;

    Tcl_VarEval(c->interp, "delete_menubar ", c->win, NULL);

    for (i = 0; i < c->num_rows; i++) {
	for (j = 0; j < c->num_columns; j++) {
	    if (c->matrix[i][j]) {
		sprintf(cmd, "update_container_menu %s %d %s", c->win, c->id,
			c->matrix[i][j]->win);
		Tcl_Eval(c->interp, cmd);
	    }
	}
    }
}
