#ifndef _CONTAINER_H_
#define _CONTAINER_H_

#include <tcl.h>
#include "canvas_box.h"

/* types of element */
#define CANVAS      0
#define RASTER      1
#define SEQ_EDITOR  2
#define RULER_AMP   3
#define RULER_LEN   4

#define ROW     (1<<0)
#define COLUMN  (1<<1)

#define NEITHER -1
#define HORIZONTAL (1<<0)
#define VERTICAL (1<<1)
#define CIRCLE (1<<2)

#define TOP_S (1<<0)
#define BOTTOM_S (1<<1)

#define SCALE_X  (1<<0)
#define SCALE_Y  (1<<1)
#define SCROLL_X (1<<0)
#define SCROLL_Y (1<<1)

/* jobs for delete_element */

#define ALL      0
#define WIN_ONLY 1

/* type of plot data */
#define GRAPH 0
#define CANV_LINE 1

/* position of element movement highlight */
#define TOP    1
#define MIDDLE 2
#define BOTTOM 3
#define LEFT   4
#define RIGHT  5

typedef struct _seq_id_dir {
    int seq_id;
    int direction; /* HORIZONTAL or VERTICAL */
} seq_id_dir;

typedef struct cursor_info_ {
    int env;
    int visible[3];
    int prev_pos;
} cursor_info;

/*
 * The cursor contains information about the current display/editing cursors.
 * Cursors can be shared. A private cursor is private for one controlling
 * display, and can only be shared with other displays that do not need
 * access to their own private cursor.
 */
typedef struct _cursor_t {
    int id;		/* Cursor identification number */
    int refs;		/* Number of concurrent uses for this cursor */
    int private;	/* Whether this is a private cursor */
    int abspos; 	/* Absolute position in sequence */
    int job;    	/* move, reused, or delete */
    char *colour;       /* colour of cursor */
    int line_width;     /* line width of cursor */
    int direction;      /* direction of cursor */
    struct _cursor_t *next;
} cursor_t;

typedef struct configs_ {
    float position; /* 0.0 to 1.0 */
    char x_direction; /* + or - (complemented) */
    char y_direction; /* + or - */
    float height;   /* <=1 for %height or >1 for pixel */
    int zoom;       /* 0=not zoomable, 1=zoom only on its own 2=zoom always */
    int scroll;    /* 0 or 1 */
} configs;

struct container_s;
struct element_s;

typedef struct plot_data_ {
    int result_id;
    cursor_t *cursor;
    int cursor_visible;
    configs **configure;      /* how to position, zoom and scroll the result */
    int n_configure;         /* number of configurations */
    double sf_m;             /* scale factor */
    double sf_c;             /* scale factor */
    int scale;               /* SCALE_X, SCALE_Y */
    int hidden;
    char *colour;
    int line_width;
    char tags[10];
    int len_ruler;          /* draw a length ruler (1) or not (0) */
    int amp_ruler;          /* draw a amplitude ruler (1) or not (0) */
    int strand;             /* TOP_S or BOTTOM_S */
    void *extra_data;       /* extra data types */
    char *name;             /* name for display in keyname, brief etc */
} plot_data;

typedef struct child_ {
    int id;
    int position;           /* -1 for left, +1 for right */
} child;

typedef struct element_s {
    int container_id;          /* container id */
    struct container_s *c;     /* container structure */

    int id;
    char *win;
    WorldPtr *world;     
    CanvasPtr *pixel;   
    StackPtr *zoom;      
    int type;        /* RULER or RASTER or CANVAS or SEQ_EDITOR */
    int orientation; /* location of length ruler HORIZONTAL or VERTICAL */
    int crosshair;   /* where to draw crosshairs: HORIZONTAL, VERTICAL */
    
    plot_data **results;
    int num_results;
    int max_results; /* currently allocated slots */

    double max_y;     /* total max amplitude */
    double min_y;     /* total min amplitude */
    int min_x;        /* total min length */
    int max_x;        /* total max length */
    ruler_s *ruler;   /* amplitude ruler configuration structure */
    int ruler_id;     /* amplitude ruler result id*/
    
    int row_index;
    int column_index;

    seq_id_dir *seqs;     /* seq id and direction array */
    int num_seqs;    
    int max_seqs;

    cursor_t **cursor;
    cursor_info *cursor_array;
    int cursor_visible;

    /* canvas specific */
    void (*scale_func)(Tcl_Interp *interp, struct element_s *e, int result_id,
		       d_box *current, CanvasPtr *canvas);
    void (*scrollregion_func)(Tcl_Interp *interp, struct element_s *e,
			      d_box *total, CanvasPtr *pixel_column,
			      CanvasPtr *pixel_row);
 
    void (*scroll_x_func)(Tcl_Interp *interp, struct element_s *e, 
			  char *command);
    void (*scroll_y_func)(Tcl_Interp *interp, struct element_s *e, 
			  char *command);
    void (*draw_crosshair_func) (Tcl_Interp *interp, struct element_s *e,
				 int pos, int orientation);
    void (*delete_crosshair_func) (Tcl_Interp *interp, struct element_s *e);
    int (*element_width) (Tcl_Interp *interp, char *win);
    int (*element_height) (Tcl_Interp *interp, char *win);
    double (*element_x) (Tcl_Interp *interp, char *win, double value);
    double (*element_y) (Tcl_Interp *interp, char *win, double value);

    /* registration specific */
    void (*replot_func)(struct element_s *e);
    void (*shutdown_func)(struct element_s *e, int job);

    int status; /* 0 - inactive, 1 - active used to prevent reentrant calls 
	* when shutting down an element */
    int seq_e_id; /* registration id */
    int parent;   /* parent id, -1 for none */

    /* FIXME: should be an array but at present only have 1 child (ruler) */
    int children_id; 
    int children_position;
    int num_children_left; 
    int num_children_right; 
} element;

typedef struct mm_ {
    double min;
    double max;
} mm;

typedef struct coord_s {
    mm visible;
    mm total;
    CanvasPtr *pixel;
    StackPtr *zoom;
    element *ruler; /* length ruler */
    int ruler_reg;  /* whether ruler has been registered */
} coord;

typedef struct container_s {
    Tcl_Interp *interp;
    char *win;
    int id;

    element ***matrix; /* matrix of elements */
    coord **row;      /* summary of row info */
    coord **column;   /* summary of column info */
    int num_rows;    /* number of rows */
    int max_rows;    /* currently allocated slots */
    int num_columns; /* number of columns */
    int max_columns; /* currently allocated slots */

    ruler_s *ruler; /* length ruler */
    tick_s *tick;   /* ruler ticks */
    int status; /* 0=inactive, 1=active used signal shutting down container */
} container;

/* transfer info between tcl keyedlist and c structure */
typedef struct element_info_s {
    int element_id;
    int container_id;
    char *element_win;
    char *container_win;
    int orientation;
    char *e_win_to;
} element_info;

void container_start_shutdown(int container_id);

int new_container(Tcl_Interp *interp, char **c_win);

int new_element(Tcl_Interp *interp, char **e_win);

container *create_container(Tcl_Interp *interp, char *c_win,
			    int container_id); 

container *get_container(int container_id);

element *create_element(Tcl_Interp *interp,
			int container_id,
			int element_id,
			char *e_win,
			int orientation,
			int crosshair);

element *get_element(int element_id);

int find_container(seq_id_dir *seq_ids,
		   int num_seqs,
		   int *direction,
		   char **e_win,
		   char **c_win);

int add_element_to_container(Tcl_Interp *interp,
			     int container_id,
			     char *c_win,
			     element *e,
			     int min_h,
			     int max_h,
			     int min_v,
			     int max_v);

int add_seq_id_to_element(element *e,
			  int seq_id,
			  int orientation);

int add_result_to_element(element *e,
			  plot_data *result,
			  double min_x,
			  double min_y,
			  double max_x,
			  double max_y,
			  int orientation,
			  int element_type);

void container_scroll_x (Tcl_Interp *interp,
			 int container_id,
			 int column,
			 char *command);

void container_scroll_y (Tcl_Interp *interp,
			 int container_id,
			 int row_num,
			 char *command);

void print_element(element *e);
void print_elements_in_container(container *c);

void container_zoom(Tcl_Interp *interp,
		    int container_id,
		    float amount,
		    int x1,
		    int y1,
		    int x2,
		    int y2);

void container_zoomback(Tcl_Interp *interp,
			int container_id);

void set_pixel_coords(double x1, double y1, double x2, double y2,
		      CanvasPtr *pixel);

d_box scale_box(element *e);
double invert_wy (element *e, double y); 
int invert_cy (element *e, int y); 
cursor_t *find_element_cursor(element *e, int seq_id, int direction);
void world_to_pixel(CanvasPtr *pixel, double wx, double wy, int *px, int *py);
void pixel_to_world(CanvasPtr *pixel, int px, int py, double *wx, double *wy);
void set_container_column_coords(Tcl_Interp *interp,
				 coord *column,
				 int x);

plot_data *find_plot_data(element *e, int result_id);

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
		     element **e_new);

int get_coord_seq_ids(container *c,
		      int index,
		      int orientation,
		      seq_id_dir **seq_ids,
		      int *num_seqs);
void delete_element_from_container(element *e);
void delete_element(element *e, int job);
void remove_result_from_element(element *e, int result_id);
void move_element_to_new_container(Tcl_Interp *interp,
				   int old_e_id,
				   int new_c_id,
				   char *new_c_win,
				   int old_c_id,
				   char *new_e_win,
				   int new_e_id,
				   int new_orientation);
int get_element_row(Tcl_Interp *interp, char *win);
int get_element_column(Tcl_Interp *interp, char *win);
int check_element_scale(element *e);
void delete_container_crosshair(Tcl_Interp *interp, int element_id);
void element_resize(Tcl_Interp *interp, int e_id);
void draw_container_crosshair(Tcl_Interp *interp, int element_id, int x, int y);
void update_container_menu(int container_id);

#endif
