#ifndef _SEQ_RASTER_H_
#define _SEQ_RASTER_H_

#include <tcl.h>
#include "seq_results.h"

/* plot_style */
#define LINE     0
#define DOT      1

/* result information jobs */
#define INPUT      0
#define OUTPUT     1
#define DIMENSIONS 2
#define INDEX      3
#define RESULT     4
#define WIN_SIZE   5
#define WIN_NAME   6

#define NEW_WINDOW 0
#define SEQ_WINDOW 1
#define MAX_NUM_SEQ 100 /* max num of seq to be displayed in single raster */

/* raster types */
#define ZOOM_BOX  0   /* sip */
#define SCALE_BAR 1   /* nip */

#define TASK_RASTER_ZOOM           1
#define TASK_RASTER_QUIT_RESULTS   2
#define TASK_RASTER_UPDATE_ALL     3
#define TASK_RASTER_UPDATE_SINGLE  4
#define TASK_RASTER_MAX_SIZE       5
#define TASK_RASTER_WINDOW_SIZE    6

#define SUPER 0
#define ADD   1
#define NEW   2


typedef struct config_ {
    float position; /* 0.0 to 1.0 */
    char x_direction; /* + or - (complemented) */
    char y_direction; /* + or - */
    float height;   /* <=1 for %height or >1 for pixel */
    int zoom;       /* 0=not zoomable, 1=zoom only on its own 2=zoom always */
    int scroll;    /* 0 or 1 */
} config;

typedef struct r_out {
    Tcl_Interp *interp;      /* tcl interptor */
    d_line max_len1;          /* min & max values for x and y in plot */
    int line_width1;          /* width of line to plot */
    int hidden;              /* boolean: whether plot is hidden or visible */
    int env_index;           /* environment for plotting, required by raster */
    char raster_win[1024];   /* name of raster window */
    int raster_id;           /* unique identifier of raster window */
    char colour1[100];        /* colour of plot */
    double max_y1;            /* max required y value of plot */
    double min_y1;           /* min required y value of plot */
    int plot_style;          /* plotting LINEs or DOTs */
    char scroll; /* whether need to scroll both 'b', x only 'x' or y only 'y'*/
    config **configure;      /* how to position, zoom and scroll the result */
    int n_configure;         /* number of configurations */
    double sf_m;             /* scale factor */
    double sf_c;             /* scale factor */
    int replot;
    char *name;              /* name to use for brief line, keybox etc */
} out_raster;

typedef struct update_struct_ {
    char *raster;
    int id;
    int old_id;
    int job;
} update_struct;

typedef struct cursor_info_ {
    int env;
    int visible[2];
    int prev_pos;
} cursor_info;

typedef struct _cursor_res_ {
    int cursor_id;
    int seq_id;
    int direction;
} cursor_res;

typedef struct _seq_id_dir {
    int seq_id;
    int direction; /* HORIZONTAL or VERTICAL */
} seq_id_dir;

typedef struct rasterresult {
    void (*op_func)(int seq_num,
		     void *obj,
		     seq_reg_data *data);
    Tcl_Interp *interp;
    char raster_win[1024];
    int id;
    seq_id_dir *seq;
    int num_seq_id;
    int num_results;
    int ed_cursor; /* FIXME - do I need this? current position of editor cursor */
    cursor_t **cursor;
    int status; /* 0 - inactive, 1 - active used to prevent reentrant calls when
		 * shutting down a raster */
    cursor_info cursor_array[NUM_CURSORS];
} RasterResult;

void init_raster_colour(Tcl_Interp *interp);
char * get_raster_colour(void);
void ReplotAllZoom(Tcl_Interp *interp, char *raster_win);

void ReplotAllCurrentZoom(Tcl_Interp *interp, char *raster_win);
int ReplotAllRasterWindow(Tcl_Interp *interp, char *raster_win);

RasterResult *raster_id_to_result(int raster_id);
RasterResult *raster_name_to_result(Tcl_Interp *interp, char * raster_win);

cursor_t *find_raster_result_cursor(RasterResult *result, int seq_id,
				    int direction);

int ReSetRasterWindowWorld(Tcl_Interp *interp, char *raster_old,
			   double orig_y1, int raster_type);

int add_seq_to_raster(RasterResult *result, int seq_id, int seq_num,
		      int direction, int line_width,
		      void (*func)(int seq_num, void *fdata, 
				   seq_reg_data *jdata));

void delete_seq_from_raster(int seq_id, 
			    int seq_num,
			    RasterResult *result,
			    void (*func)(int seq_num, void *fdata, 
					 seq_reg_data *jdata));

void raster_update_cursor(RasterResult *result, cursor_t *cursor,
			  int seq_id, Tk_Raster *raster, int cursor_show,
			  int direction);

int seq_raster_move_cursor(int raster_id, Tk_Raster *raster, int cursor_id,
			   int pos, int direction);

int seq_raster_find_edcursor(int raster_id, Tk_Raster *raster, int pos,
			     int direction, int *seq_id);

int FindRasterSize(int raster_id, d_line **max_size);

void UpdateVRuler(Tcl_Interp *interp, char *raster_win, double wy0,
		  double wy1);

void RemoveVRuler(Tcl_Interp *interp, char *raster_win, int id);

void update_raster_cursor(int new_id, int old_id);

int GetRasterZoom(int raster_id);

void UpdateScaleBars(Tcl_Interp *interp, double old_xscroll, 
		     double new_xscroll, double old_yscroll, 
		     double new_yscroll, char *raster_new, int id_old,   
		     int id_new, int job);

void remove_all_raster_cursors(Tcl_Interp *interp, Tk_Raster *raster,
			       RasterResult *result);

void raster_update_cursors(RasterResult *result, Tk_Raster *raster);

int raster_init_env(Tcl_Interp *interp, Tk_Raster *raster, cursor_t *cursor);

int raster_cursor_refresh(Tcl_Interp *interp, Tk_Raster *raster, 
			  cursor_t *changed_c, cursor_t *raster_c, int seq_id,
			  RasterResult *result, int cursor_show,
			  int direction);

int SetRasterWindowSize(Tcl_Interp *interp, char *raster_win);
int SeqReSetRasterWindowSize(Tcl_Interp *interp, char *raster_win,
			     int raster_type);
void UpdateZoomList(Tcl_Interp *interp,
		    double o_wx0, double o_wy0, double o_wx1, double o_wy1,
		    double n_wx0, double n_wy0, double n_wx1, double n_wy1,
		    char *parent, int job);
int GetWindowNumResults(Tcl_Interp *interp, char *raster_win);
char **GetRasterWindowList(Tcl_Interp *interp, char *raster_win,
			 int *num_windows);
char **GetRasterIdList(Tcl_Interp *interp, char *raster_win, int *num_windows);
void AddResultToRaster(RasterResult *result);
int DeleteResultFromRaster(RasterResult *result);
double rasterY(Tk_Raster *raster, double y);
void FindRasterResultY0(Tk_Raster *raster, int raster_id,
			config *configure, int num_results, double *y0,
			double *tick_height);
int SeqAddRasterToWindow(Tcl_Interp *interp, char *raster_win, 
			 int raster_type);
void seq_raster_callback(int seq_num, void *obj, seq_reg_data *jdata);
void SeqSuperimposeResult(Tcl_Interp *interp,char *raster_win, int result_id,
			  double o_wx0, double o_wy0, double o_wx1,
			  double o_wy1);
int raster_select_cursor_graph(int raster_id, Tk_Raster *raster,
			       char *raster_win, int pos, int max_dist);
int raster_select_cursor_dot(int raster_id, Tk_Raster *raster, 
			     char *raster_win, int rx, int ry, int max_dist,
			     int *cursor_id_h, int *cursor_id_v);
char *get_raster_frame_graph(Tcl_Interp *interp, int seq_id, int type,
			     int frame);
int get_raster_frame_dot(Tcl_Interp *interp, int seq_id_h, int seq_id_v,
			 char *raster_win);
int seq_raster_reg(Tcl_Interp *interp, char *raster_win, seq_id_dir *seq_array,
		   int num_seq_id);
int SeqDeregisterRasterWindow(int seq_id, RasterResult *old_result,
			      char *raster_win);

int init_gene_search_raster(Tcl_Interp *interp, 
			    int num_items,
			    char **win_list, 
			    char **id_list,
			    int seq_id, 
			    char **result_id,
			    char **colour_list,
			    int line_width);

int init_stick_raster(Tcl_Interp *interp, 
		      int seq_id,  
		      int result_id,
		      char *raster_win, 
		      int raster_id,
		      config *configure,
		      char *colour,
		      int line_width,
		      int tick_ht);

int init_graph_raster(Tcl_Interp *interp, 
		      int seq_id,  
		      int result_id,
		      char *raster_win, 
		      int raster_id,
		      config *configure,
		      char *colour,
		      int line_width);

int init_dot_plot(Tcl_Interp *interp, int seq_id_h, int seq_id_v, 
		  int result_id, char *name, char *raster_win, int raster_id,
		  char **opts, int n_opts, int plot_style, d_line dim);

#endif
