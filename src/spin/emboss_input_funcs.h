#ifndef _EMBOSS_INPUT_FUNCS_H
#define _EMBOSS_INPUT_FUNCS_H

typedef struct in_emboss_ {
    char *params;
} in_emboss;

int init_emboss_graph_create(Tcl_Interp *interp, 
			     int seq_id, 
			     int start, 
			     int end,
			     char *filename,
			     int *id);

int init_emboss_graph_plot(Tcl_Interp *interp, 
			   int seq_id, 
			   int result_id,
			   char *name,
			   char *raster_win,
			   int raster_id,
			   char *colour, 
			   int line_width);

int init_emboss_dot_create(Tcl_Interp *interp, 
			   int seq_id_h, 
			   int start_h, 
			   int end_h,
			   int seq_id_v, 
			   int start_v, 
			   int end_v,
			   char *filename,
			   int *id);

int init_emboss_dot_plot(Tcl_Interp *interp, 
			 int seq_id_h,
			 int seq_id_v,
			 int result_id,
			 char *name,
			 char *raster_win,
			 int raster_id,
			 char *colour, 
			 int line_width);

int init_emboss_stick_create(Tcl_Interp *interp, int seq_id, int start, 
			     int end, char *filename, int *id);
#endif
