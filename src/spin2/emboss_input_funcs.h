#ifndef _EMBOSS_INPUT_FUNCS_H
#define _EMBOSS_INPUT_FUNCS_H

typedef struct in_emboss_ {
    char *params;
} in_emboss;

typedef struct text_emboss {
    char *title;
    char *maintitle;
    char *subtitle;
    char *xtitle;
    char *ytitle;
} text_emboss;

int init_emboss_graph_create(Tcl_Interp *interp, 
			     int seq_id, 
			     int start, 
			     int end,
			     char *filename,
			     Tcl_Obj **graph_obj,
			     int *id);

int init_emboss_graph_plot(Tcl_Interp *interp,
			   int seq_id_h,
			   int seq_id_v,
			   int result_id,
			   char *e_win, 
			   int element_id,
			   char *c_win, 
			   int container_id,
			   Tcl_Obj *results,
			   int line_width,
			   char *colour,
			   char *element_type,
			   char *name);

int init_emboss_dot_create(Tcl_Interp *interp, 
			   int seq_id_h, 
			   int start_h, 
			   int end_h,
			   int seq_id_v, 
			   int start_v, 
			   int end_v,
			   char *filename,
			   Tcl_Obj **graph_obj,
			   int *id);

int init_emboss_dot_plot(Tcl_Interp *interp,
			 int seq_id_h,
			 int seq_id_v,
			 int result_id,
			 char *e_win, 
			 int element_id,
			 char *c_win, 
			 int container_id,
			 Tcl_Obj *results,
			 int line_width,
			 char *colour,
			 char *element_type,
			 char *name);

int init_emboss_stick_create(Tcl_Interp *interp, int seq_id, int start, 
			     int end, char *filename, int *id);
#endif
