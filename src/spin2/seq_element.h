#ifndef _SEQ_ELEMENT_H_
#define _SEQ_ELEMENT_H_

#include <tcl.h>
#include "container.h"
#include "seq_reg.h"
#include "tkCanvGraph.h"
#include "tclCanvGraph.h"
#include "seq_results.h"

#define MAX_NUM_SEQ 100

typedef struct seq_element_s {
    void (*op_func)(int seq_num,
		    void *obj,
		    seq_reg_data *data);
    element *e;

} seq_element;

int seq_element_reg(Tcl_Interp *interp, element *e);

element_info get_element_info(Tcl_Interp *interp, seq_id_dir *seq_ids,
			      int num_seqs, int plot_type, int frame,
			      int new_container);

void seq_find_cursor (int element_id, int direction, int *seq_id,
		      int *cursor_id);
 
void set_seq_element_funcs(element *e); 

void seq_element_cursor_notify(Tcl_Interp *interp, int seq_num, element *e,
			       cursor_t *cursor, int show_all);

int seq_canvas_move_cursor(int element_id, int cursor_id, int pos,
			   int direction);

void seq_plot_graph_func(void *obj, seq_reg_plot *plot);
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
		      char *job);
void config_result(int e_id,
		   int result_id,
		   int line_width,
		   char *colour);
void seq_ruler_reg(element *e);
int seq_find_column(Tcl_Interp *interp, int seq_id);
void seq_plot_stick_func(void *obj, seq_reg_plot *plot);
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
		      char *element_type);
int get_element_type(element *e);
int get_element_gr_type(element *e);
int compare_g_pt(const void *p1, const void *p2);
int compare_gd_line(const void *p1, const void *p2);

#endif
