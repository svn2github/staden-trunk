#ifndef _SIP_ALIGNMENT_H_
#define _SIP_ALIGNMENT_H_
#include "align.h"

/* alignment parameters */
typedef struct a_in {
    char *params;
} in_align;

int init_sip_global_align_create(Tcl_Interp *interp, int seq_id_h,
				 int seq_id_v, int start_h, int end_h, 
				 int start_v, int end_v, int match,
				 int mis_match, int start_gap, int cont_gap,
				 int strand, Tcl_Obj **graph_obj, int *id);

int init_sip_global_align_plot(Tcl_Interp *interp,
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
				  char *element_type);
#endif
