#ifndef _NIP_BASE_COMP_H
#define _NIP_BASE_COMP_H

#include "tkCanvGraph.h"
#include "tclCanvGraph.h"

typedef struct pbc_in {
    char *params;
} in_plot_base_comp;

int init_nip_base_comp_create(Tcl_Interp *interp, int seq_id, int start, 
			      int end, int win_len, int a, int c, int g, 
			      int t, int strand, Tcl_Obj **g_data, int *id);

int init_nip_base_comp_plot(Tcl_Interp *interp, int seq_id, int result_id,
			     char *e_win, char *c_win, Tcl_Obj *results,
			     int container_id, int element_id,
			     char *element_type, int line_width,
			     char *colour, int orientation);
#endif

