#ifndef _NIP_BASE_COMP_H
#define _NIP_BASE_COMP_H

typedef struct pbc_in {
    char *params;
} in_plot_base_comp;

int init_nip_base_comp_create(Tcl_Interp *interp, int seq_id, int start, 
			      int end, int win_len, int a, int c, int g, 
			      int t, int *id);

int init_nip_base_comp_plot(Tcl_Interp *interp, int seq_id, int result_id,
			    char *raster_win, int raster_id, char *colour, 
			    int line_width);
#endif

