#ifndef _SIP_QUICK_SCAN_H_
#define _SIP_QUICK_SCAN_H_

/* quick scan parameters */
typedef struct qs_in {
    char *params;
} in_quick_scan;

int init_sip_best_diagonals_create(Tcl_Interp *interp, int seq_id_h,
				   int seq_id_v, int start_h, int end_h, 
				   int start_v, int end_v, int win_len,
				   int min_match, int word_len, float min_sd,
				   int *id);

int init_sip_best_diagonals_plot(Tcl_Interp *interp, int seq_id_h, 
				 int seq_id_v, int result_id, 
				 char *raster_win, int raster_id,
				 char *colour, int line_width);

#endif
