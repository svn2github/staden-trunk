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
				 int *id);

int init_sip_global_align_plot(Tcl_Interp *interp, int seq_id_h, int seq_id_v,
			       int result_id, char *raster_win,  int raster_id,
			       char *colour, int line_width);

#endif
