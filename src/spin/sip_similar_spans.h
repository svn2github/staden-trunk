#ifndef _SIP_COMPARE_SPANS_H_
#define _SIP_COMPARE_SPANS_H_

/* compare spans parameters */
typedef struct cs_in {
    char *params;
} in_comp_spans;

typedef struct text_sim_spans_ {
    int min_score;
} text_sim_spans;

int init_sip_similar_spans_create(Tcl_Interp *interp, int seq_id_h,
				  int seq_id_v, int start_h, int end_h, 
				  int start_v, int end_v, int win_len,
				  int min_match, int *id);

int init_sip_similar_spans_plot(Tcl_Interp *interp, int seq_id_h, 
				int seq_id_v, int result_id, char *raster_win, 
				int raster_id, char *colour, int line_width);

#endif
