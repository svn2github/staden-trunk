#ifndef _SIP_SIM_H_
#define _SIP_SIM_H_

#include "array.h"
#include "seq_raster.h"
#include "align.h"

/* sim parameters */
typedef struct s_in {
    char *params;
} in_sim;

int init_sip_local_align_create(Tcl_Interp *interp, 
				int seq_id_h,
				int seq_id_v, 
				int start_h,
				int end_h, 
				int start_v,
				int end_v, 
				int num_align,
				float score_align,
				float match,
				float transition,
				float transversion,
				float start_gap,
				float cont_gap,
				int *id);

int init_sip_local_align_plot(Tcl_Interp *interp, 
			      int seq_id_h, 
			      int seq_id_v,
			      int result_id,
			      char *raster_win, 
			      int raster_id,
			      char *colour, 
			      int line_width);


void sim_align(char *seq1, char *seq2, int seq1_len, int seq2_len, 
	       int seq_type, int *num_alignments, float score_align,
	       float match, float transition,
	       float transversion, float start_gap, float cont_gap,
	       align_int **res, long *start1, long *start2, long *end1,
	       long *end2);

void sim_save(char *A, char *B, long seq1_len, long seq2_len,long M, long N,
	       align_int *S, long AP, long BP, Array ptarray, int *k);

int SaveSim(int seq1_num, int seq2_num, in_sim *input, out_raster *output,
	    Array ptarray);

void find_seq_lengths(align_int *S, long M, long N, int *len1, int *len2);
#endif
