#ifndef _SIP_COMPARE_SPANS_H_
#define _SIP_COMPARE_SPANS_H_

/* compare spans parameters */
typedef struct cs_in {
    char *params;
    int win_length;
} in_sim_spans;

typedef struct text_sim_spans_ {
    int min_score;
    int *match_score;
    int num;
} text_sim_spans;

int init_sip_similar_spans_create(Tcl_Interp *interp, int seq_id_h,
				  int seq_id_v, int start_h, int end_h, 
				  int start_v, int end_v, int win_len,
				  int min_match, int strand, int score, 
				  Tcl_Obj **graph_obj, int *id);

int init_sip_similar_spans_plot(Tcl_Interp *interp,
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
