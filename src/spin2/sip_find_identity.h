#ifndef _SIP_FIND_IDENTITY_H_
#define _SIP_FIND_IDENTITY_H_

/* #define MAX_MATCHES 1000 */

/* identities parameters */
typedef struct fi_in {
    char *params;
} in_find_identities;

typedef struct text_find_identities_ {
    int word_len;
} text_find_identities;

int init_sip_matching_words_create(Tcl_Interp *interp, int seq_id_h,
				   int seq_id_v, int start_h, int end_h, 
				   int start_v, int end_v, int word_len,
				   int strand, Tcl_Obj **graph_obj, int *id);

int init_sip_matching_words_plot(Tcl_Interp *interp,
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
