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
				   int *id);

int init_sip_matching_words_plot(Tcl_Interp *interp, int seq_id_h, 
				 int seq_id_v, int result_id, 
				 char *raster_win, int raster_id,
				 char *colour, int line_width);
#endif
