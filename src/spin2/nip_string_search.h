#ifndef _NIP_STRING_SEARCH_H_
#define _NIP_STRING_SEARCH_H_

typedef struct in_string_search_ {
    char *params;
    char *string;
} in_string_search;

int init_nip_string_search_create(int strand_sym, float match, char *string,
				  int use_iub_code, int start, int end,
				  int seq_id, Tcl_Obj **graph_obj, int *id);

int init_nip_string_search_plot(Tcl_Interp *interp, int seq_id, int result_id,
				char *e_win, char *c_win, Tcl_Obj *results,
				int container_id, int element_id,
				char *element_type, int line_width,
				char *colour, float tick_ht, int orientation);
#endif
