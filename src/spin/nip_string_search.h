#ifndef _NIP_STRING_SEARCH_H_
#define _NIP_STRING_SEARCH_H_

typedef struct in_string_search_ {
    char *params;
    char *string;
} in_string_search;

int init_nip_string_search_create(char *strand, float match, char *string,
				  int use_iub_code, int start, int end,
				  int seq_id, int *id);

int init_nip_string_search_plot(Tcl_Interp *interp, char *raster_win, 
				int raster_id, int result_id, int seq_id,
				char *colour, int line_width, int tick_ht);

#endif
