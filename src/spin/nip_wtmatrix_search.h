#ifndef _NIP_WTMATRIX_SEARCH_H_
#define _NIP_WTMATRIX_SEARCH_H_

typedef struct in_wtmatrix_search_ {
    char *params;
} in_wtmatrix_search;

typedef struct text_wtmatrix_ {
    int mark_pos;
    int length;
} text_wtmatrix;

int init_nip_wtmatrix_search_create(int start, int end, int seq_id,
				    char *wt_matrix, int *id);

int init_nip_wtmatrix_search_plot(Tcl_Interp *interp, 
				  int seq_id,  
				  int result_id,
				  char *raster_win, 
				  int raster_id,
				  char *colour,
				  int line_width,
				  int tick_ht);
#endif
