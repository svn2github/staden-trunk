#ifndef _NIP_SPLICE_SEARCH_H_
#define _NIP_SPLICE_SEARCH_H_
#include "tcl.h"
#include "splice_search.h"

typedef struct in_splice_ {
    char *params;
} in_splice;

int init_splice_search_create(int seq_id, int start, int end, char *donor,
			      char *acceptor, int strand, 
			      Tcl_Obj **graph_obj, int *id);

int init_splice_search_plot(Tcl_Interp *interp, int seq_id, int result_id,
			    char *e_win, char *c_win, Tcl_Obj *results,
			    int container_id, int element_id,
			    char *element_type, int line_width,
			    char *colour, float tick_ht, int orientation);

#endif

