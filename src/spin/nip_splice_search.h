#ifndef _NIP_SPLICE_SEARCH_H_
#define _NIP_SPLICE_SEARCH_H_
#include "tcl.h"
#include "splice_search.h"

typedef struct in_splice_ {
    char *params;
} in_splice;

int init_splice_search_create(int seq_id, int start, int end, char *donor,
			      char *acceptor, int *id);

int init_splice_search_plot(Tcl_Interp *interp, char *raster, int raster_id, 
			    char *result_id, int seq_id,  char *colours, 
			    int line_width, float tick_ht);

#endif

