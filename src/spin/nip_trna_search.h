#ifndef _NIP_TRNA_SEARCH_H_
#define _NIP_TRNA_SEARCH_H_

#include "trna_search.h"

#define MAX_TRNA_BP 44

typedef struct trna_in {
    char *params;
    TrnaSpec *t;
} in_trna_search;

typedef struct trna_ {
    TrnaRes **results;
    int n_pts;
} trna;

int init_nip_trna_search_plot(Tcl_Interp *interp, int seq_id, int result_id,
			      char *raster_win, int raster_id, char *colour,
			      int line_width, int tick_ht);

int init_nip_trna_search_create(Tcl_Interp *interp, int seq_id, int start,
				int end, int *id);
#endif
