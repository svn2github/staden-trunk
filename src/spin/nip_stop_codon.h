#ifndef _NIP_STOP_CODON_H
#define _NIP_STOP_CODON_H

#include "renz_utils.h"

#define PLUS 1
#define COMPLEMENT 2
#define BOTH 3

#define STOP_CODON  0
#define START_CODON 1

typedef struct in_s_codon_ {
    char *params;
} in_s_codon;

typedef struct s_codon_res_ {
    int *match1;
    int n_match1;
    int *match2;
    int n_match2;
    int *match3;
    int n_match3;
} s_codon_res;

int init_nip_start_codons_create(int seq_id, int start, int end,
				 char *strand_sym, int *id);

int init_nip_stop_codons_create(int seq_id, int start, int end,
				char *strand_sym, int *id);

int init_nip_stop_codons_plot(Tcl_Interp *interp, 
			      char *rasters, 
			      char *raster_ids,
			      int seq_id, 
			      char *r_id,
			      char *colours,
			      int line_width,
			      int tick_ht);

#endif

