#ifndef _NIP_GENE_SEARCH_H
#define _NIP_GENE_SEARCH_H

#include "codon_content.h"

typedef struct pgs_in {
    char *method;
    char *params;
} in_plot_gene_search;

int init_nip_codon_pref_create(Tcl_Interp *interp, int seq_id, int start,
			       int end, char *codon_table, int win_len,
			       int option, int *id);

int init_nip_gene_search_plot(Tcl_Interp *interp, char *rasters, 
			      char *raster_ids, int seq_id,  char *r_id,
			      char *colours, int line_width);

int init_nip_base_bias_create(Tcl_Interp *interp, int seq_id,int start,
			      int end, int win_len, int *id);

int init_nip_base_bias_plot(Tcl_Interp *interp, char *raster1_win, 
			    char *raster_id1, int seq_id,  char *r_id,
			    char *colour1, int line_width);

int init_nip_author_test_create(Tcl_Interp *interp, int seq_id, int start,
				int end, char *codon_table, double error,
				int *id);
#endif
