#ifndef _NIP_GENE_SEARCH_H
#define _NIP_GENE_SEARCH_H

#include "codon_content.h"

typedef struct pgs_in {
    char *method;
    char *params;
} in_plot_gene_search;

int init_nip_codon_pref_create(Tcl_Interp *interp, int seq_id, int start,
			       int end, char *codon_table, int win_len,
			       int option, int strand, Tcl_Obj **graph_obj, 
			       int *id);

int init_nip_base_bias_create(Tcl_Interp *interp, int seq_id,int start,
			      int end, int win_len, int strand, 
			      Tcl_Obj **graph_obj, int *id);

int init_nip_author_test_create(Tcl_Interp *interp, int seq_id, int start,
				int end, char *codon_table, double error,
				int strand, Tcl_Obj **graph_obj, int *id);

int init_nip_gene_search_plot(Tcl_Interp *interp,
			       int seq_id,
			       int result_id,
			       char *e_win,
			       char *c_win, 
			       Tcl_Obj *results,
			       int container_id,
			       int element_id,
			       char *element_type,
			       int line_width,
			       char *colour,
			      int orientation);
#endif
