#ifndef _TG_FUNC_H_
#define _TG_FUNC_H_

#include <tg_index.h>

int save_sequence(GapIO *io, seq_t *seq, bin_index_t *bin, range_t *r_out);

void find_pair(GapIO *io, HacheTable *pair, int recno, char *tname, bin_index_t *bin, contig_t *c,
               seq_t *seq, tg_args *a, range_t *r_out, library_t *lib);
			   
int save_range_sequence(GapIO *io, seq_t *seq, uint8_t mapping_qual, HacheTable *pair, int is_pair, char *tname, contig_t *c,
                  	 tg_args *a, int flags, library_t *lib);
					 
void create_new_contig(GapIO *io, contig_t **c, char *cname, int merge);
			 
#endif


