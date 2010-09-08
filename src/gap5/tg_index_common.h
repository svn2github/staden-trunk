#ifndef _TG_FUNC_H_
#define _TG_FUNC_H_

#include <tg_index.h>

bttmp_t *bttmp_file_open(void);

void bttmp_file_close(bttmp_t *tmp);

void bttmp_file_store(bttmp_t *tmp,  size_t name_len, char *name, tg_rec rec);

void bttmp_file_sort(bttmp_t *tmp);

char *bttmp_file_get(bttmp_t *tmp, tg_rec *rec);



tg_rec save_sequence(GapIO *io, seq_t *seq, bin_index_t *bin, range_t *r_out);

void find_pair(GapIO *io, HacheTable *pair, tg_rec recno, char *tname,
	       bin_index_t *bin, contig_t *c, seq_t *seq, tg_args *a,
	       range_t *r_out, library_t *lib); 
			   
tg_rec save_range_sequence(GapIO *io, seq_t *seq, uint8_t mapping_qual,
			   HacheTable *pair, int is_pair, char *tname,
			   contig_t *c, tg_args *a, int flags, library_t *lib);
					 
void create_new_contig(GapIO *io, contig_t **c, char *cname, int merge);

int open_tmp_file(void);
void close_tmp_file(void);
int sort_pair_file (void);
void complete_pairs(GapIO *io);			 

/*
 * Turns a comma separated list of data types into a bit-field.
 */
int parse_data_type(char *type);

int tg_index_file_type (char *fn);

void unescape_line(char *txt);

#endif


