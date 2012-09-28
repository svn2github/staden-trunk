#ifndef _FIJ_H_
#define _FIJ_H_

#include "io_utils.h"
#include "consen.h"
#include "hash_lib.h"
#include "newgap_structs.h"
#include "io_lib/hash_table.h"
#define COMPARE_ALL 0
#define COMPARE_ONE_BY_ONE 1

typedef struct contig_pair {
    tg_rec c1;
    tg_rec c2;
} contig_pair;

typedef struct fij_opts_t {
    int mask;
    int min_overlap;
    double max_percent_mismatch;
    int word_len;
    double max_prob;
    int min_match;
    int band;
    int window_len;
    int max_unknown;
    double min_conf;
    int use_conf;
    int use_hidden;
    int max_alignment;
    int fast_mode;
    double filter_words;
} fij_opts_t;

int
fij(fij_arg *fij_args,
    int num_contigs1,
    contig_list_t *contig_array1,
    int num_contigs2,
    contig_list_t *contig_array2);

int do_it_fij(fij_arg *fij_args, char *seq, int seq_len,
	      int gap_open, int gap_extend,
	      int compare_mode, HashTable *lib_hash, HashTable *links,
	      Contig_parms *contig_list1, int number_of_contigs1,
	      Contig_parms *contig_list2, int number_of_contigs2,
	      int num_shared);

void buffij(tg_rec c1, int pos1, int end1, 
	    tg_rec c2, int pos2, int end2,
	    int len, int score, double percent);

#endif /* _FIJ_H_ */
