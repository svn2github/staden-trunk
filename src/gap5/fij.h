#ifndef _FIJ_H_
#define _FIJ_H_

#include "io_utils.h"
#include "consen.h"
#include "hash_lib.h"
#define COMPARE_ALL 0
#define COMPARE_ONE_BY_ONE 1

int
fij(GapIO *io,
    int mask,
    int min_overlap,
    double max_percent_mismatch,
    int word_len,
    double max_prob,
    int min_match,
    int band,
    int window_len,
    int max_unknown,
    double min_conf,
    int use_conf,
    int use_hidden,
    int max_alignment,
    int fast_mode,
    double filter_words,
    int num_contigs1,
    contig_list_t *contig_array1,
    int num_contigs2,
    contig_list_t *contig_array2);

int do_it_fij (char seq[], int seq_len,
	       int word_len, int min_overlap,
	       double max_percent_mismatch, int compare_mode,
	       int band, int gap_open, int gap_extend, double max_prob,
	       int min_match, int max_alignment, int fast_mode,
	       double filter_words,
	       Contig_parms *contig_list1, int number_of_contigs1,
	       Contig_parms *contig_list2, int number_of_contigs2,
	       int num_shared);

#endif /* _FIJ_H_ */
