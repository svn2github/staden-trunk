#ifndef _TG_INDEX_H_
#define _TG_INDEX_H_

#include <stdio.h>

typedef struct {
    char name[L_tmpnam];
    FILE *fp;
} bttmp_t;

typedef struct {
    int append;
    int no_tree;
    int merge_contigs;
    int fmt;
    char *out_fn;
    int pair_reads;
    int min_bin_size;
    int fast_mode;
    bttmp_t *tmp;
    int reserved_seqs;
} tg_args;

#endif /* _TG_INDEX_H_ */
