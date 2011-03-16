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
    int data_type;
    int comp_mode;
    int repad;
    int store_unmapped;
    int sam_aux;
    int pair_queue;
    int store_refpos;
    int remove_dups;
} tg_args;

#define DATA_SEQ	1
#define DATA_QUAL	2
#define DATA_NAME	4
#define DATA_ANNO	8
#define DATA_ALL	15
#define DATA_BLANK	0x100 /* not even dummy seqs */

#endif /* _TG_INDEX_H_ */
