#ifndef _TG_INDEX_H_
#define _TG_INDEX_H_

typedef struct {
    int append;
    int no_tree;
    int merge_contigs;
    int fmt;
    char *out_fn;
    int pair_reads;
    int min_bin_size;
    int fast_mode;
} tg_args;

#endif /* _TG_INDEX_H_ */
