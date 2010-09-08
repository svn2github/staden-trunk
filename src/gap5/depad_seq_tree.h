#ifndef _DEPAD_SEQ_TREE_H_
#define _DEPAD_SEQ_TREE_H_

#ifndef TEST_MAIN
#    include <tg_gio.h>
#endif

#include <tree.h>

typedef struct pad_count {
    RB_ENTRY(pad_count) link;
    int pos;
    int ppos;
    int count;
} pad_count;

typedef RB_HEAD(PAD_COUNT, pad_count) pad_count_t;

RB_PROTOTYPE(PAD_COUNT, pad_count, link, pad_count_cmp);

#ifndef TEST_MAIN
pad_count_t *depad_consensus(GapIO *io, tg_rec crec);
#endif

char *repad_seq_tree(char *seq, pad_count_t *tree);
void depad_seq_tree_free(pad_count_t *tree);
int get_padded_coord(pad_count_t *tree, int unpadded);
int padtree_pad_at(pad_count_t *tree, int pos);

#endif /* _DEPAD_SEQ_TREE_H_ */
