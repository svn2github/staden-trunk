#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "depad_seq_tree.h"
#include "consensus.h"

/* Comparison function for sorting the tree */
static int pad_count_cmp(struct pad_count *p1, struct pad_count *p2) {
    return p1->pos - p2->pos;
}

RB_GENERATE(PAD_COUNT, pad_count, link, pad_count_cmp);

/*
 * Given a gapped sequence 'seq', this both depads it and also produces
 * a tree with one node per gap removed. The nodes are sorted on the padded
 * position, but also contain the orginal padded coordinate.
 */
pad_count_t *depad_seq_tree(char *seq, int offset) {
    pad_count_t *tree = malloc(sizeof(*tree));
    struct pad_count *node, *n;
    int p;
    char *in, *out;
    int count = 0;

    RB_INIT(tree);

    for (p = 0, in = out = seq; *in; in++) {
	if (*in != '*') {
	    *out++ = *in;
	    p++;
	    continue;
	}

	/* Pad, so add to the tree */
	node = (struct pad_count *)malloc(sizeof(*node));
	node->pos = p + offset;
	node->ppos = p + offset + ++count;
	node->count = 1;
	n = RB_INSERT(PAD_COUNT, tree, node);
	if (n) {
	    /* Already existed, so increment ppos */
	    n->ppos++;
	    n->count++;
	    free(node);
	}
    }
    *out = 0;

    return tree;
}

/*
 * Deallocates memory taken up by a pad_count tree.
 */
void depad_seq_tree_free(pad_count_t *tree) {
    struct pad_count *node, *next;

    if (!tree)
	return;

    for (node = RB_MIN(PAD_COUNT, tree); node; node = next) {
	next = RB_NEXT(PAD_COUNT, tree, node);
	RB_REMOVE(PAD_COUNT, tree, node);
	free(node);
    }

    free(tree);
}

/*
 * Given an unpadded sequence and a pad tree this function puts back the
 * missing pads.
 *
 * It does this by allocating and returning a new sequence buffer of the
 * appropriate length.
 */
char *repad_seq_tree(char *seq, pad_count_t *tree) {
    struct pad_count *node = RB_MAX(PAD_COUNT, tree), *next;
    int npads, count;
    size_t slen = strlen(seq);
    char *pseq, *out;
    int last = 0, i;
    
    /* Allocate the buffer to the appropriate size */
    npads = node ? node->ppos - node->pos : 0;
    if (NULL == (pseq = malloc(slen+npads+1)))
	return NULL;

    /* Put the pads back in */
    out = pseq;
    count = 0;
    for (next = RB_MIN(PAD_COUNT, tree);
	 next;
	 next = RB_NEXT(PAD_COUNT, tree, next)) {
	memcpy(out, seq, next->pos - last);
	out += next->pos - last;
	npads = next->ppos - next->pos - count;
	for (i = 0; i < npads; i++)
	    *out++ = '*';
	count += npads;
	seq += next->pos - last;
	last = next->pos;
    }
    memcpy(out, seq, slen-last);
    out += slen-last;

    *out++ = 0;

    return pseq;
}

/*
 * Converts an unpadded coordinate to a padded one.
 */
int get_padded_coord(pad_count_t *tree, int unpadded) {
    struct pad_count *node, n;

    if (!tree)
	return unpadded;

    n.pos = unpadded+1;
    node = RB_NFIND(PAD_COUNT, tree, &n);
    if (!node) {
	/* Beyond the last node */
	node = RB_MAX(PAD_COUNT, tree);
	if (!node)
	    return unpadded;
    } else {
	/* Node is now >= unpadded, so get previous and count forward */
	node = RB_PREV(PAD_COUNT, tree, node);
	if (!node)
	    return unpadded;
    }

    return node->ppos + unpadded - node->pos;
}

int padtree_pad_at(pad_count_t *tree, int pos) {
    struct pad_count n, *curr;
    n.pos = pos;
    
    curr = RB_FIND(PAD_COUNT, tree, &n);
    return curr ? curr->count : 0;
}


void padtree_dump(pad_count_t *tree) {
    struct pad_count *next;

    for (next = RB_MIN(PAD_COUNT, tree);
	 next;
	 next = RB_NEXT(PAD_COUNT, tree, next)) {
	printf("Pad at %d padded %d count=%d\n",
	       next->pos, next->ppos, next->count);
    }
}


#ifndef TEST_MAIN
pad_count_t *depad_consensus(GapIO *io, tg_rec crec) {
    contig_t *c = cache_search(io, GT_Contig, crec);
    int first = c->start;
    int last = c->end;
    int len = last-first+1;
    pad_count_t *tree;
    char *cons;

    cons = (char *)xmalloc(len+1);
    calculate_consensus_simple(io, crec, first, last, cons, NULL);
    cons[len] = 0;

    tree = depad_seq_tree(cons, first);

    free(cons);

    return tree;
}
#endif

/* ----------------------------------------------------------------------
 * A simple test application for the depad_seq_tree function.
 */
#ifdef TEST_MAIN

#include <unistd.h>
#include <fcntl.h>

static char *load(int *lenp) {
    char *data = NULL;
    int dsize = 0;
    int dcurr = 0, len;

    do {
	if (dsize - dcurr < 8192) {
	    dsize = dsize ? dsize * 2 : 8192;
	    data = realloc(data, dsize);
	}

	len = read(0, data + dcurr, 8192);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
    }

    *lenp = dcurr;
    return data;
}

int main(void) {
    int len, i;
    char *data = load(&len), *data2, *data3;
    pad_count_t *tree;

    puts(data);
    data3 = strdup(data);
    tree = depad_seq_tree(data);
    puts("");
    
    padtree_dump(tree);

    puts(data);
    data2 = repad_seq_tree(data, tree);
    puts(data2);
    puts(data3);

#if 1
    puts("Benchmarking");
    for (i = 0; i < 1000000; i++) {
	struct pad_count n, *node;

	n.pos = rand() % 2000000;
	node = RB_NFIND(PAD_COUNT, tree, &n);
	//printf("%8d %8d %6d\n", n.pos, node->pos, node->count);
    }
    puts("done");
#endif

    return 0;
}
#endif
