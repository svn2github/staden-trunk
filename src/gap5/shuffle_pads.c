/* Copyright Genome Research Limited (GRL). All rights reserved */

/*
 * TODO:
 * 
 * ----------------------------------------------------------------------
 * Allow sequences to move. Often we have alignments ending or starting like:
 *
 *  ACGGG
 *  AC*GGGTA
 *  AC*GGGTA
 *  ACGGGGTA
 *  AC*GGGTA
 *
 * The first sequence is reinforcing there being 4 Gs, but it actually only
 * has 3. The problem is that it cannot insert the pad as that changes the
 * sequence length.
 *
 * Solution. Let X be a specific base call (one of A, C, G, T, but always the
 * same member of that set).  X(n) is a run of 1 or more X.
 * Find cases where we have sequence*X(n) or X(n)*sequence.
 * Check the malign vector at the * to see if it also contains X. If so
 * trim * and X(n).
 *
 * ----------------------------------------------------------------------
 * Investigate 454 rate of miscall vs indel. Seems maybe we need to mirror
 * this and get the pad penalty much lower than a mismatch.
 *
 * ----------------------------------------------------------------------
 * Investigate the issue of reassigning confidence values during runs of
 * bases for 454 data. AGGGT may have confidence X 40 30 10 X if in the +ve
 * direction but X 10 30 40 X if in the -ve direction. After pad shuffling
 * we need to have the pads aligned against the low quality bases and not
 * the high quality ones. This means several things:
 *
 * 1. Reording the confidence of base-calls in a run
 * 2. Making sure the pads always end up at the same end (needs another
 *    algorithm after this one to do that).
 * 3. The pad confidence value cannot now just be the average of the two
 *    surrounding bases. Maybe the preceeding base confidence works.
 *
 * ----------------------------------------------------------------------
 * Remove the O(N^2) complexity code and make this O(N). The most obvious
 * case is inserting and deleting into the consensus. Currently this does
 * large scale memmoves over the entire contig, but in theory we can do
 * little more than local updates if we have the following:
 *
 * Consensus base structure:
 *     next/prev points
 *     counts[6]
 *     scores[6]
 *     orig_position
 *
 * Sequence fragment structure:
 *     Consensus base pointer (for left-most end)
 *     distance from last (relative offset rather than absolute)
 *     length
 *     sequence
 *
 * Then consensus pad insertion/deletion is just a matter of updating links.
 * Q: How do we handle removal of a base to which a sequence fragment is
 * pointing? I guess we need a list of fragments pointed to by the consensus
 * base (which is like the up/down pointers in ReAligner). Keeping this up to
 * date is a bit tricky.
 *
 * ----------------------------------------------------------------------
 */


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include "tg_gio.h"
#include "align.h"
#include "dna_utils.h"
#include "align_lib.h"
#include "text_output.h"
#include "shuffle_pads.h"

typedef struct {
    int pos;
    int size; /* +ve or -ve for ins/del */
} con_indel_t;

void print_malign(MALIGN *malign);
void print_moverlap(MALIGN *malign, MOVERLAP *o, int offset);

/*
 * Insert 'size' pads into a contig at position 'pos'.
 */
void malign_padcon(MALIGN *malign, int pos, int size, Array indels) {
    CONTIGL *cl = malign->contigl;
    con_indel_t *id;

    id = ARRP(con_indel_t, indels, ArrayMax(indels));
    id->pos = pos;
    id->size = size;

    for (; cl; cl = cl->next) {
	/* We do one of three things: nothing, insert, or shift */
	/* Nothing: */
	if (cl->mseg->offset+cl->mseg->length-1 < pos)
	    continue;

	/* Shift right: */
	if (cl->mseg->offset >= pos) {
	    cl->mseg->offset += size;
	    continue;
	}

	/* Insert */
	cl->mseg->length += size;
	cl->mseg->seq = (char *)realloc(cl->mseg->seq, cl->mseg->length+1);
	memmove(&cl->mseg->seq[pos - cl->mseg->offset + size],
		&cl->mseg->seq[pos - cl->mseg->offset],
		cl->mseg->length-size - (pos - cl->mseg->offset));
	memset(&cl->mseg->seq[pos - cl->mseg->offset], '*', size);
	cl->mseg->seq[cl->mseg->length] = 0;
    }

    malign_insert_scores(malign, pos, size);
}

/*
 * Returns the number of consensus pads added or -1 for error.
 */
int edit_mseqs(MALIGN *malign, CONTIGL *cl, MOVERLAP *o, int cons_pos,
	       Array indels) {
    int i, npads, poso;
    char *cp;

    /* Cons vector */
    npads = 0;
    for (poso = i = 0; i < o->s1_len; i++) {
	if (o->S1[i] < 0) {
	    /*printf("S1:Ins %d pads at pos %d+%d=%d\n",
	      -o->S1[i], poso, cons_pos, poso+cons_pos);*/
	    malign_padcon(malign, poso+cons_pos+npads, -o->S1[i], indels);
	    npads += -o->S1[i];
	} else {
	    poso += o->S1[i];
	}
    }

    /* sequence */
    /* Trim leading pads */
    cp = o->seq2_out;
    while(*cp == '.') {
	cp++;
	cl->mseg->offset++;
    }

    xfree(cl->mseg->seq);
    cl->mseg->seq = strdup(cp);
    for (cp = cl->mseg->seq; *cp; cp++) {
	if (*cp == '.')
	    *cp = '*';
    }

    /* Back off trailing pads */
    while (*(cp-1) == '*')
	cp--;

    cl->mseg->length = cp-cl->mseg->seq;

    /*
    printf("cl->mseg->seq=%.*s (len %d)\n",
	   cl->mseg->length, cl->mseg->seq, cl->mseg->length);
    */

    return npads;
}

/*
 * Realigning the sequences may change their start positions and hence break
 * the sorted-on-position property.
 * We make sure this is maintained here.
 */
static void resort_contigl(MALIGN *malign) {
    CONTIGL *lastl, *cl;
    int swapped;

    /* The list is very close to being sorted, so this is not too inefficient */
    do {
	swapped = 0;
	lastl = NULL;
	for (cl = malign->contigl; cl && cl->next; lastl = cl, cl = cl->next) {
	    if (cl->mseg->offset > cl->next->mseg->offset) {
		CONTIGL *n = cl->next;

		/* Swap cl & cl->next */
		if (lastl)
		    lastl->next = n;
		else
		    malign->contigl = n;
		cl->next = n->next;
		n->next = cl;
		cl = lastl ? lastl : malign->contigl;
		swapped = 1;
	    }
	}
    } while (swapped);
}

typedef struct cl_list {
    CONTIGL *cl;
    int offset;
    struct cl_list *next;
} cl_list;

static void remove_pads(GapIO *io, MALIGN *malign, contig_t *c) {
    int i, removed = 0;
    CONTIGL *cl = malign->contigl;
    cl_list *head = NULL, *c2, *last, *next;
    int npads, depth;

    for (i = 0; i < malign->length; i++) {
	/* Add new seqs to the depth array as we meet them */
	while (cl && cl->mseg->offset == i) {
	    c2 = (cl_list *)xmalloc(sizeof(cl_list));
	    c2->next = head;
	    c2->offset = 0;
	    c2->cl = cl;
	    head = c2;
	    cl = cl->next;
	}

	/* Remove any sequences we've now passed, also counting pads */
	npads = 0;
	depth = 0;
	last = NULL;
	for (c2 = head; c2; c2 = next) {
	    next = c2->next;
	    if (c2->offset == c2->cl->mseg->length) {
		if (last)
		    last->next = c2->next;
		else
		    head = c2->next;
		xfree(c2);
		continue;
	    }
	    last = c2;
	    if (c2->cl->mseg->seq[c2->offset++] == '*')
		npads++;
	    depth++;
	}

	if (npads != depth || depth == 0)
	    continue;

	/* We have a column of pads, so remove it */
	//printf("Remove pad at %d\n", i+1-removed);
	contig_delete_base(io, &c, i+1-removed);

	removed++;
    }

    malign_recalc_scores(malign, 0, malign->length-1);
}

/*
 * Iterates through all sequences in a contig realigning them against the
 * consensus vector.
 *
 * It then adds the newly aligned sequence back into the consensus, editing the
 * sequence and tag positions/lengths too.
 * To do this we may need to shuffle the start position of sequences
 * downstream, and hence also move consensus tags.
 */
MALIGN *realign_seqs(int contig, MALIGN *malign, int band, Array indels) {
    CONTIGL *lastl = NULL, *contigl;
    int nsegs, r;
    int old_start, old_end, new_start, new_end;

    for (contigl = malign->contigl, nsegs = 0; contigl; nsegs++)
	contigl = contigl->next;

    /* Loop through all sequences in the contig */
    contigl = malign->contigl;
    while (contigl) {
	int len;
	MOVERLAP *o;
	ALIGN_PARAMS *p;
	int cons_pos;
	int npads;

	/* Obtain a depadded copy of this mseg */
	len = contigl->mseg->length;


	/* Remove sequence from malign */
	malign_remove_contigl(malign, lastl, contigl);


	/* Align sequence to malign */
	p = create_align_params();
	set_align_params (p,
			  band,
			  8, /*gap_open*/
			  8, /*gap_extend*/
			  /* EDGE_GAPS_COUNT, */
			  EDGE_GAPS_ZEROX | BEST_EDGE_TRACE,
			  RETURN_EDIT_BUFFERS | RETURN_SEQ |
			  RETURN_NEW_PADS,
			  0,  /*seq1_start*/
			  0,  /*seq2_start*/
			  0,  /*old pad sym*/
			  0,  /*new pad sym*/
			  0   /*set_job*/);

	o = create_moverlap();
	init_moverlap(o, malign, contigl->mseg->seq, malign->length, len);

	cons_pos = contigl->mseg->offset;
	o->malign_len = malign->length - cons_pos;

	/* 3 bases overhang to the right */
	if (o->malign_len > contigl->mseg->length+band/2+1)
	    o->malign_len = contigl->mseg->length+band/2+1;

	/* And 3 to the left */
	if (cons_pos > band/2+1) {
	    cons_pos -= band/2+1;
	    o->malign_len += band/2+1;
	    contigl->mseg->offset -= band/2+1;
	} else {
	    o->malign_len += cons_pos;
	    contigl->mseg->offset -= cons_pos;
	    cons_pos = 0;
	}

	{
	    char *old_cons   = malign->consensus;
	    int **old_scores = malign->scores;
	    int **old_counts = malign->counts;

	    malign->consensus += cons_pos;
	    malign->counts    += cons_pos;
	    malign->scores    += cons_pos;

	    /* fixed_malign(o, p); */
	    r = realigner_malign(o, p); /* o->score = alignment score */
	    
	    /*
	    if (!r)
		print_moverlap(malign, o, cons_pos);
	    else
		puts("FAILED");
	    */

	    malign->consensus = old_cons;
	    malign->counts    = old_counts;
	    malign->scores    = old_scores;
	}

	/* Edit the sequence with the alignment */
	old_start = contigl->mseg->offset;
	old_end   = contigl->mseg->offset + contigl->mseg->length-1;
	if (r == 0 && o->S1)
	    npads = edit_mseqs(malign, contigl, o, cons_pos, indels);
	else
	    npads = 0;
	new_start = contigl->mseg->offset;
	new_end   = contigl->mseg->offset + contigl->mseg->length-1;


	/* Put sequence back */
	malign_add_contigl(malign, lastl, contigl);


	/* Update the malign structure */
	if (npads > 0) {
	    malign_recalc_scores(malign,
				 MIN(old_start, new_start),
				 MAX(old_end, new_end));
	}
	    
	/* TODO:
	 *
	 * X Realloc malign->consensus / malign->score
	 * X Move malign->consensus from here to end right by npads.
	 * X Move malign->score      " ...
	 * X Update malign->length
	 * X Recompute consensus and score over the length of this reading.
	 *
	 * If contigl was doubly linked (sorted on left and right ends
	 * separately) then we could chain left/right to only update
	 * those readings which overlap this region. For now we can
	 * just chain from left each time.  Not optimal (O(N^2) for
	 * full realignment method then) but workable perhaps.
	 *
	 * See get_malign_counts, scale_malign_scores and get_malign_consensus
	 */


	/*
	 * Check if the short-cut method gives the same result as rebuilding
	 * from scratch.
	 */
#if 0
	{
	    int i, j;
	    MALIGN *copy;
	    copy = contigl_to_malign(malign->contigl, -4, -4);

	    for (i = 0; i < copy->length; i++) {
		for (j = 0; j < copy->charset_size+2; j++) {
		    if (copy->scores[i][j] != malign->scores[i][j]) {
			printf("[%d][%d] = %d (should be %d)\n",
			       i, j,
			       malign->scores[i][j],
			       copy->scores[i][j]);
		    }
		}
	    }
	    copy->contigl = NULL;
	    destroy_malign(copy, 0);
	}
#endif

	destroy_moverlap(o);
	destroy_alignment_params(p); 

	lastl = contigl;
	contigl = contigl->next;
    }

    resort_contigl(malign);

    return malign;
}

/**
 * Builds and returns MALIGN from a Gap5 IO handle for the contig 'cnum'.
 */
MALIGN *build_malign(GapIO *io, tg_rec cnum /*, int start, int end */) {
    CONTIGL *contig, *first_contig = NULL, *last_contig = NULL;
    int i, j;
    contig_iterator *citer;
    rangec_t *r;

    /* Generate contigl linked list */
    citer = contig_iter_new(io, cnum, 1, CITER_FIRST, CITER_CSTART, CITER_CEND);
    // citer = contig_iter_new(io, cnum, 0, CITER_FIRST, start, end);
			   
    while ((r = contig_iter_next(io, citer))) {
	seq_t *s, *sorig;
	char *seq;
	int len;

	assert((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ);

	contig = create_contig_link();
	contig->id = r->rec;
	contig->mseg = create_mseg();

	sorig = s = cache_search(io, GT_Seq, r->rec);
	if ((s->len < 0) ^ r->comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	assert(s->right >= s->left);
	len = s->right - s->left + 1;
	if (NULL == (seq = malloc(len+1)))
	    return NULL;

	for (j = 0, i = s->left-1; i < s->right; i++, j++) {
	    /* Protect against the sequence containing "."; our pad sym */
	    if (s->seq[i] == '.')
		seq[j] = 'N';
	    else
		seq[j] = s->seq[i];
	}
	seq[j] = 0;

	init_mseg(contig->mseg, seq, len, r->start-1 + s->left-1);
	contig->mseg->comp = (s != sorig);

	if (last_contig) {
	    last_contig->next = contig;
	} else {
	    first_contig = contig;
	}
	last_contig = contig;

	if (s != sorig)
	    free(s);
    }
    contig_iter_del(citer);

    /* for 454 data -6 to -10 seem to work fine */
    return contigl_to_malign(first_contig, -7, -7);
}

#define LLEN 80
struct clist {
    char *seq;
    int len;
    char line[LLEN];
};

void print_malign(MALIGN *malign) {
    int i, j;
    struct clist *depth = NULL;
    int ndepth = 0;
    CONTIGL *cl = malign->contigl;

    puts("MALIGN OUTPUT");
    for (i = 0; i < malign->length; i++) {
	/* Maintain a list of CONTIGLs covering this point */

	/* ... adding new items to the list */
	while (cl && cl->mseg->offset <= i) {
	    ndepth++;
	    /* runaway loops completely kills deskpros */
	    if (ndepth > 100000)
		abort();
	    depth = (struct clist *)realloc(depth, ndepth * sizeof(*depth));
	    depth[ndepth-1].seq = cl->mseg->seq;
	    *depth[ndepth-1].seq = tolower(*depth[ndepth-1].seq);
	    depth[ndepth-1].seq[cl->mseg->length-1] =
		tolower(depth[ndepth-1].seq[cl->mseg->length-1]);
	    depth[ndepth-1].len = cl->mseg->length;
	    memset(depth[ndepth-1].line, ' ', LLEN);
	    cl = cl->next;
	}

	for (j = 0; j < ndepth; j++) {
	    depth[j].line[i%LLEN] = (depth[j].seq) ? *depth[j].seq++ : ' ';
	    if (depth[j].len > 0 && --depth[j].len == 0) {
		depth[j].seq = NULL;
	    }
	}

	/* Print line, and remove items from depth as and when needed */
	if (i%LLEN == LLEN-1) {
	    for (j = LLEN * (int)(i/LLEN); j < i; j+=10)
		printf("%10d", j+10);
	    printf("\n");
	    for (j = 0; j < ndepth; j++) {
		printf("%.*s\n", LLEN, depth[j].line);
		if (!depth[j].seq) {
		    memmove(&depth[j], &depth[j+1],
			    (ndepth-(j+1)) * sizeof(depth[j]));
		    ndepth--;
		    j--;
		}
	    }
	    printf("\n");
	}
    }

    /* Print remainder of lines */
    if ((i-1)%LLEN != LLEN-1) {
	for (j = LLEN * (int)(i/LLEN); j < i; j+=10)
	    printf("%10d", j+10);
	printf("\n");
	for (j = 0; j < ndepth; j++) {
	    printf("%.*s\n", i - LLEN * (int)(i/LLEN), depth[j].line);
	}
	printf("\n");
    }

    free(depth);
}

void print_moverlap(MALIGN *malign, MOVERLAP *o, int offset) {
    int i, j;
    struct clist *depth = NULL;
    int ndepth = 0;
    CONTIGL *cl = malign->contigl;
    int s1op = 0, s2op = 0;
    int *S1 = o->S1;
    int *S2 = o->S2;
    char *seq = o->seq2;
    int cins = 0;

    for (i = offset; i < malign->length+offset; i++) {
	/* Maintain a list of CONTIGLs covering this point */

	/* ... adding new items to the list */
	for (; cl && cl->mseg->offset+cins <= i; cl = cl->next) {
	    if (cl->mseg->offset+cins + cl->mseg->length-1 < i)
		continue;
	    ndepth++;
	    /* runaway loops completely kills deskpros */
	    if (ndepth > 1000)
		abort();
	    depth = (struct clist *)realloc(depth, ndepth * sizeof(*depth));
	    depth[ndepth-1].seq = cl->mseg->seq + i-(cl->mseg->offset+cins);
	    depth[ndepth-1].len = cl->mseg->length - (i-(cl->mseg->offset+cins));
	    memset(depth[ndepth-1].line, ' ', LLEN);
	}

	if (!s1op) {
	    s1op = *S1++;
	    if (S1-o->S1 > o->s1_len)
		break;
	}
	if (!s2op) {
	    s2op = *S2++;
	    if (S2-o->S2 > o->s2_len)
		break;
	}

	printf("%4d: ", i);

	if (s1op < 0) {
	    /* Ins to consensus */
	    s1op++;
	    printf("%c\n", *seq++);
	    cins++;
	    continue;
	} else if (s2op > 0) {
	    /* Match/mismatch */
	    printf("%c ", *seq++);
	    s2op--;
	} else if (s2op < 0) {
	    /* Ins to sequence */
	    printf("  ");
	    s2op++;
	}

	s1op--;
	for (j = 0; j < ndepth; j++) {
	    printf("%c", *depth[j].seq++);
	    if (--depth[j].len == 0) {
		depth[j].seq = NULL;
		memmove(&depth[j], &depth[j+1],
			(ndepth-(j+1)) * sizeof(depth[j]));
		ndepth--;
		j--;
	    }
	}
	printf("\n");
    }

    free(depth);
}

#include <ctype.h>
int malign_diffs(MALIGN *malign, int *tot) {
    CONTIGL *cl;
    int diff_count = 0, tot_count = 0;

    /* printf("%.*s\n", malign->length, malign->consensus); */
    for (cl = malign->contigl; cl; cl = cl->next) {
	int i;

	/*
	for (i = 0; i < cl->mseg->length; i++, end_gaps++) {
	    if (cl->mseg->seq[i] != '*')
		break;
	}
	for (i = cl->mseg->length-1; i >= 0; i--, end_gaps++) {
	    if (cl->mseg->seq[i] != '*')
		break;
	}
	*/

	for (i = 0; i < cl->mseg->length; i++) {
	    char c = toupper(malign->consensus[i+cl->mseg->offset]);
	    char s = toupper(cl->mseg->seq[i]);
	    if (c == '-')
		c = '*';

	    /*printf("%c", c==s ? '.' : s);*/
	    if (s != c)
		diff_count++;
	    tot_count++;
	}
    }

    if (tot)
	*tot = tot_count;
    return diff_count;
}

#if 0
static void update_consensus_tags(GapIO *io, int cnum, MALIGN *malign) {
    int i, last = 0;
    for (i = 0; i < malign->length; i++) {
	int p = malign->orig_pos[i];
	if (p == 0) {
	    /* Insertion */
	    shift_contig_tags(io, cnum, i+1, +1);
	} else {
	    if (p-last != 1) {
		/* Deletion */
		shift_contig_tags(io, cnum, i+1, 1-(p-last));
	    }
	    last = p;
	}
    }
}
#endif

/*
 * Moves tags on sequence 'srec' in contig 'crec' right by one if they
 * start at or beyond pos, and extends by one if they overlap pos.
 *
 * For optimisation purposes, we already know the sequence spans positions
 * start..end, so we use this for filtering our tag search.
 * However there is much else that can be optimised. We repeatedly query
 * contigs, bins, etc and almost certainly perform many contig iters over
 * the same region. (Functionality first, efficiency later.)
 */
static void tag_shift_for_insert(GapIO *io, tg_rec crec, tg_rec srec,
				 int start, int end, int pos) {
    contig_iterator *ci;
    rangec_t *r;
    contig_t *c = cache_search(io, GT_Contig, crec);;

    //printf("> tag in seq %"PRIrec" at %d+%d\n", srec, start, pos);

    ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST | CITER_ISTART,
				 start+pos, end, GRANGE_FLAG_ISANNO);
    if (!ci)
	return;

    while ((r = contig_iter_next(io, ci))) {
	range_t r2, *r_out;
	anno_ele_t *a;
	bin_index_t *bin;

	if (r->pair_rec != srec)
	    continue;

	bin_remove_item(io, &c, GT_AnnoEle, r->rec);
	r2.start    = (r->start >= start+pos) ? r->start+1 : r->start;
	r2.end      = r->end+1;
	r2.mqual    = r->mqual;
	r2.rec      = r->rec;
	r2.pair_rec = r->pair_rec;
	r2.flags    = r->flags;
	bin = bin_add_range(io, &c, &r2, &r_out, NULL, 0);

	a = cache_search(io, GT_AnnoEle, r->rec);
	if (a->bin != bin->rec ||
	    a->idx != r_out - ArrayBase(range_t, bin->rng)) {
	    /* Annotation moved bins */
	    a = cache_rw(io, a);
	    a->bin = bin->rec;
	    a->idx = r_out - ArrayBase(range_t, bin->rng);
	}
    }

    contig_iter_del(ci);
}

static void tag_shift_for_delete(GapIO *io, tg_rec crec, tg_rec srec,
				 int start, int end, int pos) {
    contig_iterator *ci;
    rangec_t *r;
    contig_t *c = cache_search(io, GT_Contig, crec);;

    //printf("< tag in seq %"PRIrec" at %d\n", srec, pos);

    ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST | CITER_ISTART,
				 start+pos, end, GRANGE_FLAG_ISANNO);
    if (!ci)
	return;

    while ((r = contig_iter_next(io, ci))) {
	range_t r2, *r_out;
	anno_ele_t *a;
	bin_index_t *bin;

	if (r->pair_rec != srec)
	    continue;

	bin_remove_item(io, &c, GT_AnnoEle, r->rec);
	r2.start    = (r->start > start+pos) ? r->start-1 : r->start;
	r2.end      = r->end-1;
	r2.mqual    = r->mqual;
	r2.rec      = r->rec;
	r2.pair_rec = r->pair_rec;
	r2.flags    = r->flags;

	if (r2.end < r2.start) {
	    /* Tag entirely removed now, it must have been on a pad */
	    a = cache_search(io, GT_AnnoEle, r->rec);
	    a = cache_rw(io, a);
	    cache_deallocate(io, a);
	    continue;
	}
	bin = bin_add_range(io, &c, &r2, &r_out, NULL, 0);

	a = cache_search(io, GT_AnnoEle, r->rec);
	if (a->bin != bin->rec ||
	    a->idx != r_out - ArrayBase(range_t, bin->rng)) {
	    /* Annotation moved bins */
	    a = cache_rw(io, a);
	    a->bin = bin->rec;
	    a->idx = r_out - ArrayBase(range_t, bin->rng);
	}
    }

    contig_iter_del(ci);
}

/*
 * Takes a multiple alignment and updates the on-disk data structures to
 * match. This needs to correct confidence values, original positions and
 * tags too.
 */
void update_io(GapIO *io, tg_rec cnum, MALIGN *malign, Array indels) {
    CONTIGL *cl;
    tg_rec rnum;
    range_t r, *r_out;
    bin_index_t *bin;
    contig_t *c = cache_search(io, GT_Contig, cnum);
    size_t i, nindel;

    cache_incr(io, c);

    /*
     * To minimise number of data modifications we use a three step approach.
     *
     * Step 1: insert columns of pads, shifting reads as appropriate.
     * Step 2: edit sequence alignments as required, possibly involving
     *         moving sequences and/or adding and removing pads.
     * Step 3: remove columns of entire pads.
     *
     * This means that when we introduce a column of pads we don't have
     * to make edits to every single read position down stream, and can
     * instead make use of the optimised recursive bin functions to do this
     * for us.
     */

    /* Step 1: make indels */
    nindel = ArrayMax(indels);
    for (i = 0; i < nindel; i++) {
	con_indel_t *id = arrp(con_indel_t, indels, i);
	int j;

	if (id->size > 0) {
	    for (j = 0; j < id->size; j++) {
		contig_insert_base(io, &c, id->pos+1, '*', -1);
	    }
	} else {
	    for (j = 0; j < -id->size; j++) {
		contig_delete_base(io, &c, id->pos+1);
	    }
	}
    }

    /* Step 2: edit alignments */
    for (cl = malign->contigl; cl; cl = cl->next) {
	seq_t *s, *sorig;
	int len, update_range;
	int shift;

	rnum = cl->id;
	
	sorig = cache_search(io, GT_Seq, rnum);
	s = dup_seq(sorig);
	if (cl->mseg->comp)
	    complement_seq_t(s);

	len = s->right - s->left + 1;

	/* Check if sequence has changed. If so assign a new one */
	if (cl->mseg->length != len ||
	    memcmp(s->seq + s->left-1, cl->mseg->seq, cl->mseg->length) != 0) {
	    int newlen = s->left-1 + ABS(s->len) - s->right + cl->mseg->length;
	    int i, j, np;
	    char *newseq  = (char *)malloc(newlen+1);
	    int1 *newconf = (int1 *)malloc(newlen+1);

	    /* Build new seq/conf arrays */
	    memcpy(newseq,  s->seq,  s->left-1);
	    memcpy(newconf, s->conf, s->left-1);

	    memcpy(&newseq[s->left-1], cl->mseg->seq, cl->mseg->length);

	    /*
	     * Step through both old and new sequences working out how
	     * they differ. This will (*should*) be entire pad movements.
	     * i = index to old seq
	     * j = index to new seq
	     * np = number of pads added minus removed from old seq.
	     */
	    np = 0;
	    for (i =j =s->left-1;
		 i < ABS(s->len) && j < s->left-1 + cl->mseg->length;
		 ) {
		/* Bases match */
		if (toupper(newseq[j]) == toupper(s->seq[i]) ||
		    (s->seq[i] == '.' && newseq[j] == 'N')) {
		    if (isupper(s->seq[i]))
			newseq[j] = toupper(newseq[j]);
		    else
			newseq[j] = tolower(newseq[j]);
		    newconf[j] = s->conf[i];
		    i++, j++;
		    continue;
		}

		/* Pad removed */
		if (s->seq[i] == '*') {
		    i++;
		    tag_shift_for_delete(io, cnum, rnum, cl->mseg->offset,
					 cl->mseg->length, i+np--);
		    /*
		    if (io_length(io, rnum) < 0) {
			tag_shift_for_delete(io, rnum, r.length - i + 1);
		    } else {
			tag_shift_for_delete(io, rnum, i+np--);
		    }
		    */
		    continue;
		}

		/* Pad created */
		if (newseq[j] == '*') {
		    int k;
		    int ql = 0, qr = 0;
		    for (k = i-1; k >= 0; k--) {
			if (s->seq[k] != '*') {
			    ql = s->conf[k];
			    break;
			}
		    }
		    for (k = i+1; k < s->right; k++) {
			if (s->seq[k] != '*') {
			    qr = s->conf[k];
			    break;
			}
		    }
		    newconf[j] = MIN(ql, qr); /* min conf of neighbours */
		    j++;
		    tag_shift_for_insert(io, cnum, rnum, cl->mseg->offset,
					 cl->mseg->length, i+ ++np);
		    /*
		    if (io_length(io, rnum) < 0) {
			tag_shift_for_insert(io, rnum, r.length - i + 1);
		    } else {
			tag_shift_for_insert(io, rnum, i+ ++np);
		    }
		    */
		    continue;
		}

		fprintf(stderr, "Alignment introduced non-pad character");
		abort();
	    }

	    /* Pads previously at the end of the reading & now removed */
	    while (i < s->right) {
		if (s->seq[i] == '*') {
		    i++;
		    tag_shift_for_delete(io, cnum, rnum, cl->mseg->offset,
					 cl->mseg->length, i+np--);
		    /*
		    if (io_length(io, rnum) < 0) {
			tag_shift_for_delete(io, rnum, r.length - i + 1);
		    } else {
			tag_shift_for_delete(io, rnum, i+np--);
		    }
		    */
		} else {
		    /* Error: clipped data that wasn't a pad */
		    abort();
		}
	    }

	    /* Should only be pads remaining in newseq, if anything */
	    s->right = j;
	    for (; j < s->left-1 + cl->mseg->length; j++) {
		if (newseq[j] != '*') {
		    fprintf(stderr, "Alignment introduced non-pad character");
		    abort();
		}
		newconf[j] = 0;
	    }

	    /* Append on the right hand cutoff data */
	    for (; i < ABS(s->len); i++, j++) {
		newseq[j]  = s->seq[i];
		newconf[j] = s->conf[i];
	    }
	    if (j != newlen) {
		abort();
	    }

	    /* Write it back out */
	    /* Copy newseq/newconf into seq_t */

	    s->seq = newseq;
	    s->conf = newconf;
	    update_range = 0;
	    if (ABS(s->len) != j) {
		/* Length change implies updating the range array too */
		s->len = s->len >= 0 ? j : -j;
		update_range = 1;
	    }

	    if (cl->mseg->comp)
		complement_seq_t(s);

	    sorig = cache_rw(io, sorig);
	    memcpy(sorig, s, sizeof(seq_t));

	    if (update_range)
		sorig = cache_item_resize(sorig, sizeof(*sorig) +
					  sequence_extra_len(sorig));

	    sequence_reset_ptr(sorig);

	    if (s->name)
		memcpy(sorig->name,       s->name,       s->name_len+1);
	    if (s->trace_name)
		memcpy(sorig->trace_name, s->trace_name, s->trace_name_len+1);
	    if (s->alignment)
		memcpy(sorig->alignment,  s->alignment,  s->alignment_len+1);
	    memcpy(sorig->seq,  s->seq,  ABS(s->len));
	    memcpy(sorig->conf, s->conf, ABS(s->len));

	    xfree(newconf);
	    xfree(newseq);
	}

	{
	    int st, en, or;
	    sequence_get_position(io, s->rec, NULL, &st, &en, &or);
	    if (or ^ (sorig->len < 0)) {
		shift = ABS(sorig->len) - sorig->right;
	    } else {
		shift = sorig->left-1;
	    }
	    st += shift;
	    if (st != cl->mseg->offset+1) {
		update_range = 1;
	    }
	}

	free(s);

	if (update_range) {
	    /* Get old range and pair data */
	    s = sorig;
	    bin = cache_search(io, GT_Bin, s->bin);
	    r = *arrp(range_t, bin->rng, s->bin_index);
	    assert(r.rec == s->rec);

	    /* Update range, tedious and slow way */
	    bin_remove_item(io, &c, GT_Seq, s->rec);
	    r.start = cl->mseg->offset + 1 - shift;
	    r.end   = r.start + ABS(s->len) - 1;
	    bin = bin_add_range(io, &c, &r, &r_out, NULL, 0);

	    /* Check if the new bin has a different complemented status too */
	    if (s->bin != bin->rec) {
		int old_comp = bin_get_orient(io, s->bin);
		int new_comp = bin_get_orient(io, bin->rec);

		if (new_comp != old_comp) {
		    int tmp;
		    s->len *= -1;
		    tmp = s->left;
		    s->left  = ABS(s->len) - (s->right-1);
		    s->right = ABS(s->len) - (tmp-1);
		}
	    }
	
	    /* Update seq bin & bin_index fields */
	    s = cache_rw(io, s);
	    s->bin = bin->rec;
	    s->bin_index = r_out - ArrayBase(range_t, bin->rng);
	}
    }

    /* Step 3 (remove pad columns) done in calling function. */

    cache_decr(io, c);
}

static int isort(const void *vp1, const void *vp2) {
    return *(const int *)vp2 - *(const int *)vp1;
}


#if 0
/*
 * Specifically for 454 data this reassigns confidence values to bases in
 * a run of the same base type.
 * It also reassigns confidence values of pads to be the minimum confidence
 * of the surrounding base call.
 */
void reassign_confidence_values(GapIO *io, int cnum) {
    GContigs c;
    GReadings r;
    int rnum;
    int scores[1000]; /* FIXME: check if we overflow! */

    contig_read(io, cnum, c);
    for (rnum = c.left; rnum; rnum = r.right) {
	char last = 0;
	char *seq;
	int1 *conf;
	int i, j, k;
	int cl, cr;

	gel_read(io, rnum, r);
	seq = TextAllocRead(io, r.sequence);
	conf = DataAllocRead(io, r.confidence, 1);

	/* Rearrange confidence in runs of bases */
	for (i = 0; i < r.length; i++) {
	    /* Find first non-pad, at 'i' */
	    while (i < r.length && seq[i] == '*')
		i++;
	    k = 0;
	    scores[k++] = conf[i];
	    last = seq[i];

	    /* Count how many there are. First diff base at 'j' */
	    j = i+1;
	    while (j < r.length && (seq[j] == '*' || seq[j] == last)) {
		if (seq[j] != '*')
		    scores[k++] = conf[j];
		j++;
	    }
		   
	    if (k != 1) {
		/* We have a run of k items (from >='i' and <'j') */
		qsort(scores, k, sizeof(*scores), isort);
		
		/* Reassign */
		j = i; k = 0;
		while (j < r.length && (seq[j] == '*' || seq[j] == last)) {
		    if (seq[j] != '*')
			conf[j] = scores[k++];
		    j++;
		}
	    }

	    i = j-1;
	}

	/* Reassign confidences to pads */
	cl = 0;
	for (i = 0; i < r.length; i++) {
	    if (seq[i] == '*') {
		for (j = i+1; j < r.length && seq[j] == '*'; j++)
		    ;
		cr = j < r.length ? conf[j] : 0;
		/* conf[i] = MIN(cl, cr); */
		conf[i] = (cl+cr)/2;
	    } else {
		cl = conf[i];
	    }
	}

	DataWrite(io, r.confidence, conf, r.length, 1);
	xfree(seq);
	xfree(conf);
    }
}
#endif


int shuffle_contigs_io(GapIO *io, int ncontigs, contig_list_t *contigs,
		       int band) {
    int i, start;
    Array indels;
    
    set_malign_lookup(5);
    /* set_alignment_matrix("/tmp/nuc_matrix", "ACGTURYMWSKDHVB-*"); */

    indels = ArrayCreate(sizeof(con_indel_t), 0);

    for (i = 0; i < ncontigs; i++) {
	tg_rec cnum = contigs[i].contig;
	int old_score, new_score, tot_score, orig_score;
	//for (start = 0; start < 1000000; start += 1000) {
	//  MALIGN *malign = build_malign(io, cnum, start, start + 1000);
	MALIGN *malign;
	int c_start, c_shift;
	bin_index_t *root_bin = NULL;

	vmessage("Shuffling pads for contig %s\n", get_contig_name(io, cnum));

	/*
	 * The shuffle pads code (malign) comes from gap4 and has lots of
	 * assumptions that the contig goes from base 1 to base N.
	 * Fixing these assumptions is a lot of work, so for now we will take
	 * the cheat route of moving the contig to ensure the assumption
	 * is valid.
	 */
	if (-1 == consensus_valid_range(io, cnum, &c_start, NULL)) {
	    verror(ERR_WARN, "shuffle_contigs_io",
		   "Failure in consensus_valid_range()");
	    return -1;
	}
	printf("Contig starts at base %d\n", c_start);
	c_shift = 1-c_start;
	if (c_shift != 0) {
	    contig_t *c = cache_search(io, GT_Contig, cnum);

	    if (!c)
		return -1;
	    if (!(root_bin = cache_search(io, GT_Bin, contig_get_bin(&c))))
		return -1;

	    printf("Shifting by %d\n", c_shift);
	    root_bin = cache_rw(io, root_bin);
	    root_bin->pos += c_shift;
	}

	malign = build_malign(io, cnum);
	resort_contigl(malign);

	ArrayMax(indels) = 0;
	orig_score = new_score = malign_diffs(malign, &tot_score);
	vmessage("Initial score %.2f%% mismatches (%d mismatches)\n",
		 (100.0 * orig_score)/tot_score, orig_score);
	UpdateTextOutput();
	do {
	    old_score = new_score;
	    malign = realign_seqs(cnum, malign, band, indels);
	    //print_malign(malign);
	    new_score = malign_diffs(malign, &tot_score);
	    vmessage("  Number of differences to consensus: %d\n", new_score);
	    UpdateTextOutput();
	} while (new_score < old_score);

	if (new_score < orig_score) {
	    //print_malign(malign);
	    update_io(io, cnum, malign, indels);
	} else {
	    vmessage("Could not reduce number of consensus differences.\n");
	}

	/* Remove pad columns */
	{
	    contig_t *c;
	    c = cache_search(io, GT_Contig, cnum);
	    cache_incr(io, c);
	    remove_pads(io, malign, c);
	    cache_decr(io, c);
	}

	destroy_malign(malign, 1);

	vmessage("Final score %.2f%% mismatches\n",
		 (100.0 * new_score)/tot_score);

	/*
	 * Sequences like
	 *   AGCT**GATGC
	 *             TGGATCGA
	 * can end up causing holes. We break the contig in this case to
	 * avoid minor database inconsistencies.
	 */
	puts("remove_contig_holes()");
	// remove_contig_holes(io, cnum);

	/* reassign_confidence_values(io, cnum); */
      //}

	/* Shift contig back */
	if (root_bin) {
	    root_bin->pos -= c_shift;
	}

	cache_flush(io);
    }

    ArrayDestroy(indels);

    return 0;
}
