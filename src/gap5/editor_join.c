/*
 * editor_join.c:
 * Contains functions for the Join Editor "align" and "join" buttons.
 */
#include <tg_gio.h>
#include <assert.h>
#include <string.h>
#include <assert.h>

#include "editor_view.h"
#include "editor_join.h"
#include "align_lib.h"
#include "hash_lib.h"
#include "align.h"
#include "dna_utils.h"

#include "dis_readings.h"

/* Add 'num' pads into the consensus for editor 'xx' at position 'pos'. */
static void add_pads(edview *xx, contig_t **ctg, int pos, int num)
{
    int i;

    if (num < 0)
	return;

    for (i = 0; i < num; i++)
	contig_insert_base(xx->io, ctg, pos, '*', -1);
}

/*
 * Converts a couple of padded sequences output from Rodger Staden's
 * alignment algorithms to an edit buffer in the same coding system as used
 * in the Myers & Miller alignment routines.
 * Basically a match is encoded as 0, N pads in sequence 1 are encoded as
 * +N, N pads in sequence 2 are encoded as -N.
 */
static int *
rsalign2myers(char *seq1, int len1, char *seq2, int len2, char pad_sym) {
    int *res = (int *)calloc(len1 + len2 + 1, sizeof(int));
    int i;
    int *S = res;
    int last_op = 0;

    for (i = 0; i < len1 && i < len2; i++) {
	if (seq1[i] != pad_sym && seq2[i] != pad_sym)
	    last_op = *S++ = 0;
	if (seq1[i] == pad_sym) {
	    if (last_op != 1)
		*S++ = 0;
	    *(S-1) += 1;
	    last_op = 1;
	}
	if (seq2[i] == pad_sym) {
	    if (last_op != 2)
		*S++ = 0;
	    *(S-1) -= 1;
	    last_op = 2;
	}
    }

    return res;
}

int align_contigs (OVERLAP *overlap, int fixed_left, int fixed_right) {

#define DO_BLOCKS 100
#define OK_MATCH 80
#define TOO_LONG_FOR_2ND_TRY 10000
    int ierr;
    ALIGN_PARAMS *params;
    Hash *h;
    int word_len, max_matches;
    int shortest_diagonal, longest_diagonal;
    int gap_open, gap_extend, band, edge_mode, job, min_match;
    int compare_method = 17;
    char *env;

    band = 1;
    gap_open = 12;
    gap_extend = 4;
    job = RETURN_SEQ | RETURN_NEW_PADS | RETURN_END_GAPS;

    edge_mode  = fixed_left  ? EDGE_GAPS_COUNT   : EDGE_GAPS_ZERO;
    edge_mode |= fixed_right ? FULL_LENGTH_TRACE : BEST_EDGE_TRACE;

    longest_diagonal = MAX(overlap->seq1_len,overlap->seq2_len);
    shortest_diagonal = MIN(overlap->seq1_len,overlap->seq2_len);
    min_match = MIN(20,(shortest_diagonal*0.1));

    if (NULL == (env = getenv("STADTABL"))) {
	verror(ERR_FATAL, "align_contigs",
	       "STADTABL environment variable is not set.");
	return -1;
    }
    else {
	char buf[1024];

	sprintf(buf, "%s/align_lib_nuc_matrix", env);
	ierr = set_alignment_matrix(buf,"ACGTURYMWSKDHVB-*");
	if (ierr) {
	    verror(ERR_FATAL, "align_contigs",
		   "%s: file not found", buf);
	    return -1;
	}
    }
    if (NULL == (params = create_align_params())) return -1;

    band = set_band_blocks(overlap->seq1_len,overlap->seq2_len);
    
    if (set_align_params (params, band, gap_open, gap_extend, edge_mode, job,
			  0, 0, 0, 0,0)) {
	destroy_alignment_params (params);
	return -1;
    };

    if ( longest_diagonal < DO_BLOCKS ) {

	ierr = affine_align(overlap,params);
	destroy_alignment_params (params);
	return ierr; /* Correct or failure */
    }

    max_matches = 100; /* dynamically grows as needed */
    word_len = 8;

    compare_method = 31;

    if ( init_hash8n ( longest_diagonal, longest_diagonal,
		     word_len, max_matches, min_match, compare_method, &h )) {
	free_hash8n(h);
	return -1;
    }

    h->seq1_len = overlap->seq1_len;
    h->seq2_len = overlap->seq2_len;
    h->seq1 = overlap->seq1;
    h->seq2 = overlap->seq2;

    if (hash_seqn(h, 1)) {
	free_hash8n(h);
      return -1;
    }
    
    if (hash_seqn(h, 2)) {
	free_hash8n(h);
      return -1;
    }

    store_hashn ( h );

    ierr = compare_b ( h, params, overlap );

    free_hash8n(h);
    if (ierr > 0) {
	/*
	 * need to check percentage.
	 * We pass the alignment either if the match is better than OK_MATCH
	 * (80% at the mo.) or if we cannot do banded alignment due to the
	 * size. So poor alignments are allowed in this case as it's not up
	 * to this code to decide what's good and what isn't.
	 */
	if(overlap->percent > OK_MATCH ||
	   longest_diagonal >= TOO_LONG_FOR_2ND_TRY) {
	    destroy_alignment_params (params);
	    return 0;
	}
    }
    /* if block alignment fails, try straight dynamic programming */


    verror(ERR_WARN, "align_contigs",
	   "Fast hashing alignment algorithm failed, "
	   "attempting full dynamic programming instead");

    if(longest_diagonal < TOO_LONG_FOR_2ND_TRY) {
	band = set_band_blocks(overlap->seq1_len,overlap->seq2_len);
	if (set_align_params (params, band, gap_open, gap_extend, edge_mode, job,
			      0, 0, 0, 0,0)) {
	    destroy_alignment_params (params);
	    return -1;
	};


        free_overlap(overlap);
	ierr = affine_align(overlap,params);

	destroy_alignment_params (params);
	return ierr;
    }
	
    verror(ERR_WARN, "align_contigs",
	   "Too large for practical use of dynamic programming");
    destroy_alignment_params (params);
    return -1;
}

/*
 * The guts of the join editor "align" function.
 * *shift returns how much we move xx1@pos1 relative to xx0@pos0. If it's
 * positive we move it right, if it's negative we move it left.
 */
static int align(edview *xx0, int pos0, int len0,
		 edview *xx1, int pos1, int len1,
		 int fixed_left, int fixed_right,
		 int *shift)
{

    char *ol0,*ol1, *cons0, *cons1;
    OVERLAP *overlap;
    char PAD_SYM = '.';
    int  *depad_to_pad0, *dp0, *depad_to_pad1, *dp1;
    int *S, *res;
    int off0 = 0, off1 = 0;
    int left0 = 0, left1 = 0;

    vfuncheader("Align contigs (join editor)");

    /* Memory allocation */
    ol0 = (char *) xmalloc(len0+1);
    ol1 = (char *) xmalloc(len1+1);
    cons0 = (char *) xmalloc(len0+1);
    cons1 = (char *) xmalloc(len1+1);
    dp0 = depad_to_pad0 = (int *)xmalloc((len0+1) * sizeof(int));
    dp1 = depad_to_pad1 = (int *)xmalloc((len1+1) * sizeof(int));

    /* Compute the consensus */
    calculate_consensus_simple(xx0->io, xx0->cnum, pos0, pos0+len0, ol0, NULL);
    calculate_consensus_simple(xx1->io, xx1->cnum, pos1, pos1+len1, ol1, NULL);

    memcpy(cons0, ol0, len0+1);
    memcpy(cons1, ol1, len1+1);

    /* Strip the pads from the consensus */
    depad_seq(ol0, &len0, depad_to_pad0);
    depad_seq(ol1, &len1, depad_to_pad1);

    if (NULL == (overlap = create_overlap())) return -1;
    init_overlap (overlap, ol0, ol1, len0, len1);

    if(-1 == align_contigs (overlap, fixed_left, fixed_right)) {
	xfree(ol0);
	xfree(ol1);
	destroy_overlap(overlap);
	return -1;
    }

    /*
    overlap->seq1_out[overlap->right+1] = 0;
    overlap->seq2_out[overlap->right+1] = 0;
    */

    S = res = rsalign2myers(overlap->seq1_out, strlen(overlap->seq1_out),
			    overlap->seq2_out, strlen(overlap->seq2_out),
			    PAD_SYM);

    /* Clip left end */
    if (*S != 0) {
	/* Pad at start, so shift contigs */
	if (*S < 0) {
	    left0 = -*S; /* used for display only */
	    depad_to_pad0 += -*S;
	    off0 = depad_to_pad0[0];
	    pos0 += off0;
	    *shift = -off0;
	    len0 -= -*S;
	} else {
	    left1 = *S; /* used for display only */
	    depad_to_pad1 += *S;
	    off1 = depad_to_pad1[0];
	    pos1 += off1;
	    *shift = off1;
	    len1 -= *S;
	}
	S++;
    } else {
	*shift = 0;
    }

    /* Clip right end */
    {
	/* FIXME: should this be shadowed or not? */
	int pos0 = 0, pos1 = 0;
	int *s = S;

	while (pos0 < len0 && pos1 < len1) {
	    if (*s < 0) {
		pos0 -= *s;
	    } else if (*s > 0) {
		pos1 += *s;
	    } else {
		pos0++;
		pos1++;
	    }

	    s++;
	}

	if (*s < 0)
	    len0 += *s;
	else if (*s > 0)
	    len1 -= *s;
    }

    /* Display the alignment. */
    {
	char *exp0, *exp1;
	int exp_len0, exp_len1;
	char name0[100];
	char name1[100];

	exp0 = (char *) xmalloc(len0+len1+1);
	exp1 = (char *) xmalloc(len0+len1+1);

	sprintf(name0, "%"PRIrec, xx0->cnum);
	sprintf(name1, "%"PRIrec, xx1->cnum);
	cexpand(ol0+left0, ol1+left1, len0, len1,
		exp0, exp1, &exp_len0, &exp_len1, 
		ALIGN_J_SSH | ALIGN_J_PADS, S);
	list_alignment(exp0, exp1, name0, name1, pos0, pos1, "");

	xfree(exp0);
	xfree(exp1);
    }

    /*************************************************************************/
    /* Now actually make the edits, keeping track of old and new pads. */
    {
	int depad_pos0 = 0, depad_pos1 = 0;
	int curr_pad0;  /* Current padded position in seq 0 */
	int curr_pad1;  /* Current padded position in seq 1 */
	int extra_pads; /* Difference between padded positions */
	int last_pad0 = -1;
	int last_pad1 = -1;
	int inserted_bases0 = 0;
	int inserted_bases1 = 0;
	contig_t *ctg0, *ctg1;

	ctg0 = cache_search(xx0->io, GT_Contig, xx0->cnum);
	cache_incr(xx0->io, ctg0);

	ctg1 = cache_search(xx1->io, GT_Contig, xx1->cnum);
	cache_incr(xx1->io, ctg1);

	while (depad_pos0 < len0 && depad_pos1 < len1) {
	    if (*S < 0) {
		depad_pos0 -= *S;
	    } else if (*S > 0) {
		depad_pos1 += *S;
	    }

	    if (depad_pos0 >= len0 || depad_pos1 >= len1)
		break;

	    curr_pad0 = depad_to_pad0[depad_pos0]-off0;
	    curr_pad1 = depad_to_pad1[depad_pos1]-off1;

	    extra_pads = (curr_pad1 - last_pad1) - (curr_pad0 - last_pad0);

	    if (extra_pads < 0) { /* Add to seq 0 */
		add_pads(xx1, &ctg1,
			 pos1 + curr_pad1 + inserted_bases1, -extra_pads);
		inserted_bases1 -= extra_pads;
	    } else if (extra_pads > 0) { /* Add to seq 1 */
		add_pads(xx0, &ctg0,
			 pos0 + curr_pad0 + inserted_bases0,  extra_pads);
		inserted_bases0 += extra_pads;
	    }
	    
	    last_pad0 = curr_pad0;
	    last_pad1 = curr_pad1;

	    if (*S == 0) {
		depad_pos0++;
		depad_pos1++;
	    }

	    S++;
	}

	cache_decr(xx0->io, ctg0);
	cache_decr(xx1->io, ctg1);
    }
    /*************************************************************************/

    xfree(res);
    xfree(ol0);
    xfree(ol1);
    xfree(dp0);
    xfree(dp1);
    destroy_overlap(overlap);

    return 0;
}

/*
 * Handles the align button in the join editor
 * Returns 0 for success
 *        -1 for failure
 */
int edJoinAlign(edview *xx, int fixed_left, int fixed_right) {
    int left0,right0;
    int left1,right1;
    int offset, ret;
    int overlapLength;
    int len0,len1;
    int shift = 0, extra;
    edview **xx2;

    int l0, l1, r0, r1; /* contig used extents */

    if (!xx->link)
	return -1;
    xx2 = xx->link->xx;
    offset = xx2[1]->displayPos - xx2[0]->displayPos;

    /* Compute overlap position and sizes */
    consensus_valid_range(xx2[0]->io, xx2[0]->cnum, &l0, &r0);
    consensus_valid_range(xx2[1]->io, xx2[1]->cnum, &l1, &r1);

    /*
     * dash => actual sequence
     * dots => contig range (start to end)
     *
     *                  l1\            /r1
     * 1:            ......------------....
     *                     ||||||||||||
     * 0:  ......-------------------------------......
     *        l0/                               \r0
     *           <--------->
     *            (-offset)
     *
     *  If l1 left of l0 then offset is +ve.
     */

    /* Set left0/left1 */
    if (fixed_left) {
	left0 = xx2[0]->cursor_apos;
	left1 = xx2[1]->cursor_apos;
    } else {
	if (offset < 0) {
	    left0 = l1-offset;
	    left1 = l1;
	} else {
	    left0 = l0;
	    left1 = l0+offset;
	}
    }

    /* Set right0/right1 */
    if (fixed_right) {
	right0 = xx2[0]->cursor_apos;
	right1 = xx2[1]->cursor_apos;
    } else {
	if (offset + r0 > r1) {
	    /* as in example above, r0 right of r1 */
	    right0 = r1-offset;
	    right1 = r1;
	} else {
	    right0 = r0;
	    right1 = r0+offset;
	}
    }

    overlapLength = right0 - left0+1;
    if (overlapLength <= 0) return 0; /* nothing to do */

    /* Add on extra data either end to allow for padding */
    extra = set_band_blocks(overlapLength, overlapLength)/2;

    if (!fixed_left) {
	left0 -= extra;
	left1 -= extra;

	if (left0 < l0) left0 = l0;
	if (left1 < l1) left1 = l1;
    }

    if (!fixed_right) {
	right0 += extra;
	right1 += extra;

	if (right0 > r0) right0 = r0;
	if (right1 > r1) right1 = r1;
    }

    len0 = right0 - left0+1;
    len1 = right1 - left1+1;

    if (len0 <= 0 || len1 <= 0)
	return 0;

    /* Do the actual alignment */
    ret = align(xx2[0], left0, len0, xx2[1], left1, len1,
		fixed_left, fixed_right, &shift);

    /* Help force full redraw */
    if (xx->r) {
	free(xx->r);
	xx->r = NULL;
    }

    if (ret)
	return ret;

    xx2[1]->displayPos = left1+shift - left0 + 1 + (xx2[0]->displayPos-1);
//  xx2[0]->displayPos = 1                       + (xx2[0]->displayPos-1);

    xx->link->lockOffset = xx2[1]->displayPos - xx2[0]->displayPos;

    xx2[0]->refresh_flags = ED_DISP_ALL;
    edview_redraw(xx2[0]);

    xx2[1]->refresh_flags = ED_DISP_ALL;
    edview_redraw(xx2[1]);

    return ret;
}

/*
 * Returns the length and number of consensus mismatches for an overlap.
 * This is copied from the start of the edJoinAlign code.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int edJoinMismatch(edview *xx, int *len, int *mismatch) {
    int left0,right0;
    int left1,right1;
    int offset;
    int overlapLength;
    int len0,len1;
    edview **xx2;
    char *c0, *c1;
    int i;

    int l0, l1, r0, r1; /* contig used extents */

    *len = 0;
    *mismatch = 0;

    if (!xx->link)
	return -1;
    xx2 = xx->link->xx;
    offset = xx2[1]->displayPos - xx2[0]->displayPos;

    /*
     * dash => actual sequence
     * dots => contig range (start to end)
     *
     *                  l1\            /r1
     * 1:            ......------------....
     *                     ||||||||||||
     * 0:  ......-------------------------------......
     *        l0/                               \r0
     *           <--------->
     *            (-offset)
     *
     *  If l1 left of l0 then offset is +ve.
     */

    /* Compute overlap position and sizes */
    consensus_valid_range(xx2[0]->io, xx2[0]->cnum, &l0, &r0);
    consensus_valid_range(xx2[1]->io, xx2[1]->cnum, &l1, &r1);

    /* Set left0/left1 */
    if (offset+l0 < l1) {
	left0 = l1-offset;
	left1 = l1;
    } else {
	left0 = l0;
	left1 = l0+offset;
    }

    if (offset+r0 > r1) {
	/* as in example above, r0 right of r1 */
	right0 = r1-offset;
	right1 = r1;
    } else {
	right0 = r0;
	right1 = r0+offset;
    }

    overlapLength = right0 - left0+1;

    /* No overlap => error */
    if (overlapLength <= 0)
	return -1;

    if (left0 < l0) left0 = l0;
    if (left1 < l1) left1 = l1;

    if (right0 > r0) right0 = r0;
    if (right1 > r1) right1 = r1;

    len0 = right0 - left0+1;
    len1 = right1 - left1+1;

    if (len0 <= 0 || len1 <= 0)
	return -1;

    assert(len0 == len1);

    c0 = xmalloc(len0+1);
    c1 = xmalloc(len1+1);
    calculate_consensus_simple(xx2[0]->io, xx2[0]->cnum, left0, right0,
			       c0, NULL);
    calculate_consensus_simple(xx2[1]->io, xx2[1]->cnum, left1, right1,
			       c1, NULL);

    *mismatch = 0;
    for (i = 0; i < len0; i++) {
	if (c0[i] != c1[i])
	    (*mismatch)++;
    }
    *len = len0;

    free(c0);
    free(c1);

    return 0;
}

#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

/*
 * Given lbin and rbin as the root bins of the left and right contig,
 * this looks for the optimal bin number to hang rbin off. If the bin
 * we return has both left and right children, then we need to create
 * a new parent containing rbin and the bin number returned, otherwise
 * use the spare child.
 *
 * Returns bin record number on success
 *         -1 on failure
 */
tg_rec find_join_bin(GapIO *io, tg_rec lbin, tg_rec rbin, int offset,
		     int offsetr, int junction) {
    bin_index_t *binl, *binr;
    int complement = 0;
    int i, f_a, f_b;
    int start, end;
    tg_rec bnum;

    binr = (bin_index_t *)cache_search(io, GT_Bin, rbin);
    binl = (bin_index_t *)cache_search(io, GT_Bin, lbin);

    start = junction + binr->pos;
    end   = junction + binr->pos + binr->size;

    /* Check who is likely to contain who */
    if (binr->size > binl->size) {
	/* Swap */
	lbin = binr->rec;
	binr = binl;
	rbin = binr->rec;
	offset = offsetr;
    }

    do {
	int child = -1;
	int off_new = offset;
	bnum = lbin;

	binl = (bin_index_t *)cache_search(io, GT_Bin, lbin);

	if (binl->flags & BIN_COMPLEMENTED) {
	    complement ^= 1;
	}

	if (complement) {
	    f_a = -1;
	    f_b = offset + binl->size-1;
	} else {
	    f_a = +1;
	    f_b = offset;
	}

	for (i =0; i < 2; i++) {
	    bin_index_t *ch;
	    if (!binl->child[i])
		continue;

	    ch = get_bin(io, binl->child[i]);

	     gio_debug(io, 1,
		       "Checking bin %"PRIrec" abs pos %d..%d vs %d..%d\n",
		       ch->rec,
		       NMIN(ch->pos, ch->pos+ch->size-1),
		       NMAX(ch->pos, ch->pos+ch->size-1),
		       start, end);

	    if (NMIN(ch->pos, ch->pos+ch->size-1) <= start &&
		NMAX(ch->pos, ch->pos+ch->size-1) >= end) {
		child = i;
		off_new = NMIN(ch->pos, ch->pos+ch->size-1);
	    }
	}

	if (child >= 0) {
	    lbin = binl->child[child];
	    offset = off_new;
	} else {
	    lbin = 0;
	}
    } while (lbin);
    
     gio_debug(io, 1, "Optimal bin to insert is above %"PRIrec"\n", bnum);

    return bnum;
}

/*
 * Invalidates cached tracks and consensus in the region covered by the
 * contig overlap.
 * 
 * Return 0 on success
 *       -1 on failure
 */
static int join_invalidate(GapIO *io, contig_t *leftc, contig_t *rightc,
			   int junction) {
    rangec_t *r;
    int i, j, nr, start, end;
    contig_t *c;

    /* Invalidate left contig */
    start = junction;
    end = leftc->end;
    c = leftc;
    for (j = 0; j < 2; j++) {
	r = contig_bins_in_range(io, &c, start, end,
				 CSIR_LEAVES_ONLY, CONS_BIN_SIZE, &nr);
	if (NULL == r) return -1;

	for (i = 0; i < nr; i++) {
	    bin_index_t *bin = cache_search(io, GT_Bin, r[i].rec);
	    if (NULL == bin) return -1;

	    if (bin->flags & BIN_CONS_VALID) {
		if (NULL == (bin = cache_rw(io, bin))) return -1;
		bin->flags |= BIN_BIN_UPDATED;
		bin->flags &= ~BIN_CONS_VALID;
	    }

	    gio_debug(io, 1, "Invalidating consensus in ctg %s, bin %"PRIrec
		      ": %d..%d (%d)\n",
		      j ? "right" : "left",
		      r[i].rec, r[i].start, r[i].end, r[i].end - r[i].start);
	}
	free(r);

	/* Invalidate right contig */
	start = rightc->start;
	end = leftc->end - junction;
	c = rightc;
    }

    return 0;
}


/*
 * Returns the last refpos marker position (passed) for a contig.
 * Returns contig_start if none present.
 */
int last_refpos (GapIO *io, contig_t *c) {
    contig_iterator *ci;
    rangec_t *r;
    int pos;

    ci = contig_iter_new_by_type(io, c->rec, 0, CITER_LAST,
				 CITER_CSTART, CITER_CEND,
				 GRANGE_FLAG_ISREFPOS);
    if (!ci)
	return contig_get_start(&c);

    r = contig_iter_next(io, ci);
    if (!r) {
	contig_iter_del(ci);
	return contig_get_start(&c);
    }

    pos = r->start;

    contig_iter_del(ci);
    return pos;
}


/*
 * Removes REFPOS markers in contig 'c' between 'from' and 'to'. We do this
 * for the overlapping contig region when joining to avoid regions of our
 * contig from being indicated as having multiple references.
 * While technically that is true after a join, it adds excessive complexity
 * and grants abilities that few people want.
 * 
 * Return 0 on success
 *       -1 on failure
 */
static int contig_remove_refpos_markers(GapIO *io, contig_t *c,
					int from, int to) {
    contig_iterator *ci;
    rangec_t *rc;

    ci = contig_iter_new_by_type(io, c->rec, 0, CITER_FIRST, from, to,
				 GRANGE_FLAG_ISREFPOS);
    if (!ci)
	return 0;

    while ((rc = contig_iter_next(io, ci))) {
	bin_index_t *bin = cache_search(io, GT_Bin, rc->orig_rec);
	range_t *r;

	if (NULL == bin) return -1;
	if (NULL == (bin = cache_rw(io, bin))) return -1;
	r = arrp(range_t, bin->rng, rc->orig_ind);
	assert((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS);

	r->flags |= GRANGE_FLAG_UNUSED;
	r->rec = bin->rng_free;
	bin->rng_free = rc->orig_ind;

	bin->flags |= BIN_BIN_UPDATED | BIN_RANGE_UPDATED;
	bin_incr_nrefpos(io, bin, -1);

	/* Otherwise adjust start/end used if we may have invalidated it */
	/* Cache till later to avoid O(N^2) complexity in worst case? */
	if (bin->start_used == r->start || bin->end_used == r->end)
	    bin_set_used_range(io, bin);
    }

    contig_iter_del(ci);
    return 0;
}

/*
 * Insert one or more new bins between bin and parent.  The location of the
 * new bin in parent is given by pos and size.  If new_bin is still
 * more than four times the size of bin, the function is called recursively
 * to add another intermediate that is approximately half the size.  This
 * should then create a resonably optimal tree structure.
 *
 * Returns a pointer to bin's new parent on success
 *         NULL on failure
 */

static bin_index_t *add_intermediate_bin(GapIO *io, bin_index_t *bin,
					 bin_index_t *parent,
					 int pos, int size) {
    tg_rec new_id;
    bin_index_t *new_bin;
    int i;

    if (NULL == (bin = cache_rw(io, bin))) return NULL;
    if (NULL == (parent = cache_rw(io, parent))) return NULL;

    new_id = bin_new(io, pos, size, bin->parent, bin->parent_type);
    if (new_id < 0) return NULL;

    gio_debug(io, 1,
	    "Adding new bin %"PRIrec" between %"PRIrec" and %"PRIrec"\n",
	    new_id, parent->rec, bin->rec);

    new_bin = (bin_index_t *) cache_search(io, GT_Bin, new_id);
    if (NULL == new_bin) return NULL;
    if (NULL == (new_bin = cache_rw(io, new_bin))) return NULL;

    new_bin->nseqs    = bin->nseqs;
    new_bin->nrefpos  = bin->nrefpos;
    new_bin->nanno    = bin->nanno;
    new_bin->child[0] = bin->rec;
    new_bin->flags   |= BIN_BIN_UPDATED;

    bin->parent      = new_id;
    bin->parent_type = GT_Bin;
    bin->pos        -= pos;
    bin->flags      |= BIN_BIN_UPDATED;

    for (i = 0; i < 2; i++) {
	if (parent->child[i] == bin->rec) parent->child[i] = new_bin->rec;
    }
    parent->flags |= BIN_BIN_UPDATED;

    if (bin->size * 4 < size) {
	/* Add more intermediates to try to keep the tree in shape */
	int sz_by_2 = size / 2;
	if (bin->pos + bin->size - sz_by_2 < sz_by_2 - bin->pos) {
	    new_bin = add_intermediate_bin(io, bin, new_bin, 0,
					   MAX(sz_by_2, bin->pos + bin->size));
	    if (NULL == new_bin) return NULL;
	} else {
	    int st = MIN(bin->pos, sz_by_2);
	    new_bin = add_intermediate_bin(io, bin, new_bin, st, size - st);
	    if (NULL == new_bin) return NULL;
	}
    }
    
    return new_bin;
}

/*
 * Ensure that bin and all of its children are big enough to reach from the
 * edge of parent to the start/end of sibling (if present).  This is necessary
 * to ensure that there are no regions in a contig where the bin tree can't
 * grow as there are no available child slots.
 *
 * If bin needs to grow by more than 50%, this is done by inserting a new bin
 * above it, thus adding a new empty child slot.  For smaller expansions,
 * bin itself is grown so that it covers the entire region.  As doing this
 * may cause the child bins to become too small, they are recursively
 * expanded as well.
 *
 * Return 0 on success
 *       -1 on failure
 */

static int recursive_grow_bins(GapIO *io, bin_index_t *bin,
			       bin_index_t *parent, bin_index_t *sibling) {
    int i, nkids;
    int shift;
    int shifted = 0;
    int new_size;
    bin_index_t *cbin[2] = { NULL, NULL };
    int comp = bin->flags & BIN_COMPLEMENTED;
    int free_start = 0;
    int free_end   = parent->size;
    int ret = -1;

    if (sibling) {
	if (sibling->pos < bin->pos) {
	    free_start = sibling->pos + sibling->size;
	} else {
	    free_end = sibling->pos;
	}
    } else {
	if (bin->pos < parent->size - (bin->pos + bin->size)) {
	    free_end = bin->pos + bin->size;
	} else {
	    free_start = bin->pos;
	}
    }
    new_size = free_end - free_start;
    shift = comp ? free_end - bin->pos - bin->size : bin->pos - free_start;

    gio_debug(io, 1, "Growing bins for %"PRIrec" %d..%d to %d..%d "
	    "parent %"PRIrec" 0..%d\n",
	    bin->rec, bin->pos, bin->pos + bin->size,
	    free_start, free_end, parent->rec, parent->size);

    if (shift == 0 && new_size == bin->size) return 0;

    if (NULL == (bin = cache_rw(io, bin))) return -1;

    if (bin->size * 3 / 2 < new_size) {
	bin_index_t *new_bin = add_intermediate_bin(io, bin, parent,
						    free_start, new_size);
	if (NULL == new_bin) return -1;

	/* Call again in case we needed to expand both ends */
	return recursive_grow_bins(io, bin, new_bin, NULL);
    }

    /* Adjust positions in the bin's range array */
    if (bin->rng && shift != 0) {
	int n = ArrayMax(bin->rng);
	for (i = 0; i < n; i++) {
	    range_t *r = arrp(range_t, bin->rng, i);
	    if (r->flags & GRANGE_FLAG_UNUSED) continue;
	    r->start += shift;
	    r->end   += shift;
	    shifted++;
	}
	if (shifted) bin->flags |= BIN_RANGE_UPDATED;
    }

    /* Collect children and update positions */
    for (i = 0, nkids = 0; i < 2; i++) {
	if (!bin->child[i]) continue;
	cbin[nkids] = get_bin(io, bin->child[i]);
	if (NULL == cbin[nkids]) goto clean;
	cache_incr(io, cbin[nkids]);
	if (shift != 0) {
	    cbin[nkids] = cache_rw(io, cbin[nkids]);
	    cbin[nkids]->pos += shift;
	    cbin[nkids]->flags |= BIN_BIN_UPDATED;
	}
	nkids++;
    }
    bin->pos  = free_start;
    bin->size = new_size;
    if (shifted) {
	bin->start_used += shift;
	bin->end_used   += shift;
    }
    bin->flags |= BIN_BIN_UPDATED;

    /* Expand children as well */
    for (i = 0; i < nkids; i++) {
	if (0 != recursive_grow_bins(io, cbin[i], bin, cbin[1 - i])) goto clean;
	cache_decr(io, cbin[i]);
	cbin[i] = NULL;
    }
    nkids = 0;
    ret = 0;
 clean:
    for (i = 0; i < nkids; i++) {
	if (NULL != cbin[i]) cache_decr(io, cbin[i]);
    }
    return ret;
}

/*
 * Transplant binr so that it becomes a new child of binl.  If binr does not
 * occupy all of the space that it should under binl, it is expanded using
 * recursive_grow_bins.
 * 
 * Returns 0 on success
 *        -1 on failure
 */

static int do_transplant(GapIO *io, bin_index_t *binl, bin_index_t *binr,
			 int comp, int pos, int child_idx,
			 bin_index_t *sibling) {
    bin_index_t *old_parent;

    gio_debug(io, 1, "Transplanting %"PRIrec" to %"PRIrec" index %d\n",
	    binr->rec, binl->rec, child_idx);
    if (NULL == (binl = cache_rw(io, binl))) return -1;
    if (NULL == (binr = cache_rw(io, binr))) return -1;

    old_parent = get_bin(io, binr->parent);
    if (NULL == old_parent) return -1;
    if (NULL == (old_parent = cache_rw(io, old_parent))) return -1;

    binl->child[child_idx] = binr->rec;
    binl->flags |= BIN_BIN_UPDATED;
    binr->pos    = pos;
    binr->parent = binl->rec;
    binr->parent_type = GT_Bin;
    binr->flags |= BIN_BIN_UPDATED;
    if (comp) binr->flags ^= BIN_COMPLEMENTED;
    if (old_parent->child[0] == binr->rec) old_parent->child[0] = 0;
    if (old_parent->child[1] == binr->rec) old_parent->child[1] = 0;
    old_parent->flags |= BIN_BIN_UPDATED;

    if (0 != bin_incr_nseq(io, binl,    binr->nseqs)) return -1;
    if (0 != bin_incr_nrefpos(io, binl, binr->nrefpos)) return -1;
    if (0 != bin_incr_nanno(io, binl,   binr->nanno)) return -1;
    if (0 != bin_incr_nseq(io, old_parent,    -binr->nseqs)) return -1;
    if (0 != bin_incr_nrefpos(io, old_parent, -binr->nrefpos)) return -1;
    if (0 != bin_incr_nanno(io, old_parent,   -binr->nanno)) return -1;
    if (0 != recursive_grow_bins(io, binr, binl, sibling)) return -1;
    return 0;
}

/*
 * Attempt to transplant binr + all of its children to underneath binl.
 * comp indicates if the bin needs to be complemented.  pos is the
 * position in the destination bin.
 *
 * Returns 0 if binr was successfully transplanted
 *         1 if there was no room for binr so nothing happened
 *        -1 if an error occurred
 */

static int transplant_bin(GapIO *io, bin_index_t *binl, bin_index_t *binr,
			  int comp, int pos) {
    /* See if binr will fit underneath binl.  If it does we can
       transplant it and all of its children directly. */
    int end = pos + binr->size;
    int free, used, ret;
    bin_index_t *ch;

    if (binr->parent_type != GT_Bin) return 1; /* Don't attempt to move root */
    if (binl->child[0] && binl->child[1]) { return 1; }
    if (!binl->child[0] && !binl->child[1]) {
	/* Both child slots free, work out best one to choose */
	int rfree = binl->size - pos - binr->size;
	return do_transplant(io, binl, binr, comp, pos,
			     pos < rfree ? 0 : 1, NULL);
    }
    free = binl->child[0] ? 1 : 0;
    used = 1 - free;
    ch = get_bin(io, binl->child[used]);
    if (NULL == ch) return -1;

    gio_debug(io, 1, "-- child #%d %"PRIrec" %d..%d\n",
	    used, ch->rec, ch->pos, ch->pos + ch->size);
    if ((end < ch->pos) || (ch->pos + ch->size <= pos)) { 
	/* Transplant to child[0] */
	cache_incr(io, ch);
	ret = do_transplant(io, binl, binr, comp, pos, free, ch);
	cache_decr(io, ch);
	return ret;
    }

    /* If we got here, there isn't enough room so give up */
    return 1;
}

static int round_up(int val) {
    int bits = 0;
    while (val) {
	val /= 2;
	bits++;
    }

    return 1<<bits;
}

/*
 * Extend the root bin of contig c so that it covers the range start..end
 * We do this by creating zero, one or two extra bins above the existing root
 * depending on if we need to extend one or both ends.  If the bin is already
 * big enough, this will do nothing.
 *
 * Returns 0 on success
 *        -1 on failure
 */

static int extend_root_bin(GapIO *io, contig_t *c, int start, int end) {
    bin_index_t *old_root, *new_root;
    tg_rec       new_id;
    int          bin_start, bin_end;

    old_root = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(&c));
    if (NULL == old_root) return -1;
    bin_start = old_root->pos;
    bin_end   = old_root->pos + old_root->size;

    if (bin_start <= start && end <= bin_end) return 0; /* Nowt to do */

    if (start < bin_start && end > bin_end) {
	/* Need to do both ends, so need two new bins */
	int ret = extend_root_bin(io, c, bin_start, end);
	if (ret) return ret;
	old_root = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(&c));
	if (NULL == old_root) return -1;
	bin_start = old_root->pos;
	bin_end   = old_root->pos + old_root->size;
    }

    /* Round the new bin up to a power of 2 in length.  With luck this will
       help stop the bin structure from becoming lop-sided through long 
       sequences of joins. */
    if (start < bin_start) start = bin_end   - round_up(bin_end - start);
    if (end   > bin_end)   end   = bin_start + round_up(end - bin_start);

    if (NULL == (old_root = cache_rw(io, old_root))) return -1;
    
    if ((new_id = bin_new(io, 0, 0, c->rec, GT_Contig)) < 0) return -1;
    new_root = (bin_index_t *)cache_search(io, GT_Bin, new_id);
    if (NULL == new_root) return -1;
    if (NULL == (new_root = cache_rw(io, new_root))) return -1;

    if (0 != contig_set_bin(io, &c, new_id)) return -1;
    gio_debug(io, 1, "Made new root bin %"PRIrec" for contig %"PRIrec"\n",
	    new_id, c->rec);

    new_root->nseqs    = old_root->nseqs;
    new_root->nrefpos  = old_root->nrefpos;
    new_root->nanno    = old_root->nanno;
    new_root->child[0] = old_root->rec;
    new_root->pos      = MIN(start, bin_start);
    new_root->size     = MAX(end, bin_end) - new_root->pos;
    new_root->flags   |= BIN_BIN_UPDATED;

    old_root->parent      = new_root->rec;
    old_root->parent_type = GT_Bin;
    old_root->pos        -= new_root->pos;
    old_root->flags      |= BIN_BIN_UPDATED;

    return 0;
}

/*
 * Attempt to join contigs by transplanting bins from one contig to the other.
 * 
 */

static int join_move_bins(GapIO *io, contig_t *cl, contig_t *cr,
			   int offset) {
    typedef struct bin_list {
        tg_rec rec;
	int parent_comp;
	int parent_offset;
	int parent_size;
        struct bin_list *next;
    } bin_list;
    bin_index_t *binl = NULL, *binr = NULL;
    bin_list *head = NULL;
    bin_list *tail = NULL;
    bin_list *item;
    bin_index_t *cbin[2] = { NULL, NULL };
    int nkids;
    int i;
    int ret = -1;

    binl = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(&cl));
    if (NULL == binl) return -1;
    cache_incr(io, binl);
    binr = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(&cr));
    if (NULL == binr) goto clean;
    cache_incr(io, binr);

    /* Make cl's root bin big enough to cover the joined contig, so
       bin_for_range will always return something.  */

    if (0 != extend_root_bin(io, cl,
			     MIN(binl->pos, binr->pos + offset),
			     MAX(binl->pos + binl->size,
				 binr->pos + binr->size + offset))) goto clean;
    cache_decr(io, binl); binl = NULL;
    cache_decr(io, binr); binr = NULL;

    gio_debug(io, 1, "cl %d..%d cr %d..%d offset %d\n",
	    contig_get_start(&cl), contig_get_end(&cl),
	    contig_get_start(&cr), contig_get_end(&cr), offset);
    if (0 != contig_set_start(io, &cl, MIN(contig_get_start(&cl),
					   contig_get_start(&cr) + offset))) {
	goto clean;
    }
    if (0 != contig_set_end(io, &cl, MAX(contig_get_end(&cl),
					 contig_get_end(&cr) + offset))) {
	goto clean;
    }

    /* Initialize bin_list */
    head = malloc(sizeof(bin_list));
    if (!head) goto clean;
    head->rec = contig_get_bin(&cr);
    head->parent_comp = 0;
    head->parent_offset = 0;
    head->parent_size = 0;
    head->next = NULL;
    tail = head;

    /* Breadth-first recursive bin mangling */
    while (head) {
	int rcomp = head->parent_comp;
	int rstart, rend, lstart = 0, lcomp = 0, new_pos;
	
	binr = get_bin(io, head->rec);
	if (NULL == binr) goto clean;
	cache_incr(io, binr);

	if (binr->flags & BIN_COMPLEMENTED) {
	    rcomp ^= 1;
	}
	rstart = (head->parent_comp
		  ? head->parent_offset+head->parent_size - binr->pos-binr->size
		  : head->parent_offset + binr->pos);
	rend = rstart + binr->size;
	binl = bin_for_range(io, &cl, rstart + offset, rend + offset - 1, 0,
			     &lstart, &lcomp);
	if (binl) {
	    cache_incr(io, binl);
	    new_pos = (lcomp
		       ? lstart + binl->size - rstart - offset - binr->size
		       : rstart + offset - lstart);

	    gio_debug(io, 1, "Trying %"PRIrec" %d..%d target %"PRIrec" %d..%d pos %d end %d %s\n",
		    binr->rec, rstart + offset, rend + offset, binl->rec,
		    lstart, lstart + binl->size,
		    new_pos, new_pos + binr->size,
		    lcomp ^ rcomp ? "comp" : "uncomp");

	    if (transplant_bin(io, binl, binr,
			       lcomp ^ head->parent_comp, new_pos)) {
		/* Didn't work; try children */
		nkids = 0;
		for (i = 0; i < 2; i++) {
		    if (!binr->child[i]) continue;
		    cbin[nkids] = get_bin(io, binr->child[i]);
		    if (NULL == cbin[nkids]) goto clean;
		    cache_incr(io, cbin[nkids]);
		    nkids++;
		}
		if (nkids == 2 && cbin[0]->pos > cbin[1]->pos) {
		    bin_index_t *tmp = cbin[0];
		    cbin[0] = cbin[1];
		    cbin[1] = tmp;
		}
		for (i = 0; i < nkids; i++) {
		    item = malloc(sizeof(bin_list));
		    if (NULL == item) goto clean;
		    item->rec = cbin[i]->rec;
		    item->parent_comp = rcomp;
		    item->parent_offset = rstart;
		    item->parent_size = binr->size;
		    item->next = NULL;
		    tail->next = item;
		    tail = item;
		    cache_decr(io, cbin[i]); cbin[i] = NULL;
		}
	    }
	    cache_decr(io, binl); binl = NULL;
	}
	
	cache_decr(io, binr); binr = NULL;
	item = head;
	head = head->next;
	free(item);
    }

    ret = 0;
 clean:
    if (NULL != binl) cache_decr(io, binl);
    if (NULL != binr) cache_decr(io, binr);
    for (i = 0; i < 2; i++) {
	if (NULL != cbin[i]) cache_decr(io, cbin[i]);
    }
    while (NULL != head) {
	item = head;
	head = head->next;
	free(item);
    }
    return 0;
}

/*
 * Move all sequences and reference position markers from bin to the
 * appropriate place in destination contig c, which should not be the one
 * that bin is in.  Any cached consensus records will be freed.
 *
 * The record numbers of the new bin for each sequence is put into seq_bins
 * for use by bin_move_annos later on.  If any annotations are seen while
 * iterating through the bin contents then *anno_seen_out is set.  This
 * allows the caller to avoid the bin_move_annos step if it is not needed.
 * 
 * The last destination bin is returned in *new_bin_out.  This is mainly so
 * that the caller knows to call bin_add_range(io, ..., -1) at the
 * end of the porcess if it needs to.
 *
 * The number of sequences moved is returned in *seqs_moved_out, and the
 * number of reference positions on *refp_moved_out.  If the range has been
 * updated, *range_changed_out is set to 1.
 *
 * Returns 0 on sucess
 *        -1 on failure
 */
static int bin_move_seqs(GapIO *io, HacheTable *seq_bins, bin_index_t *bin,
			 contig_t *c, bin_index_t **new_bin_out,
			 int f_a, int f_b, int offset,
			 int *anno_seen_out, int *seqs_moved_out,
			 int *refp_moved_out, int *range_changed_out) {
    int n = ArrayMax(bin->rng);
    int src_comp = f_a < 0 ? 1 : 0;
    bin_index_t *new_bin = NULL;
    int j;

    for (j = 0; j < n; j++) {
	HacheData hd;
	range_t *r = arrp(range_t, bin->rng, j);
	range_t *r_new = NULL;
	int      dest_comp = 0;
	int      start = r->start, end = r->end;
	seq_t   *seq;
    
	if (r->flags & GRANGE_FLAG_UNUSED) continue;
	switch (r->flags & GRANGE_FLAG_ISMASK) {
	case GRANGE_FLAG_ISSEQ:
	    r->start = (f_a > 0 ? start : end) * f_a + f_b + offset;
	    r->end   = (f_a > 0 ? end : start) * f_a + f_b + offset;
	    new_bin = bin_add_range(io, &c, r, &r_new, &dest_comp, 1);
	    if (NULL == new_bin) return -1;
	    hd.i = new_bin->rec;
	    if (NULL == HacheTableAdd(seq_bins, (char *) &r->rec,
				      sizeof(r->rec), hd, NULL)) return -1;
	    seq = cache_search(io, GT_Seq, r_new->rec);
	    if (NULL == seq) return -1;
	    if (NULL == (seq = cache_rw(io, seq))) return -1;
	    seq->bin = new_bin->rec;
	    seq->bin_index = r_new - ArrayBase(range_t, new_bin->rng);
	    if (src_comp ^ dest_comp) {
		seq->len = -seq->len;
		seq->flags ^= SEQ_COMPLEMENTED;
	    }
	    (*seqs_moved_out)++;
	    break;
	case GRANGE_FLAG_ISREFPOS:
	    r->start = (f_a > 0 ? start : end) * f_a + f_b + offset;
	    r->end   = (f_a > 0 ? end : start) * f_a + f_b + offset;
	    new_bin = bin_add_range(io, &c, r, NULL, NULL, 1);
	    if (NULL == new_bin) return -1;
	    (*refp_moved_out)++;
	    break;
	case GRANGE_FLAG_ISCONS:
	    if (0 != cache_item_remove(io, GT_Seq, r->rec)) return -1;
	    break;
	case GRANGE_FLAG_ISANNO:
	    *anno_seen_out = 1;
	    continue;
	default:
	    continue;
	}
	r->flags |= GRANGE_FLAG_UNUSED;
	r->rec = (tg_rec) bin->rng_free;
	bin->rng_free = j;
	*range_changed_out = 1;
    }

    if (NULL != new_bin) *new_bin_out = new_bin;
    return 0;
}

/*
 * Move annotations from bin to destination contig c, which should not be the
 * one bin is in.  seq_bins should contain the record numbers of the bins
 * that any sequences were moved to by bin_move_seqs so that we can put
 * sequence annotations into the correct places.
 * The last destination bin is returned in *new_bin_out mainly so that the
 * caller knows to call bin_add_range(io, ..., -1) to fix up the bin counts
 * if necessary.  It also returns the number of annotations moved 
 * in *annos_moved_out.  If the range has been updated, *range_changed_out is
 * set to 1.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int bin_move_annos(GapIO *io, HacheTable *seq_bins, bin_index_t *bin,
			  contig_t *c, bin_index_t **new_bin_out,
			  int f_a, int f_b, int offset,
			  int *annos_moved_out, int *range_changed_out) {
    bin_index_t *new_bin = NULL;
    int n = ArrayMax(bin->rng);
    int j;
    
    for (j = 0; j < n; j++) {
	range_t    *r = arrp(range_t, bin->rng, j);
	range_t    *r_new = NULL;
	anno_ele_t *anno;
	tg_rec      target = 0;
	int         start = r->start, end = r->end;
	
	if (r->flags & GRANGE_FLAG_UNUSED) continue;
	if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO) continue;
	
	anno = cache_search(io, GT_AnnoEle, r->rec);
	if (NULL == anno) return -1;
	if (anno->obj_type == GT_Seq) {
	    HacheItem *hi = HacheTableSearch(seq_bins, (char *)&r->pair_rec,
					     sizeof(r->pair_rec));
	    if (hi) {
		target = hi->data.i;
	    } else {
		/* FIXME: Tag somehow detached from its sequence */
		continue;
	    }
	}
	
	if (NULL == (anno = cache_rw(io, anno))) return -1;

	r->start = (f_a > 0 ? start : end) * f_a + f_b + offset;
	r->end   = (f_a > 0 ? end : start) * f_a + f_b + offset;
	new_bin = bin_add_to_range(io, &c, target, r, &r_new, NULL, 1);
	if (NULL == new_bin) return -1;
	anno->bin = new_bin->rec;

	r->flags |= GRANGE_FLAG_UNUSED;
	r->rec = (tg_rec) bin->rng_free;
	bin->rng_free = j;

	*range_changed_out = 1;
	(*annos_moved_out)++;
    }

    if (NULL != new_bin) *new_bin_out = new_bin;
    return 0;
}

/*
 * Move objects (sequences, reference positions and annotations) from
 * cr to cl.  This will also remove any cached consensus, as this may
 * not be correct in the joined contig anyway.  Sequence annotations
 * should end up in the same bin as the corresponding sequence at the
 * end of the process (assuming that they were in the same one to start
 * with.)  Anything else will be left behind, but for a consistent database
 * this should be nothing.
 *
 * Returns 0 on success
 *        -1 on failure 
 */
static int join_move_objects(GapIO *io, contig_t *cl, contig_t *cr,
			     int offset) {
    HacheTable  *seq_bins = HacheTableCreate(1024,
					    HASH_DYNAMIC_SIZE|HASH_POOL_ITEMS);
    int          overlap  = contig_get_end(&cl) - offset;
    rangec_t    *bin_list  = NULL;
    bin_index_t *new_bin = NULL;
    int          nbins;
    int          i;
    int          ret = -1;

    if (NULL == seq_bins) return -1;
    
    if (overlap > contig_get_end(&cr)) overlap = contig_get_end(&cr);
    bin_list = contig_bins_in_range(io, &cr, contig_get_start(&cr),
				    overlap, 0, 0, &nbins);
    if (NULL == bin_list) goto clean;

    for (i = 0; i < nbins; i++) {
	int f_a = bin_list[i].pair_start;
	int f_b = bin_list[i].pair_end;
	int seqs_moved = 0, annos_moved = 0, refp_moved = 0, range_changed = 0;
	int anno_seen = 0;
	bin_index_t *bin;

	if (NULL == (bin = get_bin(io, bin_list[i].rec))) goto clean;

	if (!bin->rng) { continue; }
	if (NULL == (bin = cache_rw(io, bin))) goto clean;

	/* Move sequences and remove cached consensus */
	if (0 != bin_move_seqs(io, seq_bins, bin, cl, &new_bin, f_a, f_b,
			       offset, &anno_seen, &seqs_moved,
			       &refp_moved, &range_changed)) goto clean;

	/* Move annotations.  Try to get them into the same bin as the
	   corresponding sequence */
	if (anno_seen) {
	    if (0 != bin_move_annos(io, seq_bins, bin, cl, &new_bin,
				    f_a, f_b, offset,
				    &annos_moved, &range_changed)) goto clean;
	}

	gio_debug(io, 1,
		"Removed %d seqs, %d tags, %d repos from bin %"PRIrec"\n",
		seqs_moved, annos_moved, refp_moved, bin_list[i].rec);
	if (0 != (bin_incr_nseq(   io, bin, -seqs_moved)))  goto clean;
	if (0 != (bin_incr_nanno(  io, bin, -annos_moved))) goto clean;
	if (0 != (bin_incr_nrefpos(io, bin, -refp_moved)))  goto clean;
	
	if (range_changed) {
	    bin->flags |= BIN_RANGE_UPDATED|BIN_BIN_UPDATED;
	    bin->flags &= ~(BIN_CONS_CACHED|BIN_CONS_VALID);
	    if (0 != bin_set_used_range(io, bin)) goto clean;
	}
	HacheTableEmpty(seq_bins, 0);
    }

    if (NULL != new_bin) {
	bin_add_range(io, NULL, NULL, NULL, NULL, -1);
    }

    ret = 0;
 clean:
    if (NULL != bin_list) free(bin_list);
    HacheTableDestroy(seq_bins, 0);
    return ret;
}

/*
 * Depth-first recurse into the bin structure, removing empty children
 * on the way back out.
 *
 * Returns 0 if bin is empty and has no children
 *         1 if bin has children or is not empty
 *        -1 on failure
 */
static int drop_empty_bins_recurse(GapIO *io, bin_index_t *bin) {
    int c;
    int child_in_use = 0;
    int ret;

    for (c = 0; c < 2; c++) {
	if (bin->child[c]) {
	    bin_index_t *cbin = cache_search(io, GT_Bin, bin->child[c]);
	    if (NULL == cbin) return -1;

	    cache_incr(io, cbin);
	    ret = drop_empty_bins_recurse(io, cbin);
	    cache_decr(io, cbin);

	    if (ret < 0) return ret;
	    if (0 == ret) {
		cache_incr(io, cbin);
		if (NULL == (bin = cache_rw(io, bin)) ||
		    NULL == (cbin = cache_rw(io, cbin))) {
		    cache_decr(io, cbin);
		    return -1;
		}
		cache_decr(io, cbin);
		if (0 != cache_deallocate(io, cbin)) return -1;
		bin->child[c] = 0;
	    } else {
		child_in_use = 1;
	    }
	}
    }

    if (!child_in_use && bin_empty(bin)) return 0;

    return 1;
}

/* 
 * Remove empty bins from contig c.  The root bin will be left even if empty
 * as having a contig without one can cause difficulties.
 * 
 * Returns 0 if all bins in contig c were empty
 *         1 if some bins were full and so remain
 *        -1 on failure
 */
static int drop_empty_bins(GapIO *io, contig_t *c) {
    bin_index_t *bin = NULL;
    int          ret = -1;

    cache_incr(io, c);
    bin = cache_search(io, GT_Bin, c->bin);
    if (NULL == bin) goto clean;
    cache_incr(io, bin);
    ret = drop_empty_bins_recurse(io, bin);

 clean:
    if (NULL != bin) cache_decr(io, bin);
    cache_decr(io, c);
    return ret;
}

/* 
 * Join contigs by overlapping bins.  This is quick but can lead to
 * a very unbalanced bin structure.
 *
 * Returns 0 on success
 *        -1 on failure
 */

int join_overlap(GapIO *io, contig_t **cl, contig_t **cr, int offset) {
    tg_rec       binp_id;
    bin_index_t *binp, *binl, *binr;
    contig_t    *cl_rw;

    if ((binp_id = bin_new(io, 0, 0, (*cl)->rec, GT_Contig)) < 0) return -1;
    binp = (bin_index_t *)cache_search(io, GT_Bin, binp_id);
    if (NULL == binp) return -1;
    if (NULL == (binp  = cache_rw(io, binp))) return -1;

    binl = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(cl));
    if (NULL == binl) return -1;
    if (NULL == (binl  = cache_rw(io, binl))) return -1;

    binr = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(cr));
    if (NULL == binr) return -1;
    if (NULL == (binr  = cache_rw(io, binr))) return -1;

    if (NULL == (cl_rw = cache_rw(io, *cl))) return -1;

    /* Update the contig links */
    if (0 != contig_set_bin(io, cl, binp_id)) return -1;
    if (0 != contig_set_start(io, cl, MIN(contig_get_start(cl),
					  contig_get_start(cr) + offset)))
	return -1;
    if (0 != contig_set_end(io, cl, MAX(contig_get_end(cl),
					contig_get_end(cr) + offset)))
	return -1;

    /* Link the new bins together */
    binp->nseqs   = binl->nseqs   + binr->nseqs;
    binp->nrefpos = binl->nrefpos + binr->nrefpos;
    binp->nanno   = binl->nanno   + binr->nanno;
    binp->child[0] = binl->rec;
    binp->child[1] = binr->rec;
    binp->pos  = MIN(binl->pos, binr->pos + offset);
    binp->size = MAX(binl->pos + binl->size, binr->pos + binr->size + offset)
	- binp->pos + 1;
    binp->flags |= BIN_BIN_UPDATED;

    binl->parent      = binp->rec;
    binl->parent_type = GT_Bin;
    binl->pos         = binl->pos - binp->pos;
    binl->flags      |= BIN_BIN_UPDATED;

    binr->parent      = binp->rec;
    binr->parent_type = GT_Bin;
    binr->pos         = binr->pos - binp->pos + offset;
    binr->flags      |= BIN_BIN_UPDATED;

    *cl = cl_rw;

    return 0;
}

/*
 * Perform the actual join process
 * Returns 0 for success
 *        -1 for failure
 */
static int do_join_contigs(GapIO *io, tg_rec clrec, tg_rec crrec, int offset,
			   int *cr_contains_stuff) {
    contig_t *cl, *cr;
    bin_index_t *binl = NULL, *binr = NULL;
    int overlap_len;
    /* Limits before we fall back to overlapping bins */
    const int MAX_OVERLAP = 1000 * MIN_BIN_SIZE;
    const int MAX_TO_MOVE = 1000000;

    *cr_contains_stuff = 1;
    if (!(cl = cache_search(io, GT_Contig, clrec))) return -1;
    if (NULL == (cl = cache_rw(io, cl))) return -1;

    if (!(cr = cache_search(io, GT_Contig, crrec))) return -1;
    if (NULL == (cr = cache_rw(io, cr))) return -1;

    /* Invalidate any cached data in the overlapping bins */
    if (0 != join_invalidate(io, cl, cr, offset)) return -1;

    if (0 != contig_remove_refpos_markers(io, cr, contig_get_start(&cr),
					  contig_get_start(&cr) +
					  last_refpos(io, cl)-offset-1))
	return -1;

    overlap_len = MIN(contig_get_end(&cl), contig_get_end(&cr) + offset)
	- MAX(contig_get_start(&cl), contig_get_start(&cr) + offset) + 1;

    if (overlap_len > 0 && overlap_len < MAX_OVERLAP) {
	if (0 != join_move_bins(io, cl, cr, offset)) return -1;
	binr = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(&cr));
	if (NULL == binr) return -1;
	if (binr->nseqs + binr->nanno + binr->nrefpos < MAX_TO_MOVE) {
	    if (0 != join_move_objects(io, cl, cr, offset)) return -1;
	    if ((*cr_contains_stuff = drop_empty_bins(io, cr)) < 0) return -1;
	}
    }

    // assert(check_contig_bin(io, clrec) == 0);

    if (*cr_contains_stuff) {
	// assert(check_contig_bin(io, crrec) == 0);
	if (0 != join_overlap(io, &cl, &cr, offset)) return -1;
    }
    
    binl = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(&cl));
    if (NULL == binl) return -1;
    cl->nseqs   = binl->nseqs;
    cl->nanno   = binl->nanno;
    cl->nrefpos = binl->nrefpos;
    return 0;
}

int join_contigs(GapIO *io, tg_rec clrec, tg_rec crrec, int offset) {
    reg_length rl;
    reg_join rj;
    contig_t *cl, *cr;
    GapIO *child_io;
    int cr_contains_stuff = 1;
    
    /* Force joins at the top-level IO */
    while (io->base)
	io = io->base;
    if (!(child_io = gio_child(io))) return -1;
    
    if (do_join_contigs(child_io, clrec, crrec, offset, &cr_contains_stuff)) {
	gio_close(child_io); /* None of this ever happened... */
	return -1;
    }

    /* Push out all the changes so far and revert to the base io */
    if (0 != cache_flush(child_io)) return -1;
    gio_close(child_io);

    /*
     * The order of notifications here is crucial.
     *
     * We need to firstly notify the right contig that it's been joined to
     * the left contig.
     *
     * Merge the contig registration lists (copy right into left).
     *
     * Then we delete the right contig.
     *
     * Finally we then tell the left contig that it's length has changed.
     */
    
    /* Notify right of join */
    rj.job = REG_JOIN_TO;
    rj.contig = clrec;
    rj.offset = offset;
    contig_notify(io, crrec, (reg_data *)&rj);

    /* Merge lists */
    if (0 != contig_register_join(io, crrec, clrec)) return -1;

    /* Destroy the old right contig */
    if (!cr_contains_stuff) {
	/* Need to get rid of the empty root bin as well */
	if (!(cr = cache_search(io, GT_Contig, crrec))) return -1;
	if (0 != cache_rec_deallocate(io, GT_Bin, contig_get_bin(&cr)))
	    return -1;
    }
    if (0 != contig_destroy(io, crrec)) return -1;

    /* Notify left of join */
    if (!(cl = cache_search(io, GT_Contig, clrec))) return -1;
    rl.job = REG_LENGTH;
    rl.length = cl->end - cl->start + 1;
    contig_notify(io, clrec, (reg_data *)&rl);

    return cache_flush(io);
}

int edJoin(edview *xx) {
    /* p=parent contig, l=left contig, r=right contig */
    tg_rec cl, cr;
    GapIO *io = xx->io;
    int offset;

    if (!xx->link)
	return -1;

    xx->link->lockOffset = xx->link->xx[1]->displayPos -
	xx->link->xx[0]->displayPos;

    if (xx->link->lockOffset > 0) {
	cl = xx->link->xx[1]->cnum;
	cr = xx->link->xx[0]->cnum;
	offset = xx->link->lockOffset;
    } else {
	cl = xx->link->xx[0]->cnum;
	cr = xx->link->xx[1]->cnum;
	offset = -xx->link->lockOffset;
    }

    cache_flush(io);

    return join_contigs(io, cl, cr, offset);
}

