/*
 * editor_join.c:
 * Contains functions for the Join Editor "align" and "join" buttons.
 */
#include <tg_gio.h>
#include <assert.h>
#include <string.h>

#include "editor_view.h"
#include "align_lib.h"
#include "hash_lib.h"
#include "align.h"
#include "dna_utils.h"

/* Add 'num' pads into the consensus for editor 'xx' at position 'pos'. */
static void add_pads(edview *xx, int pos, int num)
{
    int i, conf = 0;

    if (num < 0)
	return;

    for (i = 0; i < num; i++)
	contig_insert_base(xx->io, &xx->contig, pos, '*', conf);
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
    job = RETURN_SEQ | RETURN_NEW_PADS;

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
		add_pads(xx1, pos1 + curr_pad1 + inserted_bases1, -extra_pads);
		inserted_bases1 -= extra_pads;
	    } else if (extra_pads > 0) { /* Add to seq 1 */
		add_pads(xx0, pos0 + curr_pad0 + inserted_bases0,  extra_pads);
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
    int shift = 0;
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

    /* Add on extra data either end to allow for padding */
#define XTRA_PERC 0.30
    overlapLength = right0 - left0+1;
    if (overlapLength <= 0) return 0; /* nothing to do */

    if (!fixed_left) {
	left0 -= (int)(overlapLength * XTRA_PERC);
	left1 -= (int)(overlapLength * XTRA_PERC);

	if (left0 < l0) left0 = l0;
	if (left1 < l1) left1 = l1;
    }

    if (!fixed_right) {
	right0  += (int)(overlapLength * XTRA_PERC);
	right1  += (int)(overlapLength * XTRA_PERC);

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

    printf("*** Alignment done\n");

    if (ret)
	return ret;

    xx2[1]->displayPos = left1+shift - left0 + 1 + (xx2[0]->displayPos-1);
    xx2[0]->displayPos = 1                       + (xx2[0]->displayPos-1);

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

	    printf("Checking bin %"PRIrec" abs pos %d..%d vs %d..%d\n",
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
    
    printf("Optimal bin to insert is above %"PRIrec"\n", bnum);

    return bnum;
}

/*
 * Invalidates cached tracks and consensus in the region covered by the
 * contig overlap.
 */
int join_invalidate(GapIO *io, contig_t *leftc, contig_t *rightc,
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
	for (i = 0; i < nr; i++) {
	    bin_index_t *bin = cache_search(io, GT_Bin, r[i].rec);

	    if (bin->flags & BIN_CONS_VALID) {
		bin = cache_rw(io, bin);
		bin->flags |= BIN_BIN_UPDATED;
		bin->flags &= ~BIN_CONS_VALID;
	    }

	    printf("Invalidating consensus in ctg %s, bin %"PRIrec": %d..%d (%d)\n",
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
 * Perform the actual join process
 * Returns 0 for success
 *        -1 for failure
 */
int edJoin(edview *xx) {
    /* p=parent contig, l=left contig, r=right contig */
    contig_t *cl, *cr;
    GapIO *io = xx->io;
    bin_index_t *binp, *binl, *binr;
    int offset;
    tg_rec above, binp_id;
    reg_length rl;
    reg_join rj;

    if (!xx->link)
	return -1;

    if (xx->link->lockOffset > 0) {
	cl = xx->link->xx[1]->contig;
	cr = xx->link->xx[0]->contig;
	offset = xx->link->lockOffset;
    } else {
	cl = xx->link->xx[0]->contig;
	cr = xx->link->xx[1]->contig;
	offset = -xx->link->lockOffset;
    }

    cache_flush(io);

    /* Force joins at the top-level IO */
    while (io->base)
	io = io->base;

    /* Invalidate any cached data in the overlapping bins */
    join_invalidate(io, cl, cr, offset);

    /* Find appropriate bin to insert our new contig above */
    above = find_join_bin(io, contig_get_bin(&cl), contig_get_bin(&cr),
			  contig_offset(io, &cl), contig_offset(io, &cr),
			  offset);

    if (above == -1)
	above = contig_get_bin(&cl);

    /* Ignore this for now */
#if 0
    /* Optimisation, hang binr off binl if it fits */
    binp_id = bin_new(io, 0, 0, cl->rec, GT_Contig);
    binp = (bin_index_t *)cache_search(io, GT_Bin, binp_id);
    binl = (bin_index_t *)cache_search(io, GT_Bin, above);
    binr = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(&cr));
    cache_incr(io, binp);
    cache_incr(io, binl);
    cache_incr(io, binr);
    cache_incr(io, cl);
    cache_incr(io, cr);
    binp = cache_rw(io, binp);
    binl = cache_rw(io, binl);
    binr = cache_rw(io, binr);
    cl   = cache_rw(io, cl);

    /* Update the contig links */
    if (above == contig_get_bin(&cl)) {
	contig_set_bin  (io, &cl, binp_id);
	contig_set_start(io, &cl, MIN(contig_get_start(&cl),
				      contig_get_start(&cr)+offset));
	contig_set_end  (io, &cl, MAX(contig_get_end(&cl),
				      contig_get_end(&cr)+offset));
    } else {
	/* Or instead an internal bin */
	bin_index_t *b2 = get_bin(io, binl->parent);
	b2 = cache_rw(io, b2);
	binp->parent = b2->rec;
	binp->parent_type = GT_Bin;
	if (b2->child[0] == binl->rec) {
	    b2->child[0] = binp->rec;
	} else if (b2->child[1] == binl->rec) {
	    b2->child[1] = binp->rec;
	} else {
	    fprintf(stderr, "invalid bin relationship");
	    abort();
	}

	/* Fix offset? */
    }

    /* Link the new bins together */
    binp->child[0] = binl->rec;
    binp->child[1] = binr->rec;
    binp->pos = MIN(binl->pos, binr->pos + offset);
    binp->size = MAX(binl->pos+binl->size, binr->pos+binr->size + offset)
	- binp->pos + 1;
    binp->flags |= BIN_BIN_UPDATED;

    binl->parent = binp->rec;
    binl->parent_type = GT_Bin;
    binl->pos = binl->pos - binp->pos;
    binl->flags |= BIN_BIN_UPDATED;

    binr->parent = binp->rec;
    binr->parent_type = GT_Bin;
    binr->pos = binr->pos - binp->pos + offset;
    binr->flags |= BIN_BIN_UPDATED;

#else
    /* Object writeable copies of our objects */
    binp_id = bin_new(io, 0, 0, cl->rec, GT_Contig);
    binp = (bin_index_t *)cache_search(io, GT_Bin, binp_id);
    binl = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(&cl));
    binr = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(&cr));
    cache_incr(io, binp);
    cache_incr(io, binl);
    cache_incr(io, binr);
    cache_incr(io, cl);
    cache_incr(io, cr);

    binp = cache_rw(io, binp);
    binl = cache_rw(io, binl);
    binr = cache_rw(io, binr);
    cl   = cache_rw(io, cl);

    /* Update the contig links */
    contig_set_bin  (io, &cl, binp_id);
    contig_set_start(io, &cl, MIN(contig_get_start(&cl),
				  contig_get_start(&cr)+offset));
    contig_set_end  (io, &cl, MAX(contig_get_end(&cl),
				  contig_get_end(&cr)+offset));

    /* Link the new bins together */
    binp->nseqs = binl->nseqs + binr->nseqs;
    binp->child[0] = binl->rec;
    binp->child[1] = binr->rec;
    binp->pos = MIN(binl->pos, binr->pos + offset);
    binp->size = MAX(binl->pos+binl->size, binr->pos+binr->size + offset)
	- binp->pos + 1;
    binp->flags |= BIN_BIN_UPDATED;

    binl->parent = binp->rec;
    binl->parent_type = GT_Bin;
    binl->pos = binl->pos - binp->pos;
    binl->flags |= BIN_BIN_UPDATED;

    binr->parent = binp->rec;
    binr->parent_type = GT_Bin;
    binr->pos = binr->pos - binp->pos + offset;
    binr->flags |= BIN_BIN_UPDATED;
#endif

    cache_decr(io, binp);
    cache_decr(io, binl);
    cache_decr(io, binr);
    cache_decr(io, cl);
    cache_decr(io, cr);

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
    rj.contig = cl->rec;
    rj.offset = offset;
    contig_notify(io, cr->rec, (reg_data *)&rj);

    /* Merge lists */
    contig_register_join(io, cr->rec, cl->rec);

    /* Destroy the old right contig */
    contig_destroy(io, cr->rec);

    /* Notify left of join */
    rl.job = REG_LENGTH;
    rl.length = cl->end - cl->start + 1;
    contig_notify(io, cl->rec, (reg_data *)&rl);

    cache_flush(io);

    return 0;
}
