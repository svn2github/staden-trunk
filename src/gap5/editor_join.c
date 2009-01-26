/*
 * editor_join.c:
 * Contains functions for the Join Editor "align" and "join" buttons.
 */
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
 */
static int align(edview *xx0, int pos0, int len0,
		 edview *xx1, int pos1, int len1,
		 int fixed_left, int fixed_right)
{

    char *ol0,*ol1, *cons0, *cons1;
    OVERLAP *overlap;
    int ierr;
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

    if(-1 == (ierr =  align_contigs (overlap, fixed_left, fixed_right))) {
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
	    xx1->displayPos -= off0;
	    pos0 += off0;
	    len0 -= -*S;
	} else {
	    left1 = *S; /* used for display only */
	    depad_to_pad1 += *S;
	    off1 = depad_to_pad1[0];
	    xx0->displayPos -= off1;
	    pos1 += off1;
	    len1 -= *S;
	}
	S++;
	xx0->link->lockOffset = xx1->displayPos - xx0->displayPos;
    }

    /* Clip right end */
    {
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

	sprintf(name0, "%d", xx0->cnum);
	sprintf(name1, "%d", xx1->cnum);
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
    int left1/*,right1*/;
    int length0,length1;
    int offset, ret;
    int overlapLength;
    int len0,len1;
    int xx0_dp, xx1_dp;
    edview **xx2;

    int l0, l1, r0, r1; /* contig used extents */

    if (!xx->link)
	return -1;
    xx2 = xx->link->xx;
    offset = xx2[1]->displayPos - xx2[0]->displayPos;

    /* Compute overlap position and sizes */
    consensus_valid_range(xx2[0]->io, xx2[0]->cnum, &l0, &r0);
    consensus_valid_range(xx2[1]->io, xx2[1]->cnum, &l1, &r1);
    length0 = r0-l0+1;
    length1 = r1-l1+1;

    if (offset < 0) {
	/* -------------
	 *         |||||
	 *         --------------
	 */
	left0 = l0-offset;
	left1 = l1;
	if (fixed_left) {
	    int d = xx2[1]->cursor_apos - l0;
	    left0 += d;
	    left1 += d;
	}
    } else {
	/*         --------------
	 *         |||||
	 * -------------
	 */
	left0 = l0;
	left1 = 11+offset;
	if (fixed_left) {
	    int d = xx2[0]->cursor_apos - l0;
	    left0 += d;
	    left1 += d;
	}
    }

    if (fixed_right) {
	length0 = xx2[0]->cursor_apos;
	length1 = xx2[1]->cursor_apos;
    }


    if (offset+length0 < length1) {
	right0 = length0;
    } else {
	right0 = length1-offset;
    }
    overlapLength = right0 - left0+1;
    if (overlapLength <= 0) return -1;

    len0 = len1 = overlapLength;

    /* Add on extra data either end to allow for padding */
#define XTRA_PERC 0.30
    if (!fixed_left) {
	left0 -= (int)(overlapLength * XTRA_PERC);
	left1 -= (int)(overlapLength * XTRA_PERC);
	len0  += (int)(overlapLength * XTRA_PERC);
	len1  += (int)(overlapLength * XTRA_PERC);
    }
    if (!fixed_right) {
	len0  += (int)(overlapLength * XTRA_PERC);
	len1  += (int)(overlapLength * XTRA_PERC);
    }

    xx0_dp = xx2[0]->displayPos;
    xx1_dp = xx2[1]->displayPos;

    if (left0 < 1 && left1 < 1) {
	xx2[0]->displayPos += MAX(left0, left1)-1;
	xx2[1]->displayPos += MAX(left0, left1)-1;
    }
    if (left0 < 1) {
	len0 -= 1-left0;
	xx2[0]->displayPos += 1-left0;
	left0 = 1;
    }

    if (left1 < 1) {
	len1 -= 1-left1;
	xx2[1]->displayPos += 1-left1;
	left1 = 1;
    }

    if (len0 > length0 - left0 + 1) {
	len0 = length0 - left0 + 1;
    }
    if (len1 > length1 - left1 + 1) {
	len1 = length1 - left1 + 1;
    }

    if (len0 <= 0 || len1 <= 0)
	return -1;

    xx->link->lockOffset = xx2[1]->displayPos - xx2[0]->displayPos;

    /* Do the actual alignment */
    ret = align(xx2[0], left0, len0, xx2[1], left1, len1,
		fixed_left, fixed_right);

    if (ret) {
	/* Alignment failed - put back display positions before returning */
	xx2[0]->displayPos = xx0_dp;
	xx2[1]->displayPos = xx1_dp;
    } else {
	if (xx0_dp != xx2[0]->displayPos) {
	    xx2[0]->displayPos = xx0_dp;
	}

	if (xx1_dp != xx2[1]->displayPos) {
	    xx2[1]->displayPos = xx1_dp;
	}
    }

    xx->link->lockOffset = xx2[1]->displayPos - xx2[0]->displayPos;

    xx2[0]->refresh_flags = ED_DISP_ALL;
    edview_redraw(xx2[0]);

    xx2[1]->refresh_flags = ED_DISP_ALL;
    edview_redraw(xx2[1]);

    return ret;
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
    int offset, binp_id;

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

    //offset += contig_offset(io, &cr) - contig_offset(io, &cl);

    cache_flush(io);

    /* Force joins at the top-level IO */
    while (io->base)
	io = io->base;

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

    cache_decr(io, binp);
    cache_decr(io, binl);
    cache_decr(io, binr);
    cache_decr(io, cl);
    cache_decr(io, cr);

    /* Destroy the old right contig */
    contig_destroy(io, cr->rec);

    cache_flush(io);

    return 0;
}

#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

static int break_contig_move_bin(GapIO *io, bin_index_t *bin,
				 contig_t *cfrom, int pfrom,
				 contig_t *cto,   int pto,
				 int child_no) {
    /* Add to */
    if (pto == cto->rec) {
	/* Parent is a contig */
	if (bin->rec != cto->bin)
	    printf("Destroy old bin for contig %d (bin %d). New root=%d\n",
		   cto->rec, cto->bin, bin->rec);
	cto->bin = bin->rec;
	cto->start = 1;
	cto->end = bin->size;

	bin->parent = cto->rec;
	bin->parent_type = GT_Contig;
	bin->flags |= BIN_BIN_UPDATED;

    } else {
	/* Parent is a bin */
	bin_index_t *pb;

	if (!(pb = get_bin(io, pto)))
	    return -1;
	if (!(pb = cache_rw(io, pb)))
	    return -1;

	pb->child[child_no] = bin->rec;
	pb->flags |= BIN_BIN_UPDATED;

	bin->parent = pto;
	bin->parent_type = GT_Bin;
	bin->flags |= BIN_BIN_UPDATED;
    }

    /* Remove from: NB it may not exist? */
    if (pfrom == cfrom->rec) {
	/* Parent is a contig */
	if (cfrom->bin != bin->rec) {
	    fprintf(stderr, "pfrom incorrect\n");
	    return -1;
	}

	cfrom->bin = 0;
    } else if (pfrom > 0) {
	/* Parent is a bin */
	bin_index_t *pb;

	if (!(pb = get_bin(io, pfrom)))
	    return -1;
	if (!(pb = cache_rw(io, pb)))
	    return -1;

	if (pb->child[0] != bin->rec && pb->child[1] != bin->rec) {
	    fprintf(stderr, "pfrom incorrect\n");
	    return -1;
	}

	if (!(pb = cache_rw(io, pb)))
	    return -1;
	
	if (pb->child[0] == bin->rec)
	    pb->child[0] = 0;
	else
	    pb->child[1] = 0;
	pb->flags |= BIN_BIN_UPDATED;
    }

    return 0;
}

/*
 * Given ranges contained within a bin this makes sure that all sequences
 * referred to in these ranges have their parent listed as the new bin.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int break_contig_reparent_seqs(GapIO *io, bin_index_t *bin) {
    int i, nr = bin->rng ? ArrayMax(bin->rng) : 0;

    for (i = 0; i < nr; i++) {
	range_t *r = arrp(range_t, bin->rng, i);
	seq_t *seq = (seq_t *)cache_search(io, GT_Seq, r->rec);
	if (seq->bin != bin->rec) {
	    seq = cache_rw(io, seq);
	    seq->bin = bin->rec;
	}
    }

    return 0;
}

/*
 * A recursive break contig function.
 * bin_num	The current bin being moved or split.
 * pos		The contig break point.
 * offset	The absolute positional offset of this bin in original contig
 * pleft	The parent bin/contig record num in the left new contig
 * pright	The parent bin/contig record num in the right new contig
 * child_no     0 or 1 - whether this bin is the left/right child of its parent
 */
static int break_contig_recurse(GapIO *io,
				contig_t *cl, contig_t *cr,
				int bin_num, int pos, int offset,
				int level, int pleft, int pright,
				int child_no, int complement,
				int *left_end, int *right_start) {
    int i, f_a, f_b, rbin;
    bin_index_t *bin = get_bin(io, bin_num), *bin_dup ;
    int bin_min, bin_max;

    cache_incr(io, bin);

    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }

    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    printf("%*sBreak offset %d pos %d => test bin %d: %d..%d\n",
	   level*4, "",
	   offset, pos, bin->rec,
	   NMIN(bin->start_used, bin->end_used),
	   NMAX(bin->start_used, bin->end_used));

    bin = cache_rw(io, bin);

    bin_min = bin->rng ? NMIN(bin->start_used, bin->end_used) : offset;
    bin_max = bin->rng ? NMAX(bin->start_used, bin->end_used) : offset;

    /*
     * Add to right parent if this bin is to the right of pos,
     * or if the used portion is to the right and we have no left child.
     *
     * FIXME: Not a valid assumption!
     * The used portion of a bin is not a placeholder for the used portion
     * of all the the children beneath it. Therefore if the used portion of
     * this bin is > pos (and we have no left child) it still doesn't mean
     * that the absolute positions of the used portion of the right child
     * won't be < pos.
     */
    if (offset >= pos /*|| (bin_min >= pos && !bin->child[0])*/) {
	printf("%*sADD_TO_RIGHT pl=%d pr=%d\n", level*4, "", pleft, pright);
	if (0 != break_contig_move_bin(io, bin,
				       cl, pleft, cr, pright, 
				       child_no))
	    return -1;

	if (*right_start > bin_min)
	    *right_start = bin_min;

	cache_decr(io, bin);
	return 0;
    }

    /*
     * Add to left parent if this bin is entirely to the left of pos,
     * or if the used portion is to the left and we have no right child.
     */
    if (offset + bin->size < pos /*|| (bin_max < pos && !bin->child[1])*/) {
	printf("%*sADD_TO_LEFT\n", level*4, "");

	//if (0 != break_contig_move_bin(io, bin, cr, pright, cl, pleft, child_no))
	//return -1;
	if (*left_end < bin_max)
	    *left_end = bin_max;

	cache_decr(io, bin);
	return 0;
    }

    /*
     * Nominally the bin overlaps both left and right and so needs duplicating.
     * There are cases though at the roots of our trees where duplicating is
     * unnecessary as it leads to empty bins at the root. In this case
     * we skip creating a duplicate for the right, or alternatively steal
     * the left root bin and use that instead.
     *
     * Similarly the range_t array will either be left where it is, moved to
     * the right contig, or split in half (creating a new one for the right).
     */
    if (pright != cr->rec ||
	(bin->rng && NMAX(bin->start_used, bin->end_used) >= pos)) {
	//printf("NMAX=%d >= %d\n", NMAX(bin->start_used, bin->end_used), pos);

	rbin = 0;

	/* Possibly steal left contig's bin */
	if (pleft == cl->rec && NMIN(bin->start_used, bin->end_used) >= pos) {
#if 0
	    /* Currently this doesn't always work */
	    if (bin->child[1]) {
		bin_index_t *ch = get_bin(io, bin->child[1]);
		if (NMIN(ch->pos, ch->pos + ch->size-1) >= pos) {
		    rbin = cl->bin;
		    cl->bin = bin->child[0];
		}
	    }
#else
	    pleft = bin->rec;
#endif
	} else {
	    pleft = bin->rec;
	}

	/* Create new bin, or use root of contig if it's unused so far */
	if (!rbin && pright == cr->rec) {
	    rbin = cr->bin;
	}

	/* Otherwise we genuingly need a duplicate */
	if (!rbin)
	    rbin = bin_new(io, 0, 0, 0, GT_Bin);

	/* Initialise with duplicate values from left bin */
	bin_dup = get_bin(io, rbin);
	bin_dup = cache_rw(io, bin_dup);
	bin_dup->size = bin->size;
	bin_dup->pos = bin->pos;
	bin_dup->parent = pright;
	bin_dup->parent_type = (pright == cr->rec ? GT_Contig : GT_Bin);
	bin_dup->flags = bin->flags | BIN_BIN_UPDATED;
	bin_dup->start_used = bin->start_used;
	bin_dup->end_used = bin->end_used;

	/*
	 * Shift bin by offset if it's the contig root.
	 * It'll be shifted back by the correct amount later.
	 */
	if (pright == cr->rec) {
	    printf("moving root bin via offset=%d comp=%d\n", offset, complement);
	    bin_dup->pos += offset;
	}

	printf("%*sCreated dup for right, rec %d\n", level*4,"", bin_dup->rec);
	break_contig_move_bin(io, bin_dup, cl, 0, cr, pright, child_no);
	pright = bin_dup->rec;
    } else {
	bin_dup = NULL;
	pleft = bin->rec;
    }

    if (!bin->rng) {
	/* Empty bin */
	printf("%*sEMPTY range\n", level*4, "");
	bin->start_used = bin->end_used = 0;
	bin->flags |= BIN_BIN_UPDATED;
	if (bin_dup) {
	    bin_dup->start_used = bin_dup->end_used = 0;
	    bin_dup->flags |= BIN_BIN_UPDATED;
	}
	    
    } else if (NMIN(bin->start_used, bin->end_used) >= pos) {
	/* Move range to right contig */
	printf("%*sDUP %d, MOVE Array to right\n", level*4, "", bin_dup->rec);

	bin_dup->rng = bin->rng;
	bin_dup->rng_rec = bin->rng_rec;
	if (bin_dup->rng_rec)
	    bin_dup->flags |= BIN_RANGE_UPDATED;

	if (bin->rec != bin_dup->rec) {
	    bin->rng = NULL;
	    bin->rng_rec = 0;
	    bin->flags |= BIN_BIN_UPDATED;
	}

	if (*right_start > NMIN(bin->start_used, bin->end_used))
	    *right_start = NMIN(bin->start_used, bin->end_used);

	bin->start_used = bin->end_used = 0;
	break_contig_reparent_seqs(io, bin_dup);

    } else if (NMAX(bin->start_used, bin->end_used) < pos) {
	/* Range array already in left contig, so do nothing */
	printf("%*sMOVE Array to left\n", level*4, "");

	if (*left_end < NMAX(bin->start_used, bin->end_used))
	    *left_end = NMAX(bin->start_used, bin->end_used);

	if (bin_dup)
	    bin_dup->start_used = bin_dup->end_used = 0;

    } else {
	/* Range array covers pos, so split in two */
	int i, j, n;
	int lmin = bin->size, lmax = 0, rmin = bin->size, rmax = 0;
	printf("%*sDUP %d, SPLIT array\n", level*4, "", bin_dup->rec);

	bin->flags |= BIN_RANGE_UPDATED;
	bin_dup->flags |= BIN_RANGE_UPDATED;

	bin_dup->rng = ArrayCreate(sizeof(range_t), 0);
	n = ArrayMax(bin->rng);
	for (i = j = 0; i < n; i++) {
	    range_t *r = arrp(range_t, bin->rng, i), *r2;
	    seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	    int cstart; /* clipped sequence positions */

	    if ((s->len < 0) ^ complement) {
		cstart = NMAX(r->start, r->end) - (s->right-1);
		/* cend   = NMAX(r->start, r->end) - (s->left-1); */
	    } else {
		cstart = NMIN(r->start, r->end) + s->left-1;
		/* cend   = NMIN(r->start, r->end) + s->right-1; */
	    }
	    
	    if (cstart >= pos)  {
		r2 = (range_t *)ArrayRef(bin_dup->rng, ArrayMax(bin_dup->rng));
		*r2 = *r;
		if (rmin > r->start) rmin = r->start;
		if (rmin > r->end)   rmin = r->end;
		if (rmax < r->start) rmax = r->start;
		if (rmax < r->end)   rmax = r->end;
	    } else {
		if (lmin > r->start) lmin = r->start;
		if (lmin > r->end)   lmin = r->end;
		if (lmax < r->start) lmax = r->start;
		if (lmax < r->end)   lmax = r->end;

		if (j != i) {
		    r2 = arrp(range_t, bin->rng, j);
		    *r2 = *r;
		}
		j++;
	    }
	}

	ArrayMax(bin->rng) = j;

	printf("%d seqs => left=%d, right=%d\n",
	       n, ArrayMax(bin->rng), bin_dup ? ArrayMax(bin_dup->rng) : 0);

	if (bin_dup)
	    break_contig_reparent_seqs(io, bin_dup);

	if (lmin < lmax) {
	    bin->start_used     = lmin;
	    bin->end_used       = lmax;
	} else {
	    /* No data left in bin */
	    bin->start_used = 0;
	    bin->end_used = 0;
	}

	printf("%*sLeft=>%d..%d right=>%d..%d\n", level*4, "",
	       lmin, lmax, rmin, rmax);

	if (bin_dup) {
	    if (rmin < rmax) {
		bin_dup->start_used = rmin;
		bin_dup->end_used   = rmax;
	    } else {
		/* No data moved in bin */
		bin_dup->start_used = 0;
		bin_dup->end_used   = 0;
	    }
	}

	if (lmin < lmax) {
	    if (*left_end < NMAX(bin->start_used, bin->end_used))
		*left_end = NMAX(bin->start_used, bin->end_used);
	}
	if (bin_dup && rmin < rmax) {
	    if (*right_start > NMIN(bin_dup->start_used, bin_dup->end_used))
		*right_start = NMIN(bin_dup->start_used, bin_dup->end_used);
	}
    }


    /* Recurse */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;

	ch = get_bin(io, bin->child[i]);
	if (0 != break_contig_recurse(io, cl, cr, bin->child[i], pos,
				      NMIN(ch->pos, ch->pos + ch->size-1),
				      level+1, pleft, pright,
				      i, complement,
				      left_end, right_start))
	    return -1;
    }

    cache_decr(io, bin);
    if (bin_dup)
	cache_decr(io, bin_dup);

    return 0;
}

/*
 * Looks for redundant bins at the root containing no data and just a single
 * child.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int remove_redundant_bins(GapIO *io, contig_t *c) {
    int bnum;

    if (!(c = cache_rw(io, c)))
	return -1;

    for (bnum = c->bin; bnum;) {
	bin_index_t *bin = get_bin(io, bnum);
	if (bin->rng || (bin->child[0] && bin->child[1]))
	    break;

	/* Empty */
	c->bin = bin->child[0] ? bin->child[0] : bin->child[1];
	printf("Remove bin %d\n", bin->rec);
	bnum = c->bin;
    }

    return 0;
}

/*
 * Breaks a contig in two such that snum is the right-most reading of
 * a new contig.
 */
int break_contig(GapIO *io, int crec, int cpos) {
    contig_t *cl = (contig_t *)cache_search(io, GT_Contig, crec);
    contig_t *cr;
    int cid;
    char cname[1024], *cname_end;
    int left_end, right_start;
    bin_index_t *bin;
    int do_comp = 0;

    contig_dump_ps(io, &cl, "/tmp/tree.ps");

    strncpy(cname, contig_get_name(&cl), 1000);
    cname_end = cname + strlen(cname);
    cid = 1;
    do {
	sprintf(cname_end, "#%d", cid++);
    } while (contig_index_query(io, cname) > 0);

    if (!(cr = contig_new(io, cname)))
	return -1;
    cr = cache_rw(io, cr);
    if (0 != contig_index_update(io, cname, strlen(cname), cr->rec))
	return -1;
    printf("Break in contig %d, pos %d\n", crec, cpos);

    printf("Existing left bin = %d, right bin = %d\n",
	   cl->bin, cr->bin);

    bin = get_bin(io, cl->bin);
    do_comp = bin->flags & BIN_COMPLEMENTED;
    cache_decr(io, bin);

    left_end = 0;
    right_start = INT_MAX;
    break_contig_recurse(io, cl, cr,
			 contig_get_bin(&cl), cpos, contig_offset(io, &cl),
			 0, cl->rec, cr->rec, 0, 0, &left_end, &right_start);

    printf("New left end = %d, right start = %d\n", left_end, right_start);

    /* Ensure start/end positions of contigs work out */
    cr->start = 1;
    cr->end = cl->end - right_start + 1;
    bin = cache_rw(io, get_bin(io, cr->bin));
    bin->pos -= right_start-1;
    if ((do_comp && !(bin->flags & BIN_COMPLEMENTED)) ||
	(!do_comp && (bin->flags & BIN_COMPLEMENTED))) {
	bin->flags ^= BIN_COMPLEMENTED;
    }

    /*
    cr->start = right_start;
    cr->end = cl->end;
    bin = cache_rw(io, get_bin(io, cr->bin));
    bin->pos+=10;
    */

    cl = cache_rw(io, cl);
    cl->end = left_end;

    remove_redundant_bins(io, cl);
    remove_redundant_bins(io, cr);

    cache_flush(io);

    printf("Final left bin = %d, right bin = %d\n",
	   cl->bin, cr->bin);

    contig_dump_ps(io, &cl, "/tmp/tree_l.ps");
    contig_dump_ps(io, &cr, "/tmp/tree_r.ps");

    return 0;
}
