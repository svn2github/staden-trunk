/*
 * File:
 * Version: join.c
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: Join editor utilities
 *
 */

#include <stdio.h>

#include "edUtils.h"
#include "align.h"
#include "undo.h"
#include "xalloc.h"
#include "io_handle.h"
#include "io-reg.h"
#include "dna_utils.h"
#include "IO.h"
#include "IO2.h"
#include "contigEditor.h"
#include "gap_globals.h"
#include "align_lib.h"
#include "hash_lib.h"
#include "notes.h"

int dojoin(GapIO *io, int lnconl, int lnconr, int relx);

/* Add 'num' pads into the consensus for editor 'xx' at position 'pos'. */
#define BLOCK_SIZE 20
static void add_pads(EdStruct *xx, int pos, int num)
{
    int rem;
    int chunks;
    int k;
    char pads[BLOCK_SIZE+1] = "********************";

    if (num < 0)
	return;

    rem = num % BLOCK_SIZE;
    chunks = num / BLOCK_SIZE;
    for (k = 0; k < chunks; k++)
	insertBasesConsensus(xx, pos, BLOCK_SIZE, pads);
    if (rem)
	insertBasesConsensus(xx, pos, rem, pads);
}

int align_contigs (OVERLAP *overlap, int fixed_left, int fixed_right) {

#define DO_BLOCKS 100
#define OK_MATCH 80
#define TOO_LONG_FOR_2ND_TRY 10000
    int ierr;
    ALIGN_PARAMS *params;
    Hash *h;
    int max_seq, word_len, max_matches;
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

#if 0
static int align(EdStruct *xx0, int pos0, int len0,
		 EdStruct *xx1, int pos1, int len1)
{

    char *ol0,*ol1;
    int old_def_conf0 = xx0->default_conf_n;
    int old_def_conf1 = xx1->default_conf_n;
    OVERLAP *overlap;
    int ierr;
    int left0, left1;
    char name0[100];
    char name1[100];
    int i,pads;
    char PAD_SYM = '.';

    vfuncheader("Align contigs (join editor)");

    /* Memory allocation */
    ol0 = (char *) xmalloc(len0+1);
    ol1 = (char *) xmalloc(len1+1);

    /* Compute the consensus */
    DBcalcConsensus(xx0,pos0,len0,ol0,NULL,BOTH_STRANDS);
    DBcalcConsensus(xx1,pos1,len1,ol1,NULL,BOTH_STRANDS);

    if (NULL == (overlap = create_overlap())) return -1;
    init_overlap (overlap, ol0, ol1, len0, len1);

    if(1 > (ierr =  align_contigs (overlap))) {
	xfree(ol0);
	xfree(ol1);
	destroy_overlap(overlap);
	return -1;
    }
    
    /* Display the alignment. */

    sprintf(name0, "#%d", xx0->DBi->DB[1].number);
    sprintf(name1, "#%d", xx1->DBi->DB[1].number);
    overlap->seq1_out[overlap->right+1] = 0;
    overlap->seq2_out[overlap->right+1] = 0;
    left0 = left1 = 0;
    if(overlap->left1 == 0) left0 += overlap->left2;
    if(overlap->left2 == 0) left1 += overlap->left1;

    list_alignment(&overlap->seq1_out[overlap->left],
		   &overlap->seq2_out[overlap->left], 
		   name0, name1, pos0+left0, pos1+left1, "");

    xx1->displayPos -= left0;
    xx0->displayPos -= left1;
    xx0->link->lockOffset = xx1->displayPos - xx0->displayPos;

    /* do the edits */


    xx0->default_conf_n = -1;
    xx1->default_conf_n = -1;

    for(i=overlap->left,pads=0;i<=overlap->right;i++) {
	if(overlap->seq1_out[i] != PAD_SYM) {
	    if(pads) {
		add_pads(xx0, pos0-left1+i-pads, pads);
		pads = 0;
	    }
	}
	else {
	    pads++;
	}
    }
    if(pads) {
	add_pads(xx0, pos0-left1+i-pads-1, pads);
    }
    for(i=overlap->left,pads=0;i<=overlap->right;i++) {
	if(overlap->seq2_out[i] != PAD_SYM) {
	    if(pads) {
		add_pads(xx1, pos1-left0+i-pads, pads);
		pads = 0;
	    }
	}
	else {
	    pads++;
	}
    }
    if(pads) {
	add_pads(xx1, pos1-left0+i-pads-1, pads);
    }

    xx0->default_conf_n = old_def_conf0;
    xx1->default_conf_n = old_def_conf1;

    xfree(ol0);
    xfree(ol1);
    destroy_overlap(overlap);
    return(0);
}
#endif

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

static int align(EdStruct *xx0, int pos0, int len0,
		 EdStruct *xx1, int pos1, int len1,
		 int fixed_left, int fixed_right)
{

    char *ol0,*ol1, *cons0, *cons1;
    int old_def_conf0 = xx0->default_conf_n;
    int old_def_conf1 = xx1->default_conf_n;
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
    DBcalcConsensus(xx0,pos0,len0,ol0,NULL,BOTH_STRANDS);
    DBcalcConsensus(xx1,pos1,len1,ol1,NULL,BOTH_STRANDS);

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

	sprintf(name0, "%d", xx0->DBi->DB_contigNum);
	sprintf(name1, "%d", xx1->DBi->DB_contigNum);
	cexpand(ol0+left0, ol1+left1, len0, len1,
		exp0, exp1, &exp_len0, &exp_len1, 
		ALIGN_J_SSH | ALIGN_J_PADS, S);
	list_alignment(exp0, exp1, name0, name1, pos0, pos1, "");

	xfree(exp0);
	xfree(exp1);
    }


    /*************************************************************************/
    /* Now actually make the edits, keeping track of old and new pads. */
    openUndo(DBI(xx0));
    openUndo(DBI(xx1));

    xx0->default_conf_n = -1;
    xx1->default_conf_n = -1;
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
    xx0->default_conf_n = old_def_conf0;
    xx1->default_conf_n = old_def_conf1;
    /*************************************************************************/

    closeUndo(xx1, DBI(xx1));
    closeUndo(xx0, DBI(xx0));

    xfree(res);

    xx0->default_conf_n = old_def_conf0;
    xx1->default_conf_n = old_def_conf1;

    xfree(ol0);
    xfree(ol1);
    xfree(dp0);
    xfree(dp1);
    destroy_overlap(overlap);

    return(0);
}

#if 0
static int align_old(EdStruct *xx0, int pos0, int len0,
		 EdStruct *xx1, int pos1, int len1)
{

    char *ol0,*ol1;
    int  *depad_to_pad0_m, *depad_to_pad1_m;
    int  *depad_to_pad0,   *depad_to_pad1;
    align_int *res, *S;
    int old_def_conf0 = xx0->default_conf_n;
    int old_def_conf1 = xx1->default_conf_n;
    int off0 = 0, off1 = 0;
    int left0 = 0, left1 = 0;

    vfuncheader("Align contigs (join editor)");

    /* Memory allocation */
    ol0 = (char *) xmalloc(len0+1);
    ol1 = (char *) xmalloc(len1+1);
    depad_to_pad0 = depad_to_pad0_m = (int *)xmalloc((len0+1) * sizeof(int));
    depad_to_pad1 = depad_to_pad1_m = (int *)xmalloc((len1+1) * sizeof(int));
    S = res = (align_int *) xmalloc((len0+len1+1)*sizeof(align_int));

    /* Compute the consensus */
    DBcalcConsensus(xx0,pos0,len0,ol0,NULL,BOTH_STRANDS);
    DBcalcConsensus(xx1,pos1,len1,ol1,NULL,BOTH_STRANDS);

    /* Strip the pads from the consensus */
    depad_seq(ol0, &len0, depad_to_pad0);
    depad_seq(ol1, &len1, depad_to_pad1);

    /* Do the actual alignment */
    (void)calign(ol0, ol1, len0, len1,
		 NULL, NULL, NULL, NULL,
		 0, 0, gopenval, gextendval, 3, 0, res);

    /* Clip left end */
    if (*S != 0) {
	/* Pad at start, so shift contigs */
	if (*S < 0) {
	    left0 = -*S; /* used for display only */
	    depad_to_pad0 += -*S;
	    off0 = depad_to_pad0[0];
	    xx1->displayPos -= off0;
	    pos0 += off0;
	    len0 -= off0;
	} else {
	    left1 = *S; /* used for display only */
	    depad_to_pad1 += *S;
	    off1 = depad_to_pad1[0];
	    xx0->displayPos -= off1;
	    pos1 += off1;
	    len1 -= off1;
	}
	S++;
	xx0->link->lockOffset = xx1->displayPos - xx0->displayPos;
    }

    /* Clip right end */
    {
	int i = 0, j = 0, op;
	align_int *S2 = S;

	while (i < len0 && j < len1) {
	    if ((op = *S2++) == 0)
		i++, j++;
	    else if (op > 0)
		j += op;
	    else
		i -= op;
	}
	
	len0 = i;
	len1 = j;
    }

    /* Display the alignment. */
    {
	char *exp0, *exp1;
	int exp_len0, exp_len1;
	char name0[100];
	char name1[100];

	exp0 = (char *) xmalloc(len0+len1+1);
	exp1 = (char *) xmalloc(len0+len1+1);

	sprintf(name0, "%d", xx0->DBi->DB_contigNum);
	sprintf(name1, "%d", xx1->DBi->DB_contigNum);
	cexpand(ol0+left0, ol1+left1, len0, len1,
		exp0, exp1, &exp_len0, &exp_len1, 
		ALIGN_J_SSH | ALIGN_J_PADS, S);
	list_alignment(exp0, exp1, name0, name1, pos0, pos1, "");

	xfree(exp0);
	xfree(exp1);
    }


    /*************************************************************************/
    /* Now actually make the edits, keeping track of old and new pads. */
    openUndo(DBI(xx0));
    openUndo(DBI(xx1));

    xx0->default_conf_n = -1;
    xx1->default_conf_n = -1;
    {
	int depad_pos0 = 0, depad_pos1 = 0;
	int curr_pad0;  /* Current padded position in seq 0 */
	int curr_pad1;  /* Current padded position in seq 1 */
	int extra_pads; /* Difference between padded positions */
	int last_pad0 = -1;
	int last_pad1 = -1;
	int inserted_bases0 = 0;
	int inserted_bases1 = 0;

	while (depad_pos0 < len0 || depad_pos1 < len1) {
	    if (*S < 0) {
		depad_pos0 -= *S;
	    } else if (*S > 0) {
		depad_pos1 += *S;
	    }

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
    xx0->default_conf_n = old_def_conf0;
    xx1->default_conf_n = old_def_conf1;
    /*************************************************************************/

    closeUndo(xx1, DBI(xx1));
    closeUndo(xx0, DBI(xx0));

    xfree(res);
    xfree(ol0);
    xfree(ol1);
    xfree(depad_to_pad0_m);
    xfree(depad_to_pad1_m);

    return(0);
}
#endif



/*
 * Align the two contig editor windows
 * Returns:
 *	0 - aligned ok
 *	1 - not ok
 */
int alignOverlap(EdStruct *xx[2], int fixed_left, int fixed_right)
{
    int left0,right0;
    int left1/*,right1*/;
    int length0,length1;
    int offset = editorLockedPos(xx, 1/*force recalculation*/);
    int overlapLength;
    int len0,len1;
    int ret;
    int xx0_dp, xx1_dp;

    if (! inJoinMode(xx[0])) return 1;


    /* Compute overlap position and sizes */
    length0 = DB_Length(xx[0],0);
    length1 = DB_Length(xx[1],0);
    if (offset < 0) {
	left0 = 1-offset;
	left1 = 1;

	if (fixed_left) {
	    int d = DB_RelPos(xx[1], DBI_order(xx[1])[xx[1]->cursorSeq])-1
		+ xx[1]->cursorPos - 1;
	    left0 += d;
	    left1 += d;
	}
    } else {
	left0 = 1;
	left1 = 1+offset;

	if (fixed_left) {
	    int d = DB_RelPos(xx[0], DBI_order(xx[0])[xx[0]->cursorSeq])-1
		+ xx[0]->cursorPos - 1;
	    left0 += d;
	    left1 += d;
	}
    }

    if (fixed_right) {
	length0 = DB_RelPos(xx[0], DBI_order(xx[0])[xx[0]->cursorSeq])-1
	    + xx[0]->cursorPos;
	length1 = DB_RelPos(xx[1], DBI_order(xx[1])[xx[1]->cursorSeq])-1
	    + xx[1]->cursorPos;
    }

    if (offset+length0 < length1) {
	right0 = length0;
    } else {
	right0 = length1-offset;
    }
    overlapLength = right0 - left0+1;
    if (overlapLength <= 0) return 1;

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

    xx0_dp = xx[0]->displayPos;
    xx1_dp = xx[1]->displayPos;

    if (left0 < 1 && left1 < 1) {
	xx[0]->displayPos += MAX(left0, left1)-1;
	xx[1]->displayPos += MAX(left0, left1)-1;
    }
    if (left0 < 1) {
	len0 -= 1-left0;
	xx[0]->displayPos += 1-left0;
	left0 = 1;
    }

    if (left1 < 1) {
	len1 -= 1-left1;
	xx[1]->displayPos += 1-left1;
	left1 = 1;
    }

    if (len0 > length0 - left0 + 1) {
	len0 = length0 - left0 + 1;
    }
    if (len1 > length1 - left1 + 1) {
	len1 = length1 - left1 + 1;
    }


    xx[0]->link->lockOffset = xx[1]->displayPos - xx[0]->displayPos;

    openUndo(DBI(xx[0]));
    openUndo(DBI(xx[1]));

    /* Do the actual alignment */
    ret = align(xx[0], left0, len0, xx[1], left1, len1,
		fixed_left, fixed_right);

    if (ret) {
	/* Alignment failed - put back display positions before returning */
	xx[0]->displayPos = xx0_dp;
	xx[1]->displayPos = xx1_dp;
    } else {
	/*
	 * If displayPos has changed, put it back to the original position and
	 * adjust it once more using U_adjust_display. This will cause the undo
	 * information to be stored correctly.
	 */
	if (xx0_dp != xx[0]->displayPos) {
	    int tmp = xx[0]->displayPos - xx0_dp;
	    xx[0]->displayPos = xx0_dp;
	    U_adjust_display(xx[0], tmp);
	}

	if (xx1_dp != xx[1]->displayPos) {
	    int tmp = xx[1]->displayPos - xx1_dp;
	    xx[1]->displayPos = xx1_dp;
	    U_adjust_display(xx[1], tmp);
	}
    }

    closeUndo(xx[1], DBI(xx[1]));
    closeUndo(xx[0], DBI(xx[0]));

    return ret;

}


void joinDB(EdStruct *xx[2], GapIO *io) {
    int relx;
    int i;
    int cl, cr;
    enum States st[2];

    cl = DBI_contigNum(xx[0]);
    cr = DBI_contigNum(xx[1]);

    /*
     * Save an internal databases
     */
    for (i=0;i<2;i++)
        saveDB(xx[i], io, 0, 0 /* no notify (done by dojoin) */);

    relx = editorLockedPos(xx, 1/*force*/);


    st[0] = xx[0]->editorState; xx[0]->editorState = StateDown;
    st[1] = xx[1]->editorState; xx[1]->editorState = StateDown;

    if (relx<0) {
	dojoin(io, cl, cr, -relx);
    } else {
	dojoin(io, cr, cl, relx);
    }

    xx[0]->editorState = st[0];
    xx[1]->editorState = st[1];
}


/*
 * Joins contig 'lnconr' to contig 'lnconl' at offset 'relx', producing a
 * new contig in place of 'lnconl'.
 * The old 'lnconr' is removed.
 *
 * Returns 0 for success, -1 for error (currently never)
 */
int dojoin(GapIO *io, int lnconl, int lnconr, int relx) {
    int gl_r = io_crnbr(io, lnconl);
    int gl_l;
    int gr_l = io_clnbr(io, lnconr);
    int i, len, right = 0;
    f_int tmp;
    GContigs c;
    reg_length rl;
    reg_join rj;

    /* Shift readings in lnconr by relx */
    for (i = gr_l; i; i = io_rnbr(io, i)) {
	io_relpos(io, i) += relx;
	/* Also update their cached contig number */
	update_rnumtocnum(io, i, lnconl);
    }

    /* Linkup ends of contigs (but not yet in correct left to right order) */
    io_lnbr(io, gr_l) = gl_r;
    io_rnbr(io, gl_r) = gr_l;

    /* Call MERGE() to shuffle the links to the correct order */
    tmp = io_dbsize(io) - lnconl;
    merge_(&io_relpos(io,1), &io_length(io,1), &io_lnbr(io,1), &io_rnbr(io,1),
	   &tmp, &io_dbsize(io));

    /* MERGE() doesn't save the readings, so do that ourselves */
    for (right = i = io_clnbr(io, lnconl); i; right = i, i = io_rnbr(io, i)) {
	GReadings r;

	gel_read(io, i, r);
	r.left = io_lnbr(io, i);
	r.right = io_rnbr(io, i);
	r.position = io_relpos(io, i);
	gel_write(io, i, r);
    }

    /* merge contig annotation lists and update contig length */
    merge_contig_tags(io, lnconl, lnconr, relx);

    /* merge note lists */
    merge_contig_notes(io, lnconr, lnconl);

    contig_read(io, lnconr, c);
    len = c.length + relx;

    contig_read(io, lnconl, c);
    if (c.length < len)
	c.length = len;
    c.right = right;
    contig_write(io, lnconl, c);
    io_clength(io, lnconl) = c.length;
    io_crnbr(io, lnconl) = c.right;

    /*
     * The order of notifications here is crucial.
     *
     * We need to firstly notify the right contig that it's been joined to
     * the left contig.
     *
     * Merge the contig registration lists (copy right into left).
     *
     * Then we delete the right contig. Note that this may imply that the
     * left contig number changes, so we convert to a gel number during the
     * duration of this call.
     *
     * Finally we then tell the left contig that it's length has changed.
     */

    /* Notify right of join */
    rj.job = REG_JOIN_TO;
    rj.contig = lnconl;
    rj.offset = relx;
    contig_notify(io, lnconr, (reg_data *)&rj);

    /* Merge lists */
    contig_register_join(io, lnconr, lnconl);
        
    /* Delete old contig (lnconr) */
    gl_l = io_clnbr(io, lnconl);
    /* unlink so that io_delete_contig does not renumber, incase we're
     * deleting the last contig */
    contig_read(io, lnconr, c);
    io_clnbr(io, lnconr) = io_crnbr(io, lnconr) = 0;
    c.left = c.right = 0;
    contig_write(io, lnconr, c);
    io_delete_contig(io, lnconr);
    lnconl = rnumtocnum(io, gl_l);

    /* Notify left of join */
    rl.job = REG_LENGTH;
    rl.length = c.length;
    contig_notify(io, lnconl, (reg_data *)&rl);

    flush2t(io);

    return 0;
}
