#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "fort.h"
#include "misc.h"
#include "xalloc.h"
#include "dna_utils.h"
#include "gap_hash.h"
#include "hash_lib.h"

void set_hash8_lookup(void) {
    int i;
    /* this is to keep gap_globals happy: FIXME */
    i = 0;
}

int cmpseq_ (
		   int   *job,		/* the task to perform */
	           char *sense,		/* the orientation for seq2 */
		   int *min_mat,	/* the minimum match length */
		   int *seq1_match,	/* positions of matches in seq1 */
		   int *seq2_match,	/* positions of matches in seq2 */
		   int *len_match,	/* length of matches */
		   int *max_matches,	/* maximum number of matches */
		   char *seq1,		/* seq1 */
		   char *seq2,		/* seq2 */
		   int *seq1_len, 	/* size of seq1 and its hash array */
		   int *seq2_len,	/* size of seq2 */
	           int_fl seq1_l,       /* fortran string length for seq1 */
	           int_fl seq2_l        /* fortran string length for seq2 */

		   ) {

/*	Main interface to the sequence comparison method using hashing.

	Jobs are:
	
	1: allocate the main storage space (Generally the max required for any run)
	2: hash one sequence
	3: compare a sequence to the current hash table
	4: add a new bit of sequence to the hash table NEVER USED !!!!
	5: allocate, hash and compare 2 sequences NEVER USED !!!!
	6: deallocate the space

	Case 1 allocate all storage ( NOTE ALL CASES EXCEPT 5 ASSUME STORAGE
	                              ALREADY ALLOCATED )
	Case 2 hash seq1.
	Case 3 compares seq2 to the hash table.
	Case 4 adds seq2 at the end of the hash tables for seq1. NEVER USED!!

*/

    static Hash *h=NULL;
    switch (*job) {

    case (1):

	if ( init_hash8n ( *seq1_len, *seq2_len, 
			  8, *max_matches, *min_mat, 1, &h )) {
	    free_hash8n(h);
	    return -2;
	}

	return 0;
	break;
	
    case (2):
	
	h->seq1_len = *seq1_len;
	h->seq1 = seq1;
	if ( hash_seqn ( h, 1)) {
	    verror(ERR_WARN, "hash_seqn", "first sequence too short");
	    return -1;
	}
	(void) store_hashn ( h );
	return 0;
	break;
	
    case (3):
	
	h->seq1 = seq1;
	h->seq1_len = *seq1_len;
	h->seq2 = seq2;
	h->seq2_len = *seq2_len;
	if ( hash_seqn ( h, 2)) {
	    verror(ERR_WARN, "hash_seqn", "second sequence too short");
	    return -1;
	}
	
	return compare_seqs ( h, seq1_match, seq2_match, len_match);
	break;

    case (4):
	
	verror(ERR_WARN, "cmpseq", "illegal option 4");
	return -1;
	break;
	
    case (5):
	
	verror(ERR_WARN, "cmpseq", "illegal option 5");
	return -1;
	break;
	
    case (6):
	
	free_hash8n ( h );
	return 0;
	break;
	
	default:
	
	verror(ERR_WARN, "cmpseq", "unknown job %d", *job);
	return -2;
	break;
	
    }
}

/************************************************************/

int repeat_search (
	           int mode,		/* 1=f, 2=r, 3=b */
		   int min_match,	/* the minimum match length */
		   int **seq1_match,	/* positions of matches in seq1 */
		   int **seq2_match,	/* positions of matches in seq2 */
		   int **len_match,	/* length of matches */
		   int max_mat,		/* maximum number of matches */
		   char *seq1,		/* seq1 */
		   int seq1_len, 	/* size of seq1 and its hash array */
		   int *num_f_matches,
		   int *num_r_matches
		   ) {

    int n_matches,seq2_len,max_matches,nres;
    char *seq2,sense;
    Hash *h;
    char *depadded_seq;
    int depadded_len;
    int *depad_to_pad;
    int i;

    /* Depad sequence */
    if (NULL == (depad_to_pad = (int *)xmalloc(sizeof(int) * seq1_len)))
	return -1;
    if (NULL == (depadded_seq = (char *)xmalloc(seq1_len+1))) {
	xfree(depad_to_pad);
	return -1;
    }
    copy_seq(depadded_seq, seq1, seq1_len);
    depadded_len = seq1_len;
    depad_seq(depadded_seq, &depadded_len, depad_to_pad);
    seq1 = depadded_seq;
    seq1_len = depadded_len;

    max_matches = max_mat;
    seq2_len = seq1_len;
    seq2 = NULL;

    if ( init_hash8n ( seq1_len, seq2_len, 
		      8, max_mat, min_match, 1, &h )) {
	free_hash8n(h);
	xfree(depadded_seq);
	xfree(depad_to_pad);
	return -2;
    }
	
    h->seq1 = seq1;
    h->seq1_len = seq1_len;

    if ( hash_seqn ( h, 1)) {
	verror(ERR_WARN, "hash_seqn", "sequence too short");
	xfree(depadded_seq);
	xfree(depad_to_pad);
	return -1;
    }
    (void) store_hashn ( h );

    if ( ! (seq2 = (char *) xmalloc ( sizeof(char)*(seq1_len) ))) {
	free_hash8n ( h );
	xfree(depadded_seq);
	xfree(depad_to_pad);
	return -1;
    }
    
    (void) copy_seq ( seq2, seq1, seq1_len );
    
    h->seq2 = seq2;
    h->seq2_len = seq2_len;
    *num_f_matches = 0;
    nres = 0;
    n_matches = 0;
    if ( mode & 1 ) {

	if ( hash_seqn ( h, 2)) {
	    verror(ERR_WARN, "hash_seqn", "sequence too short");
	    free_hash8n ( h );
	    if (seq2) xfree(seq2);
	    xfree(depadded_seq);
	    xfree(depad_to_pad);
	    return -1;
	}
	sense = 'f';
	n_matches = reps ( h, seq1_match, seq2_match, len_match, 0, sense);
	*num_f_matches = n_matches;
	nres += n_matches;

    }

    *num_r_matches = 0;
    if ( mode & 2 )  {

	(void) complement_seq(seq2, seq2_len);

	if ( hash_seqn ( h, 2)) {
	    verror(ERR_WARN, "hash_seqn", "sequence too short");
	    free_hash8n ( h );
	    if (seq2) xfree(seq2);
	    xfree(depadded_seq);
	    xfree(depad_to_pad);
	    return -1;
	}

	sense = 'r';
	n_matches = reps ( h, seq1_match, seq2_match, len_match, nres, sense);
	*num_r_matches = n_matches;
	n_matches += nres;
    }

    /* Remap depadded hits to padded positions */
    for (i = 0; i < n_matches; i++) {
	int p1, p2, p1_end;
	p1 = depad_to_pad[(*seq1_match)[i]];
	p2 = depad_to_pad[(*seq2_match)[i]];
	p1_end = depad_to_pad[(*seq1_match)[i]+(*len_match)[i]-1];

	(*seq1_match)[i] = p1;
	(*seq2_match)[i] = p2;
	(*len_match) [i] = p1_end - p1 + 1;
    }

    free_hash8n ( h );
    if (seq2) xfree(seq2);
    xfree(depadded_seq);
    xfree(depad_to_pad);

    return n_matches;
}
