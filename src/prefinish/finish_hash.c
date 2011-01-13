/*
 * Includes code interfacing the 'finish' package with the sequence
 * hashing library (hash_lib.c).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "finish_hash.h"
#include "dna_utils.h"
#include "xalloc.h"
#include "misc.h"

#if 0
/*
 * Counts the number of mismatches between seq1 at position pos1 and seq2
 * position pos2. Seq2/len2 is assumed to be the shortest sequence and we
 * may step through all of this, but range checking is performed on
 * seq1 to ensure we don't over/under step the buffer.
 *
 * maxmis specifies the maximum number of mismatches we allow. If we find
 * any more than this then we bail-out early and return maxmis+1.
 * This is also returned if the seq2 isn't entirely contained within seq1.
 *
 * Returns: number of mismatches (at most maxmis+1).
 */
static int count_mismatches(int maxmis,
			    char *seq1, int len1, int pos1,
			    char *seq2, int len2, int pos2) {
    int p1, p2;
    int nmis = 0;

    /* Is seq2 contained within seq1? */
    if (pos1 - pos2 < 0 || pos1 + (len2 - pos2) >= len1)
	return maxmis + 1;

    /* Count mismatches */
    for (p1 = pos1-pos2, p2 = 0; p2 < len2; p1++, p2++) {
	if (seq1[p1] != seq2[p2]) {
	    if (++nmis > maxmis)
		return nmis;
	}
    }

    return nmis;
}
#endif

/*
 * Computes a false-priming score for this oligo, based on a weight matrix
 * used to give matches at the 3' end priority.
 *
 * End indicates which end to weight. 0 => left, 1 => right
 * If the match score is >= print_gte then the match is printed up
 * (for debugging usage.)
 *
 * Returns the match score, or 0 if none (or error).
 */
static double false_priming(int end,
			    char *seq1, int len1, int pos1,
			    char *seq2, int len2, int pos2,
			    int *perfect,
			    char *msg_buf) {
    int i;
    double score = 0;
    double max_score = 0;
    int end_len; /* length of exact match at 3' end */

    double wmatrix[50] = {
	1.2, 1.0, 1.0, 1.0, 0.9, 0.8, 0.7, 0.5, 0.5, 0.5,
	0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
	0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
	0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
	0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5
    };

    *perfect = 0;

    /*
    printf("%d/%d %d/%d\n", pos1,len1, pos2,len2);
    */

    /* If Seq2 contained within seq1? */
    pos1 -= pos2;
    seq1 += pos1;
    pos2 = 0;
    if (pos1 < 0 || pos1 + len2 >= len1)
	return 0.0;

    /* Count mismatches, ranking by position */
    if (end == 0) {
	/* left end */
	for (end_len = i = 0; i < len2; i++) {
	    double sc = wmatrix[i];
	    if (seq1[i] == seq2[i]) {
		if (i == end_len)
		    end_len++;
		score += sc;
	    }
	    max_score += sc;
	}
    } else {
	/* right end */
	for (end_len = i = len2-1; i >= 0; i--) {
	    double sc = wmatrix[len2-1-i];
	    if (seq1[i] == seq2[i]) {
		if (i == end_len)
		    end_len--;
		score += sc;
	    }
	    max_score += sc;
	}
	end_len = len2-1-end_len;
    }
    max_score += len2 * .3;
    score += end_len * .3;

    /*
     * In order to assure debugging output on different machines is the
     * same, we quantise score and max_score to 1 decimal point as the
     * remainder will be (small) rounding errors.
     */
    score = ((int)(score * 10 + 0.01)) / 10.0;
    max_score = ((int)(max_score * 10 + 0.01)) / 10.0;

    if (msg_buf) {
	sprintf(msg_buf,
		"Primer match score=%5.1f (max %5.1f) at pos %d\n"
		"    %d' %.*s %d'\n    %d' %.*s %d'\n",
		score, max_score, pos1,
		end ? 5 : 3, len2, seq1, end ? 3 : 5,
		end ? 5 : 3, len2, seq2, end ? 3 : 5);
    }

    *perfect = (score == max_score);
	

    return score;
}


extern int dna_hash8_lookup[256];

/*
 * Compares a primer sequence 'prim' of length 'lprim' against a hashed
 * sequence stored in Hash. The primer is automatically checked in both
 * directions, but the primer must consist of upper case A,C,G,T characters
 * only.
 *
 * Arguments:
 *	h		Hashed sequence
 *	prim		Primer sequence (must be uppercase)
 *	lprim		Length of primer sequence
 *	max_match	Maximum score before we reject (due to 2ndary prim)
 *	skip_self	How many matches to skip past (on skip_strand only)
 *	skip_strand	Strand (0=top,1=bot) on which skip_self counts.
 *
 * Returns:
 *	-1	Error
 *	>= 0	Score (high means strong match, low means poor match)
 */
double hash_compare_primer(Hash *h, char *prim, int lprim,
			   double max_match, int skip_self, int skip_strand) {
    int nrw, ncw, word, pw2, pw1, j;
    /* int maxmis = (1 - minmat) * lprim; */
    signed int last_pos = -1;
    int strand;
    char pcopy[FIN_MAXPRIMERLEN];
    double max_pscore = 0;
    char msg_buf[1024], best_msg_buf[1024];

    *best_msg_buf = 0;

    /* Sanity checks */
    if(h->seq1_len < h->word_length)
	return -1; 

    if(lprim < h->word_length)
	return -1; 
    
    memcpy(pcopy, prim, lprim);
    nrw = lprim - h->word_length + 1;
    
    /* Loop through both strands */
    for (strand = 0; strand < 2; strand++) {
	int self_count = strand == skip_strand ? skip_self : 0;

	/* Hash primer sequence */
	h->seq2 = pcopy;
	h->seq2_len = lprim;
	if (hash_seqn(h, 2)) {
	    fprintf(stderr, "Couldn't hash primer sequence\n");
	    return -1;
	}

	/* loop for all (nrw) complete words in values2 (primer) */
	for (pw2=0;pw2<nrw;pw2++) {
	    if ((word = h->values2[pw2]) == -1)
		continue;

	    if ((ncw = h->counts[word]) == 0)
		continue;

	    /* Check each matching word from seq1 for a 'real' match */
	    for (j=0,pw1=h->last_word[word];j<ncw;j++) {
#if 0
		int nmismatch;
#endif
		double pscore;
		int perfect;

		if (pw1 - pw2 == last_pos) {
		    pw1 = h->values1[pw1];
		    continue;
		}

		pscore = false_priming(strand ? 0 : 1,
				       h->seq1, h->seq1_len, pw1,
				       h->seq2, h->seq2_len, pw2,
				       &perfect,
				       msg_buf);

		if (pw1-pw2 != last_pos) {
		    if (self_count && perfect) {
			self_count--;
			last_pos = pw1 - pw2;
		    } else {
			/* Matches elsewhere */
			if (pscore > max_pscore) {
			    max_pscore = pscore;
			    strcpy(best_msg_buf, msg_buf);
			}
		    }
		}


#if 0
		nmismatch = count_mismatches(maxmis,
					     h->seq1, h->seq1_len, pw1,
					     h->seq2, h->seq2_len, pw2);
		if (nmismatch <= maxmis) {
		    if (fin->opts.debug[FIN_DEBUG_VPWALK] > 1)
			printf("Match strand=%d pos=%d, last_pos=%d, "
			       "self_count=%d\n",
			       strand, pw1-pw2, last_pos, self_count);
		    if (pw1-pw2 != last_pos) {
			if (self_count) {
			    self_count--;
			} else {
			    return 1;
			}
			last_pos = pw1 - pw2;
		    }
		}
#endif

		pw1 = h->values1[pw1];
	    }
	}

	/* Complement - done twice so prim is ultimately unmodified */
	complement_seq(pcopy, lprim);
    }

#if 1
    if (max_pscore >= max_match && *best_msg_buf)
	printf("%s", best_msg_buf);
#endif

    return max_pscore;
}

/*
 * This is a wrapper around hash_compare_primer for when we do not already
 * have a hashed sequence. Hashing a sequence and comparing the hashes is still
 * quicker than doing a direct comparison and it also has the benefit of
 * re-using existing code.
 *
 * Arguments:
 *	seq1		Sequence to compare against
 *	len1		Length of sequence1
 *	prim		Primer sequence (must be uppercase)
 *	lprim		Length of primer sequence
 *	max_match	Maximum score before we reject (due to 2ndary prim)
 *	skip_self	How many matches to skip past (on skip_strand only)
 *	skip_strand	Strand (0=top,1=bot) on which skip_self counts.
 *
 * Returns:
 *     -1		Error
 *	>= 0		Score for match (high means strong match)
 */
double compare_primer(char *seq1, int len1, char *prim, int lprim,
		      double max_match, int skip_self, int skip_strand) {
    Hash *h;
    double ret;
    char buf[8192], *tmp_buf;
    int allocated = 0;
    int i;

    if (len1 < 4) {
	/* shorter than word length => skip */
	return 0;
    }

    /* Pad strip */
    if (len1 > 8192) {
	if (NULL == (tmp_buf = (char *)xmalloc(len1)))
	    return -1;
	allocated = 1;
    } else {
	tmp_buf = buf;
    }

    memcpy(tmp_buf, seq1, len1);
    depad_seq(tmp_buf, &len1, NULL);

    /* Replace 'd,e,f,i with A,C,G,T again - unmasks sequence */
    for (i = 0; i < len1; i++) {
	switch(tmp_buf[i]) {
	case 'd':
	case 'D':
	    tmp_buf[i] = 'A';
	    break;
	case 'e':
	case 'E':
	    tmp_buf[i] = 'C';
	    break;
	case 'f':
	case 'F':
	    tmp_buf[i] = 'G';
	    break;
	case 'i':
	case 'I':
	    tmp_buf[i] = 'T';
	    break;
	}
    }

    if (init_hash8n(len1, lprim,
		    4 /* word length */,
		    0 /* max matches - unused */,
		    0 /* min match - unused */,
		    1 /* job */,
		    &h)) {
	fprintf(stderr, "init_hash8n failed\n");
	return -1;
    }
    
    h->seq1 = tmp_buf;
    h->seq1_len = len1;

    if (hash_seqn(h, 1)) {
	fprintf(stderr, "hash seq1 failed\n");
	return -1;
    }
    store_hashn(h);

    ret = hash_compare_primer(h, prim, lprim, max_match, skip_self,
			      skip_strand);
    
    free_hash8n(h);

    if (allocated)
	xfree(tmp_buf);

    return ret;
}

