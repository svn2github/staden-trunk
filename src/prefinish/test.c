#include <stdio.h>

#include "nw_lib.h"
#include "hash_lib.h"

#define MAXSEQ 1000000
#define MAXPRIMER 30
#define MAX_MATCHES 1000

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
int count_mismatches(int maxmis,
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
 *	minmat		Minimum match fraction for indicating a 'match'
 *	skip_self	How many matches to skip past (on skip_strand only)
 *	skip_strand	Strand (0=top,1=bot) on which skip_self counts.
 *
 * Returns:
 *	0		no matches found
 *	1		matches found
 *     -1		Error
 */
int hash_compare_primer(Hash *h, char *prim, int lprim,
			double minmat, int skip_self, int skip_strand) {
    int nrw, ncw, word, pw2, pw1, j;
    int maxmis = (1 - minmat) * lprim;
    signed int last_pos = -1;
    int strand;
    char pcopy[MAXPRIMER];

    /* Sanity checks */
    if(h->seq1_len < h->word_length)
	return 0; 

    if(lprim < h->word_length)
	return 0; 
    
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
		int nmismatch;
		nmismatch = count_mismatches(maxmis,
					     h->seq1, h->seq1_len, pw1,
					     h->seq2, h->seq2_len, pw2);
		if (nmismatch <= maxmis) {
		    if (pw1-pw2 != last_pos) {
			if (self_count)
			    self_count--;
			else
			    return 1;
			last_pos = pw1 - pw2;
		    }
		}
		pw1 = h->values1[pw1];
	    }
	}

	/* Complement - done twice so prim is ultimately unmodified */
	complement_seq(pcopy, lprim);
    }

    return 0;
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
 *	minmat		Minimum match fraction for indicating a 'match'
 *	skip_self	How many matches to skip past (on skip_strand only)
 *	skip_strand	Strand (0=top,1=bot) on which skip_self counts.
 *
 * Returns:
 *	0		no matches found
 *	1		matches found
 *     -1		Error
 */
int compare_primer(char *seq1, int len1, char *prim, int lprim,
		   double minmat, int skip_self, int skip_strand) {
    Hash *h;
    int ret;

    if (init_hash8n(len1, lprim, 4 /* word length */, 0 /* job */, &h)) {
	fprintf(stderr, "init_hash8n failed\n");
	return -1;
    }
    
    h->seq1 = seq1;
    h->seq1_len = len1;

    if (hash_seqn(h, 1)) {
	fprintf(stderr, "hash seq1 failed\n");
	return -1;
    }
    store_hashn(h);

    ret = hash_compare_primer(h, prim, lprim, minmat, skip_self, skip_strand);
    
    free_hash8n(h);

    return ret;
}

char seq[MAXSEQ];
int main(int argc, char **argv) {
    int seq_len;
    FILE *fp;
    Hash *h;
    int i;
    int iter;
    char prim[21];

    set_char_set(1);    /* 1 == DNA */
    set_dna_lookup(); 	/* general lookup and complementing */

    if (argc != 3) {
	fprintf(stderr, "Usage: %s screenfile niter\n", argv[0]);
	return 1;
    }

    fp = fopen(argv[1], "r");
    seq_len = fread(seq, 1, MAXSEQ, fp);
    fclose(fp);

    if (init_hash8n(seq_len, MAXPRIMER, 4 /* word length */,
		   0 /* job */, &h)) {
	fprintf(stderr, "init_hash8n failed\n");
	return 1;
    }

    h->seq1 = seq;
    h->seq1_len = seq_len;
    if (hash_seqn(h, 1)) {
	fprintf(stderr, "hash seq1 failed\n");
	return 1;
    }
    store_hashn(h);

    iter = atoi(argv[2]);
    srandom(15551);
    for (i = 0; i < iter; i++) {
	int j;
	for (j = 0; j < 20; j++) {
	    prim[j] = "ACGT"[random()%4];
	}
	prim[20] = 0;

	printf("Testing seq %d : ", i);

	if (hash_compare_primer(h, prim, 20, 0.77, 1, 0)) {
	    printf("1 ");
	}

	if (compare_primer(seq, seq_len, prim, 20, 0.77, 1, 0)) {
	    printf("2 ");
	}

	putchar('\n');
	fflush(stdout);
    }

    free_hash8n(h);
    return 0;
}
