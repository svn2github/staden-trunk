#include <stdio.h>

typedef void ALIGN_PARAMS;

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
 * word_length
 */
int count_mismatches(int word_length, int maxmis,
		     char *seq1, int len1, int pos1,
		     char *seq2, int len2, int pos2) {
    int p1, p2;
    int nmis = 0;

    /* FIXME: add boundary checks for , p1<0, p2<0, p1>len1, p2>len2 */

    /* Count backward matches */
    for (p1 = pos1-1, p2 = pos2-1; p2 >= 0; p1--) {
	if (p1 < 0)
	    return maxmis+1;

	if (seq1[p1] != seq2[p2]) {
	    if (++nmis > maxmis)
		return nmis;
	}
	p2--;
    }

    /* Count forward matches */
    for (p1 = pos1, p2 = pos2; p2 < len2; p1++) {
	if (p1 >= len1)
	    return maxmis+1;

	if (seq1[p1] != seq2[p2]) {
	    if (++nmis > maxmis)
		return nmis;
	}
	p2++;
    }

    return nmis;
}


extern int dna_hash8_lookup[256];

int compare_h(Hash *h, int min_match, int self_match, double perc) {
    int nrw, ncw, word, pw2, pw1, j;
    int maxmis = (1 - perc) * h->seq2_len;

    if(h->seq1_len < min_match) return -4; 
    if(h->seq2_len < min_match) return -4; 
    
    nrw = h->seq2_len - h->word_length + 1;
    
    /* 	loop for all (nrw) complete words in values2 */
    for (pw2=0;pw2<nrw;pw2++) {
 	if ((word = h->values2[pw2]) == -1)
	    continue;

	if ((ncw = h->counts[word]) == 0)
	    continue;

	for (j=0,pw1=h->last_word[word];j<ncw;j++) {
	    int nmismatch;
	    nmismatch = count_mismatches(h->word_length, maxmis,
					 h->seq1, h->seq1_len, pw1,
					 h->seq2, h->seq2_len, pw2);
	    if (nmismatch <= maxmis) {
		/* FIXME: deal with self matches */
		/* FIXME: deal with multiple instances of same match */

		/* printf("Hashm at %d (%d/%d)\n", pw1-pw2, maxmis, nmismatch); */
		return 1;
	    }
	    pw1 = h->values1[pw1];
	}
    }

    return 0;
}

#include "search_utils.h"
int compare2(char *seq1, int len1, char *seq2, int len2, double perc) {
    int nmismatch = (1 - perc) * len2;
    char *match;
    int offset = 0;

    while (match = pstrnstr_inexact(seq1, len1, seq2, len2, nmismatch)) {
	/* printf("Strnm at %d\n", match - seq1 + offset); */

	/* Only allow 1 match */
	return 1;

	/*
	offset += match-seq1+1;
	len1 -= match-seq1+1;
	seq1 = match+1;
	*/
    }

    return 0;
}

int main(int argc, char **argv) {
    char seq[MAXSEQ];
    int seq_len;
    FILE *fp;
    Hash *h;
    int i;
    int iter;
    char vec[21];

    set_char_set(1);    /* 1 == DNA */
    set_dna_lookup(); 	/* general lookup and complementing */

    if (argc != 3) {
	fprintf(stderr, "Usage: %s screenfile niter\n", argv[0]);
	return 1;
    }

    if(init_hash8n(MAXSEQ, MAXPRIMER, 4 /* word length */,
		   2 /* job */, &h)) {
	fprintf(stderr, "init_hash8n failed\n");
	return 1;
    }

    fp = fopen(argv[1], "r");
    h->seq1_len = fread(seq, 1, MAXSEQ, fp);
    h->seq1 = seq;
    fclose(fp);
    if (hash_seqn(h, 1)) {
	fprintf(stderr, "hash seq1 failed\n");
    }
    store_hashn(h);

    h->seq2_len = 20;
    h->seq2 = vec;
    iter = atoi(argv[2]);
    srandom(15551);
    for (i = 0; i < iter; i++) {
	int j;
	for (j = 0; j < 20; j++) {
	    vec[j] = "ACGT"[random()%4];
	}
	vec[20] = 0;
	if (hash_seqn(h, 2)) {
	    fprintf(stderr, "hash seq2 failed\n");
	}
    
	printf("Testing seq %d : ", i);
	if (compare_h(h, 4, 1, 0.77)) {
	    printf("1 ");
	}

	/*
	if (compare2(h->seq1, h->seq1_len, h->seq2, h->seq2_len, 0.77)) {
	    printf("2 ");
	}
	*/
	
	putchar('\n');
    }

    free_hash8n(h);
    return 0;
}
