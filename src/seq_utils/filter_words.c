#include <stdio.h>
#include <string.h>
#include "filter_words.h"
#include "dna_utils.h"

#define MATCH 1
#define MISMATCH -1.9

/*
 *      n
 * A    A
 *  C   C
 * AC   M
 *   G  G
 * A G  R
 *  CG  S
 * ACG  V
 *    T T
 * A  T W
 *  C T Y
 * AC T H
 *   GT K
 * A GT D
 *  CGT B
 * ACGT N
 */

typedef int fkey_t;

/*
 * Construct a key and a mask (allbits = 1) from a sequence.
 * If sequence is longer than Number_of_bits(fkey_t)/4 then the start of the
 * key will be clipped.
 *
 * Seq can contain the standard ambiguity codes.
 */
static fkey_t construct_key(char *seq, fkey_t *mask, int *keylen, int *keyjmp)
{
    fkey_t k = 0, m = 0;
    int seqlen = strlen(seq), i;
    int jumplen;
    char seq2[201];

    sprintf(seq2, "%.100s%.100s\n", seq, seq);
    for (i = 1; i <= seqlen; i++) {
	if (memcmp(&seq2[i], seq, seqlen) == 0) {
	    jumplen = i;
	    break;
	}
    }

    for (; *seq; seq++) {
	k = (k << 4) | ambiguity2basebit(*seq);	/* key */
	m = (m << 4) | 15;			/* mask */
    }

    if (mask)
	*mask = m;
    if (keyjmp)
	*keyjmp = jumplen;
    if (keylen)
	*keylen = seqlen;

    return k;
}

/*
 * Search for occurences of rep in seq (of length len).
 * minsize is the minimum length of consecutive 'rep' matches to use for
 * reporting.
 */
int filter_words(char *seq, char *filt, size_t len, char *rep,
		 int minsize, int maxdrop, char filter_char) {
    size_t i, j;
    fkey_t word = 0;
    int score = -1, maxscore;
    size_t start = 0, end = 0;
    fkey_t key, keymask;
    int keyjump, keylen;
    int pads = 0;

    key = construct_key(rep, &keymask, &keylen, &keyjump);

    /* Start with an entire word */
    for (i = j = 0; i < len && j < keylen-1; i++) {
	if (seq[i] == '*') {
	    pads++;
	    continue;
	}

	word = ((word << 4) | ambiguity2basebit(seq[i])) & keymask;
	/*
	 * printf("SEQ[%03d]=%c => word=%04x & %04x == %04x, %04x == %04x\n",
	 *       i, seq[i], word, key, word & key,
	 *       ~key, (signed)(word & ~key));
	 */

	j++;
    }

    /* Scan through looking for matches of a word */
    for (; i < len; i++) {
	if (seq[i] == '*'){
	    /* printf("Seq[%03d]='*'\n", i); */
	    pads++;
	    continue;
	}

	word = ((word << 4) | ambiguity2basebit(seq[i])) & keymask;
	/*
	 * printf("seq[%03d]=%c => word=%04x & %04x == %04x, %04x == %04x %c"
	 *       " score %d\n",
	 *       i, seq[i], word, key, word & key,
	 *       ~key, (signed)(word & ~key),
	 *       " M"[(word & key) && !(word & ~key)],
	 *       score);
	 */

	/*
	 * For exact matches with no ambiguity codes this can just be
	 * (word & key) == key.
	 */
	if ((word & key) && !(word & ~key)) {
	    if (score < 0) {
		pads = maxscore = score = 0;
		start = i-(keylen-1);
	    }
	    score += keyjump;
	    if (score >= maxscore) {
		maxscore = score;
		end = i;
	    }
	    for (j = 0; j < keyjump-1;) {
		if (seq[++i] == '*') {
		    /* printf("SEQ[%03d]='*'\n", i); */
		    pads++;
		    continue;
		}
		j++;
		word = ((word << 4) | ambiguity2basebit(seq[i])) & keymask;

		/*
		 * printf("SEQ[%03d]=%c => word=%04x & %04x == %04x, "
		 *        "%04x == %04x\n",
		 *        i, seq[i], word, key, word & key,
		 *        ~key, (signed)(word & ~key));
		 */
	    }
	} else {
	    if (score >= 0) {
		if (--score < 0 || score <= maxscore - maxdrop) {
		    /* printf("start=%d end=%d pads=%d minsize=%d\n",
		       start, end, pads, minsize); */
		    if ((signed)(end - start + 1 -pads) >= minsize) {
			memset(&filt[start], filter_char, end-start+1);
		    }
		    score = -1;
		    pads = 0;
		    maxscore = 0;
		}
	    } else {
		pads = 0;
		score = -1;
	    }
	}
    }

    if (score >= 0 && end-start+1-pads >= minsize) {
	memset(&filt[start], filter_char, end-start+1);
    }
    
    return 0;
}

/*
 * Search for occurences of rep in seq (of length len).
 * minsize is the minimum length of region containing 'rep' and minscore
 * is effectively the portion of minsize that is 'rep'.
 * So minsize==minscore implies 100% match.
 *    minsize==minscore+2 implies 2 bases mismatch allowed.
 */
int filter_words_local(char *seq, char *filt, size_t len, char *rep,
		       int minsize, int minscore, char filter_char) {
    size_t i, j;
    fkey_t word = 0;
    int score = -1, maxscore = 0;
    size_t start = 0, end = 0;
    fkey_t key, keymask;
    int keyjump, keylen;
    int pads = 0;

    key = construct_key(rep, &keymask, &keylen, &keyjump);
    minscore *= 100; /* precision to 2 decimal points */

    /* Start with an entire word */
    for (i = j = 0; i < len && j < keylen-1; i++) {
	if (seq[i] == '*') {
	    pads++;
	    continue;
	}

	word = ((word << 4) | ambiguity2basebit(seq[i])) & keymask;
	j++;
    }

    /* Scan through looking for matches of a word */
    for (; i < len; i++) {
	if (seq[i] == '*'){
	    pads++;
	    continue;
	}

	word = ((word << 4) | ambiguity2basebit(seq[i])) & keymask;

	/*
	 * For exact matches with no ambiguity codes this can just be
	 * (word & key) == key.
	 */
	if ((word & key) && !(word & ~key)) {
	    if (score < 0) {
		pads = maxscore = score = 0;
		start = i-(keylen-1);
	    }
	    score += keyjump * MATCH*100;
	    if (score >= maxscore) {
		maxscore = score;
		end = i;
	    }
	    for (j = 0; j < keyjump-1;) {
		if (seq[++i] == '*') {
		    pads++;
		    continue;
		}
		j++;
		word = ((word << 4) | ambiguity2basebit(seq[i])) & keymask;
	    }
	} else {
	    if (score > 0)
		score += (MISMATCH*100.0);

	    if (score <= 0) {
		if (end-start+1 >= minsize && maxscore >= minscore) {
		    memset(&filt[start], filter_char, end-start+1);
		}
		maxscore = pads = 0;
		score = -1;
	    }
	}
    }

    if (end-start+1 >= minsize && maxscore >= minscore) {
	memset(&filt[start], filter_char, end-start+1);
    }

    return 0;
}
