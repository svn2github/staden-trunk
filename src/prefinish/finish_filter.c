/* Copyright Genome Research Limited (GRL). All rights reserved */

#include <stdio.h>
#include "dust.h"
#include "finish.h"
#include "finish_filter.h"
#include "xalloc.h"
#include "dna_utils.h"

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

/*
 * Construct a key and a mask (allbits = 1) from a sequence.
 * If sequence is longer than Number_of_bits(key_t)/4 then the start of the
 * key will be clipped.
 *
 * Seq can contain the standard ambiguity codes.
 */
static key_t construct_key(char *seq, key_t *mask, int *keylen, int *keyjmp) {
    key_t k = 0, m = 0;
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
static int filter_words(char *seq, char *filt, size_t len, char *rep,
			int minsize, int maxdrop, char filter_char) {
    size_t i;
    key_t word = 0;
    int score = -1, maxscore;
    size_t start = 0, end = 0;
    int last_match_end = 0;
    key_t key, keymask;
    int keyjump, keylen;
    int pads = 0;

    key = construct_key(rep, &keymask, &keylen, &keyjump);

    /* Scan through looking for matches of a word */
    for (i = 0; i < len; i++) {
	int real_match;

	if (seq[i] == '*'){
	    pads++;
	    continue;
	}

	word = ((word << 4) | ambiguity2basebit(seq[i])) & keymask;
	real_match = 0;
	if (word & key) {
	    key_t word2 = word & key, km = keymask;

	    /* May not be a real match if not all 4-bit values also match */
	    real_match=1;
	    while (km) {
		if (!(word2 & 15)) {
		    real_match = 0;
		    break;
		}
		word2 >>= 4;
		km >>= 4;
	    }
	}
	if (real_match) {
	    if (score < 0) {
		maxscore = score = 0;
		start = i-(keylen-1);
	    }
	    /* score += keyjump; */
	    score += keylen;
	    if (score >= maxscore) {
		maxscore = score;
		end = i;
	    }
	} else {
	    if (score >= 0) {
		if (--score < 0 || score <= maxscore - maxdrop) {
		    if (end - start + 1 -pads >= minsize) {
			memset(&filt[start], filter_char, end-start+1);

			last_match_end = end;
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
 * Filter a set of known issues using my own 'words' algorithm.
 */
static void finish_filter_words(finish_t *fin) {
    int clen = io_clength(fin->io, fin->contig);

    if (fin->opts.debug[FIN_DEBUG_DUST])
	puts("Filtering using poly-* words...");

    filter_words(fin->cons, fin->filtered, clen, "AAAA", 12, 4, FILTER_POLYA);
    filter_words(fin->cons, fin->filtered, clen, "CCCC", 12, 4, FILTER_POLYC);
    filter_words(fin->cons, fin->filtered, clen, "GGGG", 12, 4, FILTER_POLYG);
    filter_words(fin->cons, fin->filtered, clen, "TTTT", 12, 4, FILTER_POLYT);
    filter_words(fin->cons, fin->filtered, clen, "KKKK", 12, 4, FILTER_POLYK);
    filter_words(fin->cons, fin->filtered, clen, "RRRR", 12, 4, FILTER_POLYR);
    filter_words(fin->cons, fin->filtered, clen, "MMMM", 12, 4, FILTER_POLYM);
    filter_words(fin->cons, fin->filtered, clen, "WWWW", 12, 4, FILTER_POLYW);
    filter_words(fin->cons, fin->filtered, clen, "YYYY", 12, 4, FILTER_POLYY);
    filter_words(fin->cons, fin->filtered, clen, "SSSS", 12, 4, FILTER_POLYS);
}

/*
 * Filter for general low complexity using dust
 */
static void finish_filter_dust(finish_t *fin) {
    int i, clen;

    if (fin->opts.debug[FIN_DEBUG_DUST])
	puts("Filtering using dust...");

    clen = io_clength(fin->io, fin->contig);
    if (NULL == (fin->filtered = (char *)xmalloc(clen)))
	return;

    /* Filter using Dust */
    memcpy(fin->filtered, fin->cons, clen);
    set_dust_level(fin->opts.dust_level);
    dust(clen, fin->filtered);

    /* Look for low complexity data with 32 of the end, if so extend to ends */
    for (i = 0; i < clen && i < 32; i++) {
	if (fin->filtered[i] == '#') {
	    for (i = 0; i < 32 && i < clen; i++)
		fin->filtered[i] = '#';
	    break;
	}
    }

    for (i = 0; clen-1-i >= 0 && i < 32; i++) {
	if (fin->filtered[clen-1-i] == '#') {
	    for (i = 0; clen-1-i >= 0 && i < 32; i++)
		fin->filtered[clen-1-i] = '#';
	    break;
	}
    }
    /*
   if (fin->opts.debug[FIN_DEBUG_DUST] > 1)
	printf("%.*s\n", clen, fin->filtered);
    */
}

#undef is_filtered
int is_filtered(char c) {
    return c == '#' || (c >=  FILTER_POLYA && c <= FILTER_POLYY);
}

void finish_filter(finish_t *fin) {
    finish_filter_dust(fin);
    finish_filter_words(fin);

    if (fin->opts.debug[FIN_DEBUG_DUST] > 1)
	printf("%.*s\n", io_clength(fin->io, fin->contig), fin->filtered);
}
