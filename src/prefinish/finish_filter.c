/* Copyright Genome Research Limited (GRL). All rights reserved */

#include <stdio.h>
#include "dust.h"
#include "finish.h"
#include "finish_filter.h"
#include "filter_words.h"
#include "xalloc.h"
#include "dna_utils.h"

/*
 * Filter a set of known issues using my own 'words' algorithm.
 */
static void finish_filter_words(finish_t *fin, char *orig,
				char *seq, int clen) {
    if (fin->opts.debug[FIN_DEBUG_DUST])
	puts("Filtering using poly-* words...");

    filter_words(orig, seq, clen, "AAAA", 12, 4, FILTER_POLYA);
    filter_words(orig, seq, clen, "CCCC", 12, 4, FILTER_POLYC);
    filter_words(orig, seq, clen, "GGGG", 12, 4, FILTER_POLYG);
    filter_words(orig, seq, clen, "TTTT", 12, 4, FILTER_POLYT);
    filter_words(orig, seq, clen, "KKKK", 12, 4, FILTER_POLYK);
    filter_words(orig, seq, clen, "RRRR", 12, 4, FILTER_POLYR);
    filter_words(orig, seq, clen, "MMMM", 12, 4, FILTER_POLYM);
    filter_words(orig, seq, clen, "WWWW", 12, 4, FILTER_POLYW);
    filter_words(orig, seq, clen, "YYYY", 12, 4, FILTER_POLYY);
    filter_words(orig, seq, clen, "SSSS", 12, 4, FILTER_POLYS);
}

/*
 * Filter for general low complexity using dust
 */
static void finish_filter_dust(finish_t *fin, char *seq, int seqlen) {
    int i;

    if (fin->opts.debug[FIN_DEBUG_DUST])
	puts("Filtering using dust...");

    /* Filter using Dust */
    set_dust_level(fin->opts.dust_level);
    dust(seqlen, seq);

    /* Look for low complexity data within 32 bases of the end,
     * if so extend to ends
     */
    for (i = 0; i < seqlen && i < 32; i++) {
	if (seq[i] == '#') {
	    for (i = 0; i < 32 && i < seqlen; i++)
		seq[i] = '#';
	    break;
	}
    }

    for (i = 0; seqlen-1-i >= 0 && i < 32; i++) {
	if (seq[seqlen-1-i] == '#') {
	    for (i = 0; seqlen-1-i >= 0 && i < 32; i++)
		seq[seqlen-1-i] = '#';
	    break;
	}
    }
    /*
   if (fin->opts.debug[FIN_DEBUG_DUST] > 1)
	printf("%.*s\n", seqlen, seq);
    */
}

#undef is_filtered
int is_filtered(char c) {
    return c == '#' || (c >=  FILTER_POLYA && c <= FILTER_POLYY);
}


/*
 * Filters 'seq' using dust and the 'words' algorithm here.
 * If 'seq' is NULL then it filters fin->cons putting the result in
 * fin->filtered.
 */
void finish_filter(finish_t *fin, char *seq, int len) {
    char *copy;

    /* Create fin->filtered if desired */
    if (!seq) {
	int clen = io_clength(fin->io, fin->contig);

	if (NULL == (fin->filtered = (char *)xmalloc(clen)))
	    return;

	/* Filter using Dust */
	memcpy(fin->filtered, fin->cons, clen);
	seq = fin->filtered;
	len = clen;
    }

    /*
     * Multiple rounds of filtering, but each needs to do a lookup against
     * the unfiltered sequence, so we filter from a temp. copy.
     */
    copy = (char *)malloc(len);
    memcpy(copy, seq, len);

    finish_filter_dust(fin, seq, len);
    finish_filter_words(fin, copy, seq, len);

    xfree(copy);

    if (fin->opts.debug[FIN_DEBUG_DUST] > 1)
	printf("filtered %.*s\n", len, seq);
}
