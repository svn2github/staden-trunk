/*
 * Other ideas:
 *
 * 1. Identify all readings on inconsistent templates and strike them out of
 *    the contig. If holes remain, then verify assembly at that point.
 *
 * 2. Look for regions of high depth. If high depth and avg length of insert
 *    spanning the gap is shorter than normal, then probable collapsed repeat.
 *
 * 3. Use the diploid scanner and haplotype splitting code too as part of
 *    the auto-break mechanism.
 *
 * 4. Spot repeats and divergence => end of repeats. Check all divergence
 *    points to ensure consistency. (Need a decent local alignment algorithm)
 *
 * 5. Emulate neighbourhoods by using a word that is (eg) 6mer, 1mer gap, 6mer.
 *    Hence with all overlapping coordinates this allows 1 base to differ
 *    while still using the memory for a 12mer and not 13mer.
 *
 * 7. Ressurrect stop-detector. Clip after stop - it's probably wrong data.
 */

/*

  Given word size N, the number of combinations of R of N being GC
  is 2^N * N! / ((N-R)! R!)
  
  Ie for N=8:
  
  0/8 256
  1/8 2048
  2/8 7168
  3/8 14336
  4/8 17920
  5/8 14336
  6/8 7168
  7/8 2048
  8/8 256

  Given a GC content of x (0 <= x <= 1) we therefore expect
  (x/2)^R * ((1-x)/2)^(N-R) * 2^N * N! / ((N-R)!R!) matches.
*/


#include "IO.h"
#include "misc.h"
#include "xalloc.h"
#include "text_output.h"
#include "filter_words.h"
#include "template.h"
#include "auto_break.h"
#include "gap_globals.h"
#include "qual.h"
#include "qualIO.h"
#include "dstring.h"

#define MIN3(a,b,c) (MIN(MIN((a),(b)),(c)))
#define MAX3(a,b,c) (MAX(MAX((a),(b)),(c)))

#ifndef WS
#    define WS 12
#endif
#define WS2 ((int)(WS/2))

#define ALLB(ws) ((1<<(2*(ws)))-1)

#define MIN_OVERLAP 10
#define MAXTSIZE 10000
#define MIN_HARD_SCORE -9
#define MIN_SOFT_SCORE 0
#define MIN_SOFT_VALID 4
/* #define NORMALISE_FOR_GC 1 */
#define NORMALISE_STR_SCORES

#define CONTIG_END_IGNORE 200

/* #define GAPPED_WORDS */

#ifdef GAPPED_WORDS
#  define WORD2GAPPED(w) ( ((w) & ALLB(WS2)) | (((w) >> 2) & ~ALLB(WS2)) )
#else
#  define WORD2GAPPED(w) (w)
#endif

typedef unsigned short count_t;
#define MAX_COUNTS ((1<<16)-1)

/* #define DEBUG */
/* #define DEBUG_SEQ */

typedef struct contig_region {
    int start;
    int end;
    int deleted;
    int rnum; /* reading at left end */
    int valid;
} contig_region_t;


typedef struct clip_pos {
    int left;
    int right;
} clip_pos_t;

void dump_gaps(Array gaps) {
    int i;
    puts("\n");
    for (i = 0; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);
	printf("Gap %d\t%d %d %d %d\n", i, 
	       gap->start, gap->end, gap->rnum, gap->deleted);
    }
}

static count_t counts[1<<(2*WS)]; /* 4^WS */
static int lookup[256];
static int clookup[256];
static double probs[4];

#if 0
/*
 * Given a template_c struct that is potentially inconsistent, this works
 * out likely real extents for the template using the max template size
 * as a guess.
 *
 * Eg.               <------ ...... ----->
 * may have extents:
 *       |<-------------------------------------------->|
 *
 * This allows us to consider this template as having an "invalidating"
 * influence on any joins over that entire extent. Ranges returned may be
 * outside the boundaries of the contig (and hence possibly negative too).
 *
 * Returns 0 for success and values in start/end
 *        -1 for failure
 */
static int template_extents(GapIO *io, template_c *temp,
			     int *start, int *end) {
    GTemplates te;
    item_t *item;
    int st = MIN(temp->start, temp->end);
    int en = MAX(temp->start, temp->end);

    if (!temp || !start || !end || temp->num <= 0)
	return -1;

    template_read(io, temp->num, te);

    for (item = head(temp->gel_cont); item; item = item->next) {
	GReadings r;
	gel_cont_t *gc = (gel_cont_t *)(item->data);
	int ptype;
	int left;
	int right;

	gel_read(io, gc->read, r);
	left  = r.position;
	right = r.position + r.sequence_length-1;;

	ptype = PRIMER_TYPE_GUESS(r);

	if ((temp->oflags & TEMP_OFLAG_IGNORE_PTYPE) ||
	    (ptype == GAP_PRIMER_FORWARD ||
	     ptype == GAP_PRIMER_CUSTFOR ||
	     ((temp->oflags & TEMP_OFLAG_IGNORE_PTYPE34) &&
	      (ptype == GAP_PRIMER_CUSTFOR ||
	       ptype == GAP_PRIMER_CUSTREV)))) {
	    if (r.sense == GAP_SENSE_ORIGINAL)
		right += te.insert_length_max;
	    else
		left -= te.insert_length_max;
	}

	if ((temp->oflags & TEMP_OFLAG_IGNORE_PTYPE) ||
	    (ptype == GAP_PRIMER_REVERSE ||
	     ptype == GAP_PRIMER_CUSTREV ||
	     ((temp->oflags & TEMP_OFLAG_IGNORE_PTYPE34) &&
	      (ptype == GAP_PRIMER_CUSTFOR ||
	       ptype == GAP_PRIMER_CUSTREV)))) {
	    if (r.sense == GAP_SENSE_ORIGINAL)
		right += te.insert_length_max;
	    else
		left -= te.insert_length_max;
	}

	if (st > left)
	    st = left;
	if (en < right)
	    en = right;
    }

    *start = st;
    *end = en;

    return 0;
}
#endif

/*
 * This checks if an inconsistent template appears to be capable of
 * spanning 'gap'.
 *
 * As there are so many reasons why a template can be inconsistent this
 * iterates through all readings on the template in 'contig' to see if the
 * range of possibilities for the "other end" could be on the opposite side
 * to the gap.
 *
 * eg where F is invalid and may span if repositioned.
 * <---F---.........---R--->   |gap|
 *  ---F--->........---R--->   |gap|
 *
 * eg where no span is possible unless both readings are misassembled.
 * <---F---........<---R---    |gap|
 *
 * Returns 1 for potentially spanning
 *         0 for not spanning
 */
static int bad_template_span(GapIO *io, int contig, template_c *t,
			     contig_region_t *gap) {
    GTemplates te;
    item_t *item;

    template_read(io, t->num, te);

    for (item = head(t->gel_cont); item; item = item->next) {
	GReadings r;
	gel_cont_t *gc = (gel_cont_t *)(item->data);

	if (gc->contig != contig)
	    continue;

	gel_read(io, gc->read, r);

	/* ---> |gap| */
	if (r.position < gap->start &&
	    r.sense == GAP_SENSE_ORIGINAL &&
	    r.position + te.insert_length_max > gap->end) {
	    return 1;
	}

	/* |gap| <--- */
	if (r.position + r.sequence_length - 1 > gap->end &&
	    r.sense != GAP_SENSE_ORIGINAL &&
	    r.position + r.sequence_length - 1 - te.insert_length_max <
	        gap->start) {
	    return 1;
	}
    }

    return 0;
}

/*
 * Checks the region from 'start' to 'end' to see whether it can be
 * validated by nearby read-pairs.
 *
 * It does this by looking up to MAXTSIZE either size of start/end
 * identifying readings that point inwards to the hole. For these reads it
 * then verifies if they validate, have no impact, or reject the region.
 *
 * Returns a score, +ve for confirmed, -ve for denied and 0 for unknown.
 */
static int check_read_pairs(GapIO *io, int contig, template_c **tarr,
			    contig_region_t *gap) {
    int i, rnum;
    int valid = 0, invalid = 0;

#ifdef DEBUG    
    printf("\n***** Contig %d, gap %d..%d, near #%d *****\n",
	   contig, gap->start, gap->end, gap->rnum);
#endif


    /* Chain left a while until we're MAXTSIZE away from the problem */
    for (rnum = gap->rnum; rnum; rnum = io_lnbr(io, rnum)) {
	if (io_relpos(io, rnum) < gap->start - MAXTSIZE)
	    break;
     }
    if (rnum == 0)
	rnum = io_clnbr(io, contig);

    /* Identify templates that appears to map close to this problem. */
    for (; rnum; rnum = io_rnbr(io, rnum)) {
	GReadings r;
	template_c *t;

	if (io_relpos(io, rnum) > gap->end + MAXTSIZE)
	    break;
	gel_read(io, rnum, r);

	/* Skip if we've validated this template already */
	t = tarr[r.template];
	if (t->num < 0)
	    continue;

#ifdef DEBUG    
	printf("=== Validating with template %d(%s)\n", r.template,
	       get_template_name(io, r.template));
#endif
	get_template_positions(io, t, contig);
	check_template_c(io, t);
	/* if (t->flags & TEMP_FLAG_SPANNING) */
#ifdef DEBUG    
	dump_template(t);
#endif
	
	/*
	 * Not all inconsistent templates have the potential to span the
	 * gap, so we check.
	 */
	if ((t->consistency &~  TEMP_CONSIST_UNKNOWN) != 0) {
	    if (bad_template_span(io, contig, t, gap)) {
#ifdef DEBUG    
		printf("Invalid %s\n", get_template_name(io, r.template));
#endif
		invalid++;
	    }

	    /* otherwise ignore for good/bad count */
	    t->num = -t->num;
	    continue;
	}
	
	t->num = -t->num;

	/*
	 * Otherwise if it spans the region, then we consider it validation.
	 * We don't just look at start/end coordinates in 't' as we can get
	 * the case where judicious quality clipping has meant that the vector
	 * is identified both sides, but the "good" quality clipped data
	 * itself is all one one side. (Phrap can sometimes do this with its
	 * methods to quality clip based on consensus similarity and hence
	 * hides the unique data away to only leave the falsely joined
	 * repetitive bit.)
	 */
	if (template_covered(io, t, contig,
			     MIN3(gap->start - 100, t->start, t->end),
			     gap->start) > 10 &&
	    template_covered(io, t, contig,
			     gap->end,
			     MAX3(gap->end + 100, t->start, t->end)) > 10) {
#ifdef DEBUG    
	    printf("Valid %s\n", get_template_name(io, r.template));
#endif
	    valid++;
	    continue;
	}
    }

    /* Re-negate the visited template numbers */
    for (i = 1; i <= Ntemplates(io); i++) {
	if (tarr[i] && tarr[i]->num < 0)
	    tarr[i]->num = -tarr[i]->num;
    }

    printf("Found %d valid, %d invalid\n", valid, invalid);

    if (valid - invalid < MIN_HARD_SCORE ||
	(valid - invalid < MIN_SOFT_SCORE &&
	 valid < MIN_SOFT_VALID))
	gap->valid = 0;
    else
	gap->valid = 1;

    return valid - invalid;
}

#undef DEBUG

/*
 * Turns an array of coverage data (0 => no coverage) to an Array of
 * contig_region_t structs.
 */
static Array coverage2contig_regions(char *valid, int clen) {
    int i, j, ngaps = 0;
    Array gaps;

    gaps = ArrayCreate(sizeof(contig_region_t), 0);
    if (NULL == gaps)
	return NULL;

#if 0
    /* Now report invalid regions, skipping the very start and end */
    for (; clen >= 1; clen--) {
	if (valid[clen])
	    break;
    }
    for (i = 1; i <= clen; i++) {
	if (valid[i])
	    break;
    }
#else
    /* Handle start and end gaps separately */
    for (i = clen; i >= 1; i--) {
	if (valid[i])
	    break;
    }
    if (clen != i) {
	contig_region_t *r;
	r = (contig_region_t *)ArrayRef(gaps, ngaps);
	r->start = i;
	r->end = clen;
	r->deleted = 0;
	r->valid = 1;
	ngaps++;
    }

    for (i = 1; i <= clen; i++) {
	if (valid[i])
	    break;
    }
    if (i != 1) {
	contig_region_t *r;
	r = (contig_region_t *)ArrayRef(gaps, ngaps);
	r->start = 1;
	r->end = i;
	r->deleted = 0;
	r->valid = 1;
	ngaps++;
    }
#endif

    for (; i <= clen; i++) {
	if (!valid[i]) {
	    contig_region_t *r;
	    for (j = i+1; j <= clen; j++) {
		if (valid[j])
		    break;
	    }
	    r = (contig_region_t *)ArrayRef(gaps, ngaps);
	    r->start = i;
	    r->end = j-1;
	    r->deleted = 0;
	    r->valid = 1;
	    ngaps++;
	    i = j-1;
	}
    }
    
    return gaps;
}

/*
 * Counts all WS-mers in the gap4 database to detect overly frequent words.
 * Only uses clipped data (to avoid vector etc).
 */
static void init_tables(void) {
    int i;

    /* Initialise lookup tables */
    for (i = 0; i < 256; i++) {
	lookup[i] = -1;
	clookup[i] = -1;
    }
    lookup['A']  = lookup['a'] = 0;
    lookup['C']  = lookup['c'] = 1;
    lookup['G']  = lookup['g'] = 2;
    lookup['T']  = lookup['t'] = 3;
#ifdef GAPPED_WORDS
    clookup['A'] = clookup['a'] = 3 << (2*WS);
    clookup['C'] = clookup['c'] = 2 << (2*WS);
    clookup['G'] = clookup['g'] = 1 << (2*WS);
    clookup['T'] = clookup['t'] = 0 << (2*WS);
#else
    clookup['A'] = clookup['a'] = 3 << (2*WS-2);
    clookup['C'] = clookup['c'] = 2 << (2*WS-2);
    clookup['G'] = clookup['g'] = 1 << (2*WS-2);
    clookup['T'] = clookup['t'] = 0 << (2*WS-2);
#endif

    memset(counts, 0, (1<<(2*WS)) * sizeof(*counts));
}

#ifdef NORMALISE_FOR_GC
static void init_gc_table(double gc) {
    probs[1] = probs[2] = gc/2;
    probs[0] = probs[3] = (1-gc)/2;
}
#endif


static char *word2str(int word) {
    static char str[WS+2];
    signed int i, j;
    
    for (j = WS, i = 0; i < WS; i++) {
#ifdef GAPPED_WORDS
	if (i == (int)(WS/2))
	    str[j--] = '.';
#endif
	str[j--] = "ACGT"[word & 3];
	word >>= 2;
    }
    str[WS+1] = 0;
    
#ifdef GAPPED_WORDS
    return str;
#else
    return str+1;
#endif
}

#if 0
static char *word2str2(int word) {
    static char str[WS+2];
    int i, j;

    for (j = WS, i = 0; i < WS+1; i++) {
	str[j--] = "ACGT"[word & 3];
	word >>= 2;
    }
    str[WS+1] = 0;

    return str;
}
#endif

double compute_prob(int word) {
    int i;
    double prob = 1;
    
    for (i = 0; i < WS; i++) {
	int base = (word >> (2*(WS-1)-2*i)) & 3;
	prob *= probs[base];
    }

    return prob;
}

/*
 * Computes the maximum theoretical redundancy of this word. Ie if we have
 * a GT repeat then word GTGTGTGTGTGT occurs 6 times more frequently within
 * a run than we'd expect.
 *
 * It does this by sliding word w along looking at the overlap. Ie:
 * ABCDEABC
 *  ABCDEABC          no match in overlap as BCDEABC != ABCDEAB
 *
 * ABCDEABC
 *      ABCDEABC      matches as ABC == ABC
 *
 * Note that on average even a random word will have 1/4 chance of the last
 * base matching the first base so claim to have a redundancy of WS/(WS-1).
 * Similarly 1/16 times it'll have redundancy of WS/(WS-2). We could
 * compensate for this by factoring this expectation into the equation, maybe
 * even utilising the GC content in there too, but for now I think it skews
 * things sufficiently little to not be too concerned.
 *
 * Returns the total number of new counts.
 */
int normalise_str_scores(void) {
    int w;
    int tc = 0;
    
    for (w = 0; w < 1<<(2*WS); w++) {
	int i, m = (1 << (2*(WS-1))) - 1;

	if (!counts[w])
	    continue;

	for (i = 1; i <= WS; i++, m >>= 2) {
	    if ( (w >> (2*i)) == (w & m) )
		break;
	}
	counts[w] /= (double)WS/i;
	if (counts[w] == 0)
	    counts[w] = 1; /* Min count of 1 */

	tc += counts[w];
    }

    return tc;
}

/*
 * Fills out the counts[] array based on all words of length WS in the
 * individual sequences. 
 *
 * Returns the total number of words indexed.
 */
int word_count(GapIO *io, double *gcp, int *depthp) {
    int rnum, i, j, tw = 0, gc = 0, at = 0;
    size_t total_contig_len = 0, total_read_len = 0;

    init_tables();

    for (i = 1; i <= NumContigs(io); i++) {
	total_contig_len += io_clength(io, i);
    }

    for (rnum = 1; rnum <= NumReadings(io); rnum++) {
	GReadings r;
	char *s, *seq;
	unsigned int word, cword;

	/* printf("Rnum #%d\n", rnum); */

	gel_read(io, rnum, r);
	if (NULL == (seq = TextAllocRead(io, r.sequence))) {
	    continue;
	}

	total_read_len += r.sequence_length;

	s = seq;
	s[r.end-1] = 0;
	s += r.start;

	cword = word = 0;

	/*
	printf("! %s ", word2str(word));
	printf("%s\n", word2str(cword));
	*/

	for (j = 0; *s; s++) {
	    if (*s == '*')
		continue;

	    switch (lookup[*s]) {
	    case -1:
		j = 0;
		break;
	    case 0:
	    case 3:
		at++;
		j++;
		word <<= 2;
		word |= lookup[*s];
		cword >>= 2;
		cword |= clookup[*s];
		break;
	    case 1:
	    case 2:
		gc++;
		j++;
		word <<= 2;
		word |= lookup[*s];
		cword >>= 2;
		cword |= clookup[*s];
		break;
	    }
	    
	    /*
	    printf("%c %s ", j >= WS ? *s : '.', word2str(word));
	    printf("%s %08x %08x\n", word2str(cword), cword, clookup[*s]);
	    */

	    if (j >= WS) {
		/* Remove middle nucleotide from word before storing */
		unsigned int w2, cw2;

		w2  = WORD2GAPPED(word);
		cw2 = WORD2GAPPED(cword);

		/*
		printf("%c %s (%d)\n",
		       *s, word2str(w2), counts[w2 & ALLB(WS)]);
		*/

		/* Count both original and complementary word */
		if (counts[ w2 & ALLB(WS)] < MAX_COUNTS)
		    counts[ w2 & ALLB(WS)]++;
		if (counts[cw2 & ALLB(WS)] < MAX_COUNTS)
		    counts[cw2 & ALLB(WS)]++;

		tw += 2;
	    }
	}

	xfree(seq);
    }

    printf("Total words = %d, GC = %5.2f%%, depth = %5.2f\n",
	   tw, (100.0*gc)/(gc+at),
	   (double)total_read_len / total_contig_len);

#ifdef NORMALISE_FOR_GC
    /* normalise counts by GC content */
    init_gc_table((double)gc/(gc+at));

    for (i = 0; i < (1<<(2*WS)); i++) {
	int c = counts[i];
	c /= (1<<(2*WS)) * compute_prob(i);
	if (c > MAX_COUNTS)
	    c = MAX_COUNTS;
	counts[i] = c;
    }
#endif

#ifdef NORMALISE_STR_SCORES
    tw = normalise_str_scores();
#endif

    if (gcp)
	*gcp = (double)gc/(gc+at);

    if (depthp)
	*depthp = (double)total_read_len / total_contig_len;

    return tw;
}

/*
 * Fills out the counts[] array based on all words of length WS in the
 * consensus sequences. 
 *
 * Returns the total number of words indexed.
 */
int word_count_cons(GapIO *io, double *gcp) {
    int cnum, j, tw = 0, gc = 0, at = 0;

    init_tables();

    for (cnum = 1; cnum <= NumContigs(io); cnum++) {
	int clen = io_clength(io, cnum);
	char *cons = (char *)malloc(clen);
	char *s;
	unsigned int word, cword;

	/* Skip single-read contigs */
	if (io_clnbr(io, cnum) == io_crnbr(io, cnum)) {
	    printf("Skipping contig =%d; singleton\n", cnum);
	    continue;
	}

	if (clen < 1630) {
	    printf("Skipping contig =%d; len %d < %d\n", cnum, clen, 1000);
	    continue;
	}

	calc_consensus(cnum, 1, clen, CON_SUM, cons, NULL, NULL, NULL,
		       consensus_cutoff, quality_cutoff,
		       database_info, (void *)io);

	/* printf("CONS =%d/%d: %.*s\n", cnum, NumContigs(io), clen, cons); */

	if (clen <= 2*CONTIG_END_IGNORE)
	    continue;

	s = cons + CONTIG_END_IGNORE;
	cons[clen-1-CONTIG_END_IGNORE] = 0;

	cword = word = 0;
	for (j = 0; *s; s++) {
	    if (*s == '*')
		continue;

	    switch (lookup[*s]) {
	    case -1:
		j = 0;
		break;
	    case 0:
	    case 3:
		at++;
		j++;
		word <<= 2;
		word |= lookup[*s];
		cword >>= 2;
		cword |= clookup[*s];
		break;
	    case 1:
	    case 2:
		gc++;
		j++;
		word <<= 2;
		word |= lookup[*s];
		cword >>= 2;
		cword |= clookup[*s];
		break;
	    }
	    
	    /*
	    printf("%c %s ", j >= WS ? *s : '.', word2str(word));
	    printf("%s %08x %08x\n", word2str(cword), cword, clookup[*s]);
	    */

	    if (j >= WS) {
		/* Remove middle nucleotide from word before storing */
		unsigned int w2, cw2;

		w2  = WORD2GAPPED(word);
		cw2 = WORD2GAPPED(cword);
		
		/* Count both original and complementary word */
		if (counts[ w2 & ALLB(WS)] < MAX_COUNTS)
		    counts[ w2 & ALLB(WS)]++;
		if (counts[cw2 & ALLB(WS)] < MAX_COUNTS)
		    counts[cw2 & ALLB(WS)]++;

		/*
		printf("%c %s (%d)",
		       *s, word2str(w2), counts[w2 & ALLB(WS)]);
		printf("\t%s (%d)\n",
		       word2str(cw2), counts[cw2 & ALLB(WS)]);
		*/
		tw += 2;
	    }
	}

	xfree(cons);
    }

    printf("Total words = %d, GC = %5.2f%%\n", tw, (100.0*gc)/(gc+at));

#ifdef NORMALISE_FOR_GC
    /* normalise counts by GC content */
    init_gc_table((double)gc/(gc+at));

    for (j = 0; j < (1<<(2*WS)); j++) {
	int c = counts[j];
	c /= (1<<(2*WS)) * compute_prob(j);
	if (c > MAX_COUNTS)
	    c = MAX_COUNTS;
	counts[j] = c;
    }
#endif

#ifdef NORMALISE_STR_SCORES
    tw = normalise_str_scores();
#endif

    if (gcp)
	*gcp = (double)gc/(gc+at);

    return tw;
}

void print_counts(double min) {
    int i;
    for (i = 0; i < 1<<(2*WS); i++) {
	if (counts[i] >= min)
	    printf("%s %d\n", word2str(i), counts[i]);
    }
}

void print_bins(void) {
    int i, j;
    int bins[10000];
    
    memset(bins, 0, sizeof(*bins)*10000);
    for (i = 0; i < (1<<(2*WS)); i++) {
	if (counts[i] < 10000)
	    bins[counts[i]]++;
    }

    for (i = 0; i < 10000; i++)
	if (bins[i])
	    break;
    for (j = 9999; j >= 0; j--)
	if (bins[j])
	    break;
    for (; i <= j; i++) {
	printf("%d %d\n", i, bins[i]);
    }
}

int filter_common_words(char *seq, char *filt, size_t len, int tw,
			double depth, double score, char filter_char,
			int debug) {
    size_t i, j;
    unsigned int word = 0;
    int pads = 0;
    double rescale = 1;

    memcpy(filt, seq, len);

    /* Start with an entire word */
    for (i = j = 0; i < WS && i < len; i++) {
	if (seq[i] == '*') {
	    pads++;
	    continue;
	}

	word = (word << 2) | lookup[seq[i]];
	j++;
    }
    
    /*
    printf("compare %f * %f (=%f) vs %f * %d / %d (=%f)\n",
	   score, depth, score * depth,
	   score, tw, 1<<(2*WS),
	   score * (double)tw / (1<<(2*WS)));
    */

    /*
     * A quick and easy hack, but not particularly robust:
     * If we've observed more words than there are buckets then we
     * expect the number of hits to be higher than the depth alone would
     * suggest.
     * Therefore we downscale a bit to compensate.
     */
    if (tw /  (1<<(2*WS)) > 1) {
	rescale = ((double)tw / (1<<(2*WS))) / depth;
    }

    /* Scan through looking for matches of a word */
    for (; i < len; i++) {
	if (seq[i] == '*'){
	    pads++;
	    continue;
	}

	word = (word << 2) | lookup[seq[i]];
	
	if (debug)
	    printf("Seq pos %ld %.*s: => %d",
		   (long)i, WS, seq+i, counts[word & ((1<<(2*WS))-1)]);
	/* if (counts[word & ((1<<(2*WS))-1)] >= score * (double)tw / (1<<(2*WS))) { */
	if ((counts[WORD2GAPPED(word) & ALLB(WS)]) / rescale
	    >= score * depth) {
	    /* FIXME: ignores pads for now */
#ifdef GAPPED_WORDS
	    memset(&filt[i-WS], filter_char, WS+1);
#else
	    memset(&filt[i-(WS-1)], filter_char, WS);
#endif
	    if (debug) putchar('*');
	}

	if (debug) putchar('\n');
    }

    /* Merge filtered blocks together if they have <= 4 bases gap */
    /* FIXME: merge if lots of valid read-pairs span gaps.
     * Eg this means that ....gap...gap.... becomes one gap if we find
     * linkage info claiming that the sequence gap...gap looks as if it
     * forms a valid contig itself.
     */
    for (i = 1; i < len; i++) {
	if (filt[i-1] == filter_char && filt[i] != filter_char) {
	    int g = i;
	    while (i < len && filt[i] != filter_char)
		i++;
	    if (i-g <= 4) {
		while (g != i && g < len)
		    filt[g++] = filter_char;
	    }
	}
    }

    return 0;
}

/*
 * Identifies regions where the input data mismatches the consensus and
 * mask out WS-1 characters either side of this location.
 * The logic is that these are base calling errors and as such will bias
 * the word-counting (as errors in repeats will make it appear to be a
 * unique word).
 */
void filter_consen_diffs(char *in, char *out, int len, char *cons) {
    int i, j;
    for (i = 0; i < len; i++) {
	if (in[i] == cons[i])
	    continue;
	for (j = i-(WS-1) >= 0 ? i-(WS-1) : 0;
	     j <= i+(WS-1) && j < len;
	     j++) {
	    out[j] = '%';
	}
    }
}

/*
 * Loops through all readings in a contig looking to see that specific low
 * complexity regions have readings that span either end. If not then we
 * label this region as a suspicious join.
 *
 * Returns: an Array of type contig_region_t on success
 *          NULL on failure
 */
static Array suspect_joins(GapIO *io, int contig, int tw, double filter_score,
			   double depth, clip_pos_t *clips) {
    int rnum, i;
    char *valid, *cp;
    int clen = io_clength(io, contig);
    char legal_chars[256];
    Array gaps;
    char *cons;
#if 0
    char *patterns[] = {"A",  "C",  "G",  "T",
			"AC", "AG", "AT", "CG", "CT", "GT"};
#endif

    if (NULL == (valid = (char *)xcalloc(clen+1, 1)))
	return NULL;

    memset(legal_chars, 0, 256);
    for (cp = "ACGTacgt"; *cp; cp++)
	legal_chars[*cp] = 1;

    /* Compute consensus */
    if (NULL == (cons = (char *)xmalloc(clen+1))) {
	xfree(valid);
	return NULL;
    }
    calc_consensus(contig, 1, clen, CON_SUM, cons, NULL, NULL, NULL,
		   consensus_cutoff, quality_cutoff,
		   database_info, (void *)io);
    
#if 0
    /*
     * Alternatively treat each motif separately and clip based on the
     * furthest left/right after all the clipping runs.
     * In realitly for speed the loop would be the other way around, but
     * conceptually it's sort of like this:
     */
    for (pnum = 0; pnum < sizeof(patterns)/sizeof(*patterns); pnum++) {
	printf("Pattern %d = '%s'\n", pnum, patterns[pnum]);
	for (rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum)) {
	    /* ... */
	    filter_words_local(seq, fseq, len, patterns[pnum], 12, 8, '#');
	    /* ... clip ... */
	}
    }
#endif

    /*
     * Loop through readings computing new clip points after masking.
     * Mark the remainder as 'valid' to identify the invalid bits.
     */
    for (rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum)) {
	char *seq, *fseq;
	size_t len;
	int unique, first, last;
	GReadings r;
	gel_read(io, rnum, r);

	if (!r.trace_name) {
	    continue;
	}

	if (NULL == (seq = TextAllocRead(io, r.sequence))) {
	    verror(ERR_WARN, "suspect_joins", "couldn't read sequence");
	    continue;
	}

	if (NULL == (fseq = strdup(seq))) {
	    xfree(seq);
	    continue;
	}
	len = strlen(seq);

	/* Filter over represented words */
	filter_common_words(seq, fseq, len, tw, depth, filter_score, '#', 0);

	/* Filter where the sequence is low qual and disagrees with consen */
	filter_consen_diffs(&seq[r.start], &fseq[r.start],
			    r.sequence_length, &cons[r.position-1]);
	
	/* Filter specific low-complexity regions */
	filter_words_local(seq, fseq, len, "A",   12, 8, '#');
	filter_words_local(seq, fseq, len, "C",   12, 8, '#');
	filter_words_local(seq, fseq, len, "G",   12, 8, '#');
	filter_words_local(seq, fseq, len, "T",   12, 8, '#');
	filter_words_local(seq, fseq, len, "AC",  12, 8, '#');
	filter_words_local(seq, fseq, len, "AG",  12, 8, '#');
	filter_words_local(seq, fseq, len, "AT",  12, 8, '#');
	filter_words_local(seq, fseq, len, "CG",  12, 8, '#');
	filter_words_local(seq, fseq, len, "CT",  12, 8, '#');
	filter_words_local(seq, fseq, len, "GT",  12, 8, '#');

#ifdef DEBUG_SEQ
	printf("SEQ #%d %d+%d %.*s\n",
	       rnum, r.position, r.sequence_length,
	       r.sequence_length, seq+r.start);
	printf("OUT #%d %d+%d %.*s\n",
	       rnum, r.position, r.sequence_length,
	       r.sequence_length, fseq+r.start);
#endif

	/* Check first good base */
	unique = 0;
	for (i = r.start+1; i < r.end && unique < MIN_OVERLAP; i++) {
#ifdef DEBUG
	    /* printf("<Pos %d: %c\n", i, fseq[i-1]); */
#endif
	    if (legal_chars[(unsigned)fseq[i-1]])
		unique++;
	}
	first = i - r.start + r.position - 1;

	/* Check for last good base */
	unique = 0;
	for (i = r.end-1; i > r.start && unique < MIN_OVERLAP; i--) {
#ifdef DEBUG
	    /* printf(">Pos %d: %c\n", i, fseq[i-1]); */
#endif
	    if (legal_chars[(unsigned)fseq[i-1]])
		unique++;
	}
	last = i - r.start + r.position - 1;

#ifdef DEBUG
	printf("Seq %d, first=%d, last=%d\n",
	       rnum, first, last);
#endif

	if (clips) {
	    clips[rnum-1].left = first;
	    clips[rnum-1].right = last;
	}

	if (last >= first)
	    memset(&valid[first-1], 1, last-first+1);

	xfree(seq);
	xfree(fseq);
    }

    gaps = coverage2contig_regions(valid, clen);
    xfree(valid);
    xfree(cons);

    return gaps;
}

/*
 * Uses read-pair information to confirm whether a gap appears to be
 * valid (or invalid).
 */
static void confirm_gaps(GapIO *io, int contig, Array gaps) {
    int i;
    template_c **tarr;

    /* Analyse templates for read-pair data */
    tarr = init_template_checks(io, 1, &contig, 1);

    for (i = 0; i <= Ntemplates(io); i++) {
	if (!tarr[i])
	    continue;
	/* ignore custom primer type info, eg w2k is either fwd or rev. */
	tarr[i]->oflags |= TEMP_OFLAG_IGNORE_PTYPE34 | TEMP_OFLAG_INTERDIST;
    }

    /* Now process gaps validating by read-pair */
    for (i = 0; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);
	int score;

	if (gap->deleted)
	    continue;

	score = check_read_pairs(io, contig, tarr, gap);
	printf("Gap %d..%d,  validity %d\n",
	       gap->start, gap->end, score);
    }

    uninit_template_checks(io, tarr);
}

/*
 * Steps through the gaps array identifying which reading numbers are
 * adjancent to the hole. This is an optimisation as it allows up to
 * randomly hop into our contig at specific positions.
 */
static void map_gaps(GapIO *io, int contig, Array gaps) {
    int rnum, i;

    rnum = io_clnbr(io, contig);
    for (i = 0; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);

	if (gap->deleted)
	    continue;

	/* left end */
	while (rnum && io_relpos(io, rnum) < gap->start)
	    rnum = io_rnbr(io, rnum);

	gap->rnum = rnum ? rnum : io_crnbr(io, contig);
    }
}

/*
 * Merges gaps if they're close together, within min_distance apart.
 */
static void merge_gaps(Array gaps, int min_distance) {
    int i;
    contig_region_t *last_gap;

    if (ArrayMax(gaps) == 0)
	return;

    last_gap = arrp(contig_region_t, gaps, 0);
    
    for (i = 1; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);

	if (gap->start - last_gap->end < min_distance) {
	    printf("Merging gap %d..%d with %d..%d\n",
		   last_gap->start, last_gap->end,
		   gap->start, gap->end);

	    last_gap->end = gap->end;
	    gap->deleted = 1;
	} else {
	    last_gap = gap;
	}
    }
}


/*
 * Analyses gaps to work out where to break. We'll produce a left contig,
 * a right contig, and possibly multiple single-read contigs for readings
 * that are contained entirely within the gap section.
 *
 * Fills out 'ds' with a list of breaks to make. Each item is in itself
 * another list with the last item 'n' being the location to break on and
 * items 1 to n-1 being readings in the middle of the gap that can be
 * moved to their own contig (either single read contigs or as a group).
 */
static void break_gaps(GapIO *io, int contig, Array gaps, clip_pos_t *clips,
		       dstring_t *ds) {
    int rnum, i;

    rnum = io_clnbr(io, contig);
    for (i = 0; i < ArrayMax(gaps); i++) {
	contig_region_t *gap = arrp(contig_region_t, gaps, i);
	int break_point = 0;

	if (gap->deleted)
	    continue;

	printf("Gap from %d to %d\n", gap->start, gap->end);

	if (gap->valid) {
	    printf("Skipping as marked as valid\n");
	    continue;
	}

	/* 
	 * Identify readings that after clipping are entirely contained
	 * within the region. These are candidates for disassembling
	 * completely
	 */
	dstring_appendf(ds, " {");
	while (rnum && io_relpos(io, rnum) < gap->end) {
	    clip_pos_t *c = &clips[rnum-1];

	    if (c->left >= gap->start && c->right <= gap->end) {
		printf("  Read #%d to self-contig\n", rnum);
		dstring_appendf(ds, " %s", io_rname(io, rnum));

	    /* first reading to overlap right edge = break point */
	    } else if (c->left >= gap->start && !break_point)
		break_point = rnum;

	    rnum = io_rnbr(io, rnum);
	}

	if (!break_point && rnum)
	    break_point = rnum;

	if (break_point) {
	    dstring_appendf(ds, " %s", io_rname(io, break_point));
	    printf("  New starting point for right contig = #%d\n",
		   break_point);
	}
	dstring_appendf(ds, "}");
    }
}

void auto_break_single_contig(GapIO *io, int contig, int start, int end,
			      int tw, double filter_score, double depth,
			      dstring_t *ds) {
    Array gaps;
    clip_pos_t *clips;

    printf("=== Checking contig %d (#%d), len %d ===\n",
	   contig, io_clnbr(io, contig), io_clength(io, contig));
    clips = (clip_pos_t *)xcalloc(NumReadings(io), sizeof(*clips));

    gaps = suspect_joins(io, contig, tw, filter_score, depth, clips);
    merge_gaps(gaps, MIN_OVERLAP*2);
    map_gaps(io, contig, gaps);
    dump_gaps(gaps);
    confirm_gaps(io, contig, gaps);
    break_gaps(io, contig, gaps, clips, ds);

    ArrayDestroy(gaps);
    xfree(clips);
}

dstring_t *auto_break_contigs(GapIO *io, int argc, contig_list_t *argv,
			      double filter_score, int by_consensus) {
    int tw, i;
    double gc;
    int depth;

    dstring_t *ds = dstring_create(NULL);

    if (by_consensus) {
	tw = word_count_cons(io, &gc); depth=1; /* Consensus based clipping */
    } else {
	tw = word_count(io, &gc, &depth); /* Sequence based clipping */
    }

    /*
    print_counts(filter_score * depth);
    print_bins();
    */

    for (i = 0; i < argc; i++) {
	auto_break_single_contig(io, argv[i].contig, argv[i].start,
				 argv[i].end, tw, filter_score, depth, ds);
    }

    return ds;
}
