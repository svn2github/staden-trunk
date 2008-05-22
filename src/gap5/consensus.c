#include <xalloc.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <float.h>

#include <tg_gio.h>

#include "consensus.h"

#define CONS_BLOCK_SIZE 1024

static unsigned char lookup[256], lookup_done = 0;

static double logodds2log[256];
static double *lo2l = logodds2log+128;

static double logodds2remainder[256];
static double *lo2r = logodds2remainder+128;

static unsigned char logodds2phred[256];
static unsigned char *lo2ph = logodds2phred+128;

static int calculate_consensus_bit(GapIO *io, int contig, int start, int end,
				   consensus_t *cons);

/*
 * The consensus calculation function - rewritten for tgap style
 * databases.
 *
 * Here con/qual holds the consensus sequence and a single phred style
 * quality value. Either can be passed in as NULL if that array is not needed,
 * but if supplied it is up to the caller to ensure they are large enough.
 *
 * start and end specify the inclusive range, counting from base 0.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
#define Q_CUTOFF -4
int calculate_consensus_simple(GapIO *io, int contig, int start, int end,
			       char *con, float *qual) {
    int i, j;
    consensus_t q[CONS_BLOCK_SIZE];
    
    /* Compute in small ranges */
    for (i = start; i < end; i += CONS_BLOCK_SIZE) {
	int st = i;
	int en = st + CONS_BLOCK_SIZE-1; /* inclusive range */
	if (en > end)
	    en = end;

	if (0 != calculate_consensus_bit(io, contig, st, en, q)) {
	    for (j = 0; j < en-st; j++) {
		if (con)
		    con[i-start+j] = 'N';
		if (qual)
		    qual[i-start+j] = 0;
	    }
	    return -1;
	}

	for (j = 0; j < en-st+1; j++) {
	    if (q[j].scores[q[j].call] > Q_CUTOFF) {
		if (con)
		    con[i-start+j] = "ACGT*"[q[j].call];
		if (qual)
		    qual[i-start+j] = q[j].scores[q[j].call];
	    } else {
		if (con)
		    con[i-start+j] = 'N';
		if (qual)
		    qual[i-start+j] = 0;
	    }
	}
    }

    return 0;
}

/*
 * The consensus calculation function - rewritten for tgap style
 * databases.
 *
 * Similar to calculate_consensus_simple, but we present the full
 * probabilities instead of just the called base.
 * "cons" is filled out with consensus_t structs. It should be passed in
 * by the caller and is expected to be the correct length.
 *
 * start and end specify the inclusive range, counting from base 0.
 *
 * Returns 0 on success,
 *        -1 on failure
 *
 */
int calculate_consensus(GapIO *io, int contig, int start, int end,
			consensus_t *cons) {
    int i;
    
    /* Compute in small ranges */
    for (i = start; i < end; i += CONS_BLOCK_SIZE) {
	int st = i;
	int en = st + CONS_BLOCK_SIZE-1;
	if (en > end)
	    en = end;

	if (0 != calculate_consensus_bit(io, contig, st, en, &cons[i-start]))
	    return -1;
    }

    return 0;
}

/*
 * The core of the consensus algorithm.
 *
 * This uses more memory than the final result in cons[] as it briefly
 * holds an array of every sequence in the start..end range, so we
 * typically wrap this up in another function to call it in smaller
 * blocks.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
static int calculate_consensus_bit(GapIO *io, int contig, int start, int end,
				   consensus_t *cons) {
    int i, j, nr;
    rangec_t *r;
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, contig);
    int len = end - start + 1;
    double (*cvec)[6];

    /* Initialise */
    if (NULL == (cvec = (double (*)[6])calloc(len, 6 * sizeof(double))))
	return -1;

    if (!lookup_done) {
	int i;

	/* Character code */
	memset(lookup, 5, 256);
	lookup_done = 1;
	lookup['A'] = lookup['a'] = 0;
	lookup['C'] = lookup['c'] = 1;
	lookup['G'] = lookup['g'] = 2;
	lookup['T'] = lookup['t'] = 3;
	lookup['*'] = lookup[','] = 4;

	/* Log odds value to log(P) */
	for (i = -128; i < 128; i++) {
	    double p = 1 / (1 + pow(10, -i / 10.0));
	    lo2l[i] = log(p);
	    lo2r[i] = log((1-p)/3);
	    lo2ph[i] = 10*log(1+pow(10, i/10.0))/log(10)+0.4999;
	}
    }

    /* Find sequences visible */
    r = contig_seqs_in_range(io, &c, start, end, &nr);

    /* Accumulate... */
    for (i = 0; i < nr; i++) {
	int sp = r[i].start;
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r[i].rec);
	seq_t *sorig = s;
	int l = s->len > 0 ? s->len : -s->len;
	unsigned char *seq;
	int left, right;
	char *conf;

	/* Complement data on-the-fly */
	if ((s->len < 0) ^ r[i].comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	if (s->len < 0) {
	    sp += s->len+1;
	}

	seq = (unsigned char *)s->seq;
	conf = s->conf;
	left = s->left;
	right = s->right;
	
	if (sp < start) {
	    seq   += start - sp;
            conf  += start - sp;
	    l     -= start - sp;
	    left  -= start - sp;
	    right -= start - sp;
	    sp = start;
	}
	if (l > end - sp + 1)
	    l = end - sp + 1;
	if (left < 1)
	    left = 1;

	for (j = left-1; j < right; j++) {
	    int q = conf[j];

	    if (sp+j > end)
		continue;
	    
	    switch(lookup[seq[j]]) {
	    case 0:
		cvec[sp-start+j][0] += lo2l[q];
		cvec[sp-start+j][1] += lo2r[q];
		cvec[sp-start+j][2] += lo2r[q];
		cvec[sp-start+j][3] += lo2r[q];
		break;

	    case 1:
		cvec[sp-start+j][0] += lo2r[q];
		cvec[sp-start+j][1] += lo2l[q];
		cvec[sp-start+j][2] += lo2r[q];
		cvec[sp-start+j][3] += lo2r[q];
		break;

	    case 2:
		cvec[sp-start+j][0] += lo2r[q];
		cvec[sp-start+j][1] += lo2r[q];
		cvec[sp-start+j][2] += lo2l[q];
		cvec[sp-start+j][3] += lo2r[q];
		break;

	    case 3:
		cvec[sp-start+j][0] += lo2r[q];
		cvec[sp-start+j][1] += lo2r[q];
		cvec[sp-start+j][2] += lo2r[q];
		cvec[sp-start+j][3] += lo2l[q];
		break;
	    }
	}
	cache_decr(io, sorig);

	if (s != sorig)
	    free(s);
    }


    /* and speculate */
    for (i = 0; i < len; i++) {
	double probs[6], tot, tot2[4], max;
	int j;

	tot = 0;
	for (j = 0; j < 4; j++) {
	    probs[j] = exp(cvec[i][j]);
	    tot += probs[j];
	    tot2[j] = 0;
	}
	for (j = 0; j < 4; j++) {
	    int k;
	    for (k = 0; k < 4; k++)
		if (j != k)
		    tot2[k] += probs[j];
	}

	/* Normalise  */
	for (j = 0; j < 4; j++) {
	    double p, omp; /* omp = 1- p */

	    p = probs[j] = tot > 0 ? probs[j] / tot : 0;
	    omp = tot2[j] / tot;

	    cons[i].scores[j] = 10 * log(p / omp) / log(10);
	}
#if 0
	/* Normalise  */
	for (j = 0; j < 4; j++) {
	    double p = probs[j] = tot > 0 ? probs[j] / tot : 0;
	    if (1-p < DBL_EPSILON)
		p = 1-DBL_EPSILON;
	    cons[i].scores[j] = 10 * log(p / (1-p)) / log(10);
	}
#endif
	
	cons[i].scores[4] = -99; /* FIXME: gap probability */

	/* And call */
	max = -4.7;
	cons[i].call = 5; /* N */
	cons[i].scores[5] = -99;
	for (j = 0; j < 5; j++) {
	    if (max < cons[i].scores[j]) {
		max = cons[i].scores[j];
		cons[i].call = j;
	    }
	}

	if (max >  99) max =  99;
	if (max < -127) max = -127;
	cons[i].phred = lo2ph[(int)max];
    }

    if (r) free(r);
    return 0;
}
