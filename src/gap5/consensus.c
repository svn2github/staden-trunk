#include <xalloc.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <float.h>

/* See consensus.txt for discussions on these algorithms */

#include <tg_gio.h>

#include "consensus.h"

#define CONS_BLOCK_SIZE 1024

/*
 * Empirically derived constants for various sequencing technologies, with
 * caveats from local alignment artifacts and gap penalties:
 *
 * Solexa overcall:  1/43406
 * Solexa undercall: 1/28135
 * Solexa subst:     1/28.74 (96.5% accurate)
 */


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
    double slx_overcall_prob  = 1.0/43500, lover,  lomover;
    double slx_undercall_prob = 1.0/28000, lunder, lomunder;

    double (*cvec)[4]; /* cvec[0-3] = A,C,G,T */
    double (*pvec)[2]; /* pvec[0] = gap, pvec[1] = base */

    /* log overcall, log one minus overcall, etc */
    lover    = log(slx_overcall_prob);
    lomover  = log(1-slx_overcall_prob);
    lunder   = log(slx_undercall_prob);
    lomunder = log(1-slx_undercall_prob);

    /* Initialise */
    if (NULL == (cvec = (double (*)[4])calloc(len, 4 * sizeof(double))))
	return -1;
    if (NULL == (pvec = (double (*)[2])calloc(len, 2 * sizeof(double))))
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
	lookup['*'] = 4;

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

    /* Accumulate... (computes the products via sum of logs) */
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
		pvec[sp-start+j][0] += lover;
		pvec[sp-start+j][1] += lomover;
		break;

	    case 1:
		cvec[sp-start+j][0] += lo2r[q];
		cvec[sp-start+j][1] += lo2l[q];
		cvec[sp-start+j][2] += lo2r[q];
		cvec[sp-start+j][3] += lo2r[q];
		pvec[sp-start+j][0] += lover;
		pvec[sp-start+j][1] += lomover;
		break;

	    case 2:
		cvec[sp-start+j][0] += lo2r[q];
		cvec[sp-start+j][1] += lo2r[q];
		cvec[sp-start+j][2] += lo2l[q];
		cvec[sp-start+j][3] += lo2r[q];
		pvec[sp-start+j][0] += lover;
		pvec[sp-start+j][1] += lomover;
		break;

	    case 3:
		cvec[sp-start+j][0] += lo2r[q];
		cvec[sp-start+j][1] += lo2r[q];
		cvec[sp-start+j][2] += lo2r[q];
		cvec[sp-start+j][3] += lo2l[q];
		pvec[sp-start+j][0] += lover;
		pvec[sp-start+j][1] += lomover;
		break;

	    case 4:
		pvec[sp-start+j][0] += lomunder;
		pvec[sp-start+j][1] += lunder;
		break;

	    default: /* 5 */
		pvec[sp-start+j][0] += lover;
		pvec[sp-start+j][1] += lomover;
		break;
	    }
	}
	cache_decr(io, sorig);

	if (s != sorig)
	    free(s);
    }


    /* and speculate */
    for (i = 0; i < len; i++) {
	double probs[6], tot2[4], max;
	double pad_prob, base_prob;
	int j;

	/* Gap or base? Work out pad probability initially */
	/* For this the sum differences is basically the log-odds score */
	cons[i].scores[4] = 10*(pvec[i][0] - pvec[i][1]);
	pad_prob = pow(10, cons[i].scores[4] / 10.0) /
	    (pow(10, cons[i].scores[4] / 10.0) + 1);
	base_prob = 1-pad_prob;

	/*
	 * And now which base type it may be.
	 * Here probs[] hold the numerators with the denominators being
	 * sum(probs[0..4]). It cancels out though so we don't need it.
	 */
	for (j = 0; j < 4; j++) {
	    probs[j] = exp(cvec[i][j]);
	    tot2[j] = 0;
	}

	/*
	 * tot2 here is for computing 1-prob without hitting rounding
	 * issues caused by 1/<some small number>.
	 *
	 * If prob is a/(a+b+c+d) then 1-a/(a+b+c+d) = (b+c+d)/(a+b+c+d).
	 * Factoring in the base_prob B then prob is B.a/(a+b+c+d).
	 * Hence 1-prob = (a.(1-B)+b+c+d)/(a+b+c+d).
	 */
	for (j = 0; j < 4; j++) {
	    int k;
	    for (k = 0; k < 4; k++)
		if (j != k)
		    tot2[k] += probs[j];
		else
		    tot2[k] += probs[j] * pad_prob;
	}

	/*
	 * Normalise.
	 * We compute p as probs[j]/total_prob.
	 * Then turn this into a log-odds via log(p/(1-p)).
	 * 1-p has some precision issues in floating point, so we use tot2
	 * as described above as an alternative way of computing it. This
	 * gives:
	 *
	 * log(p/(1-p)) = log( (probs[j]/total) / (tot2[j] / total))
	 *              = log(probs[j]/tot2[j])
	 *              = log(probs[j]) - log(tot2[j])
	 */
	for (j = 0; j < 4; j++) {
	    cons[i].scores[j] = 10 / log(10) *
		(log(base_prob * probs[j]) - log(tot2[j]));
	}

	/* And call */
	if (cons[i].scores[4] > 0) {
	    cons[i].call = 4;
	    max = cons[i].scores[4];
	} else {
	    max = -4.7; /* Consensus cutoff */
	    cons[i].call = 5; /* N */
	    for (j = 0; j < 4; j++) {
		if (max < cons[i].scores[j]) {
		    max = cons[i].scores[j];
		    cons[i].call = j;
		}
	    }
	}

	if (max >  99) max =  99;
	if (max < -127) max = -127;
	cons[i].phred = lo2ph[(int)max];
    }

    if (r) free(r);
    return 0;
}
