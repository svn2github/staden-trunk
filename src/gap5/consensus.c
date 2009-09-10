#include <xalloc.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <assert.h>

/* See consensus.txt for discussions on these algorithms */

#include <tg_gio.h>

#include "consensus.h"

#define CONS_BLOCK_SIZE 1024

#define LOG10        2.30258509299404568401
#define TENOVERLOG10 4.34294481903251827652

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
    for (i = start; i <= end; i += CONS_BLOCK_SIZE) {
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
		    con[i-start+j] = "ACGT*N"[q[j].call];
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
    for (i = start; i <= end; i += CONS_BLOCK_SIZE) {
	int st = i;
	int en = st + CONS_BLOCK_SIZE-1;
	if (en > end)
	    en = end;

	if (0 != calculate_consensus_bit(io, contig, st, en, &cons[i-start]))
	    return -1;
    }

    return 0;
}


typedef struct {
    double base_qual[4];
    double gap;
    double base;
    int depth;
} cstat;

/*
 * Like contig_seqs_in_range, but we also sum the values into the cvec array
 * at the same time.
 *
 * The purpose of this is to allow partial or complete caching of the
 * consensus data, meaning we only need to load bin tracks instead of
 * all the sequences themselves.
 */
#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))
static int contig_consensus_in_range2(GapIO *io, int brec, int start, int end,
				      int offset, cstat *cv, int complement) {
    double slx_overcall_prob  = 1.0/4350000, lover,  lomover;
    double slx_undercall_prob = 1.0/2800000, lunder, lomunder;
    int i, n, f_a, f_b;
    bin_index_t *bin = get_bin(io, brec);
    range_t *l;

    /* log overcall, log one minus overcall, etc */
    lover    = log(slx_overcall_prob);
    lomover  = log(1-slx_overcall_prob);
    lunder   = log(slx_undercall_prob);
    lomunder = log(1-slx_undercall_prob);

    cache_incr(io, bin);
    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }

    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    /*
    printf("BIN #%d %d + %d..%d (%d..%d)\n",
	   brec, offset, bin->pos, bin->pos + bin->size,
	   bin->start_used, bin->end_used);
    */
    if (!(end < NMIN(bin->start_used, bin->end_used) ||
	  start > NMAX(bin->start_used, bin->end_used))
	&& bin->rng) {
	track_t *tr = NULL;
	cstat *loc;

	tr = bin_get_track(io, bin, TRACK_CONS_ARR);

	if (tr) {
	    loc = arrp(cstat, tr->data, 0);
	} else {
	    tr = bin_create_track(io, bin, TRACK_CONS_ARR);
	    //tr = track_create_fake(TRACK_CONS_ARR, 0);
	    tr->nitems    = bin->end_used - bin->start_used + 1;
	    tr->item_size = sizeof(cstat);
	    tr->data      = ArrayCreate(tr->item_size, tr->nitems);
	    //tr->flag     |= TRACK_FLAG_FREEME;

	    bin_add_track(io, &bin, tr);
    
	    loc = arrp(cstat, tr->data, 0);
	    memset(&loc[0], 0, (bin->end_used - bin->start_used + 1) *
		   sizeof(*loc));

	    for (n = 0; n < ArrayMax(bin->rng); n++) {
		l = arrp(range_t, bin->rng, n);

		if (l->flags & (GRANGE_FLAG_UNUSED | GRANGE_FLAG_ISANNO))
		    continue;

		if (1 || (NMAX(l->start, l->end) >= start
			  && NMIN(l->start, l->end) <= end)) {
		    seq_t *s = (seq_t *)cache_search(io, GT_Seq, l->rec);
		    seq_t *sorig = s;
		    int p;

		    if ((s->len < 0)) {
			s = dup_seq(s);
			complement_seq_t(s);
		    }

		    /*
		      printf("Seq %c%d %d..%d\n",
		      s == sorig ? '>' : '<',
		      l->rec, s->left, s->right);
		    */
		    for (p = s->left-1; p < s->right; p++) {
			char base;
			double q[4];
			int j, k;

			sequence_get_base4(io, &s, p, &base, q, 0);
			/*
			  printf("  %5d+%3d %2d %c %f %f %f %f\n",
			  offset, p + l->start, p,
			  base, q[0], q[1], q[2], q[3]);
			*/
		    
			k = p + l->start - bin->start_used;
			if (k < 0 || k >= bin->end_used - bin->start_used + 1) {
			    printf("Error: start/end used summary is invalid\n");
			}
			assert(k >= 0 &&  k < bin->end_used - bin->start_used + 1);

			switch (lookup[base]) {
			case 0: case 1: case 2: case 3: /* ACGT */
			    loc[k].base_qual[0] += q[0];
			    loc[k].base_qual[1] += q[1];
			    loc[k].base_qual[2] += q[2];
			    loc[k].base_qual[3] += q[3];

			    /* Small boost for called base to resolve ties */
			    loc[k].base_qual[lookup[base]] += 1e-5;

			    /* Fall through */
			default: /* N */
			    loc[k].gap  += lover;
			    loc[k].base += lomover;
			    break;

			case 4: /* gap */
			    loc[k].gap  += lomunder;
			    loc[k].base += lunder;
			    break;
			}

			loc[k].depth++;
		    }

		    if (s != sorig)
			free(s);
		}
	    }
	}
	

	/*
	printf("Bin %d summary\n", brec);
	for (i = 0; i <= bin->end_used - bin->start_used; i++) {
	    printf("    %4d %4d %2d %f %f %f %f\n",
		   i + bin->start_used,
		   loc[i].depth,
		   loc[i].base_qual[0],
		   loc[i].base_qual[1],
		   loc[i].base_qual[2],
		   loc[i].base_qual[3]);
	}
	*/
	
	for (i = 0; i <= bin->end_used - bin->start_used; i++) {
	    int j = NORM(i + bin->start_used);
	    if (j >= start && j <= end) {
		int k = j-start;
		cv[k].base_qual[0] += loc[i].base_qual[0];
		cv[k].base_qual[1] += loc[i].base_qual[1];
		cv[k].base_qual[2] += loc[i].base_qual[2];
		cv[k].base_qual[3] += loc[i].base_qual[3];
		cv[k].gap          += loc[i].gap;
		cv[k].base         += loc[i].base;
		cv[k].depth        += loc[i].depth;
	    }
	}
	
	//free(loc);
    }

    /* Recurse down bins */
    for (i = 0; i < 2 > 0; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;
	ch = get_bin(io, bin->child[i]);
	if (end   >= NMIN(ch->pos, ch->pos + ch->size-1) &&
	    start <= NMAX(ch->pos, ch->pos + ch->size-1)) {
	    contig_consensus_in_range2(io, bin->child[i], start, end,
				       NMIN(ch->pos, ch->pos + ch->size-1),
				       cv, complement);
	}
    }

    cache_decr(io, bin);
    return 0;
}

int contig_consensus_in_range(GapIO *io, contig_t **c, int start, int end,
			      cstat *cv) {
    return contig_consensus_in_range2(io, contig_get_bin(c), start, end,
				      contig_offset(io, c), cv, 0);
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
    rangec_t *r = NULL;
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, contig);
    int len = end - start + 1;
    double slx_overcall_prob  = 1.0/4350000, lover,  lomover;
    double slx_undercall_prob = 1.0/2800000, lunder, lomunder;

    double (*cvec)[4]; /* cvec[0-3] = A,C,G,T */
    double (*pvec)[2]; /* pvec[0] = gap, pvec[1] = base */
    int *depth;
    cstat *cs;

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
    if (NULL == (depth = (int *)calloc(len, sizeof(int))))
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
	    lo2ph[i] = 10*log(1+pow(10, i/10.0))/LOG10+0.4999;
	}
    }

#if 0
    /* Accumulate... (computes the products via sum of logs) */
    /* This bit may be cached */
    cs = (cstat *)malloc((end-start+1) * sizeof(*cs));
    memset(&cs[0], 0, (end-start+1) * sizeof(*cs));
    contig_consensus_in_range(io, &c, start, end, cs);

    for (i = 0; i < end-start+1; i++) {
	cvec[i][0] = cs[i].base_qual[0];
	cvec[i][1] = cs[i].base_qual[1];
	cvec[i][2] = cs[i].base_qual[2];
	cvec[i][3] = cs[i].base_qual[3];
	pvec[i][0] = cs[i].gap;
	pvec[i][1] = cs[i].base;
	depth[i]   = cs[i].depth;
    }

    free(cs);

#else

    /* Find sequences visible */
    r = contig_seqs_in_range(io, &c, start, end, 0, &nr);

    for (i = 0; i < nr; i++) {
	int sp = r[i].start;
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, r[i].rec);
	seq_t *sorig = s;
	int l = s->len > 0 ? s->len : -s->len;
	int left, right;
	int off = 0;

	//cache_incr(io, sorig);

	/* Complement data on-the-fly */
	if ((s->len < 0) ^ r[i].comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	left = s->left;
	right = s->right;
	
	if (sp < start) {
	    off    = start - sp;
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
	    char base;
	    double q[4];

	    if (sp+j > end)
		continue;
	    
	    sequence_get_base4(io, &s, j+off, &base, q, 0);

	    switch (lookup[base]) {
	    case 0: case 1: case 2: case 3: /* ACGT */
		cvec[sp-start+j][0] += q[0];
		cvec[sp-start+j][1] += q[1];
		cvec[sp-start+j][2] += q[2];
		cvec[sp-start+j][3] += q[3];

		/* Small boost for called base to resolve ties */
		cvec[sp-start+j][lookup[base]] += 1e-5;

		/* Fall through */
	    default: /* N */
		pvec[sp-start+j][0] += lover;
		pvec[sp-start+j][1] += lomover;
		break;

	    case 4: /* gap */
		pvec[sp-start+j][0] += lomunder;
		pvec[sp-start+j][1] += lunder;
		break;
	    }

	    depth[sp-start+j]++;
	}

	//cache_decr(io, sorig);

	if (s != sorig)
	    free(s);
    }
#endif


    /* and speculate */
    for (i = 0; i < len; i++) {
	double probs[6], tot2[4], max;
	double pad_prob, base_prob;
	int j;

	/* Gap or base? Work out pad probability initially */
	/* For this the sum differences is basically the log-odds score */
	if (pvec[i][0] && pvec[i][1]) {
	    cons[i].scores[4] = 10*(pvec[i][0] - pvec[i][1]);
	    pad_prob = pow(10, cons[i].scores[4] / 10.0) /
		(pow(10, cons[i].scores[4] / 10.0) + 1);
	    base_prob = 1-pad_prob;
	} else {
	    cons[i].scores[4] = 0;
	    pad_prob = 0;
	    base_prob = 1;
	}

	/* Shift so that we never attempt to exp() of all high -ve values */
	max = cvec[i][0];
	for (j = 1; j < 4; j++) {
	    if (max < cvec[i][j])
		max = cvec[i][j];
	}
	for (j = 0; j < 4; j++) {
	    cvec[i][j] -= max;
	}

	/*
	 * And now which base type it may be.
	 * Here probs[] hold the numerators with the denominators being
	 * sum(probs[0..4]). It cancels out though so we don't need it.
	 */
	for (j = 0; j < 4; j++) {
	    probs[j] = exp(cvec[i][j]);
	    if (probs[j] == 0)
		probs[j] = DBL_MIN;
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
	cons[i].scores[5] = 0;
	for (j = 0; j < 4; j++) {
	    cons[i].scores[j] = TENOVERLOG10 *
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

	cons[i].depth = depth[i];
    }

    if (r) free(r);
    free(cvec);
    free(pvec);
    free(depth);

    return 0;
}


/*
 * Finds the portion of a contig that has non-clipped data.
 * This is a somewhat crude method by just computing consensus at the ends
 * and trimming back the zero-depth regions.
 *
 * Specify start and end as pointers for the results. Passing over NULL
 * indicates that you are not interested in that end.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
#define CONS_SZ 1024
int consensus_valid_range(GapIO *io, int contig, int *start, int *end) {
    consensus_t cons[CONS_SZ];
    signed int pos, i;
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, contig);
    if (!c)
	return -1;

    if (start) {
	pos = contig_get_start(&c);
	do {
	    if (-1 == calculate_consensus_bit(io, contig,
					      pos, pos+CONS_SZ-1, cons))
		return -1;
	    for (i = 0; i < CONS_SZ; i++, pos++) {
		if (cons[i].depth)
		    break;
	    }
	} while (i == CONS_SZ && pos <= c->end);

	*start = pos;
    }

    if (end) {
	pos = contig_get_end(&c);
	do {
	    if (-1 == calculate_consensus_bit(io, contig,
					      pos-(CONS_SZ-1), pos, cons))
		return -1;
	    for (i = CONS_SZ-1; i >= 0; i--, pos--) {
		if (cons[i].depth)
		    break;
	    }
	} while (i == 0 && pos >= c->start);

	*end = pos;
    }

    return 0;
}
