#include <xalloc.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <assert.h>

/* See consensus.txt for discussions on these algorithms */

#include <tg_gio.h>

#include "consensus.h"

#define CONS_BLOCK_SIZE 4096

#define LOG10        2.30258509299404568401
#define TENOVERLOG10 4.34294481903251827652

#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

static unsigned char lookup[256], lookup_done = 0;

static double logodds2log[256];
static double *lo2l = logodds2log+128;

static double logodds2remainder3[256];
static double *lo2r3 = logodds2remainder3+128;

static double logodds2remainder1[256];
static double *lo2r1 = logodds2remainder1+128;

static unsigned char logodds2phred[256];
static unsigned char *lo2ph = logodds2phred+128;

static int calculate_consensus_bit(GapIO *io, tg_rec contig,
				   int start, int end,
				   consensus_t *cons);
static int calculate_consensus_bit_het(GapIO *io, tg_rec contig,
				       int start, int end,
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
int calculate_consensus_simple2(GapIO *io, tg_rec contig, int start, int end,
				char *con, float *qual) {
    int i, j;
    consensus_t q[CONS_BLOCK_SIZE];
    
    /* Compute in small ranges */
    for (i = start; i <= end; i += CONS_BLOCK_SIZE) {
	int st = i;
	int en = st + CONS_BLOCK_SIZE-1; /* inclusive range */
	if (en > end)
	    en = end;

	if (0 != calculate_consensus_bit_het(io, contig, st, en, q)) {
	    for (j = 0; j < en-st; j++) {
		if (con)
		    con[i-start+j] = 'N';
		if (qual)
		    qual[i-start+j] = 0;
	    }
	    return -1;
	}

	for (j = 0; j < en-st+1; j++) {
	    if (q[j].call == 6) {
		if (con)
		    con[i-start+j] = ' ';
		if (qual)
		    qual[i-start+j] = 0;
	    } else {
		if (con)
		    con[i-start+j] = "ACGT*N"[q[j].call];
		if (qual)
		    qual[i-start+j] = q[j].scores[q[j].call];
	    }
	}
    }

    return 0;
}

int calculate_consensus_simple(GapIO *io, tg_rec contig, int start, int end,
			       char *con, float *qual) {
    int i, nr;
    rangec_t *r;
    contig_t *c;
    int left, right;

    /* Don't bother caching for temporary io structs */
    if (io->base) {
	return calculate_consensus_simple2(io, contig, start, end, con, qual);
    }

    //printf("Calculate_consensus_simple(contig=%"PRIrec", range=%d..%d)\n",
    //       contig, start, end);

    /*
     * We cache consensus in the first bin smaller than a specific size
     * (CONS_BIN_SIZE), chosen to fall between two powers of two to allow
     * a certain degree of editing without switching bins commonly.
     */
    c = (contig_t *)cache_search(io, GT_Contig, contig);
    //contig_dump_ps(io, &c, "/tmp/tree.ps");
    r = contig_bins_in_range(io, &c, start, end,
			     CSIR_SORT_BY_X | CSIR_LEAVES_ONLY,
			     CONS_BIN_SIZE, &nr);

    /* Skip through spanning bins returned */
    left = start;
    for (i = 0; i < nr; i++) {
	bin_index_t *bin;
	//int f_a, f_b;

	/* Skip entirely overlapping bins */
	if (r[i].end < left)
	    continue;

	right = MIN(r[i].start-1, end);

	bin = (bin_index_t *)cache_search(io, GT_Bin, r[i].rec);

	if (r[i].start > left) {
	    printf("Filling missing region %d..%d\n", left, right);
	    calculate_consensus_simple2(io, contig, left, right,
					con  ? &con[left-start]  : NULL,
					qual ? &qual[left-start] : NULL); 
	}

	//f_a = r[i].pair_start;
	//f_b = r[i].pair_end;

	//	if (bin->rng /* && !(end < NMIN(bin->start_used, bin->end_used) ||
	//		        start > NMAX(bin->start_used, bin->end_used)) */
	//            ) {
	if (1) {
	    seq_t *s, *dup_s = NULL, seq;
	    range_t *cons_r = NULL;
	    int n;
	    int bstart, bend; /* consensus start/end in this bin */

	    /*
	     * Fetch cached consensus sequence, or recreate it if needed.
	     * It may exist, but be invalidated too.
	     */
	    if (bin->flags & BIN_CONS_CACHED) {
		for (n = 0; bin->rng && n < ArrayMax(bin->rng); n++) {
		    range_t *l = arrp(range_t, bin->rng, n);

		    if ((l->flags&GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISCONS) {
			cons_r = l;
			break;
		    }
		}
	    }

	    bstart = r[i].start;
	    bend   = r[i].end - 1;

	    if (cons_r && (bin->flags & BIN_CONS_VALID)) {
		//bstart = cons_r->start;
		//bend   = cons_r->end;
		s = (seq_t *)cache_search(io, GT_Seq, cons_r->rec);
		cache_incr(io, s);

	    } else {
		/* Not cached, or is cached but needs updating */
		range_t r2, *r_out;
		tg_rec recno;
		int comp;
		float *tmp_qual;

		//bstart = MAX(r[i].start, c->start);
		//bend   = MIN(r[i].end-1, c->end);

		/* Load existing sequence, or fill out in-memory new struct */
		if (cons_r && !io->read_only) {
		    printf("Recomputing cached cons from %d..%d +%d\n",
			   cons_r->start, cons_r->end, r[i].start);

		    s = (seq_t *)cache_search(io, GT_Seq, cons_r->rec);
		    s = cache_rw(io, s);
		    if (cons_r->start != bstart || cons_r->end   != bend) {
			size_t extra_len =
			    (s->name       ? strlen(s->name)       : 0) +
			    (s->trace_name ? strlen(s->trace_name) : 0) +
			    (s->alignment  ? strlen(s->alignment)  : 0) +
			    (bend - bstart + 1)*2;
			s = cache_item_resize(s, sizeof(*s) + extra_len);
		    }

		    if ((s->len < 0) ^ r[i].comp)
			s->len = -s->len;
		    
		} else {
		    printf("Creating new cached cons %d..%d bin #%"PRIrec"\n",
			   bstart, bend, bin->rec);
		    memset(&seq, 0, sizeof(seq));
		    seq.pos    = bstart;
		    seq.len    = bend - bstart + 1;
		    seq.seq    = calloc(bend - bstart + 1, 1);
		    seq.conf   = calloc(seq.len, 1);
		    seq.name   = "cons";
		    seq.format = SEQ_FORMAT_CNF1;
		    seq.flags  = 0;
		    seq.left   = 1;
		    seq.right  = bend - bstart + 1;

		    seq.name_len       = 4;
		    seq.mapping_qual   = 0;
		    seq.trace_name     = NULL;
		    seq.trace_name_len = 0;
		    seq.alignment      = NULL;
		    seq.alignment_len  = 0;

		    s = &seq;
		}

		/*
		 * s is now either a cached seq_t object or the address of a
		 * local variable, which may be used to create a new sequence
		 * provided we're in read/write mode.
		 */

		if (s != &seq)
		    cache_incr(io, s);

		/* Update consensus and quality */
		tmp_qual = calloc(bend - bstart + 1, sizeof(float));

		calculate_consensus_simple2(io, contig,
					    bstart, bend,
					    s->seq, tmp_qual);

		for (n = 0; n < ABS(s->len); n++) {
		    int q = rint(tmp_qual[n]);
		    if (q < 0)
			s->conf[n] = 0;
		    else if (q > 127)
			s->conf[n] = 127;
		    else
			s->conf[n] = q;
		}			
		free(tmp_qual);


		/* If cached and range changed, move sequence */
		if (cons_r) {
		    //assert(cons_r->start == bstart);
		    //assert(cons_r->end   == bend);
		    // if (!io->read_only) { ... }
		}
		
		/* If not cached, create a new range */
		if (!io->read_only && !cons_r) {
		    r2.start    = bstart;
		    r2.end      = bend;
		    r2.rec      = 0;
		    r2.pair_rec = 0;
		    r2.mqual    = 0;
		    r2.flags    = GRANGE_FLAG_TYPE_SINGLE;
		    r2.flags   |= GRANGE_FLAG_ISCONS;
			
		    bin = bin_add_range(io, &c, &r2, &r_out, &comp, 0);
		    /*
		     * Not a valid assertion if we uncomment the bstart/bend
		     * settings above.
		     */
		    // assert(bin->rec == r[i].rec);

		    s->bin       = bin->rec;
		    s->bin_index = r_out - ArrayBase(range_t, bin->rng);
		    if (comp) {
			s->len = -s->len;
			s->flags |= SEQ_COMPLEMENTED;
		    }
		    recno = sequence_new_from(io, s);
		    r_out->rec = recno;

		    free(s->seq);
		    free(s->conf);

		    if (s != &seq)
			cache_decr(io, s);

		    s = cache_search(io, GT_Seq, recno);
		    cache_incr(io, s);
		}

		if (!io->read_only) {
		    bin = cache_rw(io, bin);
		    bin->flags |= BIN_CONS_VALID | BIN_CONS_CACHED |
			BIN_BIN_UPDATED;
		}
	    }

	    if ((s->len < 0) ^ r[i].comp) {
		if (s != &seq) {
		    cache_decr(io, s);
		    s = dup_s = dup_seq(s);
		}
		complement_seq_t(s);
	    }

	    /*
	     * Copy the cached seq to our consensus/qual buffers.
	     * Note, it may not span the entire bin if the contig is
	     * short or this bin is at the beginning or end.
	     */
	    
	    if (start <= bstart) {
		if (con) {
		    if (MIN(bend, end) + 1 > bstart) {
			memcpy(&con[bstart - start], s->seq,
			       MIN(bend, end) - bstart + 1);
		    }
		}
		if (qual) {
		    int en = MIN(bend, end) - bstart + 1;
		    int N;
		    for (n = bstart - start, N=0; N < en; n++, N++)
			qual[n] = s->conf[N];
		}
	    } else {
		if (con) {
		    if (MIN(bend, end) + 1 > start) {
			memcpy(con, &s->seq[start - bstart],
			       MIN(bend, end) - start + 1);
		    }
		}
		if (qual) {
		    int en = MIN(bend, end) - start + 1;
		    int N;
		    for (n=0, N = start-bstart; n < en; n++, N++)
			qual[n] = s->conf[N];
		}
	    }

	    if (s == &seq) {
		/* read-only cheat, so free up some data now */
		free(s->seq);
		free(s->conf);
	    } else if (dup_s) {
		free(dup_s);
	    } else {
		cache_decr(io, s);
	    }
	}

	left = r[i].end;
    }
    if (r) free(r);

    if (end > left) {
	//printf("Filling missing region %d..%d\n", left, end);
	calculate_consensus_simple2(io, contig, left, end,
				    con  ? &con[left-start]  : NULL,
				    qual ? &qual[left-start] : NULL); 
    }

    cache_flush(io);

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
int calculate_consensus(GapIO *io, tg_rec contig, int start, int end,
			consensus_t *cons) {
    int i;

    /* Compute in small ranges */
    for (i = start; i <= end; i += CONS_BLOCK_SIZE) {
	int st = i;
	int en = st + CONS_BLOCK_SIZE-1;
	if (en > end)
	    en = end;

	if (0 != calculate_consensus_bit_het(io, contig, st, en, &cons[i-start]))
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
static int calculate_consensus_bit(GapIO *io, tg_rec contig,
				   int start, int end,
				   consensus_t *cons) {
    int i, j, nr;
    rangec_t *r = NULL;
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, contig);
    int len = end - start + 1;
    static int loop = 0;

    double (*cvec)[4]; /* cvec[0-3] = A,C,G,T */
    double (*pvec)[2]; /* pvec[0] = gap, pvec[1] = base */
    char *perfect; /* quality=100 bases */
    int *depth, vst, ven;
    //cstat *cs;

    /* Initialise */
    if (NULL == (cvec = (double (*)[4])calloc(len, 4 * sizeof(double))))
	return -1;
    if (NULL == (pvec = (double (*)[2])calloc(len, 2 * sizeof(double))))
	return -1;
    if (NULL == (depth = (int *)calloc(len, sizeof(int))))
	return -1;
    if (NULL == (perfect = (char *)calloc(len, sizeof(char))))
	return -1;

    if (!lookup_done) {
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
	    lo2l[i]  = log(p);
	    lo2r3[i] = log((1-p)/3);
	    lo2r1[i] = log((1-p));
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
	    char base, base_l;
	    int qual;

	    if (sp+j > end)
		continue;

	    sequence_get_base(io, &s, j+off, &base, &qual, NULL, 0);

	    base_l = lookup[base];
	    if (base_l < 5 && qual == 100)
		perfect[sp-start+j] |= (1<<base_l);
	    
	    /* Quality 0 should never be permitted as it breaks the math */
	    if (qual < 1)
		qual = 1;

	    switch (base_l) {
	    case 0:
		cvec[sp-start+j][0] += lo2l [qual];
		cvec[sp-start+j][1] += lo2r3[qual];
		cvec[sp-start+j][2] += lo2r3[qual];
		cvec[sp-start+j][3] += lo2r3[qual];
		pvec[sp-start+j][0] += lo2r1[qual];
		pvec[sp-start+j][1] += lo2l [qual];
		break;

	    case 1:
		cvec[sp-start+j][0] += lo2r3[qual];
		cvec[sp-start+j][1] += lo2l [qual];
		cvec[sp-start+j][2] += lo2r3[qual];
		cvec[sp-start+j][3] += lo2r3[qual];
		pvec[sp-start+j][0] += lo2r1[qual];
		pvec[sp-start+j][1] += lo2l [qual];
		break;

	    case 2:
		cvec[sp-start+j][0] += lo2r3[qual];
		cvec[sp-start+j][1] += lo2r3[qual];
		cvec[sp-start+j][2] += lo2l [qual];
		cvec[sp-start+j][3] += lo2r3[qual];
		pvec[sp-start+j][0] += lo2r1[qual];
		pvec[sp-start+j][1] += lo2l [qual];
		break;

	    case 3:
		cvec[sp-start+j][0] += lo2r3[qual];
		cvec[sp-start+j][1] += lo2r3[qual];
		cvec[sp-start+j][2] += lo2r3[qual];
		cvec[sp-start+j][3] += lo2l [qual];
		pvec[sp-start+j][0] += lo2r1[qual];
		pvec[sp-start+j][1] += lo2l [qual];
		break;

	    case 5: /* N => no idea what base, but it is a base still */
		/* Is this a valid thing to do? */
		pvec[sp-start+j][0] += lo2r1[qual];
		pvec[sp-start+j][1] += lo2l [qual];
		break;

	    case 4: /* gap */
		pvec[sp-start+j][0] += lo2l [qual];
		pvec[sp-start+j][1] += lo2r1[qual];
		break;
	    }

	    depth[sp-start+j]++;
	}

	//cache_decr(io, sorig);

	if (s != sorig)
	    free(s);
    }
#endif

    /* Determine valid ranges, if necessary */
    if (depth[0] == 0 && loop == 0) {
	/* Check left edge */
	loop = 1;
	consensus_valid_range(io, contig, &vst, NULL);
	loop = 0;
	if (vst > start) {
	    vst -= start;
	} else {
	    vst = 0;
	}
    } else {
	vst = 0;
    }

    if (depth[len-1] == 0 && loop == 0) {
	/* Check right edge */
	loop = 1;
	consensus_valid_range(io, contig, NULL, &ven);
	loop = 0;
	if (ven > start) {
	    ven -= start;
	} else {
	    ven = 0;
	}
    } else {
	ven = len-1;
    }

    //    printf("range=%d..%d => i from %d..%d valid %d..%d\n",
    //	   start, end, 0, len-1, vst, ven);

    /* and speculate */
    for (i = 0; i < len; i++) {
	double probs[6], tot2[4], max;
	double pad_prob, base_prob;

	/* Perfect => manual edit at 100% */
	if (perfect[i]) {
	    cons[i].scores[0] = -127;
	    cons[i].scores[1] = -127;
	    cons[i].scores[2] = -127;
	    cons[i].scores[3] = -127;
	    cons[i].scores[4] = -127;
	    cons[i].scores[5] = 0; /* N */

	    cons[i].phred = 255;

	    switch (perfect[i]) {
	    case 1<<0:
		cons[i].call = 0;
		break;
	    case 1<<1:
		cons[i].call = 1;
		break;
	    case 1<<2:
		cons[i].call = 2;
		break;
	    case 1<<3:
		cons[i].call = 3;
		break;
	    case 1<<4:
		cons[i].call = 4;
		break;

	    default:
		/* Multiple bases marked with max quality */
		cons[i].call  = 5;
		cons[i].phred = 0;
		break;
	    }

	    cons[i].scores[cons[i].call] = cons[i].phred;
	    cons[i].depth = depth[i];

	    continue;
	}

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
	} else if (i < vst || i > ven) {
	    cons[i].call = 6; /* no base */
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
    free(perfect);

    return 0;
}


#define P_HET 1e-6

static double prior[25];     /* Sum to 1.0 */
static double lprior15[15];  /* 15 combinations of {ACGT*} */
static double acc_prior[25]; /* Cumulative values up to 32768 */

/* Precomputed matrices for the consensus algorithm */
static double pMM[101], p__[101], p_M[101];

static void consensus_init(double p_het) {
    int i;
    double acc = 0;

    memset(lookup, 5, 256*sizeof(lookup[0]));
    lookup['A'] = lookup['a'] = 0;
    lookup['C'] = lookup['c'] = 1;
    lookup['G'] = lookup['g'] = 2;
    lookup['T'] = lookup['t'] = 3;
    lookup['*'] =               4;

    for (i = 0; i < 25; i++)
	prior[i] = p_het / 20;
    prior[0] = prior[6] = prior[12] = prior[18] = prior[24] = (1-p_het)/5;

    lprior15[0]  = log(prior[0]);
    lprior15[1]  = log(prior[1]*2);
    lprior15[2]  = log(prior[2]*2);
    lprior15[3]  = log(prior[3]*2);
    lprior15[4]  = log(prior[4]*2);
    lprior15[5]  = log(prior[6]);
    lprior15[6]  = log(prior[7]*2);
    lprior15[7]  = log(prior[8]*2);
    lprior15[8]  = log(prior[9]*2);
    lprior15[9]  = log(prior[12]);
    lprior15[10] = log(prior[13]*2);
    lprior15[11] = log(prior[14]*2);
    lprior15[12] = log(prior[18]);
    lprior15[13] = log(prior[19]*2);
    lprior15[14] = log(prior[24]);

    for (i = 0; i < 25; i++) {
	acc += prior[i];
	acc_prior[i] = acc * (0xffffff+1);
    }

    for (i = 1; i < 101; i++) {
	double prob = 1 - pow(10, -i / 10.0);
	 
	pMM[i] = log(prob/5);
	p__[i] = log((1-prob)/20);
	p_M[i] = log(prob/10 + (1-prob)/40);
    }
    pMM[0] = pMM[1];
    p__[0] = p__[1];
    p_M[0] = p_M[1];
}

/*
 * As above, but allowing for the possibility of the data being
 * heterozygous.
 */
static int calculate_consensus_bit_het(GapIO *io, tg_rec contig,
				       int start, int end,
				       consensus_t *cons) {
    int i, j, nr;
    rangec_t *r = NULL;
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, contig);
    int len = end - start + 1;
    static int loop = 0;
    static int init_done =0;
    static double q2p[101];
    double min_e_exp = DBL_MIN_EXP * log(2) + 1;

    double (*scores)[15];
    double (*sumsC)[6], *sumsE;
    char *perfect; /* quality=100 bases */
    int *depth, vst, ven;

    /* Map the 15 possible combinations to 1-base or 2-base encodings */
    static int map_sing[15] = {0, 5, 5, 5, 5,
			          1, 5, 5, 5,
			             2, 5, 5,
			                3, 5,
			                   4};
    static int map_het[15] = {0,  1,  2,  3,  4,
			          6,  7,  8,  9,
			             12, 13, 14,
			                 18, 19,
			                     24};

    if (!init_done) {
	init_done = 1;
	consensus_init(P_HET);

	for (i = 0; i <= 100; i++) {
	    q2p[i] = pow(10, -i/10.0);
	}
    }

    /* Initialise */
    if (NULL == (scores = (double (*)[15])calloc(len, 15 * sizeof(double))))
	return -1;
    if (NULL == (sumsC = (double (*)[6])calloc(len, 6 * sizeof(double))))
	return -1;
    if (NULL == (sumsE = (double *)calloc(len, sizeof(double))))
	return -1;

    if (NULL == (depth = (int *)calloc(len, sizeof(int))))
	return -1;
    if (NULL == (perfect = (char *)calloc(len, sizeof(char))))
	return -1;

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
	    char base, base_l;
	    int qual;
	    double *S, MM, __, _M, qe;

	    if (sp+j > end)
		continue;

	    sequence_get_base(io, &s, j+off, &base, &qual, NULL, 0);

	    base_l = lookup[base];
	    if (base_l < 5 && qual == 100)
		perfect[sp-start+j] |= (1<<base_l);
	    
	    /* Quality 0 should never be permitted as it breaks the math */
	    if (qual < 1)
		qual = 1;
	    if (qual > 100)
		qual = 100;

	    S = scores[sp-start+j];
	    __ = p__[qual];
	    MM = pMM[qual];
	    _M = p_M[qual];

	    qe = q2p[qual];
	    sumsE[sp-start+j] += qe;
	    sumsC[sp-start+j][base_l] += 1 - qe;

	    switch (base_l) {
	    case 0:
		S[0] += MM; S[1 ]+= _M; S[2 ]+= _M; S[3 ]+= _M; S[4 ]+= _M;
		            S[5 ]+= __; S[6 ]+= __; S[7 ]+= __; S[8 ]+= __;
			                S[9 ]+= __; S[10]+= __; S[11]+= __; 
					            S[12]+= __; S[13]+= __; 
						                S[14]+= __;
		break;

	    case 1:
		S[0] += __; S[1 ]+= _M; S[2 ]+= __; S[3 ]+= __; S[4 ]+= __;
		            S[5 ]+= MM; S[6 ]+= _M; S[7 ]+= _M; S[8 ]+= _M;
			                S[9 ]+= __; S[10]+= __; S[11]+= __; 
					            S[12]+= __; S[13]+= __; 
						                S[14]+= __;
		break;

	    case 2:
		S[0] += __; S[1 ]+= __; S[2 ]+= _M; S[3 ]+= __; S[4 ]+= __;
		            S[5 ]+= __; S[6 ]+= _M; S[7 ]+= __; S[8 ]+= __;
			                S[9 ]+= MM; S[10]+= _M; S[11]+= _M; 
					            S[12]+= __; S[13]+= __; 
						                S[14]+= __;
		break;

	    case 3:
		S[0] += __; S[1 ]+= __; S[2 ]+= __; S[3 ]+= _M; S[4 ]+= __;
		            S[5 ]+= __; S[6 ]+= __; S[7 ]+= _M; S[8 ]+= __;
			                S[9 ]+= __; S[10]+= _M; S[11]+= __; 
					            S[12]+= MM; S[13]+= _M; 
						                S[14]+= __;
		break;

	    case 4:
		S[0] += __; S[1 ]+= __; S[2 ]+= __; S[3 ]+= __; S[4 ]+= _M;
		            S[5 ]+= __; S[6 ]+= __; S[7 ]+= __; S[8 ]+= _M;
			                S[9 ]+= __; S[10]+= __; S[11]+= _M; 
					            S[12]+= __; S[13]+= _M; 
						                S[14]+= MM;
		break;

	    case 5: /* N => should so something to non-pad vs pad scores */
		break;
	    }

	    depth[sp-start+j]++;
	}

	//cache_decr(io, sorig);

	if (s != sorig)
	    free(s);
    }

    /* Determine valid ranges, if necessary */
    if (depth[0] == 0 && loop == 0) {
	/* Check left edge */
	loop = 1;
	consensus_valid_range(io, contig, &vst, NULL);
	loop = 0;
	if (vst > start) {
	    vst -= start;
	} else {
	    vst = 0;
	}
    } else {
	vst = 0;
    }

    if (depth[len-1] == 0 && loop == 0) {
	/* Check right edge */
	loop = 1;
	consensus_valid_range(io, contig, NULL, &ven);
	loop = 0;
	if (ven > start) {
	    ven -= start;
	} else {
	    ven = 0;
	}
    } else {
	ven = len-1;
    }

    /* and speculate */
    for (i = 0; i < len; i++) {
	double shift, max, max_het, *S, norm[15];
	int call, het_call, ph;
	double tot1, tot2;

	/* Perfect => manual edit at 100% */
	if (perfect[i]) {
	    cons[i].scores[0] = -127;
	    cons[i].scores[1] = -127;
	    cons[i].scores[2] = -127;
	    cons[i].scores[3] = -127;
	    cons[i].scores[4] = -127;
	    cons[i].scores[5] = 0; /* N */

	    cons[i].phred = 255;

	    switch (perfect[i]) {
	    case 1<<0:
		cons[i].call = 0;
		break;
	    case 1<<1:
		cons[i].call = 1;
		break;
	    case 1<<2:
		cons[i].call = 2;
		break;
	    case 1<<3:
		cons[i].call = 3;
		break;
	    case 1<<4:
		cons[i].call = 4;
		break;

	    default:
		/* Multiple bases marked with max quality */
		cons[i].call  = 5;
		cons[i].phred = 0;
		break;
	    }

	    cons[i].scores[cons[i].call] = cons[i].phred;
	    cons[i].depth = depth[i];

	    continue;
	}

	/*
	 * Scale numbers so the maximum score is 0. This shift is essentially 
	 * a multiplication in non-log scale to both numerator and denominator,
	 * so it cancels out. We do this to avoid calling exp(-large_num) and
	 * ending up with norm == 0 and hence a 0/0 error.
	 *
	 * Can also generate the base-call here too.
	 */
	S = scores[i];
	shift = -DBL_MAX;
	max = -DBL_MAX;
	max_het = -DBL_MAX;
	for (j = 0; j < 15; j++) {
	    S[j] += lprior15[j];
	    if (shift < S[j]) {
		shift = S[j];
		//het_call = j;
	    }

	    /* Only call pure AA, CC, GG, TT, ** for now */
	    if (j != 0 && j != 5 && j != 9 && j != 12 && j != 14) {
		if (max_het < S[j]) {
		    max_het = S[j];
		    het_call = j;
		}
		continue;
	    }

	    if (max < S[j]) {
		max = S[j];
		call = j;
	    }
	}

	/*
	 * Shift and normalise.
	 * If call is, say, b we want p = b/(a+b+c+...+n), but then we do
	 * p/(1-p) later on and this has exceptions when p is very close
	 * to 1.
	 *
	 * Hence we compute b/(a+b+c+...+n - b) and
	 * rearrange (p/norm) / (1 - (p/norm)) to be p/norm2.
	 */
	for (j = 0; j < 15; j++) {
	    S[j] -= shift;
	    if (S[j] > min_e_exp) {
		S[j] = exp(S[j]);
	    } else {
		S[j] = DBL_MIN;
	    }
	    norm[j] = 0;
	}

	tot1 = tot2 = 0;
	for (j = 0; j < 15; j++) {
	    norm[j]    += tot1;
	    norm[14-j] += tot2;
	    tot1 += S[j];
	    tot2 += S[14-j];
	}

	for (j = 0; j < 15; j++) {
	    if (norm[j] == 0)
		norm[j] = DBL_MIN;
	}

	/* And store result */
	if (depth[i]) {
	    double m;

	    cons[i].depth = depth[i];

	    cons[i].call     = map_sing[call];
	    ph = -TENOVERLOG10 * log(norm[call]) + .5;
	    cons[i].phred = ph > 255 ? 255 : (ph < 0 ? 0 : ph);

	    cons[i].het_call = map_het[het_call];
	    ph = TENOVERLOG10 * (log(S[het_call]) - log(norm[het_call])) + .5;
	    cons[i].scores[6] = ph;

	    /* AA */
	    ph = TENOVERLOG10 * (log(S[0]) - log(norm[0])) + .5;
	    cons[i].scores[0] = ph;

	    /* CC */
	    ph = TENOVERLOG10 * (log(S[5]) - log(norm[5])) + .5;
	    cons[i].scores[1] = ph;

	    /* GG */
	    ph = TENOVERLOG10 * (log(S[9]) - log(norm[9])) + .5;
	    cons[i].scores[2] = ph;

	    /* TT */
	    ph = TENOVERLOG10 * (log(S[12]) - log(norm[12])) + .5;
	    cons[i].scores[3] = ph;

	    /* ** */
	    ph = TENOVERLOG10 * (log(S[14]) - log(norm[14])) + .5;
	    cons[i].scores[4] = ph;

	    /* N */
	    cons[i].scores[5] = 0; /* N */

	    /* Compute discrepancy score */
	    m = sumsC[i][0]+sumsC[i][1]+sumsC[i][2]+sumsC[i][3]+sumsC[i][4]
		- sumsC[i][cons[i].call];
	    cons[i].discrep = (m-sumsE[i])/sqrt(m+sumsE[i]);
	} else {
	    if (i < vst || i > ven)
		cons[i].call = 6; /* No base */
	    else
		cons[i].call = 5; /* N */
	    cons[i].het_call = 0;
	    cons[i].scores[0] = 0;
	    cons[i].scores[1] = 0;
	    cons[i].scores[2] = 0;
	    cons[i].scores[3] = 0;
	    cons[i].scores[4] = 0;
	    cons[i].scores[5] = 0;
	    cons[i].scores[6] = 0;
	    cons[i].phred = 0;
	    cons[i].depth = 0;
	    cons[i].discrep = 0;
	}
    }

    if (r) free(r);
    free(scores);
    free(depth);
    free(perfect);

    return 0;
}



/*
 * Returns the visible non-clipped portion of a contig. We use iterators
 * on the contig to find the first used and last used bases and return
 * these in start and end. Passing over NULL indicates that you are
 * uninterested in that end.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
int consensus_valid_range(GapIO *io, tg_rec contig, int *start, int *end) {
    contig_iterator *ci;
    rangec_t *r;

    if (start) {
	int best = INT_MAX;

	ci = contig_iter_new(io, contig, 0, CITER_FIRST | CITER_ISTART,
			     CITER_CSTART, CITER_CEND);
	
	while ((r = contig_iter_next(io, ci))) {
	    seq_t *s;
	    int left;

	    if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
		continue;

	    if (r->start > best)
		break;

	    s = cache_search(io, GT_Seq, r->rec);

	    if ((s->len < 0) ^ r->comp) {
		left = r->start + ABS(s->len) - (s->right-1) - 1;
	    } else {
		left = r->start + s->left - 1;
	    }

	    //printf("Rec %"PRIrec" start %d clipped %d\n",
	    //       r->rec, r->start, left);

	    if (best > left)
		best = left;
	}

	contig_iter_del(ci);
	*start = best != INT_MAX ? best : 0;
    }

    if (end) {
	int best = INT_MIN;

	ci = contig_iter_new(io, contig, 0, CITER_LAST | CITER_IEND,
			     CITER_CSTART, CITER_CEND);
	
	while ((r = contig_iter_prev(io, ci))) {
	    seq_t *s;
	    int right;

	    if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
		continue;

	    if (r->end < best)
		break;

	    s = cache_search(io, GT_Seq, r->rec);

	    if ((s->len < 0) ^ r->comp) {
		right = r->start + ABS(s->len) - (s->left-1) - 1;
	    } else {
		right = r->start + s->right - 1;
	    }

	    //printf("Rec %"PRIrec" right %d clipped %d\n",
	    //	   r->rec, r->end, right);

	    if (best < right)
		best = right;
	}

	contig_iter_del(ci);
	*end = best != INT_MIN ? best : 0;
    }

    return 0;
}
