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

static unsigned char lookup[256];

#if 0
static unsigned char lookup_done = 0;
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
#endif


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
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, contig);

    cache_incr(io, c);
    
    /* Compute in small ranges */
    for (i = start; i <= end; i += CONS_BLOCK_SIZE) {
	rangec_t *r = NULL;
	int nr;
	int st = i;
	int en = st + CONS_BLOCK_SIZE-1; /* inclusive range */
	if (en > end)
	    en = end;

	/* Find sequences visible */
	r = contig_seqs_in_range(io, &c, st, en, CSIR_SORT_BY_X, &nr);

	if (0 != calculate_consensus_bit_het(io, contig, st, en,
					     qual ? CONS_SCORES : 0,
					     r, nr, q)){
	    for (j = 0; j < en-st; j++) {
		if (con)
		    con[i-start+j] = 'N';
		if (qual)
		    qual[i-start+j] = 0;
	    }
	    cache_decr(io, c);
	    return -1;
	}

	if (r) free(r);

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

    cache_decr(io, c);
    return 0;
}

int calculate_consensus_simple(GapIO *io, tg_rec contig, int start, int end,
			       char *con, float *qual) {
    int i, j, nr;
    rangec_t *r;
    contig_t *c;
    int left, right;

   /* Don't bother caching for temporary io structs if short queries */
   if (io->base && end - start < 10) {
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
    cache_incr(io, c);
    //contig_dump_ps(io, &c, "/tmp/tree.ps");
    r = contig_bins_in_range(io, &c, start, end,
			     CSIR_SORT_BY_X | CSIR_LEAVES_ONLY,
			     CONS_BIN_SIZE, &nr);


    /*
     * When doing range queries, we can end up with a higher bin being
     * labelled as a leaf from contig_bins_in_range simply because we
     * pruned a child off. Eg:
     *
     *    |...|
     * [------------------]   <- Leaf
     *    |...|  [--------]
     *
     * Yet when doing full queries, the bottom node is the leaf instead.
     * This causes caching issues as we vary which nodes we'll return as
     * suitable for storing consensus in.
     *
     * We resolve this by filtering bins >= CONS_BIN_SIZE that also
     * have at least one child node.
     */
    for (i = j = 0; i < nr; i++) {
	bin_index_t *bin = cache_search(io, GT_Bin, r[i].rec);
	
	if (bin->size >= CONS_BIN_SIZE && (bin->child[0] || bin->child[1]))
	    continue;

	if (j != i) r[j] = r[i];
	j++;
    }
    nr = j;

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
	cache_incr(io, bin);
	    
	/* NB: The filtering above can cause holes too */
	if (r[i].start > left) {
	    gio_debug(io, 1, "Filling missing region %d..%d\n", left, right);
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
	    seq_t *s = NULL, *dup_s = NULL, seq;
	    range_t *cons_r = NULL;
	    int n, valid;
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

	    valid = 0;
	    if (cons_r && (bin->flags & BIN_CONS_VALID)) {
		//bstart = cons_r->start;
		//bend   = cons_r->end;
		s = (seq_t *)cache_search(io, GT_Seq, cons_r->rec);

		/* Double check sizes */
		if (s && cons_r->start == 0 && cons_r->end == bin->size-1) {
		    cache_incr(io, s);
		    valid = 1;
		} else {
		    /*
		     * This can happen when we have overlapping bins.
		     * A former call here could compute the consensus,
		     * but bin_add_range can put it in another choice.
		     * When this happens the consensus may not span the
		     * entire bin, in which case we'll recompute it
		     * again here.
		     *
		     * This has proven more reliable than attempting
		     * to coerce bin_add_range to add to specific
		     * bins.
		     */
		    //fprintf(stderr, "Cached consensus rec #%"PRIrec
		    //	    " appears to be invalid\n", s->rec);
		    valid = 0;
		    if (!s)
			cons_r = NULL;
		}
	    }
	    if (!valid) {
		/* Not cached, or is cached but needs updating */
		range_t r2, *r_out;
		tg_rec recno;
		int comp;
		float *tmp_qual;

		//bstart = MAX(r[i].start, c->start);
		//bend   = MIN(r[i].end-1, c->end);

		/* Load existing sequence, or fill out in-memory new struct */
		if (cons_r && !io->read_only) {
		    gio_debug(io, 1,
			      "Recomputing cached cons from %d..%d +%d "
			      "Bin #%"PRIrec"\n",
			      cons_r->start, cons_r->end, r[i].start,
			      bin->rec);

		    s = (seq_t *)cache_search(io, GT_Seq, cons_r->rec);
		    if (!s)
			goto load_failed;
		    s = cache_rw(io, s);

		    if (cons_r->start != 0 || cons_r->end != bin->size-1) {
			size_t extra_len =
			    (s->name       ? strlen(s->name)       : 0)+1 +
			    (s->trace_name ? strlen(s->trace_name) : 0)+1 +
			    (s->alignment  ? strlen(s->alignment)  : 0)+1 +
			    (bend - bstart + 1)*2;

			/* Reset positions/lengths */
			s->pos   = bstart;
			s->len   = (s->len > 0)
			    ?   bend - bstart + 1
			    : -(bend - bstart + 1);
			s = cache_item_resize(s, sizeof(*s) + extra_len);
			sequence_reset_ptr(s);

			bin = cache_rw(io, bin);
			bin->flags |= BIN_BIN_UPDATED | BIN_RANGE_UPDATED;
			cons_r->start = 0;
			cons_r->end = bin->size-1;
			bin->start_used = 0;
			bin->end_used = bin->size-1;
		    }

		    if ((s->len < 0) ^ r[i].comp) {
		    	s->len = -s->len;
			if (s->len < 0)
			    s->flags |= SEQ_COMPLEMENTED;
			else
			    s->flags &= ~SEQ_COMPLEMENTED;
		    }
		    
		} else {
		load_failed: /* Horrid, I know! */

		    /*
		     * If read-only we can't store this anyway, so just
		     * abort and compute the consensus we asked for instead
		     * of taking time to cache data that we immediately
		     * have to throw away.
		     */
		    if (io->read_only || end - start < 1000) {
			cache_decr(io, c);
			cache_decr(io, bin);
			if (s && s != &seq)
			    cache_decr(io, s);
			if (r)
			    free(r);
			return calculate_consensus_simple2(io, contig,
							   start, end,
							   con, qual);
		    }


		    gio_debug(io, 1, "Creating new cached cons %d..%d "
			      "bin #%"PRIrec"\n",
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
		if (!qual && io->read_only) {
		    tmp_qual = NULL;
		} else {
		    tmp_qual = calloc(bend - bstart + 1, sizeof(float));
		}

		calculate_consensus_simple2(io, contig,
					    bstart, bend,
					    s->seq, tmp_qual);

		if (tmp_qual) {
		    for (n = 0; n <= bend-bstart && n < ABS(s->len); n++) {
			int q = rint(tmp_qual[n]);
			if (q < 0)
			    s->conf[n] = 0;
			else if (q > 127)
			    s->conf[n] = 127;
			else
			    s->conf[n] = q;
		    }			
		    free(tmp_qual);
		}


		/* If cached and range changed, move sequence */
		if (cons_r) {
		    //assert(cons_r->start == bstart);
		    //assert(cons_r->end   == bend);
		    // if (!io->read_only) { ... }
		}
		
		/* If not cached, create a new range */
		if (!io->read_only && !cons_r && !io->base) {
		    r2.start    = bstart;
		    r2.end      = bend;
		    r2.rec      = 0;
		    r2.pair_rec = 0;
		    r2.mqual    = 0;
		    r2.flags    = GRANGE_FLAG_TYPE_SINGLE;
		    r2.flags   |= GRANGE_FLAG_ISCONS;
		    r2.library_rec = 0;
			
		    cache_decr(io, bin);
		    bin = bin_add_to_range(io, &c, bin->rec, &r2, &r_out,
					   &comp, 0);
		    cache_incr(io, bin);

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
		    complement_seq_t(s);
		}
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

	cache_decr(io, bin);
	left = r[i].end;
    }
    if (r) free(r);

    cache_decr(io, c);

    if (end >= left) {
	//printf("Filling missing region %d..%d\n", left, end);
	calculate_consensus_simple2(io, contig, left, end,
				    con  ? &con[left-start]  : NULL,
				    qual ? &qual[left-start] : NULL); 
    }

    if (!io->base && !io->read_only)
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
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, contig);
    cache_incr(io, c);

    /* Compute in small ranges */
    for (i = start; i <= end; i += CONS_BLOCK_SIZE) {
	rangec_t *r = NULL;
	int nr;
	int st = i;
	int en = st + CONS_BLOCK_SIZE-1;
	if (en > end)
	    en = end;

	/* Find sequences visible */
	r = contig_seqs_in_range(io, &c, st, en, 0, &nr);

	if (0 != calculate_consensus_bit_het(io, contig, st, en,
					     CONS_ALL, r, nr,
					     &cons[i-start])) {
	    if (r)
		free(r);
	    cache_decr(io, c);
	    return -1;
	}

	if (r)
	    free(r);
    }

    cache_decr(io, c);
    return 0;
}

/* Just compute the basic consensus + phred score */
int calculate_consensus_fast(GapIO *io, tg_rec contig, int start, int end,
			     consensus_t *cons) {
    int i;
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, contig);

    /* Compute in small ranges */
    for (i = start; i <= end; i += CONS_BLOCK_SIZE) {
	rangec_t *r = NULL;
	int nr;
	int st = i;
	int en = st + CONS_BLOCK_SIZE-1;
	if (en > end)
	    en = end;

	/* Find sequences visible */
	r = contig_seqs_in_range(io, &c, st, en, 0, &nr);

	if (0 != calculate_consensus_bit_het(io, contig, st, en,
					     0, r, nr, &cons[i-start])) {
	    return -1;
	    if (r)
		free(r);
	}

	if (r)
	    free(r);
    }

    return 0;
}

typedef struct {
    double base_qual[4];
    double gap;
    double base;
    int depth;
} cstat;


#if 0
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
#endif


#define P_HET 1e-6

static double prior[25];     /* Sum to 1.0 */
static double lprior15[15];  /* 15 combinations of {ACGT*} */
static double acc_prior[25]; /* Cumulative values up to 32768 */

/* Precomputed matrices for the consensus algorithm */
static double pMM[101], p__[101], p_M[101];

static double e_tab_a[1002];
static double *e_tab = &e_tab_a[500];
static double e_log[501];

static void consensus_init(double p_het) {
    int i;
    double acc = 0;

    memset(lookup, 5, 256*sizeof(lookup[0]));
    lookup['A'] = lookup['a'] = 0;
    lookup['C'] = lookup['c'] = 1;
    lookup['G'] = lookup['g'] = 2;
    lookup['T'] = lookup['t'] = 3;
    lookup['*'] =               4;

    for (i = -500; i <= 500; i++)
    	e_tab[i] = exp(i);
    for (i = 0; i <= 500; i++)
	e_log[i] = log(i);

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
 * See "A Fast, Compact Approximation of the Exponential Function"
 * by NN. Schraudolph, Neural Computation, 1999
 */
#if 0
static inline double fast_exp(double y) {
    union {
	double d;
	int i, j;
    } x;

    x.i = 0;
    x.j = 1512775 * y + 1072632447;

    return x.d;
}

static inline double fast_log(double y) {
    union {
	double d;
	int i, j;
    } x;
    
    x.d = y;
    return (x.j - 1072632447.0) / 1512775;
}
#endif

#if 1
static inline double fast_exp(double y) {
    if (y < -500)
	y = -500;
    if (y > 500)
	y = 500;

    //    printf("%f => %g %g\n", y, exp(y), e_tab[(int)y]);

    return e_tab[(int)y];
}
#endif

//#define fast_exp exp
#define fast_log log


/*
 * As above, but allowing for the possibility of the data being
 * heterozygous.
 */
int calculate_consensus_bit_het(GapIO *io, tg_rec contig,
				int start, int end, int flags,
				rangec_t *r, int nr,
				consensus_t *cons) {
    int i, j;
    int len = end - start + 1;
    static int loop = 0;
    static int init_done =0;
    static double q2p[101];
    double min_e_exp = DBL_MIN_EXP * log(2) + 1;

    double (*scores)[15];
    double (*sumsC)[6] = NULL, *sumsE = NULL;
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

    /* Currently always needed for the N vs non-N evaluation */
    flags |= CONS_COUNTS;

    /* Initialise */
    if (NULL == (scores = (double (*)[15])calloc(len, 15 * sizeof(double))))
	return -1;
    if (flags & CONS_DISCREP) {
	if (NULL == (sumsC = (double (*)[6])calloc(len, 6 * sizeof(double))))
	    return -1;
	if (NULL == (sumsE = (double *)calloc(len, sizeof(double))))
	    return -1;
    }

    if (NULL == (depth = (int *)calloc(len, sizeof(int))))
	return -1;
    if (NULL == (perfect = (char *)calloc(len, sizeof(char))))
	return -1;

    if (flags & CONS_COUNTS) {
	for (i = 0; i < len; i++) {
	    cons[i].counts[0] = 0;
	    cons[i].counts[1] = 0;
	    cons[i].counts[2] = 0;
	    cons[i].counts[3] = 0;
	    cons[i].counts[4] = 0;
	    cons[i].counts[5] = 0;
	}
    }

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

	    if (s->format != SEQ_FORMAT_CNF4) {
		/* Inline version for speed */
		base = s->seq[j+off];
		qual = s->conf[j+off];
	    } else {
		sequence_get_base(io, &s, j+off, &base, &qual, NULL, 0);
	    }

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

	    if (flags & CONS_DISCREP) {
		qe = q2p[qual];
		sumsE[sp-start+j] += qe;
		sumsC[sp-start+j][base_l] += 1 - qe;
	    }

	    if (flags & CONS_COUNTS)
		cons[sp-start+j].counts[base_l]++;


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
    if ((flags & CONS_NO_END_N) && depth[0] == 0 && loop == 0) {
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

    if ((flags & CONS_NO_END_N) && depth[len-1] == 0 && loop == 0) {
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
	int call = 0, het_call = 0, ph;
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
		//S[j] = exp(S[j]);
		S[j] = fast_exp(S[j]);
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

	/* And store result */
//	if (depth[i] &&
//	    !(fabs(norm[0] - norm[5]) < FLT_EPSILON  &&
//	      fabs(norm[0] - norm[9]) < FLT_EPSILON  &&
//	      fabs(norm[0] - norm[12]) < FLT_EPSILON &&
//	      fabs(norm[0] - norm[14]) < FLT_EPSILON)) {

	if (depth[i] && depth[i] != cons[i].counts[5] /* all N */) {
	    double m;

	    cons[i].depth = depth[i];

	    cons[i].call     = map_sing[call];
	    if (norm[call] == 0) norm[call] = DBL_MIN;
	    ph = -TENOVERLOG10 * fast_log(norm[call]) + .5;
	    cons[i].phred = ph > 255 ? 255 : (ph < 0 ? 0 : ph);

	    cons[i].het_call = map_het[het_call];
	    if (norm[het_call] == 0) norm[het_call] = DBL_MIN;
	    ph = TENOVERLOG10 * (fast_log(S[het_call]) - fast_log(norm[het_call])) + .5;
	    cons[i].scores[6] = ph;

	    if (flags & CONS_SCORES) {
		/* AA */
		if (norm[0] == 0) norm[0] = DBL_MIN;
		ph = TENOVERLOG10 * (fast_log(S[0]) - fast_log(norm[0])) + .5;
		cons[i].scores[0] = ph;

		/* CC */
		if (norm[5] == 0) norm[5] = DBL_MIN;
		ph = TENOVERLOG10 * (fast_log(S[5]) - fast_log(norm[5])) + .5;
		cons[i].scores[1] = ph;

		/* GG */
		if (norm[9] == 0) norm[9] = DBL_MIN;
		ph = TENOVERLOG10 * (fast_log(S[9]) - fast_log(norm[9])) + .5;
		cons[i].scores[2] = ph;

		/* TT */
		if (norm[12] == 0) norm[12] = DBL_MIN;
		ph = TENOVERLOG10 * (fast_log(S[12]) - fast_log(norm[12])) + .5;
		cons[i].scores[3] = ph;

		/* ** */
		if (norm[14] == 0) norm[14] = DBL_MIN;
		ph = TENOVERLOG10 * (fast_log(S[14]) - fast_log(norm[14])) + .5;
		cons[i].scores[4] = ph;

		/* N */
		cons[i].scores[5] = 0; /* N */
	    }

	    /* Compute discrepancy score */
	    if (flags & CONS_DISCREP) {
		m = sumsC[i][0]+sumsC[i][1]+sumsC[i][2]+sumsC[i][3]+sumsC[i][4]
		    - sumsC[i][cons[i].call];
		cons[i].discrep = (m-sumsE[i])/sqrt(m+sumsE[i]);
	    }
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

    free(scores);
    free(depth);
    free(perfect);
    if (sumsE)
	free(sumsE);
    if (sumsC)
	free(sumsC);

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

	ci = contig_iter_new(io, contig, 1, CITER_FIRST | CITER_ISTART |
			     CITER_SMALL_BS, CITER_CSTART, CITER_CEND);
	while (ci && (r = contig_iter_next(io, ci))) {
	    seq_t *s;
	    int left;

	    if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
		continue;

	    if (r->start > best)
		break;

	    s = cache_search(io, GT_Seq, r->rec);
	    if (!s) {
		verror(ERR_WARN, "consensus_valid_range",
		       "Failed to load seq #%"PRIrec, r->rec);
		continue;
	    }

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

	ci = contig_iter_new(io, contig, 1, CITER_LAST | CITER_IEND |
			     CITER_SMALL_BS, CITER_CSTART, CITER_CEND);
	
	while (ci && (r = contig_iter_prev(io, ci))) {
	    seq_t *s;
	    int right;

	    if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
		continue;

	    if (r->end < best)
		break;

	    s = cache_search(io, GT_Seq, r->rec);
	    if (!s) {
		verror(ERR_WARN, "consensus_valid_range",
		       "Failed to load seq #%"PRIrec, r->rec);
		continue;
	    }

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

/*
 * As per consensus_valid_range(), but we return the unclipped size of the
 * contig instead. This is somewhat easier and faster to compute.
 * Passing NULL as start or end implies you are uninterested in that end.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
int consensus_unclipped_range(GapIO *io, tg_rec contig, int *start, int *end) {
    contig_iterator *ci;
    rangec_t *r;

    if (start) {
	int best = INT_MAX;

	ci = contig_iter_new(io, contig, 1, CITER_FIRST | CITER_ISTART |
			     CITER_SMALL_BS, CITER_CSTART, CITER_CEND);

	while (ci && (r = contig_iter_next(io, ci))) {
	    if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
		continue;

	    if (r->start > best)
		break;

	    if (best > r->start)
		best = r->start;
	}

	contig_iter_del(ci);
	*start = best != INT_MAX ? best : 0;
    }

    if (end) {
	int best = INT_MIN;

	ci = contig_iter_new(io, contig, 1, CITER_LAST | CITER_IEND |
			     CITER_SMALL_BS, CITER_CSTART, CITER_CEND);
	
	while (ci && (r = contig_iter_prev(io, ci))) {
	    if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
		continue;

	    if (r->end < best)
		break;

	    if (best < r->end)
		best = r->end;
	}

	contig_iter_del(ci);
	*end = best != INT_MIN ? best : 0;
    }

    return 0;
}

/*
 * Converts a padded position into an unpadded position.
 * Returns 0 for success and writes to upos
 *        -1 for error
 */
int consensus_unpadded_pos(GapIO *io, tg_rec contig, int pos, int *upos) {
    int np, i, i_end, clipped_start;
    char *cons;
    contig_t *c;

    consensus_valid_range(io, contig, &clipped_start, NULL);

    /* For now it's the slow route: compute consensus and count '*'s */
    if (NULL == (c = cache_search(io, GT_Contig, contig)))
	return TCL_ERROR;
    if (pos <= c->start) {
	*upos = pos - clipped_start + 1;
	return 0;
    }

    if (NULL == (cons = malloc(pos - c->start + 1)))
	return -1;

    if (-1 == calculate_consensus_simple(io, contig, clipped_start,
					 pos, cons, NULL)) {
	free(cons);
	return -1;
    }

    i_end = pos - clipped_start;
    for (np = 0, i = 0; i < i_end; i++) {
	if (cons[i] == '*')
	    np++;
    }

    *upos = pos - np - clipped_start + 1;
    free(cons);

    return 0;
}

/*
 * Converts an unpadded position into a padded position.
 * Returns 0 for success and writes to pos
 *        -1 for error
 */
int consensus_padded_pos(GapIO *io, tg_rec contig, int upos, int *pos) {
    int np, i, i_end, offset, clipped_start;
    char *cons;
    contig_t *c;

    consensus_valid_range(io, contig, &clipped_start, NULL);

    /* For now it's the slow route: compute consensus and count '*'s */
    if (NULL == (c = cache_search(io, GT_Contig, contig)))
	return TCL_ERROR;
    if (upos <= 0) {
	*pos = upos + clipped_start-1;
	return 0;
    }

    /* First guess at 8k more */
    if (NULL == (cons = malloc(upos + 1 + 8192)))
	return -1;

    if (-1 == calculate_consensus_simple(io, contig,
					 clipped_start,
					 clipped_start + upos + 8192,
					 cons, NULL)) {
	free(cons);
	return -1;
    }

    i_end = upos;
    offset = clipped_start;
    np = 0;
    for(;;) {
	int sz;

	for (i = 0; i < i_end; i++) {

	    if (cons[i] == '*')
		np++;
	    if (i + offset - (clipped_start-1) >= upos + np)
		break;
	}
	
	if (i + offset - (clipped_start-1) >= upos + np)
	    break;

	/* If we haven't found the unpadded location yet, keep going */
	sz = MAX(upos + np - (i + offset), 8192);
	if (-1 == calculate_consensus_simple(io, contig,
					     i + offset, i + offset + sz,
					     cons, NULL)) {
	    free(cons);
	    return -1;
	}
					     
	offset += i;
	i_end = sz+1;
    }

    *pos = i + offset;
    free(cons);

    return 0;
}

#define WLEN 12
#define WSIZE (1<<(2*WLEN))
static unsigned char uhash[WSIZE];

static char *s(uint32_t w) {
    static char s[WLEN+1];
    int i;

    for (i = 0; i < WLEN; i++) {
	s[WLEN-i-1] = "ACGT"[w & 3];
	w >>= 2;
    }
    s[WLEN] = 0;

    return s;
}


int update_uniqueness_hash(GapIO *io) {
    int i, j, cnum;
    char *cons = NULL;
    size_t cons_l = 0;
    uint32_t w, W;

    consensus_init(P_HET);

    memset(uhash, 0, WSIZE);

    for (cnum = 0; cnum < io->db->Ncontigs; cnum++) {
	tg_rec crec = arr(tg_rec, io->contig_order, cnum);
	contig_t *c = cache_search(io, GT_Contig, crec);
	size_t len;

	if (!c) {
	    if (cons)
		xfree(cons);
	    return -1;
	}

	len = c->end - c->start + 1;
	if (cons_l < len) {
	    cons_l = len;
	    if (NULL == (cons = xrealloc(cons, cons_l)))
		return -1;
	}

	//printf("%d/%d len %d\n", cnum+1, io->db->Ncontigs, len);

	calculate_consensus_simple(io, crec, c->start, c->end, cons, NULL);
	
	for (W = w = i = j = 0; i < len && j < WLEN-1; i++) {
	    int c = lookup[cons[i]];
	    if (c < 4) {
		w = (w << 2) | c;
		W = (W << 2) | (c ^ 3);
		j++;
	    }
	}
	for (; i < len; i++) {
	    int c = lookup[cons[i]];
	    if (c < 4) {
		w = (w << 2) | c;
		w &= WSIZE-1;

		if (uhash[w] < 255)
		    uhash[w]++;

		W = (W << 2) | (c^3);
		W &= WSIZE-1;

		if (uhash[W] < 255)
		    uhash[W]++;

		//printf("%s %d\n", s(w), uhash[w]);
	    }
	}
    }

//    for (i = 0; i < WSIZE; i++) {
//	if (uhash[i])
//	    printf("%d %d\n", i, uhash[i]);
//    }

    return 0;
}

int get_uniqueness(GapIO *io, tg_rec crec, int pos) {
    char cons[WLEN*4+1];
    int i, j, w = 0;
    char *cp;

    cp = &cons[WLEN*2];
    calculate_consensus_simple(io, crec, pos-WLEN*2, pos+WLEN*2, cons, NULL);

    for (i = 0, j = WLEN/2; i >= -WLEN*2 && j > 0; i--) {
	if (lookup[cp[i] < 4])
	    j--;
    }

    if (j != 0)
	return 0;

    for (; i < WLEN*2 && j < WLEN; i++) {
	int c = lookup[cp[i]];
	if (c < 4) {
	    w = (w << 2) | c;
	    j++;
	}
    }

    printf("%s %d\n", s(w), uhash[w]);

    return uhash[w];
}

int get_uniqueness_pos(char *str, int len, int pos) {
    char cons[WLEN*4+1];
    int i, j, w = 0;

    for (i = pos, j = WLEN/2; i > 0 && j > 0; i--) {
	if (lookup[str[i]] < 4)
	    j--;
    }

    if (j != 0)
	return 0;

    for (; i < len && j < WLEN; i++) {
	int c = lookup[str[i]];
	if (c < 4) {
	    w = (w << 2) | c;
	    j++;
	}
    }

    printf("%s %d\n", s(w), uhash[w]);

    return uhash[w];
}
