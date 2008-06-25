/*
 * Computes consensus and consensus quality.
 *
 *
 * Find highest 4 chems/strands for each base type.
 * Pick out highest bases of each different chemistry, strand, base type.
 * Then perform a full bayesian check for each base type to see which is the
 * most likely.
 *
 * The method of combining the 'dependent' results is as follows:
 *
 * OLD METHOD:
 * For dependent observations ob1, ob2, .., obN
 * we take average ob "(ob1 + ob2 + ... + obN) / N" and multiply by
 * dependent_table[N], where dependent_table[] can be written to
 * give an indication of how much weighting to apply to additional results.
 * Tried using sqrt(N) for this - too optimistic. Now have a function that
 * approaches 2 as N increases, so that we never more than halve the error
 * rate.
 *
 * NEW METHOD:
 * For dependent observations ob1, ob2, .., obN
 * we find highest ob (obi) as a probability and divide this by 0.95^(N-1)
 * Ie, we take the observation with the highest probability, but we consider
 * combine with other dependent observations at a fixed .95 error rate.
 * This provides a small increase when dealing with many sequences.
 */



/*
 * Authors:
 *   jkb - James Bonfield
 *
 *    Apr 1993 jkb - created
 * 27 May 1993 jkb - change qual_cutoff to use -1 for old algorithm.
 * 01 Mar 1994 jkb - merged into gap
 * 23 Aug 1994 jkb - better (fixed size) memory usage & speed tuneups
 * 03 Feb 1999 johnt - fixed strchr overflow in qual_char
 * 14 Mar 2000 jkb - Split IO components from algorithm components. to allow
 *                   compilation of qual.c outside of Gap4.
 */
/* #define DEBUG_QUAL */


/* ----- #includes ----- */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <sys/types.h>
#include <math.h>
#include <float.h>
#include "qual.h"
#include "qualP.h"
#include "misc.h"
#include "xalloc.h"

/* ----- Standalone compilation, when used outside of Gap4 ----- */
#ifdef QUAL_STANDALONE

float consensus_cutoff;
int quality_cutoff = -1;
int chem_as_double = 1;
int consensus_mode = CONSENSUS_MODE_CONFIDENCE;

#else
#include "gap_globals.h"
#endif


/* ----- typedefs ----- */

typedef char StringBuffer[200];

/* ----- global variables ----- */
char q_lookup[3][3][2] = {
    {
	{R_NONE_NONE,	R_NONE_NONE},
	{R_NONE_GOOD,	R_NONE_GOOD},
	{R_NONE_BAD,	R_NONE_BAD}
    },
    {
	{R_GOOD_NONE,	R_GOOD_NONE},
	{R_GOOD_GOOD_NE,R_GOOD_GOOD_EQ},
	{R_GOOD_BAD,	R_GOOD_BAD}
    },
    {
	{R_BAD_NONE,	R_BAD_NONE},
	{R_BAD_GOOD,	R_BAD_GOOD},
	{R_BAD_BAD,	R_BAD_BAD}
    }
};

int consensus_iub = 0;

/* ----- `local' global variables ----- */

/*
 * qual_cutoff_def is the default quality cutoff, as specified by the
 * "set display parameters" option. 
 *
 * qual_cutoff_tmp is the value used internally by all calculations. This
 * is set to qual_cutoff_def if we pass QUAL_DEFAULT to a routine, otherwise
 * the qual_cutoff value passed to the routine.
 */
static int qual_cutoff_def = -1;
static int qual_cutoff_tmp = -1;

/* ----- local ----- */

struct con_data_t {
    char *con1, *con2;
    float *qual1, *qual2;
    float cons_cutoff;
};

struct qual_data {
    char *qual;
    float cons_cutoff;
};

typedef struct {
    int (*qual)[2]; /* base type (0-5) and base conf */
    int start;
    int end;
    int gel;
    int dir;
    int chem;
    int template;
} seq_frag;

static signed int clookup[256];
static int clookup_old[256];
static int vlookup_old[256];

static int calc_contig_info(int contig, int start, int end, int both,
			    void (*eval_func)(int (*c_qual1)[7],
					      int (*c_qual2)[7],
					      int   len,
					      int   pos,
					      void *con_data),
			    void *eval_data,
			    int (*info_func)(int         job,
					     void       *mydata,
					     info_arg_t *theirdata),
			    void *info_data);

static void quality_func(int (*c_qual1)[7], int (*c_qual2)[7], int len,
			 int pos, void *qual_data);

static void consensus_func(int (*c_qual1)[7], int (*c_qual2)[7], int len,
			   int pos, void *con_data);

static void init_clookup(void);

static int calc_contig_info_phred(int contig, int start, int end,
				  char *con, float *qual,
				  char *con2, float *qual2,
				  float cons_cutoff, int qual_cutoff,
				  void (*process_frags)(seq_frag *frag,
							int   *num_frags,
							int    from,
							int    to,
							int    start,
							char  *con1,
							float *qual1,
							char  *con2,
							float *qual2,
							float  cons_cutoff,
							int    qual_cutoff),
				  int (*info_func)(int         job,
						   void       *mydata,
						   info_arg_t *theirdata),
				  void *info_data);

static void process_frags(seq_frag *frag, int *num_frags, int from, int to,
			  int start, char *con1, float *qual1,
			  char *con2, float *qual2,
			  float cons_cutoff, int qual_cutoff);

/* ----- functions ----- */

/*
 * Our 'callback' function from calc_contig_info() when calculating the
 * quality.
 * See calc_contig_info() for comments on how and when this is called.
 */
static void quality_func(int (*c_qual1)[7], int (*c_qual2)[7], int len,
			 int pos, void *qual_data_v) {
    struct qual_data *qual_data = qual_data_v;
    register int j, t;
    register int *conq, *conq2;
    unsigned char base_p = 0, base_n = 0, qual_p, qual_n;
    char *qual;
    float cons_cutoff;

    qual = &qual_data->qual[pos];
    cons_cutoff = qual_data->cons_cutoff;

    /*
     * Collate quality data together to create quality types
     */
    for (j = 0; j < len; j++) {
	conq  = c_qual1[j];
	conq2 = c_qual2[j];
	
	if (conq[5] == -1) {
	    qual_p = Q_NO_DATA;
	} else {
	    /* check for 100% quality bases */
	    if (conq[6]) {
		switch (conq[6]) {
		case 1<<0: /* A */
		    base_p = 0;
		    break;
		case 1<<1: /* C */
		    base_p = 1;
		    break;
		case 1<<2: /* G */
		    base_p = 2;
		    break;
		case 1<<3: /* T */
		    base_p = 3;
		    break;
		case 1<<4: /* * */
		    base_p = 4;
		    break;
		default:
		    base_p = 5;
		}
		qual_p = base_p == 5 ? Q_POOR_DATA : Q_OK;

	    } else {
		/* simple check for best */
		t = conq[0]; base_p = 0;
		if (t < conq[1]) t = conq[1], base_p = 1;
		if (t < conq[2]) t = conq[2], base_p = 2;
		if (t < conq[3]) t = conq[3], base_p = 3;
		if (t < conq[4]) t = conq[4], base_p = 4;
		
		/* check whether data is good, bad or non existant */
#ifdef DEBUG_QUAL
		printf("base_p = %d, conq[base_p] = %3d, conq[5] = %3d\n",
		       base_p, conq[base_p], conq[5]);
#endif
		if (base_p == 4 /* pad */) {
		    qual_p = Q_POOR_DATA;
		} else {
		    if (conq[5]) {
			qual_p = Q_OK;
			if ((float)conq[base_p]/conq[5] - cons_cutoff
			    < -FLT_EPSILON)
			    qual_p = Q_POOR_DATA;
		    } else
			qual_p = Q_POOR_DATA;
		}
	    }
	}

	if (conq2[5] == -1) {
	    qual_n = Q_NO_DATA;
	} else {
	    if (conq2[6]) {
		switch (conq2[6]) {
		case 1<<0: /* A */
		    base_n = 0;
		    break;
		case 1<<1: /* C */
		    base_n = 1;
		    break;
		case 1<<2: /* G */
		    base_n = 2;
		    break;
		case 1<<3: /* T */
		    base_n = 3;
		    break;
		case 1<<4: /* * */
		    base_n = 4;
		    break;
		default:
		    base_n = 5;
		}
		qual_n = base_n == 5 ? Q_POOR_DATA : Q_OK;

	    } else {
		t = conq2[0]; base_n = 0;
		if (t < conq2[1]) t = conq2[1], base_n = 1;
		if (t < conq2[2]) t = conq2[2], base_n = 2;
		if (t < conq2[3]) t = conq2[3], base_n = 3;
		if (t < conq2[4]) t = conq2[4], base_n = 4;

#ifdef DEBUG_QUAL
		printf("base_n = %d, conq2[base_n] = %3d, conq2[5] = %3d\n",
		       base_n, conq2[base_n], conq2[5]);
#endif
		if (base_n == 4 /* pad */) {
		    qual_n = Q_POOR_DATA;
		} else {
		    if (conq2[5]) {
			qual_n = Q_OK;
			if ((float)conq2[base_n]/conq2[5] - cons_cutoff
			    < -FLT_EPSILON)
			    qual_n = Q_POOR_DATA;
		    } else
			qual_n = Q_POOR_DATA;
		}
	    }
	}

	/* now generate quality code based upon above measures */
	*qual++ = q_lookup[qual_p][qual_n][base_p == base_n];
#ifdef DEBUG_QUAL
	printf("%4d %s, %s + %s\n",
	       i, reasons2[qual[i]],
	       reasons[qual_p], reasons[qual_n]);
#endif
    }
}

int calc_quality(int   contig,
		 int   start,
		 int   end,
		 char *qual,
		 float cons_cutoff,
		 int   qual_cutoff,
		 int (*info_func)(int        job,
				  void       *mydata,
				  info_arg_t *theirdata),
		 void *info_data)
{
    struct qual_data qual_d;

    init_clookup();

    qual_d.qual = qual;
    qual_d.cons_cutoff = cons_cutoff;

    qual_cutoff_tmp = (qual_cutoff == QUAL_DEFAULT)
	? qual_cutoff_def : qual_cutoff;

    if (-1 == calc_contig_info(contig, start, end, 1,
			       quality_func, (void *)&qual_d,
			       info_func, info_data)) {
	return -1;
    }

    return 0;
}

/**/

/*
 * Scan along from 'position' until we find a single stranded region, or hit
 * the end of the sequence. The returned value is the new position, 0 if no
 * problem found, -1 for error.
 *
 * With contig != 0 the necessary data is initialised. At this stage, position
 * and rreg mark the left and right ends of the range of based to be scanned.
 * Setting these to 0 represents all. For the initialisation, reason and len
 * are ignored.
 *
 * Subsequent calls should be made with contig == 0 until a new contig is
 * to be used. In these cases cons_cutoff, info_func, rreg and info_data are
 * ignored. The value of rreg supplied during the initialisation is used
 * instead as this provides higher robustness.
 *
 * NOTE: This does not find the two end cases as 'holes' when both strands end
 * at the same base.
 */
int next_hole(int contig,
	      int position,
	      int rreg,
	      float cons_cutoff,
	      int   qual_cutoff,
	      char **reason,
	      int *len,
	      int (*info_func)(int        job,
			       void       *mydata,
			       info_arg_t *theirdata),
	      void *info_data)
{
    info_arg_t info;
    register int i;
    static char *quality = NULL;
    static int ll;
    static int rr;
    int iend;
    
    qual_cutoff_tmp = (qual_cutoff == QUAL_DEFAULT)
	? qual_cutoff_def : qual_cutoff;

    /* init */
    if (contig) {
	info.contig_info.contig = contig;
	info_func(GET_CONTIG_INFO, info_data, &info);

	rr = rreg ? rreg : info.contig_info.length;
	ll = position ? position : 1;

	if (quality)
	    xfree(quality);

	if (NULL == (quality = (char *)xmalloc(rr - ll + 1))) {
	    return -1;
	}

	if (calc_quality(contig, ll, rr, quality, cons_cutoff, QUAL_DEFAULT,
			 info_func, info_data) == -1) {
	    verror(ERR_FATAL, "next_hole", "Failed to calculate quality");
	    return -1;
	}

	return 0;
    }
    
    iend = rr - ll;
    for (i=position-ll; i<=iend; i++) {
	register int q = quality[i];

	/*
	 * Perhaps could be sped up by changing this if to use a precomputed
	 * array of which types are considered problems.
	 */
	if (q == R_GOOD_NONE || q == R_BAD_NONE ||
	    q == R_NONE_GOOD || q == R_NONE_BAD ||
	    q == R_NONE_NONE) {
	    int j = i + 1;

	    *reason = &quality[i];
	    
	    /*
	     * Find extent of problem. Make sure we don't scan past into
	     * another type of problem.
	     */
	    if (q == R_GOOD_NONE || q == R_BAD_NONE) {
		for (; j<=rr-ll; j++) {
		    q = quality[j];
		    if (q != R_GOOD_NONE && q != R_BAD_NONE)
			break;
		}
	    } else if (q == R_NONE_GOOD || q == R_NONE_BAD) {
		for (; j<=rr-ll; j++) {
		    q = quality[j];
		    if (q != R_NONE_GOOD && q != R_NONE_BAD)
			break;
		}
	    } else if (q == R_NONE_NONE) {
		for (; j<=rr-ll; j++) {
		    q = quality[j];
		    if (q != R_NONE_NONE)
			break;
		}
	    }

	    *len = j-i;
	    return i + ll;
	}
    }
    
    return 0;
}

/* Anything not in here (eg IUB codes) are treated as dashes */
static char qual_char[]={
    'C','T','A','G',
    'c','t','a','g',
    '*',',','-','\0'};

static char qual_ind[sizeof(qual_char)]={
    1,3,0,2,
    1,3,0,2,
    4,4,5,5};

/*
 * Note that for figures with qual_ind == 5 they will be added twice (once to
 * that qual_ind and once to index 5 as the running total). So we halve these
 * values.
 */
#define N_SCORE 5
static char qual_val[sizeof(qual_char)]={
     99, 99, 99, 99,
     99, 99, 99, 99,
     99, 99,  N_SCORE, N_SCORE};

static void init_clookup(void) {
    int i;
    static int done = 0;

    if (done)
	return;
    done = 1;

    /* Unknowns are -1 as we treat them specially (add part to all) */
    for (i=0; i<256; i++) {
	clookup[i] = -1;
	clookup_old[i] = N_SCORE;
	vlookup_old[i] = N_SCORE; /* => 2*N_SCORE as it's added twice */
    }
    for (i=0; i<9; i++)
	clookup[(unsigned char)"ACGT*acgt"[i]] = i%5;

    for (i=0; i<256; i++) {
	int c;
	char *cp;

	cp = strchr(qual_char, i);
	c = cp ? cp - qual_char : 10 /*'-'*/;
	clookup_old[i] = qual_ind[c];
	vlookup_old[i] = qual_val[c];
    }
}

/*
 * Creates and returns a `length' by 2 deep array of integers representing
 * the called sequence and it's quality.
 */
static int (*get_gel_qual(int gel, int start, int end,
			  int (*info_func)(int         job,
					   void       *mydata,
					   info_arg_t *theirdata),
			  void *info_data))[2]
{
    const char *gel_seq;
    const int1 *gel_conf;
    register int (*gel_conf2)[2];
    register int i, i_end;
    info_arg_t info;

    /* Give the C optimisers a helping hand */
    const int *clookup_t = clookup;

    /*
     * Get the sequence data
     */
    info.gel_seq.gel = gel;
    if (-1 == info_func(GET_SEQ, info_data, &info)) {
	verror(ERR_FATAL, "get_gel_qual",
	       "Failed to read sequence for gel no. %d\n", gel);
	return (int (*)[2])-1;
    }

    i_end = end - start;

    gel_conf2 = (int (*)[2])xmalloc(i_end * sizeof(int) * 2);
    gel_conf  = &info.gel_seq.gel_conf[start + info.gel_seq.gel_start];
    gel_seq   = &info.gel_seq.gel_seq[start + info.gel_seq.gel_start];

    if (consensus_mode == CONSENSUS_MODE_WEIGHTED ||
	consensus_mode == CONSENSUS_MODE_CONFIDENCE) {
	for (i=0; i < i_end; i++) {
	    gel_conf2[i][0] = clookup_t[(unsigned char)gel_seq[i]];
	    gel_conf2[i][1] = gel_conf[i];
	}
    } else {
	for (i=0; i < i_end; i++) {
	    int c = gel_seq[i];
	    gel_conf2[i][0] = clookup_old[c];
	    gel_conf2[i][1] = vlookup_old[c];
	}
    }

    /*
     * Tidy up
     */
    (void)info_func(DEL_SEQ, info_data, &info);
    
    return gel_conf2;
}

/*
 * Generates quality for a contig over a given region. If `both' is set then
 * the data is stored in two separate arrays (one for each strand), otherwise
 * it is merged into one (the first in the structure).
 *
 * In the interests of keeping memory usage to a minimum we declare our
 * 7 deep consensus array to be only twice the maximum gel length long.
 * Once we've filled one half of this we call the eval_func() to process
 * this buffer. We assume then that this function has done whatever it needs
 * to do, and will proceed to overwrite this data upon the next 'cycle' of
 * our cyclic buffer.
 *
 * eval_func() should take the following arguments:
 *
 * int (*c_qual1)[7]	The 7 deep quality array for +ve strand
 * int (*c_qual2)[7]	The 7 deep quality array for -ve strand (or NULL)
 * int len		The length of data in c_qual1 and c_qual2
 * int pos		The position in the start..end range this data is for
 * void *con_data	Private data - it's the "eval_data" argument here.
 *
 * Returns -1 for failure, 0 for success.
 */
static int calc_contig_info(int contig, int start, int end, int both,
			    void (*eval_func)(int (*c_qual1)[7],
					      int (*c_qual2)[7],
					      int len,
					      int pos,
					      void *con_data),
			    void *eval_data,
			    int (*info_func)(int         job,
					     void       *mydata,
					     info_arg_t *theirdata),
			    void *info_data)
{
    register int (*gel_qual)[2];	/* gel quality */
    register int (*contig_qual)[7];	/* contig quality info */
    int (*contig_qual2)[7];		/* contig quality info */
    int i_start, comp, not_comp;
    register int i, j, i_end, g;
    info_arg_t info;
    int max_len, max_len2, buffer_pos, buffer_ind;

    /*
     * Quick sanity checks
     */
    if (start > end)
	return -1;

    /*
     * Get contig information & allocate appropriate data.
     * We allocate a max of 2*max_gel_length
     */
    max_len = info_func(GET_GEL_LEN, info_data, &info);
    if (max_len > end - start + 1)
	max_len = end - start + 1;
    max_len2 = max_len * 2;

    info.contig_info.contig = contig;
    info_func(GET_CONTIG_INFO, info_data, &info);

    contig_qual = (int (*)[7])xmalloc(max_len2 * 7 * sizeof(int));

    /*
     * FIXME: Rewrite in a portable fashion
     *
     * (char)-1 is 0xff in 2's complement notation.
     * (int)-1 is similarly 0xffffffff (for a 4 byte int).
     *
     * memset() writes characters at a time. Hence the following ASSUMES that
     * (int)-1 is equivalen to four consecutive (char)-1's.
     */
    memset(contig_qual, -1, max_len2 * 7 * sizeof(int));

    if (both) {
	contig_qual2 = (int (*)[7])xmalloc(max_len2 * 7 * sizeof(int));
	memset(contig_qual2, -1, max_len2 * 7 * sizeof(int));
    } else {
	contig_qual2 = NULL;
    }

    
    /*
     * find left most gel in our region
     */
    info.gel_info.gel = info.contig_info.leftgel;
    do {
	info_func(GET_GEL_INFO, info_data, &info);
    } while (info.gel_info.position + info.gel_info.length < start &&
	     (info.gel_info.gel = info.gel_info.next_right));
    
    buffer_pos = start-1;
    buffer_ind = 0;

    /*
     * Yes I know this is hideous. It fixes a bug where given illegal data
     * (contigs with holes in it) we will want to run the 'eval_func' on
     * our empty buffers until we've got to the region where the next
     * sequence is visible, otherwise the first sequence (aka
     * info.gel_info.gel) may get shunted down to appear in an incorrect
     * portion of the consensus.
     *
     * Ie we want to do B AB AB AB AB AB and so we goto the B section to start
     * the ball rolling. Alternatives are simply to duplicate the code or to
     * put it in its own function, although it's have rather a lot of
     * pass-by-reference parameters.
     *
     * I also toyed with Duff's Device, but this made the code even more
     * obscure. At least this way it's self documenting :-)
     */
    goto process_buffers;

    do {
	if (info.gel_info.gel == 0)
	    break;

	i_start = info.gel_info.position < start ?
	    start - info.gel_info.position : 0;
	i_end = info.gel_info.position + info.gel_info.length > end ?
	    end - info.gel_info.position + 1: info.gel_info.length;

	if (i_end <= i_start)
	    goto fetch_next; /* ok - so it's a cheat! */

	/*
	 * Get quality info for gel
	 */
	
	if ((int (*)[2])-1 ==
	    (gel_qual = get_gel_qual(info.gel_info.gel,
				     i_start, i_end,
				     info_func, info_data))) {
	    xfree(contig_qual);
	    return -1;
	}

	/*
	 * Update running totals
	 */
	if (both && chem_as_double && info.gel_info.as_double) {
	    comp = 1;
	    not_comp = 1;
	} else {
	    comp = info.gel_info.complemented && contig_qual2;
	    not_comp = !comp;
	}
	i_end -= i_start;

	if (comp) {
	    for (i = 0, j = i_start + info.gel_info.position - start;
		 i < i_end; i++, j++) {
		register int j2 = j % max_len2;
		/*
		 * flag as used - initialise to 0 if currently -1
		 */
                if (contig_qual2[j2][0] == -1) {
		    int *speedup = contig_qual2[j2];
		    speedup[0] = speedup[1] = speedup[2] = 0;
		    speedup[3] = speedup[4] = speedup[5] = 0;
		    contig_qual2[j2][6] = 0;
                }

		if (gel_qual[i][0] != -1) {
		    if ((g=gel_qual[i][1]) >= qual_cutoff_tmp) {
			if (gel_qual[i][1] != 100) {
			    contig_qual2[j2][gel_qual[i][0]] += g;
			} else {
			    contig_qual2[j2][6] |= 1<<(gel_qual[i][0]);
			}
		    } else {
			g = 0;
		    }
		} else {
		    /* A '-', so we add 1/4 of it to all base types */
		    g = gel_qual[i][1];
		    contig_qual2[j2][0] += g/4;
		    contig_qual2[j2][1] += g/4;
		    contig_qual2[j2][2] += g/4;
		    contig_qual2[j2][3] += g/4;
		}
                
                contig_qual2[j2][5] += g;
            }
	}

	if (not_comp) {
	    for (i = 0, j = i_start + info.gel_info.position - start;
		 i < i_end; i++, j++) {
		register int j2 = j % max_len2;

                if (contig_qual[j2][0] == -1) {
		    int *speedup = contig_qual[j2];
		    speedup[0] = speedup[1] = speedup[2] = 0;
		    speedup[3] = speedup[4] = speedup[5] = 0;
		    contig_qual[j2][6] = 0;
                }

		if (gel_qual[i][0] != -1) {
		    if ((g=gel_qual[i][1]) >= qual_cutoff_tmp) {
			if (gel_qual[i][1] != 100) {
			    contig_qual[j2][gel_qual[i][0]] += g;
			} else {
			    contig_qual[j2][6] |= 1<<(gel_qual[i][0]);
			}
		    } else {
			g = 0;
		    }
		} else {
		    /* A '-', so we add 1/4 of it to all base types */
		    g = gel_qual[i][1];
		    contig_qual[j2][0] += g/4;
		    contig_qual[j2][1] += g/4;
		    contig_qual[j2][2] += g/4;
		    contig_qual[j2][3] += g/4;
		}

		contig_qual[j2][5] += g;
            }
	}

	xfree(gel_qual);

    fetch_next:
	/*
	 * fetch next gel
	 */
	info.gel_info.gel = info.gel_info.next_right;
	if (info.gel_info.gel)
	    info_func(GET_GEL_INFO, info_data, &info);

    process_buffers:
	/*
	 * If we've filled up one buffer load, then call our evaluation
	 * function.
	 */
	while (info.gel_info.position > buffer_pos + max_len) {
	    eval_func(&contig_qual[buffer_ind],
		      contig_qual2 ? &contig_qual2[buffer_ind] : NULL,
		      end - buffer_pos > max_len ? max_len : end - buffer_pos,
		      buffer_pos - (start-1), eval_data);

	    /* clear cyclic buffer section */
	    memset(&contig_qual[buffer_ind], -1, max_len * 7*sizeof(int));
	    if (contig_qual2)
		memset(&contig_qual2[buffer_ind], -1, max_len * 7*sizeof(int));

	    /* switch to opposite end of cyclic buffer */
	    buffer_ind ^= max_len;
	    buffer_pos += max_len;
	}

	/* loop until past end */

    } while (info.gel_info.position <= end);

    while (end - buffer_pos > 0) {
	eval_func(&contig_qual[buffer_ind],
		  contig_qual2 ? &contig_qual2[buffer_ind] : NULL,
		  end - buffer_pos > max_len ? max_len : end - buffer_pos,
		  buffer_pos - (start-1), eval_data);

	/*
	 * No need to clear previously used buffer as this can at most
	 * loop back once.
	 */

	buffer_ind ^= max_len;
	buffer_pos += max_len;
    }

    if (contig_qual)
	xfree(contig_qual);
    if (contig_qual2)
	xfree(contig_qual2);

    return 0;
}

/**/

/*
 * Our 'callback' function from calc_contig_info() when calculating the
 * consensus.
 * See calc_contig_info() for comments on how and when this is called.
 */
static void consensus_func(int (*c_qual1)[7], int (*c_qual2)[7], int len,
			   int pos, void *con_data_v) {
    struct con_data_t *con_data = con_data_v;
    register int j;
    register int *conq;
    int (*c_qual)[7], counts;
    float val, *qual, *discrep;
    char *con, *con1, *con2;
    unsigned char base;
    float *qual1, *qual2;
    float cons_cutoff;

    con1 = con_data->con1;
    con2 = con_data->con2;
    qual1 = con_data->qual1;
    qual2 = con_data->qual2;
    cons_cutoff = con_data->cons_cutoff;

    c_qual = c_qual1;
    con = &con1[pos];
    qual = qual1 ? &qual1[pos] : NULL;
    discrep = qual2 ? &qual2[pos] : NULL;

    counts = con2 ? 2 : 1;

    for (; counts; counts--) {
	for (j = 0; j < len; j++) {
	    register int t;

	    conq = c_qual[j];
	    
	    /*
	     * Check first for absolute 'forced' base type.
	     * Only allow one such type, otherwise it's a dash.
	     */
	    if (conq[6]) {
		switch(conq[6]) {
		case 1<<0: /* A */
		    base = 0;
		    break;
		case 1<<1: /* C */
		    base = 1;
		    break;
		case 1<<2: /* G */
		    base = 2;
		    break;
		case 1<<3: /* T */
		    base = 3;
		    break;
		case 1<<4: /* * */
		    base = 4;
		    break;
		default:
		    base = 5;
		    break;
		}
		val = base == 5 ? 0 : 100;

	    } else {
		if (conq[5] > 0) {
		    if (consensus_iub) {
			t = 0;
			if ((float)conq[0]/conq[5] - cons_cutoff
			    >= -FLT_EPSILON)
			    t |= 1;
			if ((float)conq[1]/conq[5] - cons_cutoff
			    >= -FLT_EPSILON)
			    t |= 2;
			if ((float)conq[2]/conq[5] - cons_cutoff
			    >= -FLT_EPSILON)
			    t |= 4;
			if ((float)conq[3]/conq[5] - cons_cutoff
			    >= -FLT_EPSILON)
			    t |= 8;
			base = 6 + t;
			val = 0;
		    } else {
			/* Check for best base - defaults to '-' */
			t = 0; base = 5; /* default to '-' */
			if (t < conq[0]) t = conq[0], base = 0;
			if (t < conq[1]) t = conq[1], base = 1;
			if (t < conq[2]) t = conq[2], base = 2;
			if (t < conq[3]) t = conq[3], base = 3;
			if (t < conq[4]) t = conq[4], base = 4;

			if ((val = (float)conq[base]/conq[5]) - cons_cutoff
			    < -FLT_EPSILON)
			    base = 5;
		    }
		} else {
		    val = 0;
		    base = 5;
		}
	    }

	    /* Obtaining 2nd-highest confidence */
	    if (!con2 && discrep) {
		int first = 0, second = 0;
		int k;
		for (k = 0; k < 5; k++) {
		    if (first < conq[k]) {
			second = first;
			first = conq[k];
		    } else if (second < conq[k]) {
			second = conq[k];
		    }
		}

		*discrep++ = (second * 100.0) / conq[5];
	    }
	    
	    /*
	     * We could remove this 'if' by assigning qualp to a dummy location
	     * when qualp is not set. Tmp memory vs speed tradeoff...
	     */
	    if (qual)
		*qual++ = val * 100.0;
	    
	    *con++ = "ACGT*-NACMGRSVTWYHKDBN"[base];
	}

	/*
	 * Go around for second spin if doing 2 stranded consensus
	 */
	c_qual = c_qual2;
	con = con2 ? &con2[pos] : NULL;
	qual = qual2 ? &qual2[pos] : NULL;
    }
}

/*
 * If con2 is supplied then consensus is returned in two buffers - one per
 * strand. Otherwise (con2 == NULL), consensus is returned in con
 * If qual (+ qual2) is supplied then the consensus calculations are stored in
 * that array. (Ie it is possible to see the real consensus cutoff value for
 * that base.)
 */
int calc_consensus(int   contig,
		   int   start,
		   int   end,
		   int   mode,
		   char *con,
		   char *con2,
		   float *qual,
		   float *qual2,
		   float cons_cutoff,
		   int   qual_cutoff,
		   int (*info_func)(int        job,
				    void       *mydata,
				    info_arg_t *theirdata),
		   void *info_data)
{
    register int i;
    char *conp = con, *conp2;
    int con1_extra = 0, con2_extra = 0;
    struct con_data_t con_d;

    init_clookup();

    qual_cutoff_tmp = (qual_cutoff == QUAL_DEFAULT)
	? qual_cutoff_def : qual_cutoff;

    if (CONSENSUS_MODE_CONFIDENCE == consensus_mode) {
	if (-1 ==calc_contig_info_phred(contig, start, end,
					con, qual, con2, qual2,
					cons_cutoff, qual_cutoff_tmp,
					process_frags,
					info_func, info_data)) {
	    return -1;
	}

	return 0;
    }

    /* Old algorithms - CONSENSUS_MODE_FREQ and CONSENSUS_MODE_WEIGHTED */

    /*
     * for mode == CON_WDET we allocate the extra array ourselves and
     * deallocate at the end
     */
    if (mode == CON_WDET) {
	if (NULL == (con2 = (char *)xmalloc(end - start + 1)))
	    return -1;
	con2_extra = 1;
    }

    con_d.con1 = con;
    con_d.con2 = con2;
    con_d.qual1 = qual;
    con_d.qual2 = qual2;
    con_d.cons_cutoff = cons_cutoff;

    if (-1 ==calc_contig_info(contig, start, end, con2 != NULL,
			      consensus_func, (void *)&con_d,
			      info_func, info_data)) {
	return -1;
    }

    /*
     * For the Well DETermined mode, we find places where the sequence is
     * well determined and matching on both strands, and then translate to
     * new codes. This is currently used only by the mask/marking code.
     */
    if (mode == CON_WDET) {
        unsigned char wdet_tab[256];

	/* FIXME: use Rodger's lookup table */
	for (i = 0; i < 256; i++)
	    wdet_tab[i] = i;
	wdet_tab['A'] = 'd';
	wdet_tab['C'] = 'e';
	wdet_tab['G'] = 'f';
	wdet_tab['T'] = 'i';

	conp = con;
	conp2 = con2;

	for (i = start; i <= end; i++, conp++, conp2++) {
	    if (*conp == *conp2)
		*conp = wdet_tab[(unsigned)*conp];
	    else {
		/* - x => x
		 * x y => - (where x != y && y != '-')
		 */
		if (*conp == '-')
		    *conp = *conp2;
		else if (*conp2 != '-' && *conp != *conp2)
		    *conp = '-';
	    }
/*	    printf(" %c\n",*conp); */
	}
    }

    if (con1_extra)
	xfree(con);
    if (con2_extra)
	xfree(con2);

    return 0;
}

/*
 * Set quality cutoff and return old value.
 * Should only be used as a temporary measure as the user prefers to feel in
 * control of such things themselves.
 */
int set_qual_cutoff(int new) {
    int tmp = qual_cutoff_def;

    qual_cutoff_def = new;
    return tmp;
}

/*
 * Query quality cutoff.
 */
int query_qual_cutoff(void) {
    return qual_cutoff_def;
}


/*
 * ---------------------------------------------------------------------------
 * New consensus calculation code using phred values
 * ---------------------------------------------------------------------------
 */

/*
 * A lookup table to determine how to combine multiple sequences that are
 * considered as dependent.
 * We taken the average quality value multiplied by dependent_table[count],
 * where count is the number of dependent observations.
 */

/* sqrt(N)
double dependent_table[11] =
{
    0.0,
    1.0,
    1.41421356237309504880,
    1.73205080756887729352,
    1.0,
    2.23606797749978969640,
    2.44948974278317809819,
    2.64575131106459059050,
    2.82842712474619009760,
    3.00000000000000000000,
    3.16227766016837933199,
};
*/

/* Maxes out as 2.0 */
/*
double dependent_table[11] =
{
    0.0,
    1.0,
    1.19,
    1.38,
    1.56,
    1.71,
    1.82,
    1.90,
    1.95,
    1.97,
    1.99,
};
*/

double dependent_table[11] = {
    /*  0 + 0 */ 0.0,
    /*  1 + 0 */ 0.0,
    /*  2 + 7 */ 7.0,
    /*  3 + 7 */ 14.0, 
    /*  4 + 7 */ 21.0, 
    /*  5 + 7 */ 28.0,
    /*  6 + 7 */ 35.0,
    /*  7 + 6 */ 41.0, 
    /*  8 + 4 */ 45.0,
    /*  9 + 3 */ 48.0,
    /* 10 + 1 */ 49.0,
};

double depth_scale[11] = {
    /*  0 */ 0.0,
    /*  1 */ 0.80,
    /*  2 */ 0.82,
    /*  3 */ 0.95,
    /*  4 */ 1.00,
    /*  5 */ 1.10,
    /*  6 */ 1.40,
    /*  7 */ 1.65,
    /*  8 */ 2.30,
    /*  9 */ 3.00,
    /* 10 */ 5.00
};


/*
 * A lookup table to adjust phred quality values, based on observed vs
 * expected error rates. Sampled over 14 million bases (748000 wrong).
 */
#if 0
double phred_table[] = {
    /*  0 */ 0.0000000000,
    /*  1 */ 1.0000000000,
    /*  2 */ 2.0000000000,
    /*  3 */ 3.0000000000,
    /*  4 */ 4.7588794280,
    /*  5 */ 4.6416944109,
    /*  6 */ 6.7862289566,
    /*  7 */ 7.8514577864,
    /*  8 */ 8.7907197211,
    /*  9 */ 9.4701814169,
    /* 10 */ 10.4856598779,
    /* 11 */ 11.3552237192,
    /* 12 */ 12.7612409416,
    /* 13 */ 13.7014510718,
    /* 14 */ 15.0789544175,
    /* 15 */ 16.2869136420,
    /* 16 */ 16.4517081419,
    /* 17 */ 17.8018184512,
    /* 18 */ 19.0439030456,
    /* 19 */ 19.8176078219,
    /* 20 */ 20.5708315475,
    /* 21 */ 22.3504349071,
    /* 22 */ 22.6180987466,
    /* 23 */ 24.4481565202,
    /* 24 */ 25.2954375420,
    /* 25 */ 26.0083089978,
    /* 26 */ 27.0780375368,
    /* 27 */ 28.0625536579,
    /* 28 */ 29.6375460197,
    /* 29 */ 28.7360151362,
    /* 30 */ 32.1936506907,
    /* 31 */ 32.7157501633,
    /* 32 */ 32.5052776175,
    /* 33 */ 34.1930396664,
    /* 34 */ 34.8013785119,
    /* 35 */ 31.1872519033,
    /* 36 */ 35.6909918078,
    /* 37 */ 38.4051131031,
    /* 38 */ 37.3462604393,
    /* 39 */ 37.1743518925,
    /* 40 */ 40.0000000000,
    /* 41 */ 41.0000000000,
    /* 42 */ 42.0000000000,
    /* 43 */ 43.0000000000,
    /* 44 */ 44.0000000000,
    /* 45 */ 45.0000000000,
    /* 46 */ 46.0000000000,
    /* 47 */ 47.0000000000,
    /* 48 */ 48.0000000000,
    /* 49 */ 49.0000000000,
    /* 50 */ 50.0000000000,
    /* 51 */ 51.0000000000,
    /* 52 */ 52.0000000000,
    /* 53 */ 53.0000000000,
    /* 54 */ 54.0000000000,
    /* 55 */ 55.0000000000,
    /* 56 */ 56.0000000000,
    /* 57 */ 57.0000000000,
    /* 58 */ 58.0000000000,
    /* 59 */ 59.0000000000,
    /* 60 */ 60.0000000000,
    /* 61 */ 61.0000000000,
    /* 62 */ 62.0000000000,
    /* 63 */ 63.0000000000,
    /* 64 */ 64.0000000000,
    /* 65 */ 65.0000000000,
    /* 66 */ 66.0000000000,
    /* 67 */ 67.0000000000,
    /* 68 */ 68.0000000000,
    /* 69 */ 69.0000000000,
    /* 70 */ 70.0000000000,
    /* 71 */ 71.0000000000,
    /* 72 */ 72.0000000000,
    /* 73 */ 73.0000000000,
    /* 74 */ 74.0000000000,
    /* 75 */ 75.0000000000,
    /* 76 */ 76.0000000000,
    /* 77 */ 77.0000000000,
    /* 78 */ 78.0000000000,
    /* 79 */ 79.0000000000,
    /* 80 */ 80.0000000000,
    /* 81 */ 81.0000000000,
    /* 82 */ 82.0000000000,
    /* 83 */ 83.0000000000,
    /* 84 */ 84.0000000000,
    /* 85 */ 85.0000000000,
    /* 86 */ 86.0000000000,
    /* 87 */ 87.0000000000,
    /* 88 */ 88.0000000000,
    /* 89 */ 89.0000000000,
    /* 90 */ 90.0000000000,
    /* 91 */ 91.0000000000,
    /* 92 */ 92.0000000000,
    /* 93 */ 93.0000000000,
    /* 94 */ 94.0000000000,
    /* 95 */ 95.0000000000,
    /* 96 */ 96.0000000000,
    /* 97 */ 96.9999999999,
    /* 98 */ 97.9999999999,
    /* 99 */ 99.0000000000,
    /* 100 */ 100.0000000000
};
#endif

/* NUL phred table - has no change */
double phred_table[] = {
    /*  0 */ 0,
    /*  1 */ 1,
    /*  2 */ 2,
    /*  3 */ 3,
    /*  4 */ 4,
    /*  5 */ 5,
    /*  6 */ 6,
    /*  7 */ 7,
    /*  8 */ 8,
    /*  9 */ 9,
    /* 10 */ 10,
    /* 11 */ 11,
    /* 12 */ 12,
    /* 13 */ 13,
    /* 14 */ 14,
    /* 15 */ 15,
    /* 16 */ 16,
    /* 17 */ 17,
    /* 18 */ 18,
    /* 19 */ 19,
    /* 20 */ 20,
    /* 21 */ 21,
    /* 22 */ 22,
    /* 23 */ 23,
    /* 24 */ 24,
    /* 25 */ 25,
    /* 26 */ 26,
    /* 27 */ 27,
    /* 28 */ 28,
    /* 29 */ 29,
    /* 30 */ 30,
    /* 31 */ 31,
    /* 32 */ 32,
    /* 33 */ 33,
    /* 34 */ 34,
    /* 35 */ 35,
    /* 36 */ 36,
    /* 37 */ 37,
    /* 38 */ 38,
    /* 39 */ 39,
    /* 40 */ 40,
    /* 41 */ 41,
    /* 42 */ 42,
    /* 43 */ 43,
    /* 44 */ 44,
    /* 45 */ 45,
    /* 46 */ 46,
    /* 47 */ 47,
    /* 48 */ 48,
    /* 49 */ 49,
    /* 50 */ 50,
    /* 51 */ 51,
    /* 52 */ 52,
    /* 53 */ 53,
    /* 54 */ 54,
    /* 55 */ 55,
    /* 56 */ 56,
    /* 57 */ 57,
    /* 58 */ 58,
    /* 59 */ 59,
    /* 60 */ 60,
    /* 61 */ 61,
    /* 62 */ 62,
    /* 63 */ 63,
    /* 64 */ 64,
    /* 65 */ 65,
    /* 66 */ 66,
    /* 67 */ 67,
    /* 68 */ 68,
    /* 69 */ 69,
    /* 70 */ 70,
    /* 71 */ 71,
    /* 72 */ 72,
    /* 73 */ 73,
    /* 74 */ 74,
    /* 75 */ 75,
    /* 76 */ 76,
    /* 77 */ 77,
    /* 78 */ 78,
    /* 79 */ 79,
    /* 80 */ 80,
    /* 81 */ 81,
    /* 82 */ 82,
    /* 83 */ 83,
    /* 84 */ 84,
    /* 85 */ 85,
    /* 86 */ 86,
    /* 87 */ 87,
    /* 88 */ 88,
    /* 89 */ 89,
    /* 90 */ 90,
    /* 91 */ 91,
    /* 92 */ 92,
    /* 93 */ 93,
    /* 94 */ 94,
    /* 95 */ 95,
    /* 96 */ 96,
    /* 97 */ 97,
    /* 98 */ 98,
    /* 99 */ 99,
    /* 100 */ 100
};

#define DEPEND_ADD 7

/*
 * Processes blocks of overlapping reads to produce the consensus making
 * use of the phred quality values.
 *
 * This method partitions data into sets defined by chemistry and strand
 * Ie that A(+term) and T(+term) are in the same set, but A(+term) and A(-term)
 * are in different sets.
 * Then within a set we split up into the four (five with pad) base types
 * and use Bayesian rules to calculate the relative probabilities of each
 * base type occuring within that set.
 * For multiple base types within the same set (eg A(+term) and A(+term))
 * we take the best four. This gives at most 4*5 (A,C,G,T,*) events for
 * the bayesian1 compution to compute the relative probalities for this
 * chemistry/strand set.
 * Then we use Bayesian rules again on the four possible chem/strand sets
 * to determine the final sequence.
 *
 * NB: THIS IS A WASTE OF CPU TIME!
 * It is equivalent to picking the best four from each and doing one bayesian
 * computation treating them all as independent.
 * The same applies to the plus DEPEND_ADD method below too.
 */
#if 0
static void process_frags_6(seq_frag *frag, int *num_frags, int from, int to,
			  int start, char *con1, float *qual1,
			  float cons_cutoff, int qual_cutoff)
{
    int i, j, k, l, nevents1, nevents2;
    int nf = *num_frags;
    int qual, type;
    int highest_type;
    double highest_product;
    int perfect;
    char *con = &con1[from-start];
    float *qualp = NULL;
    int qhighest[6][4][4]; /* highest 4 of each type and strand/chem */
    int chem;
    double prob, err, qnorm, sproduct;
    double bayesian1[5][5*4], bayesian2[5][5], product[5];
    int pad_present, nbase_types;
    int cons_cutoff100 = (int)(cons_cutoff * 100 + 0.001);

    if (qual1)
	qualp = &qual1[from-start];

    for (i = from; i < to; i++) {
	/* Initialise total/count arrays */
	memset(qhighest, 0, 6*4*4*sizeof(int));

	/* Fill qhighest arrays. Also remove fragments as needed */
	perfect = 0;
	pad_present = 0;
	for (j = 0; j < nf; j++) {
	    type = frag[j].qual[frag[j].start][0];
	    if (type == -1)
		type = 5;
	    /* Type: A=0, C=1, G=2, T=3, *=4, -=5 */
	    qual = frag[j].qual[frag[j].start][1];
	    chem = frag[j].dir + 2 * frag[j].chem;
	    if (qual == 100)
		perfect |= 1<<type;
	    if (qual > qhighest[type][chem][0]) {
		qhighest[type][chem][3] = qhighest[type][chem][2];
		qhighest[type][chem][2] = qhighest[type][chem][1];
		qhighest[type][chem][1] = qhighest[type][chem][0];
		qhighest[type][chem][0] = qual;
	    } else if (qual > qhighest[type][chem][1]) {
		qhighest[type][chem][3] = qhighest[type][chem][2];
		qhighest[type][chem][2] = qhighest[type][chem][1];
		qhighest[type][chem][1] = qual;
	    } else if (qual > qhighest[type][chem][2]) {
		qhighest[type][chem][3] = qhighest[type][chem][2];
		qhighest[type][chem][2] = qual;
	    } else if (qual > qhighest[type][chem][3]) {
		qhighest[type][chem][3] = qual;
	    }
	    if (type == 4)
		pad_present = 1;
	    
	    if (++frag[j].start >= frag[j].end) {
		xfree(frag[j].qual);
		memmove(&frag[j], &frag[j+1], (nf-j-1) * sizeof(*frag));
		nf--; j--;
	    }
	}

	/* Handle perfect bases (qual == 100) */
	if (perfect) {
	    char base;

	    switch(perfect) {
	    case 1<<0: base = 'A'; break;
	    case 1<<1: base = 'C'; break;
	    case 1<<2: base = 'G'; break;
	    case 1<<3: base = 'T'; break;
	    case 1<<4: base = '*'; break;
	    default:   base = '-'; break;
	    }

	    *con++ = base;
	    if (qualp)
		*qualp++ = (base == '-') ? 0 : 100;
	    
	    continue; /* to "for (i = from; i < to; i++)" loop */
	}

	/*
	 * Having found our highest values we now fill out the bayesian
	 * matrices.
	 *
	 * Ignore dash for now as we don't know what to do with it.
	 */
	nevents2 = 0;
	if (!pad_present) {
	    double add[4];
	    double tmp;
	    
	    nbase_types = 4;
	    for (k = 0; k < 4; k++) { /* strand+chem */
		for (j = 0; j < 4; j++)
		    add[j] = 0.0;

		/* bayesian1 calc on this strand+chem */
		nevents1 = 0;
		for (j = 0; j < 4; j++) { /* type */
		    for (l = 0; l < 4; l++) {
			if (qhighest[j][k][l]) {
			    tmp = (double)phred_table[qhighest[j][k][l]];
			    prob = 1 - pow(10.0, -tmp / 10.0);
			    bayesian1[0][nevents1] = (1 - prob) / 3;
			    bayesian1[1][nevents1] = (1 - prob) / 3;
			    bayesian1[2][nevents1] = (1 - prob) / 3;
			    bayesian1[3][nevents1] = (1 - prob) / 3;
			    bayesian1[j][nevents1] = prob;
			    nevents1++;
			}
		    }
		}

		/* Compute bayesian1 denominator */
		qnorm = 0;
		for (j = 0; j < 4; j++) {
		    product[j] = 1;
		    for (l = 0; l < nevents1; l++) {
			product[j] *= bayesian1[j][l];
		    }
		    qnorm += product[j];
		}

		/* Fill out bayesian2 arrays */
		if (qnorm) {
		    for (j = 0; j < 4; j++) {
			bayesian2[j][nevents2] = product[j] / qnorm;
		    }
		    nevents2++;
		}
	    }
	} else {
	    double add[5];
	    double tmp;
	    
	    nbase_types = 5;
	    for (k = 0; k < 4; k++) { /* strand+chem */
		for (j = 0; j < 5; j++)
		    add[j] = 0.0;

		/* bayesian1 calc on this strand+chem */
		nevents1 = 0;
		for (j = 0; j < 5; j++) { /* type */
		    for (l = 0; l < 4; l++) {
			if (qhighest[j][k][l]) {
			    tmp = (double)phred_table[qhighest[j][k][l]];
			    prob = 1 - pow(10.0, -tmp / 10.0);
			    bayesian1[0][nevents1] = (1 - prob) / 4;
			    bayesian1[1][nevents1] = (1 - prob) / 4;
			    bayesian1[2][nevents1] = (1 - prob) / 4;
			    bayesian1[3][nevents1] = (1 - prob) / 4;
			    bayesian1[4][nevents1] = (1 - prob) / 4;
			    bayesian1[j][nevents1] = prob;
			    nevents1++;
			}
		    }
		}

		/* Compute bayesian1 denominator */
		qnorm = 0;
		for (j = 0; j < 5; j++) {
		    product[j] = 1;
		    for (l = 0; l < nevents1; l++) {
			product[j] *= bayesian1[j][l];
		    }
		    qnorm += product[j];
		}

		/* Fill out bayesian2 arrays */
		if (qnorm) {
		    for (j = 0; j < 5; j++) {
			bayesian2[j][nevents2] = product[j] / qnorm;
		    }
		    nevents2++;
		}
	    }
	}

	/* Compute bayesian1 denominator */
	qnorm = 0;
	highest_product = 0;
	highest_type = 5;
	for (j = 0; j < nbase_types; j++) {
	    sproduct = 1;
	    for (k = 0; k < nevents2; k++) {
		sproduct *= bayesian2[j][k];
	    }
	    qnorm += sproduct;

	    if (sproduct > highest_product) {
		highest_product = sproduct;
		highest_type = j;
	    }
	}
	    
	/* Called base has highest prob. - normalise it to get new prob. */
	prob = qnorm ? highest_product / qnorm : 0;
	if (prob < 1.0)
	    err = -10.0 * log10(1-prob);
	else
	    err = 999;

	/* Now store the results in con and qual */
	if (err - cons_cutoff100 >= -FLT_EPSILON)
	    *con++ = "ACGT*-"[highest_type];
	else
	    *con++ = '-';
	if (qualp)
	    *qualp++ = err > 99 ? 99 : err;
    }

    fflush(stdout);
    *num_frags = nf;
}
#endif

/*
 * Processes blocks of overlapping reads to produce the consensus making
 * use of the phred quality values.
 *
 * This method partitions data into sets defined by chemistry and strand
 * Ie that A(+term) and T(+term) are in the same set, but A(+term) and A(-term)
 * are in different sets.
 * Then within a set we split up into the four (five with pad) base types
 * and use Bayesian rules to calculate the relative probabilities of each
 * base type occuring within that set.
 * For multiple base types within the same set (eg A(+term) and A(+term))
 * we take the highest confidence value plus DEPEND_ADD multiplied by the
 * number of other bases in this set with this base type.
 * Then we use Bayesian rules again on the four possible sets to determine
 * the final sequence.
 *
 * NB: THIS IS A WASTE OF CPU TIME!
 * It is equivalent to method 3 (below).
 * However it may provide a better basis to distinguish combining same
 * strand/chem vs differing strand/chem.
 */
#if 0
static void process_frags_5(seq_frag *frag, int *num_frags, int from, int to,
			  int start, char *con1, float *qual1,
			  float cons_cutoff, int qual_cutoff)
{
    int i, j, k, l, nevents1, nevents2;
    int nf = *num_frags;
    int qual, type;
    int highest_type;
    double highest_product;
    int perfect;
    char *con = &con1[from-start];
    float *qualp = NULL;
    int qhighest[6][4];
    int qcount[6][4];
    int chem;
    double prob, err, qnorm, sproduct;
    double bayesian1[5][5], bayesian2[5][5], product[5];
    int pad_present, nbase_types;
    int cons_cutoff100 = (int)(cons_cutoff * 100 + 0.001);

    if (qual1)
	qualp = &qual1[from-start];

    for (i = from; i < to; i++) {
	/* Initialise total/count arrays */
	memset(qhighest, 0, 6*4*sizeof(int));
	memset(qcount, 0, 6*4*sizeof(int));

	/* Fill qhighest and qcount arrays. Also remove fragments as needed */
	perfect = 0;
	pad_present = 0;
	for (j = 0; j < nf; j++) {
	    type = frag[j].qual[frag[j].start][0];
	    if (type == -1)
		type = 5;
	    /* Type: A=0, C=1, G=2, T=3, *=4, -=5 */
	    qual = frag[j].qual[frag[j].start][1];
	    chem = frag[j].dir + 2 * frag[j].chem;
	    if (qual == 100)
		perfect |= 1<<type;
	    if (qual > qhighest[type][chem])
		qhighest[type][chem] = qual;
	    qcount[type][chem]++;
	    if (type == 4)
		pad_present = 1;
	    
	    if (++frag[j].start >= frag[j].end) {
		xfree(frag[j].qual);
		memmove(&frag[j], &frag[j+1], (nf-j-1) * sizeof(*frag));
		nf--; j--;
	    }
	}

	/* Handle perfect bases (qual == 100) */
	if (perfect) {
	    char base;

	    switch(perfect) {
	    case 1<<0: base = 'A'; break;
	    case 1<<1: base = 'C'; break;
	    case 1<<2: base = 'G'; break;
	    case 1<<3: base = 'T'; break;
	    case 1<<4: base = '*'; break;
	    default:   base = '-'; break;
	    }

	    *con++ = base;
	    if (qualp)
		*qualp++ = (base == '-') ? 0 : 100;
	    
	    continue; /* to "for (i = from; i < to; i++)" loop */
	}

	/*
	 * Having found our highest values we now fill out the bayesian
	 * matrices.
	 *
	 * Ignore dash for now as we don't know what to do with it.
	 */
	nevents2 = 0;
	if (!pad_present) {
	    double add[4];
	    double tmp;
	    
	    nbase_types = 4;
	    for (k = 0; k < 4; k++) { /* strand+chem */
		for (j = 0; j < 4; j++)
		    add[j] = 0.0;

		/* bayesian1 calc on this strand+chem */
		nevents1 = 0;
		for (j = 0; j < 4; j++) { /* type */
		    if (qhighest[j][k]) {
			/* tmp = (double)phred_table[qhighest[j][k]];*/
			tmp = (double)phred_table[qhighest[j][k]] +
			    (qcount[j][k] - 1) * DEPEND_ADD;
			prob = 1 - pow(10.0, -tmp / 10.0);
			bayesian1[0][nevents1] = (1 - prob) / 3;
			bayesian1[1][nevents1] = (1 - prob) / 3;
			bayesian1[2][nevents1] = (1 - prob) / 3;
			bayesian1[3][nevents1] = (1 - prob) / 3;
			bayesian1[j][nevents1] = prob;
			nevents1++;
		    }
		}

		/* Compute bayesian1 denominator */
		qnorm = 0;
		for (j = 0; j < 4; j++) {
		    product[j] = 1;
		    for (l = 0; l < nevents1; l++) {
			product[j] *= bayesian1[j][l];
		    }
		    qnorm += product[j];
		}

		/* Fill out bayesian2 arrays */
		if (qnorm) {
		    for (j = 0; j < 4; j++) {
			bayesian2[j][nevents2] = product[j] / qnorm;
		    }
		    nevents2++;
		}
	    }
	} else {
	    double add[5];
	    double tmp;
	    
	    nbase_types = 5;
	    for (k = 0; k < 4; k++) { /* strand+chem */
		for (j = 0; j < 5; j++)
		    add[j] = 0.0;

		/* bayesian1 calc on this strand+chem */
		nevents1 = 0;
		for (j = 0; j < 5; j++) { /* type */
		    if (qhighest[j][k]) {
			/* tmp = (double)phred_table[qhighest[j][k]]; */
			tmp = (double)phred_table[qhighest[j][k]] +
			    (qcount[j][k] - 1) * DEPEND_ADD;
			prob = 1 - pow(10.0, -tmp / 10.0);
			bayesian1[0][nevents1] = (1 - prob) / 4;
			bayesian1[1][nevents1] = (1 - prob) / 4;
			bayesian1[2][nevents1] = (1 - prob) / 4;
			bayesian1[3][nevents1] = (1 - prob) / 4;
			bayesian1[4][nevents1] = (1 - prob) / 4;
			bayesian1[j][nevents1] = prob;
			nevents1++;
		    }
		}

		/* Compute bayesian1 denominator */
		qnorm = 0;
		for (j = 0; j < 5; j++) {
		    product[j] = 1;
		    for (l = 0; l < nevents1; l++) {
			product[j] *= bayesian1[j][l];
		    }
		    qnorm += product[j];
		}

		/* Fill out bayesian2 arrays */
		if (qnorm) {
		    for (j = 0; j < 5; j++) {
			bayesian2[j][nevents2] = product[j] / qnorm;
		    }
		    nevents2++;
		}
	    }
	}

	/* Compute bayesian1 denominator */
	qnorm = 0;
	highest_product = 0;
	highest_type = 5;
	for (j = 0; j < nbase_types; j++) {
	    sproduct = 1;
	    for (k = 0; k < nevents2; k++) {
		sproduct *= bayesian2[j][k];
	    }
	    qnorm += sproduct;

	    if (sproduct > highest_product) {
		highest_product = sproduct;
		highest_type = j;
	    }
	}
	    
	/* Called base has highest prob. - normalise it to get new prob. */
	prob = qnorm ? highest_product / qnorm : 0;
	if (prob < 1.0)
	    err = -10.0 * log10(1-prob);
	else
	    err = 999;

	/* Now store the results in con and qual */
	if (err - cons_cutoff100 >= -FLT_EPSILON)
	    *con++ = "ACGT*-"[highest_type];
	else
	    *con++ = '-';
	if (qualp)
	    *qualp++ = err > 99 ? 99 : err;
    }

    fflush(stdout);
    *num_frags = nf;
}
#endif

/*
 * Processes blocks of overlapping reads to produce the consensus making
 * use of the phred quality values.
 *
 * This method partitions data into sets defined by chemistry, strand, and
 * base type.
 * Ie that A(+term) and T(+term) are in different sets, but A(+term) and
 * A(+term) are in the same sets. For data within the same set we take the
 * highest confidence value plus DEPEND_ADD multiplied by the number of other
 * reads within that set.
 * Then we use Bayesian rules to determine the final sequence.
 *
 * Pads are treated as a fifth base type, but only when a column contains a
 * pad.
 *
 * If 2 consensus sequences are specified then the consensus is computed on
 * each strand independently. The confidence values are then written to qual1
 * and (if non-null) qual2.
 *
 * If 1 consensus sequence is give, but 2 quality arrays exist then the
 * second quality array is filled out to contain consensus discrepancy
 * values.
 *
 * NB: Second highest confidence is much the same as calling process_discrep
 * and using just the one confidence. Indeed the algorithms seem to be
 * identical except when pads exist. Therefore for now the 2nd highest
 * confidence code in this function is likely to be unused.
 */
static void process_frags(seq_frag *frag, int *num_frags, int from, int to,
			  int start, char *con1, float *qual1,
			  char *con2, float *qual2,
			  float cons_cutoff, int qual_cutoff)
{
    int i, j, k, nevents;
    int nf = *num_frags;
    int qual, type;
    int highest_type;
    double highest_product;
    int perfect, iub;
    char *conp = NULL;
    float *qualp = NULL;
    int qhighest[6][4];
    int qcount[6][4];
    int chem;
    double prob, err, qnorm, product;
    double bayesian[5][6*5];
    int pad_present, nbase_types;
    int cons_cutoff100 = (int)(cons_cutoff * 100 + 0.001);
    int depth;
    int num_strands = con2 ? 2 : 1, strand_num;

    for (i = from; i < to; i++) {
	for (strand_num = 0; strand_num < num_strands; strand_num++) {
	    if (strand_num == 0) {
		conp = &con1[i - start];
		if (qual1)
		    qualp = &qual1[i-start];
		else
		    qualp = NULL;
	    } else {
		conp = &con2[i - start];
		if (qual2)
		    qualp = &qual2[i-start];
		else
		    qualp = NULL;
	    }

	    /* Initialise total/count arrays */
	    memset(qhighest, 0, 6*4*sizeof(int));
	    memset(qcount, 0, 6*4*sizeof(int));
	    depth = 0;
	    
	    /* Fill qhighest and qcount arrays. Also remove frags as needed */
	    perfect = 0;
	    iub = 0;
	    pad_present = 0;
	    for (j = 0; j < nf; j++) {
		/* Type: A=0, C=1, G=2, T=3, *=4, -=5 */
		type = frag[j].qual[frag[j].start][0];
		if (type == -1)
		    type = 5;
		chem = frag[j].dir + 2 * frag[j].chem;

		/*
		 * qual == 1 implies less likely to be the called base!
		 * We map this to 2, as this is undesirable (but was the
		 * default input value for a while)
		 *
		 * qual == 0 implies "ignore this base". We do this by changing
		 * the type to dash, to force even spread of probability.
		 * Otherwise we'd actually be negatively weighting this base
		 * type.
		 */
		qual = frag[j].qual[frag[j].start][1];
		if (qual == 1)
		    qual = 2;
		else if (qual == 0)
		    type = 5;

		if (num_strands == 2) {
		    if ((frag[j].dir == 0 && strand_num == 1) ||
			(frag[j].dir == 1 && strand_num == 0))
			continue;
		}
		if (qual >= 100) {
		    perfect |= 1<<type;
		    qual = 100;
		}
		if (consensus_iub && qual >= qual_cutoff && type < 4) {
		    iub |= 1<<type;
		}
		if (qual > qhighest[type][chem])
		    qhighest[type][chem] = qual;
		qcount[type][chem]++;
		depth++;
		if (type == 4)
		    pad_present = 1;
		
		if (++frag[j].start >= frag[j].end) {
		    xfree(frag[j].qual);
		    memmove(&frag[j], &frag[j+1], (nf-j-1) * sizeof(*frag));
		    nf--; j--;
		}
	    }
	    
	    /* Handle perfect bases (qual == 100) */
	    if (perfect) {
		char base;
		
		switch(perfect) {
		case 1<<0: base = 'A'; break;
		case 1<<1: base = 'C'; break;
		case 1<<2: base = 'G'; break;
		case 1<<3: base = 'T'; break;
		case 1<<4: base = '*'; break;
		default:   base = '-'; break;
		}
		
		*conp = base;
		if (qualp)
		    *qualp = (base == '-') ? 0 : 100;
		
		continue; /* to "for (i = from; i < to; i++)" loop */
	    }
	    
	    /* Fill our bayesian matrices */
	    if (pad_present) {
		nevents = 0;
		nbase_types = 5;
		/* Dash - add to all */
		for (k = 0; k < 4; k++) {
		    if (qhighest[5][k]) {
			double tmp;
			
			/* tmp = (double)qtotal[5][k] / qcount[5][k] *
			   dependent_table[MIN(qcount[5][k],10)]; */
			
			/* tmp = (double)phred_table[qhighest[5][k]] +
			   (qcount[5][k] - 1) * DEPEND_ADD; */
			
			tmp = (double)phred_table[qhighest[5][k]] +
			    dependent_table[MIN(qcount[5][k],10)];
			
			prob = 1 - pow(10.0, -tmp / 10.0);
			
			bayesian[0][nevents] = prob / 5;
			bayesian[1][nevents] = prob / 5;
			bayesian[2][nevents] = prob / 5;
			bayesian[3][nevents] = prob / 5;
			bayesian[4][nevents] = prob / 5;
			nevents++;
		    }
		}
		/* ACGT* */
		for (j = 0; j < 5; j++) {
		    for (k = 0; k < 4; k++) {
			if (qhighest[j][k]) {
			    double tmp;
			    
			    /* tmp = (double)qtotal[j][k] / qcount[j][k] *
			       dependent_table[MIN(qcount[j][k],10)]; */
			    
			    /* tmp = (double)phred_table[qhighest[j][k]] +
			       (qcount[j][k] - 1) * DEPEND_ADD; */
			    
			    tmp = (double)phred_table[qhighest[j][k]] +
				dependent_table[MIN(qcount[j][k],10)];
			    
			    prob = 1 - pow(10.0, -tmp / 10.0);
			    
			    bayesian[0][nevents] = (1 - prob) / 4;
			    bayesian[1][nevents] = (1 - prob) / 4;
			    bayesian[2][nevents] = (1 - prob) / 4;
			    bayesian[3][nevents] = (1 - prob) / 4;
			    bayesian[4][nevents] = (1 - prob) / 4;
			    bayesian[j][nevents] = prob;
			    nevents++;
			}
		    }
		}
	    } else {
		nevents = 0;
		nbase_types = 4;
		/* Dash - add to all */
		for (k = 0; k < 4; k++) {
		    if (qhighest[5][k]) {
			double tmp;
			
			/* tmp = (double)qhighest[5][k] / qcount[5][k] *
			   dependent_table[MIN(qcount[5][k],10)]; */
			
			/* tmp = (double)phred_table[qhighest[5][k]] +
			   (qcount[5][k] - 1) * DEPEND_ADD; */

			tmp = (double)phred_table[qhighest[5][k]] +
			    dependent_table[MIN(qcount[5][k],10)];

			prob = 1 - pow(10.0, -tmp / 10.0);
			
			bayesian[0][nevents] = prob / 4;
			bayesian[1][nevents] = prob / 4;
			bayesian[2][nevents] = prob / 4;
			bayesian[3][nevents] = prob / 4;
			nevents++;
		    }
		}
		/* ACGT */
		for (j = 0; j < 4; j++) {
		    for (k = 0; k < 4; k++) {
			if (qhighest[j][k]) {
			    double tmp;
			    
			    /* tmp = (double)qhighest[j][k] / qcount[j][k] *
			       dependent_table[MIN(qcount[j][k],10)]; */
			    
			    /* tmp = (double)phred_table[qhighest[j][k]] +
			       (qcount[j][k] - 1) * DEPEND_ADD; */

			    tmp = (double)phred_table[qhighest[j][k]] +
				dependent_table[MIN(qcount[j][k],10)];

			    prob = 1 - pow(10.0, -tmp / 10.0);
			    
			    bayesian[0][nevents] = (1 - prob) / 3;
			    bayesian[1][nevents] = (1 - prob) / 3;
			    bayesian[2][nevents] = (1 - prob) / 3;
			    bayesian[3][nevents] = (1 - prob) / 3;
			    bayesian[j][nevents] = prob;
			    nevents++;
			}
		    }
		}
	    }
	    
	    qnorm = 0;
	    highest_product = 0;
	    highest_type = 5;
	    if (nevents) {    
		/* Compute denominator, used to normalise bayesian
		   probability */
		for (j = 0; j < nbase_types; j++) {
		    product = 1;
		    for (k = 0; k < nevents; k++) {
			product *= bayesian[j][k];
		    }
		    qnorm += product;
		    
		    if (product > highest_product) {
			highest_product = product;
			highest_type = j;
		    }
		}
		
		/*
		 * Called base has highest prob. - normalise it to get new
		 * prob.
		 */
		prob = qnorm ? highest_product / qnorm : 0;
		if (prob < 1.0)
		    err = MIN(-10.0 * log10(1-prob),99.4);
		else
		    err = 99.4;
		
		/* Now store the results in con and qual */
		if (err - cons_cutoff100 >= -FLT_EPSILON) {
		    *conp = "ACGT*-"[highest_type];
		} else {
		    *conp = '-';
		}

		if (consensus_iub) {
		    *conp = "NACMGRSVTWYHKDBN"[iub];
		    if (iub == 0 && highest_type == 4) {
			*conp = '*';
		    }
		}

		if (qualp)
		    *qualp = err > 99 ? 99 : err;
	    } else {
		*conp = '-';
		if (qualp)
		    *qualp = 0;
	    }

	    if (nevents && !con2 && qual2) {
		int l;
		double first=-1;
		double second=0;

		/*
		 * Loop through all base types 'l' only picking events
		 * where the called sequence matches this base type
		 * (bayesian[l][event] >= 0.25).
		 * Then compute what the consensus confidence would be
		 * from that base type alone.
		 */
		for (l = 0; l < nbase_types; l++) {
		    double prod_l = 0; 

		    qnorm = 0;
		    highest_product = 0;
		    highest_type = 5;

		    for (j = 0; j < nbase_types; j++) {
			product = 1;
			for (k = 0; k < nevents; k++) {
			    /* Only deal with bases of type 'l' */
			    if (bayesian[l][k] < 0.25)
				continue;
			    product *= bayesian[j][k];
			}
			qnorm += product;

			if (j == l)
			    prod_l = product;
		    }

		    prob = qnorm ? prod_l / qnorm : 0;
		    err = prob < 1.0 ? -10.0 * log10(1-prob) : 200;

		    if (err >= first) {
			second = first;
			first  = err;
		    } else if (err > second) {
			second = err;
		    }
		}

		/* Discrepancy = second most likely base confidence */
		qual2[i-start] = MIN(second, 199);
	    } else if (!con2 && qual2) {
		qual2[i-start] = 0;
	    }
	}
    }
    
    fflush(stdout);
    *num_frags = nf;
}

/*
 * Stirling's formula with a 1/12n correction applied to improve accuracy.
 * This seems to hold remarkably true for both low and high numbers too.
 */
double lnfact(double n) {
    /* Or Gosper's formula... 
     * return (n*ln(n) - n + ln(2*M_PI*n + M_PI/3) / 2); 
     */
    return ((n+0.5)*log(n) - n + log(2*M_PI)/2) + log(1 + 1/(12.0*n));
	/* + log(1 + 1/(288.0*n*n)); */
}

/*
 * The binomical coefficient (n,k) for n trials with k successes where
 * prob(success) = p.
 *                               k      n-k
 * P (k|n) = n! / (k! (n-k)!)   p  (1-p)
 *  p
 *
 * The coefficient we are returning here is the n! / (k! (n-k)!) bit.
 * We compute it using ln(n!) and then exp() the result back to avoid
 * excessively large numbers.
 */
double bincoef(int n, double k) {
    return exp(lnfact(n) - lnfact(k) - lnfact(n-k));
}

/*
 * Given p == 0.5 the binomial expansion simplifies a bit, so we have
 * a dedicated function for this.
 */
double binprobhalf(int n, double k) {
    return bincoef(n, k) * pow(0.5, n);
}

/*
 * This builds up discrepancy information in qual1 and qual2.
 * We use an algorithm like process_frags(), except every reading is
 * treated as independent.
 *
 * In qual1 we store the score for the second highest confidence base (with
 * all 4 base types being computed using only sequences that match that base
 * type - ie no discrepancies are involved).
 *
 * In qual2 we store a binomial coefficient computed by assuming we have a
 * population of 2 alleles with 50/50 ratio.
 */
static void process_discrep(seq_frag *frag, int *num_frags, int from, int to,
			    int start, char *con1, float *qual1,
			    char *con2, float *qual2,
			    float cons_cutoff, int qual_cutoff)
{
    int i, j, k, nevents;
    int nf = *num_frags;
    int qual, type;
    int highest_type;
    double highest_product;
    double prob, err, qnorm, product;
    double (*bayesian)[5] = NULL;
    int nbase_types;
    int count[6];

    /* Worst case is one row in bayesian[] per fragment */
    bayesian = (double (*)[5])xcalloc(*num_frags, sizeof(double)*5);

    for (i = from; i < to; i++) {
	/* Initialise total/count arrays */
	memset(count, 0, 6 * sizeof(int));
	nevents = 0;
	    
	for (j = 0; j < nf; j++) {
	    /* Type: A=0, C=1, G=2, T=3, *=4, -=5 */
	    type = frag[j].qual[frag[j].start][0];
	    if (type == -1)
		type = 5;

	    /*
	     * qual == 1 implies less likely to be the called base!
	     * We map this to 2, as this is undesirable (but was the
	     * default input value for a while)
	     *
	     * qual == 0 implies "ignore this base". We do this by changing
	     * the type to dash, to force even spread of probability.
	     * Otherwise we'd actually be negatively weighting this base
	     * type.
	     */
	    qual = frag[j].qual[frag[j].start][1];
	    if (qual == 1)
		qual = 2;
	    else if (qual == 0)
		type = 5;

	    if (qual >= qual_cutoff) {
		/* Add to bayesian array. Skip if it's "-" */
		prob = 1 - pow(10.0, -qual / 10.0);
		if (type != 5) {
		    for (k = 0; k < 5; k++) {
			bayesian[nevents][k] = (1 - prob) / 4;
		    }
		    bayesian[nevents][type] = prob;
		    count[type]++;
		    nevents++;
		}
	    }
		
	    if (++frag[j].start >= frag[j].end) {
		xfree(frag[j].qual);
		memmove(&frag[j], &frag[j+1], (nf-j-1) * sizeof(*frag));
		nf--; j--;
	    }
	}
	    
	if (nevents && qual2) {
	    int l;
	    double first=-1;
	    double second=0;
	    int first_count = 0;
	    int second_count = 0;

	    /*
	     * Loop through all base types 'l' only picking events
	     * where the called sequence matches this base type
	     * (bayesian[event][l] >= 0.25).
	     * Then compute what the consensus confidence would be
	     * from that base type alone.
	     */
	    nbase_types = 5; /* = 4;   FIXME: ignore indels for now */
	    for (l = 0; l < nbase_types; l++) {
		double prod_l = 0; 

		qnorm = 0;
		highest_product = 0;
		highest_type = 5;

		for (j = 0; j < nbase_types; j++) {
		    product = 1;
		    for (k = 0; k < nevents; k++) {
			/* Only deal with bases of type 'l' */
			if (bayesian[k][l] < 0.25)
			    continue;
			product *= bayesian[k][j];
		    }
		    if (j == l) {
			prod_l = product;
		    } else {
			qnorm += product;
		    }
		}

		prob = (qnorm + prod_l) ? qnorm / (qnorm + prod_l) : 1;
		err = prob ? -10.0 * log10(prob) : 1000;
		if (err > 1000) err = 1000;

		if (err >= first) {
		    second = first;
		    second_count = first_count;
		    first  = err;
		    first_count = count[l];
		} else if (err > second) {
		    second = err;
		    second_count = count[l];
		}
	    }

	    /* Qual1 = second highest consensus confidence */
	    if (qual1)
		qual1[i-start] = second;


	    /* Qual2 = a measure of how close we are to the 50/50 ratio */
	    if (qual2 && first_count && second_count) {
		/*
		 * Max of 10 deep to avoid oversampling and becoming
		 * too dependent on the 'optimal' 50% ratio.
		 * The reason is that we may not achieve 50% ratio even
		 * with a huge sample depth due to systematic sequencing
		 * errors.
		 */
		double cnt = first_count + second_count;
		double snd = second_count;
		if (cnt > 10) {
		    snd *= 10/cnt;
		    cnt = 10;
		}

		qual2[i-start] = binprobhalf(cnt, snd) * cnt;
	    } else if (qual2) {
		qual2[i-start] = 0;
	    }
	}
    }

    xfree(bayesian);
    
    fflush(stdout);
    *num_frags = nf;
}

/*
 * Processes blocks of overlapping reads to produce the consensus making
 * use of the phred quality values.
 *
 * This method partitions data into sets defined by chemistry, strand, and
 * base type.
 * Ie that A(+term) and T(+term) are in different sets, but A(+term) and
 * A(+term) are in the same sets. For data within the same set we take the
 * highest confidence value plus DEPEND_ADD multiplied by the number of other
 * reads within that set.
 * Then we use Bayesian rules to determine the final sequence.
 *
 * Pads are treated by working out firstly whether a base exists, and then
 * what it is (if applicable). - TODO!!!!
 */
#if 0
static void process_frags_7(seq_frag *frag, int *num_frags, int from, int to,
			    int start, char *con1, float *qual1,
			    char *con2, float *qual2,
			    float cons_cutoff, int qual_cutoff)
{
    int i, j, k, nevents;
    int nf = *num_frags;
    int qual, type;
    int highest_type;
    double highest_product;
    int perfect;
    char *conp = NULL;
    float *qualp = NULL;
    int qhighest[6][4];
    int qcount[6][4];
    int chem;
    double prob, err, qnorm, product;
    double bayesian[5][6*5];
    int pad_present;
    int cons_cutoff100 = (int)(cons_cutoff * 100 + 0.001);
    int depth;
    int num_strands = con2 ? 2 : 1, strand_num;
    double base_prob;

    for (i = from; i < to; i++) {
	for (strand_num = 0; strand_num < num_strands; strand_num++) {
	    if (strand_num == 0) {
		conp = &con1[i - start];
		if (qual1)
		    qualp = &qual1[i-start];
		else
		    qualp = NULL;
	    } else {
		conp = &con2[i - start];
		if (qual2)
		    qualp = &qual2[i-start];
		else
		    qualp = NULL;
	    }

	    /* Initialise total/count arrays */
	    memset(qhighest, 0, 6*4*sizeof(int));
	    memset(qcount, 0, 6*4*sizeof(int));
	    depth = 0;
	    
	    /* Fill qhighest and qcount arrays. Also remove frags as needed */
	    perfect = 0;
	    pad_present = 0;
	    for (j = 0; j < nf; j++) {
		type = frag[j].qual[frag[j].start][0];
		if (type == -1)
		    type = 5;
		/* Type: A=0, C=1, G=2, T=3, *=4, -=5 */
		qual = frag[j].qual[frag[j].start][1];
		chem = frag[j].dir + 2 * frag[j].chem;
		if (num_strands == 2) {
		    if ((frag[j].dir == 0 && strand_num == 1) ||
			(frag[j].dir == 1 && strand_num == 0))
			continue;
		}
		if (qual >= 100) {
		    perfect |= 1<<type;
		    qual = 100;
		}
		if (qual > qhighest[type][chem])
		    qhighest[type][chem] = qual;
		qcount[type][chem]++;
		depth++;
		if (type == 4)
		    pad_present = 1;
		
		if (++frag[j].start >= frag[j].end) {
		    xfree(frag[j].qual);
		    memmove(&frag[j], &frag[j+1], (nf-j-1) * sizeof(*frag));
		    nf--; j--;
		}
	    }
	    
	    /* Handle perfect bases (qual == 100) */
	    if (perfect) {
		char base;
		
		switch(perfect) {
		case 1<<0: base = 'A'; break;
		case 1<<1: base = 'C'; break;
		case 1<<2: base = 'G'; break;
		case 1<<3: base = 'T'; break;
		case 1<<4: base = '*'; break;
		default:   base = '-'; break;
		}
		
		*conp = base;
		if (qualp)
		    *qualp = (base == '-') ? 0 : 100;
		
		continue; /* to "for (i = from; i < to; i++)" loop */
	    }
	    
	    /* Fill our bayesian matrices */
	    if (pad_present) {
		/*
		 * Firstly work out the probability that the base exists.
		 * This is done using Bayes with 2 outcomes (gap vs base).
		 */
		nevents = 0;
		for (j = 0; j < 6; j++) { /* j==5  => pad, (ACGT*-) */
		    for (k = 0; k < 4; k++) {
			if (qhighest[j][k]) {
			    double tmp;

			    tmp = (double)phred_table[qhighest[j][k]] +
				(qcount[j][k] - 1) * DEPEND_ADD;
			    prob = 1 - pow(10.0, -tmp / 10.0);

			    if (j == 4) { /* pad */
				bayesian[0][nevents] = 1-prob;
				bayesian[1][nevents] = prob;
			    } else { /* ACGT- */
				bayesian[0][nevents] = prob;
				bayesian[1][nevents] = 1-prob;
			    }
			    nevents++;
			}
		    }
		}

		/* Compute probability of a base existing */
		qnorm = 0;
		product = 1;
		for (k = 0; k < nevents; k++) {
		    product *= bayesian[1][k]; /* pad */
		}
		qnorm += product;
		product = 1;
		for (k = 0; k < nevents; k++) {
		    product *= bayesian[0][k]; /* base */
		}
		qnorm += product;
		base_prob = qnorm ? product / qnorm : 0;
	    } else {
		base_prob = 1.0;
	    }
		
	    /* Most likely to be a base? If so, work out which */
	    if (base_prob >= 0.5) {
		/* Work out the base probabilities */
		nevents = 0;
		/* Dash - add to all */
		for (k = 0; k < 4; k++) {
		    if (qhighest[5][k]) {
			double tmp;
			
			/* tmp = (double)qhighest[5][k] / qcount[5][k] *
			   dependent_table[MIN(qcount[5][k],10)]; */
			
			tmp = (double)phred_table[qhighest[5][k]] +
			    (qcount[5][k] - 1) * DEPEND_ADD;
			prob = 1 - pow(10.0, -tmp / 10.0);
			
			bayesian[0][nevents] = prob / 4;
			bayesian[1][nevents] = prob / 4;
			bayesian[2][nevents] = prob / 4;
			bayesian[3][nevents] = prob / 4;
			nevents++;
		    }
		}
		/* ACGT */
		for (j = 0; j < 4; j++) {
		    for (k = 0; k < 4; k++) {
			if (qhighest[j][k]) {
			    double tmp;
			    
			    /* tmp = (double)qhighest[j][k] / qcount[j][k] *
			       dependent_table[MIN(qcount[j][k],10)]; */
			    
			    tmp = (double)phred_table[qhighest[j][k]] +
				(qcount[j][k] - 1) * DEPEND_ADD;
			    prob = 1 - pow(10.0, -tmp / 10.0);
			    
			    bayesian[0][nevents] = (1 - prob) / 3;
			    bayesian[1][nevents] = (1 - prob) / 3;
			    bayesian[2][nevents] = (1 - prob) / 3;
			    bayesian[3][nevents] = (1 - prob) / 3;
			    bayesian[j][nevents] = prob;
			    nevents++;
			}
		    }
		}
		
		/* Compute denominator, used to normalise bayesian probability */
		qnorm = 0;
		highest_product = 0;
		highest_type = 5;
		if (nevents) {    
		    for (j = 0; j < 4; j++) {
			product = 1;
			for (k = 0; k < nevents; k++) {
			    product *= bayesian[j][k];
			}
			qnorm += product;
			
			if (product > highest_product) {
			    highest_product = product;
			    highest_type = j;
			}
		    }
		    
		    /*
		     * Called base has highest prob. - normalise it to get new
		     * prob.
		     */
		    prob = qnorm ? highest_product / qnorm : 0;
		    if (prob < 1.0)
			err = -10.0 * log10(1-base_prob * prob);
		    else
			err = 999;
		    
		    /* Now store the results in con and qual */
		    if (err - cons_cutoff100 >= -FLT_EPSILON)
			*conp = "ACGT*-"[highest_type];
		    else
			*conp = '-';
		    if (qualp)
			*qualp = err > 99 ? 99 : err;
		} else {
		    *conp = '-';
		    if (qualp)
			*qualp = 0;
		}

	    } else { /* base_prob < 0.5 */
		err = -10.0 * log10(base_prob);

		*conp = (err - cons_cutoff100 >= -FLT_EPSILON) ? '*' : '-';
		if (qualp)
		    *qualp = err > 99 ? 99 : err;
	    }
	}
    }
    
    fflush(stdout);
    *num_frags = nf;
}
#endif

/*
 * Processes blocks of overlapping reads to produce the consensus making
 * use of the phred quality values.
 *
 * This uses Graham's notes to model valid and invalid alignments. This is
 * used to spot when we have errors introduced due to sequencing artifacts.
 */
#if 0
static void process_frags_8(seq_frag *frag, int *num_frags, int from, int to,
			    int start, char *con1, float *qual1,
			    char *con2, float *qual2,
			    float cons_cutoff, int qual_cutoff)
{
    int i, j, k, strand_num;
    int nf = *num_frags;
    char *conp = NULL;
    float *qualp = NULL;
    int num_strands = con2 ? 2 : 1;
    /* [Base_type][2*(model for strand 0)+(model for strand 1)] */
    /* 0/0=0, 0/1=1, 1/0=2, 1/1=3 */
    double product[5][4];
    int type;
    double qual, qual3;
    double model_prob[2];
    double highest_prod;
    int highest_type;
    double denominator;
    int cons_cutoff100 = (int)(cons_cutoff * 100 + 0.001);
    int perfect;

    model_prob[0] = 0.9999;
    model_prob[1] = 0.0001;

    for (i = from; i < to; i++) {
	for (strand_num = 0; strand_num < num_strands; strand_num++) {
	    if (strand_num == 0) {
		conp = &con1[i - start];
		if (qual1)
		    qualp = &qual1[i-start];
		else
		    qualp = NULL;
	    } else {
		conp = &con2[i - start];
		if (qual2)
		    qualp = &qual2[i-start];
		else
		    qualp = NULL;
	    }

	    /* Initialise products to 1.0 */
	    for (j = 0; j < 5; j++) {
		product[j][0] = model_prob[0] * model_prob[0];
		product[j][1] = model_prob[0] * model_prob[1];
		product[j][2] = model_prob[1] * model_prob[0];
		product[j][3] = model_prob[1] * model_prob[1];
	    }

	    perfect = 0;

	    /*
	     * Step through fragments multiplying the corresponding product[]
	     * element.
	     */
	    for (j = 0; j < nf; j++) {
		int i1, i2;
		
		type = frag[j].qual[frag[j].start][0];
		if (type == -1)
		    type = 5;
		qual = frag[j].qual[frag[j].start][1];
		
		if (qual == 100) {
		    /* Perfect base */
		    perfect |= 1<<type;
		}
		
		/* Skip unknown bases */
		if (type == 5)
		    goto next;
		
		qual = 1 - pow(10.0, -qual / 10.0);
		qual3 = (1-qual) / 3;
		
		if (frag[j].dir == 0) {
		    i1 = 1; i2 = 2;
		} else {
		    i1 = 2; i2 = 1;
		}

		for (k = 0; k < 5; k++) { /* k == base_type */
		    double q;
		    
		    q = (k == type) ? qual : qual3;
		    product[k][0] *= q;
		    product[k][i1] *= q;
		    product[k][i2] *= 0.25;
		    product[k][3] *= 0.25;
		}

	    next:
		if (++frag[j].start >= frag[j].end) {
		    xfree(frag[j].qual);
		    memmove(&frag[j], &frag[j+1], (nf-j-1) * sizeof(*frag));
		    nf--; j--;
		}
	    }
	    
	    /*
	     * Check for perfect bases, including disagreeing perfect bases.
	     */
	    if (perfect) {
		char base;
		
		switch(perfect) {
		case 1<<0: base = 'A'; break;
		case 1<<1: base = 'C'; break;
		case 1<<2: base = 'G'; break;
		case 1<<3: base = 'T'; break;
		case 1<<4: base = '*'; break;
		default:   base = '-'; break;
		}
		
		*conp = base;
		if (qualp)
		    *qualp = (base == '-') ? 0 : 100;
		
		continue; /* to "for (strand_num = 0"... loop */
	    }

	    /*
	     * Sum elements from each product[base_type][0..3] and store back
	     * in product[base_type][0].
	     */
	    for (k = 0; k < 5; k++) {
		product[k][0] += product[k][1] + product[k][2] + product[k][3];
	    }
	    
	    /*
	     * Normalise and find highest sum.
	     */
	    highest_type = 5;
	    highest_prod = 0;
	    denominator = 0;
	    for (k = 0; k < 5; k++) {
		if (product[k][0] > highest_prod) {
		    highest_prod = product[k][0];
		    highest_type = k;
		}
		denominator += product[k][0];
	    }
	    
	    /* Now write the consensus + confidence */
	    qual = denominator ? highest_prod / denominator : 0;
	    if (qual < 1.0)
		qual = -10 * log10(1-qual);
	    else
		qual = 99;
	    
	    if (qual - cons_cutoff100 >= -FLT_EPSILON)
		*conp = "ACGT*-"[highest_type];
	    else
		*conp = '-';
	    if (qualp)
		*qualp = qual > 99 ? 99 : qual;
	}
    }

    *num_frags = nf;
}
#endif

/*
 * Processes blocks of overlapping reads to produce the consensus making
 * use of the phred quality values.
 *
 * This uses Graham's notes to model valid and invalid alignments. This is
 * used to spot when we have errors introduced due to sequencing artifacts.
 *
 * As process_frags_8, except breaks into four groups instead of two.
 */
#if 0
static void process_frags_9(seq_frag *frag, int *num_frags, int from, int to,
			    int start, char *con1, float *qual1,
			    char *con2, float *qual2,
			    float cons_cutoff, int qual_cutoff)
{
    int i, j, k, l, strand_num;
    int nf = *num_frags;
    char *conp = NULL;
    float *qualp = NULL;
    int num_strands = con2 ? 2 : 1;
    /* [Base_type][2*(model for strand 0)+(model for strand 1)...] */
    /* 0/0=0, 0/1=1, 1/0=2, 1/1=3 */
    double product[5][16];
    int type;
    double qual, qual3;
    double mp[2];
    double highest_prod;
    int highest_type;
    double denominator;
    int cons_cutoff100 = (int)(cons_cutoff * 100 + 0.001);
    int perfect;
    int chem;

    mp[0] = 0.999;
    mp[1] = 1-mp[0];

    for (i = from; i < to; i++) {
	for (strand_num = 0; strand_num < num_strands; strand_num++) {
	    if (strand_num == 0) {
		conp = &con1[i - start];
		if (qual1)
		    qualp = &qual1[i-start];
		else
		    qualp = NULL;
	    } else {
		conp = &con2[i - start];
		if (qual2)
		    qualp = &qual2[i-start];
		else
		    qualp = NULL;
	    }

	    /*
	     * Initialise products to appropriate combinations of model prob
	     * (mp) elements 0 (conf is OK) and 1 (conf is random).
	     * This does bit twiddling such that element product[x][0] is
	     * for all chem/strand combinations being correct and product[x][1]
	     * is for all chem/strand combinations being incorrect.
	     */
	    for (j = 0; j < 5; j++) {
		for (k = 0; k < 16; k++) {
		    product[j][k] = mp[(k&1)!=0] * mp[(k&2)!=0] *
			mp[(k&4)!=0] * mp[(k&8)!=0];
		}
	    }

	    perfect = 0;

	    /*
	     * Step through fragments multiplying the corresponding product[]
	     * element.
	     */
	    for (j = 0; j < nf; j++) {
		type = frag[j].qual[frag[j].start][0];
		if (type == -1)
		    type = 5;
		qual = frag[j].qual[frag[j].start][1];
		chem = 1 << (2 * frag[j].dir + frag[j].chem);
		
		if (qual == 100) /* Perfect base */
		    perfect |= 1<<type;
		
		/* Skip unknown bases */
		if (type != 5) {
		    qual = 1 - pow(10.0, -qual / 10.0);
		    qual3 = (1-qual) / 3;
		
		    for (k = 0; k < 5; k++) { /* k == base_type */
			double q;
		    
			q = (k == type) ? qual : qual3;
			for (l = 0; l < 16; l++)
			    product[k][l] *= (l & chem) ? .25 : q;
		    }
		}

		if (++frag[j].start >= frag[j].end) {
		    xfree(frag[j].qual);
		    memmove(&frag[j], &frag[j+1], (nf-j-1) * sizeof(*frag));
		    nf--; j--;
		}
	    }
	    
	    /*
	     * Check for perfect bases, including disagreeing perfect bases.
	     */
	    if (perfect) {
		char base;
		
		switch(perfect) {
		case 1<<0: base = 'A'; break;
		case 1<<1: base = 'C'; break;
		case 1<<2: base = 'G'; break;
		case 1<<3: base = 'T'; break;
		case 1<<4: base = '*'; break;
		default:   base = '-'; break;
		}
		
		*conp = base;
		if (qualp)
		    *qualp = (base == '-') ? 0 : 100;
		
		continue; /* to "for (strand_num = 0"... loop */
	    }

	    /*
	     * Sum elements from each product[base_type][0..3] and store back
	     * in product[base_type][0].
	     */
	    for (k = 0; k < 5; k++) {
		product[k][0] += product[k][1] + product[k][2] + product[k][3];
	    }
	    
	    /*
	     * Normalise and find highest sum.
	     */
	    highest_type = 5;
	    highest_prod = 0;
	    denominator = 0;
	    for (k = 0; k < 5; k++) {
		if (product[k][0] > highest_prod) {
		    highest_prod = product[k][0];
		    highest_type = k;
		}
		denominator += product[k][0];
	    }
	    
	    /* Now write the consensus + confidence */
	    qual = denominator ? highest_prod / denominator : 0;
	    if (qual < 1.0)
		qual = -10 * log10(1-qual);
	    else
		qual = 99;
	    
	    if (qual - cons_cutoff100 >= -FLT_EPSILON)
		*conp = "ACGT*-"[highest_type];
	    else
		*conp = '-';
	    if (qualp)
		*qualp = qual > 99 ? 99 : qual;
	}
    }

    *num_frags = nf;
}
#endif

static int calc_contig_info_phred(int contig, int start, int end,
				  char *con, float *qual,
				  char *con2, float *qual2,
				  float cons_cutoff, int qual_cutoff,
				  void (*process_frags)(seq_frag *frag,
							int *num_frags,
							int from,
							int to,
							int start,
							char *con1,
							float *qual1,
							char *con2,
							float *qual2,
							float cons_cutoff,
							int qual_cutoff),
				  int (*info_func)(int         job,
						   void       *mydata,
						   info_arg_t *theirdata),
				  void *info_data)
{
    register int (*gel_qual)[2];	/* gel quality */
    int i_start;
    register int i_end;
    info_arg_t info;
    seq_frag *frag;
    int max_frags, num_frags;
    int last_pos = start;
    int frag_start, frag_end;

    /*
     * Quick sanity checks
     */
    if (start > end)
	return -1;

    info.contig_info.contig = contig;
    info_func(GET_CONTIG_INFO, info_data, &info);

    if (con)   memset(con,   '-', end-start+1);
    if (qual)  memset(qual,  0,   (end-start+1)*sizeof(float));
    if (con2)  memset(con2,  '-', end-start+1);
    if (qual2) memset(qual2, 0,   (end-start+1)*sizeof(float));

    /*
     * find left most gel in our region
     */
    info.gel_info.gel = info.contig_info.leftgel;
    do {
	info_func(GET_GEL_INFO, info_data, &info);
    } while (info.gel_info.position + info.gel_info.length < start &&
	     (info.gel_info.gel = info.gel_info.next_right));

    max_frags = 10;
    num_frags = 0;
    if (NULL == (frag = (seq_frag *)xmalloc(max_frags * sizeof(*frag))))
	return -1;

    do {
	if (info.gel_info.gel == 0)
	    break;

	i_start = info.gel_info.position < start ?
	    start - info.gel_info.position : 0;
	i_end = info.gel_info.position + info.gel_info.length > end ?
	    end - info.gel_info.position + 1: info.gel_info.length;

	if (i_end <= i_start)
	    goto fetch_next; /* ok - so it's a cheat! */

	/*
	 * Get quality info for gel
	 */
	if ((int (*)[2])-1 ==
	    (gel_qual = get_gel_qual(info.gel_info.gel,
				     i_start, i_end,
				     info_func, info_data))) {
	    xfree(frag);
	    return -1;
	}

#if PRINT_CONS
	printf("Gel %d, pos %d, length %d, right %d, istart %d iend %d\n",
	       info.gel_info.gel,
	       info.gel_info.position,
	       info.gel_info.length,
	       info.gel_info.next_right,
	       i_start,
	       i_end);
#endif


	/* Add gel to fragment list */
	if (num_frags >= max_frags) {
	    max_frags *= 2;
	    frag = (seq_frag *)xrealloc(frag, max_frags * sizeof(*frag));
	    if (frag == NULL)
		return -1;
	}

	frag[num_frags].qual = gel_qual;
	frag[num_frags].start = 0;
	frag[num_frags].end = i_end - i_start;
	frag[num_frags].gel = info.gel_info.gel;
	frag[num_frags].dir = info.gel_info.complemented;
	frag[num_frags].chem = info.gel_info.as_double ? 1 : 0;
	frag[num_frags].template = info.gel_info.template;
	num_frags++;

	fetch_next:

	last_pos = info.gel_info.position;

	/*
	 * fetch next gel
	 */
	info.gel_info.gel = info.gel_info.next_right;
	if (info.gel_info.gel) {
	    info_func(GET_GEL_INFO, info_data, &info);
	    frag_end = info.gel_info.position > end+1
		? end+1 : info.gel_info.position;
	} else {
	    frag_end = end+1;
	}
	frag_start = last_pos < start ? start : last_pos;
	if (frag_end >= frag_start)
	    process_frags(frag, &num_frags, frag_start, frag_end,
			  start, con, qual, con2, qual2,
			  cons_cutoff, qual_cutoff);

	/* loop until past end */

    } while (info.gel_info.position <= end);

    xfree(frag);
    return 0;
}


/*
 * Processes the fragments using the older consensus algorithm.
 */
#if 0
static void process_frags_old(seq_frag *frag, int *num_frags, int from, int to,
			      int start, char *con1, char *con2,
			      float *qual1, float *qual2,
			      float cons_cutoff, int qual_cutoff)
{
    int i, j;
    int nf = *num_frags;
    int conq[7];
    int qual, type;
    unsigned char base;
    int t;
    char *con = &con1[from-start];

#if PRINT_CONS
    printf("For region %d to %d. Num_frags=%d\n", from, to, nf);
    for (i = 0; i < nf; i++) {
	printf("    Fragment %d: Gel %d from %d to %d\n",
	       i, frag[i].gel, frag[i].start, frag[i].end);
    }
#endif

    for (i = from; i < to; i++) {
#if PRINT_CONS
	printf("%5d: ", i);
#endif
	conq[0] = 0;
	conq[1] = 0;
	conq[2] = 0;
	conq[3] = 0;
	conq[4] = 0;
	conq[5] = 0;
	conq[6] = 0;
	for (j = 0; j < nf; j++) {
	    type = frag[j].qual[frag[j].start][0];
	    qual = frag[j].qual[frag[j].start][1];
	    if (type != -1) {
		if (qual >= qual_cutoff) {
		    if (qual != 100)
			conq[type] += qual;
		    else
			conq[6] |= 1<<type;
		} else {
		    qual = 0;
		}
	    } else {
		conq[0] += qual/4;
		conq[1] += qual/4;
		conq[2] += qual/4;
		conq[3] += qual/4;
	    }
	    conq[5] += qual;
		
#if PRINT_CONS
	    printf("%c", "-ACGT*"[type+1]);
#endif
	    if (++frag[j].start >= frag[j].end) {
		xfree(frag[j].qual);
		memmove(&frag[j], &frag[j+1], (nf-j-1) * sizeof(*frag));
		nf--; j--;
	    }
	}

	if (conq[6]) {
	    switch(conq[6]) {
	    case 1<<0: /* A */
		base = 0;
		break;
	    case 1<<1: /* C */
		base = 1;
		break;
	    case 1<<2: /* G */
		base = 2;
		break;
	    case 1<<3: /* T */
		base = 3;
		break;
	    case 1<<4: /* * */
		base = 4;
		break;
	    default:
		base = 5;
		break;
	    }

	} else {
	    /* Check for best base - defaults to '-' */
	    t = 0; base = 5; /* default to '-' */
	    if (t < conq[0]) t = conq[0], base = 0;
	    if (t < conq[1]) t = conq[1], base = 1;
	    if (t < conq[2]) t = conq[2], base = 2;
	    if (t < conq[3]) t = conq[3], base = 3;
	    if (t < conq[4]) t = conq[4], base = 4;
	    
	    if (conq[5] > 0) {
		if ((float)conq[base]/conq[5] - cons_cutoff < -FLT_EPSILON)
		    base = 5;
	    } else {
		base = 5;
	    }
	}

	*con++ = "ACGT*-"[base];
#if PRINT_CONS
	printf(" => %c\n", "ACGT*-"[base]);
#endif
    }

    fflush(stdout);
    *num_frags = nf;
}
#endif



int calc_discrepancies(int   contig,
		       int   start,
		       int   end,
		       float *qual1,
		       float *qual2,
		       float cons_cutoff,
		       int   qual_cutoff,
		       int (*info_func)(int        job,
					void       *mydata,
					info_arg_t *theirdata),
		       void *info_data)
{
    init_clookup();

    qual_cutoff_tmp = (qual_cutoff == QUAL_DEFAULT)
	? qual_cutoff_def : qual_cutoff;

    if (-1 ==calc_contig_info_phred(contig, start, end,
				    NULL, qual1, NULL, qual2,
				    cons_cutoff, qual_cutoff_tmp,
				    process_discrep,
				    info_func, info_data)) {
	return -1;
    }

    return 0;
}
