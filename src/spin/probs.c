#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#include "sip_hash.h"
#include "readpam.h"
#include "sequence_formats.h"
#include "misc.h"
#include "dna_utils.h"
#include "probs.h"

/********************************************************************/
/*	POLYNOMIALS */

#define MAX_POLY 40000
/* #define MAX_POLY_MATRIX 2400 */
#define MAX_POLY_MATRIX 25
#define SMALL_POLY 1.0e-30
#define ZERO_SMALL(a) ( ((a) < (SMALL_POLY)) ?  (0.0) : (a) )

static double polya[MAX_POLY], polyb[MAX_POLY], polyc[MAX_POLY];
static int num_poly_terms,size_poly_step;
static double poly_matrix_scale,poly_matrix_f1;
static double poly_matrix_min;

static int poly_matrix[MAX_POLY_MATRIX][MAX_POLY_MATRIX];

/*
   num_poly_terms
   size_poly_step
       poly_matrix_min
       poly_matrix_f1
       poly_matrix_scale;
*/
/*
	MAX_POLY  max polynomial terms for polya,polyb,polyc which are
	used for multiplying out the probabilities
	MAX_POLY_MATRIX  max size of the score matrix/weight matrix that
	is scaled for use by the probability calculations. The original
	score matrices or weight matrices are supplied by the application
	programs, but the probability calculations rescale them. The
	rescaled version is stored in poly_matrix which is treated as
	an [n][m] array. Its dimensions are poly_matrix_cols, and 
	poly_matrix_rows which have to be set. We need a routine to set up:

	1. copy matrix to poly_matrix (plus other values like span length)
	2. poly_matrix_min, poly_matrix_max, poly_matrix_scale are min, max 
	and scale factors used to convert poly_matrix values to original scores
	poly_matrix_f1 is a scale factor, poly_matrix_add is used to
	rescale the matrix values. These last two need to be set up by
	the initialisation and scaling routines
	seq1_freq, seq2_freq are the base/amino acid compositions of the
	sequences being analysed. Sum (seq _freq) = 1.0

	Each matrix type - say dna score matrix, protein score matrix,
	dna weight matrix, protein weight matrix, etc needs to have its
	own initialisation routine. We also need more than one routine
	to combine the polynomials in the correct way - only score matrix
	done so far.

	The calculation puts the probability of getting score k in the 
	score matrix into polya[k]. This probability depends on the 
	chance of finding in the sequence the characters that give score k,
	which is a function of their frequencies (seq_freq). In the case
	of a score matrix the probabilities of getting score k will be of
	the form seq1_freq[i]*seq2_freq[j] where character types i and j
	give score k. There may be several sets of i and j that give score k
	so all those terms are added into polya[k]. For a weight matrix the
	probabilities of getting score k will be of the form seq1_freq[i]
	but again several values of i could give score k. 
	Once the probabilities for the first set of possible scores are in
	polya, the probabilities for the next set - say the next column in
	a weight matrix - are put in polyb. Then polya is multiplied by
	polyb and the result stored in polya (via polyc). Then the last 
	two steps ie probabilities in polyb then multiply by polya; are 
	repeated until the end of the motif is reached. At this point the
	values in polya are changed to be the cumulative probabilities.
	ie the probability in polya[k] is the probability of scoring <=k.

	Each round of polynomial multiplication increases the number of
	used elements in polya by the value of the highest value in the
	scaled score/weight matrix. The array size is handled by 
	num_poly_terms and size_poly_step, which must be initialised to
	equal the highest value in the scaled matrix.
*/


int poly_mult (void) {

    int i,j,max_terms;

    max_terms = num_poly_terms + size_poly_step;

    if ( max_terms > MAX_POLY ) return -1;

    for(i=0;i<=max_terms;i++) polyc[i] = 0.0;

    for (i=0;i<=num_poly_terms;i++) {
	for (j=0;j<=size_poly_step;j++) {
	    polyc[i+j] += polya[i] * polyb[j];
	}
    }

    num_poly_terms = max_terms;

    for(i = 0; i <= num_poly_terms; i++) 
	polya[i] = ZERO_SMALL ( polyc[i] );
    return 0;
}

int prob1 ( int job, 
	   int matrix[], 
	   int dimension1, 
	   int dimension2,
	   int span_length, 
	   double *seq1_freq, 
	   double *seq2_freq) 
{

    /* calculate probabilities for 
       span length span_length and 
       matrix matrix. matrix is sent
       as a 1 d array but we also receive
       the dimensions dimension1 and dimension2
       for use here.
    */

    int i,j,k;
    int matrix_rows, matrix_cols;
    double matrix_max,matrix_add;
    /* match_prob function needs:
       poly_matrix_min
       poly_matrix_f1
       poly_matrix_scale;
    */

    matrix_rows = dimension1;
    matrix_cols = dimension2;
    poly_matrix_f1 = span_length;
    poly_matrix_min = 99999.;
    matrix_max = -99999.;

    /* copy to 2D and get range */
    for(i = 0, k = 0; i < matrix_rows; i++) {
	for(j = 0; j < matrix_cols; j++) {
	    /*
	       printf("%d %d %d %f %f\n",i, j, poly_matrix[i][j], 
	       poly_matrix_min, matrix_max);
	       */
	    poly_matrix[i][j] = matrix[k++];
	    poly_matrix_min = MIN(poly_matrix[i][j],poly_matrix_min);
	    matrix_max = MAX(poly_matrix[i][j], matrix_max);
	}
    }

/*
    for(i = 0, k = 0; i < matrix_rows; i++) {
	for(j = 0; j < matrix_cols; j++) {
	    printf("%d ",poly_matrix[i][j]);
	}
	printf("\n");
    }
*/
    /* scale to min = 0 */

    for(i = 0; i < matrix_rows; i++) {
	for(j = 0; j < matrix_cols; j++) {
	    poly_matrix[i][j] -= poly_matrix_min;
	}
    }
    
    poly_matrix_scale = 1.0;
    matrix_add = poly_matrix_min * poly_matrix_f1;
    size_poly_step = num_poly_terms = matrix_max - poly_matrix_min; 
    /* previous line different from fortran */
    /* 
       printf(" min %f max %f scale %f terms %d\n", poly_matrix_min,
       matrix_max,
       poly_matrix_scale,num_poly_terms);
       */
    /* generate the probabilities */

    for (i=0;i<MAX_POLY;i++) polya[i] = polyb[i] = 0.0;
    for (i=0;i<matrix_rows;i++) {
	for ( j=0;j<matrix_cols;j++) {
	    k = poly_matrix[i][j];
	    polyb[k] = polya[k] += seq1_freq[i] * seq2_freq[j];
	}
    }

    for (i = 1;i < span_length;i++) {
	j = poly_mult();
	if ( j ) return j;
    }
    if ( (2 == job) || (4 == job) ) {
	for(i=num_poly_terms;i>-1;i--) polya[i] += polya[i+1];
    }
    if ( (3 == job) || (4 == job) ) {
	for(i = 0; i <= num_poly_terms; i++) {
	    polyb[i] = (i/poly_matrix_scale) + matrix_add;
	}
    }
    return 0;
}

/*	return probability for score score */

double 
match_prob (double score,
	    double min_prob) {

    int i;

    i = (score - poly_matrix_min * poly_matrix_f1) * poly_matrix_scale;
    if ( (i >= 0) && (i < MAX_POLY) && (polya[i] > min_prob) ) {
	return polya[i];
    }
    return -1.0;
}

/* return probability for score score */
double 
match_prob2 ( double score) {

    int i;

    i = (score - poly_matrix_min * poly_matrix_f1) * poly_matrix_scale;
    if ( (i >= 0) && (i < MAX_POLY) ) 
	return polya[i];
    return -1.0;
}


/*
 * FIXME - find out proper values for B and Z
 * taken from genetic_code.c but added B and Z
 */
double sip_av_protein_comp[] = { 
    8.3,/* A */
    0.0,/* B */
    1.7,/* C */
    5.3,/* D */
    6.2,/* E */
    3.9,/* F */
    7.2,/* G */
    2.2,/* H */
    5.2,/* I */
    5.7,/* K */
    9.0,/* L */
    2.4,/* M */
    4.4,/* N */
    5.1,/* P */
    4.0,/* Q */
    5.7,/* R */
    6.9,/* S */
    5.8,/* T */
    6.6,/* V */
    1.3,/* W */
    3.2,/* Y */
    0.0,/* Z */
    0.0,/* X */
    0.0,/* * */
    0.0 /* - */
};


int FindProbs(char *seq1,
	      char *seq2,
	      int seq1_lreg,      /* left cutoff of seq1 */
	      int seq1_rreg,      /* right cutoff of seq1 */
	      int seq2_lreg,      /* left cutoff of seq2 */
	      int seq2_rreg,      /* right cutoff of seq2 */
	      int span_length,
	      int seqtype,
	      int use_av_comp)
{
#define MAX_MATCHES 100000
    int i,j,k;
    int *probs_mat;
    double *seq1_freq, *seq2_freq;
    int seq1_len, seq2_len;

    set_char_set(seqtype);

    if (use_av_comp) {
	set_char_set(2); /* PROTEIN */
	if (NULL == (seq1_freq = (double *)xmalloc(char_set_size * 
						   sizeof(double))))
	    return -1;
	if (NULL == (seq2_freq = (double *)xmalloc(char_set_size * 
						   sizeof(double))))
	    return -1;

	for (i = 0; i < char_set_size; i++) {
	    seq1_freq[i] = sip_av_protein_comp[i] / 100;
	    seq2_freq[i] = sip_av_protein_comp[i] / 100;
	}
    } else {
	if (NULL == (seq1_freq = (double *)xmalloc(char_set_size * 
						   sizeof(double))))
	    return -1;
	if (NULL == (seq2_freq = (double *)xmalloc(char_set_size * 
						   sizeof(double))))
	    return -1;

	/* initialise frequency arrays */
	for (i = 0; i < char_set_size; i++) 
	    seq1_freq[i] = seq2_freq[i] = 0.0;


	seq1_len = seq1_rreg - seq1_lreg + 1;
	seq2_len = seq2_rreg - seq2_lreg + 1;
	for (i = seq1_lreg - 1; i < seq1_rreg; i++) {
	    seq1_freq[char_lookup[seq1[i]]] += 1.0;
	}
	for (i = 0; i < char_set_size; i++) {
	    seq1_freq[i] = seq1_freq[i] / seq1_len;
	}
	for (i = seq2_lreg - 1; i < seq2_rreg; i++) {
	    seq2_freq[char_lookup[seq2[i]]] += 1.0;
	}
	for (i = 0; i < char_set_size; i++) {
	    seq2_freq[i] = seq2_freq[i] / seq2_len;
	}
    }
#ifdef DEBUG
    printf("start %d %d end %d %d\n", seq1_lreg, seq2_lreg, seq1_rreg, seq2_rreg);
    for (i = 0; i < char_set_size; i++) {
	printf("seq2_freq[%d]=%f\n", i, seq2_freq[i]);
    }
#endif
    /* get matrix ready to send to probs */
    if ( ! ( probs_mat = (int *) xmalloc(sizeof(int)* 
					char_set_size * char_set_size))) { 
	return -1;
    }

    for(i = 0, k = 0; i < char_set_size; i++) {
	for(j = 0; j < char_set_size; j++) {
	    probs_mat[k++] = score_matrix[i][j];
	}
    }

    /* Ok send it */

    i = prob1 ( 4, probs_mat, char_set_size, char_set_size, span_length, 
	       seq1_freq, seq2_freq);

    set_char_set(seqtype);

    /* free up the temporary memory */
    free (probs_mat);
    xfree(seq1_freq);
    xfree(seq2_freq);
    return 0;
}

/*
 * writes a list of scores and associated probabilities over a range between
 * min_score and max_score to the output window
 */
void ListProbs(char *seq1,
	       char *seq2,
	       int seq1_lreg,      /* left cutoff of seq1 */
	       int seq1_rreg,      /* right cutoff of seq1 */
	       int seq2_lreg,      /* left cutoff of seq2 */
	       int seq2_rreg,      /* right cutoff of seq2 */
	       int span_length,
	       int seqtype,
	       int min_score,
	       int max_score,
	       int *score_hist)
{
    double prob;
    int i = min_score;
    int seq1_len, seq2_len;

    seq1_len = seq1_rreg - seq1_lreg + 1;
    seq2_len = seq2_rreg - seq2_lreg + 1;

    FindProbs(seq1, seq2, seq1_lreg, seq1_rreg, seq2_lreg, seq2_rreg, 
	      span_length, seqtype, 0);

     /* now the probs are ready to be interrogated */
    for (i = min_score; i <= max_score; i++) {
	prob =  match_prob2((double)i);
	vmessage("score %4d probability %.2e expected %12.0f observed %d\n", 
		 i, prob, (double)seq1_len * seq2_len * prob, 
		 score_hist[i-min_score]);
    }
}

/*
 * for the find matching words algorithm:
 * writes a list of scores and associated probabilities over a range between
 * min_score and max_score to the output window
 */
void ListIdentityProbs(char *seq1,
		       char *seq2,
		       int seq1_lreg,      /* left cutoff of seq1 */
		       int seq1_rreg,      /* right cutoff of seq1 */
		       int seq2_lreg,      /* left cutoff of seq2 */
		       int seq2_rreg,      /* right cutoff of seq2 */
		       int seqtype,
		       int min_score,
		       int max_score,
		       int *score_hist)
{
    double prob;
    int i = min_score;
    int seq1_len,seq2_len;

    seq1_len = seq1_rreg - seq1_lreg + 1;
    seq2_len = seq2_rreg - seq2_lreg + 1;

     /* now the probs are ready to be interrogated */
    for (i = min_score; i <= max_score; i++) {
	FindProbs(seq1, seq2, seq1_lreg, seq1_rreg, seq2_lreg, seq2_rreg, i, 
		  seqtype, 0);
	prob =  match_prob2((double)i);
	vmessage("score %4d probability %.2e expected %12.0f observed %d\n", 
		 i, prob, (double)seq1_len * seq2_len * prob, 
		 score_hist[i-min_score]);
    }
}

double FindExpectedProb(char *seq1,
			char *seq2,
			int seq1_lreg,      /* left cutoff of seq1 */
			int seq1_rreg,      /* right cutoff of seq1 */
			int seq2_lreg,      /* left cutoff of seq2 */
			int seq2_rreg,      /* right cutoff of seq2 */
			int span_length,
			int seqtype)
{
    double min_prob = 0.000000001;
    double prob;
    int seq1_len, seq2_len;

    seq1_len = seq1_rreg - seq1_lreg + 1;
    seq2_len = seq2_rreg - seq2_lreg + 1;

    FindProbs(seq1, seq2, seq1_lreg, seq1_rreg, seq2_lreg, seq2_rreg,
	      span_length, seqtype, 0);

     /* now the probs are ready to be interrogated */
    prob = match_prob((double)span_length, min_prob);

    /* prob can end up being -1, in which case, return min_prob */
    if (prob == -1) {
	return ((double)seq1_len * seq2_len * min_prob);
    } else {
	return ((double)seq1_len * seq2_len * prob);
    }

#ifdef REMOVE
    while ((prob = match_prob((double)i, min_prob)) > 0.0) {
	printf("prob %f %f\n", prob, (double)seq1_len * seq2_len * prob);
	if (i == span_length) {
	    return ((double)seq1_len * seq2_len * prob);
	}
	i++;
    }
    return -1.0;
#endif
}

/*
 * returns the score for the probability at min_prob for span_length
 */
int
FindScore(int seq1_len,
	  int seq2_len,
	  int span_length,
	  int num_matches)
{
    double min_prob = 0.000000001;
    double prob;
    int i = poly_matrix_min * span_length;
    double expected_matches;

    /* execute FindProbs from tcl */
    /* FindProbs(seq1, seq2, seq1_len, seq2_len, span_length, seqtype); */

     /* now the probs are ready to be interrogated */
    while ((prob = match_prob((double)i, min_prob)) > 0.0) {
	expected_matches = (double)seq1_len * seq2_len * prob;
#ifdef DEBUG
	printf("i %d prob %f expected %f num_matches %d\n", 
	       i, prob, expected_matches, num_matches);
#endif
	if (expected_matches < (double)num_matches){
	    return i; 
	}
	i++;
    }
    return i-1;
}
