#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <time.h>

#include "misc.h"
#include "dna_utils.h"
#include "hash_lib.h"
#include "consen.h"

#define MINMAT 12

/* Size to divide the h->diag[] into for purposes of delayed initialisation */
#define DSZ 0x800

static int dna_hash8_lookup[256];

void set_hash8_lookupn(void) {

    int i;

    for (i=0;i<256;i++) dna_hash8_lookup[i] = 4;

    dna_hash8_lookup['a'] = 0;
    dna_hash8_lookup['c'] = 1;
    dna_hash8_lookup['g'] = 2;
    dna_hash8_lookup['t'] = 3;
    dna_hash8_lookup['A'] = 0;
    dna_hash8_lookup['C'] = 1;
    dna_hash8_lookup['G'] = 2;
    dna_hash8_lookup['T'] = 3;
    dna_hash8_lookup['*'] = 0;
}

int init_hash8n (
		int max_seq,
		int max_diagonal,
		int word_length,
		int max_matches,
		int min_match,
		int job,
		Hash **h) {
    
  int size_hash;
    /* max_seq is the longest sequence, but may be either seq1 or seq2
     * max_diagonal is the length of the second longest sequence
     * if a sequence is going to be compared against itself then
     * max_seq = max_diagonal
     */

    set_hash8_lookupn ();
    
    if ( ! (*h = (Hash *) xmalloc ( sizeof(Hash) ))) return -2;

    if (word_length < 8)
	word_length = 4;
    else if (word_length < 12)
	word_length = 8;
    else if (word_length < 14)
	word_length = 12;
    else
	word_length = 14;

    size_hash = 1 << (2 * word_length);
    
    if ( HASH_JOB_BLKS & job ) {
	min_match = MAX(min_match,word_length);
    }

    (*h)->values1 = NULL;
    (*h)->values2 = NULL;
    (*h)->counts = NULL;
    (*h)->last_word = NULL;
    (*h)->diag = NULL;
    (*h)->hist = NULL;
    (*h)->expected_scores = NULL;
    (*h)->diag_match = NULL;
    (*h)->block_match = NULL;
    (*h)->max_matches = max_matches;
    (*h)->min_match = min_match;
    (*h)->matches = 0;
    (*h)->word_length = word_length;
    (*h)->size_hash = size_hash;
    (*h)->fast_mode = 0;
    (*h)->filter_words = 0;
    
    if ( ! ((*h)->values1 = (int *) xmalloc ( sizeof(int)*(max_seq) ))) 
	return -2;

    if ( ! ((*h)->values2 = (int *) xmalloc ( sizeof(int)*(max_diagonal) ))) 
	return -2;

    /* at present 3 modes are used:
     * 1) either we do a quick comparison and quick alignment using blocks
     * 2) or we do a quick comparison and slow alignment
     * 3) or simply search for repeats or fortran sequence assembly
     * for 1 need
     * HASH_JOB_DIAG   1
     * HASH_JOB_BLKS  16
     * ie job = 17
     * 
     * for 2 need
     * HASH_JOB_DIAG   1
     * HASH_JOB_HIST   2
     * HASH_JOB_EXPD   4
     * HASH_JOB_DMTCH  8
     * HASH_JOB_BLKS  16
     * ie job = 31
     *
     * for 3 need
     * HASH_JOB_DIAG   1
     * ie job = 1
     */

    if( (job != 1) && (job != 17) && (job != 31) && (job != 33)) return -2;

    if (!(job & HASH_JOB_COUNTLESS)) {
	if ( ! ((*h)->counts = (int *) xcalloc ((*h)->size_hash,sizeof(int) )))
	    return -2;
    }
  
    if ( ! ((*h)->last_word = (int *) xcalloc ( (*h)->size_hash,sizeof(int) )))
	return -2;    

    if ( HASH_JOB_DIAG & job ) {

	if ( ! ((*h)->diag = (int *) xmalloc ( sizeof(int)*(max_seq+max_diagonal+DSZ) ))) 
	    return -2;
    }

    if ( HASH_JOB_HIST & job ) {

	if(!((*h)->hist = (int *) xmalloc(sizeof(int) * (max_seq+max_diagonal+DSZ)))) 
	    return -2;
    }

    if ( HASH_JOB_EXPD & job ) {

	if(!((*h)->expected_scores = (int *) xmalloc(sizeof(int) * 
							(max_diagonal))))
	    return -2;
    }

    if ( HASH_JOB_DMTCH & job ) {

	if(!((*h)->diag_match = (Diag_Match *) xmalloc(sizeof(Diag_Match) * 
							max_matches)))
	    return -2;
	(*h)->max_matches = max_matches;
    }

    if ( HASH_JOB_BLKS & job ) {

	if(!((*h)->block_match = (Block_Match *) xmalloc(sizeof(Block_Match) * 
							max_matches)))
	    return -2;
	(*h)->max_matches = max_matches;
    }
    return 0;
}

void free_hash8n ( Hash *h ) {
    
  if ( h->values1 ) xfree ( h->values1 );
  if ( h->values2 ) xfree ( h->values2 );
  if ( h->counts ) xfree ( h->counts );
  if ( h->last_word ) xfree ( h->last_word );
  if ( h->diag )         xfree ( h->diag );
  if ( h->hist )         xfree ( h->hist );
  if ( h->expected_scores )    xfree ( h->expected_scores );
  if ( h->diag_match )    xfree ( h->diag_match );
  if ( h->block_match )    xfree ( h->block_match );
  xfree ( h );
}

static int hash_word_n(char *seq, int *start_base, int seq_len,
		       int word_length, unsigned int *uword) {
    /* 	given a sequence seq, return the hash value for the first word 
     *  after start_base that does not contain an unknown char. Tell 
     *  the caller where this is. If we reach the end of the seq set
     *  start_base and return -1.
     */
    int lstart_base = *start_base;
    int end_base    = lstart_base + word_length;
    int i;
    int base_index;
    unsigned int luword = 0;
    unsigned int mask = (1 << (word_length * 2)) - 1;

    if (seq_len < end_base) return -1;

    for (i = lstart_base; i < end_base; i++) {
	base_index = dna_hash8_lookup[(unsigned char)seq[i]];
	if (4 == base_index) {
	    /*	We've hit an unknown char, so let's start again */
	    lstart_base = i + 1;
	    end_base = lstart_base + word_length;
	    if (seq_len < end_base) {
		*start_base = lstart_base;
		return -1;
	    }
	    luword = 0;
	    i = lstart_base - 1;
	} else {
	    luword = (luword << 2) | base_index;
	}
    }
    *start_base = lstart_base;
    *uword = luword & mask;
    return 0;
}

static int hash_seq_n(char *seq, int *hash_values,
		      int seq_len, int word_length) {
    /* Given a sequence seq, return an array of hash values.
     * If we cannot find at least one word to hash on we return -1
     * otherwise we return 0.
     */
    int i, j, k, ret;
    int start_base = 0, prev_start_base, base_index;
    unsigned int uword = 0;
    unsigned int mask = (1 << (2 * word_length)) - 1;

    if ( seq_len < word_length ) return -1;

    /*	Get the hash value for the first word that contains no unknowns */
    if (hash_word_n(seq, &start_base, seq_len, word_length, &uword)) return -1;
    for (i = 0; i < start_base; i++) hash_values[i] = -1;

    /*	Now do the rest of the sequence */
    hash_values[start_base] = uword & mask;
    k = seq_len - word_length + 1;

    for (i = start_base + 1, j = start_base + word_length; i < k; i++, j++) {
	base_index = dna_hash8_lookup[(unsigned char)seq[j]];
	if (4 == base_index) {
	    /*	We've hit an unknown char, so let's start again */
	    prev_start_base = i;
	    start_base = j + 1;
	    ret = hash_word_n(seq, &start_base, seq_len, word_length, &uword);
	    for (i = prev_start_base; i < start_base; i++) hash_values[i] = -1;
	    if (ret) return 0; /* Couldn't restart */
	    hash_values[start_base] = uword & mask;
	    i = start_base;
	    j = i + word_length - 1;
	} else {
	    uword = (uword << 2) | base_index;
	    hash_values[i] = uword & mask;
	}
    }

    return 0;
}

/* ---------------------------------------------------------------------------
 * WORD SIZE = 14
 */

int hash_word14n ( char *seq, int *start_base, int seq_len, int word_length,
		unsigned int *uword) {
    return hash_word_n(seq, start_base, seq_len, word_length, uword);
}


int hash_seq14n ( char *seq, int *hash_values, int seq_len, int word_length) {
    return hash_seq_n(seq, hash_values, seq_len, word_length);
}

/* ---------------------------------------------------------------------------
 * WORD SIZE = 12
 */

static int hash_word12n ( char *seq, int *start_base, int seq_len,
			  int word_length, unsigned int *uword) {
    return hash_word_n(seq, start_base, seq_len, word_length, uword);
}


static int hash_seq12n (char *seq, int *hash_values,
			int seq_len, int word_length) {
    return hash_seq_n(seq, hash_values, seq_len, word_length);
}


/* ---------------------------------------------------------------------------
 * WORD SIZE = 8
 */

int hash_word8n ( char *seq, int *start_base, int seq_len, int word_length,
		unsigned short *uword) {
    unsigned int uiword = *uword;
    int ret = hash_word_n(seq, start_base, seq_len, word_length, &uiword);
    *uword = uiword & 0xffff;
    return ret;
}


int hash_seq8n ( char *seq, int *hash_values, int seq_len, int word_length) {
    return hash_seq_n(seq, hash_values, seq_len, word_length);
}

/* ---------------------------------------------------------------------------
 * WORD SIZE = 4
 */

int hash_word4n ( char *seq, int *start_base, int seq_len, int word_length,
		unsigned char *uword) {
    unsigned int uiword = *uword;
    int ret = hash_word_n(seq, start_base, seq_len, word_length, &uiword);
    *uword = uiword & 0xff;
    return ret;
}


int hash_seq4n ( char *seq, int *hash_values, int seq_len, int word_length) {
    return hash_seq_n(seq, hash_values, seq_len, word_length);
}

/* ---------------------------------------------------------------------------
 */

void store_hashn ( Hash *h ) {

    /*store the hash values in values: put number of occurrences of
     *each hash value in counts; put the array position of the last 
     *occurrence of each hash value in last_word, and previous
     *occurrences in values[last_word]. i.e. values is used for TWO
     *purposes: the initial indexes and then the positions.
     *Note that words containing unknown characters (like '-') are given
     *hash value -1. So we skip them here, and they are ignored.
     */


    int nw;
    register int i,j,n;

    for (i=0; i<h->size_hash; i++) {
	h->counts[i] = 0;
	h->last_word[i] = 0;
    }
    j = h->seq1_len - h->word_length + 1;
    for (i = 0; i < j; i++) {
	n = h->values1[i];
	if (-1 != n) {
	    nw = h->counts[n];
	    if (nw != 0) {
		h->values1[i] = h->last_word[n];
	    }
	    h->last_word[n] = i;
	    h->counts[n]++;
	}
    }
}

/* As per store_hashn() but no h->counts[] array. Terminate list on -1 */
void store_hashn_nocount ( Hash *h ) {
    register int i,j,n;

    for (i=0; i<h->size_hash; i++) {
	h->last_word[i] = -1;
    }
    j = h->seq1_len - h->word_length + 1;
    for (i = 0; i < j; i++) {
	n = h->values1[i];
	if (n == -1)
	    continue;
	
	h->values1[i] = h->last_word[n];
	h->last_word[n] = i;
    }
}

int diagonal_length(int seq1_len, int seq2_len, int diagonal_number) {

    /* return the length of a diagonal given the diagonal number */

    /*
       0123456789 seq1->
       1
       2
       3        diagonal numbers
       4       ^
       5     /
       6    2
       7  1
       0
       s
       e
       q
       2
       |
       v
       */

    if(diagonal_number < seq1_len)
	return MIN((diagonal_number + 1), (
	       MIN(seq1_len, seq2_len)));
    else
	return MIN((seq1_len + seq2_len - 1 - diagonal_number), (
	       MIN(seq1_len, seq2_len)));
}

void diagonal_intercepts (int diagonal_number, int seq1_len, int seq2_len,
			  int *seq1_intercept, int *seq2_intercept ) {
    /* given a diagonal number, return the intercepts on the seq1 and seq2 axes */
    
    
    if ( diagonal_number < seq1_len ) {
	*seq1_intercept = seq1_len - diagonal_number - 1;
	*seq2_intercept = 0;
    }
    else {
	*seq2_intercept = diagonal_number + 1 - seq1_len;
	*seq1_intercept = 0;
    }
    
}

/*	POLYNOMIALS */


#define MAX_POLY 20
#define SMALL_POLY 1.0e-30
#define ZERO_SMALL(a) ( ((a) < (SMALL_POLY)) ?  (0.0) : (a) )

typedef struct poly_ {
    double a[MAX_POLY];
    double b[MAX_POLY];
    double c[MAX_POLY];
    int num_terms;
    int size_step;
    int rows;
    int cols;
} Poly;

int poly_mult (Poly *poly) {
    int i,j,max_terms;

    max_terms = poly->num_terms + poly->size_step;
    if ( max_terms > MAX_POLY ) return -1;
    for(i=0;i<=max_terms;i++) poly->c[i] = 0.0;
    for (i=0;i<=poly->num_terms;i++) {
	for (j=0;j<=poly->size_step;j++) {
	    poly->c[i+j] += poly->a[i] * poly->b[j];
	}
    }
    poly->num_terms = max_terms;
    for(i = 0; i <= poly->num_terms; i++) 
	poly->a[i] = ZERO_SMALL ( poly->c[i] );
    return 0;
}

double prob_word ( int word_length, double comp[]) {
    Poly p, *poly;
    int i,j,k;

    poly = &p;
    poly->rows = poly->cols = 4;
    poly->size_step = poly->num_terms = 1;
    for (i=0;i<MAX_POLY;i++) poly->a[i] = poly->b[i] = 0.0;
    for (i=0;i<poly->rows;i++) {
	for ( j=0;j<poly->cols;j++) {
	    k = (i==j) ? 1:0;
	    poly->b[k] = poly->a[k] += comp[i] * comp[j];
	}
    }
    for (i = 1;i < word_length;i++) {
	j = poly_mult(poly);
	if ( j ) return -1.0;
    }
    /*for(i=poly->num_terms;i>-1;i--) poly->a[i] += poly->a[i+1];*/
    return poly->a[word_length];
}

void p_comp(double comp[], char *seq, int seq_len) {
    int i;
    double t;
    for(i=0;i<5;i++)comp[i]=0.0;
    if(seq_len<1)return;
    for(i=0;i<seq_len;i++) comp[dna_hash8_lookup[(unsigned)seq[i]]]++;
    for(i=0,t=0.0;i<4;i++)t+=comp[i];
    if(t>0.0) {
	for(i=0;i<4;i++)comp[i]/=t;
    }
}

int poisson_diagonals(int min_diag, int max_diag, int word_length,
		      double max_prob_in, int *expected_scores,
		      double comp[]) {
    int		diagonal_length, hits;
    double	expected_hits, sum_probs;
    double	prob_remaining;
    double	p_w;
    double limit;
    int not_found;
    double max_prob, frac, big;

    double z, emz, x;

    /* Assume the diagonal scores obey a poisson distribution
     * with an expected number of hits dependent on the diagonal
     * length and the probability of a single matching word p_w
     *
     * if Z is expected or average number of occurrences
     * then e**(-Z) {(1), (Z), (Z**2)/2!), (Z**3)/3!, ...}
     * are the probabilities of observing 0, 1, 2, 3, ... occurrences.
     *
     * We are interested in the highest scores on a diagonal and
     * this calculation seems to produce curves that fit the
     * results from simulated data reasonably well.
     * The simulation performed large numbers of comparisons 
     * of random 25% acgt composition sequences and recorded the
     * highest score obtained on each diagonal.
     * 
     * Note that for long words (eg 8) the values change very slowly
     * ie 100's of consecutive diagonals will have the same score.
     * 
     */

    for(hits = 0; hits < max_diag; hits++) {
	expected_scores[hits] = max_diag;
    }

    /* on digital alpha tcl can only manage values >1.0e-38 ! */
    limit = 1.0e-37;
    if(max_prob_in<limit)max_prob_in = limit;
    limit = 1.0e-14;
    max_prob = max_prob_in;
    if(max_prob<limit)max_prob = limit;

    /*comp[0]=comp[1]=comp[2]=comp[3]=0.25;*/
    p_w = prob_word(word_length,comp);

    if (p_w < 0.0) return -1;

    /*p_w = 1.0/pow(4.0, (float) word_length);*/
    /*printf("min_diag %d max_diag %d max_prob %e max_prob_in %e p_w %e\n",
    min_diag,max_diag,max_prob, max_prob_in, p_w);*/

    big = DBL_MAX/1000000000000.0;
    for(diagonal_length = min_diag; diagonal_length < max_diag; diagonal_length++) {

	expected_hits = (double) diagonal_length * p_w;
	/*limit = (DBL_MAX/1000000000000.0)/expected_hits;*/
	limit = big/expected_hits;

	/* sum the probabilities until a higher number of hits is
	 * sufficiently improbable
	 */

	x = 1.0;
	z = expected_hits;
	emz = exp(-1 * z);
	for(hits = 1, sum_probs = emz, not_found = 1; hits < diagonal_length; hits++) {
/*	printf("limit %e x %e sum %e rem %e hits %d\n",limit,x,sum_probs,prob_remaining,hits);*/

	    if ( x > limit ) break;
	    x *= (z/hits);
	    sum_probs += x * emz;
	    prob_remaining = 1.0 - sum_probs;
	    if(prob_remaining < max_prob) {
		expected_scores[diagonal_length] = hits;
		not_found = 0;
		break;
	    }
	}
	if ( not_found ) {
	    /* printf("not found %d %d\n",diagonal_length,hits); */
	    expected_scores[diagonal_length] = hits;
	}
    }


    if(max_prob_in<max_prob) {
	frac = 1.0+0.033*log10(max_prob/max_prob_in);
	for(hits = 0; hits < max_diag; hits++) {
	    expected_scores[hits] *= frac;
	}
    }
    /*
    for(hits = 0; hits < max_diag; hits++) {
	printf("hits %d exp %e\n",hits,expected_scores[hits]);
    }
    */
    return 0;
}

int best_intercept ( Hash *h, int *seq1_i, int *seq2_i ) {

    /* routine to find the best intercept given a set of matches
     * between a pair of sequences. This is done iteratively:
     * we find the centre of gravity, remove the outlier, and
     * repeat until the best match only remains.
     * NOTE: THIS ROUTINE ZEROES THE MATCH ARRAY
     */

    double t, sum_scores, sum_moment, c_o_g, furthest;
    int match_no, matches_left, outlier = 0;

    for(matches_left = h->matches; matches_left > 1; matches_left--) {
	sum_moment = sum_scores = 0.0;
	for(match_no = 0; match_no < h->matches; match_no++) {
	    if ( h->diag_match[match_no].prob > 0.0 ) {
		sum_moment += h->diag_match[match_no].pos *
		    h->diag_match[match_no].prob;
		sum_scores += h->diag_match[match_no].prob;
	    }
	}
	if ( sum_scores == 0 ) {
	    fprintf(stderr, "FATAL: best_intecept has sum_scores of 0\n");
	    return 0;
	}
	c_o_g = sum_moment/sum_scores;
	for(match_no = 0, furthest = 0.0; match_no < h->matches; match_no++) {
	    if ( h->diag_match[match_no].prob > 0.0 ) {
		if ((t=fabs(c_o_g - h->diag_match[match_no].pos)) > furthest) {
		    outlier = match_no;
		    furthest = t;
		}
	    }
	}
	h->diag_match[outlier].prob = 0.0;
    }
    for(match_no = 0; match_no < h->matches; match_no++) {
	if ( h->diag_match[match_no].prob > 0.0 ) {
	    diagonal_intercepts (h->diag_match[match_no].pos, h->seq1_len, 
				 h->seq2_len, seq1_i, seq2_i);
	    break;
	}
    }
    return 1;
}

int hash_seqn (Hash *h, int job) {
    assert(job == 1 || job == 2);
    assert(h->word_length >= 4 && h->word_length < 15);
    if ( job == 1 ) {
	return hash_seq_n(h->seq1, h->values1, h->seq1_len, h->word_length);
    } else {
	return hash_seq_n(h->seq2, h->values2, h->seq2_len, h->word_length);
    }
}


void remdup ( int **seq1_match, int **seq2_match, int **len_match, 
	      int offset, int *n_matches ) {
    
    /*	routine to remove duplicates from a list of repeats. It
	also removes the self match.
	Input: a list of *n_match match positions and match
	lengths.  Output: a list in which all duplicates and the selfmatch
	are removed.  *n_match is set to the new number of matches, or -1
	for error.  */
    
    register int i;
    int *index_ptr,k,keep;
    
    if ( *n_matches < 1 ) return;

    if ( ! ( index_ptr = (int *) xmalloc ( sizeof(int)*(*n_matches) ))) {
 	*n_matches = -1;
	return;
    }
    for ( i=0,k=0;i<*n_matches; i++) {
	if ( (*seq1_match)[i+offset] > (*seq2_match)[i+offset] ) {
	    index_ptr[k] = i+offset; 
	    k += 1;
	}
    }
    for ( i=0; i<k; i++) {
 	keep = index_ptr[i];

	(*seq1_match)[i+offset] = (*seq1_match)[keep];
	(*seq2_match)[i+offset] = (*seq2_match)[keep];
	(*len_match)[i+offset] = (*len_match)[keep];
    } 
    *n_matches = k;
    if ( index_ptr ) free ( index_ptr );
}

void make_reverse ( int **seq2_match, int **len_match,
		   int n_matches, int seq2_len, int offset) {
    
    int i;
    
    for (i = 0; i< n_matches; i++) {
	(*seq2_match)[i+offset] = seq2_len - (*seq2_match)[i+offset] - (*len_match)[i+offset]+2;
	
    }
}

typedef struct Edit_pair {
    int *S1;
    int *S2;
    int size;
    int next1;
    int next2;
} EDIT_PAIR;


int update_edit_pair ( EDIT_PAIR *edit_pair, OVERLAP *overlap ) {
    int i,j;

    /*
    printf(">update_edit_pair  s1 %d s2 %d\n",overlap->s1_len,overlap->s2_len);
    printf("  S1: "); for(i=0; i<overlap->s1_len; i++)
	printf(" %d", overlap->S1[i]);
    printf("\n");
    printf("  S2: "); for(i=0; i<overlap->s2_len; i++)
	printf(" %d", overlap->S2[i]);
    printf("\n");
    */

    if ( overlap->s1_len ) {
	if ( (edit_pair->size - edit_pair->next1) < overlap->s1_len ) return -1;
	for (i=edit_pair->next1,j=0;j<overlap->s1_len;i++,j++) {
	    edit_pair->S1[i] = overlap->S1[j];
	}
	edit_pair->next1 += overlap->s1_len;
	xfree ( overlap->S1 );
	overlap->S1 = NULL;
	overlap->s1_len = 0;
    }
    if ( overlap->s2_len ) {
	if ( (edit_pair->size - edit_pair->next2) < overlap->s2_len ) return -1;
	for (i=edit_pair->next2,j=0;j<overlap->s2_len;i++,j++) {
	    edit_pair->S2[i] = overlap->S2[j];
	}
	edit_pair->next2 += overlap->s2_len;
	xfree ( overlap->S2 );
	overlap->S2 = NULL;
	overlap->s2_len = 0;
    }
    return 0;
}

int block_to_edit_pair ( EDIT_PAIR *edit_pair, int length ) {
    if ( (edit_pair->size - edit_pair->next1) < 1 ) return -1;
    edit_pair->S1[edit_pair->next1++] = length;
    if ( (edit_pair->size - edit_pair->next2) < 1 ) return -1;
    edit_pair->S2[edit_pair->next2++] = length;
    return 0;
}

int align_bit ( ALIGN_PARAMS *params, OVERLAP *overlap, EDIT_PAIR *edit_pair) {

    int l1, l2;

    l1 = overlap->seq1_len;
    l2 = overlap->seq2_len;

    /*
    printf("seq1 len %d '%.*s'\n",
	   overlap->seq1_len, overlap->seq1_len, overlap->seq1);
    printf("seq2 len %d '%.*s'\n",
	   overlap->seq2_len, overlap->seq2_len, overlap->seq2);
    printf("l1=%d l2=%d band=%d\t",
	   l1, l2, params->band);
    fflush(stdout);
    */

    if (l1 == 1 && l2 == 1) {
    	edit_pair->S1[edit_pair->next1++] = 1;
    	edit_pair->S2[edit_pair->next2++] = 1;
    } else if ((l1 > 0) && (l2 > 0 )) {
	if (affine_align(overlap,params)) return -1;
	if ( update_edit_pair ( edit_pair, overlap)) return -1;
	/* printf("%f\n", overlap->score); */
    } else if (l1 > 0 ) {
	if ( edit_pair->next2 == edit_pair->size ) return -1;
	edit_pair->S2[edit_pair->next2++] = -l1;

	if ( edit_pair->next1 == edit_pair->size ) return -1;
	edit_pair->S1[edit_pair->next1++] = l1;
	/* printf("no-score\n"); */
    } else if (l2 > 0 ) {
	if ( edit_pair->next1 == edit_pair->size ) return -1;
	edit_pair->S1[edit_pair->next1++] = -l2;
	
	if ( edit_pair->next2 == edit_pair->size ) return -1;
	edit_pair->S2[edit_pair->next2++] = l2;
	/* printf("no-score\n"); */
    } else if (l1 == 0 && l2 == 0) {
	return 0;
    } else {
	printf("impossible alignment?\n");
    }

    return 0;
}

void destroy_edit_pair (EDIT_PAIR *edit_pair) {
    if (edit_pair) {
	if (edit_pair->S1) xfree (edit_pair->S1);
	if (edit_pair->S2) xfree (edit_pair->S2);
	xfree (edit_pair);
    }
}

EDIT_PAIR *create_edit_pair(int size) {
    EDIT_PAIR *edit_pair;

    if(NULL == (edit_pair = (EDIT_PAIR *) xmalloc(sizeof(EDIT_PAIR)))) {
	verror(ERR_WARN, "create_edit_pair", "xmalloc failed");
	return NULL;
    }

    if ( ! ((edit_pair->S1 = (int *) xmalloc ( sizeof(int)*(size) )))) {
	destroy_edit_pair (edit_pair);
	verror(ERR_WARN, "create_edit_pair", "xmalloc failed");
	return NULL;
    }
    if ( ! ((edit_pair->S2 = (int *) xmalloc ( sizeof(int)*(size) )))) {
	destroy_edit_pair (edit_pair);
	verror(ERR_WARN, "create_edit_pair", "xmalloc failed");
	return NULL;
    }
    edit_pair->next1 = 0;
    edit_pair->next2 = 0;
    edit_pair->size = size;
    return edit_pair;
}

/*
 * Designed with affine_align() in mind to ensure we don't go over
 * the MAX_MEMORY2 limit.
 */
int set_band_blocks(int seq1_len, int seq2_len) {
    return MIN(9990000.0/MIN(seq1_len,seq2_len),
    	       MAX(10,(MIN(seq1_len,seq2_len)*0.1)));

    return MIN(9990000.0/MIN(seq1_len,seq2_len),
    	       MAX(30,(MIN(seq1_len,seq2_len)*0.35)));
}

/* FAST MODE */
int set_band_blocks_fast(int seq1_len, int seq2_len) {
    return MIN(9990000.0/MIN(seq1_len,seq2_len),
	       MAX(10,(MIN(seq1_len,seq2_len)*0.05)));

}

/*
 * Given a bunch of hash hits we can work out what the best case scenario is
 * for the alignment:
 * +1 for every base in a matching block.
 * +1 for every base in MIN(len1,len2) in the matrix between blocks.
 * -1 for every base in ABS(len1-len2) in the matrix between blocks.
 * -1 for every base in MIN(len1,len2)/(word_len-1) between blocks.
 *
 * Ie if we assume the sections that didn't match our hash of word_len 20
 * were all {19 match, 1 mismatch}* segments then we can estimate the total
 * alignment mismatch percentage.
 * When we have very different sized axis between hashed blocks we know
 * at best we'll have at least len1-len2 pads to insert, plus whatever
 * differences the sequences have too.
 *
 * The intention here is to avoid performing alignments that we know we'll
 * reject later on anyway.
 */
int min_mismatch(Hash *h, int *mis_p, int *mat_p) {
    int i;
    int p1 = 0, p2 = 0;
    int match = 0, mismatch = 0, worst_mat=0, worst_mis=0;
    int l1, l2;

    if (h->matches == 0)
	return 100;

    /*
    printf("=== STEP 0 === %d\n",h->matches);
    for (i=0;i<h->matches;i++) {
	printf("i %d %d %d %d %d %d %d\n",i,
	       h->block_match[i].pos_seq1,
	       h->block_match[i].pos_seq2,
	       h->block_match[i].length,
	       h->block_match[i].diag,
	       h->block_match[i].best_score,
	       h->block_match[i].prev_block);
    }
    puts("");
    */

    /* First block start non-anchored while end is */
    l1 = h->block_match[0].pos_seq1;
    l2 = h->block_match[0].pos_seq2;
    mismatch = MIN(l1, l2)/h->min_match + 1;
    match =  MIN(l1, l2) - mismatch;
    match += h->block_match[0].length;

    worst_mat = h->block_match[0].length;
    worst_mis = MIN(l1, l2);

    p1 = h->block_match[0].pos_seq1 + h->block_match[0].length;
    p2 = h->block_match[0].pos_seq2 + h->block_match[0].length;

    /* Inbetween blocks => both ends anchored */
    for (i = 1; i < h->matches; i++) {
	/* Gap from p1 to [].pos_seq1 & p2 to [].pos_seq2 */
	l1 = h->block_match[i].pos_seq1 - p1;
	l2 = h->block_match[i].pos_seq2 - p2;
	match += MIN(l1, l2) - MIN(l1, l2)/h->min_match;
	mismatch += MAX(ABS(l1 - l2), MIN(l1, l2)/h->min_match+1);
	match += h->block_match[i].length;

	worst_mat += h->block_match[i].length;
	worst_mis += MAX(l1, l2);

	p1 = h->block_match[i].pos_seq1 + h->block_match[i].length;
	p2 = h->block_match[i].pos_seq2 + h->block_match[i].length;
    }

    /* Last block, only anchored at start */
    l1 = h->seq1_len - p1;
    l2 = h->seq2_len - p2;

    match += MIN(l1, l2) - (MIN(l1, l2)/h->word_length+1);
    mismatch += MIN(l1, l2)/h->word_length + 1;

    worst_mis += MIN(l1, l2);

    /*
    printf("len %d,%d match %d, mismatch %d, %5.1f%% worst %5.1f%%\n",
	   h->seq1_len, h->seq2_len, match, mismatch,
	   100.0*mismatch / (match + mismatch),
	   100.0*worst_mis / (worst_mat + worst_mis));
    */

    if (mat_p)
	*mat_p = match;
    if (mis_p)
	*mis_p = mismatch;

    return 100 * mismatch / (match + mismatch);
}

/*
 * Use the edit buffers to compute an overlap percentage mismatch
 */
void overlap_mismatch(Hash *h, EDIT_PAIR *edit_pair, OVERLAP *overlap,
		      char **seq1_p, char **seq2_p,
		      int *seq1_len_p, int *seq2_len_p,
		      int **S1_p, int **S2_p, int *n1_p, int *n2_p) {
    int i, j;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int *S1, *S2;
    int n1, n2;
    int mat, total;
    int score;
    int left1, left2, right1, right2;
    int op1, op2;

    seq1 = h->seq1;  seq1_len = h->seq1_len;
    seq2 = h->seq2;  seq2_len = h->seq2_len;
    S1 = edit_pair->S1;
    S2 = edit_pair->S2;
    n1 = edit_pair->next1;
    n2 = edit_pair->next2;

    /* Step 1: trim alignment edges */

    /* AGGAGGT... S1
     * ****GGT... S2
     */
    left1 = left2 = 0;
    if (S2[0] < 0) {
	int d = -S2[0];
	S2++; n2--;
	left2 = d;

	/* Consume first 'd' characters in seq1 */
	while (d > 0) {
	    if (S1[0] > d) {
		seq1 += d; seq1_len -= d;
		S1[0] -= d;
		d = 0;
	    } else if (S1[0] == d) {
		seq1 += d; seq1_len -= d;
		S1++; n1--;
		d = 0;
	    } else if (S1[0] >= 0) {
		seq1 += S1[0]; seq1_len -= S1[0];
		d -= S1[0];
		S1++; n1--;
	    } else { /* S1 < 0, both padded */
		d -= -S1[0];
		S1++; n1--;
	    }
	}

    /* ****GGT... S1
     * AGGAGGT... S2
     */
    } else if (S1[0] < 0) {
	int d = -S1[0];
	S1++; n1--;
	left1 = d;

	/* Consume first 'd' characters in seq2 */
	while (d > 0) {
	    if (S2[0] > d) {
		seq2 += d; seq2_len -= d;
		S2[0] -= d;
		d = 0;
	    } else if (S2[0] == d) {
		seq2 += d; seq2_len -= d;
		S2++; n2--;
		d = 0;
	    } else if (S2[0] >= 0) {
		seq2 += S2[0]; seq2_len -= S2[0]; left2 += S2[0];
		d -= S2[0];
		S2++; n2--;
	    } else { /* S2 < 0 */
		d -= -S2[0];
		S2++; n2--;
	    }
	}
    }

    right1 = right2 = MAX(left1, left2);

    /* ...GGTAGAC S1
     * ...GGT**** S2
     */
    if (S2[n2-1] < 0) {
	int d = -S2[n2-1];
	n2--;
	right1 += d;

	/* Consume last 'd' characters in seq1 */
	while (d > 0) { 
	    if (S1[n1-1] > d) {
		S1[n1-1] -= d;
		seq1_len -= d;
		d = 0;
	    } else if (S1[n1-1] == d) {
		seq1_len -= d;
		d = 0;
		n1--;
	    } else if (S1[n1-1] >= 0) {
		seq1_len -= S1[n1-1];
		d -= S1[n1-1];
		n1--;
	    } else { /* S1 < 0 */
		d -= -S1[n1-1];
		n1--;
	    }
	}

    /* ...GGTAGAC S1
     * ...GGT**** S2
     */
    } else if (S1[n1-1] < 0) {
	int d = -S1[n1-1];
	n1--;
	right2 += d;

	/* Consume last 'd' characters in seq2 */
	while (d > 0) { 
	    if (S2[n2-1] > d) {
		S2[n2-1] -= d;
		seq2_len -= d;
		d = 0;
	    } else if (S2[n2-1] == d) {
		seq2_len -= d;
		d = 0;
		n2--;
	    } else if (S2[n2-1] >= 0) {
		seq2_len -= S2[n2-1];
		d -= S2[n2-1];
		n2--;
	    } else { /* S2 < 0 */
		d -= -S2[n2-1];
		n2--;
	    }
	}
    }

    if (seq1_p)     *seq1_p = seq1;
    if (seq2_p)     *seq2_p = seq2;
    if (seq1_len_p) *seq1_len_p = seq1_len;
    if (seq2_len_p) *seq2_len_p = seq2_len;
    if (S1_p)       *S1_p = S1;
    if (S2_p)       *S2_p = S2;
    if (n1_p)       *n1_p = n1;
    if (n2_p)       *n2_p = n2;

    /* i= print_alignment(seq1, seq2, seq1_len, seq2_len, S1, S2, n1, n2, 100, stdout); */

    /* Now compute percent identity */
    i = j = 0;
    total = mat = score = 0;
    op1 = op2 = 0;
    while (i < seq1_len && j < seq2_len) {
	char c1, c2;

	while (op1 == 0) {
	    op1 = *S1++;
	}

	if (op1 < 0) {
	    c1 = '*';
	    op1++;
	} else if (op1 > 0) {
	    c1 = seq1[i++];
	    op1--;
	}

	while (op2 == 0) {
	    op2 = *S2++;
	}

	if (op2 < 0) {
	    c2 = '*';
	    op2++;
	} else if (op2 > 0) {
	    c2 = seq2[j++];
	    op2--;
	}

	right1++; right2++;
	if (c1 == c2) {
	    mat++;
	    score++;
	} else {
	    score -= 4;
	}
	total++;
    }
    right1--; right2--;

    overlap->percent = total ? 100.0 * mat/total : 0;
    overlap->length = total;
    overlap->qual = overlap->score = score;

    overlap->left1  = left1;
    overlap->right1 = right1;
    overlap->left2  = left2;
    overlap->right2 = right2;

    overlap->left = MAX(left1, left2);
    overlap->right = MIN(right1, right2);

    /* Figure out other overlap fields: direction, lo and ro */
    /* See seq_to_overlap() from align_lib.c */
    if(overlap->left1 == overlap->left2)
	overlap->direction = (overlap->right1 >= overlap->right2) ? 2 : 3;
    else if(overlap->left1 < overlap->left2)
	overlap->direction = (overlap->right1 >= overlap->right2) ? 2 : 0;
    else
	overlap->direction = (overlap->right1 <= overlap->right2) ? 3 : 1;

    /*
     * Calculate the offsets of the alignment. (i.e. the lengths of
     * the overhangs at each end.)
     *
     * Currently, if the overlap goes from a to b, then
     *    left_offset  = b.left  - a.left
     *    right_offset = b.right - a.right
     *
     * so if a containment overlap will have +ve left_offset and
     * -ve right_offset, and non-containment overlaps will have
     * both +ve.
     *
     * N.B. 'a' is not necessarily seq1 and 'b' is not necessarily
     * seq2, this depends on the direction of the overlap.
     * FIXME - maybe they should be, so that overlap info need never
     * be changed even if the readings are complemented.
     */
    switch(overlap->direction) {
    case 0: case 2:
	overlap->lo = overlap->left2 - overlap->left1;
	overlap->ro = overlap->right2 - overlap->right1;
	break;
    case 1: case 3:
	overlap->lo = overlap->left1 - overlap->left2;
	overlap->ro = overlap->right1 - overlap->right2;
	break;
    default:
	break;
    }

    assert(overlap->length == overlap->right - overlap->left + 1);
}

int align_wrap ( Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap_out) {
    int edge_mode = params->edge_mode;
    int i, s1, s2;
    int band, band_in;
    OVERLAP *overlap;
    EDIT_PAIR *edit_pair;
    int max_edit_pair;
    int max_seq;
    char NEW_PAD_SYM, OLD_PAD_SYM;
    /* int mis, mat; */

    /* Returned from overlap_mismatch */
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int *S1, *S2;
    int n1, n2;

    /*
    min_mismatch(h, &mis, &mat);
    printf("mis/mat = %d/%d = %f\n", mis, mat, 1.0*mis/(mis+mat));
    if (100.0*mis / (mis+mat) > 10)
	return -1;
    */

    NEW_PAD_SYM = params->new_pad_sym;
    OLD_PAD_SYM = params->old_pad_sym;

    /* input h which contains a list of matching blocks
     * output overlap_out which contains an overlap structure
     * ie a pair of edit buffers and a pair of aligned seqs
     *
     * strategy is to build up a pair of overall edit buffers
     * using the blocks and affine_align where necessary in between.
     * in some gaps between blocks (where 1 gap is length 0) we
     * simply add to the edit buffers. The edit buffers are stored
     * in an edit_pair structure and their contents copied to
     * overlap_out at the end. After each call to affine_align
     * we copy its edit buffers to the edit_pair structure.
     * 3 phases:
     * 1. up to first block (if there is a mismatched region)
     * 2. between blocks
     * 3. from last block to end (if there is a mismatched region)
     *
     * having built up the edit buffers, at the end we use them to
     * create a sequence alignment in order to fill out the other
     * bits of the overlap structure (would have prefered to be
     * able to avoid this and optionally only return the edit buffers,
     * but that requires extra work)
     * The alternative, which now looks more attractive, given that
     * we plan to use aligned sequences to direct editing, is to
     * rewrite all this stuff to create a pair of sequence alignments
     * as we go along, and not use edit buffers at all!!
     *
     * functions:
     * align_bit
     * sent a segment; either updates edit_pairs directly or uses
     * affine_align and then copies its edit buffers into edit_pair.
     * update_edit_pair
     * copies the edit buffers from an overlap into the edit_pair
     * frees the overlap edit buffers
     * block_to_edit_pair
     * adds a matching block to the edit_pair
     *
     * essential to make sure that the affine_align routine uses
     * the correct edge_mode
     * At the end we destroy all temporary memory: overlap and edit_pair
     */

    band = 0;
    band_in = params->band;
    max_edit_pair = MIN(h->seq1_len,h->seq2_len);

    if (NULL == (edit_pair = create_edit_pair(max_edit_pair))) {
	return -1;
    }
    if (NULL == (overlap = create_overlap())) {
	destroy_edit_pair(edit_pair);
	return -1;
    }

    init_overlap (overlap, h->seq1, h->seq2, h->seq1_len, h->seq2_len);
 
    /* align up to the first matching words,
     * align the segments between the matching words
     * align from the last matching words to the ends
     */

    /* get a start point based on the first matching word positions */

    diagonal_intercepts (h->block_match[0].diag,h->seq1_len,h->seq2_len,
			 &s1, &s2);

    overlap->seq1_len = h->block_match[0].pos_seq1;
    overlap->seq2_len = h->block_match[0].pos_seq2;
    overlap->seq1 = h->seq1;
    overlap->seq2 = h->seq2;
    
    /*
    printf("<Align pos %d+%d / %d+%d (%.10s... %.10s...)\n",
	   0, overlap->seq1_len, 0, overlap->seq2_len,
	   overlap->seq1, overlap->seq2);
    */

    params->edge_mode = (edge_mode & ~BEST_EDGE_TRACE) | FULL_LENGTH_TRACE;
    params->edge_mode &= ~EDGE_GAPS_COUNT;
    params->edge_mode |=  EDGE_GAPS_ZERO;
    if (band_in) {
	if (h->fast_mode)
	    band = set_band_blocks_fast(overlap->seq1_len,overlap->seq2_len);
	else
	    band = set_band_blocks(overlap->seq1_len,overlap->seq2_len);
    }
    set_align_params (params, band, 0,0,0,0, s1, s2,0,0,1);

    if (align_bit ( params, overlap, edit_pair)) {
 	verror(ERR_WARN, "align_wrap", "failed in align_bit");
	destroy_edit_pair(edit_pair);
	destroy_overlap(overlap);
	return -1;
    }
    free_overlap(overlap);

    /*printf("if required, done up to block 0\n");*/

    if ( block_to_edit_pair ( edit_pair, h->block_match[0].length)) {
	verror(ERR_WARN, "align_wrap", "failed in block_to_edit_pair");
	destroy_edit_pair(edit_pair);
	destroy_overlap(overlap);
	return -1;
    }
    s1 = h->block_match[0].pos_seq1 + h->block_match[0].length;
    s2 = h->block_match[0].pos_seq2 + h->block_match[0].length;
    params->edge_mode = EDGE_GAPS_COUNT | FULL_LENGTH_TRACE;
    params->edge_mode &= ~EDGE_GAPS_ZERO;
    for(i=1;i<h->matches;i++) {
	overlap->seq1_len = h->block_match[i].pos_seq1 - s1;
	overlap->seq2_len = h->block_match[i].pos_seq2 - s2;
	overlap->seq1 = &(h->seq1[s1]);
	overlap->seq2 = &(h->seq2[s2]);

	if ( MAX(overlap->seq1_len,overlap->seq2_len) > 0 ) {

	    if (band_in) {
		if (h->fast_mode)
		    band = set_band_blocks_fast(overlap->seq1_len,
						overlap->seq2_len);
		else
		    band = set_band_blocks(overlap->seq1_len,
					   overlap->seq2_len);
	    }
	    set_align_params (params, band, 0,0,0,0,0,0,0,0,1);

	    if (align_bit ( params, overlap, edit_pair)) {
		verror(ERR_WARN, "align_wrap", "failed in align_bit");
		destroy_edit_pair(edit_pair);
		destroy_overlap(overlap);
		return -1;
	    }
	    free_overlap(overlap);

	}
	s1 = h->block_match[i].pos_seq1 + h->block_match[i].length;
	s2 = h->block_match[i].pos_seq2 + h->block_match[i].length;
	if ( block_to_edit_pair ( edit_pair, h->block_match[i].length)) {
	    verror(ERR_WARN, "align_wrap", "failed in block_to_edit_pair");
	    destroy_edit_pair(edit_pair);
	    destroy_overlap(overlap);
	    return -1;
	}
    }
    /* do up to the end */
    
    overlap->seq1_len = h->seq1_len - s1;
    overlap->seq2_len = h->seq2_len - s2;
    overlap->seq1 = &(h->seq1[s1]);
    overlap->seq2 = &(h->seq2[s2]);
    
    if (band_in) {
	if (h->fast_mode)
	    band = set_band_blocks_fast(overlap->seq1_len,overlap->seq2_len);
	else
	    band = set_band_blocks(overlap->seq1_len,overlap->seq2_len);
    }
    set_align_params (params, band, 0,0,0,0, 0,0,0,0,1);
    params->edge_mode = (edge_mode & ~EDGE_GAPS_ZERO) | EDGE_GAPS_COUNT;
    params->edge_mode &= ~FULL_LENGTH_TRACE;
    params->edge_mode |= BEST_EDGE_TRACE;
    if (align_bit ( params, overlap, edit_pair)) {
	verror(ERR_WARN, "align_wrap", "failed in align_bit");
	destroy_edit_pair(edit_pair);
	destroy_overlap(overlap);
	return -1;
    }
    destroy_overlap(overlap);

    /* Compute overlap percentage mismatch without allocating huge buffers */
    if (!(params->job & RETURN_END_GAPS)) {
	overlap_mismatch(h, edit_pair, overlap_out,
			 &seq1, &seq2, &seq1_len, &seq2_len,
			 &S1, &S2, &n1, &n2);
	destroy_edit_pair(edit_pair);
	return 0;
    }

    seq1 = h->seq1;  seq1_len = h->seq1_len;
    seq2 = h->seq2;  seq2_len = h->seq2_len;
    S1 = edit_pair->S1;
    S2 = edit_pair->S2;
    n1 = edit_pair->next1;
    n2 = edit_pair->next2;


    /*
    for(i=0;i<edit_pair->next1;i++) printf("1 %d %d\n",i,edit_pair->S1[i]);
    for(i=0;i<edit_pair->next2;i++) printf("2 %d %d\n",i,edit_pair->S2[i]);
    
    i= print_alignment(h->seq1,
		    h->seq2,
		    h->seq1_len,
		    h->seq2_len,
		    edit_pair->S1,
		    edit_pair->S2,
		    edit_pair->next1,
		    edit_pair->next2,
		    100,
		    stdout);
    */
    

    /* now create the overlap from the edit buffers */
    max_seq = h->seq1_len + h->seq2_len + 1;

    /*printf("align_wrap malloc overlap_out %d\n",max_seq);*/
    if(!(overlap_out->seq1_out = (char *) xmalloc(sizeof(char) * 
						      max_seq))) {
	verror(ERR_WARN, "align_wrap", "malloc failed for seq1_out");
	destroy_edit_pair(edit_pair);
	return -1;
    }
    if(!(overlap_out->seq2_out = (char *) xmalloc(sizeof(char) * 
						      max_seq))) {
	verror(ERR_WARN, "align_wrap", "malloc failed for seq2_out");
	destroy_edit_pair(edit_pair);
	return -1;
    }

    seq_expand(seq1, overlap_out->seq1_out, &s1, S1, n1, 3, NEW_PAD_SYM);
    seq_expand(seq2, overlap_out->seq2_out, &s2, S2, n2, 3, NEW_PAD_SYM);

    overlap_out->seq_out_len = s1;

    /* save edit_pairs in overlap_out FIXME if we need them*/
    if(!(overlap_out->S1 = (int *) xmalloc(sizeof(int) * n1))) {
	verror(ERR_WARN, "align_wrap", "malloc failed for S1");
	destroy_edit_pair(edit_pair);
	return -1;
    }
    if(!(overlap_out->S2 = (int *) xmalloc(sizeof(int) * n2))) {
	verror(ERR_WARN, "align_wrap", "malloc failed for S2");
	destroy_edit_pair(edit_pair);
	return -1;
    }

    memcpy(overlap_out->S1, S1, n1*sizeof(*S1));
    memcpy(overlap_out->S2, S2, n2*sizeof(*S2));
    overlap_out->s1_len = n1;
    overlap_out->s2_len = n2;
    destroy_edit_pair(edit_pair);

    /* when we enter seq_to_overlap from here the overlap score is not set */
    overlap_out->score = 0; /* assigned by seq_to_overlap */
    if( 0 != seq_to_overlap(overlap_out,OLD_PAD_SYM,NEW_PAD_SYM)) {
	return -1;
    }
    if ( params->job & RETURN_NEW_PADS ) {
	old_pads_for_new(overlap_out->seq1_out,overlap_out->seq_out_len,
			 OLD_PAD_SYM,NEW_PAD_SYM);
	old_pads_for_new(overlap_out->seq2_out,overlap_out->seq_out_len,
			 OLD_PAD_SYM,NEW_PAD_SYM);
    }
    /*overlap_out->score = overlap_out->percent; */
    overlap_out->qual  = overlap_out->percent;

    return 0;
}

int central_diagonal ( Hash *h ) {
    int i, j;
    if(h->matches == 0) return 0;
    for(i=0,j=0;i<h->matches;i++) j+= h->block_match[i].diag;
    return j/h->matches;
}

static int c1_len, c2_len;

static int sort_func(const void *p1, const void *p2) {
    int x1,y1,x2,y2;
    Block_Match *c1 = (Block_Match *)p1;
    Block_Match *c2 = (Block_Match *)p2;

    /*
    x1 = MIN(c1->pos_seq1, c1_len - (c1->pos_seq1 + c1->length));
    y1 = MIN(c1->pos_seq2, c2_len - (c1->pos_seq2 + c1->length));
    x2 = MIN(c2->pos_seq1, c1_len - (c2->pos_seq1 + c2->length));
    y2 = MIN(c2->pos_seq2, c2_len - (c2->pos_seq2 + c2->length));
    */
    x1 = c1->pos_seq1;
    y1 = c1->pos_seq2;
    x2 = c2->pos_seq1;
    y2 = c2->pos_seq2;

    return (x1+y1) - (x2+y2);
}

int sort_blocks ( Block_Match *block_match, int matches ) {
    qsort ((void *) block_match, matches, sizeof(Block_Match), sort_func);
    return 0;
}

static int sort_pos1_func(const void *p1, const void *p2) {
    Block_Match *c1 = (Block_Match *)p1;
    Block_Match *c2 = (Block_Match *)p2;

    return c1->pos_seq1 - c2->pos_seq1;
}

int sort_pos1_blocks ( Block_Match *block_match, int matches ) {
    qsort ((void *) block_match, matches, sizeof(Block_Match), sort_pos1_func);
    return 0;
}

static int sort_len_func(const void *p1, const void *p2) {
    int x1,x2;
    Block_Match *c1 = (Block_Match *)p1;
    Block_Match *c2 = (Block_Match *)p2;

    x1 = c1->length;
    x2 = c2->length;
    return (x2-x1);
}

int sort_len_blocks ( Block_Match *block_match, int matches ) {
    qsort ((void *) block_match, matches, sizeof(Block_Match), sort_len_func);
    return 0;
}

static int sort_by_end_pos(const void *p1, const void *p2) {
    int x1,y1,x2,y2,l1,l2;
    Block_Match *c1 = (Block_Match *)p1;
    Block_Match *c2 = (Block_Match *)p2;

    x1 = c1->pos_seq1;
    y1 = c1->pos_seq2;
    l1 = c1->length;

    x2 = c2->pos_seq1;
    y2 = c2->pos_seq2;
    l2 = c2->length;

    return (x1+y1+l1+l1) - (x2+y2+l2+l2);
}
/* #define DEBUG_ALIGN_BLOCKS */
/* Returns 1 on success, <= 0 on failure */
int align_blocks ( Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap ) {
    int i,j,k,gap_pen,diag_shift,best_score,best_prev,t,tt;
    int good_blocks;
    int *index_ptr = NULL;
    double best_percent;
    int more_shuffling;

    best_score = -1000000;
    best_prev = -1;

    //UpdateTextOutput();

#ifdef DEBUG_ALIGN_BLOCKS
    printf("=== STEP 0 === %d\n",h->matches);
    for (i=0;i<h->matches;i++) {
	printf("i %d %d %d %d %d %d %d\n",i,
	       h->block_match[i].pos_seq1,
	       h->block_match[i].pos_seq2,
	       h->block_match[i].length,
	       h->block_match[i].diag,
	       h->block_match[i].best_score,
	       h->block_match[i].prev_block);
    }
    fflush(stdout);
#endif

    if ( h->matches < 1 ) return 0;

#if 0
    /* sort the blocks on length and then shrink the list
     * so that the sum of lengths is length of sequences
     */
    i = sort_len_blocks(h->block_match, h->matches);
    l = MIN(h->seq1_len,h->seq2_len);
    for (i=0,j=0; i<h->matches; i++) {
	j+=h->block_match[i].length;
	if (j>l) {
	    h->matches = i+1;
	    break;
	}
    }
#endif

    /* set each blocks score to its distance from the nearest edge
     * find the best score as this start score plus match length
     * and note the block number
     */

    for (i=0;i<h->matches;i++) {
        gap_pen = -MIN(h->block_match[i].pos_seq1, h->block_match[i].pos_seq2);
        if((t=h->block_match[i].length + gap_pen) > best_score) {
            best_score = t;
            best_prev = i;
        }

	h->block_match[i].best_score = gap_pen;
	h->block_match[i].prev_block = -1;
    }

    if (best_prev == -1) return 0; /* bail out if there is no good score */

#ifdef DEBUG_ALIGN_BLOCKS
    printf("=== STEP 1, best_score = %d,%d === %d\n",
	   best_score, best_prev, h->matches);
    for (i=0;i<h->matches;i++) {
	printf("i %d %d %d %d %d %d %d\n",i,
	       h->block_match[i].pos_seq1,
	       h->block_match[i].pos_seq2,
	       h->block_match[i].length,
	       h->block_match[i].diag,
	       h->block_match[i].best_score,
	       h->block_match[i].prev_block);
    }
    fflush(stdout);
#endif

    qsort(h->block_match, h->matches, sizeof(Block_Match), sort_by_end_pos);

#ifdef DEBUG_ALIGN_BLOCKS
    printf("=== STEP 1.5, best_score = %d,%d === %d\n",
	   best_score, best_prev, h->matches);
    for (i=0;i<h->matches;i++) {
	printf("i %d %d %d %d %d %d %d\n",i,
	       h->block_match[i].pos_seq1,
	       h->block_match[i].pos_seq2,
	       h->block_match[i].length,
	       h->block_match[i].diag,
	       h->block_match[i].best_score,
	       h->block_match[i].prev_block);
    }
    fflush(stdout);
#endif

    /* for each block look at all the previous ones to find its best
     * predecessor. This is given by: best_score + length - diag_shift
     * note the best score and block number as we proceed
     */
    for (i=1;i<h->matches;i++) {
	/* k is a counter to short-cut long searches.  This may lead to
	   sub-optimal routes being found, but in practice this rarely
	   appears to happen for good matches */
	k = 10;
	for(j=i-1; j>-1 && k; j--) {
	    int ovr;

	    /*
	     * This code aims to only stitch together blocks that
	     * overlap by <= MAX_BLOCK_OVERLAP.
	     * Fundamentallty it's the wrong thing to do though. If we limit
	     * this anywhere it should be in the computed diag_shift
	     * (consider the case of 1kb of TA repeat shifted by 20 bases,
	     * giving 980k overlap by 20 base diag).
	     *
	     * For now let's just assume it comes out in the score and
	     * always compute the score rather than pre-filtering.
	     */
	    /*
	    # define MAX_BLOCK_OVERLAP 20
	    if( ((h->block_match[j].pos_seq1 + h->block_match[j].length)
		 <= h->block_match[i].pos_seq1 + MAX_BLOCK_OVERLAP)
		&& ((h->block_match[j].pos_seq2 + h->block_match[j].length)
		 <= h->block_match[i].pos_seq2) + MAX_BLOCK_OVERLAP) {
	    */
	    if (1) {
		int len = h->block_match[j].length;
		int dist;

		diag_shift = abs(h->block_match[i].diag - 
				 h->block_match[j].diag);

		/* Reduce effective length when blocks overlap */
		ovr = 0;
		dist = h->block_match[j].pos_seq1 + h->block_match[j].length -
		    h->block_match[i].pos_seq1;
		if (dist > 0) {
		    ovr = 1;
		    len -= dist;
		}

		dist = h->block_match[j].pos_seq2 + h->block_match[j].length -
		    h->block_match[i].pos_seq2;
		if (dist > 0) {
		    ovr = 1;
		    len -= dist;
		}

		if (!ovr)
		    k--;

		/* best score does not include current match */
		t = h->block_match[j].best_score + len - diag_shift;
		/*
		 * Also compensate for the fact that stitching together
		 * blocks inherently implies we have mismatches between
		 * blocks. At the very least we have 1 per word length
		 * (and 1 per diag shift, already compensated for above).
		 */
		t -= 2*MIN(h->block_match[i].pos_seq1 -
			   (h->block_match[j].pos_seq1+
			    h->block_match[j].length),
			   h->block_match[i].pos_seq2 -
			   (h->block_match[j].pos_seq2+
			    h->block_match[j].length)) / h->min_match;
		
		if ( t > h->block_match[i].best_score ) {
		    h->block_match[i].best_score = t;
		    h->block_match[i].prev_block = j;
		}
	    }
	}
    }

#ifdef DEBUG_ALIGN_BLOCKS
    printf("=== STEP 2 === %d, best_prev=%d\n",h->matches, best_prev);
    for (i=0;i<h->matches;i++) {
	printf("i %d %d %d %d %d %d %d\n",i,
	       h->block_match[i].pos_seq1,
	       h->block_match[i].pos_seq2,
	       h->block_match[i].length,
	       h->block_match[i].diag,
	       h->block_match[i].best_score,
	       h->block_match[i].prev_block);
    }
    fflush(stdout);
#endif

    /* Find which chain gives the best score to the edge */
    best_score = INT_MIN;
    for (i=0;i<h->matches;i++) {
	Block_Match *m = &h->block_match[i];
	gap_pen = -MIN(h->seq1_len - m->pos_seq1 - m->length,
		       h->seq2_len - m->pos_seq2 - m->length);
	if (h->block_match[i].best_score + gap_pen > best_score) {
	    best_score = h->block_match[i].best_score + gap_pen;
	    best_prev = i;
	}
    }
#ifdef DEBUG_ALIGN_BLOCKS
    printf("=== STEP 3 === %d, best_score = %d, best_prev=%d\n",
	   h->matches, best_score, best_prev);
#endif
    /* best_prev is now last block */

    /* shuffle the ordered blocks to the start of the array */

    tt = h->block_match[best_prev].best_score;
    h->block_match[best_prev].best_score = -1;

    good_blocks = 1;
    for (i=best_prev,j=0;i>-1;) {
	j = h->block_match[i].prev_block;
	if (j>-1) {
	    good_blocks++;
	}
	i=h->block_match[i].prev_block;
    }

    if ( ! ( index_ptr = (int *) xmalloc ( sizeof(int)*(good_blocks) ))) {
 	return -1;
    }

    for (i=best_prev,j=good_blocks-1;i>-1;j--) {
	index_ptr[j] = i;
	i=h->block_match[i].prev_block;
    }

    h->block_match[best_prev].best_score = tt;

    for (j=0;j<good_blocks;j++) {
	i = index_ptr[j];
	if (i != j) {
	    h->block_match[j].pos_seq1   = h->block_match[i].pos_seq1;
	    h->block_match[j].pos_seq2   = h->block_match[i].pos_seq2;
	    h->block_match[j].length     = h->block_match[i].length;
	    h->block_match[j].diag       = h->block_match[i].diag;
	    h->block_match[j].best_score = h->block_match[i].best_score;
	    h->block_match[j].prev_block = h->block_match[i].prev_block;
	}
    }

#ifdef DEBUG_ALIGN_BLOCKS
    printf("=== STEP 5 === %d, %d\n",h->matches, good_blocks);
    for (i=0;i<good_blocks;i++) {
	printf("i %d %d %d %d %d %d %d\n",i,
	       h->block_match[i].pos_seq1,
	       h->block_match[i].pos_seq2,
	       h->block_match[i].length,
	       h->block_match[i].diag,
	       h->block_match[i].best_score,
	       h->block_match[i].prev_block);
    }
    fflush(stdout);
#endif

    /*
     * Clip back overlapping block elements, keeping the longest block intact
     * where possible.
     */
    do {
	more_shuffling = 0;
	for (i = 0; i < good_blocks-1; i++) {
	    int dist;
	    dist = h->block_match[i].pos_seq1 + h->block_match[i].length
		- h->block_match[i+1].pos_seq1;
	    if (dist > 0) {
		if (h->block_match[i].length > h->block_match[i+1].length) {
		    h->block_match[i+1].length   -= dist;
		    h->block_match[i+1].pos_seq1 += dist;
		    h->block_match[i+1].pos_seq2 += dist;
		} else {
		    h->block_match[i].length -= dist;
		}
	    }

	    dist = h->block_match[i].pos_seq2 + h->block_match[i].length
		- h->block_match[i+1].pos_seq2;
	    if (dist > 0) {
		if (h->block_match[i].length > h->block_match[i+1].length) {
		    h->block_match[i+1].length   -= dist;
		    h->block_match[i+1].pos_seq1 += dist;
		    h->block_match[i+1].pos_seq2 += dist;
		} else {
		    h->block_match[i].length -= dist;
		}
	    }
	}

	/* Check for blocks with length <= 0 */
	for (i = j = 0; i < good_blocks; i++) {
	    if (h->block_match[i].length > 0) {
		h->block_match[j].pos_seq1   = h->block_match[i].pos_seq1;
		h->block_match[j].pos_seq2   = h->block_match[i].pos_seq2;
		h->block_match[j].length     = h->block_match[i].length;
		h->block_match[j].diag       = h->block_match[i].diag;
		h->block_match[j].best_score = h->block_match[i].best_score;
		h->block_match[j].prev_block = h->block_match[i].prev_block;
		j++;
	    } else {
		more_shuffling = 1;
	    }
	}
	good_blocks = j;

    } while (more_shuffling);

#ifdef DEBUG_ALIGN_BLOCKS
    printf("=== STEP 6 === %d, %d\n",h->matches, good_blocks);
    for (i=0;i<good_blocks;i++) {
	printf("i %d %d %d %d %d %d %d\n",i,
	       h->block_match[i].pos_seq1,
	       h->block_match[i].pos_seq2,
	       h->block_match[i].length,
	       h->block_match[i].diag,
	       h->block_match[i].best_score,
	       h->block_match[i].prev_block);
    }
    fflush(stdout);
#endif

    if ( index_ptr ) xfree (index_ptr);
    h->matches = good_blocks;
    /*printf("returned %d matches with best score %d\n",h->matches,best_score);*/

    tt = 0;
    for (i=0;i<h->matches;i++) {
	tt += h->block_match[i].length;
    }

    i = h->matches/2;
    i = h->block_match[i].diag;
    j = diagonal_length(h->seq1_len, h->seq2_len, i);

    best_percent = 100.0 * tt / j;

    overlap->seq1 = h->seq1;
    overlap->seq2 = h->seq2;
    overlap->seq1_len = h->seq1_len;
    overlap->seq2_len = h->seq2_len;

    /* FAST MODE; Check for diagonals too far apart */
    if (h->fast_mode) {
	for (i = 0; i < h->matches-1; i++) {
	    if (ABS(h->block_match[i].diag - h->block_match[i+1].diag) > 50)
		return 0;
	}
    }

#define MAX_G 1000
    /* FAST MODE; no more than MAX_G bp without a hash match */
    if (h->fast_mode) {
	int b_st1, b_en1, b_st2, b_en2;

	for (i = 0; i < h->matches-1; i++) {
	    b_st1 = h->block_match[i+1].pos_seq1;
	    b_en1 = h->block_match[i].pos_seq1 + h->block_match[i].length;
	    if (b_st1 - b_en1 > MAX_G)
		return 0;

	    b_st2 = h->block_match[i+1].pos_seq2;
	    b_en2 = h->block_match[i].pos_seq2 + h->block_match[i].length;

	    if (b_st2 - b_en2 > MAX_G)
		return 0;
	}

	b_st1 = h->block_match[0].pos_seq1;
	b_st2 = h->block_match[0].pos_seq2;
	if (MIN(b_st1, b_st2) > MAX_G)
	    return 0;

	b_en1 = h->block_match[h->matches-1].pos_seq1 +
	    h->block_match[h->matches-1].length;
	b_en2 = h->block_match[h->matches-1].pos_seq2 +
	    h->block_match[h->matches-1].length;
	if (MIN(h->seq1_len - b_en1, h->seq2_len - b_en2) > MAX_G)
	    return 0;
    }


    /* FAST MODE */
    if ((h->fast_mode && best_percent > 30.0) || best_percent > 10.0) {
	if (align_wrap (h,params,overlap) < 0)
	    return -1;

	return 1;
    }

    return 0;
}

int align_blocks_bulk(Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap,
		      int cnum,
		      Contig_parms *contig_list, int number_of_contigs,
		      void (*add_func)(OVERLAP *overlap,
				       int cnum1,
				       int cnum2,
				       void *clientdata),
		      void *add_data) {
    int i,j,k;
    int max_mat = 0;

    /*
    printf("=== STEP -1 === %d\n",h->matches);
    for (i=0;i<h->matches;i++) {
	printf("i %d %d %d %d %d %d %d\n",i,
	       h->block_match[i].pos_seq1,
	       h->block_match[i].pos_seq2,
	       h->block_match[i].length,
	       h->block_match[i].diag,
	       h->block_match[i].best_score,
	       h->block_match[i].prev_block);
    }
    */

    if ( h->matches < 1 ) return 0;

    /*
     * Sort the blocks by position and collate into contigs
     */
    sort_pos1_blocks(h->block_match, h->matches);

    for (i = j = 0; i < h->matches; i++) {
	while (h->block_match[i].pos_seq1 > contig_list[j].contig_end_offset)
	    j++;
	h->block_match[i].contig1 = j;
    }

    /*
    printf("=== STEP 0.1 === %d\n",h->matches);
    for (i=0;i<h->matches;i++) {
	printf("i %d %d %d %d %d %d %d   contig %d\n",i,
	       h->block_match[i].pos_seq1,
	       h->block_match[i].pos_seq2,
	       h->block_match[i].length,
	       h->block_match[i].diag,
	       h->block_match[i].best_score,
	       h->block_match[i].prev_block,
	       h->block_match[i].contig1);
	if (i+1 == h->matches ||
	    h->block_match[i+1].contig1 != h->block_match[i].contig1) {
	    printf("--%d--\n", h->block_match[i].contig1);
	}
    }
    */

    /* Create a hash & overlap struct for matches from just one contig */
    max_mat = h->block_match[0].length;
    for (j=i=0;i<h->matches;i++) {
	if (i+1 == h->matches ||
	    h->block_match[i+1].contig1 != h->block_match[i].contig1) {
	    Hash h1;
	    OVERLAP o1;
	    int c = h->block_match[i].contig1;

	    if (max_mat < h->min_match) {
		goto next;
	    }

	    memcpy(&h1, h, sizeof(h1));
	    memcpy(&o1, overlap, sizeof(o1));

	    h1.seq1 = &h->seq1[contig_list[c].contig_start_offset];
	    h1.seq1_len = contig_list[c].contig_end_offset -
		contig_list[c].contig_start_offset + 1;

	    h1.matches = i+1-j;
	    /*
	    h1.block_match = malloc(h1.matches * sizeof(*h1.block_match));
	    for (k = 0; k < h1.matches; k++) {
		h1.block_match[k] = h->block_match[k+j];
		h1.block_match[k].pos_seq1 -= contig_list[c].contig_start_offset;
		h1.block_match[k].diag =
		    h1.seq1_len -
		    h1.block_match[k].pos_seq1 +
		    h1.block_match[k].pos_seq2 - 1;
		    
	    }
	    */

	    h1.block_match = &h->block_match[j];
	    for (k = 0; k < h1.matches; k++) {
		h1.block_match[k].pos_seq1 -= contig_list[c].contig_start_offset;
		h1.block_match[k].diag =
		    h1.seq1_len -
		    h1.block_match[k].pos_seq1 +
		    h1.block_match[k].pos_seq2 - 1;
		    
	    }
	    

	    o1.seq1 = h1.seq1;
	    o1.seq2 = h1.seq2;
	    o1.seq1_len = h1.seq1_len;
	    o1.seq2_len = h1.seq2_len;

	    /* Call align_blocks as if we were doing only a pair-wise search */
	    /*
	    printf("Comparing rec %"PRIrec" vs %"PRIrec"\n",
		   contig_list[c].contig_number,
		   contig_list[cnum].contig_number);
	    */
	    if (align_blocks(&h1, params, &o1)) {
		add_func(&o1, c, cnum, add_data);
	    }

	next:

	    j = i+1;
	    max_mat = 0;
	}

	if (max_mat < h->block_match[i].length)
	    max_mat = h->block_match[i].length;
    }

    return 0;
}

int compare_seqs(Hash *h, int *seq1_match_pos, int *seq2_match_pos,
		int *match_length) {
    
    int	       	ncw, nrw, word, pw1, pw2, i, j, match_size;
    int	       	diag_pos, size_hist;
    
    if(h->seq1_len < h->min_match) return -4; 
    if(h->seq2_len < h->min_match) return -4; 

    size_hist = h->seq1_len + h->seq2_len - 1;
    for(i = 0; i < size_hist; i++) h->diag[i] = -(h->word_length);
    
    nrw = h->seq2_len - h->word_length + 1;
    
    /* 	loop for all (nrw) complete words in values2 */

    h->matches = -1;
    for (pw2=0;pw2<nrw;pw2++) {
 	word = h->values2[pw2];
	if ( -1 != word ) {
	    if ( 0 != (ncw = h->counts[word]) ) {
		for (j=0,pw1=h->last_word[word];j<ncw;j++) {
		    diag_pos = h->seq1_len - pw1 + pw2 - 1;
		    if ( h->diag[diag_pos] < pw2 ) {
			if ((match_size = match_len (h->seq1,
						     pw1, h->seq1_len,
						     h->seq2,
						     pw2, h->seq2_len))
			    >= h->min_match ) {
			    h->matches++;
			    if(h->max_matches == h->matches) {
				return -5;
			    }
			    seq1_match_pos[h->matches] = pw1+1;
			    seq2_match_pos[h->matches] = pw2+1;
			    match_length[h->matches] = match_size;
			}
			h->diag[diag_pos] = pw2 + match_size;
		    }
		    pw1 = h->values1[pw1];
		}
	    }
	}
    }
    h->matches += 1;
    return h->matches;
}



int compare_a(Hash *h,
	    ALIGN_PARAMS *params, OVERLAP *overlap) {
    
    /*
     * Takes two sequences and attempts to find word matches between
     * them by hashing. The hashing is done for a certain word length
     * (8), and then the program attempts to find all cases where the
     * word can be extended to >= min_match. The position of the longest
     * matching word found during this process is extrapolated back to
     * the x or y axis of the comparison matrix. This new position is
     * used to calculate the limits of the band used by the alignment
     * routine (NB. The alignment is done externally to this function.)
     */
    
    int	       	ncw, nrw, word, pw1, pw2, i, j, match_length;
    int	       	diag_pos, size_hist;
    int	       	hist_left, hist_right;
    
    int band, band_in;
    int		match_found;
    
    if(h->seq1_len < h->word_length) return -4; 
    if(h->seq2_len < h->word_length) return -4; 

    band_in = params->band;
    size_hist = h->seq1_len + h->seq2_len - 1;
    for(i = 0; i < size_hist; i++) h->diag[i] = -(h->word_length);
    for(i = 0; i < size_hist; i++) h->hist[i] = 0;
    
    nrw = h->seq2_len - h->word_length + 1;

    /* 	loop for all (nrw) complete words in values2 */
    match_found = 0;
    for (pw2=0;pw2<nrw;pw2++) {
 	word = h->values2[pw2];
	if ( -1 != word ) {
	    if ( 0 != (ncw = h->counts[word]) ) {
		for (j=0,pw1=h->last_word[word];j<ncw;j++) {
		    diag_pos = h->seq1_len - pw1 + pw2 - 1;
		    if ( h->diag[diag_pos] < pw2 ) {
			match_length = match_len (h->seq1, pw1, h->seq1_len,
						  h->seq2, pw2, h->seq2_len);
			/* the old way:
			if ((match_length = match_len ( 
						       h->seq1, pw1, h->seq1_len,
						       h->seq2, pw2, h->seq2_len))
			    >= h->min_match ) {
			    match_found = 1;
			}
			*/
			h->hist[diag_pos] += 1 + match_length - h->word_length;
			h->diag[diag_pos] = pw2 + match_length;
		    }
		    pw1 = h->values1[pw1];
		}
	    }
	}
    }
    match_found = 1;
    if(match_found) {
	/* now find the most significant diagonals */

	hist_left = MINMAT - 1;
	hist_right = size_hist - MINMAT;
	h->matches = -1;
	for(i = hist_left; i<hist_right; i++) {
	    j = diagonal_length(h->seq1_len, h->seq2_len, i); 
	    if(h->hist[i] > h->expected_scores[j]) { 
		h->matches++;
		/*printf("match %d %d %e %e\n",i,j,h->hist[i],h->expected_scores[j]);*/
		if(h->max_matches == h->matches) {
		    h->max_matches *= 2;
		    h->diag_match =
			(Diag_Match *)xrealloc(h->diag_match,
					       sizeof(Diag_Match) *
					       h->max_matches);
		    if (NULL == h->diag_match) {
			return -5;
		    }
		}
		h->diag_match[h->matches].pos = i;
		h->diag_match[h->matches].prob = h->hist[i] / (double)j;
		/*h->diag_match[h->matches].prob = h->hist[i] / h->expected_scores[j];*/
	    }
	}
	h->matches += 1;
	if (h->matches < 1) return 0;

	/* get the best intercept defined by these diagonals */

	if (!best_intercept ( h, &pw1, &pw2))
	  return 0;
	/*printf("%d %d\n",pw1,pw2);*/

	band = 0;
	if (band_in) {
	    double perc;
	    perc = (double)band_in / 100.0;
	    j = MIN(h->seq1_len+1-pw1,h->seq2_len+1-pw2);
	    band = MAX(20,(perc*j));
	    /*printf("band_in %d j %d perc %f band %d\n",band_in,j,perc,band);*/
	}
	set_align_params (params, band, 0,0,0,0, pw1, pw2,0,0,1);
	/*print_align_params(params);*/
	i = affine_align(overlap,params);
	params->band = band_in;
	if ( i ) return -1;
	/*print_overlap_struct(overlap);*/
	return 1;
    }
    
    return 0;
}

/*
static int fast_match_len(int w,
			  const char *s1, int p1, int l1,
			  const char *s2, int p2, int l2) {
    int l = w;

    p1 += w;
    p2 += w;
    if (l1-p1 < l2-p2) {
	while (p1 < l1 && s1[p1++] == s2[p2++]) {
	    l++;
	}
    } else {
	while (p2 < l2 && s1[p1++] == s2[p2++]) {
	    l++;
	}
    }

    return l;
}
*/

/*
 * Returns match length and *bck as the distance backwards to travel.
 */
static int match_fwd_back(int word_len,
			  const char *s1, int p1, int l1,
			  const char *s2, int p2, int l2,
			  int *bck) {
    int l, i1, i2;

    /* Backwards */
    i1 = p1-1;
    i2 = p2-1;
    if (i1 < i2) {
	while (i1 >= 0 && s1[i1] == s2[i2]) {
	    i1--;
	    i2--;
	}
    } else {
	while (i2 >= 0 && s1[i1] == s2[i2]) {
	    i1--;
	    i2--;
	}
    }
    l = p1-1 - i1;
    *bck = l;

    /* Fwds */
    i1 = p1 + word_len;
    i2 = p2 + word_len;
    if (l1-i1 < l2-i2) {
	while (i1 < l1 && s1[i1] == s2[i2])
	    i1++,i2++;
    } else {
	while (i2 < l2 && s1[i1] == s2[i2])
	    i1++,i2++;
    }
    l += i1 - p1;

    return l;
}

/*
 * Returns match length and *bck as the distance backwards to travel.
 *
 * As per match_fwd_back, but we allow single bp mismatches provided they're
 * followed by 2 bases of match. The reason is that we'll *probably* end up
 * going along this diagonal in the dynamic programming step anyway.
 *
 * I say probably: in the case of overlapping blocks on close diagonals we'd
 * need to trim one. Eg:
 *
 *              /                 /                     /
 *             /                 / 	               / 
 *        /   /	            /    	              /	 
 *       /   /    =>       /         or	             /   
 *      /   	          /   	 	        /   	 
 *     /                 /         	       /         
 *
 * Both are equally optimal alignments, but if one of our diagonals went
 * through a mismatch then we may make the wrong choice. It's faster though
 * so we're testing this for efficiencies sake.
 *
 * NOTE: this forces more alignments too as previously 2 20mers with a
 * mismatch would be filtered on a minmatch=25 param, whereas now we get a
 * 41mer which is accepted.
 */
#if 0 /* UNUSED currently */
static int match_fwd_back_mm(int word_len,
			     const char *s1, int p1, int l1,
			     const char *s2, int p2, int l2,
			     int *bck, int *exact) {
    int l, i1, i2;
    int mm = word_len, max_mm = word_len;

    /* Backwards */
    i1 = p1-1;
    i2 = p2-1;
    if (i1 < i2) {
	while (i1 >= 0) {
	    if (s1[i1] == s2[i2]) {
		i1--;
		i2--;
		mm++;
	    } else {
		if (i1 >= 2 && s1[i1-1] == s2[i2-1] && s1[i1-2] == s2[i2-2]) {
		    if (max_mm < mm)
			max_mm = mm;
		    mm = 0;
		    i1 -= 3;
		    i2 -= 3;
		} else {
		    break;
		}
	    }
	}
    } else {
	while (i2 >= 0) {
	    if (s1[i1] == s2[i2]) {
		i1--;
		i2--;
		mm++;
	    } else {
		if (i2 >= 2 && s1[i1-1] == s2[i2-1] && s1[i1-2] == s2[i2-2]) {
		    i1 -= 3;
		    i2 -= 3;
		    if (max_mm < mm)
			max_mm = mm;
		    mm = 0;
		} else {
		    break;
		}
	    }
	}
    }
    l = p1-1 - i1;
    *bck = l;

    if (max_mm < mm)
	max_mm = mm;

    /* Fwds */
    i1 = p1 + word_len;
    i2 = p2 + word_len;
    if (l1-i1 < l2-i2) {
	while (i1 < l1) {
	    if (s1[i1] == s2[i2]) {
		i1++,i2++;
		mm++;
	    } else {
		if (i1+2 < l1 && s1[i1+1]==s2[i2+1] && s1[i1+2]==s2[i2+2]) {
		    i1 += 3;
		    i2 += 3;
		    if (max_mm < mm)
			max_mm = mm;
		    mm = 0;
		} else {
		    break;
		}
	    }
	}
    } else {
	while (i2 < l2) {
	    if (s1[i1] == s2[i2]) {
		i1++,i2++;
		mm++;
	    } else {
		if (i2+2 < l2 && s1[i1+1]==s2[i2+1] && s1[i1+2]==s2[i2+2]) {
		    i1 += 3;
		    i2 += 3;
		    if (max_mm < mm)
			max_mm = mm;
		    mm = 0;
		} else {
		    break;
		}
	    }
	}
    }
    l += i1 - p1;

    if (max_mm < mm)
	max_mm = mm;

    *exact = max_mm;

    return l;
}
#endif

/*
 * NOTE - test sorting step in here.
 * By sorting all hits by location we can spot blocks neighbouring hash
 * key hits and find our blocks fast.
 *
 * Eg:
 * 1st word matches positions a,b,c
 * 2nd word matches positions d,e
 * 3rd word matches positions f,g,h
 *
 * Store a-0, b-0, c-0, d-1, e-1, f-2, g-2, h-2.
 * Sort and look for runs. Runs => consecutive query words matched
 * consecutive subject words and we have a longer match.
 */
int compare_b(Hash *h,
	      ALIGN_PARAMS *params, OVERLAP *overlap) {
    
    int	       	ncw, nrw, word, pw1, pw2, i, j, match_size;
    int	       	diag_pos, size_hist;
    int job_in;
    int pw_inc;

    if(h->seq1_len < h->min_match) return 0; 
    if(h->seq2_len < h->min_match) return 0; 

    size_hist = h->seq1_len + h->seq2_len - 1;
    for(i = 0; i < size_hist; i++) h->diag[i] = -(h->word_length);
    
    nrw = h->seq2_len - h->word_length + 1;
    
    /* 	loop for all (nrw) complete words in values2 ... */
    h->matches = -1;

    /*
     * TODO: Try pw2+=min_match_len-word_len as we know it'll have to have
     * at least one word hitting then.
     *
     * If we do this, we also need to adjust match_len to search backwards
     * as well as forwards.
     */

    pw_inc = h->min_match - h->word_length + 1;
    for (pw2=0; pw2<nrw; pw2+=pw_inc) {
 	if (-1 == (word = h->values2[pw2]))
	    continue;

	/* ... and look them up on values1, extrapolating to find length */
	
	if (0 == (ncw = h->counts[word]))
	    continue;

	for (j=0,pw1=h->last_word[word]; j<ncw; j++, pw1=h->values1[pw1]) {
	    int bck;

	    diag_pos = h->seq1_len - pw1 + pw2 - 1;

	    if (h->diag[diag_pos] >= pw2)
		continue;

	    /*
	    if ((match_size = fast_match_len(h->word_length,
					     h->seq1, pw1,
					     h->seq1_len,
					     h->seq2, pw2,
					     h->seq2_len)) >= h->min_match) {
	    */
	    if ((match_size = match_fwd_back(h->word_length,
					     h->seq1, pw1,
					     h->seq1_len,
					     h->seq2, pw2,
					     h->seq2_len,
					     &bck)) >= h->min_match) {
		h->matches++;
		if(h->max_matches == h->matches) {
		    h->max_matches *= 2;
		    h->block_match = (Block_Match *)
			xrealloc(h->block_match,
				 sizeof(Block_Match) *
				 h->max_matches);
		    if (NULL == h->block_match)
			return -5;
		}

		h->block_match[h->matches].pos_seq1 = pw1-bck;
		h->block_match[h->matches].pos_seq2 = pw2-bck;
		h->block_match[h->matches].length = match_size;
		h->block_match[h->matches].diag = diag_pos;
	    }
	    h->diag[diag_pos] = pw2-bck + match_size;
	}
    }
    h->matches += 1;
    if ( h->matches < 1 ) return 0;

    job_in = params->job;
    /* Force old-style return with end-gaps in alignment */
    params->job = RETURN_SEQ | RETURN_EDIT_BUFFERS | RETURN_END_GAPS;
    pw2 = align_blocks(h, params, overlap);
    params->job = job_in;

    return pw2;
}


/*
 * As per compare_b() but we're aligning seq2 against many seq1 choices
 * contatenated together. We need to collate hits by the contig they're in
 * so we can then do alignments on each in turn.
 */
int compare_b_bulk(Hash *h,
		   ALIGN_PARAMS *params, OVERLAP *overlap, int cnum2,
		   Contig_parms *contig_list1, int number_of_contigs1,
		   int ignore_after1,
		   void (*add_func)(OVERLAP *overlap,
				    int cnum1,
				    int cnum2,
				    void *clientdata),
		   void *add_data) {
    
    int ncw, nrw, word, pw1, pw2, pw2_last, j, match_size;
    int diag_pos;
    int job_in;
    int pw_inc;
    int size_div = 1;
    char *diag_init;

    if(h->seq1_len < h->min_match) return 0; 
    if(h->seq2_len < h->min_match) return 0; 

    size_div = (int)(1 + (h->seq1_len + h->seq2_len - 1)/DSZ);
    diag_init = calloc(size_div, sizeof(char));

    nrw = h->seq2_len - h->word_length + 1;
    
    /* 	loop for all (nrw) complete words in values2 ... */
    h->matches = -1;

    pw_inc = h->min_match - h->word_length + 1;
    pw2_last = 0;
    for (pw2=0; pw2<nrw; pw2+=pw_inc) {
 	if (-1 == (word = h->values2[pw2])) {
	    /*
	     * Illegal word means we may have missed a real match by
	     * jumping forward in steps of pw_inc. In this case back up one
	     * and try again, until we hit the last valid word place we
	     * checked.
	     */
	    if (pw2 > pw2_last) {
		pw2++;
		pw2 -= pw_inc;
		continue;
	    }
	    continue;
	}

	pw2_last = pw2;

	/* ... and look them up on values1, extrapolating to find length */
	
	if (0 == (ncw = h->counts[word]))
	    continue;

	if (h->filter_words && ncw > h->filter_words)
	    continue;

	for (j=0,pw1=h->last_word[word]; j<ncw; j++, pw1=h->values1[pw1]) {
	    int bck;

	    /* Remove matches to earlier and/or self contigs
	    if (pw1 <= contig_list1[cnum2].contig_end_offset)
	    continue; */
	    if (pw1 > ignore_after1) continue;

	    diag_pos = h->seq1_len - pw1 + pw2 - 1;

	    if (!diag_init[diag_pos / DSZ]) {
		int p, q;
		diag_init[diag_pos / DSZ] = 1;
		q = (int)(diag_pos / DSZ) * DSZ;
		for (p = 0; p < DSZ; p++, q++) {
		    h->diag[q] = -h->word_length;
		}
	    }
	    if (h->diag[diag_pos] >= pw2) continue;

	    /*
	    match_size = fast_match_len(h->word_length,
					h->seq1, pw1,
					h->seq1_len,
					h->seq2, pw2,
					h->seq2_len);
	    */
	    match_size = match_fwd_back(h->word_length,
					h->seq1, pw1,
					h->seq1_len,
					h->seq2, pw2,
					h->seq2_len,
					&bck);
	    /*
	    match_size = match_fwd_back_mm(h->word_length,
					h->seq1, pw1,
					h->seq1_len,
					h->seq2, pw2,
					h->seq2_len,
					&bck, &exact);
	    if (exact >= h->min_match) {
	    */

	    if (match_size >= h->min_match) {
	    //if (match_size >= h->min_match || ncw < 3) {
		h->matches++;
		if(h->max_matches == h->matches) {
		    h->max_matches *= 2;
		    h->block_match = (Block_Match *)
			xrealloc(h->block_match,
				 sizeof(Block_Match) *
				 h->max_matches);
		    if (NULL == h->block_match)
			return -5;
		}

		h->block_match[h->matches].pos_seq1 = pw1-bck;
		h->block_match[h->matches].pos_seq2 = pw2-bck;
		h->block_match[h->matches].length = match_size;
		h->block_match[h->matches].diag = diag_pos;
	    }

	    h->diag[diag_pos] = pw2-bck + match_size;
	}
    }

    free(diag_init);

    h->matches += 1;
    if ( h->matches < 1 ) return 0;

    /* printf("Matches = %d\n", h->matches); */

    //UpdateTextOutput();

    job_in = params->job;
    params->job = 3; /* force return of edit buffers */
    pw2 = align_blocks_bulk ( h, params, overlap, cnum2,
			      contig_list1, number_of_contigs1,
			      add_func, add_data);
    params->job = job_in;

    return pw2;
}


int
gap_realloc_matches (int **seq1_match, 
		     int **seq2_match,
		     int **len_match,
		     int *max_matches) 
{
    int increment = 1000;

    *max_matches += increment;

    if (NULL == (*seq1_match = (int *)xrealloc(*seq1_match, 
					      *max_matches * sizeof(int)))) {
	return -1;
    }
    if (NULL == (*seq2_match = (int *)xrealloc(*seq2_match, 
					      *max_matches * sizeof(int)))) {
	return -1;
    }

    /* len_match is not always needed */
    if (len_match != NULL) {
	if (NULL == (*len_match = (int *)xrealloc(*len_match, 
						 *max_matches * sizeof(int)))) {
	    return -1;
	}
    }
    return 0;
}


int reps(Hash *h, 
	 int **seq1_match_pos, 
	 int **seq2_match_pos,
	 int **match_length, 
	 int offset,
	 char sense) {
    
    int ncw, nrw, word, pw1, pw2, pw2_last, i, j, match_size;
    int diag_pos, size_hist;
    int pw_inc, bck;
    /* clock_t t1, t2; */

    if(h->seq1_len < h->min_match) return -4; 
    if(h->seq2_len < h->min_match) return -4; 

    size_hist = h->seq1_len + h->seq2_len - 1;
    for(i = 0; i < size_hist; i++) h->diag[i] = -(h->word_length);

    /* if forward repeats make sure we do not bother with the main diagonal */

    if ( 'f' == sense ) h->diag[h->seq1_len-1] = h->seq1_len;
    
    nrw = h->seq2_len - h->word_length + 1;
    
    /* t1 = clock(); */

    /* 	loop for all (nrw) complete words in values2 */
    h->matches = -1;
    pw_inc = h->min_match - h->word_length + 1;
    pw2_last = 0;
    for (pw2=0;pw2<nrw;pw2+=pw_inc) {
 	if (-1 == (word = h->values2[pw2])) {
	    if (pw2 > pw2_last) {
		pw2++;
		pw2 -= pw_inc;
		continue;
	    }
	    continue;
	}

	pw2_last = pw2;

	if (0 == (ncw = h->counts[word]))
	    continue;

	for (j=0,pw1=h->last_word[word];j<ncw;j++) {
	    diag_pos = h->seq1_len - pw1 + pw2 - 1;

	    /* Remove self match */
	    if (sense == 'f' && pw1 < pw2) {
		pw1 = h->values1[pw1];
		continue;
	    }

	    if ( h->diag[diag_pos] < pw2 ) {
		if ((match_size = match_fwd_back(h->word_length,
						 h->seq1, pw1, h->seq1_len,
						 h->seq2, pw2, h->seq2_len,
						 &bck)) >= h->min_match) {
		    h->matches++;
		    if(h->max_matches == h->matches+offset) {
			
			if (-1 == (gap_realloc_matches(seq1_match_pos,
						       seq2_match_pos,
						       match_length, 
						       &h->max_matches))) {
			    return -1;
			    /* return -5; */
				    
			}
		    }
		    (*seq1_match_pos)[h->matches+offset] = pw1+1-bck;
		    (*seq2_match_pos)[h->matches+offset] = pw2+1-bck;
		    (*match_length)[h->matches+offset] = match_size;
		}
		h->diag[diag_pos] = pw2-bck + match_size;
	    }
	    pw1 = h->values1[pw1];
	}
    }
    h->matches += 1;

    /*
    t2 = clock();

    printf("Time = %f, nmatches = %d\n",
	   (double)(t2-t1) / CLOCKS_PER_SEC,
	   h->matches);
    */

    if ( h->matches ) {
	if ( sense == 'r' ) {
	    (void) make_reverse ( seq2_match_pos, match_length,
				 h->matches, h->seq2_len, offset); 	
	}
	(void) remdup ( seq1_match_pos, seq2_match_pos, match_length, offset, 
			&h->matches );
    }

    return h->matches;
}

int reps_nocount(Hash *h, 
	 int **seq1_match_pos, 
	 int **seq2_match_pos,
	 int **match_length, 
	 int offset,
	 char sense) {
    
    int nrw, word, pw1, pw2, pw2_last, i, match_size;
    int diag_pos, size_hist;
    int pw_inc, bck;
    /* clock_t t1, t2; */

    if(h->seq1_len < h->min_match) return -4; 
    if(h->seq2_len < h->min_match) return -4; 

    size_hist = h->seq1_len + h->seq2_len - 1;
    for(i = 0; i < size_hist; i++) h->diag[i] = -(h->word_length);

    /* if forward repeats make sure we do not bother with the main diagonal */

    if ( 'f' == sense ) h->diag[h->seq1_len-1] = h->seq1_len;
    
    nrw = h->seq2_len - h->word_length + 1;
    
    /* t1 = clock(); */

    /* 	loop for all (nrw) complete words in values2 */
    h->matches = -1;
    pw_inc = h->min_match - h->word_length + 1;
    pw2_last = 0;
    for (pw2=0;pw2<nrw;pw2+=pw_inc) {
 	if (-1 == (word = h->values2[pw2])) {
	    if (pw2 > pw2_last) {
		pw2++;
		pw2 -= pw_inc;
		continue;
	    }
	    continue;
	}

	pw2_last = pw2;

	for (pw1=h->last_word[word]; pw1 != -1; pw1 = h->values1[pw1]) {
	    diag_pos = h->seq1_len - pw1 + pw2 - 1;

	    /* Remove self match */
	    if (sense == 'f' && pw1 < pw2) {
		continue;
	    }

	    if ( h->diag[diag_pos] < pw2 ) {
		if ((match_size = match_fwd_back(h->word_length,
						 h->seq1, pw1, h->seq1_len,
						 h->seq2, pw2, h->seq2_len,
						 &bck)) >= h->min_match) {
		    h->matches++;
		    if(h->max_matches == h->matches+offset) {
			
			if (-1 == (gap_realloc_matches(seq1_match_pos,
						       seq2_match_pos,
						       match_length, 
						       &h->max_matches))) {
			    return -1;
			    /* return -5; */
				    
			}
		    }
		    (*seq1_match_pos)[h->matches+offset] = pw1+1-bck;
		    (*seq2_match_pos)[h->matches+offset] = pw2+1-bck;
		    (*match_length)[h->matches+offset] = match_size;
		}
		h->diag[diag_pos] = pw2-bck + match_size;
	    }
	}
    }
    h->matches += 1;

    /*
    t2 = clock();

    printf("Time = %f, nmatches = %d\n",
	   (double)(t2-t1) / CLOCKS_PER_SEC,
	   h->matches);
    */

    if ( h->matches ) {
	if ( sense == 'r' ) {
	    (void) make_reverse ( seq2_match_pos, match_length,
				 h->matches, h->seq2_len, offset); 	
	}
	(void) remdup ( seq1_match_pos, seq2_match_pos, match_length, offset, 
			&h->matches );
    }

    return h->matches;
}
