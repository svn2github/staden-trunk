/* 
   need to introduce seq1_lreg, seq1_rreg, ditto seq2 at the right
   level, then from there down seq1_len is seq1_rreg - seq1_lreg + 1
   and we send &seq1[seq1_lreg-1] 
   note positions are in "user units" ie sequence numbers from 1 
   please check the results of each with a reduced region - I have
   only checked seq1_lreg=1, seq1_rreg=seq1_len.
*/
/* 

new functions: 

1. 

int inexact_match ( char *seq, int seq_len, char *string, int string_len,
		   int min_match, int *match, int *score, int max_matches) {

     Find all positions in seq at which string has >= min_match matching
     characters.
     Requires SEQ_MISMATCH and ive included SEQ_MATCH in inexact_match for fun.

#define SEQ_MISMATCH(a,b) ( ((match_lookup[a]) != (match_lookup[b])) ? 1 : 0 )
#define SEQ_MATCH(a,b) ( ((match_lookup[a]) == (match_lookup[b])) ? 1 : 0 )

2.
   
p_compare_seqs (

     Does a hash comparison without storing the results. Assumes the word_length
     is the hash length and the minimum match. At present prints the matching
     positions but should be modified to plot them.
     Requires a change to hash_compare (maybe a job number? to tell it which of
     compare_seqs and p_compare_seqs to use)

*/

#include <math.h>

#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "sequence_formats.h"
#include "readpam.h"
#include "misc.h" /* MIN MAX */
#include "sip_results.h"  /* Point */
#include "dna_utils.h"

/* 
   General hashing routines for dna and proteins

   For dna we have the fast 8 base hashing method but this will not work for words 
   less than 8 or for proteins. So lets use the old ones from sip which are more
   general. However we can make use of many of the other routines for storing the hash
   values etc. We also have to generalise the comparisons so they know about protein 
   as well as dna. It might be convenient to allow users to use a membership of a set
   for comparing characters - ie in the table that defines identity make hydrophobics
   match etc. this could also be derived from the score matrix table: any pairs of
   characters that have some minimum score are counted as identical. In this way hashing
   is less of binary type comparison.

hash_consts [n_hash_consts] a set of constants for hashing
hash_length the hash length
char_set_size the character set size - dna or protein, 5 or 24
hash_values the hash values for a sequence
we have a table of integer values for dna (dna_lookup) and one for 
protein (protein_lookup) and may want to use either depending on the 
type of sequences being compared. Can we use a single pointer (char_lookup) for
the current type, and make it global?
Suppose we also want to compare protein against dna? I guess we make
a temporary translation of the dna and use the protein table.
*/

/* int seq_type, DNA, char_set_size, *char_lookup, *match_lookup; */

#define MAX_HASH_CONST 200

static int hash_const[MAX_HASH_CONST], word_length;

/************************************************************/

   /*
    * Set up the constants used by the general hashing routine
    * The values depend on the character set size and the hash
    * word length.
    *
    */

void set_hash_consts (void) {
    int i,j,k,m,n;

    for (i=0,k=char_set_size-1,n=0,hash_const[0]=0;i<word_length;i++) {
	m = pow( k,i );
	hash_const[0] -= hash_const[n];
	for (j=1;j<=k;j++) {
	    n++;
	    hash_const[n] = j * m;
	}
    }
}

/*
 *
 * return the hash value for the string in seq
 *
 */

int hash_value (char seq[]) {
    int i,j,k,l,n;
    for ( i=0,n=hash_const[0],j=0,k=char_set_size-1;i<word_length;i++) {
	l = char_lookup[seq[i]] + 1;
	if (char_set_size == l) return -1;
	n += hash_const[j+l];
	j += k;
    }
    return n;
}

int hash_seq ( char *seq, int *hash_values, int seq_len ) {

    int i,j,k,l,m,n,good,bad;
    memset ( hash_values, 0, sizeof(*hash_values) * seq_len);
    good = 1;
    for ( m = 0; m < seq_len-word_length + 1; m++) {
	for ( i=0,n=hash_const[0],j=0,k=char_set_size-1,bad=0;i<word_length;i++) {
	    l = char_lookup[seq[m+i]] + 1;
	    if ( l == char_set_size ) bad = 1;
	    n += hash_const[j+l];
	    j += k;
	}
	if ( bad ) {
	    hash_values[m] = -1;
	}
	else {
	    hash_values[m] = n - 1;
	    good = 0;
	}
    }
    return good;
}

/* As a special case include hashing on word length 8 for DNA */

static int dna_hash8_lookup[256];

void set_hash8_lookup(void) {

/*	hashing values */
/* 	set up table of values for permitted base characters */

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
/*    dna_hash8_lookup['*'] = 0; */

}

int hash_word8 ( char *seq, int *start_base, int seq_len,
	      unsigned short *uword) {

    /* 	given a sequence seq, return the hash value for the first word 
     *  after start_base that does not contain an unknown char. Tell 
     *  the caller where this is. If we reach the end of the seq set
     *  start_base and return -1.
     */


    register int i, word_len=8;
    register int end_base,base_index,lstart_base;
    register int unsigned short luword;

    lstart_base = *start_base;
    end_base = lstart_base + word_len;
    if ( seq_len < end_base ) return -1;

    for (i=lstart_base,luword=0,end_base=lstart_base+word_len;i<end_base;i++) {

	base_index = dna_hash8_lookup[(unsigned)seq[i]];
	if ( 4 == base_index ) {

	    /*	weve hit an unknown char, so lets start again */

	    lstart_base = i + 1;
	    end_base = lstart_base + word_len;
	    if ( seq_len < end_base ) {
		*start_base = lstart_base;
		return -1;
	    }
	    luword = 0;
	    i = lstart_base - 1;
	}
	else {
	    luword = ( luword <<2 ) | base_index;
	}
    }
    *start_base = lstart_base;
    *uword = luword;
    return 0;
}


int hash_seq8 ( char *seq, int *hash_values, int seq_len) {

/* given a sequence seq, return an array of hash values
   If we cannot find at least one word to hash on we return -1
   otherwise we return 0.
*/

    register int i,j,k,word_len=8;
    int start_base,prev_start_base,end_base,base_index;
    unsigned short uword;

    if ( seq_len < word_len ) return -1;

    /*	Get the hash value for the first word that contains no unknowns */	
    start_base = 0;
    if (hash_word8 ( seq, &start_base, seq_len, &uword)) return -1;

    for (i=0;i<start_base;i++) hash_values[i] = -1;

    /*	Now do the rest of the sequence */

    hash_values[start_base] = uword;
    k = seq_len - word_len + 1;

    for (i=start_base+1,j=start_base+word_len; i<k; i++,j++) {

	base_index = dna_hash8_lookup[(unsigned)seq[j]];
	if ( 4 == base_index ) {

	    /*	weve hit an unknown char, so lets start again */

	    prev_start_base = i;
	    start_base = j + 1;
	    if (hash_word8 ( seq, &start_base, seq_len, &uword)) {
		for (i=prev_start_base;i<start_base;i++) hash_values[i] = -1;
		return 0;
	    }

	    for (i=prev_start_base;i<start_base;i++) hash_values[i] = -1;
	    hash_values[start_base] = uword;
	    end_base = start_base + word_len;
	    i = start_base;
	    j = i + word_len - 1;

	}
	else {
	    uword = ( uword <<2 ) | base_index;
	    hash_values[i] = uword;
	}
    }
    return 0;
}



/************************************************************/

void store_hash ( 
	         int *hash_values, 	/* the hash values for each position in a seq */
		 int seq_len, 		/* size of the seq and hash array */
		 int *last_word, 	/* last occurrence of this hash value */
		 int *word_count, 	/* frequency of each hash value or word */
		 int size_wc ) {	/* number of elements in word_count and first_word */


/* 	store the hash values in hash_values: put number of occurrences of
	each hash value in word_count; put the array position of the last 
	occurrence of each hash value in last_word, and previous
	occurrences in hash_values[last_word]. 
	Note that words containing unknown characters (like '-') are given
	hash value -1. So we skip them here, and they are ignored.
*/

    int nw;
    register int i,j,n;

    for ( i=0;i<size_wc;i++ ) {
	word_count[i] = 0;
	last_word[i] = 0;
    }

    /* loop for all entries in hash_values	*/

    j = seq_len - word_length + 1;

    for ( i = 0; i < j; i++ ) {

	n = hash_values[i];

	/* is it a good value ? */

	if ( -1 != n ) {

	    nw = word_count[n];

	    /* already an entry for this word ? */

	    if ( 0 == nw ) {

		/* no, so put in last_word */

		last_word[n] = i;
		word_count[n] += 1;
	    }

	    /*	yes, so put previous last occurrence in hash_values*/

	    else {

		word_count[n] += 1;
		hash_values[i] = last_word[n];
		last_word[n] = i;
	    }
	}
    }
}

int
sip_realloc_matches (int **seq1_match, 
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

void sip_remdup ( int **seq1_match, 
		 int **seq2_match, 
		 int **len_match, 
		 int *n_matches ) 
{
    /*	routine to remove duplicates from a list of repeats * * It
	also removes the self match.
	Input: a list of *n_match match positions and match
	lengths.  * Output: a list in which all duplicates and the selfmatch
	are removed.  * *n_match is set to the new number of matches, or -1
	for error.  */

    
    register int i;
    int *index_ptr,k,keep;
    
/*	if no matches why bother */

    if ( *n_matches < 1 ) return;

    if ( ! ( index_ptr = (int *) xmalloc ( sizeof(int)*(*n_matches) ))) {
 	*n_matches = -1;	
	return;
    }
    
    /* 	Go thru list and note the index of duplicates */
    /* Every match pair, say 15,5, will have a duplicate 5,15
     * so remove all those for which the first position is < the second
     */
    
    for ( i = 0, k = 0; i < *n_matches; i++) {

	/* FIXME can we use > only here or is = required ??? */
	if ((*seq1_match)[i] >= (*seq2_match)[i]) {
	    index_ptr[k] = i; 
	    k += 1;
	}
    }
    
    for ( i = 0; i < k; i++) {
 	keep = index_ptr[i];
	(*seq1_match)[i] = (*seq1_match)[keep];
	(*seq2_match)[i] = (*seq2_match)[keep];
	if (len_match)
	    (*len_match)[i] = (*len_match)[keep];
    } 
    *n_matches = k;

    if ( index_ptr ) free ( index_ptr );

}


/************************************************************/

int compare_seqs (
		   int *hash_values1, 	/* the hash values for each position in seq1 */
		   int *last_word, 	/* last occurrence of this hash value in seq1 */
		   int *word_count, 	/* frequency of each hash value or word */
		   int *hash_values2,	/* the hash values for seq2 */
		   int min_match,	/* the minimum match length */
		   int **seq1_match,	/* positions of matches in seq1 */
		   int **seq2_match,	/* positions of matches in seq2 */
		   int **len_match,	/* length of matches */
		   int max_matches,	/* maximum number of matches */
		   char *seq1,
		   char *seq2,		
		   int seq1_len, 	/* size of seq1 */
		   int seq2_len,	/* the length of seq2 */
		   int *diag,     	/* for knowing how far weve been on a diag */
		   int seq1_lreg,
		   int seq2_lreg,
		   int same_seq

		   ) {

/*	we have a sequence seq encoded by hash tables:
	word_count[i] contains the number of occurrences of hash value i
	last_word[i] contains the seq position of the last occurrence of hash value i
	hash_values1[i] contains the hash values of previous occurrence of hash value i
	
	we have seq2 encoded in an array of hash values hash_values2[i], each i
	corresponding to the string at seq2[i]

	algorithm: for each hash value in hash_values2[i] test if it exists in word_count
	if so process all occurrences using last-word and hash_values1.
	we follow matches found this way from start to end to check they really match
	and to find their length ( the hashing routine does not distinguish "*" and "a"
	characters in the consensus, so we must do it here ).
	For dissimilar sequences we gain by checking each match start against all those
	weve already saved to see if is a continuation.
*/


    int word, match_length, n_matches;
    register int i,j,diag_pos, pw1, pw2, nrw, ncw;


    n_matches = 0;
    j = seq1_len + seq2_len;
    for (i=0;i<j;i++) diag[i] = -word_length;

    /* if looking for internal repeats skip main diagonal CHECK TIME SAVED BY THIS*/

    if ( ( same_seq ) && ( seq1_lreg == seq2_lreg ) ) diag[seq1_len - 1] = seq1_len;

    nrw = seq2_len - word_length + 1;

/* 	loop for all (nrw) complete words in hash_values2 */

    for (pw2 = 0; pw2 < nrw; pw2++) {

	word = hash_values2[pw2];

/*	legal character ? */

	if ( -1 != word ) {

/*		in consensus ?	*/

	    if ( 0 != (ncw = word_count[word]) ) {

/*		yes, so process all (ncw) of them	*/

/*		first occurrence to process is in last_word 	*/

		for (j = 0, pw1 = last_word[word]; j < ncw; j++) {

		    diag_pos = seq1_len - pw1 + pw2 - 1;
		    if ( diag[diag_pos] < pw2 ) {

/*
 *		this is not a continuation, so check its whole length 
*/

			if ( (match_length = match_len (
					     seq1, pw1, seq1_len, 
					     seq2, pw2, seq2_len )) >= min_match ) {

			    diag[diag_pos] = pw2 + match_length;
			    if ( n_matches < max_matches ) {

				(*seq1_match)[n_matches] = pw1;
				(*seq2_match)[n_matches] = pw2;
				(*len_match)[n_matches] = match_length;
				n_matches += 1;

			    }
			    else {

				/* Error: too many matches */

/*				verror(ERR_WARN, "compare_seq", "too many matches");*/
				if (-1 == (sip_realloc_matches(seq1_match, seq2_match,
							       len_match, &max_matches))) {
				    return -1;
				}
				
				(*seq1_match)[n_matches] = pw1;
				(*seq2_match)[n_matches] = pw2;
				(*len_match)[n_matches] = match_length;
				n_matches += 1;
			    }
			}
		    }

/*		get next occurrence of this word */
		    pw1 = hash_values1[pw1];

		}
	    }
	}
    }

/*	Now put the element indexes to values 1 to end */


    for ( i = 0; i<n_matches; i++) {

	(*seq1_match)[i] += seq1_lreg;
	(*seq2_match)[i] += seq2_lreg;
    }

    return n_matches;
}

/************************************************************/

int p_compare_seqs (
		     int *hash_values1, /* the hash values for each position in seq1 */
		     int *last_word, 	/* last occurrence of this hash value in seq1 */
		     int *word_count, 	/* frequency of each hash value or word */
		     int *hash_values2,	/* the hash values for seq2 */
		     char *seq1,
		     char *seq2,		
		     int seq1_len, 	/* size of seq1 */
		     int seq2_len,	/* the length of seq2 */
		     void (*pr_func)(void *data, int pw, int y),/*void (*pr_func)()*/	
		     void *data
		     ) {

/*	we have a sequence seq encoded by hash tables:
	word_count[i] contains the number of occurrences of hash value i
last_word[i] contains the seq position of the last occurrence of hash value i
	hash_values1[i] contains the hash values of previous occurrence of hash value i
	
	we have seq2 encoded in an array of hash values hash_values2[i], each i
	corresponding to the string at seq2[i]

	algorithm: for each hash value in hash_values2[i] test if it exists in word_count
	if so process all occurrences using last-word and hash_values1.


*/

    int nrw, ncw, word, pw;
    register int i,j;
    int y;
    int n_matches = 0;
    double wx0, wy0, wx1, wy1;

    RasterGetWorldScroll(data, &wx0, &wy0, &wx1, &wy1);

    nrw = seq2_len - word_length + 1;

/* 	loop for all (nrw) complete words in hash_values2 */

    for (i=0;i<nrw;i++) {

	word = hash_values2[i];

/*	legal character ? */

	if ( -1 != word ) {

/*		in consensus ?	*/

	    if ( 0 != (ncw = word_count[word]) ) {

/*		yes, so process all (ncw) of them	*/

/*		first occurrence to process is in last_word 	*/

		pw = last_word[word];

		y = (int)rasterY(data, (i+1));
		for (j=0;j<ncw;j++) {
		    /* print data directly to plotting widget */
		    pr_func(data, pw+1, y);
		    /* get next occurrence of this word */
		    pw = hash_values1[pw];
		    n_matches++;
		}
	    }
	}
    }
    return n_matches;
}

void histel_to_xy1 ( int seq1_len, int histel, int *x, int *y )
{

    if ( histel < seq1_len ) {
	*y = seq1_len - histel -1;
	*x = 0;
    }
    else {
	*x = histel - seq1_len + 1;
	*y = 0;
    }
}

void histel_to_xy ( int seq1_len, int histel, int *x, int *y )
{

    if ( histel < seq1_len ) {
	*y = seq1_len - histel -1;
	*x = 0;
    }
    else {
	*x = histel - seq1_len + 1;
	*y = 0;
    }
}

/************************************************************/
int q_compare_seqs (
		     int *hash_values1, /* the hash values for each position in seq1 */
		     int *last_word, 	/* last occurrence of this hash value in seq1 */
		     int *word_count, 	/* frequency of each hash value or word */
		     int *hash_values2,	/* the hash values for seq2 */
		     int *diag,		/* how far weve reached on each diagonal */
		     char *seq1,
		     char *seq2,		
		     int seq1_len, 	/* size of seq1 */
		     int seq2_len,	/* the length of seq2 */
		     int span_length,
		     int min_span_score,
		     double min_sd,     /* all diagonals with sd above this are plotted */
		     int max_matches,
		     int save_results,
		     int **seq1_match,	/* positions of matches in seq1 */
		     int **seq2_match,	/* positions of matches in seq2 */
	             void (*pr_func)(void *data, int pw, int y), /*void (*pr_func)()*/	
		     void *data,
		     int seq1_lreg,
		     int seq2_lreg
		   ) {

/* KATHRYN: need access to score_matrix, span_length, min_span_score and min_sd */

/*	we have a sequence seq encoded by hash tables:
	word_count[i] contains the number of occurrences of hash value i
	last_word[i] contains the seq position of the last occurrence of hash value i
	hash_values1[i] contains the hash values of previous occurrence of hash value i
	
	we have seq2 encoded in an array of hash values hash_values2[i], each i
	corresponding to the string at seq2[i]

	algorithm: for each hash value in hash_values2[i] test if it exists in word_count
	if so process all occurrences using last-word and hash_values1.

	This is quick scan so processing means accumulate the scores for each diagonal
	in a histogram, correct the histogram for the length of its diagonal, find the 
	mean and sd of the histogram, find the diagonals that exceed min_sd, scan along
	these with a window of length span_length, and plot all points above min_span_score.
	For some reason we are careful with memory and use a union for the histogram.
	seq1 on x axis. seq2 on y.

*/

    typedef union double_int
	{
	    int i;
	    double d;
	}di;
    int nrw, ncw, word, size_hist, x, y;
    int frontx, fronty, backx, backy, score, halfsp, match_length;
    register int i,j,pw1,pw2,diag_pos;
    double rj, rm, rmsq, sd, rm2, t, min_score, points, max_diagonal;
    union double_int double_int1, *hist;
    double wx0, wy0, wx1, wy1;
    int n_matches = 0;


    if (data) {
	RasterGetWorldScroll(data, &wx0, &wy0, &wx1, &wy1);
    }

    nrw = seq2_len - word_length + 1;
    size_hist = seq1_len + seq2_len;
    if ( ! (hist = (di *) xmalloc ( sizeof(double_int1)*(size_hist) ))) {
	verror(ERR_FATAL, "quick scan", "out of memory");
	goto bail_out;
    }

    j = seq1_len + seq2_len;
    for(i=0;i<j;i++) diag[i] = -word_length;

    for(i=0;i<size_hist;i++)(hist[i].i = 0);
/* 	loop for all (nrw) complete words in hash_values2 */
    for (pw2=0;pw2<nrw;pw2++) {
	word = hash_values2[pw2];
	if ( -1 != word ) {
	    if ( 0 != (ncw = word_count[word]) ) {

		for (j=0,pw1=last_word[word];j<ncw;j++) {
		    diag_pos = seq1_len - pw1 + pw2 -1;
		    if (diag[diag_pos]<pw2) {
			match_length = match_len(seq1,pw1,seq1_len,
						 seq2,pw2,seq2_len);
			diag[diag_pos] = pw2 + match_length;
			hist[diag_pos].i += match_length;
		    }
		    pw1 = hash_values1[pw1];
		}
	    }
	}
    }

    /* now do the correction for the length of each diagonal */

    max_diagonal = MIN(seq1_len,seq2_len);
    for(i=0,points=1.0;i<seq1_len;i++) {
	hist[i].d = hist[i].i / points;
	points += 1.0;
	points = MIN (max_diagonal, (points));
    }
    for(i=size_hist-1,points=1.0;i>seq1_len-1;i--) {
	hist[i].d = hist[i].i / points;
	points += 1.0;
	points = MIN (max_diagonal, (points));
    }
    rmsq = 0.0;
    rm = 0.0;
    for(i=0;i<size_hist;i++) {
	rj = hist[i].d;
	rm += rj;
	rmsq += rj*rj;
    }

    rm = rm / size_hist;
    rmsq = rmsq / size_hist;
    rm2 = rm * rm;
    sd = 0.0;
    t = rmsq - rm2;
    if (t > 0.0) sd = sqrt(t);
    min_score = rm + min_sd * sd;
    for(i=0;i<size_hist;i++) {
	hist[i].i = ( hist[i].d < min_score ) ? 0 : 1;
    }
    halfsp = span_length/2;

    for(i=span_length-1;i<size_hist-span_length+1;i++) {
	if ( hist[i].i ) {
	    (void) histel_to_xy ( seq1_len, i, &x, &y );
	    /* get started */
	    score = 0;
	    for(j=0,frontx=x, fronty=y;j<span_length;j++,frontx++,fronty++) {
		score += score_matrix[char_lookup[seq2[frontx]]][char_lookup[seq1[fronty]]];
	    }
	    if (score >= min_span_score) {
		if (save_results) {
		    if ( n_matches < max_matches ) {
			/* FIXME - is this correct? */
			(*seq1_match)[n_matches] = fronty-halfsp;
			(*seq2_match)[n_matches] = frontx-halfsp;
			n_matches += 1;
		    } else {
			if (-1 == (sip_realloc_matches(seq1_match, seq2_match,
						       NULL, &max_matches))) {
			    goto bail_out;
			}
			(*seq1_match)[n_matches] = fronty-halfsp;
			(*seq2_match)[n_matches] = frontx-halfsp;
			n_matches += 1;
		    }
		} else {
		    pr_func(data, frontx-halfsp, 
			    (int)rasterY(data, fronty-halfsp));
		}
	    }
	    /* do rest */

	    backx = x;
	    backy = y;
	    for(;frontx<seq2_len&&fronty<seq1_len;
		frontx++,fronty++,backx++,backy++) {
		score += score_matrix[char_lookup[seq2[frontx]]][char_lookup[seq1[fronty]]]
        	      -  score_matrix[char_lookup[seq2[backx]]][char_lookup[seq1[backy]]];
		if (score >= min_span_score) {
		    if (save_results) {
			if ( n_matches < max_matches ) {
			    (*seq1_match)[n_matches] = fronty-halfsp;
			    (*seq2_match)[n_matches] = frontx-halfsp;
			    n_matches += 1;
			} else {
			    if (-1 == sip_realloc_matches(seq1_match, 
							  seq2_match,
							  NULL, &max_matches)){
				goto bail_out;

			    }
			    (*seq1_match)[n_matches] = fronty-halfsp;
			    (*seq2_match)[n_matches] = frontx-halfsp;
			    n_matches += 1;
			}
		    } else {
			pr_func(data, frontx-halfsp+1, 
				(int)rasterY(data, fronty-halfsp+1));
		    }
		}
	    }
	}
    }
    free ( hist );

    for ( i = 0; i<n_matches; i++) {

	(*seq1_match)[i] += seq1_lreg;
	(*seq2_match)[i] += seq2_lreg;
    }
    return n_matches;

 bail_out:
    if (hist) free ( hist );
    return -1;
}


/************************************************************/

int hash_compare(
		 int min_match,	        /* the minimum match length */
		 int **seq1_match,	/* positions of matches in seq1 */
		 int **seq2_match,	/* positions of matches in seq2 */
		 int **len_match,	/* length of matches */
		 int max_matches,	/* maximum number of matches */
		 char *seq1,		/* seq1 */
		 char *seq2,		/* seq2 */
		 int seq1_lreg,         /* start of seq1 */
		 int seq1_rreg,         /* end of seq1 */
		 int seq2_lreg,         /* start of seq2 */
		 int seq2_rreg,         /* end of seq2 */
		 int same_seq,		/* if set then its an internal comparison */
		 void (*pr_func)(void *data, int pw, int y),     /* printing function */
		 void *data             /* printing data */
		   ) {

/*	Main interface to the sequence comparison method using hashing.

	Jobs are: allocate, hash and compare 2 sequences, deallocate

*/

    static int *hash_values1,*last_word,*word_count,*hash_values2;
    int n_matches,size_hash;
    int *diag;
    int seq1_len, seq2_len;

    hash_values1 = hash_values2 = last_word = word_count = diag = NULL;

    seq1_len = seq1_rreg - seq1_lreg + 1;
    seq2_len = seq2_rreg - seq2_lreg + 1;

    /* sanity checks */

    if ( seq1_len < MAX(word_length,min_match) ) goto bail_out;
    if ( seq2_len < MAX(word_length,min_match) ) goto bail_out;

    set_hash_consts ();

    size_hash = pow(char_set_size-1,(MIN(word_length,min_match)));


    /* allocate space for seq1	*/

    if ( ! (hash_values1 = (int *) xmalloc ( sizeof(int)*(seq1_len) ))) {
	goto bail_out;
    }

    if ( ! (last_word = (int *) xmalloc ( sizeof(int)*size_hash ))) {
	goto bail_out;
    }
    if ( ! (word_count = (int *) xmalloc ( sizeof(int)*size_hash ))) {
	goto bail_out;
    }
    
    if ( ! (hash_values2 = (int *) xmalloc ( sizeof(int)*(seq2_len) ))) {
	goto bail_out;
    }
    if ( ! (diag = (int *) xmalloc ( sizeof(int)*(seq1_len+seq2_len) ))) {
	goto bail_out;
    }

    /* hash seq1	*/
    
    if ( hash_seq ( &seq1[seq1_lreg-1], hash_values1, seq1_len )  != 0 ) {
	goto bail_out;
    }
    (void) store_hash ( hash_values1, seq1_len, last_word, word_count,
		       size_hash);

    /* hash seq2	*/

    if ( hash_seq ( &seq2[seq2_lreg-1], hash_values2, seq2_len )  != 0 ) {
	goto bail_out;
    }

    /* do the comparison	*/
    
    n_matches = compare_seqs ( hash_values1, last_word, word_count, 
				  hash_values2, min_match, seq1_match, 
				  seq2_match, len_match, max_matches, &seq1[seq1_lreg-1], 
				  &seq2[seq2_lreg-1], seq1_len, seq2_len, diag, 
			          seq1_lreg, seq2_lreg, same_seq);

    /* if looking for internal matches remove duplicates */

    /* if ( same_seq && ( n_matches > 0 ) ) { */
    /* 
     * remove duplicates if same seq, n_matches is > 0 and the matches have
     * been saved
     */
    if ( same_seq && ( n_matches > 0 ) ) {
	(void) sip_remdup ( seq1_match, seq2_match, len_match, &n_matches );
    }

    if ( hash_values1 ) free ( hash_values1 );
    if ( hash_values2 ) free ( hash_values2 );
    if ( word_count )   free ( word_count );
    if ( last_word )    free ( last_word );
    if ( diag )         free ( diag );
    return n_matches;

 bail_out:
    if ( hash_values1 ) free ( hash_values1 );
    if ( hash_values2 ) free ( hash_values2 );
    if ( word_count )   free ( word_count );
    if ( last_word )    free ( last_word );
    if ( diag )         free ( diag );
    return -1;
}

/************************************************************/

int hash_compare8(
		 int min_match,	        /* the minimum match length */
		 int **seq1_match,	/* positions of matches in seq1 */
		 int **seq2_match,	/* positions of matches in seq2 */
		 int **len_match,	/* length of matches */
		 int max_matches,	/* maximum number of matches */
		 char *seq1,		/* seq1 */
		 char *seq2,		/* seq2 */
		 int seq1_lreg,         /* start of seq1 */
		 int seq1_rreg,         /* end of seq1 */
		 int seq2_lreg,         /* start of seq2 */
		 int seq2_rreg,         /* end of seq2 */
		 int same_seq,		/* if set then its an internal comparison */
		 void (*pr_func)(void *data, int pw, int y),     /* printing function */
		 void *data             /* printing data */
		   ) {

/*	Main interface to the sequence comparison method using hashing.

	Jobs are: allocate, hash and compare 2 sequences, deallocate

*/

    static int *hash_values1,*last_word,*word_count,*hash_values2;
    int n_matches,size_hash;
    int *diag;
    int seq1_len, seq2_len;

    hash_values1 = hash_values2 = last_word = word_count = diag = NULL;

    seq1_len = seq1_rreg - seq1_lreg + 1;
    seq2_len = seq2_rreg - seq2_lreg + 1;

    /* sanity checks */

    if ( seq1_len < min_match ) goto bail_out;
    if ( seq2_len < min_match ) goto bail_out;

    set_hash8_lookup ();

    size_hash = 65536;

    /* allocate space for seq1	*/

    if ( ! (hash_values1 = (int *) xmalloc ( sizeof(int)*(seq1_len) ))) {
	goto bail_out;
    }

    if ( ! (last_word = (int *) xmalloc ( sizeof(int)*size_hash ))) {
	goto bail_out;
    }
    if ( ! (word_count = (int *) xmalloc ( sizeof(int)*size_hash ))) {
	goto bail_out;
    }
    
    if ( ! (hash_values2 = (int *) xmalloc ( sizeof(int)*(seq2_len) ))) {
	goto bail_out;
    }
    if ( ! (diag = (int *) xmalloc ( sizeof(int)*(seq1_len+seq2_len) ))) {
	goto bail_out;
    }

    /* hash seq1	*/
    
    if ( hash_seq8 ( &seq1[seq1_lreg-1], hash_values1, seq1_len )  != 0 ) {
	goto bail_out;
    }
    (void) store_hash ( hash_values1, seq1_len, last_word, word_count,
		       size_hash);

    /* hash seq2	*/

    if ( hash_seq8 ( &seq2[seq2_lreg-1], hash_values2, seq2_len )  != 0 ) {
	goto bail_out;
    }

    /* do the comparison	*/
    
    n_matches = compare_seqs ( hash_values1, last_word, word_count, 
				  hash_values2, min_match, seq1_match, 
				  seq2_match, len_match, max_matches, &seq1[seq1_lreg-1], 
				  &seq2[seq2_lreg-1], seq1_len, seq2_len, diag, 
			          seq1_lreg, seq2_lreg, same_seq);

    /* if looking for internal matches remove duplicates */

    /* if ( same_seq && ( n_matches > 0 ) ) { */
    /* 
     * remove duplicates if same seq, n_matches is > 0 and the matches have
     * been saved
     */
    if ( same_seq && ( n_matches > 0 ) ) {
	(void) sip_remdup ( seq1_match, seq2_match, len_match, &n_matches );
    }

    if ( hash_values1 ) free ( hash_values1 );
    if ( hash_values2 ) free ( hash_values2 );
    if ( word_count )   free ( word_count );
    if ( last_word )    free ( last_word );
    if ( diag )         free ( diag );
    return n_matches;

 bail_out:
    if ( hash_values1 ) free ( hash_values1 );
    if ( hash_values2 ) free ( hash_values2 );
    if ( word_count )   free ( word_count );
    if ( last_word )    free ( last_word );
    if ( diag )         free ( diag );
    return -1;
}



int sip_hash (char *seq1,
	      char *seq2,
	      int seq1_lreg, 
	      int seq1_rreg, 
	      int seq2_lreg, 
	      int seq2_rreg, 
	      int max_matches,
	      int min_match,
	      int word_len,
	      int sequence_type,
	      int same_seq,
	      int **seq1_match,
	      int **seq2_match,
	      int **len_match,
	      int *n_matches,
	      void (*pr_func)(void *data, int pw, int y),
	      void *data) 
{
    set_char_set(sequence_type);

    if ( ( DNA == sequence_type ) && ( 8 <= min_match ) ) {

	word_length = 8;

	*n_matches = hash_compare8(min_match, seq1_match, seq2_match, len_match,
			      max_matches, seq1, seq2, seq1_lreg, seq1_rreg,
			      seq2_lreg, seq2_rreg, same_seq, pr_func, data);
	return 0;
    }

    word_length = word_len; /* global assignment */

    *n_matches = hash_compare(min_match, seq1_match, seq2_match, len_match,
			      max_matches, seq1, seq2, seq1_lreg, seq1_rreg,
			      seq2_lreg, seq2_rreg, same_seq, pr_func, data);
    return 0;
}

int quick_scan(char *seq1,		/* seq1 */
	       char *seq2,		/* seq2 */
	       int seq1_lreg,           /* start of seq1 */
	       int seq1_rreg,           /* end of seq1 */
	       int seq2_lreg,           /* start of seq2 */
	       int seq2_rreg,           /* end of seq2 */
	       int sequence_type,
	       int max_matches,	        /* maximum number of matches */
	       int same_seq,		/* if set then its an internal comparison */
	       int win_length,
	       int min_match,	        /* the minimum match length */
	       int word_len,
	       double min_sd,              /* minimum std dev */
	       int save_results,       /* boolean, whether to save matches */
	       int **seq1_match,	/* positions of matches in seq1 */
	       int **seq2_match,	/* positions of matches in seq2 */
	       void (*pr_func)(void *data, int pw, int y),     /* printing function */
	       void *data             /* printing data */
	       ) {

/*	Main interface to the quick scan.

	Jobs are: allocate, hash and compare 2 sequences, deallocate

*/

/*#define SIZE_HASH 65536 */
    static int *hash_values1,*last_word,*word_count,*hash_values2;
    int n_matches,size_hash;
    int *diag;
    int seq1_len, seq2_len;
    hash_values1 = hash_values2 = last_word = word_count = diag = NULL;

    seq1_len = seq1_rreg - seq1_lreg + 1;
    seq2_len = seq2_rreg - seq2_lreg + 1;

    /* sanity checks */

    if ( seq1_len < min_match ) goto bail_out;
    if ( seq2_len < min_match ) goto bail_out;
    if ( seq1_len < win_length ) goto bail_out;
    if ( seq2_len < win_length ) goto bail_out;
    if ( seq1_len < word_len ) goto bail_out;
    if ( seq2_len < word_len ) goto bail_out;

    word_length = word_len; /* global assignment */

    set_char_set(sequence_type);

    set_hash_consts ();

    size_hash = pow(char_set_size-1, word_length);


    /* allocate space for seq1	*/
    if ( ! (hash_values1 = (int *) xmalloc ( sizeof(int)*(seq1_len) ))) {
	goto bail_out;
    }

    if ( ! (last_word = (int *) xmalloc ( sizeof(int)*size_hash ))) {
	goto bail_out;
    }
    if ( ! (word_count = (int *) xmalloc ( sizeof(int)*size_hash ))) {
	goto bail_out;
    }

    if ( ! (diag = (int *) xmalloc ( sizeof(int)*(seq1_len + seq2_len) ))) {
	goto bail_out;
    }
    
    /* hash seq1	*/
    
    if ( hash_seq ( &seq1[seq1_lreg-1], hash_values1, seq1_len )  != 0 ) {
	goto bail_out;
    }
    (void) store_hash ( hash_values1, seq1_len, last_word, word_count,
		       size_hash);
    

    /* allocate space for seq2	*/

    if ( ! (hash_values2 = (int *) xmalloc ( sizeof(int)*(seq2_len) ))) {
	goto bail_out;
    }

    /* hash seq2	*/

    if ( hash_seq ( &seq2[seq2_lreg-1], hash_values2, seq2_len )  != 0 ) {
	goto bail_out;
    }

    /* do the comparison	*/
    
    n_matches = q_compare_seqs (hash_values1, last_word, word_count, 
				hash_values2, diag, &seq1[seq1_lreg-1], 
				&seq2[seq2_lreg-1], seq1_len, seq2_len,
				win_length, min_match, min_sd,
				max_matches, save_results, seq1_match, 
				seq2_match, pr_func, data, seq1_lreg, seq2_lreg);
    if (!save_results) {
	n_matches = -1;
    }

    /* if looking for internal matches remove duplicates */

    if ( same_seq && ( n_matches > 0 ) ) {
	(void) sip_remdup ( seq1_match, seq2_match, NULL, &n_matches);
    }

    if ( hash_values1 ) free ( hash_values1 );
    if ( hash_values2 ) free ( hash_values2 );
    if ( word_count )   free ( word_count );
    if ( last_word )    free ( last_word );
    if ( diag )         free ( diag );
    return n_matches;

 bail_out:
    if ( hash_values1 ) free ( hash_values1 );
    if ( hash_values2 ) free ( hash_values2 );
    if ( word_count )   free ( word_count );
    if ( last_word )    free ( last_word );
    if ( diag )         free ( diag );
    return -1;

}
