/* vector_clip: a program for locating and marking vector segments in readings
   stored in experiment file format. It marks sequencing vector at the 5' end 3'
   ends of readings, checks for vector rearrangements, marks cloning (cosmid)
   vector, and 3' primers.

   At present the code is a bit messy - it uses alignment routines that should
   be in separate files, and some parameters passed between functions are
   probably unnecessary, but in the course of getting things to work have been
   left.

   5-8-98
   Checked out the algorithm for accumulating hits in the histogram by using
   a reading and a vector sequence composed entirely of a's. Changed it so
   instead of filling the sequence up to the cloning site with dashes, we now
   send a subsection of the sequence (from the cloning site). Now appears
   to be correct. ie all diagonals for this test data scored 1.0
   Checked out the algorithm for adding up the words around the best diagonal
   using the same data - again every position got a score of 1.0.
   Problem now is to change from a scoring system to a probability cutoff.

   Also need to deal with primer walks correctly by examining the PR record
   and not looking for the 5' end if PR = 0,3,4

*/
#include <staden_config.h>

#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#include "misc.h" 
#include "dna_utils.h"
#include <io_lib/scf.h>
#include <io_lib/expFileIO.h>
#include "getfile.h"
/*#include "vector_clip.h"*/

#define SEQUENCE_VECTOR 1  /* the jobs to do */
#define CLONING_VECTOR  2
#define HGMP            3
#define VECTOR_REARRANGEMENT 4
#define TEST_ONLY 1
#define TEST_ONLY_ALIGNMENTS 2
#define TRANSPOSON 5
#define SEQUENCE_VECTOR_PVF 6
#define SPECIAL_PRIMER 1
#define USE_QL 1
#define USE_QR 2
#define USE_SL 4
#define USE_SR 8

#define MAX_READ 4096 /* the reading */
#define MAX_VECTOR_D 100000 /* the vector sequence */
#define MIN_VECTOR 4096 /* the vector sequence */
/*#define MAX_PRIMER 100  the primer */
#define MIN_SP 40   /* minimum abs(SP) value for case when SP line 0! */
#define MINMAT_CV 20  /* how close we get to the corners of the matrix */
#define MINMAT_SV 20  /* how close we get to the corners of the matrix */
#define MIN_USEFUL 0 /* after clipping we fail those shorter than this */
#define MAX_MATCHES 1000
#define FUDGE 60	/* fudge factor to allow for crap at the 5' end of reads */

#define MAX_LINE 2048
#define VECTOR_PRIMER_NAME 0
#define F_PRIMER_SEQ 2
#define R_PRIMER_SEQ 1
#define VECTOR_FILE 6
#define MAX_VECTORS 100
#define CLONESITE_SIZE 40 /* in transposon search, alignment to vector
			     must use this length of sequence */
#define TRANS_HASH 20     /* in transposon search, alignment to vector
			     must start within the distance of the 5' end */

/* 9/1/99 johnt - must explictly import globals from DLLs with Visual C++ */
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif


/* extern int *char_lookup;  9/1/99 johnt - Not used, and doesn't link under Windows */
#define MAX_HASH_LENGTH 7
#define MAX_HASH_CONST 200
extern DLL_IMPORT int unknown_char;


int hash_const[MAX_HASH_CONST];

typedef union double_int
    {
	int i;
	double d;
    }DI;

int init_hash8 ( int seq1_len, int seq2_len,
		 int **hash_values1, int **last_word,
		 int **word_count, int **hash_values2,
		 int **diag );

int hash_seq8 ( char *seq, int *hash_values, int seq_len);

int init_hash ( int word_length, int seq1_len, int seq2_len,
	        int **hash_values1, int **last_word,
	        int **word_count, int **hash_values2,
	        int **diag, DI **hist, int *size_hash, int **line );

int hash_seq ( int word_length, char *seq, int *hash_values, int seq_len );

int get_text_seq ( char *seq, int max_len, int *seq_len, FILE *fp );

int do_hash (   int seq1_len, int seq2_len,
	        int *hash_values1, int *last_word,
	        int *word_count, int *hash_values2,
	        int *diag, int *line, int size_hash, DI *hist,
	        int word_length, char *seq1, char *seq2, 
	        double min_score, int num_diags, double *expected_scores, double max_prob,
	        int job,
	        int *x, int *y, double *score_3);

typedef struct vector_spec_ {
  char *name;
  char *f_primer_seq;
  char *r_primer_seq;
  char *file_name;
} Vector_spec;

typedef struct vector_specs_ {
  Vector_spec *vector_spec_ptr;
  int         number;
} Vector_specs;

Vector_specs* get_vector_info(FILE *fp_vf ) {

  char line[MAX_LINE], *file_name, *vector_name = NULL, *f_primer_seq = NULL,*r_primer_seq = NULL;
  char *c, *s;
  Vector_spec *v;
  static Vector_specs vs;
  int item, in_item, number;

  if (NULL == (v = (Vector_spec *)malloc((sizeof(Vector_spec ) * MAX_VECTORS)))) return NULL;

    /* format: 
     *
       name seq_r seq_f file_name 
     *
     * file_name is optional
     * file_name might contain spaces! but none of the others can.
     */ 

  number = 0;
  while (fgets(line,MAX_LINE,fp_vf) != NULL) {
    /* remove newline */
/*    line[strlen(line)-1] = '\0';*/
    for(c = line, in_item = 0, item = 0; item < 3; item++ ) {
      if ( !(in_item ) ) {
	/* looking for an item */
	for (;*c && isspace(*c); c++);
	if (*c) {
	  /* in an item */
	  in_item = 1;
	  for (s=c;*c && !isspace(*c); c++) ;
	  if (*c) {
	    /* found end of item */
	    in_item = 0;
	    *c++ = '\0';
	    if ( item == VECTOR_PRIMER_NAME ) {
	      if (NULL == (vector_name = (char *)malloc((strlen(s) + 1)*sizeof(char)))) return NULL;
	      strcpy ( vector_name, s );
	    }
	    else if ( item == F_PRIMER_SEQ ) {
	      if (NULL == (f_primer_seq = (char *)malloc((strlen(s) + 1)*sizeof(char)))) return NULL;
	      strcpy ( f_primer_seq, s );
	    }
	    else if ( item == R_PRIMER_SEQ ) {
	      if (NULL == (r_primer_seq = (char *)malloc((strlen(s) + 1)*sizeof(char)))) return NULL;
	      strcpy ( r_primer_seq, s );
	    }
	  }
	  else {
	    /* item ends oddly */
	    printf("Warning: error parsing vector_primer file\n");
	    break;
	  }
	}
	else {
	  /* item starts oddly */
	  printf("Warning: error parsing vector_primer file\n");
	  break;
	}
      }
    }
    if ( item != 3 ) continue;
    /* if the file name is present, get it */
    file_name = NULL;
    if ( strlen ( c ) ) {
      c[strlen(c)-1] = '\0';
      if (NULL == (file_name = (char *)malloc((strlen(c) + 1)*sizeof(char)))) return NULL;
      strcpy ( file_name, c );
    }

    if ( number == MAX_VECTORS ) {
      printf("Aborting: more than %d entries in vector_primer file\n",MAX_VECTORS);
      return NULL;
    }

    complement_seq ( f_primer_seq, strlen( f_primer_seq ));
    v[number].name  = vector_name;
    v[number].f_primer_seq = f_primer_seq;
    v[number].r_primer_seq = r_primer_seq;
    v[number].file_name = file_name;
    number++;
  }

  /* some of the segments between primer and cloning site may be too short
   * to give reliable matching for the 3' end searches, so here we extend them
   * FIXME given up on this!
   */

  vs.vector_spec_ptr = v;
  vs.number = number;
/*
      for(i=0;i<number;i++) {
	printf(" record %d\n",i);
	printf("%s\n",v[i].name);
	printf("%s\n",v[i].r_primer_seq);
	printf("%s\n",v[i].f_primer_seq);
	if(v[i].file_name)printf("%s\n",v[i].file_name);
      }
*/
  return &vs;
}

void write_dot(void) {
    fprintf ( stdout, "." );
    (void) fflush ( stdout );
}

int vep_error ( FILE *fp, char *file_name, int error_no ) {
    char *err_mess[] = {
/* 1 */	"Error: could not open experiment file",
/* 2 */	"Error: no sequence in experiment file",
/* 3 */	"Error: sequence too short",
/* 4 */	"Error: missing vector file name",
/* 5 */	"Error: missing cloning site",
/* 6 */	"Error: missing primer site",
/* 7 */	"Error: could not open vector file",
/* 8 */	"Error: could not write to experiment file",
/* 9 */	"Error: could not read vector file",
/* 10 */ "Error: missing primer sequence",
/* 11 */ "Error: hashing problem",
/* 12 */ "Error: alignment problem",
/* 13 */ "Error: invalid cloning site",
/* 14 */ "Error: sequence now too short",
/* 15 */ "Error: sequence entirely cloning vector",
/* 16 */ "Error: possible vector rearrangement",
/* 17 */ "Error: invalid sequence for demonstration mode"};
     

    fprintf ( stderr, "%s\n",err_mess[error_no-1] );
    if ( fp ) fprintf ( fp, "%s %s\n",file_name,err_mess[error_no-1]);
    fprintf ( stdout, "!" );
    (void) fflush ( stdout );

    return 0;
}

int poisson_diagonals ( int min_diag, int max_diag, int word_size,
		       double max_prob, double *expected_scores ) {
  int diagonal_length, hits;
  double expected_hits, prob_0_hits, fact, sum_probs;
  double prob_remaining;

  /* Assume the diagonal scores obey a poisson distribution
   * with an expected number of hits given by 
   *                                word_size
   *             diagonal_length / 4
   */

  for(hits=0;hits<min_diag;hits++) {
    expected_scores[hits] = 1.0;
  }

  /* for each diagonal length find the required score */

  for(diagonal_length=min_diag;diagonal_length<max_diag;diagonal_length++) {
    expected_hits = (double)diagonal_length/pow(4.0,word_size);
    prob_0_hits = exp(-1*expected_hits);

    /* sum the probabilities until a higher number of hits is
     * sufficiently improbable
     */

    for(hits=1,fact=1.0,sum_probs=prob_0_hits;hits<max_diag;hits++) {
      fact *= hits;
      sum_probs += prob_0_hits * pow(expected_hits,hits)/fact; 
      prob_remaining = 1.0 - sum_probs;

      if ( prob_remaining < max_prob ) {
	expected_scores[diagonal_length] = ((double)word_size*((double)hits-0.5))/(double)diagonal_length;
/*	printf("%d %f\n",diagonal_length,expected_scores[diagonal_length]); */
	break;
      }
    }
  }
  return 0;
}

int poisson_diagonal ( int diagonal_length, int word_size,
		       double score_in ) {
  int hits;
  double expected_hits, prob_0_hits, fact, sum_probs;
  double prob_remaining, min_prob=1.0e-20;
  double score_f, score;
  
  /* calculates the probability for a word length, a score and a diagonal */

  expected_hits = (double)diagonal_length/pow(4.0,word_size);
  prob_0_hits = exp(-1*expected_hits);
  score_f = (double)word_size/(double)diagonal_length;

  for(hits=1,fact=1.0,sum_probs=prob_0_hits;hits<100000;hits++) {
    fact *= hits;
    sum_probs += prob_0_hits * pow(expected_hits,hits)/fact; 
    prob_remaining = 1.0 - sum_probs;
    score = score_f*(double)(hits-0.5);
    if ( score > score_in ) {
      printf("score %f prob %E\n",score_in,prob_remaining);
      return 0;
    }
    if ( prob_remaining < min_prob ) {
      printf("error: hits %d prob %E\n",hits,prob_remaining);
      return -1;
    }
  }
  return 1;
}


void free_hash ( int *hash_values1, int *last_word,
	        int *word_count, int *hash_values2, int *diag,
	        DI *hist) {

    if ( hash_values1 ) xfree ( hash_values1 );
    if ( hash_values2 ) xfree ( hash_values2 );
    if ( word_count )   xfree ( word_count );
    if ( last_word )    xfree ( last_word );
    if ( diag )         xfree ( diag );
    if ( hist )         xfree ( hist );

}

/************************************************************/

void store_hash ( 
	         int *hash_values, 	/* the hash values for each position in a seq */
		 int seq_len, 		/* size of the seq and hash array */
		 int *last_word, 	/* last occurrence of this hash value */
		 int *word_count, 	/* frequency of each hash value or word */
		 int word_length,       /* word length */
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

/* 	zero the word counts */

    for ( i=0;i<size_wc;i++ ) {
	word_count[i] = 0;
	last_word[i] = 0;
    }

/*	loop for all entries in hash_values	*/

    j = seq_len - word_length + 1;

    for ( i = 0; i < j; i++ ) {

	n = hash_values[i];

/*	is it a good value ? */

	if ( -1 != n ) {

	    nw = word_count[n];

/* 		already an entry for this word ? */

	    if ( 0 == nw ) {

/*		no, so put in last_word */

		last_word[n] = i;
		word_count[n] += 1;
	    }

/*		yes, so put previous last occurrence in hash_values*/

	    else {

		word_count[n] += 1;
		hash_values[i] = last_word[n];
		last_word[n] = i;
	    }
	}
    }
}

/***************************************************************************/

/*
 * Needleman and Wunsch routine: give two sequences returns a full length alignment. 
 * The alignment arrays must be at least of size  seq1_len + seq2_len + 1 just in case
 * we end up with an alignment like:  abc***
 *                                    ***def
 * The output strings are terminated with '\0' and their common length is returned.
 * The routine ecconomises on space by using only two int arrays of size seq1_len
 * (as we only need to know what happened in the last row and to build the next),
 * plus a trace array of size (seq1_len * seq2_len)/4 chars. This is possible because
 * we only have 3 choices for moving: across, down or diagonally and so we encode these
 * in 2 bits. This reduces the storage by a factor of 16 over using a byte for each
 * element so is worth the complication.
 */

int nw ( char seq1[], 
	char seq2[], 
	int seq1_len, 
	int seq2_len,
	char *seq1_res,
	char *seq2_res,
	int *len_align
	) {
    int i,j, gap_pen=-2;
    int size_mat, s_c, s_r, row, column, s, e;
    int b_s, b_e, b_r, b_c, r, c,*rows1,*rows2,*rowsp1,*rowsp2;
    int r_u, r_l, r_d, t,max_seq,byte,nibble;
    char *seq1_out,*seq2_out;
    unsigned char *bit_trace, byte_across=1, byte_down=2, byte_diagonal=3;
    unsigned char mask = 3;
    int s_matrix[]={	3,0,0,0,0,
			0,3,0,0,0,
			0,0,3,0,0,
                        0,0,0,3,0,
                        0,0,0,0,0};

    max_seq = seq1_len + seq2_len + 1;
    size_mat = (seq1_len + 1) * (seq2_len + 1);

    if ( ! (rows1 = (int *) xmalloc ( sizeof(int)*(seq1_len+1) ))) return -1;
    if ( ! (rows2 = (int *) xmalloc ( sizeof(int)*(seq1_len+1) ))) {
	xfree ( rows1 );
	return -1;
    }
    if ( ! (bit_trace = (unsigned char *) xmalloc ( 1 + sizeof(char)*size_mat/4 ))) {
	xfree ( rows1 );
	xfree (rows2 );
	return -1;
    }

    /* initialise our own extra copy of the output pointers */

    seq1_out = seq1_res;
    seq2_out = seq2_res;

    j = 1 + size_mat/4;
    for(i=0;i<j;i++) {
	bit_trace[i] = 0;
    }
    for ( i=0;i<=seq1_len;i++) {
/*
	rows1[i] = gap_pen;
	rows2[i] = gap_pen;
*/

	rows1[i] = 0;
	rows2[i] = 0;
    }

    b_s = b_e = b_r = b_c = 0;
    t = 0;

    /* proceed row by row */

    for(row=1,e=seq1_len+2;row<=seq2_len;row++,e++) {

	/* switch rows: rowsp1 points to prev row, rowsp2 to new row */

	if ( t == 0 ) {
	    rowsp1 = rows1, rowsp2 = rows2+1;
	    t = 1;
	}
	else {
	    rowsp1 = rows2, rowsp2 = rows1+1;
	    t = 0;
	}
	s_r = char_lookup[seq2[row-1]];
	s_r *= 5;

	/* for each row do each column */

	for(column=1;column<=seq1_len;column++,e++,rowsp1++,rowsp2++) {

	    s_c = char_lookup[seq1[column-1]];
	    r_d = *rowsp1;
	    r_u = *(rowsp1+1);
	    r_l = *(rowsp2-1);

	    s = s_matrix[s_r+s_c];

	    /* choose a direction:

            e = element we are at
	    d = prev_row,prev_col
	    u = prev_row,this_column = d+1
	    l = prev_column,this_row = e-1
	    want to find max (r_d + s, r_u + gap, r_l + gap) 
	    */

	    byte = e/4, nibble = 2*(e%4);
	    r_d += s;
	    r_u += gap_pen;
	    r_l += gap_pen;
	    if (( r_d >= r_u ) & ( r_d >= r_l )) {
		bit_trace[byte] |= byte_diagonal << nibble;
		*rowsp2 = r_d;
		if ( r_d > b_s ) {
		    b_s = r_d;
		    b_e = e;
		    b_r = row;
		    b_c = column;
		}
	    }
	    else if ( r_u >= r_l ) {
		bit_trace[byte] |= byte_down << nibble;
		*rowsp2 = r_u;
		if ( r_u > b_s ) {
		    b_s = r_u;
		    b_e = e;
		    b_r = row;
		    b_c = column;
		}
	    }
	    else {
		bit_trace[byte] |= byte_across << nibble;
		*rowsp2 = r_l;
		if ( r_l > b_s ) {
		    b_s = r_l;
		    b_e = e;
		    b_r = row;
		    b_c = column;
		}
	    }

	}
    }

    /* now use the bit_trace to create the alignment */

    for(i=0;i<max_seq-1;i++) {
	*(seq1_out++) = '*';
	*(seq2_out++) = '*';
    }
    *seq1_out = *seq2_out = '\0';
    seq1_out--;
    seq2_out--;

    /* do any right hand end overhang */

    r = seq2_len-1;
    c = seq1_len-1;
    i = seq2_len-b_r - ( seq1_len-b_c );

    if (i > 0 ) {
	j = i;
	while (j>0) {
	    *(seq2_out--) = seq2[r];
	    seq1_out--;
	    r--;
	    j--;
	}
    }
    else if (i < 0 ) {
	j = -i;
	i = j;
	while (j>0) {
	    *(seq1_out--) = seq1[c];
	    seq2_out--;
	    c--;
	    j--;
	}
    }

    /* do up to best row and col */

    while(r>=b_r) {
	*(seq2_out--) = seq2[r];
	*(seq1_out--) = seq1[c];
	r--;
	c--;
    }

    /* follow the path for the middle section */

    r = b_r, c = b_c, e = b_e;

    while ( (r>0)&&(c>0)) {

	byte = e/4, nibble = 2*(e%4);
	if ( ((bit_trace[byte] >> nibble)&mask) == byte_diagonal ) {

	    r--,c--;
	    *(seq1_out--) = seq1[c];
	    *(seq2_out--) = seq2[r];
	}
	else if ( ((bit_trace[byte] >> nibble)&mask) == byte_down ) {

	    r--;
	    seq1_out--;
	    *(seq2_out--) = seq2[r];
	}
	else {

	    c--;
	    seq2_out--;
	    *(seq1_out--) = seq1[c];
	}
	e=(seq1_len+1)*r + c;
    }

    /* finish off the left ends dealing with any overhang */

    while(r>0) {
	r--;
	*(seq2_out--) = seq2[r];
    }

    while(c>0) {
	c--;
	*(seq1_out--) = seq1[c];
    }

    /* now, if necessary move all the data left
       to remove pads at the left end */

    c = strlen ( seq1_res );
    e = strlen ( seq2_res );
    for ( i=0;i<c;i++ ) {
	if ( (seq1_res[i] != '*') || (seq2_res[i] != '*') ) break;
    }
    r = i;
    for ( i=0,j=r;j<=c;i++,j++ ) {
	seq1_res[i] = seq1_res[j];
	seq2_res[i] = seq2_res[j];
    }
    *len_align = i - 1;
    xfree ( rows1 );
    xfree ( rows2 );
    xfree ( bit_trace );
    return 0;
}

int alignment_score ( char seq1[], 
		     char seq2[], 
		     int seq1_len, 
		     int seq2_len,
		     int job,
		     double *percent_match,
		     int *left_end,
		     int *right_end,
		     int *pad1,
		     int *pad2,
		     int *left1,
		     int *right1,
		     int *left2,
		     int *right2,
		     int tmode
		     ) {

    /* returns 4 types of alignment percentage match score and pad counts:

       **abc**
       sf*bcds
       1234567

       the match is only over the bit in the middle but sometimes we want
       to express the percentage mismatch without including the end pads.

       1 = full length , ie 1-7
       2 = unpadded part, ie 3-5
       3 = all except left padded region, ie 3-7
       4 = all except right padded region, ie 1-5
       for each the pads are counted over the same region as the percentage

    */

    int i, j, k, max_seq, len_align, ret, score, left, right, start, end;
    char *seq1_out, *seq2_out, *match;

    if ( (seq1_len < 1) || (seq2_len < 1) ) return -2;
    max_seq = seq1_len + seq2_len + 1;
    if ( ! (seq1_out = (char *) xmalloc ( sizeof(char)*max_seq ))) return -1;
    if ( ! (seq2_out = (char *) xmalloc ( sizeof(char)*max_seq ))) return -1;
    if ( ! (match = (char *) xmalloc ( sizeof(char)*81 ))) return -1;


    ret = nw ( seq1, 
	seq2, 
	seq1_len, 
	seq2_len,
	seq1_out,
	seq2_out,
	&len_align
	);

    if (  ret ) {
	xfree ( seq1_out );
	xfree ( seq2_out );
	xfree ( match );
	return ret;
    }

/*
    printf("%s\n",seq1_out);
    printf("%s\n",seq2_out);
*/

    for ( i=0;i<len_align;i++ ) {
	if ( seq1_out[i] != '*') break;
    }
    *left1 = i;
    for ( i=0;i<len_align;i++ ) {
	if ( seq2_out[i] != '*') break;
    }
    *left2 = i;
    left = MAX(*left1,*left2);

    for ( i=len_align-1;i>=0;i-- ) {
	if ( seq1_out[i] != '*') break;
    }
    *right1 = i;
    for ( i=len_align-1;i>=0;i-- ) {
	if ( seq2_out[i] != '*') break;
    }
    *right2 = i;
    right = MIN(*right1,*right2)+1;

    if ( 1 == job ) {
	start = 0, end = len_align;
    }
    else if ( 2 == job ) {
	start = left, end = right;
    }
    else if ( 3 == job ) {
	start = left, end = len_align;
    }
    else if ( 4 == job ) {
	start = 0, end = right;
    }
    else {
	return -3;
    }
    *pad1 = *pad2 = 0;

    for ( i=start,score=0;i<end;i++ ) {
	    if ( SEQ_MATCH(seq1_out[i],seq2_out[i]) ) score++;
	    if ( '*' == seq1_out[i] ) (*pad1)++;
	    if ( '*' == seq2_out[i] ) (*pad2)++;
	}

    if ( end - start ) {
	*percent_match = 100.0 * score / ( end - start );
	ret = 0;
    }
    else {
	ret = -4;
    }
    if ( tmode == TEST_ONLY_ALIGNMENTS ) {

      i = MIN((end),(start+79));
	k = MIN(79,end-start);
	match[k] = seq1_out[i] = seq2_out[i] = '\0';
	for(j=0;j<k;j++) {
	  match[j] = ' ';
	  if ( SEQ_MATCH(seq1_out[j+start],seq2_out[j+start]) ) match[j] = ':';
	}
	printf("%s\n",&seq1_out[start]);
	printf("%s\n",match);
	printf("%s\n",&seq2_out[start]);
    }
    xfree ( seq1_out );
    xfree ( seq2_out );
    xfree ( match );
    *left_end = start;
    *right_end = end;
    return ret;
}

/****************************************************************************/



int do_it_3p ( FILE *fp_i, FILE *fp_p, FILE *fp_f,
	       double cut_score,
	       int tmode) {

    char file_name[FILENAME_MAX+1],*cp, *seq;
/*    char primer[MAX_PRIMER] = "ATGAGCTCTGAGTCGTGCTATCCTG"; */
    char *primer;
    int left1,left2,right1,right2;
    Exp_info *e;
    int ql,qr,seq_length;
    char buf[10];

/* for this algorithm */

    int primer_length, x, y, ret, eret;
    int pad1, pad2;
    double percent_match;
 
    while ( fgets ( file_name, FILENAME_MAX, fp_i )) {

	if ( cp = strchr ( file_name, '\n' )) *cp = '\0';
	if ( tmode ) {
	    printf(">>>>>>>>>>>>>>>>>>>>>> %s\n", file_name );
	}

	e = exp_read_info ( file_name );
	if ( e == NULL ) {
	    eret =  vep_error ( fp_f, file_name, 1 );
	}
	else {
	    if ( exp_Nentries ( e, EFLT_SQ ) < 1 ) {
		eret =  vep_error ( fp_f, file_name, 2 );
		exp_destroy_info ( e );
		continue;
	    }
	    else {
		char *expline;

		seq = exp_get_entry ( e, EFLT_SQ );
		seq_length = strlen ( seq );

		if ( exp_Nentries ( e, EFLT_PD )) {
		    primer = exp_get_entry ( e, EFLT_PD );
		    primer_length = strlen(primer);
		}
		else {
		    eret =  vep_error ( fp_f, file_name, 10 );
		    exp_destroy_info ( e );
		    continue;
		}

		ql = 0;
		qr = seq_length;
		if ( exp_Nentries ( e, EFLT_QL )) {
		    expline = exp_get_entry ( e, EFLT_QL );
		    ql = atoi ( expline );
		}
		if ( exp_Nentries ( e, EFLT_QR )) {
		    expline = exp_get_entry ( e, EFLT_QR );
		    qr = atoi ( expline );
		}
		if ( qr < MIN_USEFUL ) {
		    eret =  vep_error ( fp_f, file_name, 3 );
		    exp_destroy_info ( e );
		    continue;
		}

		/* we are searching the whole of the reading to the right of ql
		   which slows the search but means we should find it where ever
		   it is.
		*/

		ret = alignment_score ( &seq[ql], 
		     primer, 
		     seq_length-ql,
		     primer_length,
		     2,
		     &percent_match,
		     &x,
		     &y,
		     &pad1,
                     &pad2,
		     &left1,
		     &right1,
		     &left2,
		     &right2,
		     tmode
		     );
		x += ql+1;
		if ( ret < 0 ) {
		    eret =  vep_error ( fp_f, file_name, 12 );
		    exp_destroy_info ( e );
		    continue;
		}
		if ( !tmode ) {
		    if ( percent_match >= cut_score ) {
			sprintf(buf, "%d", x);

			/* we write an sr line for this search */

			if (exp_put_str(e, EFLT_SR, buf, strlen(buf))) {
			    eret =  vep_error ( fp_f, file_name, 8 );
			    exp_destroy_info ( e );
			    continue;
			}
			else {
			    if ( fp_p ) fprintf ( fp_p, "%s\n",file_name);
			}
		    }
		    else {
			if ( fp_p ) fprintf ( fp_p, "%s\n",file_name);
		    }
		}
		else {
		    if ( percent_match >= cut_score ) {
			printf("match at %d with score %f\n",x,percent_match);
		    }
		    else {
			printf("no match\n");
		    }
		}

	    }

	}
	exp_destroy_info ( e );
	if (!(tmode)) (void) write_dot();
    }

    return 0;
}

void seq_to_dash ( char seq[], int len ) {
    int i;
    for(i=0;i<len;i++) {
	seq[i] = '-';
    }
}


int do_it_sv ( char *vector_seq, int max_vector, 
	       FILE *fp_i, FILE *fp_p, FILE *fp_f,
	       double cut_score_5, double cut_score_3, int min_left,
	       int tmode, int rmode) {

    char file_name[FILENAME_MAX+1],*cp, *seq;

    Exp_info *e;
    int ql,qr,seq_length;
    char buf[10];

/* for this algorithm */

    int vector_length = 0, x, y, ret, eret, retm;
    double score_3, score_5;
    char vector_file_name[FILENAME_MAX+1];
    char prev_vector_file_name[FILENAME_MAX+1];
    FILE *vf;
    int sp=0, sc=0; /* primer position and cloning site */
    int prev_sp, prev_sc; /* previous values */
    int new_vector;
    int sl, sr;
    int length_s;
    int length_v;
    int pad1, pad2;
    int left1,left2,right1,right2;
    int pr, primer_type;
    char *vfn;

    /* initialise this algorithm */

    prev_sp = prev_sc = 0;
    prev_vector_file_name[0] = '\0';

    while ( fgets ( file_name, FILENAME_MAX, fp_i )) {

	if ( cp = strchr ( file_name, '\n' )) *cp = '\0';

	if ( tmode ) {
	    printf(">>>>>>>>>>>>>>>>>>>>>> %s\n", file_name );
	}

	e = exp_read_info ( file_name );

	if ( e == NULL ) {
	    eret =  vep_error ( fp_f, file_name, 1 );
	    continue;
	}
	else {
	    if ( exp_Nentries ( e, EFLT_SQ ) < 1 ) {
		eret =  vep_error ( fp_f, file_name, 2 );
		exp_destroy_info ( e );
		continue;
	    }
	    else {

		char *expline;

		seq = exp_get_entry ( e, EFLT_SQ );
		seq_length = strlen ( seq );

		if ( seq_length < MIN_USEFUL ) {
		    eret =  vep_error ( fp_f, file_name, 3 );
		    exp_destroy_info ( e );
		    continue;
		}
		ql = 0;
		qr = seq_length + 1;

		if ( rmode & USE_QL ) {
		  if ( exp_Nentries ( e, EFLT_QL )) {
		    expline = exp_get_entry ( e, EFLT_QL );
		    ql = atoi ( expline );
		  }
		}
		if ( rmode & USE_QR ) {
		  if ( exp_Nentries ( e, EFLT_QR )) {
		    expline = exp_get_entry ( e, EFLT_QR );
		    qr = atoi ( expline );
		  }
		}
		sl = ql;
		sr = qr;
		if ( rmode & USE_SL ) {
		  if ( exp_Nentries ( e, EFLT_SL )) {
		    expline = exp_get_entry ( e, EFLT_SL );
		    sl = atoi ( expline );
		  }
		}
		if ( rmode & USE_SR ) {
		  if ( exp_Nentries ( e, EFLT_SR )) {
		    expline = exp_get_entry ( e, EFLT_SR );
		    sr = atoi ( expline );
		  }
		}
		if ( exp_Nentries ( e, EFLT_SC )) {
		    expline = exp_get_entry ( e, EFLT_SC );
		    sc = atoi ( expline );
		}
		else {
		    eret =  vep_error ( fp_f, file_name, 5 );
		    exp_destroy_info ( e );
		    continue;
		}

		/* internal primer ? */

		pr = 0;
		primer_type = 0;
		score_5 = 0.0;
		if ( exp_Nentries ( e, EFLT_PR )) {
		    expline = exp_get_entry ( e, EFLT_PR );
		    pr = atoi ( expline );
		    if ( (pr == 3) || (pr == 4) || (pr==0)) primer_type = SPECIAL_PRIMER;
		    /* doing this test here means PR line not obligatory
		     * but we can specify PR 0 to skip the primer search
		     * (the 3' end search is done even though we do not
		     * know the orientation of the vector! so hopefully
		     * we dont find a match! Really PR 0 should not pass
		     * through vector_clip
		     */
		}
		if ( exp_Nentries ( e, EFLT_SP )) {
		    expline = exp_get_entry ( e, EFLT_SP );
		    sp = atoi ( expline );
		}
		else {
		  if ( primer_type != SPECIAL_PRIMER ) {
		    eret =  vep_error ( fp_f, file_name, 6 );
		    exp_destroy_info ( e );
		    continue;
		  }
		}

		new_vector = 0;
		if (( sp != prev_sp) || ( sc != prev_sc)) new_vector = 1; 
		prev_sp = sp;
		prev_sc = sc;

		/* get the vector file */

		if ( exp_Nentries ( e, EFLT_SF )) {
		    vfn = exp_get_entry(e,EFLT_SF);
		    if (1 != expandpath(vfn, vector_file_name)) {
		      eret =  vep_error ( fp_f, file_name, 4 );
		      exp_destroy_info ( e );
		      continue;
		    }
		}
		else {
		    eret =  vep_error ( fp_f, file_name, 4 );
		    exp_destroy_info ( e );
		    continue;
		}

		/* we only get the vector seq if it is different from
		   the previous one, or the cloning site is different,
		   or the primer site is different (actually this only
		   matters if has moved to the other strand)
		*/

		if ( strcmp ( vector_file_name, prev_vector_file_name )) new_vector = 1;

		if ( new_vector ) {
		    vf = fopen(vector_file_name, "r");
		    if (vf == NULL ) {
			eret =  vep_error ( fp_f, file_name, 7 );
			exp_destroy_info ( e );
			continue;
		    }

		    ret = get_text_seq ( vector_seq, max_vector, &vector_length, vf);

		    fclose(vf);
		    if ( ret ) {
			eret =  vep_error ( fp_f, file_name, 9 );
			exp_destroy_info ( e );
			continue;
		    }

		    /* rotate sequence so that cloning site at vector_seq[0] */

		    ret = rotate_seq ( vector_seq, vector_length, sc+1);
		    if ( ret ) {
			eret =  vep_error ( fp_f, file_name, 13 );
			exp_destroy_info ( e );
			continue;
		    }

		    /* for forward primers we need to complement the vector */

		    if ( (sp > 0) || ( pr == 3 ) ) complement_seq ( vector_seq, vector_length );

		    strcpy( prev_vector_file_name, vector_file_name );
		}


		if ( primer_type != SPECIAL_PRIMER ) { /* do not search for primers if walking */
		  length_v = MIN((abs(sp)),vector_length);
		  length_s = MIN((abs(sp)+FUDGE),sr-sl-1);

		  ret = alignment_score ( &seq[sl], 
		     &vector_seq[vector_length-length_v], 
		     length_s,
		     length_v,
		     2,
		     &score_5,
		     &x,
		     &y,
		     &pad1,
                     &pad2,
		     &left1,
		     &right1,
		     &left2,
		     &right2,
		     tmode
		     );

		  /*printf("left1 %d right1 %d left2 %d right2 %d\n",
		       left1,right1,left2,right2);*/

		  x = MIN(right1,right2) - left1 - pad1 + 1;
		  if ( ret < 0 ) {
		    eret =  vep_error ( fp_f, file_name, 12 );
		    exp_destroy_info ( e );
		    continue;
		  }

		  if ( score_5 >= cut_score_5 ) {
		    sl += x;
		  }
		  else {
		    /* if min_left is set we use it if no match has been found.
		       ie it is not exactly a minimum but rather a default
		       Note we also allow a -ve value of min_left to mean
		       "use the primer length as a min_left value"
		     */
		    if ( min_left ) {
		      if ( -1 == min_left ) {
			sl += abs ( sp );
		      } else {
			sl += min_left;
		      }
		      /* make sure this does not go off the end! */
		      sl = MIN (sl, sr-1);
		    }
		  }
		}
		length_s = sr - sl - 1;
		length_v = MAX(abs(sp),MIN_SP);
		ret = alignment_score ( &seq[sl], vector_seq, length_s,
		      length_v, 2, &score_3, &x, &y, &pad1, &pad2,
		      &left1, &right1, &left2, &right2, tmode );

		if ( ret < 0 ) {
		  eret =  vep_error ( fp_f, file_name, 12 );
		  exp_destroy_info ( e );
		  continue;
		}

		if ( score_3 >= cut_score_3 ) sr = sl + x + 1;

		if ( !tmode ) {
		  sprintf(buf, "%d", sl);
		  if (exp_put_str(e, EFLT_SL, buf, strlen(buf))) {
		    eret =  vep_error ( fp_f, file_name, 8 );
		    exp_destroy_info ( e );
		    continue;
		  }
		  sprintf(buf, "%d", sr);
		  if (exp_put_str(e, EFLT_SR, buf, strlen(buf))) {
		    eret =  vep_error ( fp_f, file_name, 8 );
		    exp_destroy_info ( e );
		    continue;
		  }
		  if ( sr - sl < MIN_USEFUL ) {
		    eret =  vep_error ( fp_f, file_name, 14 );
		    exp_destroy_info ( e );
		    continue;
		  }
		  else {
		    if ( fp_p ) fprintf ( fp_p, "%s\n",file_name);
		  }
		}
		else {
		  retm = 0;
		  if ( score_5 >= cut_score_5 ) {
		    printf("5' match at %d with score %f\n",sl,score_5);
		    retm = 1;
		  }
		  if ( score_3 >= cut_score_3 ) {
		    printf("3' match at %d with score %f\n",sr,score_3);
		    retm = 1;
		  }
		  if ( retm == 0 ) {
		    printf("no match\n");
		  }
		}
	    }
	}
	exp_destroy_info ( e );
	if (!(tmode)) (void) write_dot();
    }
    return 0;
}


int do_it_sv_pvf ( char *vector_seq, int max_vector, 
		   FILE *fp_i, FILE *fp_p, FILE *fp_f,
		   double cut_score_5, double cut_score_3,
		   int min_left, int vector_primer_length,
		   int tmode, int rmode, Vector_specs *vs) {

    char file_name[FILENAME_MAX+1],*cp, *seq;

    Exp_info *e;
    int ql,qr,seq_length;
    char buf[10];

    /* for this algorithm */

    int x, y, ret, eret, retm;
    double score_3, score_5;
    char *seq_ptr, *vf_ptr, *sub_seq_ptr;
    char pstat[1024];
    int sp; /* primer length */
    int sl, sr;
    int primer_length_s;
    int primer_length_v;
    int primer_d;
    int pad1, pad2;
    double percent_match;
    int left1,left2,right1,right2;
    int primer_num,best_match=0,best_match_d=0, best_match_pos=0;
    int best_match_5,best_match_d5, pr_type;
    double best_percent_match;
    int length_v, length_s;

    while ( fgets ( file_name, FILENAME_MAX, fp_i )) {

	if ( cp = strchr ( file_name, '\n' )) *cp = '\0';

	if ( tmode ) {
	    printf(">>>>>>>>>>>>>>>>>>>>>> %s\n", file_name );
	}

	e = exp_read_info ( file_name );

	if ( e == NULL ) {
	    eret =  vep_error ( fp_f, file_name, 1 );
	    continue;
	}
	else {
	    if ( exp_Nentries ( e, EFLT_SQ ) < 1 ) {
		eret =  vep_error ( fp_f, file_name, 2 );
		goto next_file;
	    }
	    else {

		char *expline;

		seq = exp_get_entry ( e, EFLT_SQ );
		seq_length = strlen ( seq );
	      
		if ( seq_length < MIN_USEFUL ) {
		    eret =  vep_error ( fp_f, file_name, 3 );
		    goto next_file;
		}
		ql = 0;
		qr = seq_length + 1;
	      
		if ( rmode & USE_QL ) {
		    if ( exp_Nentries ( e, EFLT_QL )) {
			expline = exp_get_entry ( e, EFLT_QL );
			ql = atoi ( expline );
		    }
		}
		if ( rmode & USE_QR ) {
		    if ( exp_Nentries ( e, EFLT_QR )) {
			expline = exp_get_entry ( e, EFLT_QR );
			qr = atoi ( expline );
		    }
		}
		sl = ql;
		sr = qr;
		if ( rmode & USE_SL ) {
		    if ( exp_Nentries ( e, EFLT_SL )) {
			expline = exp_get_entry ( e, EFLT_SL );
			sl = atoi ( expline );
		    }
		}
		if ( rmode & USE_SR ) {
		    if ( exp_Nentries ( e, EFLT_SR )) {
			expline = exp_get_entry ( e, EFLT_SR );
			sr = atoi ( expline );
		    }
		}
	      
		/* try all primer pairs */

		pr_type = 0; /* uknown */
		score_5 = score_3 = best_percent_match = -99.0;
		best_match_5 = best_match_d5 = -9;
		vf_ptr = NULL;
		for (primer_num = 0; primer_num < vs->number; primer_num++) {
		    seq_ptr = vs->vector_spec_ptr[primer_num].f_primer_seq;
		    for (primer_d = 0; primer_d < 2; primer_d++) {
			sp = strlen (seq_ptr);
			primer_length_v = MIN(abs(sp),vector_primer_length);
			primer_length_s = MIN((abs(sp)+FUDGE),sr-sl-1);
			sub_seq_ptr = seq_ptr + strlen(seq_ptr) - primer_length_v;
			ret = alignment_score ( &seq[sl], sub_seq_ptr, 
						primer_length_s, primer_length_v,
						2, &percent_match, &x, &y, &pad1, &pad2,
						&left1, &right1, &left2, &right2, tmode );

			if ( percent_match > best_percent_match ) {
			    best_percent_match = percent_match;
			    best_match = primer_num;
			    best_match_d = primer_d;
			    best_match_pos = MIN(right1,right2) - left1 - pad1 + 1;
			}
			seq_ptr = vs->vector_spec_ptr[primer_num].r_primer_seq;
		    }
		}

		if ( best_percent_match >= cut_score_5 ) {

		    sl += best_match_pos;
		    score_5 = best_percent_match;

		    /* note the details so we can compare with 3' end! */
		
		    best_match_5 = best_match;
		    best_match_d5 = best_match_d;
		    if ( best_match_d == 0 ) pr_type = 1;
		    if ( best_match_d == 1 ) pr_type = 2;
		    vf_ptr = vs->vector_spec_ptr[best_match_5].file_name;
		    sprintf(pstat, "vector_clip: %s",
			    vs->vector_spec_ptr[best_match_5].name);
		}

		/* Weve done the 5' end so now do 3' */

		/* try all primer pairs */

		best_percent_match = -99.0;
		for (primer_num = 0; primer_num < vs->number; primer_num++) {
		    seq_ptr = vs->vector_spec_ptr[primer_num].f_primer_seq;
		    for (primer_d = 0; primer_d < 2; primer_d++) {
			length_v = strlen(seq_ptr);
			complement_seq(seq_ptr,length_v);
			length_s = sr - sl - 1;
			ret = alignment_score ( &seq[sl], seq_ptr, 
						length_s, length_v, 
						2, &percent_match, &x, &y, &pad1, &pad2,
						&left1, &right1, &left2, &right2, tmode );

			complement_seq(seq_ptr,length_v);

			if ( ret < 0 ) {
			    eret =  vep_error ( fp_f, file_name, 12 );
			    goto next_file;
			}
			if ( percent_match > best_percent_match ) {
			    best_percent_match = percent_match;
			    best_match = primer_num;
			    best_match_d = primer_d;
			    best_match_pos = x + 1;
			}
			seq_ptr = vs->vector_spec_ptr[primer_num].r_primer_seq;

		    }
		}

		if ( best_percent_match >= cut_score_3 ) {

		    score_3 = best_percent_match;
		    sr = sl + best_match_pos;
		
		    /* are the primers from the same set or in the same direction? */

		    if ( best_match_5 != -9 ) {

			if ((best_match_5 != best_match) ||
			    (best_match_d5 == best_match_d)) {
			    printf("Warning: primer pair mismatch!\n");
			}
		    }
		    /*
		     *  if (!( (best_match_5 == best_match) 
		     *	&& (best_match != best_match_d) )) {
		     *  printf("Warning: primer pair mismatch!\n");
		     *}
		     *}
		     */
		    else { /* set the pr type from the 3' end match FIXME is this OK?*/
			if ( best_match_d == 0 ) pr_type = 4;
			if ( best_match_d == 1 ) pr_type = 3;
			vf_ptr = vs->vector_spec_ptr[best_match].file_name;
			sprintf(pstat, "vector_clip: %s",
				vs->vector_spec_ptr[best_match].name);
		    }
		}

		if ( best_match_5 == -9 ) {

		    /* no 5' end match */

		    /* if min_left is set we use it if no match has been found.
		     * ie it is not exactly a minimum but rather a default
		     */
		    if ( min_left > 0 ) {
			/* make sure this does not go off the end! */
			sl = MIN (min_left, sr-1);
		    }
		}

		if ( !tmode ) {
		    sprintf(buf, "%d", sl);
		    if (exp_put_str(e, EFLT_SL, buf, strlen(buf))) {
			eret =  vep_error ( fp_f, file_name, 8 );
			goto next_file;
		    }
		    sprintf(buf, "%d", sr);
		    if (exp_put_str(e, EFLT_SR, buf, strlen(buf))) {
			eret =  vep_error ( fp_f, file_name, 8 );
			goto next_file;
		    }

		    /* decided to only write PR if we found a match */

		    if ( pr_type ) {
			sprintf(buf, "%d", pr_type);
			if (exp_put_str(e, EFLT_PR, buf, strlen(buf))) {
			    eret =  vep_error ( fp_f, file_name, 8 );
			    goto next_file;
			}
		    }
		    if ( vf_ptr ) {
			if (exp_put_str(e, EFLT_SF, vf_ptr, strlen(vf_ptr))) {
			    eret =  vep_error ( fp_f, file_name, 8 );
			    goto next_file;
			}
			if (exp_put_str(e, EFLT_PS, pstat, strlen(pstat))) {
			    eret =  vep_error ( fp_f, file_name, 8 );
			    goto next_file;
			}
		    }
		    if ( sr - sl < MIN_USEFUL ) {
			eret =  vep_error ( fp_f, file_name, 14 );
			goto next_file;
		    }
		    else {
			if ( fp_p ) fprintf ( fp_p, "%s\n",file_name);
		    }
		}
		else {
		    retm = 0;
		    if ( score_5 >= cut_score_5 ) {
			printf("5' match at %d with score %f\n",sl,score_5);
			retm = 1;
		    }
		    if ( score_3 >= cut_score_3 ) {
			printf("3' match at %d with score %f\n",sr,score_3);
			retm = 1;
		    }
		    if ( retm == 0 ) {
			printf("no match\n");
		    }
		}
	    }
	    if (!(tmode)) (void) write_dot();
	}

    next_file:
	exp_destroy_info ( e );
    }
    return 0;
}



int do_hash_tr (   int seq1_len, int seq2_len,
	        int *hash_values1, int *last_word,
	        int *word_count, int *hash_values2,
		int *diag,
	        char *seq1, char *seq2,
	        int min_match, int *x, int *y, int *score ) {

    int nrw, word, pw1, pw2, i, ncw, j, match_length, word_length = 8;
    int diag_pos;
    /* seq2 is the reading which changes each entry, but seq1 is
       the vector which may be the same for consecutive
       entries, so we use new_seq1 to tell if we need to hash seq1
    */


    *score = 0;

    if ( seq1_len < min_match ) return -4; 
    if ( seq2_len < min_match ) return -4; 

    if ( hash_seq8 ( seq2, hash_values2, seq2_len )  != 0 ) {
	return -1;
    }

    j = seq1_len + seq2_len;
    for (i=0;i<j;i++) diag[i] = -word_length;

    nrw = seq2_len - word_length + 1;

/* 	loop for all (nrw) complete words in hash_values2 */


    for (pw2=0;pw2<nrw;pw2++) {

 	word = hash_values2[pw2];

	if ( -1 != word ) {

	    if ( 0 != (ncw = word_count[word]) ) {

		for (j=0,pw1=last_word[word];j<ncw;j++) {

		    diag_pos = seq1_len - pw1 + pw2 - 1;

		    if ( diag[diag_pos] < pw2 ) {

			if ((match_length = match_len ( 
						  seq1, pw1, seq1_len,
						  seq2, pw2, seq2_len))
			    >= min_match ) {
			  if ( match_length > *score ) {
			    *score = match_length;
			    *x = pw1+1;
			    *y = pw2+1;
			  }
			}
			diag[diag_pos] = pw2 + match_length;
		    }
		    pw1 = hash_values1[pw1];
		}
	    }
	}
    }
    if ( *score > 0 ) return 1;
    return 0;
}

int do_it_tr ( char *vector_seq, int max_vector, 
	       FILE *fp_i, FILE *fp_p, FILE *fp_f,
	       double cut_score_5, double cut_score_3,
	       int min_match, int tmode, int rmode, 
	       Vector_specs *vs) {

    char file_name[FILENAME_MAX+1],*cp, *seq;

    Exp_info *e;
    int ql,qr,seq_length;
    char buf[10];

/* for this algorithm */

    int x, y, ret, eret, retm;
    double score_3, score_5;
    char *seq_ptr;
    int sp=0; /* primer position */
    int sl, sr;
    int primer_length_s;
    int primer_length_v;
    int primer_d, best_match_sp;
    int pad1, pad2;
    double percent_match;
    int left1,left2,right1,right2;
    int primer_num,best_match,best_match_d, best_match_pos=0;
    double best_percent_match;
    int length_v, length_s;

    int *hash_values1, *hash_values2, *last_word, *word_count, *diag;
    int vector_length = 0;
    char vector_file_name[FILENAME_MAX+1], *vfn;
    char prev_vector_file_name[FILENAME_MAX+1];
    FILE *vf;
    int sc, prev_sc, prev_sp;
    int new_vector;
    double score_3f, score_3r, score_vf, score_vr;
    int score_f, score_r;
    int lg, rg, xf=0, xr=0, yf=0, yr=0; 
    int match_found, score;
    int size_hash = 65536;
    int word_length = 8;

    /* This routine trims transposon data, which is complicated.
       1. The transposon ends must be stored in a vector_primer file
       2. The vector sequence should be named in the SF record
       3. Numerous scores are required.
       
       4. First we get the transposon end sequences (before we get here)
       5. Get the vector sequence and rotate it around the cloning site
       6. Search with both of the transposon end sequences and note the
          highest score. If above cut_score_5 reset SL. (dynamic program)
       7. Search the TRANS_HASH bases after SL for a match to any part of the 
          vector, on both strands (hashing).
       8. If the best match is above min_match, try to align from the 
          match point to the cloning site (dynamic program). If the 
          alignment score is >= cut_score_3 reset SL.
       9. If steps 7 and 8 fail to find a match to vector we assume that
          the transposon inserted into the target DNA and not the vector.
	  The reading could hence run into vector at its 3' end so we
	  search from SL onwards, for the sequences either side of the
	  cloning site (we do not know the orientation of the transposon
          (and hence the read) relative to the vector) (dynamic program).
	  If we find a match >= cut_score_3 reset SR.
    */
    prev_sp = prev_sc = 0;

    /* initialise this algorithm */

    if ( min_match < 8 ) min_match = 8;
    prev_vector_file_name[0] = '\0';

    if ( init_hash8 ( max_vector, MAX_READ,
		     &hash_values1, &last_word, &word_count,
		     &hash_values2, &diag ))
	return -1;

    while ( fgets ( file_name, FILENAME_MAX, fp_i )) {

	if ( cp = strchr ( file_name, '\n' )) *cp = '\0';

	if ( tmode ) {
	    printf(">>>>>>>>>>>>>>>>>>>>>> %s\n", file_name );
	}

	e = exp_read_info ( file_name );

	if ( e == NULL ) {
	    eret =  vep_error ( fp_f, file_name, 1 );
	    continue;
	}
	else {
	    if ( exp_Nentries ( e, EFLT_SQ ) < 1 ) {
		eret =  vep_error ( fp_f, file_name, 2 );
		exp_destroy_info ( e );
		continue;
	      }
	    else {

	      char *expline;

	      seq = exp_get_entry ( e, EFLT_SQ );
	      seq_length = strlen ( seq );

	      if ( seq_length < MIN_USEFUL ) {
		eret =  vep_error ( fp_f, file_name, 3 );
		exp_destroy_info ( e );
		continue;
	      }
	      ql = 0;
	      qr = seq_length + 1;
	      
	      if ( rmode & USE_QL ) {
		if ( exp_Nentries ( e, EFLT_QL )) {
		  expline = exp_get_entry ( e, EFLT_QL );
		  ql = atoi ( expline );
		}
	      }
	      if ( rmode & USE_QR ) {
		if ( exp_Nentries ( e, EFLT_QR )) {
		  expline = exp_get_entry ( e, EFLT_QR );
		  qr = atoi ( expline );
		}
	      }
	      sl = ql;
	      sr = qr;
	      if ( rmode & USE_SL ) {
		if ( exp_Nentries ( e, EFLT_SL )) {
		  expline = exp_get_entry ( e, EFLT_SL );
		  sl = atoi ( expline );
		}
	      }
	      if ( rmode & USE_SR ) {
		if ( exp_Nentries ( e, EFLT_SR )) {
		  expline = exp_get_entry ( e, EFLT_SR );
		  sr = atoi ( expline );
		}
	      }
	      if ( exp_Nentries ( e, EFLT_SC )) {
		expline = exp_get_entry ( e, EFLT_SC );
		sc = atoi ( expline );
	      }
	      else {
		eret =  vep_error ( fp_f, file_name, 5 );
		exp_destroy_info ( e );
		continue;
	      }
	      if ( exp_Nentries ( e, EFLT_SP )) {
		expline = exp_get_entry ( e, EFLT_SP );
		sp = atoi ( expline );
	      }

	      new_vector = 0;
	      if (( sp != prev_sp) || ( sc != prev_sc)) new_vector = 1; 
	      prev_sp = sp;
	      prev_sc = sc;

	      /* get the vector file */

	      if ( exp_Nentries ( e, EFLT_SF )) {
		vfn = exp_get_entry(e,EFLT_SF);
		if (1 != expandpath(vfn, vector_file_name)) {
		  eret =  vep_error ( fp_f, file_name, 4 );
		  exp_destroy_info ( e );
		  continue;
		}
	      }
	      else {
		eret =  vep_error ( fp_f, file_name, 4 );
		exp_destroy_info ( e );
		continue;
	      }

		/* we only get the vector seq if it is different from
		   the previous one, or the cloning site is different,
		   or the primer site is different (actually this only
		   matters if has moved to the other strand)
		*/

	      if ( strcmp ( vector_file_name, prev_vector_file_name )) new_vector = 1;

	      if ( new_vector ) {
		vf = fopen(vector_file_name, "r");
		if (vf == NULL ) {
		  eret =  vep_error ( fp_f, file_name, 7 );
		  exp_destroy_info ( e );
		  continue;
		}

		ret = get_text_seq ( vector_seq, max_vector, &vector_length, vf);

		fclose(vf);
		if ( ret ) {
		  eret =  vep_error ( fp_f, file_name, 9 );
		  exp_destroy_info ( e );
		  continue;
		}

		/* rotate sequence so that cloning site at vector_seq[0] */

		ret = rotate_seq ( vector_seq, vector_length, sc+1);
		if ( ret ) {
		  eret =  vep_error ( fp_f, file_name, 13 );
		  exp_destroy_info ( e );
		  continue;
		}

		strcpy( prev_vector_file_name, vector_file_name );

		if ( hash_seq8 ( vector_seq, hash_values1, vector_length )  != 0 ) {
		  eret =  vep_error ( fp_f, file_name, 11 );
		  exp_destroy_info ( e );
		  continue;
		}

		(void) store_hash ( hash_values1, vector_length, last_word, 
				   word_count, word_length, size_hash);
	      
	      }

	      /* try all primer pairs */

	      score_5 = best_percent_match = -99.0;
	      for (primer_num = 0; primer_num < vs->number; primer_num++) {
		seq_ptr = vs->vector_spec_ptr[primer_num].f_primer_seq;
		for (primer_d = 0; primer_d < 2; primer_d++) {
		  sp = strlen (seq_ptr);
		  primer_length_v = abs(sp);
		  primer_length_s = MIN((abs(sp)+FUDGE),sr-sl-1);

		  ret = alignment_score ( &seq[sl], seq_ptr, 
                     primer_length_s, primer_length_v,
		     2, &percent_match, &x, &y, &pad1, &pad2,
		     &left1, &right1, &left2, &right2, tmode );

		  /*vs->vector_spec_ptr[primer_num].f_score = percent_match;*/

		  if ( percent_match > best_percent_match ) {
		    best_percent_match = percent_match;
		    best_match = primer_num;
		    best_match_d = primer_d;
		    best_match_pos = MIN(right1,right2) - left1 - pad1 + 1;
		    best_match_sp = sp;
		  }
		  seq_ptr = vs->vector_spec_ptr[primer_num].r_primer_seq;
		}
	      }

	      if ( best_percent_match >= cut_score_5 ) {

		sl += best_match_pos;
		score_5 = best_percent_match;
	      }
		/* Weve got the 5' end so now do it again looking for vector */
	      
		/* we have to search both strands so we call do_hash
		   with the read in its original sense, then its
		   complement.
		*/

	      score_3 = 0.0;
	      match_found = 0;
	      score_f = score_r = 0;
	      /*printf("ql %d qr %d sl %d sr %d\n",ql,qr,sl,sr);*/
	      lg = MAX ( ql, sl ) - 1;
	      lg = MAX ( lg, 0 );
	      rg = MIN ( qr, sr ) - 1;
	      rg = MIN ( rg, lg + TRANS_HASH );
	      /*printf("lg %d rg %d\n",lg,rg);*/
	      
	      if ( rg - lg + 1 < min_match ) {
		eret =  vep_error ( fp_f, file_name, 3 );
		exp_destroy_info ( e );
		continue;
	      }

	      ret = do_hash_tr ( vector_length, rg-lg+1,
			        hash_values1, last_word,
			        word_count, hash_values2, diag, 
			        vector_seq, &seq[lg],
			        min_match, 
			        &x, &y, &score);
	      if ( ret < 0 ) {
		eret =  vep_error ( fp_f, file_name, 11 );
		exp_destroy_info ( e );
		continue;
	      }
	      if ( ret ) {
		xf = x;
		yf = y + lg;
		score_f = score;
		match_found = 1;
		/*printf("forward score %d at %d %d\n",score_f,xf,yf);*/
	      }


	      new_vector = 0;
	      complement_seq ( &seq[lg], rg-lg+1);
	      
	      ret = do_hash_tr ( vector_length, rg-lg+1,
				hash_values1, last_word,
				  word_count, hash_values2, diag,
				vector_seq, &seq[lg],
				min_match,
				&x, &y, &score);

	      complement_seq ( &seq[lg], rg-lg+1);

	      if ( ret < 0 ) {
		eret =  vep_error ( fp_f, file_name, 11 );
		exp_destroy_info ( e );
		continue;
	      }
	      
	      if ( ret ) {
		/*printf("x %d y %d\n",x,y);*/
		xr = vector_length - x - score + 2;
		yr = lg + y + 1;
		score_r = score;
		if ( score_r > score_f )  match_found = 2;
		/*printf("reverse score %d at %d %d\n",score,xr,yr);*/
	      }

	      score_3f = score_3r = 0.0;
	      if ( match_found == 1) { 

		match_found = 0;
		length_v = vector_length - xf;
		length_s = qr - yf - 1;
		/*printf("xf %d yf %d v %d s %d\n",xf,yf,length_v,length_s);*/
		ret = alignment_score ( &seq[yf], &vector_seq[xf], length_s,
		      length_v, 2, &score_3f, &x, &y, &pad1, &pad2,
		      &left1, &right1, &left2, &right2, tmode );
		/*printf("x %d y %d left1 %d left2 %d right1 %d right2 %d\n",
		x, y,left1,left2, right1, right2);*/
		if ( ret < 0 ) {
		  eret =  vep_error ( fp_f, file_name, 12 );
		  exp_destroy_info ( e );
		  continue;
		}

		/* got caught here: there was a repeat in the read
		   near the 5' end which also matched the start of
		   the vector: the alignment routine found an excellent
		   match to the vector, but at first I did not notice
		   that it started at 410 in the read and not near
		   the 5' end. Need to check here.
		*/

		if ( (score_3f >= cut_score_3) && ( x < TRANS_HASH) ) {
		  xf = x;
		  yf = y;
		  match_found = 1;
		}
	      }
	      else if ( match_found == 2 ) {

		match_found = 0;
		length_v = vector_length - xr;
		length_s = qr - yr - 1;
		/*printf("xr %d yr %d v %d s %d\n",xr,yr,length_v,length_s);*/
		complement_seq ( vector_seq, vector_length);
		ret = alignment_score ( &seq[yr], &vector_seq[xr], length_s,
		      length_v, 2, &score_3r, &x, &y, &pad1, &pad2,
		      &left1, &right1, &left2, &right2, tmode );
		complement_seq ( vector_seq, vector_length);
		/*printf("x %d y %d left1 %d left2 %d right1 %d right2 %d\n",
		x, y,left1,left2, right1, right2);*/

		if ( ret < 0 ) {
		  eret =  vep_error ( fp_f, file_name, 12 );
		  exp_destroy_info ( e );
		  continue;
		}
		if ( (score_3r >= cut_score_3) && ( x < TRANS_HASH) ) {
		  xr = x;
		  yr = y;
		  match_found = 2;
		}
	      }

	      if ( match_found ) {

		if ( score_3f > score_3r ) {
		  sl += yf;
		  score_5 = score_3f;
		}
		else {
		  sl += yr;
		  score_5 = score_3r;
		}
	      }
	      else { /* look for running into vector */

		/* do reverse cloning first */

		length_v = MIN(CLONESITE_SIZE,vector_length);
		length_s = sr-sl-1;

		ret = alignment_score ( &seq[sl], 
		     vector_seq,/*[vector_length-length_v], */
		     length_s,
		     length_v,
		     2,
		     &score_vr,
		     &xr,
		     &y,
		     &pad1,
                     &pad2,
		     &left1,
		     &right1,
		     &left2,
		     &right2,
		     tmode
		     );
		/*
		printf("xr %d score_vr %f\n",xr,score_vr);
		printf("left1 %d right1 %d left2 %d right2 %d\n",
		       left1,right1,left2,right2);
		*/

		if ( ret < 0 ) {
		  eret =  vep_error ( fp_f, file_name, 12 );
		  exp_destroy_info ( e );
		  continue;
		}

		/*printf("sl %d sr %d length_v %d length_s %d\n",
		       sl,sr,length_v,length_s);*/

		complement_seq ( vector_seq, vector_length);
		ret = alignment_score ( &seq[sl], vector_seq, length_s,
		      length_v, 2, &score_vf, &xf, &y, &pad1, &pad2,
		      &left1, &right1, &left2, &right2, tmode );

		complement_seq ( vector_seq, vector_length);

		/*printf("xf %d score_vf %f\n",xf,score_vf);
		printf("left1 %d right1 %d left2 %d right2 %d\n",
		       left1,right1,left2,right2);*/

		if ( ret < 0 ) {
		  eret =  vep_error ( fp_f, file_name, 12 );
		  exp_destroy_info ( e );
		  continue;
		}


		if ( score_vf > score_vr ) {
		  if ( score_vf >= cut_score_3 ) {
		    sr = sl + xf + 1;
		    score_3 = score_vf;
		  }
		}
		else if (score_vr >= cut_score_3) {
		  sr = sl + xr + 1;
		  score_3 = score_vr;
		}

	      }


	      if ( !tmode ) {
		sprintf(buf, "%d", sl);
		if (exp_put_str(e, EFLT_SL, buf, strlen(buf))) {
		  eret =  vep_error ( fp_f, file_name, 8 );
		  exp_destroy_info ( e );
		  continue;
		}
		sprintf(buf, "%d", sr);
		if (exp_put_str(e, EFLT_SR, buf, strlen(buf))) {
		  eret =  vep_error ( fp_f, file_name, 8 );
		  exp_destroy_info ( e );
		  continue;
		}
		if ( sr - sl < MIN_USEFUL ) {
		  eret =  vep_error ( fp_f, file_name, 14 );
		  exp_destroy_info ( e );
		  continue;
		}
		else {
		  if ( fp_p ) fprintf ( fp_p, "%s\n",file_name);
		}
	      }
	      else {
		retm = 0;
		if ( score_5 >= cut_score_5 ) {
		  printf("5' match at %d with score %f\n",sl,score_5);
		  retm = 1;
		}
		if ( score_3 >= cut_score_3 ) {
		  printf("3' match at %d with score %f\n",sr,score_3);
		  retm = 1;
		}
		if ( retm == 0 ) {
		  printf("no match\n");
		}
	      }
	    }
	    exp_destroy_info ( e );
	    if (!(tmode)) (void) write_dot();
	  }
      }
    
    return 0;
  }


int do_it_cv ( char *vector_seq, int max_vector,
	       FILE *fp_i, FILE *fp_p, FILE *fp_f,
	       int word_length, int num_diags, double diag_score, double max_prob,
	       int tmode) {

    char file_name[FILENAME_MAX+1],*cp, *seq;

    Exp_info *e;
    int ql,qr,seq_length;
    char buf[10];

/* for this algorithm */

    int *hash_values1, *hash_values2, *last_word, *word_count, *diag, *line;
    int size_hash;
    int vector_length = 0, x, y, ret, eret;
    DI *hist;
    char vector_file_name[FILENAME_MAX+1], *vfn;
    char prev_vector_file_name[FILENAME_MAX+1];
    FILE *vf;
    int sl, sr; /* sequencing vector left and right */
    double score_3f, score_f, score_3r, score_r;
    int lg, rg, xf=0, xr=0, cl=0, cr=0;
    int match_found;
    double *expected_scores;

#define end_fudge 4

    /* initialise this algorithm */

    if ( ! ( expected_scores = (double *) xmalloc ( sizeof(double)*
						   (max_vector+MAX_READ)/2)))
      return -2;
    if ( max_prob > 0.0 ) {
      if ( poisson_diagonals ( MINMAT_CV, MAX_READ, word_length, max_prob, 
			    expected_scores )) return -1;
    }

    prev_vector_file_name[0] = '\0';

    if ( init_hash ( word_length, max_vector, MAX_READ,
	     &hash_values1, &last_word, &word_count,
		     &hash_values2, &diag, &hist, &size_hash, &line ))
	return -1;

    while ( fgets ( file_name, FILENAME_MAX, fp_i )) {

	if ( cp = strchr ( file_name, '\n' )) *cp = '\0';

	if ( tmode ) {
	    printf(">>>>>>>>>>>>>>>>>>>>>> %s\n", file_name );
	}

	e = exp_read_info ( file_name );
	if ( e == NULL ) {
	    eret =  vep_error ( fp_f, file_name, 1 );
	}
	else {
	    if ( exp_Nentries ( e, EFLT_SQ ) < 1 ) {
		eret =  vep_error ( fp_f, file_name, 2 );
		exp_destroy_info ( e );
		continue;
	    }
	    else {

		char *expline;

		seq = exp_get_entry ( e, EFLT_SQ );
		seq_length = strlen ( seq );
		ql = 0;
		qr = seq_length + 1;

		if ( exp_Nentries ( e, EFLT_QL )) {
		    expline = exp_get_entry ( e, EFLT_QL );
		    ql = atoi ( expline );
		}
		if ( exp_Nentries ( e, EFLT_QR )) {
		    expline = exp_get_entry ( e, EFLT_QR );
		    qr = atoi ( expline );
		}

		sl = 0;
		sr = seq_length + 1;
		if ( exp_Nentries ( e, EFLT_SL )) {
		    expline = exp_get_entry ( e, EFLT_SL );
		    sl = atoi ( expline );
		}
		if ( exp_Nentries ( e, EFLT_SR )) {
		    expline = exp_get_entry ( e, EFLT_SR );
		    sr = atoi ( expline );
		}


		/* get the vector file */

		if ( exp_Nentries ( e, EFLT_CF )) {
		    vfn = exp_get_entry(e,EFLT_CF);
		    if (1 != expandpath(vfn, vector_file_name)) {
		      eret =  vep_error ( fp_f, file_name, 4 );
		      exp_destroy_info ( e );
		      continue;
		    }
		}
		else {
		    eret =  vep_error ( fp_f, file_name, 4 );
		    exp_destroy_info ( e );
		    continue;
		}


		/* we only get the vector seq if it is different from
		   the previous one.
		*/

		if ( strcmp ( vector_file_name, prev_vector_file_name )) {

		    vf = fopen(vector_file_name, "r");
		    if (vf == NULL ) {
			eret =  vep_error ( fp_f, file_name, 7 );
			exp_destroy_info ( e );
			continue;
		    }

		    ret = get_text_seq ( vector_seq, max_vector, &vector_length, vf);
		    fclose(vf);
		    if ( ret ) {
			eret =  vep_error ( fp_f, file_name, 9 );
			exp_destroy_info ( e );
			continue;
		    }
		    strcpy( prev_vector_file_name, vector_file_name );


		    if ( hash_seq ( word_length, vector_seq, hash_values1, 
				   vector_length )  != 0 ) {
		      eret =  vep_error ( fp_f, file_name, 11 );
		      exp_destroy_info ( e );
		      continue;
		    }

		    (void) store_hash ( hash_values1, vector_length, last_word,
				       word_count, word_length, size_hash);
		}

		/* we have to search both strands so we call do_hash
		   with the read in its original sense, then its
		   complement. To be safe we get the best score
		   from these two and then check against min
		*/


		/* sequence numbering is from 1 but arrays go from 0
		   and right positions can be 1 past the right hand end
		*/

		match_found = 0;
		/* to avoid the poor quality data and the sequencing vector */
		/*
		lg = MAX ( ql, sl ) - 1;
		lg = MAX ( lg, 0 );
		rg = MIN ( qr, sr ) - 1;
		rg = MIN ( rg, seq_length-1 );
		*/

		/* to avoid the sequencing vector but search within poor quality */
		lg = MAX ( 0, sl-1 );
		rg = MIN ( sr-1, seq_length-1 );
/*
		printf("ql %d sl %d qr %d sr %d lg %d rg %d\n",ql,sl,qr,sr,lg,rg);
		printf("vector_length %d rg-lg+1 %d\n",vector_length,rg-lg+1);
*/

		/* only search if at least MIN_USEFUL chars of data */

		if ( rg-lg+1 < MIN_USEFUL ) {
		    eret =  vep_error ( fp_f, file_name, 3 );
		    exp_destroy_info ( e );
		    continue;
		}
		ret = do_hash ( vector_length, rg-lg+1,
			        hash_values1, last_word,
			        word_count, hash_values2,
			        diag, line, size_hash, hist, word_length,
			        vector_seq, &seq[lg],
			        diag_score, num_diags, expected_scores, max_prob,
			        CLONING_VECTOR,
			        &x, &y, &score_3f);
		if ( ret < 0 ) {
		    eret =  vep_error ( fp_f, file_name, 11 );
		    exp_destroy_info ( e );
		    continue;
		}
		score_f = 0.0;
		if ( ret ) {
		    score_f = score_3f;
		    xf = x;
		    match_found = 1;
		}

	        complement_seq ( &seq[lg], rg-lg+1);

		ret = do_hash ( vector_length, rg-lg+1,
			        hash_values1, last_word,
			        word_count, hash_values2,
			        diag, line, size_hash, hist, word_length,
			        vector_seq, &seq[lg],
			        diag_score, num_diags, expected_scores, max_prob,
			        CLONING_VECTOR,
			        &x, &y, &score_3r);
		if ( ret < 0 ) {
		    eret =  vep_error ( fp_f, file_name, 11 );
		    exp_destroy_info ( e );
		    continue;
		}
		score_r = 0.0;
		if ( ret ) {
		    score_r = score_3r;
		    match_found = 1;
		    xr = x;
		}

		if ( match_found ) { 
		  
		    if ( ( score_f > score_r ) && ( score_f > 0.0 ) ) {
			cl = MAX ( (lg+xf-vector_length+1),sl+1);
			cr = MIN ( lg+xf,sr-1);
			/*printf("f lg %d xf %d sl %d\n",lg,xf,sl);*/
		    }
		    else if ( ( score_r > score_f ) && ( score_r > 0.0 ) ) {

			cl = MAX ( rg-xr,sl+1);
			cr = MIN ( rg-(xr-vector_length),sr-1);
			/*printf("r rg %d xr %d sl %d\n",rg,xr,sl);*/
		    }

		    if ( !tmode ) {

			sprintf(buf, "%d", cl);
			if (exp_put_str(e, EFLT_CL, buf, strlen(buf))) {
			    eret =  vep_error ( fp_f, file_name, 8 );
			    exp_destroy_info ( e );
			    continue;
			}
			sprintf(buf, "%d", cr);
			if (exp_put_str(e, EFLT_CR, buf, strlen(buf))) {
			    eret =  vep_error ( fp_f, file_name, 8 );
			    exp_destroy_info ( e );
			    continue;
			}
			if (exp_put_rng(e, EFLT_CS, &cl, &cr)) {
			    eret =  vep_error ( fp_f, file_name, 8 );
			    exp_destroy_info ( e );
			    continue;
			}
			/* check if all cloning vector */

			  /*printf("sl %d sr %d\n",sl,sr);
			  printf("cl %d cr %d\n",cl,cr);*/
			/*if (( cl == sl+1 ) && ( cr == sr-1 )) {*/
			if (abs(cl-(sl+1))<end_fudge
			    && abs(cr-(sr-1))<end_fudge) {
			    eret =  vep_error ( fp_f, file_name, 15 );
			    if (exp_put_str(e, EFLT_PS, "all cloning vector", 
					    strlen("all cloning vector"))) {
				eret =  vep_error ( fp_f, file_name, 8 );
			    }
			    exp_destroy_info ( e );
			    continue;
			}
			else {
			    if ( fp_p ) fprintf ( fp_p, "%s\n",file_name);
			}
		    }
		    else {
			printf(" mark %d to %d. score %f\n",cl,cr,MAX(score_r,score_f));
		    }

		}
		else {
		    if ( tmode ) {
			printf("no match\n");
		    }
		    else {
			if ( fp_p ) fprintf ( fp_p, "%s\n",file_name);
		    }
		}
	    }
	}
	exp_destroy_info ( e ); 
	if (!(tmode)) (void) write_dot();
    }
    free_hash ( hash_values1, last_word,
	       word_count, hash_values2,
	       diag, hist );

    xfree(expected_scores);

    return 0;
}


int dna_hash8_lookup[256];

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

int init_hash8 ( int seq1_len, int seq2_len,
	        int **hash_values1, int **last_word,
	        int **word_count, int **hash_values2,
	        int **diag ) {
    int size_hash, word_length;

    word_length = 8;

    set_hash8_lookup ();

    size_hash = 65536;


    if ( NULL == (*hash_values1 = (int *) xmalloc ( sizeof(int)*(seq1_len) ))) {
	return -2;
    }

    if ( ! (*last_word = (int *) xmalloc ( sizeof(int)*size_hash ))) {
	return -2;
    }

    if ( ! (*word_count = (int *) xmalloc ( sizeof(int)*size_hash ))) {
	return -2;
    }


    if ( ! (*hash_values2 = (int *) xmalloc ( sizeof(int)*(seq2_len) ))) {
	return -2;
    }

    if ( ! (*diag = (int *) xmalloc ( sizeof(int)*(seq2_len + seq1_len) ))) {
	return -2;
    }

    return 0;
}

void free_hash8 ( int *hash_values1, int *last_word,
	        int *word_count, int *hash_values2,
	        int *diag) {

    if ( hash_values1 ) xfree ( hash_values1 );
    if ( hash_values2 ) xfree ( hash_values2 );
    if ( word_count )   xfree ( word_count );
    if ( last_word )    xfree ( last_word );
    if ( diag )         xfree ( diag );

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

    for (i=0;i<seq_len;i++) hash_values[i] = -1;

    /*	Now do the rest of the sequence */

    hash_values[start_base] = uword;
    k = seq_len - word_len + 1;

    for (i=start_base+1,j=start_base+word_len; i<k; i++,j++) {

	base_index = dna_hash8_lookup[(unsigned)seq[j]];
	if ( 4 == base_index ) {

	    /*	weve hit an unknown char, so lets start again */

	    prev_start_base = i;
	    start_base = j + 1;
	    if (hash_word8 ( seq, &start_base, seq_len, &uword)) return 0;

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
int hash_seq8_prev ( char *seq, int *hash_values, int seq_len) {
  /* changed to new version 10.10.02 */
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

int do_hash_vr (   int seq1_len, int seq2_len,
	        int *hash_values1, int *last_word,
	        int *word_count, int *hash_values2,
		int *diag,
	        char *seq1, char *seq2, 
	        int min_match, int *x, int *y, int *score ) {

    int nrw, word, pw1, pw2, i, ncw, j, match_length, word_length = 8;
    int diag_pos;

    *score = 0;
    if ( seq1_len < min_match ) return -4; 
    if ( seq2_len < min_match ) return -4; 

    if ( hash_seq8 ( seq2, hash_values2, seq2_len )  != 0 ) {
	return -1;
    }

    j = seq1_len + seq2_len;
    for (i=0;i<j;i++) diag[i] = -word_length;

    nrw = seq2_len - word_length + 1;

/* 	loop for all (nrw) complete words in hash_values2 */


    for (pw2=0;pw2<nrw;pw2++) {

 	word = hash_values2[pw2];

	if ( -1 != word ) {

	    if ( 0 != (ncw = word_count[word]) ) {

		for (j=0,pw1=last_word[word];j<ncw;j++) {

		    diag_pos = seq1_len - pw1 + pw2 - 1;

		    if ( diag[diag_pos] < pw2 ) {

			if ((match_length = match_len ( 
						  seq1, pw1, seq1_len,
						  seq2, pw2, seq2_len))
			    >= min_match ) {

			    *score = match_length;
			    *x = pw1+1;
			    *y = pw2+1;
			    return 1;
			}
			diag[diag_pos] = pw2 + match_length;
		    }
		    pw1 = hash_values1[pw1];
		}
	    }
	}
    }
    return 0;
}

int do_it_vr ( char *vector_seq, int max_vector,
	       FILE *fp_i, FILE *fp_p, FILE *fp_f,
	       int word_length, int min_match,
	       int tmode) {

    char file_name[FILENAME_MAX+1],*cp, *seq;

    Exp_info *e;
    int ql,qr,seq_length;

/* for this algorithm */

    int *hash_values1, *hash_values2, *last_word, *word_count;
    int *diag;
    int vector_length = 0, x, y, ret, eret;
    char vector_file_name[FILENAME_MAX+1], *vfn;
    char prev_vector_file_name[FILENAME_MAX+1];
    FILE *vf;
    int sl, sr; /* sequencing vector left and right */
    int score, score_f, score_r;
    int lg, rg, xf=0, xr=0, yf=0, yr=0;
    int match_found;
    int size_hash;

    size_hash = 65536;

    /* initialise this algorithm */


    if ( min_match < 8 ) min_match = 8;
    prev_vector_file_name[0] = '\0';

    if ( init_hash8 ( max_vector, MAX_READ,
		     &hash_values1, &last_word, &word_count,
		     &hash_values2, &diag ))
	return -1;


    while ( fgets ( file_name, FILENAME_MAX, fp_i )) {

	if ( cp = strchr ( file_name, '\n' )) *cp = '\0';

	if ( tmode ) {
	    printf(">>>>>>>>>>>>>>>>>>>>>> %s\n", file_name );
	}

	e = exp_read_info ( file_name );
	if ( e == NULL ) {
	    eret =  vep_error ( fp_f, file_name, 1 );
	    exp_destroy_info ( e );
	    continue;
	}
	else {
	    if ( exp_Nentries ( e, EFLT_SQ ) < 1 ) {
		eret =  vep_error ( fp_f, file_name, 2 );
		exp_destroy_info ( e );
		continue;
	    }
	    else {

		char *expline;

		seq = exp_get_entry ( e, EFLT_SQ );
		seq_length = strlen ( seq );
		ql = 0;
		qr = seq_length + 1;

		if ( exp_Nentries ( e, EFLT_QL )) {
		    expline = exp_get_entry ( e, EFLT_QL );
		    ql = atoi ( expline );
		}
		if ( exp_Nentries ( e, EFLT_QR )) {
		    expline = exp_get_entry ( e, EFLT_QR );
		    qr = atoi ( expline );
		}

		sl = 0;
		sr = seq_length + 1;
		if ( exp_Nentries ( e, EFLT_SL )) {
		    expline = exp_get_entry ( e, EFLT_SL );
		    sl = atoi ( expline );
		}
		if ( exp_Nentries ( e, EFLT_SR )) {
		    expline = exp_get_entry ( e, EFLT_SR );
		    sr = atoi ( expline );
		}

		/* get the vector file */

		if ( exp_Nentries ( e, EFLT_SF )) {
		    vfn = exp_get_entry(e,EFLT_SF);
		    if (1 != expandpath(vfn, vector_file_name)) {
		      eret =  vep_error ( fp_f, file_name, 4 );
		      exp_destroy_info ( e );
		      continue;
		    }
		}
		else {
		  /* if the Exp file doe snot have the SF record pass it! */
		  if ( fp_p ) fprintf ( fp_p, "%s\n",file_name);
		  if ( !tmode ) {

		    if (exp_put_str(e, EFLT_CC, 
				    "no SF for vector rearrangement", 
			     strlen("no SF for vector rearrangement"))) {
		      eret =  vep_error ( fp_f, file_name, 16 );
		      exp_destroy_info ( e );
		      continue;
		    }
		  exp_destroy_info ( e );
		  continue;
		  }
		}

		/* we only get the vector seq if it is different from
		   the previous one.
		*/

		if ( strcmp ( vector_file_name, prev_vector_file_name )) {
		    vf = fopen(vector_file_name, "r");
		    if (vf == NULL ) {
			eret =  vep_error ( fp_f, file_name, 7 );
			exp_destroy_info ( e );
			continue;
		    }

		    ret = get_text_seq ( vector_seq, max_vector, &vector_length, vf);
		    fclose(vf);
/* 		    printf("vector length %d\n",vector_length); */
		    if ( ret ) {
			eret =  vep_error ( fp_f, file_name, 9 );
			exp_destroy_info ( e );
			continue;
		    }
		    strcpy( prev_vector_file_name, vector_file_name );

		    if ( hash_seq8 ( vector_seq, hash_values1, vector_length )  != 0 ) {
		      eret =  vep_error ( fp_f, file_name, 11 );
		      exp_destroy_info ( e );
		      continue;
		    }

		    (void) store_hash ( hash_values1, vector_length, last_word, 
				       word_count,
			    word_length, size_hash);
		  }

		/* we have to search both strands so we call do_hash
		   with the read in its original sense, then its
		   complement.
		*/

		match_found = 0;
		score_f = score_r = 0;
		lg = MAX ( ql, sl ) - 1;
		lg = MAX ( lg, 0 );
		rg = MIN ( qr, sr ) - 1;
		rg = MIN ( rg, seq_length-1 );
/* 		printf("lg %d rg %d\n",lg,rg); */

		if ( rg - lg + 1 < min_match ) {
		    /* Don't fail it - put it in the passed file */
		    if ( fp_p ) fprintf ( fp_p, "%s\n",file_name);
		    exp_destroy_info ( e );
		    fprintf ( stdout, "-" );
		    (void) fflush ( stdout );
		    continue;
		}
		ret = do_hash_vr ( vector_length, rg-lg+1,
			        hash_values1, last_word,
			        word_count, hash_values2, diag, 
			        vector_seq, &seq[lg],
			        min_match, 
			        &x, &y, &score);
		if ( ret < 0 ) {
		    eret =  vep_error ( fp_f, file_name, 11 );
		    exp_destroy_info ( e );
		    continue;
		}
		if ( ret ) {
		    xf = x + 1;
		    yf = y + lg + 1;
		    score_f = score;
		    match_found = 1;
/* 		    printf("forward score %d at %d %d\n",score_f,xf,yf); */
		}

		if ( !match_found ) {

		    complement_seq ( &seq[lg], rg-lg+1);

		    ret = do_hash_vr ( vector_length, rg-lg+1,
			        hash_values1, last_word,
			        word_count, hash_values2, diag,
			        vector_seq, &seq[lg],
			        min_match,
			        &x, &y, &score);
		    if ( ret < 0 ) {
			eret =  vep_error ( fp_f, file_name, 11 );
			exp_destroy_info ( e );
			continue;
		    }
		    if ( ret ) {
			xr = x + 1;
			yr = rg - lg - y + lg - score + 2;
			score_r = score;
			match_found = 1;
/*		    	printf("reverse score %d at %d %d\n",score,xr,yr); */
		    }
		}

		if ( match_found ) { 


		    /* this is left over from when I wanted to find the best match
		       but now do_hash_vr exits as soon as a long enough match is 
		       found, so this is really redundant
		       */


		    if (score_f > score_r ) {
			x = xf;
			y = yf;
			score = score_f;
		    }
		    else {
			x = xr;
			y = yr;
			score = score_r;
		    }
		    if ( !tmode ) {

			if (exp_put_str(e, EFLT_PS, "vector rearrangement", 
			strlen("vector rearrangement"))) {
			    eret =  vep_error ( fp_f, file_name, 16 );
			    exp_destroy_info ( e );
			    continue;
			}

			if ( fp_f ) fprintf ( fp_f, "%s\n",file_name);
		    }
		    else {
			printf("match length %d at %d  %d\n",score,x,y);
		    }

		}
		else {
		    if ( tmode ) {
			printf("no match\n");
		    }
		    else {
			if ( fp_p ) fprintf ( fp_p, "%s\n",file_name);
		    }
		}
	    }
	}
	exp_destroy_info ( e );
	if (!(tmode)) (void) write_dot();
    }
    free_hash8 ( hash_values1, last_word,
	         word_count, hash_values2,
	         diag );
    return 0;
}

int get_text_seq ( char *seq, int max_len, int *seq_len, FILE *fp )

/* read in a staden format (yuk) sequence file */

/* Deal with 2 special line types: comments that have ";" in column 0
   and contig consensus sequence headers that have "<----abc.00001---->"
   embedded in them */
{

#define MAX_SEQ_LINE 101

    char line[MAX_SEQ_LINE];
    int j;

    *seq_len = 0;
    while ( fgets( line,sizeof(line),fp ) != NULL ) {

	/* Check for special lines of type ";"*/

	if ( ';' != line[0] ) {

	    for (j = 0;j < MAX_SEQ_LINE && line[j]; j++) {

		if ( '<' == line[j] ) j += 20;
		if (isalpha ( (int) line[j]) || (int) line[j] == '-') {
		    if ( *seq_len >= max_len) return -1;
		     seq[*seq_len] = line[j];
		    *seq_len += 1;
		}
	    }
	}
    }
    return 0;
}


void set_hash_consts (int word_length) {
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


int hash_seq ( int word_length, char *seq, int *hash_values, int seq_len ) {

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

int diagonal_length ( int seq1_len, int seq2_len, int diagonal_number ) {

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

  if ( diagonal_number < seq1_len ) {
    return MIN ( (diagonal_number + 1), (MIN (seq1_len, seq2_len) ));
  }
  else {
    return MIN ( (seq1_len + seq2_len - 1 - diagonal_number), (MIN (seq1_len, seq2_len) ));
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
/*    printf("x %d y %d\n",*x,*y); */
}


int init_hash ( int word_length, int seq1_len, int seq2_len,
	        int **hash_values1, int **last_word,
	        int **word_count, int **hash_values2, int **diag,
	        DI **hist, int *size_hash, int **line ) {
    int size_hist;

    set_hash_consts (word_length);

    *size_hash = pow(char_set_size-1, word_length);


    if ( NULL == (*hash_values1 = (int *) xmalloc ( sizeof(int)*(seq1_len) ))) {
	return -2;
    }

    if ( ! (*last_word = (int *) xmalloc ( sizeof(int)*(*size_hash) ))) {
	return -2;
    }

    if ( ! (*word_count = (int *) xmalloc ( sizeof(int)*(*size_hash) ))) {
	return -2;
    }


    if ( ! (*hash_values2 = (int *) xmalloc ( sizeof(int)*(seq2_len) ))) {
	return -2;
    }

    size_hist = seq1_len + seq2_len - 1;
    if ( ! (*hist = (DI *) xmalloc ( sizeof(DI)*size_hist ))) {
	return -2;
    }

    if ( ! (*line = (int *) xmalloc ( sizeof(int)*(seq2_len) ))) {
	return -2;
    }

    if ( ! (*diag = (int *) xmalloc ( sizeof(int)*(seq2_len + seq1_len) ))) {
	return -2;
    }
    return 0;
}

int do_hash (   int seq1_len, int seq2_len,
	        int *hash_values1, int *last_word,
	        int *word_count, int *hash_values2, int *diag, int *line,
	        int size_hash, DI *hist, int word_length,
	        char *seq1, char *seq2,
	        double min_score, int num_diags, double *expected_scores,
	        double max_prob, int job,
	        int *x, int *y, double *score_3) {

    /* seq2 is the reading which changes each entry, but seq1 is
       the primer or the vector which may be the same for consecutive
       entries, so we use new_seq1 to tell if we need to hash seq1

       We return < 0 for errors, 0 for no match, 1 for a match

       At present the algorithms are unsatisfactory: if max_prob is set
       we find the least probable match while accumulating all matches
       that are less probable than max_prob without using the accumulated
       list! otherwise we use the old algorithm.

    */

    int size_hist;
    register int i,j,nrw, ncw, word, pw1, pw2, diag_pos, match_length;
    double points, max_diagonal;
    int score_pos, top_score_pos=0;
    double top_score;

    int max_line, line_pos, kk, half_window;
    int min_s2, max_s2, min_s1, max_s1;
    double score;
    int match_no;
    int hist_left, hist_right;
    int search_mode;
    int probs;

    typedef struct match_ {
	int pos;
	double prob;
    } Match;

    Match match[MAX_MATCHES];  
    int matches;

    probs = 1;
    search_mode = 0;
    if ( max_prob > 0.0 ) search_mode = probs;
    /* quick sanity checks */

    if ( seq1_len < word_length ) return -4; 
    if ( seq2_len < word_length ) return -4; 

    if ( hash_seq ( word_length, seq2, hash_values2, seq2_len )  != 0 ) {
	return -2;
    }


/*	now we have a sequence seq1 encoded by hash tables:
	word_count[i] contains the number of occurrences of hash value i
	last_word[i] contains the seq1 position of the last occurrence of hash value i
	hash_values1[i] contains the hash values of previous occurrence of hash value i
	
	we have seq2 encoded in an array of hash values hash_values2[i], each i
	corresponding to the string at seq2[i]

	algorithm: for each hash value in hash_values2[i] test if it exists in word_count
	if so process all occurrences using last-word and hash_values1.

	For some reason we are careful with memory and use a union for the histogram.
	seq1 on x axis. seq2 on y.

*/


    j = seq1_len + seq2_len;
    for(i=0; i<j; i++) diag[i] = -word_length;

    nrw = seq2_len - word_length + 1;
    size_hist = seq1_len + seq2_len - 1;

    for(i=0;i<size_hist;i++)(hist[i].i = 0);

    /* 	loop for all (nrw) complete words in hash_values2 */

    for (pw2=0;pw2<nrw;pw2++) {
      word = hash_values2[pw2];
      if ( -1 != word ) {
	if ( 0 != (ncw = word_count[word]) ) {
	  for (j=0,pw1=last_word[word];j<ncw;j++) {
	    diag_pos = seq1_len - pw1 + pw2 - 1;
	    if (diag[diag_pos]<pw2) {
	      match_length = match_len (seq1,pw1,seq1_len,
					seq2,pw2,seq2_len);
	      diag[diag_pos] = pw2 + match_length;
	      hist[diag_pos].i += match_length;
	    }
	    pw1 = hash_values1[pw1];
	  }
	}
      }
    }

    for(i=0;i<word_length-1;i++) {
      hist[i].d = 0.0;
    }
    for(i=size_hist-1;i>size_hist-word_length;i--) {
      hist[i].d = 0.0;
    }

    hist_left = MINMAT_CV;
    hist_right = size_hist - MINMAT_CV;

    half_window = num_diags/2;
    max_diagonal = MIN(seq1_len,seq2_len);
    matches = -1;
    *score_3 = 0.0;

    if ( search_mode == probs ) {

      /* now do the correction for the length of each diagonal
       * and find the best scores
       */

      max_diagonal = MIN(seq1_len,seq2_len);

      for(i=hist_left,j=diagonal_length(seq1_len,seq2_len,i);i<seq1_len;i++,j++) {
	points = j;
	hist[i].d = hist[i].i / points;
	j += 1;
	j = MIN (max_diagonal, (j));
	if ( hist[i].d > expected_scores[j] ) { 
	  matches++;
	  if ( MAX_MATCHES == matches ) return -5;
	/*printf("pos %d score %f exp %f\n",i,hist[i].d,expected_scores[j]);*/
	  match[matches].pos = i;
	  match[matches].prob = 100.0*(hist[i].d - expected_scores[j])/expected_scores[j];
	}
      }
      
      for(i=hist_right,j=diagonal_length(seq1_len,seq2_len,i);i>=seq1_len;i--,j++) {
	points = j;
	hist[i].d = hist[i].i / points;
	j += 1;
	j = MIN (max_diagonal, (j));
	if ( hist[i].d > expected_scores[j] ) { 
	  matches++;
	  if ( MAX_MATCHES == matches ) return -5;
	/*printf("pos %d score %f exp %f\n",i,hist[i].d,expected_scores[j]);*/
	  match[matches].pos = i;
	  match[matches].prob = 100.0*(hist[i].d - expected_scores[j])/expected_scores[j];
	}
      }
      matches += 1;
      if (matches < 1) return 0;
      
      top_score = 0.;
      for(match_no=0;match_no<matches;match_no++) {
	if ( match[match_no].prob > top_score ) {
	  top_score_pos = match[match_no].pos;
	  top_score = match[match_no].prob;
	}
      }
      *score_3 = top_score;
      *x = top_score_pos;
      return 1;
    }

    else {

      /* this is similar to the old algorithm but contains corrections*/

      for(i=hist_left,j=diagonal_length(seq1_len,seq2_len,i);i<seq1_len;i++,j++) {
	points = j;
	hist[i].d = hist[i].i / points;
	j += 1;
	j = MIN (max_diagonal, (j));
	if ( hist[i].d > min_score ) { 
	  matches++;
	  if ( MAX_MATCHES == matches ) return -5;
	  match[matches].pos = i;
	  match[matches].prob = hist[i].d;
	}
      }
      
      for(i=hist_right,j=diagonal_length(seq1_len,seq2_len,i);i>=seq1_len;i--,j++) {
	points = j;
	hist[i].d = hist[i].i / points;
	j += 1;
	j = MIN (max_diagonal, (j));
	if ( hist[i].d > min_score ) { 
	  matches++;
	  if ( MAX_MATCHES == matches ) return -5;
	  match[matches].pos = i;
	  match[matches].prob = hist[i].d;
	}
      }
      matches += 1;
      if (matches < 1) return 0;
      
      top_score = 0.;
      for(match_no=0;match_no<matches;match_no++) {
	score_pos = match[match_no].pos;

	if ( score_pos < seq1_len ) {
	  min_s2 = 0;
	  max_s2 = MIN(seq2_len,score_pos+1);
	}
	else {
	  min_s2 = score_pos - seq1_len + 1;
	  max_s2 = seq2_len;
	}

	max_line = MAX(1,max_s2-min_s2);
	min_s2 = MAX(0,min_s2-half_window);
	max_s2 = MIN(seq2_len,max_s2+half_window);
	min_s1 = MAX(0,score_pos-half_window);
	max_s1 = MIN(score_pos+half_window,seq1_len+seq2_len-1);

	for ( i=0;i<=max_line;i++) line[i] = 0;
	for (pw2=min_s2;pw2<=max_s2;pw2++) {
	  word = hash_values2[pw2];
	  if ( -1 != word ) {
	    if ( 0 != (ncw = word_count[word]) ) {
	      pw1 = last_word[word];
	      for (j=0;j<ncw;j++) {
		for(kk=0;kk<word_length;kk++) {
		  diag_pos = seq1_len - pw1 + pw2 - 1 + kk;
		  if((diag_pos>=min_s1)&&(diag_pos<=max_s1)) {
		    line_pos = pw2 - min_s2 + kk;
		    if ( (line_pos >= 0) && (line_pos <= max_line) ) line[line_pos]=1;
		  }
		}
		pw1 = hash_values1[pw1];
	      }
	    }
	  }
	}

	for (i=0,score=0.0;i<=max_line;i++) score += line[i];

	score = score/(double)max_line;

	if (score>top_score) {
	  top_score = score;
	  top_score_pos = score_pos;
	}
      }
      if ( top_score > min_score ) {
	*score_3 = top_score;
	*x = top_score_pos;
	return 1;
      }

    }
    return 0;
  }



extern DLL_IMPORT char *optarg;
extern DLL_IMPORT int optind;

void usage(int word_length, int num_diags, double diag_score, int min_match, 
	   int cut_score_5, int cut_score_3) {
    fprintf(stderr, 
	    "Usage: vector_clip [options] file_of_filenames\n"
	    "Where options are:\n"
	    "    [-s mark sequencing vector]      [-c mark cloning vector]\n"
	    "    [-h hgmp primer]                 [-r vector rearrangements]\n"
	    "    [-w word_length (%d)]             [-n num_diags (%d)]\n"
	    "    [-d diagonal score (%.2f)]       [-l minimum match (%d)]\n"
	    "    [-L minimum %% 5' match (%d)]     [-R minimum %% 3' match (%d)]\n"
	    "    [-m default 5' position]         [-t test only]\n"
	    "    [-M Max vector length (%d)]  [-P max Probability]\n"
	    "    [-v vector_primer filename]      [-i vector_primer filename]\n"
	    "    [-V vector_primer length]\n"

	    "    [-p passed fofn]                 [-f failed fofn]\n",
	     word_length, num_diags, diag_score, min_match, cut_score_5, cut_score_3,
	    MAX_VECTOR_D);
    exit(1);
}

/* program to locate and mark vector sequences in sequence readings */

int main(int argc, char **argv) {
    int c;
    int word_length, num_diags, min_match, max_vector;
    double diag_score;
    int word_length_d, num_diags_d, min_match_d, min_left;
    int vector_primer_length;
    double diag_score_d;
    double cut_score_5, cut_score_5_d;
    double cut_score_3, cut_score_3_d;
    double max_prob, max_prob_d;
    int mode,i,tmode,rmode;
    char *fofn_p, *fofn_f, *fofn_i, *vf, *vector_seq;
    char expanded_fn[FILENAME_MAX+1];
    FILE *fp_p, *fp_f, *fp_i, *fp_vf;
    Vector_specs *v = NULL;

    fofn_p = fofn_f = fofn_i = vf = NULL;
    fp_p = fp_f = fp_i = fp_vf = NULL;

    max_vector = MAX_VECTOR_D;
    word_length_d = 4;
    num_diags_d = 7;
    min_match_d = 20;
    diag_score_d = 0.35;
    cut_score_5_d = 60.0;
    cut_score_3_d = 80.0;
    word_length = -1;
    num_diags = -1;
    min_match = -1;
    diag_score = -1;
    cut_score_5 = -1.0;
    cut_score_3 = -1.0;
    vector_primer_length = 1000;
    min_left = 0;
    mode = 0;
    tmode = 0;
    rmode = 0;
    max_prob_d = -1.0;
    max_prob = max_prob_d;
    set_dna_lookup();
    set_char_set(1); /* FIXME DNA*/

    while ((c = getopt(argc, argv, "w:n:d:l:L:R:p:f:m:M:P:v:V:i:schrtT")) != -1) {
	switch (c) {
	case 'w':
	    word_length = atoi(optarg);
	    if ( word_length > MAX_HASH_LENGTH ) word_length = MAX_HASH_LENGTH;
	    break;
	case 'n':
	    num_diags = atoi(optarg);
	    break;
	case 'd':
	    diag_score = atof(optarg);
	    break;
	case 'P':
	    max_prob = atof(optarg);
	    break;
	case 'l':
	    min_match = atoi(optarg);
	    break;
	case 'V':
	    vector_primer_length = atoi(optarg);
	    break;
	case 'L':
	    cut_score_5 = atof(optarg);
	    break;
	case 'R':
	    cut_score_3 = atof(optarg);
	    break;
	case 'm':
	    min_left = atoi(optarg);
	    break;
	case 'M':
	    max_vector = atoi(optarg);
	    break;
	case 's':
	    mode = SEQUENCE_VECTOR;
	    break;
	case 'c':
	    mode = CLONING_VECTOR;
	    break;
	case 'h':
	    mode = HGMP;
	    break;
	case 'r':
	    mode = VECTOR_REARRANGEMENT;
	    break;
	case 't':
	    tmode = TEST_ONLY;
	    break;
	case 'T':
	    tmode = TEST_ONLY_ALIGNMENTS;
	    break;
	case 'p':
	    fofn_p = optarg;
	    break;
	case 'v':
	    vf = optarg;
	    mode = SEQUENCE_VECTOR_PVF;
	    break;
	case 'i':
	    vf = optarg;
	    mode = TRANSPOSON;
	    break;
	case 'f':
	    fofn_f = optarg;
	    break;
	default:
	    usage(word_length_d, num_diags_d, diag_score_d, min_match_d, (int)cut_score_5_d,(int)cut_score_3_d);
	}
    }

    if (optind == argc)
	    usage(word_length_d, num_diags_d, diag_score_d, min_match_d, (int)cut_score_5_d,(int)cut_score_3_d);

    fofn_i = argv[optind];
    if (mode == 0) mode = SEQUENCE_VECTOR;
    if ( word_length < 0 ) word_length = word_length_d;
    if ( word_length > MAX_HASH_LENGTH ) word_length = MAX_HASH_LENGTH;
    if ( num_diags < 0 ) num_diags = num_diags_d;
    if ( diag_score < 0. ) diag_score = diag_score_d;
    if ( min_match < 0 ) min_match = min_match_d;
    if ( cut_score_5 < 0.0 ) cut_score_5 = cut_score_5_d;
    if ( cut_score_3 < 0.0 ) cut_score_3 = cut_score_3_d;
    if ( max_vector < MIN_VECTOR ) max_vector = MIN_VECTOR;
    if ( vector_primer_length < 5 ) vector_primer_length = 1000;

    fp_i = fopen(fofn_i, "r");
    if (fp_i == NULL ) {
	fprintf(stderr, "Failed to open input file of file names\n");
	return -1;
    }
    if ( fofn_p ) {
	fp_p = fopen(fofn_p, "w");
	if (fp_p == NULL ) {
	    fprintf(stderr, "Failed to open output file of passed file names\n");
	    return -1;
	}
    }
    if ( fofn_f ) {
	fp_f = fopen(fofn_f, "w");
	if (fp_f == NULL ) {
	    fprintf(stderr, "Failed to open output file of failed file names\n");
	    return -1;
	}
    }

    if ( vf ) {
      if (1 != expandpath(vf, expanded_fn)) {
	fprintf(stderr, "Failed to open vector-primer pair file\n");
	return -1;
      }
      fp_vf = fopen(expanded_fn, "r");
      if (fp_vf == NULL ) {
	fprintf(stderr, "Failed to open vector-primer pair file\n");
	return -1;
      }
      if ( !(v = get_vector_info(fp_vf)) ) return -1;
      fclose ( fp_vf);
    }


    if ( ! (vector_seq = (char *) xmalloc ( sizeof(char)*max_vector ))) return -1;

    if ( mode == HGMP ) {

	i = do_it_3p ( fp_i, fp_p, fp_f, 
		       cut_score_3, tmode );
    }
    else if ( mode == SEQUENCE_VECTOR ) {

	i = do_it_sv ( vector_seq, max_vector, fp_i, fp_p, fp_f, 
		       cut_score_5, cut_score_3, min_left, tmode, rmode );
    }
    else if ( mode == SEQUENCE_VECTOR_PVF ) {

	i = do_it_sv_pvf ( vector_seq, max_vector, fp_i, fp_p, fp_f,
		       cut_score_5, cut_score_3, min_left, 
		       vector_primer_length, tmode, rmode, v);
    }
    else if ( mode == CLONING_VECTOR ) {

	i = do_it_cv ( vector_seq, max_vector, fp_i, fp_p, fp_f, word_length, num_diags, 
		       diag_score, max_prob, tmode );

    }

    else if ( mode == TRANSPOSON ) {

	i = do_it_tr ( vector_seq, max_vector, fp_i, fp_p, fp_f,
		       cut_score_5, cut_score_3, min_match, tmode, rmode,
                       v);
    }

    else if ( mode == VECTOR_REARRANGEMENT ) {

	i = do_it_vr ( vector_seq, max_vector, fp_i, fp_p, fp_f, word_length, min_match,
		       tmode );
    }

    fprintf(stdout,"\n");
    xfree ( vector_seq );
    return 0;
}
