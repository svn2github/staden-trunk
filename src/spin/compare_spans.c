/* New matching span search program */

#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "sip_hash.h"
#include "readpam.h"
#include "misc.h" /* MIN */
#include "dna_utils.h"
#include "edge.h"

void make_reverse ( int **seq2_match, int window_length, 
		   int *n_matches, int seq2_len );

#define complement_base(base) (complementary_base[base])

int compare_spans (
		    char *seq1in,	/* the vertical sequence */
		    int seq1_len,	/* length seq1 */
		    int seq1_lreg,      /* left cutoff */
		    int seq1_rreg,      /* right cutoff */
		    char *seq2in,	/* the horizontal sequence */
		    int seq2_len,	/* length seq2 */
		    int seq2_lreg,      /* left cutoff */
		    int seq2_rreg,      /* right cutoff */
		    int window_length,	/* the window length */
		    int min_score,	/* the minimum score */
		    int **seq1_match,	/* match positions for seq1 */
		    int **seq2_match,	/* match positions for seq2 */
		    int **match_score,	/* match scores */
		    int max_matches,    /* maximum number of matches */
		    int same_seq	/* for internal forward matches */
		    )

{

/* Find all spans of length window_length that reach the score min_score.

	score [ column ] 
	= score [ trailing_column ]
	- dna_matrix [ seq1 [ trailing_row ] ] [ seq2 [ trailing_column ] ]
	+ dna_matrix [ seq1 [ leading_row ] ] [ seq2 [ leading_column ] ];

	process rows from -1 to seq1_len - 1
	use trailing_row = row - 1, leading_row = row + window_length - 1

 	process columns from seq2_len - 1 to 1 - window_length
	use trailing_column = column - 1, leading_column = column + window_length - 1
	ie columns are processed right to left because we need the score for the
	previous diagonal which is stored one array element to the left.

	suppose seq1_len = 6, seq2_len = 5, window_length = 3

	then we process rows from -1 to 5; and columns  4 to -1 but we need extra
	points filled because of the leading and trailing edges, as is explained:

	we need '-'s around the outside of the sequences so that the leading and
	trailing edges add and remove appropriate values.
	The first position we process is -1, seq2_len-1 which, in the example
	is -1,4; so its leading column is at 0,6 and its trailing column at -3,3.
	The last column in the first row is -2,-2, with trailing column -3,-3. So
	we need sequence data from -3 to 6 ie -window_length to seq2_len+window_length-2.
	The last row processed is 5 ie seq1_len-1. Its leading row reaches 7,0
	ie seq1_len+window_length-2.


	-  7,-1   7,0  7,1  7,2  7,3  7,4   7,5  7,6
	-  6,-1   6,0  6,1  6,2  6,3  6,4   6,5  6,6
                  ----------------------- 
	t  5,-1 | 5,0  5,1  5,2  5,3  5,4 | 5,5  5,6
	t  4,-1 | 4,0  4,1  4,2  4,3  4,4 | 4,5  4,6
	t  3,-1 | 3,0  3,1  3,2  3,3  3,4 | 3,5  3,6
	t  2,-1 | 2,0  2,1  2,2  2,3  2,4 | 2,5  2,6
	t  1,-1 | 1,0  1,1  1,2  1,3  1,4 | 1,5  1,6
	t  0,-1 | 0,0  0,1  0,2  0,3  0,4 | 0,5  0,6
                  ----------------------- 
	- -1,-1  -1,0 -1,1 -1,2 -1,3 -1,4  -1,5 -1,6
              -    t    t    t    t    t    -    -

	allocate long_score[] for storing the scores but set 
	*score to long_score[window_length] to simplify? the addressing.

	First we fill score[] with the score for a diagonal of '-'. Then
	when we process the first row it is as though we had just completed
	a set of diagonals filled with '-'.

	At the completion of a row score[i] represents the diagonal starting 
	at column i on that row.

	Arrgh! This is all very clever but now we want to be able to do it
	over subsections of the two sequences. How?

	Well its been done so the documentation above is now out of date.
*/

    int column, trailing_column, leading_column;
    int row, trailing_row, leading_row;
    int i, *score, *long_score, min_scorem1;
    int trailing_row_index, leading_row_index, *matrix_lrow, *matrix_trow;
    int match_number = 0;
    int temp_score;
    int **row_vectors, *left_edge, *long_left_edge;
    char *tmp_seq1,*seq1;
    char *tmp_seq2,*seq2;
    int tmp_seq1_len;
    int tmp_seq2_len;
    int j,pseq1,pseq2,window,half_window;

    long_score = long_left_edge = NULL;
    tmp_seq1 = tmp_seq2 = NULL;
    row_vectors = NULL;

    /* sanity checks */

    if ( ( window_length % 2 ) == 0 ) goto bail_out;
    if ( (seq1_rreg - seq1_lreg + 1) < window_length ) goto bail_out;
    if ( (seq2_rreg - seq2_lreg + 1) < window_length ) goto bail_out;

    half_window = window_length / 2;

    tmp_seq1_len = (seq1_rreg - seq1_lreg + 1) + window_length;
    tmp_seq2_len = (seq2_rreg - seq2_lreg + 1) + window_length;

    if (NULL == (tmp_seq1 = (char *)xmalloc(tmp_seq1_len * sizeof(char))))
	goto bail_out;
    if (NULL == (tmp_seq2 = (char *)xmalloc(tmp_seq2_len * sizeof(char))))
	goto bail_out;

    if (NULL == (row_vectors = (int **)xmalloc(char_set_size * sizeof(int *))))
	goto bail_out;
    for (i = 0; i < char_set_size; i++) {
	row_vectors[i] = &score_matrix[i][0];
    }

    min_scorem1 = min_score - 1;

    if ( ! ( long_score = (int *) xmalloc ( sizeof(int) * tmp_seq2_len))) {
	goto bail_out;
    }
    if (NULL == (long_left_edge = (int *) xmalloc ( sizeof(int) * tmp_seq1_len))) {
	goto bail_out;
    }
    /* set seq1 and seq2 to be the actual sequence arrays so they can have 
       negative elements. Ditto scores.
    */

    seq1 = &tmp_seq1[half_window+1];
    seq2 = &tmp_seq2[half_window+1];

    score = &long_score[half_window+1];
    left_edge = &long_left_edge[half_window+1];

    /*  fill in temp seqs */

    for ( i=-half_window-1, j=seq1_lreg-1-half_window-1;
	  i<seq1_rreg-seq1_lreg+1+half_window;
	  i++, j++ ) {

	if ( j < 0 ) {
	    seq1[i] = char_lookup['-'];
	}
	else if ( j > seq1_len-1) {
	    seq1[i] = char_lookup['-'];
	}
	else {
	    seq1[i] = char_lookup[seq1in[j]];
	}
    }

    for ( i=-half_window-1, j=seq2_lreg-1-half_window-1;
	  i<seq2_rreg-seq2_lreg+1+half_window;
	  i++, j++ ) {

	if ( j < 0 ) {
	    seq2[i] = char_lookup['-'];
	}
	else if ( j > seq2_len-1) {
	    seq2[i] = char_lookup['-'];
	}
	else {
	    seq2[i] = char_lookup[seq2in[j]];
	}
    }


    /*	Set the scores for the first row */

    for ( column = -1; column < seq2_rreg - seq2_lreg + 1; column++ ) {
	for ( window=0, score[column]=0, pseq1=-half_window-1, pseq2=column-half_window;
	     window < window_length;
	     window++,pseq1++,pseq2++) {
	    score [ column ] += score_matrix[seq2[pseq2]][seq1[pseq1]];
	}
    }

    /*	Set the scores for the left_edge */

    for ( row = -1; row < seq1_rreg - seq1_lreg + 1; row++ ) {
	for ( window=0,left_edge[row]=0,pseq2=-half_window-1,pseq1=row-half_window;
	     window<window_length;
	     window++,pseq1++,pseq2++) {

	    left_edge [ row ] += score_matrix[seq2[pseq2]][seq1[pseq1]];
	}
    }

    /* loop for each row */

    for ( row = 0, 
	 trailing_row = row - half_window - 1,
	 leading_row = row + half_window;
	 row < seq1_rreg-seq1_lreg+1; 
	 row++, trailing_row++, leading_row++ ) {

	score [ -1 ] = left_edge [ row - 1 ];
	trailing_row_index = seq1 [ trailing_row]; 
	leading_row_index = seq1 [leading_row]; 
	matrix_trow = row_vectors [ trailing_row_index ];
	matrix_lrow = row_vectors [ leading_row_index ]; 

	/* loop right to left for each column */

	for ( column = seq2_rreg-seq2_lreg,
	     trailing_column = column -half_window - 1,
	     leading_column = column + half_window;
	     column > -1; 
	     column--, trailing_column--, leading_column-- ) {

	    temp_score =
		score [ column-1 ]
		- matrix_trow [seq2 [trailing_column]]
	      	+ matrix_lrow [seq2 [leading_column]];
	    if ( min_scorem1 < temp_score ) {

		if ( match_number == max_matches ) {

		    sip_realloc_matches(seq1_match, seq2_match,
					match_score, &max_matches);
		}

		if ( ! ( same_seq && ( row == column ) ) ) {

		    (*seq1_match) [ match_number ] = row + seq1_lreg - half_window;
		    (*seq2_match) [ match_number ] = column + seq2_lreg - half_window;
		    (*match_score) [ match_number ] = temp_score;
		    match_number++;
		}

	    }
	    score [ column ] = temp_score;
	}
    }

    if (long_score) xfree (long_score);
    if (tmp_seq1) xfree(tmp_seq1);
    if (tmp_seq2) xfree(tmp_seq2);
    if (row_vectors) xfree(row_vectors);
    if (long_left_edge) xfree(long_left_edge);
    return match_number;

 bail_out:

    if (long_score) xfree (long_score);
    if (tmp_seq1) xfree(tmp_seq1);
    if (tmp_seq2) xfree(tmp_seq2);
    if (row_vectors) xfree(row_vectors);
    if (long_left_edge) xfree(long_left_edge);
    return -1;
}
/************************************************************/

int cmpspn (
	    char *sense,	/* the orientation for seq2 */
	    int *min_mat,	/* the minimum match length */
	    int **seq1_match,	/* positions of matches in seq1 */
	    int **seq2_match,	/* positions of matches in seq2 */
	    int **match_score,	/* match scores */
	    int *max_match,	/* maximum number of matches */
	    int *window_l,	/* length of window */
	    char *seq1,		/* seq1 */
	    char *seq2,		/* seq2 */
	    int *seq1_l, 	/* size of seq1 */
	    int *seq2_l,	/* size of seq2 */
	    int seq1_lreg,      /* left cutoff of seq1 */
	    int seq1_rreg,      /* right cutoff of seq1 */
	    int seq2_lreg,      /* left cutoff of seq2 */
	    int seq2_rreg,      /* right cutoff of seq2 */
	    int same_seq        /* whether seq1 & seq2 are identical */
	    ) {



/*	Main interface to the sequence span comparison method

	Note if sense is 'r' then we must reverse and complement seq2 and
	set its results to reflect positions in the original sense.

	Note that if seq2 is NULL then we are looking for internal repeats
	and we copy seq1 to seq2 temporarily then reset seq2 to NULL.
	For the FORTRAN interface we also allow seq2_len == 0 to have
	the same effect.

*/

#define MAXERRORMESS 200 

    int start_pos,seq2NULL;
    int *start_pos_p, n_matches,min_match;
    int seq1_len, seq2_len, window_length, max_matches;
    char *error_mess_ptr, error_mess[MAXERRORMESS];

    min_match = *min_mat;
    seq1_len = *seq1_l;
    seq2_len = *seq2_l;
    window_length = *window_l;
    max_matches = *max_match;
    seq2NULL = 0;
    start_pos = 1;
    start_pos_p = &start_pos;
    min_match = *min_mat;
    error_mess_ptr = error_mess;

    /* if seq2 is NULL or seq2_len is zero then copy seq1 to seq2 */
    /* now done in sip_cmds.c */

    /* if sense is reverse then reverse and complement seq2 */

    if ( 'r' == *sense ) {
	(void) complement_seq( seq2, seq2_len);
    }
    

    /* do the comparison	*/

    n_matches = compare_spans ( seq1, seq1_len, seq1_lreg, seq1_rreg,
				seq2, seq2_len, seq2_lreg, seq2_rreg,
				window_length, min_match, 
				seq1_match, seq2_match, match_score, 
				max_matches, same_seq);
    
    /* if reverse matches make the numbers relative to original orientation */

    if ( 'r' == *sense ) {
	(void) make_reverse ( seq2_match, window_length, 
			     &n_matches, seq2_len ); 	
    }

    /* if looking for internal matches remove duplicates */

    if ( same_seq ) {
	(void) sip_remdup ( seq1_match, seq2_match, match_score, &n_matches );
	*seq2_l = 0;
    }

    return n_matches;
}

/************************************************************/

void make_reverse (int **seq2_match, 
		   int window_length,
		   int *n_matches, 
		   int seq2_len ) {
    
    int i;
    
    for (i = 0; i< *n_matches; i++) {
	(*seq2_match)[i] = seq2_len - (*seq2_match)[i] - window_length+2;
	
    }
}
