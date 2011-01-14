#include <stdio.h>
#include <float.h>
#include <math.h>
#include "misc.h"
#include "dna_utils.h"
#include "genetic_code.h"
#include "edge.h"
#include "codon_content.h"
#include "array_arith.h"
#include "trna_search.h"

/* gene search by content */


TrnaSpec *init_TrnaSpec ( 

    int max_trna_length,
    int min_trna_length,
    int max_intron_length,
    int min_intron_length,
    int max_tu_loop_length,
    int min_tu_loop_length,
    int min_aa_to_ac_length,
    int max_aa_to_ac_length,
    int min_acs_to_ace_length,
    int max_var_loop_length,
    int min_aa_score,
    int min_ac_score,
    int min_tu_score,
    int min_d_score,
    int min_total_bp_score,
    int min_total_cb_score,
    int con_base_scores[]) {
    
    TrnaSpec *t;
    int con_bases1[18] = {3,2,1,0,0,0,1,3,0,1,2,3,3,1,0,0,1,1};
    int con_bases2[18] = {3,2,3,0,2,0,3,3,0,3,2,3,3,1,2,0,3,1};
    int con_base_offsets[18] = {7,9,10,13,14,-10,1,2,6,-5,0,1,2,3,4,5,7,8};
    int i; 

    if ( ( NULL == ( t = (TrnaSpec* ) xmalloc ( sizeof (TrnaSpec) )))) return NULL;

    t->max_trna_length    = max_trna_length;
    t->min_trna_length    = min_trna_length;
    t->max_intron_length  = max_intron_length;
    t->min_intron_length  = min_intron_length;
    t->max_tu_loop_length = max_tu_loop_length;
    t->min_tu_loop_length = min_tu_loop_length;
    t->min_aa_to_ac_length = min_aa_to_ac_length;
    t->max_aa_to_ac_length = max_aa_to_ac_length;
    t->min_acs_to_ace_length = min_acs_to_ace_length;
    t->max_var_loop_length = max_var_loop_length;
    t->min_aa_score = min_aa_score;
    t->min_ac_score = min_ac_score;
    t->min_tu_score = min_tu_score;
    t->min_d_score = min_d_score;
    t->min_total_bp_score = min_total_bp_score;
    t->min_total_cb_score = min_total_cb_score;
    for ( i=0; i<18; i++ ) {
	t->con_base_scores[i] = con_base_scores[i];
	t->con_bases1[i] = con_bases1[i];
	t->con_bases2[i] = con_bases2[i];
	t->con_base_offsets[i] = con_base_offsets[i];
    }
    return t;
}


void trna_base_scores ( TrnaRes *r, TrnaSpec *t ) {

    int j,k,b;

    /* check for presence of conserved bases */

    /* fudge factors to fit fortran (done earlier)

     * r->aa_right += 1;
     * r->ac_left  += 4;
     * r->ac_right -= 4;
     * r->tu_right -= 4;
     * r->tu_left  += 4;
    */

    /*  test those from left end */

    for ( j=0, r->total_cb_score = 0; j<5; j++ ) {
	k = r->aa_left + t->con_base_offsets [ j ];
	b = char_lookup [ r->seq [ k ] ];
	if ( ( b == t->con_bases1 [ j ] ) || ( b == t->con_bases2 [ j ] ) ) {
	    r->total_cb_score += t->con_base_scores [ j ];
	}
    }

    /* those from anticodon */

    for ( j=5; j<9; j++ ) {
	k = r->ac_left + t->con_base_offsets [ j ];
	b = char_lookup [ r->seq [ k ] ];
	if ( ( b == t->con_bases1 [ j ] ) || ( b == t->con_bases2 [ j ] ) ) {
	    r->total_cb_score += t->con_base_scores [ j ];
	}
    }

    /* those from tu loop */

    for ( j=9; j<18; j++ ) {
	k = r->tu_left + t->con_base_offsets [ j ];
	b = char_lookup [ r->seq [ k ] ];
	if ( ( b == t->con_bases1 [ j ] ) || ( b == t->con_bases2 [ j ] ) ) {
	    r->total_cb_score += t->con_base_scores [ j ];
	}
    }
}

TrnaRes *init_TrnaRes (void) {
    TrnaRes *r;
    if ( ( NULL == ( r = ( TrnaRes* ) xmalloc ( sizeof (TrnaRes) )))) return NULL;
    return r;
}

TrnaRes *init_TrnaRes1 ( 
    char *seq,
    int seq_length,
    int aa_right,
    int aa_left,
    int ac_left,
    int ac_right,
    int tu_right,
    int tu_left,
    int intron_length,
    int aa_score,
    int ac_score,
    int tu_score,
    int d_score,
    int total_bp_score,
    int total_cb_score ) {
    TrnaRes *r;

    if ( ( NULL == ( r = ( TrnaRes* ) xmalloc ( sizeof (TrnaRes) )))) return NULL;

    r->seq = seq;
    r->seq_length = seq_length;
    r->aa_right = aa_right;
    r->aa_left = aa_left;
    r->ac_left = ac_left;
    r->ac_right = ac_right;
    r->tu_right = tu_right;
    r->tu_left = tu_left;
    r->intron_length = intron_length;
    r->aa_score = aa_score;
    r->ac_score = ac_score;
    r->tu_score = tu_score;
    r->d_score  = d_score;
    r->total_bp_score = total_bp_score;
    r->total_cb_score = total_cb_score;
    return r;
}



void draw_trna ( TrnaRes *r ) {

    char matrix[35][35], scores_v[5][5], scores_vh[5][5];
    int i,j,k,j2;
    int l,m,length,length_upper, length_lower, length_outer;
    int data;

    /* draw the trna cloverleaf */

    /* fudge factors to fit fortran 

    r->aa_right += 1;
    r->ac_left  += 4;
    r->ac_right -= 4;
    r->tu_right -= 4;
    r->tu_left  += 4;
*/


    vmessage("aa_right %d ac_left %d ac_right %d tu_right %d tu_left %d\n",
	   r->aa_right,r->ac_left,r->ac_right,r->tu_right,r->tu_left);


    vmessage("trna start %d\n",r->aa_left+1);
    vmessage("trna end %d\n",r->aa_right+1);
    vmessage("aa_score %d\n",r->aa_score);
    vmessage("ac_score %d\n",r->ac_score);
    vmessage("d_score %d\n",r->d_score);
    vmessage("tu_score %d\n",r->tu_score);
    vmessage("total_bp_score %d\n",r->total_bp_score);
    vmessage("total_cb_score %d\n",r->total_cb_score);

    for ( i=0; i<35; i++ ) {
	for( j=0; j<35; j++ ) {
	    matrix[i][j] = ' ';
	}
    }
    for ( i=0; i<5; i++ ) {
	for( j=0; j<5; j++ ) {
	    scores_v[i][j] = ' ';
	    scores_vh[i][j] = ' ';
	}
    }

    scores_v[0][3] = '-';
    scores_v[3][0] = '-';
    scores_v[1][2] = '-';
    scores_v[2][1] = '-';
    scores_v[2][3] = '+';
    scores_v[3][2] = '+';
    scores_vh[0][3] = '|';
    scores_vh[3][0] = '|';
    scores_vh[1][2] = '|';
    scores_vh[2][1] = '|';
    scores_vh[2][3] = '+';
    scores_vh[3][2] = '+';

    /* aa */

    if ( r->aa_right < r->seq_length ) matrix[18][5] = r->seq[r->aa_right];
    j = r->aa_right -1 ;
    k = r->aa_left;
    for ( i=6; i<=12; i++,j--,k++ ) {
	matrix[16][i] = r->seq[k];
	matrix[18][i] = r->seq[j];
	matrix[17][i] = scores_v [ char_lookup [r->seq[k]]] [ char_lookup [r->seq[j]]];
    }

    /* gap between aa and d stems */

    matrix[15][13] = r->seq[r->aa_left+7];
    matrix[14][14] = r->seq[r->aa_left+8];

    /* ac */

    k = r->ac_left - 4;
    j = r->ac_right + 4;
    for ( i=19; i<=23; i++,k++,j-- ) {
	matrix[15][i] = r->seq[k];
	matrix[17][i] = r->seq[j];
	matrix[16][i] = scores_v [ char_lookup [r->seq[k]]] [ char_lookup [r->seq[j]]];
    }

    /* ac loop */

    matrix[14][24] = r->seq[r->ac_left+1];
    matrix[14][25] = r->seq[r->ac_left+2];
    matrix[18][25] = r->seq[r->ac_left+6];
    matrix[18][24] = r->seq[r->ac_right-1];
    k = r->ac_left + 3;
    for ( i=15; i<=17; i++, k++ ) {
	matrix[i][26] = r->seq[k];
    }

    /* gap between d stem and ac loop */

    matrix[14][18] = r->seq[r->ac_left-5];

    /* tu stem */

    k = r->tu_right + 4;
    l = r->tu_left - 4;
    for ( j=19; j<=23; j++,k--,l++) {
	matrix[j][13] = r->seq[k];
	matrix[j][15] = r->seq[l];
	matrix[j][14] = scores_vh [ char_lookup [r->seq[k]]] [ char_lookup [r->seq[l]]];
    }

    /* tu loop */

    /* length of tu loop */
    length = k - l + 1;

/*    printf("%d %d %d %d\n",r->tu_right,r->tu_left,k,l); */

    if ( length > 18 ) {
#ifdef DEBUG
	printf("scream 1\n");
#endif
	return;
    }
    if ( length < 3 ) {
#ifdef DEBUG
	printf("scream 2\n");
#endif
	return;
    }
    /* number in outer */
    length_outer = ( length - 1 ) / 2;

    for ( m=24,j=1; j<=length_outer; j++, k--, l++, m++ ) {
	matrix[m][12] = r->seq[k];
	matrix[m][16] = r->seq[l];
    }
    matrix[m][15] = r->seq[l++];
    if ( !(length%2) ) matrix[m][14] = r->seq[l];

    /* d stem */

    k = r->aa_left + 9;
    l = r->ac_left - 6;
    length = l - k + 1;
    j2 = 4;
    m = 13;
    if ( length <= 10 ) {
	/* need stem of 3 not 4 */
	m = 12;
	j2 = 3;
    }
    for ( j = 1; j<= j2; j++, k++, l--, m-- ) {
	matrix[m][15] = r->seq[k];
	matrix[m][17] = r->seq[l];
	matrix[m][16] = scores_vh [ char_lookup [r->seq[k]]] [ char_lookup [r->seq[l]]];
    }

    /* d loop */

    length = l - k + 1;
    if ( length > 18 ) {
#ifdef DEBUG
	printf("scream 3\n");
#endif
	return;
    }
    if ( length < 3 ) {
#ifdef DEBUG
	printf("scream 4\n");
#endif
	return;
    }
    length_outer = ( length - 1 ) /2;
    for ( j=1; j<=length_outer; j++, k++, l--, m-- ) {
	matrix[m][14] = r->seq[k];
	matrix[m][18] = r->seq[l];
    }
    matrix[m][17] = r->seq[l--];
    if ( !(length%2) ) matrix[m][16] = r->seq[l];

    /* variable loop */

    length = r->tu_left - 5 - r->ac_right - 4;
    if ( length > 30 ) {
	vmessage("variable loop truncated\n");
	length = 30;
    }
    if ( length < 3 ) {
#ifdef DEBUG
	printf("scream 5\n");
#endif
	return;
    }
    /* number of elements in lower diagonal */
    length_lower = ( length / 2) - 1;
    if ( length_lower == 0 ) length_lower = 1;
    /* number of elements in upper diagonal */
    length_upper = ( length + 1 ) / 2;
    /* do lower */
    for ( j = 1, k = 18, l = 19, m = r->ac_right + 5; j<=length_lower; j++, k++, l++, m++ ) {
	matrix[k][l] = r->seq[m];
    }
    k++;
    l--;
    if ( length > 3 ) matrix[k][l] = r->seq[m];
    /* do upper */
    for ( j = 1, k = 19, l = 16, m = r->tu_left - 5; j<=length_upper; j++, k++, l++, m-- ) {
	matrix[k][l] = r->seq[m];
    }

    /* output */

    for ( i=0; i<35; i++ ) {
	for ( j=0,data=0; j<35; j++ ) {
	    if ( matrix[j][i] != ' ' ) data = 1;
	}
	if ( data ) {
	    for ( j=0; j<35; j++ ) {
		vmessage("%c",matrix[j][i]);
	    }
	    vmessage("\n");
	}
    }
}

void trna_draw ( TrnaRes *r ) {

    char matrix[35][35], scores_v[5][5], scores_vh[5][5];
    int i,j,k,j2;
    int l,m,length,length_upper, length_lower, length_outer;
    int data;

    /* draw the trna cloverleaf */

    /* fudge factors to fit fortran 

    r->aa_right += 1;
    r->ac_left  += 4;
    r->ac_right -= 4;
    r->tu_right -= 4;
    r->tu_left  += 4;
*/

#ifdef DEBUG
    printf("aa_right %d ac_left %d ac_right %d tu_right %d tu_left %d\n",
	   r->aa_right,r->ac_left,r->ac_right,r->tu_right,r->tu_left);

    printf("trna start %d\n",r->aa_left+1);
    printf("trna end %d\n",r->aa_right+1);
    printf("aa_score %d\n",r->aa_score);
    printf("ac_score %d\n",r->ac_score);
    printf("d_score %d\n",r->d_score);
    printf("tu_score %d\n",r->tu_score);
    printf("total_bp_score %d\n",r->total_bp_score);
    printf("total_cb_score %d\n",r->total_cb_score);
#endif

    for ( i=0; i<35; i++ ) {
	for( j=0; j<35; j++ ) {
	    matrix[i][j] = ' ';
	}
    }
    for ( i=0; i<5; i++ ) {
	for( j=0; j<5; j++ ) {
	    scores_v[i][j] = ' ';
	    scores_vh[i][j] = ' ';
	}
    }

    scores_v[0][3] = '-';
    scores_v[3][0] = '-';
    scores_v[1][2] = '-';
    scores_v[2][1] = '-';
    scores_v[2][3] = '+';
    scores_v[3][2] = '+';
    scores_vh[0][3] = '|';
    scores_vh[3][0] = '|';
    scores_vh[1][2] = '|';
    scores_vh[2][1] = '|';
    scores_vh[2][3] = '+';
    scores_vh[3][2] = '+';

    /* aa */

    if ( r->aa_right < r->seq_length ) matrix[18][5] = r->seq[r->aa_right];
    j = r->aa_right -1 ;
    k = r->aa_left;
    for ( i=6; i<=12; i++,j--,k++ ) {
	matrix[16][i] = r->seq[k];
	matrix[18][i] = r->seq[j];
	matrix[17][i] = scores_v [ char_lookup [r->seq[k]]] [ char_lookup [r->seq[j]]];
    }

    /* gap between aa and d stems */

    matrix[15][13] = r->seq[r->aa_left+7];
    matrix[14][14] = r->seq[r->aa_left+8];

    /* ac */

    k = r->ac_left - 4;
    j = r->ac_right + 4;
    for ( i=19; i<=23; i++,k++,j-- ) {
	matrix[15][i] = r->seq[k];
	matrix[17][i] = r->seq[j];
	matrix[16][i] = scores_v [ char_lookup [r->seq[k]]] [ char_lookup [r->seq[j]]];
    }

    /* ac loop */

    matrix[14][24] = r->seq[r->ac_left+1];
    matrix[14][25] = r->seq[r->ac_left+2];
    matrix[18][25] = r->seq[r->ac_left+6];
    matrix[18][24] = r->seq[r->ac_right-1];
    k = r->ac_left + 3;
    for ( i=15; i<=17; i++, k++ ) {
	matrix[i][26] = r->seq[k];
    }

    /* gap between d stem and ac loop */

    matrix[14][18] = r->seq[r->ac_left-5];

    /* tu stem */

    k = r->tu_right + 4;
    l = r->tu_left - 4;
    for ( j=19; j<=23; j++,k--,l++) {
	matrix[j][13] = r->seq[k];
	matrix[j][15] = r->seq[l];
	matrix[j][14] = scores_vh [ char_lookup [r->seq[k]]] [ char_lookup [r->seq[l]]];
    }

    /* tu loop */

    /* length of tu loop */
    length = k - l + 1;

/*    printf("%d %d %d %d\n",r->tu_right,r->tu_left,k,l); */

    if ( length > 18 ) {
#ifdef DEBUG
	printf("scream 1\n");
#endif
	return;
    }
    if ( length < 3 ) {
#ifdef DEBUG
	printf("scream 2\n");
#endif
	return;
    }
    /* number in outer */
    length_outer = ( length - 1 ) / 2;

    for ( m=24,j=1; j<=length_outer; j++, k--, l++, m++ ) {
	matrix[m][12] = r->seq[k];
	matrix[m][16] = r->seq[l];
    }
    matrix[m][15] = r->seq[l++];
    if ( !(length%2) ) matrix[m][14] = r->seq[l];

    /* d stem */

    k = r->aa_left + 9;
    l = r->ac_left - 6;
    length = l - k + 1;
    j2 = 4;
    m = 13;
    if ( length <= 10 ) {
	/* need stem of 3 not 4 */
	m = 12;
	j2 = 3;
    }
    for ( j = 1; j<= j2; j++, k++, l--, m-- ) {
	matrix[m][15] = r->seq[k];
	matrix[m][17] = r->seq[l];
	matrix[m][16] = scores_vh [ char_lookup [r->seq[k]]] [ char_lookup [r->seq[l]]];
    }

    /* d loop */

    length = l - k + 1;
    if ( length > 18 ) {
#ifdef DEBUG
	printf("scream 3\n");
#endif
	return;
    }
    if ( length < 3 ) {
#ifdef DEBUG
	printf("scream 4\n");
#endif
	return;
    }
    length_outer = ( length - 1 ) /2;
    for ( j=1; j<=length_outer; j++, k++, l--, m-- ) {
	matrix[m][14] = r->seq[k];
	matrix[m][18] = r->seq[l];
    }
    matrix[m][17] = r->seq[l--];
    if ( !(length%2) ) matrix[m][16] = r->seq[l];

    /* variable loop */

    length = r->tu_left - 5 - r->ac_right - 4;
    if ( length > 30 ) {
#ifdef DEBUG
	printf("variable loop truncated\n");
#endif
	length = 30;
    }
    if ( length < 3 ) {
#ifdef DEBUG
	printf("scream 5\n");
#endif
	return;
    }
    /* number of elements in lower diagonal */
    length_lower = ( length / 2) - 1;
    if ( length_lower == 0 ) length_lower = 1;
    /* number of elements in upper diagonal */
    length_upper = ( length + 1 ) / 2;
    /* do lower */
    for ( j = 1, k = 18, l = 19, m = r->ac_right + 5; j<=length_lower; j++, k++, l++, m++ ) {
	matrix[k][l] = r->seq[m];
    }
    k++;
    l--;
    if ( length > 3 ) matrix[k][l] = r->seq[m];
    /* do upper */
    for ( j = 1, k = 19, l = 16, m = r->tu_left - 5; j<=length_upper; j++, k++, l++, m-- ) {
	matrix[k][l] = r->seq[m];
    }

    /* output */

    for ( i=0; i<35; i++ ) {
	for ( j=0,data=0; j<35; j++ ) {
	    if ( matrix[j][i] != ' ' ) data = 1;
	}
	if ( data ) {
	    for ( j=0; j<35; j++ ) {
		printf("%c",matrix[j][i]);
	    }
	    printf("\n");
	}
    }
}

int realloc_trna(TrnaRes ***r, int *max_trna)
{
    int trna_inc = 100;
    int prev = *max_trna;
    int i;

    *max_trna += trna_inc;
    
    if (NULL == ( (*r) = (TrnaRes **)realloc((*r), *max_trna * 
					     sizeof(TrnaRes *))))
	return -1;

    for (i = prev; i < *max_trna; i++) {
	if (( (*r)[i] = init_TrnaRes ()) == NULL) 
	    return -1;
    }

    return 0;
}


int do_trna_search ( char seq[], int seq_length, int user_start, int user_end,
		     TrnaSpec *t, TrnaRes ***r, int *nmatch, 
		     int *max_total_bp_score) {

    int aa_left, aa_right, max_aa_start, aa_left_start, min_aa_end, max_aa_end;
    int aa_right_end, aa_score;
    int tu_number, tu_right, tu_left, tu_score, tu_left_match[10], tu_match_score[10];
    int tu_right_match=0, tu_match_number;
    int i,j,start,end,intron_length;
    int ac_min_start, ac_max_start, ac_left, d_left, d_right,d_score;
    int ac_right_start, ac_right, lac, rac, ac_score, ac_right_end;
    int base_pair_score [ 25 ];
    int total_base_pair;
    int max_trna = MAX_TRNA;

    *nmatch = 0;

    fill_int_array ( base_pair_score, 25, 0 );
    base_pair_score [3] = 2;
    base_pair_score [7] = 2;
    base_pair_score [11] = 2;
    base_pair_score [13] = 1;
    base_pair_score [15] = 2;
    base_pair_score [17] = 1;

    start = user_start - 1;
    end   = user_end - 1;

    /* loop for all aa stem left starts */

    max_aa_start = end - ( t->min_trna_length - 1 );


    for ( aa_left_start = start; aa_left_start <= max_aa_start; aa_left_start++ ) {

	/* loop for all aa stem right ends */

	min_aa_end = aa_left_start + t->min_trna_length - 1;
	max_aa_end = MIN ( aa_left_start + t->max_trna_length + t->max_intron_length - 1, end);

	for ( aa_right_end = min_aa_end; aa_right_end <= max_aa_end; aa_right_end++ ) {

	    /* get the aa score */
	    for ( aa_left = aa_left_start, 
		 aa_right = aa_right_end, 
		 aa_score = 0,
		 i=0; i<7; 
		 aa_left++, aa_right--, i++ )  {

		aa_score += base_pair_score [ char_lookup [ seq [ aa_left ]] 
					     + char_lookup [ seq [ aa_right ] ] * 5 ];
	    }

	    if ( aa_score >= t->min_aa_score ) {

		/* do the tu loop */


		for ( i = t->min_tu_loop_length, tu_number = 0; i <= t->max_tu_loop_length; i++ ) {
		    tu_right = aa_right;
		    tu_left  = aa_right - 9 - i;
		    tu_score = 0;
		    for ( j=0; j<5; j++, tu_left++, tu_right-- ) {

			tu_score += base_pair_score [ char_lookup [ seq [ tu_left ]] 
					     + char_lookup [ seq [ tu_right ] ] * 5 ];
		    }

		    if ( tu_score >= t->min_tu_score ) {

			tu_left_match [ tu_number ] = tu_left - 5;
			tu_match_score [ tu_number ] = tu_score;
			tu_right_match = aa_right;
			tu_number++;
		    }
		}

		/* loop for all tu stems to find ac stem */

		for ( tu_match_number = 0; tu_match_number < tu_number; tu_match_number++ ) {

		    /* try all ac left starts */

		    ac_min_start = aa_left_start + t->min_aa_to_ac_length;
		    ac_max_start = MIN ( ( tu_left_match [ tu_match_number ] -
					t->min_aa_to_ac_length ), ( aa_left_start + t->max_aa_to_ac_length ) );

		    for ( ac_left = ac_min_start; ac_left <= ac_max_start; ac_left++ ) {

			/* do the d stem first */

			d_left = aa_left_start + 8;
			d_right = ac_left -1;
			for ( i=0, d_score = 0; i<5; i++ ) {
			    d_left++;
			    d_right--;

			    d_score += base_pair_score [ char_lookup [ seq [ d_left ]] 
					     + char_lookup [ seq [ d_right ] ] * 5 ];
			}

			if ( d_score >= t->min_d_score ) {

			    /* try all ac right end positions */

			    ac_right_start = MAX ( ( ac_left + t->min_acs_to_ace_length ),
						 ( tu_left_match [ tu_match_number ] -
						  t->max_var_loop_length ));
			    ac_right_end = MIN ( ( ac_left + t->min_acs_to_ace_length + t->max_intron_length),
						( tu_left_match [ tu_match_number ] - 4 ));

			    for ( ac_right = ac_right_start; ac_right <= ac_right_end; ac_right++ ) {
				lac = ac_left - 1;
				rac = ac_right + 1;
				for ( i=0, ac_score = 0; i<5; i++ ) {
				    lac++;
				    rac--;

				    ac_score += base_pair_score [ char_lookup [ seq [ lac ]] 
					     + char_lookup [ seq [ rac ] ] * 5 ];
				}

				if ( ac_score >= t->min_ac_score ) {

				    /* we have got all stems !!! */

				    /* intron length sensisble ? */

				    intron_length = ac_right - ac_left - 16;

				    if ( ( ( intron_length == 0 ) ||
					 ( intron_length >= t->min_intron_length ) ) &&
					 ( ( aa_right_end - aa_left_start + 1 - intron_length ) <=
					  t->max_trna_length )) {

					/* high enough overall base pairing score ? */
					total_base_pair = aa_score + ac_score + 
					    d_score + tu_match_score [ tu_match_number ];
					if ( total_base_pair >= t->min_total_bp_score ) {
					    /* fudge factors to fit fortran 

					     *  r->aa_right += 1;
					     *  r->ac_left  += 4;
					     *  r->ac_right -= 4;
					     *  r->tu_right -= 4;
					     *  r->tu_left  += 4;
					     */

					    (*r)[*nmatch]->seq = seq;
					    (*r)[*nmatch]->seq_length = seq_length;
					    (*r)[*nmatch]->aa_right = aa_right_end + 1;
					    (*r)[*nmatch]->aa_left = aa_left_start;
					    (*r)[*nmatch]->ac_left = ac_left + 4;
					    (*r)[*nmatch]->ac_right = ac_right - 4;
					    (*r)[*nmatch]->tu_right = tu_right_match - 4;
					    (*r)[*nmatch]->tu_left = tu_left_match[tu_match_number] + 4;

					    /* do conserved base search in an odd place ! */

					    if ( t->min_total_cb_score ) {
					      trna_base_scores ( (*r)[*nmatch], t );
					      if ( (*r)[*nmatch]->total_cb_score < t->min_total_cb_score ) continue;
					    }
					    (*r)[*nmatch]->intron_length = intron_length;
					    (*r)[*nmatch]->aa_score = aa_score;
					    (*r)[*nmatch]->ac_score = ac_score;
					    (*r)[*nmatch]->tu_score = tu_match_score[tu_match_number];
					    (*r)[*nmatch]->d_score  = d_score;
					    (*r)[*nmatch]->total_bp_score = total_base_pair;

					    if ((*r)[*nmatch]->total_bp_score >
						*max_total_bp_score) {
						*max_total_bp_score = (*r)[*nmatch]->total_bp_score;
					    }
					    (*nmatch)++;

					    if (*nmatch >= max_trna) {
#ifdef DEBUG
						printf("REALLOC nmatch %d max_trna %d\n",
						       *nmatch, max_trna);
#endif

						if (-1 == realloc_trna(r, &max_trna))
						    return -1;

					    }
					    /* really we need to store up the results 
					    and return them. Then trna_draw is not
					    called from here */
					}
				    }
				}
			    }
			}
		    }	
		}
	    }
	}
    }
    return 0;
}


int trna_search ( char seq[], int seq_length, int user_start, int user_end,
		  TrnaRes ***results, int *nmatch, int *max_total_bp_score, 
		  TrnaSpec **t) {

    /* TrnaSpec *t; */
    /* TrnaRes **results; */
#define NIP
#ifdef NIP4
    int max_trna_length = 100; /* 70 - 130, 92 */
    int min_trna_length = 60;  /* */
    int max_intron_length =20;  /* 0 - 60, 0 */
    int min_intron_length = 0; /* 0 - 60, 0 */
    int max_tu_loop_length = 8; /* 6 - 12, 9 */
    int min_tu_loop_length = 4; /* 4 - 12, 6 */
    int min_aa_to_ac_length = 25; /* 19 in nip but this looks correct as the stem aa length forgotten*/
    int max_aa_to_ac_length = 29; /* 35 in nip but this is correct if max of 4 extra in d loop*/
    int min_acs_to_ace_length = 16;
    int max_var_loop_length = 28;
    int min_aa_score = 13;   /* 0 - 14, 11 */
    int min_ac_score = 10;  /* 0 - 10, 8 */
    int min_tu_score = 10;   /* 0 - 10, 8 */
    int min_d_score = 6;    /* 0 - 8, 3 */
    int min_total_bp_score = 0; /* 30 - 44, 30 */
    int min_total_cb_score = 18;
#endif

#ifdef NIP
    /* these are now set to find all the chmpxx trnas without introns */
    int max_trna_length = 92;
    int min_trna_length = 60;  /* */
    int max_intron_length = 0; 
    int min_intron_length = 0; 
    int max_tu_loop_length = 9; 
    int min_tu_loop_length = 6; 
    int min_aa_to_ac_length = 19;
    int max_aa_to_ac_length = 35;
    int min_acs_to_ace_length = 16;
    int max_var_loop_length = 28;
    int min_aa_score = 12;    /* this needs to be 9 for 1 chmpxx awkward customer */
    int min_ac_score = 8; 
    int min_tu_score = 9; 
    int min_d_score = 4;
    int min_total_bp_score = 36; 
    int min_total_cb_score = 16;


#endif
    int con_base_scores[18] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    int ret;
    int i;
    
    *t = init_TrnaSpec (max_trna_length,
		       min_trna_length,
		       max_intron_length,
		       min_intron_length,
		       max_tu_loop_length,
		       min_tu_loop_length,
		       min_aa_to_ac_length,
		       max_aa_to_ac_length,
		       min_acs_to_ace_length,
		       max_var_loop_length,
		       min_aa_score,
		       min_ac_score,
		       min_tu_score,
		       min_d_score,
		       min_total_bp_score,
		       min_total_cb_score,
		       con_base_scores);

    for (i = 0; i < MAX_TRNA; i++) {
	if (( (*results)[i] = init_TrnaRes ()) == NULL) 
	    return -2;
    }
    ret = do_trna_search ( seq, seq_length, user_start, user_end, *t, 
			   results, nmatch, max_total_bp_score);
    
    /* free */
    return ret;
}
