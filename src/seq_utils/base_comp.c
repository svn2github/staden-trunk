#include <stdio.h>
#include <float.h>
#include <string.h>

#include "dna_utils.h"
#include "misc.h"
#include "base_comp.h"
#include "edge.h"

int
Plot_Base_Comp(int win_len, 
	       int *score,
	       char *seq,
	       int seq_len,
	       int user_start,
	       int user_end,
	       int *match,
	       int *max_match)
{
    int i;
    int cnt = 0;
    int front, back;
    int half_span;
    int pos;

    (*max_match) = -1;
    half_span = win_len / 2;
    /*
     * do the left hand edge
     */
    for (i = 0, pos = -half_span; i < win_len; i++, pos++) {
	cnt += score[char_lookup[(unsigned int)seq[i]]];
	if (pos >= 0) 
	    match[pos] = cnt;
	if (cnt > (*max_match))
	    (*max_match) = cnt;
/*
	printf("1 i %d i %d cnt %d\n", i, i-half_span, match[i-half_span]); 
*/
    }
    /*
     * do the rest by subtracting the score from the back edge and adding the
     * score of the front edge of the sliding window
     */
    for (i = win_len; i < seq_len; i++, pos++) {
	front = i;
	back = i - win_len;
	cnt += score[char_lookup[(unsigned int)seq[front]]] -
	    score[char_lookup[(unsigned int)seq[back]]];
	match[pos] = cnt;

	/*printf("2f %d b %d pos %d cnt %d %c\n", front, back, pos, match[pos], seq[back]);*/
	if (cnt > (*max_match))
	    (*max_match) = cnt;
    }
    for (i = seq_len - win_len; i < seq_len - half_span; i++, pos++) {
	cnt -= score[char_lookup[(unsigned int)seq[i]]];
	match[pos] = cnt;
	/*printf("3f i %d pos %d cnt %d %c\n", i, pos, match[pos], seq[i]); */
    }

    return 0;
}


int 
get_base_comp_res(char seq[], int seq_length, 
		  int window_length, /* sum over all windows this length */
		  int user_start,    /* seq start numbering from 1 */
		  int user_end,      /* seq end numbering from 1 */
		  double score[],    /* score for each base type */
		  double result[],   /* put results here */
		  double *min,       /* min result */
		  double *max)       /* max result */
{

    int edge_length, middle, front_edge, back_edge, start, end, i;
    char *edge;
    *max = -1;
    *min = DBL_MAX;

    /* note well: results always start at results[0] */

    /* quick sanity checks */

    if (  ( window_length % 2 ) == 0 ) return -1;
    if ( user_start < 1 ) return -1;
    if ( user_end > seq_length ) return -1;
    if ( window_length > (user_end - user_start + 1) ) return -1;

    start = user_start - 1;
    end = user_end - 1;

    /* left end game */


    edge = seq_left_end ( seq, seq_length, start, window_length, 1);
    if ( ! edge ) return -1;

    edge_length = strlen ( edge );

    for ( i=0, result[0] = 0.0; i<window_length; i++ ) {
	result[0] += score [ char_lookup [ (unsigned int)edge [ i ] ] ];
    }

    if (result[0] > (*max)) (*max) = result[0];
    if (result[0] < (*min)) (*min) = result[0];

    for ( 
	  back_edge = 0,
	  front_edge = window_length,
	  middle = 1;
	  front_edge < edge_length;
          front_edge++,
	  back_edge++,
	  middle++) {

	  result[middle] = result[middle-1]
	    - score [ char_lookup [ (unsigned int)edge [ back_edge ] ] ]
	    + score [ char_lookup [ (unsigned int)edge [ front_edge ] ] ];
	  if (result[middle] > (*max)) (*max) = result[middle];
	  if (result[middle] < (*min)) (*min) = result[middle];
      }

    /* middle game */

    for ( 
	  back_edge = start,
	  front_edge = back_edge + window_length;
	  front_edge <= end;
          front_edge++,
	  back_edge++,
	  middle++) {

	  result[middle] = result[middle-1]
	    - score [ char_lookup [ (unsigned int)seq [ back_edge ] ] ]
	    + score [ char_lookup [ (unsigned int)seq [ front_edge ] ] ];
	  if (result[middle] > (*max)) (*max) = result[middle];
	  if (result[middle] < (*min)) (*min) = result[middle];
      }
    xfree ( edge );

    /* right end game */

    edge = seq_right_end ( seq, seq_length, end, window_length, 1);
    if ( ! (edge) ) return -1;

    edge_length = strlen ( edge );

    for ( 
	  back_edge = 0,
	  front_edge = back_edge + window_length;
	  front_edge < edge_length;
          front_edge++,
	  back_edge++,
	  middle++) {

	  result[middle] = result[middle-1]
	    - score [ char_lookup [ (unsigned int)edge [ back_edge ] ] ]
	    + score [ char_lookup [ (unsigned int)edge [ front_edge ] ] ];

	  if (result[middle] > (*max)) (*max) = result[middle];
	  if (result[middle] < (*min)) (*min) = result[middle];
      }
    xfree ( edge );

    return 0;
}

double get_base_comp_mass(int a, 
			  int c, 
			  int g,
			  int t)
{
    double mass;
    double t_mass = 304.19618;
    double c_mass = 289.18454;
    double a_mass = 313.20954;
    double g_mass = 329.20894;
    double hydrogen_mass = 1.00794;
    double phosphate_mass = 62.97256;

    mass = t * t_mass + c * c_mass + a *a_mass + g * g_mass + hydrogen_mass
	- phosphate_mass;
    return mass;
}


void get_aa_comp ( char *seq, int seq_length, double aa_comp[25] ) {
    int i;

    for ( i = 0; i < 25; i++ ) {
	aa_comp[i] = 0.0;
    }
    for ( i = 0; i < seq_length; i++ ) {
	aa_comp[protein_lookup[(unsigned int)seq[i]]] += 1.0;
    }
}

void get_aa_comp_mass(double *aa_comp, double *mass)
{
    double aa_mass[25];
    int i;

    aa_mass[0] = 71.0788;      /* A */
    aa_mass[1] = 0.0;          /* B */
    aa_mass[2] = 103.1448;     /* C */
    aa_mass[3] = 115.0886;     /* D */
    aa_mass[4] = 129.1155;     /* E */
    aa_mass[5] = 147.1766;     /* F */
    aa_mass[6] = 57.0520;      /* G */
    aa_mass[7] = 137.1412;     /* H */
    aa_mass[8] = 113.1595;     /* I */
    aa_mass[9] = 128.1742;     /* K */
    aa_mass[10] = 113.1595;    /* L */
    aa_mass[11] = 131.1986;    /* M */
    aa_mass[12] = 114.1039;    /* N */
    aa_mass[13] = 97.1167;     /* P */
    aa_mass[14] = 128.1308;    /* Q */
    aa_mass[15] = 156.1876;    /* R */
    aa_mass[16] = 87.0782;     /* S */
    aa_mass[17] = 101.1051;    /* T */
    aa_mass[18] = 99.1326;     /* V */
    aa_mass[19] = 186.2133;    /* W */
    aa_mass[20] = 163.1760;    /* Y */
    aa_mass[21] = 0.0;         /* Z */
    aa_mass[22] = 0.0;         /* X */
    aa_mass[23] = 0.0;         /* * */
    aa_mass[24] = 0.0;         /* - */

#ifdef REMOVE
    double a_mass = 71.0788;
    double b_mass = 0.0;
    double c_mass = 103.1448;
    double d_mass = 115.0886;
    double e_mass = 129.1155;
    double f_mass = 147.1766;
    double g_mass = 57.0520;
    double h_mass = 137.1412;
    double i_mass = 113.1595;
    double k_mass = 128.1742;
    double l_mass = 113.1595;
    double m_mass = 131.1986;
    double n_mass = 114.1039;
    double p_mass = 97.1167;
    double q_mass = 128.1308;
    double r_mass = 156.1876;
    double s_mass = 87.0782;
    double t_mass = 101.1051;
    double v_mass = 99.1326;
    double w_mass = 186.2133;
    double y_mass = 163.1760;
    double z_mass = 0.0;
    double x_mass = 0.0;
    double *_mass = 0.0;
    double -_mass = 0.0;
#endif

    for (i = 0; i < 25; i++) {
	mass[i] = 0.0;
    }

    for (i = 0; i < 25; i++) {
	mass[i] += aa_mass[i] * aa_comp[i];
    }
}
