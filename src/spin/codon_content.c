#include <staden_config.h>

#include <string.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "misc.h"
#include "dna_utils.h"
#include "genetic_code.h"
#include "edge.h"
#include "codon_content.h"
#include "array_arith.h"
#include "sequence_formats.h"

/* gene search by content */

CodRes *init_CodRes ( int num_results ) {

    CodRes *r;
    double *f0, *f1, *f2;
    char *top;

    if ( num_results < 1 ) return NULL;
    if ( ( NULL == ( r = ( CodRes* ) xmalloc ( sizeof (CodRes) )))) return NULL;
    if ( ( NULL == ( f0 = ( double* ) xmalloc ( sizeof (double) * num_results )))) return NULL;
    if ( ( NULL == ( f1 = ( double* ) xmalloc ( sizeof (double) * num_results )))) return NULL;
    if ( ( NULL == ( f2 = ( double* ) xmalloc ( sizeof (double) * num_results )))) return NULL ;
    if ( ( NULL == ( top = ( char* ) xmalloc ( sizeof (char) * num_results+1 )))) return NULL ;
    r->frame1 = f0;
    r->frame2 = f1;
    r->frame3 = f2;
    r->top    = top;
    r->num_results = num_results;
    r->user_start = 0;
    r->user_end = 0;
    r->max = 0.0;
    r->min = 0.0;
    r->window_length = 0;
    return r;
}

void free_CodRes ( CodRes *r ) {

    xfree (r->frame1);
    xfree (r->frame2);
    xfree (r->frame3);
    xfree (r->top);
    xfree (r);
}


void get_tops ( CodRes *results ) {

    int i, j, t;
    for(i=0,j=0,t=0;j<results->num_results;i+=3,j++) {
	if ( results->frame1[j] >= results->frame2[j] ) {
	    if ( results->frame3[j] >= results->frame1[j] ) {
		t = 3;
	    }
	    else {
		t = 1;
	    }
	}
	else {
	    if ( results->frame3[j] >= results->frame2[j] ) {
		t = 3;
	    }
	    else {
		t = 2;
	    }
	}
	results->top[j] = t;
    }
    results->top[results->num_results] = '\0';
}

void comp_from_cods ( double base_comp[], double codon_table[4][4][4] ) {

  /* make the base composition from a codon table */
  int i, j, k;
  double total;

  j = 5;
  for (i=0;i<j;i++) {
    base_comp[i] =  0.0;
  }
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      for (k=0;k<4;k++) {
	base_comp[i] += codon_table[i][j][k];
	base_comp[j] += codon_table[i][j][k];
	base_comp[k] += codon_table[i][j][k];
      }
    }
  }

  j = 5;
  for (i=0,total=0.0;i<j;i++) {
    total += base_comp[i];
  }
  if ( total > DBL_EPSILON ) {
    for (i=0;i<j;i++) {
      base_comp[i] /= total;
    }
  }
}

int set_stops_zeroes ( double freqs_64[64] ) {

  int i,j,k,i64,num_not_stop;
  double mean, total;
  char (*genetic_code)[5][5] = get_global_genetic_code();

  for ( i=0,i64=0,total=0.0,num_not_stop=0; i<4; i++ ) {
    for ( j=0; j<4; j++ ) {
      for ( k=0; k<4; k++, i64++ ) {
	if ( genetic_code[i][j][k] == '*' ) {
	  freqs_64[i64] = -1.0;
	}
	else {
	  num_not_stop++;
	  total+= freqs_64[i64];
	}
      }
      
    }
  }

  /* total codon count for non-stop codons is total */

  if ( 0.0 == total ) return -1;
  if ( 0 == num_not_stop ) return -2;

  mean = total/(double)num_not_stop;

  /* stop codon elements must be set to mean value for table */

  for (i64=0; i64<64; i64++ ) {
    if ( freqs_64[i64] < 0.0 ) freqs_64[i64] = mean;
  }

  /* zero elements must be set to 1/k where k is total codon count */

  for (i64=0; i64<64; i64++ ) {
    if ( freqs_64[i64] == 0.0 ) {
      freqs_64[i64] = 1.0/total;
    }
  }

  return 0;
}

int get_codon_scores ( char seq[], int seq_length, 
		       int window_length, /* sum over all windows this length */
		       int user_start,    /* seq start numbering from 1 */
		       int user_end,      /* seq end numbering from 1 */
		       double codon_table[4][4][4],    /* score for each codon type */
		       double result[],   /* put results here */
		       int num_results )  /* size of results array */
{
    int edge_length, middle, front_edge, back_edge, start, end, inc=3, i;
    char *edge;
    double freqs_64[64];
    double mean;
    int *genetic_code_idx = get_genetic_code_idx(0);

    /* note well: results always start at results[0] */

    /* quick sanity checks */

    if (  ( window_length % 2 ) == 0 ) return -1;
    if (  ( window_length % 3 ) != 0 ) return -1;
    if ( user_start < 1 ) return -1;
    if ( user_end > seq_length ) return -1;
    if ( window_length > (user_end - user_start + 1) ) return -1;
    set_char_set(DNA); 

    start = user_start - 1;
    end = user_end - 1;

    /* make sure we end on a codon boundary */

    i = ( end - start + 1 ) / 3;
    end = start - 1 + 3 * i;

    /****************** left end game **********************/

    codon_table_64 ( codon_table, freqs_64, TO_64 );
     mean = sum_double_array ( freqs_64, 64 ) / 64.0;
     /*     printf("mean %f\n",mean);*/
    edge = seq_left_end ( seq, seq_length, start, window_length, inc);
    if ( ! edge ) return -1;

    edge_length = strlen ( edge );

    /* printf("mean=%f\n",mean);*/

    for ( i=0, result[0] = mean; i<window_length; i+=inc ) {
      /*    printf ("legal=%d\n",legal_codon ( &edge [ i ] ));*/
	if ( legal_codon ( &edge [ i ] ) ) {
	    result[0] += 
		codon_table 
		    [ genetic_code_idx [ char_lookup [ edge [ i   ] ] ]]
		    [ genetic_code_idx [ char_lookup [ edge [ i+1 ] ] ]]
		    [ genetic_code_idx [ char_lookup [ edge [ i+2 ] ] ]];
	    /*	 printf("result1[0]=%f\n",result[0]);*/   
	}
	else {
	  result[0] += mean;
	  /*	 printf("result0[0]=%f\n",result[0]);*/ 
	}
    }

    for ( 
	  back_edge = 0,
	  front_edge = window_length,
	  middle = 1;
	  front_edge < edge_length;
          front_edge+=inc,
	  back_edge+=inc,
	  middle++) {

	result [ middle ] = result [ middle - 1 ];
	if ( legal_codon ( &edge [ back_edge ] ) ) {
	    result[middle] -= 
		codon_table 
		    [ genetic_code_idx [ char_lookup [ edge [ back_edge   ] ] ]]
		    [ genetic_code_idx [ char_lookup [ edge [ back_edge+1 ] ] ]]
		    [ genetic_code_idx [ char_lookup [ edge [ back_edge+2 ] ] ]];
	}
	else {
	  result[middle] -= mean;
	}
	if ( legal_codon ( &edge [ front_edge ] ) ) {
	    result[middle] += 
		codon_table 
		    [ genetic_code_idx [ char_lookup [ edge [ front_edge   ] ] ]]
		    [ genetic_code_idx [ char_lookup [ edge [ front_edge+1 ] ] ]]
		    [ genetic_code_idx [ char_lookup [ edge [ front_edge+2 ] ] ]];
	  }
	else {
	  result[middle] += mean;
	}
      }

    /****************** middle game **********************/

    for ( 
	  back_edge = start,
	  front_edge = back_edge + window_length;
	  front_edge <= end;
          front_edge+=inc,
	  back_edge+=inc,
	  middle++) {

	result [ middle ] = result [ middle - 1 ];
	if ( legal_codon ( &seq [ back_edge ] ) ) {
	    result[middle] -= 
		codon_table 
		    [ genetic_code_idx [ char_lookup [ seq [ back_edge   ] ] ]]
		    [ genetic_code_idx [ char_lookup [ seq [ back_edge+1 ] ] ]]
		    [ genetic_code_idx [ char_lookup [ seq [ back_edge+2 ] ] ]];
	}
	else {
	  result[middle] -= mean;
	}
	if ( legal_codon ( &seq [ front_edge ] ) ) {
	    result[middle] += 
		codon_table 
		    [ genetic_code_idx [ char_lookup [ seq [ front_edge   ] ] ]]
		    [ genetic_code_idx [ char_lookup [ seq [ front_edge+1 ] ] ]]
		    [ genetic_code_idx [ char_lookup [ seq [ front_edge+2 ] ] ]];
	}
	else {
	  result[middle] += mean;
	}
      }
    free ( edge );

    /****************** right end game **********************/

    edge = seq_right_end ( seq, seq_length, end, window_length, inc);
    if ( ! (edge) ) return -1;

    edge_length = strlen ( edge );

    for ( 
	  back_edge = 0,
	  front_edge = back_edge + window_length;
	  front_edge < edge_length;
          front_edge+=inc,
	  back_edge+=inc,
	  middle++) {

	result [ middle ] = result [ middle - 1 ];
	if ( legal_codon ( &edge [ back_edge ] ) ) {
	    result[middle] -= 
		codon_table 
		    [ genetic_code_idx [ char_lookup [ edge [ back_edge   ] ] ]]
		    [ genetic_code_idx [ char_lookup [ edge [ back_edge+1 ] ] ]]
		    [ genetic_code_idx [ char_lookup [ edge [ back_edge+2 ] ] ]];
	}
	else {
	  result[middle] -= mean;
	}
	if ( legal_codon ( &edge [ front_edge ] ) ) {
	    result[middle] += 
		codon_table 
		    [ genetic_code_idx [ char_lookup [ edge [ front_edge   ] ] ]]
		    [ genetic_code_idx [ char_lookup [ edge [ front_edge+1 ] ] ]]
		    [ genetic_code_idx [ char_lookup [ edge [ front_edge+2 ] ] ]];
	}
	else {
	  result[middle] += mean;
	}
      }
    free ( edge );
    /* for some sizes of active region we do not calculate a value for
       the last element. Here we set the last value in the array to the 
       value of the last element calculated so we are sure it is sensible */

    result[ num_results - 1 ] = result[ middle - 1 ];
    return 0;
}

/* no longer used */
#if 0
int do_pos_base_pref ( char seq[], int seq_length,
		      CodRes *results ) {

    int i, j, res;
    double m, m0, m1, m2;
    double codon_usage_table[4][4][4];

    char star[3];
    char file_name[1024];

/* before we get here we need to alloc CodRes and do 
   the following 3 initialisations:
    set_dna_lookup();
    set_char_set(DNA); 
    init_genetic_code();
*/

    file_name[0] = '\0';
    i= init_codon_pref (file_name, codon_usage_table, 0 );

    res = write_cod_table (stdout, codon_usage_table);

    res = get_codon_scores ( seq, seq_length, results->window_length, 
				results->user_start, results->user_end,
			        codon_usage_table, results->frame1,
				results->num_results );
    if ( res ) return -1;

    res = get_codon_scores ( seq, seq_length, results->window_length, 
				results->user_start+1, results->user_end,
			        codon_usage_table, results->frame2,
				results->num_results );
    if ( res ) return -1;

    res = get_codon_scores ( seq, seq_length, results->window_length, 
				results->user_start+2, results->user_end,
			        codon_usage_table, results->frame3,
				results->num_results );
    if ( res ) return -1;

    m0 = max_double_array ( results->frame1, results->num_results);
    m1 = max_double_array ( results->frame2, results->num_results);
    m2 = max_double_array ( results->frame3, results->num_results);
    m = MAX ( m0, m1 );
    m = MAX ( m, m2 );

    results->max = m;
    results->min = -m;

    /* printf("max and min %f %f\n",results->max, results->min); */
    get_tops ( results );

    for(i=0,j=0;j<results->num_results;i+=3,j++) {
	star[0] = star[1] = star[2] = ' ';
	star[results->top[j]-1] = '*';
	/*
	printf("i %d %f%c %f%c %f%c\n",i+1,results->frame1[j],star[0],results->frame2[j],star[1],results->frame3[j],star[2]);
	*/
    }

    return 0;
  }
#endif

int do_codon_pref (char seq[], int seq_length,
		   double codon_usage_table[4][4][4],
		   CodRes *results ) {

    int i, j, res;
    double m, m0, m1, m2;
    double mm, mm0, mm1, mm2;
    char star[3];

/* before we get here we need to alloc CodRes, get the codon table and do 
   the following 3 initialisations:
    set_dna_lookup();
    set_char_set(DNA); 
    init_genetic_code();
*/

    res = get_codon_scores ( seq, seq_length, results->window_length, 
			     results->user_start, results->user_end,
			     codon_usage_table, results->frame1,
			     results->num_results );
    if ( res ) return -1;

    res = get_codon_scores ( seq, seq_length, results->window_length, 
			     results->user_start+1, results->user_end,
			     codon_usage_table, results->frame2,
			     results->num_results );
    if ( res ) return -1;

    res = get_codon_scores ( seq, seq_length, results->window_length, 
			     results->user_start+2, results->user_end,
			     codon_usage_table, results->frame3,
			     results->num_results );
    if ( res ) return -1;

    m0 = max_double_array ( results->frame1, results->num_results);
    m1 = max_double_array ( results->frame2, results->num_results);
    m2 = max_double_array ( results->frame3, results->num_results);
    mm0 = min_double_array ( results->frame1, results->num_results);
    mm1 = min_double_array ( results->frame2, results->num_results);
    mm2 = min_double_array ( results->frame3, results->num_results);

    m = MAX ( m0, m1 );
    m = MAX ( m, m2 );
    mm = MIN(mm0, mm1);
    mm = MIN(mm, mm2);

    m = MAX(m, fabs(mm));

    results->max = m;
    results->min = -m;


    /* printf("max and min %f %f\n",results->max, results->min); */

    get_tops ( results );

    for(i=0,j=0;j<results->num_results;i+=3,j++) {
	star[0] = star[1] = star[2] = ' ';
	star[results->top[j]-1] = '*';
	
	/* printf("i %d %f%c %f%c %f%c\n",i+1,results->frame1[j],star[0],results->frame2[j],star[1],results->frame3[j],star[2]); */
    }

    return 0;
}

int do_author_test ( char seq[], int seq_length,
		      double codon_usage_table[4][4][4],
		      CodRes *results ) {

    int i, j, res;
    double m, m0, m1, m2;
    char star[3];

/* before we get here we need to alloc CodRes, get the codon table and do 
   the following 3 initialisations:
    set_dna_lookup();
    set_char_set(DNA); 
    init_genetic_code();
*/


    res = get_codon_scores ( seq, seq_length, results->window_length, 
				results->user_start, results->user_end,
			        codon_usage_table, results->frame1,
				results->num_results );
    if ( res ) return -1;

    res = get_codon_scores ( seq, seq_length, results->window_length, 
				results->user_start+1, results->user_end,
			        codon_usage_table, results->frame2,
				results->num_results );
    if ( res ) return -1;

    res = get_codon_scores ( seq, seq_length, results->window_length, 
				results->user_start+2, results->user_end,
			        codon_usage_table, results->frame3,
				results->num_results );
    if ( res ) return -1;

    m0 = max_double_array ( results->frame1, results->num_results);
    m1 = max_double_array ( results->frame2, results->num_results);
    m2 = max_double_array ( results->frame3, results->num_results);
    m = MAX ( m0, m1 );
    m = MAX ( m, m2 );

    results->max = m;
    m0 = min_double_array ( results->frame1, results->num_results);
    m1 = min_double_array ( results->frame2, results->num_results);
    m2 = min_double_array ( results->frame3, results->num_results);
    m = MIN ( m0, m1 );
    m = MIN ( m, m2 );

    results->min = m;

    /* printf("max and min %f %f\n",results->max, results->min); */

    get_tops ( results );

    for(i=0,j=0;j<results->num_results;i+=3,j++) {
	star[0] = star[1] = star[2] = ' ';
	star[results->top[j]-1] = '*';
	/* printf("i %d %f%c %f%c %f%c\n",i+1,results->frame1[j],star[0],results->frame2[j],star[1],results->frame3[j],star[2]); */
    }

    return 0;
}

int init_codon_pref (char *file_name, 
		     double codon_usage_table[4][4][4],
		     int option) {
    
    FILE *in_file;
    int res;
    double coding_freqs[4][4][4], non_coding_freqs[4][4][4];
    double coding_freqs_64[64], non_coding_freqs_64[64];
    double base_comp[5], total;
    int coding, non_coding;
    int i;
#define average_aa 2
#define acids_only 4
#define third_pos_only 8

    coding = non_coding = 0;
    if ( strcmp(file_name, "") != 0 ) {
      coding = 1;
      in_file = fopen(file_name,"r");
      if ( ! in_file ) {
	coding = 0;
      }

      if ( coding ) {
	res = read_cod_table ( in_file, coding_freqs );
	if ( !res ) {
	  fclose ( in_file );
	  coding = 0;
	}
      }

      non_coding = 0;
      if ( coding ) {
	res = read_cod_table ( in_file, non_coding_freqs );
	non_coding = 1;
	if ( !res ) {
	  non_coding = 0;
	}
	fclose ( in_file );
      }
      
      if ( coding ) {
	res = write_screen_cod_table(coding_freqs);
      }
      if ( non_coding ) {
	res = write_screen_cod_table(non_coding_freqs);
      }
    }

    if ( !coding ) {

      /* no codon table, so only option is original pos base pref */

      gen_cods_from_ac ( coding_freqs );
      /*res = write_cod_table (stdout, coding_freqs );*/

      /* get the base composition from the input coding table */
      
      comp_from_cods ( base_comp, coding_freqs );
      gen_cods_from_bc(non_coding_freqs, base_comp);
      codon_table_64 ( non_coding_freqs, non_coding_freqs_64, TO_64 );

      /* this code is identical to that further below for input tables*/

      /* deal with stop codons and zero values. Note zeroes set to 1/total
       * so need say total=1000 before this is done.
       */

      codon_table_64 ( coding_freqs, coding_freqs_64, TO_64 );
      scale_double_array ( coding_freqs_64, 64, 1000.0);
      
      if ( (i = set_stops_zeroes ( coding_freqs_64 )) < 0 ) return i;
      codon_table_64 ( coding_freqs, coding_freqs_64, FROM_64);
      /*res = write_cod_table (stdout, coding_freqs);*/
      total = sum_double_array ( coding_freqs_64, 64 );
      div_double_array ( coding_freqs_64, 64, total );
      
      total = sum_double_array ( non_coding_freqs_64, 64 );
      div_double_array ( non_coding_freqs_64, 64, total );
      ratio_double_arrays ( coding_freqs_64, non_coding_freqs_64, 64);
      log_double_array ( coding_freqs_64, 64 );
      scale_double_array1 ( coding_freqs_64, 64, 0.0);
      codon_table_64 ( codon_usage_table, coding_freqs_64, FROM_64 );
      /* printf("pbp\n"); */
      return 0;
    }

    if ( !non_coding ) {

      /* get the base composition from the input coding table */
      
      comp_from_cods ( base_comp, coding_freqs );
      gen_cods_from_bc(non_coding_freqs, base_comp);
    }

    if ( option & average_aa ) {

      /* printf("av aa\n"); */
      average_acid_comp ( coding_freqs );

    }

    if ( option & acids_only ) {

      /* printf("aa only\n"); */
      even_cods_per_acid ( coding_freqs );

    }


    if ( option & third_pos_only ) {

      /* printf("3 only\n"); */
      third_pos_comp( coding_freqs);
      res = write_cod_table (stdout, coding_freqs);

    }


    /* deal with stop codons and zero values. Note zeroes set to 1/total
     * so need say total=1000 before this is done.
     */

    codon_table_64 ( coding_freqs, coding_freqs_64, TO_64 );
    scale_double_array ( coding_freqs_64, 64, 1000.0);
    
    if ( (i = set_stops_zeroes ( coding_freqs_64 )) < 0 ) return i;
    /*codon_table_64 ( coding_freqs, coding_freqs_64, FROM_64);
    res = write_cod_table (stdout, coding_freqs);*/
    total = sum_double_array ( coding_freqs_64, 64 );
    div_double_array ( coding_freqs_64, 64, total );
    
    /*res = write_cod_table (stdout, non_coding_freqs);*/

    /* treat non-coding table the same */

    codon_table_64 ( non_coding_freqs, non_coding_freqs_64, TO_64 );
    scale_double_array ( non_coding_freqs_64, 64, 1000.0);
    if ( (i = set_stops_zeroes ( non_coding_freqs_64 )) < 0 ) return i;
    total = sum_double_array ( non_coding_freqs_64, 64 );
    div_double_array ( non_coding_freqs_64, 64, total );

    /* divide coding by non-coding, take logs and rescale */

    ratio_double_arrays ( coding_freqs_64, non_coding_freqs_64, 64);
    log_double_array ( coding_freqs_64, 64 );
    scale_double_array1 ( coding_freqs_64, 64, 0.0);
    codon_table_64 ( codon_usage_table, coding_freqs_64, FROM_64 );
    return 0;
  }


void get_author_weights ( double array1[], double array2[], 
			 double weights[], int num_elements ) {
    int i;
    double small = DBL_EPSILON;

    for ( i=0; i<num_elements; i++ ) {
	if ( array2[i] > small ) {
	    weights[i] = log ( array1[i] / array2[i] );
	}
	else {
	    weights[i] = 0.0;
	}
    }
}
double author_weighted_mean ( double array1[], double array2[], 
			  int num_elements ) {
    int i;
    double t=0.0;

    for ( i=0; i<num_elements; i++ ) {
	t += array1[i] * array2[i];
    }
    return t;
}

double author_variance ( double array1[], double array2[], 
			  int num_elements ) {
    int i;
    double t1=0.0,t2=0.0,t;
    double small = DBL_EPSILON;

    for ( i=0; i<num_elements; i++ ) {
	t = array1[i] * array2[i];
	t1 += t * array2[i];
	t2 += t;
    }
    t = t1 - t2 * t2;
    if ( t > small ) {
	return sqrt ( t );
    }
    else {
	return 0.0;
    }
}

double normal_x ( double prob ) {
  /* normal_x is that which is exceeded by a random variable,
   * normally distributed with zero mean and unit variance,
   * with probability prob/100.
   */
  int i, size;
  struct percent_normal {
    double p;
    double x;
  };
  struct percent_normal values[] =
    
  { 

    {40.0, 0.2533},
    {30.0, 0.5244},
    {20.0, 0.8416},
    {10.0, 1.2816},
    { 5.0, 1.6449},
    { 4.0, 1.7507},
    { 3.0, 1.8808},
    { 2.0, 2.0537},
    { 1.0, 2.3263},
    { 0.9, 2.3656},
    { 0.8, 2.4089},
    { 0.7, 2.4573},
    { 0.6, 2.5121},
    { 0.5, 2.5758},
    { 0.4, 2.6521},
    { 0.3, 2.7478},
    { 0.2, 2.8782},
    { 0.1, 3.0902},
    { 0.05, 3.2905},
    { 0.01, 3.7190},
    { 0.005, 3.8906},
    { 0.001, 4.2649},
    { 0.0005, 4.4172},
    { 0.0001, 4.753},
    { 0.00005, 4.892}
  };

  size = sizeof(values)/sizeof (struct percent_normal);
  for (i=0;i<size;i++) {
    if ( prob >= values[i].p ) return values[i].x;
  }
  return values[size-1].x;
}


       
int init_author_test (char *file_name, 
		      char *seq,
		      int seq_length,
		      double codon_usage_wm[4][4][4], 
		      double percent_error,
		      int *window_length ) 
{
  /* FIXME: need new parameter sending - percent_error
     which should lie in the range 0.0005 to 20
     */
    FILE *in_file;
    int res;
    double coding_table[4][4][4], noncoding_table[4][4][4];
    double coding_table64[64], noncoding_table64[64];
    double weights[64];
    double sum_coding, sum_noncoding, coding_wm, noncoding_wm;
    double coding_sd, noncoding_sd, fact, discrim_value, window;
    double base_comp[5];
    double kappa;

    in_file = fopen(file_name,"r");
    if ( ! in_file ) {
	return -1;
    }

    res = read_cod_table ( in_file, coding_table );

    if ( !res ) {
	return -1;
    }

    res = write_screen_cod_table(coding_table);
    if ( !res ) {
	fclose ( in_file );
	return -1;
    }

    res = read_cod_table ( in_file, noncoding_table );
    if ( !res ) {
	get_base_comp(seq, seq_length, base_comp);
	scale_double_array(base_comp, 4, 1.0);
	gen_cods_from_bc(noncoding_table, base_comp);
	scale_codon_table(noncoding_table, 1000.);
    }

    fclose ( in_file );
    res = write_screen_cod_table(noncoding_table);
    if ( !res ) {
	return -1;
    }

    codon_table_64 ( coding_table, coding_table64, TO_64 );
    codon_table_64 ( noncoding_table, noncoding_table64, TO_64 );

    reset_zeroes ( coding_table64, 64, 0.001 );
    reset_zeroes ( noncoding_table64, 64, 0.001 );

    sum_coding = sum_double_array ( coding_table64, 64 );
    sum_noncoding = sum_double_array ( noncoding_table64, 64 );

    div_double_array ( coding_table64, 64, sum_coding );
    div_double_array ( noncoding_table64, 64, sum_noncoding );

    get_author_weights ( coding_table64, noncoding_table64, weights, 64 );

    coding_wm = author_weighted_mean ( coding_table64, weights, 64 );

    noncoding_wm = author_weighted_mean ( noncoding_table64, weights, 64 );

    coding_sd = author_variance ( coding_table64, weights, 64 );
    noncoding_sd = author_variance ( noncoding_table64, weights, 64 );

    /* given the above values for coding_wm, noncoding_wm, coding_sd and
       noncoding_sd, when classifying a region into coding or noncoding,
       the probability of an error and the window_length are related by
       fact = sqrt ( window_length ) * S where S is the error size in sd's.
       The discriminating value discrim_value is 
       noncoding_wm + (fact * sqrt (window_length ) * noncoding_sd
    */

    fact = ( coding_wm - noncoding_wm ) / ( coding_sd + noncoding_sd );
    discrim_value = noncoding_wm + fact * noncoding_sd;

    /* assume an error of 1/1000 then kappa = 3.09
     *                    1/100               2.33
     *                    1/10                1.29
     */

    kappa = normal_x ( percent_error );
    window = ( kappa / fact ) * ( kappa / fact );

    *window_length = window;
    /* printf("window_length %d\n",*window_length ); */
    *window_length = 3 * ((2*(*window_length/2))+1);
    div_double_array ( weights, 64, *window_length );
    codon_table_64 ( codon_usage_wm, weights, FROM_64 );
    /* printf("window_length %d\n",*window_length ); */

    return 0;
}

CodRes1 *init_CodRes1 ( int num_results ) {

    CodRes1 *r;
    double *f0;

    if ( num_results < 1 ) return NULL;
    if ( ( NULL == ( r = ( CodRes1* ) xmalloc ( sizeof (CodRes1) )))) return NULL;
    if ( ( NULL == ( f0 = ( double* ) xmalloc ( sizeof (double) * num_results )))) return NULL;

    r->frame1 = f0;
    r->num_results = num_results;
    r->user_start = 0;
    r->user_end = 0;
    r->max = 0.0;
    r->min = 0.0;
    r->window_length = 0;
    return r;
}

void free_CodRes1 ( CodRes1 *r ) {

    xfree (r->frame1);
    xfree (r);
}

int get_pos_base_bias ( char seq[], int seq_length,
		       int user_start,
		       int user_end,
		       double result[],
		       int num_results,
		       int window_length )

{
    int edge_length, middle, front_edge, back_edge, start, end, inc=3, i;
    char *edge;
    double base_comp[5], pos_base_comp[15], pos_base_comp_bias[15];
    double small = DBL_EPSILON;
    int base_type_pos1,base_type_pos2,base_type_pos3;
    double sum_bias, t_comp;

    /* note well: results always start at results[0] */

    /* quick sanity checks */

    if (  ( window_length % 2 ) == 0 ) {
	return -1;
    }
    if (  ( window_length % 3 ) != 0 ){
	return -1;
    }
    if ( user_start < 1 ) return -1;
    if ( user_end > seq_length ) return -1;
    if ( window_length > (user_end - user_start + 1)){
	return -1;
    }
    start = user_start - 1;
    end = user_end - 1;

    /* make sure we end on a codon boundary */

    i = ( end - start + 1 ) / 3;
    end = start - 1 + 3 * i;

    /****************** left end game **********************/


    edge = seq_left_end ( seq, seq_length, start, window_length, inc);
    if ( ! edge ) return -1;

    edge_length = strlen ( edge );

    fill_double_array ( base_comp, 5, 0.0 );
    fill_double_array ( pos_base_comp, 15, 0.0 );
    fill_double_array ( pos_base_comp_bias, 15, 0.0 );

    for ( i=0; i<window_length; i+=inc ) {

	base_type_pos1 = char_lookup [ edge [ i ] ];
	base_type_pos2 = char_lookup [ edge [ i+1 ] ];
	base_type_pos3 = char_lookup [ edge [ i+2 ] ];

	base_comp [ base_type_pos1 ] += 1.0;
	base_comp [ base_type_pos2 ] += 1.0;
	base_comp [ base_type_pos3 ] += 1.0;
	pos_base_comp [ base_type_pos1 ] += 1.0;
	pos_base_comp [ base_type_pos2 + 5 ] += 1.0;
	pos_base_comp [ base_type_pos3 + 10 ] += 1.0;
    }

    for ( i=0; i<4; i++ ) {
	t_comp = base_comp [ i ] / 3.0;
	pos_base_comp_bias [ i ] = 
	    fabs ( pos_base_comp [ i ] - t_comp );
	pos_base_comp_bias [ i + 5 ] = 
	    fabs ( pos_base_comp [ i + 5 ] - t_comp );
	pos_base_comp_bias [ i + 10 ] = 
	    fabs ( pos_base_comp [ i + 10 ] - t_comp );
    }

    for ( i=0, sum_bias = 0.0; i<4; i++ ) {
	if ( base_comp [ i ] >= small ) { 
	    sum_bias += 
		( pos_base_comp_bias [ i ] / base_comp [ i ] ) +
		( pos_base_comp_bias [ i + 5 ] / base_comp [ i ] ) +
		( pos_base_comp_bias [ i + 10 ] / base_comp [ i ] );
	}
    }
    result [ 0 ] = sum_bias;

    for ( 
	  back_edge = 0,
	  front_edge = window_length,
	  middle = 1;
	  front_edge < edge_length;
          front_edge+=inc,
	  back_edge+=inc,
	  middle++) {

	base_type_pos1 = char_lookup [ edge [ back_edge ] ];
	base_type_pos2 = char_lookup [ edge [ back_edge + 1 ] ];
	base_type_pos3 = char_lookup [ edge [ back_edge + 2 ] ];

	base_comp [ base_type_pos1 ] -= 1.0;
	base_comp [ base_type_pos2 ] -= 1.0;
	base_comp [ base_type_pos3 ] -= 1.0;
	pos_base_comp [ base_type_pos1 ] -= 1.0;
	pos_base_comp [ base_type_pos2 + 5 ] -= 1.0;
	pos_base_comp [ base_type_pos3 + 10 ] -= 1.0;

	base_type_pos1 = char_lookup [ edge [ front_edge ] ];
	base_type_pos2 = char_lookup [ edge [ front_edge + 1 ] ];
	base_type_pos3 = char_lookup [ edge [ front_edge + 2 ] ];

	base_comp [ base_type_pos1 ] += 1.0;
	base_comp [ base_type_pos2 ] += 1.0;
	base_comp [ base_type_pos3 ] += 1.0;
	pos_base_comp [ base_type_pos1 ] += 1.0;
	pos_base_comp [ base_type_pos2 + 5 ] += 1.0;
	pos_base_comp [ base_type_pos3 + 10 ] += 1.0;


	for ( i=0; i<4; i++ ) {
	    t_comp = base_comp [ i ] / 3.0;
	    pos_base_comp_bias [ i ] = 
		fabs ( pos_base_comp [ i ] - t_comp );
	    pos_base_comp_bias [ i + 5 ] = 
		fabs ( pos_base_comp [ i + 5 ] - t_comp );
	    pos_base_comp_bias [ i + 10 ] = 
		fabs ( pos_base_comp [ i + 10 ] - t_comp );
	}

	for ( i=0, sum_bias = 0.0; i<4; i++ ) {
	    if ( base_comp [ i ] >= small ) { 
		sum_bias += 
		    ( pos_base_comp_bias [ i ] / base_comp [ i ] ) +
			( pos_base_comp_bias [ i + 5 ] / base_comp [ i ] ) +
			    ( pos_base_comp_bias [ i + 10 ] / base_comp [ i ] );
	    }
	}
	result [ middle ] = sum_bias;
    }
    /****************** middle game **********************/

    for ( 
	  back_edge = start,
	  front_edge = back_edge + window_length;
	  front_edge <= end;
          front_edge+=inc,
	  back_edge+=inc,
	  middle++) {

	base_type_pos1 = char_lookup [ seq [ back_edge ] ];
	base_type_pos2 = char_lookup [ seq [ back_edge + 1 ] ];
	base_type_pos3 = char_lookup [ seq [ back_edge + 2 ] ];

	base_comp [ base_type_pos1 ] -= 1.0;
	base_comp [ base_type_pos2 ] -= 1.0;
	base_comp [ base_type_pos3 ] -= 1.0;
	pos_base_comp [ base_type_pos1 ] -= 1.0;
	pos_base_comp [ base_type_pos2 + 5 ] -= 1.0;
	pos_base_comp [ base_type_pos3 + 10 ] -= 1.0;

	base_type_pos1 = char_lookup [ seq [ front_edge ] ];
	base_type_pos2 = char_lookup [ seq [ front_edge + 1 ] ];
	base_type_pos3 = char_lookup [ seq [ front_edge + 2 ] ];

	base_comp [ base_type_pos1 ] += 1.0;
	base_comp [ base_type_pos2 ] += 1.0;
	base_comp [ base_type_pos3 ] += 1.0;
	pos_base_comp [ base_type_pos1 ] += 1.0;
	pos_base_comp [ base_type_pos2 + 5 ] += 1.0;
	pos_base_comp [ base_type_pos3 + 10 ] += 1.0;

	for ( i=0; i<4; i++ ) {
	    t_comp = base_comp [ i ] / 3.0;
	    pos_base_comp_bias [ i ] = 
		fabs ( pos_base_comp [ i ] - t_comp );
	    pos_base_comp_bias [ i + 5 ] = 
		fabs ( pos_base_comp [ i + 5 ] - t_comp );
	    pos_base_comp_bias [ i + 10 ] = 
		fabs ( pos_base_comp [ i + 10 ] - t_comp );
	}

	for ( i=0, sum_bias = 0.0; i<4; i++ ) {
	    if ( base_comp [ i ] >= small ) { 
		sum_bias += 
		    ( pos_base_comp_bias [ i ] / base_comp [ i ] ) +
			( pos_base_comp_bias [ i + 5 ] / base_comp [ i ] ) +
			    ( pos_base_comp_bias [ i + 10 ] / base_comp [ i ] );
	    }
	}
	result [ middle ] = sum_bias;
    }
    xfree ( edge );

    /****************** right end game **********************/

    edge = seq_right_end ( seq, seq_length, end, window_length, inc);
    if ( ! (edge) ) return -1;

    edge_length = strlen ( edge );

    for ( 
	  back_edge = 0,
	  front_edge = back_edge + window_length;
	  front_edge < edge_length;
          front_edge+=inc,
	  back_edge+=inc,
	  middle++) {

	base_type_pos1 = char_lookup [ edge [ back_edge ] ];
	base_type_pos2 = char_lookup [ edge [ back_edge + 1 ] ];
	base_type_pos3 = char_lookup [ edge [ back_edge + 2 ] ];

	base_comp [ base_type_pos1 ] -= 1.0;
	base_comp [ base_type_pos2 ] -= 1.0;
	base_comp [ base_type_pos3 ] -= 1.0;
	pos_base_comp [ base_type_pos1 ] -= 1.0;
	pos_base_comp [ base_type_pos2 + 5 ] -= 1.0;
	pos_base_comp [ base_type_pos3 + 10 ] -= 1.0;

	base_type_pos1 = char_lookup [ edge [ front_edge ] ];
	base_type_pos2 = char_lookup [ edge [ front_edge + 1 ] ];
	base_type_pos3 = char_lookup [ edge [ front_edge + 2 ] ];

	base_comp [ base_type_pos1 ] += 1.0;
	base_comp [ base_type_pos2 ] += 1.0;
	base_comp [ base_type_pos3 ] += 1.0;
	pos_base_comp [ base_type_pos1 ] += 1.0;
	pos_base_comp [ base_type_pos2 + 5 ] += 1.0;
	pos_base_comp [ base_type_pos3 + 10 ] += 1.0;

	for ( i=0; i<4; i++ ) {
	    t_comp = base_comp [ i ] / 3.0;
	    pos_base_comp_bias [ i ] = 
		fabs ( pos_base_comp [ i ] - t_comp );
	    pos_base_comp_bias [ i + 5 ] = 
		fabs ( pos_base_comp [ i + 5 ] - t_comp );
	    pos_base_comp_bias [ i + 10 ] = 
		fabs ( pos_base_comp [ i + 10 ] - t_comp );
	}

	for ( i=0, sum_bias = 0.0; i<4; i++ ) {
	    if ( base_comp [ i ] >= small ) { 
		sum_bias += 
		    ( pos_base_comp_bias [ i ] / base_comp [ i ] ) +
			( pos_base_comp_bias [ i + 5 ] / base_comp [ i ] ) +
			    ( pos_base_comp_bias [ i + 10 ] / base_comp [ i ] );
	    }
	}
	result [ middle ] = sum_bias;
    }
    xfree ( edge );
    /* for some sizes of active region we do not calculate a value for
       the last element. Here we set the last value in the array to the 
       value of the last element calculated so we are sure it is sensible */

    result[ num_results - 1 ] = result[ middle - 1 ];
    return 0;
}

int do_pos_base_bias ( char seq[], int seq_length, CodRes1 *results ) 
{
    double m;
    int res;

    res = get_pos_base_bias ( seq, seq_length, 
			     results->user_start, 
			     results->user_end,
			     results->frame1,
			     results->num_results,
			     results->window_length);
    if ( res ) return res;
    /*
    printf(" pos_base_bias results\n");
    for ( i=0, j=0; i<results->num_results; i++,j+=3 ) {
	printf("i %d %f\n",j+1, results->frame1[i]);
    }
    */
    m = max_double_array ( results->frame1, results->num_results);
    results->max = m;
    m = min_double_array ( results->frame1, results->num_results);
    results->min = m;

    return 0;
}

void init_codon_table(double codon_table[4][4][4])
{
    int i, j, k;
    for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) {
	    for (k = 0; k < 4; k++) {
		codon_table[i][j][k] = 0.0;
	    }
	}
    }
}

void calc_codon_usage(char *seq,
		     int seq_length,
		     double codon_table[4][4][4])
		     
{
    int i;
    int inc = 3;
    int *genetic_code_idx = get_genetic_code_idx(0);

    /* make sure seq_length is a multiple of 3 */
    seq_length = (seq_length / inc) * inc;

    for (i = 0; i < seq_length; i+=inc) {
	if (legal_codon(&seq[i])) {
		codon_table 
		    [genetic_code_idx[char_lookup[seq[i  ]]]]
		    [genetic_code_idx[char_lookup[seq[i+1]]]]
		    [genetic_code_idx[char_lookup[seq[i+2]]]]++;
#ifdef DEBUG
	    printf("%c%c%c %f\n", seq[i], seq[i+1], seq[i+2], 
		   codon_table 
		   [genetic_code_idx[char_lookup[seq[i  ]]]]
		   [genetic_code_idx[char_lookup[seq[i+1]]]]
		   [genetic_code_idx[char_lookup[seq[i+2]]]]);
#endif
	}
    }
}

