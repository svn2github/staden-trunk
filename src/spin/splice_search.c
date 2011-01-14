#include <staden_config.h>

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "misc.h"
#include "array_arith.h"
#include "splice_search.h"
#include "dna_utils.h"
#include "getfile.h"

/* 7/1/99 johnt - need to explicitly import globals from dlls in Visual C++ */
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif


/* gene search by signal */

/* could put this little bit in dna_utils FIXME*/
extern DLL_IMPORT int unknown_char;

/*	Note that next two require char_match to be set */
#define SEQ_MISMATCH_ODD(a,intb) ( ( (char_match[a] < unknown_char) && (char_match[a]) == (intb)) ? 0 : 1 )

typedef struct _MatchPair {
    int  symbol;
    int  offset;
} MatchPair;	  /* within a motif, symbol symbol at offset offset
		     must match */

typedef struct _MatchMask {
    int num_elements;
    int depth;         /* depends on charset size 4 for dna, 22 for protein */
    MatchPair *pair;
} MatchMask;	/* within a motif num_elements characters are
		   100% conserved. Here we define the symbol
		   and offset pairs */

MatchMask *init_MatchMask (
    int length,
    int depth ) {

    MatchMask *m;
    MatchPair *p;

    if ( ( NULL == ( m = (MatchMask* ) xmalloc ( sizeof (MatchMask) )))) return NULL;
    if ( ( NULL == ( p = (MatchPair* ) xmalloc ( sizeof (MatchPair) * length )))) return NULL;

    m->num_elements = length;	/* set to maximum */
    m->depth = depth;           /* depends on char set size */
    m->pair = p;
    return m;
}

typedef struct _WtmatrixSpec {
    double *matrix;		/* the matrix stored as a vector */
    int    length;		/* matrix length - ied motif length */
    int    depth;		/* matrix depth - depends on charset size */
    double min;			/* cutoff score */
    double max;			/* maximum possible score */
    int    mark_pos;		/* offset to mark position. ied mark position+mark_pos */
} WtmatrixSpec;

typedef struct _WeightMatrixCounts {
  int *counts;
  int length;
  int depth;
  double min;
  double max;
  int mark_pos;
} WeightMatrixCounts;


WeightMatrixCounts *initWeightMatrixCounts (
  int length,
  int depth ) {

  WeightMatrixCounts *w;
  int *counts;

  if ( ( NULL == ( w = (WeightMatrixCounts* ) xmalloc ( sizeof (WeightMatrixCounts))))) return NULL;

  if ( ( NULL == ( counts = (int* ) xmalloc ( sizeof (int)*length*depth)))) return NULL;

  w->counts = counts;
  w->depth  = depth;
  return w;
}

void free_WeightMatrixCounts ( WeightMatrixCounts *w ) {
   if (w) {
     free ( w->counts );
     free ( w );
   }
}

void free_MatchMask ( MatchMask *m ) {

    xfree ( m->pair );
    xfree ( m );
}

WtmatrixSpec *init_Wtmatrix ( WeightMatrixCounts *c ) {

    WtmatrixSpec *w;
    double *m;
    int i, j;

    j = char_set_size * c->length;

    if ( ( NULL == ( w = (WtmatrixSpec* ) xmalloc ( sizeof (WtmatrixSpec) )))) return NULL;
    if ( ( NULL == ( m = (double* ) xmalloc ( sizeof (double) * j )))) return NULL;

    for ( i=0; i<j; i++) { 
	m[i] = 0;
    }
    w->length = c->length;
    w->depth  = char_set_size;
    w->min  = c->min;
    w->max  = c->max;
    w->mark_pos  = c->mark_pos;
    w->matrix = m;
    return w;
}

void free_WtmatrixSpec ( WtmatrixSpec *w ) {

    xfree ( w->matrix );
    xfree ( w );
}

WtmatrixRes *init_WtmatrixRes (int number_of_res,
			       WtmatrixSpec *s ) 
{

    WtmatrixRes *r;
    Wtmatch **match = NULL;

    if ( ( NULL == ( r = (WtmatrixRes* ) xmalloc ( sizeof (WtmatrixRes) )))) 
	return NULL;

    /* need to allow zero matches! */

    if ( number_of_res > 0 ) {
	if ( ( NULL == ( match = (Wtmatch** ) xmalloc ( sizeof (Wtmatch *) * number_of_res )))) return NULL;
    }
    r->match = match;
    r->number_of_res = number_of_res;
    r->length = s->length;
    r->mark_pos = s->mark_pos;
    r->min = s->min;
    r->max = s->max;
    return r;
}

void free_WtmatrixRes ( WtmatrixRes *r ) {
    int i;

    for ( i=0;i<r->number_of_res;i++ ) {
	xfree ( r->match[i] );
    }
    /* only free match if was allocated, ie if r->number_of_res > 0 */
    if (r->number_of_res > 0) 
	xfree ( r->match );

    xfree ( r );
}

void free_wt_setup ( WtmatrixSpec *w, 
	       MatchMask *wc, 
	       WtmatrixRes *results ) { 

    /* note we assume that all are initialised to NULL before malloc */
    if ( w ) free_WtmatrixSpec ( w );
    if ( wc ) free_MatchMask ( wc );
    if ( results ) free_WtmatrixRes ( results );
}

void free_splice_setup ( 
	       WeightMatrixCounts *ied_c,
	       WeightMatrixCounts *eia_c,
               WtmatrixSpec *ied, 
	       WtmatrixSpec *eia,
	       MatchMask *iedc, 
	       MatchMask *eiac,
	       WtmatrixRes *results_ied, 
	       WtmatrixRes *results_eia ) { 

    /* note we assume that all are initialised to NULL before malloc */
    if ( ied_c ) free_WeightMatrixCounts ( ied_c );
    if ( eia_c ) free_WeightMatrixCounts ( eia_c );
    if ( ied ) free_WtmatrixSpec ( ied );
    if ( eia ) free_WtmatrixSpec ( eia );
    if ( iedc ) free_MatchMask ( iedc );
    if ( eiac ) free_MatchMask ( eiac );
    if ( results_ied ) free_WtmatrixRes ( results_ied );
    if ( results_eia ) free_WtmatrixRes ( results_eia );
}


void free_splice_results2 ( SpliceResults *s ) {

    if ( s->ied_f1 ) free_WtmatrixRes ( s->ied_f1 );
    if ( s->ied_f2 ) free_WtmatrixRes ( s->ied_f2 );
    if ( s->ied_f3 ) free_WtmatrixRes ( s->ied_f3 );
    if ( s->eia_f1 ) free_WtmatrixRes ( s->eia_f1 );
    if ( s->eia_f2 ) free_WtmatrixRes ( s->eia_f2 );
    if ( s->eia_f3 ) free_WtmatrixRes ( s->eia_f3 );
}


void free_splice_results1 ( 

	       WtmatrixRes *splice_ied_f1,
	       WtmatrixRes *splice_ied_f2,
	       WtmatrixRes *splice_ied_f3,
	       WtmatrixRes *splice_eia_f1,
	       WtmatrixRes *splice_eia_f2,
	       WtmatrixRes *splice_eia_f3 ) {

    if ( splice_ied_f1 ) free_WtmatrixRes ( splice_ied_f1 );
    if ( splice_ied_f2 ) free_WtmatrixRes ( splice_ied_f2 );
    if ( splice_ied_f3 ) free_WtmatrixRes ( splice_ied_f3 );
    if ( splice_eia_f1 ) free_WtmatrixRes ( splice_eia_f1 );
    if ( splice_eia_f2 ) free_WtmatrixRes ( splice_eia_f2 );
    if ( splice_eia_f3 ) free_WtmatrixRes ( splice_eia_f3 );

}


int get_wtm_cons_chars ( int counts[], MatchMask *m ) {

    /* from counts find the 100% conserved positions */

    double *sum_column;
    int i,j,k,total,col_max,col_index=0,mask_index;

    if ( ( NULL == ( sum_column = (double* ) xmalloc ( sizeof (double) * m->num_elements )))) return -1;
    for ( i=0,mask_index=0; i<m->num_elements; i++ ) {
	for ( j=0,total=0,col_max=0; j<m->depth; j++ ) {
	    k = j * (m->num_elements) + i;
	    total += counts [ k ];
	    if ( counts [ k ] > col_max ) {
		col_max = counts [ k ];
		col_index = j; 
	    }
	}
	if ( total == col_max ) { 
	    m->pair[mask_index].offset = i;
	    m->pair[mask_index++].symbol = col_index;
	}
    }
    m->num_elements = mask_index;
    xfree ( sum_column );
    return 0;
}

int mask_match ( char *seq, int seq_len, int seq_pos, MatchMask *m ) {

    /* does MatchMask match seq starting at pos seq_pos */
    /* if so return position in seq, else return seq_len + 10 */

    int i,j,k, last_start;
    last_start = seq_len - m->pair[(m->num_elements)-1].offset - 1;
    for ( i = seq_pos; i < last_start; i++ ) {
	for ( j = 0, k = 1; j < m->num_elements; j++ ) {
	    if ( SEQ_MISMATCH_ODD( seq[i+m->pair[j].offset], m->pair[j].symbol )) {
		k = 0;
		break;
	    }
	}
	if ( k ) return i;
    }
    return seq_len + 10;
}

int do_mask_match ( char *seq, int seq_len, int user_start, int user_end, MatchMask *m ) {

    /* find all places where MatchMask matches */

    int start, end;
    int i,j, last_start;

    start = user_start - 1;
    end = user_end - 1;
    last_start = end - m->pair[(m->num_elements)-1].offset - 1;
    for ( i = start; i < last_start; i++ ) {
	if ( (j = mask_match ( seq, end, i, m )) <= end ) {
#ifdef DEBUG
	    printf("i %d j %d\n",i,j);
#endif
	    i = j;
	}
	else {
	    break;
	}
    }
    return 0;
}

int do_wt_search_cs ( char seq[], int seq_length, int user_start, int user_end,
		  WtmatrixSpec *w, MatchMask *c, WtmatrixRes *r ) {

    int start, end, max_start, left_start, num_res, matrix_pos, seq_pos, match_pos;
    double score;
    Wtmatch *m;

    start = user_start - 1;
    end   = user_end - 1;

    /* loop for all matrix left starts. Note the matrix is stored in a vector length*depth */

    max_start = end - ( w->length - 1 );

    for ( left_start = start, num_res = 0; left_start <= max_start; left_start++ ) {

	/* loop for the matrix length */

	if ( (match_pos = mask_match ( seq, end, left_start, c )) <= max_start ) {
	    left_start = match_pos;
	    for ( matrix_pos = 0, seq_pos = left_start, score = 0.0; 
		 matrix_pos < w->length && seq_pos < user_end; matrix_pos++, seq_pos++ ) {

		score += w->matrix [ char_lookup [ seq [ seq_pos ]] * (w->length) + matrix_pos ];
	    }
	    if ( score >= w->min ) {

		if ( ( NULL == ( m = (Wtmatch* ) xmalloc ( sizeof (Wtmatch ) )))) return -3;
		m->pos = left_start + w->mark_pos;
		m->score = score;
		m->seq = &seq[left_start];

		if ( num_res == r->number_of_res ) {

		    if ( (NULL == (r->match = xrealloc ( r->match, sizeof(Wtmatch *) * 
			 ( r->number_of_res + r->number_of_res/2) )))) return -2;

		    r->number_of_res = r->number_of_res + r->number_of_res/2;
		}

		r->match[num_res++] = m;
	    }
	}
	else {
	    break;
	}
    }
    r->number_of_res = num_res;

    /* set array sizes to what was used */
    if ( num_res ) {

	if ( (NULL == (r->match = xrealloc ( r->match, sizeof(Wtmatch *) * 
					    ( num_res) )))) return -3;
    }
    return 0;
}

int do_mask_match_wt ( char *seq, int seq_len, int user_start, int user_end, MatchMask *m ) {

    /* find all places where MatchMask matches */

    int start, end;
    int i,j, last_start;

    start = user_start - 1;
    end = user_end - 1;
    last_start = end - m->pair[(m->num_elements)-1].offset - 1;
    for ( i = start; i < last_start; i++ ) {
	if ( (j = mask_match ( seq, end, i, m )) < last_start ) {
#ifdef DEBUG
	    printf("i %d j %d\n",i,j);
#endif
	    i = j;
	}
	else {
	    break;
	}
    }
    return 0;
}



int get_wt_weights_old ( int counts[], WtmatrixSpec *w ) {

    /* from counts get the weights as done in nip */
    /* but note that adding 1 to each element is NOT
       the most sensitive method to use - adding 1/4
       would be better
    */

    double *sum_column;
    int i,j,k,total;

    if ( ( NULL == ( sum_column = (double* ) xmalloc ( sizeof (double) * w->length )))) return -1;
    for ( i=0; i<w->length; i++ ) {
	total = 0;
	for ( j=0; j<w->depth-1; j++ ) {
	    k = j * (w->length) + i;
	    total += counts [ k ];
	    w->matrix [ k ] =  counts [ k ];
	}

	sum_column [ i ] = (double)total;
	k = j * (w->length) + i;
	w->matrix [ k ] = sum_column [ i ] / (w->depth-1);
    }

    for ( i=0; i<w->length; i++ ) {
	for ( j=0; j<w->depth; j++ ) {
	    k = j * (w->length) + i;
	    /*printf(" i %d j %d w %f s %f ",i,j,w->matrix[k],sum_column[i]);*/
	    w->matrix[ k ] = log ( ((w->matrix[ k ] + 1) / sum_column [ i ])/0.25);
	    /*printf(" %f\n",w->matrix[k]);*/
	}
    }
    xfree ( sum_column );
    return 0;
}


int get_wt_weights ( int counts[], WtmatrixSpec *w ) {

  /* get weights as log-odds
   * assume that, for the random model, the probability of each base is 0.25
   */

  /* note we add small to avoid zero counts! */
    double *sum_column;
    /*double small = 1.0;*/
    double small;
    int i,j,k,total;

    if ( ( NULL == ( sum_column = (double* ) xmalloc ( sizeof (double) * w->length )))) return -1;

    /* loop for each column of the matrix */
    for ( i=0; i<w->length; i++ ) {
	total = 0;
	/* get the counts for each column */
	for ( j=0; j<w->depth-1; j++ ) {
	    k = j * (w->length) + i;
	    total += counts [ k ];
	}
	small = 1.0;
	if ( total ) small = 1.0/(double)total;
	sum_column [ i ] = (double)total + (w->depth-1) * small;

	/* get the counts for each character this column */
	for ( j=0; j<w->depth-1; j++ ) {
	    k = j * (w->length) + i;
	    w->matrix [ k ] =  (double)counts [ k ] + small; /* make > 0 */
	}
	k = j * (w->length) + i;
	/* set unknown character to mean count */
	w->matrix [ k ] = sum_column [ i ] / (w->depth-1);
    }

    /* loop for each column of the matrix */
    for ( i=0; i<w->length; i++ ) {
	/* set the weights for each character this column */
	for ( j=0; j<w->depth; j++ ) {
	    k = j * (w->length) + i;
	    /*printf(" i %d j %d w %f s %f ",i,j,w->matrix[k],sum_column[i]);*/
	    w->matrix [ k ] = log ( (w->matrix [ k ] / sum_column [ i ])/0.25);
	    /*printf(" %f\n",w->matrix[k]);*/
	}
    }
    xfree ( sum_column );
    return 0;
}


int do_wt_search ( char seq[], int seq_length, int user_start, int user_end,
		  WtmatrixSpec *w, WtmatrixRes *r) {

    int start, end, max_start, left_start, num_res, matrix_pos, seq_pos;
    double score;
    Wtmatch *m;

    start = user_start - 1;
    end   = user_end - 1;

    /* loop for all matrix left starts. Note the matrix is stored in a vector length*depth */

    max_start = end - ( w->length - 1 );

    for ( left_start = start, num_res = 0; left_start <= max_start; left_start++ ) {

	/* loop for the matrix length */

	for ( matrix_pos = 0, seq_pos = left_start, score = 0.0; 
	      matrix_pos < w->length; matrix_pos++, seq_pos++ ) {

	    score += w->matrix [ char_lookup [ seq [ seq_pos ]] * (w->length) + matrix_pos ];
	}
	if ( score >= w->min ) {

	    if ( ( NULL == ( m = (Wtmatch* ) xmalloc ( sizeof (Wtmatch ) )))) return -3;
	    m->pos = left_start + w->mark_pos;
	    m->score = score;
	    m->seq = &seq[left_start];

	    if ( num_res == r->number_of_res ) {

		if ( (NULL == (r->match = xrealloc ( r->match, sizeof(Wtmatch *) * 
		     ( r->number_of_res + r->number_of_res/2) )))) return -2;

		r->number_of_res = r->number_of_res + r->number_of_res/2;
	    }

	    r->match[num_res++] = m;
	}
    }
    r->number_of_res = num_res;

    /* set array sizes to what was used */
    if ( num_res ) {

	if ( (NULL == (r->match = xrealloc ( r->match, sizeof(Wtmatch *) * 
					    ( num_res) )))) return -3;
    }
    return 0;
}

WeightMatrixCounts *read_weight_matrix ( FILE *file_p, int char_set ) {

  char line[200],char_type[2];
  int NUM_COLUMNS = 20;
  int block, n_blocks, row, char_index, element;
  int wl, wmp, depth;
  double wmn, wmx;
  WeightMatrixCounts *w = NULL;

/* deal with charsetsize */
  if ( 5 == char_set ) {
    depth = 4;
  }
  else {
    depth = 22;
  }
  /* read and discard title */
      if ( ( fscanf (file_p, "%200[^\n]\n", line) == 0 ) ) return NULL;

      /* printf("%s\n",line); */

  /* read the matrix length, mark_pos, min and max */

  if ( fscanf ( file_p, "%d %d %lf %lf\n", &wl, &wmp,
		&wmn, &wmx ) != 4 ) return NULL;

  /* sanity check */

  if ( wl <= 0 ) return NULL;

  if (NULL == (w = initWeightMatrixCounts ( wl, depth ))) return NULL;
  w->length = wl;
  w->mark_pos = wmp;
  w->min = wmn;
  w->max = wmx;

  /* a table has <21 columns per block and depth rows per block
   * plus a title and 2 rows that can be ignored. Current tables have the
   * character order tcag or McLachlan order.
   * The number of blocks is 1 + (w->length - 1) / 20
   * The number of useful rows in a block is depth
   * The number of columns depends on the length and block number.
   * In the count table we store data in acgt order: all the a's
   * followed by all the c's etc and this not the order in which
   * they are read: old tables are in tcag order AND there may be
   * more than a single block. Let element be the element number
   * in counts, then 
   * element = length * char_index + block * NUM_COLUMNS + column
   *
   */

  n_blocks = 1 + ( w->length - 1 ) / NUM_COLUMNS;

  for ( block=0; block<n_blocks; block++ ) {

	
      if ( ( fscanf (file_p, "%[^\n]\n", line) == 0 ) ) return NULL;
      if ( ( fscanf (file_p, "%[^\n]\n", line) == 0 ) ) return NULL;
      for ( row=0;row<depth; row++ ) {
	  if ( ( fscanf (file_p, " %c", char_type) == 0 ) ) return NULL;
	  char_index = char_lookup [ char_type[0] ];
	  element = w->length * char_index + block * NUM_COLUMNS;
	  while ( fscanf (file_p, " %d", &(w->counts[element++]))>0);

    }
  }
  return w;
}

int splice_search (char seq[], int seq_length, int user_start, int user_end,
		   char *filename_ied, char *filename_eia,
		   SpliceResults *splice_result)
{

    WeightMatrixCounts *ied_c=NULL, *eia_c=NULL;
    WtmatrixSpec *ied=NULL, *eia=NULL;
    WtmatrixRes *results_ied=NULL, *results_eia=NULL;
    WtmatrixRes *splice_ied_f1=NULL,*splice_ied_f2=NULL,*splice_ied_f3=NULL; /* intron-exon (donors) frames 1,2,3*/
    WtmatrixRes *splice_eia_f1=NULL,*splice_eia_f2=NULL,*splice_eia_f3=NULL; /* exon-intron (acceptors) frames 1,2,3*/
    Wtmatch *m=NULL;
    MatchMask *iedc=NULL, *eiac=NULL;
    int i,num_res;
    int ret = -1;
    int j, n1, n2, n3;
    FILE *file_p;
/*
    char filename_ied[]="$STADTABL/ied.wts";
    char filename_eia[]="$STADTABL/eia.wts";
*/
    char fileexpand_ied[FILENAME_MAX+1];
    char fileexpand_eia[FILENAME_MAX+1];

    if (1 != expandpath(filename_ied, fileexpand_ied))
	return -1;
    if (1 != expandpath(filename_eia, fileexpand_eia))
	return -1;

#ifdef DEBUG    
    printf("%s %s\n", fileexpand_ied, fileexpand_eia);
#endif

    if (NULL == (file_p = fopen ( fileexpand_ied, "r"))) return -1;
    if (NULL == (ied_c = read_weight_matrix ( file_p, char_set_size ))) {
	fclose ( file_p );
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }
    fclose ( file_p );

    /* sanity check */

    if ( (user_end - user_start + 1) < ied_c->length ) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }

    if ( (ied = init_Wtmatrix ( ied_c )) == NULL ) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }
    if ( get_wt_weights ( ied_c->counts, ied ) ) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }

    if (NULL == (file_p = fopen ( fileexpand_eia, "r"))) return -1;
    if (NULL == (eia_c = read_weight_matrix ( file_p, char_set_size ))) {
	fclose ( file_p );
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }
    fclose ( file_p );

    /* sanity check */

    if ( (user_end - user_start + 1) < eia_c->length ) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }

    if ( (eia = init_Wtmatrix ( eia_c )) == NULL ) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }

    if ( get_wt_weights ( eia_c->counts, eia ) ) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }

    if (( iedc = init_MatchMask (ied_c->length,ied_c->depth)) == NULL) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }
    if (( ret = get_wtm_cons_chars ( ied_c->counts, iedc ))) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }
    if (( eiac = init_MatchMask (eia_c->length,eia_c->depth)) == NULL) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }
    if (( ret = get_wtm_cons_chars ( eia_c->counts, eiac ))) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }

    num_res = 1 + seq_length / 20; /* FIXME: should use the probability calcs to guess this */

    if (( results_ied = init_WtmatrixRes (num_res,ied)) == NULL) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }

    if (( results_eia = init_WtmatrixRes (num_res,eia)) == NULL) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }

    if (( ret = do_wt_search_cs ( seq, seq_length, user_start, user_end, 
				 ied, iedc, results_ied ))) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }

#ifdef DEBUG
    printf("number of ied matches %d\n",results_ied->number_of_res);
#endif
    /* count matches for each frame */

    n1 = n2 = n3 = 0;
    for (i=0;i<results_ied->number_of_res;i++) {
	j = results_ied->match[i]->pos%3;
	if ( j == 0 ) {
	    n3++;
	}
	else if ( j == 1 ) {
	    n1++;
	}
	else if ( j == 2 ) {
	    n2++;
	}
    }

    /* now set mark_pos to zero so it is copied as such into the new data structures */

    results_ied->mark_pos = 0;

    if (( splice_ied_f1 = init_WtmatrixRes (n1,ied)) == NULL) {

	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -2;
    }

    if (( splice_ied_f2 = init_WtmatrixRes (n2,ied)) == NULL) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -2;
    }
    if (( splice_ied_f3 = init_WtmatrixRes (n3,ied)) == NULL) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -2;
    }
    n1 = n2 = n3 = 0;
    for (i=0;i<results_ied->number_of_res;i++) {
	j = results_ied->match[i]->pos%3;
	if ( ( NULL == ( m = (Wtmatch* ) xmalloc ( sizeof (Wtmatch ) )))) {
	    free_splice_results1 (
				 splice_ied_f1,splice_ied_f2,splice_ied_f3,
				 splice_eia_f1,splice_eia_f2,splice_eia_f3 );
	    /* if ( results ) xfree ( results ); */
	    free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	    return -3;
	}
	m->pos   = results_ied->match[i]->pos;
	m->score = results_ied->match[i]->score;
	m->seq   = results_ied->match[i]->seq;

	if ( j == 0 ) {
	    splice_ied_f3->match[n3++] = m;
	}
	else if ( j == 1 ) {
	    splice_ied_f1->match[n1++] = m;
	}
	else if ( j == 2 ) {
	    splice_ied_f2->match[n2++] = m;
	}
    }

#ifdef DEBUG
    for (i=0;i<splice_ied_f1->number_of_res;i++) {
	printf("i %d pos %d 1 score %f\n",i,splice_ied_f1->match[i]->pos,splice_ied_f1->match[i]->score);
    }
    for (i=0;i<splice_ied_f2->number_of_res;i++) {
	printf("i %d pos %d 2 score %f\n",i,splice_ied_f2->match[i]->pos,splice_ied_f2->match[i]->score);
    }
    for (i=0;i<splice_ied_f3->number_of_res;i++) {
	printf("i %d pos %d 3 score %f\n",i,splice_ied_f3->match[i]->pos,splice_ied_f3->match[i]->score);
    }
#endif

#ifdef REMOVE
    if (( results_ied = init_WtmatrixRes (num_res,ied)) == NULL) return -2;
    if ( ret ) {

	free_WtmatrixSpec ( ied );
	free_WtmatrixRes ( results_ied );
	return ret;
    }
#endif
    if (( ret = do_wt_search_cs ( seq, seq_length, user_start, user_end, 
				 eia, eiac, results_eia ))) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -1;
    }

#ifdef DEBUG
    printf("number of eia matches %d\n",results_eia->number_of_res);
#endif

    /* count matches for each frame */

    n1 = n2 = n3 = 0;
    for (i=0;i<results_eia->number_of_res;i++) {
	j = results_eia->match[i]->pos%3;
	if ( j == 0 ) {
	    n3++;
	}
	else if ( j == 1 ) {
	    n1++;
	}
	else if ( j == 2 ) {
	    n2++;
	}
    }

    /* now set mark_pos to zero so it is copied as such into the new data structures */
    results_eia->mark_pos = 0;

    if (( splice_eia_f1 = init_WtmatrixRes (n1,eia)) == NULL) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -2;
    }

    if (( splice_eia_f2 = init_WtmatrixRes (n2,eia)) == NULL) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -2;
    }
    if (( splice_eia_f3 = init_WtmatrixRes (n3,eia)) == NULL) {
	free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	return -2;
    }

    n1 = n2 = n3 = 0;
    for (i=0;i<results_eia->number_of_res;i++) {
	j = results_eia->match[i]->pos%3;

	if ( ( NULL == ( m = (Wtmatch* ) xmalloc ( sizeof (Wtmatch ) )))) {

	    free_splice_results1 (
				 splice_ied_f1,splice_ied_f2,splice_ied_f3,
				 splice_eia_f1,splice_eia_f2,splice_eia_f3 );
	    /* if ( results ) xfree ( results ); */
	    free_splice_setup ( ied_c, eia_c, ied, eia,
			   iedc, eiac,
			   results_ied, results_eia );
	    return -3;
	}

	m->pos   = 1 + results_eia->match[i]->pos;
	m->score = results_eia->match[i]->score;
	m->seq   = results_eia->match[i]->seq;

	if ( j == 0 ) {
	    splice_eia_f3->match[n3++] = m;
	}
	else if ( j == 1 ) {
	    splice_eia_f1->match[n1++] = m;
	}
	else if ( j == 2 ) {
	    splice_eia_f2->match[n2++] = m;
	}
    }

#ifdef DEBUG
    for (i=0;i<splice_eia_f1->number_of_res;i++) {
	printf("i %d pos %d 1 score %f\n",i,splice_eia_f1->match[i]->pos,splice_eia_f1->match[i]->score);
    }
    for (i=0;i<splice_eia_f2->number_of_res;i++) {
	printf("i %d pos %d 2 score %f\n",i,splice_eia_f2->match[i]->pos,splice_eia_f2->match[i]->score);
    }
    for (i=0;i<splice_eia_f3->number_of_res;i++) {
	printf("i %d pos %d 3 score %f\n",i,splice_eia_f3->match[i]->pos,splice_eia_f3->match[i]->score);
    }
#endif
    splice_result->ied_f1 = splice_ied_f1;
    splice_result->ied_f2 = splice_ied_f2;
    splice_result->ied_f3 = splice_ied_f3;
    splice_result->eia_f1 = splice_eia_f1;
    splice_result->eia_f2 = splice_eia_f2;
    splice_result->eia_f3 = splice_eia_f3;

    free_splice_setup (ied_c, eia_c, ied, eia, iedc, eiac,
			results_ied, results_eia);

    return ret;
}

int weight_search ( char seq[], int seq_length, int user_start, int user_end,
                   char *filename, WtmatrixRes **weight_results) {

    WtmatrixSpec *w=NULL;
    WtmatrixRes *results=NULL;
    MatchMask *wc=NULL;
    WeightMatrixCounts *c;
    FILE *file_p;
/*    char filename[]="d.wts"; */ /* FIXME */
    int num_res;
    int ret = -1;

    if (NULL == (file_p = fopen ( filename, "r"))) return -1;
    if (NULL == (c = read_weight_matrix ( file_p, char_set_size ))) {
	fclose ( file_p );
	free_WeightMatrixCounts ( c );
	return -1;
    }
    fclose ( file_p );

    /* sanity check */

    if ( (user_end - user_start + 1) < c->length ) {
	free_WeightMatrixCounts ( c );
	return -1;
    }


    if ( (w = init_Wtmatrix ( c )) == NULL ) {
	free_wt_setup ( w, wc, results );
	free_WeightMatrixCounts ( c );
	return -1;
    }


    if ( get_wt_weights ( c->counts, w ) ) {
	free_wt_setup ( w, wc, results );
	free_WeightMatrixCounts ( c );
	return -1;
    }

    if (( wc = init_MatchMask (c->length,c->depth)) == NULL) { 
	free_wt_setup ( w, wc, results );
	return -1;
    }

    if (( ret = get_wtm_cons_chars ( c->counts, wc ))) {
	free_wt_setup ( w, wc, results );
	return -1;
    }

    free_WeightMatrixCounts ( c );

    num_res = 1 + seq_length / 20; /* FIXME: should use the probability calcs to guess this */

    if (( results = init_WtmatrixRes (num_res,w)) == NULL) {
	free_wt_setup ( w, wc, results );
	return -1;
    }

    /* if there are some conserved bases use them in the search */

    if ( wc->num_elements ) {

	if (( ret = do_wt_search_cs ( seq, seq_length, user_start, user_end, 
				 w, wc, results ))) {
	    free_wt_setup ( w, wc, results );
	    return -1;
	}
    }
    /* otherwise brute force search */
    else {
	if (( ret = do_wt_search ( seq, seq_length, user_start, user_end, 
				 w, results ))) {
	    free_wt_setup ( w, wc, results );
	    return -1;
	}
    }
#ifdef DEBUG
    printf("number of w matches %d\n",results->number_of_res);
    for (i=0;i<results->number_of_res;i++) {
	printf("i %d pos %d score %f\n",i,results->match[i]->pos,results->match[i]->score);
    }
#endif
/*    free_wt_setup ( w, wc, results );*/
    if ( w ) free_WtmatrixSpec ( w );
    if ( wc ) free_MatchMask ( wc );
    *weight_results = results;
    return ret;
}
