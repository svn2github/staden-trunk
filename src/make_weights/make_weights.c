/* make_weights: reads a set of aligned motifs, creates weights and writes
                 a weight matrix file */

#include <staden_config.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#include "misc.h" 
#include "dna_utils.h"

#define MAX_LINE 2048
#define MAX_MOTIFS 10000
#define NAME 0
#define SEQ 1
#define UNSPECIFIED 4096

/* 9/1/99 johnt - must explictly import globals from DLLs with Visual C++ */
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif


/* extern int *char_lookup;  9/1/99 johnt - Not used, and doesn't link under Windows */

typedef struct motif_spec_ {
  char *name;
  char *comment;
  char *seq;
  double score;
} Motif_spec;

typedef struct motif_specs_ {
  Motif_spec *motif_spec_ptr;
  int         number;
  int         length;
  int         name_length;
  int         depth;
  int         mark_pos;
  double      min;			/* cutoff score */
  double      max;			/* maximum possible score */
} Motif_specs;


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

  if ( ( NULL == ( counts = (int* ) xcalloc ( length*depth, sizeof (int))))) return NULL;

  w->counts = counts;
  w->depth  = depth;
  w->length  = length;
  w->mark_pos = 0;
  w->min = 0;
  w->max = 0;
  return w;
}

void free_WeightMatrixCounts ( WeightMatrixCounts *w ) {
   free ( w->counts );
   free ( w );
}

WtmatrixSpec *init_Wtmatrix ( WeightMatrixCounts *c ) {

    WtmatrixSpec *w;
    double *m;
    int i, j;

    j = c->depth * c->length;

    if ( ( NULL == ( w = (WtmatrixSpec* ) xmalloc ( sizeof (WtmatrixSpec) )))) return NULL;

    if ( ( NULL == ( m = (double* ) xmalloc ( sizeof (double) * j )))) return NULL;


    for ( i=0; i<j; i++) { 
	m[i] = 0;
    }

    w->length = c->length;
    w->depth  = c->depth;
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

Motif_specs* get_motif_info(FILE *fp_vf, int char_set ) {

  char line[MAX_LINE], *comment, *name = NULL, *seq = NULL;
  char *c, *s;
  Motif_spec *v;
  static Motif_specs vs;
  int item, in_item, number, length = 0, name_length = 0;

  if (NULL == (v = (Motif_spec *)malloc((sizeof(Motif_spec ) * MAX_MOTIFS)))) return NULL;

    /* format: 
     *
       name seq comment 
     *
     * comment is optional
     * comment might contain spaces! but none of the others can.
     */ 

  number = 0;
  while (fgets(line,MAX_LINE,fp_vf) != NULL) {

    for(c = line, in_item = 0, item = 0; item < 2; item++ ) {
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
	    /*printf("%s\n",s);*/
	    if ( item == NAME ) {
	      if (NULL == (name = (char *)malloc((strlen(s) + 1)*sizeof(char)))) return NULL;
	      strcpy ( name, s );
	      name_length = MAX ( name_length, strlen(s) );
	    }
	    else if ( item == SEQ ) {
	      if (NULL == (seq = (char *)malloc((strlen(s) + 1)*sizeof(char)))) return NULL;
	      strcpy ( seq, s );
	      length = MAX ( length, strlen(s) );
	    }
	  }
	  else {
	    /* item ends oddly */
	    printf("Warning: error parsing file\n");
	    break;
	  }
	}
	else {
	  /* item starts oddly */
	  printf("Warning: error parsing file\n");
	  break;
	}
      }
    }
    if ( item != 2 ) continue;
    /* if the comment is present, get it */
    comment = NULL;
    if ( strlen ( c ) ) {
      c[strlen(c)-1] = '\0';
      if (NULL == (comment = (char *)malloc((strlen(c) + 1)*sizeof(char)))) return NULL;
      /*printf("%s\n",c); */
      strcpy ( comment, c );
    }

    if ( number == MAX_MOTIFS ) {
      printf("Aborting: more than %d entries in file\n",MAX_MOTIFS);
      return NULL;
    }

    v[number].name  = name;
    v[number].seq = seq;
    v[number].comment = comment;
    v[number].score = -99999.0;
    number++;
  }

  vs.motif_spec_ptr = v;
  vs.number = number;
  vs.length = length;
  vs.name_length = name_length;
  vs.depth  = char_set;
  /*{
    int i;
      for(i=0;i<number;i++) {
	printf(" record %d\n",i);
	printf("%s\n",v[i].name);
	printf("%s\n",v[i].seq);
	if(v[i].comment)printf("%s\n",v[i].comment);
      }
    printf("number %d length %d\n",vs.number,vs.length);
  }*/
  return &vs;
}

int get_wt_weights ( int counts[], WtmatrixSpec *w ) {

  /* get weights as log-odds */

  /* note we add small to avoid zero counts! */
    double *sum_column;
    /*double small = 1.0;*/
    double small;
    int i,j,k,total;

    if ( ( NULL == ( sum_column = (double* ) xmalloc ( sizeof (double) * w->length )))) return -1;

    /* sum each columns counts to get total
     * set small to 1 or 1/total
     * set column total to total + small*char_set_size
     * set the weights (matrix) to counts + small
     * set unknown char (matrix) to mean for the column
     * set the weights (matrix) to log odds:
     *     assume in random model probability of each base type is 0.25
     *     so odds is p(model)/p(random) = p(model)/0.25
     *     p = matrix/total
     *     p = log ( p / 0.25 )
     *     matrix = p
    */

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
	/*printf("i %d t %d %f %f\n",i,total,sum_column[i],small);*/

	/* set the counts for each character this column */
	for ( j=0; j<w->depth-1; j++ ) {
	    k = j * (w->length) + i;
	    w->matrix [ k ] =  (double)counts [ k ] + small; /* make > 0 */
	/*printf("i %d k %d t %f %f\n",i,k,w->matrix[k],sum_column[i]);*/
	}
	k = j * (w->length) + i;
	/* set unknown character to mean count */
	w->matrix [ k ] = sum_column [ i ] / (w->depth-1);
	/*printf("n i %d t %f %d %f %f\n",i,sum_column[i],w->depth-1,w->matrix[k],t);*/
    }

    /* loop for each column of the matrix */
    for ( i=0; i<w->length; i++ ) {
	/* set the weights for each character this column */
	for ( j=0; j<w->depth; j++ ) {
	    k = j * (w->length) + i;
	    /*printf(" i %d j %d w %f s %f ",i,j,w->matrix[k],sum_column[i]);*/
	    /* old
	    p = w->matrix [ k ] / sum_column [ i ];
	    w->matrix [ k ] = log ( p / ( 1.0 - p ) );
	    */
	    w->matrix [ k ] = log ( (w->matrix [ k ] / sum_column [ i ]) / 0.25);
	    /*printf(" %f\n",w->matrix[k]);*/
	}
    }
    xfree ( sum_column );
    return 0;
}



int get_WtmatrixCounts ( Motif_specs *v, int counts[] ) {

  int motif, length, seq_pos;

  for ( motif = 0; motif < v->number; motif++ ) {
    length = strlen ( v->motif_spec_ptr[motif].seq );
    for ( seq_pos = 0; seq_pos < length; seq_pos++ ) {

      counts [ char_lookup [ v->motif_spec_ptr[motif].seq [ seq_pos ]] * length + seq_pos ]++;
    }
  }
  return 0;
}


int write_Wtmatrix ( WeightMatrixCounts *c, WtmatrixSpec *w, FILE *fp_o ) {

  char title[]="title";
  int NUM_COLUMNS = 20;
  int block, n_blocks, row, element;
  int col, last_col, i, depth, total;
  char bases[]={"acgtn"}, acids[]={"ACDEFGHIKLMNPQRSTVWY*-"};
  char *char_ptr;

  /* note that this routine fills the elements for unknown chars
   * with the total for the column before writing the file
  */


/* deal with charsetsize */
  if ( w->depth == 5 ) {
    depth = 5;
    char_ptr = bases;
  }
  else {
    depth = 22;
    char_ptr = acids;
  }

  /* put the column totals in the unused unknown char elements */

  /* loop for each column of the matrix */
  for ( col=0; col<w->length; col++ ) {
    total = 0;
    /* get the counts for each column */
    for ( row=0; row<w->depth-1; row++ ) {
      element = row * (w->length) + col;
      total += c->counts [ element ];
    }
    element = row * (w->length) + col;
    /* set unknown character to total */
    c->counts [ element ] = total;
  }

  /*  title */
  if ( ( fprintf (fp_o, "%s\n", title) <= 0 ) ) return -1;

  /* matrix length, mark_pos, min and max */

  if ( fprintf ( fp_o, "%d %d %f %f\n", w->length, w->mark_pos, w->min, w->max )
		 <= 0 ) return -1;

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
    last_col=MIN(w->length,((block+1)*NUM_COLUMNS));
    fprintf ( fp_o, "P" );
    for ( i=block*NUM_COLUMNS; i<last_col; i++ ) fprintf ( fp_o, "%6d", i);
    fprintf( fp_o,"\n");

    row = depth - 1;
    fprintf( fp_o,"%c",char_ptr[row]);
    for ( col = block*NUM_COLUMNS, element = w->length * row + block * NUM_COLUMNS; col < last_col; col++, element++ ) {
      fprintf (fp_o, "%6d", c->counts[element]);
    }
    fprintf( fp_o,"\n");

    for ( row=0;row<depth-1; row++ ) {
      fprintf( fp_o,"%c",char_ptr[row]);
      for ( col = block*NUM_COLUMNS, element = w->length * row + block * NUM_COLUMNS; col < last_col; col++, element++ ) {
	fprintf (fp_o, "%6d", c->counts[element]);
      }
      fprintf( fp_o,"\n");
    }
  }
  fclose ( fp_o );
  return 0;
}


int apply_Wtmatrix ( Motif_specs *v, WtmatrixSpec *w ) {

  int motif, length, seq_pos;
  double score, min_score, max_score;

  for ( motif = 0, min_score = 99999.0, max_score = -99999.0; motif < v->number; motif++ ) {
    length = strlen ( v->motif_spec_ptr[motif].seq );
    for ( seq_pos = 0, score = 0.0; seq_pos < length; seq_pos++ ) {
      score += w->matrix [ char_lookup [ v->motif_spec_ptr[motif].seq [ seq_pos ]] * w->length + seq_pos ];
      /*printf("%c %f\n",v->motif_spec_ptr[motif].seq[seq_pos],score);*/
    }
    v->motif_spec_ptr[motif].score = score;
    min_score = MIN(min_score,score);
    max_score = MAX(max_score,score);
  }
  w->min = min_score;
  w->max = max_score;
  return 0;
}

WeightMatrixCounts *read_weight_matrix ( FILE *file_p, int char_set ) {

  char line[120],char_type[2];
  int NUM_COLUMNS = 20;
  int block, n_blocks, row, char_index, element;
  int wl, wmp, depth;
  double wmn, wmx;
  WeightMatrixCounts *w = NULL;

/* deal with charsetsize */
  if ( 5 == char_set ) {
      depth = 4;
  } else {
      depth = 22;
  }

  /* read and discard title */
  if ( ( fscanf (file_p, "%[^\n]\n", line) == 0 ) ) return NULL;

  /* printf("%s\n",line); */

  /* read the matrix length, mark_pos, min and max */

  if ( fscanf ( file_p, "%d %d %lf %lf\n", &wl, &wmp,
		&wmn, &wmx ) != 4 ) return NULL;

  /* sanity check */

  if ( wl <= 0 ) return NULL;

  if (NULL == (w = initWeightMatrixCounts ( wl, depth+1 ))) return NULL;
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
	  /*printf(" r %d e %d\n",row,element);*/
	  while ( fscanf (file_p, " %d", &(w->counts[element++]))>0);

    }
  }
  return w;
}

int do_it ( Motif_specs *v, FILE *fp_o, FILE *fp_w ) {

  WtmatrixSpec *w=NULL;
  WeightMatrixCounts *c;
  int i;


  /* if we have an input weight matrix read it so we can use it to
   * score the aligned sequences
   * Otherwise we calculate the counts from the aligned sequences
   * In both cases we store the counts in c ready to calc the weights
   */

  if ( fp_w ) {

    if (  (c = read_weight_matrix ( fp_w, v->depth )) == NULL ) {
      fclose ( fp_w );
      free_WeightMatrixCounts ( c );
      return -1;

    }
    fclose ( fp_w );
  }
  else {

    if (NULL == (c = initWeightMatrixCounts ( v->length, v->depth ))) return -1;

    i = get_WtmatrixCounts ( v, c->counts );
  }

  /* initialise the weight matrix */

  if ( (w = init_Wtmatrix ( c )) == NULL ) {
    free_WeightMatrixCounts ( c );
    return -1;
  }

  /* get the weights from the counts */

  if ( get_wt_weights ( c->counts, w ) ) {
    free_WeightMatrixCounts ( c );
    return -1;
  }

  /* apply the weights to the input sequence alignments */

  i =  apply_Wtmatrix ( v, w );

  /* if we are not creating an output weight matrix, list the scores */

  if ( ! fp_o ) {
    for(i=0;i<v->number;i++) {
      printf("%-*s ",(int)v->name_length,v->motif_spec_ptr[i].name);
      printf("%s ",v->motif_spec_ptr[i].seq);
      printf("%d ",i);
      printf("%f ",v->motif_spec_ptr[i].score);
      if ( v->motif_spec_ptr[i].comment ) printf("%s",v->motif_spec_ptr[i].comment);
      printf("\n");
    }
    printf("min score %f max_score %f\n",w->min, w->max);
  }

  /* if the user has not specified the cutoffs, set them to the observed range */

  if ( v->min != UNSPECIFIED ) w->min = v->min;
  if ( v->max != UNSPECIFIED ) w->max = v->max;
  w->mark_pos = v->mark_pos;

  /* if required, write the output weight matrix */

  if ( fp_o ) i = write_Wtmatrix ( c, w, fp_o );
  free_WeightMatrixCounts ( c );
  
  if ( w ) free_WtmatrixSpec ( w );
  return 0;
}


extern DLL_IMPORT char *optarg;
extern DLL_IMPORT int optind;

void usage(void) {
    fprintf(stderr, 
	    "Usage: make_weights [options] input_file\n"
	    "Where options are:\n"
	    "    [-w input weights filename]      [-o output filename]\n"
	    "    [-c min score]                   [-C max score]\n"
	    "    [-m mark position]               [-v version]\n");
    exit(1);
}


int main(int argc, char **argv) {
    int i, c;
    char *fn_i, *fn_o, *fn_w;
    FILE *fp_i, *fp_o, *fp_w;
    Motif_specs *v;
    int mark_pos = -1;
    double min_score, max_score;
    fn_i = fn_o = fn_w = NULL;
    fp_i = fp_o = fp_w = NULL;

    min_score = max_score = UNSPECIFIED;
    set_dna_lookup();
    set_char_set(1); /* FIXME DNA*/

    while ((c = getopt(argc, argv, "o:w:m:c:C:vp")) != -1) {
	switch (c) {
	case 'm':
	    mark_pos = atoi(optarg);
	    break;
	case 'v':
	    printf("make_weights: version 1.0\n");
	    return 0;
	case 'c':
	    min_score = atof(optarg);
	    break;
	case 'C':
	    max_score = atof(optarg);
	    break;
	    /*
	case 'p':
	    char_set = PROTEIN;
	    break;
	    */
	case 'o':
	    fn_o = optarg;
	    break;
	case 'w':
	    fn_w = optarg;
	    break;
	default:
	    usage();
	}
    }

    if (optind == argc)
	    usage();

    fn_i = argv[optind];

    fp_i = fopen(fn_i, "r");
    if (fp_i == NULL ) {
	fprintf(stderr, "Failed to open input file\n");
	return -1;
    }
    if ( !(v = get_motif_info(fp_i,5)) ) return -1; /* DNA! */
    fclose ( fp_i);

    if ( fn_o ) {
      fp_o = fopen(fn_o, "w");
      if (fp_o == NULL ) {
	fprintf(stderr, "Failed to open output file\n");
	return -1;
      }
    }

    if ( fn_w ) {
      fp_w = fopen(fn_w, "r");
      if (fp_w == NULL ) {
	fprintf(stderr, "Failed to open input weights file\n");
	return -1;
      }
    }

    if ( mark_pos < 0 ) mark_pos = 0;
    v->mark_pos = mark_pos;
    v->min = min_score;
    v->max = max_score;
    i = do_it ( v, fp_o, fp_w );
    return i;
}
