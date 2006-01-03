#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "os.h"
#include "align_lib.h"
#include "misc.h"
#include "dna_utils.h"
#include "read_matrix.h"

/* affine_align aligns a pair of sequences.
 * int affine_align(OVERLAP *overlap,
	     ALIGN_PARAMS *params) {

 * Input: OVERLAP contains pointers to the two sequences
 *        ALIGN_PARAMS contains alignment parameters including
 *                     the band size and the fist row and column
 *                     to examine. It also allows the setting
 *                     for how gaps at the edges should be treated
 *                     and whether to produce a full length
 *                     alignment or the best overlap. Finally it
 *                     allows the choice of producing aligned
 *                     sequences, and/or edit buffers. New pads
 *                     can be distinguished from old using the
 *                     NEW_PAD_SYM and OLD_PAD_SYM variables.
 *
 *        W128 is a global score matrix
#define EDGE_GAPS_COUNT   1
#define EDGE_GAPS_ZERO    2
#define FULL_LENGTH_TRACE 4
#define BEST_EDGE_TRACE   8
 *
 * Output: OVERLAP nows contains alignments and/or edit buffers S1, S2
 * which can be used to construct the two aligned sequences.

 * Before calling affine_align:
 * 1. set the score matrix which defines the scores for each
 *    pair of characters
 * 2. initialise OVERLAP and set the pointers to the pair of sequences
 *    affine_align will allocate the required memory for the output.
 * 3. create and setup ALIGN_PARAMS
 * 4. run affine_align
 * After calling affine_align:
 * 5. use output
 * 6. OVERLAP and ALIGN_PARAMS can now be destroyed
 *
 *  setting gap modes
 *  -----------------
 * edge_mode
 *     1      gaps at left edge scored
 *     2      gaps at left edge not scored
 *     4      gaps at right edge scored ie we traceback from bottom 
 *            right corner
 *     8      gaps at right edge not scored ie we find the best
 *            score along both edges
 *
 *  setting output modes
 *  --------------------
 *    job
 *     1      return aligned sequences with new pad symbols replaced by old
 *     2      return edit buffers
 *     4      return edit buffer (not implemented)
 *     8      return aligned sequences with new pad symbols
 *
 *  example code
 *  ------------

    ALIGN_PARAMS *params;
    OVERLAP	*overlap;

    ierr = set_alignment_matrix("/home5/pubseq/share/tables/nuc_matrix",
				"ACGT");
    if(ierr) return -1;

    if (NULL == (overlap = create_overlap())) return -1;

    init_overlap (overlap, seq1, seq2, seq1_len, seq2_len);
   
    if (NULL == (params = create_align_params())) return -1;

    if (set_align_params (params, band, gap_open, gap_extend, edge_mode, job,
		      seq1_start, seq2_start, '*','!',0)) {
		      destroy_overlap (overlap);
		      destroy_alignment_params (params);
		      return -1;
		      };

    ierr = affine_align(overlap, params);

    print_overlap(overlap,stdout);

    destroy_overlap (overlap);
    destroy_alignment_params (params);
 */

#define BYTE_ACROSS	1
#define BYTE_DOWN	2
#define BYTE_DIAGONAL	3


static int malign_lookup[256];

#ifdef LOCAL_COPY_OF_W128

static int W128[128][128];


/*
 * Initialise 128x128 weight matrix from our score matrix
 */
void init_W128(int **matrix,
	       char *order, int unknown) {
    int i, j, i_end;
    unsigned char ci, cj;

    for (i = 0; i < 128; i++) {
        for (j = 0; j < 128; j++) {
	    W128[i][j] = unknown;
        }
    }

    i_end = strlen(order);
    for (i = 0; i < i_end; i++) {
	ci = order[i];
	for (j = 0; j < i_end; j++) {
	    cj = order[j];
	    W128[        ci ][        cj ] = matrix[i][j];
	    W128[tolower(ci)][        cj ] = matrix[i][j];
	    W128[        ci ][tolower(cj)] = matrix[i][j];
	    W128[tolower(ci)][tolower(cj)] = matrix[i][j];
	}
    }

}
#else
extern int W128[128][128];
#endif

void destroy_malign_counts(int **matrix, int length, int depth,
			   int *orig_matrix, int orig_length);
void scale_malign_scores(MALIGN *malign, int start, int end);

/*
 * Initialise malign matrix from our score matrix
 */

void init_malign_matrix(MALIGN *malign) {

  int i, j, ii,iii,jj,jjj;

    for(i=0;i<malign->charset_size;i++) {
      for(j=0;j<malign->charset_size;j++) {
	malign->matrix[i][j] = 0;
      }
    }
    for(i=0;i<malign->charset_size;i++) {
      ii = malign->charset[i];
      iii = malign_lookup[ii];
      for(j=0;j<malign->charset_size;j++) {
	jj = malign->charset[j];
	jjj = malign_lookup[jj];
	malign->matrix[jjj][iii] = W128[jj][ii];
      }
    }
}

int set_alignment_matrix ( char *fn, char *base_order ) {

    int **input_matrix;
    int unknown;
    int i, j, len;
    
    if ( input_matrix = create_matrix(fn, base_order)) {
	len = strlen(base_order);
	/*
	putchar(' ');
	for (j = 0; j < len; j++) {
	    printf("   %c", base_order[j]);
	}
	putchar('\n');
	for (j = 0; j < len; j++) {
	    printf("%c ", base_order[j]);
	    for (i = 0; i < len; i++) {
		printf("%3d ", input_matrix[i][j]);
	    }
	    putchar('\n');
	}
	*/
	unknown = 1000;
	for (j = 0; j < len; j++) {
	    for (i = 0; i < len; i++) {
		unknown = min(unknown,input_matrix[i][j]);
	    }
	}

	init_W128(input_matrix,base_order,unknown);
	/*
	len = 128;
	putchar('\n');
	for (j = 0; j < len; j++) {
	    for (i = 0; i < len; i++) {
		printf("%3d ", W128[i][j]);
	    }
	    putchar('\n');
	}
	*/
	free_matrix(input_matrix,base_order);
	return 0;
    }
    verror(ERR_WARN, "set_alignment_matrix", "matrix file not found");
    free_matrix(input_matrix,base_order);
    return -1;
}

int set_malign_charset(MALIGN *malign, char *charset) {
    malign->charset = strdup(charset);
    return 0;
}

void print_malign_matrix(MALIGN *malign) {
  int i,j;
  for(j=0;j<malign->charset_size;j++) {
    for(i=0;i<malign->charset_size;i++) {
      printf(" %d ",malign->matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

int set_band(int seq1_len, int seq2_len) {
    return MAX(20,(MIN(seq1_len,seq2_len)*0.2));
}


SEG *create_seg(void) {
    SEG *seg;

    if(NULL == (seg = (SEG *) xmalloc(sizeof(SEG)))) {
	verror(ERR_WARN, "create_seg", "xmalloc failed");
	return NULL;
    }

    seg->seq = NULL;
    return seg;
}

void destroy_seg (SEG *seg) {
  if ( seg ) {
    if ( seg->seq ) xfree ( seg->seq );
    xfree ( seg );
  }
}

MSEG *create_mseg(void) {
    MSEG *mseg;

    if(NULL == (mseg = (MSEG *) xmalloc(sizeof(MSEG)))) {
	verror(ERR_WARN, "create_seg", "xmalloc failed");
	return NULL;
    }

    mseg->seq = NULL;
    mseg->length = 0;
    mseg->offset = 0;
    return mseg;
}

void destroy_mseg (MSEG *mseg) {
  if ( mseg ) {
    if ( mseg->seq ) xfree ( mseg->seq );
    xfree ( mseg );
  }
}

void init_mseg (MSEG *mseg, char *seq, int length, int offset) {
    mseg->seq = seq;
    mseg->length = length;
    mseg->offset = offset;
}

CONTIGL *create_contig_link(void) {
  CONTIGL *contigl;
    if(NULL == (contigl = (CONTIGL *) xmalloc(sizeof(CONTIGL)))) {
	verror(ERR_WARN, "create_contigl", "xmalloc failed");
	return NULL;
    }
    contigl->next = NULL;
    contigl->id = 0;
    return contigl;
}

void destroy_contig_link(CONTIGL *contigl, int do_mesg) {
    if (!contigl)
	return;

    if (contigl->mseg)
	destroy_mseg(contigl->mseg);

    xfree(contigl);
    return;
}

MALIGN *create_malign(void) {
    MALIGN *malign;

    if(NULL == (malign = (MALIGN *) xmalloc(sizeof(MALIGN)))) {
	verror(ERR_WARN, "create_malign", "xmalloc failed");
	return NULL;
    }

    malign->contigl = NULL;
    malign->consensus = NULL;
    malign->orig_pos = NULL;
    malign->counts = NULL;
    malign->scores = NULL;
    malign->matrix = NULL;
    malign->charset_size = 6; /*  a,c,g,t,*,n */
    return malign;
}

void free_malign (MALIGN *malign) {
  if ( malign ) {
    if ( malign->counts ) destroy_malign_counts( malign->counts,
						 malign->length,
						 malign->charset_size,
						 malign->orig_counts,
						 malign->orig_length);
    if ( malign->scores ) destroy_malign_counts( malign->scores,
						 malign->length,
						 malign->charset_size,
						 malign->orig_scores,
						 malign->orig_length);
    if ( malign->matrix ) destroy_malign_counts( malign->matrix,
						 malign->charset_size,
						 malign->charset_size,
						 NULL, 0);
    if (malign->consensus) xfree(malign->consensus);
    if (malign->orig_pos) xfree(malign->orig_pos);
    if (malign->charset) xfree(malign->charset);
  }
  malign->contigl = NULL;
  malign->consensus = NULL;
  malign->orig_pos = NULL;
  malign->counts = NULL;
  malign->scores = NULL;
}

void destroy_malign (MALIGN *malign, int contig_links_too) {
    if ( malign ) {
	if (contig_links_too && malign->contigl) {
	    CONTIGL *cl, *next;
	    for (cl = malign->contigl; cl; cl = next) {
		next = cl->next;
		destroy_contig_link(cl, 1);
	    }
	}

	free_malign(malign);
	xfree (malign);
    }
}

void set_malign_lookup(int charset_size) {

    int i;

    for (i=0;i<256;i++) malign_lookup[i] = charset_size;

    malign_lookup['a'] = 0;
    malign_lookup['c'] = 1;
    malign_lookup['g'] = 2;
    malign_lookup['t'] = 3;
    malign_lookup['A'] = 0;
    malign_lookup['C'] = 1;
    malign_lookup['G'] = 2;
    malign_lookup['T'] = 3;
    malign_lookup['U'] = 3;
    malign_lookup['u'] = 3;
    malign_lookup['*'] = 4;
}

/*
 * WARNING: This is deeply horrid code and also technically illegal C due
 * to pointer arithmetic outside of the bounds of an array of single allocated
 * object. However practically speaking it'll work OK.
 *
 * The original malign->scores and malign->counts matrices are allocated as
 * an array of rows (int **) and a single array of int * to which all columns
 * are addresses within. This reduces malloc overhead dramatically and probably
 * increases cache hit rates.
 *
 * However for speed we then subsequently insert additional rows and allocate
 * their own private blocks of memory (see malign_insert_scores).
 * Consequently when deallocating some matrix[i]'s will be to malloced
 * blocks but most will be just pointers into various locations within
 * the original single large block. We use orig/orig_length to point to the
 * original single block (incase matrix[0] was inserted to).
 *
 * The nasties come from how to distinguish one type of matrix[i] from another.
 * My cheat here is simply to see if it is in the range of the original block.
 * If not then it's obviously been allocated later, so free it.
 */
void destroy_malign_counts(int **matrix, int length, int depth,
			   int *orig, int orig_length) {
    int i;
    if (orig) {
	for (i = 0; i < length; i++) {
	    if (matrix[i] < orig ||
		matrix[i] > &orig[orig_length * depth])
		free(matrix[i]);
	}
	free(orig);
    } else {
	free(matrix[0]);
    }
    free(matrix);
}

int **create_malign_counts(int length, int depth) {
  /* length is contig length, depth charset_size */
  int **counts;
  int i;

  counts = (int **)calloc(length,sizeof(int *));
  counts[0] = (int *)calloc(depth * length, sizeof(int));
  for(i=1;i<length;i++) {
      counts[i] = counts[0] + depth*i;
  }
  return counts;
}

/*
 * Inserts size columns into the malign->score array at position pos.
 * The columns are initialised to zeros.
 */
void malign_insert_scores(MALIGN *malign, int pos, int size) {
    int i;

    /* -ve size indicates a deletion. Can it ever occur? */

    /* printf("\nConsensus insert at %d + %d\n", pos, size); */

    if (pos >= malign->length) {
	size += pos - malign->length + 1;
	pos = malign->length-1;
    }

    /* Shuffle along counts */
    malign->counts = (int **)realloc(malign->counts,
				     (malign->length+size) * sizeof(int *));
    memmove(&malign->counts[pos+size], &malign->counts[pos],
	    sizeof(malign->counts[0]) * (malign->length - pos));
    for (i = pos; i < pos + size; i++) {
	malign->counts[i] = (int *)calloc(malign->charset_size, sizeof(int));
    }

    /* Shuffle along scores */
    malign->scores = (int **)realloc(malign->scores,
				     (malign->length+size) * sizeof(int *));
    memmove(&malign->scores[pos+size], &malign->scores[pos],
	    sizeof(malign->scores[0]) * (malign->length - pos));
    for (i = pos; i < pos + size; i++) {
	malign->scores[i] = (int *)calloc(malign->charset_size, sizeof(int));
    }

    /* Shuffle along consensus */
    malign->consensus = (char *)realloc(malign->consensus,
					malign->length+size);
    memmove(&malign->consensus[pos+size], &malign->consensus[pos],
	    malign->length - pos);
    malign->orig_pos = (int *)realloc(malign->orig_pos,
				      sizeof(int) * (malign->length+size));
    memmove(&malign->orig_pos[pos+size], &malign->orig_pos[pos],
	    sizeof(int) * (malign->length - pos));
    for (i = pos; i < pos + size; i++) {
	malign->consensus[i] = '-';
	malign->orig_pos[i] = 0;
    }

    malign->length += size;
}

/*
 * Removes 'del' from the contigl list. The previous contigl will be previous
 * or NULL if del is the first. (Knowing this info speeds up this function
 * so it is mandatory.)
 */
void malign_remove_contigl(MALIGN *malign, CONTIGL *previous, CONTIGL *del) {
    int i, j;
    int start, end;
    char *seq;

    start = del->mseg->offset;
    end = del->mseg->offset + del->mseg->length-1;
    seq = del->mseg->seq;

    if (previous)
	previous->next = del->next;
    else
	malign->contigl = del->next;

    for (j = start, i = 0; j <= end; i++, j++) {
	int c = malign_lookup[seq[i]];
	malign->counts[j][c]--;
    }

    get_malign_consensus(malign, start, end);
    scale_malign_scores(malign, start, end);
}

void malign_add_contigl(MALIGN *malign, CONTIGL *previous, CONTIGL *add) {
    int i, j;
    int start, end;

    start = add->mseg->offset;
    end = add->mseg->offset + add->mseg->length-1;

    if (previous) {
	add->next = previous->next;
	previous->next = add;
    } else {
	add->next = malign->contigl;
	malign->contigl = add;
    }

    for (j = start, i = 0; i < add->mseg->length; i++, j++) {
	int c = malign_lookup[add->mseg->seq[i]];
	malign->counts[j][c]++;
    }
    
    get_malign_consensus(malign, start, end);
    scale_malign_scores(malign, start, end);
}

void malign_recalc_scores(MALIGN *malign, int start, int end) {
    /* Recompute scores */
    get_malign_counts(malign, start, end);

    /* Recompute consensus */
    get_malign_consensus(malign, start, end);

    /* Scale scores */
    scale_malign_scores(malign, start, end);
}
 
/*
 * Fill out malign->scores[] from start to end inclusively.
 *
 * FIXME: This gets called a lot with various start and end ranges.
 * We need a way to track which sequences overlap specific regions so
 * we do not have to loop through the entire contig chaining along the
 * linked list each time. Over half of all the cpu time on a small 454
 * project (14000 reads) is spent in this function.
 *
 * We wouldn't need this for every possible value of start, just maybe
 * start/100 or something, so that we can limit the search to be within
 * a 100 bp of start rather than from the contig start. NB. How to keep this
 * up to date when moving or lengthening sequences?
 */
void get_malign_counts (MALIGN *malign, int start, int end) {
    CONTIGL *t;
    int i,j,k,l;

    /* Zero */
    for (i=start; i<=end; i++) {
	for (j=0; j<malign->charset_size; j++) {
	    malign->counts[i][j] = 0;
	}
    }

    /* Accumulate */
    for (t = malign->contigl; t && t->mseg->offset <= end; t=t->next) {
	if (t->mseg->offset + t->mseg->length-1 < start)
	    continue;
	for(j=0,k=t->mseg->offset;j<t->mseg->length;j++,k++) {
	    if (k < start)
		continue;
	    if (k > end)
		break;
	    l = malign_lookup[(int)(t->mseg->seq[j])];
	    malign->counts[k][l]++;
	}
    }
}

void get_malign_consensus(MALIGN *malign, int start, int end) {
  int i,j,k;
  char *s;

  /* must be run after get_malign_scores */

  s = malign->consensus;
  k = malign->charset_size;
  for(i=start;i<=end;i++) {
      int top = 0;
      s[i] = '-';
      for(j = 0; j < k; j++) {
	  /* Most common base call becomes the consensus */
	  int sc = malign->counts[i][j];
	  if (top < sc) {
	      top = sc;
	      s[i] = malign->charset[j];
	  }
      }
  }
}

void print_malign_scores(MALIGN *malign) {
  int i,j;
  for(i=0;i<malign->length;i++) {
      for(j=0;j<malign->charset_size;j++) {
	  printf(" %+4d ",malign->scores[i][j]);
      }
      printf("\n");
  }
  printf("\n");
}

void scale_malign_scores(MALIGN *malign, int start, int end) {
  int i,j,k;
  /* in the alignment routine all these values are added:
   * ie score is score + malign->scores[i][j];
   * even when scores[i][j] is a gap penalty.
   * to calculate the score j = seq2[row-1]
   */

  /* An update on Rodger's where the score incorporates positive from match
   * with negative from mismatches
   */
#if 0
  k = malign->matrix[0][1];
  for(i=start;i<=end;i++) {
      for (l=j=0;j<malign->charset_size;j++) {
	  l += malign->counts[i][j];
      }
      for(j=0;j<malign->charset_size;j++) {
	  if (malign->counts[i][j]) {
	      malign->scores[i][j] =
		  malign->counts[i][j] * malign->matrix[j][j] +
		  k*(l - malign->counts[i][j]);
	  } else {
	      malign->scores[i][j] = k*l;
	  }
      }
  }
#endif


  /* Simple unit based scores with 0 = perfect 1 = wrong. Based on ReAligner */
  /* Scale by 100 to fit in integers */
#if 1
  k = malign->matrix[0][1];
  for(i=start;i<=end;i++) {
      int t = 0;
      double s = 0;
      int max = 0;
      for (j = 0; j < malign->charset_size; j++) {
	  t += malign->counts[i][j];
	  if (max < malign->counts[i][j])
	      max = malign->counts[i][j];
      }
      if (t) {
	  for(j=0;j<malign->charset_size;j++) {
	      s = 0.5 * (1 - malign->counts[i][j] / (double)t);
	      if (malign->counts[i][j] != max)
		  s += 0.5;
	      malign->scores[i][j] = s * 100;
	  }
      } else {
	  for(j=0;j<malign->charset_size;j++) {
	      malign->scores[i][j] = 0;
	  }
      }
  }
#endif
}

void print_contig_links (CONTIGL *contigl) {
  CONTIGL *t;
  t = contigl;
  while(t) {
    printf("%d %d %s\n",t->mseg->length,t->mseg->offset,t->mseg->seq);
    t = t->next;
  }
}

int contigl_length (CONTIGL *contigl) {
  CONTIGL *t;
  int length;

  length = 0;
  t = contigl;
  while(t) {
    length = MAX(length,(t->mseg->length + t->mseg->offset));
    t = t->next;
  }
  return length;
}

int contigl_elements (CONTIGL *contigl) {
  CONTIGL *t;
  int elements;

  elements = 0;
  t = contigl;
  while(t) {
    elements++;
    t = t->next;
  }
  return elements;
}

void print_malign_seqs(MALIGN *malign) {
  int i;
  CONTIGL *t;
  i = 0;
  t = malign->contigl;
  while(t) {
      /*
    printf("%d %d %*.s %s\n",
	   i++,
	   t->mseg->length,
	   t->mseg->offset,
	   "                       ",
	   t->mseg->seq);
	   */
    t = t->next;
  }
}

MALIGN *contigl_to_malign(CONTIGL *contigl_in, int gap_open, int gap_extend) {
  MALIGN *malign;
  CONTIGL *contigl;
  int i;

  /* get it started */

  if(!(malign = create_malign())) {
    printf("scream contig_to_malign\n");
    return NULL;
  }
  contigl = contigl_in;
  malign->contigl = contigl;
  print_malign_seqs(malign);
  i = set_malign_charset(malign,"acgt*n");

  /*printf("charset %s\n",malign->charset);*/

  /* create the malign matrix from the external score matrix (W128) */

  malign->matrix = create_malign_counts(malign->charset_size,malign->charset_size);
  init_malign_matrix(malign);
  /*print_malign_matrix(malign);*/

  malign->length = contigl_length(contigl);

  /* get the initial scores from the alignment */

  malign->counts = create_malign_counts(malign->length,malign->charset_size);
  malign->scores = create_malign_counts(malign->length,malign->charset_size);
  malign->orig_counts = malign->counts[0];
  malign->orig_scores = malign->scores[0];
  malign->orig_length = malign->length;
  get_malign_counts(malign, 0, malign->length-1);


  //print_malign_scores(malign);

  /* make a 100% consensus for the alignment */

  malign->consensus = (char *) malloc(malign->length);
  malign->orig_pos = (int *) malloc(malign->length * sizeof(int));
  for (i = 0; i < malign->length; i++) {
      malign->orig_pos[i] = i+1;
  }
  get_malign_consensus(malign, 0, malign->length-1);
  /* printf("      %s\n",malign->consensus); */

  /* scale the scores with the gap penalties and the external score matrix */

  malign->gap_open = gap_open;
  malign->gap_extend = gap_extend;
  scale_malign_scores(malign, 0, malign->length-1);
  //print_malign_scores(malign);

  /* FIXME put in error checking for failed mallocs etc */
  return malign;
}

MOVERLAP *create_moverlap(void) {
    MOVERLAP	*moverlap;

    if(NULL == (moverlap = (MOVERLAP *) xmalloc(sizeof(MOVERLAP)))) {
	verror(ERR_WARN, "create_moverlap", "xmalloc failed");
	return NULL;
    }


    moverlap->S = NULL;
    moverlap->S1 = NULL;
    moverlap->S2 = NULL;
    moverlap->malign = NULL;
    moverlap->seq2 = NULL;
    moverlap->malign_out = NULL;
    moverlap->seq1_out = NULL;
    moverlap->seq2_out = NULL;

    return moverlap;
}

void init_moverlap (MOVERLAP *moverlap, MALIGN *malign, char *seq2, int malign_len,
		   int seq2_len) {

    moverlap->malign = malign;
    moverlap->seq2 = seq2;
    moverlap->malign_len = malign_len;
    moverlap->seq2_len = seq2_len;
    moverlap->S1  = NULL;
    moverlap->S2 = NULL;
    moverlap->S  = NULL;
    moverlap->malign_out = NULL;
    moverlap->seq1_out = NULL;
    moverlap->seq2_out = NULL;
}
   
void destroy_moverlap (MOVERLAP *moverlap) {
  if ( moverlap ) {
    if ( moverlap->S1 ) xfree ( moverlap->S1 );
    if ( moverlap->S2 ) xfree ( moverlap->S2 );
    if ( moverlap->S ) xfree ( moverlap->S );
    if ( moverlap->malign_out ) xfree ( moverlap->malign_out );
    if ( moverlap->seq1_out ) xfree ( moverlap->seq1_out );
    if ( moverlap->seq2_out ) xfree ( moverlap->seq2_out );
    xfree ( moverlap );
  }
}

void free_moverlap (MOVERLAP *moverlap) {
  if ( moverlap ) {
    if ( moverlap->S1 ) xfree ( moverlap->S1 );
    if ( moverlap->S2 ) xfree ( moverlap->S2 );
    if ( moverlap->S ) xfree ( moverlap->S );
    if ( moverlap->malign_out ) xfree ( moverlap->malign_out );
    if ( moverlap->seq1_out ) xfree ( moverlap->seq1_out );
    if ( moverlap->seq2_out ) xfree ( moverlap->seq2_out );
    moverlap->S1  = NULL;
    moverlap->S2 = NULL;
    moverlap->S  = NULL;
    moverlap->malign_out = NULL;
    moverlap->seq1_out = NULL;
    moverlap->seq2_out = NULL;  }
}

OVERLAP *create_overlap(void) {
    OVERLAP	*overlap;

    if(NULL == (overlap = (OVERLAP *) xmalloc(sizeof(OVERLAP)))) {
	verror(ERR_WARN, "create_overlap", "xmalloc failed");
	return NULL;
    }


    overlap->S = NULL;
    overlap->S1 = NULL;
    overlap->S2 = NULL;
    overlap->seq1 = NULL;
    overlap->seq2 = NULL;
    overlap->seq1_out = NULL;
    overlap->seq2_out = NULL;

    return overlap;
}

void init_overlap (OVERLAP *overlap, char *seq1, char *seq2, int seq1_len,
		   int seq2_len) {

    overlap->seq1 = seq1;
    overlap->seq2 = seq2;
    overlap->seq1_len = seq1_len;
    overlap->seq2_len = seq2_len;
    overlap->S1  = NULL;
    overlap->S2 = NULL;
    overlap->S  = NULL;
    overlap->seq1_out = NULL;
    overlap->seq2_out = NULL;
    overlap->percent = 0;
}
   
void destroy_overlap (OVERLAP *overlap) {
   if ( overlap ) {
    if ( overlap->S1 ) xfree ( overlap->S1 );
    if ( overlap->S2 ) xfree ( overlap->S2 );
    if ( overlap->S ) xfree ( overlap->S );
    if ( overlap->seq1_out ) xfree ( overlap->seq1_out );
    if ( overlap->seq2_out ) xfree ( overlap->seq2_out );
    xfree ( overlap );
   }
}

void free_overlap (OVERLAP *overlap) {
  if ( overlap ) {
    if ( overlap->S1 ) xfree ( overlap->S1 );
    if ( overlap->S2 ) xfree ( overlap->S2 );
    if ( overlap->S ) xfree ( overlap->S );
    if ( overlap->seq1_out ) xfree ( overlap->seq1_out );
    if ( overlap->seq2_out ) xfree ( overlap->seq2_out );
    overlap->S1  = NULL;
    overlap->S2 = NULL;
    overlap->S  = NULL;
    overlap->seq1_out = NULL;
    overlap->seq2_out = NULL;
  }
}

ALIGN_PARAMS *create_align_params(void) {
    ALIGN_PARAMS *params;

    if(NULL == (params = (ALIGN_PARAMS *) xmalloc(sizeof(ALIGN_PARAMS)))) {
	verror(ERR_WARN, "create_align_params", "xmalloc failed");
	return NULL;
    }
    params->gap_open = 12;
    params->gap_extend = 4;
    params->band = 0;
    params->first_row = 0;
    params->band_left = 0;
    params->band_right = 0;
    params->edge_mode = EDGE_GAPS_COUNT | BEST_EDGE_TRACE;
    params->job = RETURN_SEQ;
    params->new_pad_sym = '.';
    params->old_pad_sym = '*';
    
    return params;
}

void print_align_params(ALIGN_PARAMS *params) {
    printf("gap_open %d\ngap_extend %d\nband %d\nfirst_row %d\n"
	   "band_left %d\nband_right %d\nedge_mode %d\njob %d\n"
	   "new_pad_sym %c\nold_pad_sym %c\n",
    params->gap_open,
    params->gap_extend,
    params->band,
    params->first_row,
    params->band_left,
    params->band_right,
    params->edge_mode,
    params->job,
    params->new_pad_sym,
    params->old_pad_sym);
}
	   


int set_align_params (ALIGN_PARAMS *params, int band, int gap_open, 
		       int gap_extend, int edge_mode, int job, 
		       int seq1_start, int seq2_start, 
		      char old_pad_sym, char new_pad_sym,
		       int set_job) {


  int d, first_column;
  /*printf("band %d gap %d gap %d edge %d job %d s1 %d s2 %d pad %c\n",
	 band,gap_open,gap_extend,edge_mode,job,seq1_start,seq2_start,pad_sym);*/
  if( set_job ) {
    /* only set parameters related to band:
     * need band and start positions,
     * set band, first_row, band_left, band_right
     */
    params->band = band;
    params->first_row = 0;
    first_column = 0;
    params->band_left = 0;
    params->band_right = 0;
    if(params->band) {
      d = MIN(seq2_start, params->band);
      
      params->first_row = seq2_start - d;
      first_column = seq1_start - d;
      params->band_left = first_column - params->band;
      params->band_right = first_column + params->band;
    }
    return 0;
  }
  if ( job & RETURN_EDIT_BUFFER ) {
    verror(ERR_WARN, "affine_align", "unimplemented alignment job");
    return -1;
  }

  if ( job && !(job & RETURN_SEQ) && !(job & RETURN_EDIT_BUFFERS)) {
    verror(ERR_WARN, "affine_align", "unknown alignment job");
    return -1;
  }
  if(gap_open) params->gap_open = gap_open;
  if(gap_extend) params->gap_extend = gap_extend;
  if(edge_mode) params->edge_mode = edge_mode; 
  if(job) params->job = job;
  if ( old_pad_sym) params->old_pad_sym = old_pad_sym;
  if ( new_pad_sym) params->new_pad_sym = new_pad_sym;
  /*
     #define RETURN_SEQ 1
     #define RETURN_EDIT_BUFFERS 2
     #define RETURN_EDIT_BUFFER  4
     #define RETURN_NEW_PADS 8
  
    #define EDGE_GAPS_COUNT   1
    #define EDGE_GAPS_ZERO    2
    #define FULL_LENGTH_TRACE 4
    #define BEST_EDGE_TRACE   8
  */
  params->band = band;
  params->first_row = 0;
  first_column = 0;
  params->band_left = 0;
  params->band_right = 0;
  if(params->band) {
      d = MIN(seq2_start, params->band);
      
      params->first_row = seq2_start - d;
      first_column = seq1_start - d;
      params->band_left = first_column - params->band;
      params->band_right = first_column + params->band;
  }
  return 0;
}

void destroy_alignment_params (ALIGN_PARAMS *params) {
    if ( params ) xfree ( params );
}


void print_char_array ( FILE *file, char *array, int array_len) {
    
#define LINELENGTH 60
    
    int lines, i, i1, i2, j;
    array_len = MIN(60,array_len);
    lines = 1 + array_len/LINELENGTH;
    if ( (array_len % LINELENGTH) == 0 ) lines -= 1;
    for ( j=0; j <= lines; j++) {
	i1 =  j * LINELENGTH;
	i2 = MIN (i1+LINELENGTH-1,array_len-1);
	for ( i=i1; i <= i2; i++)
	    putc ( array[i], file );
	putc ( '\n', file );
    }
}


void seq_expand(char	  *seq,
		char	  *seq_align,
		int	  *seq_align_len,
		int *S,
		int       s_len,
		int       mode,
		char PAD_SYM) {
    
    /*
     * Note: S is a sequence edit buffer, not an
     * alignment edit buffer. See align_lib.h for
     * the difference.
     */
    
    /*
     * Expands the sequence in one of four slightly different ways,
     * depending on the value of mode:
     *
     * 0 = No asterisks at either end of the sequence returned.
     *     i.e. the length of overhang of any sequence with this
     *     sequence is not represented in the returned sequence.
     * 1 = Represent overhang at the left-hand end only.
     * 2 = Represent overhang at the right-hand end only.
     * 3 = Represent overhang at both ends.
     */
    
    int		i, j;
    int		s, s_start, s_end;
    int	l;


    /*for(i=0;i<s_len;i++)printf("%d %d\n",i,S[i]);*/
    
    if((0 == mode) || (1 == mode)) {
	/* Ignore right-hand overhang */
	for(s = s_len - 1; s >= 0; s--) {
	    if(S[s] > 0) {
		s_end = s + 1;
		break;
	    }
	}
    }

    else
	s_end = s_len;
    
    if((0 == mode) || (2 == mode)) {
	/* Ignore left-hand overhang */
	for(s = 0; s < s_len; s++) {
	    if(S[s] > 0) {
		s_start = s;
		break;
	    }
	}
    }
    else
	s_start = 0;

    *seq_align = '\0';
    for(i = 0, j = 0, s = s_start; s < s_end; s++) {
	l = S[s];
	if(l > 0) {
	    strncpy(seq_align + j, seq + i, l);
	    *(seq_align + j + l) = '\0';
	    j += l;
	    i += l;
	}
	else {
	    memset(seq_align + j, PAD_SYM, -l);
	    *(seq_align + j - l) = '\0';
	    j -= l;
	}
    }
    *(seq_align + j) = '\0';
    *seq_align_len = j;
}

int get_segment( OVERLAP *overlap, SEG *seg,
		 int job ) {

    int		len_align;
    
    int  seq_align_len;
    char PAD_SYM;

    PAD_SYM = '*';



    /* what do I need to return?
     *  1. for consen I need the righthand overhang only
     *  2. sequence excluding the overhangs
     *
     *  use job numbers:
     *  job 1: righthand overhang for seq1
     *  job 2: righthand overhang for seq2
     *  job 3: sequence excluding overhangs for seq1
     *  job 4: sequence excluding overhangs for seq2
     *
     *  note that for extending the ends of contigs the current method
     *  doe snot always provide the desired result. For example:
     *  consensus **tg*******                                       
     *  reading   cct*ttataaa                                       
     *  The code assumes that there is no gap at the left end - ie
     *  we expect there to be sufficient consensus to get the registration
     *  with the hidden data. Here there is little data and no match and it
     *  might have been better to return tttataaa - ie just remove the two
     *  chars that shoul dhave aligned with the consensus. It coul dhave been
     *  worse:
     *  consensus ********aa                                        
     *  reading   cctttataaa
     *  but it is not simply a matter of taking into account the number of
     *  pads at the left end of the consensus because in the first example
     *  shown above I would end up returning t*ttataaa
     *  I have left it for now but could sor tit if necessary FIXME!
     */

    /* OVERLAP must be set up before entry, including allocating the
     * memory for the segments
     * I think seq_expand only works for its option 3?
     */


    if ( job == 1 ) {
	seq_expand(overlap->seq1, seg->seq, &seq_align_len, 
		   overlap->S1, overlap->s1_len, 3, PAD_SYM);
	len_align = MAX(0,MAX(overlap->right1,overlap->right2) - overlap->right2);
	memmove ( seg->seq, seg->seq+overlap->right2+1, len_align );
	seg->seq[len_align] = '\0';
	seg->length = len_align;
	return 0;
    }
    if ( job == 2 ) {
	seq_expand(overlap->seq2, seg->seq, &seq_align_len, 
		   overlap->S2, overlap->s2_len, 3, PAD_SYM);
	len_align = MAX(0,MAX(overlap->right1,overlap->right2) - overlap->right1);
	memmove ( seg->seq, seg->seq+overlap->right1+1, len_align );
	seg->seq[len_align] = '\0';
	seg->length = len_align;
	return 0;
    }
    if ( job == 3 ) {
	seq_expand(overlap->seq1, seg->seq, &seq_align_len, 
		   overlap->S1, overlap->s1_len, 3, PAD_SYM);
	len_align =   overlap->length;
	memmove ( seg->seq, seg->seq+MAX(overlap->left1,overlap->left2), 
		 len_align );
	seg->seq[len_align] = '\0';
	seg->length = len_align;
	return 0;
    }
    if ( job == 4 ) {
	seq_expand(overlap->seq2, seg->seq, &seq_align_len, 
		   overlap->S2, overlap->s2_len, 3, PAD_SYM);
	len_align =   overlap->length;
	memmove ( seg->seq, seg->seq+MAX(overlap->left1,overlap->left2), 
		 len_align );
	seg->seq[len_align] = '\0';
	seg->length = len_align;
	return 0;
    }
    return -2;
}

int print_alignment(char	*seq1,
		    char	*seq2,
		    int		seq1_len,
		    int		seq2_len,
		    int	*S1,
		    int	*S2,
		    int		s1_len,
		    int		s2_len,
		    double	score,
		    FILE	*fpt) {
    
    char	*seq1_align, *seq2_align;
    int		seq1_align_len, seq2_align_len;
    char	temp_seq[51];
    int		i, j;
    int		max_out_width;
    int		max_seq;
    int		len_align;
    char PAD_SYM;
    PAD_SYM = '*';
    
    max_seq = seq1_len + seq2_len + 1;
    if(!(seq1_align = (char *) xmalloc(sizeof(char) * max_seq))) return -1;
    if(!(seq2_align = (char *) xmalloc(sizeof(char) * max_seq))) {
	xfree(seq1_align);
	return -1;
    }

    seq_expand(seq1, seq1_align, &seq1_align_len, S1, s1_len, 3, PAD_SYM);

    seq_expand(seq2, seq2_align, &seq2_align_len, S2, s2_len, 3, PAD_SYM);
    len_align = MAX(seq1_align_len, seq2_align_len);
    
    fprintf(fpt, "Alignment:\n");
    memset(temp_seq, '\0', 51);
    
    fprintf(fpt, "length = %d\n", len_align);
    fprintf(fpt, "score = %f\n", score);
    
    for(i = 0; i < len_align; i += 50) {
	fprintf(fpt, "\n     %10d%10d%10d%10d%10d\n", i + 10, i + 20, i + 30, i + 40, i + 50);

	max_out_width = MIN(len_align - i, 50);

	memset(temp_seq, ' ', 50);
	strncpy(temp_seq, seq1_align + i, max_out_width);
	fprintf(fpt, "     %-50s\n", temp_seq);

	memset(temp_seq, ' ', 50);
	strncpy(temp_seq, seq2_align + i, max_out_width);
	fprintf(fpt, "     %-50s\n", temp_seq);

	memset(temp_seq, ' ', 50);
	for(j = 0; (j < max_out_width) && (i + j < len_align); j++) {
	    *(temp_seq + j) = (toupper(*(seq1_align + i + j)) == toupper(*(seq2_align + i + j))) ? '+' : ' ';
	}
	fprintf(fpt, "     %-50s\n", temp_seq);
    }
    
    xfree(seq1_align);
    xfree(seq2_align);
    
    return 0;
}

int print_overlap(OVERLAP *overlap, FILE *fpt) {
    
    char	*seq1_align, *seq2_align;
    int		seq1_align_len, seq2_align_len;
    char	temp_seq[51];
    int		i, j;
    int		max_out_width;
    int		max_seq;
    int		len_align;
    int seq1_len, seq2_len, s1_len,s2_len, *S1, *S2;
    double score;
    char *seq1, *seq2;
    char PAD_SYM;

    PAD_SYM = '.';

    seq1 = overlap->seq1;
    seq2 = overlap->seq2;
    seq1_len = overlap->seq1_len;
    seq2_len = overlap->seq2_len;
    score = overlap->score;

    if( !overlap->seq1_out ) {
	S1 = overlap->S1;
	S2 = overlap->S2;
	s1_len = overlap->s1_len;
	s2_len = overlap->s2_len;

	max_seq = seq1_len + seq2_len + 1;
	if(!(seq1_align = (char *) xmalloc(sizeof(char) * max_seq))) return -1;
	if(!(seq2_align = (char *) xmalloc(sizeof(char) * max_seq))) {
	    xfree(seq1_align);
	    return -1;
	}
	seq_expand(seq1, seq1_align, &seq1_align_len, S1, s1_len, 3, PAD_SYM);

	seq_expand(seq2, seq2_align, &seq2_align_len, S2, s2_len, 3, PAD_SYM);
	len_align = MAX(seq1_align_len, seq2_align_len);
    }
    else {
	seq1_align = overlap->seq1_out;
	seq2_align = overlap->seq2_out;
	len_align  = overlap->seq_out_len;
    }
    
    fprintf(fpt, "Alignment:\n");
    memset(temp_seq, '\0', 51);
    
    fprintf(fpt, "length = %d\n", len_align);
    fprintf(fpt, "score = %f\n", score);
    
    for(i = 0; i < len_align; i += 50) {
	fprintf(fpt, "\n     %10d%10d%10d%10d%10d\n", i + 10, i + 20, i + 30, i + 40, i + 50);

	max_out_width = MIN(len_align - i, 50);

	memset(temp_seq, ' ', 50);
	strncpy(temp_seq, seq1_align + i, max_out_width);
	fprintf(fpt, "     %-50s\n", temp_seq);

	memset(temp_seq, ' ', 50);
	strncpy(temp_seq, seq2_align + i, max_out_width);
	fprintf(fpt, "     %-50s\n", temp_seq);

	memset(temp_seq, ' ', 50);
	for(j = 0; (j < max_out_width) && (i + j < len_align); j++) {
	    *(temp_seq + j) = (toupper(*(seq1_align + i + j)) == toupper(*(seq2_align + i + j))) ? '+' : ' ';
	}
	fprintf(fpt, "     %-50s\n", temp_seq);
    }
    
    if ( !overlap->seq1_out ) {
	xfree(seq1_align);
	xfree(seq2_align);
    }
    return 0;
}

void print_fasta(char *description, char *seq, FILE *fpt) {
    int		i;
    int		seq_len;
    char	temp[61];

    fprintf(fpt, ">%s\n", description);

    seq_len = strlen(seq);
    for(i = 0; i < seq_len; i += 60) {
	memset(temp, '\0', 61);
	strncpy(temp, seq + i, 60);
	fprintf(fpt, "%s\n", temp);
    }
}



int overlap_ends ( char *seq, int seq_len, char PAD_SYM, 
		   int *left, int *right ) {
  /* looking for new pads */
    /*
     *
     *left1       right1
     *    AACGTAAT**CGCT***
     *    01234567890123456
     *    ****TATTGCCGCTAAG
     *    left2      right2
     */

  int i,j;

  j = -1;
  for(i=0;i<seq_len;i++) {
    if ( seq[i] != PAD_SYM ) {
      j = i;
      break;
    }
  }
  if( j == -1 ) {
    return -1;
  }
  *left = j;
  j = -1;
  for(i=seq_len-1;i>-1;i--) {
    if ( seq[i] != PAD_SYM ) {
      j = i;
      break;
    }
  }
  if( j == -1 ) {
    return -1;
  }
  *right = j;
  return 0;
}

void old_pads_for_new ( char *seq, int seq_len, char OLD_PAD_SYM, char NEW_PAD_SYM) {
  int i;
  for(i=0;i<seq_len;i++) {
    if(seq[i]==NEW_PAD_SYM)seq[i]=OLD_PAD_SYM;
  }
}

int seq_to_overlap ( OVERLAP *overlap, char OLD_PAD_SYM, char NEW_PAD_SYM) {

  int i,j, score;

  if ( i = overlap_ends(overlap->seq1_out, overlap->seq_out_len,
			NEW_PAD_SYM, &overlap->left1, &overlap->right1)) {
      verror(ERR_WARN, "affine_align", "error parsing alignment");
      return -1;
  }
  if ( i = overlap_ends(overlap->seq2_out, overlap->seq_out_len,
			NEW_PAD_SYM, &overlap->left2, &overlap->right2)) {
      verror(ERR_WARN, "affine_align", "error parsing alignment");
      return -1;
  }
  overlap->left = MAX(overlap->left1, overlap->left2);
  overlap->right = MIN(overlap->right1, overlap->right2);
  
  
  /*
   * Work out the direction of the overlap:
   * 0 - suffix of seq1 overlaps with prefix of seq2
   * 1 - suffix of seq2 overlaps with prefix of seq1
   * 2 - seq1 contains seq2
   * 3 - seq2 contains seq1
   */
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
  
  overlap->length = overlap->right - overlap->left + 1;
  
  score = 0;
  for(i=overlap->left,j=0;i<=overlap->right;i++) {
      score -= 4;
      if(SEQ_MATCH((int) overlap->seq1_out[i], (int) overlap->seq2_out[i])) {
	  j++;
	  score+=5;
      }
      if( (overlap->seq1_out[i] == NEW_PAD_SYM) && (overlap->seq2_out[i] == OLD_PAD_SYM)) {
	  j++;
	  score+=5;
      }
  }

  if(overlap->length) {
      overlap->percent = 100.0 * j / overlap->length;
      overlap->score = score;
  }

  overlap->qual = overlap->score;
  return 0;
}

int seq_to_moverlap ( MOVERLAP *overlap, char OLD_PAD_SYM, char NEW_PAD_SYM) {

  int i,j;

  if ( i = overlap_ends(overlap->seq1_out, overlap->seq_out_len,
			NEW_PAD_SYM, &overlap->left1, &overlap->right1)) {
      verror(ERR_WARN, "affine_align", "error parsing alignment");
      return -1;
  }
  if ( i = overlap_ends(overlap->seq2_out, overlap->seq_out_len,
			NEW_PAD_SYM, &overlap->left2, &overlap->right2)) {
      verror(ERR_WARN, "affine_align", "error parsing alignment");
      return -1;
  }
  overlap->left = MAX(overlap->left1, overlap->left2);
  overlap->right = MIN(overlap->right1, overlap->right2);
  
  
  /*
   * Work out the direction of the overlap:
   * 0 - suffix of seq1 overlaps with prefix of seq2
   * 1 - suffix of seq2 overlaps with prefix of seq1
   * 2 - seq1 contains seq2
   * 3 - seq2 contains seq1
   */
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
  
  overlap->length = overlap->right - overlap->left + 1;
  
  for(i=overlap->left,j=0;i<=overlap->right;i++) {
      if(SEQ_MATCH((int) overlap->seq1_out[i], (int) overlap->seq2_out[i])) j++;
      if( (overlap->seq1_out[i] == NEW_PAD_SYM) && (overlap->seq2_out[i] == OLD_PAD_SYM)) j++;
  }

  if(overlap->length) {
      overlap->percent = 100.0 * j / overlap->length;
  }

  overlap->qual = overlap->score;
  return 0;
}

void print_overlap_struct(OVERLAP *overlap) {
    printf("overlap->left1 %d\n",overlap->left1);
    printf("overlap->right1 %d\n",overlap->right1);
    printf("overlap->left2 %d\n",overlap->left2);
    printf("overlap->right2 %d\n",overlap->right2);
    printf("overlap->left %d\n",overlap->left);
    printf("overlap->right %d\n",overlap->right);
    printf("overlap->length %d\n",overlap->length);
    printf("overlap->direction %d\n",overlap->direction);
    printf("overlap->lo %d\n",overlap->lo);
    printf("overlap->ro %d\n",overlap->ro);
    printf("overlap->percent %f\n",overlap->percent);
    printf("overlap->score %f\n",overlap->score);
    printf("overlap->qual %f\n",overlap->qual);
    if(overlap->seq1)printf("overlap->seq1 %p\n",overlap->seq1);
    if(overlap->seq2)printf("overlap->seq2 %p\n",overlap->seq2);
    if(overlap->seq1_out)printf("overlap->seq1_out %p\n",overlap->seq1_out);
    if(overlap->seq2_out)printf("overlap->seq2_out %p\n",overlap->seq2_out);
    if(overlap->S1)printf("overlap->S1 %p\n",(void *)overlap->S1);
    if(overlap->S2)printf("overlap->S2 %p\n",(void *)overlap->S2);
}

int seq_to_edit ( char *seq, int seq_len, int **S_out, int *S_len,
		 char PAD_SYM) {
  int i, j, gap;
  int *S;

  if(!(S = (int *) xmalloc(sizeof(int) * seq_len))) {
    verror(ERR_WARN, "affine_align", "malloc failed in seq_to_edit");
    return -1;
  }

  S[0] = 0;
  gap = 0;
  if ( seq[0] == PAD_SYM ) {
    gap = 1;
  }
  for ( i=0, j=0; i<seq_len; i++ ) {

    if ( gap ) {
      if ( seq[i] == PAD_SYM ) {
        S[j]--;
      }
      else {
	/* end of gap */
	j++;
        S[j] = 1;
        gap = 0;
      }
    }
    else {
      if ( seq[i] == PAD_SYM ) {
        /* start of gap */
	j++;
        S[j] = -1;
        gap = 1;
      }
      else {
        S[j]++;
      }
    }
  }
  *S_len = j+1;
  /*for(i=0;i<*S_len;i++)printf("%d\n",S[i]);*/
  *S_out = S;
  return 0;
}

int do_trace_back_bits ( unsigned char *bit_trace, char *seq1, char *seq2,
		   int seq1_len, int seq2_len, char **seq1_out_ret, char **seq2_out_ret,
		   int *seq_out_len, int b_r, int b_c, int b_e,
		   int band, int first_band_left, int first_row, 
		   int band_length, char PAD_SYM ) {

  int i, j, r, c, e;
  unsigned char trace_byte;
  int byte, nibble, band_left, max_seq;
  char *seq1_res, *seq2_res;
  char *seq1_out, *seq2_out;

  max_seq = seq1_len + seq2_len + 1;

  if(!(seq1_out = (char *) xmalloc(sizeof(char) * max_seq))) {
    verror(ERR_WARN, "affine_align", "malloc failed in do_trace_back");
    return -1;
  }
  if(!(seq2_out = (char *) xmalloc(sizeof(char) * max_seq))) {
    if ( seq1_out ) xfree (seq1_out);
    verror(ERR_WARN, "affine_align", "malloc failed in do_trace_back");
    return -1;
  }

  seq1_res = seq1_out;
  seq2_res = seq2_out;
		    
  for(i=0;i<max_seq-1;i++) {
    *(seq1_out++) = PAD_SYM;
    *(seq2_out++) = PAD_SYM;
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

    byte = e / 4;
    nibble = 2 * (e % 4);
    trace_byte = (bit_trace[byte] >> nibble) & 3;
      
    if(trace_byte == BYTE_DIAGONAL) {
      r--,c--;
      *(seq1_out--) = seq1[c];
      *(seq2_out--) = seq2[r];
    }
    else if(trace_byte == BYTE_DOWN) {
      r--;
      /*
       * Only gap the consensus vector if the single sequence is not a pad
       * itself. Otherwise we consider this to be a non-edit. Why? Because
       * we left seq2 padded and an insertion to seq1 is effectively removing
       * the pad from seq2. Not depadding seq2 means that realigning data
       * previously aligned and padded can be done with a narrow band around
       * the diagonal.
       */
      if (seq2[r] != '*') {
	  seq1_out--;
	  *(seq2_out--) = seq2[r];
      }
    }
    else {
      c--;
      seq2_out--;
      *(seq1_out--) = seq1[c];
    }
    if(band) {
      band_left = first_band_left + r - first_row;
      e = ((r - first_row + 1) * band_length) + (c - band_left + 1);
    }
    else {
      e = r * (seq1_len + 1) + c;
    }
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

  c = MAX(strlen ( seq1_res ),strlen ( seq2_res ));

  for ( i=0;i<c;i++ ) {
    if ( (seq1_res[i] != PAD_SYM) || (seq2_res[i] != PAD_SYM) ) break;
  }
  r = i;
  for ( i=0,j=r;j<c;i++,j++ ) {
    seq1_res[i] = seq1_res[j];
    seq2_res[i] = seq2_res[j];
  }
  seq1_res[i] = seq2_res[i] = '\0';
  *seq_out_len = i;
  seq1_out = seq1_res;
  seq2_out = seq2_res;
  *seq1_out_ret = seq1_out;
  *seq2_out_ret = seq2_out;
  return 0;
}

int do_trace_back ( unsigned char *bit_trace, char *seq1, char *seq2,
		   int seq1_len, int seq2_len, char **seq1_out_ret, char **seq2_out_ret,
		   int *seq_out_len, int b_r, int b_c, int b_e,
		   int band, int first_band_left, int first_row, 
		   int band_length, char PAD_SYM ) {

  int i, j, r, c, e;
  unsigned char trace_byte;
  int byte, band_left, max_seq;
  char *seq1_res, *seq2_res;
  char *seq1_out, *seq2_out;

  max_seq = seq1_len + seq2_len + 1;

  if(!(seq1_out = (char *) xmalloc(sizeof(char) * max_seq))) {
    verror(ERR_WARN, "affine_align", "malloc failed in do_trace_back");
    return -1;
  }
  if(!(seq2_out = (char *) xmalloc(sizeof(char) * max_seq))) {
    if ( seq1_out ) xfree (seq1_out);
    verror(ERR_WARN, "affine_align", "malloc failed in do_trace_back");
    return -1;
  }

  seq1_res = seq1_out;
  seq2_res = seq2_out;
		    
  for(i=0;i<max_seq-1;i++) {
    *(seq1_out++) = PAD_SYM;
    *(seq2_out++) = PAD_SYM;
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

    byte = e;
    trace_byte = bit_trace[byte];
      
    if(trace_byte == BYTE_DIAGONAL) {
      r--,c--;
      *(seq1_out--) = seq1[c];
      *(seq2_out--) = seq2[r];
    }
    else if(trace_byte == BYTE_DOWN) {
      r--;
      /*
       * Only gap the consensus vector if the single sequence is not a pad
       * itself. Otherwise we consider this to be a non-edit. Why? Because
       * we left seq2 padded and an insertion to seq1 is effectively removing
       * the pad from seq2. Not depadding seq2 means that realigning data
       * previously aligned and padded can be done with a narrow band around
       * the diagonal.
       */
      if (seq2[r] != '*') {
	  seq1_out--;
	  *(seq2_out--) = seq2[r];
      }
    }
    else {
      c--;
      seq2_out--;
      *(seq1_out--) = seq1[c];
    }
    if(band) {
      band_left = first_band_left + r - first_row;
      e = ((r - first_row + 1) * band_length) + (c - band_left + 1);
    }
    else {
      e = r * (seq1_len + 1) + c;
    }
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

  c = MAX(strlen ( seq1_res ),strlen ( seq2_res ));

  for ( i=0;i<c;i++ ) {
    if ( (seq1_res[i] != PAD_SYM) || (seq2_res[i] != PAD_SYM) ) break;
  }
  r = i;
  for ( i=0,j=r;j<c;i++,j++ ) {
    seq1_res[i] = seq1_res[j];
    seq2_res[i] = seq2_res[j];
  }
  seq1_res[i] = seq2_res[i] = '\0';
  *seq_out_len = i;
  seq1_out = seq1_res;
  seq2_out = seq2_res;
  *seq1_out_ret = seq1_out;
  *seq2_out_ret = seq2_out;
  return 0;
}

void destroy_af_mem ( int *F1, int *F2, int *G1, int *G2, int *H1, int *H2,
		      unsigned char *bit_trace, char *seq1_out, char *seq2_out ) {

    if ( F1 ) xfree ( F1);
    if ( G1 ) xfree ( G1);
    if ( H1 ) xfree ( H1);
    if ( F2 ) xfree ( F2);
    if ( G2 ) xfree ( G2);
    if ( H2 ) xfree ( H2);
    if ( bit_trace ) xfree ( bit_trace);
    if ( seq1_out ) xfree ( seq1_out);
    if ( seq2_out ) xfree ( seq2_out);
}



int affine_align3(OVERLAP *overlap, ALIGN_PARAMS *params) {
  /* the one using 3 tables */
  char *seq1, *seq2;
  int seq1_len, seq2_len, seq_out_len;
  int gap_open, gap_extend, edge_inc;
  int i,j;
  int s,*score_matrix_p;
  char *seq1_out, *seq2_out;
  int b_c, b_r;

  int t,big_neg,b_s,e,b_e;
  int *F1, *F2, *G1, *G2, *H1, *H2;
  int *pF1, *pF2, *pG1, *pG2, *pH1, *pH2;
  int *t_pF2, *t_pG2, *t_pH2;
  int best_F1, best_G1, best_H1, V_diag, V_extx, V_exty, V_insx, V_insy;
  int E_gap, F_gap;
  int edge_mode, best_edge_score;

  int band, band_length, two_band, band_left, band_right, first_band_left=0;
  int off_set, guard_offset, *pF_guard, *pG_guard, *pH_guard;
  int row, first_row, max_row, column, max_col;
  unsigned char *bit_trace;
  int byte, nibble, e_row, e_col, size_mat;
  char OLD_PAD_SYM, NEW_PAD_SYM;
  int gap_open_x, gap_open_y, gap_extend_x, gap_extend_y, gap_to_gap;

  /*
   *    Three possible alignment cases:
   *    IGAxi   AIGAHxi   GAxi--
   *    LGVyj   GVyj--    SLGVHyj
   *       F      G            H
   *    i.e. xi aligned with yj, xi aligned opposite a gap y,
   *    or yi aligned opposite a gap in x
   *    below these cases are contained in the recurrence relations
   *    for F, G and H respectively
   *    s(xi,yj) is score matrix
   *    d is gap_open
   *    e is gap extend
   *
   *                   F(i-1,j-1) + s(xi,yi)
   *    F(i,j)  = max  H(i-1,j-1) + s(xi,yi)      \  no gap
   *                   G(i-1,j-1) + s(xi,yi)
   *
   *                   F(i,j-1)   - d
   *    G(i,j) = max   G(i,j-1)   - e             |  gap in y
   *
   *                   F(i-1,j)   - d
   *    H(i,j) = max   H(i-1,j)   - e             -  gap in x
   *                  
   *              
   *    Find MAX(F(i,j),G(i,j),H(i,j)) and set trace accordingly:
   *                \     |      -
   *
   *    if gaps at left edge count:
   *    G(0,i) = G(i,0) = H(0,i) = H(i,0) = - d - e * i
   *    F(1,i) = F(i,1) = - d - e * i;
   *    F(0,0) = 0;
   *    otherwise all set to 0;
   *    if right end gaps count the best score is at (seq1_len,seq2_len)
   *    otherwise find the best score along the two edges
   *    trace back accordingly
   *
   *    store 2 rows for each of F, G, H
   *    use p_F1, p_G1, p_G1 to point to previous row
   *    p_F2, p_G2, p_H2 for current row being built
   *    at the start of a new row:
   *
   *    rows have length seq1_len, columns seq2_len
   *    i.e.
   *    rows: 1 - seq1_len, columns 1 - seq2_len
   *    seq1xxxxxxxxxxxxxxx
   *   s
   *   e
   *   q
   *   2
   *   y
   *   y
   *   y
   *
   *
   */

  F1 = F2 = G1 = G2 = H1 = H2 = NULL;
  bit_trace = NULL;
  seq1_out = seq2_out = NULL;
  big_neg = INT_MIN/2;
  best_edge_score = big_neg;

  seq1 = overlap->seq1;
  seq2 = overlap->seq2;
  seq1_len = overlap->seq1_len;
  seq2_len = overlap->seq2_len;

  edge_mode = params->edge_mode;
  gap_open = params->gap_open;
  gap_extend = params->gap_extend;
  OLD_PAD_SYM = params->old_pad_sym;
  NEW_PAD_SYM = params->new_pad_sym;
  band = params->band;
  first_row = params->first_row;
  band_left = params->band_left;
  band_right = params->band_right;
  edge_inc = gap_extend;
  gap_to_gap = -W128[(int)OLD_PAD_SYM][(int)OLD_PAD_SYM];
  gap_to_gap = 1;

  /* init tables */

  if(!(F1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for F1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(F2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for F2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(G1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for G1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(G2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for G2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(H1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for H1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(H2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for H2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }

  /* do recurrence */

  if ( edge_mode & EDGE_GAPS_COUNT ) {
    F1[0] = 0;
    for(i = 1,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) F1[i] = j;
    for(i = 0,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) G1[i] = j;
    for(i = 0,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) H1[i] = j;
    E_gap = -gap_open - edge_inc;
    F_gap = -gap_open;
  }
  else if ( edge_mode & EDGE_GAPS_ZERO ) {
    for(i = 0; i <= seq1_len; i++) F1[i] = 0;
    for(i = 0; i <= seq1_len; i++) G1[i] = 0;
    for(i = 0; i <= seq1_len; i++) H1[i] = 0;
    edge_inc = 0;
    E_gap = 0;
    F_gap = 0;
  }
  else {
    printf("scream: unknown gaps mode\n");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  
  /* process each row. i.e. each character of seq2 */

  b_s = big_neg;
  b_e = b_r = b_c = 0;
  t = 0;

  if ( band ) {
	/*
	 * If band is set we restrict the search to band elements either side
	 * of the main diagonal. This gives a band length of (2 * band) + 1.
	 * For the affine alignment it is necessary to know what happened in
	 * the elements to the left and above the current one, and therefore
	 * need to add another 2 to the band length.
	 */

	band_length = (2 * band) + 3;
	/*
	 * Need to add one to first_row, first_column, band_left and
	 * band_right because the alignment matrix has an extra row
	 * and an extra column over the number given by the lengths
	 * of the two sequences.
	 */
	size_mat = (MIN(seq1_len - band_left, seq2_len - first_row) + 1) 
	            * band_length;

	if(!(bit_trace = (unsigned char *) 
	    xmalloc(1 + sizeof(char) * size_mat / 4))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, 
			     seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat / 4;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	first_row++;
	band_left++;
	band_right++;
	first_band_left = band_left;

	two_band = band * 2;
	max_row = MIN(seq2_len, first_row + seq1_len - band_left + 1);

	for(row = first_row; row <= max_row; row++, band_left++, band_right++) {

	  guard_offset = band_left + two_band;

	  if(t == 0) {
	    pF1   = F1;
	    pF2   = F2;
	    pG1   = G1;
	    pG2   = G2;
	    pH1   = H1;
	    pH2   = H2;
	    pF_guard = F1 + guard_offset;
	    pG_guard = G1 + guard_offset;
	    pH_guard = H1 + guard_offset;
	    F2[0] = F_gap;
	    H2[0] = E_gap;
	    G2[0] = E_gap;
	    t = 1;
	  } else {
	    pF1   = F2;
	    pF2   = F1;
	    pG1   = G2;
	    pG2   = G1;
	    pH1   = H2;
	    pH2   = H1;
	    pF_guard = F2 + guard_offset;
	    pG_guard = G2 + guard_offset;
	    pH_guard = H2 + guard_offset;
	    F1[0] = F_gap;
	    H1[0] = E_gap;
	    G1[0] = E_gap;
	    t = 0;
	  }
	  if ( (off_set = band_left - 1 ) > 0 ) {
	    pF1 += off_set;
	    pF2 += off_set;
	    pG1 += off_set;
	    pG2 += off_set;
	    pH1 += off_set;
	    pH2 += off_set;
	    *pF2 = big_neg;
	    *pG2 = big_neg;
	    *pH2 = big_neg;
	  }
	  t_pF2 = pF2;
	  t_pG2 = pG2;
	  t_pH2 = pH2;

	  if ( band_right <= seq1_len ) {
	    *pF_guard = big_neg;
	    *pG_guard = big_neg;
	    *pH_guard = big_neg;
	  }
	  E_gap -= edge_inc;
	  F_gap -= edge_inc;

	  score_matrix_p = W128[(int)seq2[row-1]];
	  if ( seq2[row-1] != OLD_PAD_SYM ) {
	    gap_open_y = gap_open;
	    gap_extend_y = gap_extend;
	  }
	  else {
	    gap_open_y = gap_extend_y = gap_to_gap;
	  }
	    /*
	     * Need to prevent the possibility of certain transitions
	     * at the ends of each band; column-1 should not be allowed at
	     * the left-hand end of a band, and row-1 should not be
	     * allowed at the right-hand end. Such transitions would
	     * take the trace-back out of the the parts of the 
	     * alignment matrix that are covered by the band.
	     *
	     *      0 1 2 3 4 5 6
	     *    0   _ _ _           For example, shouldn't consider
	     *    1  |_|_|_|_         column-1 at position (row = 3, col = 3)
	     *    2    |_|_|_|_       so give (3, 2) a very low score.
	     *    3      |_|_|_|_     Likewise, shouldn't consider row-1
	     *    4        |_|_|_|    at position (3, 5), so give (2, 5)
	     *    5          |_|_|    a very low score.
	     * this is done using the guard pointers at the right edge
	     * and initialising the left edges where required.   
	     */

	  /* process each column. i.e. each character of seq1 */
	  
	  max_col = MIN(seq1_len, band_right);
	  for(column = MAX(1, band_left);
		column <= max_col;
		column++, pF1++, pF2++, pG1++, pG2++, pH1++, pH2++) {

	    /* move diagonally? */
	    s = score_matrix_p[(int)seq1[column-1]];
	    if ( seq1[column-1] != OLD_PAD_SYM ) {
	      gap_open_x = gap_open;
	      gap_extend_x = gap_extend;
	    }
	    else {
	      gap_open_x = gap_extend_x = gap_to_gap;
	    }

	    V_diag = *pF1 + s;
	    V_insx = *pH1 + s;
	    V_insy = *pG1 + s;
	    if ( V_insx > V_diag ) {
	      if ( V_insx > V_insy ) {
		best_F1 = V_insx;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    else {
	      if ( V_diag > V_insy ) {
		best_F1 = V_diag;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    *(pF2+1) = best_F1;

	    /* gap in x? */
	    V_diag =  *pF2 - gap_open_x;
	    V_extx =  *pH2 - gap_extend_x;
	    if ( V_diag > V_extx ) {
	      best_H1 = V_diag;
	    }
	    else {
	      best_H1 = V_extx;
	    }
	    *(pH2+1) = best_H1;

	    /* gap in y? */
	    V_diag = *(pF1+1) - gap_open_y;
	    V_exty = *(pG1+1) - gap_extend_y;
	    if ( V_diag > V_exty ) {
	      best_G1 = V_diag;
	    }
	    else {
	      best_G1 = V_exty;
	    }
	    *(pG2+1) = best_G1;

	    e_row = (row - first_row + 1) * band_length;
	    e_col = column - band_left + 1;
	    e = e_row + e_col;
	    byte = e / 4;
	    nibble = 2 * (e % 4);
	    
	    /* find the best move */
	    if ( best_H1 > best_F1 ) {
	      if ( best_H1 > best_G1 ) {
		bit_trace[byte] |= BYTE_ACROSS << nibble;
		b_s = best_H1;
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
	      }
	    }
	    else {
	      if ( best_F1 > best_G1 ) {
		bit_trace[byte] |= BYTE_DIAGONAL << nibble;
		b_s = best_F1;
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
	      }
	    }
	  }
	  
	  if ( column > seq1_len ) {
	    if ( edge_mode & BEST_EDGE_TRACE ) {
	      best_H1 = MAX(best_H1,best_G1);
	      best_F1 = MAX(best_H1,best_F1);
	      if ( best_F1 > best_edge_score ) {
		best_edge_score = best_F1;
		b_r = row;
		b_e = ((row - first_row + 1) * band_length) + 
		       (seq1_len - band_left + 1);
	      }
	    }
	  }
      }
	
	
        row--;
	band_left--;
	band_right--;
	if ( edge_mode & BEST_EDGE_TRACE ) {
	  b_c = seq1_len;

	  pF2 = t_pF2+1;
	  pG2 = t_pG2+1;
	  pH2 = t_pH2+1;
	  max_col = MIN(seq1_len, band_right);
	  for(column = MAX(1, band_left);
	      column <= max_col;
	      column++, pF2++, pG2++, pH2++) {
	    best_F1 = *pF2;
	    best_G1 = *pG2;
	    best_H1 = *pH2;
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_c = column;
	      b_r = seq2_len;
	      b_e = ((row - first_row + 1) * band_length) + 
		     (column - band_left + 1);
	    }
	  }
      }
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	    b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = ((b_r - first_row + 1) * band_length) + 
		 (b_c - band_left + 1);
	  best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = ((b_r - first_row + 1) * band_length) + 
		 (b_c - band_left + 1);
	  best_edge_score = b_s;
        }

  }
  else {

	/* Initialise the bit trace */
	size_mat = (seq1_len + 1) * (seq2_len + 1);
	if(!(bit_trace = (unsigned char *) xmalloc(1 + sizeof(char) * size_mat / 4))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat / 4;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	for(row = 1, e = seq1_len + 2; row <= seq2_len; row++, e++) {

	  if(t == 0) {
	    pF1   = F1;
	    pF2   = F2;
	    pG1   = G1;
	    pG2   = G2;
	    pH1   = H1;
	    pH2   = H2;
	    F2[0] = F_gap;
	    H2[0] = E_gap;
	    G2[0] = E_gap;
	    t = 1;
	  } else {
	    pF1   = F2;
	    pF2   = F1;
	    pG1   = G2;
	    pG2   = G1;
	    pH1   = H2;
	    pH2   = H1;
	    F1[0] = F_gap;
	    H1[0] = E_gap;
	    G1[0] = E_gap;
	    t = 0;
	  }
	  
	  E_gap -= edge_inc;
	  F_gap -= edge_inc;

	  score_matrix_p = W128[(int)seq2[row-1]];
	  if ( seq2[row-1] != OLD_PAD_SYM ) {
	    gap_open_y = gap_open;
	    gap_extend_y = gap_extend;
	  }
	  else {
	    gap_open_y = gap_extend_y = gap_to_gap;
	  }
	  
	  /* process each column. i.e. each character of seq1 */
	  
	  for(column = 1; column <= seq1_len; column++, e++) {
	    /* move diagonally? */
	    s = score_matrix_p[(int)seq1[column-1]];
	    if ( seq1[column-1] != OLD_PAD_SYM ) {
	      gap_open_x = gap_open;
	      gap_extend_x = gap_extend;
	    }
	    else {
	      gap_open_x = gap_extend_x = gap_to_gap;
	    }
	    V_diag = pF1[column-1] + s;
	    V_insx = pH1[column-1] + s;
	    V_insy = pG1[column-1] + s;
	    if ( V_insx > V_diag ) {
	      if ( V_insx > V_insy ) {
		best_F1 = V_insx;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    else {
	      if ( V_diag > V_insy ) {
		best_F1 = V_diag;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    pF2[column] = best_F1;

	    printf("%3d %3d %3d %3d %3d ",row,column,V_diag,V_insx,V_insy);
	    
	    /* gap in x? */
	    V_diag =  pF2[column-1] - gap_open_x;
	    V_extx =  pH2[column-1] - gap_extend_x;
	    if ( V_diag > V_extx ) {
	      best_H1 = V_diag;
	    }
	    else {
	      best_H1 = V_extx;
	    }
	    pH2[column] = best_H1;

	    printf("%3d %3d ",V_diag,V_extx);
	    
	    /* gap in y? */
	    V_diag =  pF1[column] - gap_open_y;
	    V_exty = pG1[column] - gap_extend_y;
	    if ( V_diag > V_exty ) {
	      best_G1 = V_diag;
	    }
	    else {
	      best_G1 = V_exty;
	    }
	    pG2[column] = best_G1;

	    printf("%3d %3d %3d %3d %3d ",V_diag,V_exty,s,gap_open_x,gap_open_y);

	    byte = e / 4;
	    nibble = 2 * (e % 4);

	    /* choose best move */
	    if ( best_H1 > best_F1 ) {
	      if ( best_H1 > best_G1 ) {
		bit_trace[byte] |= BYTE_ACROSS << nibble;
		b_s = best_H1;
		printf(" -");
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
		printf(" |");
	      }
	    }
	    else {
	      if ( best_F1 > best_G1 ) {
		bit_trace[byte] |= BYTE_DIAGONAL << nibble;
		b_s = best_F1;
		printf(" \\");
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
		printf(" |");
	      }
	    }
	    printf("\n");

	  }
	  
	  if ( edge_mode & BEST_EDGE_TRACE ) {
	    /* best right edge score */
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_r = row;
	      b_e = (row + 1) * (seq1_len + 1) - 1;
	    }
	  }
	}
	
	if ( edge_mode & BEST_EDGE_TRACE ) {
	  /* best bottom edge score */
	  b_c = seq1_len;
	  for(column = 1; column <= seq1_len; column++) {
	    best_F1 = pF2[column];
	    best_G1 = pG2[column];
	    best_H1 = pH2[column];
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_c = column;
	      b_r = seq2_len;
	      b_e = (row - 1) * (seq1_len + 1) + column;
	    }
	  }
	}
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = seq2_len * (seq1_len + 1) + seq1_len;
	  best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = seq2_len * (seq1_len + 1) + seq1_len;
	  best_edge_score = b_s;
        }
    }


  /* do traceback */

  overlap->score = best_edge_score;

  if( i = do_trace_back_bits ( bit_trace, seq1, seq2, seq1_len, seq2_len,
		      &seq1_out, &seq2_out, &seq_out_len, b_r, b_c, b_e,
		      band, first_band_left, first_row, band_length, NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
      return -1;
  }

  overlap->seq1_out = seq1_out;
  overlap->seq2_out = seq2_out;
  overlap->seq_out_len = seq_out_len;

  if ( i = seq_to_overlap (overlap, OLD_PAD_SYM, NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
  }

  if ( params->job & RETURN_EDIT_BUFFERS ) {
    if (seq_to_edit ( seq1_out,seq_out_len,&overlap->S1,&overlap->s1_len,NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
    }
    if (seq_to_edit ( seq2_out,seq_out_len,&overlap->S2,&overlap->s2_len,NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
    }
  }

  if ( params->job & RETURN_SEQ ) {
    if ( !(params->job & RETURN_NEW_PADS) ) {
      old_pads_for_new(seq1_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
      old_pads_for_new(seq2_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
    }
    seq1_out = seq2_out = NULL; /* stop them being freed! */
  }
  else {
    overlap->seq1_out = overlap->seq2_out = NULL;
    /* ie we let destroy_af_mem free the memory, but we must
     * ensure that othr routines do not try to free it too 
     */
  }
  destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );

  return 0;
}

int affine_align_bits(OVERLAP *overlap, ALIGN_PARAMS *params) {
  /* the one using 3 tables */
  char *seq1, *seq2;
  int seq1_len, seq2_len, seq_out_len;
  int gap_open, gap_extend, edge_inc;
  int i,j;
  int s,*score_matrix_p;
  char *seq1_out, *seq2_out;
  int b_c, b_r;

  int t,big_neg,b_s,e,b_e;
  int *F1, *F2, *G1, *G2, *H1, *H2;
  int *pF1, *pF2, *pG1, *pG2, *pH1, *pH2;
  int *t_pF2, *t_pG2, *t_pH2;
  int best_F1, best_G1, best_H1, V_diag, V_extx, V_exty, V_insx, V_insy;
  int F_gap, start_edge_pens_x, start_edge_pens_y, G_gap, H_gap, g;
  int edge_mode, best_edge_score;

  int band, band_length, two_band, band_left, band_right, first_band_left=0;
  int off_set, guard_offset, *pF_guard, *pG_guard, *pH_guard;
  int row, first_row, max_row, last_row, column, max_col, last_column;
  unsigned char *bit_trace;
  int byte, nibble, e_row, e_col, size_mat;
  char OLD_PAD_SYM, NEW_PAD_SYM;
  int gap_open_x, gap_open_y, gap_extend_x, gap_extend_y, gap_to_gap;

  /*
   *    Three possible alignment cases:
   *    IGAxi   AIGAHxi   GAxi--
   *    LGVyj   GVyj--    SLGVHyj
   *       F      G            H
   *    i.e. xi aligned with yj, xi aligned opposite a gap y,
   *    or yi aligned opposite a gap in x
   *    below these cases are contained in the recurrence relations
   *    for F, G and H respectively
   *    s(xi,yj) is score matrix
   *    d is gap_open
   *    e is gap extend
   *
   *                   F(i-1,j-1) + s(xi,yi)
   *    F(i,j)  = max  H(i-1,j-1) + s(xi,yi)      \  no gap
   *                   G(i-1,j-1) + s(xi,yi)
   *
   *                   F(i,j-1)   - d
   *    G(i,j) = max   G(i,j-1)   - e             |  gap in y
   *
   *                   F(i-1,j)   - d
   *    H(i,j) = max   H(i-1,j)   - e             -  gap in x
   *                  
   *              
   *    Find MAX(F(i,j),G(i,j),H(i,j)) and set trace accordingly:
   *                \     |      -
   *
   *    if gaps at left edge count:
   *    G(0,i) = G(i,0) = H(0,i) = H(i,0) = - d - e * i
   *    F(1,i) = F(i,1) = - d - e * i;
   *    F(0,0) = 0;
   *    otherwise all set to 0;
   *    if right end gaps count the best score is at (seq1_len,seq2_len)
   *    otherwise find the best score along the two edges
   *    trace back accordingly
   *
   *    store 2 rows for each of F, G, H
   *    use p_F1, p_G1, p_G1 to point to previous row
   *    p_F2, p_G2, p_H2 for current row being built
   *    at the start of a new row:
   *
   *    rows have length seq1_len, columns seq2_len
   *    i.e.
   *    rows: 1 - seq1_len, columns 1 - seq2_len
   *    seq1xxxxxxxxxxxxxxx
   *   s
   *   e
   *   q
   *   2
   *   y
   *   y
   *   y
   *
   *
   */

  F1 = F2 = G1 = G2 = H1 = H2 = NULL;
  bit_trace = NULL;
  seq1_out = seq2_out = NULL;
  big_neg = INT_MIN/2;
  best_edge_score = big_neg;

  seq1 = overlap->seq1;
  seq2 = overlap->seq2;
  seq1_len = overlap->seq1_len;
  seq2_len = overlap->seq2_len;

  edge_mode = params->edge_mode;
  gap_open = params->gap_open;
  gap_extend = params->gap_extend;
  OLD_PAD_SYM = params->old_pad_sym;
  NEW_PAD_SYM = params->new_pad_sym;
  band = params->band;
  first_row = params->first_row;
  band_left = params->band_left;
  band_right = params->band_right;
  edge_inc = gap_extend;
  gap_to_gap = -W128[(int)OLD_PAD_SYM][(int)OLD_PAD_SYM];
  gap_to_gap = 1;

  /* init tables */

  if(!(F1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for F1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(F2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for F2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(G1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for G1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(G2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for G2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(H1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for H1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(H2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for H2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }

  /* do recurrence */


  if ( edge_mode & EDGE_GAPS_COUNT ) {
    F1[0] = 0;

    /* deal with pads at start of seq1: if present we set the edge
     * scores for G */

    for(i = 0; i < seq1_len; i++) {
	if (seq1[i] != OLD_PAD_SYM) break;
    }
    start_edge_pens_y = i;

    for(j=0,g=-gap_to_gap;j<start_edge_pens_y;j++,g-=gap_to_gap) G1[j] = g;

    g = start_edge_pens_y ? g+gap_to_gap : -gap_open;
    for(; i <= seq1_len; i++,g-=edge_inc) G1[i] = g;

    /* deal with pads at start of seq2: if present we set the edge
     * scores for H */

    for(i = 0; i < seq2_len; i++) {
	if (seq2[i] != OLD_PAD_SYM) break;
    }
    start_edge_pens_x = i;

    g = start_edge_pens_x ? -gap_to_gap : -gap_open;
    H_gap = g;

    for(i = 0,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) H1[i] = j;

    for(i = 1,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) F1[i] = j;

    /*
    for(i = 1,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) F1[i] = j;
    for(i = 0,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) G1[i] = j;
    for(i = 0,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) H1[i] = j;
    */
    G_gap = -gap_open - edge_inc;
    F_gap = -gap_open;
  }
  else if ( edge_mode & EDGE_GAPS_ZERO ) {
    for(i = 0; i <= seq1_len; i++) F1[i] = 0;
    for(i = 0; i <= seq1_len; i++) G1[i] = -gap_open;
    for(i = 0; i <= seq1_len; i++) H1[i] = -gap_open;
    edge_inc = 0;
    F_gap = 0;
    G_gap = -gap_open;
    H_gap = -gap_open;
    start_edge_pens_x = -1;
  }
  else {
    printf("scream: unknown gaps mode\n");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  
  /* process each row. i.e. each character of seq2 */

  b_s = big_neg;
  b_e = b_r = b_c = 0;
  t = 0;

  if ( band ) {
	/*
	 * If band is set we restrict the search to band elements either side
	 * of the main diagonal. This gives a band length of (2 * band) + 1.
	 * For the affine alignment it is necessary to know what happened in
	 * the elements to the left and above the current one, and therefore
	 * need to add another 2 to the band length.
	 */

	band_length = (2 * band) + 3;
	/*
	 * Need to add one to first_row, first_column, band_left and
	 * band_right because the alignment matrix has an extra row
	 * and an extra column over the number given by the lengths
	 * of the two sequences.
	 */
	size_mat = (MIN(seq1_len - band_left, seq2_len - first_row) + 1) 
	            * band_length;

	if(!(bit_trace = (unsigned char *) 
	    xmalloc(1 + sizeof(char) * size_mat / 4))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, 
			     seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat / 4;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	first_row++;
	band_left++;
	band_right++;
	first_band_left = band_left;

	two_band = band * 2;
	max_row = MIN(seq2_len, first_row + seq1_len - band_left + 1);
	last_row = 0;

	for(row = first_row; row <= max_row; row++, band_left++, band_right++) {

	  guard_offset = band_left + two_band;

	  if(t == 0) {
	    pF1   = F1;
	    pF2   = F2;
	    pG1   = G1;
	    pG2   = G2;
	    pH1   = H1;
	    pH2   = H2;
	    pF_guard = F1 + guard_offset;
	    pG_guard = G1 + guard_offset;
	    pH_guard = H1 + guard_offset;
	    F2[0] = F_gap;
	    G2[0] = G_gap;
	    H2[0] = H_gap;
	    F_gap -= edge_inc;
	    G_gap -= edge_inc;
	    if ( row > start_edge_pens_x ) {
		H_gap -= edge_inc;
	    }
	    else {
	      H_gap -= gap_to_gap;
	    }
	    t = 1;
	  } else {
	    pF1   = F2;
	    pF2   = F1;
	    pG1   = G2;
	    pG2   = G1;
	    pH1   = H2;
	    pH2   = H1;
	    pF_guard = F2 + guard_offset;
	    pG_guard = G2 + guard_offset;
	    pH_guard = H2 + guard_offset;
	    F1[0] = F_gap;
	    G1[0] = G_gap;
	    H1[0] = H_gap;
	    F_gap -= edge_inc;
	    G_gap -= edge_inc;
	    if ( row > start_edge_pens_x ) {
		H_gap -= edge_inc;
	    }
	    else {
	      H_gap -= gap_to_gap;
	    }
	    t = 0;
	  }
	  if ( (off_set = band_left - 1 ) > 0 ) {
	    pF1 += off_set;
	    pF2 += off_set;
	    pG1 += off_set;
	    pG2 += off_set;
	    pH1 += off_set;
	    pH2 += off_set;
	    *pF2 = big_neg;
	    *pG2 = big_neg;
	    *pH2 = big_neg;
	  }
	  t_pF2 = pF2;
	  t_pG2 = pG2;
	  t_pH2 = pH2;

	  if ( band_right <= seq1_len ) {
	    *pF_guard = big_neg;
	    *pG_guard = big_neg;
	    *pH_guard = big_neg;
	  }

	  score_matrix_p = W128[(int)seq2[row-1]];
	  if ( seq2[row-1] != OLD_PAD_SYM ) {
	    gap_open_y = gap_open;
	    gap_extend_y = gap_extend;
	  }
	  else {
	    gap_open_y = gap_extend_y = gap_to_gap;
	  }
	    /*
	     * Need to prevent the possibility of certain transitions
	     * at the ends of each band; column-1 should not be allowed at
	     * the left-hand end of a band, and row-1 should not be
	     * allowed at the right-hand end. Such transitions would
	     * take the trace-back out of the the parts of the 
	     * alignment matrix that are covered by the band.
	     *
	     *      0 1 2 3 4 5 6
	     *    0   _ _ _           For example, shouldn't consider
	     *    1  |_|_|_|_         column-1 at position (row = 3, col = 3)
	     *    2    |_|_|_|_       so give (3, 2) a very low score.
	     *    3      |_|_|_|_     Likewise, shouldn't consider row-1
	     *    4        |_|_|_|    at position (3, 5), so give (2, 5)
	     *    5          |_|_|    a very low score.
	     * this is done using the guard pointers at the right edge
	     * and initialising the left edges where required.   
	     */

	  /* process each column. i.e. each character of seq1 */
	  
	  max_col = MIN(seq1_len, band_right);
	  for(column = MAX(1, band_left);
		column <= max_col;
		column++, pF1++, pF2++, pG1++, pG2++, pH1++, pH2++) {

	      last_row = MAX(row,last_row);

	    /* move diagonally? */
	    s = score_matrix_p[(int)seq1[column-1]];
	    if ( seq1[column-1] != OLD_PAD_SYM ) {
	      gap_open_x = gap_open;
	      gap_extend_x = gap_extend;
	    }
	    else {
	      gap_open_x = gap_extend_x = gap_to_gap;
	    }

	    V_diag = *pF1 + s;
	    V_insx = *pH1 + s;
	    V_insy = *pG1 + s;
	    if ( V_insx > V_diag ) {
	      if ( V_insx > V_insy ) {
		best_F1 = V_insx;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    else {
	      if ( V_diag > V_insy ) {
		best_F1 = V_diag;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    *(pF2+1) = best_F1;

	    /* gap in x? */
	    V_diag =  *pF2 - gap_open_x;
	    V_extx =  *pH2 - gap_extend_x;
	    if ( V_diag > V_extx ) {
	      best_H1 = V_diag;
	    }
	    else {
	      best_H1 = V_extx;
	    }
	    *(pH2+1) = best_H1;

	    /* gap in y? */
	    V_diag = *(pF1+1) - gap_open_y;
	    V_exty = *(pG1+1) - gap_extend_y;
	    if ( V_diag > V_exty ) {
	      best_G1 = V_diag;
	    }
	    else {
	      best_G1 = V_exty;
	    }
	    *(pG2+1) = best_G1;

	    e_row = (row - first_row + 1) * band_length;
	    e_col = column - band_left + 1;
	    e = e_row + e_col;
	    byte = e / 4;
	    nibble = 2 * (e % 4);
	    
	    /* find the best move */
	    if ( best_H1 > best_F1 ) {
	      if ( best_H1 > best_G1 ) {
		bit_trace[byte] |= BYTE_ACROSS << nibble;
		b_s = best_H1;
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
	      }
	    }
	    else {
	      if ( best_F1 > best_G1 ) {
		bit_trace[byte] |= BYTE_DIAGONAL << nibble;
		b_s = best_F1;
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
	      }
	    }
	  }
	  
	  if ( column > seq1_len ) {
	    if ( edge_mode & BEST_EDGE_TRACE ) {
	      best_H1 = MAX(best_H1,best_G1);
	      best_F1 = MAX(best_H1,best_F1);
	      if ( best_F1 > best_edge_score ) {
		best_edge_score = best_F1;
		b_r = row;
		b_e = ((row - first_row + 1) * band_length) + 
		       (seq1_len - band_left + 1);
	      }
	    }
	  }
      }
	
	
	last_column = max_col;
	row = max_row;
	band_left--;
	band_right--;
	if ( edge_mode & BEST_EDGE_TRACE ) {
	  b_c = last_column;

	  pF2 = t_pF2+1;
	  pG2 = t_pG2+1;
	  pH2 = t_pH2+1;
	  max_col = MIN(seq1_len, band_right);
	  for(column = MAX(1, band_left);
	      column <= max_col;
	      column++, pF2++, pG2++, pH2++) {
	    best_F1 = *pF2;
	    best_G1 = *pG2;
	    best_H1 = *pH2;
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_c = column;
	      b_r = last_row;
	      b_e = ((row - first_row + 1) * band_length) + 
		     (column - band_left + 1);
	    }
	  }
      }
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	    b_r = last_row;
	  b_c = last_column;
	  b_e = ((b_r - first_row + 1) * band_length) + 
		 (b_c - band_left + 1);
	  best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	  b_r = last_row;
	  b_c = last_column;
	  b_e = ((b_r - first_row + 1) * band_length) + 
		 (b_c - band_left + 1);
	  best_edge_score = b_s;
        }

  }
  else {

	/* Initialise the bit trace */
	size_mat = (seq1_len + 1) * (seq2_len + 1);
	if(!(bit_trace = (unsigned char *) xmalloc(1 + sizeof(char) * size_mat / 4))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat / 4;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	for(row = 1, e = seq1_len + 2; row <= seq2_len; row++, e++) {

	  if(t == 0) {
	    pF1   = F1;
	    pF2   = F2;
	    pG1   = G1;
	    pG2   = G2;
	    pH1   = H1;
	    pH2   = H2;
	    F2[0] = F_gap;
	    G2[0] = G_gap;
	    H2[0] = H_gap;
	    F_gap -= edge_inc;
	    G_gap -= edge_inc;
	    if ( row > start_edge_pens_x ) {
		H_gap -= edge_inc;
	    }
	    else {
	      H_gap -= gap_to_gap;
	    }
	    t = 1;
	  } else {
	    pF1   = F2;
	    pF2   = F1;
	    pG1   = G2;
	    pG2   = G1;
	    pH1   = H2;
	    pH2   = H1;
	    F1[0] = F_gap;
	    G1[0] = G_gap;
	    H1[0] = H_gap;
	    F_gap -= edge_inc;
	    G_gap -= edge_inc;
	    if ( row > start_edge_pens_x ) {
		H_gap -= edge_inc;
	    }
	    else {
	      H_gap -= gap_to_gap;
	    }
	    t = 0;
	  }

	  score_matrix_p = W128[(int)seq2[row-1]];
	  if ( seq2[row-1] != OLD_PAD_SYM ) {
	    gap_open_y = gap_open;
	    gap_extend_y = gap_extend;
	  }
	  else {
	    gap_open_y = gap_extend_y = gap_to_gap;
	  }
	  
	  /* process each column. i.e. each character of seq1 */
	  
	  for(column = 1; column <= seq1_len; column++, e++) {
	    /* move diagonally? */
	    s = score_matrix_p[(int)seq1[column-1]];
	    if ( seq1[column-1] != OLD_PAD_SYM ) {
	      gap_open_x = gap_open;
	      gap_extend_x = gap_extend;
	    }
	    else {
	      gap_open_x = gap_extend_x = gap_to_gap;
	    }
	    V_diag = pF1[column-1] + s;
	    V_insx = pH1[column-1] + s;
	    V_insy = pG1[column-1] + s;
	    if ( V_insx > V_diag ) {
	      if ( V_insx > V_insy ) {
		best_F1 = V_insx;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    else {
	      if ( V_diag > V_insy ) {
		best_F1 = V_diag;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    pF2[column] = best_F1;

	    /*printf("%3d %3d %3d %3d %3d ",row,column,V_diag,V_insx,V_insy);*/
	    
	    /* gap in x? */
	    V_diag =  pF2[column-1] - gap_open_x;
	    V_extx =  pH2[column-1] - gap_extend_x;
	    if ( V_diag > V_extx ) {
	      best_H1 = V_diag;
	    }
	    else {
	      best_H1 = V_extx;
	    }
	    pH2[column] = best_H1;

	    /*printf("%3d %3d ",V_diag,V_extx);*/
	    
	    /* gap in y? */
	    V_diag =  pF1[column] - gap_open_y;
	    V_exty = pG1[column] - gap_extend_y;
	    if ( V_diag > V_exty ) {
	      best_G1 = V_diag;
	    }
	    else {
	      best_G1 = V_exty;
	    }
	    pG2[column] = best_G1;

	    /*printf("%3d %3d %3d %3d %3d ",V_diag,V_exty,s,gap_open_x,gap_open_y);*/

	    byte = e / 4;
	    nibble = 2 * (e % 4);

	    /* choose best move */
	    if ( best_H1 > best_F1 ) {
	      if ( best_H1 > best_G1 ) {
		bit_trace[byte] |= BYTE_ACROSS << nibble;
		b_s = best_H1;
		/*printf(" -");*/
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
		/*printf(" |");*/
	      }
	    }
	    else {
	      if ( best_F1 > best_G1 ) {
		bit_trace[byte] |= BYTE_DIAGONAL << nibble;
		b_s = best_F1;
		/*printf(" \\");*/
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
		/*printf(" |");*/
	      }
	    }
	    /*printf("\n");*/

	  }
	  
	  if ( edge_mode & BEST_EDGE_TRACE ) {
	    /* best right edge score */
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_r = row;
	      b_e = (row + 1) * (seq1_len + 1) - 1;
	    }
	  }
	}
	
	if ( edge_mode & BEST_EDGE_TRACE ) {
	  /* best bottom edge score */
	  b_c = seq1_len;
	  for(column = 1; column <= seq1_len; column++) {
	    best_F1 = pF2[column];
	    best_G1 = pG2[column];
	    best_H1 = pH2[column];
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_c = column;
	      b_r = seq2_len;
	      b_e = (row - 1) * (seq1_len + 1) + column;
	    }
	  }
	}
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = seq2_len * (seq1_len + 1) + seq1_len;
	  best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = seq2_len * (seq1_len + 1) + seq1_len;
	  best_edge_score = b_s;
        }
    }


  /* do traceback */

  overlap->score = best_edge_score;

  if( i = do_trace_back_bits ( bit_trace, seq1, seq2, seq1_len, seq2_len,
		      &seq1_out, &seq2_out, &seq_out_len, b_r, b_c, b_e,
		      band, first_band_left, first_row, band_length, NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
      return -1;
  }

  overlap->seq1_out = seq1_out;
  overlap->seq2_out = seq2_out;
  overlap->seq_out_len = seq_out_len;

  if ( i = seq_to_overlap (overlap, OLD_PAD_SYM, NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
  }

  if ( params->job & RETURN_EDIT_BUFFERS ) {
    if (seq_to_edit ( seq1_out,seq_out_len,&overlap->S1,&overlap->s1_len,NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
    }
    if (seq_to_edit ( seq2_out,seq_out_len,&overlap->S2,&overlap->s2_len,NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
    }
  }

  if ( params->job & RETURN_SEQ ) {
    if ( !(params->job & RETURN_NEW_PADS) ) {
      old_pads_for_new(seq1_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
      old_pads_for_new(seq2_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
    }
    seq1_out = seq2_out = NULL; /* stop them being freed! */
  }
  else {
    overlap->seq1_out = overlap->seq2_out = NULL;
    /* ie we let destroy_af_mem free the memory, but we must
     * ensure that othr routines do not try to free it too 
     */
  }
  destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );

  return 0;
}


int affine_align_big(OVERLAP *overlap, ALIGN_PARAMS *params) {
  /* the one using 3 tables */
  char *seq1, *seq2;
  int seq1_len, seq2_len, seq_out_len;
  int gap_open, gap_extend, edge_inc;
  int i,j;
  int s,*score_matrix_p;
  char *seq1_out, *seq2_out;
  int b_c, b_r;

  int t,big_neg,b_s,e,b_e;
  int *F1, *F2, *G1, *G2, *H1, *H2;
  int *pF1, *pF2, *pG1, *pG2, *pH1, *pH2;
  int *t_pF2, *t_pG2, *t_pH2;
  int best_F1, best_G1, best_H1, V_diag, V_extx, V_exty, V_insx, V_insy;
  int F_gap, start_edge_pens_x, start_edge_pens_y, G_gap, H_gap, g;
  int edge_mode, best_edge_score;

  int band, band_length, two_band, band_left, band_right, first_band_left=0;
  int off_set, guard_offset, *pF_guard, *pG_guard, *pH_guard;
  int row, first_row, max_row, last_row, column, max_col, last_column;
  unsigned char *bit_trace;
  int e_row, e_col, size_mat;
  char OLD_PAD_SYM, NEW_PAD_SYM;
  int gap_open_x, gap_open_y, gap_extend_x, gap_extend_y, gap_to_gap;

  /*
   *    Three possible alignment cases:
   *    IGAxi   AIGAHxi   GAxi--
   *    LGVyj   GVyj--    SLGVHyj
   *       F      G            H
   *    i.e. xi aligned with yj, xi aligned opposite a gap y,
   *    or yi aligned opposite a gap in x
   *    below these cases are contained in the recurrence relations
   *    for F, G and H respectively
   *    s(xi,yj) is score matrix
   *    d is gap_open
   *    e is gap extend
   *
   *                   F(i-1,j-1) + s(xi,yi)
   *    F(i,j)  = max  H(i-1,j-1) + s(xi,yi)      \  no gap
   *                   G(i-1,j-1) + s(xi,yi)
   *
   *                   F(i,j-1)   - d
   *    G(i,j) = max   G(i,j-1)   - e             |  gap in y
   *
   *                   F(i-1,j)   - d
   *    H(i,j) = max   H(i-1,j)   - e             -  gap in x
   *                  
   *              
   *    Find MAX(F(i,j),G(i,j),H(i,j)) and set trace accordingly:
   *                \     |      -
   *
   *    if gaps at left edge count:
   *    G(0,i) = G(i,0) = H(0,i) = H(i,0) = - d - e * i
   *    F(1,i) = F(i,1) = - d - e * i;
   *    F(0,0) = 0;
   *    otherwise all set to 0;
   *    if right end gaps count the best score is at (seq1_len,seq2_len)
   *    otherwise find the best score along the two edges
   *    trace back accordingly
   *
   *    store 2 rows for each of F, G, H
   *    use p_F1, p_G1, p_G1 to point to previous row
   *    p_F2, p_G2, p_H2 for current row being built
   *    at the start of a new row:
   *
   *    rows have length seq1_len, columns seq2_len
   *    i.e.
   *    rows: 1 - seq1_len, columns 1 - seq2_len
   *    seq1xxxxxxxxxxxxxxx
   *   s
   *   e
   *   q
   *   2
   *   y
   *   y
   *   y
   *
   *
   */

  F1 = F2 = G1 = G2 = H1 = H2 = NULL;
  bit_trace = NULL;
  seq1_out = seq2_out = NULL;
  big_neg = INT_MIN/2;
  best_edge_score = big_neg;

  seq1 = overlap->seq1;
  seq2 = overlap->seq2;
  seq1_len = overlap->seq1_len;
  seq2_len = overlap->seq2_len;

  edge_mode = params->edge_mode;
  gap_open = params->gap_open;
  gap_extend = params->gap_extend;
  OLD_PAD_SYM = params->old_pad_sym;
  NEW_PAD_SYM = params->new_pad_sym;
  band = params->band;
  first_row = params->first_row;
  band_left = params->band_left;
  band_right = params->band_right;
  edge_inc = gap_extend;
  gap_to_gap = 1;

  /* init tables */

  if(!(F1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for F1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(F2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for F2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(G1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for G1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(G2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for G2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(H1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for H1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(H2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for H2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }

  /* do recurrence */

  if ( edge_mode & EDGE_GAPS_COUNT ) {
    F1[0] = 0;

    /* deal with pads at start of seq1: if present we set the edge
     * scores for G */

    for(i = 0; i < seq1_len; i++) {
	if (seq1[i] != OLD_PAD_SYM) break;
    }
    start_edge_pens_y = i;

    for(j=0,g=-gap_to_gap;j<start_edge_pens_y;j++,g-=gap_to_gap) G1[j] = g;

    g = start_edge_pens_y ? g+gap_to_gap : -gap_open;
    for(; i <= seq1_len; i++,g-=edge_inc) G1[i] = g;

    /* deal with pads at start of seq2: if present we set the edge
     * scores for H */

    for(i = 0; i < seq2_len; i++) {
	if (seq2[i] != OLD_PAD_SYM) break;
    }
    start_edge_pens_x = i;

    g = start_edge_pens_x ? -gap_to_gap : -gap_open;
    H_gap = g;

    for(i = 0,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) H1[i] = j;

    for(i = 1,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) F1[i] = j;

    /*
    for(i = 1,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) F1[i] = j;
    for(i = 0,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) G1[i] = j;
    for(i = 0,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) H1[i] = j;
    */
    G_gap = -gap_open - edge_inc;
    F_gap = -gap_open;
  }
  else if ( edge_mode & EDGE_GAPS_ZERO ) {
    for(i = 0; i <= seq1_len; i++) F1[i] = 0;
    for(i = 0; i <= seq1_len; i++) G1[i] = -gap_open;
    for(i = 0; i <= seq1_len; i++) H1[i] = -gap_open;
    edge_inc = 0;
    F_gap = 0;
    H_gap = -gap_open;
    G_gap = -gap_open;
    start_edge_pens_x = -1;
  }
  else {
    printf("scream: unknown gaps mode\n");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  
  /* process each row. i.e. each character of seq2 */

  b_s = big_neg;
  b_e = b_r = b_c = 0;
  t = 0;

  if ( band ) {
	/*
	 * If band is set we restrict the search to band elements either side
	 * of the main diagonal. This gives a band length of (2 * band) + 1.
	 * For the affine alignment it is necessary to know what happened in
	 * the elements to the left and above the current one, and therefore
	 * need to add another 2 to the band length.
	 */

	band_length = (2 * band) + 3;
	/*
	 * Need to add one to first_row, first_column, band_left and
	 * band_right because the alignment matrix has an extra row
	 * and an extra column over the number given by the lengths
	 * of the two sequences.
	 */
	size_mat = (MIN(seq1_len - band_left, seq2_len - first_row) + 1) 
	            * band_length;

	if(!(bit_trace = (unsigned char *) 
	    xmalloc(1 + sizeof(char) * size_mat))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, 
			     seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	first_row++;
	band_left++;
	band_right++;
	first_band_left = band_left;

	two_band = band * 2;
	max_row = MIN(seq2_len, first_row + seq1_len - band_left + 1);

	last_row = 0;

	for(row = first_row, e_row = band_length; row <= max_row; row++, band_left++, band_right++, e_row+=band_length) {

	  guard_offset = band_left + two_band;

	  if(t == 0) {
	    pF1   = F1;
	    pF2   = F2;
	    pG1   = G1;
	    pG2   = G2;
	    pH1   = H1;
	    pH2   = H2;
	    pF_guard = F1 + guard_offset;
	    pG_guard = G1 + guard_offset;
	    pH_guard = H1 + guard_offset;
	    F2[0] = F_gap;
	    G2[0] = G_gap;
	    H2[0] = H_gap;
	    F_gap -= edge_inc;
	    G_gap -= edge_inc;
	    if ( row > start_edge_pens_x ) {
		H_gap -= edge_inc;
	    }
	    else {
	      H_gap -= gap_to_gap;
	    }
	    t = 1;
	  } else {
	    pF1   = F2;
	    pF2   = F1;
	    pG1   = G2;
	    pG2   = G1;
	    pH1   = H2;
	    pH2   = H1;
	    pF_guard = F2 + guard_offset;
	    pG_guard = G2 + guard_offset;
	    pH_guard = H2 + guard_offset;
	    F1[0] = F_gap;
	    G1[0] = G_gap;
	    H1[0] = H_gap;
	    F_gap -= edge_inc;
	    G_gap -= edge_inc;
	    if ( row > start_edge_pens_x ) {
		H_gap -= edge_inc;
	    }
	    else {
	      H_gap -= gap_to_gap;
	    }
	    t = 0;
	  }
	  if ( (off_set = band_left - 1 ) > 0 ) {
	    pF1 += off_set;
	    pF2 += off_set;
	    pG1 += off_set;
	    pG2 += off_set;
	    pH1 += off_set;
	    pH2 += off_set;
	    *pF2 = big_neg;
	    *pG2 = big_neg;
	    *pH2 = big_neg;
	  }
	  t_pF2 = pF2;
	  t_pG2 = pG2;
	  t_pH2 = pH2;

	  if ( band_right <= seq1_len ) {
	    *pF_guard = big_neg;
	    *pG_guard = big_neg;
	    *pH_guard = big_neg;
	  }

	  score_matrix_p = W128[(int)seq2[row-1]];
	  if ( seq2[row-1] != OLD_PAD_SYM ) {
	    gap_open_y = gap_open;
	    gap_extend_y = gap_extend;
	  }
	  else {
	    gap_open_y = gap_extend_y = gap_to_gap;
	  }
	    /*
	     * Need to prevent the possibility of certain transitions
	     * at the ends of each band; column-1 should not be allowed at
	     * the left-hand end of a band, and row-1 should not be
	     * allowed at the right-hand end. Such transitions would
	     * take the trace-back out of the the parts of the 
	     * alignment matrix that are covered by the band.
	     *
	     *      0 1 2 3 4 5 6
	     *    0   _ _ _           For example, shouldn't consider
	     *    1  |_|_|_|_         column-1 at position (row = 3, col = 3)
	     *    2    |_|_|_|_       so give (3, 2) a very low score.
	     *    3      |_|_|_|_     Likewise, shouldn't consider row-1
	     *    4        |_|_|_|    at position (3, 5), so give (2, 5)
	     *    5          |_|_|    a very low score.
	     * this is done using the guard pointers at the right edge
	     * and initialising the left edges where required.   
	     */

	  /* process each column. i.e. each character of seq1 */
	  
	  max_col = MIN(seq1_len, band_right);
	  for(column = MAX(1, band_left),
	      e_col = MAX(1, band_left) - band_left + 1;
		column <= max_col;
		column++, pF1++, pF2++, pG1++, pG2++, pH1++, pH2++, e_col++) {

	      /* need to know last row we actually reach */

	      last_row = MAX(row,last_row);

	      /* move diagonally? */
	    s = score_matrix_p[(int)seq1[column-1]];
	    if ( seq1[column-1] != OLD_PAD_SYM ) {
	      gap_open_x = gap_open;
	      gap_extend_x = gap_extend;
	    }
	    else {
	      gap_open_x = gap_extend_x = gap_to_gap;
	    }

	    V_diag = *pF1 + s;
	    V_insx = *pH1 + s;
	    V_insy = *pG1 + s;
	    if ( V_insx > V_diag ) {
	      if ( V_insx > V_insy ) {
		best_F1 = V_insx;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    else {
	      if ( V_diag > V_insy ) {
		best_F1 = V_diag;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    *(pF2+1) = best_F1;

	    /* gap in x? */
	    V_diag =  *pF2 - gap_open_x;
	    V_extx =  *pH2 - gap_extend_x;
	    if ( V_diag > V_extx ) {
	      best_H1 = V_diag;
	    }
	    else {
	      best_H1 = V_extx;
	    }
	    *(pH2+1) = best_H1;

	    /* gap in y? */
	    V_diag = *(pF1+1) - gap_open_y;
	    V_exty = *(pG1+1) - gap_extend_y;
	    if ( V_diag > V_exty ) {
	      best_G1 = V_diag;
	    }
	    else {
	      best_G1 = V_exty;
	    }
	    *(pG2+1) = best_G1;

/*	    e_row = (row - first_row + 1) * band_length;*/
/*	    e_col = column - band_left + 1;*/
	    e = e_row + e_col;
	    
	    /* find the best move */
	    if ( best_H1 > best_F1 ) {
	      if ( best_H1 > best_G1 ) {
		bit_trace[e] = BYTE_ACROSS;
		b_s = best_H1;
	      }
	      else {
		bit_trace[e] = BYTE_DOWN;
		b_s = best_G1;
	      }
	    }
	    else {
	      if ( best_F1 > best_G1 ) {
		bit_trace[e] = BYTE_DIAGONAL;
		b_s = best_F1;
	      }
	      else {
		bit_trace[e] = BYTE_DOWN;
		b_s = best_G1;
	      }
	    }
	  }
	  
	  if ( column > seq1_len ) {
	    if ( edge_mode & BEST_EDGE_TRACE ) {
	      best_H1 = MAX(best_H1,best_G1);
	      best_F1 = MAX(best_H1,best_F1);
	      if ( best_F1 > best_edge_score ) {
		best_edge_score = best_F1;
		b_r = row;
		b_e = ((row - first_row + 1) * band_length) + 
		       (seq1_len - band_left + 1);
	      }
	    }
	  }
      }
	
	last_column = max_col;
        row = max_row;
	band_left--;
	band_right--;
	if ( edge_mode & BEST_EDGE_TRACE ) {
	  b_c = last_column;

	  pF2 = t_pF2+1;
	  pG2 = t_pG2+1;
	  pH2 = t_pH2+1;
	  max_col = MIN(seq1_len, band_right);
	  for(column = MAX(1, band_left);
	      column <= max_col;
	      column++, pF2++, pG2++, pH2++) {
	    best_F1 = *pF2;
	    best_G1 = *pG2;
	    best_H1 = *pH2;
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_c = column;
	      b_r = last_row;
	      b_e = ((row - first_row + 1) * band_length) + 
		     (column - band_left + 1);
	    }
	  }
      }
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	    b_r = last_row;
	    b_c = last_column;
	  b_e = ((b_r - first_row + 1) * band_length) + 
		 (b_c - band_left + 1);
	  best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	  b_r = last_row;
	  b_c = last_column;
	  b_e = ((b_r - first_row + 1) * band_length) + 
		 (b_c - band_left + 1);
	  best_edge_score = b_s;
        }

  }
  else {

	/* Initialise the bit trace */
	size_mat = (seq1_len + 1) * (seq2_len + 1);
	if(!(bit_trace = (unsigned char *) xmalloc(1 + sizeof(char) * size_mat))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	for(row = 1, e = seq1_len + 2; row <= seq2_len; row++, e++) {

	  if(t == 0) {
	    pF1   = F1;
	    pF2   = F2;
	    pG1   = G1;
	    pG2   = G2;
	    pH1   = H1;
	    pH2   = H2;
	    F2[0] = F_gap;
	    G2[0] = G_gap;
	    F_gap -= edge_inc;
	    G_gap -= edge_inc;
	    if ( row > start_edge_pens_x ) {
		H2[0] = H_gap;
		H_gap -= edge_inc;
	    }
	    else {
	      H2[0] = H_gap;
	      H_gap -= gap_to_gap;
	    }
	    t = 1;
	  } else {
	    pF1   = F2;
	    pF2   = F1;
	    pG1   = G2;
	    pG2   = G1;
	    pH1   = H2;
	    pH2   = H1;

	    F1[0] = F_gap;
	    G1[0] = G_gap;
	    F_gap -= edge_inc;
	    G_gap -= edge_inc;
	    if ( row > start_edge_pens_x ) {
		H1[0] = H_gap;
		H_gap -= edge_inc;
	    }
	    else {
	      H1[0] = g;
	      g -= gap_to_gap;
	    }
	    t = 0;
	  }
	  

	  score_matrix_p = W128[(int)seq2[row-1]];
	  if ( seq2[row-1] != OLD_PAD_SYM ) {
	    gap_open_y = gap_open;
	    gap_extend_y = gap_extend;
	  }
	  else {
	    gap_open_y = gap_extend_y = gap_to_gap;
	  }
	  
	  /* process each column. i.e. each character of seq1 */
	  
	  for(column = 1; column <= seq1_len; column++, e++) {

	    /* move diagonally? */
	    s = score_matrix_p[(int)seq1[column-1]];
	    if ( seq1[column-1] != OLD_PAD_SYM ) {
	      gap_open_x = gap_open;
	      gap_extend_x = gap_extend;
	    }
	    else {
	      gap_open_x = gap_extend_x = gap_to_gap;
	    }
	    V_diag = pF1[column-1] + s;
	    V_insx = pH1[column-1] + s;
	    V_insy = pG1[column-1] + s;
	    if ( V_insx > V_diag ) {
	      if ( V_insx > V_insy ) {
		best_F1 = V_insx;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    else {
	      if ( V_diag > V_insy ) {
		best_F1 = V_diag;
	      }
	      else {
		best_F1 = V_insy;
	      }
	    }
	    pF2[column] = best_F1;

	    /*printf("%3d %3d %3d %3d %3d ",row,column,V_diag,V_insx,V_insy);*/
	    
	    /* gap in x? */
	    V_diag =  pF2[column-1] - gap_open_x;
	    V_extx =  pH2[column-1] - gap_extend_x;
	    if ( V_diag > V_extx ) {
	      best_H1 = V_diag;
	    }
	    else {
	      best_H1 = V_extx;
	    }
	    pH2[column] = best_H1;

	    /*printf("%3d %3d ",V_diag,V_extx);*/
	    
	    /* gap in y? */
	    V_diag =  pF1[column] - gap_open_y;
	    V_exty = pG1[column] - gap_extend_y;
	    if ( V_diag > V_exty ) {
	      best_G1 = V_diag;
	    }
	    else {
	      best_G1 = V_exty;
	    }
	    pG2[column] = best_G1;

	    /*printf("%3d %3d %3d %3d %3d ",V_diag,V_exty,s,gap_open_x,gap_open_y);*/


	    /* choose best move */
	    if ( best_H1 > best_F1 ) {
	      if ( best_H1 > best_G1 ) {
		bit_trace[e] = BYTE_ACROSS;
		b_s = best_H1;
		/*printf("%3d -",b_s);*/
	      }
	      else {
		bit_trace[e] = BYTE_DOWN;
		b_s = best_G1;
		/*printf("%3d |",b_s);*/
	      }
	    }
	    else {
	      if ( best_F1 > best_G1 ) {
		bit_trace[e] = BYTE_DIAGONAL;
		b_s = best_F1;
		/*printf("%3d \\",b_s);*/
	      }
	      else {
		bit_trace[e] = BYTE_DOWN;
		b_s = best_G1;
		/*printf("%3d |",b_s);*/
	      }
	    }
	    /*printf("\n");*/

	  }
	  
	  if ( edge_mode & BEST_EDGE_TRACE ) {
	    /* best right edge score */
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_r = row;
	      b_e = (row + 1) * (seq1_len + 1) - 1;
	    }
	  }
	}
	
	if ( edge_mode & BEST_EDGE_TRACE ) {
	  /* best bottom edge score */
	  b_c = seq1_len;
	  for(column = 1; column <= seq1_len; column++) {
	    best_F1 = pF2[column];
	    best_G1 = pG2[column];
	    best_H1 = pH2[column];
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_c = column;
	      b_r = seq2_len;
	      b_e = (row - 1) * (seq1_len + 1) + column;
	    }
	  }
	}
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = seq2_len * (seq1_len + 1) + seq1_len;
	  best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = seq2_len * (seq1_len + 1) + seq1_len;
	  best_edge_score = b_s;
        }
    }


  /* do traceback */

  overlap->score = best_edge_score;

  if( i = do_trace_back ( bit_trace, seq1, seq2, seq1_len, seq2_len,
		      &seq1_out, &seq2_out, &seq_out_len, b_r, b_c, b_e,
		      band, band ? first_band_left : 0, first_row,
			  band ? band_length : 0, NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
      return -1;
  }

  overlap->seq1_out = seq1_out;
  overlap->seq2_out = seq2_out;
  overlap->seq_out_len = seq_out_len;

  if ( i = seq_to_overlap (overlap, OLD_PAD_SYM, NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
  }

  if ( params->job & RETURN_EDIT_BUFFERS ) {
    if (seq_to_edit ( seq1_out,seq_out_len,&overlap->S1,&overlap->s1_len,NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
    }
    if (seq_to_edit ( seq2_out,seq_out_len,&overlap->S2,&overlap->s2_len,NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
    }
  }

  if ( params->job & RETURN_SEQ ) {
    if ( !(params->job & RETURN_NEW_PADS) ) {
      old_pads_for_new(seq1_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
      old_pads_for_new(seq2_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
    }
    seq1_out = seq2_out = NULL; /* stop them being freed! */
  }
  else {
    overlap->seq1_out = overlap->seq2_out = NULL;
    /* ie we let destroy_af_mem free the memory, but we must
     * ensure that othr routines do not try to free it too 
     */
  }
  destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );

  return 0;
}

int affine_align2_bits(OVERLAP *overlap, ALIGN_PARAMS *params) {

  char *seq1, *seq2;
  int seq1_len, seq2_len, seq_out_len;
  int gap_open, gap_extend, edge_inc;
  int i,j;
  int s,*score_matrix_p;
  char *seq1_out, *seq2_out;
  int b_c, b_r;

  int t,big_neg,b_s,e,b_e;
  int *F1, *F2, *G1, *G2, *H1, *H2;
  int *pF1, *pF2, *pG1, *pG2;
  int *t_pF2, *t_pG2;
  int best_F1, best_G1, best_H1, FV_diag, GV_diag,FV_insx, GV_insx, FV_insy, GV_insy;
  int E_gap, F_gap;
  int edge_mode, best_edge_score;

  int band, band_length, two_band, band_left, band_right, first_band_left=0;
  int off_set, guard_offset, *pF_guard, *pG_guard;
  int row, first_row, max_row, column, max_col;
  unsigned char *bit_trace;
  int byte, nibble, e_row, e_col, size_mat;
  char OLD_PAD_SYM, NEW_PAD_SYM;
  int gap_open_x, gap_open_y, gap_extend_x, gap_extend_y, gap_to_gap;

  /*
   *    Three possible alignment cases:
   *    IGAxi   AIGAHxi   GAxi--
   *    LGVyj   GVyj--    SLGVHyj
   *       F      G            H
   *    i.e. xi aligned with yj, xi aligned opposite a gap y,
   *    or yi aligned opposite a gap in x
   *    below these cases are contained in the recurrence relations
   *    for F, G and H respectively
   *    s(xi,yj) is score matrix
   *    d is gap_open
   *    e is gap extend
   *
   *                   F(i-1,j-1) + s(xi,yi)
   *    F(i,j)  = max  G(i-1,j-1) + s(xi,yi)      \  no gap
   *                   FV_diag
   *                   GV_diag
   *
   *                   F(i,j-1)   - d
   *    G(i,j) = max   G(i,j-1)   - e             |  gap in y
   *                   F(i-1,j)   - d
   *                   G(i-1,j)   - e             -  gap in x
   *                   FV_insy
   *                   GV_insy
   *                   FV_insx
   *                   GV_insx
   *    Find MAX(F(i,j),G(i,j)) and set trace accordingly:
   *                \     |      -
   *
   *    if gaps at left edge count:
   *    G(0,i) = G(i,0) = - d - e * i
   *    F(1,i) = F(i,1) = - d - e * i;
   *    F(0,0) = 0;
   *    otherwise all set to 0;
   *    if right end gaps count the best score is at (seq1_len,seq2_len)
   *    otherwise find the best score along the two edges
   *    trace back accordingly
   *
   *    store 2 rows for each of F, G, H
   *    use p_F1, p_G1, p_G1 to point to previous row
   *    p_F2, p_G2, p_H2 for current row being built
   *    at the start of a new row:
   *
   *    rows have length seq1_len, columns seq2_len
   *    i.e.
   *    rows: 1 - seq1_len, columns 1 - seq2_len
   *    seq1xxxxxxxxxxxxxxx
   *   s
   *   e
   *   q
   *   2
   *   y
   *   y
   *   y
   *
   *
   */

  F1 = F2 = G1 = G2 = H1 = H2 = NULL;
  bit_trace = NULL;
  seq1_out = seq2_out = NULL;
  big_neg = INT_MIN/2;
  best_edge_score = big_neg;

  seq1 = overlap->seq1;
  seq2 = overlap->seq2;
  seq1_len = overlap->seq1_len;
  seq2_len = overlap->seq2_len;

  edge_mode = params->edge_mode;
  gap_open = params->gap_open;
  gap_extend = params->gap_extend;
  OLD_PAD_SYM = params->old_pad_sym;
  NEW_PAD_SYM = params->new_pad_sym;
  band = params->band;
  first_row = params->first_row;
  band_left = params->band_left;
  band_right = params->band_right;
  edge_inc = gap_extend;
  gap_to_gap = -W128[(int)OLD_PAD_SYM][(int)OLD_PAD_SYM];
  gap_to_gap = 1;

  /* init tables */

  if(!(F1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for F1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(F2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for F2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(G1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for G1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(G2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for G2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }

  /* do recurrence */

  if ( edge_mode & EDGE_GAPS_COUNT ) {
    F1[0] = 0;
    for(i = 1,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) F1[i] = j;
    for(i = 0,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) G1[i] = j;
    E_gap = -gap_open - edge_inc;
    F_gap = -gap_open;
  }
  else if ( edge_mode & EDGE_GAPS_ZERO ) {
    for(i = 0; i <= seq1_len; i++) F1[i] = 0;
    for(i = 0; i <= seq1_len; i++) G1[i] = 0;
    edge_inc = 0;
    E_gap = 0;
    F_gap = 0;
  }
  else {
    printf("scream: unknown gaps mode\n");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  
  /* process each row. i.e. each character of seq2 */

  b_s = big_neg;
  b_e = b_r = b_c = 0;
  t = 0;

  if ( band ) {
	/*
	 * If band is set we restrict the search to band elements either side
	 * of the main diagonal. This gives a band length of (2 * band) + 1.
	 * For the affine alignment it is necessary to know what happened in
	 * the elements to the left and above the current one, and therefore
	 * need to add another 2 to the band length.
	 */

	band_length = (2 * band) + 3;
	/*
	 * Need to add one to first_row, first_column, band_left and
	 * band_right because the alignment matrix has an extra row
	 * and an extra column over the number given by the lengths
	 * of the two sequences.
	 */
	size_mat = (MIN(seq1_len - band_left, seq2_len - first_row) + 1) 
	            * band_length;

	if(!(bit_trace = (unsigned char *) 
	    xmalloc(1 + sizeof(char) * size_mat / 4))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, 
			     seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat / 4;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	first_row++;
	band_left++;
	band_right++;
	first_band_left = band_left;

	two_band = band * 2;
	max_row = MIN(seq2_len, first_row + seq1_len - band_left + 1);

	for(row = first_row; row <= max_row; row++, band_left++, band_right++) {

	  guard_offset = band_left + two_band;

	  if(t == 0) {
	    pF1   = F1;
	    pF2   = F2;
	    pG1   = G1;
	    pG2   = G2;
	    pF_guard = F1 + guard_offset;
	    pG_guard = G1 + guard_offset;
	    F2[0] = F_gap;
	    G2[0] = E_gap;
	    t = 1;
	  } else {
	    pF1   = F2;
	    pF2   = F1;
	    pG1   = G2;
	    pG2   = G1;
	    pF_guard = F2 + guard_offset;
	    pG_guard = G2 + guard_offset;
	    F1[0] = F_gap;
	    G1[0] = E_gap;
	    t = 0;
	  }
	  if ( (off_set = band_left - 1 ) > 0 ) {
	    pF1 += off_set;
	    pF2 += off_set;
	    pG1 += off_set;
	    pG2 += off_set;
	    *pF2 = big_neg;
	    *pG2 = big_neg;
	  }
	  t_pF2 = pF2;
	  t_pG2 = pG2;

	  if ( band_right <= seq1_len ) {
	    *pF_guard = big_neg;
	    *pG_guard = big_neg;
	  }
	  E_gap -= edge_inc;
	  F_gap -= edge_inc;

	  score_matrix_p = W128[(int)seq2[row-1]];
	  if ( seq2[row-1] != OLD_PAD_SYM ) {
	    gap_open_y = gap_open;
	    gap_extend_y = gap_extend;
	  }
	  else {
	    gap_open_y = gap_extend_y = gap_to_gap;
	  }
	    /*
	     * Need to prevent the possibility of certain transitions
	     * at the ends of each band; column-1 should not be allowed at
	     * the left-hand end of a band, and row-1 should not be
	     * allowed at the right-hand end. Such transitions would
	     * take the trace-back out of the the parts of the 
	     * alignment matrix that are covered by the band.
	     *
	     *      0 1 2 3 4 5 6
	     *    0   _ _ _           For example, shouldn't consider
	     *    1  |_|_|_|_         column-1 at position (row = 3, col = 3)
	     *    2    |_|_|_|_       so give (3, 2) a very low score.
	     *    3      |_|_|_|_     Likewise, shouldn't consider row-1
	     *    4        |_|_|_|    at position (3, 5), so give (2, 5)
	     *    5          |_|_|    a very low score.
	     * this is done using the guard pointers at the right edge
	     * and initialising the left edges where required.   
	     */

	  /* process each column. i.e. each character of seq1 */
	  
	  max_col = MIN(seq1_len, band_right);
	  for(column = MAX(1, band_left);
		column <= max_col;
		column++, pF1++, pF2++, pG1++, pG2++) {

	    /* move diagonally? */
	    s = score_matrix_p[(int)seq1[column-1]];
	    if ( seq1[column-1] != OLD_PAD_SYM ) {
	      gap_open_x = gap_open;
	      gap_extend_x = gap_extend;
	    }
	    else {
	      gap_open_x = gap_extend_x = gap_to_gap;
	    }

	    FV_diag = *pF1 + s;
	    GV_diag = *pG1 + s;
	    if ( GV_diag > FV_diag ) {
		best_F1 = GV_diag;
	    }
	    else {
	      best_F1 = FV_diag;
	    }
	    *(pF2+1) = best_F1;

	    /* gap in x? */
	    FV_insx =  *pF2 - gap_open_x;
	    GV_insx =  *pG2 - gap_extend_x;
	    /* gap in y? */
	    FV_insy = *(pF1+1) - gap_open_y;
	    GV_insy = *(pG1+1) - gap_extend_y;
	    best_H1 = MAX(FV_insx,GV_insx);
	    best_G1 = MAX(FV_insy,GV_insy);
	    *(pG2+1) = MAX(best_H1,best_G1);

	    e_row = (row - first_row + 1) * band_length;
	    e_col = column - band_left + 1;
	    e = e_row + e_col;
	    byte = e / 4;
	    nibble = 2 * (e % 4);
	    
	    /* find the best move */
	    if ( best_H1 > best_F1 ) {
	      if ( best_H1 > best_G1 ) {
		bit_trace[byte] |= BYTE_ACROSS << nibble;
		b_s = best_H1;
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
	      }
	    }
	    else {
	      if ( best_F1 > best_G1 ) {
		bit_trace[byte] |= BYTE_DIAGONAL << nibble;
		b_s = best_F1;
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
	      }
	    }
	  }
	  
	  if ( column > seq1_len ) {
	    if ( edge_mode & BEST_EDGE_TRACE ) {
	      best_H1 = MAX(best_H1,best_G1);
	      best_F1 = MAX(best_H1,best_F1);
	      if ( best_F1 > best_edge_score ) {
		best_edge_score = best_F1;
		b_r = row;
		b_e = ((row - first_row + 1) * band_length) + 
		       (seq1_len - band_left + 1);
	      }
	    }
	  }
      }
	
	
        row--;
	band_left--;
	band_right--;
	if ( edge_mode & BEST_EDGE_TRACE ) {
	  b_c = seq1_len;

	  pF2 = t_pF2+1;
	  pG2 = t_pG2+1;
	  max_col = MIN(seq1_len, band_right);
	  for(column = MAX(1, band_left);
	      column <= max_col;
	      column++, pF2++, pG2++) {
	    best_F1 = *pF2;
	    best_G1 = *pG2;
	    best_F1 = MAX(best_G1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_c = column;
	      b_r = seq2_len;
	      b_e = ((row - first_row + 1) * band_length) + 
		     (column - band_left + 1);
	    }
	  }
      }
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	    b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = ((b_r - first_row + 1) * band_length) + 
		 (b_c - band_left + 1);
	  best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = ((b_r - first_row + 1) * band_length) + 
		 (b_c - band_left + 1);
	  best_edge_score = b_s;
        }

  }
  else {

	/* Initialise the bit trace */
	size_mat = (seq1_len + 1) * (seq2_len + 1);
	if(!(bit_trace = (unsigned char *) xmalloc(1 + sizeof(char) * size_mat / 4))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat / 4;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	for(row = 1, e = seq1_len + 2; row <= seq2_len; row++, e++) {

	  if(t == 0) {
	    pF1   = F1;
	    pF2   = F2;
	    pG1   = G1;
	    pG2   = G2;
	    F2[0] = F_gap;
	    G2[0] = E_gap;
	    t = 1;
	  } else {
	    pF1   = F2;
	    pF2   = F1;
	    pG1   = G2;
	    pG2   = G1;
	    F1[0] = F_gap;
	    G1[0] = E_gap;
	    t = 0;
	  }
	  
	  E_gap -= edge_inc;
	  F_gap -= edge_inc;

	  score_matrix_p = W128[(int)seq2[row-1]];
	  if ( seq2[row-1] != OLD_PAD_SYM ) {
	    gap_open_y = gap_open;
	    gap_extend_y = gap_extend;
	  }
	  else {
	    gap_open_y = gap_extend_y = gap_to_gap;
	  }
	  
	  /* process each column. i.e. each character of seq1 */
	  
	  for(column = 1; column <= seq1_len; column++, e++) {
	    /* move diagonally? */
	    s = score_matrix_p[(int)seq1[column-1]];
	    if ( seq1[column-1] != OLD_PAD_SYM ) {
	      gap_open_x = gap_open;
	      gap_extend_x = gap_extend;
	    }
	    else {
	      gap_open_x = gap_extend_x = gap_to_gap;
	    }
	    FV_diag = pF1[column-1] + s;
	    GV_diag = pG1[column-1] + s;
	    if ( GV_diag > FV_diag ) {
		best_F1 = GV_diag;
	    }
	    else {
	      best_F1 = FV_diag;
	    }
	    pF2[column] = best_F1;
	    
	    /* gap in x? */
	    FV_insx =  pF2[column-1] - gap_open_x;
	    GV_insx =  pG2[column-1] - gap_extend_x;
	    /* gap in y? */
	    FV_insy = pF1[column] - gap_open_y;
	    GV_insy = pG1[column] - gap_extend_y;
	    best_H1 = MAX(FV_insx,GV_insx);
	    best_G1 = MAX(FV_insy,GV_insy);

	    pG2[column] = MAX(best_H1,best_G1);

	    byte = e / 4;
	    nibble = 2 * (e % 4);

	    /* choose best move */
	    if ( best_H1 > best_F1 ) {
	      if ( best_H1 > best_G1 ) {
		bit_trace[byte] |= BYTE_ACROSS << nibble;
		b_s = best_H1;
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
	      }
	    }
	    else {
	      if ( best_F1 > best_G1 ) {
		bit_trace[byte] |= BYTE_DIAGONAL << nibble;
		b_s = best_F1;
	      }
	      else {
		bit_trace[byte] |= BYTE_DOWN << nibble;
		b_s = best_G1;
	      }
	    }
	  }
	  
	  if ( edge_mode & BEST_EDGE_TRACE ) {
	    /* best right edge score */
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_r = row;
	      b_e = (row + 1) * (seq1_len + 1) - 1;
	    }
	  }
	}
	
	if ( edge_mode & BEST_EDGE_TRACE ) {
	  /* best bottom edge score */
	  b_c = seq1_len;
	  for(column = 1; column <= seq1_len; column++) {
	    best_F1 = pF2[column];
	    best_G1 = pG2[column];
	    best_F1 = MAX(best_G1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_c = column;
	      b_r = seq2_len;
	      b_e = (row - 1) * (seq1_len + 1) + column;
	    }
	  }
	}
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = seq2_len * (seq1_len + 1) + seq1_len;
	  best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = seq2_len * (seq1_len + 1) + seq1_len;
	  best_edge_score = b_s;
        }
    }


  /* do traceback */

  overlap->score = best_edge_score;

  if( i = do_trace_back_bits ( bit_trace, seq1, seq2, seq1_len, seq2_len,
		      &seq1_out, &seq2_out, &seq_out_len, b_r, b_c, b_e,
		      band, first_band_left, first_row, band_length, NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
      return -1;
  }

  overlap->seq1_out = seq1_out;
  overlap->seq2_out = seq2_out;
  overlap->seq_out_len = seq_out_len;

  if ( i = seq_to_overlap (overlap, OLD_PAD_SYM, NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
  }

  if ( params->job & RETURN_EDIT_BUFFERS ) {
    if (seq_to_edit ( seq1_out,seq_out_len,&overlap->S1,&overlap->s1_len,NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
    }
    if (seq_to_edit ( seq2_out,seq_out_len,&overlap->S2,&overlap->s2_len,NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
    }
  }

  if ( params->job & RETURN_SEQ ) {
    if ( !(params->job & RETURN_NEW_PADS) ) {
      old_pads_for_new(seq1_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
      old_pads_for_new(seq2_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
    }
    seq1_out = seq2_out = NULL; /* stop them being freed! */
  }
  else {
    overlap->seq1_out = overlap->seq2_out = NULL;
    /* ie we let destroy_af_mem free the memory, but we must
     * ensure that othr routines do not try to free it too 
     */
  }
  destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );

  return 0;
}

int affine_align(OVERLAP *overlap, ALIGN_PARAMS *params) {

    /* decide which algorithm to use */
#define MAX_MEMORY 10000000

    int mem;

    if (params->band) {
	mem = 2 * params->band * MIN(overlap->seq1_len,overlap->seq2_len);
    }
    else {
	mem = overlap->seq1_len * overlap->seq2_len;
    }
    if (mem > MAX_MEMORY) {
	return affine_align_bits(overlap,params);
    }
    return affine_align_big(overlap,params);
}

int affine_align2_big(OVERLAP *overlap, ALIGN_PARAMS *params) {

  char *seq1, *seq2;
  int seq1_len, seq2_len, seq_out_len;
  int gap_open, gap_extend, edge_inc;
  int i,j;
  int s,*score_matrix_p;
  char *seq1_out, *seq2_out;
  int b_c, b_r;

  int t,big_neg,b_s,e,b_e;
  int *F1, *F2, *G1, *G2, *H1, *H2;
  int *pF1, *pF2, *pG1, *pG2;
  int *t_pF2, *t_pG2;
  int best_F1, best_G1, best_H1, FV_diag, GV_diag,FV_insx, GV_insx, FV_insy, GV_insy;
  int E_gap, F_gap;
  int edge_mode, best_edge_score;

  int band, band_length, two_band, band_left, band_right, first_band_left=0;
  int off_set, guard_offset, *pF_guard, *pG_guard;
  int row, first_row, max_row, column, max_col;
  unsigned char *bit_trace;
  int e_row, e_col, size_mat;
  char OLD_PAD_SYM, NEW_PAD_SYM;
  int gap_open_x, gap_open_y, gap_extend_x, gap_extend_y, gap_to_gap;

  /*
   *    Three possible alignment cases:
   *    IGAxi   AIGAHxi   GAxi--
   *    LGVyj   GVyj--    SLGVHyj
   *       F      G            H
   *    i.e. xi aligned with yj, xi aligned opposite a gap y,
   *    or yi aligned opposite a gap in x
   *    below these cases are contained in the recurrence relations
   *    for F, G and H respectively
   *    s(xi,yj) is score matrix
   *    d is gap_open
   *    e is gap extend
   *
   *                   F(i-1,j-1) + s(xi,yi)
   *    F(i,j)  = max  G(i-1,j-1) + s(xi,yi)      \  no gap
   *                   FV_diag
   *                   GV_diag
   *
   *                   F(i,j-1)   - d
   *    G(i,j) = max   G(i,j-1)   - e             |  gap in y
   *                   F(i-1,j)   - d
   *                   G(i-1,j)   - e             -  gap in x
   *                   FV_insy
   *                   GV_insy
   *                   FV_insx
   *                   GV_insx
   *    Find MAX(F(i,j),G(i,j)) and set trace accordingly:
   *                \     |      -
   *
   *    if gaps at left edge count:
   *    G(0,i) = G(i,0) = - d - e * i
   *    F(1,i) = F(i,1) = - d - e * i;
   *    F(0,0) = 0;
   *    otherwise all set to 0;
   *    if right end gaps count the best score is at (seq1_len,seq2_len)
   *    otherwise find the best score along the two edges
   *    trace back accordingly
   *
   *    store 2 rows for each of F, G, H
   *    use p_F1, p_G1, p_G1 to point to previous row
   *    p_F2, p_G2, p_H2 for current row being built
   *    at the start of a new row:
   *
   *    rows have length seq1_len, columns seq2_len
   *    i.e.
   *    rows: 1 - seq1_len, columns 1 - seq2_len
   *    seq1xxxxxxxxxxxxxxx
   *   s
   *   e
   *   q
   *   2
   *   y
   *   y
   *   y
   *
   *
   */

  F1 = F2 = G1 = G2 = H1 = H2 = NULL;
  bit_trace = NULL;
  seq1_out = seq2_out = NULL;
  big_neg = INT_MIN/2;
  best_edge_score = big_neg;

  seq1 = overlap->seq1;
  seq2 = overlap->seq2;
  seq1_len = overlap->seq1_len;
  seq2_len = overlap->seq2_len;

  edge_mode = params->edge_mode;
  gap_open = params->gap_open;
  gap_extend = params->gap_extend;
  OLD_PAD_SYM = params->old_pad_sym;
  NEW_PAD_SYM = params->new_pad_sym;
  band = params->band;
  first_row = params->first_row;
  band_left = params->band_left;
  band_right = params->band_right;
  edge_inc = gap_extend;
  gap_to_gap = -W128[(int)OLD_PAD_SYM][(int)OLD_PAD_SYM];
  gap_to_gap = 1;

  /* init tables */

  if(!(F1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for F1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(F2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for F2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(G1 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for G1");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  if(!(G2 = (int *) xmalloc(sizeof(int) * (seq1_len + 2)))) {
    verror(ERR_WARN, "affine_align", "xmalloc failed for G2");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }

  /* do recurrence */

  if ( edge_mode & EDGE_GAPS_COUNT ) {
    F1[0] = 0;
    for(i = 1,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) F1[i] = j;
    for(i = 0,j=-gap_open; i <= seq1_len; i++,j-=edge_inc) G1[i] = j;
    E_gap = -gap_open - edge_inc;
    F_gap = -gap_open;
  }
  else if ( edge_mode & EDGE_GAPS_ZERO ) {
    for(i = 0; i <= seq1_len; i++) F1[i] = 0;
    for(i = 0; i <= seq1_len; i++) G1[i] = 0;
    edge_inc = 0;
    E_gap = 0;
    F_gap = 0;
  }
  else {
    printf("scream: unknown gaps mode\n");
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
    return -1;
  }
  
  /* process each row. i.e. each character of seq2 */

  b_s = big_neg;
  b_e = b_r = b_c = 0;
  t = 0;

  if ( band ) {
	/*
	 * If band is set we restrict the search to band elements either side
	 * of the main diagonal. This gives a band length of (2 * band) + 1.
	 * For the affine alignment it is necessary to know what happened in
	 * the elements to the left and above the current one, and therefore
	 * need to add another 2 to the band length.
	 */

	band_length = (2 * band) + 3;
	/*
	 * Need to add one to first_row, first_column, band_left and
	 * band_right because the alignment matrix has an extra row
	 * and an extra column over the number given by the lengths
	 * of the two sequences.
	 */
	size_mat = (MIN(seq1_len - band_left, seq2_len - first_row) + 1) 
	            * band_length;

	if(!(bit_trace = (unsigned char *) 
	    xmalloc(1 + sizeof(char) * size_mat))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, 
			     seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	first_row++;
	band_left++;
	band_right++;
	first_band_left = band_left;

	two_band = band * 2;
	max_row = MIN(seq2_len, first_row + seq1_len - band_left + 1);

	for(row = first_row, e_row = band_length; row <= max_row; row++, band_left++, band_right++, e_row+=band_length) {

	  guard_offset = band_left + two_band;

	  if(t == 0) {
	    pF1   = F1;
	    pF2   = F2;
	    pG1   = G1;
	    pG2   = G2;
	    pF_guard = F1 + guard_offset;
	    pG_guard = G1 + guard_offset;
	    F2[0] = F_gap;
	    G2[0] = E_gap;
	    t = 1;
	  } else {
	    pF1   = F2;
	    pF2   = F1;
	    pG1   = G2;
	    pG2   = G1;
	    pF_guard = F2 + guard_offset;
	    pG_guard = G2 + guard_offset;
	    F1[0] = F_gap;
	    G1[0] = E_gap;
	    t = 0;
	  }
	  if ( (off_set = band_left - 1 ) > 0 ) {
	    pF1 += off_set;
	    pF2 += off_set;
	    pG1 += off_set;
	    pG2 += off_set;
	    *pF2 = big_neg;
	    *pG2 = big_neg;
	  }
	  t_pF2 = pF2;
	  t_pG2 = pG2;

	  if ( band_right <= seq1_len ) {
	    *pF_guard = big_neg;
	    *pG_guard = big_neg;
	  }
	  E_gap -= edge_inc;
	  F_gap -= edge_inc;

	  score_matrix_p = W128[(int)seq2[row-1]];
	  if ( seq2[row-1] != OLD_PAD_SYM ) {
	    gap_open_y = gap_open;
	    gap_extend_y = gap_extend;
	  }
	  else {
	    gap_open_y = gap_extend_y = gap_to_gap;
	  }
	    /*
	     * Need to prevent the possibility of certain transitions
	     * at the ends of each band; column-1 should not be allowed at
	     * the left-hand end of a band, and row-1 should not be
	     * allowed at the right-hand end. Such transitions would
	     * take the trace-back out of the the parts of the 
	     * alignment matrix that are covered by the band.
	     *
	     *      0 1 2 3 4 5 6
	     *    0   _ _ _           For example, shouldn't consider
	     *    1  |_|_|_|_         column-1 at position (row = 3, col = 3)
	     *    2    |_|_|_|_       so give (3, 2) a very low score.
	     *    3      |_|_|_|_     Likewise, shouldn't consider row-1
	     *    4        |_|_|_|    at position (3, 5), so give (2, 5)
	     *    5          |_|_|    a very low score.
	     * this is done using the guard pointers at the right edge
	     * and initialising the left edges where required.   
	     */

	  /* process each column. i.e. each character of seq1 */
	  
	  max_col = MIN(seq1_len, band_right);
	  for(column = MAX(1, band_left), 
	      e_col = MAX(1, band_left) - band_left + 1;
		column <= max_col;
		column++, pF1++, pF2++, pG1++, pG2++, e_col++) {

	    /* move diagonally? */
	    s = score_matrix_p[(int)seq1[column-1]];
	    if ( seq1[column-1] != OLD_PAD_SYM ) {
	      gap_open_x = gap_open;
	      gap_extend_x = gap_extend;
	    }
	    else {
	      gap_open_x = gap_extend_x = gap_to_gap;
	    }

	    FV_diag = *pF1 + s;
	    GV_diag = *pG1 + s;
	    if ( GV_diag > FV_diag ) {
		best_F1 = GV_diag;
	    }
	    else {
	      best_F1 = FV_diag;
	    }
	    *(pF2+1) = best_F1;

	    /* gap in x? */
	    FV_insx =  *pF2 - gap_open_x;
	    GV_insx =  *pG2 - gap_extend_x;
	    /* gap in y? */
	    FV_insy = *(pF1+1) - gap_open_y;
	    GV_insy = *(pG1+1) - gap_extend_y;
	    best_H1 = MAX(FV_insx,GV_insx);
	    best_G1 = MAX(FV_insy,GV_insy);
	    *(pG2+1) = MAX(best_H1,best_G1);

/*	    e_row = (row - first_row + 1) * band_length;*/
/*	    e_col = column - band_left + 1;*/
	    e = e_row + e_col;
	    
	    /* find the best move */
	    if ( best_H1 > best_F1 ) {
	      if ( best_H1 > best_G1 ) {
		bit_trace[e] = BYTE_ACROSS;
		b_s = best_H1;
	      }
	      else {
		bit_trace[e] = BYTE_DOWN;
		b_s = best_G1;
	      }
	    }
	    else {
	      if ( best_F1 > best_G1 ) {
		bit_trace[e] = BYTE_DIAGONAL;
		b_s = best_F1;
	      }
	      else {
		bit_trace[e] = BYTE_DOWN;
		b_s = best_G1;
	      }
	    }
	  }
	  
	  if ( column > seq1_len ) {
	    if ( edge_mode & BEST_EDGE_TRACE ) {
	      best_H1 = MAX(best_H1,best_G1);
	      best_F1 = MAX(best_H1,best_F1);
	      if ( best_F1 > best_edge_score ) {
		best_edge_score = best_F1;
		b_r = row;
		b_e = ((row - first_row + 1) * band_length) + 
		    (seq1_len - band_left + 1);
	    }
	  }
	}
      }
	
	
        row--;
	band_left--;
	band_right--;
	if ( edge_mode & BEST_EDGE_TRACE ) {
	  b_c = seq1_len;

	  pF2 = t_pF2+1;
	  pG2 = t_pG2+1;
	  max_col = MIN(seq1_len, band_right);
	  for(column = MAX(1, band_left);
	      column <= max_col;
	      column++, pF2++, pG2++) {
	    best_F1 = *pF2;
	    best_G1 = *pG2;
	    best_F1 = MAX(best_G1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_c = column;
	      b_r = seq2_len;
	      b_e = ((row - first_row + 1) * band_length) + 
		     (column - band_left + 1);
	    }
	  }
      }
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	    b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = ((b_r - first_row + 1) * band_length) + 
		 (b_c - band_left + 1);
	  best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = ((b_r - first_row + 1) * band_length) + 
		 (b_c - band_left + 1);
	  best_edge_score = b_s;
        }

  }
  else {

	/* Initialise the bit trace */
	size_mat = (seq1_len + 1) * (seq2_len + 1);
	if(!(bit_trace = (unsigned char *) xmalloc(1 + sizeof(char) * size_mat))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	for(row = 1, e = seq1_len + 2; row <= seq2_len; row++, e++) {

	  if(t == 0) {
	    pF1   = F1;
	    pF2   = F2;
	    pG1   = G1;
	    pG2   = G2;
	    F2[0] = F_gap;
	    G2[0] = E_gap;
	    t = 1;
	  } else {
	    pF1   = F2;
	    pF2   = F1;
	    pG1   = G2;
	    pG2   = G1;
	    F1[0] = F_gap;
	    G1[0] = E_gap;
	    t = 0;
	  }
	  
	  E_gap -= edge_inc;
	  F_gap -= edge_inc;

	  score_matrix_p = W128[(int)seq2[row-1]];
	  if ( seq2[row-1] != OLD_PAD_SYM ) {
	    gap_open_y = gap_open;
	    gap_extend_y = gap_extend;
	  }
	  else {
	    gap_open_y = gap_extend_y = gap_to_gap;
	  }
	  
	  /* process each column. i.e. each character of seq1 */
	  
	  for(column = 1; column <= seq1_len; column++, e++) {
	    /* move diagonally? */
	    s = score_matrix_p[(int)seq1[column-1]];
	    if ( seq1[column-1] != OLD_PAD_SYM ) {
	      gap_open_x = gap_open;
	      gap_extend_x = gap_extend;
	    }
	    else {
	      gap_open_x = gap_extend_x = gap_to_gap;
	    }
	    FV_diag = pF1[column-1] + s;
	    GV_diag = pG1[column-1] + s;
	    if ( GV_diag > FV_diag ) {
		best_F1 = GV_diag;
	    }
	    else {
	      best_F1 = FV_diag;
	    }
	    pF2[column] = best_F1;
	    
	    /* gap in x? */
	    FV_insx =  pF2[column-1] - gap_open_x;
	    GV_insx =  pG2[column-1] - gap_extend_x;
	    /* gap in y? */
	    FV_insy = pF1[column] - gap_open_y;
	    GV_insy = pG1[column] - gap_extend_y;
	    best_H1 = MAX(FV_insx,GV_insx);
	    best_G1 = MAX(FV_insy,GV_insy);

	    pG2[column] = MAX(best_H1,best_G1);

	    /* choose best move */
	    if ( best_H1 > best_F1 ) {
	      if ( best_H1 > best_G1 ) {
		bit_trace[e] = BYTE_ACROSS;
		b_s = best_H1;
	      }
	      else {
		bit_trace[e] = BYTE_DOWN;
		b_s = best_G1;
	      }
	    }
	    else {
	      if ( best_F1 > best_G1 ) {
		bit_trace[e] = BYTE_DIAGONAL;
		b_s = best_F1;
	      }
	      else {
		bit_trace[e] = BYTE_DOWN;
		b_s = best_G1;
	      }
	    }
	  }
	  
	  if ( edge_mode & BEST_EDGE_TRACE ) {
	    /* best right edge score */
	    best_H1 = MAX(best_H1,best_G1);
	    best_F1 = MAX(best_H1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_r = row;
	      b_e = (row + 1) * (seq1_len + 1) - 1;
	    }
	  }
	}
	
	if ( edge_mode & BEST_EDGE_TRACE ) {
	  /* best bottom edge score */
	  b_c = seq1_len;
	  for(column = 1; column <= seq1_len; column++) {
	    best_F1 = pF2[column];
	    best_G1 = pG2[column];
	    best_F1 = MAX(best_G1,best_F1);
	    if ( best_F1 > best_edge_score ) {
	      best_edge_score = best_F1;
	      b_c = column;
	      b_r = seq2_len;
	      b_e = (row - 1) * (seq1_len + 1) + column;
	    }
	  }
	}
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = seq2_len * (seq1_len + 1) + seq1_len;
	  best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	  b_r = seq2_len;
	  b_c = seq1_len;
	  b_e = seq2_len * (seq1_len + 1) + seq1_len;
	  best_edge_score = b_s;
        }
    }


  /* do traceback */

  overlap->score = best_edge_score;

  if( i = do_trace_back ( bit_trace, seq1, seq2, seq1_len, seq2_len,
		      &seq1_out, &seq2_out, &seq_out_len, b_r, b_c, b_e,
		      band, first_band_left, first_row, band_length, NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
      return -1;
  }

  overlap->seq1_out = seq1_out;
  overlap->seq2_out = seq2_out;
  overlap->seq_out_len = seq_out_len;

  if ( i = seq_to_overlap (overlap, OLD_PAD_SYM, NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
  }

  if ( params->job & RETURN_EDIT_BUFFERS ) {
    if (seq_to_edit ( seq1_out,seq_out_len,&overlap->S1,&overlap->s1_len,NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
    }
    if (seq_to_edit ( seq2_out,seq_out_len,&overlap->S2,&overlap->s2_len,NEW_PAD_SYM)) {
      destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
      return -1;
    }
  }

  if ( params->job & RETURN_SEQ ) {
    if ( !(params->job & RETURN_NEW_PADS) ) {
      old_pads_for_new(seq1_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
      old_pads_for_new(seq2_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
    }
    seq1_out = seq2_out = NULL; /* stop them being freed! */
  }
  else {
    overlap->seq1_out = overlap->seq2_out = NULL;
    /* ie we let destroy_af_mem free the memory, but we must
     * ensure that othr routines do not try to free it too 
     */
  }
  destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );

  return 0;
}

#if 0
/*
 * Uses fixed weights instead of affine weights.
 * FIXME: banded mode not updated in this code yet
 */
int fixed_malign(MOVERLAP *moverlap, ALIGN_PARAMS *params) {

    char *seq1, *seq2;
    int seq1_len, seq2_len, seq_out_len;
    int edge_inc;
    int i,j;
    int s;
    char *seq1_out, *seq2_out;
    int b_c, b_r;

    int t,big_neg,b_s,e,b_e;
    int *F1, *F2, *G1, *G2, *H1, *H2;
    int *pF1, *pF2, *pG1, *pG2, *pH1, *pH2;
    int *t_pF2, *t_pG2, *t_pH2;
    int best_F1, best_G1, best_H1, V_diag, V_extx, V_insx, V_insy;
    int E_gap, F_gap;
    int edge_mode, best_edge_score;

    int band, band_length, two_band, band_left, band_right, first_band_left=0;
    int off_set, guard_offset, *pF_guard;
    int row, first_row, max_row, column, max_col;
    unsigned char *bit_trace;
    int byte, nibble, size_mat;
    char OLD_PAD_SYM, NEW_PAD_SYM;
    int gap_open_x, gap_extend_x;
    int row_index, gap_open_index, gap_extend_index, gap_match_index;
    int gap_open_score, gap_extend_score, gap_pen;
    int **scores;
    MALIGN *malign;


    /*
     *    Three possible alignment cases:
     *    IGAxi   AIGAHxi   GAxi--
     *    LGVyj   GVyj--    SLGVHyj
     *       F      G            H
     *    i.e. xi aligned with yj, xi aligned opposite a gap y,
     *    or yi aligned opposite a gap in x
     *    below these cases are contained in the recurrence relations
     *    for F, G and H respectively
     *    s(xi,yj) is score matrix
     *    d is gap_open
     *    e is gap extend
     *
     *                   F(i-1,j-1) + s(xi,yi)
     *    F(i,j)  = max  H(i-1,j-1) + s(xi,yi)      \  no gap
     *                   G(i-1,j-1) + s(xi,yi)
     *
     *                   F(i,j-1)   - d
     *    G(i,j) = max   G(i,j-1)   - e             |  gap in y
     *
     *                   F(i-1,j)   - d
     *    H(i,j) = max   H(i-1,j)   - e             -  gap in x
     *                  
     *              
     *    Find MAX(F(i,j),G(i,j),H(i,j)) and set trace accordingly:
     *                \     |      -
     *
     *    if gaps at left edge count:
     *    G(0,i) = G(i,0) = H(0,i) = H(i,0) = - d - e * i
     *    F(1,i) = F(i,1) = - d - e * i;
     *    F(0,0) = 0;
     *    otherwise all set to 0;
     *    if right end gaps count the best score is at (seq1_len,seq2_len)
     *    otherwise find the best score along the two edges
     *    trace back accordingly
     *
     *    store 2 rows for each of F, G, H
     *    use p_F1, p_G1, p_G1 to point to previous row
     *    p_F2, p_G2, p_H2 for current row being built
     *    at the start of a new row:
     *
     *    rows have length seq1_len, columns seq2_len
     *    i.e.
     *    rows: 1 - seq1_len, columns 1 - seq2_len
     *    seq1xxxxxxxxxxxxxxx
     *   s
     *   e
     *   q
     *   2
     *   y
     *   y
     *   y
     *
     *   multiple sequence alignment version.
     * 
     *   compare a multiple sequence alignment (malign) against a single seq (seq)
     *
     *   Want the scoring to be compatible with the standard scoring matrix used
     *   for aligning pairs of sequences and to use affine gap weighting scheme.
     *   The malign can have ragged ends - ie different depths and the scores
     *   should reflect this.
     *
     *   Let gap_open = d, gap_extend = e, depth of malign at i is n(i)
     *   character type k has count C(i,k) at position i
     *   score matrix value for character types j,k = M(j,k)
     *
     *   score for character types j,k at i = n(i).C(i,j).M(j,k) if C(i,j) !=0
     *                                      = n(i).M(j,k)        if C(i,j) =0
     *   ie if the character in seq occurs in malign +ve score
     *      if not mismatch score, in both cases score multiplied by depth.
     *
     *   score for introducing gap in seq   = n(i).C(i,j).M(j,k) if C(i,'*') !=0
     *                                      = n(i).d             if C(i,'*') =0
     *   score for extending gap in seq     = n(i).C(i,j).M(j,k) if C(i,'*') !=0
     *                                      = n(i).e             if C(i,'*') =0
     *   ditto for gaps in malign.
     *   NB i thought that gaps in malign should be weighted to reflect the
     *   depth so that we would favour gaps in malign over gaps in seq, but
     *   this resulted in too many gaps in seq:
     *   eg match = 4, mismatch -8, d -12, e -4, 2 sequences in malign
     *   1 gap  in seq(-12) > 1 mismatch(-16) > 1 gap in malign(-24)
     *   6 gaps in seq(-32) = 2 mismatch(-32) < 1 gap in malign(-24)
     *   so they gaps in seq and malign are scored the same
     *
     *   Implementation
     *   
     *   precompute the malign scores (scores[][]) for every position
     *   add two extra rows to this matrix: one for gap_open, one for gap_extend
     *   these are precomputed. scores[malign_length][charset_size+2]
     *   If edge gaps count, fill the edges with the values from the gap_open
     *   and gap_extend scores. ie correctly reflect the depth of the malign
     *   throughout its length.
     *
     *   gap_open_index and gap_extend_index are indices to the gap_open and
     *   gap_extend elements of the malign score matrix, gap_match_index to
     *   the elements containing the scores for padding characters
     *
     *   In the recurrence, for each row, we check the sequence character type:
     *   if it is not a pad we set
     *              gap_open_x   = gap_open_index
     *              gap_extend_x = gap_extend_index
     *          else{
     *              gap_open_x = gap_extend_x = gap_match_index
     *
     *   this enables us to set the appropriate costs for introducing or aligning
     *   gaps: if the seq contains a pad it gets a scores dependent on the number
     *   of pads in the malign at that point; if it does not contain a pad it gets
     *   the appropriate gap penalties.
     *
     *   For each column:
     *
     *          gap_open_score   = scores[column-1][gap_open_x]
     *          gap_extend_score = scores[column-1][gap_extend_x]
     */

    malign = moverlap->malign;
    scores = malign->scores;
    gap_open_index = malign->charset_size;
    gap_extend_index = gap_open_index + 1;
    gap_match_index = malign->charset_size - 2;

    F1 = F2 = G1 = G2 = H1 = H2 = NULL;
    bit_trace = NULL;
    seq1_out = seq2_out = NULL;
    big_neg = INT_MIN/2;
    best_edge_score = big_neg;

    seq1 = moverlap->malign->consensus;
    seq2 = moverlap->seq2;
    seq1_len = moverlap->malign_len;
    seq2_len = moverlap->seq2_len;

    edge_mode = params->edge_mode;

    OLD_PAD_SYM = params->old_pad_sym;
    NEW_PAD_SYM = params->new_pad_sym;
    band = params->band;
    first_row = params->first_row;
    band_left = params->band_left;
    band_right = params->band_right;


    /* init tables */

    if(!(F1 = (int *) xcalloc(seq1_len + 2, sizeof(int)))) {
	verror(ERR_WARN, "affine_align", "xmalloc failed for F1");
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }
    if(!(F2 = (int *) xcalloc(seq1_len + 2, sizeof(int)))) {
	verror(ERR_WARN, "affine_align", "xmalloc failed for F2");
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }
    if(!(G1 = (int *) xcalloc(seq1_len + 2, sizeof(int)))) {
	verror(ERR_WARN, "affine_align", "xmalloc failed for G1");
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }
    if(!(G2 = (int *) xcalloc(seq1_len + 2, sizeof(int)))) {
	verror(ERR_WARN, "affine_align", "xmalloc failed for G2");
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }
    if(!(H1 = (int *) xcalloc(seq1_len + 2, sizeof(int)))) {
	verror(ERR_WARN, "affine_align", "xmalloc failed for H1");
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }
    if(!(H2 = (int *) xcalloc(seq1_len + 2, sizeof(int)))) {
	verror(ERR_WARN, "affine_align", "xmalloc failed for H2");
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }

    /* do recurrence */

    if ( edge_mode & EDGE_GAPS_COUNT ) {
	F1[0] = 0;
	E_gap = scores[0][gap_open_index];
	H1[0] = G1[0] = E_gap;
	for(i = 1; i < seq1_len; i++) {
	    F1[i] = E_gap;
	    E_gap += scores[i][gap_extend_index];
	    H1[i] = E_gap;
	    G1[i] = E_gap;
	}
	E_gap = scores[0][gap_open_index] + scores[0][gap_extend_index];
	/* F_gap = F1[0]; */
	F_gap = scores[0][gap_open_index];
	edge_inc = 1;
    }
    else if ( edge_mode & EDGE_GAPS_ZERO ) {
	for(i = 0; i <= seq1_len; i++) F1[i] = 0;
	for(i = 0; i <= seq1_len; i++) G1[i] = 0;
	for(i = 0; i <= seq1_len; i++) H1[i] = 0;
	edge_inc = 0;
	E_gap = 0;
	F_gap = 0;

	/* FIXME - so we insert pads in X but not Y */
	edge_inc = 1;
	E_gap = scores[0][gap_open_index] + scores[0][gap_extend_index];
	F_gap = scores[0][gap_open_index];
    }
    else {
	printf("scream: unknown gaps mode\n");
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }
  
    /* process each row. i.e. each character of seq2 */

    b_s = big_neg;
    b_e = b_r = b_c = 0;
    t = 0;

    if ( band ) {
	/*
	 * If band is set we restrict the search to band elements either side
	 * of the main diagonal. This gives a band length of (2 * band) + 1.
	 * For the affine alignment it is necessary to know what happened in
	 * the elements to the left and above the current one, and therefore
	 * need to add another 2 to the band length.
	 */

	band_length = (2 * band) + 3;
	/*
	 * Need to add one to first_row, first_column, band_left and
	 * band_right because the alignment matrix has an extra row
	 * and an extra column over the number given by the lengths
	 * of the two sequences.
	 */
	size_mat = (MIN(seq1_len - band_left, seq2_len - first_row) + 1) 
	    * band_length;

	if(!(bit_trace = (unsigned char *) 
	     xmalloc(1 + sizeof(char) * size_mat / 4))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, 
			     seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat / 4;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	first_row++;
	band_left++;
	band_right++;
	first_band_left = band_left;

	two_band = band * 2;
	max_row = MIN(seq2_len, first_row + seq1_len - band_left + 1);

	for(row = first_row; row <= max_row; row++, band_left++, band_right++) {

	    guard_offset = band_left + two_band;

	    if(t == 0) {
		pF1   = F1;
		pF2   = F2;
		pF_guard = F1 + guard_offset;
		t = 1;
	    } else {
		pF1   = F2;
		pF2   = F1;
		pF_guard = F2 + guard_offset;
		t = 0;
	    }
	    pF2[0] = F_gap;

	    if ( (off_set = band_left - 1 ) > 0 ) {
		pF1 += off_set;
		pF2 += off_set;
		*pF2 = big_neg;
	    }
	    t_pF2 = pF2;

	    if ( band_right <= seq1_len ) {
		*pF_guard = big_neg;
	    }
	    gap_pen = scores[MIN(row-1,seq1_len-1)][gap_extend_index];
	    if(edge_inc) {
		E_gap += gap_pen;
		F_gap += gap_pen;
	    }
	    row_index = malign_lookup[(int)seq2[row-1]];

	    /* FIXME: got the axes switched for the single sequence methods
	     * what is the correct thing to do here????
	     */

	    /* Always use the score matrix for pads */
	    if ( seq2[row-1] != OLD_PAD_SYM && 0) {
		gap_open_x = gap_open_index;
		gap_extend_x = gap_extend_index;
	    }
	    else {
		gap_open_x = gap_extend_x = gap_match_index;
	    }
	    /*
	     * Need to prevent the possibility of certain transitions
	     * at the ends of each band; column-1 should not be allowed at
	     * the left-hand end of a band, and row-1 should not be
	     * allowed at the right-hand end. Such transitions would
	     * take the trace-back out of the the parts of the 
	     * alignment matrix that are covered by the band.
	     *
	     *      0 1 2 3 4 5 6
	     *    0   _ _ _           For example, shouldn't consider
	     *    1  |_|_|_|_         column-1 at position (row = 3, col = 3)
	     *    2    |_|_|_|_       so give (3, 2) a very low score.
	     *    3      |_|_|_|_     Likewise, shouldn't consider row-1
	     *    4        |_|_|_|    at position (3, 5), so give (2, 5)
	     *    5          |_|_|    a very low score.
	     * this is done using the guard pointers at the right edge
	     * and initialising the left edges where required.   
	     */

	    /* process each column. i.e. each character of seq1 */
	  
	    max_col = MIN(seq1_len, band_right);
	    for(column = MAX(1, band_left);
		column <= max_col;
		column++, pF1++, pF2++, pG1++, pG2++, pH1++, pH2++) {

		byte = e / 4;
		nibble = 2 * (e % 4);

		s = scores[column-1][row_index];


		V_diag = *pF1 + s;
		V_insx = *pH1 + s;
		V_insy = *pG1 + s;

		V_diag = pF1[0] + s;
		V_insx = pF2[0] + scores[column-1][gap_match_index];
		V_insy = pF1[-1] + scores[column-1][gap_open_index];

		if ( V_insx > V_diag ) {
		    if ( V_insx > V_insy ) {
			best_F1 = V_insx;
			bit_trace[byte] |= BYTE_ACROSS << nibble;
			b_s = V_insx;
		    }
		    else {
			best_F1 = V_insy;
			bit_trace[byte] |= BYTE_DOWN << nibble;
			b_s = V_insy;
		    }
		}
		else {
		    if ( V_diag > V_insy ) {
			best_F1 = V_diag;
			bit_trace[byte] |= BYTE_DIAGONAL << nibble;
			b_s = V_diag;
		    }
		    else {
			best_F1 = V_insy;
			bit_trace[byte] |= BYTE_DOWN << nibble;
			b_s = V_insy;
		    }
		}
		*(pF2+1) = best_F1;

		/* gap in x? */
		V_diag =  *pF2 + gap_open_score;
		V_extx =  *pH2 + gap_extend_score;
		if ( V_diag > V_extx ) {
		    best_H1 = V_diag;
		}
		else {
		    best_H1 = V_extx;
		}
		*(pH2+1) = best_H1;
	    }
	  
	    if ( column > seq1_len ) {
		if ( edge_mode & BEST_EDGE_TRACE ) {
		    best_H1 = MAX(best_H1,best_G1);
		    best_F1 = MAX(best_H1,best_F1);
		    if ( best_F1 > best_edge_score ) {
			best_edge_score = best_F1;
			b_r = row;
			b_e = ((row - first_row + 1) * band_length) + 
			    (seq1_len - band_left + 1);
		    }
		}
	    }
	    /*    
		  for(q=0;q<seq1_len;q++)printf(" %3d ",pF2[q]);
		  printf("\n");
		
		  for(q=0;q<seq1_len;q++)printf(" %3d ",pG2[q]);
		  printf("\n");
		  for(q=0;q<seq1_len;q++)printf(" %3d ",pH2[q]);
		  printf("\n");
	    */
	  
	}
	
	
        row--;
	band_left--;
	band_right--;
	if ( edge_mode & BEST_EDGE_TRACE ) {
	    b_c = seq1_len;

	    pF2 = t_pF2+1;
	    pG2 = t_pG2+1;
	    pH2 = t_pH2+1;
	    max_col = MIN(seq1_len, band_right);
	    for(column = MAX(1, band_left);
		column <= max_col;
		column++, pF2++, pG2++, pH2++) {
		best_F1 = *pF2;
		best_G1 = *pG2;
		best_H1 = *pH2;
		best_H1 = MAX(best_H1,best_G1);
		best_F1 = MAX(best_H1,best_F1);
		if ( best_F1 > best_edge_score ) {
		    best_edge_score = best_F1;
		    b_c = column;
		    b_r = seq2_len;
		    b_e = ((row - first_row + 1) * band_length) + 
			(column - band_left + 1);
		}
	    }
	}
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	    b_r = seq2_len;
	    b_c = seq1_len;
	    b_e = ((b_r - first_row + 1) * band_length) + 
		(b_c - band_left + 1);
	    best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	    b_r = seq2_len;
	    b_c = seq1_len;
	    b_e = ((b_r - first_row + 1) * band_length) + 
		(b_c - band_left + 1);
	    best_edge_score = b_s;
        }

    } else {
	/* Initialise the bit trace */
	size_mat = (seq1_len + 1) * (seq2_len + 1);
	if(!(bit_trace = (unsigned char *)xmalloc(1 +  size_mat / 4))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2,
			     bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat / 4;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}

	/* Step through rows, where a row is a consensus vector */
	for(row = 1, e = seq1_len + 2; row <= seq2_len; row++, e++) {
	    if(t == 0) {
		pF1 = F1;
		pF2 = F2;
		t = 1;
	    } else {
		pF1 = F2;
		pF2 = F1;
		t = 0;
	    }

	    pF2[0] = F_gap;
	  
	    if (0) {
		int i;
		printf("F Row %d:", row-1);
		for (i = 0; i <= seq1_len; i++) {
		    printf(" %+4d", pF1[i]);
		}
		puts("");
	    }

	    gap_pen = scores[MIN(row-1,seq1_len-1)][gap_extend_index];
	    if(edge_inc) {
		F_gap += gap_pen;
	    }

	    row_index = malign_lookup[(int)seq2[row-1]];
	  
	    /* process each column. i.e. each character of seq1 */
	    for(column = 1; column <= seq1_len; column++, e++) {
		byte = e / 4;
		nibble = 2 * (e % 4);

		s = scores[column-1][row_index];

		V_diag = pF1[column-1] + s;
		V_insx = pF2[column-1] + scores[column-1][gap_match_index];
		V_insy = pF1[column] + scores[column-1][gap_open_index];

		if ( V_insx > V_diag ) {
		    if ( V_insx > V_insy ) {
			best_F1 = V_insx;
			bit_trace[byte] |= BYTE_ACROSS << nibble;
		    } else {
			best_F1 = V_insy;
			bit_trace[byte] |= BYTE_DOWN << nibble;
		    }
		} else {
		    if ( V_diag > V_insy ) {
			best_F1 = V_diag;
			bit_trace[byte] |= BYTE_DIAGONAL << nibble;
		    } else {
			best_F1 = V_insy;
			bit_trace[byte] |= BYTE_DOWN << nibble;
		    }
		}
		pF2[column] = best_F1;
	    }
	  
	    if ( edge_mode & BEST_EDGE_TRACE ) {
		/* best right edge score */
		if ( best_F1 > best_edge_score ) {
		    best_edge_score = best_F1;
		    b_r = row;
		    b_e = (row + 1) * (seq1_len + 1) - 1;
		}
	    }
	}
	if (0) {
	    int i;
	    printf("F Row %d:", row-1);
	    for (i = 0; i <= seq1_len; i++) {
		printf(" %+4d", pF2[i]);
	    }
	    puts("");
	}
	
	if ( edge_mode & BEST_EDGE_TRACE ) {
	    /* FIXME: only find exit point along bottom and not right edge */
	    best_edge_score = big_neg;

	    /* best bottom edge score */
	    b_c = seq1_len;
	    for(column = 1; column <= seq1_len; column++) {
		best_F1 = pF2[column];
		if ( best_F1 > best_edge_score ) {
		    best_edge_score = best_F1;
		    b_c = column;
		    b_r = seq2_len;
		    b_e = (row - 1) * (seq1_len + 1) + column;
		}
	    }
	}
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	    b_r = seq2_len;
	    b_c = seq1_len;
	    b_e = seq2_len * (seq1_len + 1) + seq1_len;
	    best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	    b_r = seq2_len;
	    b_c = seq1_len;
	    b_e = seq2_len * (seq1_len + 1) + seq1_len;
	    best_edge_score = b_s;
	}
    }


    /* do traceback */

    moverlap->score = best_edge_score;

    if( i = do_trace_back_bits ( bit_trace, seq1, seq2, seq1_len, seq2_len,
				 &seq1_out, &seq2_out, &seq_out_len, b_r, b_c, b_e,
				 band, first_band_left, first_row, band_length, NEW_PAD_SYM)) {
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }
    /*
      printf("%s\n",seq1_out);
      printf("%s\n",seq2_out);
    */
    moverlap->seq1_out = seq1_out;
    moverlap->seq2_out = seq2_out;
    moverlap->seq_out_len = seq_out_len;

    if ( i = seq_to_moverlap (moverlap, OLD_PAD_SYM, NEW_PAD_SYM)) {
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
	return -1;
    }

    if ( params->job & RETURN_EDIT_BUFFERS ) {
	if (seq_to_edit ( seq1_out,seq_out_len,&moverlap->S1,&moverlap->s1_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
	    return -1;
	}
	if (seq_to_edit ( seq2_out,seq_out_len,&moverlap->S2,&moverlap->s2_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, NULL, NULL );
	    return -1;
	}
    }

    if ( params->job & RETURN_SEQ ) {
	if ( !(params->job & RETURN_NEW_PADS) ) {
	    old_pads_for_new(seq1_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
	    old_pads_for_new(seq2_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
	}
	seq1_out = seq2_out = NULL; /* stop them being freed! */
    }
    else {
	moverlap->seq1_out = moverlap->seq2_out = NULL;
	/* ie we let destroy_af_mem free the memory, but we must
	 * ensure that othr routines do not try to free it too 
	 */
    }
    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );

    return 0;
}
#endif

/*
 * Realigns an MALIGN using ReAligner (Anson & Myers) scoring functions.
 * This basically scores 0 for a match or a number <= 1 for a mismatch or
 * pad. Pad costs are looked up in the score matrix just like a mismatch is
 * as costs are position specific.
 *
 * This code also scores zero for introducing a pad in the consensus vector
 * such that it will align against an existing pad in the sequence. The reason
 * for this is so that we can leave pads in the sequence (ie do not depad it)
 * and if sequence was already aligned against the vector then the band can
 * be centred perfectly on the main diagonal of the matrix. This allows us to
 * get away with very small bands in order to realign sequences with the
 * minimal amount of effort.
 */
int realigner_malign(MOVERLAP *moverlap, ALIGN_PARAMS *params) {
    char *seq1, *seq2;
    int seq1_len, seq2_len, seq_out_len;
    int edge_inc;
    int i,j;
    char *seq1_out, *seq2_out;
    int b_c, b_r;

    int t,big_neg,b_s,e,b_e;
    int *F1, *F2;
    int *pF1, *pF2;
    int *t_pF2;
    int E_gap, edge_mode, best_edge_score;

    int band, band_length, two_band, band_left, band_right, first_band_left=0;
    int row, first_row, max_row, column, max_col;
    unsigned char *bit_trace;
    int byte, nibble, size_mat;
    char OLD_PAD_SYM, NEW_PAD_SYM;
    int row_index, gap_match_index;
    int **scores;
    MALIGN *malign;


    malign = moverlap->malign;
    scores = malign->scores;
    gap_match_index = malign_lookup['*'];

    F1 = F2 = NULL;
    bit_trace = NULL;
    seq1_out = seq2_out = NULL;
    big_neg = INT_MAX/2;
    best_edge_score = big_neg;

    seq1 = moverlap->malign->consensus;
    seq2 = moverlap->seq2;
    seq1_len = moverlap->malign_len;
    seq2_len = moverlap->seq2_len;

    edge_mode = params->edge_mode;

    OLD_PAD_SYM = params->old_pad_sym;
    NEW_PAD_SYM = params->new_pad_sym;
    band = params->band;
    first_row = params->first_row;
    band_left = params->band_left;
    band_right = params->band_right;


    /* init tables */
    if(!(F1 = (int *) xcalloc(seq1_len + 2, sizeof(int)))) {
	verror(ERR_WARN, "affine_align", "xmalloc failed for F1");
	destroy_af_mem ( F1, F2, 0, 0, 0, 0, bit_trace, seq1_out, seq2_out );
	return -1;
    }
    if(!(F2 = (int *) xcalloc(seq1_len + 2, sizeof(int)))) {
	verror(ERR_WARN, "affine_align", "xmalloc failed for F2");
	destroy_af_mem ( F1, F2, 0, 0, 0, 0, bit_trace, seq1_out, seq2_out );
	return -1;
    }

    /* do recurrence */
    if ( edge_mode & EDGE_GAPS_COUNT ) {
	for(E_gap = 0, i = 0; i < seq1_len; i++) {
	    F1[i] = E_gap;
	    E_gap += 100;
	}
	E_gap = 100;
	edge_inc = 1;
    } else if ( edge_mode & EDGE_GAPS_ZEROX ) {
	/* ZEROX means we penalise in Y but not X */
	for(i = 0; i <= seq1_len; i++)
	    F1[i] = 0;
	edge_inc = 1;
	E_gap = 100;
    } else {
	printf("scream: unknown gaps mode\n");
	destroy_af_mem ( F1, F2, 0, 0, 0, 0, bit_trace, seq1_out, seq2_out );
	return -1;
    }
  
    /* process each row. i.e. each character of seq2 */
    b_s = big_neg;
    b_e = b_r = b_c = 0;
    t = 0;

    if (band) {
	int off_set, guard_offset, *pF_guard;

	/*
	 * If band is set we restrict the search to band elements either side
	 * of the main diagonal. This gives a band length of (2 * band) + 1.
	 * It is necessary to know what happened in the elements to the left
	 * and above the current one, and therefore need to add another 2 to
	 * the band length.
	 */
	band_length = (2 * band) + 3;

	/*
	 * Need to add one to first_row, first_column, band_left and
	 * band_right because the alignment matrix has an extra row
	 * and an extra column over the number given by the lengths
	 * of the two sequences.
	 */
	size_mat = (MIN(seq1_len - band_left, seq2_len - first_row) + 1) 
	    * band_length;

	if(!(bit_trace = (unsigned char *) 
	     xmalloc(1 + sizeof(char) * size_mat / 4))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, 0, 0, 0, 0, bit_trace, 
			     seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat / 4;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}
	first_row++;
	band_left++;
	band_right++;
	first_band_left = band_left;

	two_band = band * 2;
	max_row = MIN(seq2_len, first_row + seq1_len - band_left + 1);

	for(row = first_row; row <= max_row; row++, band_left++, band_right++) {
	    int e_row, e_col, insy_cost;

	    guard_offset = band_left + two_band;

	    if(t == 0) {
		pF1 = F1;
		pF2 = F2;
	    } else {
		pF1 = F2;
		pF2 = F1;
	    }
	    t ^= 1;
	    pF_guard = pF1 + guard_offset;

	    /*
	     * Need to prevent the possibility of certain transitions
	     * at the ends of each band; column-1 should not be allowed at
	     * the left-hand end of a band, and row-1 should not be
	     * allowed at the right-hand end. Such transitions would
	     * take the trace-back out of the the parts of the 
	     * alignment matrix that are covered by the band.
	     *
	     *      0 1 2 3 4 5 6
	     *    0   _ _ _           For example, shouldn't consider
	     *    1  |_|_|_|_         column-1 at position (row = 3, col = 3)
	     *    2    |_|_|_|_       so give (3, 2) a very bad score.
	     *    3      |_|_|_|_     Likewise, shouldn't consider row-1
	     *    4        |_|_|_|    at position (3, 5), so give (2, 5)
	     *    5          |_|_|    a very bad score.
	     * this is done using the guard pointers at the right edge
	     * and initialising the left edges where required.   
	     */
	    pF2[0] = E_gap;
	    if ( (off_set = band_left - 1 ) >= 0 ) {
		pF1 += off_set;
		pF2 += off_set;
		*pF2 = big_neg;
	    }
	    t_pF2 = pF2;

	    if ( band_right <= seq1_len ) {
		*pF_guard = big_neg;
	    }

	    if(edge_inc) {
		E_gap += scores[MIN(row-1,seq1_len-1)][gap_match_index];
	    }
	    row_index = malign_lookup[(int)seq2[row-1]];

	    /* process each column. i.e. each character of seq1 */
	    max_col = MIN(seq1_len, band_right);
	    column = MAX(1, band_left);
	    e_row = (row - first_row + 1) * band_length;
	    e_col = column - band_left + 1;
	    e = e_row + e_col;
	    insy_cost = (seq2[row-1] == '*' ? 0 : 100);

	    for(; column <= max_col; e++, column++, pF1++, pF2++) {
		int V_diag, V_insx, V_insy;

		byte = e / 4;
		nibble = 2 * (e % 4);

		V_diag = pF1[0] + scores[column-1][row_index];
		V_insx = pF2[0] + scores[column-1][gap_match_index];
		V_insy = pF1[1] + insy_cost;

		if (V_diag <= V_insx && V_diag <= V_insy) {
		    b_s = V_diag;
		    bit_trace[byte] |= BYTE_DIAGONAL << nibble;
		} else if (V_insx <= V_insy) {
		    b_s = V_insx;
		    bit_trace[byte] |= BYTE_ACROSS << nibble;
		} else {
		    b_s = V_insy;
		    bit_trace[byte] |= BYTE_DOWN << nibble;
		}

		*(pF2+1) = b_s;
	    }
	}
	
	
        row--;
	band_left--;
	band_right--;
	if ( edge_mode & BEST_EDGE_TRACE ) {
	    /* Only find exit point along bottom and not right edge */
	    best_edge_score = big_neg;

	    /* best bottom edge score */
	    b_c = seq1_len;
	    pF2 = t_pF2+1;
	    max_col = MIN(seq1_len, band_right);
	    for(column = MAX(1, band_left);
		column <= max_col;
		column++, pF2++) {
		if (best_edge_score > *pF2) {
		    best_edge_score = *pF2;
		    b_c = column;
		    b_r = seq2_len;
		    b_e = ((row - first_row + 1) * band_length) + 
			(column - band_left + 1);
		}
	    }
	}

	if ( edge_mode & FULL_LENGTH_TRACE ) {
	    b_r = seq2_len;
	    b_c = seq1_len;
	    b_e = ((b_r - first_row + 1) * band_length) + 
		(b_c - band_left + 1);
	    best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	    b_r = seq2_len;
	    b_c = seq1_len;
	    b_e = ((b_r - first_row + 1) * band_length) + 
		(b_c - band_left + 1);
	    best_edge_score = b_s;
        }
    } else {
	/* NON-BANDED code */

	/* Initialise the bit trace */
	size_mat = (seq1_len + 1) * (seq2_len + 1);
	if(!(bit_trace = (unsigned char *)xmalloc(1 +  size_mat / 4))) {
	    verror(ERR_WARN, "affine_align", "xmalloc failed for bit_trace");
	    destroy_af_mem ( F1, F2, 0, 0, 0, 0,
			     bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	j = 1 + size_mat / 4;
	for(i = 0; i < j; i++) {
	    bit_trace[i] = 0;
	}

	/* Step through rows, where a row is a consensus vector */
	for(row = 1, e = seq1_len + 2; row <= seq2_len; row++, e++) {
	    int insy_cost;

	    if(t == 0) {
		pF1 = F1;
		pF2 = F2;
	    } else {
		pF1 = F2;
		pF2 = F1;
	    }
	    t ^= 1;

	    pF2[0] = E_gap;
	  
	    if(edge_inc) {
		E_gap += scores[MIN(row-1,seq1_len-1)][gap_match_index];
	    }
	    row_index = malign_lookup[(int)seq2[row-1]];

	    insy_cost = (seq2[row-1] == '*' ? 0 : 100);
	  
	    /* process each column. i.e. each character of seq1 */
	    for(column = 1; column <= seq1_len; column++, e++) {
		int V_diag, V_insx, V_insy;

		byte = e / 4;
		nibble = 2 * (e % 4);


		V_diag = pF1[column-1] + scores[column-1][row_index];
		V_insx = pF2[column-1] + scores[column-1][gap_match_index];
		V_insy = pF1[column] + insy_cost;

		if (V_diag <= V_insx && V_diag <= V_insy) {
		    b_s = V_diag;
		    bit_trace[byte] |= BYTE_DIAGONAL << nibble;
		} else if (V_insx <= V_insy) {
		    b_s = V_insx;
		    bit_trace[byte] |= BYTE_ACROSS << nibble;
		} else {
		    b_s = V_insy;
		    bit_trace[byte] |= BYTE_DOWN << nibble;
		}

		pF2[column] = b_s;
	    }
	}
	
	if ( edge_mode & BEST_EDGE_TRACE ) {
	    /* Only find exit point along bottom and not right edge */
	    best_edge_score = big_neg;

	    /* best bottom edge score */
	    b_c = seq1_len;
	    for(column = 1; column <= seq1_len; column++) {
		if (best_edge_score > pF2[column]) {
		    best_edge_score = pF2[column];
		    b_c = column;
		    b_r = seq2_len;
		    b_e = (row - 1) * (seq1_len + 1) + column;
		}
	    }
	}
	if ( edge_mode & FULL_LENGTH_TRACE ) {
	    b_r = seq2_len;
	    b_c = seq1_len;
	    b_e = seq2_len * (seq1_len + 1) + seq1_len;
	    best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	    b_r = seq2_len;
	    b_c = seq1_len;
	    b_e = seq2_len * (seq1_len + 1) + seq1_len;
	    best_edge_score = b_s;
	}
    }


    /* do traceback */

    moverlap->score = best_edge_score;

    if( i = do_trace_back_bits ( bit_trace, seq1, seq2, seq1_len, seq2_len,
				 &seq1_out, &seq2_out, &seq_out_len, b_r, b_c, b_e,
				 band, first_band_left, first_row, band_length, NEW_PAD_SYM)) {
	destroy_af_mem ( F1, F2, 0,0,0,0, bit_trace, seq1_out, seq2_out );
	return -1;
    }
    /*
      printf("%s\n",seq1_out);
      printf("%s\n",seq2_out);
    */
    moverlap->seq1_out = seq1_out;
    moverlap->seq2_out = seq2_out;
    moverlap->seq_out_len = seq_out_len;

    if ( i = seq_to_moverlap (moverlap, OLD_PAD_SYM, NEW_PAD_SYM)) {
	destroy_af_mem ( F1, F2, 0,0,0,0, bit_trace, NULL, NULL);
	return -1;
    }

    if ( params->job & RETURN_EDIT_BUFFERS ) {
	if (seq_to_edit ( seq1_out,seq_out_len,&moverlap->S1,&moverlap->s1_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, 0,0,0,0, bit_trace, NULL, NULL );
	    return -1;
	}
	if (seq_to_edit ( seq2_out,seq_out_len,&moverlap->S2,&moverlap->s2_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, 0,0,0,0, bit_trace, NULL, NULL );
	    return -1;
	}
    }

    if ( params->job & RETURN_SEQ ) {
	if ( !(params->job & RETURN_NEW_PADS) ) {
	    old_pads_for_new(seq1_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
	    old_pads_for_new(seq2_out,seq_out_len,OLD_PAD_SYM,NEW_PAD_SYM);
	}
	seq1_out = seq2_out = NULL; /* stop them being freed! */
    } else {
	moverlap->seq1_out = moverlap->seq2_out = NULL;
	/* ie we let destroy_af_mem free the memory, but we must
	 * ensure that othr routines do not try to free it too 
	 */
    }
    destroy_af_mem ( F1, F2, 0,0,0,0, bit_trace, seq1_out, seq2_out );

    return 0;
}
