#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <staden.h>
#include <sp_alignment_structs.h>
#include <sp_align_lib.h>
#include <sp_hash_lib.h>
#include <sp_alignment.h>
#include <read_matrix.h>


namespace sp {


/* so far this is all for DNA so lets initialise the tables here */

/**
  * initialise dna lookup and dna character set
  */
void init_DNA_lookup() {
    set_dna_lookup();
    set_char_set(1);
}

/**
  * create overlap structure and set pointers to NULL
  */

OVERLAP *create_overlap(void) {
    OVERLAP *overlap;

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

/**
  * initialise overlap structure:
  * set pointers to seqs and set their lengths
  */

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
    overlap->percent = 0.0;
    overlap->left1 = 0;
    overlap->left2 = 0;
    overlap->right1 = 0;
    overlap->right2 = 0;
    overlap->left = 0;
    overlap->right = 0;
    overlap->length = 0;
    overlap->direction = 0;
    overlap->lo = 0;
    overlap->ro = 0;
    overlap->score = 0.0;
    overlap->qual = 0.0;
}


/**
  * destroy overlap structure and all its contents 
  */

void destroy_overlap (OVERLAP *overlap) {
   if ( overlap ) {
    if ( overlap->S1 ) xfree ( overlap->S1 );
    if ( overlap->S2 ) xfree ( overlap->S2 );
    if ( overlap->S ) xfree ( overlap->S );
    if ( overlap->seq1_out ) xfree ( overlap->seq1_out );
    if ( overlap->seq2_out ) xfree ( overlap->seq2_out );
    xfree ( overlap );
   }
   overlap = NULL;
}

/**
  * free overlap structure output arrays and set to NULL
  * ie S1, S2, S, seq1_out, seq2_out
  */

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

/**
  * add seq1 to overlap
  *
  */

int set_overlap_seq1(OVERLAP *overlap, char *seq, int seq_len) {

    /* paranoia */

    if (NULL == overlap) return -1;

    overlap->seq1 = seq;
    overlap->seq1_len = seq_len;
    return 0;
}

/**
  * add seq2 to overlap
  *
  */

int set_overlap_seq2(OVERLAP *overlap, char *seq, int seq_len) {

    /* paranoia */

    if (NULL == overlap) return -1;

    overlap->seq2 = seq;
    overlap->seq2_len = seq_len;
    return 0;
}

/**
  * destroy alignment parameters and all its contents
  * ie all the hashing arrays and the align structure itself
  */

void destroy_align_params (ALIGN_PARAMS *params) {
    if (params) {
   destroy_hash8n(params->hash);
   xfree (params);
    }
}

/**
  *create align params structure and initialise what we can
  * ie gap_open, gap_extend, band, edge mode, return job and word length
  */

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
    params->edge_mode = SP_ALIGNMENT_LEFT_EDGE_GAPS_COUNT | SP_ALIGNMENT_BEST_RIGHT_EDGE;
    params->return_job = SP_ALIGNMENT_RETURN_SEQ;
    params->seq1_start = 0;
    params->seq2_start = 0;
    params->seq1_end = 0;
    params->seq2_end = 0;
    params->new_pad_sym = '.';
    params->old_pad_sym = '*';
    params->algorithm = 0;
    params->score_matrix = NULL;
    params->hash = NULL;
    params->word_length = 8;
    params->min_match = 0;
    params->max_prob = 0.0;
    
    return params;
}

/**
  * print align_params structure
  */

void print_align_params(ALIGN_PARAMS *params) {
    printf("gap_open %d\ngap_extend %d\nband %d\nfirst_row %d\n"
      "band_left %d\nband_right %d\nedge_mode %d\njob %d\n"
      "seq1_start %d\nseq2_start %d\nseq1_end %d\nseq2_end %d\n"
      "new_pad_sym %c\nold_pad_sym %c\n"
      "algorithm %d\nword_length %d\n"
      "min_match %d\nmax_prob %e\n",
    params->gap_open,
    params->gap_extend,
    params->band,
    params->first_row,
    params->band_left,
    params->band_right,
    params->edge_mode,
    params->return_job,
    params->seq1_start,
    params->seq2_start,
    params->seq1_end,
    params->seq2_end,
    params->new_pad_sym,
    params->old_pad_sym,
    params->algorithm,
    params->word_length,
    params->min_match,
    params->max_prob);

}
      
/**
  * set up all the alignment parameters (except hashing arrays) if non-zero
  */
#ifdef DYNMAT
int set_align_params (ALIGN_PARAMS *params, int band, int gap_open, 
             int gap_extend, int return_job, 
             int seq1_start, int seq2_start, 
             char old_pad_sym, char new_pad_sym,
             int seq1_end, int seq2_end,
             int algorithm, int word_length, 
             int min_match, int user_edge_mode,
             double max_prob, int** score_matrix )
#else
int set_align_params (ALIGN_PARAMS *params, int band, int gap_open, 
             int gap_extend, int return_job, 
             int seq1_start, int seq2_start, 
             char old_pad_sym, char new_pad_sym,
             int seq1_end, int seq2_end,
             int algorithm, int word_length, 
             int min_match, int user_edge_mode,
             double max_prob, int (*score_matrix)[128][128]) 
#endif
{


  int d, first_column;

  if(seq1_start > 0) params->seq1_start = seq1_start;
  if(seq2_start > 0) params->seq2_start = seq2_start;
  if(seq1_end > 0) params->seq1_end = seq1_end;
  if(seq2_end > 0) params->seq2_end = seq2_end;

  if ( return_job & SP_ALIGNMENT_RETURN_EDIT_BUFFER ) {
    verror(ERR_WARN, "affine_align", "unimplemented alignment job");
    return -1;
  }

  if ( return_job && !(return_job & SP_ALIGNMENT_RETURN_SEQ) && !(return_job & SP_ALIGNMENT_RETURN_EDIT_BUFFERS)) {
    verror(ERR_WARN, "affine_align", "unknown alignment job");
    return -1;
  }
  if(gap_open) params->gap_open = gap_open;
  if(gap_extend) params->gap_extend = gap_extend;

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

  if(return_job) params->return_job = return_job;
  if (old_pad_sym) params->old_pad_sym = old_pad_sym;
  if (new_pad_sym) params->new_pad_sym = new_pad_sym;
  to_internal_edges(user_edge_mode, &params->edge_mode);
  if(score_matrix) params->score_matrix = score_matrix;
  if(algorithm) params->algorithm = algorithm;
  if(word_length) params->word_length = word_length;
  if(min_match) params->min_match = min_match;
  if(max_prob > 0.0) params->max_prob = max_prob;
      
  return 0;
}
   
/**
  * set the align params band value only
  */

int set_align_params_band_size (ALIGN_PARAMS *params, int band) {
    params->band = band;
    return 0;
}

    /**
     * set 0 end points to full length and do range checking
     * assumes overlap already set correctly
     */

int set_align_params_range (ALIGN_PARAMS *params, OVERLAP *overlap,
             int seq1_start, int seq1_end,
             int seq2_start, int seq2_end) {

    if(seq1_start <= 0) seq1_start = 0;
    if(seq2_start <= 0) seq2_start = 0;
    if(seq1_end <= 0) seq1_end = overlap->seq1_len - 1;
    if(seq2_end <= 0) seq2_end = overlap->seq2_len - 1;
    if(seq1_end > overlap->seq1_len - 1) seq1_end = overlap->seq1_len - 1;
    if(seq2_end > overlap->seq2_len - 1) seq2_end = overlap->seq2_len - 1;
    
    params->seq1_start = seq1_start;
    params->seq1_end = seq1_end;
    params->seq2_start = seq2_start;
    params->seq2_end = seq2_end;
    return 0;
}


    /** set parameters related to band:
     * need band and start positions,
     * set band, first_row, band_left, band_right
     */

int set_align_params_banding (ALIGN_PARAMS *params, int band, int seq1_start,
               int seq2_start) {

  int d, first_column;

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

/**
  * convert user edge mode to align params edge mode
  */

int set_align_params_edge_mode (ALIGN_PARAMS *params, int user_edge_mode) {
    to_internal_edges(user_edge_mode, &params->edge_mode);
    return 0;
}

/**
  * set align params algorithm only
  */

int set_align_params_algorithm (ALIGN_PARAMS *params, int algorithm) {
    params->algorithm = algorithm;
    return 0;
}

/**
  * get alignment score matrix from file
  */

#ifdef DYNMAT
int get_alignment_matrix( int** matrix_128, char *fn, char *base_order )
#else
int get_alignment_matrix ( int matrix_128[128][128], char *fn, char *base_order )
#endif
{
    int **input_matrix;
    int unknown;
    int i, j, len;
    
    if ( input_matrix = create_matrix(fn, base_order)) {
   len = strlen(base_order);
   unknown = 1000;
   for (j = 0; j < len; j++) {
       for (i = 0; i < len; i++) {
      unknown = MIN(unknown,input_matrix[i][j]);
       }
   }

   to_128(matrix_128,input_matrix,base_order,unknown);
   free_matrix(input_matrix,base_order);
   return 0;
    }
    verror(ERR_WARN, "get_alignment_matrix", "matrix file not found");
    free_matrix(input_matrix,base_order);
    return -1;
}


/**
  * print a 128 by 128 score matrix
  */

void print_128(int matrix_128[128][128]) {

    int i,j,len = 128;
    putchar('\n');
    for (j = 0; j < len; j++) {
   for (i = 0; i < len; i++) {
       printf("%3d ", matrix_128[i][j]);
   }
   putchar('\n');
    }
}


/**
  * print an alignment from an overlap structure
  * either from edit buffers or from its output arrays
  */
int print_overlap(OVERLAP *overlap, FILE *fpt) {
    
    char *seq1_align, *seq2_align;
    int     seq1_align_len, seq2_align_len;
    char temp_seq[51];
    int     i, j;
    int     max_out_width;
    int     max_seq;
    int     len_align;
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

/**
  * print alignment edit buffers for debugging
  */

void print_edit_buffers(OVERLAP *overlap) {

  int i;

  for(i=0;i<overlap->s1_len;i++) {
    printf("1 %d\n",overlap->S1[i]);
  }

  for(i=0;i<overlap->s2_len;i++) {
    printf("2 %d\n",overlap->S2[i]);
  }
}

/**
  * print the position of an overlap and its percentage match
  */

void print_overlap_posn(OVERLAP *overlap) {

    printf("start %d end %d percent %f\n",overlap->left, overlap->right, overlap->percent);
}

/**
  * print the contents of an overlap structure
  */

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

/**
  * initialise all hashing arrays for all possible algorithms
  * defined by job
  */

int init_hash8n (
      int max_seq,
      int max_diagonal,
      int word_length,
      int max_matches,
      int min_match,
      int job,
      Hash **h) {
/* FIXME why is this in this file?? */    
  int size_hash;
    /* max_seq is the longest sequence, but may be either seq1 or seq2
     * max_diagonal is the length of the second longest sequence
     * if a sequence is going to be compared against itself then
     * max_seq = max_diagonal
     */

    set_hash8_lookupn ();
    
    if ( ! (*h = (Hash *) xmalloc ( sizeof(Hash) ))) return -2;

    if ( (word_length != 8) && (word_length != 4)) {

   if (word_length < 4) word_length = 4;
   if (word_length > 4) word_length = 8;
    }
    size_hash = (int)pow(4.0, (double) word_length);
    
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
    
    if ( ! ((*h)->values1 = (int *) xmalloc ( sizeof(int)*(max_seq) ))) 
   return -2;

    if ( ! ((*h)->values2 = (int *) xmalloc ( sizeof(int)*(max_diagonal) ))) 
   return -2;

    /* at present 3 modes are used:
     * 1) either we do a quick comparison and quick alignment using blocks
     * 2) or we do a quick comparison and slow alignment
     * 3) or simply search for repeats or fortran sequence assembly
     * for 1 need
     * HASH_JOB_DIAG 1
     * HASH_JOB_BLKS  16
     * ie job = 17;
     * 
     * for 2 need
     * HASH_JOB_DIAG 1
     * HASH_JOB_HIST 2
     * HASH_JOB_EXPD 4
     * HASH_JOB_DMTCH 8
     * ie job = 31
     *
     * for 3 need HASH_JOB_DIAG 1 only
     *
     * for 4 need
     * HASH_JOB_DIAG 1
     * HASH_JOB_DMTCH 8
     * ie job = 9
     */

    if( (job != 1) && (job != 17) && (job != 31) && (job != 9) ) return -2;

  if ( ! ((*h)->counts = (int *) xmalloc ( sizeof(int)*((*h)->size_hash) ))) 
      return -2;
  
  if ( ! ((*h)->last_word = (int *) xmalloc ( sizeof(int)*((*h)->size_hash) ))) 
      return -2;    

    if ( HASH_JOB_DIAG & job ) {

   if ( ! ((*h)->diag = (int *) xmalloc ( sizeof(int)*(max_seq+max_diagonal) ))) 
       return -2;
    }

    if ( HASH_JOB_HIST & job ) {

   if(!((*h)->hist = (int *) xmalloc(sizeof(int) * (max_seq+max_diagonal)))) 
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

/**
  * destroy hash structure
  */

void destroy_hash8n ( Hash *h ) {
  /* this should be used instead of free_hash8n */
  if (h) {
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
}

/**
  * free the contents of the hash arrays and set them to NULL
  */
void free_hash8n1 ( Hash *h ) {
/* FIXME why is this here, what uses it??? */    
  if ( h->values1 ) xfree ( h->values1 );
  if ( h->values2 ) xfree ( h->values2 );
  if ( h->counts ) xfree ( h->counts );
  if ( h->last_word ) xfree ( h->last_word );
  if ( h->diag )         xfree ( h->diag );
  if ( h->hist )         xfree ( h->hist );
  if ( h->expected_scores )    xfree ( h->expected_scores );
  if ( h->diag_match )    xfree ( h->diag_match );
  if ( h->block_match )    xfree ( h->block_match );

  h->values1         = NULL;
  h->values2 = NULL;
  h->counts = NULL;
  h->last_word = NULL;
  h->diag   = NULL;
  h->hist = NULL;
  h->expected_scores = NULL;
  h->diag_match      = NULL;
 h->block_match     = NULL;
}


    /**store the hash values in values: put number of occurrences of
     *each hash value in counts; put the array position of the last 
     *occurrence of each hash value in last_word, and previous
     *occurrences in values[last_word]. i.e. values is used for TWO
     *purposes: the initial indexes and then the positions.
     *Note that words containing unknown characters (like '-') are given
     *hash value -1. So we skip them here, and they are ignored.
     */


void store_hashn ( Hash *h ) {

    int nw;
    register int i,j,n;

    for ( i=0;i<h->size_hash;i++ ) {
   h->counts[i] = 0;
   h->last_word[i] = 0;
    }
    j = h->seq1_len - h->word_length + 1;
    for ( i = 0; i < j; i++ ) {
   n = h->values1[i];
   if ( -1 != n ) {
       nw = h->counts[n];
       if ( 0 == nw ) {
      h->last_word[n] = i;
      h->counts[n] += 1;
       }
       else {
      h->counts[n] += 1;
      h->values1[i] = h->last_word[n];
      h->last_word[n] = i;
       }
   }
    }
}

/**
  * hash sequence defined by job using preset parameters
  */

int hash_seqn (Hash *h, int job) {

    if ( job == 1 ) {
   if (h->word_length == 8 ) {
       if ( hash_seq8n ( h->seq1, h->values1, 
             h->seq1_len, h->word_length ) != 0 ) {
      return -1;
       }
   }
   else {
       if ( hash_seq4n ( h->seq1, h->values1, 
             h->seq1_len, h->word_length ) != 0 ) {
      return -1;
       }
   }
   return 0;
    }
    else if ( job == 2 ) {
   if (h->word_length == 8 ) {
       if ( hash_seq8n ( h->seq2, h->values2, 
             h->seq2_len, h->word_length ) != 0 ) {
      return -1;
       }
   }
   else {
       if ( hash_seq4n ( h->seq2, h->values2, 
             h->seq2_len, h->word_length ) != 0 ) {
      return -1;
       }
   }
   return 0;
    }
    return -2;
}

/**
  * return the composition of seq in comp
  * sum(comp) = 1.0
  */

void p_comp(double comp[], char *seq, int seq_len) {
    int i;
    double t;

    for(i=0;i<5;i++)comp[i]=0.0;
    if(seq_len<1)return;
    for ( i = 0; i < seq_len; i++ ) {
        comp[ char_lookup [ (unsigned int)seq[i] ] ] += 1.0;
    }
    for(i=0,t=0.0;i<4;i++)t+=comp[i];
    if(t>0.0) {
   for(i=0;i<4;i++)comp[i]/=t;
    }
}

/**
  * calc expected scores for probability max_prob_in
  * for word_length for all diagonal lengths from min_diag
  * to max_diag and store them in expected_scores
  * probabities of each word depend on the supplied composition
  */

int poisson_diagonals(int min_diag, int max_diag, int word_length,
            double max_prob_in, int *expected_scores,
            double comp[]) {
    int     diagonal_length, hits;
    double  expected_hits, sum_probs;
    double  prob_remaining;
    double  p_w;
    double limit;
    int not_found;
    double max_prob, frac, big;

    double z, emz, x;

    /* FIXME should this be in this file? */

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
/* printf("limit %e x %e sum %e rem %e hits %d\n",limit,x,sum_probs,prob_remaining,hits);*/

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
       printf("not found %d %d\n",diagonal_length,hits);
       expected_scores[diagonal_length] = hits;
   }
    }


    if(max_prob_in<max_prob) {
   frac = 1.0+0.033*log10(max_prob/max_prob_in);
   for(hits = 0; hits < max_diag; hits++) {
       expected_scores[hits] = (int)(expected_scores[hits] * frac);
   }
    }

    /*
       for(hits = 0; hits < max_diag; hits++) {
   printf("hits %d exp %e\n",hits,expected_scores[hits]);
   }
   */

    return 0;
}

/**
  * set up all the hasing arrays hidden in params
  * used by programs that call the aligner once FIXME is this true?
  */

int construct_hash_all(ALIGN_PARAMS *params, OVERLAP *overlap) {

    Hash *h;
    int hash_len1, hash_len2, max_seq, max_matches, longest_diagonal;
    double comp[5];

    hash_len1 = params->seq1_end - params->seq1_start + 1;
    hash_len2 = params->seq2_end - params->seq2_start + 1;

    max_seq = MAX(hash_len1, hash_len2);
    max_matches = max_seq;

    longest_diagonal = MAX(hash_len1, hash_len2);
    max_matches = MIN(10000,longest_diagonal);

    if (init_hash8n(max_seq, longest_diagonal,
          params->word_length, max_matches, params->min_match, 
          params->algorithm, &h )) {
   destroy_hash8n(h);
   return -1;
    }

    h->seq1_len = hash_len1;
    h->seq2_len = hash_len2;
    h->seq1 = overlap->seq1 + params->seq1_start;
    h->seq2 = overlap->seq2 + params->seq2_start;

    if (hash_seqn(h, 1)) {
   destroy_hash8n(h);
      return -1;
    }
    
    if (hash_seqn(h, 2)) {
   destroy_hash8n(h);
      return -1;
    }

    store_hashn ( h );

    if (params->algorithm == SP_ALIGNMENT_ALGORITHM_C) {
        p_comp(comp,overlap->seq1,overlap->seq1_len);
   if(poisson_diagonals(params->min_match, longest_diagonal,
              h->word_length, params->max_prob, 
              h->expected_scores, comp)) {
       destroy_hash8n(h);
       return -1;
   }
    }
    params->hash = h;
    return 0;
}

/**
  * stick seq1 into params, hash and store 
  * FIXME could rename this: prepare_seq1_for_aligner
  * and have it return 0 if not a hashing algorithm
  */

int hash_seq1(ALIGN_PARAMS *params, char *seq, int seq_len) {


    /* paranoia */

    if ((params->algorithm != SP_ALIGNMENT_ALGORITHM_C) &&
       (params->algorithm != SP_ALIGNMENT_ALGORITHM_B)) return -1;
    if (NULL == params->hash) return -1;

    params->hash->seq1 = seq;
    params->hash->seq1_len = seq_len;

    if (hash_seqn(params->hash, 1)) {
      return -1;
    }
    
    store_hashn ( params->hash);
    return 0;
}

/**
  * stick seq2 into params and hash
  * FIXME could rename this: prepare_seq2_for_aligner
  * and have it return 0 if not a hashing algorithm
  */

int hash_seq2(ALIGN_PARAMS *params, char *seq, int seq_len) {


    /* paranoia */

    if ((params->algorithm != SP_ALIGNMENT_ALGORITHM_C) &&
       (params->algorithm != SP_ALIGNMENT_ALGORITHM_B)) return -1;
    if (NULL == params->hash) return -1;

    params->hash->seq2 = seq;
    params->hash->seq2_len = seq_len;

    if (hash_seqn(params->hash, 2)) {
      return -1;
    }

    return 0;
}

/**
  * set the poisson values based on this seq and the contents of params
  * ie params must have been set beforehand
  */

int set_align_params_poisson(ALIGN_PARAMS *params, char *seq, int seq_len) {


    double comp[5];

    if (params->algorithm == SP_ALIGNMENT_ALGORITHM_C) {
        p_comp(comp,seq,seq_len);
   if(poisson_diagonals(params->min_match, 
              MAX(params->hash->seq1_len,params->hash->seq2_len),
              params->word_length, params->max_prob, 
              params->hash->expected_scores, comp)) {
       return -1;
   }
   return 0;
    }
    return -1;
}


/**
  * initialise the hash tables with size max_seq1, max_seq2
  */

int init_align_params_hashing(ALIGN_PARAMS *params, int max_seq1, int max_seq2) {

    Hash *h;
    int hash_len1, hash_len2, max_seq, max_matches, longest_diagonal;

    hash_len1 = max_seq1;
    hash_len2 = max_seq2;

    max_seq = MAX(hash_len1, hash_len2);
    max_matches = max_seq;

    longest_diagonal = MAX(hash_len1, hash_len2);
    max_matches = MIN(10000,longest_diagonal);

    if (init_hash8n(max_seq, longest_diagonal,
          params->word_length, max_matches, params->min_match, 
          params->algorithm, &h )) {
   destroy_hash8n(h);
   return -1;
    }
    params->hash = h;
    return 0;
}

/**
  * having set up overlap, set up align params ready for aligner
  * used by programs that do one off alignments
  * FIXME could rename this: prepare_seqs_for_aligner
  * then it would be compatible with prepare_seq(1,2)_for_aligner
  *
  */

int prepare_for_aligner(ALIGN_PARAMS *params, OVERLAP *overlap) {

    Hash *h;
    int hash_len1, hash_len2, max_seq, max_matches, longest_diagonal;
    double comp[5];

    /* paranoia */

    if (NULL == params) return -2;
    if (NULL == overlap) return -2;

    /* nothing to do? */


    if((params->algorithm != SP_ALIGNMENT_ALGORITHM_C) &&
       (params->algorithm != SP_ALIGNMENT_ALGORITHM_B)) return 0;

    hash_len1 = params->seq1_end - params->seq1_start + 1;
    hash_len2 = params->seq2_end - params->seq2_start + 1;

    max_seq = MAX(hash_len1, hash_len2);
    max_matches = max_seq;

    longest_diagonal = MAX(hash_len1, hash_len2);
    max_matches = MIN(10000,longest_diagonal);

    if (init_hash8n(max_seq, longest_diagonal,
          params->word_length, max_matches, params->min_match, 
          params->algorithm, &h )) {
   destroy_hash8n(h);
   return -1;
    }

    h->seq1_len = hash_len1;
    h->seq2_len = hash_len2;
    h->seq1 = overlap->seq1 + params->seq1_start;
    h->seq2 = overlap->seq2 + params->seq2_start;

    if (hash_seqn(h, 1)) {
   destroy_hash8n(h);
      return -1;
    }
    
    if (hash_seqn(h, 2)) {
   destroy_hash8n(h);
      return -1;
    }

    store_hashn ( h );

    if (params->algorithm == SP_ALIGNMENT_ALGORITHM_C) {
        p_comp(comp,overlap->seq1,overlap->seq1_len);
   if(poisson_diagonals(params->min_match, longest_diagonal,
              h->word_length, params->max_prob, 
              h->expected_scores, comp)) {
       destroy_hash8n(h);
       return -1;
   }
   
    }
    params->hash = h;
    return 0;
}


/**
  * interface to alignment by all alignment routines
  * params and overlap need to be fully set up beforehand
  * (including all hash tables)
  */

int aligner(ALIGN_PARAMS *params, OVERLAP *overlap) {

    if ( params->algorithm == SP_ALIGNMENT_ALGORITHM_A) {
   return affine_align(overlap,params);
    }

    if (params->algorithm == SP_ALIGNMENT_ALGORITHM_B) {
   return compare_b ( params->hash, params, overlap );
    }

    if (params->algorithm == SP_ALIGNMENT_ALGORITHM_C) {
   return compare_c ( params->hash, params, overlap );
    }

    return -1;
}


}
