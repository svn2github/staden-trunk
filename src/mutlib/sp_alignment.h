#if !defined(SP_ALIGNMENT_H)
#define SP_ALIGNMENT_H


#include <stdio.h>    /* For FILE structure */


namespace sp {


#define SP_ALIGNMENT_RETURN_SEQ          1
#define SP_ALIGNMENT_RETURN_EDIT_BUFFERS 2
#define SP_ALIGNMENT_RETURN_EDIT_BUFFER  4
#define SP_ALIGNMENT_RETURN_NEW_PADS     8

#define SP_ALIGNMENT_LEFT_EDGE_GAPS_COUNT   1
#define SP_ALIGNMENT_BEST_RIGHT_EDGE        2

#define EDGE_GAPS_COUNT   1
#define EDGE_GAPS_ZERO    2
#define FULL_LENGTH_TRACE 4
#define BEST_EDGE_TRACE   8

#define HASH_JOB_DIAG 1
#define HASH_JOB_HIST 2
#define HASH_JOB_EXPD 4
#define HASH_JOB_DMTCH 8
#define HASH_JOB_BLKS  16

void init_DNA_lookup(void);


#ifdef DYNMAT
int get_alignment_matrix( int** score_matrix, char *fn, char *base_order);
#else
int get_alignment_matrix (int score_matrix[128][128], char *fn, char *base_order);
#endif


void print_128(int matrix_128[128][128]);

OVERLAP *create_overlap(void);

void init_overlap (OVERLAP *overlap, char *seq1, char *seq2, int seq1_len, int seq2_len);
void destroy_overlap (OVERLAP *overlap);

void free_overlap (OVERLAP *overlap);

int set_overlap_seq1(OVERLAP *overlap, char *seq, int seq_len);
int set_overlap_seq2(OVERLAP *overlap, char *seq, int seq_len);

ALIGN_PARAMS *create_align_params(void);

void destroy_align_params (ALIGN_PARAMS *params);

#ifdef DYNMAT
int set_align_params (ALIGN_PARAMS *params, int band, int gap_open, 
             int gap_extend, int return_job, 
             int seq1_start, int seq2_start, 
             char old_pad_sym, char new_pad_sym,
             int seq1_end, int seq2_end,
             int algorithm, int word_length, int min_match, 
             int user_edge_mode,
             double max_prob, int** score_matrix );
#else
int set_align_params (ALIGN_PARAMS *params, int band, int gap_open, 
             int gap_extend, int return_job, 
             int seq1_start, int seq2_start, 
             char old_pad_sym, char new_pad_sym,
             int seq1_end, int seq2_end,
             int algorithm, int word_length, int min_match, 
             int user_edge_mode,
             double max_prob, int (*score_matrix)[128][128]);
#endif

int prepare_for_aligner(ALIGN_PARAMS *params, OVERLAP *overlap);

int aligner(ALIGN_PARAMS *params, OVERLAP *overlap);


int set_align_params_band_size (ALIGN_PARAMS *params, int band);


int set_align_params_banding (ALIGN_PARAMS *params, int band, int seq1_start,
               int seq2_start);

int set_align_params_range (ALIGN_PARAMS *params, OVERLAP *overlap,
             int seq1_start, int seq1_end,
             int seq2_start, int seq2_end);

int init_align_params_hashing(ALIGN_PARAMS *params, int max_seq1, int max_seq2);

int set_align_params_edge_mode (ALIGN_PARAMS *params, int user_edge_mode);

int affine_align(OVERLAP *overlap, ALIGN_PARAMS *params);

int print_overlap(OVERLAP *overlap, FILE *fpt);

void print_edit_buffers(OVERLAP *overlap);

void print_overlap_posn(OVERLAP *overlap);

void print_overlap_struct(OVERLAP *overlap);

void print_align_params(ALIGN_PARAMS *params);

void set_hash8_lookupn(void);

int init_hash8n (int max_seq, int max_diagonal, int word_length, 
       int max_matches,
       int min_match,
       int job,
       Hash **h);

void destroy_hash8n ( Hash *h );

void free_hash8n1 ( Hash *h );

void store_hashn ( Hash *h );

void p_comp(double comp[], char *seq, int seq_len);

int poisson_diagonals(int min_diag, int max_diag, int word_length,
            double max_prob_in, int *expected_scores,
            double comp[]);

int hash_seqn (Hash *h, int job);

int hash_seq1(ALIGN_PARAMS *params, char *seq, int seq_len);
int hash_seq2(ALIGN_PARAMS *params, char *seq, int seq_len);


}


#endif

