#if !defined(SP_HASH_LIB_H)
#define SP_HASH_LIB_H
#include "sp_alignment_structs.h"
#include "sp_alignment.h"
#include "sp_align_lib.h"


namespace sp {


#define MINMAT 20

typedef struct Edit_pair {
    int *S1;
    int *S2;
    int size;
    int next1;
    int next2;
} EDIT_PAIR;

#define MAX_POLY 20

typedef struct poly_ {
    double a[MAX_POLY];
    double b[MAX_POLY];
    double c[MAX_POLY];
    int num_terms;
    int size_step;
    int rows;
    int cols;
} Poly;


void print_h (Hash *h);

int hash_word8n ( char *seq, int *start_base, int seq_len, int word_length, unsigned short *uword);

int hash_seq8n ( char *seq, int *hash_values, int seq_len, int word_length);

int hash_word4n ( char *seq, int *start_base, int seq_len, int word_length, unsigned char *uword);

int hash_seq4n ( char *seq, int *hash_values, int seq_len, int word_length);


int diagonal_length(int seq1_len, int seq2_len, int diagonal_number);

void diagonal_intercepts (int diagonal_number, int seq1_len, int seq2_len,int *seq1_intercept, int *seq2_intercept );

int poly_mult (Poly *poly);

double prob_word ( int word_length, double comp[]);

int best_intercept ( Hash *h, int *seq1_i, int *seq2_i );

void remdup (int *seq1_match, int *seq2_match, int *len_match, int *n_matches);

void make_reverse (int *seq2_match, int *len_match, int n_matches, int seq2_len);
    
int update_edit_pair ( EDIT_PAIR *edit_pair, OVERLAP *overlap);


int block_to_edit_pair ( EDIT_PAIR *edit_pair, int length );

int align_bit ( ALIGN_PARAMS *params, OVERLAP *overlap, EDIT_PAIR *edit_pair);

void destroy_edit_pair (EDIT_PAIR *edit_pair);

EDIT_PAIR *create_edit_pair(int size);

int set_band_blocks(int seq1_len, int seq2_len);

void left_edit_buffer (OVERLAP *overlap, ALIGN_PARAMS *params, int *S1_next, int *S2_next);

void right_edit_buffer (OVERLAP *overlap, ALIGN_PARAMS *params, int *S1_next, int *S2_next);

int align_wrap ( Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap_out);

int central_diagonal ( Hash *h );

int align_blocks ( Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap);

int compare_seqs(Hash *h, int *seq1_match_pos, int *seq2_match_pos, int *match_length);

int compare_d(Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap);

int compare_b(Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap);
    
int compare_c(Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap);

int reps(Hash *h, int *seq1_match_pos, int *seq2_match_pos, int *match_length, char sense);


}


#endif

