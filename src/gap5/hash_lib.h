#if !defined(HASH_H)
#define HASH_H
#include "align_lib.h"
#include "consen.h"
#include "io_lib/hash_table.h"

#define MINMAT 12

#define HASH_JOB_DIAG       1
#define HASH_JOB_HIST       2
#define HASH_JOB_EXPD       4
#define HASH_JOB_DMTCH      8
#define HASH_JOB_BLKS      16
#define HASH_JOB_COUNTLESS 32

typedef struct block_match_ {
    int pos_seq1;
    int pos_seq2;
    int diag;
    int length;
    int best_score;
    int prev_block;
    int next_block;
    int contig1;
} Block_Match;

typedef struct diag_match_ {
    int pos;
    double prob;
} Diag_Match;

typedef struct hash_ {
  int word_length;
  int size_hash;
  int seq1_len;
  int seq2_len;
  int *values1;
  int *values2;
  int *counts;
  int *last_word;
  int *diag;
  int *hist;
  char *seq1;
  char *seq2;
  int *expected_scores;
  Diag_Match *diag_match;
  Block_Match *block_match;
  int max_matches;
  int matches;
  int min_match;
  int fast_mode;
  int filter_words;
} Hash;

void set_hash8_lookupn(void);

int init_hash8n (
		int max_seq, 
		int max_diagonal,
		int word_length,
		int max_matches,
		int min_match,
		int job,
		Hash **h);


void free_hash_pos ( Hash *h );

void free_hash8n ( Hash *h );

int hash_word8n ( char *seq, int *start_base, int seq_len, int word_length,
	      unsigned short *uword);

int hash_seq8n ( char *seq, int *hash_values, int seq_len, int word_length);

int hash_word4n ( char *seq, int *start_base, int seq_len, int word_length,
	      unsigned char *uword);

int hash_seq4n ( char *seq, int *hash_values, int seq_len, int word_length);

void store_hashn ( Hash *h );
void store_hashn_nocount ( Hash *h );

int hash_seqn (Hash *h, int job);

int diagonal_length(int seq1_len, int seq2_len, int diagonal_number);

void diagonal_intercepts (int diagonal_number, int seq1_len, int seq2_len,
			  int *seq1_intercept, int *seq2_intercept );

void p_comp(double comp[], char *seq, int seq_len);

int poisson_diagonals(int min_diag, int max_diag, int word_size,
		      double max_prob, int *expected_scores, double comp[]);

int best_intercept ( Hash *h, int *seq1_i, int *seq2_i );

int compare_b_bulk(Hash *h,
		   ALIGN_PARAMS *params, OVERLAP *overlap,
		   int cnum, tg_rec crec2,
		   Contig_parms *contig_list, int number_of_contigs,
		   int ignore_after1, HashTable *links,
		   void (*add_func)(OVERLAP *overlap,
				    int cnum1,
				    int cnum2,
				    void *clientdata),
		   void *add_data);
int compare_b(Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap);
int compare_a(Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap);

int compare_seqsb(Hash *h);

int align_blocks(Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap);

int compare_seqs(Hash *h, int *seq1_match_pos, int *seq2_match_pos, int *match_length);

int reps (Hash *h, int **seq1_match_pos, int **seq2_match_pos, int **match_len,
	  int offset, char sense);
int reps_nocount(Hash *h, int **seq1_match_pos, int **seq2_match_pos,
		 int **match_length, int offset, char sense);

int set_band_blocks(int seq1_len, int seq2_len);

#endif
