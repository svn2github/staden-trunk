#if !defined(SP_ALIGNMENT_STRUCTS_H)
#define SP_ALIGNMENT_STRUCTS_H


namespace sp {


#define SP_ALIGNMENT_ALGORITHM_A 1  /* dynamic programming (affine_align) */
#define SP_ALIGNMENT_ALGORITHM_B 17 /* block alignment (compare_b) */
#define SP_ALIGNMENT_ALGORITHM_C 31 /* poisson scores  (compare_c) */
#define SP_ALIGNMENT_ALGORITHM_D 17 /* unused (no alignment) (compare_d) */


#ifndef DYNMAT
typedef int W128_D[128][128];
/*typedef int (*W128_P)[128][128];*/
typedef W128_D *W128_P;
#endif


typedef struct Overlap {
    double  percent;
    int     length;
    int     direction;
    int     lo, ro; /* Left and right offsets */
    int     left1, left2, left, right1, right2, right;
    double  score;
    double  qual;
    int    *S; /* Alignment edit buffer */
    int     s_len;
    int    *S1, *S2; /* Sequence edit buffers for seq1 and seq2 */
    int     s1_len, s2_len; /* lengths of sequence edit buffers */
    int     seq1_len, seq2_len; /* lengths of the original sequences */
    char   *seq1, *seq2; /* the original sequences */
    char   *seq1_out;
    char   *seq2_out;
    int     seq_out_len;
    /*
     * example overlap: 
     *
     * 01234567890123456
     * AACGTAAT**CGCT***
     * ****TATTGCCGCTAAG
     *
     * left1, right1, left2, right2, left, right are all positions in the alignment
     * left1  = first non-pad in seq1 (0)
     * right1 = last non-pad in seq1 (13)
     * left2  = first non-pad in seq2 (4)
     * right2 = last non-pad in seq2 (16)
     * left   = position to right of left overhang ie max(left1,left2) (4)
     * right  = position to left of right overhang ie min(right1,right2) (13)
     * length = length of overlapping section of the overlap = right - left + 1 (10)
     *
     *left1       right1
     *    AACGTAAT**CGCT***
     *    01234567890123456
     *    ****TATTGCCGCTAAG
     *    left2      right2
     *        <-length->
     *     left    right
     *
     * direction = describes the type of overlap:
     * 0 - suffix of seq1 overlaps with prefix of seq2
     * 11111111
     *     2222222
     * 1 - suffix of seq2 overlaps with prefix of seq1
     *     111111111
     * 22222222
     * 2 - seq1 contains seq2
     * 11111111111  111111111111  11111111111
     *     2222     222222        22222222222
     * 3 - seq2 contains seq1
     *     11111     11111111111
     * 22222222222   2222222222222
     *
     * for example overlap direction (0)
     * if(overlap->left1 == overlap->left2)
     *   overlap->direction = (overlap->right1 >= overlap->right2) ? 2 : 3;
     * else if(overlap->left1 < overlap->left2)
     *   overlap->direction = (overlap->right1 >= overlap->right2) ? 2 : 0;
     * else
     *   overlap->direction = (overlap->right1 <= overlap->right2) ? 3 : 1;
     *
     * lo and ro (left and right offsets): the lengths of the overhangs at each end:
     * depends on the direction and can be +ve or -ve:
     * a containment overlap will have +ve left_offset and -ve right_offset, 
     * a non-containment overlap will have +ve left_offset and +ve right_offset
     * 
     * 11111111
     *     2222222
     * lo (4) ro (3)
     *     111111111
     * 22222222
     * lo (4) ro (5)
     * 11111111111  lo (4) ro (3)  111111111111 lo (0) ro (-6) 11111111111 lo (0) ro (0)
     *     2222                    222222                      22222222222
     *
     *     11111  lo (4) ro (-2)   11111111111    lo (0) ro (-2)
     * 22222222222                  2222222222222
     *
     * switch(overlap->direction) {
     *  case 0: case 2:
     *    overlap->lo = overlap->left2 - overlap->left1;
     *    overlap->ro = overlap->right2 - overlap->right1;
     *    break;
     *  case 1: case 3:
     *    overlap->lo = overlap->left1 - overlap->left2;
     *    overlap->ro = overlap->right1 - overlap->right2;
     *    break;
     *  default:
     *    break;
     *
     * The formats of the edit buffers are as follows, using this
     * alignment as an example:
     *
     * AACGTAAT**CGCT***
     * ****TATTGCCGCTAAG
     *
     *    There are two sequence edit buffers, one for each sequence
     *    in the alignment. A sequence edit buffer contains +n for
     *    n bases from the sequence, and -n for n pads in the sequence.
     *    This would give S1 = {+8, -2, +4, -3} and S2 = {-4, 13} for
     *    the example alignment.
     *
     *    The alignment edit buffer contains '0' when two bases are
     *    aligned, +n for n pads in seq1, and -n for n pads in seq2.
     *    This would give S[] = {-4, 0, 0, 0, 0, +2, 0, 0, 0, 0, +3}
     *    for the example alignment. They are included only for
     *    compatibility with other programs in the package and are not
     *    used (or created).
     *
     */
} OVERLAP;

typedef struct diag_match_ {
    int pos;
    double prob;
} Diag_Match;

typedef struct block_match_ {
    int pos_seq1;
    int pos_seq2;
    int diag;
    int length;
    int best_score;
    int prev_block;
} Block_Match;

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
} Hash;

typedef struct Align_params {
    int band;
    int gap_open;
    int gap_extend;
    int edge_mode;
    int return_job;
    int seq1_start;
    int seq2_start;
    int seq1_end;
    int seq2_end;
    int first_row;
    int band_left;
    int band_right;
    char old_pad_sym;
    char new_pad_sym;
    int algorithm;
    int word_length;
    int min_match;
    double max_prob;
#ifdef DYNMAT
    int** score_matrix;
#else
    W128_P score_matrix;
#endif
    Hash *hash;
} ALIGN_PARAMS;
/*    int (*score_matrix)[128][128];*/



}



#endif

