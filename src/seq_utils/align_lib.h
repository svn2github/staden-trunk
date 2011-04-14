#if !defined(ALIGN_LIB_H)
#define ALIGN_LIB_H

typedef struct Mseg {
  char *seq;
  int length;
  int offset;
} MSEG;

typedef struct Contigl {
  MSEG *mseg;
  struct Contigl *next;
  int id;
}CONTIGL;

CONTIGL *create_contig_link(void);

void print_contig_links(CONTIGL *contigl);

/*
 * FIXME: See shuffle_pads.c header for possibile ways to improve the
 * multiple alignment code here.
 */

typedef struct Malign {
    char *charset;
    int charset_size;
    int  length;
    int  **matrix;
    CONTIGL *contigl;
    char *consensus;
    int *orig_pos;
    int **counts;
    int **scores;
    int *orig_counts;
    int *orig_scores;
    int orig_length;
    int gap_open;
    int gap_extend;
} MALIGN;

MSEG **get_malign_segs(CONTIGL *contigl);

void print_malign_seqs(MALIGN *malign);

int set_malign_charset(MALIGN *malign, char *charset);

int **create_malign_counts(int length, int depth);

void get_malign_counts (MALIGN *malign, int start, int end);

void print_malign_counts(MALIGN *malign);

void print_malign_matrix(MALIGN *malign);

void print_malign_scores(MALIGN *malign);

void print_contig_links (CONTIGL *contigl);

int contigl_length (CONTIGL *contigl);

MALIGN *contigl_to_malign(CONTIGL *contigl, int gap_open, int gap_extend);

void destroy_malign (MALIGN *malign, int contigl_too);

MSEG *create_mseg(void);

void init_mseg (MSEG *mseg, char *seq, int length, int offset);

void destroy_mseg (MSEG *mseg);

void get_malign_consensus(MALIGN *malign, int start, int end);

void set_malign_lookup (int charset_size);

void malign_insert_scores(MALIGN *malign, int pos, int size);

void malign_remove_contigl(MALIGN *malign, CONTIGL *previous, CONTIGL *del);

void malign_add_contigl(MALIGN *malign, CONTIGL *previous, CONTIGL *add);

void malign_recalc_scores(MALIGN *malign, int start, int end);

typedef struct Moverlap {
    double	percent;
    int		length;
    int 	direction;
    int		lo, ro; /* Left and right offsets */
    int		left1, left2, left, right1, right2, right;
    double 	score;
    double 	qual;
    int	*S; /* Alignment edit buffer */
    int		s_len;
    int	*S1, *S2; /* Sequence edit buffers for seq1 and seq2 */
    int		s1_len, s2_len; /* lengths of sequence edit buffers */
    int		malign_len, seq2_len; /* lengths of the original sequences */
  MALIGN *malign;
  MALIGN *malign_out;
    char        *seq2; /* the original sequence */
    char        *seq1_out, *seq2_out; 
    int    seq_out_len;
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
     *     11111  lo (-4) ro (-2)   11111111111    lo (0) ro (-2)
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
} MOVERLAP;

MOVERLAP *create_moverlap(void);

void init_moverlap (MOVERLAP *moverlap, MALIGN *malign, char *seq2, 
int malign_len, int seq2_len);

void destroy_moverlap (MOVERLAP *moverlap);


typedef struct Overlap {
    double	percent;
    int		length;
    int 	direction;
    int		lo, ro; /* Left and right offsets */
    int		left1, left2, left, right1, right2, right;
    double 	score;
    double 	qual;
    int	*S; /* Alignment edit buffer */
    int		s_len;
    int	*S1, *S2; /* Sequence edit buffers for seq1 and seq2 */
    int		s1_len, s2_len; /* lengths of sequence edit buffers */
    int		seq1_len, seq2_len; /* lengths of the original sequences */
    char        *seq1, *seq2; /* the original sequences */
    char *seq1_out;
    char *seq2_out;
    int    seq_out_len;
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
     *     11111  lo (-4) ro (-2)   11111111111    lo (0) ro (-2)
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


#define RETURN_SEQ 1
#define RETURN_EDIT_BUFFERS 2
#define RETURN_EDIT_BUFFER  4
#define RETURN_NEW_PADS 8
#define RETURN_END_GAPS 16

/* Score gaps at the start of seqs */
#define EDGE_GAPS_COUNT   1

/* No edge gaps for start.  *Must* set _COUNT or _ZERO, despite bit field */
#define EDGE_GAPS_ZERO    2

/* Whether to only go from (0,0)-(N,M) or end at any point on edges */
#define FULL_LENGTH_TRACE 4

/* Don't score gaps at seq ends */
#define BEST_EDGE_TRACE   8

/* Penalise in Y, but not in X */
#define EDGE_GAPS_ZEROX   16

typedef struct Align_params {
    int gap_open;
    int gap_extend;
    int band;
    int first_row;
    int band_left;
    int band_right;
    int edge_mode;
    int job;
    char old_pad_sym;
    char new_pad_sym;
    int seq1_start;
    int seq2_start;
} ALIGN_PARAMS;


typedef struct Seg {
    int length;
    char *seq;
} SEG;

SEG *create_seg(void);

void destroy_seg (SEG *seg);


int get_segment( OVERLAP *overlap, SEG *seg,
		 int job );

int set_band(int seq1_len, int seq2_len);

OVERLAP *create_overlap(void);

void init_overlap (OVERLAP *overlap, char *seq1, char *seq2, int seq1_len,
		   int seq2_len);

void destroy_overlap (OVERLAP *overlap);

void free_overlap (OVERLAP *overlap);

ALIGN_PARAMS *create_align_params(void);

int set_align_params (ALIGN_PARAMS *params, int band, int gap_open, 
		       int gap_extend, int edge_mode, int job, int seq1_start, int seq2_start, char old_pad_sym, char new_pad_sym, int set_job);

void destroy_alignment_params (ALIGN_PARAMS *params);

void print_char_array ( FILE *file, char *array, int array_len);

int affine_align(OVERLAP *overlap,
	     ALIGN_PARAMS *params);

int affine_malign(MOVERLAP *moverlap,
	     ALIGN_PARAMS *params);

int realigner_malign(MOVERLAP *moverlap,
		     ALIGN_PARAMS *params);

int seq_to_edit ( char *seq, int seq_len, int **S_out, int *S_len,
		 char PAD_SYM);
void old_pads_for_new ( char *seq, int seq_len, char OLD_PAD_SYM, char NEW_PAD_SYM);

int seq_to_overlap ( OVERLAP *overlap, char OLD_PAD_SYM, char NEW_PAD_SYM);

int seq_to_moverlap ( MOVERLAP *moverlap, char OLD_PAD_SYM, char NEW_PAD_SYM);

int do_trace_back ( unsigned char *bit_trace, char *seq1, char *seq2,
		   int seq1_len, int seq2_len, char **seq1_out_ret, char **seq2_out_ret,
		   int *seq_out_len, int b_r, int b_c, int b_e,
		   int band, int first_band_left, int first_row, 
		   int band_length, char PAD_SYM );


void seq_expand(char	  *seq,
		char	  *seq_align,
		int	  *seq_align_len,
		int *S,
		int       s_len,
		int       mode,
		char PAD_SYM);

int print_alignment(char	*seq1,
		    char	*seq2,
		    int		seq1_len,
		    int		seq2_len,
		    int	*S1,
		    int	*S2,
		    int		s1_len,
		    int		s2_len,
		    double	score,
		    FILE	*fpt);

int print_overlap ( OVERLAP *overlap, FILE *fp);

void print_fasta(char *description, char *seq, FILE *fpt);

void init_W128(int **matrix,
	       char *order, int unknown);

int set_alignment_matrix ( char *fn, char *base_order );

#endif
