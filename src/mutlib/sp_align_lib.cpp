#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <os.h>
#include <staden.h>
#include <sp_alignment_structs.h>
#include <sp_alignment.h>
#include <sp_align_lib.h>
#include <read_matrix.h>


namespace sp {


int SIZE_MAT;

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
 OVERLAP *overlap;

 ierr = set_alignment_matrix("/home5/pubseq/share/tables/nuc_matrix",
 "ACGT");
 if(ierr) return -1;

 if (NULL == (overlap = create_overlap())) return -1;

 init_overlap (overlap, seq1, seq2, seq1_len, seq2_len);
   
 if (NULL == (params = create_align_params())) return -1;

 if (set_align_params (params, band, gap_open, gap_extend, edge_mode, job,
 seq1_start, seq2_start, '*','!',0)) {
 destroy_overlap (overlap);
 destroy_align_params (params);
 return -1;
 };

 ierr = affine_align(overlap, params);

 print_overlap(overlap,stdout);

 destroy_overlap (overlap);
 destroy_align_params (params);
*/

#define BYTE_ACROSS  1
#define BYTE_DOWN 2
#define BYTE_DIAGONAL   3


/**
 * convert a user edge mode to an internal mode used by alignment routines
 */

void to_internal_edges(int user, int *edge_mode) {

    int temp = 0;

    if (0 == user) {
	*edge_mode = EDGE_GAPS_ZERO | FULL_LENGTH_TRACE;
	return;
    }

    if(user & SP_ALIGNMENT_LEFT_EDGE_GAPS_COUNT) {
	temp |= EDGE_GAPS_COUNT;
    }
    else {
	temp |= EDGE_GAPS_ZERO;
    }
    if(user & SP_ALIGNMENT_BEST_RIGHT_EDGE) {
	temp |= BEST_EDGE_TRACE;
    }
    else {
	temp |= FULL_LENGTH_TRACE;
    }
    *edge_mode = temp;

}

/**
 * converts an input matrix of size strlen(order) x strlen(order)
 * to 128 x 128 and returns it in matrix_128
 */

#ifdef DYNMAT
void to_128( int** matrix_128, int **matrix, char *order, int unknown )
#else
    void to_128(int matrix_128[128][128], int **matrix, char *order, int unknown)
#endif
{
    int i, j, i_end;
    unsigned char ci, cj;

    for (i = 0; i < 128; i++) {
        for (j = 0; j < 128; j++) {
	    matrix_128[i][j] = unknown;
        }
    }

    i_end = strlen(order);
    for (i = 0; i < i_end; i++) {
	ci = order[i];
	for (j = 0; j < i_end; j++) {
	    cj = order[j];
	    matrix_128[        ci ][        cj ] = matrix[i][j];
	    matrix_128[tolower(ci)][        cj ] = matrix[i][j];
	    matrix_128[        ci ][tolower(cj)] = matrix[i][j];
	    matrix_128[tolower(ci)][tolower(cj)] = matrix[i][j];
	}
    }
    /*print_128(matrix_128);*/

}

/**
 * malloc a SEG structure
 */

SEG *create_seg(void) {
    SEG *seg;

    if(NULL == (seg = (SEG *) xmalloc(sizeof(SEG)))) {
	verror(ERR_WARN, "create_seg", "xmalloc failed");
	return NULL;
    }

    seg->seq = NULL;
    return seg;
}

/**
 * destroy a SEG structure
 */

void destroy_seg (SEG *seg) {
    if ( seg ) {
	if ( seg->seq ) xfree ( seg->seq );
	xfree ( seg );
    }
}

/**
 * given an edit buffer and an original sequence
 * return a padded sequence
 */

void seq_expand(char   *seq,
		char    *seq_align,
		int     *seq_align_len,
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
    
    int     i, j;
    int     s, s_start, s_end;
    int  l;


    /*for(i=0;i<s_len;i++)printf("%d %d\n",i,S[i]);*/
    
    s_end = s_len;
    if((0 == mode) || (1 == mode)) {
	/* Ignore right-hand overhang */
	for(s = s_len - 1; s >= 0; s--) {
	    if(S[s] > 0) {
		s_end = s + 1;
		break;
	    }
	}
    }
    
    s_start = 0;
    if((0 == mode) || (2 == mode)) {
	/* Ignore left-hand overhang */
	for(s = 0; s < s_len; s++) {
	    if(S[s] > 0) {
		s_start = s;
		break;
	    }
	}
    }

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

/**
 * sum the scores for the aligned characters in an overlap output
 */

#ifdef DYNMAT
int overlap_score (OVERLAP *overlap, int** score_matrix )
{
    int i,s;

    for(s=0,i=overlap->left;i<=overlap->right;i++) {
	s += score_matrix[(int)overlap->seq1_out[i]][(int)overlap->seq2_out[i]];
    }
    return s;
}
#else
int overlap_score (OVERLAP *overlap, W128_P score_matrix)
{
    int i,s;

    for(s=0,i=overlap->left;i<=overlap->right;i++) {
	s += (*score_matrix)[(int)overlap->seq1_out[i]][(int)overlap->seq2_out[i]];
    }
    return s;
}
#endif

int get_segment( OVERLAP *overlap, SEG *seg,
		 int job ) {

    int     len_align;
    
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

/**
 * prints the alignment from an overlaps edit buffers
 */

int print_alignment(char   *seq1,
		    char *seq2,
		    int     seq1_len,
		    int     seq2_len,
		    int  *S1,
		    int  *S2,
		    int     s1_len,
		    int     s2_len,
		    double  score,
		    FILE *fpt) {
    
    char *seq1_align, *seq2_align;
    int     seq1_align_len, seq2_align_len;
    char temp_seq[51];
    int     i, j;
    int     max_out_width;
    int     max_seq;
    int     len_align;
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

/**
 * prints a sequence in fasta format
 */

void print_fasta(char *description, char *seq, FILE *fpt) {
    int     i;
    int     seq_len;
    char temp[61];

    fprintf(fpt, ">%s\n", description);

    seq_len = strlen(seq);
    for(i = 0; i < seq_len; i += 60) {
	memset(temp, '\0', 61);
	strncpy(temp, seq + i, 60);
	fprintf(fpt, "%s\n", temp);
    }
}


/**
 * find the ends of a sequence seq defined as being the first
 * non pad positions in from each end
 * these are returned in left and right which are offsets from 0
 */

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

/**
 * swaps old pad symbols for new ones
 */

void old_pads_for_new ( char *seq, int seq_len, char OLD_PAD_SYM, char NEW_PAD_SYM) {
    int i;
    for(i=0;i<seq_len;i++) {
	if(seq[i]==NEW_PAD_SYM)seq[i]=OLD_PAD_SYM;
    }
}

/**
 * from a sequence alignment stored in an overlap work out and fill in
 * all the other overlap values
 */

int seq_to_overlap ( OVERLAP *overlap, char OLD_PAD_SYM, char NEW_PAD_SYM) {

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

/**
 * given an edit buffer and its length
 * where possible combine adjacent edits into one
 * and return the new length
 */

void shrink_edit_buffer(int *S, int *s_len) {

    int len, value, v, sign, new_sign, i, k;

    len = *s_len;
    value = S[0];
    sign = (value > 0) ? 1 : 0;
    for(i=1,k=0;i<len;i++) {
	v = S[i];
	new_sign = (v > 0) ? 1 : 0;
	if(new_sign == sign) {
	    value += v;
	}
	else {
	    S[k++] = value;
	    value = v;
	    sign = new_sign;
	}
    }
    S[k++] = value;
    *s_len = k;
}

/**
 * given an overlap shrink its edit buffers and reset their lengths
 */
void shrink_edit_buffers(OVERLAP *overlap) {

    shrink_edit_buffer(overlap->S1, &overlap->s1_len);
    shrink_edit_buffer(overlap->S2, &overlap->s2_len);
}

/**
 * from a padded sequence calculate its edit buffer
 * returned in S_out with length S_len
 */
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

/**
 * given a bit trace from an alignment routine create the sequence alignment
 * returned in seq1_out_ret and seq2_out_ret with length seq_out_len
 */

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
	    seq1_out--;
	    *(seq2_out--) = seq2[r];
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

/**
 * given an unsigned char trace from an alignment routine 
 * create the sequence alignment
 * returned in seq1_out_ret and seq2_out_ret with length seq_out_len
 */
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
	if((byte<0)||(byte>=SIZE_MAT))
	    printf("SCREAM trace SIZE_MAT %d byte %d seq1_len %d seq2_len %d fbl %d band %d bl %d fr %d\n",SIZE_MAT,byte,seq1_len,seq2_len,first_band_left,band,band_length,first_row);
	trace_byte = bit_trace[byte];
      
	if(trace_byte == BYTE_DIAGONAL) {
	    r--,c--;
	    *(seq1_out--) = seq1[c];
	    *(seq2_out--) = seq2[r];
	}
	else if(trace_byte == BYTE_DOWN) {
	    r--;
	    seq1_out--;
	    *(seq2_out--) = seq2[r];
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

/**
 * destroy all the memory used by affine alignment routines
 * and set seq1_out and seq2_out to NULL so callers know not
 * to use them!
 */

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
    seq1_out = NULL;
    seq2_out = NULL;
}

/**
 * dynamic programming routine using 3 tables
 */

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
#ifndef DYNMAT
    W128_P W128_p;
    W128_p = params->score_matrix;
#endif


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

#ifdef DYNMAT
	    score_matrix_p = params->score_matrix[ seq2[row-1] ];
#else
	    score_matrix_p = (int *) (*W128_p + (int)seq2[row-1]);
#endif

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

#ifdef DYNMAT
	    score_matrix_p = params->score_matrix[ seq2[row-1] ];
#else
	    score_matrix_p = (int *) (*W128_p + (int)seq2[row-1]);
#endif

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
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_EDIT_BUFFERS ) {
	if (seq_to_edit ( seq1_out,seq_out_len,&overlap->S1,&overlap->s1_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	if (seq_to_edit ( seq2_out,seq_out_len,&overlap->S2,&overlap->s2_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_SEQ ) {
	if ( !(params->return_job & SP_ALIGNMENT_RETURN_NEW_PADS) ) {
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

/**
 * dynamic programming routine using 3 tables
 * and bit patterns to store the traceback
 * suitable for long sequences (ie less memory) but slower
 * than affine_align_big
 */
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
#ifndef DYNMAT
    W128_P W128_p;
    W128_p = params->score_matrix;
#endif


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


#ifdef DYNMAT
	    score_matrix_p = params->score_matrix[ seq2[row-1] ];
#else
	    score_matrix_p = (int *) (*W128_p + (int)seq2[row-1]);
#endif

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
	row = last_row;
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
	if(best_edge_score < 0) printf("scream: best_edge_score not set\n");
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

#ifdef DYNMAT
	    score_matrix_p = params->score_matrix[ seq2[row-1] ];
#else
	    score_matrix_p = (int *) (*W128_p + (int)seq2[row-1]);
#endif

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
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_EDIT_BUFFERS ) {
	if (seq_to_edit ( seq1_out,seq_out_len,&overlap->S1,&overlap->s1_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	if (seq_to_edit ( seq2_out,seq_out_len,&overlap->S2,&overlap->s2_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_SEQ ) {
	if ( !(params->return_job & SP_ALIGNMENT_RETURN_NEW_PADS) ) {
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


/**
 * dynamic programming routine using 3 tables
 * and unsigned char to store the traceback
 * suitable for not too long sequences (else too much memory) and faster
 * than affine_align_bits (which should be used for long sequences
 * what is long??
 */
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
    int MAX_E = 0;
    int best_max = 0;
    int best_max_row = 0;
    int best_max_col = 0;
#ifndef DYNMAT
    W128_P W128_p;
    W128_p = params->score_matrix;
#endif


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

    /*
      print_align_params(params);
      print_overlap_struct(overlap);
      printf("seq1_len %d seq2_len %d\n",overlap->seq1_len,overlap->seq2_len);
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
	/*
	  printf("size_mat %d band %d band_left %d first_row %d band_length %d\n",size_mat,band,band_left,first_row,band_length);
	*/
	SIZE_MAT = size_mat + 1;
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

	    /*score_matrix_p = W128[(int)seq2[row-1]];*/
	    /*score_matrix_p = *(W128_p)[(int)seq2[row-1]];*/

#ifdef DYNMAT
	    score_matrix_p = params->score_matrix[ seq2[row-1] ];
#else
	    score_matrix_p = (int *) (*W128_p + (int)seq2[row-1]);
#endif

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
		/*
		  printf("row %d column %d\n",row,column);
		*/

		/* have to make sure which row we actually reached ! */

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

		/*     e_row = (row - first_row + 1) * band_length;*/
		/*     e_col = column - band_left + 1;*/
		e = e_row + e_col;

		/*
		  printf("row %d column %d first_row %d band_left %d e %d e_row %d e_col %d MAX_E %d max_row %d\n",
		  row,column,first_row,band_left,e,e_row,e_col,MAX_E,max_row);
		*/

		MAX_E = MAX(e,MAX_E);
       
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

		{
		    int best = 0;
		    if (best < best_F1)
			best = best_F1;
		    if (best < best_G1)
			best = best_G1;
		    if (best < best_H1)
			best = best_H1;
		    if (best_max < best) {
			best_max = best;
			best_max_row = row;
			best_max_col = column;
		    }
		}

	    }
	    if((e<0)||(e>=SIZE_MAT))

		printf("SCREAM 1trace SIZE_MAT %d e %d seq1_len %d seq2_len %d fbl %d band %d bl %d fr %d\n",SIZE_MAT,e,seq1_len,seq2_len,first_band_left,band,band_length,first_row);
     
	    if ( column > seq1_len ) {
		if ( edge_mode & BEST_EDGE_TRACE ) {
		    best_H1 = MAX(best_H1,best_G1);
		    best_F1 = MAX(best_H1,best_F1);
		    if ( best_F1 > best_edge_score ) {
			best_edge_score = best_F1;
			b_r = row;
			b_e = ((row - first_row + 1) * band_length) + 
			    (seq1_len - band_left + 1);
			if((b_e<0)||(b_e>=SIZE_MAT))

			    printf("SCREAM 22trace SIZE_MAT %d b_e %d seq1_len %d seq2_len %d fbl %d band %d bl %d fr em %d%d\n",SIZE_MAT,b_e,seq1_len,seq2_len,first_band_left,band,band_length,first_row,edge_mode);
		    }
		}
	    }
	}

	/*
	  printf("row %d max_col %d first_row %d band_left %d max_row %d\n",row,max_col,first_row,band_left,max_row);
	*/

	last_column = max_col;
	row = last_row;
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
		    if((b_e<0)||(b_e>=SIZE_MAT))

			printf("SCREAM 0trace SIZE_MAT %d b_e %d seq1_len %d seq2_len %d fbl %d band %d bl %d fr %d edge_mode %d\n",SIZE_MAT,b_e,seq1_len,seq2_len,first_band_left,band,band_length,first_row,edge_mode);
		}
	    }
	}
	/*
	  printf("row %d last_column %d b_e %d\n",row,last_column,
	  ((last_row - 1 - first_row + 1) * band_length) + 
	  (last_column - band_left + 1));
	*/

	if ( edge_mode & FULL_LENGTH_TRACE ) {
	    b_r = last_row;
	    /*b_c = seq1_len;*/
	    b_c = last_column;
	    b_e = ((b_r - first_row + 1) * band_length) + 
		(b_c - band_left + 1);
	    if((b_e<0)||(b_e>=SIZE_MAT))
		printf("SCREAM 1 trace SIZE_MAT %d b_e %d seq1_len %d seq2_len %d fbl %d band %d bl %d fr %d edge_mode %d\n",SIZE_MAT,b_e,seq1_len,seq2_len,first_band_left,band,band_length,first_row,edge_mode);
	    best_edge_score = b_s;
	}
	if ( !(edge_mode & BEST_EDGE_TRACE) ) { /* fall back */
	    b_r = last_row;
	    /*b_c = seq1_len;*/
	    b_c = last_column;
	    b_e = ((b_r - first_row + 1) * band_length) + 
		(b_c - band_left + 1);
	    if((b_e<0)||(b_e>=SIZE_MAT))
		printf("SCREAM 2 trace SIZE_MAT %d b_e %d seq1_len %d seq2_len %d fbl %d band %d bl %d fr %d edge_mode %d\n",SIZE_MAT,b_e,seq1_len,seq2_len,first_band_left,band,band_length,first_row,edge_mode);
	    best_edge_score = b_s;
        }

	if(best_edge_score < -100000) printf("SCREAM best edge not set\n");
    }
    else {

	/* Initialise the bit trace */
	size_mat = (seq1_len + 1) * (seq2_len + 1);
	SIZE_MAT = size_mat + 1;
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

#ifdef DYNMAT
	    score_matrix_p = params->score_matrix[ seq2[row-1] ];
#else
	    score_matrix_p = (int *) (*W128_p + (int)seq2[row-1]);
#endif

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
			    band, first_band_left, first_row, band_length, NEW_PAD_SYM)) {
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }

    overlap->seq1_out = seq1_out;
    overlap->seq2_out = seq2_out;
    overlap->seq_out_len = seq_out_len;

    if ( i = seq_to_overlap (overlap, OLD_PAD_SYM, NEW_PAD_SYM)) {
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_EDIT_BUFFERS ) {
	if (seq_to_edit ( seq1_out,seq_out_len,&overlap->S1,&overlap->s1_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	if (seq_to_edit ( seq2_out,seq_out_len,&overlap->S2,&overlap->s2_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_SEQ ) {
	if ( !(params->return_job & SP_ALIGNMENT_RETURN_NEW_PADS) ) {
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

/**
 * unused alignment routine using only 2 tables and bit pattern trace back
 * meant to be OK for reasonable similar sequences but not good enough
 * for our padded ones at least
 */

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
#ifndef DYNMAT
    W128_P W128_p;
    W128_p = params->score_matrix;
#endif


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

#ifdef DYNMAT
	    score_matrix_p = params->score_matrix[ seq2[row-1] ];
#else
	    score_matrix_p = (int *) (*W128_p + (int)seq2[row-1]);
#endif

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

#ifdef DYNMAT
	    score_matrix_p = params->score_matrix[ seq2[row-1] ];
#else
	    score_matrix_p = (int *) (*W128_p + (int)seq2[row-1]);
#endif

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
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_EDIT_BUFFERS ) {
	if (seq_to_edit ( seq1_out,seq_out_len,&overlap->S1,&overlap->s1_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	if (seq_to_edit ( seq2_out,seq_out_len,&overlap->S2,&overlap->s2_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_SEQ ) {
	if ( !(params->return_job & SP_ALIGNMENT_RETURN_NEW_PADS) ) {
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

/**
 * the proper route into the affine alignment routines
 * uses faster method for sequences of reasonable length
 * and slower for longer ones (also takes into account band size)
 * memory calculated as 2 * band * MIN(seq1_len,seq2_len)
 * or seq1_len * seq2_len
 * use fast method if memory < MAX_MEMORY (10,000,000)
 * overlap and params must be set up beforehand
 * alignment returned in overlap
 */

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

/**
 * unused alignment routine using only 2 tables and unsigned char trace back
 * meant to be OK for reasonable similar sequences but not good enough
 * for our padded ones at least
 */

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
#ifndef DYNMAT
    W128_P W128_p;
    W128_p = params->score_matrix;
#endif


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


#ifdef DYNMAT
	    score_matrix_p = params->score_matrix[ seq2[row-1] ];
#else
	    score_matrix_p = (int *) (*W128_p + (int)seq2[row-1]);
#endif

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

		/*     e_row = (row - first_row + 1) * band_length;*/
		/*     e_col = column - band_left + 1;*/
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

#ifdef DYNMAT
	    score_matrix_p = params->score_matrix[ seq2[row-1] ];
#else
	    score_matrix_p = (int *) (*W128_p + (int)seq2[row-1]);
#endif

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
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_EDIT_BUFFERS ) {
	if (seq_to_edit ( seq1_out,seq_out_len,&overlap->S1,&overlap->s1_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	if (seq_to_edit ( seq2_out,seq_out_len,&overlap->S2,&overlap->s2_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_SEQ ) {
	if ( !(params->return_job & SP_ALIGNMENT_RETURN_NEW_PADS) ) {
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


int malign_lookup[256];

int set_malign_charset(MALIGN *malign, const char *charset) {

    if(NULL == (malign->charset = (char *) xmalloc(sizeof(charset)+1))) {
	verror(ERR_WARN, "set_malign_charset", "xmalloc failed");
	return -1;
    }
    strcpy(malign->charset,charset);
    return 0;
}

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
	    /* FIXME malign->matrix[jjj][iii] = W128[jj][ii]; */
	}
    }
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
    return contigl;
}

MALIGN *create_malign(void) {
    MALIGN *malign;

    if(NULL == (malign = (MALIGN *) xmalloc(sizeof(MALIGN)))) {
	verror(ERR_WARN, "create_malign", "xmalloc failed");
	return NULL;
    }

    malign->contigl = NULL;
    malign->msegs = NULL;
    malign->nseqs = 0;
    malign->consensus = NULL;
    malign->scores = NULL;
    malign->matrix = NULL;
    malign->charset_size = 6; /*  a,c,g,t,*,n */
    return malign;
}

void destroy_malign (MALIGN *malign) {
    if ( malign ) {
	if ( malign->contigl ) xfree ( malign->contigl );
	if ( malign->msegs ) xfree ( malign->msegs );
	if ( malign->consensus ) xfree ( malign->consensus );
	if ( malign->scores ) xfree ( malign->scores );
	if ( malign->matrix ) xfree ( malign->matrix );
	xfree ( malign );
    }
}
void free_malign (MALIGN *malign) {
    if ( malign ) {
	if ( malign->contigl ) xfree ( malign->contigl );
	if ( malign->msegs ) xfree ( malign->msegs );
	if ( malign->consensus ) xfree ( malign->consensus );
	if ( malign->scores ) xfree ( malign->scores );
    }
    malign->contigl = NULL;
    malign->msegs = NULL;
    malign->consensus = NULL;
    malign->scores = NULL;
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

int **create_malign_counts(int length, int depth) {
    /* length is contig length, depth charset_size */
    int **counts;
    int i;

    counts = (int **)malloc(length*sizeof(int *));
    for(i=0;i<length;i++) {
	counts[i] = (int *)calloc(depth,sizeof(int));
    }
    return counts;
}
 
void get_malign_counts (MALIGN *malign) {
    CONTIGL *t;
    int i,j,k,l;

    t = malign->contigl;
    while(t) {
	for(j=0,k=t->mseg->offset;j<t->mseg->length;j++,k++) {
	    l = malign_lookup[(int)(t->mseg->seq[j])];
	    malign->scores[k][l]++;
	}
	t = t->next;
    }
    /* put the totals in the rows for gap open and gap extend */
    k = malign->charset_size;
    l = k + 1;
    for(i=0;i<malign->length;i++) {
	for(j=0;j<malign->charset_size;j++) {
	    malign->scores[i][k] += malign->scores[i][j];
	    malign->scores[i][l] += malign->scores[i][j];
	}
    }
}

void get_malign_consensus(MALIGN *malign) {
    int i,j,k;
    char *s;

    /* must be run after get_malign_scores
     * and before scale_malign_scores ! 
     */

    s = malign->consensus;
    k = malign->charset_size;
    for(i=0;i<malign->length;i++) {
	s[i] = '-';
	for(j=0;j<malign->charset_size;j++) {
	    if(malign->scores[i][j] == malign->scores[i][k]) {
		s[i] = malign->charset[j];
		break;
	    }
	}
    }
}

void print_malign_scores(MALIGN *malign) {
    int i,j;
    for(j=0;j<malign->charset_size+2;j++) {
	for(i=0;i<malign->length;i++) {
	    printf(" %d ",malign->scores[i][j]);
	}
	printf("\n");
    }
    printf("\n");
}

void scale_malign_scores(MALIGN *malign, int gap_open, int gap_extend) {
    int i,j,k,l;
    /* in the alignment routine all these values are added:
     * ie score is score + malign->scores[i][j];
     * even when scores[i][j] is a gap penalty.
     * to calculate the score j = seq2[row-1]
     */

    /* do the characters which are non-zero */
    for(i=0;i<malign->length;i++) {
	for(j=0;j<malign->charset_size;j++) {
	    malign->scores[i][j] *= malign->matrix[j][j];
	}
    }
    /* do the characters which are zero */
    k = malign->matrix[0][1];
    for(i=0;i<malign->length;i++) {
	l = malign->scores[i][malign->charset_size];
	for(j=0;j<malign->charset_size;j++) {
	    if (!(malign->scores[i][j])) malign->scores[i][j] = k*l;
	}
    }
    /* set the gap_open and gap_extend elements */
    j = malign->charset_size;
    for(i=0;i<malign->length;i++) {
	malign->scores[i][j]*=gap_open;
	malign->scores[i][j+1]*=gap_extend;
    }
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

MSEG **get_malign_segs(CONTIGL *contigl) {

    CONTIGL *t;
    MSEG **msegs;
    MSEG *mseg;

    int i,nseqs;

    nseqs = contigl_elements(contigl);
    msegs = (MSEG **)malloc(nseqs*sizeof(MSEG *));
    t = contigl;
    i = 0;
    while(t) {
	mseg = create_mseg();
	init_mseg(mseg,t->mseg->seq,t->mseg->length,t->mseg->offset);
	msegs[i++] = mseg;
	t = t->next;
    }
    return msegs;
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
	  space,
	  t->mseg->seq);
	*/
	t = t->next;
    }
}

MALIGN *contigl_to_malign(CONTIGL *contigl_in) {

    MALIGN *malign;
    CONTIGL *contigl;
    int i;
    int gap_open, gap_extend;
    gap_open = -12;
    gap_extend = -4;

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
    malign->nseqs = contigl_elements(contigl);

    /* get the initial scores from the alignment */

    malign->scores = create_malign_counts(malign->length,malign->charset_size+2);
    get_malign_counts(malign);


    /*print_malign_scores(malign);*/

    /* make a 100% consensus for the alignment */

    malign->consensus = (char *) malloc(malign->length);
    get_malign_consensus(malign);
    printf("      %s\n",malign->consensus);

    /* scale the scores with the gap penalties and the external score matrix */

    scale_malign_scores(malign,gap_open,gap_extend);
    print_malign_scores(malign);

    /* FIXME put in error checking for failed mallocs etc */
    return malign;
}

MOVERLAP *create_moverlap(void) {
    MOVERLAP   *moverlap;

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

int affine_malign(MOVERLAP *moverlap, ALIGN_PARAMS *params) {

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
    int best_F1, best_G1, best_H1, V_diag, V_extx, V_exty, V_insx, V_insy;
    int E_gap, F_gap;
    int edge_mode, best_edge_score;

    int band, band_length, two_band, band_left, band_right, first_band_left=0;
    int off_set, guard_offset, *pF_guard, *pG_guard, *pH_guard;
    int row, first_row, max_row, column, max_col;
    unsigned char *bit_trace;
    int byte, nibble, e_row, e_col, size_mat;
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
	E_gap = scores[0][gap_open_index];
	H1[0] = G1[0] = E_gap;
	for(i = 1; i < seq1_len; i++) {
	    F1[i] = E_gap;
	    E_gap += scores[i][gap_extend_index];
	    H1[i] = E_gap;
	    G1[i] = E_gap;
	}
	E_gap = scores[0][gap_open_index] + scores[0][gap_extend_index];
	F_gap = F1[0];
	edge_inc = 1;
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
	    gap_pen = scores[MIN(row-1,seq1_len-1)][gap_extend_index];
	    if(edge_inc) {
		E_gap += gap_pen;
		F_gap += gap_pen;
	    }
	    row_index = malign_lookup[(int)seq2[row-1]];

	    /* FIXME: got the axes switched for the single sequence methods
	     * what is the correct thing to do here????
	     */

	    if ( seq2[row-1] != OLD_PAD_SYM ) {
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

		/* move diagonally? */
		s = scores[column-1][row_index];
		gap_open_score = scores[column-1][gap_open_x];
		gap_extend_score = scores[column-1][gap_extend_x];

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
		V_diag =  *pF2 + gap_open_score;
		V_extx =  *pH2 + gap_extend_score;
		if ( V_diag > V_extx ) {
		    best_H1 = V_diag;
		}
		else {
		    best_H1 = V_extx;
		}
		*(pH2+1) = best_H1;

		/* gap in y? */
		V_diag = *(pF1+1) + gap_open_score;
		V_exty = *(pG1+1) + gap_extend_score;
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
     
	    gap_pen = scores[MIN(row-1,seq1_len-1)][gap_extend_index];
	    if(edge_inc) {
		E_gap += gap_pen;
		F_gap += gap_pen;
	    }

	    row_index = malign_lookup[(int)seq2[row-1]];

     
	    if ( seq2[row-1] != OLD_PAD_SYM ) {
		gap_open_x = gap_open_index;
		gap_extend_x = gap_extend_index;
	    }
	    else {
		gap_open_x = gap_extend_x = gap_match_index;
	    }


	    /* process each column. i.e. each character of seq1 */
     
	    for(column = 1; column <= seq1_len; column++, e++) {
		/* move diagonally? */
		s = scores[column-1][row_index];
		gap_open_score = scores[column-1][gap_open_x];
		gap_extend_score = scores[column-1][gap_extend_x];

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

		/* gap in x? */
		V_diag =  pF2[column-1] + gap_open_score;
		V_extx =  pH2[column-1] + gap_extend_score;
		if ( V_diag > V_extx ) {
		    best_H1 = V_diag;
		}
		else {
		    best_H1 = V_extx;
		}
		pH2[column] = best_H1;
       
		/* gap in y? */
		V_diag =  pF1[column] + gap_open_score;
		V_exty = pG1[column] + gap_extend_score;
		if ( V_diag > V_exty ) {
		    best_G1 = V_diag;
		}
		else {
		    best_G1 = V_exty;
		}
		pG2[column] = best_G1;

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
		/*
		  if (( best_H1 > best_F1 ) && ( best_H1 > best_G1 )) {
		  bit_trace[byte] |= BYTE_ACROSS << nibble;
		  b_s = best_H1;
		  }
		  else if (( best_G1 > best_F1 ) && ( best_G1 > best_H1 )) {
		  bit_trace[byte] |= BYTE_DOWN << nibble;
		  b_s = best_G1;
		  }
		  else {
		  bit_trace[byte] |= BYTE_DIAGONAL << nibble;
		  b_s = best_F1;
		  }
		*/
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
	destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	return -1;
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_EDIT_BUFFERS ) {
	if (seq_to_edit ( seq1_out,seq_out_len,&moverlap->S1,&moverlap->s1_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
	if (seq_to_edit ( seq2_out,seq_out_len,&moverlap->S2,&moverlap->s2_len,NEW_PAD_SYM)) {
	    destroy_af_mem ( F1, F2, G1, G2, H1, H2, bit_trace, seq1_out, seq2_out );
	    return -1;
	}
    }

    if ( params->return_job & SP_ALIGNMENT_RETURN_SEQ ) {
	if ( !(params->return_job & SP_ALIGNMENT_RETURN_NEW_PADS) ) {
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

} /* namespace */
