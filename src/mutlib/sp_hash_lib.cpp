#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <staden.h>
#include <sp_alignment_structs.h>
#include <sp_hash_lib.h>


namespace sp {


#define MINMAT 20


int dna_hash8_lookup[256];

/**
  * set the lookup values for hashing
  */

void set_hash8_lookupn(void) {

    int i;

    for (i=0;i<256;i++) dna_hash8_lookup[i] = 4;

    dna_hash8_lookup['a'] = 0;
    dna_hash8_lookup['c'] = 1;
    dna_hash8_lookup['g'] = 2;
    dna_hash8_lookup['t'] = 3;
    dna_hash8_lookup['A'] = 0;
    dna_hash8_lookup['C'] = 1;
    dna_hash8_lookup['G'] = 2;
    dna_hash8_lookup['T'] = 3;
    dna_hash8_lookup['*'] = 0;
}

void print_h (Hash *h) {
    printf("word_length %d size_hash %d seq1_len %d seq2_len %d\n",
           h->word_length,h->size_hash,h->seq1_len,h->seq2_len);
}

    /** given a sequence seq, return the hash value for the first word 
     *  after start_base that does not contain an unknown char. Tell 
     *  the caller where this is. If we reach the end of the seq set
     *  start_base and return -1.
     */
    
int hash_word8n ( char *seq, int *start_base, int seq_len, int word_length,
                unsigned short *uword) {
    
    register int i;
    register int end_base,base_index,lstart_base;
    register int unsigned short luword;
    
    lstart_base = *start_base;
    end_base = lstart_base + word_length;
    if ( seq_len < end_base ) return -1;
    
    for (i=lstart_base,luword=0,end_base=lstart_base+word_length;i<end_base;i++) {

        base_index = dna_hash8_lookup[(unsigned)seq[i]];
        if ( 4 == base_index ) {

            /*  weve hit an unknown char, so lets start again */

            lstart_base = i + 1;
            end_base = lstart_base + word_length;
            if ( seq_len < end_base ) {
                *start_base = lstart_base;
                return -1;
            }
            luword = 0;
            i = lstart_base - 1;
        }
        else {
            luword = ( luword <<2 ) | base_index;
        }
    }
    *start_base = lstart_base;
    *uword = luword;
    return 0;
}


    /** given a sequence seq, return an array of hash values
      * If we cannot find at least one word to hash on we return -1
      * otherwise we return 0.
      */

int hash_seq8n ( char *seq, int *hash_values, int seq_len, int word_length) {


    register int i,j,k;
    int start_base,prev_start_base,end_base,base_index;
    unsigned short uword;

    if ( seq_len < word_length ) return -1;

    /*  Get the hash value for the first word that contains no unknowns */      
    start_base = 0;
    if (hash_word8n ( seq, &start_base, seq_len, word_length, &uword)) return -1;

    for (i=0;i<start_base;i++) hash_values[i] = -1;

    /*  Now do the rest of the sequence */

    hash_values[start_base] = uword;
    k = seq_len - word_length + 1;

    for (i=start_base+1,j=start_base+word_length; i<k; i++,j++) {

        base_index = dna_hash8_lookup[(unsigned)seq[j]];
        if ( 4 == base_index ) {

            /*  weve hit an unknown char, so lets start again */

            prev_start_base = i;
            start_base = j + 1;
            if (hash_word8n ( seq, &start_base, seq_len, word_length, &uword)) {
                for (i=prev_start_base;i<start_base;i++) hash_values[i] = -1;
                return 0;
            }

            for (i=prev_start_base;i<start_base;i++) hash_values[i] = -1;
            hash_values[start_base] = uword;
            end_base = start_base + word_length;
            i = start_base;
            j = i + word_length - 1;

        }
        else {
            uword = ( uword <<2 ) | base_index;
            hash_values[i] = uword;
        }
    }
    return 0;
}


    /** given a sequence seq, return the hash value for the first word 
     *  after start_base that does not contain an unknown char. Tell 
     *  the caller where this is. If we reach the end of the seq set
     *  start_base and return -1.
     */
    
int hash_word4n ( char *seq, int *start_base, int seq_len, int word_length,
                unsigned char *uword) {
    
    
    register int i;
    register int end_base,base_index,lstart_base;
    register char unsigned luword;
    
    lstart_base = *start_base;
    end_base = lstart_base + word_length;
    if ( seq_len < end_base ) return -1;
    
    for (i=lstart_base,luword=0,end_base=lstart_base+word_length;i<end_base;i++) {

        base_index = dna_hash8_lookup[(unsigned)seq[i]];
        if ( 4 == base_index ) {

            /*  weve hit an unknown char, so lets start again */

            lstart_base = i + 1;
            end_base = lstart_base + word_length;
            if ( seq_len < end_base ) {
                *start_base = lstart_base;
                return -1;
            }
            luword = 0;
            i = lstart_base - 1;
        }
        else {
            luword = ( luword <<2 ) | base_index;
        }
    }
    *start_base = lstart_base;
    *uword = luword;
    return 0;
}


    /** given a sequence seq, return an array of hash values
      * If we cannot find at least one word to hash on we return -1
      * otherwise we return 0.
      */

int hash_seq4n ( char *seq, int *hash_values, int seq_len, int word_length) {

    register int i,j,k;
    int start_base,prev_start_base,end_base,base_index;
    unsigned char uword;

    if ( seq_len < word_length ) return -1;

    /*  Get the hash value for the first word that contains no unknowns */      
    start_base = 0;
    if (hash_word4n ( seq, &start_base, seq_len, word_length, &uword)) return -1;

    for (i=0;i<start_base;i++) hash_values[i] = -1;

    /*  Now do the rest of the sequence */

    hash_values[start_base] = uword;
    k = seq_len - word_length + 1;

    for (i=start_base+1,j=start_base+word_length; i<k; i++,j++) {

        base_index = dna_hash8_lookup[(unsigned)seq[j]];
        if ( 4 == base_index ) {

            /*  weve hit an unknown char, so lets start again */

            prev_start_base = i;
            start_base = j + 1;
            if (hash_word4n ( seq, &start_base, seq_len, word_length, &uword)) {
                for (i=prev_start_base;i<start_base;i++) hash_values[i] = -1;
                return 0;
            }

            for (i=prev_start_base;i<start_base;i++) hash_values[i] = -1;
            hash_values[start_base] = uword;
            end_base = start_base + word_length;
            i = start_base;
            j = i + word_length - 1;

        }
        else {
            uword = ( uword <<2 ) | base_index;
            hash_values[i] = uword;
        }
    }
    return 0;
}

    /** return the length of a diagonal given the diagonal number
      *
      * 0123456789 seq1->
      * 1
      * 2
      * 3        diagonal numbers
      * 4       ^
      * 5     /
      * 6    2
      * 7  1
      * 0
      * s
      * e
      * q
      * 2
      * |
      * v
      */

int diagonal_length(int seq1_len, int seq2_len, int diagonal_number) {

    if(diagonal_number < seq1_len)
        return MIN((diagonal_number + 1), (
               MIN(seq1_len, seq2_len)));
    else
        return MIN((seq1_len + seq2_len - 1 - diagonal_number), (
               MIN(seq1_len, seq2_len)));
}

    /** given a diagonal number, return the intercepts on the 
      *seq1 and seq2 axes
      */

void diagonal_intercepts (int diagonal_number, int seq1_len, int seq2_len,
                          int *seq1_intercept, int *seq2_intercept ) {

    if ( diagonal_number < seq1_len ) {
        *seq1_intercept = seq1_len - diagonal_number - 1;
        *seq2_intercept = 0;
    }
    else {
        *seq2_intercept = diagonal_number + 1 - seq1_len;
        *seq1_intercept = 0;
    }
    
}


#define SMALL_POLY 1.0e-30
#define ZERO_SMALL(a) ( ((a) < (SMALL_POLY)) ?  (0.0) : (a) )

/**
  * multiplication
  */

int poly_mult (Poly *poly) {
    int i,j,max_terms;

    max_terms = poly->num_terms + poly->size_step;
    if ( max_terms > MAX_POLY ) return -1;
    for(i=0;i<=max_terms;i++) poly->c[i] = 0.0;
    for (i=0;i<=poly->num_terms;i++) {
        for (j=0;j<=poly->size_step;j++) {
            poly->c[i+j] += poly->a[i] * poly->b[j];
        }
    }
    poly->num_terms = max_terms;
    for(i = 0; i <= poly->num_terms; i++) 
        poly->a[i] = ZERO_SMALL ( poly->c[i] );
    return 0;
}

/**
  * calculate the probability for a word of length word_length
  */

double prob_word ( int word_length, double comp[]) {
    Poly p, *poly;
    int i,j,k;

    poly = &p;
    poly->rows = poly->cols = 4;
    poly->size_step = poly->num_terms = 1;
    for (i=0;i<MAX_POLY;i++) poly->a[i] = poly->b[i] = 0.0;
    for (i=0;i<poly->rows;i++) {
        for ( j=0;j<poly->cols;j++) {
            k = (i==j) ? 1:0;
            poly->b[k] = poly->a[k] += comp[i] * comp[j];
        }
    }
    for (i = 1;i < word_length;i++) {
        j = poly_mult(poly);
        if ( j ) return -1.0;
    }
    /*for(i=poly->num_terms;i>-1;i--) poly->a[i] += poly->a[i+1];*/
    return poly->a[word_length];
}

    /** routine to find the best intercept given a set of matches
     * between a pair of sequences. This is done iteratively:
     * we find the centre of gravity, remove the outlier, and
     * repeat until the best match only remains.
     * NOTE: THIS ROUTINE ZEROES THE MATCH ARRAY
     */
int best_intercept ( Hash *h, int *seq1_i, int *seq2_i ) {


    double t, sum_scores, sum_moment, c_o_g, furthest;
    int match_no, matches_left, outlier;
    sum_scores = sum_moment = 0.0;

    for(matches_left = h->matches; matches_left > 1; matches_left--) {
        sum_moment = sum_scores = 0.0;
        for(match_no = 0; match_no < h->matches; match_no++) {
            if ( h->diag_match[match_no].prob > 0.0 ) {
                sum_moment += h->diag_match[match_no].pos *
                    h->diag_match[match_no].prob;
                sum_scores += h->diag_match[match_no].prob;
            }
        }
        if ( sum_scores == 0 ) {
            fprintf(stderr, "FATAL: best_intecept has sum_scores of 0\n");
            return 0;
        }
        c_o_g = sum_moment/sum_scores;
        outlier = -1;
        for(match_no = 0, furthest = 0.0; match_no < h->matches; match_no++) {
            if ( h->diag_match[match_no].prob > 0.0 ) {
                if ((t=fabs(c_o_g - h->diag_match[match_no].pos)) > furthest) {
                    outlier = match_no;
                    furthest = t;
                }
            }
        }
        if(outlier == -1) {
            /* must have matches on same diagonal whcih is cog */
            for(match_no = 0, furthest = -1.0; match_no < h->matches; match_no++) {
                if ( h->diag_match[match_no].prob > 0.0 ) {
                    if ((t=fabs(c_o_g - h->diag_match[match_no].pos)) > furthest) {
                        outlier = match_no;
                        furthest = t;
                    }
                }
            }
        }
        h->diag_match[outlier].prob = 0.0;
    }
    for(match_no = 0; match_no < h->matches; match_no++) {
        if ( h->diag_match[match_no].prob > 0.0 ) {
            diagonal_intercepts (h->diag_match[match_no].pos, h->seq1_len, 
                                 h->seq2_len, seq1_i, seq2_i);
            break;
        }
    }
    return 1;
    }

    /** routine to remove duplicates from a list of repeats. It
      * also removes the self match.
      * Input: a list of *n_match match positions and match
      * lengths.  Output: a list in which all duplicates and the selfmatch
      * are removed.  *n_match is set to the new number of matches, or -1
      * for error.
      */
    
void remdup ( int *seq1_match, int *seq2_match, int *len_match, int
             *n_matches ) {
    
    register int i;
    int *index_ptr,k,keep;
    
    if ( *n_matches < 1 ) return;

    if ( ! ( index_ptr = (int *) xmalloc ( sizeof(int)*(*n_matches) ))) {
        *n_matches = -1;
        return;
    }
    for ( i=0,k=0;i<*n_matches; i++) {
        if ( seq1_match[i] > seq2_match[i] ) {
            index_ptr[k] = i; 
            k += 1;
        }
    }
    for ( i=0; i<k; i++) {
        keep = index_ptr[i];
        seq1_match[i] = seq1_match[keep];
        seq2_match[i] = seq2_match[keep];
        len_match[i] = len_match[keep];
    } 
    *n_matches = k;
    if ( index_ptr ) free ( index_ptr );
}

/**
  * set the match positions from a complemented sequence
  * to be relative to the original orientation
  */

void make_reverse ( int *seq2_match, int *len_match,
                   int n_matches, int seq2_len ) {
    
    int i;
    
    for (i = 0; i< n_matches; i++) {
        seq2_match[i] = seq2_len - seq2_match[i] - len_match[i]+2;
        
    }
}
/**
  * add edit pair data from an alignment of a segment (stored in overlap)
  * to the end of an existing edit buffer edit_pair
  */

int update_edit_pair ( EDIT_PAIR *edit_pair, OVERLAP *overlap ) {
    int i,j;
    /*printf("s1 %d s2 %d\n",overlap->s1_len,overlap->s2_len);*/
    if ( overlap->s1_len ) {
        if ( (edit_pair->size - edit_pair->next1) < overlap->s1_len ) return -1;
        for (i=edit_pair->next1,j=0;j<overlap->s1_len;i++,j++) {
            edit_pair->S1[i] = overlap->S1[j];
        }
        edit_pair->next1 += overlap->s1_len;
        xfree ( overlap->S1 );
        overlap->S1 = NULL;
        overlap->s1_len = 0;
    }
    if ( overlap->s2_len ) {
        if ( (edit_pair->size - edit_pair->next2) < overlap->s2_len ) return -1;
        for (i=edit_pair->next2,j=0;j<overlap->s2_len;i++,j++) {
            edit_pair->S2[i] = overlap->S2[j];
        }
        edit_pair->next2 += overlap->s2_len;
        xfree ( overlap->S2 );
        overlap->S2 = NULL;
        overlap->s2_len = 0;
    }
    return 0;
}

/**
  * add a block of alignment (ie a single perfect match) to
  * an existing edit_pair buffer
  */

int block_to_edit_pair ( EDIT_PAIR *edit_pair, int length ) {

  if ( (edit_pair->size - edit_pair->next1) < 1 ) return -1;
  edit_pair->S1[edit_pair->next1++] = length;
  if ( (edit_pair->size - edit_pair->next2) < 1 ) return -1;
  edit_pair->S2[edit_pair->next2++] = length;
  return 0;
}

/**
  * sent a pair of sequence segments to align (in overlap)
  * if they are both non-zero uses affine_align to align them
  * then sticks their edit buffers onto a growing edit_buffer
  * else if one is zero length it updates the edit_pair buffer accordingly
  */

int align_bit ( ALIGN_PARAMS *params, OVERLAP *overlap, EDIT_PAIR *edit_pair) {

    int l1, l2, ret;

    l1 = overlap->seq1_len;
    l2 = overlap->seq2_len;
    /*
    printf("align_bit l1 l2 %d %d\n",l1,l2);
    */
    if ((l1 > 0) && (l2 > 0 )) {
        if(ret = affine_align(overlap,params)) return -1;
        if ( update_edit_pair ( edit_pair, overlap)) return -1;
    }
    else {
        if (l1 > 0 ) {
            if ( edit_pair->next2 == edit_pair->size ) return -1;
            edit_pair->S2[edit_pair->next2++] = -l1;

            if ( edit_pair->next1 == edit_pair->size ) return -1;
            edit_pair->S1[edit_pair->next1++] = l1;
        }
        if (l2 > 0 ) {
            if ( edit_pair->next1 == edit_pair->size ) return -1;
            edit_pair->S1[edit_pair->next1++] = -l2;

            if ( edit_pair->next2 == edit_pair->size ) return -1;
            edit_pair->S2[edit_pair->next2++] = l2;
        }
    }
    return 0;
}

/**
  * destroy an edit_pair
  */

void destroy_edit_pair (EDIT_PAIR *edit_pair) {
    if (edit_pair) {
        if (edit_pair->S1) xfree (edit_pair->S1);
        if (edit_pair->S2) xfree (edit_pair->S2);
        xfree (edit_pair);
    }
}

/**
  * create and return an edit_pair of size size
  */

EDIT_PAIR *create_edit_pair(int size) {
    EDIT_PAIR *edit_pair;

    if(NULL == (edit_pair = (EDIT_PAIR *) xmalloc(sizeof(EDIT_PAIR)))) {
        verror(ERR_WARN, "create_edit_pair", "xmalloc failed");
        return NULL;
    }

    if ( ! ((edit_pair->S1 = (int *) xmalloc ( sizeof(int)*(size) )))) {
        destroy_edit_pair (edit_pair);
        verror(ERR_WARN, "create_edit_pair", "xmalloc failed");
        return NULL;
    }
    if ( ! ((edit_pair->S2 = (int *) xmalloc ( sizeof(int)*(size) )))) {
        destroy_edit_pair (edit_pair);
        verror(ERR_WARN, "create_edit_pair", "xmalloc failed");
        return NULL;
    }
    edit_pair->next1 = 0;
    edit_pair->next2 = 0;
    edit_pair->size = size;
    return edit_pair;
}

/**
  * returns a band_size dependent on the length of the sequences to be aligned
  * max(30,35%(min))
  */

int set_band_blocks(int seq1_len, int seq2_len) {
    int band;
/*    return MAX(30,(MIN(seq1_len,seq2_len)*0.35));*/
    /* want no more that seq2_len but up to 35% of seq1_len */
    band = (int)MIN(((seq1_len+1)/2),(seq2_len*0.35));
    printf("seq1_len %d seq2_len %d band %d\n", seq1_len,seq2_len,band);
    return band;
}


    /** fiddle left end of edit buffer for cases when the sequences aligned
     *  are segments which align up to the start point or where one is
     *  a segment and the other starts somewhere to the right
     *
     *  if seq1_start == 0, seq2_start != 0 we must add an extra edit pair
     *  if both are != 0 and != we add 2 edit pairs to the shorter, 1 to
     *  the longer, if == add 1 pair, if both 0 add nothing
     */

void left_edit_buffer (OVERLAP *overlap, ALIGN_PARAMS *params, 
                       int *S1_next, int *S2_next) {

    int j1, j2;

    j1 = j2 = 0;

    if ((params->seq1_start > 0) && (params->seq2_start > 0)) {
        if (params->seq1_start > params->seq2_start) {
            overlap->S1[0] = params->seq1_start;
            overlap->S2[0] = params->seq2_start - params->seq1_start;
            overlap->S2[1] = params->seq2_start;
            j1 = 1;
            j2 = 2;
        }
        else if (params->seq2_start > params->seq1_start) {
            overlap->S2[0] = params->seq2_start;
            overlap->S1[0] = params->seq1_start - params->seq2_start;
            overlap->S1[1] = params->seq1_start;
            j1 = 2;
            j2 = 1;
        }
        else if (params->seq1_start == params->seq2_start) {
            overlap->S1[0] = params->seq1_start;
            overlap->S2[0] = params->seq2_start;
            j1 = j2 = 1;
        }
    }
    else {
        if (params->seq1_start > 0) {
            overlap->S1[0] = params->seq1_start;
            overlap->S2[0] = -params->seq1_start;
            j1 = j2 = 1;
        }
        if (params->seq2_start > 0) {
            overlap->S2[0] = params->seq2_start;
            overlap->S1[0] = -params->seq2_start;
            j1 = j2 = 1;
        }
    }
    /*
    printf("left_edit %d %d\n",j1,j2);
    */
    *S1_next = j1;
    *S2_next = j2;
}

    /** fiddle right end of edit buffer for cases when the sequences aligned
     *  are segments which align after the end point or where one is
     *  a segment and the other ends somewhere to the left
     *  as we assume they align at the ends of the segments, we just
     *  have to add on the ends, then pad out the shorter one (if necessary)
     *
     *  if seq1_end == 0, seq2_end != 0 we must add an extra edit pair
     *  if both != 0 we must add 2 edit pairs to the shorter, 1 to the longer
     *  if == add 1 pair if both 0 do nothing
     */
void right_edit_buffer (OVERLAP *overlap, ALIGN_PARAMS *params, 
                        int *S1_next, int *S2_next) {

    int i, i1, i2, j1, j2, k1, k2;

    
    if ((params->seq1_end == 0) && (params->seq2_end == 0)) return;
    i1 = *S1_next;
    i2 = *S2_next;
    
    
    for(i=0,j1=0;i<i1;i++) {
        if(overlap->S1[i]>0)j1+=overlap->S1[i];
    }
    
    
    for(i=0,j2=0;i<i2;i++) {
        if(overlap->S2[i]>0)j2+=overlap->S2[i];
    }
    
    
    k1 = overlap->seq1_len - params->seq1_end - 1;
    k2 = overlap->seq2_len - params->seq2_end - 1;
    /*
    printf("start 1 %d start 2 %d\n",params->seq1_start,params->seq2_start);
    printf("end 1 %d end 2 %d\n",params->seq1_end,params->seq2_end);
    printf("len %d len %d\n",overlap->seq1_len,overlap->seq2_len);
    printf("k1 %d k2 %d\n",k1,k2);
    */

    if ((params->seq1_end < overlap->seq1_len-1) && 
        (params->seq2_end < overlap->seq2_len-1)) {
        if (params->seq1_end > params->seq2_end) {
            overlap->S1[i1++] = k1;
            overlap->S2[i2++] = k2;
            overlap->S2[i2++] = k2 - k1;
        }
        else if (params->seq2_end > params->seq1_end) {
            overlap->S2[i2++] = k2;
            overlap->S1[i1++] = k1;
            overlap->S1[i1++] = k1 - k2;
        }
        else if (params->seq1_end == params->seq2_end) {
            overlap->S1[i1++] = k1;
            overlap->S2[i2++] = k2;
        }
    }
    else {
        if (params->seq1_end < overlap->seq1_len - 1) {
            overlap->S1[i1++] = k1;
            overlap->S2[i2++] = -k1;
        }
        if (params->seq2_end < overlap->seq2_len - 1) {
            overlap->S2[i2++] = k2;
            overlap->S1[i1++] = -k2;
        }
    }
    
    *S1_next = i1;
    *S2_next = i2;
}

    /** input h which contains a list of matching blocks
     * output overlap_out which contains an overlap structure
     * ie a pair of edit buffers and a pair of aligned seqs
     *
     * strategy is to build up a pair of overall edit buffers
     * using the blocks and affine_align where necessary in between.
     * in some gaps between blocks (where 1 gap is length 0) we
     * simply add to the edit buffers. The edit buffers are stored
     * in an edit_pair structure and their contents copied to
     * overlap_out at the end. After each call to affine_align
     * we copy its edit buffers to the edit_pair structure.
     * 3 phases of block alignment:
     * 1. up to first block (if there is a mismatched region)
     * 2. between blocks
     * 3. from last block to end (if there is a mismatched region)
     *
     * NB the algorithm can also handle the alignment of segments
     * and then add the surrounding bits of sequence afterwards by
     * use of the edit buffers and the full length sequences stored
     * in the overlap structure.
     * Explanation:
     * 1. the seq1_start and seq2_start values in ALIGN_PARAMS are set
     * 2. if seq1_end!=0 || seq2_end!=0 values in ALIGN_PARAMS are set
     * 3. h->seq1_len = seq1_end - seq1_start + 1
     * 4. h->seq1 = overlap->seq1 + seq1_start
     * 5. ditto seq2
     * 6. left_edit_buffer and right_edit_buffer are used to update
     *    the edit buffers at each end with values that add in the 
     *    external sequence
     * 7. the edit buffers and the overlap->seq arrays are used to
     *    create the complete alignment in overlap->seq_out
     * 8. NOTE that this only works if either one sequence lies
     *    entirely within the other (ie only the longer sequence is
     *    a segment), or the regions outside the segment being aligned
     *    align exactly at the boundaries (this is assumed by the routines
     *    that add on the sequences outside the segments).
     *
     * having built up the edit buffers, at the end we use them to
     * create a sequence alignment in order to fill out the other
     * bits of the overlap structure (would have prefered to be
     * able to avoid this and optionally only return the edit buffers,
     * but that requires extra work)
     * The alternative, which now looks more attractive, given that
     * we plan to use aligned sequences to direct editing, is to
     * rewrite all this stuff to create a pair of sequence alignments
     * as we go along, and not use edit buffers at all!!
     *
     * functions:
     * align_bit
     * sent a segment; either updates edit_pairs directly or uses
     * affine_align and then copies its edit buffers into edit_pair.
     * update_edit_pair
     * copies the edit buffers from an overlap into the edit_pair
     * frees the overlap edit buffers
     * block_to_edit_pair
     * adds a matching block to the edit_pair
     *
     * essential to make sure that the affine_align routine uses
     * the correct edge_mode
     * At the end we destroy all temporary memory: overlap and edit_pair
     */

int align_wrap ( Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap_out) {


    int i, j, s1, s2;
    int S1_next, S2_next;
    int len_seq;
    int band, band_in;
    OVERLAP *overlap;
    EDIT_PAIR *edit_pair;
    int max_edit_pair;
    int max_seq;
    char NEW_PAD_SYM, OLD_PAD_SYM;

    NEW_PAD_SYM = params->new_pad_sym;
    OLD_PAD_SYM = params->old_pad_sym;

    band = 0;
    band_in = params->band;
    max_edit_pair = MIN(h->seq1_len,h->seq2_len);

    if (NULL == (edit_pair = create_edit_pair(max_edit_pair))) {
        return -1;
    }
    if (NULL == (overlap = create_overlap())) {
        destroy_edit_pair(edit_pair);
        return -1;
    }

    init_overlap (overlap, h->seq1, h->seq2, h->seq1_len, h->seq2_len);
 
    /* align up to the first matching words,
     * align the segments between the matching words
     * align from the last matching words to the ends
     */

    /* get a start point based on the first matching word positions */

    diagonal_intercepts (h->block_match[0].diag,h->seq1_len,h->seq2_len,
                         &s1, &s2);

    overlap->seq1_len = h->block_match[0].pos_seq1;
    overlap->seq2_len = h->block_match[0].pos_seq2;
    overlap->seq1 = h->seq1;
    overlap->seq2 = h->seq2;
    len_seq = MAX(overlap->seq1_len,overlap->seq2_len);
    
    /*
    printf("intercepts s1 %d s2 %d len_to_align %d\n",s1,s2,len_seq);
    printf("block 0 pos1 %d pos2 %d len %d\n",h->block_match[0].pos_seq1,
           h->block_match[0].pos_seq2,h->block_match[0].length);
    */
    params->edge_mode = 6;
    if ( band_in) band = set_band_blocks(overlap->seq1_len,overlap->seq2_len);
    set_align_params_banding (params, band, s1, s2);

    if (align_bit ( params, overlap, edit_pair)) {
        verror(ERR_WARN, "align_wrap", "failed in align_bit");
        destroy_edit_pair(edit_pair);
        destroy_overlap(overlap);
        return -1;
    }
    free_overlap(overlap);

    /*printf("if required, done up to block 0\n");*/

    if ( block_to_edit_pair ( edit_pair, h->block_match[0].length)) {
        verror(ERR_WARN, "align_wrap", "failed in block_to_edit_pair");
        destroy_edit_pair(edit_pair);
        destroy_overlap(overlap);
        return -1;
    }
    s1 = h->block_match[0].pos_seq1 + h->block_match[0].length;
    s2 = h->block_match[0].pos_seq2 + h->block_match[0].length;
    params->edge_mode = 5;
    for(i=1;i<h->matches;i++) {
        overlap->seq1_len = h->block_match[i].pos_seq1 - s1;
        overlap->seq2_len = h->block_match[i].pos_seq2 - s2;
        overlap->seq1 = &(h->seq1[s1]);
        overlap->seq2 = &(h->seq2[s2]);
        len_seq = MAX(overlap->seq1_len,overlap->seq2_len);

        /*
        printf("s1 %d s2 %d len_to_align %d\n",s1,s2,len_seq);
        printf("block %d pos1 %d pos2 %d len %d\n",i,h->block_match[i].pos_seq1,
               h->block_match[i].pos_seq2,h->block_match[i].length);
        */
        if ( len_seq > 0 ) {

            if(band_in)band = set_band_blocks(overlap->seq1_len,overlap->seq2_len);
            set_align_params_banding (params, band, 0, 0);
            if (align_bit ( params, overlap, edit_pair)) {
                verror(ERR_WARN, "align_wrap", "failed in align_bit");
                destroy_edit_pair(edit_pair);
                destroy_overlap(overlap);
                return -1;
            }
            free_overlap(overlap);

        }
        s1 = h->block_match[i].pos_seq1 + h->block_match[i].length;
        s2 = h->block_match[i].pos_seq2 + h->block_match[i].length;
        if ( block_to_edit_pair ( edit_pair, h->block_match[i].length)) {
            verror(ERR_WARN, "align_wrap", "failed in block_to_edit_pair");
            destroy_edit_pair(edit_pair);
            destroy_overlap(overlap);
            return -1;
        }
    }
    /* do up to the end */
    
    overlap->seq1_len = h->seq1_len - s1;
    overlap->seq2_len = h->seq2_len - s2;
    /*printf("at end s1 %d s2 %d len1 %d len2 %d\n",s1,s2,overlap->seq1_len,overlap->seq2_len);*/
    overlap->seq1 = &(h->seq1[s1]);
    overlap->seq2 = &(h->seq2[s2]);
    
    if(band_in)band = set_band_blocks(overlap->seq1_len,overlap->seq2_len);
    set_align_params_banding (params, band, 0, 0);
    params->edge_mode = 9;
    if (align_bit ( params, overlap, edit_pair)) {
        verror(ERR_WARN, "align_wrap", "failed in align_bit");
        destroy_edit_pair(edit_pair);
        destroy_overlap(overlap);
        return -1;
    }
    destroy_overlap(overlap);
    /*
    for(i=0;i<edit_pair->next1;i++) printf("1 %d %d\n",i,edit_pair->S1[i]);
    for(i=0;i<edit_pair->next2;i++) printf("2 %d %d\n",i,edit_pair->S2[i]);
    
    
    i= print_alignment(h->seq1,
                    h->seq2,
                    h->seq1_len,
                    h->seq2_len,
                    edit_pair->S1,
                    edit_pair->S2,
                    edit_pair->next1,
                    edit_pair->next2,
                    100,
                    stdout);
    */
    /* now create the overlap from the edit buffers */

    max_seq = overlap_out->seq1_len + overlap_out->seq2_len + 1;

    if(!(overlap_out->seq1_out = (char *) xmalloc(sizeof(char) * 
                                                      max_seq))) {
        verror(ERR_WARN, "align_wrap", "malloc failed for seq1_out");
        destroy_edit_pair(edit_pair);
        return -1;
    }
    if(!(overlap_out->seq2_out = (char *) xmalloc(sizeof(char) * 
                                                      max_seq))) {
        verror(ERR_WARN, "align_wrap", "malloc failed for seq2_out");
        destroy_edit_pair(edit_pair);
        return -1;
    }

    shrink_edit_buffer(edit_pair->S1,&edit_pair->next1);
    shrink_edit_buffer(edit_pair->S2,&edit_pair->next2);

    /*
    for(i=0;i<edit_pair->next1;i++) printf("1 %d %d\n",i,edit_pair->S1[i]);
    for(i=0;i<edit_pair->next2;i++) printf("2 %d %d\n",i,edit_pair->S2[i]);
    */

    /* save edit_pairs in overlap_out */

    if(!(overlap_out->S1 = (int *) xmalloc(sizeof(int) * 
                                                      edit_pair->next1+4))) {
        verror(ERR_WARN, "align_wrap", "malloc failed for S1");
        destroy_edit_pair(edit_pair);
        return -1;
    }
    if(!(overlap_out->S2 = (int *) xmalloc(sizeof(int) * 
                                                      edit_pair->next2+4))) {
        verror(ERR_WARN, "align_wrap", "malloc failed for S2");
        destroy_edit_pair(edit_pair);
        return -1;
    }

    /* need to deal with case where we are aligning a subsection
     *
     * note that above we have allocated enough space for 4 extra elements
     */

    left_edit_buffer(overlap_out,params,&S1_next,&S2_next);

    j = S1_next;
    for (i=0;i<edit_pair->next1;i++,j++) overlap_out->S1[j] = edit_pair->S1[i];
    S1_next = j;
    overlap_out->s1_len = S1_next;
    j = S2_next;
    for (i=0;i<edit_pair->next2;i++,j++) overlap_out->S2[j] = edit_pair->S2[i];
    S2_next = j;
    overlap_out->s2_len = S2_next;
    /*
    print_edit_buffers(overlap_out);
    */
    right_edit_buffer(overlap_out,params,&S1_next,&S2_next);
    
    overlap_out->s1_len = S1_next;
    overlap_out->s2_len = S2_next;
    destroy_edit_pair(edit_pair);
    /*
    print_edit_buffers(overlap_out);
    */
    shrink_edit_buffers(overlap_out);

    seq_expand(overlap_out->seq1, overlap_out->seq1_out, &s1, overlap_out->S1, 
               overlap_out->s1_len, 3, NEW_PAD_SYM);
    seq_expand(overlap_out->seq2, overlap_out->seq2_out, &s2, overlap_out->S2, 
               overlap_out->s2_len, 3, NEW_PAD_SYM);

    overlap_out->seq_out_len = s1;
    /* when we enter seq_to_overlap from here the overlap score is not set */
    overlap_out->score = 0;
    if( i = seq_to_overlap(overlap_out,OLD_PAD_SYM,NEW_PAD_SYM)) {
        return -1;
    }
    if ( params->return_job & RETURN_NEW_PADS ) {
        old_pads_for_new(overlap_out->seq1_out,overlap_out->seq_out_len,
                         OLD_PAD_SYM,NEW_PAD_SYM);
        old_pads_for_new(overlap_out->seq2_out,overlap_out->seq_out_len,
                         OLD_PAD_SYM,NEW_PAD_SYM);
    }
    i = overlap_score(overlap_out,params->score_matrix);
    overlap_out->score = overlap_out->qual = i;
    return 0;
}

int central_diagonal ( Hash *h ) {
    int i, j;
    if(h->matches == 0) return 0;
    for(i=0,j=0;i<h->matches;i++) j+= h->block_match[i].diag;
    return j/h->matches;
}

/**
  * comparison function for sorting aligned blocks on position
  */

extern "C" int sort_func(const void *p1, const void *p2) {
    int x1,y1,x2,y2;
    Block_Match *c1 = (Block_Match *)p1;
    Block_Match *c2 = (Block_Match *)p2;

    x1 = c1->pos_seq1;
    y1 = c1->pos_seq2;
    x2 = c2->pos_seq1;
    y2 = c2->pos_seq2;
    return (x1+y1) - (x2+y2);
}

/**
  * sort block matches on distance from edge of sequences
  */

int sort_blocks ( Block_Match *block_match, int matches ) {
    qsort ((void *) block_match, matches, sizeof(Block_Match), sort_func);
    return 0;
}

/**
  * comparison function for sorting aligned blocks on length
  */
extern "C" int sort_len_func(const void *p1, const void *p2) {
    int x1,x2;
    Block_Match *c1 = (Block_Match *)p1;
    Block_Match *c2 = (Block_Match *)p2;

    x1 = c1->length;
    x2 = c2->length;
    return (x2-x1);
}

/**
  * sort block matches on length
  */
int sort_len_blocks ( Block_Match *block_match, int matches ) {
    qsort ((void *) block_match, matches, sizeof(Block_Match), sort_len_func);
    return 0;
}

/**
  * from an input set of matching blocks (in h) produce a sequence alignment
  * output in overlap
  * algorithm:
  * 1. sort the blocks on length and then shrink the list
  *    so that the sum of lengths is length of sequences
  * 2. sort the blocks on distance from starts of sequences
  * 3. set each blocks score to its distance from the nearest edge
  * 4. find the best score as this start score plus match length
  * 5. and note the block number
  * 6. for each block look at all the previous ones to find its best
  *    predecessor. This is given by: best_score + length - diag_shift
  *    note the best score and block number as we proceed
  *    This leaves us with a linked list of the best blocks to use
  *    to produce a complete alignment.
  * 7. if the sum of the lengths of these blocks is > 20% of the diagonal
  *    length use align_wrap to produce a complete alignment
  */

int align_blocks ( Hash *h, ALIGN_PARAMS *params, OVERLAP *overlap ) {
    int i,j,l,gap_pen,diag_shift,best_score,best_prev,t,tt;
    int good_blocks, first_block;
    int *index_ptr = NULL;
    double best_percent;

    gap_pen = -1;
    best_score = -1000000;
    best_prev = -1;

    if ( h->matches < 1 ) return 0;

    /* sort the blocks on length and then shrink the list
     * so that the sum of lengths is length of sequences
     */

    i = sort_len_blocks(h->block_match, h->matches);
    l = MIN(h->seq1_len,h->seq2_len);
    for(i=0,j=0;i<h->matches;i++) {
        j+=h->block_match[i].length;
        if(j>l) {
            h->matches = i+1;
            break;
        }
    }

    /* sort the blocks on distance from starts of sequences */

    i = sort_blocks(h->block_match, h->matches);

    /* set each blocks score to its distance from the nearest edge
     * find the best score as this start score plus match length
     * and note the block number
     */

    for (i=0;i<h->matches;i++) {
        gap_pen = 
        -MIN(h->block_match[i].pos_seq1,h->block_match[i].pos_seq2);
        if((t=h->block_match[i].length + gap_pen) > best_score) {
            best_score = t;
            best_prev = i;
        }
        h->block_match[i].best_score = gap_pen;
        h->block_match[i].prev_block = -1;
    }

    if (best_prev == -1) return 0; /* bail out if there is no good score */
    /*
    printf("\n before %d\n",h->matches);
    for (i=0;i<h->matches;i++) {
        printf("i %d %d %d %d %d %d %d\n",i,
               h->block_match[i].pos_seq1,
               h->block_match[i].pos_seq2,
               h->block_match[i].length,
               h->block_match[i].diag,
               h->block_match[i].best_score,
               h->block_match[i].prev_block);

    }
    */

    /* for each block look at all the previous ones to find its best
     * predecessor. This is given by: best_score + length - diag_shift
     * note the best score and block number as we proceed
     */

    for (i=1;i<h->matches;i++) {
        for(j=i-1;j>-1;j--) {

            if( ((h->block_match[j].pos_seq1 + h->block_match[j].length)
                 <= h->block_match[i].pos_seq1)
                && ((h->block_match[j].pos_seq2 + h->block_match[j].length)
                 <= h->block_match[i].pos_seq2)) {

                diag_shift = abs(h->block_match[i].diag - 
                                 h->block_match[j].diag);

                /* best score doe snot include current match */
                if ( (t = h->block_match[j].best_score +
                      h->block_match[j].length - diag_shift) >
                    h->block_match[i].best_score ) {
                    h->block_match[i].best_score = t;
                    h->block_match[i].prev_block = j;
                    tt = t+h->block_match[i].length;
                    /*printf("i %d j %d t %d tt %d\n",i,j,t,tt);*/
                    if (tt>best_score) {
                        best_score = tt;
                        best_prev = i;
                    }
                }
            }
        }
    }

    /* best_prev is now last block */

    /* shuffle the ordered blocks to the start of the array */

    tt = h->block_match[best_prev].best_score;
    h->block_match[best_prev].best_score = -1;

    good_blocks = 1;
    first_block = 0;
    for (i=best_prev,j=0;i>-1;) {
        j = h->block_match[i].prev_block;
        if (j>-1) {
            good_blocks++;
            first_block = j;
        }
        i=h->block_match[i].prev_block;
    }

    if ( ! ( index_ptr = (int *) xmalloc ( sizeof(int)*(good_blocks) ))) {
        return -1;
    }

    for (i=best_prev,j=good_blocks-1;i>-1;j--) {
        index_ptr[j] = i;
        i=h->block_match[i].prev_block;
    }

    h->block_match[best_prev].best_score = tt;

    for (j=0;j<good_blocks;j++) {
        i = index_ptr[j];
        if (i != j) {
            h->block_match[j].pos_seq1   = h->block_match[i].pos_seq1;
            h->block_match[j].pos_seq2   = h->block_match[i].pos_seq2;
            h->block_match[j].length     = h->block_match[i].length;
            h->block_match[j].diag       = h->block_match[i].diag;
            h->block_match[j].best_score = h->block_match[i].best_score;
            h->block_match[j].prev_block = h->block_match[i].prev_block;
        }
    }

    if ( index_ptr ) xfree (index_ptr);
    h->matches = good_blocks;
    /*printf("returned %d matches with best score %d\n",h->matches,best_score);*/

    /*
    tt = 0;
    for (i=0;i<h->matches;i++) {
        tt += h->block_match[i].length;
        printf("i %d %d %d %d %d %d %d\n",i,
               h->block_match[i].pos_seq1,
               h->block_match[i].pos_seq2,
               h->block_match[i].length,
               h->block_match[i].diag,
               h->block_match[i].best_score,
               h->block_match[i].prev_block);
    }
    */
    
           /*
              could set 2 scores:
              the % of the diagonal which is covered by blocks
              the % of the region between and including the two
              outermost blocks ( ie a local or repeat score)
              local is best score + distance to first block
              global is best score - distance to right edge
              */

    i = h->matches/2;
    i = h->block_match[i].diag;
    j = diagonal_length(h->seq1_len, h->seq2_len, i);
    best_percent = 
           100.0*(double)(best_score - h->block_match[0].best_score) /
           (double)j;
    /*
    printf("local %d global %d %f\n",
           best_score - h->block_match[0].best_score,
           best_score - 
           MIN( (h->seq1_len - (h->block_match[h->matches-1].pos_seq1 +
                 h->block_match[h->matches-1].length)),
                (h->seq2_len - (h->block_match[h->matches-1].pos_seq2 +
                 h->block_match[h->matches-1].length))),best_percent);
                 */
/*
    overlap->seq1 = h->seq1;
    overlap->seq2 = h->seq2;
    overlap->seq1_len = h->seq1_len;
    overlap->seq2_len = h->seq2_len;
*/
    if (best_percent > 20.0) {
        i = align_wrap (h,params,overlap);
        if (i) return i;
        /*if (!i) print_overlap(overlap,stdout);*/
    }
    else {
        return 0;
    }
    return 1;
}

/**
  * compare a pair of sequences using hashing
  */

int compare_seqs(Hash *h, int *seq1_match_pos, int *seq2_match_pos,
                int *match_length) {
    
    int         ncw, nrw, word, pw1, pw2, i, j, match_size;
    int         diag_pos, size_hist;
    
    if(h->seq1_len < h->min_match) return -4; 
    if(h->seq2_len < h->min_match) return -4; 

    size_hist = h->seq1_len + h->seq2_len - 1;
    for(i = 0; i < size_hist; i++) h->diag[i] = -(h->word_length);
    
    nrw = h->seq2_len - h->word_length + 1;
    
    /*  loop for all (nrw) complete words in values2 */

    h->matches = -1;
    for (pw2=0;pw2<nrw;pw2++) {
        word = h->values2[pw2];
        if ( -1 != word ) {
            if ( 0 != (ncw = h->counts[word]) ) {
                for (j=0,pw1=h->last_word[word];j<ncw;j++) {
                    diag_pos = h->seq1_len - pw1 + pw2 - 1;
                    if ( h->diag[diag_pos] < pw2 ) {
                        if ((match_size = match_len ( 
                                                       h->seq1, pw1, h->seq1_len,
                                                       h->seq2, pw2, h->seq2_len))
                            >= h->min_match ) {
                            h->matches++;
                            if(h->max_matches == h->matches) {
                                return -5;
                            }
                            seq1_match_pos[h->matches] = pw1+1;
                            seq2_match_pos[h->matches] = pw2+1;
                            match_length[h->matches] = match_size;
                        }
                        h->diag[diag_pos] = pw2 + match_size;
                    }
                    pw1 = h->values1[pw1];
                }
            }
        }
    }
    h->matches += 1;
    return h->matches;
}



    /**Accumulates all matches between two sequences and then assess their
      * diagonal scores to find the most significant. Finds the centre
      * of gravity of these and uses that as the focus of a dynamic 
      * programming alignment using affine_align
      * used for finding poor matches in fij
     */
    
int compare_c(Hash *h,
            ALIGN_PARAMS *params, OVERLAP *overlap) {
    
    int         ncw, nrw, word, pw1, pw2, i, j, match_length;
    int         diag_pos, size_hist;
    int         hist_left, hist_right;
    
    int band, band_in;
    int         match_found;
    
    if(h->seq1_len < h->word_length) return -4; 
    if(h->seq2_len < h->word_length) return -4; 

    band_in = params->band;
    size_hist = h->seq1_len + h->seq2_len - 1;
    for(i = 0; i < size_hist; i++) h->diag[i] = -(h->word_length);
    for(i = 0; i < size_hist; i++) h->hist[i] = 0;
    
    nrw = h->seq2_len - h->word_length + 1;

    /*  loop for all (nrw) complete words in values2 */
    match_found = 0;
    for (pw2=0;pw2<nrw;pw2++) {
        word = h->values2[pw2];
        if ( -1 != word ) {
            if ( 0 != (ncw = h->counts[word]) ) {
                for (j=0,pw1=h->last_word[word];j<ncw;j++) {
                    diag_pos = h->seq1_len - pw1 + pw2 - 1;
                    if ( h->diag[diag_pos] < pw2 ) {
                        match_length = match_len (h->seq1, pw1, h->seq1_len,
                                                  h->seq2, pw2, h->seq2_len);
                        /* the old way:
                        if ((match_length = match_len ( 
                                                       h->seq1, pw1, h->seq1_len,
                                                       h->seq2, pw2, h->seq2_len))
                            >= h->min_match ) {
                            match_found = 1;
                        }
                        */
                        h->hist[diag_pos] += 1 + match_length - h->word_length;
                        h->diag[diag_pos] = pw2 + match_length;
                    }
                    pw1 = h->values1[pw1];
                }
            }
        }
    }
    match_found = 1;
    if(match_found) {
        /* now find the most significant diagonals */

        hist_left = MINMAT - 1;
        hist_right = size_hist - MINMAT;
        h->matches = -1;
        for(i = hist_left; i<hist_right; i++) {
            j = diagonal_length(h->seq1_len, h->seq2_len, i); 
            if(h->hist[i] > h->expected_scores[j]) { 
                h->matches++;
                /*printf("match %d %d %e %e\n",i,j,h->hist[i],h->expected_scores[j]);*/
                if(h->max_matches == h->matches) {
                    printf("too many matches %d\n",h->max_matches);
                    return -5;
                }
                h->diag_match[h->matches].pos = i;
                h->diag_match[h->matches].prob = h->hist[i] / (double)j;
                /*h->diag_match[h->matches].prob = h->hist[i] / h->expected_scores[j];*/
            }
        }
        h->matches += 1;
        if (h->matches < 1) return 0;

        /* get the best intercept defined by these diagonals */

        if (!best_intercept ( h, &pw1, &pw2))
          return 0;

        band = 0;
        if (band_in) {
            double perc;
            perc = (double)band_in / 100.0;
            j = MIN(h->seq1_len+1-pw1,h->seq2_len+1-pw2);
            band = (int) MAX(20.0,(perc*j));
            /*printf("band_in %d j %d perc %f band %d\n",band_in,j,perc,band);*/
        }
        set_align_params_banding (params, band, pw1, pw2);
        /*print_align_params(params);*/
        i = affine_align(overlap,params);
        params->band = band_in;
        if ( i ) return -1;
        /*print_overlap_struct(overlap);*/
        return 1;
    }
    
    return 0;
}

/**
  * sequence alignment using block matching with dynamic programming
  * for the bits between the blocks
  * good for similar sequences and very fast
  */

int compare_b(Hash *h,
            ALIGN_PARAMS *params, OVERLAP *overlap) {
    
    /* block alignment algorithm */

    int         ncw, nrw, word, pw1, pw2, i, j, match_size;
    int         diag_pos, size_hist;
    int job_in;

    if(h->seq1_len < h->min_match) return -4; 
    if(h->seq2_len < h->min_match) return -4; 

    size_hist = h->seq1_len + h->seq2_len - 1;
    for(i = 0; i < size_hist; i++) h->diag[i] = -(h->word_length);
    
    nrw = h->seq2_len - h->word_length + 1;
    
    /*  loop for all (nrw) complete words in values2 */

    h->matches = -1;
    for (pw2=0;pw2<nrw;pw2++) {
        word = h->values2[pw2];
        if ( -1 != word ) {
            if ( 0 != (ncw = h->counts[word]) ) {
                for (j=0,pw1=h->last_word[word];j<ncw;j++) {
                    diag_pos = h->seq1_len - pw1 + pw2 - 1;
                    if ( h->diag[diag_pos] < pw2 ) {
                        if ((match_size = match_len ( 
                                                       h->seq1, pw1, h->seq1_len,
                                                       h->seq2, pw2, h->seq2_len))
                            >= h->min_match ) {
                            h->matches++;
                            if(h->max_matches == h->matches) {
                                return -5;
                            }

                            h->block_match[h->matches].pos_seq1 = pw1;
                            h->block_match[h->matches].pos_seq2 = pw2;
                            h->block_match[h->matches].length = match_size;
                            h->block_match[h->matches].diag = diag_pos;
                        }
                        h->diag[diag_pos] = pw2 + match_size;
                    }
                    pw1 = h->values1[pw1];
                }
            }
        }
    }
    h->matches += 1;
    if ( h->matches < 1 ) return 0;
    job_in = params->return_job;
    params->return_job = 3; /* force return of edit buffers */
    /*printf(" blocks %d\n",h->matches);*/
    pw2 = align_blocks ( h, params, overlap );
    params->return_job = job_in;
    return pw2;
}


    /**Takes two sequences and attempts to find word matches between
     * them by hashing. The hashing is done for a certain word length
     * and then the program attempts to find all cases where the
     * word can be extended to >= min_match. The position of the centre
     * of gravity of all the matching blocks above min_match is found
     * and extrapolated back to the x or y axis of the comparison matrix. 
     * This is used to calculate the limits of the band used by the alignment
     * routine which are returned in params
     */
    
int compare_d(Hash *h,
            ALIGN_PARAMS *params, OVERLAP *overlap) {
    
    int         ncw, nrw, word, pw1, pw2, i, j, match_size;
    int         diag_pos, size_hist;
    int band_in;


    if(h->seq1_len < h->min_match) return -4; 
    if(h->seq2_len < h->min_match) return -4; 
    band_in = params->band;
    size_hist = h->seq1_len + h->seq2_len - 1;
    for(i = 0; i < size_hist; i++) h->diag[i] = -(h->word_length);
    
    nrw = h->seq2_len - h->word_length + 1;
    
    /*  loop for all (nrw) complete words in values2 */

    h->matches = -1;
    for (pw2=0;pw2<nrw;pw2++) {
        word = h->values2[pw2];
        if ( -1 != word ) {
            if ( 0 != (ncw = h->counts[word]) ) {
                for (j=0,pw1=h->last_word[word];j<ncw;j++) {
                    diag_pos = h->seq1_len - pw1 + pw2 - 1;
                    if ( h->diag[diag_pos] < pw2 ) {
                        if ((match_size = match_len ( 
                                                       h->seq1, pw1, h->seq1_len,
                                                       h->seq2, pw2, h->seq2_len))
                            >= h->min_match ) {
                            h->matches++;
                            if(h->max_matches == h->matches) {
                                return -5;
                            }
                            h->diag_match[h->matches].pos = diag_pos;
                            h->diag_match[h->matches].prob = 
                                (double)match_size/(double)diagonal_length(h->seq1_len, h->seq2_len, diag_pos); 
                        }
                        h->diag[diag_pos] = pw2 + match_size;
                    }
                    pw1 = h->values1[pw1];
                }
            }
        }
    }
    h->matches += 1;
    if ( h->matches < 1 ) return 0;
    if (!best_intercept ( h, &pw1, &pw2)) return 0;

    set_align_params_banding (params, band_in, pw1, pw2);

    return 1;
}

/**
  * finds all matching words between two sequences and returns them in an array
  * for foward matches removes the main diagonal
  */

int reps(Hash *h, int *seq1_match_pos, int *seq2_match_pos,
                int *match_length, char sense) {
    
    int         ncw, nrw, word, pw1, pw2, i, j, match_size;
    int         diag_pos, size_hist;
    

    if(h->seq1_len < h->min_match) return -4; 
    if(h->seq2_len < h->min_match) return -4; 

    size_hist = h->seq1_len + h->seq2_len - 1;
    for(i = 0; i < size_hist; i++) h->diag[i] = -(h->word_length);

    /* if forward repeats make sure we do not bother with the main diagonal */

    if ( 'f' == sense ) h->diag[h->seq1_len-1] = h->seq1_len;
    
    nrw = h->seq2_len - h->word_length + 1;
    
    /*  loop for all (nrw) complete words in values2 */

    h->matches = -1;
    for (pw2=0;pw2<nrw;pw2++) {
        word = h->values2[pw2];
        if ( -1 != word ) {
            if ( 0 != (ncw = h->counts[word]) ) {
                for (j=0,pw1=h->last_word[word];j<ncw;j++) {
                    diag_pos = h->seq1_len - pw1 + pw2 - 1;
                    if ( h->diag[diag_pos] < pw2 ) {
                        if ((match_size = match_len ( 
                                                       h->seq1, pw1, h->seq1_len,
                                                       h->seq2, pw2, h->seq2_len))
                            >= h->min_match ) {
                            h->matches++;
                            if(h->max_matches == h->matches) {
                                return -5;
                            }
                            seq1_match_pos[h->matches] = pw1+1;
                            seq2_match_pos[h->matches] = pw2+1;
                            match_length[h->matches] = match_size;
                        }
                        h->diag[diag_pos] = pw2 + match_size;
                    }
                    pw1 = h->values1[pw1];
                }
            }
        }
    }
    h->matches += 1;
    if ( h->matches ) {
        if ( sense == 'r' ) {
            (void) make_reverse ( seq2_match_pos, match_length,
                                 h->matches, h->seq2_len );     
        }
        (void) remdup ( seq1_match_pos, seq2_match_pos, match_length, &h->matches );
    }
    return h->matches;
}


}
