#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>


#include "misc.h" 
#include "dna_utils.h"

char *seq_left_end ( char seq[], int seq_length, int start, int window_length,
		    int inc) {

    int num_elements, middle_start, half_window, left_start, right_end;
    char *edge_seq;
    int i, j;

    /* work in array element nos.
       middle_start first element that caller will remove from left edge
       left_start is first element we need to fill with dash or seq
       right_end is last element we need to fill with seq
    */


    middle_start = start;
    half_window = window_length/2;
    left_start = middle_start - inc * ( half_window / inc );
    right_end = middle_start + window_length - 1;
    num_elements = window_length + inc * ( half_window / inc );

    /* quick sanity checks */

    if ( right_end > seq_length ) return NULL;
    if ( start < 0 ) return NULL;

    if ( ( NULL == ( edge_seq = ( char* ) xmalloc ( sizeof (char) * (num_elements+1))))) return NULL;


    /* if required fill start with dashes */
    edge_seq[num_elements] = '\0';
    for ( i = 0, j = left_start; j < 0; i++, j++ ) {
	edge_seq[i] = '-';
    }

    /* fill rest with real sequence */

    for ( ; j <= right_end; i++, j++ ) {
	edge_seq[i] = seq[j];
    }
    return edge_seq;
}

char *seq_right_end ( char seq[], int seq_length, int end, int window_length,
		     int inc ) {

    int num_elements, left_start, right_end;
    char *edge_seq;
    int i, j, jend;

    /* work in array element nos.
       return the last base which the middle game used
       - so its score can be removed from the left edge
    */

    /* quick sanity checks */


    if ( end >= seq_length ) return NULL;
    if ( window_length > seq_length ) return NULL;

    left_start = end - window_length + 1;
    right_end  = end + window_length/2;
    if ( inc == 3 ) right_end += 1;
    num_elements = right_end - left_start + 1;

    if ( ( NULL == ( edge_seq = ( char* ) xmalloc ( sizeof (char) * (num_elements+1))))) return NULL;

    /* fill start with real sequence */

    edge_seq[num_elements] = '\0';
    jend = MAX(end+1,seq_length);
    for ( i = 0, j = left_start; (j < jend)&&( i<num_elements); i++, j++ ) {
	edge_seq[i] = seq[j];
    }

    /* if required fill rest with dashes */

    for ( ; j <= right_end; i++, j++ ) {
	edge_seq[i] = '-';
    }
    return edge_seq;
}


void get_base_comp ( char seq[], int seq_length, double base_comp[5] ) {
    int i;
    for ( i = 0; i < 5; i++ ) {
	base_comp[i] = 0.0;
    }
    for ( i = 0; i < seq_length; i++ ) {
	base_comp[ char_lookup [ (unsigned int)seq[i] ] ] += 1.0;
    }
}



