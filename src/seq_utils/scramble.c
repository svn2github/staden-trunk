#include <stdlib.h>

#include "sequence_formats.h"
/* routines for scrambling a sequence */
/*
   send seq of length seq_len
   put in a scramble structure
   generate seq_len random numbers and 
   put them into the scramble structure
   sort the scramble structure
   copy the scrambled sequence to seq
*/
struct char_int_pair {
    char character;
    int  number;
};

int compare_pair(const void *vp1, const void *vp2) {
    struct char_int_pair *p1 = (struct char_int_pair *)vp1;
    struct char_int_pair *p2 = (struct char_int_pair *)vp2;

    if (p1->number < p2->number)
	return -1;
    else if (p1->number == p2->number)
	return 0;
    else
	return 1;
}
 
int scramble_seq ( char seq[], int seq_len, int seed ) {

    int i;
    struct char_int_pair *scrambler;

    if ( ! ( scrambler = (struct char_int_pair *) malloc ( sizeof(struct char_int_pair)*(seq_len) ))) {
        return -1;
    }

    srand(seed);
    for(i=0;i<seq_len;i++) {
	scrambler[i].character = seq[i];
	scrambler[i].number = rand();
    }

    qsort ((void *) scrambler, seq_len, sizeof(struct char_int_pair),
	   compare_pair);

    for(i=0;i<seq_len;i++) {
	seq[i] = scrambler[i].character;
    }

    free ( scrambler );
    return 0;
}
