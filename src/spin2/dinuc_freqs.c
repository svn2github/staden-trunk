#include <ctype.h>
#include <stdio.h>

#include "dna_utils.h"

void calc_dinuc_freqs ( char seq[],
		       int user_start, int user_end,
		       double freqs[5][5]) {
    /* calculate dinucleotide frequencies from user_start to user_end */

    int i,j;

    for ( i=0;i<5;i++ ) {
	for ( j=0;j<5;j++ ) {
	    freqs[i][j] = 0.0;
	}
    }

    /* quick sanity check */

    if ( user_end - user_start < 1 ) return;

    for ( i=user_start-1; i< user_end-1; i++ ) {
	freqs[dna_lookup[seq[i]]][dna_lookup[seq[i+1]]]+=1.0;
    }

    for ( i=0;i<5;i++ ) {
	for ( j=0;j<5;j++ ) {
	    freqs[i][j] /= (user_end - user_start) / 100.0;
	}
    }

}

void calc_expected_dinuc_freqs ( char seq[],
		       int user_start, int user_end,
		       double freqs[5][5]) {

    /* calculate expected dinucleotide frequencies 
     * from user_start to user_end 
     * assuming random base order
     */

    int i,j;
    double comp[5];

    for ( i=0;i<5;i++ ) {
	comp[i] = 0.0;
	for ( j=0;j<5;j++ ) {
	    freqs[i][j] = 0.0;
	}
    }

    /* quick sanity check */

    if ( user_end - user_start < 1 ) return;

    for ( i=user_start-1; i< user_end-1; i++ ) {
	comp[dna_lookup[seq[i]]]+=1.0;
    }

    for ( i=0;i<5;i++ ) {
	comp[i] /= (double) (user_end - user_start);
    }

    for ( i=0;i<5;i++ ) {
	for ( j=0;j<5;j++ ) {
	    freqs[i][j] = comp[i] * comp[j] * 100.0;
	}
    }

}

int print_dinuc_table ( FILE *fp, double freqs[5][5] ) {
    char base[]="acgt";
    int i,j;
    if ( fprintf(fp, "       a       c       g       t\n") < 0 ) return 1;
    for ( i=0;i<4;i++ ) {
	if ( fprintf(fp, " %c",base[i]) < 0 ) return 1;
	for ( j=0;j<4;j++ ) {
	    if ( fprintf(fp, " %7.2f",freqs[i][j]) < 0 ) return 1;
	}
	if ( fprintf(fp, "\n") < 0 ) return 1;
    }
    /* caller does the close */
    return 0;
}

#ifdef REMOVE
main(){

    double freqs[5][5];
    char seq[]="-aaaaattttcccgg-";
    int i;
    set_dna_lookup();
    calc_dinuc_freqs ( seq, 2,15,freqs);
    i = print_dinuc_table ( stdout, freqs );
    calc_expected_dinuc_freqs ( seq, 2,15,freqs);
    i = print_dinuc_table ( stdout, freqs );
}
#endif
