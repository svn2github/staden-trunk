#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <os.h>
#include <misc.h>
#include "assert.h"
#include "consen.h"


/****************************************************************************/

/* Return 1 if base is not (a,c,g,t,A,C,G,T), ELSE 0 */

int unknown_base(char base) {

    static char known[] = {"acgtACGT"};
    int i,j;

    j = strlen ( known );
    for ( i=0; i<j; i++) {
	if ( base == known[i] )
	    return 0;
    }
    return 1;
}


/**********************************************************************/

int rotate_buf_old ( int buf_pos, int buf_size ) {

    int t;

    t = buf_pos + 1;
    if ( t == buf_size ) t = 0;
    return t;
}

#define rotate_buf(p,s) ((p+1)%s)

/**********************************************************************/

/* Return base number of last base of good data. ie index+1 of last
   base before the sequence contains too many unknown characters.

   Algorithm

   Use a rotating array to store the positions of the last several unknown 
   bases. Also note the element numbers in this array of the leftmost and
   rightmost unknown bases. If we have stored the positions of at least the 
   minimum number of unknown bases, and their distance apart is less than
   the window length, return the position of the current leftmost unknown
   base.

*/


int bad_data_start ( char *seq, int window_len, int max_unknown,
		    int seq_length, int dir) {

    int max_unknownp1, *unknown_ptr, leftu, rightu, num_bad, i;
    int istart, iend, cycle;

    cycle = max_unknownp1 = max_unknown + 1;
    leftu = 0; 			/* index of leftmost unknown in unknown_ptr */
    rightu = -1;		/* index of rightmost unknown in unknown_ptr */
    num_bad = 0;		/* count of unknowns in unknown_ptr */

    /* allocate space for rotating array */
    unknown_ptr = (int *) malloc ( cycle * sizeof ( int ) );
    if ( NULL == unknown_ptr ) return 0;

    if (dir == 1) {
	istart = 0;
	iend = seq_length;
    } else {
	istart = seq_length-1;
	iend = -1;
    }

    for (i = istart; i != iend; i += dir) {

	if ( unknown_base ( seq[i] ) ) {

	    if (dir == -1 && i <= window_len) {
		max_unknownp1 = max_unknown * ((float)i / window_len) + 1;
	    }

	    num_bad += 1;
	    rightu = ( rightu + 1 ) % cycle;
	    unknown_ptr[rightu] = i;

	    if ( num_bad >= max_unknownp1 ) {
		/*
		 * Got enough bad base positions stored in rotating buffer.
		 * Are the first and last within window_len of one another ? 
		 */
		/* assert((rightu+1)%cycle == leftu);*/
		if ( ABS(unknown_ptr[rightu] - unknown_ptr[leftu])
		    < window_len ) {
		    int ret = unknown_ptr[leftu];
		    free ( unknown_ptr );
		    return ret;
		}

		leftu = ( leftu + 1 ) % cycle;
	    }
	}
    }
    free ( unknown_ptr );
    return dir == 1 ? seq_length : -1;
}



/************************************************************/

int end_of_good ( char *seq, 
		 int start,
		 int window_len1, 
		 int max_unknown1, 
		 int window_len2,
		 int max_unknown2) {

    int window_len, max_unknown, jstart, bad_start;

    window_len = window_len1;
    max_unknown = max_unknown1;

    jstart = MIN ( start, strlen( seq ) );

    /* 	Do the first search from start on */
    bad_start = bad_data_start(&seq[jstart], window_len, max_unknown,
			       strlen(&seq[jstart]), 1);

    /* 	If required do another search from here on */
    if ( window_len2 ) {

	jstart = bad_start + jstart + window_len1/2;
	jstart = MIN ( jstart, strlen( seq ) );
	bad_start = bad_data_start(&seq[jstart], window_len2, max_unknown2,
				   strlen(&seq[jstart]), 1);

    }

    return bad_start + jstart;
}

int start_of_good ( char *seq, 
		   int start,
		   int window_len1, 
		   int max_unknown1, 
		   int window_len2,
		   int max_unknown2) {

    int window_len, max_unknown, jstart, bad_start;

    window_len = window_len1;
    max_unknown = max_unknown1;

    jstart = MAX ( start, 1 );

    /* 	Do the first search from start on */
    bad_start = bad_data_start(seq, window_len, max_unknown, jstart, -1)+2; 


    /* 	If required do another search from here on */
    if ( window_len2 ) {
	jstart = bad_start - window_len1/2;
	jstart = MAX ( jstart, 1 );
	bad_start = bad_data_start(seq, window_len2, max_unknown2, jstart, -1)+2;
    }

    return bad_start;
}
