/*	routines for finding open reading frames

*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "dna_utils.h"
#include "array_arith.h"
#include "genetic_code.h"
#include "text_output.h"
#include "misc.h"


/* PUT ME IN array_arith */

int minimum_element ( int list[], int num_elements ) {

    /* return element number of first element in list that contains
     * the lists lowest value
     */

    int i, min;

    for ( i=1,min=list[0];i<num_elements;i++) {
	min = MIN((list[i]),min);
    }
    for ( i=0;i<num_elements;i++) {
	if ( min == list[i] ) return i;
    }

    return 0;
}

/* PUT ME in genetic_code */

int write_seq_lines ( FILE *fp, char *seq, int seq_len ) {
#define LINELENGTH 60

    /* writes seq to fp, 60 characters per line */

    int i;
    for(i=0;i<seq_len;i++) {
	if ( i && (i%LINELENGTH==0)) if ( fprintf(fp,"\n") < 0 ) return 1;
	if ( fprintf(fp,"%c",seq[i]) < 0 ) return 1;
    }
    if ( fprintf(fp,"\n") < 0 ) return 1;
    return 0;
}

int write_screen_seq_lines (char *seq, int seq_len) {
#define LINELENGTH 60

    /* writes seq to text output window, 60 characters per line */

    int i;
    for(i=0;i<seq_len;i++) {
	if ( i && (i%LINELENGTH==0)) vmessage("\n");
	vmessage("%c",seq[i]);
    }
    vmessage("\n");
    return 0;
}

char *orf_protein_seqf ( char *dna, int num_bases ) {

    /* return a string containing the translation up to the next
     * stop codon or end of data.
     */
    char *s;
    int i,j;

    if ( ( NULL == ( s = ( char* ) malloc ( sizeof (char) * (num_bases))))) return NULL;

    for ( i=0,j=0; i+2<num_bases; i+=3,j++ ) {
	if ( ( '*' == (s[j] = codon_to_acid1(&dna[i])))) {
	    j++;
	    break;
	}
    }

    if ( j > 0 && '*' != s[j-1] ) {       /* always return * as last symbol */
	s[j++] = '*';
    }

    s[j++] = '\0';

    /* Shrink memory usage (if realloc will allow it) */
    if ( ( NULL == ( s = ( char* ) realloc ( s, sizeof (char) * (j+1)))))
	return NULL;

    return s;
}

int write_open_frames_f ( FILE *fp, char *dna, int num_bases, 
			  int user_start,
			  int user_end, 
			  int min_open ) {

    /* finds open reading frames of length >= min_open and writes them
     * to fp in fasta format
     */

    int starts[3], frame, len_open;
    char *protein, line[80];

    /* min_open -= 1; */
    starts[0] = user_start-1;
    starts[1] = user_start;
    starts[2] = user_start+1;
    frame = 0;

    while ( starts[frame] < user_end - 3 * min_open ) {

	if ( min_open < (len_open = strlen ( protein = (orf_protein_seqf ( 
	           &dna[starts[frame]], user_end-starts[frame] ))))) {
	    memset(line,' ',80);
	    sprintf(line, ">%d", starts[frame]+1);
	    line[strlen(line)] = ' ';
	    sprintf(&line[21], "%d..%d", starts[frame]+1, starts[frame]+len_open*3-3);
	    if ( fprintf( fp,"%s\n",line) < 0 ) goto bail_out;
	    if (write_seq_lines ( fp, protein, len_open )) goto bail_out;
	}
	starts[frame] += 3 * (len_open) ;
	frame = minimum_element(starts,3);
	free ( protein );
    }
    return 0;
 bail_out:
    free ( protein );
    return 1;
}

int write_screen_open_frames_f ( char *dna, int num_bases, 
				 int user_start,
				 int user_end, 
				 int min_open ) {

    /* finds open reading frames of length >= min_open and writes them
     * to fp in fasta format
     */

    int starts[3], frame, len_open;
    char *protein, line[80];

    /* min_open -= 1; */
    starts[0] = user_start-1;
    starts[1] = user_start;
    starts[2] = user_start+1;
    frame = 0;

    while ( starts[frame] < user_end - 3 * min_open ) {

	if ( min_open < (len_open = strlen ( protein = (orf_protein_seqf ( 
	           &dna[starts[frame]], user_end-starts[frame] ))))) {
	    memset(line,' ',80);
	    sprintf(line, ">%d", starts[frame]+1);
	    line[strlen(line)] = ' ';
	    sprintf(&line[21], "%d..%d", starts[frame]+1, starts[frame]+len_open*3-3);
	    vmessage("%s\n", line);
	    if (write_screen_seq_lines ( protein, len_open )) goto bail_out;
	}
	starts[frame] += 3 * (len_open) ;
	frame = minimum_element(starts,3);
	free ( protein );
    }
    return 0;
 bail_out:
    free ( protein );
    return 1;
}

void write_open_frames_f_ft ( FILE *fp, char *dna, int num_bases, 
			      int user_start,
			      int user_end, 
			      int min_open ) {
    
    /* find open reading frames (plus strand) and write them to fp in 
     * embl ft format.
     */


    int starts[3], frame, len_open;
    char *protein, line[80];

    /* min_open -= 1; */
    starts[0] = user_start-1;
    starts[1] = user_start;
    starts[2] = user_start+1;
    frame = 0;

    while ( starts[frame] < user_end - 3 * min_open ) {

	if ( min_open < (len_open = strlen ( protein = (orf_protein_seqf ( 
	           &dna[starts[frame]], user_end-starts[frame] ))))) {
	    memset(line,' ',80);
	    strcpy( line, "FT   CDS");
	    line[8] = ' ';
	    sprintf(&line[21], "%d..%d", starts[frame]+1, starts[frame]+len_open*3-3);
	    if ( fprintf(fp, "%s\n",line) < 0 ) goto bail_out;
	}
	starts[frame] += 3 * (len_open) ;
	frame = minimum_element(starts,3);
	free ( protein );
    }
    return;
 bail_out:
    free ( protein );
}

void write_screen_open_frames_f_ft (char *dna, int num_bases, 
				    int user_start,
				    int user_end, 
				    int min_open ) {
    
    /* find open reading frames (plus strand) and write them to output window
     * in embl ft format.
     */


    int starts[3], frame, len_open;
    char *protein, line[80];

    /*  min_open -= 1; */
    starts[0] = user_start-1;
    starts[1] = user_start;
    starts[2] = user_start+1;
    frame = 0;

    while ( starts[frame] < user_end - 3 * min_open ) {

	if ( min_open < (len_open = strlen ( protein = (orf_protein_seqf ( 
	           &dna[starts[frame]], user_end-starts[frame] ))))) {
	    memset(line,' ',80);
	    strcpy( line, "FT   CDS");
	    line[8] = ' ';
	    sprintf(&line[21], "%d..%d", starts[frame]+1, starts[frame]+len_open*3-3);
	    vmessage("%s\n", line); 
	}
	starts[frame] += 3 * (len_open) ;
	frame = minimum_element(starts,3);
	free ( protein );
    }
    return;
}

char *orf_protein_seq_r ( char *dna, int num_bases ) {

    /* return a string containing the translation up to the next
     * stop codon or end of data.
     */

    char *s;
    int i,j;

    if ( ( NULL == ( s = ( char* ) malloc ( sizeof (char) * (num_bases))))) return NULL;

    for ( i=0,j=0; i+2<num_bases; i+=3,j++ ) {
	if ( ( '*' == (s[j] = codon_to_cacid1(&dna[i])))) {
	    j++;
	    break;
	}
    }

    if ( j > 0 && '*' != s[j-1] ) {       /* always return * as last symbol */
	s[j++] = '*';
    }
    reverse_dna ( s, j-1 );
    s[j++] = '\0';

    if ( ( NULL == ( s = ( char* ) realloc ( s, sizeof (char) * (j+1)))))
	return NULL;

    return s;
}

int write_open_frames_r ( FILE *fp, char *dna, int num_bases, 
			    int user_start,
			    int user_end, 
			 int min_open ) {

    /* finds open reading frames of length >= min_open and writes them
     * to fp in fasta format  (minus strand)
     */

    int starts[3], frame, len_open;
    char *protein, line[80];

    /* min_open -= 1; */
    starts[0] = user_start-1;
    starts[1] = user_start;
    starts[2] = user_start+1;
    frame = 0;

    while ( starts[frame] < user_end - 3 * min_open ) {

	if ( min_open < (len_open = strlen ( protein = (orf_protein_seq_r ( 
	           &dna[starts[frame]], user_end-starts[frame] ))))) {
	    memset(line,' ',80);
	    sprintf(line, ">%d", starts[frame]+1);
	    line[strlen(line)] = ' ';
	    sprintf(&line[21], "complement(%d..%d)", starts[frame]+1, starts[frame]+len_open*3-3);
	    if ( fprintf(fp, "%s\n",line) < 0 ) goto bail_out;
	    if (write_seq_lines ( fp, protein, len_open )) goto bail_out;
	}
	starts[frame] += 3 * (len_open) ;
	frame = minimum_element(starts,3);
	free ( protein );
    }
    return 0;
 bail_out:
    free ( protein );
    return 1;
}
int write_screen_open_frames_r ( char *dna, int num_bases, 
				 int user_start,
				 int user_end, 
				 int min_open ) {

    /* finds open reading frames of length >= min_open and writes them
     * to text output window in fasta format  (minus strand)
     */

    int starts[3], frame, len_open;
    char *protein, line[80];

    /* min_open -= 1; */
    starts[0] = user_start-1;
    starts[1] = user_start;
    starts[2] = user_start+1;
    frame = 0;

    while ( starts[frame] < user_end - 3 * min_open ) {

	if ( min_open < (len_open = strlen ( protein = (orf_protein_seq_r ( 
	           &dna[starts[frame]], user_end-starts[frame] ))))) {
	    memset(line,' ',80);
	    sprintf(line, ">%d", starts[frame]+1);
	    line[strlen(line)] = ' ';
	    sprintf(&line[21], "complement(%d..%d)", starts[frame]+1, starts[frame]+len_open*3-3);
	    vmessage("%s\n",line);
	    if (write_screen_seq_lines (protein, len_open )) goto bail_out;
	}
	starts[frame] += 3 * (len_open) ;
	frame = minimum_element(starts,3);
	free ( protein );
    }
    return 0;
 bail_out:
    free ( protein );
    return 1;
}

void write_open_frames_r_ft ( FILE *fp, char *dna, int num_bases, 
			      int user_start,
			      int user_end, 
			      int min_open ) {
    
    /* find open reading frames (minus strand) and write them to fp in 
     * embl ft format.
     */

    int starts[3], frame, len_open;
    char *protein, line[80];

    /* min_open -= 1; */
    starts[0] = user_start-1;
    starts[1] = user_start;
    starts[2] = user_start+1;

    frame = 0;

    while ( starts[frame] < user_end - 3 * min_open ) {

	if ( min_open < (len_open = strlen ( protein = (orf_protein_seq_r ( 
	           &dna[starts[frame]], user_end-starts[frame] ))))) {
	    memset(line,' ',80);
	    strcpy( line, "FT   CDS");
	    line[8] = ' ';
	    sprintf(&line[21], "complement(%d..%d)", starts[frame]+1, starts[frame]+len_open*3-3);
	    if ( fprintf(fp,"%s\n",line) < 0 ) goto bail_out;
	}
	starts[frame] += 3 * (len_open) ;
	frame = minimum_element(starts,3);
	free ( protein );
    }
    return;
 bail_out:
    free ( protein );
}

void write_screen_open_frames_r_ft (char *dna, int num_bases, 
				    int user_start,
				    int user_end, 
				    int min_open ) {
    
    /* find open reading frames (minus strand) and write them to text output
     * window in embl ft format.
     */

    int starts[3], frame, len_open;
    char *protein, line[80];

    /* min_open -= 1; */
    starts[0] = user_start-1;
    starts[1] = user_start;
    starts[2] = user_start+1;

    frame = 0;

    while ( starts[frame] < user_end - 3 * min_open ) {

	if ( min_open < (len_open = strlen ( protein = (orf_protein_seq_r ( 
	           &dna[starts[frame]], user_end-starts[frame] ))))) {
	    memset(line,' ',80);
	    strcpy( line, "FT   CDS");
	    line[8] = ' ';
	    sprintf(&line[21], "complement(%d..%d)", starts[frame]+1, starts[frame]+len_open*3-3);
	    vmessage("%s\n",line);
	}
	starts[frame] += 3 * (len_open) ;
	frame = minimum_element(starts,3);
	free ( protein );
    }
    return;
}

#ifdef REMOVE
main () {
/*                         1  2  3  4  5  6  7  8
               12345678901234567890123456789012345678901234567890123456*/
    char *dna="atcgcttaaacgatgttatagacagataacagatacagatagacagatacagtaaa";
    char *protein;
    int len_open, min_len_open;
    int i, num_bases;
    int lreg,rreg;

    set_dna_lookup();
    init_genetic_code ();
    printf("%s\n",dna);
    num_bases = strlen ( dna );
    printf("num_bases %d\n",num_bases);

    lreg = 1;
    rreg = num_bases;
    min_len_open = 6;
    write_open_frames_f_ft ( stdout, dna, num_bases, lreg, rreg, min_len_open );
    write_open_frames_r_ft ( stdout, dna, num_bases, lreg, rreg, min_len_open );
    write_open_frames_f ( stdout, dna, num_bases, lreg, rreg, min_len_open );
    write_open_frames_r ( stdout, dna, num_bases, lreg, rreg, min_len_open );
}
#endif
