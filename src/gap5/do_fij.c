#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "fij.h"
#include "misc.h" 
#include "dna_utils.h"
#include "align_lib.h"
/*#include "hash_lib.h"*/

/* 1/6/98 johnt - need to explicitly import globals from DLLs with Visual C++ */
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif

extern DLL_IMPORT int unknown_char;

void buffij(int left_gel1, int seq2_start, int left_gel2, int seq1_start,
	    int len_align, int score, double percent_mismatch);

int do_it_fij ( char seq[], int seq_len,
		int word_len, int min_overlap,
		double max_percent_mismatch, int compare_mode,
		int band, int gap_open, int gap_extend, double max_prob,
		int min_match, int max_alignment,
		Contig_parms *contig_list, int number_of_contigs) {

    int ret, i, j, longest_diagonal;
    int max_contig;
    int seq1_len, seq2_len, contig1_num, contig2_num;
    double percent_mismatch;
    static char buf[1024],name1[10],name2[10];
    int max_seq;
    int seq1_start_f, seq2_start_f, seq1_start_r, seq1_end_r, seq2_start_r;
    double comp[5];
    int max_matches;
    char *depad_seq1    = NULL, *depad_seq2    = NULL;
    int  *depad_to_pad1 = NULL, *depad_to_pad2 = NULL;
    int edge_mode, job, seq1_start, seq2_start;
    int compare_method;
    Hash *h;
    OVERLAP	*overlap;
    ALIGN_PARAMS *params;

    if (NULL == (params = create_align_params())) return -1;

    edge_mode = 10;
    seq1_start = 0;
    seq2_start = 0;
    job = 11;

    if (set_align_params (params, band, gap_open, gap_extend, edge_mode, job,
			  seq1_start, seq2_start,0,0,0)) {
		      destroy_alignment_params (params);
		      return -1;
		      };


    if (NULL == (overlap = create_overlap())) {
	destroy_alignment_params (params);
	return -1;
    }
    init_overlap (overlap, seq, seq, seq_len, seq_len);


    /* 
     * find longest contig 
     * initialise hashing
     * compare each contig in each orientation
     */

    if ( 2 > number_of_contigs ) {
	destroy_alignment_params (params);
	destroy_overlap (overlap);
	return -2;
    }
    for ( i=0,j=0,max_contig=0;i<number_of_contigs;i++) {
	if ( contig_list[i].contig_end_offset
	    - contig_list[i].contig_start_offset > max_contig ) {
	    j = i;
	    max_contig = contig_list[i].contig_end_offset
		       - contig_list[i].contig_start_offset;
	}
    }
    max_contig += 1;

    /*
    for ( i=0,longest_diagonal=0;i<number_of_contigs;i++) {
	if ( i != j ) {
	    if ( contig_list[i].contig_end_offset
		- contig_list[i].contig_start_offset > longest_diagonal ) {
		longest_diagonal = contig_list[i].contig_end_offset
		    - contig_list[i].contig_start_offset;
	    }
	}
    }
    */
    longest_diagonal = max_contig;
    max_seq = 2 * max_contig + 1;
    if (NULL == (depad_seq1 = (char *) xmalloc(sizeof(char) * max_contig))
	|| NULL == (depad_seq2 = (char *) xmalloc(sizeof(char) * max_contig))
	|| NULL == (depad_to_pad1 = (int *) xmalloc(sizeof(int) * max_contig))
	|| NULL == (depad_to_pad2 = (int *) xmalloc(sizeof(int) * max_contig))) {

	if (depad_seq1) xfree(depad_seq1);
	if (depad_seq2) xfree(depad_seq2);
	if (depad_to_pad1) xfree(depad_to_pad1);
	if (depad_to_pad2) xfree(depad_to_pad2);
	destroy_alignment_params (params);
	destroy_overlap (overlap);

	return -2;
    }

    max_matches = longest_diagonal;

    if ( min_match ) {
	compare_method = 17;
    }
    else {
	compare_method = 31;
    }
    if ( init_hash8n ( max_contig, longest_diagonal,
		     word_len, max_matches, min_match, compare_method, &h )) {
	free_hash8n(h);
	destroy_alignment_params (params);
	destroy_overlap (overlap);
	if (depad_seq1) xfree(depad_seq1);
	if (depad_seq2) xfree(depad_seq2);
	if (depad_to_pad1) xfree(depad_to_pad1);
	if (depad_to_pad2) xfree(depad_to_pad2);
	return -1;
    }
    
    if ( HASH_JOB_EXPD & compare_method ) {

        p_comp(comp,seq,seq_len);

	if(poisson_diagonals(MINMAT, longest_diagonal,
			     h->word_length, max_prob, h->expected_scores, comp)) {
	    free_hash8n(h);
	    destroy_alignment_params (params);
	    destroy_overlap (overlap);
	    if (depad_seq1) xfree(depad_seq1);
	    if (depad_seq2) xfree(depad_seq2);
	    if (depad_to_pad1) xfree(depad_to_pad1);
	    if (depad_to_pad2) xfree(depad_to_pad2);
	    return -1;
	}
    }

    for ( contig1_num = 0; contig1_num < number_of_contigs; contig1_num++ ) {
	seq1_len = contig_list[contig1_num].contig_end_offset 
	    - contig_list[contig1_num].contig_start_offset + 1;

	if ( seq1_len >= min_overlap ) {

	    /* Depad seq1 */
	    
	    overlap->seq1 = &seq[(contig_list[contig1_num].contig_start_offset)];
	    copy_seq(depad_seq1, overlap->seq1, seq1_len);
	    depad_seq(depad_seq1, &seq1_len, depad_to_pad1);
	    overlap->seq1 = h->seq1 = depad_seq1;
	    h->seq1_len = overlap->seq1_len = seq1_len;

	    if (hash_seqn(h, 1)) {
		verror(ERR_WARN, "find internal joins",
		       "hashing 1st sequence (#%d)",
		       contig_list[contig1_num].contig_left_gel); 
		continue;
	    }
    
	    (void) store_hashn ( h );

	    for ( contig2_num = contig1_num+1; contig2_num < number_of_contigs; 
		 contig2_num++ ) {

		seq2_len = contig_list[contig2_num].contig_end_offset 
		    - contig_list[contig2_num].contig_start_offset + 1;

		if ( seq2_len >= min_overlap ) {
		    overlap->seq2 = &seq[(contig_list[contig2_num].contig_start_offset)];
		    /* depad seq2.  This isn't really very efficient as we will
		       depad the sequences lots of times, but it uses less memory */

		    copy_seq(depad_seq2, overlap->seq2, seq2_len);
		    depad_seq(depad_seq2, &seq2_len, depad_to_pad2);
		    overlap->seq2 = h->seq2 = depad_seq2;
		    h->seq2_len = overlap->seq2_len = seq2_len;

		    if (hash_seqn(h, 2)) {
			verror(ERR_WARN, "find internal joins",
			       "hashing 2nd sequence (#%d)",
			       contig_list[contig2_num].contig_left_gel); 
			continue;
		    }

		    /*printf("%d %d\n",contig_list[contig1_num].contig_left_gel,contig_list[contig2_num].contig_left_gel);*/
				    
		    ret = 0;
		    if ( compare_method != 17 ) {
			/*
			 * Accept that some matches may just be too slow
			 * or memory hungry for sensitive mode, so we fall
			 * back to fast methods when this fails.
			 */
			ret = compare_a ( h, params, overlap );
			if (ret < 0) {
			    verror(ERR_WARN, "find internal joins",
				   "alignment too large for sensitive mode;"
				   " falling back to quick alignment"); 
			}
		    }
		    if ( ret < 0 || compare_method == 17 ) {
			ret = compare_b ( h, params, overlap );
		    }

		    if ( ret < 0 ) {
			verror(ERR_WARN, "find internal joins", "hashing" ); 
			continue;
		    }
		    
		    if ( ret ) {

			percent_mismatch = 100.0 - overlap->percent;

			if ( (overlap->length >= min_overlap) && 
			    (percent_mismatch <= max_percent_mismatch )) {

			    /* note conversion depadded to padded coordinates */

			    seq1_start_f = depad_to_pad1[overlap->left2]
				- contig_list[contig1_num].contig_left_extension + 1 ; 
			    /* add 1 to get base number */
			    seq2_start_f = depad_to_pad2[overlap->left1]
				- contig_list[contig2_num].contig_left_extension + 1 ; 
			
			    sprintf(name1,"%d",
				    contig_list[contig1_num].contig_left_gel);
			    sprintf(name2,"%d",
				    contig_list[contig2_num].contig_left_gel);
			    sprintf(buf,
				    " Possible join between contig %d "
				    "in the + sense and contig %d\n"
				    " Length %d",
				    contig_list[contig1_num].contig_left_gel,
				    contig_list[contig2_num].contig_left_gel,
				    overlap->length);
			    
			    /* Oops.  The initial coordinates in list_alignment are
			       padded, but then the alignment is depadded.  Hopefully
			       no-one will notice ! */

			    overlap->seq1_out[overlap->right+1] = '\0';
			    overlap->seq2_out[overlap->right+1] = '\0';

			    if (overlap->length <= max_alignment) {
				ret = list_alignment(&overlap->seq1_out[overlap->left],
						     &overlap->seq2_out[overlap->left],
						     name1,name2,seq1_start_f,seq2_start_f,buf);
			    } else {
				vmessage("%s\n", buf);
				vmessage(" Percentage mismatch %5.1f\n\n",
					 percent_mismatch);
			    }

			    buffij(
				   contig_list[contig1_num].contig_left_gel,seq2_start_f,
				   contig_list[contig2_num].contig_left_gel,seq1_start_f,
				   overlap->length, (int)overlap->score,
				   percent_mismatch);
			}
		    }
		}
		free_overlap(overlap);
	    }
	}
	if ( compare_mode == COMPARE_SINGLE ) break;
    }

    /* now do complementary strand */

    for ( contig1_num = 0; contig1_num < number_of_contigs; contig1_num++ ) {
	seq1_len = contig_list[contig1_num].contig_end_offset 
	    - contig_list[contig1_num].contig_start_offset + 1;

	if ( seq1_len >= min_overlap ) {

	    /* Depad seq1 */
	    
	    overlap->seq1 = &seq[(contig_list[contig1_num].contig_start_offset)];
	    copy_seq(depad_seq1, overlap->seq1, seq1_len);
	    complement_seq(depad_seq1, seq1_len);
	    depad_seq(depad_seq1, &seq1_len, depad_to_pad1);
	    overlap->seq1 = h->seq1 = depad_seq1;
	    h->seq1_len = overlap->seq1_len = seq1_len;

	    /*printf("seq1 %s\n",overlap->seq1);*/
	    if (hash_seqn(h, 1)) {
		verror(ERR_WARN, "find internal joins",
		       "hashing 1st sequence (#%d)",
		       contig_list[contig1_num].contig_left_gel); 
		continue;
	    }
    
	    (void) store_hashn ( h );
    
	    for ( contig2_num = contig1_num+1; contig2_num < number_of_contigs; 
		 contig2_num++ ) {

		seq2_len = contig_list[contig2_num].contig_end_offset 
		    - contig_list[contig2_num].contig_start_offset + 1;

		if ( seq2_len >= min_overlap ) {


		    overlap->seq2 = &seq[(contig_list[contig2_num].contig_start_offset)];
		    /* depad seq2.  This isn't really very efficient as we will
		       depad the sequences lots of times, but it uses less memory */

		    copy_seq(depad_seq2, overlap->seq2, seq2_len);
		    depad_seq(depad_seq2, &seq2_len, depad_to_pad2);
		    overlap->seq2 = h->seq2 = depad_seq2;
		    h->seq2_len = overlap->seq2_len = seq2_len;

		    if (hash_seqn(h, 2)) {
			verror(ERR_WARN, "find internal joins", 
			       "hashing 2nd sequence (#%d)",
			       contig_list[contig2_num].contig_left_gel); 
			continue;
		    }

		    
		    /*printf("%d %d\n",contig_list[contig1_num].contig_left_gel,contig_list[contig2_num].contig_left_gel);*/

		    if ( compare_method == 17 ) {
			ret = compare_b ( h, params, overlap );
		    }
		    else {
			ret = compare_a ( h, params, overlap );
		    }

		    if ( ret < 0 ) {
			verror(ERR_WARN, "find internal joins", "hashing" ); 
			continue;
		    }
		    
		    if ( ret ) {

			percent_mismatch = 100.0 - overlap->percent;

			if ( (overlap->length >= min_overlap) && 
			    (percent_mismatch <= max_percent_mismatch )) {

			    /* note conversion depadded to padded coordinates */

                            seq1_start_r = depad_to_pad1[overlap->left2]
			        - contig_list[contig1_num].contig_right_extension + 1 ; 
			    {
				int p = overlap->left2 + overlap->length - 1;
				int diff = 0;
				if (p > seq1_len-1) {
				    diff = p - (seq1_len-1);
				    p = seq1_len-1;
				}
				seq1_end_r = depad_to_pad1[p] + diff;
			    }
			    
		            seq2_start_r = depad_to_pad2[overlap->left1]
			        - contig_list[contig2_num].contig_left_extension + 1 ; 
			
			    sprintf(name1,"%d",
				    contig_list[contig1_num].contig_left_gel);
			    sprintf(name2,"%d",
				    contig_list[contig2_num].contig_left_gel);
			    sprintf(buf," Possible join between contig %d "
				    "in the - sense and contig %d\n"
				    " Length %d",
				    contig_list[contig1_num].contig_left_gel,
				    contig_list[contig2_num].contig_left_gel,
				    overlap->length);
			    overlap->seq1_out[overlap->right+1] = '\0';
			    overlap->seq2_out[overlap->right+1] = '\0';

			    if (overlap->length <= max_alignment) {
				ret = list_alignment(&overlap->seq1_out[overlap->left],
						     &overlap->seq2_out[overlap->left],
						     name1,name2,seq1_start_r,seq2_start_r,buf);
			    } else {
				vmessage("%s\n", buf);
				vmessage(" Percentage mismatch %5.1f\n\n",
					 percent_mismatch);
			    }
			
			    buffij(
				   -contig_list[contig1_num].contig_left_gel,
				   seq2_start_r,
				   contig_list[contig2_num].contig_left_gel,
				   seq1_end_r,
				   overlap->length, (int)overlap->score,
				   percent_mismatch);
			}
		    }
		}
		free_overlap(overlap);
	    }
	}
	if ( compare_mode == COMPARE_SINGLE ) break;
    }

    xfree(depad_seq1);
    xfree(depad_seq2);
    xfree(depad_to_pad1);
    xfree(depad_to_pad2);
    free_hash8n ( h );
    destroy_alignment_params (params);
    destroy_overlap (overlap);
    return 0;
}
