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

void buffij(tg_rec left_gel1, int seq2_start, tg_rec left_gel2, int seq1_start,
	    int len_align, int score, double percent_mismatch);

typedef struct {
    Contig_parms *contig_list;
    int number_of_contigs;
    int *depad_to_pad1;
    int *depad_to_pad2;
    int min_overlap;
    int max_percent_mismatch;
    int max_alignment;
    int seq2_len;
    int one_by_one;
} add_fij_t;

static void add_fij_overlap(OVERLAP *overlap, int contig1_num,
			    int contig2_num, void *clientdata) {
    double percent_mismatch;
    int seq1_start_f, seq2_start_f;
    static char buf[1024],name1[10],name2[10];
    add_fij_t *cd = (add_fij_t *)clientdata;
    int c1_s;

    c1_s = cd->one_by_one
	? 0 
	: cd->contig_list[contig1_num].contig_start_offset;

    percent_mismatch = 100.0 - overlap->percent;

    if ( (overlap->length >= cd->min_overlap) && 
	 (percent_mismatch <= cd->max_percent_mismatch )) {

	/* note conversion depadded to padded coordinates */
	seq1_start_f = cd->depad_to_pad1[overlap->left2 + c1_s]
	    - cd->depad_to_pad1[c1_s]
	    - cd->contig_list[contig1_num].contig_left_extension + 1 ; 

	/* add 1 to get base number */
	seq2_start_f = cd->depad_to_pad2[overlap->left1]
	    - cd->contig_list[contig2_num].contig_left_extension + 1 ; 
			
	sprintf(name1,"%"PRIrec,
		cd->contig_list[contig1_num].contig_left_gel);
	sprintf(name2,"%"PRIrec,
		cd->contig_list[contig2_num].contig_left_gel);
	sprintf(buf,
		" Possible join between contig %"PRIrec
		" in the + sense and contig %"PRIrec"\n"
		" Length %d",
		cd->contig_list[contig1_num].contig_left_gel,
		cd->contig_list[contig2_num].contig_left_gel,
		overlap->length);
			    
#if 0
	/* Oops.  The initial coordinates in list_alignment are
	   padded, but then the alignment is depadded.  Hopefully
	   no-one will notice ! */

	overlap->seq1_out[overlap->right+1] = '\0';
	overlap->seq2_out[overlap->right+1] = '\0';

	if (overlap->length <= cd->max_alignment) {
	    ret = list_alignment(&overlap->seq1_out[overlap->left],
				 &overlap->seq2_out[overlap->left],
				 name1,name2,seq1_start_f,seq2_start_f,buf);
	} else {
	    vmessage("%s\n", buf);
	    vmessage(" Percentage mismatch %5.1f\n\n",
		     percent_mismatch);
	}
#else
	vmessage("%s\n", buf);
	vmessage(" Percentage mismatch %5.1f\n\n",
		 percent_mismatch);
#endif
	    
	buffij(cd->contig_list[contig1_num].contig_left_gel,
	       seq2_start_f,
	       cd->contig_list[contig2_num].contig_left_gel,
	       seq1_start_f,
	       overlap->length, (int)overlap->score,
	       percent_mismatch);
    }

    //free_overlap(overlap);
}

static void add_fij_overlap_r(OVERLAP *overlap, int contig1_num,
			      int contig2_num, void *clientdata) {
    double percent_mismatch;
    int seq1_start_r, seq2_end_r, seq2_start_r;
    static char buf[1024],name1[10],name2[10];
    add_fij_t *cd = (add_fij_t *)clientdata;
    int c1_s;

    c1_s = cd->one_by_one
	? 0 
	: cd->contig_list[contig1_num].contig_start_offset;

    percent_mismatch = 100.0 - overlap->percent;

    if ( (overlap->length >= cd->min_overlap) && 
	 (percent_mismatch <= cd->max_percent_mismatch )) {

	/* note conversion depadded to padded coordinates */
	seq1_start_r = cd->depad_to_pad1[overlap->left2 + c1_s]
	    - cd->depad_to_pad1[c1_s]
	    - cd->contig_list[contig1_num].contig_left_extension + 1 ; 

	/* add 1 to get base number */
	seq2_start_r = cd->depad_to_pad2[overlap->left1]
	    - cd->contig_list[contig2_num].contig_right_extension + 1 ; 

	{
	    int p = overlap->left1 + overlap->length - 1;
	    int diff = 0;
	    if (p > cd->seq2_len-1) {
		diff = p - (cd->seq2_len-1);
		p = cd->seq2_len-1;
	    }
	    seq2_end_r = cd->depad_to_pad2[p] + diff;
	}
			
	sprintf(name1,"%"PRIrec,
		cd->contig_list[contig1_num].contig_left_gel);
	sprintf(name2,"%"PRIrec,
		cd->contig_list[contig2_num].contig_left_gel);
	sprintf(buf,
		" Possible join between contig %"PRIrec
		" in the - sense and contig %"PRIrec"\n"
		" Length %d",
		cd->contig_list[contig1_num].contig_left_gel,
		cd->contig_list[contig2_num].contig_left_gel,
		overlap->length);
			    
	/* Oops.  The initial coordinates in list_alignment are
	   padded, but then the alignment is depadded.  Hopefully
	   no-one will notice ! */

#if 0
	overlap->seq1_out[overlap->right+1] = '\0';
	overlap->seq2_out[overlap->right+1] = '\0';

	if (overlap->length <= cd->max_alignment) {
	    ret = list_alignment(&overlap->seq1_out[overlap->left],
				 &overlap->seq2_out[overlap->left],
				 name1,name2,seq1_start_r,seq2_start_r,buf);
	} else {
	    vmessage("%s\n", buf);
	    vmessage(" Percentage mismatch %5.1f\n\n",
		     percent_mismatch);
	}
#else
	vmessage("%s\n", buf);
	vmessage(" Percentage mismatch %5.1f\n\n",
		 percent_mismatch);
#endif

	buffij(cd->contig_list[contig1_num].contig_left_gel,
	       seq2_end_r,
	       -cd->contig_list[contig2_num].contig_left_gel,
	       seq1_start_r,
	       overlap->length, (int)overlap->score,
	       percent_mismatch);
    }

    //free_overlap(overlap);
}

int do_it_fij ( char seq[], int seq_len,
		int word_len, int min_overlap,
		double max_percent_mismatch, int compare_mode,
		int band, int gap_open, int gap_extend, double max_prob,
		int min_match, int max_alignment, int fast_mode,
		double filter_words,
		Contig_parms *contig_list, int number_of_contigs) {

    int ret, i, j, longest_diagonal;
    int max_contig;
    int seq1_len, seq2_len, contig1_num, contig2_num;
    double comp[5];
    int max_matches;
    char *depad_seq1    = NULL, *depad_seq2    = NULL;
    int  *depad_to_pad1 = NULL, *depad_to_pad2 = NULL;
    int edge_mode, job, seq1_start, seq2_start;
    int compare_method;
    Hash *h;
    OVERLAP *overlap;
    ALIGN_PARAMS *params;
    Contig_parms *contig_list_depadded = NULL;
    add_fij_t add_fij_cd;
    int one_by_one;
    int c1_loop_size;

    if (NULL == (params = create_align_params())) return -1;

    edge_mode = 10;
    seq1_start = 0;
    seq2_start = 0;
    job = RETURN_NEW_PADS | RETURN_EDIT_BUFFERS | RETURN_SEQ;

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

    max_contig = seq_len;

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
    if (NULL == (depad_seq1   = (char *)xmalloc(sizeof(char) * max_contig))||
	NULL == (depad_seq2   = (char *)xmalloc(sizeof(char) * max_contig))||
	NULL == (depad_to_pad1   = (int *)xmalloc(sizeof(int) * max_contig)) ||
	NULL == (depad_to_pad2   = (int *)xmalloc(sizeof(int) * max_contig))) {

	if (depad_seq1)      xfree(depad_seq1);
	if (depad_seq2)      xfree(depad_seq2);
	if (depad_to_pad1)   xfree(depad_to_pad1);
	if (depad_to_pad2)   xfree(depad_to_pad2);
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
	destroy_alignment_params (params);
	destroy_overlap (overlap);
	if (depad_seq1)      xfree(depad_seq1);
	if (depad_seq2)      xfree(depad_seq2);
	if (depad_to_pad1)   xfree(depad_to_pad1);
	if (depad_to_pad2)   xfree(depad_to_pad2);
	return -1;
    }
    h->fast_mode = fast_mode;
    h->filter_words =  0;
    
    if ( HASH_JOB_EXPD & compare_method ) {

        p_comp(comp,seq,seq_len);

	if(poisson_diagonals(MINMAT, longest_diagonal,
			     h->word_length, max_prob, h->expected_scores,
			     comp)) {
	    free_hash8n(h);
	    destroy_alignment_params (params);
	    destroy_overlap (overlap);
	    if (depad_seq1)      xfree(depad_seq1);
	    if (depad_seq2)      xfree(depad_seq2);
	    if (depad_to_pad1)   xfree(depad_to_pad1);
	    if (depad_to_pad2)   xfree(depad_to_pad2);
	    return -1;
	}
    }

    add_fij_cd.contig_list = contig_list;
    add_fij_cd.number_of_contigs = number_of_contigs;
    add_fij_cd.depad_to_pad1 = depad_to_pad1;
    add_fij_cd.depad_to_pad2 = depad_to_pad2;
    add_fij_cd.min_overlap = min_overlap;
    add_fij_cd.max_percent_mismatch = max_percent_mismatch;
    add_fij_cd.max_alignment = max_alignment;

    /*
     * TODO: For now we have compare_b_bulk(), but no compare_a_bulk().
     * So if we're using the sensitive method then for now we need to revert
     * back to comparing individual contigs rather than a contig against
     * the entire combined consensus.
     */
    one_by_one = (compare_method != 17 || compare_mode == COMPARE_SINGLE)
	? 1 : 0;
    add_fij_cd.one_by_one = one_by_one;

    if (one_by_one) {
	c1_loop_size = number_of_contigs;
    } else {
	c1_loop_size = 1;
    }
    
    for ( contig1_num = 0; contig1_num < c1_loop_size; contig1_num++ ) {
	if (one_by_one) {
	    seq1_len = contig_list[contig1_num].contig_end_offset 
		- contig_list[contig1_num].contig_start_offset + 1;

	    if (seq1_len < min_overlap)
		continue;

	    /* Depad seq1 */
	    overlap->seq1 = &seq[(contig_list[contig1_num].contig_start_offset)];
	    copy_seq(depad_seq1, overlap->seq1, seq1_len);
	    depad_seq(depad_seq1, &seq1_len, depad_to_pad1);
	    overlap->seq1 = h->seq1 = depad_seq1;
	    h->seq1_len = overlap->seq1_len = seq1_len;
	} else {
	    seq1_len = seq_len;

	    /* Depad seq1 */
	    overlap->seq1 = seq;
	    copy_seq(depad_seq1, overlap->seq1, seq1_len);
	    depad_seq(depad_seq1, &seq1_len, depad_to_pad1);
	    overlap->seq1 = h->seq1 = depad_seq1;
	    h->seq1_len = overlap->seq1_len = seq1_len;
	}

	/* Convert padded to depadded coords in contig_list */
	if (!one_by_one) {
	    int p, cnum;
	    contig_list_depadded = malloc(number_of_contigs *
					  sizeof(Contig_parms));
	    memcpy(contig_list_depadded,
		   contig_list,
		   number_of_contigs * sizeof(Contig_parms));


	    for (cnum = p = 0; p < seq1_len; p++) {
		if (depad_to_pad1[p] == contig_list[cnum].contig_start_offset) 
		    contig_list_depadded[cnum].contig_start_offset = p;

		if (depad_to_pad1[p] >= contig_list[cnum].contig_end_offset)
		    contig_list_depadded[cnum++].contig_end_offset = p;
	    }

	    add_fij_cd.contig_list = contig_list_depadded;
	} else {
	    add_fij_cd.contig_list = contig_list;
	}

	vmessage("Hashing sequences... ");
	//UpdateTextOutput();
	if (hash_seqn(h, 1)) {
	    verror(ERR_WARN, "find internal joins",
		   "hashing 1st sequence (#%"PRIrec")",
		   contig_list[contig1_num].contig_left_gel); 
	    continue;
	}
	vmessage("done\n\n");
	//UpdateTextOutput();
    
	(void) store_hashn ( h );

	if (filter_words > 0) {
	    /* Compute expected mean number of hits per hash key */
	    double m = (seq1_len - 20*number_of_contigs)
		/ pow(4, h->word_length);
	    h->filter_words = filter_words * m + 0.5;
	    if (h->filter_words < 2)
		h->filter_words = 2;
	} else if (filter_words < 0) {
	    h->filter_words = -filter_words;
	    if (h->filter_words < 2)
		h->filter_words = 2;
	} else {
	    h->filter_words = 0;
	}
	if (h->filter_words)
	    vmessage("Filtering words occuring more than %d times.\n",
		     h->filter_words);

	for (contig2_num = one_by_one ? contig1_num+1 : 0;
	     contig2_num < number_of_contigs; 
	     contig2_num++) {

	    if (one_by_one) {
		if (contig_list[contig2_num].contig_number ==
		    contig_list[contig1_num].contig_number)
		    continue;
	    }

	    //UpdateTextOutput();

	    seq2_len = contig_list[contig2_num].contig_end_offset 
		- contig_list[contig2_num].contig_start_offset + 1;

	    if (seq2_len < min_overlap)
		continue;

	    /* --- FWD strand --- */
	    /* Depad seq2 */
	    overlap->seq2 =
		&seq[(contig_list[contig2_num].contig_start_offset)];
	    copy_seq(depad_seq2, overlap->seq2, seq2_len);
	    depad_seq(depad_seq2, &seq2_len, depad_to_pad2);
	    overlap->seq2 = h->seq2 = depad_seq2;
	    h->seq2_len = overlap->seq2_len = seq2_len;

	    fflush(stdout);
	    if (hash_seqn(h, 2)) {
		verror(ERR_WARN, "find internal joins",
		       "hashing 2nd sequence (#%"PRIrec")",
		       contig_list[contig2_num].contig_left_gel); 
		continue;
	    }

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
		if (ret) {
		    add_fij_overlap(overlap, contig1_num, contig2_num,
				      &add_fij_cd);
		}
	    }
	    if ( ret < 0 || compare_method == 17 ) {
		if (one_by_one) {
		    ret = compare_b ( h, params, overlap );

		    if ( ret < 0 ) {
			verror(ERR_WARN, "find internal joins", "hashing" ); 
			continue;
		    }
		    
		    if ( ret ) {
			add_fij_overlap(overlap, contig1_num, contig2_num,
					&add_fij_cd);
		    }
		} else {
		    ret = compare_b_bulk ( h, params, overlap,
					   contig2_num,
					   contig_list_depadded,
					   number_of_contigs,
					   add_fij_overlap,
					   &add_fij_cd);
		}
	    }
	    free_overlap(overlap);

	    /* --- REV strand --- */
	    seq2_len = contig_list[contig2_num].contig_end_offset 
		- contig_list[contig2_num].contig_start_offset + 1;

	    /* Complement and depad seq2 */
	    overlap->seq2 =
		&seq[(contig_list[contig2_num].contig_start_offset)];
	    copy_seq(depad_seq2, overlap->seq2, seq2_len);
	    complement_seq(depad_seq2, seq2_len);
	    depad_seq(depad_seq2, &seq2_len, depad_to_pad2);
	    overlap->seq2 = h->seq2 = depad_seq2;
	    h->seq2_len = overlap->seq2_len = seq2_len;

	    fflush(stdout);
	    if (hash_seqn(h, 2)) {
		verror(ERR_WARN, "find internal joins",
		       "hashing 2nd sequence (#%"PRIrec")",
		       contig_list[contig2_num].contig_left_gel); 
		continue;
	    }

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
		if (one_by_one) {
		    ret = compare_b ( h, params, overlap );

		    if ( ret < 0 ) {
			verror(ERR_WARN, "find internal joins", "hashing" ); 
			continue;
		    }
		    
		    if ( ret ) {
			add_fij_overlap_r(overlap, contig1_num, contig2_num,
					  &add_fij_cd);
		    }
		} else {
		    add_fij_cd.seq2_len = seq2_len;
		    ret = compare_b_bulk ( h, params, overlap,
					   contig2_num,
					   contig_list_depadded,
					   number_of_contigs,
					   add_fij_overlap_r,
					   &add_fij_cd);
		}
	    }
	    free_overlap(overlap);
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
    if (contig_list_depadded)
	xfree(contig_list_depadded);

    return 0;
}
