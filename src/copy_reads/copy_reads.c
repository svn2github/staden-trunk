#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "gap_globals.h"
#include "copy_reads_globals.h"
#include "misc.h" 
#include "dna_utils.h"
#include "align_lib.h"
#include "fij.h"
#include "assemble_direct.h"
#include "seqInfo.h"
#include "tagUtils.h"
#include "complement.h"
#include "copy_reads.h"
#include "gap_cli_arg.h"
#include "active_tags.h"
#include "list_proc.h"
#include "hash_lib.h"


#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif
extern DLL_IMPORT int unknown_char;



/* #define DEBUG 1 */

int create_si_struct(GapIO *io,
		     Tcl_Interp *interp,
		     GReadings r,
		     int r_num,
		     SeqInfo **si_p)
{
    char *buf = NULL;
    int1 *conf = NULL;
    int2 *opos = NULL;
    char *seq;
    SeqInfo *si;
    char tag_str[1050];
    char comment[100];
    char *tag_name, *tag_text;
    int buf_size;
    int opos_size;

    /* The largest thing in buf is opos array */
    buf_size = r.length + 1;
    opos_size = r.length * 11 + 1;
    if (NULL == (buf = (char *)xmalloc(buf_size))) 
	return -1;
    
    /* need to fill out an SeqInfo structure for reading */
    if ((si = allocSeqInfo()) == NULL) {
	return -1;
    }
    si->e = exp_create_info();
    
    /* +1 will allow it to be NULL terminated by TextRead */
    seq  = (char *)xmalloc((r.length+1) * sizeof(*seq));
    TextRead(io, r.sequence, seq, r.length+1);
    
    strcpy(buf, io_rname(io, r_num));
    
    exp_set_entry(si->e, EFLT_ID, buf);
    exp_set_entry(si->e, EFLT_EN, buf);

    if (r.trace_name) {
	  TextRead(io, r.trace_name, buf, buf_size);
	  exp_set_entry(si->e, EFLT_LN, buf);
    }
    if (r.trace_type) {
	  TextRead(io, r.trace_type, buf, buf_size);
	  exp_set_entry(si->e, EFLT_LT, buf);
    }

    if (r.confidence) {
	if (NULL == (conf = (int1 *) xmalloc(sizeof(int1) * opos_size))) {
	    return -1;
	}
	
	DataRead(io, r.confidence, conf, r.length * sizeof(*conf), 
		 sizeof(*conf));
    }

    if (r.orig_positions) {
	if (NULL == (opos = (int2 *) xmalloc(sizeof(int2) * opos_size))) {
	    return -1;
	}
	DataRead(io, r.orig_positions, opos, r.length * sizeof(*opos), 
		 sizeof(*opos));
    }
    
    si->length = r.length;

#ifdef DEBUG
    printf("si struct r.sense %d\n", r.sense);
#endif
    if (r.sense == GAP_SENSE_ORIGINAL) {

	si->start = r.start;
	si->end = r.end;
	exp_set_entry(si->e, EFLT_SQ, seq);
    } else {
	int s, e, l;

	s = r.start;
	e = r.end;
	l = r.length;
	io_complement_seq(&l, &s, &e, seq, conf, opos);
	
	exp_set_entry(si->e, EFLT_SQ, seq);
	si->start = s;
	si->end = e;
    }

    tag_text = get_default_string(interp, copy_reads_defs, 
				  "COPY_READS.TAG_DESTINATION.TEXT");
    sprintf(comment, "%s %s", tag_text, io_name(io));
    tag_name = get_default_string(interp, copy_reads_defs, 
				  "COPY_READS.TAG_DESTINATION.NAME");

    values2tag(tag_str, tag_name, si->start+1, si->end-1, r.strand, comment);
    
    exp_set_entry(si->e, EFLT_TG, tag_str);

    si->confidence = conf;
    si->origpos = opos;
    *si_p = si;

    xfree(buf);
    xfree(seq);
    return 0;
}

double calc_average_read_quality(GapIO *io, 
				 GReadings r,
				 int1 *conf)
{
    int i;
    int conf_sum = 0;

    if (!r.confidence) {
	return 0;
    }
    DataRead(io, r.confidence, conf, r.length * sizeof(*conf), 
	     sizeof(*conf));

    /* only do visible part of reading */
    for (i = r.start; i < r.start + r.sequence_length; i++) {
	conf_sum += conf[i];
    }
    return (conf_sum / r.sequence_length);
}

int copy_reads(Tcl_Interp *interp,
	       GapIO *io_from,
	       Contig_parms contig_from,
	       GapIO *io_to,
	       Contig_parms contig_to,
	       int start,
	       int end,
	       int pos_to,
	       int complement,
	       double max_mism,
	       double min_average_qual,
	       int display,
	       Tcl_DString *copied_reads)
{
    int i;
    GReadings reading, rl;
    char *contig_name_to;
    consen_info *ci = NULL;
    align_info *ai;
    int ierr;
    int old_clen;
    int rnum;
    int name;
    int offset = 0;
    int c_num;
    SeqInfo *si = NULL;
    int direction = 0;
    int first_read = 1;
    int stop = 0;
    double average_qual;
    int1 *conf = NULL;
    char tag_str[1050];
    char comment[100];
    char *tag_name, *tag_text;
    int conf_alloc = 10000; /* realloced as needed */

    int tol;
    int maxpads = 25; /* hard coded in assemble_direct.c */
    
    contig_name_to = get_contig_name(io_to, contig_to.contig_number);

#ifdef DEBUG
    printf("READ_RAID contig to number %d %d %d %d\n", contig_to.contig_number,
	   contig_to.contig_left_gel, io_clength(io_to, contig_to.contig_number), complement);
    printf("start %d end %d pos_to %d\n", start, end, pos_to);
#endif

    /* Calculate consensus */
    if (NULL == (ci = all_consensus(io_to, gap4_global_get_consensus_cutoff())))
	return -1;

    /* 
     * tolerance is equal to the max pads that can be entered into an
     * alignment 
     */
    tol = 2*maxpads;
    /*
     * not very efficient but makes everything so much easier!
     * NB requires the FROM db to have write permissions
     */
    if (complement == GAP_STRAND_REVERSE) {
	complement_contig(io_from, contig_from.contig_number);
    }
    complement = GAP_STRAND_FORWARD;

    if (NULL == (conf = (int1 *)xmalloc(conf_alloc)))
	return -1;

    for (i = io_clnbr(io_from, contig_from.contig_number); i; 
	 i = io_rnbr(io_from, i)) {
	gel_read(io_from, i, reading);

	/* Realloc conf if necessary */
	if (reading.length > conf_alloc) {
	    int1 *conf2;
	    conf_alloc = reading.length + 100;
	    conf2 = (int1 *)xrealloc(conf, conf_alloc);
	    if (!conf2) {
		xfree(conf);
		return -1;
	    }
	    conf = conf2;
	}


	vmessage("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
	vmessage("Doing reading %s (#%d)\n", io_rname(io_from, i), i);
#ifdef DEBUG
	printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
	printf("Doing reading %s (#%d)\n", io_rname(io_from, i), i);
#endif
	/*
	 *  reject if reading is already in TO database, which will happen
	 * if a reading is entered and then you do the complement
	 */
	if (get_gel_num(io_to, io_rname(io_from, i), GGN_NAME) > 0) {
	    vmessage("  Reading is already in the database\n");
	    if (first_read) {
		continue;
	    } else {
		stop = 1;
	    }
	}

	/* 
	 * reject if reading start is beyond the end of the overlap or
	 * reading end is before the start of the overlap
	 */
	if ((reading.position > end) || ((reading.position + reading.sequence_length) < start)) { 
#ifdef DEBUG
	    printf("reading not in overlap region pos %d end %d len %d start %d \n", reading.position, end, reading.sequence_length, start);
#endif
	    vmessage("  Reading is not in the overlap region\n");
#ifdef DEBUG
	    vmessage("  start %d end %d pos %d seq_len %d\n", start, end,
		     reading.position, reading.sequence_length);
	    printf("start %d end %d pos %d seq_len %d\n", start, end,
		     reading.position, reading.sequence_length);
#endif
	    if (first_read) {
		continue;
	    } else {
		stop = 1;
	    }
	} 

	/* 
	 * reject if average quality is < min_average_quality
	 */
	if (!stop && reading.confidence) {
	    average_qual = calc_average_read_quality(io_from, reading, conf);
	    if (average_qual < min_average_qual) {
		vmessage("  Average quality of reading (%f) is too low\n", average_qual);
#ifdef DEBUG
		printf("  Average quality of reading (%f) is too low\n", average_qual);
#endif
		if (first_read) {
		    continue;
		} else {
		    stop = 1;
		}
	    } 
	}

	if (reading.left > 0 && !first_read) {
	    gel_read(io_from, reading.left, rl);
	    offset = reading.position - rl.position;
#ifdef DEBUG
	    printf("r.pos %d rl.pos %d offset %d\n", reading.position,
		   rl.position, offset);
#endif
	} else {
#ifdef DEBUG	    
	    printf("First read before pos_to %d start %d\n", pos_to, start);
#endif
	    pos_to = (pos_to + (reading.position - start));
	}

#ifdef DEBUG
	printf("pos_to %d offset %d pos_to %d\n", pos_to, offset, pos_to+offset);
#endif
	pos_to += offset;

	/* reading shouldn't be entered but I needed the offset to be calc */
	if (stop) {
	    stop = 0;
	    continue;
	}
#ifdef DEBUG
	/* if dealing with overlap region only */
	printf("name %s start %d seq_len %d end %d len %d position %d\n", 
	       io_rname(io_from, i), reading.start, 
	       reading.sequence_length, reading.end, reading.length,
	       reading.position);

	printf("start %d end %d pos_to %d\n", start, end, pos_to);
#endif
	/* need to store reading in it's original orientation */
	create_si_struct(io_from, interp, reading, i, &si);

	direction = reading.sense;
#ifdef DEBUG
	printf("sense %d complement %d direction %d\n", 
	       reading.sense, complement, direction);
#endif
	ai = assemble_align(io_to, si, ci, contig_to.contig_number, &pos_to, 
			    direction, tol, display, maxpads, max_mism, &ierr);
	
	name = get_gel_num(io_to, contig_name_to, GGN_NAME);
	if (NULL == ai) {
	    if (ierr == 2) {
		vmessage("  Percentage mismatch is too high\n");
#ifdef DEBUG
		printf("  Percentage mismatch is too high\n");
#endif
		name = 0;
	    } else if (ierr == 3) {
		vmessage("  Reading does not overlap\n");
#ifdef DEBUG
		printf("  Reading does not overlap\n");
#endif
		name = 0;
	    } else if (ierr == 4) {
		vmessage("  Reading does not overlap "
			 "within tolerance\n");
#ifdef DEBUG
		printf("  Reading does not overlap "
			 "within tolerance\n");
#endif
		name = 0;
	    } else {
		verror(ERR_WARN, "copy reads",
		       "failed to align reading '%s'", io_rname(io_from, i));
		name = 0;
	    }
	}

	c_num = contig_to.contig_number;

	if (name != 0) {
	    old_clen = io_clength(io_to, c_num);

	    if (-1 == (rnum = enter_reading(io_to, si, direction, ai, 
					    c_num, pos_to)))
		return -1;
	    
	    if (-1 == link_reading(io_to, name, rnum, c_num, pos_to))
		return -1;
	
	    if (-1 == recalc_consensus(io_to, ci, c_num, pos_to,
				       ABS(io_length(io_to, rnum)),
				       old_clen, 
				       io_clength(io_to, c_num))) {
		verror(ERR_WARN, "copy reads",
		       "failed to recalculate consensus - quitting");
		return -1;
	    }

	    vmessage("Entered reading %s\n\n", io_rname(io_from, i));
#ifdef DEBUG
	    printf("Entered reading %s\n\n", io_rname(io_from, i));
#endif
	    Tcl_DStringAppendElement(copied_reads, io_rname(io_from, i));

	    /* add tag to FROM database */
	    tag_text = get_default_string(interp, copy_reads_defs, 
					  "COPY_READS.TAG_SOURCE.TEXT");
	    sprintf(comment, "%s %s", tag_text, io_name(io_to));
	    tag_name = get_default_string(interp, copy_reads_defs, 
					  "COPY_READS.TAG_SOURCE.NAME");

	    values2tag(tag_str, tag_name, si->start+1, si->end-1, 
		       reading.strand, comment);
	    create_tag_for_gel(io_from, i, si->length, tag_str, NULL, 0, NULL,
			       0);
    

	    /* need to reset pos_to when it becomes -ve */
	    if (pos_to <= 0) {
		pos_to = 1;
	    }

	    /* 
	     * if this was the first read to be entered, set first_read to
	     * be false
	     */
	    if (first_read)
		first_read = 0;

	} else {
	    vmessage("Didn't enter reading %s\n\n", io_rname(io_from, i));
	}
	freeSeqInfo(si);
	if (ai) {
	    xfree(ai->res);
	    xfree(ai);
	}
    }

    /*
     * complement contig back again
     */
    if (complement == GAP_STRAND_REVERSE) {
	complement_contig(io_from, contig_from.contig_number);	
    }
    
    
    if (conf)
	xfree(conf);
    
    flush2t(io_to);
    free_all_consensus(ci);
    return 0;
}

int check_cons_match(char *seq1,
		     char *seq2,
		     double max_percent_mismatch)
{
    int win_len = 100;
    int len;
    int i, j, mism = 0;
    int max_mismatch;

    len = strlen(seq1);

    if (len < win_len)
	win_len = len;

    max_mismatch = win_len * max_percent_mismatch / 100;

    for (i = 0, j = 0; i < win_len; i++, j++) {
	if (!same_char(seq1[i], seq2[j])) {
	    mism++;
	}
    }

    if (mism > max_mismatch) {

	vmessage("\nLocal mismatch of %f found at postion %d over a window length of %d\n",
		 (float)mism / win_len * 100, i-win_len+1, win_len);
	vmessage("Aborting this match\n\n");
	return -1;
    }

    do {
	mism -= !same_char(seq1[i++-win_len], seq2[j++-win_len]);
	if (i < len-2) {
	    mism += !same_char(seq1[i], seq2[j]);
	}
	if (mism > max_mismatch) {
	    vmessage("\nLocal mismatch of %f found at postion %d over a window length of %d\n",
		     (float)mism / win_len * 100, i-win_len, win_len);
	    vmessage("Aborting this match\n\n");
	    return -1;
	}
    } while (i < len - 1);

    return 0;
}

void compare_consensus(Tcl_Interp *interp,
		       char *seq,
		       Contig_parms *contig_list,
		       int number_of_contigs,		
		       GapIO *io_from,
		       Contig_parms contig_from,
		       GapIO *io_to,
		       int min_overlap,
		       double max_percent_mismatch, 
		       double align_max_mism,
		       OVERLAP	*overlap,
		       ALIGN_PARAMS *params,
		       int seq2_len,
		       char *depad_seq2,
		       int *depad_to_pad1,
		       int *depad_to_pad2,
		       Hash *h,
		       int complement,
		       double min_average_qual,
		       int display_cons,
		       int display_seq,
		       Tcl_DString *copied_reads)
{
    int contig2_num;
    int ret;
    double percent_mismatch;
    int seq1_start_f, seq2_start_f;
    int seq1_end_f, seq2_end_f;
    static char buf[80],name1[10],name2[10];
    int overlap_len;

    for ( contig2_num = 0; contig2_num < number_of_contigs; 
	  contig2_num++ ) {
	vmessage("Comparing source contig %s (#%d) with destination contig %s (#%d)\n",
		 io_rname(io_from, contig_from.contig_left_gel),
		 contig_from.contig_left_gel,
		 io_rname(io_to, contig_list[contig2_num].contig_left_gel),
		 contig_list[contig2_num].contig_left_gel);

	seq2_len = contig_list[contig2_num].contig_end_offset 
	    - contig_list[contig2_num].contig_start_offset + 1;

	if ( seq2_len >= min_overlap ) {
	    overlap->seq2 = &seq[(contig_list[contig2_num].contig_start_offset)];
	    /* depad seq2.  This isn't really very efficient as we will
	     * depad the sequences lots of times, but it uses less memory 
	     */
	    
	    copy_seq(depad_seq2, overlap->seq2, seq2_len);
	    depad_seq(depad_seq2, &seq2_len, depad_to_pad2);
	    overlap->seq2 = h->seq2 = depad_seq2;
	    h->seq2_len = overlap->seq2_len = seq2_len;
	    
	    if (hash_seqn(h, 2)) {
		verror(ERR_WARN, "copy reads", "hashing 2" ); 
		continue;
	    }

	    /* always use compare_method 17 (quick method of fij) */
	    ret = compare_b ( h, params, overlap );

	    if ( ret < 0 ) {
		verror(ERR_WARN, "copy reads", "hashing" ); 
		continue;
	    }
	    
	    if ( ret ) {
		
		percent_mismatch = 100.0 - overlap->percent;

		if ( (overlap->length >= min_overlap) && 
		     (percent_mismatch <= max_percent_mismatch )) {
		    
		    /* note conversion depadded to padded coordinates */
		    
		    seq1_start_f = depad_to_pad1[overlap->left2]
			- contig_from.contig_left_extension + 1 ; 
		    
		    /* add 1 to get base number */
		    seq2_start_f = depad_to_pad2[overlap->left1]
			- contig_list[contig2_num].contig_left_extension + 1 ; 

		    overlap_len = overlap->right - overlap->left;

		    vmessage("Overlap found at position %d of contig #%d and position %d of contig #%d of length %d\n",
			     seq1_start_f, contig_from.contig_left_gel,
			     seq2_start_f, contig_list[contig2_num].contig_left_gel, overlap_len);
		    
		    overlap->seq1_out[overlap->right+1] = '\0';
		    overlap->seq2_out[overlap->right+1] = '\0';

		    if (display_cons) {
			sprintf(name2,"%d",
				contig_list[contig2_num].contig_left_gel);
			sprintf(buf," Possible join between contig in the + sense and contig %d",
				contig_list[contig2_num].contig_left_gel);
			
			/* Oops.  The initial coordinates in list_alignment are
			   padded, but then the alignment is depadded.  Hopefully
			   no-one will notice ! */
			ret = list_alignment(&overlap->seq1_out[overlap->left],
					     &overlap->seq2_out[overlap->left],
					     name1,name2,seq1_start_f,seq2_start_f,buf);
		    }
		    
		    if (-1 == (check_cons_match(&overlap->seq1_out[overlap->left],
						&overlap->seq2_out[overlap->left],
						max_percent_mismatch)))
			{
			    continue;
			}

		    /* 
		     * must -1 from overlap_len because depad_to_pad 
		     * starts at index 0 so if overlap_len == seq_len,
		     * need to index overlap_len-1
		     */
		    seq1_end_f = depad_to_pad1[overlap_len-1+overlap->left2];
		    seq2_end_f = depad_to_pad2[overlap_len-1+overlap->left1];

#ifdef DEBUG
		    printf("left1 %d left2 %d left %d right1 %d right2 %d right %d\n", overlap->left1, overlap->left2, overlap->left, overlap->right1, overlap->right2, overlap->right);
		    
		    printf("r1 %d l1 %d r2 %d l2 %d\n", 
			   depad_to_pad1[overlap->right],
			   depad_to_pad1[overlap->left],
			   depad_to_pad2[overlap->right],
			   depad_to_pad2[overlap->left]);
		    
		    printf("end seq1 %d seq2 %d\n",
			   seq1_end_f, seq2_end_f);
		    
		    printf("overlap coords are s1 %d %d s2 %d %d\n",
			   seq1_start_f, seq1_end_f, 
			   seq2_start_f, seq2_end_f);
#endif
		    copy_reads(interp, io_from, contig_from, io_to,
			       contig_list[contig2_num],
			       seq1_start_f, 
			       seq1_end_f,
			       seq2_start_f, complement,
			       align_max_mism, min_average_qual,
			       display_seq, copied_reads);
		}
	    }
	}
	free_overlap(overlap);
    }
}

int hash_consensus( Tcl_Interp *interp,
		    char *seq,    /* all contigs in TO db */
		    int seq_len,
		    char *seq1,  /* single contig in FROM db */
		    int seq1_len_orig,
		    int word_len, 
		    int min_overlap,
		    double max_percent_mismatch, 
		    int band, 
		    int gap_open, 
		    int gap_extend, 
		    int min_match,
		    double align_max_mism,
		    double min_average_qual,
		    int display_cons,
		    int display_seq,
		    Contig_parms *contig_list, 
		    int number_of_contigs,
		    GapIO *io_from,
		    Contig_parms contig_from,
		    GapIO *io_to,
		    Tcl_DString *copied_reads) {

    int i, longest_diagonal;
    int max_contig;
    int seq2_len;
    int max_seq;
    int max_matches;
    char *depad_seq1,    *depad_seq2;
    int  *depad_to_pad1, *depad_to_pad2;
    int edge_mode, job, seq1_start, seq2_start;
    int compare_method;
    int seq1_len = seq1_len_orig;

    Hash *h;
    
    OVERLAP	*overlap;
    ALIGN_PARAMS *params;

    if (NULL == (params = create_align_params())) return -1;

    edge_mode = 10;
    seq1_start = 0;
    seq2_start = 0;
    seq2_len = 0;  /* sent to compare_consensus but reset there before use! */
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
    init_overlap (overlap, seq1, seq, seq1_len, seq_len);

    /* 
     * find longest contig 
     * initialise hashing
     * compare each contig in each orientation
     */
    for (i = 0, max_contig = 0; i < number_of_contigs; i++) {
	if (contig_list[i].contig_end_offset
	    - contig_list[i].contig_start_offset > max_contig ) {
	    max_contig = contig_list[i].contig_end_offset
		       - contig_list[i].contig_start_offset;
	}
    }
    max_contig += 1;

    /* need to check against single comparing contig */
    if (max_contig < seq1_len) 
	max_contig = seq1_len;

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

    /* always use the quick method of fij */
    compare_method = 17;

    if ( init_hash8n (max_contig, longest_diagonal,
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
    
    if ( seq1_len >= min_overlap ) {

	/* Depad seq1 */
	
	overlap->seq1 = seq1;
	copy_seq(depad_seq1, overlap->seq1, seq1_len);
	depad_seq(depad_seq1, &seq1_len, depad_to_pad1);
	overlap->seq1 = h->seq1 = depad_seq1;
	h->seq1_len = overlap->seq1_len = seq1_len;

	if (hash_seqn(h, 1)) {
	    verror(ERR_WARN, "copy reads", "hashing 1" ); 
	    return -1;
	}
	
	(void) store_hashn ( h );

	compare_consensus(interp, seq, contig_list, number_of_contigs, io_from, 
			  contig_from, io_to, min_overlap, 
			  max_percent_mismatch, align_max_mism, overlap, 
			  params, seq2_len, depad_seq2, depad_to_pad1, 
			  depad_to_pad2, h, GAP_STRAND_FORWARD, 
			  min_average_qual, display_cons, display_seq,
			  copied_reads);
    }

    /* now do complementary strand */
    vmessage("Complementing contig %d\n", contig_from.contig_left_gel);
    seq1_len = seq1_len_orig;
    if ( seq1_len >= min_overlap ) {

	/* Depad seq1 */
	
	overlap->seq1 = seq1;
	copy_seq(depad_seq1, overlap->seq1, seq1_len);
	complement_seq(depad_seq1, seq1_len);
	depad_seq(depad_seq1, &seq1_len, depad_to_pad1);
	overlap->seq1 = h->seq1 = depad_seq1;
	h->seq1_len = overlap->seq1_len = seq1_len;

	if (hash_seqn(h, 1)) {
	    verror(ERR_WARN, "copy reads", "hashing 1" ); 
	    return -1;
	}
    
	(void) store_hashn ( h );
	compare_consensus(interp, seq, contig_list, number_of_contigs, io_from, 
			  contig_from, io_to, min_overlap, 
			  max_percent_mismatch, align_max_mism, overlap, 
			  params, seq2_len, depad_seq2, depad_to_pad1, 
			  depad_to_pad2, h, GAP_STRAND_REVERSE, 
			  min_average_qual, display_cons, display_seq,
			  copied_reads);
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

int init_copy_reads(Tcl_Interp *interp,
		    GapIO *io_from,
		    GapIO *io_to,
		    int compare_mode, /* could be sorted out to save work below*/
		    int mask,
		    int min_overlap,
		    double max_percent_mismatch,
		    int word_len,
		    int min_match,
		    int band,
		    double align_max_mism,
		    double min_average_qual,
		    int display_cons,
		    int display_seq,
		    int min_contig_len,
		    int num_contigs_from,
		    contig_list_t *contig_array_from,
		    int num_contigs_to,
		    contig_list_t *contig_array_to,
		    Tcl_DString *copied_reads)
{
    char *consensus_from, *consensus_to, *consensus;
    int i;
    int ret;
    int max_consensus;
    Contig_parms *contig_list_from, *contig_list_to;
    int database_size_from, database_size_to;
    int task_mask;
    int consensus_length_from, consensus_length_to;
    int max_read_length_from, max_read_length_to;
    Hidden_params p;
    int contig_len;

    /* FIXME: get these from elsewhere */

    int gap_open, gap_extend;

#if 0
    gap_open = gap4_global_get_gopenval();
    gap_extend = gap4_global_get_gextendval();
#endif

    /* 
     * hard code these since we are not using the Millers and Myers table 
     * but using the tables/align_lib_nuc_matrix 
     */
    gap_open = 12;
    gap_extend= 4;

    p.min = p.max = p.verbose = 0;
    p.do_it = 0;
    p.use_conf = 1;
    p.test_mode = 0;
    p.start = 0;
    p.lwin1 = 0;
    p.lcnt1 = 0;
    p.rwin1 = 0;
    p.rcnt1 = 0;
    p.qual_val = 0;
    p.window_len = 0;
    p.gap_open = gap4_global_get_gopenval();
    p.gap_extend = gap4_global_get_gextendval();
    p.band = 30; /* FIXME: hardwired 30 bases band for aligning hidden data */
    max_consensus = gap4_global_get_maxseq();

    max_read_length_from = max_gel_len(io_from);
    max_read_length_to = max_gel_len(io_to);

    if ((consensus_from = (char *)xmalloc(max_consensus * sizeof(char)))==NULL){
	return(-1);
    }
    if ((consensus_to = (char *)xmalloc(max_consensus * sizeof(char)))==NULL){
	return(-1);
    }
    if ((consensus = (char *)xmalloc(max_consensus * sizeof(char)))==NULL){
	return(-1);
    }

    database_size_from = io_dbsize(io_from);
    if ( ! (contig_list_from = get_contig_list (database_size_from, io_from,
						num_contigs_from, 
						contig_array_from))) {
	xfree (consensus_from);
	return -5;
    }

    database_size_to = io_dbsize(io_to);
    if ( ! (contig_list_to = get_contig_list (database_size_to, io_to,
					      num_contigs_to, 
					      contig_array_to))) {
	xfree (consensus_to);
	return -5;
    }

    task_mask = ADDTITLE | NORMALCONSENSUS;
    if ( mask == 2 ) task_mask |= MARKING; 
    if ( mask == 3 ) task_mask |= MASKING;

    consensus_length_from = 0;
    if ( ret = make_consensus ( task_mask, io_from, consensus_from, NULL,
				contig_list_from, num_contigs_from,
				&consensus_length_from,
				max_read_length_from,
				max_consensus,
				p,
				gap4_global_get_consensus_cutoff() ) ) {

	xfree ( consensus_from );
	xfree(contig_list_from);
	return -1;
    }
    /*
     * loop through the FROM database, comparing each contig with all contigs
     * in the TO database. Use the COMPARE_SINGLE mode of fij which requires
     * that the first element in the contig array should be the single contig
     * to compare
     */

    for (i = 0; i < num_contigs_from; i++) {
	contig_len = contig_list_from[i].contig_end_offset
	    - contig_list_from[i].contig_start_offset + 1;

	/* only use contigs >= min_contig_len */
	if (contig_len < min_contig_len) {
	    continue;
	}

	strncpy(consensus, 
		&consensus_from[contig_list_from[i].contig_start_offset], 
		contig_len);
#ifdef DEBUG
	printf("contig from number %d start offset %d end %d len %d\n\n\n", 
	       contig_list_from[i].contig_number, 
	       contig_list_from[i].contig_start_offset,
	       contig_list_from[i].contig_end_offset, contig_len);
#endif
	/* 
	 * need to do this each time because adding readings to the TO db
	 * will change it's length and sequence
	 */
	consensus_length_to = 0;
	if ( ret = make_consensus ( task_mask, io_to, consensus_to, NULL,
				    contig_list_to, num_contigs_to,
				    &consensus_length_to,
				    max_read_length_to,
				    max_consensus,
				    p,
				    gap4_global_get_consensus_cutoff() ) ) {

	    xfree ( consensus_to );
	    xfree(contig_list_to);
	    return -1;
	}

	if ( ret = hash_consensus (interp, consensus_to, consensus_length_to,
				   consensus, contig_len,
				   word_len, min_overlap, 
				   max_percent_mismatch, band, gap_open, 
				   gap_extend, min_match,
				   align_max_mism, min_average_qual,
				   display_cons, display_seq,
				   contig_list_to, num_contigs_to,
				   io_from, contig_list_from[i], io_to,
				   copied_reads)) {
	    
	    xfree(contig_list_from);
	    xfree (consensus_from);
	    xfree(contig_list_to);
	    xfree (consensus_to);
	    xfree(consensus);
	    return -1;
	}
    }

    xfree(contig_list_from);
    xfree (consensus_from);
    xfree(contig_list_to);
    xfree (consensus_to);
    xfree(consensus);

    return 0;
}


int tcl_copy_reads(ClientData clientData,
		   Tcl_Interp *interp, 
		   int objc, 
		   Tcl_Obj *CONST objv[])
{
    contig_list_t *contig_array_from = NULL;
    contig_list_t *contig_array_to = NULL;
    int num_contigs_from = 0;
    int num_contigs_to = 0;
    GapIO *io_from;
    GapIO *io_to;
    int mode = 0, mask = 0;
    Tcl_DString copied_reads;

    typedef struct {
	int handle_from;
	int handle_to;
	int word_len;
	int min_overlap;
	float max_mis;
	char *mode;
	char *inlist_from;
	char *inlist_to;
	int band;
	int min_match;
	char *mask;
	char *tag_list;
	float align_max_mism;
	int min_contig_len;
	float min_average_qual;
	int display_cons;
	int display_seq;
    } raid_arg;

    raid_arg args;
    cli_args a[] = {
	{"-io_from",	  ARG_INT,  1, NULL,      offsetof(raid_arg, handle_from)},
	{"-io_to",	  ARG_INT,  1, NULL,      offsetof(raid_arg, handle_to)},
	{"-word_length",  ARG_INT,  1, "8",	  offsetof(raid_arg, word_len)},
	{"-min_overlap",  ARG_INT,  1, "20",   offsetof(raid_arg, min_overlap)},
	{"-max_pmismatch",ARG_FLOAT,1, "30.0",    offsetof(raid_arg, max_mis)},
	{"-contigs_from", ARG_STR,  1, NULL,      offsetof(raid_arg, inlist_from)},
	{"-contigs_to",	  ARG_STR,  1, NULL,      offsetof(raid_arg, inlist_to)},
	{"-band",         ARG_INT,  1, "1",      offsetof(raid_arg, band)},
	{"-min_match",    ARG_INT,  1, "20",	 offsetof(raid_arg, min_match)},
	{"-mask",	  ARG_STR,  1, "none",    offsetof(raid_arg, mask)},
	{"-tag_types",	  ARG_STR,  1, "",        offsetof(raid_arg, tag_list)},
	{"-align_max_mism", ARG_FLOAT,  1, "10.0",offsetof(raid_arg, align_max_mism)},
	{"-min_contig_len", ARG_INT,  1, "2000",offsetof(raid_arg, min_contig_len)},
	{"-min_average_qual", ARG_FLOAT,  1, "30.0",offsetof(raid_arg, min_average_qual)},
	{"-display_cons", ARG_INT,  1, "0",offsetof(raid_arg, display_cons)},
	{"-display_seq",  ARG_INT,  1, "0",offsetof(raid_arg, display_seq)},
	{NULL,	0,	0, NULL, 0}
    };

    vfuncheader("copy reads");

    if (-1 == gap_parse_obj_args(a, &args, objc, objv)) {
	vmessage("Error in parsing arguments\n");
	return TCL_ERROR;
    }
    /* Parse mode and mask */
    if (strcmp(args.mask, "none") == 0)
	mask = 1;
    else if (strcmp(args.mask, "mark") == 0)
	mask = 2;
    else if (strcmp(args.mask, "mask") == 0)
	mask = 3;
    else {
	Tcl_SetResult(interp, "invalid mask mode", TCL_STATIC);
	return TCL_ERROR;
    }

    /* set this always */
    mode = COMPARE_SINGLE;

    if ( (io_from = io_handle(&args.handle_from)) == NULL){	
	verror(ERR_FATAL, "copy_reads", "invalid io handle %d", args.handle_from);
	return TCL_OK;
    }
    if ( (io_to = io_handle(&args.handle_to)) == NULL){	
	verror(ERR_FATAL, "copy_reads", "invalid io handle");
	return TCL_OK;
    }

    /* create contig name array */
    active_list_contigs(io_from, args.inlist_from, &num_contigs_from, &contig_array_from);
    active_list_contigs(io_to, args.inlist_to, &num_contigs_to, &contig_array_to);

    if (SetActiveTags(args.tag_list) == -1) {
	return TCL_OK;	
    } 

    Tcl_DStringInit(&copied_reads);

    if (init_copy_reads(interp, io_from, io_to, mode, mask, args.min_overlap, 
			(double)args.max_mis, args.word_len, 
			args.min_match, args.band,
			(double)args.align_max_mism, 
			(double)args.min_average_qual,
			args.display_cons, args.display_seq,
			args.min_contig_len,
			num_contigs_from, contig_array_from,
			num_contigs_to, contig_array_to, &copied_reads) < 0 ) {
	verror(ERR_WARN, "copy reads", "Failure in Copy Reads");
	SetActiveTags("");
	return TCL_OK;
    }

    xfree(contig_array_from);
    xfree(contig_array_to);
    SetActiveTags("");	
    Tcl_DStringResult(interp, &copied_reads);

    return TCL_OK;

}

/*
 * ---------------------------------------------------------------------------
 * Tcl initialisation and top-level registered functions
 * ---------------------------------------------------------------------------
 */

/*
 * This is called when the library is dynamically linked in with the calling
 * program. Use it to initialise any tables and to register the necessary
 * commands.
 */
int Copy_reads_Init(Tcl_Interp *interp) {
    init_globals(interp);
    init_copy_reads_globals(interp);

    if (NULL == Tcl_CreateObjCommand(interp, "copy_reads", tcl_copy_reads,
				     (ClientData)NULL,
				     (Tcl_CmdDeleteProc *) NULL))
	return TCL_ERROR;

    return TCL_OK;
}

int Copy_reads_SafeInit(Tcl_Interp *interp) {
    return Copy_reads_Init(interp);
}

int Copy_reads_Unload(Tcl_Interp *interp, int flags) {
    Tcl_SetResult(interp, "Pkg_Unload() function not implemented",
		  TCL_STATIC);
    return TCL_ERROR;
}

int Copy_reads_SafeUnload(Tcl_Interp *interp, int flags) {
    Tcl_SetResult(interp, "Pkg_SafeUnload() function not implemented",
		  TCL_STATIC);
    return TCL_ERROR;
}
