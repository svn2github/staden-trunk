/* routine to handle all consensus calculations for the assembly program

   There are several different types of "consensus" required:

   1. all contigs or one contig
   2. contigs sorted on length
   3. contigs with added hidden data
   4. contigs masked
   5. single stranded
   6. normal consensus
   7. quality codes
   8. with title
   9. without an added title

   Do them all in one loop by putting the contig numbers in a list first.

   Tell consen what to do by using variable "task_mask" which has the
   appropriate bits from the list above set.

   add_contig_title adds a contig title to an array of chars
   get_hidden gets some hidden data to stick on the end of a contig
   mask_contig sent the consensususes returns it masked
   sort_contigs sorts a list of contigs on length

*/

#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <fcntl.h>
#include <unistd.h>
#include "io_utils.h"
#include <io_lib/expFileIO.h>
#include "misc.h"
#include "consen.h"
#include "consensus.h"
#include "gap_globals.h"
#include "dna_utils.h"
#include "list_proc.h"
#include "extract.h"
#include "align_lib.h"
#include "gap4_compat.h"
#include "qual.h"

/*#define MAXGEL 1024*/
#define MAXGEL_PLUS 1024

/************************************************************/

static char standard_to_masked[256];
static char standard_to_marked[256];
static char marked_to_masked[256];
static char masked_to_marked[256];


/************************************************************/

void set_mask_lookup() {

/* 	set up table of values for masking and marking base characters */


    int i;
    for (i=0;i<256;i++) {
	standard_to_masked[i] = i;
	standard_to_marked[i] = i;
	masked_to_marked[i]   = i;
	marked_to_masked[i]   = i;
    }

    standard_to_masked['A'] = 'd';
    standard_to_masked['C'] = 'e';
    standard_to_masked['G'] = 'f';
    standard_to_masked['T'] = 'i';

    standard_to_marked['A'] = 'a';
    standard_to_marked['C'] = 'c';
    standard_to_marked['G'] = 'g';
    standard_to_marked['T'] = 't';

    masked_to_marked['d'] = 'a';
    masked_to_marked['e'] = 'c';
    masked_to_marked['f'] = 'g';
    masked_to_marked['i'] = 't';

    marked_to_masked['a'] = 'd';
    marked_to_masked['c'] = 'e';
    marked_to_masked['g'] = 'f';
    marked_to_masked['t'] = 'i';

}


/*
 * Add title to a consensus. 
 *
 * Inputs: 
 * 
 * 1. consensus:     20 char chunk of consensus array
 * 2. project_name:  project name terminated by .
 * 3. left_gelnumber: the contig left gel number
 *
 * Output:
 *
 * 1. Title is 20 chars in length of form: 
 *
 * <project_name.left_gelnumber-->
 */
void add_contig_title(char *consensus, char *project_name,
		      tg_rec left_gelnumber)
{

    int  plen, rlen;
    static char dashes[]="----------------";
    char buf[50];
    char *cp;

    /* Find out how many digits for reading name */
    rlen = sprintf(buf, "%"PRIrec, left_gelnumber);

    /* Length of project name */
    if (cp = strchr(project_name, '.'))
	plen = cp - project_name;
    else
	plen = strlen(project_name);
    if (plen + rlen + 3 /* < . > */ > 20) {
	plen = 20 - rlen - 3;
    }
    
    sprintf(consensus, "<%.*s.%.*"PRIrec"%.*s>",
	    plen, project_name,
	    rlen, left_gelnumber,
	    17-plen-rlen, dashes);
}

/*int comparecontigs ( int *l1, int *l2 );*/

/************************************************************/


int comparecontigs ( const void *v1, const void *v2 ) {


    Contig_parms *c1, *c2;
    int l1, l2;

    c1 = (Contig_parms *)v1;
    c2 = (Contig_parms *)v2;
    l1 = c1->contig_end - c1->contig_start;
    l2 = c2->contig_end - c2->contig_start;


    return l2 - l1;

}


/************************************************************/


int sort_contigs ( Contig_parms *contig_list, 
		  int number_of_contigs) {

/*	sort contigs on length

	input:  a list of contig numbers.
	output: the sorted list.
*/


    qsort ( ( void * ) contig_list, number_of_contigs, 
	       sizeof (Contig_parms), comparecontigs );

    return 0;
}




/************************************************************/
/*
 * Given an array of contig numbers (1..Num) or a single contig
 * ((8000-1)..(8000-Num)), generate a Contig_parms array containing
 * contig number, start, end, and left gel.
 */
Contig_parms *get_contig_list (int database_size, GapIO *io,
			       int number_of_contigs,
			       contig_list_t *contig_array) {
    int i;
    Contig_parms *contig_list;

    if (!contig_array)
	number_of_contigs = NumContigs(io);

    if (0 == number_of_contigs)
	return NULL;

    if ( ! ( contig_list = (Contig_parms *)xmalloc(sizeof(Contig_parms) *
						   number_of_contigs) )) {
	return NULL;
    }

    /* 30-12-97 rs change for new fij: initialise extensions to 0 */
    /* 17-12-98 rs change for find_repeats bug: add offsets to structure */

    /* if contig_array is NULL, all contigs are used to make the list */
    for ( i=0;i<number_of_contigs;i++) {
	if (contig_array) {
	    contig_list[i].contig_number = contig_array[i].contig;
	    contig_list[i].contig_start  = contig_array[i].start;
	    contig_list[i].contig_end    = contig_array[i].end;
	    contig_list[i].contig_left_extension  = 0;
	    contig_list[i].contig_right_extension = 0;
	    contig_list[i].contig_start_offset = 0;
	    contig_list[i].contig_end_offset = 0;
	} else {
	    contig_list[i].contig_number = arr(tg_rec, io->contig_order, i);
	    contig_list[i].contig_start  = 1;
	    contig_list[i].contig_end    = ABS(io_clength(io, i+1));
	    contig_list[i].contig_left_extension  = 0;
	    contig_list[i].contig_right_extension = 0;
	    contig_list[i].contig_start_offset = 0;
	    contig_list[i].contig_end_offset = 0;
	}
	contig_list[i].contig_left_gel =
	    io_clnbr(io, contig_list[i].contig_number);
    }

    return contig_list;
}

/* given a position in a consensus sequence return the contig left gel num
 * ASSUME that contig_start contains array element in consensus of left end 
 * of contig!! and that the contig_list is otherwise correctly filled in.
 */

int contig_listel_from_con_pos ( Contig_parms *contig_list, 
				 int number_of_contigs, int pos_in_contig ) {

    int i, hi, lo;

    if ( 0 == number_of_contigs ) return -1;
    if ( 1 == number_of_contigs ) return 0;

    lo = 0;
    hi = number_of_contigs - 1;

    while (lo < hi) {
	i = (hi + lo) / 2;
	if (pos_in_contig < contig_list[i].contig_start_offset) {
	    hi = i;
	} else if (pos_in_contig >= contig_list[i + 1].contig_start_offset) {
	    lo = i + 1;
	} else {
	    return i;
	}
    }
    /* If here, pos_in_contig is either before the first contig or after 
       the last */
    if (pos_in_contig < contig_list[0].contig_start_offset) return 0;
    return number_of_contigs - 1;
}

	
/**********************************************************************/


void maskit ( char *seq, int seq_length, int job) {

    /*	routine to do masking and marking	*/

    int i;

    switch (job) {

    case (STANDARD_TO_MASKED):


	for ( i=0;i<seq_length;i++ )
	    seq[i] = standard_to_masked[ (unsigned)seq[i] ];
	break;

    case (MASKED_TO_MARKED):

	for ( i=0;i<seq_length;i++ )
	    seq[i] = masked_to_marked[ (unsigned)seq[i] ];
	break;

    case (MARKED_TO_MASKED):

	for ( i=0;i<seq_length;i++ )
	    seq[i] = marked_to_masked[ (unsigned)seq[i] ];
	break;

    case (STANDARD_TO_MARKED):

	for ( i=0;i<seq_length;i++ )
	    seq[i] = standard_to_marked[ (unsigned)seq[i] ];
	break;

    default:

	verror(ERR_WARN, "maskit", "unknown job no. %d", job);
	break;
    }
}

/****************************************************************************/
#if 0
int mask_consensus(GapIO *io, char *consensus, int contig, int lreg, int rreg, 
		   int job) {
    GAnnotations *ap;
    GContigs c;
    GReadings r;
    int gel;
    extern char **active_tag_types;
    extern int number_of_active_tags;

    /* Routine to mask regions of a consensus	*/

    /* *consensus		the consensus
       *active_tag_types	the list of tag types
       number_of_active_tags	the number of tag types in the list
       lreg, rreg		the start and end points for the consensus
                                note consensus[0] corresponds to lreg
                                and  consensus[rreg-lreg] to rreg

       Deal with tags on the consensus (send -contig to vtagget) and on
       the individual reads. Mask_job = 1, mark_job = 2;
       Masking and marking are done by changing the bases to new character sets.
       The algorithm may mask the same bases several times if they are
       covered by several tags that appear on the list, but it is the
       simplest thing to do.
       */

    /* Is there anything to do ? */

    if ( number_of_active_tags == 0 )
	return 0;

    /* init */

    if (0 != contig_read(io, contig, c))
	return -1;

    if (!lreg)
	lreg = 1;
    if (!rreg)
	rreg= c.length;
    
    /* do the tags on readings first	*/

    gel = c.left;

    while (gel) {
	gel_read(io, gel, r);

	if ( r.position <= rreg ) {
	    /* init vtagget() */
	    ap = vtagget(io, gel, number_of_active_tags, active_tag_types);


	    while (ap && ap != (GAnnotations *)-1) {
		int e;

		/* Normalise tags if necessary */
		if (r.sense)
		    ap->position = r.length - ap->position - ap->length + 2;

		/* 100% cutoff data - reject tag */
		if (ap->position + ap->length - 1 <= r.start ||
		    ap->position >= r.end)
		    goto next;

		/* overlap cutoff with used - clip appropriately */
		if (ap->position <= r.start) {
		    ap->length   -= r.start - ap->position + 1;
		    ap->position += r.start - ap->position + 1;
		}

		e = r.position - r.start + ap->position - 1;

		if ((e + ap->length > lreg) && ( e <= rreg ) ) {
		    if (e < lreg) {
			ap->length -= lreg-e;
			e = lreg;
		    }

		    if (e <= rreg && e + ap->length - 1 > rreg) {
			ap->length = rreg - e + 1;
		    }

		    e = e - lreg + 1;
		    (void) maskit ( &consensus[e-1], ap->length, job);

		}

	    next:
		ap = vtagget(io, 0, number_of_active_tags, active_tag_types);
	    }
	}
	gel = r.right;
    }


    /* now do the tags on the consensus	*/


	/* init vtagget() */

    gel = -contig;

    ap = vtagget(io, gel, number_of_active_tags, active_tag_types);

    while (ap && ap != (GAnnotations *)-1 && ap->position <= rreg) {
	int e;

	e = ap->position;
	    
	if (e + ap->length >= lreg) {
	    if (e < lreg) {
		ap->length -= lreg-e;
		e = lreg;
	    }

	    if (e <= rreg && e + ap->length - 1 > rreg) {
		ap->length = e + ap->length - 1 - rreg;
	    }
	    
	    (void) maskit ( &consensus[e-1], ap->length, job);

	}

	ap = vtagget(io, 0, number_of_active_tags, active_tag_types);

    }
    return 0;

}
#endif

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
		  int max_unknown1) {


    int window_len, max_unknown, jstart, bad_start;

    window_len = window_len1;
    max_unknown = max_unknown1;

    jstart = MIN ( start, (int)strlen( seq ) );

    /* 	Do the first search from start on */
    bad_start = bad_data_start(&seq[jstart], window_len, max_unknown,
			       strlen(&seq[jstart]), 1);

    /* 	If required do another search from here on 
     *  if ( window_len2 ) {
     *
     *	jstart = bad_start + jstart + window_len1/2;
     *  jstart = MIN ( jstart, strlen( seq ) );
     *	bad_start = bad_data_start(&seq[jstart], window_len2, max_unknown2,
     *			   strlen(&seq[jstart]), 1);
     *
     *  }
     */

    return bad_start + jstart;
}

/*
 * Scans rightwards in a quality buffer until the average quality drops below
 * a specific threshold, defined by the global qual_val and window_len
 * parameters.
 *
 * Having found this window, the procedure repeats with successively smaller
 * windows until the exact base is identified.
 */
int scan_right(Hidden_params p, int1 *conf, int start_pos, int len) {
    int i, total, rclip;
    int lowest_total;
    int win_len = p.window_len;

    do {
	lowest_total = p.qual_val * win_len;

	total = 0;
	for (i = start_pos; i < start_pos + win_len && i < len; i++)
	    total += conf[i];

	if (i + win_len < len) {
	    i = start_pos;
	    do {
		total = total - conf[i] + conf[i+win_len];
		i++;
	    } while (i <= (len - win_len - 1) && total >= lowest_total);
	}

	start_pos = i-1;
    } while (--win_len > 0);

    //rclip = i == len ? len + 1 : i;
    rclip = i == len ? len : i + 1;
    if (p.verbose)
	printf("    right clip = %d of %d\n", rclip, len);

    return rclip;
}

#if 0
int get_hidden_seq (GapIO *io, tg_rec read_number,
		   char *hidden_seq,
		   int *length_hidden) {

    /* get the hidden sequence for read read_number in hidden_seq
       return its length in length_hidden (which is the max length
       of hidden_seq on entry). Return 1 for success 0 for error.
       */

    int len;
    int comp;

    if ( io_get_extension(io,read_number,hidden_seq,*length_hidden,
				 &len,&comp) == 0) {
	*length_hidden = len;

	if (comp) {
	    complement_seq(hidden_seq, len);
	}
	return 1;
    }
    return 0;
}
#endif

void revconf( int1 *conf, int len ) {
    register int i, middle, j;
    int1 temp;

    middle = len/2;
    for ( i = 0, j = len-1; i < middle; i++, j--) {
	temp = conf[i];
	conf[i] = conf[j];
	conf[j] = temp;
    }
}

/************************************************************/

#if 0
int get_hidden(GapIO *io, int contig, int max_read_length, 
	       int end, Hidden_params p,
	       char *hidden_seq, char *t_hidden_seq, char *consensus) {

    /* Routine to get the hidden data near to the ends of a contig
       so that we can add it to its ends.
       */

    int read_number, right_most, left_most, length_hidden;
    GContigs c;
    GReadings r;
    int i, l, new_end, length_extension;
    int start_hidden, best_length_hidden;

    int band, ret;
    OVERLAP	*overlap;
    ALIGN_PARAMS params;
    SEG Seg, *seg;
    double percent_mismatch, max_percent_mismatch;

    int1 *conf = NULL;
    int use_conf;

    if (NULL == (overlap = create_overlap())) return -1;

    seg = &Seg;

    max_percent_mismatch = 101.0;

    band = p.band;

    params.gap_open = p.gap_open;
    params.gap_extend = p.gap_extend;
    params.edge_mode = 9;
    params.job = 2;
    params.old_pad_sym = '*';
    params.new_pad_sym = '.';
    params.first_row = 0;

    use_conf = p.use_conf;

    /* set the output pointer */
    seg->seq = hidden_seq;

    /* Get the contig info	*/

    if (0 != contig_read(io, contig, c)) {
	destroy_overlap(overlap);
	return -1;
    }

    length_extension = 0;
    start_hidden = best_length_hidden = 0; /* only set to avoid warnings */

    switch ( end ) {

    case ( RIGHT_END ):

	/* set the position of the furthest right weve reached	*/

      right_most = c.length;
      read_number = c.right;

      while ( read_number ) {

	gel_read(io, read_number, r);

	/* if we are max_read_length from the end we cannot get 
	   a better extension so bail out	*/
	
	if ( right_most - r.position > max_read_length )
	  break;
	
	/* only interested in original sense reads */

	if ( ! ( r.sense ) ) {

	    length_hidden = MAXGEL_PLUS-1;
	    if ( get_hidden_seq 
		(io, read_number, t_hidden_seq,&length_hidden) 
		&& length_hidden) {

		if (use_conf) {
		    if (0 != io_aread_seq(io, read_number, NULL, NULL, NULL,
					  NULL, &conf, NULL, 0)) {
			use_conf = 0;
			if (conf) {
			    xfree(conf);
			    conf = NULL;
			}
		    }
		}
		if (use_conf) {
		    for (i = r.end; i < r.length; i++)
			if (conf[i] != 0)
			    break;
		    if (i == r.length) {
			if (p.verbose)
			    puts("    Confidence values are all zero - using sequence");
			use_conf = 0;
		    }
		}
		if (use_conf) {

		    length_hidden = MIN
			(length_hidden,scan_right
			 (p, conf, r.end, r.length)-r.end+1);
		}
		
		if (!use_conf) {

		    length_hidden = end_of_good ( t_hidden_seq, p.start,
						 p.rwin1, p.rcnt1);
		}

		/* longer extension ?	*/
	    
		new_end = r.position + r.sequence_length + length_hidden - 1;
		l = MAX ( right_most, new_end );
		
		if ( l > right_most ) {
		    right_most = l;
		    length_extension = new_end - c.length;
		    start_hidden = r.position + r.sequence_length;
		    best_length_hidden = length_hidden;
		    (void) memcpy ( hidden_seq, 
				   t_hidden_seq,
				   best_length_hidden);
		}
	    }
	}
	read_number = r.left;
	if (conf) {
	    xfree(conf);
	    conf = NULL;
	}
      }

      if (conf) {
	  xfree(conf);
	  conf = NULL;
      }


      /* when we get here if length_extension != 0 we have the
       * longest segment of hidden data in hidden_seq
       * It starts at start_hidden (-1 for element) in the consensus
       * and is of length best_length_hidden.
       */

      if ( length_extension != 0 ) {
	
	  if ( start_hidden <= c.length ) {
	      /* make seq1 the consensus, seq2 the hidden data */

	      (void) memcpy (t_hidden_seq, 
			     hidden_seq,
			     best_length_hidden);
	      init_overlap (overlap, 
	      &consensus[start_hidden-1],
	      t_hidden_seq,
	      c.length-start_hidden+1,
	      best_length_hidden);
	      /*
	      overlap.seq1 = &consensus[start_hidden-1];
	      overlap.seq1_len = c.length-start_hidden+1;
	      overlap.seq2 = t_hidden_seq;
	      overlap.seq2_len = best_length_hidden;
	      */
	      params.band = MIN(band,overlap->seq1_len);
	      params.band_left = -params.band;
	      params.band_right = params.band;

	      ret = affine_align(overlap, &params);

	      if ( ret < 0 ) {
		  destroy_overlap(overlap);
		  return -2;
	      }
	      percent_mismatch = 100.0 - overlap->percent;

	      if (percent_mismatch < max_percent_mismatch ) {
		  ret = get_segment(overlap,seg,2);
		  if ( ret < 0 ) {
		      destroy_overlap(overlap);
		      return -2;
		  }
		  length_extension = seg->length;
	      }
	      else {
		  length_extension = 0;
		  destroy_overlap(overlap);
		  return -3; 
	      }
	  }
      }
      destroy_overlap(overlap);
      return length_extension;
      break;

  case ( LEFT_END ):

      /* set the position of the furthest left weve reached	*/

      left_most = 1;
      read_number = c.left;

      while ( read_number ) {

	  gel_read(io, read_number, r);
	  
	  /* if we are max_read_length from the end we cannot get 
	     a better extension so bail out	*/
	
	  if ( r.position > max_read_length )
	      break;
	
	  /* only interested in reverse sense reads */
	  
	  if ( r.sense ) {
	      
	      length_hidden = MAXGEL_PLUS-1;
	      if ( get_hidden_seq 
		  (io, read_number, t_hidden_seq,&length_hidden) 
		  && length_hidden) {

		  if (use_conf) {
		      if (0 != io_aread_seq(io, read_number, NULL, NULL, NULL,
					    NULL, &conf, NULL, 0)) {
			  use_conf = 0;
			  if (conf) {
			      xfree(conf);
			      conf = NULL;
			  }
		      }
		  }
		  if (use_conf) {
		      /* reverse the order of the hidden confs */
		      revconf(conf,r.length);
		      for (i = r.end; i < r.length; i++)
			  if (conf[i] != 0)
			      break;
		      if (i == r.length) {
			  if (p.verbose)
			      puts("    Confidence values are all zero - using sequence");
			  use_conf = 0;
		      }
		  }
		  if (use_conf) {
		      
		      length_hidden = MIN
			  (length_hidden,scan_right
			   (p, conf, r.length-r.start+1, r.length)
			   -(r.length-r.start));
		  }
		  
		  if (!use_conf) {
		      
		      length_hidden = end_of_good ( t_hidden_seq, p.start,
						   p.rwin1, p.rcnt1);
		  }
		  
		  /* Longer extension ?	*/
		  
		  new_end = r.position - length_hidden;
		  l = MIN ( left_most, new_end );
		  
		  if ( l < left_most ) {
		      left_most = l;
		      length_extension = length_hidden - r.position + 1;
		      
		      start_hidden = r.position - 1;
		      best_length_hidden = length_hidden;
		      (void) memcpy ( hidden_seq, 
				     t_hidden_seq,
				     best_length_hidden);
		  }
	      }
	  }
	  read_number = r.right;

	  if (conf) {
	      xfree(conf);
	      conf = NULL;
	  }
      }

      if (conf) {
	  xfree(conf);
	  conf = NULL;
      }
      
      /* when we get here if length_extension != 0 we have the
       * longest segment of hidden data in hidden_seq
       * in the correct orientation for alignment
       * It starts at start_hidden (-1 for element) in the consensus
       * and is of length best_length_hidden.
       */

      if ( length_extension != 0 ) {
	  
	  if ( start_hidden > 0 ) {

	      complement_seq(consensus, start_hidden);

	      (void) memcpy (t_hidden_seq, 
			     hidden_seq,
			     best_length_hidden);
	      init_overlap (overlap, 
	      consensus,
	      t_hidden_seq,
	      start_hidden,
	      best_length_hidden);
	      /*
	      overlap.seq1 = consensus;
	      overlap.seq1_len = start_hidden;
	      overlap.seq2 = t_hidden_seq;
	      overlap.seq2_len = best_length_hidden;
	      */
	      params.band = MIN(band,overlap->seq1_len);
	      params.band_left = -params.band;
	      params.band_right = params.band;

	      ret = affine_align(overlap, &params);

	      if ( ret < 0 ) {
		  destroy_overlap(overlap);
		  return -2;
	      }
	      percent_mismatch = 100.0 - overlap->percent;

	      if (percent_mismatch < max_percent_mismatch ) {

		  ret = get_segment(overlap,seg,2);

		  if ( ret < 0 ) {
		      destroy_overlap(overlap);
		      return -2;
		  }
		  length_extension = seg->length;
	      }
	      else {
		  length_extension = 0;
		  destroy_overlap(overlap);
		  return -3;
	      }
	      complement_seq(consensus, start_hidden);
	  }
	  if(length_extension>0)complement_seq(hidden_seq, length_extension);
      }
      destroy_overlap(overlap);
      return length_extension;
      break;

  default:

      verror(ERR_WARN, "get_hidden", "unknown end");
      destroy_overlap(overlap);
      return -2;
      break;
  }
}
#endif

static int get_hidden_start(GapIO *io, tg_rec contig, Hidden_params p,
			    char *cons, char *hidden_seq) {
    int start;
    rangec_t *r;
    int nr, i, max_ext = 0;
    contig_t *c;

    if (consensus_valid_range(io, contig, &start, NULL) != 0) {
	return -1;
    }

    c = cache_search(io, GT_Contig, contig);
    cache_incr(io, c);

    if (!(r = contig_seqs_in_range(io, &c, start-1, start-1, 0, &nr))) {
	cache_decr(io, c);
	return -1;
    }

    //printf("Contig #%"PRIrec" start has %d seqs\n", contig, nr);

    /* All these reads extend the contig by at least one base */
    for (i = 0; i < nr; i++) {
	seq_t *s = cache_search(io, GT_Seq, r[i].rec);
	seq_t *sorig = s;
	int cstart, cend;
	int j, k, slen, lclip, ext;

	s = dup_seq(s);
	if ((s->len < 0) ^ r[i].comp) {
	    complement_seq_t(s);
	}

	/* Contig start/end coords */
	cstart = r[i].start + s->left-1;
	cend   = r[i].start + s->right-1;

	revconf(s->conf, ABS(s->len));
	lclip = scan_right(p, s->conf, ABS(s->len) - s->left + 1, ABS(s->len));
	lclip = ABS(s->len) - lclip+1;
	ext = start - (r[i].start + lclip-1);

	//printf("#%"PRIrec", %d %.*s\n", 
	//       s->rec,
	//       r[i].start + s->left-1 - start,
	//       ext, s->seq+lclip-1);

	if (ext > 0 && r[i].start + s->left-1 > start && 0) {
	    OVERLAP *overlap;
	    ALIGN_PARAMS *params;
	    int j, k, jx, kx, j_end;
	    char rseq[MAXGEL_PLUS], rcons[MAXGEL_PLUS];
	    int rseq_len, rcons_len;

	    //printf("Hidden_seq for #%"PRIrec" = %.*s\n",
	    //	   r[i].rec, ABS(s->len)-s->right, &s->seq[s->right]);

	    if (NULL == (params = create_align_params()) ||
		NULL == (overlap = create_overlap())) {
		if (params)
		    destroy_alignment_params(params);
		cache_decr(io, c);
		return -1;
	    }

	    set_align_params (params,
			      20, /* band */
			      12, /* gap_open */
			      4,  /* gap_extend */
			      EDGE_GAPS_COUNT | BEST_EDGE_TRACE,
			      /* RETURN_SEQ | */ RETURN_EDIT_BUFFERS,
			      0, 0, 0, 0,0);

	    /*
	     * I've no idea why, but the consensus algorithm has big issues
	     * anchoring itself on the right end, but anchored left with
	     * floating right works fine. So we reverse complement to get
	     * the sequences in a form that works. Yuk!
	     */
	    if (r[i].start + s->left-1 - start > MAXGEL_PLUS) {
		memcpy(rcons, cons, MAXGEL_PLUS);
		rcons_len = MAXGEL_PLUS;
	    } else {
		memcpy(rcons, cons, r[i].start + s->left-1 - start);
		rcons_len = r[i].start + s->left-1 - start;
	    }
	    complement_seq(rcons, rcons_len);

	    if (s->left-1 > MAXGEL_PLUS) {
		memcpy(rseq, s->seq + s->left-1 - MAXGEL_PLUS, MAXGEL_PLUS);
		rseq_len = MAXGEL_PLUS;
	    } else {
		memcpy(rseq, s->seq, s->left-1);
		rseq_len = s->left-1;
	    }
	    complement_seq(rseq, rseq_len);

	    init_overlap(overlap, rcons, rseq, rcons_len, rseq_len);

	    affine_align(overlap, params);

	    //puts(overlap->seq1_out);
	    //puts(overlap->seq2_out);

	    j = k = jx = kx = 0;
	    /* Trim pads on end of consensus */
	    if (overlap->S1[overlap->s1_len-1] < 0)
		overlap->s1_len--;

	    while (jx < overlap->s1_len && kx < overlap->s2_len) {
		char bc, br;

		/* Cons base */
		if (overlap->S1[jx] > 0) {
		    if (--overlap->S1[jx] == 0)
			jx++;
		    bc = overlap->seq1[j++];
		} else {
		    if (++overlap->S1[jx] == 0)
			jx++;
		    bc = '*';
		}

		/* Read base */
		if (overlap->S2[kx] > 0) {
		    if (--overlap->S2[kx] == 0)
			kx++;
		    br = overlap->seq2[k++];
		} else {
		    if (++overlap->S2[kx] == 0)
			kx++;
		    br = '*';
		}

		//printf("%3d %c/%c %3d\n", j, bc, br, k);
	    }

	    /* Complement back again */
	    k = rseq_len - k - 1;
	    ext = k - lclip + 1;

	    destroy_alignment_params(params);
	    destroy_overlap(overlap);
	}

	if (max_ext < ext) {
	    max_ext = ext;
	    memcpy(hidden_seq, &s->seq[lclip-1], ext);
	    hidden_seq[ext] = 0;
	    //printf("New hidden_seq=%s\n", hidden_seq);
	}

	if (sorig != s)
	    free(s);
    }

    cache_decr(io, c);
    return max_ext;
}

static int get_hidden_end(GapIO *io, tg_rec contig, Hidden_params p,
			  char *cons, char *hidden_seq) {
    int end;
    rangec_t *r;
    int nr, i, max_ext = 0;
    contig_t *c;

    if (consensus_valid_range(io, contig, NULL, &end) != 0) {
	return -1;
    }

    c = cache_search(io, GT_Contig, contig);
    cache_incr(io, c);

    if (!(r = contig_seqs_in_range(io, &c, end+1, end+1, 0, &nr))) {
	cache_decr(io, c);
	return -1;
    }

    //printf("Contig #%"PRIrec" end has %d seqs\n", contig, nr);

    /* All these reads extend the contig by at least one base */
    for (i = 0; i < nr; i++) {
	seq_t *s = cache_search(io, GT_Seq, r[i].rec);
	seq_t *sorig = s;
	int cstart, cend;
	int j, k, slen, rclip, rfrom, ext;

	if ((s->len < 0) ^ r[i].comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	/* Contig start/end coords */
	cstart = r[i].start + s->left-1;
	cend   = r[i].start + s->right-1;

	/*
	 * We may have a sequence that isn't at the end. Eg:
	 *
	 * Cons: AGACTTAGCGTAGACCC
	 * Read: AGAgttagccgtagacccgttagcga
	 *
	 * It's not in alignment perfectly, so we align to work
	 * out the precise piece to attach to the end; "gttagcga"
	 * in this example.
	 */
	rclip = scan_right(p, s->conf, s->left, ABS(s->len));
	ext = rclip-1 + r[i].start - end;

	if (ext > 0 && end > r[i].start + s->right-1) {
	    OVERLAP *overlap;
	    ALIGN_PARAMS *params;
	    int j, k, jx, kx;

	    //printf("Hidden_seq for #%"PRIrec" = %.*s\n",
	    //	   r[i].rec, ABS(s->len)-s->right, &s->seq[s->right]);

	    if (NULL == (params = create_align_params()) ||
		NULL == (overlap = create_overlap())) {
		if (params)
		    destroy_alignment_params(params);
		cache_decr(io, c);
		return -1;
	    }

	    set_align_params (params,
			      20, /* band */
			      12, /* gap_open */
			      4,  /* gap_extend */
			      /* left is fixed, right not */
			      /* vs left is floating, right fixed */
			      /* EDGE_GAPS_ZERO | FULL_LENGTH_TRACE */
			      EDGE_GAPS_COUNT | BEST_EDGE_TRACE,
			      /*RETURN_SEQ | */RETURN_EDIT_BUFFERS,
			      0, 0, 0, 0,0);

	    init_overlap(overlap,
			 &cons[r[i].start + s->right-1],
			 &s->seq[s->right],
			 end - (r[i].start + s->right-1),
			 rclip - s->right);

	    affine_align(overlap, params);

	    //puts(overlap->seq1_out);
	    //puts(overlap->seq2_out);

	    j = k = jx = kx = 0;
	    /* Trim pads on end of consensus */
	    if (overlap->S1[overlap->s1_len-1] < 0)
		overlap->s1_len--;

	    while (jx < overlap->s1_len && kx < overlap->s2_len) {
		char bc, br;

		/* Cons base */
		if (overlap->S1[jx] > 0) {
		    if (--overlap->S1[jx] == 0)
			jx++;
		    bc = overlap->seq1[j++];
		} else {
		    if (++overlap->S1[jx] == 0)
			jx++;
		    bc = '*';
		}

		/* Read base */
		if (overlap->S2[kx] > 0) {
		    if (--overlap->S2[kx] == 0)
			kx++;
		    br = overlap->seq2[k++];
		} else {
		    if (++overlap->S2[kx] == 0)
			kx++;
		    br = '*';
		}

		//printf("%3d %c/%c %3d\n", j, bc, br, k);
	    }

	    /*
	     * At this stage j should be the end of the consensus
	     * and k should be the offset into the sequence for the first
	     * base beyond the consensus.
	     */
	    rfrom = k + s->right;
	    //printf("rfrom=%d rclip=%d, sright=%d, old ext=%d, new=%d\n",
	    //	   rfrom, rclip, s->right, ext, rclip-rfrom);
	    ext = rclip - rfrom + 1;

	    destroy_alignment_params(params);
	    destroy_overlap(overlap);
	} else {
	    rfrom = s->right;
	}

	if (max_ext < ext) {
	    max_ext = ext;
	    memcpy(hidden_seq, &s->seq[rfrom], ext);
	    hidden_seq[ext] = 0;
	    //printf("New hidden_seq=%s\n", hidden_seq);
	}

	if (sorig != s)
	    free(s);
    }

    cache_decr(io, c);
    return max_ext;
}

/*
 * Looks at the clipped/hidden data off either the start (end == LEFT_END) or
 * end (end == RIGHT_END) of the contig and chooses one that extends the
 * contig the furthest via good quality data.
 *
 * During this process we may need to perform sequence alignments against
 * the existing contig if the best candidate for extension is not the same
 * as the current leftmost or rightmost sequence. The 'p' and 'consensus'
 * parameters govern this behaviour.
 *
 * Returns the number of bases extended, with the sequence in hidden_seq,
 *         or 0 for no extension
 *         or -1 for failure.
 */
static int get_hidden(GapIO *io, tg_rec contig, int end,
		      Hidden_params p, char *consensus,
		      char *hidden_seq) {
    if (end == LEFT_END)
	return get_hidden_start(io, contig, p, consensus, hidden_seq);
    else if (end == RIGHT_END)
	return get_hidden_end(io, contig, p, consensus, hidden_seq);
    else
	return -1;
}

/************************************************************/

int make_consensus( int task_mask, GapIO *io,
		   char **consensus2, float *quality,
		   Contig_parms *contig_list, int number_of_contigs,
		   int *consensus_length, int max_read_length,
		   Hidden_params p, float percd ) {
		   
    tg_rec contig, left_gel_number;
    int i,j, start, end;
    int contig_length, consensus_start, contig_start;
    int left_extension, right_extension, max_consensus;
    char *hidden_seq,*t_hidden_seq;
    /*char hidden_seq[MAXGEL_PLUS],t_hidden_seq[MAXGEL_PLUS];*/
    char *project_name, *consensus = NULL;

/* routine to handle all consensus calculations for the assembly program

   There are several different types of "consensus" required:

   1. all contigs or one contig
   2. contigs sorted on length
   3. contigs with added hidden data
   4. contigs masked
   5. single stranded
   6. normal consensus
   7. quality codes
   8. with title
   9. without an added title

   Do them all in one loop by putting the contig numbers in a list first.

   Tell consen what to do by using variable "task_mask" which has the
   appropriate bits from the list above set.

   If we pass in 'quality' as non NULL, it is filled in with the quality
   return values from calc_consensus. These are either log error rates or
   consensus cutoff values depending on the global consensus_mode variable.
   'quality' is assumed to be the same size as 'consensus' and all data
   copies are kept in sync. Eg ADD_TITLE to consensus will leave a title
   lengthed hole in the quality array.

   Return values: 0 success
                 -1 insufficient space
                 -2 masking failed
		 -3 add hidden failed
*/

    project_name = io_name(io);

    hidden_seq = t_hidden_seq = NULL;
    
    max_read_length = find_max_gel_len(io, 0, 0);
    if ( task_mask & ADDHIDDENDATA ) {
      /* need to allocate twice allowed length to satisfy the alignment routines*/
	if ((hidden_seq = (char *)xmalloc(2*MAXGEL_PLUS * sizeof(char)))==NULL){
	    return(-1);
	}
	if ((t_hidden_seq = (char *)xmalloc(2*MAXGEL_PLUS * sizeof(char)))==NULL){
	    if (hidden_seq) xfree(hidden_seq);
	    return(-1);

	}
    }
    /*	number_of_contigs is the number of consensus sequences to produce */
    /*	we have a list of contig numbers in contig_list	*/


    /* Is there enough space in consensus? We cannot check for ADDHIDDENDATA!*/

    for ( i=0,j=0;i<number_of_contigs;i++ ) {
	j += contig_list[i].contig_end - contig_list[i].contig_start + 1;
    }
    if ( task_mask & ADDTITLE ) {
	j += number_of_contigs * 20;
    }

    if ( task_mask & ADDHIDDENDATA ) {
	/* Worst case, but may need to shrink & fail */
	j += number_of_contigs * MAXGEL_PLUS*2;
    }
    max_consensus = j;
    if (NULL == (consensus = (char *)malloc(j+1)))
	return -1;

    if ( task_mask & SORTCONTIGS ) {
	if ( i = sort_contigs ( contig_list, number_of_contigs ) ) {
	    if (consensus) xfree(consensus);
	    return -2;
	}
    }

    consensus_start = *consensus_length; 

    /*	loop for all contigs	*/

    for ( i=0;i<number_of_contigs;i++ ) {

	contig = contig_list[i].contig_number;
	left_extension = 0;
	right_extension = 0;
	start  = contig_list[i].contig_start;
	end    = contig_list[i].contig_end;
	contig_length = end - start + 1;

	if ( task_mask & ADDTITLE ) {

	    /* is there enough space left ? */
	    if ( (*consensus_length + 20) > max_consensus ) {
		/* This shouldn't occur now due to precomputing this */
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		if (consensus) xfree(consensus);
		return -1;
	    }
	    left_gel_number = contig_list[i].contig_left_gel;

	    (void) add_contig_title ( &consensus[*consensus_length], 
				     project_name, left_gel_number);
	    *consensus_length += 20;
	}
	contig_start = *consensus_length;
	/* note the element number of the first base each contig */

	contig_list[i].contig_start_offset = *consensus_length - consensus_start;

	if ( task_mask & NORMALCONSENSUS ) {
	    /* is there enough space left ? */

	    if ( (*consensus_length + contig_length) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		if (consensus) xfree(consensus);
		return -1;
	    }
	    calc_consensus(contig, start, end, CON_SUM,
			   &consensus[*consensus_length], NULL,
			   quality ? &quality[*consensus_length] : NULL, NULL,
			   percd, quality_cutoff, database_info, (void *)io);
	    *consensus_length += contig_length;
	}

	else if ( task_mask & SINGLESTRANDED ) {
	    /* is there enough space left ? */

	    if ( (*consensus_length + contig_length) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		if (consensus) xfree(consensus);
		return -1;
	    }
	    calc_consensus(contig, start, end, CON_WDET,
			   &consensus[*consensus_length], NULL,
			   quality ? &quality[*consensus_length] : NULL, NULL,
			   percd, quality_cutoff, database_info, (void *)io);
	    *consensus_length += contig_length;
	}

	else if ( task_mask & QUALITYCODES ) {
	    /* is there enough space left ? */

	    if ( (*consensus_length + contig_length) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		if (consensus) xfree(consensus);
		return -1;
	    }
	    /*
	     * R_GOOD_GOOD_EQ	'a'
	     * R_GOOD_BAD	'b'
	     * R_BAD_GOOD	'c'
	     * R_GOOD_NONE	'd'
	     * R_NONE_GOOD	'e'
	     * R_BAD_BAD 	'f'
	     * R_BAD_NONE	'g'
	     * R_NONE_BAD	'h'
	     * R_GOOD_GOOD_NE	'i'
	     * R_NONE_NONE	'j'
	     */

	    calc_quality(contig, start, end, 
			   &consensus[*consensus_length], 
			   percd, quality_cutoff, database_info, (void *)io);

	    *consensus_length += contig_length;
	}

	if ( task_mask & ADDHIDDENDATA ) {
	    assert(hidden_seq);

	    left_extension = get_hidden(io, contig, LEFT_END,
					p, &consensus[contig_start],
					hidden_seq);
	    //printf("L %d %.*s\n", left_extension, left_extension, hidden_seq);
	    /* check for normal operation */
	    if ( left_extension < 0 ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		if (consensus) xfree(consensus);
		return -3;
	    }
	    /* is there enough space left ? */
	    if ( (*consensus_length + left_extension) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		if (consensus) xfree(consensus);
		return -1;
	    }
	    (void) memmove ( &consensus[contig_start+left_extension],
			   &consensus[contig_start], contig_length);

	    (void) memcpy ( &consensus[contig_start], hidden_seq,
			   left_extension);
	    if (quality) {
		int i;
		for (i = 0; i < left_extension; i++)
		    quality[contig_start + i] = 0.0;
	    }
	    *consensus_length += left_extension;
	    contig_list[i].contig_left_extension = left_extension;

	    right_extension = get_hidden(io, contig, RIGHT_END,
					 p, &consensus[contig_start +
						       left_extension],
					 hidden_seq);
	    //printf("R %d %.*s\n", right_extension, right_extension, hidden_seq);
	    /* check for normal operation */

	    if ( right_extension < 0 ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		if (consensus) xfree(consensus);
		return -3;
	    }
	    /* is there enough space left ? */

	    if ( (*consensus_length + right_extension) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		if (consensus) xfree(consensus);
		return -1;
	    }
	    (void) memcpy ( &consensus[*consensus_length], hidden_seq,
				       right_extension);
	    if (quality) {
		int i;
		for (i = 0; i < right_extension; i++)
		    quality[*consensus_length + i] = 0.0;
	    }
	    *consensus_length += right_extension;
	    contig_list[i].contig_right_extension = right_extension;

	}
#if 0
	if ( task_mask & MASKING ) {
/*	    printf("do masking\n");*/
            if ( mask_consensus(io, 
			       &consensus[*consensus_length - contig_length - 
					  right_extension], 
			       contig, start, end, 1) ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		if (consensus) xfree(consensus);
		return -2;
	    }
	}

	if ( task_mask & MARKING ) {
/*	    printf("do masking\n");*/
            if ( mask_consensus(io, 
			       &consensus[*consensus_length - contig_length - 
					  right_extension], 
			       contig, start, end, 4) ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		if (consensus) xfree(consensus);
		return -2;
	    }
	}
#endif
	/* note the element number of the last base each contig */

	contig_list[i].contig_end_offset = *consensus_length - consensus_start - 1;

    }
    if (hidden_seq) xfree(hidden_seq);
    if (t_hidden_seq) xfree(t_hidden_seq);
    *consensus2 = consensus;
    return 0;
}


/************************************************************/
/* 
 * read through the consensus and find the ends of the contigs.
 * store their positions in contig_ends[] and 
 * their left gel numbers in contig_numbers[]
 *
 * NOTE: contig_ends[] and contig_numbers[] should be allocated to
 * the number of contigs PLUS ONE.
 */
int find_contig_ends ( char *seq, int seq_len, 
		       int *contig_ends, tg_rec *contig_numbers ) {


    tg_rec left_gel;
    int i, contig_index;
    char *dot;

    contig_index = 0;
    for ( i = 0; i < seq_len; ) {

	if ( seq[i] == '<' ) {

	    /* found the next contig */
	    if ( !(dot = strchr ( &seq[i],'.' ))) {

		/* whoops, lets bail out */
		return 0;
	    }
	    left_gel = atoi ( ++dot );
	    contig_ends [ contig_index ] = i;
	    contig_numbers [ contig_index ] = left_gel;
/*	    printf("left gel %d left end %d\n",
		   contig_numbers [ contig_index ],
		   contig_ends [ contig_index ]);
*/
	    contig_index++;
	    i += 20;
	}
	else i++;
    }

    /* stick position of last char+1 in next element */
    contig_ends [ contig_index ] = seq_len;
    return contig_index;
}


/* writes an array of chars, 60 per line */
int plain_fmt_output( FILE *fp, char *seq, int seq_len, int nopads)
{
#define LINELENGTH 60
    int i, j;

    j = i = 0;
    while (i < seq_len) {
	for (j = 0; j < LINELENGTH;) {
	    if (nopads == 0 || seq[i] != '*') {
		if (fprintf ( fp, "%c", seq[i]) < 0) return 1;
		j++;
	    }
	    if (++i >= seq_len) {
		break;
	    }
	}
	if (fprintf ( fp, "\n") < 0) return 1;
    }
    return 0;
}


char *set_fasta_table(void) {

/* set up table of values for converting to fasta characters */
/* This is not the strict fasta alphabet, but this is better for real work. */


    int i;
    char *fasta_base_table;
    char *codes = "ACGTRYMWSKDHVBDEFI";
    int codes_l = strlen(codes);

    if ( NULL == (fasta_base_table = (char *) malloc ( 257 * sizeof ( char ) )))
	return NULL;

    for (i=0;i<256;i++) fasta_base_table[i] = 'n';

    for (i = 0; i < codes_l; i++) {
	int lower = tolower(codes[i]);
	fasta_base_table[codes[i]] = lower;
	fasta_base_table[lower] = lower;
    }

    fasta_base_table['*'] = '*';

    return fasta_base_table;
}

/************************************************************/

int convert_to_fasta ( char *seq, int *seq_len_p, int nopads ) {

    int i, j;
    char *fasta_base_table;
    int seq_len = *seq_len_p;

    if ( NULL == (fasta_base_table = set_fasta_table()) ) return 1;

    if (nopads) {
	for ( i = j = 0; i < seq_len; i++){
	    if (seq[i] != '*')
		seq[j++] = fasta_base_table [ (unsigned int)seq[i] ];
	}
	seq[j] = 0;
	*seq_len_p = j;
    } else {
	for ( i = 0; i < seq_len; i++){
	    seq[i] = fasta_base_table [ (unsigned int)seq[i] ];
	}
    }
    free ( fasta_base_table );
    return 0;
}


int fasta_fmt_output ( FILE *fp, char *seq, int seq_len, char *entry_name,
		      int nopads, char *title) {

    /* do fasta format output */

    if ( convert_to_fasta ( seq, &seq_len, nopads )) return 1;

/*  entry_name is project name for single contig but left read otherwise
    title is always left read number
*/

    fprintf(fp,">%s %s \n",entry_name, title);

    /* do the sequence */
    if ( plain_fmt_output( fp, seq, seq_len, nopads ) != 0 ) return 1;
    return 0;
}

#if 0
/*
 * Creates a consensus experiment file containing tags
 * 'gel_anno' specifies whether to output annotations from gel readings.
 * 'truncate' specifies whether to allow tags within the cutoff data to be
 * output. 
 *
 * Returns 0 for success, -1 for failure.
 */

int expt_fmt_output(GapIO *io, mFILE *fp, char *seq, float *qual,
		    tg_rec left_read, int lreg, int rreg,
		    tg_rec gel_anno, int truncate, int gel_notes, int nopads) {
    GContigs c;
    GReadings r;
    Exp_info *e;
    int err = 0, gel, contig;
    char *name;
    char *new_seq = seq;
    int new_rreg = rreg;
    int *pads = NULL;
    int1 *new_qual = NULL;

    if (nopads) {
	int i, j;

	/* Strip pads from seq and quality */
	if (NULL == (new_seq = xmalloc(rreg - lreg + 2)))
	    return -1;

	if (NULL == (pads = (int *)xcalloc((rreg - lreg + 2), sizeof(int)))) {
	    xfree(new_seq);
	    return -1;
	}
	
	if (qual) {
	    if (NULL == (new_qual = (int1 *)xmalloc(4 * (rreg -lreg + 2) *
						    sizeof(int1)))) {
		xfree(new_seq);
		xfree(pads);
		return -1;
	    }

	    for (i = j = 0; i <= rreg - lreg; i++) {
		if (seq[i] != '*') {
		    new_qual[j++] = (int1)(qual[i]+0.499);
		}
	    }
	}

	for (i = j = 0; i <= rreg - lreg; i++) {
	    if (seq[i] != '*') {
		new_seq[j++] = seq[i];
	    }
	    pads[i] = i+1-j;
	}
	new_seq[j] = 0;
	new_rreg = j-1 + lreg;
    } else {
	/* Convert quality from float to int1 */
	if (qual) {
	    int i, j;

	    if (NULL == (new_qual = (int1 *)xmalloc(4 * (rreg -lreg + 2) *
						    sizeof(int1)))) {
		return -1;
	    }

	    for (i = j = 0; i <= rreg - lreg; i++, j++) {
		new_qual[j] = (int1)(qual[i]+0.499);
	    }

	    new_rreg = rreg;
	}
    }

    

    /* allocate experiment file info */
    e = exp_create_info();
    e->fp = fp;

    if ( (contig = rnumtocnum ( io, left_read )) == -1 ) return -1;
    contig_read(io, contig, c);
    name = io_rname(io, c.left);

    err |= exp_put_str(e, EFLT_ID, name, strlen(name));
    err |= exp_put_str(e, EFLT_EN, name, strlen(name));

    /* Contig annotations */
    err |= output_annotations(io, e, c.annotations, 1-lreg, 1, 0, 0, 1,
			      lreg-1, rreg+1, "Contig Annotations",
			      pads, rreg - lreg + 1);

    /* Contig notes */
    err |= output_notes(io, e, c.notes, GT_Contigs, contig);

    if (gel_notes) {
	/* Reading notes */
	for (gel = c.left; gel; gel = r.right) {
	    
	    gel_read(io, gel, r);
	    if (r.notes) {
		err |= output_notes(io, e, r.notes, GT_Readings, gel);
	    }
	}
    }

    if (gel_anno) {
	/* Reading annotations  */
	for (gel = c.left; gel; gel = r.right) {
	    char buf[100], *name;
	    
	    gel_read(io, gel, r);

	    if (r.position + r.sequence_length < lreg)
		continue;

	    if (r.position > rreg)
		break;

	    if (r.annotations) {
		name = io_rname(io, gel);
		sprintf(buf, "Annotations for reading %s", name);

		err |= output_annotations(io, e, r.annotations,
					  r.position - r.start - lreg,
					  0, r.sense, r.length, 0,
					  truncate ? r.start : 0,
					  truncate ? r.end : 0, buf,
					  pads, rreg - lreg + 1);
	    } else {
		/*
		 * Output blank CC line anyway so that the consenus file can
		 * be parsed to obtain a list of all readings, not just those
		 * containing tags. Requested by Gos Micklem.
		 */
		name = io_rname(io, gel);
		sprintf(buf, "Annotations for reading %s", name);
		exp_put_str(e, EFLT_CC, buf, strlen(buf));
	    }
	}
    }

    if (qual) {
	char *buf;

	buf = (char *)xmalloc(5 * (new_rreg-lreg+1));

	if (buf && new_qual) {
	    conf2str(new_qual, new_rreg-lreg+1, buf);
	    err |= exp_put_str(e, EFLT_AV, buf, strlen(buf));
	    xfree(buf);
	}
    }
    err |= exp_put_str(e, EFLT_SQ, new_seq, rreg-lreg+1);

    /* deallocate experiment file info - will close fp for us unless
     we set it to NULL */
    e->fp = NULL;
    exp_destroy_info(e);

    if (err == -1) {
	verror(ERR_WARN,"extract_consensus","Writing experiment file failed.");
    }

    if (nopads) {
	xfree(new_seq);
	xfree(pads);
    }
    if (new_qual)
	xfree(new_qual);

    return err ? -1 : 0;
}
#endif


#define STADEN_FORMAT 1
#define FASTA_FORMAT 2
#define EXPT_FORMAT 3
#define INTERNAL_FORMAT 4
#if 0    
int write_consensus (GapIO *io, FILE *fp, 
		     int output_format, char *seq, float *qual, int seq_len,
		     int max_contigs,
		     int gel_anno, int truncate, int gel_notes,
		     int num_contigs, contig_list_t *contig_array,
		     int nopads, int name_format) {

    int contig_index, number_of_contigs, *contig_ends;
    tg_rec *contig_numbers;
    mFILE *mf = NULL;

    contig_ends    = (int *)    malloc ( max_contigs * sizeof ( int ) );
    contig_numbers = (tg_rec *) malloc ( max_contigs * sizeof ( tg_rec ) );

    number_of_contigs =
	find_contig_ends (seq, seq_len, contig_ends, contig_numbers);

    if (number_of_contigs != num_contigs) {
	verror(ERR_WARN, "write_consensus", "number of contigs is unknown!");
	return 1;
    }

    /* Now loop for each contig writing out the data in the selected format */
    for ( contig_index = 0; contig_index < num_contigs; contig_index++) {

	if ( output_format == STADEN_FORMAT ) {

	    /* do the horrible header */
	    {
		/* Work around DU V4.0 printf bug - %.20s calls strlen() */
		char buf[21];
		memcpy(buf, &seq[contig_ends[contig_index]], 20);
		buf[20] = 0;

		if (fprintf(fp,"%.20s\n", buf) < 0)
		    goto error;
	    }

	    /* do the sequence */
	    if (plain_fmt_output(fp, &seq[contig_ends[contig_index]+20], 
				 contig_ends[contig_index+1] - 
				 (contig_ends[contig_index]+20), nopads) != 0)
		goto error;

	} else if ( output_format == FASTA_FORMAT ) {
	    char title[1024], *entry_name_ptr;
	    char tname[DB_NAMELEN+1];

	    sprintf(title, "%s.%d",
		    io_name(io), contig_numbers[contig_index]);
	    if (name_format == CONS_NAME_LREADING) {
		entry_name_ptr = io_rname(io, contig_numbers[contig_index]);
	    } else /* if (name_format == CONS_NAME_LTEMPLATE) */ { 
		GReadings r;
		GTemplates t;

		gel_read(io, contig_numbers[contig_index], r);
		if (r.template) {
		    template_read(io, r.template, t);
		    TextRead(io, t.name, tname, DB_NAMELEN);
		    tname[DB_NAMELEN] = 0;
		    entry_name_ptr = tname;
		} else {
		    entry_name_ptr = "?";
		}
	    }
	    if ( fasta_fmt_output (fp, 
				   &seq[contig_ends[contig_index]+20], 
				   contig_ends[contig_index+1] - 
				   (contig_ends[contig_index]+20),
				   entry_name_ptr, nopads, title))
		goto error;

	} else if ( output_format == EXPT_FORMAT ) {
	    if (!mf)
		mf = mfreopen(NULL, "w", fp);
	    if( expt_fmt_output( io, mf,
				&seq[contig_ends[contig_index]+20],
				qual ? &qual[contig_ends[contig_index]+20]
				     : NULL,
				contig_numbers[contig_index],
				contig_array[contig_index].start,
				contig_array[contig_index].end,
				gel_anno, truncate, gel_notes, nopads ) )
		goto error;
	}
    }

    if (mf) {
	mfflush(mf);
	mfdestroy(mf);
    }
    
    free ( contig_ends );
    free ( contig_numbers );
    return 0;
 error:
    free ( contig_ends );
    free ( contig_numbers );
    return 1;
}
#endif
