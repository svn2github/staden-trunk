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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <fcntl.h>
#include <unistd.h>
#include "io_utils.h"
#include "IO.h"
#include "fort.h"
#include "FtoC.h"
#include "expFileIO.h"
#include "misc.h"
#include "edUtils.h"
#include "tagUtils.h"
#include "consen.h"
#include "gap_globals.h"
#include "dna_utils.h"
#include "list_proc.h"
#include "extract.h"
#include "align_lib.h"

/*#define MAXGEL 1024*/
#define MAXGEL_PLUS 1024

/************************************************************/

/*
 * get a list of unattatched readings ie the names of the
 * readings that are in contigs with only one reading
 */
char *unattached_reads (GapIO *io) {

    int i;
    GContigs c;
    char *name;
    int number_of_contigs = NumContigs(io);
    void *dl = alloc_dlist();
    char *res;
    
    for ( i=0;i<number_of_contigs;i++ ) {

	GT_Read(io, arr(GCardinal, io->contigs, i),
		&c, sizeof(c), GT_Contigs);
	if ( c.left == c.right ) {
	    name = io_rname(io, c.left);
	    add_to_dlist(dl, name);
	}
    }

    res = strdup(read_dlist(dl));
    free_dlist(dl);

    return res;
}

/* Search for a minimal set of reads covering a contig */

int rr_read (GapIO *io, 
	     int first_read,
	     int max_length)

{

/* routine given first_read, finds the read that overlaps first_read,
   within a distance max_length, and extends furthest right. Returns
   the read number or 0 for end reached. If we do not find one within
   max_length we return the first one found.

   Keep reading until we go passed the end of the first read or further
   away than max_length, or reach the end of the contig. For each read
   test the position of its right end. If it is further out than the
   current furthest, then note its number (rightmost_read) and the position
   it reaches (furthest_right). We do not worry about strands and are
   assuming that the best reads to choose are those that gave the longest
   visible data first time around */

    GReadings r;
    int furthest_right, rightmost_read, current_read, not_at_end;
    int current_range, end_position;

    rightmost_read = 0;
    current_read = first_read;
    not_at_end = 1;

    /* get position for end of search (end_position) */

    gel_read(io, current_read, r);
    end_position = r.position + MIN(max_length, r.sequence_length) - 1;
    furthest_right = end_position;
    current_read = r.right;
    if ( r.right == 0 ) not_at_end = 0;

    while ( not_at_end ) {

	/* get next reading */

	gel_read(io, current_read, r);
	    
	/* set flag if passed end of allowed range */

	if ( r.position >= end_position )
	    not_at_end = 0;

	if (( not_at_end ) || ( rightmost_read == 0 )) {
	    current_range = r.position +
		MIN(r.sequence_length, max_length) - 1;
	    if ( current_range > furthest_right ) {
		furthest_right = current_range;
		rightmost_read = current_read;
	    }
	    if ( r.right != 0 ) {

		current_read = r.right;
	    }
	    else 
		/* set flag for end of contig */

		not_at_end = 0;

	}

    }

    return rightmost_read;
}

    

char *minimal_coverage (GapIO *io, int num_c, contig_list_t *contigs) {
    char *name;
    int i, max_length, current_read;
    GContigs c;
    void *dl;
    char *out;

    max_length = 30000; /* FIXME: work out what to do with this in future.
			 * Currently it's set high to disable it.
			 */

    dl = alloc_dlist();

    for ( i=0;i<num_c;i++ ) {

	GT_Read(io, arr(GCardinal, io->contigs, contigs[i].contig-1),
		&c, sizeof(c), GT_Contigs);
	current_read = c.left;

	name = io_rname(io, current_read);
	add_to_dlist(dl, name);

	while ( current_read = rr_read( io, current_read, max_length) ) {
	    name = io_rname(io, current_read);
	    add_to_dlist(dl, name);
	}
    }

    out = strdup(read_dlist(dl));
    free_dlist(dl);

    return out;
}



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
void add_contig_title(char *consensus, char *project_name, int left_gelnumber)
{

    int  plen, rlen;
    static char dashes[]="----------------";
    char buf[50];
    char *cp;

    /* Find out how many digits for reading name */
    rlen = sprintf(buf, "%d", left_gelnumber);

    /* Length of project name */
    if (cp = strchr(project_name, '.'))
	plen = cp - project_name;
    else
	plen = strlen(project_name);
    if (plen + rlen + 3 /* < . > */ > 20) {
	plen = 20 - rlen - 3;
    }
    
    sprintf(consensus, "<%.*s.%.*d%.*s>",
	    plen, project_name,
	    rlen, left_gelnumber,
	    17-plen-rlen, dashes);
}

f_proc_ret cadtit_(char *consensus, char *project_name, f_int *leftgel,
		   f_implicit consensus_l, f_implicit name_l) {
    /*
     * Don't need to Fstr2Cstr as we never check the end of project_name,
     * we only look for the dot.
     */
    add_contig_title(consensus, project_name, *leftgel);
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
	    contig_list[i].contig_number = i+1;
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

    int i;

    if ( 0 == number_of_contigs ) return -1;
    if ( 1 == number_of_contigs ) return 0;

    for ( i=1; i<number_of_contigs; i++ ) {

	if ( pos_in_contig <= contig_list[i].contig_start_offset ) {
	    return i - 1;
	}
    }
    /* must be in last contig */
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
/**********************************************************************/

void maskc_ (char *seq, f_int *seq_len, f_int *jobin, f_implicit seq_l) {

    int seq_length, job;


    seq_length = *seq_len;
    job = *jobin;

    (void) maskit ( seq, seq_length, job );
}




/****************************************************************************/

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

    if (-1 == GT_Read(io, arr(GCardinal, io->contigs, contig-1),
		      &c, sizeof(c), GT_Contigs))
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

    rclip = i == len ? len + 1 : i;
    if (p.verbose)
	printf("    right clip = %d\n", rclip);

    return rclip;
}


int get_hidden_seq (GapIO *io, int read_number,
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

    if (-1 == GT_Read(io, arr(GCardinal, io->contigs, contig-1),
		      &c, sizeof(c), GT_Contigs)) {
	destroy_overlap(overlap);
	return -1;
    }

    length_extension = 0;

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

/************************************************************/

int make_consensus( int task_mask, GapIO *io,
		   char *consensus, float *quality,
		   Contig_parms *contig_list, int number_of_contigs,
		   int *consensus_length, int max_read_length,
		   int max_consensus, Hidden_params p, float percd ) {
		   
    int contig,left_gel_number,i,j, start, end;
    int contig_length, consensus_start, contig_start;
    int left_extension, right_extension;
    char *hidden_seq,*t_hidden_seq;
    /*char hidden_seq[MAXGEL_PLUS],t_hidden_seq[MAXGEL_PLUS];*/
    char *project_name;

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

    if ( j > max_consensus ) {
	if (hidden_seq) xfree(hidden_seq);
	if (t_hidden_seq) xfree(t_hidden_seq);
	return -1;
    }

    if ( task_mask & SORTCONTIGS ) {
	if ( i = sort_contigs ( contig_list, number_of_contigs ) ) return -2;
    }

    consensus_start = *consensus_length; 

    /*	loop for all contigs	*/

    for ( i=0;i<number_of_contigs;i++ ) {

	contig = contig_list[i].contig_number;
	left_extension = 0;
	right_extension = 0;

	if ( task_mask & ADDTITLE ) {

	    /* is there enough space left ? */
	    if ( (*consensus_length + 20) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
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
	    start  = contig_list[i].contig_start;
	    end    = contig_list[i].contig_end;
	    contig_length = end - start + 1;

	    /* is there enough space left ? */

	    if ( (*consensus_length + contig_length) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		return -1;
	    }
	    calc_consensus(contig, start, end, CON_SUM,
			   &consensus[*consensus_length], NULL,
			   quality ? &quality[*consensus_length] : NULL, NULL,
			   percd, quality_cutoff, database_info, (void *)io);
	    *consensus_length += contig_length;
	}

	else if ( task_mask & SINGLESTRANDED ) {
	    start  = contig_list[i].contig_start;
	    end    = contig_list[i].contig_end;
	    contig_length = end - start + 1;

	    /* is there enough space left ? */

	    if ( (*consensus_length + contig_length) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		return -1;
	    }
	    calc_consensus(contig, start, end, CON_WDET,
			   &consensus[*consensus_length], NULL,
			   quality ? &quality[*consensus_length] : NULL, NULL,
			   percd, quality_cutoff, database_info, (void *)io);
	    *consensus_length += contig_length;
	}

	else if ( task_mask & QUALITYCODES ) {

	    start  = contig_list[i].contig_start;
	    end    = contig_list[i].contig_end;
	    contig_length = end - start + 1;

	    /* is there enough space left ? */

	    if ( (*consensus_length + contig_length) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
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


	    left_extension = get_hidden(io, contig, max_read_length,LEFT_END,
					p, hidden_seq,
					t_hidden_seq, &consensus[contig_start]);
	    /* check for normal operation */
	    if ( left_extension < 0 ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		return -3;
	    }
	    /* is there enough space left ? */
	    if ( (*consensus_length + left_extension) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
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

	    right_extension = get_hidden(io, contig, max_read_length,RIGHT_END,
					p, hidden_seq, t_hidden_seq, 
					 &consensus[contig_start+left_extension]);
	    /* check for normal operation */

	    if ( right_extension < 0 ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
		return -3;
	    }
	    /* is there enough space left ? */

	    if ( (*consensus_length + right_extension) > max_consensus ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
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

	if ( task_mask & MASKING ) {
/*	    printf("do masking\n");*/
            if ( mask_consensus(io, 
			       &consensus[*consensus_length - contig_length - 
					  right_extension], 
			       contig, start, end, 1) ) {
		if (hidden_seq) xfree(hidden_seq);
		if (t_hidden_seq) xfree(t_hidden_seq);
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
		return -2;
	    }
	}

	/* note the element number of the last base each contig */

	contig_list[i].contig_end_offset = *consensus_length - consensus_start - 1;

    }
    if (hidden_seq) xfree(hidden_seq);
    if (t_hidden_seq) xfree(t_hidden_seq);
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
		       int *contig_ends, int *contig_numbers ) {


    int i,left_gel,contig_index;
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
    if ( NULL == (fasta_base_table = (char *) malloc ( 257 * sizeof ( char ) )))
	return NULL;

    for (i=0;i<256;i++) fasta_base_table[i] = 'n';

    fasta_base_table['a'] = 'a';
    fasta_base_table['c'] = 'c';
    fasta_base_table['g'] = 'g';
    fasta_base_table['t'] = 't';
    fasta_base_table['A'] = 'a';
    fasta_base_table['C'] = 'c';
    fasta_base_table['G'] = 'g';
    fasta_base_table['T'] = 't';

    fasta_base_table['d'] = 'd';
    fasta_base_table['e'] = 'e';
    fasta_base_table['f'] = 'f';
    fasta_base_table['i'] = 'i';
    fasta_base_table['D'] = 'd';
    fasta_base_table['E'] = 'e';
    fasta_base_table['F'] = 'f';
    fasta_base_table['I'] = 'i';

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

/*
 * Creates a consensus experiment file containing tags
 * 'gel_anno' specifies whether to output annotations from gel readings.
 * 'truncate' specifies whether to allow tags within the cutoff data to be
 * output. 
 *
 * Returns 0 for success, -1 for failure.
 */

int expt_fmt_output(GapIO *io, FILE *fp, char *seq, float *qual,
		    int left_read, int lreg, int rreg,
		    int gel_anno, int truncate, int gel_notes, int nopads) {
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
    GT_Read(io, arr(GCardinal, io->contigs, contig-1),
	    &c, sizeof(c), GT_Contigs);

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


    
int write_consensus (GapIO *io, FILE *fp, 
		     int output_format, char *seq, float *qual, int seq_len,
		     int max_contigs,
		     int gel_anno, int truncate, int gel_notes,
		     int num_contigs, contig_list_t *contig_array,
		     int nopads, int name_format) {

#define STADEN_FORMAT 1
#define FASTA_FORMAT 2
#define EXPT_FORMAT 3
#define INTERNAL_FORMAT 4

    int contig_index, number_of_contigs, *contig_ends, *contig_numbers;

    contig_ends    = (int *) malloc ( max_contigs * sizeof ( int ) );
    contig_numbers = (int *) malloc ( max_contigs * sizeof ( int ) );

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
	    
	    if( expt_fmt_output( io, fp,
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
    free ( contig_ends );
    free ( contig_numbers );
    return 0;
 error:
    free ( contig_ends );
    free ( contig_numbers );
    return 1;
}
/*
 routine to calculate and write consensus sequence(s) to a file
   in one of three formats: staden, fasta or expt. Several types
   of consensus are possible:

   job 1: normal
       2: extended - ie with hidden data added to ends
       3: unfinished - ie only regions that are not covered
                       by good data on both strands are written as a,c,g,t
		       the rest is "-".
       4: quality codes - ie a consensus of quality codes. Previously these 
                          have been 0,1,2,3,4,5,5,6,7, etc but for 
			  compatibility with most? sequence reading routines 
			  it is better to switch to a new character set of 
			  a,b,c,d,e,f, etc. Note for fasta output we are
			  currently changing characters to a,c,g,t,n!

   Algorithm: get consensus for all contigs in internal format 
   - ie <---lambda.0001--->acagatatattat< etc
   Then process each contig and call the appropriate output routines. Note
   that for expt output we may need to add in the tags etc.

*/

int make_consensus_files ( int task_mask, int output_format, int gel_anno,
			  int truncate, int gel_notes, FILE *fp, 
			  GapIO *io, char *consensus, float *quality,
			  int database_size, int nconts,
			  int *consensus_len, int max_read_length, 
			  int max_consensus, Hidden_params p, float percd,
			  int num_contigs, contig_list_t *contig_array,
			  int nopads, int name_format)

{

    int consensus_length;
    int max_contigs;
    int success;
    Contig_parms *contig_list;

    /* fiddle the index to the consensus array: 
       if consensus_length != 0 set it one back, otherwise leave it */

    consensus_length = 0;
    contig_list =  get_contig_list (database_size, io, num_contigs,
				    contig_array);

    if ( success =  make_consensus(task_mask, io, 
				   consensus, quality, contig_list,
				   num_contigs, 
				   &consensus_length,
				   max_read_length, max_consensus, p, percd)) {

	free ( contig_list );
	return success;
    }
    max_contigs = num_contigs + 1;

    success = write_consensus (io, fp, output_format, 
			       consensus, quality, consensus_length,
			       max_contigs, gel_anno, truncate, gel_notes,
			       num_contigs, contig_array, nopads,
			       name_format);
    free ( contig_list );
    *consensus_len = consensus_length;
    return success;
}



/************************************************************/

void precon_ (char *consensus, char *project_name, 
	      float *percd_in, int *idbsiz,
	      int *num_contigs, contig_list_t *contig_array,
	      int *task, int *idevr,
	      int *consensus_len, int *maxgel, int *maxcon,
	      int *window, int *nbad, int *iladd, int *iradd, 
	      int *iok)

{

    int i, consensus_length;
    int number_of_contigs, task_mask, max_read_length;
    int max_consensus, success = 1, database_size;
    GapIO *io;
    Contig_parms *contig_list;
    Hidden_params p;

    max_read_length = *maxgel;
    task_mask = *task;
    database_size = *idbsiz;

    consensus_length = *consensus_len;

    /* fiddle the index to the consensus array: 
       if consensus_length != 0 set it one back, otherwise leave it */

    consensus_length = MAX ( 0, consensus_length - 1 );

    max_consensus = *maxcon;

    /*get a pointer to the GapIO structure for the lower routines */
    if ( ( io = io_handle ( idevr ) ) == NULL ) {
	*iok = success;
	return;
    }

    p.min = p.max = p.verbose = p.use_conf = p.qual_val = p.window_len =0;
    p.test_mode = 0;
    p.start = 0;
    p.lwin1 = 0;
    p.lcnt1 = 0;
    p.rwin1 = *window;
    p.rcnt1 = *nbad;

    number_of_contigs = *num_contigs;
    contig_list =  get_contig_list(database_size, io,
				   number_of_contigs, contig_array);


    success =  make_consensus ( task_mask, io, consensus, NULL,
			       contig_list, number_of_contigs, 
			       &consensus_length, max_read_length, 
			       max_consensus, p, *percd_in);


    /* if required copy over the extension lengths to iladd and iradd */

    if ( task_mask & ADDHIDDENDATA ) {

	for ( i = 0;i<number_of_contigs;i++ ) {
	    iladd[i] = contig_list[i].contig_left_extension;
	    iradd[i] = contig_list[i].contig_right_extension;
	}
    }
    free ( contig_list );
    *consensus_len = consensus_length;
    *iok = success;
}

/*
 * As precon - but for a single contig only
 */
void precn1_ (char *consensus, char *project_name, 
	      float *percd_in, int *idbsiz,
	      int *contig, int *lreg, int *rreg,
	      int *task, int *idevr,
	      int *consensus_len, int *maxgel, int *maxcon,
	      int *window, int *nbad, int *iladd, int *iradd, 
	      int *iok)

{
    int one = 1;
    contig_list_t cl;

    cl.contig = *idbsiz - *contig;
    cl.start = *lreg;
    cl.end = *rreg;

    precon_(consensus, project_name, percd_in, idbsiz, &one, &cl,
	    task, idevr, consensus_len, maxgel, maxcon,
	    window, nbad, iladd, iradd, iok);
}

int
consensus_dialog (GapIO *io,
		  int mask_or_mark,
		  int consensus_type,
		  int output_format,
		  int gel_anno,
		  int truncate,
		  int gel_notes,
		  int use_conf,
		  int min_conf,
		  int win_size,
		  int dash,
		  char *out_file,
		  int num_contigs,
		  contig_list_t *contig_array,
		  int nopads,
		  int name_format) 
{
    int task_mask;
    int database_size, c_nconts;
    int max_read_length, c_maxseq;
    int success, consensus_length;
    float c_percd;
    int nconts;
    char *seq1; /* MAXSEQ */
    float *qual; /* MAXSEQ */
    FILE *fp;
    Hidden_params p;

    if ((seq1 = (char *)xmalloc(maxseq * sizeof(char)))==NULL){
	return(-1);
    }
    if (output_format == EXPT_FORMAT) {
	if ((qual = (float *)xmalloc(maxseq * sizeof(float)))==NULL) {
	    xfree(seq1);
	    return(-1);
	}
    } else {
	qual = NULL;
    }

    nconts = NumContigs(io);

/*  the dialogue should supply: 
    mask_or_mark 	1 means masking, 2 means mark, 0 means neither
    consensus_type	1 normal, 2 extended, 3 unfinished, 4 quality
    output_format	1 staden, 2 fasta, 3 experiment, 4 internal
    gel_anno		1 for yes
    truncate		1 for yes

    if extended then win_size and dash
    other items are set here by borrowing code from fij
*/

    p.min = p.max = p.verbose = p.qual_val = p.window_len =0;
    p.test_mode = 0;
    p.start = 0;
    p.lwin1 = 0;
    p.lcnt1 = 0;
    p.rwin1 = 0;
    p.rcnt1 = 0;

    p.use_conf = use_conf;
    p.qual_val = min_conf;
    p.rwin1    = win_size;
    p.window_len=win_size;
    p.rcnt1    = dash;

    p.band = 30; /* FIXME */

    task_mask = 1;
    if (mask_or_mark ) {
	if ( 1 == mask_or_mark ) task_mask += MASKING;
	else task_mask += MARKING;
    }
    if ( 1 == consensus_type ) task_mask += NORMALCONSENSUS;
    if ( 2 == consensus_type ) task_mask += ADDHIDDENDATA + NORMALCONSENSUS;
    if ( 3 == consensus_type ) task_mask += SINGLESTRANDED; /* = unfinished */
    if ( 4 == consensus_type ) task_mask += QUALITYCODES;

    database_size = io_dbsize(io);
    c_nconts = nconts;
    max_read_length = find_max_gel_len(io, 0, 0) + 1;
    c_maxseq = maxseq;
    c_percd = consensus_cutoff;

    if (NULL == (fp = fopen(out_file, "w"))) {
	verror(ERR_WARN, "consensus_dialogue", "%s: %s", out_file,
	       strerror(errno));
	xfree(seq1);
	if (qual) xfree(qual);
	return -1;
    }

    /* FIXME - file handling errors to be checked for */
    success = make_consensus_files ( task_mask, output_format,
				     gel_anno, truncate, gel_notes,
				     fp, io,
				     seq1, qual,
				     database_size,
				     c_nconts,
				     &consensus_length,
				     max_read_length,
				     c_maxseq,
				     p, c_percd,
				     num_contigs, contig_array,
				     nopads,
				     name_format);

    if (0 != success)
	verror(ERR_WARN, "consensus_dialog",
	       "couldn't create consensus: code %d", success);

    /* now close the file */
    fclose(fp);

    xfree(seq1);
    if (qual)
	xfree(qual);

    return 0;
}
