#ifndef _CONSEN_H_
#define _CONSEN_H_

#include "io_utils.h"

/*****************************************************************************
 * Macros
 ****************************************************************************/

/* Jobs for maskit */
#define STANDARD_TO_MASKED 1
#define MASKED_TO_MARKED   2
#define MARKED_TO_MASKED   3
#define STANDARD_TO_MARKED 4

/* Jobs for precon */
#define ADDTITLE 1
#define ADDHIDDENDATA 2
#define NORMALCONSENSUS 4
#define SINGLESTRANDED 8
#define QUALITYCODES 16
#define MASKING 32
#define MARKING 64
#define SORTCONTIGS 128

#define RIGHT_END 1
#define LEFT_END 2

#define CONS_NAME_LREADING 1
#define CONS_NAME_LTEMPLATE 2

/*****************************************************************************
 * Typedefs & structures
 ****************************************************************************/

typedef struct contig_parms_ {
    tg_rec contig_number;
    int contig_start;
    int contig_end;
    tg_rec contig_left_gel;
    int contig_left_extension;
    int contig_right_extension;
    int contig_start_offset;
    int contig_end_offset;
} Contig_parms;


typedef struct {

    /* do it or not? */
    int do_it;

    /* For both clipping methods */
    int min;		/* minimum 5' clip point */
    int max;		/* maximum 3' clip point */
    int verbose;
    int use_conf;	/* which method to use, 1 => confidence */
    int test_mode;	/* 1 => do not write out changes */

    /* For N-count clipping */
    int start;		/* Start point for scanning left/right. */
    int lwin1, lcnt1;	/* 1st left clip window length and number of N's */
    int rwin1, rcnt1;	/* 1st right clip window length and number of N's */

    /* For confidence value clipping */
    int qual_val;	/* average quality value */
    int window_len;	/* over this window length */

    /* For alignment */
    int gap_open;
    int gap_extend;
    int band;
} Hidden_params;

/*****************************************************************************
 * Function prototypes
 ****************************************************************************/

/* Setup masking functions */
void set_mask_lookup(void);



/*
   Add title to a consensus. 

   Inputs: 
   
   1. consensus:     20 char chunk of consensus array
   2. project_name:  project name terminated by .
   3. left_gelnumber: the contig left gel number

   Output:

   1. Title is 20 chars in length of form: 

   <project_name.left_gelnumber-->
*/
void add_contig_title(char *consensus, char *project_name, tg_rec left_gelnumber);



/*	sort contigs on length

	input:  a list of contig numbers.
	output: the sorted list.
*/
int sort_contigs ( Contig_parms *contig_list, 
		  int number_of_contigs);



/*
 * Given an array of contig numbers (1..Num) or a single contig
 * ((8000-1)..(8000-Num)), generate a Contig_parms array containing
 * contig number, start, end, and left gel.
 */
Contig_parms *get_contig_list (int database_size, GapIO *io,
			       int number_of_contigs,
			       contig_list_t *contig_array);



/*	routine to do masking and marking	*/
void maskit ( char *seq, int seq_length, int job);
void maskc_ (char *seq, f_int *seq_len, f_int *jobin, f_implicit seq_l);



    /* Routine to mask regions of a consensus	*/

    /* *consensus		the consensus
       *tag_type_list		the list of tag types
       number_of_tag_types	the number of tag types in the list
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
int mask_consensus(GapIO *io, char *consensus, int contig, int lreg, int rreg, 
		   int job);



void precon_ ( char *consensus, char *project_name, 
	     float *percd,
	     int *idbsiz, int *num_contigs, contig_list_t *clist,
	      int *task, int *idevr,
	      int *consensus_len, int *maxgel, int *maxcon,
	      int *window, int *nbad, int *iladd, int *iradd, 
	      int *iok);

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
		  int name_format); 

char *unattached_reads(GapIO *io);

char *minimal_coverage(GapIO *io, int num_c, contig_list_t *contigs);

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

   Return values: 0 success
                 -1 insufficient space
                 -2 masking failed
		 -3 add hidden failed
*/
int make_consensus( int task_mask, GapIO *io,
		   char **consensus, float *quality,
		   Contig_parms *contig_list, int number_of_contigs,
		   int *consensus_length, int max_read_length,
		   Hidden_params p, float percd );


/* 
 * read through the consensus and find the ends of the contigs.
 * store their positions in contig_ends[] and 
 * their left gel numbers in contig_numbers[]
 *
 * NOTE: contig_ends[] and contig_numbers[] should be allocated to
 * the number of contigs PLUS ONE.
 */

int find_contig_ends ( char *seq, int seq_len, 
		       int *contig_ends, tg_rec *contig_numbers );


int end_of_good ( char *seq, 
		 int start,
		 int window_len1, 
		 int max_unknown1);

/* given a position in a consensus sequence return the contig left gel num
 * ASSUME that contig_start contains array element in consensus of left end 
 * of contig!! and that the contig_list is otherwise correctly filled in.
 */

int contig_listel_from_con_pos ( Contig_parms *contig_list, 
				 int number_of_contigs, int pos_in_contig );

#endif
