#ifndef _FIND_FRAGMENTS_H
#define _FIND_FRAGMENTS_H

#include "qual.h"

/*
 * Structure type used in find_fragments() below.
 */
typedef struct {
    int num;		/* Reading number */
    int abs_start;	/* Abs_start and abs_end are an inclusive range */
    int abs_end;	/*   of absolute positions in the contig (1..n) */
    int seq_start;	/* Seq_start and seq_end are an inclusive range */
    int seq_end;	/*   relative to the each sequence. */
    int cutoff_len;	/* Length of left hand cutoff data */
} seq_frag;

/*
 * find_fragments()
 *
 * This breaks down the assembly between start and end (inclusive, starting
 * from 1) into small fragments where each fragment consists of a set of
 * sequences which all span the entire fragment size.
 *
 * Or to put it another way:
 *
 * Seq1  ----------------------------------
 * Seq2      ------------------------
 * Seq3           --------------------------------
 *
 * Gives:
 *
 * Seq1  |----|-----|-------------------|------|       |
 * Seq2  |    |-----|-------------------|      |       |
 * Seq3  |    |     |-------------------|------|-------|
 *
 * This greatly simplifies the process of certain types of depth analysis.
 * See find_fragments.c for an example usage.
 *
 * Arguments:
 *	io		The Gap IO handle
 *	contig		The contig number (1..NumContigs)
 *	start		The first base to include (1..ContigLen)
 *	end		The last base to include (start..ContigLen)
 *	info_func	A function to call when querying contig and seq info
 *	info_data	Client specific data passed into info_func
 *	clientfunc	A function to call with each fragment set
 *	clientdata	Client specific data passed into callback
 *
 * Returns:
 *	 0 for success
 *	-1 for failure
 */
int find_fragments(GapIO *io,
		   int contig,
		   int start,
		   int end,
		   int (*info_func)(int         job,
				    void       *mydata,
				    info_arg_t *theirdata),
		   void *info_data,
		   void (*clientfunc)(GapIO *io,
				      int contig,
				      int start,
				      int end,
				      seq_frag *frag,
				      int num_frags,
				      void *clientdata),
		   void *clientdata);

#endif /* _FIND_FRAGMENTS_H */
