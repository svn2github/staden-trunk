#include <stdlib.h>
#include <limits.h>

#include "IO.h"
#include "misc.h"
#include "find_fragments.h"

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
 *
 * Arguments:
 *	io		The Gap IO handle
 *	contig		The contig number (1..NumContigs)
 *	start		The first base to include (1..ContigLen)
 *	end		The last base to include (start..ContigLen)
 *	info_func	A function to call when querying contig and seq info
 *	info_data	Client specific data passed into info_func
 *	clientfunc	A function to call with each fragment set
 *	clientdata	Client specific data passed into clientfunc
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
		   void *clientdata) {
    seq_frag *frag;
    int max_frags, num_frags;
    int last_pos = start;
    int next_end = INT_MAX;
    info_arg_t info;

    /*
     * Firstly, find the left most gel in our region.
     * This leaves the first GReadings spanning this region in "r".
     */
    info.contig_info.contig = contig;
    info_func(GET_CONTIG_INFO, info_data, &info);
    info.gel_info.gel = info.contig_info.leftgel;
    do {
	info_func(GET_GEL_INFO, info_data, &info);
    } while (info.gel_info.position + info.gel_info.length < start &&
	     (info.gel_info.gel = info.gel_info.next_right));

    /*
     * Initialise our "frag" array.
     * See qual.c:calc_contig_info_phred() for similar code structuring.
     */
    max_frags = 10;
    num_frags = 0;
    if (NULL == (frag = (seq_frag *)xmalloc(max_frags * sizeof(*frag))))
	return -1;

    /*
     * Now we scan through sequences finding regions of identical depth.
     * We produce a series of sequence fragments within this region.
     * This makes non-trivial depth-based analysis (such as depth of
     * templates) much easier.
     */
    do {
	int frag_start, frag_end;

	if (info.gel_info.gel == 0)
	    break;

	/* Keep track of the nearest sequence end point */
	if (next_end > info.gel_info.position + info.gel_info.length - 1) {
	    next_end = info.gel_info.position + info.gel_info.length - 1;
	}

	/* Add gel to fragment list */
	if (num_frags >= max_frags) {
	    max_frags *= 2;
	    frag = (seq_frag *)xrealloc(frag, max_frags * sizeof(*frag));
	    if (frag == NULL)
		return -1;
	}
	
	frag[num_frags].num = info.gel_info.gel;
	frag[num_frags].abs_start = info.gel_info.position;
	frag[num_frags].abs_end =
	    info.gel_info.position + info.gel_info.length - 1;
	frag[num_frags].cutoff_len = info.gel_info.start;
	num_frags++;


	/*
	 * Candidate end position - where the next sequence starts.
	 */
	last_pos = info.gel_info.position;
	info.gel_info.gel = info.gel_info.next_right;
	if (info.gel_info.gel) {
	    info_func(GET_GEL_INFO, info_data, &info);
	    frag_end = MIN(info.gel_info.position-1, end);
	} else {
	    frag_end = end;
	}

	frag_start = MAX(last_pos, start);

	if (frag_end >= frag_start) {
	    int i;
	    int tmp_end;

	    /*
	     * We have a fragment, but we may need to break it up further.
	     * Specifically, our frag_start to frag_end at this point will
	     * indicate the region between two sequence start points, but
	     * one or more sequences could end within this range. We keep
	     * track of the next end position and break this fragment into
	     * sub fragments.
	     */
	    do {
		tmp_end = MIN(frag_end, next_end);

		if (tmp_end >= frag_start) {
		    /*
		     * At long last, we now have our fragment information.
		     * Compute the start and end points within the individual
		     * sequences and call our callback function.
		     */
		    for (i = 0; i < num_frags; i++) {
			frag[i].seq_start =frag_start - frag[i].abs_start +
			    frag[i].cutoff_len;
			frag[i].seq_end   = tmp_end - frag[i].abs_start +
			    frag[i].cutoff_len;
		    }

		    clientfunc(io, contig, frag_start, tmp_end,
			       frag, num_frags, clientdata);
		}
		
		/*
		 * Check if a sequence ends - if so, remove it from our frag
		 * array.
		 * This also recomputes the new next_end value.
		 */
		frag_start = MAX(start, next_end+1);
		next_end = INT_MAX-1;
		for (i = 0; i < num_frags; i++) {
		    if (frag[i].abs_end <= tmp_end) {
			memmove(&frag[i], &frag[i+1],
				(num_frags-i-1) * sizeof(*frag));
			num_frags--; i--;
		    } else if (next_end > frag[i].abs_end) {
			next_end = frag[i].abs_end;
		    }
		}
	    } while (frag_start <= frag_end);
	}

    } while (info.gel_info.position <= end);

    xfree(frag);
    return 0;
}

#if 0
/*
 * Here's an example usage of find_fragments.
 */
static void find_fragments_callback(GapIO *io, int contig, int start, int end,
				    seq_frag *frag, int num_frags,
				    void *clientdata) {
    int i, j;
    GReadings r;
    char *seq;

    printf("Fragment from %d to %d\n", start, end);
    for (i = 0; i < num_frags; i++) {
	printf("\tSeq %4d, pos %4d to %4d ",
	       frag[i].num, frag[i].seq_start, frag[i].seq_end);

	gel_read(io, frag[i].num, r);
	seq = TextAllocRead(io, r.sequence);
	putchar('\t');
	for (j = frag[i].seq_start; j <= frag[i].seq_end; j++)
	    putchar(seq[j]);
	putchar('\n');
	xfree(seq);
    }
}

int find_fragments_main(GapIO *io, int contig, int start, int end) {
    int r;

    r = find_fragments(io, contig, start, end, database_info, io,
		       find_fragments_callback, NULL);
    return r;
}
#endif
