#ifndef _PROBE_H_
#define _PROBE_H_

#include <tcl.h>
#include "IO.h"

/*
 * Find oligos suitable for the 'probe' sequencing strategy. The strategy
 * is to find oligos near the end of the contigs pointing towards the centre
 * of the contig. These oligos are then used in a screening process to find
 * templates containing these oligos. The forward readings for each template
 * is then sequenced and assembled. Note that the oligos chosen therefore need
 * to be unique in the consensus sequence and not found in the vector
 * sequences.
 *
 * This routine lists suitable oligos. Arguments are:
 *
 * min_size, max_size
 *	The range of lengths of suitable oligos
 *
 * match_perc
 *	Oligos matching other sequence with >= match_perc are considered
 *	non-unique.
 *
 * from, to
 *	Which area of contig consensus to search for oligos in. Specified as
 *	offsets from the contig ends.
 *
 * vectors
 *	A file of filenames of vector sequences to screen against.
 *
 * Returns:  0 for success
 *           -1 for failure
 */
int find_probes(GapIO *io, int num_contigs, int *contig_arr,
		int min_size, int max_size, float match_perc,
		int from, int to, char *vectors, char *primer_defs,
		Tcl_DString *dstr);

#endif
