#ifndef _CHECK_ASSEMBLY_H_
#define _CHECK_ASSEMBLY_H_

/*
 * Scans through one or more contigs checking each reading for correct
 * assembly.
 * This is done by obtaining the cutoff data and aligning this with
 * the consensus sequence.
 *
 * Returns -1 for failure, 0 for success.
 */
int check_assembly(GapIO *io, int num_contigs, int *contigs, int cutoff,
		   int minlen, int winsize, int maxdash, float maxperc);

#endif
