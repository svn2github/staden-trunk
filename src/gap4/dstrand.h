#ifndef _DSTRAND_H_
#define _DSTRAND_H_

void dbl_complement(GapIO *io, int *lreg, int *rreg, int contig);

void double_strand_list(GapIO *io, int num_contigs,
			contig_list_t *contigs,
			int miscount, float misperc);

void double_strand_single(GapIO *io, int contig, int lreg, int rreg,
			  int miscount, int misperc);

#endif
