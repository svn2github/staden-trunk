#ifndef _QUALIO_H
#define _QUALIO_H

int *count_confidence(GapIO *io, int contig, int start, int end);
int list_confidence(int *freqs, int length);
int get_base_confidences(GapIO *io, int contig,
			 int *match_freqs, int *mismatch_freqs);
extern double list_base_confidence(int *matfreqs, int *misfreqs);

#endif
