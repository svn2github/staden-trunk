#ifndef _QUALIO_H
#define _QUALIO_H

int *count_confidence(GapIO *io, tg_rec contig, int start, int end);
int list_confidence(int *freqs, int length);
int get_base_confidences(GapIO *io, tg_rec contig, int start, int end,
			 int *match_freqs, int *mismatch_freqs,
			 long matrix[6][6]);
extern double list_base_confidence(int *matfreqs, int *misfreqs,
				   long matrix[6][6]);

#endif
