#ifndef _QUALIO_H
#define _QUALIO_H

int *count_confidence(GapIO *io, int contig, int start, int end);
int list_confidence(int *freqs, int length);

#endif
