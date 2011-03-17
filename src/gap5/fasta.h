#ifndef _FASTA_H_
#define _FASTA_H_

#include <tg_index.h>

int parse_fasta_or_fastq(GapIO *io, char *fn, tg_args *a, int format);

#endif /* _FASTA_H_ */
