#ifndef _PCR_PRIMERS_H
#define _PCR_PRIMERS_H_

#include "IO.h"
#include "primlib.h"
#include "dstring.h"
#include "finish.h"

#define MAX_PRIMER_LEN 50

/* Gap4 wrapper around the primer3 primer_pair struct */
typedef struct {
    primer_pair *pair;
    int contig[2];	/* contig number for left/right primer */
    int pos[2];		/* padded contig position of left end */
    int len[2];		/* padded primer length */
    char seq[2][MAX_PRIMER_LEN+1];	/* unpadded primer sequence */
} g4_primer_pair;

dstring_t *finish_pcr_primers(finish_t *fin, char *pdefs,
			      contig_list_t *contigs, int ncontigs);

#endif /* _PCR_PRIMERS_H */
