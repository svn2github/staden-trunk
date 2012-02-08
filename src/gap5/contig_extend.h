#ifndef _CONTIG_EXTEND_H_
#define _CONTIG_EXTEND_H_

#include <tg_gio.h>

int contig_extend(GapIO *io, tg_rec *contigs, int ncontigs, int min_depth,
		  int match_score, int mismatch_score);

#endif /* _CONTIG_EXTEND_H_ */

