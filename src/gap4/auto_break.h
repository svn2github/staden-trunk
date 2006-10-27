#ifndef AUTO_BREAK_H
#define AUTO_BREAK_H

#include "dstring.h"

dstring_t *auto_break_contigs(GapIO *io, int argc, contig_list_t *argv,
			      double filter_score, int by_consensus);

#endif /* AUTO_BREAK_H */
