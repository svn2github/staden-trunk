#ifndef _SHUFFLE_PADS_H
#define _SHUFFLE_PADS_H

#include <tg_gio.h>
#include "io_utils.h"

int shuffle_contigs_io(GapIO *io, int ncontigs, contig_list_t *contigs,
		       int band);


#endif /* _SHUFFLE_PADS_H */
