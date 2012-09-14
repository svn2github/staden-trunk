#ifndef _SHUFFLE_PADS_H
#define _SHUFFLE_PADS_H

#include <tg_gio.h>
#include "io_utils.h"

int shuffle_contigs_io(GapIO *io, int ncontigs, contig_list_t *contigs,
		       int band, int flush);

int remove_pad_columns(GapIO *io, int rargc, contig_list_t *rargv,
		       int percent_pad, int quiet);

#endif /* _SHUFFLE_PADS_H */
