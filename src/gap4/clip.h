#ifndef _CLIP_H_
#define _CLIP_H_

#include "IO.h"
#include "io_utils.h"

void quality_clip(GapIO *io, int num_contigs, contig_list_t *cl, int qual_avg);
void difference_clip(GapIO *io, int num_contigs, contig_list_t *cl,
		     int add_tags);

#endif /* _CLIP_H_ */
