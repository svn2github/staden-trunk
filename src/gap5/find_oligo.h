#include "io_utils.h"

int find_oligos(GapIO *io,
		int num_contigs,
		contig_list_t *contig_array,
		float mis_match,
		char *string,
		int consensus_only,
		int in_cutoff);

int
find_oligo_file(GapIO *io,
		int num_contigs,
		contig_list_t *contig_array,
		float mis_match,
		char *file,
		int consensus_only,
		int in_cutoff);
