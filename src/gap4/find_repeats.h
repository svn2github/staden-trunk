#include "io_utils.h"
#include "IO.h"

int
find_repeats(GapIO *io,
	     int handle,
             int idir,
             int minmat,
	     f_int mask,
	     float percd,
	     int num_contigs,
	     contig_list_t *clist,
	     char *outfile);

