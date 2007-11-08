#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <io_lib/hash_table.h>
#include <io_lib/srf.h>

int main(int argc, char **argv) {
    srf_t *srf;
    char *archive, *trace;
    uint64_t cpos, hpos, dpos;
    int i;

    /*
     * We accept multiple arguments, but it's largely pointless as we just
     * concatenate the output together forming an invalid ZTR file.
     *
     * It's here for testing purposes basically so we can extract everything
     * and sum the output for validity checking.
     */
    if (argc < 3) {
	fprintf(stderr, "Usage: srf_extract archive_name trace_name ...\n");
	return 1;
    }
    archive = argv[1];

    srf = srf_open(archive, "r");

#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif

    for (i = 2; i < argc; i++) {
	trace = argv[i];

	/* Search index */
	switch (srf_find_trace(srf, trace, &cpos, &hpos, &dpos)) {
	case -1:
	    fprintf(stderr, "Malformed or missing index. "
		    "Consider running srf_index\n");
	    return 1;

	case -2:
	    fprintf(stderr, "%s: not found\n", trace);
	    break;
	
	default:
	    /* The srf object holds the latest data and trace header blocks */
	    fwrite(srf->th.trace_hdr, 1, srf->th.trace_hdr_size, stdout);
	    fwrite(srf->tb.trace,     1, srf->tb.trace_size,     stdout);
	    break;
	}
    }
	
    return 0;
}
