#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "Read.h"
#include "conf.h"
#include "open_trace_file.h"

/*
 * infp and outfp maybe the same FILE *, so we need to fseek between reading
 * and writing.
 */
static int do_it(FILE *infp, FILE *outfp, int in_f, int out_f, char *fn,
		 int phred_scale) {
    Read *r;

    if (NULL == (r = fread_reading(infp, fn, in_f))) {
	fprintf(stderr, "Couldn't read reading file\n");
	return 1;
    }

    calc_conf_values(r, phred_scale);

    rewind(outfp);
    if (r->format == TT_CTF || r->format == TT_ZTR)
	out_f = r->format;

    if (-1 == (fwrite_reading(outfp, r, out_f))) {
	fprintf(stderr, "Couldn't write reading file\n");
	read_deallocate(r);
	return 1;
    }

    ftruncate(fileno(outfp), ftell(outfp));
    read_deallocate(r);

    return 0;
}

static void usage(void) {
    fprintf(stderr, "Usage: eba [trace_file]\n");

    exit(1);
}

int main(int argc, char **argv) {
    FILE *ifp = stdin;
    FILE *ofp = stdout;
    char *fn;
    int in_type = TT_ANY, out_type = TT_ANY;
    int phred_scale = 1;
    int a = 1;

    while (a < argc) {
	if (strcmp(argv[a], "-phred_scale") == 0)
	    phred_scale = 1;
	else if (strcmp(argv[a], "-old_scale") == 0)
	    phred_scale = 0;
	else
	    break;
	a++;
    }

    switch (argc-a) {
    case 0: {
	fn = "(stdin)";
	break;
    }
    case 1: {
	ifp = fopen_compressed(argv[a], &ofp);
	if (ifp == NULL) {
	    perror(argv[a]);
	    return 1;
	}
	fn = argv[a];
	break;
    }
    default:
	usage();
	return 1;
    }

    return do_it(ifp, ofp, in_type, out_type, fn, phred_scale);
}
