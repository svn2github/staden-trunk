#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "Read.h"
#include "conf.h"
#include "open_trace_file.h"
#include "misc.h"


/*
 * infp and outfp maybe the same FILE *, so we need to fseek between reading
 * and writing.
 */
static int do_it(FILE *infp, FILE *outfp, int in_f, int out_f, char *fn,
		 int phred_scale, int avg_qual,
		 int filtered, int non_filtered) {
    Read *r, *rf;

    if (NULL == (r = fread_reading(infp, fn, in_f))) {
	fprintf(stderr, "Couldn't read reading file\n");
	return 1;
    }

    if (filtered) {
	rf = read_dup(r, "?");
	calc_conf_values(rf, phred_scale, avg_qual, 1);
	if (phred_scale) {
	    rescale_scores(rf);
	}
    }
    if (non_filtered)
	calc_conf_values(r, phred_scale, avg_qual, 0);

    /* Average confidence values */
    if (avg_qual)
	average_conf(r);

    /* Move to phred scale. Only tuned on non-filtered non-'cosa' values */
    if (phred_scale) {
	rescale_scores(r);
    }

    /* Merge filtered and non-filtered, or just copy filtered confidences */
    if (filtered) {
	int i;
	char *nf_conf[4], *f_conf[4], base[256];

	nf_conf[0] = rf->prob_A;
	nf_conf[1] = rf->prob_C;
	nf_conf[2] = rf->prob_G;
	nf_conf[3] = rf->prob_T;
	f_conf[0] = r->prob_A;
	f_conf[1] = r->prob_C;
	f_conf[2] = r->prob_G;
	f_conf[3] = r->prob_T;

	memset(base, 5, 256);
	base['A'] = 0;
	base['a'] = 0;
	base['C'] = 1;
	base['c'] = 1;
	base['G'] = 2;
	base['g'] = 2;
	base['T'] = 3;
	base['t'] = 3;

	for (i = 0; i < r->NBases; i++) {
	    char c1, c2;
	    int bind;

	    bind = base[r->base[i]];
	    if (bind == 5) {
		f_conf[0][i] = f_conf[1][i] = f_conf[2][i] = f_conf[3][i] = 0;
		continue;
	    }

	    c1 = nf_conf[bind][i];
	    if (non_filtered) {
		c2 =  f_conf[bind][i];
		f_conf[bind][i] = MIN(c1,c2);
	    } else {
		f_conf[bind][i] = c1;
	    }
	}

	read_deallocate(rf);
    }

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
    int avg_qual = 1;
    int filtered = 0;
    int non_filtered = 1;

    while (a < argc) {
	if (strcmp(argv[a], "-phred_scale") == 0)
	    phred_scale = 1;
	else if (strcmp(argv[a], "-old_scale") == 0)
	    phred_scale = 0;
	else if (strcmp(argv[a], "-average") == 0)
	    avg_qual = atoi(argv[++a]);
	else if (strcmp(argv[a], "-non_filtered") == 0)
	    non_filtered = atoi(argv[++a]);
	else if (strcmp(argv[a], "-filtered") == 0)
	    filtered = atoi(argv[++a]);
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

    return do_it(ifp, ofp, in_type, out_type, fn, phred_scale,
		 avg_qual, filtered, non_filtered);
}
