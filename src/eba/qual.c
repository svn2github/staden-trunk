#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Read.h"
#include "conf.h"
#include "open_trace_file.h"
#include "misc.h"

static int combined[] = {
    /*   0 */ 3,
    /*   1 */ 3,
    /*   2 */ 3,
    /*   3 */ 3,
    /*   4 */ 3,
    /*   5 */ 3,
    /*   6 */ 3, 
    /*   7 */ 3,
    /*   8 */ 3,
    /*   9 */ 4, 
    /*  10 */ 5, 
    /*  11 */ 5, 
    /*  12 */ 6, 
    /*  13 */ 6, 
    /*  14 */ 6, 
    /*  15 */ 6, 
    /*  16 */ 7, 
    /*  17 */ 6, 
    /*  18 */ 7, 
    /*  19 */ 7, 
    /*  20 */ 8, 
    /*  21 */ 8, 
    /*  22 */ 8, 
    /*  23 */ 9, 
    /*  24 */ 9, 
    /*  25 */ 10,
    /*  26 */ 10,
    /*  27 */ 11,
    /*  28 */ 11,
    /*  29 */ 11,
    /*  30 */ 12,
    /*  31 */ 12,
    /*  32 */ 13,
    /*  33 */ 13,
    /*  34 */ 14,
    /*  35 */ 14,
    /*  36 */ 15,
    /*  37 */ 15,
    /*  38 */ 15,
    /*  39 */ 16,
    /*  40 */ 17,
    /*  41 */ 17,
    /*  42 */ 18,
    /*  43 */ 18,
    /*  44 */ 19,
    /*  45 */ 19,
    /*  46 */ 20,
    /*  47 */ 21,
    /*  48 */ 21,
    /*  49 */ 22,
    /*  50 */ 23,
    /*  51 */ 23,
    /*  52 */ 24,
    /*  53 */ 24,
    /*  54 */ 25,
    /*  55 */ 26,
    /*  56 */ 26,
    /*  57 */ 27,
    /*  58 */ 28,
    /*  59 */ 28,
    /*  60 */ 29,
    /*  61 */ 29,
    /*  62 */ 30,
    /*  63 */ 31,
    /*  64 */ 32,
    /*  65 */ 32,
    /*  66 */ 33,
    /*  67 */ 34,
    /*  68 */ 35,
    /*  69 */ 36,
    /*  70 */ 37,
    /*  71 */ 37,
    /*  72 */ 38,
    /*  73 */ 38,
    /*  74 */ 38,
    /*  75 */ 39,
    /*  76 */ 40,
    /*  77 */ 40,
    /*  78 */ 41,
    /*  79 */ 41,
    /*  80 */ 42,
    /*  81 */ 43,
    /*  82 */ 43,
    /*  83 */ 43,
    /*  84 */ 43,
    /*  85 */ 43,
    /*  86 */ 43,
    /*  87 */ 43,
    /*  88 */ 43,
    /*  89 */ 43,
    /*  90 */ 43,
    /*  91 */ 43,
    /*  92 */ 43,
    /*  93 */ 43,
    /*  94 */ 43,
    /*  95 */ 43,
    /*  96 */ 43,
    /*  97 */ 43,
    /*  98 */ 43,
    /*  99 */ 43,
    /* 100 */ 43
};


/*
 * infp and outfp maybe the same FILE *, so we need to fseek between reading
 * and writing.
 */
static int do_it(FILE *infp, FILE *outfp, int in_f, int out_f, char *fn,
		 int phred_scale, int avg_qual, int filtered,
		 int non_filtered, int offset, int dump) {
    Read *r, *rf;

    if (NULL == (r = fread_reading(infp, fn, in_f))) {
	fprintf(stderr, "Couldn't read reading file\n");
	return 1;
    }

    if (filtered) {
	rf = read_dup(r, "?");
	filter_trace(rf);
	calc_conf_values(rf, phred_scale, 1, offset);

	if (phred_scale) {
	    rescale_scores(rf, 1);
	}
    }
    if (non_filtered) {
	calc_conf_values(r, phred_scale, 0, offset);

	/* Average confidence values */
	if (avg_qual)
	    average_conf(r);

	if (phred_scale) {
	    rescale_scores(r, 0);
	}
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
		double p1, p2, p;
		c2 =  f_conf[bind][i];
		/* Tried MIN(c1,c2) */
		/* f_conf[bind][i] = (c1 * c2) / 100;*/
		p1 = 1-pow(10,(double)c1/-10.0);
		p2 = 1-pow(10,(double)c2/-10.0);
		p = (p1*p2) / (p1*p2 + (1-p1)*(1-p2));
		f_conf[bind][i] = combined[f_conf[bind][i]];
	    } else {
		f_conf[bind][i] = c1;
	    }
	}

	read_deallocate(rf);
    }

    rewind(outfp);
    if (r->format == TT_CTF || r->format == TT_ZTR)
	out_f = r->format;

    if (dump) {
	int i, count = 20;
	switch(dump) {
	case 2:
	    printf(">%s\n", fn);
	    break;
	case 3:
	    printf("BaseQuality : %s\n", fn);
	    break;
	}

	for (i = 0; i < r->NBases; i++) {
	    char conf;
	    switch(r->base[i]) {
	    case 'A': case 'a':
		conf = r->prob_A[i];
		break;
	    case 'C': case 'c':
		conf = r->prob_C[i];
		break;
	    case 'G': case 'g':
		conf = r->prob_G[i];
		break;
	    case 'T': case 't':
		conf = r->prob_T[i];
		break;
	    default:
		conf = 0;
	    }
	    if (dump == 1) {
		putchar(conf);
	    } else {
		if (dump == 4)
		    printf("%d%c", i, r->base[i]);
		printf("%d%c", conf, --count == 0 ? '\n' : ' ');
		if (!count)
		    count = 20;
	    }
	}

	if (count != 20)
	    putchar('\n');

	if (dump == 3)
	    putchar('\n');

    } else {
	if (-1 == (fwrite_reading(outfp, r, out_f))) {
	    fprintf(stderr, "Couldn't write reading file\n");
	    read_deallocate(r);
	    return 1;
	}

	ftruncate(fileno(outfp), ftell(outfp));
    }

    read_deallocate(r);

    return 0;
}

static void usage(void) {
    fprintf(stderr, "Usage: eba [options] [trace_file]\n");
    fprintf(stderr, "   -phred_scale                Use phred log scale\n");
    fprintf(stderr, "   -old_scale                  Use S/N ratios\n");
    fprintf(stderr, "   -average 0/1                Whether to avg non-filtered results\n");
    fprintf(stderr, "   -non_filtered 0/1           Compute S/N on non-filtered traces\n");
    fprintf(stderr, "   -filtered 0/1               Compute S/N on filtered traces\n");
    fprintf(stderr, "   -offset value               Add value to denominator in S/N calc.\n");
    fprintf(stderr, "   -dump raw/fasta/caf         Output quality to stdout instead of new trace.\n");
    fprintf(stderr, "\n  eg. eba -phred_scale -non_filtered 1 -average 1 -filtered 1 -offset 50 a.scf\n");

    exit(1);
}

/*
 * Experimentation on /nfs/repository/p444/dJ1050E16/ (approx 2000 sequences
 * forming a finished 90Kb contig) show that combining averaged non_filtered
 * (ie old eba) with filtered (new eba) gives the best amount of
 * discrepancy between good and bad base calls. Compared to phred this seems
 * to discriminate better for poor data and not so well on very good data.
 */
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
    int offset = 20;
    int dump = 0;

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
	else if (strcmp(argv[a], "-offset") == 0)
	    offset = atoi(argv[++a]);
	else if (strcmp(argv[a], "-dump") == 0) {
	    ++a;
	    if (strcmp(argv[a], "raw") == 0)
		dump = 1;
	    else if (strcmp(argv[a], "fasta") == 0)
		dump = 2;
	    else if (strcmp(argv[a], "caf") == 0)
		dump = 3;
	    else if (strcmp(argv[a], "debug") == 0)
		dump = 4;
	    else
		usage();
	}
	else if (strcmp(argv[a], "-h") == 0)
	    usage();
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
		 avg_qual, filtered, non_filtered, offset, dump);
}
