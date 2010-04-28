#include <staden_config.h>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "seqInfo.h"
#include "consen.h"
#include "xalloc.h"

/* johnt 1/6/99 must explicitly import globals from DLLs with Visual C++*/
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif
 
extern DLL_IMPORT char *optarg;
extern DLL_IMPORT int optind;

typedef struct {
    /* For both clipping methods */
    int min;		/* minimum 5' clip point */
    int max;		/* maximum 3' clip point */
    int verbose;
    int use_conf;	/* which method to use, 1 => confidence */
    int test_mode;	/* 1 => do not write out changes */
    int min_len;	/* minimum length */

    /* For N-count clipping */
    int start;		/* Start point for scanning left/right. */
    int lwin1, lcnt1;	/* 1st left clip window length and number of N's */
    int lwin2, lcnt2;	/* 2nd left clip window length and number of N's */
    int rwin1, rcnt1;	/* 1st right clip window length and number of N's */
    int rwin2, rcnt2;	/* 2nd right clip window length and number of N's */

    /* For confidence value clipping */
    int qual_val;	/* average quality value */
    int window_len;	/* over this window length */
} params;


/*
 * Scans through a quality buffer finding the highest average block of length
 * window_len.
 */
static int find_highest_conf(params p, int1 *conf, int len) {
    int i, total, best_total, best_pos;

    if (p.window_len >= len)
	return len/2;

    for (total = i = 0; i < p.window_len; i++) {
	total += conf[i];
    }
    best_total = total;
    best_pos = 0;

    for (i = 0; i < len - p.window_len; i++) {
	total = total - conf[i] + conf[i+p.window_len];
	if (total > best_total) {
	    best_total = total;
	    best_pos = i;
	}
    }

    if (p.verbose)
	printf("    Start position = %d (total %d)\n",
	       best_pos+1, best_total);

    return best_pos+1;
}


/*
 * Scans leftwards in a quality buffer until the average quality drops below
 * a specific threshold, defined by the global qual_val and window_len
 * parameters.
 *
 * Having found this window, the procedure repeats with successively smaller
 * windows until the exact base is identified.
 */
int scan_left(params p, int1 *conf, int start_pos) {
    int i, total, lclip;
    int lowest_total;
    int win_len = p.window_len;

    do {
	lowest_total = p.qual_val * win_len;

	total = 0;
	for (i = start_pos; i < start_pos + win_len; i++)
	    total += conf[i];

	i = start_pos;
	do {
	    i--;
	    total = total + conf[i] - conf[i+win_len];
	} while (i > 0 && total >= lowest_total);

	start_pos = i+2;
    } while (--win_len > 0);

    lclip = i;
    if (p.verbose)
	printf("    left clip = %d\n", lclip);

    return lclip;
}


/*
 * Scans rightwards in a quality buffer until the average quality drops below
 * a specific threshold, defined by the global qual_val and window_len
 * parameters.
 *
 * Having found this window, the procedure repeats with successively smaller
 * windows until the exact base is identified.
 */
int scan_right(params p, int1 *conf, int start_pos, int len) {
    int i, total, rclip;
    int lowest_total;
    int win_len = p.window_len;

    do {
	lowest_total = p.qual_val * win_len;

	total = 0;
	for (i = start_pos; i < start_pos + win_len; i++)
	    total += conf[i];

	i = start_pos;
	do {
	    total = total - conf[i] + conf[i+win_len];
	    i++;
	} while (i <= (len - win_len - 1) && total >= lowest_total);

	start_pos = i-1;
    } while (--win_len > 0);

    rclip = i == len ? len + 1 : i;
    if (p.verbose)
	printf("    right clip = %d\n", rclip);

    return rclip;
}


/*
 * Quality clips file 'file', updating the QL and QR records in the process.
 *
 * Returns 0 for success.
 *        -1 for failure.
 */
static int qclip(char *file, params p) {
    SeqInfo *si = NULL;
    int start_pos, right_pos, left_pos;
    FILE *fp;

    if (p.verbose)
	printf("Clipping file %s\n", file);

    /* Read the sequence and confidence */
    if (NULL == (si = read_sequence_details(file, 0))) {
	fprintf(stderr, "Failed to read file '%s'\n", file);
	return -1;
    }

    if (p.use_conf) {
	int i;
	int1 *conf;

	conf = xmalloc(si->length * sizeof(*conf));
	if (SeqInfo_conf(si, conf, si->length) == -1) {
	    if (p.verbose)
		puts("    Could not load confidence values - using sequence");
	    p.use_conf = 0;
	}

	for (i = 0; i < si->length; i++)
	    if (conf[i] != 0)
		break;
	if (i == si->length) {
	    if (p.verbose)
		puts("    Confidence values are all zero - using sequence");
	    p.use_conf = 0;
	}

	if (p.use_conf) {
	    /* Identify the best location to start from */
	    start_pos = find_highest_conf(p, conf, si->length);

	    /* Scan left, and scan right */
	    left_pos = scan_left(p, conf, start_pos)+1;
	    right_pos = scan_right(p, conf, start_pos, si->length)+1;
	}
	xfree(conf);
    }

    if (!p.use_conf) {
	char *seq = exp_get_entry(si->e, EFLT_SQ);

	left_pos = start_of_good(seq, p.start, p.lwin1, p.lcnt1,
				 p.lwin2, p.lcnt2) - 1;
	right_pos = end_of_good(seq, p.start, p.rwin1, p.rcnt1,
				p.rwin2, p.rcnt2) + 1;
    }
    
    if (left_pos < p.min)
	left_pos = p.min;
    if (right_pos > p.max)
	right_pos = p.max;
    if (left_pos >= right_pos)
	left_pos = right_pos-1;

    if (right_pos - left_pos < p.min_len) {
	fprintf(stderr, "Sequence too short (length=%d)\n",
		right_pos - left_pos);
	freeSeqInfo(si);
	return -1;
    }

    /* Append details onto the end of the Exp File */
    if (!p.test_mode) {
	if (NULL == (fp = fopen(file, "a"))) {
	    fprintf(stderr, "Failed to write file '%s'\n", file);
	    freeSeqInfo(si);
	    return -1;
	}
	fprintf(fp, "QL   %d\n", left_pos);
	fprintf(fp, "QR   %d\n", right_pos);
	fclose(fp);
    } else {
	printf("%-30s QL %4d            QR %4d\n", file, left_pos, right_pos);
    }

    freeSeqInfo(si);
    return 0;
}

static void usage(void) {
    fprintf(stderr,
	"Usage for using confidence codes (default mode):\n"
	"       qclip [-c] [-vt] [-m min 5' cutoff] [-M max 3' cutoff] [-x min_length]\n"
	"                  [-w window_len(30)] [-q average_quality (10)] file ...\n\n"
	"Usage for using sequence only:\n"
	"       qclip -n   [-vt] [-m min 5' cutoff] [-M max 3' cutoff] [-x min_length]\n"
	"                  [-s start_offset(70)]\n"
	"                  [-L left_window_len(20)] [-l left_N_count(3)]\n"
	"                  [-R right_window_len(100)] [-r right_N_count(5)] file ...\n"
	    );
    exit(1);
}

int main(int argc, char **argv) {
    int c, i, ret = 0;
    params p;

    /* Defaults */
    p.min = 0;
    p.max = 10000000;
    p.min_len = 0;
    p.verbose = 0;
    p.use_conf = 1;
    p.start = 70;
    p.lwin1 = 20;
    p.lcnt1 = 3;
    p.lwin2 = 0;
    p.lcnt2 = 0;
    p.rwin1 = 100;
    p.rcnt1 = 5;
    p.rwin2 = 0;
    p.rcnt2 = 0;
    p.qual_val = 10;
    p.window_len = 30;
    p.test_mode = 0;

    while ((c = getopt(argc, argv, "q:w:vtncm:M:R:r:L:l:s:x:")) != -1) {
	switch (c) {
	    /* Both methods */
	case 'v':
	    p.verbose = 1;
	    break;

	case 't':
	    p.test_mode = 1;
	    break;

	case 'm':
	    p.min = atoi(optarg);
	    break;

	case 'M':
	    p.max = atoi(optarg);
	    break;

	case 'x':
	    p.min_len = atoi(optarg);
	    break;

	case 'n':
	    p.use_conf = 0;
	    break;

	case 'c':
	    p.use_conf = 1;
	    break;

	    /* New method */
	case 'q':	    
	    p.qual_val = atoi(optarg);
	    break;

	case 'w':
	    p.window_len = atoi(optarg);
	    break;

	    /* Old method */
	case 'R':
	    p.rwin1 = atoi(optarg);
	    break;

	case 'r':
	    p.rcnt1 = atoi(optarg);
	    break;

	case 'L':
	    p.lwin1 = atoi(optarg);
	    break;

	case 'l':
	    p.lcnt1 = atoi(optarg);
	    break;

	case 's':
	    p.start = atoi(optarg);
	    break;

	default:
	    usage();
	}
    }

    if (optind == argc)
	usage();

    for (i = optind; i < argc; i++) {
	int ret_val;

	ret_val = qclip(argv[i], p);
	if (p.verbose)
	    printf("    qclip() returned %d\n", ret_val);

	ret |= ret_val;
    }

    return ret ? 1 : 0;
}
