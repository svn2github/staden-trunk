#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#include "Read.h"
#include "misc.h"

/* johnt 1/6/99 must explicitly import globals from DLLs with Visual C++*/
#ifdef _MSC_VER
#  define DLL_IMPORT __declspec(dllimport)
#else
#  define DLL_IMPORT
#endif
 
extern DLL_IMPORT char *optarg;
extern DLL_IMPORT int optind;

typedef struct {
    int verbose;
    int signal_val;	/* average signal strength */
    int window_len;	/* over this window length */
    double threshold;
} params;


/*
 * Scans through a quality buffer finding the highest average block of length
 * window_len.
 */
static int find_highest_peak(params *p, int *peaks, int len, int *avg_height) {
    int i, total, best_total, best_pos;

    if (p->window_len >= len)
	return len/2;

    for (total = i = 0; i < p->window_len; i++) {
	total += peaks[i];
    }
    best_total = total;
    best_pos = 0;

    for (i = 0; i < len - p->window_len; i++) {
	total = total - peaks[i] + peaks[i+p->window_len];
	if (total > best_total) {
	    best_total = total;
	    best_pos = i;
	}
    }

    if (p->verbose)
	printf("    Start position = %d (total %d)\n",
	       best_pos+1, best_total);

    if (avg_height)
	*avg_height = best_total / (double)p->window_len + 0.5;

    return best_pos+1;
}

#define PI 3.1415926535897932
int scan_right(char *name, params *p, int *peaks, int start_pos,
	       int highest, int len, int level) {
    int i;
    int lowest_total;
    int win_len2 = p->window_len/2;
    int total_l = 0, total_r = 0;
    double best_val = 0;
    int peakpos[10];
    double peakval[10];
    int npp;
    double costab[100], area;

    /* Check boundary cases */
    if (len < p->window_len)
	return 0;
    if (start_pos < 0)
	start_pos = 0;
    start_pos += win_len2;
    if (start_pos + win_len2 >= len)
	return 0;

    /* Produce cosine matrix spread over +1 to 0 in winlen/2 steps */
    for (area = i = 0; i < win_len2; i++) {
	costab[i] = cos(i*(PI/2)/win_len2);
	area += costab[i];
    }

    /* When to bail out due to lack of signal strength */
    lowest_total = p->signal_val * highest / 100.0 * area;

    peakpos[npp = 0] = 0;
    do {
	double ratio;

	total_l = 0;
	total_r = 0;
	for (i = 0; i < win_len2; i++) {
	    total_l += peaks[start_pos-1-i] * costab[i];
	    total_r += peaks[start_pos+i] * costab[i];
	}

	ratio = total_r ? (double)total_l/total_r : 0.0;

	/*
	if (total_r)
	    printf("%d: %d %d %d %f\n", level, start_pos+1, total_l, total_r, ratio);
	*/

	/* Abort when signal drops too low */
	if (total_l < lowest_total) {
	    /*
	      printf("bailout %d < %d\n",
		   total_l, lowest_total);
	    */
	    break;
	}

	/* Mark current highest value */
	if (ratio > p->threshold && best_val < ratio) {
	    int j;
	    double avg_local = 0;
	    double avg_long = 0;

	    /*
	     * First check whether over a longer range the peaks increase
	     * again. If so it's probably a poor stretch (poly-X perhaps)
	     * and not a real 'stop'.
	     */
	    for (j = start_pos; j < start_pos+win_len2; j++)
		avg_local += peaks[j];

	    for (avg_long = avg_local; j < start_pos+200 && j < len; j++)
		avg_long += peaks[j];

	    if (avg_long/(j-start_pos) < 1.2*avg_local/win_len2 ||
		avg_long/(j-start_pos) < (p->signal_val/2) * highest) {
		best_val = ratio;
		peakpos[npp] = start_pos;
		peakval[npp] = ratio;
	    } else {
		printf("Skipped peak at %d due to later peak increases "
		       "%f -> %f, high=%d\n",
		       start_pos,
		       avg_local/win_len2,
		       avg_long/(j-start_pos),
		       highest);
	    }
	}

	/*
	 * If we're more than 'N' away from a peak without being higher than
	 * it and we've dropped below the threshold then mark the highest
	 * point as a peak and initial a new peak search.
	 */
	if (peakpos[npp] && start_pos-peakpos[npp] > 25 &&
	    ratio <= p->threshold * 0.7) {
	    if (++npp >= 10)
		break;
	    peakpos[npp] = 0;
	    best_val = 0;
	}
    } while (++start_pos <= (len - win_len2 - 1));
    if (peakpos[npp])
	npp++;

    if (level == 0) {
	for (i = 0; i < npp; i++) {
	    int newwin = 0.6 * p->window_len;
	    int oldwin = p->window_len;
	    int refinedpos;

	    /*
	     * Refine position by using a smaller window length and searching
	     * again around startpos.
	     */
	    p->window_len = newwin;
	    refinedpos = scan_right(name, p, peaks,
				    peakpos[i]-oldwin/2,
				    highest, peakpos[i]+oldwin/2, 1);
	    p->window_len = oldwin;
	    if (!refinedpos)
		refinedpos = peakpos[i];

	    printf("%s Peak %d at %d / %d height %f\n",
		   name, i, refinedpos+1, peakpos[i]+1, peakval[i]);

	}
    }

    return peakpos[0];
}

/*
 * Allocates and initialises a peak-height array from a Read.
 * The buffer returned should be freed using xfree().
 * Returns NULL on failure.
 */
int *create_peaks(Read *r) {
    int i;
    int *peaks;
    int tot = 0;

#define CAP_WIN 50

    if (!(peaks = (int *)xcalloc(r->NBases, sizeof(int))))
	return NULL;

    for (i = 0; i < r->NBases; i++) {
	int bpos = r->basePos[i];
	int peak = 0;

	if (peak < r->traceA[bpos])
	    peak = r->traceA[bpos];
	if (peak < r->traceC[bpos])
	    peak = r->traceC[bpos];
	if (peak < r->traceG[bpos])
	    peak = r->traceG[bpos];
	if (peak < r->traceT[bpos])
	    peak = r->traceT[bpos];

	/* Cap any abnormally high peaks so as not to bias drop search */
	if (i < CAP_WIN) {
	    tot += peak;
	} else {
	    if (peak*CAP_WIN > tot*2) {
		/*printf("Capping peak %d from %d to %d\n",
		       i, peak, tot*2/CAP_WIN);
		*/
		peak = tot*2/CAP_WIN;
	    }
	    tot = tot - peaks[i-CAP_WIN] + peak;
	}

	peaks[i] = peak;
    }

    return peaks;
}

/*
 * Searches for a likely stop point in trace 'r'.
 *
 * Returns first base position for success.
 *        -1 for failure.
 */
static int find_stops(Read *r, params *p) {
    int *peaks;
    int pos, maxheight;

    if (NULL == (peaks = create_peaks(r)))
	return -1;

    pos = find_highest_peak(p, peaks, r->NBases, &maxheight);
    if (p->verbose)
	printf("max height (avg over winlen) = %d\n", maxheight);
    pos = scan_right(r->trace_name, p, peaks, pos, maxheight, r->NBases, 0);

    xfree(peaks);

    return pos;
}

static void usage(void) {
    fprintf(stderr,
	    "Usage: stops [-v]\n"
	    "             [-w window_len(30)]\n"
	    "             [-s signal_strength(5)]\n"
	    "             [-t threshold(3.0)]\n"
	    "             file ...\n");
    fprintf(stderr,
	    "  ( Eg. stops -v -w 50 -t 2.5 *SCF )\n");
    exit(1);
}

int main(int argc, char **argv) {
    int c, i, ret = 0;
    params p;

    /* Defaults */
    p.verbose = 0;
    p.signal_val = 5;
    p.window_len = 50;
    p.threshold = 3.0;

    while ((c = getopt(argc, argv, "vw:s:t:")) != -1) {
	switch (c) {
	case 'v':
	    p.verbose = 1;
	    break;

	case 's':
	    p.signal_val = atoi(optarg);
	    break;

	case 'w':
	    p.window_len = atoi(optarg);
	    break;

	case 't':
	    p.threshold = atof(optarg);
	    break;

	default:
	    usage();
	}
    }

    if (optind == argc)
	usage();

    for (i = optind; i < argc; i++) {
	int ret_val;
	Read *r;

	if (p.verbose)
	    printf("Processing %s\n", argv[i]);

	if (!(r = read_reading(argv[i], TT_ANY))) {
	    fprintf(stderr, "Couldn't read '%s'\n", argv[i]);
	    continue;
	}

	ret_val = find_stops(r, &p);
	if (ret_val >= 0)
	    ret_val = 0;

	read_deallocate(r);

	ret |= ret_val;
    }

    return ret ? 1 : 0;
}
