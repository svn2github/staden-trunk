#include <xalloc.h>
#include "tg_gio.h"
#include "read_depth.h"

#define WIN_ITEMS 1024
#define WIN_SZ2 7 /* half window size, +/- this value */

/*
 * Computes the reading depth over a range within the contig. Note that 
 * the resolution of data returned may not be at the per-base level
 * depending on the size of the range queried. The actual resolution
 * of the data is samples from start_out to end_out in steps of
 * inc_out.
 *
 * The array is xmalloced() and it is up to the caller to xfree() it.
 *
 * Returns average sequence depth on success (plus first/res filled out)
 *         NULL on failure (first/res undefined).
 */
int *avg_sequence_depth(GapIO *io, tg_rec cnum, int start, int end,
			int *start_out, int *end_out, int *inc_out) {
    int len = end-start+1, len2;
    int rshift;
    int *depth;
    int i, j;
    rangec_t *r;
    int nr;
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, cnum);

    if (NULL == c) return NULL;

    /* Compute resolution such that length is less than WIN_ITEMS */
    rshift = 0;
    len2 = len;
    while (len2 > WIN_ITEMS) {
	len2>>=1;
	rshift++;
    }
    
    start = start & ~((1<<rshift)-1);
    end = (end & ~((1<<rshift)-1)) + 1;

    *start_out = start;
    *end_out = end;
    *inc_out = 1<<rshift;

    depth = xcalloc(len2+1, sizeof(*depth));
    if (NULL == depth) return NULL;

    r = contig_seqs_in_range(io, &c, start, end, 0, &nr);
    if (NULL == r) {
	free(depth);
	return NULL;
    }

    for (i = 0; i < nr; i++) {
	/* FIXME: care about clipped depth? */
	for (j = r[i].start; j <= r[i].end; j++) {
	    int p = j-start;
	    if (p >= 0 && p < len)
		depth[p >> rshift]++;
	}
    }

    for (i = 0; i < len2; i++) {
	depth[i] /= 1<<rshift;
    }

    free(r);
    return depth;
}

#if 0
/*
 * As per avg_sequence_depth but this returns a min_max_avg_t structure
 * containing minimum, maximum and average depth values.
 *
 * Returns average sequence depth on success
 *         NULL on failure
 */
min_max_avg_t *sequence_depth(GapIO *io, tg_rec cnum, int start, int end,
			      int *start_out, int *end_out, int *inc_out) {
    int len = end-start+1, len2;
    int rshift;
    int *depth;
    min_max_avg_t *mma_depth;
    int i, j;
    rangec_t *r;
    int nr;
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, cnum);

    /* Compute resolution such that length is less than WIN_ITEMS */
    rshift = 0;
    len2 = len;
    while (len2 > WIN_ITEMS) {
	len2>>=1;
	rshift++;
    }
    
    start = start & ~((1<<rshift)-1);
    end = (end & ~((1<<rshift)-1)) + 1;

    *start_out = start;
    *end_out = end;
    *inc_out = 1<<rshift;

    depth = xcalloc(len+1, sizeof(*depth));
    mma_depth = xcalloc(len2+1, sizeof(*mma_depth));

    r = contig_seqs_in_range(io, &c, start, end, 0, &nr);

    /* Accumulate */
    for (i = 0; i < nr; i++) {
	/* FIXME: care about clipped depth? */
	for (j = r[i].start; j <= r[i].end; j++) {
	    int p = j-start;
	    if (p >= 0 && p < len)
		depth[p]++;
	}
    }

    /* Average */
    for (i = 0; i < len; i += 1<<rshift) {
	uint64_t min, max, tot;

	min = INT_MAX; max = 0; tot = 0;
	for (j = 0; j < 1<<rshift; j++) {
	    if (min > depth[i+j])
		min = depth[i+j];
	    if (max < depth[i+j])
		max = depth[i+j];
	    tot += depth[i+j];
	}
	mma_depth[i >> rshift].avg = tot / (1<<rshift);
	mma_depth[i >> rshift].min = min;
	mma_depth[i >> rshift].max = max;
    }

    /*
     * Right idea I think, but should be moved to the viewer portion.
     */
    if (1) {
	int tot = 0;
	int *avg2 = xcalloc(len2+1, sizeof(*avg2));
	for (i = 0; i <= WIN_SZ2*2; i++)
	    tot += mma_depth[i].avg;
	for (i = WIN_SZ2; i < len2-WIN_SZ2; i++) {
	    avg2[i] = tot/(WIN_SZ2*2+1);
	    tot += mma_depth[i+WIN_SZ2+1].avg - mma_depth[i-WIN_SZ2].avg;
	}
	for (i = 0; i < len2; i++)
	    mma_depth[i].avg = avg2[i];
	xfree(avg2);
    }

    free(r);
    free(depth);

    return mma_depth;
}
#endif

/*
 * As per avg_sequence_depth but this returns a min_max_avg_t structure
 * containing minimum, maximum and average depth values.
 *
 * Returns average sequence depth on success
 *         NULL on failure
 */
min_max_avg_t *sequence_depth(GapIO *io, tg_rec cnum, int start, int end,
			      int *start_out, int *end_out, int *inc_out) {
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, cnum);
    int len = end-start+1, i;
    double bpv;
    min_max_avg_t *mma_depth;
    track_t *track;

    bpv = (double)len / WIN_ITEMS;
    mma_depth = xcalloc(WIN_ITEMS, sizeof(*mma_depth));

    track = contig_get_track(io, &c, start, end, TRACK_READ_DEPTH, bpv);
    for (i = 0; i < WIN_ITEMS; i++) {
	mma_depth[i].avg = arr(int, track->data, i);
	mma_depth[i].min = mma_depth[i].avg;
	mma_depth[i].max = mma_depth[i].avg;
    }

    *start_out = start;
    *end_out = end;
    *inc_out = bpv;

    track_free(track);

    return mma_depth;
}
