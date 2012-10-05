/*
 * Implements the various tracks (read depth, GC content, etc) that we
 * may wish to implement. The actual track_t object itself is in tg_track.c
 * (singular).
 */

#include "tg_gio.h"
#include "tg_tracks.h"
#include "xalloc.h"

/*
 * Computed the read depth at a resolution of 1 base per value. This is
 * the starting point for all averaged lower-resolution plots.
 *
 * Returns a malloced track of read depth covering 'bin' on success.
 *         NULL on failure
 */
int *track_read_depth_r1(GapIO *io, bin_index_t *bin) {
    int i, j;
    int *depth;
    rangec_t *r;
    tg_rec cnum;
    int nr, pos;
    contig_t *c;
    
    depth = (int *)xcalloc(bin->size, sizeof(*depth));

    if (!bin->rng)
	return depth;

    if (-1 == bin_get_position(io, bin, &cnum, &pos, NULL))
	return NULL;
    c = (contig_t *)cache_search(io, GT_Contig, cnum);
    if (NULL == c) return NULL;
    r = contig_seqs_in_range(io, &c, pos, pos+bin->size-1, 0, &nr);
    if (NULL == r) return NULL;

    /* Accumulate */
    for (i = 0; i < nr; i++) {
        /* FIXME: care about clipped depth? If so need to read l->rec */
        for (j = r[i].start; j <= r[i].end; j++) {
	    int p = j - pos;
	    if (p >= 0 && p < bin->size)
		depth[p]++;
        }
    }

    free(r);

    return depth;
}
