#ifndef _TG_TRACKS_H_
#define _TG_TRACKS_H_

/*
 * Computed the read depth at a resolution of 1 base per value. This is
 * the starting point for all averaged lower-resolution plots.
 *
 * Returns a malloced track of read depth covering 'bin' on success.
 *         NULL on failure
 */
int *track_read_depth_r1(GapIO *io, bin_index_t *bin);

#endif /* _TG_TRACK_H_ */
