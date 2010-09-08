#ifndef _READ_DEPTH_H
#define _READ_DEPTH_H

typedef struct {
    int min, max, avg;
} min_max_avg_t;

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
			int *start_out, int *end_out, int *inc_out);

/*
 * As per avg_sequence_depth but this returns a min_max_avg_t structure
 * containing minimum, maximum and average depth values.
 *
 * Returns average sequence depth on success
 *         NULL on failure
 */
min_max_avg_t *sequence_depth(GapIO *io, tg_rec cnum, int start, int end,
			      int *start_out, int *end_out, int *inc_out);

#endif /* _READ_DEPTH_H */
