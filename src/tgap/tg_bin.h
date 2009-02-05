#ifndef _TG_BIN_H_
#define _TG_BIN_H_

/* #define MIN_BIN_SIZE 256 */
#define MIN_BIN_SIZE 4096

/* Size of bin to use for depth track */
#define RD_ELEMENTS 8192

int bin_new(GapIO *io, int pos, int sz, int parent, int parent_type);
bin_index_t *bin_add_range(GapIO *io, contig_t **c, range_t *r,
			   range_t **r_out);

/*
 * A bit like bin_get_track, but this is designed to auto-generate and
 * update the track as desired. The expectation is that this will always
 * succeed and anything else is a fatal error.
 */
track_t *bin_query_track(GapIO *io, bin_index_t *bin, int type);

/*
 * Creates a fake track struct, to be freed with track_free
 *
 * Returns track_t on success
 *         NULL on failure
 */
track_t *track_create_fake(int type, int size);

/*
 * Some tracks are "official" cached objects and are deallocated as part of
 * the tg_cache system. Others are temporary on-the-fly structs generated
 * for the purpose of temporary display. These we need to free.
 *
 * This function knows which is which and does the appropriate thing.
 */
void track_free(track_t *t);

/*
 * This finds a bin suitable for a given range. If such a bin doesn't
 * exist it can optionally create one if the 'extend' flag is true.
 *
 * When a bin is found the absolute offset of that bin is returned
 * in 'offset_r' (may be NULL).
 *
 * Returns the bin pointer on success
 *         NULL on failure
 */
bin_index_t *bin_for_range(GapIO *io, contig_t **c,
			   int start, int end, int extend,
			   int *offset_r);

/*
 * Finds the contig number and position of the start of a bin.
 * The right position is obviously this + bin->size (regardless of
 * whether it has been complemented).
 *
 * Returns 0 on success (plus *contig & *pos)
 *        -1 on failure
 */
int bin_get_position(GapIO *io, bin_index_t *bin, int *contig, int *pos);

#endif /* _TG_CONTIG_H_ */
