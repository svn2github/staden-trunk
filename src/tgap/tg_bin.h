#ifndef _TG_BIN_H_
#define _TG_BIN_H_

/* Default, overridden in tg_index */
#define MIN_BIN_SIZE 4096

/* Size of bin to use for depth track */
#define RD_ELEMENTS 8192

int bin_new(GapIO *io, int pos, int sz, int parent, int parent_type);
bin_index_t *bin_add_range(GapIO *io, contig_t **c, range_t *r,
			   range_t **r_out, int *complemented);

/*
 * A bit like bin_get_track, but this is designed to auto-generate and
 * update the track as desired. The expectation is that this will always
 * succeed and anything else is a fatal error.
 */
track_t *bin_query_track(GapIO *io, bin_index_t *bin, int type);

/*
 * Finds the track of a given type for this bin.
 *
 * Returns track_t pointer on success (do not free)
 *         NULL on failure (eg bin too small)
 */
track_t *bin_get_track(GapIO *io, bin_index_t *bin, int type);

/*
 * Invalidates a track.
 * Returns 0 on success
 *        -1 on failure.
 */
int bin_invalidate_track(GapIO *io, bin_index_t *bin, int type);

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
 * The relative orientation of this bin to the contig may be returned in
 * complemented, if non NULL. Note this is not the same as the
 * BIN_COMPLEMENTED bit in bin->flags, as that is in indication of whether
 * this bin is complemented relative to its parent. Instead complemented is
 * an XORed product of all bins to this point.
 *
 * Returns the bin pointer on success
 *         NULL on failure
 */
bin_index_t *bin_for_range(GapIO *io, contig_t **c,
			   int start, int end, int extend,
			   int *offset_r, int *complemented);

/*
 * Finds the contig number and position of the start of a bin.
 * The right position is obviously this + bin->size (regardless of
 * whether it has been complemented).
 *
 * Returns 0 on success (plus *contig & *pos)
 *        -1 on failure
 */
int bin_get_position(GapIO *io, bin_index_t *bin, int *contig, int *pos);

/*
 * Removes a record referenced from a known bin
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_item_from_bin(GapIO *io, contig_t **c, bin_index_t **binp,
			     int type, int rec);

/*
 * Removes a record referenced from a bin. As above but we only know the
 * record and not which bin it's part of
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_item(GapIO *io, contig_t **c, int type, int rec);

/*
 * Finds the contig number and position of a record number.
 *
 * If non-NULL r_out is filled with the associated range_t struct.
 *
 * If non-NULL i_out is filled with a pointer to the object referred to
 * by type/rec. (This is just for minor optimisations.) If returned it
 * will have had cache_incr() run on it, so the caller should
 * use cache_decr() to permit deallocation.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_get_item_position(GapIO *io, int type, GRec rec,
			  int *contig, int *start, int *end, int *orient,
			  int *bin, range_t *r_out, void **i_out);

#endif /* _TG_CONTIG_H_ */
