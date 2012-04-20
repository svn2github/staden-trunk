#ifndef _TG_BIN_H_
#define _TG_BIN_H_

/* Default, overridden in tg_index */
#define MIN_BIN_SIZE 4096

/* Size of bin to use for depth track */
#define RD_ELEMENTS 8192

tg_rec bin_new(GapIO *io, int pos, int sz, tg_rec parent, int parent_type);

/*
 * Adds a range to the contig.
 *
 * Delay_nseq controls whether we update bin_incr_nseq() immediately or
 * whether to delay until later.
 * A value of 0 will update nseq on each call.
 * A value of 1 will update nseq whenever we happen to write to a different
 * bin than before.
 * A value of -1 will update any pending nseqs and do no other action (so
 * c, r, r_out, etc are all ignored and NULL is returned).
 *
 * Returns the bin we added the range to on success
 *         NULL on failure
 */
bin_index_t *bin_add_range(GapIO *io, contig_t **c, range_t *r,
			   range_t **r_out, int *complemented,
			   int delay_nseq);

/*
 * As per bin_add_range() but add the range to a specific bin.
 * We use this for annotations to keep them in the same bin as their
 * associated sequence.
 *
 * Returns the bin we added the range to on success
 *         NULL on failure
 */
bin_index_t *bin_add_to_range(GapIO *io, contig_t **c, tg_rec brec, range_t *r,
			      range_t **r_out, int *complemented,
			      int delay_nseq);

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
int bin_get_position(GapIO *io, bin_index_t *bin, tg_rec *contig, int *pos,
		     int *comp);

/*
 * Removes a record referenced from a known bin
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_item_from_bin(GapIO *io, contig_t **c, bin_index_t **binp,
			     int type, tg_rec rec);

/*
 * Removes a record referenced from a bin. As above but we only know the
 * record and not which bin it's part of
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_item(GapIO *io, contig_t **c, int type, tg_rec rec);

/*
 * A faster alternative to bin_remove_item_from_bin(). Note this does not
 * ensure the contig is consistent afterwards or that the bin ranges are
 * valid. These checks need to be done afterwards by calling
 * bin_set_used_range() and possibly consensus_unclipped_range().
 *
 * If bin_idx is already known, pass it in. Otherwise pass in -1.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int fast_remove_item_from_bin(GapIO *io, contig_t **c, bin_index_t **binp,
			      int type, tg_rec rec, int bin_idx);

/*
 * Removes the refpos marker at a specific position in the contig.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_refpos(GapIO *io, tg_rec crec, int pos);

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
int bin_get_item_position(GapIO *io, int type, tg_rec rec,
			  tg_rec *contig, int *start, int *end, int *orient,
			  tg_rec *bin, range_t *r_out, void **i_out);

/*
 * Computes the bin orientation with respect to the contig.
 * Returns 1 for complemented
 *         0 for uncomplemented.
 */
int bin_get_orient(GapIO *io, tg_rec rec);

/*
 * Adds 'n' to the nseq counter for a bin and all parent bins chaining up
 * to the root node. 'n' may be negative too.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_incr_nseq(GapIO *io, bin_index_t *bin, int n);
int bin_incr_nrefpos(GapIO *io, bin_index_t *bin, int n);
int bin_incr_nanno(GapIO *io, bin_index_t *bin, int n);


/*
 * Clears the BIN_CONS_VALID flag on any bins that may contain
 * consensus data over the region contig:start..end.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_invalidate_consensus(GapIO *io, tg_rec contig, int start, int end);

/*
 * Returns whether a bin is empty (has no ranges). This isn't as trivial
 * as it sounds as we may have bin->rng allocated and containing data, but
 * with all data elements being UNUSED.
 *
 * Returns 1 for empty
 *         0 if not
 */
int bin_empty(bin_index_t *bin);

/*
 * Call after updating range array to ensure that the bin start_used and
 * end_used values are correct.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_set_used_range(GapIO *io, bin_index_t *bin);

#endif /* _TG_CONTIG_H_ */
