#ifndef _DIS_READINGS_H_
#define _DIS_READINGS_H_

#include <tg_gio.h>

/**
 * Removes a set of readings from either the contig or the database.
 *
 * When removing from the database we need to delete everything related to
 * that reading (annotations, template if not used elsewhere, etc). When
 * moving to a new contig, we prefer to keep all readings clustered together.
 *
 * Ie if we remove A & B (overlapping) from one contig and C from another
 * then we create two new contigs containing A & B in one and C in the other.
 *
 * When creating new contigs, we have the option of copying over any
 * overlapping consensus tags to the new contigs. This choice only refers to
 * consensus tags; reading tags are always copied.
 */
int disassemble_readings(GapIO *io, tg_rec *rnums, int nreads, int move,
                         int remove_holes, int duplicate_tags);

/*
 * As per disassemble readings, but removes entire contigs.
 *
 * This is substantially faster as it doesn't need to track a lot of the
 * changes to contig dimensions.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int disassemble_contigs(GapIO *io, tg_rec *cnums, int ncontigs);

void bin_destroy_recurse(GapIO *io, tg_rec rec);

/*
 * Looks for contig gaps between start..end in contig and if it finds them,
 * breaking the contig in two.
 */
int remove_contig_holes(GapIO *io, tg_rec contig, int start, int end,
			int empty_contigs_only);

#endif
