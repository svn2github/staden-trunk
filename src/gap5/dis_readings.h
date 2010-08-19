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
int disassemble_readings(GapIO *io, int *rnums, int nreads, int move,
                         int remove_holes, int duplicate_tags);

#endif
