#ifndef _DIS_READINGS_H_
#define _DIS_READINGS_H_

/*
 * Calls remove_contig_holes for all contigs in the DB.
 * Returns the same codes.
 */
int remove_contig_holes_all(GapIO *io);

/*
 * remove_contig_holes - checks for gaps in a contig at the start, end or
 * internal. Internal holes require splitting the contig in two, while start
 * and end gaps simply require adjustment of the length and shuffling sequence
 * positions.
 *
 * This function obtains all information from the disk (or cached) database 
 * structures rather than the internal arrays incase these are not
 * consistent (yet). It does however assume that the reading positions
 * are sorted left to right.
 *
 * Returns 0 for success
 *        -1 for error
 */
int remove_contig_holes(GapIO *io, int cnum);

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
			 int duplicate_tags);

/**
 * Deletes an entire contig icluding all the readings, annotations, notes, etc
 * This is a wrapper around disassemble readings.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int delete_contig(GapIO *io, int contig);

#endif
