#ifndef _BREAK_CONTIG_H_
#define _BREAK_CONTIG_H_

#include "tg_gio.h"

/*
 * Breaks a contig in two such that snum is the right-most reading of
 * a new contig.
 *
 * If break_holes is true then any gaps in coverage caused by breaking will
 * then attempt to break the contig further.
 *
 * Returns the new (right) contig record on success
 *        -1 on failure
 */
tg_rec break_contig(GapIO *io, tg_rec crec, int cpos, int break_holes);

/*
 * Compute the visible statr position of a contig. This isn't just the extents
 * of start_used / end_used in the bins as this can included invisible
 * data such as cached consensus sequences.
 *
 * 'from' is the starting point to search from in a contig iterator.
 * Specify CITER_CSTART if you don't know what this is, otherwise if the
 * contig has been edited (eg shrunk) then you may want to specify an older
 * start coord in order to correctly trim all data.
 */
int contig_visible_start(GapIO *io, tg_rec crec, int from);

/*
 * Compute the visible end position of a contig. This isn't just the extents
 * of start_used / end_used in the bins as this can included invisible
 * data such as cached consensus sequences.
 */
int contig_visible_end(GapIO *io, tg_rec crec, int to);

#endif
