#ifndef _BREAK_CONTIG_H_
#define _BREAK_CONTIG_H_

#include "tg_gio.h"

/*
 * Breaks a contig in two such that snum is the right-most reading of
 * a new contig.
 *
 * If break_holes is true then any gaps in coverage caused by breaking will
 * then attempt to break the contig further.
 */
int break_contig(GapIO *io, tg_rec crec, int cpos, int break_holes);

#endif
