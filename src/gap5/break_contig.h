#ifndef _BREAK_CONTIG_H_
#define _BREAK_CONTIG_H_

#include "tg_gio.h"

/*
 * Breaks a contig in two such that snum is the right-most reading of
 * a new contig.
 */
int break_contig(GapIO *io, int crec, int cpos);

#endif
