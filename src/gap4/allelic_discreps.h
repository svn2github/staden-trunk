#ifndef _ALLELIC_DISCREPS_H
#define _ALLELIC_DISCREPS_H

#include "dstring.h"
/*
 * MAIN ENTRY POINT.
 *
 * Given a contig containing a mixed assembly from two alleles, this
 * function attempts to split templates into two sets corresponding
 * to each of the two alleles.
 */
dstring_t *allelic_discreps(GapIO *io, int contig);


#endif
