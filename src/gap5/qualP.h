#ifndef _QUALP_H
#define _QUALP_H

#include "qual.h"

/*
 * This should only be included by qual.c
 */

#ifdef DEBUG_QUAL
static char *reasons[] = {
    "No data",
    "OK",
    "Poor data"
    };

static char *reasons2[] = {
    "Good == Good",
    "Good != Good",
    "Good == Bad ",
    "Good != Bad ",
    "Bad  == Good",
    "Bad  != Good",
    "Bad  == Bad ",
    "Bad  != Bad ",
    "Good    None",
    "Bad     None",
    "None    Good",
    "None    Bad ",
    "None    None"
    };
#endif

#define Q_NO_DATA	0	/* Nothing to declare */
#define Q_OK		1	/* Good data */
#define Q_POOR_DATA	2	/* Quality of base is <= qual_cutoff,  */
                                /* or data doesn't meet conflict %age tests */

/* Index with Quality (NODATA, OK, POOR), equal (2 consensuses) */
extern char q_lookup[3][3][2]; /* defined in qual.c */

#endif
