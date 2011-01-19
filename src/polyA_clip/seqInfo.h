/*
 * Copyright (c) Medical Research Council 1998. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies and that credit is given
 * where due.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * Kathryn Beal, as part of the Staden Package at the MRC Laboratory of
 * Molecular Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

#ifndef _SEQINFO_H_
#define _SEQINFO_H_

#include "misc.h"

#include <io_lib/array.h>
#include <io_lib/expFileIO.h>

#define DB_NAMELEN 32

#define ERR_WARN 0
#define ERR_FATAL 1

/*
 * The following structure contains all the useful information
 * extracted from the experiment/staden file
 */
typedef struct {

    /* details from the staden/experiment file */
    Exp_info *e;		/* almost everything what we want is here */

    /* derived values */
    int length;
    int start;
    int end;

    /* other values */
    int1 *confidence;
    int2 *origpos;

} SeqInfo;

extern SeqInfo *allocSeqInfo(void);
extern void freeSeqInfo(SeqInfo *si);

extern SeqInfo *read_sequence_details(char *filename, int ignore_vec);

extern Exp_info *exp_read_staden_info(char *filename);
/*
 * Read a staden file into an Exp_info data structure
 */

char *read_sequence_name(SeqInfo *si);

int SeqInfo_conf(SeqInfo *si, int1 *conf, int length);
void SeqInfo_opos(SeqInfo *si, int2 *opos, int length);

#endif /*_SEQINFO_H_*/
