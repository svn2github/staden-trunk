/*
 * File: seqInfo.h
 * Version:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: misc routines in a temp file
 *
 * Created:
 * Updated:
 *
 */

#ifndef _SEQINFO_H_
#define _SEQINFO_H_

#include "os.h"
#include "array.h"

#include <io_lib/expFileIO.h>

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

/*
 * Read a staden file into an Exp_info data structure
 */

char *read_sequence_name(SeqInfo *si);

void SeqInfo_conf(SeqInfo *si, int1 *conf, int length);
void SeqInfo_opos(SeqInfo *si, int2 *opos, int length);
char *SeqInfo_used_seq(SeqInfo *si, int *length);

#endif /*_SEQINFO_H_*/
