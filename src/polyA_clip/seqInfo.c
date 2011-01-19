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

/*
 * This file contains the various utility functions required by exp_file.c
 * It primarily acts as an interface between exp_file.c and the expFileIO.c
 * and scf_extras.c files.
 * It is derived from the Gap4 source copy of seqInfo.c, but is much
 * reduced in size.
 */

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <stdlib.h>

#include "seqInfo.h"

#include <os.h>
#include <io_lib/expFileIO.h>
#include <io_lib/scf_extras.h>

#include "xalloc.h"

/*
 * Usage: verror(priority, name, format, args...);
 * NB: don't pass more than 8K per call
 */
/* ARGSUSED */
__PRINTF_FORMAT__(3,4)
void verror(int priority, const char *name, const char *fmt, ...) { 
    va_list args;

    fprintf(stderr, "%s: ", name);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
}

/*************************************************************
 * Create and free routines
 *************************************************************/
SeqInfo *allocSeqInfo(void)
/*
 *
 */
{
    SeqInfo *si;
    if ( (si = (SeqInfo *)xmalloc(sizeof(SeqInfo))) != NULL) {
	si->confidence = NULL;
	si->origpos = NULL;
	si->length = -1;
	si->start = -1;
	si->end = -1;
    }

    return si;
}



void freeSeqInfo(SeqInfo *si)
/*
 *
 */
{
#define FREE(F) \
    if (si->F != NULL) {\
        xfree(si->F);\
	si->F = NULL;\
    }

    if (si != NULL) {
	if(si->e != NULL) { exp_destroy_info(si->e); si->e = NULL; }
	FREE(confidence);
	FREE(origpos);
	xfree(si);
    }
}

/*************************************************************
 * Utilities
 *************************************************************/

static void determine_active_region(Exp_info *e, int *left, int *right,
				    int ignore_vec)
{
    /*
     * need to determine active region
     */
    int SQlen;
    int CSfrom, CSto, SL, SR, QL, QR;
    int l,r;
    int noCS; /* if there is no CS line */

    SQlen = strlen(exp_get_entry(e,EFLT_SQ));
    noCS = exp_get_rng(e,EFLT_CS,&CSfrom,&CSto);
    if (exp_get_int(e,EFLT_SL,&SL)) SL = 0;
    if (exp_get_int(e,EFLT_SR,&SR)) SR = SQlen+1;
    if (exp_get_int(e,EFLT_QL,&QL)) QL = 0;
    if (exp_get_int(e,EFLT_QR,&QR)) QR = SQlen+1;

    if (ignore_vec) {
	*left = QL;
	*right = QR;
	return;
    }

    /* find where the good bits are */
    if (SL>QL) l = SL; else l = QL;
    if (SR<QR) r = SR; else r = QR;

    /*
     * NOTE: tmp fix to remove CS clipping.
     */
    noCS = -1;
    if (!noCS) {
	/*
	 * Six posibilities:
	 *           l                r
	 * 1  -----  >                <
	 * 2  -------|------->        <
	 * 3  -------|<---------------|-------
	 * 4         >    <--------   |
	 * 5         >    <-----------|-------
	 * 6         >                <  -----
	 *
	 * key: [-] CS [>] new left c.o [<] new right c.o [|] orig c.o
	 *
	 * Coded for understanding not efficiency!!!
	 */
	if (/*1*/ CSfrom <= l && CSto <= l)
	    l = l; /* null statement */
	else if (/*2*/ CSfrom <= l+1 && CSto < r)
	    l = CSto;
	else if (/*3*/ CSfrom <= l+1 && CSto >= r)
	    r = l+1;
	else if (/*4*/ CSfrom < r && CSto < r)
	    r = CSfrom;
	else if (/*5*/ CSfrom < r && CSto >= r)
	    r = CSfrom;
	else if (/*6*/ CSfrom >= r && CSto >= r)
	    l = l; /* null statement */
    }
    if (l>r) l = r-1; /* just a quick consistency check */

    *left = l;
    *right = r;
}


/*
 * Reads in sequence from file filname.
 * Format of this file must be either:
 * 	staden (80 char lines, title lines prefixed with ";")
 * or   experiment file
 *
 * "ignore_vec" specifies whether to ignore the vector sequence lines (SL,
 * SR, CL, CL, CS) when determining the active region. This is necessary
 * when inputting extracted data to preassembly.
 *
 * returns length read in
 */
SeqInfo *read_sequence_details(char *filename, int ignore_vec)
{
    SeqInfo *si = NULL;
    Exp_info *e;

    /*
     * read sequence details into experiment file format
     */
    e = exp_read_info(filename);

    if ( e != NULL ) {
	exp_close(e);
	if (e->Nentries[EFLT_SQ] == 0 || ( si = allocSeqInfo() ) == NULL ) {
	    exp_destroy_info(e);
	    return NULL;
	} else {
	    si->e = e;
	    si->length = strlen(exp_get_entry(e,EFLT_SQ));
	    determine_active_region(e, &si->start, &si->end, ignore_vec);
	}

	/* orig pos and conf. values */
	if (e->Nentries[EFLT_ON]) {
	    int2 *opos;

	    if (NULL != (opos = (int2 *)xmalloc((si->length+1)*
						sizeof(int2)))) {
		if (str2opos(opos, si->length+1, exp_get_entry(e, EFLT_ON))
		    != si->length)
		    verror(ERR_WARN, "read_sequence_details",
			   "Experiment file %s - 'ON' line has wrong number "
			   "of items", filename);
		si->origpos = opos;
	    } else
		si->origpos = NULL;
	}

	if (e->Nentries[EFLT_AV]) {
	    int1 *conf;

	    if (NULL != (conf = (int1 *)xmalloc((si->length+1)*
						sizeof(int1)))) {
		if (str2conf(conf, si->length+1, exp_get_entry(e, EFLT_AV))
		    != si->length)
		    verror(ERR_WARN, "read_sequence_details",
			   "Experiment file %s - 'AV' line has wrong number "
			   "of items", filename);
		si->confidence = conf;
	    } else
		si->confidence = NULL;
	}
    }

    return si;
}


/*
 * Reads an ID line from an experiment file, which has already been loaded
 * into a SeqInfo structure.
 *
 * Returns name for success
 *         NULL for failure
 *
 * The name returned is valid only until the next call to idline().
 */
char *read_sequence_name(SeqInfo *si) {
    static char name[DB_NAMELEN+1];
    char *namep;
    int i;

    if (exp_Nentries(si->e, EFLT_ID) < 1) {
	verror(ERR_WARN, "read_sequence_name", "No ID line in experiment file");
	if (exp_Nentries(si->e, EFLT_EN) < 1) {
	    verror(ERR_WARN, "read_sequence_name", "Not even an EN line!");
	    return NULL;
	} else {
	    namep = exp_get_entry(si->e, EFLT_EN);
	}
    } else {
	namep = exp_get_entry(si->e, EFLT_ID);
    }

    /* Copy up to the first white space */
    i = 0;
    do {
	name[i++] = *namep++;
    } while (i < DB_NAMELEN &&
	     *namep != ' ' && *namep != '\t' &&
	     *namep != '\n' && *namep != '\r' &&
	     *namep != '\0');
    name[i] = 0;

    return name;
}


/*
 * Fills a buffer 'conf' of length 'length' with the confidence values.
 *
 * These are either read from si, read from the trace file
 * If no confidence values available they are guessed and we return -1.
 * Otherwise we return 0
 */
int SeqInfo_conf(SeqInfo *si, int1 *conf, int length) {
    if (si->confidence) {
	memcpy(conf, si->confidence, sizeof(int1) * length);
    } else {
	int guess = 0;

	/* Read from SCF file */
	if (exp_Nentries(si->e,EFLT_LT) && exp_Nentries(si->e,EFLT_LN)) {
	    if (0 != get_read_conf(si->e, length, NULL, conf))
		guess = 1;
	} else {
	    guess = 1;
	}

	if (guess) {
	    int i;

	    for (i = 0; i < length; i++)
		conf[i] = 2;

	    return -1;
	}
    }
    return 0;
}

/*
 * Fills a buffer 'opos' of length 'length' with the original positions.
 *
 * These are either read from si, or guessed.
 */
void SeqInfo_opos(SeqInfo *si, int2 *opos, int length) {
    if (si->origpos) {
        memcpy(opos, si->origpos, sizeof(int2) * length);
    } else {
	int i, j;
	char *seq;

	seq = exp_get_entry(si->e, EFLT_SQ);

	for (i = j = 0; i < length; i++)
	    if (seq[i] != '*')
		opos[i] = ++j;
	    else
		opos[i] = 0;
    }
}

