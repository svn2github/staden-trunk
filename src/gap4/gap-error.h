/*
 * File:
 * Version:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */

#ifndef _GAP_ERROR_H_
#define _GAP_ERROR_H_

#include "misc.h"

extern void error_sig(int sig);

/* NOT FATAL */
extern void GAP_ERROR(char *reason, ...) __PRINTF_FORMAT__(1,2);


/* FATAL */
extern void GAP_ERROR_FATAL(char *reason, ...) __PRINTF_FORMAT__(1,2);

extern void set_gap_fatal_errors(int val);

extern char *GapErrorString(int err);

/* GAPERR values must not intersect with GERR values */
#define GAPERR_BASE		1000
#define NUM_GAPERRS		2

#define GAPERR_NO_ERROR		(GAPERR_BASE+0)
#define GAPERR_INVALID_TYPE	(GAPERR_BASE+1)
#define GAPERR_NOT_FOUND	(GAPERR_BASE+2)
#define GAPERR_TRUSTME		(GAPERR_BASE+3)

#define gaperr_set(e) ((e)?xerr_set((e), GapErrorString((e))):0)

#endif /*_GAP_ERROR_H_*/
