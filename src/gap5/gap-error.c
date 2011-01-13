/*
 * File: gap-error.c
 * Version:
 *
 * Orig. Author: Simon Dear
 *               MRC Laboratory of Molecular Biology
 *	         Hills Road
 *	         Cambridge CB2 2QH
 *	         United Kingdom
 *
 * Description:
 *
 * Created:
 * Updated: jkb 16/8/94 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <signal.h>

/* 11/1/99 johnt - defines for Windows support */
#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif

#include "misc.h"
#include "gap-error.h"
#include "g-error.h"
#include "array.h"
#include "bitmap.h"
#include "stack_dump.h"

char *GapErrorString(int err) {
    if (err >= GAPERR_BASE) {
	switch(err) {
	case GAPERR_NO_ERROR:		return "no error";
	case GAPERR_INVALID_TYPE:	return "invalid type";
	case GAPERR_NOT_FOUND:		return "does not exist";
	case GAPERR_TRUSTME:		return "you just can't!";
	default:			return "unknown error";
	}
    }

    if (err >= ARRAY_ERR_START)
	return ArrayErrorString(err);

    if (err >= BITMAP_ERR_START) {
	return BitmapErrorString(err);
    }

    return gerrors[err];
}

int gap_fatal_errors = 1;

/*************************************************************
 * Error reporting routines
 *************************************************************/

void error_sig(int sig) {
    verror(ERR_FATAL, "signal_handler",
	   "Program terminated unexpectedly with signal %d.", sig);
    if (sig != SIGINT
#if defined(SIGQUIT)
	&& sig != SIGQUIT
#endif
	) {
	verror(ERR_FATAL, "signal_handler",
	       "This is probably a bug.");
	verror(ERR_FATAL, "signal_handler",
	       "Please report all bug reports at "
	       "https://sourceforge.net/projects/staden/");
	/*
	 * Force crash.
	 */
/* 11/1/99 johnt - no SIGBUS under WINNT */
	signal(SIGSEGV, SIG_DFL);
#if defined(SIGBUS)
	signal(SIGBUS, SIG_DFL);
#endif
	stack_trace();
	abort();
    } else {
	exit(1);
    }
}

static void xperror_out_func(char *name, char *str) {
    verror(ERR_FATAL, name, str);
}

/* NOT FATAL */
void GAP_ERROR(char *reason, ...) {
    char buf[8192];
    va_list args;

    va_start(args, reason);
    vsprintf(buf, reason, args);

    xperror(buf, xperror_out_func);
}

/* FATAL */
void GAP_ERROR_FATAL(char *reason, ...) {
    char buf[8192];
    va_list args;

    va_start(args, reason);
    vsprintf(buf, reason, args);

    xperror(buf, xperror_out_func);

    if (gap_fatal_errors) {
#ifdef _WIN32
	/* 11/1/99 johnt - WINNT will not have stdout/err defined unless running in console mode
	 * so use a message box
	 */
	if( fileno(stderr) == -1 ){
	    MessageBox(NULL,buf,"Gap4 Error",MB_OK|MB_ICONERROR|MB_TASKMODAL);
	    return;
	}
#endif
	fputs("Gap4 has found an unrecoverable error - These are usually bugs.\nPlease submit all errors at https://sourceforge.net/projects/staden/\n", stderr);
	signal(SIGSEGV, SIG_DFL);
/* 11/1/99 johnt - No SIGBUS on WINNT */
#if defined(SIGBUS)
	signal(SIGBUS, SIG_DFL);
#endif
	stack_trace();
	*((int *)0) = 88;
	abort();
    } else {
	verror(ERR_FATAL, "NOTE  ", "Continue at own risk!");
    }
}

void set_gap_fatal_errors(int val) {
    gap_fatal_errors = val ? 1 : 0;
}

