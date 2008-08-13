/*
 * File: g-error.c
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: routines for servicing gap server errors
 *
 * Created: prior to 18 September 1992
 */

#include <stdio.h>

#include "g-error.h"
#include "g-misc.h" /* IMPORT: G_Number */

/* error list */
char *gerrors[] = {
    "", /* NONE */
    "unknown error",
    "operation not implemented",
    "file name too long",
    "file exists",
    "cannot create file",
    "cannot open file",
    "cannot set or upgrade lock",
    "operation would block",
    "operation not allowed",
    "out of memory",
    "file is full",
    "invalid arguments",
    "invalid record identifier",
    "read error",
    "write error",
    "seek error",
    "fatal corruption detected",
    "something rather pathetic has happened",
    "maximum number of clients connected",
    "client already connected",
    "couldn't sync data",

    /* Tree code */
    "item not found in free-tree",
    "overlap within free-tree",

    "Incompatible database bit-size"
};

#ifdef gerr_set
#    undef gerr_set
#endif

int gerr_set(int errcode) {
    if (errcode)
	return xerr_set(errcode, gerrors[errcode]);
    else
	return 0;
}

int gerr_set_lf(int errcode, int line, char *file) {
    if (errcode)
	return xerr_set_globals(errcode, gerrors[errcode], line, file);
    else
	return 0;
}

