/*
 * File: g-error.h
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: header file for gap server error servicing routines
 *
 * Created: prior to 28 September 1992
 */

#ifndef _G_ERROR_H_
#define _G_ERROR_H_

#include "xerror.h"

/*
 * About all we can say regarding the official ANSI errno.h is that it
 * defines EDOM and ERANGE. C implementations typically define many more
 * such as EPERM, ENOENT, etc. POSIX.1 defines many more, but not all
 * the ones we want (such as EWOULDBLOCK)
 */

/*
 * The G lib error messages. Numbers larger than (or equal to) GERR_START
 * are internal error messages; often being more general such as
 * 'write failed'. Other messages are synonyms for POSIX errors.
 *
 * To print an error both the gerrnum and errno are checked.
 */
/* General IO */
#define GERR_NONE			0
#define GERR_UNKNOWN			1
#define GERR_NOT_IMPLEMENTED		2
#define GERR_NAME_TOO_LONG		3 /* ENAMETOOLONG */
#define GERR_FILE_EXISTS		4 /* EEXIST */
#define GERR_CANT_CREATE		5
#define GERR_OPENING_FILE		6
#define GERR_CANT_LOCK			7
#define GERR_WOULD_BLOCK		8
#define GERR_PERMISSION			9  /* EACCESS */
#define GERR_OUT_OF_MEMORY		10 /* ENOMEM */
#define GERR_FILE_FULL			11
#define GERR_INVALID_ARGUMENTS		12
#define GERR_INVALID_RECORD_ID		13
#define GERR_READ_ERROR			14
#define GERR_WRITE_ERROR		15
#define GERR_SEEK_ERROR			16
#define GERR_MAX_CLIENTS		19
#define GERR_ALREADY_CONNECTED  	20
#define GERR_SYNC			21

/* Free Tree code */
#define TREE_NO_ERROR 			GERR_NONE
#define TREE_OUT_OF_MEMORY		GERR_OUT_OF_MEMORY
#define TREE_OUT_OF_SPACE		GERR_FILE_FULL
#define TREE_NOT_FOUND			22
#define TREE_OVERLAP			23

#define GERR_WRONG_BITSIZE              24

/*
 * A "g" interface to the xerr_set() code. All we need is an error number.
 * This is then translated into both a unique number and a string to pass
 * to the xerr_set routine.
 */
extern char *gerrors[];
/*int gerr_set(int errcode); */
#define gerr_set(e) (gerr_set_lf((e),__LINE__,__FILE__))
int gerr_set_lf(int errcode, int line, char *file);

#endif /*_G_ERROR_H_*/
