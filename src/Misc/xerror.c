/*
 * File: error.c
 *
 * Error handling. The basic method is when an error should be marked,
 * call xerr_set(errnum, string), which is a macro, and return -1, NULL, etc.
 * To print the error then call xperror("reason"). This'll print the
 * system error message (errno), the local error message (xerrno), and the
 * specified reason. Note that the system error message may not be valid,
 * but it's printed anyway.
 */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "xerror.h"
#include "misc.h"

static int   xerrnum;	/* Last error number */
static char *xerrstr;	/* Last error string */
static int   xerrline;	/* Error line number */
static char *xerrfile;	/* Error file */

/* 7/1/99 johnt - changed local variable errno to errnum throughout this file
                  as errno is a macro in Visual C++ !
*/

/*
 * Yet again on SunOS 4.1 an ANSI defined function is missing. We'll
 * write our own strerror using Sun's sys_errlist[] array (which doesn't
 * appear to be prototyped anywhere).
 */
#ifdef NOSTRERROR
char *strerror(int errnum) {
    extern char *sys_errlist[];
    extern int sys_nerr;

    if (errnum >= 0 && errnum < sys_nerr) {
	return sys_errlist[errnum];
    } else {
	return "unknown error";
    }
}
#endif

void xperror(char *reason, void (*out_func)(char *name, char *str)) {
    char buf[1024];
    sprintf(buf, "%s [%d]", strerror(errno), errno);
    out_func("SYSMSG ", buf);
    sprintf(buf, "%s [%d]", xerrstr, xerrnum); 
    out_func("ERROR  ", buf);
    sprintf(buf, "%s",      reason);
    out_func("COMMENT", buf);
    sprintf(buf, "%s:%d",   xerrfile, xerrline);
    out_func("FILE   ", buf);
}

int xerr_set_globals(int errnum, char *string, int line, char *file) {
    xerrnum  = errnum;
    xerrstr  = string;
    xerrline = line;
    xerrfile = file;

    return errnum;
}

int get_xerrnum(void) {
    return xerrnum;
}
