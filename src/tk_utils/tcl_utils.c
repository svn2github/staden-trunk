#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <tcl.h>

#include "tcl_utils.h"
#include "text_output.h"
#include "misc.h"
#include "vlen.h"
#include "tclXkeylist.h"

/*
 *-----------------------------------------------------------------------------
 * Result handling mechanisms
 *-----------------------------------------------------------------------------
 */

/*
 * We use sprintf to generate a result and then pass it onto Tcl_SetResult().
 *
 * You maybe wondering why we don't simply sprintf into interp->result. The
 * complete answer is long and tedious, but it involves requiring usage of
 * the Tcl_ResetResult function. The simplest answer is to avoid it
 * completely.
 */
__PRINTF_FORMAT__(2,3)
void vTcl_SetResult(Tcl_Interp *interp, char *fmt, ...) {
    static char string[8192];
    char *stringp = string;
    int len;
    va_list args;
    
    va_start(args, fmt);

    /* Use static buffer for small output */
    if ((len = vflen(fmt, args)) > 8192) {
	if (NULL == (stringp = (char *)xmalloc(len))) {
	    verror(ERR_FATAL, "vTcl_SetResult", "out of memory");
	    return;
	}
    }

    vsprintf(stringp, fmt, args);
    va_end(args);

    Tcl_SetResult(interp, stringp, TCL_VOLATILE);

    if (stringp != string)
	xfree(stringp);
}


__PRINTF_FORMAT__(2,3)
char *vTcl_DStringAppend(Tcl_DString *dsPtr, char *fmt, ...) {
    char string[8192], *stringp = string;
    va_list args;
    int len;
    char *ret;

    va_start(args, fmt);

    /* Use static buffer for small output */
    if ((len = vflen(fmt, args)) > 8192) {
	if (NULL == (stringp = (char *)xmalloc(len))) {
	    verror(ERR_FATAL, "vTcl_DStringAppend", "out of memory");
	    return NULL;
	}
    }

    vsprintf(stringp, fmt, args);
    ret = Tcl_DStringAppend(dsPtr, stringp, -1);
    va_end(args);

    if (stringp != string)
	xfree(stringp);

    return ret;
}

__PRINTF_FORMAT__(2,3)
char *vTcl_DStringAppendElement(Tcl_DString *dsPtr, char *fmt, ...) {
    char string[8192], *stringp = string;
    va_list args;
    int len;
    char *ret;

    va_start(args, fmt);

    /* Use static buffer for small output */
    if ((len = vflen(fmt, args)) > 8192) {
	if (NULL == (stringp = (char *)xmalloc(len))) {
	    verror(ERR_FATAL, "vTcl_DStringAppend", "out of memory");
	    return NULL;
	}
    }

    vsprintf(stringp, fmt, args);
    ret = Tcl_DStringAppendElement(dsPtr, stringp);
    va_end(args);

    if (stringp != string)
	xfree(stringp);

    return ret;
}

/*
 * Turns a non writable string to a writable string. This is required
 * for passing strings to many tcl and keyed list routines (although I
 * consider this a bug of tcl). This routine name is intentionally short
 * to keep the impact of this untidyness to a minimum.
 *
 * NB: Only handles up to 8K long strings, which is ok for the Tcl purposes.
 */
char *w(char *str) {
    static char buf[8192];

    strcpy(buf, str);
    return buf;
}

__PRINTF_FORMAT__(1,2)
char *vw(char *fmt, ...) {
    static char buf[8192];
    va_list args;

    va_start(args, fmt);

    vsprintf(buf, fmt, args);
    return buf;
}


/*
 *-----------------------------------------------------------------------------
 * Pointer to Tcl and Tcl to pointer conversions.
 *-----------------------------------------------------------------------------
 */

/*
 * These allow C to pass pointers into Tcl, and for Tcl to pass them back
 * again. Tcl itself cannot ever make use of the pointers directly, but it
 * provides an easy method of passing C data through Tcl back into C.
 * Note that the pointer must not be reallocated after being converted to Tcl.
 *
 * WARNING: If you change how these work, be warned that the contig_selector
 * assumes that they start with "p_".
 */

/*
 * C pointer to Tcl string.
 * The returned value is 'owned' by this function and is valid until the next
 * call.
 */
char *CPtr2Tcl(void *ptr) {
    static char buf[1024];
    sprintf(buf, "p_%p", ptr);
    return buf;
}

/*
 * Tcl string to C pointer
 */
void *TclPtr2C(char *str) {
    void *ptr;

/* SunOS4.x sscanf can't always cope with %p */
#ifndef BUGGY_SSCANF
    sscanf(str, "p_%p", &ptr);
#else
    sscanf(str, "p_%lx", &ptr);
#endif
    return ptr;
}

/*
 *-----------------------------------------------------------------------------
 * "Defaults" reading
 *-----------------------------------------------------------------------------
 */

char *get_default_string(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *name) {
    Tcl_Obj *val;
    char *str;

    TclX_KeyedListGet(interp, defs_ptr, name, &val);
    if (val) {
	str = Tcl_GetStringFromObj(val, NULL);
	return str;
    } else {
	fprintf(stderr, "Invalid key '%s'\n", name);
	return NULL;
    }
}

char *get_default_astring(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *name) {
    Tcl_Obj *val;
    char *str;

    TclX_KeyedListGet(interp, defs_ptr, name, &val);
    if (val) {
	str = Tcl_GetStringFromObj(val, NULL);
	return strdup(str);
    } else {
	fprintf(stderr, "Invalid key '%s'\n", name);
	return NULL;
    }
}

int get_default_int(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *name) {
    Tcl_Obj *val;
    int i;

    TclX_KeyedListGet(interp, defs_ptr, name, &val);
    if (val) {
	(void)Tcl_GetIntFromObj(interp, val, &i);
	return i;
    } else {
	fprintf(stderr, "Invalid key '%s'\n", name);
	return -1;
    }
}

double get_default_double(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *name) {
    Tcl_Obj *val;
    double i;

    TclX_KeyedListGet(interp, defs_ptr, name, &val);
    if (val) {
	(void)Tcl_GetDoubleFromObj(interp, val, &i);
	return i;
    } else {
	fprintf(stderr, "Invalid key '%s'\n", name);
	return -1.0;
    }
}

Tcl_Obj *get_default_obj(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *name) {
    Tcl_Obj *val;

    TclX_KeyedListGet(interp, defs_ptr, name, &val);
    return val;
}
