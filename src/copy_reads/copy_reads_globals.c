/*
 * Initialises global variables for use in copy_reads.
 *
 */
#include <tcl.h>

Tcl_Obj *copy_reads_defs = NULL;
static Tcl_Obj *defs_name2;

static char *copy_reads_defs_trace(ClientData cd, Tcl_Interp *interp,
				 char *n1, char *n2, int flags);

/* Main global setup function */
int init_copy_reads_globals(Tcl_Interp *interp) {

    /* do the same for copy_reads */
    {
	Tcl_Obj *val;
	defs_name2 = Tcl_NewStringObj("copy_reads_defs", -1); /* global */

	val = Tcl_ObjGetVar2(interp, defs_name2, NULL, TCL_GLOBAL_ONLY);
	if (NULL == val)
	    val = Tcl_NewStringObj("", -1);

	copy_reads_defs = Tcl_ObjSetVar2(interp, defs_name2, NULL, val,
				  TCL_GLOBAL_ONLY);
	Tcl_TraceVar(interp, "copy_reads_defs", TCL_TRACE_WRITES | TCL_GLOBAL_ONLY,
		     copy_reads_defs_trace, NULL);
    }
    return TCL_OK;
}

static char *copy_reads_defs_trace(ClientData cd, Tcl_Interp *interp,
				 char *n1, char *n2, int flags) {
    copy_reads_defs = Tcl_ObjGetVar2(interp, defs_name2, NULL, TCL_GLOBAL_ONLY);
    return NULL;
}
