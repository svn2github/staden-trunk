#include <string.h>
#include <tcl.h>

#include "dna_utils.h"
#include "sequence_formats.h"
#include "misc.h"
#include "genetic_code.h"
#include "sip_results.h"

Tcl_Obj *sip_defs = NULL;
static Tcl_Obj *defs_name;

static char *sip_defs_trace(ClientData cd, Tcl_Interp *interp,
			    char *n1, char *n2, int flags);

/* Main global setup function */
int sip_init_globals(Tcl_Interp *interp) {
    /*
     * sip_defs (a Tcl_Obj pointer)
     *
     * We keep this up to date by creating a write trace on the object and
     * doing an ObjGetVar2 when it changes. This way the object is always
     * valid.
     * Firstly we have to create sip_defs though as initially it doesn't
     * exist.
     */
    {
	Tcl_Obj *val;

	defs_name = Tcl_NewStringObj("sip_defs", -1); /* global */

	/*	printf ("Sip'default name%s\n", defs_name);*/

	val = Tcl_ObjGetVar2(interp, defs_name, NULL, TCL_GLOBAL_ONLY);
	if (NULL == val)
	    val = Tcl_NewStringObj("", -1);

	sip_defs = Tcl_ObjSetVar2(interp, defs_name, NULL, val,
				     TCL_GLOBAL_ONLY);
	Tcl_TraceVar(interp, "sip_defs", TCL_TRACE_WRITES | TCL_GLOBAL_ONLY,
		     sip_defs_trace, NULL);
    }

    /* set up lookup tables in readpam.c */
    set_dna_lookup();
    set_protein_lookup();
    init_genetic_code();

    /* set up identity protein and dna matrices */
    set_matrix_identity(PROTEIN);
    set_matrix_identity(DNA);

    return TCL_OK;

}

static char *sip_defs_trace(ClientData cd, Tcl_Interp *interp,
			    char *n1, char *n2, int flags) {
    sip_defs = Tcl_ObjGetVar2(interp, defs_name, NULL, TCL_GLOBAL_ONLY);
    return NULL;
}

