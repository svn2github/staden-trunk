#include <string.h>

#include "sequence_formats.h"  /* DNA PROTEIN */
#include "dna_utils.h"
#include "nip_globals.h"
#include "genetic_code.h"

/* char *nip_defs = NULL; */
Tcl_Obj *nip_defs = NULL;
static Tcl_Obj *defs_name;

static char *nip_defs_trace(ClientData cd, Tcl_Interp *interp,
			    char *n1, char *n2, int flags);

/* Main global setup function */
int nip_init_globals(Tcl_Interp *interp) {

    set_dna_lookup();
    set_char_set(DNA);
    set_iubc_lookup();
    init_genetic_code();
    
    /*
     * nip_defs (a Tcl_Obj pointer)
     *
     * We keep this up to date by creating a write trace on the object and
     * doing an ObjGetVar2 when it changes. This way the object is always
     * valid.
     * Firstly we have to create nip_defs though as initially it doesn't
     * exist.
     */
    {
	Tcl_Obj *val;

	defs_name = Tcl_NewStringObj("nip_defs", -1); /* global */
     

	val = Tcl_ObjGetVar2(interp, defs_name, NULL, TCL_GLOBAL_ONLY);
	if (NULL == val)
	    val = Tcl_NewStringObj("", -1);

	nip_defs = Tcl_ObjSetVar2(interp, defs_name, NULL, val,
				     TCL_GLOBAL_ONLY);
	Tcl_TraceVar(interp, "nip_defs", TCL_TRACE_WRITES | TCL_GLOBAL_ONLY,
		     nip_defs_trace, NULL);
    }
    
    return TCL_OK;

}

static char *nip_defs_trace(ClientData cd, Tcl_Interp *interp,
			    char *n1, char *n2, int flags) {
    nip_defs = Tcl_ObjGetVar2(interp, defs_name, NULL, TCL_GLOBAL_ONLY);
    return NULL;
}
