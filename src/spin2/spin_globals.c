#include <string.h>
#include "spin_globals.h"

Tcl_Obj *spin_defs = NULL;
static Tcl_Obj *defs_name;

int cutoff_align1 = 3;
int cutoff_align2 = 0;
int cutoff_align3 = -3;
char *symbol_align0 =NULL;
char *symbol_align1 =NULL;
char *symbol_align2 =NULL;
char *symbol_align3 =NULL;

static char *spin_defs_trace(ClientData cd, Tcl_Interp *interp,
				     char *n1, char *n2, int flags);

/* Main global setup function */
int spin_init_globals(Tcl_Interp *interp) {

    /*
     * spin_defs (a Tcl_Obj pointer)
     *
     * We keep this up to date by creating a write trace on the object and
     * doing an ObjGetVar2 when it changes. This way the object is always
     * valid.
     * Firstly we have to create nip_defs though as initially it doesn't
     * exist.
     */
    {
	Tcl_Obj *val;

	defs_name = Tcl_NewStringObj("spin_defs", -1); /* global */

	val = Tcl_ObjGetVar2(interp, defs_name, NULL, TCL_GLOBAL_ONLY);
	if (NULL == val)
	    val = Tcl_NewStringObj("", -1);

	spin_defs = Tcl_ObjSetVar2(interp, defs_name, NULL, val,
					   TCL_GLOBAL_ONLY);
	Tcl_TraceVar(interp, "spin_defs", 
		     TCL_TRACE_WRITES | TCL_GLOBAL_ONLY,

		     spin_defs_trace, NULL);

        symbol_align0 = Tcl_Alloc(2);
	strcpy(symbol_align0, "*");  
	symbol_align1 = Tcl_Alloc(2);
	strcpy(symbol_align1, "|");
        symbol_align2 = Tcl_Alloc(2);
	strcpy(symbol_align2, ":");
        symbol_align3 = Tcl_Alloc(2);
	strcpy(symbol_align3, ".");

Tcl_LinkVar(interp, "cutoff_align1", (char *)&cutoff_align1, TCL_LINK_INT);
Tcl_LinkVar(interp, "cutoff_align2", (char *)&cutoff_align2, TCL_LINK_INT);
Tcl_LinkVar(interp, "cutoff_align3", (char *)&cutoff_align3, TCL_LINK_INT);
Tcl_LinkVar(interp, "symbol_align0", (char *)&symbol_align0, TCL_LINK_STRING);
Tcl_LinkVar(interp, "symbol_align1", (char *)&symbol_align1, TCL_LINK_STRING);
Tcl_LinkVar(interp, "symbol_align2", (char *)&symbol_align2, TCL_LINK_STRING);
Tcl_LinkVar(interp, "symbol_align3", (char *)&symbol_align3, TCL_LINK_STRING);

 }
    
    return TCL_OK;

}

static char *spin_defs_trace(ClientData cd, Tcl_Interp *interp,
				     char *n1, char *n2, int flags) {
    spin_defs = Tcl_ObjGetVar2(interp, defs_name, NULL, 
				       TCL_GLOBAL_ONLY);
    return NULL;
}

