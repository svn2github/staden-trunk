#include <assert.h>
#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include "user_defaults.h"
#include "licence.h"
#include "tclXkeylist.h"
#include "tclCanvGraph.h"
#include "tkCanvGraph.h"

extern int TclXKeylist_Init(Tcl_Interp *);
extern int Raster_Init(Tcl_Interp *);
extern int Tk_utils_Misc_Init(Tcl_Interp *);
extern int TextOutput_Init(Tcl_Interp *);
extern int Trace_Init(Tcl_Interp *);
extern int Sheet_Init(Tcl_Interp *);
extern int Container_Init(Tcl_Interp *);

static char *tk_utils_defs_trace(ClientData cd, Tcl_Interp *interp,
				 char *n1, char *n2, int flags);

Tcl_Obj *tk_utils_defs = NULL;    /* Todo: Add accessor function and make static */
static Tcl_Interp* our_interp;
static Tcl_Obj*    defs_name;


Tcl_Interp* GetInterp( void )
{
    assert(our_interp);
    return our_interp;
}

char* GetInterpResult( void )
{
    assert(our_interp);
    return Tcl_GetStringResult(GetInterp());
}

int Tk_utils_Init(Tcl_Interp *interp) {
    char *s, c[20], *lib = NULL, buf[1024];

    check_licence();
    c[0] = get_licence_type();
    c[1] = '\0';

    Tcl_SetVar2(interp, "licence","type", c, TCL_GLOBAL_ONLY);
    c[0] = get_licence_os();
    Tcl_SetVar2(interp, "licence","os", c, TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, "licence","hostid", get_licence_id(), TCL_GLOBAL_ONLY);
    sprintf(c, "%d", get_licence_users());
    Tcl_SetVar2(interp, "licence", "users", c, TCL_GLOBAL_ONLY);
    sprintf(c, "%d", get_licence_expire());
    Tcl_SetVar2(interp, "licence", "expire", c, TCL_GLOBAL_ONLY);

    our_interp = interp;

    /* Keyed lists from tclX */
    TclX_KeyedListInit(interp);
 
    /* Our updated Raster widget */
    Raster_Init(interp);

    /* Our own widgets and commands */
    Tk_utils_Misc_Init(interp);
    TextOutput_Init(interp);
    Trace_Init(interp);
    Sheet_Init(interp);

    Container_Init(interp);

    /* Tk_CreateItemType(&tkGraphType); */
    Tcl_GraphInit(interp);

    /* SeqReg_Init(interp); */

    /*
     * The auto_path.
     */
    if (NULL == (lib = getenv("STASH_LIB"))) {
	if (lib = getenv("STADENROOT")) {
	    sprintf(buf, "%s/lib/tk_utils", lib);
	    lib = buf;
	}
    }

    if (lib) {
	char *argv[3];
	int argc = 3;
	char *merged;
	argv[0] = "lappend";
	argv[1] = "auto_path";
	argv[2] = lib;
	Tcl_Eval(interp, merged = Tcl_Merge(argc, argv));
	Tcl_Free(merged);
    }

    /*
     * Set packages(name). This is done to prevent subsequent reloading
     * of this library (for efficiency reasons). The only reason that this
     * is necessary is that currently gap4 dynamically links with some
     * libraries at link time. When they're all at run time this won't
     * be necessary.
     */
    if (s = Tcl_GetVar2(interp, "packages", "tk_utils", TCL_GLOBAL_ONLY))
	sprintf(c, "%d", atoi(s)|2);
    else
	strcpy(c, "2");
    Tcl_SetVar2(interp, "packages", "tk_utils", c, TCL_GLOBAL_ONLY);

    /*
     * tk_utils_defs (a Tcl_Obj pointer)
     *
     * We keep this up to date by creating a write trace on the object and
     * doing an ObjGetVar2 when it changes. This way the object is always
     * valid.
     * Firstly we have to create tk_utils_defs though as initially it doesn't
     * exist.
     */
    {
	Tcl_Obj *val = Tcl_NewStringObj("", -1);

	defs_name = Tcl_NewStringObj("tk_utils_defs", -1); /* global */
	tk_utils_defs = Tcl_ObjSetVar2(interp, defs_name, NULL, val,
				       TCL_GLOBAL_ONLY);
	Tcl_TraceVar(interp, "tk_utils_defs",
		     TCL_TRACE_WRITES | TCL_GLOBAL_ONLY,
		     tk_utils_defs_trace, NULL);
    }

    return Tcl_PkgProvide(interp, "tk_utils", "1.0");
}

int Tk_utils_SafeInit(Tcl_Interp *interp) {
    return Tk_utils_Init(interp);
}

static char *tk_utils_defs_trace(ClientData cd, Tcl_Interp *interp,
				 char *n1, char *n2, int flags) {
    tk_utils_defs = Tcl_ObjGetVar2(interp, defs_name, NULL, TCL_GLOBAL_ONLY);
    return NULL;
}
