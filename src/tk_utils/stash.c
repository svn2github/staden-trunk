/*
 * A replacement for the Tcl_AppInit function that has delayed initialisation
 * of Tk_Init and adds on Tk_utils_Init.
 */

#include <tcl.h>
#include <tk.h>
#include "tkCanvGraph.h"

extern int Tk_Init(Tcl_Interp *interp);
extern int Tk_utils_Init(Tcl_Interp *interp);

static int tkinit(ClientData clientData,
		  Tcl_Interp *interp,
		  int argc,
		  char **argv) {
  int ret = Tk_Init(interp);
  Tk_CreateItemType(&tkGraphType);
  return ret;
}

int Stash_AppInit(Tcl_Interp *interp) {
    if (Tcl_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }

    /* Delay initialising Tk package */
    Tcl_CreateCommand(interp, "tkinit", tkinit,
		      (ClientData)NULL, NULL);

    /* Tk_utils initialisation */
    if (Tk_utils_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }

    Tcl_SetVar(interp, "tcl_rcFileName", "~/.stashrc", TCL_GLOBAL_ONLY);
    return TCL_OK;
}

