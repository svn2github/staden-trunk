#include <string.h>
#include <tcl.h>
#include <stdlib.h>

#include "tkSeqed.h"
#include "tkSeqedNames.h"
#include "nip_cmds.h"
#include "sequtils_cmds.h"

int Nip_Init(Tcl_Interp *interp) {
    char *s, c[20];

    /*
     * Set packages(name). This is done to prevent subsequent reloading
     * of this library (for efficiency reasons). The only reason that this
     * is necessary is that currently gap4 dynamically links with some
     * libraries at link time. When they're all at run time this won't
     * be necessary.
     */
    if (s = Tcl_GetVar2(interp, "packages", "nip", TCL_GLOBAL_ONLY))
	sprintf(c, "%d", atoi(s)|2);
    else
	strcpy(c, "2");
    Tcl_SetVar2(interp, "packages", "nip", c, TCL_GLOBAL_ONLY);

    if (Seqed_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }

    if (SeqedNames_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }

    if (NipCmds_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }

    return TCL_OK;
}
