/* 
 * tkAppInit.c --
 *
 *	Provides a default version of the Tcl_AppInit procedure for
 *	use in wish and similar Tk-based applications.
 *
 * Copyright (c) 1993 The Regents of the University of California.
 * Copyright (c) 1994 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

/* static char sccsid[] = "@(#) tkAppInit.c 1.12 94/12/17 16:30:56"; */

#include <stdlib.h>
#include <string.h>
#include "tk.h"

#include "newgap_cmds.h"
#include "gap-tcl.h"

/*
 * The following variable is a special hack that is needed in order for
 * Sun shared libraries to be used for Tcl.
 */

#ifdef NEED_MATHERR
extern int matherr();
int *tclDummyMathPtr = (int *) matherr;
#endif

/*
 *----------------------------------------------------------------------
 *
 * main --
 *
 *	This is the main program for the application.
 *
 * Results:
 *	None: Tk_Main never returns here, so this procedure never
 *	returns either.
 *
 * Side effects:
 *	Whatever the application does.
 *
 *----------------------------------------------------------------------
 */

int
main(int argc, char **argv) {
    extern char *prog_name;

    prog_name = argv[0];

    Tcl_Main(argc, argv, Tcl_AppInit);
    return 0;			/* Needed only to prevent compiler warning. */
}

/*
 *----------------------------------------------------------------------
 *
 * Tcl_AppInit --
 *
 *	This procedure performs application-specific initialization.
 *	Most applications, especially those that incorporate additional
 *	packages, will have their own version of this procedure.
 *
 * Results:
 *	Returns a standard Tcl completion code, and leaves an error
 *	message in interp->result if an error occurs.
 *
 * Side effects:
 *	Depends on the startup script.
 *
 *----------------------------------------------------------------------
 */


int tkinit(ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {
    return Tk_Init(interp);
}

int
Tcl_AppInit(Tcl_Interp *interp)
{
   extern int Tk_utils_Init(Tcl_Interp *);
   extern int Gap_Init(Tcl_Interp *interp);
   char *lib, buf[1025];

    /*
     * Redefine TCL_LIBRARY and TK_LIBRARY to the staden package versions
     * in $STADLIB/{tcl,tk}.
     */
    if (NULL != (lib = getenv("STADLIB"))) {
        sprintf(buf, "TCL_LIBRARY=%s/tcl", lib);
        Tcl_PutEnv(buf);
        sprintf(buf, "TK_LIBRARY=%s/tk", lib);
        Tcl_PutEnv(buf);
    }

    if (Tcl_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
/*
    if (Tk_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
*/
    Tcl_CreateCommand(interp, "tkinit", tkinit,
		      (ClientData)NULL, NULL);

    /*
     * Call the init procedures for included packages.  Each call should
     * look like this:
     *
     * if (Mod_Init(interp) == TCL_ERROR) {
     *     return TCL_ERROR;
     * }
     *
     * where "Mod" is the name of the module.
     */

    /* Library init routines */
    if (Tk_utils_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }

    /* Gap4 init routines */
    if (Gap_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }

    /*
     * Call Tcl_CreateCommand for application-specific commands, if
     * they weren't already created by the init procedures called above.
     */

    {
        char *s, c[10];
	/*
	 * Set packages(name). This is done to prevent subsequent reloading
	 * of this library (for efficiency reasons). The only reason that this
	 * is necessary is that currently gap4 dynamically links with some
	 * libraries at link time. When they're all at run time this won't
	 * be necessary.
	 */
	if (NULL != (s = Tcl_GetVar2(interp, "packages", "gap", TCL_GLOBAL_ONLY)))
	    sprintf(c, "%d", atoi(s)|2);
	else
	    strcpy(c, "2");
	Tcl_SetVar2(interp, "packages", "gap", c, TCL_GLOBAL_ONLY);
    }

    /*
     * Specify a user-specific startup file to invoke if the application
     * is run interactively.  Typically the startup file is "~/.apprc"
     * where "app" is the name of the application.  If this line is deleted
     * then no user-specific startup file will be run under any conditions.
     */

    Tcl_SetVar(interp, "tcl_rcFileName", "~/.wishrc", TCL_GLOBAL_ONLY);
    return TCL_OK;
}
