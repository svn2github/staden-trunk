#include <tcl.h>
#include <stdio.h>
#include <stdlib.h>

int
Tcl_AppInit(Tcl_Interp *interp)
{
   char *lib, buf[1025];

    /*
     * Redefine TCL_LIBRARY and TK_LIBRARY to the staden package versions
     * in $STADLIB/{tcl,tk}.
     */
#ifndef _WIN32 /* 11/1/98 johnt - not required for WIN32 */
    if (NULL != (lib = getenv("STADLIB"))) {
        sprintf(buf, "TCL_LIBRARY=%s/tcl", lib);
        Tcl_PutEnv(buf);
        sprintf(buf, "TK_LIBRARY=%s/tk", lib);
        Tcl_PutEnv(buf);
    }
#endif

#ifdef TRAP_SIGNALS
    /*
     * Use the BSD signal() command to trap probable program crashes.
     * This then adds a debug message to make sure that we tell people
     * to email us.
     */
#if defined(SIGBUS) /* 11/1/99 johnt - SIGBUS not defined under WINNT */
    signal(SIGBUS,  error_sig);
#endif
    signal(SIGSEGV, error_sig);
    signal(SIGILL,  error_sig);
    signal(SIGFPE,  error_sig);
#if defined(SIGSYS)
    signal(SIGSYS,  error_sig);
#endif /* SIGSYS */
#endif /* TRAP_SIGNALS */

    if (Tcl_Init(interp) == TCL_ERROR) {
	WishPanic(interp->result);
	return TCL_ERROR;
    }

#ifdef _WIN32

    if( needConsole ){
	/* running in Windows mode, so initialise TK, and Init console */
	if (Tk_Init(interp) == TCL_ERROR) {
	    WishPanic(interp->result);
	    return TCL_ERROR;
	}
	Tcl_StaticPackage(interp, "Tk", Tk_Init, Tk_SafeInit);

        if (TkConsoleInit(interp) == TCL_ERROR) {
	    WishPanic(interp->result);
	    return TCL_ERROR;
	}
    }

    /* modified Windows Save Dialog, to set save type, based on user selection and file extension */
    Tcl_CreateCommand(interp,"sp_getSaveFile",Tk_CustomGetSaveFileCmd,(ClientData)NULL,NULL);
    /* modified Windows Open Dialog, to allow multiple file opens */
    Tcl_CreateCommand(interp,"sp_getOpenFile",Tk_CustomGetOpenFileCmd,(ClientData)NULL,NULL);

#endif


    /*
    Tcl_CreateCommand(interp, "tkinit", tkinit,
		      (ClientData)NULL, NULL);
    */

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
	WishPanic(interp->result);
	return TCL_ERROR;
    }

    /*
     * Call Tcl_CreateCommand for application-specific commands, if
     * they weren't already created by the init procedures called above.
     */

    /*
     * Specify a user-specific startup file to invoke if the application
     * is run interactively.  Typically the startup file is "~/.apprc"
     * where "app" is the name of the application.  If this line is deleted
     * then no user-specific startup file will be run under any conditions.
     */

    Tcl_SetVar(interp, "tcl_rcFileName", "~/.stashrc", TCL_GLOBAL_ONLY);
    return TCL_OK;
}
