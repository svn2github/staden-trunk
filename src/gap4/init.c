#include <tcl.h>

#include "tkEditor.h"
#include "tkEdNames.h"
#include "tkSheet.h"
#include "contigEditor.h"
#include "newgap_cmds.h"
#include "gap-tcl.h"
#include "active_tags.h"

/*
#ifndef _WIN32
#  define TRAP_SIGNALS
#endif
*/

#ifdef TRAP_SIGNALS
#include <signal.h>
#include "gap-error.h"
#endif

#ifdef TRAP_SIGNALS
void trap_signals(void) {
    /*
     * Use the BSD signal() command to trap probable program crashes.
     * This then adds a debug message to make sure that we tell people
     * to email us.
     */
#if defined(SIGBUS) /* 11/1/99 johnt - SIGBUS not defined on WINNT */
    signal(SIGBUS,  error_sig);
#endif
    signal(SIGSEGV, error_sig);
    signal(SIGILL,  error_sig);
    signal(SIGFPE,  error_sig);
    signal(SIGINT,  error_sig);
#if defined(SIGQUIT)
    signal(SIGQUIT, error_sig);
#endif
#if defined(SIGSYS)
    signal(SIGSYS,  error_sig);
#endif /* SIGSYS */
}
#endif /* TRAP_SIGNALS */

int Gap_Init(Tcl_Interp *interp) {
#ifdef TRAP_SIGNALS
    trap_signals();
#endif

    Editor_Init(interp);
    EdNames_Init(interp);
    Sheet_Init(interp);
    Ced_Init(interp);
    NewGap_Init(interp);
    Db_Init(interp);

    get_tag_types();

    return Tcl_PkgProvide(interp, "gap4", "1.0");
}

int Gap_SafeInit(Tcl_Interp *interp) {
    return Gap_Init(interp);
}

int Gap_Unload(Tcl_Interp *interp, int flags) {
    Tcl_SetResult(interp, "Pkg_Unload() function not implemented",
		  TCL_STATIC);
    return TCL_ERROR;
}

int Gap_SafeUnload(Tcl_Interp *interp, int flags) {
    Tcl_SetResult(interp, "Pkg_SafeUnload() function not implemented",
		  TCL_STATIC);
    return TCL_ERROR;
}
