#ifndef _NIP_GLOBALS_H_
#define _NIP_GLOBALS_H_

#include <tcl.h>

extern Tcl_Obj *nip_defs;

/* Main global setup function */
int nip_init_globals(Tcl_Interp *interp);

#endif
