#ifndef _SIP_GLOBALS_H_
#define _SIP_GLOBALS_H_

#include <tcl.h>

extern Tcl_Obj *sip_defs;		/* "defs" */

/* Main global setup function */
int sip_init_globals(Tcl_Interp *interp);

#endif
