#ifndef _SEQED_GLOBALS_H_
#define _SEQED_GLOBALS_H_

#include <tcl.h>

extern Tcl_Interp *our_interp;
void init_global_interp(Tcl_Interp *interp);

#endif
