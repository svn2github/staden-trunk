#ifndef _CAPTURE_H_
#define _CAPTURE_H_

#include <tcl.h>

int tcl_capture(ClientData clientData, Tcl_Interp *interp,
		int argc, char **argv);

#endif
