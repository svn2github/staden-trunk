#ifndef TCL_IO_LIB_H
#define TCL_IO_LIB_H

#include <io_lib/Read.h>
#include <tcl.h>

int tcl_read_seq_trace(ClientData clientData,
		       Tcl_Interp *interp,
		       int objc,
		       Tcl_Obj *CONST objv[]);


#endif /* TCL_IO_LIB_H */
