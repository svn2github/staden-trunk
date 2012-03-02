#ifndef _TG_TCL_H_
#define _TG_TCL_H_

#include <tcl.h>

char *io_obj_as_string(GapIO *io);
GapIO *io_from_obj(Tcl_Obj *obj);

int G5_Init(Tcl_Interp *interp);
int G5_SafeInit(Tcl_Interp *interp);
int G5_Unload(Tcl_Interp *interp, int flags);
int G5_SafeUnload(Tcl_Interp *interp, int flags);

#endif /* _TG_TCL_H_ */

