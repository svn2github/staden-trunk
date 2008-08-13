#ifndef _TG_TCL_H_
#define _TG_TCL_H_

#include <tcl.h>

char *io_obj_as_string(GapIO *io);
GapIO *io_from_obj(Tcl_Obj *obj);

#endif /* _TG_TCL_H_ */

