#ifndef _IMPORT_GFF_H_
#define _IMPORT_GFF_H_

#include <tcl.h>

int tcl_import_gff(ClientData clientData, Tcl_Interp *interp,
		   int objc, Tcl_Obj *CONST objv[]);

#endif /* _IMPORT_GFF_H_ */
