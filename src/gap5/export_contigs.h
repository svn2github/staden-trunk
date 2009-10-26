#ifndef _EXPORT_CONTIGS_H_
#define _EXPORT_CONTIGS_H_

#include <tcl.h>

int tcl_export_contigs(ClientData clientData, Tcl_Interp *interp,
		       int objc, Tcl_Obj *CONST objv[]);

int tcl_export_tags(ClientData clientData, Tcl_Interp *interp,
		    int objc, Tcl_Obj *CONST objv[]);

#endif /* _EXPORT_CONTIGS_H_ */
