#ifndef _TCLXKEYLIST_H_
#define _TCLXKEYLIST_H_

#include <tcl.h>

/*
 * Exported keyed list object manipulation functions.
 */
extern Tcl_Obj *
TclX_NewKeyedListObj _ANSI_ARGS_((void));

extern int
TclX_KeyedListGet _ANSI_ARGS_((Tcl_Interp *interp,
                               Tcl_Obj    *keylPtr,
                               char       *key,
                               Tcl_Obj   **valuePtrPtr));

extern int
TclX_KeyedListSet _ANSI_ARGS_((Tcl_Interp *interp,
                               Tcl_Obj    *keylPtr,
                               char       *key,
                               Tcl_Obj    *valuePtr));

extern int
TclX_KeyedListDelete _ANSI_ARGS_((Tcl_Interp *interp,
                                  Tcl_Obj    *keylPtr,
                                  char       *key));

extern int
TclX_KeyedListGetKeys _ANSI_ARGS_((Tcl_Interp *interp,
                                   Tcl_Obj    *keylPtr,
                                   char       *key,
                                   Tcl_Obj   **listObjPtrPtr));

extern void
TclX_KeyedListInit _ANSI_ARGS_((Tcl_Interp *interp));

#endif
