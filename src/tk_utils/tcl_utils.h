#ifndef _TCL_UTILS_H
#define _TCL_UTILS_H

#include <tcl.h>
#include "misc.h" /* for __PRINTF_FORMAT__ */

#ifdef _MSC_VER
#  ifdef BUILDING_TK_UTILS_DLL
#    define TK_UTILS_EXPORT __declspec(dllexport)
#  else
#    define TK_UTILS_EXPORT __declspec(dllimport)
#  endif
#else
#  define TK_UTILS_EXPORT
#endif


extern TK_UTILS_EXPORT Tcl_Obj *tk_utils_defs;   /* Todo: Make this a function */


Tcl_Interp* GetInterp( void );
char* GetInterpResult( void );
void dump_tcl_stack(void);
void vTcl_SetResult(Tcl_Interp *interp, char *fmt, ...) __PRINTF_FORMAT__(2,3);
char *vTcl_DStringAppend(Tcl_DString *dsPtr, char *fmt, ...) __PRINTF_FORMAT__(2,3);
char *vTcl_DStringAppendElement(Tcl_DString *dsPtr, char *fmt, ...) __PRINTF_FORMAT__(2,3);
char *w(char *str);
char *vw(char *fmt, ...) __PRINTF_FORMAT__(1,2);
void *TclPtr2C(char *str);
char *CPtr2Tcl(void *ptr);
char *get_default_string(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *name);
char *get_default_astring(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *name);
int get_default_int(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *name);
double get_default_double(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *name);
Tcl_Obj *get_default_obj(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *name);


#endif
