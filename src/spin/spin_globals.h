#ifndef _SPIN_GLOBALS_H_
#define _SPIN_GLOBALS_H_
#include <tcl.h>

/* Main global setup function */
int spin_init_globals(Tcl_Interp *interp);
extern Tcl_Obj *spin_defs;

extern int cutoff_align1;
extern int cutoff_align2;
extern int cutoff_align3;
extern char *symbol_align0;
extern char *symbol_align1;
extern char *symbol_align2;
extern char *symbol_align3;


#endif
