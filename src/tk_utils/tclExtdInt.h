#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <tcl.h>
#include <stdarg.h>
#include <assert.h>

#define FALSE 0
#define TRUE 1

/*
 * Macro that behaves like strdup, only uses ckalloc.  Also macro that does the
 * same with a string that might contain zero bytes,
 */
#define ckstrdup(sourceStr) \
  (strcpy (ckalloc (strlen (sourceStr) + 1), sourceStr))

#define ckbinstrdup(sourceStr, length) \
  ((char *) memcpy (ckalloc (length + 1), sourceStr, length + 1))


/*
 * Prototypes for TCL extension utilities
 */
int  TclX_IsNullObj( Tcl_Obj *objPtr );
void TclX_AppendObjResult( Tcl_Interp *interp, char *arg1, ... );
int  TclX_WrongArgs( Tcl_Interp *interp, Tcl_Obj *commandNameObj, char *string );
