/*
 * tclXutil.c
 *
 * Utility functions for Extended Tcl.
 *-----------------------------------------------------------------------------
 * Copyright 1991-1997 Karl Lehenbauer and Mark Diekhans.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose and without fee is hereby granted, provided
 * that the above copyright notice appear in all copies.  Karl Lehenbauer and
 * Mark Diekhans make no representations about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *-----------------------------------------------------------------------------
 * $Id: tclXutil.c,v 1.1.1.1 2003-06-09 11:24:38 jkb Exp $
 *-----------------------------------------------------------------------------
 */

#include "tclExtdInt.h"

/*
 * Used to return argument messages by most commands.
 */
char *tclXWrongArgs = "wrong # args: ";


/*-----------------------------------------------------------------------------
 * TclX_WrongArgs --
 *
 *   Easily create "wrong # args" error messages.
 *
 * Parameters:
 *   o commandNameObj - Object containing name of command (objv[0])
 *   o string - Text message to append.
 * Returns:
 *   TCL_ERROR
 *-----------------------------------------------------------------------------
 */
int TclX_WrongArgs( Tcl_Interp *interp, Tcl_Obj *commandNameObj, char *string )
{
    char    *commandName;
    Tcl_Obj *resultPtr = Tcl_GetObjResult (interp);
    int      commandLength;

    commandName = Tcl_GetStringFromObj (commandNameObj, &commandLength);

    Tcl_AppendStringsToObj (resultPtr,
			    tclXWrongArgs,
			    commandName,
			    (char *)NULL);

    if (*string != '\0') {
	Tcl_AppendStringsToObj (resultPtr, " ", string, (char *)NULL);
    }
    return TCL_ERROR;
}


/*-----------------------------------------------------------------------------
 * TclX_AppendObjResult --
 *
 *   Append a variable number of strings onto the object result already
 * present for an interpreter.
 *
 * Parameters:
 *   o interp - Interpreter to set the result in.
 *   o args - Strings to append, terminated by a NULL.
 *-----------------------------------------------------------------------------
 */
void TclX_AppendObjResult( Tcl_Interp *interp, char *arg1, ... )
{
    Tcl_Obj *resultPtr;
    va_list argList;
    char *string;

    va_start(argList, arg1);
    resultPtr = Tcl_GetObjResult (interp);

    while (1) {
        string = va_arg(argList, char *);
        if (string == NULL) {
            break;
        }
        Tcl_AppendToObj (resultPtr, string, -1);
    }
    va_end(argList);
}


/*-----------------------------------------------------------------------------
 * TclX_IsNullObj --
 *
 *   Check if an object is {}, either in list or zero-lemngth string form, with
 * out forcing a conversion.
 *
 * Parameters:
 *   o objPtr - Object to check.
 * Returns:
 *   True if NULL, FALSE if not.
 *-----------------------------------------------------------------------------
 */
int TclX_IsNullObj( Tcl_Obj *objPtr )
{
    static Tcl_ObjType *listType = NULL, *stringType = NULL;
    int length;

    /*
     * Only get types once, as they must be static.
     */
    if (listType == NULL) {
        listType = Tcl_GetObjType ("list");
        stringType = Tcl_GetObjType ("string");
    }

    if (objPtr->typePtr == NULL) {
        return (objPtr->length == 0);
    } else {
        if (objPtr->typePtr == listType) {
            Tcl_ListObjLength (NULL, objPtr, &length);
            return (length == 0);
        } else if (objPtr->typePtr == stringType) {
            Tcl_GetStringFromObj (objPtr, &length);
            return (length == 0);
        }
    }
    Tcl_GetStringFromObj (objPtr, &length);
    return (length == 0);
}
