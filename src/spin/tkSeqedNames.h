#ifndef _TK_EDNAMES_H
#define _TK_EDNAMES_H

#include <tk.h>
#include "sheet.h"


/* We don't need the definition here, but just need to know it exists */
/* struct _EdStruct; */

typedef struct {
#   include "tkSheet_struct.h"
    
    /* And our data to edit */
    /* struct _EdStruct *xx; */
} seqedNames;
int SeqedNames_Init(Tcl_Interp *interp);

#endif
