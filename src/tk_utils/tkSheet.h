#ifndef _TK_SHEET_H
#define _TK_SHEET_H

#include <tk.h>
#include "sheet.h"

typedef struct _tkSheet {
#   include "tkSheet_struct.h"
} tkSheet;

extern tkSheet *LastSheet(void);

int Sheet_Init(Tcl_Interp *interp);

#endif
