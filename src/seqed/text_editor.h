#if !defined(TEXT_EDITOR_H)
#define TEXT_EDITOR_H

#include "editor.h"


int TextEditor_Init (Tcl_Interp *interp);
int remove_edit(EDITS *edits);

extern  EDITOR_RECORDS *editor_records;

#endif
