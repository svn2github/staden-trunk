#if !defined(HANG_EDITOR_H)
#define HANG_EDITOR_H

#include "editor.h"

int Hang_editor_Init(Tcl_Interp *interp);
int un_trim_sequence (EDITOR_RECORD *er, EDIT *edit, int add_to_undo);
int un_fill_sequence (EDITOR_RECORD *er, EDIT *edit, int add_to_undo);

extern  EDITOR_RECORDS *editor_records;

#endif
