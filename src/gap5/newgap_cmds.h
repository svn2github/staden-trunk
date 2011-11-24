#ifndef _NEWGAP_CMDS_H
#define _NEWGAP_CMDS_H

#include <stdio.h>
#include <tk.h>
#include "tg_gio.h"
#include "list.h"
#include "gap_globals.h"

/* a view is defined as a distinct window */
extern int current_view; /* global value indicating the current view */

int NewGap_Init(Tcl_Interp *interp);

#endif /* _NEWGAP_CMDS_H */
