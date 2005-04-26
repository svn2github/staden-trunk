#ifndef _TK_EDITOR_H
#define _TK_EDITOR_H

#include <tk.h>
#include "tkSheet.h"
#include "tkSheet_common.h"
#include "sheet.h"


/* We don't need the definition here, but just need to know it exists */
struct _EdStruct;

typedef struct {
#   include "tkSheet_struct.h"

    XColor *qual_bg[10];
    XColor *qual_below;
    XColor *diff_bg;
    XColor *edit_bg[4];
    XColor *tmpl_bg[6];
    XColor *set_bg[10];
    char *xScrollCmd;
    char *yScrollCmd;
    char *highlight_cmd;
    int max_height;

    /* And our data to edit */
    struct _EdStruct *xx;
} Editor;

#define TKSHEET(ed)   ((tkSheet *)(ed))
#define EDTKWIN(ed)   ((ed)->sw.tkwin)
#define EDINTERP(ed)  (TKSHEET(ed)->interp)
#define EDDISPLAY(ed) ((ed)->sw.display)

int Editor_Init(Tcl_Interp *interp);
void ed_set_slider_pos(struct _EdStruct *xx, int pos);
int edGetSelection(ClientData clientData, int offset, char *buffer,
		   int bufsize);

#endif
