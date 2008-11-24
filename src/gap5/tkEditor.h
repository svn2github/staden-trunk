#ifndef _TK_EDITOR_H
#define _TK_EDITOR_H

#include <tk.h>
#include "tkSheet.h"
#include "tkSheet_common.h"
#include "sheet.h"


/* We don't need the definition here, but just need to know it exists */
struct _edview;

typedef struct {
#   include "tkSheet_struct.h"

    XColor *qual_bg[10];
    XColor *qual_below;
    XColor *diff1_bg;
    XColor *diff2_bg;
    XColor *diff1_fg;
    XColor *diff2_fg;
    XColor *edit_bg[4];
    XColor *tmpl_bg[6];
    XColor *set_bg[10];
    char *xScrollCmd;
    char *yScrollCmd;
    char *highlight_cmd;
    int max_height;
    
    /* And our data to edit */
    struct _edview *xx;

    /* Display options */
    int display_cutoffs;
    int display_quality;
    int display_mapping_quality;
    int display_differences; /* 0 = no, 1, 2, 3 = display modes */
    int display_differences_case;
    int display_differences_qual;
    int consensus_at_top; /* 0 => normal, 1 => flipped for join editor */
} Editor;

#define TKSHEET(ed)   ((tkSheet *)(ed))
#define EDTKWIN(ed)   ((ed)->sw.tkwin)
#define EDINTERP(ed)  (TKSHEET(ed)->interp)
#define EDDISPLAY(ed) ((ed)->sw.display)

#include "editor_view.h"

int Editor_Init(Tcl_Interp *interp);
void ed_set_slider_pos(struct _edview *xx, int pos);
int edGetSelection(ClientData clientData, int offset, char *buffer,
		   int bufsize);

#endif
