#ifndef _CONSISTENCY_CANVAS_BOX_H
#define _CONSISTENCY_CANVAS_BOX_H

#include "consistency_display.h"
#include "canvas_box.h"

int consistency_cursor_refresh(Tcl_Interp *interp, GapIO *io, 
			       obj_consistency_disp *c, int cnum,
			       cursor_t *changed_c, cursor_t *canvas_c,
			       CanvasPtr *canvas, win **win_list, int num_wins,
			       int reg_id, int offset, int *visible,
			       WorldPtr *world, int cursor_show);
#endif

