#ifndef _NIP_CANVAS_BOX_H_
#define _NIP_CANVAS_BOX_H_

#include "editor_reg.h"

typedef struct {
    int seq_id;
    int start;
    int end;
    char *renz_start;
    char *renz_end;
} seq_copy_arg;

/*folowing struct come from "seq_reg_structs.h", due to same struct has been
  defined in "nip_structs.h", so "seq_reg_structs.h" can not be included.*/
typedef struct {
    int seq_num;
    int id;
    int pos;
    int direction;
} ed_cursor_notify_arg;

typedef struct {
    int cursorid;
    int seq_num;
} eqc_arg;


int ed_nip_canvas_cursor_refresh(Tcl_Interp *interp, 
			      int seq_num,
			      cursor_e *changed_c, 
			      cursor_e *canvas_c,
			      CanvasPtr *canvas, win **win_list, int num_wins,
			      int reg_id, int *visible, WorldPtr *world,
			      int cursor_show);

int RenzCanvas_Init(Tcl_Interp *interp);

#endif
