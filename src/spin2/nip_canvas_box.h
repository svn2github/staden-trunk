#ifndef _NIP_CANVAS_BOX_H_
#define _NIP_CANVAS_BOX_H_
int nip_canvas_cursor_refresh(Tcl_Interp *interp, int seq_num,
			      cursor_t *changed_c, cursor_t *canvas_c,
			      CanvasPtr *canvas, win **win_list, int num_wins,
			      int reg_id, int *visible, WorldPtr *world,
			      int cursor_show);
#endif
