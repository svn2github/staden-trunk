#ifndef _GAP_CANVAS_BOX_H
#define _GAP_CANVAS_BOX_H

#include "tg_gio.h"
#include "canvas_box.h"

/* TASK_CANVAS_WORLD structure */
typedef struct {
    int canvasx;  /* input */
    int cnum;     /* input */
    double basex; /* returned */
} task_world_t;

int canvas_cursor_refresh(Tcl_Interp *interp, GapIO *io, int cnum,
			  cursor_t *changed_c, cursor_t *canvas_c,
			  CanvasPtr *canvas, win **win_list, int num_wins,
			  int reg_id, int offset, int *visible,
			  WorldPtr *world, int cursor_show);

int find_cursor_contig(GapIO *io, int id, c_offset *contig_offset,
		       int *contig_array, int num_contigs, double wx);

/*
 * Deletes an existing cursor - calls canvas_cursor_delete in tcl
 */
void canvas_cursor_delete(Tcl_Interp *interp, GapIO *io, cursor_t *cursor,
			  CanvasPtr *canvas, win **win_list, int num_wins);

#endif
