#ifndef _CONSISTENCY_DISPLAY_H
#define _CONSISTENCY_DISPLAY_H

#include "canvas_box.h"
#include "io-reg.h"

/* configure array for what to display on the consistency display */
#define NUM_CONS_CONFIGS 6
#define CONS_RULER       0
#define CONS_TICKS       1
#define CONFIDENCE       2
#define READING_COV      3
#define READPAIR_COV     4
#define CONSISTENCY      5

typedef struct {
    Tcl_Interp *interp;
    c_offset *contig_offset;
    int *contigs;
    int num_contigs;
    int start;
    int end;
    char frame[100];
    int id;
    ruler_s *ruler;
    cursor_s xhair;
    win **win_list;
    int num_wins;
    d_box *orig_total;
    PlotRec *ruler_coord;
    int do_update;
    int buffer_count;
    cursor_t **cursor;
    int *cursor_visible;
    int configs[NUM_CONS_CONFIGS];
} obj_consistency_disp;

int add_consistency_window(GapIO *io,obj_consistency_disp *c,
			   char *window, char scroll, int id,
			   double x1, double y1, double x2, double y2);
void delete_consistency_window(obj_consistency_disp *c, int win_num);
int get_consistency_win_num(obj_consistency_disp *c, int id);
void clear_consistency(GapIO *io, obj_consistency_disp *c);
int consistency_reg(GapIO *io,Tcl_Interp *interp, int *contig_array,
		    int num_contigs, int start, int end, char *frame, 
		    ruler_s *ruler, cursor_s xhair);
void create_consistency_ruler(GapIO *io, obj_consistency_disp *c);
void consistency_update_cursors(GapIO *io, obj_consistency_disp *c,
				int show_cursor);
void consistency_shutdown(GapIO *io, obj_consistency_disp *c);

void consistency_canvasScrollX(Tcl_Interp *interp,
			       obj_consistency_disp *c,
			       win **win_list,
			       int num_wins,
			       char *scroll_args);

/* only scroll single window */
void consistencyScrollY(Tcl_Interp *interp,
			char *window,
			int scroll,
			d_box *visible,
			CanvasPtr *canvas,
			char *scroll_args);

#endif
