#ifndef _NIP_RESTRICTION_ENZYMES_H_
#define _NIP_RESTRICTION_ENZYMES_H_

#include "seq_reg.h"

typedef struct renz_res_ {
    int num_enzymes;
    R_Enz *r_enzyme;
    int num_match;
    R_Match *match;

    int yoffset;
    tick_s *tick;
    cursor_s cursor;
    int text_offset;
    char *text_colour;

    int sequence_type; /* HACK to do something with */
    char re_win[100];
    char names_win[100];
    char frame[100];
    ruler_s *ruler; 
    int sequence_len;
    win **win_list;
    int num_wins;
    WorldPtr *world;
    CanvasPtr *canvas;
    StackPtr *zoom;
} renz_res;

int nip_renz_reg(Tcl_Interp *interp, int seq_id, out_canvas *output, 
		 char *filename,
		 char *frame, char *names_win, char *re_win, char *inlist,
		 int num_items, int start, int end, int text_offset,
		 char *text_fill, tick_s *tick, int yoffset, ruler_s *ruler,
		 cursor_s cursor);

void nip_renz_info(int seq_num, renz_res *data, int start, int print_opt);

#endif
