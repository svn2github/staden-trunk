#ifndef _STOP_CODON_H_
#define _STOP_CODON_H_

#include "canvas_box.h"

typedef struct obj_codon_t {
    Tcl_Interp *interp;
    int num_codons;
    char **codon;
    int num_match;
    s_match c_match;
    int sequence_type;
    int strand;
    int start;
    int end;
    int yoffset;
    tick_s *tick;
    cursor_s xhair;
    int text_offset;
    int editor_avail;

    int id;
    char window[100];    /* stop codon plot window */
    char names_win[100]; /* frame name window */
    char frame[100];     /* parent frame */
    ruler_s *ruler;
    win **win_list;
    int num_wins;
    WorldPtr *world;
    CanvasPtr *canvas;
    StackPtr *zoom;
    cursor_t *cursor;
    int cursor_visible;
} obj_codon;


int
codon_reg(Tcl_Interp *interp,
	  int strand,
	  GapIO *io,
	  char *frame,
	  char *names,
	  char *plot,
	  int contig_num,
	  int lreg,
	  int rreg,
	  tick_s *tick,
	  int yoffset,
	  ruler_s *ruler,
	  cursor_s cursor);

int
stop_codon_replot(Tcl_Interp *interp,
		  GapIO *io, 
		  obj_codon *s,
		  char *seq);

#endif
