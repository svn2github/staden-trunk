#if !defined(RENZYME_H)
#define RENZYME_H

#include "canvas_box.h"
#include "cli_arg.h"
#include "renzyme_box.h"
#include "renzyme_search.h"
#include "editor.h"

/* this comes from SPIN out_canvas */
typedef struct out_canvas_e_ {
  Tcl_Interp *interp;
  cursor_e *cursor;
  int cursor_visible;
} out_canvas_e;

typedef struct {
    char *frame;
    char *win_name;
    char *plot;
    char *win_ruler;
    char *inlist;
    int  num_items;
    int  text_offset;
    char *text_fill;
    int  tick_ht;
    int  tick_wd;
    char *tick_fill;
    int  cursor_wd;
    char *cursor_fill;
    int  yoffset;
    int  seq_id;
    int  start;
    int  end;
} ed_renz_arg;

typedef struct ed_renz_res_ {
    int num_enzymes;
    RENZYMES *r_enzyme;
    int num_match;
    R_MATCH *match;
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
    int seq_id;
    SELECTION *sel;
} ed_renz_res;


int EditorPlotRenz(ClientData cData, Tcl_Interp *interp, int argc, char **argv);
int RenzymeCmds_Init(Tcl_Interp *interp);
RENZYME *renzyme_copy ( RENZYME *r);
int find_matches(RENZYME *rs,                                          /* in */
		 char *sequence,                                       /* in */
		 int sequence_len,                                     /* in */
		 int sequence_type,                                    /* in */
		 int rs_id,                                            /* in */
		 R_MATCH **match,                                      /* out */
		 int *total_matches);                                  /* out */

int find_matches_all (RENZYMES *rss,                                    /* in */
		      char *sequence,                                   /* in */
		      int sequence_len,                                 /* in */
		      int sequence_type,                                /* in */
		      R_MATCH **match,                                  /* out */
		      int *total_matches);                              /* out */

#endif
