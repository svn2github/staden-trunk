#ifndef _FT_VIEWER_H_
#define _FT_VIEWER_H_

#include "seq_reg.h"
#include "canvas_box.h"

#define FT_RULER   0
#define FT_FORWARD 1
#define FT_REVERSE 2

#define FT_SINGLE 1  /* draw all features on single line */
#define FT_FEATURE 2 /* draw same features on same line */
#define FT_ALL 3     /* draw all features on separate lines */

typedef struct ft_viewer_in {
    char *params;
    int start;
    int end;
} in_ft_viewer;

typedef struct feats_ {
    double origin;
    int display;
} feats;

typedef struct feat_s_ {
    char *win;    /* window name */
    int display;  /* display mode */
    double offset; /* offset of separated features */
    int min; /* min height */
    int max; /* max height */
} feat_s;

typedef struct ft_viewer_res_ {
    int plot_type; /* 0 = horizontal, 1 = vertical, 2 = circular */
    int ft_offset;
    int ft_height;
    double origin;           /* origin of circle */
    feat_s forward;
    feat_s reverse;

    int sequence_type; /* HACK to do something with */
    char *frame;
    cursor_s cursor;
    ruler_s *ruler; 
    int sequence_len;
    win **win_list;
    int num_wins;
    WorldPtr *world;
    CanvasPtr *canvas;
    StackPtr *zoom;
} ft_viewer_res;

typedef struct ft_key_s_ {
    int index;
    char *bg_colour;
    char *shape;
} ft_key_s;

int ft_viewer_reg(Tcl_Interp *interp,
		  int seq_id,
		  char *frame,
		  char *forward,
		  char *reverse,
		  int start,
		  int end,
		  int plot_type,
		  int ft_offset,
		  int ft_height,
		  ruler_s *ruler,
		  cursor_s cursor,
		  double origin);

int init_ft_viewer_create(Tcl_Interp *interp,
			  int seq_id,
			  int start,
			  int end,
			  int strand,
			  int orientation,
			  Tcl_Obj **graph_obj,
			  int *id);

int init_ft_viewer_plot(Tcl_Interp *interp,
			int seq_id,
			int result_id,
			char *e_win,
			char *c_win, 
			Tcl_Obj *results,
			int container_id,
			int element_id,
			char *colour,
			int orientation,
			int ft_height,
			float ft_offset,
			char *element_type,
			double origin,
			int display_mode);
#endif

