#ifndef _RESTRICTION_ENZYME_MAP_H_
#define _RESTRICTION_ENZYME_MAP_H_

#include "canvas_box.h"
#include "renz_utils.h"

void renz_shutdown(R_Enz *r_enzyme, 
		   int num_enzymes,
		   R_Match *match,
		   CanvasPtr *canvas,
		   WorldPtr *world,
		   StackPtr *zoom);

void
PlotStickMap(Tcl_Interp *interp,
	     char *win_name,
	     int unpadded_cut_site,
	     int cut_site,
	     int xoffset,
	     int yoffset,
	     int tick_ht,
	     int tick_wd,
	     char *colour,
	     int index,
	     int from,
	     int to);

void plot_renz_matches(Tcl_Interp *interp,
		       char *re_win,
		       char *names_win,
		       int text_offset,
		       char *t_colour,
		       int yoffset,
		       int num_enzymes,
		       R_Enz *r_enzyme,
		       ruler_s *ruler,
		       int sequence_len,
		       int num_matches,
		       R_Match *match,
		       tick_s *tick,
		       char *frame,
		       WorldPtr *world,
		       CanvasPtr *canvas,
		       win **win_list,
		       int num_wins,
		       StackPtr **zoom);
#endif
