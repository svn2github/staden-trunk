#ifndef _RULER_DISPLAY_H
#define _RULER_DISPLAY_H

#include "canvas_box.h"
#include "IO.h"

int display_ruler(Tcl_Interp *interp,  GapIO *io, CanvasPtr *canvas,
		  c_offset *contig_offset, int *contig_array,
		  int num_contigs, int disp_ruler, int disp_ticks,
		  ruler_s *ruler, char *frame, PlotRec **ruler_coord);

int
display_single_ruler(Tcl_Interp *interp, GapIO *io, CanvasPtr *canvas,    
		     int contig_array, int disp_ruler, int disp_ticks,
		     ruler_s *ruler);

void display_ruler_v(Tcl_Interp *interp, CanvasPtr *canvas, ruler_s *ruler,
		     double wy1, double wy2);
#endif
