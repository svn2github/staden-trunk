#include <tk.h>
#include "canvas_box.h"

void display_ruler_ticks(Tcl_Interp *interp, CanvasPtr *canvas, int x_offset,
			 int y_offset, ruler_s *ruler, int start, int end);

void display_ruler_ticks_v(Tcl_Interp *interp, CanvasPtr *canvas, 
			   ruler_s *ruler, double wy1, double wy2);

void display_ruler_ticks_c(Tcl_Interp *interp,  CanvasPtr *canvas, 
			   ruler_s *ruler, int wx1, int wx2, double origin,
			   cir_s circle);

void ruler_ticks(double ruler_min, double ruler_max, int def_num_ticks,
		 double *firstTick, double *step, int *numTicks);
