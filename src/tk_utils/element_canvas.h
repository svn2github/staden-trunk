#ifndef _ELEMENT_CANVAS_H_
#define _ELEMENT_CANVAS_H_

#include <tcl.h>
#include "container.h"
#include "tclCanvGraph.h"

double canvas_x(Tcl_Interp *interp, char *win, double val);
double canvas_y(Tcl_Interp *interp, char *win, double val);

int canvas_width(Tcl_Interp *interp, char *win);
int canvas_height(Tcl_Interp *interp, char *win);


void canvas_scroll_x(Tcl_Interp *interp, element *e, char *command);

void canvas_scroll_y(Tcl_Interp *interp, element *e, char *command);

int create_graph(Tcl_Interp *interp, char *e_win, Tcl_Obj *results,
		 int width, char *foreground, char *tags, int orientation);
void create_canv_line(Tcl_Interp *interp,
		      char *e_win,
		      Graph *graph,plot_data *result,
		      int width,
		      char *foreground,
		      char *tags,
		      int orientation);

void canvas_scale(Tcl_Interp *interp, element *e, int result_id, 
		  d_box *current, CanvasPtr *canvas);

void canvas_scrollregion(Tcl_Interp *interp, element *e, d_box *total,
			 CanvasPtr *pixel_column, CanvasPtr *pixel_row);

void canvas_zoom(Tcl_Interp *interp, element **e_list, int num_e);

void draw_canvas_crosshair(Tcl_Interp *interp, element *e, int pos,
			   int orientation);

void delete_canvas_crosshair(Tcl_Interp *interp, element *e);
void canvas_move(Tcl_Interp *interp,
		 element *e,
		 int result_id,
		 double x_amount,
		 double y_amount);
#endif
