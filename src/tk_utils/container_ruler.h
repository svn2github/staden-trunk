#include <tcl.h>
#include "container.h"

int add_length_ruler(Tcl_Interp *interp,
		     container *c,
		     int row_index,
		     int column_index,		     
		     int row_num,
		     int column_num,
		     int orientation);

int add_element_ruler(Tcl_Interp *interp,
		      container *c,
		      element *e,
		      int orientation);

void update_length_ruler(Tcl_Interp *interp,
			 container *c,
			 coord *coords);

