#include <tk.h>

/*****************************************************************************/
/*                            C2TclArray                                     */
/*****************************************************************************/
/* convert C array into tcl array */
/* elment_array: C array */
/* num_elements: number of elements in C array */
/* tcl_array: tcl array name */
/* tcl_num_ele: tcl variable to contain number of elements in array */

int
C2TclArray(Tcl_Interp *interp,                                    /* in, out */
	   char *element_array[],                                      /* in */
	   int num_elements,                                           /* in */
	   char *tcl_array,                                            /* in */
	   char *tcl_num_ele);                                         /* in */
