#include "gap_array.h"
#include "active_tags.h"
#include "tcl_utils.h"

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
	   char *tcl_num_ele)                                          /* in */
{
  int item;
  char item_str[1024];
  char num_ele_str[1024];
  
  /* convert num_elements into a string */
  sprintf(num_ele_str, "%d", num_elements);

  /* create tcl variable passed as tcl_num_ele to hold num_elements */
  if (Tcl_SetVar(interp, tcl_num_ele, num_ele_str, 0) == NULL) {

      return TCL_ERROR; 
      
  } /* end if */

  for (item = 0; item < num_elements; item++) {

      /* convert each element of the array into a string for Tcl_SetVar */
      sprintf(item_str, "%d", item);

      /* create tcl array and set the corresponding C item of each element of
       * the array 
       */
      if (Tcl_SetVar2(interp, tcl_array, item_str, element_array[item], 0) 
	  ==  NULL) {  
	  Tcl_SetResult(interp, "Tcl_SetVar2 failed", TCL_STATIC);
	  return TCL_ERROR; 

      } /* end if */

  } /* end for */

  return TCL_OK;

} /* end C2TclArray */


