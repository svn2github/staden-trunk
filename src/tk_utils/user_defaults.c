/* NB: Not used; not finished; needs updated to tclX v8 */

#include <tcl.h>
#include "tclXkeylist.h"

static char *def_list = NULL;

int tcl_set_def(ClientData clientData, Tcl_Interp *interp,
		int argc, char **argv) {
    char *new_list;

    /*
     * proc set_def {a1 a2} {
     *    global defn
     *    global $defn
     *    keylset $defn $a1 $a2
     * }
     */

    if (argc != 3) {
	Tcl_AppendResult(interp, "wrong # args: ", argv[0],
			 " name value", NULL);
	return TCL_ERROR;
    }

    new_list = Tcl_SetKeyedListField(interp, argv[1], argv[2], def_list);
    if (NULL == new_list) {
	Tcl_AppendResult(interp, "out of memory", NULL);
	return TCL_ERROR;
    }

    if (def_list)
	xfree(def_list);
    def_list = new_list;

    return TCL_OK;
}
