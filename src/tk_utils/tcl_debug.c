#include <tcl.h>
#include <tcl_utils.h>		/* For GetInterp() */


/*
 *-----------------------------------------------------------------------------
 * Some handy debugging routines
 *-----------------------------------------------------------------------------
 */

/*
 * Displays the Tcl stack frame - doesn't work always unfortunately.
 */
void dump_tcl_stack(void) {

    /* We use VarEval as it can take a non writable string, unlike Eval */

    Tcl_VarEval( GetInterp(),
		"for {set i [info level]} {$i > 0} {incr i -1} {"
		"    puts \"Level $i: [info level $i]\""
		"}",
		NULL);
}
