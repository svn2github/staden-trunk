#ifndef _GAP_CLI_ARG_H_
#define _GAP_CLI_ARG_H_

#include "cli_arg.h"
#define ARG_OBJ 7
#define ARG_DBL ARG_DOUBLE

/*
 * Parses command line arguments.
 * 'args' specifies an array cli_args. Each cli_arg contains a default. Each
 * argument with a default of NULL are non optional. This function will return
 * -1 if any of these arguments aren't specified, or if unknown arguments
 * are specified. Returns 0 otherwise.
 *
 * The _args functions set the default values into the store and then run
 * the _config functions, which do the main work of parsing the
 * arguments and filling out the relevant store bits.
 * Ie _args is primarily for argument parsing to functions and _config 
 * is for subsequent configuration of data objects. The _config
 * function ignores the .def section of cli_args.
 */
int gap_parse_args(cli_args *args, void *store, int argc, char **argv);
int gap_parse_config(cli_args *args, void *store, int argc, char **argv);
int gap_parse_obj_args(cli_args *args, void *store, int objc,
		       Tcl_Obj * const objv[]);
int gap_parse_obj_config(cli_args *args, void *store, int objc,
			 Tcl_Obj * const objv[]);

#endif
