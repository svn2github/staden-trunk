#ifndef _SPIN_CLI_ARG_H_
#define _SPIN_CLI_ARG_H_

#include "cli_arg.h"
#define ARG_OBJ 7

/*
 * Parses command line arguments.
 * 'args' specifies an array cli_args. Each cli_arg contains a default. Each
 * argument with a default of NULL are non optional. This function will return
 * -1 if any of these arguments aren't specified, or if unknown arguments
 * are specified. Returns 0 otherwise.
 */
int spin_parse_args(cli_args *args, void *store, int argc, char **argv);
int spin_parse_obj_args(cli_args *args, void *store, int objc,
		       Tcl_Obj * const objv[]);

#endif
