#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "misc.h"
#include "cli_arg.h"
#include "tk.h"

static void parse_args_set(cli_args *a, char *store, char *val) {
    if (a->type == ARG_STR) {
	*((char **) &store[a->offset]) = val;
    } else if (a->type == ARG_ARR) {
	strncpy((char *) &store[a->offset], val, a->value-1);
    } else if (a->type == ARG_INT) {
	*((int *) &store[a->offset]) = atoi(val);
    } else if (a->type == ARG_FLOAT) {
	*((float *) &store[a->offset]) = atof(val);
    } else { /* ARG_DOUBLE */
	*((double *) &store[a->offset]) = atof(val);
    }

    a->def = ""; /* mark as used */
}

/*
 * Parses command line arguments.
 * 'args' specifies an array cli_args. Each cli_arg contains a default. Each
 * argument with a default of NULL are non optional. This function will return
 * -1 if any of these arguments aren't specified, or if unknown arguments
 * are specified. Returns 0 otherwise.
 */

int parse_args(cli_args *args, void *store, int argc, char **argv) {
    int i, ret = 0;
    cli_args *a;

    /* Initialise defaults */
    for (a = args; a->command; a++)
	if (a->def)
	    parse_args_set(a, store, a->def);
	else if (a->type == ARG_ARR)
	    memset((char *)&((char *)store)[a->offset], 0, a->value); /* YUK */

    for (i = 1; i < argc; i++) {
	for (a = args; a->command; a++) {
	    if (strcmp(a->command, argv[i])== 0) {
		if (a->value) {
		    if (i == argc - 1) {
			/*
			  verror(ERR_WARN, "parse_args",
			  "No argument given for option '%s'\n",
			  argv[i]);
			*/
			ret = -1;

			break;
		    }
		    
		    parse_args_set(a, store, argv[++i]);
		} else
		    parse_args_set(a, store, "1");

		break;
	    }
	}

	if (!a->command) {
	    /*
	      verror(ERR_WARN,"parse_args", "Unknown option '%s'\n", argv[i]);
	    */
	    ret = -1;
	}
    }

    /* Check all required options were set */
    for (a = args; a->command; a++) {
	if (!a->def)
	    return -1;
    }
    return ret;
}
