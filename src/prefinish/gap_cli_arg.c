#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>

#include "misc.h"
#include "io_utils.h"	/* io_handle() */

/*
typedef int GapIO;
GapIO *io_handle(int *x) {}
*/

#include "gap_cli_arg.h"

static void parse_args_set(cli_args *a, char *store, char *val) {
    if (a->type == ARG_STR) {
	*((char **)&store[a->offset]) = val;
    } else if (a->type == ARG_ARR) {
	strncpy((char *)&store[a->offset], val, a->value-1);
    } else if (a->type == ARG_IO) {
	int handle = atoi(val);
	GapIO *io = io_handle(&handle);

	if (io)
	    *((GapIO **)&store[a->offset]) = io;
	else
	    return; /* Doesn't set a->def */
    } else if (a->type == ARG_INT) {
	*((int *)&store[a->offset]) = atoi(val);
    } else if (a->type == ARG_OBJ) {
	/* Cannot initialise objs */
	*((Tcl_Obj **)&store[a->offset]) = NULL;
    } else if (a->type == ARG_FLOAT) {
	*((float *)&store[a->offset]) = atof(val);
    } else if (a->type == ARG_DOUBLE) {
	*((double *)&store[a->offset]) = atof(val);
    } else {
	fprintf(stderr, "Unknown argument type %d\n", a->type);
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
int gap_parse_args(cli_args *args, void *store, int argc, char **argv) {
    cli_args *a;
    int ret;

    /* Initialise defaults */
    for (a = args; a->command; a++)
	if (a->def)
	    parse_args_set(a, store, a->def);
	else if (a->type == ARG_ARR)
	    memset((char *)&((char *)store)[a->offset], 0, a->value); /* YUK */

    ret = gap_parse_config(args, store, argc, argv);
    
    /* Check all required options were set */
    for (a = args; a->command; a++)
	if (!a->def) {
	    verror(ERR_WARN, "parse_args",
		   "No argument given for option '%s'\n",
		   a->command);
	    return -1;
	}

    return ret;
}

int gap_parse_config(cli_args *args, void *store, int argc, char **argv) {
    int i, ret = 0;
    cli_args *a;
    char *aname;

    for (i = 1; i < argc; i++) {
	aname = argv[i];

	for (a = args; a->command; a++) {
	    if (strcmp(a->command, aname)== 0) {
		if (a->value) {
		    if (i == argc - 1) {
			verror(ERR_WARN, "parse_args",
			       "No argument given for option '%s'\n",
			       aname);
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
	    verror(ERR_WARN,"parse_args", "Unknown option '%s'\n", aname);
	    ret = -1;
	}
    }

    return ret;
}


static void parse_args_obj_set(cli_args *a, char *store, Tcl_Obj *val) {
    if (a->type == ARG_OBJ) {
	*((Tcl_Obj **)&store[a->offset]) = val;
    } else if (a->type == ARG_STR) {
	*((char **)&store[a->offset]) = Tcl_GetStringFromObj(val, NULL);
    } else if (a->type == ARG_IO) {
	int handle;
	GapIO *io;

	Tcl_GetIntFromObj(NULL, val, &handle);
	io = io_handle(&handle);

	if (io)
	    *((GapIO **)&store[a->offset]) = io;
	else
	    return; /* Doesn't set a->def */
    } else if (a->type == ARG_INT) {
	int i;

	Tcl_GetIntFromObj(NULL, val, &i);
	*((int *)&store[a->offset]) = i;
    } else if (a->type == ARG_FLOAT) {
	double d;

	Tcl_GetDoubleFromObj(NULL, val, &d);
	*((float *)&store[a->offset]) = d;
    } else if (a->type == ARG_DOUBLE) {
	double d;

	Tcl_GetDoubleFromObj(NULL, val, &d);
	*((double *)&store[a->offset]) = d;
    } else {
	fprintf(stderr, "Unknown argument type %d\n", a->type);
    }

    a->def = ""; /* mark as used */
}

int gap_parse_obj_args(cli_args *args, void *store, int objc,
		       Tcl_Obj * const objv[]) {
    cli_args *a;
    int ret;

    /* Initialise defaults */
    for (a = args; a->command; a++)
	if (a->def)
	    parse_args_set(a, store, a->def);
	else if (a->type == ARG_ARR)
	    memset((char *)&((char *)store)[a->offset], 0, a->value); /* YUK */

    ret = gap_parse_obj_config(args, store, objc, objv);

    /* Check all required options were set */
    for (a = args; a->command; a++)
	if (!a->def)
	    return -1;

    return ret;
}

int gap_parse_obj_config(cli_args *args, void *store, int objc,
			 Tcl_Obj * const objv[]) {
    int i, ret = 0;
    cli_args *a;
    Tcl_Obj *one_str;

    one_str = Tcl_NewStringObj("1", -1);

    for (i = 1; i < objc; i++) {
	char *aname = Tcl_GetStringFromObj(objv[i], NULL);
	for (a = args; a->command; a++) {
	    if (strcmp(a->command, aname) == 0) {
		if (a->value) {
		    if (i == objc - 1) {
			verror(ERR_WARN, "parse_args",
			       "No argument given for option '%s'\n",
			       aname);
			ret = -1;

			break;
		    }
		    
		    parse_args_obj_set(a, store, objv[++i]);
		} else
		    parse_args_obj_set(a, store, one_str);

		break;
	    }
	}

	if (!a->command) {
	    verror(ERR_WARN,"parse_args", "Unknown option '%s'\n", aname);
	    ret = -1;
	}
    }

    Tcl_DecrRefCount(one_str);

    return ret;
}
