#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>

#include "misc.h"
#include "tg_gio.h"

#include "gap_cli_arg.h"

/* Uncomment the next line to turn on logging for most of the tcl -> c 
   interface.  Probably too much information for normal use. */
/* #define LOG_FILE */
#ifdef LOG_FILE
#   include "text_output.h"
#endif

static void parse_args_set(cli_args *a, char *store, char *val) {
    if (a->type == ARG_STR) {
	*((char **)&store[a->offset]) = val;
    } else if (a->type == ARG_ARR) {
	strncpy((char *)&store[a->offset], val, a->value-1);
    } else if (a->type == ARG_IO) {
	return; /* no default possible */
    } else if (a->type == ARG_INT) {
	*((int *)&store[a->offset]) = atoi(val);
    } else if (a->type == ARG_REC) {
	*((tg_rec *)&store[a->offset]) = atorec(val);
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

static void parse_args_obj_set(cli_args *a, char *store, Tcl_Obj *val) {
    if (a->type == ARG_OBJ) {
	*((Tcl_Obj **)&store[a->offset]) = val;
    } else if (a->type == ARG_STR) {
	*((char **)&store[a->offset]) = Tcl_GetStringFromObj(val, NULL);
    } else if (a->type == ARG_IO) {
	GapIO *io = io_from_obj(val);
	*((GapIO **)&store[a->offset]) = io;
    } else if (a->type == ARG_INT) {
	int i;

	if (Tcl_GetIntFromObj(NULL, val, &i) == TCL_OK) {	
	    *((int *)&store[a->offset]) = i;
	} else {
	    *((int *)&store[a->offset]) =
		atoi(Tcl_GetStringFromObj(val, NULL));
	}
    } else if (a->type == ARG_REC) {
	Tcl_WideInt i;

	if (Tcl_GetWideIntFromObj(NULL, val, &i) == TCL_OK) {	
	    *((tg_rec *)&store[a->offset]) = i;
	} else {
	    *((tg_rec *)&store[a->offset]) =
		atorec(Tcl_GetStringFromObj(val, NULL));
	}
    } else if (a->type == ARG_FLOAT) {
	double d;

	if (Tcl_GetDoubleFromObj(NULL, val, &d) == TCL_OK) {
	    *((float *)&store[a->offset]) = d;
	} else {
	    *((float *)&store[a->offset]) =
		atof(Tcl_GetStringFromObj(val, NULL));
	}
    } else if (a->type == ARG_DOUBLE) {
	double d;

	if (Tcl_GetDoubleFromObj(NULL, val, &d) == TCL_OK) {
	    *((double *)&store[a->offset]) = d;
	} else {
	    *((double *)&store[a->offset]) =
		atof(Tcl_GetStringFromObj(val, NULL));
	}
    } else {
	fprintf(stderr, "Unknown argument type %d\n", a->type);
    }

    a->def = ""; /* mark as used */
}

int gap_parse_obj_args(cli_args *args, void *store, int objc,
		       Tcl_Obj * const objv[]) {
    cli_args *a;
    int ret;

#ifdef LOG_FILE
    tcl_log_str(NULL, NULL, objc, objv);
#endif
    /* Initialise defaults */
    for (a = args; a->command; a++)
	if (a->def) /* String default values */
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
