#include <tcl.h>
#include <tk.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "misc.h"
#include "getfile.h"
#include "tcl_utils.h"
#include "capture.h"
#include "traceType.h"
#include "cli_arg.h"
#include "renz_utils.h"

/* 7/1/99 johnt - added S_ISDIR and S_ISREG definitions to WINNT support */
#ifndef S_ISDIR
#define S_ISDIR(m)      (((m)&S_IFMT) == S_IFDIR)
#endif /*!S_ISDIR*/
#ifndef S_ISREG
#define S_ISREG(m)      (((m)&S_IFMT) == S_IFREG)
#endif /*!S_ISREG*/


int tcl_expandpath(ClientData clientData, Tcl_Interp *interp,
		   int argc,  char **argv) {
    char exp_filename[FILENAME_MAX+1];

    if (argc != 2)
	return TCL_ERROR;

    if (0 == expandpath(argv[1], exp_filename)) {
	return TCL_ERROR;
    }

    vTcl_SetResult(interp, "%s", exp_filename);
    return TCL_OK;
}

int tcl_createconsole(ClientData clientData, Tcl_Interp *interp,
		      int argc,  char **argv) {
    Tk_CreateConsoleWindow(interp);
    return TCL_OK;
}

int tcl_mkdir(ClientData clientData, Tcl_Interp *interp,
	      int argc,  char **argv) {
    struct stat statbuf;

    if (argc != 2)
	return TCL_ERROR;

    /* open output directory */
    /* check for directory existance, and create if needed */
    if (-1 == stat(argv[1], &statbuf)) {
	if (-1 == mkdir(argv[1], 0777)) {
	    perror(argv[1]);
	    verror(ERR_WARN, "mkdir", "cannot create directory %s", argv[1]);
	    return TCL_ERROR;
	}
    } else {
	if (!S_ISDIR(statbuf.st_mode)) {
	    verror(ERR_WARN, "mkdir", "%s already exists and is not a "
		   "directory", argv[1]);
	    return TCL_ERROR;
	}
    }
    return TCL_OK;
}

int tcl_dir_or_file(ClientData clientData, Tcl_Interp *interp,
		    int argc,  char **argv) {
    Tcl_DString files;
    Tcl_DString dirs;
    Tcl_DString result;
    int i;
    struct stat st;

    if (argc < 2) {
	Tcl_SetResult(interp,
		      "wrong # args: should be \"tcl_dir_or_file "
		      "filename ...\"\n", TCL_STATIC);
	return TCL_ERROR;
    }
    Tcl_DStringInit(&files);
    Tcl_DStringInit(&dirs);
    Tcl_DStringInit(&result);

    for (i=1; i<argc; i++) {
	if (stat(argv[i], &st) != -1) {
	    if (S_ISDIR(st.st_mode))
		Tcl_DStringAppendElement(&dirs, argv[i]);
	    else
		Tcl_DStringAppendElement(&files, argv[i]);
	}
    }

    Tcl_DStringAppendElement(&result, Tcl_DStringValue(&dirs));
    Tcl_DStringAppendElement(&result, Tcl_DStringValue(&files));

    Tcl_DStringFree(&dirs);
    Tcl_DStringFree(&files);

    Tcl_DStringResult(interp, &result);

    return TCL_OK;
}

int tcl_trace_type(ClientData clientData, Tcl_Interp *interp,
		   int argc,  char **argv) {
    if (argc != 2) {
	Tcl_SetResult(interp,
		      "wrong # args: should be \"trace_type filename\"\n",
		      TCL_STATIC);
	return TCL_ERROR;
    }

    vTcl_SetResult(interp, "%s",
		   trace_type_int2str(determine_trace_type(argv[1])));
    return TCL_OK;
}

int Tk_utils_Misc_Init(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "expandpath", tcl_expandpath,
		      (ClientData) NULL, NULL);

    Tcl_CreateCommand(interp, "create_console", tcl_createconsole,
		      (ClientData) NULL, NULL);

    Tcl_CreateCommand(interp, "mkdir", tcl_mkdir,
		      (ClientData) NULL, NULL);

    Tcl_CreateCommand(interp, "dir_or_file", tcl_dir_or_file,
		      (ClientData)NULL, NULL);

    Tcl_CreateCommand(interp, "trace_type", tcl_trace_type,
		      (ClientData)NULL, NULL);

    /* from capture.c */
    Tcl_CreateCommand(interp, "capture", tcl_capture,
		      (ClientData)NULL, NULL);

    return TCL_OK;
}

