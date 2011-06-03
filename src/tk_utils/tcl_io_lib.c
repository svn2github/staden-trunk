#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "tcl_io_lib.h"

/*
 * A tcl interface to the staden-io_lib package, for reading sequences
 * from trace files or URLs to remote trace servers.
 */
int tcl_read_seq_trace(ClientData clientData,
		       Tcl_Interp *interp,
		       int objc,
		       Tcl_Obj *CONST objv[]) {
    Tcl_Obj *l, *o;
    int i, j;
    static char sanitized[256];
    static int init = 0;

    /* Initialise sanitization array. BWA can't handle non-ACGT for example */
    if (!init) {
	memset(sanitized, 'N', 256);
	for (i = 0; i < 8; i++) {
	    sanitized["ACGTacgt"[i]] = "ACGTACGT"[i];
	}
	init = 1;
    }

    l = Tcl_NewListObj(0, NULL);
    for (i = 1; i < objc; i++) {
	char *n = Tcl_GetString(objv[i]);
	Read *r;
	char *qual;

	r = read_reading(n, TT_ANY);
	if (!r) {
	    Tcl_SetResult(interp, "Failed to read trace\n", TCL_STATIC);
	    return TCL_ERROR;
	}


	/* Seq */
	for (j = 0; j < r->NBases; j++) {
	    r->base[j] = sanitized[r->base[j]];
	}
	o = Tcl_NewStringObj(r->base, r->NBases);
	Tcl_ListObjAppendElement(interp, l, o);


	/* Qual */
	if (NULL == (qual = malloc(r->NBases)))
	    return TCL_ERROR;
	
	for (j = 0; j < r->NBases; j++) {
	    switch (r->base[j]) {
	    case 'a':
	    case 'A':
		qual[j] = r->prob_A[j] + '!'; /* fastq format */
		break;

	    case 'c':
	    case 'C':
		qual[j] = r->prob_C[j] + '!';
		break;

	    case 'g':
	    case 'G':
		qual[j] = r->prob_G[j] + '!';
		break;

	    case 't':
	    case 'T':
	    case 'u':
	    case 'U':
		qual[j] = r->prob_T[j] + '!';
		break;

	    default:
		qual[j] = '!';
		break;
	    }
	}

	o = Tcl_NewStringObj(qual, r->NBases);
	Tcl_ListObjAppendElement(interp, l, o);
	
	free(qual);
	read_deallocate(r);
    }

    Tcl_SetObjResult(interp, l);
    return TCL_OK;
}
