#ifndef _COMPLEMENT_CONTIG_H_
#define _COMPLEMENT_CONTIG_H_

/*
 * Complements a contig. A C interface to the Fortran CMPLMT function.
 *
 * Returns 0 for success, -1 for failure
 */
int complement_contig(GapIO *io, int contig);

int tk_complement_contig(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]);

#endif
