#include <tk.h>

#include "IO.h"
#include "io_handle.h"
#include "complement.h"
#include "gap_cli_arg.h"
#include "io-reg.h"
#include "misc.h"
#include "tcl_utils.h"
#include "fort.h"
#include "list_proc.h"

/*
 * Complements a contig. A C interface to the Fortran CMPLMT function.
 *
 * Returns 0 for success, -1 for failure
 */
int complement_contig(GapIO *io, int contig) {
    f_int *handle;
    f_int llino;
    f_int lincon;
    char *gel;
    reg_complement rc;
    int max_len = find_max_gel_len(io, contig, 0)+1;

    if (contig_lock_write(io, contig) == -1) {
	verror(ERR_WARN, "complement_contig", "Contig is busy");
	return -1;
    }

    if (NULL == (gel = (char *)xmalloc(max_len)))
	return -1;

    handle = handle_io(io);
    lincon = io_dbsize(io) - contig;
    llino  = io_clnbr(io, contig);

    /*
     * SUBROUTINE CMPLMT(RELPG,LNGTHG,LNBR,RNBR,NGELS,NCONTS,
     * +LINCON,LLINO,GEL,IDBSIZ,IDEVR,MAXGEL)
     * INTEGER RELPG(IDBSIZ)
     * INTEGER LNGTHG(IDBSIZ),LNBR(IDBSIZ),RNBR(IDBSIZ)
     * CHARACTER GEL(MAXGEL)
     */
    cmplmt_(&io_relpos(io,1), &io_length(io,1), &io_lnbr(io,1), &io_rnbr(io,1),
	    &NumReadings(io), &NumContigs(io), &lincon, &llino, gel,
	    &io_dbsize(io), handle, &max_len,
	    /* implicit args */ max_len);

    xfree(gel);

    /*
     * Notify complement - this should be done in CMPLMT really, as should
     * the busy checks etc.
     */
    rc.job = REG_COMPLEMENT;

    contig_notify(io, contig, (reg_data *)&rc);
    
    flush2t(io);

    return 0;
}


/*
 * Tk interface to complement_contig
 */
typedef struct {
    GapIO *io;
    char *contigs;
} com_args;

int tk_complement_contig(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) {
    com_args args;
    int num_contigs;
    contig_list_t *contig_array = NULL;
    int *contigs;
    int err = 0, i;

    cli_args a[] = {
	{"-io",      ARG_IO, 1, NULL, offsetof(com_args, io)},
	{"-contigs", ARG_STR,1, NULL, offsetof(com_args, contigs)},
	{NULL,      0,      0, NULL, 0}
    };

    vfuncheader("complement contig");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /* create contig name array */
    active_list_contigs(args.io, args.contigs, &num_contigs, &contig_array);
    if (num_contigs == 0) {
	xfree(contig_array);
	return TCL_OK;
    }
    contigs = to_contigs_only(num_contigs, contig_array);
    xfree(contig_array);

    for (i = 0; i < num_contigs; i++)
	if (-1 == complement_contig(args.io, contigs[i]))
	    err |= 1;

    xfree(contigs);

    if (err)
	Tcl_SetResult(interp, "1", TCL_STATIC);
    else 
	Tcl_SetResult(interp, "0", TCL_STATIC);

    return TCL_OK;
}
