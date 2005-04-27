#include <tk.h>

#include "IO.h"
#include "io_handle.h"
#include "complement.h"
#include "gap_cli_arg.h"
#include "io-reg.h"
#include "misc.h"
#include "tcl_utils.h"
#include "fort.h"
#include "tagUtils.h"
#include "list_proc.h"
#include "dna_utils.h"

typedef struct {
    int position;
    int rnum;
} read_details;

static int rd_sort_func(const void *p1, const void *p2) {
    return ((const read_details *)p1)->position
	- ((const read_details *)p2)->position;
}

/*
 * Complements a contig. This is a C rewrite of the old CMPLMT fortran
 * function.
 *
 * Initially we resort the sequence coordinates so that they are sorted
 * decrementally by their end coordinates. Then we effectively complement
 * those coordinates. This gets written back to the reading
 * structures/contig left/right ends.
 *
 * Secondly we reverse/complement the sequences themselves.
 *
 * Finally we update the tags on the consensus.
 *
 * Returns 0 for success, -1 for failure
 */
int complement_contig(GapIO *io, int contig) {
    read_details *rd = NULL;
    int rnum, clen = io_clength(io, contig);
    int count, i, err = 0;
    GContigs c;
    reg_complement rc;

    if (contig_lock_write(io, contig) == -1) {
	verror(ERR_WARN, "complement_contig", "Contig is busy");
	return -1;
    }

    /* Fill our read_details array with sequence numbers and lengths */
    for (count = 0, rnum =  io_clnbr(io, contig);
	 rnum;
	 rnum = io_rnbr(io, rnum))
	count++;
    if (NULL == (rd = (read_details *)malloc(sizeof(*rd) * count)))
	return -1;

    for (count = 0, rnum =  io_clnbr(io, contig);
	 rnum;
	 rnum = io_rnbr(io, rnum), count++) {
	rd[count].rnum = rnum;
	rd[count].position = io_relpos(io, rnum)+ABS(io_length(io, rnum))-1;
    }

    /*
     * Sort incrementally.
     * This list is mostly sorted already, so maybe a bubble sort
     * would be preferable.
     */
    qsort(rd, count, sizeof(*rd), rd_sort_func);

    /*
     * Reproduce the contig linked list from the sorted array.
     */
    io_crnbr(io, contig) = rd[0].rnum;
    for (i = 0; i < count; i++) {
	io_rnbr(io, rd[i].rnum)   = i > 0 ? rd[i-1].rnum : 0;
	io_lnbr(io, rd[i].rnum)   = i < count-1 ? rd[i+1].rnum : 0;
	io_relpos(io, rd[i].rnum) = clen+1-rd[i].position;
	io_length(io, rd[i].rnum) = -io_length(io, rd[i].rnum);
    }
    io_clnbr(io, contig) = rd[count-1].rnum;

    /*
     * Update the GReadings and GContigs structures.
     * The r.start, r.end and r.length will be updated later by the
     * call to io_complement_seq().
     */
    for (i = 0; i < count; i++) {
	GReadings r;
	rnum = rd[i].rnum;
	gel_read(io, rnum, r);
	r.left     = io_lnbr(io, rnum);
	r.right    = io_rnbr(io, rnum);
	r.position = io_relpos(io, rnum);
	r.sense   ^= 1;
	gel_write(io, rnum, r);
    }
    contig_read(io, contig, c);
    c.left = io_clnbr(io, contig);
    c.right = io_crnbr(io, contig);
    contig_write(io, contig, c);

    /* Complement the sequence data */
    for (i = 0; i < count; i++) {
	int length, start, end;
	char *seq = NULL;
	int1 *conf = NULL;
	int2 *opos = NULL;
	if (0 != io_aread_seq(io, rd[i].rnum,
			      &length, &start, &end,
			      &seq, &conf, &opos, 0)) {
	    err = 1;
	    continue;
	}
	io_complement_seq(&length, &start, &end, seq, conf, opos);
	io_write_seq(io, rd[i].rnum, &length, &start, &end, seq, conf, opos);
	if (seq)  xfree(seq);
	if (conf) xfree(conf);
	if (opos) xfree(opos);
    }
    xfree(rd);

    /* Complement the contig tags */
    complement_contig_tags(io, contig);
    
    flush2t(io);

    /*
     * Notify complement - this should be done in CMPLMT really, as should
     * the busy checks etc.
     */
    rc.job = REG_COMPLEMENT;

    contig_notify(io, contig, (reg_data *)&rc);

    return err;
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

f_proc_ret cmplmt_(f_int *RELPG, f_int *LNGTHG, f_int *LNBR,   f_int *RNBR,
		   f_int *NGELS, f_int *NCONTS, f_int *LINCON, f_int *LLINO,
		   char  *GEL,   f_int *IDBSIZ, f_int *IDEVR,  f_int *MAXGEL,
		   f_implicit gel_l) {
    GapIO *io;
    int contig;

    if ( (io = io_handle(IDEVR)) == NULL) f_proc_return();
    contig = io_dbsize(io) - *LINCON;
    complement_contig(io, contig);
    f_proc_return();
}

f_proc_ret sqcom_(char *SEQ, f_int *IDIM, f_implicit gel_l) {
    complement_dna(SEQ, *IDIM);
    f_proc_return();
}

f_proc_ret sqrev_(char *SEQ, f_int *IDIM, f_implicit gel_l) {
    reverse_dna(SEQ, *IDIM);
    f_proc_return();
}
