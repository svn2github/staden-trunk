#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "list_proc.h"
#include "IO.h"
#include "misc.h"
#include "gap_globals.h"
#include "fort.h"
#include "complement.h"

/*
 * Loop through the input files finding the longest one.
 * This is a worst-case scenario for the memory allocation in assembly
 * (provided they're not compressed).
 */
static int largest_file(void) {
    char *item;
    struct stat statbuf;
    int largest = 20000;
    /* default needed when reading via RAWDATA so no local files found */

    while (item = get_active_list_item()) {
	if (stat(item, &statbuf) == -1) {
	    /* perror(item); */
	    continue;
	}
	if (largest < statbuf.st_size)
	    largest = statbuf.st_size;
    }

    rewind_active_list();

    return largest;
}

/*
 * This checks for mutation-study type contigs where we expect to get all the
 * forward sequences in one orientation and all the reverse sequences in
 * another orientation. In this case there's a 50-50 chance of having the
 * wrong 'logical' orientation (reverses should be complemented). This
 * function will correct this.
 *
 * In order to not complement existing contigs we need a way of determining
 * which contigs are new. The number of contigs could shrink, but still
 * containing new contigs (if joins were made), so we define new contigs as
 * those consisting entirely of new readings. For this we need to know the
 * number of readings in the database prior to assembly.
 *
 * Arguments:
 *	io			GapIO
 *	prev_nreadings		Number of readings prior to assembly.
 */
static void check_contig_orientations(GapIO *io, int prev_nreadings) {
    int nc = NumContigs(io);
    int cnum, rnum;
    int pos_count, neg_count, other_count;

    for (cnum = 1; cnum <= nc; cnum++) {
	/* Is this a new contig? */
	for (rnum = io_clnbr(io, cnum); rnum; rnum = io_rnbr(io, rnum)) {
	    if (rnum < prev_nreadings)
		break;
	}
	if (rnum)
	    continue;

	/* Yes, count number of pos and neg orientation contigs */
	pos_count = neg_count = other_count = 0;
	for (rnum = io_clnbr(io, cnum); rnum; rnum = io_rnbr(io, rnum)) {
	    GReadings r;
	    gel_read(io, rnum, r);

	    switch (STRAND(r)) {
	    case GAP_STRAND_FORWARD:
		if (r.sense == GAP_SENSE_ORIGINAL)
		    pos_count++;
		else
		    neg_count++;
		break;

	    case GAP_STRAND_REVERSE:
		if (r.sense == GAP_SENSE_REVERSE)
		    pos_count++;
		else
		    neg_count++;
		break;

	    default:
		other_count++;
	    }
	}

	/* Complement it if neg_count is 90% or more of the total */
	if (pos_count + neg_count >= 2 &&
	    neg_count / (pos_count + neg_count + other_count) >= 0.9) {
	    complement_contig(io, cnum);
	}
    }
}

char *
auto_assemble(f_int handle,                                            /* in */
	      char *inlist,                                            /* in */
              f_int iopt,                                              /* in */
	      f_int entry,                                             /* in */
	      f_int disp_mode,                                         /* in */
	      f_int min_mat,                                           /* in */
	      f_int min_ovr,                                           /* in */
	      f_int max_pad,                                           /* in */
	      float max_mis,                                           /* in */
	      f_int align,                                             /* in */
	      int joins,                                               /* in */
	      f_int fail_mode,                                         /* in */
	      f_int win_size,                                          /* in */
	      f_int max_dash,                                          /* in */
	      f_int ignore_prev,                                       /* in */
	      float percd)                                             /* in */
{
    /* dbas.f dbauto */
    GapIO *io;
    f_int maxcon;
    f_int maxsav;
    f_int maxglm;
    char *seq1;
    char *seq2;
    char *seq3;
    char *seq4;
    char *seq5;
    char *seqc2;
    char *seqc3;
    char *seqg2;
    char *seqg3;
    char *rnames;
    f_int *sav1;
    f_int *sav2;
    f_int *sav3;
    f_int *cends;
    f_int *nends;
    f_int iok;
    char namarc[F_NAMLEN];
    char *dbname;
    f_int nconts, ngels;
    f_int *clist = NULL;
    void *dl;
    char *ret;
    int max_gel_len;
    int prev_nreadings;

    if ( (io = io_handle(&handle)) == NULL){	
	verror(ERR_FATAL, "auto_assemble", "invalid io handle");
	return NULL;
    }

    invalidate_rnumtocnum(io, 1);

    prev_nreadings = NumReadings(io);

    /*
     * NB use a local copy of ngels and nconts here - DO NOT pass
     *  NumReadings(io) & NumContigs(io) 
     */
    if (-1 == set_active_list(inlist)) {
	return NULL;
    }

    /*
     * The maximum space to allocation for a single reading is the MAX
     * of the current reading length in the database plus the maximum
     * file length, with an extra amount to allow for padding.
     * (This is just a heuristic hack as a workaround for the fixed memory
     * model of fortran. It will go away when we rewrite the assembly
     * engine.)
     */
    {
	int l1, l2;
	l1 = largest_file();
	l2 = find_max_gel_len(io, 0, 0);
	
	max_gel_len = 100 + 1.5 * MAX(l1, l2);
    }

    dbname = io_name(io);

    maxcon = io_dbsize(io);
    maxglm = 2 * max_gel_len;
    /* maxsav = io_dbsize(io); */
    maxsav = maxglm;

    nconts = NumContigs(io);
    ngels = NumReadings(io);

    if ((cends = (f_int *)xmalloc(maxcon * sizeof(f_int)))==NULL){
	return NULL;
    }
    if ((nends = (f_int *)xmalloc(maxcon * sizeof(f_int)))==NULL){
	return NULL;
    }
    if ((seq1 = (char *)xmalloc(maxseq * sizeof(char)))==NULL){
	return NULL;
    }
    if ((seq2 = (char *)xmalloc(maxglm * sizeof(char)))==NULL){
	return NULL;
    }
    if ((seq3 = (char *)xmalloc(maxglm * sizeof(char)))==NULL){
	return NULL;
    }
    if ((seq4 = (char *)xmalloc(maxglm * sizeof(char)))==NULL){
	return NULL;
    }
    if ((seq5 = (char *)xmalloc(maxglm * sizeof(char)))==NULL){
	return NULL;
    }
    if ((seqc2 = (char *)xmalloc(maxglm * 2 * sizeof(char)))==NULL){
	return NULL;
    }
    if ((seqg2 = (char *)xmalloc(maxglm * 2 * sizeof(char)))==NULL){
	return NULL;
    }
    if ((seqg3 = (char *)xmalloc(maxglm * sizeof(char)))==NULL){
	return NULL;
    }
    if ((seqc3 = (char *)xmalloc(maxglm * sizeof(char)))==NULL){
	return NULL;
    }
    if ((rnames = (char *)xmalloc(io_maxdbsize(io) * DB_NAMELEN *
				   sizeof(char)))
        ==NULL){
	return NULL;
    }
    if ((sav1 = (f_int *)xmalloc(maxsav * sizeof(f_int)))==NULL){
	return NULL;
    }
    if ((sav2 = (f_int *)xmalloc(maxsav * sizeof(f_int)))==NULL){
	return NULL;
    }
    if ((sav3 = (f_int *)xmalloc(maxsav * sizeof(f_int)))==NULL){
	return NULL;
    }


    entry = 1 - entry;
    joins = 1 - joins;
    tolist_(NULL, NULL, 1, 0);
    dbauto_(&io_relpos(io,1), &io_length(io,1), &io_lnbr(io,1),
             &io_rnbr(io,1), &io_maxdbsize(io), &io_dbsize(io),
             &ngels, &nconts, clist, &max_gel_len,
             seq1, seq2, seq3, seq4, seq5, seqc2, seqg2, seqg3, seqc3, rnames,
             &maxseq, &maxglm,
             sav1, sav2, sav3, &maxsav, cends, nends, &maxcon,
             &handle, namarc, dbname, &percd, 
             &entry, &disp_mode, &min_mat, &max_pad, &max_mis, &align, 
             &iopt, &ignore_prev, &joins, &fail_mode, &win_size, 
	     &max_dash, "dummy", &iok, &min_ovr,
             maxseq, maxglm, maxglm, maxglm, maxglm, 2*maxglm, 2*maxglm, 
             maxglm, maxglm, DB_NAMELEN,
             F_NAMLEN, DB_FILELEN, strlen("dummy")); 
    flush2t(io);

    dl = tolist_(NULL, NULL, 0, 0);
    if (dl) {
	ret = strdup(read_dlist(dl));
	free_dlist(dl);
    } else {
	ret = strdup("");
    }

    check_contig_orientations(io, prev_nreadings);

    invalidate_rnumtocnum(io, 0);

    xfree(cends);
    xfree(nends);
    xfree(seq1);
    xfree(seq2);
    xfree(seq3);
    xfree(seq4);
    xfree(seq5);
    xfree(seqc2);
    xfree(seqg2);
    xfree(seqc3);
    xfree(seqg3);
    xfree(rnames);
    xfree(sav1);
    xfree(sav2);
    xfree(sav3);

    return ret;

} /* end auto_assemble */


