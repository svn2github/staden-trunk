#include <string.h>
#include <ctype.h>

#include <tg_gio.h>
#include <xalloc.h>

#include "gap4_compat.h"
#include "misc.h"
#include "io_utils.h"

#define unimp(m) (printf("%s(): unimplemented function\n", #m), -1)
#define outmoded(m) printf("%s(): outmoded concept\n", #m);

/* ----------------------------------------------------------------------
 * Emulating read access to the old FORTRAN arrays
 * Fortunately these had already been wrapped up in C via macros to look like
 * functions, so we can replace them easily here.
 */

/* Contig length */
int io_clength(GapIO *io, int cnum) {
    contig_t *c;

    c = (contig_t *)cache_search(io, GT_Contig, cnum);
    return contig_get_end(&c) - contig_get_start(&c) + 1;
}

/* Left most reading in a contig - outmoded */
int io_clnbr(GapIO *io, int cnum) {
    int rec;

    contig_iterator *ci = contig_iter_new(io, cnum, 1, CITER_FIRST,
					  CITER_CSTART, CITER_CEND);

    rec = contig_iter_next(io, ci)->rec;
    contig_iter_del(ci);

    return rec;
}

/* Right most reading in a contig - outmoded */
int io_crnbr(GapIO *io, int cnum) {
    int rec;
    contig_iterator *ci = contig_iter_new(io, cnum, 1, CITER_LAST,
					  CITER_CSTART, CITER_CEND);

    rec = contig_iter_prev(io, ci)->rec;
    contig_iter_del(ci);

    return rec;
}

/* Reading length */
int io_length(GapIO *io, int rnum) {
    seq_t *s;

    s = (seq_t *)cache_search(io, GT_Seq, rnum);
    return sequence_get_len(&s);
}

/* Reading position */
int io_relpos(GapIO *io, int rnum) {
    int cnum, pos;

    sequence_get_position(io, rnum, &cnum, &pos);
    return pos;
}

/* Left neighbour */
int io_lnbr(GapIO *io, int rnum) {
    static contig_iterator *ci = NULL;
    static int last_rnum = -1;
    static GapIO *last_io;

    rangec_t *r;

    if (ci == NULL || rnum != last_rnum || io != last_io) {
	int cnum, pos;

	sequence_get_position(io, rnum, &cnum, &pos);

	if (ci)
	    contig_iter_del(ci);
	ci = contig_iter_new(io, cnum, 1, CITER_LAST, CITER_CSTART, pos);
	if (!ci)
	    return 0;

	/* Now get tot he correct point in the iterator: rnum */
	do {
	    r = contig_iter_prev(io, ci);
	} while (r && r->rec != rnum);

	if (!r)
	    return 0;

	last_io = io;
    }

    /*
     * At this stage we have an iterator and the last item used it in was
     * 'rnum', so we can simply call prev.
     */
    r = contig_iter_prev(io, ci);
    if (r) {
	last_rnum = r->rec;
	return r->rec;
    }

    return 0;
}

/* Right neightbour */
int io_rnbr(GapIO *io, int rnum) {
    static contig_iterator *ci = NULL;
    static int last_rnum = -1;
    static GapIO *last_io;
    rangec_t *r;

    if (ci == NULL || rnum != last_rnum || io != last_io) {
	int cnum, pos;
	sequence_get_position(io, rnum, &cnum, &pos);

	if (ci)
	    contig_iter_del(ci);
	ci = contig_iter_new(io, cnum, 1, CITER_LAST, pos, CITER_CEND);
	if (!ci)
	    return 0;

	/* Now get tot he correct point in the iterator: rnum */
	do {
	    r = contig_iter_next(io, ci);
	} while (r && r->rec != rnum);

	if (!r)
	    return 0;

	last_io = io;
    }

    /*
     * At this stage we have an iterator and the last item used it in was
     * 'rnum', so we can simply call prev.
     */
    r = contig_iter_next(io, ci);
    if (r) {
	last_rnum = r->rec;
	return r->rec;
    }

    return 0;
}

/* ---------------------------------------------------------------------- */

int io_write_seq(GapIO *io,	/*  */
		 int N,		/* record/gel number  */
		 int  *length,	/* length of complete string */
		 int  *start,	/* start */
		 int  *end,	/* end */
		 char *seq,	/* complete sequence */
		 int1 *conf,	/* confidence vals */
		 int2 *opos)	/* original pos */
{
    return unimp(io_write_seq);
}

int io_write_rd(GapIO *io,	/*  */
		int N,		/* record/gel number  */
		char *file,	/* trace file */
		int filelen,
		char *type,	/* trace file type */
		int typelen)
{
    return unimp(io_write_rd);
}

int io_write_free_annotation(GapIO *io, int *f)
{
    return unimp(io_write_free_annotation);
}

int io_write_annotation(GapIO *io,
			int N,
			int *anno)
{
    return unimp(io_write_annotation);
}

/*
 * Obtains a contig name (leftmost reading name) from a contig number.
 * The returned value is only valid until the next call of get_contig_name().
 */
char *get_contig_name(GapIO *io, int number) {
    static char name[DB_NAMELEN+1];
    contig_t *c;

    c = (contig_t *)cache_search(io, GT_Contig, number);
    strcpy(name, contig_get_name(&c));

    return name;
}

/*
 * Converts a template name to a template number.
 *
 * Arguments:
 *     io	- GapIO *
 *     tname    - the string described above
 *
 * Returns:
 *    0 for failure, otherwise the template number
 */
int template_name_to_number(GapIO *io, char *tname) {
    return unimp(template_name_to_number);
}

/*
 * Converts a reading name to a reading number.
 *
 * Arguments:
 *     io	- GapIO *
 *     rname    - the string described above
 *
 * Returns:
 *    0 for failure, otherwise the gel number
 */
int read_name_to_number(GapIO *io, char *gel_name) {
    int n = sequence_index_query(io, gel_name);
    return n > 0 ? n : 0;
}

/*
 * Converts a contig name to a contig number.
 * WARNING: slow for now.
 *
 * Returns:
 *     contig record number for success
 *     0 for failure
 */
int contig_name_to_number(GapIO *io, char *name) {
    int n = contig_index_query(io, name);
    return n > 0 ? n : 0;
}

/*
 * Converts a reading identifier string into a gel number.
 * "#[0-9]+" represents a gel number in string form.
 * "=.*"     represents a contig identifier "=[0-9]+"
 * ".*"      represents a gel name, as given to read_name_to_number().
 *
 * Arguments:
 *     io	- GapIO *
 *     gel_name - the string described above
 *     is_name  - whether to assume this is a text name rather than an ID.
 *		  not really needed now as we've removed the / anyway.
 *
 * Returns:
 *    -1 for failure, otherwise the gel number
 */
int get_gel_num(GapIO *io, char *gel_name, int is_name) {
    int n;

    if (*gel_name == '#')
	return atoi(gel_name+1);

    if (*gel_name == '=')
	return io_clnbr(io, atoi(gel_name+1));

    n = read_name_to_number(io, gel_name);
    return n ? n : -1;
}

/*
 * Converts a gel number to a contig number. Left-chaining is done
 * automatically.
 *
 * Arguments:
 *     io	- GapIO *
 *     gel	- Any gel number
 *
 * Returns:
 *    -1 for failure, otherwise the contig number
 */
int rnumtocnum(GapIO *io, int gel) {
    return sequence_get_contig(io, gel);
}

void update_rnumtocnum(GapIO *io, int gel, int contig) {
    unimp(update_rnumtocnum);
}

void invalidate_rnumtocnum(GapIO *io, int disable) {
    unimp(invalidate_rnumtocnum);
}

/*
 * Get contig number - converts a identifier string to a contig number by
 * combining the above routines.
 *
 * Arguments:
 *    (see get_gel_num)
 *
 * Returns:
 *    -1 for failure, otherwise the contig number
 */
int get_contig_num(GapIO *io, char *gel_name, int is_name) {
    int gel;

    /* Look for contig name first */
    gel = contig_name_to_number(io, gel_name);
    if (gel)
	return gel;

    /* Otherwise we assume it is a reading name */
    gel = get_gel_num(io, gel_name, is_name);
    if (gel == -1)
	return -1;
    
    /* And return it's contig number */
    return rnumtocnum(io, gel);
}

f_int nameno_(char *name, f_int *ngels, f_int *handle, char name_l) {
    return unimp(nameno_);
}

/*
 * Chains left along a gel list to find the leftmost gel in a contig.
 * This function is database corrupt resistant.
 *
 * Arguments:
 *     io	- GapIO *
 *     gel	- Any gel number
 *
 * Returns:
 *    -1 for failure, otherwise the leftmost gel number
 */
int chain_left(GapIO *io, int gel) {
    int cnum;
    
    outmoded(chain_left);
    cnum = rnumtocnum(io, gel);
    return io_clnbr(io, cnum);
}

/*
 * Given a list of reading names we return a list of reading numbers. This is
 * done with efficiency in mind - we prefer to do fewer disk accesses at the
 * expense of more memory accesses.
 *
 * Returns: 0 for success
 *         -1 for failure
 * and the number and values of readings in argc and argv.
 */
int lget_gel_num(GapIO *io, int listArgc, char **listArgv, /* INPUT list  */
		 int *rargc, int **rargv) {                /* OUTPUT list */
    int i,j,count=0;

    if (NULL == (*rargv = (int *)xmalloc(listArgc * sizeof(int))))
	return -1;

    /* Translate =cnum and #rnum into numbers */
    for (j=0; j<listArgc; j++) {
	if (listArgv[j][0] == '#') {
	    (*rargv)[j] = atoi(&listArgv[j][1]);
	    count++;
	} else if (listArgv[j][0] == '=') {
	    int num = atoi(&listArgv[j][1]);
	    (*rargv)[j] = num ? io_clnbr(io, num) : 0;
	    count++;
	} else {
	    (*rargv)[j] = 0;
	}
    }

    /* With tcl hash tables - no need to optimise this code anyway */
    for (i = 0; i < listArgc; i++) {
	if ((*rargv)[i] == 0) {
	    int rnum = get_gel_num(io, listArgv[i], GGN_ID);
	    if (rnum != -1) {
		(*rargv)[i] = rnum;
		count++;
	    }
	}
    }

    if (count != listArgc) {
	for (i=j=0; j<listArgc; j++) {
	    if ((*rargv)[j] != 0) {
		(*rargv)[i++] = (*rargv)[j];
	    }
	}
    }

    *rargc = count;
    return 0;
}


/*
 * Given a list of reading names we return a list of contig numbers. This is
 * done with efficiency in mind - we prefer to do fewer disk accesses at the
 * expense of more memory accesses.
 *
 * Returns: 0 for success
 *         -1 for failure
 * and the number and values of contigs in argc and argv.
 */
int lget_contig_num(GapIO *io, int listArgc, char **listArgv, /* INPUT list  */
		    int *rargc, contig_list_t **rargv) {      /* OUTPUT list */
    int i, j, count=0;
    char *p;

    /* allocate: +1 ensures that we never alloc zero bytes */
    if (NULL == (*rargv = (contig_list_t *)xmalloc(1 + listArgc *
						   sizeof(contig_list_t))))
	return -1;

    /* Scan through list finding the reading indentifiers and start/ends */
    for (j = 0; j < listArgc; j++) {
	/* Find end of identifier */
	for (p = listArgv[j]; *p && !isspace(*p); p++)
	    ;

	if (*p) {
	    *p = 0;
	    p++;
	    /* look for start and end */
	    switch(sscanf(p, "%d %d",
			  &(*rargv)[j].start,
			  &(*rargv)[j].end)) {
	    case 1:
		(*rargv)[j].end = 0;
		break;
	    case 0:
	    case -1:
		(*rargv)[j].start = 1;
		(*rargv)[j].end = 0;
	    }
	} else {
	    *p = 0;
	    (*rargv)[j].start = 1;
	    (*rargv)[j].end = 0;
	}

	/* If identifier is #num, then translation is trivial */
    }

    /* Translate #cnum and #rnum into numbers */
    for (j=0; j<listArgc; j++) {
	if (listArgv[j][0] == '#') {
	    int num = atoi(&listArgv[j][1]);
	    if (num > 0) {
		(*rargv)[j].contig = rnumtocnum(io, num);
		count++;
	    } else {
		(*rargv)[j].contig = 0;
	    }
	} else if (listArgv[j][0] == '=') {
	    int num = atoi(&listArgv[j][1]);
	    if (num > 0) {
		(*rargv)[j].contig = num;
		count++;
	    } else {
		(*rargv)[j].contig = 0;
	    }
	} else {
	    (*rargv)[j].contig = 0;
	}
    }

    /* Convert names to numbers */
    /* Translate reading names to reading numbers */
    for (j = 0; j < listArgc; j++) {
	if (0 == (*rargv)[j].contig) {
	    int cnum = contig_name_to_number(io, listArgv[j]);
	    if (cnum) {
		(*rargv)[j].contig = cnum;
		count++;
	    }
	}
    }

    /* Remove duplicates? */
    //fprintf(stderr, "FIXME: remove duplicates in lget_contig_num");
    /* FIXME: to do */

    /*
     * Handle case when we've failed to find some; we just shuffle down to
     * fill any holes in the rargv structure.
     */
    for (i=j=0; j<listArgc; j++) {
	if ((*rargv)[j].contig != 0) {
	    (*rargv)[i++] = (*rargv)[j];
	}
    }


    /* Set ranges */
    for (j=0; j<listArgc; j++) {
	if ((*rargv)[j].end == 0)
	    (*rargv)[j].end = ABS(io_clength(io, (*rargv)[j].contig));
    }

    /* Check for failures */
    for (i=j=0; j<listArgc; j++) {
	if ((*rargv)[j].contig > 0) {
	    (*rargv)[i++] = (*rargv)[j];
	}
    }

    *rargc = i;
    return 0;
}

char *get_read_name(GapIO *io, int number) {
    seq_t *s;

    s = (seq_t *)cache_search(io, GT_Seq, number);
    return sequence_get_name(&s);
}

/* ----------------------------------------------------------------------
 */
GapIO *open_db(char *project, char *version, int *status, int create,
	       int read_only) {
    char name[1024];
    sprintf(name, "%s.%s", project, version);

    *status = 0;
    return gio_open(name, read_only, create);
}

int close_db(GapIO *io) {
    gio_close(io);
    return 0;
}


/* ----------------------------------------------------------------------
 * The old <datatype>_read and <datatype>_write functions.
 */
#if 0
int x_gel_read     (GapIO *io, int n, void *v) { unimp("gel_read"); }
int x_contig_read  (GapIO *io, int n, void *v) { unimp("contig_read"); }
int x_tag_read     (GapIO *io, int n, void *v) { unimp("tag_read"); }
int x_template_read(GapIO *io, int n, void *v) { unimp("template_read"); }
int x_vector_read  (GapIO *io, int n, void *v) { unimp("vector_read"); }
int x_clone_read   (GapIO *io, int n, void *v) { unimp("clone_read"); }

int x_gel_write     (GapIO *io, int n, void *v) { unimp("gel_write"); }
int x_contig_write  (GapIO *io, int n, void *v) { unimp("contig_write"); }
int x_tag_write     (GapIO *io, int n, void *v) { unimp("tag_write"); }
int x_template_write(GapIO *io, int n, void *v) { unimp("template_write"); }
int x_vector_write  (GapIO *io, int n, void *v) { unimp("vector_write"); }
int x_clone_write   (GapIO *io, int n, void *v) { unimp("clone_write"); }
#endif

/* ----------------------------------------------------------------------
 */
int CalcTotalContigLen(GapIO *io) {
    int len = 0, i;
    contig_t *c;

    for (i = 0; i < NumContigs(io); i++) {
	int crec = arr(GCardinal, io->contig_order, i);
	c = (contig_t *)cache_search(io, GT_Contig, crec);
	len += c->end - c->start + 1;
    }

    return len;
}

void bell(void) {
    extern Tcl_Interp *GetInterp(void);

    Tcl_Eval(GetInterp(), "bell");
}

void complement_contig(void) {
    puts("FIXME: complement_contig unimplemented");
}
