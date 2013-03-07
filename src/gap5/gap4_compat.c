#include <string.h>
#include <ctype.h>

#include <tg_gio.h>
#include <xalloc.h>

#include "gap4_compat.h"
#include "misc.h"
#include "io_utils.h"
#include "consensus.h"
#include "tg_scaffold.h"
#include "io_lib/hash_table.h"

#define unimp(m) (printf("%s(): unimplemented function\n", #m), -1)
#define outmoded(m) printf("%s(): outmoded concept\n", #m);

/* ----------------------------------------------------------------------
 * Emulating read access to the old FORTRAN arrays
 * Fortunately these had already been wrapped up in C via macros to look like
 * functions, so we can replace them easily here.
 */

/*
 * Contig length, but not quite in the tradtional 1..length sense of Gap4.
 * Rather this is end-start+1 as start can be anywhere.
 */
int io_clength(GapIO *io, tg_rec cnum) {
    contig_t *c;

    c = (contig_t *)cache_search(io, GT_Contig, cnum);
    if (!c) {
	verror(ERR_FATAL, "io_clength()",
	       "Failed to load contig #%"PRIrec" in io_clength()",
	       cnum);
	return 0;
    }
    return contig_get_end(&c) - contig_get_start(&c) + 1;
}

/* Clipped contig length */
int io_cclength(GapIO *io, tg_rec cnum) {
    int start, end;
    if (-1 == consensus_valid_range(io, cnum, &start, &end)) {
	verror(ERR_FATAL, "io_cclength()",
	       "Failed to load contig #%"PRIrec" in io_clength()",
	       cnum);
	return 0;
    }

    return end - start + 1;
}

/* Left most reading in a contig - outmoded */
tg_rec io_clnbr(GapIO *io, tg_rec cnum) {
    tg_rec rec;
    rangec_t *r;
    contig_iterator *ci = contig_iter_new(io, cnum, 1, CITER_FIRST,
					  CITER_CSTART, CITER_CEND);
    if (!ci)
	return 0;

    r = contig_iter_next(io, ci);
    rec = r ? r->rec : 0;
    contig_iter_del(ci);

    return rec;
}

/* Right most reading in a contig - outmoded */
tg_rec io_crnbr(GapIO *io, tg_rec cnum) {
    tg_rec rec;
    rangec_t *r;
    contig_iterator *ci = contig_iter_new(io, cnum, 1, CITER_LAST,
					  CITER_CSTART, CITER_CEND);
    if (!ci)
	return 0;

    r = contig_iter_next(io, ci);
    rec = r ? r->rec : 0;
    contig_iter_del(ci);

    return rec;
}

/* Reading length */
int io_length(GapIO *io, tg_rec rnum) {
    seq_t *s;

    s = (seq_t *)cache_search(io, GT_Seq, rnum);
    return sequence_get_len(&s);
}

/* Reading position */
int io_relpos(GapIO *io, tg_rec rnum) {
    tg_rec cnum;
    int pos;

    sequence_get_position(io, rnum, &cnum, &pos, NULL, NULL);
    return pos;
}

/* Left neighbour */
tg_rec io_lnbr(GapIO *io, tg_rec rnum) {
    static contig_iterator *ci = NULL;
    static tg_rec last_rnum = -1;
    static GapIO *last_io;

    rangec_t *r;

    if (ci == NULL || rnum != last_rnum || io != last_io) {
	tg_rec cnum;
	int pos;

	sequence_get_position(io, rnum, &cnum, &pos, NULL, NULL);

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
tg_rec io_rnbr(GapIO *io, tg_rec rnum) {
    static contig_iterator *ci = NULL;
    static tg_rec last_rnum = -1;
    static GapIO *last_io;
    rangec_t *r;

    if (ci == NULL || rnum != last_rnum || io != last_io) {
	tg_rec cnum;
	int pos;
	sequence_get_position(io, rnum, &cnum, &pos, NULL, NULL);

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

/*
 * Obtains a contig name (leftmost reading name) from a contig number.
 * The returned value is only valid until the next call of get_contig_name().
 */
char *get_contig_name(GapIO *io, tg_rec number) {
    static char name[DB_NAMELEN+1];
    contig_t *c;

    c = (contig_t *)cache_search(io, GT_Contig, number);
    if (c) {
	strncpy(name, contig_get_name(&c), DB_NAMELEN);
	name[DB_NAMELEN] = 0;
    } else {
	strcpy(name, "(unknown contig)");
    }

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
tg_rec template_name_to_number(GapIO *io, char *tname) {
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
tg_rec read_name_to_number(GapIO *io, char *gel_name) {
    tg_rec n = sequence_index_query(io, gel_name);
    return n > 0 ? n : 0;
}

/*
 * Converts a scaffold name to a reading number.
 *
 * Arguments:
 *     io	- GapIO *
 *     rname    - the string described above
 *
 * Returns:
 *    0 for failure, otherwise the scaffold record number
 */
tg_rec scaffold_name_to_number(GapIO *io, char *name) {
    tg_rec n = 0;

    /* Check numeric values first */
    if (*name == '=' || *name == '#') {
	n = atorec(name+1);

	if (cache_exists(io, GT_Scaffold, n)) {
	    return n;
	} else if (cache_exists(io, GT_Contig, n)) {
	    contig_t *c = cache_search(io, GT_Contig, n);
	    return c->scaffold;
	}
    }

    /* Also check the name index */
    n = scaffold_index_query(io, name);
    return n > 0 ? n : 0;
}

/*
 * Converts a contig name to a contig number.
 * Name can be =contig_rec, #contig_rec, #reading_rec, contig_name
 * or (slower) reading_name.
 *
 * Returns:
 *     contig record number for success
 *     0 for failure
 */
tg_rec contig_name_to_number(GapIO *io, char *name) {
    tg_rec n = 0;

    /* Check numeric values first */
    if (*name == '=' || *name == '#') {
	n = atorec(name+1);

	if (cache_exists(io, GT_Contig, n)) {
	    return n;
	} else if (cache_exists(io, GT_Seq, n)) {
	    if ((n = rnumtocnum(io, n)) > 0)
		return n;
	}
    }

    /* Also check contig names and sequence names */
    if ((n = contig_index_query(io, name)) > 0)
	return n;

    if ((n = read_name_to_number(io, name)) > 0)
	n = rnumtocnum(io, n);

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
tg_rec get_gel_num(GapIO *io, char *gel_name, int is_name) {
    tg_rec n;

    if (*gel_name == '#')
	return atorec(gel_name+1);

    if (*gel_name == '=')
	return io_clnbr(io, atorec(gel_name+1));

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
tg_rec rnumtocnum(GapIO *io, tg_rec gel) {
    if (cache_exists(io, GT_Seq, gel))
	return sequence_get_contig(io, gel);
    else if (cache_exists(io, GT_Contig, gel))
	return gel;
    else
	return -1;
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
tg_rec get_contig_num(GapIO *io, char *gel_name, int is_name) {
    tg_rec gel;

    /* Look for contig name first */
    gel = contig_name_to_number(io, gel_name);
    if (gel)
	return gel;

    /* Otherwise we assume it is a reading name */
    gel = get_gel_num(io, gel_name, is_name);
    if (gel <= 0)
	return -1;
    
    /* And return it's contig number */
    return rnumtocnum(io, gel);
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
tg_rec chain_left(GapIO *io, tg_rec gel) {
    tg_rec cnum;
    
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
		 int *rargc, tg_rec **rargv) {             /* OUTPUT list */
    int i,j,count=0;

    if (NULL == (*rargv = (tg_rec *)xmalloc(listArgc * sizeof(tg_rec))))
	return -1;

    /* Translate =cnum and #rnum into numbers */
    for (j=0; j<listArgc; j++) {
	if (listArgv[j][0] == '#') {
	    (*rargv)[j] = atorec(&listArgv[j][1]);
	    count++;
	} else if (listArgv[j][0] == '=') {
	    tg_rec num = atorec(&listArgv[j][1]);
	    (*rargv)[j] = num ? io_clnbr(io, num) : 0;
	    count++;
	} else {
	    (*rargv)[j] = 0;
	}
    }

    /* With tcl hash tables - no need to optimise this code anyway */
    for (i = 0; i < listArgc; i++) {
	if ((*rargv)[i] == 0) {
	    tg_rec rnum = get_gel_num(io, listArgv[i], GGN_ID);
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
static int lget_contig_num_base(GapIO *io,
				int listArgc, char **listArgv,       /* IN  */
				int *rargc, contig_list_t **rargv) { /* OUT */
    int i, j, count=0;
    HashTable *h = HashTableCreate(1024, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);

    if (NULL == h) return -1;

    /* allocate: +1 ensures that we never alloc zero bytes */
    if (NULL == (*rargv = (contig_list_t *)xmalloc(1 + listArgc *
						   sizeof(contig_list_t)))) {
	HashTableDestroy(h, 0);
	return -1;
    }

    /* Scan through list finding the reading indentifiers and start/ends */
    for (i = j = 0; j < listArgc; j++) {
	/* See if list item can be split into {contig start end} */
	int iargc;
	char **iargv = NULL;
	char *ctg_name;
	tg_rec num, rnum;
	HashData hd;
	int new;
	contig_t *c;

	if (TCL_OK == Tcl_SplitList(NULL, listArgv[j], &iargc, &iargv)
	    && iargc > 0) {
	    ctg_name = iargv[0];
	    if (!(iargc > 2 && 1 == sscanf(iargv[2], "%d", &(*rargv)[i].end))) {
		(*rargv)[i].end = INT_MAX;
	    }
	    if (!(iargc > 1 && 1 == sscanf(iargv[1], "%d", &(*rargv)[i].start))) {
		(*rargv)[i].start = INT_MAX;
	    }
	} else {
	    ctg_name = listArgv[j];
	    (*rargv)[i].start = INT_MAX;
	    (*rargv)[i].end = INT_MAX;
	}
    
	/* Translate =cnum, #rnum and contig names into numbers */
	if (ctg_name[0] == '#'
	    && (rnum = atorec(&ctg_name[1])) > 0
	    && (num = rnumtocnum(io, rnum)) > 0) {
	    (*rargv)[i].contig = num;
	    count++;
	} else if (ctg_name[0] == '='
		   && (num = atorec(&ctg_name[1])) > 0
		   && (c = cache_search(io, GT_Contig, num)) != NULL) {
	    (*rargv)[i].contig = num;
	    count++;
	} else if ((num = contig_name_to_number(io, ctg_name)) > 0) {
	    (*rargv)[i].contig = num;
	    count++;
	} else {
	    verror(ERR_WARN, "lget_contig_num_base()",
		   "Unknown contig name %s", ctg_name);
	    (*rargv)[i].contig = 0;
	}
	if (iargv) Tcl_Free((char *) iargv);

	/* Skip over failures */
	if ((*rargv)[i].contig <= 0)
	    continue;

	/* Check for duplicates */
	hd.i = 1;
	if (NULL == HashTableAdd(h, (char *)&(*rargv)[i].contig,
				 sizeof(tg_rec), hd, &new)) {
	    /* Out of memory... */
	    xfree(*rargv);
	    HashTableDestroy(h, 0);
	    return -1;
	}
	if (new) i++;  /* Not seen already, so keep it */
    }

    HashTableDestroy(h, 0);

    *rargc = i;
    return 0;
}

/* Full extent of contigs including cutoff data */
int lget_contig_num2(GapIO *io, int listArgc, char **listArgv,
		     int *rargc, contig_list_t **rargv) {
    int ret, i;

    ret = lget_contig_num_base(io, listArgc, listArgv, rargc, rargv);
    if (ret != 0)
	return ret;

    /* Set ranges */
    for (i=0; i<*rargc; i++) {
	contig_t *c;
	int st, en;

	c = cache_search(io, GT_Contig, (*rargv)[i].contig);
	st = c->start;
	en = c->end;

	if ((*rargv)[i].start == INT_MAX || (*rargv)[i].start < st)
	    (*rargv)[i].start = st;

	if ((*rargv)[i].end == INT_MAX || (*rargv)[i].end > en)
	    (*rargv)[i].end = en;

	if ((*rargv)[i].start > en)
	    (*rargv)[i].start = en;

	if ((*rargv)[i].end < st)
	    (*rargv)[i].end = st;
    }

    return 0;
}

/* As above, but clipped to only the used portion of a contig */
int lget_contig_num(GapIO *io, int listArgc, char **listArgv, /* INPUT list  */
		    int *rargc, contig_list_t **rargv) {      /* OUTPUT list */
    int ret, i;

    ret = lget_contig_num_base(io, listArgc, listArgv, rargc, rargv);
    if (ret != 0)
	return ret;

    /* Set ranges */
    for (i=0; i<*rargc; i++) {
	int st, en;

	consensus_valid_range(io, (*rargv)[i].contig, &st, &en);

	if ((*rargv)[i].start == INT_MAX || (*rargv)[i].start < st)
	    (*rargv)[i].start = st;

	if ((*rargv)[i].end == INT_MAX || (*rargv)[i].end > en)
	    (*rargv)[i].end = en;

	if ((*rargv)[i].start > en)
	    (*rargv)[i].start = en;

	if ((*rargv)[i].end < st)
	    (*rargv)[i].end = st;
    }

    return 0;
}

/*
 * As per lget_contig_num_base, but with scaffold IDs instead.
 * Returns 0 on success
 *        -1 on failure
 */
int lget_scaffold_num(GapIO *io, int listArgc, char **listArgv, /* IN  list  */
		      int *rargc, tg_rec **rargv) {             /* OUT list */
    int i, j, count=0;
    char *p;
    HashTable *h;

    /* allocate: +1 ensures that we never alloc zero bytes */
    if (NULL == (*rargv = (tg_rec *)xmalloc(1 + listArgc *
					    sizeof(tg_rec))))
	return -1;

    /* Scan through list finding the scaffold indentifiers */
    for (j = 0; j < listArgc; j++) {
	/* Find end of identifier */
	for (p = listArgv[j]; *p && !isspace(*p); p++)
	    ;

	*p = 0;

	/* If identifier is #num, then translation is trivial */
    }

    /* Translate #cnum and #rnum into numbers */
    for (j=0; j<listArgc; j++) {
	if (listArgv[j][0] == '#' || listArgv[j][0] == '=') {
	    tg_rec num = atorec(&listArgv[j][1]);
	    if (num > 0) {
		(*rargv)[j] = num;
		count++;
	    } else {
		(*rargv)[j] = 0;
	    }
	} else {
	    (*rargv)[j] = 0;
	}
    }

    /* Convert names to numbers */
    for (j = 0; j < listArgc; j++) {
	if (0 == (*rargv)[j]) {
	    tg_rec srec = scaffold_index_query(io, listArgv[j]);
	    if (srec) {
		(*rargv)[j] = srec;
		count++;
	    } else {
		verror(ERR_WARN, "scaffold_index_query()",
		       "Unknown scaffold name %s", listArgv[j]);
	    }
	}
    }


    /* Remove duplicates */
    h = HashTableCreate(1024, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);
    for (i = j = 0; j < listArgc; j++) {
	HashData hd;
	int new;

	if ((*rargv)[j] == 0)
	    continue;

	hd.i = 1;
	HashTableAdd(h, (char *)&(*rargv)[j], sizeof(tg_rec),
		     hd, &new);

	if (new) {
	    (*rargv)[i++] = (*rargv)[j];
	}
    }
    HashTableDestroy(h, 0);
    listArgc = i;

    /*
     * Handle case when we've failed to find some; we just shuffle down to
     * fill any holes in the rargv structure.
     */
    for (i=j=0; j<listArgc; j++) {
	if ((*rargv)[j] != 0) {
	    (*rargv)[i++] = (*rargv)[j];
	}
    }


    /* Check for failures */
    for (i=j=0; j<listArgc; j++) {
	if ((*rargv)[j] > 0) {
	    (*rargv)[i++] = (*rargv)[j];
	}
    }

    *rargc = i;
    return 0;
}

char *get_read_name(GapIO *io, tg_rec number) {
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
#if 0
/* Unclipped length */
int64_t CalcTotalContigLen(GapIO *io) {
    int64_t len = 0;
    int i;
    contig_t *c;

    for (i = 0; i < NumContigs(io); i++) {
	tg_rec crec = arr(tg_rec, io->contig_order, i);
	c = (contig_t *)cache_search(io, GT_Contig, crec);
	len += c ? c->end - c->start + 1 : 0;
    }

    return len;
}
#else
/* Clipped length */
int64_t CalcTotalContigLen(GapIO *io) {
    int64_t len = 0;
    int i;

    for (i = 0; i < NumContigs(io); i++) {
	tg_rec crec = arr(tg_rec, io->contig_order, i);
	len += io_cclength(io, crec);
    }

    return len;
}
#endif

void bell(void) {
    extern Tcl_Interp *GetInterp(void);

    Tcl_Eval(GetInterp(), "bell");
}

/*
 * Complements an individual contig.
 * Returns 0 for success
 *        -1 for failure
 */
int complement_contig(GapIO *io, tg_rec crec) {
    contig_t *c;
    bin_index_t *b;
    int len;
    int ustart, uend, delta;
    reg_complement rc;

    if (contig_lock_write(io, crec) == -1) {
	verror(ERR_WARN, "complement_contig", "Contig is busy");
	return -1;
    }

    if (!(c = (contig_t *)cache_search(io, GT_Contig, crec)))
	return -1;

    cache_incr(io, c);

    /*
     * Attempt to have symmetry in coordinates such that a contig
     * with used (unclipped) data from position A to B will be
     * complemented with the used portion of data still going from A to B.
     */
    consensus_valid_range(io, crec, &ustart, &uend);
    delta = (ustart - c->start) - (c->end - uend);

    /* Empty contig is special case */
    if (!contig_get_bin(&c)) {
	cache_decr(io, c);
	return 0;
    }

    if (!(b = (bin_index_t *)cache_search(io, GT_Bin, contig_get_bin(&c)))) {
	cache_decr(io, c);
	return -1;
    }

    if (!(b = cache_rw(io, b))) {
	cache_decr(io, c);
	return -1;
    }
    if (!(c = cache_rw(io, c))) {
	cache_decr(io, c);
	return -1;
    }

    b->flags ^= BIN_COMPLEMENTED;
    b->flags |= BIN_BIN_UPDATED;

    len = c->end - c->start + 1;
    b->pos -= b->size - ((c->start-b->pos)*2 + len);

    /* Shift by the difference between clipped vs unclipped coords */
    b->pos   += delta;
    c->start += delta;
    c->end   += delta;
    c->timestamp = io_timestamp_incr(io);

    cache_flush(io);
    
    rc.job = REG_COMPLEMENT;
    contig_notify(io, crec, (reg_data *)&rc);

    cache_decr(io, c);

    return 0;
}
