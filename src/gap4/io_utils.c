#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "io_utils.h"
#include "misc.h"
#include "FtoC.h"


/*
 * ---------------------------------------------------------------------------
 * Name to number type conversions.
 */

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
#if GAP_CACHE==2
    /* Using tcl hash tables */
    char cname[128];
    Tcl_HashEntry *hash;

    /* Shouldn't be needed now, but we may have old databases */
    Fstr2Cstr(tname, DB_NAMELEN, cname, 128);
    
    if (!(hash = Tcl_FindHashEntry(&io->tname_hash, tname))) {
	return 0;
    }
    return (int)Tcl_GetHashValue(hash);

#else
    /* Via a linear search through all names */
    char cname[128];
    int i;
    GTemplates t;
    char buf[128];
	
    /* Shouldn't be needed now, but we may have old databases */
    Fstr2Cstr(tname, 128, cname, 128);
    
    for (i = 1; i <= Ntemplates(io); i++) {
	/* read template record */
	template_read(io, i, t);
	/* read template name */
	TextRead(io, t.name, buf, sizeof(buf));
	if (strcmp(buf, tname)==0)
	    return i;
    }
    return 0;
#endif
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
#if GAP_CACHE==2
    /* Using tcl hash tables */
    char cname[DB_NAMELEN+1];
    Tcl_HashEntry *hash;

    /* Shouldn't be needed now, but we may have old databases */
    Fstr2Cstr(gel_name, DB_NAMELEN, cname, DB_NAMELEN+1);
    
    if (!(hash = Tcl_FindHashEntry(&io->rname_hash, cname))) {
	return 0;
    }
    return (int)Tcl_GetHashValue(hash);

#else
    /* Via a search through all names (which may be cached in memory) */
    char cname[DB_NAMELEN+1];
    int z;
    int from, to;
    static int last = 0;
    int i;

    /* Shouldn't be needed now, but we may have old databases */
    Fstr2Cstr(gel_name, DB_NAMELEN, cname, DB_NAMELEN+1);
    
    /* Reset last when we've reopened a database */
    if (last > NumReadings(io))
	last = NumReadings(io);

    for (z = 0; z < 2; z++) {
	from = z ? 0    : last;
	to   = z ? last : NumReadings(io);

	for (i = from; i < to; i++) {
	    char *name = get_read_name(io, i+1);
	    if (*name) {
		if (strncmp(name, cname, DB_NAMELEN) == 0) {
		    return (last = i)+1;
		    break;
		}
	    }
	}
    }
    return 0;
#endif
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
    int gel_num = -1;
     
    if (is_name == GGN_ID) {
	/* Maybe specifying by contig id? */
	if (*gel_name == '=') {
	    int num = atoi(gel_name+1);
		
	    if (!num)
		return -1;

	    return io_clnbr(io, num);
	}
    }

    /* simplest form - a number */
    if (is_name == GGN_ID && *gel_name == '#') {
	gel_num = atoi(gel_name+1);
    } else {
	gel_num = read_name_to_number(io, gel_name);
    }

    /*
     * Check gel_num is valid - we may have typed "999999"
     */
    return (gel_num < 1 || gel_num > NumReadings(io)) ? -1 : gel_num;
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
    char *visited;
    int l;

    /*
     * Optimise: if we have the contig number cached in rnum2cnum then
     * we can shortcut this by jumping straight to the first reading number.
     */
    if (io->cached_rnum2cnum && (l = arr(int, io->rnum2cnum, gel-1))) {
	return io_clnbr(io, l);
    }

    if (NULL == (visited = (char *)xcalloc(NumReadings(io)+1, 1)))
	return -1;

    while (l = io_lnbr(io,gel)) {
	if (visited[l]) {
	    verror(ERR_FATAL, "chain_left",
		   "Loop detected: %d found previously\n", l);
	    xfree(visited);
	    return -1;
	}
	visited[gel = l] = 1;
    }

    xfree(visited);
    return gel;
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
    int i, lgel, cnum;

    if (io->cached_rnum2cnum) {
	cnum = arr(int, io->rnum2cnum, gel-1);
	if (cnum) {
	    return cnum;
	}
    }

    /* First chain left */
    if (-1 == (lgel = chain_left(io, gel)))
	return -1;

    /* Then search for the contig with this left neighbour */
    for (cnum = 1; cnum <= NumContigs(io); cnum++) {
	if (io_clnbr(io, cnum) == lgel) {
	    if (io->cached_rnum2cnum) {
		/* Now cache details for this entire contig */
		for (i = lgel; i; i = io_rnbr(io, i))
		    arr(int, io->rnum2cnum, i-1) = cnum;
	    }
	    return cnum;
	}
    }

    return -1;
}

/*
 * Updates the reading number to contig number cache.
 * Specify contig as zero to remove this reading from the cache.
 */
void update_rnumtocnum(GapIO *io, int gel, int contig) {
    ArrayRef(io->rnum2cnum, gel-1);
    arr(int, io->rnum2cnum, gel-1) = contig;
}

/*
 * Invalidates the entire reading number to contig number cache.
 */
void invalidate_rnumtocnum(GapIO *io, int disable) {
    int i, nr = NumReadings(io);

    ArrayRef(io->rnum2cnum, nr-1);
    for (i = 1; i <= nr; i++) {
	arr(int, io->rnum2cnum, i-1) = 0;
    }
    io->cached_rnum2cnum = !disable;
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

    /* Get the gel number */
    gel = get_gel_num(io, gel_name, is_name);
    if (gel == -1)
	return -1;
    
    /* And return it's contig number */
    return rnumtocnum(io, gel);
}


/*
 * Fortran interface (INTEGER NAMENO) to get_gel_num()
 */
f_int nameno_(char *name, f_int *ngels, f_int *handle, char name_l) {
    GapIO *io;
    int ret;

    if ((io = io_handle(handle)) == NULL)
	return 0;

    ret = get_gel_num(io, name, GGN_NAME);
    return ret > 0 ? ret : 0;
}

/*
 * ---------------------------------------------------------------------------
 * List versions of get_gel_num and get_contig_num that are more efficient
 * when dealing with several strings at a time.
 */

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

#if GAP_CACHE==2
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

#else

    /* Translate reading names to reading numbers */
    for (n=NumReadings(io), i=1; count < listArgc && i<=n; i++) {
	char *name;

	name = get_read_name(io, i);
	for (j=0; j<listArgc; j++) {
	    if (strncmp(listArgv[j], name, DB_NAMELEN) == 0) {
		(*rargv)[j] = i;
		count++;
		/* break; */
	    }
	}
    }
#endif

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
    int i, j, n, count=0;
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
	    if (num > 0 && num <= NumReadings(io)) {
		(*rargv)[j].contig = num;
		count++;
	    } else {
		(*rargv)[j].contig = 0;
	    }
	} else if (listArgv[j][0] == '=') {
	    int num = atoi(&listArgv[j][1]);
	    if (num > 0 && num <= NumContigs(io)) {
		(*rargv)[j].contig = io_clnbr(io, num);
		count++;
	    } else {
		(*rargv)[j].contig = 0;
	    }
	} else {
	    (*rargv)[j].contig = 0;
	}
    }

    /* Translate reading names to reading numbers */
    for (j = 0; j < listArgc; j++) {
	if (0 == (*rargv)[j].contig) {
	    int rnum = read_name_to_number(io, listArgv[j]);
	    if (rnum) {
		(*rargv)[j].contig = rnum;
		count++;
	    }
	}
    }

    /*
     * Handle case when we've failed to find some; we just shuffle down to
     * fill any holes in the rargv structure.
     */
    if (count != listArgc) {
	for (i=j=0; j<listArgc; j++) {
	    if ((*rargv)[j].contig != 0) {
		(*rargv)[i++] = (*rargv)[j];
	    }
	}
    }


    /*
     * Chain left and negate to distinguish those we've done and those that we
     * haven't. Also set end value correctly if it's now 0.
     */
    for (j = 0; j < count; j++) {
	(*rargv)[j].contig = -chain_left(io, (*rargv)[j].contig);
    }

    /* Translate reading numbers to contig numbers, removing duplicates */
    count = 0;
    for (n=NumContigs(io), i=1; i<=n; i++) {
	int to_find = 1;

	for (j=0; j<listArgc; j++) {
	    if (io_clnbr(io, i) == -(*rargv)[j].contig && to_find) {
		(*rargv)[j].contig = i;
		if ((*rargv)[j].end == 0)
		    (*rargv)[j].end = ABS(io_clength(io, i));
		count++;
		to_find = 0;
	    }
	}

	if (count == listArgc)
	    break;
    }

    /* Check for failures */
    if (count != listArgc) {
	for (i=j=0; j<listArgc; j++) {
	    if ((*rargv)[j].contig > 0) {
		(*rargv)[i++] = (*rargv)[j];
	    }
	}
    }

    *rargc = count;
    return 0;
}

/*
 * Convert a "contig_list_t *" array into an "int *" array containing only
 * the contig elements.
 * This is useful for routines that only work on whole contigs.
 */
int *to_contigs_only(int num_contigs, contig_list_t *cl) {
    int *contigs;
    int i;

    /* +1 to always ensure we never malloc 0 bytes */
    if (NULL == (contigs = (int *)xmalloc(num_contigs * sizeof(int) + 1)))
	return NULL;

    for (i = 0; i < num_contigs; i++)
	contigs[i] = cl[i].contig;

    return contigs;
}

/*
 * ---------------------------------------------------------------------------
 * Number to name type conversions.
 */

void cache_template_name(GapIO *io, int number, char *name) {
#if GAP_CACHE==2
    if (*name) {
	Tcl_HashEntry *hash;
	int new;
	hash = Tcl_CreateHashEntry(&io->tname_hash, name, &new);
	Tcl_SetHashValue(hash, (ClientData)number);
    }
#endif
}


void cache_read_name(GapIO *io, int number, char *name) {
#if GAP_CACHE==1
    /* Update cached string */
    strncpy(arr(name_t, io->read_names, number-1).name, name, DB_NAMELEN+1);
#elif GAP_CACHE==2
    /* Update hash tables */
    if (*name) {
	Tcl_HashEntry *hash;
	int new;
	hash = Tcl_CreateHashEntry(&io->rname_hash, name, &new);
	arr(Tcl_HashEntry *, io->read_names, number-1) = hash;
	Tcl_SetHashValue(hash, (ClientData)number);
    } else {
	arr(Tcl_HashEntry *, io->read_names, number-1) = NULL;
    }
#else
    /* Do nothing */
#endif
}

void cache_delete_read_name(GapIO *io, int number) {
    char *name;

#if GAP_CACHE==1
    /* Update cached string */
    *(arr(name_t, io->read_names, number-1).name) = 0;
#elif GAP_CACHE==2
    /* Update hash tables */
    name = get_read_name(io, number);
    if (name && *name) {
	Tcl_HashEntry *hash;

	hash = Tcl_FindHashEntry(&io->rname_hash, name);
	if (hash)
	    Tcl_DeleteHashEntry(hash);

	arr(Tcl_HashEntry *, io->read_names, number-1) = NULL;
    }
#else
    /* Do nothing */
#endif
}

/*
 * Obtains a reading name from a reading number.
 * The returned value is only valid until the next call of get_read_name().
 */
char *get_read_name(GapIO *io, int number) {
#if GAP_CACHE==1
    /* Read from io->read_names cache */
    if (!arr(name_t, io->read_names, number-1).name[0]) {
	GReadings r;
	char name[DB_NAMELEN+1];

	gel_read(io, number, r);
	TextRead(io, r.name, name, DB_NAMELEN+1);
	Fstr2Cstr(name, DB_NAMELEN, name, DB_NAMELEN+1);
	strcpy(arr(name_t, io->read_names, number-1).name, name);
    }
    return arr(name_t, io->read_names, number-1).name;

#elif GAP_CACHE==2
    Tcl_HashEntry *h;
    /* Read from io->read_names cache */
    if ((h = arr(Tcl_HashEntry *, io->read_names, number-1))) {
	return Tcl_GetHashKey(&io->rname_hash, h);
    } else {
	/* Read from disk */
	static char name[DB_NAMELEN+1];
	GReadings r;

	gel_read(io, number, r);
	TextRead(io, r.name, name, DB_NAMELEN);
	Fstr2Cstr(name, DB_NAMELEN, name, DB_NAMELEN+1);

	cache_read_name(io, number, name);

	return name;
    }
#else
    /* Read from disk */
    static char name[DB_NAMELEN+1];
    GReadings r;

    gel_read(io, number, r);
    TextRead(io, r.name, name, DB_NAMELEN);
    Fstr2Cstr(name, DB_NAMELEN, name, DB_NAMELEN+1);

    return name;
#endif
}

/*
 * Obtains a contig name (leftmost reading name) from a contig number.
 * The returned value is only valid until the next call of get_contig_name().
 */
char *get_contig_name(GapIO *io, int number) {
    static char name[DB_NAMELEN+1];

    strcpy(name, get_read_name(io, io_clnbr(io, number)));

    return name;
}


char *get_vector_name(GapIO *io, int vector) {
    static char name[1025];
    GVectors v;

    if (vector > Nvectors(io))
	return "???";

    vector_read(io, vector, v);
    if (TextRead(io, v.name, name, sizeof(name)-1))
	return "???";
    name[sizeof(name)-1] = 0;

    return name;
}

char *get_template_name(GapIO *io, int template) {
    static char name[1025];
    GTemplates t;

    if (template > Ntemplates(io))
	return "???";

    template_read(io, template, t);
    if (TextRead(io, t.name, name, sizeof(name)-1))
	return "???";
    name[sizeof(name)-1] = 0;

    return name;
}

char *get_clone_name(GapIO *io, int clone) {
    static char name[1025];
    GClones c;

    if (clone > Nclones(io))
	return "???";

    clone_read(io, clone, c);
    if (TextRead(io, c.name, name, sizeof(name)-1))
	return "???";
    name[sizeof(name)-1] = 0;

    return name;
}
