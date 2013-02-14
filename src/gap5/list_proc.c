#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include <io_lib/hash_table.h>
#include <ctype.h>

#include "io_utils.h"
#include "FtoC.h"
#include "newgap_cmds.h"
#include "misc.h"
#include "gap4_compat.h" /* read_name_to_number() */

static char **active_list_argv = NULL;
static int active_list_argc = 0;
static int active_list_index;

/*
 * Gets the active tcl list to 'list'.
 *
 * returns  0 for success
 * returns -1 for failure
 */
int set_active_list(char *list) {
    if (active_list_argv)
	Tcl_Free((char *)active_list_argv);

    if (Tcl_SplitList(GetInterp(), list,
		      &active_list_argc, &active_list_argv) != TCL_OK) {
	active_list_argv = NULL;
	active_list_argc = 0;
	return -1;
    }

    active_list_index = 0;

    return 0;
}

/*
 * return the size of the active list
 */
int active_list_size(void) {
    return active_list_argc;
}

/*
 * Find the next string in the current active list.
 *
 * returns item for success
 * returns NULL for failure
 */
char *get_active_list_item(void)
{
    if (active_list_index < active_list_argc)
	return active_list_argv[active_list_index++];
    else
	return NULL;
}


/*
 * Reset the list index to the beginning again.
 */
int rewind_active_list(void) {
    active_list_index = 0;

    return 0;
}

/* --------------------------------------------------------------------------
 * Lists menu interaction via the Tcl NGList array.
 */

/*
 * Adds to an list. Creates the list if it doesn't exist
 */
void add_to_list(char *name, char *item) {
    static char last_list[100];
    static int created = 0;

    if (!GetInterp())
	return;

    /*
     * Check and create the list, only check if we have really to.
     */
    if (strncmp(last_list, name, 100) != 0 || created == 0) {
	created = 1;
	strncpy(last_list, name, 100);

	if (Tcl_GetVar2(GetInterp(), "NGList", name, TCL_GLOBAL_ONLY) == NULL)
	    Tcl_VarEval(GetInterp(), "ListCreate2 ", name, " \"\" ",
			"SEQID", NULL);
    }

    Tcl_SetVar2(GetInterp(), "NGList", name, item,
		TCL_GLOBAL_ONLY | TCL_APPEND_VALUE | TCL_LIST_ELEMENT);
}

/*
 * Reads and returns a list contents.
 */
char *read_list(char *listname) {
    return Tcl_GetVar2(GetInterp(), "NGList", listname, TCL_GLOBAL_ONLY);
}

/*
 * Resets a Tcl NGList list
 */
void clear_list(char *listname, char *tag) {
    if (Tcl_GetVar2(GetInterp(), "NGList", listname, TCL_GLOBAL_ONLY))
	Tcl_SetVar2(GetInterp(), "NGList", listname, "", TCL_GLOBAL_ONLY);
}

/*
 * Removes duplicates from an NGList
 */
void list_remove_duplicates(char *listname) {
    Tcl_VarEval(GetInterp(), "ListRemoveDuplicates ", listname, NULL);
}

/* ------------------------------------------------------------------------ */

/*
 * Access to the Tcl dynamic string mechanism to hold a list in C memory
 * rather than our "list" idea in Tcl. The routines here are to abstract away
 * the use of Tcl for storing this data. We have alloc, add and free functions.
 */
void *alloc_dlist(void) {
    Tcl_DString *ds;

    ds = (Tcl_DString *)xmalloc(sizeof(*ds));
    if (ds == NULL)
	return NULL;

    Tcl_DStringInit(ds);

    return ds;
}

char *add_to_dlist(void *dl, char *item) {
    Tcl_DString *ds = (Tcl_DString *)dl;

    return Tcl_DStringAppendElement(ds, item);
}

char *read_dlist(void *dl) {
    Tcl_DString *ds = (Tcl_DString *)dl;

    return ds->string;
}

void free_dlist(void *dl) {
    Tcl_DString *ds = (Tcl_DString *)dl;

    Tcl_DStringFree(ds);
    xfree(ds);
}

/*
 *----------------------------------------------------------------------------
 * Simple utility functions to combine routines
 *----------------------------------------------------------------------------
 */
int active_list_readings(GapIO *io, char *list, int *argc, tg_rec **argv) {
    if (-1 == set_active_list(list))
	return -1;

    return lget_gel_num(io, active_list_argc, active_list_argv, argc, argv);
}

int active_list_contigs(GapIO *io, char *list,
			int *argc, contig_list_t **argv) {
    if (-1 == set_active_list(list))
	return -1;

    return lget_contig_num(io, active_list_argc, active_list_argv, argc, argv);
}

/* As above, but include cutoff data at ends of contig */
int active_list_contigs_extended(GapIO *io, char *list,
				 int *argc, contig_list_t **argv) {
    if (-1 == set_active_list(list))
	return -1;

    return lget_contig_num2(io, active_list_argc, active_list_argv,
			    argc, argv);
}

int active_list_scaffold(GapIO *io, char *list,
			 int *argc, tg_rec **argv) {
    if (-1 == set_active_list(list))
	return -1;

    return lget_scaffold_num(io, active_list_argc, active_list_argv,
			     argc, argv);
}

/*
 *----------------------------------------------------------------------------
 * Our interfaces into handy Tcl functions.
 * Used to keep direct Tcl usage in our C files to a minimum.
 *----------------------------------------------------------------------------
 */

/*
 * Splits a list into argc & argv.
 * Returns 0 for success, -1 for failure.
 */
int SplitList(char *list, int *argc, char ***argv) {
    return (Tcl_SplitList(GetInterp(), list, argc, argv) == TCL_OK) ? 0 : -1;
}


/* --------------------------------------------------------------------------
 * Read pairing.
 * Given a list of reading names or #nums, this returns a list of reading
 * numbers (given names are not unique) containing the original read and its
 * pair. Sets *nr to be the number of returned items.
 *
 * Returns malloced array of record numbers on success (caller frees)
 *         NULL on failure
 */
tg_rec *pair_readings(GapIO *io, char *inlist, int *nr) {
    int nread = 0, nalloc = 0;
    tg_rec *r = NULL, rec;
    char *cp;
    HashTable *h;
    HashItem *hi;
    HashIter *iter;
    
    h = HashTableCreate(1024, HASH_DYNAMIC_SIZE | HASH_POOL_ITEMS);

    while (*inlist && isspace(*inlist))
	inlist++;

    while (*inlist) {
	char tmp;

	for (cp = inlist; *cp && !isspace(*cp); cp++)
	    ;
	tmp = *cp;
	*cp = 0;

	rec = 0;
	if (*inlist == '#') {
	    char *endp;
	    rec = strtol(inlist+1, &endp, 10);
	    if (*endp != 0)
		rec = 0;
	}

	if (!rec)
	    rec = read_name_to_number(io, inlist);

	/* Insert rec and rec's pair */
	if (rec > 0) {
	    seq_t *s = cache_search(io, GT_Seq, rec);

	    if (s) {
		HashData hd;
		hd.i = 0;
		HashTableAdd(h, (char *)&rec, sizeof(rec), hd, NULL);

		rec = sequence_get_pair(io, s);
		if (rec > 0 && (s = cache_search(io, GT_Seq, rec))) {
		    HashTableAdd(h, (char *)&rec, sizeof(rec), hd, NULL);
		}
	    }
	}

	*cp = tmp;
	while (*cp && isspace(*cp))
	    cp++;

	inlist = cp;
    }

    /* We used a hash as it automatically drops duplicates.
     * Now iterate over the hair to produce a linear list for returning.
     */
    iter= HashTableIterCreate();
    while ((hi = HashTableIterNext(h, iter))) {
	if (nread >= nalloc) {
	    nalloc = nalloc ? nalloc * 2 : 256;
	    r = realloc(r, nalloc * sizeof(*r));
	    if (!r)
		return NULL;
	}
	r[nread++] = *(tg_rec *)hi->key;
    }

    HashTableIterDestroy(iter);
    HashTableDestroy(h, 0);

    *nr = nread;

    return r;
}
