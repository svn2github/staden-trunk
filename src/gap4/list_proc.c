#include <tcl.h>
#include <stdlib.h>
#include <string.h>

#include "io_utils.h"
#include "FtoC.h"
#include "newgap_cmds.h"
#include "misc.h"

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


/*
 * Adds to an list. Creates the list if it doesn't exist
 */
void add_to_list(char *name, char *item) {
    static char last_list[100];
    static int created = 0;

    /*
     * Check and create the list, only check if we have really to.
     */
    if (strncmp(last_list, name, 100) != 0 || created == 0) {
	created = 1;
	strncpy(last_list, name, 100);

	if (Tcl_GetVar2(GetInterp(), "NGList", name, TCL_GLOBAL_ONLY) == NULL)
	    Tcl_VarEval(GetInterp(), "ListCreate2 ", name, " \"\"", NULL);
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
 * Fortran interfaces
 *----------------------------------------------------------------------------
 */

/*
 * Get next reading (an interface to get_alist_item).
 */
int gnread_(char *ITEM, f_implicit ITEM_l)
{
    char *citem;

    if (NULL == (citem = get_active_list_item()))
	return 1;

    Cstr2Fstr(citem, ITEM, ITEM_l);

    return 0;
}

/* fortran interface to add_to_list() */
void *tolist_(char *name_p, char *info_p,
	      f_implicit name_l, f_implicit info_l) {
    char name[256], info[256];
    static char *last = NULL;
    static void *dl = NULL;

    /* Init */
    if (name_p == NULL && name_l != 0) {
	dl = NULL;
	last = NULL;
	return NULL;
    }

    /* Return list */
    if (name_p == NULL && name_l == 0) {
	return dl;
    }

    if (last != name_p) {
	last = name_p;

	dl = alloc_dlist();
    }

    Fstr2Cstr(name_p, name_l, name, 255);
    Fstr2Cstr(info_p, info_l, info, 255);

    add_to_dlist(dl, info);

    return NULL;
}

/*
 *----------------------------------------------------------------------------
 * Simple utility functions to combine routines
 *----------------------------------------------------------------------------
 */
int active_list_readings(GapIO *io, char *list, int *argc, int **argv) {
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
