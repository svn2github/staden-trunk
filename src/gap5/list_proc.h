#ifndef _LIST_PROC_H_
#define _LIST_PROC_H_

#include "io_utils.h"

/*
 * Gets the active tcl list to 'list'.
 *
 * returns  0 for success
 * returns -1 for failure
 */
int set_active_list(char *list);

/*
 * Find the next string in the current active list.
 *
 * returns item for success
 * returns NULL for failure
 */
char *get_active_list_item(void);

/*
 * Reset the list index to the beginning again.
 */
int rewind_active_list(void);

int active_list_readings(GapIO *io, char *list, int *argc, tg_rec **argv);

int active_list_contigs(GapIO *io, char *list,
			int *argc, contig_list_t **argv);

int active_list_scaffold(GapIO *io, char *list,
			 int *argc, tg_rec **argv);

/* As above, but include cutoff data at ends of contig */
int active_list_contigs_extended(GapIO *io, char *list,
				 int *argc, contig_list_t **argv);

/*
 * Splits a list into argc & argv.
 * Returns 0 for success, -1 for failure.
 */
int SplitList(char *list, int *argc, char ***argv);


/*
 * Adds to an _existing_ list. ListCreate should be called prior to this.
 */
void add_to_list(char *name, char *item);

/*
 * Reads and returns a list contents.
 */
char *read_list(char *listname);

/*
 * Resets a Tcl NGList list
 */
void clear_list(char *listname);

/*
 * Removes duplicates from an NGList
 */
void list_remove_duplicates(char *listname);

/*
 * Access to the Tcl dynamic string mechanism to hold a list in C memory
 * rather than our "list" idea in Tcl. The routines here are to abstract away
 * the use of Tcl for storing this data. We have alloc, add and free functions.
 */
void *alloc_dlist(void);
char *add_to_dlist(void *dl, char *item);
void free_dlist(void *dl);
char *read_dlist(void *dl);

/*
 * Given a list of reading names or #nums, this returns a list of reading
 * numbers (given names are not unique) containing the original read and its
 * pair. Sets *nr to be the number of returned items.
 *
 * Returns malloced array of record numbers on success (caller frees)
 *         NULL on failure
 */
tg_rec *pair_readings(GapIO *io, char *inlist, int *nr);

#endif
