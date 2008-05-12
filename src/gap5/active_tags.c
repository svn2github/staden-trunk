#include <stdio.h>

#include "tagdb.h"
#include "list_proc.h"
#include "misc.h"

/*
 *-----------------------------------------------------------------------------
 * Handle active tag list. We need to initialise it and then allow the user
 * to set his own selections. 
 *-----------------------------------------------------------------------------
 */
char **active_tag_types = NULL;
int number_of_active_tags = 0;

/*
 * Get the possible tag types from tagdb and set their booleans to 0.
 * Read the tag db into tag_db, which sets their number to tag_db_count
 */
void get_tag_types (void) {
    static int done = 0;

    if (!done) {
	readInTagDB();
	done = 1;
    }
}

/*
 * Sets an active tag array from the 'list' string.
 * An empty 'list' sets the array to contain nothing, but a NULL list sets
 * the array to the default - all.
 */
int SetActiveTags2 (char *list, int *num, char ***types) {
    if (*types)
	Tcl_Free((char *)*types);
 
    if (list) {
	if (SplitList(list, num, types) == -1) {
	    *types = NULL;
	    *num = 0;
	    return -1;
	}
    } else {
	int i;

	if (NULL == (*types = (char **)Tcl_Alloc(tag_db_count * sizeof(char *)))){
	    *num = 0;
	    return -1;
	}

	for (i = 0; i < tag_db_count; i++) {
	    (*types)[i] = tag_db[i].id;
	}
	*num = tag_db_count;
    }

    return 0;
}

int SetActiveTags (char *list) {
    return SetActiveTags2(list, &number_of_active_tags, &active_tag_types);
}


/*
 * Should use a hash table when we've got _large_ numbers of tags (which we
 * currently don't.
 */
int idToIndex(char *id)
{
    int i;
    if (id==NULL) return 0;
    for (i=0; i<tag_db_count; i++) {
        if (strncmp(id,tag_db[i].id,4)==0)
            return i;
    }
    return 0;
}

/*
 * Convert a tag id to a tag type string
 */
char *indexToId(int index) {
    return tag_db[index].search_id;
}
