#ifndef _ACTIVE_TAGS_H_
#define _ACTIVE_TAGS_H_

#include "tagdb.h"

/*
 * Sets an active tag array from the 'list' string.
 * An empty 'list' sets the array to contain nothing, but a NULL list sets
 * the array to the default - all.
 */
int SetActiveTags2 (char *list, int *num, char ***types);
int SetActiveTags(char *list);

/*
 * Get the possible tag types from tagdb and set their booleans to 0.
 * Read the tag db into tag_db, which sets their number to tag_db_count
 */
void get_tag_types (void);

int idToIndex(char *id);
char *indexToId(int index);

extern char **active_tag_types;
extern int number_of_active_tags;

#endif /* _ACTIVE_TAGS_H_ */

