/*
 * File: list.h
 * Version: 1.0
 *
 * Original Author: James Bonfield
 *
 * Description:
 *     Yet another generalised list package. (It's simpler to write one than
 *     hunt around for a free system.)
 *
 * Created: 02/08/94
 */

#ifndef _LIST_H
#define _LIST_H

#include <errno.h>

/*
 * ---------------------------------------------------------------------------
 * Data types
 * ---------------------------------------------------------------------------
 */

/*
 * A singly linked list of items
 */
typedef struct _item {
    struct _item	*next;
    void		*data;
    /* int		id; */
} item_t;

/*
 * The list header - also treated as an item under special circumstances!
 * Probably easiest to do away with last entirely and always add to start of
 * list. This also allows for removal of 'data' from this structure.
 */
typedef struct _list {
    item_t *first;
    void   *data;   /* always NULL - a dummy record to match struct _item */
    item_t *last;
} list_t;

/*
 * ---------------------------------------------------------------------------
 * Macros
 * ---------------------------------------------------------------------------
 */

/* list accessing */
#define head(list) ((list)->first)
#define tail(list) ((list)->last)

/* Create / delete items */
#define new_item() ((item_t *)xmalloc(sizeof(item_t)))
#define free_item(item) \
    do { \
	if (item) \
	    xfree(item); \
    } while (0)

/*
 * ---------------------------------------------------------------------------
 * Functions
 * ---------------------------------------------------------------------------
 */

/*
 * Creates a new list
 *
 * Returns:
 *    A list pointer, or NULL for error.
 */
extern list_t *new_list(void);

/*
 * Frees a list and all items within that list.
 * Can optionally also free all the data within that list
 *
 * Arguments:
 *    lp	- the list to free
 *    func	- a function to call if the data in the list is to be freed.
 *		  The function is called with the data as the argument.
 *		  Specify as NULL if data shouldn't be freed.
 */
extern void free_list(list_t *lp, void (*func)(void *data));

/*
 * Adds an entry to a list
 *
 * Arguments:
 *     lp	- a list pointer
 *     dp	- the pointer to data to add
 *
 * Returns:
 *     0 for success
 *    -1 for failure (error reporting performed)
 */
int add_item(list_t *lp, void *dp);

/*
 * Deletes data from a list.
 *
 * Arguments:
 *    lp	- list pointer
 *    dp	- data pointer
 *    func	- Optional function to call when removing the data. Set to NULL
 *		  if data shouldn't be freed.
 *    del_all	- 'set' if we wish to delete all occurances of dp in lp.
 *		  Otherwise only delete the first.
 *
 * Returns:
 *     0 for success; -1 for failure (dp not found in lp).
 */
int delete_item(list_t *lp, void *dp, void (*free)(void *data),
		int (*comp)(void *d1, void *d2), int del_all);

#endif /* _LIST_H */
