/*
 * File: list.c
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

#include <stdio.h>
#include <stdlib.h>
#include "list.h"
#include "xalloc.h"

/* static long item_id = 0; */

/* item 'id' allocation */
/* #define next_id() (++item_id) */

/*
 * Creates a new list
 *
 * Returns:
 *    A list pointer, or NULL for error.
 */
list_t *new_list() {
    list_t *lp;

    if (NULL == (lp = (list_t *)xcalloc(1, sizeof(list_t)))) {
	return NULL;
    }

    return lp;
}

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
void free_list(list_t *lp, void (*func)(void *data)) {
    item_t *i = head(lp), *in;

    if (func) {
	while(i) {
	    in = i->next;
	    func(i->data);
	    xfree(i);
	    i = in;
	}
    } else {
	while(i) {
	    in = i->next;
	    xfree(i);
	    i = in;
	}
    }

    xfree(lp);
}

/*
 * Adds an entry to a list
 *
 * Arguments:
 *     lp	- a list pointer
 *     dp	- the pointer to data to add
 *
 * Returns:
 *    >0 for success, value == item 'id'.
 *    -1 for failure (error reporting performed)
 */
int add_item(list_t *lp, void *dp) {
    item_t *i = new_item();

    if (!i) {
	return -1;
    }

    i->next = NULL;
    i->data = dp;
    /* i->id   = next_id(); */

    if (tail(lp)) {
	tail(lp)->next = i;
	tail(lp) = i;
    } else {
	head(lp) = tail(lp) = i;
    }

    return 0 /* i->id */ ;
}

/*
 * **INTERNAL**
 * Finds item previous to that containing 'data' within a list, starting at
 * item 'ip'.
 *
 * Arguments:
 *    ip	- item pointer
 *    dp	- data pointer
 *    comp	- comparitive function, as in strcmp().
 *
 * Returns:
 *    Item with 'next' pointing to the item containing dp.
 */
item_t *find_item(item_t *ip, void *dp, int (*comp)(void *d1, void *d2)) {
    while (ip->next) {
	if (comp(ip->next->data, dp) == 0)
	    return ip;
	
	ip = ip->next;
    }

    return NULL;
}

/*
 * Deletes data from a list.
 *
 * Arguments:
 *    lp	- list pointer
 *    dp	- data pointer
 *    del	- Optional function to call when removing the data. Set to NULL
 *		  if data shouldn't be freed.
 *    comp	- Comparitive function, eg strcmp().
 *    del_all	- 'set' if we wish to delete all occurances of dp in lp.
 *		  Otherwise only delete the first.
 *
 * Returns:
 *     0 for success; -1 for failure (dp not found in lp).
 */
int delete_item(list_t *lp, void *dp, void (*del)(void *data),
		int (*comp)(void *d1, void *d2), int del_all) {
    item_t *i = (item_t *)lp, *tmp;
    int found = -1;

    do {
	i = find_item(i, dp, comp);
	if (!i)
	    break;

	found = 0;

	if (i->next->data && del)
	    del(i->next->data);
    
	tmp = i->next;
	i->next = i->next->next;

	if (tail(lp) == tmp) {
	    tail(lp) = i->next ? i->next : head(lp);
	}

	xfree(tmp);
    } while (del_all);

    return found;
}

/*
 * Debugging - simply prints the list and pointers
 */
void dump_list(list_t *lp) {
    item_t *i = head(lp);

    printf("%p LIST %p %p\n", (void *)lp, (void *)head(lp), (void *)tail(lp));
    while(i) {
	printf("%p %p(%ld)\n", (void *)i, i->data, (long)(i->data));
	i = i->next;
    }
}

