#include <stdio.h>
#include <stdlib.h>
#include "hash.h"
#include "xalloc.h"

#define TRUE 1
#define FALSE 0

void
InitialiseTable(TablePtr T[])
{
    int i;
        
    for (i = 0; i < HASHMODULUS; i++) {

	T[i] = NULL;

    } /* end for */

} /* end InitialiseTable */

int
Hash(int key) 
{
    return (key % HASHMODULUS);

} /* end Hash */

void
ChainInsert(TablePtr T[],
	    int newkey,
	    ItemType newinfo)
{

    int H_index;
    TablePtr new_node;
    
    H_index = Hash(newkey);
    new_node = (Table *)xmalloc(sizeof(Table));

    /* check able to allocate memory. If not, then exit */
    if (new_node == NULL){
	return;

    } /* end if */

    new_node->key = newkey;
    new_node->info = newinfo;
    new_node->next = T[H_index];
 
    T[H_index] = new_node;
    
} /* end ChainInsert */

void
ChainSearch(TablePtr T[],
	    int search_key,
	    int *found,
	    ItemType *search_info)
{
    int H_index;
    TablePtr node_ptr;

    *found = FALSE;
    H_index = Hash(search_key);

    node_ptr = T[H_index];

    while ((node_ptr != NULL) && (! *found)) {
	if (node_ptr->key == search_key) {
	    *found = TRUE;
	} else {
	    node_ptr = node_ptr->next;
	} 
    }
    if (*found) {
	*search_info = node_ptr->info;
    }

} /* end ChainSearch */

void 
ChainDelete(TablePtr T[], int search_key) {
    int H_index = Hash(search_key);
    TablePtr node_ptr = T[H_index], last_ptr = NULL;
    
    /* Find node */
    for (; node_ptr != NULL; last_ptr = node_ptr, node_ptr = node_ptr->next) {
	if (node_ptr->key == search_key)
	    break;
    }

    if (node_ptr) {
	/* Unlink */
	if (last_ptr) {
	    last_ptr->next = node_ptr->next;
	} else {
	    T[H_index] = node_ptr->next;
	}

	/* and free */
	xfree(node_ptr);
    }
} /* end ChainDelete */
