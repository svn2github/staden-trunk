#ifndef _FREE_TREE_H_
#define _FREE_TREE_H_

#include "g-error.h"
#include "g-os.h"

typedef struct free_tree_n_ {
    /* Tree node pointers */
    struct free_tree_n_ *left; 
    struct free_tree_n_ *right;
    struct free_tree_n_ *parent;

    /* Linked lists threaded through the tree; groups nodes of similar size */
    struct free_tree_n_ *next;
    struct free_tree_n_ *prev;

    /* AVL balance: -1, 0, +1 for left heavy, balanced, right heavy */
    int balance;

    /* The node contents themselves */
    GImage pos;
    GImage len;
} free_tree_n;

typedef struct free_tree_ {
    free_tree_n *root;
    free_tree_n *wilderness;

    free_tree_n **node_blocks;  /* node block allocation system */
    int nblocks;
    free_tree_n *free_nodes;    /* follow left pointer for linked list */

    /*
     * Maintain a series of linked lists containing blocks of specific ranges.
     * In malloc library terms this is sometimes referred to as a
     * "segregated fit".
     * The first 57 lists are from 0 to 4080 inclusive in steps of
     * 16 (initially), 32 and 64 at the end. 
     * The remainder are 4096, 8192, 16384, ... in powers of 2.
     */
#define FBLOCKM 4095
#define FBSZ 16
#define NFBLOCK 57
#define NFBLOCK2 (NFBLOCK+64)
    int listsz[FBLOCKM/FBSZ+2];
    free_tree_n *lists[NFBLOCK2];
} free_tree;


/*
 * Allocates a free tree, initialising it to a single free node starting
 * at pos, extending to pos+len.
 *
 * Returns a free_tree pointer for success
 *         NULL for failure
 */
free_tree *freetree_create(GImage pos, GImage len);

/*
 * Deallocates memory used by a free_tree
 */
void freetree_destroy(free_tree *t);

/*
 * Registers a (pos+len) pair with the free tree, effectively excluding that
 * space from the tree.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int freetree_register(free_tree *t, GImage pos, GImage len);

/*
 * Allocates a space from the free tree. If this doesn't exist, then space
 * from the first node (aka end of the disk file) is used.
 * 
 * Returns offset for success
 *        -1 for failure
 */
GImage freetree_allocate(free_tree *t, GImage len);

/*
 * Adds some space (pos+len) back to the free_tree.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int freetree_unregister(free_tree *t, GImage pos, GImage len);

#define FREETREE_BLOCK(len, block_size) (\
    ((len) % (block_size)) \
    ? ((len) - ((len) % (block_size)) + (block_size)) \
    : (len))

#ifdef CACHE_FREETREE
/*
 * Saves a freetree to disk at the end of the file pointed to by 'fd'.
 * To verify validity for subsequent reads, we also store the current global
 * time stamp along with a 32-bit CRC of the entire tree.
 * (If the tree becomes corrupted then we would inevitably corrupt the entire
 * database.)
 *
 * Returns:
 *	0 for success
 *	-1 for failure
 */
int freetree_save_int4(free_tree *t, int fd, int last_time);
int freetree_save_int8(free_tree *t, int fd, int last_time);

/*
 * Loads a freetree from the current location in the file pointed to by fd.
 * We check that the last_time matches the one passed in here and that the
 * 32-bit CRC is correct.
 * If any of these operations fail we return NULL and the tree will be
 * recomputed from scratch.
 *
 * Returns:
 *	A free_tree pointer for success.
 *	NULL for failure.
 */
free_tree *freetree_load_int4(int fd, int last_time);
free_tree *freetree_load_int8(int fd, int last_time);
#endif

#endif /* _FREE_TREE_H_ */
