#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#define NDEBUG
#include <assert.h>
#include "freetree.h"
#include "xalloc.h"
#include "g-os.h"

/*
 * TODO:
 * - Use a block allocating system. According to dmalloc the memory usage
 *   of this module is 50% administration (aka malloc headers) so we could
 *   feasibly save nearly half of our memory.
 *
 * - Use an AVL tree height balancing system to reduce the time taken by
 *   the tree_find_pos algorithm.
 *
 * - Add check to tree_find_pos to error on cases where the block isn't there.
 *   This will spot the 'overlapping image' cases, I hope.
 *
 * - Remove requirement for parent? Or maybe use a threaded tree. To remove
 *   the parent we need to change some algorithms to be recursive. These
 *   appear to be not often ran algorithms except for the length search. This
 *   may be better done by searching the allocated blocks space as mentioned
 *   in the first optimisation above. The tree_delete_node() routine has to
 *   modify the left or right pointers of it's parent node, but this node is
 *   always known before calling tree_delete_node() so it isn't a problem.
 *   (We'd have to change the search functions to return the parent as well
 *   as the current node.)
 *
 * - Remove hack of having the root node always containing upto INT_MAX in
 *   space. We should simply recode to keep track of the last value in the
 *   tree somewhere else. This would simplify things in the long run as
 *   this node at present would cause problems with an AVL balancing system.
 */

/*
 * Allocates a new free_tree_n node.
 * Because we may allocate many nodes we have our own allocator which
 * allocates blocks of 100 nodes at a time. This is primarily to reduce the
 * malloc overhead.
 *
 * Returns free_tree pointer for success
 *         NULL for failure
 */
#define BLOCK_FACTOR 100
static free_tree_n *new_node(free_tree *t) {
    free_tree_n *n;

    /* Is there a free node? If not add a new block */
    if (!t->free_nodes) {
	int i;
	free_tree_n *prev;

	t->nblocks++;
	t->node_blocks = (free_tree_n **)xrealloc(t->node_blocks, t->nblocks *
						  sizeof(free_tree_n *));
	n = t->node_blocks[t->nblocks-1] =
	    (free_tree_n *)xmalloc(BLOCK_FACTOR * sizeof(free_tree_n));

	prev = NULL;
	for (i = 0; i < BLOCK_FACTOR; i++) {
	    n[i].left = prev;
	    n[i].right = NULL;
	    n[i].parent = NULL;
	    n[i].pos = 0;
	    n[i].len = 0;
	    prev = &n[i];
	}
	t->free_nodes = &n[BLOCK_FACTOR-1];
    }
	
    /* pick one off the free list */
    n = t->free_nodes;
    t->free_nodes = n->left;

    n->left = NULL;
    n->right = NULL;
    n->parent = NULL;
    n->pos = 0;
    n->len = 0;

    return n;
}

static void del_node(free_tree *t, free_tree_n *n) {
    n->left = t->free_nodes;
    n->right = NULL;
    n->parent = NULL;
    t->free_nodes = n;
}

/*
 * Search for a tree element in which pos is a member of. We scan the tree
 * in an inorder fasion.
 *
 * Complexity: O(log(N)) for a balanced tree (which these aren't).
 *
 * Always returns a free_tree_n pointer.
 */
static free_tree_n *tree_find_pos(free_tree_n *t, int pos) {
    for (;;) {
	if (pos < t->pos) {
	    /* left child */
	    if (t->left) {
		t = t->left;
	    } else {
		return t;
	    }

	} else if (t->pos + t->len > pos) {
	    /* This segment */
	    return t;

	} else {
	    /* right child */
	    if (t->right) {
		t = t->right;
	    } else {
		return t;
	    }
	}
    }
}


/*
 * Search for a tree element in which the (pos,len) pair butt up against, or
 * NULL if there are none.
 *
 * Complexity: O(log(N)) for a balanced tree (which these aren't).
 *
 * Returns the node next to (pos,len) for success
 *         NULL for not found.
 *
 * status is filled with -1 for error (overlap found)
 *                        0 for no butt (inbetween free blocks)
 *                        1 for butt with left end
 *                        2 for butt with right end
 */
static free_tree_n *tree_find_pos_len(free_tree_n *t, int pos, int len,
				      int *status) {
    int end;

    *status = 0;
    for (;;) {
	if (pos < t->pos) {
	    /* NB: should check and error when pos + len > t->pos */
	    assert(pos + len <= t->pos);
	    if (pos + len == t->pos) {
		*status = 1;
		return t;
	    }

	    /* left child */
	    if (t->left) {
		t = t->left;
	    } else {
		/* Not found */
		return t;
	    }

	} else {
	    end = t->pos + t->len;
	    if (pos > end) {
		/* right child */
		if (t->right) {
		    t = t->right;
		} else {
		    /* Not found */
		    return t;
		}
	    } else {
		if (pos == end) {
		    *status = 2;
		    return t;
		}
		

		/* Overlap */
		*status = -1;
		return NULL;
	    }
	}
    }
}


/*
 * Search the tree for the first block greater than a particular length.
 * We scan the tree in an preorder fashion as it will on average visit less
 * nodes before finding a gap and the position ordering is not relevant.
 * We use a 'roving pointer' to remember the last place we searched in order
 * to try and reduce small blocks accumulating at the start.
 *
 * Complexity: O(N/f)? (f is some fragmentation value)
 *
 * Returns a free_tree_n pointer when one is found.
 *           NULL for failure
 */
static free_tree_n *tree_find_len(free_tree *tr, int len) {
    free_tree_n *root = tr->tree->left;
    free_tree_n *rover = tr->rover;
    free_tree_n *t;
    int i;

    /* Do we know the largest block? */
    if (tr->largest) {
	if (len <= tr->largest->len) {
	    free_tree_n *ret = tr->largest;
	    tr->largest = NULL;
	    return tr->rover = ret;
	} else {
	    goto quick;
	}
    }

    /* puts("============== TREE_FIND_LEN ================"); */

    if (NULL == rover)
	rover = tr->rover = root;

    t = rover;

    for (i=0; i < 2; i++) {
	for (; t->parent;) {
	    if (t == rover && i == 1)
		break;

	    /* printf("pos=%d,len=%d\n", t->pos, t->len); */
	    /* Check node */
	    if (t->len >= len)
		return tr->rover=t;

	    /* Traverse left */
	    if (t->left) {
		t = t->left;
		continue;
	    }
	    
	    /* Traverse right */
	    if (t->right) {
		t = t->right;
		continue;
	    }
	    
	    /* Backtrack */
	    for (;t->parent;) {
		/* from left, go to right */
		if (t->parent->left == t && t->parent->right) {
		    t = t->parent->right;
		    break;
		}
		
		/* from left but no right, or from right */
		t = t->parent;
	    }
	}
	
	/*
	 * If we've fallen out of the above loop then it means we've hit the
	 * parent (NOTE: there is no right neighbour to the parent) or
	 * the rover (second loop around).
	 *
	 * In the first case we scan the bit that starting from the rover
	 * missed by now starting from the root node and rechecking everything
	 * again (this time stopping at the rover).
	 */
	if (i == 0) {
	    if (t->left == NULL)
		break;
	    t = t->left;
	}

	if (rover == root)
	    break;
    }

 quick:

    fflush(stdout);
    tr->rover = root;
    return root->len >= len ? root : NULL;
}


/*
 * Walks 'left' along a tree one node in an inorder order.
 *
 * Returns the predecessor node when successful
 *             NULL for failure
 */
static free_tree_n *tree_walk_left(free_tree_n *n) {
    if (n->left) {
	/* go left and chain right */
	for (n = n->left; n->right; n = n->right)
	    ;
	return n;
    }

    /*
     * otherwise we backtrack (possibly zero steps) until we're the
     * right child of our parent.
     */
    for (; n->parent && n->parent->left == n; n = n->parent)
	;

    return n->parent;
}


/*
 * Walks 'right' along a tree one node in an inorder order.
 *
 * Returns the successor node when successful
 *             NULL for failure
 */
static free_tree_n *tree_walk_right(free_tree_n *n) {
    if (n->right) {
	/* go right and and chain left */
	for (n = n->right; n->left; n = n->left)
	    ;
	return n;
    }

    /*
     * otherwise we backtrack (possibly zero steps) until we're the
     * left child of our parent.
     */
    for (; n->parent && n->parent->right == n; n = n->parent)
	;

    return n->parent;
}


/*
 * Deletes a node - we DON'T check if the node to delete is the root, but
 * this would cause an invalid tree anyway.
 */
static void tree_delete_node(free_tree *t, free_tree_n *node) {
    free_tree_n **rover = &t->rover;

    if (node == *rover)
	*rover = NULL;

    if (node == t->largest)
	t->largest = NULL;

    if (node->left == NULL && node->right == NULL) {
	/* Leaf node */
	if (node->parent->left == node)
	    node->parent->left = NULL;
	else
	    node->parent->right = NULL;

    } else if (node->left == NULL) {
	/* Only right child - replace by our right */
	if (node->parent->left == node)
	    node->parent->left = node->right;
	else
	    node->parent->right = node->right;
	node->right->parent = node->parent;

    } else if (node->right == NULL) {
	/* Only left child - replace by our left */
	if (node->parent->left == node)
	    node->parent->left = node->left;
	else
	    node->parent->right = node->left;
	node->left->parent = node->parent;

    } else {
	free_tree_n *succ;
	/*
	 * This is an interior node - replace with next ordered item in
	 * the tree. We can do this by taking the nodes right child and
	 * following left.
	 */
	for (succ = node->right; succ->left; succ=succ->left)
	    ;

	/* Make sure we did actually manage to go left... */
	if (succ->parent != node) {
	    succ->parent->left = succ->right;
	    if (succ->right)
		succ->right->parent = succ->parent;
	}

	/* Unlink node */
	if (node->parent->left == node) {
	    node->parent->left = succ;
	} else {
	    node->parent->right = succ;
	}
	succ->left = node->left;
	succ->left->parent = succ;
	if (node->right != succ) {
	    succ->right = node->right;
	    succ->right->parent = succ;
	}
	succ->parent = node->parent;
    }

    del_node(t, node);
}

static int last_end;
static int tree_sum;

/*
 * Prints a tree to stdout in the "inorder" order.
 */
static void tree_print_r(free_tree_n *t, int indent) {
    if (t->left) {
	/*putchar('l');*/
	assert(t->left->parent == t);
	tree_print_r(t->left, indent+1);
    }
    printf("pos=%d, len=%d\n", t->pos, t->len);
    tree_sum += t->pos;
    assert(t->pos > last_end);
    last_end = t->pos + t->len;
    if (t->right) {
	/*putchar('r');*/
	assert(t->right->parent == t);
	tree_print_r(t->right, indent-1);
    }
    /*putchar('b');*/
}

void tree_print(free_tree *t) {
    puts("============== TREE ============");
    last_end = -1;
    tree_sum = 0;
    tree_print_r(t->tree, 0);
    printf("Tree sum=%d\n", tree_sum);
}


/*
 * Generates postscript to print the tree. It's pretty rough and ready, with
 * constant scale factors (hand edit the postscript to fit!), but serves
 * well for quick debugging and analysis of how the tree is balanced.
 */
void tree_print_ps(free_tree_n *t) {
    float step=65536;
    float yinc = 10000;
    float degrade = .9;

    puts("%!PS");
    puts("500 380 translate 90 rotate newpath 0 0 moveto .003 .002 scale");

    for (; t->parent;) {
	if (!t->left) {
	    if (!t->right) {
		for (;t->parent; step*=2, yinc/=degrade) {
		    if (t->parent->left == t && t->parent->right) {
			t = t->parent->right;
			printf("%f %f rmoveto\n", step*2, -yinc/degrade);
			printf("%f %f rlineto\n", step*2, yinc/degrade);
			break;
		    }
		    if (t->parent->left == t)
			printf("%f %f rmoveto\n", step*2, -yinc/degrade);
		    else
			printf("%f %f rmoveto\n", -step*2, -yinc/degrade);
		    t=t->parent;
		}
	    } else {
		printf("%f %f rlineto\n", step, yinc);
		step/=2;
		yinc*=degrade;
		t = t->right;
	    }
	} else {
	    printf("%f %f rlineto\n", -step, yinc);
	    step/=2;
	    yinc*=degrade;
	    t = t->left;
	}
    }

    puts("stroke showpage");
}


int tree_total_depth(free_tree_n *t) {
    int depth = 0, total_depth = 0;
    int nodes = 0;

    for (; t->parent;) {
	/* Traverse left */
	if (t->left) {
	    depth++;
	    t = t->left;
	    continue;
	}
	
	/* Traverse right */
	if (t->right) {
	    depth++;
	    t = t->right;
	    continue;
	}
	    
	/* Backtrack */
	for (;t->parent;) {
	    /* from left, go to right */
	    if (t->parent->left == t && t->parent->right) {
		total_depth += depth;
		nodes++;
		t = t->parent->right;
		break;
	    }
	    
	    /* from left but no right, or from right */
	    total_depth += depth;
	    nodes++;
	    t = t->parent;
	    depth--;
	}
    }
    printf("total depth=%d, nodes=%d\n", total_depth, nodes);

    return total_depth;
}


/*
 * Rotate left the tree around node 't'. Returns the new root (the new 't').
 */
free_tree_n *tree_rotate_left(free_tree_n *t) {
    free_tree_n *r;

    r = t->right;
    t->right = r->left;
    t->right->parent = t;
    r->left = t;
    r->parent = t->parent;
    t->parent = r;

    return r;
}

/*
 * Rotate right the tree around node 't'. Returns the new root (the new 't').
 */
free_tree_n *tree_rotate_right(free_tree_n *t) {
    free_tree_n *l;

    l = t->left;
    t->left = l->right;
    t->left->parent = t;
    l->right = t;
    l->parent = t->parent;
    t->parent = l;

    return l;
}

/*
 *-----------------------------------------------------------------------------
 * "Freelist" management (historical name).
 * From here on things get a little backwards. For instance adding an item
 * to the free list can result in removing an item from the tree. Hence the
 * clear split between the tree code and the freelist code.
 *-----------------------------------------------------------------------------
 */


/*
 * Allocates a free tree, initialising it to a single free node starting
 * at pos, extending to pos+len.
 *
 * Returns a free_tree pointer for success
 *         NULL for failure
 */
free_tree *freetree_create(int pos, int len) {
    free_tree *t;

    if (NULL == (t = (free_tree *)xmalloc(sizeof(free_tree)))) {
	return NULL;
    }

    t->node_blocks = NULL;
    t->nblocks = 0;
    t->free_nodes = NULL;

    if (NULL == (t->tree = new_node(t))) {
	xfree(t);
	return NULL;
    }
    t->tree->pos = pos;
    t->tree->len = len;

    t->rover = NULL;
    t->largest = NULL;

    return t;
}


/*
 * Deallocates memory used by a free_tree
 */
void freetree_destroy(free_tree *t) {
    int i;

    if (!t)
	return;

    if (t->node_blocks) {
	for (i = 0; i < t->nblocks; i++) {
	    if (t->node_blocks[i])
		xfree(t->node_blocks[i]);
	}
    
	xfree(t->node_blocks);
    }

    xfree(t);
}


/*
 * Registers a (pos+len) pair with the free tree, effectively excluding that
 * space from the tree.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int freetree_register(free_tree *t, int pos, int len) {
    free_tree_n *n, *lnode;

    /* Find node for this position */
    n = tree_find_pos(t->tree, pos);

    /*
     * Our node found could be fitting one of the following cases
     * (order of checks matters).
     *
     * 1. pos == n->pos && len == n->len
     *    entire node is filled - so remove it.
     * 2. pos == n->pos
     *    overlaps start of node, update n->pos and n->len
     * 3. pos + len == n->pos + n->len
     *    overlaps end of node, update n->len
     * 4. else...
     *    contained within node - so split it.
     */
    if (pos == n->pos && len == n->len) {
	tree_delete_node(t, n);
    } else if (pos == n->pos) {
	n->len -= len;
	n->pos += len;
	assert(n->len > 0);
    } else if (pos + len == n->pos + n->len) {
	n->len -= len;
	assert(n->len > 0);
    } else {
	if (NULL == (lnode = new_node(t))) {
	    (void)gerr_set(TREE_OUT_OF_MEMORY);
	    return -1;
	}
	
	{
	    free_tree_n *tmp;

	    lnode->pos = n->pos;
	    lnode->len = pos - n->pos;
	    lnode->left = NULL;
	    lnode->right = NULL;
	    assert(lnode->pos >= 0);
	    assert(lnode->len > 0);

	    n->len = n->pos + n->len - (pos + len);
	    n->pos = pos + len;
	    assert(n->pos >= 0);
	    assert(n->len > 0);

	    if (n->left) {
		for (tmp = n->left; tmp->right; tmp=tmp->right)
		    ;

		tmp->right = lnode;
		lnode->parent = tmp;
	    } else {
		lnode->parent = n;
		n->left = lnode;
	    }
	}
    }

    return 0;
}

/*
 * Allocates a space from the free tree. If this doesn't exist, then space
 * from the first node (aka end of the disk file) is used.
 * 
 * Returns offset for success
 *        -1 for failure
 */
int freetree_allocate(free_tree *t, int len) {
    free_tree_n *n;
    int pos;

    /* Find free block */
    if (t->tree->left) {
	if (NULL == (n = tree_find_len(t, len))) {
	    /* None found in tree - use root node */
	    n = t->tree;
	    if (n->len < len) {
		(void)gerr_set(TREE_OUT_OF_SPACE);
		return -1;
	    }
	}
    } else {
	n = t->tree;
    }

    /*
     * Cases:
     * 1. entire block used - remove it.
     * 2. part of block used - use start.
     */
    pos = n->pos;
    if (len == n->len) {
	tree_delete_node(t, n);
    } else {
	if (n == t->largest)
	    t->largest = NULL;
	n->pos += len;
	n->len -= len;
	assert(n->len > 0);
    }

    /* printf("Allocated at position %d+%d\n", pos, len); */
    return pos;
}

/*
 * Reallocates an existing block to be a larger block. Speed tests show
 * that this doesn't seem to improve over a straight alloc+free system 
 * by any more than 2% or so. Additionally it actually seems to increase
 * (by a tiny amount) the size of files and number of free items. I can't
 * see why; perhaps just a quirk of the data being used by gap.
 *
 * Returns offset for success
 *        -1 for failure
 */
int freetree_reallocate(free_tree *t, int pos, int old_len, int new_len) {
    free_tree_n *n;
    int new_pos;

    /*
     * Search tree for the position at the end of this block. This will
     * tell us whether there is room to grow the block by.
     */
    n = tree_find_pos(t->tree, pos + old_len);
    if (n->pos == pos + old_len && n->pos + n->len >= pos + new_len) {

	/*
	 * We can extend this existing block. Two cases to handle.
	 * 1. extenstion covers entire free block.
	 * 2. part of free block used.
	 */
	if (n->pos + n->len == pos + new_len) {
	    tree_delete_node(t, n);
	    /* printf("Reallocated position %d by deleting node\n", pos); */
	} else {
	    int diff = pos + new_len - n->pos;

	    if (n == t->largest)
		t->largest = NULL;

	    n->pos += diff;
	    n->len -= diff;
	    assert(n->len > 0);
	    /* printf("Reallocated position %d by extending node\n", pos); */
	}
	new_pos = pos;

    } else {
	/*
	 * We can't extend, so simply reallocate and unregister.
	 */
	if (-1 != (new_pos = freetree_allocate(t, new_len)))
	    freetree_unregister(t, pos, old_len);
	/* printf("Reallocated position %d by alloc()+unreg()\n", pos); */
    }

    return new_pos;
}

/*
 * Adds some space (pos+len) back to the free_tree.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int freetree_unregister(free_tree *t, int pos, int len) {
    free_tree_n *n, *l, *r;
    int status;

    /*
     * Search tree for (pos,len) pair. This function will tell us whether
     * the (pos,len) pair butt up to a free block or is between blocks.
     */
    n = tree_find_pos_len(t->tree, pos, len, &status);
    switch(status) {
    default:
	/* printf("Dealloc %d+%d: OVERLAP\n", pos, len); */
	(void)gerr_set(TREE_OVERLAP);
	return -1;

    case 0:
	/*
	 * In between blocks, so create a new one. The node returned by
	 * tree_find_pos_len() will conveniently be the parent of the new
	 * node to create. We create either a right child or a left child
	 * depending on the relationship between pos and t->pos.
	 */
	if (NULL == (l = new_node(t))) {
	    (void)gerr_set(TREE_OUT_OF_MEMORY);
	    return -1;
	}
	l->left  = NULL;
	l->right = NULL;
	l->pos   = pos;
	l->len   = len;
	l->parent= n;

	if (!t->largest || l->len >= t->largest->len)
	    t->largest = l;

	if (n->pos > pos)
	    n->left = l;
	else
	    n->right = l;

	break;

    case 1:
	/*
	 * Butt with left end. We need to check that it doesn't also
	 * butt with the right end of another, in which case we'd like
	 * to join blocks.
	 *
	 * Find predecessor in tree to see whether this 
	 * requires joining two nodes. Actually this is impossible due to the
	 * design of the tree_find_pos_len() routine (it returns '2' in this
	 * case), but the code here shouldn't "know" the internal algorithm
	 * of the tree_find_pos_len() routine.
	 */
	l = tree_walk_left(n);

	if (l && l->pos + l->len == pos) {
	    /*
	     * To join nodes we merge the left into the right and
	     * delete the left. We don't need to worry about keeping
	     * position order within the tree as this must already be the
	     * case - we can never join two nodes when there's another node
	     * between them or we'd have detected an overlap.
	     */

	    n->len = n->pos + n->len - l->pos;
	    n->pos = l->pos;

	    if (!t->largest || n->len >= t->largest->len)
		t->largest = n;

	    assert(n->len > 0);
	    assert(n->pos >= 0);
	    tree_delete_node(t, l);
	} else {
	    /*
	     * Otherwise we simply update the len or pos+len pair to include
	     * the newly unregistered item.
	     */
	    n->pos -= len;
	    n->len += len;

	    if (!t->largest || n->len >= t->largest->len)
		t->largest = n;
	}

	break;

    case 2:
	/*
	 * Butt with right end. We need to check that it doesn't also
	 * butt with the left end of another, in which case we'd like
	 * to join blocks.
	 *
	 * Find successor in tree to see whether this 
	 * requires joining two nodes.
	 */
	r = tree_walk_right(n);

	if (r && r->pos == pos + len) {
	    /*
	     * To join nodes we merge the left into the right and
	     * delete the left. We don't need to worry about keeping
	     * position order within the tree as this must already be the
	     * case - we can never join two nodes when there's another node
	     * between them or we'd have detected an overlap.
	     */

	    r->len = r->pos + r->len - n->pos;
	    r->pos = n->pos;

	    if (!t->largest || r->len >= t->largest->len)
		t->largest = r;

	    assert(r->len > 0);
	    assert(r->pos >= 0);
	    tree_delete_node(t, n);
	} else {
	    /*
	     * Otherwise we simply update the len or pos+len pair to include
	     * the newly unregistered item.
	     */
	    n->len += len;

	    if (!t->largest || n->len >= t->largest->len)
		t->largest = n;
	}

	break;
    }

    /* printf("Unregistered at position %d+%d\n", pos, len); */
    return 0;
}

/*
 * ---------------------------------------------------------------------------
 * The code below is for saving and loading free_trees to disk.
 * This is used for fast database opening as it avoids the (expensive)
 * freetree_register function.
 */

static int find_node_addr(free_tree *t, free_tree_n *n) {
    int k;
    int ind;

    if (!n)
	return -1;

    ind = -1;
    for (k = 0; k < t->nblocks; k++) {
	if (n >= &t->node_blocks[k][0] &&
	    n <  &t->node_blocks[k][BLOCK_FACTOR]) {
	    ind = n - t->node_blocks[k];
	    break;
	}
    }
    if (-1 == ind) {
	fprintf(stderr, "Fatal: couldn't save tree\n");
	return -1;
    }

    return ind + k * BLOCK_FACTOR;
}

static int node_ind2addr(free_tree *t, int ind, free_tree_n **nptr) {
    int major;
    int minor;

    if (ind == -1)
	return -1;

    major = ind / BLOCK_FACTOR;
    minor = ind % BLOCK_FACTOR;

    if (major < 0 || major >= t->nblocks)
	return -1;

    if (minor < 0)
	return -1;

    *nptr = &t->node_blocks[major][minor];
    return 0;
}

/*
 * Calculates a 32-bit CRC for a block of data.
 */
#define POLY 0xedb88320
static int calc_crc(char *data, int len) {
    unsigned int table[256], crc, i, j;

    for (i = 0; i < 256; i++) {
        crc = i;
        for (j = 8; j > 0; j--)
            crc = (crc & 1) ? (crc >> 1) ^ POLY : (crc >> 1);
        table[i] = crc;
    }

    crc = 0xffffffff;
    for (i = 0; i < (unsigned int)len; i++)
        crc = (crc >> 8) ^ table[(crc ^ data[i]) & 0xff];

    return crc;
}

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
int freetree_save(free_tree *t, int fd, int last_time) {
    int i, j;
    free_tree_n *n;
    int4 *data;
    int data_ind = 0;
    int endian = 1;

    /* Allocate a block of data to write our tree into */
    data = (int4 *)xmalloc((t->nblocks * BLOCK_FACTOR * 5 + 6) * sizeof(int4));
    if (!data)
	return -1;

    /* Main tree structure */
    data[data_ind++] = t->nblocks;
    data[data_ind++] = last_time;
    data[data_ind++] = find_node_addr(t, t->free_nodes);
    data[data_ind++] = find_node_addr(t, t->tree);
    data[data_ind++] = find_node_addr(t, t->rover);

    for (i = 0; i < t->nblocks; i++) {
	for (j = 0; j < BLOCK_FACTOR; j++) {
	    n = &t->node_blocks[i][j];

	    data[data_ind++] = find_node_addr(t, n->left);
	    data[data_ind++] = find_node_addr(t, n->right);
	    data[data_ind++] = find_node_addr(t, n->parent);
	    data[data_ind++] = n->pos;
	    data[data_ind++] = n->len;
	}
    }

    if (*(char *)&endian) {
	int i;
	for (i = 0; i < data_ind; i++)
	    swap_int4(data[i], data[i]);
    }

    data[data_ind] = calc_crc((char *)data, data_ind * 4);
    data_ind++;

    if (*(char *)&endian) {
	swap_int4(data[data_ind-1], data[data_ind-1]);
    }

    /* LOW LEVEL IO HERE */
    write(fd, data, data_ind * sizeof(int4));

    xfree(data);

    return 0;
}

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
free_tree *freetree_load(int fd, int last_time) {
    free_tree *t = NULL;
    int4 *data = NULL;
    int4 nblocks;
    int endian = 1;
    int sz;
    int data_ind = 0;
    int last_block;
    int i, j;

    /* Load the first 5 blocks */
    sz = 5 * sizeof(int4);
    if (!(data = (int4 *)xmalloc(sz)))
	goto error;

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (20 != read(fd, data, sz))
	goto error;

    /* Check time stamp */
    if (be_int4(data[1]) != last_time) {
	fprintf(stderr, "Incorrect tree timestamp - rebuilding tree\n");
	goto error;
    }

    /* Load the other blocks */
    nblocks = be_int4(data[0]);
    last_block = nblocks * BLOCK_FACTOR * 5 + 5;
    sz =  (last_block+1) * 4;
    if (!(data = (int4 *)xrealloc(data, sz)))
	goto error;

    sz -= 5*4; /* first 5 already loaded */
    
    /* LOW LEVEL IO HERE */
    errno = 0;
    if (sz != read(fd, &data[5], sz)) {
	fprintf(stderr, "Tree too short\n");
	goto error;
    }


    /* Check CRC */
    if (*(char *)&endian) {
	swap_int4(data[last_block], data[last_block]);
    }
    if (data[last_block] != calc_crc((char *)data, last_block*4)) {
	fprintf(stderr, "Invalid tree CRC\n");
	goto error;
    }

    /* Byte-swap blocks */
    if (*(char *)&endian) {
	int i;

	for (i = 0; i < last_block ; i++)
	    swap_int4(data[i], data[i]);
    }

    /* Construct the tree top */
    if (!(t = (free_tree *)xmalloc(sizeof(free_tree))))
	goto error;
    t->tree = NULL;
    t->rover = NULL;
    t->node_blocks = NULL;
    t->nblocks = 0;
    t->free_nodes = NULL;
    t->largest = NULL;

    /* Allocate tree memory */
    if (!(t->node_blocks = (free_tree_n **)xmalloc(nblocks *
						   sizeof(free_tree_n *))))
	goto error;
    for (i = 0; i < nblocks; i++)
	t->node_blocks[i] = NULL;
    for (i = 0; i < nblocks; i++) {
	free_tree_n *n;

	n = t->node_blocks[i] =
	    (free_tree_n *)xmalloc(BLOCK_FACTOR * sizeof(free_tree_n));
	if (!n)
	    goto error;
	for (j = 0; j < BLOCK_FACTOR; j++) {
	    n[j].left = NULL;
	    n[j].right = NULL;
	    n[j].parent = NULL;
	    n[j].pos = 0;
	    n[j].len = 0;
	}
    }

    /* Decode blocks into tree */
    t->nblocks    = nblocks;
    data_ind = 2;
    if (-1 == node_ind2addr(t, data[data_ind++], &t->free_nodes))
	goto error;

    if (-1 == node_ind2addr(t, data[data_ind++], &t->tree))
	goto error;

    if (-1 == node_ind2addr(t, data[data_ind++], &t->rover)) {
	t->rover = NULL;
    }

    for (i = 0; i < t->nblocks; i++) {
	for (j = 0; j < BLOCK_FACTOR; j++) {
	    free_tree_n *n = &t->node_blocks[i][j];
	    if (-1 ==  node_ind2addr(t, data[data_ind++], &n->left))
		n->left = NULL;

	    if (-1 == node_ind2addr(t, data[data_ind++], &n->right))
		n->right = NULL;

	    if (-1 == node_ind2addr(t, data[data_ind++], &n->parent))
		n->right = NULL;

	    n->pos = data[data_ind++];
	    n->len = data[data_ind++];
	}
    }

    if (data)
	xfree(data);

    return t;

 error:
    if (data)
	xfree(data);
    if (t)
	freetree_destroy(t);
    return NULL;
}
