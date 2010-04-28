#include <staden_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
// #define NDEBUG
#include <assert.h>
#include "freetree.h"
#include "xalloc.h"
#include "g-os.h"
#include "misc.h"

/*
 * Apr 2005: This code now uses AVL trees for maintaining a mostly balanced
 * freetree. Implementing AVL balancing isn't trivial, but the following
 * web resources helped greatly.
 *
 * GNU libavl docs:
 *    http://www.stanford.edu/~blp/avl/libavl.html/Rebalancing-AVL-Trees.html
 *
 * A C++ AVL tree (less detailed, but still with some useful hints):
 *    http://www.cmcrossroads.com/bradapp/ftp/src/libs/C++/AvlTrees.html
 */

/*
 * TODO:
 * - Add check to tree_find_pos to error on cases where the block isn't there.
 *   This will spot the 'overlapping image' cases, I hope.
 *
 * - Use hints from various malloc systems for improved performance.
 *   "Dynamic Storage Allocation: A Survey and Critical Review" is a good
 *   review article on such things.
 *
 * - Test use of rovers with the segregated freelist. This is "next fit" which
 *   may turn out to cause fragmentation, but in theory it should be faster.
 */

/* Turns size 'l' into an index to the free_tree->lists array. */
#define SEGLISTNUM(t, l) (((l)<=FBLOCKM) \
    ? (t)->listsz[((l)/FBSZ)] \
    : (ilog2((l))-12+NFBLOCK))

/*
 * Returns the integer component of log2(len).
 */
static int ilog2(GImage len) {
    int l = 0;
    while (len >>= 1) l++;
    return l;
}



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
	    n[i].prev = NULL;
	    n[i].next = NULL;
	    n[i].pos = 0;
	    n[i].len = 0;
	    n[i].balance = 0;
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
    n->balance = 0;
    n->prev = NULL;
    n->next = NULL;

    return n;
}

static void del_node(free_tree *t, free_tree_n *n) {
    n->left = t->free_nodes;
    n->right = NULL;
    n->parent = NULL;
    n->prev = NULL;
    n->next = NULL;
    t->free_nodes = n;
}

/*
 * Search for a tree element in which pos is a member of. We scan the tree
 * in an inorder fasion.
 *
 * Complexity: O(log(N)) for a balanced tree.
 *
 * Always returns a free_tree_n pointer.
 */
static free_tree_n *tree_find_pos(free_tree *ft, GImage pos) {
    free_tree_n *t = ft->root;

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
 * Complexity: O(log(N)) for a balanced tree.
 *
 * Returns the node next to (pos,len) for success
 *         NULL for not found.
 *
 * status is filled with -1 for error (overlap found)
 *                        0 for no butt (inbetween free blocks)
 *                        1 for butt with left end
 *                        2 for butt with right end
 */
static free_tree_n *tree_find_pos_len(free_tree *ft, GImage pos,
				      GImage len, int *status) {
    GImage end;
    free_tree_n *t = ft->root;

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


static free_tree_n *tree_find_len(free_tree *t, GImage len) {
    int sz;
    free_tree_n *n;
    /* int count = 0; */

    /*
     * Loop through the free lists, starting with optimal size and working up.
     */
    for (sz = SEGLISTNUM(t, len); sz < NFBLOCK2; sz++) {
	/* count++; */
	for (n = t->lists[sz]; n; n = n->next) {
	    if (n->len >= len) {
		/* printf("tree_find_len %"PRIGImage" => %"PRIGImage
		          " (%d searched)\n", len, n->len, count);*/
		return n;
	    }
	    /* count++; */
	}
    }

    /*
    printf("tree_find_len %"PRIGImage" => %"PRIGImage
	   " (%d searched, got wilderness)\n",
	   len, t->wilderness->len, count);
    */

    return t->wilderness->len >= len ? t->wilderness : NULL;
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


static GImage last_end;
static GImage tree_sum;

static int tree_height(free_tree_n *t) {
    int lh = 0, rh = 0;

    if (!t) return 0;

    lh = 1 + tree_height(t->left);
    rh = 1 + tree_height(t->right);

    return MAX(lh, rh);
}

/*
 * Prints a tree to stdout in the "inorder" order.
 */
static int tree_print_r(free_tree_n *t) {
    int hl, hr;
    int err = 0;
    if (t->left) {
	assert(t->left->parent == t);
	err |= tree_print_r(t->left);
    }
    hl = tree_height(t->left);
    hr = tree_height(t->right);
    printf("%.*s%p: l=%p r=%p p=%p hgt=%d/%d/%d %"
	   PRIGImage"+%"PRIGImage"\n",
	   ABS(t->balance), t->balance+"--=+++"+2,
	   t, t->left, t->right, t->parent,
	   tree_height(t), hr - hl, t->balance,
	   t->pos, t->len);

    err |= (hr-hl != t->balance);

    tree_sum += t->pos;
    assert(t->pos > last_end);
    assert(t->len > 0);
    last_end = t->pos + t->len;
    if (t->right) {
	assert(t->right->parent == t);
	err |= tree_print_r(t->right);
    }
    return err;
}

void tree_print_lists(free_tree *t) {
    int i;
    free_tree_n *n;

    printf("============== FREETREE THREADED LISTS ============\n");
    for (i = 0; i < NFBLOCK2; i++) {
	int count;
	printf(">>> list %d\n", i);
	count = 0;
	for (n = t->lists[i]; n; n = n->next) {
	    printf("%"PRIGImage" ", n->len);
	    assert(n->len > 0);
	    if (n->prev == NULL)
		assert(n == t->lists[i]);
	    if (n->prev)
		assert(n->prev->next == n);
	    if (n->next)
		assert(n->next->prev == n);
	    count++;
	}
	printf(" => %d entries\n", count);
    }
}

void tree_print(free_tree *t) {
    int err;
    printf("============== TREE %p, root %p ============\n", t, t->root);
    last_end = -1;
    tree_sum = 0;
    err = tree_print_r(t->root);
    if (err) {
	printf("@@@@@@ INVALID BALANCES @@@@@@\n");
    }
    assert(err == 0);
    printf("Tree sum=%"PRIGImage"\n", tree_sum);
}

static void tree_print_dot_r(FILE *fp, free_tree_n *t) {
    if (t->left)
	fprintf(fp, "edge [color=\"#00BB00\"] a%p -> a%p\n", t, t->left);
    if (t->right)
	fprintf(fp, "edge [color=\"#2020FF\"] a%p -> a%p\n", t, t->right);

    if (t->left) {
	tree_print_dot_r(fp, t->left);
    }

    if (t->right) {
	tree_print_dot_r(fp, t->right);
    }
}

void tree_print_dot(free_tree *t) {
    FILE *fp = fopen("freetree.dot", "w");
    fprintf(fp, "digraph G {\n");
    fprintf(fp, "rotate=90 size=\"11,3\" ratio=fill\n");
    tree_print_dot_r(fp, t->root);
    fprintf(fp, "}\n");
    fclose(fp);
}

/*
 * Generates postscript to print the tree. It's pretty rough and ready, with
 * constant scale factors (hand edit the postscript to fit!), but serves
 * well for quick debugging and analysis of how the tree is balanced.
 */
void tree_print_ps(free_tree_n *t) {
    float step=65536;
    float yinc = 10000;
    float yred = 0.98;
    float xred = 0.7;
    int i;
    int depth = 0, maxdepth = 0;

    puts("%!PS");
    puts("500 380 translate 90 rotate newpath 0 0 moveto .003 .0005 scale");

    for (; t->parent;) {
	if (!t->left) {
	    if (!t->right) {
		for (;t->parent; step/=xred, yinc/=yred) {
		    if (t->parent->left == t && t->parent->right) {
			t = t->parent->right;
			printf("%f %f rmoveto\n", step/xred, -yinc/yred);
			printf("%f %f rlineto\n", step/xred, yinc/yred);
			break;
		    }
		    if (t->parent->left == t)
			printf("%f %f rmoveto\n", step/xred, -yinc/yred);
		    else
			printf("%f %f rmoveto\n", -step/xred, -yinc/yred);
		    t=t->parent;
		    depth--;
		}
	    } else {
		printf("%f %f rlineto\n", step, yinc);
		step*=xred;
		yinc*=yred;
		depth++;
		t = t->right;
	    }
	} else {
	    printf("%f %f rlineto\n", -step, yinc);
	    step*=xred;
	    yinc*=yred;
	    depth++;
	    t = t->left;
	}
	if (maxdepth < depth)
	    maxdepth = depth;
    }
    puts("stroke");
    {
	double ypos = 10000;
	for (yinc=10000,i=0; i < maxdepth; i++) {
	    printf("-100000 %f moveto 100000 %f lineto\n", ypos, ypos);
	    yinc *= yred;
	    ypos += yinc;
	}
    }
    
    puts("stroke showpage");
}


/*
 * Rotate left the tree around node 'A'. Returns the new root ('B').
 *
 *    A              B
 *   / \            / \
 *  1   B    =>    A   3
 *     / \        / \
 *    2   3      1   2
 */
free_tree_n *tree_rotate_left(free_tree_n *A) {
    free_tree_n *B = A->right;

    B->parent = A->parent;
    if ((A->right = B->left))
	A->right->parent = A;
    B->left = A;
    A->parent = B;

    A->balance = - --B->balance;

    return B;
}

/*
 * Rotate right the tree around node 'A'. Returns the new root ('B').
 *
 *      A          B
 *     / \        / \
 *    B   3  =>  1   A
 *   / \            / \
 *  1   2          2   3
 */
free_tree_n *tree_rotate_right(free_tree_n *A) {
    free_tree_n *B = A->left;

    B->parent = A->parent;
    if ((A->left = B->right))
	A->left->parent = A;
    B->right = A;
    A->parent = B;

    A->balance = - ++B->balance;
    return B;
}

/*
 * Rotate double-left around node 'A'. Returns the new root ('B').
 *
 *         A              B
 *        / \           /   \
 *       1   C         A     C
 *          / \   =>  / \   / \
 *         B   4     1   2 3   4
 *        / \
 *       2   3
 */
free_tree_n *tree_rotate_left2(free_tree_n *A) {
    free_tree_n *C = A->right;
    free_tree_n *B = C->left;

    B->parent = A->parent;
    if ((A->right = B->left))
	A->right->parent = A;
    if ((C->left = B->right))
	C->left->parent = C;
    B->left = A;
    A->parent = B;
    B->right = C;
    C->parent = B;

    B->left->balance  = -MAX(B->balance, 0);
    B->right->balance = -MIN(B->balance, 0);
    B->balance = 0;

    return B;
}

/*
 * Rotate double-right around node 'A'. Returns the new root (the new 'B').
 *
 *         A              B
 *        / \           /   \
 *       C   4         C     A
 *      / \       =>  / \   / \
 *     1   B         1   2 3   4
 *        / \
 *       2   3
 */
free_tree_n *tree_rotate_right2(free_tree_n *A) {
    free_tree_n *C = A->left;
    free_tree_n *B = C->right;

    B->parent = A->parent;
    if ((C->right = B->left))
	C->right->parent = C;
    if ((A->left = B->right))
	A->left->parent = A;
    B->left = C;
    C->parent = B;
    B->right = A;
    A->parent = B;

    B->left->balance  = -MAX(B->balance, 0);
    B->right->balance = -MIN(B->balance, 0);
    B->balance = 0;

    return B;
}

/*
 * Rebalance the (sub)tree rooted at node 'n' if required. Called after
 * an insertion.
 */
void tree_rebalance(free_tree *t, free_tree_n *n) {
    free_tree_n *parent = n->parent;
    free_tree_n *new_root = NULL;

    switch (n->balance) {
    case -2:
	switch (n->left->balance) {
	case -1:
	    new_root = tree_rotate_right(n);
	    break;
	case +1:
	    new_root = tree_rotate_right2(n);
	    break;
	default:
	    abort();
	}
	break;

    case +2:
	switch (n->right->balance) {
	case +1:
	    new_root = tree_rotate_left(n);
	    break;
	case -1:
	    new_root = tree_rotate_left2(n);
	    break;
	default:
	    abort();
	}
	break;
    }

    if (!new_root)
	return;

    if (parent) {
	if (parent->left == n)
	    parent->left = new_root;
	else
	    parent->right = new_root;
    } else {
	t->root = new_root;
    }
}

/*
 * Insert node 'n' to the appropriate linked list blocked by n->len.
 * It does this by inserting to the head of the list, as this means that reuse
 * of blocks can be made to try and reuse blocks most recently freed (so that
 * disk caching is optimised).
 */
static void list_insert_node(free_tree *t, free_tree_n *n) {
    int sz = SEGLISTNUM(t, n->len);
    assert(sz >= 0 && sz < NFBLOCK2);
    if ((n->next = t->lists[sz]))
	n->next->prev = n;
    n->prev = NULL;
    t->lists[sz] = n;
}

/*
 * Deletes node 'n' from the linked lists, updating the master list pointer
 * and linking left/right neighbours as appropriate.
 */
static void list_delete_node(free_tree *t, free_tree_n *n) {
    int sz = SEGLISTNUM(t, n->len);
    assert(sz >= 0 && sz < NFBLOCK2);

    if (n->prev)
	n->prev->next = n->next;
    if (n->next)
	n->next->prev = n->prev;
    if (t->lists[sz] == n)
	t->lists[sz] = n->next;
    n->next = n->prev = NULL;
}

/*
 * Handles resizing of node 'n' from size 'oldsz' to 'newsz' with regards to
 * the linked list.
 * This does NOT actually resize the node, but instead it makes sure that it
 * is in the correct linked list.
 */
static void list_resize_node(free_tree *t, free_tree_n *n,
			     GImage oldsz, GImage newsz) {
    int bo = SEGLISTNUM(t, oldsz);
    int bn = SEGLISTNUM(t, newsz);

    assert(bo >= 0 && bo < NFBLOCK2);
    assert(bn >= 0 && bn < NFBLOCK2);

    /* Nothing to do? */
    if (bo == bn)
	return;

    /* Otherwise remove from one list and insert to the other */
    /* Remove */
    if (n->prev)
	n->prev->next = n->next;
    if (n->next)
	n->next->prev = n->prev;
    if (t->lists[bo] == n)
	t->lists[bo] = n->next;

    /* Reinsert */
    if ((n->next = t->lists[bn]))
	n->next->prev = n;
    n->prev = NULL;
    t->lists[bn] = n;
}

/*
 * Inserts 'node' to be a child of 'parent' to either the left (direction = -1)
 * or right (direction = +1) child.
 */
static void tree_insert_node(free_tree *t,
			     free_tree_n *parent,
			     free_tree_n *node,
			     int direction) {
    /*
    printf("Insert node %p to %p dir %d\n",
	   node, parent, direction);
    */

    /* Update the threaded linked-lists */
    list_insert_node(t, node);

    /* Link it in */
    if (direction == -1) {
	assert(parent->left == NULL);
	parent->left = node;
	node->parent = parent;
    } else {
	assert(parent->right == NULL);
	parent->right = node;
	node->parent = parent;
    }

    /* Correct balancing */
    /* Propogate change upwards when parent's balance is zero */
    while (parent && parent->balance == 0) {
	parent->balance += (parent->left == node) ? -1 : +1;
	node = parent;
	parent = node->parent;
    }
    if (!parent)
	return;
    /*
     * Now we have the other two cases, where parent's balance is either
     * the opposite of direction (in which case we are balancing it) or
     * the same (in which case it becomes too unbalanced and rotations are
     * requried).
     * 
     * So we apply the balance and if it's -2 or +2 we rebalance.
     */
    parent->balance += (parent->left == node) ? -1 : +1;

    tree_rebalance(t, parent);
}

/*
 * Deletes a node, maintaining AVL balance factors.
 */
/*static*/ void tree_delete_node(free_tree *t, free_tree_n *n) {
    free_tree_n *p = n->parent, *r = n->right, *s;
    int dir, sdir;

    /* printf("Delete_node %p\n", n); */

    /* Update the threaded linked-lists */
    list_delete_node(t, n);

    if (p) {
	dir = (p->left == n) ? -1 : +1;
    } else {
	dir = 0;
    }

    if (r == NULL) {
	/* No right child, so move n->left to the place of n */
	if (dir == -1) {
	    if ((p->left = n->left))
		p->left->parent = p;
	} else if (dir == +1) {
	    if ((p->right = n->left))
		p->right->parent = p;
	} else {
	    t->root = n->left;
	    n->left->parent = NULL;
	}

	s = p;
	sdir = dir;

    } else if (r->left == NULL) {
	/* We have a right child 'r', but it does not have a left child */
	if ((r->left = n->left)) {
	    r->left->parent = r;
	}
	r->parent = p;
	if (dir == -1) {
	    p->left = r;
	} else if (dir == +1) {
	    p->right = r;
	} else {
	    t->root = r;
	}

	r->balance = n->balance;
	s = r;
	sdir = +1;

    } else {
	/*
	 * We have a right child 'r' which also has a left child. Here we do
	 * an inorder successor 's' to find the node that should replace 'n'
	 */
	for (s = r->left; s->left; r = s, s = r->left);

	/* Move s */
	if ((s->left = n->left))
	    s->left->parent = s;
	if ((r->left = s->right))
	    r->left->parent = r;
	s->right = n->right;
	s->right->parent = s;
	s->parent = p;
	if (dir == -1) {
	    p->left = s;
	} else if (dir == +1) {
	    p->right = s;
	} else {
	    t->root = s;
	}

	s->balance = n->balance;
	s = r; /* first node to backtrack from when balancing */
	sdir = -1;
    }


    /*
     * Now rebalance, starting from 's'and working upwards.
     * 'sdir' holds the direction from 's' that we came from (or initially
     * from where the node was deleted).
     */
    while (s) {
	if (sdir == -1) {
	    s->balance++;
	    /*
	     * Case 1: deleting the node balances the tree at 's',
	     * but the change in height means theparent may need updating
	     * (so we keep looping).
	     */
	    /* s->balance == 0; do nothing */

	    /*
	     * Case 2: deleting the node unbalanced the tree at 's',
	     * but the total height of s will not change
	     * (was h/h, now h-1/h, so still max height of h).
	     * => update and break
	     */
	    if (s->balance == +1) {
		break;

	    /*
	     * Case 3: deleting the node causes the already unbalanced tree
	     * at 's' to become overly unbalanced, requiring a rotation.
	     */
	    } else if (s->balance > +1) {
		free_tree_n *x = s->right;
		free_tree_n *p = s->parent;
		int bal = x->balance;
		if (bal < 0) {
		    r = tree_rotate_left2(s);
		} else {
		    r = tree_rotate_left(s);
		}

		if (p) {
		    if (p->left == s)
			p->left = r;
		    else
			p->right = r;
		} else {
		    t->root = r;
		}

		if (bal == 0)
		    break;

		s = r;
	    }
	} else if (sdir == +1) {
	    s->balance--;
	    if (s->balance == -1) {
		break;
	    } else if (s->balance < -1) {
		free_tree_n *x = s->left;
		free_tree_n *p = s->parent;
		int bal = x->balance;
		if (bal > 0) {
		    r = tree_rotate_right2(s);
		} else {
		    r = tree_rotate_right(s);
		}

		if (p) {
		    if (p->left == s)
			p->left = r;
		    else
			p->right = r;
		} else {
		    t->root = r;
		}

		if (bal == 0)
		    break;

		s = r;
	    }
	}

	if (s->parent)
	    sdir = (s->parent->left == s) ? -1 : +1;
	s = s->parent;
    }

    /* Deallocate the deleted node */
    del_node(t, n);
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
free_tree *freetree_create(GImage pos, GImage len) {
    free_tree *t;
    int i, j;

    if (NULL == (t = (free_tree *)xmalloc(sizeof(free_tree)))) {
	return NULL;
    }

    t->node_blocks = NULL;
    t->nblocks = 0;
    t->free_nodes = NULL;

    if (NULL == (t->root = new_node(t))) {
	xfree(t);
	return NULL;
    }
    t->root->pos = pos;
    t->root->len = len;

    t->wilderness = t->root;

    for (i = 0; i < NFBLOCK2; i++) {
	t->lists[i] = NULL;
    }
    /*
     * For the first 4Kb, we use direct indexed lookups of all possible
     * sizes of blocks. After that we use powers of 2.
     */
    for (i = j = 0; i <= 256/FBSZ; j++) {
	t->listsz[i++] = j;
    }
    for (; i <= 512/FBSZ; j++) {
	t->listsz[i++] = j;
	t->listsz[i++] = j;
    }

    for (; i <= 1024/FBSZ; j++) {
	t->listsz[i++] = j;
	t->listsz[i++] = j;
	t->listsz[i++] = j;
	t->listsz[i++] = j;
    }

    for (; i <= FBLOCKM/FBSZ; j++) {
	t->listsz[i++] = j;
	t->listsz[i++] = j;
	t->listsz[i++] = j;
	t->listsz[i++] = j;
	t->listsz[i++] = j;
	t->listsz[i++] = j;
	t->listsz[i++] = j;
	t->listsz[i++] = j;
    }

    /*
    for (i = 0; i <= FBLOCKM/FBSZ; i++) {
	t->listsz[i] = i;
    }
    */

    return t;
}


/*
 * Deallocates memory used by a free_tree
 */
void freetree_destroy(free_tree *t) {
    int i;

    /*
      tree_print(t);
      tree_print_lists(t);
    */

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
int freetree_register(free_tree *t, GImage pos, GImage len) {
    free_tree_n *n, *lnode;

    /* Find node for this position */
    n = tree_find_pos(t, pos);

    /*
    printf("Register %ld+%ld => found %p (%ld+%ld)\n",
	   (long)pos, (long)len, n, (long)n->pos, (long)n->len);
    */

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
	list_resize_node(t, n, n->len, n->len-len);
	n->len -= len;
	n->pos += len;
	assert(n->len > 0);
    } else if (pos + len == n->pos + n->len) {
	list_resize_node(t, n, n->len, n->len-len);
	n->len -= len;
	assert(n->len > 0);
    } else {
	free_tree_n *tmp;
	GImage old_size;

	if (NULL == (lnode = new_node(t))) {
	    (void)gerr_set(TREE_OUT_OF_MEMORY);
	    return -1;
	}

	/* Create new node (will be left hand end) */
	lnode->pos = n->pos;
	lnode->len = pos - n->pos;
	lnode->left = NULL;
	lnode->right = NULL;
	assert(lnode->pos >= 0);
	assert(lnode->len > 0);

	/* Shrink old node (right hand end) */
	old_size = n->len;
	n->len = n->pos + n->len - (pos + len);
	n->pos = pos + len;
	assert(n->pos >= 0);
	assert(n->len > 0);
	list_resize_node(t, n, old_size, n->len);

	/* Insert to the tree */
	if (n->left) {
	    for (tmp = n->left; tmp->right; tmp=tmp->right)
		;
	    
	    tree_insert_node(t, tmp, lnode, +1);
	} else {
	    tree_insert_node(t, n, lnode, -1);
	}
    }

    return 0;
}

/*
 * Allocates a space from the free tree. If this doesn't exist, then space
 * from the wildeness node (aka end of the disk file) is used.
 * 
 * Returns offset for success
 *        -1 for failure
 */
GImage freetree_allocate(free_tree *t, GImage len) {
    free_tree_n *n;
    GImage pos;

    /* Find free block */
    if (NULL == (n = tree_find_len(t, len))) {
	/* None found in tree - use root node */
	n = t->wilderness;
	if (n->len < len) {
	    (void)gerr_set(TREE_OUT_OF_SPACE);
	    return -1;
	}
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
	list_resize_node(t, n, n->len, n->len - len);
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
int freetree_reallocate(free_tree *t, GImage pos, GImage old_len,
			GImage new_len) {
    free_tree_n *n;
    GImage new_pos;

    /*
     * Search tree for the position at the end of this block. This will
     * tell us whether there is room to grow the block by.
     */
    n = tree_find_pos(t, pos + old_len);
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
	    GImage diff = pos + new_len - n->pos;

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
int freetree_unregister(free_tree *t, GImage pos, GImage len) {
    free_tree_n *n, *l, *r;
    int status;

    /*
     * Search tree for (pos,len) pair. This function will tell us whether
     * the (pos,len) pair butt up to a free block or is between blocks.
     */
    n = tree_find_pos_len(t, pos, len, &status);
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
	assert(l->len > 0);

	tree_insert_node(t, n, l, n->pos > pos ? -1 : +1);

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
	    GImage old_len;

	    /*
	     * To join nodes we merge the left into the right and
	     * delete the left. We don't need to worry about keeping
	     * position order within the tree as this must already be the
	     * case - we can never join two nodes when there's another node
	     * between them or we'd have detected an overlap.
	     */

	    old_len = n->len;
	    n->len = n->pos + n->len - l->pos;
	    n->pos = l->pos;
	    list_resize_node(t, n, old_len, n->len);

	    assert(n->len > 0);
	    assert(n->pos >= 0);
	    tree_delete_node(t, l);
	} else {
	    /*
	     * Otherwise we simply update the len or pos+len pair to include
	     * the newly unregistered item.
	     */
	    list_resize_node(t, n, n->len, n->len + len);
	    n->pos -= len;
	    n->len += len;
	    assert(n->len > 0);
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
	    GImage old_len;
	    /*
	     * To join nodes we merge the left into the right and
	     * delete the left. We don't need to worry about keeping
	     * position order within the tree as this must already be the
	     * case - we can never join two nodes when there's another node
	     * between them or we'd have detected an overlap.
	     */

	    old_len = r->len;
	    r->len = r->pos + r->len - n->pos;
	    r->pos = n->pos;
	    list_resize_node(t, r, old_len, r->len);

	    assert(r->len > 0);
	    assert(r->pos >= 0);
	    tree_delete_node(t, n);
	} else {
	    /*
	     * Otherwise we simply update the len or pos+len pair to include
	     * the newly unregistered item.
	     */
	    list_resize_node(t, n, n->len, n->len + len);
	    n->len += len;
	    assert(n->len > 0);
	}

	break;
    }

    /* printf("Unregistered at position %d+%d\n", pos, len); */
    return 0;
}

/*
 * NOTE: For now the code to read/write the cached freetree has been disabled.
 * This is for several main reasons.
 *
 * 1. Old freetrees will be unbalanced, so loading from scratch with the AVL
 *    balancing will improve performance.
 *
 * 2. Benchmarking shows that using the AVL code to initiate a new tree from
 *    scratch is very close to the performance of loading a cached tree.
 *
 * 3. The format of the tree would need updating too, to store the balance
 *    factors and the segregated lists.
 */

#if CACHE_FREETREE
/*
 * ---------------------------------------------------------------------------
 * The code below is for saving and loading free_trees to disk.
 * This is used for fast database opening as it avoids the (expensive)
 * freetree_register function.
 */

#define NUM(n) (n ? (int)find_node_addr(g,n) : 0)

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

#define TYPE int4
#include "freetree-io.h"
#undef TYPE

#define TYPE int8
#include "freetree-io.h"
#undef TYPE

free_tree *freetree_load_int4(int fd, int last_time) {
    return NULL;
}
free_tree *freetree_load_int8(int fd, int last_time) {
    return NULL;
}
int freetree_save_int4(free_tree *t, int fd, int last_time) {
    return -1;
}
int freetree_save_int8(free_tree *t, int fd, int last_time) {
    return -1;
}
#endif
