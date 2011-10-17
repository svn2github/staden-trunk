/*
 * As b+tree.c but we now do all tree accesses via method calls to a tree
 * record number, instead of holding direct pointers in the btree_node_t.
 *
 * This is more analogous to how it'd work in practice on disk.
 */

#include <staden_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include "b+tree2.h"
#include "tg_utils.h"

char *btree_check(btree_t *t, btree_node_t *n, char *pleaf);

btree_node_t *btree_new_node(void) {
    btree_node_t *n = malloc(sizeof(*n));
    int i;
    for (i = 0; i <= BTREE_MAX; i++) {
	n->keys[i] = NULL;
	n->chld[i] = 0;
    }
    n->leaf = 1;
    n->used = 0;
    n->parent = 0;
    n->next = 0;
    n->cache = NULL;

    return n;
}

void btree_del_node(btree_node_t *n) {
    int i;

    for (i = 0; i < n->used; i++) {
	if (n->keys[i])
	    free(n->keys[i]);
    }
    free(n);
}

btree_t *btree_new(void *cd, BTRec root) {
    btree_t *t = malloc(sizeof(*t));
    t->cd = cd;
    if (root)
	t->root = btree_node_get(t->cd, root);
    else
	t->root = btree_node_new(t->cd);

    if (!t->root) {
	free(t);
	return NULL;
    }
    
    btree_inc_ref(cd, t->root);

    return t;
}

static void btree_del_recurse(btree_t *t, btree_node_t *n) {
    int i;

    for (i = 0; i < n->used; i++) {
	if (!n->leaf && n->chld[i]) {
	    btree_del_recurse(t, btree_node_get(t->cd, n->chld[i]));
	}
    }
    btree_node_del(t->cd, n);
}

void btree_del(btree_t *t) {
    btree_del_recurse(t, t->root);
    free(t);
}


static int btree_find_key(btree_node_t *n, char *str) {
    int low, high, mid;

#if 0
    int i;
    /* Simple brute force search */
    for (i = 0; i < n->used; i++) {
	if (strcmp(n->keys[i], str) >= 0)
	    break;
    }
    return i;
#endif

    /* Binary search between low..high) */
    low = 0; high = n->used-1; mid = 0;
    while (low < high) {
	mid = (high-low+1)/2 + low;
	if (strcmp(n->keys[mid], str) >= 0)
	    high = mid-1;
	else
	    low = mid+1;
    }

    if (n->keys[low] && strcmp(n->keys[low], str) < 0)
	low++;

    return low;
}

/* Set parent index 'ind' to have key 'str' */
static int btree_set_parent_key(btree_t *t, btree_node_t *p,
				int ind, char *str, int recurse) {
    btree_inc_ref(t->cd, p);

    if (recurse && p->parent && p->keys[ind]) {
	btree_node_t *par = btree_node_get(t->cd, p->parent);
	int ind2 = btree_find_key(par, p->keys[ind]);
	if (par->keys[ind2] &&
	    0 == strcmp(par->keys[ind2], p->keys[ind]))
	    btree_set_parent_key(t, par, ind2, str, 1);
    }

    if (p->keys[ind])
	free(p->keys[ind]);
    p->keys[ind] = strdup(str);

    btree_node_put(t->cd, p);
    btree_dec_ref(t->cd, p);

    return 0;
}

static int btree_insert_key(btree_t *t, btree_node_t *n, int ind,
			    char *str, BTRec child) {
    int i, j;
    btree_node_t *par = NULL;

    if (n->parent)
	if (NULL == (par = btree_node_get(t->cd, n->parent)))
	    return -1;

    i = ind == -1 ? btree_find_key(n, str) : ind;
    //i = btree_find_key(n, str);
    
    /* Insert to this node. Does it have too many items yet? */
    for (j = n->used; j > i; j--) {
	n->keys[j] = n->keys[j-1];
	n->chld[j] = n->chld[j-1];
    }
    n->keys[i] = strdup(str);
    n->chld[i] = child;

    if (n->used == ind && par) {
	/* Inserted to end of list, so update parent too */
	ind = btree_find_key(par, n->keys[n->used-1]);
	btree_set_parent_key(t, par, ind, str, 1);
    }

    if (++n->used > BTREE_MAX) {
	btree_node_t *n2;

	/*
	 * At this point we split our tree in two.
	 * The copy the right hand half into a new node 'n2' and insert
	 * the highest key in the left hand int the parent (it should
	 * already have the highest from n2 if needed).
	 */
	n2 = btree_node_new(t->cd);
	n2->leaf = n->leaf;
	for (i = BTREE_MAX/2, j = 0; i <= BTREE_MAX; i++, j++) {
	    n2->keys[j] = n->keys[i]; n->keys[i] = NULL;
	    n2->chld[j] = n->chld[i]; n->chld[i] = 0;
	    if (!n2->leaf && n2->chld[j]) {
		btree_node_t *c = btree_node_get(t->cd, n2->chld[j]);
		c->parent = n2->rec;
		btree_node_put(t->cd, c);
	    }
	}
	n->used = BTREE_MAX/2;
	n2->used = j;
	n2->parent = n->parent;
	if (n->leaf) {
	    n2->next = n->next;
	    n->next = n2->rec;
	}

	btree_inc_ref(t->cd, n);
	btree_inc_ref(t->cd, n2);

	if (n->parent) {
	    ind = btree_find_key(par, n->keys[n->used-1]);
	    btree_set_parent_key(t, par, ind, n->keys[n->used-1], 0);
	    btree_insert_key(t, par, -1, n2->keys[n2->used-1], n2->rec);
	} else {
	    /* We need to make a new root */
	    btree_dec_ref(t->cd, t->root);
	    t->root = btree_node_new(t->cd);
	    btree_inc_ref(t->cd, t->root);
	    t->root->leaf = 0;
	    t->root->chld[0] = n->rec;
	    t->root->chld[1] = n2->rec;
	    t->root->keys[0] = strdup(n->keys[n->used-1]);
	    t->root->keys[1] = strdup(n2->keys[n2->used-1]);
	    t->root->used = 2;
	    n->parent = t->root->rec;
	    n2->parent = t->root->rec;

	    btree_node_put(t->cd, t->root);
	}
	
	btree_node_put(t->cd, n);
	btree_node_put(t->cd, n2);

	btree_dec_ref(t->cd, n);
	btree_dec_ref(t->cd, n2);
    } else {
	btree_node_put(t->cd, n);
    }

    return 0;
}

/*
 * Searches for 'str' in a btree and returns the node found.
 * 'index' : n->keys[index1] == str (str match, not ptr)
 */
static btree_node_t *btree_find_recurse(btree_t *t, char *str, int *index) {
    btree_node_t *n = t->root;
    int i;

    /* Find the leaf node to insert into */
    for (;;) {
	i = btree_find_key(n, str);

	if (n->leaf)
	    break;
	else
	    n = n->chld[i]
		? btree_node_get(t->cd, n->chld[i])
		: btree_node_get(t->cd, n->chld[i-1]);

	if (!n)
	    return NULL;
    }

    if (index)
	*index = i;

    return n;
}

int btree_insert(btree_t *t, char *str, BTRec value) {
    btree_node_t *n;
    int ind;

    n = btree_find_recurse(t, str, &ind);
    if (!n || !n->keys[ind] || 0 != strcmp(n->keys[ind], str))
	return btree_insert_key(t, n, ind, str, value);
    /* else already existed */

    //printf("Dup key %s\n", str);
    return btree_insert_key(t, n, ind, str, value);
}


BTRec btree_search(btree_t *t, char *str, int prefix) {
    int ind;
    btree_node_t *n = btree_find_recurse(t, str, &ind);

    if (prefix) {
	if (!n || !n->keys[ind] || 0 != strncmp(n->keys[ind], str,strlen(str)))
	    return -1;
    } else {
	if (!n || !n->keys[ind] || 0 != strcmp(n->keys[ind], str))
	    return -1;
    }

    return n->chld[ind];
}

BTRec *btree_search_all(btree_t *t, char *str, int prefix, int *nhits) {
    int ind;
    btree_node_t *n = btree_find_recurse(t, str, &ind);
    BTRec *results = NULL;
    size_t res_size = 0;
    size_t res_alloc = 0;

    if (!n || !n->keys[ind]) {
	*nhits = 0;
	return NULL;
    }

    /*
     * n->keys[ind] is the first hit, or the nearest one preceeding it.
     * We now iterate from here on.
     */
    do {
	if (( prefix && strncmp(n->keys[ind], str, strlen(str)) == 0) ||
	    (!prefix && strcmp(n->keys[ind], str) == 0)) {
	    /* Store a match */
	    if (res_size >= res_alloc) {
		res_alloc = res_alloc ? res_alloc*2 : 16;
		results = realloc(results, res_alloc * sizeof(*results));
		if (!results)
		    return NULL;
	    }
	    results[res_size++] = n->chld[ind];

	    /* Increment to next potential one */
	    if (++ind >= n->used) {
		n = n->next ? btree_node_get(t->cd, n->next) : NULL;
		ind = 0;
	    }
	} else {
	    break;
	}
    } while (n);

    *nhits = res_size;
    return results;
}

/*
 * Returns a new btree iterator, optionally starting with the first match
 * on or after prefix 'str'. If str is NULL we start from the first node
 * in the tree.
 *
 * Call btree_next() to fetch data and btree_iter_del() to destroy an
 * iterator. Updating a btree will make operations on an active iterator
 * undefined.
 *
 * Returns iterator on success
 *         NULL on failure
 */
btree_iter_t *btree_iter_new(btree_t *t, char *str) {
    btree_iter_t *iter = malloc(sizeof(*iter));
    if (!iter)
	return NULL;

    iter->ind = 0;
    iter->t = t;
    iter->n = btree_find_recurse(t, str ? str : "", &iter->ind);
    if (!iter->n || !iter->n->keys[iter->ind]) {
	free(iter);
	return NULL;
    }

    return iter;
}

void btree_iter_del(btree_iter_t *iter) {
    if (iter)
	free(iter);
}

/*
 * Iterates through a btree return key,value pairs.
 * Note: the key memory is only guaranteed to be valid until the next
 * (tg_)cache query.
 * 
 * Returns the key and sets *rec on success.
 *         NULL on failure (eot)
 */
char *btree_next(btree_iter_t *iter, BTRec *rec) {
    if (!iter || !iter->n)
	return NULL;

    while (iter->ind >= iter->n->used) {
	if (!iter->n->next)
	    return NULL;

	iter->n = btree_node_get(iter->t->cd, iter->n->next);
	iter->ind = 0;
    }

    if (rec) {
	*rec = iter->n->chld[iter->ind];
    }

    return iter->n->keys[iter->ind++];
}

/* Redistributes data between 'from' and 'to' evenly */
static int redist(btree_t *t, btree_node_t *to, btree_node_t *from,
		  int to_ind, int from_ind, int front) {
    int i, j, c;
    int nmove = (to->used + from->used)/2 - to->used;

    if (front) {
	/* Shuffle up to make room */
	for (i = to->used-1; i >= 0; i--) {
	    to->keys[i+nmove] = to->keys[i];
	    to->chld[i+nmove] = to->chld[i];
	}
	i = 0;
	j = from->used - nmove;
    } else {
	i = to->used;
	j = 0;
    }

    /* Move the keys over */
    for (c = 0; c < nmove; c++, i++, j++) {
	to->keys[i] = from->keys[j]; from->keys[j] = NULL;
	to->chld[i] = from->chld[j]; from->chld[j] = 0;
	if (!to->leaf && to->chld[i]) {
	    btree_node_t *n = btree_node_get(t->cd, to->chld[i]);
	    n->parent = to->rec;
	    btree_node_put(t->cd, n);
	}
    }

    to->used += nmove;
    from->used -= nmove;

    if (!front) {
	/* Shuffle down 'from' */
	for (i = 0, j = nmove; i < from->used; i++, j++) {
	    from->keys[i] = from->keys[j]; from->keys[j] = NULL;
	    from->chld[i] = from->chld[j]; from->chld[j] = 0;
	}
    }

    btree_node_put(t->cd, from);
    btree_node_put(t->cd, to);

    btree_set_parent_key(t, btree_node_get(t->cd, to->parent),
			 to_ind, to->keys[to->used-1], 1);
    btree_set_parent_key(t, btree_node_get(t->cd, from->parent),
			 from_ind, from->keys[from->used-1], 1);

    return 0;
}

static int btree_delete_key(btree_t *t, btree_node_t *n, int ind, char *str) {
    int i, j;
    btree_node_t *left, *right;

    do {
	btree_node_t *par;

	btree_inc_ref(t->cd, n);

	/* Remove index 'ind' */
	if (n->keys[ind])
	    free(n->keys[ind]);

	for (i = ind+1; i < n->used; i++) {
	    n->keys[i-1] = n->keys[i];
	    n->chld[i-1] = n->chld[i];
	}
	n->keys[i-1] = NULL;
	n->chld[i-1] = 0;
	n->used--;

	/* If root and one node left, collapse root */
	if (t->root == n && n->used == 1 && !n->leaf && n->chld[0]) {
	    btree_dec_ref(t->cd, t->root);
	    btree_node_del(t->cd, t->root);
	    t->root = btree_node_get(t->cd, n->chld[0]);
	    btree_inc_ref(t->cd, t->root);
	    t->root->parent = 0;
	    btree_node_put(t->cd, t->root);
	    btree_dec_ref(t->cd, n);
	    return 0;
	}

	/* Replace 'ind' with our parent's index to us */
	par = n->parent ? btree_node_get(t->cd, n->parent) : NULL;
	if (par) {
	    btree_inc_ref(t->cd, par);
	    ind = btree_find_key(par, str);

	    if (par->keys[ind] && 0 == strcmp(par->keys[ind], str))
		btree_set_parent_key(t, par, ind, n->keys[n->used-1], 1);
	}
	 
	if (n->used >= BTREE_MIN || !par) {
	    btree_node_put(t->cd, n);
	    btree_dec_ref(t->cd, n);
	    return 0;
	}

	/*
	 * We've underflowed so we have MAX/2-1 entries.
	 * We know our neighbour has >= MAX/2 entries.
	 *
	 * Either this + neighbour >= MAX items (=> redistribute)
	 * or     this + neighbour <  MAX items (=> merge and recurse up to
	 *                                          parent).
	 */
	left  = (ind-1 >= 0)
	    ? btree_node_get(t->cd, par->chld[ind-1])
	    : NULL;
	if (left) btree_inc_ref(t->cd, left);

	right = (ind+1 <= par->used-1)
	    ? btree_node_get(t->cd, par->chld[ind+1])
	    : NULL;
	if (right) btree_inc_ref(t->cd, right);


	/* Can we redistribute? */
	if (left && left->used  > BTREE_MIN) {
	    redist(t, n, left, ind, ind-1, 1);
	    btree_dec_ref(t->cd, n);
	    btree_dec_ref(t->cd, left);
	    if (right) btree_dec_ref(t->cd, right);
	    break;

	} else if (right && right->used > BTREE_MIN) {
	    redist(t, n, right, ind, ind+1, 0);
	    btree_dec_ref(t->cd, n);
	    if (left) btree_dec_ref(t->cd, left);
	    btree_dec_ref(t->cd, right);
	    break;

	} else {
	    /* Merge with one of our neighbours */
	    if (left) {
		/* left = left; */
		if (right) btree_dec_ref(t->cd, right);
		right = n;
		btree_inc_ref(t->cd, right);
		ind--;
	    } else {
		if (left) btree_dec_ref(t->cd, left);
		left = n;
		btree_inc_ref(t->cd, left);
		/* right = right; */
	    }

	    if (left == right) {
		/* we're already a one-node item, recurse up again */
		btree_dec_ref(t->cd, n);
		n = btree_node_get(t->cd, left->parent);
	    } else {
		/* Merge right into left */
		for (j = left->used, i = 0; i < right->used; i++, j++) {
		    left->keys[j] = right->keys[i]; right->keys[i] = NULL;
		    left->chld[j] = right->chld[i]; right->chld[i] = 0;
		    if (!left->leaf && left->chld[j]) {
			btree_node_t *c = btree_node_get(t->cd, left->chld[j]);
			c->parent = left->rec;
			btree_node_put(t->cd, c);
		    }
		}
		left->used = j;
		right->used = 0;
		left->next = right->next;
		btree_node_put(t->cd, left);
		
		/* Recurse */
		btree_dec_ref(t->cd, n);
		n = btree_node_get(t->cd, left->parent);
		btree_inc_ref(t->cd, n);
		btree_set_parent_key(t, n, ind, left->keys[left->used-1], 1);
		ind++; /* right node */
		str = left->keys[left->used-1];

		btree_dec_ref(t->cd, right);
		btree_node_del(t->cd, right);
		right = NULL;
		btree_dec_ref(t->cd, n);
	    }

	    if (left)  btree_dec_ref(t->cd, left);
	    if (right) btree_dec_ref(t->cd, right);
	}

    } while (n);

    return 0;
}

/*
 * Deletes the first item it finds with string 'str'
 */
int btree_delete(btree_t *t, char *str) {
    btree_node_t *n;
    int ind;

    n = btree_find_recurse(t, str, &ind);
    if (!n || !n->keys[ind] || 0 != strcmp(n->keys[ind], str))
	return 0; /* No match */

    return btree_delete_key(t, n, ind, str);
}

/*
 * Deletes a specific pairing of (str,rec).
 */
int btree_delete_rec(btree_t *t, char *str, BTRec rec) {
    btree_node_t *n;
    int ind;

    n = btree_find_recurse(t, str, &ind);
    if (!n || !n->keys[ind] || 0 != strcmp(n->keys[ind], str))
	return 0; /* No match */

    /* n->keys[ind] is a hit, so search from here on finding rec */
    while (strcmp(n->keys[ind], str) == 0) {
	/* Found it? If so remove */
	if (n->chld[ind] == rec)
	    return btree_delete_key(t, n, ind, str);

	/* Otherwise increment to next potential match */
	if (++ind >= n->used) {
	    n = n->next ? btree_node_get(t->cd, n->next) : NULL;
	    ind = 0;
	}
    }

    return 0;
}

#define INDENT (depth ? printf("%*c", depth, ' ') : 0)
void btree_print(btree_t *t, btree_node_t *n, int depth) {
    int i;

    if (depth == 0)
	puts("");
    INDENT; printf("Node %"PRIbtr", leaf=%d, parent %"PRIbtr
		   ", next %"PRIbtr", used %d\n",
		   n->rec, n->leaf, n->parent, n->next, n->used);
    btree_inc_ref(t->cd, n);
    for (i = 0; i < n->used; i++) {
	INDENT; printf("key %d = %s val %"PRIbtr"\n",
		       i, n->keys[i] ? n->keys[i] : "-", n->chld[i]);
	if (!n->leaf && n->chld[i]) {
	    btree_print(t, btree_node_get(t->cd, n->chld[i]), depth+2);
	}
    }
    btree_dec_ref(t->cd, n);
}

char *btree_check(btree_t *t, btree_node_t *n, char *pleaf) {
    int i;
    char *str = NULL;
    char *prev = pleaf;

    btree_inc_ref(t->cd, n);
    for (i = 0; i < n->used; i++) {
	assert(n->keys[i]);
	assert(strcmp(n->keys[i], prev) > 0);
	prev = n->keys[i];
	if (n->leaf) {
	    //assert(n->chld[i] == 0);
	    pleaf = str = n->keys[i];

	    if (n->next && i == n->used-1)
		assert(strcmp(n->keys[i],
			      btree_node_get(t->cd, n->next)->keys[0]) < 0);
	} else {
	    btree_node_t *c = btree_node_get(t->cd, n->chld[i]);
	    assert(c);
	    assert(c->parent == n->rec);
	    str = btree_check(t, c, pleaf);
	    assert(strcmp(n->keys[i], str) == 0);
	}
    }
    btree_dec_ref(t->cd, n);

    return str;
}

#define FRONT_COMPRESSION
int btree_size(btree_t *t, btree_node_t *n) {
    int sz = 0, sz2 = 0;
    int i;
    unsigned char c;
    char *last = NULL;
    
    sz++;  /* nused */
    sz++; /* leaf */
    c = n->used; write(1, &c, 1);
    c = n->leaf; write(1, &c, 1);
    for (i = 0; i < n->used; i++) {
#ifdef FRONT_COMPRESSION
	if (!last) {
	    sz += 1+strlen(n->keys[i])+1;
	    c = 0; write(1, &c, 1);
	    write(1, n->keys[i], strlen(n->keys[i])+1);
	} else {
	    char *cp1 = n->keys[i], *cp2 = last;
	    while (*cp1 == *cp2)
		cp1++, cp2++;
	    c = cp2-last;
	    write(1, &c, 1), sz++;

	    //cp1[strlen(cp1)-1] |= 0x80;
	    write(1, cp1, strlen(cp1)+1);
	    sz += strlen(cp1)+1;
	}
	last = n->keys[i];
#else
	sz += strlen(n->keys[i])+1;
	write(1, n->keys[i], strlen(n->keys[i])+1);
#endif
	write(1, &n->chld[i], 4);
	sz += 4; /* rec.no */
	if (!n->leaf && n->chld[i]) {
	    sz2 += btree_size(t, btree_node_get(t->cd, n->chld[i]));
	}
    }

    //fprintf(stderr, "Size of this tree = %d\n", sz);

    return sz + sz2;
}

void btree_list(btree_t *t, char *prefix) {
    btree_node_t *n;
    int ind, i;
    size_t l = strlen(prefix);

    n = btree_find_recurse(t, prefix, &ind);
    for (i = ind; n;) {
	if (0 != strncmp(prefix, n->keys[i], l))
	    break;

	printf("%s\n", n->keys[i]);
	if (++i == n->used) {
	    n = btree_node_get(t->cd, n->next);
	    i = 0;
	}
    }
}

/*
 * A btree node is encoded on disk as follows:
 * 1x 1byte: Flags (bit 0 => leaf flag, bits 1-7 reserverd)
 * 1x 1byte: Number of used keys/children (N).
 * 1x 4byte: Parent node record (0 if none)
 * 1x 4byte: Next node record (0 if none, eg internal)
 * Nx 4byte: N children record nums.
 *           Either other node records (leaf==0) or the data itself (leaf==1)
 * Nx {
 *    1byte: Shared key prefix length (common data from previous key)
 *    ?byte: Key suffix
 *    1byte: '\0' - key nul termination char
 * }
 */


/*
 * Decodes the on-disk btree format into an in-memory C struct.
 *
 * Returns allocated btree_node_t on success
 *         NULL on failure
 */
btree_node_t *btree_node_decode(unsigned char *buf) {
    btree_node_t *n;
    unsigned char *bufp;
    int i;
    char *last;

    if (NULL == (n = btree_new_node()))
	return NULL;

    /* Static data */
    n->leaf = buf[0];
    n->used = buf[1];

    n->parent =
	(buf[2] << 24) |
	(buf[3] << 16) |
	(buf[4] <<  8) |
	(buf[5] <<  0);

    n->next =
	(buf[6] << 24) |
	(buf[7] << 16) |
	(buf[8] <<  8) |
	(buf[9] <<  0);

    /* Array of children */
    for (i = 0; i < n->used; i++) {
	n->chld[i] =
	    (buf[10+i*4] << 24) |
	    (buf[11+i*4] << 16) |
	    (buf[12+i*4] <<  8) |
	    (buf[13+i*4] <<  0);
    }

    /* And finally the variable sizes elements - the keys */
    last = "";
    bufp = &buf[10 + n->used*4];
    for (i = 0; i < n->used; i++) {
	int dist = *bufp++;
	size_t l = strlen((char *)bufp);
	n->keys[i] = (char *)malloc(dist + l + 1);
	if (dist)
	    strncpy(n->keys[i], last, dist);
	strcpy(n->keys[i]+dist, (char *)bufp);
	bufp += l+1;

	last = n->keys[i];
    }

    return n;
}

/*
 * Converts an in-memory btree_node_t struct to a serialised character stream
 * suitable for storing on disk.
 *
 * Returns malloced char* on success, also setting *size to the length used.
 *         NULL on failure
 */
unsigned char *btree_node_encode(btree_node_t *n, size_t *size) {
    size_t alloc;
    unsigned char *buf, *bufp;
    int i;
    char *last;

    /* Build static component */
    alloc = 1+1+4+4 + 4*n->used + 8*n->used /* guesstimate */;
    if (NULL == (buf = (unsigned char *)malloc(alloc)))
	return NULL;

    assert(n->used <= 255);
    buf[0] = n->leaf;
    buf[1] = n->used;

    buf[2] = ((int)n->parent >> 24) & 0xff;
    buf[3] = ((int)n->parent >> 16) & 0xff;
    buf[4] = ((int)n->parent >>  8) & 0xff;
    buf[5] = ((int)n->parent >>  0) & 0xff;

    buf[6] = ((int)n->next >> 24) & 0xff;
    buf[7] = ((int)n->next >> 16) & 0xff;
    buf[8] = ((int)n->next >>  8) & 0xff;
    buf[9] = ((int)n->next >>  0) & 0xff;

    for (i = 0; i < n->used; i++) {
	buf[10+i*4] = ((int)n->chld[i] >> 24) & 0xff;
	buf[11+i*4] = ((int)n->chld[i] >> 16) & 0xff;
	buf[12+i*4] = ((int)n->chld[i] >>  8) & 0xff;
	buf[13+i*4] = ((int)n->chld[i] >>  0) & 0xff;
    }

    /* And add in the variable sized element - key strings */
    bufp = &buf[10 + n->used*4];
    last = "";
    for (i = 0; i < n->used; i++) {
	char *cp1, *cp2;

	/* Determine the common prefix length with last key */
	cp1 = n->keys[i];
	cp2 = last;

	while (*cp1 == *cp2 && *cp2)
	    cp1++, cp2++;

	/* Realloc buf as needed, ensuring bufp updates too */
	while (bufp + strlen(cp1) + 2 - buf >= alloc) {
	    size_t diff = bufp - buf;
	    alloc += 1000;
	    buf = (unsigned char *)realloc(buf, alloc);
	    bufp = buf + diff;
	}

	/* Copy it over */
	*bufp++ = cp2-last;
	while((*bufp++ = *cp1++))
	    ;

	last = n->keys[i];
    }

    *size = bufp - buf;

    return buf;
}

/*
 * Decodes the on-disk btree format into an in-memory C struct.
 *
 * Returns allocated btree_node_t on success
 *         NULL on failure
 */
btree_node_t *btree_node_decode2(unsigned char *buf, int fmt) {
    btree_node_t *n;
    unsigned char *bufp, *bufp2, *bufp3;
    int i;
    char *last;

    if (NULL == (n = btree_new_node()))
	return NULL;

    /* Static data */
    n->leaf = buf[0];
    n->used = (buf[1] << 8) | (buf[2] << 0);

    if (fmt == 1) {
	n->parent =
	    (buf[4] << 24) |
	    (buf[5] << 16) |
	    (buf[6] <<  8) |
	    (buf[7] <<  0);

	n->next =
	    (buf[8]  << 24) |
	    (buf[9]  << 16) |
	    (buf[10] <<  8) |
	    (buf[11] <<  0);

	bufp = &buf[12];
	for (i = 0; i < n->used; i++) {
	    uint32_t i4;
	    bufp += u72int(bufp, &i4); n->chld[i] = i4;
	}
    } else {
	uint64_t i8;

	bufp = &buf[4];
	bufp += u72intw(bufp, &i8); n->parent = i8;
	bufp += u72intw(bufp, &i8); n->next = i8;

	for (i = 0; i < n->used; i++) {
	    bufp += u72intw(bufp, &i8); n->chld[i] = i8;
	}
    }

    /* And finally the variable sizes elements - the keys */
    last = "";
    bufp2 = bufp;
    bufp3 = bufp2 + n->used;
    bufp  = bufp3 + n->used;

    for (i = 0; i < n->used; i++) {
	int dist = *bufp2++;
	size_t l = *bufp3++;
	n->keys[i] = (char *)malloc(dist + l + 1);
	if (dist)
	    strncpy(n->keys[i], last, dist);
	strncpy(n->keys[i]+dist, (char *)bufp, l);
	*(n->keys[i]+dist+l) = 0;
	bufp += l;

	last = n->keys[i];
    }
    for (;i < BTREE_MAX; i++) {
	n->keys[i] = NULL;
	n->chld[i] = 0;
    }

    return n;
}

/*
 * Converts an in-memory btree_node_t struct to a serialised character stream
 * suitable for storing on disk.
 *
 * The encoded stream comes in several distinct parts with size
 * returned in parts. These consist of:
 *
 * 1) A small 12-byte header + the child record numbers
 * 2) The size of the common prefix between each key string.
 * 3) The size of the key suffixes (bit that differs to previous key)
 * 4) the key suffixes themselves.
 *
 * Returns malloced char* on success, also setting size & parts[0..3]
 *         NULL on failure
 */
unsigned char *btree_node_encode2(btree_node_t *n, size_t *size,
				  size_t *parts, int fmt) {
    size_t alloc;
    unsigned char *buf, *bufp, *bufp2, *bufp3;
    int i;
    char *last;

    /* Build static component */
    alloc = 1+2+1+4+4 + 4*n->used + 8*n->used /* guesstimate */;
    if (NULL == (buf = (unsigned char *)malloc(alloc)))
	return NULL;

    /*
    printf("Encode node %d used %d leaf %d\n",
	   n->rec, n->used, n->leaf);
    for (i = 0; i < n->used; i++) {
	printf("   %4d %10d %s\n", i, n->chld[i],
	       n->keys[i] ? n->keys[i] : "");
    }
    */

    assert(n->used <= 65535);
    buf[0] = n->leaf;
    buf[1] = (n->used >> 8) & 0xff;
    buf[2] = (n->used >> 0) & 0xff;
    buf[3] = 0;

    if (fmt == 1) {
	buf[4] = ((int)n->parent >> 24) & 0xff;
	buf[5] = ((int)n->parent >> 16) & 0xff;
	buf[6] = ((int)n->parent >>  8) & 0xff;
	buf[7] = ((int)n->parent >>  0) & 0xff;

	buf[8] = ((int)n->next >> 24) & 0xff;
	buf[9] = ((int)n->next >> 16) & 0xff;
	buf[10]= ((int)n->next >>  8) & 0xff;
	buf[11]= ((int)n->next >>  0) & 0xff;

	bufp = &buf[12];
	for (i = 0; i < n->used; i++) {
	    bufp += int2u7((int)n->chld[i], bufp);
	}
    } else {
	bufp = &buf[4];
	bufp += intw2u7(n->parent, bufp);
	bufp += intw2u7(n->next, bufp);

	for (i = 0; i < n->used; i++) {
	    bufp += intw2u7(n->chld[i], bufp);
	}
    }


    /* And add in the variable sized element - key strings */

    if (parts) {
	parts[0] = bufp - buf;
	parts[1] = n->used;
	parts[2] = n->used;
    }

    last = "";
    bufp2 = bufp;
    bufp3 = bufp2 + n->used;
    bufp  = bufp3 + n->used;

    for (i = 0; i < n->used; i++) {
	char *cp1, *cp2;
	int len;

	/* Determine the common prefix length with last key */
	cp1 = n->keys[i];
	cp2 = last;

	while (*cp1 == *cp2 && *cp2)
	    cp1++, cp2++;

	/* Realloc buf as needed, ensuring bufp* updates too */
	while (bufp + strlen(cp1) + 2 - buf >= alloc) {
	    size_t diff  = bufp  - buf;
	    size_t diff2 = bufp2 - buf;
	    size_t diff3 = bufp3 - buf;
	    alloc += 1000;
	    buf = (unsigned char *)realloc(buf, alloc);
	    bufp  = buf + diff;
	    bufp2 = buf + diff2;
	    bufp3 = buf + diff3;
	}

	/* Copy it over */
	*bufp2++ = cp2-last;      /* prefix length */
	len = 0;
	while((*bufp++ = *cp1++)) /* suffix itself */
	    len++;
	bufp--;
	*bufp3++ = len;           /* suffix length */

	last = n->keys[i];
    }

    *size = bufp - buf;
    if (parts)
	parts[3] = bufp - buf - parts[0] - parts[1] - parts[2];

    return buf;
}

#ifdef TEST_MAIN
/* -------------------------------------------------------------------------
 * Test code. This needs to supply the basic btree_node_* functions
 * (get, put, new and del) allow with a main() test function.
 */

/*
 * Mapping of record numbers to pointers. This is used in simulating
 * btree_node_get and btree_node_put method calls.
 */
#define MAX_REC 400000
btree_node_t *rec_map[MAX_REC];

btree_node_t *btree_node_get(void *cd, BTRec r) {
    return (r >= 0 && r < MAX_REC) ? rec_map[r] : NULL;
}

int btree_node_put(void *cd, btree_node_t *n) {
    return 0;
}

btree_node_t *btree_node_new(void *cd) {
    static int i = 1;
    int j;

    for (j = 0; j < MAX_REC; j++) {
	if (!rec_map[i]) {
	    rec_map[i] = btree_new_node();
	    rec_map[i]->rec = i;
	    return rec_map[i];
	}

	if (++i == MAX_REC)
	    i = 1;
    }

    fprintf(stderr, "Ran out of nodes\n");

    return NULL;
}

void btree_node_del(void *cd, btree_node_t *n) {
    rec_map[n->rec] = NULL;
    btree_del_node(n);
}

void btree_inc_ref(void *cd, btree_node_t *n) {}
void btree_dec_ref(void *cd, btree_node_t *n) {}

char *lines[1000000];
int main(int argc, char **argv) {
    char line[1024];
    btree_t *tree = btree_new(NULL, 0);
    FILE *fp = argc == 1 ? stdin : fopen(argv[1], "r");
    int nlines = 0;
    int i;

    while (fgets(line, 1024, fp)) {
	char *cp;
	if ((cp = strchr(line, '\n')))
	    *cp = 0;

	//puts(line);
	btree_insert(tree, line, nlines);
	//btree_print(tree, tree->root, 0);
	//btree_check(tree, tree->root, "");
	
	lines[nlines++] = strdup(line);
	if (nlines == 1000000)
	    break;
    }

    btree_check(tree, tree->root, "");
    //btree_print(tree, tree->root, 0);
    //return 0;

    if (argc >= 3)
	srandom(atoi(argv[2]));
#if 1
    for (i = 0; i < 1000000; i++) {
	int n = random() % nlines;
	if (random()&1) {
	    //fprintf(stdout, "Insert %s %d\n", lines[n], n);
	    btree_insert(tree, lines[n], n);
	} else {
	    //fprintf(stdout, "Delete %s\n", lines[n]);
	    btree_delete(tree, lines[n]);
	}

	//btree_print(tree, tree->root, 0);

	if ((i % 10000) == 0) {
	    putchar('.');
	    fflush(stdout);
	    btree_check(tree, tree->root, "");
	}
	//btree_check(tree, tree->root, "");
    }
    puts("");
#endif

    btree_check(tree, tree->root, "");
    btree_print(tree, tree->root, 2);
    btree_list(tree, "yt");

    //fprintf(stderr, "btree size = %d\n", btree_size(tree->root));

    btree_del(tree);

    return 0;
}

#endif
