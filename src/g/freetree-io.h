/*
 * A very simple implementation of C++ style templates.
 *
 * Before #including this file define the TYPE macro.
 */

/*
 * This contains 32-bit and 64-bit functions to load and save the
 * freetree data.
 * The freetree is saved to the end of the AUX file and has the following
 * format:
 *
 * HEADER: (1 of)
 *    No. blocks
 *    Last time
 *    Free nodes number
 *    Tree root node number
 *    Rover node number
 *
 * NODES: (No.blocks * BLOCK_FACTOR copies)
 *    Left node number
 *    Right node number
 *    Parent node number
 *    Position
 *    Length
 *
 * CRC32 (1 of)
 *    32-bit CRC of all the above data
 */

#define PASTE(a,b) a##b
#define PASTY(a,b) PASTE(a,b)

int PASTY(freetree_save_,TYPE)(free_tree *t, int fd, int last_time) {
    int i, j;
    free_tree_n *n;
    TYPE *data;
    int data_ind = 0;
    int endian = 1;

    /* Allocate a block of data to write our tree into */
    data = (TYPE *)xmalloc((t->nblocks * BLOCK_FACTOR * 5 + 6) * sizeof(TYPE));
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
	    /*
	     * n->len may be MAXINT8, but strip this down to the max allowed
	     * within this data type
	     */
	    data[data_ind++] = n->len & (((uint8)1 << (sizeof(TYPE)*8-1))-1);
	}
    }

    if (*(char *)&endian) {
	int i;
	for (i = 0; i < data_ind; i++)
	    PASTY(swap_,TYPE)(data[i], data[i]);
    }

    data[data_ind] = calc_crc((char *)data, data_ind * sizeof(TYPE));
    data_ind++;

    if (*(char *)&endian) {
	PASTY(swap_,TYPE)(data[data_ind-1], data[data_ind-1]);
    }

    /* LOW LEVEL IO HERE */
    write(fd, data, data_ind * sizeof(TYPE));

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
free_tree *PASTY(freetree_load_,TYPE)(int fd, int last_time) {
    free_tree *t = NULL;
    TYPE *data = NULL;
    int4 nblocks;
    int endian = 1;
    int sz;
    int data_ind = 0;
    int last_block;
    int i, j;

    /* Load the first 5 blocks */
    sz = 5 * sizeof(TYPE);
    if (!(data = (TYPE *)xmalloc(sz)))
	goto error;

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (sz != read(fd, data, sz))
	goto error;

    /* Check time stamp */
    if (PASTY(be_,TYPE)(data[1]) != last_time) {
	fprintf(stderr, "Incorrect tree timestamp - rebuilding tree\n");
	goto error;
    }

    /* Load the other blocks */
    nblocks = PASTY(be_,TYPE)(data[0]);
    last_block = nblocks * BLOCK_FACTOR * 5 + 5;
    sz =  (last_block+1) * sizeof(TYPE);
    if (!(data = (TYPE *)xrealloc(data, sz)))
	goto error;

    sz -= 5*sizeof(TYPE); /* first 5 already loaded */

    /* LOW LEVEL IO HERE */
    errno = 0;
    if (sz != read(fd, &data[5], sz)) {
	fprintf(stderr, "Tree too short\n");
	goto error;
    }

    /* Check CRC */
    if (*(char *)&endian) {
	PASTY(swap_,TYPE)(data[last_block], data[last_block]);
    }
    if (data[last_block] != calc_crc((char *)data, last_block*sizeof(TYPE))) {
	fprintf(stderr, "Invalid tree CRC\n");
	goto error;
    }

    /* Byte-swap blocks */
    if (*(char *)&endian) {
	int i;

	for (i = 0; i < last_block ; i++)
	    PASTY(swap_,TYPE)(data[i], data[i]);
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
