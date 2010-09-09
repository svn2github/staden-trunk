#ifndef _BTREE2_H_
#define _BTREE2_H_

/* The order of the tree. Keep even for now */
//#define BTREE_MAX 100
#define BTREE_MAX 4000
#define BTREE_MIN (BTREE_MAX/2)

typedef int64_t BTRec;
#define PRIbtr PRId64

typedef struct btree_node {
    char *keys[BTREE_MAX+1];
    BTRec rec; /* this record number */
    BTRec chld[BTREE_MAX+1]; /* zero => none */
    BTRec parent;
    BTRec next;
    int leaf;
    int used;
    void *cache; /* Used for linking back to the GView info */
} btree_node_t;

typedef struct {
    void *cd; /* arbitrary client data to pass back to user funcs. */
    btree_node_t *root;
} btree_t;

btree_t *btree_new(void *cd, BTRec root);
void btree_del(btree_t *t);
btree_node_t *btree_new_node(void);
void btree_del_node(btree_node_t *n);

/*
 * Converts an in-memory btree_node_t struct to a serialised character stream
 * suitable for storing on disk.
 *
 * Returns malloced char* on success, also setting *size to the length used.
 *         NULL on failure
 */
unsigned char *btree_node_encode(btree_node_t *n, size_t *size);
unsigned char *btree_node_encode2(btree_node_t *n, size_t *size, size_t *prts,
				  int fmt);

/*
 * Decodes the on-disk btree format into an in-memory C struct.
 *
 * Returns allocated btree_node_t on success
 *         NULL on failure
 */
btree_node_t *btree_node_decode(unsigned char *buf);
btree_node_t *btree_node_decode2(unsigned char *buf, int fmt);

int btree_insert(btree_t *t, char *str, BTRec value);
int btree_delete(btree_t *t, char *str);
BTRec btree_search(btree_t *t, char *str, int prefix);
BTRec *btree_search_all(btree_t *t, char *str, int prefix, int *nhits);

/* Access methods to be supplied by whatever is calling the btree code */
extern btree_node_t *btree_node_get(void *cd, BTRec r);
extern int btree_node_put(void *cd, btree_node_t *n);
extern btree_node_t *btree_node_new(void *cd);
extern void btree_node_del(void *cd, btree_node_t *n);
extern void btree_inc_ref(void *cd, btree_node_t *n);
extern void btree_dec_ref(void *cd, btree_node_t *n);

#endif /* _BTREE2_H_ */
