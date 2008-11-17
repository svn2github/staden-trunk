#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <fcntl.h>

#include "g.h"
#include "misc.h"
#include "tg_iface_g.h"
#include "tg_utils.h"
#include "b+tree2.h"
/* #include "io_lib/deflate_interlaced.h" */

/* #define INDEX_NAMES */

static iface iface_g;

/* ------------------------------------------------------------------------
 * This file houses the Gap5 interface to the "g library".
 * No knowledge of the g library should be visible outside of this one file.
 *
 * See tg_face.h for more information on what the database interface layer
 * needs to provide.
 *
 * None of this code should be externally visible except for the one single
 * exported function listed in tg_iface_g.h.
 */


#define DATABASE_RECORD 0

/*
 * Internal g-library connection information.
 * This struct is never visible outside of this file.
 */
typedef struct {
    GDB *gdb;
    GClient client;
    GLock mode;
    HacheTable *seq_name_hash;
    btree_t *seq_name_tree;
    HacheTable *contig_name_hash;
    btree_t *contig_name_tree;
} g_io;


#if 0
/* ------------------------------------------------------------------------ */
/*
 * Data compression routines. These hook into the interlaced-deflate
 * algorithm using in SRF/ZTR.
 */

#define SYM_EOF 256
static huffman_code_t codes_slx_qual[] = {
    { 40,  1}, { 10,  5}, {  3,  6}, {  4,  6}, {  5,  6}, { 11,  6},
    { 12,  6}, { 13,  6}, { 14,  6}, { 15,  6}, { 16,  6}, { 17,  6},
    { 18,  6}, { 19,  6}, { 20,  6}, { 21,  6}, { 22,  6}, { 23,  6},
    { 24,  6}, { 25,  6}, { 26,  6}, { 27,  6}, { 29,  6}, { 30,  6},
    { 31,  6}, {  6,  7}, {  7,  7}, {  8,  7}, {  9,  7}, { 28,  7},
    { 32,  7}, { 33,  7}, { 34,  7}, { 35,  7}, { 36,  7}, { 37,  7},
    { 38,  7}, { 39,  7}, {  2, 10}, {229, 13}, {254, 13}, {  0, 15},
    {  1, 15}, { 41, 15}, { 42, 15}, { 43, 15}, { 44, 15}, { 45, 15},
    { 46, 15}, { 47, 15}, {'0', 15}, {'1', 15}, {'2', 15}, {'3', 15},
    {'4', 15}, {'5', 15}, {'6', 15}, {'7', 15}, {'8', 15}, {'9', 15},
    { 58, 15}, { 59, 15}, { 60, 15}, { 61, 15}, { 62, 15}, { 63, 15},
    { 64, 15}, {'A', 15}, {'B', 15}, {'C', 15}, {'D', 15}, {'E', 15},
    {'F', 15}, {'G', 15}, {'H', 15}, {'I', 15}, {'J', 15}, {'K', 15},
    {'L', 15}, {'M', 15}, {'N', 15}, {'O', 15}, {'P', 15}, {'Q', 15},
    {'R', 15}, {'S', 15}, {'T', 15}, {'U', 15}, {'V', 15}, {'W', 15},
    {'X', 15}, {'Y', 15}, {'Z', 15}, { 91, 15}, { 92, 15}, { 93, 15},
    { 94, 15}, { 95, 15}, { 96, 15}, {'a', 15}, {'b', 15}, {'c', 15},
    {'d', 15}, {'e', 15}, {'f', 15}, {'g', 15}, {'h', 15}, {'i', 15},
    {'j', 15}, {'k', 15}, {'l', 15}, {'m', 15}, {'n', 15}, {'o', 15},
    {'p', 15}, {'q', 15}, {'r', 15}, {'s', 15}, {'t', 15}, {'u', 15},
    {'v', 15}, {'w', 15}, {'x', 15}, {'y', 15}, {'z', 15}, {123, 15},
    {124, 15}, {125, 15}, {126, 15}, {127, 15}, {128, 15}, {129, 15},
    {130, 15}, {131, 15}, {132, 15}, {133, 15}, {134, 15}, {135, 15},
    {136, 15}, {137, 15}, {138, 15}, {139, 15}, {140, 15}, {141, 15},
    {142, 15}, {143, 15}, {144, 15}, {145, 15}, {146, 15}, {147, 15},
    {148, 15}, {149, 15}, {150, 15}, {151, 15}, {152, 15}, {153, 15},
    {154, 15}, {155, 15}, {156, 15}, {157, 15}, {158, 15}, {159, 15},
    {160, 15}, {161, 15}, {162, 15}, {163, 15}, {164, 15}, {165, 15},
    {166, 15}, {167, 15}, {168, 15}, {169, 15}, {170, 15}, {171, 15},
    {172, 15}, {173, 15}, {174, 15}, {175, 15}, {176, 15}, {177, 15},
    {178, 15}, {179, 15}, {180, 15}, {181, 15}, {182, 15}, {183, 15},
    {184, 15}, {185, 15}, {186, 15}, {187, 15}, {188, 15}, {189, 15},
    {190, 15}, {191, 15}, {192, 15}, {193, 15}, {194, 15}, {195, 15},
    {196, 15}, {197, 15}, {198, 15}, {199, 15}, {200, 15}, {201, 15},
    {202, 15}, {203, 15}, {204, 15}, {205, 15}, {206, 15}, {207, 15},
    {208, 15}, {209, 15}, {210, 15}, {211, 15}, {212, 15}, {213, 15},
    {214, 15}, {215, 15}, {216, 15}, {217, 15}, {218, 15}, {219, 15},
    {220, 15}, {221, 15}, {222, 15}, {223, 15}, {224, 15}, {225, 15},
    {226, 15}, {227, 15}, {228, 15}, {230, 15}, {231, 15}, {232, 15},
    {233, 15}, {234, 15}, {235, 15}, {236, 15}, {237, 15}, {238, 15},
    {239, 15}, {240, 15}, {241, 15}, {242, 15}, {243, 15}, {244, 15},
    {245, 15}, {246, 15}, {247, 15}, {248, 15}, {249, 15}, {250, 15},
    {251, 15}, {252, 15}, {253, 15}, {255, 15}, {SYM_EOF, 15},
};

/*
 * Returns block_t pointer on success
 *         NULL on failure.
 */
block_t *ideflate_encode(unsigned char *data, int len,
			 int code_set, int rec_size) {
    /* Encoding */
    int blk_size = 1<<20;
    unsigned char *d2 = data;
    block_t *blk;
    huffman_codeset_t *cs;
    static int init_code_sets = 0;
    huffman_codeset_t *cuser0;

    if (!init_code_sets) {
	cuser0 = codes2codeset(codes_slx_qual, 
			       sizeof(codes_slx_qual)/sizeof(huffman_code_t),
			       CODE_USER);
	init_code_sets = 1;
    }

    blk = block_create(NULL, 8192);

    do {
	int l2 = len > blk_size ? blk_size : len;

	if (code_set != 0)
	    l2 = len; /* predefined code-sets have final-block bit set */
	cs = generate_code_set(code_set, rec_size,
			       d2, l2, /* Data and length */
			       1,      /* eof */
			       MAX_CODE_LEN,
			       0);     /* all codes */
	if (!cs)
	    return NULL;

	store_codes(blk, cs, l2 == len);
	if (code_set != 0) {
	    blk->data[blk->byte = 0] = 0;  /* starting bit no. preseved */
	} else {
	    /*
	      fprintf(stderr, "codes finished at %d bytes, %d bits\n",
		    blk->byte, blk->bit);
	    */
	}

	huffman_multi_encode(blk, cs, code_set, d2, l2);

	huffman_codeset_destroy(cs);
	len -= l2;
	d2  += l2;

    } while (len > 0);

    return blk;
}

/*
 * Returns block_t ptr on success and sets new_data to point to the byte
 *         immediately following the compressed input.
 *         NULL for failure.
 */
block_t *ideflate_decode(unsigned char *data, int len, int code_set,
			 unsigned char **new_data) {
    huffman_codeset_t *cs = NULL;
    block_t *blk_in;
    int bfinal;
    block_t *out;

    blk_in = block_create(NULL, 1000 + len);

    /* Inefficient; see ztr's compression.c for caching this */
    if (code_set != 0) {
	cs = generate_code_set(code_set, 1,  /* no. codes */
			       NULL, 0,      /* data + size */
			       1,	     /* eof */
			       MAX_CODE_LEN,
			       0);	     /* all_codes */
	store_codes(blk_in, cs, 1);
    }

    if (blk_in->bit != 0) {
	blk_in->data[blk_in->byte] |= data[0];
	memcpy(&blk_in->data[blk_in->byte+1], data+1, len-1);
    } else {
	memcpy(&blk_in->data[blk_in->byte], data, len);
    }

    /* Do the decoding */
    do {
	if (!cs)
	    cs = restore_codes(blk_in, &bfinal);
	out = huffman_multi_decode(blk_in, cs);
	if (!bfinal) {
	    fprintf(stderr, "Multiple huffman blocks needed in "
		    "ideflate_decode: not catered for.\n");
	    block_destroy(out, 0);
	    return NULL;
	}
	huffman_codeset_destroy(cs);
	cs = NULL;

	break;
    } while (!bfinal);

    if (new_data)
	*new_data = data + blk_in->byte + (blk_in->bit != 0);

    return out;
}
#endif

/* ------------------------------------------------------------------------ */
/*
 * Simple interfaces to the underlying g-library. These are very basic, but
 * sufficient. Will need to rework this with proper lock utilisation schemes
 * for Gap4 proper.
 */

/* Hacky allocation scheme - always increment the record number */
static int allocate(g_io *io) {
    /* FIXME: we should track free records using a freerecs bitmap
     * as is done in Gap4.
     */
    static int record = -1;

    if (record == -1)
	record = io->gdb->gfile->header.num_records;

    return record++;
}

static GView lock(g_io *io, int rec, int mode) {
    if (!mode)
	mode = G_LOCK_EX;

    return g_lock_N_(io->gdb, io->client, 0, rec, mode);
}

static int unlock(g_io *io, GView v) {
    return g_unlock_(io->gdb, io->client, v);
}

static int g_write(g_io *io, GView v, void *buf, size_t len) {
    return g_write_(io->gdb, io->client, v, buf, len);
}

/*
 * Writes buf as a series of little-endian int32_t elements, byte
 * swapping if appropriate.
 */
static int g_write_le4(g_io *io, GView v, void *buf, size_t len) {
    static GCardinal le_a[256];
    int32_t *le = le_a;
    int ret;
    size_t i, j;

    if (len % sizeof(int32_t) != 0)
	return -1;

    if (len > 256 * sizeof(int32_t))
	le = (int32_t *)malloc(len * sizeof(int32_t));

    for (i = j = 0; i < len; i+=sizeof(int32_t), j++) {
	le[j] = le_int4(((int32_t *)buf)[j]);
    }

    ret = g_write_(io->gdb, io->client, v, le, len);

    if (le != le_a)
	free(le);
    return ret;
}

static int g_read(g_io *io, GView v, void *buf, size_t len) {
    return g_read_(io->gdb, io->client, v, buf, len);
}

static void *g_read_alloc(g_io *io, GView v, size_t *len) {
    GViewInfo vi;
    void *buf;

    g_view_info_(io->gdb, io->client, v, &vi);
    if (len)
	*len = vi.used;

    if (NULL == (buf = malloc(vi.used)))
	return NULL;
    
    if (g_read_(io->gdb, io->client, v, buf, vi.used) == 0)
	return buf;
    else {
	free(buf);
	return NULL;
    }
}

/*
 * Reads buf as a series of little-endian int32_t elements, byte
 * swapping if appropriate.
 */
static int g_read_le4(g_io *io, GView v, void *buf, size_t len) {
    int ret;
    size_t i, j;
    int32_t *le = (int32_t *)buf;

    if (len % sizeof(int32_t) != 0)
	return -1;

    ret = g_read_(io->gdb, io->client, v, buf, len);

    for (i = j = 0; i < len; i+=sizeof(int32_t), j++) {
	le[j] = le_int4(((int32_t *)buf)[j]);
    }

    return ret;
}

static int g_flush(g_io *io, GView v) {
    return g_flush_(io->gdb, io->client, v);
}


/* ------------------------------------------------------------------------
 * Generic io_database functions
 */


/* Generic functions - shared by all objects within the g library */
static GRec io_generic_create(void *dbh, void *unused) {
    return allocate((g_io *)dbh);
}

static int io_generic_destroy(void *dbh, GRec r) {
    fprintf(stderr, "io_generic_free unimplemented yet\n");
    return -1;
}

static GView io_generic_lock(void *dbh, GRec r, int mode) {
    return lock((g_io *)dbh, r, mode);
}

static int io_generic_unlock(void *dbh, GView v) {
    g_flush((g_io *)dbh, v);
    return unlock((g_io *)dbh, v);
}

static int io_generic_upgrade(void *dbh, GView v, int mode) {
    g_io *gio = (g_io *)dbh;
    return g_upgrade_(gio->gdb, gio->client, v, mode);
}

static int io_generic_abandon(void *dbh, GView v) {
    g_io *gio = (g_io *)dbh;
    return g_abandon_(gio->gdb, gio->client, v);
}

static cached_item *io_generic_read(void *dbh, GRec rec, int type) {
    GView v;
    void *buf;
    size_t buf_len;
    cached_item *ci;

    /* Load from disk */
    if (-1 == (v = io_generic_lock(dbh, rec, G_LOCK_RO)))
	return NULL;

    buf = g_read_alloc((g_io *)dbh, v, &buf_len);

    if (!(ci = cache_new(type, rec, v, NULL, buf_len))) {
	free(buf);
        return NULL;
    }

    memcpy(&ci->data, buf, buf_len);
    free(buf);

    ci->data_size = buf_len;
    return ci;
}

static int io_generic_write(void *dbh, cached_item *ci) {
    int ret;

    assert(ci->lock_mode >= G_LOCK_RW);

    ret = g_write((g_io *)dbh, ci->view, &ci->data, ci->data_size);
    g_flush((g_io *)dbh, ci->view); /* Should we auto-flush? */
    return ret;
}

static int io_generic_info(void *dbh, GView v, GViewInfo *vi) {
    g_io *io = (g_io *)dbh;
    return g_view_info_(io->gdb, io->client, v, vi);
}

/* ------------------------------------------------------------------------
 * The B+Tree cache methods
 *
 * Our B+Tree keys are sequence names, contig names, etc.
 * The data attached to these is the corresponding record number for the
 * object.
 *
 * However do not confuse those with the HacheTable keys and data.
 * They represent the record number (key) for a BTree block and the
 * in-memory representation of that btree block.
 */

typedef struct {
    g_io *io;
    HacheTable *h;
} btree_query_t;

static HacheData *btree_load_cache(void *clientdata, char *key, int key_len,
				   HacheItem *hi) {
    btree_query_t *bt = (btree_query_t *)clientdata;
    g_io *io = bt->io;
    HacheTable *h = bt->h;
    cached_item *ci;
    GView v;
    size_t len;
    char *buf;
    btree_node_t *n;
    BTRec rec = *((BTRec *)key);
    static HacheData hd;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    if (NULL == (buf = g_read_alloc(io, v, &len))) {
	return NULL;
    }

    /* Decode the btree element */
    n = btree_node_decode(buf);
    n->rec = rec;

    ci = cache_new(GT_BTree, rec, v, NULL, sizeof(btree_node_t *));
    ci->data = n;
    n->cache = ci;

    //printf("btree_load_cache(%d)\n", n->rec);

    free(buf);

    hd.p = ci;
    ci->hi = hi;

    HacheTableDecRef(h, hi);

    return &hd;
}

static void btree_del_cache(void *clientdata, HacheData hd) {
    g_io *io = (g_io *)clientdata;
    cached_item *ci = hd.p;
    btree_node_t *n = (btree_node_t *)ci->data;

    assert(ci->updated == 0);
    
    unlock(io, ci->view);
    free(ci);
    //printf("btree_del_cache(%d)\n", n->rec);
    btree_del_node(n);
}

static int btree_write(g_io *io, btree_node_t *n) {
    int ret;
    size_t len;
    GView v;
    char *data = (char *)btree_node_encode(n, &len);
    cached_item *ci = n->cache;

    if (ci) {
	ret = g_write(io, ci->view, data, len);
	g_flush(io, ci->view);
    } else {
	if (-1 == (v = lock(io, n->rec, G_LOCK_EX))) {
	    fprintf(stderr, "Failed to lock btree node %d\n", n->rec);
	    return -1;
	}
	ret = g_write(io, v, data, len);
	unlock(io, v);
    }

    free(data);

    if (ret) {
	fprintf(stderr, "Failed to write btree node %d\n", n->rec);
    }

    return ret;
}


/*
 * FIXME: These should be passed in as pointers to functions when creating
 * a btree so we can have multiple btrees using different access functions.
 * For now this suffices though.
 */
btree_node_t *btree_node_get(void *clientdata, BTRec r) {
    btree_query_t *bt = (btree_query_t *)clientdata;
    HacheItem *hi = HacheTableSearch(bt->h, (char *)&r, sizeof(r));

    return hi
	? (btree_node_t *)(((cached_item *)hi->data.p)->data)
	: (fprintf(stderr, "Failed to load btree %d\n", r), NULL);
}

/*
 * Create a node and add it to our internal cache too
 *
 * Returns the record number on success
 *         -1 on failure
 */
int btree_node_create(g_io *io, HacheTable *h) {
    GRec rec;
    btree_node_t *n;
    cached_item *ci;
    GView v;
    static HacheData hd;
    HacheItem *hi;

    /* Allocate a new record */
    rec = allocate(io);
    n = btree_new_node();
    n->rec = rec;

    /* Lock it and populate our hash */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return -1;

    ci = cache_new(GT_BTree, rec, v, NULL, sizeof(btree_node_t *));
    ci->data = n;
    n->cache = ci;

    hd.p = ci;

    hi = HacheTableAdd(h, (char *)&rec, sizeof(rec), hd, NULL);
    ci->hi = hi;

    HacheTableDecRef(h, hi);

    return rec;
}

btree_node_t *btree_node_new(void *cd) {
    btree_query_t *bt = (btree_query_t *)cd;
    int rec = btree_node_create(bt->io, bt->h);
    
    return btree_node_get(cd, rec);
}

int btree_node_put(void *cd, btree_node_t *n) {
    btree_query_t *bt = (btree_query_t *)cd;
    cached_item *ci = n->cache;

    if (!ci->updated) {
	if (-1 == g_upgrade_(bt->io->gdb, bt->io->client, ci->view, G_LOCK_RW))
	    return -1;
	ci->updated = 1;
	HacheTableIncRef(bt->h, ci->hi);
    }

    return 0;
}

void btree_node_del(void *cd, btree_node_t *n) {
    /* FIXME: deallocate disk storage space too */
    btree_del_node(n);
}

void btree_flush(g_io *io, HacheTable *h) {
    int i;

    if (!h)
	return;

    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    cached_item *ci = hi->data.p;
	    btree_node_t *n;

	    if (!ci->updated)
		continue;
	    
	    n = (btree_node_t *)(ci->data);
	    if (0 == btree_write(io, n)) {
		ci->updated = 0;
		/*
		 * FIXME - may need to move this out of the loop by creating
		 * a list of items to DecRef here and then working through
		 * them.
		 *
		 * Reason? Because DecRef itself has the capability of
		 * removing items and hence can alter the hi->next list
		 * we're iterating through.
		 */
		HacheTableDecRef(h, hi);
	    }
	}
    }

    //HacheOrderPurge(h);
}

void btree_destroy(g_io *io, HacheTable *h) {
    int i;

    if (!h)
	return;

    puts("\n=== btree_hash ===");
    HacheTableStats(h, stdout);

    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    btree_del_cache(io, hi->data);
	}
    }

    HacheTableDestroy(h, 0);
}

void btree_inc_ref(void *cd, btree_node_t *n) {
    btree_query_t *bt = (btree_query_t *)cd;
    cached_item *ci = n->cache;

    HacheTableIncRef(bt->h, ci->hi);
}

void btree_dec_ref(void *cd, btree_node_t *n) {
    btree_query_t *bt = (btree_query_t *)cd;
    cached_item *ci = n->cache;

    HacheTableDecRef(bt->h, ci->hi);
}

/* ------------------------------------------------------------------------
 * GDatabase access methods
 */
static int bitsize = G_64BIT;
static int io_database_create_files(char *fn) {
    char auxfn[1024];
    int fd;
    AuxHeader auxheader;
    GCardinal i;
    int (*(*low_level_vector))(int fd, void *x, int num);
    int endian = 1;
    dheap_t *h;

    /*
     * Determine the default vectors for creating a database.
     * This is temporary as from here on the g-library will auto-sense when
     * it opens the database files.
     */
    if ( *(char *)&endian ) {
	low_level_vector = (bitsize == G_64BIT)
	    ? low_level_vectors_swapped64
	    : low_level_vectors_swapped32;
    } else {
	low_level_vector = (bitsize == G_64BIT)
	    ? low_level_vectors64
	    : low_level_vectors32;
    }

    /* check file name isn't too long */
    if ( strlen(fn) + strlen(G_AUX_SUFFIX) >= sizeof(auxfn) ) return gerr_set(GERR_NAME_TOO_LONG);
	
    strcpy(auxfn,fn);
    strcat(auxfn,G_AUX_SUFFIX);

    /* create files */
    /* LOW LEVEL IO HERE */
    if (NULL == (h = heap_create(fn)))
	return gerr_set(GERR_CANT_CREATE);
    // heap_destroy(h);

    /* LOW LEVEL IO HERE */
    if ( (fd = creat(auxfn,G_DEF_PERMS)) == -1 ) return gerr_set(GERR_CANT_CREATE);

    /* initialise header */
    auxheader.file_size = 0;
    auxheader.block_size = (int32_t)8;
    auxheader.num_records = 0;
    auxheader.max_records = 100; /* dynamically grows anyway */
    auxheader.last_time = G_YEAR_DOT;
    auxheader.flags = (GFlags) 0;
    auxheader.spare1 = (GFlags) 0;
    auxheader.free_time = G_YEAR_DOT;
    for (i=G_Number(auxheader.spare)-1;i>=0;i--) auxheader.spare[i]=0;
    auxheader.format = bitsize;

    /* write(fd,&auxheader,sizeof(auxheader)); */
    (void) (low_level_vector[GOP_WRITE_AUX_HEADER])(fd,&auxheader,1);

    /* LOW LEVEL IO HERE */
    close(fd);

    return 0;
}

static void *io_database_connect(char *dbname, int ro) {
    g_io *io;
    btree_query_t *bt;

    if (NULL == (io = (g_io *)calloc(1, sizeof(*io))))
	return NULL;

    io->gdb = g_open_database_(&dbname, 1, ro);
    io->client = g_connect_client_(io->gdb, 0, G_LOCK_EX, &io->mode);
#ifdef INDEX_NAMES
    io->seq_name_hash = HacheTableCreate(10000,
					 HASH_DYNAMIC_SIZE | HASH_OWN_KEYS);

    if (NULL == (bt = (btree_query_t *)malloc(sizeof(*bt))))
	return NULL;
    bt->io = io;
    bt->h = io->seq_name_hash;
    io->seq_name_hash->clientdata = bt;
    io->seq_name_hash->load = btree_load_cache;
    io->seq_name_hash->del  = btree_del_cache;
#else
    io->seq_name_hash = NULL;
#endif

    io->contig_name_hash = HacheTableCreate(1000,
					    HASH_DYNAMIC_SIZE | HASH_OWN_KEYS);

    if (NULL == (bt = (btree_query_t *)malloc(sizeof(*bt))))
	return NULL;
    bt->io = io;
    bt->h = io->contig_name_hash;
    io->contig_name_hash->clientdata = bt;
    io->contig_name_hash->load = btree_load_cache;
    io->contig_name_hash->del  = btree_del_cache;

    io->seq_name_tree = NULL; /* Initialised when reading GDatabase */
    io->contig_name_tree = NULL;

    return io;
}

static int io_database_commit(void *dbh) {
    g_io *io = (g_io *)dbh;

    btree_flush(io, io->seq_name_hash);
    btree_flush(io, io->contig_name_hash);
    
    return 0;
}

static int io_database_disconnect(void *dbh) {
    g_io *io = (g_io *)dbh;

    //btree_print(io->seq_name_tree, io->seq_name_tree->root, 0);

    io_database_commit(dbh);
    btree_destroy(io, io->seq_name_hash);
    btree_destroy(io, io->contig_name_hash);

    g_disconnect_client_(io->gdb, io->client);
    g_shutdown_database_(io->gdb);

    free(io);

    return 0;
}

static GRec io_database_create(void *dbh, void *from) {
    g_io *io = (g_io *)dbh;
    GCardinal db_rec = allocate(io);
    GView v;
    GDatabase db;

    /* init_db is only called on a blank database => first record is 0 */
    assert(db_rec == 0);

    db.Ncontigs = 0;
    db.version = 0;

    /* Contig order */
    db.contig_order = allocate(io); /* contig array */
    v = lock(io, db.contig_order, G_LOCK_EX);
    if (-1 == g_write_le4(io, v, NULL, 0))
	return -1;
    g_flush(io, v);
    unlock(io, v);

#ifdef INDEX_NAMES
    /* Sequence name btree */
    db.seq_name_index = btree_node_create(io, io->seq_name_hash);
    assert(db.seq_name_index > 0);
#else
    db.seq_name_index = 0;
#endif
    db.contig_name_index = btree_node_create(io, io->contig_name_hash);

    /* Database struct itself */
    v = lock(io, db_rec, G_LOCK_EX);
    if (-1 == g_write_le4(io, v, (void *)&db, sizeof(db)))
	return -1;
    g_flush(io, v);
    unlock(io, v);
   
    return 0;
}

static cached_item *io_database_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    GDatabase *db;
    cached_item *ci;

    if (NULL == (ci = io_generic_read(dbh, rec, GT_Database)))
	return NULL;
    db = (GDatabase *)&ci->data;

#ifdef INDEX_NAMES
    /* Initialise the seq_name btree if needed */
    if (io->seq_name_tree)
	return ci;

    /* Read the root */
    if (db->seq_name_index) {
	btree_query_t *bt = (btree_query_t *)malloc(sizeof(*bt));
	bt->io = io;
	bt->h = io->seq_name_hash;
	io->seq_name_tree = btree_new(bt, db->seq_name_index);
	assert(io->seq_name_tree);
	assert(io->seq_name_tree->root);
    }

    //printf("seq_name_hash=%p\n", io->seq_name_hash);
#endif

    if (db->contig_name_index) {
	btree_query_t *bt = (btree_query_t *)malloc(sizeof(*bt));
	bt->io = io;
	bt->h = io->contig_name_hash;
	io->contig_name_tree = btree_new(bt, db->contig_name_index);
	assert(io->contig_name_tree);
	assert(io->contig_name_tree->root);
    }

    return ci;

}

/* ------------------------------------------------------------------------
 * contig_t access methods
 */
static cached_item *io_contig_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    GContig_header *ch;
    contig_t *c;
    size_t len, slen;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    if (NULL == (ch = g_read_alloc(io, v, &len))) {
	free(ci);
	return NULL;
    }

    /* Generate in-memory data structure */
    slen = sizeof(*ch) + sizeof(char *) + ((unsigned char *)(ch+1))[0];
    if (!(ci = cache_new(GT_Contig, rec, v, NULL, slen)))
	return NULL;

    c = (contig_t *)&ci->data;
    c->rec    = rec;
    c->start  = le_int4(ch->start);
    c->end    = le_int4(ch->end);
    c->bin    = le_int4(ch->bin);
    c->name   = (char *)(&c->name + 1);
    strncpy(c->name, ((unsigned char *)ch) + sizeof(GContig_header) + 1,
	    ((unsigned char *)ch)[sizeof(GContig_header)]);
    c->name[((unsigned char *)ch)[sizeof(GContig_header)]] = 0;

    free(ch);

    return ci;
}

static int io_contig_write_view(g_io *io, contig_t *c, GView v) {
    GContig_header *ch;
    size_t len;

    /* Construct on-disc representation */
    len = sizeof(*ch) + 1 + (c->name ? strlen(c->name) : 0);
    ch = (GContig_header *)malloc(len);
    ch->start = le_int4(c->start);
    ch->end   = le_int4(c->end);
    ch->bin = le_int4(c->bin);
    ((unsigned char *)ch)[sizeof(GContig_header)] =
	c->name ? strlen(c->name) : 0;
    if (c->name)
	strncpy(&((char *)ch)[sizeof(GContig_header)+1], c->name,
		((unsigned char *)ch)[sizeof(GContig_header)]);

    /* Write the data */
    if (-1 == g_write(io, v, (char *)ch, len)) {
	free(ch);
	return -1;
    }

    g_flush(io, v);
    free(ch);
    return 0;
}

static int io_contig_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    contig_t *c = (contig_t *)&ci->data;

    assert(ci->lock_mode >= G_LOCK_RW);
    return io_contig_write_view(io, c, ci->view);
}

static GRec io_contig_create(void *dbh, void *vfrom) {
    g_io *io = (g_io *)dbh;
    GRec rec;
    GView v;
    contig_t *from = (contig_t *)vfrom;

    rec = allocate(io);
    v = lock(io, rec, G_LOCK_EX);
    
    if (from) {
	io_contig_write_view(io, from, v);
    } else {
	contig_t c;
	c.rec = rec;
	c.start = c.end = 0;
	c.bin = 0;
	c.name = NULL;
	io_contig_write_view(io, &c, v);
    }
    unlock(io, v);

    return rec;
}

static GRec io_contig_index_query(void *dbh, char *name) {
    g_io *io = (g_io *)dbh;
    
    if (!io->contig_name_tree)
	return -1;

    return btree_search(io->contig_name_tree, name);
}

static int io_contig_index_add(void *dbh, char *name, GRec rec) {
    g_io *io = (g_io *)dbh;
    
    if (!io->contig_name_tree)
	return -1;

    btree_insert(io->contig_name_tree, name, rec);
    return io->contig_name_tree->root->rec;
}

/* ------------------------------------------------------------------------
 * Array access methods
 */
static cached_item *io_array_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    GViewInfo vi;
    ArrayStruct *ar;
    int ar_sz;

    if (!(ci = cache_new(GT_RecArray, rec, 0, NULL, sizeof(*ar))))
	return NULL;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    g_view_info_(io->gdb, io->client, v, &vi);
    ar_sz = vi.used / sizeof(GCardinal);

    ar = ArrayCreate(sizeof(GCardinal), ar_sz);
    ArrayRef(ar, ar_sz);
    g_read(io, v, ArrayBase(GCardinal, ar), ar_sz * sizeof(GCardinal));

    memcpy(&ci->data, ar, sizeof(*ar));
    free(ar);

    ci->view = v;
    return ci;
}

static int io_array_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    Array ar;
    int ret;

    assert(ci->lock_mode >= G_LOCK_RW);
    ar = (Array)&ci->data;
    ret = g_write(io, ci->view, ArrayBase(GCardinal, ar), 
		   ArrayMax(ar) * sizeof(GCardinal));
    g_flush(io, ci->view);

    return ret;
}


/* ------------------------------------------------------------------------
 * anno access methods
 */
static cached_item *io_anno_read(void *dbh, GRec rec) {
    return io_generic_read(dbh, rec, 0);
}


/* ------------------------------------------------------------------------
 * dnasrc access methods
 */
static cached_item *io_dnasrc_read(void *dbh, GRec rec) {
    return io_generic_read(dbh, rec, GT_DNASource);
}


/* ------------------------------------------------------------------------
 * vector access methods
 */
static cached_item *io_vector_read(void *dbh, GRec rec) {
    return io_generic_read(dbh, rec, 0);
}


/* ------------------------------------------------------------------------
 * bin access methods
 */
static cached_item *io_bin_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    GBin *b;
    bin_index_t *bin;
    void *buf;
    size_t buf_len;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    if (NULL == (buf = g_read_alloc(io, v, &buf_len)))
	return NULL;
    if (buf_len < sizeof(GBin)) {
	buf = realloc(buf, sizeof(GBin));
	memset(buf + buf_len, 0, sizeof(GBin) - buf_len);
    }
    b = (GBin *)buf;

    /* Allocate our overlapping data objects */
    if (!(ci = cache_new(GT_Bin, rec, v, NULL, sizeof(*bin))))
	return NULL;
    bin = (bin_index_t *)&ci->data;

    /* Construct bin */
    bin->rec         = rec;
    bin->pos         = b->pos;
    bin->size        = b->size;
    bin->start_used  = b->start;
    bin->end_used    = b->end;
    bin->child[0]    = b->child[0];
    bin->child[1]    = b->child[1];
    bin->bin_id      = b->id;
    bin->flags       = b->flags;
    bin->parent      = b->parent;
    bin->parent_type = b->parent_type;
    bin->rng_rec     = b->range;
    bin->rng         = NULL;
    bin->track_rec   = b->track;
    bin->track       = NULL;

    /* Load ranges */
    if (b->range) {
	GViewInfo vi;
	int nranges;

	v = lock(io, b->range, G_LOCK_RO);
	g_view_info_(io->gdb, io->client, v, &vi);
	nranges = vi.used / sizeof(GRange);

	bin->rng = ArrayCreate(sizeof(GRange), nranges);
	ArrayRef(bin->rng, nranges-1);
	g_read(io, v, ArrayBase(GRange, bin->rng),
	       nranges * sizeof(GRange));
	unlock(io, v);
    }

    /* Load tracks */
    if (b->track) {
	GViewInfo vi;
	int ntracks;

	v = lock(io, b->track, G_LOCK_RO);
	g_view_info_(io->gdb, io->client, v, &vi);
	ntracks = vi.used / sizeof(GBinTrack);

	bin->track = ArrayCreate(sizeof(GBinTrack), ntracks);
	ArrayRef(bin->track, ntracks-1);
	g_read(io, v, ArrayBase(GBinTrack, bin->track),
	       ntracks * sizeof(GBinTrack));
	unlock(io, v);
    }

    free(buf);
    return ci;
}

static int io_bin_write_view(g_io *io, bin_index_t *bin, GView v) {
    GBin g;
    int err = 0;

    /* Ranges */
    if (bin->flags & BIN_RANGE_UPDATED) {
	GView v;

	bin->flags &= ~BIN_RANGE_UPDATED;

	if (!bin->rng_rec) {
	    bin->rng_rec = allocate(io);
	    bin->flags |= BIN_BIN_UPDATED;
	}

	v = lock(io, bin->rng_rec, G_LOCK_EX);
	err |= g_write(io, v, ArrayBase(GRange, bin->rng),
		       sizeof(GRange) * ArrayMax(bin->rng));
	err |= unlock(io, v);
    }

    /* Tracks */
    if (bin->flags & BIN_TRACK_UPDATED) {
	GView v;

	bin->flags &= ~BIN_TRACK_UPDATED;

	if (!bin->track_rec) {
	    bin->track_rec = allocate(io);
	    bin->flags |= BIN_BIN_UPDATED;
	}

	v = lock(io, bin->track_rec, G_LOCK_EX);
	err |= g_write(io, v, ArrayBase(GBinTrack, bin->track),
		       sizeof(GBinTrack) * ArrayMax(bin->track));
	err |= unlock(io, v);
    }

    /* Bin struct itself */
    if (bin->flags & BIN_BIN_UPDATED) {
	bin->flags &= ~BIN_BIN_UPDATED;
	g.pos         = bin->pos;
	g.size        = bin->size;
	g.start       = bin->start_used;
	g.end         = bin->end_used;
	g.id          = bin->bin_id;
	g.flags       = bin->flags;
	g.parent      = bin->parent;
	g.parent_type = bin->parent_type;
	g.child[0]    = bin->child[0];
	g.child[1]    = bin->child[1];
	g.range       = bin->rng_rec;
	g.track       = bin->track_rec;

	err |= g_write(io, v, &g, sizeof(g));
	g_flush(io, v);
    }

    return err;
}

static int io_bin_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    bin_index_t *bin = (bin_index_t *)&ci->data;

    assert(ci->lock_mode >= G_LOCK_RW);
    return io_bin_write_view(io, bin, ci->view);
}

static int io_bin_create(void *dbh, void *vfrom) {
    bin_index_t *from = vfrom;
    g_io *io = (g_io *)dbh;
    GRec rec;
    GView v;

    rec = allocate(io);
    v = lock(io, rec, G_LOCK_EX);
    
    if (from) {
	io_bin_write_view(io, from, v);
    } else {
	bin_index_t b;
	b.pos         = 0;
	b.size        = 0;
	b.start_used  = 0;
	b.end_used    = 0;
	b.parent      = 0;
	b.parent_type = 0;
	b.child[0]    = 0;
	b.child[1]    = 0;
	b.bin_id      = rec;
	b.rng	      = NULL;
	b.rng_rec     = 0;
	b.track	      = NULL;
	b.track_rec   = 0;
	b.flags       = 0;
	io_bin_write_view(io, &b, v);
    }
    unlock(io, v);

    return rec;
}


/* ------------------------------------------------------------------------
 * track access methods
 */
static cached_item *io_track_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    GTrack_Header *t;
    track_t *track;
    char *buf;
    size_t buf_len;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    if (NULL == (buf = g_read_alloc(io, v, &buf_len)))
	return NULL;
    t = (GTrack_Header *)buf;

    /* Allocate our overlapping data objects */
    if (!(ci = cache_new(GT_Track, rec, v, NULL, sizeof(*track) + buf_len)))
	return NULL;
    track = (track_t *)&ci->data;

    /* Construct track_t */
    track->rec       = rec;
    track->type      = t->type;
    track->flag      = t->flags;
    track->item_size = t->item_size;
    track->nitems    = t->nitems;
    track->data      = ArrayCreate(track->item_size, track->nitems);
    memcpy(ArrayBase(char, track->data), buf + sizeof(GTrack_Header),
	   t->item_size * t->nitems);

    free(buf);
    return ci;
}

static int io_track_write_view(g_io *io, track_t *track, GView v) {
    GTrack_Header *h;
    char *data;
    int err = 0;

    data = (char *)malloc(sizeof(*h) + track->item_size * track->nitems);
    if (!data)
	return -1;
    h = (GTrack_Header *)data;

    /* The structure header fields */
    h->type      = track->type;
    h->flags     = track->flag & ~TRACK_FLAG_FREEME; /* not fake */
    h->item_size = track->item_size;
    h->nitems    = track->data ? track->nitems : 0;
    
    /* The array */
    if (h->nitems)
	memcpy(data + sizeof(*h), ArrayBase(char, track->data),
	       track->item_size * track->nitems);
    
    err |= g_write(io, v, h, sizeof(*h) + track->item_size * track->nitems);
    g_flush(io, v);
    free(data);

    return err;
}

static int io_track_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    track_t *track = (track_t *)&ci->data;

    printf("Write track %d\n", ci->rec);

    assert(ci->lock_mode >= G_LOCK_RW);
    return io_track_write_view(io, track, ci->view);
}

static int io_track_create(void *dbh, void *vfrom) {
    track_t *from = vfrom;
    g_io *io = (g_io *)dbh;
    GRec rec;
    GView v;

    rec = allocate(io);
    v = lock(io, rec, G_LOCK_EX);
    
    if (from) {
	io_track_write_view(io, from, v);
    } else {
	track_t t;
	t.type = TRACK_UNKNOWN;
	t.flag = 0;
	t.item_size = 0;
	t.nitems = 0;
	t.rec = 0;
	t.data = NULL;
	io_track_write_view(io, &t, v);
    }
    unlock(io, v);

    return rec;
}


/* ------------------------------------------------------------------------
 * seq_t access methods
 */

/*
 * Decodes an on-disk sequence structure into a malloced seq_t struct.
 *
 * On disc struct:
 * ? bytes bin record no.
 * ? byte 'left clip'
 * ? byte 'right clip'
 * ? byte sequence length
 * ? byte other_end (rec.no. of paired read)
 * ? byte parent_rec
 * 1 byte parent_type
 * 1 byte seq_tech (3 bottom bits)
 *      + flags (3 next bits)
 *      + format (2 top bits)
 * 1 byte mapping_quality
 * ? bytes name, nul terminated
 * ? byte trace name, nul terminated
 * ? bytes alignment, nul terminated
 * remainder is seq/qual (various formats).
 *
 * Returns a pointer to a seq_t struct
 *      or NULL on failure.
 */
static cached_item *seq_decode(unsigned char *buf, size_t len) {
    cached_item *ci;
    unsigned char *cp, flags, mapping_qual, seq_tech, format;
    size_t slen;
    signed int i, j;
    seq_t *seq;
    uint32_t left, right, bin, seq_len;
    int parent_type, parent_rec, other_end;

    cp = buf;
    cp += u72int(cp, &bin);
    cp += u72int(cp, &left);
    cp += u72int(cp, &right);
    cp += u72int(cp, &seq_len);
    cp += u72int(cp, &other_end);
    cp += u72int(cp, &parent_rec);
    parent_type = *cp++;
    format = *cp++;
    seq_tech = format & ((1<<3)-1);
    format >>= 3;
    flags = format & ((1<<3)-1);
    format >>= 3;
    mapping_qual = *cp++;
    /* cp is now the variable sized section starting with reading name */

    /* Generate in-memory data structure */
    slen = sizeof(seq_t) + len - (cp-buf) +
	seq_len * (1 + (format == SEQ_FORMAT_CNF4 ? 4 : 1));

    if (!(ci = cache_new(GT_Seq, 0, 0, NULL, slen)))
        return NULL;
    seq = (seq_t *)&ci->data;

    seq->bin          = bin;
    seq->left         = left;
    seq->right        = right;
    seq->parent_type  = parent_type;
    seq->parent_rec   = parent_rec;
    seq->other_end    = other_end;
    seq->seq_tech     = seq_tech;
    seq->flags        = flags;
    seq->format       = format;
    seq->mapping_qual = mapping_qual;
    seq->len          = (seq->flags & SEQ_COMPLEMENTED)
	                    ? -seq_len : seq_len;

    memcpy(&seq->data, cp, len - (cp-buf));

    /* Name */
    seq->name = (char *)&seq->data;
    seq->name_len = strlen(seq->name);
    cp += seq->name_len + 1;

    /* Trace name */
    seq->trace_name = seq->name + seq->name_len + 1;
    seq->trace_name_len = strlen(seq->trace_name);
    cp += seq->trace_name_len + 1;

    /* Alignment */
    seq->alignment = seq->trace_name + seq->trace_name_len + 1;
    seq->alignment_len = strlen(seq->alignment);
    cp += seq->alignment_len + 1;

    /* Seq/Qual */
    seq->seq = seq->alignment + seq->alignment_len + 1;
    seq->conf = seq->seq + seq_len;

    /* cp is now at the encoded seq/qual area */
    switch (seq->format) {
    case SEQ_FORMAT_MAQ:
	for (i = 0; i < seq_len; i++) {
	    seq->conf[i] = cp[i] >> 2;
	    seq->seq[i] = seq->conf[i] ? "ACGT"[cp[i] & 3] : 'N';
	}
	break;

    case SEQ_FORMAT_CNF1:
#if 1
	for (i = j = 0; i < seq_len;) {
	    unsigned char c = cp[j++];
	    seq->seq[i++] = "ACGTN*"[c%6]; c /= 6;
	    if (i >= seq_len) break;
	    seq->seq[i++] = "ACGTN*"[c%6]; c /= 6;
	    if (i >= seq_len) break;
	    seq->seq[i++] = "ACGTN*"[c%6];
	}
	cp += j;
	for (i = j = 0; i < seq_len;) {
	    seq->conf[i++] = *cp++;
	    if (i+1 < seq_len && cp[-1] == cp[0]) {
		unsigned char dup = *cp++;
		for (j = *cp++; j >= 0; j--)
		    seq->conf[i++] = dup;
	    }
	}
#else
	for (i = 0; i < seq_len; i++) {
	    seq->seq[i] = cp[i];
	}
	cp += seq_len;
	for (i = 0; i < seq_len; i++) {
	    seq->conf[i] = cp[i];
	}
#endif
	break;

    case SEQ_FORMAT_CNF4:
	for (i = 0; i < seq_len; i++) {
	    seq->seq[i] = cp[i];
	}
	cp += seq_len;
	for (i = 0; i < seq_len*4; i++) {
	    seq->conf[i] = cp[i];
	}
	break;

    default:
	fprintf(stderr, "Unknown sequence format '%d'\n", seq->format);
	break;
    }


    return ci;
}

static cached_item *io_seq_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    void *bloc;
    size_t bloc_len;
    GView v;
    cached_item *ci;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    bloc = g_read_alloc(io, v, &bloc_len);

    if (!bloc)
	return NULL;

    ci = seq_decode(bloc, bloc_len);
    free(bloc);

    if (!ci)
	return NULL;

    ci->view = v;
    ci->rec = rec;
    return ci;
}

/* See seq_decode for the storage format */
static int io_seq_write_view(g_io *io, seq_t *seq, GView v, GRec rec) {
    int err = 0;
    int i, j, seq_len, name_len, trace_name_len, data_len;
    unsigned char block[1024], *cp = block, *cpstart;
    static unsigned char base2val_maq[256] = { /* ACGT => 0123, else 9 */
     /* 0 1 2 3 4 5 6 7 8 9 A B C D E F */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* 00 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* 10 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* 20 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* 30 */
	9,0,9,1,9,9,9,2,9,9,9,9,9,9,9,9, /* 40 */
	9,9,9,9,3,3,9,9,9,9,9,9,9,9,9,9, /* 50 */
	9,0,9,1,9,9,9,2,9,9,9,9,9,9,9,9, /* 60 */
	9,9,9,9,3,3,9,9,9,9,9,9,9,9,9,9, /* 70 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* 80 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* 20 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* A0 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* B0 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* C0 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* D0 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, /* E0 */
	9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9  /* F0 */
    };
    static unsigned char base2val_cnf1[256] = { /* ACGTN* => 012345, else 4 */
     /* 0 1 2 3 4 5 6 7 8 9 A B C D E F */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, /* 00 */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, /* 10 */
	4,4,4,4,4,4,4,4,4,4,5,4,4,4,4,4, /* 20 */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, /* 30 */
	4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, /* 40 */
	4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4, /* 50 */
	4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4, /* 60 */
	4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4, /* 70 */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, /* 80 */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, /* 20 */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, /* A0 */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, /* B0 */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, /* C0 */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, /* D0 */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, /* E0 */
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4  /* F0 */
    };

    name_len = seq->name_len;
    if (name_len > 255)
	name_len = 255;

    trace_name_len = seq->trace_name_len;
    if (trace_name_len > 255)
	trace_name_len = 255;

    seq_len = ABS(seq->len);

    /* Auto-detect format where possible */
    if (seq->format == 0) {
	seq->format = SEQ_FORMAT_MAQ;
	for (i = 0; i < seq_len; i++) {
	    if (seq->seq[i] != 'A' && seq->seq[i] != 'C' &&
		seq->seq[i] != 'G' && seq->seq[i] != 'T') {
		seq->format = SEQ_FORMAT_CNF1;
		break;
	    }
	}
    }

    /* Worst case */
    data_len =
	  5 /* bin rec.no */
	+ 5 /* left clip */
	+ 5 /* right clip */
	+ 5 /* seq len */
	+ 5 /* other_end */
	+ 5 /* parent_rec */
	+ 1 /* parent_type */
	+ 1 /* seq_tech/flags/format */
	+ 1 /* mapping_quality */
	+ name_len + 1 /* name */
	+ trace_name_len + 1 /* trace name */
	+ seq_len*5 + 2; /* deflate worst expansion? */
    if (data_len > 1024) {
	if (NULL == (cp = (unsigned char *)malloc(data_len)))
	    return -1;
    }

    /* Clips */
    cpstart = cp;
    cp += int2u7(seq->bin, cp);
    cp += int2u7(seq->left, cp);
    cp += int2u7(seq->right, cp);
    cp += int2u7(seq_len, cp);

    /* Read-pair info */
    cp += int2u7(seq->other_end, cp);
    cp += int2u7(seq->parent_rec, cp);
    *cp++ = seq->parent_type;

    /* flags & m.quality */
    *cp++ = (seq->format << 6) | (seq->flags << 3) | seq->seq_tech;
    *cp++ = seq->mapping_qual;

    /* Name */
    strcpy(cp, seq->name);
    cp += seq->name_len + 1;

    /* Trace name */
    strcpy(cp, seq->trace_name);
    cp += seq->trace_name_len + 1;

    /* Alignment */
#if 0
  {
    char *al = seq->alignment;
    /* Squashed format */
    j = 0;
    while (*al) {
	unsigned int u = 0;
	unsigned int v = 0;
	switch (*al) {
	case 'M':
	    u = 0;
	    break;
	case 'I':
	    u = 1;
	    break;
	case 'D':
	    u = 2;
	    break;
	case 'E':
	    u = 3;
	    break;
	default:
	    fprintf(stderr, "Unknown alignment type '%c'\n",
		    seq->alignment[i]);
	    u = 3;
	}
	v = strtol(al+1, &al, 10);
	u |= v << 2;
	cp += int2u7(u, cp);
    }
    cp += int2u7(0, cp); /* match of length zero */
  }
#else
    strcpy(cp, seq->alignment);
    cp += seq->alignment_len + 1;
#endif

    /* Seq/Conf */
    switch (seq->format) {
    case SEQ_FORMAT_MAQ:
	for (i = 0; i < seq_len; i++) {
	    unsigned char v = base2val_maq[((unsigned char *)seq->seq)[i]];
	    if (v != 9) {
		if (seq->conf[i] <= 0)
		    v = 0;
		else if (seq->conf[i] >= 64)
		    v |= 63 << 2;
		else
		    v |= seq->conf[i] << 2;
	    } else {
		v = 0;
	    }
	    *cp++ = v;
	}
	break;

    case SEQ_FORMAT_CNF1:
#if 1
	for (i = j = 0; i < seq_len; i+=3, j++) {
	    unsigned char c1 = seq->seq[i];
	    unsigned char c2 = i+1 < seq_len ? seq->seq[i+1] : 0;
	    unsigned char c3 = i+2 < seq_len ? seq->seq[i+2] : 0;
	    *cp++ = base2val_cnf1[c1] + base2val_cnf1[c2]*6 + base2val_cnf1[c3]*36;
	}
	if (1) {
	    for (i = 0; i < seq_len; i = j) {
		j = i+1;
		while (j-i<257 && j < seq_len && seq->conf[i]==seq->conf[j])
		    j++;

		if (j-i == 1) {
		    *cp++ = seq->conf[i];
		} else {
		    *cp++ = seq->conf[i];
		    *cp++ = seq->conf[i];
		    *cp++ = j-i-2;
		}
	    }  
	}
#else
	for (i = 0; i < seq_len; i++) {
	    *cp++ = seq->seq[i];
	}
	for (i = 0; i < seq_len; i++) {
	    *cp++ = seq->conf[i];
	}
#endif
	break;

    case SEQ_FORMAT_CNF4:
	for (i = 0; i < seq_len; i++) {
	    *cp++ = seq->seq[i];
	}
	for (i = 0; i < seq_len*4; i++) {
	    *cp++ = seq->conf[i];
	}
	break;

    default:
	fprintf(stderr, "Unknown sequence format '%d'\n", seq->format);
	break;
    }

    //    printf("Write rec %d len %d %.*s:%d\n", rec, cp-cpstart,
    //	   seq->name_len, seq->name, seq->mapping_qual);
    err |= g_write(io, v, (void *)cpstart, cp-cpstart);
    g_flush(io, v);

    if (cpstart != block)
	free(cpstart);

    return err ? -1 : 0;
}

static int io_seq_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    seq_t *seq = (seq_t *)&ci->data;

    assert(ci->lock_mode >= G_LOCK_RW);
    return io_seq_write_view(io, seq, ci->view, ci->rec);
}

static int io_seq_create(void *dbh, void *vfrom) {
    seq_t *from = vfrom;
    g_io *io = (g_io *)dbh;
    GRec rec;
    GView v;

    rec = allocate(io);
    v = lock(io, rec, G_LOCK_EX);
    
    if (from) {
	io_seq_write_view(io, from, v, rec);
    } else {
	seq_t s;
	s.bin = 0;
	s.pos = 0;
	s.len = 0;
	s.left = s.right = 0;
	s.seq_tech = 0;
	s.flags = 0;
	s.parent_rec = 0;
	s.parent_type = 0;
	s.other_end = 0;
	s.mapping_qual = 0;
	s.name_len = 0;
	s.name = NULL;
	s.seq = NULL;
	s.conf = NULL;
	io_seq_write_view(io, &s, v, rec);
    }
    unlock(io, v);

    return rec;
}

static GRec io_seq_index_query(void *dbh, char *name) {
    g_io *io = (g_io *)dbh;
    
    if (!io->seq_name_tree)
	return -1;

    return btree_search(io->seq_name_tree, name);
}

static int io_seq_index_add(void *dbh, char *name, GRec rec) {
    g_io *io = (g_io *)dbh;
    
    if (!io->seq_name_tree)
	return -1;

    btree_insert(io->seq_name_tree, name, rec);
    return io->seq_name_tree->root->rec;
}

/* ------------------------------------------------------------------------
 * The externally visible object interfaces themselves.
 */
static iface iface_g = {
    io_database_create_files,
    io_database_connect,
    io_database_disconnect,
    io_database_commit,

    {
	/* Generic array */
	io_generic_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_array_read,
	io_array_write,
	io_generic_info,
    },

    {
	/* Database */
	io_database_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_database_read,
	io_generic_write,
	io_generic_info,
    },

    {
	/* Contig */
	io_contig_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_contig_read,
	io_contig_write,
	io_generic_info,
	io_contig_index_query,
	io_contig_index_add,
    },

    {
	/* Seq */
	io_seq_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_seq_read,
	io_seq_write,
	io_generic_info,
	io_seq_index_query,
	io_seq_index_add,
    },

    {
	/* Anno */
	io_generic_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_anno_read,
	io_generic_write,
	io_generic_info,
    },

    {
	/* DNASource */
	io_generic_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_dnasrc_read,
	io_generic_write,
	io_generic_info,
    },

    {
	/* Vector */
	io_generic_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_vector_read,
	io_generic_write,
	io_generic_info,
    },

    {
	/* Bin */
	io_bin_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_bin_read,
	io_bin_write,
	io_generic_info,
    },

    {
	/* Track */
	io_track_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_track_read,
	io_track_write,
	io_generic_info,
    },
};

iface *get_iface_g(void) {
    return &iface_g;
}
