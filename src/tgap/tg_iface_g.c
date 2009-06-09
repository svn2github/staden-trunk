#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <fcntl.h>

#include "g.h"
#include "misc.h"
#include "tg_iface_g.h"
#include "tg_utils.h"
#include "b+tree2.h"
#include "io_lib/deflate_interlaced.h"
#include "dna_utils.h"
#include "tg_gio.h"

static int wrstats[100];
static int wrcounts[100];
static int rdstats[100];
static int rdcounts[100];

#define INDEX_NAMES

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


/* ------------------------------------------------------------------------ */
/*
 * Data compression routines using zlib.
 */
#include <zlib.h>
static char *mem_deflate(char *data, size_t size, size_t *cdata_size) {
    z_stream s;
    char *cdata = NULL; /* Compressed output */
    int cdata_alloc = 0;
    int cdata_pos = 0;
    int err;

    cdata = malloc(cdata_alloc = size*1.05+10);
    cdata_pos = 0;

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)data;
    s.avail_in = size;
    s.total_in = 0;
    s.next_out  = cdata;
    s.avail_out = cdata_alloc;
    s.total_out = 0;
    s.data_type = Z_BINARY;

    err = deflateInit2(&s, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15,
		       9, Z_DEFAULT_STRATEGY);

    /* Encode to 'cdata' array */
    for (;s.avail_in;) {
	s.next_out = (unsigned char *)&cdata[cdata_pos];
	s.avail_out = cdata_alloc - cdata_pos;
	if (cdata_alloc - cdata_pos <= 0) {
	    fprintf(stderr, "Deflate produced larger output than expected. Abort\n"); 
	    return NULL;
	}
	err = deflate(&s, Z_NO_FLUSH);
	cdata_pos = cdata_alloc - s.avail_out;
	if (err != Z_OK) {
	    fprintf(stderr, "zlib deflate error: %s\n", s.msg);
	    break;
	}
    }
    if (deflate(&s, Z_FINISH) != Z_STREAM_END) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    *cdata_size = s.total_out;

    if (deflateEnd(&s) != Z_OK) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    return cdata;
}

static char *mem_deflate_parts(char *data,
			       size_t *part_size, int nparts,
			       size_t *cdata_size) {
    z_stream s;
    char *cdata = NULL; /* Compressed output */
    int cdata_alloc = 0;
    int cdata_pos = 0;
    int err, i;
    size_t size = 0;

    for (i = 0; i < nparts; i++)
	size += part_size[i];

    cdata = malloc(cdata_alloc = size*1.05+256*nparts);
    cdata_pos = 0;

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)data;
    s.total_in = 0;
    s.next_out  = cdata;
    s.avail_out = cdata_alloc;
    s.total_out = 0;
    s.data_type = Z_BINARY;

    err = deflateInit2(&s, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15,
		       9, Z_FILTERED);

    /* Encode to 'cdata' array */
    for (i = 0; i < nparts; i++) {
	s.avail_in = part_size[i];
	for (;s.avail_in;) {
	    s.next_out = (unsigned char *)&cdata[cdata_pos];
	    s.avail_out = cdata_alloc - cdata_pos;
	    if (cdata_alloc - cdata_pos <= 0) {
		fprintf(stderr, "Deflate produced larger output than expected. Abort\n"); 
		return NULL;
	    }
	    err = deflate(&s, Z_SYNC_FLUSH); // also try Z_FULL_FLUSH
	    cdata_pos = cdata_alloc - s.avail_out;
	    if (err != Z_OK) {
		fprintf(stderr, "zlib deflate error: %s\n", s.msg);
		break;
	    }
	}
    }
    if (deflate(&s, Z_FINISH) != Z_STREAM_END) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    *cdata_size = s.total_out;

    if (deflateEnd(&s) != Z_OK) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    return cdata;
}

static char *mem_deflate_lparts(char *data,
				size_t *part_size, int *level, int nparts,
				size_t *cdata_size) {
    z_stream s;
    char *cdata = NULL; /* Compressed output */
    int cdata_alloc = 0;
    int cdata_pos = 0;
    int err, i;
    size_t size = 0;

    for (i = 0; i < nparts; i++)
	size += part_size[i];

    cdata = malloc(cdata_alloc = size*1.05+256*nparts);
    cdata_pos = 0;

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)data;
    s.total_in = 0;
    s.next_out  = cdata;
    s.avail_out = cdata_alloc;
    s.total_out = 0;
    s.data_type = Z_BINARY;

    err = deflateInit2(&s, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15,
		       9, Z_FILTERED);

    /* Encode to 'cdata' array */
    for (i = 0; i < nparts; i++) {
	s.avail_in = part_size[i];
	for (;s.avail_in;) {
	    s.next_out = (unsigned char *)&cdata[cdata_pos];
	    s.avail_out = cdata_alloc - cdata_pos;
	    if (cdata_alloc - cdata_pos <= 0) {
		fprintf(stderr, "Deflate produced larger output than expected. Abort\n"); 
		return NULL;
	    }
	    deflateParams(&s, level[i], Z_DEFAULT_STRATEGY);
	    err = deflate(&s, Z_SYNC_FLUSH); // also try Z_FULL_FLUSH
	    cdata_pos = cdata_alloc - s.avail_out;
	    if (err != Z_OK) {
		fprintf(stderr, "zlib deflate error: %s\n", s.msg);
		break;
	    }
	}
    }
    if (deflate(&s, Z_FINISH) != Z_STREAM_END) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    *cdata_size = s.total_out;

    if (deflateEnd(&s) != Z_OK) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    return cdata;
}

static char *mem_inflate(char *cdata, size_t csize, size_t *size) {
    z_stream s;
    char *data = NULL; /* Uncompressed output */
    int data_alloc = 0;
    int err;

    /* Starting point at uncompressed size, 4x compressed */
    data = malloc(data_alloc = csize*4+10);

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)cdata;
    s.avail_in = csize;
    s.total_in = 0;
    s.next_out  = data;
    s.avail_out = data_alloc;
    s.total_out = 0;

    err = inflateInit(&s);

    /* Decode to 'data' array */
    for (;s.avail_in;) {
	s.next_out = (unsigned char *)&data[s.total_out];
	err = inflate(&s, Z_NO_FLUSH);
	if (err == Z_STREAM_END)
	    break;

	if (err != Z_OK) {
	    fprintf(stderr, "zlib inflate error: %s\n", s.msg);
	    break;
	}

	/* More to come, so realloc */
	data = realloc(data, data_alloc += s.avail_in*4 + 10);
	s.avail_out += s.avail_in*4+10;
    }
    inflateEnd(&s);

    *size = s.total_out;
    return data;
}


#if 0
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
    return ret ? -1 : 0;
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

    return ret ? -1 : 0;
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
    return ret ? -1 : 0;
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

    rdstats[GT_BTree] += len;
    rdcounts[GT_BTree]++;

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
    btree_query_t *bt = (btree_query_t *)clientdata;
    g_io *io = bt->io;
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
	wrstats[GT_BTree] += len;
	wrcounts[GT_BTree]++;
	ret = g_write(io, ci->view, data, len);
	g_flush(io, ci->view);
    } else {
	if (-1 == (v = lock(io, n->rec, G_LOCK_EX))) {
	    fprintf(stderr, "Failed to lock btree node %d\n", n->rec);
	    return -1;
	}
	wrstats[GT_BTree] += len;
	wrcounts[GT_BTree]++;
	ret = g_write(io, v, data, len);
	unlock(io, v);
    }

    free(data);

    if (ret) {
	fprintf(stderr, "Failed to write btree node %d\n", n->rec);
    }

    return ret ? -1 : 0;
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
	    cached_item *ci = hi->data.p;
	    btree_node_t *n = (btree_node_t *)ci->data;
	    assert(ci->updated == 0);
	    unlock(io, ci->view);
	    free(ci);
	    btree_del_node(n);
	}
    }

    if (h->clientdata)
	free(h->clientdata);

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
    heap_destroy(h, 1);

    /* LOW LEVEL IO HERE */
    if ( (fd = creat(auxfn,G_DEF_PERMS)) == -1 ) return gerr_set(GERR_CANT_CREATE);

    /* initialise header */
    auxheader.file_size = 0;
    auxheader.block_size = (int32_t)8; /* unused now? */
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
    if (!io->gdb)
	return NULL;
    io->client = g_connect_client_(io->gdb, 0, G_LOCK_EX, &io->mode);
    if (io->client == -1)
	return NULL;
    
    //g_lock_file_N_(io->gdb, io->client, 0);

#ifdef INDEX_NAMES
    io->seq_name_hash = HacheTableCreate(1024,
					 HASH_DYNAMIC_SIZE | HASH_OWN_KEYS);
    io->seq_name_hash->name = "io->seq_name_hash";

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

    io->contig_name_hash = HacheTableCreate(1024,
					    HASH_DYNAMIC_SIZE | HASH_OWN_KEYS);
    io->contig_name_hash->name = "io->contig_name_hash";

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

int io_database_lock(void *dbh) {
    g_io *io = (g_io *)dbh;
    g_lock_file_N_(io->gdb, io->client, 0);
    return 0;
}

int io_database_unlock(void *dbh) {
    g_io *io = (g_io *)dbh;
    g_unlock_file_N_(io->gdb, io->client, 0);
    return 0;
}

static int io_database_commit(void *dbh) {
    g_io *io = (g_io *)dbh;

    btree_flush(io, io->seq_name_hash);
    btree_flush(io, io->contig_name_hash);
    
    //g_unlock_file_N_(io->gdb, io->client, 0);
    //g_lock_file_N_(io->gdb, io->client, 0);

    return 0;
}

static int io_database_disconnect(void *dbh) {
    g_io *io = (g_io *)dbh;

    //btree_print(io->seq_name_tree, io->seq_name_tree->root, 0);

    io_database_commit(dbh);

    if (io->seq_name_hash) {
	btree_destroy(io, io->seq_name_hash);
	if (io->seq_name_tree) {
	    free(io->seq_name_tree->cd);
	    free(io->seq_name_tree);
	}
    }

    if (io->contig_name_hash) {
	btree_destroy(io, io->contig_name_hash);
	if (io->contig_name_tree) {
	    free(io->contig_name_tree->cd);
	    free(io->contig_name_tree);
	}
    }

    g_disconnect_client_(io->gdb, io->client);
    g_shutdown_database_(io->gdb);

    free(io);

    printf("\n*** I/O stats (type, read, write) ***\n");
    printf("GT_RecArray     \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_RecArray],     wrstats[GT_RecArray],
	   rdcounts[GT_RecArray],     rdstats[GT_RecArray]);
    printf("GT_Bin          \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_Bin],          wrstats[GT_Bin],
	   rdcounts[GT_Bin],          rdstats[GT_Bin]);
    printf("GT_Range        \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_Range],        wrstats[GT_Range],
	   rdcounts[GT_Range],        rdstats[GT_Range]);
    printf("GT_BTree        \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_BTree],        wrstats[GT_BTree],
	   rdcounts[GT_BTree],        rdstats[GT_BTree]);
    printf("GT_Track        \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_Track],        wrstats[GT_Track],
	   rdcounts[GT_Track],        rdstats[GT_Track]);
    printf("GT_Contig       \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_Contig],       wrstats[GT_Contig],
	   rdcounts[GT_Contig],       rdstats[GT_Contig]);
    printf("GT_Seq          \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_Seq],          wrstats[GT_Seq],
	   rdcounts[GT_Seq],          rdstats[GT_Seq]);
    printf("GT_Anno         \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_Anno],         wrstats[GT_Anno],
	   rdcounts[GT_Anno],         rdstats[GT_Anno]);
    printf("GT_AnnoEle      \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_AnnoEle],      wrstats[GT_AnnoEle],
	   rdcounts[GT_AnnoEle],      rdstats[GT_AnnoEle]);
    printf("GT_SeqBlock     \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_SeqBlock],     wrstats[GT_SeqBlock],
	   rdcounts[GT_SeqBlock],     rdstats[GT_SeqBlock]);
    printf("GT_AnnoEleBlock \t%7d\t%14d\t%7d\t%14d\n",
	   wrcounts[GT_AnnoEleBlock], wrstats[GT_AnnoEleBlock],
	   rdcounts[GT_AnnoEleBlock], rdstats[GT_AnnoEleBlock]);

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

    /* Libraries */
    db.Nlibraries = 0;
    db.library = allocate(io);
    v = lock(io, db.library, G_LOCK_EX);
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

    io_database_commit(io);
   
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
	return NULL;
    }

    rdstats[GT_Contig] += len;
    rdcounts[GT_Contig]++;

    /* Generate in-memory data structure */
    slen = sizeof(*c) + sizeof(char *) + ((unsigned char *)(ch+1))[0] + 1;
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
    ch->bin   = le_int4(c->bin);
    ((unsigned char *)ch)[sizeof(GContig_header)] =
	c->name ? strlen(c->name) : 0;
    if (c->name)
	strncpy(&((char *)ch)[sizeof(GContig_header)+1], c->name,
		((unsigned char *)ch)[sizeof(GContig_header)]);

    /* Write the data */
    wrstats[GT_Contig] += len;
    wrcounts[GT_Contig]++;
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

    rdstats[GT_RecArray] += vi.used;
    rdcounts[GT_RecArray]++;

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
    wrstats[GT_RecArray] += ArrayMax(ar) * sizeof(GCardinal);
    wrcounts[GT_RecArray]++;
    ret = g_write(io, ci->view, ArrayBase(GCardinal, ar), 
		   ArrayMax(ar) * sizeof(GCardinal));
    g_flush(io, ci->view);

    return ret ? -1 : 0;
}


/* ------------------------------------------------------------------------
 * anno element access methods
 */

/*
 * Reads an anno_ele_t struct and returns a cached_item containing this.
 * In memory we have the structure and the comment packed together into
 * one malloc call.
 *
 * On disc we have:
 * ? byte bin record
 * ? byte tag_type
 * ? byte obj_type
 * ? byte obj_rec
 * ? byte anno_rec
 * ? byte comment length (L)
 * L bytes of comment
 */
static cached_item *io_anno_ele_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    void *bloc;
    char *cp;
    size_t bloc_len;
    GView v;
    cached_item *ci;
    int anno_rec, tag_type, obj_type, obj_rec, comment_len, bin_rec;
    anno_ele_t *e;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    bloc = g_read_alloc(io, v, &bloc_len);

    rdstats[GT_AnnoEle] += bloc_len;
    rdcounts[GT_AnnoEle]++;

    if (!bloc)
	return NULL;


    /* Decode it */
    cp = bloc; 
    cp += u72int(cp, &bin_rec);
    cp += u72int(cp, &tag_type);
    cp += u72int(cp, &obj_type);
    cp += u72int(cp, &obj_rec);
    cp += u72int(cp, &anno_rec);
    cp += u72int(cp, &comment_len);
    ci = cache_new(GT_AnnoEle, rec, v, NULL, sizeof(anno_ele_t) +
		   comment_len + 1);
    e = (anno_ele_t *)&ci->data;
    e->bin      = bin_rec;
    e->tag_type = tag_type;
    e->obj_type = obj_type;
    e->obj_rec  = obj_rec;
    e->anno_rec = anno_rec;
    e->rec = rec;
    if (comment_len) {
	e->comment = (char *)&e->data;
	memcpy(e->comment, cp, comment_len);
	e->comment[comment_len] = 0;
    } else {
	e->comment = NULL;
    }
    

    /* And tidy up */
    free(bloc);
    ci->view = v;
    ci->rec = rec;
    return ci;

}

static int io_anno_ele_write_view(g_io *io, anno_ele_t *e, GView v) {
    int data_len, err = 0;
    unsigned char block[1024], *cp = block, *cpstart;
    int comment_len;

    /* Allocate memory if needed */
    comment_len = e->comment ? strlen(e->comment) : 0;
    data_len = 5 + 5 + 5 + 5 + 5 + comment_len;
    if (data_len > 1024) {
	if (NULL == (cp = (unsigned char *)malloc(data_len)))
	    return -1;
    }

    /* Encode */
    cpstart = cp;
    cp += int2u7(e->bin, cp);
    cp += int2u7(e->tag_type, cp);
    cp += int2u7(e->obj_type, cp);
    cp += int2u7(e->obj_rec, cp);
    cp += int2u7(e->anno_rec, cp);
    cp += int2u7(comment_len, cp);
    if (comment_len) {
	memcpy(cp, e->comment, comment_len);
	cp += comment_len;
    }

    /* Write */
    wrstats[GT_AnnoEle] += cp-cpstart;
    wrcounts[GT_AnnoEle]++;
    err |= g_write(io, v, (void *)cpstart, cp-cpstart);
    g_flush(io, v);

    if (cpstart != block)
	free(cpstart);

    return err ? -1 : 0;
}

static int io_anno_ele_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    anno_ele_t *e = (anno_ele_t *)&ci->data;

    assert(ci->lock_mode >= G_LOCK_RW);
    return io_anno_ele_write_view(io, e, ci->view);
}

static int io_anno_ele_create(void *dbh, void *vfrom) {
    g_io *io = (g_io *)dbh;
    anno_ele_t *from = (anno_ele_t *)vfrom;
    GRec rec;
    GView v;

    rec = allocate(io);
    v = lock(io, rec, G_LOCK_EX);
    
    if (from) {
	io_anno_ele_write_view(io, from, v);
    } else {
	anno_ele_t e;
	memset(&e, 0, sizeof(e));
	io_anno_ele_write_view(io, &e, v);
    }

    unlock(io, v);

    return rec;
}


/* ------------------------------------------------------------------------
 * anno access methods
 */
static cached_item *io_anno_read(void *dbh, GRec rec) {
    return io_generic_read(dbh, rec, 0);
}

static int io_anno_write(void *dbh, cached_item *ci) {
    return io_generic_write(dbh, ci);
}


/* ------------------------------------------------------------------------
 * dnasrc access methods
 */
static cached_item *io_library_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    library_t *lib;
    char *ch, *zpacked;
    size_t len, ssz;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    ch = g_read_alloc(io, v, &len);

    if (ch && len) {
	zpacked = mem_inflate(ch, len, &ssz);
	free(ch);
	len = ssz;
	ch = zpacked;
    }

    /* Generate in-memory data structure */
    if (!(ci = cache_new(GT_Library, rec, v, NULL, sizeof(*lib))))
	return NULL;

    lib = (library_t *)&ci->data;
    lib->rec = rec;
    if (ch == NULL || len == 0) {
	lib->insert_size[0] = 0;
	lib->insert_size[1] = 0;
	lib->insert_size[2] = 0;
	lib->sd[0] = 0;
	lib->sd[1] = 0;
	lib->sd[2] = 0;
	lib->machine = 0;
	lib->lib_type = 0;
	memset(lib->size_hist, 0, 3 * LIB_BINS * sizeof(lib->size_hist[0][0]));
    } else {
	int i, j, tmp;
	char *cp = ch;

	cp += u72int(cp, &lib->insert_size[0]);
	cp += u72int(cp, &lib->insert_size[1]);
	cp += u72int(cp, &lib->insert_size[2]);
	cp += u72int(cp, &tmp); lib->sd[0] = tmp/100.0;
	cp += u72int(cp, &tmp); lib->sd[1] = tmp/100.0;
	cp += u72int(cp, &tmp); lib->sd[2] = tmp/100.0;
	cp += u72int(cp, &lib->machine);
	cp += u72int(cp, &lib->lib_type);
	
	for (j = 0; j < 3; j++) {
	    int last = 0;
	    for (i = 0; i < LIB_BINS; i++) {
		cp += s72int(cp, &lib->size_hist[j][i]);
		lib->size_hist[j][i] += last;
		last = lib->size_hist[j][i];
	    }
	}
    }

    if (ch)
	free(ch);

    ci->view = v;
    ci->rec = rec;

    return ci;
}

static int io_library_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    library_t *lib = (library_t *)&ci->data;
    char cpstart[LIB_BINS*5*3+100], *cp = cpstart;
    int tmp, i, j, err;
    char *gzout;
    size_t ssz;

    assert(ci->lock_mode >= G_LOCK_RW);

    cp += int2u7(lib->insert_size[0], cp);
    cp += int2u7(lib->insert_size[1], cp);
    cp += int2u7(lib->insert_size[2], cp);
    tmp = lib->sd[0] * 100; cp += int2u7(tmp, cp);
    tmp = lib->sd[1] * 100; cp += int2u7(tmp, cp);
    tmp = lib->sd[2] * 100; cp += int2u7(tmp, cp);
    cp += int2u7(lib->machine, cp);
    cp += int2u7(lib->lib_type, cp);

    for (j = 0; j < 3; j++) {
	int last = 0;
	for (i = 0; i < LIB_BINS; i++) {
	    cp += int2s7(lib->size_hist[j][i] - last, cp);
	    last = lib->size_hist[j][i];
	}
    }

    /* Compress it */
    gzout = mem_deflate(cpstart, cp-cpstart, &ssz);
    //err = g_write(io, ci->view, cpstart, cp-cpstart);
    err = g_write(io, ci->view, gzout, ssz);
    free(gzout);
    g_flush(io, ci->view);

    if (0) {
	int i, j;

	puts("\nlibrary\n");
	for (j = 0; j < 3; j++) {
	    printf("  type %d\n", j);
	    for (i = 0; i < 1792; i++) {
		if (!lib->size_hist[j][i]) continue;
		printf("        %d\t%f\n", ibin2isize(i),
		       (double)lib->size_hist[j][i] / ibin_width(i));
	    }
	}
	printf("Rec %d view %d\n", ci->rec, ci->view);
    }

    return err;
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
static char *pack_rng_array(GRange *rng, int nr, int *sz) {
    int i;
    size_t part_sz[7];
    GRange last, last_tag;
    char *cp[6], *cp_orig[6], *out, *out_orig;
    //char *cpt, *cpt_orig;
    //HacheTable *h = HacheTableCreate(16, HASH_DYNAMIC_SIZE);
    //int ntags;

    memset(&last, 0, sizeof(last));
    memset(&last_tag, 0, sizeof(last_tag));

#if 0
    /* Allocate an index to all the tag types concerned */
    for (i = 0; i < nr; i++) {
	GRange r = rng[i];
	HacheData hd;

	if (!(r.flags & GRANGE_FLAG_ISANNO))
	    continue;

	HacheTableAdd(h, &r.mqual, sizeof(r.mqual), hd, NULL);
    }
    for (ntags = i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    hi->data.i = ntags++;
	}
    }

    /* First write out tag type dictionary */
    cpt = cpt_orig = malloc(5+ntags*4);
    cpt += int2u7(ntags, cpt);
    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    memcpy(cpt, hi->key, 4);
	    cpt += 4;
	}
    }
    /* FIXME: add cpt dict to 'out' later on */
#endif

    /* Pack the 6 structure elements to their own arrays */
    for (i = 0; i < 6; i++)
	cp[i] = cp_orig[i] = malloc(nr * 5);

    for (i = 0; i < nr; i++) {
	GRange r = rng[i];

	if (r.flags & GRANGE_FLAG_ISANNO) {
	    r.end   -= r.start;
	    r.start -= last_tag.start;
	    r.rec   -= last_tag.rec;
	} else {
	    r.end   -= r.start;
	    r.start -= last.start;
	    r.rec   -= last.rec;
	}

	cp[0] += int2u7(r.start, cp[0]);
	cp[1] += int2u7(r.end,   cp[1]);
	cp[2] += int2u7(r.rec,   cp[2]);
	cp[3] += int2u7(r.mqual, cp[3]);
	cp[4] += int2u7(r.flags, cp[4]);

	if (r.flags & GRANGE_FLAG_ISANNO) {
	    if (!(r.flags & GRANGE_FLAG_TYPE_SINGLE))
		cp[5] += int2s7(r.pair_rec - last_tag.pair_rec, cp[5]);
	    last_tag = rng[i];
	} else {
	    if (!(r.flags & GRANGE_FLAG_TYPE_SINGLE))
		cp[5] += int2s7(r.pair_rec - last.pair_rec, cp[5]);
	    last = rng[i];
	}
    }

    for (i = 0; i < 6; i++)
	part_sz[i+1] = cp[i]-cp_orig[i];

    /* Construct a header with nr and the size of the 6 packed struct fields */
    *sz = 7*5 + part_sz[1] + part_sz[2] + part_sz[3] +
	part_sz[4] + part_sz[5] + part_sz[6];

    out = out_orig = (char *)malloc(*sz);
    out += int2u7(nr, out);
    for (i = 0; i < 6; i++)
	out += int2u7(cp[i]-cp_orig[i], out);
    part_sz[0] = out-out_orig;

    /* Followed by the serialised 6 packed fields themselves */
    for (i = 0; i < 6; i++) {
	int len = cp[i]-cp_orig[i];
#if 0
	/* Run length encoding - doesn't gain us much so skip */
	int last_sym = 0, count = 0, j;

	for (j = 0; j < len; j++) {
	    *out++ = cp_orig[i][j];
	    if (cp_orig[i][j] == last_sym)
		count++;
	    else
		count = 0;

	    if (count >= 2) {
		int k = 0;
		while (cp_orig[i][++j] == last_sym && k < 255)
		    k++;
		*out++ = k;
		if (cp_orig[i][j] != last_sym)
		    j--;
	    }
	    last_sym = cp_orig[i][j];
	}
#else
	memcpy(out, cp_orig[i], len);
	out += len;
#endif
	free(cp_orig[i]);
    }

    *sz = out-out_orig;

    /* Gzip it too */
    {
    	char *gzout;
	size_t ssz;

	if (*sz < 512)
	    gzout = mem_deflate(out_orig, *sz, &ssz);
	else
	    gzout = mem_deflate_parts(out_orig, part_sz, 7, &ssz);
	*sz = ssz;

    	free(out_orig);
    	out_orig = gzout;
    }

    //write(2, out_orig, *sz);

    //HacheTableDestroy(h, 0);

    return out_orig;
}

static GRange *unpack_rng_array(char *packed, int packed_sz, int *nr) {
    int i, off[6];
    char *cp[6], *zpacked = NULL;
    GRange last, last_tag, *r, *ls = &last, *lt = &last;
    size_t ssz;

    /* First of all, inflate the compressed data */
    zpacked = packed = mem_inflate(packed, packed_sz, &ssz);
    packed_sz = ssz;

    /* Unpack number of ranges */
    packed += u72int(packed, nr);

    /* Unpack offsets of the 6 range components */
    for (i = 0; i < 6; i++) 
	packed += u72int(packed, &off[i]);
    cp[0] = packed;
    for (i = 1; i < 6; i++)
	cp[i] = cp[i-1] + off[i-1];

    r = (GRange *)malloc(*nr * sizeof(*r));
    memset(ls, 0, sizeof(*ls));
    memset(lt, 0, sizeof(*lt));

    /* And finally unpack from the 6 components in parallel for each struct */
    for (i = 0; i < *nr; i++) {
	cp[0] += u72int(cp[0], &r[i].start);
	cp[1] += u72int(cp[1], &r[i].end);
	cp[2] += u72int(cp[2], &r[i].rec);
	cp[3] += u72int(cp[3], &r[i].mqual);
	cp[4] += u72int(cp[4], &r[i].flags);
	if (r[i].flags & GRANGE_FLAG_ISANNO) {
	    if (!(r[i].flags & GRANGE_FLAG_TYPE_SINGLE)) {
		int32_t pr;
		cp[5] += s72int(cp[5], &pr);
		r[i].pair_rec = pr + lt->pair_rec;
	    } else {
		r[i].pair_rec = 0;
	    }
	    
	    r[i].rec += lt->rec;
	    r[i].start += lt->start;
	    r[i].end += r[i].start;

	    lt = &r[i];
	} else {
	    if (!(r[i].flags & GRANGE_FLAG_TYPE_SINGLE)) {
		int32_t pr;
		cp[5] += s72int(cp[5], &pr);
		r[i].pair_rec = pr + ls->pair_rec;
	    } else {
		r[i].pair_rec = 0;
	    }
	    
	    r[i].rec += ls->rec;
	    r[i].start += ls->start;
	    r[i].end += r[i].start;

	    ls = &r[i];
	}
    }

    assert(cp[5] - zpacked == packed_sz);

    if (zpacked)
	free(zpacked);

    return r;
}

/* Used in the on-disc encoding */
#define BIN_COMPLEMENTED  (1<<0)
#define BIN_NO_RANGE      (1<<1)
#define BIN_NO_TRACK      (1<<2)
#define BIN_NO_LCHILD     (1<<3)
#define BIN_NO_RCHILD     (1<<4)
#define BIN_SIZE_EQ_POS   (1<<5)
#define BIN_POS_ZERO      (1<<6)
#define BIN_ROOT_NODE     (1<<7)

static cached_item *io_bin_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    GBin *b, g;
    bin_index_t *bin;
    char *buf, *cp;
    size_t buf_len;
    int bflag;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    if (NULL == (buf = g_read_alloc(io, v, &buf_len)))
	return NULL;
    if (buf_len < sizeof(GBin)) {
	buf = realloc(buf, sizeof(GBin));
	memset(buf + buf_len, 0, sizeof(GBin) - buf_len);
    }
    //b = (GBin *)buf;
    b = &g;
    cp = buf;

    rdstats[GT_Bin] += buf_len;
    rdcounts[GT_Bin]++;

    cp += u72int(cp, &bflag);
    g.flags = (bflag & BIN_COMPLEMENTED) ? BIN_COMPLEMENTED : 0;
    g.parent_type = (bflag & BIN_ROOT_NODE) ? GT_Contig : GT_Bin;

    if (bflag & BIN_POS_ZERO)
	g.pos = 0;
    else
	cp += s72int(cp, &g.pos);

    if (bflag & BIN_SIZE_EQ_POS)
	g.size = g.pos;
    else
	cp += u72int(cp, &g.size);

    if (bflag & BIN_NO_RANGE) {
	g.range = 0;
	g.start = 0;
	g.end   = 0;
    } else {
	cp += u72int(cp, &g.start);
	cp += u72int(cp, &g.end);
	g.end += g.start;
	cp += u72int(cp, &g.range);
    }

    if (bflag & BIN_NO_LCHILD)
	g.child[0] = 0;
    else
	cp += u72int(cp, &g.child[0]);

    if (bflag & BIN_NO_RCHILD)
	g.child[1] = 0;
    else
	cp += u72int(cp, &g.child[1]);

    if (bflag & BIN_NO_TRACK)
	g.track = 0;
    else
	cp += u72int(cp, &g.track);

    cp += u72int(cp, &g.parent);
    cp += u72int(cp, &g.nseqs);

    /*
    printf("<%d / p=%d+%d, %d..%d p=%d/%d, ch=%d/%d, id=%d, f=%d t=%d ns=%d r=%d\n",
	   rec,
	   g.pos, g.size, g.start, g.end,
	   g.parent_type, g.parent, g.child[0], g.child[1],
	   g.id, g.flags, g.track, g.nseqs, g.range);
    */
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
    bin->nseqs       = b->nseqs;

    /* Load ranges */
    if (b->range) {
	GViewInfo vi;
	int nranges;
	GRange *r;
	char *buf;

	v = lock(io, b->range, G_LOCK_RO);
	g_view_info_(io->gdb, io->client, v, &vi);
	
	buf = (char *)malloc(vi.used);
	g_read(io, v, buf, vi.used);
	r = unpack_rng_array(buf, vi.used, &nranges);
	free(buf);

	rdstats[GT_Range] += vi.used;
	rdcounts[GT_Range]++;

	//printf("Unpacked %d ranges from %d bytes\n", nranges, vi.used);

	bin->rng = ArrayCreate(sizeof(GRange), nranges);
	if (ArrayBase(GRange, bin->rng))
	    free(ArrayBase(GRange, bin->rng));
	bin->rng->base = r;
	ArrayRef(bin->rng, nranges-1);

	//	g_read(io, v, ArrayBase(GRange, bin->rng),
	//	       nranges * sizeof(GRange));
	unlock(io, v);
    }

    /* Load tracks */
    if (b->track) {
	GViewInfo vi;
	int ntracks;

	v = lock(io, b->track, G_LOCK_RO);
	g_view_info_(io->gdb, io->client, v, &vi);
	ntracks = vi.used / sizeof(GBinTrack);

	rdstats[GT_Track] += vi.used;
	rdcounts[GT_Track]++;

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
    unsigned int bflag;

    /* Ranges */
    if (bin->flags & BIN_RANGE_UPDATED) {
	GView v;
	char *cp;
	int sz;

	bin->flags &= ~BIN_RANGE_UPDATED;

	if (!bin->rng_rec) {
	    bin->rng_rec = allocate(io);
	    bin->flags |= BIN_BIN_UPDATED;
	}

	cp = pack_rng_array(ArrayBase(GRange, bin->rng), ArrayMax(bin->rng),
			    /* bin->start_used, */ &sz);
	//printf("Packed %d ranges in %d bytes\n", ArrayMax(bin->rng), sz);

#ifdef DEBUG
	{
	    char fn[1024];
	    sprintf(fn, "/tmp/jkb/rng.%d", bin->rec);
	    int fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC, 0666);
	    if (fd == -1) {
		perror (fn);
	    } else {
		write(fd, cp, sz);
		close(fd);
	    }
	}
#endif

	v = lock(io, bin->rng_rec, G_LOCK_EX);
	//	err |= g_write(io, v, ArrayBase(GRange, bin->rng),
	//	       sizeof(GRange) * ArrayMax(bin->rng));
	wrstats[GT_Range] += sz;
	wrcounts[GT_Range]++;
	err |= g_write(io, v, cp, sz);
	free(cp);
	err |= unlock(io, v);
    }

    /* Tracks */
    if (bin->flags & BIN_TRACK_UPDATED) {
	GView v;

	bin->flags &= ~BIN_TRACK_UPDATED;

	if (bin->track) {
	    if (!bin->track_rec) {
		bin->track_rec = allocate(io);
		bin->flags |= BIN_BIN_UPDATED;
	    }

	    v = lock(io, bin->track_rec, G_LOCK_EX);
	    wrstats[GT_Track] += sizeof(GBinTrack) * ArrayMax(bin->track);
	    wrcounts[GT_Track]++;
	    err |= g_write(io, v, ArrayBase(GBinTrack, bin->track),
			   sizeof(GBinTrack) * ArrayMax(bin->track));
	    err |= unlock(io, v);
	}
    }

    /* Bin struct itself */
    if (bin->flags & BIN_BIN_UPDATED) {
	char cpstart[12*5], *cp = cpstart;

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
	g.nseqs       = bin->nseqs;

#ifdef DEBUG
	{
	    char fn[1024];
	    sprintf(fn, "/tmp/jkb/bin.%d", bin->rec);
	    int fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC, 0666);
	    if (fd == -1) {
		perror (fn);
	    } else {
		write(fd, &g, sizeof(g));
		close(fd);
	    }
	}
#endif
	/*
	printf(">%d / p=%d+%d, %d..%d p=%d/%d, ch=%d/%d, f=%d t=%d ns=%d r=%d\n",
	       bin->rec,
	       g.pos, g.size, g.start, g.end,
	       g.parent_type, g.parent, g.child[0], g.child[1],
	       g.flags, g.track, g.nseqs, g.range);
	*/

	/* Encode the bin in a more efficient struct */
	bflag = 0;
	if (g.flags & BIN_COMPLEMENTED) bflag |= BIN_COMPLEMENTED;
	if (!g.child[0])                bflag |= BIN_NO_LCHILD;
	if (!g.child[1])                bflag |= BIN_NO_RCHILD;
	if (!g.range)                   bflag |= BIN_NO_RANGE;
	if (!g.track)                   bflag |= BIN_NO_TRACK;
	assert(g.parent_type == GT_Bin || g.parent_type == GT_Contig);
	if (g.parent_type == GT_Contig) bflag |= BIN_ROOT_NODE;
	if (g.pos == 0)                 bflag |= BIN_POS_ZERO;
	if (g.size == g.pos)            bflag |= BIN_SIZE_EQ_POS;

	cp += int2u7(bflag, cp);
	if (!(bflag & BIN_POS_ZERO))     cp += int2s7(g.pos, cp);
	if (!(bflag & BIN_SIZE_EQ_POS))  cp += int2u7(g.size, cp);
	if (!(bflag & BIN_NO_RANGE)) {
	    cp += int2u7(g.start, cp);
	    cp += int2u7(g.end - g.start, cp);
	    cp += int2u7(g.range, cp);
	}
	if (!(bflag & BIN_NO_LCHILD))    cp += int2u7(g.child[0], cp);
	if (!(bflag & BIN_NO_RCHILD))    cp += int2u7(g.child[1], cp);
	if (!(bflag & BIN_NO_TRACK))     cp += int2u7(g.track, cp);
	cp += int2u7(g.parent, cp);
	cp += int2u7(g.nseqs, cp);

	wrstats[GT_Bin] += cp-cpstart;
	//wrstats[GT_Bin] += sizeof(g);
	wrcounts[GT_Bin]++;
	err |= g_write(io, v, cpstart, cp - cpstart);
	//err |= g_write(io, v, &g, sizeof(g));
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
    //bin_index_t *from = vfrom;
    g_io *io = (g_io *)dbh;
    GRec rec;

    rec = allocate(io);
    return rec;

    /*
     * Delay this now until after we've fully initialised it.
     * This was previously causing a double write as we want the ability
     * to allocate a record before we actually write data to it.
     *
     * See changes to bin_new() in tg_bin.c
     */
#if 0
    GView v;
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
	b.nseqs	      = 0;
	io_bin_write_view(io, &b, v);
    }
    unlock(io, v);

    return rec;
#endif
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

    rdstats[GT_Track] += buf_len;
    rdcounts[GT_Track]++;

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
    
    wrstats[GT_Track] += sizeof(*h) + track->item_size * track->nitems;
    wrcounts[GT_Track]++;
    err |= g_write(io, v, h, sizeof(*h) + track->item_size * track->nitems);
    g_flush(io, v);
    free(data);

    return err;
}

static int io_track_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    track_t *track = (track_t *)&ci->data;

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
 * ? byte bin record no.
 * ? byte bin index
 * ? byte 'left clip'
 * ? byte 'right clip'
 * ? byte sequence length
 * ? byte parent_rec;
 * 1 byte parent_type;
 * 1 byte seq_tech (3 bottom bits)
 *      + flags (3 next bits)
 *      + format (2 top bits)
 * 1 byte mapping_quality
 * ? byte Nanno (num. annotations)
 * ? byte Annotation record numbers (Nanno copies of ? byte) 
 * ? bytes name, nul terminated
 * ? byte trace name, nul terminated
 * ? bytes alignment, nul terminated
 * remainder is seq/qual (various formats).
 *
 * Returns a pointer to a seq_t struct
 *      or NULL on failure.
 */
static cached_item *seq_decode(unsigned char *buf, size_t len, int rec) {
    cached_item *ci;
    unsigned char *cp, flags, mapping_qual, seq_tech, format;
    size_t slen;
    signed int i, j;
    seq_t *seq;
    uint32_t left, right, bin, seq_len;
    int parent_rec, parent_type, bin_index;
    Array anno;

    if (len) {
	int Nanno;
	cp = buf;
	cp += u72int(cp, &bin);
	cp += u72int(cp, &bin_index);
	cp += u72int(cp, &left);
	cp += u72int(cp, &right);
	cp += u72int(cp, &seq_len);
	cp += u72int(cp, &parent_rec);
	parent_type = *cp++;
	format = *cp++;
	seq_tech = format & ((1<<3)-1);
	format >>= 3;
	flags = format & ((1<<3)-1);
	format >>= 3;
	mapping_qual = *cp++;
	cp += u72int(cp, &Nanno);
	if (Nanno) {
	    anno = ArrayCreate(sizeof(int), Nanno);
	    for (i = 0; i < Nanno; i++) {
		cp += u72int(cp, arrp(int, anno, i));
	    }
	} else {
	    anno = NULL;
	}
	/* cp is now the variable sized section starting with reading name */
    } else {
	/* new sequence */
	bin = 0;
	bin_index = 0;
	left = 0;
	right = 0;
	seq_len = 0;
	parent_rec = 0;
	parent_type = 0;
	seq_tech = 0;
	flags = 0;
	mapping_qual = 0;
	format = 0;
	anno = NULL;
	cp = buf = "\0\0\0"; /* seq name, trace name, alignment */
	len = 3;
    }

    /* Generate in-memory data structure */
    slen = sizeof(seq_t) + len - (cp-buf) +
	seq_len * (1 + (format == SEQ_FORMAT_CNF4 ? 4 : 1));

    if (!(ci = cache_new(GT_Seq, 0, 0, NULL, slen)))
        return NULL;
    seq = (seq_t *)&ci->data;

    seq->rec          = rec;
    seq->bin          = bin;
    seq->bin_index    = bin_index;
    seq->left         = left;
    seq->right        = right;
    seq->parent_type  = parent_type;
    seq->parent_rec   = parent_rec;
    seq->seq_tech     = seq_tech;
    seq->flags        = flags;
    seq->format       = format;
    seq->mapping_qual = mapping_qual;
    seq->anno         = anno;
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

/* basic MTF code */
static unsigned char order[256];
static void reorder_init(void) {
    int i;
    for (i = 0; i < 256; i++)
	order[i] = i;
}
static int reorder(int c) {
    int i;
    for (i = 0; i < 256; i++) {
	if (order[i] == c)
	    break;
    }
    if (i == 256)
	abort();

    memmove(&order[1], &order[0], i);
    order[0] = c;

    return i;
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

    rdstats[GT_Seq] += bloc_len;
    rdcounts[GT_Seq]++;

    if (!bloc)
	return NULL;

    ci = seq_decode(bloc, bloc_len, rec);
    free(bloc);

    printf("Read rec %d => %p\n", rec, ci);

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

    seq->format = SEQ_FORMAT_CNF1;

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
	+ 5 /* bin index */
	+ 5 /* left clip */
	+ 5 /* right clip */
	+ 5 /* seq len */
	+ 5 /* parent_rec */
	+ 1 /* parent_type */
	+ 1 /* seq_tech/flags/format */
	+ 1 /* mapping_quality */
	+ 5 /* Nanno */
	+ 5*(seq->anno ? ArrayMax(seq->anno) : 0) /* anno rec.nos */
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
    cp += int2u7(seq->bin_index, cp);
    cp += int2u7(seq->left, cp);
    cp += int2u7(seq->right, cp);
    cp += int2u7(seq_len, cp);

    /* Read-pair info */
    cp += int2u7(seq->parent_rec, cp);
    *cp++ = seq->parent_type;

    /* flags & m.quality */
    *cp++ = (seq->format << 6) | (seq->flags << 3) | seq->seq_tech;
    *cp++ = seq->mapping_qual;

    /* Annotations */
    if (seq->anno) {
	cp += int2u7(ArrayMax(seq->anno), cp);
	for (i = 0; i < ArrayMax(seq->anno); i++) {
	    cp += int2u7(arr(int, seq->anno, i), cp);
	}
    } else {
	cp += int2u7(0, cp);
    }

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
	for (i = 0; i < seq_len; i++) {
	    *cp++ = seq->conf[i];
	}
#else
#if 1
    {
	char tmp[1024], *t = tmp;
	static huffman_codeset_t *cs = NULL;
	block_t *blk;

	if (!cs) {
	    huffman_code_t codes[] = {
	{  0,  2}, {  1,  2}, {  2,  3}, {  3,  3}, {  4,  6}, { 28,  6},
	{ 29,  6}, { 30,  6}, { 31,  6}, { 32,  6}, {'T',  6}, {  5,  7},
	{  6,  7}, { 25,  7}, { 26,  7}, { 27,  7}, { 33,  7}, { 34,  7},
	{'C',  7}, {'D',  7}, {'G',  7}, {'H',  7}, {  7,  8}, {  8,  8},
	{ 21,  8}, { 22,  8}, { 23,  8}, { 24,  8}, {'A',  8}, {'B',  8},
	{'E',  8}, {  9,  9}, { 10,  9}, { 11,  9}, { 16,  9}, { 17,  9},
	{ 18,  9}, { 19,  9}, { 20,  9}, { 12, 10}, { 13, 10}, { 15, 11},
	{'m', 12}, {'q', 13}, {'z', 13}, { 14, 15}, { 35, 15}, { 37, 15},
	{ 43, 15}, { 45, 15}, { 46, 15}, {'0', 15}, {'1', 15}, {'2', 15},
	{'3', 15}, {'4', 15}, {'5', 15}, {'6', 15}, {'8', 15}, {'L', 15},
	{'P', 15}, {'a', 15}, {'c', 15}, {'d', 15}, {'e', 15}, {'f', 15},
	{'g', 15}, {'h', 15}, {'i', 15}, {'n', 15}, {'o', 15}, {'p', 15},
	{'r', 15}, {'s', 15}, {'t', 15}, {'u', 15}, {256 /* EOF */, 15},
	    };

	    cs = codes2codeset(codes, sizeof(codes)/sizeof(huffman_code_t), 99);
	}

	reorder_init();
	for (i = 0; i < seq_len; i++) {
	    *t++ = reorder(seq->seq[i]);
	}
	reorder_init();
	for (i = 0; i < seq_len; i++) {
	    *t++ = reorder(seq->conf[i]);
	}

	blk = block_create(NULL, 8192);
	huffman_multi_encode(blk, cs, 99, tmp, seq_len*2);

	//	write(1, blk->data, blk->byte + (blk->bit != 0));
	//	write(2, tmp, seq_len*2);
	memcpy(cp, blk->data, blk->byte + (blk->bit != 0));
	cp += blk->byte + (blk->bit != 0);
	block_destroy(blk, 0);
    }
#else
	for (i = 0; i < seq_len; i++) {
	    *cp++ = seq->seq[i];
	}
	for (i = 0; i < seq_len; i++) {
	    *cp++ = seq->conf[i];
	}
#endif
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
    wrstats[GT_Seq] += cp-cpstart;
    wrcounts[GT_Seq]++;
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

#if 0
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
	s.bin_index = 0;
	s.pos = 0;
	s.len = 0;
	s.left = s.right = 0;
	s.seq_tech = 0;
	s.flags = 0;
	s.parent_rec = 0;
	s.parent_type = 0;
	s.mapping_qual = 0;
	s.anno = NULL;
	s.name_len = 0;
	s.name = NULL;
	s.seq = NULL;
	s.conf = NULL;
	io_seq_write_view(io, &s, v, rec);
    }
    unlock(io, v);

    return rec;
}
#else
static int io_seq_create(void *dbh, void *vfrom) {
    g_io *io = (g_io *)dbh;
    return allocate(io);
}
#endif

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
 * seq_block access methods
 */
static cached_item *io_seq_block_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    GView v;
    cached_item *ci;
    seq_block_t *b;
    unsigned char *buf, *cp;
    size_t buf_len;
    seq_t in[SEQ_BLOCK_SZ];
    int i, last;

    set_dna_lookup();

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    if (!(ci = cache_new(GT_SeqBlock, rec, v, NULL, sizeof(*b))))
	return NULL;

    b = (seq_block_t *)&ci->data;
    cp = buf = (unsigned char *)g_read_alloc((g_io *)dbh, v, &buf_len);

    rdstats[GT_SeqBlock] += buf_len;
    rdcounts[GT_SeqBlock]++;

    /* Ungzip it too */
    if (1) {
	size_t ssz;
	buf = mem_inflate(buf, buf_len, &ssz);
	free(cp);
	cp = buf;
	buf_len = ssz;
    }
    b->est_size = buf_len;

    if (!buf_len) {
	memset(&b->rec[0], 0, SEQ_BLOCK_SZ*sizeof(b->rec[0]));
	memset(&b->seq[0], 0, SEQ_BLOCK_SZ*sizeof(b->seq[0]));
	free(buf);
	return ci;
    }

    /* Decode the fixed size components of our sequence structs */
    /* Bin */
    for (i = 0; i < SEQ_BLOCK_SZ; i++)
	cp += u72int(cp, &in[i].bin);

    /* Bin index */
    for (last = i = 0; i < SEQ_BLOCK_SZ; i++) {
	int32_t bi;
	if (!in[i].bin) continue;
	cp += s72int(cp, &bi);
	in[i].bin_index = last + bi;
	last = in[i].bin_index;
    }

    /* left clip */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, &in[i].left);
    }

    /* right clip */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, &in[i].right);
    }

    /* length */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, &in[i].len);
    }

    /* parent rec */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, &in[i].parent_rec);
    }

    /* parent type */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	in[i].parent_type = *cp++;
    }

    /* flags */
    for (i = 0 ; i < SEQ_BLOCK_SZ; i++) { 
	if (!in[i].bin) continue;
	unsigned char f = *cp++;
	in[i].seq_tech = f & 7;
	in[i].flags = (f >> 3) & 7;
	in[i].format = (f >> 6) & 3;
	if (in[i].flags & SEQ_COMPLEMENTED)
	    in[i].len = -in[i].len;
    }

    /* mapping quality */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	in[i].mapping_qual = *cp++;
    }

    /* name length */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, &in[i].name_len);
    }

    /* trace name length */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, &in[i].trace_name_len);
    }

    /* alignment length */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, &in[i].alignment_len);
    }


    /* Convert our static structs to cached_items */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (in[i].bin) {
	    cached_item *si;
	    size_t extra_len;

	    extra_len =
		sizeof(seq_t) + 
		in[i].name_len +
		in[i].trace_name_len +
		in[i].alignment_len + 
		ABS(in[i].len) +
		ABS(in[i].len) * (in[i].format == SEQ_FORMAT_CNF4 ? 4 : 1);
	    if (!(si = cache_new(GT_Seq, 0, 0, NULL, extra_len)))
		return NULL;

	    b->rec[i] = (rec << SEQ_BLOCK_BITS) + i;
	    b->seq[i] = (seq_t *)&si->data;
	    in[i].rec = (rec << SEQ_BLOCK_BITS) + i;
	    in[i].anno = NULL;
	    *b->seq[i] = in[i];
	    b->seq[i]->block = b;
	    b->seq[i]->idx = i;
	} else {
	    b->rec[i] = 0;
	    b->seq[i] = NULL;
	}
    }

    /* Decode variable sized components */
    /* Names */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!b->seq[i]) continue;
	b->seq[i]->name = (char *)&b->seq[i]->data;
	memcpy(b->seq[i]->name, cp, b->seq[i]->name_len);
	cp += b->seq[i]->name_len;
	b->seq[i]->name[b->seq[i]->name_len] = 0;
    }

    /* Trace names, delta from seq name */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!b->seq[i]) continue;
	b->seq[i]->trace_name = b->seq[i]->name + b->seq[i]->name_len + 1;
	if (b->seq[i]->trace_name_len) {
	    int tlen = *cp++;
	    memcpy(b->seq[i]->trace_name, b->seq[i]->name, tlen);
	    memcpy(&b->seq[i]->trace_name[tlen], cp,
		   b->seq[i]->trace_name_len - tlen);
	    cp += b->seq[i]->trace_name_len - tlen;
	}
	b->seq[i]->trace_name[b->seq[i]->trace_name_len] = 0;
    }

    /* Alignment strings (unused at present) */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!b->seq[i]) continue;
	b->seq[i]->alignment = b->seq[i]->trace_name + b->seq[i]->trace_name_len + 1;
	memcpy(b->seq[i]->alignment, cp, b->seq[i]->alignment_len);
	cp += b->seq[i]->alignment_len;
	b->seq[i]->alignment[b->seq[i]->alignment_len] = 0;
    }

    /* Sequence */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!b->seq[i]) continue;
	b->seq[i]->seq = b->seq[i]->alignment + b->seq[i]->alignment_len + 1;
	memcpy(b->seq[i]->seq, cp, ABS(b->seq[i]->len));
	if (b->seq[i]->len < 0)
	    complement_seq(b->seq[i]->seq, -b->seq[i]->len);
	cp += ABS(b->seq[i]->len);
    }

    /* Quality */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!b->seq[i]) continue;
	b->seq[i]->conf = b->seq[i]->seq + ABS(b->seq[i]->len);
	memcpy(b->seq[i]->conf, cp, ABS(b->seq[i]->len));
	cp += ABS(b->seq[i]->len);
    }

    assert(cp - buf == buf_len);
    free(buf);

    return ci;
}

static int io_seq_block_write(void *dbh, cached_item *ci) {
    int err;
    g_io *io = (g_io *)dbh;
    seq_block_t *b = (seq_block_t *)&ci->data;
    int i, last_index;
    char *cp, *cp_start;
    char *out[17], *out_start[17];
    size_t out_size[17], total_size;
    int level[17];

    set_dna_lookup();

    /* Compute worst-case sizes, for memory allocation */
    for (i = 0; i < 17; i++) {
	out_size[i] = 0;
    }
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	seq_t *s = b->seq[i];
	if (!s) {
	    out_size[0]++;
	    continue;
	}

	out_size[0] += 5; /* bin */
	out_size[1] += 5; /* bin_index */
	out_size[2] += 5; /* left */
	out_size[3] += 5; /* right */
	out_size[4] += 5; /* len */
	out_size[5] += 5; /* parent_rec */
	out_size[6] ++;   /* parent type */
	out_size[7] ++;   /* format */
	out_size[8] ++;   /* m.qual */
	out_size[9] += 5; /* name len */
	out_size[10]+= 5; /* trace name len */
	out_size[11]+= 5; /* alignment len */
	out_size[12]+= s->name_len+1;
	out_size[13]+= s->trace_name_len;
	out_size[14]+= s->alignment_len;
	out_size[15]+= ABS(s->len); /* seq */
	out_size[16]+= ABS(s->len) * (s->format == SEQ_FORMAT_CNF4 ? 4 : 1);
    }
    for (i = 0; i < 17; i++)
	out_start[i] = out[i] = malloc(out_size[i]+1);


    /* serialised sequences, separated by type of data */
    last_index = 0;
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	int delta;
	seq_t *s = b->seq[i];

	if (!s) {
	    *out[0]++= 0; /* signifies sequence not present */
	    continue;
	}

	out[0] += int2u7(s->bin, out[0]);
	delta = s->bin_index - last_index;
	last_index = s->bin_index;
	out[1] += int2s7(delta, out[1]);
	out[2] += int2u7(s->left, out[2]);
	out[3] += int2u7(s->right, out[3]);
	out[4] += int2u7(ABS(s->len), out[4]);

	out[5] += int2u7(s->parent_rec, out[5]);
	*out[6]++ = s->parent_type;
	
	/* flags & m.quality */
	s->format = SEQ_FORMAT_CNF1;
	*out[7]++ = (s->format << 6) | (s->flags << 3) | s->seq_tech;
	*out[8]++ = s->mapping_qual;
	
	/* Name */
	out[9] += int2u7(s->name_len, out[9]);
	memcpy(out[12], s->name, s->name_len);
	out[12] += s->name_len;

	/* Trace name */
	if (s->trace_name_len == 0 ||
	    (s->trace_name_len == s->name_len &&
	     0 == memcpy(s->trace_name, s->name, s->name_len))) {
	    /* Trace name and seq name are identical, so just skip it */
	    *out[10]++ = 0;
	} else {
	    out[10] += int2u7(s->trace_name_len, out[10]);

	    if (s->trace_name_len) {
		/* Delta to name */
		int j;
		for (j = 0; j < s->name_len && j < s->trace_name_len; j++)
		    if (s->trace_name[j] != s->name[j])
			break;
		*out[13]++ = j;
		memcpy(out[13], &s->trace_name[j], s->trace_name_len-j);
		out[13] += s->trace_name_len - j;
	    }
	}

	/* Alignment */
	s->alignment_len = 0; /* this is unused for now */
	out[11] += int2u7(s->alignment_len, out[11]);
	memcpy(out[14], s->alignment, s->alignment_len);
	out[14] += s->alignment_len;

	/* Sequences - store in alignment orientation for better compression */
	if (s->len < 0) {
	    complement_seq(s->seq, ABS(s->len));
	    memcpy(out[15], s->seq,  ABS(s->len)); out[15] += ABS(s->len);
	    complement_seq(s->seq, ABS(s->len));
	} else {
	    memcpy(out[15], s->seq,  ABS(s->len)); out[15] += ABS(s->len);
	}

	/* Quality */
	/* Reversing this too would be consistent with the sequence and it
	 * saves a tiny amount of space, but it's only about 0.3% and costs
	 * cpu
	 */
	memcpy(out[16], s->conf, ABS(s->len)); out[16] += ABS(s->len);
    }


    /* Concatenate data types together and adjust out_size to actual usage */
    for (total_size = i = 0; i < 17; i++) {
	out_size[i] = out[i] - out_start[i];
	total_size += out_size[i];
    }
    cp = cp_start = (char *)malloc(total_size+1);
    for (i = 0; i < 17; i++) {
	memcpy(cp, out_start[i], out_size[i]);
	cp += out_size[i];
	free(out_start[i]);
	level[i] = 5; /* default gzip level */
    }
    assert(cp - cp_start == total_size);

    level[12] = 7; /* name, 8 is approx 1% better, but 10% slower */
    level[13] = 7; /* trace name */
    level[16] = 6; /* conf */

    /* NB: should name and trace name be compressed together? They're
     * probably the same or highly related?
     */

    //printf("Block %d, est.size %d, actual size %d, diff=%d %f gzip=",
    //       ci->rec, b->est_size, cp-cp_start,
    //       cp-cp_start - b->est_size, (cp-cp_start + 0.0) / b->est_size);
    
    /* Gzip it too */
    if (1) {
	char *gzout;
	size_t ssz;

	//gzout = mem_deflate(cp_start, cp-cp_start, &ssz);
	gzout = mem_deflate_lparts(cp_start, out_size, level, 17, &ssz);
	free(cp_start);
	cp_start = gzout;
	cp = cp_start + ssz;
    }
    //printf("%d\n", cp-cp_start);

    /* Finally write the serialised data block */
    assert(ci->lock_mode >= G_LOCK_RW);
    wrstats[GT_SeqBlock] += cp-cp_start;
    wrcounts[GT_SeqBlock]++;
    err = g_write(io, ci->view, cp_start, cp - cp_start);
    g_flush(io, ci->view);
   
    free(cp_start);

    return err ? -1 : 0;
}

static GRec io_seq_block_create(void *dbh, void *vfrom) {
    g_io *io = (g_io *)dbh;
    GRec rec;
    GView v;

    rec = allocate(io);
    v = lock(io, rec, G_LOCK_EX);

    /* Write blank data here? */

    unlock(io, v);
    
    return rec;
}


/* ------------------------------------------------------------------------
 * seq_block access methods
 */
static cached_item *io_anno_ele_block_read(void *dbh, GRec rec) {
    g_io *io = (g_io *)dbh;
    GView v;
    cached_item *ci;
    anno_ele_block_t *b;
    unsigned char *buf, *cp;
    size_t buf_len;
    anno_ele_t in[ANNO_ELE_BLOCK_SZ];
    int i, last;
    int comment_len[ANNO_ELE_BLOCK_SZ];

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    if (!(ci = cache_new(GT_AnnoEleBlock, rec, v, NULL, sizeof(*b))))
	return NULL;

    b = (anno_ele_block_t *)&ci->data;
    cp = buf = (unsigned char *)g_read_alloc((g_io *)dbh, v, &buf_len);

    rdstats[GT_AnnoEleBlock] += buf_len;
    rdcounts[GT_AnnoEleBlock]++;

    /* Ungzip it too */
    if (1) {
	size_t ssz;
	buf = mem_inflate(buf, buf_len, &ssz);
	free(cp);
	cp = buf;
	buf_len = ssz;
    }
    b->est_size = buf_len;

    if (!buf_len) {
	memset(&b->rec[0], 0, ANNO_ELE_BLOCK_SZ*sizeof(b->rec[0]));
	memset(&b->ae[0],  0, ANNO_ELE_BLOCK_SZ*sizeof(b->ae[0]));
	free(buf);
	return ci;
    }

    /* Decode the fixed size components of our sequence structs */
    /* Bin */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++)
	cp += u72int(cp, &in[i].bin);

    /* Tag type */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, &in[i].tag_type);
    }

    /* Obj type */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, &in[i].obj_type);
    }
    
    /* Obj record */
    for (last = i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	int32_t tmp;
	if (!in[i].bin) continue;
	cp += s72int(cp, &tmp);
	in[i].obj_rec = last + tmp;
	tmp = in[i].obj_rec;
	//cp += u72int(cp, &in[i].obj_rec);
    }

    /* Anno record */
    for (last = i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	int32_t tmp;
	if (!in[i].bin) continue;
	cp += s72int(cp, &tmp);
	in[i].anno_rec = last + tmp;
	tmp = in[i].anno_rec;
	//cp += u72int(cp, &in[i].anno_rec);
    }

    /* Comment length */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, &comment_len[i]);
    }


    /* Convert our static structs to cached_items */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	if (in[i].bin) {
	    cached_item *si;
	    size_t extra_len;

	    extra_len = sizeof(anno_ele_t) + comment_len[i];
	    if (!(si = cache_new(GT_Seq, 0, 0, NULL, extra_len)))
		return NULL;

	    b->rec[i] = (rec << ANNO_ELE_BLOCK_BITS) + i;
	    b->ae[i]  = (anno_ele_t *)&si->data;
	    in[i].rec = (rec << ANNO_ELE_BLOCK_BITS) + i;
	    *b->ae[i] = in[i];
	    b->ae[i]->block = b;
	    b->ae[i]->idx = i;
	    b->ae[i]->comment = (char *)&b->ae[i]->data;
	} else {
	    b->rec[i] = 0;
	    b->ae[i] = NULL;
	}
    }


    /* Decode variable sized components */
    /* Comments */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	memcpy(b->ae[i]->comment, cp, comment_len[i]);
	b->ae[i]->comment[comment_len[i]] = 0;
	cp += comment_len[i];
    }

    assert(cp - buf == buf_len);
    free(buf);

    return ci;
}

static int io_anno_ele_block_write(void *dbh, cached_item *ci) {
    int err;
    g_io *io = (g_io *)dbh;
    anno_ele_block_t *b = (anno_ele_block_t *)&ci->data;
    int i, last_obj_rec, last_anno_rec;
    char *cp, *cp_start;
    char *out[7], *out_start[7];
    size_t out_size[7], total_size;
    int level[7];

    /* Compute worst-case sizes, for memory allocation */
    for (i = 0; i < 7; i++) {
	out_size[i] = 0;
    }
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	anno_ele_t *e = b->ae[i];
	if (!e) {
	    out_size[0]++;
	    continue;
	}

	out_size[0] += 5; /* bin */
	out_size[1] += 5; /* tag type */
	out_size[2] += 5; /* obj type */
	out_size[3] += 5; /* obj record */
	out_size[4] += 5; /* anno record */
	out_size[5] += 5; /* comment length */
	out_size[6] += e->comment ? strlen(e->comment) : 0; /* comments */
    }
    for (i = 0; i < 7; i++)
	out_start[i] = out[i] = malloc(out_size[i]+1);
    
    /* serialised annotations, separated by type of data */
    last_obj_rec = last_anno_rec = 0;
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	int delta;
	anno_ele_t *e = b->ae[i];

	if (!e) {
	    *out[0]++=0; /* signifies annotation not present */
	    continue;
	}

	out[0] += int2u7(e->bin, out[0]);
	out[1] += int2u7(e->tag_type, out[1]);
	out[2] += int2u7(e->obj_type, out[2]);

	delta = e->obj_rec - last_obj_rec;
	last_obj_rec = e->obj_rec;
	out[3] += int2s7(delta, out[3]);
	//out[3] += int2u7(e->obj_rec, out[3]);

	delta = e->anno_rec - last_anno_rec;
	last_anno_rec = e->anno_rec;
	out[4] += int2s7(delta, out[4]);
	//out[4] += int2u7(e->anno_rec, out[4]);

	if (e->comment) {
	    int comment_len = strlen(e->comment);
	    out[5] += int2u7(comment_len, out[5]);
	    memcpy(out[6], e->comment, comment_len);
	    out[6] += comment_len;
	} else {
	    out[5] += int2u7(0, out[5]);
	}
    }

    /* Concatenate data types together and adjust out_size to actual usage */
    for (total_size = i = 0; i < 7; i++) {
	out_size[i] = out[i] - out_start[i];
	total_size += out_size[i];
    }
    cp = cp_start = (char *)malloc(total_size+1);
    for (i = 0; i < 7; i++) {
	memcpy(cp, out_start[i], out_size[i]);
	cp += out_size[i];
	free(out_start[i]);
	level[i] = 6; /* default gzip level */
    }
    assert(cp - cp_start == total_size);

    level[6] = 8; /* comments */

    /* Gzip it too */
    if (1) {
	char *gzout;
	size_t ssz;

	//gzout = mem_deflate(cp_start, cp-cp_start, &ssz);
	gzout = mem_deflate_lparts(cp_start, out_size, level, 7, &ssz);
	free(cp_start);
	cp_start = gzout;
	cp = cp_start + ssz;
    }

    /* Finally write the serialised data block */
    assert(ci->lock_mode >= G_LOCK_RW);
    wrstats[GT_AnnoEleBlock] += cp-cp_start;
    wrcounts[GT_AnnoEleBlock]++;
    err = g_write(io, ci->view, cp_start, cp - cp_start);
    g_flush(io, ci->view);
    
    free(cp_start);

    return err ? -1 : 0;
}

static GRec io_anno_ele_block_create(void *dbh, void *vfrom) {
    g_io *io = (g_io *)dbh;
    GRec rec;
    GView v;

    rec = allocate(io);
    v = lock(io, rec, G_LOCK_EX);

    /* Write blank data here? */

    unlock(io, v);
    
    return rec;
}

/* ------------------------------------------------------------------------
 * The externally visible object interfaces themselves.
 */
static iface iface_g = {
    io_database_create_files,
    io_database_connect,
    io_database_disconnect,
    io_database_commit,
    io_database_lock,
    io_database_unlock,

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
	/* AnnoEle */
	io_anno_ele_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_anno_ele_read,
	io_anno_ele_write,
	io_generic_info,
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
	io_anno_write,
	io_generic_info,
    },

    {
	/* Library */
	io_generic_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_library_read,
	io_library_write,
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
	/* Seq_block */
	io_seq_block_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_seq_block_read,
	io_seq_block_write,
	io_generic_info,
    },

    {
	/* Anno_ele_block */
	io_anno_ele_block_create,
	io_generic_destroy,
	io_generic_lock,
	io_generic_unlock,
	io_generic_upgrade,
	io_generic_abandon,
	io_anno_ele_block_read,
	io_anno_ele_block_write,
	io_generic_info,
    },
};

iface *get_iface_g(void) {
    return &iface_g;
}
