#include <staden_config.h>
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
    int comp_mode;
} g_io;


/*
 * Dummy nul-compression functions. These do nothing except satisfy
 * the internal API requirements. We use them simply to compute the
 * base-line so we can estimate the overhead due to zlib vs overhead
 * due to lzma.
 */
static char *nul_mem_deflate(char *data, size_t size, size_t *cdata_size) {
    char *out = malloc(size);
    memcpy(out, data, size);
    *cdata_size = size;
    return out;
}

static char *nul_mem_deflate_parts(char *data,
				   size_t *part_size, int nparts,
				   size_t *cdata_size) {
    size_t tot_size = 0;
    int i;
    for (i = 0; i < nparts; i++)
	tot_size += part_size[i];
    return nul_mem_deflate(data, tot_size, cdata_size);
}

static char *nul_mem_deflate_lparts(char *data,
				    size_t *part_size, int *level, int nparts,
				    size_t *cdata_size) {
    return nul_mem_deflate_parts(data, part_size, nparts, cdata_size);
}

static char *nul_mem_inflate(char *cdata, size_t csize, size_t *size) {
    char *out = malloc(csize);
    memcpy(out, cdata, csize);
    *size = csize;
    return out;
}

#ifdef HAVE_LIBLZMA
/* ------------------------------------------------------------------------ */
/*
 * Data compression routines using liblzma (xz)
 *
 * On a test set this shrunk the main db from 136157104 bytes to 114796168, but
 * caused tg_index to grow from 2m43.707s to 15m3.961s. Exporting as bfastq
 * went from 18.3s to 36.3s. So decompression suffers too, but not as bad
 * as compression times.
 *
 * For now we disable this functionality. If it's to be reenabled make sure you
 * improve the mem_inflate implementation as it's just a test hack at the
 * moment.
 */
#include <lzma.h>

/* Fast mode, but not too bad on compression still */
#define LZMA_LEVEL 3

static char *lzma_mem_deflate(char *data, size_t size, size_t *cdata_size) {
    char *out;
    size_t out_size = lzma_stream_buffer_bound(size);
    *cdata_size = 0;

    out = malloc(out_size);

    /* Single call compression */
    if (LZMA_OK != lzma_easy_buffer_encode(LZMA_LEVEL,
					   LZMA_CHECK_CRC32, NULL,
					   (uint8_t *)data, size,
					   (uint8_t *)out, cdata_size,
					   out_size))
    	return NULL;

    return out;
}

static char *lzma_mem_deflate_parts(char *data,
				    size_t *part_size, int nparts,
				    size_t *cdata_size) {
    size_t tot_size = 0;
    int i;
    for (i = 0; i < nparts; i++)
	tot_size += part_size[i];
    return lzma_mem_deflate(data, tot_size, cdata_size);
}

static char *lzma_mem_deflate_lparts(char *data,
				     size_t *part_size, int *level, int nparts,
				     size_t *cdata_size) {
    return lzma_mem_deflate_parts(data, part_size, nparts, cdata_size);
}

static char *lzma_mem_inflate(char *cdata, size_t csize, size_t *size) {
    lzma_stream strm = LZMA_STREAM_INIT;
    size_t out_size = 0, out_pos = 0;
    char *out = NULL;
    int r;

    /* Initiate the decoder */
    if (LZMA_OK != lzma_stream_decoder(&strm, 50000000, 0))
	return NULL;

    /* Decode loop */
    strm.avail_in = csize;
    strm.next_in = (uint8_t *)cdata;

    for (;strm.avail_in;) {
	if (strm.avail_in > out_size - out_pos) {
	    out_size += strm.avail_in * 4 + 32768;
	    out = realloc(out, out_size);
	}
	strm.avail_out = out_size - out_pos;
	strm.next_out = (uint8_t *)&out[out_pos];

	r = lzma_code(&strm, LZMA_RUN);
	if (LZMA_OK != r && LZMA_STREAM_END != r) {
	    fprintf(stderr, "r=%d\n", r);
	    fprintf(stderr, "mem=%"PRId64"d\n", (int64_t)lzma_memusage(&strm));
	    return NULL;
	}

	out_pos = strm.total_out;

	if (r == LZMA_STREAM_END)
	    break;
    }

    /* finish up any unflushed data; necessary? */
    r = lzma_code(&strm, LZMA_FINISH);
    if (r != LZMA_OK && r != LZMA_STREAM_END) {
	fprintf(stderr, "r=%d\n", r);
	return NULL;
    }

    out = realloc(out, strm.total_out);
    *size = strm.total_out;

    lzma_end(&strm);

    return out;
}
#endif

/* ------------------------------------------------------------------------ */
/*
 * Data compression routines using zlib.
 */
#include <zlib.h>
static char *zlib_mem_deflate(char *data, size_t size, size_t *cdata_size) {
    z_stream s;
    unsigned char *cdata = NULL; /* Compressed output */
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
	s.next_out = &cdata[cdata_pos];
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
    return (char *)cdata;
}

static char *zlib_mem_deflate_parts(char *data,
				    size_t *part_size, int nparts,
				    size_t *cdata_size) {
    z_stream s;
    unsigned char *cdata = NULL; /* Compressed output */
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
	    s.next_out = &cdata[cdata_pos];
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
    return (char *)cdata;
}

static char *zlib_mem_deflate_lparts(char *data,
				     size_t *part_size, int *level, int nparts,
				     size_t *cdata_size) {
    z_stream s;
    unsigned char *cdata = NULL; /* Compressed output */
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
	    s.next_out = &cdata[cdata_pos];
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
    return (char *)cdata;
}

static char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size) {
    z_stream s;
    unsigned char *data = NULL; /* Uncompressed output */
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
	s.next_out = &data[s.total_out];
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
    return (char *)data;
}

static char *mem_deflate(int mode,
			 char *data, size_t size, size_t *cdata_size) {
    switch (mode) {
    case COMP_MODE_NONE:
	return nul_mem_deflate (data, size, cdata_size);
    case COMP_MODE_ZLIB:
	return zlib_mem_deflate(data, size, cdata_size);
#ifdef HAVE_LIBLZMA
    case COMP_MODE_LZMA:
	return lzma_mem_deflate(data, size, cdata_size);
#endif	
    }
    
    return NULL;
}

static char *mem_deflate_parts(int mode, char *data,
			       size_t *part_size, int nparts,
			       size_t *cdata_size) {
    switch (mode) {
    case COMP_MODE_NONE:
	return nul_mem_deflate_parts (data, part_size, nparts, cdata_size);
    case COMP_MODE_ZLIB:
	return zlib_mem_deflate_parts(data, part_size, nparts, cdata_size);
#ifdef HAVE_LIBLZMA
    case COMP_MODE_LZMA:
	return lzma_mem_deflate_parts(data, part_size, nparts, cdata_size);
#endif	
    }
    
    return NULL;
}

static char *mem_deflate_lparts(int mode, char *data,
				size_t *part_size, int *level, int nparts,
				size_t *cdata_size) {
    switch (mode) {
    case COMP_MODE_NONE:
	return nul_mem_deflate_lparts (data, part_size, level, nparts, cdata_size);
    case COMP_MODE_ZLIB:
	return zlib_mem_deflate_lparts(data, part_size, level, nparts, cdata_size);
#ifdef HAVE_LIBLZMA
    case COMP_MODE_LZMA:
	return lzma_mem_deflate_lparts(data, part_size, level, nparts, cdata_size);
#endif	
    }
    
    return NULL;
}

static char *mem_inflate(int mode,
			 char *cdata, size_t csize, size_t *size) {
    switch (mode) {
    case COMP_MODE_NONE:
	return nul_mem_inflate (cdata, csize, size);
    case COMP_MODE_ZLIB:
	return zlib_mem_inflate(cdata, csize, size);
#ifdef HAVE_LIBLZMA
    case COMP_MODE_LZMA:
	return lzma_mem_inflate(cdata, csize, size);
#endif	
    }
    
    return NULL;
}

/* ------------------------------------------------------------------------ */
/*
 * Simple interfaces to the underlying g-library. These are very basic, but
 * sufficient. Will need to rework this with proper lock utilisation schemes
 * for Gap4 proper.
 */

/* Hacky allocation scheme - always increment the record number */
#if 1
static int other_record = 0;
static int other_record_start = 0;
void set_reserved_seqs(int rseqs) {
    other_record = rseqs / SEQ_BLOCK_SZ;
    other_record_start = other_record;
}

static tg_rec allocate(g_io *io, int type) {
    static tg_rec record = -1;
    tg_rec r;
    
    if (record == -1)
	record = io->gdb->gfile->header.num_records;

    /*
     * Temporary hack to prevent lots of small records (bins, contigs, etc)
     * from pushing the number of bits for the record number beyond
     * 31-SEQ_BLOCK_BITS (ie currently record ~2 million).
     *
     * FIXME: the correct solution is to enable record numbers to be 64-bit.
     */
    if (other_record) {
	switch (type) {
	case GT_SeqBlock:
	case GT_AnnoEleBlock:
	case GT_Database:
	    r = record++;
	    if (r == other_record_start) {
		fprintf(stderr, "\n*** Ran out of seq/anno record numbers.\n");
		fprintf(stderr, "*** Please use a higher value in the -r "
			"option of tg_index.\n");
		exit(1);
	    }
	    break;

	default:
	    r = other_record++;
	}
    } else {
	if ((r = g_free_rec_(io->gdb, io->client, 0)) != G_NO_REC)
	    return r;

	r = record++;
#if 0
	if (type == GT_SeqBlock && r >= (1<<(31-SEQ_BLOCK_BITS))) {
	    fprintf(stderr, "\n*** Too many database records to cope with the"
		    " sequence 'blocking factor'.");
	    fprintf(stderr, "*** Please rerun tg_index using the -r option "
		    "to reserve record space for\n    more sequences\n");
	    exit(1);
	} else if (type == GT_AnnoEleBlock &&
		   r >= (1<<(31-ANNO_ELE_BLOCK_BITS))) {
	    fprintf(stderr, "\n*** Too many database records to cope with the"
		    " annotation 'blocking factor'.\n");
	    fprintf(stderr, "*** Please rerun tg_index using the -r option "
		    "to reserve record space for\n    more annotations\n");
	    exit(1);
	}
#endif
    }

    return r;
}
#else
static int allocate(g_io *io, int type) {
    /*
     * FIXME, the header should point to a record free-list, to allow
     * reclaiming of deallocated records.
     */
    static int record = -1;

    if (record == -1)
	record = io->gdb->gfile->header.num_records;

    return record++;
}
#endif

static int deallocate(g_io *io, tg_rec rec, GView v) {
    return g_remove_(io->gdb, io->client, v);
}

static GView lock(g_io *io, tg_rec rec, int mode) {
    if (!mode)
	mode = G_LOCK_EX;

    return g_lock_N_(io->gdb, io->client, 0, (GRec)rec, mode);
}

static int unlock(g_io *io, GView v) {
    return g_unlock_(io->gdb, io->client, v);
}

static int g_write(g_io *io, GView v, void *buf, size_t len) {
    return g_write_(io->gdb, io->client, v, buf, len);
}

static int g_writev(g_io *io, GView v, GIOVec *vec, GCardinal vcnt) {
    return g_writev_(io->gdb, io->client, v, vec, vcnt);
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

    if (!vi.used)
	return NULL;

    if (NULL == (buf = malloc(vi.used)))
	return NULL;

    if (g_read_(io->gdb, io->client, v, buf, vi.used) == 0)
	return buf;
    else {
	free(buf);
	return NULL;
    }
}

static int g_flush(g_io *io, GView v) {
    return g_flush_(io->gdb, io->client, v);
}


/* ------------------------------------------------------------------------
 * Generic io_database functions
 */


/* Generic functions - shared by all objects within the g library */
static tg_rec io_generic_create(void *dbh, void *unused) {
    return allocate((g_io *)dbh, GT_Generic);
}

static int io_generic_destroy(void *dbh, tg_rec r, GView v) {
    return deallocate((g_io *)dbh, r, v);
}

static GView io_generic_lock(void *dbh, tg_rec r, int mode) {
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


/*
 * Generic reading and writing of N 32-bit integer values
 */
static GCardinal *io_generic_read_i4(g_io *io, GView v, int type,
				     size_t *nitems) {
    unsigned char *buf, *cp;
    size_t buf_len;
    uint32_t i, ni;
    GCardinal *card;

    /* Load from disk */
    cp = buf = g_read_alloc(io, v, &buf_len);
    if (buf_len < 2) {
	*nitems = 0;
	return NULL;
    }

    assert(cp[0] == type);
    assert((cp[1] & 0x3f) == 0); /* initial format */
    cp += 2;
    cp += u72int(cp, &ni);
    *nitems = ni;

    if (NULL == (card = (GCardinal *)malloc(*nitems * sizeof(GCardinal)))) {
	free(buf);
        return NULL;
    }

    for (i = 0; i < ni; i++)
	cp += u72int(cp, (uint32_t *)&card[i]);

    assert(cp-buf == buf_len);
    free(buf);

    return card;
}

/*
 * Generic reading and writing of N tg_rec integer values
 */
static tg_rec *io_generic_read_rec(g_io *io, GView v, int type,
				   size_t *nitems) {
    unsigned char *buf, *cp;
    size_t buf_len;
    uint64_t i, ni, i64;
    uint32_t i32;
    tg_rec *recs;
    int fmt;

    /* Load from disk */
    cp = buf = g_read_alloc(io, v, &buf_len);
    if (buf_len < 2) {
	*nitems = 0;
	return NULL;
    }

    assert(cp[0] == type);
    fmt = (cp[1] & 0x3f);
    assert(fmt <= 1);

    cp += 2;
    if (fmt == 0) {
	cp += u72intw(cp, &i64);
	ni = i64;
    } else {
	cp += u72int(cp, &i32);
	ni = i32;
    }
    *nitems = ni;

    if (NULL == (recs = (tg_rec *)malloc(*nitems * sizeof(*recs)))) {
	free(buf);
        return NULL;
    }

    if (fmt == 0) {
	for (i = 0; i < ni; i++) {
	    cp += u72int(cp, &i32);
	    recs[i] = i32;
	}
    } else {
	for (i = 0; i < ni; i++) {
	    cp += u72intw(cp, &i64);
	    recs[i] = i64;
	}
    }

    assert(cp-buf == buf_len);
    free(buf);

    return recs;
}

static cached_item *io_generic_read(void *dbh, tg_rec rec, int type) {
    GView v;
    unsigned char *buf, *cp;
    size_t buf_len;
    cached_item *ci;
    uint32_t nitems, i;
    GCardinal *card;

    /* Load from disk */
    if (-1 == (v = io_generic_lock(dbh, rec, G_LOCK_RO)))
	return NULL;

    cp = buf = g_read_alloc((g_io *)dbh, v, &buf_len);
    if (buf_len < 2)
	return NULL;

    assert(cp[0] == type);
    assert((cp[1] & 0x3f) == 0); /* initial format */
    cp += 2;
    cp += u72int(cp, &nitems);

    if (!(ci = cache_new(type, rec, v, NULL, nitems * sizeof(GCardinal)))) {
	free(buf);
        return NULL;
    }
    ci->data_size = nitems * sizeof(GCardinal);
    card = (GCardinal *)&ci->data;

    for (i = 0; i < nitems; i++)
	cp += u72int(cp, (uint32_t *)&card[i]);

    assert(cp-buf == buf_len);
    free(buf);

    return ci;
}

/*
 * Returns number of bytes written on success
 *         -1 on failure
 */
static int io_generic_write_i4(g_io *io, GView v, int type,
			       void *buf, size_t len) {
    int ret, i;
    unsigned char *cp_start, *cp;
    GCardinal *card = (GCardinal *)buf;
    uint32_t nitems = len / sizeof(GCardinal);

    /* Allocate memory based on worst case sizes */
    if (NULL == (cp = cp_start = malloc(5 * nitems + 2 + 5)))
	return -1;

    *cp++ = type;
    *cp++ = 0;
    cp += int2u7(nitems, cp);
    for (i = 0; i < nitems; i++) {
	cp += int2u7((uint32_t)card[i], cp);
    }

    ret = g_write(io, v, cp_start, cp - cp_start);
    g_flush(io, v); /* Should we auto-flush? */

    free(cp_start);
    return ret ? -1 : cp - cp_start;
}

/*
 * Generic reading and writing of N tg_rec integer values
 */
static int io_generic_write_rec(g_io *io, GView v, int type,
				tg_rec *recs, size_t nrec) {
    int ret, i;
    unsigned char *cp_start, *cp;

    /* Allocate memory based on worst case sizes */
    if (NULL == (cp = cp_start = malloc(10 * nrec + 2 + 10)))
	return -1;

    *cp++ = type;
    *cp++ = 1; /* fmt=0 for 32-bit only, fmt=1 for 64-bit recs */
    cp += intw2u7(nrec, cp);
    for (i = 0; i < nrec; i++) {
	cp += intw2u7((uint64_t)recs[i], cp);
    }

    ret = g_write(io, v, cp_start, cp - cp_start);
    g_flush(io, v); /* Should we auto-flush? */

    free(cp_start);
    return ret ? -1 : cp - cp_start;
}

static int io_generic_write(void *dbh, cached_item *ci) {
    assert(ci->lock_mode >= G_LOCK_RW);

    return io_generic_write_i4((g_io *)dbh, ci->view,
			       ci->type, &ci->data, ci->data_size) ? 0 : -1;
}

static int io_generic_info(void *dbh, GView v, GViewInfo *vi) {
    g_io *io = (g_io *)dbh;
    return g_view_info_(io->gdb, io->client, v, vi);
}

/*
 * Returns 1 if the record / type combination is valid.
 *         0 if not.
 */
static int io_rec_exists(void *dbh, int type, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    GView v;
    unsigned char buf;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return 0;
    
    buf = 0;
    if (0 != g_read(dbh, v, &buf, 1))
	return 0;

    return buf == type ? 1 : 0;
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
    char *buf, *buf2;
    btree_node_t *n;
    BTRec rec = *((BTRec *)key);
    static HacheData hd;
    int fmt;
    int comp_mode;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    if (NULL == (buf = g_read_alloc(io, v, &len))) {
	return NULL;
    }

    assert(buf[0] == GT_BTree);
    fmt = buf[1] & 0x3f;
    assert(fmt <= 2); /* format number */
    comp_mode = ((unsigned char)buf[1]) >> 6;

    if (fmt >= 1) {
	size_t ssz;
	char *unpacked = mem_inflate(comp_mode, buf+2, len-2, &ssz);

	free(buf);
	buf2 = buf = unpacked;
	len = ssz+2;
    } else {
	buf2 = buf+2;
    }

    rdstats[GT_BTree] += len;
    rdcounts[GT_BTree]++;

    /* Decode the btree element */
    switch (fmt) {
    case 0:
	n = btree_node_decode((unsigned char *)buf2);
	break;
    case 1:
    case 2:
	n = btree_node_decode2((unsigned char *)buf2, fmt);
	break;
    default:
	abort();
    }
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
    size_t len, gzlen;
    size_t parts[4];
    GView v;
    int bfmt = 2;
    char *data = (char *)btree_node_encode2(n, &len, parts, bfmt);
    char fmt[2];
    GIOVec vec[2];
    cached_item *ci = n->cache;
    char *gzout;

    /* Set up data type and version */
    fmt[0] = GT_BTree;
    fmt[1] = bfmt | (io->comp_mode << 6);
    vec[0].buf = fmt;  vec[0].len = 2;

    gzout = mem_deflate_parts(io->comp_mode, data, parts, 4, &gzlen);
    free(data); data = gzout; len = gzlen;

    vec[1].buf = data; vec[1].len = len;

    if (ci) {
	wrstats[GT_BTree] += len;
	wrcounts[GT_BTree]++;
	//ret = g_write(io, ci->view, b2, len);
	ret = g_writev(io, ci->view, vec, 2);
	g_flush(io, ci->view);
    } else {
	if (-1 == (v = lock(io, n->rec, G_LOCK_EX))) {
	    fprintf(stderr, "Failed to lock btree node %"PRIbtr"\n", n->rec);
	    return -1;
	}
	wrstats[GT_BTree] += len;
	wrcounts[GT_BTree]++;
	//ret = g_write(io, v, b2, len);
	ret = g_writev(io, v, vec, 2);
	//unlock(io, v);
    }

    free(data);

    if (ret) {
	fprintf(stderr, "Failed to write btree node %"PRIbtr"\n", n->rec);
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

    if (hi) {
	return (btree_node_t *)(((cached_item *)hi->data.p)->data);
    } else {
	fprintf(stderr, "Failed to load btree %"PRIbtr"\n", r);
	return NULL;
    }
}

/*
 * Create a node and add it to our internal cache too
 *
 * Returns the record number on success
 *         -1 on failure
 */
tg_rec btree_node_create(g_io *io, HacheTable *h) {
    tg_rec rec;
    btree_node_t *n;
    cached_item *ci;
    GView v;
    static HacheData hd;
    HacheItem *hi;

    /* Allocate a new record */
    rec = allocate(io, GT_BTree);
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
    tg_rec rec = btree_node_create(bt->io, bt->h);
    
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
    auxheader.free_record = G_NO_REC;
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
    io->comp_mode = COMP_MODE_ZLIB;

    return io;
}

int io_database_setopt(void *dbh, io_opt opt, int val) {
    g_io *io = (g_io *)dbh;

    switch (opt) {
    case OPT_COMP_MODE:
	io->comp_mode = val;
	return 0;

    default:
	fprintf(stderr, "Unknown io_option: %d\n", (int)opt);
    }

    return -1;
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
	    free(io->seq_name_tree);
	}
    }

    if (io->contig_name_hash) {
	btree_destroy(io, io->contig_name_hash);
	if (io->contig_name_tree) {
	    free(io->contig_name_tree);
	}
    }

    g_disconnect_client_(io->gdb, io->client);
    g_shutdown_database_(io->gdb);

    free(io);

    printf("\n*** I/O stats (type, write count/size read count/size) ***\n");
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


static cached_item *io_database_read(void *dbh, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    database_t *db;
    cached_item *ci;
    GView v;
    unsigned char *buf, *cp;
    size_t buf_len;
    uint32_t nitems, i32;
    int fmt;

    /* Load from disk */
    if (-1 == (v = io_generic_lock(dbh, rec, G_LOCK_RO)))
	return NULL;

    cp = buf = g_read_alloc((g_io *)dbh, v, &buf_len);
    if (buf_len < 2)
	return NULL;

    fmt = cp[1] & 0x3f;
    assert(cp[0] == GT_Database);
    assert(fmt <= 1); /* initial format */
    cp += 2;

    if (fmt == 0) {
	cp += u72int(cp, &nitems); /* ignore now */
    }

    /* Generate in-memory data structure */
    if (!(ci = cache_new(GT_Database, rec, v, NULL, sizeof(*db))))
	return NULL;
    db = (database_t *)&ci->data;

    cp += u72int(cp, (uint32_t *)&db->version);
    cp += u72int(cp, (uint32_t *)&db->Ncontigs);
    cp += u72int(cp, &i32); db->contig_order = i32;
    cp += u72int(cp, (uint32_t *)&db->Nlibraries);
    cp += u72int(cp, &i32); db->library = i32;
    cp += u72int(cp, &i32); db->seq_name_index = i32;
    cp += u72int(cp, &i32); db->contig_name_index = i32;

    assert(cp-buf == buf_len);
    free(buf);

#ifdef INDEX_NAMES
    /* Initialise the seq_name btree if needed */
    if (io->seq_name_tree)
	return ci;

    /* Read the root */
    if (db->seq_name_index) {
	btree_query_t *bt = (btree_query_t *)io->seq_name_hash->clientdata;
	bt->io = io;
	bt->h = io->seq_name_hash;
	io->seq_name_tree = btree_new(bt, db->seq_name_index);
	if (!io->seq_name_tree || !io->seq_name_tree->root)
	    return NULL;
    }

    //printf("seq_name_hash=%p\n", io->seq_name_hash);
#endif

    if (db->contig_name_index) {
	btree_query_t *bt = (btree_query_t *)io->contig_name_hash->clientdata;
	bt->io = io;
	bt->h = io->contig_name_hash;
	io->contig_name_tree = btree_new(bt, db->contig_name_index);
	if (!io->contig_name_tree || !io->contig_name_tree->root)
	    return NULL;
    }

    return ci;

}


static int io_database_write_view(g_io *io, database_t *db, GView v) {
    unsigned char buf[sizeof(*db)*2], *cp = buf;

    /* Construct the on-disc format */
    *cp++ = GT_Database;
    *cp++ = 1; /* format */

    cp += int2u7(db->version, cp);
    cp += int2u7(db->Ncontigs, cp);
    cp += intw2u7(db->contig_order, cp);
    cp += int2u7(db->Nlibraries, cp);
    cp += intw2u7(db->library, cp);
    cp += intw2u7(db->seq_name_index, cp);
    cp += intw2u7(db->contig_name_index, cp);
    
    /* Write it out */
    if (-1 == g_write(io, v, buf, cp-buf))
	return -1;

    g_flush(io, v);
    return 0;
}

/*
 * Writes a database_t record.
 * Returns 0 on success
 *        -1 on failure
 */
static int io_database_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    database_t *db = (database_t *)&ci->data;;

    assert(ci->lock_mode >= G_LOCK_RW);
    return io_database_write_view(io, db, ci->view);
}


static tg_rec io_database_create(void *dbh, void *from) {
    g_io *io = (g_io *)dbh;
    tg_rec db_rec = allocate(io, GT_Database);
    GView v;
    database_t db;

    /* init_db is only called on a blank database => first record is 0 */
    assert(db_rec == 0);

    db.Ncontigs = 0;
    db.version = DB_VERSION;

    /* Contig order */
    db.contig_order = allocate(io, GT_RecArray); /* contig array */
    v = lock(io, db.contig_order, G_LOCK_EX);
    //    if (-1 == g_write_le4(io, v, NULL, 0))
    //	return -1;
    g_flush(io, v);
    unlock(io, v);

    /* Libraries */
    db.Nlibraries = 0;
    db.library = allocate(io, GT_RecArray);
    v = lock(io, db.library, G_LOCK_EX);
    //    if (-1 == g_write_le4(io, v, NULL, 0))
    //	return -1;
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
    if (-1 == io_database_write_view(io, &db, v))
	return -1;

    unlock(io, v);

    io_database_commit(io);
   
    return 0;
}

/*
 * Creates an B+Tree index, but does not attach it to a specific object.
 */
int io_database_create_index(void *dbh, cached_item *ci, int type) {
    g_io *io = (g_io *)dbh;
    HacheTable *h = HacheTableCreate(1024, HASH_DYNAMIC_SIZE | HASH_OWN_KEYS);
    btree_query_t *bt;
    database_t *db = (database_t *)&ci->data;

    if (NULL == (bt = (btree_query_t *)malloc(sizeof(*bt))))
	return -1;

    bt->io = io;
    bt->h  = h;

    h->clientdata = bt;
    h->load       = btree_load_cache;
    h->del        = btree_del_cache;

    switch(type) {
    case DB_INDEX_NAME:
	if (db->seq_name_index)
	    return -1; /* already exists */

	io->seq_name_hash = h;
	h->name = "io->seq_name_hash";
	db->seq_name_index = btree_node_create(io, h);
	io->seq_name_tree = btree_new(bt, db->seq_name_index);

	assert(io->seq_name_tree);
	assert(io->seq_name_tree->root);
	break;
	
    case DB_INDEX_CONTIG:
	if (db->contig_name_index)
	    return -1; /* already exists */

	io->contig_name_hash = h;
	h->name = "io->contig_name_hash";
	db->contig_name_index = btree_node_create(io, h);
	io->contig_name_tree = btree_new(bt, db->contig_name_index);

	assert(io->contig_name_tree);
	assert(io->contig_name_tree->root);
	break;

    default:
	return -1;
    }

    io_database_commit(io);

    return 0;
}

/* ------------------------------------------------------------------------
 * contig_t access methods
 */
static cached_item *io_contig_read(void *dbh, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    contig_t *c;
    size_t len, slen;
    unsigned char *ch, *cp;
    int32_t start, end;
    uint32_t bin, nlen;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    if (NULL == (cp = ch = g_read_alloc(io, v, &len))) {
	return NULL;
    }

    if (len < 2)
	return NULL;

    assert(cp[0] == GT_Contig);
    assert(cp[1] == 0);
    cp += 2;

    rdstats[GT_Contig] += len;
    rdcounts[GT_Contig]++;

    /* Decode the fixed size bits */
    cp += s72int(cp, &start);
    cp += s72int(cp, &end);
    cp += u72int(cp, &bin);
    cp += u72int(cp, &nlen);

    /* Generate in-memory data structure */
    slen = sizeof(*c) + sizeof(char *) + nlen + 1;
    if (!(ci = cache_new(GT_Contig, rec, v, NULL, slen)))
	return NULL;

    c = (contig_t *)&ci->data;
    c->rec    = rec;
    c->start  = start;
    c->end    = end;
    c->bin    = bin;
    c->name   = (char *)&c->data;
    memcpy(c->name, cp, nlen);
    c->name[nlen] = 0;

    free(ch);

    return ci;
}

static int io_contig_write_view(g_io *io, contig_t *c, GView v) {
    size_t len;
    unsigned char *cp, *buf;
    int nlen;

    /* Estimate worst-case memory requirements */
    nlen = c->name ? strlen(c->name) : 0;
    len = 2 + 5+5+5 + 5+nlen;
    if (NULL == (cp = buf = malloc(len)))
	return -1;


    /* Construct on-disc representation */
    *cp++ = GT_Contig;
    *cp++ = 0; /* format version */
    cp += int2s7(c->start, cp);
    cp += int2s7(c->end, cp);
    cp += int2u7((int32_t)c->bin, cp);
    cp += int2u7(nlen, cp);
    if (c->name) {
	memcpy(cp, c->name, nlen);
	cp += nlen;
    }
    len = cp-buf; /* Actual length */

    /* Write the data */
    wrstats[GT_Contig] += len;
    wrcounts[GT_Contig]++;
    if (-1 == g_write(io, v, (char *)buf, len)) {
	free(buf);
	return -1;
    }

    g_flush(io, v);
    free(buf);
    return 0;
}

static int io_contig_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    contig_t *c = (contig_t *)&ci->data;

    assert(ci->lock_mode >= G_LOCK_RW);
    assert(ci->rec > 0);
    return io_contig_write_view(io, c, ci->view);
}

static tg_rec io_contig_create(void *dbh, void *vfrom) {
    g_io *io = (g_io *)dbh;
    tg_rec rec;
    GView v;
    contig_t *from = (contig_t *)vfrom;

    rec = allocate(io, GT_Contig);
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

static tg_rec io_contig_index_query(void *dbh, char *name, int prefix) {
    g_io *io = (g_io *)dbh;
    
    if (!io->contig_name_tree)
	return -1;

    return btree_search(io->contig_name_tree, name, prefix);
}

static int io_contig_index_add(void *dbh, char *name, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    
    if (!io->contig_name_tree)
	return -1;

    btree_insert(io->contig_name_tree, name, rec);
    return io->contig_name_tree->root->rec;
}

static int io_contig_index_del(void *dbh, char *name) {
    g_io *io = (g_io *)dbh;
    
    if (!io->contig_name_tree)
	return -1;

    return btree_delete(io->contig_name_tree, name);
}

/* ------------------------------------------------------------------------
 * Array access methods
 */
static cached_item *io_array_read(void *dbh, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    GViewInfo vi;
    ArrayStruct *ar;

    if (!(ci = cache_new(GT_RecArray, rec, 0, NULL, sizeof(*ar))))
	return NULL;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    g_view_info_(io->gdb, io->client, v, &vi);
    rdstats[GT_RecArray] += vi.used;
    rdcounts[GT_RecArray]++;

    ar = ArrayCreate(sizeof(tg_rec), 0);
    if (ar->base) free(ar->base);
    ar->base = io_generic_read_rec(io, v, GT_RecArray, &ar->dim);
    ar->max = ar->dim;

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
    assert(ci->rec > 0);
    ar = (Array)&ci->data;
    ret = io_generic_write_rec(io, ci->view, GT_RecArray,
			       ArrayBase(tg_rec, ar),
			       ArrayMax(ar));

    wrstats[GT_RecArray] += ret;
    wrcounts[GT_RecArray]++;

    return ret >= 0 ? 0 : -1;
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
static cached_item *io_anno_ele_read(void *dbh, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    void *bloc;
    unsigned char *cp;
    size_t bloc_len;
    GView v;
    cached_item *ci;
    uint32_t anno_rec, tag_type, obj_type, obj_rec, comment_len, bin_rec;
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
    assert(cp[0] == GT_AnnoEle);
    assert(cp[1] == 0); /* format */
    cp += 2;
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
    data_len = 2 + 5 + 5 + 5 + 5 + 5 + comment_len;
    if (data_len > 1024) {
	if (NULL == (cp = (unsigned char *)malloc(data_len)))
	    return -1;
    }

    /* Encode */
    cpstart = cp;
    *cp++ = GT_AnnoEle;
    *cp++ = 0; /* format */ 
    cp += int2u7((int32_t)e->bin, cp);
    cp += int2u7(e->tag_type, cp);
    cp += int2u7(e->obj_type, cp);
    cp += intw2u7(e->obj_rec, cp);
    cp += intw2u7(e->anno_rec, cp);
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
    assert(ci->rec > 0);
    return io_anno_ele_write_view(io, e, ci->view);
}

static tg_rec io_anno_ele_create(void *dbh, void *vfrom) {
    g_io *io = (g_io *)dbh;
    anno_ele_t *from = (anno_ele_t *)vfrom;
    tg_rec rec;
    GView v;

    rec = allocate(io, GT_AnnoEle);
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
static cached_item *io_anno_read(void *dbh, tg_rec rec) {
    return io_generic_read(dbh, rec, GT_Anno);
}

static int io_anno_write(void *dbh, cached_item *ci) {
    assert(ci->rec > 0);
    return io_generic_write(dbh, ci);
}


/* ------------------------------------------------------------------------
 * dnasrc access methods
 */
static cached_item *io_library_read(void *dbh, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    library_t *lib, l;
    char *ch, *zpacked;
    size_t len, ssz;
    int fmt = -1, comp_mode;
    char *name = NULL;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    ch = g_read_alloc(io, v, &len);

    if (ch && len) {
	assert(ch[0] == GT_Library);
	fmt = ch[1] & 0x3f;
	comp_mode = ((unsigned char)ch[1]) >> 6;
	assert(fmt >= 0 && fmt <= 1); /* format */

	zpacked = mem_inflate(comp_mode, ch+2, len-2, &ssz);
	free(ch);
	len = ssz;
	ch = zpacked;
    }

    /* Generate a static in-memory data structure */
    if (ch == NULL || len == 0) {
	l.insert_size[0] = 0;
	l.insert_size[1] = 0;
	l.insert_size[2] = 0;
	l.sd[0] = 0;
	l.sd[1] = 0;
	l.sd[2] = 0;
	l.machine = 0;
	l.lib_type = 0;
	l.name = NULL;
	memset(l.size_hist, 0, 3 * LIB_BINS * sizeof(l.size_hist[0][0]));
    } else {
	int i, j;
	uint32_t tmp;
	unsigned char *cp = (unsigned char *)ch;

	cp += u72int(cp, (uint32_t *)&l.insert_size[0]);
	cp += u72int(cp, (uint32_t *)&l.insert_size[1]);
	cp += u72int(cp, (uint32_t *)&l.insert_size[2]);
	cp += u72int(cp, &tmp); l.sd[0] = tmp/100.0;
	cp += u72int(cp, &tmp); l.sd[1] = tmp/100.0;
	cp += u72int(cp, &tmp); l.sd[2] = tmp/100.0;
	cp += u72int(cp, (uint32_t *)&l.machine);
	cp += u72int(cp, (uint32_t *)&l.lib_type);
	
	for (j = 0; j < 3; j++) {
	    int last = 0;
	    for (i = 0; i < LIB_BINS; i++) {
		cp += s72int(cp, &l.size_hist[j][i]);
		l.size_hist[j][i] += last;
		last = l.size_hist[j][i];
	    }
	}

	if (fmt && *cp) {
	    name = (char *)cp;
	}
    }

    /* Copy over to a dynamically allocated cache item */
    if (!(ci = cache_new(GT_Library, rec, v, NULL, sizeof(*lib) +
			 (name ? strlen(name)+1 : 0))))
	return NULL;
    lib = (library_t *)&ci->data;
    memcpy(lib, &l, sizeof(l));
    lib->rec = rec;

    if (name) {
	lib->name = (char *)&lib->data;
	strcpy(lib->name, name);
    } else {
	lib->name = NULL;
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
    unsigned char cpstart[LIB_BINS*5*3+100], *cp = cpstart;
    int tmp, i, j, err;
    char *gzout;
    size_t ssz;
    char fmt[2];
    GIOVec vec[2];

    assert(ci->lock_mode >= G_LOCK_RW);
    assert(ci->rec > 0);

    fmt[0] = GT_Library;
    fmt[1] = (lib->name ? 1 : 0) | (io->comp_mode << 6);

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
    if (lib->name) {
	strcpy((char *)cp, lib->name);
	cp += strlen(lib->name)+1;
    }

    /* Compress it */
    gzout = mem_deflate(io->comp_mode, (char *)cpstart, cp-cpstart, &ssz);
    //err = g_write(io, ci->view, cpstart, cp-cpstart);
    vec[0].buf = fmt;   vec[0].len = 2;
    vec[1].buf = gzout; vec[1].len = ssz;
    err = g_writev(io, ci->view, vec, 2);
    free(gzout);
    g_flush(io, ci->view);

    return err;
}

/* ------------------------------------------------------------------------
 * vector access methods
 */
static cached_item *io_vector_read(void *dbh, tg_rec rec) {
    return io_generic_read(dbh, rec, 0 /*GT_Vector*/);
}


/* ------------------------------------------------------------------------
 * bin access methods
 */
static char *pack_rng_array(int comp_mode, int fmt,
			    GRange *rng, int nr, int *sz) {
    int i;
    size_t part_sz[7];
    GRange last, last_tag;
    unsigned char *cp[6], *cp_orig[6], *out;
    char *out_orig;
    //char *cpt, *cpt_orig;
    //HacheTable *h = HacheTableCreate(16, HASH_DYNAMIC_SIZE);
    //int ntags;

    memset(&last, 0, sizeof(last));
    memset(&last_tag, 0, sizeof(last_tag));

    /* Pack the 6 structure elements to their own arrays */
    for (i = 0; i < 6; i++)
	cp[i] = cp_orig[i] = malloc(nr * 10);

    for (i = 0; i < nr; i++) {
	GRange r = rng[i];

	if (r.flags & GRANGE_FLAG_UNUSED) {
	    if (fmt == 0)
		cp[2] += int2u7((int32_t)r.rec, cp[2]);
	    cp[4] += int2u7(GRANGE_FLAG_UNUSED, cp[4]);
	    continue;
	}
	
	//printf("%04d %d\t%d\t", i, r.rec, r.start);

	if ((r.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    r.end   -= r.start;
	    r.start -= last_tag.start;
	    r.rec   -= last_tag.rec;
	} else {
	    r.end   -= r.start;
	    r.start -= last.start;
	    r.rec   -= last.rec;
	}

	if (fmt == 0) {
	    cp[0] += int2u7(r.start, cp[0]);
	    cp[1] += int2u7(r.end,   cp[1]);
	    cp[2] += int2u7((int32_t)r.rec,   cp[2]);
	    cp[3] += int2u7(r.mqual, cp[3]);
	    cp[4] += int2u7(r.flags, cp[4]);

	    if ((r.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
		if (!(r.flags & GRANGE_FLAG_TYPE_SINGLE))
		    cp[5] += int2s7(r.pair_rec - last_tag.pair_rec, cp[5]);
		last_tag = rng[i];
	    } else {
		if (!(r.flags & GRANGE_FLAG_TYPE_SINGLE))
		    cp[5] += int2s7(r.pair_rec - last.pair_rec, cp[5]);
		last = rng[i];
	    }
	} else {
	    cp[0] += int2s7 (r.start, cp[0]);
	    cp[1] += int2u7 (r.end,   cp[1]);
	    cp[2] += intw2s7(r.rec,   cp[2]);
	    cp[3] += int2u7 (r.mqual, cp[3]);
	    cp[4] += int2u7 (r.flags, cp[4]);

	    if ((r.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
		if (!(r.flags & GRANGE_FLAG_TYPE_SINGLE))
		    cp[5] += intw2s7(r.pair_rec - last_tag.pair_rec, cp[5]);
		last_tag = rng[i];
	    } else {
		if (!(r.flags & GRANGE_FLAG_TYPE_SINGLE))
		    cp[5] += intw2s7(r.pair_rec - last.pair_rec, cp[5]);
		last = rng[i];
	    }
	}
    }

    for (i = 0; i < 6; i++)
	part_sz[i+1] = cp[i]-cp_orig[i];

    /* Construct a header with nr and the size of the 6 packed struct fields */
    *sz =  7*5 + part_sz[1] + part_sz[2] + part_sz[3] +
	part_sz[4] + part_sz[5] + part_sz[6];

    out = malloc(*sz);
    out_orig = (char *)out;
    out += int2u7(nr, out);
    for (i = 0; i < 6; i++)
	out += int2u7(cp[i]-cp_orig[i], out);
    part_sz[0] = (char *)out-out_orig;

    /* Followed by the serialised 6 packed fields themselves */
    for (i = 0; i < 6; i++) {
	int len = cp[i]-cp_orig[i];
	memcpy(out, cp_orig[i], len);
	out += len;
	free(cp_orig[i]);
    }

    *sz = (char *)out-out_orig;

    /* Gzip it too */
    {
    	char *gzout;
	size_t ssz;

	if (*sz < 512)
	    gzout = mem_deflate(comp_mode, out_orig, *sz, &ssz);
	else
	    gzout = mem_deflate_parts(comp_mode, out_orig, part_sz, 7, &ssz);
	*sz = ssz;

    	free(out_orig);
    	out_orig = gzout;
    }

    //write(2, out_orig, *sz);

    //HacheTableDestroy(h, 0);

    return out_orig;
}

static GRange *unpack_rng_array(int comp_mode, int fmt,
				unsigned char *packed,
				int packed_sz, int *nr) {
    uint32_t i, off[6];
    int32_t i32;
    unsigned char *cp[6], *zpacked = NULL;
    GRange last, *r, *ls = &last, *lt = &last;
    size_t ssz;

    /* First of all, inflate the compressed data */
    zpacked = packed = (unsigned char *)mem_inflate(comp_mode,
						    (char *)packed,
						    packed_sz, &ssz);
    packed_sz = ssz;

    /* Unpack number of ranges */
    packed += u72int(packed, (uint32_t *)nr);

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
	r[i].y = 0;

	if (fmt == 0) {
	    cp[2] += u72int(cp[2], &i32); r[i].rec = i32;
	}
	cp[4] += u72int(cp[4], (uint32_t *)&r[i].flags);
	if (r[i].flags & GRANGE_FLAG_UNUSED) {
	    if (fmt != 0)
		r[i].rec = 0;
	    r[i].start = 0;
	    r[i].end = 0;
	    r[i].mqual = 0;
	    r[i].pair_rec = 0;
	    continue;
	}
	
	if (fmt == 0) {
	    cp[0] += u72int(cp[0], (uint32_t *)&r[i].start);
	    cp[1] += u72int(cp[1], (uint32_t *)&r[i].end);
	    cp[3] += u72int(cp[3], (uint32_t *)&r[i].mqual);
	    if ((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
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
	} else {
	    int64_t rec_tmp;
	    cp[0] += s72int (cp[0], (int32_t *)&r[i].start);
	    cp[1] += u72int (cp[1], (uint32_t *)&r[i].end);
	    cp[2] += s72intw(cp[2], (int64_t *)&rec_tmp);
	    cp[3] += u72int (cp[3], (uint32_t *)&r[i].mqual);
	    r[i].rec = rec_tmp;
	    if ((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
		if (!(r[i].flags & GRANGE_FLAG_TYPE_SINGLE)) {
		    int64_t pr;
		    cp[5] += s72intw(cp[5], &pr);
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
		    int64_t pr;
		    cp[5] += s72intw(cp[5], &pr);
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
#define BIN_CONS_CACHED_  (1<<8)
#define BIN_CONS_VALID_   (1<<9)

static cached_item *io_bin_read(void *dbh, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    GBin g, *b = &g;
    bin_index_t *bin;
    unsigned char *buf, *cp;
    size_t buf_len;
    uint32_t bflag;
    int version;
    int comp_mode;
    uint64_t i64;

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    buf = g_read_alloc(io, v, &buf_len);
    if (!buf && buf_len)
	return NULL;
    if (buf_len == 0) {
	/* Blank bin */
	g.pos = g.size = 0;
	g.start = g.end = 0;
	g.child[0] = g.child[1] = 0;
	g.id = 0;
	g.flags = 0;
	g.parent = g.parent_type = 0;
	g.range = 0;
	g.track = 0;
	g.nseqs = 0;
	g.rng_free = -1;
	goto empty_bin;
    }
    cp = buf;

    rdstats[GT_Bin] += buf_len;
    rdcounts[GT_Bin]++;

    assert(cp[0] == GT_Bin);
    version = cp[1];
    assert(version <= 1); /* format */
    cp += 2;
    cp += u72int(cp, &bflag);
    g.flags = (bflag & BIN_COMPLEMENTED) ? BIN_COMPLEMENTED : 0;
    if (bflag & BIN_CONS_CACHED_)
	g.flags |= BIN_CONS_CACHED;
    if (bflag & BIN_CONS_VALID_)
	g.flags |= BIN_CONS_VALID;
    g.parent_type = (bflag & BIN_ROOT_NODE) ? GT_Contig : GT_Bin;

    if (bflag & BIN_POS_ZERO)
	g.pos = 0;
    else
	cp += s72int(cp, &g.pos);

    if (bflag & BIN_SIZE_EQ_POS)
	g.size = g.pos;
    else
	cp += u72int(cp, (uint32_t *)&g.size);

    if (bflag & BIN_NO_RANGE) {
	g.range = 0;
	g.start = 0;
	g.end   = 0;
    } else {
	cp += u72int(cp, (uint32_t *)&g.start);
	cp += u72int(cp, (uint32_t *)&g.end);
	g.end += g.start;
	cp += u72intw(cp, &i64); g.range = i64;
    }

    if (bflag & BIN_NO_LCHILD) {
	g.child[0] = 0;
    } else {
	cp += u72intw(cp, &i64); g.child[0] = i64;
    }

    if (bflag & BIN_NO_RCHILD) {
	g.child[1] = 0;
    } else {
	cp += u72intw(cp, &i64); g.child[1] = i64;
    }

    if (bflag & BIN_NO_TRACK) {
	g.track = 0;
    } else {
	cp += u72intw(cp, &i64); g.track = i64;
    }

    cp += u72intw(cp, &i64); g.parent = i64;
    cp += u72int(cp, (uint32_t *)&g.nseqs);

    if (version > 0) {
	cp += s72int(cp, &g.rng_free);
    } else {
	g.rng_free = -1;
    }

 empty_bin:
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
    bin->flags       = b->flags;
    bin->parent      = b->parent;
    bin->parent_type = b->parent_type;
    bin->rng_rec     = b->range;
    bin->rng         = NULL;
    bin->track_rec   = b->track;
    bin->track       = NULL;
    bin->nseqs       = b->nseqs;
    bin->rng_free    = b->rng_free;

    /* Load ranges */
    if (b->range) {
	GViewInfo vi;
	int nranges;
	GRange *r;
	unsigned char *buf;

	v = lock(io, (int)b->range, G_LOCK_RO);
	g_view_info_(io->gdb, io->client, v, &vi);
	
	if (vi.used) {
	    unsigned int fmt;

	    buf = malloc(vi.used);
	    g_read(io, v, buf, vi.used);
	    fmt = buf[1] & 0x3f;
	    assert(buf[0] == GT_Range);
	    assert(fmt <= 1);
	    comp_mode = ((unsigned char)buf[1]) >> 6;
	    r = unpack_rng_array(comp_mode, fmt, buf+2, vi.used-2, &nranges);
	    free(buf);

	    rdstats[GT_Range] += vi.used;
	    rdcounts[GT_Range]++;

	    //printf("Unpacked %d ranges from %d bytes\n", nranges, vi.used);

	    bin->rng = ArrayCreate(sizeof(GRange), nranges);
	    if (ArrayBase(GRange, bin->rng))
		free(ArrayBase(GRange, bin->rng));
	    bin->rng->base = r;
	    ArrayRef(bin->rng, nranges-1);
	} else {
	    bin->rng = NULL;
	}

	//	g_read(io, v, ArrayBase(GRange, bin->rng),
	//	       nranges * sizeof(GRange));
	unlock(io, v);
    }

    /* Load tracks */
    if (b->track) {
	GViewInfo vi;
	size_t nitems, i;
	GBinTrack *bt;

	if (-1 == (v = lock(io, (int)b->track, G_LOCK_RO)))
	    return NULL;

	g_view_info_(io->gdb, io->client, v, &vi);
	rdstats[GT_Track] += vi.used;
	rdcounts[GT_Track]++;

	bt = (GBinTrack *)io_generic_read_i4(io, v, GT_RecArray, &nitems);
	nitems /= sizeof(GBinTrack) / sizeof(GCardinal);
	bin->track = ArrayCreate(sizeof(bin_track_t), nitems);
	bin->track->max = bin->track->dim = nitems;
	for (i = 0; i < nitems; i++) {
	    bin_track_t *t = arrp(bin_track_t, bin->track, i);
	    t->type  = bt->type;
	    t->flags = bt->flags;
	    t->rec   = bt->rec;
	    t->track = NULL; /* cached ptr for temporary tracks */
	}
	free(bt);
	unlock(io, v);
    }

    free(buf);
    return ci;
}

static int io_bin_write_view(g_io *io, bin_index_t *bin, GView v) {
    GBin g;
    int err = 0;
    uint32_t bflag;

    /* Ranges */
    if (bin->flags & BIN_RANGE_UPDATED) {
	GView v;
	char *cp, fmt[2];
	int sz;
	GIOVec vec[2];
	int fmt2 = 1;

	fmt[0] = GT_Range;
	fmt[1] = fmt2 | (io->comp_mode << 6);

	bin->flags &= ~BIN_RANGE_UPDATED;

	if (!bin->rng_rec) {
	    bin->rng_rec = allocate(io, GT_Range);
	    bin->flags |= BIN_BIN_UPDATED;
	}

	cp = pack_rng_array(io->comp_mode, fmt2, ArrayBase(GRange, bin->rng),
			    ArrayMax(bin->rng),
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

	v = lock(io, (int)bin->rng_rec, G_LOCK_EX);
	//	err |= g_write(io, v, ArrayBase(GRange, bin->rng),
	//	       sizeof(GRange) * ArrayMax(bin->rng));
	wrstats[GT_Range] += sz+2;
	wrcounts[GT_Range]++;
	vec[0].buf = fmt;   vec[0].len = 2;
	vec[1].buf = cp;    vec[1].len = sz;
	err |= g_writev(io, v, vec, 2);
	free(cp);
	err |= unlock(io, v);
    }

    /* Tracks */
    if (bin->flags & BIN_TRACK_UPDATED) {
	GView v;
	size_t nb, i, j;
	GBinTrack *bt;
	
	bin->flags &= ~BIN_TRACK_UPDATED;

	if (bin->track) {
	    bt = (GBinTrack *)malloc(sizeof(GBinTrack) * ArrayMax(bin->track));

	    /* Create GBinTrack array from non-temporary tracks */
	    for (i = j = 0; i < ArrayMax(bin->track); i++) {
		bin_track_t *from = arrp(bin_track_t, bin->track, i);

		if (from->flags & TRACK_FLAG_FREEME) /* temp. */
		    continue;

		bt[j].type  = from->type;
		bt[j].flags = from->flags;
		bt[j].rec   = (GCardinal)from->rec;
		j++;
	    }

	    /* Write them out, if we have any left */
	    if (j) {
		int o;

		if (!bin->track_rec) {
		    bin->track_rec = allocate(io, GT_Track);
		    bin->flags |= BIN_BIN_UPDATED;
		}

		v = lock(io, (int)bin->track_rec, G_LOCK_EX);

		if ((o = io_generic_write_i4(io, v, GT_RecArray, bt,
					     j * sizeof(GBinTrack))) < 0)
		    err |= 1;
		nb = o;
		err |= unlock(io, v);

		wrstats[GT_Track] += nb;
		wrcounts[GT_Track]++;
	    }

	    free(bt);
	}
    }

    /* Bin struct itself */
    if (bin->flags & BIN_BIN_UPDATED) {
	unsigned char cpstart[12*5+2], *cp = cpstart;

	bin->flags &= ~BIN_BIN_UPDATED;
	g.pos         = bin->pos;
	g.size        = bin->size;
	g.start       = bin->start_used;
	g.end         = bin->end_used;
	g.id          = 0; /* unused */
	g.flags       = bin->flags;
	g.parent      = bin->parent;
	g.parent_type = bin->parent_type;
	g.child[0]    = bin->child[0];
	g.child[1]    = bin->child[1];
	g.range       = bin->rng_rec;
	g.track       = bin->track_rec;
	g.nseqs       = bin->nseqs;
	g.rng_free    = bin->rng_free;

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
	if (g.flags & BIN_CONS_CACHED)  bflag |= BIN_CONS_CACHED_;
	if (g.flags & BIN_CONS_VALID)   bflag |= BIN_CONS_VALID_;

	*cp++ = GT_Bin;
	*cp++ = g.rng_free == -1 ? 0 : 1; /* Format */

	cp += int2u7(bflag, cp);
	if (!(bflag & BIN_POS_ZERO))     cp += int2s7(g.pos, cp);
	if (!(bflag & BIN_SIZE_EQ_POS))  cp += int2u7(g.size, cp);
	if (!(bflag & BIN_NO_RANGE)) {
	    cp += int2u7(g.start, cp);
	    cp += int2u7(g.end - g.start, cp);
	    cp += intw2u7(g.range, cp);
	}
	if (!(bflag & BIN_NO_LCHILD))    cp += intw2u7(g.child[0], cp);
	if (!(bflag & BIN_NO_RCHILD))    cp += intw2u7(g.child[1], cp);
	if (!(bflag & BIN_NO_TRACK))     cp += intw2u7(g.track, cp);
	cp += intw2u7(g.parent, cp);
	cp += int2u7(g.nseqs, cp);
	if (g.rng_free != -1)
	    cp += int2s7(g.rng_free, cp);

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
    assert(ci->rec > 0);
    return io_bin_write_view(io, bin, ci->view);
}

static tg_rec io_bin_create(void *dbh, void *vfrom) {
    //bin_index_t *from = vfrom;
    g_io *io = (g_io *)dbh;
    tg_rec rec;

    rec = allocate(io, GT_Bin);
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
	b.rng_free    = -1;
	io_bin_write_view(io, &b, v);
    }
    unlock(io, v);

    return rec;
#endif
}

static int io_bin_destroy(void *dbh, tg_rec r, GView v) {
    g_io *io = (g_io *)dbh;
    GBin g;
    unsigned char *buf, *cp;
    size_t buf_len;
    uint32_t bflag;

    /* Load from disk */
    buf = g_read_alloc(io, v, &buf_len);
    if (!buf) {
	/* Can happen when we create a new record and immediately try to
	 * destroy it, which break contig sometimes does.
	 */
	return deallocate(io, r, v);
    }

    /* Decode */
    cp = buf;

    rdstats[GT_Bin] += buf_len;
    rdcounts[GT_Bin]++;

    assert(cp[0] == GT_Bin);
    assert(cp[1] <= 1); /* format */
    cp += 2;
    cp += u72int(cp, &bflag);
    if (!(bflag & BIN_POS_ZERO))
	cp += s72int(cp, &g.pos);

    if (!(bflag & BIN_SIZE_EQ_POS))
	cp += u72int(cp, (uint32_t *)&g.size);

    if (!(bflag & BIN_NO_RANGE)) {
	GView rv;
	uint64_t i64;

	cp += u72int(cp, (uint32_t *)&g.start);
	cp += u72int(cp, (uint32_t *)&g.end);
	g.end += g.start;
	cp += u72intw(cp, &i64); g.range = i64;

	rv = lock(io, (int)g.range, G_LOCK_RW);
	deallocate(io, (int)g.range, rv);
	unlock(io, rv);
    }

    free(buf);
    return deallocate(io, r, v);
}

/* ------------------------------------------------------------------------
 * track access methods
 */

typedef struct {
    double base_qual[4];
    double gap;
    double base;
    int depth;
} cstat;

static cached_item *io_track_read(void *dbh, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    cached_item *ci;
    GView v;
    track_t *track;
    size_t buf_len, i;
    unsigned char *buf, *cp;
    uint32_t type, flags, item_size, nitems;
    
    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    buf = g_read_alloc(io, v, &buf_len);
    if (!buf && buf_len)
	return NULL;
    cp = buf;

    rdstats[GT_Track] += buf_len;
    rdcounts[GT_Track]++;

    assert(cp[0] == GT_Track);
    assert(cp[1] == 0);
    cp += 2;

    /* Decode fixed size portions */
    cp += u72int(cp, &type);
    cp += u72int(cp, &flags);
    cp += u72int(cp, &item_size);
    cp += u72int(cp, &nitems);

    /* Allocate our overlapping data objects */
    if (!(ci = cache_new(GT_Track, rec, v, NULL, sizeof(*track) + 
			 item_size * nitems)))
	return NULL;
    track = (track_t *)&ci->data;

    /* Construct track_t */
    track->rec       = rec;
    track->type      = type;
    track->flag      = flags;
    track->item_size = item_size;
    track->nitems    = nitems;
    track->data      = ArrayCreate(item_size, nitems);

    switch (type) {
    case TRACK_CONS_ARR:
	for (i = 0; i < track->nitems; i++) {
	    cstat *c = arrp(cstat, track->data, i);
	    int32_t i4;
	    uint32_t u4;

	    cp += s72int(cp, &i4); c->base_qual[0] = i4;
	    cp += s72int(cp, &i4); c->base_qual[1] = i4;
	    cp += s72int(cp, &i4); c->base_qual[2] = i4;
	    cp += s72int(cp, &i4); c->base_qual[3] = i4;
	    cp += s72int(cp, &i4); c->gap = i4;
	    cp += s72int(cp, &i4); c->base = i4;
	    cp += u72int(cp, &u4); c->depth = u4;
	}
	break;

    default:
	assert(buf_len - (cp-buf) == track->item_size * track->nitems);
	memcpy(ArrayBase(char, track->data), cp,
	       track->item_size * track->nitems);
    }

    free(buf);
    return ci;
}

static int io_track_write_view(g_io *io, track_t *track, GView v) {
    unsigned char *data, *cp;
    int err = 0, i;

    cp = data = malloc(2 + 4*5 + track->item_size * track->nitems);
    if (!data)
	return -1;

    /* Encode the fixed portions */
    *cp++ = GT_Track;
    *cp++ = 0; /* format */
    cp += int2u7(track->type, cp);
    cp += int2u7(track->flag & ~TRACK_FLAG_FREEME, cp);
    cp += int2u7(track->item_size, cp);
    cp += int2u7(track->data ? track->nitems : 0, cp);

    /* The array - try compressing this */
    switch(track->type) {
    case TRACK_CONS_ARR:
	for (i = 0; i < track->nitems; i++) {
	    cstat *c = arrp(cstat, track->data, i);
	    cp += int2s7(c->base_qual[0] + 0.5, cp);
	    cp += int2s7(c->base_qual[1] + 0.5, cp);
	    cp += int2s7(c->base_qual[2] + 0.5, cp);
	    cp += int2s7(c->base_qual[3] + 0.5, cp);
	    cp += int2s7(c->gap + 0.5, cp);
	    cp += int2s7(c->base + 0.5, cp);
	    cp += int2u7(c->depth + 0.5, cp);
	}
	break;

    default:
	if (track->nitems) {
	    memcpy(cp, ArrayBase(char, track->data),
		   track->item_size * track->nitems);
	    cp += track->item_size * track->nitems;
	}
    }
    
    wrstats[GT_Track] += cp-data;
    wrcounts[GT_Track]++;
    err |= g_write(io, v, data, cp-data);
    g_flush(io, v);
    free(data);

    return err;
}

static int io_track_write(void *dbh, cached_item *ci) {
    g_io *io = (g_io *)dbh;
    track_t *track = (track_t *)&ci->data;

    assert(ci->lock_mode >= G_LOCK_RW);
    assert(ci->rec > 0);
    return io_track_write_view(io, track, ci->view);
}

static tg_rec io_track_create(void *dbh, void *vfrom) {
    track_t *from = vfrom;
    g_io *io = (g_io *)dbh;
    tg_rec rec;
    GView v;

    rec = allocate(io, GT_Track);
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
static cached_item *seq_decode(unsigned char *buf, size_t len, tg_rec rec) {
    cached_item *ci;
    unsigned char *cp, flags, mapping_qual, seq_tech, format;
    size_t slen;
    signed int i, j;
    seq_t *seq;
    uint32_t left, right, bin, seq_len;
    uint32_t parent_rec, parent_type, bin_index;
    Array anno;

    if (len) {
	uint32_t Nanno;
	cp = buf;
	assert(cp[0] == GT_Seq);
	assert(cp[1] == 0);
	cp += 2;
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
		cp += u72int(cp, arrp(uint32_t, anno, i));
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
	/* seq name, trace name, alignment */
	cp = buf = (unsigned char *)"\0\0\0";
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

    /* SAM Aux - not supported in old single-seq mode */
    seq->aux_len = 0;
    seq->sam_aux = NULL;

    /* cp is now at the encoded seq/qual area */
    switch (seq->format) {
    case SEQ_FORMAT_MAQ:
	for (i = 0; i < seq_len; i++) {
	    seq->conf[i] = cp[i] >> 2;
	    seq->seq[i] = seq->conf[i] ? "ACGT"[cp[i] & 3] : 'N';
	}
	break;

    case SEQ_FORMAT_CNF1:
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

static cached_item *io_seq_read(void *dbh, tg_rec rec) {
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

    printf("Read rec %"PRIrec" => %p\n", rec, ci);

    if (!ci)
	return NULL;

    ci->view = v;
    ci->rec = rec;
    return ci;
}

/* See seq_decode for the storage format */
static int io_seq_write_view(g_io *io, seq_t *seq, GView v, tg_rec rec) {
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
	  2
	+ 5 /* bin rec.no */
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
    *cp++ = GT_Seq;
    *cp++ = 0; /* format */
    cp += intw2u7(seq->bin, cp);
    cp += int2u7(seq->bin_index, cp);
    cp += int2u7(seq->left, cp);
    cp += int2u7(seq->right, cp);
    cp += int2u7(seq_len, cp);

    /* Read-pair info */
    cp += intw2u7(seq->parent_rec, cp);
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
    strcpy((char *)cp, seq->name);
    cp += seq->name_len + 1;

    /* Trace name */
    strcpy((char *)cp, seq->trace_name);
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
    strcpy((char *)cp, seq->alignment);
    cp += seq->alignment_len + 1;
#endif

    /* Seq/Conf */
    switch (seq->format) {
    case SEQ_FORMAT_MAQ:
	for (i = 0; i < seq_len; i++) {
	    unsigned char qv = base2val_maq[((unsigned char *)seq->seq)[i]];
	    if (qv != 9) {
		if (seq->conf[i] <= 0)
		    qv = 0;
		else if (seq->conf[i] >= 64)
		    qv |= 63 << 2;
		else
		    qv |= seq->conf[i] << 2;
	    } else {
		qv = 0;
	    }
	    *cp++ = qv;
	}
	break;

    case SEQ_FORMAT_CNF1:
	for (i = j = 0; i < seq_len; i+=3, j++) {
	    unsigned char c1 = seq->seq[i];
	    unsigned char c2 = i+1 < seq_len ? seq->seq[i+1] : 0;
	    unsigned char c3 = i+2 < seq_len ? seq->seq[i+2] : 0;
	    *cp++ = base2val_cnf1[c1] + base2val_cnf1[c2]*6 + base2val_cnf1[c3]*36;
	}
	for (i = 0; i < seq_len; i++) {
	    *cp++ = seq->conf[i];
	}
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
    assert(ci->rec > 0);
    return io_seq_write_view(io, seq, ci->view, (GRec)ci->rec);
}

static tg_rec io_seq_create(void *dbh, void *vfrom) {
    g_io *io = (g_io *)dbh;
    return allocate(io, GT_Seq);
}

static tg_rec io_seq_index_query(void *dbh, char *name, int prefix) {
    g_io *io = (g_io *)dbh;
    
    if (!io->seq_name_tree)
	return -1;

    return btree_search(io->seq_name_tree, name, prefix);
}

static tg_rec *io_seq_index_query_all(void *dbh, char *name, int prefix,
				      int *nrecs) {
    g_io *io = (g_io *)dbh;
    
    if (!io->seq_name_tree)
	return NULL;

    return btree_search_all(io->seq_name_tree, name, prefix, nrecs);
}

static int io_seq_index_add(void *dbh, char *name, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    
    if (!io->seq_name_tree)
	return -1;

    btree_insert(io->seq_name_tree, name, rec);
    return io->seq_name_tree->root->rec;
}

static int io_seq_index_del(void *dbh, char *name) {
    g_io *io = (g_io *)dbh;
    
    if (!io->seq_name_tree)
	return -1;

    return btree_delete(io->seq_name_tree, name);
}

/* ------------------------------------------------------------------------
 * seq_block access methods
 */

#define REORDER_BY_READ_GROUP
/* NB: Do not undef now unless you remove sam_aux support too */

/*
 * Define REORDER_BY_READ_GROUP if you wish to experiment with sorting data
 * by their read group.
 *
 * The theory (and practice) is that when using mixed libraries, such as
 * many of the 1000Genomes project bam files, we have different profiles
 * for read names and quality values. It reduced the space taken up in
 * names by 17% and quality values in 5.6% - overall coming out at 6-7%.
 * This is dramatically increased if we alter the SEQ_BLOCK_BITS parameter,
 * with 8192 seqs stored in upto 1Mb chunks giving 14% savings (19% when
 * using lzma which can make better use of larger blocks). Increasing
 * SEQ_BLOCK_BITS has other, detrimental, effects though.
 *
 * Either way the reading code will handle it as the first format byte
 * is adjusted to indicate whether reordering took place.
 */

static cached_item *io_seq_block_read(void *dbh, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    GView v;
    cached_item *ci;
    seq_block_t *b;
    unsigned char *buf, *cp;
    size_t buf_len;
    seq_t in[SEQ_BLOCK_SZ];
    int i, j, k, last;
    int reorder_by_read_group = 0;
    int sam_aux = 0;
    int first_seq = 0;
    int wide_recs = 0;
    uint32_t i32;
    uint64_t i64;

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

    if (!buf_len) {
	b->est_size = 0;
	memset(&b->seq[0], 0, SEQ_BLOCK_SZ*sizeof(b->seq[0]));
	free(buf);
	return ci;
    }

    assert(buf[0] == GT_SeqBlock);
    assert((buf[1] & 0x3f) <= 7); /* format */
    if ((buf[1] & 0x3f) >= 1)
	reorder_by_read_group = 1;
    if ((buf[1] & 0x3f) & 2)
	sam_aux = 1;
    if ((buf[1] & 0x3f) & 4)
	wide_recs = 1;

    /* Ungzip it too */
    if (1) {
	size_t ssz;
	int comp_mode = ((unsigned char)buf[1]) >> 6;
	buf = (unsigned char *)mem_inflate(comp_mode, 
					   (char *)buf+2, buf_len-2, &ssz);
	free(cp);
	cp = buf;
	buf_len = ssz;
    }
    b->est_size = buf_len;

    /* Decode the fixed size components of our sequence structs */
    /* Bin */
    if (wide_recs) {
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    cp += u72intw(cp, &i64);
	    in[i].bin = i64;
	}
    } else {
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    cp += u72int(cp, &i32);
	    in[i].bin = i32;
	}
    }

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
	cp += u72int(cp, (uint32_t *)&in[i].left);
    }

    /* right clip */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, (uint32_t *)&in[i].right);
    }

    /* length */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, (uint32_t *)&in[i].len);
    }

    /* parent rec */
    if (wide_recs) {
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    if (!in[i].bin) continue;
	    cp += u72intw(cp, &i64);
	    in[i].parent_rec = i64;
	}
    } else {
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    if (!in[i].bin) continue;
	    cp += u72int(cp, &i32);
	    in[i].parent_rec = i32;
	}
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
	if (in[i].bin) {
	    first_seq = i;
	    break;
	}
    }

    if (reorder_by_read_group) {
	int p1 = 1;

	for (i = k = 0; i < SEQ_BLOCK_SZ; i++) {
	    if (!in[i].bin) continue;

	    /*
	     * If this is one we've already processed, then continue.
	     * We mark "already processed" by negating the parent_rec.
	     * For parent_rec == 0 this obviously doesn't work, so we
	     * always deal with those in the first pass.
	     */
	    if (in[i].parent_rec < 0 || (in[i].parent_rec == 0 && !p1)) {
		in[i].parent_rec = -in[i].parent_rec;
		continue;
	    }

	    cp += u72int(cp, (uint32_t *)&in[i].name_len);
	    for (j = i+1; j < SEQ_BLOCK_SZ; j++) {
		if (!in[j].bin) continue;
		if ((in[j].parent_rec && in[j].parent_rec != in[i].parent_rec)
		    || (!p1 && !in[j].parent_rec))
		    continue;
		
		cp += u72int(cp, (uint32_t *)&in[j].name_len);
		in[j].parent_rec = -in[j].parent_rec;
	    }

	    p1 = 0;
	}
    } else {
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    if (!in[i].bin) continue;
	    cp += u72int(cp, (uint32_t *)&in[i].name_len);
	}
    }

    /* trace name length */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, (uint32_t *)&in[i].trace_name_len);
    }

    /* alignment length */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, (uint32_t *)&in[i].alignment_len);
    }

    /* sam aux length */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	if (sam_aux)
	    cp += u72int(cp, (uint32_t *)&in[i].aux_len);
	else
	    in[i].aux_len = 0;
    }

    /* Convert our static structs to cached_items */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	if (in[i].bin) {
	    cached_item *si;
	    size_t extra_len;

	    extra_len =
		sizeof(seq_t) + 
		in[i].name_len + 1 +
		in[i].trace_name_len + 1 + 
		in[i].alignment_len + 1 + 
		in[i].aux_len + 1 + 
		ABS(in[i].len) +
		ABS(in[i].len) * (in[i].format == SEQ_FORMAT_CNF4 ? 4 : 1);
	    if (!(si = cache_new(GT_Seq, 0, 0, NULL, extra_len)))
		return NULL;

	    b->seq[i] = (seq_t *)&si->data;
	    in[i].rec = ((tg_rec)rec << SEQ_BLOCK_BITS) + i;
	    in[i].anno = NULL;
	    *b->seq[i] = in[i];
	    b->seq[i]->block = b;
	    b->seq[i]->idx = i;
	} else {
	    b->seq[i] = NULL;
	}
    }

    /* Decode variable sized components */
    /* Names */
    if (reorder_by_read_group) {
	for (i = k = 0; i < SEQ_BLOCK_SZ; i++) {
	    seq_t *s = b->seq[i];
	    if (!s) continue;

	    if (s->parent_rec < 0 || (i > first_seq && !s->parent_rec)) {
		s->parent_rec = -s->parent_rec;
		continue;
	    }

	    while (!b->seq[k]) k++;

	    s->name = (char *)&s->data;
	    //s->name_len = in[k++].name_len;
	    memcpy(s->name, cp, s->name_len);
	    cp += s->name_len;
	    s->name[s->name_len] = 0;

	    for (j = i+1; j < SEQ_BLOCK_SZ; j++) {
		seq_t *s2 = b->seq[j];

		if (!s2)
		    continue;

		if ((s2->parent_rec && s2->parent_rec != s->parent_rec) ||
		    (i > first_seq && !s2->parent_rec))
		    continue;

		while (!b->seq[k]) k++;

		s2->name = (char *)&s2->data;
		//s2->name_len = in[k++].name_len;
		memcpy(s2->name, cp, s2->name_len);
		cp += s2->name_len;
		s2->name[s2->name_len] = 0;

		s2->parent_rec = -s2->parent_rec;
	    }
	}
    } else {
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    if (!b->seq[i]) continue;
	    b->seq[i]->name = (char *)&b->seq[i]->data;
	    memcpy(b->seq[i]->name, cp, b->seq[i]->name_len);
	    cp += b->seq[i]->name_len;
	    b->seq[i]->name[b->seq[i]->name_len] = 0;
	}
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
    if (reorder_by_read_group) {
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    seq_t *s = b->seq[i];
	    if (!s) continue;

	    if (s->parent_rec < 0 || (i > first_seq && !s->parent_rec)) {
		s->parent_rec = -s->parent_rec;
		continue;
	    }

	    s->conf = s->seq + ABS(s->len);
	    memcpy(s->conf, cp, ABS(s->len));
	    cp += ABS(s->len);

	    for (j = i+1; j < SEQ_BLOCK_SZ; j++) {
		seq_t *s2 = b->seq[j];

		if (!s2) continue;

		if ((s2->parent_rec && s2->parent_rec != s->parent_rec) ||
		    (i > first_seq && !s2->parent_rec))
		    continue;

		s2->conf = s2->seq + ABS(s2->len);
		memcpy(s2->conf, cp, ABS(s2->len));
		cp += ABS(s2->len);

		s2->parent_rec = -s2->parent_rec;
	    }
	}
    } else {
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    if (!b->seq[i]) continue;
	    b->seq[i]->conf = b->seq[i]->seq + ABS(b->seq[i]->len);
	    memcpy(b->seq[i]->conf, cp, ABS(b->seq[i]->len));
	    cp += ABS(b->seq[i]->len);
	}
    }

    /* Sam auxillary records */
    if (sam_aux) {
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    if (!b->seq[i]) continue;
	}
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    if (!b->seq[i]) continue;
	    b->seq[i]->sam_aux = b->seq[i]->conf + 
		(b->seq[i]->format == SEQ_FORMAT_CNF4 ? 4 : 1)
		* ABS(b->seq[i]->len);
	    memcpy(b->seq[i]->sam_aux, cp, b->seq[i]->aux_len);
	    cp += b->seq[i]->aux_len;
	}
    }

    assert(cp - buf == buf_len);
    free(buf);

    return ci;
}

static int io_seq_block_write(void *dbh, cached_item *ci) {
    int err;
    g_io *io = (g_io *)dbh;
    seq_block_t *b = (seq_block_t *)&ci->data;
    int i, j, last_index;
    unsigned char *cp, *cp_start;
    unsigned char *out[19], *out_start[19], *out_malloc;
    size_t out_size[19], total_size;
    int level[19];
    GIOVec vec[2];
    char fmt[2];
    int nb = 0;
    int have_sam_aux = 0;
    int nparts = 19;
    int first_seq = -1;
    int wide_recs = sizeof(tg_rec) > sizeof(uint32_t);

    assert(ci->lock_mode >= G_LOCK_RW);
    assert(ci->rec > 0);

    set_dna_lookup();

    /* Compute worst-case sizes, for memory allocation */
    for (i = 0; i < 19; i++) {
	out_size[i] = 0;
    }
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	seq_t *s = b->seq[i];
	if (!s) {
	    out_size[0]++;
	    continue;
	}

	out_size[0] += 10;/* bin */
	out_size[1] += 5; /* bin_index */
	out_size[2] += 5; /* left */
	out_size[3] += 5; /* right */
	out_size[4] += 5; /* len */
	out_size[5] += 10;/* parent_rec */
	out_size[6] ++;   /* parent type */
	out_size[7] ++;   /* format */
	out_size[8] ++;   /* m.qual */
	out_size[9] += 5; /* name len */
	out_size[10]+= 5; /* trace name len */
	out_size[11]+= 5; /* alignment len */
	out_size[12]+= s->name_len+1;
	out_size[13]+= s->trace_name_len+1;
	out_size[14]+= s->alignment_len;
	out_size[15]+= ABS(s->len); /* seq */
	out_size[16]+= ABS(s->len) * (s->format == SEQ_FORMAT_CNF4 ? 4 : 1);
	out_size[17]+= 5; /* aux_len - moved before [12] later */
	out_size[18]+= s->aux_len;
	nb += ABS(s->len);
    }
    for (total_size = i = 0; i < 19; i++)
	total_size += out_size[i]+1;
    out_malloc = malloc(total_size);
    for (total_size = i = 0; i < 19; i++) {
	out_start[i] = out[i] = &out_malloc[total_size];
	total_size += out_size[i]+1;
    }

    /* serialised sequences, separated by type of data */
    last_index = 0;
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	int delta;
	seq_t *s = b->seq[i];

	if (!s) {
	    *out[0]++= 0; /* signifies sequence not present */
	    continue;
	}

	if (first_seq == -1)
	    first_seq = i;

	if (wide_recs) {
	    out[0] += intw2u7(s->bin, out[0]);
	} else {
	    out[0] += int2u7((int32_t)s->bin, out[0]);
	}
	delta = s->bin_index - last_index;
	last_index = s->bin_index;
	out[1] += int2s7(delta, out[1]);
	out[2] += int2u7(s->left, out[2]);
	out[3] += int2u7(s->right, out[3]);
	out[4] += int2u7(ABS(s->len), out[4]);

	if (wide_recs) {
	    out[5] += intw2u7(s->parent_rec, out[5]);
	} else {
	    out[5] += int2u7((int32_t)s->parent_rec, out[5]);
	}
	*out[6]++ = s->parent_type;
	
	/* flags & m.quality */
	s->format = SEQ_FORMAT_CNF1;
	*out[7]++ = (s->format << 6) | (s->flags << 3) | s->seq_tech;

	/* Duplicated in range, but adds about 1% on test bam input */
	*out[8]++ = s->mapping_qual;
	
#ifndef REORDER_BY_READ_GROUP
	/* Name */
	out[9] += int2u7(s->name_len, out[9]);
	memcpy(out[12], s->name, s->name_len);
	out[12] += s->name_len;
#endif

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

	/* Sam auxillary records */
	out[17] += int2u7(s->aux_len, out[17]);
	if (s->aux_len) {
	    memcpy(out[18], s->sam_aux, s->aux_len);
	    out[18] += s->aux_len;
	    have_sam_aux = 1;
	}

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
#ifndef REORDER_BY_READ_GROUP
	memcpy(out[16], s->conf, ABS(s->len)); out[16] += ABS(s->len);
#endif
    }

#if 0
    /* rotated quality */
    /* Generally this isn't an improvement */
    unsigned char *tmp = out[16];
    for (j = 0; nb; j++) {
	for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	    int delta;
	    seq_t *s = b->seq[i];

	    if (!s) {
		continue;
	    }

	    if (j > ABS(s->len))
		continue;

	    *out[16]++ = s->conf[j];
	    nb--;
	}
    }
#endif

#ifdef REORDER_BY_READ_GROUP
    /* quality, name sorted by library */
    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	seq_t *s = b->seq[i];

	if (!s) {
	    continue;
	}

	/*
	 * First time through, with i==0, we permit s->parent_rec == 0 to
	 * be written out (mixed in with the true parent_rec of seq[0]).
	 *
	 * Subsequent times we only right out +ve parent_rec values.
	 *
	 * FIXME: Use qsort instead, as per sam_comp2.c. It's faster
	 * and less confusing to understand, although a *different*
	 * format due to the order changing.
	 */
	if (s->parent_rec < 0 || (i > first_seq && !s->parent_rec)) {
	    s->parent_rec = -s->parent_rec;
	    continue;
	}

	memcpy(out[16], s->conf, ABS(s->len)); out[16] += ABS(s->len);

	out[9] += int2u7(s->name_len, out[9]);
	memcpy(out[12], s->name, s->name_len);
	out[12] += s->name_len;

	for (j = i+1; j < SEQ_BLOCK_SZ; j++) {
	    seq_t *s2 = b->seq[j];

	    if (!s2)
		continue;

	    if ((s2->parent_rec && s2->parent_rec != s->parent_rec) ||
		(i > first_seq && !s2->parent_rec))
		continue;

	    memcpy(out[16], s2->conf, ABS(s2->len)); out[16] += ABS(s2->len);

	    out[9] += int2u7(s2->name_len, out[9]);
	    memcpy(out[12], s2->name, s2->name_len);
	    out[12] += s2->name_len;

	    s2->parent_rec = -s2->parent_rec;
	}
    }
#endif

    /* Concatenate data types together and adjust out_size to actual usage */
    if (have_sam_aux) {
	/* Reorder as we need aux_len before the variable sized portions */
	int oz = out_size[17];
	unsigned char *o = out[17], *os = out_start[17];
	for (i = 16; i >= 12; i--) {
	    out[i+1] = out[i];
	    out_start[i+1] = out_start[i];
	    out_size[i+1] = out_size[i];
	}
	out[12] = o;
	out_size[12] = oz;
	out_start[12] = os;
    } else {
	nparts = 17;
    }
    for (total_size = i = 0; i < nparts; i++) {
	out_size[i] = out[i] - out_start[i];
	total_size += out_size[i];
    }

    cp = cp_start = malloc(total_size+1);
    for (i = 0; i < nparts; i++) {
	memcpy(cp, out_start[i], out_size[i]);
	cp += out_size[i];
	level[i] = 5; /* default gzip level */
    }
    assert(cp - cp_start == total_size);

    level[12] = 7; /* name, 8 is approx 1% better, but 10% slower */
    level[13] = 7; /* trace name */
    level[16] = 6; /* conf */
    level[18] = 7;

    /* NB: should name and trace name be compressed together? They're
     * probably the same or highly related?
     */

    //printf("Block %d, est.size %d, actual size %d, diff=%d %f gzip=",
    //       ci->rec, b->est_size, cp-cp_start,
    //       cp-cp_start - b->est_size, (cp-cp_start + 0.0) / b->est_size);
    
    /* Gzip it too */
    if (1) {
	unsigned char *gzout;
	size_t ssz;

	//gzout = mem_deflate(cp_start, cp-cp_start, &ssz);
	gzout = (unsigned char *)mem_deflate_lparts(io->comp_mode,
						    (char *)cp_start,
						    out_size, level,
						    nparts, &ssz);
	free(cp_start);
	cp_start = gzout;
	cp = cp_start + ssz;
    }

    /* Finally write the serialised data block */
    /* Format:
     * bit 0     0 => orig format (1 => reorder by RG, also any other bit)
     * bit 1     0 => no sam_aux, 1 => have them
     * bit 2     0 => 32-bit rec, 1 => 64-bit rec
     * bit 3-5   reserved (0)
     * bit 6-7   Compression method
     *
     * NB: any bit 0-5 set also implies reorder by RG. Ie format 2 originally
     * meant reorder + sam_aux (and so still does).
     */
    fmt[0] = GT_SeqBlock;
    fmt[1] = (have_sam_aux ? 2 : 1) | (io->comp_mode << 6); /* format */
    if (wide_recs)
	fmt[1] |= (1<<2);
    vec[0].buf = fmt;      vec[0].len = 2;
    vec[1].buf = cp_start; vec[1].len = cp - cp_start;
    
    assert(ci->lock_mode >= G_LOCK_RW);
    wrstats[GT_SeqBlock] += cp-cp_start + 2;
    wrcounts[GT_SeqBlock]++;
    err = g_writev(io, ci->view, vec, 2);
    g_flush(io, ci->view);
   
    free(cp_start);
    free(out_malloc);

    return err ? -1 : 0;
}

static tg_rec io_seq_block_create(void *dbh, void *vfrom) {
    g_io *io = (g_io *)dbh;
    tg_rec rec;
    GView v;

    rec = allocate(io, GT_SeqBlock);
    v = lock(io, rec, G_LOCK_EX);

    /* Write blank data here? */

    unlock(io, v);
    
    return rec;
}


/* ------------------------------------------------------------------------
 * seq_block access methods
 */
static cached_item *io_anno_ele_block_read(void *dbh, tg_rec rec) {
    g_io *io = (g_io *)dbh;
    GView v;
    cached_item *ci;
    anno_ele_block_t *b;
    unsigned char *buf, *cp;
    size_t buf_len;
    anno_ele_t in[ANNO_ELE_BLOCK_SZ];
    int i;
    uint64_t last, i64;
    int comment_len[ANNO_ELE_BLOCK_SZ];

    /* Load from disk */
    if (-1 == (v = lock(io, rec, G_LOCK_RO)))
	return NULL;

    if (!(ci = cache_new(GT_AnnoEleBlock, rec, v, NULL, sizeof(*b))))
	return NULL;

    b = (anno_ele_block_t *)&ci->data;
    cp = buf = (unsigned char *)g_read_alloc((g_io *)dbh, v, &buf_len);

    if (!buf_len) {
	b->est_size = 0;
	memset(&b->ae[0],  0, ANNO_ELE_BLOCK_SZ*sizeof(b->ae[0]));
	free(buf);
	return ci;
    }

    assert(buf[0] == GT_AnnoEleBlock);
    assert((buf[1] & 0x3f) == 0); /* Format */

    rdstats[GT_AnnoEleBlock] += buf_len;
    rdcounts[GT_AnnoEleBlock]++;

    /* Ungzip it too */
    if (1) {
	size_t ssz;
	int comp_mode = ((unsigned char)buf[1]) >> 6;
	buf = (unsigned char *)mem_inflate(comp_mode,
					   (char *)buf+2, buf_len-2, &ssz);
	free(cp);
	cp = buf;
	buf_len = ssz;
    }
    b->est_size = buf_len;

    /* Decode the fixed size components of our sequence structs */
    /* Bin */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	cp += u72intw(cp, &i64);
	in[i].bin = i64;
    }

    /* Tag type */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, (uint32_t *)&in[i].tag_type);
    }

    /* Obj type */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, (uint32_t *)&in[i].obj_type);
    }
    
    /* Obj record */
    for (last = i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	int64_t tmp;
	if (!in[i].bin) continue;
	cp += s72intw(cp, &tmp);
	in[i].obj_rec = last + tmp;
	last = in[i].obj_rec;
	//cp += u72int(cp, (uint32_t *)&in[i].obj_rec);
    }

    /* Anno record */
    for (last = i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	int64_t tmp;
	if (!in[i].bin) continue;
	cp += s72intw(cp, &tmp);
	in[i].anno_rec = last + tmp;
	tmp = in[i].anno_rec;
	//cp += u72int(cp, (uint32_t *)&in[i].anno_rec);
    }

    /* Comment length */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	if (!in[i].bin) continue;
	cp += u72int(cp, (uint32_t *)&comment_len[i]);
    }


    /* Convert our static structs to cached_items */
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	if (in[i].bin) {
	    cached_item *si;
	    size_t extra_len;

	    extra_len = sizeof(anno_ele_t) + comment_len[i];
	    if (!(si = cache_new(GT_AnnoEle, 0, 0, NULL, extra_len)))
		return NULL;

	    b->ae[i]  = (anno_ele_t *)&si->data;
	    in[i].rec = ((tg_rec)rec << ANNO_ELE_BLOCK_BITS) + i;
	    *b->ae[i] = in[i];
	    b->ae[i]->block = b;
	    b->ae[i]->idx = i;
	    b->ae[i]->comment = (char *)&b->ae[i]->data;
	} else {
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
    int i;
    tg_rec last_obj_rec, last_anno_rec;
    unsigned char *cp, *cp_start;
    unsigned char *out[7], *out_start[7];
    size_t out_size[7], total_size;
    int level[7];
    GIOVec vec[2];
    char fmt[2];

    assert(ci->lock_mode >= G_LOCK_RW);
    assert(ci->rec > 0);

    /* Compute worst-case sizes, for memory allocation */
    for (i = 0; i < 7; i++) {
	out_size[i] = 0;
    }
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	anno_ele_t *e = b->ae[i];
	if (!e) {
	    out_size[0]++;
	    continue;
	}

	out_size[0] += 10;/* bin */
	out_size[1] += 5; /* tag type */
	out_size[2] += 5; /* obj type */
	out_size[3] += 10;/* obj record */
	out_size[4] += 10;/* anno record */
	out_size[5] += 5; /* comment length */
	out_size[6] += e->comment ? strlen(e->comment) : 0; /* comments */
    }
    for (i = 0; i < 7; i++)
	out_start[i] = out[i] = malloc(out_size[i]+1);
    
    /* serialised annotations, separated by type of data */
    last_obj_rec = last_anno_rec = 0;
    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	tg_rec delta;
	anno_ele_t *e = b->ae[i];

	if (!e) {
	    *out[0]++=0; /* signifies annotation not present */
	    continue;
	}

	out[0] += intw2u7(e->bin, out[0]);
	out[1] += int2u7(e->tag_type, out[1]);
	out[2] += int2u7(e->obj_type, out[2]);

	delta = e->obj_rec - last_obj_rec;
	last_obj_rec = e->obj_rec;
	out[3] += intw2s7(delta, out[3]);
	//out[3] += int2u7(e->obj_rec, out[3]);

	delta = e->anno_rec - last_anno_rec;
	last_anno_rec = e->anno_rec;
	out[4] += intw2s7(delta, out[4]);
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
    cp = cp_start = malloc(total_size+1);
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
	unsigned char *gzout;
	size_t ssz;

	//gzout = mem_deflate(cp_start, cp-cp_start, &ssz);
	gzout = (unsigned char *)mem_deflate_lparts(io->comp_mode,
						    (char *)cp_start,
						    out_size, level, 7, &ssz);
	free(cp_start);
	cp_start = gzout;
	cp = cp_start + ssz;
    }

    /* Finally write the serialised data block */
    fmt[0] = GT_AnnoEleBlock;
    fmt[1] = 0 | (io->comp_mode << 6); /* format */
    vec[0].buf = fmt;      vec[0].len = 2;
    vec[1].buf = cp_start; vec[1].len = cp - cp_start;

    assert(ci->lock_mode >= G_LOCK_RW);
    wrstats[GT_AnnoEleBlock] += cp-cp_start + 2;
    wrcounts[GT_AnnoEleBlock]++;
    err = g_writev(io, ci->view, vec, 2);
    g_flush(io, ci->view);
    
    free(cp_start);

    return err ? -1 : 0;
}

static tg_rec io_anno_ele_block_create(void *dbh, void *vfrom) {
    g_io *io = (g_io *)dbh;
    tg_rec rec;
    GView v;

    rec = allocate(io, GT_AnnoEleBlock);
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
    io_database_setopt,
    io_rec_exists,

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
	io_database_write,
	io_generic_info,
	io_database_create_index,
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
	io_contig_index_del,
    },

    {
	/* Bin */
	io_bin_create,
	io_bin_destroy,
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
	io_seq_index_query_all,
	io_seq_index_add,
	io_seq_index_del,
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
