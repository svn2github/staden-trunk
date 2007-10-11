#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/hash_table.h>

#include "srf.h"

/*
 * ---------------------------------------------------------------------------
 * Object creation / destruction.
 */

/*
 * Allocates and returns an srf_t structure.
 * fp is the file point associated with this SRF file. Pass as NULL if unknown
 * at this stage.
 *
 * Returns malloced memory on success
 *         NULL on failure
 */
srf_t *srf_create(FILE *fp) {
    srf_t *srf = (srf_t *)calloc(1, sizeof(*srf));

    if (srf)
	srf->fp = fp;
    return srf;
}

/*
 * Opens an SRF archive for reading or writing as defined by 'mode'. Mode is
 * passed directly on to fopen() and so uses the same flags.
 *
 * Returns a srf_t struct pointer on success.
 *         NULL on failure
 */
srf_t *srf_open(char *fn, char *mode) {
    FILE *fp;

    return (fp = fopen(fn, mode)) ? srf_create(fp) : NULL;
}

/*
 * Deallocates an srf_t struct. If auto_close is true then it also closes
 * any associated FILE pointer.
 */
void srf_destroy(srf_t *srf, int auto_close) {
    if (!srf)
	return;
    
    if (srf->th.trace_hdr) {
	free(srf->th.trace_hdr);
	srf->th.trace_hdr = NULL;
    }

    if (auto_close && srf->fp)
	fclose(srf->fp);

    free(srf);
}

/*
 * ---------------------------------------------------------------------------
 * Base data type I/O.
 */

/* 
 * Writes a null-terminated C string in pascal-string form.
 * Returns the number of bytes written or
 *         -1 for failure.
 */
int srf_write_pstring(srf_t *srf, char *str) {
    size_t l = str ? strlen(str) : 0;
    if (l > 255)
	return -1;

    if (l)
	return fprintf(srf->fp, "%c%s", (int)l, str);
    else
	return fprintf(srf->fp, "%c", (int)l);
}

/*
 * Reads a pascal-style string from the srf file.
 * 'str' passed in needs to be at least 256 bytes long. The string read
 * will be stored there and nul terminated.
 *
 * Returns the length of the string read (minus nul char),
 *         -1 for failure
 */
int srf_read_pstring(srf_t *srf, char *str) {
    int len;

    if (EOF == (len = fgetc(srf->fp)))
	return -1;
    if (len != fread(str, 1, len, srf->fp))
	return -1;
    str[len] = '\0';

    return len;
}

/*
 * Read/write unsigned 32-bit and 64-bit values in big-endian format
 * (ie the same as ZTR and SCF endian uses).
 * Functions all return 0 for success, -1 for failure.
 */
int srf_read_uint32(srf_t *srf, uint32_t *val) {
    unsigned char d[4];
    if (1 != fread(d, 4, 1, srf->fp))
	return -1;

    *val = (d[0] << 24) | (d[1] << 16) | (d[2] << 8) | (d[3] << 0);
    return 0;
}

int srf_write_uint32(srf_t *srf, uint32_t val) {
    unsigned char d[4];
    d[0] = (val >> 24) & 0xff;
    d[1] = (val >> 16) & 0xff;
    d[2] = (val >>  8) & 0xff;
    d[3] = (val >>  0) & 0xff;

    return fwrite(d, 4, 1, srf->fp) ? 0 : -1;
}

int srf_read_uint64(srf_t *srf, uint64_t *val) {
    unsigned char d[8];
    if (1 != fread(d, 8, 1, srf->fp))
	return -1;

    *val = ((uint64_t)d[0] << 56)
	 | ((uint64_t)d[1] << 48)
	 | ((uint64_t)d[2] << 40)
	 | ((uint64_t)d[3] << 32)
	 | ((uint64_t)d[4] << 24)
	 | ((uint64_t)d[5] << 16)
	 | ((uint64_t)d[6] <<  8)
	 | ((uint64_t)d[7] <<  0);
    return 0;
}

int srf_write_uint64(srf_t *srf, uint64_t val) {
    unsigned char d[8];
    d[0] = (val >> 56) & 0xff;
    d[1] = (val >> 48) & 0xff;
    d[2] = (val >> 40) & 0xff;
    d[3] = (val >> 32) & 0xff;
    d[4] = (val >> 24) & 0xff;
    d[5] = (val >> 16) & 0xff;
    d[6] = (val >>  8) & 0xff;
    d[7] = (val >>  0) & 0xff;

    return fwrite(d, 8, 1, srf->fp) ? 0 : -1;
}

/*
 * ---------------------------------------------------------------------------
 * Mid level I/O - srf block type handling
 */


/*
 * Allocates and initialises a container header structure. An existing
 * srf_cont_hdr_t may be passed in to avoid allocation of a new object.
 *
 * Returns: allocated cont_header on success, to be freed using
 *          srf_destroy_cont_hdr() (if allocated here),
 *          NULL on failure.
 */
srf_cont_hdr_t *srf_construct_cont_hdr(srf_cont_hdr_t *ch,
				       char *bc,
				       char *bc_version) {
    if (!ch) {
	if (NULL == (ch = (srf_cont_hdr_t *)calloc(1, sizeof(*ch))))
	    return NULL;
    }

    ch->next_block_offset = 0;
    ch->block_type = SRFB_CONTAINER;
    strcpy(ch->version, SRF_VERSION);
    ch->container_type = 'Z';
    strncpy(ch->base_caller, bc, 255);
    strncpy(ch->base_caller_version, bc_version, 255);

    return ch;
}

/*
 * Deallocates an srf_cont_hdr_t constructed by srf_construct_cont_hdr().
 */
void srf_destroy_cont_hdr(srf_cont_hdr_t *ch) {
    if (ch)
	free(ch);
}

/*
 * Reads a container header and stores the result in 'ch'.
 * Returns 0 for success
 *        -1 for failure
 */
int srf_read_cont_hdr(srf_t *srf, srf_cont_hdr_t *ch) {
    char magic[3];

    if (!ch)
	return -1;

    /* Check block type */
    if (EOF == (ch->block_type = fgetc(srf->fp)))
	return -1;
    if (ch->block_type != SRFB_CONTAINER)
	return -1;

    /* Check magic number && version */
    if (3 != fread(magic, 1, 3, srf->fp))
	return -1;
    if (srf_read_pstring(srf, ch->version) < 0)
	return -1;
    if (strncmp(magic, "SRF", 3) || strcmp(ch->version, SRF_VERSION))
	return -1;
    
    /* Containter type, base caller bits */
    if (EOF == (ch->container_type = fgetc(srf->fp)) ||
	srf_read_pstring(srf, ch->base_caller) < 0||
	srf_read_pstring(srf, ch->base_caller_version) < 0)
	return -1;

    /* Offsets */
    if (0 != srf_read_uint32(srf, &ch->next_block_offset))
	return -1;

    return 0;
}

/*
 * Writes a container header to disk.
 *
 * 4 Block type + magic number ("SSRF")
 * 1+n: pString for version number ("1.0")
 * 1: "Z" => container type is ZTR
 * 1+n: pString for base-caller ("eg Bustard")
 * 1+n: pString for base-caller version (eg "1.8.28")
 * 4: uint32 distance from start of header to first data block.
 * 4: uint32 distance from start of block to index block (0). FIXME
 *
 * Returns 0 for success
 *        -1 for failure
 */
int srf_write_cont_hdr(srf_t *srf, srf_cont_hdr_t *ch) {
    uint32_t sz = 0;
    int32_t z;

    if (!ch)
	return -1;

    /* Magic number && version */
    if (4 != fwrite(SRF_MAGIC, 1, 4, srf->fp))
	return -1;
    sz += 4;
    if ((z = srf_write_pstring(srf, ch->version)) < 0)
	return -1;
    sz += z;
    
    /* Containter type, base caller bits */
    if (EOF == fputc(ch->container_type, srf->fp))
	return -1;
    sz++;
    if ((z = srf_write_pstring(srf, ch->base_caller)) < 0)
	return -1;
    sz += z;
    if ((z = srf_write_pstring(srf, ch->base_caller_version)) < 0)
	return -1;
    sz += z;

    /* No XML data - so sz stays as-is for now */
    sz += 0;

    /* Offsets */
    sz += 8;
    ch->next_block_offset = sz;
    if (0 != srf_write_uint32(srf, ch->next_block_offset))
	return -1;

    return 0;
}

/*
 * Initialises a srf_trace_header_t and inserts some passed in values.
 * If the supplied th is NULL then a new structure is allocated.
 *
 * Returns a pointer to the dh passed in (or allocated if NULL) on success
 *         NULL on failure
 */
srf_trace_hdr_t *srf_construct_trace_hdr(srf_trace_hdr_t *th,
					 char *prefix,
					 unsigned char *header,
					 uint32_t header_sz) {
    if (!th) {
	if (NULL == (th = (srf_trace_hdr_t *)calloc(1, sizeof(*th))))
	    return NULL;
    }

    th->block_type = SRFB_TRACE_HEADER;
    strncpy(th->id_prefix, prefix, 255);
    th->trace_hdr_size = header_sz;
    th->trace_hdr = header;

    return th;
}

/*
 * Deallocates a srf_trace_hdr_t if allocated by
 * srf_construct_trace_hdr().
 * Do not use this if you passed in a static srf_trace_hdr to the construct
 * function.
 */
void srf_destroy_trace_hdr(srf_trace_hdr_t *th) {
    if (th)
	free(th);
}

/*
 * Reads a data header and stores the result in 'th'.
 * Returns 0 for success
 *        -1 for failure
 */
int srf_read_trace_hdr(srf_t *srf, srf_trace_hdr_t *th) {
    /* Check block type */
    if (EOF == (th->block_type = fgetc(srf->fp)))
	return -1;
    if (th->block_type != SRFB_TRACE_HEADER)
	return -1;

    /* Read-id prefix */
    if (srf_read_pstring(srf, th->id_prefix) < 0)
	return -1;

    /* The data header itself */
    if (0 != srf_read_uint32(srf, &th->trace_hdr_size))
	return -1;
    if (th->trace_hdr_size) {
	if (NULL == (th->trace_hdr = malloc(th->trace_hdr_size)))
	    return -1;
	if (th->trace_hdr_size != fread(th->trace_hdr, 1,
					  th->trace_hdr_size, srf->fp)) {
	    free(th->trace_hdr);
	    return -1;
	}
    } else {
	th->trace_hdr = NULL;
    }

    return 0;
}

/*
 * Writes a srf_trace_hdr_t structure to disk.
 * Returns 0 for sucess
 *        -1 for failure
 */
int srf_write_trace_hdr(srf_t *srf, srf_trace_hdr_t *th) {
    if (!srf->fp)
	return -1;

    if (EOF == fputc(th->block_type, srf->fp))
	return -1;

    if (-1 == srf_write_pstring(srf, th->id_prefix))
	return -1;

    /* The ztr header blob itself... */
    if (-1 == srf_write_uint32(srf, th->trace_hdr_size))
	return -1;
    if (th->trace_hdr_size !=
	fwrite(th->trace_hdr, 1, th->trace_hdr_size, srf->fp))
	return -1;

    return 0;
}


srf_trace_body_t *srf_construct_trace_body(srf_trace_body_t *tb,
					   char *suffix,
					   unsigned char *body,
					   uint32_t body_size) {
    if (!tb) {
	if (NULL == (tb = (srf_trace_body_t *)calloc(1, sizeof(*tb))))
	    return NULL;
    }
    tb->block_type = SRFB_TRACE_BODY;
    strncpy(tb->read_id, suffix, 255);
    tb->trace = body;
    tb->trace_size = body_size;
    tb->flags = 0;

    return tb;
}

void srf_destroy_trace_body(srf_trace_body_t *tb) {
    if (tb)
	free(tb);
}

/*
 * Writes a new trace body.
 *
 * Returns: 0 on success
 *          -1 on failure
 */
int srf_write_trace_body(srf_t *srf, srf_trace_body_t *tb) {
    if (!srf->fp)
	return -1;

    if (EOF == fputc(tb->block_type, srf->fp))
	return -1;
    if (-1 == srf_write_pstring(srf, tb->read_id))
	return -1;

    if (EOF == (fputc(tb->flags, srf->fp)))
	return -1;

    /* Tbe ztr footer blob itself... */
    if (0 != srf_write_uint32(srf, tb->trace_size))
	return -1;
    if (tb->trace_size != fwrite(tb->trace, 1, tb->trace_size, srf->fp))
	return -1;

    return 0;
}

/*
 * Reads a trace header + trace 'blob' and stores the result in 'th'
 * If no_trace is true then it skips loading the trace data itself.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int srf_read_trace_body(srf_t *srf, srf_trace_body_t *tb, int no_trace) {
    int c;

    /* Check block type */
    if (EOF == (tb->block_type = fgetc(srf->fp)))
	return -1;
    if (tb->block_type != SRFB_TRACE_BODY)
	return -1;

    /* Read-id suffix */
    if (srf_read_pstring(srf, tb->read_id) < 0)
	return -1;

    if (EOF == (c = fgetc(srf->fp)))
	return -1;
    else
	tb->flags = c;

    /* The trace data itself */
    if (0 != srf_read_uint32(srf, &tb->trace_size))
	return -1;
    if (!no_trace) {
	if (tb->trace_size) {
	    if (NULL == (tb->trace = malloc(tb->trace_size)))
		return -1;
	    if (tb->trace_size != fread(tb->trace, 1, tb->trace_size,
					srf->fp)) {
		free(tb->trace);
		return -1;
	    }
	} else {
	    tb->trace = NULL;
	}
    } else {
	tb->trace = NULL;
    }

    return 0;
}


/*
 * Reads a SRF index header. See srf_write_index_hdr for the format.
 *
 * Returns 0 on success and fills out *hdr
 *         -1 on failure
 */
int srf_read_index_hdr(srf_t *srf, srf_index_hdr_t *hdr) {
    /* Load footer */
    if (0 != fseeko(srf->fp, -16, SEEK_END))
	return -1;

    if (4 != fread(hdr->magic,   1, 4, srf->fp))
	return -1;
    if (4 != fread(hdr->version, 1, 4, srf->fp))
	return -1;
    if (0 != srf_read_uint64(srf, &hdr->size))
	return -1;

    /* Check for validity */
    if (memcmp(hdr->magic,   SRF_INDEX_MAGIC,   4) ||
	memcmp(hdr->version, SRF_INDEX_VERSION, 4))
	return -1;

    /* Seek to index header and re-read */
    if (0 != fseeko(srf->fp, -hdr->size, SEEK_END))
	return -1;
    
    if (4 != fread(hdr->magic,   1, 4, srf->fp))
	return -1;
    if (4 != fread(hdr->version, 1, 4, srf->fp))
	return -1;
    if (0 != srf_read_uint64(srf, &hdr->size))
	return -1;
    if (0 != srf_read_uint32(srf, &hdr->n_container))
	return -1;
    if (0 != srf_read_uint32(srf, &hdr->n_data_block_hdr))
	return -1;
    if (0 != srf_read_uint32(srf, &hdr->n_buckets))
	return -1;
    if (1 != fread(&hdr->hash_func, 1, 1, srf->fp))
	return -1;

    /* Check once more */
    if (memcmp(hdr->magic,   SRF_INDEX_MAGIC,   4) ||
	memcmp(hdr->version, SRF_INDEX_VERSION, 4))
	return -1;

    return 0;
}

/*
 * Writes a SRF index header.
 *
 * Header:
 *   x4    magic number, starting with 'I'.
 *   x4    version code (eg "1.00")
 *   x8    index size
 *   x4    number of containers
 *   x4    number of DBHs
 *   x4    number of hash buckets (~10 billion traces per file is enough).
 *   x1    hash function (see hash_table.h)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int srf_write_index_hdr(srf_t *srf, srf_index_hdr_t *hdr) {
    if (4 != fwrite(hdr->magic,   1, 4, srf->fp))
	return -1;
    if (4 != fwrite(hdr->version, 1, 4, srf->fp))
	return -1;
    if (0 != srf_write_uint64(srf, hdr->size))
	return -1;
    if (0 != srf_write_uint32(srf, hdr->n_container))
	return -1;
    if (0 != srf_write_uint32(srf, hdr->n_data_block_hdr))
	return -1;
    if (0 != srf_write_uint32(srf, hdr->n_buckets))
	return -1;
    if (1 != fwrite(&hdr->hash_func, 1, 1, srf->fp))
	return -1;

    return 0;
}


/*
 * ---------------------------------------------------------------------------
 * Higher level I/O functions
 */

/*
 * Fetches the next trace from an SRF container.
 * Name, if defined (which should be a buffer of at least 512 bytes long)
 * will be filled out to contain the read name.
 * 
 * Returns mFILE containing trace on success
 *         NULL on failure.
 */
mFILE *srf_next_trace(srf_t *srf, char *name) {
    do {
	int type;

	switch(type = srf_next_block_type(srf)) {
	case -1:
	    /* EOF */
	    return NULL;

	case SRFB_CONTAINER:
	    if (0 != srf_read_cont_hdr(srf, &srf->ch))
		return NULL;
	    break;

	case SRFB_TRACE_HEADER:
	    if (srf->th.trace_hdr) {
		free(srf->th.trace_hdr);
		srf->th.trace_hdr = NULL;
	    }

	    if (0 != srf_read_trace_hdr(srf, &srf->th))
		return NULL;

	    break;

	case SRFB_TRACE_BODY: {
	    mFILE *mf = mfcreate(NULL, 0);
	    srf_trace_body_t tb;

	    if (!mf || 0 != srf_read_trace_body(srf, &tb, 0))
		return NULL;

	    if (name)
		sprintf(name, "%s%s", srf->th.id_prefix, tb.read_id);

	    if (srf->th.trace_hdr_size)
		mfwrite(srf->th.trace_hdr, 1, srf->th.trace_hdr_size, mf);
	    if (tb.trace_size)
		mfwrite(tb.trace, 1, tb.trace_size, mf);
	    if (tb.trace)
		free(tb.trace);

	    mrewind(mf);
	    return mf;
	}

	case SRFB_INDEX: {
	    off_t pos = ftell(srf->fp);
	    srf_index_hdr_t hdr;
	    srf_read_index_hdr(srf, &hdr);

	    /* Skip the index body */
	    fseeko(srf->fp, pos + hdr.size, SEEK_SET);
	    break;
	}

	default:
	    fprintf(stderr, "Block of unknown type '%c'. Aborting\n", type);
	    return NULL;
	}
    } while (1);

    return NULL;
}

/*
 * Returns the type of the next block.
 * -1 for none (EOF)
 */
int srf_next_block_type(srf_t *srf) {
    int c = fgetc(srf->fp);
    if (c == EOF)
	return -1;

    ungetc(c, srf->fp);

    return c;
}

/*
 * Reads the next SRF block from an archive and returns the block type.
 * If the block is a trace it'll return the full trace name too (maximum
 * 512 bytes).
 *
 * Returns block type on success, writing to pos and name as appropriate
 *         -1 on EOF
 *         -2 on failure
 */
int srf_next_block_details(srf_t *srf, uint64_t *pos, char *name) {
    int type;
    *pos = ftell(srf->fp);

    switch(type = srf_next_block_type(srf)) {
    case -1:
	/* EOF */
	return -1;

    case SRFB_CONTAINER:
	if (0 != srf_read_cont_hdr(srf, &srf->ch))
	    return -2;
	break;

    case SRFB_TRACE_HEADER:
	if (0 != srf_read_trace_hdr(srf, &srf->th))
	    return -2;
	
	break;

    case SRFB_TRACE_BODY: {
	srf_trace_body_t tb;

	/* Inefficient, but it'll do for testing purposes */
	if (0 != srf_read_trace_body(srf, &tb, 1))
	    return -2;

	if (name)
	    sprintf(name, "%s%s", srf->th.id_prefix, tb.read_id);
	
	break;
    }

    case SRFB_INDEX: {
	srf_index_hdr_t hdr;
	srf_read_index_hdr(srf, &hdr);

	/* Skip the index body */
	fseeko(srf->fp, *pos + hdr.size, SEEK_SET);
	break;
    }


    default:
	fprintf(stderr, "Block of unknown type '%c'. Aborting\n", type);
	return -2;
    }

    return type;
}

/*
 * Searches through 'nitems' 8-byte values stored in 'srf' at file offset
 * 'start' onwards for the closest value <= 'query'.
 *
 * Returns 0 on success, setting *res
 *        -1 on failure
 */
static int binary_scan(srf_t *srf, int nitems, uint64_t start, uint64_t query,
		       uint64_t *res) {
    int min = 0;
    int max = nitems;
    int guess, i;
    uint64_t pos = 0, best = 0;

    if (nitems <= 0)
	return -1;

    /* Binary search on disk for approx location */
    while (max - min > 100) {
	guess = (max - min) / 2 + min;

	if (guess == max)
	    guess = max-1;

	if (-1 == fseeko(srf->fp, guess * 8 + start, SEEK_SET))
	    return -1;
	if (0 != srf_read_uint64(srf, &pos))
	    return -1;
	if (pos > query) {
	    max = guess;
	} else {
	    min = guess;
	}
    }

    /* Within a small distance => linear scan now to avoid needless disk IO */
    if (-1 == fseeko(srf->fp, min * 8 + start, SEEK_SET))
	return -1;
    for (i = min; i < max; i++) {
	if (0 != srf_read_uint64(srf, &pos))
	    return -1;
	if (pos > query) {
	    break;
	} else {
	    best = pos;
	}
    }

    assert(best <= query);
    *res = best;

    return 0;
}

/*
 * Searches in an SRF index for a trace of a given name.
 * If found it sets the file offsets for the container (cpos), data block
 * header (hpos) and data block (dpos).
 *
 * On a test with 2 containers and 12 headers this averaged at 6.1 reads per
 * trace fetch and 8.0 seeks.
 *
 * Returns 0 on success
 *        -1 on failure (eg no index)
 *        -2 on trace not found in index.
 */
int srf_find_trace(srf_t *srf, char *tname,
		   uint64_t *cpos, uint64_t *hpos, uint64_t *dpos) {
    srf_index_hdr_t hdr;
    uint64_t hval, bnum;
    uint32_t bucket_pos;
    off_t ipos, skip;

    /* Check for valid index */
    if (0 != srf_read_index_hdr(srf, &hdr)) {
	return -1;
    }
    ipos = ftello(srf->fp);
    skip = hdr.n_container * 8 + hdr.n_data_block_hdr * 8;

    /* Hash and load the bucket */
    hval = hash64(hdr.hash_func, (unsigned char *)tname, strlen(tname));
    bnum = hval & (hdr.n_buckets - 1);
    if (-1 == fseeko(srf->fp, ipos + skip + bnum * 4, SEEK_SET))
	return -1;

    if (0 != srf_read_uint32(srf, &bucket_pos))
	return -1;
    if (!bucket_pos)
	return -2;

    /* Secondary hash is the top 7-bits */
    hval >>= 57;

    /* Jump to the item list */
    if (-1 == fseeko(srf->fp, ipos-SRF_INDEX_HDR_SIZE + bucket_pos, SEEK_SET))
	return -1;
    for (;;) {
	int h = fgetc(srf->fp);
	off_t saved_pos;
	
	if ((h & 0x7f) != hval) {
	    if (h & 0x80)
		return -2; /* end of list and not found */
	    /*
	     * fseeko(srf->fp, 8, SEEK_CUR);
	     * Use fread instead as it's likely already cached and linux
	     * fseeko involves a real system call (lseek).
	     */
	    fread(dpos, 1, 8, srf->fp);
	    continue;
	}

	/* Potential hit - investigate to see if it's the real one: */
	/* Seek to dpos and get trace id suffix. Compare to see if valid */
	if (0 != srf_read_uint64(srf, dpos))
	    return -1;
	saved_pos = ftello(srf->fp);
	if (-1 == fseeko(srf->fp, (off_t)*dpos, SEEK_SET))
	    return -1;
	if (0 != srf_read_trace_body(srf, &srf->tb, 0))
	    return -1;
	if (strcmp(tname + strlen(tname) - strlen(srf->tb.read_id),
		   srf->tb.read_id)) {
	    if (-1 == fseeko(srf->fp, saved_pos, SEEK_SET))
		return -1;
	    continue;
	}

	/* Binary search the index to identify the matching cpos & hpos */
	if (0 != binary_scan(srf, hdr.n_container,
			     ipos, *dpos, cpos))
	    return -1;

	if (0 != binary_scan(srf, hdr.n_data_block_hdr,
			     ipos + hdr.n_container * 8,
			     *dpos, hpos))
	    return -1;

	/* Check trace id prefix matches */
	if (-1 == fseeko(srf->fp, *hpos, SEEK_SET))
	    return -1;
	if (0 != srf_read_trace_hdr(srf, &srf->th))
	    return -1;
	if (strncmp(tname, srf->th.id_prefix, strlen(srf->th.id_prefix))) {
	    if (-1 == fseeko(srf->fp, saved_pos, SEEK_SET))
		return -1;
	    continue;
	}

	/* FIXME: what to do with base-caller and cpos */

	/* Found it! */
	break;
    }

    return 0;
}
