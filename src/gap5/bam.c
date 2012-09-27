/*
 * Handles reading (only for now) from a BAM format file.
 * We use zlib to read the data in and fetch data one read at a time.
 */

#include <staden_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <fcntl.h>
#include <zlib.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>

#include "bam.h"
#include "os.h"

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif


static int bam_more_input(bam_file_t *b);
static int bam_more_output(bam_file_t *b);
static int reg2bin(int start, int end);
static int bgzf_write(int fd, const void *buf, size_t count);

/*
 * Allow for unaligned memory access. This is used in cigar string processing
 * as the packed BAM data struct has cigar after read name instead of before.
 */
#define ALLOW_UAC

/*
 * Reads len bytes from fp into data.
 *
 * Returns the number of bytes read.
 *         0 for eof.
 *        -1 for failure.
 */
static int bam_read(bam_file_t *b, void *data, size_t len) {
    int nb = 0, n;
    unsigned char *cdata = data;

    while (len) {
	/* Consume any available uncompressed output */
	if (b->out_sz) {
	    size_t l = MIN(b->out_sz, len);
	    memcpy(cdata, b->out_p, l);
	    b->out_p += l;
	    b->out_sz -= l;
	    cdata += l;
	    len -= l;
	    nb += l;

	    if (!len)
		return nb;
	}

	if (!b->gzip) {
	    /* Already uncompressed, so easy to deal with */
	    if (!b->in_sz)
		if (-1 == bam_more_input(b))
		    return nb ? nb : 0;
		    
	    b->out_p  = b->in_p;
	    b->out_sz = b->in_sz;
	    b->in_sz  = 0;
	    continue;
	}

	n = bam_more_output(b);
	if (n == -1)
	    return -1;
	if (n == 0)
	    return nb;
    }

    return nb;
}

/*
 * Reads a line of text of unknown length.
 * 'str' is both input and output. If *str == NULL then memory is allocated
 * for the line. If *str != NULL then it is expected to point to an existing
 * block of memory that we can write into and realloc as required.
 *
 * Similarly *len is both input and output. It is expected to hold the
 * current allocated size of *str. It is modified if we realloc it.
 *
 * Lines have the \n removed and will be null terminated.
 *
 * Returns actual line length used (note note the same as *len) on success
 *        -1 on failure
 */
int bam_get_line(bam_file_t *b, unsigned char **str, size_t *len) {
    unsigned char *buf = *str;
    int used_l = 0;
    size_t alloc_l = *len;
    int next_condition;

    while (b->out_sz || bam_more_output(b) > 0) {
	int tmp;
	unsigned char *from = b->out_p;
	unsigned char *to   = &buf[used_l];

	/*
	 * Next condition is the number of loop iterations before something
	 * has to be done - either getting more uncompressed output or
	 * resizing the buffer. We don't care which, but it allows us to
	 * have just one check per loop instead of two. Once out of the loop
	 * we can then afford to determine which case is and deal with it.
	 */
	tmp = next_condition = MIN(b->out_sz, alloc_l-used_l);

	while (next_condition-- > 0) { /* these 3 lines are 50% of SAM cpu */
	    if (*from != '\n') {
		*to++ = *from++;
	    } else {
		b->out_p = from;
		used_l = to-buf;
		b->out_p++;
		buf[used_l] = 0;
		*str = buf;
		*len = alloc_l;
		b->out_sz -= (tmp - next_condition);
		return used_l;
	    }
	}

	used_l = to-buf;
	b->out_p = from;
	b->out_sz -= tmp;

	if (used_l >= alloc_l) {
	    alloc_l = alloc_l ? alloc_l * 2 : 1024;
	    if (NULL == (buf = realloc(buf, alloc_l)))
		return -1;
	}
    }

    return b->out_sz ? -1 : 0;
}

int load_bam_header(bam_file_t *b) {
    char magic[4];
    int i;

    if (4 != bam_read(b, magic, 4))
	return -1;
    if (memcmp(magic, "BAM\x01",4) != 0)
	return -1;
    if (4 != bam_read(b, &b->header_len, 4))
	return -1;
    b->header_len = le_int4(b->header_len);
    if (!(b->header = malloc(b->header_len)))
	return -1;
    if (b->header_len != bam_read(b, b->header, b->header_len))
	return -1;
    if (4 != bam_read(b, &b->nref, 4))
	return -1;
    b->nref = le_int4(b->nref);

    b->ref = malloc(b->nref * sizeof(*b->ref));
    for (i = 0; i < b->nref; i++) {
	uint32_t nlen;
	if (4 != bam_read(b, &nlen, 4))
	    return -1;
	nlen = le_int4(nlen);

	b->ref[i].name = calloc(1, nlen);
	if (nlen != bam_read(b, b->ref[i].name, nlen))
	    return -1;
	if (4 != bam_read(b, &b->ref[i].len, 4))
	    return -1;
	b->ref[i].len = le_int4(b->ref[i].len);
    }

    for (i = 0; i < b->header_len; i++) {
	if (b->header[i] == '\n')
	    b->line++;
    }

    return 0;
}

int load_sam_header(bam_file_t *b) {
    unsigned char *str = NULL;
    size_t alloc = 0, len;
    int header_pos = 0;

    b->header = NULL;
    b->header_len = 0;
    
    b->ref  = NULL;
    b->nref = 0;

    while ((b->out_sz > 0 || bam_more_output(b) > 0) && *b->out_p == '@') {
	b->line++;
	if ((len = bam_get_line(b, &str, &alloc)) == -1)
	    return -1;

	/* Add header lines to b->header */
	if (header_pos + len + 1>= b->header_len) {
	    b->header_len += 8192;
	    b->header = realloc(b->header, b->header_len);
	    if (!b->header)
		return -1;
	}
	memcpy(&b->header[header_pos], str, len);
	b->header[header_pos+len] = '\n';
	header_pos += len+1;

	/* Also parse any reference lines while reading the header */
	if (str[1] == 'S' && str[2] == 'Q') {
	    int rlen = -1;
	    char *rname = NULL;
	    unsigned char *cp = str+3;
	    HashData hd;

	    //printf("line=%ld/%ld/'%s'\n", (long)len, (long)strlen(str), str);

	    if (*cp != '\t') {
		fprintf(stderr, "Malformed header line '%s'\n", str);
		return -1;
	    }
	    cp++;

	    /* Tokenise line into key/value pairs */
	    while (*cp) {
		char *key = (char *)cp;
		char *val;

		if (!key[0] || !key[1] || key[2] != ':') {
		    fprintf(stderr, "Malformed header line '%s'\n", str);
		    return -1;
		}
		key[2] = '\0';

		cp = (unsigned char *)(val = key+3);
		while (*cp && *cp != '\t')
		    cp++;

		if (*cp)
		    *cp++ = 0;

		if (0 == strcmp(key, "LN")) {
		    rlen = atoi(val);
		} else if (0 == strcmp(key, "SN")) {
		    rname = val;
		}
	    }

	    if (!rname || rlen == -1) {
		fprintf(stderr, "No SN or LN value in @SQ line\n");
		return -1;
	    }

	    b->ref = realloc(b->ref, (b->nref+1)*sizeof(*b->ref));
	    if (!b->ref)
		return -1;
	    b->ref[b->nref].len  = rlen;
	    b->ref[b->nref].name = strdup(rname);

	    hd.i = b->nref;
	    HashTableAdd(b->ref_hash, b->ref[b->nref].name, 0, hd, NULL);

	    b->nref++;
	}
    }

    /* Blank header if none supplied */
    if (!b->header)
	b->header = strdup("");

    if (header_pos)
	b->header[header_pos] = '\0';
    b->header_len = header_pos;

    return 0;
}

/*
 * Extracts @RG records from the header and places them in a hash table.
 * Returns 0 on success
 *        -1 on failure
 */
int bam_parse_header(bam_file_t *b) {
    int i, j, lno = 0;
    int ntabs, ref_len;
    tag_list_t *tags;

    if (!b->header)
	return -1;

    if (!b->rg_hash)
	b->rg_hash = HashTableCreate(4, HASH_FUNC_HSIEH |
				     HASH_DYNAMIC_SIZE |
				     HASH_NONVOLATILE_KEYS);
    if (!b->ref_hash)
	b->ref_hash = HashTableCreate(16, HASH_FUNC_HSIEH |
				      HASH_DYNAMIC_SIZE |
				      HASH_NONVOLATILE_KEYS);

    if (!b->rg_hash || !b->ref_hash)
	return -1;

    /* Extract the data into RG and SQ hashes */
    b->nref = 0;
    for (i = 0; i < b->header_len; i++) {
	char *id = NULL;
	int id_len;
	HashData hd;
	int i_start = i, i_len;
	int is_rg;

	for (j = i; b->header[j] && j < b->header_len; j++)
	    if (b->header[j] == '\n')
		break;
	i_len = j - i_start;
	lno++;

	if (b->header[i] != '@') {
	    fprintf(stderr, "Header line does not start with '@' at line %d:\n"
		    "\"%.*s\"\n",
		    lno, i_len, &b->header[i_start]);
            i = i_start + i_len;
            continue;
	}

	/* RG or SQ only */
	if (!((b->header[i+1] == 'R' && b->header[i+2] == 'G') ||
	      (b->header[i+1] == 'S' && b->header[i+2] == 'Q'))) {
	    i = i_start + i_len;
	    continue;
	}

	is_rg = (b->header[i+1] == 'R' && b->header[i+2] == 'G');

	/* Tokenise, not fast but doesn't matter */
	if (is_rg) {
	    for (ntabs = 0, j = i; j < b->header_len && b->header[j] != '\n'; j++)
		if (b->header[j] == '\t')
		    ntabs++;

	    if (ntabs == 0) {
		fprintf(stderr, "Missing tab in header line %d:\n\"%.*s\"\n",
			lno, i_len, &b->header[i_start]);
		i = i_start + i_len;
		continue;
	    }

	    if (NULL == (tags = malloc((ntabs+1) * sizeof(*tags))))
		return -1;

	    ntabs = 0;
	} else {
	    ref_len = 0;
	}
	
	while (i < b->header_len && b->header[i] != '\n') {
	    if (b->header[i] == '\t')
		break;
	    i++;
	}

	while (i < b->header_len && b->header[i] != '\n') {
	    if (b->header[i] == '\t') {
		int key;
		char *value;

		if (!b->header[i+1] || !b->header[i+2] ||
		    b->header[i+3] != ':') {
		    fprintf(stderr,
			    "Malformed key:value pair in header line %d:\n"
			    "\"%.*s\"\n",
			    lno, i_len, &b->header[i_start]);
		    i = i_start + i_len;
		    continue;
		}

		key = (b->header[i+1]<<8) + b->header[i+2]; 
		i+=4;
		
		value = &b->header[i];
		while (i < b->header_len && b->header[i] != '\n') {
		    if (b->header[i] == '\t')
			break;
		    i++;
		}

		if (is_rg) {
		    if (key == (('I'<<8) | 'D')) {
			id = value;
			id_len = &b->header[i] - value;
		    } else {
			tags[ntabs].value  = value;
			tags[ntabs].key    = key;
			tags[ntabs].length = &b->header[i] - value;
			ntabs++;
		    }
		} else {
		    if (key == (('S'<<8) | 'N')) {
			id = value;
			id_len = &b->header[i] - value;
		    }

		    if (key == (('L'<<8) | 'N')) {
			ref_len = atoi(value);
		    }
		}
	    }
	}
	if (is_rg) {
	    tags[ntabs].value = NULL;
	    
	    if (id) {
		int new = 0;
		hd.p = tags;
		HashTableAdd(b->rg_hash, id, id_len, hd, &new);
		if (!new)
		    free(tags);
	    } else {
		fprintf(stderr, "No ID record in @RG line\n");
	    }
	} else {
	    if (id) {
		char tmp;

		hd.i = b->nref++;
		HashTableAdd(b->ref_hash, id, id_len, hd, NULL);
		b->ref = realloc(b->ref, b->nref * sizeof(*b->ref));
		tmp = id[id_len]; id[id_len] = 0;
		b->ref[b->nref-1].name = strdup(id);
		id[id_len] = tmp;
		b->ref[b->nref-1].len = ref_len;
	    } else {
		fprintf(stderr, "No SN record in @SQ line\n");
	    }
	}

	if (0 && is_rg) {
	    printf("RG ID=%.*s\n", id_len, id);
	    tag_list_t *t = tags;
	    while (t->value) {
		printf("%c%c -> '%.*s'\n",
		       t->key >> 8, t->key & 0xff,
		       t->length,
		       t->value);
		t++;
	    }
	}
    }

    return 0;
}

tag_list_t *bam_find_rg(bam_file_t *b, char *id) {
    HashItem *hi = HashTableSearch(b->rg_hash, id, 0);
    return hi ? (tag_list_t *)hi->data.p : NULL;
}

#ifndef O_BINARY
#    define O_BINARY 0
#endif

bam_file_t *bam_open(char *fn, char *mode) {
    bam_file_t *b = calloc(1, sizeof *b);
    
    if (!b)
	return NULL;

    b->in_p     = b->in;
    b->in_sz    = 0;
    b->out_p    = b->out;
    b->out_sz   = 0;
    b->next_len = -1;
    b->bs       = NULL;
    b->bs_size  = 0;
    b->z_finish = 1;
    b->bgzf     = 0;
    b->no_aux   = 0;
    b->line     = 0;
    b->binary   = 0;

    /* Creation */
    if (*mode == 'w') {
	b->mode = O_WRONLY | O_TRUNC | O_CREAT;
	if (mode[1] == 'b') {
	    b->mode |= O_BINARY;
	    b->binary = 1;
	}

	if (strcmp(fn, "-") == 0) {
	    b->fd = 1; /* Stdout */
	} else {
	    if (-1 == (b->fd = open(fn, b->mode, 0666)))
		goto error;
	}

	return b;
    }

    if (*mode != 'r')
	return NULL;

    if (strcmp(mode, "rb") == 0) {
	b->mode = O_RDONLY | O_BINARY;
    } else {
	b->mode = O_RDONLY;
    }
    if (-1 == (b->fd = open(fn, b->mode, 0)))
	goto error;

    b->ref_hash = HashTableCreate(16, HASH_FUNC_HSIEH |
				      HASH_DYNAMIC_SIZE |
				      HASH_NONVOLATILE_KEYS);

    /* Load first block so we can check */
    bam_more_input(b);
    if (b->in_sz < 2)
	return NULL;
    if (b->in_p[0] == 31 && b->in_p[1] == 139)
	b->gzip = 1;
    else
	b->gzip = 0;

    if (b->gzip) {
	/* Set up zlib */
	b->s.zalloc    = NULL;
	b->s.zfree     = NULL;
	b->s.opaque    = NULL;
	inflateInit2(&b->s, -15);
    }

    /* Load header */
    if (strcmp(mode, "rb") == 0) {
	if (-1 == load_bam_header(b))
	    goto error;
	b->bam = 1;
    } else {
	if (-1 == load_sam_header(b))
	    goto error;
	b->bam = 0;
    }

    /* Parse RG & SQ records in header */
    if (-1 == bam_parse_header(b))
	return NULL;

    return b;

 error:
    if (b) {
	if (b->header)
	    free(b->header);
	free(b);
    }

    return NULL;
}

void bam_close(bam_file_t *b) {
    if (b->mode & O_WRONLY) {
	if (b->binary) {
	    if (bgzf_write(b->fd, b->out, b->out_p - b->out)) {
		fprintf(stderr, "Write failed in bam_close()\n");
	    }

	    /* Output a blank BGZF block too to mark EOF */
	    if (28 != write(b->fd, "\037\213\010\4\0\0\0\0\0\377\6\0\102\103"
			    "\2\0\033\0\3\0\0\0\0\0\0\0\0\0", 28)) {
		fprintf(stderr, "Write failed in bam_close()\n");
	    }
	} else {
	    if (b->out_p - b->out != write(b->fd, b->out, b->out_p - b->out)) {
		fprintf(stderr, "Write failed in bam_close()\n");
	    }
	}
    }

    if (!b)
	return;

    if (b->bs)
	free(b->bs);
    if (b->header)
	free(b->header);
    if (b->ref) {
	int i;
	for (i = 0; i < b->nref; i++)
	    free(b->ref[i].name);
	free(b->ref);
    }

    if (b->gzip)
	inflateEnd(&b->s);

    if (b->ref_hash)
	HashTableDestroy(b->ref_hash, 0);

    if (b->rg_hash)
	HashTableDestroy(b->rg_hash, 1);

    close(b->fd);

    free(b);
}

/*
 * Loads more data into the input (compressed) buffer.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
static int bam_more_input(bam_file_t *b) {
    size_t l;

    if (b->in != b->in_p) {
	memmove(b->in, b->in_p, b->in_sz);
	b->in_p = b->in;
    }

    l = read(b->fd, &b->in[b->in_sz], Z_BUFF_SIZE - b->in_sz);
    if (l <= 0)
	return -1;
    
    b->in_sz += l;
    return 0;
}

/*
 * Converts compressed input to the uncompressed output buffer
 *
 * Returns number of additional output bytes on success
 *         0 on eof
 *        -1 on failure.
 */
static int bam_more_output(bam_file_t *b) {
    int err = Z_OK;
    unsigned char *bgzf;
    int xlen, bsize;

    assert(b->out_sz == 0);

    if (!b->gzip) {
	/* Already uncompressed, so easy to deal with */
	if (!b->in_sz)
	    if (-1 == bam_more_input(b))
		return 0;
		    
	b->out_p  = b->in_p;
	b->out_sz = b->in_sz;
	b->in_sz  = 0;
	return b->out_sz;
    }

    /* Uncompress another BGZF block */
    /* BGZF header */
    if (b->in_sz < 18) {
	if (-1 == bam_more_input(b))
	    return 0;
	if (b->in_sz < 18)
	    return -1;
	if (b->in_sz == 0)
	    return 0; /* eof */
    }
	
    if (b->z_finish) {
	/*
	 * BGZF header is gzip + extra fields.
	 */
	bgzf = b->in_p;
	b->in_p += 10; b->in_sz -= 10;

	if (bgzf[0] != 31 || bgzf[1] != 139)
	    return -1; /* magic number failure */
	if ((bgzf[3] & 4) == 4) {
	    /* has extra fields, eg BGZF */
	    xlen = bgzf[10] + bgzf[11]*256;
	    b->in_p += 2; b->in_sz -= 2;
	} else {
	    xlen = 0;
	}
    } else {
	/* Continuing with an existing data stream */
	xlen = 0;
    }

    /* BGZF */
    if (xlen == 6) {
	b->bgzf = 1;
	b->in_p += 6; b->in_sz -= 6;

	if (bgzf[12] != 'B' || bgzf[13] != 'C' ||
	    bgzf[14] !=  2  || bgzf[15] !=  0)
	    return -1;
	bsize = bgzf[16] + bgzf[17]*256;
	bsize -= 6+19;

	/* Inflate */
	if (b->in_sz < bsize + 8) {
	    do {
		if (bam_more_input(b) == -1)
		    return -1; /* Truncated */
	    } while (b->in_sz < bsize + 8);
	}
	b->s.avail_in  = bsize;
	b->s.next_in   = b->in_p;
	b->s.avail_out = Z_BUFF_SIZE;
	b->s.next_out  = b->out;
	b->s.total_out = 0;
	    
	inflateReset(&b->s);
	err = inflate(&b->s, Z_FINISH);
	if (err != Z_STREAM_END) {
	    fprintf(stderr, "Inflate returned error code %d\n", err);
	    return -1;
	}
	b->z_finish = 1;

	b->in_p   += bsize + 8; /* crc & isize */
	b->in_sz  -= bsize + 8;
	b->out_sz  = b->s.total_out;
	b->out_p   = b->out;
    } else {
	/* Some other gzip variant, but possibly still having xlen */
	while (xlen) {
	    int d = MIN(b->in_sz, xlen);
	    xlen     -= d;
	    b->in_p  += d;
	    b->in_sz -= d;
	    if (b->in_sz == 0)
		bam_more_input(b);
	    if (b->in_sz == 0)
		return -1; /* truncated file */
	}
	    
	b->s.avail_in  = b->in_sz;
	b->s.next_in   = b->in_p;
	b->s.avail_out = Z_BUFF_SIZE;
	b->s.next_out  = b->out;
	b->s.total_out = 0;
	if (b->z_finish)
	    inflateReset(&b->s);
	    
	do {
	    err = inflate(&b->s, Z_BLOCK);
	    //printf("err=%d\n", err);

	    if (err == Z_OK || err == Z_STREAM_END) {
		b->in_p  += b->in_sz - b->s.avail_in;
		b->in_sz  = b->s.avail_in;
		b->out_sz = b->s.total_out;
		b->out_p  = b->out;
	    }

	    if (err == Z_STREAM_END) {
		b->z_finish = 1;
		/* Consume (ignore) CRC & ISIZE */
		if (b->in_sz < 8)
		    bam_more_input(b);

		if (b->in_sz < 8)
		    return -1; /* truncated file */

		b->in_sz -= 8;
		b->in_p  += 8;
	    } else {
		b->z_finish = 0;
	    }
	} while (err != Z_OK && err != Z_STREAM_END);
    }

    return b->out_sz;
}

/*
 * Decodes the next line of SAM into a bam_seq_t struct.
 */
int sam_next_seq(bam_file_t *b, bam_seq_t **bsp) {
    static unsigned char *str = NULL;
    static size_t alloc_l = 0;
    int used_l, n, sign;
    unsigned char *cpf, *cpt, *cp;
    int cigar_len;
    bam_seq_t *bs;
    HashItem *hi;
    int start, end;
    static int lookup[256] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 00 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 10 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 20 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 30 */
	0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0, /* 40 */
	0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0, /* 50 */
	0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0, /* 60 */
	0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0, /* 70 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 80 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 90 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* a0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* b0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* c0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* d0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* e0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};/* f0 */

    /* Fetch a single line */
    if ((used_l = bam_get_line(b, &str, &alloc_l)) <= 0) {
	return used_l;
    }

    used_l *= 4; // FIXME

    /* Over sized memory, for worst case? FIXME: cigar can break this! */
    if (!*bsp || used_l + sizeof(*bs) > (*bsp)->alloc) {
	if (!(*bsp = realloc(*bsp, used_l + sizeof(*bs))))
	    return -1;
	(*bsp)->alloc = used_l + sizeof(*bs);
	(*bsp)->blk_size = 0; /* compute later */
    }

    bs = *bsp;
    bs->flag_nc = 0;
    bs->bin_mq_nl = 0;
    
    /* Decode line */
    cpf = str;
    cpt = (unsigned char *)&bs->data;
    
    /* Name */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	*cpt++ = *cpf++;
    *cpt++ = 0;
    if (!*cpf++) return -1;
    bam_set_name_len(bs, cpf-cp);

    /* flag */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bam_set_flag(bs, n);

    /* ref */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	cpf++;
    if (*cp == '*') {
	/* Unmapped */
	bs->ref = -1;
    } else {
	hi = HashTableSearch(b->ref_hash, (char *)cp, cpf-cp);
	if (!hi) {
	    HashData hd;

	    fprintf(stderr, "Reference seq %.*s unknown\n", (int)(cpf-cp), cp);

	    /* Fabricate it instead */
	    b->ref = realloc(b->ref, (b->nref+1)*sizeof(*b->ref));
	    if (!b->ref)
		return -1;
	    b->ref[b->nref].len  = 0; /* Unknown value */
	    b->ref[b->nref].name = malloc(cpf-cp+1);
	    memcpy(b->ref[b->nref].name, cp, cpf-cp);
	    b->ref[b->nref].name[cpf-cp] = 0;

	    hd.i = b->nref;
	    hi = HashTableAdd(b->ref_hash, b->ref[b->nref].name, 0, hd, NULL);
	    b->nref++;
	}
	bs->ref = hi->data.i;
    }
    cpf++;

    /* Pos */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bs->pos = n-1;
    start = end = bs->pos;

    /* map qual */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bam_set_map_qual(bs, n);

    /* cigar */
    n = 0;
    cigar_len = 0;
    if (*cpf == '*') {
	cpf++;
    } else {
	//while (*cpf && *cpf != '\t') {
	while (*cpf > '\t') {
	    if (isdigit(*cpf)) {
		n = n*10 + *cpf++ - '0';
	    } else {
		unsigned char op;
		union {
		    unsigned char c[4];
		    uint32_t i;
		} c4i;

		switch (*cpf++) {
		case 'M': op=BAM_CMATCH;         end+=n; break;
		case 'I': op=BAM_CINS;                   break;
		case 'D': op=BAM_CDEL;           end+=n; break;
		case 'N': op=BAM_CREF_SKIP;      end+=n; break;
		case 'S': op=BAM_CSOFT_CLIP;             break;
		case 'H': op=BAM_CHARD_CLIP;             break;
		case 'P': op=BAM_CPAD;                   break;
		case '=': op=BAM_CBASE_MATCH;    end+=n; break;
		case 'X': op=BAM_CBASE_MISMATCH; end+=n; break;
		default:
		    fprintf(stderr, "Unknown cigar opcode '%c'\n", cpf[-1]);
		    return -1;
		}

		c4i.i = (n << 4) | op;
		*cpt++ = c4i.c[0];
		*cpt++ = c4i.c[1];
		*cpt++ = c4i.c[2];
		*cpt++ = c4i.c[3];

		n = 0;
		cigar_len++;
	    }
	}
    }
    bam_set_cigar_len(bs, cigar_len);
    //printf("pos %d, %d..%d => bin %d\n", bs->pos, start, end, reg2bin(start, end));
    bam_set_bin(bs, reg2bin(start,end));
    if (!*cpf++) return -1;
    
    /* mate ref name */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	cpf++;
    if (*cp == '*' && cp[1] == '\t') {
	bs->mate_ref = -1;
    } else if (*cp == '=' && cp[1] == '\t') {
	bs->mate_ref = bs->ref;
    } else {
	hi = HashTableSearch(b->ref_hash, (char *)cp, cpf-cp);
	if (!hi) {
	    HashData hd;

	    fprintf(stderr, "Mate ref seq \"%.*s\" unknown\n", (int)(cpf-cp), cp);

	    /* Fabricate it instead */
	    b->ref = realloc(b->ref, (b->nref+1)*sizeof(*b->ref));
	    if (!b->ref)
		return -1;
	    b->ref[b->nref].len  = 0; /* Unknown value */
	    b->ref[b->nref].name = malloc(cpf-cp+1);
	    memcpy(b->ref[b->nref].name, cp, cpf-cp);
	    b->ref[b->nref].name[cpf-cp] = 0;

	    hd.i = b->nref;
	    hi = HashTableAdd(b->ref_hash, b->ref[b->nref].name, 0, hd, NULL);
	    b->nref++;
	}
	bs->mate_ref = hi->data.i;
    }
    if (!*cpf++) return -1;

    /* mate pos */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bs->mate_pos = n-1;

    /* insert size */
    n = 0;
    if (*cpf == '-') {
	sign = -1;
	cpf++;
    } else {
	sign = 1;
    }
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bs->ins_size = n*sign;

    /* seq */
    cp = cpf;
    n = 0;
    //while (*cpf && *cpf != '\t') {
    if (cpf[0] == '*' && cpf[1] == '\t') {
	cpf++;
	bs->len = 0;
    } else {
	while (*cpf > '\t') {
	    if (n == 0) {
		*cpt = lookup[*cpf]<<4;
		n = 1;
	    } else {
		n = 0;
		*cpt++ |= lookup[*cpf];
	    }
	    cpf++;
	}
	if (n == 1)
	    cpt++;
	bs->len = cpf-cp;
    }
    if (!*cpf++) return -1;

    /* qual */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    if (cpf[0] == '*' && (cpf[1] == '\0' || cpf[1] == '\t')) {
	/* no qual */
	memset(cpt, '\xFF', bs->len);
	cpt += bs->len;
	cpf++;
    } else {
	while (*cpf > '\t')
	    *cpt++ = *cpf++ - '!';
    }

    assert((char *)cpt == (char *)(bam_aux(bs)));

    if (!*cpf++ || b->no_aux) goto skip_aux;

    /* aux */
    while (*cpf) {
	unsigned char *key = cpf, *value;
	if (!(key[0] && key[1] && key[2] == ':' && key[3] && key[4] == ':'))
	    return -1;
	cpf += 5;

	value = cpf;
	while (*cpf && *cpf != '\t')
	    cpf++;

	*cpt++ = key[0];
	*cpt++ = key[1];

	switch(key[3]) {
	case 'A':
	    *cpt++ = 'A';
	    *cpt++ = *value;
	    break;

	case 'i':
	    n = atoi((char *)value);
	    if (n >= 0) {
		if (n < 256) {
		    *cpt++ = 'C';
		    *cpt++ = n;
		} else if (n < 65536) {
		    *cpt++ = 'S';
		    *cpt++ = n & 0xff;
		    *cpt++ = (n >> 8) & 0xff;
		} else {
		    *cpt++ = 'I';
		    *cpt++ = n & 0xff;
		    *cpt++ = (n >> 8) & 0xff;
		    *cpt++ = (n >>16) & 0xff;
		    *cpt++ = (n >>24) & 0xff;
		}
	    } else {
		if (n >= -128 && n < 128) {
		    *cpt++ = 'c';
		    *cpt++ = n;
		} else if (n >= -32768 && n < 32768) {
		    *cpt++ = 's';
		    *cpt++ = n & 0xff;
		    *cpt++ = (n >> 8) & 0xff;
		} else {
		    *cpt++ = 'i';
		    *cpt++ = n & 0xff;
		    *cpt++ = (n >> 8) & 0xff;
		    *cpt++ = (n >>16) & 0xff;
		    *cpt++ = (n >>24) & 0xff;
		}
	    }
	    break;

	case 'f': {
	    union {
		float f;
		int i;
	    } u;
	    u.f = atof((char *)value);
	    *cpt++ = 'f';
	    *cpt++ = (u.i) & 0xff;
	    *cpt++ = (u.i >> 8) & 0xff;
	    *cpt++ = (u.i >>16) & 0xff;
	    *cpt++ = (u.i >>24) & 0xff;
	    break;
	}

	case 'Z':
	    *cpt++ = 'Z';
	    while (value != cpf)
		*cpt++=*value++;
	    *cpt++ = 0;
	    break;

	case 'H':
	    *cpt++ = 'H';
	    while (value != cpf)
		*cpt++=*value++;
	    *cpt++ = 0;
	    break;

	case 'B': {
	    char subtype = *value++;
	    unsigned char *sz;
	    int count = 0;

	    *cpt++ = 'B';
	    *cpt++ = subtype;
	    sz = cpt; cpt += 4; /* Fill out later */

	    while (*value == ',') {
		value++;
		switch (subtype) {
		case 'c': case 'C':
		    *cpt++ = strtol((char *)value, (char **)&value, 10);
		    break;

		case 's': case 'S':
		    n = strtol((char *)value, (char **)&value, 10);
		    *cpt++ = n & 0xff;
		    *cpt++ = (n >> 8) & 0xff;
		    break;
		    
		case 'i': case 'I':
		    n = strtoll((char *)value, (char **)&value, 10);
		    *cpt++ = n & 0xff;
		    *cpt++ = (n >> 8) & 0xff;
		    *cpt++ = (n >>16) & 0xff;
		    *cpt++ = (n >>24) & 0xff;
		    break;
		    
		case 'f': {
		    union {
			float f;
			int i;
		    } u;

		    u.f = strtod((char *)value, (char **)&value);
		    *cpt++ = (u.i) & 0xff;
		    *cpt++ = ((u.i) >> 8) & 0xff;
		    *cpt++ = ((u.i) >>16) & 0xff;
		    *cpt++ = ((u.i) >>24) & 0xff;
		    break;
		}
		}
		count++;
	    }
	    if (value != cpf) {
		fprintf(stderr, "Malformed %c%c:B:... auxiliary field\n",
			key[0], key[1]);
		value = cpf;
	    }
	    *sz++ = count & 0xff;
	    *sz++ = (count >> 8) & 0xff;
	    *sz++ = (count >>16) & 0xff;
	    *sz++ = (count >>24) & 0xff;
	    break;
	}

	default:
	    fprintf(stderr, "Unknown aux format code '%c'\n", key[3]);
	    break;
	}

	if (*cpf == '\t')
	    cpf++;
    }

 skip_aux:
    *cpt++ = 0;

    bs->blk_size = (unsigned char *)cpt - (unsigned char *)&bs->ref - 1;
    if (bs->blk_size >= bs->alloc)
	abort();

    return 1;
}

int bam_name2ref(bam_file_t *b, char *ref) {
    HashItem *hi;

    if (!b->ref_hash || !ref)
	return -1;

    hi = HashTableSearch(b->ref_hash, ref, strlen(ref));
    return hi ? hi->data.i : -1;
}

/*
 * Fills out the next bam_seq_t struct.
 * bs must be non-null, but *bs may be NULL or an existing bam_seq_t pointer.
 * This function will alloc and/or grow the memory accordingly, allowing for
 * efficient reuse.
 *
 * Returns 1 on success
 *         0 on eof
 *        -1 on error
 */
#ifdef ALLOW_UAC
int bam_next_seq(bam_file_t *b, bam_seq_t **bsp) {
    int32_t blk_size, blk_ret;
    bam_seq_t *bs;

    b->line++;

    if (!b->bam)
	return sam_next_seq(b, bsp);

    if (b->next_len > 0) {
	blk_size = b->next_len;
    } else {
	if (4 != bam_read(b, &blk_size, 4))
	    return 0;
	blk_size = le_int4(blk_size);
    }

    if (!*bsp || blk_size > (*bsp)->blk_size) {
	/* 20 extra is for bs->alloc to bs->cigar_len plus next_len */
	if (!(*bsp = realloc(*bsp, blk_size+20)))
	    return -1;
	(*bsp)->alloc = blk_size+20;
	(*bsp)->blk_size = blk_size;
    }
    bs = *bsp;
    
    if ((blk_ret = bam_read(b, &bs->ref, blk_size+4)) == 0)
	return 0;

    if (blk_size+4 != blk_ret) {
	if (blk_size != blk_ret) {
	    return -1;
	} else {
	    b->next_len = 0;
	    ((char *)(&bs->ref))[blk_size] = 0;
	}
    } else {
	memcpy(&b->next_len, &((char *)(&bs->ref))[blk_size], 4);
	((char *)(&bs->ref))[blk_size] = 0;
    }
    b->next_len = le_int4(b->next_len);

    bs->blk_size  = blk_size;
    bs->ref       = le_int4(bs->ref);
    bs->pos       = le_int4(bs->pos);
    bs->bin_mq_nl = le_int4(bs->bin_mq_nl);
    bs->flag_nc   = le_int4(bs->flag_nc);
    bs->len       = le_int4(bs->len);
    bs->mate_ref  = le_int4(bs->mate_ref);
    bs->mate_pos  = le_int4(bs->mate_pos);
    bs->ins_size  = le_int4(bs->ins_size);

    /* Unpack flag_nc into separate flag & cigar_len */
    bs->cigar_len = bs->flag_nc & 0xffff;
    
    if (10 == be_int4(10)) {
	int i, cigar_len = bam_cigar_len(bs);
	uint32_t *cigar = bam_cigar(bs);
	for (i = 0; i < cigar_len; i++) {
	    cigar[i] = le_int4(cigar[i]);
	}
    }

    return 1;
}

#else

int bam_next_seq(bam_file_t *b, bam_seq_t **bsp) {
    int32_t blk_size, blk_ret;
    bam_seq_t *bs;

    if (!b->bam)
	return sam_next_seq(b, bsp);

    if (b->next_len > 0) {
	blk_size = b->next_len;
    } else {
	if (4 != bam_read(b, &blk_size, 4))
	    return 0;
	blk_size = le_int4(blk_size);
    }

    if (!*bsp || blk_size > (*bsp)->blk_size) {
	if (!(*bsp = realloc(*bsp, blk_size+24)))
	    return -1;
	(*bsp)->alloc = blk_size+24;
	(*bsp)->blk_size = blk_size;
    }
    bs = *bsp;

    /* The fixed-sized fields */
    if ((blk_ret = bam_read(b, &bs->ref, 32)) == 0)
	return 0;

    if (blk_ret != 32)
	return -1;

    //bs->blk_size  = blk_size;
    bs->ref       = le_int4(bs->ref);
    bs->pos       = le_int4(bs->pos);
    bs->bin_mq_nl = le_int4(bs->bin_mq_nl);
    bs->flag_nc   = le_int4(bs->flag_nc);
    bs->len       = le_int4(bs->len);
    bs->mate_ref  = le_int4(bs->mate_ref);
    bs->mate_pos  = le_int4(bs->mate_pos);
    bs->ins_size  = le_int4(bs->ins_size);

    /* Unpack flag_nc into separate flag & cigar_len */
    bs->cigar_len = bs->flag_nc & 0xffff;

    /* Name */
    if (bam_read(b, &bs->data, bam_name_len(bs)) != bam_name_len(bs))
	return -1;

    /* Pad name out to end on a word-aligned boundary */
    blk_ret = blk_size - 32 - bam_name_len(bs);
    bam_set_name_len(bs, round4(bam_name_len(bs)));

    /* The remainder, word aligned */
    if (bam_read(b, &bs->data + bam_name_len(bs), blk_ret) != blk_ret)
	return -1;

    if (10 == be_int4(10)) {
	int i, cigar_len = bam_cigar_len(bs);
	uint32_t *cigar = bam_cigar(bs);
	for (i = 0; i < cigar_len; i++) {
	    bs[i] = le_int4(bs[i]);
	}
    }

    return 1;
}
#endif

/*
 * Looks for aux field 'key' and returns the value.
 * Returns NULL if not found.
 */
char *bam_aux_find(bam_seq_t *b, char *key) {
    char *cp = bam_aux(b);
    static int type_size[256] = {
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  4,  0,  4,  0,  0,  0,  0,  0,  7,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  4, 11,  0,  7,  0,  0,  7,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

    while (*cp) {
	int sz;

	//printf("%c%c:%c:?\n", cp[0], cp[1], cp[2]);

	if (cp[0] == key[0] && cp[1] == key[1])
	    return cp+3;
	
	if ((sz = type_size[cp[2]])) {
	    /* Fixed length fields */
	    cp += sz;
	} else {
	    /* Variable length, null terminated */
	    if (cp[2] == 'Z' || cp[2] == 'H') {
		cp += 3;
		while (*cp++)
		    ;
	    } else {
		fprintf(stderr, "Unknown aux type code %c(%d) in seq %s\n",
			cp[2], cp[2], bam_name(b));
		return NULL;
	    }
	}
    }

    return NULL;
}


/*
 * An iterator on bam aux fields. NB: This code is not reentrant or multi-
 * thread capable. The values returned are valid until the next call to
 * this function.
 * key:  points to an array of 2 characters (eg "RG", "NM")
 * type: points to an address of 1 character (eg 'Z', 'i')
 * val:  points to an address of a bam_aux_t union.
 *
 * Pass in *iter_handle as NULL to initialise the search and then
 * pass in the modified value on each subsequent call to continue the search.
 *
 * Returns 0 if the next value is valid, setting key, type and val.
 *        -1 when no more found.
 */
int bam_aux_iter(bam_seq_t *b, char **iter_handle,
		 char *key, char *type, bam_aux_t *val) {
    char *s;

    if (!iter_handle || !*iter_handle) {
	s = (char *)bam_aux(b);
    } else {
	s = *iter_handle;
    }

    /* We null terminate our aux list for ease */
    if (s[0] == 0)
	return -1;

    key[0] = s[0];
    key[1] = s[1];
    
    switch (s[2]) {
    case 'A':
	if (type) *type = 'A';
	if (val) val->i = *(s+3);
	s+=4;
	break;

    case 'C':
	if (type) *type = 'i';
	if (val) val->i = *(uint8_t *)(s+3);
	s+=4;
	break;

    case 'c':
	if (type) *type = 'i';
	if (val) val->i = *(int8_t *)(s+3);
	s+=4;
	break;

    case 'S':
	if (type) *type = 'i';
	if (val)
	    val->i = (uint16_t)((((unsigned char *)s)[3]<< 0) +
				(((unsigned char *)s)[4]<< 8));
	s+=5;
	break;

    case 's':
	if (type) *type = 'i';
	if (val)
	    val->i = (int16_t)((((unsigned char *)s)[3]<< 0) +
			       (((unsigned char *)s)[4]<< 8));
	s+=5;
	break;

    case 'I':
	if (type) *type = 'i';
	if (val)
	    val->i = (uint32_t)((((unsigned char *)s)[3]<< 0) +
				(((unsigned char *)s)[4]<< 8) +
				(((unsigned char *)s)[5]<<16) +
				(((unsigned char *)s)[6]<<24));
	s+=7;
	break;

    case 'i':
	if (type) *type = 'i';
	if (val)
	    val->i = (int32_t)((((unsigned char *)s)[3]<< 0) +
			       (((unsigned char *)s)[4]<< 8) +
			       (((unsigned char *)s)[5]<<16) +
			       (((unsigned char *)s)[6]<<24));
	s+=7;
	break;

    case 'f':
	if (type) *type = 'f';
	if (val) memcpy(&val->f, s+3, 4);
	s+=7;
	break;

    case 'd':
	if (type) *type = 'd';
	if (val) memcpy(&val->d, s+3, 8);
	s+=11;
	break;

    case 'Z': case 'H':
	if (type) *type = s[2];
	s+=3;
	if (val) val->s = s;
	while (*s++);
	break;

    case 'B': {
	uint32_t count;
	if (type) *type = 'B';
	count = (unsigned int)((((unsigned char *)s)[4]<< 0) +
			       (((unsigned char *)s)[5]<< 8) +
			       (((unsigned char *)s)[6]<<16) +
			       (((unsigned char *)s)[7]<<24));

	if (val) {
	    val->B.n = count;
	    val->B.t = s[3];
	    val->B.s = s+8;
	}
	s+=8;

	switch(val->B.t) {
	case 'c': case 'C': s +=   count; break;
	case 's': case 'S': s += 2*count; break;
	case 'i': case 'I': s += 4*count; break;
	case 'f':           s += 4*count; break;
	default:
	    fprintf(stderr, "Unknown sub-type '%c' for aux type 'B'\n",
		    val->B.t);
	    return -1;
	}
	break;
    }

    default:
	fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	return -1;
    }

    if (iter_handle)
	*iter_handle = s;

    return 0;
}

static int reg2bin(int start, int end) {
    if (end>start) end--;
    if ((start>>14) == (end>>14)) return ((1<<15)-1)/7 + (start>>14);
    if ((start>>17) == (end>>17)) return ((1<<12)-1)/7 + (start>>17);
    if ((start>>20) == (end>>20)) return ((1<<9 )-1)/7 + (start>>20);
    if ((start>>23) == (end>>23)) return ((1<<6 )-1)/7 + (start>>23);
    if ((start>>26) == (end>>26)) return ((1<<3 )-1)/7 + (start>>26);
    return 0;
}

/*
 * Constructs a bam_seq_t from components.
 * Ignores auxiliary tags for now.
 *
 * Returns -1 on error
 *          number of bytes written to bam_seq_t on success (ie tag offset)
 */
int bam_construct_seq(bam_seq_t *b, int s_size,
		      char *qname, size_t qname_len,
		      int flag,
		      int rname,      // Ref ID
		      int pos,
		      int start, int end, // aligned start/end coords
		      int mapq,
		      int ncigar, uint32_t *cigar,
		      int mrnm,       // Mate Ref ID
		      int mpos,
		      int isize,
		      int len,
		      char *seq,
		      char *qual) {
    char *cp;
    int i;
    static char L[256];

    if (L[0] == 0) {
	cp = "=ACMGRSVTWYHKDBN";
	memset(L, 15, 256);
	for (i = 0; i < 16; i++) {
	    L[cp[i]] = L[tolower(cp[i])] = i;
	}
    }

    b->ref = rname;
    b->pos = pos-1;
    bam_set_bin(b, reg2bin(start-1,end-1));
    bam_set_map_qual(b, mapq);
    bam_set_name_len(b, qname_len+1);
    bam_set_flag(b, flag);
    bam_set_cigar_len(b, ncigar);
    b->len = len;
    b->mate_ref = mrnm;
    b->mate_pos = mpos-1;
    b->ins_size = isize;

    if (s_size < sizeof(b) + 4*ncigar + (len+1)/2 + len)
	return -1;

    cp = bam_name(b);
    memcpy(cp, qname, qname_len);
    cp[qname_len] = 0;
    cp += qname_len+1;

    /* Cigar */
    for (i = 0; i < ncigar; i++) {
	uint32_t i32 = cigar[i];
	*cp++ = i32&0xff;
	i32>>=8;
	*cp++ = i32&0xff;
	i32>>=8;
	*cp++ = i32&0xff;
	i32>>=8;
	*cp++ = i32&0xff;
    }

    /* Seq */
    for (i = 0; i < len; i += 2) {
	uint8_t b2 = L[seq[i]]<<4;
	if (i+1 < len) b2 += L[seq[i+1]];
	*cp++ = b2;
    }

    /* Qual */
    if (qual) {
	memcpy(cp, qual, len);
	cp += len;
    } else {
	for (i = 0; i < len; i++) {
	    *cp++ = '*';
	}
    }

    *cp = 0; /* terminate aux list, for ease of parsing later */

    /* cp now points to the auxiliary tags if required */
    b->blk_size = (int)(cp-(char *)&b->ref);
    return (int)(cp-(char *)b);
}
		      

static char *append_int(char *cp, int32_t i) {
    int32_t j;

    if (i < 0) {
	*cp++ = '-';
	i = -i;

	if (i < 0) {
	    *cp++ = '2'; *cp++ = '1'; *cp++ = '4'; *cp++ = '7';
	    *cp++ = '4'; *cp++ = '8'; *cp++ = '3'; *cp++ = '6';
	    *cp++ = '4'; *cp++ = '8';
	    return cp;
	}
    } else if (i == 0) {
	*cp++ = '0';
	return cp;
    }

    //if (i < 10)         goto b0;
    if (i < 100)        goto b1;
    //if (i < 1000)       goto b2;
    if (i < 10000)      goto b3;
    //if (i < 100000)     goto b4;
    if (i < 1000000)    goto b5;
    //if (i < 10000000)   goto b6;
    if (i < 100000000)  goto b7;

 b9: if ((j = i / 1000000000)) {*cp++ = j + '0'; i -= j*1000000000; goto x8;}
 b8: if ((j = i / 100000000))  {*cp++ = j + '0'; i -= j*100000000;  goto x7;}
 b7: if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
 b6: if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5: if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
 b4: if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3: if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
 b2: if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1: if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
 b0: if (i)                     *cp++ = i + '0';
    return cp;

 x8: *cp++ = i / 100000000 + '0', i %= 100000000;
 x7: *cp++ = i / 10000000  + '0', i %= 10000000;
 x6: *cp++ = i / 1000000   + '0', i %= 1000000;
 x5: *cp++ = i / 100000    + '0', i %= 100000;
 x4: *cp++ = i / 10000     + '0', i %= 10000;
 x3: *cp++ = i / 1000      + '0', i %= 1000;
 x2: *cp++ = i / 100       + '0', i %= 100;
 x1: *cp++ = i / 10        + '0', i %= 10;
 x0: *cp++ = i             + '0';

    return cp;
}

/*
 * Unsigned version of above.
 * Only differs when the int has the top bit set (~2.15 billion and above).
 */
static char *append_uint(char *cp, uint32_t i) {
    uint32_t j;

    if (i == 0) {
	*cp++ = '0';
	return cp;
    }

    //if (i < 10)         goto b0;
    if (i < 100)        goto b1;
    //if (i < 1000)       goto b2;
    if (i < 10000)      goto b3;
    //if (i < 100000)     goto b4;
    if (i < 1000000)    goto b5;
    //if (i < 10000000)   goto b6;
    if (i < 100000000)  goto b7;

 b9: if ((j = i / 1000000000)) {*cp++ = j + '0'; i -= j*1000000000; goto x8;}
 b8: if ((j = i / 100000000))  {*cp++ = j + '0'; i -= j*100000000;  goto x7;}
 b7: if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
 b6: if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5: if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
 b4: if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3: if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
 b2: if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1: if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
 b0: if (i)                     *cp++ = i + '0';
    return cp;

 x8: *cp++ = i / 100000000 + '0', i %= 100000000;
 x7: *cp++ = i / 10000000  + '0', i %= 10000000;
 x6: *cp++ = i / 1000000   + '0', i %= 1000000;
 x5: *cp++ = i / 100000    + '0', i %= 100000;
 x4: *cp++ = i / 10000     + '0', i %= 10000;
 x3: *cp++ = i / 1000      + '0', i %= 1000;
 x2: *cp++ = i / 100       + '0', i %= 100;
 x1: *cp++ = i / 10        + '0', i %= 10;
 x0: *cp++ = i             + '0';

    return cp;
}

/*
 * This is set up so that count should never be more than BGZF_BUFF_SIZE. 
 * This has been chosen to deliberately be small enough such that the
 * bgzf header/footer + worst-case expansion (deflateBound() func) of 'buf'
 * are <= 65536, thus ensuring BGZF BSIZE is always 16-bit.
 */
static int bgzf_write(int fd, const void *buf, size_t count) {
    unsigned char blk[Z_BUFF_SIZE];
    z_stream s;
    int cdata_pos;
    int cdata_size;
    int cdata_alloc;
    int err;
    uint32_t crc;

    /* Initialise zlib stream */
    cdata_pos = 18;
    cdata_alloc = Z_BUFF_SIZE;
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)buf;
    s.avail_in = count;
    s.total_in = 0;
    s.next_out  = blk + cdata_pos;
    s.avail_out = cdata_alloc;
    s.total_out = 0;
    s.data_type = Z_BINARY;

    /* Compress it */
    err = deflateInit2(&s, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8,
		       Z_DEFAULT_STRATEGY);

    if (err != Z_OK) {
	fprintf(stderr, "zlib deflateInit2 error: %s\n", s.msg);
	return -1;
    }

    /* Encode to 'cdata' array */
    for (;s.avail_in;) {
	s.next_out = blk + cdata_pos;
	s.avail_out = cdata_alloc - cdata_pos;
	if (cdata_alloc - cdata_pos <= 0) {
	    fprintf(stderr, "Deflate produced larger output than expected. Abort\n"); 
	    return -1;
	}
	err = deflate(&s, Z_NO_FLUSH); // or Z_FINISH?
	cdata_pos = cdata_alloc - s.avail_out;
	if (err != Z_OK) {
	    fprintf(stderr, "zlib deflate error: %s\n", s.msg);
	    break;
	}
    }
    if (deflate(&s, Z_FINISH) != Z_STREAM_END) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    cdata_size = s.total_out;

    if (deflateEnd(&s) != Z_OK) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }

    assert(cdata_size <= 65536);

    /* Fill out gzip header */
    blk[ 0] =  31; // ID1
    blk[ 1] = 139; // ID2
    blk[ 2] =   8; // CM (deflate)
    blk[ 3] =   4; // FLAGS (FEXTRA)
    blk[ 4] = blk[5] = blk[6] = blk[7] = 0; // MTIME
    blk[ 8] =   0; // XFL
    blk[ 9] = 255; // OS (unknown)
    blk[10] =   6; // XLEN
    blk[11] =   0; // XLEN

    /* Extra BGZF fields */
    blk[12] =  66; // SI1
    blk[13] =  67; // SI2
    blk[14] =   2; // SLEN
    blk[15] =   0; // SLEN
    blk[16] = ((cdata_size + 25) >> 0) & 0xff;
    blk[17] = ((cdata_size + 25) >> 8) & 0xff;

    crc = crc32(0L, NULL, 0L);
    crc = crc32(crc, (unsigned char *)buf, count);
    blk[18+cdata_size+0] = (crc >> 0) & 0xff;
    blk[18+cdata_size+1] = (crc >> 8) & 0xff;
    blk[18+cdata_size+2] = (crc >>16) & 0xff;
    blk[18+cdata_size+3] = (crc >>24) & 0xff;
    blk[18+cdata_size+4] = (count >> 0) & 0xff;
    blk[18+cdata_size+5] = (count >> 8) & 0xff;
    blk[18+cdata_size+6] = (count >>16) & 0xff;
    blk[18+cdata_size+7] = (count >>24) & 0xff;

    //printf("count=%d/%x, cdata_size=%d/%x\n", count, count, cdata_size, cdata_size);
    if (18+cdata_size+8 != write(fd, blk, 18+cdata_size+8))
	return -1;

    return 0;
}

/*
 * Writes a single bam sequence object.
 * Returns 0 on success
 *        -1 on failure
 */
int bam_put_seq(bam_file_t *fp, bam_seq_t *b) {
    static unsigned char *buf = NULL;
    static size_t buf_len = 0;
    size_t len;
    char *auxh, aux_key[2], type;
    bam_aux_t val;

    if (!fp->binary) {
	/* SAM */
	unsigned char *end = fp->out + BGZF_BUFF_SIZE, *dat;
	int sz, i, n;

#define BF_FLUSH() \
	do {			  \
	    if (fp->out_p - fp->out != \
		write(fp->fd, fp->out, fp->out_p - fp->out)) \
		return -1;				     \
	    fp->out_p=fp->out;				\
	} while(0)

	/* QNAME */
	if (end - fp->out_p < (sz = bam_name_len(b))) BF_FLUSH();
	memcpy(fp->out_p, bam_name(b), sz-1); fp->out_p += sz-1;
	*fp->out_p++ = '\t';

	/* FLAG */
	if (end-fp->out_p < 5) BF_FLUSH();
	fp->out_p = append_int(fp->out_p, bam_flag(b)); *fp->out_p++ = '\t';

	/* RNAME */
	if (b->ref != -1) {
	    size_t l = strlen(fp->ref[b->ref].name);
	    if (end-fp->out_p < l+1) BF_FLUSH();
	    memcpy(fp->out_p, fp->ref[b->ref].name, l);
	    fp->out_p += l;
	} else {
	    if (end-fp->out_p < 2) BF_FLUSH();
	    *fp->out_p++ = '*';
	}
	*fp->out_p++ = '\t';

	/* POS */
	if (end-fp->out_p < 12) BF_FLUSH();
	fp->out_p = append_int(fp->out_p, b->pos+1); *fp->out_p++ = '\t';

	/* MAPQ */
	if (end-fp->out_p < 5) BF_FLUSH();
	fp->out_p = append_int(fp->out_p, bam_map_qual(b)); *fp->out_p++ = '\t';

	/* CIGAR */
	n = bam_cigar_len(b); dat = (char *)bam_cigar(b);
	for (i = 0; i < n; i++, dat+=4) {
	    uint32_t c = dat[0] + (dat[1]<<8) + (dat[2]<<16) + (dat[3]<<24);
	    if (end-fp->out_p < 13) BF_FLUSH();
	    fp->out_p = append_int(fp->out_p, c>>4);
	    *fp->out_p++="MIDNSHP=X"[c&15];
	}
	if (n==0) {
	    if (end-fp->out_p < 2) BF_FLUSH();
	    *fp->out_p++='*';
	}
	*fp->out_p++='\t';

	/* NRNM */
	if (b->mate_ref != -1) {
	    if (b->mate_ref == b->ref) {
		if (end-fp->out_p < 2) BF_FLUSH();
		*fp->out_p++ = '=';
	    } else {
		size_t l = strlen(fp->ref[b->mate_ref].name);
		if (end-fp->out_p < l+1) BF_FLUSH();
		memcpy(fp->out_p, fp->ref[b->mate_ref].name, l);
		fp->out_p += l;
	    }
	} else {
	    if (end-fp->out_p < 2) BF_FLUSH();
	    *fp->out_p++ = '*';
	}
	*fp->out_p++ = '\t';

	/* MPOS */
	if (end-fp->out_p < 12) BF_FLUSH();
	fp->out_p = append_int(fp->out_p, b->mate_pos+1); *fp->out_p++ = '\t';

	/* ISIZE */
	if (end-fp->out_p < 12) BF_FLUSH();
	fp->out_p = append_int(fp->out_p, b->ins_size); *fp->out_p++ = '\t';

	/* SEQ */
	n = (b->len+1)/2;
	dat = bam_seq(b);

	/* BAM encoding */
//	while (n) {
//	    int l = end-fp->out_p < n ? end-fp->out_p : n;
//	    memcpy(fp->out_p, dat, l); fp->out_p += l;
//	    n -= l; dat += l;
//	    if (end == fp->out_p) BF_FLUSH();
//	}
	if (b->len != 0) {
	    for (i = 0; i < b->len; i++) {
		if (end - fp->out_p < 3) BF_FLUSH();
		*fp->out_p++ = "=ACMGRSVTWYHKDBN"[*dat >> 4];
		if (++i >= b->len)
		    break;
		*fp->out_p++ = "=ACMGRSVTWYHKDBN"[*dat++ & 15];
	    }
	} else {
	    if (end - fp->out_p < 2) BF_FLUSH();
	    *fp->out_p++ = '*';
	}
	*fp->out_p++ = '\t';

	/* QUAL */
	n = b->len;
	dat = bam_qual(b);
	/* BAM encoding */
//	while (n) {
//	    int l = end-fp->out_p < n ? end-fp->out_p : n;
//	    memcpy(fp->out_p, dat, l); fp->out_p += l;
//	    n -= l; dat += l;
//	    if (end == fp->out_p) BF_FLUSH();
//	}
	if (b->len != 0) {
	    if (*dat == 0xff) {
		if (end - fp->out_p < 2) BF_FLUSH();
		*fp->out_p++ = '*';
		dat += b->len;
	    } else {
		for (i = 0; i < b->len; i++) {
		    if (end - fp->out_p < 3) BF_FLUSH();
		    *fp->out_p++ = *dat++ + '!';
		}
	    }
	} else {
	    if (end - fp->out_p < 2) BF_FLUSH();
	    *fp->out_p++ = '*';
	}

	/* Auxiliary tags */
	auxh = NULL;
	while (0 == bam_aux_iter(b, &auxh, aux_key, &type, &val)) {
	    if (end - fp->out_p < 20) BF_FLUSH();
	    *fp->out_p++ = '\t';
	    *fp->out_p++ = aux_key[0];
	    *fp->out_p++ = aux_key[1];
	    *fp->out_p++ = ':';
	    *fp->out_p++ = type;
	    *fp->out_p++ = ':';
	    switch(type) {
	    case 'A':
		*fp->out_p++ = val.i;
		break;

	    case 'C':
		fp->out_p = append_int(fp->out_p, (uint8_t)val.i);
		break;

	    case 'c':
		fp->out_p = append_int(fp->out_p, (int8_t)val.i);
		break;

	    case 'S':
		fp->out_p = append_int(fp->out_p, (uint16_t)val.i);
		break;

	    case 's':
		fp->out_p = append_int(fp->out_p, (int16_t)val.i);
		break;

	    case 'I':
		fp->out_p = append_int(fp->out_p, (uint32_t)val.i);
		break;

	    case 'i':
		fp->out_p = append_int(fp->out_p, (int32_t)val.i);
		break;

	    case 'f':
		fp->out_p += sprintf(fp->out_p, "%g", val.f);
		break;

	    case 'd':
		fp->out_p += sprintf(fp->out_p, "%g", val.d);
		break;

	    case 'Z':
	    case 'H': {
		size_t l = strlen(val.s);
		if (end - fp->out_p < l+2) BF_FLUSH();
		memcpy(fp->out_p, val.s, l);
		fp->out_p += l;
		break;
	    }

	    case 'B': {
		uint32_t count = val.B.n, sz, j, i_end;
		unsigned char *s = val.B.s;
		*fp->out_p++ = val.B.t;

		/*
		 * Chew through count items 4000 at a time.
		 * This is because 4000*14 (biggest %g output plus comma?)
		 * is just shy of 64k, so we avoid buffer overflows.
		 */
		switch (val.B.t) {
		case 'C': case 'c': sz = 4; break;
		case 'S': case 's': sz = 6; break;
		default:            sz = 14; break;
		}

		for (j = 0; j < count; j += 4000) {
		    int i_start = j;
		    int i_end = j + 4000 < count ? j + 4000 : count;

		    if (end - fp->out_p < 5+(i_end-i_start)*sz) BF_FLUSH();

		    switch (val.B.t) {
			int i;
		    case 'C':
			for (i = i_start; i < i_end; i++, s++) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_int(fp->out_p, (uint8_t)s[0]);
			}
			break;

		    case 'c':
			for (i = i_start; i < i_end; i++, s++) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_int(fp->out_p, (int8_t)s[0]);
			}
			break;

		    case 'S':
			for (i = i_start; i < i_end; i++, s+=2) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_int(fp->out_p,
						   (uint16_t)((s[0] << 0) +
							      (s[1] << 8)));
			}
			break;

		    case 's':
			for (i = i_start; i < i_end; i++, s+=2) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_int(fp->out_p,
						   (int16_t)((s[0] << 0) +
							     (s[1] << 8)));
			}
			break;

		    case 'I':
			for (i = i_start; i < i_end; i++, s+=4) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_uint(fp->out_p,
						    (uint32_t)((s[0] << 0) +
							       (s[1] << 8) +
							       (s[2] <<16) +
							       (s[3] <<24)));
			}
			break;

		    case 'i':
			for (i = i_start; i < i_end; i++, s+=4) {
			    *fp->out_p++ = ',';
			    fp->out_p = append_int(fp->out_p,
						   (int32_t)((s[0] << 0) +
							     (s[1] << 8) +
							     (s[2] <<16) +
							     (s[3] <<24)));
			}
			break;

		    case 'f': {
			union {
			    float f;
			    unsigned char c[4];
			} u;
			for (i = i_start; i < i_end; i++, s+=4) {
			    *fp->out_p++ = ',';
			    u.c[0] = s[0];
			    u.c[1] = s[1];
			    u.c[2] = s[2];
			    u.c[3] = s[3];
			    fp->out_p += sprintf(fp->out_p, "%g", u.f);
			}
			break;
		    }

		    default:
			fprintf(stderr, "Unhandled sub-type of aux type B\n");
		    }
		}
		break;
	    }

	    default:
		fprintf(stderr, "Unhandled auxilary type '%c' in "
			"bam_put_seq()\n", type);
	    }
	}

	*fp->out_p++ = '\n';

    } else {
	/* BAM */
	unsigned char *end = fp->out + BGZF_BUFF_SIZE, *ptr;
	size_t to_write;

#define CF_FLUSH() \
	do {			  \
	    if (bgzf_write(fp->fd, fp->out, fp->out_p - fp->out))	\
		return -1;				     \
	    fp->out_p=fp->out;				\
	} while(0)


	if (end - fp->out_p < b->blk_size+4) CF_FLUSH();
	*fp->out_p++ = (b->blk_size >> 0) & 0xff;
	*fp->out_p++ = (b->blk_size >> 8) & 0xff;
	*fp->out_p++ = (b->blk_size >>16) & 0xff;
	*fp->out_p++ = (b->blk_size >>24) & 0xff;

#if 0
	memcpy(fp->out_p, (unsigned char *)&b->ref, b->blk_size);
	fp->out_p += b->blk_size;
#else

	to_write = b->blk_size;
	ptr = (unsigned char *)&b->ref;
	//printf("blk_size=%d\n", b->blk_size);
	do {
	    size_t blk_len = MIN(to_write, end - fp->out_p);
	    memcpy(fp->out_p, ptr, blk_len);
	    fp->out_p += blk_len;
	    to_write -= blk_len;
	    ptr += blk_len;

	    if (to_write) {
		//printf("flushing %d+%d\n",
		//       (int)(ptr-(unsigned char *)&b->ref),
		//       (int)(fp->out_p-fp->out));
		CF_FLUSH();
	    }
	} while(to_write > 0);
#endif

	//len = 4*9 + bam_name_len(b) + 1 + 4*bam_cigar_len(b) +
	//    (bam_seq_len(b)+1)/2 + bam_seq_len(b);
    }

    return 0;
}

/*
 * Writes a header block.
 * Returns 0 for success
 *        -1 for failure
 */
int bam_write_header(bam_file_t *out) {
    char *header, *hp;
    size_t hdr_size;
    int i;

    hdr_size = 12 + out->header_len+1;
    for (i = 0; i < out->nref; i++) {
	hdr_size += strlen(out->ref[i].name)+1 + 8;
    }
    if (NULL == (hp = header = malloc(hdr_size)))
	return -1;

    if (out->binary) {
	*hp++ = 'B'; *hp++ = 'A'; *hp++ = 'M'; *hp++ = 1;
	*hp++ = (out->header_len >> 0) & 0xff;
	*hp++ = (out->header_len >> 8) & 0xff;
	*hp++ = (out->header_len >>16) & 0xff;
	*hp++ = (out->header_len >>24) & 0xff;
    }
    memcpy(hp, out->header, out->header_len);
    hp += out->header_len;

    if (out->binary) {
	int i;

	*hp++ = (out->nref >> 0) & 0xff;
	*hp++ = (out->nref >> 8) & 0xff;
	*hp++ = (out->nref >>16) & 0xff;
	*hp++ = (out->nref >>24) & 0xff;

	for (i = 0; i < out->nref; i++) {
	    size_t l = strlen(out->ref[i].name)+1;
	    *hp++ = (l>> 0) & 0xff;
	    *hp++ = (l>> 8) & 0xff;
	    *hp++ = (l>>16) & 0xff;
	    *hp++ = (l>>24) & 0xff;

	    strcpy(hp, out->ref[i].name);
	    hp += l;

	    l = out->ref[i].len;
	    *hp++ = (l>> 0) & 0xff;
	    *hp++ = (l>> 8) & 0xff;
	    *hp++ = (l>>16) & 0xff;
	    *hp++ = (l>>24) & 0xff;
	}
    }

    if (out->binary) {
	int len = hp-header;
	char *cp = header;

	while (len) {
	    int sz = BGZF_BUFF_SIZE < len ? BGZF_BUFF_SIZE : len;
	    if (bgzf_write(out->fd, cp, sz))
		return -1;
	    cp  += sz;
	    len -= sz;
	}
    } else {
	if (hp-header != write(out->fd, header, hp-header))
	    return -1;
    }

    free(header);

    return 0;
}


#ifdef TEST_SAM_TO_BAM
/*
 * (cd gap5;gcc -O2  -I/software/badger/opt/xz_utils/include -DUSE_NON_CONST  -I/nfs/users/nfs_j/jkb/staden/trunk/src/build.seq1/../gap5 -I/nfs/users/nfs_j/jkb/staden/trunk/src/build.seq1/../Misc -I"/usr/include/tcl8.4/tk-private/generic" -I"/usr/include/tcl8.4/tk-private/unix" -I"/usr/include/tcl8.4/tcl-private/generic" -I"/usr/include/tcl8.4/tcl-private/unix" -I/software/badger/include -I/nfs/users/nfs_j/jkb/staden/trunk/src/build.seq1/../tk_utils -I/nfs/users/nfs_j/jkb/staden/trunk/src/build.seq1/../seq_utils -I/nfs/users/nfs_j/jkb/staden/trunk/src/build.seq1/../primer3/src   -I/software/badger/opt/xz_utils/include -I/nfs/users/nfs_j/jkb/staden/trunk/src/build.seq1  -DSVN_VERSION="3063" -fPIC  /nfs/users/nfs_j/jkb/staden/trunk/src/build.seq1/../gap5/bam.c -o ../bam_main -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DTEST_SAM_TO_BAM -L /nfs/sam_scratch/jkb/staden.i686/lib -lstaden-read) 
 */
int main(int argc, char **argv) {
    size_t i;
    bam_file_t *in, *out;
    bam_seq_t *s;
    int nr = 0;
    int nm[2];
    char *mode = "rb";


    if (argc != 3) {
	fprintf(stderr, "Usage: %s sam/bam-in-file bam-out-file\n", argv[0]);
	return 1;
    }
    
    {
	size_t l = strlen(argv[1]);
	if (l >= 4 && argv[1][l-3] == 's')
	    mode = "r";
    }

    in = bam_open(argv[1], mode);
    if (!in) {
	fprintf(stderr, "Failed to open bam file %s\n", argv[1]);
	return 1;
    }

    {
	size_t l = strlen(argv[2]);
	if (l >= 4 && argv[2][l-3] == 's')
	    mode = "w";
	else
	    mode = "wb";
    }

    out = bam_open(argv[2], mode);
    if (!out) {
	fprintf(stderr, "Failed to open bam file %s\n", argv[2]);
	return 1;
    }

    /* Copy header and refs from in to out, for writing purposes */
    out->header     = in->header;
    out->header_len = in->header_len;
    out->ref        = in->ref;
    out->nref       = in->nref;

    if (in->header) {
	if (bam_write_header(out))
	    return 1;
    }

    nm[0] = nm[1] = 0;
    s = NULL;
    while (bam_next_seq(in, &s) > 0) {
	out->ref = in->ref; // tmp hack
	if (-1 == bam_put_seq(out, s))
	    return 1;
	nr++;
	nm[(bam_flag(s) & BAM_FUNMAP) > 0]++;
    }

    printf("Mapped      = %d\n",nm[0]);
    printf("Unmpped     = %d\n",nm[1]);
    printf("Total reads = %d\n",nr);

    out->header = NULL;
    out->ref    = NULL;

    bam_close(in);
    bam_close(out);

    return 0;
}
#endif

#ifdef TEST_BAM_MAIN

//#define RG_COUNT

#ifdef RG_COUNT
#include "io_lib/hash_table.h"
#endif

int main(int argc, char **argv) {
    size_t i;
    bam_file_t *b;
    bam_seq_t *s;
    int nr = 0;
    int nm[2];
    char *mode = "rb";

#ifdef RG_COUNT
    HashTable *h = HashTableCreate(4, HASH_FUNC_HSIEH | HASH_DYNAMIC_SIZE);
#endif

    if (argc != 2) {
	fprintf(stderr, "Usage: %s bam_file\n", argv[0]);
	return 1;
    }
    
    {
	size_t l = strlen(argv[1]);
	if (l >= 4 && argv[1][l-3] == 's')
	    mode = "r";
    }

    b = bam_open(argv[1], mode);
    if (!b) {
	fprintf(stderr, "Failed to open bam file\n");
	return 1;
    }

    //b->no_aux = 1;

    nm[0] = nm[1] = 0;
    s = NULL;
    while (bam_next_seq(b, &s) > 0) {
	nr++;
	nm[(bam_flag(s) & BAM_FUNMAP) > 0]++;
#ifdef RG_COUNT
	//printf("%s RG='%s'\n", bam_name(s), bam_aux_find(s, "RG"));
	{
	    char *rg = bam_aux_find(s, "RG");
	    HashItem *hi;
	    if (!rg)
		rg = "unknown";

	    hi = HashTableSearch(h, rg, 0);
	    if (!hi) {
		HashData hd;
		hd.p = calloc(2,sizeof(int));
		hi = HashTableAdd(h, rg, 0, hd, NULL);
	    }
	    ((int *)hi->data.p)[(bam_flag(s) & BAM_FUNMAP) > 0]++;
	}
#endif
    }

    printf("Mapped      = %d\n",nm[0]);
    printf("Unmpped     = %d\n",nm[1]);
    printf("Total reads = %d\n",nr);

#ifdef RG_COUNT
    {
	HashIter *iter = HashTableIterCreate();
	HashItem *hi;
	
	printf("Mapped\tUnmap\tRead-group\n");
	while (hi = HashTableIterNext(h, iter)) {
	    printf("%d\t%d\t%.*s\n",
		   ((int *)hi->data.p)[0],
		   ((int *)hi->data.p)[1],
		   hi->key_len, hi->key);
	}
	HashTableIterDestroy(iter);
	HashTableDestroy(h, 1);
    }
#endif

    bam_close(b);

    return 0;
}

#endif /* TEST_MAIN */
