/*
 * This file contains the heap allocation system for Gap5 - a disk
 * based block allocator.  It replaces the AVL trees used in Gap4.
 * NB: This is basically a malloc() implementation.
 *
 * The main database file is a heap of allocated objects using
 * boundary data. The boundary data is the key change from Gap4 and
 * the simplification that allows us to move away from AVL trees to a
 * more malloc-style system.
 *
 * Here we implement a segregated fit allocator with many doubly
 * linked rings, using a next-fit strategy (via a roving pointer)
 * and block coalescing on free.
 *
 * All allocation requests are rounded to the next 8 byte boundary.
 * We then return from the smallest available pool, splitting free
 * blocks if required.
 * Pools are of sizes 16 to 1024 in steps of 8 inclusive, and then
 * every 16, 32, 64, 128, 256, 512 ... bytes in size from then on,
 * giving a total of 155 pools for the full 32-bit size range.
 *
 * Every allocated block consists of:
 *
 *   31bits   1bit   8*N bits             ? bits    1 bit
 * +--------+------+----------------...-+---...---+------+
 * | Length | Free | Allocated data ... | Padding | Free |
 * +--------+------+----------------...-+---...---+------+
 *
 * Here Padding is extra storage to ensure the entire total length is
 * a multiple of 64-bits.
 *
 * A free block consists of:
 *
 *   31bits   1bit  64bits         64bits  31bits   1bit
 * +--------+------+------+       +------+--------+------+
 * | Length | Free | Prev |  ...  | Next | Length | Free |
 * +--------+------+------+       +------+--------+------+
 *
 * The Prev and Next values are pointers - offsets into the previous
 * or next free block. The ends of the linked list are linked together
 * to form a ring.
 *
 * Here we note that the minimum value of N supported is 16 bytes
 * (64-bits) in order to allow for Prev and Next. This yields a
 * Padding value of 31 bits in order to round the entire allocated
 * storage including boundary tags to the next multiple of 64-bits (24
 * bytes).
 *
 * The minimum overhead for an allocated block is 5 bytes and averages
 * at 9 bytes.
 */

/*
 * To-do:
 * Each pool should remember the largest item found within it to shortcut
 * the linear searching through it.
 */

#include <staden_config.h>

#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <string.h>

#include <sys/time.h>
#include <time.h>

#include "g-alloc.h"
#include "os.h"

#define MIN_LEN 24

//#define DEBUG
#define DEBUG_FP stdout
//#define VALGRIND

/* An allocated block on the heap */
typedef struct {
    uint64_t pos;
    uint32_t len;
    uint32_t blen; /* length of block immediately before us */
    uint64_t prev; /* previous/next pointers to form a ring */
    uint64_t next;
    char free;     /* non-zero => free */
    char bfree;    /* free status for block immediately before us */
} block_t;

/* Rounds a length to the next multiple of 8 */
#define SIZE_ROUND(l) (((l) & 7) ? ((l) & ~7) + 8 : (l))

/* Identifies which pool to use given a specific length */
static inline int pool(uint32_t l) {
    int pool;

    if (l < 16)
	return 0;

    if (l <= 1024)
	return (l>>3) - 2;
    
    /* Slow method. For faster way, see:
     * http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup
     */
    l-=1024-8;
    l >>= 3;
    pool = 126;
    while (l >>= 1)
	pool++;

    return pool;
}

static int get_block(dheap_t *h, int64_t offset, block_t *b) {
    union {
	unsigned char c[24];
	int32_t i32[6];
	int64_t i64[3];
    } data;

    /* Attempt to read the footer from the previous block too */
    if (offset >= 4) {
	if (-1 == lseek(h->fd, offset-4, SEEK_SET))
	    return -1;
	if (24 != read(h->fd, data.c, 24))
	    return -1;
    } else {
	if (-1 == lseek(h->fd, offset, SEEK_SET))
	    return -1;
	if (20 != read(h->fd, &data.c[4], 20))
	    return -1;
	data.c[0] = data.c[1] = data.c[2] = data.c[3] = 0;
    }

    b->pos = offset;
    
    //b->blen = be_int4(*(int32_t *)&data[0]);
    b->blen = be_int4(data.i32[0]);
    b->bfree = b->blen & 1;
    b->blen &= ~1;

    //b->len = be_int4(*(int32_t *)&data[4]);
    b->len = be_int4(data.i32[1]);
    b->free = b->len & 1;
    b->len &= ~1;

    if (b->free) {
	//b->prev = be_int8(*(int64_t *)&data[8]);
	//b->next = be_int8(*(int64_t *)&data[16]);
	b->prev = be_int8(data.i64[1]);
	b->next = be_int8(data.i64[2]);
    } else {
	b->prev = b->next = 0;
    }

    return 0;
}

/*
 * Destructively writes a heap block, overwriting the contents of the block
 * itself unless header_only is specified. With header_only true we also
 * do not write the 'length|free' footer data either.
 *
 * Zero indicates whether to zero the block. We don't bother doing this if
 * we're shrinking an already free block. (Zeroing the data is relevant for
 * allocation only, and even then it's questionable.)
 * 
 * Returns 0 for sucess
 *        -1 for failure
 */
#define PB_SIZE 1024
static int put_block(dheap_t *h, block_t *b, int header_only, int zero) {
    unsigned char header[PB_SIZE];
    uint32_t len;
    int32_t i32;
    int64_t i64;

    if (-1 == lseek(h->fd, b->pos, SEEK_SET))
	return -1;

    len = b->len | b->free;

    //*(uint32_t *)&header[0] = be_int4(len);
    i32 = be_int4(len);
    header[0] = ((unsigned char *)&i32)[0];
    header[1] = ((unsigned char *)&i32)[1];
    header[2] = ((unsigned char *)&i32)[2];
    header[3] = ((unsigned char *)&i32)[3];

    //*(uint64_t *)&header[4] = be_int8(b->prev);
    i64 = be_int8(b->prev);
    header[ 4] = ((unsigned char *)&i64)[0];
    header[ 5] = ((unsigned char *)&i64)[1];
    header[ 6] = ((unsigned char *)&i64)[2];
    header[ 7] = ((unsigned char *)&i64)[3];
    header[ 8] = ((unsigned char *)&i64)[4];
    header[ 9] = ((unsigned char *)&i64)[5];
    header[10] = ((unsigned char *)&i64)[6];
    header[11] = ((unsigned char *)&i64)[7];
    
    //*(uint64_t *)&header[12] = be_int8(b->next);
    i64 = be_int8(b->next);
    header[12] = ((unsigned char *)&i64)[0];
    header[13] = ((unsigned char *)&i64)[1];
    header[14] = ((unsigned char *)&i64)[2];
    header[15] = ((unsigned char *)&i64)[3];
    header[16] = ((unsigned char *)&i64)[4];
    header[17] = ((unsigned char *)&i64)[5];
    header[18] = ((unsigned char *)&i64)[6];
    header[19] = ((unsigned char *)&i64)[7];

    if (header_only) {
	/* Header only is used when updating the next/prev pointers */
	if (20 != write(h->fd, header, 20))
	    return -1;

    } else {
#ifdef VALGRIND
	/* For valgrind we always zero data */
	zero = 1;
#endif
	
	/*
	 * We zero the data block. This has two effects.
	 *
	 * 1. It appears that writing the header + gap + footer in a single
	 *    write is faster, especially over NFS to a BlueArc. Sometimes
	 *    20 fold faster!
	 *
	 * 2. It means in g-request.c we no longer need to call the
	 *    write_zeros function which was doing this anyway.
	 */
	if (b->len <= PB_SIZE) {
	    assert(b->len >= MIN_LEN);
	    if (zero)
		memset(&header[20], 0, b->len-4-20);
	    *(uint32_t *)&header[b->len-4] = be_int4(len);
	    if (b->len != write(h->fd, header, b->len))
		return -1;
	} else {
	    if (zero) {
		/* A long block, so calloc it instead */
		unsigned char *h2 = calloc(1, b->len);
		if (!h2)
		    return -1;
		memcpy(h2, header, 20);
		*(uint32_t *)&h2[b->len-4] = be_int4(len);
		
		if (b->len != write(h->fd, h2, b->len))
		    return -1;
		
		free(h2);
	    } else {
		/* A long block, so we write the header and footer only */
		uint32_t be_len = be_int4(len);
		if (20 != write(h->fd, header, 20))
		    return -1;
		if (-1 == lseek(h->fd, b->len-4-20, SEEK_CUR))
		    return -1;
		if (4 != write(h->fd, &be_len, 4))
		    return -1;
	    }
	}
    }

    return 0;
}

static int write_pool(dheap_t *h, int pool) {
    int64_t pbe = be_int8(h->pool[pool]);

    if (-1 == lseek(h->fd, 8*pool, SEEK_SET))
	return -1;

    if (8 != write(h->fd, &pbe, 8))
	return -1;

    return 0;
}

static int unlink_block(dheap_t *h, block_t *b) {
    int p = pool(b->len);

    if (b->next != b->pos) {
	block_t nb;

	if (-1 == get_block(h, b->next, &nb))
	    return -1;
	nb.prev = b->prev;
	put_block(h, &nb, 1, 0);
    }

    if (b->prev != b->pos) {
	block_t pb;

	if (-1 == get_block(h, b->prev, &pb))
	    return -1;
	pb.next = b->next;
	put_block(h, &pb, 1, 0);
    }

    if (h->pool[p] == b->pos) {
	/* Pool points to this item */
	if (b->next == b->pos) {
	    h->pool[p] = 0;
	    h->maxsz[p] = 0;
	} else {
	    /* maxsz is still a valid upper-bound */
	    h->pool[p] = b->next;
	}

	write_pool(h, p);
    }

    return 0;
}

static int link_block(dheap_t *h, block_t *b) {
    int p = pool(b->len);

    b->free = 1;

    if (h->maxsz[p] < b->len)
	h->maxsz[p] = b->len;

    if (h->pool[p]) {
	/* Link into existing ring */
	block_t pb, nb;

	if (-1 == get_block(h, h->pool[p], &pb))
	    return -1;
	if (-1 == get_block(h, pb.next, &nb))
	    return -1;
	b->prev = h->pool[p];
	b->next = pb.next;
	pb.next = b->pos;
	if (pb.pos != nb.pos)
	    nb.prev = b->pos;
	else
	    pb.prev = b->pos;

	put_block(h, &pb, 1, 0);
	put_block(h, b,   0, 0);

	if (pb.pos != nb.pos)
	    put_block(h, &nb, 0, 0);
    } else {
	/* No previous, so create a new ring */
	h->timer++;
	h->pool[p] = b->pos;
	b->prev = b->pos;
	b->next = b->pos;
	put_block(h, b, 0, 0);
    }

    write_pool(h, p);

    return 0;
}

/*
 * Load the heap meta-data from an existing heap-file.
 * Returns dheap_t pointer on success.
 *         NULL on failure
 */
dheap_t *heap_fdload(int fd) {
    dheap_t *h;
    int i;
    struct stat sb;

    if (NULL == (h = malloc(sizeof(*h))))
	return NULL;
    h->fd = fd;

    if (8*NPOOLS != read(h->fd, &h->pool[0], 8*NPOOLS))
	return NULL;

    for (i = 0; i < NPOOLS; i++)
	h->pool[i] = be_int8(h->pool[i]);

    if (-1 == fstat(h->fd, &sb))
	return NULL;

    h->wilderness = sb.st_size;

    for (i = 0; i < NPOOLS; i++) {
	h->next_free_pool[i] = 0;
	h->next_free_time[i] = 0;
	h->maxsz[i] = 0;
    }
    h->timer = 1;

    /* heap_check(h); */

    return h;
}

dheap_t *heap_load(char *file, int mode) {
    int fd;

    if (-1 == (fd = open(file, mode))) {
	return NULL;
    }

    return heap_fdload(fd);
}


/*
 * Create a new heap in 'file'.
 * Returns dheap_t pointer on success.
 *         NULL on failure
 */
dheap_t *heap_create(char *file) {
    int fd, i;
    uint64_t pool[NPOOLS];

    if (-1 == (fd = open(file, O_RDWR|O_CREAT|O_TRUNC, 0666)))
	return NULL;

    for (i = 0; i < NPOOLS; i++) {
	pool[i] = 0;
    }

    write(fd, pool, 8*NPOOLS);
    close(fd);

    return heap_load(file, O_RDWR);
}

/*
 * Deallocates memory used by the heap, optionally closing the internal
 * filedescriptor too if close_me is true.
 */
void heap_destroy(dheap_t *h, int close_me) {
    if (!h)
	return;

    /* heap_check(h); */
    if (close_me)
	close(h->fd);

    free(h);
}

static int64_t wilderness_allocate(dheap_t *h, uint32_t length) {
    block_t b;

    b.free = 0;
    b.len = length;
    b.pos = h->wilderness;
    b.prev = 0;
    b.next = 0;

    h->wilderness += length;
    put_block(h, &b, 0, 1);
    
#ifdef DEBUG
    static int count = 0;
    fprintf(DEBUG_FP, "%ld (W) %d\n", b.pos+4, count++);
#endif
    return b.pos + 4;
}

/*
 * Allocates an object of a specific length from the disk heap 'h'.
 * Returns the file offset on success
 *        -1 on failure
 */
int64_t heap_allocate(dheap_t *h, uint32_t length, uint32_t *allocated) {
    int p, pred, porig;
    uint64_t rover, orig;

    /* Round size to the next multiple of 8 after boundary tags */
    length = SIZE_ROUND(length+5);
    if (length < MIN_LEN)
	length = MIN_LEN;

    if (allocated)
	*allocated = length-5;

#ifdef DEBUG
    fprintf(DEBUG_FP, "H-%p: heap_allocate %d => ", h, length);
    fflush(DEBUG_FP);
#endif

    /* Search for the first item in the pool that's big enough */
    porig = p = pool(length);
    if (!h->pool[p] &&
	h->next_free_pool[p] &&
	h->next_free_time[p] == h->timer)
	p = h->next_free_pool[p];
    pred = p;

    assert(pred >= porig);

    //    if (length > 70000)
    //	printf("  largest=%d  ", heap_largest_check(h));

    for (; p < NPOOLS && p-pred < 75; p++) {
	uint64_t ms = 0;
	orig = rover = h->pool[p];

	if (!rover)
	    continue;

	if (h->maxsz[p] && h->maxsz[p] < length)
	    continue;

	/* Found a pool with free blocks, so search for one large enough */
	do {
	    block_t b;

	    assert(rover != 0);
	    if (-1 == get_block(h, rover, &b))
		return -1;

	    if (b.len >= length) {
		assert(p >= pred);

		unlink_block(h, &b);
		b.free = 0;

		/*
		 * Min. free block size = MIN_LEN.
		 * If the remainder is less than that we'll just use
		 * the entire thing.
		 */
		if (b.len - length < MIN_LEN) {
		    put_block(h, &b, 0, 1);
		    if (!h->pool[porig])
			h->next_free_pool[porig] = p;
		} else {
		    block_t b2 = b;
		    b.len = length;

		    put_block(h, &b, 0, 1);

		    b2.len -= length;
		    if (b2.len >= MIN_LEN) {
			b2.pos += length;
			b2.free = 1;
			link_block(h, &b2);
			put_block(h, &b2, 0, 0);
		    }
		    if (!h->pool[porig])
			h->next_free_pool[porig] = pool(b2.len);
		    if (h->next_free_pool[porig] <= porig)
			h->next_free_pool[porig] = 0;
		}

#ifdef DEBUG
		fprintf(DEBUG_FP, "%ld\n", rover+4);
#endif
		h->next_free_time[porig] = h->timer;

		return rover + 4;
	    }

	    if (ms < b.len)
		ms = b.len;
	    
	    rover = b.next;
	} while (rover != orig);
	h->maxsz[p] = ms;

	/* If still here then none was found, we try next pool */
    }
    
    if (!h->pool[porig]) {
	h->next_free_pool[porig] = NPOOLS;
	h->next_free_time[porig] = h->timer;
    }

    /*
     * We've now searched all pools and none free. 
     * So instead we allocate off the wilderness block.
     */
    return wilderness_allocate(h, length);
}

/*
 * Frees the object at 'pos' from the disk heap h.
 * Returns 0 for success
 *        -1 for failure (shouldn't happen!)
 */
int heap_free(dheap_t *h, int64_t pos) {
    block_t b, pb, nb;

    //h->timer++;

    if (-1 == get_block(h, pos-4, &b))
	return -1;

#ifdef DEBUG
    fprintf(DEBUG_FP, "H-%p: heap_free %ld <= %ld\n", h, b.len, b.pos);
#endif

    if (b.pos + b.len == h->wilderness) {
	//h->wilderness = b.pos;
	//return 0;
	return link_block(h, &b);
    }

    if (-1 == get_block(h, b.pos + b.len, &nb)) {
	return -1;
    }

    assert(b.free == 0);

    if (!b.bfree && !nb.free) {
	/* [prev-free] ... [us] ... [next-free] */
	if (-1 == link_block(h, &b))
	    return -1;
	
    } else if (!nb.free) {
	/* [prev-free][us] ... [next-free] */
	if (-1 == get_block(h, b.pos - b.blen, &pb))
	    return -1;
	unlink_block(h, &pb);
	pb.len += b.len;
	link_block(h, &pb);

    } else if (!b.bfree) {
	/* [prev-free] ... [us][next-free] */
	unlink_block(h, &nb);
	b.len += nb.len;
	link_block(h, &b);

    } else {
	/* [prev-free][us][next-free] */
	unlink_block(h, &nb);
	if (-1 == get_block(h, b.pos - b.blen, &pb))
	    return -1;
	unlink_block(h, &pb);
	pb.len += b.len + nb.len;
	link_block(h, &pb);
    }

    return 0;
}

/*
 * Brute force check on the heap validity.
 * It's pretty inefficient and could be sped up, but for now this is just
 * a debugging tool rather than expecting to use it as a full integrity
 * checker.
 */
typedef struct free_bit_struct {
    uint64_t pos;
    uint64_t prev;
    uint64_t next;
    uint32_t len;
    struct free_bit_struct *n;
} free_bit;

int heap_largest_check(dheap_t *h) {
    free_bit *head = NULL, *fb = NULL;
    int largest = 0;

    /* Check the in-memory pool data matches on disk */
    {
	uint64_t p[NPOOLS];
	int i;
	lseek(h->fd, 0, SEEK_SET);
	read(h->fd, &p[0], 8*NPOOLS);
	for (i = 1; i < NPOOLS-1; i++) {
	    int pmin, pmax;

	    p[i] = be_int8(p[i]);
	    assert(p[i] == h->pool[i]);

	    if (i < 1) {
		pmin = 0;
		pmax = 16;
	    } else if (i < 126) {
		pmin = (i+2)<<3;
		pmax = (i+3)<<3;
	    } else {
		pmin = 1024-8 + (8 << (i-126));
		pmax = 1024-8 + (8 << (i+1-126));
	    }

	    pmax--;

	    assert(pool(pmin) == i);
	    assert(pool(pmax) == i);
	    assert(pool(pmin-1) == i-1);
	    assert(pool(pmax+1) == i+1);
	}
    }

    /* Now read each block in turn */
    {
	uint32_t len, len2;
	uint64_t prev, next;
	uint64_t offset = 8*NPOOLS;

	while (4 == read(h->fd, &len, 4)) {
	    read(h->fd, &prev, 8);
	    read(h->fd, &next, 8);

	    len = be_int4(len);
	    prev = be_int8(prev);
	    next = be_int8(next);

	    if (len & 1) {
		if (largest < (len & ~1))
		    largest = (len & ~1);
	    }

	    assert(len < 10000000);
	    assert((len & ~1) > 0);

	    if (len & 1) {
		fb = (free_bit *)calloc(1, sizeof(*fb));
		fb->pos = offset;
		fb->len = len;
		fb->prev = prev;
		fb->next = next;
		fb->n = head;
		head = fb;

		assert(fb->prev);
		assert(fb->next);
	    }

	    offset += len & ~1;
	    lseek(h->fd, offset-4, SEEK_SET);
	    read(h->fd, &len2, 4);
	    len2 = be_int4(len2);

	    if (len&1) {
		assert(len == len2);
	    }
	    assert((len&1) == (len2&1));
	}
    }

    /* Check all pools point to valid items of the correct size */
    {
	int i;

	/* Very inefficient, but ok for debugging */
	for (i = 0; i < NPOOLS; i++) {
	    int pmin, pmax;
	    uint64_t pos, more;
	    free_bit *last;

	    if (h->pool[i] == 0)
		continue;

	    if (i < 1) {
		pmin = 0;
		pmax = 16;
	    } else if (i < 126) {
		pmin = (i+2)<<3;
		pmax = (i+3)<<3;
	    } else {
		pmin = 1024-8 + (8 << (i-126));
		pmax = 1024-8 + (8 << (i+1-126));
	    }
	    pmax--;

	    pos = h->pool[i];
	    last = NULL;
	    more = 2;
	    do {
		/* Find free_bit struct */
		fb = head;
		while (fb) {
		    free_bit *next;
		    next = fb->n;
		    if (pos == fb->pos)
			break;
		    fb = next;
		}

		assert(fb);
		if (more == 2) {
		    assert(fb->len != 0);
		    assert(fb->len >= pmin);
		    assert(fb->len <= pmax);
		}
		fb->len = 0; /* Mark it as "been here" */

		if (fb->pos == h->pool[i])
		    more--;

		if (last)
		    assert(last->next == fb->pos);

		last = fb;
		pos = fb->next;
	    } while (more);
	}
    }

    /* Free, checking we've visited all too */
    fb = head;
    while (fb) {
	free_bit *next;
	next = fb->n;
	assert(fb->len == 0);
	free(fb);
	fb = next;
    }

    return largest;
}

void heap_check(dheap_t *h) {
    free_bit *head = NULL, *fb = NULL;

    puts("-- heap check --");
    /* Check the in-memory pool data matches on disk */
    {
	uint64_t p[NPOOLS];
	int i;
	lseek(h->fd, 0, SEEK_SET);
	read(h->fd, &p[0], 8*NPOOLS);
	for (i = 1; i < NPOOLS-1; i++) {
	    int pmin, pmax;

	    p[i] = be_int8(p[i]);
	    assert(p[i] == h->pool[i]);

	    if (i < 1) {
		pmin = 0;
		pmax = 16;
	    } else if (i < 126) {
		pmin = (i+2)<<3;
		pmax = (i+3)<<3;
	    } else {
		pmin = 1024-8 + (8 << (i-126));
		pmax = 1024-8 + (8 << (i+1-126));
	    }

	    pmax--;

	    if (p[i])
		printf(" pool(%d) = %"PRIu64"d (%d..%d)\n",
		       i, p[i], pmin, pmax);

	    assert(pool(pmin) == i);
	    assert(pool(pmax) == i);
	    assert(pool(pmin-1) == i-1);
	    assert(pool(pmax+1) == i+1);
	}
    }

    /* Now read each block in turn */
    {
	uint32_t len, len2;
	uint64_t prev, next;
	uint64_t offset = 8*NPOOLS;

	while (4 == read(h->fd, &len, 4)) {
	    read(h->fd, &prev, 8);
	    read(h->fd, &next, 8);

	    len = be_int4(len);
	    prev = be_int8(prev);
	    next = be_int8(next);

	    if (len & 1) {
		printf("%5"PRIu64"d+%4"PRIu32"d free prev=%5"PRIu64"d "
		       "next=%5"PRIu64"d\n",
		       offset, len & ~1, prev, next);
	    } else {
		printf("%5"PRIu64"d+%4"PRIu32"d used\n", offset, len & ~1);
	    }

	    assert(len < 10000000);
	    assert((len & ~1) > 0);

	    if (len & 1) {
		fb = (free_bit *)calloc(1, sizeof(*fb));
		fb->pos = offset;
		fb->len = len;
		fb->prev = prev;
		fb->next = next;
		fb->n = head;
		head = fb;

		assert(fb->prev);
		assert(fb->next);
	    }

	    offset += len & ~1;
	    lseek(h->fd, offset-4, SEEK_SET);
	    read(h->fd, &len2, 4);
	    len2 = be_int4(len2);

	    if (len&1) {
		assert(len == len2);
	    }
	    assert((len&1) == (len2&1));
	}
    }

    /* Check all pools point to valid items of the correct size */
    {
	int i;

	/* Very inefficient, but ok for debugging */
	for (i = 0; i < NPOOLS; i++) {
	    int pmin, pmax;
	    uint64_t pos, more;
	    free_bit *last;

	    if (h->pool[i] == 0)
		continue;

	    if (i < 1) {
		pmin = 0;
		pmax = 16;
	    } else if (i < 126) {
		pmin = (i+2)<<3;
		pmax = (i+3)<<3;
	    } else {
		pmin = 1024-8 + (8 << (i-126));
		pmax = 1024-8 + (8 << (i+1-126));
	    }
	    pmax--;

	    pos = h->pool[i];
	    last = NULL;
	    more = 2;
	    do {
		/* Find free_bit struct */
		fb = head;
		while (fb) {
		    free_bit *next;
		    next = fb->n;
		    if (pos == fb->pos)
			break;
		    fb = next;
		}

		assert(fb);
		if (more == 2) {
		    assert(fb->len != 0);
		    assert(fb->len >= pmin);
		    assert(fb->len <= pmax);
		}
		fb->len = 0; /* Mark it as "been here" */

		if (fb->pos == h->pool[i])
		    more--;

		if (last)
		    assert(last->next == fb->pos);

		last = fb;
		pos = fb->next;
	    } while (more);
	}
    }

    /* Free, checking we've visited all too */
    fb = head;
    while (fb) {
	free_bit *next;
	next = fb->n;
	assert(fb->len == 0);
	free(fb);
	fb = next;
    }
}

#ifdef HEAP_CHECKER
int main(int argc, char **argv) {
    dheap_t *h;

    if (argc != 2) {
	fprintf(stderr, "Usage: heap_check filename\n");
	return 1;
    }

    h = heap_load(argv[1], O_RDONLY);
    heap_check(h);

    return 0;
}
#endif

#ifdef TEST_MAIN
void write_data(int fd, uint64_t pos, int val) {
    char data[256];

    memset(data, val & 0xff, val & 0xff);
    lseek(fd, pos, SEEK_SET);
    write(fd, data, val & 0xff);
}

#define NP 50
int test1(int argc, char **argv) {
    dheap_t *h = heap_create("foo");
    int i;
    int64_t pos[NP];

    srand(atoi(argv[1]));
    memset(pos, 0, NP * sizeof(*pos));

    for (i = 0; i < 10000000; i++) {
	int p = rand() % NP;
	int j;

#if 0
	puts("");
	for (j = 0; j < NP; j++) {
	    if (pos[j])
		printf("%d:%ld\n", j, pos[j]);
	}
#endif

	if (rand()%5 != 0) {
	    int len;
	    int lenr;

	    if (rand()%10 >= 9)
		len = 1000 + rand() & 0xfff;
	    else if (rand()%10 >= 6)
		len = 100 + rand() & 0xff;
	    else
		len = 100 + rand() & 0x7f;

	    /* Alloc */
	    if (!pos[p]) {
		pos[p] = heap_allocate(h, len, &lenr);
		//printf("Alloc %4d => %-6ld\t(-4..+%d)\n", len, pos[p], lenr);
		//write_data(h->fd, pos[p], len);
		//heap_check(h);
	    }
	} else {
	    if (pos[p]) {
		/* Free */
		//printf("Free %ld\n", pos[p]);
		heap_free(h, pos[p]);
		pos[p] = 0;
		//heap_check(h);
	    }
	}

#if 0
	if (i % 1000 == 0) {
	    heap_check(h);
	}
#endif
    }

    heap_check(h);

    return 0;
}

int test2(int argc, char **argv) {
    dheap_t *h = heap_create("/tmp/foo");
    int i;
    int64_t pos;

    for (i = 0; i < 1000000; i++) {
	pos = heap_allocate(h, 10000000, NULL);
	printf("%8d: pos = %ld\n", i, (long)pos);
	assert(pos >= 0);
    }
    return 0;
}

int main(int argc, char **argv) {
    return test2(argc, argv);
}

#endif
