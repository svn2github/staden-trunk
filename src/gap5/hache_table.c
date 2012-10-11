#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

#include "os.h"
#include "hache_table.h"

/* #define DEBUG */

/* =========================================================================
 * TCL's hash function. Basically hash*9 + char.
 * =========================================================================
 */

uint32_t HacheTcl(uint8_t *data, int len) {
    uint32_t hash = 0;
    int i;

    for (i = 0; i < len; i++) {
	hash += (hash<<3) + data[i];
    }

    return hash;
}

/* =========================================================================
 * Paul Hsieh's hash function
 * http://www.azillionmonkeys.com/qed/hash.html
 * =========================================================================
 */

#undef get16bits
#if defined(__i386__) || defined(__i686__) || defined(__x86_64__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((const uint8_t *)(d))[1] << 8UL)\
                      +((const uint8_t *)(d))[0])
#endif

static uint32_t HacheHsieh(uint8_t *data, int len) {
    uint32_t hash = 0, tmp;
    int rem;

    if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
    case 3: hash += get16bits (data);
	hash ^= hash << 16;
	hash ^= data[sizeof (uint16_t)] << 18;
	hash += hash >> 11;
	break;
    case 2: hash += get16bits (data);
	hash ^= hash << 11;
	hash += hash >> 17;
	break;
    case 1: hash += *data;
	hash ^= hash << 10;
	hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
/*
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 2;
    hash += hash >> 15;
    hash ^= hash << 10;
*/
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}

static uint32_t HacheHsieh16(uint8_t *data) {
    uint32_t hash = 0, tmp;

    hash  += get16bits (data);
    tmp    = (get16bits (data+2) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    hash  += hash >> 11;

    hash  += get16bits (data+4);
    tmp    = (get16bits (data+6) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    hash  += hash >> 11;

    hash  += get16bits (data+8);
    tmp    = (get16bits (data+10) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    hash  += hash >> 11;

    hash  += get16bits (data+12);
    tmp    = (get16bits (data+14) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    hash  += hash >> 11;

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}

static uint32_t HacheHsieh8(uint8_t *data) {
    uint32_t hash = 0, tmp;

    hash  += get16bits (data);
    tmp    = (get16bits (data+2) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    hash  += hash >> 11;

    hash  += get16bits (data+4);
    tmp    = (get16bits (data+6) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    hash  += hash >> 11;

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}

/* =========================================================================
 * Bob Jenkins' hash function
 * http://burtleburtle.net/bob/hash/doobs.html
 * =========================================================================
 */

#define hashsize(n) ((uint32_t)1<<(n))
#define hashmask(n) (hashsize(n)-1)

/*
--------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.
For every delta with one or two bits set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() was built out of 36 single-cycle latency instructions in a 
  structure that could supported 2x parallelism, like so:
      a -= b; 
      a -= c; x = (c>>13);
      b -= c; a ^= x;
      b -= a; x = (a<<8);
      c -= a; b ^= x;
      c -= b; x = (b>>13);
      ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage 
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
--------------------------------------------------------------------
*/
#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

/*
--------------------------------------------------------------------
hash() -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  len     : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 6*len+35 instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (uint8_t **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);

By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

See http://burtleburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------
*/

static uint32_t HacheJenkins(uint8_t *k, int length /*, uint32_t initval */)
{
   register uint32_t a,b,c,len;

   /* Set up the internal state */
   len = length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = 0; /* initval; */        /* the previous hash value */

   /*---------------------------------------- handle most of the key */
   while (len >= 12)
   {
      a += (k[0] +((uint32_t)k[1]<<8) +((uint32_t)k[2]<<16) +((uint32_t)k[3]<<24));
      b += (k[4] +((uint32_t)k[5]<<8) +((uint32_t)k[6]<<16) +((uint32_t)k[7]<<24));
      c += (k[8] +((uint32_t)k[9]<<8) +((uint32_t)k[10]<<16)+((uint32_t)k[11]<<24));
      mix(a,b,c);
      k += 12; len -= 12;
   }

   /*------------------------------------- handle the last 11 bytes */
   c += length;
   switch(len)              /* all the case statements fall through */
   {
   case 11: c+=((uint32_t)k[10]<<24);
   case 10: c+=((uint32_t)k[9]<<16);
   case 9 : c+=((uint32_t)k[8]<<8);
      /* the first byte of c is reserved for the length */
   case 8 : b+=((uint32_t)k[7]<<24);
   case 7 : b+=((uint32_t)k[6]<<16);
   case 6 : b+=((uint32_t)k[5]<<8);
   case 5 : b+=k[4];
   case 4 : a+=((uint32_t)k[3]<<24);
   case 3 : a+=((uint32_t)k[2]<<16);
   case 2 : a+=((uint32_t)k[1]<<8);
   case 1 : a+=k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}

/*
 * An interface to the above hash functions.
 * Returns:
 *    A 32-bit hash key, suitable for masking down to smaller bit sizes
 */
uint32_t hache(int func, uint8_t *key, int key_len) {
    switch (func) {
    case HASH_FUNC_HSIEH:
	return HacheHsieh(key, key_len);

    case HASH_FUNC_TCL:
	return HacheTcl(key, key_len);
	
    case HASH_FUNC_JENKINS:
	return HacheJenkins(key, key_len);
	
    case HASH_FUNC_NULL:
    	return *(int *)key;
    }
    
    return 0;
}

/*
 * Special case of above for 16-byte keys
 */
static uint32_t hache16(int func, uint8_t *key) {
    switch (func) {
    case HASH_FUNC_HSIEH:
	return HacheHsieh16(key);

    case HASH_FUNC_TCL:
	return HacheTcl(key, 16);
	
    case HASH_FUNC_JENKINS:
	return HacheJenkins(key, 16);
	
    case HASH_FUNC_NULL:
    	return *(int *)key;
    }
    
    return 0;
}

static uint32_t hache8(int func, uint8_t *key) {
    switch (func) {
    case HASH_FUNC_HSIEH:
	return HacheHsieh8(key);

    case HASH_FUNC_TCL:
	return HacheTcl(key, 8);
	
    case HASH_FUNC_JENKINS:
	return HacheJenkins(key, 8);
	
    case HASH_FUNC_NULL:
    	return *(int *)key;
    }
    
    return 0;
}

/* =========================================================================
 * Hache Table handling code
 * =========================================================================
 */

static char *hname(HacheTable *h) {
    static char name[100];
    if (h->name)
	return h->name;

    sprintf(name, "%p", h);
    return name;
}


/* Multiplicative factors indicating when to grow or shrink the hash table */
#define HASH_TABLE_RESIZE 3

/*
 * Creates a HacheItem for use with HacheTable h.
 *
 * Returns:
 *    A pointer to new HacheItem on success
 *    NULL on failure.
 */
static HacheItem *HacheItemCreate(HacheTable *h) {
    HacheItem *hi;

    hi = (h->options & HASH_POOL_ITEMS ? 
    	pool_alloc(h->hi_pool) : malloc(sizeof(*hi)));

    if (NULL == hi) return NULL;

    hi->data.p    = NULL;
    hi->data.i    = 0;
    hi->next      = NULL;
    hi->key       = NULL;
    hi->key_len   = 0;
    hi->ref_count = 1;
    hi->order     = -1;
    hi->h         = h;
    hi->in_use_next = NULL;
    hi->in_use_prev = NULL;

    h->nused++;

    //printf("Hash %p item %p\n", h, hi);

    return hi;
}

/*
 * Deallocates a HacheItem created via HacheItemCreate.
 *
 * This function will not remove the item from the HacheTable so be sure to
 * call HacheTableDel() first if appropriate.
 */
static void HacheItemDestroy(HacheTable *h, HacheItem *hi, int deallocate_data) {
    assert(hi->h == h);

    if (!(h->options & HASH_NONVOLATILE_KEYS) || (h->options & HASH_OWN_KEYS))
	if (hi->key)
	    free(hi->key);

    if (deallocate_data) {
	if (h->del) {
	    h->del(h->clientdata, hi->data);
	} else if (hi->data.p) {
	    free(hi->data.p);
	}
    }

    if (hi->in_use_next)
	hi->in_use_next->in_use_prev = hi->in_use_prev;
    if (hi->in_use_prev)
	hi->in_use_prev->in_use_next = hi->in_use_next;
    if (h->in_use == hi)
	h->in_use = hi->in_use_next;

    
    if (h->options & HASH_POOL_ITEMS) 
    	pool_free(h->hi_pool, hi);
    else if (hi)
	free(hi);

    h->nused--;
}

/*
 * Creates a new HacheTable object. Size will be rounded up to the next
 * power of 2. It is a starting point and hash tables may be grown or shrunk
 * as needed (if HASH_DYNAMIC_SIZE is used).
 *
 * Options are as defined in the header file (see HASH_* macros).
 *
 * Returns:
 *    A pointer to a HacheTable on success
 *    NULL on failue
 */
HacheTable *HacheTableCreate(int size, int options) {
    HacheTable *h;
    int i, bits;
    uint32_t mask;
    int osize = size;

    if (!(h = (HacheTable *)malloc(sizeof(*h))))
	return NULL;

    if (options & HASH_POOL_ITEMS) {
        h->hi_pool = pool_create(sizeof(HacheItem));
	if (NULL == h->hi_pool) {
	    free(h);
	    return NULL;
	}
    } else {
        h->hi_pool = NULL;
    }

    if (size < 4)
	size = 4; /* an inconsequential minimum size */

    /* Round the requested size to the next power of 2 */
    bits = 0;
    size--;
    while (size) {
	size /= 2;
	bits++;
    }
    size = 1<<bits;
    mask = size-1;

    h->nbuckets = size;
    h->mask = mask;
    h->options = options;
    h->bucket = (HacheItem **)malloc(sizeof(*h->bucket) * size);
    h->nused = 0;
    h->searches = 0;
    h->hits = 0;
    h->cache_size = osize;
    //h->cache_size = 1000; /* Randomly picked default! */
    //h->cache_size = 2; /* Minimum; below this and foo=foo->next() breaks */
    //h->cache_size = 10; /* Small enough to test valid incr/decr */
    h->ordering = (HacheOrder *)malloc(h->cache_size * sizeof(*h->ordering));
    h->head = h->tail = -1;
    h->free = 0;
    for (i = 0; i < h->cache_size; i++) {
	h->ordering[i].hi = NULL;
	h->ordering[i].next = i == h->cache_size-1 ? -1 : i+1;
	h->ordering[i].prev = i-1;
    }
    h->clientdata = NULL;
    h->load = NULL;
    h->del = NULL;
    h->in_use = NULL;
    h->name = NULL;

    for (i = 0; i < size; i++) {
	h->bucket[i] = NULL;
    }

    return h;
}

/*
 * Deallocates a HacheTable object (created by HacheTableCreate).
 *
 * The deallocate_data parameter is a boolean to indicate whether the
 * data attached to the hash table should also be free()d. DO NOT USE
 * this if the HacheData attached was not a pointer allocated using
 * malloc().
 */
void HacheTableDestroy(HacheTable *h, int deallocate_data) {
    int i;

    if (!h)
	return;

    //HacheTableRefInfo(h, stdout);

    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi = h->bucket[i], *next = NULL;
	for (hi = h->bucket[i]; hi; hi = next) {
	    assert(hi->h == h);
	    next = hi->next;
	    HacheItemDestroy(h, hi, deallocate_data);
	}
    }
    
    if (h->hi_pool) pool_destroy(h->hi_pool);

    if (h->bucket)
	free(h->bucket);

    if (h->ordering)
	free(h->ordering);

    free(h);
}


/*
    Deletes all entries from the HacheTable while still leaving the
    HacheTable pointer.
*/

int HacheTableEmpty(HacheTable *h, int deallocate_data) {
    int i;
    
    if (!h) return -1;
    if (h->nbuckets == 0) return 0;
    
    // the destruction
    
    for (i = 0; i < h->nbuckets; i++) {
    	HacheItem *hi = h->bucket[i], *next = NULL;
	
	for (hi = h->bucket[i]; hi; hi = next) {
	    assert(hi->h == h);
	    next = hi->next;
	    HacheItemDestroy(h, hi, deallocate_data);
	}
    }
    
    if (h->bucket) free(h->bucket);
    
    if (h->ordering) free(h->ordering);
    
    // and a bit of creation

    if (h->hi_pool) {
    	pool_destroy(h->hi_pool);
	h->hi_pool = pool_create(sizeof(HacheItem));
	
	if (NULL == h->hi_pool) return -1;
    }
    
    // the creation proper
    h->bucket = (HacheItem **)malloc(sizeof(*h->bucket) * h->nbuckets);
    h->mask = h->nbuckets - 1;
    h->nused = 0;
    h->searches = 0;
    h->hits = 0;
    h->ordering = (HacheOrder *)malloc(h->cache_size * sizeof(*h->ordering));
    h->head = h->tail = -1;
    h->free = 0;

    for (i = 0; i < h->cache_size; i++) {
	h->ordering[i].hi = NULL;
	h->ordering[i].next = i == h->cache_size-1 ? -1 : i+1;
	h->ordering[i].prev = i-1;
    }

    h->clientdata = NULL;
    h->load = NULL;
    h->del = NULL;
    h->in_use = NULL;

    for (i = 0; i < h->nbuckets; i++) {
	h->bucket[i] = NULL;
    }
    
    return 0;
}



/*
 * Resizes a HacheTable to have 'newsize' buckets.
 * This is called automatically when adding or removing items so that the
 * hash table keeps at a sensible scale.
 *
 * FIXME: Halving the size of the hash table is simply a matter of coaelescing
 * every other bucket. Instead we currently rehash (which is slower).
 * Doubling the size of the hash table currently requires rehashing, but this
 * too could be optimised by storing the full 32-bit hash of the key along
 * with the key itself. This then means that it's just a matter of seeing what
 * the next significant bit is. It's a memory vs speed tradeoff though and
 * re-hashing is pretty quick.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int HacheTableResize(HacheTable *h, int newsize) {
    HacheTable *h2;
    int i;

#ifdef DEBUG
    fprintf(stdout, "Resizing HacheTable %s to %d\n", hname(h), newsize);
#endif

    /* Create a new hash table and rehash everything into it */
    h2 = HacheTableCreate(newsize, h->options);
    
    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi, *next;
	for (hi = h->bucket[i]; hi; hi = next) {
	    assert(hi->h == h);
	    uint32_t hv = hache(h2->options & HASH_FUNC_MASK,
			       (uint8_t *)hi->key, hi->key_len) & h2->mask;
	    next = hi->next;
	    hi->next = h2->bucket[hv];
	    h2->bucket[hv] = hi;
	}
    }

    /* Swap the links over & free */
    free(h->bucket);
    h->bucket   = h2->bucket;
    h->nbuckets = h2->nbuckets;
    h->mask     = h2->mask;
    if (h2->ordering)
	free(h2->ordering);
    free(h2);
    
    return 0;
}

/*
 * Expands the cache size in a hash table. This is only typically required due
 * to too many items being marked as in-use.
 */
int HacheTableExpandCache(HacheTable *h) {
    HacheOrder *newo;
    int i, j = h->cache_size;

    fprintf(stderr, "Cache order %s full, doubling size (%d)!\n",
	    hname(h), h->cache_size*2);

    /* Double size */
    newo = (HacheOrder *)realloc(h->ordering, h->cache_size*2 * sizeof(*newo));
    if (NULL == newo)
	return -1;

    h->cache_size *= 2;
    h->ordering = newo;

    /* Add free entries */
    for (i = j; i < h->cache_size; i++) {
	h->ordering[i].hi = NULL;
	h->ordering[i].next = i == h->cache_size-1 ? -1 : i+1;
	h->ordering[i].prev = i-1;
    }
    if (-1 != h->free) {
	h->ordering[h->cache_size-1].next = h->free;
	h->ordering[h->free].prev = h->cache_size-1;
    }

    h->ordering[j].prev = -1;
    h->free = j;

    return 0;
}

/*
 * Moves a HacheItem in the ordering array to the end of the array, updating
 * head/tail appropriately.
 */
void HacheOrderAccess(HacheTable *h, HacheItem *hi) {
    int i = hi->order;

    assert(hi->h == h);

    if (i == -1)
	return; /* presumably this is a locked item */

    if (h->tail == i)
	return; /* nothing to do */

    /* Remove links to us */
    if (h->ordering[i].next != -1)
	h->ordering[h->ordering[i].next].prev = h->ordering[i].prev;
    if (h->ordering[i].prev != -1)
	h->ordering[h->ordering[i].prev].next = h->ordering[i].next;

    /* Fix head */
    if (h->head == i)
	h->head = h->ordering[i].next;

    /* Add to end */
    h->ordering[h->tail].next = i;
    h->ordering[i].next = -1;
    h->ordering[i].prev = h->tail;

    /* Fix tail */
    h->tail = i;
}

/*
 * Removes a HacheItem from the ordering array
 */
void HacheOrderRemove(HacheTable *h, HacheItem *hi) {
    int i = hi->order;
    
    assert(hi->h == h);

    if (hi->order == -1)
	return;

    /* Unlink */
    if (h->ordering[i].next != -1)
	h->ordering[h->ordering[i].next].prev = h->ordering[i].prev;
    if (h->ordering[i].prev != -1)
	h->ordering[h->ordering[i].prev].next = h->ordering[i].next;

    /* Update head/tails */
    if (h->head == i)
	h->head = h->ordering[i].next;
    if (h->tail == i)
	h->tail = h->ordering[i].prev;

    /* Add to free list */
    h->ordering[i].hi = NULL;
    h->ordering[i].prev = -1;
    h->ordering[i].next = h->free;
    h->free = i;
}

/*
 * Adds a HacheItem to the ordering array, freeing an old item if required.
 * Returns the order index of the newly added item.
 */
int HacheOrderAdd(HacheTable *h, HacheItem *hi) {
    int i;

    assert(hi->h == h);

    /* Free the oldest unused item if it's full */
    if (h->free == -1) {
	if (h->head != -1) {
	    HacheTableDel(h, h->ordering[h->head].hi, 1);
	}
    }

    if (h->free == -1) {
	if (-1 == HacheTableExpandCache(h)) {
	    fprintf(stderr, "Failed to expand\n");
	    return -1;
	}
    }

    /* Add the latest item to the end */
    i = h->free;
    h->free = h->ordering[i].next;
    if (-1 != h->free)
	h->ordering[h->free].prev = -1;

    h->ordering[i].hi = hi;
    h->ordering[i].prev = h->tail;
    h->ordering[i].next = -1;

    if (-1 != h->tail)
	h->ordering[h->tail].next = i;
    h->tail = i;

    if (-1 == h->head)
	h->head = i;

    /*    
    printf("hi=%p %p->ordering[%d]={%p,%d,%d}\n",
	   hi, h, i,
	   h->ordering[i].hi,
	   h->ordering[i].next,
	   h->ordering[i].prev);
    */

    return i;
}

/* Purges the LRU cache to check for cache coherency bugs */
void HacheOrderPurge(HacheTable *h) {
    int i, in;
    for (i = h->head; i != -1; i = in) {
	in = h->ordering[i].next;
	HacheTableDel(h, h->ordering[i].hi, 1);
    }
}

void HacheTableIncRef(HacheTable *h, HacheItem *hi) {
    assert(hi->h == h);

    hi->ref_count++;
    if (hi->order != -1) {
	HacheOrderRemove(h, hi);
	hi->order = -1;
    }

    if (!(h->in_use == hi || hi->in_use_prev || hi->in_use_next)) {
	hi->in_use_next = h->in_use;
	if (h->in_use)
	    h->in_use->in_use_prev = hi;
	hi->in_use_prev = NULL;
	h->in_use = hi;
    }
}

void HacheTableDecRef(HacheTable *h, HacheItem *hi) {
    assert(hi->h == h);

    if (hi->ref_count <= 0) {
	fprintf(stderr, "WARNING: attempting to decrement reference count "
		"on hache item %p when ref_count is already <= 0\n", hi);
    }

    if (hi && hi->ref_count > 0) {
	if (--hi->ref_count <= 0) {
	    hi->order = HacheOrderAdd(h, hi);

	    if (hi->in_use_next)
		hi->in_use_next->in_use_prev = hi->in_use_prev;
	    if (hi->in_use_prev)
		hi->in_use_prev->in_use_next = hi->in_use_next;
	    if (h->in_use == hi)
		h->in_use = hi->in_use_next;
	    hi->in_use_next = NULL;
	    hi->in_use_prev = NULL;
	}
    }
}

/*
 * Adds a HacheData item to HacheTable h with a specific key. Key can be binary
 * data, but if key_len is passed as zero then strlen() will be used to
 * determine the key length.
 *
 * The "new" pointer may be passed as NULL. When not NULL it is filled out
 * as a boolean to indicate whether the key is already in this hash table.
 *
 * The HASH_ALLOW_DUP_KEYS option (specified when using HacheTableCreate)
 * will allow duplicate keys to be stored, and hence *new is also zero.
 * By default duplicate keys are disallowed.
 *
 * Keys are considered to be volatile memory (ie temporary storage) and so the
 * hash table takes separate copies of them. To avoid this use the
 * HASH_NONVOLATILE_KEYS option.
 *
 * If the HASH_OWN_KEYS option was specified when creating the table then
 * keys will be considered to be owned by the hash table. In this case
 * the key will be freed when the table is destroyed regardless of
 * whether the HASH_NONVOLATILE_KEYS option was used to allocate its
 * own private copy.
 *
 * Returns:
 *    The HacheItem created (or matching if a duplicate) on success
 *    NULL on failure
 */
HacheItem *HacheTableAdd(HacheTable *h, char *key, int key_len, HacheData data,
			 int *new) {
    uint32_t hv;
    HacheItem *hi;

#ifdef DEBUG
    printf("HacheTableAdd %s %d data.p %p\n", hname(h), *(int *)key, data.p);
#endif

    if (!key_len)
	key_len = strlen(key);

    hv = hache(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len) & h->mask;
    
    /* Already exists? */
    if (!(h->options & HASH_ALLOW_DUP_KEYS)) {
	for (hi = h->bucket[hv]; hi; hi = hi->next) {
	    if (key_len == hi->key_len &&
		memcmp(key, hi->key, key_len) == 0) {
		if (new) *new = 0;
		return hi;
	    }
	}
    }

    /* No, so create a new one and link it in */
    if (NULL == (hi = HacheItemCreate(h)))
	return NULL;

    if (h->options & HASH_NONVOLATILE_KEYS)
	hi->key = key;
    else {
	hi->key = (char *)malloc(key_len+1);
	memcpy(hi->key, key, key_len);
	hi->key[key_len] = 0; /* null terminate incase others print keys */
    }
    hi->key_len = key_len;
    hi->data = data;
    hi->next = h->bucket[hv];
    h->bucket[hv] = hi;

    if ((h->options & HASH_DYNAMIC_SIZE) &&
	h->nused > HASH_TABLE_RESIZE * h->nbuckets)
	HacheTableResize(h, h->nbuckets*4);

    if (new) *new = 1;

    return hi;
}

/*
 * Given a new key this rehashes an existing item in the hash table.
 * This is unlikely to fail except in low memory conditions or if the
 * new key already exists and we have disallowed duplicate keys.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int HacheTableRehash(HacheTable *h, HacheItem *hi, char *key, int key_len) {
    uint32_t hv, hv_old;
    HacheItem *ti, *last, *next;

    assert(hi->h == h);

    hv = hache(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len) & h->mask;
    hv_old = hache(h->options & HASH_FUNC_MASK,
		   (uint8_t *)hi->key, hi->key_len) & h->mask;

    /* Already exists? */
    if (!(h->options & HASH_ALLOW_DUP_KEYS)) {
	for (ti = h->bucket[hv]; ti; ti = ti->next) {
	    if (key_len == ti->key_len &&
		memcmp(key, ti->key, key_len) == 0) {
		return -1;
	    }
	}
    }    

    /* Re-key the item */
    if (h->options & HASH_NONVOLATILE_KEYS) {
	hi->key = key;
    } else {
	char *cp = (char *)malloc(key_len+1);
	if (!cp)
	    return -1;
	free(hi->key);
	hi->key = cp;
	memcpy(hi->key, key, key_len);
	hi->key[key_len] = 0; /* null terminate incase others print keys */
    }
    hi->key_len = key_len;

    /* Remove from old loc */
    for (last = NULL, next = h->bucket[hv_old]; next;
	 last = next, next = next->next) {
	if (next == hi) {
	    /* Link last to hi->next */
	    if (last)
		last->next = hi->next;
	    else
		h->bucket[hv_old] = hi->next;
	}
    }

    /* Add to new loc */
    hi->next = h->bucket[hv];
    h->bucket[hv] = hi;

    return 0;
}


/*
 * Removes a specified HacheItem from the HacheTable. (To perform this it needs
 * to rehash based on the hash key as hash_item only has a next pointer and
 * not a previous pointer.)
 * 
 * The HacheItem itself is also destroyed (by an internal call to
 * HacheItemDestroy). The deallocate_data parameter controls whether the data
 * associated with the HacheItem should also be free()d.
 *
 * See also the HacheTableRemove() function to remove by key instead of
 * HacheItem.
 *
 * Returns 0 on success
 *        -1 on failure (eg HacheItem not in the HacheTable);
 */
int HacheTableDel(HacheTable *h, HacheItem *hi, int deallocate_data) {
    uint32_t hv;
    HacheItem *next, *last;

    assert(hi->h == h);

#ifdef DEBUG
    printf("HacheTableDel %s %p\n", hname(h), hi->data.p);
#endif

    hv = hache(h->options & HASH_FUNC_MASK,
	      (uint8_t *)hi->key, hi->key_len) & h->mask;

    for (last = NULL, next = h->bucket[hv]; next;
	 last = next, next = next->next) {
	if (next == hi) {
	    /* Link last to hi->next */
	    if (last)
		last->next = hi->next;
	    else
		h->bucket[hv] = hi->next;

	    HacheOrderRemove(h, hi);
	    HacheItemDestroy(h, hi, deallocate_data);

	    return 0;
	}
    }

    return -1;
}


/*
 * Searches the HacheTable for the data registered with 'key' and removes
 * these items from the HacheTable. In essence this is a combination of
 * HacheTableSearch and HacheTableDel functions.
 *
 * If HASH_ALLOW_DUP_KEYS is used this will remove all items matching 'key',
 * otherwise just a single item will be removed.
 *
 * If 'deallocate_data' is true the data associated with the HacheItem will
 * be free()d.
 *
 * Returns
 *    0 on success (at least one item found)
 *   -1 on failure (no items found).
 */
int HacheTableRemove(HacheTable *h, char *key, int key_len,
		     int deallocate_data) {
    uint32_t hv;
    HacheItem *last, *next, *hi;
    int retval = -1;

    if (!key_len)
	key_len = strlen(key);

    hv = hache(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len) & h->mask;

    last = NULL;
    next = h->bucket[hv];

    while (next) {
	hi = next;
	if (key_len == hi->key_len &&
	    memcmp(key, hi->key, key_len) == 0) {
	    /* An item to remove, adjust links and destroy */
	    if (last)
		last->next = hi->next;
	    else
		h->bucket[hv] = hi->next;

	    next = hi->next;
	    HacheOrderRemove(h, hi);
	    HacheItemDestroy(h, hi, deallocate_data);

	    retval = 0;
	    if (!(h->options & HASH_ALLOW_DUP_KEYS))
		break;

	} else {
	    /* We only update last when it's something we haven't destroyed */
	    last = hi;
	    next = hi->next;
	}
    }

    return retval;
}

/*
 * As HacheTableSearch below, but this won't force a load or update access
 * patterns. It's simply a means to check if an item exists. One intended
 * use is to query and item and if present decrement the reference.
 */
HacheItem *HacheTableQuery(HacheTable *h, char *key, int key_len) {
    uint32_t hv;
    HacheItem *hi;

#ifdef DEBUG
    printf("HacheTableQuery %s %d\n", hname(h), *(int *)key);
#endif

    h->searches++;

    if (!key_len)
	key_len = strlen(key);

    /* Return if present */
    if (key_len == 16) {
	/* cache_search key size is most common */
	hv = hache16(h->options & HASH_FUNC_MASK,
		     (uint8_t *)key) & h->mask;
	for (hi = h->bucket[hv]; hi; hi = hi->next) {
	    if (hi->key_len == 16) {
		uint32_t *a = (uint32_t *)key;
		uint32_t *b = (uint32_t *)hi->key;
		int m =
		    (a[0] - b[0] == 0) &&
		    (a[1] - b[1] == 0) &&
		    (a[2] - b[2] == 0) &&
		    (a[3] - b[3] == 0);

		if (m) {
		    h->hits++;
		    HacheOrderAccess(h, hi);
#ifdef DEBUG
		    printf("\thit = %p\n", hi->data.p);
#endif
		    return hi;
		}
	    }
	}
    }

    hv = hache(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len) & h->mask;
    for (hi = h->bucket[hv]; hi; hi = hi->next) {
	if (key_len == hi->key_len &&
	    memcmp(key, hi->key, key_len) == 0) {
	    h->hits++;
	    HacheOrderAccess(h, hi);
#ifdef DEBUG
	    printf("\thit = %p\n", hi->data.p);
#endif
	    return hi;
	}
    }

    /* Unable to load, so it's really not present */
    return NULL;
}

/*
 * Searches the HacheTable for the data registered with 'key'.
 * If HASH_ALLOW_DUP_KEYS is used this will just be the first one found.
 * You will then need to use HacheTableNext to iterate through the matches.
 *
 * Returns
 *    HacheItem if found
 *    NULL if not found
 */
HacheItem *HacheTableSearch(HacheTable *h, char *key, int key_len) {
    uint32_t hv;
    HacheItem *hi;

#ifdef DEBUG
    printf("HacheTableSearch %s %d\n", hname(h), *(int *)key);
#endif

    h->searches++;

    if (!key_len)
	key_len = strlen(key);

    /* Return if present */
    hv = hache(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len) & h->mask;
    for (hi = h->bucket[hv]; hi; hi = hi->next) {
	if (key_len == hi->key_len &&
	    memcmp(key, hi->key, key_len) == 0) {
	    h->hits++;
	    HacheOrderAccess(h, hi);
#ifdef DEBUG
	    printf("\thit = %p\n", hi->data.p);
#endif
	    return hi;
	}
    }

    /* Attempt a load & store if not */
    if (h->load) {
	HacheData dummy, *hd;
	dummy.p = NULL;

	/* Allocate storage in hash table */
	hi = HacheTableAdd(h, key, key_len, dummy, NULL);
	if (!hi)
	    return NULL;

	/* Load the object and link to store */
	hd = h->load(h->clientdata, key, key_len, hi);
	if (hd) {
	    hi->data = *hd;
	    return hi;
	} else {
	    HacheTableDel(h, hi, 0);
	    return NULL;
	}
    }

    /* Unable to load, so it's really not present */
    return NULL;
}

/*
 * Find the next HacheItem (starting from 'hi') to also match this key.
 * This is only valid when the HASH_ALLOW_DUP_KEYS is in use.
 *
 * Returns
 *    HacheItem if found
 *    NULL if not found
 */
HacheItem *HacheTableNext(HacheItem *hi, char *key, int key_len) {
    if (!hi)
	return NULL;

    for (hi = hi->next; hi; hi = hi->next) {
	if (key_len == hi->key_len &&
	    memcmp(key, hi->key, key_len) == 0) {
	    return hi;
	}
    }

    return NULL;
}

/*
 * Reverses the order of items in a buckets of the hash table.
 * HacheTableAdd has the effect of reversing the list as it always adds to the
 * start (for efficiency), but sometimes we wish to keep this order consistent.
 */
void HacheTableReverse(HacheTable *h) {
    int i;
    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi = h->bucket[i], *last = NULL, *next = NULL;

	if (!hi)
	    continue;

	while (hi) {
	    next = hi->next;
	    hi->next = last;
	    last = hi;
	    hi = next;
	}

	h->bucket[i] = last;
    }
}

/*
 * Dumps a textual represenation of the hash table to stdout.
 */
void HacheTableDump(HacheTable *h, FILE *fp) {
    int i;
    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    int j, printable=1;
	    for (j = 0; j < hi->key_len; j++) {
		if (!isprint(hi->key[j])) {
		    printable=0;
		    break;
		}
	    }

	    if (printable) {
		fprintf(fp, "%.*s\n", hi->key_len, hi->key);
	    } else {
		if (hi->key_len == 4) {
		    fprintf(fp, "%d\n", *(int *)(hi->key));
		} else {
		    fprintf(fp, "%p ", hi->key);
		    for (j = 0; j < hi->key_len; j++)
			fprintf(fp, "%02x ", (unsigned char)(hi->key[j]));
		    fprintf(fp, "\n");
		}
	    }
	}
    }
}

/*
 * Reports how much data in the hache table is in use (reference count > 0)
 * and how much is just cached data to be freed as and when needed.
 */
void HacheTableRefInfo(HacheTable *h, FILE *fp) {
    int nr = 0, nu = 0, no = 0, nf = 0;

    if (!fp) fp = stdout;
    
    int i;
    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    if (hi->ref_count)
		nr++;
	    else
		nu++;
	    if (hi->order != -1)
		no++;
	}
    }

    if (h->cache_size) {
	for (i = h->free; i != -1; i = h->ordering[i].next)
	    nf++;
    }

    fprintf(fp, "Hache Table %s\n", hname(h));
    fprintf(fp, "    Cache size       %d\n", h->cache_size);
    fprintf(fp, "    Refcount > 0     %d\n", nr);
    fprintf(fp, "    Refcount = 0     %d\n", nu);
    fprintf(fp, "    Items with order %d\n", no);
    fprintf(fp, "    Items to reuse   %d\n", nf);
    assert(no + nf == h->cache_size);
    assert(no == nu);
}

/*
 * Produces some simple statistics on the hash table population.
 */
void HacheTableStats(HacheTable *h, FILE *fp) {
    int i;
    double avg = (double)h->nused / h->nbuckets;
    double var = 0;
    int maxlen = 0;
    int filled = 0;
    int clen[51];
    int count1, count2;

    if (!fp) fp = stdout;

    for (i = 0; i <= 50; i++)
	clen[i] = 0;

    for (i = 0; i < h->nbuckets; i++) {
	int len = 0;
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    assert(hi->h == h);
	    len++;
	}
	if (len > 0) {
	    filled++;
	    if (len > maxlen)
		maxlen = len;
	}
	clen[len <= 50 ? len : 50]++;
	var += (len-avg) * (len-avg);
    }
    var /= h->nbuckets;
    /* sd = sqrt(var); */

    fprintf(fp, "Nbuckets  = %u\n", h->nbuckets);
    fprintf(fp, "Nused     = %d\n", h->nused);
    fprintf(fp, "Avg chain = %f\n", avg);
    fprintf(fp, "Chain var.= %f\n", var);
    fprintf(fp, "%%age full = %f\n", (100.0*filled)/h->nbuckets);
    fprintf(fp, "max len   = %d\n", maxlen);
    fprintf(fp, "cache_size= %d\n", h->cache_size);

    for (count1 = count2 = i = 0; i < h->cache_size; i++) {
	if (h->ordering[i].hi) {
	    count1++;
	    if (h->ordering[i].hi->ref_count)
		count2++;
	}
    }
    fprintf(fp, "N.cached  = %d\n", count1);
    fprintf(fp, "N.locked  = %d\n", count2);
    fprintf(fp, "Searches  = %d\n", h->searches);
    fprintf(fp, "Cache hits= %d (%6.2f)%%\n",
	    h->hits, 100.0*h->hits/h->searches);

    /* Reset counters for next stats call */
    h->hits = h->searches = 0; 

    for (i = 0; i <= maxlen; i++) {
	fprintf(fp, "Chain %2d   = %d\n", i, clen[i]);
    }

    //HacheTableLeakCheck(h);
}

/*
 * For debugging purposes only. This function severs all links to items
 * with a reference count > 0. Obviously this breaks the HacheTable in
 * various ways, but the purpose is that it should be immediately followed
 * up by an exit() call and an analysis of memory leaks.
 *
 * In theory all items with a reference count > 0 will have pointers to them
 * in other pieces of code. If they do not then we know we leaked memory
 * somewhere by virtue of incrementing the reference count and not
 * decrementing it again before losing our pointer to the object.
 */
void HacheTableLeakCheck(HacheTable *h, FILE *fp) {
    int i;
    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi, *next, *last = NULL;
	for (hi = h->bucket[i]; hi; last = hi, hi = next) {
	    assert(hi->h == h);
	    next = hi->next;

	    if (!hi->ref_count)
		continue;

	    //	    printf("Has ref count %d: %.*s\n",
	    //		   hi->ref_count, hi->key_len, hi->key);

	    /* Remove all memory links to haches, keys and data */
	    if (last)
		last->next = next;
	    else
		h->bucket[i] = next;
	    hi->next = NULL;
	    hi->h = NULL;
	    hi->key = NULL;
	    hi->data.p = 0;
	    if (hi->in_use_next) {
		hi->in_use_next->in_use_prev = NULL;
		hi->in_use_next = NULL;
	    }
	    if (hi->in_use_prev) {
		hi->in_use_prev->in_use_next = NULL;
		hi->in_use_prev = NULL;
	    }
	}
    }
}

/*
 * Iterates through members of a hash table returning items sequentially.
 *
 * Returns the next HacheItem on success
 *         NULL on failure.
 */
HacheItem *HacheTableIterNext(HacheTable *h, HacheIter *iter) {
    do {
	if (iter->hi == NULL) {
	    if (++iter->bnum >= h->nbuckets)
		break;
	    iter->hi = h->bucket[iter->bnum];
	} else {
	    iter->hi = iter->hi->next;
	}
    } while (!iter->hi);
    
    return iter->hi;
}

void HacheTableIterReset(HacheIter *iter) {
    if (iter) {
	iter->bnum = -1;
	iter->hi = NULL;
    }
}

HacheIter *HacheTableIterCreate(void) {
    HacheIter *iter = (HacheIter *)malloc(sizeof(*iter));

    HacheTableIterReset(iter);
    return iter;
}

void HacheTableIterDestroy(HacheIter *iter) {
    if (iter)
	free(iter);
}
