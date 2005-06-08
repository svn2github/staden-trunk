#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash_table.h"
#include "os.h"

/*
 * TO DO
 *
 * - HashItemRemove
 *   Remove single or remove all if HASH_ALLOW_DUP_KEYS is used.
 */


/* =========================================================================
 * TCL's hash function. Basically hash*9 + char.
 * =========================================================================
 */

uint32_t HashTcl(uint8_t *data, int len) {
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
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((const uint8_t *)(d))[1] << 8UL)\
                      +((const uint8_t *)(d))[0])
#endif

uint32_t HashHsieh(uint8_t *data, int len) {
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
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 2;
    hash += hash >> 15;
    hash ^= hash << 10;

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

uint32_t HashJenkins(uint8_t *k, int length /*, uint32_t initval */)
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
uint32_t hash(int func, uint8_t *key, int key_len) {
    switch (func) {
    case HASH_FUNC_HSIEH:
	return HashHsieh(key, key_len);

    case HASH_FUNC_TCL:
	return HashTcl(key, key_len);
	
    case HASH_FUNC_JENKINS:
	return HashJenkins(key, key_len);
    }
    
    return 0;
}

/* =========================================================================
 * Hash Table handling code
 * =========================================================================
 */

/* Multiplicative factors indicating when to grow or shrink the hash table */
#define HASH_TABLE_RESIZE 3

/*
 * Creates a HashItem for use with HashTable h.
 *
 * Returns:
 *    A pointer to new HashItem on success
 *    NULL on failure.
 */
static HashItem *HashItemCreate(HashTable *h) {
    HashItem *hi;

    if (!(hi = (HashItem *)malloc(sizeof(*hi))))
	return NULL;

    hi->data.p    = NULL;
    hi->data.i    = 0;
    hi->next      = NULL;
    hi->key       = NULL;
    hi->key_len   = 0;

    h->nused++;
    
    return hi;
}

/*
 * Deallocates a HashItem created via HashItemCreate.
 */
static void HashItemDestroy(HashTable *h, HashItem *hi) {
    if (!(h->options & HASH_NONVOLATILE_KEYS) && hi->key)
	free(hi->key);

    if (hi)
	free(hi);

    h->nused--;
}

/*
 * Creates a new HashTable object. Size will be rounded up to the next
 * power of 2. It is a starting point and hash tables may be grown or shrunk
 * as needed (if HASH_DYNAMIC_SIZE is used).
 *
 * Options are as defined in the header file (see HASH_* macros).
 *
 * Returns:
 *    A pointer to a HashTable on success
 *    NULL on failue
 */
HashTable *HashTableCreate(int size, int options) {
    HashTable *h;
    int i, bits;
    uint32_t mask;

    if (!(h = (HashTable *)malloc(sizeof(*h))))
	return NULL;

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
    h->bucket = (HashItem **)malloc(sizeof(*h->bucket) * size);
    h->nused = 0;

    for (i = 0; i < size; i++) {
	h->bucket[i] = NULL;
    }

    return h;
}

/*
 * Deallocates a HashTable object (created by HashTableCreate).
 */
void HashTableDestroy(HashTable *h) {
    int i;

    if (!h)
	return;

    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi = h->bucket[i], *next = NULL;
	for (hi = h->bucket[i]; hi; hi = next) {
	    next = hi->next;
	    HashItemDestroy(h, hi);
	}
    }

    if (h->bucket)
	free(h->bucket);

    free(h);
}

/*
 * Resizes a HashTable to have 'newsize' buckets.
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
int HashTableResize(HashTable *h, int newsize) {
    HashTable *h2;
    int i;

    fprintf(stderr, "Resizing to %d\n", newsize);

    /* Create a new hash table and rehash everything into it */
    h2 = HashTableCreate(newsize, h->options);
    
    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi, *next;
	for (hi = h->bucket[i]; hi; hi = next) {
	    uint32_t hv = hash(h2->options & HASH_FUNC_MASK,
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
    free(h2);

    return 0;
}

/*
 * Adds a HashData item to HashTable h with a specific key. Key can be binary
 * data, but if key_len is passed as zero then strlen() will be used to
 * determine the key length.
 *
 * The "new" pointer may be passed as NULL. When not NULL it is filled out
 * as a boolean to indicate whether the key is already in this hash table.
 *
 * The HASH_ALLOW_DUP_KEYS option (specified when using HashTableCreate)
 * will allow duplicate keys to be stored, and hence *new is also zero.
 * By default duplicate keys are disallowed.
 *
 * Keys are considered to be volatile memory (ie temporary storage) and so the
 * hash table takes separate copies of them. To avoid this use the
 * HASH_NONVOLATILE_KEYS option.
 *
 * Returns:
 *    The HashItem created (or matching if a duplicate) on success
 *    NULL on failure
 */
HashItem *HashTableAdd(HashTable *h, char *key, int key_len, HashData data,
		       int *new) {
    uint32_t hv;
    HashItem *hi;

    if (!key_len)
	key_len = strlen(key);

    hv = hash(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len) & h->mask;

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
    if (NULL == (hi = HashItemCreate(h)))
	return NULL;

    if (h->options & HASH_NONVOLATILE_KEYS)
	hi->key = key;
    else {
	hi->key = strdup(key);
    }
    hi->key_len = key_len;
    hi->data = data;
    hi->next = h->bucket[hv];
    h->bucket[hv] = hi;

    if ((h->options & HASH_DYNAMIC_SIZE) &&
	h->nused > HASH_TABLE_RESIZE * h->nbuckets)
	HashTableResize(h, h->nbuckets*4);

    if (new) *new = 1;

    return hi;
}

/*
 * Searches the HashTable for the data registered with 'key'.
 * If HASH_ALLOW_DUP_KEYS is used this will just be the first one found.
 * You will then need to use HashTableNext to iterate through the matches.
 *
 * Returns
 *    HashItem if found
 *    NULL if not found
 */
HashItem *HashTableSearch(HashTable *h, char *key, int key_len) {
    uint32_t hv;
    HashItem *hi;

    if (!key_len)
	key_len = strlen(key);

    hv = hash(h->options & HASH_FUNC_MASK, (uint8_t *)key, key_len) & h->mask;
    for (hi = h->bucket[hv]; hi; hi = hi->next) {
	if (key_len == hi->key_len &&
	    memcmp(key, hi->key, key_len) == 0)
	    return hi;
    }

    return NULL;
}

/*
 * Find the next HashItem (starting from 'hi') to also match this key.
 * This is only valid when the HASH_ALLOW_DUP_KEYS is in use.
 *
 * Returns
 *    HashItem if found
 *    NULL if not found
 */
HashItem *HashTableNext(HashItem *hi, char *key, int key_len) {
    if (!hi)
	return NULL;

    for (hi = hi->next; hi; hi = hi->next) {
	if (key_len == hi->key_len &&
	    memcmp(key, hi->key, key_len) == 0)
	    return hi;
    }

    return NULL;
}

/*
 * Dumps a textual represenation of the hash table to stdout.
 */
void HashTableDump(HashTable *h, FILE *fp) {
    int i;
    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    fprintf(fp, "%.*s %"PRId64" %"PRId32"\n",
		    hi->key_len, hi->key,
		    hi->data.idx.pos, hi->data.idx.len);
	}
    }
}

/*
 * Produces some simple statistics on the hash table population.
 */
void HashTableStats(HashTable *h, FILE *fp) {
    int i;
    double avg = (double)h->nused / h->nbuckets;
    double var = 0;
    int maxlen = 0;
    int filled = 0;
    int clen[51];

    for (i = 0; i <= 50; i++)
	clen[i] = 0;

    for (i = 0; i < h->nbuckets; i++) {
	int len = 0;
	HashItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
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

    fprintf(fp, "Nbuckets  = %d\n", h->nbuckets);
    fprintf(fp, "Nused     = %d\n", h->nused);
    fprintf(fp, "Avg chain = %f\n", avg);
    fprintf(fp, "Chain var.= %f\n", var);
    fprintf(fp, "%%age full = %f\n", (100.0*filled)/h->nbuckets);
    fprintf(fp, "max len   = %d\n", maxlen);
    for (i = 0; i <= maxlen; i++) {
	fprintf(fp, "Chain %2d   = %d\n", i, clen[i]);
    }
}

/*
 * --------------------------------------------------------------------
 * Below we have a specialisation of the HashTable code where the data
 * attached to the hash table is a position,size pair. This allows for the
 * hash table to encode positions and sizes of items within a file archive.
 * --------------------------------------------------------------------
 */

/*
 * Writes the HashTable structures to 'fp'.
 * This is a specialisation of the HashTable where the HashData is a
 * position,size tuple.
 *
 * This consists of the following format:
 * Header:
 *    ".hsh" (magic numebr)
 *    x3     (3-bytes of version code, eg "000")
 *    x1     (HASH_FUNC_? function used)
 *    x4     (4-bytes big-endian; number of bytes in stored hash, inc. header)
 *    x4     (4-bytes big-endian; number of hash buckets)
 * Archive name:
 *    x1     (length, zero => no name)
 *    ?      (archive filename)
 * Buckets (multiples of)
 *    x4     (4-byte offset of linked list pos,  rel. to the start of the hdr)
 * Items (per bucket chain, not written if Bucket[?]==0)
 *    x1     (key length, zero => end of chain)
 *    ?      (key)
 *    x8     (position)
 *    x4     (size)
 *
 * It is designed such that on-disk querying of the hash table can be done
 * purely by forward seeks. (This is generally faster due to pre-fetching of
 * the subsequent blocks by many disk controllers.)
 */
void HashFileSave(HashTable *h, FILE *fp, char *archive) {
    int i;
    HashItem *hi;
    HashFileHeader hh;
    uint64_t *bucket_pos;
    uint64_t offset = 0;

   
    /* Compute the coordinates of the hash items */
    offset = sizeof(hh) + 1 + (archive?strlen(archive):0) + h->nbuckets * 4;
    bucket_pos = (uint64_t *)calloc(h->nbuckets, sizeof(uint64_t));
    for (i = 0; i < h->nbuckets; i++) {
	bucket_pos[i] = offset;

	if (!(hi = h->bucket[i]))
	    continue;
	for (; hi; hi = hi->next) {
	    offset += 1 + hi->key_len + 8 + 4;
	}
	offset++;
    }

    /* Write the header: */
    hh.magic[0] = '.';
    hh.magic[1] = 'h';
    hh.magic[2] = 's';
    hh.magic[3] = 'h';
    hh.vers[0]  = '0';
    hh.vers[1]  = '0';
    hh.vers[2]  = '0';
    hh.hfunc    = h->options & HASH_FUNC_MASK;
    hh.nbuckets = be_int4(h->nbuckets);
    hh.size     = be_int4(offset);
    fwrite(&hh, sizeof(hh), 1, fp);

    /* Write the archive filename, if known */
    if (archive) {
	fputc(strlen(archive), fp);
	fputs(archive, fp);
    } else {
	fputc(0, fp);
    }

    /* Write out hash buckets */
    for (i = 0; i < h->nbuckets; i++) {
	uint32_t zero = 0;
	uint32_t be32;

	if (!(hi = h->bucket[i])) {
	    fwrite(&zero, 4, 1, fp);
	    continue;
	}

	be32 = be_int4(bucket_pos[i]);
	fwrite(&be32, 4, 1, fp);
    }
    free(bucket_pos);

    /*
     * Finally write the hash_item linked lists. The first item is the
     * hash key length. We append a zero to the end of the list so we
     * can check this key length to determine the end of this hash
     * item list.
     */
    for (i = 0; i < h->nbuckets; i++) {
	if (!(hi = h->bucket[i]))
	    continue;
	for (; hi; hi = hi->next) {
	    uint64_t be64;
	    uint32_t be32;
	    fprintf(fp, "%c%.*s", hi->key_len,
		    hi->key_len, hi->key);
	    be64 = be_int8(hi->data.idx.pos);
	    fwrite(&be64, 8, 1, fp);
	    be32 = be_int4(hi->data.idx.len);
	    fwrite(&be32, 4, 1, fp);
	}
        fputc(0, fp);
    }
}

/*
 * Reads an entire HashTable from fp.
 *
 * Returns:
 *    A filled out HashTable pointer on success
 *    NULL on failure   
 */
HashTable *HashFileLoad(FILE *fp) {
    int i;
    HashTable *h;
    HashItem *hi;
    HashFileHeader hh;
    uint64_t *bucket_pos;
    unsigned char *htable;
    int htable_pos;
    char fname[256];
    int fnamelen;

    if (NULL == (htable = (unsigned char *)malloc(20)))
	return NULL;
    if (sizeof(hh) != fread(htable, 1, sizeof(hh), fp))
	return NULL;

    /* Read and create the hash table header */
    memcpy(&hh, htable, sizeof(hh));
    hh.nbuckets = be_int4(hh.nbuckets);
    hh.size = be_int4(hh.size);
    h = HashTableCreate(hh.nbuckets, hh.hfunc);
    bucket_pos = (uint64_t *)calloc(h->nbuckets, sizeof(uint64_t));

    /* Load the archive filename */
    fnamelen = fgetc(fp);
    if (fnamelen) {
	fread(fname, 1, fnamelen, fp);
    }
    fname[fnamelen] = 0;

    /* Load the rest of the hast table to memory */
    htable_pos = sizeof(hh) + fnamelen + 1;
    if (NULL == (htable = (unsigned char *)realloc(htable, hh.size)))
	return NULL;
    if (hh.size != fread(&htable[htable_pos], 1, hh.size, fp))
	return NULL;

    /* Identify the "bucket pos". Detemines which buckets have data */
    for (i = 0; i < h->nbuckets; i++) {
	memcpy(&bucket_pos[i], &htable[htable_pos], 4);
	bucket_pos[i] = be_int4(bucket_pos[i]);
	htable_pos += 4;
    }

    /* Read the hash table items */
    for (i = 0; i < h->nbuckets; i++) {
	if (!bucket_pos[i])
	    continue;
	for (;;) {
	    int c;
	    char key[256];
	    uint64_t pos;
	    uint32_t size;

	    c = htable[htable_pos++];
	    if (c == EOF || !c)
		break;

	    memcpy(key, &htable[htable_pos], c);
	    htable_pos += c;
	    memcpy(&pos, &htable[htable_pos], 8);
	    htable_pos += 8;
	    pos = be_int8(pos);
	    memcpy(&size, &htable[htable_pos], 4);
	    htable_pos += 4;
	    size = be_int4(size);

	    hi = HashItemCreate(h);
	    hi->next = h->bucket[i];
	    h->bucket[i] = hi;
	    hi->key_len = c;
	    hi->key = (char *)malloc(c+1);
	    memcpy(hi->key, key, c);
	    hi->key[c] = 0; /* For debugging convenience only */
	    hi->data.idx.pos = pos;
	    hi->data.idx.len = size;
	}
    }

    fprintf(stderr, "done\n");
    fflush(stderr);
    free(bucket_pos);

    return h;
}

/*
 * Opens a stored hash table file. It also internally keeps an open file to
 * hash and the archive files.
 *
 * Returns the HashFile pointer on success
 *         NULL on failure
 */
HashFile *HashFileOpen(char *fname) {
    HashFile *hf = (HashFile *)malloc(sizeof(*hf));
    char archive[256];
    int archive_len;

    hf->hfp = hf->afp = NULL;

    /* Open the hash and read the header */
    if (NULL == (hf->hfp = fopen(fname, "rb")))
	return NULL;

    if (sizeof(hf->h) != fread(&hf->h, 1, sizeof(hf->h), hf->hfp)) {
	HashFileClose(hf);
	return NULL;
    }
    hf->h.nbuckets = be_int4(hf->h.nbuckets);

    /* Load the main archive filename */
    if ((archive_len = fgetc(hf->hfp)))
	fread(archive, 1, archive_len, hf->hfp);
    archive[archive_len] = 0;

    hf->header_size = sizeof(hf->h) + 1 + archive_len;

    /* Open the main archive too */
    if (NULL == (hf->afp = fopen(archive, "rb"))) {
	HashFileClose(hf);
	return NULL;
    }

    return hf;
}

/*
 * Closes a and deallocates a HashFile opened via HashFileOpen()
 */
void HashFileClose(HashFile *hf) {
    if (!hf)
	return;

    if (hf->hfp)
	fclose(hf->hfp);
    if (hf->afp)
	fclose(hf->afp);
    free(hf);
}

/*
 * Searches the named HashFile for a specific key.
 * When found it returns the position and size of the object in pos and size.
 *
 * Returns
 *    0 on success (pos & size updated)
 *   -1 on failure
 */
int HashFileQuery(HashFile *hf, uint8_t *key, int key_len,
		  uint64_t *r_pos, uint32_t *r_size) {
    uint32_t hval;
    uint32_t pos;
    int klen;
    int cur_offset = 0;

    /* Hash 'key' to compute the bucket number */
    hval = hash(hf->h.hfunc, key, key_len) & (hf->h.nbuckets-1);

    /* Read the bucket to find the first linked list item location */
    if (-1 == fseek(hf->hfp, 4*hval + hf->header_size, SEEK_SET))
	return -1;
    if (4 != fread(&pos, 1, 4, hf->hfp))
	return -1;
    pos = be_int4(pos);
    cur_offset = 4*hval + 4 + hf->header_size;

    /* Jump to the HashItems list and look through for key */
    if (-1 == fseek(hf->hfp, pos - cur_offset, SEEK_CUR))
	return -1;

    for (klen = fgetc(hf->hfp); klen; klen = fgetc(hf->hfp)) {
	char k[256];
	uint64_t pos;
	uint32_t size;

	fread(k, klen, 1, hf->hfp);
	fread(&pos, 8, 1, hf->hfp);
	pos = be_int8(pos);
	fread(&size, 4, 1, hf->hfp);
	size = be_int4(size);
	if (klen == key_len && 0 == memcmp(key, k, key_len)) {
	    *r_pos = pos;
	    *r_size = size;
	    return 0;
	}
    }

    return -1;
}
