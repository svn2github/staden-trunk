/*
 * This adds a hash table index (".hsh" v1.00 format) to an SRF archive.
 * It does this either inline on the file itself (provided it doesn't already
 * have an index) or by producing an external index file.
 */

/* ---------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "hash_table.h"
#include "os.h"
#include "array.h"
#include "srf.h"


/*
 * Writes the HashTable structures to 'fp'.
 * This is a specialisation of the HashTable where the HashData is a
 * position,size tuple.
 *
 * Header:
 *   x4    magic number, starting with 'I'.
 *   x4    version code (eg "1.00")
 *   x8    index size (should be x4 as we assume bucket locs are x4?)
 *   x4    number of containers
 *   x4    number of DBHs
 *   x4    number of hash buckets (~10 billion traces per file is enough).
 *   x1    hash function (should be a 64-bit one)
 *
 * Containers: (1 entry per container)
 *   x8    file position of container header
 *
 * Data Block Headers: (1 entry per DBH)
 *   x8    file position of container header
 *
 * Buckets: (1 entry per bucket)
 *   x4    4-byte offset of linked list pos,  rel. to the start of the hdr
 *
 * Items: (1 per trace)
 *   x1    name disambiguation hash, top-most bit set => last item in list
 *   x8    data position
 *
 * Footer:
 *   x4    magic number
 *   x4    version
 *   x8    index size
 *
 * It is designed such that on-disk querying of the hash table can be done
 * purely by forward seeks. (This is generally faster due to pre-fetching of
 * the subsequent blocks by many disk controllers.)
 *
 * Returns: the number of bytes written on success
 *         -1 for error
 */
int HFSave(Array ch_pos, Array th_pos, HashTable *h, srf_t *srf) {
    unsigned int i, j;
    srf_index_hdr_t hdr;
    uint64_t *bucket_pos;

    /* Compute index size and bucket offsets */
    hdr.size = SRF_INDEX_HDR_SIZE;
    hdr.size += ArrayMax(ch_pos) * 8;
    hdr.size += ArrayMax(th_pos) * 8;
    hdr.size += h->nbuckets * 4;
    hdr.hash_func = h->options & HASH_FUNC_MASK;
    if (NULL == (bucket_pos = (uint64_t *)calloc(h->nbuckets,
						 sizeof(*bucket_pos))))
	return -1;
    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	if (!(hi = h->bucket[i]))
	    continue;
	bucket_pos[i] = hdr.size;
	for (; hi; hi = hi->next)
	    hdr.size += 9;
    }
    hdr.size += 16; /* footer */

    /* Construct and write out the index header */
    memcpy(hdr.magic,   SRF_INDEX_MAGIC,   4);
    memcpy(hdr.version, SRF_INDEX_VERSION, 4);
    hdr.n_container = ArrayMax(ch_pos);
    hdr.n_data_block_hdr = ArrayMax(th_pos);
    hdr.n_buckets = h->nbuckets;
    if (0 != srf_write_index_hdr(srf, &hdr))
	return -1;

    /* Write the container and data block header arrays */
    j = ArrayMax(ch_pos);
    for (i = 0; i < j; i++) {
	if (0 != srf_write_uint64(srf, arr(uint64_t, ch_pos, i)))
	    return -1;
    }

    j = ArrayMax(th_pos);
    for (i = 0; i < j; i++) {
	if (0 != srf_write_uint64(srf, arr(uint64_t, th_pos, i)))
	    return -1;
    }

    /* Write out buckets */
    for (i = 0; i < h->nbuckets; i++) {
	if (0 != srf_write_uint32(srf, bucket_pos[i]))
	    return -1;
    }

    /* Write out the trace locations themselves */
    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	/*
	fprintf(stderr, "Bucket %d offset %lld vs %lld\n",
		i, ftell(srf->fp), bucket_pos[i]);
	*/
	if (!(hi = h->bucket[i]))
	    continue;
	for (; hi; hi = hi->next) {
	    uint64_t pos = hi->data.i;
	    uint32_t h7;

	    /* Rehash key in 7 bits;  */
	    h7 = hash64(hdr.hash_func, (uint8_t *)hi->key, hi->key_len) >> 57;
	    if (!hi->next)
		h7 |= 0x80;
	    /*
	    fprintf(stderr, "\t%.*s => %x @ %lld\n",
		    hi->key_len, hi->key, h7, pos);
	    */
	    if (fputc(h7, srf->fp) < 0)
		return -1;
	    if (0 != srf_write_uint64(srf, pos))
		return -1;
	}
    }

    /* Footer */
    if (4 != fwrite(hdr.magic,   1, 4, srf->fp))
	return -1;
    if (4 != fwrite(hdr.version, 1, 4, srf->fp))
	return -1;
    if (0 != srf_write_uint64(srf, hdr.size))
	return -1;

    return 0;
}

int main(int argc, char **argv) {
    srf_t *srf;
    uint64_t pos;
    char name[512];
    int type;
    char *archive;
    HashTable *db_hash;
    Array ch_pos, th_pos;
    
    if (argc != 2) {
	fprintf(stderr, "Usage: hash_srf srf_file\n");
	return 1;
    }
    archive = argv[1];

    if (NULL == (srf = srf_open(archive, "r+b"))) {
	perror(argv[1]);
	return 1;
    }

    /* Create the arrays and hash table */
    if (!(ch_pos = ArrayCreate(sizeof(uint64_t), 0)))
	return 1;

    if (!(th_pos = ArrayCreate(sizeof(uint64_t), 0)))
	return 1;

    if (!(db_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS3)))
	return 1;

    /* Scan through file gathering the details to index in memory */
    while ((type = srf_next_block_details(srf, &pos, name)) >= 0) {
	HashData hd;

	switch (type) {
	case SRFB_CONTAINER:
	    ARR(uint64_t, ch_pos, ArrayMax(ch_pos)) = pos;
	    break;

	case SRFB_TRACE_HEADER:
	    ARR(uint64_t, th_pos, ArrayMax(th_pos)) = pos;
	    break;

	case SRFB_TRACE_BODY:
	    hd.i = pos;
	    HashTableAdd(db_hash, name, strlen(name), hd, NULL);
	    break;

	default:
	    abort();
	}
    }

    /* Write out the index */
    HashTableStats(db_hash, stderr);
    HFSave(ch_pos, th_pos, db_hash, srf);

    HashTableDestroy(db_hash, 0);
    ArrayDestroy(ch_pos);
    ArrayDestroy(th_pos);
    srf_destroy(srf, 1);
    
    return 0;
}
