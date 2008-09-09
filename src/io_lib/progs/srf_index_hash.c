/*
 * This adds a hash table index (".hsh" v1.01 format) to an SRF archive.
 * It does this either inline on the file itself (provided it doesn't already
 * have an index) or by producing an external index file.
 */

/* ---------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <io_lib/hash_table.h>
#include <io_lib/os.h>
#include <io_lib/array.h>
#include <io_lib/srf.h>

typedef struct {
    uint64_t pos;
    uint32_t dbh;
} pos_dbh;

/*
 * Writes the HashTable structures to 'fp'.
 * This is a specialisation of the HashTable where the HashData is a
 * position,size tuple.
 *
 * Header:
 *   x4    magic number, starting with 'I'.
 *   x4    version code (eg "1.00")
 *   x8    index size (should be x4 as we assume bucket locs are x4?)
 *
 *   x1    index type ('E' normally)
 *   x1    dbh_pos_stored_sep (indicates if the item list contains the
 *         "data block header" index number).
 *
 *   x4    number of containers
 *   x4    number of DBHs
 *   x8    number of hash buckets
 *
 *   x*    dbhFile  p-string (NULL if held within the same file)
 *   x*    contFile p-string (NULL if held within the same file)
 *
 * Containers: (1 entry per container)
 *   x8    file position of container header
 *
 * Data Block Headers: (1 entry per DBH)
 *   x8    file position of container header
 *
 * Buckets: (1 entry per bucket)
 *   x8    8-byte offset of linked list pos,  rel. to the start of the hdr
 *
 * Items: (1 per trace)
 *   x1    name disambiguation hash, top-most bit set => last item in list
 *   x8    data position
 *  (x4)  (dbh_index - optional; present if dbh_pos_stored_sep is 1)
 *
 * Footer:
 *   x4    magic number
 *   x4    version
 *   x8    index size
 *
 * Returns: the number of bytes written on success
 *         -1 for error
 */
int HFSave(char *ch_file, Array ch_pos,
	   char *th_file, Array th_pos,
	   int dbh_pos_stored_sep,
	   HashTable *h, srf_t *srf) {
    unsigned int i, j;
    srf_index_hdr_t hdr;
    uint64_t *bucket_pos;
    int item_sz;

    /* Option: whether to store dbh positions directly in the index */
    hdr.dbh_pos_stored_sep = dbh_pos_stored_sep;

    /* Compute index size and bucket offsets */
    hdr.size = 34 +
	1 + (ch_file ? strlen(ch_file) : 0) +
	1 + (th_file ? strlen(th_file) : 0);
    hdr.size += 8*(ArrayMax(ch_pos) + ArrayMax(th_pos) + h->nbuckets);
    if (NULL == (bucket_pos = (uint64_t *)calloc(h->nbuckets,
						 sizeof(*bucket_pos))))
	return -1;
    for (i = 0; i < h->nbuckets; i++) {
	HashItem *hi;
	if (!(hi = h->bucket[i]))
	    continue;
	bucket_pos[i] = hdr.size;
	item_sz = 1 + 8 + (hdr.dbh_pos_stored_sep ? 4 : 0);
	for (; hi; hi = hi->next)
	    hdr.size += item_sz;
    }
    hdr.size += 16; /* footer */

    /* Construct and write out the index header */
    memcpy(hdr.magic,   SRF_INDEX_MAGIC,   4);
    memcpy(hdr.version, SRF_INDEX_VERSION, 4);
    hdr.index_type = 'E';
    hdr.n_container = ArrayMax(ch_pos);
    hdr.n_data_block_hdr = ArrayMax(th_pos);
    hdr.n_buckets = h->nbuckets;
    if (th_file)
	strncpy(hdr.dbh_file,  th_file, 255);
    else
	hdr.dbh_file[0] = 0;
    if (ch_file)
	strncpy(hdr.cont_file, ch_file, 255);
    else
	hdr.cont_file[0] = 0;
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
	if (0 != srf_write_uint64(srf, bucket_pos[i]))
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
	    uint64_t pos;
	    uint32_t dbh;
	    uint32_t h7;

	    if (hdr.dbh_pos_stored_sep) {
		pos_dbh *pdbh = (pos_dbh *)hi->data.p;
		pos = pdbh->pos;
		dbh = pdbh->dbh;
	    } else {
		pos = hi->data.i;
	    }

	    /* Rehash key in 7 bits;  */
	    h7 = hash64(h->options & HASH_FUNC_MASK,
			(uint8_t *)hi->key, hi->key_len) >> 57;
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

	    if (hdr.dbh_pos_stored_sep)
		if (0 != srf_write_uint32(srf, dbh))
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

/* ------------------------------------------------------------------------ */
void usage(int code) {
    printf("Usage: srf_index_hash [-c] srf_file\n");
    printf(" Options:\n");
    printf("    -c       check an existing index, don't re-index\n");
    exit(code);
}

int main(int argc, char **argv) {
    srf_t *srf;
    uint64_t pos;
    char name[512];
    int i, type, new;
    char *archive;
    HashTable *db_hash;
    Array ch_pos, th_pos;
    int dbh_pos_stored_sep = 0;
    pos_dbh *pdbh;
    int check = 0;
    off_t old_index = 0;
    
    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;
	} else if (!strcmp(argv[i], "-c")) {
	    check = 1;
	} else if (!strcmp(argv[i], "-h")) {
	    usage(0);
	} else {
	    usage(1);
	}
    }

    if (argc != (i+1)) {
      usage(1);
    }

    archive = argv[i];

    if( check ){
        srf = srf_open(archive, "rb");
    } else {
        srf = srf_open(archive, "r+b");
    }
    if (NULL == srf ){
 	perror(argv[i]);
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

	/* Only want this set if the last block in the file is an index */
	old_index = 0;

	switch (type) {
	case SRFB_CONTAINER:
	    ARR(uint64_t, ch_pos, ArrayMax(ch_pos)) = pos;
	    break;

	case SRFB_TRACE_HEADER:
	    ARR(uint64_t, th_pos, ArrayMax(th_pos)) = pos;
	    break;

	case SRFB_TRACE_BODY:
	    if (dbh_pos_stored_sep) {
		if (NULL == (pdbh = (pos_dbh *)malloc(sizeof(*pdbh))))
		    return 1;
		pdbh->pos = pos;
		pdbh->dbh = ArrayMax(th_pos);
		hd.p = pdbh;
	    } else {
		hd.i = pos;
	    }
	    HashTableAdd(db_hash, name, strlen(name), hd, &new);
            if (0 == new) {
              fprintf(stderr, "duplicate read name %s\n", name);
              return 1;
            }
	    break;

	case SRFB_INDEX:
	    /* An old index */
	    old_index = pos;
	    break;

        case SRFB_NULL_INDEX:
            break;

	default:
	    abort();
	}
    }

    /* the type should be -1 (EOF) */
    if( type != -1 )
        abort();

    /* are we really at the end of the srf file */
    pos = ftell(srf->fp);
    fseek(srf->fp, 0, SEEK_END);
    if( pos != ftell(srf->fp) ){
        fprintf(stderr, "srf file is corrupt\n");
	return 1;
    }
    
    if( check ){
      ArrayDestroy(th_pos);
      srf_destroy(srf, 1);
      return 0;
    }

    /* Write out the index */
    if (old_index)
	fseeko(srf->fp, old_index, SEEK_SET);

    HashTableStats(db_hash, stderr);
    HFSave(NULL, ch_pos,
	   NULL, th_pos,
	   dbh_pos_stored_sep,
	   db_hash, srf);

    /* Truncate incase we've somehow overwritten an old longer index */
    ftruncate(fileno(srf->fp), ftello(srf->fp));

    HashTableDestroy(db_hash, 0);
    ArrayDestroy(ch_pos);
    ArrayDestroy(th_pos);
    srf_destroy(srf, 1);
    
    return 0;
}
