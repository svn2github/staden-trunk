#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>

#include "tg_gio.h"
#include "misc.h"

static int load_counts[100];
static int unload_counts[100];
static int write_counts[100];

/*
 * This module implements a layer of caching on top of the underlying
 * IO interface.
 *
 * The cache is basically a hash table with the hash key being type/rec_num
 * and the 'payload' being the object itself. For efficiencies sake this
 * isn't the on-disk representation but rather a decoded data structure
 * holding a more C-friendly representation. Hence this code also contains
 * the encoding/decoding functions.
 *
 * Objects are kept in cache until either there is no more room or until
 * they are written to disk.
 *
 * When reading we create a RO view (unless explicitly asked for otherwise).
 *
 * When writing we upgrade the view to RW and increment it's reference count
 * to prevent cache expiry.
 */

/*
 * FIXME: the default refcount should be zero, so if you just loop through
 * reading data it's automatically deallocated.
 *
 * If you need to keep it around for a while then the code should
 * explicitly increment the reference count.
 */

/* The cache key - a combination of record number and data type. */
typedef struct {
    GRec rec;
    char type;
} cache_key_t;

/*
 * Allocates and initialises a new cached_item object from a type, view and
 * HacheItem.
 * 'e_len' is the amount of extra storage needed to house the gap object
 * itself.
 *
 * Returns ptr on success.
 *         NULL on failure
 */
cached_item *cache_new(int type, GRec rec, GView v,
		       HacheItem *hi, size_t e_len) {
    cached_item *ci;

    if (NULL == (ci = (cached_item *)malloc(sizeof(cached_item) + e_len)))
	return NULL;

    ci->view = v;
    ci->rec = rec;
    ci->lock_mode = G_LOCK_RO;
    ci->type = type;
    ci->hi = hi;
    ci->updated = 0;
    ci->data_size = e_len;

    return ci;
}

/* ----------------------------------------------------------------------
 * Callback functions used by the Hache table.
 */


/*
 * A forced unload of a sequence.
 * This is called by the cache when it has insufficient space to load new
 * data into. It will only be called on objects with a zero reference count,
 * which in the context of this code means they have not been locked R/W.
 */
static void seq_unload(GapIO *io, cached_item *ci) {
    seq_t *s = (seq_t *)&ci->data;
    io->iface->seq.unlock(io->dbh, ci->view);

    if (s->anno)
	ArrayDestroy(s->anno);

    free(ci);
}

/*
 * Writes a sequence to disc.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
static int seq_write(GapIO *io, cached_item *ci) {
    return io->iface->seq.write(io->dbh, ci);
}


/*
 * A forced unload of a sequence block.
 *
 * This is called by the cache when it has insufficient space to load new
 * data into. It will only be called on objects with a zero reference count,
 * which in the context of this code means they have not been locked R/W.
 */
static void seq_block_unload(GapIO *io, cached_item *ci) {
    int i;
    seq_block_t *b = (seq_block_t *)&ci->data;

    io->iface->seq_block.unlock(io->dbh, ci->view);

    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
	seq_t *s = b->seq[i];
	cached_item *si;

	if (s) {
	    if (s->anno)
		ArrayDestroy(s->anno);
	    si = ci_ptr(s);
	    free(si);
	}
    }
    free(ci);
}

/*
 * Writes a sequence block to disc.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
static int seq_block_write(GapIO *io, cached_item *ci) {
    return io->iface->seq_block.write(io->dbh, ci);
}


/*
 * A forced unload of a track.
 */
static void track_unload(GapIO *io, cached_item *ci) {
    track_t *track = (track_t *)&ci->data;

    io->iface->track.unlock(io->dbh, ci->view);

    if (track->data)
	ArrayDestroy(track->data);
    free(ci);
}

/*
 * Writes a track to disc.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
static int track_write(GapIO *io, cached_item *ci) {
    return io->iface->track.write(io->dbh, ci);
}


/*
 * A forced unload of a bin.
 */
static void bin_unload(GapIO *io, cached_item *ci) {
    bin_index_t *bin = (bin_index_t *)&ci->data;

    io->iface->bin.unlock(io->dbh, ci->view);

    if (bin->rng)
	ArrayDestroy(bin->rng);
    if (bin->track)
	ArrayDestroy(bin->track);
    free(ci);
}

/*
 * Writes a bin to disc.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
static int bin_write(GapIO *io, cached_item *ci) {
    return io->iface->bin.write(io->dbh, ci);
}


static void contig_unload(GapIO *io, cached_item *ci) {
    io->iface->seq.unlock(io->dbh, ci->view);
    free(ci);
}

/*
 * Writes a contig structure to disc.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
static int contig_write(GapIO *io, cached_item *ci) {
    return io->iface->contig.write(io->dbh, ci);
}

static int array_write(GapIO *io, cached_item *ci) {
    return io->iface->array.write(io->dbh, ci);
}

static void array_unload(GapIO *io, cached_item *ci) {
    Array ar = (Array)&ci->data;
    //ArrayDestroy(ar);
    free(ar->base);

    io->iface->seq.unlock(io->dbh, ci->view);
    free(ci);
}

static void database_unload(GapIO *io, cached_item *ci) {
    io->iface->database.unlock(io->dbh, ci->view);
    free(ci);
}

/*
 * Writes an annotation / annotation_element to disc.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
static int anno_ele_write(GapIO *io, cached_item *ci) {
    return io->iface->anno_ele.write(io->dbh, ci);
}
static int anno_write(GapIO *io, cached_item *ci) {
    return io->iface->anno.write(io->dbh, ci);
}

static void anno_ele_unload(GapIO *io, cached_item *ci) {
    io->iface->anno_ele.unlock(io->dbh, ci->view);
    free(ci);
}
static void anno_unload(GapIO *io, cached_item *ci) {
    io->iface->anno.unlock(io->dbh, ci->view);
    free(ci);
}


/*
 * Writes an library to disk.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
static int library_write(GapIO *io, cached_item *ci) {
    return io->iface->library.write(io->dbh, ci);
}

static void library_unload(GapIO *io, cached_item *ci) {
    io->iface->library.unlock(io->dbh, ci->view);
    free(ci);
}


/*
 * Resize a cached_item object to a new size.
 * The object passed in and returned is the object itself in the cache, but
 * the memory we realloc is the cached_item struct holding that.
 *
 * This may need to reset the HacheTable pointer back to this object too
 * has it holds (key,value) pairings where value is the pointer to this
 * item.
 *
 * Returns new void* item pointer on success
 *        NULL on failure
 */
void *cache_item_resize(void *item, size_t size) {
    cached_item *ci = ci_ptr(item);
    cached_item *new = (cached_item *)realloc(ci, size + sizeof(*ci));

    if (NULL == new)
	return NULL;

    if (ci == new)
	return item;

    if (new->hi) {
	assert(new->hi->data.p == ci);
	new->data_size = size;
	new->hi->data.p = new;
    }

    if (new->type == GT_Seq) {
	seq_t *s = (seq_t *)&new->data;
	s->block->seq[s->idx] = s;
    }

    return &new->data;
}

/*
 * Called when attempting to load a new record.
 * 'key' here is a 4-byte record number followed by a 1 byte type field -
 * ie a cache_key_t struct.
 */
static HacheData *cache_load(void *clientdata, char *key, int key_len,
			    HacheItem *hi) {
    GapIO *io = (GapIO *)clientdata;
    cached_item *ci;

    cache_key_t *k = (cache_key_t *)key;
    static HacheData hd;

    //printf("Cache load %d type %d\n", k->rec, k->type);

    load_counts[k->type]++;

    switch (k->type) {
    case GT_Database:
	ci = io->iface->database.read(io->dbh, k->rec);
	break;
	
    case GT_Seq:
	ci = io->iface->seq.read(io->dbh, k->rec);
	printf("Query seq %d => %p\n", k->rec, ci);
	break;

    case GT_SeqBlock:
	/* FIXME
	 * At this point we need to reparent the unpacked sequences
	 * and populate as objects in the cache too.
	 */
	ci = io->iface->seq_block.read(io->dbh, k->rec);
	break;

    case GT_Bin:
	ci = io->iface->bin.read(io->dbh, k->rec);
	break;

    case GT_Track:
	ci = io->iface->track.read(io->dbh, k->rec);
	break;

    case GT_Contig:
	ci = io->iface->contig.read(io->dbh, k->rec);
	break;

    case GT_RecArray:
	ci = io->iface->array.read(io->dbh, k->rec);
	break;

    case GT_AnnoEle:
	ci = io->iface->anno_ele.read(io->dbh, k->rec);
	break;

    case GT_Anno:
	ci = io->iface->anno.read(io->dbh, k->rec);
	break;

    case GT_Library:
	ci = io->iface->library.read(io->dbh, k->rec);
	break;

    default:
	return NULL;
    }

    if (!ci)
	return NULL;

    hd.p = ci;
    ci->hi = hi;

    /*
     * After thinking long and hard I've decided that all items we load
     * should start off with a zero reference count. This means tight loops
     * are fine (eg bin=get_bin(bin->parent)) and also forces the number
     * of cache_incr to be the same as cache_decr which is good for bug
     * spotting.
     *
     * Finally it means we do not get issues when saving data as to whether
     * we need to decrement the reference count.
     */
    HacheTableDecRef(io->cache, hi);

    return &hd;
}

/* Callback from Hache */
static void cache_unload(void *clientdata, HacheData hd) {
    GapIO *io = (GapIO *)clientdata;
    cached_item *ci = hd.p;

    //    printf("Cache unload %d\n", ci->rec);

    assert(ci->updated == 0);

    unload_counts[ci->type]++;

    switch (ci->type) {
    case GT_Seq:
	seq_unload(io, ci);
	break;

    case GT_SeqBlock:
	/* FIXME:
	 * At this stage we need to decrement the reference counts on
	 * the child GT_Seq objects and unload those too.
	 */
	seq_block_unload(io, ci);
	break;

    case GT_Bin:
	bin_unload(io, ci);
	break;

    case GT_Track:
	track_unload(io, ci);
	break;

    case GT_Contig:
	contig_unload(io, ci);
	break;

    case GT_RecArray:
	array_unload(io, ci);
	break;

    case GT_Database:
	database_unload(io, ci);
	break;

    case GT_AnnoEle:
	anno_ele_unload(io, ci);
	break;

    case GT_Anno:
	anno_unload(io, ci);
	break;

    case GT_Library:
	library_unload(io, ci);
	break;
    }
}



/* ----------------------------------------------------------------------
 * External interfaces - create/destroy/query/update etc
 */

/*
 * Creates the IO Hache
 *
 * Returns 0 on success
 *         -1 failure.
 */
int cache_create(GapIO *io) {
    HacheTable *h;

    //    if (NULL == (h = HacheTableCreate(131072, HASH_DYNAMIC_SIZE|HASH_OWN_KEYS)))
    if (NULL == (h = HacheTableCreate(8192, HASH_DYNAMIC_SIZE|HASH_OWN_KEYS)))
	return -1;

    h->clientdata = io;
    h->load = cache_load;
    h->del  = cache_unload;

    io->cache = h;
    return 0;
}

void cache_destroy(GapIO *io) {
    int i;
    HacheTable *h = io->cache;

    if (!h)
	return;

    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    cache_unload(io, hi->data);
	}
    }

    HacheTableDestroy(io->cache, 0);
}

/*
 * Ensure key is initialised correctly, making sure that it is blank
 * in the gaps between structure elements and the padding at the end.
 */
cache_key_t construct_key(GRec rec, int type) {
    cache_key_t k;
    memset(&k, 0, sizeof(k));
    k.rec = rec;
    k.type = type;
    return k;
}


/*
 * Upgrades a lock on a view to a higher level, eg from read-only to
 * read-write or exclusive access.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
int cache_upgrade(GapIO *io, cached_item *ci, int mode) {
    int ret;

    switch(ci->type) {
    case GT_Database:
	ret = io->iface->database.upgrade(io->dbh, ci->view, mode);
	ci->lock_mode = mode;
	break;
	
    case GT_Seq:
	ret = io->iface->seq.upgrade(io->dbh, ci->view, mode);
	ci->lock_mode = mode;
	break;

    case GT_SeqBlock:
	ret = io->iface->seq_block.upgrade(io->dbh, ci->view, mode);
	ci->lock_mode = mode;
	break;

    case GT_Contig:
	ret = io->iface->contig.upgrade(io->dbh, ci->view, mode);
	ci->lock_mode = mode;
	break;
	
    case GT_Bin:
	ret = io->iface->bin.upgrade(io->dbh, ci->view, mode);
	ci->lock_mode = mode;
	break;
	
    case GT_Track:
	ret = io->iface->track.upgrade(io->dbh, ci->view, mode);
	ci->lock_mode = mode;
	break;
	
    case GT_RecArray:
	ret = io->iface->array.upgrade(io->dbh, ci->view, mode);
	ci->lock_mode = mode;
	break;
 
    case GT_AnnoEle:
	ret = io->iface->anno_ele.upgrade(io->dbh, ci->view, mode);
	ci->lock_mode = mode;
	break;

    case GT_Anno:
	ret = io->iface->anno.upgrade(io->dbh, ci->view, mode);
	ci->lock_mode = mode;
	break;

    case GT_Library:
	ret = io->iface->library.upgrade(io->dbh, ci->view, mode);
	ci->lock_mode = mode;
	break;

    default:
	return -1;
    }

    /*
     * Bump reference count to indicate we cannot expire this record.
     */
    //    if (ret == 0)
    //	HacheTableIncRef(io->cache, ci->hi);

    return ret;
}

/*
 * Upgrades a lock on a view to a higher level, eg from read-only to
 * read-write or exclusive access.
 *
 * For overlayed GapIOs we make sure the locked item is in this layer.
 *
 * Returns locked object pointer on success
 *        NULL for failure.
 */
void *cache_lock(GapIO *io, int type, GRec rec, int mode) {
    cache_key_t key = construct_key(rec, type);
    cached_item *ci;
    HacheTable *h = io->cache;
    HacheItem *hi = HacheTableSearch(h, (char *)&key, sizeof(key));
    
    if (!hi)
	return NULL;

    if (NULL == (ci = hi->data.p))
	return NULL;

    if (cache_upgrade(io, ci, mode) == 0)
	return &ci->data;
    else
	return NULL;
}

/*
 * qsort callback.
 * Sorts cached_item * on view number.
 */
int qsort_ci_view(const void *p1, const void *p2) {
    cached_item **c1 = (cached_item **)p1;
    cached_item **c2 = (cached_item **)p2;

    /* FIXME: use c2-c1 if lock_file_N is used as it reverses the order */
    return (*c2)->view - (*c1)->view;
}

/*
 * Flushes changes in the cache back to disk.
 */
int cache_flush(GapIO *io) {
    int i, ret = 0;
    HacheTable *h = io->cache;
    HacheItem *hi;
    Array to_flush;
    int nflush = 0;

    //printf(">>> cache flush <<<\n");
    //HacheTableRefInfo(io->cache, stdout);

    /*
     * If this is a derived io then we need to pass these items up to
     * our parent, remove them from this io and call cache_flush on the
     * parent instead (possibly recursively).
     */
    if (io->base) {
	for (i = 0; i < h->nbuckets; i++) {
	    HacheItem *hi, *next;
	    for (hi = h->bucket[i]; hi; hi = next) {
		HacheData data;
		cached_item *ci = hi->data.p;

		next = hi->next;

		if (!ci->updated)
		    continue;

		/*
		 * For blocked data structures, merge with base copy first.
		 */
		switch (ci->type) {
		case GT_SeqBlock: {
		    HacheItem *htmp;
		    seq_block_t *bn = (seq_block_t *)&ci->data;
		    seq_block_t *bo;
		    int i;

		    htmp = HacheTableSearch(io->base->cache,
					    hi->key, hi->key_len);
		    bo = (seq_block_t *)&((cached_item *)htmp->data.p)->data;

		    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
			if (!bn->seq[i]) {
			    bn->seq[i] = bo->seq[i];
			    bn->seq[i]->block = bn;
			}
		    }

		    break;
		}
		}

		/* Purge from parent */
		HacheTableRemove(io->base->cache,
				 ci->hi->key, ci->hi->key_len, 0);

		/* FIXME: search for parent ci first and deallocate after */

		/* Move this item to parent */
		data.p = ci;
		ci->hi = HacheTableAdd(io->base->cache, 
				       ci->hi->key, ci->hi->key_len,
				       data, NULL);
		HacheTableIncRef(ci->hi->h, ci->hi);

		/* Remove from this hache */
		HacheTableRemove(io->cache, ci->hi->key, ci->hi->key_len, 0);
	    }
	}

	return cache_flush(io->base);
    }

    to_flush = ArrayCreate(sizeof(cached_item *), 8192);

    /* Identify the list of items that need flushing */
    //fprintf(stderr, "\n");
    for (hi = h->in_use; hi; hi = hi->in_use_next) {
	cached_item *ci = hi->data.p;
	//	fprintf(stderr, "Item %d, type %d, updated %d\n",
	//		ci->view, ci->type, ci->updated);
	if (ci->updated) {
	    ARR(cached_item *, to_flush, nflush++) = ci;
	}
    }

    /* Sort them by view number, which is likely to be approx on-disk order */
    qsort(ArrayBase(cached_item *, to_flush), nflush, sizeof(cached_item *),
	  qsort_ci_view);


    io->iface->lock(io->dbh);

    /* Flush them out */
    for (i = 0; i < nflush; i++) {
	cached_item *ci = arr(cached_item *, to_flush, i);
	//	struct timeval tp1, tp2;
	//	long d1, d2;
	//	gettimeofday(&tp1, NULL);

	write_counts[ci->type]++;

	switch (ci->type) {
	case GT_Database:
	    ret = io->iface->database.write(io->dbh, ci);
	    break;

	case GT_Contig:
	    ret = contig_write(io, ci);
	    break;

	case GT_Seq:
	    ret = seq_write(io, ci);
	    break;

	case GT_SeqBlock:
	    ret = seq_block_write(io, ci);
	    break;

	case GT_Bin:
	    ret = bin_write(io, ci);
	    break;

	case GT_Track:
	    ret = track_write(io, ci);
	    break;

	case GT_RecArray:
	    ret = array_write(io, ci);
	    break;

	case GT_AnnoEle:
	    ret = anno_ele_write(io, ci);
	    break;

	case GT_Anno:
	    ret = anno_write(io, ci);
	    break;

	case GT_Library:
	    ret = library_write(io, ci);
	    break;

	default:
	    fprintf(stderr, "Unable to write object of type %d\n",
		    ci->type);
	    ret = -1;
	}

	//	gettimeofday(&tp2, NULL);
	//	d1 = (tp2.tv_sec - tp1.tv_sec)*1000000 + tp2.tv_usec - tp1.tv_usec;

	if (ret == 0) {
	    ci->updated = 0;
	    HacheTableDecRef(io->cache, ci->hi);
	}

	//	gettimeofday(&tp1, NULL);
	//	d2 = (tp1.tv_sec - tp2.tv_sec)*1000000 + tp1.tv_usec - tp2.tv_usec;
	//	printf("Flush %d type %d, time %ld %ld\n",
	//	       i, ci->type, d1, d2);
    }

    ArrayDestroy(to_flush);

    io->iface->unlock(io->dbh);
    io->iface->commit(io->dbh);

    //printf(">>> flush done <<<\n");
    //HacheTableRefInfo(io->cache, stdout);

#if 0
    printf("\nType       Load\tUnload\tWrite\n");
    printf("----------------------------------\n");
    printf("RecArray     %d\t%d\t%d\n",
	   load_counts[GT_RecArray],
	   unload_counts[GT_RecArray],
	   write_counts[GT_RecArray]);
    printf("Bin          %d\t%d\t%d\n",
	   load_counts[GT_Bin],
	   unload_counts[GT_Bin],
	   write_counts[GT_Bin]);
    printf("Range        %d\t%d\t%d\n",
	   load_counts[GT_Range],
	   unload_counts[GT_Range],
	   write_counts[GT_Range]);
    printf("BTree        %d\t%d\t%d\n",
	   load_counts[GT_BTree],
	   unload_counts[GT_BTree],
	   write_counts[GT_BTree]);
    printf("Database     %d\t%d\t%d\n",
	   load_counts[GT_Database],
	   unload_counts[GT_Database],
	   write_counts[GT_Database]);
    printf("Contig       %d\t%d\t%d\n",
	   load_counts[GT_Contig],
	   unload_counts[GT_Contig],
	   write_counts[GT_Contig]);
    printf("Seq          %d\t%d\t%d\n",
	   load_counts[GT_Seq],
	   unload_counts[GT_Seq],
	   write_counts[GT_Seq]);
    printf("DNASource    %d\t%d\t%d\n",
	   load_counts[GT_DNASource],
	   unload_counts[GT_DNASource],
	   write_counts[GT_DNASource]);
    printf("Track        %d\t%d\t%d\n",
	   load_counts[GT_Track],
	   unload_counts[GT_Track],
	   write_counts[GT_Track]);
    memset(load_counts, 0, 100*sizeof(int));
    memset(unload_counts, 0, 100*sizeof(int));
    memset(write_counts, 0, 100*sizeof(int));
#endif

    return ret;
}

/*
 * Returns true or false depending on whether the cache has unsaved data in
 * it.
 */
int cache_updated(GapIO *io) {
    int i;
    HacheTable *h = io->cache;

    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    cached_item *ci = hi->data.p;
	    if (ci->updated) {
		return 1;
	    }
	}
    }

    return 0;
}

/*
 * Loads an in item into the cache (if not present) and returns it.
 * The query parameters are the object type (GT_*) and the record number.
 * For "io overlays" this may just return the object in the base io
 * instead.
 *
 * Returns a pointer to the object on success
 *         NULL on failure
 */
void *cache_search(GapIO *io, int type, GRec rec) {
    int sub_rec = 0;
    int otype = type, orec = rec;
    cache_key_t k;
    HacheItem *hi;
    
    if (type == GT_Seq) {
	sub_rec = rec & (SEQ_BLOCK_SZ-1);
	rec >>= SEQ_BLOCK_BITS;
	type = GT_SeqBlock;
	//	printf("Converting seq %d to seq_block %d,%d\n", orec, rec, sub_rec);
    }

    k = construct_key(rec, type);
    hi = HacheTableQuery(io->cache, (char *)&k, sizeof(k));

    /* Pass one layer up if we're an overlay on top of another GapIO */
    if (!hi && io->base) {
	return cache_search(io->base, otype, orec);
    } else if (!hi) {
	/* Otherwise if it's not found, force a load */
	hi = HacheTableSearch(io->cache, (char *)&k, sizeof(k));
    }

    if (!hi)
	return NULL;

    if (otype == GT_Seq) {
	seq_block_t *b = (seq_block_t *)&((cached_item *)hi->data.p)->data;
	HacheData hd;

	/*
	 * If this is a child I/O then it's possible this block has partial
	 * data. Hence we look both here and also the parent I/O.
	 */
	if (!b->seq[sub_rec] && io->base) {
	    return cache_search(io->base, otype, orec);
	} else {
	    return b->seq[sub_rec];
	}

	/*
	k = construct_key(b->rec[sub_rec], GT_Seq);
	if (b->seq[sub_rec]) {
	    //	    printf("Adding seq %d to cache\n", b->rec[sub_rec]);
	    hd.p = ci_ptr(b->seq[sub_rec]);
	    hi = HacheTableAdd(io->cache, (char *)&k, sizeof(k), hd, NULL);
	    ci_ptr(b->seq[sub_rec])->hi = hi;
	    return b->seq[sub_rec];
	} else {
	    return NULL;
	}
	*/
    }

    return &((cached_item *)hi->data.p)->data;
}

/*
 * Creates a new item.
 */
int cache_item_create(GapIO *io, int type, void *from) {
    static int brec = 0;
    static int sub_rec = SEQ_BLOCK_SZ;
    seq_block_t *b;

    if (type != GT_Seq) {
	fprintf(stderr,
		"cache_item_create only implemented for GT_Seq right now\n");
	return -1;
    }

    if (sub_rec == SEQ_BLOCK_SZ) {
	sub_rec = 0;
	brec = io->iface->seq_block.create(io->dbh, NULL);
    }

    b = (seq_block_t *)cache_search(io, GT_SeqBlock, brec);

    /* Start new blocks if they contain too much data too */
    if (b->est_size > 150000) {
	//printf("New sub block after %d/%d seqs\n", sub_rec, SEQ_BLOCK_SZ);
	sub_rec = 0;
	brec = io->iface->seq_block.create(io->dbh, NULL);
	b = (seq_block_t *)cache_search(io, GT_SeqBlock, brec);
    }

    b = cache_rw(io, b);
    b->rec[sub_rec] = (brec << SEQ_BLOCK_BITS) + sub_rec;

    /* FIXME: move this somewhere sensible */
    {
	seq_t *s, *f = (seq_t *)from;
	int slen = sizeof(seq_t) + f->name_len+1 + f->trace_name_len+1 +
	    f->alignment_len+1 +
	    ABS(f->len) * (1 + (f->format == SEQ_FORMAT_CNF4 ? 4 : 1));
	cached_item *ci = cache_new(GT_Seq, 0, 0, NULL, slen);

	s = (seq_t *)&ci->data;
	*s = *f;

	s->name = (char *)&s->data;
	strcpy(s->name, f->name ? f->name : "");
	s->name_len = strlen(s->name);

	s->trace_name = s->name + s->name_len + 1;
	strcpy(s->trace_name, f->trace_name ? f->trace_name : "");
	s->trace_name_len = strlen(s->trace_name);

	s->alignment = s->trace_name + s->trace_name_len + 1;
	strcpy(s->alignment, f->alignment ? f->alignment : "");
	s->alignment_len = strlen(s->alignment);

	s->seq = s->alignment + s->alignment_len + 1;
	memcpy(s->seq, f->seq, ABS(f->len));

	s->conf = s->seq + ABS(s->len);
	memcpy(s->conf, f->conf, ABS(f->len)*
	       (f->format == SEQ_FORMAT_CNF4 ? 4 : 1));

	s->block = b;
	s->idx = sub_rec;
	b->seq[sub_rec] = s;
	b->est_size += 2*ABS(s->len) + s->name_len + s->alignment_len +
	    s->trace_name_len + 15;
    }

    //printf("Create seq %d,%d => %d\n", brec, sub_rec, b->rec[sub_rec]);

    return (brec << SEQ_BLOCK_BITS) + sub_rec++;
}

/*
 * Obtains the parent cached_item for this sub-record, or returns this
 * record if there is no parent.
 *
 * This is used for going from sequences to sequence-blocks, the latter of
 * which is the actual granularity for cache storage and locking.
 */
cached_item *cache_master(cached_item *ci) {
    seq_t *s;

    if (!ci || ci->type != GT_Seq)
	return ci;

    s = (seq_t *)&ci->data;
    if (!s->block)
	return ci;

    return ci_ptr(s->block);
}
void cache_incr(GapIO *io, void *data) {
    cached_item *ci = cache_master(ci_ptr(data));
    HacheTableIncRef(ci->hi->h, ci->hi);
}

void cache_decr(GapIO *io, void *data) {
    cached_item *ci = cache_master(ci_ptr(data));
    HacheTableDecRef(ci->hi->h, ci->hi);
}

/*
 * Takes an item from a base io and duplicates into to this io.
 * This is tricky to get right, but basically we need to have the same
 * underlying "View" as our client still only has one connection to the
 * database, but we produce a separate cached_item with duplicated data
 * so it can safely be modified.
 *
 * NB: Once we flush the modified cached_item back we will need to migrate
 * the changes over to the corresponding cached_item in the io->base
 * cache.
 */
cached_item *cache_dup(GapIO *io, cached_item *sub_ci) {
    cached_item *ci = cache_master(sub_ci);
    HacheItem *hi_old = ci->hi;
    HacheItem *hi_new;
    cached_item *ci_new = NULL;

    /* Force load into this io, if not loaded already */
    hi_new = HacheTableQuery(io->cache, hi_old->key, hi_old->key_len);
    if (!hi_new) {
	HacheData hd;

	printf("Cache_dup ci type %d/%d, sub_ci type %d/%d\n", 
	       ci->type, ci->rec, sub_ci->type, sub_ci->rec);

	/* Duplicate the cached_item into this io, keeping the same view */
	ci_new = (cached_item *)malloc(sizeof(*ci) + ci->data_size);
	memcpy(ci_new, ci, sizeof(*ci) + ci->data_size);
	hd.p = ci_new;
	hi_new = HacheTableAdd(io->cache, hi_old->key, hi_old->key_len,
			       hd, NULL);
	ci_new->hi = hi_new;

	switch(ci_new->type) {
	case GT_SeqBlock: {
	    int i;
	    /* Duplicate our arrays */
	    seq_block_t *ob = (seq_block_t *)&ci->data;
	    seq_block_t *b = (seq_block_t *)&ci_new->data;
	    memcpy(b->rec, ob->rec, SEQ_BLOCK_SZ * sizeof(*b->rec));

	    /*
	     * When duplicating a seq_block we know that we're in a child
	     * I/O rather than the master one. In this case we just leave
	     * our seq entries as NULL and use copy-on-write semantics for 
	     * migrating the sequences as and when desired.
	     */
	    for (i = 0; i < SEQ_BLOCK_SZ; i++) {
		b->seq[i] = NULL;
	    }
	    break;
	}

	case GT_Seq: {
	    /* reset internal name/seq/conf pointers */
	    seq_t *os = (seq_t *)&ci->data;
	    seq_t *s  = (seq_t *)&ci_new->data;
	    s->name = (char *)&s->data;
	    s->trace_name = s->name + s->name_len + 1;
	    s->alignment = s->trace_name + s->trace_name_len + 1;
	    s->seq = s->alignment + s->alignment_len + 1;
	    s->conf = s->seq + ABS(s->len);

	    /* Duplicate the annotation Array */
	    /* FIXME: update to new format, whatever that is! */
	    if (s->anno) {
		s->anno = ArrayCreate(sizeof(int), ArrayMax(os->anno));
		memcpy(ArrayBase(int, s->anno),
		       ArrayBase(int, os->anno),
		       ArrayMax(os->anno) * sizeof(int));
	    }
	    break;
	}

	case GT_Bin: {
	    /* Duplicate the Array */
	    bin_index_t *ob = (bin_index_t *)&ci->data;
	    bin_index_t *nb = (bin_index_t *)&ci_new->data;
	    if (ob->rng) {
		nb->rng = ArrayCreate(sizeof(range_t), ArrayMax(ob->rng));
		ArrayMax(nb->rng) = ArrayMax(ob->rng);
		memcpy(ArrayBase(range_t, nb->rng),
		       ArrayBase(range_t, ob->rng),
		       ArrayMax(ob->rng) * sizeof(range_t));
	    }
	    if (ob->track) {
		nb->track = ArrayCreate(sizeof(track_t), ArrayMax(ob->track));
		ArrayMax(nb->track) = ArrayMax(ob->track);
		memcpy(ArrayBase(GBinTrack, nb->track),
		       ArrayBase(GBinTrack, ob->track),
		       ArrayMax(ob->track) * sizeof(GBinTrack));
	    }
	    break;
	}

	case GT_Track: {
	    /* Duplicate the Array */
	    track_t *ob = (track_t *)&ci->data;
	    track_t *nb = (track_t *)&ci_new->data;
	    if (ob->data) {
		nb->data = ArrayCreate(nb->item_size, ArrayMax(ob->data));
		ArrayMax(nb->data) = ArrayMax(ob->data);
		memcpy(ArrayBase(char, nb->data),
		       ArrayBase(char, ob->data),
		       ArrayMax(ob->data) * nb->item_size);
	    }
	    break;
	}

	case GT_AnnoEle: {
	    /* reset internal comment pointer */
	    anno_ele_t *oe = (anno_ele_t *)&ci->data;
	    anno_ele_t *ne = (anno_ele_t *)&ci_new->data;
	    if (oe->comment)
		ne->comment = (char *)&ne->data;
	    else
		ne->comment = NULL;
	    break;
	}

	case GT_Anno: {
	    anno_t *oa = (anno_t *)&ci->data;
	    anno_t *na = (anno_t *)&ci_new->data;
	    
	    na->key   = oa->key   ? strdup(oa->key)   : NULL;
	    na->value = oa->value ? strdup(oa->value) : NULL;
	    if (oa->ele) {
		na->ele = ArrayCreate(sizeof(int), ArrayMax(oa->ele));
		ArrayMax(na->ele) = ArrayMax(oa->ele);
		memcpy(ArrayBase(char, na->ele),
		       ArrayBase(char, oa->ele),
		       ArrayMax(na->ele) * sizeof(int));
	    }
	    break;
	}

	case GT_Library: {
	    puts("FIXME: implement library_dup");
	    break;
	}

	}
    } else {
	/* Already existed, so just return it */
	ci_new = (cached_item *)hi_new->data.p;
    }

    /*
     * Partial duplicate of above.
     * When ci and sub_ci differ (ie ci is a sub-record of a block) then we
     * also need to populate the entry sub_ci within our duplicated ci.
     */
    if (ci != sub_ci) {
	cached_item *sub_new;

	/* Duplicate */
	switch (sub_ci->type) {
	case GT_Seq: {
	    seq_block_t *b = (seq_block_t *)&ci_new->data;
	    seq_t *os  = (seq_t *)&sub_ci->data, *s;
	    int sub_rec = os->rec & (SEQ_BLOCK_SZ-1);

	    /* Already duplicated? */
	    if (b->seq[os->idx]) {
		sub_new = sub_ci;
		break;
	    }
	    
	    printf("Duplicate seq %d in block %d\n", sub_rec, ci->rec);

	    sub_new = (cached_item *)malloc(sizeof(*ci) + sub_ci->data_size);
	    memcpy(sub_new, sub_ci, sizeof(*ci) + sub_ci->data_size);
	    s = (seq_t *)&sub_new->data;

	    s->name = (char *)&s->data;
	    s->trace_name = s->name + s->name_len + 1;
	    s->alignment = s->trace_name + s->trace_name_len + 1;
	    s->seq = s->alignment + s->alignment_len + 1;
	    s->conf = s->seq + ABS(s->len);

	    /* Duplicate the annotation Array */
	    /* FIXME: update to new format, whatever that is! */
	    if (s->anno) {
		s->anno = ArrayCreate(sizeof(int), ArrayMax(os->anno));
		memcpy(ArrayBase(int, s->anno),
		       ArrayBase(int, os->anno),
		       ArrayMax(os->anno) * sizeof(int));
	    }

	    s->block = b;
	    b->seq[s->idx] = s;

	    break;
	}
	}

	ci_new = sub_new;
    }

    return ci_new;
}

/*
 * Locks a cached item for read-write access, returning a new pointer.
 * For 'derived' GapIOs this locks the original but returns the duplicated
 * (copy on write) version. Otherwise it locks in situ and returns the
 * same address.
 *
 * Returns new ptr on success (maybe same or different to data)
 *         NULL on failure
 */
void *cache_rw(GapIO *io, void *data) {
    cached_item *ci = ci_ptr(data);
    cached_item *mi = cache_master(ci);

    if (io->base) {
	ci = cache_dup(io, ci);
	mi = cache_master(ci);
	/*
	if (ci->type == GT_Seq) {
	    seq_t *s = (seq_t *)&ci->data;
	    seq_block_t *b = (seq_block_t *)&mi->data;
	    ci = ci_ptr(b->seq[s->idx]);
	}
	*/
	data = &ci->data;
    }

    /* Ensure it's locked RW */
    if (mi->lock_mode < G_LOCK_RW) {
	if (-1 == cache_upgrade(io, mi, G_LOCK_RW)) {
	    fprintf(stderr, "lock denied for rec %d\n", mi->rec);
	    return NULL;
	}
    }

    /* Also set updated flag to 1 and bump ref count if needed */
    if (!mi->updated) {
	mi->updated = 1;
	HacheTableIncRef(mi->hi->h, mi->hi);
    }

    return data;
}
