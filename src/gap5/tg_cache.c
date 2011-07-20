#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>

#ifdef CACHE_REF_DEBUG
/*  So tg_gio.h doesn't redefine prototypes */
#   define WAS_CACHE_REF_DEBUG
#   undef CACHE_REF_DEBUG
#endif
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
    tg_rec rec;
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
cached_item *cache_new(int type, tg_rec rec, GView v,
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
    ci->forgetme = 0;
    ci->data_size = e_len;

    return ci;
}

static int chksum(cached_item *ci) {
    return HacheTcl((uint8_t *)&ci->data, ci->data_size);
}

/*
 * Frees a cached item, for use after appropriate <TYPE>_unload function 
 * has been called. This is usually nothing more than a free() call, but
 * centralised here to provide extra debugging functionality when needed.
 */
static void cache_free(cached_item *ci) {
    if (ci->chk_sum && chksum(ci) != ci->chk_sum && ci->lock_mode < G_LOCK_RW){
	printf("Chksum differs on ci for rec %"PRIrec"\n", ci->rec);
	abort();
    }

#ifdef WAS_CACHE_REF_DEBUG
    /*
     * Also memset the data so it's blatantly obvious when we try to use it.
     * This helps to detect possible reference count issues where we should
     * have done a cache_incr() in order to prevent something from being
     * purged.
     */
    memset(ci, 'z', ci->data_size + sizeof(*ci));
#endif

    free(ci);
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
static void seq_unload(GapIO *io, cached_item *ci, int unlock) {
    seq_t *s = (seq_t *)&ci->data;
    if (unlock)
	io->iface->seq.unlock(io->dbh, ci->view);

    if (s->anno)
	ArrayDestroy(s->anno);

    cache_free(ci);
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
static void seq_block_unload(GapIO *io, cached_item *ci, int unlock) {
    int i;
    seq_block_t *b = (seq_block_t *)&ci->data;

    if (unlock)
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
    cache_free(ci);
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
static void track_unload(GapIO *io, cached_item *ci, int unlock) {
    track_t *track = (track_t *)&ci->data;

    if (unlock)
	io->iface->track.unlock(io->dbh, ci->view);

    if (track->data)
	ArrayDestroy(track->data);
    cache_free(ci);
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
static void bin_unload(GapIO *io, cached_item *ci, int unlock) {
    bin_index_t *bin = (bin_index_t *)&ci->data;

    if (ci->forgetme)
	io->iface->bin.destroy(io->dbh, ci->rec, ci->view);

    if (unlock)
	io->iface->bin.unlock(io->dbh, ci->view);

    if (bin->rng)
	ArrayDestroy(bin->rng);
    if (bin->track)
	ArrayDestroy(bin->track);
    cache_free(ci);
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


static void contig_unload(GapIO *io, cached_item *ci, int unlock) {
    if (ci->forgetme)
	io->iface->contig.destroy(io->dbh, ci->rec, ci->view);

    if (unlock)
	io->iface->contig.unlock(io->dbh, ci->view);
    cache_free(ci);
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

static void array_unload(GapIO *io, cached_item *ci, int unlock) {
    Array ar = (Array)&ci->data;
    //ArrayDestroy(ar);
    if (ar->base)
	free(ar->base);

    if (unlock)
	io->iface->seq.unlock(io->dbh, ci->view);
    cache_free(ci);
}

static void database_unload(GapIO *io, cached_item *ci, int unlock) {
    if (unlock)
	io->iface->database.unlock(io->dbh, ci->view);
    cache_free(ci);
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
static int anno_ele_block_write(GapIO *io, cached_item *ci) {
    return io->iface->anno_ele_block.write(io->dbh, ci);
}

static void anno_ele_unload(GapIO *io, cached_item *ci, int unlock) {
    if (unlock)
	io->iface->anno_ele.unlock(io->dbh, ci->view);
    cache_free(ci);
}
static void anno_unload(GapIO *io, cached_item *ci, int unlock) {
    if (unlock)
	io->iface->anno.unlock(io->dbh, ci->view);
    cache_free(ci);
}
static void anno_ele_block_unload(GapIO *io, cached_item *ci, int unlock) {
    int i;
    anno_ele_block_t *b = (anno_ele_block_t *)&ci->data;

    if (unlock)
	io->iface->anno_ele_block.unlock(io->dbh, ci->view);

    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
	if (b->ae[i]) {
	    cached_item *si = ci_ptr(b->ae[i]);
	    if (si)
		free(si);
	}
    }
    cache_free(ci);
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

static void library_unload(GapIO *io, cached_item *ci, int unlock) {
    if (unlock)
	io->iface->library.unlock(io->dbh, ci->view);

    /* Disable checking as even in r/o mode we update insert_size on the fly */
    ci->chk_sum = 0;

    cache_free(ci);
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

    new->data_size = size;

    if (ci == new)
	return item;

    if (new->hi) {
	assert(new->hi->data.p == ci);
	new->hi->data.p = new;
    }

    switch (new->type) {
    case GT_Seq: {
	seq_t *s = (seq_t *)&new->data;
	s->block->seq[s->idx] = s;
	sequence_reset_ptr(s);
	break;
    }

    case GT_AnnoEle: {
	anno_ele_t *e = (anno_ele_t *)&new->data;
	e->block->ae[e->idx] = e;
	break;
    }
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
	break;

    case GT_SeqBlock:
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

    case GT_AnnoEleBlock:
	ci = io->iface->anno_ele_block.read(io->dbh, k->rec);
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

    ci->chk_sum = chksum(ci);

    return &hd;
}

/*
 * Ensure key is initialised correctly, making sure that it is blank
 * in the gaps between structure elements and the padding at the end.
 */
static void construct_key(tg_rec rec, int type, cache_key_t *k) {
    memset(k, 0, sizeof(*k));
    k->rec = rec;
    k->type = type;
}

/* Callback from Hache */
static void cache_unload(void *clientdata, HacheData hd) {
    GapIO *io = (GapIO *)clientdata;
    cached_item *ci = hd.p;
    int unlock = 1;

    //    printf("Cache unload %d\n", ci->rec);

    assert(io->base || ci->updated == 0);

    /*
     * If we're a child I/O, do not unlock the view of this item if the
     * base I/O also has a view open on it.
     */
    if (io->base) {
	HacheItem *hi_base;
	cache_key_t k;

	construct_key(ci->rec, ci->type, &k);
	hi_base = HacheTableQuery(io->base->cache, (char *)&k, sizeof(k));

	if (hi_base)
	    unlock = 0;
    }

    unload_counts[ci->type]++;

    switch (ci->type) {
    case GT_Seq:
	seq_unload(io, ci, unlock);
	break;

    case GT_SeqBlock:
	seq_block_unload(io, ci, unlock);
	break;

    case GT_Bin:
	bin_unload(io, ci, unlock);
	break;

    case GT_Track:
	track_unload(io, ci, unlock);
	break;

    case GT_Contig:
	contig_unload(io, ci, unlock);
	break;

    case GT_RecArray:
	array_unload(io, ci, unlock);
	break;

    case GT_Database:
	database_unload(io, ci, unlock);
	break;

    case GT_AnnoEle:
	anno_ele_unload(io, ci, unlock);
	break;

    case GT_AnnoEleBlock:
	anno_ele_block_unload(io, ci, unlock);
	break;

    case GT_Anno:
	anno_unload(io, ci, unlock);
	break;

    case GT_Library:
	library_unload(io, ci, unlock);
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
#ifdef WAS_CACHE_REF_DEBUG
    /* Test smaller cache to stress-test ref counting bugs */
    if (NULL == (h = HacheTableCreate(64, HASH_DYNAMIC_SIZE|HASH_OWN_KEYS)))
	return -1;
#else
    if (NULL == (h = HacheTableCreate(1024, HASH_DYNAMIC_SIZE|HASH_OWN_KEYS)))
	return -1;
#endif
    h->name = "tg_cache";

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

    HacheTableStats(h, stderr);

    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    cache_unload(io, hi->data);
	}
    }

    HacheTableDestroy(io->cache, 0);
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

    case GT_AnnoEleBlock:
	ret = io->iface->anno_ele_block.upgrade(io->dbh, ci->view, mode);
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
void *cache_lock(GapIO *io, int type, tg_rec rec, int mode) {
    cache_key_t key;
    cached_item *ci;
    HacheTable *h = io->cache;
    HacheItem *hi;

    construct_key(rec, type, &key);    
    hi = HacheTableSearch(h, (char *)&key, sizeof(key));
    
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

    //printf("\n>>> cache flush <<<\n");
    //HacheTableRefInfo(io->cache, stdout);

    /*
     * If this is a derived io then we need to pass these items up to
     * our parent, remove them from this io and call cache_flush on the
     * parent instead (possibly recursively).
     */
    if (io->base) {
	for (i = 0; i < h->nbuckets; i++) {
	    HacheItem *next;
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
		    int j;

		    htmp = HacheTableSearch(io->base->cache,
					    hi->key, hi->key_len);
		    bo = (seq_block_t *)&((cached_item *)htmp->data.p)->data;

		    for (j = 0; j < SEQ_BLOCK_SZ; j++) {
			if (!bn->seq[j]) {
			    bn->seq[j] = bo->seq[j];
			    if (bn->seq[j])
				bn->seq[j]->block = bn;
			}
		    }

		    break;
		}

		case GT_AnnoEleBlock: {
		    HacheItem *htmp;
		    anno_ele_block_t *bn = (anno_ele_block_t *)&ci->data;
		    anno_ele_block_t *bo;
		    int j;

		    htmp = HacheTableSearch(io->base->cache,
					    hi->key, hi->key_len);
		    bo = (anno_ele_block_t *)&((cached_item *)htmp->data.p)->data;

		    for (j = 0; j < ANNO_ELE_BLOCK_SZ; j++) {
			if (!bn->ae[j]) {
			    bn->ae[j] = bo->ae[j];
			    if (bn->ae[j])
				bn->ae[j]->block = bn;
			}
		    }
		    break;
		}
		}



		/*
		 * Synchronise references.
		 *
		 * When we use the child io all reads look in the child first
		 * and then base if not found in child. Once we run cache_rw
		 * on the child, we perform a copy-on-write system (cache_dup)
		 * and create a duplicate of the object in the child cache.
		 *
		 * FIXME: at this stage, if anything else decides to hold a
		 * long-term reference to an object then it could be
		 * overwritten when we flush the child IO. We maybe need
		 * some registration system, or at least a standard protocol
		 * for how to program long-duration windows. (Don't trust
		 * cache_incr.)
		 */
		/* FIXME: not implemented! Need to think on this more. */

		/* Purge from parent */
		HacheTableRemove(io->base->cache,
				 ci->hi->key, ci->hi->key_len, 0);

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
	//printf("View %d, type %d, rec %d, updated %d, forgetme %d\n",
	//	ci->view, ci->type, ci->rec, ci->updated, ci->forgetme);
	if (ci->updated) {
	    ARR(cached_item *, to_flush, nflush++) = ci;
	    //printf("To flush %d, updated\n", ci->rec);
	} else {
	    if (ci->forgetme) {
		ARR(cached_item *, to_flush, nflush++) = ci;
	    }
	}
    }


#if 0
    /*
     * DEBUG: double check all items to make sure that we're not going
     * to miss something which should be flushed but is somehow not
     * labelled as in_use.
     */
    {
	HacheIter *iter = HacheTableIterCreate();
	while (hi = HacheTableIterNext(h, iter)) {
	    cached_item *ci = hi->data.p;

	    if (ci->updated || ci->forgetme) {
		int i;
		for (i = 0; i < nflush; i++) {
		    if (arr(cached_item *, to_flush, i) == ci)
			break;
		}
		assert(i < nflush);
	    }
	}
	HacheTableIterDestroy(iter);
    }
#endif


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

	/*
	 * Should we delay this until we've written data, so that these
	 * records aren't reused during this self-same flush? It may
	 * make it impossible to roll back as it currently stands, but
	 * maybe this is impossible already!? (We have no tool to attempt
	 * it anyway currently.)
	 */
	if (ci->forgetme) {
	    ci->updated = 0;
	    HacheTableDecRef(io->cache, ci->hi);
	    HacheTableDel(io->cache, ci->hi, 1);
	    continue;
	}

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

	case GT_AnnoEleBlock:
	    ret = anno_ele_block_write(io, ci);
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
    printf("Library    %d\t%d\t%d\n",
	   load_counts[GT_Library],
	   unload_counts[GT_Library],
	   write_counts[GT_Library]);
    printf("Track        %d\t%d\t%d\n",
	   load_counts[GT_Track],
	   unload_counts[GT_Track],
	   write_counts[GT_Track]);
    memset(load_counts, 0, 100*sizeof(int));
    memset(unload_counts, 0, 100*sizeof(int));
    memset(write_counts, 0, 100*sizeof(int));
#endif

#ifdef WAS_CACHE_REF_DEBUG
    /* Count updated & locked records, to check for cache ref-count leaks */
    {
	int total[100];
	int updated[100];
	int ref_count[100];
	int in_use[100];
	size_t sz[100], tsz, tc;

	memset(total,     0, 100*sizeof(int));
	memset(updated,   0, 100*sizeof(int));
	memset(ref_count, 0, 100*sizeof(int));
	memset(in_use,    0, 100*sizeof(int));
	memset(sz,        0, 100*sizeof(size_t));

	for (i = 0; i < h->nbuckets; i++) {
	    HacheItem *next;
	    for (hi = h->bucket[i]; hi; hi = next) {
		HacheData data;
		cached_item *ci = hi->data.p;

		next = hi->next;

		if (ci->type != GT_Library &&
		    chksum(ci) != ci->chk_sum && ci->lock_mode < G_LOCK_RW) {
		    printf("Chksum differs on ci for rec %"PRIrec"\n",
			   ci->rec);
		    abort();
		}

		if (ci->updated) updated[ci->type]++;
		if (hi->ref_count) ref_count[ci->type]++;
		if (hi->in_use_prev || hi->in_use_next || h->in_use == hi)
		    in_use[ci->type]++;
		total[ci->type]++;

		sz[ci->type] += sizeof(*ci) + ci->data_size + 8;

		switch(ci->type) {
		    int j;
		    seq_block_t *sb;
		    anno_ele_block_t *eb;
		    cached_item *si;

		case GT_SeqBlock:
		    sb = (seq_block_t *)&ci->data;
		    for (j = 0; j < SEQ_BLOCK_SZ; j++) {
			seq_t *s = sb->seq[j];
			if (s) {
			    si = ci_ptr(s);
			    sz[ci->type] += sizeof(*si) + si->data_size + 8;
			}
		    }
		    break;

		case GT_AnnoEleBlock:
		    eb = (anno_ele_block_t *)&ci->data;
		    for (j = 0; j < SEQ_BLOCK_SZ; j++) {
			anno_ele_t *e = eb->ae[j];
			if (e) {
			    si = ci_ptr(e);
			    sz[ci->type] += sizeof(*si) + si->data_size + 8;
			}
		    }
		    break;
		}
	    }
	}

	for (tc = tsz = i = 0; i < 100; i++) {
	    tsz += sz[i];
	    tc  += total[i];
	}
	printf("After flush: %ld items, %ld bytes in cache\n", tc, tsz);
	for (i = 0; i < 100; i++) {
	    char *s[] = {
		"", "", "", "RecArray", "", "Bin", "Range", "BTree",
		"", "", "", "", "", "", "", "", "Database", "Contig",
		"Seq", "Library", "Track", "AnnoEle", "Anno",
		"SeqBlock", "AnnoEleBlock"
	    };
	    if (!total[i]) continue;
	    printf("  Type %d %s\n", i, s[i]);
	    printf("    Size      = %ld (%ld/item)\n", sz[i], sz[i]/total[i]);
	    printf("    Total rec = %d\n", total[i]);
	    printf("    locked    = %d\n", ref_count[i]);
	    if (in_use[i] != ref_count[i])
		printf("****'in_use'  = %d\n", in_use[i]);
	    if (updated[i])
		printf("****updated   = %d\n", updated[i]);
	}

	cache_ref_debug_dump();
    }
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
 * Dumps information about the cache, for debugging purposes.
 */
static int interested_rec = -1;
void cache_dump(GapIO *io) {
    int i;
    HacheTable *h = io->cache;
    int nused = 0, nlocked = 0, nupdated = 0;

    printf("\nCheck for io = %p (%s)\n", io,
	   io->base ? "child" : "base");

    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    cached_item *ci = hi->data.p;

	    if (interested_rec != -1 && ci->rec != interested_rec)
		continue;

	    printf("  rec=%"PRIrec"\tv=%d\tlock=%d\ttype=%d\tci=%p\trc=%d\n",
		   ci->rec, ci->view, ci->lock_mode, ci->type, ci,
		   hi->ref_count);

	    nused++;
	    if (ci->lock_mode >= G_LOCK_RW)
		nlocked++;
	    if (ci->updated)
		nupdated++;

	    assert(ci->updated == 0 || ci->lock_mode >= G_LOCK_RW);
	    assert(ci->hi == hi);
	    assert(hi->h == io->cache);
	}
    }
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
void *cache_search(GapIO *io, int type, tg_rec rec) {
    int sub_rec = 0;
    int otype = type;
    tg_rec orec = rec;
    cache_key_t k;
    HacheItem *hi;
    
    switch (type) {
    case GT_Seq:
	sub_rec = rec & (SEQ_BLOCK_SZ-1);
	rec >>= SEQ_BLOCK_BITS;
	type = GT_SeqBlock;
	break;

    case GT_AnnoEle:
	sub_rec = rec & (ANNO_ELE_BLOCK_SZ-1);
	rec >>= ANNO_ELE_BLOCK_BITS;
	type = GT_AnnoEleBlock;
	break;
    }

    construct_key(rec, type, &k);
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

    switch (otype) {
    case GT_Seq: {
	seq_block_t *b = (seq_block_t *)&((cached_item *)hi->data.p)->data;

	/*
	 * If this is a child I/O then it's possible this block has partial
	 * data. Hence we look both here and also the parent I/O.
	 */
	if (!b->seq[sub_rec] && io->base) {
	    return cache_search(io->base, otype, orec);
	} else {
	    return b->seq[sub_rec];
	}
	break;
    }

    case GT_AnnoEle: {
	anno_ele_block_t *b = (anno_ele_block_t *)&((cached_item *)hi->data.p)->data;

	/*
	 * If this is a child I/O then it's possible this block has partial
	 * data. Hence we look both here and also the parent I/O.
	 */
	if (!b->ae[sub_rec] && io->base) {
	    return cache_search(io->base, otype, orec);
	} else {
	    return b->ae[sub_rec];
	}
	break;
    }
    }

    return &((cached_item *)hi->data.p)->data;
}


/*
 * As per cache_search, but do not load the item if it's not already in the
 * cache.
 *
 * Returns a pointer to the object on success
 *         NULL on failure or if not in cache
 */
void *cache_search_no_load(GapIO *io, int type, tg_rec rec) {
    int sub_rec = 0;
    int otype = type;
    tg_rec orec = rec;
    cache_key_t k;
    HacheItem *hi;
    
    switch (type) {
    case GT_Seq:
	sub_rec = rec & (SEQ_BLOCK_SZ-1);
	rec >>= SEQ_BLOCK_BITS;
	type = GT_SeqBlock;
	break;

    case GT_AnnoEle:
	sub_rec = rec & (ANNO_ELE_BLOCK_SZ-1);
	rec >>= ANNO_ELE_BLOCK_BITS;
	type = GT_AnnoEleBlock;
	break;
    }

    construct_key(rec, type, &k);
    hi = HacheTableQuery(io->cache, (char *)&k, sizeof(k));

    if (!hi && io->base)
	return cache_search_no_load(io->base, otype, orec);

    if (!hi)
	return NULL;

    switch (otype) {
    case GT_Seq:
	{
	    seq_block_t *b = (seq_block_t *)&((cached_item *)hi->data.p)->data;
	    return b->seq[sub_rec];
	}

    case GT_AnnoEle:
	{
	    anno_ele_block_t *b = (anno_ele_block_t *)
		((cached_item *)hi->data.p)->data;
	    return b->ae[sub_rec];
        }
    }

    return &((cached_item *)hi->data.p)->data;
}


/*
 * Returns whether (type,rec) is a legal combination and exists.
 * 1 = yes, 0 = no.
 */
int cache_exists(GapIO *io, int type, int rec) {
    switch (type) {
    case GT_Seq:
	return io->iface->exists(io->dbh, GT_SeqBlock, rec >> SEQ_BLOCK_BITS);

    case GT_AnnoEle:
	return io->iface->exists(io->dbh, GT_AnnoEleBlock,
				 rec >> ANNO_ELE_BLOCK_BITS);

    default:
	return io->iface->exists(io->dbh, type, rec);
    }
}

/*
 * Creates a new seq_t item.
 */
static int cache_item_init_seq(GapIO *io, void *from, tg_rec rec) {
    seq_t *s, *f = (seq_t *)from;
    size_t slen = sizeof(seq_t) + sequence_extra_len(f);
    cached_item *ci = cache_new(GT_Seq, 0, 0, NULL, slen);
    tg_rec brec, sub_rec;
    seq_block_t *b;

    sub_rec = rec & (SEQ_BLOCK_SZ-1);
    brec    = rec >> SEQ_BLOCK_BITS;

    s = (seq_t *)&ci->data;
    if (sequence_copy(s, f) == -1)
	return -1;

    b = (seq_block_t *)cache_search(io, GT_SeqBlock, brec);

    s->rec = rec;
    s->block = b;
    s->idx = (int)sub_rec;
    b->seq[sub_rec] = s;
    b->est_size += 15 + sequence_extra_len(s);

    return 0;
}

static tg_rec cache_item_create_seq(GapIO *io, void *from) {
    static tg_rec brec = 0;
    static tg_rec sub_rec = SEQ_BLOCK_SZ;
    seq_block_t *b;

    if (sub_rec == SEQ_BLOCK_SZ) {
	sub_rec = 0;
	brec = io->iface->seq_block.create(io->dbh, NULL);
    }

    b = (seq_block_t *)cache_search(io, GT_SeqBlock, brec);

    /* Start new blocks if they contain too much data too */
    if (b->est_size > 250000) {
	//printf("New sub block after %d/%d seqs\n", sub_rec, SEQ_BLOCK_SZ);
	sub_rec = 0;
	brec = io->iface->seq_block.create(io->dbh, NULL);
	b = (seq_block_t *)cache_search(io, GT_SeqBlock, brec);
    }

    b = cache_rw(io, b);

    /* FIXME: move this somewhere sensible */
    if (from) {
	if (cache_item_init_seq(io, from, (brec << SEQ_BLOCK_BITS) + sub_rec))
	    return -1;
    }

    //assert(brec < (1<<(31-SEQ_BLOCK_BITS)));

    return (brec << SEQ_BLOCK_BITS) + sub_rec++;
}

/*
 * Creates a new anno_ele_t item.
 */
static int cache_item_init_anno_ele(GapIO *io, void *from, tg_rec rec) {
    anno_ele_t *t, *f = (anno_ele_t *)from;
    int slen = sizeof(anno_ele_t) +
	(f->comment ? strlen(f->comment) : 0)+1;
    cached_item *ci = cache_new(GT_AnnoEle, 0, 0, NULL, slen);
    tg_rec brec, sub_rec;
    anno_ele_block_t *b;

    sub_rec = rec & (ANNO_ELE_BLOCK_SZ-1);
    brec    = rec >> ANNO_ELE_BLOCK_BITS;

    t = (anno_ele_t *)&ci->data;
    *t = *f;
    t->comment = (char *)&t->data;
    strcpy(t->comment, f->comment ? f->comment : "");

    b = (anno_ele_block_t *)cache_search(io, GT_AnnoEleBlock, brec);

    t->rec = (brec << ANNO_ELE_BLOCK_BITS) + sub_rec;
    t->block = b;
    t->idx = (int)sub_rec;
    b->ae[sub_rec] = t;
    b->est_size += strlen(t->comment) + 10;

    return 0;
}

static tg_rec cache_item_create_anno_ele(GapIO *io, void *from) {
    static tg_rec brec = 0;
    static tg_rec sub_rec = ANNO_ELE_BLOCK_SZ;
    anno_ele_block_t *b;

    if (sub_rec == ANNO_ELE_BLOCK_SZ) {
	sub_rec = 0;
	brec = io->iface->anno_ele_block.create(io->dbh, NULL);
    }

    b = (anno_ele_block_t *)cache_search(io, GT_AnnoEleBlock, brec);

    /* Start new blocks if they contain too much data too */
    if (b->est_size > 150000) {
	sub_rec = 0;
	brec = io->iface->anno_ele_block.create(io->dbh, NULL);
	b = (anno_ele_block_t *)cache_search(io, GT_AnnoEleBlock, brec);
    }

    b = cache_rw(io, b);

    if (from) {
	if (cache_item_init_anno_ele(io, from, (brec << ANNO_ELE_BLOCK_BITS) + sub_rec))
	    return -1;
    }

    //assert(brec < (1<<(31-ANNO_ELE_BLOCK_BITS)));

    return (brec << ANNO_ELE_BLOCK_BITS) + sub_rec++;
}

/*
 * Creates a new item.
 *
 * Note this is using io->iface and io->dbh so it directly creates records
 * using the lower level API, even if this is a child IO. (Both io->dbh and
 * io->base->dbh will be the same.)
 *
 * Therefore our child I/O isn't quite the copy-on-write layer we think it
 * is for item creation or destruction. It's sufficient for our current
 * purposes though and the failure mode is harmless - we run the risk of,
 * say, creating a tag in the editor, quitting without saving it, using up
 * a record for the GT_AnnoEleBlock but never writing one.
 */
tg_rec cache_item_create(GapIO *io, int type, void *from) {
    switch (type) {
    case GT_Seq:
	return cache_item_create_seq(io, from);

    case GT_AnnoEle:
	return cache_item_create_anno_ele(io, from);

    default:
	fprintf(stderr,
		"cache_item_create only implemented for GT_Seq/GT_AnnoEle right now\n");
    }
    
    return -1;
}


/*
 * Removes an item from one of the blocked structures.
 * The caller is expected to deallocate the seq/anno record too.
 */
int cache_item_remove(GapIO *io, int type, tg_rec rec) {
    int sub_rec = 0;
    seq_block_t *sb;
    anno_ele_block_t *ab;

    switch (type) {
    case GT_Seq:
	sub_rec = rec & (SEQ_BLOCK_SZ-1);
	rec >>= SEQ_BLOCK_BITS;
	sb = cache_search(io, GT_SeqBlock, rec);
	sb = cache_rw(io, sb);
	sb->seq[sub_rec] = NULL;
	break;

    case GT_AnnoEle:
	sub_rec = rec & (ANNO_ELE_BLOCK_SZ-1);
	rec >>= ANNO_ELE_BLOCK_BITS;
	ab = cache_search(io, GT_AnnoEleBlock, rec);
	ab = cache_rw(io, ab);
	ab->ae[sub_rec] = NULL;
	break;

    default:
	fprintf(stderr, "cache_item_remove only implemented for "
		"GT_Seq/GT_AnnoEle.\n");
	return -1;
    }

    return 0;
}


/*
 * Initialises an item, for use when cache_item_create was used with 
 * from == NULL. This allows us to pre-allocate a record number and then
 * initialise it once we know the full details.
 */
int cache_item_init(GapIO *io, int type, void *from, tg_rec rec) {
    switch (type) {
    case GT_Seq:
	return cache_item_init_seq(io, from, rec);

    case GT_AnnoEle:
	return cache_item_init_anno_ele(io, from, rec);

    default:
	fprintf(stderr,
		"cache_item_init only implemented for GT_Seq/GT_AnnoEle right now\n");
    }
    
    return -1;
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
    anno_ele_t *e;

    if (!ci)
	return NULL;

    switch (ci->type) {
    case GT_Seq:
	s = (seq_t *)&ci->data;
	if (!s->block)
	    return ci;
	else
	    return ci_ptr(s->block);

    case GT_AnnoEle:
	e = (anno_ele_t *)&ci->data;
	if (!e->block)
	    return ci;
	else
	    return ci_ptr(e->block);
    }

    return ci;
}

void cache_incr(GapIO *io, void *data) {
    cached_item *ci = cache_master(ci_ptr(data));
    HacheTableIncRef(ci->hi->h, ci->hi);
}

void cache_decr(GapIO *io, void *data) {
    cached_item *ci = cache_master(ci_ptr(data));
    HacheTableDecRef(ci->hi->h, ci->hi);

    assert(ci->hi->ref_count >= 0);
    assert(ci->updated == 0 || ci->hi->ref_count > 0);
}

static HacheTable *ref_debug = NULL;
void cache_incr_debug(GapIO *io, void *data, char *where) {
    char key[100];
    cached_item *ci = cache_master(ci_ptr(data));
    HacheData hd;

    if (!ref_debug)
	ref_debug = HacheTableCreate(1024, HASH_DYNAMIC_SIZE);

    sprintf(key, "%p-%d", data, ci->hi->ref_count - ci->updated);
    hd.p = strdup(where);
    HacheTableAdd(ref_debug, key, 0, hd, NULL);

    //fprintf(stderr, "INCR %s %s\n", key, where);

    cache_incr(io, data);
}

void cache_decr_debug(GapIO *io, void *data, char *where) {
    char key[100];
    cached_item *ci = cache_master(ci_ptr(data));
    sprintf(key, "%p-%d", data, ci->hi->ref_count-1 - ci->updated);

    if (HacheTableRemove(ref_debug, key, 0, 1) != 0) {
	fprintf(stderr, "Failed to remove %s - not in hash table?\n",
		key);
    }

    //fprintf(stderr, "DECR %s %s\n", key, where);

    cache_decr(io, data);
}

void cache_ref_debug_dump(void) {
    HacheIter *iter = HacheTableIterCreate();
    HacheItem *hi;
    HacheTable *h = HacheTableCreate(16, HASH_DYNAMIC_SIZE);

    while (hi = HacheTableIterNext(ref_debug, iter)) {
	HacheData hd;
	HacheItem *hi2;
	hd.i = 0;
	hi2 = HacheTableAdd(h, hi->data.p, 0, hd, NULL);
	hi2->data.i++;
    }
    HacheTableIterDestroy(iter);

    iter = HacheTableIterCreate();
    while (hi = HacheTableIterNext(h, iter)) {
	printf("%"PRId64"\t%s\n", hi->data.i, hi->key);
    }
    HacheTableIterDestroy(iter);

    HacheTableDestroy(h, 0);
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

    /* Ensure base io copy cannot go out of scope */
    HacheTableIncRef(ci->hi->h, ci->hi);

    /* Force load into this io, if not loaded already */
    hi_new = HacheTableQuery(io->cache, hi_old->key, hi_old->key_len);
    if (!hi_new) {
	HacheData hd;

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
	    seq_block_t *b = (seq_block_t *)&ci_new->data;

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

	    sequence_reset_ptr(s);

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
		memcpy(ArrayBase(bin_track_t, nb->track),
		       ArrayBase(bin_track_t, ob->track),
		       ArrayMax(ob->track) * sizeof(bin_track_t));
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

	case GT_AnnoEleBlock: {
	    int i;
	    /* Duplicate our arrays */
	    anno_ele_block_t *b = (anno_ele_block_t *)&ci_new->data;

	    for (i = 0; i < ANNO_ELE_BLOCK_SZ; i++) {
		b->ae[i] = NULL;
	    }
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
	    //int sub_rec = os->rec & (SEQ_BLOCK_SZ-1);

	    /* Already duplicated? */
	    if (b->seq[os->idx]) {
		sub_new = sub_ci;
		break;
	    }
	    
	    //printf("Duplicate seq %d in block %d\n", sub_rec, ci->rec);

	    sub_new = (cached_item *)malloc(sizeof(*ci) + sub_ci->data_size);
	    memcpy(sub_new, sub_ci, sizeof(*ci) + sub_ci->data_size);
	    s = (seq_t *)&sub_new->data;

	    sequence_reset_ptr(s);

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

	    /* Bump reference count of master again */
	    HacheTableIncRef(ci_new->hi->h, ci_new->hi);

	    break;
	}

	case GT_AnnoEle: {
	    anno_ele_block_t *b = (anno_ele_block_t *)&ci_new->data;
	    anno_ele_t *oe  = (anno_ele_t *)&sub_ci->data, *e;
	    //int sub_rec = oe->rec & (ANNO_ELE_BLOCK_SZ-1);

	    /* Already duplicated? */
	    if (b->ae[oe->idx]) {
		sub_new = sub_ci;
		break;
	    }
	    
	    //printf("Duplicate anno %d in block %d\n", sub_rec, ci->rec);

	    sub_new = (cached_item *)malloc(sizeof(*ci) + sub_ci->data_size);
	    memcpy(sub_new, sub_ci, sizeof(*ci) + sub_ci->data_size);
	    e = (anno_ele_t *)&sub_new->data;

	    e->comment = (char *)&e->data;
	    e->block = b;
	    b->ae[e->idx] = e;

	    /* Bump reference count of master again */
	    HacheTableIncRef(ci_new->hi->h, ci_new->hi);

	    break;
	}

	default:
	    sub_new = NULL;
	    break;
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
    
    if (io->read_only)
	return NULL;

    if (io->base && mi->lock_mode < G_LOCK_RW) {
	ci = cache_dup(io, ci);
	mi = cache_master(ci);
	data = &ci->data;
    }

    /* Ensure it's locked RW */
    if (mi->lock_mode < G_LOCK_RW) {
	if (-1 == cache_upgrade(io, mi, G_LOCK_RW)) {
	    fprintf(stderr, "lock denied for rec %"PRIrec"\n", mi->rec);
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

/*
 * Deallocates a record from disc, freeing up the record number in the process
 * too. For now the record we deallocate shouldn't have been modified (ie
 * made r/w). There is no specific need for this constraint except for the
 * fact that I cannot see why we would modify something and then want to
 * throw away both the changes and completely destroy the original record
 * in the same step. So for now we check as it may spot errors. Relax this
 * if needed.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cache_deallocate(GapIO *io, void *data) {
    cached_item *ci = ci_ptr(data);

    /* Mark it for removal */
    ci->forgetme = 1;

    return 0;
}

/*
 * As above, but destroys a type/record instead of a cache object we
 * have a pointer to.
 */
int cache_rec_deallocate(GapIO *io, int type, tg_rec rec) {
    void *v = cache_search(io, type, rec);
    cached_item *ci;

    if ((ci = ci_ptr(v))) {
	/* Need write access to remove a record */
	if (ci->lock_mode < G_LOCK_RW) {
	    if (-1 == cache_upgrade(io, ci, G_LOCK_RW)) {
		fprintf(stderr, "lock denied for rec %"PRIrec"\n", ci->rec);
		return -1;
	    }
	}

	/* Mark for removal */
	ci->forgetme = 1;

	/* Actual deallocate will be done at next flush */
	HacheTableIncRef(ci->hi->h, ci->hi);

	return 0;
    } else {
	return -1;
    }
}
