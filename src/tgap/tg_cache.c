#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "tg_gio.h"
#include "misc.h"

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
    io->iface->seq.unlock(io->dbh, ci->view);
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

    assert(new->hi->data.p == ci);
    new->data_size = size;
    new->hi->data.p = new;

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

    switch (k->type) {
    case GT_Database:
	ci = io->iface->database.read(io->dbh, k->rec);
	break;
	
    case GT_Seq:
	ci = io->iface->seq.read(io->dbh, k->rec);
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

    default:
	return NULL;
    }

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

    assert(ci->updated == 0);

    switch (ci->type) {
    case GT_Seq:
	seq_unload(io, ci);
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

    if (NULL == (h = HacheTableCreate(32768, HASH_DYNAMIC_SIZE|HASH_OWN_KEYS)))
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

    return (*c1)->view - (*c2)->view;
}

/*
 * Flushes changes in the cache back to disk.
 */
int cache_flush(GapIO *io) {
    int i, ret = 0;
    HacheTable *h = io->cache;
    Array to_flush;
    int nflush = 0;

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

		/* Purge from parent */
		HacheTableRemove(io->base->cache,
				 ci->hi->key, ci->hi->key_len, 0);
		/* FIXME: search for parent ci first and deallocate after */

		/* Move this item to parent */
		data.p = ci;
		ci->hi = HacheTableAdd(io->base->cache, 
				       ci->hi->key, ci->hi->key_len,
				       data, NULL);

		/* Remove from this hache */
		HacheTableRemove(io->cache, ci->hi->key, ci->hi->key_len, 0);
	    }
	}

	return cache_flush(io->base);
    }

    to_flush = ArrayCreate(sizeof(cached_item *), 8192);

    /* Identify the list of items that need flushing */
    //fprintf(stderr, "\n");
    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    cached_item *ci = hi->data.p;
	    //fprintf(stderr, "Item %d, type %d, updated %d\n",
	    //    ci->view, ci->type, ci->updated);
	    if (ci->updated) {
		arr(cached_item *, to_flush, nflush++) = ci;
	    }
	}
    }

    /* Sort them by view number, which is likely to be approx on-disk order */
    qsort(ArrayBase(cached_item *, to_flush), nflush, sizeof(cached_item *),
	  qsort_ci_view);

    /* Flush them out */
    for (i = 0; i < nflush; i++) {
	cached_item *ci = arr(cached_item *, to_flush, i);
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

	case GT_Bin:
	    ret = bin_write(io, ci);
	    break;

	case GT_Track:
	    ret = track_write(io, ci);
	    break;

	case GT_RecArray:
	    ret = array_write(io, ci);
	    break;

	default:
	    fprintf(stderr, "Unable to write object of type %d\n",
		    ci->type);
	    ret = -1;
	}

	if (ret == 0) {
	    ci->updated = 0;
	    HacheTableDecRef(io->cache, ci->hi);
	}
    }

    ArrayDestroy(to_flush);

    io->iface->commit(io->dbh);

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
    cache_key_t k = construct_key(rec, type);
    HacheItem *hi = HacheTableQuery(io->cache, (char *)&k, sizeof(k));

    /* Pass one layer up if we're an overlay on top of another GapIO */
    if (!hi && io->base) {
	return cache_search(io->base, type, rec);
    } else if (!hi) {
	/* Otherwise if it's not found, force a load */
	hi = HacheTableSearch(io->cache, (char *)&k, sizeof(k));
    }

    if (!hi)
	return NULL;

    return &((cached_item *)hi->data.p)->data;
}

void cache_incr(GapIO *io, void *data) {
    cached_item *ci = ci_ptr(data);
    HacheTableIncRef(io->cache, ci->hi);
}

void cache_decr(GapIO *io, void *data) {
    cached_item *ci = ci_ptr(data);
    HacheTableDecRef(io->cache, ci->hi);
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
cached_item *cache_dup(GapIO *io, cached_item *ci) {
    HacheItem *hi_old = ci->hi;
    HacheItem *hi_new;

    /* Force load into this io, if not loaded already */
    hi_new = HacheTableQuery(io->cache, hi_old->key, hi_old->key_len);
    if (!hi_new) {
	cached_item *ci_new;
	HacheData hd;

	/* Duplicate the cached_item into this io, keeping the same view */
	ci_new = (cached_item *)malloc(sizeof(*ci) + ci->data_size);
	memcpy(ci_new, ci, sizeof(*ci) + ci->data_size);
	hd.p = ci_new;
	hi_new = HacheTableAdd(io->cache, hi_old->key, hi_old->key_len,
			       hd, NULL);
	ci_new->hi = hi_new;

	switch(ci_new->type) {
	case GT_Seq: {
	    /* reset internal name/seq/conf pointers */
	    seq_t *s = (seq_t *)&ci_new->data;
	    s->name = (char *)&s->data;
	    s->seq = s->name + s->name_len + 1;
	    s->conf = s->seq + ABS(s->len);
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
	}
	}
    }

    if (!hi_new)
	return NULL;

    return (cached_item *)hi_new->data.p;
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

    if (io->base) {
	ci = cache_dup(io, ci);
	data = &ci->data;
    }

    /* Ensure it's locked RW */
    if (ci->lock_mode < G_LOCK_RW) {
	if (-1 == cache_upgrade(io, ci, G_LOCK_RW)) {
	    fprintf(stderr, "lock denied for rec %d\n", ci->rec);
	    return NULL;
	}
    }

    /* Also set updated flag to 1 and bump ref count if needed */
    if (!ci->updated) {
	ci->updated = 1;
	HacheTableIncRef(io->cache, ci->hi);
    }

    return data;
}
