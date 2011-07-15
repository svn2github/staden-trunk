#ifndef _TG_CACHE_ITEM_H_
#define _TG_CACHE_ITEM_H_

#include "hache_table.h"

/*
 * A cached item is some ancillary house-keeping components plus an
 * in-memory representation of the object itself. When querying the cache
 * we return pointer to the object, but this isn't a legal address for
 * free(). That works in our favour infact as it can be trapped by Valgrind
 * and we only intend for allocation and deallocation to occur inside
 * this file.
 *
 * Note that as a confusion the data component is not a real pointer. Instead
 * it's a placeholder for the data itself - ie the object is at
 * &cached_item.data and not pointed to by cached_item.data. This removes
 * the need for an additional malloc and also serves as a check against
 * other code freeing our data structures (as the object we actually return
 * isn't at the start of a malloced block so valgrind etc would complain).
 *
 * The HacheItem pointer is here so we can realloc this object if we need
 * to grow it and update the HacheTable pointer to the new object.
 */
typedef struct cached_item_s {
    GView view;
    int8_t lock_mode;
    uint8_t updated;
    uint8_t forgetme;
    int8_t type;
    tg_rec rec;
    HacheItem *hi;
    size_t data_size;
    int chk_sum;
    void *data;
} cached_item;

/*
 * A hideous hack to go from a cached data item to the cached_item
 * itself.
 */
#define ci_ptr(c) \
    ((cached_item *)((char *)(c) - offsetof(cached_item, data)))

cached_item *cache_new(int type, tg_rec rec, GView v,
		       HacheItem *hi, size_t e_len);

#endif /* _TG_CACHE_ITEM_H_ */
