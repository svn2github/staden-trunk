#ifndef _G_ALLOC_H_
#define _G_ALLOC_H_

#include <inttypes.h>

#define NPOOLS 155

/* Our disk-heap struct */
typedef struct {
    int fd;                     /* File descriptor to heap on disk */
    uint64_t pool[NPOOLS];      /* A pointer to a free node, 0 if none */
    uint64_t maxsz[NPOOLS];	/* Maximum size of data in this pool */
    int next_free_pool[NPOOLS]; /* A skip to the next pool with data */
    int next_free_time[NPOOLS]; /* next_free_pool[x]==timer => cache valid */
    int timer;			/* ++ when pool[] changes; invalidate cache */
    uint64_t wilderness;        /* End of file marker */
} dheap_t;


/*
 * Load the heap meta-data from an existing heap-file.
 * Returns dheap_t pointer on success.
 *         NULL on failure
 */
dheap_t *heap_fdload(int fd);
dheap_t *heap_load(char *file, int mode);

/*
 * Create a new heap in 'file'.
 * Returns dheap_t pointer on success.
 *         NULL on failure
 */
dheap_t *heap_create(char *file);

/*
 * Deallocates memory used by the heap, optionally closing the internal
 * filedescriptor too if close_me is true.
 */
void heap_destroy(dheap_t *h, int close_me);

/*
 * Allocates an object of a specific length from the disk heap 'h'.
 * Returns the file offset on success and size actually allocated.
 *        -1 on failure
 */
int64_t heap_allocate(dheap_t *h, uint32_t length, uint32_t *alloc);

/*
 * Frees the object at 'pos' from the disk heap h.
 * Returns 0 for success
 *        -1 for failure (shouldn't happen!)
 */
int heap_free(dheap_t *h, int64_t pos);

/*
 * Brute force check on the heap validity.
 * It's pretty inefficient and could be sped up, but for now this is just
 * a debugging tool rather than expecting to use it as a full integrity
 * checker.
 */
void heap_check(dheap_t *h);

#endif /* _G_ALLOC_H_ */
