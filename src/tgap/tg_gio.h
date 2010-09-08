#ifndef _TG_GIO_H_
#define _TG_GIO_H_

#include "tg_struct.h"
#include "tg_iface.h"
#include "hache_table.h"

/* ------------------------------------------------------------------------- */
/* The GapIO structure itself - the starting point of all calls */

/*
 * This implements the GapIO structure and ancillary functions
 * The structure itself is the primary I/O mechanism for manipulating Gap5
 * objects.
 *
 * It contains a default database interface pointer too, but this can
 * optionally be passed over as non-null to higher level functions
 * that want to support multiple interfaces (such as the consensus
 * algorithm).
 *
 * The 'base' pointer allows for GapIO overlays. Eg the contig editor
 * can have its own cache, but in reality it simply forwards the vast
 * majority of queries through to the base cache instead.
 */

#define COMP_MODE_ZLIB 0
#define COMP_MODE_NONE 1
#define COMP_MODE_LZMA 2

/* GapIO itself */
typedef struct GapIO {
    /* --- Valid always --- */
    /* DB object cache */
    HacheTable *cache;
    struct GapIO *base;

    /* --- Valid if base == NULL --- */
    /* I/O interface */
    iface *iface;
    void *dbh; /* Database handle to pass into iface functions */

    /* Cached components of the database */
    database_t *db; /* Cached database and view */

    /* Contig order array, also maps contig number to rec.num */
    ArrayStruct *contig_order;
    ArrayStruct *library;

    /* View information - FIXME: move elsewhere*/
    int contig_num;

    /* Contig registration scheme hooks */
    HacheTable *contig_reg;     /* Registration arrays for each contig */
    HacheTable *contig_cursor;	/* Hash of cursor_t lists */

    /* Minimum size for newly created bins */
    int min_bin_size;

    /* Whether we opened in read-write or read-only */
    int read_only;

    /* DB name */
    char *name;
} GapIO;

GapIO *gio_open(char *fn, int ro, int create);
void gio_close(GapIO *io);
GapIO *gio_child(GapIO *io_p);
int gio_read_contig(GapIO *io, int cnum, contig_t **c);


/* ------------------------------------------------------------------------- */
/* The cache aspect of the I/O mechanism */
int cache_create(GapIO *io);
void cache_destroy(GapIO *io);
int cache_flush(GapIO *io);
int cache_updated(GapIO *io);
void *cache_search(GapIO *io, int type, tg_rec rec);
void *cache_search_no_load(GapIO *io, int type, tg_rec rec);
int cache_upgrade(GapIO *io, cached_item *ci, int mode);
void *cache_lock(GapIO *io, int type, tg_rec rec, int mode);
void *cache_item_resize(void *item, size_t size);

void *cache_rw(GapIO *io, void *data);

/* New preferred cache incr/decr functions */
void cache_incr(GapIO *io, void *data);
void cache_decr(GapIO *io, void *data);

/* Create and initialise a new item. */
tg_rec cache_item_create(GapIO *io, int type, void *from);
int cache_item_init(GapIO *io, int type, void *from, tg_rec rec);
int cache_item_remove(GapIO *io, int type, tg_rec rec);
int cache_rec_deallocate(GapIO *io, int type, tg_rec rec);

/* ------------------------------------------------------------------------- */
/* And now the object specific defintions */
#include "tg_register.h"
#include "tg_contig.h"
#include "tg_sequence.h"
#include "tg_bin.h"
#include "tg_track.h"
#include "tg_anno.h"
#include "tg_library.h"
#include "tg_utils.h"
#include "tg_tcl.h"

#endif /* _TG_GIO_H_ */
