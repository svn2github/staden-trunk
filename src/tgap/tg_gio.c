#include <assert.h>
#include <string.h>

#include "tg_gio.h"
#include "misc.h"

#include "tg_iface_g.h"

/* ------------------------------------------------------------------------
 * This is the primary IO layer that the rest of Gap5 uses. It consists of a
 * GapIO struct with functions to create, destroy and some basic common
 * manipulations.
 *
 * Internally it uses an implementation agnostic database interface held
 * in the void *dbh element.
 */

/* ------------------------------------------------------------------------ */
/* Open/close/commit */

/*
 * Open a database, optionally in read-only mode and creating if desired too.
 *
 * Returns GapIO pointer to DB on success.
 *         NULL on failure.
 */
GapIO *gio_open(char *fn, int ro, int create) {
    GapIO *io = (GapIO *)calloc(1, sizeof(*io));
    char *cp;

    io->iface = get_iface_g();
    if (create) {
	if (-1 == io->iface->create(fn))
	    return NULL;
    }

    io->min_bin_size = MIN_BIN_SIZE; /* default */

    /* Initialise the cache */
    cache_create(io);

    if (NULL == (io->dbh = io->iface->connect(fn, ro))) {
	if (!ro) {
	    ro = 1;
	    if (NULL == (io->dbh = io->iface->connect(fn, ro)))
		return NULL;
	} else {
	    return NULL;
	}
    }

    io->read_only = ro;

    if (create) {
	io->iface->database.create(io->dbh, NULL);
    }

    /* Cache the GDatabase struct */
    if (NULL == (io->db = cache_search(io, GT_Database, 0)))
	return NULL;
    cache_incr(io, io->db);

    if (io->db->version > DB_VERSION) {
	verror(ERR_WARN, "Open Database",
	       "Database version %d is too new for this version of gap5",
	       io->db->version);
	return NULL;
    }

    /* Load the contigs array */
    io->contig_order = cache_search(io, GT_RecArray, io->db->contig_order);
    cache_incr(io, io->contig_order);

    /* Load the library array */
    io->library = cache_search(io, GT_RecArray, io->db->library);
    cache_incr(io, io->library);

    /* Initialise the contig and cursor registration hashes */
    contig_register_init(io);

    io->iface->setopt(io->dbh, OPT_COMP_MODE, COMP_MODE_ZLIB);

    /* Copy the name */
    if (NULL == (cp = strrchr(fn, '/')))
	cp = fn;
    io->name = strdup(cp);

    return io;
}


/*
 * Closes a database, automatically committing any unsaved changes.
 * Also frees any associated memory.
 */
void gio_close(GapIO *io) {
    if (io->base) {
	cache_destroy(io);
	free(io);
	return;
    }

    cache_decr(io, io->db);
    cache_decr(io, io->contig_order);
    cache_decr(io, io->library);

    cache_flush(io);
    cache_destroy(io);

    contig_register_destroy(io);

    io->iface->commit(io->dbh);
    io->iface->disconnect(io->dbh);

    if (io->name)
	free(io->name);

    free(io);
}


/*
 * Creates a child of a GapIO. The child can basically be considered as a
 * read-through but write-sticky cache. By that I mean that using the child
 * GapIO works as if the parent was used instead, but writes to the
 * child GapIO are not propagated through to the parent until an
 * explicit request is made.
 *
 * The purpose for this mechanism is to allow tools like the contig editor
 * be able to read and write to a local cached copy of the database with
 * the ability to cancel and effectively rollback all local changes made.
 */
GapIO *gio_child(GapIO *io_p) {
    GapIO *io = (GapIO *)calloc(1, sizeof(*io));

    io->iface = get_iface_g();
    cache_create(io);
    
    io->base = io_p;
    io->dbh = io->base->dbh;
    io->read_only = io->base->read_only;
    io->min_bin_size = io->base->min_bin_size;
    return io;
}

/* ------------------------------------------------------------------------ */
/* 'Binning' callback or utilisation functions */


/* ------------------------------------------------------------------------ */
/* Other DB structures */

/*
 * Reads a given GContig number to the GContig pointer.
 * Here cnum is a contig number from 0 to N-1 indicating the record number
 * in the cnum-th element of the contig_order array.
 *
 * Returns 0 on success.
 *        -1 on failure.
 */
int gio_read_contig(GapIO *io, int cnum, contig_t **c) {
    GRec crec;

    if (!io->contig_order)
	return -1;

    crec = arr(GCardinal, io->contig_order, cnum);
    *c = (contig_t *)cache_search(io, GT_Contig, crec);

    cache_decr(io, *c);

    return 0;
}
