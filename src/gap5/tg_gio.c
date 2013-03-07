#if !defined(__MINGW32__) && !defined(_MSC_VER)
/* Uncomment to enable .log files */
/* #   define DO_LOGGING */
#endif
#include <assert.h>
#include <string.h>
#ifdef DO_LOGGING
#   include <unistd.h>     /* For getuid */
#   include <sys/types.h>
#   include <pwd.h>        /* For getpwuid */
#   include "text_output.h"
#endif
#include "tg_gio.h"
#include "misc.h"
#include "actf.h"

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

void xperror_fatal(char *name, char *str) {
    verror(ERR_FATAL, name, str);
}

void xperror_warn(char *name, char *str) {
    verror(ERR_WARN, name, str);
}

static int db_version = DB_VERSION;
int gio_set_db_version(int vers) {
    if (vers < 0)
	return -1;
    if (vers > DB_VERSION)
	return -1;

    db_version = vers;
    return 0;
}

#ifdef DO_LOGGING
/* Start logging */
static void open_log_file(GapIO *io, char *fn) {
    char log_buf[1024];
    char *user;
    char *logfn;
    struct passwd *pw;
    int uid = getuid();
    size_t name_len = strlen(fn) + 5;

    pw = getpwuid(uid);
    user = pw ? pw->pw_name : "unknown";
    snprintf(log_buf, sizeof(log_buf), "opening %s r%c by %s(%d)",
	     fn, io->read_only ? 'o' : 'w', user, uid);
    logfn = malloc(name_len);
    if (NULL != logfn) {
	snprintf(logfn, name_len, "%s%s", fn, ".log");
	log_file(logfn, log_buf);
	free(logfn);
    }
}
#endif

/*
 * Open a database, optionally in read-only mode and creating if desired too.
 *
 * Returns GapIO pointer to DB on success.
 *         NULL on failure.
 */
GapIO *gio_open(char *fn, int ro, int create) {
    GapIO *io = (GapIO *)calloc(1, sizeof(*io));
    char *cp;
    int lock_err;

    /* Check locks */
    lock_err = actf_lock(ro, fn, create);
    if (!create && (lock_err == 3 || lock_err == 5)) {
	vmessage("Opening database in read only mode instead.\n");
	ro = 1;
	lock_err = actf_lock(ro, fn, create);
    }
    if (lock_err != 0) {
	verror(ERR_WARN, "Open Database",
	       "Unable to lock and/or open the database.");
	return NULL;
    }

    io->iface = get_iface_g();
    if (create) {
	if (0 != io->iface->create(fn)) {
	    xperror("In tg_gio.c:gio_open()", xperror_fatal);
	    return NULL;
	}
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
	io->iface->database.create(io->dbh, NULL, db_version);
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

    /* Load the scaffold array */
    if (io->db->scaffold) {
	io->scaffold =
	    cache_search(io, GT_RecArray, io->db->scaffold);
	cache_incr(io, io->scaffold);
    } else {
	/* FIXME: create a dummy order of 1 per contig? */
	io->scaffold = 0;
    }
    
    /* Load the library array */
    io->library = cache_search(io, GT_RecArray, io->db->library);
    cache_incr(io, io->library);

    /* Initialise the contig and cursor registration hashes */
    contig_register_init(io);

    io->iface->setopt(io->dbh, OPT_COMP_MODE, COMP_MODE_ZLIB);

    /* Copy the name */
    if ((cp = strrchr(fn, '/')))
	cp++;
    else
	cp = fn;
    io->name = strdup(cp);

    io->last_bin = 0;
    io->incr_svalue = io->incr_rvalue = io->incr_avalue = 0;

    io->max_template_size = 0;

    io->debug_level = 0;
    io->debug_fp = stderr;

    //update_uniqueness_hash(io);

#ifdef DO_LOGGING
    open_log_file(io, fn);
#endif

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

#ifdef DO_LOGGING
    {
	char buf[256];
	snprintf(buf, sizeof(buf), "closing database %s ...",
		 io->name ? io->name : "");
	log_file(NULL, buf);
    }
#endif

    cache_decr(io, io->db);
    cache_decr(io, io->contig_order);
    if (io->scaffold)
	cache_decr(io, io->scaffold);
    cache_decr(io, io->library);

    cache_flush(io);
    cache_destroy(io);

    contig_register_destroy(io);

    io->iface->commit(io->dbh);
    io->iface->disconnect(io->dbh);

    actf_unlock(io->read_only, io->name);

#ifdef DO_LOGGING
    log_file(NULL, "...closed.");
#endif

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

    assert(0 == io_p->last_bin); /* No pending updates in bin_add_to_range */

    io->iface = get_iface_g();
    cache_create(io);
    
    io->base = io_p;
    io->dbh = io->base->dbh;
    io->read_only = io->base->read_only;
    io->min_bin_size = io->base->min_bin_size;
    io->debug_level = io->base->debug_level;
    io->debug_fp = io->base->debug_fp;
    io->last_bin = 0;
    io->max_template_size = io->base->max_template_size;
    return io;
}

/*
 * Returns the original GapIO regardless of how many children deep we are.
 */
GapIO *gio_base(GapIO *io) {
    while (io->base)
	io = io->base;
    return io;
}

int io_timestamp_incr(GapIO *io) {
    io = gio_base(io);
    io->db = cache_rw(io, io->db);
    return ++io->db->timestamp;
}

/*
 * Sets debugging to a specific level.
 * Level 0 turns off debugging output.
 * Returns the previous debug level.
 */
int gio_debug_level(GapIO *io, int level) {
    int r = io->debug_level;
    io->debug_level = level;

    if (io->iface)
	io->iface->setopt(io->dbh, OPT_DEBUG_LEVEL, level);

    return r;
}

__PRINTF_FORMAT__(3,4)
void gio_debug(GapIO *io, int level, char *fmt, ...) {
    va_list args;

    if (io->debug_level < level)
	return;

    va_start(args, fmt);
    vfprintf(io->debug_fp, fmt, args);
    va_end(args);
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
    tg_rec crec;

    if (io->base)
	io = io->base;

    if (!io->contig_order)
	return -1;

    crec = arr(tg_rec, io->contig_order, cnum);
    *c = (contig_t *)cache_search(io, GT_Contig, crec);

    return 0;
}
