/*
 * File: tg_register.c: derived from gap4/io-reg.c
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *         Hills Road
 *         Cambridge CB2 2QH
 *         United Kingdom
 *
 * Description: Handles registration of gap data with display routines.
 *
 * Created: 9 March 1995
 */

/*
 * This allows arbitrary code to register for contig-level events.
 *
 * We implement this by storing a function pointer and void* data in a
 * hash table keyed on the contig number. The hash table is setup to allow
 * duplicate entries so we can have multiple displays registered for
 * a single contig. (The use of a hash table also acts as a sparse array, so
 * the 'key' here is the real contig record number and not an index
 * into a 1..NContigs array.)
 *
 * We also register data by their given ID in the same hache table.
 * This is achieved by negating the ID and ensuring it starts at 1 (ie -1).
 * Hence we can quickly obtain the function/data pairing from a given ID
 * or alternatively we can search through a contig instead.
 */

#include <string.h>
#include <time.h>
#include <assert.h>

#include "tg_gio.h"
#include "hache_table.h"
#include "xalloc.h"

/*
 *-----------------------------------------------------------------------------
 * (De)registration functions
 *-----------------------------------------------------------------------------
 */

/* Debugging aid */
//#define LOG_FILE
static void log_file(FILE *fp, char *buf) {
    puts(buf);
}

/*
 * Initialise the contig register lists
 *
 * Returns 0 on success and -1 for error;
 */
int contig_register_init(GapIO *io) {
    io_contig_reg(io) = NULL;
    io_cursor_reg(io) = NULL;

    /* Create blank hash tables */
    if (NULL == (io_contig_reg(io) = HacheTableCreate(16384,
						      HASH_DYNAMIC_SIZE |
						      HASH_OWN_KEYS |
						      HASH_ALLOW_DUP_KEYS)))
	return -1;
    if (NULL == (io_cursor_reg(io) = HacheTableCreate(16384,
						      HASH_DYNAMIC_SIZE |
						      HASH_OWN_KEYS)))

	return -1;

    io_contig_reg(io)->name = "io_contig_reg(io)";
    io_cursor_reg(io)->name = "io_cursor_reg(io)";

    io_contig_reg(io)->load = NULL;
    io_contig_reg(io)->del  = NULL;
    io_cursor_reg(io)->load = NULL;
    io_cursor_reg(io)->del  = NULL;

    return 0;
}

/*
 * Deallocates memory used by the contig registration scheme.
 */
void contig_register_destroy(GapIO *io) {
    if (io_contig_reg(io))
	HacheTableDestroy(io_contig_reg(io), 0);
    if (io_cursor_reg(io))
	HacheTableDestroy(io_cursor_reg(io), 0);

    io_contig_reg(io) = NULL;
    io_cursor_reg(io) = NULL;
}

/*
 * Debugging output
 */
static void contig_register_dump(GapIO *io) {
    HacheTable *h = io_contig_reg(io);
    HacheIter *iter;
    HacheItem *hi;

    puts("====contig_register_dump====");
    iter = HacheTableIterCreate();
    while (hi = HacheTableIterNext(h, iter)) {
	contig_reg_t *r = (contig_reg_t *)hi->data.p;
	printf("%p(%p) %p %2d %d*fn %p/%p {%p,%p}\n",
	       hi, hi->h, r,
	       *(int *)hi->key,
	       r->ref_count,
	       r->func, r->fdata,
	       r->hi[0], r->hi[1]);

	assert(h == hi->h);
	assert(hi == r->hi[0] || hi == r->hi[1]);
    }
    puts("");
}

/*
 * Allocates new register id numbers as requires. No protection for wrap-
 * around, but this ought to take ages to ever happen.
 */
int register_id() {
    static int id = 0;

    return ++id;
}


/*
 * An internal function to handle removal of the contig_reg_t from our
 * hache tables. This code is shared by multiple places; basically whenever
 * the reference count hits zero.
 *
 * Iter and next maybe NULL, but if non-null they are used to make
 * sure that the current next element in the iterator is still
 * valid. As we are removing data, we may need to step on to the next
 * iteration if we're invalidating our current next candidate.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int contig_reg_remove(GapIO *io, contig_reg_t *r,
			     HacheIter *iter, HacheItem **next) {
    int i;

    if (!r)
	return -1;

    /* Mark for removal */
    r->flags |= REG_FLAG_INACTIVE;

    /* But possibly still in use, so delay */
    if (!r->ref_count == 0)
	return 0;

    /* Yes - so remove from both locations in hache */
    for (i = 0; i < 2; i++) {
	if (!r->hi[i])
	    continue;

	if (next && *next == r->hi[i]) {
	    if (iter)
		*next = HacheTableIterNext(r->hi[i]->h, iter);
	    else
		*next = NULL;
	}

	if (HacheTableDel(io_contig_reg(io), r->hi[i], 0) != 0)
	    return -1;
    }

    free(r);

    return 0;
}


/*
 * As per send_event below, but we broadcast to every contig and to contig
 * 0 too.
 */
void broadcast_event(GapIO *io, HacheTable *h,
		     reg_data *msg, int except) {
    HacheIter *iter;
    int bitcheck = msg->job;
    HacheItem *hi, *next;

    /* Incr ref counts to ensure they're not removed, yet */
    iter = HacheTableIterCreate();
    while (hi = HacheTableIterNext(h, iter)) {
	contig_reg_t *cr = (contig_reg_t *)hi->data.p;
	cr->ref_count++;
    }

    /* Do the notification */
    HacheTableIterReset(iter);
    while (hi = HacheTableIterNext(h, iter)) {
	contig_reg_t *cr = (contig_reg_t *)hi->data.p;
	int key = *(int *)hi->key;

	/*
	 * One item may register to multiple contigs, so iterate through
	 * registration IDs only - aka -ve indices in the hash.
	 */
	if (key >= 0)
	    continue;

	if (cr->flags & REG_FLAG_INACTIVE)
	    continue;

	if (cr->flags & bitcheck && cr->id != except)
	    cr->func(io, 0, cr->fdata, msg);
    }

    /* Decr ref counts, removing if appropriate */
    HacheTableIterReset(iter);
    for (hi = HacheTableIterNext(h, iter); hi; hi = next) {
	contig_reg_t *cr = (contig_reg_t *)hi->data.p;
	int key = *(int *)hi->key;
	next = HacheTableIterNext(h, iter);

	/* Skip 'contig' root list entries */
	if (key >= 0)
	    continue;

	if (--cr->ref_count == 0) {
	    puts("delete me");
	    contig_reg_remove(io, cr, iter, &next);
	}
    }
}

/*
 * Sends message 'msg' to 'contig' if they've registered for message
 * types msg->job. A single ID can be excepted from getting an event.
 * This is typically the window that's generating the event in the first
 * place. Set except to -1 to send to all registered functions.
 *
 * If contig is negative it implies send to item registered with contig 0,
 * but claiming to be for contig '-contig' instead.
 *
 * We need to be particularly careful here about what happens if the
 * callback function triggers new items to be registered or removed.
 */
static void send_event(GapIO *io, HacheTable *h, int contig,
		       reg_data *msg, int except) {
    HacheItem *hi;
    int bitcheck = msg->job;
    int ocontig = contig;

    if (!h)
	return;

    if (contig < 0) {
	ocontig = -contig;
	contig = 0;
    }

    /* Incr ref counts to ensure they're not removed */
    hi = HacheTableSearch(h, (char *)&contig, sizeof(contig));
    while (hi) {
	contig_reg_t *cr = (contig_reg_t *)hi->data.p;
	cr->ref_count++;
	hi = HacheTableNext(hi, (char *)&contig, sizeof(contig));
    }

    /* And now do the actual notification */
    hi = HacheTableSearch(h, (char *)&contig, sizeof(contig));
    while (hi) {
	contig_reg_t *cr = (contig_reg_t *)hi->data.p;
	int key = *(int *)hi->key;
	hi = HacheTableNext(hi, (char *)&contig, sizeof(contig));

	/* Skip 'contig' root list entries */
	if (contig < 0 && key >= 0)
	    continue;

	if (cr->flags & REG_FLAG_INACTIVE)
	    continue;

	if (cr->flags & bitcheck && cr->id != except)
	    cr->func(io, ocontig, cr->fdata, msg);
    }

    /* Decr ref counts and remove if appropriate */
    hi = HacheTableSearch(h, (char *)&contig, sizeof(contig));
    while (hi) {
	HacheItem *next;
	contig_reg_t *cr = (contig_reg_t *)hi->data.p;

	next = HacheTableNext(hi, (char *)&contig, sizeof(contig));

	if (--cr->ref_count == 0) {
	    puts("delete me");
	    contig_reg_remove(io, cr, NULL, &next);
	}

	hi = next;
    }
}


/*
 * Registers func(io, contig, fdata, jdata) with contig 'contig'.
 * Doesn't check if the (func,fdata) pair are already existant.
 *
 * NB. DONT CALL THIS UNTIL YOU'VE INITIALISED ALL DATA. An example is
 * with the template display and the editor. Contig_register sends a
 * REGISTER notification to the editor, which sends out a CURSOR_NOTIFY,
 * which gets back to the routine calling contig_register so that it
 * can draw the editor cursor. So make sure the calling routine has
 * finalised all updates before registering as notifications will
 * start arriving IMMEDIATELY.
 *
 * Returns 0 on success, and -1 for error. 
 */
int contig_register(GapIO *io, int contig,
		    void (*func)(GapIO *io, int contig, void *fdata,
				 reg_data *jdata),
		    void *fdata,
		    int id, int flags, int type) {
    contig_reg_t *r;
    reg_register reg;
    static int uid = 0;
    HacheData hd;
    HacheItem *hi;
    int nid;

    /* Allocate new item */
    if (NULL == (r = (contig_reg_t *)calloc(1, sizeof(*r))))
	return -1;
    hd.p = r;
    hi = HacheTableAdd(io_contig_reg(io), (char *)&contig, sizeof(contig),
		       hd, NULL);
    HacheTableIncRef(io_contig_reg(io), hi);
    r->hi[0] = hi;

    nid = -id;
    hi = HacheTableAdd(io_contig_reg(io), (char *)&nid, sizeof(nid),
		       hd, NULL);
    HacheTableIncRef(io_contig_reg(io), hi);
    r->hi[1] = hi;

#ifdef LOG_FILE    
  {
    static int last_id = -1;
    /* Log file */
    if (id != last_id || 1) {
	char buf[1024], buf2[1024];
	reg_query_name qn;

	qn.job = REG_QUERY_NAME;
	qn.line = buf;
	buf[0] = 0;
	func(io, contig, fdata, (reg_data *)&qn);
	sprintf(buf2, "> Register id=%d cnum=%d func=%p data=%p :%.900s",
		id, contig, func, fdata, buf);
	log_file(NULL, buf2);
	last_id = id;
    }
  }
#endif

    /* Copy over our data */
    r->func = func;
    r->fdata = fdata;
    r->id = id;
    r->time = time(NULL);
    r->flags = flags;
    r->type = type;
    r->uid = ++uid;
    r->ref_count = 1;

    //contig_register_dump(io);

    /*
     * Notify registration to those that are interested.
     * Send notifications to both contig 'n' and contig 0.
     */
    reg.job = REG_REGISTER;
    reg.contig = contig;
    reg.id = id;
    reg.type = type;
    send_event(io, io_contig_reg(io),  contig, (reg_data *)&reg, -1);
    send_event(io, io_contig_reg(io), -contig, (reg_data *)&reg, -1);
    return 0;
}


/*
 * Deregisters func(io, contig, data, jdata) from contig 'contig'.
 * Contig 0 represents the global state which gets all messages.
 * Contig maybe negative in which case it indicates a specific result id,
 * but we still need to find the positive version too.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int contig_deregister(GapIO *io, int contig,
		      void (*func)(GapIO *io, int contig, void *fdata,
				   reg_data *jdata),
		      void *fdata) {
    contig_reg_t *r = NULL;
    reg_register reg;
    HacheItem *hi, *next;
    HacheIter *iter;

    //contig_register_dump(io);

    /*
     * NB: we may have multiple registrations for this func/fdata pair.
     * We may also have been given a registration ID rather than a contig
     * number. For now we take the brute force approach of just searching
     * the entire hache table. It's not ideal, but is still bounded by
     * the number of plots (mostly).
     */
    iter = HacheTableIterCreate();

    next = HacheTableIterNext(io_contig_reg(io), iter);
    while (hi = next) {
	next = HacheTableIterNext(io_contig_reg(io), iter);
	r = (contig_reg_t *)hi->data.p;
	if (r->func != func ||
	    r->fdata != fdata)
	    continue;

	/*
	 * Found a candidate.
	 * We decrement the reference count if active and mark it as inactive.
	 *
	 * If this makes the reference count hit zero then we remove it,
	 * otherwise we leave it and hope it's removed later once the reference
	 * count hits zero.
	 */
	if (!(r->flags & REG_FLAG_INACTIVE)) {
	    r->flags |= REG_FLAG_INACTIVE;
	}

	/* Notify the deregistration to those that are interested */
	reg.job = REG_DEREGISTER;
	reg.contig = contig;
	reg.id = r->id;
	reg.type = r->type;

	/*
	 * We may need to delay these until after removal incase we generate
	 * more removals.
	 */
	send_event(io, io_contig_reg(io),  contig, (reg_data *)&reg, -1);
	send_event(io, io_contig_reg(io), -contig, (reg_data *)&reg, -1);
   
	if (--r->ref_count == 0)
	    contig_reg_remove(io, r, iter, &next);
    }

    HacheTableIterDestroy(iter);

    //contig_register_dump(io);

    return 0;
}

/*
 * Event handling for destroying a contig.
 * Call this just after we delete it. It sends out notification events and
 * updated the contig_reg hash tables.
 */
void contig_register_delete(GapIO *io, int contig) {
    HacheTable *h = io_contig_reg(io);
    HacheItem *hi;
    reg_delete rd;

    while (io->base)
	io = io->base;

    /* Send the event */
    rd.job = REG_DELETE;
    contig_notify(io, contig, (reg_data *)&rd);

    /*
     * Remove items from the hash table incase the objects registered with
     * the contig haven't done this themselves.
     */
    hi = HacheTableSearch(h, (char *)&contig, sizeof(contig));
    while (hi) {
	contig_reg_t *cr = (contig_reg_t *)hi->data.p;
	hi = HacheTableNext(hi, (char *)&contig, sizeof(contig));

	if (--cr->ref_count == 0)
	    contig_reg_remove(io, cr, NULL, NULL);
    }
}

/*
 *-----------------------------------------------------------------------------
 * Cursor functions
 *-----------------------------------------------------------------------------
 */
static int cursor_id = 0;

/*
 * Returns the head of the cursor list for 'contig'
 *         NULL if none found.
 */
static cursor_t *io_cursor_get(GapIO *io, int contig) {
    HacheItem *hi;
    HacheTable *h = io_cursor_reg(io);

    if (!h)
	return NULL;

    hi = HacheTableSearch(io_cursor_reg(io), (char *)&contig, sizeof(contig));
    if (!hi)
	return NULL;

    return (cursor_t *)hi->data.p;
}

/* Sets the cursor list to 'gc' */
static void io_cursor_set(GapIO *io, int contig, cursor_t *gc) {
    
    /* Remove old incase gc == NULL or we have an old anyway */
    HacheTableRemove(io_cursor_reg(io), (char *)&contig, sizeof(contig), 0);

    /* Set the new value */
    if (gc) {
	HacheData hd;
	hd.p = gc;
	HacheTableAdd(io_cursor_reg(io), (char *)&contig, sizeof(contig),
		      hd, NULL);
    }
    return;
}

/*
 * Create a cursor for this contig.
 * If private == 1 then create a new cursor if no non-private ones can be
 * found. Otherwise use an existing one if available.
 * Ie, private cursors can be used more than once, but only by one private
 * user at any one time.
 * Returns the cursor pointer, or NULL for failure.
 */
cursor_t *create_contig_cursor(GapIO *io, int contig, int private, int sent_by)
{
    cursor_t *gc, *gcend;
    reg_cursor_notify cn;

    /* If private, look for a non private cursor to use */
    if (private) {
	for (gc = io_cursor_get(io, contig); gc; gc = gc->next)
	    if (!gc->private)
		break;
	if (gc) {
	    gc->private = private;
	    gc->refs++;
	    goto notify;
	}
    }

    /* If not private, use the first cursor we find */
    if (!private && (gc = io_cursor_get(io, contig))) {
	gc->refs++;
	goto notify;
    }

    /* Otherwise, create our own */
    if (NULL == (gc = xmalloc(sizeof(*gc))))
	return NULL;
    gc->id = cursor_id++;
    gc->refs = 1;
    gc->seq = 0;
    gc->pos = 1;
    gc->abspos = 1;
    gc->private = private;
    gc->sent_by = sent_by;
    gc->next = 0;

    /* Find last existing cursor */
    gcend = io_cursor_get(io, contig);
    while(gcend && gcend->next)
	gcend = gcend->next;

    /* Add gc to it */
    if (gcend)
	gcend->next = gc;
    else
	io_cursor_set(io, contig, gc);

 notify:
    cn.cursor = gc;
    cn.job = REG_CURSOR_NOTIFY;
    gc->job = CURSOR_MOVE | CURSOR_INCREMENT;
    contig_notify(io, contig, (reg_data *)&cn);
#ifdef DEBUG
    printf("CREATE_CURSOR contig %d id %d refs %d private %d\n", 
	   contig, gc->id, gc->refs, gc->private);
#endif

    return gc;
}

/*
 * Given a cursor identifier, return the cursor structure or NULL if not
 * found. We only look in the specific contig number.
 */
cursor_t *find_contig_cursor(GapIO *io, int contig, int id)
{
    cursor_t *gc;

    for (gc = io_cursor_get(io, contig); gc && gc->id != id; gc = gc->next);

    return (gc && gc->id == id) ? gc : NULL;
}

/*
 * Deletes a contig cursor for this option. If the cursor is in use
 * more than once then this simply decrements the reference count.
 * 'private' indicates whether this cursor was "your private one".
 *
 * Always sends a notification about this cursor, either to say it is
 * being destroyed (abspos = -1), or to say that the reference count
 * (and possible private status) have changed.
 */
void delete_contig_cursor(GapIO *io, int contig, int id, int private)
{
    cursor_t *gc, *gcp;
    int c = contig;
    reg_cursor_notify cn;

    if (!(gc = find_contig_cursor(io, c, id)))
	return;

    if (private)
	gc->private = 0;

    /* Decr ref count */
    gc->job = CURSOR_DECREMENT;
    if (--gc->refs <= 0) {
	gc->job |= CURSOR_DELETE;
    }

    /* Notify of the impending demise or of new reference count */
    cn.cursor = gc;
    cn.job = REG_CURSOR_NOTIFY;
    contig_notify(io, c, (reg_data *)&cn);

    if (gc->refs > 0)
	return;

    /* If first in contig, skip to next */
    if (io_cursor_get(io, c) == gc) {
	io_cursor_set(io, c, gc->next);
	xfree(gc);
	return;
    }

    /* Otherwise, find previous and link to next*/
    for (gcp = io_cursor_get(io, c); gcp && gcp->next != gc; gcp = gcp->next);
    if (gcp && gcp->next == gc) {
	gcp->next = gc->next;
	xfree(gc);
    }

    return;
}


/*
 *-----------------------------------------------------------------------------
 * Contig level actions
 *-----------------------------------------------------------------------------
 */

/*
 * Uses the register list for a given contig to call a particular job.
 * Contig 0 is a special registration list for windows that want to track
 * all contigs, so we always duplicate data there too.
 */
void contig_notify(GapIO *io, int contig, reg_data *jdata) {
    while (io->base)
	io = io->base;

    send_event(io, io_contig_reg(io), contig, jdata, -1);
    if (contig)
	send_event(io, io_contig_reg(io), -contig, jdata, -1);
}

void contig_notify_except(GapIO *io, int contig, reg_data *jdata, int id) {
    while (io->base)
	io = io->base;

    send_event(io, io_contig_reg(io), contig, jdata, id);
    if (contig)
	send_event(io, io_contig_reg(io), -contig, jdata, id);
}

/*
 * Joins two registers lists. This checks for duplicate joins and so won't
 * register something again if it was already registered to cto.
 * Objects previously registered with cfrom and now being moved to cto
 * will be registered generating a new REG_REGISTER event (for whoever is
 * listening). Objects that were already registered in cto are not registered.
 * Finally we drop objects from the cfrom list *without* sending
 * corresponding REG_DEREGISTER events.
 *
 * This also merges the cursor lists.
 *
 * Returns 0 for success, and -1 for error.
 */
int contig_register_join(GapIO *io, int cfrom, int cto) {
    HacheTable *h = io_contig_reg(io);
    HacheItem *hif, *hit;
    contig_reg_t *rf, *rt;
    cursor_t *gc;
    int offset = 0; /* HACK for now */
    
    /* Copy lists */
    hif = HacheTableSearch(h, (char *)&cfrom, sizeof(cfrom));
    while (hif) {
	HacheItem *next;

	rf = (contig_reg_t *)hif->data.p;

	/* Inefficient, but check if already in cto */
	hit = HacheTableSearch(h, (char *)&cto, sizeof(cto));
	while (hit) {
	    rt = (contig_reg_t *)hit->data.p;
	    if (rf->id == rt->id)
		break;
	    hit = HacheTableNext(hit, (char *)&cfrom, sizeof(cto));
	}

	next = HacheTableNext(hif, (char *)&cfrom, sizeof(cfrom));

	if (!hit) {
	    /* Id not found, so move hif to new contig */
	    if (!HacheTableRehash(h, hif, (char *)&cto, sizeof(cto))) {
		fprintf(stderr, "Failed to rehash hi=%p\n", hif);
	    }
	} else {
	    /* Dup, so silently deregister from the cfrom list */
	    HacheTableDel(h, hif, 0);
	}

	hif = next;
    }

    /*
     * Update cursor lists.
     * This basically just concatenates the two cursors together.
     * If the sequence number for a cursor in 'cfrom' is 0 (consensus),
     * we also update the position.
     */
    for (gc = io_cursor_get(io, cto); gc && gc->next; gc = gc->next)
	;
    if (gc) {
	gc->next = io_cursor_get(io, cfrom);
    } else {
	io_cursor_set(io, cto, io_cursor_get(io, cfrom));
    }
    for (gc = io_cursor_get(io, cfrom); gc; gc = gc->next) {
	if (gc->seq == cfrom || gc->seq == cto || gc->seq == 0) {
	    gc->pos += offset;
	    gc->abspos = gc->pos;
	} else {
	    int cnum, pos;
	    sequence_get_position(io, gc->seq, &cnum, &pos, NULL, NULL);
	    gc->abspos = pos + gc->pos;
	}
    }
	
    io_cursor_set(io, cfrom, NULL);

#ifdef LOG_FILE    
    {
	char buf[1024];
	sprintf(buf, "> Register_join done");
	log_file(NULL, buf);
    }
#endif

    return 0;
}

/*
 * Iterates through all contigs_reg_t registered with 'contig' and having
 * a id of 'id'. start_from maybe NULL, but if not then it's used to
 * store the last hache table position to facilitate iteration.
 *
 * An 'id' of zero implies no filtering on ID is needed - return all.
 *
 * Returns contig_reg_t* on success
 *         NULL on failure
 */
contig_reg_t *get_reg_by_contig_id(GapIO *io, int contig, int id,
				   HacheItem **start_from) {
    HacheItem *hi;
    HacheTable *h;
    contig_reg_t *r;

    h = io_contig_reg(io);
    hi = start_from && *start_from
	? HacheTableNext(*start_from, (char *)&contig, sizeof(contig))
	: HacheTableSearch(h, (char *)&contig, sizeof(contig));
    while (hi) {
	r = (contig_reg_t *)hi->data.p;
	if (!id || r->id == id) {
	    if (start_from)
		*start_from = hi;
	    return r;
	}
	hi = HacheTableNext(hi, (char *)&contig, sizeof(contig));
    }

    if (start_from)
	*start_from = NULL;

    return NULL;
}

contig_reg_t *get_reg_by_id(GapIO *io, int id,
			    HacheItem **start_from) {
    HacheItem *hi;
    HacheTable *h;
    int nid = -id;

    h = io_contig_reg(io);
    hi = start_from && *start_from
	? HacheTableNext(*start_from, (char *)&nid, sizeof(nid))
	: HacheTableSearch(h, (char *)&nid, sizeof(id));

    if (start_from)
	*start_from = hi;

    return hi ? (contig_reg_t *)hi->data.p : NULL;
}

contig_reg_t **get_reg_by_type(GapIO *io, int type, int *nresult) {
    HacheTable *h;
    contig_reg_t **res = NULL;
    int nres = 0, nalloc = 0, i;

    h = io_contig_reg(io);
    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	contig_reg_t *r;

	/* Loop through all +ve keys => contig registrations */
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    int key = *(int *)hi->key;
	    if (key < 0)
		continue;

	    if (nres >= nalloc) {
		nalloc += 10;
		res = (contig_reg_t **)realloc(res, nalloc*sizeof(*res));
	    }

	    r = (contig_reg_t *)hi->data.p;
	    if (r->type != type)
		continue;

	    res[nres++] = r;
	}	   
    }

    *nresult = nres;
    return res;
}

/*
 * As per get_reg_by_id iterator, but returning as an array instead.
 * The last element of the array will be NULL.
 */
contig_reg_t **result_to_regs(GapIO *io, int id) {
    contig_reg_t **rl, *r;
    int size = 8, count = 0;
    HacheItem *hi = NULL;

    if (NULL == (rl = (contig_reg_t **)xmalloc(size * sizeof(contig_reg_t *))))
        return NULL;

    while ((r = get_reg_by_id(io, id, &hi))) {
	rl[count++] = r;
	if (count >= size-1) {
	    size *= 2;
	    rl = (contig_reg_t **)xrealloc(rl, size * sizeof(contig_reg_t *));
	    if (NULL == rl) {
		return NULL;
	    }
	}
    }

    /* We've always guaranteed one element in rl to be left free. */
    rl[count] = NULL;

    return rl;
}

/*
 * Fills out a result_name_t struct for every result registered with the
 * GapIO, primarily containing the result name but also a reference to how
 * to contact it more directly.
 *
 * On return *nresult holds the size of the returned array.
 * The caller should deallocate this array using free().
 */
result_name_t *result_names(GapIO *io, int *nresult) {
    result_name_t *res = NULL;
    int nres = 0, nalloc = 0;
    HacheTable *h = io_contig_reg(io);
    int i;

    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi;
	reg_query_name qn;
	contig_reg_t *r;

	/* Loop through all +ve keys => contig registrations */
	for (hi = h->bucket[i]; hi; hi = hi->next) {
	    int key = *(int *)hi->key;
	    if (key < 0)
		continue;
	    
	    if (nres >= nalloc) {
		nalloc += 10;
		res = (result_name_t *)realloc(res, nalloc*sizeof(*res));
	    }

	    r = (contig_reg_t *)hi->data.p;
	    if (r->flags & REG_FLAG_INACTIVE)
		continue;

	    qn.job = REG_QUERY_NAME;
	    qn.line = res[nres].name;
	    r->func(io, 0, r->fdata, (reg_data *)&qn);

	    res[nres].id = r->id;
	    res[nres].contig = key;
	    res[nres].r = r;
	    nres++;
	}
    }

    *nresult = nres;
    return res;
}


char *result_time(GapIO *io, int id) {
    contig_reg_t *r = get_reg_by_id(io, id, NULL);
    static char buf[80];

    if (!r)
	return "unknown";

    /* %r doesn't work for windows ! */
    strftime(buf, sizeof(buf)-1, "%a %I:%M:%S %p", localtime(&r->time));

    return buf;
}

/*
 * Returns the data component of a contig_reg_t for a specific id.
 *
 * We return only the first data for id found, or NULL if none found.
 */
void *result_data(GapIO *io, int id) {
    contig_reg_t *r = get_reg_by_id(io, id, NULL);

    return r ? r->fdata : NULL;
}


/*
 * Uses the register list for a given result to call a particular job.
 * As per contig_notify except on a result basis rather than contig basis.
 * 'all' declares whether to send message to all registrations with this id
 * or just the first found.
 */
void result_notify(GapIO *io, int id, reg_data *jdata, int all) {
    contig_reg_t *r;
    HacheItem *hi = NULL;

    if (jdata->job != REG_GENERIC)
	printf("result_notify(id=%d, jdata->job=%d)\n",
	       id, jdata->job);

    while ((r = get_reg_by_id(io, id, &hi))) {
	if ((r->flags & jdata->job) && !(r->flags & REG_FLAG_INACTIVE)) {
	    r->func(io, 0, r->fdata, jdata);
	    if (!all)
		return;
	}
    }
}

/*
 * Notifies all (or the first) registered items of a given type.
 * Works across all contigs.
 *
 * Returns 0 for success, -1 for failure.
 */
int type_notify(GapIO *io, int type, reg_data *jdata) {
    contig_reg_t **res;
    int nres, i, ret = -1, changed;

    do {
	if (NULL == (res = get_reg_by_type(io, type, &nres)))
	    return ret;

	changed = 0;
	for (i = 0; i < nres; i++) {
	    if ((res[i]->flags & jdata->job) &&
		!(res[i]->flags & REG_FLAG_INACTIVE)) {
		res[i]->func(io, 0, res[i]->fdata, jdata);
		/* The callback function may have changed our arrays */
		changed = 1;
		break;
	    }
	}

	ret = 0;
	free(res);
    } while (changed);

    return 0;
}

/*
 * Returns the id for a registered item of a particular type, If contig is
 * non zero, we search only this contig. Otherwise we scan all.
 *
 * We return only the first data for id found, or 0 if none found.
 */
int type_to_result(GapIO *io, int type, int contig) {
    int nres;
    contig_reg_t **res = get_reg_by_type(io, type, &nres);
    int id = -1;

    if (nres) {
	id = res[0]->id;
	free(res);
    }

    return id;
}


/*
 *-----------------------------------------------------------------------------
 * Lock management
 *-----------------------------------------------------------------------------
 */

void busy_dialog(GapIO *io, int contig) {
    char buf[1024];

    sprintf(buf, "tk_messageBox \
			-icon warning \
			-title {Contig is busy} \
			-message {The contig is busy, probably due to an "
	                         "editor in use for this contig. Changes will "
	                         "not be made for this contig.} \
                        -type ok");

    Tcl_Eval(GetInterp(), buf);
}

/*
 * Attempts to lock a contig for exclusive write access.
 * No record of the lock is kept; that's done implicitly by registration.
 *
 * Returns 0 for success and -1 for failure.
 */
int contig_lock_write(GapIO *io, int contig) {
    reg_get_lock lg;
    reg_set_lock ls;

    /*
     * Notify a lock request. If a view is busy and requires exclusive access
     * rights then it will clear lg.lock
     */
    lg.job = REG_GET_LOCK;
    lg.lock = REG_LOCK_WRITE;

    contig_notify(io, contig, (reg_data *)&lg);

    /*
     * If lg.lock remains set then nothing has objected. We notify our
     * decision to go ahead.
     */
    if (lg.lock & REG_LOCK_WRITE) {
	ls.job = REG_SET_LOCK;
	ls.lock = REG_LOCK_WRITE;

	contig_notify(io, contig, (reg_data *)&ls);

	return 0;
    }

    busy_dialog(io, contig);
    return -1;
}
