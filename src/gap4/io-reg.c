/*
 * File: io-reg.c
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

#include <string.h>
#include <time.h>

#include "IO.h"
#include "array.h"
#include "io-reg.h"
#include "misc.h"
#include "tk-io-reg.h"
#include "text_output.h"


/*
 *-----------------------------------------------------------------------------
 * (De)registration functions
 *-----------------------------------------------------------------------------
 */

/*
 * Initialise the contig register lists
 *
 * Returns 0 on success and -1 for error;
 */
int contig_register_init(GapIO *io) {
    int i;

    /* Allocate Arrays */
    if (NULL == (io_contig_reg(io) = ArrayCreate(sizeof(Array),
						 Ncontigs(io)+1)))
	return -1;
    if (NULL == (io_cursor_reg(io) = ArrayCreate(sizeof(Array),
						 Ncontigs(io))))
	return -1;

    /* Create empty Arrays for each contig */
    for (i = 0; i <= Ncontigs(io); i++) {
	if (NULL == (io_reg(io, i) = ArrayCreate(sizeof(contig_reg_t), 0)))
	    return -1;
	io_Nreg(io, i) = 0;
	if (i > 0)
	    io_cursor(io, i) = NULL;
    }

    return 0;
}

/*
 * Deallocates memory used by the contig registration scheme.
 */
void contig_register_destroy(GapIO *io) {
    int i;

    for (i = 0; i <= Ncontigs(io); i++)
	ArrayDestroy(io_reg(io, i));

    ArrayDestroy(io_contig_reg(io));
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
    int i, n;
    reg_register reg;
    static int uid = 0;
    static int last_id = -1;

    /* Allocate new item */
    /* ArrayRef already increments io_Nreg, so no need to do that here */
    if (NULL == (r = (contig_reg_t *)ArrayRef(io_reg(io, contig),
					      io_Nreg(io, contig))))
	return -1;

    /* Log file */
    if (id != last_id) {
	char buf[1024];
	sprintf(buf, "> Register id=%d cnum=%d", id, contig);
	log_file(NULL, buf);
	last_id = id;
    }

    /* Copy over our data */
    r->func = func;
    r->fdata = fdata;
    r->id = id;
    r->time = time(NULL);
    r->flags = flags;
    r->type = type;
    r->uid = ++uid;

    /*
     * Notify registration to those that are interested.
     * Send notifications to both contig 'n' and contig 0.
     */
    r = ArrayBase(contig_reg_t, io_reg(io, contig));
    n = io_Nreg(io, contig) - 1;
    reg.job = REG_REGISTER;
    reg.contig = contig;
    reg.id = id;
    reg.type = type;
    for (i=0; i < n; i++) {
	if (r[i].flags & REG_REGISTER) {
	    r[i].func(io, contig, r[i].fdata, (reg_data *)&reg);
	}
    }

    r = ArrayBase(contig_reg_t, io_reg(io, 0));
    n = io_Nreg(io, 0) - 1;
    reg.job = REG_REGISTER;
    reg.contig = contig;
    reg.id = id;
    reg.type = type;
    for (i=0; i < n; i++) {
	if (r[i].flags & REG_REGISTER) {
	    r[i].func(io, contig, r[i].fdata, (reg_data *)&reg);
	}
    }

    /*
     * Result manager - could be done by acknowledging REG_REGISTER requests,
     * but as it needs to monitor for all contigs it's easiest like this.
     */
    update_results(io);

    return 0;
}


/*
 * Deregisters func(io, contig, data, jdata) from contig 'contig'.
 *
 * Returns 0 for success, and -1 for error.
 * If contig==0 then it searches for the (func,jdata) pair in all contigs.
 */
int contig_deregister(GapIO *io, int contig,
		      void (*func)(GapIO *io, int contig, void *fdata,
				   reg_data *jdata),
		      void *fdata) {
    contig_reg_t *r;
    int i, n;
    reg_register reg;

    /* If needed, find contig. Crude, but it works */
    if (contig == 0) {
	int c;

	contig = -1;

	for (c = 0; c <= NumContigs(io) && contig == -1; c++) {
	    for (i = 0; i < io_Nreg(io, c); i++) {
		contig_reg_t *r = &arr(contig_reg_t, io_reg(io, c), i);
		if (r->func == func && r->fdata == fdata) {
		    contig = c;
		    break;
		}
	    }
	}

	if (contig == -1)
	    return -1;
    }

    i = io_Nreg(io, contig);
    r = ArrayBase(contig_reg_t, io_reg(io, contig));
    /* Search for element in the array */
    for (; i > 0; i--) {
	if (r[i-1].func == func && r[i-1].fdata == fdata)
	    break;
    }

    if (i > 0) {
	/* Match found */
	reg_query_name qn;
	char buf[1024], buf2[1024];
	static int last_id = -1;

	/* Find name, for log file */
	if (r[i-1].id != last_id) {
	    qn.job = REG_QUERY_NAME;
	    qn.line = buf;
	    buf[0] = 0;
	    r[i-1].func(io, contig, r[i-1].fdata, (reg_data *)&qn);
	    sprintf(buf2, "> Deregister id=%d cnum=%d: %.900s",
		    r[i-1].id, contig, buf);
	    log_file(NULL, buf2);
	    last_id = r[i-1].id;
	}

	/* Shuffle everything else down */
	reg.job = REG_DEREGISTER;
	reg.contig = contig;
	reg.id = r[i-1].id;
	reg.type = r[i-1].type;
	memmove(&r[i-1], &r[i], (io_Nreg(io, contig) - i) * sizeof(r[i]));

	/* Decrement count */
	io_Nreg(io, contig)--;

	/*
	 * Notify registration to those that are interested
	 * WARNING: This will break if one of these notifications also
	 * causes a deregister.
	 */
	r = ArrayBase(contig_reg_t, io_reg(io, contig));
	n = io_Nreg(io, contig);
	for (i=0; i < n; i++) {
	    if (r[i].flags & REG_DEREGISTER) {
		r[i].func(io, contig, r[i].fdata, (reg_data *)&reg);
	    }
	}

	r = ArrayBase(contig_reg_t, io_reg(io, 0));
	n = io_Nreg(io, 0);
	for (i=0; i < n; i++) {
	    if (r[i].flags & REG_DEREGISTER) {
		r[i].func(io, contig, r[i].fdata, (reg_data *)&reg);
	    }
	}
    }

    /*
     * Result manager - could be done by acknowledging REG_REGISTER requests,
     * but as it needs to monitor for all contigs it's easiest like this.
     */
    update_results(io);

    return 0;
}


/*
 *-----------------------------------------------------------------------------
 * Contig level actions
 *-----------------------------------------------------------------------------
 */

/*
 * Uses the register list for a given contig to call a particular job.
 *
 * Contig 0 is for all contigs. NB: may not work if this causes the number
 *    of contigs to change!
 *
 * except_id is an id value to skip in notifications. Use -1 if you don't
 * have one.
 */
static void contig_notify_common(GapIO *io, int contig, reg_data *jdata,
				 int except_id) {
    contig_reg_t *r;
    int i, j, k, n;
    int *uids;
    signed int send_to_zero;
    int orig_contig;

    send_to_zero = (contig > 0);
    if (contig < 0)
	contig = -contig;
    orig_contig = contig;

    if (contig == 0) {
	reg_buffer_start rs;
	reg_buffer_end re;

	rs.job = REG_BUFFER_START;
	for (i = 1; i <= NumContigs(io); i++)
	    contig_notify_common(io, -i, (reg_data *)&rs, except_id);
	for (i = 1; i <= NumContigs(io); i++)
	    contig_notify_common(io, -i, jdata, except_id);
	re.job = REG_BUFFER_END;
	for (i = 1; i <= NumContigs(io); i++)
	    contig_notify_common(io, -i, (reg_data *)&re, except_id);

	/* Plus flow through for contig 0 */
    }

    /*
     * Notify two lists - one for this contig, and one items registered
     * with contig 0.
     */
    do {
	n = io_Nreg(io, contig);
	r = ArrayBase(contig_reg_t, io_reg(io, contig));

	if (0 == n)
	    return;

	/*
	 * Take a snap shot of all the registrations that need notifying.
	 * We do this by remembering their unique uids.
	 */
	if (NULL == (uids = (int *)xmalloc(n * sizeof(int))))
	    return;

	for (k = 0; k < n; k++)
	    uids[k] = r[k].uid;

	/* Now loop through all uids sending the notification */
	for (j = i = 0; j < k; i++, j++) {
	    /*
	     * Don't optimise this outside the loop - it may not be a constant
	     */
	    n = io_Nreg(io, contig);
	    
	    /*
	     * It's possible that the next registration isn't our next uid.
	     * This occurs when the previous r[i].func() modified the
	     * registration list for this contig (eg by deregistering
	     * something). So we scan along to find the correct 'i' in this
	     * case remembering that it may not even exist any more.
	     */
	    if (i >= n || r[i].uid != uids[j]) {
		for (i = 0; i < n; i++) {
		    if (r[i].uid == uids[j])
			break;
		}
		if (i == n)
		    continue;
	    }

	    /* Send the notification request itself */
	    if (r[i].flags & jdata->job && r[i].id != except_id)
		r[i].func(io, orig_contig, r[i].fdata, jdata);
	}

	xfree(uids);

	contig = 0;
    } while (--send_to_zero >= 0);
}

/*
 * Uses the register list for a given contig to call a particular job.
 *
 * Contig 0 is for all contigs. NB: may not work if this causes the number
 *    of contigs to change!
 */
void contig_notify(GapIO *io, int contig, reg_data *jdata) {
    contig_notify_common(io, contig, jdata, -1);
}


void contig_notify_except(GapIO *io, int contig, reg_data *jdata,
			  int id) {
    contig_notify_common(io, contig, jdata, id);
}


/*
 * Joins two registers lists. This checks for duplicate joins and so won't
 * register something again if it was already registered to cto.
 * This also merges the cursor lists.
 *
 * Returns 0 for success, and -1 for error.
 */
int contig_register_join(GapIO *io, int cfrom, int cto) {
    Array a, ato;
    int n, nto, i, j, found;
    contig_reg_t *r, *rto;
    cursor_t *gc;
    int offset = 0; /* HACK for now */
    
    /* copy list */
    a = io_reg(io, cfrom);
    ato = io_reg(io, cto);
    nto = io_Nreg(io, cto);

    for (i = 0, n = io_Nreg(io, cfrom); i < n; i++) {
	r = &arr(contig_reg_t, a, i);

	/* Check if this registration already exists in cto */
	found = 0;
	for (j = 0; j < nto; j++) {
	    rto = &arr(contig_reg_t, ato, j);
	    if (rto->func == r->func && rto->fdata == r->fdata) {
		found = 1;
		break;
	    }
	}

	if (!found)
	    (void)contig_register(io, cto, r->func, r->fdata, r->id,
				  r->flags, r->type);
    }

    /*
     * NOTE: The following line has been commented and uncommented many times.
     * Each time seemed things had changed slightly for reversing the comments
     * to appear to fix a bug. The current form (uncommented) should be
     * final. The change.log entry dated 22/04/96 states:
     *
     *	io-reg.c: Uncommented our "clear from list" bit of
     *	contig_register_join(). Note: (the following are NOT obvious!)
     *
     *	a) This is necessary in order to not free memory for plots that
     *	are now part of the 'to' contig when deleting the 'from' contig.?
     *	An example of this is when joining contig A with B when contig A
     *	or B is also shown in a stop codon plot, restriction enzyme
     *	plot, or (theoretically) any other non 2D plot.
     *
     *	b) A 2D plot is (and MUST be) always registered with all contigs.
     *	Hence when contig F is joined to contig T plots registered with
     *	contig N (not F or T) must still be replotted.
     *
     *	c) To achieve b) we make sure that creating a new contig will
     *	always call quit_displays first, which avoids the problem of
     *	having to reregister all the 2D plots.
     */
    /* clear cfrom list */
    io_Nreg(io, cfrom) = 0;
    
    /*
     * Update cursor lists.
     * This basically just concatenates the two cursors together.
     * If the sequence number for a cursor in 'cfrom' is 0 (consensus),
     * we also update the position.
     */
    for (gc = io_cursor(io, cto); gc && gc->next; gc = gc->next)
	;
    if (gc) {
	gc->next = io_cursor(io, cfrom);
    } else {
	io_cursor(io, cto) = io_cursor(io, cfrom);
    }
    for (gc = io_cursor(io, cfrom); gc; gc = gc->next) {
	if (gc->seq == 0) {
	    gc->pos += offset;
	    gc->abspos = gc->pos;
	} else {
	    gc->abspos = io_relpos(io, gc->seq) + gc->pos;
	}
    }
    io_cursor(io, cfrom) = NULL;

    /* Flag update for result-manager */
    update_results(io);

    return 0;
}


/*
 * Dump lists
 */
void contig_register_dump(GapIO *io) {
    contig_reg_t *r;
    int i, n, c;

    for (c = 0; c <= NumContigs(io); c++) {
	n = io_Nreg(io, c);
	r = ArrayBase(contig_reg_t, io_reg(io, c));
	
	printf("Contig %d\n", c);

	for (i = 0; i < n; i++)
	    printf("    Function 0x%p      Data 0x%p\n",
		   (void *)(r[i].func), r[i].fdata);
    }
}


/*
 *-----------------------------------------------------------------------------
 * Result (id) level actions
 *-----------------------------------------------------------------------------
 */

/*
 * Converts a registration id to an array of contig_reg_t structures,
 * suitable for further interrogation.
 *
 * Returns NULL for failure or a NULL terminated 'contig_reg_t *' array which
 * is expected to be xfree()d by the calling function.
 */
contig_reg_t **result_to_regs(GapIO *io, int id) {
    contig_reg_t **rl;
    int size = 8, count = 0, i, c;

    if (NULL == (rl = (contig_reg_t **)xmalloc(size * sizeof(contig_reg_t *))))
	return NULL;

    for (c = 0; c <= NumContigs(io); c++) {
	for (i = 0; i < io_Nreg(io, c); i++) {
	    contig_reg_t *r = &arr(contig_reg_t, io_reg(io, c), i);

	    if (r->id == id) {
		rl[count++] = r;
		if (count >= size-1) {
		    size *= 2;
		    if (NULL == (rl = (contig_reg_t **)xrealloc(rl,
					size * sizeof(contig_reg_t *)))) {
			xfree(rl);
			return NULL;
		    }
		}
	    }
	}
    }

    /* We've always guaranteed one element in rl to be left free. */
    rl[count] = NULL;

    return rl;
}


/*
 * Generates description of functions registered with a particular contig.
 * If contig -1 is specified then all are listed.
 * 'contig' is modified to return the contig number this result was from
 * (useful when sending contig -1), as is 'reg' to return the index into
 * the registration array for this contig. This (contig,reg) pair specifies
 * a particular result without the need for remembering pointers.
 * 'id' contains a unique id number for this result.
 *
 * Sending contig as 0 will list functions registered for all contigs. NB this
 * is not the same as contig -1, which lists all registrations.
 */
char *result_names(GapIO *io, int *contig, int *reg, int *id, int first) {
    reg_query_name qn;
    static char buf[80];
    static int n;
    static int c;
    contig_reg_t *r;

    qn.job = REG_QUERY_NAME;
    qn.line = buf;

    if (first) {
	n = 0;
	c = (*contig != -1) ? *contig : 0;
    } else {
	n++;
    }

    while (n >= io_Nreg(io, c)) {
	if (*contig == -1) {
	    n = 0;
	    if (++c > NumContigs(io))
		return NULL;
	} else
	    return NULL;
    }

    *qn.line = '\0';
    r = &arr(contig_reg_t, io_reg(io, c), n);
    if (r->flags & REG_FLAG_INVIS)
	*buf = 0;
    else
	r->func(io, c, r->fdata, (reg_data *)&qn);
    *contig = c;
    if (reg) *reg = n;
    if (id) *id = r->id;

    return qn.line;
}

/*
 * Returns a time (string) that a given id was registered. Assumes all
 * contig registrations for a particular id are registered together.
 * 'contig' isn't really needed here, but it's currently known by tcl
 * and speeds up our search.
 */
char *result_time(GapIO *io, int contig, int id) {
    int i, n;
    contig_reg_t *r;
    static char buf[80];
    
    n = io_Nreg(io, contig);
    r = ArrayBase(contig_reg_t, io_reg(io, contig));
    
    for (i = 0; i < n && r[i].id != id; i++)
	;

    if (i == n)
	return "unknown";

    /* %r doesn't work for windows ! */
    strftime(buf, sizeof(buf)-1, "%a %I:%M:%S %p", localtime(&r[i].time));
    return buf;
}


/*
 * Uses the register list for a given result to call a particular job.
 * As per contig_notify except on a result basis rather than contig basis.
 * 'all' declares whether to send message to all registrations with this id
 * or just the first found.
 */
void result_notify(GapIO *io, int id, reg_data *jdata, int all) {
    int i, c;

    for (c = 0; c <= NumContigs(io); c++) {
	for (i = 0; i < io_Nreg(io, c); i++) {
	    contig_reg_t *r = &arr(contig_reg_t, io_reg(io, c), i);

	    if (r->id == id && (r->flags & jdata->job)) {
		r->func(io, c, r->fdata, jdata);
		if (!all)
		    return;
	    }
	}
    }
}


/*
 * Returns the data component of a contig_reg_t for a specific id.
 * If contig is non zero, we search only this contig. Otherwise we scan
 * all.
 *
 * We return only the first data for id found, or NULL if none found.
 */
void *result_data(GapIO *io, int id, int contig) {
    int i, c, cend, cst;
    
    if (contig) {
	cend = contig;
	cst  = contig;
    } else {
	cst  = 0;
	cend = NumContigs(io);
    }

    for (c = cst; c <= cend; c++) {
	for (i = 0; i < io_Nreg(io, c); i++) {
	    contig_reg_t *r = &arr(contig_reg_t, io_reg(io, c), i);

	    if (r->id == id)
		return r->fdata;
	}
    }

    return NULL;
}


/*
 *-----------------------------------------------------------------------------
 * Type level actions
 *-----------------------------------------------------------------------------
 */

/*
 * Returns the id for a registered item of a particular type, If contig is
 * non zero, we search only this contig. Otherwise we scan all.
 *
 * We return only the first data for id found, or 0 if none found.
 */
int type_to_result(GapIO *io, int type, int contig) {
    int i, c, cend, cst;
    
    if (contig) {
	cend = contig;
	cst  = contig;
    } else {
	cst  = 0;
	cend = NumContigs(io);
    }

    for (c = cst; c <= cend; c++) {
	for (i = 0; i < io_Nreg(io, c); i++) {
	    contig_reg_t *r = &arr(contig_reg_t, io_reg(io, c), i);

	    if (r->type == type)
		return r->id;
	}
    }

    return 0;
}


/*
 * Notifies all (or the first) registered items of a given type.
 * Works across all contigs.
 *
 * Returns 0 for success, -1 for failure.
 */
int type_notify(GapIO *io, int type, reg_data *jdata, int all) {
    int i, j, k, contig, ret = -1;
    int *uids = NULL, n;
    contig_reg_t *r;

    for (contig = 0; contig <= NumContigs(io); contig++) {
	n = io_Nreg(io, contig);
	r = ArrayBase(contig_reg_t, io_reg(io, contig));

	if (0 == n)
	    continue;

	/*
	 * Take a snap shot of all the registrations that need notifying.
	 * We do this by remembering their unique uids.
	 */
	if (NULL == (uids = (int *)xrealloc(uids, n * sizeof(int))))
	    return -1;

	for (k = 0; k < n; k++)
	    uids[k] = r[k].uid;

	/* Now loop through all uids sending the notification */
	for (j = i = 0; j < k; i++, j++) {
	    /*
	     * Don't optimise this outside the loop - it may not be a constant
	     */
	    n = io_Nreg(io, contig);
	    
	    /*
	     * It's possible that the next registration isn't our next uid.
	     * This occurs when the previous r[i].func() modified the
	     * registration list for this contig (eg by deregistering
	     * something). So we scan along to find the correct 'i' in this
	     * case remembering that it may not even exist any more.
	     */
	    if (i >= n || r[i].uid != uids[j]) {
		for (i = 0; i < n; i++) {
		    if (r[i].uid == uids[j])
			break;
		}
		if (i == n)
		    continue;
	    }

	    /* Send the notification request itself */
	    if (r[i].type == type && (r[i].flags & jdata->job)) {
		r[i].func(io, contig, r[i].fdata, jdata);
		if (!all) {
		    xfree(uids);
		    return 0;
		} else {
		    ret=0;
		}
	    }
	}
    }

    xfree(uids);

    return ret;
}


/*
 * Notifies all (or the first) registered items of a given type within
 * a specified contig.
 *
 * Returns 0 for success, -1 for failure.
 */
int type_contig_notify(GapIO *io, int contig, int type,
			reg_data *jdata, int all) {
    contig_reg_t *r;
    int i, n, ret = -1;

    n = io_Nreg(io, contig);
    r = ArrayBase(contig_reg_t, io_reg(io, contig));
    
    for (i = 0; i < n; i++) {
	if (r[i].type == type && (r[i].flags & jdata->job)) {
	    r[i].func(io, contig, r[i].fdata, jdata);
	    if (!all)
		return 0;
	    else
		ret=0;
	}
    }

    return ret;
}


/*
 *-----------------------------------------------------------------------------
 * Lock management
 *-----------------------------------------------------------------------------
 */

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

/*
 *-----------------------------------------------------------------------------
 * Cursor functions
 *-----------------------------------------------------------------------------
 */
static int cursor_id = 0;

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
	for (gc = io_cursor(io, contig); gc; gc = gc->next)
	    if (!gc->private)
		break;
	if (gc) {
	    gc->private = private;
	    gc->refs++;
	    goto notify;
	}
    }

    /* If not private, use the first cursor we find */
    if (!private && io_cursor(io, contig)) {
	io_cursor(io, contig)->refs++;
	gc = io_cursor(io, contig);
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
    gcend = io_cursor(io, contig);
    while(gcend && gcend->next)
	gcend = gcend->next;

    /* Add gc to it */
    if (gcend)
	gcend->next = gc;
    else
	io_cursor(io, contig) = gc;

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
 * found. If contig != NULL, only look in this contig. Otherwise look at all.
 * If found and contig != NULL, fill with the correct contig number.
 *
 * Always sends a notification about this cursor to inform of the new
 * reference count, and possible private status.
 */
cursor_t *find_contig_cursor(GapIO *io, int *contig, int id)
{
    cursor_t *gc;
    int c;

    if (contig && *contig != 0) {
	for (gc = io_cursor(io, *contig); gc && gc->id != id; gc = gc->next);
	if (gc && gc->id == id)
	    return gc;
	return NULL;
    }

    for (c = 1; c <= NumContigs(io); c++) {
	if (contig)
	    *contig = c;
	for (gc = io_cursor(io, c); gc && gc->id != id; gc = gc->next);
	if (gc && gc->id == id)
	    return gc;
    }
    return NULL;
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

    if (!(gc = find_contig_cursor(io, &c, id)))
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
    if (io_cursor(io, c) == gc) {
	io_cursor(io, c) = gc->next;
	xfree(gc);
	return;
    }

    /* Otherwise, find previous and link to next*/
    for (gcp = io_cursor(io, c); gcp && gcp->next != gc; gcp = gcp->next);
    if (gcp && gcp->next == gc) {
	gcp->next = gc->next;
	xfree(gc);
    }

    return;
}
