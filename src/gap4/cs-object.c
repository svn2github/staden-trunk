#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xalloc.h"
#include "cs-object.h"
#include "newgap_cmds.h" /* GetInterp() */
#include "misc.h"
#include "text_output.h"
#include "contig_selector.h"
#include "tcl_utils.h"
#include "tk-io-reg.h"

/*
 * ============================================================================
 * Hashing for translation of tk canvas IDs to generic objects
 * ============================================================================
 */

#define Hash(key) ((key) % HASHMODULUS)

void InitialiseHash(HTablePtr T[]) {
    int i;

    for (i = 0; i < HASHMODULUS; i++)
	T[i] = NULL;

} /* end InitialiseTable */


void HashInsert(HTablePtr T[], int newkey, HItemType newinfo) {
    int H_index;
    HTablePtr new_node;

    H_index = Hash(newkey);
    /* Allocate and check memory */
    if (NULL == (new_node = (HTableItem *)xmalloc(sizeof(HTableItem))))
	return;

    new_node->key = newkey;
    new_node->info = newinfo;
    new_node->next = T[H_index];

    T[H_index] = new_node;

} /* end HashInsert */


void HashDelete(HTablePtr T[], int search_key) {
    int H_index = Hash(search_key);
    HTablePtr node_ptr = T[H_index], last_ptr = NULL;

    /* Find node */
    for (; node_ptr != NULL; last_ptr = node_ptr, node_ptr = node_ptr->next) {
	if (node_ptr->key == search_key)
	    break;
    }

    if (node_ptr) {
	/* Unlink */
	if (last_ptr) {
	    last_ptr->next = node_ptr->next;
	} else {
	    T[H_index] = node_ptr->next;
	}

	/* and free */
	xfree(node_ptr);
    }
} /* end HashDelete */


HItemType *HashSearch(HTablePtr T[], int search_key) {
    HTablePtr node_ptr = T[Hash(search_key)];

    for (; node_ptr != NULL; node_ptr = node_ptr->next) {
	if (node_ptr->key == search_key)
	    break;
    }

    return node_ptr ? node_ptr->info : NULL;

} /* end HashSearch */


void csmatch_reset_hash(HTablePtr T[], mobj_repeat *r) {
    int i;
    for (i = 0; i < r->num_match; i++) {
	HashDelete(csplot_hash, r->match[i].inum);
	HashInsert(csplot_hash, r->match[i].inum, &r->match[i]);
    }
}

/*
 * ============================================================================
 * Manipulation of entire sets of results.
 * ============================================================================
 */
void DeleteRepeats(Tcl_Interp *interp, mobj_repeat *r, char *csplot_name,
		   HTablePtr T[]) {
    int i;

    /* Loop through each item removing from the hash table */
    for (i = r->num_match-1; i >= 0; i--) {
	HashDelete(T, r->match[i].inum);
    }

    /* Remove from canvas */
    Tcl_VarEval(interp, csplot_name, " delete ", r->tagname, NULL);
}


/*
 * ============================================================================
 * Object callback mechanism - 'dehashes' and invokes the appropriate jobs
 * ============================================================================
 */
char *obj_get_ops(int inum) {
    extern HTablePtr csplot_hash[HASHMODULUS];
    obj_generic *obj;
    char *ops;

    if (NULL == (obj = (obj_generic *)HashSearch(csplot_hash, inum))) {
	verror(ERR_FATAL, "obj_get_ops",
	       "Unknown canvas item number! (%d)\n", inum);
	return NULL;
    }

    if (!obj->call.func)
	return NULL;

    ops = obj->call.func(OBJ_LIST_OPERATIONS, NULL, obj, obj->call.data);

/*
    printf("\n==== Operations allowed ====\n");
    while (l = strlen(ops)) {
	printf("%d. %s\n", op++, ops);
	ops += l+1;
    }
    printf("\n");
*/
    return ops;
}

void obj_invoke_op(int inum, int op) {
    extern HTablePtr csplot_hash[HASHMODULUS];
    obj_generic *obj;

    if (NULL == (obj = (obj_generic *)HashSearch(csplot_hash, inum))) {
	verror(ERR_FATAL, "obj_invoke_op",
	       "Unknown canvas item number! (%d)\n", inum);
	return;
    }

    if (obj->call.func)
	obj->call.func(OBJ_INVOKE_OPERATION, &op, obj, obj->call.data);
}

char *obj_get_brief(int inum) {
    extern HTablePtr csplot_hash[HASHMODULUS];
    obj_generic *obj;

    if (NULL == (obj = (obj_generic *)HashSearch(csplot_hash, inum))) {
	verror(ERR_FATAL, "obj_get_brief",
	       "Unknown canvas item number! (%d)\n", inum);
	return NULL;
    }

    if (obj->call.func)
	return obj->call.func(OBJ_GET_BRIEF, NULL, obj, obj->call.data);
    else
	return NULL;
}


/*
 * Hides a match object from a repeat style metaobject
 */
void obj_hide(Tcl_Interp *interp, char *cs_plot, obj_match *obj,
	      mobj_repeat *r, HTablePtr T[]) {
    /*
     * Set this object to be hidden, and replot them all. We could replot
     * only this one, but we don't know its item number and it's fiddly
     * to find out.
     */
    if (obj->flags & OBJ_FLAG_HIDDEN)
	return;

    obj->flags |= OBJ_FLAG_HIDDEN;

    DeleteRepeats(interp, r, cs_plot, T);
    PlotRepeats(r->io, r);
}

/*
 * Shows a match object (after a hide) from a repeat style metaobject.
 * NB: Never used at present as we can't highlight this object to reveal!
 */
void obj_reveal(Tcl_Interp *interp, char *cs_plot, obj_match *obj,
		mobj_repeat *r, HTablePtr T[]) {
    /*
     * Set this object to be visible, and replot them all. We could replot
     * only this one, but we don't know its item number and it's fiddly
     * to find out.
     */
    if (!(obj->flags & OBJ_FLAG_HIDDEN))
	return;

    obj->flags &= ~OBJ_FLAG_HIDDEN;

    DeleteRepeats(interp, r, cs_plot, T);
    PlotRepeats(r->io, r);
}

/*
 * Removes a single match object and replots.
 */
void obj_remove(Tcl_Interp *interp, char *cs_plot, obj_match *obj,
		mobj_repeat *r, HTablePtr T[]) {
    DeleteRepeats(interp, r, cs_plot, T);

    /* Sorry about the arithmetic! I'm feeling terse today */
    memmove(obj, obj+1, (--r->num_match - (obj - r->match)) * sizeof(*obj));

    if (r->num_match > 0) {
	PlotRepeats(r->io, r);
    } else {
	csmatch_remove(r->io, cs_plot, r, T);
    }
}

void obj_invoke_next(mobj_repeat *mobj) {
    csmatch_invoke_next(mobj);
}

int obj_get_next(mobj_repeat *mobj) {
    int next;

    next = csmatch_get_next(mobj);
    if (next != -1)
	return mobj->match[next].inum;
    else
	return -1;
}

/*
 * ============================================================================
 * Manipulation of contig selector objects for various REG operations
 * ============================================================================
 */


/*
 * Handle a REG_JOIN_TO request for match objects
 * 'r' isn't always a mobj_repeat, but sometimes a mobj_template or
 * mobj_fij. For the time being though these are all equivalent.
 */
void csmatch_join_to(GapIO *io, int contig, reg_join *j, mobj_repeat *r,
		     HTablePtr T[], char *cs_plot) {
    int i;
    /*
       printf("Joining %d to %d at offset %d\n",
       contig, j->contig, j->offset);
       */
    for (i = 0; i < r->num_match; i++) {
	if (abs(r->match[i].c1) == contig) {
	    r->match[i].pos1 += j->offset;
	    r->match[i].c1 = r->match[i].c1 > 0
		? j->contig : -j->contig;
	}

	if (abs(r->match[i].c2) == contig) {
	    r->match[i].pos2 += j->offset;
	    r->match[i].c2 = r->match[i].c2 > 0
		? j->contig : -j->contig;
	}

	/* For FIJ: remove match if moved onto diagonal */
	if (r->match_type == REG_TYPE_FIJ &&
	    r->match[i].c1 == r->match[i].c2) {
	    obj_match *o = &r->match[i];

	    if (i <= r->current) r->current--;
	    i--;
	    memmove(o, o+1, (--r->num_match - (o - r->match)) * sizeof(*o));
	}
    }

    if (r->num_match > 0) {
	DeleteRepeats(GetInterp(), r, cs_plot, T);
	PlotRepeats(io, r);
    } else {
	csmatch_remove(io, cs_plot, r, T);
    }

    return;
}

/*
 * Handle a REG_COMPLEMENT request for match objects
 */
void csmatch_complement(GapIO *io, int contig, mobj_repeat *r,
			HTablePtr T[], char *cs_plot) {
    int n, i;
    int clen = io_clength(io, contig);

    n = r->num_match;
    for (i = 0; i < n; i++) {
	if (abs(r->match[i].c1) == contig) {
	    r->match[i].pos1 = clen+1 - (r->match[i].pos1 +
					 r->match[i].length - 1);

	    r->match[i].c1 = -r->match[i].c1;
	}

	if (abs(r->match[i].c2) == contig) {
	    r->match[i].pos2 = clen+1 - (r->match[i].pos2 +
					 r->match[i].length - 1);

	    r->match[i].c2 = -r->match[i].c2;
	}
    }

    /* FIXME */
    DeleteRepeats(GetInterp(), r, cs_plot, T);

    PlotRepeats(io, r);

    return;
}


/*
 * Pop up a configuration window for adjusting a plot.
 *
 * NOTE: This is ghastly. We write the mobj_repeat pointer as a text string
 * using the %p printf format. This is then passed into Tcl, which faithfully
 * passes it around not knowing what it means until it gets back into C
 * at the tk_matchresult_configure function, whereupon it is converted back.
 * NB: This relies on mobj_repeat never been reallocated.
 */
void csmatch_configure(GapIO *io, char *cs_plot, mobj_repeat *r) {
    char *tclptr;

    tclptr = CPtr2Tcl(r);
    if (TCL_OK != Tcl_VarEval(GetInterp(), "cs_config ", cs_plot, " ",
			      tclptr, NULL)) {
	puts(GetInterpResult());
    }
}


/*
 * Remove a result from the csplot window
 */
void csmatch_remove(GapIO *io, char *cs_plot,
		    mobj_repeat *reg_dat,
		    HTablePtr T[]) {
    int c;

    /* Delete from the canvas and hash table */
    DeleteRepeats(GetInterp(), reg_dat, cs_plot, T);

    /*
     * Remove from the registration lists.
     * Loop through all contigs for time being.
     */
    for (c = 1; c <= NumContigs(io); c++)
	contig_deregister(io, c, reg_dat->reg_func, reg_dat);

    /*
     * Pop down configuration window if visible
     */
    if (TCL_OK != Tcl_VarEval(GetInterp(), "cs_config_quit ", cs_plot, " ",
			      reg_dat->tagname, NULL)) {
	puts(GetInterpResult());
    }

    /* Inform contig selector next button */
    Tcl_VarEval(GetInterp(), "CSLastUsedFree ", CPtr2Tcl(reg_dat), NULL);

    /* Free memory */
    if (reg_dat->match)
	xfree(reg_dat->match);
    if (reg_dat->params)
	xfree(reg_dat->params);
    xfree(reg_dat);
}


/*
 * Provide information about a set of matches - the meta object.
 */
void csmatch_info(mobj_repeat *r, char *name) {
    int i;

    vfuncheader("%s result list", name);

    vmessage("Number of matches = %d\n", r->num_match);
    vmessage("Display colour = %s\n", r->colour);

    for (i = 0; i < r->num_match; i++) {
	obj_match *obj = &r->match[i];
	int op = -1; /* For match operations, -1 is always "Information" */

	vmessage("\nMatch %d:\n", i);
	obj->func(OBJ_INVOKE_OPERATION, &op, obj, obj->data);
    }
}

/*
 * Invoke the default operation on the next match in a meta-object.
 */
void csmatch_invoke_next(mobj_repeat *r) {
    int op = -2; /* Default */
    int next;
    void bell(void);

    next = csmatch_get_next(r);
    if (next == -1) {
	bell();
	return;
    }

    r->current = next;
    r->match[r->current].flags |= OBJ_FLAG_VISITED;
    r->match[r->current].func(OBJ_INVOKE_OPERATION, &op,
			      &r->match[r->current], r);
}

int csmatch_get_next(mobj_repeat *r) {
    int next = r->current, countdown = r->num_match;

    do {
	next++;
	if (next >= r->num_match)
	    next = 0;

	countdown--;
    } while ((r->match[next].flags & OBJ_FLAG_VISITED) && countdown >= 0);

    if (countdown < 0)
	return -1;
    else
	return next;
}

/*
 * Reveals all match objects in a repeat style metaobject. This is currently
 * the only way to reshow objects as we haven't got a user interface to
 * obj_reveal().
 */
void csmatch_reveal(Tcl_Interp *interp, char *cs_plot, mobj_repeat *r,
		    HTablePtr T[]) {
    int i;

    for (i = 0; i < r->num_match; i++)
	r->match[i].flags &= ~OBJ_FLAG_HIDDEN;

    DeleteRepeats(interp, r, cs_plot, T);
    PlotRepeats(r->io, r);

    r->all_hidden = 0;
    update_results(r->io);
}


/*
 * Hides all match objects in a repeat style metaobject without removing
 * them. Reveal later using csmatch_reveal().
 */
void csmatch_hide(Tcl_Interp *interp, char *cs_plot, mobj_repeat *r,
		  HTablePtr T[]) {
    int i;

    for (i = 0; i < r->num_match; i++)
	r->match[i].flags |= OBJ_FLAG_HIDDEN;

    DeleteRepeats(interp, r, cs_plot, T);
    PlotRepeats(r->io, r);

    /*
     * We also need to shut down the configure window as this doesn't work
     * with everything hidden.
     */
    Tcl_VarEval(interp, "cs_config_quit ", cs_plot, " ", r->tagname, NULL);

    r->all_hidden = 1;
    update_results(r->io);
}

/*
 * Handle a REG_NUMBER_CHANGE request for match objects
 */
void csmatch_renumber(GapIO *io, int old_contig, int new_contig,
		      mobj_repeat *r, HTablePtr T[], char *cs_plot) {
    int n, i;

    n = r->num_match;
    for (i = 0; i < n; i++) {
	if (abs(r->match[i].c1) == old_contig)
	    r->match[i].c1 = r->match[i].c1 > 0 ? new_contig : -new_contig;

	if (abs(r->match[i].c2) == old_contig)
	    r->match[i].c2 = r->match[i].c2 > 0 ? new_contig : -new_contig;
    }

    /* FIXME */
    DeleteRepeats(GetInterp(), r, cs_plot, T);

    PlotRepeats(io, r);

    return;
}

/*
 * Replot match objects, eg when a contig length has changed.
 */
void csmatch_replot(GapIO *io, mobj_repeat *r, HTablePtr T[], char *cs_plot) {
    /* FIXME */
    DeleteRepeats(GetInterp(), r, cs_plot, T);
    PlotRepeats(io, r);

    return;
}

/*
 * Handle a REG_DELETE request for match objects
 */
void csmatch_contig_delete(GapIO *io, mobj_repeat *r, int contig,
			   char *cs_plot, HTablePtr T[]) {
    int i, n;

    n = r->num_match;
    for (i = 0; i < n; i++) {
	if (abs(r->match[i].c1) == contig || abs(r->match[i].c2 == contig)) {
	    /*
	     * Found a match to be removed - we copy the last match in our
	     * array to this match and decrement our num_match count
	     */
	    if (i < n-1) {
		memcpy(&r->match[i], &r->match[n-1], sizeof(r->match[i]));
		i--;
	    }
	    n--;
	}
    }
    r->num_match = n;

    /* FIXME */
    DeleteRepeats(GetInterp(), r, cs_plot, T);

    PlotRepeats(io, r);

    return;
}



/*
 * Reset the visited flags - used by the 'next' mechanism
 */
void csmatch_reset_next(mobj_repeat *r) {
    int i;

    for (i = 0; i < r->num_match; i++) {
	r->match[i].flags &= ~OBJ_FLAG_VISITED;
    }
    r->current = -1;
}
