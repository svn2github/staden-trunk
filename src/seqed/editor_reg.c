#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "array.h"
#include "seq_reg.h"
#include "editor_reg.h"
#include "misc.h"
#include "tkRaster.h"
#include "seq_raster.h"
#include "spin_globals.h"
#include "tcl_utils.h"


static Array editor_reg;
static Array editor_cursor_reg;/*FIXME*/

extern void init_cursor_colour(Tcl_Interp *interp);
extern int add_cursor_free_array(int id);
extern int get_cursor_id(void);
extern char *get_cursor_colour(int id);

/* Initialise the editor registration lists */
int editor_register_init (Tcl_Interp *interp) {

    /* Allocate Array of arrays */
    if (NULL == (editor_reg = ArrayCreate(sizeof(Array), 0)))
	return -1;
    if (NULL == (editor_cursor_reg = ArrayCreate(sizeof(Array), 0)))
	return -1;
    init_cursor_colour(interp);
    return 0;
}

/*
 * Add a new sequence to the editor registration list 
 */
int add_editor_reg(int index) {
    
    /* return ref to memory and allocate more memory if necessary */
    (void)ArrayRef(editor_reg, index);
    (void)ArrayRef(editor_cursor_reg, index);

    if (NULL == (editor_func_array(index) = ArrayCreate(sizeof(edit_reg), 0)))
    	return -1;
    editor_Nfuncs(index) = 0;
    editor_cursor(index) = NULL;
    return 0;
}

void remove_sequence_from_registration (int index) {

    /* remove sequence from array structure */
    if (index < (ArrayMax(editor_reg) - 1)) {
	memmove(arrp(Array, editor_reg, index),
		arrp(Array, editor_reg, index+1),
		(ArrayMax(editor_reg) - index - 1) * sizeof(Array));
	memmove(arrp(Array, editor_cursor_reg, index),
		arrp(Array, editor_cursor_reg, index+1),
		(ArrayMax(editor_cursor_reg) - index - 1) * sizeof(Array));
    }
    ArrayMax(editor_reg)--;
    ArrayMax(editor_cursor_reg)--;

}

int get_editor_reg_id() {

    static int ed_reg_id = 0;
    return(ed_reg_id++);
}

int editor_register (int seq_num,
		  void (*func)(int seq_num, void *fdata, editor_reg_data *jdata),
		  void *fdata,
		  int type,
		  int id)
{
    edit_reg *r;
    int i;

    /* Check if this (func,fdata) pair already exists. */
    r = ArrayBase(edit_reg, editor_func_array(seq_num));
    for (i = 0; i < editor_Nfuncs(seq_num); i++) {
	if (r[i].func == func && r[i].fdata == fdata) {
	    return 0;
	}
    }
    /* Allocate new function item */
    if (NULL == (r = (edit_reg *)ArrayRef(editor_func_array(seq_num), 
					 editor_Nfuncs(seq_num))))
	return -1;
    /* Copy over our data */
    r->func = func;
    r->fdata = fdata;
    r->time = time(NULL);
    r->type = type;
    r->id = id;

    return 0;
}
/*
 * Deregisters func(seq_num, data, jdata) from sequence 'seq_num'.
 *
 * Returns 0 for success, and -1 for error.
 */
int 
editor_deregister(int seq_num,
		  void (*func)(int seq_num, void *fdata, editor_reg_data *jdata),
		  void *fdata) 
{
    edit_reg *r;
    int i;
    int num_funcs;

#ifdef DEBUG
    printf("DEREGISTER seq %d \n", seq_num);
    /*seq_register_dump();*/
#endif
    num_funcs = editor_Nfuncs(seq_num);
    r = ArrayBase(edit_reg, editor_func_array(seq_num));
    /* Search for element in the array */

    /* 
     * can have the same function registered twice with the same sequence eg
     * sip4 with the same sequence. Therefore must go through all the functions
     * registered and not stop when find the first one
     */
    for (i = 0; i < num_funcs; i++) {
	if (r[i].func == func && r[i].fdata == fdata) {

	    /* Match found - shuffle everything else down */
	    memmove(&r[i], &r[i+1], (editor_Nfuncs(seq_num) - i - 1) * sizeof(r[i]));
	    
	    /* Decrement count */
	    editor_Nfuncs(seq_num)--;

	    /* decrement i and num_funcs as moved everything up in the list */
	    i--;
	    num_funcs--;
	}
    }

#ifdef REMOVE
    for (i = 0; i < num_funcs; i++) {
	if (r[i].func == func && r[i].fdata == fdata)
	    break;
    }

    if (i < num_funcs) {
	/* Match found - shuffle everything else down */
	memmove(&r[i], &r[i+1], (editor_Nfuncs(seq_num) - i - 1) * sizeof(r[i]));

	/* Decrement count */
	editor_Nfuncs(seq_num)--;
    }
#endif
    return 0;
}

/*
 * Uses the register list for a given seq_num to call a particular job.
 */
void editor_notify(int seq_num, editor_reg_data *jdata) {

    edit_reg *r;
    int i, j, k;
    int num_funcs;
    int *ids;

#ifdef DEBUG
    seq_register_dump();
#endif

    num_funcs = editor_Nfuncs(seq_num);
    r = ArrayBase(edit_reg, editor_func_array(seq_num));

    if (0 == num_funcs)
	return;

    /*
     * Take a snap shot of all the registrations that need notifying.
     * We do this by remembering their unique uids.
     */
    if (NULL == (ids = (int *)xmalloc(num_funcs * sizeof(int))))
	return;

    for (k = 0; k < num_funcs; k++)
	ids[k] = r[k].id; 

     /* Now loop through all uids sending the notification */
    for (j = i = 0; j < k; i++, j++) {
	/* Don't optimise this outside the loop - it may not be a constant */
	num_funcs = editor_Nfuncs(seq_num);
     
	/*
	 * It's possible that the next registration isn't our next id.
	 * This occurs when the previous r[i].func() modified the registration
	 * list for this sequence (eg by deregistering something).
	 * So we scan along to find the correct 'i' in this case remembering
	 * that it may not even exist any more.
	 */
	if (i >= num_funcs || r[i].id != ids[j]) {
	    for (i = 0; i < num_funcs; i++) {
		if (r[i].id == ids[j])
		    break;
	    }
	    if (i == num_funcs) {
		continue;
	    }	
	}
	/* Send the notification request itself */
	r[i].func(seq_num, r[i].fdata, jdata);
    }
    xfree(ids);
}

/*
 * Uses the register list for a given result to call a particular job.
 */
void 
editor_result_notify(int id, editor_reg_data *jdata, int all) {
    
    int i, j, k, s, n;
    int *ids;
    edit_reg *r;

    for (s = 0; s < ArrayMax(editor_reg); s++) {
	n = editor_Nfuncs(s);
	r = ArrayBase(edit_reg, editor_func_array(s));
	
#ifdef DEBUG
	printf("seq_result_notify s %d n %d\n", s, n);
#endif
	if (0 == n)
	    continue;
	
	/*
	 * Take a snap shot of all the registrations that need notifying.
	 * We do this by remembering their unique uids.
	 */
	if (NULL == (ids = (int *)xmalloc(n * sizeof(int))))
	    return;

	for (k = 0; k < n; k++)
	    ids[k] = r[k].id;
	
	/* Now loop through all uids sending the notification */
	for (j = i = 0; j < k; i++, j++) {
	    /* Don't optimise this outside the loop - it may not be a constant */
	    n = editor_Nfuncs(s);

	    /*
	     * It's possible that the next registration isn't our next id.
	     * This occurs when the previous r[i].func() modified the registration
	     * list for this sequence (eg by deregistering something).
	     * So we scan along to find the correct 'i' in this case remembering
	     * that it may not even exist any more.
	     */
	    if (i >= n || r[i].id != ids[j]) {
		for (i = 0; i < n; i++) {
		    if (r[i].id == ids[j])
			break;
		}
		if (i == n)
		    continue;
	    }
	    
	    /* Send the notification request itself */
	    if (r[i].id == id) {
		r[i].func(s, r[i].fdata, jdata);
		if (!all) {
		    xfree(ids);
		    return;
		}

	    }
	}
	xfree(ids);
    }

#ifdef REMOVE
    for (s = 0; s < ArrayMax(sequence_reg); s++) {
	n = editor_Nfuncs(s);
	for (i = 0; i < n; i++) {
	    edit_reg *r = &arr(edit_reg, editor_func_array(s), i);
	    if (r->id == id) {
		r->func(s, r->fdata, jdata);
		if (!all) {
		    return;
		}
		
	    }
	}
    }
#endif
}


/*
 * Create a cursor for this sequence.
 * If private == 1 then create a new cursor if no non-private ones can be
 * found. Otherwise use an existing one if available.
 * Ie, private cursors can be used more than once, but only by one private
 * user at any one time.
 * Returns the cursor pointer, or NULL for failure.
 */
cursor_e *editor_create_cursor(int seq_num, 
			       int private, 
			       char *colour,
			       int line_width,
			       int cursor_num,
			       int direction)
{
    cursor_e *gc, *gcend;
    seq_reg_cursor_notify cn;

    /* If private, look for a non private cursor to use */
    if (private) {	
	for (gc = editor_cursor(seq_num); gc; gc = gc->next)	 
	if (!gc->private && gc->direction == direction)	   
		break;	
	if (gc) {
	    gc->private = private;
	    gc->refs++;
	    goto notify;
	}
    }

    /* If not private, use the first cursor we find */
#if REMOVE
    if (!private && editor_cursor(seq_num)) {
	editor_cursor(seq_num)->refs++;
	gc = editor_cursor(seq_num);
	goto notify;
    }
#endif
    if (!private) {
	for (gc = editor_cursor(seq_num); gc; gc = gc->next) {
	    if (gc->direction == direction) {
		--cursor_num; /* cursor_num? what is useful */
	    }
	    if (cursor_num <= 0) 
		break;
	}
	if (gc /* && gc->direction == direction*/) {
	    gc->refs++;
	    goto notify;
	}
    }

    /* Otherwise, create our own */
    if (NULL == (gc = xmalloc(sizeof(*gc))))
	return NULL;
    gc->id = get_cursor_id();
    if (gc->id >= NUM_CURSORS) {
	verror(ERR_WARN, "create cursor", "Too many cursors\n");
	return NULL;
    }

    gc->refs = 1;
    gc->abspos = 1;
    gc->posy = 1;
    gc->private = private;
    gc->next = NULL;
    
    /* allocate colour */
    if (!colour) {
	gc->colour = strdup(get_cursor_colour(gc->id));
    } else {
	gc->colour = strdup(colour);
    }
    gc->line_width = line_width;
    gc->direction = direction;
    gc->sent_by = 0;
    gc->job = 0;

    /* Find last existing cursor */
    gcend = editor_cursor(seq_num);
    while(gcend && gcend->next)
	gcend = gcend->next;

    /* Add gc to it */
    if (gcend)
	gcend->next = gc;
    else
	editor_cursor(seq_num) = gc;

 notify:
    cn.cursor = gc;
    cn.job = SEQ_CURSOR_NOTIFY;
    gc->job = CURSOR_MOVE | CURSOR_INCREMENT;

#ifdef DEBUG
    printf("CREATE_CURSOR num %d id %d refs %d dir %d private %d cursor_num %d\n", 
	   seq_num, gc->id, gc->refs, gc->direction, gc->private, cursor_num);
#endif
    /*
      HACK - is this needed? Causes problems with redrawing raster cursors
      */
    editor_notify(seq_num, (editor_reg_data *)&cn);
    return gc;
}


/*
 * Given a cursor identifier, return the cursor structure or NULL if not
 * found. If seq_num != NULL,only look in this seq_num. Otherwise look at all.
 * If found and seq_num != NULL, fill with the correct seq number.
 *
 * Always sends a notification about this cursor to inform of the new
 * reference count, and possible private status.
 */
cursor_e *ed_find_cursor(int *seq_num, int id, int direction)
{
    cursor_e *gc;
    int s;

    if (seq_num && *seq_num != -1) {
	for (gc = editor_cursor(*seq_num); gc; gc = gc->next) {
	    if (gc && gc->id == id) {
		if (direction == -1) {
		    return gc;
		} else if (gc->direction == direction) {
		    return gc;
		}
	    }
	}
	return NULL;
    }

    for (s = 0; s < ArrayMax(editor_reg); s++) {
	if (seq_num)
	    *seq_num = s;
	for (gc = editor_cursor(s); gc; gc = gc->next) {
	    if (gc && gc->id == id) {
		if (direction == -1) {
		    return gc;
		} else if (gc->direction == direction) {
		    return gc;
		}
	    }
	}

    }
    return NULL;
}

/*
 * Deletes a sequence cursor for this option. If the cursor is in use
 * more than once then this simply decrements the reference count.
 * 'private' indicates whether this cursor was "your private one".
 *
 * Always sends a notification about this cursor, either to say it is
 * being destroyed (abspos = -1), or to say that the reference count
 * (and possible private status) have changed.
 */

void ed_delete_cursor(int seq_num, int id, int private)
{
    cursor_e *gc, *gcp;
    int s = seq_num;
    seq_reg_cursor_notify cn;
#ifdef DEBUG    
    printf("SEQ delete_cursor seq_num %d id %d private %d\n",
	   seq_num, id, private);
#endif
    if (!(gc = ed_find_cursor(&s, id, -1)))
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
    cn.job = SEQ_CURSOR_NOTIFY;
    editor_notify(s, (editor_reg_data *)&cn);

    if (gc->refs > 0)
	return;

    /* If first in seq_num, skip to next */
    if (editor_cursor(s) == gc) {
	editor_cursor(s) = gc->next;
	add_cursor_free_array(gc->id);
	free(gc->colour);
	xfree(gc);
	return;
    }

    /* Otherwise, find previous and link to next*/
    for (gcp = editor_cursor(s); gcp && gcp->next != gc; gcp = gcp->next);
    if (gcp && gcp->next == gc) {
	gcp->next = gc->next;
	add_cursor_free_array(gc->id);
	free(gc->colour);
	xfree(gc);
    }
    return;
}

text_editor_result *init_text_editor_result (void) {

    text_editor_result *r;
    
    if (NULL == (r = (text_editor_result *)xmalloc(sizeof(text_editor_result))))
	return NULL;
    r->op_func = NULL;
    r->seq_id = 0;
    r->ed_id = 0;
    r->interp = NULL;
    r->colour = NULL;
    r->index = 0; /* event register ID */
    r->cursor_id = 0; /*cursor register ID */
    return r;
}

/*seq_reg_changed *init_sequence_changed (void) {



}*/
