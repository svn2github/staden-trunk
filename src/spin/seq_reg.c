#include <stdio.h>
#include <string.h>
#include <time.h>

#include "array.h"
#include "seq_reg.h"
#include "misc.h"
#include "tkRaster.h"
#include "seq_raster.h"
#include "spin_globals.h"
#include "tcl_utils.h"

/* 
 * overall structure:
 * array of size num sequences each of which has an array of associated
 * functions eg
 * seq[0] may have compare spans
 *                 hashing
 * seq[1] may have compare spans
 * seq[2] may have hashing
 *
 */

static Array sequence_reg;
static Array seq_cursor_reg;

#define NUM_CURSOR_COL 10
static char *cursor_colour[NUM_CURSOR_COL];

void init_cursor_colour(Tcl_Interp *interp)
{
    cursor_colour[0] = get_default_string(interp, spin_defs, "CURSOR.COLOUR.0");
    cursor_colour[1] = get_default_string(interp, spin_defs, "CURSOR.COLOUR.1");
    cursor_colour[2] = get_default_string(interp, spin_defs, "CURSOR.COLOUR.2");
    cursor_colour[3] = get_default_string(interp, spin_defs, "CURSOR.COLOUR.3");
    cursor_colour[4] = get_default_string(interp, spin_defs, "CURSOR.COLOUR.4");
    cursor_colour[5] = get_default_string(interp, spin_defs, "CURSOR.COLOUR.5");
    cursor_colour[6] = get_default_string(interp, spin_defs, "CURSOR.COLOUR.6");
    cursor_colour[7] = get_default_string(interp, spin_defs, "CURSOR.COLOUR.7");
    cursor_colour[8] = get_default_string(interp, spin_defs, "CURSOR.COLOUR.8");
    cursor_colour[9] = get_default_string(interp, spin_defs, "CURSOR.COLOUR.9");
}



void seq_register_dump(void);

/*
 *-----------------------------------------------------------------------------
 * (De)registration functions
 *-----------------------------------------------------------------------------
 */

/*
 * Initialise the sequence register lists
 *
 * Returns 0 on success and -1 for error;
 */
int seq_register_init(Tcl_Interp *interp) {

    /* Allocate Array of arrays */
    if (NULL == (sequence_reg = ArrayCreate(sizeof(Array), 0)))
	return -1;

    if (NULL == (seq_cursor_reg = ArrayCreate(sizeof(Array), 0)))
	return -1;

    init_cursor_colour(interp);
    init_raster_colour(interp);
    return 0;
}

/*
 * Add a new sequence to the registration
 */
int
add_reg_seq(int index)
{
    
    /* return ref to memory and allocate more memory if necessary */
    (void)ArrayRef(sequence_reg, index);
    (void)ArrayRef(seq_cursor_reg, index);

    if (NULL == (seq_func_array(index) = ArrayCreate(sizeof(seq_reg), 0)))
	return -1;
    seq_Nfuncs(index) = 0;
    seq_cursor(index) = NULL;
    return 0;
}

/*
 * Delete a sequence from the registration
 */
void
delete_reg_seq(int index)
{
    seq_reg_plot jdata;

    /* remove registered functions */
    jdata.job = SEQ_DELETE;
    seq_notify(index, (seq_reg_data *)&jdata);

    
#ifdef DEBUG
    printf("delete_reg_seq \n");
    seq_register_dump();
#endif
    /* Deallocate */
    ArrayDestroy(seq_func_array(index));

    /* remove sequence from array structure */
    if (index < (ArrayMax(sequence_reg) - 1)) {
	memmove(arrp(Array, sequence_reg, index),
		arrp(Array, sequence_reg, index+1),
		(ArrayMax(sequence_reg) - index - 1) * sizeof(Array));
	memmove(arrp(Array, seq_cursor_reg, index),
		arrp(Array, seq_cursor_reg, index+1),
		(ArrayMax(seq_cursor_reg) - index - 1) * sizeof(Array));
    }
    ArrayMax(sequence_reg)--;
    ArrayMax(seq_cursor_reg)--;
#ifdef DEBUG
    printf("delete_reg_seq \n");
    seq_register_dump();
#endif
}

/*
 * return a unique identifier for each set of results
 */
int 
get_reg_id()
{
    static int id = 0;
    return(id++);
}

/*
 * Registers func(seq_num, fdata, jdata) with sequence 'seq_num'.
 * Doesn't check if the (func,fdata) pair are already existant.
 *
 * Returns 0 on success, and -1 for error. 
 */
int 
seq_register(int seq_num,
	     void (*func)(int seq_num, void *fdata, seq_reg_data *jdata),
	     void *fdata,
	     int type,
	     int id)
{
    seq_reg *r;
    int i;

    /* Check if this (func,fdata) pair already exists. */
    r = ArrayBase(seq_reg, seq_func_array(seq_num));
    for (i = 0; i < seq_Nfuncs(seq_num); i++) {
	if (r[i].func == func && r[i].fdata == fdata) {
	    return 0;
	}
    }

    /* Allocate new function item */
    if (NULL == (r = (seq_reg *)ArrayRef(seq_func_array(seq_num), 
					 seq_Nfuncs(seq_num))))
	return -1;

    /* Copy over our data */
    r->func = func;
    r->fdata = fdata;
    r->time = time(NULL);
    r->type = type;
    r->id = id;

    /* update results manager with new function */
    /* update_results(); */

    return 0;
}

/*
 * Deregisters func(seq_num, data, jdata) from sequence 'seq_num'.
 *
 * Returns 0 for success, and -1 for error.
 */
int 
seq_deregister(int seq_num,
	       void (*func)(int seq_num, void *fdata, seq_reg_data *jdata),
	       void *fdata) 
{
    seq_reg *r;
    int i;
    int num_funcs;

#ifdef DEBUG
    printf("DEREGISTER seq %d \n", seq_num);
    seq_register_dump();
#endif
    num_funcs = seq_Nfuncs(seq_num);
    r = ArrayBase(seq_reg, seq_func_array(seq_num));
    /* Search for element in the array */

    /* 
     * can have the same function registered twice with the same sequence eg
     * sip4 with the same sequence. Therefore must go through all the functions
     * registered and not stop when find the first one
     */
    for (i = 0; i < num_funcs; i++) {
	if (r[i].func == func && r[i].fdata == fdata) {
	    /* Match found - shuffle everything else down */
	    memmove(&r[i], &r[i+1], (seq_Nfuncs(seq_num) - i - 1) * sizeof(r[i]));
	    
	    /* Decrement count */
	    seq_Nfuncs(seq_num)--;

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
	memmove(&r[i], &r[i+1], (seq_Nfuncs(seq_num) - i - 1) * sizeof(r[i]));

	/* Decrement count */
	seq_Nfuncs(seq_num)--;
    }
#endif
    return 0;
}

/*
 * check to see if unique identifier id has already been registered
 */
int
is_seq_reg(int id) {
    int s, i, n;
    seq_reg *r;

    for (s = 0; s < ArrayMax(sequence_reg); s++) {
	n = seq_Nfuncs(s);
	for (i = 0; i < n; i++) {
	    r = &arr(seq_reg, seq_func_array(s), i);

	    if (r->id == id) {
		return 1;
	    }
	}
    }
    return 0;
}
/*
 * Uses the register list for a given seq_num to call a particular job.
 */
void 
seq_notifyOLD(int seq_num, 
           seq_reg_data *jdata) 
{
    seq_reg *r;
    int i;
    int num_funcs;
    
#ifdef DEBUG
    seq_register_dump();
#endif

    num_funcs = seq_Nfuncs(seq_num);
    r = ArrayBase(seq_reg, seq_func_array(seq_num));

    if (0 == num_funcs)
        return;

    /* 
     * Now loop through all funcs sending the notification
     * MUST go backwards to deal with case of deletion 
     */
    for (i = num_funcs - 1; i >= 0; i--) {
        /* Send the notification request itself */
        r[i].func(seq_num, r[i].fdata, jdata);
    }
}


/*
 * Uses the register list for a given seq_num to call a particular job.
 */
void seq_notify(int seq_num, 
		seq_reg_data *jdata) 
{
    seq_reg *r;
    int i, j, k;
    int num_funcs;
    int *ids;

#ifdef DEBUG
    seq_register_dump();
#endif

    num_funcs = seq_Nfuncs(seq_num);
    r = ArrayBase(seq_reg, seq_func_array(seq_num));

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
	num_funcs = seq_Nfuncs(seq_num);
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
	    if (i == num_funcs)
		continue;
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
seq_result_notify(int id, seq_reg_data *jdata, int all) {
    int i, j, k, s, n;
    int *ids;
    seq_reg *r;

    for (s = 0; s < ArrayMax(sequence_reg); s++) {

	n = seq_Nfuncs(s);
	r = ArrayBase(seq_reg, seq_func_array(s));

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
	    n = seq_Nfuncs(s);

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
	n = seq_Nfuncs(s);
	for (i = 0; i < n; i++) {
	    seq_reg *r = &arr(seq_reg, seq_func_array(s), i);
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
 * Uses the register list for all results to call a particular job.
 */
void 
seq_result_notify_all(seq_reg_data *jdata) {
    int i, s, n, j, k;
    int *uids;
    seq_reg *r;

    /* loop through sequences */
    for (s = 0; s < ArrayMax(sequence_reg); s++) {

	/* number of functions associated with each sequence */
	n = seq_Nfuncs(s);

	if (n == 0)
	    continue;

	r = ArrayBase(seq_reg, seq_func_array(s));
	if (NULL == (uids = (int *)xmalloc(n * sizeof(int))))
	    return;

	/*
	 * Take a snap shot of all the registrations that need notifying.
	 * We do this by remembering their unique uids.
	 */
	for (k = 0; k < n; k++) {
	    uids[k] = r[k].id;
	}

	/* Now loop through all uids sending the notification */
	for (j = i = 0; j < k; i++, j++) {
	    /* Don't optimise this outside the loop -
	       it may not be a constant */
	    
	    n = seq_Nfuncs(s);
	    /*
	     * It's possible that the next registration isn't our next uid.
	     * This occurs when the previous r[i].func() modified the 
	     * registration list for this contig (eg by deregistering 
	     * something). So we scan along to find the correct 'i' in this 
	     * case remembering that it may not even exist any more.
	     */
	    if (i >= n || r[i].id != uids[j]) {
		for (i = 0; i < n; i++) {
		    if (r[i].id == uids[j])
			break;
		}
		if (i == n)
		    continue;
	    }
	    /* Send the notification request itself */
	    r[i].func(s, r[i].fdata, jdata);
	}
	xfree(uids);
    }
   
}

/*
 * Uses the register list for all results of a given type to call a 
 * particular job.
 */
void seq_type_notify(int seq_num, int type, seq_reg_data *jdata) {
    int i;
    int num_funcs;
    seq_reg *r;

    num_funcs = seq_Nfuncs(seq_num);
    r = ArrayBase(seq_reg, seq_func_array(seq_num));

    if (0 == num_funcs)
	return;

    /* 
     * Now loop through all funcs sending the notification
     * MUST go backwards to deal with case of deletion 
     */
    for (i = num_funcs - 1; i >= 0; i--) {
	/* Send the notification request itself */
	if (r[i].type == type) {
	    r[i].func(seq_num, r[i].fdata, jdata);
	}
    }
}

/*
 * return the type of result
 */
int seq_get_type(int id) {
    int s, i, n;
    seq_reg *r;

    for (s = 0; s < ArrayMax(sequence_reg); s++) {
	n = seq_Nfuncs(s);
	for (i = 0; i < n; i++) {
	    r = &arr(seq_reg, seq_func_array(s), i);

	    if (r->id == id) {
		return (r->type);
	    }
	}
    }
    return -1;
}

/*
 * Generates description of functions registered with a particular sequence.
 * 'reg' to return the index into
 * the registration array for this sequence. This (contig,reg) pair specifies
 * a particular result without the need for remembering pointers.
 * 'id' contains a unique id number for this result.
 */
seq_reg_name * seq_result_names(int *num_elements) 
{
    seq_reg_query_name qn;
    static char buf[80];
    int n, i, j, k;
    seq_reg *r;
    int *id_array;
    int num_results;
    int found = 0;
    int num_ids = 0;
    seq_reg_name *data;
    int cnt = 0;

    num_results = seq_num_results();
    if (num_results == 0) {
	return NULL;
    }
    id_array = (int *)xmalloc(num_results * sizeof(int));
    data = (seq_reg_name *)xmalloc(num_results * sizeof(seq_reg_name));

    /* initialise id_array */
    for (i = 0; i < num_results; i++) {
	id_array[i] = -1;
	data[i].line = (char *)xmalloc(100 * sizeof(char));
	data[i].time = (char *)xmalloc(100 * sizeof(char));
    }

    qn.job = SEQ_QUERY_NAME;
    qn.line = buf;
    *qn.line = '\0';

    for (i = 0; i < ArrayMax(sequence_reg); i++) {
	n = seq_Nfuncs(i);

	for (j = 0; j < n; j++) {
	    r = &arr(seq_reg, seq_func_array(i), j);
	    found = 0;                            /* reset found to be false */

	    /* want to check if already done result */
	    for (k = 0; k < num_ids; k++) {
		if (id_array[k] == r->id){
		    found = 1;
		    break;
		}
	    }
	    /* if not in id_array
	     * add to results array 
	     */
	    if (!found) {
		r->func(i, r->fdata, (seq_reg_data *)&qn);
		strcpy(data[cnt].line, qn.line);
		data[cnt].id = r->id;
		strcpy(data[cnt].time, seq_result_time(i, r->id));
#ifdef DEBUG
		printf("line %s %s id %d time %s\n", qn.line, data[cnt].line, data[cnt].id, data[cnt].time);
#endif
		cnt++;
		id_array[num_ids++] = r->id;
	    }
	}
    }
    *num_elements = cnt;

    xfree(id_array);
    return data;
}

/*
 * Returns a time (string) that a given id was registered. Assumes all
 * contig registrations for a particular id are registered together.
 * 'contig' isn't really needed here, but it's currently known by tcl
 * and speeds up our search.
 */
char *seq_result_time(int seq_num, int id) {
    int i, n;
    seq_reg *r;
    static char buf[80];
    
    n = seq_Nfuncs(seq_num);
    r = ArrayBase(seq_reg, seq_func_array(seq_num));
    
    for (i = 0; i < n && r[i].id != id; i++)
	;

    if (i == n)
	return "unknown";

    /* %r doesn't work on Windows ! */
    strftime(buf, sizeof(buf)-1, "%a %I:%M:%S %p", localtime(&r[i].time));
    return buf;
}

/*
 * Dump lists
 */
void seq_register_dump(void) {
    seq_reg *r;
    int i, n, c;

    for (c = 0; c < ArrayMax(sequence_reg); c++) {

	printf("sequence %d\n", c);
	printf("num funcs!! %d \n", (int)seq_Nfuncs(c));
	n = seq_Nfuncs(c);
	r = ArrayBase(seq_reg, seq_func_array(c));

	for (i = 0; i < n; i++) 
	  printf("    Function 0x%p      Data 0x%p ID %d \n",
		   r[i].func, r[i].fdata, r[i].id);
    }
}

/*
 * returns the number of sequences registered
 */
int
seq_num_seqs()
{
    return ArrayMax(sequence_reg);
}

/*
 * total number of results registered
 */
int
seq_num_results()
{
    int s;
    int max_seqs;
    int num_funcs = 0;
  
    max_seqs = ArrayMax(sequence_reg);
    /* sequence array */
    for (s = 0; s < max_seqs; s++) {
	num_funcs += seq_Nfuncs(s);
    }
    return (num_funcs);
}

/*
 * creates an array of data that succeed the comparison function
 */
int
search_reg_data(int (*comparison)(void *fdata, int type),
		void **array,
		int *num_elements)
{
    seq_reg *r;
    int i, n, j, k;
    int cnt = 0;
    int *id_array;
    int num_results;
    int found = 0;
    int num_ids = 0;

#ifdef START_DEBUG  
    printf("start SEARCH-REG-DATA\n");
    seq_register_dump();
#endif

    num_results = seq_num_results();
    if (num_results == 0) {
	*num_elements = 0;
	return -1;
    }
    id_array = (int *)xmalloc(num_results * sizeof(int));
    for (i = 0; i < num_results; i++)
	id_array[i] = -1;

    for (i = 0; i < ArrayMax(sequence_reg); i++) {
	n = seq_Nfuncs(i);
	r = ArrayBase(seq_reg, seq_func_array(i));

	for (j = 0; j < n; j++) {
	    found = 0;
	    /* want to check if already done result */
	    for (k = 0; k < num_ids; k++) {
/*
		printf("i %d j %d k %d num_ids %d, id1 %d id2 %d\n", 
		       i, j, k, num_ids, id_array[k], r[j].id);
*/
		if (id_array[k] == r[j].id){
		    found = 1;
		    break;
		}
	    }
	    /* if not in id_array, and comparison is true, 
	     * add to results array 
	     */
/*
	    printf("found %d i %d j %d k %d num_ids %d id %d\n", 
		   found, i, j, k, num_ids, r[j].id);
*/
	    if (!found && comparison(r[j].fdata, r[j].type)) {
		array[cnt++] = r[j].fdata;
		id_array[num_ids++] = r[j].id;
	    }
	}
    }
    *num_elements = cnt;

#ifdef DEBUG
    for ( i = 0; i < cnt; i++) 
	printf("i %d id %d \n", i, array[i]);
#endif
    xfree(id_array);
    return 0;
}

/*
 * Converts a registration id to an array of seq_nums
 *
 * Returns NULL for failure or a NULL terminated 'int *' array which
 * is expected to be xfree()d by the calling function.
 */
int *result_to_seq_nums(int id, int *num_seqs) {
    int *rl;
    int i, j, n;
    int count = 0;
    seq_reg *r;

    if (NULL == (rl = (int *)xmalloc((ArrayMax(sequence_reg)+1) * 
				      sizeof(int))))
	return NULL;

    for (i = 0; i < ArrayMax(sequence_reg); i++) {
	n = seq_Nfuncs(i);

	for (j = 0; j < n; j++) {
	    r = &arr(seq_reg, seq_func_array(i), j);
	    if (r->id == id) {
		rl[count++] = i;
	    }
	}
    }

    *num_seqs = count;
    return rl;
}


/*
 * Returns the data component of a contig_reg_t for a specific id.
 *
 * We return only the first data for id found, or NULL if none found.
 */
void *result_data(int id, int seq_num) {
    int i, j, end, st, n;
    seq_reg *r;
    
    if (seq_num > -1) {
	end = seq_num;
	st  = seq_num;
    } else {
	st  = 1;
	end = ArrayMax(sequence_reg);
    }

    for (i = st; i <= end; i++) {
	n = seq_Nfuncs(i);
	for (j = 0; j < n; j++) {
	    r = &arr(seq_reg, seq_func_array(i), j);
	    if (r->id == id)
		return (r->fdata);
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
 * Returns the id for a registered item of a particular type, If sequence is
 * > -1, we search only this seq_num. Otherwise we scan all.
 *
 * We return only the first data for id found, or -1 if none found.
 */
int type_to_result(int type, int seq_num) {
    int i, j, end, st, n;
    seq_reg *r;
    
    if (seq_num > -1) {
	end = seq_num;
	st  = seq_num;
    } else {
	st  = 1;
	end = ArrayMax(sequence_reg);
    }

    for (i = st; i <= end; i++) {
	n = seq_Nfuncs(i);
	for (j = 0; j < n; j++) {
	    r = &arr(seq_reg, seq_func_array(i), j);
	    if (r->type == type)
		return r->id;
	}
    }

    return -1;
}

/*
 *-----------------------------------------------------------------------------
 * Cursor functions
 *-----------------------------------------------------------------------------
 */
static int cursor_id = 0;
static int *cursor_free_array = NULL;
static int size_free_array = 0;
static int num_cursors = 0;

int get_num_cursors(void)
{
    return num_cursors;
}

int get_cursor_id(void)
{
    int id;

    num_cursors++;
    if (size_free_array == 0) {
	return cursor_id++;
    }
    size_free_array--;
    id = cursor_free_array[0];
    memmove(&cursor_free_array[0], &cursor_free_array[1], 
	    size_free_array * sizeof(int));

    return id;
}

int add_cursor_free_array(int id) 
{
    static int array_size = 0;

    if (size_free_array >= array_size) {
	array_size += 10;
	if (NULL == (cursor_free_array = (int *)xrealloc(cursor_free_array, 
				      array_size * sizeof(int)))) {
	    xfree(cursor_free_array);
	    return -1;
	}
    }
    cursor_free_array[size_free_array++] = id;
    num_cursors--;
    return 0;
}

char *get_cursor_colour(int id)
{
    return(cursor_colour[id%NUM_CURSOR_COL]);
}

/*
 * Create a cursor for this sequence.
 * If private == 1 then create a new cursor if no non-private ones can be
 * found. Otherwise use an existing one if available.
 * Ie, private cursors can be used more than once, but only by one private
 * user at any one time.
 * Returns the cursor pointer, or NULL for failure.
 */
cursor_t *create_cursor(int seq_num, 
			int private, 
			char *colour,
			int line_width,
			int cursor_num,
			int direction)
{
    cursor_t *gc, *gcend;
    seq_cursor_notify cn;

    /* If private, look for a non private cursor to use */
    if (private) {
	for (gc = seq_cursor(seq_num); gc; gc = gc->next)
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
    if (!private && seq_cursor(seq_num)) {
	seq_cursor(seq_num)->refs++;
	gc = seq_cursor(seq_num);
	goto notify;
    }
#endif
    if (!private) {
#ifdef REMOVE
	for (gc = seq_cursor(seq_num); gc && --cursor_num > 0; gc = gc->next)
	    ;
#endif
	for (gc = seq_cursor(seq_num); gc; gc = gc->next) {
	    if (gc->direction == direction) {
		--cursor_num;
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
    /* gc->id = cursor_id++; */
    gc->id = get_cursor_id();
    if (gc->id >= NUM_CURSORS) {
	verror(ERR_WARN, "create cursor", "Too many cursors\n");
	return NULL;
    }

    gc->refs = 1;
    gc->abspos = 1;
    gc->private = private;
    gc->next = 0;
    
    /* allocate colour */
    if (!colour) {
	gc->colour = strdup(get_cursor_colour(gc->id));
    } else {
	gc->colour = strdup(colour);
    }
    gc->line_width = line_width;
    gc->direction = direction;

    /* Find last existing cursor */
    gcend = seq_cursor(seq_num);
    while(gcend && gcend->next)
	gcend = gcend->next;

    /* Add gc to it */
    if (gcend)
	gcend->next = gc;
    else
	seq_cursor(seq_num) = gc;

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
    seq_notify(seq_num, (seq_reg_data *)&cn);
    return gc;
}

int find_nearest_cursor(Tk_Raster *raster,
			int seq_num, 
			int pos,
			int max_dist,
			int direction,
			int *cursor_pos)                             /* out */
{
    cursor_t *gc;
    int closest = INT_MAX;
    int cursor_id = -1;
    int diff = INT_MAX;
    int rx, ry;
    double wx0, wy0, wx1, wy1;

#ifdef DEBUG
    printf("find_nearest_cursor %d\n", seq_num);
#endif

    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);
 
    for (gc = seq_cursor(seq_num); gc; gc = gc->next) {
	WorldToRaster(raster, gc->abspos, rasterY(raster, gc->abspos), &rx, &ry);

	if (direction == HORIZONTAL && gc->direction == HORIZONTAL) {
	    diff = abs(rx - pos);
#ifdef DEBUG
	    printf("H abspos %d rx %d pos %d diff %d closest %d id %d direct %d\n", 
		   gc->abspos, rx, pos, diff, closest, gc->id, gc->direction);
#endif
	} else if (direction == VERTICAL && gc->direction == VERTICAL) {
	    diff = abs(ry - pos);
#ifdef DEBUG
	    printf("V abspos %d ry %d pos %d diff %d closest %d id %d direct %d\n", 
		   gc->abspos, ry, pos, diff, closest, gc->id, gc->direction);
#endif
	} else {
	    diff = INT_MAX;
	}

	if (diff < closest) {
	    closest = diff;
	    cursor_id = gc->id;
	    if (direction == HORIZONTAL) {
		*cursor_pos = rx;
	    } else {
		*cursor_pos = ry;
	    }
	}
#ifdef DEBUG
	printf("closest %d max_dist %d\n", closest, max_dist);
#endif
    }
    if (closest > max_dist) {
	return -1;
    } else {
	return cursor_id;
    }
}

/*
 * Given a cursor identifier, return the cursor structure or NULL if not
 * found. If seq_num != NULL,only look in this seq_num. Otherwise look at all.
 * If found and seq_num != NULL, fill with the correct seq number.
 *
 * Always sends a notification about this cursor to inform of the new
 * reference count, and possible private status.
 */
cursor_t *find_cursor(int *seq_num, int id, int direction)
{
    cursor_t *gc;
    int s;

    if (seq_num && *seq_num != -1) {
	for (gc = seq_cursor(*seq_num); gc; gc = gc->next) {
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

    for (s = 0; s < ArrayMax(sequence_reg); s++) {
	if (seq_num)
	    *seq_num = s;
	for (gc = seq_cursor(s); gc; gc = gc->next) {
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

int seq_num_cursors(void)
{
    int s;
    int num = 0;
    cursor_t *gc;

    for (s = 0; s < ArrayMax(sequence_reg); s++) {
	for (gc = seq_cursor(s); gc; gc = gc->next) {
	    num++;
	}
    }

   return num;
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
void delete_cursor(int seq_num, int id, int private)
{
    cursor_t *gc, *gcp;
    int s = seq_num;
    seq_cursor_notify cn;
#ifdef DEBUG    
    printf("SEQ delete_cursor seq_num %d id %d private %d\n",
	   seq_num, id, private);
#endif
    if (!(gc = find_cursor(&s, id, -1)))
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
    seq_notify(s, (seq_reg_data *)&cn);

    if (gc->refs > 0)
	return;

    /* If first in seq_num, skip to next */
    if (seq_cursor(s) == gc) {
	seq_cursor(s) = gc->next;
	add_cursor_free_array(gc->id);
	free(gc->colour);
	xfree(gc);
	return;
    }

    /* Otherwise, find previous and link to next*/
    for (gcp = seq_cursor(s); gcp && gcp->next != gc; gcp = gcp->next);
    if (gcp && gcp->next == gc) {
	gcp->next = gc->next;
	add_cursor_free_array(gc->id);
	free(gc->colour);
	xfree(gc);
    }
    return;
}
