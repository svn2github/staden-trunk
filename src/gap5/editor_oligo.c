/*
 * The "select_oligos" function of the contig editor.
 *
 * FIXME: move this outside of the editor and use it as a general
 * command operating on a contig region? (Or contigs - for PCR.)
 */

#include <tcl.h>

#include "editor_view.h"
#include "primlib.h"
#include "dna_utils.h"
#include "gap4_compat.h" /* bell */

/*
 *----------------------------------------------------------------------------
 * Publically callable functions below here
 *----------------------------------------------------------------------------
 */

/*
 * Find suitable oligos using Primer3 with current parameter settings.
 * Return a Tcl list object on success.
 *          NULL on failure.
number of oligos found (or -1 for error).
 */
Tcl_Obj *edSelectOligoGenerate(edview *xx, int is_fwds, int bkwd_width,
			       int fwd_width, int avg_read_len,
			       char *primer_defs) {
    Tcl_Obj *lobj;
    int apos = xx->cursor_apos;

    int i, j;
    primlib_state *state;
    primlib_args *args;

    char *consensus;
    int *opos;
    int left, right;
    int cleft, cright;
    int consensusLength;


    /* Create an initialise primlib state */
    state = primlib_create();
    if (NULL == (args = primlib_str2args(primer_defs)))
	return NULL;

    primlib_set_args(state, args);
    free(args);

    /*
     * Conceptually we select off consensus.
     * Determine consensus around point...
     *    the oligo will be selected from this region
     */
    if (is_fwds) {
	left  = apos - bkwd_width;
	right = apos + fwd_width;
    } else {
	left  = apos - fwd_width;
	right = apos + bkwd_width;
    }
    if (0 == consensus_valid_range(xx->io, xx->cnum, &cleft, &cright)) {
	if (left < cleft)
	    left = cleft;
	if (right > cright)
	    right = cright;
    } else {
	contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);

	if (left  < c->start)
	    left  = c->start;
	if (right > c->end)
	    right = c->end;
    }

    /*
     * Allocated and calculate the consensus sequence.
     */
    consensusLength = right - left + 1;
    if (NULL == (consensus = (char *)xmalloc(consensusLength+1)) ||
	NULL == (opos = (int *)xmalloc((consensusLength+1)*sizeof(int)))) {
	bell();
	return NULL;
    }
    calculate_consensus_simple(xx->io, xx->cnum, left, right, consensus, NULL);
    consensus[consensusLength] = 0;

    /*
     * Complement if necessary
     */
    if (!is_fwds) {
	complement_seq(consensus, consensusLength);
    }

    /* Depad it */
    for (j = i = 0; i < consensusLength; i++) {
	opos[i] = j;
	if (consensus[i] != '*') {
	    consensus[j] = consensus[i];
	    j++;
	}
    }
    consensus[j] = 0;

    /* Pick the actual oligos */
    if (-1 == primlib_choose(state, consensus) || state->nprimers == 0) {
	xfree(opos);
	xfree(consensus);
	primlib_destroy(state);

	return NULL;
    }


    /*
     * Generate a Tcl list of the results.
     *
     * This consists of start, end, seq, quality, gc_content, melting temp.
     */
    lobj = Tcl_NewListObj(0, NULL);
    for (i = 0; i < state->nprimers; i++) {
	Tcl_Obj *l = Tcl_NewListObj(0, NULL);
	int st = state->primers[i].start, st2 = st;
	int en = st + state->primers[i].length-1, en2 = en;
	for (j = st; j < consensusLength; j++) {
	    if (is_fwds) {
		if (opos[j] == st)
		    st2 = j;
		if (opos[j] == en)
		    en2 = j;
	    } else {
		if (opos[j] == st)
		    en2 = consensusLength - j - 1;
		if (opos[j] == en)
		    st2 = consensusLength - j - 1;
	    }
	}

	Tcl_ListObjAppendElement(xx->interp, l, 
				 Tcl_NewStringObj("start", -1));
	Tcl_ListObjAppendElement(xx->interp, l,
				 Tcl_NewIntObj(st2+left));

	Tcl_ListObjAppendElement(xx->interp, l, 
				 Tcl_NewStringObj("end", -1));
	Tcl_ListObjAppendElement(xx->interp, l,
				 Tcl_NewIntObj(en2+left));

	Tcl_ListObjAppendElement(xx->interp, l, 
				 Tcl_NewStringObj("sequence", -1));
	Tcl_ListObjAppendElement(xx->interp, l,
				 Tcl_NewStringObj(&consensus[st], en-st+1));

	Tcl_ListObjAppendElement(xx->interp, l, 
				 Tcl_NewStringObj("quality", -1));
	Tcl_ListObjAppendElement(xx->interp, l,
				 Tcl_NewDoubleObj(state->primers[i].quality));

	Tcl_ListObjAppendElement(xx->interp, l, 
				 Tcl_NewStringObj("GC", -1));
	Tcl_ListObjAppendElement(xx->interp, l,
				 Tcl_NewDoubleObj(state->primers[i].gc_content));

	Tcl_ListObjAppendElement(xx->interp, l, 
				 Tcl_NewStringObj("temperature", -1));
	Tcl_ListObjAppendElement(xx->interp, l,
				 Tcl_NewDoubleObj(((int)(state->primers[i].temp * 100))/100.0));

	Tcl_ListObjAppendElement(xx->interp, l, 
				 Tcl_NewStringObj("end_stability", -1));
	Tcl_ListObjAppendElement(xx->interp, l,
				 Tcl_NewDoubleObj(state->primers[i].end_stability));

	Tcl_ListObjAppendElement(xx->interp, l, 
				 Tcl_NewStringObj("self_any", -1));
	Tcl_ListObjAppendElement(xx->interp, l,
				 Tcl_NewDoubleObj(state->primers[i].self_any / 100.0));

	Tcl_ListObjAppendElement(xx->interp, l, 
				 Tcl_NewStringObj("self_end", -1));
	Tcl_ListObjAppendElement(xx->interp, l,
				 Tcl_NewDoubleObj(state->primers[i].self_end / 100.0));

	Tcl_ListObjAppendElement(xx->interp, l, 
				 Tcl_NewStringObj("self_end", -1));
	Tcl_ListObjAppendElement(xx->interp, l,
				 Tcl_NewDoubleObj(state->primers[i].self_end / 100.0));

	Tcl_ListObjAppendElement(xx->interp, lobj, l);
    }

    xfree(opos);
    xfree(consensus);
    primlib_destroy(state);

    return lobj;
}

