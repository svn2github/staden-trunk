#include <tcl.h>

#include "edUtils.h"
#include "tman_interface.h"
#include "Read.h"
#include "misc.h"
#include "gap_globals.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "tkTraceIO.h"
#include "mutlib.h"

static Read *diff_readings(EdStruct *xx,
			   Read *r1, int seq1, int off1,
			   Read *r2, int seq2, int off2,
			   int *startp, int *start1p) {
    int start1, end1, start2, end2, start, end;
    Read *r;
    tracediff_t td;

    /* One consensus trace works well, but not two. */
    if (!seq1 && !seq2)
	return NULL;

    /* Sequences must be in the same orientation */
    if (DB_Comp(xx, seq1) != DB_Comp(xx, seq2))
	return NULL;

    if (r1 == NULL || r2 == NULL)
	return NULL;

    /*
     * Compute start and end of the overlap point in each reading. This is
     * done using the non cutoff data only. We take 1 off the end position
     * to ensure that we always have a next base to find the position of.
     */

    /* This sets start and end to be positions in the consensus */
    if (xx->diff_trace_size) {
	start = positionInContig(xx, xx->cursorSeq, xx->cursorPos)
	    - xx->diff_trace_size;
	start1 = MAX(start, DB_RelPos(xx, seq1) - DB_Start(xx, seq1)-1);
	start2 = MAX(start, DB_RelPos(xx, seq2) - DB_Start(xx, seq2)-1);
	end = positionInContig(xx, xx->cursorSeq, xx->cursorPos)
	    + xx->diff_trace_size;
	end1 = MIN(end, DB_RelPos(xx, seq1) - DB_Start(xx, seq1)-1
		   + DB_Length2(xx, seq1) - 1);
	end2 = MIN(end, DB_RelPos(xx, seq2) - DB_Start(xx, seq2)-1
		   + DB_Length2(xx, seq2) - 1);
    } else {
	if (xx->reveal_cutoffs) {
	    start1 = DB_RelPos(xx, seq1) - DB_Start(xx, seq1)-1;
	    start2 = DB_RelPos(xx, seq2) - DB_Start(xx, seq2)-1;
	    end1   = DB_RelPos(xx, seq1) - DB_Start(xx, seq1)-1
		+ DB_Length2(xx, seq1) - 1;
	    end2   = DB_RelPos(xx, seq2) - DB_Start(xx, seq2)-1
		+ DB_Length2(xx, seq2) - 1;
	} else {
	    start1 = DB_RelPos(xx, seq1);
	    start2 = DB_RelPos(xx, seq2);
	    end1   = DB_RelPos(xx, seq1) + DB_Length(xx, seq1) - 1;
	    end2   = DB_RelPos(xx, seq2) + DB_Length(xx, seq2) - 1;
	}
    }
    start  = MAX(start1, start2);
    end    = MIN(end1, end2);
    start = MAX(start, 1);
    end   = MAX(end, 1);
    start = MIN(start, DB_Length(xx, 0));
    end   = MIN(end, DB_Length(xx, 0));
    if (end <= start)
	return NULL;

    /* Now we convert these to positions in the full-length editor sequences */
    start1 = start - (DB_RelPos(xx, seq1)-1) + DB_Start(xx, seq1);
    start2 = start - (DB_RelPos(xx, seq2)-1) + DB_Start(xx, seq2);
    end1   = end   - (DB_RelPos(xx, seq1)-1) + DB_Start(xx, seq1);
    end2   = end   - (DB_RelPos(xx, seq2)-1) + DB_Start(xx, seq2);

    /* Change from first/last used base to last/first clipped base */
    start1--;
    start2--;
    end1++;
    end2++;

    /*
     * And now we convert these to positions in the trace (orig orientation)
     * When comparing against a consensus trace we just use the entire lot
     * as this is how it's been generated in the first place (we hope).
     */
    if (seq1) {
	start1 = origpos(xx, seq1, start1);
	end1   = origpos(xx, seq1, end1);
    } else {
	end1  -= start1;
	start1 = 0;
    }
    if (seq2) {
	start2 = origpos(xx, seq2, start2);
	end2   = origpos(xx, seq2, end2);
    } else {
	end2  -= start2;
	start2 = 0;
    }


    /*
     * If complemented, change the start and end positions so that they map
     * on to a trace counting with base 1 at the left, as the tracediff
     * library does not have access to the trace display widget details (which
     * contains the real numbering order).
     */
    if (start1 > end1) {
	start1 = r1->NBases - start1 + 1;
	end1 = r1->NBases - end1 + 1;
    }
    if (start2 > end2) {
	start2 = r2->NBases - start2 + 1;
	end2 = r2->NBases - end2 + 1;
    }

    *startp = start;

    /* Initialise Mark's trace diff code */
    TraceDiffInit(&td);
    if( xx->compare_trace_yscale )
    	TraceDiffSetParameter( &td, TRACEDIFF_PARAMETER_YSCALE, 1.0 );
    TraceDiffSetReference(&td, r2, MUTLIB_STRAND_FORWARD, start2, end2);
    TraceDiffSetInput(&td, r1, MUTLIB_STRAND_FORWARD, start1, end1);

    /* Do the difference without analysis of the results */
    TraceDiffExecute(&td, TRACEDIFF_ALGORITHM_DEFAULT_DIFFERENCE_ONLY);
    if (TraceDiffGetResultCode(&td)) {
	verror(ERR_WARN, "diff_readings", "%s", TraceDiffGetResultString(&td));
	return NULL;
    }

    /* Get a copy of the result and then destroy the TraceDiff instance */
    r = TraceDiffGetDifference(&td,start1p, NULL);

    if (!seq1) {
	*start1p += start2-1;
    }

    if (r) {
	r = read_dup(r, NULL);
	/* set baseline and maxTraceVal */
	/* diff_reset_zero(r); */
    }
    TraceDiffDestroy(&td);

    return r;
}

/* ------------------------------------------------------------------------- */

/*
 * Given two already loaded traces (ed1 and ed2), produce a new trace display
 * containing the traces differences over the region that the traces overlap.
 */
DisplayContext *diff_edc_traces(EdStruct *xx, tman_dc *ed1, tman_dc *ed2) {
    Tcl_CmdInfo info;
    Read *r1, *r2, *r;
    char *pname;
    Tcl_Interp *interp = EDINTERP(xx->ed);
    tman_dc *ed;
    DisplayContext *dc;
    static int diff_counter = 0;
    int exists;
    char buf[1024], name[1024];
    int start, start1;

    /* Get the two read structures */
    Tcl_GetCommandInfo(interp, ed1->dc->path, &info);
    r1 = ((DNATrace *)info.clientData)->read;

    Tcl_GetCommandInfo(interp, ed2->dc->path, &info);
    r2 = ((DNATrace *)info.clientData)->read;

    /* Produce a diff read structure */
    r = diff_readings(xx, r1, ed1->seq, ed1->pos, r2, ed2->seq, ed2->pos,
		      &start, &start1);
    if (r == NULL) {
	bell();
	return NULL;
    }

    /* Create a trace display */
    pname = get_default_string(interp, gap_defs, "TRACE_DISPLAY.WIN");
    if (ed1->seq)
	sprintf(name, " {diffs: #%d #%d}",
		DB_Number(xx, ed1->seq), DB_Number(xx, ed2->seq));
    else
	sprintf(name, " {diffs: =%d #%d}",
		-DB_Number(xx, ed1->seq), DB_Number(xx, ed2->seq));
    Tcl_VarEval(interp, "trace_create ",
		Tk_PathName(EDTKWIN(xx->ed)), pname, " ",
		Tk_PathName(EDTKWIN(xx->ed)),
		name, NULL);
    pname = interp->result;

    /* Fill out the tman_dc and DisplayContext structures */
    sprintf(buf, "Diffs %d", diff_counter++);
    dc = getTDisplay(xx, buf, 0, 0, &exists);
    strcpy(dc->path, pname);
    ed = find_free_edc();
    ed->dc = dc;
    ed->pos = start-1;
    ed->xx = xx;
    ed->type = TRACE_TYPE_DIFF;
    ed->derivative_seq = ed1->seq ? ed1->seq : ed2->seq;
    ed->derivative_offset = start1;

    /* Add the Read to the trace widget */
    Tcl_GetCommandInfo(interp, interp->result, &info);
    trace_memory_load((DNATrace *)info.clientData, r);
    dc->tracePtr = (DNATrace *)info.clientData;

    /* Adjust position */
    {
	int num, end;
	int pos = positionInContig(xx, xx->cursorSeq, xx->cursorPos);

	num = tman_get_trace_position(xx, ed, pos, &end);
	repositionSeq(xx, dc, num);
    }

    return dc;
}
