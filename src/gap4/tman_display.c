#include <tk.h>

#include "misc.h"
#include "edStructs.h"
#include "tman_display.h"
#include "tman_interface.h"
#include "gap_globals.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"

static int init_done = 0;
static DisplayContext contexts[MAXCONTEXTS];
static int context_list[MAXCONTEXTS];

/*
 * Converts a trace path name to a display context */
DisplayContext *trace_path_to_dc(char *path) {
    int i;

    for (i = 0; i < MAXCONTEXTS; i++) {
	if (context_list[i] >= 0 && strncmp(contexts[context_list[i]].path,
					    path, 1024) == 0) {
	    return &contexts[context_list[i]];
	}
    }

    return NULL;
}

/*
 * 'Flash / Flick' the trace. Added to counteract the useful feature that the
 * earlier buggy drawCursor had (it redrew half the trace!) - easy to spot the
 * trace being redrawn
 */
static void flashCursor(EdStruct *xx, DisplayContext *dc) {
    Tcl_VarEval(EDINTERP(xx->ed), "trace_highlight [winfo parent ",
		dc->path, "] 1", NULL);
}

void deleteTraceDisplay(EdStruct *xx, DisplayContext *dc) {
    char buf[1024];
    tman_dc *edc;
    int i, num = -1;
    int mini_trace;

    if (!dc)
	return;

    for (i = 0; i < MAXCONTEXTS; i++) {
	if (context_list[i] >= 0 && &contexts[context_list[i]] == dc) {
	    num = i;
	    break;
	}
    }
    
    mini_trace = dc->mini_trace;

    /* Remove num and shuffle remaining items down */
    if ((edc = find_edc(dc)) && !mini_trace)
	tman_unhighlight(edc);

    dc->used = 0;
    strcpy(buf, dc->path);

    /*
     * This order is important. If we destroy the widget before removing
     * it from the context list then it will be destroyed twice, as there is
     * a <Destroy> binding on the widget would call (eventually) this code
     * again.
     */
    if (num < MAXCONTEXTS-1) {
	memmove(&context_list[num], &context_list[num+1],
		sizeof(int) * (MAXCONTEXTS-1 - num));
    }
    context_list[MAXCONTEXTS-1] = -1;

    if (mini_trace) {
	/* Mini traces are just a dnatrace widget */
	Tcl_VarEval(EDINTERP(xx->ed), "destroy ", buf, NULL);
    } else {
	/*
	 * Full traces are a complex of windows, of with dc->path is a child.
	 * So we destroy the parent instead.
	 */
	Tcl_VarEval(EDINTERP(xx->ed), "dnatrace_remove ", buf, NULL);
    }
}

void freeTDisplay(char *file) {
    int i;

    /* Search for context with this path and return if found */
    for (i = 0; i < MAXCONTEXTS; i++) {
	if (context_list[i] >= 0 && strncmp(contexts[context_list[i]].file,
					    file, FILE_NAME_LENGTH) == 0) {
	    if (i != MAXCONTEXTS-1) {
		memmove(&context_list[i], &context_list[i+1],
			sizeof(int) * (MAXCONTEXTS-1 - i));
	    }
	    context_list[MAXCONTEXTS-1] = -1;
	    return;
	}
    }
}

DisplayContext *getTDisplay(EdStruct *xx, char *file, int allow_dup,
			    int mini_trace, int *exists) {
    int i, fre = -1;
    DisplayContext *dc;

    /* Initialise if required */
    if (!init_done) {
	for (i = 0; i < MAXCONTEXTS; i++) {
	    context_list[i] = -1;
	    contexts[i].used = 0;
	    contexts[i].complemented = 0;
	    contexts[i].mini_trace = 0;
	}
	init_done = 1;
    }

    if (!allow_dup && !mini_trace) {
	/* Search for context with this path and return if found */
	for (i = 0; i < MAXCONTEXTS; i++) {
	    if (context_list[i] >= 0 &&
		strncmp(contexts[context_list[i]].file,
			file, FILE_NAME_LENGTH) == 0 &&
		contexts[context_list[i]].mini_trace == mini_trace) {
		*exists = 1;
		return &contexts[context_list[i]];
	    }
	}
    }

    *exists = 0;
    /* Make a free context */
    for (i = 0; i < MAXCONTEXTS; i++) {
	if (context_list[i] == -1)
	    break;
    }

    if (i == MAXCONTEXTS) {
	/*
	 * None free, so remove the first.
	 * FIXME: we could change this by having 'timestamps' on the 'used'
	 * element and searching for the oldest one. Does this really matter
	 * though?
	 */
	deleteTraceDisplay(xx, &contexts[context_list[0]]);
	i--;
    }
    fre = i;

    /* 'fre' now holds the context number allocated */
    for (i = 0; i < MAXCONTEXTS; i++)
	if (!contexts[i].used)
	    break;
    context_list[fre] = i;
    dc = &contexts[i];
    strncpy(dc->file, file, FILE_NAME_LENGTH);
    dc->used = 1;
    dc->complemented = 0;
    dc->mini_trace = mini_trace;

    return dc;
}

DisplayContext *manageTrace(EdStruct *xx,
			    char *format,
			    char *rawDataFile,
			    int baseNum,
			    int leftCutOff,
			    int cutLength,
			    int complemented,
			    int baseSpacing,
			    char *traceTitle,
			    int allow_dup,
			    int small_seq
			    ) {
    char *traceName;
    DisplayContext *dc;
    int exists;
    Tcl_Interp *interp = EDINTERP(xx->ed);
    char buf[1024];
    char *pname;
    Tcl_CmdInfo info;
    char *edpath;
    char seqbuf[1024];

    if ((traceName=(char *)strrchr(rawDataFile,'/'))==NULL)
        traceName = rawDataFile;
    else
        traceName++;

    dc = getTDisplay(xx, traceName, allow_dup, small_seq, &exists);
    if (exists) {
	repositionSeq(xx, dc, baseNum);
	flashCursor(xx, dc);
	return dc;
    }

    pname = get_default_string(interp, gap_defs, "TRACE_DISPLAY.WIN");

    /*
     * If we're the bottom half of a join editor, combine traces with the
     * top half.
     */
    if (inJoinMode(xx) && xx->link && xx == xx->link->xx[1] && !small_seq) {
	edpath = Tk_PathName(EDTKWIN(xx->link->xx[0]->ed));
    } else {
	edpath = Tk_PathName(EDTKWIN(xx->ed));
    }

    if (small_seq) {
	/* Mini-traces embedded in the editor */
	sprintf(seqbuf, "%d %d", small_seq, xx->lines_per_seq-1);
	if (TCL_OK != Tcl_VarEval(interp, "trace_small_add ",
				  edpath, pname, " {", rawDataFile, "} {",
				  edpath, "} ", seqbuf, NULL)) {
	    freeTDisplay(traceName);
	    puts(Tcl_GetStringResult(interp));
	    return NULL;
	}
    } else {
	/* The full-blown trace display. */
	if (TCL_OK != Tcl_VarEval(interp, "trace_add ",
				  edpath, pname, " {", rawDataFile, "} {",
				  edpath, "} {", traceTitle, "}", NULL)) {
	    freeTDisplay(traceName);
	    return NULL;
	}
    }
    strcpy(dc->path, Tcl_GetStringResult(interp));

    /* Get Trace widget pointer */
    if (-1 == Tcl_GetCommandInfo(interp, Tcl_GetStringResult(interp), &info)) {
	freeTDisplay(traceName);
	return NULL;
    }
    dc->tracePtr = (DNATrace *)info.clientData;

    /* Set orientation and cutoffs */
    if (complemented)
	Tcl_VarEval(interp, dc->path, " complement", NULL);
    dc->complemented = complemented;

    if (complemented) {
	leftCutOff = dc->tracePtr->Ned - (leftCutOff - 1);
	cutLength = 2-cutLength;
    }

    sprintf(buf, "%s left_cutoff %d", dc->path, leftCutOff);
    Tcl_Eval(interp, buf);

    sprintf(buf, "%s right_cutoff %d", dc->path, leftCutOff + cutLength);
    Tcl_Eval(interp, buf);

    /* Tcl_VarEval(interp, "update idletasks", NULL); */

    /* Adjust position */
    repositionSeq(xx, dc, baseNum);

    return dc;
}

void repositionSeq(EdStruct *xx, DisplayContext *dc, int baseNum) {
    char buf[1024];
    int offset;
    double f1;

/* New style syntax */
#if 1
    /* Convert base number into centred sample number */
    if (dc->tracePtr->comp)
	baseNum = dc->tracePtr->Ned - baseNum - 1;
    offset = trace_get_pos(dc->tracePtr, baseNum) -
	TKSHEET(xx->ed)->sw.font_width / 2;
    offset -= dc->tracePtr->disp_width / 2;
    f1 = offset / (double)dc->tracePtr->read->NPoints;

    sprintf(buf, "%s xview moveto %g;%s icursor %d\n",
	    dc->path, f1,
	    dc->path, baseNum);

#else
    sprintf(buf, "%s xview C%%%d;%s icursor %d\n",
	    dc->path, baseNum,
	    dc->path, baseNum);
#endif
    Tcl_Eval(EDINTERP(xx->ed), buf);
}

void deleteTrace(EdStruct *xx, char *path) {
    deleteTraceDisplay(xx, trace_path_to_dc(path));
}

void diffTrace(EdStruct *xx, char *path1, char *path2) {
    signed int t1 = -1, t2 = -1;
    tman_dc *edc1, *edc2;
    int i;

    /*
     * Ensure that we've got a blank trace slot available. Not doing this may
     * cause the traces to be shuffled up when we bring up the differences,
     * which due to the use of edc pointers causes confusion.
     */
    for (i = 0; i < MAXCONTEXTS; i++) {
	if (context_list[i] == -1)
	    break;
    }
    if (i == MAXCONTEXTS)
	deleteTraceDisplay(xx, &contexts[context_list[0]]);

    for (i = 0; i < MAXCONTEXTS; i++) {
	if (context_list[i] >= 0 &&
	    strncmp(contexts[context_list[i]].path, path1, 1024) == 0) {
	    t1 = i;
	    if (t2 != -1)
		break;
	}

	if (context_list[i] >= 0 &&
	    strncmp(contexts[context_list[i]].path, path2, 1024) == 0) {
	    t2 = i;
	    if (t1 != -1)
		break;
	}
    }

    if (t1 == -1 || t2 == -1 || t1 == t2) {
	bell();
	return;
    }

    if (NULL == (edc1 = find_edc(&contexts[context_list[t1]])) ||
	NULL == (edc2 = find_edc(&contexts[context_list[t2]]))) {
	bell();
	return;
    }

    if (DBI_contigNum(edc1->xx) != DBI_contigNum(edc2->xx)) {
	bell();
	return;
    }

    diff_edc_traces(xx, edc1, edc2);
}
