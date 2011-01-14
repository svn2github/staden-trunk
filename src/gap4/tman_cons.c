#include <ctype.h>
#include <tcl.h>

#include "tman_interface.h"
#include "edUtils.h"
#include <io_lib/Read.h>
#include "misc.h"
#include "gap_globals.h"
#include "tcl_utils.h"
#include <io_lib/traceType.h>
#include "tclXkeylist.h"
#include "tkTraceIO.h"

/* Standard number of samples per base in consensus trace. Must be even */
#define TWIDTH 12

typedef struct {
    Read *r;
    char *seq;
    int  *opos;
} diff_cons_seq;

typedef struct {
    signed int opos;
    double back;
    double forw;
} diff_cons_trace;

#define strand_matches(xx,seq,strand) \
    ((strand) == DB_Comp((xx), (seq)))

#if 0
static int cons_matches(EdStruct *xx, char *con, int start, int seq,
			int pos) {
    char *s;

    s = DB_Seq(xx, seq) + DB_Start(xx, seq) + pos - DB_RelPos(xx, seq);
    con = con + pos - start;

    if (toupper(s[0]) == con[0])
	return 1;
    else
	return 0;
}
#endif

#define cons_matches(xx,con,start,seq,pos) \
    ((toupper(*(DB_Seq((xx), (seq)) + DB_Start((xx), (seq)) + \
		(pos) - DB_RelPos((xx), (seq)))) == \
      *((con) + (pos) - (start))) ? 1 : 0)

/*
 * This function produces an array of trace sample positions for a given
 * sequence and range using its trace data and original position information.
 * Bases without original positions are assigned new ones by looking at the
 * positions of neighbouring bases.
 */
static int *get_trace_pos(Read *r, EdStruct *xx, int seq, int pos,
			  int from, int to, char *s, int do_comp) {
    int2 *op, *opos = NULL;
    int i, last_val = 0, last_pos = 0, next_pos, next_val;
    double avg_space = (r->basePos[r->NBases-1] - r->basePos[0]) /
	(double)(r->NBases - 1);
    int next_p;
    int *trace_pos;

    if (NULL == (trace_pos = (int *)xmalloc((to - from + 2) * sizeof(int))))
	return NULL;

    /* Copy integer DB_Opos to opos, complementing as required */
    if (seq != 0) {
	op = DB_Opos(xx, seq);
	if (DB_Comp(xx, seq) == COMPLEMENTED && do_comp) {
	    opos = (int2 *)xmalloc((to - from + 2) * sizeof(*opos));
	    if (NULL == opos) {
		xfree(trace_pos);
		return NULL;
	    }

	    for (i = 0; op[i] == 0 && i < DB_Length2(xx, seq); i++);
	    if (i != DB_Length2(xx, seq))
		last_val = op[i];
	    else
		last_val = DB_Length2(xx, seq);
	    for (i = from; i < to; i++) {
		if (op[i] != 0)
		    opos[i-from] = last_val - op[i] + 1;
		else
		    opos[i-from] = 0;
	    }
	    op = opos;
	} else {
	    op += from;
	}
    } else {
	/* Consensus sequences have no opos[] values, so we compute some */
	opos = (int2 *)xmalloc((to - from + 2) * sizeof(*opos));
	if (NULL == opos) {
	    xfree(trace_pos);
	    return NULL;
	}

	for (i = from; i <= to+1; i++) {
	    opos[i-from] = i+1-pos;
	}
	opos[to-from] = 0;
	opos[to-from+1] = 0;
	op = opos;
    }

    /*
     * Now expand 0s to midpoints between surrounding bases.
     * We have special cases for the boundary conditions.
     */
    if (op[0] == 0) {
	for (i = from; i <= to && op[i-from] == 0; i++);
	if (op[i-from]) {
	    last_pos = i-1;
	    last_val = op[i-from]-1;
	    if (last_val < 0)
		last_val = 0;
	} else {
	    last_pos = -1;
	    last_val = 0;
	}
    }
    for (i = from; i < to; i++) {
	/* A non zero position is fine, but be careful of the right end */
	if (op[i-from]) {
	    last_pos = i;
	    last_val = op[i-from]-1;

	    if (op[i-from] < r->NBases) {
		next_p = r->basePos[op[i-from]];
	    } else {
		next_p = r->basePos[op[i-from]-1] + avg_space;
		if (next_p > r->NPoints)
		    next_p = r->NPoints;
	    }

	    /*
	     * printf("%c %3d %3d %5d+%2d\n",
	     *        s[i], i, op[i-from]-1, r->basePos[op[i-from]-1],
	     *        next_p - r->basePos[op[i-from]-1] + 1);
	     */
	    trace_pos[i-from] = r->basePos[op[i-from]-1];

	} else {
	    /*
	     * For a zero position we have to scan rightwards to find the
	     * next base with a non zero original position. We then use
	     * this position to interpolate the positions of the bases
	     * inbetween.
	     */
	    double lval = last_pos >= 0 ? r->basePos[last_val] : 0;
	    double rval;

	    /* Scan right, checking for the end of the sequence */
	    for (next_pos = i;
		 next_pos < to && op[next_pos-from]==0;
		 next_pos++);
	    if (next_pos >= to) {
		next_pos = to-1;
		next_val = last_val + 1;
		rval = r->basePos[last_val+1];
	    } else {
		next_val = op[next_pos-from]-1;

		if (next_val < r->NBases) {
		    rval = r->basePos[next_val];
		} else {
		    rval = r->basePos[next_val-1] + avg_space;
		    if (rval > r->NPoints)
			rval = r->NPoints;
		}
	    }

	    /*
	     * So we now have lval and rval as the previous and next trace
	     * positions for bases that have original positions, and
	     * last_pos and next_pos specify the base numbers for the bases
	     * used to find lval and rval.
	     * So do some linear interpolation to get trace positions for
	     * all of the bases between last_pos and next_pos.
	     */
	    do {
		double tmp = next_pos - last_pos;
		if (tmp)
		    trace_pos[i-from] =
			(int)(lval + (i-last_pos)/tmp*(rval-lval) +.5);
		else
		    trace_pos[i-from] = lval;
	    } while (i+1 < next_pos && i+1 < to && i++);
	}
    }
    trace_pos[to-from] = trace_pos[to-from-1] + avg_space;
    if (trace_pos[to-from] >= r->NPoints)
	trace_pos[to-from] = r->NPoints-1;
#if 0
    /* Shift positions to be between bases */
    last = trace_pos[0] - avg_space;
    for (i = 0; i <= to-from; i++) {
	int tmp;

	tmp = trace_pos[i];
	trace_pos[i] = (trace_pos[i] + last) / 2;
	last = tmp;
    }
    if (trace_pos[0] < 0) trace_pos[0] = 0;
#endif

    if (opos) xfree(opos);

    return trace_pos;
}

static int cons_realloc_trace(Read *r, int *max_points, int new_points) {
    if (new_points > *max_points) {
	int old_max = *max_points;

	*max_points = new_points * 1.5;
	r->traceA = xrealloc(r->traceA, *max_points * sizeof(*r->traceA));
	r->traceC = xrealloc(r->traceC, *max_points * sizeof(*r->traceC));
	r->traceG = xrealloc(r->traceG, *max_points * sizeof(*r->traceG));
	r->traceT = xrealloc(r->traceT, *max_points * sizeof(*r->traceT));

	if (NULL == r->traceA ||
	    NULL == r->traceC ||
	    NULL == r->traceG ||
	    NULL == r->traceT)
	    return -1;

	memset(&r->traceA[old_max], 0,
	       (*max_points - old_max) * sizeof(*r->traceA));
	memset(&r->traceC[old_max], 0,
	       (*max_points - old_max) * sizeof(*r->traceC));
	memset(&r->traceG[old_max], 0,
	       (*max_points - old_max) * sizeof(*r->traceG));
	memset(&r->traceT[old_max], 0,
	       (*max_points - old_max) * sizeof(*r->traceT));
    }

    return 0;
}

static int do_empty_cons_base(EdStruct *xx, char *cons, int pos, int start,
			      Read *r, int off, int *max_points) {
    int i, width = TWIDTH, istart;

    if (-1 == cons_realloc_trace(r, max_points, off + width/2 + 1))
	return -1;

    if (off-width/2 < 0)
	istart = -(off-width/2);
    else
	istart = 0;

    for (i = istart; i < width; i++) {
	r->traceA[off+i-width/2] = 0;
	r->traceC[off+i-width/2] = 0;
	r->traceG[off+i-width/2] = 0;
	r->traceT[off+i-width/2] = 0;
    }
    r->base[pos-start] = cons[pos-start];
    r->basePos[pos-start] = off + width/2;
    r->prob_A[pos-start] = 0;
    r->prob_C[pos-start] = 0;
    r->prob_G[pos-start] = 0;
    r->prob_T[pos-start] = 0;

    return width;
}

static int do_cons_base(EdStruct *xx, char *cons, int pos, int start,
			int count, int *seqList, diff_cons_seq *rlist,
			Read *r, int off, int match, int *max_points) {
    int i, j, istart;
    int width;
    static diff_cons_trace *tr = NULL;
    static int diff_count = 0;
    int used_count;
    double avg_back, avg_wid;

    if (count == 0)
	return do_empty_cons_base(xx, cons, pos, start, r, off, max_points);

    if (count > diff_count) {
	tr = xrealloc(tr, count * 2 * sizeof(*tr));
	if (NULL == tr)
	    return -1;
	diff_count = count * 2;
    }

    /* printf("--- At position %3d ---", pos); */
    width = 0;
    used_count = 0;
    avg_back = 0;
    avg_wid = 0;
    for (i = 0; i < count; i++) {
	int seq = seqList[i];
	int p = pos-DB_RelPos(xx, seq);

	if (match && cons_matches(xx, cons, start, seq, pos) == 0) {
	    tr[i].opos = -1;
	    continue;
	}

	/*
	 * printf(" %d:%c",
	 *        DB_Number(xx, seq),
	 *        rlist[i].seq[p]);
	 */

	tr[i].opos = rlist[i].opos[p];
	if (p > 0)
	    tr[i].back = (double)(rlist[i].opos[p-1] + rlist[i].opos[p])/2.0;
	else
	    tr[i].back = tr[i].opos;
	if (p + 1 < DB_Length(xx, seq))
	    tr[i].forw  = (double)(rlist[i].opos[p] + rlist[i].opos[p+1])/2.0;
	else
	    tr[i].forw = tr[i].opos;
	width += tr[i].forw - tr[i].back;
	avg_back += tr[i].opos - tr[i].back;
	avg_wid += tr[i].forw - tr[i].back;
	used_count++;
    }
    /* putchar('\n'); */

    if (used_count == 0)
	return do_empty_cons_base(xx, cons, pos, start, r, off, max_points);

    /* width = ABS((int)((double)width / used_count + 0.5)); */
    /* width = 2 * (int)((double)width / used_count / 2 + 0.5); */

    /*
     * We have to set the width to a fixed amount to get the consensus trace
     * to be smooth. I'm certain that variable widths will work provided that
     * we iron out the rounding problems involved with this.
     */
    width = TWIDTH;

    if (-1 == cons_realloc_trace(r, max_points, off + width/2 + 1))
	return -1;

    if (off-width/2 < 0)
	istart = -(off-width/2);
    else
	istart = 0;
    for (i = istart; i < width; i++) {
	double scale_a = 0.0, scale_c = 0.0, scale_g = 0.0, scale_t = 0.0;
	double ratio;
	double posd, posm;
	int posi;
	int i1, i2;

	for (j = 0; j < count; j++) {
	    if (tr[j].opos == -1)
		continue;

	    if (DB_Comp(xx, seqList[j]) == UNCOMPLEMENTED) {
		ratio = (double)(tr[j].forw - tr[j].back) / width;
		posd = (double)i * ratio;
		posi = (signed int)posd;
		posm = posd - posi;

		i1 = (int)(tr[j].back+0.5) + posi;
		i2 = i1+1;
		scale_a += (rlist[j].r->traceA[i2] - rlist[j].r->traceA[i1])
		    * posm + rlist[j].r->traceA[i1];

		scale_c += (rlist[j].r->traceC[i2] - rlist[j].r->traceC[i1])
		    * posm + rlist[j].r->traceC[i1];

		scale_g += (rlist[j].r->traceG[i2] - rlist[j].r->traceG[i1])
		    * posm + rlist[j].r->traceG[i1];

		scale_t += (rlist[j].r->traceT[i2] - rlist[j].r->traceT[i1])
		    * posm + rlist[j].r->traceT[i1];
	    } else {
		ratio = (double)(tr[j].forw - tr[j].back) / width;
		posd = (double)i * ratio;
		posi = (signed int)posd;
		posm = posi - posd;

		i1 = (int)(tr[j].back+0.5) + posi;
		i2 = i1-1;
		scale_t += (rlist[j].r->traceA[i2] - rlist[j].r->traceA[i1])
		    * posm + rlist[j].r->traceA[i1];

		scale_g += (rlist[j].r->traceC[i2] - rlist[j].r->traceC[i1])
		    * posm + rlist[j].r->traceC[i1];

		scale_c += (rlist[j].r->traceG[i2] - rlist[j].r->traceG[i1])
		    * posm + rlist[j].r->traceG[i1];

		scale_a += (rlist[j].r->traceT[i2] - rlist[j].r->traceT[i1])
		    * posm + rlist[j].r->traceT[i1];
	    }
	}

	/* printf("i=%d a=%f c=%f g=%f t=%f\n",
	 *      i, scale_a, scale_c, scale_g, scale_t);
	 */

	r->traceA[off+i-width/2] = (TRACE)(scale_a / used_count);
	r->traceC[off+i-width/2] = (TRACE)(scale_c / used_count);
	r->traceG[off+i-width/2] = (TRACE)(scale_g / used_count);
	r->traceT[off+i-width/2] = (TRACE)(scale_t / used_count);
    }

    /*
     * Reposition base positions to the centre of peaks. This is vitally
     * important when combined with the trace_diff code.
     */
    {
	signed int tmp;
	if (avg_wid)
	    if (avg_back)
		tmp = off-width/2 + avg_back / avg_wid * width + 0.5;
	    else
		tmp = off;
	else
	    tmp = off-width/2 + 0.5;
	if (tmp < 0) tmp = 0;
	r->base[pos-start] = cons[pos-start];
	r->basePos[pos-start] = tmp;
	r->prob_A[pos-start] = 0;
	r->prob_C[pos-start] = 0;
	r->prob_G[pos-start] = 0;
	r->prob_T[pos-start] = 0;
    }

    return width;
}

static int tidy_up(Read *r, int nbases, int npoints) {
    int max;
    int i;

    r->NBases = nbases;
    r->NPoints = npoints;

    for (max = 1, i = 0; i < npoints; i++) {
	if (r->traceA[i] > max) max = r->traceA[i];
	if (r->traceC[i] > max) max = r->traceC[i];
	if (r->traceG[i] > max) max = r->traceG[i];
	if (r->traceT[i] > max) max = r->traceT[i];
    }

    r->maxTraceVal = max;

    return 0;
}

/* ------------------------------------------------------------------------- */

/*
 * Produce a consensus trace from a specific region of this contig.
 */
Read *cons_trace(EdStruct *xx, int start, int end, int strand,
		 int match, int exception) {
    int *seqList, i, j, count, next;
    Read *r;
    int max_points = 10000;
    char *con = NULL;
    diff_cons_seq *rlist = NULL;
    char fileName[256];
    char t_type[5];
    int form;
    int offset = 0, w;

    /* Get the consensus sequence */
    if (NULL == (con = (char *)xmalloc(end - start + 2)))
	goto error;
    DBcalcConsensus(xx, start, end - start + 1, con, NULL, BOTH_STRANDS);

    /* Allocate a list of read pointers and positions */
    if (NULL == (rlist = (diff_cons_seq *)xcalloc(DBI_gelCount(xx),
						  sizeof(*rlist))))
	goto error;

    /* Allocate a read structure */
    if (NULL == (r = read_allocate(max_points, end - start + 1)))
	goto error;

    /* Derive the initial list of sequences covering the start point */
    count = 0;
    seqList = DBI_list(xx);
    for (i = 1;
	 i <= DBI_gelCount(xx) && DB_RelPos(xx, DBI_order(xx)[i]) <= start;
	 i++) {
	int seq = DBI_order(xx)[i];
	DBgetSeq(DBI(xx), seq);
	if (DB_RelPos(xx, seq) + DB_Length(xx, seq) > start &&
	    strand_matches(xx, seq, strand) &&
	    seq != exception) {
	    if (get_trace_path(xx, seq, fileName, t_type) == 0) {
		form = trace_type_str2int(t_type);
		rlist[count].r = read_reading(fileName, form);
		if (rlist[count].r) {
		    rlist[count].seq = DBgetSeq(DBI(xx), seq);
		    rlist[count].opos =
			get_trace_pos(rlist[count].r, xx, seq, 0,
				      DB_Start(xx, seq),
				      DB_Start(xx, seq) + DB_Length(xx, seq),
				      DB_Seq(xx, seq), 0);

		    seqList[count++] = seq;
		}
	    }
	}
    }
    if (i <= DBI_gelCount(xx))
	next = i;
    else
	next = 0;

    /*
     * Loop along the sequence updating seqList as we go.
     * At each point we know how many sequences there are so we can
     * produce the consensus from these sequences.
     */
    for (i = start; i <= end; i++) {
	w = do_cons_base(xx, con, i, start, count, seqList, rlist, r, offset,
			 match, &max_points);
	if (w == -1)
	    goto error;
	offset += w;

	/* Update seqList for the next position */
	if (i < end) {
	    /* Remove sequences */
	    for (j = 0; j < count; j++) {
		int seq = seqList[j];
		if (DB_RelPos(xx, seq) + DB_Length(xx, seq) - 1 <= i) {
		    read_deallocate(rlist[j].r);
		    xfree(rlist[j].opos);
		    memmove(&seqList[j], &seqList[j+1],
			    (count-1-j) * sizeof(*seqList));
		    memmove(&rlist[j], &rlist[j+1],
			    (count-1-j) * sizeof(*rlist));
		    count--;
		    j--;
		}
	    }

	    /* Add sequences */
	    while (next && DB_RelPos(xx, next) <= i+1) {
		/* printf("next=%d %d %d\n",
		       next, DB_RelPos(xx, next), i+1); */
		DBgetSeq(DBI(xx), next);
		if (strand_matches(xx, next, strand) &&
		    get_trace_path(xx, next, fileName, t_type) == 0) {
		    form = trace_type_str2int(t_type);
		    rlist[count].r = read_reading(fileName, form);
		    if (rlist[count].r) {
			rlist[count].seq = DBgetSeq(DBI(xx), next);
			rlist[count].opos =
			    get_trace_pos(rlist[count].r, xx, next, 0,
					  DB_Start(xx, next),
					  DB_Start(xx,next)+DB_Length(xx,next),
					  DB_Seq(xx, next), 0);

			seqList[count++] = next;
		    }
		}
		if (++next > DBI_gelCount(xx))
		    next = 0;
	    }
	}
    }

    for (i = 0; i < count; i++) {
	read_deallocate(rlist[i].r);
	xfree(rlist[i].opos);
    }

    tidy_up(r, end-start + 1, offset);

    xfree(con);
    xfree(rlist);
    return r;

 error:
    if (con) xfree(con);
    if (rlist) xfree(rlist);
    return NULL;
}

/*
 * Add a consensus trace to the trace display.
 */
void cons_edc_trace(EdStruct *xx, int start, int end, int strand, int match,
		    int exception) {
    Read *r;
    char *pname;
    char buf[1024];
    Tcl_Interp *interp = EDINTERP(xx->ed);
    int exists;
    tman_dc *ed;
    DisplayContext *dc;
    static int cons_counter = 0;
    Tcl_CmdInfo info;
    int pos;

    /* Produce the read structure */
    if (NULL == (r = cons_trace(xx, start, end, strand, match, exception))) {
	bell();
	return;
    }

    /* Create a trace display */
    pname = get_default_string(interp, gap_defs, "TRACE_DISPLAY.WIN");
    Tcl_VarEval(interp, "trace_create ",
		Tk_PathName(EDTKWIN(xx->ed)), pname, " ",
		Tk_PathName(EDTKWIN(xx->ed)),
		" consensus", NULL);
    pname = Tcl_GetStringResult(interp);

    /* Fill out the tman_dc and DisplayContext structures */
    sprintf(buf, "Cons %d", cons_counter++);
    dc = getTDisplay(xx, buf, 0, 0, &exists);
    strcpy(dc->path, pname);
    ed = find_free_edc();
    ed->dc = dc;
    ed->pos = start-1;
    ed->xx = xx;
    ed->seq = 0;
    ed->type = TRACE_TYPE_CON;

    /* Add the Read to the trace widget */
    Tcl_GetCommandInfo(interp, Tcl_GetStringResult(interp), &info);
    trace_memory_load((DNATrace *)info.clientData, r);
    dc->tracePtr = (DNATrace *)info.clientData;

    /* Adjust position */
    Tcl_Eval(interp, "update idletasks");
    pos = positionInContig(xx, xx->cursorSeq, xx->cursorPos) - start;
    repositionSeq(xx, dc, pos);
}
