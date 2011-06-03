/*
 * File: tman_interface.c
 * Version: 2.0
 *
 * This version interfaces between the contig editor trace requests
 * (edUtils) and the trace display mechanism (currently tcl/tk).
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: An standard interface to trace manipulation (except for the
 * manageTrace() call in tman_main.c. Designed to allow interaction
 * between the contig editor and the tman module.
 *
 * This should be the only code that understands both the contig editor and the
 * trace manager details.
 *
 * Created: 05 April 1994 (1.0)
 * Updated: 08 Feb   1995 (2.0)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <tk.h>

#include "edUtils.h"
#include "tman_interface.h"
#include "tman_display.h"
#include "contigEditor.h"
#include "edStructs.h"
#include <io_lib/Read.h>
#include "os.h"
#include "misc.h"
#include "array.h"
#include "dstring.h"

static tman_dc edc[MAXCONTEXTS]; /* initialise to null? */

/*
 * Look for first edc element with a NULL dc field (ie not used)
 */
tman_dc *find_free_edc(void) {
    int i;

    for (i = 0; i < MAXCONTEXTS && edc[i].dc != NULL; i++)
	;

    if (i == MAXCONTEXTS) {
	fprintf(stderr, "WARNING - Reusing an old trace! "
		"This should never happen.\n");
	i = 0; /* cheat - shouldn't happen */
    }

    edc[i].derivative_seq = 0;
    edc[i].derivative_offset = 0;

    return &edc[i];
}

/*
 * ---------------------------------------------------------------------------
 * Public callable functions.
 */

/*
 * Search for dc in edc array.
 * If we find it then we know which element this dc corresponds to.
 */
tman_dc *find_edc(DisplayContext *dc) {
    int i;

    for (i = 0; i < MAXCONTEXTS && edc[i].dc != dc; i++)
	;

    return (i == MAXCONTEXTS) ? NULL : &edc[i];
}


/*
 * Brings up a trace window.
 * Will highlight the appropriate reading.
 * Also unhighlights a reading if necessary (ie that we already had MAXCONTEXTS
 * traces displayed.
 */
DisplayContext *tman_manage_trace(
				  char *format,
				  char *rawDataFile,
				  int baseNum,
				  int leftCutOff,
				  int cutLength,
				  int complemented,
				  int baseSpacing,
				  char *traceTitle,
				  EdStruct *xx,
				  int seq,
				  int allow_dup,
				  int mini_trace
				  ) {
    tman_dc *ed;
    DisplayContext *dc;

    dc = manageTrace(xx, format, rawDataFile, baseNum, leftCutOff, cutLength,
		     complemented, baseSpacing, traceTitle, allow_dup,
		     mini_trace ? seq : 0);

    if (!dc)
	return NULL;

    /*
     * Is this a reuse of an existant context? If so then we need to notify
     * the editor of its removal.
     */
    if (ed = find_edc(dc)) {
	tman_unhighlight(ed);
    } else {
	ed = find_free_edc();
    }

    ed->dc = dc;
    ed->seq = seq;
    ed->pos = 0;
    ed->type = mini_trace ? TRACE_TYPE_MINI : TRACE_TYPE_SEQ;
    ed->xx = xx;
    ed->derivative_seq = 0;
    ed->derivative_offset = 0;

    if (!mini_trace)
	tman_highlight(ed);

    return dc;
}

/*
 * Converts editor consensus position to trace coordinate (bases)
 */
int tman_get_trace_position(EdStruct *xx, tman_dc *dc, int pos, int *end) {
    int seq, p, num;

    seq = dc->derivative_seq ? dc->derivative_seq : dc->seq;
    
    if (!dc->derivative_seq)
	(void) DBgetSeq(DBI(xx), seq); /* force sequence to be read */

    if (DB_Length2(xx, seq) == 0)
	return 0;

    p = pos - DB_RelPos(xx, seq) + DB_Start(xx, seq) + 1;
    if (p < 1) {
	/* best guess */
	return tman_get_trace_position(xx, dc, 1 + DB_RelPos(xx,seq)
				       - DB_Start(xx, seq) - 1, end) -
	    (DB_Comp(xx, seq) != COMPLEMENTED ? 1-p : p-1);
    } else if (p > DB_Length2(xx, seq)) {
	/* best guess */
	int diff = p - DB_Length2(xx, seq);
	pos = DB_Length2(xx, seq) + DB_RelPos(xx, seq) - DB_Start(xx, seq) - 1;
	return tman_get_trace_position(xx, dc, pos, end) +
	    (DB_Comp(xx, seq) != COMPLEMENTED ? diff : -diff);
    }
    
    num = origpos(xx, seq, p) - 1;

    if (dc->derivative_offset && dc->derivative_seq) {
	if (DB_Comp(xx, seq) == COMPLEMENTED) {
	    /* Just trust me on this one! */
	    num = DB_Length2(xx, seq) - num -2;
	}
	num -= dc->derivative_offset;
    }
    
    if (end)
	*end = DB_Length2(xx, seq);

    return num;
}

/*
 * Repositions displayed traces to be centred upon position 'pos'.
 * 'pos' is a contig offset.
 */
void tman_reposition_traces(EdStruct *xx, int pos, int mini_trace) {
    int i, num, end;

    if (!xx->trace_lock)
	return;

    for (i = 0; i < MAXCONTEXTS; i++) {
	if (edc[i].dc) {
	    EdStruct *yy = edc[i].xx;

	    switch (edc[i].type) {
	    case TRACE_TYPE_MINI:
		if (xx != yy || !mini_trace)
		    continue;

		num = tman_get_trace_position(xx, &edc[i], pos, &end);
		break;

	    case TRACE_TYPE_SEQ:
	    case TRACE_TYPE_POS_CONTROL:
	    case TRACE_TYPE_NEG_CONTROL:
	    case TRACE_TYPE_DIFF:
		if (xx != yy || mini_trace)
		    continue;

		num = tman_get_trace_position(xx, &edc[i], pos, &end);
		break;

	    case TRACE_TYPE_CON:
		num = pos - edc[i].pos - 1;
		end = 999999;
		break;

	    default:
		continue;
	    }


#if 0
	    /*
	     * Remove sequence? if it is out of range
	     */
	    if (num < 0 || num >= end) {
		continue;
	    }
#endif

	    repositionSeq(xx, edc[i].dc, num);
	}
    }
}

/*
 * Notifies the contig editor of the removal of a trace.
 * This simply unhighlights the relevent line on the sequence name display.
 */
void tman_unhighlight(tman_dc *edc) {
    EdStruct *xx = edc->xx;
/*
    if (!xx || xx->editorState == StateDown)
	return;
*/
    DBsetFlags(xx, edc->seq, DB_Flags(xx, edc->seq) & ~DB_FLAG_TRACE_SHOWN);

    edc->dc = NULL;

    RedisplayName(xx, edc->seq);
    redisplaySequences(xx, 1);
}

/*
 * Highlights a trace in the contig editor.
 */
void tman_highlight(tman_dc *edc) {
    EdStruct *xx = edc->xx;

    if (!xx || xx->editorState == StateDown)
	return;

    DBsetFlags(xx, edc->seq, DB_Flags(xx, edc->seq) | DB_FLAG_TRACE_SHOWN);

    RedisplayName(xx, edc->seq);
    redisplaySequences(xx, 1);
}

/*
 * gets the locking action
 */
int tman_get_lock(EdStruct *xx) {
    return xx->trace_lock;
}

/*
 * Sets/clears the locking action
 */
void tman_set_lock(EdStruct *xx, int val) {
    if (inJoinMode(xx) && xx->link) {
	xx->link->xx[0]->trace_lock = val;
	xx->link->xx[1]->trace_lock = val;
    } else {
	xx->trace_lock = val;
    }
}

/*
 * Acknowledges a shutdown of the contig editor window.
 * As this may still leave the traces window visible, we need to clear the
 * xx structures so that we do not try to access this free()d memory.
 */
void tman_notify_exit(void) {
    int i;

    for (i = 0; i < MAXCONTEXTS; i++) {
	edc[i].xx = (EdStruct *)0;
    }
}

/*
 * Shuts down all traces in the current trace display.
 *
 * limit_to is 0 for no limits (all traces), 1 for mini_only, 2 for
 * full only.
 */
void tman_shutdown_traces(EdStruct *xx, int limit_to) {
    int i;

    /* Shut down existing traces */
    for (i = 0; i < MAXCONTEXTS; i++) {
	if (edc[i].dc != NULL && edc[i].xx == xx) {
	    if (limit_to == 1 && edc[i].dc->mini_trace == 0)
		continue;
	    if (limit_to == 2 && edc[i].dc->mini_trace != 0)
		continue;
	    deleteTrace(xx, edc[i].dc->path);
	    edc[i].dc = NULL;
	}
    }
}

/*
 * Which trace types to bring up when auto-displaying traces.
 * See the code in tman_problem_traces() and tman_init_problem_traces() for
 * what these numbers mean.
 */
static int auto_traces[100] = {0, 3, 6, -1};

/*
 * Parses the specification of which traces to display at a problem, and
 * puts them in the auto_traces[] array. This must be kept in sync with
 * the tman_problem_traces() function.
 */
void tman_init_problem_traces(char *spec_orig) {
    int i = 0;
    char *token;
    char *spec;

    spec = strdup(spec_orig);

    token = strtok(spec, "\t ,/:");
    while (token) {
	int plus = 0, minus = 0;
	int ind = 0;

	if (*token == '+') {
	    plus = 1;
	    token++;
	} else if (*token == '-') {
	    minus = 1;
	    token++;
	}
	if (*token == '2') {
	    ind += 10;
	    token++;
	}

	switch (*token) {
	case 'd':
	case 'D':
	    auto_traces[i] = ind + (plus ? 1 : (minus ? 2 : 0));
	    break;

	case 'p':
	case 'P':
	    auto_traces[i] = ind + (plus ? 4 : 7);
	    break;

	case 't':
	case 'T':
	    auto_traces[i] = ind + (plus ? 5 : 8);
	    break;

	default:
	    auto_traces[i] = ind + (plus ? 3 : 6);
	}

	i++;
	token = strtok(NULL, "\t ,/:");
    }
    auto_traces[i] = -1;

    xfree(spec);
}

/*
 * Shuts down all currently displayed traces and displays those which are
 * suitable for solving a problem at the specified position.
 */
void tman_problem_traces(EdStruct *xx, int pos) {
    int *seqList, i, tmp;
    char c;
    struct {
	int seq;
	int score;
    } found[20] = {{0,-1}, {0,-1}, {0,-1}, {0,-1}, {0,-1},
		   {0,-1}, {0,-1}, {0,-1}, {0,-1}, {0,-1},
		   {0,-1}, {0,-1}, {0,-1}, {0,-1}, {0,-1},
		   {0,-1}, {0,-1}, {0,-1}, {0,-1}, {0,-1}};

    tman_shutdown_traces(xx, 2);

    /* Gather the list of traces covering this point */
    tmp = xx->reveal_cutoffs;
    xx->reveal_cutoffs = 0;
    seqList = sequencesInRegion(xx, pos, 1);
    xx->reveal_cutoffs = tmp;

    /*
     * What's the consensus character at this point?
     * Run this at 1% to ensure we always get a good estimate of the
     * correct base, rather than simply getting dashes all the while.
     */
    {
	float tmp = xx->con_cut;

	xx->con_cut = 0.01;
	DBcalcConsensus(xx, pos, 1, &c,  NULL, BOTH_STRANDS);
	xx->con_cut = tmp;
    }

    /*
     * Find suitable readings.
     * We're looking for any traces from this set:
     *     found[0]	1st reading conflicting with the consensus (any)
     *     found[1]	1st +ve reading conflicting with the consensus
     *     found[2]	1st -ve reading conflicting with the consensus
     *     found[3]	1st good +ve sense (any) (use quality & consensus)
     *     found[4]	1st good +ve sense primer (use quality & consensus)
     *     found[5]	1st good +ve sense terminator (use quality & consensus)
     *     found[6]	1st good -ve sense (any) (use quality & consensus)
     *     found[7]	1st good -ve sense primer (use quality & consensus)
     *     found[8]	1st good -ve sense terminator (use quality & consensus)
     *     ...
     *     found[10]	2nd reading conflicting with the consensus (any)
     *     found[11]	2nd +ve reading conflicting with the consensus
     *     found[12]	2nd -ve reading conflicting with the consensus
     *     found[13]	2nd good +ve sense (any) (use quality & consensus)
     *     found[14]	2nd good +ve sense primer (use quality & consensus)
     *     found[15]	2nd good +ve sense terminator (use quality & consensus)
     *     found[16]	2nd good -ve sense (any) (use quality & consensus)
     *     found[17]	2nd good -ve sense primer (use quality & consensus)
     *     found[18]	2nd good -ve sense terminator (use quality & consensus)
     */
    for (i=0; seqList[i]; i++) {
	int seq = seqList[i], p = pos - DB_RelPos(xx, seqList[i]) + 1;
	int score;
	int term = (DB_Flags(xx, seq) & DB_FLAG_TERMINATOR) > 0;

	score = getQual(xx, seq, p);

	if (DB_WSeq(xx, seq)[p-1] == c) {
	    /* An agreeing reading */
	    if (DB_Comp(xx, seq) == UNCOMPLEMENTED) {
		if (score > found[4+term].score) {
		    found[14+term] = found[4+term];
		    found[4+term].score = score;
		    found[4+term].seq = seq;
		} else if (score > found[14+term].score) {
		    found[14+term].score = score;
		    found[14+term].seq = seq;
		}
		if (score > found[3].score) {
		    found[13] = found[3];
		    found[3].score = score;
		    found[3].seq = seq;
		} else if (score > found[13].score) {
		    found[13].score = score;
		    found[13].seq = seq;
		}
	    } else {
		if (score > found[7+term].score) {
		    found[17+term] = found[7+term];
		    found[7+term].score = score;
		    found[7+term].seq = seq;
		} else if (score > found[17+term].score) {
		    found[17+term].score = score;
		    found[17+term].seq = seq;
		}
		if (score > found[6].score) {
		    found[16] = found[6];
		    found[6].score = score;
		    found[6].seq = seq;
		} else if (score > found[16].score) {
		    found[16].score = score;
		    found[16].seq = seq;
		}
	    }
	} else {
	    /* A disagreeing reading */
	    if (score > found[1+term].score) {
		found[11+term] = found[1+term];
		found[1+term].score = score;
		found[1+term].seq = seq;
	    } else if (score > found[11+term].score) {
		found[11+term].score = score;
		found[11+term].seq = seq;
	    }
	    if (score > found[0].score) {
		found[10] = found[0];
		found[0].score = score;
		found[0].seq = seq;
	    } else if (score > found[10].score) {
		found[10].score = score;
		found[10].seq = seq;
	    }
	}
    }

    /* Display the traces */
    i = 0;
    while (auto_traces[i] != -1) {
	if (found[auto_traces[i]].seq) {
	    int tmp = xx->compare_trace;
	    xx->compare_trace = -1;

	    showTrace(xx, found[auto_traces[i]].seq,
		      pos - DB_RelPos(xx, found[auto_traces[i]].seq) + 1,
		      xx->fontWidth * 2, 0, 0);

	    xx->compare_trace = tmp;
	}
	i++;
    }
}

/*
 * Remove all traces displayed for a particular DBinfo. The reason for this
 * is that the edc[i].seq numbers will change when one DBinfo is merged with
 * another due to joining contigs. It is possible to recalculate, but quite
 * fiddle to get the timing correct.
 */
void tman_handle_join(DBInfo *old, DBInfo *new) {
    int i;

    for (i = 0; i < MAXCONTEXTS; i++) {
	if (edc[i].dc == NULL || DBI(edc[i].xx) != old)
	    continue;

	deleteTrace(edc[i].xx, edc[i].dc->path);
    }
}

/*
 * Diff two traces. If they're not already displayed them make them so.
 * Returns 0 for success, -1 for failure.
 */
DisplayContext *diff_traces(EdStruct *xx, int seq1, int seq2, int pos) {
    int i;
    tman_dc *edc1 = NULL, *edc2 = NULL;

    tman_shutdown_traces(xx, 2);

    /* The first trace could be the consensus */
    if (seq1 != 0)
	showTrace(xx, seq1, pos - DB_RelPos(xx, seq1) + 1,
		  xx->fontWidth * 2, 1, 0);
    else
	cons_edc_trace(xx, DB_RelPos(xx, seq2),
		       DB_RelPos(xx, seq2) + DB_Length(xx, seq2) - 1,
		       DB_Comp(xx, seq2),
		       xx->compare_trace_match,
		       xx->compare_trace_select ? seq2 : 0);

    showTrace(xx, seq2, pos - DB_RelPos(xx, seq2) + 1,
	      xx->fontWidth * 2, 1, 0);

    for (i = 0; i < MAXCONTEXTS; i++) {
	if (edc[i].dc && edc[i].seq == seq1)
	    edc1 = &edc[i];
	if (edc[i].dc && edc[i].seq == seq2)
	    edc2 = &edc[i];
    }

    if (NULL == edc1 || NULL == edc2) {
	bell();
	return NULL;
    }

    return diff_edc_traces(xx, edc1, edc2);
}

/*
 * Given a sequence 'seq', find a sequence which is from the same template
 * but is on the opposite strand. If more than one exists then we want the one
 * with the largest overlap. seqlist is a list of sequence to scan through as
 * we've already narrowed down the field somewhat.
 *
 * Returns an editor sequence number or 0 for failure.
 */
static int auto_diff_sibling(EdStruct *xx, int seq, int *seqlist) {
    int i, sibling_seq = -1, sibling_overlap = 0;

    for (i = 0; seqlist[i]; i++) {
	int left, right;

	int j = seqlist[i];
	if (DB_Comp(xx, seq) == DB_Comp(xx, j))
	    continue;

	if (DBI_DB(xx)[seq].template != DBI_DB(xx)[j].template)
	    continue;

	left = MAX(DB_RelPos(xx, seq), DB_RelPos(xx, j));
	right = MIN(DB_RelPos(xx, seq) + DB_Length(xx, seq)-1,
		    DB_RelPos(xx, j) + DB_Length(xx, j)-1);

	if (sibling_overlap < right - left) {
	    sibling_overlap = right - left;
	    sibling_seq = j;
	}
    }

    return sibling_seq;
}

/*
 * Given sequence 'seq', identify the positive and negative reference
 * traces, chosing from 'seqlist' possibilities.
 *
 * Returns, plus and minus as editor sequence numbers. Zero indicates
 * none found.
 */
static void auto_diff_references(EdStruct *xx, int seq, int *seqlist,
				 int *plus, int *minus) {
    int seq_5; /* pos of 5' end */
    int closest_seq_p = 0;
    int closest_dist_p = INT_MAX;
    int closest_seq_m = 0;
    int closest_dist_m = INT_MAX;
    int i;

    /* Find 5' position of our comparison sequence */
    seq_5 = DB_Comp(xx, seq) == UNCOMPLEMENTED
	? DB_RelPos(xx, seq)
	: DB_RelPos(xx, seq) + DB_Length(xx, seq) - 1;

    /* Identify which reference traces are the closest */
    for (i = 0; seqlist[i]; i++) {
	int j = seqlist[i];
	int pos_5;

	/* Must be a reference trace */
	if (!(DB_Flags(xx, j) & DB_FLAG_REFTRACE))
	    continue;

	/* Must be on the same strand as 'seq' */
	if (DB_Comp(xx, j) != DB_Comp(xx, seq))
	    continue;

	pos_5 = DB_Comp(xx, j) == UNCOMPLEMENTED
	? DB_RelPos(xx, j)
	: DB_RelPos(xx, j) + DB_Length(xx, j) - 1;

	if (DB_Flags(xx, j) & DB_FLAG_REFTRACE_NEG) {
	    if (closest_dist_m > ABS(pos_5 - seq_5)) {
		closest_dist_m = ABS(pos_5 - seq_5);
		closest_seq_m = j;
	    }
	}

	if (DB_Flags(xx, j) & DB_FLAG_REFTRACE_POS) {
	    if (closest_dist_p > ABS(pos_5 - seq_5)) {
		closest_dist_p = ABS(pos_5 - seq_5);
		closest_seq_p = j;
	    }
	}
    }

    *plus = closest_seq_p;
    *minus = closest_seq_m;
}


/*
 * Attempts to guess positive and negative references by looking for good
 * quality traces that agree and disagree with the specified sequence 'seq'.
 * The sequences chosen will match the same strand as 'seq'.
 * We only pick new references if the passed in 'plus' and 'minus' sequence
 * numbers are zero.
 *
 * Returns: new plus/minus sequences, or zero if that type is not found.
 */
typedef struct {
    int read;
    int template;
    int conf;
    char base;
    int strand;
} read_template_pair;

static int read_template_sort(const void *p1, const void *p2) {
    const read_template_pair *rt1 = (read_template_pair *)p1;
    const read_template_pair *rt2 = (read_template_pair *)p2;
    return rt1->template - rt2->template;
}

static void guess_references(EdStruct *xx, int seqtop, int seqbot, int pos,
			     int *seqlist, int *rtop, int *rbot) {
    int i, j, count;
    int seq_pos;
    char seq_base1, seq_base2;
    int best_conf_match = 5, best_conf_mismatch = 5;
    int mismatch_seq_top = 0, mismatch_seq_bot = 0;
    read_template_pair *tlist; /* templates in seqlist[] */

    /* The base call for the query sequence */
    DBgetSeq(DBI(xx), seqtop);
    seq_pos = pos - DB_RelPos(xx, seqtop) + DB_Start(xx, seqtop);
    seq_base1 = DB_Seq(xx, seqtop)[seq_pos];
    DBgetSeq(DBI(xx), seqbot);
    seq_pos = pos - DB_RelPos(xx, seqbot) + DB_Start(xx, seqbot);
    seq_base2 = DB_Seq(xx, seqbot)[seq_pos];

    /* Identify which traces are close & good quality, store in tlist[] */
    for (i = 0; seqlist[i]; i++)
	;
    tlist = xcalloc(i, sizeof(*tlist));
    for (i = count = 0; seqlist[i]; i++) {
	int j = seqlist[i];
	char base;
	int conf;

	if (j == seqtop || j == seqbot)
	    continue;

	/* Determine if this is a better pos/neg base */
	DBgetSeq(DBI(xx), j);
	seq_pos = pos - DB_RelPos(xx, j) + DB_Start(xx, j);
	if (seq_pos < DB_Start(xx, j) || seq_pos+2 > DB_End(xx, j))
	    continue;
	base = DB_Seq(xx, j)[seq_pos];
	conf = DB_Conf(xx, j)[seq_pos];

	tlist[count].read = j;
	tlist[count].template = DBI_DB(xx)[j].template;
	tlist[count].base = base;
	tlist[count].conf = conf;
	tlist[count].strand = DB_Comp(xx, j);
	count++;
    }

    /* Resort list */
    qsort(tlist, count, sizeof(*tlist), read_template_sort);
    
    /* Remove duplicates where a template occurs more than one per strand */
    for (i = 0; i < count; i++) {
	for (j = i+1; j < count && tlist[j].template==tlist[i].template; j++) {
	    if (tlist[j].strand == tlist[i].strand) {
		/* duplicates on one strand */
		if (tlist[j].conf > tlist[i].conf) {
		    tlist[i].conf = 0;
		} else {
		    tlist[j].conf = 0;
		}
	    }
	}
    }

    /* Remove cases where the basecalls on top/bot strand disagree */
    for (i = 0; i < count; i++) {
	if (tlist[i].conf == 0)
	    continue;
	for (j = i+1; j < count && tlist[j].template==tlist[i].template; j++) {
	    if (tlist[j].conf == 0)
		continue;

	    if (tlist[i].base != tlist[j].base) {
		tlist[i].conf = 0;
		tlist[j].conf = 0;
	    }
	}
    }

    /* Now combine +/- strand scores together & pick best scores */
    for (i = 0; i < count; i++) {
	int score = tlist[i].conf;
	if (score == 0)
	    continue;
	for (j = i+1; j < count && tlist[j].template==tlist[i].template; j++) {
	    score += tlist[j].conf;
	}
	for (; i < j; i++) {
	    tlist[i].conf = score;

	    if (best_conf_match <= score &&
		(tlist[i].base == seq_base1 ||
		 tlist[i].base == seq_base2)) {
		best_conf_match = score;
		
	    }
	    if (best_conf_mismatch <= score && 
		(tlist[i].base != seq_base1 ||
		 tlist[i].base != seq_base2)) {
		best_conf_mismatch = score;
		if (tlist[i].strand == UNCOMPLEMENTED)
		    mismatch_seq_top = tlist[i].read;
		else
		    mismatch_seq_bot = tlist[i].read;
	    }
	}
	i--;
    }

    *rtop = mismatch_seq_top;
    *rbot = mismatch_seq_bot;
}

static void trace_columns(EdStruct *xx, int cols) {
    char *edpath;
    Tcl_Interp *interp = EDINTERP(xx->ed);
    char buf[10];

    /*
     * If we're the bottom half of a join editor, combine traces with the
     * top half.
     */
    if (inJoinMode(xx) && xx->link && xx == xx->link->xx[1]) {
	edpath = Tk_PathName(EDTKWIN(xx->link->xx[0]->ed));
    } else {
	edpath = Tk_PathName(EDTKWIN(xx->ed));
    }

    if (cols < 1) cols = 1;
    if (cols > 4) cols = 4;
    sprintf(buf, "%d", cols);
    /*
    Tcl_SetVar(interp, "update idletasks; trace_columns",
	       buf, TCL_GLOBAL_ONLY);
    */
    Tcl_SetVar(interp, "trace_columns",
	       buf, TCL_GLOBAL_ONLY);
}

static tman_dc *seq2edc(int seq) {
    int i;
    for (i = 0; i < MAXCONTEXTS; i++) {
	if (edc[i].dc && edc[i].seq == seq && edc[i].type != TRACE_TYPE_MINI)
	    return &edc[i];
    }

    return NULL;
}

/*
 * Called when the user has selected to automatically perform trace
 * differencing on their reference traces (both the wildtype/negative control
 * and optionally a positive control).
 * This may also be called in a limited mode where we wish to observe the
 * fwd/rev traces for a template, but not to automatically show
 * references and differences. The code looks at xx->diff_traces to determine
 * which mode to use.
 *
 * 'contig_pos' is a contig position rather than a sequence position.
 *
 * 'dcs' is a 2 dimensional array of Display Context pointers. If non-null it
 * will be filled out with the DisplayContexts for the displayed traces.
 *
 * Returns 0 for success
 *	  -1 for failure
 */
int auto_diff(EdStruct *xx, int seq, int contig_pos) {
    int seq_5;
    int *seqlist;
    int ref_pos_top, ref_pos_bot;
    int ref_neg_top, ref_neg_bot;
    int sibling;
    int top_seq, bot_seq;
    int columns;
    int pos;
    DisplayContext *ref_pos_top_dc = NULL, *ref_pos_bot_dc = NULL;
    DisplayContext *ref_neg_top_dc = NULL, *ref_neg_bot_dc = NULL;
    DisplayContext *top_seq1_dc    = NULL, *top_seq2_dc    = NULL;
    DisplayContext *bot_seq1_dc    = NULL, *bot_seq2_dc    = NULL;

    /* Find sequences 1000 base pairs either side of the 5' position */
    seq_5 = DB_Comp(xx, seq) == UNCOMPLEMENTED
	? DB_RelPos(xx, seq)
	: DB_RelPos(xx, seq) + DB_Length(xx, seq) - 1;
    seqlist = sequencesInRegion(xx, seq_5 - 1000, 2000);
    if (!seqlist) {
	bell();
	return -1;
    }

    /* Find top and bottom strand mutant sequences */
    top_seq = bot_seq = 0;
    sibling = auto_diff_sibling(xx, seq, seqlist);

    if (DB_Comp(xx, seq) == UNCOMPLEMENTED) {
	top_seq = seq;
	bot_seq = sibling != -1 ? sibling : 0;
    } else {
	top_seq = sibling != -1 ? sibling : 0;
	bot_seq = seq;
    }

    /* Find positive and negative controls */
    ref_pos_top = ref_pos_bot = 0;
    ref_neg_top = ref_neg_bot = 0;
    if (xx->diff_traces) {
	if (top_seq) {
	    auto_diff_references(xx, top_seq, seqlist, &ref_pos_top, &ref_neg_top);
	}
	if (bot_seq) {
	    auto_diff_references(xx, bot_seq, seqlist, &ref_pos_bot, &ref_neg_bot);
	}

	/* If negative references do not exist, pick something good look instead */
	if (!ref_neg_top && !ref_neg_bot)
	    guess_references(xx, top_seq, bot_seq, contig_pos, seqlist, 
			     &ref_neg_top, &ref_neg_bot);
    }


    /* No references found - revert to just top/bot sequence */
    if (!xx->diff_traces ||
	(!ref_pos_top && !ref_pos_bot && !ref_neg_top && !ref_neg_bot)) {
	if (top_seq) {
	    pos = contig_pos - DB_RelPos(xx, top_seq) + 1;
	    showTrace(xx, top_seq, pos, xx->fontWidth * 2, 1, 0);
	}
	if (bot_seq) {
	    pos = contig_pos - DB_RelPos(xx, bot_seq) + 1;
	    showTrace(xx, bot_seq, pos, xx->fontWidth * 2, 1, 0);
	}
	return 0;
    }

    columns = 0;
    if (top_seq && ref_neg_top) columns++;
    if (bot_seq && ref_neg_bot) columns++;
    if (top_seq && ref_pos_top) columns++;
    if (bot_seq && ref_pos_bot) columns++;

    /* Shutdown any existing traces and display N columns of new traces */
    tman_shutdown_traces(xx, 2);
    trace_columns(xx, columns);

    if (top_seq && ref_neg_top) {
	pos = contig_pos - DB_RelPos(xx, top_seq) + 1;
	top_seq1_dc = showTrace(xx, top_seq, pos, xx->fontWidth * 2, 1, 0);
	if (!top_seq1_dc) {
	    bell();
	    top_seq = 0;
	}
    }
    if (bot_seq && ref_neg_bot) {
	pos = contig_pos - DB_RelPos(xx, bot_seq) + 1;
	bot_seq1_dc = showTrace(xx, bot_seq, pos, xx->fontWidth * 2, 1, 0);
	if (!bot_seq1_dc) {
	    bell();
	    bot_seq = 0;
	}
    }
    if (top_seq && ref_pos_top) {
	pos = contig_pos - DB_RelPos(xx, top_seq) + 1;
	top_seq2_dc = showTrace(xx, top_seq, pos, xx->fontWidth * 2, 1, 0);
	if (!top_seq2_dc) {
	    bell();
	    top_seq = 0;
	}
    }
    if (bot_seq && ref_pos_bot) {
	pos = contig_pos - DB_RelPos(xx, bot_seq) + 1;
	bot_seq2_dc = showTrace(xx, bot_seq, pos, xx->fontWidth * 2, 1, 0);
	if (!bot_seq2_dc) {
	    bell();
	    bot_seq = 0;
	}
    }
    
    if (top_seq && ref_neg_top) {
	pos = contig_pos - DB_RelPos(xx, ref_neg_top) + 1;
	ref_neg_top_dc = showTrace(xx, ref_neg_top, pos, xx->fontWidth * 2, 1,
				   0);
	if (!ref_neg_top_dc) {
	    deleteTraceDisplay(xx, top_seq1_dc);
	    bell();
	    columns--;
	    ref_neg_top = 0;
	} else {
	    tman_dc *edc = find_edc(ref_neg_top_dc);
	    edc->type = TRACE_TYPE_NEG_CONTROL;
	}
    }
    if (bot_seq && ref_neg_bot) {
	pos = contig_pos - DB_RelPos(xx, ref_neg_bot) + 1;
	ref_neg_bot_dc = showTrace(xx, ref_neg_bot, pos, xx->fontWidth * 2, 1,
				   0);
	if (!ref_neg_bot_dc) {
	    deleteTraceDisplay(xx, bot_seq1_dc);
	    bell();
	    columns--;
	    ref_neg_bot = 0;
	} else {
	    tman_dc *edc = find_edc(ref_neg_bot_dc);
	    edc->type = TRACE_TYPE_NEG_CONTROL;
	}
    }
    if (top_seq && ref_pos_top) {
	pos = contig_pos - DB_RelPos(xx, ref_pos_top) + 1;
	ref_pos_top_dc = showTrace(xx, ref_pos_top, pos, xx->fontWidth * 2, 1,
				   0);
	if (!ref_pos_top_dc) {
	    deleteTraceDisplay(xx, top_seq2_dc);
	    bell();
	    columns--;
	    ref_pos_top = 0;
	} else {
	    tman_dc *edc = find_edc(ref_pos_top_dc);
	    edc->type = TRACE_TYPE_POS_CONTROL;
	}
    }
    if (bot_seq && ref_pos_bot) {
	pos = contig_pos - DB_RelPos(xx, ref_pos_bot) + 1;
	ref_pos_bot_dc = showTrace(xx, ref_pos_bot, pos, xx->fontWidth * 2, 1,
				   0);
	if (!ref_pos_bot_dc) {
	    deleteTraceDisplay(xx, bot_seq2_dc);
	    bell();
	    columns--;
	    ref_pos_top = 0;
	} else {
	    tman_dc *edc = find_edc(ref_pos_bot_dc);
	    edc->type = TRACE_TYPE_POS_CONTROL;
	}
    }

    if (top_seq && ref_neg_top) {
	diff_edc_traces(xx, seq2edc(top_seq), seq2edc(ref_neg_top));
    }
    if (bot_seq && ref_neg_bot) {
	diff_edc_traces(xx, seq2edc(bot_seq), seq2edc(ref_neg_bot));
    }
    if (top_seq && ref_pos_top) {
	diff_edc_traces(xx, seq2edc(top_seq), seq2edc(ref_pos_top));
    }
    if (bot_seq && ref_pos_bot) {
	diff_edc_traces(xx, seq2edc(bot_seq), seq2edc(ref_pos_bot));
    }

    trace_columns(xx, columns);

    return 0;
}

/*
 * From a single scroll of trace named 'path', this scrolls all other
 * traces in the same trace display. Ie a real trace->trace lock mode.
 */
void edScrollTraces(EdStruct *xx, char *path, char *command) {
    DisplayContext *dc;
    int i, orig;
    DNATrace *t;
    Read *r;
    tman_dc *ed;
    int2 *opos;
    int last = 0;
    int comp;
    int scroll_dir;
    double gpos;
    int pos;
    int com_argc;
    char **com_argv;
    char *com_argv_tmp[5];
    double f1;
    int count;
    int seq;
    int scroll_mode = TK_SCROLL_UNITS;

    /* Find trace number */
    dc = trace_path_to_dc(path);
    ed = find_edc(dc);
    t = dc->tracePtr;
    r = t->read;
    xx = ed->xx;


    /* Process scrollbar callback */
    if (!strchr(command, ' ')) {
	/* Old syntax */
	pos = atoi(command);
    } else {
	/* New syntax */
	if (Tcl_SplitList(EDINTERP(xx->ed), command, &com_argc, &com_argv)
	    != TCL_OK)
	    return;

	com_argv_tmp[0] = "a";
	com_argv_tmp[1] = "b";
	com_argv_tmp[2] = com_argv[0];
	com_argv_tmp[3] = com_argv[1];
	com_argv_tmp[4] = com_argv[2];
	com_argc+=2;

	scroll_mode = Tk_GetScrollInfo(EDINTERP(xx->ed), com_argc,
				       com_argv_tmp, &f1, &count);
	switch (scroll_mode) {
	default:
	case TK_SCROLL_ERROR:
	    pos = -1;
	    break;

	case TK_SCROLL_MOVETO:
	    pos = f1 * r->NPoints;
	    break;

	case TK_SCROLL_PAGES:
	    pos = t->disp_offset + count * t->disp_width * 0.9;
	    break;

	case TK_SCROLL_UNITS:
	    pos = t->disp_offset + count;
	    /* FIXME: This is just a hack to workaround the FIXME below */
	    if (ed->derivative_seq) {
		int x;
		for (x = 0; x < count; x++)
		    edCursorRight(ed->xx);
		for (x = 0; x < -count; x++)
		    edCursorLeft(ed->xx);
		return;
	    }
	}
	Tcl_Free((char *)com_argv);
	if (pos == -1)
	    return;
    }

    if (!xx->trace_lock) {
	/* Just scroll the one trace without moving the editor cursor */
	char buf[2000];

	gpos = pos / (double)r->NPoints;
	sprintf(buf, "%s xview moveto %g", path, gpos);
	Tcl_Eval(EDINTERP(xx->ed), buf);
	return;
    }

    if (scroll_mode == TK_SCROLL_UNITS) {
	if (pos >= t->disp_offset)
	    scroll_dir = 1;
	else
	    scroll_dir = -1;
    } else {
	scroll_dir = 0; /* avoid gcc whinge */
    }

    pos += t->disp_width/2;

    /* Convert the pos to a trace base coordinate */
    if (t->comp) {
	for (i = 0; i < t->Ned; i++) {
	    if (r->basePos[t->edPos[i]-1] <= pos)
		break;
	}
	if (scroll_mode == TK_SCROLL_UNITS) {
	    if (i > 0 && scroll_dir == -1 && r->basePos[t->edPos[i-1]-1] > pos)
		i++;
	}
    } else {
	for (i = 0; i < t->Ned; i++) {
	    if (r->basePos[t->edPos[i]-1] >= pos)
		break;
	}
	if (scroll_mode == TK_SCROLL_UNITS) {
	    if (i > 0 && scroll_dir == 1 && r->basePos[t->edPos[i-1]-1] < pos)
		i++;
	}
    }

    seq = ed->derivative_seq ? ed->derivative_seq : ed->seq;
    if (ed->derivative_seq) {
	/*
	 * If this is a trace difference, then we now do all maths on the
	 * sequence it was "derived" from, as it has real original position
	 * arrays (etc).
	 */
	seq = ed->derivative_seq;
	comp = DB_Comp(xx, seq) == COMPLEMENTED;

	if (comp) {
	    opos = DB_Opos(xx, seq);
	    
	    /*
	     * FIXME: not correct for scrolling by arrow heads, but close
	     * enough to be workable using the main scrollbar.
	     * This is hard to fix - it is best solved by starting again with
	     * the difference trace being the same orientation of the sequence
	     * it is derived from!.
	     */
	    i = origpos(xx, seq, i + ed->derivative_offset);
	} else {
	    i += ed->derivative_offset;
	}

	pos = 0;
    } else {
	seq = ed->seq;
	pos = ed->pos;
	comp = DB_Comp(xx, seq) == COMPLEMENTED;
    }

    orig = i;

    /* Now convert the original trace base into an edited sequence base no. */
    if (opos = DB_Opos(xx, seq)) {
	if (comp) {
	    for (i = DB_Length2(xx, seq)-1; i >= 0; i--) {
		if (opos[i] >= orig)
		    break;
		if (opos[i])
		    last = i;
	    }
	    if (i < 0)
		i = 0;
	} else {
	    for (i = 0; i < DB_Length2(xx, seq); i++) {
		if (opos[i] >= orig)
		    break;
		if (opos[i])
		    last = i;
	    }
	}
	if (i < DB_Length2(xx, seq) &&
	    opos[i] != orig &&
	    (last && i && opos[last] && opos[i])) {
	    /* Interpolate */
	    i = ((double)(orig - opos[last]) / (opos[i] - opos[last])) *
		(i - last) + last;
	}
    } else {
	i = orig;
    }

    i = i - DB_Start(xx, seq) + 1;

    setCursorPosSeq(xx, i + pos, seq);
    showCursor(xx,xx->cursorSeq, xx->cursorPos);
    repositionTraces(xx);
}


/*
 * Part of the save_trace_images() code below.
 */
static void trace_report_table(dstring_t *html, EdStruct *xx,
			       dstring_t *seqname,
			       int seq, int pos,
			       mutation_t *mut, int ref_pos,
			       tman_dc *traces[3][2], int ncolumns,
			       int page_break,
			       char *extra_heading,
			       int skip_heading) {
    int row, col;
    char *user_comment[4];

    if (!skip_heading) {
	if (page_break) {
	    dstring_append(html, "<h2 style=\"page-break-before: always\">");
	} else {
	    dstring_append(html, "<h2>");
	}
	dstring_appendf(html, "Sample %s, position %d</h2>\n",
			dstring_str(seqname), ref_pos);
    }
    if (extra_heading) {
	if (skip_heading && page_break)
	    dstring_appendf(html,
			    "<h3 style=\"page-break-before: always\">"
			    "%s</h3>\n",
			    extra_heading);
	else
	    dstring_appendf(html,
			    "<h3>%s</h3>\n",
			    extra_heading);
    }

    dstring_appendf(html, 
		"<table border=1 cellspacing=0 cellpadding=5 width=\"50%%\">\n"
		"  <tbody>\n");

    /* Now step through the 2D array of traces producing an HTML table */
    for (row = 0; row < 3; row++) {
	/* First line of tag comments; stats from mutscan */
	if (row == 0) {
	    dstring_append(html, "    <tr>\n");
	    for (col = 0; col < ncolumns; col++) {
		char *tag_text = NULL;
		char *tag_type = NULL;
		int strand = 0;
		int trace_pos = 0;
		
		if (mut->seq_top == traces[row][col]->seq) {
		    tag_type = mut->tag_type_top;
		    tag_text = mut->tag_text_top;
		    strand = 1;
		}
		if (mut->seq_bot == traces[row][col]->seq) {
		    tag_type = mut->tag_type_bot;
		    tag_text = mut->tag_text_bot;
		    strand = 2;
		}
		dstring_append(html, "      <td valign=top align=center>\n");
		Tcl_VarEval(EDINTERP(xx->ed), traces[row][col]->dc->path,
			    " icursor", NULL);
		trace_pos = atoi(Tcl_GetStringResult(EDINTERP(xx->ed)));
		if (traces[row][col]->dc->complemented) {
		    /* Negate position */
		    trace_pos = traces[row][col]->dc->tracePtr->Ned
			- trace_pos;
		} else {
		    trace_pos++; /* Count from 1 and not zero */
		}
		dstring_appendf(html, "      <small>%.4s %c %d(%d)</small>\n",
				tag_type ? tag_type : "????",
				"?+-"[strand],
				ref_pos,
				trace_pos);

		user_comment[col] = 0;
		if (tag_text) {
		    char *next_line = strchr(tag_text, '\n');
		    size_t text_len;
		    if (next_line) {
			text_len = next_line - tag_text;
			user_comment[col] = next_line+1;
		    } else {
			text_len = strlen(tag_text);
		    }
		    dstring_appendf(html, "        <br><small>%.*s</small>\n",
				    (int)text_len, tag_text);
		} else {
		    dstring_append(html, "        <br><small>Mutation not tagged</small>\n");
		}
		dstring_append(html, "      </td>\n");
	    }
	    dstring_append(html, "    </tr>\n");
	}

	/* Images */
	dstring_append(html, "    <tr>\n");
	for (col = 0; col < ncolumns; col++) {
	    dstring_t *ds = dstring_create("");

	    if (traces[row][col]->type == TRACE_TYPE_NEG_CONTROL) {
		dstring_appendf(ds, "%s (-ve control)",
				traces[row][col]->dc->file);
	    } else if (traces[row][col]->type == TRACE_TYPE_POS_CONTROL) {
		dstring_appendf(ds, "%s (+ve control)",
				traces[row][col]->dc->file);
	    } else if (traces[row][col]->type == TRACE_TYPE_DIFF) {
		/* Assumption that row 0 and 1 are seq and control */
		dstring_appendf(ds, "(%s - %s)",
				traces[0][col]->dc->file,
				traces[1][col]->dc->file);
	    } else {
		dstring_appendf(ds, "%s", traces[row][col]->dc->file);
	    }

	    dstring_appendf(html,
			    "      <td valign=top align=center width=\"%d%%\">\n"
			    "         <img src=\"seq%d_%d_%d.png\" "
			    "alt=\"seq%d_%d_%d trace\"><br>\n"
			    "         <small>%s</small><br>\n"
			    "      </td>\n",
			    (int)(100/ncolumns),
			    seq, pos, col*4+row,
			    seq, pos, col*4+row,
			    dstring_str(ds));
	    dstring_destroy(ds);
	}
	dstring_append(html, "    </tr>\n");

	/* Sequence */
	switch (traces[row][0]->type) {
	case TRACE_TYPE_SEQ:
	    dstring_append(html, "   <tr>\n");
	    for (col = 0; col < ncolumns; col++) {
		DNATrace *t = traces[row][col]->dc->tracePtr;

		/* Position in sequence */
		int spos;
		int start, end;

		/* Inclusive base sequence ranges, +- 19 from cursor pos */
		spos = t->cursor_pos;
		start = MAX(spos-19, 0);
		end = MIN(spos+19, t->Ned-1);
		dstring_appendf(html,
				"    <td valign=top align=center"
				" style=\"white-space: nowrap\">\n"
				"      <small><tt>%.*s"
				"<span style=\"background-color: #ffAAAA\">"
				"%c</span>%.*s</small></tt>\n"
				"    </td>\n",
				spos-start, &t->edBases[start],
				t->edBases[spos],
				end-spos, &t->edBases[spos+1]);

	    }
	    dstring_append(html, "   </tr>\n");
	    break;

	case TRACE_TYPE_POS_CONTROL:
	case TRACE_TYPE_NEG_CONTROL:
	    for (col = 0; col < ncolumns; col++) {
		char *bases;
		int spos, start, end;

		/* Identify column span */
		int type = traces[row][col]->type;
		int ncols = 1;
		if (type == TRACE_TYPE_NEG_CONTROL) {
		    if (col+1 < ncolumns && traces[row][col+1]->type == type) {
			col++;
			ncols++;
		    }
		}

		if (type == TRACE_TYPE_NEG_CONTROL) {
		    spos = pos - DB_RelPos(xx, DBI(xx)->reference_seq);
		    bases = DBgetSeq(DBI(xx), DBI(xx)->reference_seq);
		    start = MAX(spos-19, 0);
		    end = MIN(spos+19, DB_Length(xx, DBI(xx)->reference_seq));
		} else {
		    int xseq = traces[row][col]->seq;
		    DBgetSeq(DBI(xx), xseq);
		    bases = _DB_Seq(DBI(xx), xseq);
		    spos = pos - (DB_RelPos(xx, xseq)-1) + DB_Start(xx, xseq);
		    start = MAX(spos-19, 0);
		    end = MIN(spos+19, DB_Length(xx, xseq));
		}

		dstring_appendf(html,
				"    <td valign=top align=center colspan=%d"
				" style=\"white-space: nowrap\">\n"
				"      <small><tt>%.*s"
				"<span style=\"background-color: #ffAAAA\">"
				"%c</span>%.*s</small></tt>\n"
				"    </td>\n",
				ncols,
				spos-start, &bases[start],
				bases[spos],
				end-spos, &bases[spos+1]);
	    }
	}

	/* Any tag user-added comments */
	if (row == 2) {
	    /*
	     * Need to firstly check if any comments exist, as some browsers
	     * create visible blank rows when no table cells are filled.
	     */
	    int doit = 0;
	    for (col = 0; col < ncolumns; col++) {
		if (user_comment[col]) {
		    doit = 1;
		    break;
		}
	    }

	    if (doit) {
		dstring_append(html, "    <tr>\n");
		for (col = 0; col < ncolumns; col++) {
		    dstring_append(html, "      <td>\n");	
		    if (user_comment[col]) {
			dstring_t *com = dstring_create(user_comment[col]);
			dstring_to_html(com);
			dstring_appendf(html,
					"        <small>%s</small>\n",
					dstring_str(com));
			dstring_destroy(com);
		    } else {
			dstring_append(html, "        <br>\n");
		    }
		    dstring_append(html, "      </td>\n");
		}
		dstring_append(html, "    </tr>\n");
	    }
	}
    }
    dstring_append(html,
		   "  </tbody>\n"
		   "</table>\n");
}

/*
 * Calls auto_diff() to bring up the 2x3 or 4x3 difference plots for a 'seq' 
 * at a given consensus 'pos'. 
 * The plots are then saved as PNG files with wrapping HTML code.
 *
 * Side effects: This leaves no trace displays open upon exit, and will shut
 * down any existing ones.
 *
 * Returns 0 for success,
 *        -1 for failure
 */
int save_trace_images(dstring_t *html, EdStruct *xx, int seq, int pos,
		      mutation_t *mut, int ref_pos, int page_break,
		      char *dir) {
    int i;
    int tmp_pos;
    int tmp_dts;
    int tmp_seq;
    int tmp_lps;
    int tmp_diff;
    int ncolumns;
    int row, col;
    dstring_t *seqname;
    tman_dc *traces[3][4];
    tman_dc *traces_neg[3][2];
    tman_dc *traces_pos[3][2];
    int ncolumns_neg, ncolumns_pos;
    int skip_heading;

    seqname = dstring_create(io_rname(DBI_io(xx), DB_Number(xx, seq)));
    dstring_escape_html(seqname);

    /*
     * Call auto_diff to display the traces.
     * Note that the difference traces appear to assume that cursorPos
     * is the correct location to start up at, so we temporarily move to
     * there first.
     */
    tman_shutdown_traces(xx, 2);
    tmp_pos = xx->cursorPos;
    tmp_seq = xx->cursorSeq;
    xx->cursorPos = pos;
    xx->cursorSeq = 0;
    tmp_dts = xx->diff_trace_size;
    tmp_diff = xx->diff_traces;
    xx->diff_trace_size = 10; /* speeds up differencing */
    xx->diff_traces = 1;
    tmp_lps = xx->lines_per_seq;
    xx->lines_per_seq = 1;
    auto_diff(xx, seq, pos);
    tman_reposition_traces(xx, pos, 0);
    xx->cursorPos = tmp_pos;
    xx->cursorSeq = tmp_seq;
    xx->diff_trace_size = tmp_dts;
    xx->lines_per_seq = tmp_lps;
    xx->diff_traces = tmp_diff;

    /* Count how many traces there are, to work out the number of columns */
    ncolumns = 0;
    for (i = 0; i < MAXCONTEXTS; i++) {
	if (edc[i].dc && edc[i].type != TRACE_TYPE_MINI)
	    ncolumns++;
    }

    if (ncolumns%3 != 0) {
	/* missing trace - no way to tell which one, so abort. */
	return -1;
    }

    ncolumns /= 3;

    dstring_appendf(html, "<a name=\"Sample_%s_%d\"></a>\n",
		    dstring_str(seqname), ref_pos);

    /* Loop through all traces producing PNG images */
    row = col = 0;
    memset(&traces[0][0], 0, 3*4*sizeof(traces[0][0]));
    for (i = 0; i < MAXCONTEXTS && row < 3; i++) {
	char buf[1024];
	Tcl_DString ds;

	if (!edc[i].dc || edc[i].type == TRACE_TYPE_MINI)
	    continue;

	Tcl_DStringInit(&ds);
	Tcl_DStringAppendElement(&ds, edc[i].dc->path);
	Tcl_DStringAppendElement(&ds, "save_image");
	sprintf(buf, "%s/seq%d_%d_%d.png", dir, seq, pos, col*4+row);
	Tcl_DStringAppendElement(&ds, buf);
	Tcl_DStringAppendElement(&ds, "150");

	Tcl_Eval(EDINTERP(xx->ed), Tcl_DStringValue(&ds));
	Tcl_DStringFree(&ds);

	traces[row][col] = &edc[i];
	if (++col == ncolumns) {
	    col = 0;
	    row++;
	}
    }

    /*
     * Identify which trace columns compare against -ve controls and
     * which compare against +ve controls.
     *
     * If we have more than 2 columns, then we split the table into two
     * based on whether it's a diff against the -ve or +ve control.
     */
    ncolumns_pos = ncolumns_neg = 0;
    if (ncolumns <= 2) {
	for (row = 0; row < 3; row++) {
	    for (col = 0; col < ncolumns; col++) {
		traces_neg[row][col] = traces[row][col];
	    }
	}
	ncolumns_neg = ncolumns;
    } else {
	for (col = 0; col < ncolumns; col++) {
	    if (traces[1][col]->type == TRACE_TYPE_POS_CONTROL) {
		for (row = 0; row < 3; row++)
		    traces_pos[row][ncolumns_pos] = traces[row][col];
		ncolumns_pos++;
	    } else {
		for (row = 0; row < 3; row++)
		    traces_neg[row][ncolumns_neg] = traces[row][col];
		ncolumns_neg++;
	    }
	}
    }

    /* Now generate the html for the table(s) */
    skip_heading = 0;
    if (ncolumns_neg) {
	trace_report_table(html, xx, seqname, seq, pos, mut, ref_pos,
			   traces_neg, ncolumns_neg, page_break,
			   ncolumns > 2
			     ? "Difference vs. negative control"
			     : "",
			   skip_heading);
	page_break = 1;
	skip_heading = 1;
    }

    if (ncolumns_pos) {
	trace_report_table(html, xx, seqname, seq, pos, mut, ref_pos,
			   traces_pos, ncolumns_pos, page_break,
			   ncolumns > 2
			     ? "Difference vs. positive control"
			     : "",
			   skip_heading);
    }

    tman_shutdown_traces(xx, 2);

    dstring_destroy(seqname);

    return 0;
}
