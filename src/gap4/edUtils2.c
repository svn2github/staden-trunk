#include <math.h>
#include <stddef.h> /* 11/1/99 johnt - for size_t on WINNT */
#include <ctype.h>

#include "qual.h"
#include "gap-dbstruct.h"
#include "edUtils.h"
#include "IO.h"
#include "io-reg.h"
#include "undo.h"
#include "edStructs.h"
#include "contigEditor.h"
#include "tagUtils.h"
#include "select.h"
#include "extend.h"
#include "fortran.h"
#include "fort.h"
#include "misc.h"
#include "tman_interface.h"
#include "locks.h"
#include "io_utils.h"
#include "IO2.h"
#include "xalloc.h"
#include "io_handle.h"
#include "active_tags.h"
#include "gap_globals.h"
#include "notes.h"

/* ------ Internal function declarations ----- */
static void initEdStruct(EdStruct *xx, int flags, int displayWidth);
static void updateJoinCursor(EdStruct *xx, int seq, int pos);

#ifdef CACHE_CONSENSUS_2
/*
 * Simple consensus caching scheme.
 *
 * Basically, the consensus is either entirely known, or not at all.
 * We call valid_consensus() to determine in a region (which is ignored) of
 * the consensus is valid. If 'xx->valid_consensus' is true then it is.
 *
 * Whenever we do _any_ change to the consensus, we invalidate the entire
 * thing. This even includes changing the consensus or quality cutoff
 * parameters.
 *
 * To consensus (and quality) itself is stored in the EdStruct and is
 * updated only when we ask for the entire consensus - such as when computing
 * the padded to unpadded base positions.
 */

/*
 * Resizes the consensus buffer.
 */
void resize_consensus(EdStruct *xx, int len) {
    len++; /* Just to be careful! */
    if (len > xx->consensus_len) {
	xx->consensus = xrealloc(xx->consensus, (size_t)(len * 1.2));
	xx->quality = xrealloc(xx->quality, (size_t)(len * 1.2 *sizeof(xx->quality[0])));
	xx->consensus_len = len * 1.2;
	xx->valid_consensus = 0;
    }
}

void invalidate_consensus(EdStruct *xx) {
    xx->valid_consensus = 0;
}

/*ARGSUSED*/
void invalidate_consensus_region(EdStruct *xx, int from, int to) {
    xx->valid_consensus = 0;
}

int valid_consensus(EdStruct *xx, int from, int to) {
    if (from >= 1 && to <= DB_Length(xx, 0))
	return xx->valid_consensus;
    else
	return 0;
}

#endif

#if !defined(CACHE_CONSENSUS_2)
void resize_consensus(EdStruct *xx, int len) {}
void invalidate_consensus(EdStruct *xx) {}
void invalidate_consensus_region(EdStruct *xx, int from, int to) {}
int valid_consensus(EdStruct *xx, int from, int to) { return 0; }
#endif

/*
 *----------------------------------------------------------------------------
 * EdStruct functions
 *---------------------------------------------------------------------------
 */

/*
 * Static variables
 * defining the state of the contig editor
 */
int edused[MAXEDSTATES] = {0};
EdStruct edstate[MAXEDSTATES] = {{0}};

EdStruct *editor_id_to_edstruct(int id) {
    return &edstate[id];
}

int editor_available(int contig, int nojoin) {
    int j;

    for (j=0; j<MAXEDSTATES; j++) {
	if (edused[j] && edstate[j].DBi &&
	    DBI_contigNum(&edstate[j]) == contig) {
	    if (nojoin == 0 || edstate[j].link == NULL)
		return j;
	}
    }

    return -1;
}

int move_editor(int id, int seq, int pos) {
    EdStruct *xx = &edstate[id];
    int i;

    for (i = 1; i <= DBI_gelCount(xx); i++) {
	if (DB_Number(xx, i) == seq) {
	    seq = i;
	    break;
	}
    }
    setCursorPosSeq(xx, pos, seq);
    redisplayWithCursor(xx);
    front_editor(xx);

    return 0;
}


int editor_select_region(int id, int seq, int pos, int len) {
    EdStruct *xx = &edstate[id];
    int i;

    for (i = 1; i <= DBI_gelCount(xx); i++) {
	if (DB_Number(xx, i) == seq) {
	    seq = i;
	    break;
	}
    }

    /* Display the selection */
    _select_region(xx, seq, pos, len);

    return 0;
}


/*
 * Get the next free EdStruct
 */
EdStruct *getFreeEdStruct(GapIO *io, int contig,
			  void (*dispFunc)(void *, int, int, int, void *))
{
    EdStruct *xx;
    int i, j, inuse = 0, flags;

    /* Find an unused one */
    for (i=0; i<MAXEDSTATES; i++)
	if (!edused[i])
	    break;

    if (i == MAXEDSTATES)
	return NULL;

    /* Check for whether this contig is in use */
    for (j=0; j<MAXEDSTATES; j++) {
	if (edused[j] && edstate[j].DBi &&
	    DBI_contigNum(&edstate[j]) == contig) {
	    inuse = 1;
	    break;
	}
    }

    xx = &edstate[i];
    edused[i] = 1;

    if (inuse)
	xx->DBi = edstate[j].DBi;
    else
	xx->DBi = NULL;

    flags = DB_DELAYED_READ | DB_DATA_TYPE_DNA;
    if (!io_rdonly(io))
	flags |= DB_ACCESS_UPDATE;
    initEdStruct(xx, flags, DEFAULT_DISPLAY_WIDTH);

    DBI_dispFunc(xx)[DBI_nextDisp(xx)] = dispFunc;
    DBI_dispData(xx)[DBI_nextDisp(xx)] = xx;
    DBI_nextDisp(xx)++;

    semaphoreGrab(activeLock);

    return xx;
}

/*
 * Removes an EdStruct. We only free memory here local to the edStruct. Any
 * DBInfo data maybe required by other edStructs. If not, then this will
 * already have been freed for us.
 */
static void destroyEdStruct(EdStruct *xx) {
    int i;

    /* find edstate number */
    for (i=0; i<MAXEDSTATES; i++) {
	if (xx == &edstate[i])
	    break;
    }

    /* Mark edstate as free */
    edused[i] = 0;
    edstate[i].DBi = NULL;

    /* Free memory */
    if (xx->tag_list)
	xfree(xx->tag_list);
#ifdef CACHE_CONSENSUS_2
    if (xx->consensus)
	xfree(xx->consensus);
    if (xx->quality)
	xfree(xx->quality);
#endif
    if (xx->set)
	xfree(xx->set);
    if (xx->set_collapsed)
	xfree(xx->set_collapsed);

    semaphoreRelease(activeLock);
}

void initEdStruct(EdStruct *xx, int flags, int displayWidth)
{
    static int editor_id = 0;
    int i;

    xx->editor_id = editor_id++;

    if (!xx->DBi) {
	if ((xx->DBi = (DBInfo *)xcalloc(sizeof(DBInfo),1)) == NULL)
	    return;

	xx->DBi->DB_flags = flags;
	xx->DBi->DB_gelCount = 0;
	xx->DBi->DB_contigNum = 0;
	xx->DBi->DBlist = NULL;
	xx->DBi->DBorder = NULL;
	xx->DBi->DB = NULL;
	xx->DBi->next_display = 0;
	xx->DBi->last_undo = 0;
	xx->DBi->num_undo = 0;
	xx->DBi->discarded_undos = 0;
	xx->DBi->edits_made = 0;
	xx->DBi->since_auto_save = 0;
	xx->DBi->store_undo = 1;
	xx->DBi->open_undo_count = 0;
	xx->DBi->registration_id = 0;
	xx->DBi->reference_seq = 0;
	xx->DBi->templates = NULL;
    }

    xx->displayPos = 1;
    xx->displayYPos = 0;
    xx->displayWidth = displayWidth;
    xx->displayHeight = 0;
    xx->totalHeight = 0;
    xx->cursorPos = 1; /* Don't use setCursorPos or setCursorSeq here as */
    xx->cursorSeq = 0; /* xx->DBi->DB hasn't been initialised yet */
    xx->rulerDisplayed = 1;
    xx->consensusDisplayed = 1;
    xx->fontWidth = 0;
    xx->fontHeight = 0;
    xx->displayedConsensus[0] = '\0';
    xx->reveal_cutoffs = 0;
    xx->showDifferences = 0;
    xx->compare_strands = 0;
    xx->display_traces = 0;
    xx->diff_traces = 0;
    xx->read_pair_traces = 0;
    xx->auto_save = 0;
    xx->link = NULL;
    xx->editorState = StateDown;
    xx->con_cut = 0.0;
    xx->qual_cut = 0;
    xx->consensus_mode = consensus_mode;
    xx->select_made = 0;
    xx->select_seq = 0;
    xx->select_start_pos = 0;
    xx->select_end_pos = 0;
    xx->select_tag = NULL;
    xx->insert_mode = 0;
    xx->super_edit = 0;
    xx->ed = NULL;
    for (i=0; i<10; i++)
	xx->qual_bg[i] = 0;
    xx->qual_below = 0;
    xx->extent_left = 0;
    xx->extent_right = 0;
    xx->sel_oli = NULL;
    xx->trans_mode = 1;
    xx->trace_lock = 0;
    xx->show_qual = 0;
    xx->show_cons_qual = 0;
    xx->tmp_tag = NULL;
    xx->refresh_flags = 0;
    xx->refresh_seq = 0;
    xx->refresh_pos = 0;
    xx->compare_trace = -1;
    xx->compare_trace_match = 1;
    xx->compare_trace_select = 1;
    xx->compare_trace_algorithm = 1;
    xx->compare_trace_yscale = 1;
    xx->group_mode = POSITION;
    xx->show_edits = 0;
    for (i=0; i<4; i++)
	xx->edit_bg[i] = 0;
    for (i=0; i<6; i++)
	xx->tmpl_bg[i] = 0;
    for (i=0; i<10; i++)
	xx->set_bg[i] = 0;
    xx->names_xpos = 0;
    xx->default_conf_r = 100;
    xx->default_conf_n = 100;
#ifdef CACHE_CONSENSUS_2
    xx->consensus = NULL;
    xx->quality = NULL;
    xx->valid_consensus = 0;
    xx->consensus_len = 0;
#endif
    xx->unpadded_ruler = 0;
    xx->status_lines = NULL;
    xx->status_depth = 0;
    xx->lines_per_seq = 1;
    xx->diff_trace_size = 0;
    xx->diff_qual = 0;
    xx->set = NULL;
    xx->set_collapsed = NULL;
    xx->nsets = 0;

    /* Set all tags to be displayed by default */
    if (NULL == (xx->tag_list = (int *)xmalloc(tag_db_count * sizeof(int))))
	return;

    for (i = 0; i < tag_db_count; i++) {
	xx->tag_list[i] = 1;
    }

    for (i = 0; i < MAX_STATUS_LINES; i++) {
	xx->status[i] = 0;
    }
}


/*
 * Converts a reading number to an editor sequence number.
 * Returns the sequence number or -1 if not found.
 */
int rnum_to_edseq(EdStruct *xx, int rnum) {
    int i;

    for (i=1; i<=DBI_gelCount(xx); i++) {
	if (DB_Number(xx, i) == rnum)
	    return i;
    }
    return -1;
}

/*
 * Create the initial sequence display
 */
int createEdDisplay(EdStruct *xx, int seq, int pos)
{
    /*
     * Initial position on screen
     */
    xx->displayPos = 1;
    setCursorPosSeq(xx, pos, 0);
    
    if (-1 != (seq = rnum_to_edseq(xx, seq))) {
	if (xx->cursorPos >= 1 &&
	    xx->cursorPos <= DB_Length(xx, seq)) {
	    setCursorPosSeq(xx, pos, seq);
	}
    }
    
    /*
     * Display
     * set xx->displayPos to force repositioning of cursor 
     */
    xx->displayPos = positionInContig(xx,xx->cursorSeq,xx->cursorPos) +
	2*xx->displayWidth;
    redisplayWithCursor(xx);
    return 0;
}


EdLink *CreateEdLink(EdStruct *xx0, EdStruct *xx1) {
    EdLink *el;

    if (NULL == (el = (EdLink *)xmalloc(sizeof(EdLink))))
	return NULL;

    el->xx[0] = xx0;
    el->xx[1] = xx1;
    el->lockOffset = 0;
    el->locked = 0;
    el->diffs = NULL;

    xx0->link = el;
    xx1->link = el;

    return el;
}

void DestroyEdLink(EdLink *el) {

    el->xx[0]->link = el->xx[1]->link = NULL;
    xfree(el);
}


/*
 * Call a function several times (once for each registered EdStruct) of
 * this DBInfo.
 */
void DBI_callback(DBInfo *dbi, int type, int seq, int pos, void *pointer) {
    static void (*display_func[MAX_DISP_PROCS])
	(void *xx, int type, int seq, int pos, void *);
    static void *display_data[MAX_DISP_PROCS];
    int x, y;
    
    /*
     * Take a temporary copy of the _DBI_dispFunc prior to calling them
     * incase the act of calling (as with DBCALL_QUIT) changes the
     * arrays. Ok, so it's not a truely reentrant routine, but it only matters
     * for these nasty cases (DBCALL_QUIT) which don't call DBI_callback
     * themselves.
     */
    for (x = y = 0; x < MAX_DISP_PROCS; x++) {
	if (_DBI_dispFunc(dbi)[x]) {
	    display_func[y] = _DBI_dispFunc(dbi)[x];
	    display_data[y] = _DBI_dispData(dbi)[x];
	    y++;
	}
    }

    for (x = 0; x < y; x++)
	display_func[x](display_data[x], type, seq, pos, pointer);
}


/*
 * Get an EdStruct from a DBInfo and editor_id.
 *
 * The abstraction breaks down a little bit here as we're assuming this
 * is an edstruct - our original scheme allows for any 'view' to be attached
 * to the DBInfo structure.
 */
EdStruct *DBI_to_EdStruct(DBInfo *dbi, int ed_id) {
    int x;

    for (x = 0; x < MAX_DISP_PROCS; x++) {
	if (_DBI_dispFunc(dbi)[x] &&
	    ((EdStruct *)_DBI_dispData(dbi)[x])->editor_id == ed_id)
	    return (EdStruct *)_DBI_dispData(dbi)[x];
    }

    return NULL;
}

/*
 *----------------------------------------------------------------------------
 * Consenus bits and pieces
 *---------------------------------------------------------------------------
 */

/*
 * Calculate dynamic consensus length
 */
int calculate_consensus_length(EdStruct *xx)
{
    int clen;
    int sequenceEnd,i;
    
    /* set initial consensus length */
    i = DBI_gelCount(xx);
    sequenceEnd = DB_RelPos(xx,DBI_order(xx)[i]) +
	DB_Length(xx,DBI_order(xx)[i]) - 1;
    clen = sequenceEnd;

    for(i--; i>=1; i--) {
	sequenceEnd = DB_RelPos(xx,DBI_order(xx)[i]) +
	    DB_Length(xx,DBI_order(xx)[i]) - 1;

	if (clen < sequenceEnd)
	    clen = sequenceEnd;
    }

    return clen;
}

/*
 * Set dynamic consensus length
 */
void calculateConsensusLength(EdStruct *xx)
{
    DBsetLength(xx,0,calculate_consensus_length(xx));
    DBsetLength2(xx, 0, DB_Length(xx, 0));
}


/*
 * `Info' function. Obtains various information from either the database
 * or from the contig editor structures. For this reason a pointer to the
 * relevant routine is supplied by the calling code.
 * This is the contig editor version.
 */
int contEd_info(int job, void *mydata, info_arg_t *theirdata) {
    EdStruct *xx = (EdStruct *)mydata;
    DBInfo *db = DBI(xx);

    switch (job) {
    case GET_SEQ:
	{
	    gel_seq_t *gel_seq = &theirdata->gel_seq;
	    register int gel = _DBI_order(db)[gel_seq->gel];

	    (void)DBgetSeq(db, gel);
	    gel_seq->gel_seq     = _DB_Seq    (db, gel);
	    /*
	     * Reference sequences have confidence 100 so that the consensus
	     * becomes the reference. We fabricate this on the fly (rather
	     * than editing the real data, which makes 'undo' trickier).
	     */
	    if (gel == db->reference_seq) {
		gel_seq->gel_conf = (int1 *)xmalloc(_DB_Length2(db, gel));
		memset(gel_seq->gel_conf, 100, _DB_Length2(db, gel));
	    } else {
		gel_seq->gel_conf    = _DB_Conf   (db, gel);
	    }
	    gel_seq->gel_opos    = _DB_Opos   (db, gel);
	    gel_seq->gel_length  = _DB_Length2(db, gel);
	    gel_seq->gel_end	 = _DB_End    (db, gel);
	    gel_seq->gel_start   = _DB_Start  (db, gel);

	    return 0;
	}
    case DEL_SEQ:
	{
	    gel_seq_t *gel_seq = &theirdata->gel_seq;
	    register int gel = _DBI_order(db)[gel_seq->gel];

	    /* Free reference seq confidence array */
	    if (gel == db->reference_seq) {
		xfree(gel_seq->gel_conf);
		gel_seq->gel_conf = NULL;
	    }
	    /* No need to free anything else */
	    return 0;
	}
    case GET_CONTIG_INFO:
	{
	    contig_info_t *contig_info = &theirdata->contig_info;
	    int gel = 1;
	    
	    contig_info->length = _DB_Length(db, 0);

	    contig_info->leftgel = 0;
	    for (gel = 1; gel <= _DBI_gelCount(db); gel = gel+1) {
		/* Skip if it's "hidden" for disassembly */
		if (_DB_Flags(db, _DBI_order(db)[gel]) & DB_FLAG_INVIS) {
		    continue;
		}

		/* Skip if we're showing a set and it's not in this one */
		if (xx->set && xx->curr_set &&
		    xx->set[_DBI_order(db)[gel]] != xx->curr_set) {
		    continue;
		}

		contig_info->leftgel = gel;
		break;
	    }

	    return 0;
	}
    case DEL_CONTIG_INFO:
	{
	    return 0;
	}
    case GET_GEL_INFO:
	{
	    gel_info_t *gel_info = &theirdata->gel_info;
	    register int gel = gel_info->gel;
	    register int gelo = _DBI_order(db)[gel];

	    gel_info->length       = _DB_Length(db, gelo);
	    gel_info->unclipped_len= _DB_Length2(db, gelo);
	    gel_info->complemented = _DB_Comp  (db, gelo) == COMPLEMENTED;
	    gel_info->position     = _DB_RelPos(db, gelo);
	    gel_info->as_double    = _DB_Flags (db, gelo) & DB_FLAG_TERMINATOR;
	    gel_info->start	   = _DB_Start (db, gelo);
	    gel_info->template     = _DB_Template(db, gelo);

	    gel_info->next_right = 0;
	    for (gel++; gel <= _DBI_gelCount(db); gel = gel+1) {
		/* Skip if it's "hidden" for disassembly */
		if (_DB_Flags(db, _DBI_order(db)[gel]) & DB_FLAG_INVIS) {
		    continue;
		}

		/* Skip if we're showing a set and it's not in this one */
		if (xx->set && xx->curr_set &&
		    xx->set[_DBI_order(db)[gel]] != xx->curr_set) {
		    continue;
		}

		gel_info->next_right = gel;
		break;
	    }

	    return 0;
	}
    case DEL_GEL_INFO:
	{
	    return 0;
	}
    case GET_GEL_LEN:
	{
	    return dbi_max_gel_len(db, 1);
	}
    default:
	verror(ERR_FATAL, "contEd_info", "Unknown job number (%d)", job);
	return -1;
    }
}


/*
 * calculate the consensus for position `pos' in contig,
 * for `width' characters
 */
void DBcalcConsensus (EdStruct *xx,int pos, int width, char *str,
		      float *qual, int mode)
{
    int tmp_mode;

    if (xx->compare_strands) {
	char *str2 = (char *)xmalloc(width+1);
	float *qual2;
	int i;
	
	if (qual) {
	    qual2 = (float *)xmalloc((width+1) * sizeof(float));
	    if (qual2 == NULL)
		return;
	} else
	    qual2 = NULL;

	if (str2 == NULL)
	    return;
	
	tmp_mode = consensus_mode;
	consensus_mode = xx->consensus_mode;
	calc_consensus(0, pos, pos + width - 1, CON_SUM,
		       str, str2, qual, qual2,
		       xx->con_cut,
		       xx->consensus_mode ? xx->qual_cut : -1,
		       contEd_info,
		       (void *)xx);
	consensus_mode = tmp_mode;

	for (i=0; i<width; i++) {
	    if (str[i] != str2[i]) {
		str[i] = '-';
		if (qual)
		    qual[i] = 0;
	    } else if (qual) {
		if (consensus_mode == CONSENSUS_MODE_CONFIDENCE) {
		    double e1, e2, e;

		    if (qual[i] != 100 || qual2[i] != 100) {
			/* Bayesian statistics */
			e1 = pow(10.0, -qual[i]/10.0);
			e2 = pow(10.0, -qual2[i]/10.0);
			e = 1-(1-e1)*(1-e2)/((1-e1)*(1-e2)+e1*e2/3);
			e = e ? -10*log10(e) : 99;
			qual[i] = MIN(e, 99);

			/*
			 * A quicker approximation would be
			 * qual[i] = MIN(qual[i] + qual2[i] + 4, 99);
			 */
		    } else {
			qual[i] = 100;
		    }
		} else {
		    qual[i] = (qual[i] + qual2[i]) / 2;
		}
	    }
	}

	if (qual2)
	    xfree(qual2);
	xfree(str2);
    } else {

	if (mode == BOTH_STRANDS) {
#ifdef CACHE_CONSENSUS_2
	    if (valid_consensus(xx, pos, pos + width - 1)) {
		memcpy(str, &xx->consensus[pos-1], width);
		if (qual)
		    memcpy(qual, &xx->quality[pos-1],
			   width * sizeof(*qual));
	    } else {
		if (pos == 1 && width == DB_Length(xx, 0)) {
		    resize_consensus(xx, width);
		    tmp_mode = consensus_mode;
		    consensus_mode = xx->consensus_mode;
		    calc_consensus(0, 1, width, CON_SUM,
				   xx->consensus, NULL, xx->quality, NULL,
				   xx->con_cut,
				   xx->consensus_mode ? xx->qual_cut : -1,
				   contEd_info,
				   (void *)xx);
		    consensus_mode = tmp_mode;
		    memcpy(str, xx->consensus, width);
		    if (qual)
			memcpy(qual, xx->quality, width * sizeof(*qual));
		    xx->valid_consensus = 1;
		} else {
#endif
		    tmp_mode = consensus_mode;
		    consensus_mode = xx->consensus_mode;
		    calc_consensus(0, pos, pos + width - 1, CON_SUM,
				   str, NULL, qual, NULL,
				   xx->con_cut,
				   xx->consensus_mode ? xx->qual_cut : -1,
				   contEd_info,
				   (void *)xx);
		    consensus_mode = tmp_mode;
#ifdef CACHE_CONSENSUS_2
		}
	    }
#endif
	} else {
	    /*
	     * Calculate only on a single strand. So pass over both arrays.
	     * Temp allocate second array and throw away afterwards.
	     */
	    
	    char *dummy1 = (char *)xmalloc(width+1);
	    float *dummy2 = NULL;

	    if (dummy1 == NULL)
		return;

	    if (qual) {
		dummy2 = (float *)xmalloc((width+1) * sizeof(float));
		if (dummy2 == NULL)
		    return;
	    }

	    if (mode == COMPLEMENTED) {
		tmp_mode = consensus_mode;
		consensus_mode = xx->consensus_mode;
		calc_consensus(0, pos, pos + width - 1, CON_SUM,
			       dummy1, str, dummy2, qual,
			       xx->con_cut,
			       xx->consensus_mode ? xx->qual_cut : -1,
			       contEd_info, (void *)xx);
		consensus_mode = tmp_mode;
	    } else {
		tmp_mode = consensus_mode;
		consensus_mode = xx->consensus_mode;
		calc_consensus(0, pos, pos + width - 1, CON_SUM,
			       str, dummy1, qual, dummy2,
			       xx->con_cut,
			       xx->consensus_mode ? xx->qual_cut : -1,
			       contEd_info, (void *)xx);
		consensus_mode = tmp_mode;
	    }
	    
	    xfree(dummy1);
	    if (dummy2)
		xfree(dummy2);
	}
    }
}

/*
 * As for DBcalcConsensus, this computes the consensus but for a specific
 * set only.
 */
void DBcalcSetConsensus(EdStruct *xx,int pos, int width, int set,
			char *str, float *qual) {
    int tmp_curr = xx->curr_set;
    xx->curr_set = set;
    DBcalcConsensus(xx, pos, width, str, qual, BOTH_STRANDS);
    xx->curr_set = tmp_curr;
}

/*
 * calculate the consensus discrepancy figures for position `pos' in contig,
 * for `width' characters
 */
void DBcalcDiscrepancies(EdStruct *xx,int pos, int width, float *qual)
{
    int tmp_mode;
    char *str = xmalloc(width);

    tmp_mode = xx->consensus_mode;
    consensus_mode = CONSENSUS_MODE_CONFIDENCE;
    calc_consensus(0, pos, pos + width - 1, CON_SUM,
		   str, NULL, NULL, qual,
		   xx->con_cut, 0,
		   contEd_info,
		   (void *)xx);
    consensus_mode = tmp_mode;
    xfree(str);
}


/*
 *----------------------------------------------------------------------------
 * Database IO
 *---------------------------------------------------------------------------
 */

/*
 * Read in the sequence for 'seq' and return it.
 */
char *DBgetSeq(DBInfo *db, int seq)
{
    int_f i = _DB_Number(db,seq);
    int length;
    
    /* already in memory? */
    if (!seq || _DB_Flags(db,seq) & DB_FLAG_SEQ_IN_MEMORY)
	return _DB_Seq(db,seq)+_DB_Start(db,seq);
    
    /*
     * Find total length of sequence.
     */
    (void)get_read_info(_DBI_io(db),
			i,
			NULL, 0, /* clone */
			NULL, 0, /* vec */
			NULL, 0, /* template */
			NULL, 0, /* template vec */
			&length,
			NULL, NULL, /* insert min/max */
			NULL,    /* dir */
			NULL,    /* strands */
			NULL,    /* primer */
			NULL, NULL, NULL, NULL /* ids for top four */
			);

    /*
     * allocate memory. We only allocate as much as we need plus a little
     * more to allow for edits.
     */
    
    length += SEQ_LENGTH_INC + 0.1 * length;

    /* force reading */
    {
	int tmp_len;
	GReadings r;

	(void) io_aread_seq(_DBI_io(db),
			    (int) i,
			    &tmp_len,
			    &_DBI_DB(db)[seq].gap_start,
			    &_DBI_DB(db)[seq].gap_end,
			    &_DBI_DB(db)[seq].sequence,
			    &_DBI_DB(db)[seq].gap_conf,
			    &_DBI_DB(db)[seq].gap_opos,
			    0);
	_DBI_DB(db)[seq].gap_length = tmp_len;
	_DBsetAlloced(db, seq, tmp_len);
	gel_read(_DBI_io(db), i, r);
	_DBI_DB(db)[seq].template = r.template;
    }
    
    /* mark as read */
    _DBsetFlags(db,seq,_DB_Flags(db,seq)|DB_FLAG_SEQ_IN_MEMORY);
    
    return _DB_Seq(db,seq)+_DB_Start(db,seq);
}


/*
 * Force tags into memory
 */
tagStruct *DBgetTags (DBInfo *db, int seq)
{
    int i;
    
    /* already in memory? */
    if (_DB_Flags(db,seq) & DB_FLAG_TAG_IN_MEMORY)
	return (tagStruct *)_DB_Tags(db,seq);

    /* ensure seq is loaded */
    DBgetSeq(db, seq);
    
    /* read in tag list */
    i = _DB_Number(db,seq);
    _DBsetTags(db,seq,readTagList(db,i,seq));
    
    /* mark as read */
    _DBsetFlags(db,seq,_DB_Flags(db,seq)|DB_FLAG_TAG_IN_MEMORY);
    
    return (tagStruct *)_DB_Tags(db,seq);
}


/*
 * Force reading in the sequence for seq
 */
char *DBgetName(DBInfo *db, int seq)
{
    int_f i;
    char buf[NAMELEN+1];
    
    /* already in memory? */
    if (!seq || _DB_Flags(db, seq) & DB_FLAG_NAME_IN_MEMORY)
	return _DB_Name(db, seq);
    
    /* allocate memory */
    if ((_DBsetName(db, seq,(char *)xmalloc(sizeof(char)*(NAMELEN+1))))==NULL)
	return NULL;
    
    /* force reading */
    i = _DB_Number(db, seq);
    readn_(handle_io(_DBI_io(db)), &i, buf, DB_NAMELEN); /* FORIO */
    buf[DB_NAMELEN]='\0';
    sprintf(_DB_Name(db, seq),"%+*d %-*s",
	    DB_GELNOLEN,
	    (_DB_Comp(db, seq) == COMPLEMENTED)
	      ? -_DB_Number(db, seq)
	      : _DB_Number(db, seq),
	    DB_NAMELEN,
	    buf);
    
    /* mark as read */
    _DBsetFlags(db, seq, _DB_Flags(db, seq) | DB_FLAG_NAME_IN_MEMORY);
    
    return _DB_Name(db, seq);
}


/*
 * Get the template name for sequence 'seq'.
 * Static storage, so non-reentrant.
 */
char *DBgetTemplateName(DBInfo *db, int seq)
{
    char *rname = DBgetName(db, seq);
    char tname[DB_NAMELEN+1];
    static char name[NAMELEN+1];
    int i;
    GReadings r;
    GTemplates t;
    
    /* Read template name */
    i = _DB_Number(db, seq);
    if (i <= 0) {
	return rname;
    }
    gel_read(_DBI_io(db), i, r);
    if (!r.template) {
	strcpy(tname, "(unknown)");
    } else {
	template_read(_DBI_io(db), r.template, t);
	if (!t.name) {
	    strcpy(tname, "(unknown)");
	} else {
	    TextRead(_DBI_io(db), t.name, tname, DB_NAMELEN);
	    tname[DB_NAMELEN] = 0;
	}
    }

    /* Combine template name with gel number */
    sprintf(name, "%.*s %-*s", DB_GELNOLEN, rname, DB_NAMELEN, tname);
    return name;
}

/*
 * Get part of a sequence from its `pos' base for `width' bases
 * Bases number from 0?
 */
void DBgetSequence(EdStruct *xx, int seq, int pos, int width, char *str)
{
    char *src;
    int length = DB_Length(xx,seq);
    int i;
    
    src = DBgetSeq(DBI(xx),seq);
    
    /* Lefthand cut off */
    if (pos<0) {
	i = (width<-pos)?width:-pos;
	getLCut(xx,seq, -pos, i, str);
    } else
	i=0;
    
    /*copy sequence*/
    for (; i<width && (pos+i)<length; i++) {
	str[i]=src[pos+i]; 
    }
    
    /* Righthand cut off */
    if (i < width) {
	getRCut(xx,seq, pos+i-length, width-i, &str[i]);
    }
    
    str[width]='\0';
}

/*
 *----------------------------------------------------------------------------
 * Database manipulation
 *---------------------------------------------------------------------------
 */


/*
 * Callback function from the contig registration scheme
 */
/*ARGSUSED*/
void DBi_reg(GapIO *io, int contig, void *fdata, reg_data *jdata) {
    DBInfo *db = (DBInfo *)fdata;
    static char params[100];

    switch(jdata->job) {
    case REG_QUERY_NAME:
	{
	    if (_DBI_order(db)) {
		sprintf(jdata->name.line, "Contig editor @ %d",
			_DB_Number(db, _DBI_order(db)[1]));
	    } else {
		sprintf(jdata->name.line, "Contig editor @ =%d",
			_DBI_contigNum(db));
	    }

	    return;
	}

    case REG_LENGTH:
	{
	    if (_DBI_flags(db) & DB_NO_REGS)
		return;

	    if (_editsMade(db)) {
		verror(ERR_FATAL, "contig_editor",
		       "Cannot update as data is unsaved, yet changed\n");
		return;
	    }

	    /* printf("New length is %d\n", jdata->length.length); */
	    
	    /* Reinitialise data */
	    contig_deregister(io, _DBI_contigNum(db), DBi_reg, db);
	    db->registration_id = -db->registration_id;
	    DBI_callback(db, DBCALL_REINIT, 0, 0, NULL);

	    break;
	}

    case REG_JOIN_TO:
	{
	    int id;

	    if (_editsMade(db)) {
		verror(ERR_FATAL, "contig_editor",
		       "Cannot update as data is unsaved, yet changed\n");
		return;
	    }

	    /* printf("Joining %d, to new contig %d at %d\n",
		   _DBI_contigNum(db), jdata->join.contig, jdata->join.offset);
	     */

	    /*
	     * We must inform the shift of display position before the
	     * change of contig number/length as we'd end up adversely
	     * affecting our 'left' contig in the join.
	     *
	     * The contig_deregister here is because we're joining to an
	     * already registered contig (hopefully). If it doesn't then it
	     * implies a REG_JOIN_TO has been done by something other than
	     * the join editor. This is _currently_ impossible.
	     *
	     * The DBInfo structure for what we're joining to should exist
	     * at this stage (we check to make sure though). We can free
	     * our current DBInfo and update the edStructs to point to
	     * the new value.
	     */
	    contig_deregister(io, _DBI_contigNum(db), DBi_reg, db);

	    DBI_callback(db, DBCALL_JOIN_SHIFT, 0, jdata->join.offset, NULL);

	    /* Find new DBInfo */
	    if (id = type_to_result(io, REG_TYPE_EDITOR, jdata->join.contig)) {
		tman_handle_join(db, (DBInfo *)result_data(io, id,
				     jdata->join.contig));
		DBI_callback(db, DBCALL_RELINK, 0, 0,
			     result_data(io, id, jdata->join.contig));
	    }

	    /* At this point db is now invalid (it's been freed). */
	    
	    break;
	}

    case REG_GET_LOCK:
	{
	    int x, unsaved = 0;

	    /*
	     * Once we've modified data, we need exclusive access, so
	     * clear any write lock.
	     */
	    if (jdata->glock.lock & REG_LOCK_WRITE) {
		if (_editsMade(db)) {
		    jdata->glock.lock &= ~REG_LOCK_WRITE;
		    break;
		}

		/* Find edstruct and check if it's a join editor */
		for (x = 0; x < MAX_DISP_PROCS; x++) {
		    if (_DBI_dispFunc(db)[x] == db_callback_tk) {
		        EdStruct *xx = (EdStruct *)(_DBI_dispData(db)[x]);
			if (xx->link &&
			    (DBI_edits_made(xx->link->xx[0]) ||
			     DBI_edits_made(xx->link->xx[1]))) {
			    unsaved = 1;
			}
		    }
		}

		if (unsaved)
		    jdata->glock.lock &= ~REG_LOCK_WRITE;
	    }

	    break;
	}

    case REG_SET_LOCK:
	{
	    /*
	     * We need exclusive access, so if someone is locking for
	     * write this should only occur when we haven't made
	     * edits. Otherwise scream!
	     */
	    if (jdata->slock.lock & REG_LOCK_WRITE) {
		if (!_editsMade(db))
		    DBI_callback(db, DBCALL_QUIT, 0, 0, NULL);
		else
		    verror(ERR_FATAL, "contig_editor", "HELP - Lock ignored!");
	    }

	    break;
	}

    case REG_QUIT:
	{
	    /*
	     * We are being asked to quit. We can only allow this is we
	     * haven't made changes.
	     */
	    if (_editsMade(db)) {
		jdata->glock.lock &= ~REG_LOCK_WRITE;
	    } else {
	        int x, unsaved = 0;
	        /*
		 * This is icky - if we've got another editor window
		 * linked to this as a join editor and this editor
		 * is not saved, then we also clear the write lock.
		 */
		for (x = 0; x < MAX_DISP_PROCS; x++) {
		    if (_DBI_dispFunc(db)[x] == db_callback_tk) {
		        EdStruct *xx = (EdStruct *)(_DBI_dispData(db)[x]);
			if (xx->link &&
			    (DBI_edits_made(xx->link->xx[0]) ||
			     DBI_edits_made(xx->link->xx[1]))) {
			    unsaved = 1;
			}
		    }
		}

		if (unsaved)
		    jdata->glock.lock &= ~REG_LOCK_WRITE;
		else
		    DBI_callback(db, DBCALL_QUIT, 0, 0, NULL);
	    }

	    break;
	}

    case REG_GET_OPS:
	{
	    /* jdata->get_ops.ops = "Information\0"; */
	    break;
	}

    case REG_PARAMS:
	{
	    sprintf(params, "Contig: %d",
		    _DB_Number(db, _DBI_order(db)[1]));
	    jdata->params.string = params;
	    break;
	}

    case REG_NUMBER_CHANGE:
	{
	    _DBI_contigNum(db) = jdata->number.number;
	    break;
	}

    case REG_GENERIC:
	{
#if 0
	    /*
	     * An example of translating an editior id to and EdStruct
	     * from within the registration scheme.
	     */
	    EdStruct *xx = DBI_to_EdStruct(db, *((int *)jdata->generic.data));

	    if (!xx)
		break;

	    switch(jdata->generic.task) {
	    case 0:
		edCursorLeft(xx);
		break;
	    case 1:
		edCursorRight(xx);
		break;
	    case 2:
		edCursorUp(xx);
		break;
	    case 3:
		edCursorDown(xx);
		break;
	    }
#endif

	    switch(jdata->generic.task) {

	    case TASK_EDITOR_GETCON: {
		/* Find edstruct */
		EdStruct *xx = NULL;
		int x;
		task_editor_getcon *tc;

		for (x = 0; x < MAX_DISP_PROCS; x++) {
		    if (_DBI_dispFunc(db)[x] == db_callback_tk) {
		        xx = (EdStruct *)(_DBI_dispData(db)[x]);
			break;
		    }
		}
		if (!xx)
		    break;

		/*
		 * Returns con == NULL for failure.
		 * Otherwise returns the excerpt of consensus between lreg
		 * and rreg. (If lreg and rreg are both 0 then this is all
		 * of the contig.)
		 */
		tc = (task_editor_getcon *)jdata->generic.data;

		if (!tc->lreg && !tc->rreg) {
		    tc->lreg = 1;
		    tc->rreg = _DB_Length(db, 0);
		}

		if (NULL == (tc->con = (char *)xmalloc(tc->rreg-tc->lreg + 2)))
		    return;

		calc_consensus(0, tc->lreg, tc->rreg, CON_SUM, tc->con, NULL,
			       NULL, NULL, tc->con_cut, tc->qual_cut,
			       contEd_info, (void *)xx);
		tc->con[tc->rreg] = 0;
	    }
	    }
	}

    case REG_HIGHLIGHT_READ:
	{
	    int seq, flags, oflags;

	    if ((seq = NumberToSeq(db, jdata->highlight.seq)) != -1) {
		oflags = flags = _DB_Flags(db, seq);
		if (jdata->highlight.val) {
		    flags |= DB_FLAG_SELECTED;
		} else {
		    flags &= ~DB_FLAG_SELECTED;
		}
		_DBsetFlags(db, seq, flags);
		if (flags != oflags) {
		    int i;
		    for (i = 0; i < MAX_DISP_PROCS; i++) {
			if (_DBI_dispFunc(db)[i])
			    RedisplayName(((EdStruct *)_DBI_dispData(db)[i]),
					  seq);
		    }
		    redisplayDBSequences(db, 1);
		}
	    }
	    break;
	}

    case REG_REGISTER:
	{
	    /*
	     * We're listening for new registrations. When we see one we
	     * send out our current cursor coordinates for each display
	     * [(0,0) expands to these]. The purpose of this is, for example,
	     * to enable the template display to start up and to then
	     * immediately get told by the editor where it's cursor
	     * position is. It's easier to do this than to code the opposite
	     * scheme of the template display sending a "where's the cursor?"
	     * request to each editor in turn.
	     */
	    DBI_callback(db, DBCALL_CURSOR_NOTIFY, 0, 0, NULL);
	    break;
	}

    case REG_CURSOR_NOTIFY:
	{
	    int i, seq, pos;

	    seq = jdata->cursor_notify.cursor->seq;
	    if (seq != -1 && seq != 0) {
		for (i = _DBI_gelCount(db); i > 0; i--) {
		    if (_DB_Number(db, i) == seq) {
			break;
		    }
		}
		seq = i;
		pos = jdata->cursor_notify.cursor->pos;
	    } else {
		if (seq == -1) {
		    seq = 0;
		    pos = jdata->cursor_notify.cursor->abspos;
		} else {
		    pos = jdata->cursor_notify.cursor->pos;
		}
	    }
	    
	    /* Redisplay editor with this cursor id */
	    for (i = 0; i < MAX_DISP_PROCS; i++) {
		if (_DBI_dispFunc(db)[i] &&
		    ((EdStruct *)_DBI_dispData(db)[i])->cursor ==
		    jdata->cursor_notify.cursor) {
		    EdStruct *xx = (EdStruct *)_DBI_dispData(db)[i];
		    if (xx->cursorSeq != seq || xx->cursorPos != pos) {
			xx->cursorSeq = seq;
			xx->cursorPos = pos;
			updateJoinCursor(xx, seq, pos);
			xx->refresh_flags |= ED_DISP_CURSOR;
			redisplayWithCursor(xx);
			repositionTraces(xx);
		    }
		}
	    }
	    break;
	}
    }

    return;
}

/*
 * Parses the REFS notes for sequence 'seqnum', starting with note
 * number 'notenum'.
 */
static void edParseNote(GapIO *io, EdStruct *xx, int seqnum, int notenum) {
    GNotes n;
    int refs = str2type("REFS");
    int reft = str2type("REFT");
    char *text;

    /* Look for a REFS or REFT note */
    for (; notenum; notenum = n.next) {
	note_read(io, notenum, n);
	if ((n.type != refs && n.type != reft) || n.annotation == 0) {
	    continue;
	}

	/* Found one, so now parse it */

	/* REFS */
	if (NULL == (text = TextAllocRead(io, n.annotation)))
	    continue;

	if (n.type == refs) {
	    int offset;
	    int length;

	    if (sscanf(text, "sequence %d %d", &offset, &length) != 2) {
		length = 0;
		if (sscanf(text, "sequence %d", &offset) != 1) {
		    verror(ERR_WARN, "contig_editor", "Invalid REFS note");
		    xfree(text);
		    return;
		}
	    }
	    DBaddFlags(xx, seqnum, DB_FLAG_REFSEQ);
	    DBI(xx)->reference_seq = seqnum;
	    DBI(xx)->reference_len = length;
	    DBI(xx)->reference_offset = offset;
	}

	/* REFT */
	if (n.type == reft) {
	    if (strncmp(text, "control -ve", 11) == 0) {
		DBclearFlags(xx, seqnum, DB_FLAG_REFTRACE);
		DBaddFlags(xx, seqnum, DB_FLAG_REFTRACE_NEG);
	    } else if (strncmp(text, "control +ve", 11) == 0) {
		DBclearFlags(xx, seqnum, DB_FLAG_REFTRACE);
		DBaddFlags(xx, seqnum, DB_FLAG_REFTRACE_POS);
	    } else {
		verror(ERR_WARN, "contig_editor", "Invalid REFT note");
	    }
	}
    
	xfree(text);
    }

    return;
}

/*
 * Create an internal database and read all relevant data into it
 */
int initialiseDB(/* FORIO */
		 EdStruct *xx,
		 GapIO *io,
		 int cnum, 	/* contig number */
		 int idbsiz,	/* size of database */
		 int llino	/* left-most gel in contig */
		 )
{
    int i,c;
    
    DBI_contigNum(xx) = cnum;
    DBI_io(xx) = io;

    /*
     * Register data
     */
    if (!DBI_registration_id(xx))
	DBI_registration_id(xx) = register_id();
    contig_register(io, cnum, DBi_reg, DBI(xx), DBI_registration_id(xx),
		    REG_REQUIRED | REG_DATA_CHANGE | REG_LOCKS |
		    REG_NUMBER_CHANGE | REG_GENERIC | REG_HIGHLIGHT_READ |
		    REG_FLAG_INVIS | REG_REGISTER |
		    REG_CURSOR_NOTIFY,
		    REG_TYPE_EDITOR);
    
    /* Free old memory */
    if (DBI_DB(xx) != NULL) {
	for (i=0; i <= DBI_gelCount(xx); i++) {
	    if (DB_Name(xx, i))
		xfree(DB_Name(xx,i));
	    if (DB_Seq (xx, i))
		xfree(DB_Seq (xx,i));
	    if (DB_Conf(xx, i))
		xfree(DB_Conf(xx,i));
	    if (DB_Opos(xx, i))
		xfree(DB_Opos(xx,i));
	    destroyTagList(DB_Tags(xx,i));
	}
	
	xfree(DBI_DB(xx));
    }
    if (DBI_list(xx))  xfree(DBI_list(xx));
    if (DBI_order(xx)) xfree(DBI_order(xx));
    
    /*
     * count number of gels in contig
     */
    for (DBI_gelCount(xx)=1,i=llino;
	 DBI_gelCount(xx) < idbsiz && (int)io_rnbr(io,i);
	 DBI_gelCount(xx)++,i=(int)io_rnbr(io,i))
	;
    
    if ((DBI_DB(xx) = (DBStruct *)xcalloc(DBI_gelCount(xx)+1,
					  sizeof(DBStruct)))==NULL)
	goto disaster;
    
    if ((DBI_list(xx) =  (int *)xcalloc(DBI_gelCount(xx) + 1,
					sizeof(int)))==NULL)
	goto disaster;
    
    if ((DBI_order(xx) = (int *)xcalloc(DBI_gelCount(xx)+1,
					sizeof(int)))==NULL)
	goto disaster;
    
    /*
     * read information into local database
     */
    for (c=1,i=llino; c < idbsiz && i; c++,i=(int)io_rnbr(io,i)) {
	GReadings r;
	
	DBsetRelPos(xx,c,io_relpos(io,i));
	DBsetLength(xx,c,abs(io_length(io,i)));
	DBsetComp(xx,c,(io_length(io,i)<0) ? -1 : 1);
	DBsetNumber(xx,c,i);

	gel_read(io, i, r);

	if (r.chemistry & GAP_CHEM_TERMINATOR)
	    DBsetFlags(xx,c,DB_FLAG_TERMINATOR);
	else
	    DBsetFlags(xx,c,DB_FLAG_NONE);
	
	DBI_order(xx)[c] = c;
	
	/* See if it has a reference note */
	edParseNote(io, xx, c, r.notes);

	if (DBI_flags(xx) & DB_STORAGE_INTERNAL) {
	    if (DBgetSeq(DBI(xx),c)==NULL)
		goto disaster;
	    if (DBgetName(DBI(xx),c)==NULL)
		goto disaster;
	    (void)DBgetTags(DBI(xx),c);
	}
    }

    if (DBI(xx)->reference_seq < 0) {
	verror(ERR_WARN, "contig_editor",
	       "Reference seq listed in REFS note is not in this contig");
	DBI(xx)->reference_seq = 0;
    }

    /* Clear selected tag - it's been freed */
    xx->select_tag = NULL;
    
    /*
     * Set up consensus
     */
    {
        DBsetRelPos(xx,0,1);
	DBsetComp(xx,0,UNCOMPLEMENTED);
	calculateConsensusLength(xx);
        if ((DBsetSeq(xx,0,(char *)xmalloc(MAX_DISPLAY_WIDTH+1)))==NULL)
	    goto disaster;
	if ((DBsetName(xx,0,(char *)xmalloc(sizeof(char)*(NAMELEN+1))))==NULL)
	    goto disaster;
	sprintf(DB_Name(xx,0),"%*s %-*s",
		DB_GELNOLEN," ",
		DB_NAMELEN, "CONSENSUS");
        DBI_order(xx)[0] = 0;
	DBsetNumber(xx, 0, -DBI_contigNum(xx));
    }

    /*
     * Compute template information.
     */
    if (xx->DBi->templates)
	uninit_template_checks(DBI_io(xx), xx->DBi->templates);
    xx->DBi->templates = init_template_checks(DBI_io(xx), 1, &cnum, 1);
    template_check_set_flags(DBI_io(xx), xx->DBi->templates,
			     TEMP_OFLAG_INTERDIST, 0);
    /* xx->DBi->templates = init_template_checks(DBI_io(xx), 0, NULL, 0);*/
    check_all_templates(DBI_io(xx), xx->DBi->templates);

    /*
     * Redisplay. Important as we've possibly changed the ordering.
     */
    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 0);
    

    return 0;
    
 disaster:
    
    freeDB(xx, 1);
    return 1;
}

/*
 * Free an internal database
 */
void freeDB(EdStruct *xx, int delete_ed)
{
    int i, j, inuse = 0, this_one = -1;
    
    /* Check how many times this contig is in use */
    for (j=0; j<MAXEDSTATES; j++) {
	if (edused[j] && edstate[j].DBi && DBI_DB(&edstate[j]) &&
	    DBI_DB(&edstate[j]) == DBI_DB(xx))
	    inuse++;
    }

    /* Find it in the DBI_dispData array */
    for (j = 0; j < MAX_DISP_PROCS; j++)
	if (DBI_dispData(xx)[j] == xx)
	    this_one = j;

    /* Remove display function and data from this DBi */
    if (this_one != -1) {
	for (i = this_one; i < MAX_DISP_PROCS-1; i++) {
	    DBI_dispFunc(xx)[i] = DBI_dispFunc(xx)[i+1];
	    DBI_dispData(xx)[i] = DBI_dispData(xx)[i+1];
	}
	DBI_dispFunc(xx)[i] = NULL;
	DBI_dispData(xx)[i] = NULL;
	
	DBI_nextDisp(xx)--;
    }
    
    if (inuse < 2) {
	/* Deregister */
	contig_deregister(DBI_io(xx), DBI_contigNum(xx), DBi_reg, DBI(xx));

	if (DBI_DB(xx) != NULL) {
	    for (i=0; i <= DBI_gelCount(xx); i++) {
		if (DB_Name(xx, i))
		    xfree(DB_Name(xx,i));
		if (DB_Seq (xx, i))
		    xfree(DB_Seq (xx,i));
		if (DB_Conf(xx, i))
		    xfree(DB_Conf(xx,i));
		if (DB_Opos(xx, i))
		    xfree(DB_Opos(xx,i));
		destroyTagList(DB_Tags(xx,i));
	    }
	
	    xfree(DBI_DB(xx));
	}

	xfree(DBI_list(xx));
	xfree(DBI_order(xx));

	DBI_DB(xx) = NULL;
	DBI_list(xx) = NULL;
	DBI_order(xx) = NULL;

	destroyFreeTagList();

	xfree(DBI(xx));
    }

    if (delete_ed)
	destroyEdStruct(xx);
}



/*
 * Writes reference notes
 */
static int writeRefNotes(EdStruct *xx, int seq) {
    GReadings r;
    GNotes n;
    GapIO *io = DBI_io(xx);
    int nnote;
    int refs_type = str2type("REFS");
    int reft_type = str2type("REFT");

    /* Find existing REF[ST] notes and delete them */
    gel_read(io, DB_Number(xx, seq), r);
    for (nnote = r.notes; nnote; nnote = n.next) {
	note_read(io, nnote, n);
	if ((n.type != refs_type && n.type != reft_type))
	    continue;

	delete_note(io, nnote);
    }

    if (DB_Flags(xx, seq) & DB_FLAG_REFSEQ) {
	/* Sequence is reference sequence => create REFS note */
	nnote = new_note(io, refs_type, GT_Readings, DB_Number(xx, seq));
	if (nnote) {
	    char comment[1024];
	    if (DBI(xx)->reference_len)
		sprintf(comment, "sequence %d %d",
			DBI(xx)->reference_offset,
			DBI(xx)->reference_len);
	    else
		sprintf(comment, "sequence %d",
			DBI(xx)->reference_offset);
	    edit_note(io, nnote, NULL, comment);
	}
    }

    if (DB_Flags(xx, seq) & DB_FLAG_REFTRACE) {
	/* Sequence is a reference trace => create REFT note */
	nnote = new_note(io, reft_type, GT_Readings, DB_Number(xx, seq));
	if (nnote) {
	    if (DB_Flags(xx, seq) & DB_FLAG_REFTRACE_NEG)
		edit_note(io, nnote, NULL, "control -ve");
	    else if (DB_Flags(xx, seq) & DB_FLAG_REFTRACE_POS)
		edit_note(io, nnote, NULL, "control +ve");
	    else
		verror(ERR_WARN, "writeRefNotes",
		       "Unknown reference trace type");
	}
    }
    
    return 0;
}

/*
 * Save an internal database
 */
void saveDB(EdStruct *xx, GapIO *io, int auto_save, int notify) {
    int i;
    int_f N,leftN,rightN,cn;
    int flag;
    reg_length jl;

    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return;
    }

    for (i=1; i<=DBI_gelCount(xx); i++) {
	flag = DB_Flags(xx,DBI_order(xx)[i]);
	N = DB_Number(xx,DBI_order(xx)[i]);

	/* if (flag & DB_FLAG_REL_MODIFIED) { */

	/*
	 * DB_FLAG_REL_MODIFIED is set when the left-right neighbours of
	 * this reading change. However the code to set this at present does
	 * not reset DB_FLAG_REL_MODIFIED for its two neighbours.
	 */
	if (1) {
	    /*
	     * update relationships
	     */
	    io_relpos(io,N) = DB_RelPos(xx,DBI_order(xx)[i]);
	    io_length(io,N) = (DB_Comp(xx,DBI_order(xx)[i]) == COMPLEMENTED)
		? -DB_Length(xx,DBI_order(xx)[i])
		    : DB_Length(xx,DBI_order(xx)[i]);

	    if (i==1)
		leftN = 0;
	    else
		leftN = DB_Number(xx,DBI_order(xx)[i-1]);

	    if (i==DBI_gelCount(xx))
		rightN = 0;
	    else
		rightN = DB_Number(xx,DBI_order(xx)[i+1]);
	    
	    io_lnbr(io,N) = leftN;
	    io_rnbr(io,N) = rightN;

	    /* update gel line */ /* FORIO */
	    writeg_(handle_io(io),&N,&io_relpos(io,N),&io_length(io,N),
		    &io_lnbr(io,N),&io_rnbr(io,N));
	}

	/*
	 * update working versions
	 */
	if (flag&DB_FLAG_SEQ_IN_MEMORY && flag&DB_FLAG_SEQ_MODIFIED) {
	    int tmp_len = DBI_DB(xx)[DBI_order(xx)[i]].gap_length;

	    (void) io_write_seq(io,
				(int) N,
				&tmp_len,
				&DBI_DB(xx)[DBI_order(xx)[i]].gap_start,
				&DBI_DB(xx)[DBI_order(xx)[i]].gap_end,
				DBI_DB(xx)[DBI_order(xx)[i]].sequence,
				DBI_DB(xx)[DBI_order(xx)[i]].gap_conf,
				DBI_DB(xx)[DBI_order(xx)[i]].gap_opos
				);
	}

	/*
	 * Update reference notes
	 */
	if (flag & DB_FLAG_NOTE_MODIFIED) {
	    writeRefNotes(xx, DBI_order(xx)[i]);
	}

	/*
	 * update tag list
	 */
	if (flag&DB_FLAG_TAG_IN_MEMORY && flag&DB_FLAG_TAG_MODIFIED) {
            writeTagList(xx, DBI_order(xx)[i]);
        }
	
        /*
	 * Clear sequence-edited flags
	 */
	DBsetFlags(xx,DBI_order(xx)[i],
		   flag & ~(DB_FLAG_SEQ_MODIFIED |
			    DB_FLAG_REL_MODIFIED |
			    DB_FLAG_TAG_MODIFIED));
    }
    
    /*
     * update contig relationships
     */
    calculateConsensusLength(xx);
    cn = DBI_contigNum(xx);
    io_clength(io,cn) = DB_Length(xx,0);
    io_clnbr  (io,cn) = DB_Number(xx,DBI_order(xx)[1]);
    io_crnbr  (io,cn) = DB_Number(xx, DBI_order(xx)[DBI_gelCount(xx)]);
    /* update contig line */

    /* FORIO */
    writec_(handle_io(io),
	    &cn,
	    &io_clength(io,cn),
	    &io_clnbr  (io,cn),
	    &io_crnbr  (io,cn)
	    );

    /*
     * Update consensus tag list
     */
    flag = DB_Flags(xx, 0);
    if (flag & DB_FLAG_TAG_IN_MEMORY && flag & DB_FLAG_TAG_MODIFIED) {
	writeTagList(xx, 0);
    }

    if (auto_save)
	resetEdits(xx);
    else
	freeAllUndoLists(xx);
    
    /* FORIO */
    flush2t(io);

    /* Send notification to other registered functions */
    if (notify) {
	/* DBI_flags(xx) |= DB_NO_REGS; - search for NO_REGS in change.log */
	jl.job = REG_LENGTH;
	jl.length = DB_Length(xx, 0);
	contig_notify(io, cn, (reg_data *)&jl);
	/* DBI_flags(xx) &= ~DB_NO_REGS; */
    }
}

/*
 *----------------------------------------------------------------------------
 * Obtaining region/position information
 *---------------------------------------------------------------------------
 */

/*
 * Return number of sequences on screen
 */
int linesInRegion(EdStruct *xx, int pos, int width)
{
    int i, count;
    int *sets = NULL;

    sets = (int *)xcalloc(xx->nsets+1, sizeof(int));
    
    /* i = posToIndex(xx,pos - dbi_max_gel_len(DBI(xx),0)); */
    i = 1;

    if (xx->reveal_cutoffs) {
	for (count=0 ; i && i <= DBI_gelCount(xx); i++) {
	    int p;
	    int snum;

	    p = DB_RelPos(xx, DBI_order(xx)[i]) - DB_Start(xx, DBI_order(xx)[i]);
	    snum = xx->set ? xx->set[DBI_order(xx)[i]] : 0;

	    /*
	    if (p >= pos + width + MAX_CUTOFF)
		break;
	    */

	    if (p + DB_Length2(xx,DBI_order(xx)[i]) > pos &&
		p < pos + width && DB_Length(xx, DBI_order(xx)[i]) &&
		(!xx->set || xx->curr_set == 0 || snum == xx->curr_set)) {
		if (xx->set_collapsed && xx->set_collapsed[snum] && sets[snum])
		    continue;
		sets[snum]++;
		count++;
	    }
	}
    } else {
	for (count=0 ;
	     i &&
	     i <= DBI_gelCount(xx) &&
	     DB_RelPos(xx,DBI_order(xx)[i]) < (pos+width) ;
	     i++) {
	    int snum = xx->set ? xx->set[DBI_order(xx)[i]] : 0;
	    if (DB_RelPos(xx,DBI_order(xx)[i]) +
		     DB_Length(xx,DBI_order(xx)[i]) > pos &&
		DB_Length(xx,DBI_order(xx)[i]) &&
		(!xx->set || xx->curr_set == 0 || snum == xx->curr_set)) {
		if (xx->set_collapsed && xx->set_collapsed[snum] && sets[snum])
		    continue;
		sets[snum]++;
		count++;
	    }
	}
    }

    count += xx->consensusDisplayed;
    xfree(sets);

    return count;
}


/*
 * Return number of sequences on screen
 */
int linesOnScreen (EdStruct *xx, int pos, int width)
{
    int i, count;
    int *sets = NULL;
    sets = (int *)xcalloc(xx->nsets+1, sizeof(int));
    
    /* i = posToIndex(xx,pos - dbi_max_gel_len(DBI(xx),0) - MAX_CUTOFF); */
    i = 1;
    for (count=0; i && i<=DBI_gelCount(xx); i++) {
	int relPos, length;
	int len_lcut, len_rcut;
	int snum;
	
	if (xx->reveal_cutoffs) {
	    len_lcut = lenLCut(xx,DBI_order(xx)[i]);
	    len_rcut = lenRCut(xx,DBI_order(xx)[i]);
	} else {
	    if (DB_RelPos(xx, DBI_order(xx)[i]) > pos + width)
		break;
	    len_lcut = len_rcut = 0;
	}

	relPos = DB_RelPos(xx,DBI_order(xx)[i]) - len_lcut;
	length = DB_Length(xx,DBI_order(xx)[i]) + len_lcut + len_rcut;

	snum = xx->set ? xx->set[DBI_order(xx)[i]] : 0;

        if (relPos < pos+width && relPos+length>pos) {
	    if (!xx->set || xx->curr_set == 0 || snum == xx->curr_set) {
		if (xx->set_collapsed && xx->set_collapsed[snum] && sets[snum])
		    continue;
		sets[snum]++;
		count++;
	    }
	}
    }

    count += xx->consensusDisplayed;
    xfree(sets);
    
    return count;
}


/*
 * Sorts a list of sequences so that sequences with the same template
 * are next to each other.
 * We use bubble sort here, despite inefficiencies, as it is a nice and
 * easy "stable" sort. This allows us to firstly sort by strand so that
 * the sequences within a template are always shown in +/- order.
 */
static void sort_seq_by_template(EdStruct *xx, int *list, int count) {
    int i, swaps;
    
    /* Trivial bubble sort */
    do {
	swaps = 0;
	for (i = 0; i < count-1; i++) {
	    if (DBI_DB(xx)[list[i]].template < DBI_DB(xx)[list[i+1]].template)
	    {
		int t;
		t = list[i];
		list[i] = list[i+1];
		list[i+1] = t;
		swaps = 1;
	    }
	}
    } while (swaps);
}

/*
 * Sorts a list of sequences so that sequences with the same grouping
 * are next to each other.
 * We use bubble sort here, despite inefficiencies, as it is a nice and
 * easy "stable" sort. This allows us to firstly sort by strand so that
 * the sequences within a set are also sorted by template and by +/-
 * order.
 */
static void sort_seq_by_set(EdStruct *xx, int *list, int count) {
    int i, swaps;

    if (!xx->set)
	return;


    /* Trivial bubble sort */
    do {
	swaps = 0;
	for (i = 0; i < count-1; i++) {
	    if (xx->set[list[i]] > xx->set[list[i+1]])
	    {
		int t;
		t = list[i];
		list[i] = list[i+1];
		list[i+1] = t;
		swaps = 1;
	    }
	}
    } while (swaps);
}


/*
 * Sort sequence list by strand using a stable sort.
 */
static void sort_seq_by_strand(EdStruct *xx, int *list, int count) {
    int i, swaps;

    /* Trivial bubble sort */
    do {
	swaps = 0;
	for (i = 0; i < count-1; i++) {
	    if (DB_Comp(xx, list[i])   == COMPLEMENTED &&
		DB_Comp(xx, list[i+1]) != COMPLEMENTED) {
		int t;
		t = list[i];
		list[i] = list[i+1];
		list[i+1] = t;
		swaps = 1;
	    }
	}
    } while (swaps);
}


/*
 * Sort sequence list by clone using a stable sort.
 * Not desparately efficient as the clone information is not cached, instead
 * it is read from the database on-the-fly.
 */
static void sort_seq_by_clone(EdStruct *xx, int *list, int count) {
    int i, swaps;
    GapIO *io = DBI_io(xx);
    char (*clones)[DB_NAMELEN+1];
    char **clones2;

    /* Precache clone names */
    clones = (char (*)[DB_NAMELEN+1])xmalloc(count * sizeof(*clones));
    clones2 = (char **)xmalloc(count * sizeof(*clones2));
    if (!clones || !clones2)
	return;

    for (i = 0; i < count; i++) {
	GReadings r;
	GTemplates t;
	GClones c;
	
	/* clones2 holds ptrs to clones, which we swap during sorting */
	clones2[i] = clones[i];

	clones[i][0] = '\0';
	gel_read(io, DB_Number(xx, list[i]), r);
	if (!r.template)
	    continue;

	template_read(io, r.template, t);
	if (!t.clone)
	    continue;

	clone_read(io, t.clone, c);
	if (!c.name)
	    continue;

	TextRead(io, c.name, clones[i], DB_NAMELEN);
	clones[i][DB_NAMELEN] = 0;
    }

    /* Trivial bubble sort */
    do {
	swaps = 0;
	for (i = 0; i < count-1; i++) {
	    if (strcmp(clones2[i], clones2[i+1]) < 0) {
		int t;
		char *c;

		t = list[i];
		list[i] = list[i+1];
		list[i+1] = t;

		c = clones2[i];
		clones2[i] = clones2[i+1];
		clones2[i+1] = c;

		swaps = 1;
	    }
	}
    } while (swaps);

    xfree(clones);
    xfree(clones2);
}


/*
 * Sort functions by sequence name and number - these do not need to use
 * a stable sort as the names and numbers are unique.
 */
static EdStruct *tmp_xx; /* Global as qsort can't pass clientdata in */
static int qsort_seq_by_alpha(const void *ps1, const void *ps2) {
    int s1 = *(int *)ps1, s2 = *(int *)ps2;
    return strcmp(DBgetName(DBI(tmp_xx), s1) + DB_GELNOLEN+1,
		  DBgetName(DBI(tmp_xx), s2) + DB_GELNOLEN+1);
}

static int qsort_seq_by_numeric(const void *ps1, const void *ps2) {
    int s1 = *(int *)ps1, s2 = *(int *)ps2;
    return DB_Number(tmp_xx, s1) - DB_Number(tmp_xx, s2);
}

static int qsort_seq_by_chemistry(const void *ps1, const void *ps2) {
    int s1 = *(int *)ps1, s2 = *(int *)ps2;
    GReadings r1, r2;
    if (DB_Number(tmp_xx, s1) <= 0 || DB_Number(tmp_xx, s2) <= 0)
	return 0;

    gel_read(DBI_io(tmp_xx), DB_Number(tmp_xx, s1), r1);
    gel_read(DBI_io(tmp_xx), DB_Number(tmp_xx, s2), r2);
    if (r1.chemistry < r2.chemistry) {
	return -1;
    } else if (r1.chemistry > r2.chemistry) {
	return 1;
    } else {
	return DB_RelPos(tmp_xx, s1) - DB_RelPos(tmp_xx, s2);
    }
}

/*
 * Sort sequence list in-line by xx->group_mode.
 *
 * Assumption: list is generated in POSITIONal sort.
 */
static void sort_seq_list(EdStruct *xx, int *list, int count) {
    static int *last_list = NULL;
    static int *last_sorted = NULL;
    static int last_count = 0;
    static int last_mode  = 0;

    /* Check if it's the same list as before */
    if (last_list && last_sorted &&
	last_count == count &&
	last_mode == xx->group_mode &&
	memcmp(list, last_list, count * sizeof(*list)) == 0) {
	memcpy(list, last_sorted, count * sizeof(*list));
	return;
    }
	       
    /* Cache this input */
    last_list = xrealloc(last_list, count * sizeof(*list));
    if (last_list)
	memcpy(last_list, list, count * sizeof(*list));
    last_count = count;
    last_mode  = xx->group_mode;

    switch(xx->group_mode) {
    case POSITION:
	break; /* Starting point */

    case TEMPLATE:
	sort_seq_by_strand(xx, list, count);
	sort_seq_by_template(xx, list, count);
	break;

    case STRAND:
	sort_seq_by_strand(xx, list, count);
	break;

    case CLONE:
	sort_seq_by_clone(xx, list, count);
	break;

    case ALPHA:
	tmp_xx = xx;
	qsort(list, count, sizeof(*list), qsort_seq_by_alpha);
	break;

    case NUMERIC:
	tmp_xx = xx;
	qsort(list, count, sizeof(*list), qsort_seq_by_numeric);
	break;

    case CHEMISTRY:
	tmp_xx = xx;
	qsort(list, count, sizeof(*list), qsort_seq_by_chemistry);
	break;

    case SET:
	/*
	sort_seq_by_strand(xx, list, count);
	sort_seq_by_template(xx, list, count);
	*/
	sort_seq_by_set(xx, list, count);
	break;
    }

    /*
     * Move any REFSEQ tagged sequences to the statr of the list.
     * This cannot be done for POSITION collating order as there is code
     * (in searchUtils.c for example) that relies on the position sorting
     * to be absolutely correct for optimisation purposes.
     */
    if (xx->group_mode != POSITION) {
	int i;
	for (i = 1; i < count; i++) {
	    if (DB_Flags(xx, list[i]) & DB_FLAG_REFSEQ) {
		/* move to first pos */
		int rs = list[i];
		memmove(&list[1], &list[0], i*sizeof(*list));
		list[0] = rs;
	    }
	}
    }

    /* Cache this output */
    last_sorted = xrealloc(last_sorted, count * sizeof(*list));
    if (last_sorted)
	memcpy(last_sorted, list, count * sizeof(*list));
}


/*
 * Return a pointer to list of sequences in region of contig
 */
int *sequencesInRegion(EdStruct *xx,int pos, int width)
{
    int i, count;
    int *sets = (int *)xcalloc(xx->nsets+1, sizeof(int));
    
    /* i = posToIndex(xx,pos - dbi_max_gel_len(DBI(xx),0)); */
    i = 1;

    if (xx->reveal_cutoffs) {
	for (count=0 ; i && i <= DBI_gelCount(xx); i++) {
	    int p;
	    int snum = xx->set ? xx->set[DBI_order(xx)[i]] : 0;

	    p = DB_RelPos(xx, DBI_order(xx)[i]) - DB_Start(xx, DBI_order(xx)[i]);
	    /*
	    if (p >= pos + width + MAX_CUTOFF)
		break;
	    */

	    if (p + DB_Length2(xx,DBI_order(xx)[i]) > pos &&
		p < pos + width && DB_Length(xx, DBI_order(xx)[i]) &&
		(!xx->set || xx->curr_set == 0 || snum == xx->curr_set)) {
		if (xx->set_collapsed && xx->set_collapsed[snum] && sets[snum])
		    continue;
		DBI_list(xx)[count++]=DBI_order(xx)[i];
		sets[snum]++;
	    }
	}
    } else {
	for (count=0 ;
	     i &&
	     i <= DBI_gelCount(xx) &&
	     DB_RelPos(xx, DBI_order(xx)[i]) < (pos+width) ;
	     i++) {
	    int snum = xx->set ? xx->set[DBI_order(xx)[i]] : 0;

	    if (DB_RelPos(xx,DBI_order(xx)[i]) +
		DB_Length(xx,DBI_order(xx)[i]) > pos &&
		DB_Length(xx,DBI_order(xx)[i]) &&
		(!xx->set || xx->curr_set == 0 || snum == xx->curr_set)) {
		if (xx->set_collapsed && xx->set_collapsed[snum] && sets[snum])
		    continue;
		DBI_list(xx)[count++]=DBI_order(xx)[i];
		sets[snum]++;
	    }
	}
    }

    if (xx->group_mode)
	sort_seq_list(xx, DBI_list(xx), count);
    sort_seq_by_set(xx, DBI_list(xx), count);

    if (xx->consensusDisplayed) DBI_list(xx)[count++] = 0;
    xfree(sets);
    
    return DBI_list(xx);
}


/*
 * Return a pointer to list of sequences on screen
 */
int *sequencesOnScreen(EdStruct *xx,int pos, int width)
{
    int i, count;
    int *sets = (int *)xcalloc(xx->nsets+1, sizeof(int));
    
    /* i = posToIndex(xx,pos - dbi_max_gel_len(DBI(xx),0) - MAX_CUTOFF); */
    i = 1;

    for (count=0; i && i<=DBI_gelCount(xx); i++) {

	int relPos, length;
	int len_lcut, len_rcut;
	int snum = xx->set ? xx->set[DBI_order(xx)[i]] : 0;
	
	if (xx->reveal_cutoffs) {
	    len_lcut = lenLCut(xx,DBI_order(xx)[i]);
	    len_rcut = lenRCut(xx,DBI_order(xx)[i]);
	} else {
	    if (DB_RelPos(xx, DBI_order(xx)[i]) > pos + width)
		break;
	    len_lcut = len_rcut = 0;
	}

	relPos = DB_RelPos(xx,DBI_order(xx)[i]) - len_lcut;
	length = DB_Length(xx,DBI_order(xx)[i]) + len_lcut + len_rcut;

        if (relPos < pos+width && relPos+length>pos) {
	    if (!xx->set || xx->curr_set == 0 || snum == xx->curr_set) {
		if (xx->set_collapsed && xx->set_collapsed[snum] && sets[snum])
		    continue;
		sets[snum]++;
		DBI_list(xx)[count++]=DBI_order(xx)[i];
	    }
	}
    }

    if (xx->group_mode)
	sort_seq_list(xx, DBI_list(xx), count);
    sort_seq_by_set(xx, DBI_list(xx), count);

    if (xx->consensusDisplayed) DBI_list(xx)[count++] = 0;
    xfree(sets);
   
    return DBI_list(xx);
}


/*
 * returns relative position in a sequence as an 
 * absolute position in the contig
 */
int positionInContig(EdStruct *xx, int seq, int pos)
{
    return DB_RelPos(xx,seq) + pos - 1;
}


/*
 * returns true if base in `seq' at position `pos' is currently
 * being displayed on screen 
 */
int onScreen (EdStruct *xx, int seq, int pos, int *wrong_x)
{
    int posInContig, visible = 1;
    int *seqList;
    int screenRow;

    posInContig = positionInContig(xx,seq,pos);

    /*
     * Scan through sequences in the displayed Y region to see whether
     * it's visible. The consensus is always visible.
     */
    seqList = sequencesOnScreen(xx,xx->displayPos, xx->displayWidth);
    for(screenRow=xx->displayYPos;
	screenRow < xx->displayHeight/xx->lines_per_seq + xx->displayYPos-2 &&
	seqList[screenRow] != seq;
	screenRow++);
    if (seqList[screenRow] != seq && seq != 0) {
	visible = 0;
    }

    if (wrong_x)
	*wrong_x = !(posInContig >= xx->displayPos &&
		     posInContig < xx->displayPos + xx->displayWidth);

    return (posInContig >= xx->displayPos &&
            posInContig < xx->displayPos + xx->displayWidth &&
	    visible);
}


/*
 * Get maximum extents of sequence, taking into account cutoffs.
 */
void extents(EdStruct *xx, int *left, int *right)
{
    if (xx->reveal_cutoffs) {
	int eleft, eright;
	int maxgel;
	int i;
	
	/*
	 * Determine for left end
	 */
	eleft = 1;
	
	for (i=1; i<=DBI_gelCount(xx); i++) {
	    int thisleft;
	    
	    thisleft = DB_RelPos(xx,DBI_order(xx)[i]) -
		lenLCut(xx,DBI_order(xx)[i]);
	    
	    if (eleft > thisleft)
		eleft = thisleft;
	}
	
	
	eright = DB_Length(xx,0);
	maxgel = dbi_max_gel_len(DBI(xx),0);
	for (i = DBI_gelCount(xx); i >= 1; i--) {
	    int thisright;
	    
	    thisright =  DB_RelPos(xx,DBI_order(xx)[i]) +
		DB_Length(xx,DBI_order(xx)[i]) + lenRCut(xx,DBI_order(xx)[i]) - 1;
	    
	    if (eright < thisright)
		eright = thisright;
	}
	
	*left = eleft;
	*right = eright;
	
    } else {
	*left = 1;
	*right = DB_Length(xx,0);
    }
}




/*
 * get information about relative positions of two joined contigs
 */
static void joinedExtents(EdStruct *xx, int *leftPos, int* rightPos)
{
    int offset = editorLockedPos(xx->link->xx, 0/*don't force recalculation*/);
    EdStruct *otherxx;
    int left,right;
    int otherleft,otherright;
    
    otherxx = xx->link->xx[0];
    extents(xx,&left,&right);

    if (otherxx==xx) {
	otherxx = xx->link->xx[1];
	extents(otherxx,&otherleft,&otherright);
        *leftPos = min(left,otherleft-offset);
	*rightPos = max(right,otherright-offset);
    } else {
	extents(otherxx,&otherleft,&otherright);
        *leftPos = min(left,otherleft+offset);
	*rightPos = max(right,otherright+offset);
    }
}




/*
 * Get maximum extents of sequence allowing for cutoff and join mode.
 */
void getExtents(EdStruct *xx)
{
    if (inJoinMode(xx) && editorLocked(xx))
	joinedExtents(xx, &xx->extent_left, &xx->extent_right);
    else {
	extents(xx, &xx->extent_left, &xx->extent_right);
	if (inJoinMode(xx)  && !editorLocked(xx)) {
	    xx->extent_right += xx->displayWidth - 2;
	    xx->extent_left  -= xx->displayWidth - 1;
	}
    }
}



/*
 *----------------------------------------------------------------------------
 * Scrolling
 *---------------------------------------------------------------------------
 */

/*
 * Increase the leftmost base position on the screen by a symbolic amount
 */
static void incDisplayPosP (EdStruct *xx, int distance)
{
    switch (distance) {
    case D_screen     : xx->displayPos += xx->displayWidth; break;
    case D_halfScreen : xx->displayPos += xx->displayWidth/2; break;
    case D_character  : xx->displayPos += 1; break;
    }

    if (xx->displayPos > xx->extent_right + 2 - xx->displayWidth)
	xx->displayPos = xx->extent_right + 2 - xx->displayWidth;

    xx->refresh_flags |= ED_DISP_SCROLL;
    redisplaySequences(xx, 0);
}


/*
 * Decrease the leftmost base position on the screen by a symbolic ammount
 */
static void decDisplayPosP (EdStruct *xx, int distance)
{
    switch (distance) {
    case D_screen     : xx->displayPos -= xx->displayWidth; break;
    case D_halfScreen : xx->displayPos -= xx->displayWidth/2; break;
    case D_character  : xx->displayPos -= 1; break;
    }

    if (xx->displayPos < xx->extent_left) xx->displayPos = xx->extent_left;

    xx->refresh_flags |= ED_DISP_SCROLL;
    redisplaySequences(xx, 0);
}


/*
 * Increase the leftmost base position on the screen by a symbolic ammount
 */
void incDisplayPos (EdStruct *xx, int distance)
{
    if (editorLocked(xx)) {
	incDisplayPosP(xx->link->xx[0], distance);
	incDisplayPosP(xx->link->xx[1], distance);
    } else
	incDisplayPosP(xx, distance);
    
    redisplayDisagreement(xx);
}


/*
 * Decrease the leftmost base position on the screen by a symbolic ammount
 */
void decDisplayPos (EdStruct *xx, int distance)
{
    if (editorLocked(xx)) {
	decDisplayPosP(xx->link->xx[0], distance);
	decDisplayPosP(xx->link->xx[1], distance);
    } else
	decDisplayPosP(xx, distance);
    
    redisplayDisagreement(xx);
}



/*
 * Adjust displayed region
 */
void setDisplayPosP(EdStruct *xx, int pos)
{
    if (editorLocked(xx)) {
        int offset = editorLockedPos(xx->link->xx, 1/*don't force recalc*/);
	EdStruct *otherxx;

	otherxx = xx->link->xx[0];
	if (otherxx == xx) {
	    otherxx = xx->link->xx[1];
	    otherxx->displayPos = pos + offset;
	} else {
	    otherxx->displayPos = pos - offset;
	}

	ed_set_slider_pos(otherxx, otherxx->displayPos);
	otherxx->refresh_flags |= ED_DISP_SCROLL;
	redisplaySequences(otherxx, 0);
    }
    
    xx->displayPos = pos;
}

/*
 * centralise pos on screen
 */
void setDisplayPos(EdStruct *xx, int pos)
{
    setDisplayPosP(xx, pos);
    decDisplayPos(xx, D_halfScreen);
}


/*
 * Set left of screen to pos
 */
void setDisplayPos2(EdStruct *xx, int pos)
{
    /* Set NO_DIFFS to prevent too many diff-line redraws (flicker) */
    if (xx->link) {
	xx->link->xx[0]->refresh_flags |= ED_DISP_NO_DIFFS;
	xx->link->xx[1]->refresh_flags |= ED_DISP_NO_DIFFS;
    } else {
	xx->refresh_flags |= ED_DISP_NO_DIFFS;
    }
    setDisplayPosP(xx, pos);

    xx->refresh_flags |= ED_DISP_SCROLL;
    redisplaySequences(xx, 0);

    /* Clear NO_DIFFS, and redraw diffs, once only */
    if (xx->link) {
	xx->link->xx[0]->refresh_flags &= ~ED_DISP_NO_DIFFS;
	xx->link->xx[1]->refresh_flags &= ~ED_DISP_NO_DIFFS;
    } else {
	xx->refresh_flags &= ~ED_DISP_NO_DIFFS;
    }
    redisplayDisagreement(xx);
}

/*
 * Set the y display position (vert. scrolling)
 */
void setDisplayYPos(EdStruct *xx, int pos)
{
    xx->displayYPos = pos;
    xx->refresh_flags |= ED_DISP_SCROLL;
    redisplaySequences(xx, 0);
}

/*
 *----------------------------------------------------------------------------
 * Display code
 *---------------------------------------------------------------------------
 */


/*
 * Calls the display callback.
 * If update_all is set the other views of this data are redrawn as well.
 * For efficiency it's often not required to do this (eg when scrolling one
 * view).
 */
void redisplaySequences(EdStruct *xx, int update_all) {
    int i;

    /*
     * We need to copy the redisplay optimisation information from this view
     * to other views.
     */
    if (update_all) {
	int flags = xx->refresh_flags;
	int seq = xx->refresh_seq;

	for (i = 0; i < MAX_DISP_PROCS; i++) {
	    if (DBI_dispFunc(xx)[i]) {
		EdStruct *xx2 = (EdStruct *)(DBI_dispData(xx)[i]);
		
		xx2->refresh_flags |= flags;
		xx2->refresh_seq = seq;
		DBI_dispFunc(xx)[i](DBI_dispData(xx)[i], DBCALL_REDISPLAY,
				    0, 0, NULL);
	    }
	}
    } else {
	for (i = 0; i < MAX_DISP_PROCS; i++) {
	    if (DBI_dispData(xx)[i] == xx)
		DBI_dispFunc(xx)[i](DBI_dispData(xx)[i], DBCALL_REDISPLAY,
				    0, 0, NULL);
	}
    }
}

/*
 * As redisplaySequences except with a db instead and no 'update_all'
 */
void redisplayDBSequences(DBInfo *db, int names_only) {
    int i;

    for (i = 0; i < MAX_DISP_PROCS; i++) {
	if (_DBI_dispFunc(db)[i])
	    _DBI_dispFunc(db)[i](_DBI_dispData(db)[i], DBCALL_REDISPLAY,
				 0, names_only, NULL);
    }
}

/*
 * ensure that the cursor is visible on the screen
 */
void showCursor(EdStruct *xx, int seq, int pos)
{
    int wrong_x;

    if (onScreen(xx, seq, pos, &wrong_x)) {
        positionCursor(xx,seq,pos);
    } else {
	int screenRow;
	int *seqList;

	if (wrong_x)
	    setDisplayPos(xx,positionInContig(xx,seq,pos));

	/*
	 * Scan through the displayed sequences at this point to see
	 * whether the sequence is above or below the visable range.
	 */
	seqList = sequencesOnScreen(xx,xx->displayPos, xx->displayWidth);

	/* Above */
	for (screenRow = 0; screenRow < xx->displayYPos &&
	     seqList[screenRow] != seq;
	     screenRow++);
	if (seqList[screenRow] == seq) {
	    xx->displayYPos = screenRow;

	} else {
	    /* Below */
	    for (screenRow = xx->displayHeight/xx->lines_per_seq
		     + xx->displayYPos - 1;
		 seqList[screenRow] != 0 && seqList[screenRow] != seq;
		 screenRow++);
	    if (seqList[screenRow] != 0) {
		xx->displayYPos =
		    (screenRow - (xx->displayHeight/xx->lines_per_seq - 2));
	    }
	}

	xx->refresh_flags |= ED_DISP_YSCROLL;
	redisplaySequences(xx, 0);
    }
}


/*
 * Redisplay screen, ensuring cursor display
 */
void redisplayWithCursor(EdStruct *xx)
{
    xx->refresh_flags |= ED_DISP_CURSOR;
    if (onScreen(xx, xx->cursorSeq, xx->cursorPos, NULL)) {
	redisplaySequences(xx, 1);
    } else
	showCursor(xx,xx->cursorSeq,xx->cursorPos);
}


/*
 * Keeps the 'opposite' editor in a join editor in cursor sync.
 */
static void updateJoinCursor(EdStruct *xx, int seq, int pos) {
    EdStruct *oxx;
    int pos1, pos2;
    int seq2, tmp;
    int valid;

    /* Skip if we're not in a locked join editor */
    if (!(inJoinMode(xx) && editorLocked(xx) && xx->link))
	return;

    /* Find consensus position in other editor */
    oxx = xx->link->xx[0] == xx ? xx->link->xx[1] : xx->link->xx[0];
    pos1 = positionInContig(xx, xx->cursorSeq, xx->cursorPos);
    seq2 = oxx->cursorSeq;

    if (xx == xx->link->xx[1]) {
	pos2 = pos1 - xx->link->lockOffset;
    } else {
	pos2 = pos1 + xx->link->lockOffset;
    }

    /* Are we still in a valid position for our original sequence? */
    valid = 1;
    tmp = pos2 - DB_RelPos(oxx, seq2) + 1;
    if (oxx->reveal_cutoffs) {
	if (tmp <= -DB_Start(oxx, seq2) ||
	    tmp > DB_Length2(oxx, seq2) - DB_Start(oxx, seq2) + 1)
	    valid = 0;
    } else {
	if (tmp < 1 || tmp > DB_Length(oxx, seq2)+1)
	    valid = 0;
    }

    if (valid) {
	/* Yes... so keep that sequence */
	oxx->cursorPos = tmp;
	oxx->cursorSeq = seq2;
    } else {
	/* Invalid in old sequence, so move to consensus */

	/* Range checking on consensus */
	if (pos2 < 1)
	    pos2 = 1;
	if (pos2 > DB_Length(oxx, 0)+1)
	    pos2 = DB_Length(oxx, 0)+1;

	/* Due to the range checking, we may now be valid again - check! */
	valid = 1;
	tmp = pos2 - DB_RelPos(oxx, seq2) + 1;
	if (oxx->reveal_cutoffs) {
	    if (tmp <= -DB_Start(oxx, seq2) ||
		tmp > DB_Length2(oxx, seq2) - DB_Start(oxx, seq2) + 1)
		valid = 0;
	} else {
	    if (tmp < 1 || tmp > DB_Length(oxx, seq2)+1)
		valid = 0;
	}

	/* And now set the cursor */
	if (valid) {
	    oxx->cursorPos = tmp;
	    oxx->cursorSeq = seq2;
	} else {
	    oxx->cursorPos = pos2;
	    oxx->cursorSeq = 0;
	}
    }

    /* Redraw the cursor, but only if still visible */
    if (onScreen(oxx, oxx->cursorSeq, oxx->cursorPos, NULL))
	showCursor(oxx, oxx->cursorSeq, oxx->cursorPos);

    /* Update trace positions */
    repositionTraces(oxx);

    /*
     * Only perform callback for this display, rather than all displays of
     * this data.
     */
    db_callback_tk(oxx, DBCALL_CURSOR_NOTIFY,
		   oxx->cursorSeq ? DB_Number(oxx, oxx->cursorSeq) : 0,
		   oxx->cursorPos, NULL);
}


/*
 * Sets the cursor position
 */
void setCursorPos(EdStruct *xx, int pos) {
    xx->refresh_flags |= ED_DISP_CURSOR;
    xx->cursorPos = pos;

    updateJoinCursor(xx, xx->cursorSeq, pos);

    /*
     * Only perform callback for this display, rather than all displays of
     * this data.
     */
    db_callback_tk(xx, DBCALL_CURSOR_NOTIFY,
		   xx->cursorSeq ? DB_Number(xx, xx->cursorSeq) : 0,
		   xx->cursorPos, NULL);
}


/*
 * Sets the cursor sequence
 */
void setCursorSeq(EdStruct *xx, int seq) {
    xx->refresh_flags |= ED_DISP_CURSOR;
    xx->cursorSeq = seq;

    /*
     * Only perform callback for this display, rather than all displays of
     * this data.
     */
    db_callback_tk(xx, DBCALL_CURSOR_NOTIFY,
		   xx->cursorSeq ? DB_Number(xx, xx->cursorSeq) : 0,
		   xx->cursorPos, NULL);
}


/*
 * Sets both cursor position and sequence (more efficient as it only
 * generates one CURSOR_NOTIFY request).
 */
void setCursorPosSeq(EdStruct *xx, int pos, int seq) {
    xx->refresh_flags |= ED_DISP_CURSOR;
    xx->cursorPos = pos;
    xx->cursorSeq = seq;

    updateJoinCursor(xx, xx->cursorSeq, pos);

    /*
     * Only perform callback for this display, rather than all displays of
     * this data.
     */
    db_callback_tk(xx, DBCALL_CURSOR_NOTIFY,
		   xx->cursorSeq ? DB_Number(xx, xx->cursorSeq) : 0,
		   xx->cursorPos, NULL);
}


/*
 *----------------------------------------------------------------------------
 * Displaying traces
 *---------------------------------------------------------------------------
 */

int get_trace_path(EdStruct *xx, int seq, char *fileName, char *t_type) {
    char t_fname[FILE_NAME_LENGTH+1]; /* file name of trace */
    int i;
    
    memset(t_fname, 0, FILE_NAME_LENGTH+1);

    if (readrd(handle_io(DBI_io(xx)), DB_Number(xx, seq),
	       t_type, t_fname,
	       sizeof(t_type), sizeof(t_fname)) != 0)
	return 1;
    
    t_type[4] = '\0';

    /* convert fortran string to c string */
    for (i=FILE_NAME_LENGTH-1;
	 i>=0 && (!t_fname[i] || isspace((int)t_fname[i]));
	 i--)
	;

    t_fname[++i] = '\0';
    
    /* skip if no raw data file for trace */
    if (t_fname[0] == '\0') return 1;
    
    /* don't check for existance - this is done later... */
    strcpy(fileName, t_fname);

    return 0;
}

DisplayContext *showTrace(EdStruct *xx, int seq, int pos, int baseSpacing,
			  int differencing, int mini_trace)
{
    int baseNum;
    int t_lcut;
    int t_rcut;
    char fileName[256];
    char t_type[5];
    int allow_dup = 0;
    DisplayContext *dc;

    if (!mini_trace) {
	int skip = 0;

	if (!differencing && (xx->diff_traces || xx->read_pair_traces)) {
	    if (0 == auto_diff(xx, seq, positionInContig(xx, seq, pos)))
		return NULL;
	    else
		skip = 1;
	} else if (differencing) {
	    allow_dup = 1;
	} else {
	    allow_dup = 0;
	}
	
	if (!skip && !differencing && xx->compare_trace >= 0) {
	    return diff_traces(xx, xx->compare_trace, seq,
			       positionInContig(xx, seq, pos));
	}
    }
    
    if (get_trace_path(xx, seq, fileName, t_type))
	return NULL;

    baseNum = origpos(xx, seq, DB_Start(xx,seq) + pos);
    t_lcut  = origpos(xx, seq, DB_Start(xx,seq));
    t_rcut  = origpos(xx, seq, DB_Start(xx,seq) + DB_Length(xx, seq));
    
    if (DB_Comp(xx,seq) == UNCOMPLEMENTED) {
	dc = tman_manage_trace(t_type, fileName, baseNum -1 , t_lcut,
			       t_rcut - t_lcut + 1, /*not complemented*/0,
			       baseSpacing, DBgetName(DBI(xx),seq),
			       xx, seq, allow_dup, mini_trace);
    } else {
	dc = tman_manage_trace(t_type, fileName, baseNum -1, t_lcut,
			       t_rcut - t_lcut + 1, /*complemented*/1,
			       baseSpacing, DBgetName(DBI(xx),seq),
			       xx, seq, allow_dup, mini_trace);
    }
    
    return dc;
}


void edInvokeTrace(EdStruct *xx) {
    int baseSpacing = xx->fontWidth * 2;
    int *slist;
    
    if (xx->cursorSeq) {
	showTrace(
		  xx,
		  xx->cursorSeq,
		  xx->cursorPos,
		  baseSpacing,
		  0,
		  0 /* full */);
    } else {
	int *seqList;
	int i, j, tmpcmp, tmpdiff, tmprp;

        seqList = sequencesInRegion(xx,xx->cursorPos,1);
	tmpcmp  = xx->compare_trace;
	tmpdiff = xx->diff_traces;
	tmprp   = xx->read_pair_traces;

	xx->compare_trace = -1;
	xx->diff_traces = 0;
	xx->read_pair_traces = 0;

	/* copy list */
	for (i=0; seqList[i]; i++);
	slist = (int *)xcalloc(i+1, sizeof(int));
	memcpy(slist, seqList, i * sizeof(int));

	/* Shut down existing traces */
	tman_shutdown_traces(xx, 2);

	for (i=0, j=MAXCONTEXTS; j && slist[i]; i++) {

	    /*
	     * Must only invoke traces covering the cursor position
	     */
	    if (xx->cursorPos - DB_RelPos(xx, slist[i]) +
		DB_Start(xx, slist[i]) < 0)
		continue;

	    showTrace(
		      xx,
		      slist[i],
		      xx->cursorPos-DB_RelPos(xx,slist[i])+1,
		      baseSpacing,
		      0,
		      0 /* full */);
	    j--;
	}
	
	xfree(slist);
	xx->compare_trace = tmpcmp;
	xx->diff_traces = tmpdiff;
	xx->read_pair_traces = tmprp;
    }
}


void repositionTraces(EdStruct *xx) {
    /* Full trace display only */
    tman_reposition_traces(xx,
			   xx->cursorPos + DB_RelPos(xx, xx->cursorSeq) - 1,
	                   0);
}



/*
 *----------------------------------------------------------------------------
 * Sequence manipulation - Low level routines
 *
 * These routines manipulate the memory directly and don't care about side
 * effects.
 *
 * All return:
 *    0 - Succes
 *    1 - Failure
 *----------------------------------------------------------------------------
 */


/*
 * Set cursor sequence and position for all displays of this info. Not too
 * useful at present!
 */
int _cursor_set(DBInfo *db, int seq, int pos) {
    DBI_callback(db, DBCALL_CURSOR, seq, pos, NULL);
    return 0;
}

/*
 * Adjust the relative position of a gel
 */
int _shift_right(DBInfo *db, int seq, int num_bases, int seq_flags)
{
    /* adjust relative position */
    _DBsetRelPos(db, seq, _DB_RelPos(db,seq) + num_bases);
    
    /* set flags */
    _DBsetFlags(db, seq, seq_flags);
    
    return 0;
}


/*
 * Adjust the relative position of a gel
 */
int _shift_left(DBInfo *db, int seq, int num_bases, int seq_flags)
{
    /* adjust relative position */
    _DBsetRelPos(db, seq, _DB_RelPos(db, seq) - num_bases);
    
    /* set flags */
    _DBsetFlags(db, seq, seq_flags);
    
    return 0;
}


/*
 * Insert bases into a sequence before base at position pos. We are flagged
 * whether or not we insert into the cutoff. This is due to a nasty special
 * case for undoing deletions of the base at position 0. We wish the insert in
 * this case to add to the cutoff data. Normal inserts to position zero add to
 * the used data.
 */
int _insert_bases(DBInfo *db, int seq, int pos, int num_bases, char *bases,
		  int1 *conf, int2 *opos, int seq_flags, int cutoff)
{
    int len = _DB_Length(db, seq), i;

    (void) DBgetSeq(db, seq); /* force sequence to be read */

    for (i = 0; i < num_bases; i++) {
	DBI_callback(db, DBCALL_INSERT, seq, pos, NULL);
    }

    {
	int tmp_len = _DB_Length2(db, seq);

	(void) io_insert_seq(&tmp_len,
			     &_DB_Start(db,seq),
			     &_DB_End(db,seq),
			     _DB_Seq(db,seq),
			     _DB_Conf(db,seq),
			     _DB_Opos(db,seq),
			     /**********/
			     pos+_DB_Start(db,seq),
			     bases,
			     conf,
			     opos,
			     num_bases
			     );
	_DB_Length2(db, seq) = tmp_len;
    }

    /* handle special case for a cutoff insert */
    if (cutoff && pos == 1) {
	_DBsetStart(db, seq, _DB_Start(db, seq) + num_bases);
    } else if (cutoff && pos == len + 1) {
	_DBsetEnd(db, seq, _DB_End(db, seq) - num_bases);
    } else {
	/* adjust length, only if we're inserting into the non cutoff data */
	if (pos <= _DB_Length(db, seq)+1 && pos > 0)
	    _DBsetLength(db,seq,len+num_bases);
    }

    /* set flags */
    _DBsetFlags(db,seq,seq_flags);

    return 0;
}


/*
 * Delete bases from a sequence before base at position pos
 */
int _delete_bases(DBInfo *db, int seq, int pos, int num_bases, int seq_flags)
{
    int len = _DB_Length(db,seq), i;

    (void) DBgetSeq(db,seq); /* force sequence to be read */
    
    for (i = 0; i < num_bases; i++) {
	DBI_callback(db, DBCALL_DELETE, seq, pos, NULL);
    }

    {
	int tmp_len = _DB_Length2(db, seq);

	(void) io_delete_seq(&tmp_len,
			     &_DB_Start(db,seq),
			     &_DB_End(db,seq),
			     _DB_Seq(db,seq),
			     _DB_Conf(db,seq),
			     _DB_Opos(db,seq),
			     /**********/
			     pos+_DB_Start(db,seq),
			     num_bases
			     );
	_DB_Length2(db,seq) = tmp_len;
    }

    /* adjust length, only if we're deleting from the non cutoff data */
    if (pos <= _DB_Length(db, seq) && pos > 0)
	_DBsetLength(db,seq,len-num_bases);
    
    /* set flags */
    _DBsetFlags(db,seq,seq_flags);

    return 0;
}


int _adjust_base_conf(DBInfo *db, int seq, int pos, int val, int opos,
		      int flags) {
    (_DB_Conf(db, seq) + _DB_Start(db, seq) - 1)[pos] = val;
    (_DB_Opos(db, seq) + _DB_Start(db, seq) - 1)[pos] = opos;

    /* set flags */
    _DBsetFlags(db, seq, flags);

    return 0;
}


/*
 * Replace bases of a sequence starting at base numbered pos
 */
int _replace_bases(DBInfo *db, int seq, int pos, int num_bases, char *bases,
		   int1 *conf, int2 *opos, int seq_flags,
		   int diff_only, int conf_only)
{
    (void) DBgetSeq(db, seq); /* force sequence to be read */

    {
	int tmp_len = _DB_Length2(db, seq);

	(void) io_replace_seq(&tmp_len,
			      &_DB_Start(db, seq),
			      &_DB_End(db, seq),
			      _DB_Seq(db, seq),
			      _DB_Conf(db, seq),
			      _DB_Opos(db, seq),
			      /**********/
			      pos + _DB_Start(db, seq),	
			      bases,
			      conf,
			      opos,
			      num_bases,
			      diff_only,
			      conf_only
			      );
	_DB_Length2(db, seq) = tmp_len;
    }

    /* set flags */
    _DBsetFlags(db, seq, seq_flags);
    
    return 0;
}


/*
 * Transpose two characters
 */
int _transpose_bases(DBInfo *db, int seq, int pos, int seq_flags) {
    char *seqp, tmp_seqp;
    int1 *conf, tmp_conf;
    int2 *opos, tmp_opos;

    (void)DBgetSeq(db, seq); /* force read */
    seqp = _DB_Seq(db, seq);
    conf = _DB_Conf(db, seq);
    opos = _DB_Opos(db, seq);

    if (pos < 0 || pos >= _DB_Length2(db, seq)) {
	bell();
	return -1;
    }

    pos += _DB_Start(db, seq);

    tmp_seqp = seqp[pos]; seqp[pos] = seqp[pos+1]; seqp[pos+1] = tmp_seqp;
    tmp_conf = conf[pos]; conf[pos] = conf[pos+1]; conf[pos+1] = tmp_conf;
    tmp_opos = opos[pos]; opos[pos] = opos[pos+1]; opos[pos+1] = tmp_opos;

    /* set flags */
    _DBsetFlags(db,seq,seq_flags);

    return 0;
}

/*
 * Set an individual sequence flag
 */
int _set_flags(DBInfo *db, int seq, int seq_flags) {
    /* set flags */
    _DBsetFlags(db,seq,seq_flags);

    return 0;
}

/*
 * Mark a sequence as the reference sequence
 */
int _set_reference_seq(DBInfo *db, int seq, int flags,
		       int refseq, int length, int offset){
    /* set flags */
    _DBsetFlags(db,seq,flags);

    /* set reference */
    db->reference_seq = refseq;
    db->reference_len = length;
    db->reference_offset = offset;

    return 0;
}

/*
 * Move the sequence seq into its proper position
 */
int _reorder_seq(DBInfo *db, int seq, int old_index, int new_index,
		 int seq_flags)
{
    int i;
    
    if (old_index < new_index) {
	/* moving right - so shuffle everything else left */
	for (i = old_index+1; i<= new_index; i++)
	    _DBI_order(db)[i-1] = _DBI_order(db)[i];
    } else {
	/* moving left - so shuffle everything else right */
	for (i = old_index-1; i>= new_index; i--)
	    _DBI_order(db)[i+1] = _DBI_order(db)[i];
    }
    _DBI_order(db)[new_index] = seq;
    
    
    /* set flags */
    _DBsetFlags(db,seq,seq_flags);
    
    return 0;
}



/*
 *----------------------------------------------------------------------------
 * Sequence manipulation - Mid level routines
 *
 * These are interfaces to the low level routines. They ensure the update
 * information is kept and allow for undo.
 *----------------------------------------------------------------------------
 */

int U_shift_right(DBInfo *db, int seq, int num_bases)
{
    
    int flags, old_flags;
    UndoStruct *u;
    
    old_flags = _DB_Flags(db,seq);
    
    /* UNDO - _shift_left(db,seq,num_bases,old_flags) */
    if ( (u = newUndoStruct(db)) != NULL ) {
	u->db = db;
	u->command = UndoShiftLeft;
	u->sequence = seq;
	u->info.shift_left.num_bases = num_bases;
	u->info.shift_left.flags = old_flags;
	recordUndo(db,u);
    }
    
    flags = old_flags | DB_FLAG_REL_MODIFIED;
    return _shift_right(db,seq,num_bases,flags);
}



int U_shift_left(DBInfo *db, int seq, int num_bases)
{
    int flags, old_flags;
    UndoStruct *u;
    
    old_flags = _DB_Flags(db,seq);
    /* UNDO - _shift_right(db),seq,num_bases,old_flags) */
    if ( (u = newUndoStruct(db)) != NULL ) {
	u->db = db;
	u->command = UndoShiftRight;
	u->sequence = seq;
	u->info.shift_right.num_bases = num_bases;
	u->info.shift_right.flags = old_flags;
	recordUndo(db,u);
    }
    
    flags = old_flags | DB_FLAG_REL_MODIFIED;
    return _shift_left(db,seq,num_bases,flags);
}



int U_insert_bases(EdStruct *xx, int seq, int pos, int num_bases, char *bases)
{
    int flags, old_flags, ret;
    UndoStruct *u;
    int1 new_conf[100], *new_confp, *free_me;
    int1 *old_conf;
    int i;
    
    if (num_bases <= 100) {
	new_confp = new_conf;
	free_me = NULL;
    } else {
	free_me = new_confp = (int1 *)xmalloc(num_bases * sizeof(int1));
    }

    old_flags = DB_Flags(xx,seq);
    (void) DBgetSeq(DBI(xx),seq);
    old_conf = DB_Conf(xx,seq) + DB_Start(xx,seq);

    /* UNDO - _delete_bases(DBI(xx),seq,pos,num_bases,old_flags) */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoDeleteBases;
	u->sequence = seq;
	u->info.delete_bases.position = pos > 0 ? pos : pos-num_bases;
	u->info.delete_bases.num_bases = num_bases;
	u->info.delete_bases.flags = old_flags;
	recordUndo(DBI(xx),u);
    }
    
    flags = old_flags | DB_FLAG_REL_MODIFIED | DB_FLAG_SEQ_MODIFIED;
    
    for (i = 0; i < num_bases; i++)
	new_confp[i] = bases[i] == '-' ? 0 : xx->default_conf_n;

    ret = _insert_bases(DBI(xx), seq, pos, num_bases, bases,
			new_confp, NULL, flags, 0);

    RedisplaySeq(xx, seq);
    if (pos < 1)
	U_adjust_cursor(xx, -num_bases);

    if (free_me) xfree(free_me);
    return ret;
}


/*
 * YUK!
 * Every time this routine is called from this file, tagDeleteBases()
 * also needs to be called.
 * When this is called from elsewhere, things are ok
 */
int U_delete_bases(EdStruct *xx, int seq, int pos, int num_bases)
{
    int flags, old_flags, ret;
    char *old_bases;
    int1 *old_conf;
    int2 *old_opos;
    UndoStruct *u;
    
    old_flags = DB_Flags(xx,seq);
    (void) DBgetSeq(DBI(xx),seq);
    old_bases = DB_Seq(xx,seq) + DB_Start(xx,seq);
    old_conf = DB_Conf(xx,seq) + DB_Start(xx,seq);
    old_opos = DB_Opos(xx,seq) + DB_Start(xx,seq);
    /* UNDO -
       _insert_bases(DBI(xx),seq,pos,num_bases,bases,old_flags,1)
       
       where:
       s     <- DBgetSeq(DBI(xx),seq);
       bases <- strncpy(bases,&seq[(pos-num_bases+1)-1],num_bases);
       */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoInsertBases;
	u->sequence = seq;
	u->info.insert_bases.position = pos > 0 ? pos : pos + num_bases;
	if (pos == 0 ||
	    (pos == DB_End(xx,seq) - DB_Start(xx,seq)))
	    u->info.insert_bases.cutoff = 1;
	else
	    u->info.insert_bases.cutoff = 0;
	u->info.insert_bases.num_bases = num_bases;
	u->info.insert_bases.flags = old_flags;
	/* copy sequence */
	packBCO(&u->info.insert_bases.b_c_o,
		&old_bases[pos-1],
		&old_conf[pos-1],
		&old_opos[pos-1],
		num_bases);
	recordUndo(DBI(xx),u);
    }
    
    flags = old_flags | DB_FLAG_REL_MODIFIED | DB_FLAG_SEQ_MODIFIED;
    ret = _delete_bases(DBI(xx),seq,pos,num_bases,flags);
    RedisplaySeq(xx, seq);
    if (pos < 1)
	U_adjust_cursor(xx, num_bases);

    return ret;
}


int U_replace_bases(EdStruct *xx, int seq, int pos, int num_bases, char *bases,
		    int diff_only)
{
    int flags, old_flags;
    char *old_bases;
    int1 *old_conf;
    int2 *old_opos;
    UndoStruct *u;
    int1 new_conf[100], *new_confp, *free_me;
    int rval;

    if (num_bases <= 100) {
	new_confp = new_conf;
	free_me = NULL;
    } else {
	free_me = new_confp = (int1 *)xmalloc(num_bases * sizeof(int1));
    }

    
    if (!xx->reveal_cutoffs && pos > DB_Length(xx, seq))
	return 0;

    old_flags = DB_Flags(xx,seq);
    (void) DBgetSeq(DBI(xx),seq);
    old_bases = DB_Seq(xx,seq) + DB_Start(xx,seq);
    old_conf = DB_Conf(xx,seq) + DB_Start(xx,seq);
    old_opos = DB_Opos(xx,seq) + DB_Start(xx,seq);
    /* UNDO -
       _replace_bases(DBI(xx),seq,pos,num_bases,old_bases,old_flags)
       
       where:
       s         <- DBgetSeq(DBI(xx),seq);
       old_bases <- strncpy(bases,&seq[(pos-num_bases+1)-1],num_bases);
       */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoReplaceBases;
	u->sequence = seq;
	u->info.replace_bases.position = pos;
	u->info.replace_bases.num_bases = num_bases;
	u->info.replace_bases.flags = old_flags;
	if (diff_only)
	    u->info.replace_bases.flags |= DB_FLAG_TMP_DIFF_ONLY;
	/* copy sequence */
	packBCO(&u->info.insert_bases.b_c_o,
		&old_bases[(pos)-1],
		&old_conf[(pos)-1],
		&old_opos[(pos)-1],
		num_bases);
	recordUndo(DBI(xx),u);
    }
    
    flags = old_flags | DB_FLAG_REL_MODIFIED | DB_FLAG_SEQ_MODIFIED;

    if (xx->default_conf_r != -1) {
	int i;
	for (i = 0; i < num_bases; i++)
	    new_confp[i] = bases[i] == '-' ? 0 : xx->default_conf_r;
    } else {
	new_confp = &old_conf[pos-1];
    }

    if (0 != _replace_bases(DBI(xx), seq, pos, num_bases, bases,
			    new_confp, NULL, flags, diff_only, 0))
	rval = 0;
    else
	rval = num_bases;

    if (free_me) xfree(free_me);
    return rval;
}


/*
 * Changes the confidence values of bases that disagree with the sequence.
 */
int U_replace_conf(EdStruct *xx, int seq, int pos, int num_bases, char *bases)
{
    int flags, old_flags;
    char *old_bases;
    int1 *old_conf;
    int2 *old_opos;
    UndoStruct *u;
    int1 *conf;
    int r;
    
    if (!xx->reveal_cutoffs && pos > DB_Length(xx, seq))
	return 0;

    if (NULL == (conf = (int1 *)xcalloc(num_bases, sizeof(*conf))))
	return 0;

    old_flags = DB_Flags(xx,seq);
    (void) DBgetSeq(DBI(xx),seq);
    old_bases = DB_Seq(xx,seq) + DB_Start(xx,seq);
    old_conf = DB_Conf(xx,seq) + DB_Start(xx,seq);
    old_opos = DB_Opos(xx,seq) + DB_Start(xx,seq);
    /* UNDO -
       _replace_bases(DBI(xx),seq,pos,num_bases,old_bases,old_conf,NULL,
                      old_flags)
       
       where:
       s         <- DBgetSeq(DBI(xx),seq);
       old_bases <- strncpy(bases,&seq[(pos-num_bases+1)-1],num_bases);
       */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoReplaceBases;
	u->sequence = seq;
	u->info.replace_bases.position = pos;
	u->info.replace_bases.num_bases = num_bases;
	u->info.replace_bases.flags = old_flags
	    | DB_FLAG_TMP_CONF_ONLY | DB_FLAG_TMP_DIFF_ONLY;
	/* copy sequence */
	packBCO(&u->info.insert_bases.b_c_o,
		bases,
		&old_conf[(pos)-1],
		&old_opos[(pos)-1],
		num_bases);
	recordUndo(DBI(xx),u);
    }
    
    flags = old_flags | DB_FLAG_REL_MODIFIED | DB_FLAG_SEQ_MODIFIED;

    /*
     * YUK! adjust tags - Do this here?
     */
    if (0 != _replace_bases(DBI(xx),seq,pos,num_bases,bases,conf,NULL,flags,
			    1 /*diff_only*/, 1 /*conf_only*/))
	r = 0;
    else
	r = num_bases;

    xfree(conf);
    return r;
}



int U_reorder_seq(EdStruct *xx, int seq, int old_index, int new_index)
{
    int flags, old_flags;
    UndoStruct *u;
    
    old_flags = DB_Flags(xx,seq);
    /* UNDO - _reorder_seq(DBI(xx), seq, new_index, old_index, old_flags) */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoReorderSeq;
	u->sequence = seq;
	u->info.reorder_seq.old_id = new_index;
	u->info.reorder_seq.new_id = old_index;
	u->info.reorder_seq.flags = old_flags;
	recordUndo(DBI(xx),u);
    }
    
    flags = old_flags | DB_FLAG_REL_MODIFIED;
    xx->refresh_flags |= ED_DISP_SEQS | ED_DISP_NAMES;
    return _reorder_seq(DBI(xx), seq, old_index, new_index, flags);
}




void U_change_consensus_length(EdStruct *xx, int len)
{
    int old_len;
    UndoStruct *u;
    
    old_len = DB_Length(xx,0);
    /* UNDO - DBsetLength(xx,0, old_len) */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoConsensusLength;
	u->info.consensus_length.num_bases = old_len;
	recordUndo(DBI(xx),u);
    }
    
    DBsetLength(xx,0,len);
    DBsetLength2(xx,0,len);
    xx->refresh_flags |= ED_DISP_CONS | ED_DISP_RULER | ED_DISP_STATUS | \
	ED_DISP_SCROLL;
}



void U_adjust_cursor(EdStruct *xx, int n)
{
    UndoStruct *u;
    int old_pos;
    
    old_pos = xx->cursorPos;
    /* UNDO - xx->cursorPos = old_pos; */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoAdjustCursor;
	u->sequence = xx->cursorSeq;
	u->info.adjust_cursor.xx = xx;
	u->info.adjust_cursor.position = old_pos;
	u->info.adjust_cursor.editor_id = xx->editor_id;
	recordUndo(DBI(xx),u);
    }
    
    setCursorPos(xx, xx->cursorPos + n);
}


void U_adjust_base_conf(EdStruct *xx, int seq, int pos, int val) {
    int flags, old_flags;
    int1 *old_conf;
    int2 *old_opos;
    UndoStruct *u;
    int2 new_opos;

    old_flags = DB_Flags(xx,seq);
    old_conf = DB_Conf(xx,seq) + DB_Start(xx,seq);
    old_opos = DB_Opos(xx,seq) + DB_Start(xx,seq);
    
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoAdjustBaseConf;
	u->sequence = seq;
	u->info.adjust_base_conf.position = pos;
	u->info.adjust_base_conf.conf = old_conf[pos-1];
	u->info.adjust_base_conf.flags = old_flags;
	u->info.adjust_base_conf.opos = old_opos[pos-1];
	recordUndo(DBI(xx), u);
    }

    flags = old_flags | DB_FLAG_SEQ_MODIFIED | DB_FLAG_REL_MODIFIED;

    new_opos = 0;
    _adjust_base_conf(DBI(xx), seq, pos, val, new_opos, flags);
}

void U_transpose_bases(EdStruct *xx, int seq, int pos) {
    int flags, old_flags;
    UndoStruct *u;

    old_flags = DB_Flags(xx, seq);

    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoTransposeBases;
	u->sequence = seq;
	u->info.transpose_bases.position = pos;
	u->info.transpose_bases.flags = old_flags;
	recordUndo(DBI(xx), u);
    }

    flags = old_flags | DB_FLAG_SEQ_MODIFIED | DB_FLAG_REL_MODIFIED;

    _transpose_bases(DBI(xx), seq, pos, flags);
}

void U_set_flags(EdStruct *xx, int seq, int flags) {
    int old_flags;
    UndoStruct *u;

    old_flags = DB_Flags(xx, seq);

    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoSetFlags;
	u->sequence = seq;
	u->info.set_flags.flags = old_flags;
	recordUndo(DBI(xx), u);
    }

    if ( (flags & DB_FLAG_REFTRACE) != (old_flags & DB_FLAG_REFTRACE) )
	flags |= DB_FLAG_NOTE_MODIFIED;

    _set_flags(DBI(xx), seq, flags);
}

void U_set_reference_seq(EdStruct *xx, int seq, int refseq, int length,
			 int offset) {
    int flags, old_flags;
    UndoStruct *u;

    old_flags = DB_Flags(xx, seq);

    /* Remove old reference seq if appropriate */
    if (refseq && DBI(xx)->reference_seq) {
	int old_ref = DBI(xx)->reference_seq;
	DBI(xx)->reference_seq = 0;
	U_set_reference_seq(xx, old_ref, 0, 0, 0);
    }

    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoSetReferenceSeq;
	u->sequence = seq;
	u->info.set_reference_seq.flags = old_flags;
	u->info.set_reference_seq.refseq = DBI(xx)->reference_seq;
	u->info.set_reference_seq.length = DBI(xx)->reference_len;
	u->info.set_reference_seq.offset = DBI(xx)->reference_offset;
	recordUndo(DBI(xx), u);
    }

    old_flags |= DB_FLAG_NOTE_MODIFIED;
    if (refseq) {
	flags = old_flags | DB_FLAG_REFSEQ;
    } else { 
	flags = old_flags &~ DB_FLAG_REFSEQ;
    }
    _set_reference_seq(DBI(xx), seq, flags, refseq, length, offset);
}


/*
 *----------------------------------------------------------------------------
 * Sequence manipulation - High level routines
 *
 * These are interfaces to the mid level routines. Each routine performs a
 * logical function
 *----------------------------------------------------------------------------
 */

/***********************************
 ** Routines for simplifying code **
 ***********************************/

/*
 * Handle U_delete_bases() and tagDeleteBases()
 */
int handle_delete_bases(EdStruct *xx, int seq, int pos, int num_bases)
{
    char *bases;
    /* delete tag bases first */
    bases = DBgetSeq(DBI(xx),seq);
    tagDeleteBases(xx,seq,pos,num_bases);
    /* delete bases */
    return U_delete_bases(xx,seq,pos,num_bases);
}




/*
 * Handle both U_insert_bases() and tagInsertBases()
 */
int handle_insert_bases(EdStruct *xx, int seq, int pos, int num_bases,
			char *bases)
{
    /* insert tag bases first */
    tagInsertBases(xx,seq,pos,num_bases);
    /* insert bases */
    return U_insert_bases(xx,seq,pos,num_bases,bases);
}



/********************
 ** The real stuff **
 ********************/


/*
 * Adjust the relative position of a gel
 */
int shiftRight(EdStruct *xx, int seq, int num_bases)
{
    int i;
    int seq_index;
    int new_seq_index;
    int left_enders;
    int new_left_index;
    int pos;
    int end;
    
    if (seq==0) {
	/* we can't shift the consensus right! */
	return 1;
    }
    
    seq_index = 0;	/* init. to keep compiler happy - not really needed */
    new_left_index = 0;	/* init. to keep compiler happy - not really needed */
    left_enders = 0;

    if (DB_RelPos(xx,seq)==1) {
	/*
	 * we must consider a special condition here:
	 * when this is the left-most gel, shifting this gel right
	 * has the effect of shifting everything else left!
	 */
	new_left_index = 0;
	for (i=1;
	     i<=DBI_gelCount(xx) && DB_RelPos(xx,DBI_order(xx)[i]) <= num_bases;
	     i++) {
	    if (DB_RelPos(xx,DBI_order(xx)[i]) == 1)
		left_enders++;

	    /* keep a record of index of first reading not seq */
	    if (DBI_order(xx)[i] != seq && ! new_left_index)
		new_left_index = i;

	    /* keep a record of seq's index */
	    if (DBI_order(xx)[i] == seq)
		seq_index = i;
	}

	new_seq_index = i-1;
    } else {
	/* determine shuffling */
	pos = DB_RelPos(xx,seq);
	seq_index = seqToIndex(xx,seq);
	for (i = seq_index;
	     i <= DBI_gelCount(xx)
	         && DB_RelPos(xx,DBI_order(xx)[i]) < pos+num_bases;
	     i++);
	new_seq_index = i-1;
    }

    if (left_enders == 1) {
	int all_left, this_right;

	/* special condition met! */
	if (new_left_index) {
	    all_left = DB_RelPos(xx,DBI_order(xx)[new_left_index]) - 1;
	    this_right = num_bases - all_left;
	} else {
	    all_left = num_bases;
	    this_right = 0;
	}
	
	if (this_right)
	    U_shift_right(DBI(xx),seq,this_right);
	
	if (all_left) {
	    /* move all other sequences all_right */
	    for (i=1; i<seq;i++)
		U_shift_left(DBI(xx),i,all_left);
	    for (i=seq+1; i<=DBI_gelCount(xx);i++)
		U_shift_left(DBI(xx),i,all_left);
	}
	
    } else {
	/* only need to adjust this sequence */
	U_shift_right(DBI(xx),seq,num_bases);
    }
    RedisplaySeq(xx, seq);
    
    if (new_seq_index != seq_index) {
	/* shuffling required */
	U_reorder_seq(xx, seq, seq_index, new_seq_index);
    }
    
    /* adjust consensus length */
    if (DB_RelPos(xx, seq) <= num_bases + 1 ||
	DB_RelPos(xx, seq) + DB_Length(xx, seq) + num_bases + 1
	>= DB_Length(xx, 0)) {
	end = calculate_consensus_length(xx); 

	if (DB_Length(xx,0) != end) {
	    /* xx->displayPos += num_bases; */
	    U_change_consensus_length(xx,end);
	    U_adjust_cursor(xx, 0);
	}
    }
    invalidate_consensus(xx);
    
    return 0;
}



/*
 * Adjust the relative position of a gel
 */
int shiftLeft(EdStruct *xx, int seq, int num_bases)
{
    int i;
    int seq_index;
    int new_seq_index;
    int pos;
    int end;
    
    if (seq==0) {
	/* we can't shift the consensus right! */
	return 1;
    }
    
    /* determine shuffling */
    pos = DB_RelPos(xx,seq);
    seq_index = seqToIndex(xx,seq);

    for (i=seq_index; i>0 && DB_RelPos(xx,DBI_order(xx)[i])>pos-num_bases;i--)
	;
    new_seq_index = i+1;
    
    
    if (DB_RelPos(xx,seq) <= num_bases) {
	int all_right, this_left;
	
	this_left = DB_RelPos(xx,seq) - 1;
	all_right = num_bases - this_left;
	
	if (this_left)
	    U_shift_left(DBI(xx),seq,this_left);
	
	if (all_right) {
	    /* move all other sequences all_right */
	    for (i=1; i<seq;i++)
		U_shift_right(DBI(xx),i,all_right);

	    for (i=seq+1; i<=DBI_gelCount(xx);i++)
		U_shift_right(DBI(xx),i,all_right);
	}
	
    } else {
	/* only need to adjust this sequence */
	U_shift_left(DBI(xx),seq,num_bases);
    }
    
    if (new_seq_index != seq_index) {
	/* shuffling required */
	U_reorder_seq(xx, seq, seq_index, new_seq_index);
    }
    RedisplaySeq(xx, seq);
    
    /* adjust consensus length */
    if (DB_RelPos(xx, seq) <= num_bases + 1 ||
	DB_RelPos(xx, seq) + DB_Length(xx, seq) + num_bases + 1
	>= DB_Length(xx, 0)) {
	end = calculate_consensus_length(xx);
	if (DB_Length(xx,0) != end) {
	    /* xx->displayPos -= num_bases; */
	    U_change_consensus_length(xx,end);
	    U_adjust_cursor(xx, 0);
	}
    }
    invalidate_consensus(xx);
    
    return 0;
}



/*
 * Insert the bases into the sequence, before position pos
 *
 * Returns: number of bases inserted
 */
int insertBases(EdStruct *xx, int seq, int pos, int num_bases, char *bases)
{
    int len;
    int max_len;
    int overflow;
    int end;
    
    /* should use insertBasesConsensus instead */
    if (seq==0) return 0;
    
    (void) DBgetSeq(DBI(xx),seq);

    /* check if we'll need to extend our allocated memory for the seqs */
    len = DB_Length2(xx,seq);
    max_len = DB_Alloced(xx, seq);
    overflow = len + num_bases - max_len;

    if (overflow > 0) {
	size_t length = len + num_bases + SEQ_LENGTH_INC + len * 0.1;

	/*
	 * Realloc more space.
	 * What do we do if there isn't any more!? Worry about it later...
	 */
	DBsetSeq (xx, seq, (char *)xrealloc(DB_Seq (xx, seq),
					    length));
	DBsetOpos(xx, seq, (int2 *)xrealloc(DB_Opos(xx, seq),
					    length * sizeof(int2)));
	DBsetConf(xx, seq, (int1 *)xrealloc(DB_Conf(xx, seq),
					    length * sizeof(int1)));
	DBsetAlloced(xx, seq, length);
    }

    if (num_bases) {
	/* handle insert */
	handle_insert_bases(xx,seq,pos,num_bases,bases);
	
	/* may need to adjust consensus length */
	end = DB_RelPos(xx,seq) + DB_Length(xx,seq) - 1;

	if (DB_Length(xx,0) < end) {
	    U_change_consensus_length(xx,end);
	}
    }

    invalidate_consensus(xx);

    return num_bases;
}




/*
 * Insert the bases into the sequence, before position pos
 *
 * Returns: number of bases deleted
 */
int deleteBases(EdStruct *xx, int seq, int pos, int num_bases)
{
    int len;
    int end;
    int old_end;
    
    /* should use deleteBasesConsensus instead */
    if (seq==0) return 0;
    
    /* check if we are going to run out of gel to delete */
    len = DB_Length(xx,seq);
    if (num_bases > len) num_bases = len;
    
    if (num_bases) {
	handle_delete_bases(xx,seq,pos,num_bases);
	
	/* may need to adjust consensus length */
	old_end = DB_RelPos(xx,seq) + len - 1;
	if (old_end == DB_Length(xx,0)) {
	    end = calculate_consensus_length(xx);
	    U_change_consensus_length(xx,end);
	}
    }

    invalidate_consensus(xx);
    
    return num_bases;
}


/*
 * Replace the bases into the sequence, before position pos
 *
 * Returns: the number of bases replaced
 */
int replaceBases(EdStruct *xx, int seq, int pos, int num_bases, char *bases)
{
    int len;
    
    /* should use replaceBasesConsensus instead */
    if (seq==0) return 0;
    
    /* check if we are going to overrun maximum gel length */
    len = DB_Length2(xx,seq) - DB_Start(xx, seq);
    
    /* if we are we should trim num_bases accordingly */
    if (num_bases > (len-pos+1))
	num_bases = (len-pos+1);
    
    if (num_bases > 0) {
	num_bases = U_replace_bases(xx,seq,pos,num_bases,bases, 0);
	RedisplaySeq(xx, seq);
#ifdef CACHE_CONSENSUS_2
	/* Update consensus buffers */
	invalidate_consensus(xx);
#endif
    } else {
	bell();
    }
    
    /* no need to adjust consensus length */
    
    return num_bases;
}


/*
 * Returns:
 *	0 - OK
 *	1 - Error
 */
int replaceBasesConsensus(EdStruct *xx, int pos, int num_bases, char *bases)
{
    int ind;
    int i;
    /* following are for gels being modified */
    int gel_rel_pos;
    int gel_pos;
    int gel_len;
    int short_fall, over_shoot;
    int this_num_bases;
    char *this_bases;
    
    /* find the first sequence that overlaps pos */
    /*
    int max_len;
    int try_pos;
    max_len = dbi_max_gel_len(DBI(xx),0);
    try_pos = pos - max_len;
    ind = posToIndex(xx,try_pos);
    */
    ind = 1;
    
    /* didn't find index */
    if (!ind) return 1;
    
    /*
     * replace in all sequences that overlap
     */
    for (i=ind;
	 i <= DBI_gelCount(xx) && 
	 DB_RelPos(xx,DBI_order(xx)[i]) <= (pos + num_bases - 1) ;
	 i++) {
	gel_rel_pos = DB_RelPos(xx,DBI_order(xx)[i]);
	gel_pos = pos - gel_rel_pos + 1;
	gel_len =  DB_Length(xx,DBI_order(xx)[i]);
	
	this_num_bases = num_bases;
	this_bases = bases;
	
	/* check for falling short of start */
	short_fall = gel_rel_pos - pos;
	if (short_fall > 0) {
	    if (short_fall >= num_bases) {
		/* no overlap - skip */
		continue;
	    } else {
		gel_pos = 1;
		this_bases = bases+short_fall;
		this_num_bases = num_bases-short_fall;
	    }
	}
	
	/* check for overshooting the end */
	over_shoot = (this_num_bases + pos - 1) -
	    (gel_rel_pos + gel_len - 1);
	
	if (over_shoot > 0) {
	    if (over_shoot >= num_bases) {
		/* no overlap - skip */
		continue;
	    } else {
		this_num_bases = num_bases - over_shoot;
	    }
	}
	
	if (this_num_bases) {
	    char *old_bases;
	    int j;
	    
	    /*
	     * To optimise this we only make changes, and only store
	     * undo info, for places where we have at least one edit
	     * to do for this reading.
	     */
	    (void) DBgetSeq(DBI(xx),DBI_order(xx)[i]);
	    old_bases = DB_Seq(xx,DBI_order(xx)[i]) +
		DB_Start(xx,DBI_order(xx)[i]);
	    
	    for (j = 0; j < this_num_bases; j++) {
		if (toupper(this_bases[j]) != 
		    toupper(old_bases[j+gel_pos-1]))
		    break;
	    }
	    
	    if (j != this_num_bases) {
		if (xx->super_edit & SUPEREDIT_MODIFY_CONF)
		    U_replace_conf(xx, DBI_order(xx)[i], gel_pos,
				   this_num_bases, this_bases);
		else
		    U_replace_bases(xx, DBI_order(xx)[i], gel_pos,
				    this_num_bases, this_bases, 1);
	    }
	}
    }

#ifdef CACHE_CONSENSUS_2
    /* Update consensus buffers */
    invalidate_consensus(xx);
#endif

    xx->refresh_flags |= ED_DISP_SEQS | ED_DISP_STATUS;
    
    /* No need to adjust consensus length */
    return 0;
}


/*
 * Returns:
 *	0 - OK
 *	1 - Error
 */
int insertBasesConsensus(EdStruct *xx, int pos, int num_bases, char *bases)
{
    int ind;
    int end;
    int i;
    /* following are for gels being modified */
    int gel_rel_pos;
    int gel_pos;
    int gel_len;
    
    /* find the first sequence that overlaps pos */
    /*
    int max_len;
    int try_pos;
    max_len = dbi_max_gel_len(DBI(xx),0);
    try_pos = pos - max_len;
    ind = posToIndex(xx,try_pos);
    */
    ind = 1;
    
    /* didn't find index */
    if (!ind) return 1;
    
    /*
     * insert into all sequences
     */
    for (i=ind; i <= DBI_gelCount(xx); i++)
	{
	    gel_rel_pos = DB_RelPos(xx,DBI_order(xx)[i]);
	    
	    if (gel_rel_pos > pos) {
		/*
		 * Readings falling to the right are shifted
		 */
		/*
		 * fixme!
		 * To optimise this, have one routine to shift all readings
		 * right!
		 */
		U_shift_right(DBI(xx),DBI_order(xx)[i],num_bases);
	    } else {
		/*
		 * if hit any part of sequence then we should insert
		 */
		gel_pos = pos - gel_rel_pos + 1;
		gel_len =  DB_Length(xx,DBI_order(xx)[i]);
		if ( (gel_rel_pos + gel_len - 1) + 1 >= pos )
		    /* should determine how many bases to insert etc */
		    insertBases(xx,DBI_order(xx)[i],gel_pos, num_bases, bases);
	    }
	}
    
    /*
     * Notify for the consensus insert. This makes sure that the selection
     * is kept correct if it's on the consensus
     */
    for (i = 0; i < num_bases; i++)
	DBI_callback(DBI(xx), DBCALL_INSERT, 0, pos, NULL);

    /* adjust consensus length */
    end = calculate_consensus_length(xx);
    if (DB_Length(xx,0) != end) {
	U_change_consensus_length(xx,end);
    }

    /* adjust annotations */
    tagInsertBases(xx, 0, pos, num_bases);

    xx->refresh_flags |= ED_DISP_SEQS | ED_DISP_SCROLL;

    invalidate_consensus(xx);
    
    return 0;
}


/*
 * Returns:
 *	0 - OK
 *	1 - Error
 */
int deleteBasesConsensus(EdStruct *xx, int pos, int num_bases)
{
    int ind;
    int i;
    int end;
    /* following are for gels being modified */
    int gel_rel_pos;
    int gel_pos;
    int gel_len;
    int this_num_bases, this_shift_left;
    
    /* find the first sequence that overlaps pos */
    if (pos <= num_bases) {
	/* deleting close to start */
	num_bases = pos;
	ind = 1;
    } else {
	/*
	int max_len;
	int try_pos;
	max_len = dbi_max_gel_len(DBI(xx),0);
	try_pos = pos - max_len - num_bases;
	ind = posToIndex(xx,try_pos);
	*/
	ind = 1;
	/* didn't find index */
	if (!ind) return 1;
    }
    
    /*
     * insert into all sequences
     */
    for (i=ind; i <= DBI_gelCount(xx); i++)
	{
	    gel_rel_pos = DB_RelPos(xx,DBI_order(xx)[i]);
	    gel_pos = pos - gel_rel_pos + 1;
	    gel_len =  DB_Length(xx,DBI_order(xx)[i]);
	    
	    this_num_bases = num_bases;
	    this_shift_left = 0;
	    
	    /* four special cases to consider */
	    if (gel_pos >= gel_len + num_bases) {
		/* case 1 - this gel is not affected */
		this_num_bases = 0;
	    }

	    if (gel_pos > gel_len) {
		/* case 2 - some trimming of right end of reading needed */
		this_num_bases = gel_len + 1 - (gel_pos - this_num_bases);
		gel_pos = gel_len + 1;
	    }

	    if (gel_pos > 0 && gel_pos < this_num_bases) {
		/* case 3 - will consume left end - and more! */
		this_num_bases = gel_pos - 1;
		this_shift_left = num_bases - this_num_bases;
	    }

	    if (gel_pos < 1) {
		/* case 4 - this gel needs shifting left */
		/* fixme
		 * Optimise by shift all readings left in one call
		 */
		this_shift_left = num_bases;
		this_num_bases = 0;
	    }
	    
	    if (this_num_bases > 0) {
		/* delete if necessary */
		handle_delete_bases(xx,DBI_order(xx)[i],gel_pos,this_num_bases);
	    }
	    
	    if (this_shift_left > 0) {
		/* shift if necessary */
		U_shift_left(DBI(xx),DBI_order(xx)[i],this_shift_left);
	    }
	    
	}
    
    /*
     * Notify for the consensus delete. This makes sure that the selection
     * is kept correct if it's on the consensus
     */
    for (i = 0; i < num_bases; i++)
	DBI_callback(DBI(xx), DBCALL_DELETE, 0, pos, NULL);

    invalidate_consensus(xx);
    
    /* adjust consensus length */
    end = calculate_consensus_length(xx);
    if (DB_Length(xx,0) != end) {
	U_change_consensus_length(xx,end);
    }
    
    /* adjust annotations */
    tagDeleteBases(xx, 0, pos, num_bases);

    xx->refresh_flags |= ED_DISP_SEQS | ED_DISP_SCROLL;
    
    return 0;
}



/*
 * Returns:
 * 0 = Success
 * 1 = Error
 */
int adjustBaseConf(EdStruct *xx, int seq, int pos, int val, int move) {
    if (seq == 0)
	return 1;

    openUndo(DBI(xx));

    U_adjust_base_conf(xx, seq, pos, val);
    if (move)
	U_adjust_cursor(xx,1);

    closeUndo(xx, DBI(xx));
    RedisplaySeq(xx, seq);

#ifdef CACHE_CONSENSUS_2
    /* Update consensus buffers */
    invalidate_consensus(xx);
#endif

    redisplayWithCursor(xx);

    return 0;
}


/*
 * Switches two bases around maintaining original position information
 */
int transpose(EdStruct *xx, int seq, int pos, int dir, int num_bases) {
    int i;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }

    /* does not work for the consensus */
    if (seq == 0) {
	return 1;
    }
    
    /* can only shuffle pads when not in super edit */
    if (DBgetSeq(DBI(xx), seq)[pos-1] != '*' &&
	!(xx->super_edit & SUPEREDIT_TRANSPOSE_ANY))
	return 1;

    /* doit */
    openUndo(DBI(xx));

    for (i=0; i<num_bases; i++) {
	U_transpose_bases(xx, seq, pos-1 - (dir==-1));
	U_adjust_cursor(xx, dir);
    }

    closeUndo(xx, DBI(xx));

    invalidate_consensus(xx);

    redisplayWithCursor(xx);

    return 0;
}

/*
 * Sets a reference trace.
 */
int set_reference_trace(EdStruct *xx, int seq, int flag) {
    int flags;

    if (seq == 0)
	return 1;

    flags = DB_Flags(xx, seq);
    flags &= ~DB_FLAG_REFTRACE;
    flags |= flag;

    openUndo(DBI(xx));
    U_set_flags(xx, seq, flags);
    closeUndo(xx, DBI(xx));

    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 1);

    return 0;
}

/*
 * Sets a reference sequence.
 */
int set_reference_seq(EdStruct *xx, int seq, int refseq,
		      int length, int offset) {
    openUndo(DBI(xx));
    U_set_reference_seq(xx, seq, refseq, length, offset);
    closeUndo(xx, DBI(xx));

    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 1);

    return 0;
}


/*
 *----------------------------------------------------------------------------
 * Miscellaneous utility routines
 *---------------------------------------------------------------------------
 */

/**
 * Counts the number of good and bad spanning templates based on their
 * lengths. Inconsistencies due to orientation are "sort of" worked out,
 * but only if they are inconsistent internally.
 */
static void spanning_template_stats(EdStruct *xx[2], int overlapLength,
				    int *ptgood, int *ptbad) {
    int i;
    int tgood = 0, tbad = 0;
    int c1, c2;
    int numt = Ntemplates(DBI_io(*xx));
    int offset = editorLockedPos(xx, 1/*force recalculation*/);

    c1 = offset < 0 ? DBI_contigNum(xx[0]) : DBI_contigNum(xx[1]);
    c2 = offset < 0 ? DBI_contigNum(xx[1]) : DBI_contigNum(xx[0]);

    for (i = 1; i <= numt; i++) {
	template_c *t = DBI(*xx)->templates[i];

	if (!t || !(t->flags & TEMP_FLAG_SPANNING))
	    continue;

	check_template_length_overlap(DBI_io(*xx), t,
				      c1, c2, /*overlapLength,*/
				      ABS(offset));

	/* Spanning, but not between both c1 and c2 */
	if (!t->computed_length)
	    continue;

	if (t->consistency)
	    tbad++;
	else
	    tgood++;
    }

    *ptgood = tgood;
    *ptbad = tbad;
}

void countDisagreements(EdStruct *xx[2], int *overlapLength, int *wingeCount,
			int *ptgood, int *ptbad)
{
    int left0,right0;
    int left1;
    int length0,length1;
    int offset = editorLockedPos(xx, 1/*force recalculation*/);
    int i;
    char *ol0,*ol1;
    
    *ptgood = 0;
    *ptbad = 0;

    if (offset < 0) {
	left0 = 1-offset;
	left1 = 1;
    } else {
	left0 = 1;
	left1 = 1+offset;
    }
    length0 = DB_Length(xx[0],0);
    length1 = DB_Length(xx[1],0);
    right0 = (offset+length0 < length1)
	? length0
	: length1-offset;
    *overlapLength = right0 - left0+1;
    *wingeCount  = 0;
    
    if (*overlapLength > 0) {
	ol0 = (char *)xmalloc(*overlapLength+1);
	ol1 = (char *)xmalloc(*overlapLength+1);
	DBcalcConsensus(xx[0],left0,*overlapLength,ol0,NULL,BOTH_STRANDS);
	DBcalcConsensus(xx[1],left1,*overlapLength,ol1,NULL,BOTH_STRANDS);
	for (i=0;i<*overlapLength;i++) if(ol0[i]!=ol1[i])(*wingeCount)++;
	xfree(ol0);
	xfree(ol1);
    }

    spanning_template_stats(xx, *overlapLength, ptgood, ptbad);
}


/*
 * Find the first sequence that starts at or to the right of a
 * given position using a binary search
 */
int posToIndex(EdStruct *xx, int pos)
{
    int Min, Max, Mid;
    
    /* binary search */
    /* Min, Max, Mid refer to pairs of numbers: ie MAX --> [MAX-1],[MAX] */
    Min = 1;
    Max = DBI_gelCount(xx) + 1;
    
    do {
	int r1,r2;
	
	Mid = (Max+Min)/2;
	
	/* compare */
	r1 = (Mid==1)?(pos-1):DB_RelPos(xx,DBI_order(xx)[Mid-1]);
	r2 = (Mid==DBI_gelCount(xx)+1)?(pos+1):DB_RelPos(xx,DBI_order(xx)[Mid]);
	
	if (r1 < pos && r2 >= pos) 
	    return (Mid==DBI_gelCount(xx)+1)?(Mid-1):Mid;
	
	if (r1 < pos)
	    Min = Mid+1;
	else
	    Max = Mid-1;
	
    } while (Max>=Min);
    
    return 0;
}


/*
 * Find the first sequence that starts at or to the right of a
 * given position
 */
int posToSeq(EdStruct *xx, int pos)
{
    int ind;
    
    ind = posToIndex(xx,pos);
    if (ind)
	return DB_RelPos(xx,DBI_order(xx)[ind]);
    else
	return 0;
    
}


int seqToIndex(EdStruct *xx, int seq)
{
    int i;
    int ind;
    
    ind = posToIndex(xx,DB_RelPos(xx,seq));
    
    if (ind) {
	for (i=ind; i<=DBI_gelCount(xx) && DBI_order(xx)[i]!=seq;i++)
	    ;
	
	if (i<=DBI_gelCount(xx))
	    return i;
    }
    
    return 0;
}

/*
 * Converts a gel reading number to an editor 'seq' number. -1 for unknown
 */
int NumberToSeq(DBInfo *db, int number) {
    int i;

    for (i=1; i <= _DBI_gelCount(db); i++) {
	if (_DB_Number(db, i) == number)
	    return i;
    }

    return -1;
}

/*
 * Determine position in original sequence corresponding to pos. Pos is the
 * offset from the start of the entire sequence, not just the used sequence.
 */
int origpos(EdStruct *xx, int seq, int pos)
{
    int2 *opos;
    int i;
    int left=0, right=0;
    
    /*
     * ensure sequence and related information is read into memory
     */
    if (DBgetSeq(DBI(xx),seq) == NULL) return 0;

    if (pos < 1)
	pos = 1;
    if (pos > DB_Length2(xx, seq))
	pos = DB_Length2(xx, seq);

    opos = DB_Opos(xx,seq);

    /*
     * we have an immediate answer
     */
    if (opos[pos-1])
	return opos[pos-1];
    
    /*
     * A zero value means the base has been edited at some stage
     * We should look left and right for the first unedited base,
     * and average the original base positions of these.
     */
    for(i=pos-1;i>0 && (left=opos[i-1])==0;i--);
    for(i=pos+1;i<=DB_Length2(xx,seq) && (right=opos[i-1])==0;i++);

    if (!left)
	left = right;

    if (!right)
	right = left;

    /*
     * For consistency of trace display, round down for uncomplemented
     * sequences and up for complemented sequences.
     */
    if (DB_Comp(xx, seq) == UNCOMPLEMENTED) {
	return (left + right) / 2;
    } else {
	return (left + right) / 2.0 + 0.5;
    }
}



/*
 *----------------------------------------------------------------------------
 * General cutoff obtaining
 *---------------------------------------------------------------------------
 */

/*
 * get *width* right most characters of left cutoff
 */
void getLeftCutOff(EdStruct *xx,int seq, int width, char *str)
{
    /*
     * ensure sequence and related information is read into memory
     */
    if (DBgetSeq(DBI(xx),seq) == NULL) return;
    
    
    if (xx->reveal_cutoffs && width >0 ) {
	char *s = DB_Seq(xx,seq);
	
	if (s != NULL) {
	    int l = DB_Start(xx,seq);
	    for (;l<width;width--)*str++=' ';
	    strncpy(str,&s[l-width],width);
	    return;
	}
    }
    
    for(;width>0;width--)*str++=' ';
}


/*
 * get *width* characters starting at position in left cutoff
 */
void getLCut(EdStruct *xx,int seq, int pos, int width, char *str)
{
    /*
     * ensure sequence and related information is read into memory
     */
    if (DBgetSeq(DBI(xx),seq) == NULL) return;
    
    if (xx->reveal_cutoffs && width >0 ) {
	char *s = DB_Seq(xx,seq);
	
	if (s != NULL) {
	    int l = DB_Start(xx,seq);
	    for (;l<pos;pos--,width--)*str++=' ';
	    strncpy(str,&s[l-pos],width);
	    return;
	}
    }
    
    for(;width>0;width--)*str++=' ';
}





/*
 * get *width* right most characters of left cutoff
 */
void getRightCutOff(EdStruct *xx,int seq, int width, char *str)
{
    /*
     * ensure sequence and related information is read into memory
     */
    if (DBgetSeq(DBI(xx),seq) == NULL) return;
    
    if (xx->reveal_cutoffs && width >0 ) {
	char *s = DB_Seq(xx,seq) + DB_End(xx,seq);
	
	if (s != NULL) {
	    int l = DB_Length2(xx,seq) - DB_End(xx,seq) + 1;
	    for (;l<width;width--)str[width-1]=' ';
	    strncpy(str,s,width);
	    return;
	}
    }
    
    for(;width>0;width--)*str++=' ';
}


/*
 * get *width* characters starting at position in left cutoff
 */
void getRCut(EdStruct *xx,int seq, int pos, int width, char *str)
{
    /*
     * ensure sequence and related information is read into memory
     */
    if (DBgetSeq(DBI(xx),seq) == NULL) return;
    
    if (xx->reveal_cutoffs && width >0 ) {
	char *s = DB_Seq(xx,seq) + DB_End(xx,seq) - 1;
	
	if (s != NULL) {
	    int l = DB_Length2(xx,seq) - DB_End(xx,seq) + 1;
	    for (;l<pos+width;width--)str[width-1]=' ';
	    strncpy(str,&s[pos],width);
	    return;
	}
    }
    
    for(;width>0;width--)*str++=' ';
}


int lenRCut(EdStruct *xx, int seq)
{
    
    /*
     * ensure sequence and related information is read into memory
     */
    if (DBgetSeq(DBI(xx),seq) == NULL) return 0;
    
    return DB_Length2(xx,seq) - DB_End(xx,seq) + 1;
}


int lenLCut(EdStruct *xx, int seq)
{
    /*
     * ensure sequence and related information is read into memory
     */
    if (DBgetSeq(DBI(xx),seq) == NULL) return 0;
    
    return DB_Start(xx,seq);
}

/*
 * Finds the quality value for a base, allowing for this base being a pad
 * (whereby we average from the surrounding bases).
 */
int getQual(EdStruct *xx, int seq, int pos) {
    int score = -1, i, iend;
    char *s = DBgetSeq(DBI(xx), seq);
    int1 *conf = DB_Conf(xx, seq) + DB_Start(xx, seq);

    if (s[pos-1] != '*')
	return conf[pos-1];

    /* average left and right bases */
    for (i = pos-2, iend = -DB_Start(xx, seq); i >= iend; i--) {
	if (s[i] != '*') {
	    score = conf[i];
	    break;
	}
    }

    for (i = pos, iend = DB_Length(xx, seq) - DB_Start(xx, seq);
	 i < iend; i++) {
	if (s[i] != '*') {
	    if (score != -1)
		score = (score + conf[i]) / 2;
	    else
		score = conf[i];
	    break;
	}
    }

    return score;
}

/*
 * Find the current maximum sequence length in the editor.
 * If 'clipped' is true then we only look at the used quality-clipped portion,
 * otherwise we consider the hidden data too.
 */
int dbi_max_gel_len(DBInfo *db, int clipped) {
    int i, end = _DBI_gelCount(db);
    int max_len = 0;

    if (clipped) {
	for (i = 1; i <= end; i++) {
	    if (max_len < _DB_Length(db, i))
		max_len = _DB_Length(db, i);
	}
    } else {
	for (i = 1; i <= end; i++) {
	    if (max_len < _DB_Length2(db, i))
		max_len = _DB_Length2(db, i);
	}
    }

    return max_len;
}

/* Returns the character underneath the editor cursor */
int edGetChar(EdStruct *xx) {
    char *seq = DBgetSeq(DBI(xx), xx->cursorSeq);
    return seq ? seq[xx->cursorPos-1] : 0;
}

