/*
 *----------------------------------------------------------------------------
 * The interface routines - hooks for the outside world
 *
 * These routines are mappings of events to edit functions
 *
 * All return:
 *    0 - OK
 *    1 - Error
 *----------------------------------------------------------------------------
 */

#include <string.h>

#include "gap-dbstruct.h"
#include "edUtils.h"
#include "contigEditor.h"
#include "extend.h"
#include "active_tags.h"
#include "misc.h"
#include "io-reg.h"
#include "template_display.h" /* IMPORT: update_reading_list */
#include "gap_globals.h"
#include "qualIO.h"
#include "notes.h"
#include "tman_interface.h"

#define REPLACE (0L)
#define INSERT  (1L << 0)
#define CONSENSUS (0L)
#define READING (1L << 1)

/*
 *----------------------------------------------------------------------------
 * Locally required functions
 *---------------------------------------------------------------------------
 */

static void caretDown2(EdStruct *xx, int i, int seqCount, int *seqp, int *posp,
		       int *seqList, int posInContig) {
    int pos, seq;

    if (xx->editorState == StateDown)
	return;

    i++;
    if (i == seqCount)
	i = 0;

    
    seq = seqList[i];
    pos = posInContig - DB_RelPos(xx,seq) + 1;

    if (pos < -DB_Start(xx, seq) + 1) {
	caretDown2(xx, i, seqCount, seqp, posp, seqList, posInContig);
	return;
    } else if (pos > DB_Length2(xx, seq) - DB_Start(xx, seq) + 1) {
	caretDown2(xx, i, seqCount, seqp, posp, seqList, posInContig);
	return;
    }

    *seqp = seq;
    *posp = pos;
}

static void caretUp2(EdStruct *xx, int i, int seqCount, int *seqp, int *posp,
		       int *seqList, int posInContig) {
    int pos, seq;

    if (xx->editorState == StateDown)
	return;

    if (i==0) i = seqCount;
    i--;
    
    seq = seqList[i];
    pos = posInContig - DB_RelPos(xx,seq) + 1;

    if (pos < -DB_Start(xx, seq) + 1) {
	caretUp2(xx, i, seqCount, seqp, posp, seqList, posInContig);
	return;
    } else if (pos > DB_Length2(xx, seq) - DB_Start(xx, seq) + 1) {
	caretUp2(xx, i, seqCount, seqp, posp, seqList, posInContig);
	return;
    }

    *seqp = seq;
    *posp = pos;
}


static int validKey(EdStruct *xx, char key)
{
/*
 * Staden Codes:
 *
 * static char validKeys1[] = "CcTtAaGg1234DdVvBbHhKkLlMmNnRrYy5678-*";
 * static char validKeys2[] = "ctag1234dvbhklmnry5678-*";
 */

    /* IUB codes */
    static char validKeys1[] = "CcTtUuAaGg-*RrYyMmKkSsWwBbDdHhVvNn";
    static char validKeys2[] = "ctuag-*rymkswbdhvn";
    char *validKeys;

    if (DBI_io(xx)->db.data_class == GAP_PROTEIN) {
	return (key!=' ');
    } else {
	validKeys = (xx->super_edit & SUPEREDIT_UPPERCASE)
	    ? validKeys1 : validKeys2;
	return (strchr(validKeys,key) != NULL && key != '\0');
    }
}


/*
 *----------------------------------------------------------------------------
 * Globally accessable functions
 *---------------------------------------------------------------------------
 */

/* ---- from edUtils.c ---- */

/*
 * swap two characters and move the cursor left
 */
int edTransposeLeft(EdStruct *xx, int num_bases) {
    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    RedisplaySeq(xx, xx->cursorSeq);
    return transpose(xx, xx->cursorSeq, xx->cursorPos, -1, num_bases);
}

/*
 * swap two characters and move the cursor right
 */
int edTransposeRight(EdStruct *xx, int num_bases) {
    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    RedisplaySeq(xx, xx->cursorSeq);
    return transpose(xx, xx->cursorSeq, xx->cursorPos, 1, num_bases);
}

/*
 * Handle a key press - typically anything that's not handled elsewhere,
 * which boils down to edits.
 *
 * If nomove is set, we do not move the cursor to the next base once an
 * edit is made.
 */
int edKeyPress(EdStruct *xx, char key, int nomove)
{
    int err;
    int mode;
    int bases;
    
    int seq, pos;
    
    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    seq = xx->cursorSeq;
    pos = xx->cursorPos;
    if (! onScreen(xx, seq, pos, NULL)) {
	/*
	 * When cursor is currently off screen,
	 * give a warning tone, and recentre screen on cursor
	 */
	showCursor(xx, seq, pos);
	return 1;
    }
    
    /* determine mode */
    mode = ((seq) ? READING : CONSENSUS)
	| ((editModeIsInsert(xx)) ? INSERT : REPLACE);
    
    /* spaces have special meaning */
    if (key == ' ') {
	/* spaces inserted at start of a sequence mean shiftRight */
	if (seq && 
	    ((xx->reveal_cutoffs && pos == -DB_Start(xx, seq) + 1) ||
	     (!xx->reveal_cutoffs && pos == 1)) &&
	    xx->super_edit & SUPEREDIT_SHIFT_READ) {
	    openUndo(DBI(xx));
	    err = shiftRight(xx,seq, 1);
	    closeUndo(xx, DBI(xx));

	    getExtents(xx);
	    if (err) {
		return 1;
	    } else {
		redisplayWithCursor(xx);
		return 0;
	    }
	} else {
	    return 1;
	}
    }
    
    /* is mode compatible with update privilages */
    if (mode == (READING|INSERT) && !(xx->super_edit & SUPEREDIT_INS_READ)) {
	return 1;
    }
    
    /* otherwise - check for valid key */
    if (! validKey(xx, key)) {
	return 1;
    }
    
    /* start recording undos for command */
    openUndo(DBI(xx));
    
    bases = 0;/* init to keep to the compiler happy - not really needed */
    switch (mode) {
    case (CONSENSUS|INSERT):
	if (key != '*' && !(xx->super_edit & SUPEREDIT_INS_ANY_CON)) {
	    err = 1;
	    break;
	}
	err = insertBasesConsensus(xx,pos,1,&key);
	bases = 1; nomove = 0;
	break;
    case (CONSENSUS|REPLACE):
	if (key == '*') {
	    err = insertBasesConsensus(xx,pos,1,&key);
	    bases = 1; nomove = 0;
	} else {
	    if (!(xx->super_edit & SUPEREDIT_REPLACE_CON)) {
		err = 1;
		break;
	    }
	    err = replaceBasesConsensus(xx,pos,1,&key);
	    bases = 1;
	}
	break;
    case (READING|INSERT):
	bases = insertBases(xx,seq,pos,1,&key);
	err = (bases==0);
	nomove = 0;
	break;
    case (READING|REPLACE):
	bases = replaceBases(xx,seq,pos,1,&key);
	err = (bases==0);
	break;
    default:
	err = 1;
	break;
    }
    
    
    if (err) {
	closeUndo(xx, DBI(xx));
	return 1;
    }
    
    if (!nomove)
	U_adjust_cursor(xx,bases);
    
    /* stop recording undos for command */
    closeUndo(xx, DBI(xx));
    
    redisplayWithCursor(xx);

    getExtents(xx);
    return 0;
}


/*
 * Handle a delete key press
 */
int edKeyDeleteCommon(EdStruct *xx, int to_left)
{
    int err;
    int mode;
    int bases;
    char consensus[2];
    
    int seq, pos;
    
    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    seq = xx->cursorSeq;
    pos = xx->cursorPos;
    if (! onScreen(xx, seq, pos, NULL)) {
	/*
	 * When cursor is currently off screen,
	 * give a warning tone, and recentre screen on cursor
	 */
	showCursor(xx, seq, pos);
	return 1;
    }
    
    /* deletions at the left end of a reading have special meaning */
    if ((xx->reveal_cutoffs && pos == -DB_Start(xx, seq) + 1) ||
	(!xx->reveal_cutoffs && pos == 1) && 
	(xx->super_edit & SUPEREDIT_SHIFT_READ)) {
	/* can only shift left a gel reading - not the consensus */
	if (seq) {
	    openUndo(DBI(xx));
	    err = shiftLeft(xx,seq, 1);
	    closeUndo(xx, DBI(xx));
	    getExtents(xx);
	    
	    if (err) {
		return 1;
	    } else {
		redisplayWithCursor(xx);
		return 0;
	    }
	} else {
	    return 1;
	}
    }
    
    mode = ((seq)? READING : CONSENSUS );
    
    /* start recording undos for command */
    openUndo(DBI(xx));
    
    switch (mode) {
    case (CONSENSUS):
	/* currently only allow deletion of pads */
	DBcalcConsensus(xx,pos-1,1,consensus,NULL,BOTH_STRANDS);
	if (*consensus != '*') {
	    if (!(xx->super_edit & SUPEREDIT_DEL_ANY_CON)) {
		if (*consensus != '-' ||
		    !(xx->super_edit & SUPEREDIT_DEL_DASH_CON)) {
		    err = 1;
		    break;
		}
	    }
	}
	err = deleteBasesConsensus(xx,pos-1,1);
	break;
    case (READING):
	if (xx->super_edit & SUPEREDIT_DEL_READ &&
	    (xx->reveal_cutoffs || pos > 1)) {
	    bases = deleteBases(xx,seq,pos-1,1);
	    err = (bases==0);
	} else
	    err = 1;
	break;
    default:
	err = 1;
	break;
    }
    
    if (err) {
	closeUndo(xx, DBI(xx));
	return 1;
    }

    if (to_left && seq && pos > 1)
      shiftRight(xx,seq, 1);

    U_adjust_cursor(xx,-1);
	
    /* stop recording undo for command */
    closeUndo(xx, DBI(xx));
	
    redisplayWithCursor(xx);
    getExtents(xx);

    return 0;
}

/*
 * Delete the base to the left of the cursor, shifting the data to the
 * rightwards left by one base to keep the reading left end at the
 * same position.
 */
int edKeyDelete(EdStruct *xx)
{
    return edKeyDeleteCommon(xx, 0);
}

/*
 * Delete the base to the left of the cursor, shifting the reading right
 * by one.
 */
int edKeyDeleteLeft(EdStruct *xx)
{
    return edKeyDeleteCommon(xx, 1);
}


/*
 * Sets the confidence value to be 100%
 */
int edConf100(EdStruct *xx) {
    int seq, pos;

    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    seq = xx->cursorSeq;
    pos = xx->cursorPos;

    if (! onScreen(xx, seq, pos, NULL)) {
	/*
	 * When cursor is currently off screen,
	 * give a warning tone, and recentre screen on cursor
	 */
	showCursor(xx, seq, pos);
	return 0;
    }

    if (pos > DB_Length2(xx, seq) - DB_Start(xx, seq) ||
	pos <= -DB_Start(xx, seq)) {
	bell();
	return 1;
    }

    if (adjustBaseConf(xx, seq, pos, 100, 1)) {
	bell();
	return 1;
    }
    
    return 0;
}


/*
 * Sets the confidence value to be 0%
 */
int edConf0(EdStruct *xx) {
    int seq, pos;

    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    seq = xx->cursorSeq;
    pos = xx->cursorPos;

    if (! onScreen(xx, seq, pos, NULL)) {
	/*
	 * When cursor is currently off screen,
	 * give a warning tone, and recentre screen on cursor
	 */
	showCursor(xx, seq, pos);
	return 0;
    }

    if (pos > DB_Length2(xx, seq) - DB_Start(xx, seq) ||
	pos <= -DB_Start(xx, seq)) {
	bell();
	return 1;
    }

    if (adjustBaseConf(xx, seq, pos, 0, 1)) {
	bell();
	return 1;
    }

    return 0;
}


/*
 * Increments the confidence of a base by a specific amount. Use negatives
 * to decrememnt
 */
int edConfIncr(EdStruct *xx, int amount) {
    int seq, pos;
    int1 *old_conf;
    signed int cval;

    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    seq = xx->cursorSeq;
    pos = xx->cursorPos;

    if (!seq) {
	bell();
	return 1;
    }

    if (! onScreen(xx, seq, pos, NULL)) {
	/*
	 * When cursor is currently off screen,
	 * give a warning tone, and recentre screen on cursor
	 */
	showCursor(xx, seq, pos);
	return 0;
    }

    if (pos > DB_Length2(xx, seq) - DB_Start(xx, seq) ||
	pos <= -DB_Start(xx, seq)) {
	bell();
	return 1;
    }

    (void) DBgetSeq(DBI(xx),seq);
    old_conf = DB_Conf(xx,seq) + DB_Start(xx,seq);
    cval = old_conf[pos-1];
    if ((cval == 100 && amount > 0) ||
	(cval == 0 && amount < 0)) {
	bell();
	return 1;
    }
    cval += amount;
    if (cval > 100) cval = 100;
    if (cval < 0) cval = 0;
    if (adjustBaseConf(xx, seq, pos, cval, 0)) {
	bell();
	return 1;
    }

    return 0;
}

/*
 * updates the global Readings list with selections made in the contig editor
 */
void
update_reading_list(GapIO *io,
		    int r_num,
		    int high_light)
{
    char cmd[1024];
    char *r_name;

    /* can't do the consensus */
    if (r_num > 0) {
	r_name = get_read_name(io, r_num);
	sprintf(cmd, "UpdateReadingListItem %s %d", r_name, high_light);
	Tcl_Eval(GetInterp(), cmd);
    }
}

/*
 * State is -1 to toggle
 *           0 to clear
 *           1 to set
 */
int edSelectRead(EdStruct *xx, int seq, int state) {
    int flag;
    reg_highlight_read rh;

    if (xx->editorState == StateDown)
	return 1;

    flag = DB_Flags(xx, seq);
    if (state == -1) {
	flag ^= DB_FLAG_SELECTED;
    } else if (state == 0) {
	flag &= ~DB_FLAG_SELECTED;
    } else {
	flag |= DB_FLAG_SELECTED;
    }

    DBsetFlags(xx, seq ,flag);
    RedisplayName(xx, seq);
    redisplaySequences(xx, 1);

    rh.job = REG_HIGHLIGHT_READ;
    rh.seq = DB_Number(xx, seq);
    rh.val = (flag & DB_FLAG_SELECTED) ? 1 : 0;

    update_reading_list(DBI_io(xx), rh.seq, rh.val);
    contig_notify(DBI_io(xx), DBI_contigNum(xx), (reg_data *)&rh);
    
    return 0;
}


int edStartRead(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return 1;

    setCursorPos(xx, 1);
    showCursor(xx, xx->cursorSeq, xx->cursorPos);
    repositionTraces(xx);

    return 0;
}


int edEndRead(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return 1;

    setCursorPos(xx, DB_Length(xx, xx->cursorSeq) + 1);
    showCursor(xx, xx->cursorSeq, xx->cursorPos);
    repositionTraces(xx);

    return 0;
}


int edStartRead2(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return 1;

    if (xx->reveal_cutoffs) {
	setCursorPos(xx, 1-DB_Start(xx, xx->cursorSeq));
    } else {
	setCursorPos(xx, 1);
    }

    showCursor(xx, xx->cursorSeq, xx->cursorPos);
    repositionTraces(xx);

    return 0;
}


int edEndRead2(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return 1;

    if (xx->reveal_cutoffs) {
	setCursorPos(xx, DB_Length2(xx, xx->cursorSeq) -
		     DB_Start(xx, xx->cursorSeq) + 1);
    } else {
	setCursorPos(xx, DB_Length(xx, xx->cursorSeq) + 1);
    }

    showCursor(xx, xx->cursorSeq, xx->cursorPos);
    repositionTraces(xx);

    return 0;
}


int edStartContig(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return 1;

    setCursorPosSeq(xx, 1, 0);

    showCursor(xx, xx->cursorSeq, xx->cursorPos);
    repositionTraces(xx);

    return 0;
}


int edEndContig(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return 1;

    setCursorPosSeq(xx, DB_Length(xx, 0) + 1, 0);

    showCursor(xx, xx->cursorSeq, xx->cursorPos);
    repositionTraces(xx);

    return 0;
}


/*
 * Move cursor right
 */
int edCursorRight(EdStruct *xx)
{
    if (xx->editorState == StateDown)
	return 1;

    if (xx->cursorPos <= DB_Length(xx,xx->cursorSeq) ||
	(xx->reveal_cutoffs && xx->cursorPos + DB_Start(xx, xx->cursorSeq) <=
	 DB_Length2(xx, xx->cursorSeq))) {
	setCursorPos(xx, xx->cursorPos+1);
        showCursor(xx,xx->cursorSeq, xx->cursorPos);
	repositionTraces(xx);
	return 0;
    } else {
	bell();
	return 1;
    }
}


/*
 * Move cursor left
 */
int edCursorLeft(EdStruct *xx)
{
    if (xx->editorState == StateDown)
	return 1;

    if (xx->cursorPos > 1 ||
	(xx->reveal_cutoffs && xx->cursorPos> -DB_Start(xx, xx->cursorSeq)+1)){
	setCursorPos(xx, xx->cursorPos-1);
        showCursor(xx,xx->cursorSeq, xx->cursorPos);
	repositionTraces(xx);
	return 0;
    } else {
	bell();
	return 1;
    }
}


/*
 * Move cursor down,
 * cycle if necessary
 */
int edCursorDown(EdStruct *xx)
{
    int *seqList,seqCount, cseq, cpos;
    int posInContig;
    int i;
    
    if (xx->editorState == StateDown)
	return 1;

    posInContig = positionInContig(xx,xx->cursorSeq,xx->cursorPos);
    seqList = sequencesInRegion(xx,posInContig-1,2);
    seqCount = linesInRegion(xx,posInContig-1,2);

    if (seqCount == 1)
	return 0;

    for(i=0;
	i<seqCount && seqList[i]!=xx->cursorSeq;
	i++);
    
    cseq = xx->cursorSeq;
    cpos = xx->cursorPos;
    caretDown2(xx, i, seqCount, &cseq, &cpos, seqList, posInContig);
    if (xx->cursorSeq != cseq || xx->cursorPos != cpos)
	setCursorPosSeq(xx, cpos, cseq);

    showCursor(xx, xx->cursorSeq, xx->cursorPos);

    return 0;
}


/*
 * Move cursor up,
 * cycle if necessary
 */
int edCursorUp(EdStruct *xx)
{
    int *seqList,seqCount, cseq, cpos;
    int posInContig;
    int i;
    
    if (xx->editorState == StateDown)
	return 1;

    posInContig = positionInContig(xx,xx->cursorSeq,xx->cursorPos);
    seqList = sequencesInRegion(xx,posInContig-1,2);
    seqCount = linesInRegion(xx,posInContig-1,2);
    if (seqCount == 1)
	return 0;
    for(i=0;
	i<seqCount && seqList[i]!=xx->cursorSeq;
	i++);
    
    cseq = xx->cursorSeq;
    cpos = xx->cursorPos;
    caretUp2(xx, i, seqCount, &cseq, &cpos, seqList, posInContig);
    if (xx->cursorSeq != cseq || xx->cursorPos != cpos)
	setCursorPosSeq(xx, cpos, cseq);

    showCursor(xx,xx->cursorSeq, xx->cursorPos);

    return 0;
}    


/* ---- from extend.c ---- */

/*
 * Extends the left cutoff by one base
 */
int edExtendLeft(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    RedisplaySeq(xx, xx->cursorSeq);
    return meta_arrow(xx, EXTEND_LEFT);
}


/*
 * Extends the right cutoff by one base
 */
int edExtendRight(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    RedisplaySeq(xx, xx->cursorSeq);
    return meta_arrow(xx, EXTEND_RIGHT);
}


/*
 * Sets the left cutoff position
 */
int edZapLeft(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    RedisplaySeq(xx, xx->cursorSeq);
    zap_Left(xx);
    
    return 0;
}


/*
 * Sets the right cutoff position
 */
int edZapRight(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    RedisplaySeq(xx, xx->cursorSeq);
    zap_Right(xx);

    return 0;
}


/*
 * Set the cursor position from an X,Y position within the sheet widget
 */
int edSetCursor(EdStruct *xx, int x, int y) {
    int *seqList;

    if (xx->editorState == StateDown)
	return 1;

    if (y < 0 || y >= xx->displayHeight || x < 0 || x >= xx->displayWidth)
	return 1;

    /* Special case for cons. when we've got more sequences 'underneath' it */
    y /= xx->lines_per_seq;
    if (y == (xx->displayHeight-1)/xx->lines_per_seq)
	y = (xx->totalHeight-1)/xx->lines_per_seq;
    else
	y += xx->displayYPos;

    seqList = sequencesOnScreen(xx, xx->displayPos, xx->displayWidth);
    setCursorPosSeq(xx,
		    xx->displayPos - DB_RelPos(xx,seqList[y]) + x + 1,
		    seqList[y]);

    if (!(xx->set &&
	  xx->set_collapsed &&
	  xx->set[seqList[y]] &&
	  xx->set_collapsed[xx->set[seqList[y]]])) {
	if (xx->reveal_cutoffs) {
	    if (xx->cursorPos < -DB_Start(xx, xx->cursorSeq) + 1)
		setCursorPos(xx, -DB_Start(xx, xx->cursorSeq) + 1);
	    else
		if (xx->cursorPos + DB_Start(xx, xx->cursorSeq) >
		    DB_Length2(xx,xx->cursorSeq))
		    setCursorPos(xx, DB_Length2(xx,xx->cursorSeq)
				 - DB_Start(xx, xx->cursorSeq) + 1);
	} else {
	    if (xx->cursorPos<1)
		setCursorPos(xx, 1);
	    else
		if (xx->cursorPos > DB_Length(xx,xx->cursorSeq)+1)
		    setCursorPos(xx, DB_Length(xx,xx->cursorSeq)+1);
	}
    }

    positionCursor(xx, xx->cursorSeq, xx->cursorPos);
    repositionTraces(xx);

    return 0;
}

/*
 * Sets the cursor to an absolute position on the consensus
 */
void edSetCursorConsensus(EdStruct *xx, int pos) {
    if (xx->editorState == StateDown)
	return;

    if (pos < 1)
	pos = 1;
    else if (pos > DB_Length(xx, 0) + 1)
	pos = DB_Length(xx, 0) + 1;

    setCursorPosSeq(xx, pos, 0);
	
    positionCursor(xx, 0, pos);
    showCursor(xx, 0, pos);
    repositionTraces(xx);

    return;
}


/* ---- from contigEditor.c ---- */

/*
 * Set reveal cutoffs, -1 is toggle
 */
int edSetRevealCutoffs(EdStruct *xx, int reveal) {
    if (xx->editorState == StateDown)
	return 1;

    if (reveal == -1)
	xx->reveal_cutoffs ^= 1;
    else
	xx->reveal_cutoffs = reveal;

    /*
     * Move cursor onto consensus if we're removing cutoff displays and
     * the cursor is currently in the cutoff.
     */
    if (!xx->reveal_cutoffs && 
	(xx->cursorPos < 1 ||
	 xx->cursorPos > DB_Length(xx, xx->cursorSeq) + 1)) {
	setCursorPosSeq(xx,
			positionInContig(xx, xx->cursorSeq, xx->cursorPos),
			0);

	if (xx->cursorPos < 1)
	    setCursorPos(xx, 1);
	else if (xx->cursorPos > DB_Length(xx, 0) + 1)
	    setCursorPos(xx, DB_Length(xx, 0) + 1);
    }
    
    getExtents(xx);
    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 0);
    ed_set_slider_pos(xx, xx->displayPos);

    return 0;
}


/*
 * Set insert mode. -1 is toggle
 */
int edSetInsertMode(EdStruct *xx, int insert) {
    if (insert == -1)
	xx->insert_mode ^= 1;
    else {
	if (insert)
	    xx->insert_mode = 1;
	else
	    xx->insert_mode = 0;
    }

    return 0;
}


/*
 * Set ruler display mode - padded(0) vs unpadded(1).
 * ispadded == -1 implies toggle.
 */
int edSetRulerMode(EdStruct *xx, int ispadded) {
    if (ispadded == -1)
	xx->unpadded_ruler ^= 1;
    else {
	xx->unpadded_ruler = ispadded;
    }

    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 0);

    return 0;
}


/*
 * Controls whether reading names (0) or template names (1) are shown
 * in the edNames component.
 */
int edSetTemplateNameMode(EdStruct *xx, int val) {
    xx->template_names = val;

    xx->refresh_flags |= ED_DISP_NAMES;
    redisplaySequences(xx, 0);

    return 0;
}
int edGetTemplateNameMode(EdStruct *xx) {
    return xx->template_names;
}


/*
 * Get gel number from position. Returns the gel number or -1 for failure.
 */
int edGetGelNumber(EdStruct *xx, int x, int y) {
    int *seqList;

    if (xx->editorState == StateDown)
	return -1;

    if (y < 0 || y >= xx->displayHeight ||
	x < 0 || x >= xx->displayWidth)
	return -1;

    /* Special case for cons. when we've got more sequences 'underneath' it */
    if (y == xx->displayHeight - 1)
	return 0;

    y /= xx->lines_per_seq;
    y += xx->displayYPos;

    seqList = sequencesOnScreen(xx, xx->displayPos, xx->displayWidth);
    return seqList[y];
}


/*
 * Get gel name from number
 */
char *edGetGelName(EdStruct *xx, int number) {
    if (xx->editorState == StateDown)
	return NULL;

    return DBgetName(DBI(xx), number);
}

/*
 * Get all gel names from a number and sequences to the right of this number.
 *
 * Returns a dstring holding the list
 *        or NULL for failure
 */
dstring_t *edGetGelNamesToRight(EdStruct *xx, int number) {
    dstring_t *ds = NULL;
    int seq;
    int pos = DB_RelPos(xx, number);

    ds = dstring_create(NULL);
    for (seq = 1; seq <= DBI_gelCount(xx); seq++) {
	if (DB_RelPos(xx, seq) > pos ||
	    (DB_RelPos(xx, seq) == pos && seq >= number)) {
	    dstring_appendf(ds, "{%s} ", DBgetName(DBI(xx), seq));
	}
    }

    return ds;
}


/*
 * Set consensus mode
 */
void edSetCMode(EdStruct *xx, int value) {
    if (xx->editorState == StateDown)
	return;

    xx->consensus_mode = value;
    xx->refresh_flags |= ED_DISP_CONS | ED_DISP_STATUS;
    invalidate_consensus(xx);
    redisplaySequences(xx, 0);
}


/*
 * Set consensus cutoff
 */
void edSetCCutoff(EdStruct *xx, int value) {
    if (xx->editorState == StateDown)
	return;

    xx->con_cut = ((float)value)/100;
    xx->refresh_flags |= ED_DISP_CONS | ED_DISP_STATUS;
    invalidate_consensus(xx);
    redisplaySequences(xx, 0);
}


/*
 * Set quality cutoff
 */
void edSetQCutoff(EdStruct *xx, int value) {
    if (xx->editorState == StateDown)
	return;

    xx->qual_cut = value;
    xx->refresh_flags |= ED_DISP_SEQS | ED_DISP_STATUS;
    invalidate_consensus(xx);
    redisplaySequences(xx, 0);
}


/*
 * Set highlight disagreements quality cutoff
 */
void edSetDifferenceQuality(EdStruct *xx, int value) {
    if (xx->editorState == StateDown)
	return;

    xx->diff_qual = value;
    xx->refresh_flags |= ED_DISP_SEQS;
    redisplaySequences(xx, 0);
}


/*
 * Set or clear the locking of join editors. Value of 1 means they are locked.
 * 0 means unlock.
 */
int edSetJoinLock(EdStruct *xx, int value) {
    EdStruct *xx0 = xx->link->xx[0];
    EdStruct *xx1 = xx->link->xx[1];

    if (xx->editorState == StateDown)
	return -1;

    if (xx->editorMode != JOINMODE)
	return -1;

    if (xx->link->locked = value) {
	xx->link->lockOffset = xx1->displayPos - xx0->displayPos;
    }

    getExtents(xx0);
    getExtents(xx1);

    xx0->refresh_flags |= ED_DISP_ALL;
    xx1->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx0, 0);
    redisplaySequences(xx1, 0);

    return 0;
}


/*
 * Make a join
 */
void edJoin(EdStruct *xx) {
    if (xx->editorState == StateDown)
	return;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return;
    }
    
    if (xx->editorMode != JOINMODE || !xx->link)
	return;

    joinDB(xx->link->xx, DBI_io(xx));
}


/*
 * Perform an alignment in the join editor
 */
void edJoinAlign(EdStruct *xx, int fixed_left, int fixed_right) {
    if (xx->editorState == StateDown)
	return;

    if (!xx->link)
	return;

    if (alignOverlap(xx->link->xx, fixed_left, fixed_right)) {
	bell();
    } else {
	xx->link->xx[0]->refresh_flags |= ED_DISP_ALL;
	xx->link->xx[1]->refresh_flags |= ED_DISP_ALL;

	/* Force the cursor to be in the correct position */
	xx->link->lockOffset = xx->link->xx[1]->displayPos -
	    xx->link->xx[0]->displayPos;
	setCursorPos(xx, xx->cursorPos);

	redisplaySequences(xx->link->xx[0], 1);
	redisplaySequences(xx->link->xx[1], 1);
    }
}


/*
 * Sets the active annotation types for the contig editor
 */
void edSetActiveAnnos(EdStruct *xx, int argc, char **argv) {
    int i;

    if (xx->editorState == StateDown)
	return;

    for (i = 0; i < tag_db_count; i++)
	xx->tag_list[i] = 0;

    for (i = 0; i < argc; i++) {
	xx->tag_list[idToIndex(argv[i])] = 1;
    }

    edSelectClear(xx);
    xx->refresh_flags |= ED_DISP_SEQS;
    redisplaySequences(xx, 0);
}

/*
 * Returns the active annotation types
 */
char *edGetActiveAnnos(EdStruct *xx) {
    static char types[8192], *cp;
    int i;

    *types = 0;

    if (xx->editorState == StateDown)
	return types;

    cp = types;
    for (i = 0; i < tag_db_count; i++) {
	if (xx->tag_list[i]) {
	    strcpy(cp, indexToId(i));
	    cp[4] = ' ';
	    cp += 5;
	}
    }
    *cp = 0;

    return types;
}


/* ----- status line manipulation ----- */


/*
 * Turns on an editor status. If it's already on then we do nothing.
 */
void edStatusAdd(EdStruct *xx, int mode) {
    if (xx->editorState == StateDown)
	return;

    if (mode < 0 || mode >= MAX_STATUS_LINES)
	return;

    xx->status[mode] = 1;

    xx->refresh_flags |= ED_DISP_STATUS | ED_DISP_SCROLL | ED_DISP_HEIGHT;
    redisplaySequences(xx, 0);
}


/*
 * Turns off an editor status. If it's already off then we do nothing.
 */
void edStatusDelete(EdStruct *xx, int mode) {
    if (xx->editorState == StateDown)
	return;

    if (mode < 0 || mode >= MAX_STATUS_LINES)
	return;

    xx->status[mode] = 0;

    xx->refresh_flags |= ED_DISP_STATUS | ED_DISP_SCROLL | ED_DISP_HEIGHT;
    redisplaySequences(xx, 0);
}


/*
 * Sets the codon translation mode. Set to either 1 or 3 characters.
 */
void edStatusTransMode(EdStruct *xx, int mode) {
    if (xx->editorState == StateDown)
	return;

    xx->trans_mode = mode == 3 ? 3 : 1;

    xx->refresh_flags |= ED_DISP_STATUS;
    redisplaySequences(xx, 0);
}

/*
 * Controls whether to show the quality of bases as grey scale levels
 */
void edShowQuality(EdStruct *xx, int mode) {
    if (xx->editorState == StateDown)
	return;

    xx->show_qual = mode;

    xx->refresh_flags |= ED_DISP_READS;
    redisplaySequences(xx, 0);
}

void edShowCQuality(EdStruct *xx, int mode) {
    if (xx->editorState == StateDown)
	return;

    xx->show_cons_qual = mode;

    xx->refresh_flags |= ED_DISP_CONS;
    redisplaySequences(xx, 0);
}


/*
 * Controls whether to show the edits as pseudo-tags
 */
void edShowEdits(EdStruct *xx, int mode) {
    if (xx->editorState == StateDown)
	return;

    xx->show_edits = mode;

    xx->refresh_flags |= ED_DISP_READS;
    redisplaySequences(xx, 0);
}


/**
 * Hides a reading from the editor and consensus.
 * If seq is negative then this implies hide +seq and all sequences to the
 * right of this point.
 * Otherwise we just hide the single seq itself.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int edHideRead(EdStruct *xx, int xseq, int check_cursor) {
    int seq = ABS(xseq);

    if (check_cursor && ! onScreen(xx, seq, xx->cursorPos, NULL)) {
	/*
	 * When cursor is currently off screen,
	 * give a warning tone, and recentre screen on cursor
	 */
	showCursor(xx, seq, xx->cursorPos);
	bell();
	return 1;
    }

    if (seq == 0)
        return 1;

    if (xseq < 0) {
	int pos = DB_RelPos(xx, seq);
	for (seq = 1; seq <= DBI_gelCount(xx); seq++) {
	    if (DB_RelPos(xx, seq) > pos ||
		(DB_RelPos(xx, seq) == pos && seq >= -xseq)) {
		DBsetFlags(xx, seq, DB_Flags(xx, seq) ^ DB_FLAG_INVIS);
	    }
	}
	xx->refresh_flags |= ED_DISP_ALL;
    } else {
	DBsetFlags(xx, seq, DB_Flags(xx, seq) ^ DB_FLAG_INVIS);
    }

    RedisplayName(xx, seq);
    xx->refresh_flags |= ED_DISP_CONS | ED_DISP_STATUS | ED_DISP_SELECTION;
    redisplaySequences(xx, 1);

    return 0;
}


/*
 * Returns the list of hidden readings
 */
int *edGetHiddenReads(EdStruct *xx) {
    int *seqList, i, count = 0;

    if (xx->editorState == StateDown)
        return NULL;
  
    for (i = 1; i <= DBI_gelCount(xx); i++) {
        if (DB_Flags(xx, i) & DB_FLAG_INVIS)
	    count++;
    }

    if (NULL == (seqList = (int *)xmalloc(++count * sizeof(int))))
        return NULL;

    count = 0;
    for (i = 1; i <= DBI_gelCount(xx); i++) {
        if (DB_Flags(xx, i) & DB_FLAG_INVIS)
	    seqList[count++] = DB_Number(xx, i);
    }
    seqList[count] = 0;

    return seqList;
}

/*
 * Sets the trace for comparing against.
 */
void edSetTraceComparator(EdStruct *xx, int seq) {
    if (xx->compare_trace != -1 &&
	(DB_Flags(xx, xx->compare_trace) & DB_FLAG_SELECTED)) {
	edSelectRead(xx, xx->compare_trace, -1);
    }

    xx->compare_trace = seq;
    if (seq != -1 && !(DB_Flags(xx, seq) & DB_FLAG_SELECTED))
	edSelectRead(xx, seq, -1);
}


/* ----- 'brief' line manipulation ----- */

static void add_number(char *buf, int *j, int l1, int l2, int val) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*d", l1, l2, val);
	else
	    *j += sprintf(buf + *j, "%*d", l1, val);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*d", l2, val);
	else
	    *j += sprintf(buf + *j, "%d", val);
}

static void add_double(char *buf, int *j, int l1, int l2, double val) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*f", l1, l2, val);
	else
	    *j += sprintf(buf + *j, "%*f", l1, val);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*f", l2, val);
	else
	    *j += sprintf(buf + *j, "%f", val);
}

static void add_string(char *buf, int *j, int l1, int l2, char *str) {
    if (l1)
	if (l2)
	    *j += sprintf(buf + *j, "%*.*s", l1, l2, str);
	else
	    *j += sprintf(buf + *j, "%*s", l1, str);
    else
	if (l2)
	    *j += sprintf(buf + *j, "%.*s", l2, str);
	else
	    *j += sprintf(buf + *j, "%s", str);
}

/*
 * Formats tag information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 *
 * %%	Single % sign
 * %p	Tag position
 * %d	Tag direction (+/-/=, Raw: 0/1/2)
 * %D	Tag direction (----->/<-----/<---->, Raw: 0/1/2)
 * %t	Tag type (always 4 characters)
 * %l	Tag length
 * %#	Tag number (0 if unknown)
 * %c	Tag comment
 *
 * Additionally, some formats (p, l, n and c) can be specified as
 * %<number><format> (eg %80c) to allow AT MOST that many characters.
 */
int edSetBriefTag(EdStruct *xx, int seq, tagStruct *t, char *format) {
    char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    
    for (i = j = 0; format[i]; i++) {
	if (format[i] != '%') {
	    status_buf[j++] = format[i];
	    continue;
	}

	l1 = strtol(&format[++i], &cp, 10);
	i = cp - format;
	if (format[i] == '.') {
	    l2 = strtol(&format[++i], &cp, 10);
	    i = cp - format;
	} else {
	    l2 = 0;
	}
	if (format[i] == 'R') {
	    raw = 1;
	    i++;
	} else {
	    raw = 0;
	}

	switch(format[i]) {
	case '%':
	    status_buf[j++] = '%';
	    break;

	case 'd': /* Direction */
	case 'D':
	    {
		int sense = normaliseSense(xx, seq, t->tagrec.sense);
		char *dirs[] = {"----->", "<-----", "<---->"};

		if (raw)
		    add_number(status_buf, &j, l1, l2, sense);
		else
		    if (format[i] == 'd')
			j += sprintf(&status_buf[j], "%c", "+-="[sense]);
		    else
			add_string(status_buf, &j, l1, l2, dirs[sense]);
	    }
	    break;

	case 't': /* Type */
	    status_buf[j++] = t->tagrec.type.c[0];
	    status_buf[j++] = t->tagrec.type.c[1];
	    status_buf[j++] = t->tagrec.type.c[2];
	    status_buf[j++] = t->tagrec.type.c[3];
	    break;

	case 'p': /* Position */
	    add_number(status_buf, &j, l1, l2, t->tagrec.position);
	    break;

	case 'l': /* Length */
	    add_number(status_buf, &j, l1, l2, t->tagrec.length);
	    break;

	case '#': /* Number */
	    add_number(status_buf, &j, l1, l2, t->original_tag_id);
	    break;


	case 'c': /* Comment */
	    force_comment(DBI_io(xx), t);
	    add_string(status_buf, &j, l1, l2,
		       t->newcomment ? t->newcomment : "(no comment)");
	    break;

	default:
	    status_buf[j++] = format[i];
	}
    }
    status_buf[j] = 0;

    return tk_update_brief_line(xx, status_buf);
}

/*
 * Formats reading information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 *
 * %%	Single % sign
 * %n	Reading name (Raw: number)
 * %#	Reading number
 * %t	Trace name
 * %p	Position
 * %l	Clipped length
 * %L	Total length
 * %s	Start of clip
 * %e	End of clip
 * %S   Sense (whether complemented, +/-, Raw: 0/1)
 * %a	Chemistry (primer/terminator, Raw: integer)
 * %d	Strand (+/-, Raw 0/1)
 * %c	Reading freetext comment
 * %P	Primer (unknown/forward universal/reverse universal/forward custom/
 *              reverse custom,  Raw: 0/1/2/3/4)
 * %t   Trace name
 * %Tn	Template name (Raw: template number)
 * %T#	Template number
 * %Tv	Template vector (Raw: template vector number)
 * %Tc	Template consistency (Raw: as a number)
 * %Ti	Template insert size
 * %Cn	Clone name (Raw: clone number)
 * %C#	Clone number
 * %Cv	Clone vector (Raw: clone vector number)
 *
 * Additionally specifying %<number><format> forces AT MOST that many
 * characters to be displayed.
 * Specifying %R<format> (or %<number>R<format>) requests the raw data to
 * be displayed. This only works for some fields. Eg %Rp displays 0 to 4, but
 * %p displays, for instance, "forward universal"
 */
static int edSetBriefRead(EdStruct *xx, int seq, char *format) {
    char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    int rnum = DB_Number(xx, seq);
    GapIO *io = DBI_io(xx);
    
    for (i = j = 0; format[i]; i++) {
	if (format[i] != '%') {
	    status_buf[j++] = format[i];
	    continue;
	}

	l1 = strtol(&format[++i], &cp, 10);
	i = cp - format;
	if (format[i] == '.') {
	    l2 = strtol(&format[++i], &cp, 10);
	    i = cp - format;
	} else {
	    l2 = 0;
	}
	if (format[i] == 'R') {
	    raw = 1;
	    i++;
	} else {
	    raw = 0;
	}

	switch(format[i]) {
	case '%':
	    status_buf[j++] = '%';
	    break;

	case '#':
	    add_number(status_buf, &j, l1, l2, DB_Number(xx, seq));
	    break;

	case 'n':
	    if (raw)
		add_number(status_buf, &j, l1, l2, DB_Number(xx, seq));
	    else
		add_string(status_buf, &j, l1, l2, io_rname(io, rnum));
	    break;

	case 'p':
	    add_number(status_buf, &j, l1, l2, DB_RelPos(xx, seq));
	    break;

	case 'l':
	    add_number(status_buf, &j, l1, l2, DB_Length(xx, seq));
	    break;

	case 'L':
	    add_number(status_buf, &j, l1, l2, DB_Length2(xx, seq));
	    break;

	case 's':
	    add_number(status_buf, &j, l1, l2, DB_Start(xx, seq));
	    break;

	case 'e':
	    add_number(status_buf, &j, l1, l2, DB_End(xx, seq));
	    break;

	case 'S':
	    if (raw)
		add_number(status_buf, &j, l1, l2,
			   DB_Comp(xx, seq) == COMPLEMENTED ? 1 : 0);
	    else
		add_string(status_buf, &j, l1, l2,
			   DB_Comp(xx, seq) == COMPLEMENTED ? "+" : "-");
	    break;

	case 'a':
	    {
		GReadings r;
		gel_read(io, rnum, r);

		if (raw)
		    add_number(status_buf, &j, l1, l2, r.chemistry);
		else {
		    char *type_name;

		    switch (r.chemistry & GAP_CHEM_TYPE_MASK) {
		    case GAP_CHEM_TYPE_UNKNOWN:
			type_name = "unknown";
			break;
			
		    case GAP_CHEM_TYPE_ABI_RHOD:
			type_name = "rhodamine";
			break;
			
		    case GAP_CHEM_TYPE_ABI_DRHOD:
			type_name = "dRhodamine";
			break;
			
		    case GAP_CHEM_TYPE_BIGDYE1:
			type_name = "BigDyeV1";
			break;

		    case GAP_CHEM_TYPE_BIGDYE2:
			type_name = "BigDyeV2";
			break;
			
		    case GAP_CHEM_TYPE_BIGDYE3:
			type_name = "BigDyeV3";
			break;
			
		    case GAP_CHEM_TYPE_ET:
			type_name = "ET";
			break;
			
		    case GAP_CHEM_TYPE_MB_ET:
			type_name = "MegaBACE-ET";
			break;
			
		    case GAP_CHEM_TYPE_LICOR:
			type_name = "LiCor";
			break;

		    case GAP_CHEM_TYPE_SOLEXA:
			type_name = "Solexa";
			break;

		    case GAP_CHEM_TYPE_SOLID:
			type_name = "SOLiD";
			break;

		    case GAP_CHEM_TYPE_454:
			type_name = "454";
			break;

		    case GAP_CHEM_TYPE_OX_NANO:
			type_name = "Nanopore";
			break;


		    default:
			type_name = "undefined";
			break;
		    }
		    add_string(status_buf, &j, l1, l2, type_name);

		    if ((r.chemistry & GAP_CHEM_TYPE_MASK) 
			< GAP_CHEM_TYPE_SOLEXA) {
			add_string(status_buf, &j, l1, l2,
				   (r.chemistry & GAP_CHEM_TERMINATOR)
				   ? " terminator" : " primer");
		    }
		}
	    }
	    break;

	case 't':
	    {
		char path[FILE_NAME_LENGTH];
		char type[5];

		io_read_rd(io, rnum, path, FILE_NAME_LENGTH, type, 5);
		path[FILE_NAME_LENGTH-1] = 0;
		type[4] = 0;
		
		add_string(status_buf, &j, l1, l2, path);
	    }
	    break;

	case 'P':
	    {
		GReadings r;
		int primer;

		gel_read(io, rnum, r);
		primer = PRIMER_TYPE(r);

		if (raw)
		    add_number(status_buf, &j, l1, l2, primer);
		else {
		    char *str;
		    if      (primer == 1) str = "forward universal";
		    else if (primer == 2) str = "reverse universal";
		    else if (primer == 3) str = "forward custom";
		    else if (primer == 4) str = "reverse custom";
		    else                  str = "unknown";
		    add_string(status_buf, &j, l1, l2, str);
		}
	    }
	    break;

	case 'd':
	    {
		GReadings r;
		int strand;

		gel_read(io, rnum, r);
		strand = STRAND(r);

		if (raw)
		    add_number(status_buf, &j, l1, l2, strand);
		else {
		    char *str;
		    if      (strand == 0) str = "+";
		    else if (strand == 1) str = "-";
		    else                  str = "?";
		    add_string(status_buf, &j, l1, l2, str);
		}
	    }
	    break;

	case 'c': /* INFO note */
	    {
		GNotes n;
		GReadings r;
		int nnote;
		int ref_type = str2type("INFO");

		gel_read(io, rnum, r);
		for (nnote = r.notes; nnote; nnote = n.next) {
		    note_read(io, nnote, n);
		    if (n.type == ref_type)
			break;
		}
		if (nnote && n.annotation) {
		    char *text;
		    if (text = TextAllocRead(io, n.annotation)) {
			add_string(status_buf, &j, l1, l2, text);
			xfree(text);			
		    }
		}
	    }
	    break;

	case 'T': /* Template - subfields */
	    {
		GTemplates t;
		GReadings r;
		char name[DB_NAMELEN+1];
		template_c *tc;

		gel_read(io, rnum, r);
		if (0 == r.template && format[i+1] != '#') {
		    i++;
		    break;
		} else {
		    template_read(io, r.template, t);
		}

		if (r.template && DBI(xx)->templates[r.template]) {
		    tc = DBI(xx)->templates[r.template];
		} else {
		    tc = NULL;
		}

		switch (format[++i]) {
		case '#':
		    add_number(status_buf, &j, l1, l2, r.template);
		    break;

		case 'n':
		    if (raw) {
			add_number(status_buf, &j, l1, l2, r.template);
			break;
		    }

		    if (t.name) {
			TextRead(io, t.name, name, DB_NAMELEN);
			name[DB_NAMELEN] = 0;
		    } else {
			strcpy(name, "(unknown)");
		    }
		    add_string(status_buf, &j, l1, l2, name);
		    break;

		case 'i':
		    {
			char buf[1024];
			sprintf(buf, "%d..%d",
				t.insert_length_min, t.insert_length_max);
			add_string(status_buf, &j, l1, l2, buf);
		    }
		    break;

		case 'c':
		    {
			int tstat, guessed, spanned;
			char buf[10];

			tstat = tc ? tc->consistency : TEMP_CONSIST_UNKNOWN;
			guessed = tc->flags & (TEMP_FLAG_GUESSED_START |
					       TEMP_FLAG_GUESSED_END);
			spanned = tc->flags & TEMP_FLAG_SPANNING;

			if (raw) {
			    add_number(status_buf, &j, l1, l2, tstat);
			    break;
			}
			*buf = 0;

			if (tstat & TEMP_CONSIST_DIST) {
			    if ((tc->start > tc->end && tc->direction == 0) ||
				(tc->start < tc->end && tc->direction == 1))
				strcat(buf, "D");
			    else
				strcat(buf, "d");
			}
			if (tstat & TEMP_CONSIST_PRIMER)
			    strcat(buf, "P");
			if (tstat & TEMP_CONSIST_STRAND)
			    strcat(buf, "S");
			if (tstat & TEMP_CONSIST_UNKNOWN)
			    strcat(buf, "?");
			if (guessed)
			    strcat(buf, "E");
			if (tstat & TEMP_CONSIST_INTERDIST)
				strcat(buf, "I");
			if (spanned)
			    strcat(buf, "O");
			if (!tstat && !guessed)
			    strcat(buf, "ok");
			add_string(status_buf, &j, l1, l2, buf);
		    }
		    break;

		case 'v':
		    {
			GVectors v;

			if (raw) {
			    add_number(status_buf, &j, l1, l2, 0);
			    break;
			}

			if (t.vector) {
			    vector_read(io, t.vector, v);
			} else {
			    add_string(status_buf, &j, l1, l2, "(unknown)");
			}

			if (v.name) {
			    TextRead(io, v.name, name, DB_NAMELEN);
			    name[DB_NAMELEN] = 0;
			} else {
			    strcpy(name, "(unknown)");
			}
			add_string(status_buf, &j, l1, l2, name);
		    }
		    break;
		}
	    }
	    break;

	case 'C': /* Clone - subfields */
	    {
		GTemplates t;
		GClones c;
		GReadings r;
		char name[DB_NAMELEN+1];

		gel_read(io, rnum, r);
		if (r.template == 0) {
		    i++;
		    break;
		}
		template_read(io, r.template, t);
		if (0 == t.clone && format[i+1] != '#') {
		    i++;
		    break;
		} else {
		    clone_read(io, t.clone, c);
		}

		switch (format[++i]) {
		case '#':
		    add_number(status_buf, &j, l1, l2, t.clone);
		    break;

		case 'n':
		    if (raw) {
			add_number(status_buf, &j, l1, l2, t.clone);
			break;
		    }

		    if (c.name) {
			TextRead(io, c.name, name, DB_NAMELEN);
			name[DB_NAMELEN] = 0;
		    } else {
			strcpy(name, "(unknown)");
		    }
		    add_string(status_buf, &j, l1, l2, name);
		    break;

		case 'v':
		    {
			GVectors v;

			if (raw) {
			    add_number(status_buf, &j, l1, l2, c.vector);
			    break;
			}

			if (c.vector) {
			    vector_read(io, c.vector, v);
			} else {
			    add_string(status_buf, &j, l1, l2, "(unknown)");
			}

			if (v.name) {
			    TextRead(io, v.name, name, DB_NAMELEN);
			    name[DB_NAMELEN] = 0;
			} else {
			    strcpy(name, "(unknown)");
			}
			add_string(status_buf, &j, l1, l2, name);
		    }
		    break;
		}
	    }
	    break;

	default:
	    status_buf[j++] = format[i];
	}
    }
    status_buf[j] = 0;

    return tk_update_brief_line(xx, status_buf);
}

/*
 * Formats contig information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 *
 * %%	Single % sign
 * %n	Left most reading name (Raw: number)
 * %s	(As %n)
 * %e	Right most reading name (Raw: number)
 * %#   Contig number
 * %l	Length
 * %c	Contig freetext comment
 * %E   Expected number of errors
 */
static int edSetBriefContig(EdStruct *xx, int seq, char *format) {
    char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    GapIO *io = DBI_io(xx);
    int recompute = 0;
    int lc;
    
    for (i = j = 0; format[i]; i++) {

	if (format[i] != '%') {
	    status_buf[j++] = format[i];
	    continue;
	}

	l1 = strtol(&format[++i], &cp, 10);
	i = cp - format;
	if (format[i] == '.') {
	    l2 = strtol(&format[++i], &cp, 10);
	    i = cp - format;
	} else {
	    l2 = 0;
	}
	if (format[i] == 'R') {
	    raw = 1;
	    i++;
	} else {
	    raw = 0;
	}

	switch(format[i]) {
	case '%':
	    status_buf[j++] = '%';
	    break;

	case 'n':
	case 's':
	    {
		GContigs c;

		contig_read(io, -DB_Number(xx, seq), c);
		if (raw)
		    add_number(status_buf, &j, l1, l2, c.left);
		else
		    add_string(status_buf, &j, l1, l2, io_rname(io, c.left));
	    }
	    break;

	case 'e':
	    {
		GContigs c;

		contig_read(io, -DB_Number(xx, seq), c);
		if (raw)
		    add_number(status_buf, &j, l1, l2, c.right);
		else
		    add_string(status_buf, &j, l1, l2, io_rname(io, c.right));
	    }
	    break;

	case '#':
	    add_number(status_buf, &j, l1, l2, -DB_Number(xx, seq));
	    break;

	case 'l':
	case 'L':
	    add_number(status_buf, &j, l1, l2, DB_Length(xx, seq));
	    break;

	case 'c': /* INFO note */
	    {
		GNotes n;
		GContigs c;
		int nnote;
		int ref_type = str2type("INFO");

		contig_read(io, -DB_Number(xx, seq), c);
		nnote = c.notes;
		for (nnote = c.notes; nnote; nnote = n.next) {
		    note_read(io, nnote, n);
		    if (n.type == ref_type)
			break;
		}
		if (nnote && n.annotation) {
		    char *text;
		    if (text = TextAllocRead(io, n.annotation)) {
			add_string(status_buf, &j, l1, l2, text);
			xfree(text);			
		    }
		}
	    }
	    break;

	case 'E':
	    {
		char *con;
		float *qual;
		int i;
		double nerrs;

		if (consensus_mode != CONSENSUS_MODE_CONFIDENCE) {
		    add_string(status_buf, &j, l1, l2, "unknown");
		    break;
		}

		con = (char *)xmalloc(DB_Length(xx, seq));
		qual = (float *)xmalloc(DB_Length(xx, seq) * sizeof(float));
		if (NULL == con || NULL == qual) {
		    add_string(status_buf, &j, l1, l2, "unknown");
		    break;
		}

		calc_consensus(0, 1, DB_Length(xx, seq), CON_SUM,
			       con, NULL, qual, NULL,
			       xx->con_cut, xx->qual_cut,
			       contEd_info, (void *)xx);

		nerrs = 0.0;
		for (i = 0; i < DB_Length(xx, seq); i++) {
		    nerrs += pow(10.0, -qual[i]/10.0);
		}

		add_double(status_buf, &j, l1, l2, nerrs);

		xfree(con);
		xfree(qual);
		recompute = 1;
	    }

	default:
	    status_buf[j++] = format[i];
	}
    }
    status_buf[j] = 0;

    lc = tk_update_brief_line(xx, status_buf);
    return recompute ? -1 : lc;
}

/*
 * Set the brief status line with whatever is relevant for the data underneath
 * the current cursor position.
 */
int edSetBriefSeqStatus(EdStruct *xx, int x, int y) {
    int seq, pos;
    tagStruct *t;
    static int last_seq = -1, last_pos = -1, last_code = -1;

    /* Set 'seq' and 'pos' in the same manner as the edSetCursor function */
    if (-1 == (seq = edGetGelNumber(xx, x, y)))
	return -1;
    pos = xx->displayPos - DB_RelPos(xx, seq) + x + 1;

    if (xx->reveal_cutoffs) {
	if (pos < -DB_Start(xx, seq) + 1)
	    pos = -DB_Start(xx, seq) + 1;
	else
	    if (pos + DB_Start(xx, seq) > DB_Length2(xx, seq))
		pos = DB_Length2(xx, seq) - DB_Start(xx, seq) + 1;
    } else {
	if (pos < 1)
	    pos = 1;
	else
	    if (pos > DB_Length(xx, seq) + 1)
		pos = DB_Length(xx, seq) + 1;
    }

    /*
     * Is this the same as last time? If so only redraw if the brief line has
     * changed for some other reason.
     */
    if (seq == last_seq && pos == last_pos &&
	tk_update_brief_line(xx, NULL) == last_code)
	return 0;

    /* Is there a tag at this point? */
    if (NULL == (t = findTag(xx, seq, pos + DB_Start(xx, seq))))
	return 0;

    force_comment(DBI_io(xx), t);

    last_code = edSetBriefTag(xx, seq, t,
			      get_default_string(EDINTERP(xx->ed), gap_defs,
						 "TAG_BRIEF_FORMAT"));

    last_seq = seq;
    last_pos = pos;

    return 0;
}

/*
 * Set the brief status line with whatever is relevant for the data underneath
 * the current cursor position.
 */
int edSetBriefNameStatus(EdStruct *xx, int x, int y) {
    int seq;
    static int last_seq = -1, last_code = -1;

    /* Set 'seq' and 'pos' in the same manner as the edSetCursor function */
    if (-1 == (seq = edGetGelNumber(xx, x, y)))
	return -1;

    /* Optimisation - no change to status line required */
    if (seq == last_seq && tk_update_brief_line(xx, NULL) == last_code)
	return 0;

    if (seq) {
	last_code = edSetBriefRead(xx, seq,
				   get_default_string(EDINTERP(xx->ed),
				       gap_defs, "READ_BRIEF_FORMAT"));
    } else {
	last_code = edSetBriefContig(xx, seq,
				   get_default_string(EDINTERP(xx->ed),
				       gap_defs, "CONTIG_BRIEF_FORMAT"));
    }

    last_seq = seq;

    return 0;
}

/*
 * Counts the number of pads between the start of the contig and 'pos'
 * This is used to provide an unpadded base position.
 */
static int count_pads(EdStruct *xx, int pos) {
    char *con;
    int i, pads;

    /* We cannot have any pads in the cutoff data */
    if (pos > DB_Length(xx, 0))
	pos = DB_Length(xx, 0);

    if (pos <= 0 || NULL == (con = xmalloc(DB_Length(xx, 0)+1)))
	return 0;

    /* Ask for the whole consensus - this way it will be cached */
    DBcalcConsensus(xx, 1, DB_Length(xx, 0), con, NULL, BOTH_STRANDS);

    for (pads = i = 0; i < pos; i++)
	if (con[i] == '*')
	    pads++;

    xfree(con);

    return pads;
}


/*
 * Calculates the unpadded base number for position 'pos'.
 * This function is typically called many times in quick succession with pos,
 * pos+1, pos+2 etc.
 *
 * First call this function with cache_count set to a positive number to
 * indicate how many positions to calculate.
 *
 * Next call this function several times (in theory, 'cache_count' times).
 * cache_count should be set to 0.
 *
 * Finally call again with cache_count set to -1 to free up the memory.
 *
 * Returns the base position.
 */
int edUnpaddedBaseNumber(EdStruct *xx, int pos, int cache_count) {
    static int last_pos = 0, last_pads;
    static char *con = NULL;
    int pads, i;

    if (cache_count > 0) {
	if (pos + cache_count >= 0) {
	    /* Ask for the whole consensus - this way it will be cached */
	    if (NULL == (con = xmalloc(DB_Length(xx, 0)+1)))
		return 0;
	    DBcalcConsensus(xx, 1, DB_Length(xx, 0), con, NULL,
			    BOTH_STRANDS);
	}

	pads = 0;
	if (pos >= 0) {
	    for (i = 0; i < pos-1 && i < DB_Length(xx, 0)-1; i++) {
		if (con[i] == '*')
		    pads++;
	    }
	}
	last_pos = pos-1;
	last_pads = pads;

    } else if (cache_count == -1) {
	if (con)
	    xfree(con);
	con = NULL;
	return 0;
    } else {
	if (pos == last_pos + 1) {
	    if (pos >= 0) {
		if (pos < DB_Length(xx, 0)-1 && con[pos-1] == '*')
		    last_pads++;
	    } else {
		last_pads = 0;
	    }
	    pads = last_pads;
	    last_pos++;
	} else {
	    /* Getting here implies incorrect usage of this function */
	    return 0;
	}
    }

    return pos - pads;
}

/*
 * Set the brief status line with information about an individual base.
 *
 * 'mode' is 1 or 2 for 1st and 2nd mouse buttons. This allows us to bind
 * different information to each mouse button. The 1st (left) button is the
 * default mode, and is used after editor searches.
 *
 * %%   Single % sign
 * %c   Confidence value (phred style)
 * %p   Confidence value (as probability)
 * %P   Padded consensus base position
 * %U   Unpadded consensus base position
 */
int edSetBriefSeqBase(EdStruct *xx, int x, int y, int mode) {
    char status_buf[8192]; /* NB: no bounds checking! */
    int i, j, l1, l2, raw;
    char *cp;
    int1 *conf;
    float c2[2] = {0.0, 0.0};
    int seq, pos;
    static char *format1 = NULL;
    static char *format2 = NULL;
    char *format;

    /* Do nothing is quality cutoff is -1 when left button is used. */
    /*
    if (mode == 1 && xx->qual_cut < 0)
	return 0;
    */

    /* Get the format */
    if (mode == 2) {
	if (!format2) {
	    format2 = get_default_string(EDINTERP(xx->ed), gap_defs,
					 "BASE_BRIEF_FORMAT2");
	}
	format = format2;
    } else {
	if (!format1) {
	    format1 = get_default_string(EDINTERP(xx->ed), gap_defs,
					 "BASE_BRIEF_FORMAT1");
	}
	format = format1;
    }
    
    /* Set 'seq' and 'pos' in the same manner as the edSetCursor function */
    if (x == -1 && y == -1) {
	seq = xx->cursorSeq;
	pos = xx->cursorPos;
    } else {
	if (-1 == (seq = edGetGelNumber(xx, x, y)))
	    return -1;
	pos = xx->displayPos - DB_RelPos(xx, seq) + x + 1;
    }

    if (xx->reveal_cutoffs) {
	if (pos < -DB_Start(xx, seq) + 1)
	    pos = -DB_Start(xx, seq) + 1;
	else
	    if (pos + DB_Start(xx, seq) > DB_Length2(xx, seq))
		pos = DB_Length2(xx, seq) - DB_Start(xx, seq) + 1;
    } else {
	if (pos < 1)
	    pos = 1;
	else
	    if (pos > DB_Length(xx, seq) + 1)
		pos = DB_Length(xx, seq) + 1;
    }
    
    (void) DBgetSeq(DBI(xx), seq);
    if (seq) {
	conf = DB_Conf(xx, seq) + DB_Start(xx, seq);
    } else {
	char dummy[2];
	conf = NULL;
	DBcalcConsensus(xx, pos, 1, dummy, c2, BOTH_STRANDS);
    }

    /* Update the format line */
    for (i = j = 0; format[i]; i++) {
	if (format[i] != '%') {
	    status_buf[j++] = format[i];
	    continue;
	}

	l1 = strtol(&format[++i], &cp, 10);
	i = cp - format;
	if (format[i] == '.') {
	    l2 = strtol(&format[++i], &cp, 10);
	    i = cp - format;
	} else {
	    l2 = 0;
	}
	if (format[i] == 'R') {
	    raw = 1;
	    i++;
	} else {
	    raw = 0;
	}

	switch(format[i]) {
	case '%':
	    status_buf[j++] = '%';
	    break;

	case 'c':
	    if (conf && pos <= DB_Length2(xx, seq) - DB_Start(xx, seq))
		add_number(status_buf, &j, l1, l2, conf[pos-1]);
	    else
		add_number(status_buf, &j, l1, l2, c2[0]+0.499);
	    /*
		add_string(status_buf, &j, l1, l2, "/");
		add_double(status_buf, &j, l1, l2, c2[0]);
	    */
	    break;

	case 'p':
	    if (xx->consensus_mode == CONSENSUS_MODE_CONFIDENCE) {
		double p;
		if (conf && pos <= DB_Length2(xx, seq) - DB_Start(xx, seq))
		    p = 1 - pow(10, conf[pos-1]/-10.0);
		else
		    p = 1 - pow(10, c2[0]/-10.0);
		add_double(status_buf, &j, l1, l2, p);
	    } else {
		add_string(status_buf, &j, l1, l2, "-");
	    }
	    break;

	case 'P':
	    {
		int pos, ref;

		pos = DB_RelPos(xx, xx->cursorSeq) + xx->cursorPos - 1;

		if (ref = DBI(xx)->reference_seq) {
		    char *bases = DBgetSeq(DBI(xx), DBI(xx)->reference_seq);
		    int npads, i;

		    pos -= DB_RelPos(xx, ref) - 1;

		    if (DB_Comp(xx, ref) == UNCOMPLEMENTED) {
			/* Standard orientation */
			for (npads = i = 0;
			     i < pos && i < DB_Length2(xx, ref);
			     i++) {
			    if (bases[i] == '*')
				npads++;
			}
			pos = pos - npads + DBI(xx)->reference_offset-1;
		    } else {
			/* Reverse orientation */
			for (npads = 0, i = DB_Length(xx, ref)-1;
			     i >= pos-1 && i >= 0;
			     i--) {
			    if (bases[i] == '*')
				npads++;
			}
			pos = DB_Length(xx, ref) - pos - npads +
			    DBI(xx)->reference_offset;
		    }
		    if (DBI(xx)->reference_len) {
			pos = ((pos - 1) % DBI(xx)->reference_len) + 1;
			while (pos < 1)
			    pos += DBI(xx)->reference_len;
		    }
		}
		add_number(status_buf, &j, l1, l2, pos);
	    }
	    break;

	case 'U':
	    {
		int pos;

		pos = DB_RelPos(xx, xx->cursorSeq) + xx->cursorPos - 1;

		if (DBI(xx)->reference_seq) {
		    char *bases = DBgetSeq(DBI(xx), DBI(xx)->reference_seq);
		    int npads, i;

		    for (npads = i = 0; i < pos &&
			 i < DB_Length2(xx, DBI(xx)->reference_seq); i++) {
			if (bases[i] == '*')
			    npads++;
		    }
		    pos = pos - npads + DBI(xx)->reference_offset-1;
		    if (DBI(xx)->reference_len) {
			pos = ((pos - 1) % DBI(xx)->reference_len) + 1;
			while (pos < 0)
			    pos += DBI(xx)->reference_len;
		    }
		} else {
		    pos -= count_pads(xx, pos);
		}
		add_number(status_buf, &j, l1, l2, pos);
	    }
	    break;

	default:
	    status_buf[j++] = format[i];
	}
    }
    status_buf[j] = 0;

    return tk_update_brief_line(xx, status_buf);
}

/*
 * List the confidence values in a tabular form. Based on the list_confidence
 * command in newgap_cmds.c.
 */
int edListConfidence(EdStruct *xx, int start, int end, int info_only) {
    char status_buf[8192];
    char *con;
    float *qual;
    double err, err_rate;
    int i;
    double total_errs;
    int freqs[101];

    for (i = 0; i <= 100; i++)
	freqs[i] = 0;

    qual = (float *)xmalloc((end - start + 1) * sizeof(*qual));
    con = (char *)xmalloc((end - start + 1) * sizeof(*con));
    if (!qual || !con)
	return -1;
    
    /* Get the consensus along with the confidence values */
    calc_consensus(0, start, end, CON_SUM, con, NULL, qual, NULL,
		   xx->con_cut, xx->qual_cut,
		   contEd_info, (void *)xx);

    for (i = 0; i < end - start + 1; i++) {
	if (qual[i] < 0) qual[i] = 0;
	if (qual[i] > 100) qual[i] = 100;
	freqs[(int)(qual[i] + 0.499)]++;
    }

    xfree(qual);
    xfree(con);

    if (!info_only) {
	total_errs = list_confidence(freqs, end - start + 1);
    }

    total_errs = 0.0;
    for (i = 0; i <= 99; i++) {
	err = freqs[i] * pow(10.0, -i/10.0);
	total_errs += err;
    }
    err_rate = total_errs ? (end - start + 1) / total_errs : 0;

    sprintf(status_buf, "Expected no. of errors between %d and %d is %.2f. "
	    "Error rate = 1/%.0f",
	    start, end, total_errs, err_rate);
    tk_update_brief_line(xx, status_buf);
    
    return 0;
}


#if 0
/*
 * Mark the reference sequence.
 * Returns 0 for success, -1 for failure.
 */
int edSetReferenceSeq(EdStruct *xx, int seq, int length, int offset) {
    GReadings r;
    GNotes n;
    int nnote;
    int ref_type = str2type("REFS");
    GapIO *io = DBI_io(xx);
    char comment[1024];

    /* Find existing REFS notes and delete them */
    if (DBI(xx)->reference_seq) {
	gel_read(io, DB_Number(xx, seq), r);
	for (nnote = r.notes; nnote; nnote = n.next) {
	    note_read(io, nnote, n);
	    if (n.type == ref_type) {
		delete_note(io, nnote);
	    }
	}
    }

    /* Update the cached editor info */
    DBI(xx)->reference_seq = seq;
    DBI(xx)->reference_len = length;
    DBI(xx)->reference_offset = offset;

    /* Create a new note on the new sequence */
    nnote = new_note(io, ref_type, GT_Readings, DB_Number(xx, seq));
    if (!nnote)
	return -1;

    /* Update the note comment */
    seq = DB_Number(xx, seq);
    sprintf(comment, "sequence %d %d", offset, length);
    if (0 != edit_note(io, nnote, NULL /* type */, comment))
	return -1;

    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 1);

    return 0;
}
#endif

/*
 * Turns on or off the mini traces. The height is how many extra lines
 * to reserve for the traces. Height of zero disables them.
 */
void edSetMiniTraces(EdStruct *xx, int height) {
    /* Shutdown any old ones first */
    if (xx->lines_per_seq != 1) {
	tman_shutdown_traces(xx, 1);
    }

    if (height != 0) {
	int *seqList;
	int i;

	seqList = sequencesInRegion(xx, xx->displayPos, xx->displayWidth);

	/* Do this now as the height is needed in showTrace */
	xx->lines_per_seq = height+1;

	for (i = 0; seqList[i]; i++) {
	    showTrace(xx, seqList[i],
		      xx->displayPos + xx->displayWidth/2 -
		          (DB_RelPos(xx, seqList[i])-1),
		      xx->fontWidth, 0, 1);
	}
    }

    /*
     * Change the lines per seq and redisplay
     * The low level sheet calls are to force a complete clear and redraw and
     * hence avoid the optimisation of only redrawing the sequence lines that
     * have changed (and not blanking the lines inbetween).
     */
    sheet_clear(&TKSHEET(xx->names)->sw);
    sheet_clear(&TKSHEET(xx->ed)->sw);
    TKSHEET(xx->names)->flags |= SHEET_REDRAW_ALL;
    TKSHEET(xx->ed)->flags |= SHEET_REDRAW_ALL;
    xx->lines_per_seq = height+1;
    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 0);
    SheetDisplay((ClientData)xx->names);
    SheetDisplay((ClientData)xx->ed);
}


/*
 * Finds the next difference between editors in the join editor and positions
 * the cursor there.
 */
void edNextDifference(EdStruct *xx, int fwd) {
    char *con1 = NULL, *con2 = NULL;
    int pos1, pos2;
    int clen1, clen2;

    if (!inJoinMode(xx))
	goto error;
    
    /* Get the consensi */
    clen1 = DB_Length(xx->link->xx[0], 0);
    clen2 = DB_Length(xx->link->xx[1], 0);
    if (!(con1 = xmalloc(clen1+1)))
	goto error;
    if (!(con2 = xmalloc(clen2+1)))
	goto error;
    DBcalcConsensus(xx->link->xx[0], 1, clen1, con1, NULL, BOTH_STRANDS);
    DBcalcConsensus(xx->link->xx[1], 1, clen2, con2, NULL, BOTH_STRANDS);

    /* Find the current cursor positions */
    pos1 = positionInContig(xx->link->xx[0],
			    xx->link->xx[0]->cursorSeq,
			    xx->link->xx[0]->cursorPos);
    pos2 = pos1 + xx->link->lockOffset;

    if (fwd) {
	while (++pos1 <= clen1 && ++pos2 <= clen2) {
	    if (con1[pos1-1] != con2[pos2-1])
		break;
	}
    } else {
	while (--pos1 > 0 && --pos2 > 0) {
	    if (con1[pos1-1] != con2[pos2-1])
		break;
	}
    }

    setCursorPosSeq(xx->link->xx[0], pos1, 0);
    setCursorPosSeq(xx->link->xx[1], pos2, 0);
    redisplayWithCursor(xx->link->xx[0]);
    redisplayWithCursor(xx->link->xx[1]);

 error:
    if (con1)
	xfree(con1);
    if (con2)
	xfree(con2);
}

void edViewSet(EdStruct *xx, int set) {
    xx->curr_set = set;
    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 0);
}

void edMoveSet(EdStruct *xx, int set_num, int nseqs, char **seqs) {
    int i;

    if (!xx->set) {
	xx->set = (int *)xcalloc(DBI(xx)->DB_gelCount+1, sizeof(int));
    }

    for (i = 0; i < nseqs; i++) {
	int rnum = get_gel_num(DBI_io(xx), seqs[i], GGN_ID);
	if (rnum > 0)
	    rnum = rnum_to_edseq(xx, rnum);
	if (rnum > 0) {
	    xx->set[rnum] = set_num;
	}
    }

    /* Incase this expands (or shrinks) our set count, we recalloc */
    if (set_num > xx->nsets) {
	xx->set_collapsed = (int *)xrealloc(xx->set_collapsed,
					    (set_num+1)*sizeof(int));
	for (i = xx->nsets+1; i <= set_num; i++)
	    xx->set_collapsed[i] = 0;
	xx->nsets = set_num;
    }

    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 0);
}

/*
 * Set the set collapsed status.
 * 0  => not collapsed
 * 1  => collapsed
 * -1 => toggle
 *
 * Returns the new status
 */
int edCollapseSet(EdStruct *xx, int set, int mode) {
    switch (mode) {
    case 0:
	xx->set_collapsed[set] = 0;
	break;
    case 1:
	xx->set_collapsed[set] = 1;
	break;
    case -1:
	xx->set_collapsed[set] ^= 1;
    }
    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 0);

    return xx->set_collapsed[set];
}

/*
 * Returns 0 if seq is not in a set.
 * Returns <0 for the set number if seq is the "topmost" visible in that
 * set.
 * Returns >0 for the set number if seq is in a set, but isn't the topmost.
 */
int edFindSet(EdStruct *xx, int seq) {
    int *seqList;
    int pos = xx->displayPos;
    int width = xx->displayWidth;
    int lastset = 0, set;
    int i, k;

    /* Find out what's visible */
    seqList = sequencesOnScreen(xx,pos, width);

    for (i = xx->lines_per_seq-1; i < xx->displayHeight + xx->lines_per_seq-1;
	 i+=xx->lines_per_seq) {
	if (i < xx->displayHeight-1)
	    k = i/xx->lines_per_seq + xx->displayYPos;
	else
	    k = (xx->totalHeight-1)/xx->lines_per_seq;


	set = xx->set ? xx->set[seqList[k]] : 0;

	if (set != lastset) {
	    if (seqList[k] == seq)
		return -set;
	}
	if (seqList[k] == seq)
	    return set;
	lastset = set;
    }

    return 0;
}
