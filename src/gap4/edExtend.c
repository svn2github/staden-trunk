/*
 * File: extend.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: For moving cutoffs
 *
 * Created: ?
 * Updated: 16 April 1993
 *
 */

/*
 * Change log:
 *
 *   1/10/91 SD  Added calculateConsensusLength to extend, unextend, undo_unextend
 *   29/4/92 SD  Changes related to general speed up in edUtils.c
 *   18/5/92 SD  Construct (*nc=*++nc) not liked by dec
 */

#include <stdio.h>
#include <stdlib.h>

#include "edUtils.h"
#include "fortran.h"
#include "tagUtils.h"
#include "contigEditor.h"
#include "undo.h"
#include "select.h"
#include "extend.h"

#define OTHER_END(E) ((E)==EXTEND_LEFT ? EXTEND_RIGHT : EXTEND_LEFT)

int adjustMark(EdStruct *xx, int seq, int num_bases, int direction, int end);

/*
 * Handle cut-off adjust
 */
int meta_arrow (EdStruct *xx, int key)
{
    
    int seq,pos;
    int seq_length;
    int end;
    
    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	bell();
	return 1;
    }
    
    seq = xx->cursorSeq;
    pos = xx->cursorPos;
    seq_length = DB_Length(xx,seq);
    
    /*
     * determine which operation is to take place
     */
    if (!seq) {
	/* no operation on consensus */
	if (key==EXTEND_LEFT)
	    edCursorLeft(xx);
	else
	    edCursorRight(xx);
    } else {
	/*
	 * determine which end we are operating on
	 */
	if (seq_length == 0)
	    /*
	     * a special case to make things work nicely
	     */
	    end = key;
	else if (pos == 1)
	    end = EXTEND_LEFT;
	else if (pos == seq_length+1) 
	    end = EXTEND_RIGHT;
	else
	    end = 0;

	/*
	 * May need to redraw numbering even if consensus length doesn't
	 * change
	 */
	if (seq == DBI(xx)->reference_seq)
	    xx->refresh_flags |= ED_DISP_RULER;

	if (end) {
	    if (! adjustMark(xx,seq,1,key,end)) {
		redisplayWithCursor(xx);
	    } else
		bell();
	} else {
	    if (key==EXTEND_LEFT)
		edCursorLeft(xx);
	    else
		edCursorRight(xx);
	}
	
    }

    getExtents(xx);
    return 0;
}




/*
 * Zap (unextend) to right end
 */
void zap_Right (EdStruct *xx)
{
    int seq,pos;
    int seq_length;
    int num_bases;
    
    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE) || xx->cursorSeq == 0) {
	bell();
	return;
    }
    
    seq = xx->cursorSeq;
    pos = xx->cursorPos;
    seq_length = DB_Length(xx,seq);

    /* Disallow in left cutoff region */
    if (xx->cursorPos <= 0) {
	bell();
	return;
    }

    num_bases = seq_length - pos + 1;
    setCursorPos(xx, seq_length+1);

    /*
     * May need to redraw numbering even if consensus length doesn't
     * change
     */
    if (seq == DBI(xx)->reference_seq)
	xx->refresh_flags |= ED_DISP_RULER;
    
    if (num_bases>0) {
	if (! adjustMark(xx,seq,num_bases,EXTEND_LEFT,EXTEND_RIGHT)) {
	    redisplayWithCursor(xx);
	} else
	    bell();
    } else {
	if (! adjustMark(xx,seq,-num_bases,EXTEND_RIGHT,EXTEND_RIGHT)) {
	    redisplayWithCursor(xx);
	} else
	    bell();
    }

    getExtents(xx);
}




/*
 * Zap (unextend) to left end
 */
void zap_Left (EdStruct *xx)
{
    int seq,pos;
    int seq_length;
    int num_bases;
    
    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE) || xx->cursorSeq == 0) {
	bell();
	return;
    }
    
    seq = xx->cursorSeq;
    pos = xx->cursorPos;
    seq_length = DB_Length(xx,seq);

    /* Disallow in right cutoff region */
    if (xx->cursorPos > seq_length+1) {
	bell();
	return;
    }

    num_bases = pos - 1;
    setCursorPos(xx, 1);

    if (num_bases>0) {
	if (! adjustMark(xx,seq,num_bases,EXTEND_RIGHT,EXTEND_LEFT)) {
	    redisplayWithCursor(xx);
	} else
	    bell();
    } else {
	if (! adjustMark(xx,seq,-num_bases,EXTEND_LEFT,EXTEND_LEFT)) {
	    redisplayWithCursor(xx);
	} else
	    bell();
    }

    getExtents(xx);
}




/*************************************************************/
/*************************************************************/
/*************************************************************/
/*************************************************************/
/*************************************************************/
/*************************************************************/


/*************************************************************
 ** Low level
 *************************************************************/


int _adjust_ends(DBInfo *db, int seq, int start_bases, int end_bases,
		 int seq_flags)
{
    _DBsetStart(db,seq,_DB_Start(db,seq)+start_bases);
    _DBsetEnd(db,seq,_DB_End(db,seq)+end_bases);
    _DBsetLength(db,seq,_DB_Length(db,seq)-start_bases+end_bases);
    _DBsetFlags(db,seq,seq_flags);

    return 0;
}




/*************************************************************
 ** Undo level
 *************************************************************/

int U_adjust_ends(EdStruct *xx, int seq, int start_bases, int end_bases)
{
    UndoStruct *u;
    int seq_flags, old_seq_flags;
    
    old_seq_flags = DB_Flags(xx,seq);
    /* UNDO - _left_cutoff_add_bases(xx,seq,tag,num_bases,bases) */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoAdjustEnds;
	u->sequence = seq;
	u->info.adjust_ends.start_bases = -start_bases;
	u->info.adjust_ends.end_bases   = -end_bases;
	u->info.adjust_ends.seq_flags   = old_seq_flags;
	recordUndo(DBI(xx),u);
    }
    
    seq_flags = old_seq_flags | DB_FLAG_SEQ_MODIFIED | DB_FLAG_REL_MODIFIED;
    return _adjust_ends(DBI(xx), seq, start_bases, end_bases, seq_flags);
}


void U_adjust_display(EdStruct *xx, int n)
{
    UndoStruct *u;
    
    /* UNDO - xx->displayPos += n; */
    if ( (u = newUndoStruct(DBI(xx))) != NULL ) {
	u->db = DBI(xx);
	u->command = UndoAdjustDisplay;
	u->info.adjust_display.xx = xx;
	u->info.adjust_display.position = -n;
	u->info.adjust_display.editor_id = xx->editor_id;
	recordUndo(DBI(xx),u);
    }
    
    xx->displayPos += n;
    xx->refresh_flags |= ED_DISP_RULER | ED_DISP_SCROLL;
}






/*************************************************************
 ** Functional level
 *************************************************************/












/*
 * Extend a reading by num_bases in a particular direction
 */
int adjustMark(EdStruct *xx, int seq, int num_bases, int direction, int end)
{
    int max_bases;
    int extend;
    int orig_end;
    int comp;
    int old_len, new_len;
    
    /*
     * Can't adjust consensus!
     */
    if (seq == 0) return 1;
    
    old_len = DB_Length(xx, 0);

    /*
     * Determine if we are extending or unextending
     * Determine whether sequence is in the original sense
     * Determine original end we're dealing with
     */
    extend   = (end==direction);
    comp     = (DB_Comp(xx,seq)==COMPLEMENTED);
    orig_end = comp?OTHER_END(end):end;



    /*
     * Don't take more than is available
     */
    if (extend) {
	/* rules for extending */
	if (end == EXTEND_LEFT)
	    max_bases = DB_Start(xx,seq);
	else
	    max_bases = DB_Length2(xx,seq) - DB_End(xx,seq) + 1;
    } else {
	/* rules for unextending */
	max_bases = DB_Length(xx,seq) - 1;
    }
    if (num_bases > max_bases) num_bases = max_bases;

    /* no more room */
    if (num_bases <= 0) return 1;
    
    /* START UNDO */
    openUndo(DBI(xx));


    /*
     * Phase 1 - adjust marks
     */
    if (extend) {
	/* extending */
	if (end == EXTEND_LEFT) {
	    U_adjust_ends(xx,seq,-num_bases,0);
	} else {
	    U_adjust_ends(xx,seq,0,num_bases);
	}
	    
    } else {
	/* unextending */
	if (end == EXTEND_LEFT) {
	    U_adjust_ends(xx,seq,num_bases,0);
	} else {
	    U_adjust_ends(xx,seq,0,-num_bases);
	}
    }
    

    /*
     * Phase 1b - when dealing with left end we need to adjust position in gel
     */
    if (end == EXTEND_LEFT) {
	if (extend)
	    shiftLeft(xx,seq,num_bases);
	else
	    shiftRight(xx,seq,num_bases);
    }
    

    /*
     * Tags are now in absolute positions
     */
#if 0
    /*
     * Phase 2 - adjust tag positions
     */
    /* complemented? current_gel_end key_direction */
    if (extend) {
	/* determine current position of left end */
	if (direction==EXTEND_LEFT)
	    pos = 1;			/* fwd left  left / rev left  left */
	else
	    pos = DB_Length(xx,seq) + 1;/* rev right right/ fwd right right */
	tagInsertBases(xx,seq,pos,num_bases);
    } else {
	if (orig_end == EXTEND_LEFT) {
	    if (direction==EXTEND_LEFT)
		pos = DB_Length(xx,seq) + 1 - num_bases; /* rev right left */
	    else
		pos = 1;				 /* fwd left  right */
	} else {
	    if (direction==EXTEND_LEFT)
		pos = DB_Length(xx,seq) + 2 - num_bases; /* fwd right left */
	    else
		pos = 0;				 /* rev left  right */
	}
	tagDeleteBases(xx,seq,pos,num_bases);
    }
#endif


    /*
     * Phase 3 - adjust cursor position
     */
    if (end == EXTEND_RIGHT) {
	if (extend)
	    U_adjust_cursor(xx,num_bases);
	else
	    U_adjust_cursor(xx,-num_bases);
    } else {
	/* YUK! force cursor position on undo */
	U_adjust_cursor(xx,0);
    }

    
    /*
     * Phase 4 - adjust consensus length = Lengths are a changing
     */

    if (end == EXTEND_LEFT && DB_Length(xx, 0) != old_len) {
	int move_by = DB_Length(xx, 0) - old_len;
	U_adjust_display(xx, move_by);
    }


    if (xx->link) {
	xx->link->lockOffset = xx->link->xx[1]->displayPos -
	    xx->link->xx[0]->displayPos;
	setCursorPos(xx, xx->cursorPos);
    }

#if 0
    if (end == EXTEND_LEFT && DB_Length(xx, 0) != old_len) {
	if (DB_Length(xx, 0) > old_len) {
	    xx->displayPos += num_bases;
	    /* decDisplayPos(xx, num_bases); */
	} else {
	    xx->displayPos -= num_bases;
	    /* incDisplayPos(xx, num_bases); */
	}
    }
#endif

    new_len = calculate_consensus_length(xx);

    /*
     * This bit doesn't appear to be needed as it's apparently done elsewhere.
     * But I certainly can't see where! - jkb
     */
    if (DB_Length(xx, 0) != new_len) {
	U_change_consensus_length(xx, new_len);
    }


    /*
     * Phase 5 - adjust consensus tags
     */
    if (DB_Length(xx, 0) != old_len) {
	if (end == EXTEND_LEFT) {
	    if (DB_Length(xx, 0) > old_len)
		tagInsertBases(xx, 0, 1, DB_Length(xx, 0) - old_len);
	    else
		tagDeleteBases(xx, 0, old_len - DB_Length(xx, 0),
			       old_len - DB_Length(xx, 0));
	} else {
	    if (DB_Length(xx, 0) < old_len)
		tagDeleteBases(xx, 0, old_len, old_len - DB_Length(xx, 0));
	}
    }
    
#if 0
    if (DB_Length(xx, 0) != old_len) {
	if (extend && end == EXTEND_LEFT) {
		tagInsertBases(xx, 0, 0, DB_Length(xx, 0) - old_len);
	} else if (end == EXTEND_LEFT) {
	    tagDeleteBases(xx, 0, 0,
			   old_len - DB_Length(xx, 0));
	}
    }
#endif

    closeUndo(xx, DBI(xx));
    invalidate_consensus(xx);

    return 0;
    
}
