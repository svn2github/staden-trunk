/*
 * File: undo.c
 * Version:
 *
 * Author: 
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "undo.h"
#include "edUtils.h"
#include "tagUtils.h"
#include "extend.h"
#include "contigEditor.h"
#include "xalloc.h"
#include "text_output.h"


/* static void DumpUndo(UndoStruct *u); */

/*************************************************************
 ** Packs String, int1 and int2 into a single buffer.
 *************************************************************/

/*
 * Copy Base, Conf and Opos arrays into a packed structure.
 */
void packBCO(PackedBCO *ps, char *bases, int1 *conf,
	     int2 *opos, int len) {
    size_t size;

    size = 4 * len; /* len * (int2 + int1 + int1) */
    ps->len = len;

    if (size > sizeof(void *)) {
	if (!(ps->data.ptr = (void *)xmalloc(size)))
	    return;

	memcpy(ps->data.ptr, opos, len * 2);
	memcpy(ps->data.ptr + 2*len, bases, len);
	memcpy(ps->data.ptr + 3*len, conf, len);
    } else {
	memcpy(ps->data.array, opos, len * 2);
	memcpy(ps->data.array + 2*len, bases, len);
	memcpy(ps->data.array + 3*len, conf, len);
    }
}

/*
 * Recover space
 */
void unpackBCO(PackedBCO *ps) {
    if ((size_t)(4 * ps->len) > sizeof(char *))
	xfree(ps->data.ptr);
    ps->len = 0;
}

/*
 * Return pointers to the packed strings
 */
void getBCO(PackedBCO *ps, char **bases, int1 **conf, int2 **opos) {
    if ((size_t)(4 * ps->len) > sizeof(char *)) {
	*opos  = (int2 *)ps->data.ptr;
	*bases = ps->data.ptr + 2 * ps->len;
	*conf  = (int1 *)(ps->data.ptr + 3 * ps->len);
    } else {
	*opos  = (int2 *)ps->data.array;
	*bases = ps->data.array + 2 * ps->len;
	*conf  = (int1 *)(ps->data.array + 3 * ps->len);
    }
}





/*************************************************************
 *
 *************************************************************/

/*
 * Number of undos since the last auto save. Must be <= MAX_SAVE_EDITS
 */
#define AUTO_SAVE_EDITS 50

#define FREE_RECLAIM 1
#define FREE_UNDO    2


/*
 *    0 - no edits made
 *    1 - edits made
 */
int _editsMade(DBInfo *db)
{
    return (_DBI_edits_made(db) || _DBI_since_auto_save(db)) &&
	(_DBI_num_undo(db) || _DBI_discarded_undos(db));
}

/*
 *    0 - no edits made
 *    1 - edits made
 */
int editsMade(EdStruct *xx)
{
    return _editsMade(DBI(xx));
}


UndoStruct *newUndoStruct(DBInfo *db)
/*
 * Allocate a new undo struct
 */
{
    UndoStruct *u;

    if (!_DBI_store_undo(db))
	return NULL;

    u = (UndoStruct *)xmalloc(sizeof(UndoStruct));
    if (u != NULL) {
	u->command = UndoNullCommand;
	u->info.insert_bases.num_bases = 0;
	u->next = NULL;
    }
    
    return u;
}


void freeUndoStruct(UndoStruct *u, int free_op)
/*
 * Free an undo struct
 */
{
    if (u == NULL) return;
    
    /*
     * Free any allocated string space
     */
    switch (u->command) {
    case UndoInsertBases:
	unpackBCO(&u->info.insert_bases.b_c_o);
	break;
    case UndoReplaceBases:
	unpackBCO(&u->info.insert_bases.b_c_o);
	break;
	/*
	 * YUK! in some circumstances we want to free the tag,
	 * but most the time wo do not
	 */
    case UndoInsertAnnotation:
	if (free_op == FREE_RECLAIM)
	    freeTag(u->info.insert_annotation.new_tag);
	break;
    default:
	break;
    }
    
    xfree(u);
}



void freeUndoList(UndoStruct *u,int free_op)
/*
 * Free one command (an undo list)
 */
{
    UndoStruct *next;
    while (u != NULL) {
	next = u->next;
	freeUndoStruct(u,free_op);
	u = next;
    }
    
}


void freeUndoLists(EdStruct *xx, int free_op)
/*
 * Free all undos
 */
{
    for ( ; DBI_num_undo(xx) > 0 ; DBI_num_undo(xx)--,
	 DBI_last_undo(xx) = (DBI_last_undo(xx) + MAX_SAVE_EDITS - 1) %
	 MAX_SAVE_EDITS)
	freeUndoList(DBI_undo_lists(xx)[DBI_last_undo(xx)],free_op);
}


void resetEdits(EdStruct *xx) {
    DBI_since_auto_save(xx) = DBI_edits_made(xx) = 0;
}

void freeAllUndoLists(EdStruct *xx)
/*
 * Free all undos
 */
{
    freeUndoLists(xx, FREE_RECLAIM);
    DBI_discarded_undos(xx) = 0;
    resetEdits(xx);
}


void openUndo(DBInfo *db)
/*
 * Start recording undo for next command
 */
{
    if (++_DBI_open_undo_count(db) > 1)
	return;

    if (!_DBI_store_undo(db)) {
	/* Do not store undo info, but mark edits as having been made */
	_DBI_edits_made(db)++;
	_DBI_discarded_undos(db)++;
	return;
    }

    _DBI_last_undo(db) = (_DBI_last_undo(db) + 1) % MAX_SAVE_EDITS;
    
    if (_DBI_num_undo(db) == MAX_SAVE_EDITS) {
	/* must free up last_undo */
	freeUndoList(_DBI_undo_lists(db)[_DBI_last_undo(db)],FREE_RECLAIM);
	_DBI_discarded_undos(db) = 1;
    } else
	_DBI_num_undo(db)++;

    _DBI_edits_made(db)++;

    /* ensure this is null */
    _DBI_undo_lists(db)[_DBI_last_undo(db)] = NULL;
    
}


void closeUndo(EdStruct *xx, DBInfo *db)
/*
 * Start recording undo for next command
 */
{
    if (--_DBI_open_undo_count(db) != 0)
	return;

    if (!_DBI_store_undo(db))
	return;

    if (_DBI_num_undo(db) > 0 &&
	_DBI_undo_lists(db)[_DBI_last_undo(db)] == NULL){
	_DBI_last_undo(db) = (_DBI_last_undo(db) + MAX_SAVE_EDITS - 1) %
	    MAX_SAVE_EDITS;
	_DBI_num_undo(db)--;
	if (--_DBI_edits_made(db) < 0) {
	    _DBI_edits_made(db) = 0;
	}
    }
    
    if (xx->auto_save && _DBI_edits_made(db) >= AUTO_SAVE_EDITS) {
	vmessage("Contig Editor: auto saving\n");
	UpdateTextOutput();
	saveDB(xx, _DBI_io(db), 1 /* an autosave */, 1 /* notify others */);
    }
}


void recordUndo(DBInfo *db, UndoStruct *u)
/*
 * Record an edit operation for current command
 */
{
    if (!_DBI_store_undo(db))
	return;

    u->next = _DBI_undo_lists(db)[_DBI_last_undo(db)];
    _DBI_undo_lists(db)[_DBI_last_undo(db)] = u;
}

/*
 * Stop recording undo information
 */
void stopUndo(EdStruct *xx) {
    int edits_made = DBI_edits_made(xx);
    int num_undo = DBI_num_undo(xx);
    DBI_store_undo(xx) = 0;
    freeAllUndoLists(xx);
    DBI_edits_made(xx) = edits_made;
    DBI_discarded_undos(xx) = num_undo;
}

void startUndo(EdStruct *xx) {
    DBI_store_undo(xx) = 1;
}





/*************************************************************
 *
 *************************************************************/



void undoEdit(UndoStruct *u)
/*
 * Undo an edit of composite command
 */
{
    if (u==NULL) return;

    switch(u->command) {
    case UndoShiftRight:
	_shift_right(
		     u->db,
		     u->sequence,
		     u->info.shift_right.num_bases,
		     u->info.shift_right.flags
		     );
	break;
    case UndoShiftLeft:
	_shift_left(
		    u->db,
		    u->sequence,
		    u->info.shift_left.num_bases,
		    u->info.shift_left.flags
		    );
	break;
    case UndoDeleteBases:
	_delete_bases(
		      u->db,
		      u->sequence,
		      u->info.delete_bases.position,
		      u->info.delete_bases.num_bases,
		      u->info.delete_bases.flags
		      );
	break;
    case UndoReplaceBases: {
	char *bases;
	int1 *conf;
	int2 *opos;

	getBCO(&u->info.replace_bases.b_c_o, &bases, &conf, &opos);
	_replace_bases(
		       u->db,
		       u->sequence,
		       u->info.replace_bases.position,
		       u->info.replace_bases.num_bases,
		       bases,
		       conf,
		       opos,
		       u->info.replace_bases.flags & ~DB_FLAG_TMP,
		       u->info.replace_bases.flags & DB_FLAG_TMP_DIFF_ONLY,
		       u->info.replace_bases.flags & DB_FLAG_TMP_CONF_ONLY
		       );
	break;
    }
    case UndoInsertBases: {
	char *bases;
	int1 *conf;
	int2 *opos;

	getBCO(&u->info.replace_bases.b_c_o, &bases, &conf, &opos);

	_insert_bases(
		      u->db,
		      u->sequence,
		      u->info.insert_bases.position,
		      u->info.insert_bases.num_bases,
		      bases,
		      conf,
		      opos,
		      u->info.insert_bases.flags,
		      u->info.insert_bases.cutoff);
	break;
    }
    case UndoReorderSeq:
	_reorder_seq(
		     u->db,
		     u->sequence,
		     u->info.reorder_seq.old_id,
		     u->info.reorder_seq.new_id,
		     u->info.reorder_seq.flags
		     );
	break;
    case UndoConsensusLength:
	_DBsetLength(
		     u->db,
		     0,
		     u->info.consensus_length.num_bases
		     );
	_DBsetLength2(
		      u->db,
		      0,
		      u->info.consensus_length.num_bases
		      );
	break;
    case UndoAdjustCursor:
	/* Only do if the EdStruct hasn't been reallocated */
	if (u->info.adjust_cursor.xx->editor_id ==
	    u->info.adjust_cursor.editor_id) {
	    setCursorPosSeq(u->info.adjust_cursor.xx,
			    u->info.adjust_cursor.position,
			    u->sequence);
	}
	break;
	

    case UndoAdjustDisplay:
	/* Only do if the EdStruct hasn't been reallocated */
	if (u->info.adjust_display.xx->editor_id ==
	    u->info.adjust_display.editor_id) {
	    u->info.adjust_display.xx->displayPos +=
		u->info.adjust_display.position;
	}
	break;
	
	/*
	 * Annotation manipulation
	 */
    case UndoAdjustPositionAnnotation:
	_adjust_position_annotation(
				    u->db,
				    u->sequence,
				    u->info.adjust_position_annotation.tag,
				    u->info.adjust_position_annotation.position,
				    u->info.adjust_position_annotation.seq_flags,
				    u->info.adjust_position_annotation.tag_flags
				    );
	break;
    case UndoAdjustLengthAnnotation:
	_adjust_length_annotation(
				  u->db,
				  u->sequence,
				  u->info.adjust_length_annotation.tag,
				  u->info.adjust_length_annotation.length,
				  u->info.adjust_length_annotation.seq_flags,
				  u->info.adjust_length_annotation.tag_flags
				  );
	break;
	
    case UndoInsertAnnotation:
	_insert_annotation(
			   u->db,
			   u->sequence,
			   u->info.insert_annotation.tag,
			   u->info.insert_annotation.new_tag,
			   u->info.insert_annotation.seq_flags
			   );
	/* to unsure the following is not reclaimed during garbage collect */
	u->info.insert_annotation.new_tag = NULL;
	break;
	
    case UndoDeleteAnnotation:
	_delete_annotation(
			   u->db,
			   u->sequence,
			   u->info.delete_annotation.tag,
			   u->info.delete_annotation.seq_flags
			   );
	/* garbage collect */
	freeTag(u->info.delete_annotation.old_tag);
	break;
	
    case UndoDestroyAnnotation:
	_destroy_annotation(
			    u->db,
			    u->info.destroy_annotation.xx,
			    u->sequence,
			    u->info.destroy_annotation.tag,
			    u->info.destroy_annotation.seq_flags
			    );
	break;
	/*
	 * Edit manipulation
	 */
	
	
	
    case UndoAdjustEnds:
	_adjust_ends(
		     u->db,
		     u->sequence,
		     u->info.adjust_ends.start_bases,
		     u->info.adjust_ends.end_bases,
		     u->info.adjust_ends.seq_flags
		     );
	break;
	
    case UndoAdjustBaseConf:
	_adjust_base_conf(
			  u->db,
			  u->sequence,
			  u->info.adjust_base_conf.position,
			  u->info.adjust_base_conf.conf,
			  u->info.adjust_base_conf.opos,
			  u->info.adjust_base_conf.flags
			  );
	break;

    case UndoTransposeBases:
	_transpose_bases(
			 u->db,
			 u->sequence,
			 u->info.transpose_bases.position,
			 u->info.transpose_bases.flags
			 );
	break;

    case UndoSetFlags:
	_set_flags(
		   u->db,
		   u->sequence,
		   u->info.set_flags.flags
		   );
	break;

    case UndoSetReferenceSeq:
	_set_reference_seq(
		   u->db,
		   u->sequence,
		   u->info.set_reference_seq.flags,
		   u->info.set_reference_seq.refseq,
		   u->info.set_reference_seq.length,
		   u->info.set_reference_seq.offset
		   );
	break;

    default:
	break;
    }
}






void undoLastCommand(EdStruct *xx)
/*
 * Undo last keypress 
 */
{
    UndoStruct *u;
    
    if (xx->editorState == StateDown)
	return;

    if (DBI_num_undo(xx) == 0) {
	/* no undo available */
	bell();
	return;
    }
    
    u = DBI_undo_lists(xx)[DBI_last_undo(xx)];

    if (u == NULL || !(_DBI_flags(u->db) & DB_ACCESS_UPDATE)) {
	bell();
	return;
    }

    if (--DBI_edits_made(xx) < 0) {
	DBI_edits_made(xx) = 0;
	DBI_since_auto_save(xx) = 1;
    }
    
    while ( u != NULL ) {
	/* fprintf(stderr, "Undoing edit on db %p\n", u->db); */
	/* DumpUndo(u); */
	undoEdit(u);
	u = u->next;
    }
    
    
    /* remove last undo */
    freeUndoList(DBI_undo_lists(xx)[DBI_last_undo(xx)],FREE_UNDO);
    DBI_undo_lists(xx)[DBI_last_undo(xx)] = NULL;
    DBI_last_undo(xx) = (DBI_last_undo(xx) + MAX_SAVE_EDITS - 1) %
	MAX_SAVE_EDITS;
    DBI_num_undo(xx)--;

#ifdef CACHE_CONSENSUS_2
    /* Update consensus buffers */
    invalidate_consensus(xx);
#endif

    /* Not optimal, but works */
    xx->refresh_flags = ED_DISP_ALL;
    redisplayWithCursor(xx);
}






#if 0
static void DumpUndo(UndoStruct *u)
/*
 * Dump undo
 */
{
    switch(u->command) {
	
    case UndoShiftRight:
	fprintf(stderr,"   _shift_right(seq=%d,num_bases=%d)\n",
		u->sequence,
		u->info.shift_right.num_bases
		);
	break;
    case UndoShiftLeft:
	fprintf(stderr,"   _shift_left(seq=%d,num_bases=%d)\n",
		u->sequence,
		u->info.shift_left.num_bases
		);
	break;
    case UndoDeleteBases:
	fprintf(stderr,"   _delete_bases(seq=%d,pos=%d,num_bases=%d)\n",
		u->sequence,
		u->info.delete_bases.position,
		u->info.delete_bases.num_bases
		);
	break;
    case UndoReplaceBases:
	fprintf(stderr,"   _replace_bases(seq=%d,pos=%d,bases=\"%.*s\")\n",
		u->sequence,
		u->info.replace_bases.position,
		u->info.replace_bases.num_bases,
		/* getChar(&u->info.replace_bases.bases)*/ "unknown"
		);
	break;
    case UndoInsertBases:
	fprintf(stderr,"   _insert_bases(seq=%d,pos=%d,bases=\"%.*s\")\n",
		u->sequence,
		u->info.insert_bases.position,
		u->info.insert_bases.num_bases,
		/* getChar(&u->info.insert_bases.bases) */ "unknown"
		);
	break;
    case UndoReorderSeq:
	fprintf(stderr,"   _reorder_seq(seq=%d,old_id=%d,new_id=%d)\n",
		u->sequence,
		u->info.reorder_seq.old_id,
		u->info.reorder_seq.new_id
		);
	break;
    case UndoConsensusLength:
	fprintf(stderr,"   consensus length=%d\n",
		u->info.consensus_length.num_bases
		);
	break;
    case UndoAdjustCursor:
	fprintf(stderr,"   adjust_cursor pos=%d\n",
		u->info.adjust_cursor.position
		);
	break;
    case UndoAdjustDisplay:
	fprintf(stderr,"   adjust_display(editor_id=%d,position=%d)\n",
		u->info.adjust_display.editor_id,
		u->info.adjust_display.position
		);
    case UndoAdjustPositionAnnotation:
	fprintf(stderr,"   _adjust_position_annotation(tag=%p,pos=%d)\n",
		(void *)u->info.adjust_position_annotation.tag,
		u->info.adjust_position_annotation.position
		);
	break;
    case UndoAdjustLengthAnnotation:
	fprintf(stderr,"   _adjust_length_annotation(tag=%p,len=%d)\n",
		(void *)u->info.adjust_length_annotation.tag,
		u->info.adjust_length_annotation.length
		);
	break;
    case UndoInsertAnnotation:
	fprintf(stderr,"   _insert_annotation(seq=%d,tag=%p,newtag=%p)\n",
		u->sequence,
		(void *)u->info.insert_annotation.tag,
		(void *)u->info.insert_annotation.new_tag
		);
	break;
    case UndoDeleteAnnotation:
	fprintf(stderr,"   _delete_annotation(seq=%d,tag=%p,old_tag=%p)\n",
		u->sequence,
		(void *)u->info.delete_annotation.tag,
		(void *)u->info.delete_annotation.old_tag
		);
	break;
    case UndoDestroyAnnotation:
	fprintf(stderr,"   _destroy_annotation(seq=%d,tag=%p)\n",
		u->sequence,
		(void *)u->info.destroy_annotation.tag
		);
	break;
	
	
	/*
	 * Modify raw data parameters
	 */
    case UndoAdjustEnds:
	fprintf(stderr,"   _adjust_ends(start=%d,end=%d)\n",
		u->info.adjust_ends.start_bases,
		u->info.adjust_ends.end_bases
		);
	
	break;
	
    case UndoAdjustBaseConf:
	fprintf(stderr, "  _adjust_base_conf(seq=%d,pos=%d,val=%d)\n",
		u->sequence,
		u->info.adjust_base_conf.position,
		u->info.adjust_base_conf.conf
		);
	break;

    case UndoTransposeBases:
	fprintf(stderr,"   transposebases(seq=%d,pos=%d)\n",
		u->sequence,
		u->info.transpose_bases.position
		);
	break;

    default:
	fprintf(stderr,"   unknown command\n");
	break;
    }
}
#endif



/* YUK! name archaic */
void cleanUpStacks(EdStruct *xx)
/*
 * free all stacks
 */
{
    freeAllUndoLists(xx);
}

/* YUK! name archaic */
void cleanUpAllStacks(EdStruct *xx)
/*
 * free all stacks
 */
{
    cleanUpStacks(xx);
}



