/*
 * File: undo.h
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

#ifndef _UNDO_H_
#define _UNDO_H_

#include "edUtils.h"


enum UndoCommands {
    UndoNullCommand,
    UndoShiftRight,
    UndoShiftLeft,
    UndoDeleteBases,
    UndoReplaceBases,
    UndoInsertBases,
    UndoReorderSeq,
    UndoConsensusLength,
    UndoAdjustCursor,
    UndoAdjustPositionEdit,
    UndoAdjustPositionAnnotation,
    UndoAdjustLengthAnnotation,
    UndoDeleteAnnotation,
    UndoInsertAnnotation,
    UndoDestroyAnnotation,
    UndoAdjustEnds,
    UndoAdjustBaseConf,
    UndoTransposeBases,
    UndoAdjustDisplay,
    UndoSetFlags,
    UndoSetReferenceSeq
};

typedef struct {
    union {
	char array[sizeof(char *)];
	char *ptr;
    } data;
    int len;
} PackedBCO;

extern void packBCO(PackedBCO *ps, char *bases, int1 *conf,
		    int2 *opos, int len);


typedef struct bnode *undoStructP;

typedef struct undoNode {
    DBInfo *db;
    struct undoNode *next;
    int command;
    int sequence;

    union {

	struct {
	    int flags;
	    int num_bases;
	} shift_left, shift_right;

	struct {
	    int flags;
	    int position;
	    int num_bases;
	} delete_bases;

	struct {
	    PackedBCO b_c_o;
	    int flags;
	    int position;
	    int num_bases;
	} replace_bases;

	struct {
	    PackedBCO b_c_o;
	    int flags;
	    int position;
	    int num_bases;
	    int cutoff;
	} insert_bases;

	struct {
	    int flags;
	    int old_id;
	    int new_id;
	} reorder_seq;

	struct {
	    int num_bases;
	} consensus_length;

	struct {
	    EdStruct *xx;
	    int position;
	    int editor_id;
	} adjust_cursor, adjust_display;


	struct {
	    tagStruct *tag;
	    EdStruct *xx;
	    int seq_flags;
	} destroy_annotation;

	struct {
	    tagStruct *tag;
	    int position;
	    int tag_flags;
	    int seq_flags;
	} adjust_position_annotation;


	struct {
	    tagStruct *tag;
	    int length;
	    int tag_flags;
	    int seq_flags;
	} adjust_length_annotation;

	struct {
	    tagStruct *tag;
	    tagStruct *new_tag;
	    int seq_flags;
	} insert_annotation;

	struct {
	    tagStruct *tag;
	    tagStruct *old_tag;
	    int seq_flags;
	} delete_annotation;

	struct {
	    int start_bases;
	    int end_bases;
	    int seq_flags;
	} adjust_ends;

	struct {
	    int position;
	    int flags;
	    int1 conf;
	    int2 opos;
	} adjust_base_conf;

	struct {
	    int position;
	    int flags;
	} transpose_bases;

	struct {
	    int flags;
	} set_flags;
	
	struct {
	    int flags;
	    int refseq;
	    int length;
	    int offset;
	} set_reference_seq;

    } info;

} UndoStruct;

extern UndoStruct *newUndoStruct(DBInfo *db);
extern void recordUndo(DBInfo *db, UndoStruct *u);
int _editsMade(DBInfo *db);
int editsMade(EdStruct *xx);
void undoLastCommand(EdStruct *xx);
void resetEdits(EdStruct *xx);
void freeAllUndoLists(EdStruct *xx);
void openUndo(DBInfo *db);
void closeUndo(EdStruct *xx, DBInfo *db);
void stopUndo(EdStruct *xx);
void startUndo(EdStruct *xx);

#endif  /*_UNDO_H_*/
