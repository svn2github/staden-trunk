/*  Last edited: Jan  7 10:34 2004 (mng) */
/*
 * File: edUtils.h
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

#ifndef _edUtils_h
#define _edUtils_h

#include "qual.h"
#include "edStructs.h"
#include "tagUtils.h"
#include "IO.h"
#include "tman_display.h"
#include "dstring.h"

/*
 * Useful distances
 * (treat as symbolic rather than actual distances)
 */
#define D_screen     80
#define D_halfScreen 40
#define D_character   1

/*
 * Flags for the internal database
 */
/* for whole database */
#define DB_ACCESS            (1L<<0)
#define DB_ACCESS_READONLY   (0L)
#define DB_ACCESS_UPDATE     (DB_ACCESS)

#define DB_DATA_TYPE         (1L<<1)
#define DB_DATA_TYPE_DNA     (0L)
#define DB_DATA_TYPE_PROTEIN (DB_DATA_TYPE)

#define DB_STORAGE           (1L<<2)
#define DB_STORAGE_DISK      (0L)
#define DB_STORAGE_INTERNAL  (DB_STORAGE)
#define DB_DELAYED_READ      (0L)

#define DB_NO_REGS	     (1L<<3)


/* for each sequence */
#define DB_FLAG_NONE             (0L)
#define DB_FLAG_IN_MEMORY        (1L<<0)
#define DB_FLAG_SEQ_MODIFIED     (1L<<1)
#define DB_FLAG_REL_MODIFIED     (1L<<2)
#define DB_FLAG_TAG_MODIFIED     (1L<<3)
#define DB_FLAG_SELECTED         (1L<<4)
#define DB_FLAG_TAG_IN_MEMORY    (1L<<5)
#define DB_FLAG_SEQ_IN_MEMORY    (1L<<0)
#define DB_FLAG_NAME_IN_MEMORY   (1L<<6)
#define DB_FLAG_TRACE_SHOWN	 (1L<<7)
#define DB_FLAG_TERMINATOR	 (1L<<8)
#define DB_FLAG_INVIS		 (1L<<9)
#define DB_FLAG_REFTRACE_NEG	 (1L<<10)
#define DB_FLAG_REFTRACE_POS	 (1L<<11)
#define DB_FLAG_REFTRACE	 (DB_FLAG_REFTRACE_NEG | DB_FLAG_REFTRACE_POS)
#define DB_FLAG_REFSEQ		 (1L<<12)
#define DB_FLAG_NOTE_MODIFIED	 (1L<<13)

#define DB_FLAG_TMP_DIFF_ONLY	 (1L<<20)
#define DB_FLAG_TMP_CONF_ONLY	 (1L<<21)
#define DB_FLAG_TMP		 (DB_FLAG_TMP_DIFF_ONLY|DB_FLAG_TMP_CONF_ONLY)

/*
 * Callback types for DBI_dispFunc()
 */
#define DBCALL_REDISPLAY	0
#define DBCALL_INSERT		1
#define DBCALL_DELETE		2
#define DBCALL_CURSOR		3
#define DBCALL_ADJUST_START	4
#define DBCALL_REINIT		5
#define DBCALL_JOIN_SHIFT	6
#define DBCALL_QUIT		7
#define DBCALL_CURSOR_NOTIFY	8
#define DBCALL_RELINK	        9

/*
 * Superedit flags
 */
#define SUPEREDIT_INS_READ	(1<<0)
#define SUPEREDIT_DEL_READ	(1<<1)
#define SUPEREDIT_INS_ANY_CON	(1<<2)
#define SUPEREDIT_DEL_DASH_CON	(1<<3)
#define SUPEREDIT_DEL_ANY_CON	(1<<4)
#define SUPEREDIT_REPLACE_CON	(1<<5)
#define SUPEREDIT_MODIFY_CONF	(1<<6)
#define SUPEREDIT_TRANSPOSE_ANY	(1<<7)
#define SUPEREDIT_SHIFT_READ	(1<<8)
#define SUPEREDIT_UPPERCASE	(1<<9)

/*
 * Display update flags - see xx->refresh_flags
 */
#define ED_DISP_NAMES		(1<<0)
#define ED_DISP_READS		(1<<1)
#define ED_DISP_CONS		(1<<2)
#define ED_DISP_STATUS		(1<<4)
#define ED_DISP_RULER		(1<<5)
#define ED_DISP_CURSOR		(1<<6)
#define ED_DISP_SELECTION	(1<<7)
#define ED_DISP_HEIGHT		(1<<8)
#define ED_DISP_YSCROLL		(1<<9)
#define ED_DISP_SCROLL		((1<<3) | ED_DISP_YSCROLL)
/* Sequence and consensus */
#define ED_DISP_SEQS		(ED_DISP_READS | ED_DISP_CONS)
/* Everything */
#define ED_DISP_ALL		(ED_DISP_SEQS | ED_DISP_NAMES | \
				 ED_DISP_SCROLL | ED_DISP_STATUS | \
				 ED_DISP_RULER | ED_DISP_CURSOR | \
				 ED_DISP_SELECTION | ED_DISP_HEIGHT)
/* Single items only - see xx->refresh_seq */
#define ED_DISP_NAME		(1<<10)
#define ED_DISP_SEQ		(1<<11)
/* For optimising redisplayDisagreements to remove flicker */
#define ED_DISP_NO_DIFFS	(1<<12)


/****************/
#define DB_RelPos(X,A)   ((X)->DBi->DB[(A)].relPos)
#define DB_Length(X,A)   ((X)->DBi->DB[(A)].length) /* clipped length */
#define DB_Number(X,A)   ((X)->DBi->DB[(A)].number)
#define DB_Comp(X,A)     ((X)->DBi->DB[(A)].complemented)
#define DB_Name(X,A)     ((X)->DBi->DB[(A)].name)
#define DB_Seq(X,A)      ((X)->DBi->DB[(A)].sequence)
#define DB_Flags(X,A)    ((X)->DBi->DB[(A)].flags)
#define DB_Tags(X,A)     ((X)->DBi->DB[(A)].tagList)
#define DB_Opos(X,A)     ((X)->DBi->DB[(A)].gap_opos)
#define DB_Conf(X,A)     ((X)->DBi->DB[(A)].gap_conf)
#define DB_Start(X,A)    ((X)->DBi->DB[(A)].gap_start)
#define DB_End(X,A)      ((X)->DBi->DB[(A)].gap_end)
#define DB_Length2(X,A)  ((X)->DBi->DB[(A)].gap_length) /* full length */
#define DB_Alloced(X,A)  ((X)->DBi->DB[(A)].alloc_length)
#define DB_WSeq(X,A)     (DB_Seq((X),(A))+DB_Start((X),(A)))

#define _DB_RelPos(DBi,A)   ((DBi)->DB[(A)].relPos)
#define _DB_Length(DBi,A)   ((DBi)->DB[(A)].length)
#define _DB_Number(DBi,A)   ((DBi)->DB[(A)].number)
#define _DB_Comp(DBi,A)     ((DBi)->DB[(A)].complemented)
#define _DB_Name(DBi,A)     ((DBi)->DB[(A)].name)
#define _DB_Seq(DBi,A)      ((DBi)->DB[(A)].sequence)
#define _DB_Flags(DBi,A)    ((DBi)->DB[(A)].flags)
#define _DB_Tags(DBi,A)     ((DBi)->DB[(A)].tagList)
#define _DB_Opos(DBi,A)     ((DBi)->DB[(A)].gap_opos)
#define _DB_Conf(DBi,A)     ((DBi)->DB[(A)].gap_conf)
#define _DB_Start(DBi,A)    ((DBi)->DB[(A)].gap_start)
#define _DB_End(DBi,A)      ((DBi)->DB[(A)].gap_end)
#define _DB_Length2(DBi,A)  ((DBi)->DB[(A)].gap_length)
#define _DB_Alloced(DBi,A)  ((DBi)->DB[(A)].alloc_length)
#define _DB_WSeq(DBi,A)     (_DB_Seq((DBi),(A))+_DB_Start((DBi),(A)))

#define DBI(X)		  ((X)->DBi)
#define DBI_DB(X)         ((X)->DBi->DB)
#define DBI_flags(X)      ((X)->DBi->DB_flags)
#define DBI_gelCount(X)   ((X)->DBi->DB_gelCount)
#define DBI_contigNum(X)  ((X)->DBi->DB_contigNum)
#define DBI_list(X)       ((X)->DBi->DBlist)
#define DBI_order(X)      ((X)->DBi->DBorder)
#define DBI_dispFunc(X)   ((X)->DBi->display_func)
#define DBI_dispData(X)   ((X)->DBi->display_data)
#define DBI_nextDisp(X)   ((X)->DBi->next_display)
#define DBI_io(X)	  ((X)->DBi->io)
#define DBI_undo_lists(X) ((X)->DBi->undo_lists)
#define DBI_last_undo(X)	((X)->DBi->last_undo)
#define DBI_num_undo(X)		((X)->DBi->num_undo)
#define DBI_discarded_undos(X)	((X)->DBi->discarded_undos)
#define DBI_edits_made(X)	((X)->DBi->edits_made)
#define DBI_since_auto_save(X)	((X)->DBi->since_auto_save)
#define DBI_store_undo(X)	((X)->DBi->store_undo)
#define DBI_open_undo_count(X)	((X)->DBi->open_undo_count)
#define DBI_registration_id(X)	((X)->DBi->registration_id)

#define _DBI_DB(DBi)         ((DBi)->DB)
#define _DBI_flags(DBi)      ((DBi)->DB_flags)
#define _DBI_gelCount(DBi)   ((DBi)->DB_gelCount)
#define _DBI_contigNum(DBi)  ((DBi)->DB_contigNum)
#define _DBI_list(DBi)       ((DBi)->DBlist)
#define _DBI_order(DBi)      ((DBi)->DBorder)
#define _DBI_dispFunc(DBi)   ((DBi)->display_func)
#define _DBI_dispData(DBi)   ((DBi)->display_data)
#define _DBI_nextDisp(DBi)   ((DBi)->next_display)
#define _DBI_io(DBi)         ((DBi)->io)
#define _DBI_undo_lists(DBi) ((DBi)->undo_lists)
#define _DBI_last_undo(DBi)		((DBi)->last_undo)
#define _DBI_num_undo(DBi)		((DBi)->num_undo)
#define _DBI_discarded_undos(DBi)	((DBi)->discarded_undos)
#define _DBI_edits_made(DBi)		((DBi)->edits_made)
#define _DBI_since_auto_save(DBi)	((DBi)->since_auto_save)
#define _DBI_store_undo(DBi)		((DBi)->store_undo)
#define _DBI_open_undo_count(DBi)	((DBi)->open_undo_count)
#define _DBI_registration_id(DBi)	((DBi)->registration_id)

/****************/
#define COMPLEMENTED -1
#define BOTH_STRANDS 0
#define UNCOMPLEMENTED 1

#define DBgetGelName(xx,i) ( &( DBgetName(DBI(xx),i) )[DB_GELNOLEN+1] )

#define DBsetRelPos(X,A,B)  (X)->DBi->DB[(A)].relPos = (B)
#define DBsetLength(X,A,B)  (X)->DBi->DB[(A)].length = (B)
#define DBsetNumber(X,A,B)  (X)->DBi->DB[(A)].number = (B)
#define DBsetComp(X,A,B)    (X)->DBi->DB[(A)].complemented = (B)
#define DBsetName(X,A,B)    (X)->DBi->DB[(A)].name = (B)
#define DBsetSeq(X,A,B)     (X)->DBi->DB[(A)].sequence = (B)
#define DBsetFlags(X,A,B)   (X)->DBi->DB[(A)].flags = (B)
#define DBaddFlags(X,A,B)   (X)->DBi->DB[(A)].flags |= (B)
#define DBclearFlags(X,A,B) (X)->DBi->DB[(A)].flags &= ~(B)
#define DBsetTags(X,A,B)    (X)->DBi->DB[(A)].tagList = (B)
#define DBsetOpos(X,A,B)    (X)->DBi->DB[(A)].gap_opos = (B)
#define DBsetConf(X,A,B)    (X)->DBi->DB[(A)].gap_conf = (B)
#define DBsetStart(X,A,B)   (X)->DBi->DB[(A)].gap_start = (B)
#define DBsetEnd(X,A,B)     (X)->DBi->DB[(A)].gap_end = (B)
#define DBsetLength2(X,A,B) (X)->DBi->DB[(A)].gap_length = (B)
#define DBsetAlloced(X,A,B) (X)->DBi->DB[(A)].alloc_length = (B)

#define _DBgetGelName(DBi,i) ( &( DBgetName((DBi),i) )[DB_GELNOLEN+1] )

#define _DBsetRelPos(DBi,A,B)  (DBi)->DB[(A)].relPos = (B)
#define _DBsetLength(DBi,A,B)  (DBi)->DB[(A)].length = (B)
#define _DBsetNumber(DBi,A,B)  (DBi)->DB[(A)].number = (B)
#define _DBsetComp(DBi,A,B)    (DBi)->DB[(A)].complemented = (B)
#define _DBsetName(DBi,A,B)    (DBi)->DB[(A)].name = (B)
#define _DBsetSeq(DBi,A,B)     (DBi)->DB[(A)].sequence = (B)
#define _DBsetFlags(DBi,A,B)   (DBi)->DB[(A)].flags = (B)
#define _DBsetTags(DBi,A,B)    (DBi)->DB[(A)].tagList = (B)
#define _DBsetOpos(DBi,A,B)    (DBi)->DB[(A)].gap_opos = (B)
#define _DBsetConf(DBi,A,B)    (DBi)->DB[(A)].gap_conf = (B)
#define _DBsetStart(DBi,A,B)   (DBi)->DB[(A)].gap_start = (B)
#define _DBsetEnd(DBi,A,B)     (DBi)->DB[(A)].gap_end = (B)
#define _DBsetLength2(DBi,A,B) (DBi)->DB[(A)].gap_length = (B)
#define _DBsetAlloced(DBi,A,B) (DBi)->DB[(A)].alloc_length = (B)

/*
 * Useful macros
 */
#define normalisePos(X,S,P,L) \
    ( (DB_Comp((X),(S))==UNCOMPLEMENTED) ? (P) : (DB_Length((X),(S)) - (P) - (L) + 2) )

#define normalisePos2(X,S,P,L) \
    ( (DB_Comp((X),(S))==UNCOMPLEMENTED) ? (P) : (DB_Length2((X),(S)) - (P) - (L) + 2) )

#define normaliseSense(X,S,D) \
    ( (DB_Comp((X),(S))==UNCOMPLEMENTED) ? (D) : ((D) == 2 ? 2 : ((D)==0)) )

#define db_version (DBI_io(ss)->db.version)

#define editModeIsInsert(xx) (xx->insert_mode)

#define RedisplaySeq(xx,seq) \
((xx->refresh_seq > 0 && xx->refresh_seq != seq) ? \
 (xx->refresh_flags |= ED_DISP_SEQS | ED_DISP_STATUS) : \
 (xx->refresh_seq = seq, \
  xx->refresh_flags |= (ED_DISP_SEQ | ED_DISP_CONS | ED_DISP_STATUS))) \

#define RedisplayName(xx,seq) \
((xx->refresh_seq > 0 && xx->refresh_seq != seq) ? \
 (xx->refresh_flags = ED_DISP_NAMES) : \
 (xx->refresh_seq = seq, xx->refresh_flags |= ED_DISP_NAME)) \

/*----------------------------------------------------------------------------
 * Prototypes
 *---------------------------------------------------------------------------*/

extern void resize_consensus(EdStruct *xx, int len);
extern void invalidate_consensus(EdStruct *xx);
extern void invalidate_consensus_region(EdStruct *xx, int from, int to);
extern int valid_consensus(EdStruct *xx, int from, int to);

extern EdStruct *getFreeEdStruct(GapIO *io, int contig,
				 void (*dispFunc)(void *, int, int,
						  int, void *));
extern int createEdDisplay(EdStruct *xx, int seq, int pos);
extern int calculate_consensus_length(EdStruct *xx);
extern void calculateConsensusLength(EdStruct *xx);
extern int contEd_info(int job, void *mydata, info_arg_t *theirdata);
extern void DBcalcConsensus (EdStruct *xx,int pos, int width, char *str,
			     float *qual, int mode);
extern char *DBgetSeq(DBInfo *db, int seq);
extern tagStruct *DBgetTags (DBInfo *db, int seq);
extern char *DBgetName(DBInfo *db, int seq);
extern void DBgetSequence(EdStruct *xx, int seq, int pos, int width, char *str);
extern int initialiseDB/* FORIO */(
		 EdStruct *xx,
	         GapIO *io,
		 int cnum,	/* contig number */
		 int idbsiz,	/* size of database */
		 int llino	/* left-most gel in contig */
		 );
extern void freeDB(EdStruct *xx, int delete_ed);
extern void saveDB(
	    EdStruct *xx,
	    GapIO *io,
	    int   auto_save,
	    int notify)
;
extern void joinDB(
	    EdStruct *xx[2],
	    GapIO *io
	    );
extern int *sequencesInRegion(EdStruct *xx,int pos, int width);
extern int *sequencesOnScreen(EdStruct *xx,int pos, int width);
extern int positionInContig(EdStruct *xx, int seq, int pos);
extern void incDisplayPos (EdStruct *xx, int distance);
extern void decDisplayPos (EdStruct *xx, int distance);
extern void setDisplayPosPercent (EdStruct *xx, float percent);
extern void setDisplayPos(EdStruct *xx, int pos);
extern void setDisplayYPos(EdStruct *xx, int pos);
extern void redisplayWithCursor(EdStruct *xx);
void setCursorPos(EdStruct *xx, int pos);
void setCursorSeq(EdStruct *xx, int seq);
void setCursorPosSeq(EdStruct *xx, int pos, int seq);

extern void edInvokeTrace(EdStruct *xx);
extern void repositionTraces(EdStruct *xx);
extern int _cursor_set(DBInfo *db, int seq, int pos);
extern int _shift_right(DBInfo *db, int seq, int num_bases, int seq_flags);
extern int _shift_left(DBInfo *db, int seq, int num_bases, int seq_flags);
extern int _insert_bases(DBInfo *db, int seq, int pos, int num_bases,
			 char *bases,
			 int1 *conf, int2 *opos, int seq_flags, int cutoff);
extern int _delete_bases(DBInfo *db, int seq, int pos, int num_bases,
			 int seq_flags);
extern int _adjust_base_conf(DBInfo *db, int seq, int pos, int val, int opos,
			     int flags);
extern int _replace_bases(DBInfo *db, int seq, int pos, int num_bases,
			  char *bases, int1 *conf, int2 *opos, int seq_flags,
			  int diff_only, int conf_only);
extern int _transpose_bases(DBInfo *db, int seq, int pos, int seq_flags);
extern int _reorder_seq(DBInfo *db, int seq, int old_index, int new_index,
			int seq_flags);
extern int _set_reference_seq(DBInfo *db, int seq, int flags,
			      int refseq, int length, int offset);
extern int _set_flags(DBInfo *db, int seq, int seq_flags);
extern int U_shift_right(DBInfo *db, int seq, int num_bases);
extern int U_shift_left(DBInfo *db, int seq, int num_bases);
extern int U_insert_bases(EdStruct *xx, int seq, int pos, int num_bases, char *bases);
extern int U_delete_bases(EdStruct *xx, int seq, int pos, int num_bases);
extern int U_replace_bases(EdStruct *xx, int seq, int pos, int num_bases, char *bases, int diff_only);
extern int U_reorder_seq(EdStruct *xx, int seq, int old_index, int new_index);
extern void U_change_consensus_length(EdStruct *xx, int len);
extern void U_adjust_cursor(EdStruct *xx, int n);
extern void U_adjust_base_conf(EdStruct *xx, int seq, int pos, int val);
extern void U_transpose_bases(EdStruct *xx, int seq, int pos);
extern void U_adjust_display(EdStruct *xx, int n);
extern int handle_delete_bases(EdStruct *xx, int seq, int pos, int num_bases);
extern int handle_insert_bases(EdStruct *xx, int seq, int pos, int num_bases,
			char *bases);
extern int shiftRight(EdStruct *xx, int seq, int num_bases);
extern int shiftLeft(EdStruct *xx, int seq, int num_bases);
extern int insertBases(EdStruct *xx, int seq, int pos, int num_bases, char *bases);
extern int deleteBases(EdStruct *xx, int seq, int pos, int num_bases);
extern int replaceBases(EdStruct *xx, int seq, int pos, int num_bases, char *bases);
extern int replaceBasesConsensus(EdStruct *xx, int pos, int num_bases, char *bases);
extern int insertBasesConsensus(EdStruct *xx, int pos, int num_bases, char *bases);
extern int deleteBasesConsensus(EdStruct *xx, int pos, int num_bases);
extern int adjustBaseConf(EdStruct *xx, int seq, int pos, int val, int move);
extern int transpose(EdStruct *xx, int seq, int pos, int dir, int num_bases);
extern int set_reference_trace(EdStruct *xx, int seq, int flags);
extern int set_reference_seq(EdStruct *xx, int seq,
			     int refseq, int length, int offset);
extern int edTransposeLeft(EdStruct *xx, int num_bases);
extern int edTransposeRight(EdStruct *xx, int num_bases);
extern int edKeyPress(EdStruct *xx, char key, int nomove);
extern int edKeyDelete(EdStruct *xx);
extern int edKeyDeleteLeft(EdStruct *xx);
extern int edConf100(EdStruct *xx);
extern int edConf0(EdStruct *xx);
extern int edConfIncr(EdStruct *xx, int amount);
extern int edSelectRead(EdStruct *xx, int seq, int state);
extern void edSelectClear(EdStruct *xx);
extern void edSelectFrom(EdStruct *xx, int pos);
extern void edSelectTo(EdStruct *xx, int pos);
extern int edStartRead(EdStruct *xx);
extern int edEndRead(EdStruct *xx);
extern int edStartRead2(EdStruct *xx);
extern int edEndRead2(EdStruct *xx);
extern int edStartContig(EdStruct *xx);
extern int edEndContig(EdStruct *xx);
extern int edCursorRight(EdStruct *xx);
extern int edCursorLeft(EdStruct *xx);
extern int edCursorDown(EdStruct *xx);
extern int edCursorUp(EdStruct *xx);
extern void countDisagreements(EdStruct *xx[2], int *overlapLength,
			       int *wingeCount, int *ptgood, int *ptbad);
extern int posToIndex(EdStruct *xx, int pos);
extern int posToSeq(EdStruct *xx, int pos);
extern int seqToIndex(EdStruct *xx, int seq);
extern int origpos(EdStruct *xx, int seq, int pos);
extern void getLeftCutOff(EdStruct *xx,int seq, int width, char *str);
extern void getLCut(EdStruct *xx,int seq, int pos, int width, char *str);
extern void getRightCutOff(EdStruct *xx,int seq, int width, char *str);
extern void getRCut(EdStruct *xx,int seq, int pos, int width, char *str);
extern int lenRCut(EdStruct *xx, int seq);
extern int lenLCut(EdStruct *xx, int seq);
extern int linesOnScreen (EdStruct *xx, int pos, int width);
extern int onScreen (EdStruct *xx, int seq, int pos, int *wrong_x);
extern EdLink *CreateEdLink(EdStruct *xx0, EdStruct *xx1);
extern void redisplaySequences(EdStruct *xx, int update_all);
extern void redisplayDBSequences(DBInfo *db, int names_only);
extern int edSetRevealCutoffs(EdStruct *xx, int reveal);
extern int edSetInsertMode(EdStruct *xx, int insert);
extern int edSetSuperedit(EdStruct *xx, int superedit);
extern int edGetGelNumber(EdStruct *xx, int x, int y);
extern char *edGetGelName(EdStruct *xx, int number);
extern dstring_t *edGetGelNamesToRight(EdStruct *xx, int number);
extern int edGetStatusNumber(EdStruct *xx, int x, int y);
extern void edSetCMode(EdStruct *xx, int value);
extern void edSetCCutoff(EdStruct *xx, int value);
extern void edSetQCutoff(EdStruct *xx, int value);
extern void edSetDifferenceQuality(EdStruct *xx, int value);
extern int edSetJoinLock(EdStruct *xx, int value);
extern void edJoin(EdStruct *xx);
extern void edJoinAlign(EdStruct *xx);
extern void edSetActiveAnnos(EdStruct *xx, int argc, char **argv);
extern void edSetMiniTraces(EdStruct *xx, int height);
extern void edNextDifference(EdStruct *xx, int fwd);
extern int getQual(EdStruct *xx, int seq, int pos);
extern int edDoSearch(EdStruct *xx, int forwards, int strand, char
		      *type, char *value);
extern int edHideRead(EdStruct *xx, int seq, int check_cursor);

void status_strand(EdStruct *xx, int pos, int width,
		   XawSheetInk *splodge, char *sline,
		   char *name);

void front_editor(EdStruct *xx);
int move_editor(int id, int seq, int pos);
int editor_available(int contig, int nojoin);
EdStruct *editor_id_to_edstruct(int id);

void setDisplayPosP(EdStruct *xx, int pos);
void setDisplayPos(EdStruct *xx, int pos);
void setDisplayPos2(EdStruct *xx, int pos);
void setDisplayYPos(EdStruct *xx, int pos);
int rnum_to_edseq(EdStruct *xx, int rnum);
int inJoinMode(EdStruct *xx);
int editorLocked(EdStruct *xx);
void ed_set_slider_pos(EdStruct *xx, int pos);
void ed_set_yslider_pos(EdStruct *xx, int pos);
void ed_set_slider_pos(EdStruct *xx, int pos);
void bell(void);
void edStatusAdd(EdStruct *xx, int mode);
void edStatusDelete(EdStruct *xx, int mode);
void edStatusTransMode(EdStruct *xx, int mode);
void edShowQuality(EdStruct *xx, int mode);
void edShowCQuality(EdStruct *xx, int mode);
void edShowEdits(EdStruct *xx, int mode);
void edSetTraceComparator(EdStruct *xx, int seq);
int edSetCursor(EdStruct *xx, int x, int y);
void edSetCursorConsensus(EdStruct *xx, int pos);
int edZapLeft(EdStruct *xx);
int edZapRight(EdStruct *xx);
int edExtendLeft(EdStruct *xx);
int edExtendRight(EdStruct *xx);
void delete_edStruct(EdStruct *xx);
void extents(EdStruct *xx, int *left, int *right);
int NumberToSeq(DBInfo *db, int number);
int redisplayDisagreement(EdStruct *xx);
void positionCursor(EdStruct *xx, int seq, int pos);
DisplayContext *showTrace(EdStruct *xx, int seq, int pos, int baseSpacing,
	      int differencing, int mini_trace);
int get_trace_path(EdStruct *xx, int seq, char *fileName, char *t_type);
void editor_send_pos(int contig);
void DestroyEdLink(EdLink *el);
void setNamePos(EdStruct *xx, int pos);

/*
 * Set the brief status line with whatever is relevant for the data underneath
 * the current cursor position.
 */
int edSetBriefSeqStatus(EdStruct *xx, int x, int y);

/*
 * Set the brief status line with whatever is relevant for the data underneath
 * the current cursor position.
 */
int edSetBriefNameStatus(EdStruct *xx, int x, int y);

/*
 * List the confidence values in a tabular form. Based on the list_confidence
 * command in newgap_cmds.c.
 */
int edListConfidence(EdStruct *xx, int start, int end, int info_only);

/*
 * ensure that the cursor is visible on the screen
 */
void showCursor(EdStruct *xx, int seq, int pos);

/*
 * Get maximum extents of sequence allowing for cutoff and join mode.
 */
void getExtents(EdStruct *xx);

/*
 * Return number of sequences on screen
 */
int linesInRegion(EdStruct *xx, int pos, int width);

/*
 * Return number of sequences on screen
 */
int linesOnScreen (EdStruct *xx, int pos, int width);

/*
 * -- from join.c --
 * Align the two contig editor windows
 * Returns:
 *	0 - aligned ok
 *	1 - not ok
 */
int alignOverlap(EdStruct *xx[2]);

/*
 * -- from tkEdUtils.c --
 * Creates a selection spanning the length of the tag underneath the cursor
 */
void _select_tag(EdStruct *xx, int seq, tagStruct *t);

/*
 * -- from tkEdUtils.c --
 * Updates the 'brief' (NB: not a 'status' line in the same sense as the
 * above uses) line.
 *
 * To optimise usage of this, it returns a 'code number' which is updated
 * each time a new message is displayed. Passing 'msg' as NULL will return
 * the current code. Therefore calling routines can optimise by knowing
 * whether or not an update is required.
 */
int tk_update_brief_line(EdStruct *xx, char *msg);

/*
 * -- from tkEdUtils.c --
 * Automatically called when X wishes to obtain a selection. We register
 * this procedure in the initialise code of tkEditor.c.
 *
 * Return codes expected:
 *    -1  Failure
 *    >0  Number of bytes
 */
int edGetSelection(ClientData clientData, int offset, char *buffer,
		   int bufsize);

/*
 * -- from tkEdUtils.c --
 */
int tk_edid_to_editor(ClientData clientData, Tcl_Interp *interp, 
		       int argc, char *argv[]);
void tk_redisplaySequences(EdStruct *xx);
void ed_set_slider_pos(EdStruct *xx, int pos);
void ed_set_nslider_pos(EdStruct *xx, int pos);

/*
 * -- from edInterface.c --
 */
int edSetBriefSeqBase(EdStruct *xx, int x, int y, int mode);

/*
 * -- from edInterface.c --
 * Set ruler display mode - padded(0) vs unpadded(1).
 * ispadded == -1 implies toggle.
 */
int edSetRulerMode(EdStruct *xx, int ispadded);

/*
 * -- from edInterface.c --
 * Returns the active annotation types
 */
char *edGetActiveAnnos(EdStruct *xx);

/*
 * -- from edInterface.c --
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
int edUnpaddedBaseNumber(EdStruct *xx, int pos, int cache_count);

/*
 * -- from edInterface.c --
 * Formats tag information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 *
 * %%	Single % sign
 * %p	Tag position
 * %d	Tag direction (+/-/=, Raw: 0/1/2)
 * %D	Tag direction (----->/<-----/<---->, Raw: 0/1/2)
 * %t	Tag type (always 4 characters)
 * %l	Tag length
 * %n	Tag number (0 if unknown)
 * %c	Tag comment
 *
 * Additionally, some formats (p, l, n and c) can be specified as
 * %<number><format> (eg %80c) to allow AT MOST that many characters.
 */
int edSetBriefTag(EdStruct *xx, int seq, tagStruct *t, char *format);

/*
 * Find the current maximum sequence length in the editor.
 * If 'clipped' is true then we only look at the used quality-clipped portion,
 * otherwise we consider the hidden data too.
 */
int dbi_max_gel_len(DBInfo *db, int clipped);

/*
 * Redisplays the editor status lines (strands, translations)
 */
void tk_redisplaySeqStatusCompute(EdStruct *xx, int pos, int width);


void find_exons(EdStruct *xx, int pos, int width, int generate);

/*
 * Reports mutations occuring inside and outside of exons. When in an exon
 * we also report if it is a silent mutation.
 *
 * diffs_tagged is true when we only want to report mutations where HETE and
 * MUTA tags exist, otherwise base call differences are used instead.
 *
 * sort_by_position, when true, outputs mutations column-by-column instead of
 * sequence by sequence.
 *
 * Returns: HTML report of the mutations as a dstring pointer. This should
 *		be freed by the caller.
 *	    NULL if failed, and sets err_msg pointer to an error string.
 */
dstring_t *report_mutations(EdStruct *xx, int diffs_tagged,
			    int sort_by_position, char *dir, int detailed,
			    char **err_msg);

/* Returns the character underneath the editor cursor */
int edGetChar(EdStruct *xx);

#endif
