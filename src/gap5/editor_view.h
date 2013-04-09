#ifndef _EDITOR_VIEW_H_
#define _EDITOR_VIEW_H_

#include "tg_gio.h"
#include "tkSheet.h"
#include "tkEditor.h"
#include "tkEdNames.h"
#include "consensus.h"

#define MAX_DISPLAY_WIDTH 1000
#define MAX_NAME_WIDTH    256
#define WIN_NAME_SIZE 100

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
#define ED_DISP_XSCROLL		(1<<3)
#define ED_DISP_YSCROLL		(1<<9)
/* Sequence and consensus */
#define ED_DISP_SEQS		(ED_DISP_READS | ED_DISP_CONS)
/* Everything */
#define ED_DISP_ALL		(ED_DISP_SEQS | ED_DISP_NAMES | \
				 ED_DISP_XSCROLL | ED_DISP_YSCROLL | \
                                 ED_DISP_STATUS | \
				 ED_DISP_RULER | ED_DISP_CURSOR | \
				 ED_DISP_SELECTION | ED_DISP_HEIGHT)
/* Single items only - see xx->refresh_seq */
#define ED_DISP_NAME		(1<<10)
#define ED_DISP_SEQ		(1<<11)
/* For optimising redisplayDisagreements to remove flicker */
#define ED_DISP_NO_DIFFS	(1<<12)

/* Flags not implemented yet */
#define DB_Flags(xx, seq) (0)
#define DB_setFlags(xx, seq, flags) (0)

struct _EdLink;
enum States {StateDown=0,StateUp};

typedef struct _edview {
    /* A derived IO struct */
    GapIO *io;
    tg_rec cnum;

    /* Necessary Tcl/Tk bits and bobs */
    Tcl_Interp *interp;
    int editor_id;
    char edname[20];
    char seq_win[WIN_NAME_SIZE];
    char name_win[WIN_NAME_SIZE];
    Editor *ed;
    edNames *names;

    /* Display window size and location */
    int displayPos;
    int displayYPos;
    int displayWidth;
    int displayHeight;
    int names_xPos;

    void (*dispFunc)(void *, int, int, int, void *);
    enum States editorState;
    int refresh_flags;
    tg_rec refresh_seq;

    /* Cached data */
    char        displayedConsensus[MAX_DISPLAY_WIDTH];
    consensus_t cachedConsensus[MAX_DISPLAY_WIDTH];

    /* Cursor coordinates, sequence number and position in seq */
    int    cursor_type;
    tg_rec cursor_rec;
    int    cursor_pos;
    int    cursor_apos; /* absolute position in contig */
    cursor_t *cursor;
    int reg_id; /* registration id */

    /* Y coordinates for the elements in the editor */
    int y_seq_start;
    int y_seq_end;
    int y_cons;
    int y_numbers;

    /* Join editor bits */
    struct _EdLink *link;

    int trace_lock;

    /* Cached range query & Y position */
    rangec_t *r;
    int nr;
    int max_height;
    int r_start, r_end;
    /* FIXME: add cached index into r[] foreach row[y], as it's sorted on y */

    /* Maps r[i].anno.obj_rec to i */
    HacheTable *anno_hash;

    /* Maps record numbers to r[i] indices */
    HacheTable *rec_hash;

    /* Selection */
    int    select_made;
    tg_rec select_seq;
    int    select_start;
    int    select_end;

    /* Traces recently loaded */
    HacheTable *trace_hash;
    
    /* Seq sort settings */
    seq_sort_t sort_settings;
} edview;

typedef struct _EdLink {
    edview *xx[2];		/* A pair of linked editors */
    int locked;			/* Set if these are locked together */
    int lockOffset;		/* And their relative distances */
    tkSheet *diffs;
} EdLink;

/*
 * A C interface to the edit_contig and join_contig Tcl functions.
 */
int edit_contig(GapIO *io, tg_rec cnum, tg_rec rnum, int pos);
int join_contig(GapIO *io, tg_rec cnum[2], tg_rec rnum[2], int pos[2]);

/*
 * Allocates and initialises a new edview
 */
edview *edview_new(GapIO *io, tg_rec contig, tg_rec crec, int cpos,
		   Editor *ed, edNames *names,
		   void (*dispFunc)(void *, int, int, int, void *),
		   Tcl_Interp *interp);
/*
 * Deallocates an edview
 */
void edview_destroy(edview *xx);

/*
 * Changes the contig number in the edview struct.
 */
void edview_renumber(edview *xx, tg_rec new_contig);

/*
 * Finds an existing editor widget for a specific contig.
 * Returns edview on success
 *         NULL on failure
 */
edview *edview_find(GapIO *io, tg_rec contig);

/*
 * The main editor redraw function
 */
int edview_redraw(edview *xx);

/* Update X scrollbar of names display */
void ed_set_nslider_pos(edview *xx, int pos);
int set_displayPos(edview *xx, int pos);

/*
 * Formats tag information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 */
char *edGetBriefTag(edview *xx, tg_rec anno_ele, char *format);

/*
 * Formats reading information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 */
char *edGetBriefSeq(edview *xx, tg_rec seq, int pos, char *format);

/*
 * Formats consensus information for the status line.
 * This is done using a format string where certain % rules are replaced by
 * appropriate components.
 */
char *edGetBriefCon(edview *xx, tg_rec crec, int pos, char *format);

/*
 * Given an X,Y coordinate return the reading id under this position.
 *
 * Returns record number on success
 *         -1 on failure.
 */
int edGetGelNumber(edview *xx, int x, int y);

/*
 * Identifies the type of object underneath a specific row and column.
 * 'name' is a boolean which when true indicates the row,col are in the
 * names panel instead of the sequence panel.
 * 'seq_only' forces the item to be a sequence or consensus, and not
 * an object on them (eg annotation).
 *
 * Returns the item type GT_* on success and the record/pos in *rec, *pos
 *         -1 on failure (eg numbers, off screen, etc)
 */
int edview_item_at_pos(edview *xx, int row, int col, int name, int exact,
		       int seq_only, tg_rec *rec, int *pos);

/* Cursor movement control */
void edSetApos(edview *xx);
int edSetCursorPos(edview *xx, int type, tg_rec rec, int pos, int visible);
int edCursorUp(edview *xx);
int edCursorDown(edview *xx);
int edCursorLeft(edview *xx);
int edCursorRight(edview *xx);
int edReadStart(edview *xx);
int edReadStart2(edview *xx);
int edReadEnd(edview *xx);
int edReadEnd2(edview *xx);
int edContigStart(edview *xx);
int edContigEnd(edview *xx);

/*
 * Force the cursor to be visible. If x_safe or y_safe are true then
 * we omit some of the searching and assume there is no reason to check
 * that x or y is still visible.
 *
 * Returns 1 if redraw has taken place
 *         0 if not
 */
int showCursor(edview *xx, int x_safe, int y_safe);


int inJoinMode(edview *xx);
void edDisplayTrace(edview *xx);


/*
 * Given a sequence record number this identifies all other sequence
 * records from the same template. The returned array is malloced and should
 * be freed by the caller once finished with.
 *
 * Returns pointer to array of records of size *nrec on success
 *         NULL on failure (or zero found)
 */
tg_rec *edGetTemplateReads(edview *xx, tg_rec seqrec, int *nrec);

/*
 * Handles the align button in the join editor
 * Returns 0 for success
 *        -1 for failure
 */
int edJoinAlign(edview *xx, int fixed_left, int fixed_right);

/*
 * Perform the actual join process
 * Returns 0 for success
 *        -1 for failure
 */
int edJoin(edview *xx);

/*
 * Populates the cache of visible items in xx->r and xx->nr.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int edview_visible_items(edview *xx, int start, int end);

/*
 * X11 Selection handling
 */
int edSelectClear(edview *xx);
void edSelectFrom(edview *xx, int pos);
void edSelectTo(edview *xx, int pos);
void edSelectSet(edview *xx, tg_rec rec, int start, int end);

/*
 * Searching - see edview_search.c
 */
int edview_search(edview *xx, int forwards, int strand,
		  char *type, char *value);

/*
 * Convert from a record and position to a window X,Y coordinate in sheet
 * units.
 *
 * Returns 0 on success and stores via x and y pointers.
 *        -1 on failure (rec/pos not visible).
 */
int edGetXY(edview *xx, int rec_type, tg_rec rec, int pos, int *x, int *y);


/* In editor_join.c */
int edJoinMismatch(edview *xx, int *len, int *mismatch);

/*
 * Finds the next/previous difference between a pair of join editors.
 *
 * Returns 0 on sucess
 *        -1 on failure
 */
int edNextDifference(edview *xx);
int edPrevDifference(edview *xx);

/* Set the group_by/sort settings */
void edview_set_sort_order(edview *xx);
void ed_set_base_sort_point(edview *xx);
void ed_set_sequence_sort(edview *xx);

/* Compute original positions array via alignments */
int origpos(edview *xx, tg_rec srec, int pos);

/* Find records between from_rec and to_rec in the names layout. */
Array edview_items_between(edview *xx, tg_rec from_rec, tg_rec to_rec);

#endif /* _EDITOR_VIEW_H_ */
