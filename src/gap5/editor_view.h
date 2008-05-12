#ifndef _EDITOR_VIEW_H_
#define _EDITOR_VIEW_H_

#include "tg_gio.h"
#include "tkSheet.h"
#include "tkEditor.h"
#include "tkEdNames.h"

#define MAX_DISPLAY_WIDTH 1000
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

/* Orig positions not implemented yet */
#define origpos(xx, seq, pos) ((pos)+1)

/* Flags not implemented yet */
#define DB_Flags(xx, seq) (0)
#define DB_setFlags(xx, seq, flags) (0)

struct _EdLink;
enum States {StateDown=0,StateUp};

typedef struct _edview {
    /* A derived IO struct */
    GapIO *io;
    int cnum;
    contig_t *contig;

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
    int refresh_seq;

    /* Cached data */
    char displayedConsensus[MAX_DISPLAY_WIDTH];

    /* Cursor coordinates, sequence number and position in seq */
    int cursor_type;
    int cursor_rec;
    int cursor_pos;
    int cursor_apos; /* absolute position in contig */

    /* Y coordinates for the elements in the editor */
    int y_seq_start;
    int y_seq_end;
    int y_cons;
    int y_numbers;

    /* Join editor bits */
    struct _EdLink *link;

    int trace_lock;
} edview;

typedef struct _EdLink {
    edview *xx[2];		/* A pair of linked editors */
    int locked;			/* Set if these are locked together */
    int lockOffset;		/* And their relative distances */
    tkSheet *diffs;
} EdLink;

/*
 * Allocates and initialises a new edview
 */
edview *edview_new(GapIO *io, int contig,
		   Editor *ed, edNames *names,
		   void (*dispFunc)(void *, int, int, int, void *));

/*
 * Deallocates an edview
 */
int edview_redraw(edview *xx);

/* Update X scrollbar of names display */
void ed_set_nslider_pos(edview *xx, int pos);
int set_displayPos(edview *xx, int pos);

/*
 * Formats reading information for the status line. This is done using a format
 * string where certain % rules are replaced by appropriate components.
 */
char *edGetBriefSeq(edview *xx, int seq, int pos, char *format);

/*
 * Given an X,Y coordinate return the reading id under this position.
 *
 * Returns record number on success
 *         -1 on failure.
 */
int edGetGelNumber(edview *xx, int x, int y);

/*
 * Identifies the type of object underneath a specific row and column.
 *
 * Returns the item type GT_* on success and the record/pos in *rec, *pos
 *         -1 on failure (eg numbers, off screen, etc)
 */
int edview_item_at_pos(edview *xx, int row, int col, int *rec, int *pos);

/* Cursor movement control */
int edSetCursorPos(edview *xx, int type, int rec, int pos);
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

int inJoinMode(edview *xx);
void edDisplayTrace(edview *xx);


#endif /* _EDITOR_VIEW_H_ */
