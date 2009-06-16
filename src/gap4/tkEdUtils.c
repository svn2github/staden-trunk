#include <ctype.h>
#include <stdlib.h>
#include <X11/Xatom.h> /* XA_PRIMARY - included in Tk distribrution */
#include <math.h>
#include <time.h>

#include "gap_globals.h"
#include "edUtils.h"
#include "edStructs.h"
#include "tkEditor.h"
#include "tkEdNames.h"
#include "tkSheet.h"
#include "genetic_code.h"
#include "tcl_utils.h"
#include "misc.h"
#include "dna_utils.h"
#include "notedb.h"
#include "tman_interface.h"
#include "dstring.h"

#ifdef _WIN32
/* 6/1/99 johnt - defines various X calls and other stuff for Windows */
#  include "tkWinPort.h"
#endif

#define SHOW_QUAL

/*
 * Returns 'type' or 0 to indicate whether position 'pos' of sequence 'seq' has
 * been edited. We are scanning in direction 'dir' and 'comp' tells us whether
 * the sequence has been complimented or not.
 *
 * type == 1 for base deleted
 * type == 2 for base changed/inserted
 * type == 3 for pad changed/inserted
 * type == 4 for confidence
 *
 * Call will 'xx == NULL' to initialise.
 */
static int is_edit(EdStruct *xx, int seq, int pos, int dir, int comp) {
    static int last, last_pos;
    int c1, o1, o2, tmpl = last, tmpp = last_pos;
    char b1;
    
    if (!xx) {
	last = 0;
	last_pos = 0;
	return 0;
    }

    /*
     * At present the opos array is 16 bit, so we do not report edits in
     * longer sequences.
     */
    if (ABS(DB_Length2(xx, seq)) >= 32768)
	return 0;

    if (pos >= 0 && pos < DB_Length2(xx, seq)) {
	c1 = DB_Conf(xx, seq)[pos];

	/* remember last non zero original position */
	o1 = DB_Opos(xx, seq)[pos];
	if (o1) {
	    last = o1;
	    last_pos = pos;
	}
	
	b1 = DB_Seq(xx, seq)[pos];
	
	if (pos + dir >= 0 && pos + dir < DB_Length2(xx, seq)) {
	    o2 = DB_Opos(xx, seq)[pos+dir];
    
	    /* Hmm! Checks for a deletion */
	    if (tmpp && (comp * ((last-tmpl) - comp*(tmpp-last_pos))) < 0) {
		return 1;
	    }
	}
	    
	/* treat pads as a different status to other bases */
	if (!o1)
	    return (b1 == '*') ? 3 /* pad */ : 2; /* base */
    }

    /* all clear */
    return 0;
}

static void redisplaySelection(EdStruct *xx);

/*
 *---------------------------------------------------------------------------
 * Functions to call directly from C.
 * These take EdStruct arguments and do appropriate translations to tcl
 * calls.
 *---------------------------------------------------------------------------
 */

void front_editor(EdStruct *xx) {
    char *path = Tk_PathName(EDTKWIN(xx->ed));
    if (!path)
	return;

    Tcl_VarEval(EDINTERP(xx->ed), "set x [winfo toplevel ", path,
		"];wm deiconify $x; raise $x", NULL);
}

void ed_set_slider_pos(EdStruct *xx, int pos) {
    Editor *ed = xx->ed;
    char buf[100];
    int len = xx->displayWidth;
    int total = xx->extent_right - xx->extent_left + 1;

    if (!ed || xx->editorState == StateDown)
	return;

    if (!ed->xScrollCmd)
	return;

    pos -= xx->extent_left;

    /* total width start end */
    sprintf(buf, " %.20f %.20f",
	    pos / (double)total,
	    (pos + len) / (double)total);

    /*
     * Ignore errors as the widget may not exist (eg when shutting
     * down the editor). Hacky, but simple.
     */
    if (Tcl_VarEval(EDINTERP(ed), ed->xScrollCmd, buf, NULL)
	!= TCL_OK) {
	Tcl_AddErrorInfo(EDINTERP(ed),
			 "\n    (xscrollcommand executed by Editor)");
	Tcl_BackgroundError(EDINTERP(ed));
    }

    /* Editor 'mini-traces' only */
    tman_reposition_traces(xx, xx->displayPos + xx->displayWidth/2, 1);
}

static void shuffle_traces(EdStruct *xx) {
    char *edpath = Tk_PathName(EDTKWIN(xx->ed));
    int *seqList;
    int i;
    int pos = xx->displayPos;
    int width = xx->displayWidth;
    char cmd[1024];

    if (xx->lines_per_seq == 1)
	return;

    seqList = sequencesOnScreen(xx,pos, width);

    /* Initialise font height */
    if (!xx->fontHeight) {
	Tk_FontMetrics fm;
	Tk_Font fn = Tk_GetFont(EDINTERP(xx->ed), EDTKWIN(xx->ed),
				"sheet_font");

	if (fn == NULL) {
	    fprintf(stderr, "Font sheet_font not found\n");
	    return;
	}

	Tk_GetFontMetrics(fn, &fm);
	xx->fontHeight = fm.linespace;
    }

    /*
     * We need this as traces that were beneath the consensus (due to scrolling
     * from a deep to a shallow section of the contig) need to be hidden so
     * that they do not overlap the status lines.
     */
    sprintf(cmd, "foreach t [winfo children %s] {place forget $t}",
	    edpath);
    if (TCL_OK != Tcl_Eval(EDINTERP(xx->ed), cmd)) {
	puts(Tcl_GetStringResult(EDINTERP(xx->ed)));
    }
    
    for (i=0 ; i < xx->displayHeight - xx->lines_per_seq;
	 i+=xx->lines_per_seq ) {
	int k;

	/* Set 'k' to an index into the visible sequences */
	if (i != xx->displayHeight-1)
	    k = i/xx->lines_per_seq + xx->displayYPos;
	else
	    k = (xx->totalHeight-1)/xx->lines_per_seq;

	if (seqList[k] == 0)
	    break;

	sprintf(cmd, "place %s.trace_%d -y %d; raise %s.trace_%d",
		edpath, seqList[k], (i+1) * xx->fontHeight,
		edpath, seqList[k]);
	if (TCL_OK != Tcl_Eval(EDINTERP(xx->ed), cmd)) {
	    /*
	     * This can fail if a trace is not visible because it has not
	     * been loaded. Currently we take the easy solution to this
	     * by simply recomputing the visible traces.
	     */
	    if (xx->lines_per_seq > 1) {
		/* So load it and try again... */
		showTrace(xx, seqList[k],
			  xx->displayPos + xx->displayWidth/2 -
			  (DB_RelPos(xx, seqList[k])-1),
			  xx->fontWidth, 0, 1);
		Tcl_Eval(EDINTERP(xx->ed), cmd);
	    }
	}
    }
}

void ed_set_yslider_pos(EdStruct *xx, int pos) {
    Editor *ed = xx->ed;
    char buf[100];

    if (!ed || xx->editorState == StateDown)
	return;

    if (!ed->yScrollCmd)
	return;

    shuffle_traces(xx); 

    pos -= xx->extent_left - 1;

    /* total width start end */
    sprintf(buf, " %.20f %.20f",
	    (xx->displayYPos*xx->lines_per_seq) / (double)xx->totalHeight,
	    (xx->displayYPos*xx->lines_per_seq + xx->displayHeight) /
	    (double)xx->totalHeight);

    /*
     * Ignore errors as the widget may not exist (eg when shutting
     * down the editor). Hacky, but simple.
     */
    if (Tcl_VarEval(EDINTERP(ed), ed->yScrollCmd, buf, NULL)
	!= TCL_OK) {
	Tcl_AddErrorInfo(EDINTERP(ed),
			 "\n    (yscrollcommand executed by Editor)");
	Tcl_BackgroundError(EDINTERP(ed));
    }
}

/* Update X scrollbar of names display */
void ed_set_nslider_pos(EdStruct *xx, int pos) {
    edNames *en = xx->names;
    char buf[1024];

    if (!en || xx->editorState == StateDown)
	return;

    if (en->xScrollCmd) {
	double fract1, fract2;
	fract1 = pos / (double)DB_NAMELEN;
	fract2 = (pos + en->sw.columns - (DB_GELNOLEN + 2)) /
	    (double)DB_NAMELEN;
	sprintf(buf, " %.20f %.20f", fract1, fract2);
	if (Tcl_VarEval(EDINTERP(en), en->xScrollCmd, buf, NULL) != TCL_OK) {
	    printf("Error in editor names scroll: %s\n", Tcl_GetStringResult(EDINTERP(en)));
	}
    }
}

void positionCursor(EdStruct *xx, int seq, int pos) {
    if (onScreen(xx, seq, pos, NULL)) {
        int screenRow,screenColumn;
        int *seqList;
        
        screenColumn = positionInContig(xx,seq,pos) - xx->displayPos;
        seqList = sequencesOnScreen(xx,xx->displayPos, xx->displayWidth);

	if (seq == 0) {
	    /* Consensus */
	    screenRow = xx->displayHeight - 1;
	} else {
	    /* Sequence */
	    for(screenRow=xx->displayYPos;
		screenRow < xx->displayHeight/xx->lines_per_seq +
		    xx->displayYPos && seqList[screenRow] != seq;
		screenRow++)
		;
	    if (seqList[screenRow] != seq) {
		XawSheetDisplayCursor(&xx->ed->sw,False);
		return;
	    }

	    screenRow = (screenRow - xx->displayYPos) * xx->lines_per_seq +
		(xx->lines_per_seq-1);
	}

        XawSheetDisplayCursor(&xx->ed->sw, True);
        XawSheetPositionCursor(&xx->ed->sw, screenColumn,
			       screenRow + xx->rulerDisplayed);
    } else {
        XawSheetDisplayCursor(&xx->ed->sw,False);
    }
}

#if 0
#define COL_TO_PIXEL(W,C) \
    ((W)->font->max_bounds.width * (C) + (W)->border_width)
#define ROW_TO_PIXEL(W,R) \
    (fontHeight((W)->font) * ((R)+1) + (W)->border_width)
#define GET_ARRAY_CELL(A,R,C)\
    ( &A->base[(R * A->cols + C)*A->size] )

void sheet_copy(Sheet *sw, int srcx, int srcy, int width, int height,
		int dstx, int dsty) {
    int r, c;
    sheet_paper paper_base = (sheet_paper *)sw->paper->base;
    sheet_ink ink_base = (sheet_ink *)sw->ink->base;

    /* Copy raw X display */

    XCopyArea(sw->display, sw->window, sw->window, sw->sparegc,
	      (int)COL_TO_PIXEL(sw, srcx),  (int)ROW_TO_PIXEL(sw, srcy),
	      (int)COL_TO_PIXEL(sw, width), (int)ROW_TO_PIXEL(sw, height),
	      (int)COL_TO_PIXEL(sw, dstx),  (int)ROW_TO_PIXEL(sw, dsty));

    /* Copy paper and ink details */
    if (srcx > dstx) {
	memcpy(paper_base, &paper_base[srcx-dstx],
	       sizeof(*paper_base) * (sw->columns * (height+1) - (srcx-dstx)));
	memcpy(ink_base, &ink_base[srcx-dstx],
	       sizeof(*ink_base) * (sw->columns * (height+1) - (srcx-dstx)));
    } else {
	memcpy(&paper_base[dstx-srcx], paper_base,
	       sizeof(*paper_base) * (sw->columns * (height+1) - (dstx-srcx)));
	memcpy(&ink_base[dstx-srcx], ink_base,
	       sizeof(*ink_base) * (sw->columns * (height+1) - (dstx-srcx)));
    }
#if 0
    if (srcx < dstx) {
	for (r = srcy+1; r < srcy+1+height; r++) {
	    sheet_ink ink_base;
	    sheet_paper paper_base;

	    ink_base = (sheet_ink)GET_ARRAY_CELL(sw->ink, r, srcx);
	    paper_base = (sheet_paper)GET_ARRAY_CELL(sw->paper, r, srcx);

	    for (c = srcx+width-1; c >= srcx; c--, ink_base--, paper_base--) {
		ink_base[c] = ink_base[c-1];
		paper_base[c] = paper_base[c-1];
	    }
	}
    } else {
	for (r = srcy+1; r < srcy+1+height; r++) {
	    sheet_ink ink_base;
	    sheet_paper paper_base;

	    ink_base = (sheet_ink)GET_ARRAY_CELL(sw->ink, r, dstx-1);
	    paper_base = (sheet_paper)GET_ARRAY_CELL(sw->paper, r, dstx-1);

	    for (c = srcx; c < srcx+width; c++, ink_base++, paper_base++) {
		ink_base[c] = ink_base[c+1];
		paper_base[c] = paper_base[c+1];
	    }
	}
    }
#endif
}

void tk_editor_scroll(EdStruct *xx, int orig_pos, int new_pos) {
    if (orig_pos < new_pos) {
	sheet_copy(&(TKSHEET(xx->ed)->sw),
		   new_pos - orig_pos, -1,
		   xx->displayWidth - (new_pos - orig_pos), xx->displayHeight,
		   0, -1);
    } else {
	sheet_copy(&(TKSHEET(xx->ed)->sw),
		   0, -1,
		   xx->displayWidth - (orig_pos - new_pos), xx->displayHeight,
		   orig_pos - new_pos , -1);
    }
}
#endif

/*
 * Calculates the numbers for the contig editor ruler line.
 * The return value is the index into this ruler buffer to plot.
 */
static int generate_ruler(EdStruct *xx, char *ruler, XawSheetInk *ink,
			  int pos, int width) {
    char *k = ruler;
    int j;

    int padded_pos[MAX_DISPLAY_WIDTH+21];

    memset(ruler, ' ', MAX_DISPLAY_WIDTH+21);
    /*
     * for (j = 0; j < MAX_DISPLAY_WIDTH+21; j++)
     *     ink[j].sh = sh_default;
     */

    if (DBI(xx)->reference_seq) {
	/* Number relative to a specific sequence number */
	char *bases = DBgetSeq(DBI(xx), DBI(xx)->reference_seq);
	int reflen = DB_Length(xx, DBI(xx)->reference_seq);
	int rp = DB_RelPos(xx, DBI(xx)->reference_seq);

	/* 
	 * Compute unpadded base positions for this window.
	 * WARNING: This code is full of warts. I carefully worked it out, but
	 * it was not quite perfect. Working out again gave a slightly
	 * different, but equally imperfect answer. So I applied the next
	 * logical step - I hope you like fudge!
	 *
	 * If the change the code - check for these cases and combinations:
	 *   Position from left end when uncomplemented
	 *   Position from right end when complemented
	 *   Negative offsets
	 *   When relpos of ref is 1 and when it is not, both when refseq
	 *     is complemented and when it is not.
	 *   Effect of pads in ref seq
	 *   That base numbers are never positioned above the pads
	 *   Circular sequences (both orientations)
	 *   Offset base numbers
	 */
	if ((DB_Comp(xx, DBI(xx)->reference_seq) == UNCOMPLEMENTED)) {
	    int unpadded = MIN(0, pos);

	    for (j = unpadded; j < pos + width + 9; j++) {
		if (j >= pos)
		    padded_pos[j-pos] = unpadded-rp;
		if (j-rp < 0 || j-rp >= reflen || bases[j-rp] != '*')
		    unpadded++;
	    }
	} else {
	    int unpadded = reflen - MAX(pos + width + 8, reflen);

	    for (j = MAX(pos + width + 8, reflen); j >= pos; j--) {
		if (j <= pos + width + 8)
		    padded_pos[j-pos] = unpadded+(rp-1);
		if (j-rp-1 < 0 || j-rp-1 >= reflen || bases[j-rp-1] != '*')
		    unpadded++;
	    }
	}

	for (j = pos; j < pos + width + 9; j++, k++, ink++) {
	    int unpadded2 = (padded_pos[j-pos] + DBI(xx)->reference_offset);

	    if (DBI(xx)->reference_len) {
		unpadded2 %= DBI(xx) -> reference_len;
		while (unpadded2 < 0) {
		    unpadded2 += DBI(xx)->reference_len;
		}
		if (unpadded2 == 0)
		    continue; /* Don't display 0 for circular seqs */
	    }
	    if (!(unpadded2 % 10)) {
		sprintf(k, "%10d", unpadded2);
	    }
	}

	/*
	 * Pads in the consensus can sometimes leave nulls in the buffer,
	 * which turn into square boxes on some fonts. Replace these with
	 * spaces.
	 */
	for (j = 0; j < MAX_DISPLAY_WIDTH+21; j++)
	    if (ruler[j] == '\0')
		ruler[j] = ' ';

	return 9;

    } else if (!xx->unpadded_ruler) {
	/* Basic numbering */
	int lower,times;
	lower = (pos - pos%10);
	times = width/10 + 3;
	for (j=0;j<times;j++,k+=10,lower+=10)
	    sprintf(k,"%10d",lower);
	return 9+pos%10;
	
    } /* else */ {
	/* Number by unpadded base calls, using the consensus */
	int unpadded;
	
	edUnpaddedBaseNumber(xx, pos, width + 9);
	for (j = pos; j < pos + width + 9; j++, k++) {
	    unpadded = edUnpaddedBaseNumber(xx, j, 0);
	    if (unpadded%10)
		continue;
	    sprintf(k, "%10d", unpadded);
	    if (k+10-ruler < MAX_DISPLAY_WIDTH+21)
		k[10]=' ';
	}
	edUnpaddedBaseNumber(xx, pos, -1);

	return 9;
    }
}

/*
 * Deals with redisplay after editor scroll and the associated height
 * changes that may occur. This may change xx->refresh_flags.
 */
static void tk_redisplaySeqScroll(EdStruct *xx, int status_depth_changed) {
    int pos = xx->displayPos;
    int width = xx->displayWidth;
    int i, h_changed;

    if (xx->refresh_flags & ED_DISP_SCROLL) {
	/*puts("set slider pos");*/
	ed_set_slider_pos(xx, pos);
	ed_set_nslider_pos(xx, xx->names_xpos);
    }
    
    /*
     * Set Up Text Window sizes
     */
    i = xx->totalHeight;
    xx->totalHeight = (linesOnScreen(xx,pos,width)-xx->consensusDisplayed)
	*xx->lines_per_seq + xx->consensusDisplayed;
    h_changed = (xx->totalHeight != i) || status_depth_changed;
    xx->displayHeight = xx->totalHeight > xx->ed->max_height
	? xx->ed->max_height : xx->totalHeight;

    if (xx->displayYPos > (xx->totalHeight - xx->displayHeight +
			   xx->lines_per_seq-1) / xx->lines_per_seq) {
	xx->displayYPos = (xx->totalHeight - xx->displayHeight +
			   xx->lines_per_seq-1) / xx->lines_per_seq;
	xx->refresh_flags |= ED_DISP_YSCROLL;
    }

    if (xx->refresh_flags & ED_DISP_YSCROLL) {
	ed_set_yslider_pos(xx, xx->displayYPos);

	/* Vertical scroll implies change of readings and names only */
	xx->refresh_flags |= ED_DISP_READS | ED_DISP_NAMES;
    }

    if (xx->refresh_flags & ED_DISP_SCROLL) {
	/* It's simplest to redraw everything when a scroll
	 * occurrs. We cannot rely on checking the height is the
	 * same and simply copying areas of screen/memory around
	 * as it's possible to add one sequence and remove one
	 * sequence to keep the same height.
	 */
	xx->refresh_flags |= ED_DISP_ALL;
    }
    xx->refresh_pos = xx->displayPos;

    if ((xx->refresh_flags & ED_DISP_SCROLL && h_changed) ||
	(xx->refresh_flags & ED_DISP_HEIGHT)) {
	TKSHEET(xx->ed)->divider=0;
	TKSHEET(xx->names)->divider=0;
	sheet_set_display_height(TKSHEET(xx->ed),
				 xx->displayHeight + 1 + xx->status_depth);
	sheet_set_display_height(TKSHEET(xx->names),
				 xx->displayHeight + xx->status_depth);
    }

    xx->refresh_flags &= ~ED_DISP_SCROLL;
}

static void note_ink(XawSheetInk *ink, int db, int sense) {
    if (sense == 0) {
	if (note_db[db].fg_colour) {
	    ink->sh |= sh_fg;
	    ink->fg = note_db[db].fg_pixel;
	}
	if (note_db[db].bg_colour) {
	    ink->sh |= sh_bg;
	    ink->bg = note_db[db].bg_pixel;
	}
    } else {
	if (note_db[db].gf_colour) {
	    ink->sh |= sh_fg;
	    ink->fg = note_db[db].gf_pixel;
	}
	if (note_db[db].gb_colour) {
	    ink->sh |= sh_bg;
	    ink->bg = note_db[db].gb_pixel;
	}
    }
}

/*
 * Returns the a character and colour indicating which, if any, notes
 * are present for this sequence (or contig).
 */
static char tk_redisplaySeqNotes(EdStruct *xx, int seq, XawSheetInk *ink) {
    char c = ' ';
    int db;

    if (seq) {
	if (DB_Flags(xx, seq) & DB_FLAG_REFSEQ) {
	    c = 'S';
	    db = note_id2index("REFS");
	    note_ink(ink, db == -1 ? 0 : db, 0);
	}
	if (DB_Flags(xx, seq) & DB_FLAG_REFTRACE_POS) {
	    c = (DB_Comp(xx, seq) == UNCOMPLEMENTED) ? 'f' : 'r';
	    db = note_id2index("REFT");
	    note_ink(ink, db == -1 ? 0 : db, 0);
	}
	if (DB_Flags(xx, seq) & DB_FLAG_REFTRACE_NEG) {
	    c = (DB_Comp(xx, seq) == UNCOMPLEMENTED) ? 'F' : 'R';
	    db = note_id2index("REFT");
	    note_ink(ink, db == -1 ? 0 : db, 1);
	}

	if (c == ' ') {
	    GReadings r;
	    gel_read(DBI_io(xx), DB_Number(xx, seq), r);
	    if (r.notes)
		c = '*';
	}

    } else {
	GContigs co;
	contig_read(DBI_io(xx), DBI_contigNum(xx), co);
	if (co.notes)
	    c = '*';
    }

    return c;
}

/*
 * Handles redisplaying of the sequence names
 */
static void tk_redisplaySeqNames(EdStruct *xx, int *seqList) {
    XawSheetInk nsplodge[NAMELEN];
    int i, j, search;
    int lastset = 0, set;

    /* all names, or just xx->refresh_seq ? */
    search = (xx->refresh_flags & ED_DISP_NAMES) ? 0 : 1;

    for (i = xx->lines_per_seq-1; i < xx->displayHeight + xx->lines_per_seq-1;
	 i+=xx->lines_per_seq) {
	int k;
	char *name, buf[NAMELEN+2];

	if (i < xx->displayHeight-1)
	    k = i/xx->lines_per_seq + xx->displayYPos;
	else
	    k = (xx->totalHeight-1)/xx->lines_per_seq;

	if (seqList[k] == 0) {
	    if (i < xx->displayHeight-1) {
		/* blank the gap */
		memset(buf, ' ', NAMELEN+1);
		XawSheetPutText(&xx->names->sw, 0, i, NAMELEN+1, buf);
	    }
	    i = xx->displayHeight-1;
	}

	for (j=0 ; j < NAMELEN ; j++) {
	    nsplodge[j].sh = sh_default;
	    nsplodge[j].fg = 0; /* not used, but prevents warnings */
	    nsplodge[j].bg = 0; /* not used, but prevents warnings */
	}

	/* "Set" colouring */
	set = xx->set ? xx->set[seqList[k]] : 0;
	*buf = ' ';
	if (seqList[k] != 0) {
	    for (j=DB_GELNOLEN+2 ; j < NAMELEN ; j++) {
		nsplodge[j].sh |= sh_bg;
		nsplodge[j].bg = xx->set_bg[set % 10];
	    }
#if 0
	    if (set != lastset) {
		*buf = "-+"[xx->set_collapsed && xx->set_collapsed[set]];
	    }
#endif
	    lastset = set;
	}

	if (search && seqList[k] != xx->refresh_seq)
	    continue;

	if (DB_Flags(xx, seqList[k]) & DB_FLAG_SELECTED) {
	    for (j=DB_GELNOLEN+2 ; j < NAMELEN ; j++)
		nsplodge[j].sh |= sh_inverse;
	}

	if (DB_Flags(xx, seqList[k]) & DB_FLAG_TRACE_SHOWN) {
	    for (j=0 ; j < NAMELEN ; j++) {
		nsplodge[j].sh |= sh_light;
	    }
	}

	if (DB_Flags(xx, seqList[k]) & DB_FLAG_INVIS) {
	    for (j=1; j <= DB_GELNOLEN; j++) {
		nsplodge[j].sh |= sh_bg;
		nsplodge[j].bg = xx->qual_bg[0];
	    }
	}

	/*
	 * Colour the base between number and name if the template for this
	 * sequence is inconsistent.
	 */
	{
	    int tnum = DBI(xx)->DB[seqList[k]].template;
	    if (DBI(xx)->templates[tnum]) {
		if (DBI(xx)->templates[tnum]->consistency) {
		    nsplodge[DB_GELNOLEN+1].sh |= sh_bg;
		    switch (DBI(xx)->templates[tnum]->consistency) {
		    case TEMP_CONSIST_DIST:
			nsplodge[DB_GELNOLEN+1].bg = xx->tmpl_bg[0];
			break;
		    case TEMP_CONSIST_STRAND:
			nsplodge[DB_GELNOLEN+1].bg = xx->tmpl_bg[1];
			break;
		    case TEMP_CONSIST_PRIMER:
			nsplodge[DB_GELNOLEN+1].bg = xx->tmpl_bg[2];
			break;
		    default:
			/* Multiple problems */
			nsplodge[DB_GELNOLEN+1].bg = xx->tmpl_bg[3];
			break;
		    }
		} else if (DBI(xx)->templates[tnum]->flags &
			   (TEMP_FLAG_GUESSED_END |
			    TEMP_FLAG_GUESSED_START)) {
		    nsplodge[DB_GELNOLEN+1].sh |= sh_bg;
		    if (DBI(xx)->templates[tnum]->flags & TEMP_FLAG_SPANNING)
			nsplodge[DB_GELNOLEN+1].bg = xx->tmpl_bg[5];
		    else
			nsplodge[DB_GELNOLEN+1].bg = xx->tmpl_bg[4];
		}
	    }
	}

	if (xx->template_names) {
	    name = DBgetTemplateName(DBI(xx), seqList[k]);
	} else {
	    name = DBgetName(DBI(xx), seqList[k]);
	}

	/* Indicate whether notes are present */
	buf[0] = tk_redisplaySeqNotes(xx, seqList[k],&nsplodge[0]);
	XawSheetPutJazzyText(&xx->names->sw, 0, i, 1, buf, nsplodge);

	/* Display the reading number and name */
	strncpy(buf, name, DB_GELNOLEN+1);
	strcpy(buf + DB_GELNOLEN+1, name + DB_GELNOLEN+1 + xx->names_xpos);

	if (seqList[k] == 0) {
	    strcpy(buf, name);
	    if (xx->displayYPos > 0) {
		buf[0] = '<';
	    }
	    if (xx->displayYPos*xx->lines_per_seq + xx->displayHeight
		< xx->totalHeight) {
		buf[1] = '>';
	    }
	}

	XawSheetPutJazzyText(&xx->names->sw, 1, i, NAMELEN, buf, nsplodge+1);
    }
}

/*
 * Redraws the emacs-style edit status (-**-, etc)
 */
static void tk_redisplaySeqEditStatus(EdStruct *xx) {
    char status1 = '-', status2 = '-';
    char buf[4];
    int consensus_line = xx->displayHeight - 1;

    if (!DBI_flags(xx) & DB_ACCESS_UPDATE)
	status1 = status2 = '%';
    if (DBI_edits_made(xx)) {
	status1 = '*';
	if (status2 == '-')
	    status2 = '*';
    }
    buf[0] = '-';
    buf[1] = status2;
    buf[2] = status1;
    buf[3] = '-';
    
    XawSheetPutText(&xx->names->sw, xx->names->sw.columns-4,
		    consensus_line, 4, buf);
}

/*
 * Colours sequence bases according to edit information
 */
static void tk_redisplaySeqEdits(EdStruct *xx, XawSheetInk *splodge,
				 int seq, int pos, int width) {
    int j, m, dir, type;
    int istart = 0, iend = width, cstart, cend;
		
    if (xx->reveal_cutoffs) {
	cstart = pos - DB_RelPos(xx, seq) +
	    DB_Start(xx, seq);
	cend = cstart + width;
		    
	if (cstart < 0) {
	    istart = -cstart;
	    cstart = 0;
	}
	if (cend-1 > DB_Length2(xx, seq)) {
	    iend -= (cend-1 - DB_Length2(xx, seq));
	}
		    
	m = cstart;
    } else {
	cstart = pos - DB_RelPos(xx, seq);
	cend = cstart + width;
		    
	if (cstart < 0) {
	    istart = -cstart;
	    cstart = 0;
	}
	if (cend-1 > DB_Length(xx, seq)) {
	    iend -= (cend-1 - DB_Length(xx, seq));
	}
		    
	m = cstart + DB_Start(xx, seq);
    }
		
    is_edit(0, 0, 0, 0, 0);
    dir = DB_Comp(xx, seq) == -1 ? 1 : -1;
    /* FIXME - slow */
    for (j = istart; j < iend; j++, m++) {
	if (type = is_edit(xx, seq, m, 1, dir)) {
	    splodge[j].sh |= sh_bg;
	    splodge[j].bg = xx->edit_bg[type-1];
	}
    }
}

/*
 * Colouration for highlight disagreements
 */
static void tk_redisplaySeqDiffs(EdStruct *xx, XawSheetInk *splodge,
				 char *seq, int pos, int width,
				 int seqn) {
    int j;
    int1 *conf;
    int cstart;

    /* Get confidence */
    cstart = pos - DB_RelPos(xx, seqn) + DB_Start(xx, seqn);
    conf = DB_Conf(xx, seqn);

    if ((xx->showDifferences&3) == 1) {
	/* By dots */
	if (xx->showDifferences&4) {
	    int j;
	    for (j=0; j<width; j++)
		if (seq[j]==xx->displayedConsensus[j])
		    seq[j]='.';
		else if (seq[j] != ' ' && conf[cstart + j] < xx->diff_qual)
		    seq[j]=':';
	} else {
	    int j;
	    for (j=0; j<width; j++)
		if (tolower(seq[j]) == tolower(xx->displayedConsensus[j]))
		    seq[j]='.';
		else if (seq[j] != ' ' && conf[cstart + j] < xx->diff_qual)
		    seq[j]=':';
	}

    } else if (xx->showDifferences&4) {
	/* By background */
	for (j = 0; j < width; j++) {
	    if (seq[j] != xx->displayedConsensus[j] &&
		seq[j] != ' ' &&
		conf[cstart + j] >= xx->diff_qual) {
		if ((xx->showDifferences&3) == 3) {
		    splodge[j].sh |= sh_bg;
		    splodge[j].bg = xx->diff_bg;
		} else {
		    splodge[j].sh |= sh_fg;
		    splodge[j].fg = xx->diff_bg;
		}
	    }
	}

    } else {
	/* By foreground */
	for (j = 0; j < width; j++) {
	    if (tolower(seq[j]) !=
		tolower(xx->displayedConsensus[j]) &&
		seq[j] != ' ' &&
		conf[cstart + j] >= xx->diff_qual) {
		if ((xx->showDifferences&3) == 3) {
		    splodge[j].sh |= sh_bg;
		    splodge[j].bg = xx->diff_bg;
		} else {
		    splodge[j].sh |= sh_fg;
		    splodge[j].fg = xx->diff_bg;
		}
	    }
	}
    }
}

/*
 * Colours sequences according to their confidence values
 */
static void tk_redisplaySeqReadQual(EdStruct *xx, XawSheetInk *splodge,
				    int seq, int pos, int width) {
    int q = xx->qual_cut;
    int istart = 0, iend = width, cstart, cend, j, m;
    int1 *conf;
	    
    if (xx->reveal_cutoffs) {
	cstart = pos - DB_RelPos(xx, seq) + DB_Start(xx, seq);
	cend = cstart + width;
		
	if (cstart < 0) {
	    istart = -cstart;
	    cstart = 0;
	}
	if (cend > DB_Length2(xx, seq)) {
	    iend -= (cend - DB_Length2(xx, seq));
	}
		
	conf =  &DB_Conf(xx, seq)[cstart];
    } else {
	cstart = pos - DB_RelPos(xx, seq);
	cend = cstart + width;
		
	if (cstart < 0) {
	    istart = -cstart;
	    cstart = 0;
	}
	if (cend > DB_Length(xx, seq)) {
	    iend -= (cend - DB_Length(xx, seq));
	}
		
	conf =  &DB_Conf(xx, seq)[cstart + DB_Start(xx, seq)];
    }
	    
    for (m = 0, j = istart; j < iend; j++, m++) {
	if (conf[m] < q || (conf[m] == 0 && q != -1)) {
	    if (xx->show_qual) {
		if (!(splodge[j].sh & sh_bg))
		    splodge[j].bg = xx->qual_bg[0];
		splodge[j].sh |= sh_default | sh_bg | sh_fg;
	    } else {
		splodge[j].sh |= sh_default | sh_fg;
	    }
	    splodge[j].fg = xx->qual_below;
	} else {
	    int col;
		    
	    if (!xx->show_qual || splodge[j].sh & sh_bg)
		continue;
		    
	    col = (10 * (conf[m] - q + 1)) / (100 - q + 1);
	    if (col >= 10)
		col = 9;
	    if (col < 0)
		col = 0;
		    
	    splodge[j].sh |= sh_default | sh_bg;
	    splodge[j].bg = xx->qual_bg[col];
	}
    }
}

/*
 * Colours "set" consensus according to their confidence values
 */
static void tk_redisplaySetQual(EdStruct *xx, XawSheetInk *splodge,
				float *conf, int width) {
    int q = xx->qual_cut;
    int i;
	    
    for (i = 0; i < width; i++) {
	if (conf[i] < q || (conf[i] == 0 && q != -1)) {
	    if (xx->show_qual) {
		if (!(splodge[i].sh & sh_bg))
		    splodge[i].bg = xx->qual_bg[0];
		splodge[i].sh |= sh_default | sh_bg | sh_fg;
	    } else {
		splodge[i].sh |= sh_default | sh_fg;
	    }
	    splodge[i].fg = xx->qual_below;
	} else {
	    int col;
		    
	    if (!xx->show_qual || splodge[i].sh & sh_bg)
		continue;
		    
	    col = (10 * (conf[i] - q + 1)) / (100 - q + 1);
	    if (col >= 10)
		col = 9;
	    if (col < 0)
		col = 0;
		    
	    splodge[i].sh |= sh_default | sh_bg;
	    splodge[i].bg = xx->qual_bg[col];
	}
    }
}

/*
 * Colours the consensus according to the confidence values.
 */
static void tk_redisplaySeqConsQual(EdStruct *xx, XawSheetInk *splodge,
				    int seq, int pos, int width) {
    int j, col, q = xx->qual_cut;
	    
    for (j = 0; j < width; j++) {
	if (xx->displayedConfidence[j] < q ||
	    xx->displayedConfidence[j] == 0) {
	    if (!(splodge[j].sh & sh_bg))
		splodge[j].bg = xx->qual_bg[0];
	    splodge[j].fg = xx->qual_below;
	    splodge[j].sh |= sh_default | sh_bg | sh_fg;
	} else {
	    col = (10 * (xx->displayedConfidence[j] - q + 1)) / (100 - q + 1);
	    if (col >= 10)
		col = 9;
	    if (col < 0)
		col = 0;
	    if (!(splodge[j].sh & sh_bg))
		splodge[j].bg = xx->qual_bg[col];
	    splodge[j].sh |= sh_default | sh_bg;
	}
    }
}

/*
 * Redraws the main sequences panel, including all the various highlights
 * such as confidence values, tags, show edits, etc.
 */
static void tk_redisplaySeqSequences(EdStruct *xx, int *seqList) {
    int i, search;
    char seq_str[MAX_DISPLAY_WIDTH+21];
    int pos = xx->displayPos;
    int width = xx->displayWidth;
    int diff = 0;
    float set_qual[MAX_DISPLAY_WIDTH+21];
    int set_cons;

    /* Get the consensus sequence */
    if (xx->refresh_flags & ED_DISP_CONS) {
	DBcalcConsensus(xx,pos,width,xx->displayedConsensus,
			xx->displayedConfidence, BOTH_STRANDS);
    }

    /* Draw ruler */
    if (xx->rulerDisplayed &&
	(xx->refresh_flags & ED_DISP_RULER || xx->unpadded_ruler)) {
	int seq_offset;
	
	seq_offset = generate_ruler(xx, seq_str, NULL, pos, width);
	XawSheetPutText(&xx->ed->sw,diff, 0, (Dimension)width,
			&seq_str[seq_offset]);
    }
    xx->refresh_flags &= ~ED_DISP_RULER;
    
    /* Are we looking to refresh just a single sequence? */
    search = ((xx->refresh_flags & ED_DISP_SEQ) &&
	      !(xx->refresh_flags & ED_DISP_READS) &&
	      (xx->showDifferences == 0)) ? 1 : 0;
    
    /* Loop through all visible sequences */
    for (i=xx->lines_per_seq-1 ; i < xx->displayHeight + xx->lines_per_seq-1;
	 i+=xx->lines_per_seq ) {
	int k;
	char * ptr;
	XawSheetInk splodge[MAX_DISPLAY_WIDTH+1];

	/* Initialise to prevent workshop warnings */
	memset(splodge, 0, (MAX_DISPLAY_WIDTH+1) * sizeof(splodge[0]));
	
	/* Set 'k' to an index into the visible sequences */
	if (i < xx->displayHeight-1)
	    k = i/xx->lines_per_seq + xx->displayYPos;
	else
	    k = (xx->totalHeight-1)/xx->lines_per_seq;
	
	/*
	 * Make sure the consensus is always visible and at the bottom of the
	 * screen, even if it means leaving a gap.
	 */
	if (seqList[k] == 0) {
	    if (i < xx->displayHeight-1) {
		/* blank the gap */
		char blank[MAX_DISPLAY_WIDTH+1];
		memset(blank, ' ', MAX_DISPLAY_WIDTH);
		XawSheetPutText(&xx->ed->sw, diff,
				i+xx->rulerDisplayed, (Dimension)width,
				blank);
	    }
	    i = xx->displayHeight-1;
	}

	/* Get relevant portion of the sequence or consensus */
	set_cons = 0;
	if (seqList[k] == 0) {
	    if (!(xx->refresh_flags & ED_DISP_CONS))
		continue;
	    
	    ptr = xx->displayedConsensus;
	} else {
	    int set = xx->set ? xx->set[seqList[k]] : 0;
	    if (search && seqList[k] != xx->refresh_seq)
		continue;
	    
	    DBgetSequence(xx, seqList[k], pos-DB_RelPos(xx,seqList[k]),
			  width, seq_str);
	    if (set) {
		if (xx->set_collapsed && xx->set_collapsed[set]) {
		    DBcalcSetConsensus(xx, pos, width, set, seq_str, set_qual);
		    set_cons = 1;
		}
	    }
	    ptr = seq_str;
	}

	/* Colour it according to the displayed tags */
	getTagSplodge(xx, seqList[k],pos - DB_RelPos(xx,seqList[k]),
		      width, splodge);
	
	/* Show edits */
	if (xx->show_edits && seqList[k] != 0) {
	    tk_redisplaySeqEdits(xx, splodge, seqList[k], pos, width);
	}
	
	/* Colour diffs if required */
	if (seqList[k] != 0 && xx->showDifferences) {
	    tk_redisplaySeqDiffs(xx, splodge, ptr, pos, width, seqList[k]);
	}
	
	/* Remove cutoff tags displayed when not in 'reveal cutoffs' mode. */
	if (!xx->reveal_cutoffs) {
	    int j;
	    
	    for (j=0; j<width; j++) {
		if (ptr[j] == ' ')
		    splodge[j].sh &= ~(sh_fg | sh_bg);
	    }
	}
	
#ifdef SHOW_QUAL
	if (seqList[k] != 0) {
	    if (set_cons) {
		tk_redisplaySetQual(xx, splodge, set_qual, width);
	    } else {
		tk_redisplaySeqReadQual(xx, splodge, seqList[k], pos, width);
	    }
	} else {
	    if (xx->show_cons_qual)
		tk_redisplaySeqConsQual(xx, splodge, seqList[k],pos,width);
	}
#endif /* SHOW_QUAL */
	
	/* And finally draw it */
	XawSheetPutJazzyText(&xx->ed->sw,diff, i+xx->rulerDisplayed,
			     (Dimension)width,ptr,splodge);

	/* Display joined disagreements whenever we redraw the consensus */
	if (seqList[k] == 0 && inJoinMode(xx) &&
	    !(xx->refresh_flags & ED_DISP_NO_DIFFS))
	    redisplayDisagreement(xx);
    }
    
    xx->refresh_flags &= ~(ED_DISP_SEQS | ED_DISP_SEQ);
}


/*
 * Redisplays the editor status lines (strands, translations)
 */
void tk_redisplaySeqStatusCompute(EdStruct *xx, int pos, int width) {
    xx->status_depth = 0;

    if (xx->status[EDITOR_SL_STRAND]) {
	/* Strands */
	int l = xx->status_depth++;
	xx->status_lines =
	    (EdStatus *)xrealloc(xx->status_lines,
				 xx->status_depth * sizeof(EdStatus));
	status_strand(xx, pos, width,
		      xx->status_lines[l].colours,
		      xx->status_lines[l].line,
		      xx->status_lines[l].name);
    }

    if (xx->status[EDITOR_SL_AUTOTRANSLATE]) {
	/* Automatic translations */
	find_exons(xx, pos, width, 0);

    } else if (xx->status[EDITOR_SL_FRAME1p] ||
	       xx->status[EDITOR_SL_FRAME2p] ||
	       xx->status[EDITOR_SL_FRAME3p] ||
	       xx->status[EDITOR_SL_FRAME1c] ||
	       xx->status[EDITOR_SL_FRAME2c] ||
	       xx->status[EDITOR_SL_FRAME3c]) {
	/* Translate selected frames */
	find_exons(xx, pos, width, 1);
    }
}

/*
 * Redisplays the editor status lines (strands, translations)
 */
static void tk_redisplaySeqStatusDisplay(EdStruct *xx) {
    int same_cons;
    int i, diff = 0;
    XawSheetInk nsplodge[NAMELEN];
    int width = xx->displayWidth;

    for (i = 0; i < NAMELEN; i++)
	nsplodge[i].sh = sh_default;

    /*
     * Add status lines
     */
    same_cons = 0;
    for (i = 0; i < xx->status_depth; i++) {
	XawSheetPutText(&xx->names->sw,
			0,
			i + xx->rulerDisplayed + xx->displayHeight - 1,
			NAMELEN,
			xx->status_lines[i].name);

	XawSheetPutJazzyText(&xx->ed->sw,
			     diff,
			     i + xx->rulerDisplayed + xx->displayHeight,
			     (Dimension)width,
			     xx->status_lines[i].line,
			     xx->status_lines[i].colours);
    }

    if (xx->status_depth) {
	sheet_draw_separator(TKSHEET(xx->names),
			     xx->rulerDisplayed + xx->displayHeight - 2);
	sheet_draw_separator(TKSHEET(xx->ed),
			     xx->rulerDisplayed + xx->displayHeight - 1);
    }
}

/*
 * The main editor redraw function.
 */
void tk_redisplaySequences(EdStruct *xx) {
    int *seqList;
    int pos = xx->displayPos;
    int width = xx->displayWidth;
    int cur_depth;

    if (!xx->ed || xx->editorState == StateDown)
	return;

    /* Work out the status line; this may control window height */
    cur_depth = xx->status_depth;
    if (xx->refresh_flags & (ED_DISP_STATUS | ED_DISP_SCROLL)) {
	tk_redisplaySeqStatusCompute(xx, xx->displayPos, xx->displayWidth);
    }

    /* Deal with ED_DISP_SCROLL and ED_DISP_YSCROLL events */
    tk_redisplaySeqScroll(xx, cur_depth != xx->status_depth);

    /* Find out what's visible */
    seqList = sequencesOnScreen(xx,pos, width);

    /* Update sequence names list */
    if (xx->refresh_flags & (ED_DISP_NAMES | ED_DISP_NAME)) {
	tk_redisplaySeqNames(xx, seqList);
    }

    /* Consensus emacs-style edit-status ----, -%%-, -**- */
    tk_redisplaySeqEditStatus(xx);

    /* Redraw the main sequences section, including consensus */
    if (xx->refresh_flags & (ED_DISP_SEQS | ED_DISP_SEQ)) {
	tk_redisplaySeqSequences(xx, seqList);
    }

    /* We've already computed them, but now we actually display them */
    if (xx->refresh_flags & ED_DISP_STATUS) {
	tk_redisplaySeqStatusDisplay(xx);
    }

    /* Editor cursor position */
    if (xx->refresh_flags & ED_DISP_CURSOR) {
	positionCursor(xx, xx->cursorSeq, xx->cursorPos);
    }

    /* Underlining for current selection */
    if (xx->refresh_flags & ED_DISP_SELECTION) {
	redisplaySelection(xx);
    }

    xx->refresh_flags = 0;
    xx->refresh_seq = 0;
}

int redisplayDisagreement(EdStruct *xx) {
    char spare[MAX_DISPLAY_WIDTH+1];
    int i;
    EdStruct *xx0, *xx1;
    
    if (!xx->ed || xx->editorState == StateDown)
	return 0;

    if (inJoinMode(xx) && xx->link) {
        xx0 = xx->link->xx[0];
        xx1 = xx->link->xx[1];

        for (i=0;i<xx->displayWidth;i++) {
            if (i + xx0->displayPos <= 0 ||
                i + xx1->displayPos <= 0 ||
                i + xx0->displayPos > DB_Length(xx0, 0) ||
                i + xx1->displayPos > DB_Length(xx1, 0))
                spare[i] = ' ';
            else
                spare[i]=(xx0->displayedConsensus[i] ==
                          xx1->displayedConsensus[i])?' ':'!';
        }
        XawSheetPutText(&xx->link->diffs->sw, 0, 0,
			(Dimension)xx->displayWidth, spare);
    }

    return 0;
}


/*
 *---------------------------------------------------------------------------
 * Selection handling
 *---------------------------------------------------------------------------
 */


/*
 * Draws the selection
 */
static void toggle_select(EdStruct *xx, int seq, int from_pos, int to_pos) {
    int *seqList;
    int s_from,s_to;
    int screenRow;

    if (from_pos > to_pos) {
	int temp = from_pos;
	from_pos = to_pos;
	to_pos = temp;
    }

    /* clip to screen */
    s_from = positionInContig(xx,seq,from_pos) - xx->displayPos
	- DB_Start(xx, seq);
    if (s_from>=xx->displayWidth) return;
    if (s_from<0) s_from=0;
    s_to = positionInContig(xx,seq,to_pos) - xx->displayPos
	- DB_Start(xx, seq);
    if (s_to<0) return;
    if (s_to>=xx->displayWidth) s_to = xx->displayWidth-1;

    seqList = sequencesOnScreen(xx,xx->displayPos, xx->displayWidth);
    if (seq != 0) {
	for (screenRow = 0;
	     seqList[screenRow] && seqList[screenRow] != seq;
	     screenRow++)
	    ;
	if (seqList[screenRow] != seq)
	    return;

	screenRow = (screenRow - xx->displayYPos) * xx->lines_per_seq;
    } else {
	screenRow = xx->displayHeight - 1;
    }
    XawSheetOpHilightText(&xx->ed->sw,
			  s_from,screenRow + xx->rulerDisplayed,
			  (Dimension)(s_to-s_from+1), sh_select, HOP_TOG);

}


static void redisplaySelection(EdStruct *xx) {
    if (!xx->ed || xx->editorState == StateDown)
	return;

    if (!xx->select_made || xx->select_start_pos == xx->select_end_pos)
	return;

    if (xx->select_start_pos < xx->select_end_pos) {
	toggle_select(xx, xx->select_seq, xx->select_start_pos,
		    xx->select_end_pos-1);
    } else {
	toggle_select(xx, xx->select_seq, xx->select_end_pos,
		    xx->select_start_pos-1);
    }
}

/*
 * Called from tk when our selection has been taken from us
 */
static void EdSelectionLost(ClientData cd) {
    EdStruct *xx = (EdStruct *)cd;

    /* Undisplay the selection */
    redisplaySelection(xx);

    /* And reinitialise */
    xx->select_made = 0;
    xx->select_seq = 0;
    xx->select_start_pos = 0;
    xx->select_end_pos = 0;
}

void edSelectClear(EdStruct *xx) {
    if (EDTKWIN(xx->ed))
	Tk_ClearSelection(EDTKWIN(xx->ed), XA_PRIMARY);
    EdSelectionLost((ClientData)xx);
}

void edSelectFrom(EdStruct *xx, int pos) {
    /* Undisplay an old selection */
    if (xx->select_made)
	redisplaySelection(xx);
    else
	xx->select_made = 1;

    xx->select_seq = xx->cursorSeq;

    pos = xx->displayPos - DB_RelPos(xx,xx->select_seq) + pos + 1
        + DB_Start(xx, xx->select_seq);
    
    if (xx->reveal_cutoffs) {
	if (pos < 1)
	    pos = 1;
	else if (pos > DB_Length2(xx,xx->select_seq)+1)
	    pos = DB_Length2(xx,xx->select_seq)+1;
    } else {
	if (pos < DB_Start(xx, xx->select_seq) + 1)
	    pos = DB_Start(xx, xx->select_seq) + 1;
	else if (pos > DB_Length(xx,xx->select_seq) +
		 DB_Start(xx, xx->select_seq) + 1)
	    pos = DB_Length(xx,xx->select_seq) +
		 DB_Start(xx, xx->select_seq) + 1;
    }

    xx->select_start_pos = pos;
    xx->select_end_pos = xx->select_start_pos;
    xx->select_tag = NULL;

    Tk_OwnSelection(EDTKWIN(xx->ed), XA_PRIMARY, EdSelectionLost,
		    (ClientData)xx);

    /* Display new selection */
    redisplaySelection(xx);
}

void edSelectTo(EdStruct *xx, int pos) {
    /* Perform a select_from if we don't own the selection */
    if (!xx->select_made) {
	/* edSelectFrom(xx, pos); */
	return;
    }

    /* Undisplay old selection */
    redisplaySelection(xx);

    pos = xx->displayPos - DB_RelPos(xx,xx->select_seq) + pos + 1
        + DB_Start(xx, xx->select_seq);
    
    if (xx->reveal_cutoffs) {
	if (pos<1)
	    pos = 1;
	else if (pos > DB_Length2(xx,xx->select_seq)+1)
	    pos = DB_Length2(xx,xx->select_seq)+1;
    } else {
	if (pos < DB_Start(xx, xx->select_seq) + 1)
	    pos = DB_Start(xx, xx->select_seq) + 1;
	else if (pos > DB_Length(xx,xx->select_seq) +
		 DB_Start(xx, xx->select_seq) + 1)
	    pos = DB_Length(xx,xx->select_seq) +
		 DB_Start(xx, xx->select_seq) + 1;
    }

    xx->select_end_pos = pos;

    /* display new selection */
    redisplaySelection(xx);
}

/*
 * Automatically called when X wishes to obtain a selection. We register
 * this procedure in the initialise code of tkEditor.c.
 *
 * Return codes expected:
 *    -1  Failure
 *    >0  Number of bytes
 */
int edGetSelection(ClientData clientData, int offset, char *buffer,
		   int bufsize) {
    Editor *ed = (Editor *)clientData;
    int start, end, len;
    EdStruct *xx = ed->xx;

    /* Do we have a selection? */
    if (!xx->select_made)
	return -1;

    start = xx->select_start_pos + offset;
    end   = xx->select_end_pos;

    if (start > end) {
	len = start;
	start = end;
	end = len;
    }

    len = end - start > bufsize ? bufsize : end - start;

    if (len) {
	if (xx->select_seq) {
	    start -= DB_Start(xx, xx->select_seq) + 1;
	    DBgetSequence(xx, xx->select_seq, start, len, buffer);
	} else {
	    DBcalcConsensus(xx, start - DB_Start(xx, xx->select_seq),
			    len, buffer, NULL, BOTH_STRANDS);
	}
    }

    return len;
}


/*
 * Fills in information about a selection.
 * Returns 0 if no selection made, otherwise returns 1.
 */
int getSelection(EdStruct *xx, int *seq, int *start, int *length,
		 tagStruct **t)
{
    if (! xx->select_made) return 0;

    if (xx->select_start_pos <= xx->select_end_pos) {
        *seq = xx->select_seq;
        *start = xx->select_start_pos;
        *length = xx->select_end_pos - xx->select_start_pos;
	if (t)
	    *t = xx->select_tag;
    } else {
        *seq = xx->select_seq;
        *start = xx->select_end_pos;
        *length = xx->select_start_pos - xx->select_end_pos;
	if (t)
	    *t = xx->select_tag;
    }
    return 1;
}


/*
 * Creates a selection spanning the length of the tag underneath the cursor
 */
void _select_tag(EdStruct *xx, int seq, tagStruct *t) {
    if (t==NULL) return;

    /* Undisplay an old selection */
    if (xx->select_made)
	redisplaySelection(xx);
    else
        xx->select_made = 1;

    xx->select_seq = seq;
    xx->select_start_pos = normalisePos2(xx, seq, t->tagrec.position,
					 t->tagrec.length);
    xx->select_end_pos = xx->select_start_pos + t->tagrec.length;
    xx->select_tag = t;

    Tk_OwnSelection(EDTKWIN(xx->ed), XA_PRIMARY, EdSelectionLost,
		    (ClientData)xx);
    xx->refresh_flags |= ED_DISP_SELECTION;

    /* Update the brief line */
    edSetBriefTag(xx, seq, t, get_default_string(EDINTERP(xx->ed), gap_defs,
						 "TAG_BRIEF_FORMAT"));
    /*
    {
    char status_buf[1024];
    sprintf(status_buf,
	    "Tag #%d: Position %d+%d %c%c%c%c %d \"%.80s\"",
	    t->original_tag_id, t->tagrec.position, t->tagrec.length,
	    t->tagrec.type.c[0], t->tagrec.type.c[1],
	    t->tagrec.type.c[2], t->tagrec.type.c[3],
	    normaliseSense(xx, seq, t->tagrec.sense),
	    t->newcomment ? t->newcomment : "(no comment)");
    tk_update_brief_line(xx, status_buf);
    }
    */

    redisplaySelection(xx);
}


/*
 * Creates a selection spanning the a specific sequence region
 */
void _select_region(EdStruct *xx, int seq, int pos, int length) {
    /* Undisplay an old selection */
    if (xx->select_made)
	redisplaySelection(xx);
    else
        xx->select_made = 1;

    xx->select_seq = seq;
    xx->select_start_pos = pos;
    xx->select_end_pos = pos + length;
    xx->select_tag = NULL;

    Tk_OwnSelection(EDTKWIN(xx->ed), XA_PRIMARY, EdSelectionLost,
		    (ClientData)xx);
    /*
     * FIXME:
     * ED_DISP_SELECTION should be all we need to update, but at present
     * this does not seem to be redrawing properly.
     * ED_DISP_ALL is an easy workaround for now.
     */
    xx->refresh_flags |= ED_DISP_ALL;

    redisplaySelection(xx);
}


/*
 * Adjusts the current selection when an insert is performed. (We either move
 * the selection or extend it's length.)
 */
void selectInsertBase(EdStruct *xx, int seq, int pos) {
    pos += DB_Start(xx, seq);

        if (xx->select_made && xx->select_seq==seq) {
        int inverted=(xx->select_end_pos < xx->select_start_pos);
        int start,end;

        if (inverted) {
            start=xx->select_end_pos;
            end  =xx->select_start_pos;
        } else {
            end  =xx->select_end_pos;
            start=xx->select_start_pos;
        }

        if (pos <= start) {
            xx->select_start_pos++;
            xx->select_end_pos++;
        } else if (pos < end) {
            if (inverted)
                xx->select_start_pos++;
            else
                xx->select_end_pos++;
        }
    }
}


/*
 * Adjusts the current selection when a delete is performed. (We either move
 * the selection or extend it's length.)
 */
void selectDeleteBase(EdStruct *xx, int seq, int pos) {
    pos += DB_Start(xx, seq);

    if (xx->select_made && xx->select_seq==seq) {
        int inverted=(xx->select_end_pos < xx->select_start_pos);
        int start,end;

        if (inverted) {
            start=xx->select_end_pos;
            end  =xx->select_start_pos;
        } else {
            end  =xx->select_end_pos;
            start=xx->select_start_pos;
        }

        if (pos < start) {
            xx->select_start_pos--;
            xx->select_end_pos--;
        } else if (pos < end) {
            if (inverted)
                xx->select_start_pos--;
            else
                xx->select_end_pos--;
        }
    }
}


/*
 *---------------------------------------------------------------------------
 * Status line calculations
 *---------------------------------------------------------------------------
 */
void status_strand(EdStruct *xx, int pos, int width,
		   XawSheetInk *splodge, char *sline,
		   char *name) {
    char qual[MAX_DISPLAY_WIDTH];
    int j;

    calc_quality(0, pos, pos + width - 1, qual,
		 xx->con_cut, xx->qual_cut,
		 contEd_info, (void *)xx);

    for (j = 0; j < width; j++) {
	splodge[j].sh = sh_default;
	
	switch (qual[j]) {
	case R_GOOD_GOOD_NE:
	    sline[j] = '!';
	    break;

	case R_GOOD_GOOD_EQ:
	case R_GOOD_BAD:
	case R_BAD_GOOD:
	case R_BAD_BAD:
	    sline[j] = '=';
	    break;
	    
	case R_GOOD_NONE:
	case R_BAD_NONE:
	    sline[j] = '+';
	    break;
	    
	case R_NONE_GOOD:
	case R_NONE_BAD:
	    sline[j] = '-';
	    break;
	    
	default: /* R_NONE_NONE */
	    sline[j] = ' ';
	}
    }

    sprintf(name, " %*s %-*s", DB_GELNOLEN, " ", DB_NAMELEN, "Strands");
}

/*
 * Updates the 'brief' (NB: not a 'status' line in the same sense as the
 * above uses) line.
 *
 * To optimise usage of this, it returns a 'code number' which is updated
 * each time a new message is displayed. Passing 'msg' as NULL will return
 * the current code. Therefore calling routines can optimise by knowing
 * whether or not an update is required.
 */
int tk_update_brief_line(EdStruct *xx, char *msg) {
    static int last_code = 0;
    Tcl_DString ds;
    int i, diff;
    static char tmp_msg[1025];

    if (NULL == msg)
	return last_code;

    if (NULL == xx->ed->highlight_cmd)
	return 0;

    /* Copy, compare and convert msg */
    diff = 0;
    for (i = 0; msg[i] && i < 1024; i++) {
	char c;
	c = (msg[i] == '\n' || msg[i] == '\r') ? ' ' : msg[i];
	if (diff == 0 && tmp_msg[i] != c)
	    diff = 1;
	tmp_msg[i] = c;
    }
    if (tmp_msg[i] != 0)
	diff = 1;
    tmp_msg[i] = 0;

    /* No difference from last line, so do nothing */
    if (!diff)
	return last_code;

    Tcl_DStringInit(&ds);
    Tcl_DStringAppend(&ds, xx->ed->highlight_cmd, -1);
    Tcl_DStringAppend(&ds, " ", 1);
    Tcl_DStringAppendElement(&ds, tmp_msg);
    if (TCL_OK != Tcl_Eval(EDINTERP(xx->ed), Tcl_DStringValue(&ds))) {
	fprintf(stderr, "Tcl_Eval: %s\n", Tcl_GetStringResult(EDINTERP(xx->ed)));
    }

    Tcl_DStringFree(&ds);

    return ++last_code;
}

/*
 *---------------------------------------------------------------------------
 * Miscellaneous
 *---------------------------------------------------------------------------
 */

void delete_edStruct(EdStruct *xx) {
    extern int edused[MAXEDSTATES];
    extern EdStruct edstate[MAXEDSTATES];
    int j, inuse = 0;

            
    for (j = 0; j < MAXEDSTATES; j++) {
	if (edused[j] && edstate[j].DBi == DBI(xx))
	    inuse++;
    }
    
    if (inuse == 1) {
	freeAllUndoLists(xx);
    }

    if (EDTKWIN(xx->ed))
	Tk_ClearSelection(EDTKWIN(xx->ed), XA_PRIMARY);
    XSync(Tk_Display(Tk_MainWindow(EDINTERP(xx->ed))), True);

    if (inJoinMode(xx) && xx->link) {
	DestroyEdLink(xx->link);
    }

    freeDB(xx, 1);
}


int tk_edid_to_editor(ClientData clientData, Tcl_Interp *interp, 
		       int argc, char *argv[]) {
    extern int edused[MAXEDSTATES];
    extern EdStruct edstate[MAXEDSTATES];
    int i, edid;

    if (argc != 2) {
	return TCL_ERROR;
    }
    edid = atoi(argv[1]);

    for (i = 0; i < MAXEDSTATES; i++) {
	if (edused[i] && edstate[i].editor_id == edid)
	    break;
    }

    if (i != MAXEDSTATES) {
	vTcl_SetResult(interp, "%s", Tk_PathName(EDTKWIN(edstate[i].ed)));
    }

    return TCL_OK;
}

/*
 * Set the name X display position
 */
void setNamePos(EdStruct *xx, int pos)
{
    xx->refresh_flags |= ED_DISP_NAMES;
    xx->names_xpos = pos;

    redisplaySequences(xx, 0);
}
