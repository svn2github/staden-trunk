#include <staden_config.h>

#include <tk.h>
#include <string.h>
#include <X11/Xlib.h>

#include "tkSeqed.h"
#include "tkSeqedUtils.h"
#include "tkSheet_common.h"
#include "dna_utils.h"
#include "genetic_code.h"
#include "misc.h"
#include "seqed_translate.h"
#include "seq_reg.h"
#include "sheet.h"
#include "seqed_restriction_enzymes.h"
#include "seqed_write.h"
#include "seqed_search.h"
#include "tkSeqedUtils.h"
#include "seq_results.h"

static region *exons;
static int num_exons;

void
initSeqed(tkSeqed *se)
{
    int i;
    se->rid = 0; /* registration id */
    se->seq_id = 0; 
    se->sequence = NULL;
    se->xScrollCmd = NULL;
    se->yScrollCmd = NULL;
    se->displayWidth = 0;
    se->displayHeight = 0;
    se->displayPos = 1;
    se->displayYPos = 0;
    se->cursorPos = 1;
    se->cursorSeq = 0;
    se->extent_left = 0;
    se->extent_right = 0;
    se->num_lines = 0;
    se->rulerDisplayed = 1;
    se->complementDisplayed = 0;
    se->translationDisplayed = 0;
    se->autoDisplayed = 0;
    se->renzDisplayed = 0;
    se->renz_lines = 0;
    se->trans_mode = 1;
    se->trans_lines = 0;
    se->auto_trans_lines = 0;
    for (i = 0; i < MAX_TRANS_LINES; i++) {
	se->trans[i] = -1;
    }
    se->trans_lines_p = 0;
    se->trans_lines_c = 0;
    for (i = 0; i < MAX_LINES; i++) {
	se->lines[i] = 0;
    }
}

void
setDimensions(tkSeqed *se)
{
    se->displayWidth = se->sw.columns;
    se->displayHeight = se->sw.rows;
    /* 
     * set SEQ to be at a constant height given by se->anchor_pos. Initialise
     * it to be 3 up from the bottom - to allow room for ruler and complement
     */ 
    se->anchor_pos = se->displayHeight - 3; 
#ifdef DEBUG
    printf("setDimensions anchor %d height %d\n",  
	   se->anchor_pos, se->displayHeight);
#endif
}

void
set_extents(tkSeqed *se, int length)
{
/*
    se->extent_left = 0;
    se->extent_right = length - 1;
*/
    se->extent_left = 1;       /* start at 1 */
    se->extent_right = length;
}


/*
 * positions the lines within the widget
 */
void set_lines(tkSeqed *se,
	       int set_display_pos,
	       int keep_x_pos)
{
    int i, j;
    int cnt;
    int start;
    static int prev_top = INT_MAX;

#ifdef DEBUG
    printf("\nSTART setlines ypos %d\n", se->displayYPos);
#endif

    /* reinitialise line numbers to -1 */
    for (i = 0; i < num_exons; i++) {
	exons[i].line_num = -1;
    }

    /* initialise previous y position ie top sequence displayed */
    if (prev_top == INT_MAX) {
	prev_top = se->num_lines;
    }
#ifdef DEBUG
    printf("prev %d top %d keep_x %d\n", prev_top, se->num_lines, keep_x_pos);
#endif
    /* sheet_clear(&se->sw);  */

    /* 
     * number of displayed lines has decreased 
     * NB this is out by one call because num_lines is updated at the end of 
     * this routine ie num_lines is prev and prev is prev_prev really!
     * but don't know what I should do anyway with the scrollbar 
     */

    if (prev_top > se->num_lines) {
	se->displayYPos = se->num_lines - (prev_top - se->displayYPos);

	if (se->displayYPos < 0) 
	    se->displayYPos = 0;
#ifdef DEBUG
	printf("CHANGED pos %d\n", se->displayYPos);
#endif
    }


    if (keep_x_pos) {
	if (set_display_pos) {
	    se->displayYPos = set_display_pos;
	} else {
	    se->displayYPos = 0;
	}
    }

    cnt = 0 - se->displayYPos;
    start = cnt;

#ifdef DEBUG
    printf("IN SET LINES %d ypos %d num lines %d\n", 
	   cnt, se->displayYPos, se->num_lines);
#endif
    /*
     * the following assignments MUST be in the order the lines are to appear
     * in the final display
     */

    if (se->renzDisplayed) {
	/* cnt = se->renz_lines; */
	se->lines[RENZ] = cnt;
	cnt += se->renz_lines;
    }

    if (se->translationDisplayed) {
	/* forward translation */
	for (j = 0; j < se->trans_lines; j++) {
	    if (se->trans[j] < 4)
		se->lines[se->trans[j]] = cnt++;
	}
    }
    
    if (se->autoDisplayed) {
	se->auto_trans_lines = find_auto_lines(&exons, num_exons, 
					       se->displayWidth,
					       se->displayPos-1, 0);
	se->lines[AUTO] = cnt;
	cnt += se->auto_trans_lines;
    }

    se->cursorSeq = cnt;
    se->lines[SEQ] = cnt++;

    if (se->rulerDisplayed)
	se->lines[RULER] = cnt++;

    if (se->complementDisplayed)
	se->lines[COMP] = cnt++;
    
    if (se->autoDisplayed) {
	se->auto_c_trans_lines = find_auto_lines(&exons, num_exons, 
					       se->displayWidth,
					       se->displayPos-1, 1);
	se->lines[AUTO_C] = cnt;
	cnt += se->auto_c_trans_lines;
    }

    if (se->translationDisplayed) {
	/* complement translation */
	for (j = 0; j < se->trans_lines; j++) {
	    if (se->trans[j] > 3)
		se->lines[se->trans[j]] = cnt++;
	}
    }

    /* 
     * difference between initial cnt and final cnt is the number of lines
     * to be displayed
     */
#ifdef DEBUG
    printf("cnt %d start %d\n", cnt, start);
#endif
    prev_top = se->num_lines;
    se->num_lines = cnt - start;

#ifdef DEBUG
    printf("MIN %d MAX %d\n", se->heightmin, se->heightmax);
#endif

    /*
     * reduce display height if number of lines is less than the maximum
     * display height else set at maximum display height
     */
/*    if (se->num_lines <= se->displayHeight) { */
#ifdef DEBUG
    printf("num_lines %d heightmin %d heightmax %d\n",
	   se->num_lines, se->heightmin, se->heightmax);
#endif
    if (se->num_lines >= se->heightmin && se->num_lines <= se->heightmax) {

#ifdef DEBUG
       printf("RESET DISPLAY HEIGHT cnt %d ht %d top %d ypos %d\n", cnt, se->displayHeight, se->num_lines, se->displayYPos);
#endif

      /* sheet_set_display_height(TKSHEET(se), se->num_lines); */
    } else {
	/* sheet_set_display_height(TKSHEET(se), se->displayHeight); */
    }
    /*
     * need to update vertical scrollbar position 
     */
    seqed_set_v_sb_pos(se, se->displayYPos);
}

void
seqed_redisplay_translation(tkSeqed *se, 
			    int pos)
{
    int j;
    char sline[MAX_DISPLAY_WIDTH+1];
    /*
    XawSheetInk splodge[MAX_DISPLAY_WIDTH+1];
    char name[10];
    */
    for (j = 0; j < se->trans_lines; j++) {
	/*
	seqed_translate_frame(se, &(se->sequence[pos-1]), pos, 
			      se->displayWidth, se->trans[j], sline, 
			      name, se->trans_mode, splodge); 
        XawSheetPutJazzyText(&se->sw, 0, se->lines[se->trans[j]], 
			     se->displayWidth, sline, splodge);
			      */
	seqed_write_translation(&se->sequence[pos-1], se->trans[j], 
				se->trans_mode, pos, 
				se->displayWidth, 0, sline);
	XawSheetPutText(&se->sw, 0, se->lines[se->trans[j]], 
			se->displayWidth, sline);

    }
}

int
seqed_redisplay_auto_translation(tkSeqed *se, 
				 int pos)
{
    int i;
    char **sline;
    char name[10];
    XawSheetInk **splodge;
    int num_lines = -1;
    int complement;
    int max_lines;

#ifdef DEBUG
    printf("start seqed_redisplay_auto_translation\n");
#endif

    /*
     * allocate memory to store display lines 
     */
    max_lines = se->auto_trans_lines + se->auto_c_trans_lines;
    if (max_lines == 0)
	return 1;

    if (NULL == (splodge = (XawSheetInk **)xmalloc(max_lines * 
						   sizeof(XawSheetInk*))))
	return TCL_OK;

    if (NULL == (sline = (char **)xmalloc(max_lines * sizeof(char*))))
	return TCL_OK;

    for (i = 0; i < max_lines; i++) {
	if (NULL == (sline[i] = (char *)xmalloc((MAX_DISPLAY_WIDTH+1) * 
						sizeof(char))))
	    return TCL_OK;
	if (NULL == (splodge[i] = (XawSheetInk *)xmalloc((MAX_DISPLAY_WIDTH+1)
						       * sizeof(XawSheetInk))))
	    return TCL_OK;
    }

    /* complement */
    for (complement = 0; complement < 2; complement++) {

	if (complement) 
	    num_lines = se->auto_c_trans_lines;
	else
	    num_lines = se->auto_trans_lines;
	
	if (num_lines == 0)
	    continue;
#ifdef DEBUG
	printf("num_lines %d j %d comp %d\n", num_lines, j, complement);
#endif	
	/* 
	 * loop through exons visible on display
	 */
	for (i = 0; i < num_exons; i++) {
	    
	    /* printf("translate i %d line %d\n", i, exons[i].line_num); */
	    if (exons[i].line_num > -1 && exons[i].complement == complement) {

/*
		printf("pos %d i %d line %d start %d end %d\n", 
		       pos, i, exons[i].line_num, exons[i].start, exons[i].end);
*/
		seqed_auto_translate(se, &(se->sequence[pos-1]), pos, 
				     se->displayWidth, sline[exons[i].line_num], name, 
				     splodge[exons[i].line_num], 
				     se->trans_mode, exons, i, exons[i].start,
				     exons[i].end, exons[i].line_num, 
				     complement);
	    }
	}
	/*
	 * print out lines to display
	 */
	for (i = 0; i < num_lines; i++) {
	    if (complement) {
/*
		printf("line number %d\n", se->lines[AUTO_C]);
		printf("sline %d %.*s\n", i, se->displayWidth, sline[i]); 
*/
		XawSheetPutJazzyText(&se->sw, 0, se->lines[AUTO_C]+i,
				     se->displayWidth, sline[i], splodge[i]);
	    } else {
		/* printf("sline %d %.*s\n", i, se->displayWidth, sline[i]); */
		XawSheetPutJazzyText(&se->sw, 0, se->lines[AUTO]+i,
				     se->displayWidth, sline[i], splodge[i]);
	    }
	}
    }

    for (i = 0; i < max_lines; i++) {
	xfree(sline[i]);
	xfree(splodge[i]);
    }
    xfree(sline);
    xfree(splodge);
    return TCL_OK;
}

void seqed_redisplay_sequence(tkSeqed *se, 
			      int pos) 
{
    char sline[MAX_DISPLAY_WIDTH+1];
    int displayWidth;
    int x = 0;
     
    memset(sline, ' ', MAX_DISPLAY_WIDTH);

    if (se->seq_len < se->displayWidth) {
	displayWidth = se->seq_len;
    } else {
	displayWidth = se->displayWidth;
    }
    seqed_write_sequence(&se->sequence[pos+1], pos+1, displayWidth, sline);
    XawSheetPutText(&se->sw, x, se->lines[SEQ], displayWidth, sline);
}

void
seqed_redisplay_complement(tkSeqed *se, 
			   int pos) 
{
    int x = 0;
    char sline[MAX_DISPLAY_WIDTH+1];
    int displayWidth;

    memset(sline, ' ', MAX_DISPLAY_WIDTH);
    if (se->seq_len < se->displayWidth) {
	displayWidth = se->seq_len;
    } else {
	displayWidth = se->displayWidth;
    }
    seqed_write_complement(&se->sequence[pos+1], pos+1, displayWidth, sline);
    XawSheetPutText(&se->sw, x, se->lines[COMP], displayWidth, sline);
}

void
seqed_redisplay_ruler(tkSeqed *se, 
		      int pos) 
{
    int x = 0;
    char sline[MAX_DISPLAY_WIDTH+1];
    int displayWidth;

    if (se->seq_len < se->displayWidth) {
	displayWidth = se->seq_len;
    } else {
	displayWidth = se->displayWidth;
    }
    seqed_write_ruler(pos, displayWidth, sline);
    XawSheetPutText(&se->sw, x, se->lines[RULER], displayWidth, sline);
}

void
seqed_redisplay_seq(tkSeqed *se, 
		    int pos, /* first base of sequence to print (1=first) */
		    int keep_x_pos)
{
    int display_pos;

    pos -= se->extent_left - 1;

    /*
     * bit of a HACK because I don't know how many lines I require to display
     * the restriction enzymes until I have done the calculation and therefore
     * I call set_lines within the seqed_redisplay_renzyme function
     * 
     */
/*    set_lines(se); */
    sheet_clear(&se->sw);

    if (se->renzDisplayed) {
	if (-1 == seqed_redisplay_renzyme(se, pos, keep_x_pos)) {
	    verror(ERR_WARN, "sequence editor", "failed to redisplay restriction enzyme plot \n");
	}
    } else {
	set_lines(se, 0, keep_x_pos);
	display_pos = se->lines[SEQ] - se->anchor_pos;
	set_lines(se, display_pos, keep_x_pos);
    }

/*
    for (i = 0; i < num_exons; i++) {
	printf("set lines i %d line %d\n", i, exons[i].line_num);
    }
*/
    if (se->rulerDisplayed)
	seqed_redisplay_ruler(se, pos);

    if (se->complementDisplayed)
	seqed_redisplay_complement(se, pos);
	
    if (se->autoDisplayed)
	seqed_redisplay_auto_translation(se, pos);

    if (se->translationDisplayed)
	seqed_redisplay_translation(se, pos);

    seqed_redisplay_sequence(se, pos);

    seqed_positionCursor(se, se->cursorSeq, se->cursorPos);

    /* callback to raster to ensure the cursor is plotted there aswell */
    /* can't see the purpose of this; commented out on 2/10/97 */
    /* seqed_setCursorPos(se, se->cursorPos); */

    /* set the scrollbar */
    /*  seqed_set_h_sb_pos(se, se->cursorPos); */
    seqed_set_h_sb_pos(se, pos);

    se->flags |= SHEET_REDRAW_TEXT;
    if (!(se->flags & SHEET_REDRAW_PENDING)) {
	se->flags |= SHEET_REDRAW_PENDING;
	Tcl_DoWhenIdle(SheetDisplay, (ClientData)se);
    }
}

/*
 * add ruler to Seqed widget
 */
void
seqed_add_ruler(tkSeqed *se, 
		int length, 
		char *sequence,
		int x,
		int y) 
{
    se->rulerDisplayed = 1;

    seqed_redisplay_ruler(se, 0);
    
}

/*
 * add sequence to Seqed widget
 * save sequence and it's length in structure
 */
int
seqed_add_sequence(tkSeqed *se, 
		   int length, 
		   char *sequence,
		   char *seq_name,
		   int sequence_type,
		   int seq_id,
		   int x,
		   int y) 
{
    char *new_seq;

    if (NULL == (new_seq = (char *)xmalloc((length+5) * sizeof(char))))
	return TCL_OK;
    
    se->sequence = new_seq;
    se->sequence[0] = ' ';
    se->sequence[1] = ' ';

    strcpy(&(se->sequence[2]), sequence);
    se->sequence[length+2] = ' ';
    se->sequence[length+3] = ' ';
    se->sequence[length+4] = '\0';
    se->seq_len = length;
    se->seq_name = seq_name;
    se->sequence_type = sequence_type;
    se->seq_id = seq_id;

    set_extents(se, length);

    sheet_set_display_height(TKSHEET(se), se->displayHeight);
    seqed_redisplay_seq(se, 1, 1);
    seqed_showCursor(se, se->cursorSeq, se->cursorPos);
    return TCL_OK;
}

void reset_anchor(tkSeqed *se)
{
    
    se->anchor_pos = se->displayHeight - 3; 
#ifdef DEBUG
    printf("reset_anchor pos %d height %d numlines %d\n", 
	   se->anchor_pos, se->displayHeight, se->num_lines);
#endif
}

void seqed_set_h_sb_pos(tkSeqed *se, 
			int pos) 
{
    char buf[100];
    int len = se->displayWidth;
    int total = se->extent_right - se->extent_left + 1;
    
    if (se->seq_len < se->displayWidth) {
	len = se->seq_len;
    } else {
	len = se->displayWidth;
    }

#ifdef DEBUG
    printf("seqed_set_h_sb_pos %d\n", pos);
#endif
/*
    if (!ed || sw->editorState == StateDown)
	return;
*/
    if (se->xScrollCmd) {
	 pos -= se->extent_left - 1;

	/* total width start end */
	 sprintf(buf, " %g %g", pos / (double)total,
		 (pos + se->displayWidth) / (double)total);
	/*
	 * Ignore errors as the widget may not exist (eg when shutting
	 * down the editor). Hacky, but simple.
	 */

        if (Tcl_VarEval(se->interp, se->xScrollCmd, buf, NULL) != TCL_OK) {
            Tcl_AddErrorInfo(se->interp,
                             "\n    (xscrollcommand executed by Sheet)");
            Tcl_BackgroundError(se->interp);
        }
    }
}

void seqed_set_v_sb_pos(tkSeqed *se, 
			int pos) 
{
    double start, end, theight;
    int dheight = se->displayHeight;
    int lines = se->num_lines;
    char buf[100];

    if (!se->yScrollCmd)
	return;

    /*
     * Using the new scrollbar syntax, we set the start and end of the scroll
     * bar slider in fractions from 0 to 1. The total height is:
     *
     *  Max(P+H, N) - MIN(0,P)
     *
     * Therefore, start =       P - min(0, P)
     *                      ----------------------
     *                      Max(P+H, N) - MIN(0,P)
     *
     * And end = start +              h
     *                      ----------------------
     *                      Max(P+H, N) - MIN(0,P)
     *
     * Where P = display position
     *       H = window height
     *       N = number of lines in output
     */

    theight = MAX(pos + dheight, lines) - MIN(0, pos);
    start = (pos - MIN(0, pos)) / theight;
    end = start + dheight / theight;

#ifdef DEBUG
    printf("Start=%f End=%f\n", start, end);
#endif

    sprintf(buf, " %f %f", start, end); 
    /*
     * Ignore errors as the widget may not exist (eg when shutting
     * down). Hacky, but simple.
     */
    if (Tcl_VarEval(se->interp, se->yScrollCmd, buf, NULL) != TCL_OK) {
	Tcl_AddErrorInfo(se->interp,
			 "\n    (yscrollcommand executed by Sheet)");
	Tcl_BackgroundError(se->interp);
    }
}


/*
 * set translation mode to either 1 letter or 3 letter
 */
void
seqedTransMode(tkSeqed *se, 
	       int mode)
{
    se->trans_mode = mode == 3 ? 3 : 1;
    seqed_redisplay_seq(se, se->displayPos, 1);
}

/* 7/1/99 johnt - re implemented ColourNameToPixel to make portable to Windows */
Pixel ColourNameToPixel(Tcl_Interp *interp, Tk_Window tkwin, char *name)
{
    XColor *cp = Tk_GetColor(interp, tkwin, name);

    if (cp) {
	return cp->pixel;
    } else {
	verror(ERR_WARN, "ColourNameToPixel", "Colourmap is full");
	return 0;
    }
}

Pixel get_new_colour(Tcl_Interp *interp,
		     int num_colours)
{
    static int col = 0;
    char colour[20];

    if (col == 0) 
	strcpy(colour, "blue");
    if (col == 1) 
	strcpy(colour, "red");
    if (col == 2) 
	strcpy(colour, "green");
    if (col == 3) 
	strcpy(colour, "purple");
    if (col == 4) 
	strcpy(colour, "brown");
    if (col == 5) 
	strcpy(colour, "yellow");
    if (col == 6) 
	strcpy(colour, "cyan");
    if (col == 7) 
	strcpy(colour, "hotpink");
    if (col == 8) 
	strcpy(colour, "orange");
    if (col == 9) 
	strcpy(colour, "yellowgreen");
    if (col == 10) 
	strcpy(colour, "coral");

    col++;
    return ColourNameToPixel(interp,Tk_MainWindow(interp),colour); /* 7/1/99 johnt - use new format ColourNameToPixel */
}

/*
 * assume exons are in order at this point
 */  
#define G4
int
parse_feature_table(Tcl_Interp *interp,
		    tkSeqed *se)
{
    int i;

#ifdef HSPROPG    
    num_exons = 9;
#endif

#ifdef G4
    num_exons = 11;
#endif
    
#ifdef PS2CG
    num_exons = 10;
#endif 

    if (NULL == (exons = (region *)xmalloc((num_exons) * sizeof(region))))
	return TCL_OK;

#ifdef HSPROPG    
    exons[0].start = 2206;
    exons[0].end = 2281;
    exons[0].join = -1;
    exons[0].complement = 0;
    exons[1].start = 2376;
    exons[1].end = 2526;
    exons[2].start = 3774;
    exons[2].end = 3949;
    exons[3].start = 4410;
    exons[3].end = 4580;
    exons[4].start = 4719;
    exons[4].end = 4910;
    exons[5].start = 5105;
    exons[5].end = 5278;
    exons[6].start = 5532;
    exons[6].end = 5723;
    exons[7].start = 5882;
    exons[7].end = 5993;
    exons[8].start = 7610;
    exons[8].end = 7775;

    for (i = 1; i < num_exons; i++) {
	exons[i].join = i-1;
	exons[i].complement = 0;
    }
#endif

#ifdef G4
    exons[0].start = 59;
    exons[0].end = 1723;
    exons[1].start = 698;
    exons[1].end = 1720;
    exons[2].start = 1276;
    exons[2].end = 1638;
    exons[3].start = 1638;
    exons[3].end = 1808;
    exons[4].start = 1720;
    exons[4].end = 1974;
    exons[5].start = 1976;
    exons[5].end = 2434;
    exons[6].start = 2154;
    exons[6].end = 2444;
    exons[7].start = 2477;
    exons[7].end = 2554;
    exons[8].start = 2600;
    exons[8].end = 3883;
    exons[9].start = 4020;
    exons[9].end = 4553;
    exons[10].start = 4564;
    exons[10].end = 5577;
    for (i = 0; i < num_exons; i++) {
	exons[i].join = -1;
	exons[i].complement = 0;
    }

#endif

#ifdef PS2CG
    exons[0].start = 175;
    exons[0].end = 748;
    exons[0].join = -1;
    exons[0].complement = 0;
    exons[1].start = 175;
    exons[1].end = 748;
    exons[1].join = -1;
    exons[1].complement = 0;
    exons[2].start = 175;
    exons[2].end = 411;
    exons[2].join = -1;
    exons[2].complement = 0;
    exons[3].start = 797;
    exons[3].end = 810;
    exons[3].join = 0;
    exons[3].complement = 0;
    exons[4].start = 797;
    exons[4].end = 2917;
    exons[4].join = 2;
    exons[4].complement = 0;
    exons[5].start = 811;
    exons[5].end = 1502;
    exons[5].join = 1;
    exons[5].complement = 0;
    exons[6].start = 2925;
    exons[6].end = 4076;
    exons[6].join = -1;
    exons[6].complement = 1;
    exons[7].start = 4045;
    exons[7].end = 5004;
    exons[7].join = -1;
    exons[7].complement = 1;
    exons[8].start =4045 ;
    exons[8].end = 4659;
    exons[8].join = -1;
    exons[8].complement = 1;
    exons[9].start = 5285;
    exons[9].end = 5297;
    exons[9].join = -1;
    exons[9].complement = 1;
#endif

    exons[0].num_char = exons[0].end % 3;
    exons[0].colour = se->sw.foreground;

    for (i = 1; i < num_exons; i++) {
	if (exons[i].join != -1) {
	    /* found a join */
	    exons[i].num_char = ( exons[i].end - (exons[i].start - 
					exons[exons[i].join].num_char)+1)%3;

	    exons[i].colour = get_new_colour(interp, num_exons);
	    exons[exons[i].join].colour = exons[i].colour;
	} else {
	    exons[i].num_char = ((exons[i].end - exons[i].start) + 1) % 3;

	    exons[i].colour = se->sw.foreground;
	} 
    } 
#ifdef DEBUG
    for (i = 0; i < num_exons; i++)
	printf("EXON %d %d\n", i, exons[i].num_char);
#endif
    return TCL_OK;
}

void
seqedTranslateAdd(Tcl_Interp *interp,
		  tkSeqed *se,
		  int mode)
{
    int i;

    if (mode == AUTO) {
	se->autoDisplayed = 1;
	parse_feature_table(interp, se);
	seqed_redisplay_seq(se, se->displayPos, 1);
	return;
    } 

    se->translationDisplayed = 1;
    
    /* Delete if already existing (so we can move it to the end) */
    for (i = 0; i < se->trans_lines; i++) {
	if (se->trans[i] == mode) {
	    memmove(se->trans + i, se->trans + i + 1,
		    (MAX_TRANS_LINES - i - 1) * sizeof(*se->trans));
	    se->trans_lines--;
	    break;
	}
    }

    /* Add mode */
    if (se->trans_lines < MAX_TRANS_LINES) {
	se->trans[se->trans_lines++] = mode;
    }

    seqed_redisplay_seq(se, se->displayPos, 1);

}
/*
 * Turns off an seqed translation . If it's already off then we do nothing.
 */
void seqedTranslateDelete(tkSeqed *se,
			  int mode) 
{
    int i;

    if (mode == AUTO) {
	se->autoDisplayed = 0;
	seqed_redisplay_seq(se, se->displayPos, 1);
	return;
    }

    for (i = 0; i < se->trans_lines; i++) {
	if (se->trans[i] == mode) {
	    memmove(se->trans + i, se->trans + i + 1,
		    (MAX_TRANS_LINES - i - 1) * sizeof(*se->trans));
	    se->trans_lines--;
	    break;
	}
    }
    

    if (se->trans_lines == 0 ) 
	se->translationDisplayed = 0;
    else
	se->translationDisplayed = 1;

    reset_anchor(se);
    seqed_redisplay_seq(se, se->displayPos, 1);
}

void
seqed_add_renzyme(tkSeqed *se,
		  char *filename,
		  char *list,
		  int num_items)
{
    seqedREnzyme(se, filename, list, num_items, se->displayWidth);

    se->renzDisplayed = 1;
    seqed_redisplay_seq(se, se->displayPos, 1);
}

void
seqedSeqInfo(tkSeqed *se,
	     int x,
	     int y)
{


}

/*
 *----------------------------------------------------------------------------
 * Scrolling
 *---------------------------------------------------------------------------
 */

/*
 * Increase the leftmost base position on the screen by a symbolic amount
 */
static void seqed_incDisplayPosP (tkSeqed *se, int distance)
{
    switch (distance) {
    case D_screen     : se->displayPos += se->displayWidth; break;
    case D_halfScreen : se->displayPos += se->displayWidth/2; break;
    case D_character  : se->displayPos += 1; break;
    }

    if (se->displayPos > se->extent_right + 2 - se->displayWidth)
	se->displayPos = se->extent_right + 2 - se->displayWidth;

    seqed_redisplay_seq(se, se->displayPos, 1); 
}


/*
 * Decrease the leftmost base position on the screen by a symbolic ammount
 */
static void seqed_decDisplayPosP (tkSeqed *se, int distance)
{
    switch (distance) {
    case D_screen     : se->displayPos -= se->displayWidth; break;
    case D_halfScreen : se->displayPos -= se->displayWidth/2; break;
    case D_character  : se->displayPos -= 1; break;
    }

    if (se->displayPos < se->extent_left) se->displayPos = se->extent_left;
    
    if (se->displayPos > se->extent_right + 2 - se->displayWidth)
	se->displayPos = se->extent_right + 2 - se->displayWidth;

    seqed_redisplay_seq(se, se->displayPos, 1); 
}


/*
 * Increase the leftmost base position on the screen by a symbolic ammount
 */
void seqed_incDisplayPos (tkSeqed *se, int distance)
{
/*
    if (editorLocked(xx)) {
	incDisplayPosP(xx->link->xx[0], distance);
	incDisplayPosP(xx->link->xx[1], distance);
    } else
*/
	seqed_incDisplayPosP(se, distance);
    
    /* redisplayDisagreement(xx); */
}


/*
 * Decrease the leftmost base position on the screen by a symbolic ammount
 */
void seqed_decDisplayPos (tkSeqed *se, int distance)
{
/*
    if (editorLocked(xx)) {
	decDisplayPosP(xx->link->xx[0], distance);
	decDisplayPosP(xx->link->xx[1], distance);
    } else
*/
	seqed_decDisplayPosP(se, distance);
   
    /* redisplayDisagreement(xx); */
}



/*
 * Adjust displayed region
 */
void seqed_setDisplayPosP(tkSeqed *se, int pos)
{
#ifdef REMOVE
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

	ed_set_h_sb_pos(otherxx, otherxx->displayPos);
	redisplaySequences(otherxx, 0);
    }
#endif   

    se->displayPos = pos;
}

/*
 * centralise pos on screen
 */
void seqed_setDisplayPos(tkSeqed *se, int pos)
{
    seqed_setDisplayPosP(se, pos);
    seqed_decDisplayPos(se, D_halfScreen);
}


/*
 * Set left of screen to pos
 */
void seqed_setDisplayPos2(tkSeqed *se, int pos, int keep_x_pos)
{
    seqed_setDisplayPosP(se, pos);
    seqed_redisplay_seq(se, pos, keep_x_pos); 
}

/*
 * Sets the cursor position
 */
void seqed_setCursorPos(tkSeqed *se, int pos) {
    seq_cursor_notify cn;

    /* HACK - is this correct? */
    /* if cursor is not visible, don't bother doing anything */
#ifdef HACK
    if (!seqed_onScreen(se, pos)) {
	se->cursorPos = pos;
	return;
    }
#endif    

    se->prev_pos = se->cursor->abspos;
    se->cursor->abspos = pos;
    se->cursor->job = CURSOR_MOVE;

    se->cursorPos = pos;

    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = se->cursor;
    cn.cursor->job = CURSOR_MOVE;

#ifdef DEBUG
    printf("!!!!!!seqed_setCursorPos %d abspos %d prev %d\n", pos,
	   se->cursor->abspos, se->prev_pos);
#endif
     /* seq_result_notify(se->rid, (seq_reg_data *)&info, 0);  */

/*    seq_result_notify(raster_id, (seq_reg_data *)&info2, 0); */

    seq_notify(GetSeqNum(se->seq_id), (seq_reg_data *)&cn); 
}

/*
 * returns true if base in `seq' at position `pos' is currently
 * being displayed on screen 
 */
int seqed_onScreen (tkSeqed *se, int pos)
{

#ifdef DEBUG
    printf("onscreen pos %d dispos %d \n", pos, se->displayPos);
#endif

    return (pos >= se->displayPos &&
            pos < se->displayPos + se->displayWidth);
}


void seqed_positionCursor(tkSeqed *se, int seq, int pos) {

    int screenRow,screenColumn;
        
    screenColumn = pos - se->displayPos;
    screenRow = seq;
/*
    int *seqList;
    seqList = sequencesOnScreen(xx,xx->displayPos, xx->displayWidth);
    for(screenRow=0;
	screenRow<xx->displayHeight && seqList[screenRow] != seq;
	screenRow++);
*/
#ifdef DEBUG
    printf("seqed_positionCursor\n");
#endif
    if (seqed_onScreen(se, pos)) {
	XawSheetDisplayCursor(&se->sw, True);
	XawSheetPositionCursor(&se->sw, screenColumn, screenRow);
    } else {
        XawSheetDisplayCursor(&se->sw, False);
    }

/*
    if (onScreen(se, seq, pos)) {
        int screenRow,screenColumn;
        int *seqList;
        
        screenColumn = positionInContig(xx,seq,pos) - xx->displayPos;
        seqList = sequencesOnScreen(xx,xx->displayPos, xx->displayWidth);
        for(screenRow=0;
            screenRow<xx->displayHeight && seqList[screenRow] != seq;
            screenRow++);
        XawSheetDisplayCursor(&xx->ed->sw,True);
        XawSheetPositionCursor(&xx->ed->sw,screenColumn,
			       screenRow+xx->rulerDisplayed);
    } else
        XawSheetDisplayCursor(&xx->ed->sw,False);
*/
}

/*
 * ensure that the cursor is visible on the screen
 */
void seqed_showCursor(tkSeqed *se, int seq, int pos)
{

#ifdef DEBUG
    printf("showcursor pos %d\n", pos);
#endif
    if (seqed_onScreen(se, pos)) {
        seqed_positionCursor(se, seq, pos);
   } else {
	seqed_setDisplayPos(se, pos);
    }

}

void seqed_save(tkSeqed *se,
		char *filename,
		int from,
		int to,
		int line_length)
{
    FILE *fp;
    if ((fp = fopen(filename,"w")) != NULL ) {
	seqed_write(se, fp, from, to, line_length);
	fclose(fp);
    }
}

int seqed_search(tkSeqed *se,
		 char *string,
		 int direction,
		 int strand,
		 double per_match,
		 int new_search,
		 int use_iub_code)
{
    if (new_search) {
	seqed_string_search(&se->sequence[2], se->seq_len, se->seq_name,
			    string, direction, strand, per_match, 
			    se->cursorPos, use_iub_code);
	seqed_next_string(se, direction);
    } else {
	seqed_next_string(se, direction);
    }
    return 0;
    
}

