/*
 * The Sheet widget for tk. This is simply a tk interface to the sheet.c file.
 */

#include <tk.h>
#include <string.h>
#include "tkSeqed.h"
#include "tk_defs.h"
#include "xalloc.h"
#include "tcl_utils.h"
#include "tkSheet_common.h"
#include "tkSeqedUtils.h"
#include "misc.h"
#include "seq_reg.h"
#include "dna_utils.h"
#include "seqedInterface.h"
#include "seqed_restriction_enzymes.h"
#include "seqed_search.h"

/* ---- Local defines ---- */
#define offset(field) Tk_Offset(tkSeqed, field)

/* ---- Local function prototypes ---- */
int SeqedCmd(ClientData clientData, Tcl_Interp *interp,
	     int argc, char **argv);

static int SeqedConfigure(Tcl_Interp *interp, tkSeqed *se,
			  int argc, char **argv, int flags);

static int SeqedWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv);

static void SeqedSheetExtension(ClientData clientData, int job, void *data);


/* ---- Local static data ---- */
static Tk_ConfigSpec configSpecs[] = {
#   include "tkSheet_config.h"
    {TK_CONFIG_STRING,
         "-xscrollcommand", "xScrollCommand", "ScrollCommand",
         "", offset(xScrollCmd),        TK_CONFIG_NULL_OK},
    {TK_CONFIG_STRING,
         "-yscrollcommand", "yScrollCommand", "ScrollCommand",
         "", offset(yScrollCmd),        TK_CONFIG_NULL_OK},
    {TK_CONFIG_BOOLEAN, "-ruler", "ruler", "Ruler", "1",
	 Tk_Offset(tkSeqed, rulerDisplayed), 0},
    {TK_CONFIG_BOOLEAN, "-complement", "complement", "Complement", "0",
	 Tk_Offset(tkSeqed, complementDisplayed), 0},
    {TK_CONFIG_INT, "-heightmin", "heightmin", "HeightMin", "0",
	 Tk_Offset(tkSeqed, heightmin), 0},
    {TK_CONFIG_INT, "-heightmax", "heightmax", "HeightMax", "0",
	 Tk_Offset(tkSeqed, heightmax), 0},
    {TK_CONFIG_INT, "-pos", "pos", "Pos", "1",
	 Tk_Offset(tkSeqed, cursorPos), 0},

    {TK_CONFIG_END,
	 (char *)NULL,	(char *)NULL,	(char *)NULL,	(char *) NULL,
         0,	0,	NULL},
};

#undef offset

/* ---- Globally callable functions ---- */

/*
 * Init.
 */
int Seqed_Init(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "seqed", SeqedCmd,
                      NULL, (Tcl_CmdDeleteProc *)NULL);

    return TCL_OK;
}

/* ---- Locally callable functions ---- */

/*
 * Our class command procedure; should only be used by tkEditor.c
 */
int SeqedCmd(ClientData clientData, Tcl_Interp *interp,
	     int argc, char **argv) {
    Tk_Window tkwin;
    tkSeqed *se;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " pathName ?options?\"", (char *)NULL);
        return TCL_ERROR;
    }


    /*
     * Allocate
     */
    if (NULL == (se = (tkSeqed *)xmalloc(sizeof(tkSeqed)))) {
        return TCL_ERROR;
    }

    /*
     * Call common sheet initialisation code
     */
    tkwin = SheetCmdCommon(interp, Tk_MainWindow(interp), 
			   configSpecs, (tkSheet *)se,
			   argv[1], "Seqed");
    if (NULL == tkwin) {
	xfree(se);
	return TCL_ERROR;
    }

    /*
     * some initialisation procedures
     */
 
    TKSHEET(se)->extensionData = (ClientData)se;
    TKSHEET(se)->extension = SeqedSheetExtension;

    initSeqed(se);
    /* initSeqed(se, atoi(DEF_SHEET_WIDTH), 5); */
    set_dna_lookup();

    /*
     * Register our instance of the widget class
     */
    Tcl_CreateCommand(interp, Tk_PathName(tkwin),
                      SeqedWidgetCmd, (ClientData)se,
                      (Tcl_CmdDeleteProc *)NULL);


    /*
     * And process our arguments - send them to Configure
     */
    if (SeqedConfigure(interp, se, argc-2, argv+2, 0) != TCL_OK) {
        Tk_DestroyWindow(tkwin);
        return TCL_ERROR;
    }

    /*
     * set the displayHeight and displayWidth from the -height -width values
     */
    setDimensions(se);
    Tcl_SetResult(interp, Tk_PathName(tkwin), TCL_VOLATILE);
    return TCL_OK;
}


/*
 * Configure.
 */
static int SeqedConfigure(Tcl_Interp *interp, tkSeqed *se,
			  int argc, char **argv, int flags) {

    if (TCL_ERROR == SheetConfigureCommon(interp, (tkSheet *)se,
					  argc, argv, flags))
	return TCL_ERROR;
    return TCL_OK;
}


/*
 * widget command
 */
static int SeqedWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv) {
    tkSeqed *se = (tkSeqed *)clientData;
    int result = TCL_OK;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " option ?arg arg ...?\"",
                         (char *) NULL);
        return TCL_ERROR;
    }

    Tcl_Preserve((ClientData)se);

    if (strcmp(argv[1], "configure") == 0) {
	result = SheetWidgetCmdConfig(interp, (tkSheet *)se, argc, argv);

    } else if ('x' == *argv[1] && strcmp(argv[1], "xview") == 0) {
	int pos;
	double f1, f2;
	int count;
	double d_extent;
	int displayWidth;

	d_extent = se->extent_right - se->extent_left + 1;
	pos = se->displayPos;

	if (argc == 2) {
	    /* Query */
	    f1 = pos / d_extent;
	    f2 = (pos + se->displayWidth) / d_extent;
	    vTcl_SetResult(interp, "%g %g", f1, f2);

	} else {
	    int scroll = 1;
	    switch (Tk_GetScrollInfo(interp, argc, argv, &f1, &count)) {
	    case TK_SCROLL_ERROR:
		scroll = 0;
		break;
	    case TK_SCROLL_MOVETO:
		pos = f1 * d_extent + se->extent_left;
		break;
	    case TK_SCROLL_PAGES:
		pos += count * 0.9 * se->displayWidth;
		break;
	    case TK_SCROLL_UNITS:
		pos += count;
		break;
	    }

	    if (scroll) {
		displayWidth = MIN(se->seq_len+1, se->displayWidth);

		if (pos < se->extent_left)
		    pos = se->extent_left;
		if (pos > se->extent_right + 2 - displayWidth)
		    pos = se->extent_right + 2 - displayWidth;
	    
		seqed_setDisplayPos2(se, pos, 1);
		seqed_set_h_sb_pos(se, pos); 
	    }
	}

    } else if ('y' == *argv[1] && strcmp(argv[1], "yview") == 0) {
	int type, from, to, offset = 0, count;
	double fraction, oldoffset;

	if (argc == 2) {
	    
	} else {
	    from = MIN(0, se->displayYPos);
	    to = MAX(se->displayYPos + se->displayHeight, se->num_lines);

	    type = Tk_GetScrollInfo(interp, argc, argv, &fraction, &count);
	    switch (type) {
		case TK_SCROLL_ERROR:
		    goto fail;

		case TK_SCROLL_MOVETO:
		    oldoffset = (double)(se->displayYPos -
					 MIN(0, se->displayYPos)) /
			(MAX(se->displayYPos+se->displayHeight, se->num_lines)-
			 MIN(0, se->displayYPos));
		    offset = fraction * (to - from) + from;
#ifdef DEBUG
		    printf("TK_SCROLL_MOVETO %f(%f) %d\n",
			   fraction, oldoffset, offset);
#endif		    
		    count = fraction > oldoffset ? 1 : -1;

		    break;

		case TK_SCROLL_PAGES:
		    offset = se->displayYPos + count * se->displayHeight;
		    break;

		case TK_SCROLL_UNITS:
		    offset = se->displayYPos + count;
		    break;
		}
#ifdef DEBUG
	    printf("count=%d, offset=%d, +h = %d, nlines=%d\n",
		   count, offset, offset + se->displayHeight, se->num_lines);
#endif
	    if (count > 0) {
		if (offset + se->displayHeight > se->num_lines) {
		    offset = se->num_lines - se->displayHeight;
		    if (offset < 0)
			offset = 0;
		}
	    } else {
		if (offset < 0)
		    offset = 0;
	    }
#ifdef DEBUG
	    printf("ypos = %d\n", offset);
#endif
	    se->anchor_pos -= (offset - se->displayYPos);
	    se->displayYPos = offset;
	    seqed_redisplay_seq(se, se->displayPos, 0);
	    seqed_set_v_sb_pos(se, offset);
	}
    } else if (strcmp(argv[1], "add") == 0) {
	int x, y;

	if (argc != 5) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " add x y text\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(se->interp, argv[2], &x);
	Tcl_GetInt(se->interp, argv[3], &y);

	XawSheetPutText(&se->sw, x, y, strlen(argv[4]), argv[4]);
	se->flags |= SHEET_REDRAW_TEXT;
	if (!(se->flags & SHEET_REDRAW_PENDING)) {
	    se->flags |= SHEET_REDRAW_PENDING;
	    Tcl_DoWhenIdle(SheetDisplay, (ClientData)se);
	}
    } else if ('r' == *argv[1] && strcmp(argv[1], "ruler")== 0) {
	int mode;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " ruler args\"",
			     (char *) NULL);
	    goto fail;
	}
	Tcl_GetBoolean(interp, argv[2], &mode);
	se->rulerDisplayed = mode;
	if (se->sequence) 
	    seqed_redisplay_seq(se, se->displayPos, 1);

    } else if ('c' == *argv[1] && strcmp(argv[1], "complement")== 0) {

	int mode;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " complement args\"",
			     (char *) NULL);
	    goto fail;
	}
	Tcl_GetBoolean(interp, argv[2], &mode);
	se->complementDisplayed = mode;	
	if (se->sequence) 
	    seqed_redisplay_seq(se, se->displayPos, 1);
	
    } else if ('r' == *argv[1] && strcmp(argv[1], "restriction_enzyme")== 0) {

	if (argc != 3 && argc != 6) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " restriction_enzyme command args\"",
			     (char *) NULL);
	    goto fail;
	}
	if (argc == 3) {
	    if (strcmp(argv[2], "delete") != 0) {
		Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " restriction_enzyme delete\"",
				 (char *) NULL);
		goto fail;
	    }
	    seqed_delete_renzyme(se);
	} 
	if (argc == 6) {
	    if (strcmp(argv[2], "add") != 0) {
		Tcl_AppendResult(interp, "wrong # args: should be \"",
				 argv[0], " restriction_enzyme add filename list length_of_list\"",
			     (char *) NULL);
		goto fail;
	    } 
	    seqed_add_renzyme(se, argv[3], argv[4], atoi(argv[5]));
	}
    } else if ('s' == *argv[1] && strcmp(argv[1], "seq_info")== 0) {
	int x, y;
#ifdef DEBUG
	printf("seq_info argc %d\n", argc);
#endif
	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " seq_info x y\"",
			     (char *) NULL);
	    goto fail;
	}
	
	sheet_arg_x(TKSHEET(se), argv[2], &x);
	sheet_arg_y(TKSHEET(se), argv[3], &y);
	seqedSeqInfo(se, x, y);

    } else if ('t' == *argv[1] && strcmp(argv[1], "translate")== 0) {

	if (argc != 4 ||
	    (strcmp(argv[2], "add") != 0 && strcmp(argv[2], "delete") != 0)) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " translate command type\"",
			     (char *) NULL);
	    goto fail;
	}
	
	if (strcmp(argv[2], "add") == 0) {
	    seqedTranslateAdd(interp, se, atoi(argv[3]));
	} else if (strcmp(argv[2], "delete") == 0) {
	    seqedTranslateDelete(se, atoi(argv[3]));
	}

    } else if ('t' == *argv[1] && strcmp(argv[1], "translation_mode")== 0) {
	
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " translation_mode mode\"",
			     (char *) NULL);
	    goto fail;
	}

	seqedTransMode(se, atoi(argv[2]));

    } else if ('q' == *argv[1] && strcmp(argv[1], "quit")== 0) {
	seq_reg_quit info;

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " quit\"",
			     (char *) NULL);
	    goto fail;
	}

	info.job = SEQ_QUIT;
	seq_result_notify(se->rid, (seq_reg_data *)&info, 0);

    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_left") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_left\"",
			     (char *) NULL);
	    goto fail;
	}

	seqedCursorLeft(se);

    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_right") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_right\"",
			     (char *) NULL);
	    goto fail;
	}

	seqedCursorRight(se);
    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_up") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_up\"",
			     (char *) NULL);
	    goto fail;
	}

	seqedCursorUp(se);

    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_down") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_down\"",
			     (char *) NULL);
	    goto fail;
	}

	seqedCursorDown(se);

    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_set") == 0) {
	int x,y;

	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_set xpos ypos\"",
			     (char *) NULL);
	    goto fail;
	}
	
	sheet_arg_x(TKSHEET(se), argv[2], &x);
	sheet_arg_y(TKSHEET(se), argv[3], &y);
	/* sheet_arg_y subtracts 1 from y which I don't want it to! */
	y++;
	
	if (seqedSetCursor(se, x, y)) {
	    bell();
	    Tcl_SetResult(interp, "-1", TCL_STATIC);
	} else
	    Tcl_SetResult(interp, "0", TCL_STATIC);

    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_sequence") == 0) {
	/* 
	 * if argc == 2 return the current position of the cursor 
	 * if argc == 3 set position of cursor to be argv[2] 
	 */
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_sequence ?xpos?\"",
			     (char *) NULL);
	    goto fail;
	}
	
	if (argc == 3) {

	    /* move sheet to cursor position */
	    seqed_showCursor(se, se->cursorSeq, atoi(argv[2]));

	    /* draw cursor */
	    seqed_positionCursor(se, se->cursorSeq, atoi(argv[2]));

	    /* callback to position cursor in all relevant displays */
	    seqed_setCursorPos(se, atoi(argv[2]));
	}

	vTcl_SetResult(interp, "%d", se->cursorPos);
    } else if ('s' == *argv[1] && strcmp(argv[1], "save") == 0) {
	int from, to, line_length;

	if (argc != 6) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " save filename from to line_length\"",
			     (char *) NULL);
	    goto fail;
	}
	Tcl_GetInt(interp, argv[3], &from);
	Tcl_GetInt(interp, argv[4], &to);
	Tcl_GetInt(interp, argv[5], &line_length);
	
	seqed_save(se, argv[2], from, to, line_length);

    } else if ('s' == *argv[1] && strcmp(argv[1], "search") == 0) {
	double per_match;
	int strand, direction, new_search, use_iub_code;
	char *string;
	
	if (argc != 8) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " search string direction strand percentage_match new_search\"",
			     (char *) NULL);
	    goto fail;
	}
	if (strcmp(argv[3], "+") == 0) {
	    direction = 0;
	} else {
	    direction = 1;
	}

	if (strcmp(argv[4], "+") == 0) {
	    strand = 0;
	} else {
	    strand = 1;
	}
	string = strdup(argv[2]);
	Tcl_GetDouble(interp, argv[5], &per_match);
	Tcl_GetInt(interp, argv[6], &new_search);
	Tcl_GetInt(interp, argv[7], &use_iub_code);

	seqed_search(se, string, direction, strand, per_match, 
		     new_search, use_iub_code);
	free(string);
    } else if ('s' == *argv[1] && strcmp(argv[1], "search_destroy") == 0) {

	seqed_string_search_free();
    } else {
	Tcl_AppendResult(interp, "bad option \"", argv[1], "\": must be "
			 "configure or add", NULL);
	result = TCL_ERROR;
    }

    Tcl_Release((ClientData)se);
    return result;

 fail:
    Tcl_Release((ClientData)se);
    return TCL_ERROR;
}

void delete_seqed(tkSeqed *se) {

    xfree(se->sequence);
}
/*
 * Callbacks from the sheet widget for our registered extension function.
 * See tkSheet.h for job numbers available and tkSheet.c for where they're
 * called.
 */
static void SeqedSheetExtension(ClientData clientData, int job, void *data) {
    tkSeqed *se = (ClientData)clientData;

    switch (job) {
    case SHEET_JOB_RESIZE:

	if (se->sw.columns > MAX_DISPLAY_WIDTH) {
	    int font_width = Tk_TextWidth(se->sw.font, "0", 1);

	    se->sw.columns = MAX_DISPLAY_WIDTH;
	    se->sw.width_in_pixels = se->sw.columns * font_width +
		2 * (se->sw.border_width /* + se->sw.pad_x */);
	}

	se->displayWidth = se->sw.columns;
	se->displayHeight = se->sw.rows;
#ifdef DEBUG
	printf("SHEET_JOB_RESIZE: height %d\n", se->sw.rows);
#endif
	seqed_redisplay_seq(se, se->displayPos, 1);
	break;

    case SHEET_JOB_DESTROY:
	/*
	if (ed->xx) {
	    delete_edStruct(ed->xx);
	    ed->xx = NULL;
	}
	*/ 
	if (se) {
	    delete_seqed(se);
	}
	se = NULL;
	break;
    }
}

