/*
 * The Editor tk widget.
 * We employ a little black magic here by "inheriting" the tkSheet command
 * here. This is done by providing the required tk interface functions here
 * that call the tkSheet functions and then add their own bits and pieces.
 *
 * We also reconfigure the class from Sheet to Editor thus allowing easy use
 * of bind to add the editing facilities.
 */

#include <stdlib.h>
#include <X11/Xatom.h> /* XA_STRING */

#include "tkEditor.h"
#include "edUtils.h"
#include "tkSheet_common.h"
#include "tk_defs.h"
#include "undo.h"
#include "tman_interface.h"
#include "oligo.h"
#include "xalloc.h"
#include "tcl_utils.h"
#include "edCommands.h"
#include "contigEditor.h"
#include "gap_globals.h"
#include "dstring.h"

/* ---- Local function prototypes ---- */
static int EditorConfigure(Tcl_Interp *interp, Editor *ed,
			  int argc, char **argv, int flags);

static int EditorCmd(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv);

static int EditorWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, char **argv);

static void EditorSheetExtension(ClientData clientData, int job, void *data);

/* ---- Local static data ---- */

#define offset(field) Tk_Offset(Editor, field)
static Tk_ConfigSpec configSpecs[] = {
#   include "tkSheet_config.h"
    /* Only used by the Editor */
    {TK_CONFIG_STRING,
         "-xscrollcommand", "xScrollCommand", "ScrollCommand",
         "", offset(xScrollCmd),        TK_CONFIG_NULL_OK, NULL},
    {TK_CONFIG_STRING,
         "-yscrollcommand", "yScrollCommand", "ScrollCommand",
         "", offset(yScrollCmd),        TK_CONFIG_NULL_OK, NULL},
    {TK_CONFIG_STRING,
         "-highlightcommand", "highlightCommand", "HighlightCommand",
         "", offset(highlight_cmd),     TK_CONFIG_NULL_OK, NULL},
    {TK_CONFIG_INT,
         "-max_height", "maxHeight",	"MaxHeight",	"20",
	 offset(max_height),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour0","qualColour0",	"Background",	"#494949",
	 offset(qual_bg[0]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour1","qualColour1",	"Background",	"#595959",
	 offset(qual_bg[1]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour2","qualColour2",	"Background",	"#696969",
	 offset(qual_bg[2]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour3","qualColour3",	"Background",	"#797979",
	 offset(qual_bg[3]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour4","qualColour4",	"Background",	"#898989",
	 offset(qual_bg[4]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour5","qualColour5",	"Background",	"#999999",
	 offset(qual_bg[5]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour6","qualColour6",	"Background",	"#a9a9a9",
	 offset(qual_bg[6]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour7","qualColour7",	"Background",	"#b9b9b9",
	 offset(qual_bg[7]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour8","qualColour8",	"Background",	"#c9c9c9",
	 offset(qual_bg[8]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour9","qualColour9",	"Background",	"#d9d9d9",
	 offset(qual_bg[9]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qual_fg","qualForeground",	"Foreground","#ff5050",
	 offset(qual_below),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-diff_bg","diffBackground",	"Background","green3",
	 offset(diff_bg),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-editcolour0","qualColour0",	"Background",	"red",
	 offset(edit_bg[0]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-editcolour1","qualColour1",	"Background",	"green",
	 offset(edit_bg[1]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-editcolour2","qualColour2",	"Background",	"pink",
	 offset(edit_bg[2]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-editcolour3","qualColour3",	"Background",	"lightblue",
	 offset(edit_bg[3]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour0","tmplColour0",	"Background",	"lightblue",
	 offset(tmpl_bg[0]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour1","tmplColour1",	"Background",	"lightblue",
	 offset(tmpl_bg[1]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour2","tmplColour2",	"Background",	"lightblue",
	 offset(tmpl_bg[2]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour3","tmplColour3",	"Background",	"lightblue",
	 offset(tmpl_bg[3]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour4","tmplColour4",	"Background",	"lightblue",
	 offset(tmpl_bg[4]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour5","tmplColour5",	"Background",	"lightblue",
	 offset(tmpl_bg[5]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-setcolour0","setColour0",	"Background",	"white",
	 offset(set_bg[0]),		0, NULL},
    {TK_CONFIG_COLOR,
         "-setcolour1","setColour1",	"Background",	"#c8bfff",
	 offset(set_bg[1]),		0, NULL},
    {TK_CONFIG_COLOR,
         "-setcolour2","setColour2",	"Background",	"#c8bf98",
	 offset(set_bg[2]),		0, NULL},
    {TK_CONFIG_COLOR,
         "-setcolour3","setColour3",	"Background",	"#c6f3ad",
	 offset(set_bg[3]),		0, NULL},
    {TK_CONFIG_COLOR,
         "-setcolour4","setColour4",	"Background",	"#dbc6c2",
	 offset(set_bg[4]),		0, NULL},
    {TK_CONFIG_COLOR,
         "-setcolour5","setColour5",	"Background",	"#ffffa4",
	 offset(set_bg[5]),		0, NULL},
    {TK_CONFIG_COLOR,
         "-setcolour6","setColour6",	"Background",	"#add6da",
	 offset(set_bg[6]),		0, NULL},
    {TK_CONFIG_COLOR,
         "-setcolour7","setColour7",	"Background",	"#c1edf8",
	 offset(set_bg[7]),		0, NULL},
    {TK_CONFIG_COLOR,
         "-setcolour8","setColour8",	"Background",	"#b8cb7d",
	 offset(set_bg[8]),		0, NULL},
    {TK_CONFIG_COLOR,
         "-setcolour9","setColour9",	"Background",	"#a8a8a4",
	 offset(set_bg[9]),		0, NULL},
    {TK_CONFIG_END,
	 (char *)NULL,	(char *)NULL,	(char *)NULL,	(char *) NULL,
         0,	0,	NULL},
};
#undef offset

/* ---- Global functions ---- */

/*
 * Init.
 */
int Editor_Init(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "editor", EditorCmd,
                      NULL, (Tcl_CmdDeleteProc *)NULL);

    return TCL_OK;
}


/* ---- Local functions ---- */

/*
 * Our class command procedure. This is different from usual class command
 * procs as we actually use the SheetCmd() procedure and then carefully
 * subvert their actions to our own uses.
 */
static int EditorCmd(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv) {
    Tk_Window tkwin;
    Editor *ed;
    int i;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " pathName ?options?\"", (char *)NULL);
        return TCL_ERROR;
    }

    /*
     * Allocate
     */
    if (NULL == (ed = (Editor *)xmalloc(sizeof(Editor)))) {
        return TCL_ERROR;
    }


    /*
     * Call common sheet initialisation code
     */
    tkwin = SheetCmdCommon(interp, Tk_MainWindow(interp),
			   configSpecs, (tkSheet *)ed,
			   argv[1], "Editor");
    if (NULL == tkwin) {
	xfree(ed);
	return TCL_ERROR;
    }


    /*
     * Initialise tag colours - FIXME, do in one place.
     */
    setUpColourMap(interp, Tk_MainWindow(interp));


    /*
     * Initialised rest of Editor structure.
     */
    TKSHEET(ed)->extensionData = (ClientData)ed;
    TKSHEET(ed)->extension = EditorSheetExtension;
    for (i = 0; i < 10; i++)
	ed->qual_bg[i] = NULL;
    for (i = 0; i < 4; i++)
	ed->edit_bg[i] = NULL;
    for (i = 0; i < 6; i++)
	ed->tmpl_bg[i] = NULL;
    for (i = 0; i < 10; i++)
	ed->set_bg[i] = NULL;
    ed->qual_below = NULL;
    ed->diff_bg    = NULL;
    ed->xScrollCmd = NULL;
    ed->yScrollCmd = NULL;
    ed->highlight_cmd = NULL;
    ed->max_height = 0;
    ed->xx = NULL;
    /*
     * Ideally we want gridding enabled so that "half characters" are not
     * visible, but with causes major headaches with metacity (the default
     * gnome2 window manager).
     */
    ed->grid = 0;


    /*
     * Register our instance of the widget class
     */
    Tcl_CreateCommand(interp, Tk_PathName(tkwin),
                      EditorWidgetCmd, (ClientData)ed,
                      (Tcl_CmdDeleteProc *)NULL);

    /*
     * Add a selection handler
     */
    Tk_CreateSelHandler(tkwin, XA_PRIMARY, XA_STRING, edGetSelection,
			(ClientData)ed, XA_STRING);

    /*
     * And process our arguments - send them to Configure
     */
    if (EditorConfigure(interp, ed, argc-2, argv+2, 0) != TCL_OK) {
        Tk_DestroyWindow(tkwin);
        return TCL_ERROR;
    }


    Tcl_SetResult(interp, Tk_PathName(tkwin), TCL_VOLATILE);
    return TCL_OK;
}



/*
 * Configure.
 */
static int EditorConfigure(Tcl_Interp *interp, Editor *ed,
			  int argc, char **argv, int flags) {
    if (TCL_ERROR == SheetConfigureCommon(interp, (tkSheet *)ed,
					  argc, argv, flags))
	return TCL_ERROR;

    return TCL_OK;
}


/*
 * Editor widget command
 */
static int EditorWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			       int argc, char **argv) {
    Editor *ed = (Editor *)clientData;
    int result = TCL_OK;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " option ?arg arg ...?\"",
                         (char *) NULL);
        return TCL_ERROR;
    }

    Tcl_Preserve((ClientData)TKSHEET(ed));

    if ('c' == *argv[1] && strcmp(argv[1], "configure") == 0) {
	result = SheetWidgetCmdConfig(interp, TKSHEET(ed), argc, argv);

    } else if (!ed->xx) {
	ed->xx = 0; /* Do nothing and flow through to 'success' */

    } else if ('u' == *argv[1] && strcmp(argv[1], "update_brief") == 0) {
	int x, y;

	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " update_brief xpos ypos\"",
			     (char *) NULL);
	    goto fail;
	}

	sheet_arg_x(TKSHEET(ed), argv[2], &x);
	sheet_arg_y(TKSHEET(ed), argv[3], &y);

	edSetBriefSeqStatus(ed->xx, x, y);

    } else if ('u' == *argv[1] && strcmp(argv[1], "update_brief_base") == 0) {
	int x, y;
	int mode = 1;

	if (argc < 1 || argc > 5) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " update_brief_base ?-mode2? ?xpos ypos?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc > 2 && strcmp(argv[2], "-mode2") == 0) {
	    mode = 2;
	    argc--;
	    argv++;
	}

	if (argc == 4) {
	    sheet_arg_x(TKSHEET(ed), argv[2], &x);
	    sheet_arg_y(TKSHEET(ed), argv[3], &y);
	} else {
	    x = y = -1; /* use cursorPos/Seq */
	}

	edSetBriefSeqBase(ed->xx, x, y, mode);

    } else if ('x' == *argv[1] && strcmp(argv[1], "xview") == 0) {
	int offset;
	double d_extent;
	double f1, f2;
	int type, count;
	char buf[1024];
	int scroll = 1;

	getExtents(ed->xx);
	d_extent = ed->xx->extent_right - ed->xx->extent_left + 1;
	offset = ed->xx->displayPos + ed->xx->extent_left;

	if (argc == 2) {
	    /* Query */
	    f1 = offset / d_extent;
	    f2 = (offset + ed->xx->displayWidth) / d_extent;
	    sprintf(buf, "%.20f %.20f", f1, f2);
	    Tcl_SetResult(interp, buf, TCL_VOLATILE);
	    scroll = 0;

	} else if (argc == 3) {
	    /* Old syntax */
	    Tcl_GetInt(interp, argv[2], &offset);

	} else {
	    /* New syntax */
	    type = Tk_GetScrollInfo(interp, argc, argv, &f1, &count);
	    switch (type) {
	    case TK_SCROLL_ERROR:
		scroll = 0;
		break;
	    case TK_SCROLL_MOVETO:
		offset = f1 * d_extent + ed->xx->extent_left;
		break;
	    case TK_SCROLL_PAGES:
		offset = ed->xx->displayPos +
		    count * 0.9 * ed->xx->displayWidth;
		break;
	    case TK_SCROLL_UNITS:
		offset = ed->xx->displayPos + count;
		break;
	    }
	}

	if (scroll) {
	    /* Check validity */
	    if (offset < ed->xx->extent_left)
		offset = ed->xx->extent_left;
	    if (offset > ed->xx->extent_right + 2 - ed->xx->displayWidth)
		offset = ed->xx->extent_right + 2 - ed->xx->displayWidth;

	    setDisplayPos2(ed->xx, offset);
	    ed_set_slider_pos(ed->xx, offset);
	}

    } else if ('y' == *argv[1] && strcmp(argv[1], "yview") == 0) {
	double f1, f2;
	int scroll = 1;
	int offset = ed->xx->displayYPos * ed->xx->lines_per_seq;
	int type;
	int count;

	if (argc == 2) {
	    /* Query */
	    f1 = ed->xx->displayYPos * ed->xx->lines_per_seq /
		(double)ed->xx->totalHeight;
	    f2 = (ed->xx->displayYPos * ed->xx->lines_per_seq +
		  ed->xx->displayHeight) /
		(double)ed->xx->totalHeight;
	    vTcl_SetResult(interp, "%g %g", f1, f2);
	    scroll = 0;

	} else if (argc == 3) {
	    /* Old syntax */
	    Tcl_GetInt(interp, argv[2], &offset);

	} else {
	    /* New syntax */
	    type = Tk_GetScrollInfo(interp, argc, argv, &f1, &count);
	    switch (type) {
	    case TK_SCROLL_ERROR:
		scroll = 0;
		break;
	    case TK_SCROLL_MOVETO:
		offset = f1 * ed->xx->totalHeight;
		break;
	    case TK_SCROLL_PAGES:
		offset += count * 0.9 * ed->xx->displayHeight;
		break;
	    case TK_SCROLL_UNITS:
		offset += count * ed->xx->lines_per_seq;
		break;
	    }
	}

	offset /= ed->xx->lines_per_seq;

	if (scroll) {
	    /* Check validity */
	    if (offset < 0)
		offset = 0;
	    if (offset > ed->xx->totalHeight - ed->xx->displayHeight)
		offset = ed->xx->totalHeight - ed->xx->displayHeight;

	    setDisplayYPos(ed->xx, offset);
	    ed_set_yslider_pos(ed->xx, offset);
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_reveal") == 0) {
	int cutoff;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_reveal ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2)
	    edSetRevealCutoffs(ed->xx, -1);
	else {
	    Tcl_GetInt(interp, argv[2], &cutoff);
	    edSetRevealCutoffs(ed->xx, cutoff);
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_grouping") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_grouping ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3) {
	    int v;
	    Tcl_GetInt(interp, argv[2], &v);
	    if (v >= 0)
		ed->xx->group_mode = v;
	    else
		ed->xx->group_mode ^= -v;

	    ed->xx->refresh_flags |= ED_DISP_ALL;
	    redisplaySequences(ed->xx, 0);
	}

	vTcl_SetResult(interp, "%d", ed->xx->group_mode);

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_insert") == 0) {
	int val;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_insert ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2)
	    edSetInsertMode(ed->xx, -1);
	else {
	    Tcl_GetInt(interp, argv[2], &val);
	    edSetInsertMode(ed->xx, val);
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_unpadded_ruler") == 0) {
	int val;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_unpadded_ruler ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2)
	    edSetRulerMode(ed->xx, -1);
	else {
	    Tcl_GetInt(interp, argv[2], &val);
	    edSetRulerMode(ed->xx, val);
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "show_template_names") == 0) {
	int val;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " show_template_names ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3) {
	    Tcl_GetInt(interp, argv[2], &val);
	    if (-1 == val)
		val = edGetTemplateNameMode(ed->xx) ^ 1;
	    edSetTemplateNameMode(ed->xx, val);
	}
	vTcl_SetResult(interp, "%d", edGetTemplateNameMode(ed->xx));

    } else if ('s' == *argv[1] && strcmp(argv[1], "superedit") == 0) {
	int s_argc;
	char **s_argv;
	int mode;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " superedit modes\"",
			     (char *) NULL);
	    goto fail;
	}

	if (Tcl_SplitList(interp, argv[2], &s_argc, &s_argv) != TCL_OK)
	    goto fail;

	if (s_argc != 10)
	    goto fail;

	mode = 0;
	mode |= atoi(s_argv[0]) * SUPEREDIT_INS_READ;
	mode |= atoi(s_argv[1]) * SUPEREDIT_DEL_READ;
	mode |= atoi(s_argv[2]) * SUPEREDIT_INS_ANY_CON;
	mode |= atoi(s_argv[3]) * SUPEREDIT_DEL_DASH_CON;
	mode |= atoi(s_argv[4]) * SUPEREDIT_DEL_ANY_CON;
	mode |= atoi(s_argv[5]) * SUPEREDIT_REPLACE_CON;
	mode |= atoi(s_argv[6]) * SUPEREDIT_SHIFT_READ;
	mode |= atoi(s_argv[7]) * SUPEREDIT_TRANSPOSE_ANY;
	mode |= atoi(s_argv[8]) * SUPEREDIT_UPPERCASE;
	mode |= atoi(s_argv[9]) * SUPEREDIT_MODIFY_CONF;

	ed->xx->super_edit = mode;

	Tcl_Free((char *)s_argv);

    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_left") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_left\"",
			     (char *) NULL);
	    goto fail;
	}

	edCursorLeft(ed->xx);

    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_right") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_right\"",
			     (char *) NULL);
	    goto fail;
	}

	edCursorRight(ed->xx);
    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_up") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_up\"",
			     (char *) NULL);
	    goto fail;
	}

	edCursorUp(ed->xx);

    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_down") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_down\"",
			     (char *) NULL);
	    goto fail;
	}

	edCursorDown(ed->xx);

    } else if ('t' == *argv[1] && strcmp(argv[1], "transpose_left") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " transpose_left\"",
			     (char *) NULL);
	    goto fail;
	}

	if (edTransposeLeft(ed->xx, 1))
	    bell();

    } else if ('t' == *argv[1] && strcmp(argv[1], "transpose_right") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " transpose_right\"",
			     (char *) NULL);
	    goto fail;
	}

	if (edTransposeRight(ed->xx, 1))
	    bell();

    } else if ('d' == *argv[1] && strcmp(argv[1], "delete_key") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " delete_key\"",
			     (char *) NULL);
	    goto fail;
	}

	if (edKeyDelete(ed->xx))
	    bell();

    } else if ('d' == *argv[1] && strcmp(argv[1], "delete_left_key") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " delete_left_key\"",
			     (char *) NULL);
	    goto fail;
	}

	if (edKeyDeleteLeft(ed->xx))
	    bell();

    } else if ('e' == *argv[1] && strcmp(argv[1], "edit_key") == 0) {
	char key, nomove = 0;

	if (argc != 3 && argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " edit_key ?-nomove? character\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 4) {
	    nomove = 1;
	    key = argv[3][0];
	} else {
	    key = argv[2][0];
	}

	if (edKeyPress(ed->xx, key, nomove))
	    bell();

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_confidence") == 0) {
	int val;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_confidence value\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &val);
	if (100 == val)
	    edConf100(ed->xx);
	else if (0 == val)
	    edConf0(ed->xx);
	else
	    bell();

    } else if ('i' == *argv[1] && strcmp(argv[1], "increment_confidence")==0) {
	int val;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " increment_confidence amount\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &val);
	edConfIncr(ed->xx, val);

    } else if ('r' == *argv[1] && strcmp(argv[1], "read_start") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " read_start\"",
			     (char *) NULL);
	    goto fail;
	}

	edStartRead(ed->xx);

    } else if ('r' == *argv[1] && strcmp(argv[1], "read_end") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " read_end\"",
			     (char *) NULL);
	    goto fail;
	}

	edEndRead(ed->xx);

    } else if ('r' == *argv[1] && strcmp(argv[1], "read_start2") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " read_start2\"",
			     (char *) NULL);
	    goto fail;
	}

	edStartRead2(ed->xx);

    } else if ('r' == *argv[1] && strcmp(argv[1], "read_end2") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " read_end2\"",
			     (char *) NULL);
	    goto fail;
	}

	edEndRead2(ed->xx);

    } else if ('c' == *argv[1] && strcmp(argv[1], "contig_start") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " contig_start\"",
			     (char *) NULL);
	    goto fail;
	}

	edStartContig(ed->xx);

    } else if ('c' == *argv[1] && strcmp(argv[1], "contig_end") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " contig_end\"",
			     (char *) NULL);
	    goto fail;
	}

	edEndContig(ed->xx);

    } else if ('e' == *argv[1] && strcmp(argv[1], "extend_left") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " extend_left\"",
			     (char *) NULL);
	    goto fail;
	}

	edExtendLeft(ed->xx);

    } else if ('e' == *argv[1] && strcmp(argv[1], "extend_right") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " extend_right\"",
			     (char *) NULL);
	    goto fail;
	}

	edExtendRight(ed->xx);

    } else if ('z' == *argv[1] && strcmp(argv[1], "zap_left") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " zap_left\"",
			     (char *) NULL);
	    goto fail;
	}

	edZapLeft(ed->xx);

    } else if ('z' == *argv[1] && strcmp(argv[1], "zap_right") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " zap_right\"",
			     (char *) NULL);
	    goto fail;
	}

	edZapRight(ed->xx);

    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_set") == 0) {
	int x,y;

	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_set xpos ypos\"",
			     (char *) NULL);
	    goto fail;
	}

	sheet_arg_x(TKSHEET(ed), argv[2], &x);
	sheet_arg_y(TKSHEET(ed), argv[3], &y);

	if (edSetCursor(ed->xx, x, y)) {
	    bell();
	    Tcl_SetResult(interp, "-1", TCL_STATIC);
	} else
	    Tcl_SetResult(interp, "0", TCL_STATIC);

    } else if ('c' == *argv[1] && strcmp(argv[1], "cursor_consensus") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " cursor_consensus ?xpos?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3)
	    edSetCursorConsensus(ed->xx, atoi(argv[2]));

	vTcl_SetResult(interp, "%d",
		       DB_RelPos(ed->xx, ed->xx->cursorSeq)-1 +
		       ed->xx->cursorPos);

    } else if ('u' == *argv[1] && strcmp(argv[1], "undo") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " undo\"",
			     (char *) NULL);
	    goto fail;
	}

	undoLastCommand(ed->xx);

    } else if ('s' == *argv[1] && strcmp(argv[1], "select") == 0) {
	if (argc < 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " select option ?arg?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (strcmp(argv[2], "clear") == 0) {
	    if (argc != 3) {
		Tcl_AppendResult(interp, "wrong # args: should be \"",
				 argv[0], " select clear\"",
				 (char *) NULL);
		goto fail;
	    }

	    edSelectClear(ed->xx);

	} else {
	    int arg;

	    if (argc != 4) {
		Tcl_AppendResult(interp, "wrong # args: should be \"",
				 argv[0], " select option arg\"",
				 (char *) NULL);
		goto fail;
	    }

	    sheet_arg_x(TKSHEET(ed), argv[3], &arg);

	    if (strcmp(argv[2], "from") == 0) {
		edSelectFrom(ed->xx, arg);

	    } else if (strcmp(argv[2], "to") == 0) {
		edSelectTo(ed->xx, arg);

	    } else if (strcmp(argv[2], "adjust") == 0) {
		edSelectTo(ed->xx, arg);

	    } else {
		Tcl_AppendResult(interp, "bad select option \"", argv[2],
				 "\": must be adjust, clear, from, or to",
				 (char *) NULL);
	    }
	}

    } else if ('d' == *argv[1] && strcmp(argv[1], "delete_anno") == 0) {
	tagStruct *t = NULL;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " delete_anno ?tagptr?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3)
	    sscanf(argv[2], "%p", &t);

	deleteAnnotation(ed->xx, t);

    } else if ('c' == *argv[1] && strcmp(argv[1], "create_anno") == 0) {
	char *res;

	if (argc != 2 && argc != 5) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " create_anno ?type text strand?\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_ResetResult(interp);
	if (argc == 5) {
	    /* Create a predefined tag at the current position */
	    saveAnnotation(ed->xx, argv[2], argv[3], atoi(argv[4]));
	} else {
	    /* Invoke tag editor */
	    if (NULL != (res = createAnnotation(ed->xx)))
		Tcl_SetResult(interp, res, TCL_VOLATILE);
	}

    } else if ('e' == *argv[1] && strcmp(argv[1], "edit_anno") == 0) {
	tagStruct *t = NULL;
	char *res;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " edit_anno ?tagptr?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3)
	    sscanf(argv[2], "%p", &t);

	Tcl_ResetResult(interp);
	if (NULL != (res = editAnnotation(ed->xx, t)))
	    Tcl_SetResult(interp, res, TCL_VOLATILE);

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_anno") == 0) {
	tagStruct *t = NULL;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_anno tagptr\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3)
	    sscanf(argv[2], "%p", &t);

	force_comment(DBI_io(ed->xx), t);
	Tcl_SetResult(interp, t->newcomment, TCL_VOLATILE);

    } else if ('l' == *argv[1] && strcmp(argv[1], "list_anno") == 0) {
	dstring_t *ds;

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " list_anno\"",
			     (char *) NULL);
	    goto fail;
	}

	if (NULL != (ds = listAnnotation(ed->xx))) {
	    Tcl_SetResult(interp, dstring_str(ds), TCL_VOLATILE);
	    dstring_destroy(ds);
	}

    } else if ('c' == *argv[1] && strcmp(argv[1], "create_tmp_anno") == 0) {
	if (argc != 8) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " create_anno seq pos len type text strand\"",
			     (char *) NULL);
	    goto fail;
	}

	/* Create a predefined tag at the current position */
	createTmpAnnotation(ed->xx,
			    atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
			    argv[5], argv[6], atoi(argv[7]));

    } else if ('s' == *argv[1] && strcmp(argv[1], "save") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " save\"",
			     (char *) NULL);
	    goto fail;
	}

	saveDB(ed->xx, DBI_io(ed->xx), 0 /* autosave */, 1 /* notify others*/);
	redisplaySequences(ed->xx, 1);

    } else if ('a' == *argv[1] && strcmp(argv[1], "auto_save") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " auto_save ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    /* toggle */
	    ed->xx->auto_save ^= 1;
	} else {
	    /* set */
	    Tcl_GetInt(interp, argv[2], &ed->xx->auto_save);
	}

	saveDB(ed->xx, DBI_io(ed->xx), 0, 1);

    } else if ('s' == *argv[1] && strcmp(argv[1], "store_undo") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " store_undo ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    /* Report */
	    vTcl_SetResult(interp, "%d", DBI_store_undo(ed->xx));
	} else {
	    int i = atoi(argv[2]);
	    if (i == -1) { /* toggle */
		DBI_store_undo(ed->xx) ^= 1;
	    } else { /* set */
		DBI_store_undo(ed->xx) = i;
	    }

	    if (DBI_store_undo(ed->xx))
		startUndo(ed->xx);
	    else
		stopUndo(ed->xx);
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "search") == 0) {
	char *value = "";
	int found;

	if (argc < 5 || argc > 6) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " search direction strand type ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 6)
	    value = argv[5];

	found = edDoSearch(ed->xx, (strcmp("backward", argv[2]) != 0),
			   *argv[3], argv[4], value);

	vTcl_SetResult(interp, "%d", found);

    } else if ('s' == *argv[1] && strcmp(argv[1], "show_differences") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " show_differences ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    vTcl_SetResult(interp, "%d", ed->xx->showDifferences);
	} else {
	    Tcl_GetInt(EDINTERP(ed), argv[2], &ed->xx->showDifferences);
	    ed->xx->refresh_flags |= ED_DISP_READS;
	    redisplaySequences(ed->xx, 0);
	}

    } else if ('r' == *argv[1] && strcmp(argv[1], "replace_confidence") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " replace_confidence ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    vTcl_SetResult(interp, "%d", ed->xx->default_conf_r);
	} else {
	    Tcl_GetInt(EDINTERP(ed), argv[2], &ed->xx->default_conf_r);
	}

    } else if ('i' == *argv[1] && strcmp(argv[1], "insert_confidence") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " insert_confidence ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    vTcl_SetResult(interp, "%d", ed->xx->default_conf_n);
	} else {
	    Tcl_GetInt(EDINTERP(ed), argv[2], &ed->xx->default_conf_n);
	}

    } else if ('c' == *argv[1] && strcmp(argv[1], "compare_strands") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " compare_strands ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2)
	    ed->xx->compare_strands ^= 1;
	else
	    Tcl_GetInt(EDINTERP(ed), argv[2], &ed->xx->compare_strands);

	ed->xx->refresh_flags |= ED_DISP_CONS | ED_DISP_STATUS;
	redisplaySequences(ed->xx, 0);

    } else if ('a' == *argv[1] && strcmp(argv[1], "autodisplay_traces") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " autodisplay_traces ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2)
	    ed->xx->display_traces ^= 1;
	else
	    Tcl_GetInt(EDINTERP(ed), argv[2], &ed->xx->display_traces);

    } else if ('a' == *argv[1] && strcmp(argv[1], "autodiff_traces") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " autodiff_traces ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2)
	    ed->xx->diff_traces ^= 1;
	else
	    Tcl_GetInt(EDINTERP(ed), argv[2], &ed->xx->diff_traces);

    } else if ('r' == *argv[1] && strcmp(argv[1], "read_pair_traces") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " read_pair_traces ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2)
	    ed->xx->read_pair_traces ^= 1;
	else
	    Tcl_GetInt(EDINTERP(ed), argv[2], &ed->xx->read_pair_traces);

    } else if ('d' == *argv[1] && strcmp(argv[1], "dump_contig") == 0) {
	int left, right, llength, nwidth;

	if (argc != 7) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " dump_contig file from to "
			     "line_len name_width\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[3], &left);
	Tcl_GetInt(interp, argv[4], &right);
	Tcl_GetInt(interp, argv[5], &llength);
	Tcl_GetInt(interp, argv[6], &nwidth);

	dumpContig(ed->xx, argv[2], left, right, llength, nwidth);

    } else if ('c' == *argv[1] && strcmp(argv[1], "consensus_trace") == 0) {
	int left, right, strand, matching;

	if (argc != 7) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " consensus_trace file from to strand matching\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[3], &left);
	Tcl_GetInt(interp, argv[4], &right);
	Tcl_GetInt(interp, argv[5], &strand);
	Tcl_GetInt(interp, argv[6], &matching);

	save_consensus_trace(ed->xx, argv[2], left, right, strand, matching);

    } else if ('a' == *argv[1] && strcmp(argv[1], "align") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " align\"",
			     (char *) NULL);
	    goto fail;
	}

	align_read(ed->xx);

    } else if ('s' == *argv[1] && strcmp(argv[1], "select_oligo") == 0) {
	/*
	 * select_oligos generate sense fwd back avg_readlen
	 * select_oligos next
	 * select_oligos prev
	 * select_oligos accept template
	 * select_oligos quit
	 *
	 */
	if (argc == 8 && strcmp(argv[2], "generate") == 0) {
	    vTcl_SetResult(interp, "%d",
			   edSelectOligoGenerate(ed->xx,
						 atoi(argv[3]) /* sense */,
						 atoi(argv[4]) /* fwd */,
						 atoi(argv[5]) /* bwd */,
						 atoi(argv[6]) /* readlen */,
						      argv[7]  /* p3_def */));
	} else if (argc == 3 && strcmp(argv[2], "next") == 0) {
	    char *ret = edSelectOligoNext(ed->xx);

	    if (ret) {
		Tcl_SetResult(interp, ret, TCL_VOLATILE);
		free(ret);
	    }
	} else if (argc == 3 && strcmp(argv[2], "prev") == 0) {
	    char *ret = edSelectOligoPrev(ed->xx);

	    if (ret) {
		Tcl_SetResult(interp, ret, TCL_VOLATILE);
		free(ret);
	    }
	} else if (argc == 4 && strcmp(argv[2], "accept") == 0) {
	    Tcl_SetResult(interp, edSelectOligoAccept(ed->xx, argv[3]),
			  TCL_STATIC);
	} else if (argc == 3 && strcmp(argv[2], "quit") == 0) {
	    edSelectOligoQuit(ed->xx);
	} else {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " select_oligo command ?options?\"",
			     (char *) NULL);
	    goto fail;
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "strip_pads") == 0) {
	int sp_consensus_mode;
	float sp_consensus_cutoff;

	if (argc != 4 && argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " strip_pads ?consensus_mode "
			     "consensus_cutoff(fraction)?\"", (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    sp_consensus_mode = CONSENSUS_MODE_FREQ;
	    sp_consensus_cutoff = 1.0;
	} else {
	    double tmp;
	    Tcl_GetInt(interp, argv[2], &sp_consensus_mode);
	    Tcl_GetDouble(interp, argv[3], &tmp);
	    sp_consensus_cutoff = tmp;
	}
	strip_pads(ed->xx, sp_consensus_mode, sp_consensus_cutoff);

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_number") == 0) {
	char buf[10];
	int x, y, num;

	if (argc != 2 && argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_number ?xpos ypos?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 4) {
	    sheet_arg_x(TKSHEET(ed), argv[2], &x);
	    sheet_arg_y(TKSHEET(ed), argv[3], &y);

	    if (-1 != (num = edGetGelNumber(ed->xx, x, y))) {
		sprintf(buf, "%d", num);
		Tcl_AppendResult(interp, buf, NULL);
	    } /* otherwise return a blank */
	} else {
	    sprintf(buf, "%d", ed->xx->cursorSeq);
	    Tcl_AppendResult(interp, buf, NULL);
	}

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_read_number") == 0) {
	char buf[10];
	int x, y, num;

	if (argc != 2 && argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_read_number ?xpos ypos?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 4) {
	    sheet_arg_x(TKSHEET(ed), argv[2], &x);
	    sheet_arg_y(TKSHEET(ed), argv[3], &y);

	    if (-1 != (num = DB_Number(ed->xx,edGetGelNumber(ed->xx, x, y)))) {
		sprintf(buf, "%d", num);
		Tcl_AppendResult(interp, buf, NULL);
	    } /* otherwise return a blank */
	} else {
	    sprintf(buf, "%d", DB_Number(ed->xx, ed->xx->cursorSeq));
	    Tcl_AppendResult(interp, buf, NULL);
	}

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_name") == 0) {
	char *name;
	int num;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_name ?gel_number?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3)
	    Tcl_GetInt(interp, argv[2], &num);
	else
	    num = ed->xx->cursorSeq;

	Tcl_ResetResult(interp);
	if (NULL != (name = edGetGelName(ed->xx, num))) {
	    Tcl_AppendResult(interp, name, NULL);
	} /* otherwise return a blank */

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_names_to_right") == 0) {
	int num;
	dstring_t *ds = NULL;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_name ?gel_number?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3)
	    Tcl_GetInt(interp, argv[2], &num);
	else
	    num = ed->xx->cursorSeq;

	if (NULL != (ds = edGetGelNamesToRight(ed->xx, num))) {
	    Tcl_SetResult(interp, dstring_str(ds), TCL_VOLATILE);
	    dstring_destroy(ds);
	} else {
	    /* otherwise return a blank */
	    Tcl_ResetResult(interp);
	}

    } else if ('f' == *argv[1] && strcmp(argv[1], "find_read") == 0) {
	int num;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " find_read ?identifier?\"",
			     (char *) NULL);
	    goto fail;
	}

	num = get_gel_num(DBI_io(ed->xx), argv[2], GGN_ID);
	num = rnum_to_edseq(ed->xx, num);
	vTcl_SetResult(interp, "%d", num);

    } else if ('i' == *argv[1] && strcmp(argv[1], "invoke_trace") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " invoke_trace\"",
			     (char *) NULL);
	    goto fail;
	}

	edInvokeTrace(ed->xx);

    } else if ('d' == *argv[1] && strcmp(argv[1], "delete_trace") == 0) {
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " delete_trace path\"",
			     (char *) NULL);
	    goto fail;
	}

	deleteTrace(ed->xx, argv[2]);

    } else if ('d' == *argv[1] && strcmp(argv[1], "diff_trace") == 0) {
	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " delete_trace path1 path2\"",
			     (char *) NULL);
	    goto fail;
	}

	diffTrace(ed->xx, argv[2], argv[3]);

    } else if ('t' == *argv[1] && strcmp(argv[1], "trace_redraw") == 0) {
	/* repositionTraces(ed->xx); */

    } else if ('t' == *argv[1] && strcmp(argv[1], "trace_config") == 0) {
	if (argc != 6 && argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " trace_config ?cons_match cons_select "
			     "algorithm yscale?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    vTcl_SetResult(interp, "%d %d %d %d",
			   ed->xx->compare_trace_match,
			   ed->xx->compare_trace_select,
			   ed->xx->compare_trace_algorithm,
			   ed->xx->compare_trace_yscale);
	} else {
	    ed->xx->compare_trace_match     = atoi(argv[2]);
	    ed->xx->compare_trace_select    = atoi(argv[3]);
	    ed->xx->compare_trace_algorithm = atoi(argv[4]);
	    ed->xx->compare_trace_yscale    = atoi(argv[5]);
	}

    } else if ('t' == *argv[1] && strcmp(argv[1], "trace_comparator") == 0) {
	int num;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " trace_comparator ?identifier?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3) {
	    if (strcmp(argv[2], "0") == 0) {
		num = 0; /* Consensus */
	    } else {
		num = get_gel_num(DBI_io(ed->xx), argv[2], GGN_ID);
		num = rnum_to_edseq(ed->xx, num);
	    }

	    edSetTraceComparator(ed->xx, num);
	} else {
	    edSetTraceComparator(ed->xx, -1);
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_trace_lock") == 0) {
	int lock;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_trace_lock ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    tman_set_lock(ed->xx, !tman_get_lock(ed->xx));
	} else {
	    Tcl_GetBoolean(interp, argv[2], &lock);
	    tman_set_lock(ed->xx, lock);
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_ccutoff") == 0) {
	int cutoff;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_ccutoff ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    char buf[10];
	    sprintf(buf, "%d", (int)(100 * ed->xx->con_cut + 0.1));
	    Tcl_AppendResult(interp, buf, NULL);
	} else {
	    Tcl_GetInt(interp, argv[2], &cutoff);
	    edSetCCutoff(ed->xx, cutoff);
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_qcutoff") == 0) {
	int cutoff;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_qcutoff ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    char buf[10];
	    sprintf(buf, "%d", ed->xx->qual_cut);
	    Tcl_AppendResult(interp, buf, NULL);
	} else {
	    Tcl_GetInt(interp, argv[2], &cutoff);
	    edSetQCutoff(ed->xx, cutoff);
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_disagreement_cutoff") == 0) {
	int cutoff;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_disagreement_cutoff ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    char buf[10];
	    sprintf(buf, "%d", ed->xx->diff_qual);
	    Tcl_AppendResult(interp, buf, NULL);
	} else {
	    Tcl_GetInt(interp, argv[2], &cutoff);
	    edSetDifferenceQuality(ed->xx, cutoff);
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_consensus_mode") == 0) {
	int mode;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_consensus_mode ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 2) {
	    char buf[10];
	    sprintf(buf, "%d", ed->xx->consensus_mode);
	    Tcl_AppendResult(interp, buf, NULL);
	} else {
	    Tcl_GetInt(interp, argv[2], &mode);
	    edSetCMode(ed->xx, mode);
	}

    } else if ('j' == *argv[1] && strcmp(argv[1], "join_lock") == 0) {
	int lock;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " join_lock ?value?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3) {
	    Tcl_GetInt(interp, argv[2], &lock);
	    edSetJoinLock(ed->xx, lock);
	} else {
	    edSetJoinLock(ed->xx, !editorLocked(ed->xx));
	}

    } else if ('j' == *argv[1] && strcmp(argv[1], "join_mode") == 0) {
	char buf[10];

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " join_mode\"",
			     (char *) NULL);
	    goto fail;
	}

	sprintf(buf, "%d",
		inJoinMode(ed->xx) &&
		DBI_contigNum(ed->xx->link->xx[0]) !=
		DBI_contigNum(ed->xx->link->xx[1]));
	Tcl_AppendResult(interp, buf, NULL);

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_extents") == 0) {
	int left,right;

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_extents\"",
			     (char *) NULL);
	    goto fail;
	}
	extents(ed->xx, &left, &right);

	vTcl_SetResult(interp, "%d %d", left, right);

	
    } else if ('g' == *argv[1] && strcmp(argv[1], "get_contig_num") == 0) {
	vTcl_SetResult(interp, "%d", DBI_contigNum(ed->xx));
 
    } else if ('j' == *argv[1] && strcmp(argv[1], "join_percentage") == 0) {
	char buf[100];
	int overlapLength, wingeCount;
	float perc = -1.0;
	int tgood, tbad;

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " join_percentage\"",
			     (char *) NULL);
	    goto fail;
	}

	if (ed->xx->link) {
	    countDisagreements(ed->xx->link->xx, &overlapLength, &wingeCount,
			       &tgood, &tbad);
	    if (overlapLength > 0)
		perc = (float)(100 * wingeCount) / (float)overlapLength;
	}

	sprintf(buf, "%.20f %d %d", perc, tgood, tbad);
	Tcl_AppendResult(interp, buf, NULL);

    } else if ('j' == *argv[1] && strcmp(argv[1], "join_align") == 0) {
	int fixed_left = 0, fixed_right = 0;
	if (argc != 2 && argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " join_align ?fixed-left fixed-right?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 4) {
	    Tcl_GetInt(interp, argv[2], &fixed_left);
	    Tcl_GetInt(interp, argv[3], &fixed_right);
	}

	edJoinAlign(ed->xx, fixed_left, fixed_right);

    } else if ('j' == *argv[1] && strcmp(argv[1], "join") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " join\"",
			     (char *) NULL);
	    goto fail;
	}

	edJoin(ed->xx);

    } else if ('e' == *argv[1] && strcmp(argv[1], "edits_made") == 0) {
	char buf[10];

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " edits_made\"",
			     (char *) NULL);
	    goto fail;
	}

	sprintf(buf, "%d", editsMade(ed->xx));
	Tcl_AppendResult(interp, buf, NULL);

    } else if ('q' == *argv[1] && strcmp(argv[1], "quit") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " quit\"",
			     (char *) NULL);
	    goto fail;
	}

	db_callback_tk(ed->xx, DBCALL_QUIT, 0, 0, NULL);

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_displayed_annos")== 0) {
	if (argc < 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_displayed_annos ?anno_type? ...\"",
			     (char *) NULL);
	    goto fail;
	}

	edSetActiveAnnos(ed->xx, argc-2, &argv[2]);

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_displayed_annos")== 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_displayed_annos\"",
			     (char *) NULL);
	    goto fail;
	}

	vTcl_SetResult(interp, "%s", edGetActiveAnnos(ed->xx));

    } else if ('s' == *argv[1] && strcmp(argv[1], "status")== 0) {
	if (argc != 4 ||
	    (strcmp(argv[2], "add") != 0 && strcmp(argv[2], "delete") != 0)) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " status command type\"",
			     (char *) NULL);
	    goto fail;
	}


	if (strcmp(argv[2], "add") == 0) {
	    edStatusAdd(ed->xx, atoi(argv[3]));
	} else if (strcmp(argv[2], "delete") == 0) {
	    edStatusDelete(ed->xx, atoi(argv[3]));
	}

    } else if ('t' == *argv[1] && strcmp(argv[1], "translation_mode")== 0) {
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " translation_mode mode\"",
			     (char *) NULL);
	    goto fail;
	}

	edStatusTransMode(ed->xx, atoi(argv[2]));

    } else if ('s' == *argv[1] &&
	       (strcmp(argv[1], "show_quality")== 0 ||
		strcmp(argv[1], "show_reading_quality")== 0)) {
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " show_quality value\"",
			     (char *) NULL);
	    goto fail;
	}

	edShowQuality(ed->xx, atoi(argv[2]));

    } else if ('s' == *argv[1] &&
	       strcmp(argv[1], "show_consensus_quality")== 0) {
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " show_consensus_quality value\"",
			     (char *) NULL);
	    goto fail;
	}

	edShowCQuality(ed->xx, atoi(argv[2]));

    } else if ('s' == *argv[1] && strcmp(argv[1], "show_edits")== 0) {
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " show_edits value\"",
			     (char *) NULL);
	    goto fail;
	}

	edShowEdits(ed->xx, atoi(argv[2]));

    } else if ('h' == *argv[1] && strcmp(argv[1], "hide_read")== 0) {
	int num;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " hide_read ?seq_number?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3)
	    Tcl_GetInt(interp, argv[2], &num);
	else
	    num = ed->xx->cursorSeq;

	edHideRead(ed->xx, num, argc != 3);

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_hidden_reads")== 0) {
        int *reads, i;
	extern int *edGetHiddenReads(EdStruct *xx); /* FIXME */

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_hidden_reads\"",
			     (char *) NULL);
	    goto fail;
	}

	if (NULL != (reads = edGetHiddenReads(ed->xx))) {
	    i = 0;
	    while (reads[i]) {
	        char buf[20];

		sprintf(buf, "#%d ", reads[i++]);
	        Tcl_AppendResult(interp, buf, NULL);
	    }

	    xfree(reads);
	}

    } else if ('i' == *argv[1] && strcmp(argv[1], "io")== 0) {
        char buf[20];

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " io\"",
			     (char *) NULL);
	    goto fail;
	}

	sprintf(buf, "%d", *handle_io(DBI_io(ed->xx)));
	Tcl_AppendResult(interp, buf, NULL);

    } else if ('w' == *argv[1] && strcmp(argv[1], "write_mode") == 0) {
        if (ed->xx) {
            if (argc > 2) {
                int i = atoi(argv[2]);
                if (i == -1)
                    DBI_flags(ed->xx) ^= DB_ACCESS_UPDATE;
                else if (i == 1)
                    DBI_flags(ed->xx) |= DB_ACCESS_UPDATE;
                else
                    DBI_flags(ed->xx) &= ~DB_ACCESS_UPDATE;
                redisplaySequences(ed->xx, 1);
            }
            vTcl_SetResult(interp, "%d",
                           (int)(DBI_flags(ed->xx) & DB_ACCESS_UPDATE));
        }
    } else if ('l' == *argv[1] && strcmp(argv[1], "list_confidence") == 0) {
	int left, right, info_only;

	if (argc != 5) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " list_confidence from to info_only\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &left);
	Tcl_GetInt(interp, argv[3], &right);
	Tcl_GetInt(interp, argv[4], &info_only);

	edListConfidence(ed->xx, left, right, info_only);

    } else if ('t' == *argv[1] && strcmp(argv[1], "trace_scroll") == 0) {
	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " trace_scroll path scroll_cmd\"",
			     (char *) NULL);
	    goto fail;
	}

	edScrollTraces(ed->xx, argv[2], argv[3]);

    } else if ('n' == *argv[1] && strcmp(argv[1], "number_of_views") == 0) {
	char buf[20];
	int nviews = 0;
	int i;

	if (argc != 2) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " number_of_views\"",
			     (char *) NULL);
	    goto fail;
	}

	for (i = 0; i < MAX_DISP_PROCS; i++) {
	    if (DBI_dispFunc(ed->xx)[i]) {
		nviews++;
	    }
	}

	sprintf(buf, "%d", nviews);
	Tcl_AppendResult(interp, buf, NULL);

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_reference_seq") == 0) {
	int seq, length, offset;

	if (argc != 5 && argc != 3) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_reference_seq seq ?length offset?\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &seq);
	if (argc == 5) {
	    /* Set */
	    Tcl_GetInt(interp, argv[3], &length);
	    Tcl_GetInt(interp, argv[4], &offset);

	    vTcl_SetResult(interp, "%d",
			   set_reference_seq(ed->xx, seq, seq, length, offset));
	} else {
	    /* Clear */
	    vTcl_SetResult(interp, "%d",
			   set_reference_seq(ed->xx, seq, 0, 0, 0));
	}

    } else if ('s' == *argv[1] && strcmp(argv[1], "set_reference_trace") == 0){
	int seq, control;
	int flags;

	if (argc != 4) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_reference_trace seq -1/0/1",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &seq);
	Tcl_GetInt(interp, argv[3], &control);

	flags = 0;
	if (control < 0)
	    flags = DB_FLAG_REFTRACE_NEG;
	else if (control > 0)
	    flags = DB_FLAG_REFTRACE_POS;

	set_reference_trace(ed->xx, seq, flags);

    } else if ('r' == *argv[1] && strcmp(argv[1], "report_mutations") == 0) {
	int mode = 1; /* use tags */
	int sort_by = 1;
	char *err = NULL;
	dstring_t *html;
	int detail = 2; /* 0 == no html, 1 == basic html, 2 == full html */
	char *dir = NULL;
	int minconf = 15;

	if (argc >= 3)
	    Tcl_GetInt(interp, argv[2], &mode);
	if (argc >= 4)
	    Tcl_GetInt(interp, argv[3], &sort_by);
	if (argc >= 5)
	    dir = argv[4];
	if (argc >= 6)
	    Tcl_GetInt(interp, argv[5], &detail);
	if (argc >= 7)
	    Tcl_GetInt(interp, argv[6], &minconf);

	Tcl_ResetResult(interp);
	html = report_mutations(ed->xx, mode, sort_by, dir, detail,
				minconf, &err);
	if (err) {
	    Tcl_SetResult(interp, err, TCL_VOLATILE);
	    goto fail;
	}
	
	if (html) {
	    Tcl_SetResult(interp, dstring_str(html), TCL_VOLATILE);
	    dstring_destroy(html);
	}

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_flags") == 0) {
	int seq;
	int flags;
	char flag_str[1024];

	if (argc != 3) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_flags seq",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &seq);
	flags = DB_Flags(ed->xx, seq);

	*flag_str = 0;
	if (flags & DB_FLAG_REFTRACE_NEG) {
	    strcat(flag_str, "REFTRACE_NEG ");
	}
	if (flags & DB_FLAG_REFTRACE_POS) {
	    strcat(flag_str, "REFTRACE_POS ");
	}
	if (flags & DB_FLAG_REFSEQ) {
	    strcat(flag_str, "REFSEQ ");
	}
	Tcl_ResetResult(interp);
	Tcl_SetResult(interp, flag_str, TCL_VOLATILE);

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_template_seqs") == 0) {
	int seq;
	int tnum;
	dstring_t *ds = NULL;
	template_c *tarr;
	item_t *item;

	if (argc != 3) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_template_seqs seq",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_ResetResult(interp);

	Tcl_GetInt(interp, argv[2], &seq);
	tnum = DBI(ed->xx)->DB[seq].template;
	tarr = DBI(ed->xx)->templates[tnum];

	if (!tarr)
	    goto fail;
	
	ds = dstring_create(NULL);
	for (item = head(tarr->gel_cont); item; item = item->next) {
	    gel_cont_t *gc = (gel_cont_t *)item->data;
	    char *rname = get_read_name(DBI_io(ed->xx), gc->read);
	    dstring_appendf(ds, "%s %d %d ", rname, gc->contig,
			    io_relpos(DBI_io(ed->xx), gc->read));
	}

	Tcl_SetResult(interp, dstring_str(ds), TCL_VOLATILE);
	dstring_destroy(ds);

    } else if ('s' == *argv[1] && strcmp(argv[1], "show_mini_traces") == 0) {
	int height;

	if (argc != 3) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " show_mini_traces height",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &height);
	edSetMiniTraces(ed->xx, height);

    } else if ('n' == *argv[1] && strcmp(argv[1], "next_difference") == 0) {

	edNextDifference(ed->xx, 1);

    } else if ('p' == *argv[1] && strcmp(argv[1], "prev_difference") == 0) {

	edNextDifference(ed->xx, 0);

    } else if ('v' == *argv[1] && strcmp(argv[1], "view_set") == 0) {
	int set;

	if (argc != 3) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " view_set set_number\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &set);
	edViewSet(ed->xx, set);

    } else if ('m' == *argv[1] && strcmp(argv[1], "move_to_set") == 0) {
	int set;
	char **s_argv;
	int s_argc;

	if (argc < 3) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " move_to_set set_number ?seq ...?\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &set);
	if (Tcl_SplitList(interp, argv[3], &s_argc, &s_argv) != TCL_OK)
	    goto fail;

	edMoveSet(ed->xx, set, s_argc, s_argv);

	Tcl_Free((char *)s_argv);

    } else if ('c' == *argv[1] && strcmp(argv[1], "collapse_set") == 0) {
	int set;

	if (argc != 4 && argc != 3) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " collapse_set set_number ?1/0/-1?\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &set);
	if (argc == 3) {
	    vTcl_SetResult(interp, "%d",
			   ed->xx->set_collapsed
			   ? ed->xx->set_collapsed[set]
			   : 0);
	} else {
	    int mode;
	    Tcl_GetInt(interp, argv[3], &mode);
	    vTcl_SetResult(interp, "%d", edCollapseSet(ed->xx, set, mode));
	}

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_set") == 0) {
	int seq;

	if (argc != 3) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_set seq_number\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &seq);
	vTcl_SetResult(interp, "%d", edFindSet(ed->xx, seq));

    } else if ('s' == *argv[1] && strcmp(argv[1], "show_problem_traces") == 0) {
	int pos;

	if (argc != 3) {
	    Tcl_ResetResult(interp);
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " show_problem_traces position\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &pos);
	tman_problem_traces(ed->xx, pos);
	Tcl_SetResult(interp, "0", TCL_STATIC);

    } else {
	Tcl_AppendResult(interp, "bad option \"", argv[1], "\": must be FIXME",
			 NULL);
	goto fail;
    }

    Tcl_Release((ClientData)TKSHEET(ed));
    return result;

 fail:
    Tcl_Release((ClientData)TKSHEET(ed));
    return TCL_ERROR;
}


/*
 * Callbacks from the sheet widget for our registered extension function.
 * See tkSheet.h for job numbers available and tkSheet.c for where they're
 * called.
 */
static void EditorSheetExtension(ClientData clientData, int job, void *data) {
    Editor *ed = (Editor *)clientData;

    if (!ed->xx)
	return;

    switch (job) {
    case SHEET_JOB_RESIZE:
	if (ed->sw.columns > MAX_DISPLAY_WIDTH) {
	    int font_width = Tk_TextWidth(ed->sw.font, "0", 1);

	    ed->sw.columns = MAX_DISPLAY_WIDTH;
	    ed->sw.width_in_pixels = ed->sw.columns * font_width +
		2 * (ed->sw.border_width /* + ed->sw.pad_x */);
	}

	ed->xx->displayWidth = ed->sw.columns;
	ed->xx->refresh_flags |= ED_DISP_ALL | ED_DISP_HEIGHT;
	redisplaySequences(ed->xx, 0);
	break;

    case SHEET_JOB_DESTROY:
	if (ed->xx) {
	    delete_edStruct(ed->xx);
	    ed->xx = NULL;
	}

	break;
    }
}
