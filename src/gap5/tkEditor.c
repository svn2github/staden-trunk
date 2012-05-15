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
#include <string.h>

#include "tkEditor.h"
#include "tkSheet_common.h"
#include "tk_defs.h"
//#include "undo.h"
#include "tman_interface.h"
#include "editor_oligo.h"
#include "xalloc.h"
#include "tcl_utils.h"
//#include "edCommands.h"
//#include "contigEditor.h"
#include "gap_globals.h"
#include "dstring.h"
#include "gap4_compat.h"
#include "tagdb.h"
#include "notedb.h"

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
	 "-qualcolour0","selQualColour0",	"Background",	"#2929ff",
	 offset(qual_bg2[0]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour1","selQualColour1",	"Background",	"#3939ff",
	 offset(qual_bg2[1]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour2","selQualColour2",	"Background",	"#4949ff",
	 offset(qual_bg2[2]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour3","selQualColour3",	"Background",	"#5959ff",
	 offset(qual_bg2[3]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour4","selQualColour4",	"Background",	"#6969ff",
	 offset(qual_bg2[4]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour5","selQualColour5",	"Background",	"#7979ff",
	 offset(qual_bg2[5]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour6","selQualColour6",	"Background",	"#8989ff",
	 offset(qual_bg2[6]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour7","selQualColour7",	"Background",	"#9999ff",
	 offset(qual_bg2[7]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour8","selQualColour8",	"Background",	"#a9a9ff",
	 offset(qual_bg2[8]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qualcolour9","selQualColour9",	"Background",	"#b9b9ff",
	 offset(qual_bg2[9]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-qual_fg","qualForeground",	"Foreground","#ff5050",
	 offset(qual_below),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-diff1_bg","diff1Background",	"Background","green3",
	 offset(diff1_bg),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-diff2_bg","diff2Background",	"Background","green3",
	 offset(diff2_bg),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-diff1_fg","diff1Foreground",	"Foreground","green3",
	 offset(diff1_fg),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-diff2_fg","diff2Foreground",	"Foreground","green3",
	 offset(diff2_fg),		0, NULL},
    {TK_CONFIG_COLOR,
         "-stripe_bg","stripeBackground", "Background","#f0f0f0",
	 offset(stripe_bg),		0, NULL},
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
	 "-tmplcolour0","tmplColour0",	"Background",	"#800000",
	 offset(tmpl_bg[0]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour1","tmplColour1",	"Background",	"#c0c0ff",
	 offset(tmpl_bg[1]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour2","tmplColour2",	"Background",	"white",
	 offset(tmpl_bg[2]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour3","tmplColour3",	"Background",	"#888888",
	 offset(tmpl_bg[3]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour4","tmplColour4",	"Background",	"#ff2020",
	 offset(tmpl_bg[4]),		0, NULL},
    {TK_CONFIG_COLOR,
	 "-tmplcolour5","tmplColour5",	"Background",	"orange",
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
    {TK_CONFIG_BOOLEAN,
         "-display_cutoffs", "displayCutoffs", "DisplayCutoffs",
         "0", offset(display_cutoffs), 0, NULL},
    {TK_CONFIG_BOOLEAN,
         "-display_quality", "displayQuality", "DisplayQuality",
         "0", offset(display_quality), 0, NULL},
    {TK_CONFIG_BOOLEAN,
         "-display_mapping_quality", "displayMappingQuality",
         "DisplayMappingQuality",
         "0", offset(display_mapping_quality), 0, NULL},
    {TK_CONFIG_INT,
         "-display_differences", "displayDifferences", "DisplayDifferences",
         "0", offset(display_differences), 0, NULL},
    {TK_CONFIG_INT,
         "-display_differences_quality", "displayDifferencesQuality",
         "DisplayDifferencesQuality",
         "0", offset(display_differences_qual), 0, NULL},
    {TK_CONFIG_INT,
         "-consensus_at_top", "consensusAtTop", "ConsensusAtTop",
         "0", offset(consensus_at_top), 0, NULL},
    {TK_CONFIG_BOOLEAN,
         "-differences_case_sensitive", "differencesCaseSensitive",
         "DifferencesCaseSensitive",
         "0", offset(display_differences_case), 0, NULL},
    {TK_CONFIG_INT,
         "-stripe_mode", "stripeMode", "StripeMode",
         "0", offset(stripe_mode), 0, NULL},
    {TK_CONFIG_INT,
     "-stack_mode", "stackMode", "StackMode",
         "0", offset(stack_mode), 0, NULL},
    {TK_CONFIG_INT,
         "-hide_anno", "hideAnno", "HideAnno",
         "0", offset(hide_annos), 0, NULL},
    {TK_CONFIG_INT,
         "-pos_type", "posType", "PosType",
         "0", offset(pos_type), 0, NULL},
    {TK_CONFIG_INT,
         "-group_by_primary", "groupByPrimary", "GroupByPrimary",
         "0", offset(group_primary), 0, NULL},
    {TK_CONFIG_INT,
         "-group_by_secondary", "groupBySecondary", "GroupBySecondary",
         "0", offset(group_secondary), 0, NULL},
    {TK_CONFIG_STRING,
         "-output_list", "outputList", "OutputList",
         "readings", offset(output_list), TK_CONFIG_NULL_OK, NULL},
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

static int attach_edview(Tcl_Interp *interp, Editor *ed,
			 int argc, char **argv) {
    edNames *name_sheet;
    Tcl_CmdInfo cmdinfo;
    edview *xx;
    GapIO *io;
    Tcl_Obj *io_str;
    Tk_Window tkwin = ed->sw.tkwin;

    if (argc != 7) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], " init io contig_num cursor_rec cursor_pos name_path\"",
			 (char *) NULL);
	return TCL_ERROR;
    }


    /* Convert the ed->io into a real io object */
    io_str = Tcl_NewStringObj(argv[2], -1);
    io = io_from_obj(io_str);
    Tcl_DecrRefCount(io_str);
    if (!io)
	return TCL_ERROR;

    /* Obtain seq and name sheet widget from the window pathnames */
    if (0 == Tcl_GetCommandInfo(interp, argv[6], &cmdinfo)) {
	return TCL_ERROR;
    }
    name_sheet = (edNames *)cmdinfo.clientData;

    xx = edview_new(io, atorec(argv[3]), atorec(argv[4]), atoi(argv[5]),
		    ed, name_sheet, NULL, interp);
    if (!xx) {
	return TCL_ERROR;
    }

    strncpy(xx->seq_win,  Tk_PathName(tkwin), WIN_NAME_SIZE);
    strncpy(xx->name_win, argv[4], WIN_NAME_SIZE);
    
    ed->xx = xx;
    name_sheet->xx = xx;
    return TCL_OK;
}

static Pixel ColourNameToPixel(Tcl_Interp *interp, Tk_Window tkwin, char *name)
{
    XColor *cp = Tk_GetColor(interp, tkwin, name);

    if (cp) {
	return cp->pixel;
    } else {
	verror(ERR_WARN, "ColourNameToPixel", "Colourmap is full");
	return 0;
    }
}

static void setUpColourMap(Tcl_Interp *interp, Tk_Window tkwin)
{
    static int done = 0;
    int i;

    if (done)
	return;
    else
	done = 1;

    for (i=0;i<tag_db_count;i++) {
        tag_db[i].fg_pixel =  (tag_db[i].fg_colour == NULL) ?
            1 : ColourNameToPixel(interp, tkwin, tag_db[i].fg_colour);
        tag_db[i].bg_pixel =  (tag_db[i].bg_colour == NULL) ?
            0 : ColourNameToPixel(interp, tkwin, tag_db[i].bg_colour);
        tag_db[i].gf_pixel =  (tag_db[i].gf_colour == NULL) ?
            1 : ColourNameToPixel(interp, tkwin, tag_db[i].gf_colour);
        tag_db[i].gb_pixel =  (tag_db[i].gb_colour == NULL) ?
            0 : ColourNameToPixel(interp, tkwin, tag_db[i].gb_colour);
    }

    for (i = 0; i < note_db_count; i++) {
	note_db[i].fg_pixel = note_db[i].fg_colour
	    ? ColourNameToPixel(interp, tkwin, note_db[i].fg_colour) : 1;
	note_db[i].bg_pixel = note_db[i].bg_colour
	    ? ColourNameToPixel(interp, tkwin, note_db[i].bg_colour) : 0;
	note_db[i].gf_pixel = note_db[i].gf_colour
	    ? ColourNameToPixel(interp, tkwin, note_db[i].gf_colour) : 1;
	note_db[i].gb_pixel = note_db[i].gb_colour
	    ? ColourNameToPixel(interp, tkwin, note_db[i].gb_colour) : 0;
    }
}

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
     * Initialise tag colours.
     */
    setUpColourMap(interp, Tk_MainWindow(interp));


    /*
     * Initialised rest of Editor structure.
     */
    TKSHEET(ed)->extensionData = (ClientData)ed;
    TKSHEET(ed)->extension = EditorSheetExtension;
    for (i = 0; i < 10; i++) {
	ed->qual_bg[i] = NULL;
	ed->qual_bg2[i] = NULL;
    }
    for (i = 0; i < 4; i++)
	ed->edit_bg[i] = NULL;
    for (i = 0; i < 6; i++)
	ed->tmpl_bg[i] = NULL;
    for (i = 0; i < 10; i++)
	ed->set_bg[i] = NULL;
    ed->qual_below = NULL;
    ed->diff1_bg    = NULL;
    ed->diff2_bg    = NULL;
    ed->diff1_fg    = NULL;
    ed->diff2_fg    = NULL;
    ed->stripe_bg   = NULL;
    ed->xScrollCmd = NULL;
    ed->yScrollCmd = NULL;
    ed->highlight_cmd = NULL;
    ed->output_list = NULL;
    ed->max_height = 0;
    ed->xx = NULL;
    /*
     * Ideally we want gridding enabled so that "half characters" are not
     * visible, but with causes major headaches with metacity (the default
     * gnome2 window manager).
     */
    ed->grid = 1;


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


static int EditorWidgetCmd(ClientData clientData, Tcl_Interp *interp,	
			   int argc, char **argv) {
    Editor *ed = (Editor *)clientData;
    int result = TCL_OK;
    int index;
    Tcl_Obj *s_obj;
    edview *xx = ed->xx;

    static char *optionStrings[] = {
	"configure",     "init",         "io",	         "redraw",
	"xview",         "yview",        "get_number",   "get_seq_status",
	"set_cursor",    "contig_start", "contig_end",   "contig_rec",
	"display_trace", "delete_trace", "trace_scroll", "save",
	"cursor_up",     "cursor_down",  "cursor_left",  "cursor_right",
	"read_start",    "read_start2",  "read_end",     "read_end2",
	"get_template_seqs", "edits_made", "link_to",    "lock",
	"join_align",	 "join_mismatch","join",         "join_offset",
	"select",	 "edit_annotation",  "clear_visibility_cache",
	"cursor_id",     "get_cursor",	  "search",      "get_xy",
	"change_contig", "select_oligo","show_cursor",
	"reference_pos", "next_difference", "prev_difference",
	"set_sort_order", "set_trace_lock", "set_base_sort_point", 
	"set_sequence_sort", NULL
    };
    enum options {
	_CONFIGURE,      _INIT,          _IO,            _REDRAW,
	_XVIEW,          _YVIEW,         _GET_NUMBER,    _GET_SEQ_STATUS,
	_SET_CURSOR,     _CONTIG_START,  _CONTIG_END,    _CONTIG_REC,
	_DISPLAY_TRACE,  _DELETE_TRACE,  _TRACE_SCROLL,  _SAVE,
	_CURSOR_UP,      _CURSOR_DOWN,   _CURSOR_LEFT,   _CURSOR_RIGHT,
	_READ_START,     _READ_START2,   _READ_END,      _READ_END2,
	_GET_TEMPLATE_SEQS, _EDITS_MADE, _LINK_TO,       _LOCK,
	_JOIN_ALIGN,     _JOIN_MISMATCH, _JOIN,          _JOIN_OFFSET,
	_SELECT,	 _EDIT_ANNOTATION,  _CLEAR_VISIBILITY_CACHE,
	_CURSOR_ID,      _GET_CURSOR,	 _SEARCH,	 _GET_XY,
	_CHANGE_CONTIG,  _SELECT_OLIGO,  _SHOW_CURSOR,
	_REFERENCE_POS,  _NEXT_DIFFERENCE, _PREV_DIFFERENCE,
	_SET_SORT_ORDER, _SET_TRACE_LOCK, _SET_BASE_SORT_POINT,
	_SET_SEQUENCE_SORT
    };

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " option ?arg arg ...?\"",
                         (char *) NULL);
	goto fail;
    }

    Tcl_Preserve((ClientData)TKSHEET(ed));

    s_obj = Tcl_NewStringObj(argv[1], -1);
    Tcl_IncrRefCount(s_obj);
    if (Tcl_GetIndexFromObj(interp, s_obj, optionStrings, "option", 0,
	    &index) != TCL_OK) {
	Tcl_DecrRefCount(s_obj);
	goto fail;
    }
    Tcl_DecrRefCount(s_obj);

    switch ((enum options)index) {
    default:
	goto fail;

    /* Reconfigure widget */
    case _CONFIGURE:
	result = SheetWidgetCmdConfig(interp, TKSHEET(ed), argc, argv);
	break;

    /* Initialise by creating an edview struct */
    case _INIT:
	result = attach_edview(interp, ed, argc, argv);
	break;

    /* GapIO */
    case _IO:
	if (argc == 3) {
	    /* Update the io we're using, eg due to a contig join */
	    Tcl_Obj *io_str;
	    GapIO *io;

	    io_str = Tcl_NewStringObj(argv[2], -1);
	    io = io_from_obj(io_str);
	    Tcl_DecrRefCount(io_str);
	    if (!io)
		return TCL_ERROR;
	    xx->io = io;
	}
	Tcl_SetResult(interp, io_obj_as_string(xx->io) , TCL_VOLATILE);
	break;
	
    /* Contig record number */
    case _CONTIG_REC:
	Tcl_SetObjResult(interp, Tcl_NewWideIntObj(xx->cnum));
	break;

    /* Save contig */
    case _SAVE: {
	reg_length rl;
	contig_t *c;

	/*
	 * We shouldn't have to do the decr, flush, search & incr here.
	 * However cache_flush on a child I/O is currently buggy in that
	 * it moves the pointer to the parent object (OK) and in doing
	 * so fails to correct the reference counting hash (not OK).
	 *
	 * This is a convenient work around until we update the caching
	 * layer.
	 */
	result = cache_flush(xx->io);
	c = cache_search(xx->io, GT_Contig, xx->cnum);

	rl.job = REG_LENGTH;
	rl.length = c->end - c->start + 1;
	contig_notify(xx->io->base, c->rec, (reg_data *)&rl);

	vTcl_SetResult(interp, "%d", result);
	break;
    }

    /* Query: have we been editing? */
    case _EDITS_MADE:
	result = cache_updated(xx->io);
	vTcl_SetResult(interp, "%d", result);
	if (result != -1)
	    result = 0;
	break;

    /* Refresh display after external changes */
    case _REDRAW: {
	int reload_seq = 0;
	
	if (argc == 3) {
	    reload_seq = atoi(argv[2]);
	}

	ed->xx->refresh_flags = ED_DISP_ALL;
	if (reload_seq && ed->xx->r) {
	    free(ed->xx->r);
	    ed->xx->r = NULL;
	    ed->xx->nr = 0;
	}

	edview_redraw(ed->xx);
	break;
    }
	
    /* X Scroll */
    case _XVIEW: {
	double f1;
	int type, count, offset;
	contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);

	if (argc == 2) {
	    /* editor xview */
	    vTcl_SetResult(interp, "%d", xx->displayPos);
	    break;

	} else if (argc == 3) {
	    /* editor xview <position> */
	    offset = atoi(argv[2]);

	} else {
	    /* else 'moveto ?' or 'scroll ? units' or 'scroll ? pages' */

	    type = Tk_GetScrollInfo(interp, argc, argv, &f1, &count);
	    switch (type) {
	    default:
		goto fail;
		
	    case TK_SCROLL_MOVETO:
		offset = f1 * contig_get_length(&c) + contig_get_start(&c);
		break;

	    case TK_SCROLL_PAGES:
		offset = xx->displayPos + count * 0.9 * xx->displayWidth;
		break;

	    case TK_SCROLL_UNITS:
		offset = xx->displayPos + count;
		break;
	    }
	}

	/* Bounds checking */
	if (offset < c->start - xx->displayWidth + 1)
	    offset = c->start - xx->displayWidth + 1;
	if (offset > c->end)
	    offset = c->end;

	set_displayPos(xx, offset);
	break;
    }

    /* Y Scroll */
    case _YVIEW: {
	double f1;
	int type, count, offset;

	/* Cache elsewhere - ie last number of seqs visible? */
	edview_visible_items(xx, xx->displayPos,
			     xx->displayPos + xx->displayWidth);

	//type = Tk_GetScrollInfoObj(interp, objc, objv, &f1, &count);
	type = Tk_GetScrollInfo(interp, argc, argv, &f1, &count);
	if (!xx->ed->consensus_at_top)
	    count = -count;

	if (f1 < 0)
	    f1 = 0;

	switch (type) {
	default:
	    goto fail;

	case TK_SCROLL_MOVETO:
	    if (xx->ed->consensus_at_top)
		offset = f1 * xx->max_height + 0.5;
	    else
		offset = (1-f1) * xx->max_height - (xx->displayHeight-1) + 0.5;
	    break;

	case TK_SCROLL_PAGES:
	    offset = xx->displayYPos + count * 0.9 * xx->displayHeight;
	    break;

	case TK_SCROLL_UNITS:
	    offset = xx->displayYPos + count;
	    break;
	}

	if (offset >= 0)
	    xx->displayYPos = offset;
	xx->refresh_flags = ED_DISP_YSCROLL;
	edview_redraw(xx);
	break;
    }

    /* Get reading number under at x,y coord (eg mouse pointer), or
     * the current editior cursor if no x,y specified.
     */
    case _GET_NUMBER: {
	char buf[100];
	int x, y, type, pos;
	tg_rec rec;

	if (argc != 2 && !(argc >= 4 && argc <= 6)) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_number ?xpos ypos ?exact ?seq_only???\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc >= 4) {
	    int exact = argc == 5 ? atoi(argv[4]) : 1;
	    int seq_only = argc == 6 ? atoi(argv[5]) : 0;
	    sheet_arg_x(TKSHEET(ed), argv[2], &x); /* cell coordinates */
	    sheet_arg_y(TKSHEET(ed), argv[3], &y); y++;

	    if (-1 != (type = edview_item_at_pos(ed->xx, y, x, 0, exact,
						 seq_only, &rec, &pos))) {
		sprintf(buf, "%d %"PRIrec" %d", type, rec, pos);
		Tcl_AppendResult(interp, buf, NULL);
	    } /* otherwise return a blank */
	} else {
	    sprintf(buf, "%d %"PRIrec" %d",
		    ed->xx->cursor_type,
		    ed->xx->cursor_rec,
		    ed->xx->cursor_pos);
	    Tcl_AppendResult(interp, buf, NULL);
	}
	break;
    }

    /* Given a record number and offset, turn this into an x,y coord in
     * "sheet" units. This is useful if we want to go from editor cursor
     * record/pos to x/y and back to top-most record/pos (eg tag on top
     * of current sequence) via get_number.
     */
    case _GET_XY: {
	int type, pos;
	tg_rec rec;
	int x,y;

	if (argc != 2 && argc != 5) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_xy ?rec_type rec_no position?"
			     " ?make_visible?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 5) {
	    type = atoi(argv[2]);
	    rec  = atorec(argv[3]);
	    pos  = atoi(argv[4]);
	} else {
	    type = ed->xx->cursor_type;
	    rec  = ed->xx->cursor_rec;
	    pos  = ed->xx->cursor_pos;
	}

	if (edGetXY(xx, type, rec, pos, &x, &y) == 0) {
	    vTcl_SetResult(interp, "%d %d", x, y);
	} else {
	    Tcl_ResetResult(interp);
	}

	break;
    }

    /* Sets the editor cursor position */
    case _SET_CURSOR: {
	int visible = 1;
	if (argc != 5 && argc != 6) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " set_cursor rec_type rec_no position"
			     " ?make_visible?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc >= 6)
	    Tcl_GetInt(interp, argv[5], &visible);
	edSetCursorPos(ed->xx, atoi(argv[2]), atorec(argv[3]), atoi(argv[4]),
		       visible);
	break;
    }

    /* Gets the cursor position, either absolute or relative */
    case _GET_CURSOR: {
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_cursor {absolute|relative}\"",
			     (char *) NULL);
	    goto fail;
	}

	if (*argv[2] == 'a')
	    vTcl_SetResult(interp, "%d %"PRIrec" %d",
			   GT_Contig, xx->cnum, xx->cursor_apos);
	else
	    vTcl_SetResult(interp, "%d %"PRIrec" %d",
			   xx->cursor_type, xx->cursor_rec, xx->cursor_pos);
	
	break;
    }

    /* Forces the cursor to be visible, returns 1 if it wasn't before. */
    case _SHOW_CURSOR:
	vTcl_SetResult(interp, "%d", showCursor(xx, 0, 0));
	break;

    /* Cursor movement operations */
    case _CURSOR_UP:
	if (xx->ed->consensus_at_top)
	    edCursorUp(ed->xx);
	else
	    edCursorDown(ed->xx);
	break;

    case _CURSOR_DOWN:
	if (xx->ed->consensus_at_top)
	    edCursorDown(ed->xx);
	else
	    edCursorUp(ed->xx);
	break;

    case _CURSOR_LEFT:
	edCursorLeft(ed->xx);
	break;

    case _CURSOR_RIGHT:
	edCursorRight(ed->xx);
	break;

    case _READ_START:
	edReadStart(ed->xx);
	break;

    case _READ_START2:
	edReadStart2(ed->xx);
	break;

    case _READ_END:
	edReadEnd(ed->xx);
	break;

    case _READ_END2:
	edReadEnd2(ed->xx);
	break;

    case _CONTIG_START:
	edContigStart(ed->xx);
	break;

    case _CONTIG_END:
	edContigEnd(ed->xx);
	break;

    case _DISPLAY_TRACE:
	edDisplayTrace(ed->xx);
	break;

    case _DELETE_TRACE:
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " delete_trace path\"",
			     (char *) NULL);
	    goto fail;
	}

	deleteTrace(ed->xx, argv[2]);
	break;

    case _TRACE_SCROLL:
	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " trace_scroll path scroll_cmd\"",
			     (char *) NULL);
	    goto fail;
	}

	edScrollTraces(ed->xx, argv[2], argv[3]);
	break;

    case _GET_SEQ_STATUS: {
	char *msg;

	if (argc != 6) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_seq_status rec_type rec_no"
			              " position format\"",
			     (char *) NULL);
	    goto fail;
	}

	switch (atoi(argv[2])) {
	case GT_Seq:
	    msg = edGetBriefSeq(ed->xx, atorec(argv[3]), atoi(argv[4]), argv[5]);
	    break;

	case GT_Contig:
	    msg = edGetBriefCon(ed->xx, atorec(argv[3]), atoi(argv[4]), argv[5]);
	    break;

	case GT_AnnoEle:
	    msg = edGetBriefTag(ed->xx, atorec(argv[3]), argv[5]);
	    break;

	default:
	    msg = "";
	    break;
	}
	Tcl_SetResult(interp, msg, TCL_VOLATILE);
	break;
    }

    case _GET_TEMPLATE_SEQS: {
	int nrec, i;
	tg_rec *rec;
	dstring_t *ds = dstring_create(NULL);

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_template_seqs rec.no.\"",
			     (char *) NULL);
	    goto fail;
	}

	rec = edGetTemplateReads(ed->xx, atorec(argv[2]), &nrec);	
	for (i = 0; i < nrec; i++) {
	    dstring_appendf(ds, "%"PRIrec" ", rec[i]);
	}
	Tcl_SetResult(interp, dstring_str(ds), TCL_VOLATILE);
	dstring_destroy(ds);
	break;
    }

    case _LINK_TO: {
	Editor *ed1 = ed, *ed2;
	tkSheet *diffs;
	Tcl_CmdInfo cmdinfo;
	
	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " link_to ed2 diffs.\"",
			     (char *) NULL);
	    goto fail;
	}

	if (0 == Tcl_GetCommandInfo(interp, argv[2], &cmdinfo)) {
	    return TCL_ERROR;
	}
	ed2 = (Editor *)cmdinfo.clientData;
	
	if (0 == Tcl_GetCommandInfo(interp, argv[3], &cmdinfo)) {
	    return TCL_ERROR;
	}
	diffs = (tkSheet *)cmdinfo.clientData;

	if (NULL == (ed1->xx->link = (EdLink *)malloc(sizeof(EdLink))))
	    return TCL_ERROR;
	ed2->xx->link = ed1->xx->link;
	ed1->xx->link->xx[0] = ed1->xx;
	ed1->xx->link->xx[1] = ed2->xx;
	ed1->xx->link->diffs = diffs;
	ed1->xx->link->locked = 1;
	ed1->xx->link->lockOffset = ed2->xx->displayPos - ed1->xx->displayPos;
	
	break;
    }

    case _LOCK: {
	int val;
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " lock value\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &val);
	if (ed->xx->link) {
	    ed->xx->link->locked = val;
	    ed->xx->link->lockOffset = ed->xx->link->xx[1]->displayPos
		- ed->xx->link->xx[0]->displayPos;
	}
	break;
    }

    case _JOIN_ALIGN: {
	int fixed_left = 0;
	int fixed_right = 0;
	if (argc != 4 && argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " join_align ?fixed_left fixed_right?\"",
			     (char *) NULL);
	    goto fail;
	}
	if (argc == 4) {
	    Tcl_GetInt(interp, argv[2], &fixed_left);
	    Tcl_GetInt(interp, argv[3], &fixed_right);
	}

	result = (0 == edJoinAlign(ed->xx, fixed_left, fixed_right))
	    ? TCL_OK
	    : TCL_ERROR;
	break;
    }

    case _JOIN_MISMATCH: {
	int len, mismatch;

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " join_mismatch\"",
			     (char *) NULL);
	    goto fail;
	}

	result = (0 == edJoinMismatch(ed->xx, &len, &mismatch))
	    ? TCL_OK
	    : TCL_ERROR;
	vTcl_SetResult(interp, "%d %d", len, mismatch);
	break;
    }

    case _JOIN:
	result = edJoin(ed->xx) == 0 ? TCL_OK : TCL_ERROR;
	break;

    case _JOIN_OFFSET:
	if (!xx->link) {
	    result = TCL_ERROR;
	    break;
	}

	Tcl_SetObjResult(interp, Tcl_NewIntObj(xx->link->lockOffset));
	break;

    case _SELECT: {
	char *select_options[] = {"clear", "from", "to", "get", "set", NULL};
	enum  select_options     { CLEAR,   FROM,   TO,   GET,   SET};
	Tcl_Obj *obj;

	if (argc < 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " select sub-cmd ?options?\"",
			     (char *) NULL);
	    goto fail;
	}

	obj = Tcl_NewStringObj(argv[2], -1);
	Tcl_IncrRefCount(obj);
	if (Tcl_GetIndexFromObj(interp, obj, select_options, "option", 0,
				&index) != TCL_OK) {
	    Tcl_DecrRefCount(obj);
	    goto fail;
	}
	Tcl_DecrRefCount(obj);

	if (index == CLEAR) {
	    if (argc != 3) {
		Tcl_AppendResult(interp, "wrong # args: should be \"",
				 argv[0], " select clear\"",
				 (char *) NULL);
		goto fail;
	    }
	    edSelectClear(ed->xx);
	} else if (index == TO || index == FROM) {
	    int arg;
	    if (argc != 4) {
		Tcl_AppendResult(interp, "wrong # args: should be \"",
				 argv[0], " select ",
				 index == TO ? "to" : "from",
				 " arg\"",
				 (char *) NULL);
		goto fail;
	    }

	    sheet_arg_x(TKSHEET(ed), argv[3], &arg);
	    if (index == TO) {
		edSelectTo(ed->xx, arg);
	    } else {
		edSelectFrom(ed->xx, arg);
	    }
	    
	} else if (index == GET) {
	    if (ed->xx->select_seq) {
		vTcl_SetResult(interp, "%d %"PRIrec" %d %d",
			       ed->xx->select_seq == ed->xx->cnum
			       ? GT_Contig
			       : GT_Seq,
			       ed->xx->select_seq,
			       ed->xx->select_start,
			       ed->xx->select_end);
	    } else {
		/* The base under the editing cursor */
		vTcl_SetResult(interp, "%d %"PRIrec" %d %d",
			       ed->xx->cursor_rec == ed->xx->cnum
			       ? GT_Contig
			       : GT_Seq,
			       ed->xx->cursor_rec,
			       ed->xx->cursor_pos,
			       ed->xx->cursor_pos);
	    }

	} else if (index == SET) {
	    tg_rec rn = ed->xx->cnum;
	    int bs = 0;
	    int be = 0;

	    if (argc != 5 && argc != 6) {
		Tcl_AppendResult(interp, "wrong # args: should be \"",
				 argv[0], " select set ",
				 " ?rec_num? base_start base_end\"",
				 (char *) NULL);
		goto fail;
	    }

	    if (argc == 5) {
		bs = atoi(argv[3]);
		be = atoi(argv[4]);
	    } else {
		rn = atorec(argv[4]);
		bs = atoi(argv[5]);
		be = atoi(argv[6]);
	    }

	    edSelectSet(ed->xx, rn, bs, be);

	} else {
	    Tcl_AppendResult(interp, "wrong sub-command: should be "
			     "clear, from, to or get",
			     (char *) NULL);
	    goto fail;
	}
	break;
    }

    case _EDIT_ANNOTATION:
	//tagEditor(ed->xx);
	break;

    case _CURSOR_ID:
	Tcl_SetObjResult(interp, Tcl_NewIntObj(xx && xx->cursor
					       ? xx->cursor->id
					       : -1));
	break;

    case _SEARCH: {
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
	
	found = edview_search(ed->xx, (strcmp("backward", argv[2]) != 0),
			      *argv[3], argv[4], value);
	vTcl_SetResult(interp, "%d", found);

	break;
    }

    case _CHANGE_CONTIG:
	if (argc == 3)
	    ed->xx->cnum = atorec(argv[2]);
	break;

    case _SELECT_OLIGO: {
	Tcl_Obj *obj;

	if (argc != 7) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " select_oligo is_fwds dist_fwd "
			     "dist_back avg_readlen primer3_params\"",
			     (char *) NULL);
	    goto fail;
	}

	obj = edSelectOligoGenerate(ed->xx,
				    atoi(argv[2]) /* is_fwds? */,
				    atoi(argv[3]) /* fwd */,
				    atoi(argv[4]) /* bwd */,
				    atoi(argv[5]) /* readlen */,
				    argv[6]  /* p3_def */);

	if (obj)
	    Tcl_SetObjResult(interp, obj);
	else
	    Tcl_ResetResult(interp);
		       
	break;
    }

    case _CLEAR_VISIBILITY_CACHE:
	if (xx->r)
	    free(xx->r);
	xx->r = NULL;

	break;

    case _REFERENCE_POS: {
	int rpos, dir, rid;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " reference_pos padded_coord\"",
			     (char *) NULL);
	    goto fail;
	}
	
	rpos = padded_to_reference_pos(ed->xx->io, ed->xx->cnum,
				       atoi(argv[2]), &dir, &rid);

	vTcl_SetResult(interp, "%d %d %d", rpos, dir, rid);
	break;
    }

    case _NEXT_DIFFERENCE:
	edNextDifference(ed->xx);
	break;

    case _PREV_DIFFERENCE:
	edPrevDifference(ed->xx);
	break;
	
    case _SET_SORT_ORDER:
	edview_set_sort_order(ed->xx);
	break;
	
    case _SET_BASE_SORT_POINT:
    	ed_set_base_sort_point(ed->xx);
	break;

    case _SET_TRACE_LOCK: {
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
    }
    
    case _SET_SEQUENCE_SORT:
    	ed_set_sequence_sort(ed->xx);
    	break;
    
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
	ed->xx->displayHeight = ed->sw.rows;
	ed->xx->refresh_flags |= ED_DISP_ALL | ED_DISP_HEIGHT;
	edview_redraw(ed->xx);
	break;

    case SHEET_JOB_DESTROY:
	if (ed->xx) {
	    edview_destroy(ed->xx);
	    ed->xx = NULL;
	}

	break;
    }
}
