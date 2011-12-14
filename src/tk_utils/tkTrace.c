/*
 * The Trace widget.
 *
 * Bugs: -font only works for initial creation of widget.
 */


#include <staden_config.h>
#include <stdlib.h>
#include <stdio.h>          /* For sprintf() */
#include <tk.h>
#include <math.h>
#include <X11/Xutil.h>
#ifdef HAVE_PNG
#    include <png.h>
#endif
#include <os.h>
#include <io_lib/Read.h>
#include <io_lib/traceType.h>

#include "tk_defs.h"

#include "os.h"
#include "tkTrace.h"
#include "tkTraceIO.h"
#include "trace_print.h"
#include "misc.h"
#include "split.h"
#include "tcl_utils.h"

#define DEF_TRACE_HEIGHT "220"
#define DEF_TRACE_WIDTH  "700"
/* #define DEF_TRACE_FONT   "Fixed -12 roman" */
#define DEF_TRACE_FONT   "trace_font"
#define DEF_TRACE_CONF_FONT   "trace_conf_font"

/* Stipple definitions */
#define stripe4_width 16
#define stripe4_height 16
static unsigned char stripe4_bits[] = {
   0xc3, 0xc3, 0x87, 0x87, 0x0f, 0x0f, 0x1e, 0x1e, 0x3c, 0x3c, 0x78, 0x78,
   0xf0, 0xf0, 0xe1, 0xe1, 0xc3, 0xc3, 0x87, 0x87, 0x0f, 0x0f, 0x1e, 0x1e,
   0x3c, 0x3c, 0x78, 0x78, 0xf0, 0xf0, 0xe1, 0xe1};

static Tk_ConfigSpec configSpecs[] = {
    {TK_CONFIG_BORDER, "-background", "background", "Background", NORMAL_BG,
    Tk_Offset(DNATrace, border), TK_CONFIG_COLOR_ONLY,
	 (Tk_CustomOption *) NULL},
    {TK_CONFIG_BORDER, "-background", "background", "Background", "white",
    Tk_Offset(DNATrace, border), TK_CONFIG_MONO_ONLY,
	 (Tk_CustomOption *) NULL},
    {TK_CONFIG_SYNONYM, "-bd", "borderWidth", (char *) NULL, (char *) NULL,
	 0, 0, (Tk_CustomOption *) NULL},
    {TK_CONFIG_SYNONYM, "-bg", "background", (char *) NULL,  (char *) NULL,
	 0, 0, (Tk_CustomOption *) NULL},
    {TK_CONFIG_PIXELS, "-borderwidth", "borderWidth", "BorderWidth", "0.5m",
    Tk_Offset(DNATrace, borderWidth), 0, (Tk_CustomOption *) NULL},
    {TK_CONFIG_PIXELS, "-height", "height", "Height", DEF_TRACE_HEIGHT,
    Tk_Offset(DNATrace, height), 0, NULL},
    {TK_CONFIG_RELIEF, "-relief", "relief", "Relief", "raised",
    Tk_Offset(DNATrace, relief), 0, (Tk_CustomOption *) NULL},
    {TK_CONFIG_PIXELS, "-width",  "width",  "Width",  DEF_TRACE_WIDTH,
    Tk_Offset(DNATrace, width), 0, NULL},
    {TK_CONFIG_INT, "-cursor_pos", "cursorPos", "CursorPos", "0",
    Tk_Offset(DNATrace, cursor_pos), 0, NULL},
    {TK_CONFIG_COLOR, "-colour1", "colour1", "Colour1", "#00c500",
    Tk_Offset(DNATrace, Acol), 0, NULL},
    {TK_CONFIG_COLOR, "-colour2", "colour2", "Colour2", "blue",
    Tk_Offset(DNATrace, Ccol), 0, NULL},
    {TK_CONFIG_COLOR, "-colour3", "colour3", "Colour3", "black",
    Tk_Offset(DNATrace, Gcol), 0, NULL},
    {TK_CONFIG_COLOR, "-colour4", "colour4", "Colour4", "orangered",
    Tk_Offset(DNATrace, Tcol), 0, NULL},
    {TK_CONFIG_COLOR, "-cursor_colour", "cursorColour", "CursorColour", "blue",
    Tk_Offset(DNATrace, CursorCol), 0, NULL},
    {TK_CONFIG_COLOR, "-qual_colour", "qualColour", "QualColour", "grey",
    Tk_Offset(DNATrace, CutoffCol), 0, NULL},
    {TK_CONFIG_COLOR, "-vector_colour", "vectorColour", "VectorColour", "pink",
    Tk_Offset(DNATrace, VectorCol), 0, NULL},
    {TK_CONFIG_COLOR, "-conf_colour", "confColour", "ConfColour", "skyblue3",
    Tk_Offset(DNATrace, ConfCol), 0, NULL},
    {TK_CONFIG_COLOR, "-conf_neg_colour", "confColour", "ConfColour","pink1",
    Tk_Offset(DNATrace, ConfNegCol), 0, NULL},
    {TK_CONFIG_STRING, "-xscrollcommand", "xScrollCommand", "ScrollCommand",
    "", Tk_Offset(DNATrace, xScrollCmd), TK_CONFIG_NULL_OK, NULL},
    {TK_CONFIG_BOOLEAN, "-shownumbers", "showNumbers", "ShowNumbers", "1",
    Tk_Offset(DNATrace, show_numbers), 0, NULL},
    {TK_CONFIG_BOOLEAN, "-showsequence", "showSequence", "ShowSequence", "1",
    Tk_Offset(DNATrace, show_sequence ), 0, NULL},
    {TK_CONFIG_BOOLEAN, "-showedits", "showEdits", "ShowEdits", "1",
    Tk_Offset(DNATrace, show_edits), 0, NULL},
    {TK_CONFIG_BOOLEAN, "-showtrace", "showTrace", "ShowTrace", "1",
    Tk_Offset(DNATrace, show_trace), 0, NULL},
    {TK_CONFIG_BOOLEAN, "-showconf", "showConf", "ShowConf", "1",
    Tk_Offset(DNATrace, show_conf), 0, NULL},
    {TK_CONFIG_FONT, "-font", "font", "Font", DEF_TRACE_FONT,
    Tk_Offset(DNATrace, font), 0, NULL},
    {TK_CONFIG_FONT, "-conf_font", "confFont", "ConfFont", DEF_TRACE_CONF_FONT,
    Tk_Offset(DNATrace, conf_font), 0, NULL},
    {TK_CONFIG_INT, "-showends", "showEnds", "ShowEnds", "1",
    Tk_Offset(DNATrace, show_ends), 0, NULL},
    {TK_CONFIG_INT, "-line_width", "lineWidth", "LineWidth", "0",
         Tk_Offset(DNATrace, line_width), 0, NULL},
    {TK_CONFIG_DOUBLE, "-xmag", "xmag", "Xmag", "1.58",
    Tk_Offset(DNATrace, scale_x), 0, NULL},
    {TK_CONFIG_DOUBLE, "-conf_scale", "confScale", "ConfScale", "1.00",
    Tk_Offset(DNATrace, scale_conf), 0, NULL},
    {TK_CONFIG_INT, "-trace_scale", "traceScale", "traceScale", "0",
    Tk_Offset(DNATrace, trace_scale), 0, NULL},
    {TK_CONFIG_INT, "-style", "style", "Style", "0",
    Tk_Offset(DNATrace, style), 0, NULL},
    {TK_CONFIG_END, (char *) NULL, (char *) NULL, (char *) NULL, (char *) NULL,
	 0, 0, NULL},
};


static int TraceCmd(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv);
void TraceEventProc(ClientData clientData, XEvent *eventPtr);
static int TraceConfigure(Tcl_Interp *interp, DNATrace *tracePtr,
			  int argc, char **argv, int flags);
static int TraceWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv);
void TraceDisplay(ClientData clientData);
void TraceDestroy(char *clientData);
static void TraceScroll(DNATrace *t, char *str);
static int set_visible_region(DNATrace *t);
static int resample(DNATrace *t, int spacing);
static int save_image(DNATrace *t, char *fname, int height,
		      int start, int end);

/*
 * Trace init
 */
int Trace_Init(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "dnatrace", TraceCmd,
		      NULL, (Tcl_CmdDeleteProc *)NULL);
    return TCL_OK;
}

/*
 * Our class command procedure.
 */
static int TraceCmd(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {
    DNATrace *tracePtr;
    Tk_Window tkwin;

    if (argc < 2) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], " pathName ?options?\"", (char *)NULL);
	return TCL_ERROR;
    }

    if (NULL == (tkwin = Tk_CreateWindowFromPath(interp, Tk_MainWindow(interp),
						 argv[1], (char *)NULL))) {
	return TCL_ERROR;
    }
    Tk_SetClass(tkwin, "DnaTrace");

    if (NULL == (tracePtr = (DNATrace *) xmalloc(sizeof(DNATrace)))) {
	Tk_DestroyWindow(tkwin);
	return TCL_ERROR;
    }

    /* printf("Create  0x%p\n", tracePtr); */

    tracePtr->tkwin = tkwin;
    tracePtr->display = Tk_Display(tkwin);
    tracePtr->interp = interp;
    tracePtr->width = 0;
    tracePtr->height = 0;
    tracePtr->cursor_pos = 0;
    tracePtr->cursor_pos_old = 0;
    tracePtr->read = NULL;
    tracePtr->Agc = NULL;
    tracePtr->Cgc = NULL;
    tracePtr->Ggc = NULL;
    tracePtr->Tgc = NULL;
    tracePtr->CursorGC = NULL;
    tracePtr->CutoffGC = NULL;
    tracePtr->VectorGC = NULL;
    tracePtr->QualVecGC = NULL;
    tracePtr->ConfGC = NULL;
    tracePtr->ConfNegGC = NULL;
    tracePtr->Acol = NULL;
    tracePtr->Ccol = NULL;
    tracePtr->Gcol = NULL;
    tracePtr->Tcol = NULL;
    tracePtr->CursorCol = NULL;
    tracePtr->CutoffCol = NULL;
    tracePtr->VectorCol = NULL;
    tracePtr->ConfCol = NULL;
    tracePtr->ConfNegCol = NULL;
    tracePtr->border = NULL;
    tracePtr->borderWidth = 0;
    tracePtr->relief = TK_RELIEF_FLAT;
    tracePtr->disp_width = 300;
    tracePtr->disp_offset = 0;
    tracePtr->disp_offset_old = -1;
    tracePtr->last_disp_width = -1;
    tracePtr->last_disp_offset = -1;
    tracePtr->last_disp_height[0] = -1;
    tracePtr->last_disp_height[1] = -1;
    tracePtr->last_disp_height[2] = -1;
    tracePtr->last_disp_height[3] = -1;
    tracePtr->scale_x = 1.58;
    tracePtr->scale_y = 1.0;
    tracePtr->scale_conf = 1.0;
    tracePtr->trace_scale = 0;
    tracePtr->xScrollCmd = NULL;
    tracePtr->old_scroll_mode = 0;
    tracePtr->flags = 0;
    tracePtr->show_numbers = 1;
    tracePtr->show_sequence  = 1;
    tracePtr->show_edits = 1;
    tracePtr->show_trace = 1;
    tracePtr->show_conf = 1;
    tracePtr->font = NULL;
    tracePtr->conf_font = NULL;
    tracePtr->Ned = 0;
    tracePtr->MaxNed = 0;
    tracePtr->edBases = NULL;
    tracePtr->edPos = NULL;
    tracePtr->edConf = NULL;
    tracePtr->pn = 0;
    tracePtr->ps = 0;
    tracePtr->pe = 0;
    tracePtr->pt = 0;
    tracePtr->comp = 0;
    tracePtr->show_ends = 1;
    tracePtr->line_width = 0; /* default == fast */
    tracePtr->leftVector = -1;
    tracePtr->rightVector = -1;
    tracePtr->tracePos = NULL;
    tracePtr->tracePosE = NULL;
    tracePtr->visible = VisibilityFullyObscured;
    tracePtr->yticks = 0;

    tracePtr->ps_options.orientation = NULL;
    tracePtr->ps_options.font_type = NULL;
    tracePtr->ps_trace.basePos_lookup = NULL;
    tracePtr->ps_trace.options_A.colour_str = NULL;
    tracePtr->ps_trace.options_A.XCol = NULL;
    tracePtr->ps_trace.options_A.dash_str = NULL;
    tracePtr->ps_trace.options_A.dash = NULL;
    tracePtr->ps_trace.options_C.colour_str = NULL;
    tracePtr->ps_trace.options_C.XCol = NULL;
    tracePtr->ps_trace.options_C.dash_str = NULL;
    tracePtr->ps_trace.options_C.dash = NULL;
    tracePtr->ps_trace.options_G.colour_str = NULL;
    tracePtr->ps_trace.options_G.XCol = NULL;
    tracePtr->ps_trace.options_G.dash_str = NULL;
    tracePtr->ps_trace.options_G.dash = NULL;
    tracePtr->ps_trace.options_T.colour_str = NULL;
    tracePtr->ps_trace.options_T.XCol = NULL;
    tracePtr->ps_trace.options_T.dash_str = NULL;
    tracePtr->ps_trace.options_T.dash = NULL;
    tracePtr->ps_trace.options_N.colour_str = NULL;
    tracePtr->ps_trace.options_N.XCol = NULL;
    tracePtr->ps_trace.options_N.dash_str = NULL;
    tracePtr->ps_trace.options_N.dash = NULL;

    Tk_CreateEventHandler(tkwin, ExposureMask | StructureNotifyMask |
			  VisibilityChangeMask,
			  TraceEventProc, (ClientData)tracePtr);
    Tcl_CreateCommand(interp, Tk_PathName(tkwin),
		      TraceWidgetCmd, (ClientData)tracePtr,
		      (Tcl_CmdDeleteProc *)NULL);
    if (TraceConfigure(interp, tracePtr, argc-2, argv+2, 0) != TCL_OK) {
	Tk_DestroyWindow(tracePtr->tkwin);
	return TCL_ERROR;
    }

    Tcl_SetResult(interp, Tk_PathName(tkwin), TCL_VOLATILE);
    return TCL_OK;
}

/*
 * Configure.
 */
static int TraceConfigure(Tcl_Interp *interp, DNATrace *tracePtr,
			  int argc, char **argv, int flags) {
    Tk_FontMetrics fm;

    if (Tk_ConfigureWidget(interp, tracePtr->tkwin, configSpecs,
			   argc, argv, (char *)tracePtr, flags) != TCL_OK) {
	return TCL_ERROR;
    }
    Tk_GetFontMetrics(tracePtr->font, &tracePtr->fm);
    tracePtr->font_width = Tk_TextWidth(tracePtr->font, "0", 1);

    Tk_GetFontMetrics(tracePtr->conf_font, &fm);
    tracePtr->conf_font_width = Tk_TextWidth(tracePtr->conf_font, "00", 2);
    tracePtr->conf_font_height = fm.ascent;

    Tk_SetWindowBackground(tracePtr->tkwin,
			   Tk_3DBorderColor(tracePtr->border)->pixel);

    if (tracePtr->Agc != NULL) {
	Tk_FreeGC(tracePtr->display, tracePtr->Agc);
	Tk_FreeGC(tracePtr->display, tracePtr->Cgc);
	Tk_FreeGC(tracePtr->display, tracePtr->Ggc);
	Tk_FreeGC(tracePtr->display, tracePtr->Tgc);
	Tk_FreeGC(tracePtr->display, tracePtr->CursorGC);
	Tk_FreeGC(tracePtr->display, tracePtr->CutoffGC);
	Tk_FreeGC(tracePtr->display, tracePtr->VectorGC);
	Tk_FreeGC(tracePtr->display, tracePtr->QualVecGC);
	Tk_FreeGC(tracePtr->display, tracePtr->ConfGC);
	Tk_FreeGC(tracePtr->display, tracePtr->ConfNegGC);
	tracePtr->Agc = NULL;
    }

    if (tracePtr->Agc == NULL) {
	XGCValues gcv;
	long mask = GCFont | GCFunction | GCForeground | GCBackground
	    | GCLineWidth;

	gcv.function = GXcopy;
	gcv.background = WhitePixelOfScreen(Tk_Screen(tracePtr->tkwin));
	gcv.foreground = tracePtr->Acol->pixel;
	gcv.font = Tk_FontId(tracePtr->font);
	gcv.line_width = tracePtr->line_width;
	tracePtr->Agc = Tk_GetGC(tracePtr->tkwin, mask, &gcv);

	gcv.foreground = tracePtr->Ccol->pixel;
	tracePtr->Cgc = Tk_GetGC(tracePtr->tkwin, mask, &gcv);

	gcv.foreground = tracePtr->Gcol->pixel;
	tracePtr->Ggc = Tk_GetGC(tracePtr->tkwin, mask, &gcv);

	gcv.foreground = tracePtr->Tcol->pixel;
	tracePtr->Tgc = Tk_GetGC(tracePtr->tkwin, mask, &gcv);

	gcv.foreground = tracePtr->CursorCol->pixel;
	tracePtr->CursorGC = Tk_GetGC(tracePtr->tkwin, mask, &gcv);

	mask &= ~GCLineWidth;

	gcv.foreground = tracePtr->CutoffCol->pixel;
	tracePtr->CutoffGC = Tk_GetGC(tracePtr->tkwin, mask, &gcv);

	gcv.font = Tk_FontId(tracePtr->conf_font);
	gcv.foreground = tracePtr->ConfCol->pixel;
	gcv.line_width = 0;
	tracePtr->ConfGC = Tk_GetGC(tracePtr->tkwin, mask, &gcv);

	gcv.font = Tk_FontId(tracePtr->conf_font);
	gcv.foreground = tracePtr->ConfNegCol->pixel;
	gcv.line_width = 0;
	tracePtr->ConfNegGC = Tk_GetGC(tracePtr->tkwin, mask, &gcv);

	gcv.foreground = tracePtr->VectorCol->pixel;
	gcv.background = tracePtr->CutoffCol->pixel;
	tracePtr->VectorGC = Tk_GetGC(tracePtr->tkwin, mask, &gcv);

	tracePtr->stipple = XCreateBitmapFromData(
	        Tk_Display(tracePtr->tkwin),
		RootWindowOfScreen(Tk_Screen(tracePtr->tkwin)),
                (char *)stripe4_bits, stripe4_width, stripe4_height);
	gcv.stipple = tracePtr->stipple;
	gcv.fill_style = FillStippled;
	mask |= GCStipple | GCFillStyle;
	tracePtr->QualVecGC = Tk_GetGC(tracePtr->tkwin, mask, &gcv);
    }

    Tk_GeometryRequest(tracePtr->tkwin, tracePtr->width, tracePtr->height);

    tracePtr->pos[TRACEP_N].bd = 0;
    tracePtr->pos[TRACEP_S].bd = 0;
    tracePtr->pos[TRACEP_E].bd = 0;
    tracePtr->pos[TRACEP_T].bd = 0;

    tracePtr->pos[TRACEP_N].h =
	tracePtr->show_numbers  ? tracePtr->fm.linespace : 0;
    tracePtr->pos[TRACEP_S].h =
	tracePtr->show_sequence ? tracePtr->fm.linespace : 0;
    tracePtr->pos[TRACEP_E].h =
	tracePtr->show_edits    ? tracePtr->fm.linespace : 0;

    tracePtr->flags |= TRACE_REDRAW | TRACE_BORDER | TRACE_RESCALE;
    if (!(tracePtr->flags & TRACE_WAITING)) {
	tracePtr->flags |= TRACE_WAITING;
	Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
    }

    /* calculate new display offset/width if xmag changes */
    {
	int width = tracePtr->width - 2*tracePtr->borderWidth;

	tracePtr->disp_width = (int)(width / tracePtr->scale_x + .999);
    }


    return TCL_OK;
}

/*
 * Trace widget command.
 */
static int TraceWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv) {
    DNATrace *tracePtr = (DNATrace *)clientData;
    int result = TCL_OK;

    if (argc < 2) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], " option ?arg arg ...?\"",
			 (char *) NULL);
	return TCL_ERROR;
    }

    Tcl_Preserve((ClientData) tracePtr);

    if (strcmp(argv[1], "configure") == 0 ||
	strcmp(argv[1], "config") == 0) {
	if (argc == 2) {
	    result = Tk_ConfigureInfo(interp, tracePtr->tkwin,
				      configSpecs, (char *) tracePtr,
				      (char *) NULL, 0);
	} else if (argc == 3) {
	    result = Tk_ConfigureInfo(interp, tracePtr->tkwin,
				      configSpecs, (char *) tracePtr,
				      argv[2], 0);
	} else {
	    result = TraceConfigure(interp, tracePtr, argc-2,
				     argv+2, TK_CONFIG_ARGV_ONLY);
	}
    } else if (strcmp(argv[1], "load") == 0) {
	const char *nativepath;
	Tcl_DString dstr;

	if (argc != 3 && argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " load filename ?format?\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	tracePtr->disp_offset = 0;
	Tcl_DStringInit(&dstr);
	nativepath = Tcl_UtfToExternalDString(NULL, argv[2], -1, &dstr);

	if (-1 == trace_load(tracePtr, (char *)nativepath,
			     argc == 4 ? argv[3] : "anytr")) {
	    Tcl_AppendResult(interp, "Failed to load trace", NULL);

	    Tcl_DStringFree(&dstr);
	    result = TCL_ERROR;
	    goto release;
	}

	Tcl_DStringFree(&dstr);

	tracePtr->flags |= TRACE_REDRAW;
	if (!(tracePtr->flags & TRACE_WAITING)) {
	    tracePtr->flags |= TRACE_WAITING;
	    Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
	}

    } else if (strcmp(argv[1], "save") == 0) {
	const char *nativepath;
	Tcl_DString dstr;

	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " save filename format\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	Tcl_DStringInit(&dstr);
	nativepath = Tcl_UtfToExternalDString(NULL, argv[2], -1, &dstr);

	if (-1 == trace_save(tracePtr, (char *)nativepath, argv[3])) {
	    Tcl_AppendResult(interp, "Failed to load trace", NULL);

	    Tcl_DStringFree(&dstr);
	    result = TCL_ERROR;
	    goto release;
	}

	Tcl_DStringFree(&dstr);

    } else if (strcmp(argv[1], "loaded") == 0) {
	char buf[10];

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " loaded\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	sprintf(buf, "%d", tracePtr->read ? 1 : 0);
	Tcl_AppendResult(interp, buf, NULL);

    } else if (strcmp(argv[1], "info") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " info\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if (tracePtr->read ) {
	    char  ibuf[256];
	    Read* rs = tracePtr->read;
	    if( rs->info )
		Tcl_AppendResult(interp, rs->info, NULL);
	    sprintf( ibuf, "POINTS=%d\nBASES=%d\nBASELINE=%d\nMAXVAL=%d",
		rs->NPoints, rs->NBases, rs->baseline, (int)rs->maxTraceVal );
            Tcl_AppendResult( interp, ibuf, NULL );
	}
    } else if (strcmp(argv[1], "xview") == 0) {
	int offset = 0, centre = 0;
	double f1, f2;
	int type, count;

	if (!tracePtr->read)
	    goto release;

	if (argc == 2) { /* Query */
	    if (tracePtr->disp_width == 0) {
		Tcl_SetResult(interp, "0 1", TCL_STATIC);
	    } else {
		char buf[1024];
		f1 = tracePtr->disp_offset / (double)tracePtr->read->NPoints;
		f2 = (tracePtr->disp_offset + tracePtr->disp_width) /
		    (double)tracePtr->read->NPoints;
		sprintf(buf, "%g %g", f1, f2);
		Tcl_SetResult(interp, buf, TCL_VOLATILE);
	    }
	    goto release;

	} else if (argc == 3) { /* old syntax */
	    char *ap = argv[2];

	    if (*ap == 'C') {
		centre = 1;
		ap++;
	    }
	    if (*ap == '%') {
		offset = trace_get_pos(tracePtr, atoi(ap+1));
	    } else {
		offset = atoi(ap);
	    }
	    if (centre)
		offset -= tracePtr->disp_width / 2;
	    tracePtr->old_scroll_mode = 1;

	} else { /* new syntax */
	    type = Tk_GetScrollInfo(interp, argc, argv, &f1, &count);
	    switch (type) {
	    case TK_SCROLL_ERROR:
		goto release;
	    case TK_SCROLL_MOVETO:
		offset = (int) (f1 * tracePtr->read->NPoints + 0.5);
		break;
	    case TK_SCROLL_PAGES:
		offset = tracePtr->disp_offset +
		    count * 0.9 * tracePtr->disp_width;
		break;
	    case TK_SCROLL_UNITS:
		offset = tracePtr->disp_offset + count * 5;
		break;
	    }
	    tracePtr->old_scroll_mode = 0;
	}

	if (!tracePtr->show_ends) {
	    if (offset < 0)
		offset = 0;
	    else if (offset > tracePtr->read->NPoints - tracePtr->disp_width)
		offset = tracePtr->read->NPoints - tracePtr->disp_width;
	}

	if (tracePtr->disp_offset != offset) {
	    if (tracePtr->disp_offset_old == -1)
		tracePtr->disp_offset_old = tracePtr->disp_offset;
	    tracePtr->disp_offset = offset;
	    tracePtr->flags |= TRACE_SCROLL;
	    if (!(tracePtr->flags & TRACE_WAITING)) {
		tracePtr->flags |= TRACE_WAITING;
		Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
	    }
	}

    } else if (strcmp(argv[1], "xmag") == 0) {
	double new;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " xmag scale\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if (tracePtr->read) {
	    int width = Tk_Width(tracePtr->tkwin) - 2*tracePtr->borderWidth;
	    new = atof(argv[2]);

	    if (new < 0)
		new = 0;
	    else if (new > 10000)
		new = 10000;

	    if (tracePtr->style == STYLE_PYRO)
		new *= 10;
	    else if (tracePtr->style == STYLE_STICK)
		new *= 30;

	    /* scale_x size when all samples are visible */
	    /*
	     * min is the minimum pixels/sample ratio.
	     * The magic 10 below is maximum pixels/sample which corresponds
	     * to 1000 on our 0-1000 scale.
	     */
	    /*
	    min = ((double)width) / tracePtr->read->NPoints;
	    new = new * (10 - min) / 1000 + min;
	    */
	    new = (new+100)/(100+100);

	    /* Shift disp_offset to keep centre of screen stationary */
	    tracePtr->disp_offset -= ((int)((float)(width / new) + .999) -
				      tracePtr->disp_width) / 2;
	    /* Now update tracePtr records */
	    tracePtr->disp_width =  (int)((float)(width / new) + .999);
	    tracePtr->scale_x = new;

	    if (!tracePtr->show_ends) {
		if (tracePtr->disp_offset + tracePtr->disp_width
		    > tracePtr->read->NPoints)
		    tracePtr->disp_offset =
			tracePtr->read->NPoints - tracePtr->disp_width;
		if (tracePtr->disp_offset < 0)
		    tracePtr->disp_offset = 0;
	    }

	    tracePtr->flags |= TRACE_REDRAW;
	    if (!(tracePtr->flags & TRACE_WAITING)) {
		tracePtr->flags |= TRACE_WAITING;
		Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
	    }
	}

    } else if (strcmp(argv[1], "ymag") == 0) {
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " ymag scale\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	tracePtr->scale_y = atof(argv[2]);
	tracePtr->flags |= TRACE_REDRAW;
	if (!(tracePtr->flags & TRACE_WAITING)) {
	    tracePtr->flags |= TRACE_WAITING;
	    Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
	}

    } else if (strcmp(argv[1], "icursor") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " icursor index\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if (argc == 2) {
	    char buf[10];

	    sprintf(buf, "%d", tracePtr->cursor_pos);
	    Tcl_AppendResult(interp, buf, NULL);
	} else {
	    if (tracePtr->read == NULL) {
		result = TCL_OK;
		goto release;
	    }

	    tracePtr->cursor_pos_old = tracePtr->cursor_pos;

	    if (argv[2][0] == '@') {
		tracePtr->cursor_pos = pixel_to_base(tracePtr,
						     atoi(&argv[2][1]), 0);
	    } else {
		tracePtr->cursor_pos = atoi(argv[2]);
	    }

	    if (tracePtr->cursor_pos < 0)
		tracePtr->cursor_pos = 0;
	    if (tracePtr->cursor_pos > tracePtr->Ned)
		tracePtr->cursor_pos = tracePtr->Ned;

	    tracePtr->flags |= TRACE_CURSOR;
	    TraceDisplay((ClientData)tracePtr);
/*
	    if (!(tracePtr->flags & TRACE_WAITING)) {
		tracePtr->flags |= TRACE_WAITING;
		Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
	    }
*/
	}

    } else if (strcmp(argv[1], "insert") == 0) {
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " insert base\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	trace_insert(tracePtr, tracePtr->cursor_pos, argv[2][0]);

	tracePtr->flags |= TRACE_REDRAW;
	if (!(tracePtr->flags & TRACE_WAITING)) {
	    tracePtr->flags |= TRACE_WAITING;
	    Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
	}

    } else if (strcmp(argv[1], "delete") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " delete\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}
	trace_delete(tracePtr, tracePtr->cursor_pos);

	tracePtr->flags |= TRACE_REDRAW;
	if (!(tracePtr->flags & TRACE_WAITING)) {
	    tracePtr->flags |= TRACE_WAITING;
	    Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
	}

    } else if (strcmp(argv[1], "left_cutoff") == 0) {
	if (argc != 3 && argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " left_cutoff ?index?\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if (tracePtr->read) {
	    if (argc == 2) {
		vTcl_SetResult(interp, "%d", tracePtr->read->leftCutoff);
	    } else {
		if (argv[2][0] == '@') {
		    tracePtr->read->leftCutoff =
			pixel_to_base(tracePtr, atoi(&argv[2][1]), 1);
		} else {
		    tracePtr->read->leftCutoff = atoi(argv[2]);
		}

		if (tracePtr->read->leftCutoff >= tracePtr->read->rightCutoff
		    && tracePtr->read->rightCutoff)
		    tracePtr->read->leftCutoff = tracePtr->read->rightCutoff-1;

		tracePtr->flags |= TRACE_REDRAW;
		if (!(tracePtr->flags & TRACE_WAITING)) {
		    tracePtr->flags |= TRACE_WAITING;
		    Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
		}
	    }
	}

    } else if (strcmp(argv[1], "resample") == 0) {
	int spacing;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " resample samples_per_base\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	Tcl_GetInt(interp, argv[2], &spacing);
	resample(tracePtr, spacing);

    } else if (strcmp(argv[1], "right_cutoff") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " right_cutoff ?index?\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if (tracePtr->read) {
	    if (argc == 2) {
		vTcl_SetResult(interp, "%d", tracePtr->read->rightCutoff);
	    } else {
		if (argv[2][0] == '@') {
		    tracePtr->read->rightCutoff =
			pixel_to_base(tracePtr,atoi(&argv[2][1]), 1) + 1;
		} else {
		    tracePtr->read->rightCutoff = atoi(argv[2]);
		}

		if (tracePtr->read->rightCutoff <= tracePtr->read->leftCutoff
		    && tracePtr->read->leftCutoff
		    && tracePtr->read->rightCutoff)
		    tracePtr->read->rightCutoff = tracePtr->read->leftCutoff+1;
		else if (tracePtr->read->rightCutoff == tracePtr->Ned + 1)
		    tracePtr->read->rightCutoff = 0;

		tracePtr->flags |= TRACE_REDRAW;
		if (!(tracePtr->flags & TRACE_WAITING)) {
		    tracePtr->flags |= TRACE_WAITING;
		    Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
		}
	    }
	}

    } else if (strcmp(argv[1], "left_vector") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " left_vector ?index?\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if (tracePtr->read) {
	    if (argc == 2) {
		vTcl_SetResult(interp, "%d", tracePtr->leftVector);
	    } else {
		if (argv[2][0] == '@') {
		    tracePtr->leftVector =
			pixel_to_base(tracePtr, atoi(&argv[2][1]), 1);
		} else {
		    tracePtr->leftVector = atoi(argv[2]);
		}

		if (tracePtr->leftVector >= tracePtr->rightVector &&
		    tracePtr->rightVector != -1)
		    tracePtr->leftVector = tracePtr->rightVector - 1;

		tracePtr->flags |= TRACE_REDRAW;
		if (!(tracePtr->flags & TRACE_WAITING)) {
		    tracePtr->flags |= TRACE_WAITING;
		    Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
		}
	    }
	}

    } else if (strcmp(argv[1], "right_vector") == 0) {
	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " right_vector ?index?\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if (tracePtr->read) {
	    if (argc == 2) {
		vTcl_SetResult(interp, "%d", tracePtr->rightVector);
	    } else {
		if (argv[2][0] == '@') {
		    tracePtr->rightVector =
			pixel_to_base(tracePtr, atoi(&argv[2][1]), 1) + 1;
		} else {
		    tracePtr->rightVector = atoi(argv[2]) + 1;
		    if (tracePtr->rightVector == 0)
			tracePtr->rightVector = -1;
		}

		if (tracePtr->rightVector <= tracePtr->leftVector &&
		    tracePtr->leftVector != -1 && tracePtr->rightVector != -1)
		    tracePtr->rightVector = tracePtr->leftVector + 1;
		else if (tracePtr->rightVector >= tracePtr->Ned + 2)
		    tracePtr->rightVector = tracePtr->Ned + 1;

		tracePtr->flags |= TRACE_REDRAW;
		if (!(tracePtr->flags & TRACE_WAITING)) {
		    tracePtr->flags |= TRACE_WAITING;
		    Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
		}
	    }
	}

    } else if (strcmp(argv[1], "complement") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " complement\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	complement_trace(tracePtr);

	tracePtr->flags |= TRACE_REDRAW;
	if (!(tracePtr->flags & TRACE_WAITING)) {
	    tracePtr->flags |= TRACE_WAITING;
	    Tcl_DoWhenIdle(TraceDisplay, (ClientData)tracePtr);
	}

    } else if (strcmp(argv[1], "position") == 0) {
	char buf[20];

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " position type\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	/* Displays the display offset and width in sample points or
	 * bases.
	 */
	if (argv[2][0] == 's') {
	    sprintf(buf, "%d %d", tracePtr->disp_offset, tracePtr->disp_width);
	    Tcl_AppendResult(interp, buf, NULL);
	} else if (argv[2][0] == 'o') {
	    if (tracePtr->read) {
		int st, en, off;

		off = tracePtr->disp_offset >= 0 ? tracePtr->disp_offset : 0;
		st = tracePtr->tracePos[off];
		en = tracePtr->tracePos[off + tracePtr->disp_width];
		sprintf(buf, "%d %d", st, en - st);
		Tcl_AppendResult(interp, buf, NULL);
	    }
	} else if (argv[2][0] == 'e') {
	    if (tracePtr->read) {
		int st, en, off;

		off = tracePtr->disp_offset >= 0 ? tracePtr->disp_offset : 0;
		st = tracePtr->tracePosE[off];
		en = off  + tracePtr->disp_width;
		if (en >= tracePtr->read->NPoints)
		    en = tracePtr->read->NPoints - 1;
		en = tracePtr->tracePosE[en];
		sprintf(buf, "%d %d", st, en - st);
		Tcl_AppendResult(interp, buf, NULL);
	    }
	}

    } else if (strcmp(argv[1], "sequence") == 0) {
	char tmp[4096];

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " sequence\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}
/*
	strncpy(tmp, tracePtr->read->base, tracePtr->read->NBases);
	tmp[tracePtr->read->NBases] = 0;
*/
	if (tracePtr->edBases) {
	    strncpy(tmp, tracePtr->edBases, tracePtr->Ned);
	    tmp[tracePtr->Ned] = 0;

	    Tcl_AppendResult(interp, tmp, NULL);
	} else {
	    Tcl_ResetResult(interp);
	}

    } else if (strcmp(argv[1], "format") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " format\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	Tcl_SetResult(interp, trace_type_int2str(tracePtr->read->format),
		      TCL_STATIC);

    } else if (strcmp(argv[1], "orig_format") == 0) {
	int fm;

	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " orig_format\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	fm = tracePtr->read->orig_trace_format;
	if (TT_UNK == fm)
	    fm = tracePtr->read->format;
	Tcl_SetResult(interp, trace_type_int2str(fm), TCL_STATIC);

    } else if (strcmp(argv[1], "flash") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " flash\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	trace_flash(tracePtr);

    } else if (strcmp(argv[1], "ps_configure") == 0) {
	if(TCL_ERROR == ps_configure(&(tracePtr->ps_options), argc - 1, argv + 1)) {
	    result = TCL_ERROR;
	    goto release;
	}

    } else if (strcmp(argv[1], "ps_configure_trace") == 0) {
	if(TCL_ERROR == ps_configure_trace(tracePtr, argc - 1, argv + 1)) {
	    result = TCL_ERROR;
	    goto release;
	}

    } else if (strcmp(argv[1], "ps_configure_trace_line") == 0) {
	if (argc < 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " ps_configure_trace_line [A|C|G|T] <command> <arguments>\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if(TCL_ERROR == ps_configure_trace_line(tracePtr, argc - 2, argv + 2)) {
	    result = TCL_ERROR;
	    goto release;
	}

    } else if (strcmp(argv[1], "print") == 0) {
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " print file_name\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if(-1 == trace_print(tracePtr, argv[2])) {
	    result = TCL_ERROR;
	    goto release;
	}

    } else if (strcmp(argv[1], "visible_region") == 0) {
	return set_visible_region(tracePtr);

    } else if (strcmp(argv[1], "save_image") == 0) {
	int start_sample, end_sample;
	int height = tracePtr->pos[TRACEP_T].h;

	if (argc != 6 && argc != 3 && argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " save_image filename "
			     "?height ?start_sample end_sample??\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if (argc >= 4)
	    Tcl_GetInt(interp, argv[3], &height);

	if (argc == 6) {
	    Tcl_GetInt(interp, argv[4], &start_sample);
	    Tcl_GetInt(interp, argv[5], &end_sample);
	} else {
	    double avg_spacing =
		(tracePtr->read->basePos[tracePtr->read->NBases-1] -
		 tracePtr->read->basePos[0]) /
		(double)tracePtr->read->NBases;
	    int base_sample;

	    if (tracePtr->cursor_pos > 0 &&
		tracePtr->cursor_pos < tracePtr->Ned)
		base_sample = trace_get_pos(tracePtr, tracePtr->cursor_pos);
	    else
		base_sample = tracePtr->disp_offset + tracePtr->disp_width/2;

	    start_sample = base_sample - 5 * avg_spacing;
	    end_sample = base_sample + 5 * avg_spacing;
	}

	if (-1 == save_image(tracePtr, argv[2], height,
			     start_sample, end_sample)) {
	    result = TCL_ERROR;
	    goto release;
	}

    } else if (strcmp(argv[1], "base_info") == 0) {
	int index, conf, base;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " base_info index\"",
			     (char *)NULL);
	    result = TCL_ERROR;
	    goto release;
	}

	if (tracePtr->read == NULL) {
	    result = TCL_OK;
	    goto release;
	}

	if (argv[2][0] == '@') {
	    int pos = atoi(&argv[2][1]);
	    index = pixel_to_base(tracePtr, pos, 0);

	    /*
	     * Correct for the fact that pixel_to_base doesn't give the
	     * nearest base, but the next base.
	     */
	    pos = pixel_to_point(tracePtr, pos - tracePtr->borderWidth - 1);
	    if (index < tracePtr->read->NBases-1 && index > 0) {
		int left = abs(pos - tracePtr->read->basePos[index-1]);
		int mid  = abs(pos - tracePtr->read->basePos[index]);
		int right= abs(pos - tracePtr->read->basePos[index+1]);

		if (left < mid && left < right) {
		    index--;
		} else if (right < mid && right < left) {
		    index++;
		}
	    } else if (index <= 0) {
		index = 0;
	    } else  /* index >= tracePtr->read->NBases - 1 */ {
		index = tracePtr->read->NBases;
	    }

	} else {
	    index = atoi(argv[2]);
	}

	if (index < 0 || index >= tracePtr->read->NBases) {
	    result = TCL_OK;
	    goto release;
	}

	base = tracePtr->read->base[index];
	switch (base) {
	case 'a':
	case 'A':
	    conf = tracePtr->read->prob_A[index];
	    break;
	case 'c':
	case 'C':
	    conf = tracePtr->read->prob_C[index];
	    break;
	case 'g':
	case 'G':
	    conf = tracePtr->read->prob_G[index];
	    break;
	case 't':
	case 'T':
	    conf = tracePtr->read->prob_T[index];
	    break;
	case '-':
	case 'N':
	    conf = (tracePtr->read->prob_A[index] +
		    tracePtr->read->prob_C[index] +
		    tracePtr->read->prob_G[index] +
		    tracePtr->read->prob_T[index]) / 4;
	    break;
	default:
	    conf = 0;
	    break;
	}

	vTcl_SetResult(interp, "Position %4d, Base %c, Confidence %3d (1/%g)",
		       index+1, base, conf, pow(10, -conf/10.0));

    } else {
	Tcl_AppendResult(interp, "unknown option \"", argv[1], "\".", NULL);
	result = TCL_ERROR;
    }

 release:
    Tcl_Release((ClientData)tracePtr);

    return result;
}

/*
 * Event handling
 */
void TraceEventProc(ClientData clientData, XEvent *eventPtr) {
    DNATrace *tracePtr = (DNATrace *) clientData;

    if (eventPtr->type == Expose) {
	tracePtr->flags |= TRACE_REDRAW | TRACE_BORDER;
	if (!(tracePtr->flags & TRACE_WAITING)) {
	    tracePtr->flags |= TRACE_WAITING;
	    Tcl_DoWhenIdle(TraceDisplay, (ClientData) tracePtr);
	}

    } else if (eventPtr->type == ConfigureNotify) {
	/* resize */
	int width = Tk_Width(tracePtr->tkwin) - 2*tracePtr->borderWidth;

	tracePtr->disp_width =  (int)((float)(width/tracePtr->scale_x) + .999);
	if (tracePtr->read && tracePtr->disp_offset + tracePtr->disp_width >
	    tracePtr->read->NPoints && !tracePtr->show_ends) {
	    tracePtr->disp_offset = tracePtr->read->NPoints -
		tracePtr->disp_width;
	    if (tracePtr->disp_offset < 0)
		tracePtr->disp_offset = 0;
	}

#if 0
	/*
	 * It seems we cope with display width showing more samples than the
	 * trace has, so disable this for now.
	 * (It causes problems with the positioning of automatically
	 * generated difference traces that are invoked from the editor.)
	 */
	if (tracePtr->read &&
	    tracePtr->disp_width >= tracePtr->read->NPoints) {
	    tracePtr->disp_width = tracePtr->read->NPoints -1;
	    tracePtr->scale_x = ((double)width)	/ tracePtr->disp_width;
	}
#endif

	tracePtr->flags |= TRACE_REDRAW | TRACE_BORDER | TRACE_RESCALE;
	if (!(tracePtr->flags & TRACE_WAITING)) {
	    tracePtr->flags |= TRACE_WAITING;
	    Tcl_DoWhenIdle(TraceDisplay, (ClientData) tracePtr);
	}

    } else if (eventPtr->type == DestroyNotify) {
	Tcl_DeleteCommand(tracePtr->interp,
			  Tk_PathName(tracePtr->tkwin));
	tracePtr->tkwin = NULL;
	if (tracePtr->flags & TRACE_WAITING)
	    Tcl_CancelIdleCall(TraceDisplay,
			       (ClientData) tracePtr);
	Tcl_EventuallyFree((ClientData) tracePtr, TraceDestroy);

    } else if (eventPtr->type == VisibilityNotify) {
	if (tracePtr->visible == VisibilityFullyObscured) {
	    /*
	     * Ensure it's redrawn. An expose event may will not be
	     * generated if we have backing store for the window.
	     */
	    tracePtr->flags |= TRACE_REDRAW | TRACE_BORDER;
	    if (!(tracePtr->flags & TRACE_WAITING)) {
		tracePtr->flags |= TRACE_WAITING;
		Tcl_DoWhenIdle(TraceDisplay, (ClientData) tracePtr);
	    }
	}
	tracePtr->visible = eventPtr->xvisibility.state;
	/* Reset visibility optimisation */
	tracePtr->disp_offset_old = -1;
    }
}

static void get_pixmap(Display *d, Drawable w, DNATrace *tracePtr, int num,
		       Pixmap p, int disp_offset, int disp_width,
		       Tk_Window tkwin) {
    int ql, qr, l0, r0, t, sl, sr;

    sl = ql = l0 = point_to_pixel(tracePtr, disp_offset-.3);
    sr = qr = r0 = point_to_pixel(tracePtr, disp_offset + disp_width+.3);

    if (tracePtr->read->leftCutoff) {
	if (tracePtr->read->leftCutoff == tracePtr->Ned)
	    t = tracePtr->read->NPoints;
	else
	    t = (trace_get_pos(tracePtr, tracePtr->read->leftCutoff - 1) +
		 trace_get_pos(tracePtr, tracePtr->read->leftCutoff)) / 2;

	if (t > disp_offset) {
	    ql = point_to_pixel(tracePtr, t-.3);

	    if (ql > r0)
		ql = r0;
	}
    }

    if (tracePtr->read->rightCutoff) {
	if (tracePtr->read->rightCutoff == 1)
	    t = 0;
	else
	    t = ((trace_get_pos(tracePtr, tracePtr->read->rightCutoff-1) +
		  trace_get_pos(tracePtr, tracePtr->read->rightCutoff-2)) / 2);
	if (t < disp_offset + disp_width) {
	    qr = point_to_pixel(tracePtr, t+.3);

	    if (qr < l0)
		qr = l0;
	    else if (qr > r0)
		qr = r0;
	}
    }

    if (tracePtr->leftVector != -1) {
	if (tracePtr->leftVector == tracePtr->Ned)
	    t = tracePtr->read->NPoints;
	else
	    t = (trace_get_pos(tracePtr, tracePtr->leftVector - 1) +
		 trace_get_pos(tracePtr, tracePtr->leftVector)) / 2;

	if (t > disp_offset) {
	    sl = point_to_pixel(tracePtr, t-.3);

	    if (sl > r0)
		sl = r0;
	}
    }

    if (tracePtr->rightVector != -1) {
	if (tracePtr->rightVector == 1)
	    t = 0;
	else
	    t = ((trace_get_pos(tracePtr, tracePtr->rightVector-1) +
		  trace_get_pos(tracePtr, tracePtr->rightVector-2)) / 2);
	if (t < disp_offset + disp_width) {
	    sr = point_to_pixel(tracePtr, t+.3);

	    if (sr < l0)
		sr = l0;
	    else if (sr > r0)
		sr = r0;
	}
    }

    /*
     * Clear the background to a colour based on the clip information.
     * This is horrendously complicated due to the need to stipple when
     * quality and vector clips overlap. No guarantees the logic is correct
     * here - it's pretty tricky to dream up all the possible cases.
     */
    {
	int hstart1=l0, hend1=l0;
	int hstart2=r0, hend2=r0;
	int hstart3=r0, hend3=r0;
	int qstart1=l0, qend1=l0;
	int qstart2=r0, qend2=r0;
	int vstart1=l0, vend1=l0;
	int vstart2=r0, vend2=r0;
	int bstart1=l0, bend1=r0;

	if (ql != l0) {
	    qstart1 = l0;
	    qend1 = ql;
	    bstart1 = ql;
	}

	if (qr != r0) {
	    qstart2 = qr;
	    qend2 = r0;
	    bend1 = qr;
	}

	if (sl != l0) {
	    if (ql > sl) {
		hend1 = sl;
		qstart1 = sl;
	    } else {
		bstart1 = sl;
		hend1 = ql;
		qend1 = l0;
		vstart1 = ql;
		if (sl < qr) {
		    vend1 = sl;
		} else {
		    vend1 = qr;
		    hstart2 = qr;
		    hend2 = sl;
		    qstart2 = sl;
		}
	    }
	}

	if (sr != r0) {
	    if (sr > qr) {
		if (hstart2 < qr && hend2 > qr) {
		    hend2 = r0;
		    qstart2 = r0;
		} else {
		    hstart3 = sr;
		    qend2 = sr;
		}
	    } else {
		bend1 = sr;
		hstart3 = qr;
		qstart2 = r0;
		vend2 = qr;
		if (sr >= ql) {
		    vstart2 = sr;
		} else {
		    vstart2 = ql;
		    hstart2 = sr;
		    hend2 = ql;
		    qend1 = sr;
		}
	    }
	}

/*
	printf("\nsr=%3d,qr=%3d,r0=%3d sl=%3d,ql=%3d,l0=%3d\n",
	       sr, qr, r0, sl, ql, l0);
	printf("h1 %03d-%03d\n", hstart1, hend1);
	printf("h2 %03d-%03d\n", hstart2, hend2);
	printf("h3 %03d-%03d\n", hstart3, hend3);
	printf("v1 %03d-%03d\n", vstart1, vend1);
	printf("v2 %03d-%03d\n", vstart2, vend2);
	printf("q1 %03d-%03d\n", qstart1, qend1);
	printf("q2 %03d-%03d\n", qstart2, qend2);
	printf("b1 %03d-%03d\n", bstart1, bend1);
*/

	/*
	 * Now we have several ranges of possible states.
	 * qstart/qend detail where quality colour needs to be drawn.
	 * vstart/vend detail where vector colour needs to be drawn.
	 * hstart/hend detail where we need to stipple both colours.
	 * finally bstart/bend details whatever is left.
	 */

	if (tracePtr->pos[num].h > 0) {
	    if (hstart1 != hend1) {
		XFillRectangle(d, p, tracePtr->CutoffGC,
			       hstart1, 0,
			       hend1 - hstart1, tracePtr->pos[num].h);
		XFillRectangle(d, p, tracePtr->QualVecGC,
			       hstart1, 0,
			       hend1 - hstart1, tracePtr->pos[num].h);
	    }

	    if (hstart2 != hend2) {
		XFillRectangle(d, p, tracePtr->CutoffGC,
			       hstart2, 0,
			       hend2 - hstart2, tracePtr->pos[num].h);
		XFillRectangle(d, p, tracePtr->QualVecGC,
			       hstart2, 0,
			       hend2 - hstart2, tracePtr->pos[num].h);
	    }

	    if (hstart3 != hend3) {
		XFillRectangle(d, p, tracePtr->CutoffGC,
			       hstart3, 0,
			       hend3 - hstart3, tracePtr->pos[num].h);
		XFillRectangle(d, p, tracePtr->QualVecGC,
			       hstart3, 0,
			       hend3 - hstart3, tracePtr->pos[num].h);
	    }

	    if (qstart1 != qend1) {
		XFillRectangle(d, p, tracePtr->CutoffGC,
			       qstart1, 0,
			       qend1 - qstart1, tracePtr->pos[num].h);
	    }

	    if (qstart2 != qend2) {
		XFillRectangle(d, p, tracePtr->CutoffGC,
			       qstart2, 0,
			       qend2 - qstart2, tracePtr->pos[num].h);
	    }

	    if (vstart1 != vend1) {
		XFillRectangle(d, p, tracePtr->VectorGC,
			       vstart1, 0,
			       vend1 - vstart1, tracePtr->pos[num].h);
	    }

	    if (vstart2 != vend2) {
		XFillRectangle(d, p, tracePtr->VectorGC,
			       vstart2, 0,
			       vend2 - vstart2, tracePtr->pos[num].h);
	    }

	    if (bstart1 < bend1) {
		Tk_Fill3DRectangle(tkwin, p, tracePtr->border,
				   bstart1, 0, bend1 - bstart1,
				   tracePtr->pos[num].h, 0, 0);
	    }
	}
    }
}



static void set_pixmap(Display *d, Pixmap p, Drawable w, DNATrace *tracePtr,
		       int num) {
    if (tracePtr->pos[num].h > 0)
	XCopyArea(d, p, w, tracePtr->Agc,
		  0, 0, Tk_Width(tracePtr->tkwin) - 2*tracePtr->borderWidth,
		  tracePtr->pos[num].h, tracePtr->borderWidth,
		  tracePtr->pos[num].y);
}


static void resize_pixmaps(DNATrace *tracePtr) {
    int width;

    if (tracePtr->pn) {
	Tk_FreePixmap(tracePtr->display, tracePtr->pn);
	tracePtr->pn = 0;
    }

    if (tracePtr->ps) {
	Tk_FreePixmap(tracePtr->display, tracePtr->ps);
	tracePtr->ps = 0;
    }

    if (tracePtr->pe) {
	Tk_FreePixmap(tracePtr->display, tracePtr->pe);
	tracePtr->pe = 0;
    }

    if (tracePtr->pt) {
	Tk_FreePixmap(tracePtr->display, tracePtr->pt);
	tracePtr->pt = 0;
    }

    tracePtr->width = Tk_Width(tracePtr->tkwin);
    tracePtr->height = Tk_Height(tracePtr->tkwin);

    width = tracePtr->width - 2 * tracePtr->borderWidth;

    if (tracePtr->show_numbers && tracePtr->pos[TRACEP_N].h > 0)
	tracePtr->pn = Tk_GetPixmap(tracePtr->display,
				     Tk_WindowId(tracePtr->tkwin),
				     width,
				     tracePtr->pos[TRACEP_N].h,
				     Tk_Depth(tracePtr->tkwin));

    if (tracePtr->show_sequence && tracePtr->pos[TRACEP_S].h > 0)
	tracePtr->ps = Tk_GetPixmap(tracePtr->display,
				     Tk_WindowId(tracePtr->tkwin),
				     width,
				     tracePtr->pos[TRACEP_S].h,
				     Tk_Depth(tracePtr->tkwin));

    if (tracePtr->show_edits && tracePtr->pos[TRACEP_E].h > 0)
	tracePtr->pe = Tk_GetPixmap(tracePtr->display,
				     Tk_WindowId(tracePtr->tkwin),
				     width,
				     tracePtr->pos[TRACEP_E].h,
				     Tk_Depth(tracePtr->tkwin));

    if (tracePtr->show_trace && tracePtr->pos[TRACEP_T].h > 0)
	tracePtr->pt = Tk_GetPixmap(tracePtr->display,
				     Tk_WindowId(tracePtr->tkwin),
				     width,
				     tracePtr->pos[TRACEP_T].h,
				     Tk_Depth(tracePtr->tkwin));
}

/*
 * Display
 */
void TraceDisplay(ClientData clientData) {
    DNATrace *tracePtr = (DNATrace *)clientData;
    int width, height, bd = tracePtr->borderWidth;
    Display *d;
    Drawable w;
    int disp_width;
    int disp_offset;

    tracePtr->flags &= ~TRACE_WAITING;

    if (tracePtr->visible == VisibilityFullyObscured) {
	/* Delay redraw until some later time */
	return;
    }

    if (!Tk_IsMapped(tracePtr->tkwin)) {
	Tcl_CancelIdleCall(TraceDisplay, (ClientData)tracePtr);
	return;
    }

    if (!tracePtr->font)
	return;

    width = Tk_Width(tracePtr->tkwin);
    height = Tk_Height(tracePtr->tkwin);
    d = tracePtr->display;
    w = Tk_WindowId(tracePtr->tkwin);
    disp_width = tracePtr->disp_width;
    disp_offset = tracePtr->disp_offset;

    /* Rescale */
    if (tracePtr->flags & TRACE_RESCALE) {
	int i, t = 0;

	/* Calculate new trace display size */
	for( i=0; i<TRACEP_T; i++ )
	    t += tracePtr->pos[i].h;

	tracePtr->pos[TRACEP_T].h = Tk_Height(tracePtr->tkwin) - t - 2*bd;

	tracePtr->pos[TRACEP_T].bd = 0;
	if (tracePtr->pos[TRACEP_T].h <= 0)
	    tracePtr->pos[TRACEP_T].h = 0;

	tracePtr->pos[TRACEP_T].y = bd;
	tracePtr->pos[TRACEP_S].y = tracePtr->pos[TRACEP_T].y + tracePtr->pos[TRACEP_T].h;
	tracePtr->pos[TRACEP_E].y = tracePtr->pos[TRACEP_S].y + tracePtr->pos[TRACEP_S].h;
	tracePtr->pos[TRACEP_N].y = tracePtr->pos[TRACEP_E].y + tracePtr->pos[TRACEP_E].h;

	/* Resize the pixmaps */
	resize_pixmaps(tracePtr);

	/* Force redraw of everything else */
	tracePtr->flags &= ~TRACE_RESCALE;
	tracePtr->flags |= (TRACE_REDRAW | TRACE_BORDER);
    }

    if (tracePtr->flags & TRACE_BORDER) {
	/* Draw the border, and redraw all */
	Tk_Fill3DRectangle(tracePtr->tkwin, w, tracePtr->border, 0, 0,
			   width, height, tracePtr->borderWidth, tracePtr->relief);
	tracePtr->flags |= TRACE_REDRAW;
	tracePtr->flags &= ~TRACE_BORDER;
    }

    width -= 2*tracePtr->borderWidth;

    if (tracePtr->flags & TRACE_SCROLL && !(tracePtr->flags & TRACE_REDRAW)) {
	int old_x = point_to_pixel(tracePtr, tracePtr->disp_offset_old);
	int new_x = point_to_pixel(tracePtr, tracePtr->disp_offset);

	if (ABS(new_x - old_x) < width && tracePtr->disp_offset_old != -1) {
	    int from_x, to_x;

	    if (new_x > old_x) {
		/* scroll right */
		int new_o;

		from_x = new_x - old_x;
		to_x = 0;
		width = width - (new_x - old_x) + 10;

		new_o = disp_offset + disp_width -
		    (disp_offset - tracePtr->disp_offset_old) - 1;
		disp_width -= (new_o - disp_offset);
		disp_offset = new_o;

	    } else {
		/* scroll left */
		from_x = 0;
		to_x = old_x - new_x;

		width = width - (old_x - new_x) + 10;
		disp_width = tracePtr->disp_offset_old -
		    tracePtr->disp_offset + 1;
	    }

	    /* Scroll remaining screen */
	    if (tracePtr->show_numbers && tracePtr->pn)
		XCopyArea(d, tracePtr->pn, tracePtr->pn, tracePtr->Agc,
			  from_x, 0, width, tracePtr->pos[TRACEP_N].h,
			  to_x, 0);

	    if (tracePtr->show_sequence && tracePtr->ps)
		XCopyArea(d, tracePtr->ps, tracePtr->ps, tracePtr->Agc,
			  from_x, 0, width, tracePtr->pos[TRACEP_S].h,
			  to_x, 0);

	    if (tracePtr->show_edits && tracePtr->pe)
		XCopyArea(d, tracePtr->pe, tracePtr->pe, tracePtr->Agc,
			  from_x, 0, width, tracePtr->pos[TRACEP_E].h,
			  to_x, 0);

	    if (tracePtr->show_trace && tracePtr->pt)
		XCopyArea(d, tracePtr->pt, tracePtr->pt, tracePtr->Agc,
			  from_x, 0, width, tracePtr->pos[TRACEP_T].h,
			  to_x, 0);
	}

	tracePtr->flags &= ~TRACE_SCROLL;
	tracePtr->flags |= TRACE_REDRAW;
    }

    if (tracePtr->flags & TRACE_REDRAW) {
	if (tracePtr->read) {
	    int doff, dwid, doit;

	    /* We draw a little bit to the left and right of our actual region
	     * just to ensure we overlap the correct data when scrolling.
	     */
	    doff = disp_offset > 0 ? disp_offset - 1 : 0;
	    if (doff + disp_width + 2 > tracePtr->read->NPoints) {
		dwid = tracePtr->read->NPoints - doff;
	    } else {
		dwid = disp_width + 2;
	    }

	    /*
	     * If we're scrolling at the end of a trace then there's no need
	     * to redraw anything, but we do need to copy the pixmaps still.
	     */
	    doit = (dwid > 0 && doff < tracePtr->read->NPoints);

	    if (tracePtr->show_numbers && tracePtr->pos[TRACEP_N].h > 0) {
		if (tracePtr->last_disp_offset != disp_offset ||
		    tracePtr->last_disp_width != disp_width ||
		    tracePtr->last_disp_height[TRACEP_N] !=
		    tracePtr->pos[TRACEP_N].h) {

		    get_pixmap(d, w, tracePtr, TRACEP_N, tracePtr->pn,
			       disp_offset, disp_width, tracePtr->tkwin);
		    if (doit && tracePtr->Ned > 0)
			trace_draw_numbers(tracePtr, d, tracePtr->pn,
					   doff, disp_width + 2,
					   0, tracePtr->pos[TRACEP_N].h);
		    tracePtr->last_disp_height[TRACEP_N] =
			tracePtr->pos[TRACEP_N].h;
		}
		set_pixmap(d, tracePtr->pn, w, tracePtr, TRACEP_N);
	    }

	    if (tracePtr->show_sequence && tracePtr->pos[TRACEP_S].h > 0) {
		if (tracePtr->last_disp_offset != disp_offset ||
		    tracePtr->last_disp_width != disp_width ||
		    tracePtr->last_disp_height[TRACEP_S] !=
		    tracePtr->pos[TRACEP_S].h) {

		    get_pixmap(d, w, tracePtr, TRACEP_S, tracePtr->ps,
			       disp_offset, disp_width, tracePtr->tkwin);
		    if (doit && tracePtr->Ned > 0)
			trace_draw_sequence(tracePtr, d, tracePtr->ps,
					    doff, disp_width + 2,
					    0, tracePtr->pos[TRACEP_S].h);
		    tracePtr->last_disp_height[TRACEP_S] =
			tracePtr->pos[TRACEP_S].h;
		}
		set_pixmap(d, tracePtr->ps, w, tracePtr, TRACEP_S);
	    }

	    if (tracePtr->show_edits && tracePtr->pos[TRACEP_E].h > 0) {
		if (tracePtr->last_disp_offset != disp_offset ||
		    tracePtr->last_disp_width != disp_width ||
		    tracePtr->last_disp_height[TRACEP_E] !=
		    tracePtr->pos[TRACEP_E].h) {

		    get_pixmap(d, w, tracePtr, TRACEP_E, tracePtr->pe,
			       disp_offset, disp_width, tracePtr->tkwin);
		    if (doit && tracePtr->Ned > 0)
			trace_draw_edits(tracePtr, d, tracePtr->pe,
					 doff, disp_width + 2,
					 0, tracePtr->pos[TRACEP_E].h);
		    tracePtr->last_disp_height[TRACEP_E] =
			tracePtr->pos[TRACEP_E].h;
		}
		set_pixmap(d, tracePtr->pe, w, tracePtr, TRACEP_E);
	    }

	    if (tracePtr->show_trace && tracePtr->pos[TRACEP_T].h > 0) {
		if (tracePtr->last_disp_offset != disp_offset ||
		    tracePtr->last_disp_width != disp_width ||
		    tracePtr->last_disp_height[TRACEP_T] !=
		    tracePtr->pos[TRACEP_T].h) {

		    get_pixmap(d, w, tracePtr, TRACEP_T, tracePtr->pt,
			       disp_offset, disp_width, tracePtr->tkwin);
		    if (doit)
			trace_draw_trace(tracePtr, d, tracePtr->pt,
					 doff, dwid,
					 0, tracePtr->pos[TRACEP_T].h);
		    tracePtr->last_disp_height[TRACEP_T] =
			tracePtr->pos[TRACEP_T].h;
		}
		set_pixmap(d, tracePtr->pt, w, tracePtr, TRACEP_T);
	    }

	    /*
	     * Disable this optimisation at present as it doesn't appear
	     * to work 100% of the time. WHY!?
	     *
	     * tracePtr->last_disp_offset = disp_offset;
	     * tracePtr->last_disp_width = disp_width;
	     */

	    tracePtr->flags &= ~TRACE_REDRAW;
	}

    }

    /* Cursor change */
    if (tracePtr->flags & TRACE_CURSOR) {
	int off;

	if (tracePtr->show_edits) {
	    /*
	     * Horizontal line in the edits display.
	     * This extends from 4 pixels onwards of the previous base
	     * for 8 pixels. We'll convert from trace to pixel and back to
	     * trace again to compute the region to redraw.
	     */
	    int pos, wid;

	    if (tracePtr->cursor_pos_old > 0)
		off = trace_get_pos(tracePtr, tracePtr->cursor_pos_old-1);
	    else
		off = 0;

	    pos = point_to_pixel(tracePtr, off) + 4;
	    off = pixel_to_point(tracePtr, pos) - 1;
	    wid = pixel_to_point(tracePtr, pos+8) - off + 2;

	    /* Cursor is in the sequence line */
	    get_pixmap(d, w, tracePtr, TRACEP_E, tracePtr->pe,
		       off, wid, tracePtr->tkwin);
	    trace_draw_edits(tracePtr, d, tracePtr->pe,
			     off, wid,
			     0, tracePtr->pos[TRACEP_E].h);
	    set_pixmap(d, tracePtr->pe, w, tracePtr, TRACEP_E);

	} else {
	    /*
	     * Vertical line in the trace display.
	     * The physical cursor is drawn 1 pixel before and after the trace
	     * cursor position. Hence we sub/add 1 to the pixel position and
	     * then convert to the prev/next trace pos. This handles both
	     * cases where pixels/point > 1 and < 1.
	     */
	    int pos, wid;

	    off = trace_get_pos(tracePtr, tracePtr->cursor_pos_old);
	    pos = point_to_pixel(tracePtr, off);
	    off = pixel_to_point(tracePtr, pos - 1) - 1;
	    wid = pixel_to_point(tracePtr, pos + 1) - off + 1;

	    get_pixmap(d, w, tracePtr, TRACEP_T, tracePtr->pt,
		       off, wid, tracePtr->tkwin);
	    trace_draw_trace(tracePtr, d, tracePtr->pt,
			     off-1, wid+2,
			     0, tracePtr->pos[TRACEP_T].h);
	    set_pixmap(d, tracePtr->pt, w, tracePtr, TRACEP_T);
	}

	tracePtr->flags &= ~TRACE_CURSOR;
    }

    /* Update scrollbar */
    tracePtr->disp_offset_old = -1;
    TraceScroll(tracePtr, "Display");
}

/*
 * Destroy
 */

void TraceDestroy(char *clientData) {
    DNATrace *tracePtr = (DNATrace *) clientData;

    /* printf("Destroy 0x%p\n", tracePtr); */

    Tcl_CancelIdleCall(TraceDisplay, (ClientData)tracePtr);

    if (tracePtr->pn)
	Tk_FreePixmap(tracePtr->display, tracePtr->pn);
    if (tracePtr->ps)
	Tk_FreePixmap(tracePtr->display, tracePtr->ps);
    if (tracePtr->pe)
	Tk_FreePixmap(tracePtr->display, tracePtr->pe);
    if (tracePtr->pt)
	Tk_FreePixmap(tracePtr->display, tracePtr->pt);

    Tk_FreeOptions(configSpecs, (char *) tracePtr,
		   tracePtr->display, 0);

    if (tracePtr->Agc != NULL)
	Tk_FreeGC(tracePtr->display, tracePtr->Agc);

    if (tracePtr->Cgc != NULL)
	Tk_FreeGC(tracePtr->display, tracePtr->Cgc);

    if (tracePtr->Ggc != NULL)
	Tk_FreeGC(tracePtr->display, tracePtr->Ggc);

    if (tracePtr->Tgc != NULL)
	Tk_FreeGC(tracePtr->display, tracePtr->Tgc);

    if (tracePtr->CursorGC != NULL)
	Tk_FreeGC(tracePtr->display, tracePtr->CursorGC);

    if (tracePtr->CutoffGC != NULL)
	Tk_FreeGC(tracePtr->display, tracePtr->CutoffGC);

    if (tracePtr->VectorGC != NULL)
	Tk_FreeGC(tracePtr->display, tracePtr->VectorGC);

    if (tracePtr->QualVecGC != NULL)
	Tk_FreeGC(tracePtr->display, tracePtr->QualVecGC);

    if (tracePtr->ConfGC != NULL)
	Tk_FreeGC(tracePtr->display, tracePtr->ConfGC);

    if (tracePtr->ConfNegGC != NULL)
	Tk_FreeGC(tracePtr->display, tracePtr->ConfNegGC);

    trace_unload(tracePtr);

    xfree(tracePtr);
}


static void TraceScroll(DNATrace *t, char *str) {
    char string[100];

    if (t->xScrollCmd && t->read) {
	if (t->old_scroll_mode) {
	    sprintf(string, " %d %d %d %d", t->read->NPoints, t->disp_width,
		    t->disp_offset, t->disp_offset + t->disp_width);
	} else {
	    double start, end;
	    start = t->disp_offset / (double)t->read->NPoints;
	    end = (t->disp_offset + t->disp_width) / (double)t->read->NPoints;
	    sprintf(string, " %g %g", start, end);
	}

	if (Tcl_VarEval(t->interp, t->xScrollCmd, string, NULL) != TCL_OK) {
	    Tcl_AddErrorInfo(t->interp,
			     "\n    (xscrollcommand executed by trace)");
	    Tcl_BackgroundError(t->interp);
	}
    }
}

int *trace_index_to_basePos(uint_2 *basePos, int NBases, int NPoints) {
    /*
     * Returns an array of size NPoints. If i is a value of an element in
     * basePos, then basePos_lookup[i] has the value of the index of the
     * basePos element.
     * All other other elements of basePos_lookup are set equal to -1.
     */
    int	i;
    int	*basePos_lookup;

    if (!NPoints)
	return NULL;

    basePos_lookup = (int *) xmalloc(NPoints * sizeof(int));
    if(NULL == basePos_lookup) {
	return NULL;
    }

    for(i = 0; i < NPoints; i++) {
	basePos_lookup[i] = -1;
    }
    for(i = 0; i < NBases; i++) {
	basePos_lookup[MIN(basePos[i], NPoints - 1)] = i;
    }

    return basePos_lookup;
}

static int set_visible_region(DNATrace *t) {
    int min, max, first, last;

    min = 0;
    max =  t->read->NBases - 1;
    first = t->disp_offset;
    first = (first < 0) ? min : t->tracePos[first];
    last = t->disp_offset + t->disp_width - 1;
    last = (last > (t->read->NPoints - 1)) ? max : t->tracePos[last];

    vTcl_SetResult(t->interp, "%d %d %d %d", min, max, first, last);

    return TCL_OK;
}

/*
 * Resamples in X such that each base call consists of 'spacing' trace slices.
 * Assumes that bases are centred on each peak.
 *
 * Returns 0 for success
 *        -1 for failure
 */
static int resample(DNATrace *t, int spacing) {
    Read *r = t->read;
    int nbases = r->NBases;
    TRACE *traceA = NULL;
    TRACE *traceC = NULL;
    TRACE *traceG = NULL;
    TRACE *traceT = NULL;
    int base, nsamples;

    traceA = (TRACE *)xcalloc(nbases * spacing, sizeof(TRACE));
    traceC = (TRACE *)xcalloc(nbases * spacing, sizeof(TRACE));
    traceG = (TRACE *)xcalloc(nbases * spacing, sizeof(TRACE));
    traceT = (TRACE *)xcalloc(nbases * spacing, sizeof(TRACE));
    if (!traceA || !traceC || !traceG || !traceT) {
	goto error;
    }

    /* Foreach base find the peak, centred on that base */
    nsamples = 0;
    for (base = 0; base < r->NBases; base++) {
	double left_samp, right_samp;
	double ratio, xpos;
	int sample;

	left_samp  = (base == 0)
	    ? 0
	    : (r->basePos[base] + r->basePos[base-1])/2.0;
	right_samp = (base == r->NBases - 1)
	    ? r->NPoints-2
	    : (r->basePos[base] + r->basePos[base+1])/2.0;

	if (right_samp > r->NPoints-2)
	    right_samp = r->NPoints-2;

	/* Iterate around samples for this peak */
	if (!r->flow) {
	    xpos = left_samp;
	    ratio = (right_samp - left_samp) / (spacing);
	    for (sample = 0; sample < spacing; sample++) {
		int l;
		double fl, fr;
		l = (int)xpos;
		fl = 1-(xpos-l);
		fr = 1-(l+1-xpos);
		
		traceA[nsamples] = fl * r->traceA[l] + fr * r->traceA[l+1];
		traceC[nsamples] = fl * r->traceC[l] + fr * r->traceC[l+1];
		traceG[nsamples] = fl * r->traceG[l] + fr * r->traceG[l+1];
		traceT[nsamples] = fl * r->traceT[l] + fr * r->traceT[l+1];
		nsamples++;
		xpos += ratio;
	    }
	} else {
	    for (sample = 0; sample < spacing; sample++) {
		int bp = r->basePos[base];
		if (sample == spacing/2) {
		    traceA[nsamples] = bp ? r->traceA[bp] : 0;
		    traceC[nsamples] = bp ? r->traceC[bp] : 0;
		    traceG[nsamples] = bp ? r->traceG[bp] : 0;
		    traceT[nsamples] = bp ? r->traceT[bp] : 0;
		} else {
		    traceA[nsamples] = 0;
		    traceC[nsamples] = 0;
		    traceG[nsamples] = 0;
		    traceT[nsamples] = 0;
		}
		nsamples++;
	    }
	}
    }

    /* Update read structure to point to new arrays */
    for (base = 0; base < r->NBases; base++) {
	r->basePos[base] = base*spacing + spacing/2;
    }

    xfree(r->traceA); r->traceA = traceA;
    xfree(r->traceC); r->traceC = traceC;
    xfree(r->traceG); r->traceG = traceG;
    xfree(r->traceT); r->traceT = traceT;
    r->NPoints = nsamples;

    t->tracePos = (uint_2 *)xrealloc(t->tracePos, r->NPoints * sizeof(uint_2));
    t->tracePosE= (uint_2 *)xrealloc(t->tracePosE,r->NPoints * sizeof(uint_2));

    return 0;

 error:
    if (traceA) xfree(traceA);
    if (traceC) xfree(traceC);
    if (traceG) xfree(traceG);
    if (traceT) xfree(traceT);
    return -1;
}

#ifdef HAVE_PNG
/*
 * Takes an X "drawable" and saves it as a PNG format file to 'fp'.
 * (xoff,yoff)-(xoff+width-1,yoff+height-1) specifies the inclusive area
 * of the drawable to use.
 *
 * The colour palette is chosen from the GCs allocated in the DNATrace
 * pointer.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
int drawable_to_png(DNATrace *t, FILE *fp, Display *disp, Drawable d,
		    int xoff, int yoff, int width, int height) {
    XImage *i = NULL;
    int x, y;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_colorp palette = NULL;
    png_byte *row = NULL;

    /*
     * The pixmap is in the server-side. Create a client side XImage and
     * copy the pixmap over to it.
     */
#ifdef _WIN32
    /*
     * Tk's tkWinImage.c implementation of XGetImage only supports
     * plane_mask == 1. Consequently we have to put up with black and
     * white images under windows.
     */
    i = XGetImage(disp, d, xoff, yoff, width, height, 1, XYPixmap);
#else
    i = XGetImage(disp, d, xoff, yoff, width, height, AllPlanes, XYPixmap);
#endif
    if (!i)
	goto error;

    /* Allocate PNG structs */
    if (!(png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
					    NULL, NULL, NULL)))
	goto error;
    if (!(info_ptr = png_create_info_struct(png_ptr)))
	goto error;

    /* PNG requires setjmp/longjmp to be used */
    if (setjmp(png_jmpbuf(png_ptr)))
	goto error;

    /* Create the info bits */
    png_init_io(png_ptr, fp);

    png_set_IHDR(png_ptr,
		 info_ptr,
		 width,
		 height,
		 4 /* bit_depth */,
		 PNG_COLOR_TYPE_PALETTE,
		 PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_DEFAULT,
		 PNG_FILTER_TYPE_DEFAULT);

    /* Create the palette */
    palette = (png_colorp)png_malloc(png_ptr, 8 * sizeof(*palette));
#define SET_PALETTE(palette,xcol) \
	(palette).red   = (xcol)->red   / 256; \
	(palette).green = (xcol)->green / 256; \
	(palette).blue  = (xcol)->blue  / 256;
    palette[0].red   = 255;
    palette[0].green = 255;
    palette[0].blue  = 255;
    SET_PALETTE(palette[1], t->Acol);
    SET_PALETTE(palette[2], t->Ccol);
    SET_PALETTE(palette[3], t->Gcol);
    SET_PALETTE(palette[4], t->Tcol);
    SET_PALETTE(palette[5], t->CursorCol);
    SET_PALETTE(palette[6], t->ConfCol);
    SET_PALETTE(palette[7], t->ConfNegCol);
    png_set_PLTE(png_ptr, info_ptr, palette, 8);

    /* Write PNG info, including palette */
    png_write_info(png_ptr, info_ptr);

    /*
     * Convert the XImage to a series of PNG rows
     * Lookup colours by checking against our GCs we're drawing with.
     */
    row = (png_byte *)xmalloc((width/2+1) * sizeof(*row));
    for (y = 0; y < height; y++) {
	memset(row, 0, width/2);
	for (x = 0; x < width; x++) {
	    long pixel = XGetPixel(i, x, y);
	    int value;

	    if (pixel == t->Acol->pixel) {
		value = 1;
	    } else if (pixel == t->Ccol->pixel) {
		value = 2;
	    } else if (pixel == t->Gcol->pixel) {
		value = 3;
	    } else if (pixel == t->Tcol->pixel) {
		value = 4;
	    } else if (pixel == t->CursorCol->pixel) {
		value = 5;
	    } else if (pixel == t->ConfCol->pixel) {
		value = 6;
	    } else if (pixel == t->ConfNegCol->pixel) {
		value = 7;
	    } else {
		value = 0;
	    }

	    if (x%2 == 0) {
		row[x/2] = 16*value;
	    } else {
		row[x/2] += value;
	    }
	}
	png_write_row(png_ptr, row);
    }
    xfree(row);
    png_write_end(png_ptr, NULL);

    png_free(png_ptr, palette);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    XDestroyImage(i);

    return 0;

 error:
    if (i)
	XDestroyImage(i);

    if (png_ptr)
	png_destroy_write_struct(&png_ptr, &info_ptr);

    return -1;
}
#endif


/*
 * Generates an image in PPM format of the trace from sample coordinates
 * 'start' to 'end' inclusive and saves it to 'filename'.
 *
 * The string returned is allocated by Tcl_Alloc and deallocation is the
 * responsibility of the caller.
 *
 * Returns 0 for success
 *	  -1 for failure
 */
static int save_image(DNATrace *t, char *filename, int trace_height,
		      int start, int end) {
#ifdef HAVE_PNG
    Pixmap p = 0;
    int width;
    int height;
    int tmp_disp_offset;
    FILE *fp = NULL;
    Drawable rootwin;
    int depth;

    if (NULL == (fp = fopen(filename, "wb")))
	goto error;

    width = point_to_pixel(t, end) - point_to_pixel(t, start) + 1;
    height = trace_height + t->pos[TRACEP_S].h + 2;

    /*
     * Allocate a pixmap and draw to it.
     * Normally we'd use t->tkwin and Tk_Depth(t->tkwin) in these arguments,
     * but we need this to work when the window has not yet been mapped, so
     * we compute the drawable and depth from the root window.
     */
    rootwin = RootWindow(t->display, DefaultScreen(t->display));
    depth = DefaultDepth(t->display, DefaultScreen(t->display));
    p = Tk_GetPixmap(t->display,
		     rootwin,
		     width,
		     height,
		     depth);
    if (!p)
	goto error;

    /* Draw into our pixmap, temporarily moving the disp offset */
    XFillRectangle(t->display, p, t->CutoffGC, 0, 0, width, height);
    tmp_disp_offset = start;
    t->disp_offset = start;
    trace_draw_trace(t, t->display, p, start, end-start+1,
		     0, trace_height);
    trace_draw_sequence(t, t->display, p, start, end-start+1,
			trace_height + 2, t->pos[TRACEP_S].h);
    t->disp_offset = tmp_disp_offset;


    if (-1 == drawable_to_png(t, fp, t->display, p, 0, 0, width, height))
	goto error;

    fclose(fp);
    Tk_FreePixmap(t->display, p);

    return 0;

 error:
    if (fp)
	fclose(fp);

    if (p)
	Tk_FreePixmap(t->display, p);
#endif

    return -1;
}


