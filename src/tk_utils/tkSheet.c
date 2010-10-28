/*
 * The Sheet widget for tk. This is simply a tk interface to the sheet.c file.
 */

#include <string.h>
#include <tk.h>
#include "tkSheet.h"
#include "tk_defs.h"
#include "xalloc.h"
#include "tcl_utils.h"
#include "tkSheet_common.h"

/* ---- Local defines ---- */
#define offset(field) Tk_Offset(tkSheet, field)

/* ---- Local function prototypes ---- */
int SheetCmd(ClientData clientData, Tcl_Interp *interp,
	     int argc, char **argv);

static int SheetConfigure(Tcl_Interp *interp, tkSheet *sw,
			  int argc, char **argv, int flags);

static int SheetWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv);

/* ---- Local static data ---- */
static Tk_ConfigSpec configSpecs[] = {
#   include "tkSheet_config.h"
    {TK_CONFIG_END,
	 (char *)NULL,	(char *)NULL,	(char *)NULL,	(char *) NULL,
         0,	0,	NULL},
};

#undef offset

static tkSheet *last_sheet = NULL;

/* ---- Globally callable functions ---- */

/*
 * Init.
 */
int Sheet_Init(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "sheet", SheetCmd,
                      NULL, (Tcl_CmdDeleteProc *)NULL);

    return TCL_OK;
}


/*
 * Provides an interface for C to grab the recently allocated tkSheet pointer.
 */
tkSheet *LastSheet() {
    return last_sheet;
}


/* ---- Locally callable functions ---- */

/*
 * Our class command procedure; should only be used by tkEditor.c
 */
int SheetCmd(ClientData clientData, Tcl_Interp *interp,
	     int argc, char **argv) {
    Tk_Window tkwin;
    tkSheet *sw;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " pathName ?options?\"", (char *)NULL);
        return TCL_ERROR;
    }


    /*
     * Allocate
     */
    if (NULL == (sw = (tkSheet *)xmalloc(sizeof(tkSheet)))) {
        return TCL_ERROR;
    }
    last_sheet = sw;


    /*
     * Call common sheet initialisation code
     */
    tkwin = SheetCmdCommon(interp, Tk_MainWindow(interp), configSpecs, sw,
			   argv[1], "Sheet");
    if (NULL == tkwin) {
	xfree(sw);
	return TCL_ERROR;
    }


    /*
     * Register our instance of the widget class
     */
    Tcl_CreateCommand(interp, Tk_PathName(tkwin),
                      SheetWidgetCmd, (ClientData)sw,
                      (Tcl_CmdDeleteProc *)NULL);


    /*
     * And process our arguments - send them to Configure
     */
    if (SheetConfigure(interp, sw, argc-2, argv+2, 0) != TCL_OK) {
        Tk_DestroyWindow(tkwin);
        return TCL_ERROR;
    }


    Tcl_SetResult(interp, Tk_PathName(tkwin), TCL_VOLATILE);
    return TCL_OK;
}


/*
 * Configure.
 */
static int SheetConfigure(Tcl_Interp *interp, tkSheet *sw,
			  int argc, char **argv, int flags) {
    if (TCL_ERROR == SheetConfigureCommon(interp, sw,
					  argc, argv, flags))
	return TCL_ERROR;

    return TCL_OK;
}


/*
 * widget command
 */
static int SheetWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv) {
    tkSheet *sw = (tkSheet *)clientData;
    int result = TCL_OK;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " option ?arg arg ...?\"",
                         (char *) NULL);
        return TCL_ERROR;
    }

    Tcl_Preserve((ClientData)sw);

    if (strcmp(argv[1], "configure") == 0 ||
	strcmp(argv[1], "config") == 0) {
	result = SheetWidgetCmdConfig(interp, sw, argc, argv);

    } else if (strcmp(argv[1], "add") == 0) {
	int x, y;

	if (argc != 5) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " add x y text\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(sw->interp, argv[2], &x);
	Tcl_GetInt(sw->interp, argv[3], &y);

	XawSheetPutText(&sw->sw, x, y, (Dimension)strlen(argv[4]), argv[4]);
	sw->flags |= SHEET_REDRAW_TEXT;
	if (!(sw->flags & SHEET_REDRAW_PENDING)) {
	    sw->flags |= SHEET_REDRAW_PENDING;
	    Tcl_DoWhenIdle(SheetDisplay, (ClientData)sw);
	}

    } else {
	Tcl_AppendResult(interp, "bad option \"", argv[1], "\": must be "
			 "configure or add", NULL);
	result = TCL_ERROR;
    }

    Tcl_Release((ClientData)sw);
    return result;

 fail:
    Tcl_Release((ClientData)sw);
    return TCL_ERROR;
}
