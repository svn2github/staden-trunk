/*
 * The EdNames tk widget.
 * As with the Editor widget this is an inheritance of the tkSheet widget.
 *
 * See tkEditor for better commenting of this method.
 */

#include <tk.h>
#include <string.h>
#include "tk_defs.h"
#include "tkSheet_common.h"
#include "tkSeqedNames.h"
#include "tcl_utils.h"
#include "xalloc.h"

#undef TKSHEET
#define TKSHEET(en)   ((tkSheet *)(en))

/* ---- Local function prototypes ---- */
static int SeqedNamesCmd(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv);
static int SeqedNamesWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, char **argv);
static int SeqedNamesConfigure(Tcl_Interp *interp, seqedNames *sn,
			  int argc, char **argv, int flags);

#define offset(field) Tk_Offset(seqedNames, field)
static Tk_ConfigSpec configSpecs[] = {
#   include "tkSheet_config.h"
    {TK_CONFIG_END,
	 (char *)NULL,	(char *)NULL,	(char *)NULL,	(char *) NULL,
         0,	0,	NULL},
};
#undef offset

/* ---- Global functions ---- */

/*
 * Init.
 */
int SeqedNames_Init(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "seqednames", SeqedNamesCmd,
                      NULL,
                      (Tcl_CmdDeleteProc *)NULL);

    return TCL_OK;
}

/* ---- Local functions ---- */

/*
 * Our class command procedure. This is different from usual class command
 * procs as we actually use the SheetCmd() procedure and then carefully
 * subvert their actions to our own uses.
 */
static int SeqedNamesCmd(ClientData clientData, Tcl_Interp *interp,
			 int argc, char **argv) {
    Tk_Window tkwin;
    seqedNames *sn;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " pathName ?options?\"", (char *)NULL);
        return TCL_ERROR;
    }

    /*
     * Allocate
     */
    if (NULL == (sn = (seqedNames *)xmalloc(sizeof(seqedNames)))) {
        return TCL_ERROR;
    }

    /*
     * Call common sheet initialisation code
     */
    tkwin = SheetCmdCommon(interp, Tk_MainWindow(interp),
			   configSpecs, (tkSheet *)sn,
			   argv[1], "SeqedNames");
    if (NULL == tkwin) {
	xfree(sn);
	return TCL_ERROR;
    }

    /*
     * Register our instance of the widget class
     */
    Tcl_CreateCommand(interp, Tk_PathName(tkwin),
                      SeqedNamesWidgetCmd, (ClientData)sn,
                      (Tcl_CmdDeleteProc *)NULL);


    /*
     * And process our arguments - send them to Configure
     */
    if (SeqedNamesConfigure(interp, sn, argc-2, argv+2, 0) != TCL_OK) {
        Tk_DestroyWindow(tkwin);
        return TCL_ERROR;
    }


    Tcl_SetResult(interp, Tk_PathName(tkwin), TCL_VOLATILE);
    return TCL_OK;
}





/*
 * Configure.
 */
static int SeqedNamesConfigure(Tcl_Interp *interp, seqedNames *sn,
			  int argc, char **argv, int flags) {

    if (TCL_ERROR == SheetConfigureCommon(interp, (tkSheet *)sn,
					  argc, argv, flags))
	return TCL_ERROR;

    return TCL_OK;
}


/*
 * seqedNames widget command
 */
static int SeqedNamesWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			       int argc, char **argv) {
    seqedNames *sn = (seqedNames *)clientData;
    int result = TCL_OK;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " option ?arg arg ...?\"",
                         (char *) NULL);
        return TCL_ERROR;
    }

    Tcl_Preserve((ClientData)TKSHEET(sn));

    if ('c' == *argv[1] && strcmp(argv[1], "configure") == 0) {
	result = SheetWidgetCmdConfig(interp, TKSHEET(sn), argc, argv);
    } else {
	Tcl_AppendResult(interp, "bad option \"", argv[1], "\": must be FIXME",
			 NULL);
	goto fail;
    }

    Tcl_Release((ClientData)TKSHEET(sn));
    return result;

 fail:
    Tcl_Release((ClientData)TKSHEET(sn));
    return TCL_ERROR;
}
