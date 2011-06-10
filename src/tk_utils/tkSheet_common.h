#ifndef _TK_SHEET_COMMON_H
#define _TK_SHEET_COMMON_H

#include "tkSheet.h"

/*
 * JOB numbers for the extension function
 */
#define SHEET_JOB_RESIZE 0
#define SHEET_JOB_DESTROY 1

#define SHEET_REDRAW_PENDING 	0x01
#define SHEET_REDRAW_BORDER	0x02
#define SHEET_REDRAW_TEXT	0x04
#define SHEET_REDRAW_ALL 	(SHEET_REDRAW_BORDER | SHEET_REDRAW_TEXT)
#define SHEET_DESTROYED 0x08

#define DEF_SHEET_WIDTH  "80"
#define DEF_SHEET_HEIGHT "1"
/* #define DEF_SHEET_FONT   "Fixed -20 roman" */
#define DEF_SHEET_FONT   "sheet_font"
#define DEF_SHEETB_FONT  "sheet_bold_font"

/*
 * Initialise and create the window. This routine should be used from within
 * the initial SheetCmd type of procedures.
 *
 * Returns a Tk_Window upon success,
 *           NULL for failure
 */
Tk_Window SheetCmdCommon(Tcl_Interp *interp, Tk_Window main,
			 Tk_ConfigSpec configSpecs[],
			 tkSheet *sw, char *path, char *class);

/*
 * The common sheet Configure procedure. This calls Tk_ConfigureWidget with
 * the appropriate Tk_ConfigSpecs and then updates any necessary elements of
 * the sheet widget in question.
 *
 * Returns TCL_OK or TCL_ERROR
 */
int SheetConfigureCommon(Tcl_Interp *interp, tkSheet *sw,
			 int argc, char **argv, int flags);

/*
 * Handles the configure and cget widget command options.
 *
 * Returns TCL_OK or TCL_ERROR
 */
int SheetWidgetCmdConfig(Tcl_Interp *interp, tkSheet *sw,
			 int argc, char **argv);

/*
 * Redisplay the sheet widget
 */
void SheetDisplay(ClientData clientData);


/*
 * Parsing arguments dealing with sheet coordinates. '@val' is in pixels,
 * otherwise 'val' is in character units.
 */
void sheet_arg_x(tkSheet *sw, char *arg, int *val);
void sheet_arg_y(tkSheet *sw, char *arg, int *val);


/*
 * Reconfigures the height.
 */
void sheet_set_display_height(tkSheet *tsw, int height);


/*
 * Drawing the beveled separator
 */
void sheet_draw_separator(tkSheet *sw, int position);

#endif
