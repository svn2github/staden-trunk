/*
 * Common utilities to be called by all sheet based widgets.
 */

#include <math.h>
#include <tk.h>
#include "tkSheet.h"
#include "tk_defs.h"
#include "xalloc.h"
#include "tkSheet_common.h"

/* ---- Local defines ---- */

static void SheetEventProc(ClientData clientData, XEvent *eventPtr);

static void SheetDestroy(char *clientData);

/* ---- Globally callable functions ---- */

/*
 * Initialise and create the window. This routine should be used from within
 * the initial SheetCmd type of procedures.
 *
 * Returns a Tk_Window upon success,
 *           NULL for failure
 */
Tk_Window SheetCmdCommon(Tcl_Interp *interp, Tk_Window main,
			 Tk_ConfigSpec configSpecs[],
			 tkSheet *sw, char *path, char *class) {
    Tk_Window tkwin;
    
    /* Create the window */
    tkwin = Tk_CreateWindowFromPath(interp, main, path, (char *)NULL);
    if (NULL == tkwin)
	return NULL;

    /* Set its class */
    Tk_SetClass(tkwin, class);


    /*
     * Initialise variables so that configure works correctly.
     */
    sw->configSpecs = configSpecs;
    sw->sw.tkwin = tkwin;
    sw->sw.display = Tk_Display(tkwin);
    sw->sw.font = NULL;
    sw->sw.bold_font = NULL;
    sw->interp = interp;
    sw->border = NULL;
    sw->fg = NULL;
    sw->light = NULL;
    sw->indel_fg = NULL;
    sw->flags = 0;
    sw->initialised = 0;
    sw->extension = NULL;
    sw->extensionData = NULL;
    sw->divider = 0;
    sw->grid = 0;
#ifdef COUNT_REDRAW_REQUESTS
    sw->count = 0;
#endif

    /*
     * Request events
     */
    Tk_CreateEventHandler(tkwin, ExposureMask | StructureNotifyMask,
                          SheetEventProc, (ClientData)sw);

    return tkwin;
}


/*
 * The common sheet Configure procedure. This calls Tk_ConfigureWidget with
 * the appropriate Tk_ConfigSpecs and then updates any necessary elements of
 * the sheet widget in question.
 *
 * Returns TCL_OK or TCL_ERROR
 */
int SheetConfigureCommon(Tcl_Interp *interp, tkSheet *sw,
			 int argc, char **argv, int flags) {
    int font_width, font_height;

    /*
     * Standard configure options.
     */
    if (Tk_ConfigureWidget(interp, sw->sw.tkwin, sw->configSpecs,
                           argc, argv, (char *)sw, flags) != TCL_OK) {
        return TCL_ERROR;
    }
    Tk_GetFontMetrics(sw->sw.font, &sw->sw.fm);

    Tk_SetBackgroundFromBorder(sw->sw.tkwin, sw->border);

    sw->sw.font_width = font_width = Tk_TextWidth(sw->sw.font, "0", 1);
    font_height = sw->sw.fm.linespace;
    sw->sw.width_in_pixels = sw->sw.columns * font_width +
	2 * (sw->sw.border_width /* + sw->sw.pad_x */);
    sw->sw.height_in_pixels = sw->sw.rows * font_height +
	2 * (sw->sw.border_width /* + sw->sw.pad_y */);
    
    Tk_GeometryRequest(sw->sw.tkwin, sw->sw.width_in_pixels,
		       sw->sw.height_in_pixels);
    Tk_SetInternalBorder(sw->sw.tkwin, sw->sw.border_width);

    /*
     * Are we initialising things? If so we've got things to setup
     * and memory to allocate.
     */
    if (!sw->initialised) {
	sheet_create(&sw->sw, sw->light->pixel, sw->fg->pixel,
		     Tk_3DBorderColor(sw->border)->pixel,
		     sw->indel_fg->pixel);
	sw->initialised = 1;

	if (sw->grid)
	    Tk_SetGrid(sw->sw.tkwin, sw->sw.columns, sw->sw.rows,
		       sw->sw.font_width, sw->sw.fm.linespace);
    } else {
	int old_rows = sw->sw.rows, old_columns = sw->sw.columns;
	sheet_resize(&sw->sw, old_rows, old_columns);
        sheet_config(&sw->sw, sw->light->pixel, sw->fg->pixel,
		     Tk_3DBorderColor(sw->border)->pixel,
		     sw->indel_fg->pixel);
    }

    sw->flags |= SHEET_REDRAW_ALL;
    if (!(sw->flags & SHEET_REDRAW_PENDING)) {
	sw->flags |= SHEET_REDRAW_PENDING;
	Tcl_DoWhenIdle(SheetDisplay, (ClientData)sw);
    }

    return TCL_OK;
}


/*
 * Handles the configure and cget widget command options.
 *
 * Returns TCL_OK or TCL_ERROR
 */
int SheetWidgetCmdConfig(Tcl_Interp *interp, tkSheet *sw,
			 int argc, char **argv) {
    if (argc == 2) {
	return Tk_ConfigureInfo(interp, sw->sw.tkwin,
				sw->configSpecs, (char *)sw,
				(char *) NULL, 0);
    } else if (argc == 3) {
	return Tk_ConfigureInfo(interp, sw->sw.tkwin,
				sw->configSpecs, (char *)sw,
				argv[2], 0);
    } else {
	return SheetConfigureCommon(interp, sw, argc-2,
				    argv+2, TK_CONFIG_ARGV_ONLY);
    }
}


/*
 * Parsing arguments dealing with sheet coordinates. '@val' is in pixels,
 * otherwise 'val' is in character units.
 */
void sheet_arg_x(tkSheet *sw, char *arg, int *val) {
    if (*arg == '@') {
	Tcl_GetInt(sw->interp, &arg[1], val);
	*val = PIXEL_TO_COL(&sw->sw, *val);
    } else {
	Tcl_GetInt(sw->interp, arg, val);
    }
}

void sheet_arg_y(tkSheet *sw, char *arg, int *val) {
    if (*arg == '@') {
	Tcl_GetInt(sw->interp, &arg[1], val);


	*val = PIXEL_TO_ROW(&sw->sw, *val);
    } else {
	Tcl_GetInt(sw->interp, arg, val);
    }
#ifdef DEBUG
    printf("val %d \n", *val);
#endif
    (*val)--;
}



/*
 *---------------------------------------------------------------------------
 * Now that the tk bits are out the way...
 * Some utilities for manipulating the sheet widget from C.
 *---------------------------------------------------------------------------
 */

/*
 * Reconfigures the height.
 */
void sheet_set_display_height(tkSheet *tsw, int height) {
    int font_height, new_height;

    if (height == tsw->sw.rows)
	return;

    font_height = tsw->sw.fm.linespace;
    new_height = height * font_height + 2 * (tsw->sw.border_width );

    Tk_GeometryRequest(tsw->sw.tkwin, tsw->sw.width_in_pixels, new_height);
    Tk_SetInternalBorder(tsw->sw.tkwin, tsw->sw.border_width);
    if (tsw->grid) {
	Tk_UnsetGrid(tsw->sw.tkwin);
	Tk_SetGrid(tsw->sw.tkwin, tsw->sw.columns, height,
		   tsw->sw.font_width, tsw->sw.fm.linespace);
    }
    
    /*
     * We need to resize here anyway (although not display) as the code that
     * calls this will immediately assume the sheet is large enough, even
     * though it hasn't been accept (ie ConfigureNotify) yet. So we resize
     * here to the max of current and new.
     */
    if (tsw->sw.rows < height) {
	int old_rows = tsw->sw.rows;

	tsw->sw.rows = height;
	tsw->sw.height_in_pixels = tsw->sw.rows * font_height +
	    2 * (tsw->sw.border_width);

	sheet_resize(&tsw->sw, old_rows, tsw->sw.columns);
    }
}


/*
 * Drawing the beveled separator
 */
void sheet_draw_separator(tkSheet *sw, int position) {
    int y;
    Drawable w;

    if (!(sw->divider = position)) {
	return;
    }

    y = ROW_TO_PIXEL(&sw->sw, position);
    w = Tk_WindowId(sw->sw.tkwin);

    Tk_3DHorizontalBevel(sw->sw.tkwin, w, sw->border,
			 0, y,
			 sw->sw.width_in_pixels, 2,
			 0, 1, 1, /* leftIn, rightIn, topBevel */
			 sw->relief);
    Tk_3DHorizontalBevel(sw->sw.tkwin, w, sw->border,
			 0, y + 1,
			 sw->sw.width_in_pixels, 1,
			 0, 0, 0, /* leftIn, rightIn, topBevel */
			 sw->relief);
}


/*
 * Display
 */
void SheetDisplay(ClientData clientData) {
    tkSheet *sw = (tkSheet *)clientData;
    Drawable w;

    sw->flags &= ~SHEET_REDRAW_PENDING;

    if (sw->flags & SHEET_DESTROYED)
	return;

    if (!(sw->sw.tkwin)) /* destroyed */
	return;
    
    if (!(w = Tk_WindowId(sw->sw.tkwin))) /* not mapped */
	return;

    /* Border */
    if (sw->flags & SHEET_REDRAW_BORDER) {
	sw->flags &= ~SHEET_REDRAW_BORDER;

	Tk_Draw3DRectangle(sw->sw.tkwin, w, sw->border,
			   0, 0,
			   sw->sw.width_in_pixels,
			   sw->sw.height_in_pixels,
			   sw->sw.border_width, sw->relief);

    }

    if (sw->flags & SHEET_REDRAW_TEXT) {
	sw->flags &= ~SHEET_REDRAW_TEXT;

	sheet_display(&sw->sw);

	if (sw->divider) {
	    int y = ROW_TO_PIXEL(&sw->sw, sw->divider);

	    Tk_3DHorizontalBevel(sw->sw.tkwin, w, sw->border,
				 0, y,
				 sw->sw.width_in_pixels, 1,
				 0, 1, 1, /* leftIn, rightIn, topBevel */
				 sw->relief);
	    Tk_3DHorizontalBevel(sw->sw.tkwin, w, sw->border,
				 0, y + 1,
				 sw->sw.width_in_pixels, 1,
				 0, 0, 0, /* leftIn, rightIn, topBevel */
				 sw->relief);
	}
    }
}


/* ---- Static functions ---- */

/*
 * Event handling
 */
static void SheetEventProc(ClientData clientData, XEvent *eventPtr) {
    tkSheet *sw = (tkSheet *)clientData;

    if (!sw->initialised)
	return;

    switch(eventPtr->type) {
    case Expose:

	sw->flags |= SHEET_REDRAW_ALL;
	if (!(sw->flags & SHEET_REDRAW_PENDING)) {
	    sw->flags |= SHEET_REDRAW_PENDING;
	    Tcl_DoWhenIdle(SheetDisplay, (ClientData)sw);
	}

	break;

    case ConfigureNotify: {
	int old_width, old_height, font_width, font_height;
	double dcol, drow;

#ifdef COUNT_REDRAW_REQUESTS
	sw->count = 0;
#endif

	old_width = sw->sw.columns;
	old_height= sw->sw.rows;
	font_width = Tk_TextWidth(sw->sw.font, "0", 1);
	font_height = sw->sw.fm.linespace;

	dcol = (eventPtr->xconfigure.width - 2*sw->sw.border_width)
	    / (double)font_width;
	drow = (eventPtr->xconfigure.height - 2*sw->sw.border_width)
	    / (double)font_height;

	sw->sw.columns =  ceil(dcol);
	sw->sw.rows    =  ceil(drow);

#if 0
	if (sw->sw.columns == old_width && sw->sw.rows == old_height) {
	    break;
	}
#endif

	sw->sw.width_in_pixels = sw->sw.columns * font_width +
	    2 * (sw->sw.border_width);
	sw->sw.height_in_pixels = sw->sw.rows * font_height +
	    2 * (sw->sw.border_width);

#ifdef DEBUG
	printf("old_ht %d font_ht %d drow %f rows %d pixels %d\n", 
	       old_height, font_height,
	       drow, sw->sw.rows, sw->sw.height_in_pixels );
#endif
	sheet_resize(&sw->sw, old_height, old_width);

	    sheet_clear(&sw->sw);
#if 0
	if (fmod(drow, 1.0) != 0.0) {
	    sheet_clear(&sw->sw);
	    Tk_GeometryRequest(sw->sw.tkwin, sw->sw.width_in_pixels,
	    		       sw->sw.height_in_pixels);
	}
#endif

	sw->flags |= SHEET_REDRAW_ALL;
	SheetDisplay((ClientData)sw);

	if (sw->extension)
	    sw->extension(sw->extensionData, SHEET_JOB_RESIZE, NULL);

	break;
    }

    case DestroyNotify:

        Tcl_DeleteCommand(sw->interp, Tk_PathName(sw->sw.tkwin));
        sw->sw.tkwin = NULL;
	sw->flags |= SHEET_DESTROYED;
        if (sw->flags & SHEET_REDRAW_PENDING)
            Tcl_CancelIdleCall(SheetDisplay,
			       (ClientData)sw);
        Tcl_EventuallyFree((ClientData)sw, SheetDestroy);

	break;
/*
    default:
	printf("Unrequested event %d\n", eventPtr->type);
*/
    }
}


/*
 * Destroy
 */
static void SheetDestroy(char *clientData) {
    tkSheet *sw = (tkSheet *)clientData;

    if (sw->extension)
	sw->extension(sw->extensionData, SHEET_JOB_DESTROY, NULL);

    Tk_FreeOptions(sw->configSpecs, (char *)sw, sw->sw.display, 0);

    Tcl_CancelIdleCall(SheetDisplay, (ClientData)sw);

    sheet_destroy(&sw->sw);
    xfree(sw);
}
