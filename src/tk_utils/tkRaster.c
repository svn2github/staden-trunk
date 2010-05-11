/*
 * tkRaster.c --
 *
 *	This module implements "raster" widgets.  A "raster" is
 *      a rectangular array of pixels visible through a window.
 *      Commands are provided for drawing simple shapes on the
 *      raster. A raster may also have additional named pixmaps that
 *      may be used to move rectangular chunks of pixels to/from
 *      the displayed pixmap through the use of bitblt operations.
 *
 * Acknowlegment:
 *
 *   	This stuff has been mostly copied from file tkSquare.c of
 *      tk's distribution.
 *
 * Disclaimer:
 *
 *      Use this at your own risk.
 */
/* HACK - put somewhere else */
#define ROUND(x)   ((x) < 0.0 ? ceil((x) - 0.5) : floor((x) + 0.5))

#define DEF_RASTER_HEIGHT "300"
#define DEF_RASTER_WIDTH "400"
#define DEF_RASTER_BORDER_WIDTH "2"
#define DEF_RASTER_GEOMETRY (char *)NULL
#define DEF_RASTER_X_SCROLL_CMD (char*)NULL
#define DEF_RASTER_Y_SCROLL_CMD (char*)NULL
#define DEF_RASTER_STIPPLE None
#define DEF_RASTER_SCROLL_INCREMENT "10"

#define RASTER_DRAW_PIXMAP 1
#define RASTER_DRAW_BORDER 2
#define RASTER_DRAW_ALL    3


#include <stdio.h>
#include <memory.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <tcl.h>
#include <tk.h>
#include <X11/Xlib.h>

#define HIGH INT_MAX
#define LOW INT_MIN


/* For the prototypes for our implementations of X functions */
#ifdef _WIN32
#include "tkWinX.h"
#endif


#include "tkRaster.h"
#include "tkRasterBuiltIn.h"
#include "xalloc.h"
#include "tcl_utils.h"

/*
 * A Drawing Environment is a set of arguments that affect the
 * drawing of primitives (color, line thickness, fill pattern etc). It
 * is implemented as a  XGCValues structure plus a bitmask. A raster may
 * have many of these loaded at the same time
 */
/*
typedef struct {
   XGCValues gcValues;
   unsigned long valMask;
   XColor * fgColor, * bgColor;
   int index;
} DrawEnvironment;
*/
/*
 * A data structure of the following type is kept for each Raster
 * widget managed by this file:
 */

typedef struct Raster_t {
    Tk_Window tkwin;		/* Window that embodies the Raster.  NULL
				 * means window has been deleted but
				 * widget record hasn't been cleaned up yet. */
    Display *display;		/* X's token for the window's display. */
    Tcl_Interp *interp;		/* Interpreter associated with widget. */
    int x, y;			/* Position of viewable pixmap's upper-left
				 * corner within the widget's area */

    int new_x, new_y;
    int initialised;		/* Startup: has configure finished yet */
    char * geometry;		/* This tells us the geometry that should
				   be requested for the window. If null,
				   window should be big enough to display
				   raster */

    char *xScrollCmd;		/* Command prefix for communicating with
				 * horizontal scrollbar.  NULL means no
				 * horizontal scrollbar.  Malloc'ed*/

    char *yScrollCmd;		/* Command prefix for communicating with
				 * vertical scrollbar.  NULL means no
				 * vertical scrollbar.  Malloc'ed*/

    int scrollIncrement;	/* Scroll increment units */
    int winWidth;		/* Width of window. */
    int winHeight;		/* Height of window. */

    int width;			/* Width of raster. */
    int height;                 /* Height of raster. */

    Cursor cursor;              /* Current cursor for window, or None. */

    double wx_start;            /* extents of scrolling region */
    double wy_start;            /* used in x and y scrolling */
    double wx_end;
    double wy_end;

    /*
     * GC related options
     */

    GC drawGC;			/* Graphics context for drawing */
    Tcl_HashTable drawEnvTable; /* Set of Drawing Environments currently
				   defined for this Raster */
    int drawEnvCount;		/* Number of Drawing Environments created
				   for this raster so far */
    DrawEnvironment* currentDrawEnv;
                                /* Drawing Environment currently loaded
				   in drawGC */

    /*
     * Information used when displaying widget:
     */

    int borderWidth;		/* Width of 3-D border around whole widget. */
    Tk_3DBorder bgBorder;	/* Used for drawing background. */
    XColor *fgColor;		/* For drawing visible pixmap. */
    int relief;			/* Indicates whether window as a whole is
				 * raised, sunken, or flat. */
    GC copyGC;        	        /* Graphics context for copying to screen */
    Pixmap pm;			/* Pixmap for storing the raster */
    int doubleBuffer;		/* Non-zero means double-buffer redisplay
				 * with pixmap;  zero means draw straight
				 * onto the display. */
    int updatePending;		/* Non-zero means a call to DisplayRaster
				 * has already been scheduled. */
    int x0, y0, x1, y1;		/* coordinates of rectangle of the
				   pixmap that has to be redisplayed */
    int doclip;                 /* Tells the display routine whether
				   updating the window can be done using
				   a clip mask or not. If so, the 'clip'
				   rectangle is used. */
    int cx0,cy0,cx1,cy1;	/* Clip Rectangle */

    /*
     *  Values used for World-To-Raster-And-Back transformations:
     */

    double ax, bx;     		/* rx = (int) ax * (wx - bx) */
    double ay, by;		/* ry = (int) ay * (wy - by) */
    double wx0, wy0, wx1, wy1;  /* Coords of upperleft corner and bottom right
				   corners of the "world" */

    /*
     *  Area modified when drawing a primitive
     */

    int px0, py0, px1, py1;
    int knownarea;

    /*
     *  Private storage needed for each primitive
     */

    ClientData *primitiveData;

    /*void (*plot_func)(struct Raster_t *raster, int x0, int y0, int x1, int y1); */
    void (*plot_func)(Tk_Raster *raster, char *raster_win, int job,
		      int x0, int y0, int x1, int y1);
} Raster;



/*
 * A Tcl Hash table is used to parse the built-in and user-specified
 * raster primitives. An entry value in this hash table is a pointer
 * to a RasterImplement structure.
 */

typedef struct {
   char* primitive;
   int primitiveno;
   RasterPrimDrawProc * draw;
   RasterPrimInitProc * init;
   RasterPrimFreeProc * free;
} RasterImplement;

static Tcl_HashTable PrimitiveTable;  /* Hash table for storing primitive
					 implementations */
static int PrimitiveCount = 0;	      /* Number of primitives presently
					 supported */

/*
 * Forward declarations for procedures defined later in this file:
 */

static int  ConfigureRaster _ANSI_ARGS_((Tcl_Interp *, Raster *,
					int argc, char **argv, int flags));
static void DestroyRaster _ANSI_ARGS_((char *clientData));
static void DisplayRaster _ANSI_ARGS_((ClientData clientData));
static void RasterEventProc _ANSI_ARGS_((ClientData clientData,
					 XEvent *eventPtr));
static int  RasterWidgetCmd _ANSI_ARGS_((ClientData clientData, Tcl_Interp *,
					int argc, char **argv));
static int  RasterDraw _ANSI_ARGS_((Tcl_Interp*, Raster*, RasterImplement*,
				   int argc, char ** argv));
static void arrangeDisplay _ANSI_ARGS_((Raster*, int x0, int y0,
					int x1, int y1, int type));
static void arrangeExpose _ANSI_ARGS_((Raster*, int x0, int y0,
				       int x1, int y1));
static int  myOptionParse _ANSI_ARGS_((ClientData, Tcl_Interp *,
				       Tk_Window, char*, char*, int));
static char * myOptionPrint _ANSI_ARGS_((ClientData, Tk_Window, char*,
					 int, Tcl_FreeProc ** ));

static int ConfigDrawEnv _ANSI_ARGS_((Tcl_Interp*, Raster*, DrawEnvironment*,
				      int argc, char * argv []));

static int CreateDrawEnv _ANSI_ARGS_((Tcl_Interp *, Raster*,
				      int argc,  char* argv []));
static void DestroyDrawEnv _ANSI_ARGS_((Raster*, DrawEnvironment*));
static void RasterScrollX(Raster *raster, int x, int old);
static void RasterScrollY(Raster *raster, int y, int old);
static void RasterClear(Raster *RasterPtr);

/*
 * Not implemented (by Tk) configuration options used in this widgets
 * are supported via custom options (see Tk_ConfigureWidget).
 * The parsing and printing of these options are done through procs
 * myOptionParse and myOptionPrint respectively. Those functions
 * require as a ClientData argument a table where each entry is
 * composed of a string and an associated integer value (AttrNameValue
 * below).
 */

typedef struct {
   char * optionname;
   int optionvalue;
} AttrNameValue;

AttrNameValue FillStyleTable [] = {
   { "solid", FillSolid },
   { "opaquestippled", FillOpaqueStippled },
   { "stippled", FillStippled },
   { (char*) NULL, 0 }
};

AttrNameValue LineStyleTable [] = {
   { "solid", LineSolid },
   { "doubledash", LineDoubleDash },
   { "onoffdash", LineOnOffDash },
   { (char*) NULL, 0 }
};

AttrNameValue FunctionTable [] = {
   { "clear", GXclear },
   { "and", GXand },
   { "andreverse", GXandReverse },
   { "copy", GXcopy },
   { "andinverted", GXandInverted },
   { "noop", GXnoop },
   { "xor", GXxor },
   { "or", GXor },
   { "nor", GXnor },
   { "equiv", GXequiv },
   { "invert", GXinvert },
   { "orreverse", GXorReverse },
   { "copyinverted", GXcopyInverted },
   { "orinverted", GXorInverted },
   { "nand", GXnand },
   { "set", GXset },
   { (char*) NULL, 0 }
};

Tk_CustomOption FillStyleOption = {
   myOptionParse,
   myOptionPrint,
   (ClientData) FillStyleTable
};

Tk_CustomOption LineStyleOption = {
   myOptionParse,
   myOptionPrint,
   (ClientData) LineStyleTable
};

Tk_CustomOption FunctionOption = {
   myOptionParse,
   myOptionPrint,
   (ClientData) FunctionTable
};


/*
 * Information used for argv parsing: Some config items correspond to
 * setting values in the Graphics Context structure (drawGC of Raster).
 * Depending on the type of feature being drawn (Line, Rect, etc),
 * some GC values may be reset at drawing time and reset to the
 * old values after the drawing is finished by passing appropriate
 * config options. See the code for rasterDraw below.
 */

static Tk_ConfigSpec configSpecs[] = {
    {TK_CONFIG_BORDER, "-background", "background",  "Background",
	"#d9d9d9", Tk_Offset(Raster, bgBorder),  TK_CONFIG_COLOR_ONLY},
    {TK_CONFIG_BORDER, "-background", "background", "Background",
	"white", Tk_Offset(Raster, bgBorder), TK_CONFIG_MONO_ONLY},
    {TK_CONFIG_COLOR, "-foreground", "foreground", "Foreground",
	"brown", Tk_Offset(Raster, fgColor), TK_CONFIG_COLOR_ONLY},
    {TK_CONFIG_COLOR, "-foreground", "foreground", "Foreground",
	"black", Tk_Offset(Raster, fgColor), TK_CONFIG_MONO_ONLY},
    {TK_CONFIG_SYNONYM, "-bg", "background", (char *) NULL,
	(char *) NULL, 0, 0},
    {TK_CONFIG_SYNONYM, "-fg", "foreground", (char *) NULL,
	(char *) NULL, 0, 0},
    {TK_CONFIG_SYNONYM, "-bd", "borderWidth", (char *) NULL,
	(char *) NULL, 0, 0},
    {TK_CONFIG_PIXELS, "-borderwidth", "borderWidth", "BorderWidth",
	DEF_RASTER_BORDER_WIDTH, Tk_Offset(Raster, borderWidth), 0},
    {TK_CONFIG_ACTIVE_CURSOR, "-cursor", "cursor", "Cursor",
        (char*) NULL, Tk_Offset(Raster, cursor), TK_CONFIG_NULL_OK},
    {TK_CONFIG_INT, "-dbl", "doubleBuffer", "DoubleBuffer",
	"0", Tk_Offset(Raster, doubleBuffer), 0},
    {TK_CONFIG_INT, "-x", "x", "X",
	"0", Tk_Offset(Raster, x), 0},
    {TK_CONFIG_INT, "-y", "y", "Y",
	"0", Tk_Offset(Raster, y), 0},
    {TK_CONFIG_STRING, "-geometry", "geometry", "Geometry",
	DEF_RASTER_GEOMETRY, Tk_Offset(Raster, geometry), TK_CONFIG_NULL_OK},
    {TK_CONFIG_RELIEF, "-relief", "relief", "Relief",
	"raised", Tk_Offset(Raster, relief), 0},
    {TK_CONFIG_PIXELS, "-width", "width", "Width",
	DEF_RASTER_WIDTH, Tk_Offset(Raster, width), 0},
    {TK_CONFIG_PIXELS, "-height", "height", "Height",
	DEF_RASTER_HEIGHT, Tk_Offset(Raster, height), 0},
    {TK_CONFIG_STRING, "-xscrollcommand", "xScrollCommand", "ScrollCommand",
	DEF_RASTER_X_SCROLL_CMD, Tk_Offset(Raster, xScrollCmd),
	TK_CONFIG_NULL_OK},
    {TK_CONFIG_STRING, "-yscrollcommand", "yScrollCommand", "ScrollCommand",
	DEF_RASTER_Y_SCROLL_CMD, Tk_Offset(Raster, yScrollCmd),
	TK_CONFIG_NULL_OK},
    {TK_CONFIG_PIXELS, "-scrollincrement", "scrollIncrement","ScrollIncrement",
	DEF_RASTER_SCROLL_INCREMENT, Tk_Offset(Raster, scrollIncrement), 0},
    {TK_CONFIG_END, (char *) NULL, (char *) NULL, (char *) NULL,
	(char *) NULL, 0, 0}
};



/*----------------------------------------------------------------------
 *  RasterInit --
 *
 *  	Initializes the Raster Widget module. This involves basically
 *  	setting up the built in raster primitives.
 *
 *  Results:
 *  	A standard Tcl result indicating the success of the initialization
 *
 *  Side effects:
 *	A new hash table is created.
 *----------------------------------------------------------------------
 */
int Raster_Init (Tcl_Interp *interp)
{
   Tcl_InitHashTable (&PrimitiveTable, TCL_STRING_KEYS);
   Tcl_CreateCommand(interp, "raster", RasterCmd,
		     (ClientData)Tk_MainWindow(interp),
		     NULL);
  return RasterBuiltInInit (interp);
}


/*
 *--------------------------------------------------------------
 *
 * RasterCmd --
 *
 *	This procedure is invoked to process the "Raster" Tcl
 *	command.  It creates a new "Raster" widget.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	A new widget is created and configured.
 *
 *--------------------------------------------------------------
 */

int
RasterCmd(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)
{
    Tk_Window tkwin;
    Raster *RasterPtr;
    ClientData newenv;
    ClientData * dataPtr;
    Tcl_HashSearch search;
    Tcl_HashEntry* entryPtr;
    RasterImplement* primitivePtr;
    int result;

    if (argc < 2) {
	Tcl_AppendResult(interp, "wrong # args:  should be \"",
		argv[0], " pathName ?options?\"", (char *) NULL);
	return TCL_ERROR;
    }

    tkwin = Tk_CreateWindowFromPath(interp, Tk_MainWindow(interp), argv[1],
				    (char *) NULL);
    if (tkwin == NULL) {
	return TCL_ERROR;
    }
    Tk_SetClass(tkwin, "Raster");

    /*
     * Allocate and initialize the widget record.
     */

    RasterPtr = (Raster *) ckalloc(sizeof(Raster));
    RasterPtr->tkwin = tkwin;
    RasterPtr->display = Tk_Display(tkwin);
    RasterPtr->interp = interp;
    RasterPtr->x = 0;
    RasterPtr->y = 0;
    RasterPtr->new_x = 0;
    RasterPtr->new_y = 0;
    RasterPtr->width = None;
    RasterPtr->height = None;
    RasterPtr->geometry = NULL;
    RasterPtr->winWidth = None;
    RasterPtr->winHeight = None;
    RasterPtr->borderWidth = 0;
    RasterPtr->bgBorder = NULL;
    RasterPtr->fgColor = NULL;
    RasterPtr->relief = TK_RELIEF_FLAT;
    RasterPtr->cursor = None;
    RasterPtr->copyGC = None;
    RasterPtr->drawGC = None;
    RasterPtr->pm = None;
    RasterPtr->doubleBuffer = 0;
    RasterPtr->updatePending = 0;
    RasterPtr->ax = RasterPtr->ay = 1.0;
    RasterPtr->bx = RasterPtr->by = 0.0;
    RasterPtr->wx0 = RasterPtr->wx1 =
    RasterPtr->wy0 = RasterPtr->wy1 = 0.0;
    RasterPtr->xScrollCmd = NULL;
    RasterPtr->yScrollCmd = NULL;
    RasterPtr->scrollIncrement = 1;
    RasterPtr->drawEnvCount = 0;
    RasterPtr->currentDrawEnv = (DrawEnvironment*) NULL;
    RasterPtr->px0 = HIGH;
    RasterPtr->py0 = HIGH;
    RasterPtr->px1 = LOW;
    RasterPtr->py1 = LOW;
    RasterPtr->knownarea = 0;
    RasterPtr->initialised = 0;
    RasterPtr->wx_start = DBL_MAX/2;
    RasterPtr->wy_start = DBL_MAX/2;
    RasterPtr->wx_end = -DBL_MAX/2;
    RasterPtr->wy_end = -DBL_MAX/2;
    RasterPtr->plot_func = NULL;
    RasterPtr->x0 = HIGH;
    RasterPtr->y0 = HIGH;
    RasterPtr->x1 = LOW;
    RasterPtr->y1 = LOW;

    Tcl_InitHashTable (&(RasterPtr->drawEnvTable), TCL_ONE_WORD_KEYS);

    /* A new drawing environment numbered "0" must be created
       and installed with the following command */
    if (CreateDrawEnv (interp, RasterPtr, 0, (char**) NULL) != TCL_OK ||
	DrawEnvIndex (interp, RasterPtr, 0, &newenv) != TCL_OK) {
       return TCL_ERROR;
    }
    RasterPtr->currentDrawEnv = (DrawEnvironment*) newenv;

    Tk_CreateEventHandler(RasterPtr->tkwin, ExposureMask|StructureNotifyMask,
	    RasterEventProc, (ClientData) RasterPtr);
    Tcl_CreateCommand(interp, Tk_PathName(RasterPtr->tkwin), RasterWidgetCmd,
	    (ClientData) RasterPtr, NULL);
    if (ConfigureRaster(interp, RasterPtr, argc-2, argv+2, 0) != TCL_OK ||
	SetDrawEnv (interp, RasterPtr, newenv) != TCL_OK) {
	Tk_DestroyWindow(RasterPtr->tkwin);
	return TCL_ERROR;
    }

    /* Initialize primitives */
    RasterPtr->primitiveData = (ClientData*) malloc (sizeof(ClientData) *
						     PrimitiveCount);
    result = TCL_OK;
    for (entryPtr = Tcl_FirstHashEntry (&PrimitiveTable, &search);
	 entryPtr != NULL;
	 entryPtr = Tcl_NextHashEntry (&search)) {
       primitivePtr = (RasterImplement*) Tcl_GetHashValue (entryPtr);
       dataPtr = &(RasterPtr->primitiveData[primitivePtr->primitiveno]);
       if (primitivePtr->init != NULL) {
	  if ((*(primitivePtr->init)) (interp, RasterPtr, dataPtr) != TCL_OK)
	     result = TCL_ERROR;
       } else {
	   *dataPtr = (ClientData)0;
       }
    }
    if (result != TCL_OK) {
       Tk_DestroyWindow(RasterPtr->tkwin);
       return TCL_ERROR;
    }

    Tcl_SetResult(interp, Tk_PathName(RasterPtr->tkwin), TCL_STATIC);
    return TCL_OK;
}

/*----------------------------------------------------------------------
 * RasterAddPrimitive --
 *
 *  	Implements a new geometric primitive to be used in connection with
 *  	a raster
 *
 * Results:
 *
 *  	A standard tcl result
 *
 * Side effects:
 *
 *	Another entry in the PrimitiveTable is created
 *
 *----------------------------------------------------------------------
 */
int RasterAddPrimitive (interp, primitivename, drawproc, initproc, freeproc)
     Tcl_Interp * interp;
     char * primitivename;
     RasterPrimDrawProc * drawproc;
     RasterPrimInitProc * initproc;
     RasterPrimFreeProc * freeproc;
{
   int new;
   Tcl_HashEntry * entryPtr;
   RasterImplement *implemPtr;
   entryPtr = Tcl_CreateHashEntry (&PrimitiveTable, primitivename, &new);
   if (!new) {
      Tcl_AppendResult (interp, primitivename,
		       " could not be installed", (char*)NULL);
      return TCL_ERROR;
   }
   implemPtr = (RasterImplement*) malloc (sizeof (RasterImplement));
   Tcl_SetHashValue (entryPtr, implemPtr);
   implemPtr->primitive = primitivename;
   implemPtr->primitiveno = PrimitiveCount++;
   implemPtr->init = initproc;
   implemPtr->draw = drawproc;
   implemPtr->free = freeproc;
   return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * RasterWidgetCmd --
 *
 *	This procedure is invoked to process the Tcl command
 *	that corresponds to a widget managed by this module.
 *	See the user documentation for details on what it does.
 *
 * Results:
 *	A standard Tcl result.
 *
 * Side effects:
 *	See the user documentation.
 *
 *--------------------------------------------------------------
 */

static int
RasterWidgetCmd(clientData, interp, argc, argv)
    ClientData clientData;		/* Information about Raster widget. */
    Tcl_Interp *interp;			/* Current interpreter. */
    int argc;				/* Number of arguments. */
    char **argv;			/* Argument strings. */
{
    Raster *RasterPtr = (Raster *) clientData;
    int result = TCL_OK;
    int length;
    char c;
    int index;
    int rx, ry;
    double x0, x1, y0, y1;
    Tcl_HashEntry *entryPtr;
    RasterImplement * implemPtr;
    DrawEnvironment * drawEnvPtr;

    if (argc < 2) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
		argv[0], " option ?arg arg ...?\"", (char *) NULL);
	return TCL_ERROR;
    }

    Tcl_Preserve((ClientData) RasterPtr);
    c = argv[1][0];
    length = strlen(argv[1]);

    if ((c == 'c') && (strncmp(argv[1], "cget", length) == 0)) {
	/* CGET COMMAND */
	if (argc != 3) {
            Tcl_AppendResult(interp, "wrong # args: should be \"",
                    argv[0], " cget option\"",
                    (char *) NULL);
            goto error;
        }
        result = Tk_ConfigureValue(interp, RasterPtr->tkwin, configSpecs,
                (char *) RasterPtr, argv[2], 0);
    }
    else if ((c == 'c') && (strncmp(argv[1], "configure", length) == 0)) {
        /* CONFIGURE COMMAND */
	if (argc == 2) {
	    result = Tk_ConfigureInfo(interp, RasterPtr->tkwin, configSpecs,
		    (char *) RasterPtr, (char *) NULL, 0);
	} else if (argc == 3) {
	    result = Tk_ConfigureInfo(interp, RasterPtr->tkwin, configSpecs,
		    (char *) RasterPtr, argv[2], 0);
	} else {
	    result = ConfigureRaster(interp, RasterPtr, argc-2, argv+2,
		    TK_CONFIG_ARGV_ONLY);
	}
	if (result != TCL_OK) goto error;
	arrangeDisplay (RasterPtr, LOW, LOW, HIGH, HIGH, RASTER_DRAW_ALL);
    }
    else if ((c == 'c') && (strncmp (argv [1], "clear", length) == 0)) {
       /* CLEAR COMMAND */
       if (argc != 2) {
	  Tcl_AppendResult (interp, "wrong # args: should be \"",
			    argv[0], " clear", (char *) NULL);
	  goto error;
       }
       XFillRectangle (RasterPtr->display, RasterPtr->pm, RasterPtr->copyGC,
		       0, 0, RasterPtr->width, RasterPtr->height);
       arrangeDisplay (RasterPtr, LOW, LOW, HIGH, HIGH, RASTER_DRAW_PIXMAP);
    }
    else if ((c == 'x') &&  (strncmp(argv[1], "xview", length) == 0)) {
       /* XVIEW COMMAND */
	int count, type;
	double fraction;
	double world_length = RasterPtr->wx_end - RasterPtr->wx_start;
	int win_width = Tk_Width(RasterPtr->tkwin) -
	    (2 * RasterPtr->borderWidth);

	if (argc == 2) {
/*
	    PrintScrollFractions(canvasPtr->xOrigin + canvasPtr->inset,
		    canvasPtr->xOrigin + Tk_Width(canvasPtr->tkwin)
		    - canvasPtr->inset, canvasPtr->scrollX1,
		    canvasPtr->scrollX2, interp->result);
*/
	} else {
	    type = Tk_GetScrollInfo(interp, argc, argv, &fraction, &count);
#ifdef DEBUG
	    printf("XVIEW\n");
	    printf("FRACTION %f\n", fraction);
#endif

	    switch (type) {
		case TK_SCROLL_ERROR:
		    goto error;
		case TK_SCROLL_MOVETO: {
		    double wdx;

		    //		    if (fraction <= 0)
		    //			fraction = 0;

		    wdx = RasterPtr->wx1 - RasterPtr->wx0;

#ifdef DEBUG
		    printf("world_length %f wdx %f %f\n", world_length, wdx,
			   (world_length - wdx) / world_length);
#endif
		    //		    if (fraction >= (world_length - wdx) / world_length)
		    //			fraction = (world_length - wdx) / world_length;
		    WorldToRaster(RasterPtr, (world_length * fraction +
					      RasterPtr->wx_start), 0,
				  &rx, &ry);
		    RasterPtr->new_x = rx + RasterPtr->x;
#ifdef DEBUG
		    printf("MOVETO %f %d rx %d\n", fraction,
			   RasterPtr->new_x, rx);
#endif
		    break;
		}
	    case TK_SCROLL_PAGES: {
		int rx0, ry0, rx1, ry1;
#ifdef DEBUG
		printf("scroll pages \n");
		printf("count %d \n", count);

#endif
		WorldToRaster(RasterPtr, RasterPtr->wx_start,
			      RasterPtr->wy_start, &rx0, &ry0);
		WorldToRaster(RasterPtr, RasterPtr->wx_end,
			      RasterPtr->wy_end, &rx1, &ry1);

		RasterPtr->new_x += count * .9 * win_width;

		if (rx0 - (count*.9* win_width) >= 0) {
		    RasterPtr->new_x += rx0 - (count * .9 * win_width);
		}

		if (rx1-(count*.9*win_width) <= win_width) {
		    RasterPtr->new_x += (rx1 - (count*.9*win_width) - win_width);
		}

		break;
	    }
	    case TK_SCROLL_UNITS: {
		int rx0, ry0, rx1, ry1;

		    WorldToRaster(RasterPtr, RasterPtr->wx_start,
				  RasterPtr->wy_start, &rx0, &ry0);
		    WorldToRaster(RasterPtr, RasterPtr->wx_end,
				  RasterPtr->wy_end, &rx1, &ry1);

		    if (RasterPtr->scrollIncrement > 0) {
			RasterPtr->new_x += count*RasterPtr->scrollIncrement;
		    } else {
			RasterPtr->new_x += count * .1 *
			    (Tk_Width(RasterPtr->tkwin));
		    }
#ifdef DEBUG
		    printf("SCROLL %d x %d new %d rx0 %d rx1 %d width %d\n",
			   count*RasterPtr->scrollIncrement,
			   RasterPtr->x, RasterPtr->new_x, rx0, rx1,
			   win_width);
		    printf("wx0 %f wx1 %f\n", RasterPtr->wx0, RasterPtr->wx1);
#endif
		    if (rx0 - (count*RasterPtr->scrollIncrement) >= 0) {
			RasterPtr->new_x += rx0-(count*RasterPtr->scrollIncrement);
		    }

		    if (rx1-(count*RasterPtr->scrollIncrement) <= win_width) {
			RasterPtr->new_x += (rx1 - (count*RasterPtr->scrollIncrement) - win_width);
		    }

		    break;
		}
	    }
	}

	arrangeDisplay (RasterPtr, LOW, LOW, HIGH, HIGH, RASTER_DRAW_PIXMAP);
    }
    else if ((c == 'y') &&  (strncmp(argv[1], "yview", length) == 0)) {
	/* YVIEW COMMAND */
	int count, type;
	double fraction;
	double world_height = RasterPtr->wy_end - RasterPtr->wy_start;
	int win_height = Tk_Height(RasterPtr->tkwin) -
	    (2 * RasterPtr->borderWidth);

	if (argc == 2) {
	    /*
	    PrintScrollFractions(canvasPtr->yOrigin + canvasPtr->inset,
	    canvasPtr->yOrigin + Tk_Height(canvasPtr->tkwin)
	    - canvasPtr->inset, canvasPtr->scrollY1,
	    canvasPtr->scrollY2, interp->result);
	    */

	} else {
	    type = Tk_GetScrollInfo(interp, argc, argv, &fraction, &count);
	    switch (type) {
	    case TK_SCROLL_ERROR:
		goto error;
	    case TK_SCROLL_MOVETO: {
		double wdy;

		if (fraction <= 0)
		    fraction = 0.0;

		wdy = RasterPtr->wy1 - RasterPtr->wy0;

#ifdef DEBUG
		printf("world_height %f wdy %f %f\n", world_height, wdy,
		       (world_height - wdy) / world_height);
#endif
		if (fraction >= (world_height - wdy) / world_height) {
		    fraction = (world_height - wdy) / world_height;
		}
#ifdef DEBUG
		printf("fraction %f \n", fraction);
#endif
		WorldToRaster(RasterPtr, 0, (world_height * fraction +
					     RasterPtr->wy_start),
			      &rx, &ry);
		RasterPtr->new_y = ry + RasterPtr->y;
#ifdef DEBUG
		printf("MOVETO %f %d ry %d RasterPtr->y %d\n", fraction,
		       RasterPtr->new_y, ry, RasterPtr->y);
#endif
		break;
	    }
	    case TK_SCROLL_PAGES:
		{
		    int rx0, ry0, rx1, ry1;
#ifdef DEBUG
		    printf("scroll pages \n");
#endif
		    WorldToRaster(RasterPtr, RasterPtr->wx_start,
				  RasterPtr->wy_start, &rx0, &ry0);
		    WorldToRaster(RasterPtr, RasterPtr->wx_end,
				  RasterPtr->wy_end, &rx1, &ry1);

		    RasterPtr->new_y += count * .9 * win_height;

		    if (ry0 - (count*.9* win_height) >= 0) {
			RasterPtr->new_y += ry0 - (count * .9 * win_height);
		    }

		    if (ry1-(count*.9*win_height) <= win_height) {
			RasterPtr->new_y += (ry1 - (count*.9*win_height) - win_height);
		    }
		    break;
		}
	    case TK_SCROLL_UNITS: {
		int rx0, ry0, rx1, ry1;

		WorldToRaster(RasterPtr, RasterPtr->wx_start,
			      RasterPtr->wy_start, &rx0, &ry0);
		WorldToRaster(RasterPtr, RasterPtr->wx_end,
			      RasterPtr->wy_end,
			      &rx1, &ry1);

		if (RasterPtr->scrollIncrement > 0) {
		    RasterPtr->new_y += count*RasterPtr->scrollIncrement;
		} else {
		    RasterPtr->new_y += count * .1 *
			(Tk_Height(RasterPtr->tkwin));
		}
		printf("SCROLL %d y %d new %d ry0 %d ry1 %d height %d\n",
		       count*RasterPtr->scrollIncrement,
		       RasterPtr->y, RasterPtr->new_y, ry0, ry1,
		       win_height);
		printf("wy0 %f wy1 %f\n", RasterPtr->wy0, RasterPtr->wy1);
		if (ry0 - (count*RasterPtr->scrollIncrement) >= 0) {
		    RasterPtr->new_y += ry0-(count*RasterPtr->scrollIncrement);
		    printf("NEW_Y1 %d\n", RasterPtr->new_y);
		}

		if (ry1-(count*RasterPtr->scrollIncrement) <= win_height) {
		    RasterPtr->new_y += (ry1 - (count*RasterPtr->scrollIncrement) - win_height);
		    printf("NEW_Y2 %d\n", RasterPtr->new_y);
		}

		break;
	    }
	    }
	}

	arrangeDisplay (RasterPtr, LOW, LOW, HIGH, HIGH, RASTER_DRAW_PIXMAP);
    }
    else if ((c == 'e') && (strncmp (argv [1], "envset", length) == 0)) {
       /* ENVSET COMMAND */
       if (argc > 2) {
	  if (Tcl_GetInt(RasterPtr->interp, argv[2], &index) != TCL_OK) {
	     goto error;
	  }
	  if (DrawEnvIndex (interp, RasterPtr, index,
			    (ClientData*)&drawEnvPtr) != TCL_OK){
	     goto error;
	  }
	  result = SetDrawEnv (interp, RasterPtr, drawEnvPtr);
       }
       else {
	   vTcl_SetResult(interp, "%d", RasterPtr->currentDrawEnv->index);
       }
    }
    else if ((c == 'e') && (strncmp (argv [1], "envcreate", length) == 0)) {
       /* ENVCREATE COMMAND */
       result = CreateDrawEnv (interp, RasterPtr, argc-2, argv+2);
    }
    else if ((c == 'e') && (strncmp (argv [1], "envconfigure", length) == 0)) {
       /* ENVCONFIG COMMAND */
       if (argc > 2 && isdigit(argv [2][0])) {
	  if (Tcl_GetInt(RasterPtr->interp, argv[2], &index) != TCL_OK) {
	     goto error;
	  }
	  argc -= 3;
	  argv += 3;
       }
       else {
	  index = 0;
	  argc -= 2;
	  argv += 2;
       }
       if (DrawEnvIndex (interp, RasterPtr, index,
			 (ClientData*) &drawEnvPtr) != TCL_OK) {
	  goto error;
       }
       if (argc < 2) {
	  result = ConfigInfoDrawEnv (interp, RasterPtr, drawEnvPtr,
				      argc, argv);
       } else {
	  result = ConfigDrawEnv (interp, RasterPtr, drawEnvPtr, argc, argv);
       }
    }
    else if ((c == 'e') && (strncmp (argv [1], "envdelete", length) == 0)) {
       /* ENVDELETE COMMAND */
       if (argc > 2) {
	  if (Tcl_GetInt(RasterPtr->interp, argv[2], &index) != TCL_OK) {
	     goto error;
	  }
	  if (DrawEnvIndex (interp, RasterPtr, index,
			    (ClientData*)&drawEnvPtr) != TCL_OK){
	     goto error;
	  }
       }
       else {
	  drawEnvPtr = RasterPtr->currentDrawEnv;
	  index = drawEnvPtr->index;
       }
       if (index == 0) {
	  Tcl_AppendResult (interp, "Cannot delete drawing environment 0",
			    (char*) NULL);
	  goto error;
       }
       DestroyDrawEnv (RasterPtr, drawEnvPtr);
       if (drawEnvPtr == RasterPtr->currentDrawEnv) {
	  DrawEnvIndex (interp, RasterPtr, 0, (ClientData*) &drawEnvPtr);
	  SetDrawEnv (interp, RasterPtr, drawEnvPtr);
       }
    }
    else if ((c == 'w') && (strncmp (argv [1], "world", length) == 0)) {
       /* WORLD COMMAND */
       if (argc == 2) {
	  /*return world coordinates */
	  vTcl_SetResult(interp, "%.9g %.9g %.9g %.9g",
			 RasterPtr->wx0, RasterPtr->wy0, RasterPtr->wx1,
			 RasterPtr->wy1);
       }
       else {
	  /* set world coordinates */
	  if (argc < 6) {
	     Tcl_AppendResult (interp, "wrong # args: ",
			       " expected 4 world coordinates", (char*)NULL);
	     goto error;
	  }
	  if (Tcl_GetDouble (interp, argv [2], &x0) != TCL_OK ||
	      Tcl_GetDouble (interp, argv [3], &y0) != TCL_OK ||
	      Tcl_GetDouble (interp, argv [4], &x1) != TCL_OK ||
	      Tcl_GetDouble (interp, argv [5], &y1) != TCL_OK) {
	     goto error;
	  }
	  if (x1 == x0 || y1 == y0) {
	     Tcl_AppendResult (interp, "coordinates must define a rectangle",
			       (char*) NULL);
	     goto error;
	  }
#ifdef DEBUG
	  printf("*************************WORLD*******************\n");
#endif
	  SetRasterCoords (RasterPtr, x0, y0, x1, y1);
       }
    }
    else if ((c == 'w') && (strncmp (argv [1], "world_scroll", length) == 0)) {
       /* WORLD COMMAND */
       if (argc == 2) {
	  /*return world coordinates */
	  vTcl_SetResult(interp, "%.9g %.9g %.9g %.9g",
			 RasterPtr->wx0, RasterPtr->wy0, RasterPtr->wx1,
			 RasterPtr->wy1);
       }
       else {
	  /* set world coordinates */
	  if (argc < 6) {
	     Tcl_AppendResult (interp, "wrong # args: ",
			       " expected 4 world coordinates", (char*)NULL);
	     goto error;
	  }
	  if (Tcl_GetDouble (interp, argv [2], &x0) != TCL_OK ||
	      Tcl_GetDouble (interp, argv [3], &y0) != TCL_OK ||
	      Tcl_GetDouble (interp, argv [4], &x1) != TCL_OK ||
	      Tcl_GetDouble (interp, argv [5], &y1) != TCL_OK) {
	     goto error;
	  }
	  if (x1 == x0 || y1 == y0) {
	     Tcl_AppendResult (interp, "coordinates must define a rectangle",
			       (char*) NULL);
	     goto error;
	  }
#ifdef DEBUG
	  printf("*************************WORLD*******************\n");
#endif
	  RasterSetWorldScroll (RasterPtr, x0, y0, x1, y1);
       }
    }
    else if ((c == 'w') && (strncmp (argv [1], "world_size", length) == 0)) {
	/* WORLD_SIZE COMMAND */
	/* return the total world coords of raster */

	vTcl_SetResult(interp, "%.9g %.9g %.9g %.9g",
		       RasterPtr->wx_start, RasterPtr->wy_start,
		       RasterPtr->wx_end, RasterPtr->wy_end);

    }
    else if ((c == 't') && (strncmp (argv [1], "toworld", length) == 0)) {
       /* TOWORLD COMMAND */
       if (argc < 4) {
	  Tcl_AppendResult (interp, "wrong # args: ",
			    " expected 2 raster pixmap coordinates",
			    (char*)NULL);
	  goto error;
       }
       if (Tcl_GetInt (interp, argv [2], &rx) != TCL_OK ||
	   Tcl_GetInt (interp, argv [3], &ry) != TCL_OK) {
	  goto error;
       }
/*       RasterToWorld (RasterPtr, rx - RasterPtr->x, ry - RasterPtr->y,
		      &x0, &y0);
*/
       RasterToWorld (RasterPtr, rx, ry, &x0, &y0);
#ifdef DEBUG
       printf("TOWORLD %d %d %f %f\n", rx, ry, x0, y0);
#endif
       vTcl_SetResult(interp, "%.9g %.9g", x0, y0);
    }
    else if ((c == 't') && (strncmp (argv [1], "topixmap", length) == 0)) {
       /* TOPIXMAP COMMAND */
       if (argc < 4) {
	  Tcl_AppendResult (interp, "wrong # args: ",
			    " expected 2 world coordinates",
			    (char*)NULL);
	  goto error;
       }
       if (Tcl_GetDouble (interp, argv [2], &x0) != TCL_OK ||
	   Tcl_GetDouble (interp, argv [3], &y0) != TCL_OK) {
	  goto error;
       }
       WorldToRaster (RasterPtr, x0, y0, &rx, &ry);
#ifdef DEBUG
       printf("TOPIXMAP %f %f %d %d\n", x0, y0, rx, ry);
#endif
       vTcl_SetResult(interp, "%d %d", rx, ry);
    }
    else if ((entryPtr = Tcl_FindHashEntry (&PrimitiveTable, argv[1]))
	       != NULL) {
        implemPtr = (RasterImplement*) Tcl_GetHashValue (entryPtr);
        result = RasterDraw (interp, RasterPtr, implemPtr, argc, argv);
    }
    else if ((c == 'x') && (strncmp (argv [1], "xmag", length) == 0)) {
       /* XMAG COMMAND */
	int value;
	double middle, bases;
	double wx0, wx1;
	int min_bases; /* minimum number of bases allowed */
	double base_length;
	int scale_min, scale_max;
	double scroll_length;

	scale_min = get_default_int(interp, tk_utils_defs,
				    "RASTER.SCALEX.MIN");
	scale_max = get_default_int(interp, tk_utils_defs,
				    "RASTER.SCALEX.MAX");
	min_bases = get_default_int(interp, tk_utils_defs,
				    "RASTER.SCALEX.MIN_BASES");

	if (argc < 3) {
	    Tcl_AppendResult (interp, "wrong # args: ",
			      " expected 1 scale coord",
			      (char*)NULL);
	    goto error;
	}
	if (Tcl_GetInt (interp, argv [2], &value) != TCL_OK)
	    goto error;

	if (RasterPtr->wx_end == -DBL_MAX/2 ||
	    RasterPtr->wx_start == DBL_MAX/2) {
	    Tcl_Release((ClientData) RasterPtr);
	    return result;
	}

	wx0 = RasterPtr->wx0;
	wx1 = RasterPtr->wx1;
	if (RasterPtr->wx0 < RasterPtr->wx_start) {
	    scroll_length = (double)(RasterPtr->wx1 - RasterPtr->wx0);
	    wx0 = RasterPtr->wx_start;
	    wx1 = wx0 + scroll_length;
	}

	if (RasterPtr->wx1 > RasterPtr->wx_end) {
	    scroll_length = (double)(RasterPtr->wx1 - RasterPtr->wx0);
	    wx1 = RasterPtr->wx_end;
	    wx0 = wx1 - scroll_length;
	}

	/* REAL HACK - think about this! */
	if (value == scale_min) {
	    wx1 = RasterPtr->wx_end;
	    wx0 = RasterPtr->wx_start;
	}

#ifdef DEBUG
	printf("start %f end %f wx0 %f wx1 %f\n",
	       RasterPtr->wx0, RasterPtr->wx1, wx0, wx1);
#endif

	/* middle = (double)(RasterPtr->wx1 + RasterPtr->wx0) / 2; */
	middle = (wx1 + wx0) / 2;
	base_length = RasterPtr->wx_end - RasterPtr->wx_start;
	/* base_length = wx1 - wx0; */

	bases = RasterSetXMag(interp, base_length, value);

#ifdef DEBUG
	printf("base_length %f\n", base_length);
	printf("bases %f value %d\n", bases, value);
#endif

	wx0 = middle - (bases / 2);
	wx1 = middle + (bases / 2);

#ifdef DEBUG
	printf("value %d middle %f bases %f wx0 %f wx1 %f \n",
	       value, middle, bases, wx0, wx1);
#endif
	SetRasterCoords (RasterPtr, wx0, RasterPtr->wy0, wx1, RasterPtr->wy1);
	if (RasterPtr->plot_func) {
	    /* need to round up wx1 */
	    RasterPtr->plot_func(RasterPtr, Tk_PathName(RasterPtr->tkwin),
				 RASTER_REPLOT_ALL,
				 (int)wx0, (int)RasterPtr->wy0,
				 (int)wx1+1, (int)RasterPtr->wy1);
	}
    }
    else if ((c == 'y') && (strncmp (argv [1], "ymag", length) == 0)) {
       /* YMAG COMMAND */
	double value;
	double middle, scale;
	double wy0, wy1;
	int scale_min = 1;
	double scale_length, scroll_length;

	if (argc < 3) {
	    Tcl_AppendResult (interp, "wrong # args: ",
			      " expected 1 scale coord",
			      (char*)NULL);
	    goto error;
	}
	if (Tcl_GetDouble (interp, argv [2], &value) != TCL_OK)
	    goto error;

	if (RasterPtr->wy_end == -DBL_MAX/2 ||
	    RasterPtr->wy_start == DBL_MAX/2) {
	    Tcl_Release((ClientData) RasterPtr);
	    return result;
	}

	wy0 = RasterPtr->wy0;
	wy1 = RasterPtr->wy1;
	if (RasterPtr->wy0 < RasterPtr->wy_start) {
	    scroll_length = RasterPtr->wy1 - RasterPtr->wy0;
	    wy0 = RasterPtr->wy_start;
	    wy1 = wy0 + (scroll_length);
	}

	if (RasterPtr->wy1 > RasterPtr->wy_end) {
	    scroll_length = RasterPtr->wy1 - RasterPtr->wy0;
	    wy1 = RasterPtr->wy_end;
	    wy0 = wy1 - (scroll_length);
	}

	/* REAL HACK - think about this! */
	if (value == scale_min) {
	    wy1 = RasterPtr->wy_end;
	    wy0 = RasterPtr->wy_start;
	}

#ifdef DEBUG
	printf("wy0 %f wy1 %f wy_start %f wy_end %f \n",
	        RasterPtr->wy0,  RasterPtr->wy1,  RasterPtr->wy_start,
	       RasterPtr->wy_end);

	printf("wy0 %f wy1 %f\n", wy0, wy1);
#endif
	/* middle = (double)(RasterPtr->wy1 + RasterPtr->wy0) / 2; */

	middle = (wy1 + wy0) / 2;

	scale_length = RasterPtr->wy_end - RasterPtr->wy_start;
	scale = RasterSetYMag(interp, scale_length, value);

	wy0 = middle - (scale / 2);
	wy1 = middle + (scale / 2);

#ifdef DEBUG
	printf("value %f middle %f scale %f wy0 %f wy1 %f \n",
	       value, middle, scale, wy0, wy1);
#endif

	SetRasterCoords(RasterPtr, RasterPtr->wx0, wy0, RasterPtr->wx1, wy1);

	if (RasterPtr->plot_func) {
	    RasterPtr->plot_func(RasterPtr,
				 Tk_PathName(RasterPtr->tkwin),
				 RASTER_REPLOT_ALL,
				 (int)RasterPtr->wx0, (int)wy0,
				 (int)RasterPtr->wx1, (int)wy1);
	}
    } else {
        /* OTHERS */
	Tcl_AppendResult(interp, "bad option \"", argv[1],
		"\"", (char *) NULL);
	goto error;
    }

    Tcl_Release((ClientData) RasterPtr);
    return result;

    error:
    Tcl_Release((ClientData) RasterPtr);
    return TCL_ERROR;
}

/*
 *______________________________________________________________________
 *
 * RasterDraw --
 *
 *	This procedure is invoked by RasterWidgetCmd to draw shapes
 *      on the pixmap.
 *
 * Results:
 *	A standart Tcl result
 *
 * Side effects:
 *	An appropriate thing is drawn
 *
 *----------------------------------------------------------------------
 */

static int
RasterDraw (interp, RasterPtr, implemPtr, argc, argv)
     Tcl_Interp * interp;
     Raster* RasterPtr;
     RasterImplement* implemPtr;
     int argc;
     char ** argv;
{
   RasterPtr->px0 = HIGH;
   RasterPtr->py0 = HIGH;
   RasterPtr->px1 = LOW;
   RasterPtr->py1 = LOW;
   RasterPtr->knownarea = 0;

   if ((*(implemPtr->draw)) (interp, RasterPtr,
			     RasterPtr->primitiveData [implemPtr->primitiveno],
			     argc, argv) != TCL_OK) {
      return TCL_ERROR;
   }

   if (RasterPtr->knownarea) {
      /* Primitive informed modified area of pixmap */
      if (RasterPtr->px1 >= 0 && RasterPtr->py1 >= 0) {
	 /* No point in redisplaying an empty area */
	 arrangeDisplay (RasterPtr, RasterPtr->px0, RasterPtr->py0,
			 RasterPtr->px1, RasterPtr->py1,
			 RASTER_DRAW_PIXMAP);
      }
   }
   else {
      /* Redisplay the whole pixmap */
      arrangeDisplay(RasterPtr, 0, 0, RasterPtr->width-1, RasterPtr->height-1,
		     RASTER_DRAW_PIXMAP);
   }

   return TCL_OK;
}

static void
ResizeRaster(Raster *RasterPtr,
	     int wid,
	     int hgt)
{

    Pixmap newpm;
#ifdef DEBUG
    printf("ResizeRaster\n");
#endif
    /* 7/1/99 johnt - changed X Pixmap functions to TK functions for Windows portability */
    newpm = Tk_GetPixmap (RasterPtr->display,
			   RootWindowOfScreen (Tk_Screen (RasterPtr->tkwin)),
			   RasterPtr->width, RasterPtr->height,
			   DefaultDepthOfScreen(Tk_Screen(RasterPtr->tkwin)));

    XFillRectangle (RasterPtr->display, newpm, RasterPtr->copyGC,
		    0, 0, RasterPtr->width, RasterPtr->height);
    XCopyArea (RasterPtr->display, RasterPtr->pm, newpm,
	       RasterPtr->copyGC, 0, 0, wid - RasterPtr->borderWidth,
	       hgt - RasterPtr->borderWidth, 0, 0);
    Tk_FreePixmap (RasterPtr->display, RasterPtr->pm);
    RasterPtr->pm = newpm;

    if (RasterPtr->plot_func) {
	/* HACK - is this the best place to put these? */
	RasterClear(RasterPtr);

	/* set world coords because raster height and width have changed */
	SetRasterCoords(RasterPtr, RasterPtr->wx0, RasterPtr->wy0,
			RasterPtr->wx1, RasterPtr->wy1);

	RasterPtr->plot_func(RasterPtr, Tk_PathName(RasterPtr->tkwin),
			     RASTER_REPLOT_SLIVER,
			     RasterPtr->wx_start, RasterPtr->wy_start,
			     RasterPtr->wx_end, RasterPtr->wy_end);
    }
}

/*
 *----------------------------------------------------------------------
 *
 * ConfigureRaster --
 *
 *	This procedure is called to process an argv/argc list in
 *	conjunction with the Tk option database to configure (or
 *	reconfigure) a Raster widget.
 *
 * Results:
 *	The return value is a standard Tcl result.  If TCL_ERROR is
 *	returned, then interp->result contains an error message.
 *
 * Side effects:
 *	Configuration information, such as colors, border width,
 *	etc. get set for RasterPtr;  old resources get freed,
 *	if there were any.
 *
 *----------------------------------------------------------------------
 */

static int
ConfigureRaster(interp, RasterPtr, argc, argv, flags)
    Tcl_Interp *interp;			/* Used for error reporting. */
    Raster *RasterPtr;			/* Information about widget. */
    int argc;				/* Number of valid entries in argv. */
    char **argv;			/* Arguments. */
    int flags;				/* Flags to pass to
					 * Tk_ConfigureWidget. */
{
   XGCValues gcValues;
   unsigned long valuemask;

    if (Tk_ConfigureWidget(interp, RasterPtr->tkwin, configSpecs,
	    argc, argv, (char *) RasterPtr, flags) != TCL_OK) {
	return TCL_ERROR;
    }

    /*
     *  Check width and height for reasonable values
     */

    if (RasterPtr->width < 1 || RasterPtr->width > 4096) {
       Tcl_AppendResult (interp, "Illegal raster width", (char*)NULL);
       return TCL_ERROR;
    }

    if (RasterPtr->height < 1 || RasterPtr->height > 4096) {
       Tcl_AppendResult (interp, "Illegal raster height", (char*)NULL);
       return TCL_ERROR;
    }

    /*
     * Set the background for the window and create a graphics context
     * for use during redisplay and another for drawing on the pixmap
     */

    Tk_SetWindowBackground (RasterPtr->tkwin,
			    Tk_3DBorderColor(RasterPtr->bgBorder)->pixel);


    if (RasterPtr->copyGC == None) {
       valuemask = GCFunction|GCGraphicsExposures|GCForeground;
       gcValues.function = GXcopy;
       gcValues.graphics_exposures = False;
       gcValues.foreground = Tk_3DBorderColor(RasterPtr->bgBorder)->pixel;
       RasterPtr->copyGC = Tk_GetGC(RasterPtr->tkwin,valuemask, &gcValues);
    }

    if (RasterPtr->drawGC == None) {
       valuemask = GCFunction|GCGraphicsExposures|GCForeground|GCBackground;
       gcValues.function = GXcopy;
       gcValues.background = Tk_3DBorderColor(RasterPtr->bgBorder)->pixel;
       gcValues.foreground = RasterPtr->fgColor->pixel;
       gcValues.graphics_exposures = False;
       RasterPtr->drawGC = XCreateGC (RasterPtr->display,
				      RootWindowOfScreen (
					     Tk_Screen (RasterPtr->tkwin)),
				      valuemask, &gcValues);
    }

    /*
     * Create a pixmap for storing the raster
     */
    if (RasterPtr->pm == None) {
       /* 7/1/99 johnt - changed X Pixmap functions to TK functions for Windows portability */
       RasterPtr->pm = Tk_GetPixmap (RasterPtr->display,
				      RootWindowOfScreen (
					    Tk_Screen (RasterPtr->tkwin)),
				      RasterPtr->width, RasterPtr->height,
				      DefaultDepthOfScreen(Tk_Screen(
						      RasterPtr->tkwin)));
       XFillRectangle (RasterPtr->display, RasterPtr->pm, RasterPtr->copyGC,
		       0, 0, RasterPtr->width, RasterPtr->height);
    }
    else {
       /* See if pixmap needs to be resized */
       unsigned int wid, hgt, udummy;
       int idummy;
       Window wdummy;
       XGetGeometry (RasterPtr->display, RasterPtr->pm, &wdummy,
		     &idummy, &idummy, &wid, &hgt, &udummy, &udummy);
       if (wid != RasterPtr->width || hgt != RasterPtr->height) {
#ifdef DEBUG
	  printf("Resize to %d x %d width %d height %d\n", wid, hgt,
		 RasterPtr->width, RasterPtr->height);
#endif
	  ResizeRaster(RasterPtr, wid, hgt);
/*
	  Pixmap newpm = XCreatePixmap (RasterPtr->display,
	      RootWindowOfScreen (Tk_Screen (RasterPtr->tkwin)),
	      RasterPtr->width, RasterPtr->height,
	      DefaultDepthOfScreen(Tk_Screen(RasterPtr->tkwin)));
	  XFillRectangle (RasterPtr->display, newpm, RasterPtr->copyGC,
			  0, 0, RasterPtr->width, RasterPtr->height);
	  XCopyArea (RasterPtr->display, RasterPtr->pm, newpm,
		     RasterPtr->copyGC, 0, 0, wid, hgt, 0, 0);
	  XFreePixmap (RasterPtr->display, RasterPtr->pm);
	  RasterPtr->pm = newpm;
*/
       }
    }


    /*
     * If geometry was specified, take this as the intended size for the
     * window through which the raster will be viewed. otherwise, try to
     * get a window large enough for the raster plus the border
     */
    if (RasterPtr->geometry == NULL) {
       RasterPtr->winWidth = RasterPtr->width+2*RasterPtr->borderWidth;
       RasterPtr->winHeight = RasterPtr->height+2*RasterPtr->borderWidth;
    }
    else {
       int h, w;
       if (sscanf (RasterPtr->geometry, "%dx%d", &h, &w) != 2) {
	  Tcl_AppendResult(interp, "bad geometry \"", RasterPtr->geometry,
			   "\": expected widthxheight", (char *) NULL);
	  return TCL_ERROR;
       }
       RasterPtr->winWidth = w;
       RasterPtr->winHeight = h;
    }

    RasterPtr->initialised = 1;

    /*
     * Register the desired geometry for the window.  Then arrange for
     * the window to be redisplayed.
     */
    Tk_GeometryRequest(RasterPtr->tkwin,
		       RasterPtr->winWidth,
		       RasterPtr->winHeight);
    Tk_SetInternalBorder(RasterPtr->tkwin, RasterPtr->borderWidth);
    arrangeDisplay (RasterPtr, LOW, LOW, HIGH, HIGH, RASTER_DRAW_ALL);

#ifdef DEBUG
    printf("winheight %d \n", RasterPtr->winHeight);
#endif
    return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * RasterEventProc --
 *
 *	This procedure is invoked by the Tk dispatcher for various
 *	events on Rasters.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	When the window gets deleted, internal structures get
 *	cleaned up.  When it gets exposed, it is redisplayed.
 *
 *--------------------------------------------------------------
 */

static void
RasterEventProc(clientData, eventPtr)
    ClientData clientData;	/* Information about window. */
    XEvent *eventPtr;		/* Information about event. */
{
   Raster *RasterPtr = (Raster *) clientData;

   if (!RasterPtr->initialised)
       return;

   if (eventPtr->type == Expose) {
      if (eventPtr->xexpose.count >= 0) {
	 arrangeExpose (RasterPtr, eventPtr->xexpose.x, eventPtr->xexpose.y,
			eventPtr->xexpose.x+eventPtr->xexpose.width,
			eventPtr->xexpose.y+eventPtr->xexpose.height);
      }
   } else if (eventPtr->type == ConfigureNotify) {
/*
       int wid = RasterPtr->width;
       int hgt = RasterPtr->height;
*/
       unsigned int wid, hgt, udummy;
       int idummy;
       Window wdummy;
       XGetGeometry (RasterPtr->display, RasterPtr->pm, &wdummy,
		     &idummy, &idummy, &wid, &hgt, &udummy, &udummy);
#ifdef DEBUG
       printf("old %d\n", RasterPtr->height);
#endif
#if 1
       RasterPtr->width = Tk_Width(RasterPtr->tkwin);
       RasterPtr->height = Tk_Height(RasterPtr->tkwin);
#ifdef DEBUG
       printf("Resize to %d x %d (%d x %d)\n", wid, hgt,
	    RasterPtr->width, RasterPtr->height);
#endif
       ResizeRaster(RasterPtr, wid, hgt);
#else
       RasterPtr->width = wid;
       RasterPtr->height = hgt;
       ResizeRaster(RasterPtr,
		    Tk_Width(RasterPtr->tkwin),
		    Tk_Height(RasterPtr->tkwin));
#endif

#ifdef DEBUG
       printf("new %d\n", RasterPtr->height);
#endif
      arrangeDisplay (RasterPtr, LOW, LOW, HIGH, HIGH, RASTER_DRAW_ALL);
   } else if (eventPtr->type == DestroyNotify) {
      Tcl_DeleteCommand(RasterPtr->interp, Tk_PathName(RasterPtr->tkwin));
      RasterPtr->tkwin = NULL;
      if (RasterPtr->updatePending) {
	 Tcl_CancelIdleCall(DisplayRaster, (ClientData) RasterPtr);
      }
      Tcl_EventuallyFree((ClientData) RasterPtr, DestroyRaster);
   }
}

/*
 *--------------------------------------------------------------
 *
 * DisplayRaster --
 *
 *	This procedure redraws the contents of a Raster window.
 *	It is invoked as a do-when-idle handler, so it only runs
 *	when there's nothing else for the application to do.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Information appears on the screen.
 *
 *--------------------------------------------------------------
 */

static void DisplayRaster(ClientData clientData)
{
    Raster *RasterPtr = (Raster *) clientData;
    Tk_Window tkwin = RasterPtr->tkwin;
    Pixmap pm = None;
    Drawable d;
    int result, winwid, winhgt, wid, hgt;
    XRectangle clip;
    int mode = RasterPtr->updatePending;
    int scrolling = 0;

#ifdef DEBUG
    printf("DisplayRaster update %d\n", RasterPtr->updatePending);
    printf("0:old=%d,%d new=%d,%d\n",
	   RasterPtr->x,RasterPtr->y,
	   RasterPtr->new_x,RasterPtr->new_y);
#endif
    RasterPtr->updatePending = 0;
    if (tkwin == NULL || !Tk_IsMapped(tkwin)) {
	return;
    }

#if 0
    /*
     * Analyze the Raster's redisplay rectangle. If we're not asked to
     * redisplay the whole pixmap then do a quick and dirty job.
     */
    if (RasterPtr->x0 > 0 || RasterPtr->y0 > 0 ||
	RasterPtr->x1 <= RasterPtr->width ||
	RasterPtr->y1 <= RasterPtr->height) {

       winx = RasterPtr->x + RasterPtr->borderWidth + RasterPtr->x0;
       winy = RasterPtr->y + RasterPtr->borderWidth + RasterPtr->y0;

       if (winx < RasterPtr->borderWidth) {
	  RasterPtr->x0 += RasterPtr->borderWidth - winx;
	  winx = RasterPtr->borderWidth;
       }
       if (winy < RasterPtr->borderWidth) {
	  RasterPtr->y0 += RasterPtr->borderWidth - winy;
	  winy = RasterPtr->borderWidth;
       }
       winwid = Tk_Width (tkwin) - 2*RasterPtr->borderWidth;
       winhgt = Tk_Height (tkwin) - 2*RasterPtr->borderWidth;
       if (RasterPtr->x1 + RasterPtr->x > winwid) {
	  RasterPtr->x1 = winwid - RasterPtr->x;
       }
       if (RasterPtr->y1 + RasterPtr->y > winhgt) {
	  RasterPtr->y1 = winhgt - RasterPtr->y;
       }
       wid = RasterPtr->x1 - RasterPtr->x0;
       hgt = RasterPtr->y1 - RasterPtr->y0;
       if (wid > 0 && hgt > 0) {
	  XCopyArea (Tk_Display(tkwin), RasterPtr->pm, Tk_WindowId(tkwin),
		     RasterPtr->copyGC,
		     RasterPtr->x0, RasterPtr->y0, wid, hgt,
		     RasterPtr->x0, RasterPtr->y0);
       }
       RasterPtr->x0 = HIGH;
       RasterPtr->y0 = HIGH;
       RasterPtr->x1 = LOW;
       RasterPtr->y1 = LOW;
       return;
    }
#endif

    /*
     * Create a pixmap for double-buffering, if necessary.
     */

    if (RasterPtr->doubleBuffer) {
	/* 7/1/99 johnt - changed X Pixmap functions to TK functions for Windows portability */
	pm = Tk_GetPixmap(Tk_Display(tkwin), Tk_WindowId(tkwin),
			   Tk_Width(tkwin), Tk_Height(tkwin),
			   DefaultDepthOfScreen(Tk_Screen(tkwin)));
	d = pm;
    } else {
	d = Tk_WindowId(tkwin);
    }

    /*
     * Set a clipping region if applicable
     */

    if (RasterPtr->doclip) {
       clip.x = RasterPtr->cx0;
       clip.y = RasterPtr->cy0;
       clip.width = RasterPtr->cx1-RasterPtr->cx0;
       clip.height = RasterPtr->cy1-RasterPtr->cy0;
       XSetClipRectangles (Tk_Display(tkwin), RasterPtr->copyGC, 0, 0,
			   &clip, 1, Unsorted);
    }

    if (mode & RASTER_DRAW_BORDER) {
	/*
	 * Redraw the widget's background and border.
	 */
	Tk_Draw3DRectangle(tkwin, d, RasterPtr->bgBorder,
			   0, 0, Tk_Width(tkwin), Tk_Height(tkwin),
			   RasterPtr->borderWidth, RasterPtr->relief);
    }

    if (mode & RASTER_DRAW_PIXMAP) {
	/*
	 * Display the Raster.
	 */
	wid = RasterPtr->width;
	hgt = RasterPtr->height;

	/* scrolling */
	if (RasterPtr->new_x != RasterPtr->x ||
	    RasterPtr->new_y != RasterPtr->y) {
	    int dx, dy;

	    /* initialise */
	    if (RasterPtr->plot_func) {
		RasterPtr->plot_func(RasterPtr, Tk_PathName(RasterPtr->tkwin),
				     RASTER_INIT,
				     RasterPtr->wx_start, RasterPtr->wy_start,
				     RasterPtr->wx_end, RasterPtr->wy_end);
	    }

	    scrolling = 1;
	    dx = RasterPtr->new_x - RasterPtr->x;
	    dy = RasterPtr->new_y - RasterPtr->y;

#ifdef DEBUG
	    printf("1:old=%d,%d new=%d,%d\n",
		   RasterPtr->x,RasterPtr->y,
		   RasterPtr->new_x,RasterPtr->new_y);
	    printf("dx %d wid %d wid-dx %d\n", dx, wid, wid-dx);
#endif
	    /* FIXME: - overflow in XCopyArea? */
	    /* width & height (wid-dx, hgt) are of type unsigned int */
	    /* only scroll in x if need to */
	    if (RasterPtr->new_x != RasterPtr->x) {
		if (RasterPtr->new_x > RasterPtr->x) {
		    if (wid > dx) {
#ifdef DEBUG
			printf("XCopyArea %d %d %d %d\n", dx, 0, wid-dx, hgt);
#endif
			/* Copying right chunk to left, clear right edge */
			XCopyArea (Tk_Display(tkwin), RasterPtr->pm,
				   RasterPtr->pm, RasterPtr->copyGC,
				   dx, 0,
				   wid - dx, hgt,
				   0, 0);
		    }
#ifdef DEBUG
		    printf("XFillRec %d %d %d %d\n",
			   wid - dx < 0 ? 0 : wid - dx, 0, wid, hgt);
#endif
		    XFillRectangle (RasterPtr->display, RasterPtr->pm,
				    RasterPtr->copyGC,
				    wid - dx < 0 ? 0 : wid - dx, 0,
				    wid, hgt);
		} else if (RasterPtr->new_x <= RasterPtr->x){

		    /* Copying left chunk to right, clear left edge */
		    if (wid + dx > 0) {
			XCopyArea (Tk_Display(tkwin), RasterPtr->pm,
				   RasterPtr->pm, RasterPtr->copyGC,
				   0, 0,
				   wid + dx, hgt,
				   -dx, 0);
		    }
		    XFillRectangle (RasterPtr->display, RasterPtr->pm,
				    RasterPtr->copyGC,
				    0, 0,
				    -dx > wid ? wid : -dx, hgt);
		}
		RasterScrollX(RasterPtr, RasterPtr->new_x, RasterPtr->x);
	    }
	    /* only scroll in y if need to */
	    if (RasterPtr->new_y != RasterPtr->y) {
		if (RasterPtr->new_y > RasterPtr->y) {
		    /* Copying bottom chunk to top, clear bottom edge */
		    XCopyArea (Tk_Display(tkwin), RasterPtr->pm, RasterPtr->pm,
			       RasterPtr->copyGC,
			       0, dy,
			       wid, hgt - dy,
			       0, 0);

		    XFillRectangle (RasterPtr->display, RasterPtr->pm,
				    RasterPtr->copyGC,
				    0, hgt - dy < 0 ? 0 : hgt - dy,
				    wid, hgt);
		} else if (RasterPtr->new_y <= RasterPtr->y){
		    /* Copying top chunk to bottom, clear top edge */
		    XCopyArea (Tk_Display(tkwin), RasterPtr->pm, RasterPtr->pm,
			       RasterPtr->copyGC,
			       0, 0,
			       wid, hgt + dy,
			       0, -dy);

		    XFillRectangle (RasterPtr->display, RasterPtr->pm,
				    RasterPtr->copyGC,
				    0, 0,
				    wid, -dy > hgt ? hgt : -dy);
		}
		RasterScrollY(RasterPtr, RasterPtr->new_y, RasterPtr->y);
#ifdef DEBUG
		printf("after scrolling \n");
#endif
	    }
	}
	winwid = Tk_Width (tkwin) - 2*RasterPtr->borderWidth;
	winhgt = Tk_Height (tkwin) - 2*RasterPtr->borderWidth;
	if (wid > 0 && hgt > 0) {
	    XCopyArea (Tk_Display(tkwin), RasterPtr->pm, d, RasterPtr->copyGC,
		       0, 0, winwid, winhgt,
		       RasterPtr->borderWidth, RasterPtr->borderWidth);
	}

	/*
	 * If double-buffered, copy to the screen and release the pixmap.
	 */

	if (RasterPtr->doubleBuffer) {
	    XCopyArea(Tk_Display(tkwin), pm, Tk_WindowId(tkwin),
		      RasterPtr->copyGC,
		      0, 0, Tk_Width(tkwin), Tk_Height(tkwin), 0, 0);
	    /* 7/1/99 johnt - changed X Pixmap functions to TK functions for Windows portability */
	    Tk_FreePixmap(Tk_Display(tkwin), pm);
	}
    }

    RasterPtr->x0 = HIGH;
    RasterPtr->y0 = HIGH;
    RasterPtr->x1 = LOW;
    RasterPtr->y1 = LOW;

    /*
     * Reset the clipping stuff
     */
    if (RasterPtr->doclip) {
       RasterPtr->doclip = 0;
       XSetClipMask (Tk_Display (tkwin), RasterPtr->copyGC, None);
    }

    RasterPtr->x = RasterPtr->new_x;
    RasterPtr->y = RasterPtr->new_y;

    if (/*scrolling &&*/ RasterPtr->xScrollCmd != NULL) {
       char args [200];
       double world_length = RasterPtr->wx_end - RasterPtr->wx_start;

#ifdef DEBUG_DISPLAY
       printf("scroll wx0 %f wx1 %f %f %f %f\n", RasterPtr->wx0,
	      RasterPtr->wx1,
	      (RasterPtr->wx0 - RasterPtr->wx_start)/ world_length,
	      (RasterPtr->wx1 - RasterPtr->wx_start)/ world_length,
	      world_length);
#endif
       sprintf(args, " %.20f %.20f",
	       (RasterPtr->wx0 - RasterPtr->wx_start) / world_length,
	       (RasterPtr->wx1 - RasterPtr->wx_start) / world_length);

       result = Tcl_VarEval (RasterPtr->interp, RasterPtr->xScrollCmd,
			     args, (char*) NULL);

       if (result != TCL_OK) {
	  Tcl_BackgroundError (RasterPtr->interp);
       }
       Tcl_ResetResult (RasterPtr->interp);

    }

    if (/*scrolling &&*/ RasterPtr->yScrollCmd != NULL) {
       char args [200];
       double world_height = RasterPtr->wy_end - RasterPtr->wy_start;

#ifdef DEBUG_DISPLAY
       printf("start %f end %f height %f\n", RasterPtr->wy_start, RasterPtr->wy_end,
	      world_height);

       printf("scroll wy0 %f wy1 %f %f %f %f \n",
	      RasterPtr->wy0, RasterPtr->wy1,
	      (RasterPtr->wy0- RasterPtr->wy_start) / world_height,
	      (RasterPtr->wy1- RasterPtr->wy_start) / world_height,
	      world_height);
#endif
       sprintf(args, " %.20f %.20f",
	       (RasterPtr->wy0 - RasterPtr->wy_start) / world_height,
	       (RasterPtr->wy1 - RasterPtr->wy_start) / world_height);


       result = Tcl_VarEval (RasterPtr->interp, RasterPtr->yScrollCmd,
			     args, (char*) NULL);
       if (result != TCL_OK) {
	  Tcl_BackgroundError (RasterPtr->interp);
       }
       Tcl_ResetResult (RasterPtr->interp);
    }

#ifdef DEBUG
    printf("2:old=%d,%d new=%d,%d\n",
	   RasterPtr->x,RasterPtr->y,
	   RasterPtr->new_x,RasterPtr->new_y);
#endif
}

/*
 *----------------------------------------------------------------------
 *
 * DestroyRaster --
 *
 *	This procedure is invoked by Tcl_EventuallyFree or Tcl_Release
 *	to clean up the internal structure of a Raster at a safe time
 *	(when no-one is using it anymore).
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Everything associated with the Raster is freed up.
 *
 *----------------------------------------------------------------------
 */

static void
DestroyRaster(char *clientData)
{
    Raster *RasterPtr = (Raster *) clientData;
    ClientData * dataPtr;
    Tcl_HashSearch search;
    Tcl_HashEntry* entryPtr;
    RasterImplement* primitivePtr;
    DrawEnvironment* drawEnvPtr;

    /* Free primitives */
    for (entryPtr = Tcl_FirstHashEntry (&PrimitiveTable, &search);
	 entryPtr != NULL;
	 entryPtr = Tcl_NextHashEntry (&search)) {
       primitivePtr = (RasterImplement*) Tcl_GetHashValue (entryPtr);
       if (primitivePtr->free != NULL) {
	  dataPtr = &(RasterPtr->primitiveData[primitivePtr->primitiveno]);
	  (*(primitivePtr->free)) (RasterPtr, *dataPtr);
       }
    }
    free (RasterPtr->primitiveData);

    /* Free DrawEnvironments */
    for (entryPtr = Tcl_FirstHashEntry (&RasterPtr->drawEnvTable, &search);
	 entryPtr != NULL;
	 entryPtr = Tcl_NextHashEntry (&search)) {
       drawEnvPtr = (DrawEnvironment*) Tcl_GetHashValue (entryPtr);
       DestroyDrawEnv (RasterPtr, drawEnvPtr);
    }
    Tcl_DeleteHashTable (&RasterPtr->drawEnvTable);

    /* Free Conf Options */
    Tk_FreeOptions(configSpecs, (char *) RasterPtr, RasterPtr->display, 0);

    /* Free X resources */
    if (RasterPtr->copyGC != None) {
	Tk_FreeGC(RasterPtr->display, RasterPtr->copyGC);
	XFreeGC (RasterPtr->display, RasterPtr->drawGC);
    }
    if (RasterPtr->pm != None) {
       /* 7/1/99 johnt - changed X Pixmap functions to TK functions for Windows portability */
       Tk_FreePixmap (RasterPtr->display, RasterPtr->pm);
    }
    ckfree((char *) RasterPtr);
}


/*=============================================================================
 */

void RasterRefresh(Raster *RasterPtr) {
   if (RasterPtr->knownarea) {
      /* Primitive informed modified area of pixmap */
      if (RasterPtr->px1 >= 0 && RasterPtr->py1 >= 0) {
	 /* No point in redisplaying an empty area */
	 arrangeDisplay (RasterPtr, RasterPtr->px0, RasterPtr->py0,
			 RasterPtr->px1, RasterPtr->py1,
			 RASTER_DRAW_PIXMAP);
      }
   }
   else {
      /* Redisplay the whole pixmap */
      arrangeDisplay(RasterPtr, 0, 0, RasterPtr->width-1, RasterPtr->height-1,
		     RASTER_DRAW_PIXMAP);
   }

   RasterPtr->px0 = HIGH;
   RasterPtr->py0 = HIGH;
   RasterPtr->px1 = LOW;
   RasterPtr->py1 = LOW;
   RasterPtr->knownarea = 0;
}

void tk_RasterRefresh(Tk_Raster *rasterptr)
{
   Raster* RasterPtr = (Raster*) rasterptr;

   RasterRefresh(RasterPtr);
}

/*
 *______________________________________________________________________
 *
 * arrangeDisplay --
 *
 *	This procedure is invoked whenever we discover that a part of
 *      the pixmap has been modified and/or a part of the window must
 *      be redrawn.
 *
 * Results:
 *
 *	None.
 *
 * Side effects:
 *
 *      A Call to DisplayRaster is eventually scheduled. The Raster
 *      structure pointed to by RasterPtr is updatedd to reflect the
 *      area that must be redisplayed.
 *
 *----------------------------------------------------------------------
 */

static void
arrangeDisplay (RasterPtr, x0, y0, x1, y1, type)
     Raster * RasterPtr;
     int x0, y0, x1, y1, type;
{
#ifdef DEBUG
    printf("arrangeDisplay update %d type %d\n", RasterPtr->updatePending,
	   type);
#endif
   RasterPtr->doclip = 0;
   if (x0 < RasterPtr->x0) RasterPtr->x0 = x0;
   if (y0 < RasterPtr->y0) RasterPtr->y0 = y0;
   if (x1 > RasterPtr->x1) RasterPtr->x1 = x1;
   if (y1 > RasterPtr->y1) RasterPtr->y1 = y1;
   if (!RasterPtr->updatePending) {
       Tcl_DoWhenIdle(DisplayRaster, (ClientData) RasterPtr);
   }
   RasterPtr->updatePending |= type;
}

/*
 *______________________________________________________________________
 *
 * arrangeExpose --
 *
 *	This procedure is invoked when we get an Expose event. It
 * 	arranges for the redisplay routine to be invoked using a
 *  	clip region constraining it to the exposed part of the window.
 *
 * Results:
 *
 *	None.
 *
 * Side effects:
 *
 *      A Call to DisplayRaster is eventually scheduled. The Raster
 *      structure pointed to by RasterPtr is updated to reflect the
 *      area that must be redisplayed.
 *
 *----------------------------------------------------------------------
 */

static void
arrangeExpose (RasterPtr, x0, y0, x1, y1)
     Raster * RasterPtr;
     int x0, y0, x1, y1;
{

#ifdef DEBUG
    printf("arrangeExpose update %d clip %d \n", RasterPtr->updatePending,
	   RasterPtr->doclip);
#endif
    if (RasterPtr->updatePending) {
	if (RasterPtr->doclip) {
	    /* Extend the clipping rectangle */
	    if (x0 < RasterPtr->cx0) RasterPtr->cx0 = x0;
	    if (y0 < RasterPtr->cy0) RasterPtr->cy0 = y0;
	    if (x1 > RasterPtr->cx1) RasterPtr->cx1 = x1;
	    if (y1 > RasterPtr->cy1) RasterPtr->cy1 = y1;
	}
	RasterPtr->x0 = RasterPtr->y0 = LOW;
	RasterPtr->x1 = RasterPtr->y1 = HIGH;
    } else {
	RasterPtr->doclip = 1;
	RasterPtr->cx0 = x0;
	RasterPtr->cy0 = y0;
	RasterPtr->cx1 = x1;
	RasterPtr->cy1 = y1;
	RasterPtr->x0 = RasterPtr->y0 = LOW;
	RasterPtr->x1 = RasterPtr->y1 = HIGH;
	Tcl_DoWhenIdle(DisplayRaster, (ClientData) RasterPtr);
    }
    RasterPtr->updatePending = RASTER_DRAW_ALL;
}

/*--------------------------------------------------------------------------
 *
 *  Configuration options for a Raster "Drawing Environment"
 */

Tk_ConfigSpec DrawEnvSpecs [] = {

    /*
    {TK_CONFIG_INT, "-bgpixel", "bgpixel",  "Bgpixel",
     (char*) NULL, Tk_Offset(DrawEnvironment, bgColor->pixel), 0},
    {TK_CONFIG_INT, "-fgpixel", "fgpixel",  "Fgpixel",
     (char*) NULL, Tk_Offset(DrawEnvironment, fgColor->pixel), 0},
     */
    {TK_CONFIG_COLOR, "-background", "background",  "Background",
	NULL, Tk_Offset(DrawEnvironment, bgColor),  TK_CONFIG_NULL_OK},
    {TK_CONFIG_COLOR, "-foreground", "foreground", "Foreground",
	NULL, Tk_Offset(DrawEnvironment, fgColor), TK_CONFIG_NULL_OK},
    {TK_CONFIG_PIXELS, "-linewidth", "linewidth", "LineWidth", (char*) NULL,
	Tk_Offset(DrawEnvironment, gcValues.line_width), 0 },
    {TK_CONFIG_CUSTOM, "-linestyle", "linestyle", "LineStyle", "solid",
	Tk_Offset(DrawEnvironment, gcValues.line_style),
	TK_CONFIG_DONT_SET_DEFAULT, &LineStyleOption },
    {TK_CONFIG_CAP_STYLE, "-capstyle", "capstyle", "CapStyle",(char*) NULL,
	Tk_Offset(DrawEnvironment, gcValues.cap_style), 0},
    {TK_CONFIG_JOIN_STYLE, "-joinstyle", "joinstyle","JoinStyle",(char*) NULL,
	Tk_Offset(DrawEnvironment, gcValues.join_style), 0},
    {TK_CONFIG_CUSTOM, "-function", "function", "Function", "copy",
	Tk_Offset(DrawEnvironment, gcValues.function),
	TK_CONFIG_DONT_SET_DEFAULT, &FunctionOption },
    {TK_CONFIG_BITMAP, "-stipple", "stipple", "Stipple", (char*) NULL,
	Tk_Offset(DrawEnvironment, gcValues.stipple), 0},
    {TK_CONFIG_CUSTOM, "-fillstyle", "fillstyle", "FillStyle", "solid",
	Tk_Offset(DrawEnvironment, gcValues.fill_style),
	TK_CONFIG_DONT_SET_DEFAULT,
	&FillStyleOption },
    {TK_CONFIG_SYNONYM, "-bg", "background", (char *) NULL,
	(char *) NULL, 0, 0},
    {TK_CONFIG_SYNONYM, "-fg", "foreground", (char *) NULL,
	(char *) NULL, 0, 0},
    {TK_CONFIG_END, (char *) NULL, (char *) NULL, (char *) NULL,
	(char *) NULL, 0, 0}
};

int
GetFgPixel(Tcl_Interp* interp,
	   Tk_Raster* rasterptr,
	   int index)
{
    DrawEnvironment* drawEnv;
    Raster* RasterPtr = (Raster*) rasterptr;

    if (DrawEnvIndex (interp, RasterPtr, index,
		      (ClientData*) &drawEnv) != TCL_OK) {
	return -1;
    }
    return drawEnv->gcValues.foreground;
}


int GetBgPixel(Tk_Raster* rasterptr)
{
    Raster* RasterPtr = (Raster*) rasterptr;
    return (Tk_3DBorderColor(RasterPtr->bgBorder)->pixel);
}

int
SetFgPixel(Tcl_Interp* interp,
	   Tk_Raster* rasterptr,
	   int index,
	   int fg_pixel)
{
    DrawEnvironment* drawEnv;
    Raster* RasterPtr = (Raster*) rasterptr;

    if (DrawEnvIndex (interp, RasterPtr, index,
		      (ClientData*) &drawEnv) != TCL_OK) {
	return -1;
    }
    drawEnv->gcValues.foreground = fg_pixel;
    return 0;
}

static int ConfigDrawEnv (interp, RasterPtr, drawEnv, argc, argv)
/*         -------------
 *
 *  Processes configuration options related to Drawing Environments.
 *  valMask and gcValues entries in the given DrawEnviroment structure
 *  are modified according to the configuration options. Returns
 *  a standard Tcl result.
 */
     Tcl_Interp* interp;
     Raster* RasterPtr;
     DrawEnvironment* drawEnv;
     int argc;
     char * argv [];
{
   if (Tk_ConfigureWidget(interp, RasterPtr->tkwin, DrawEnvSpecs,
      argc, argv, (char *) drawEnv, TK_CONFIG_ARGV_ONLY) != TCL_OK) {
      return TCL_ERROR;
   }

   if (drawEnv->bgColor) {
      /* Background was changed */
      drawEnv->valMask |= GCBackground;
      drawEnv->gcValues.background = drawEnv->bgColor->pixel;
   }
   if (drawEnv->fgColor) {
      /* Foreground was changed */
      drawEnv->valMask |= GCForeground;
      drawEnv->gcValues.foreground = drawEnv->fgColor->pixel;
   }
   if (DrawEnvSpecs [4].specFlags & TK_CONFIG_OPTION_SPECIFIED) {
      /* Line width was changed */
      drawEnv->valMask |= GCLineWidth;
   }
   if (DrawEnvSpecs [5].specFlags & TK_CONFIG_OPTION_SPECIFIED) {
      /* Line width was changed */
      drawEnv->valMask |= GCLineStyle;
   }
   if (DrawEnvSpecs [6].specFlags & TK_CONFIG_OPTION_SPECIFIED) {
      /* Cap Style was changed */
      drawEnv->valMask |= GCCapStyle;
   }
    if (DrawEnvSpecs [7].specFlags & TK_CONFIG_OPTION_SPECIFIED) {
      /* Join Style was changed */
      drawEnv->valMask |= GCJoinStyle;
   }
   if (DrawEnvSpecs [8].specFlags & TK_CONFIG_OPTION_SPECIFIED) {
      /* BITBLT function was changed */
      drawEnv->valMask |= GCFunction;
   }
   if (DrawEnvSpecs [9].specFlags & TK_CONFIG_OPTION_SPECIFIED) {
      /* Stipple bitmap was changed */
      drawEnv->valMask |= GCStipple;
   }
   if (DrawEnvSpecs [10].specFlags & TK_CONFIG_OPTION_SPECIFIED) {
      /* FillStyle was changed */
      drawEnv->valMask |= GCFillStyle;
   }
   
   drawEnv->drawGC = XCreateGC(RasterPtr->display, 
	                      RootWindowOfScreen(Tk_Screen (RasterPtr->tkwin)),
			              drawEnv->valMask, &(drawEnv->gcValues));

   if (drawEnv == RasterPtr->currentDrawEnv) {
      return SetDrawEnv (interp, RasterPtr, drawEnv);
   }

   return TCL_OK;
}

/* removed static kfs 10.05.00 */
int ConfigInfoDrawEnv (Tcl_Interp* interp,
		       Tk_Raster* raster,
		       DrawEnvironment*drawEnv,
		       int argc,
		       char **argv)
/*         -----------------
 *
 *  Returns a string describing configuration options of a drawing
 *  environment.
 */
{
    Raster* RasterPtr = (Raster*) raster;

   if (argc == 0) {
      return Tk_ConfigureInfo (interp, RasterPtr->tkwin,
			       DrawEnvSpecs, (char*) drawEnv, (char*)NULL, 0);
   }
   return Tk_ConfigureInfo (interp, RasterPtr->tkwin,
			    DrawEnvSpecs, (char*)drawEnv, argv [0], 0);
}

int DrawEnvIndex (interp, rasterptr, drawenvno, drawenvptrptr)
/*  ------------
 *
 *  Searches rasterptr's drawEnvTable for a DrawEnv for a drawing environment
 *  with the given 'drawenvno'. If successful, puts in *drawenvptr a pointer
 *  to that DrawEnv. Returns a standard TCL result
 */
     Tcl_Interp * interp;
     Tk_Raster* rasterptr;
     int drawenvno;
     ClientData* drawenvptrptr;
{
   Raster* RasterPtr = (Raster*) rasterptr;
   DrawEnvironment** drawEnvPtrPtr = (DrawEnvironment**) drawenvptrptr;
   Tcl_HashEntry *entryPtr;

   entryPtr = Tcl_FindHashEntry (&(RasterPtr->drawEnvTable), (char*)drawenvno);
   if (entryPtr == NULL) {
      Tcl_AppendResult (interp, "Given drawing environment does not exist",
			(char*) NULL);
      return TCL_ERROR;
   }
   *drawEnvPtrPtr = (DrawEnvironment*) Tcl_GetHashValue (entryPtr);
   return TCL_OK;
}


void GetDrawEnv (rasterptr, drawenvptr)
/*   ----------
 *
 *  Returns in * drawenvptr a pointer to the raster's current drawing env
 *
 */
     Tk_Raster* rasterptr;
     ClientData* drawenvptr;
{
   *drawenvptr = (ClientData)((Raster*) rasterptr)->currentDrawEnv;
}

int SetDrawEnv (interp, rasterptr, drawenvptr)
/*  ----------
 *
 *  Change the 'drawGC' of 'rasterptr' to reflect the drawing environment
 *  given by drawenvptr.
 *
 */
     Tcl_Interp * interp;
     Tk_Raster* rasterptr;
     ClientData drawenvptr;
{
   Raster* RasterPtr = (Raster*) rasterptr;
   DrawEnvironment* drawEnvPtr = (DrawEnvironment*) drawenvptr;

   if (drawEnvPtr->valMask != 0) {
        RasterPtr->drawGC = drawEnvPtr->drawGC;
   }

   RasterPtr->currentDrawEnv = drawEnvPtr;

   return TCL_OK;
}


static int CreateDrawEnv (interp, RasterPtr, argc, argv)
/*  -------------
 *
 *  Creates a new drawing environment for the given raster, setting it
 *  according to the options given in argc and argv. If successful,
 *  a pointer to the just created drawing environment is put in *drawenvptrptr.
 *  Returns a standard tcl result
 *
 */
     Tcl_Interp * interp;
     Raster* RasterPtr;
     int argc;
     char* argv [];
{
   Tcl_HashEntry *entryPtr;
   DrawEnvironment *drawEnvPtr;
   int new;

   drawEnvPtr = (DrawEnvironment*) malloc (sizeof (DrawEnvironment));
   drawEnvPtr->index = RasterPtr->drawEnvCount;
   drawEnvPtr->fgColor = Tk_GetColor (interp, RasterPtr->tkwin, "brown");
   drawEnvPtr->bgColor = Tk_GetColor (interp, RasterPtr->tkwin, "#d9d9d9");
   drawEnvPtr->gcValues.background = drawEnvPtr->bgColor->pixel;
   drawEnvPtr->gcValues.foreground = drawEnvPtr->fgColor->pixel;
   drawEnvPtr->gcValues.stipple = None;
   drawEnvPtr->gcValues.line_width = 1;
   drawEnvPtr->gcValues.line_style = LineSolid;
   drawEnvPtr->gcValues.cap_style = CapButt;
   drawEnvPtr->gcValues.fill_style = FillSolid;
   drawEnvPtr->gcValues.join_style = JoinRound;
   drawEnvPtr->gcValues.function = GXcopy;
   drawEnvPtr->valMask = (GCBackground|GCForeground|GCLineWidth|
			  GCLineStyle|GCCapStyle|GCJoinStyle|
			  GCFunction|GCFillStyle);
   if (ConfigDrawEnv (interp, RasterPtr, drawEnvPtr, argc, argv) != TCL_OK) {
      free (drawEnvPtr);
      return TCL_ERROR;
   }

   entryPtr = Tcl_CreateHashEntry (&(RasterPtr->drawEnvTable),
				   (char*) (RasterPtr->drawEnvCount),
				   &new);

   if (!new) {
      Tcl_AppendResult (interp, "Could not create a Drawing Environment",
			(char*)NULL);
      free (drawEnvPtr);
      return TCL_ERROR;
   }

   Tcl_SetHashValue (entryPtr, drawEnvPtr);

   vTcl_SetResult(interp, "%d", RasterPtr->drawEnvCount++);
   return TCL_OK;
}

static void DestroyDrawEnv (RasterPtr, drawenvptr)
/*          --------------
 *
 * Destroys the given drawing environment freeing its options
 */
     Raster* RasterPtr;
     DrawEnvironment*  drawenvptr;
{
   Tcl_HashEntry *entryPtr;
   Tk_FreeOptions (DrawEnvSpecs, (char*) drawenvptr, RasterPtr->display, 0);

   entryPtr = Tcl_FindHashEntry (&(RasterPtr->drawEnvTable),
				 (char*)(drawenvptr->index));
   Tcl_DeleteHashEntry (entryPtr);
   free (drawenvptr);
}

/*========================================================================
 *
 *  Procedures myOptionParse and myOptionPrint implement all custom options
 *  used for configuring drawing environments.
 *
 */

static int
myOptionParse (clientData, interp, tkwin, value, widgRec, offset)
     ClientData clientData;
     Tcl_Interp * interp;
     Tk_Window tkwin;
     char* value;
     char* widgRec;
     int offset;
{
   AttrNameValue* table = (AttrNameValue*) clientData;
   AttrNameValue* ptr = table;

   while (ptr->optionname != (char*) NULL) {
      if (strcmp (value, ptr->optionname) == 0) {
	 * (int *) (widgRec+offset) = ptr->optionvalue;
	 return TCL_OK;
      }
      ptr++;
   }

   Tcl_AppendResult (interp, "wrong arg \"",
		     value, "\": should be one of ", (char*) NULL);

   ptr = table;
   while (ptr->optionname != (char*)NULL) {
      Tcl_AppendResult (interp, "\"", ptr->optionname,
			(ptr+1)->optionname == (char*)NULL ? "\"." : "\",",
			(char*) NULL);
      ptr++;
   }

   return TCL_ERROR;
}

static char *
myOptionPrint (clientData, tkwin, widgRec, offset, freeProcPtr)
     ClientData clientData;
     Tk_Window tkwin;
     char* widgRec;
     int offset;
     Tcl_FreeProc ** freeProcPtr;
{
   AttrNameValue* ptr = (AttrNameValue*) clientData;
   int value = *(int*) (widgRec+offset);
   while (ptr->optionname != (char*)NULL) {
      if (ptr->optionvalue == value) return ptr->optionname;
      ptr++;
   }
   return (char*) NULL;
}

/*======================================================================
 *
 *  Utility functions for handling world-to-raster-and-back  transformations
 */

void SetRasterCoords (Tk_Raster* raster,
		      double x0,
		      double y0,
		      double x1,
		      double y1)
/*
 *  (x0,y0) are the coordinates of the upper-left corner and
 *  (x1,y1) are the coordinates of the lower-right corner
 */
{
   Raster* RasterPtr = (Raster*) raster;

   /* prevent any assignments if values are both 0 */
   if ((x1 - x0 == 0) || (y1-y0 == 0))
       return;

   RasterPtr->ax = RasterPtr->width / (x1-x0);
   RasterPtr->ay = RasterPtr->height / (y1-y0);

/*
   printf("SetRasterCoords width %d height %d \n", RasterPtr->width,
	  RasterPtr->height);
*/
#ifdef KFS_REMOVED
   double ratio;

   ratio = RasterPtr->ax / RasterPtr->ay;
   if (ratio < 0.0) ratio = -ratio;
   if (ratio > 1.0) {
      RasterPtr->ax /= ratio;
      RasterPtr->bx = x0 - (RasterPtr->width / RasterPtr->ax + x0 - x1) / 2;
      RasterPtr->by = y0;
   }
   else {
      RasterPtr->ay *= ratio;
      RasterPtr->bx = x0;
      RasterPtr->by = y0 - (RasterPtr->height / RasterPtr->ay + y0 - y1) / 2;
   }
#endif
   /* KFS: must set bx and by */
   RasterPtr->bx = x0;
   RasterPtr->by = y0;

   RasterPtr->wx0 = x0;
   RasterPtr->wy0 = y0;
   RasterPtr->wx1 = x1;
   RasterPtr->wy1 = y1;

#ifdef DEBUG
   printf("*******SetRasterCoords x0 %f x1 %f y0 %f y1 %f\n", x0, x1, y0, y1);
   printf("set to *%f+%f\n", RasterPtr->ax, RasterPtr->bx);
   printf("set to *%f+%f\n", RasterPtr->ay, RasterPtr->by);
#endif
}

void GetWorldToRasterConversion(Tk_Raster *raster, double *ax, double *ay, double *bx, double *by) {
/*
 * get the conversion factors for use outside
 * for speed rather than anything else
 */
    Raster* RasterPtr = (Raster*) raster;
    *ax = RasterPtr->ax;
    *ay = RasterPtr->ay;
    *bx = RasterPtr->bx;
    *by = RasterPtr->by;
}     


void WorldToRaster (raster, wx, wy, rx, ry)
/*
 *  (wx,wy) are  world coordinates
 *  (*rx,*ry) are set to the corresponding raster coordinates
 */
     Tk_Raster * raster;
     double wx, wy;
     int * rx, * ry;
{
   Raster* RasterPtr = (Raster*) raster;

    *rx = (wx - RasterPtr->bx) * RasterPtr->ax;
    *ry = (wy - RasterPtr->by) * RasterPtr->ay;
#ifdef DEBUG
   printf("WorldToRaster ax %f bx %f ay %f by %f wx %f wy %f rx %d ry %d\n",
	  RasterPtr->ax, RasterPtr->bx, RasterPtr->ay, RasterPtr->by, wx, wy,
	  *rx, *ry);
#endif
}

void RasterToWorld (raster, rx, ry, wx, wy)
/*
 *  (rx,ry) are raster coordinates
 *  (*wx,*wy) are set to the corresponding world coordinates
 */
     Tk_Raster* raster;
     int rx, ry;
     double * wx, * wy;
{
   Raster* RasterPtr = (Raster*) raster;

   *wx = rx / RasterPtr->ax + RasterPtr->bx;
   *wy = ry / RasterPtr->ay + RasterPtr->by;
#ifdef DEBUG
   printf("ry %d ay %f by %f wy %f \n",
	  RasterPtr->height, ry, RasterPtr->ay, RasterPtr->by, *wy);
#endif
}

/*======================================================================
 *
 *  SetRasterModifiedArea may be used by primitive implementations
 *  to let the raster know what part of the pixmap was modified.
 *  If not called, raster will assume that all pixmap was modified.
 *  If implementing a raster command that does not actually modify the
 *  pixmap, this routine should be called with 0 0 0 0 as args.
 */

void SetRasterModifiedArea (raster, rx0, ry0, rx1, ry1)
     Tk_Raster* raster;
     int rx0, ry0;
     int rx1, ry1;
{
   Raster* RasterPtr = (Raster*) raster;
   int lw, tmp;

   if (rx0 > rx1) { tmp = rx0; rx0 = rx1; rx1 = tmp; }
   if (ry0 > ry1) { tmp = ry0; ry0 = ry1; ry1 = tmp; }
   RasterPtr->knownarea = 1;
   if (rx1 == 0 && rx0 == 0 && ry1 == 0 && ry0 == 0) return;

   lw = RasterPtr->currentDrawEnv->gcValues.line_width;
   rx0 -= lw; if (rx0 < 0) rx0 = 0;
   rx1 += lw; if (rx1 >= RasterPtr->width) rx1 = RasterPtr->width-1;
   ry0 -= lw; if (ry0 < 0) ry0 = 0;
   ry1 += lw; if (ry1 >= RasterPtr->height) ry1 = RasterPtr->height-1;

   if (rx0 < RasterPtr->px0) RasterPtr->px0 = rx0;
   if (ry0 < RasterPtr->py0) RasterPtr->py0 = ry0;
   if (rx1 > RasterPtr->px1) RasterPtr->px1 = rx1;
   if (ry1 > RasterPtr->py1) RasterPtr->py1 = ry1;
}

/*======================================================================
 *
 *  The following are utility functions to access elements of a Tk_Raster.
 *  They may be used to implement drawing primitives which have to
 *  access X functions directly
 */

Drawable GetRasterDrawable (raster)
     Tk_Raster* raster;
{
   return ((Raster*)raster)->pm;
}

Display* GetRasterDisplay (raster)
     Tk_Raster* raster;
{
   return ((Raster*)raster)->display;
}

Tk_Window GetRasterTkWin (raster)
     Tk_Raster* raster;
{
   return ((Raster*)raster)->tkwin;
}

GC GetRasterGC (raster)
     Tk_Raster* raster;
{
   return ((Raster*)raster)->drawGC;
}


/*=======================================================================
 *
 *  Below are the functions that implement the built-in drawing primitives
 */


void RasterDrawPoints (raster, coord, npts)
/*   ----------------
 *
 *  Draws the points (coord [2i], coord [2i+1]) for 0 <= i < npts
 */
    Tk_Raster * raster;
     double* coord;
     int npts;
{
   int pointwid = ((Raster*) raster)->currentDrawEnv->gcValues.line_width ;
   GC gc = GetRasterGC (raster);
   Drawable d = GetRasterDrawable (raster);
   Display * dsp = GetRasterDisplay (raster);
   int i;
   int n = 2 * npts;
   XPoint *pt, *ptptr;
   int rx, ry;
   int minx = HIGH, miny = HIGH, maxx = LOW, maxy = LOW;

   if (npts < 1) return;

   for (i = 0, ptptr = pt = (XPoint*) malloc (sizeof (XPoint)*npts);
	i < n; i+=2, ptptr++) {
      WorldToRaster (raster, coord [i], coord [i+1], &rx, &ry);
      if (rx < minx) minx = rx;
      if (rx > maxx) maxx = rx;
      if (ry < miny) miny = ry;
      if (ry > maxy) maxy = ry;
      ptptr->x = rx;
      ptptr->y = ry;
   }

   if (pointwid >= 2) {
      int halfwid = pointwid/2;
      for (i=0, ptptr = pt; i < npts; i++, ptptr++) {
	 XFillArc (dsp, d, gc, ptptr->x-halfwid, ptptr->y-halfwid,
		   pointwid, pointwid, 0, 360*64);
      }
   } else {
      XDrawPoints (dsp, d, gc, 	pt, npts, CoordModeOrigin);
   }
   free (pt);

   SetRasterModifiedArea (raster, minx, miny, maxx, maxy);
}

void RasterDrawPoint (Tk_Raster *raster,
		      int x,
		      int y)
/*   ---------------
 *
 *  Draws the point x, y
 */
{
   int pointwid = ((Raster*) raster)->currentDrawEnv->gcValues.line_width ;
   GC gc = GetRasterGC (raster);
   Drawable d = GetRasterDrawable (raster);
   Display * dsp = GetRasterDisplay (raster);
   int rx, ry;
   int minx = HIGH, miny = HIGH, maxx = LOW, maxy = LOW;

   WorldToRaster (raster, x, y, &rx, &ry);
   if (rx < minx) minx = rx;
   if (rx > maxx) maxx = rx;
   if (ry < miny) miny = ry;
   if (ry > maxy) maxy = ry;

   if (pointwid >= 2) {
      int halfwid = pointwid/2;
      XFillArc (dsp, d, gc, rx - halfwid, ry - halfwid,
		pointwid, pointwid, 0, 360*64);
   } else {
      XDrawPoint (dsp, d, gc, rx, ry);
   }

#ifdef DEBUG
   printf("x %d y %d rx %d ry %d \n", x, y, rx, ry);
#endif
   /* XDrawPoint (dsp, d, gc, rx, ry); */

   SetRasterModifiedArea (raster, minx, miny, maxx, maxy);
}

void
RasterDrawLine (Tk_Raster *raster,
		int x1, double y1, int x2, double y2)
/*   ---------------
 *
 *  Draws a single line segments connecting x1 y1 with x2 y2
 */
{
   int rx1, ry1, rx2, ry2;
   int minx = HIGH, miny = HIGH, maxx = LOW, maxy = LOW;

   WorldToRaster (raster, x1, y1, &rx1, &ry1);
   WorldToRaster (raster, x2, y2, &rx2, &ry2);

#ifdef DEBUG
   printf("line (%d,%f)-(%d,%f) [%d,%d]-[%d,%d]\n",
	  x1, y1, x2, y2,
	  rx1, ry1, rx2, ry2);
#endif
   if (rx1 < minx) minx = rx1;
   if (rx1 > maxx) maxx = rx1;
   if (ry1 < miny) miny = ry1;
   if (ry1 > maxy) maxy = ry1;

   if (rx2 < minx) minx = rx2;
   if (rx2 > maxx) maxx = rx2;
   if (ry2 < miny) miny = ry2;
   if (ry2 > maxy) maxy = ry2;

   XDrawLine (GetRasterDisplay (raster),
	      GetRasterDrawable (raster),
	      GetRasterGC (raster),
	      rx1, ry1, rx2, ry2);

   /* printf("minx %d miny %d maxx %d maxy %d\n", minx, miny, maxx, maxy); */
   SetRasterModifiedArea (raster, minx, miny, maxx, maxy);
}

void RasterDrawLines (raster, coord, npts)
/*   ---------------
 *
 *  Draws the npts-1 line segments connecting the successive
 *  'npts' points (coord [2i], coord [2i+1]), for 0 <= i < npts.
 */
     Tk_Raster * raster;
     double* coord;
     int npts;
{
   int i, j, num;
   int n = npts*2;
   XPoint *pt, *ptptr;
   int rx, ry;
   int minx = HIGH, miny = HIGH, maxx = LOW, maxy = LOW;

   if (npts < 1) return;

   for (i = 0, ptptr = pt = (XPoint*) malloc (sizeof (XPoint)*npts);
	i < n; i+=2, ptptr++) {
      WorldToRaster (raster, coord [i], coord [i+1], &rx, &ry);
      if (rx < minx) minx = rx;
      if (rx > maxx) maxx = rx;
      if (ry < miny) miny = ry;
      if (ry > maxy) maxy = ry;
      ptptr->x = rx;
      ptptr->y = ry;
   }

   /*
    * fix to deal with x servers which can't cope with more than 2^16 lines
    */
   if (npts < 65000) {
       XDrawLines (GetRasterDisplay (raster),
		   GetRasterDrawable (raster),
		   GetRasterGC (raster),
		   pt, npts, CoordModeOrigin);
   } else {
       j = 1;
       for (i = 0; i < npts; i+=65000, j = i) {
	   if (npts < i + 65000) {
	       num = npts - i + 1;
	   } else {
	       num = 65000;
	   }

	   XDrawLines (GetRasterDisplay (raster),
		       GetRasterDrawable (raster),
		       GetRasterGC (raster),
		       &pt[j-1], num, CoordModeOrigin);
       }
   }
   /*
   XDrawLines (GetRasterDisplay (raster),
	       GetRasterDrawable (raster),
	       GetRasterGC (raster),
	       pt, npts, CoordModeOrigin);
   */
   free (pt);

   SetRasterModifiedArea (raster, minx, miny, maxx, maxy);
}

void RasterDrawSegments (Tk_Raster * raster,
			 double* coord,
			 int nsegs)
/*   ---------------
 *
 *  Draws the nsegs line segments connecting pairs of successive
 *  points. coord contains nsegs*4 values.
 */
{
   int i, num;
   int n = nsegs*4;
   XSegment *seg, *segptr;
   int rx1, ry1, rx2, ry2;
   int minx = HIGH, miny = HIGH, maxx = LOW, maxy = LOW;

   if (nsegs < 1) return;

   for (i = 0, segptr = seg = (XSegment*) malloc (sizeof (XSegment)*(nsegs));
	i < n; i+=4, segptr++) {
      WorldToRaster (raster, coord [i], coord [i+1], &rx1, &ry1);
      WorldToRaster (raster, coord [i+2], coord [i+3], &rx2, &ry2);
      if (rx1 < minx) minx = rx1;
      if (rx1 > maxx) maxx = rx1;
      if (ry1 < miny) miny = ry1;
      if (ry1 > maxy) maxy = ry1;
      if (rx2 < minx) minx = rx2;
      if (rx2 > maxx) maxx = rx2;
      if (ry2 < miny) miny = ry2;
      if (ry2 > maxy) maxy = ry2;
      segptr->x1 = rx1;
      segptr->y1 = ry1;
      segptr->x2 = rx2;
      segptr->y2 = ry2;
   }

   /*
    * fix to deal with x servers which can't cope with more than 2^16 lines
    */
   if (nsegs < 32000) {
       XDrawSegments (GetRasterDisplay (raster),
		   GetRasterDrawable (raster),
		   GetRasterGC (raster),
		   seg, nsegs);
   } else {
       for (i = 0; i < nsegs; i+=32000) {
	   if (nsegs < i + 32000) {
	       num = nsegs - i;
	   } else {
	       num = 32000;
	   }

	   XDrawSegments (GetRasterDisplay (raster),
		       GetRasterDrawable (raster),
		       GetRasterGC (raster),
		       &seg[i], num);
       }
   }

   free (seg);

   SetRasterModifiedArea (raster, minx, miny, maxx, maxy);
}

void RasterDrawRectangles (raster, coord, nrects)
/*   --------------------
 *
 *  Draws rectangles with sides aligned with the coordinate axes and
 *  with diagonal corners at (coord [4i], coord [4i+1]) and
 *  (coord [4i+2], coord [4i+3]) for  0 <= i < nrects,
 */
     Tk_Raster * raster;
     double* coord;
     int nrects;
{
   int i;
   int n = nrects*4;
   XRectangle *rect, *rectptr;
   int rx, ry, tmp;
   int minx = HIGH, miny = HIGH, maxx = LOW, maxy = LOW;

   if (nrects < 1) return;

   for (i = 0, rectptr= rect= (XRectangle*) malloc(sizeof (XRectangle)*nrects);
	i < n; i+=4, rectptr++) {
      WorldToRaster (raster, coord [i], coord [i+1], &rx, &ry);
      if (rx < minx) minx = rx;
      if (rx > maxx) maxx = rx;
      if (ry < miny) miny = ry;
      if (ry > maxy) maxy = ry;
      rectptr->x = rx;
      rectptr->y = ry;
      WorldToRaster (raster, coord [i+2], coord [i+3], &rx, &ry);
      if (rx < minx) minx = rx;
      if (rx > maxx) maxx = rx;
      if (ry < miny) miny = ry;
      if (ry > maxy) maxy = ry;
      if (rectptr->x > rx) { tmp = rectptr->x; rectptr->x = rx; rx = tmp; }
      if (rectptr->y > ry) { tmp = rectptr->y; rectptr->y = ry; ry = tmp; }
      rectptr->width = rx - rectptr->x;
      rectptr->height = ry - rectptr->y;
   }

   XDrawRectangles (GetRasterDisplay (raster),
		    GetRasterDrawable (raster),
		    GetRasterGC (raster),
		    rect, nrects);
   free (rect);

   SetRasterModifiedArea (raster, minx, miny, maxx, maxy);
}

void RasterFillRectangles (raster, coord, nrects)
/*   --------------------
 *
 *  Draws filled rectangles with sides aligned with the coordinate axes and
 *  with diagonal corners at (coord [4i], coord [4i+1]) and
 *  (coord [4i+2], coord [4i+3]) for  0 <= i < nrects,
 */
     Tk_Raster * raster;
     double* coord;
     int nrects;
{
   int i;
   int n = nrects*4;
   XRectangle *rect, *rectptr;
   int rx, ry, tmp;
   int minx = HIGH, miny = HIGH, maxx = LOW, maxy = LOW;

   if (nrects < 1) return;

   for (i = 0, rectptr= rect= (XRectangle*) malloc(sizeof (XRectangle)*nrects);
	i < n; i+=4, rectptr++) {
      WorldToRaster (raster, coord [i], coord [i+1], &rx, &ry);
      if (rx < minx) minx = rx;
      if (rx > maxx) maxx = rx;
      if (ry < miny) miny = ry;
      if (ry > maxy) maxy = ry;
      rectptr->x = rx;
      rectptr->y = ry;
      WorldToRaster (raster, coord [i+2], coord [i+3], &rx, &ry);
      if (rx < minx) minx = rx;
      if (rx > maxx) maxx = rx;
      if (ry < miny) miny = ry;
      if (ry > maxy) maxy = ry;
      if (rectptr->x > rx) { tmp = rectptr->x; rectptr->x = rx; rx = tmp; }
      if (rectptr->y > ry) { tmp = rectptr->y; rectptr->y = ry; ry = tmp; }
      rectptr->width = rx - rectptr->x;
      rectptr->height = ry - rectptr->y;
   }

   XFillRectangles (GetRasterDisplay (raster),
		    GetRasterDrawable (raster),
		    GetRasterGC (raster),
		    rect, nrects);
   free (rect);

   SetRasterModifiedArea (raster, minx, miny, maxx, maxy);
}

void RasterFillPolygon (raster, coord, npts)
/*   -----------------
 *
 *  Fills the polygon bounded by the line segments connecting the successive
 *  'npts' points (coord [2i], coord [2i+1]), for 0 <= i < npts.
 *  If the last point is not equal to the first point, a line segment
 *  connecting them is added.
 */
     Tk_Raster * raster;
     double* coord;
     int npts;
{
   int i;
   int n = npts*2;
   XPoint *pt, *ptptr;
   int rx, ry;
   int minx = HIGH, miny = HIGH, maxx = LOW, maxy = LOW;

   if (npts < 1) return;

   for (i = 0, ptptr = pt = (XPoint*) malloc (sizeof (XPoint)*npts);
	i < n; i+=2, ptptr++) {
      WorldToRaster (raster, coord [i], coord [i+1], &rx, &ry);
      if (rx < minx) minx = rx;
      if (rx > maxx) maxx = rx;
      if (ry < miny) miny = ry;
      if (ry > maxy) maxy = ry;
      ptptr->x = rx;
      ptptr->y = ry;
   }

   XFillPolygon (GetRasterDisplay (raster),
		 GetRasterDrawable (raster),
		 GetRasterGC (raster),
		 pt, npts, Complex, CoordModeOrigin);
   free (pt);

   SetRasterModifiedArea (raster, minx, miny, maxx, maxy);
}

/*
 * C interface to tcl routine envcreate
 */
int
CreateDrawEnviron(Tcl_Interp *interp,
		  Tk_Raster *RasterPtr,
		  int argc, char **argv)
{
    int result;

    CreateDrawEnv(interp, RasterPtr, argc, argv);
    result = atoi(Tcl_GetStringResult(interp));
    return result;
}

/*
 * C interface to tcl routine envset
 */
int
SetDrawEnviron (Tcl_Interp *interp,
		Tk_Raster *RasterPtr,
		int index)
{
    int result;
    DrawEnvironment *drawEnvPtr;

    if (DrawEnvIndex (interp, RasterPtr, index,
		      (ClientData*)&drawEnvPtr) != TCL_OK){
#ifdef DEBUG
	printf("ERROR1: failed DrawEnvIndex\n");
#endif
	return -1;
    }
    result = SetDrawEnv (interp, RasterPtr, drawEnvPtr);
    return result;
}

char *
GetRasterColour (Tcl_Interp *interp,
		 Tk_Raster *RasterPtr,
		 int index)
{
    int result;
    DrawEnvironment *drawEnvPtr;
    static char colour[10];

    if (DrawEnvIndex (interp, RasterPtr, index,
		      (ClientData*)&drawEnvPtr) != TCL_OK){
#ifdef DEBUG
	printf("ERROR2: failed DrawEnvIndex %d\n", index);
#endif
    }
    result = SetDrawEnv (interp, RasterPtr, drawEnvPtr);
    sprintf(colour, "#%02x%02x%02x", drawEnvPtr->fgColor->red/256,
	    drawEnvPtr->fgColor->green/256, drawEnvPtr->fgColor->blue/256);
    return colour;
}

int GetRasterLineWidth (Tcl_Interp *interp,
		       Tk_Raster *RasterPtr,
		       int index)
{
    int result;
    DrawEnvironment *drawEnvPtr;

    if (DrawEnvIndex (interp, RasterPtr, index,
		      (ClientData*)&drawEnvPtr) != TCL_OK){

    }
    result = SetDrawEnv (interp, RasterPtr, drawEnvPtr);

    return (drawEnvPtr->gcValues.line_width);
}


static void
RasterClear(Raster *RasterPtr)
{
    if (RasterPtr->plot_func) {
	RasterPtr->plot_func(RasterPtr, Tk_PathName(RasterPtr->tkwin),
			     RASTER_INIT, 0, 0, 0, 0);
    }

    /* set the crosshair flag to be "not visible" */
    Tcl_VarEval(RasterPtr->interp, "unset_raster_xh ",
		Tk_PathName(RasterPtr->tkwin), NULL);
    XFillRectangle (RasterPtr->display, RasterPtr->pm, RasterPtr->copyGC,
		       0, 0, RasterPtr->width, RasterPtr->height);
    arrangeDisplay (RasterPtr, LOW, LOW, HIGH, HIGH, RASTER_DRAW_PIXMAP);
}

void
tk_RasterClear(Tk_Raster *rasterptr)
{
    Raster* RasterPtr = (Raster*) rasterptr;

    RasterClear(RasterPtr);
}

void
RasterClearArea(Raster *RasterPtr,
		int x0,	int y0, int x1, int y1)
{

    XFillRectangle (RasterPtr->display, RasterPtr->pm, RasterPtr->copyGC,
		       x0, y0, x1, y1);
    arrangeDisplay (RasterPtr, LOW, LOW, HIGH, HIGH, RASTER_DRAW_PIXMAP);

}

void GetRasterCoords(Tk_Raster *rasterptr,
		     double *wx0,
		     double *wy0,
		     double *wx1,
		     double *wy1)
{
    Raster* raster = (Raster*) rasterptr;

    *wx0 = raster->wx0;
    *wy0 = raster->wy0;
    *wx1 = raster->wx1;
    *wy1 = raster->wy1;
}

static void
RasterScrollX(Raster *raster,
	      int new_x, /* pixel coords */
	      int old_x) /* pixel coords */
{
    double max_x, min_x;
    double wx0, wx1;
    double wnew_x, wold_x, wnew_y, wold_y;
    double dx;

    /* find the portion of the raster (left or right) to be replotted */
    RasterToWorld(raster, new_x-old_x, 0, &wnew_x, &wnew_y);
    RasterToWorld(raster, 0, 0, &wold_x, &wold_y);

    min_x = raster->wx0; /* min world base on raster */
    max_x = raster->wx1; /* max world base on raster */
    dx = wnew_x - wold_x; /* amount of scrolling */
    if (dx > 0) {
	/* fill right edge */
	wx1 = max_x + dx;
	if (dx > max_x - min_x)
	    wx0 = min_x + dx;
	else
	    wx0 = max_x;
    } else {
	/* fill left edge */
	wx0 = min_x + dx;
	if (-dx > max_x - min_x)
	    wx1 = max_x + dx;
	else
	    wx1 = min_x;
    }

#ifdef SCROLL
    printf("woldx %f wnewx %f dx %f wx0 %f wx1 %f max_x %f min_x %f\n", wold_x, wnew_x, dx, wx0, wx1, max_x, min_x);
    printf("raster->wx0=%f,dx=%f, +=%f\n",
	   raster->wx0, dx, raster->wx0+dx);
#endif
    /*
     * set up new world coords so that the replotting results appear in the
     * correct place
     */
    SetRasterCoords (raster, wnew_x, raster->wy0,
		     wnew_x + (raster->wx1 - raster->wx0), raster->wy1);

#ifdef SCROLL
/*    printf("SetRasterCoords wx0 %f wx1 %f\n", wx0, wx1);*/
#endif
    if (raster->plot_func) {
	raster->plot_func(raster, Tk_PathName(raster->tkwin),
			  RASTER_REPLOT_SLIVER,
			  (int)wx0, 0, (int)(wx1+1), 0);
    }
}

static void
RasterScrollY(Raster *raster,
	      int new_y, /* pixel coords */
	      int old_y) /* pixel coords */
{
    double max_y, min_y;
    double wy0, wy1;
    double wnew_x, wold_x, wnew_y, wold_y;
    double dy = 0;

    /* find the portion of the raster (top or bottom) to be replotted */
    RasterToWorld(raster, 0, new_y-old_y, &wnew_x, &wnew_y);
    RasterToWorld(raster, 0, 0, &wold_x, &wold_y);

    min_y = raster->wy0; /* min world base on raster */
    max_y = raster->wy1; /* max world base on raster */
    dy = wnew_y - wold_y; /* amount of scrolling */
    if (dy < 0) {
	/* fill top edge */
	wy0 = min_y;
	wy1 = min_y + dy;
    } else {
	/* fill bottom edge */
	wy0 = max_y;
	wy1 = max_y + dy;
    }

#ifdef SCROLL
    printf("woldy %f wnewy %f dy %f wy0 %f wy1 %f max_y %f min_y %f\n", wold_y, wnew_y, dy, wy0, wy1, max_y, min_y);
    printf("raster->wy0=%f,dy=%f, +=%f\n",
	   raster->wy0, dy, raster->wy0+dy);
#endif

    /*
     * set up new world coords so that the replotting results appear in the
     * correct place
     */
    SetRasterCoords (raster, raster->wx0, wnew_y,
		     raster->wx1, wnew_y + (raster->wy1 - raster->wy0));

#ifdef SCROLL
    printf("SetRasterCoords wy0 %f wy1 %f\n", wy0, wy1);
    printf("start %f end %f \n", raster->wy_start, raster->wy_end);

#endif
    /*
     * because we display the results "upside down", need to subtract the
     * max y in world coords. Need to add 1 because I plot to max y + 1
     */
    if (raster->plot_func) {
	raster->plot_func(raster, Tk_PathName(raster->tkwin),
			  RASTER_REPLOT_ALL,
			  (int)raster->wx0,
			  (int)(raster->wy_end - wy0)+1,
			  (int)raster->wx1, (int)(raster->wy_end - wy1) + 2);
    }
}

void RasterResetWorldScroll(Tk_Raster *rasterptr)
{
    Raster* raster = (Raster*) rasterptr;

#ifdef DEBUG
    printf("RasterResetWorldScroll\n");
#endif

    raster->wx_start = DBL_MAX/2;
    raster->wy_start = DBL_MAX/2;
    raster->wx_end = -DBL_MAX/2;
    raster->wy_end = -DBL_MAX/2;
}


/*
 * set the scroll region for a single raster
 * return 1 if the scroll region was changed, else return 0
 */
int
RasterSetWorldScroll(Tk_Raster *rasterptr,
		     double x0, double y0,
		     double x1, double y1)
{
    Raster* raster = (Raster*) rasterptr;
    int flag = 0;

#ifdef DEBUG
    printf("RasterSetWorldScroll %s x0 %f x1 %f wx0 %f wx1 %f\n",
	   Tk_PathName(raster->tkwin), x0, x1, raster->wx_start, raster->wx_end);
    printf("RasterSetWorldScroll y0 %f y1 %f wy0 %f wy1 %f\n",
	   y0, y1, raster->wy_start, raster->wy_end);
#endif

    if (x0 != raster->wx_start) {
	raster->wx_start = x0;
	flag = 1;
    }
    if (y0 != raster->wy_start) {
	raster->wy_start = y0;
	flag = 1;
    }
    if (x1 != raster->wx_end) {
	raster->wx_end = x1;
	flag = 1;
    }
    if (y1 != raster->wy_end) {
	raster->wy_end = y1;
	flag = 1;
    }

    /*
     * need some checks here to ensure that the start and ends are never
     * the identical for the same direction which could lead to a divide
     * by 0 error
     */

    if (raster->wx_start == raster->wx_end) {
	raster->wx_start -= 1e-10;
	raster->wx_end += 1e-10;
    }

    if (raster->wy_start == raster->wy_end) {
	raster->wy_start -= 1e-10;
	raster->wy_end += 1e-10;
    }

#ifdef DEBUG
    printf("RasterSetWorldScroll %s x0 %f x1 %f y0 %f y1 %f \n",
	   Tk_PathName(raster->tkwin), raster->wx_start, raster->wx_end,
	   raster->wy_start, raster->wy_end);
#endif

    return flag;
}

void
RasterGetWorldScroll(Tk_Raster *rasterptr,
		     double *x0, double *y0,
		     double *x1, double *y1)
{
    Raster* raster = (Raster*) rasterptr;

    *x0 = raster->wx_start;
    *y0 = raster->wy_start;
    *x1 = raster->wx_end;
    *y1 = raster->wy_end;

#ifdef DEBUG
    printf("RasterGetWorldScroll %s %f %f %f %f \n",
	   Tk_PathName(raster->tkwin), *x0, *y0, *x1, *y1);
#endif
}

void
RasterWinSize(Tk_Raster *rasterptr,
	      int *width, int *height)
{
    Raster* raster = (Raster*) rasterptr;

    *width = raster->width;
    *height = raster->height;
}

void
RasterInitPlotFunc(Tk_Raster *rasterptr,
		   void (*plot_func)(Tk_Raster *raster, char *raster_win,
				     int job, int x0, int y0, int x1,
				     int y1))
{

    Raster* raster = (Raster*) rasterptr;
    raster->plot_func = plot_func;
}

void RasterDeletePlotFunc(Tk_Raster *rasterptr)
{

    Raster* raster = (Raster*) rasterptr;
    raster->plot_func = NULL;
}

/*
 * find the new_range given an old_range and a mag value
 * from a scale range of say, 1 to 100, a value of 1 gives a mag of 100/100,
 * a value of 2 gives 99/100 and a value of 100 depicts 10 bases
 * using: scale_min   old_range
 *        scale_max   min_bases
 *  scale_min = m * old_range + c
 *  scale_max = m * min_bases + c
 */
double RasterSetXMag(Tcl_Interp *interp,
		     double old_range,
		     int value)
{
    double m, c;
    double new_range;
    int scale_min, scale_max;
    int min_bases;

    scale_min = get_default_int(interp, tk_utils_defs, "RASTER.SCALEX.MIN");
    scale_max = get_default_int(interp, tk_utils_defs, "RASTER.SCALEX.MAX");
    min_bases = get_default_int(interp, tk_utils_defs, "RASTER.SCALEX.MIN_BASES");

    m = (old_range - min_bases) / (scale_min - scale_max);
    c = min_bases - (m * scale_max);
    new_range = m * value + c;

#ifdef DEBUG
    printf("RasterSetXMag value %d old %f min_bases %d m %f c %f new %f\n",
	   value, old_range, min_bases, m, c, new_range);
#endif
    return new_range;
}

/*
 * find the mag value to convert old_range into new_range
 */
double RasterGetXMag(Tcl_Interp *interp,
		     double old_range,
		     double new_range)
{
    int scale_min, scale_max;
    double value;
    double m, c;
    int min_bases;

    scale_min = get_default_int(interp, tk_utils_defs, "RASTER.SCALEX.MIN");
    scale_max = get_default_int(interp, tk_utils_defs, "RASTER.SCALEX.MAX");
    min_bases = get_default_int(interp, tk_utils_defs, "RASTER.SCALEX.MIN_BASES");

    m = (old_range - min_bases) / (scale_min - scale_max);
    c = min_bases - (m * scale_max);
    value = (new_range - c) / m;
    return value;
}

double RasterSetYMag(Tcl_Interp *interp,
		     double old_range,
		     int value)
{
    double new_range;
    int scale_min, scale_max;
    int scale_range;

    scale_min = get_default_int(interp, tk_utils_defs, "RASTER.SCALEY.MIN");
    scale_max = get_default_int(interp, tk_utils_defs, "RASTER.SCALEY.MAX");

    scale_range = scale_max - scale_min + 1;
    /* percentage decrement */
    new_range = (double)(scale_range - value + 1) / scale_range * old_range;
#ifdef DEBUG
    printf("RasterSetYMag: scale %d value %d old %f new %f\n", scale_range, value,
	   old_range, new_range);
#endif
    return new_range;
}

double RasterGetYMag(Tcl_Interp *interp,
		     double old_range,
		     double new_range)
{
    double value;
    int scale_min, scale_max;
    int scale_range;

    scale_min = get_default_int(interp, tk_utils_defs, "RASTER.SCALEY.MIN");
    scale_max = get_default_int(interp, tk_utils_defs, "RASTER.SCALEY.MAX");
    scale_range = scale_max - scale_min + 1;

    value = scale_range - (new_range * scale_range / old_range) + 1;

#ifdef DEBUG
    printf("RasterGetYMag old %f new %f value %f\n", old_range, new_range, value);
#endif
    return value;
}

void RasterCallPlotFunc(Tk_Raster *rasterptr,
			int job,
			int x0, int y0, int x1, int y1)
{
    Raster* raster = (Raster*) rasterptr;

    if (raster->plot_func) {
	raster->plot_func(raster, Tk_PathName(raster->tkwin), job, x0, y0, x1, y1);
    }
}
