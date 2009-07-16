#include <stdio.h>
#include <float.h>
#include "os.h"
#include "tkInt.h"
#include "tkPort.h"
#include "tkCanvas.h"
#include "tkCanvGraph.h"
#include "tclCanvGraph.h"
#include "matrix.h"

#ifdef __MINGW32__
extern int		Tk_CanvasTagsParseProc _ANSI_ARGS_((
				ClientData clientData, Tcl_Interp * interp, 
				Tk_Window tkwin, CONST char * value, 
				char * widgRec, int offset));
extern char *		Tk_CanvasTagsPrintProc _ANSI_ARGS_((
				ClientData clientData, Tk_Window tkwin, 
				char * widgRec, int offset, 
				Tcl_FreeProc ** freeProcPtr));
#endif

/* HACK - put somewhere else */
#define ROUND(x)   ((x) < 0.0 ? ceil((x) - 0.5) : floor((x) + 0.5))

/*
 * The structure below defines the record for each graph item.
 */

typedef struct GraphItem  {
    Tk_Item header;		/* Generic stuff that's the same for all
				 * types.  MUST BE FIRST IN STRUCTURE. */
    Tk_Anchor anchor;		/* Where to anchor bitmap relative to
				 * (x,y). */
    double an_x, an_y;		/* Coordinates of positioning point for
				 * graph. */
    int x, y;		        /* Position of viewable pixmap's upper-left
				 * corner within the widget's area */
    int new_x, new_y;           /* how much scrolling has occurred */

    Graph *graph;		/* string format graph to display in window. */
    double min_x, max_x, min_y, max_y; /* min and max x and y pixel values */
    GC copyGC;        	        /* Graphics context for copying to screen */
    GC paintGC;                 /* Graphics context for painting to canvas */
    GC drawGC;			/* Graphics context to use for drawing
				 * graph on screen. */
    Pixmap pm;			/* Pixmap for storing the raster */

    int x0, y0, x1, y1;		/* coordinates of rectangle of the
				   pixmap that has to be redisplayed */
    double line_width;          /* width to use for graph */
    int cap_style;		/* Cap style for line. */
    int join_style;		/* Join style for line. */
    int line_style;
    int fill_style;
    XColor *fgColor;		/* Foreground color to use for graph. */
    Tk_3DBorder bgBorder;	/* Background color to use for graph. */

    int invertx;                /* non-zero means reverse x values */
    int inverty;                /* non-zero means reverse y values */

    Tcl_Interp *interp;
    int redraw;                 /* whether to redraw the pixmap data */
    double originX, originY;		/* Origin about which to scale rect. */
    double M[3][3];             /* scaling matrix */
    int prev_width;             /* store previous width/height to determine */
    int prev_height;            /* window has resized */
    int vertical;               /* orientation of window */
} GraphItem;

/*
 * Prototypes for procedures defined in this file:
 */

static int		GraphCoords _ANSI_ARGS_((Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr, int argc,
			    Tcl_Obj *CONST argv[]));
static int		GraphToArea _ANSI_ARGS_((Tk_Canvas canvas,
			    Tk_Item *itemPtr, double *rectPtr));
static double		GraphToPoint _ANSI_ARGS_((Tk_Canvas canvas,
			    Tk_Item *itemPtr, double *coordPtr));
static int		GraphToPostscript _ANSI_ARGS_((Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr, int prepass));
static void		ComputeGraphBbox _ANSI_ARGS_((Tk_Canvas canvas,
			    GraphItem *graphPtr));
static int		ConfigureGraph _ANSI_ARGS_((Tcl_Interp *interp,
			    Tk_Canvas canvas, Tk_Item *itemPtr, int argc,
			    Tcl_Obj *CONST argv[], int flags));
static int		CreateGraph _ANSI_ARGS_((Tcl_Interp *interp,
			    Tk_Canvas canvas, struct Tk_Item *itemPtr,
			    int argc, Tcl_Obj *CONST argv[]));
static void		DeleteGraph _ANSI_ARGS_((Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display));
static void		DisplayGraph _ANSI_ARGS_((Tk_Canvas canvas,
			    Tk_Item *itemPtr, Display *display, Drawable dst,
			    int x, int y, int width, int height));
static void		ScaleGraph _ANSI_ARGS_((Tk_Canvas canvas,
			    Tk_Item *itemPtr, double originX, double originY,
			    double scaleX, double scaleY));
static void		TranslateGraph _ANSI_ARGS_((Tk_Canvas canvas,
			    Tk_Item *itemPtr, double deltaX, double deltaY));

int pixel_to_canvas(Tcl_Interp *interp,
		    Tk_Window tkwin,
		    char direction,
		    int pixel);

void GraphPlotFunc(Tk_Canvas canvas,
		   Display *display,
		   GraphItem *graphPtr,
		   int x1,
		   int x2);

void ResizeGraph(Tk_Canvas canvas,
		 Display *display,
		 Tk_Window tkwin,
		 GraphItem *graphPtr,
		 int old_width,
		 int old_height,
		 int new_width,
		 int new_height);

void canvas_array_line(GraphItem *graphPtr,
		       gd_line array,
		       g_pt *start,
		       g_pt *end);

void canvas_array_pt(GraphItem *graphPtr,
		     g_pt array,
		     g_pt *pt);

void find_nearest_line(GraphItem *graphPtr,
		       gd_line *array,
		       int num_lines,
		       int index,
		       g_pt pt,
		       int start_x,
		       int counter,
		       double *min);

void find_nearest_cline(GraphItem *graphPtr,
			g_pt *array,
			int n_pts,
			int index,
			g_pt pt,
			int start_x,
			int counter,
			double *min);

void find_nearest_dot(GraphItem *graphPtr,
		      g_pt *array,
		      int n_pts,
		      int index,
		      g_pt pt,
		      int start_x,
		      int counter,
		      double *min);

int Tk_CanvasGraphParseProc(ClientData clientData,
                            Tcl_Interp *interp,
                            Tk_Window tkwin,
                            CONST char *value,
                            char *widgRec,
                            int offset);
char *
Tk_CanvasGraphPrintProc(ClientData clientData,
                       Tk_Window tkwin,
                       char *widgRec,
                       int offset,
                        Tcl_FreeProc **freeProcPtr);

/* From tkUtil.c */
static char *
TkOrientPrintProc_(
		  ClientData clientData,      /* Ignored. */
		  Tk_Window tkwin,            /* Window containing canvas widget. */
		  char *widgRec,              /* Pointer to record for item. */
		  int offset,                 /* Offset into item. */
    Tcl_FreeProc **freeProcPtr) /* Pointer to variable to fill in with
                                 * information about how to reclaim storage
                                 * for return string. */
{
    register int *statePtr = (int *) (widgRec + offset);

    if (*statePtr) {
        return "vertical";
    } else {
        return "horizontal";
    }
}

/* From tkUtil.c */
static char *
TkPixelPrintProc_(
		 ClientData clientData,      /* not used */
		 Tk_Window tkwin,            /* not used */
		 char *widgRec,              /* Widget structure record */
		 int offset,                 /* Offset of tile in record */
		 Tcl_FreeProc **freeProcPtr) /* not used */
{
    double *doublePtr = (double *) (widgRec + offset);
    char *p = (char *) ckalloc(24);

    Tcl_PrintDouble(NULL, *doublePtr, p);
    *freeProcPtr = TCL_DYNAMIC;
    return p;
}

static int
TkOrientParseProc_(
		  ClientData clientData,      /* some flags.*/
		  Tcl_Interp *interp,         /* Used for reporting errors. */
		  Tk_Window tkwin,            /* Window containing canvas widget. */
		  const char *value,          /* Value of option. */
		  char *widgRec,              /* Pointer to record for item. */
		  int offset)                 /* Offset into item. */
{
    int c;
    size_t length;

    register int *orientPtr = (int *) (widgRec + offset);

    if(value == NULL || *value == 0) {
        *orientPtr = 0;
        return TCL_OK;
    }

    c = value[0];
    length = strlen(value);

    if ((c == 'h') && (strncmp(value, "horizontal", length) == 0)) {
        *orientPtr = 0;
        return TCL_OK;
    }
    if ((c == 'v') && (strncmp(value, "vertical", length) == 0)) {
        *orientPtr = 1;
        return TCL_OK;
    }
    Tcl_AppendResult(interp, "bad orientation \"", value,
            "\": must be vertical or horizontal", NULL);
    *orientPtr = 0;
    return TCL_ERROR;
}

static int
TkGetDoublePixels_(
		  Tcl_Interp *interp,         /* Use this for error reporting. */
		  Tk_Window tkwin,            /* Window whose screen determines conversion
					       * from centimeters and other absolute
					       * units. */
		  CONST char *string,         /* String describing a number of pixels. */
		  double *doublePtr)          /* Place to store converted result. */
{
    char *end;
    double d;

    d = strtod((char *) string, &end);
    if (end == string) {
    error:
        Tcl_AppendResult(interp, "bad screen distance \"", string, "\"", NULL);
        return TCL_ERROR;
    }
    while ((*end != '\0') && isspace(UCHAR(*end))) {
        end++;
    }
    switch (*end) {
    case 0:
        break;
    case 'c':
        d *= 10*WidthOfScreen(Tk_Screen(tkwin));
        d /= WidthMMOfScreen(Tk_Screen(tkwin));
        end++;
        break;
    case 'i':
        d *= 25.4*WidthOfScreen(Tk_Screen(tkwin));
        d /= WidthMMOfScreen(Tk_Screen(tkwin));
        end++;
        break;
    case 'm':
        d *= WidthOfScreen(Tk_Screen(tkwin));
        d /= WidthMMOfScreen(Tk_Screen(tkwin));
        end++;
        break;
    case 'p':
        d *= (25.4/72.0)*WidthOfScreen(Tk_Screen(tkwin));
        d /= WidthMMOfScreen(Tk_Screen(tkwin));
        end++;
        break;
    default:
        goto error;
    }
    while ((*end != '\0') && isspace(UCHAR(*end))) {
        end++;
    }
    if (*end != 0) {
        goto error;
    }
    *doublePtr = d;
    return TCL_OK;
}

static int
TkPixelParseProc_(
    ClientData clientData,      /* If non-NULL, negative values are allowed as
                                 * well */
    Tcl_Interp *interp,         /* Interpreter to send results back to */
    Tk_Window tkwin,            /* Window on same display as tile */
    const char *value,          /* Name of image */
    char *widgRec,              /* Widget structure record */
    int offset)                 /* Offset of tile in record */
{
    double *doublePtr = (double *) (widgRec + offset);
    int result;

    result = TkGetDoublePixels_(interp, tkwin, value, doublePtr);

    if ((result == TCL_OK) && (clientData == NULL) && (*doublePtr < 0.0)) {
        Tcl_AppendResult(interp, "bad screen distance \"", value, "\"", NULL);
        return TCL_ERROR;
    }
    return result;
}

static Tk_CustomOption pixelOption = {
    (Tk_OptionParseProc *) TkPixelParseProc_,
    TkPixelPrintProc_, (ClientData) NULL
};

static Tk_CustomOption graphOption = {
    (Tk_OptionParseProc *) Tk_CanvasGraphParseProc,
    Tk_CanvasGraphPrintProc, (ClientData) NULL
};

static Tk_CustomOption tagsOption = {
    (Tk_OptionParseProc *) Tk_CanvasTagsParseProc,
    Tk_CanvasTagsPrintProc, (ClientData) NULL
};

static Tk_CustomOption orientOption = {
    (Tk_OptionParseProc *) TkOrientParseProc_,
    TkOrientPrintProc_,
    (ClientData) NULL
};

static Tk_ConfigSpec configSpecs[] = {
    {TK_CONFIG_ANCHOR, "-anchor", (char *) NULL, (char *) NULL,
        "center", Tk_Offset(GraphItem, anchor), TK_CONFIG_DONT_SET_DEFAULT},

    {TK_CONFIG_CUSTOM, "-graph", (char *) NULL, (char *) NULL,
        (char *) NULL, Tk_Offset(GraphItem, graph), TK_CONFIG_NULL_OK,
        &graphOption},

    {TK_CONFIG_COLOR, "-fill", (char *) NULL, (char *) NULL,
	"black", Tk_Offset(GraphItem, fgColor), TK_CONFIG_NULL_OK},

    /*
    {TK_CONFIG_COLOR, "-foreground", "foreground", "Foreground",
        "black", Tk_Offset(GraphItem, fgColor), TK_CONFIG_NULL_OK},
    {TK_CONFIG_SYNONYM, "-fg", "foreground", (char *) NULL,
        (char *) NULL, 0, 0},
    */
    {TK_CONFIG_CUSTOM, "-width", (char *) NULL, (char *) NULL,
        "1.0", Tk_Offset(GraphItem, line_width),
        TK_CONFIG_DONT_SET_DEFAULT, &pixelOption},

    {TK_CONFIG_CUSTOM, "-tags", (char *) NULL, (char *) NULL,
	(char *) NULL, 0, TK_CONFIG_NULL_OK, &tagsOption},

    {TK_CONFIG_BOOLEAN, "-invertx", (char *) NULL, (char *) NULL,
	"", Tk_Offset(GraphItem, invertx), 0},

    {TK_CONFIG_BOOLEAN, "-inverty", (char *) NULL, (char *) NULL,
	"", Tk_Offset(GraphItem, inverty), 0},

    {TK_CONFIG_CUSTOM, "-orient", (char *) NULL, (char *) NULL,
	"horizontal", Tk_Offset(GraphItem, vertical), TK_CONFIG_DONT_SET_DEFAULT,
	&orientOption},

    {TK_CONFIG_END, (char *) NULL, (char *) NULL, (char *) NULL,
        (char *) NULL, 0, 0}
};


/*
 * The structures below defines the graph item type in terms of
 * procedures that can be invoked by generic item code.
 */
Tk_ItemType tkGraphType = {
    "graph",                            /* name */
    sizeof(GraphItem),                  /* itemSize */
    CreateGraph,                        /* createProc */
    configSpecs,                        /* configSpecs */
    ConfigureGraph,                     /* configureProc */
    GraphCoords,                        /* coordProc */
    DeleteGraph,			/* deleteProc */
    DisplayGraph,			/* displayProc */
    TK_CONFIG_OBJS,			/* flags */
    GraphToPoint,			/* pointProc */
    GraphToArea,			/* areaProc */
    GraphToPostscript,			/* postscriptProc */
    ScaleGraph,			        /* scaleProc */
    TranslateGraph,			/* translateProc */
    (Tk_ItemIndexProc *) NULL,		/* indexProc */
    (Tk_ItemCursorProc *) NULL,		/* icursorProc */
    (Tk_ItemSelectionProc *) NULL,	/* selectionProc */
    (Tk_ItemInsertProc *) NULL,		/* insertProc */
    (Tk_ItemDCharsProc *) NULL,		/* dTextProc */
    (Tk_ItemType *) NULL,		/* nextPtr */
};

#define MAX_STATIC_POINTS 200

int
Tk_CanvasGraphParseProc(ClientData clientData,
                        Tcl_Interp *interp,
                        Tk_Window tkwin,
                        CONST char *value,
                        char *widgRec,
                        int offset)
{
    /*
    Graph *graph = (Graph *)(widgRec + offset);
    printf("n_pts %d\n", graph->n_pts);
    */
#ifdef DEBUG
    printf("Tk_CanvasGraphParseProc\n");
#endif
    return TCL_OK;
}

char *
Tk_CanvasGraphPrintProc(ClientData clientData,
                       Tk_Window tkwin,
                       char *widgRec,
                       int offset,
                       Tcl_FreeProc **freeProcPtr)
{

    return NULL;
}

void InitGraphDim(Tcl_Interp *interp,
		  Tk_Canvas canvas,
		  Tk_Item *itemPtr)
{
    GraphItem *graphPtr = (GraphItem *) itemPtr;
    Graph *graph;
    m_coords new_pt, old_pt;

    graph = graphPtr->graph;

#ifdef DEBUG
    printf("DIMENSIONS %f %f %f %f\n", graph->dim.x0, graph->dim.y0,
	   graph->dim.x1, graph->dim.y1);
    printf("MIN0 %f %f MAX %f %f\n", graphPtr->min_x, graphPtr->min_y,
	   graphPtr->max_x, graphPtr->max_y);
#endif

    if (graphPtr->vertical) {
	make_coord(&old_pt, graph->dim.y0, graph->dim.x0);
    } else {
	make_coord(&old_pt, graph->dim.x0, graph->dim.y0);
    }
    multiply_coords_by_matrix(&new_pt, graphPtr->M, &old_pt);

    if (new_pt.x < graphPtr->min_x)
	graphPtr->min_x = new_pt.x;
    if (new_pt.y < graphPtr->min_y)
	graphPtr->min_y = new_pt.y;

    if (graphPtr->vertical) {
	make_coord(&old_pt, graph->dim.y1, graph->dim.x1);
    } else {
	make_coord(&old_pt, graph->dim.x1, graph->dim.y1);
    }
    multiply_coords_by_matrix(&new_pt, graphPtr->M, &old_pt);

    if (new_pt.x > graphPtr->max_x)
	graphPtr->max_x = new_pt.x;

    if (new_pt.y > graphPtr->max_y)
	graphPtr->max_y = new_pt.y;

#ifdef DEBUG
    printf("MIN1 %f %f MAX %f %f\n", graphPtr->min_x, graphPtr->min_y,
	   graphPtr->max_x, graphPtr->max_y);
#endif
    ComputeGraphBbox(canvas, graphPtr);
}



/*
 *--------------------------------------------------------------
 *
 * CreateGraph --
 *
 *	This procedure is invoked to create a new graph
 *	item in a canvas.
 *
 * Results:
 *	A standard Tcl return value.  If an error occurred in
 *	creating the item, then an error message is left in
 *	the interp's result;  in this case itemPtr is left uninitialized,
 *	so it can be safely freed by the caller.
 *
 * Side effects:
 *	A new graph item is created.
 *
 *--------------------------------------------------------------
 */

static int
CreateGraph(Tcl_Interp *interp,
	    Tk_Canvas canvas,
	    Tk_Item *itemPtr,
	    int argc,
	    Tcl_Obj *CONST *argv)
{
    GraphItem *graphPtr = (GraphItem *) itemPtr;
    int i;

#ifdef DEBUG
    printf("&&&&&&&&&&&&&&&&&&&&CreateGraph %p\n", graphPtr);
#endif
    if (argc==1) {
	i = 1;
    } else {
	char *arg = Tcl_GetStringFromObj(argv[1], NULL);
	if (((argc>1) && (arg[0] == '-')
		&& (arg[1] >= 'a') && (arg[1] <= 'z'))) {
	    i = 1;
	} else {
	    i = 2;
	}
    }

    if (argc < i) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
		Tk_PathName(Tk_CanvasTkwin(canvas)), " create ",
		itemPtr->typePtr->name, " x y ?options?\"",
		(char *) NULL);
	return TCL_ERROR;
    }

    /*
     * Initialize item's record.
     */

    graphPtr->anchor = TK_ANCHOR_CENTER;
    graphPtr->fgColor = NULL;
    graphPtr->copyGC = None;
    graphPtr->drawGC = None;
    graphPtr->paintGC = None;
    graphPtr->line_width = 1;
    graphPtr->line_style = LineSolid;
    graphPtr->cap_style = CapButt;
    graphPtr->fill_style = FillSolid;
    graphPtr->join_style = JoinRound;
    graphPtr->graph = NULL;
    graphPtr->invertx = 0;
    graphPtr->inverty = 0;
    graphPtr->min_x = DBL_MAX/2;
    graphPtr->min_y = DBL_MAX/2;
    graphPtr->max_x = -DBL_MAX/2;
    graphPtr->max_y = -DBL_MAX/2;
    graphPtr->pm = None;
    graphPtr->x = 0;
    graphPtr->y = 0;
    graphPtr->new_x = 0;
    graphPtr->new_y = 0;
    graphPtr->x0 = INT_MAX;
    graphPtr->y0 = INT_MAX;
    graphPtr->x1 = INT_MIN;
    graphPtr->y1 = INT_MIN;
    graphPtr->interp = NULL;
    graphPtr->originX = 0.0;
    graphPtr->originY = 0.0;
    graphPtr->prev_width = 0;
    graphPtr->prev_height = 0;
    graphPtr->vertical = 0;

    make_identity_matrix(graphPtr->M);
    /*
     * Process the arguments to fill in the item record.
     */

    if ((GraphCoords(interp, canvas, itemPtr, i, argv) != TCL_OK)) {
	goto error;
    }

    if (ConfigureGraph(interp, canvas, itemPtr, argc-i, argv+i, 0) == TCL_OK) {
	InitGraphDim(interp, canvas, itemPtr);
	return TCL_OK;
    }


    error:
    DeleteGraph(canvas, itemPtr, Tk_Display(Tk_CanvasTkwin(canvas)));
    return TCL_ERROR;
}

/*
 *--------------------------------------------------------------
 *
 * GraphCoords --
 *
 *	This procedure is invoked to process the "coords" widget
 *	command on graph items.  See the user documentation for
 *	details on what it does.
 *
 * Results:
 *	Returns TCL_OK or TCL_ERROR, and sets the interp's result.
 *
 * Side effects:
 *	The coordinates for the given item may be changed.
 *
 *--------------------------------------------------------------
 */

static int
GraphCoords(interp, canvas, itemPtr, argc, argv)
     Tcl_Interp *interp;			/* Used for error reporting. */
     Tk_Canvas canvas;			/* Canvas containing item. */
     Tk_Item *itemPtr;			/* Item whose coordinates are to be
					 * read or modified. */
     int argc;				/* Number of coordinates supplied in
					 * argv. */
     Tcl_Obj *CONST argv[];		/* Array of coordinates: x1, y1,
					 * x2, y2, ... */
{
    GraphItem *graphPtr = (GraphItem *) itemPtr;
    Graph *graph;

    graph = graphPtr->graph;

#ifdef DEBUG
    printf("\n\n\nGraphCoords %d\n\n\n\n", argc);
#endif
    if (argc == 0) {
	Tcl_Obj *obj = Tcl_NewObj();
	Tcl_Obj *subobj = Tcl_NewDoubleObj(graph->dim.x0);
	Tcl_ListObjAppendElement(interp, obj, subobj);
	subobj = Tcl_NewDoubleObj(graph->dim.y0);
	Tcl_ListObjAppendElement(interp, obj, subobj);
	subobj = Tcl_NewDoubleObj(graph->dim.x1);
	Tcl_ListObjAppendElement(interp, obj, subobj);
	subobj = Tcl_NewDoubleObj(graph->dim.y1);
	Tcl_ListObjAppendElement(interp, obj, subobj);
	Tcl_SetObjResult(interp, obj);
    } else if (argc <3) {
	if (argc==1) {
	    if (Tcl_ListObjGetElements(interp, argv[0], &argc,
		    (Tcl_Obj ***) &argv) != TCL_OK) {
		return TCL_ERROR;
	    } else if (argc != 2) {
		char buf[64 + TCL_INTEGER_SPACE];

		sprintf(buf, "wrong # coordinates: expected 2, got %d", argc);
		Tcl_SetResult(interp, buf, TCL_VOLATILE);
		return TCL_ERROR;
	    }
	}
	if ((Tk_CanvasGetCoordFromObj(interp, canvas, argv[0],
				      &graphPtr->an_x) != TCL_OK)
	    || (Tk_CanvasGetCoordFromObj(interp, canvas, argv[1],
					 &graphPtr->an_y)
 		    != TCL_OK)) {
	    return TCL_ERROR;
	}
#ifdef DEBUG
	printf("ANCHOR %f %f\n", graphPtr->an_x, graphPtr->an_y);
#endif
	ComputeGraphBbox(canvas, graphPtr);
    } else {
	char buf[64 + TCL_INTEGER_SPACE];

	sprintf(buf, "wrong # coordinates: expected 0 or 2, got %d", argc);
	Tcl_SetResult(interp, buf, TCL_VOLATILE);
	return TCL_ERROR;
    }
    return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * ConfigureGraph --
 *
 *	This procedure is invoked to configure various aspects
 *	of a graph item, such as its anchor position.
 *
 * Results:
 *	A standard Tcl result code.  If an error occurs, then
 *	an error message is left in the interp's result.
 *
 * Side effects:
 *	Configuration information may be set for itemPtr.
 *
 *--------------------------------------------------------------
 */

static int
ConfigureGraph(Tcl_Interp *interp,
	       Tk_Canvas canvas,
	       Tk_Item *itemPtr,
	       int argc,
	       Tcl_Obj *CONST *argv,
	       int flags)
{
    GraphItem *graphPtr = (GraphItem *) itemPtr;
    XGCValues gcValues;
    Tk_Window tkwin;
    unsigned long mask;
    int i;
    int width, height;
    int clear_pm = 0;

#ifdef DEBUG
    printf("ConfigureGraph %d\n", argc);
#endif

    tkwin = Tk_CanvasTkwin(canvas);

    {
	Tcl_Obj *obj;

	for (i = 0; i < argc; i++) {
	    obj = (Tcl_Obj *)argv[i];
#ifdef DEBUG
	    printf("tcl_obj refcount %d length %d ",
		   obj->refCount, obj->length);
	    if (obj->bytes)
		printf("bytes %s\n", obj->bytes);
	    if (obj->typePtr)
		printf("name %s\n", obj->typePtr->name);
	    else
		printf("\n");
#endif
	}
    }

    if (Tk_ConfigureWidget(interp, tkwin, configSpecs, argc, (char **) argv,
	    (char *) graphPtr, flags|TK_CONFIG_OBJS) != TCL_OK) {

	printf("ERROR %s\n", Tcl_GetStringResult(interp));
	return TCL_ERROR;
    }

    /* MEGA HACK to deal with graph object directly, not via string conversion */
    for (i = 0; i < argc; i++) {
	Tcl_Obj *obj = (Tcl_Obj *)argv[i];
	if (obj->typePtr && (strcmp(obj->typePtr->name, "graph") == 0)) {
		graphPtr->graph = Tcl_GetGraphFromObj(obj);
	}
    }

    /*
     * A few of the options require additional processing, such as those
     * that determine the graphics context.
     */

    width = Tk_Width(tkwin);
    height = Tk_Height(tkwin);
#ifdef DEBUG
    printf("width %d height %d\n", width, height);
#endif

    if (graphPtr->pm == None) {
	graphPtr->pm = Tk_GetPixmap (Tk_Display(tkwin),
				     Tk_WindowId(tkwin),
				     width, height,
				     1);
	clear_pm = 1; /* flag since need create pixmap in order to create gc */
    }

    if (graphPtr->copyGC == None) {
	mask = GCFunction|GCGraphicsExposures|GCForeground;
	gcValues.function = GXcopy;
	gcValues.graphics_exposures = False;
	gcValues.foreground = 0;

	graphPtr->copyGC = XCreateGC(Tk_Display(tkwin), graphPtr->pm, mask, &gcValues);
    }

    if (graphPtr->drawGC == None) {
	mask = GCForeground|GCCapStyle|GCJoinStyle|GCLineWidth|GCLineStyle|GCFillStyle|GCFunction|GCGraphicsExposures;

	gcValues.foreground = 1;
	gcValues.cap_style = graphPtr->cap_style;
	gcValues.join_style = graphPtr->join_style;
	gcValues.line_width = (int) (graphPtr->line_width + 0.5);
	gcValues.line_style = graphPtr->line_style;
	gcValues.fill_style = graphPtr->fill_style;
	gcValues.function = GXcopy;
	gcValues.graphics_exposures = False;
	graphPtr->drawGC = XCreateGC(Tk_Display(tkwin), graphPtr->pm, mask, &gcValues);
    } else {
	mask = GCCapStyle|GCJoinStyle|GCLineWidth|GCLineStyle|GCFillStyle;

	gcValues.cap_style = graphPtr->cap_style;
	gcValues.join_style = graphPtr->join_style;
	gcValues.line_width = (int) (graphPtr->line_width + 0.5);
	gcValues.line_style = graphPtr->line_style;
	gcValues.fill_style = graphPtr->fill_style;
	XChangeGC(Tk_Display(tkwin), graphPtr->drawGC, mask, &gcValues);
    }

    mask = GCFunction|GCGraphicsExposures|GCForeground;
    gcValues.function = GXcopy;
    gcValues.graphics_exposures = False;
    gcValues.foreground = graphPtr->fgColor->pixel;

    if (graphPtr->paintGC == None) {
	graphPtr->paintGC = Tk_GetGC(tkwin, mask, &gcValues);
    } else {
	Tk_FreeGC(Tk_Display(tkwin), graphPtr->paintGC);
	graphPtr->paintGC = Tk_GetGC(tkwin, mask, &gcValues);
    }

    if (clear_pm) {
	XFillRectangle (Tk_Display(tkwin), graphPtr->pm, graphPtr->copyGC,
			0, 0, width, height);
    }

    graphPtr->redraw = 1;

    graphPtr->interp = interp;

    /* set scroll extents */
    graphPtr->x0 = pixel_to_canvas(graphPtr->interp, tkwin, 'x', 0);
    graphPtr->x1 = pixel_to_canvas(graphPtr->interp, tkwin, 'x', Tk_Width(tkwin));

    graphPtr->y0 = pixel_to_canvas(graphPtr->interp, tkwin, 'y', 0);
    graphPtr->y1 = pixel_to_canvas(graphPtr->interp, tkwin, 'y', Tk_Height(tkwin));

    ComputeGraphBbox(canvas, graphPtr);
    return TCL_OK;
}

/*
 *--------------------------------------------------------------
 *
 * DeleteGraph --
 *
 *	This procedure is called to clean up the data structure
 *	associated with a graph item.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Resources associated with itemPtr are released.
 *
 *--------------------------------------------------------------
 */

static void
DeleteGraph(Tk_Canvas canvas,
	    Tk_Item *itemPtr,
	    Display *display)
{
    GraphItem *graphPtr = (GraphItem *) itemPtr;

#ifdef DEBUG
    printf("DeleteGraph\n");
#endif

#ifdef TODO
    if (graphPtr->graph != None) {
	Tk_FreeBitmap(display, graphPtr->graph);
    }
#endif

    if (graphPtr->fgColor != NULL) {
	Tk_FreeColor(graphPtr->fgColor);
    }
    if (graphPtr->copyGC != None) {
	XFreeGC(display, graphPtr->copyGC);
    }
    if (graphPtr->drawGC != None) {
	XFreeGC(display, graphPtr->drawGC);
    }

    if (graphPtr->paintGC != None) {
	Tk_FreeGC(display, graphPtr->paintGC);
    }
}

/*
 *--------------------------------------------------------------
 *
 * ComputeGraphBbox --
 *
 *	This procedure is invoked to compute the bounding box of
 *	all the pixels that may be drawn as part of a graph item.
 *	This procedure is where the child graph's placement is
 *	computed.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The fields x1, y1, x2, and y2 are updated in the header
 *	for itemPtr.
 *
 *--------------------------------------------------------------
 */

	/* ARGSUSED */
static void
ComputeGraphBbox(Tk_Canvas canvas, GraphItem *graphPtr)
{
    Tk_Window tkwin;
    int x, y;
    int width, height;

    tkwin = Tk_CanvasTkwin(canvas);

    /*
     * Store the information in the item header.
     */

    x = (int) (graphPtr->an_x + ((graphPtr->an_x >= 0) ? 0.5 : - 0.5));
    y = (int) (graphPtr->an_y + ((graphPtr->an_y >= 0) ? 0.5 : - 0.5));

    /*
     * Compute location and size of graph, using anchor information.
     */
    width = graphPtr->max_x - graphPtr->min_x + 1;
    height = graphPtr->max_y - graphPtr->min_y + 1;

    switch (graphPtr->anchor) {
	case TK_ANCHOR_N:
	    x -= width/2;
	    break;
	case TK_ANCHOR_NE:
	    x -= width;
	    break;
	case TK_ANCHOR_E:
	    x -= width;
	    y -= height/2;
	    break;
	case TK_ANCHOR_SE:
	    x -= width;
	    y -= height;
	    break;
	case TK_ANCHOR_S:
	    x -= width/2;
	    y -= height;
	    break;
	case TK_ANCHOR_SW:
	    y -= height;
	    break;
	case TK_ANCHOR_W:
	    y -= height/2;
	    break;
	case TK_ANCHOR_NW:
	    break;
	case TK_ANCHOR_CENTER:
	    x -= width/2;
	    y -= height/2;
	    break;
    }

    /*
     * Store the information in the item header.
     */
    /*
    graphPtr->header.x1 = x;
    graphPtr->header.y1 = y;
    graphPtr->header.x2 = x + width;
    graphPtr->header.y2 = y + height;
    */

    graphPtr->header.x1 = graphPtr->min_x + x;
    graphPtr->header.y1 = graphPtr->min_y + y;
    graphPtr->header.x2 = graphPtr->max_x + x;
    graphPtr->header.y2 = graphPtr->max_y + y;

#ifdef DEBUG
    printf("anchor x %f y %f\n", graphPtr->an_x, graphPtr->an_y);
    printf("ComputeGraphBbox x %d y %d width %d height %d min_x %f max_x %f\n",
	   x, y, width, height, graphPtr->min_x, graphPtr->max_x);

    printf("ComputeGraphBbox %d %d %d %d\n", graphPtr->header.x1,
	   graphPtr->header.y1, graphPtr->header.x2, graphPtr->header.y2);
#endif


}

void GraphScrollX(Tk_Canvas canvas,
		  Display *display,
		  GraphItem *graphPtr,
		  int new_x,
		  int old_x)
{
    int dx;
    Tk_Window tkwin = Tk_CanvasTkwin(canvas);
    int min_x, max_x, wx0, wx1;
    double y0, y1;

    min_x = graphPtr->x0;
    max_x = graphPtr->x1;

    if (graphPtr->vertical) {
	graphPtr->x0 = pixel_to_canvas(graphPtr->interp, tkwin, 'x', 0);
	graphPtr->x1 = pixel_to_canvas(graphPtr->interp, tkwin, 'x', Tk_Width(tkwin));
	if (graphPtr->inverty) {
	    y0 = graphPtr->max_y - graphPtr->y1 + graphPtr->min_y;
	    y1 = graphPtr->max_y - graphPtr->y0 + graphPtr->min_y;
	} else {
	    y0 = graphPtr->y0;
	    y1 = graphPtr->y1;
	}

	GraphPlotFunc(canvas, display, graphPtr, y0, y1);

    } else {
	dx = new_x - old_x; /* amount of scrolling */
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
	printf("oldx %d newx %d dx %d wx0 %d wx1 %d max_x %d min_x %d\n", old_x, new_x, dx, wx0, wx1, max_x, min_x);
#endif

	graphPtr->x0 = pixel_to_canvas(graphPtr->interp, tkwin, 'x', 0);
	graphPtr->x1 = pixel_to_canvas(graphPtr->interp, tkwin, 'x', Tk_Width(tkwin));

	GraphPlotFunc(canvas, display, graphPtr, wx0, wx1+1);
    }
}

void GraphScrollY(Tk_Canvas canvas,
		  Display *display,
		  GraphItem *graphPtr,
		  int new_y,
		  int old_y)
{
    Tk_Window tkwin = Tk_CanvasTkwin(canvas);
    int dy;
    int min_y, max_y, wy0, wy1, tmp;

    if (graphPtr->vertical) {

	min_y = graphPtr->y0;
	max_y = graphPtr->y1;

	dy = new_y - old_y; /* amount of scrolling */
	if (dy > 0) {
	    /* fill top edge */
	    wy1 = max_y + dy;
	    if (dy > max_y - min_y)
		wy0 = min_y + dy;
	    else
		wy0 = max_y;
	} else {
	    /* fill bottom edge */
	    wy0 = min_y + dy;
	    if (-dy > max_y - min_y)
		wy1 = max_y + dy;
	    else
		wy1 = min_y;
	}

#ifdef SCROLL
	printf("oldy %d newy %d dy %d wy0 %d wy1 %d max_y %d min_y %d\n", old_y, new_y, dy, wy0, wy1, max_y, min_y);

#endif

	if (graphPtr->inverty) {
	    tmp = wy0;
	    wy0 = graphPtr->max_y - wy1 + graphPtr->min_y;
	    wy1 = graphPtr->max_y - tmp + graphPtr->min_y;
	}

	graphPtr->y0 = pixel_to_canvas(graphPtr->interp, tkwin, 'y', 0);
	graphPtr->y1 = pixel_to_canvas(graphPtr->interp, tkwin, 'y', Tk_Height(tkwin));
	GraphPlotFunc(canvas, display, graphPtr, wy0, wy1+1);

    } else {

	graphPtr->y0 = pixel_to_canvas(graphPtr->interp, tkwin, 'y', 0);
	graphPtr->y1 = pixel_to_canvas(graphPtr->interp, tkwin, 'y', Tk_Height(tkwin));

	GraphPlotFunc(canvas, display, graphPtr, graphPtr->x0, graphPtr->x1);
	/*
	XDrawLine(display, graphPtr->pm, graphPtr->drawGC, 0, 500, 500, 0);
	*/
    }
}

static int b_search_pt(int key,
		       g_pt *array,
		       int num_pts)
{
    int found = 0;
    int middle, low, high;

    low = 0;
    high = num_pts-1;

    while (!found) {

	middle = (low + high) / 2;

	if (key < array[middle].x) {
	    high = middle - 1;
	} else if (key > array[middle].x) {
	    low = middle + 1;
	} else {
	    found = 1;
	    return middle;
	}
	if (low > high) {
	    break;
	}
    }
    return middle;
}

int b_search_nearest_line(int key,
			  gd_line *array,
			  int num_pts)
{
    int found = 0;
    int middle, low, high, index = 0;
    int min = INT_MAX;
    int dist;

    low = 0;
    high = num_pts-1;

    while (!found) {

	middle = (low + high) / 2;
	dist = abs(key - array[middle].x0);
#ifdef DEBUG
	printf("middle %d key %d array %f dist %d\n",
	       middle, key, array[middle].x0, dist);
#endif
	if (key < array[middle].x0) {
	    high = middle - 1;
	} else if (key > array[middle].x0) {
	    low = middle + 1;

	} else {
	    found = 1;
	    return middle;
	}
	if (dist < min) {
	    index = middle;
	    min = dist;
	}
	if (low > high) {
	    break;
	}

    }
#ifdef DEBUG
    printf("index %d middle %d\n", index, middle);
#endif
    return index;
}

int b_search_nearest_pt(int key,
			g_pt *array,
			int num_pts)
{
    int found = 0;
    int middle, low, high;

    low = 0;
    high = num_pts-1;

    while (!found) {

	middle = (low + high) / 2;

	if (key < array[middle].x) {
	    high = middle - 1;
	} else if (key > array[middle].x) {
	    low = middle + 1;

	} else {
	    found = 1;
	    return middle;
	}
	if (low > high) {
	    break;
	}

    }
    return middle;
}

static int b_search_line(int key,
			 gd_line *array,
			 int num_pts)
{
    int found = 0;
    int middle, low, high;

    low = 0;
    high = num_pts-1;

    while (!found) {

	middle = (low + high) / 2;

	if (key < array[middle].x0) {
	    high = middle - 1;
	} else if (key > array[middle].x0) {
	    low = middle + 1;

	} else {
	    found = 1;
	    return middle;
	}
	if (low > high) {
	    break;
	}

    }
    return middle;
}

void find_start_end_line(gd_line *array,
			 int n_lines,
			 int x0,
			 int x1,
			 int *First,
			 int *Last)
{
    int first = -1;
    int last = -1;

    /*
     * it isn't necessarily the case that I will have as many points as
     * sequence bases so need to check if each base is within my range
     */

    /* need to find array indices for first and last points */
    if (array[0].x0 >= x0) {
	first = 0;
    } else {
	first = b_search_line(x0, array, n_lines);
	if (first > 0)
	    first--;
    }

    last = b_search_line(x1, array, n_lines);

#ifdef DEBUG
    printf("LAST x1 %d %d\n", x1, last);
#endif
    last+=2;
    if (last > n_lines)
	last = n_lines;
#ifdef DEBUG
    printf("first %d last %d\n", first, last);
#endif

    *First = first;
    *Last = last;
}

void find_start_end_pt(g_pt *array,
		       int n_lines,
		       int x0,
		       int x1,
		       int *First,
		       int *Last)
{
    int first = -1;
    int last = -1;

    /*
     * it isn't necessarily the case that I will have as many points as
     * sequence bases so need to check if each base is within my range
     */

    /* need to find array indices for first and last points */
    if (array[0].x >= x0) {
	first = 0;
    } else {
	first = b_search_pt(x0, array, n_lines);
	if (first > 0)
	    first--;
    }
#ifdef DEBUG
    printf("first %f %d %d\n", array[0].x, x0, first);
#endif
    last = b_search_pt(x1, array, n_lines);

#ifdef DEBUG
    printf("LAST x1 %d %d\n", x1, last);
#endif
    last+=2;
    if (last > n_lines)
	last = n_lines;
#ifdef DEBUG
    printf("first %d last %d\n", first, last);
#endif
    *First = first;
    *Last = last;

}

double pixelx_to_world(double px,
		       GraphItem *graphPtr)
{
    double wx;

    wx = (px - graphPtr->M[0][1] - graphPtr->M[0][2]) / graphPtr->M[0][0];

#ifdef DEBUG
    print_matrix(graphPtr->M);
    printf("pixelx to world %f %f\n", px, wx);
#endif
    return wx;

}

double pixely_to_world(double py,
		       GraphItem *graphPtr)
{
    double wy;

    wy = (py - graphPtr->M[1][0] - graphPtr->M[1][2]) / graphPtr->M[1][1];

#ifdef DEBUG
    print_matrix(graphPtr->M);
    printf("pixely to world %f %f\n", py, wy);
#endif
    return wy;

}

/*
 * plot contiguous array of lines ie XDrawLines
 */
void Graph_Plot_Func(Tk_Canvas canvas,
		     Display *display,
		     GraphItem *graphPtr,
		     parray p,
		     int x0,  /* pixels */
		     int x1)  /* pixels */
{
    int i;
    int first = -1;
    int last = -1;
    XPoint staticPoints[MAX_STATIC_POINTS];
    XPoint *pointPtr;
    XPoint *pPtr;
    double wx0, wx1;
    int free = 0;
    g_pt pt;
    int height;
    Tk_Window tkwin = Tk_CanvasTkwin(canvas);

    if (p.n_pts == 0)
	return;

#ifdef DEBUG
    printf("Graph_Plot_Func x0 %d x1 %d\n", x0, x1);
#endif

#ifdef DEBUG
    printf("y0 %d max_y %f min_y %f\n", graphPtr->y0, graphPtr->max_y, graphPtr->min_y);
#endif

    if (graphPtr->vertical) {
	wx0 = pixely_to_world(x0, graphPtr);
	wx1 = pixely_to_world(x1, graphPtr);
    } else {
	wx0 = pixelx_to_world(x0, graphPtr);
	wx1 = pixelx_to_world(x1, graphPtr);
    }
    if (p.n_pts <= MAX_STATIC_POINTS) {
	pointPtr = staticPoints;
    } else {
	pointPtr = (XPoint *) ckalloc((unsigned) (p.n_pts *
						  sizeof(XPoint)));
	free++;
    }

    height = Tk_Height(tkwin);

    if (p.n_pts == 1) {
#ifdef TODO
	/* RasterDrawPoints(raster, coords, 1); */
#endif
    } else {
	/*
	 * it isn't necessarily the case that I will have as many points as
	 * sequence bases so need to check if each base is within my range
	 */
	find_start_end_pt(p.p_array, p.n_pts, wx0, wx1, &first, &last);

#ifdef DEBUG
	printf("first %d %f last %d %f\n", first, p.p_array[first].x,
	       last, p.p_array[last-1].x);
	printf("sx0 %d sy0 %d header.x1 %d header.x2 %d\n",
	       graphPtr->x0, graphPtr->y0,
	       graphPtr->header.x1, graphPtr->header.x2);
#endif

#ifdef DEBUG
	printf("max_x %f min_x %f x0 %d max_y %f min_y %f y0 %d\n",
	       graphPtr->max_x, graphPtr->min_x, graphPtr->x0,
	       graphPtr->max_y, graphPtr->min_y, graphPtr->y0);
#endif
	for (i = first, pPtr = pointPtr; i < last;  i++, pPtr++) {

	    canvas_array_pt(graphPtr, p.p_array[i], &pt);

	    /* FIXME need to check if need x0 offset in invertx/inverty
	      for graphs that do not start at the origin */
	    if (graphPtr->invertx) {
		pPtr->x = graphPtr->max_x - (pt.x + graphPtr->x0) + graphPtr->min_x;
	    } else {
		pPtr->x = pt.x - graphPtr->x0;
	    }

	    if (graphPtr->inverty) {


		pPtr->y = graphPtr->max_y - (pt.y + graphPtr->y0) + graphPtr->min_y;

#if 0
		pPtr->y = height + graphPtr->min_y - (pt.y + graphPtr->y0);
#endif
#ifdef DEBUG
		printf("max_y %f y %f y0 %d min_y %f pPtry %f\n",
		       graphPtr->max_y, y, graphPtr->y0, graphPtr->min_y, pPtr->y);
#endif
	    } else {
		pPtr->y = pt.y - graphPtr->y0;
	    }

	    pPtr->y += graphPtr->an_y;

#ifdef PRINT
	    printf("%d %f %f %d %d %f\n", i, p.p_array[i].x, p.p_array[i].y,
		   pPtr->x, pPtr->y, graphPtr->an_y);
#endif
	}
#ifdef DEBUG
	for (i = 0; i < last-first; i++) {
	    printf("%d %d\n", pointPtr[i].x, pointPtr[i].y);
	}
#endif
	XDrawLines(display, graphPtr->pm, graphPtr->drawGC, pointPtr,
		   last-first, CoordModeOrigin);
    }
    if (free)
	ckfree((char *)pointPtr);
}

/*
 * plot discrete array of lines ie XDrawSegments
 */
void Graph_Segment_Func(Tk_Canvas canvas,
			Display *display,
			GraphItem *graphPtr,
			darray d,
			int px0,  /* pixels */
			int px1)  /* pixels */
{
    XSegment *seg, *segptr;
    int i, num;
    int first, last, len;
    double wx0, wx1;
    m_coords *new_pt, *old_pt;
    g_pt start, end;
    int width, height;
    Tk_Window tkwin = Tk_CanvasTkwin(canvas);

#ifdef DEBUG
    printf("&&&&&&&&&&&&&&&&&Graph_Segment_Func %p %d %f x0 %d x1 %d\n", graphPtr, d.n_dlines, graphPtr->an_y, px0, px1);
#endif
    if (d.n_dlines == 0)
	return;

    old_pt = (m_coords *)ckalloc(sizeof(m_coords));
    new_pt = (m_coords *)ckalloc(sizeof(m_coords));

    segptr = seg = (XSegment*) ckalloc (sizeof (XSegment)*(d.n_dlines));

    if (graphPtr->vertical) {
	wx0 = pixely_to_world(px0, graphPtr);
	wx1 = pixely_to_world(px1, graphPtr);
    } else {
	wx0 = pixelx_to_world(px0, graphPtr);
	wx1 = pixelx_to_world(px1, graphPtr);
    }

    find_start_end_line(d.d_array, d.n_dlines, wx0, wx1, &first, &last);

    width = Tk_Width(tkwin);
    height = Tk_Height(tkwin);

    for (i = first, segptr = segptr; i < last; i++, segptr++) {

	canvas_array_line(graphPtr, d.d_array[i], &start, &end);

	if (graphPtr->invertx) {
	    segptr->x1 = graphPtr->max_x - (start.x + graphPtr->x0) + graphPtr->min_x;
	    segptr->x2 = graphPtr->max_x - (end.x + graphPtr->x0) + graphPtr->min_x;
	} else {
	    segptr->x1 = start.x - graphPtr->x0;
	    segptr->x2 = end.x - graphPtr->x0;
	}

	if (graphPtr->inverty) {


	    segptr->y1 = graphPtr->max_y - (start.y + graphPtr->y0) + graphPtr->min_y;
	    segptr->y2 = graphPtr->max_y - (end.y + graphPtr->y0) + graphPtr->min_y;
	    /*
	    segptr->y1 = height + graphPtr->min_y - (start.y + graphPtr->y0);
	    segptr->y2 = height + graphPtr->min_y - (end.y + graphPtr->y0);
	    */
#ifdef DEBUG
	    printf("SEGMENTY max_y %f min_y %f start.y %f end.y %f y0 %d y1 %d y2 %d\n",
		   graphPtr->max_y,  graphPtr->min_y, start.y, end.y, graphPtr->y0, segptr->y1, segptr->y2);
#endif
	} else {
	    segptr->y1 = start.y - graphPtr->y0;
	    segptr->y2 = end.y - graphPtr->y0;
	}

	segptr->y1 += graphPtr->an_y;
	segptr->y2 += graphPtr->an_y;

#ifdef DEBUG
	printf("dlines %d %f %f %f %f ", i, d.d_array[i].x0, d.d_array[i].y0,
	       d.d_array[i].x1, d.d_array[i].y1);

       printf("segment %d %d %d %d\n", segptr->x1, segptr->y1, segptr->x2, segptr->y2);
#endif
    }

    len = last-first;

    /*
     * fix to deal with x servers which can't cope with more than 2^16 lines
     */
    if (len < 32000) {
	XDrawSegments (display,
		       graphPtr->pm,
		       graphPtr->drawGC,
		       seg, len);
    } else {
	for (i = 0; i < len; i+=32000) {
	    if (len < i + 32000) {
		num = len - i;
	    } else {
		num = 32000;
	    }

	    XDrawSegments (display,
			   graphPtr->pm,
			   graphPtr->drawGC,
			   &seg[i], num);
	}
    }
    ckfree((char *)seg);
    ckfree((char *)old_pt);
    ckfree((char *)new_pt);
}

/*
 * plot discrete dots
 */
void Graph_Dot_Func(Tk_Canvas canvas,
		    Display *display,
		    GraphItem *graphPtr,
		    parray p,
		    int x0,  /* pixels */
		    int x1)  /* pixels */
{
    int i;
    XPoint staticPoints[MAX_STATIC_POINTS];
    XPoint *pointPtr;
    XPoint *pPtr;
    int pointwid = graphPtr->line_width ;
    int free_ptr = 0;
    int first, last;
    g_pt pt;
    double wx0, wx1;
    int height;
    Tk_Window tkwin = Tk_CanvasTkwin(canvas);

#ifdef DEBUG
    printf("Graph_Dot_Func\n");
#endif
    if (p.n_pts == 0)
	return;

    if (p.n_pts <= MAX_STATIC_POINTS) {
	pPtr = pointPtr = staticPoints;
    } else {
	pPtr = pointPtr = (XPoint *) ckalloc((unsigned) (p.n_pts * sizeof(XPoint)));
	free_ptr = 1;
    }

    height = Tk_Height(tkwin);

    wx0 = pixelx_to_world(x0, graphPtr);
    wx1 = pixelx_to_world(x1, graphPtr);

    find_start_end_pt(p.p_array, p.n_pts, wx0, wx1, &first, &last);

#ifdef DEBUG
    printf("y0 %d max_y %f min_y %f\n", graphPtr->y0, graphPtr->max_y, graphPtr->min_y);
#endif

    for (i = first, pPtr = pointPtr; i < last; i++, pPtr++) {

	canvas_array_pt(graphPtr, p.p_array[i], &pt);

	if (graphPtr->invertx) {
	    pPtr->x = graphPtr->max_x - (pt.x + graphPtr->x0) + graphPtr->min_x;
	} else {
	    pPtr->x = pt.x - graphPtr->x0;
	}

	if (graphPtr->inverty) {

	    pPtr->y = graphPtr->max_y - (pt.y + graphPtr->y0) + graphPtr->min_y;
#if 0
	    pPtr->y = height + graphPtr->min_y - (pt.y + graphPtr->y0);
#endif
	} else {
	    pPtr->y = pt.y - graphPtr->y0;
	}

#ifdef DEBUG
	printf("dot %f %f %d %d\n", p.p_array[i].x, p.p_array[i].y, pPtr->x, pPtr->y);
#endif
    }

    if (pointwid >= 2) {
	int halfwid = pointwid/2;
	for (i=0, pPtr = pointPtr; i < last-first; i++, pPtr++) {
#ifdef DEBUG
	    printf("xfillarc %d %d %d\n", halfwid, pPtr->x, pPtr->y);
#endif
	    XFillArc (display, graphPtr->pm, graphPtr->drawGC,
		      pPtr->x-halfwid, pPtr->y-halfwid,
		      pointwid, pointwid, 0, 360*64);
	    }
    } else {
#ifdef DEBUG
	printf("xdrawpoints \n");
#endif
	XDrawPoints (display, graphPtr->pm, graphPtr->drawGC, pointPtr,
		     last-first, CoordModeOrigin);
    }
    if (free_ptr)
	ckfree ((char *)pointPtr);

}

/*
 * plot array of filled rectangles
 */
void Graph_FilledRectangle_Func(Tk_Canvas canvas,
				Display *display,
				GraphItem *graphPtr,
				darray d,
				int px0,  /* pixels */
				int px1)  /* pixels */
{
    XRectangle *rect, *rectptr;
    int i;
    int first, last;
    double wx0, wx1;
    g_pt start, end;
    int width, height;
    Tk_Window tkwin = Tk_CanvasTkwin(canvas);

#ifdef DEBUG
    printf("&&&&&&&&&&&&&&&&&Graph_FilledRectanle_Func %p %d %f x0 %d x1 %d\n", graphPtr, d.n_dlines, graphPtr->an_y, px0, px1);
#endif
    if (d.n_dlines == 0)
	return;

    if (graphPtr->vertical) {
	wx0 = pixely_to_world(px0, graphPtr);
	wx1 = pixely_to_world(px1, graphPtr);
    } else {
	wx0 = pixelx_to_world(px0, graphPtr);
	wx1 = pixelx_to_world(px1, graphPtr);
    }

    find_start_end_line(d.d_array, d.n_dlines, wx0, wx1, &first, &last);

    width = Tk_Width(tkwin);
    height = Tk_Height(tkwin);

    rectptr = rect = (XRectangle*) ckalloc (sizeof (XRectangle)*(d.n_dlines));

    for (i = first ; i < last; i++, rectptr++) {

	canvas_array_line(graphPtr, d.d_array[i], &start, &end);
	if (start.x < end.x) {
	    rectptr->x = start.x - graphPtr->x0;
	    rectptr->width = abs(end.x - graphPtr->x0 - rectptr->x);
	} else {
	    rectptr->x = end.x - graphPtr->x0;
	    rectptr->width = abs(start.x - graphPtr->x0 - rectptr->x);
	}
	if (start.y < end.y) {
	    rectptr->y = start.y - graphPtr->y0;
	    rectptr->height = abs(end.y - graphPtr->y0 - rectptr->y);
	} else {
	    rectptr->y = end.y - graphPtr->y0;
	    rectptr->height = abs(start.y - graphPtr->y0 - rectptr->y);
	}

#ifdef DEBUG
	printf("dlines %d %f %f %f %f ", i, d.d_array[i].x0, d.d_array[i].y0,
	       d.d_array[i].x1, d.d_array[i].y1);

	printf("x %d y %d width %d height %d\n", rectptr->x, rectptr->y,
	       rectptr->width,  rectptr->height);

#endif
    }
    XFillRectangles(display, graphPtr->pm, graphPtr->drawGC, rect, d.n_dlines);
}

void GraphPlotFunc(Tk_Canvas canvas,
		   Display *display,
		   GraphItem *graphPtr,
		   int x1,
		   int x2)
{
    int i;
    Tk_Window tkwin = Tk_CanvasTkwin(canvas);

#ifdef DEBUG
    printf("GraphPlotFunc x1 %d x2 %d\n", x1, x2);

    printf("ndarrays %d nparrays %d\n",
	   graphPtr->graph->n_darrays, graphPtr->graph->n_parrays);
#endif


    /*
     * seems an awful HACK - need to reset these here for new plots created
     * when a previous plot has been scrolled
     */
    graphPtr->x0 = pixel_to_canvas(graphPtr->interp, tkwin, 'x', 0);
    graphPtr->x1 = pixel_to_canvas(graphPtr->interp, tkwin, 'x', Tk_Width(tkwin));

    graphPtr->y0 = pixel_to_canvas(graphPtr->interp, tkwin, 'y', 0);
    graphPtr->y1 = pixel_to_canvas(graphPtr->interp, tkwin, 'y', Tk_Height(tkwin));

    for (i = 0; i < graphPtr->graph->n_darrays; i++) {
	if (graphPtr->graph->d_arrays[i].type == G_LINES) {
	    /* discrete lines */
	    Graph_Segment_Func(canvas, display, graphPtr,
			       graphPtr->graph->d_arrays[i], x1, x2);
	} else if (graphPtr->graph->d_arrays[i].type == G_RECTANGLEFILL) {
	    Graph_FilledRectangle_Func(canvas, display, graphPtr,
				       graphPtr->graph->d_arrays[i], x1, x2);
	}
    }

   for (i = 0; i < graphPtr->graph->n_parrays; i++) {
       if (graphPtr->graph->p_arrays[i].type == G_LINE) {
	   /* contiguous line */
	   Graph_Plot_Func(canvas, display, graphPtr,
			   graphPtr->graph->p_arrays[i], x1, x2);
       } else if (graphPtr->graph->p_arrays[i].type == G_DOT) {
	   Graph_Dot_Func(canvas, display, graphPtr,
			  graphPtr->graph->p_arrays[i], x1, x2);
       }
   }
}

/*
 *--------------------------------------------------------------
 *
 * DisplayGraph --
 *
 *	This procedure is invoked to draw a graph item in a given
 *	drawable.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	ItemPtr is drawn in drawable using the transformation
 *	information in canvas.
 *
 *--------------------------------------------------------------
 */
static void DisplayGraph(Tk_Canvas canvas,
			 Tk_Item *itemPtr,
			 Display *display,
			 Drawable drawable,
			 int x,
			 int y,
			 int width,
			 int height)
{
    GraphItem *graphPtr = (GraphItem *) itemPtr;
    Tk_Window tkwin = Tk_CanvasTkwin(canvas);
    int win_width=0, win_height=0;
    short drawableX=0, drawableY=0;
    int graphX=0, graphY=0, graphWidth=0, graphHeight=0;

    if (x > graphPtr->header.x1) {
	graphX = x - graphPtr->header.x1;
	graphWidth = graphPtr->header.x2 - x;
    } else {
	graphX = 0;
	if ((x+width) < graphPtr->header.x2) {
	    graphWidth = x + width - graphPtr->header.x1;
	} else {
	    graphWidth = graphPtr->header.x2 - graphPtr->header.x1;
	}
    }

    if (y > graphPtr->header.y1) {
	graphY = y - graphPtr->header.y1;
	graphHeight = graphPtr->header.y2 - y;
    } else {
	graphY = 0;
	if ((y+height) < graphPtr->header.y2) {
	    graphHeight = y + height - graphPtr->header.y1;
	} else {
	    graphHeight = graphPtr->header.y2 - graphPtr->header.y1;
	}
    }

#ifdef DEBUG
    printf("DISPLAY GRAPH x %d y %d w %d h %d\n", x, y, width, height);
    printf("DISPLAY GRAPH graphx %d y %d w %d h %d\n", graphX, graphY, graphWidth, graphHeight);
#endif

    graphX = x;
    graphY = y;

    graphWidth = width;
    graphHeight = height;
    Tk_CanvasDrawableCoords(canvas, (double) x, (double) y,
	    &drawableX, &drawableY);

#ifdef DEBUG
    printf("drawableX %d drawableY %d\n", drawableX, drawableY);
#endif

    graphPtr->new_x = graphX;
    graphPtr->new_y = graphY;
    win_width = Tk_Width(tkwin);
    win_height = Tk_Height(tkwin);

    if (graphPtr->prev_width == 0)
	graphPtr->prev_width = win_width;

    if (graphPtr->prev_height == 0)
	graphPtr->prev_height = win_height;

    /* check if window has been resized */
    if (graphPtr->prev_width != win_width) {
#ifdef DEBUG
	printf("RESIZE\n");
#endif
	ResizeGraph(canvas, display, tkwin, graphPtr, graphPtr->prev_width,
		    graphPtr->prev_height, width, height);
	graphPtr->prev_width = win_width;
    }
    if (graphPtr->prev_height != win_height) {
#ifdef DEBUG
	printf("RESIZE hgt %d height %d\n", graphPtr->prev_height, win_height);
#endif
	ResizeGraph(canvas, display, tkwin, graphPtr, graphPtr->prev_width,
		    graphPtr->prev_height,  width, height);
	graphPtr->prev_height = win_height;
    }

    /* redraw everything */
    if (graphPtr->redraw) {
#ifdef DEBUG
	printf("REDRAW %d %d\n", graphX, graphX+graphWidth);
#endif
	/* clear before doing drawing */
	XFillRectangle (display, graphPtr->pm, graphPtr->copyGC,
			0, 0, win_width, win_height);

	if (graphPtr->vertical) {
	    GraphPlotFunc(canvas, display, graphPtr, graphY, graphY + graphHeight);
	} else {
	    GraphPlotFunc(canvas, display, graphPtr, graphX, graphX + graphWidth);
	}
	graphPtr->redraw = 0;
    } else if (((TkCanvas *)canvas)->flags & UPDATE_SCROLLBARS) {
	/* scroll event */
	int dx, dy;

#ifdef DEBUG
	printf("SCROLL\n");
#endif
	dx = graphPtr->new_x - graphPtr->x;
	dy = graphPtr->new_y - graphPtr->y;

#ifdef DEBUG
	printf("1:old=%d,%d new=%d,%d\n",
	       graphPtr->x, graphPtr->y,
	       graphPtr->new_x,graphPtr->new_y);
	printf("dx %d win_width %d win_width-dx %d\n", dx, win_width, win_width-dx);
	printf("dy %d win_height %d win_height-dy %d\n", dy, win_height, win_height-dy);
#endif

	/* FIXME: - overflow in XCopyArea? */
	/* width & height (wid-dx, hgt) are of type unsigned int */
	/* only scroll in x if need to */
	if (graphPtr->new_x != graphPtr->x) {

	    if (graphPtr->vertical) {
		XFillRectangle (display, graphPtr->pm,
				graphPtr->copyGC,
				0, 0,
				win_width, win_height);

		if (graphPtr->new_x > graphPtr->x) {
		    /* Copying bottom chunk to top, clear bottom edge */
		    XCopyArea (Tk_Display(tkwin), graphPtr->pm, graphPtr->pm,
			       graphPtr->copyGC,
			       dx, 0,
			       win_width - dx, win_height,
			       0, 0);
		} else if (graphPtr->new_x <= graphPtr->x){
		    /* Copying top chunk to bottom, clear top edge */
		    XCopyArea (Tk_Display(tkwin), graphPtr->pm, graphPtr->pm,
			       graphPtr->copyGC,
			       0, 0,
			       win_width + dx, win_height,
			       -dx, 0);
		}


	    } else {
		if (graphPtr->new_x > graphPtr->x) {
		    if (win_width > dx) {
#ifdef DEBUG
			printf("XCopyArea %d %d %d %d\n", dx, 0, win_width-dx, win_height);
#endif
#ifdef DEBUG
			printf("Copying right chunk to left, clear right edge\n");
#endif
			/* Copying right chunk to left, clear right edge */
			XCopyArea (Tk_Display(tkwin), graphPtr->pm,
				   graphPtr->pm, graphPtr->copyGC,
				   dx, 0,
				   win_width - dx, win_height,
				   0, 0);
		    }
#ifdef DEBUG
		    printf("XFillRec %d %d %d %d\n",
			   win_width - dx < 0 ? 0 : win_width - dx, 0, win_width, win_height);
#endif
		    XFillRectangle (display, graphPtr->pm,
				    graphPtr->copyGC,
				    win_width - dx < 0 ? 0 : win_width - dx, 0,
				win_width, win_height);
		} else if (graphPtr->new_x <= graphPtr->x){
		    /* Copying left chunk to right, clear left edge */
#ifdef DEBUG
		    printf("Copying left chunk to right, clear left edge\n");
#endif
		    if (win_width + dx > 0) {
			XCopyArea (Tk_Display(tkwin), graphPtr->pm,
				   graphPtr->pm, graphPtr->copyGC,
				   0, 0,
				   win_width + dx, win_height,
				   -dx, 0);
		    }
		    XFillRectangle (display, graphPtr->pm,
				    graphPtr->copyGC,
				    0, 0,
				    -dx > win_width ? win_width : -dx, win_height);
		}
	    }
#ifdef DEBUG
	    printf("graphscrollx new %d x %d\n", graphPtr->new_x, graphPtr->x);
#endif
	    GraphScrollX(canvas, display, graphPtr, graphPtr->new_x, graphPtr->x);
	}
	/* only scroll in y if need to */
	if (graphPtr->new_y != graphPtr->y) {

	    if (graphPtr->vertical) {
		if (graphPtr->new_y > graphPtr->y) {
		    if (win_height > dy) {
			XCopyArea (Tk_Display(tkwin), graphPtr->pm,
				   graphPtr->pm, graphPtr->copyGC,
				   0, dy,
				   width, win_height - dy,
				   0, 0);
		    }
		    XFillRectangle (display, graphPtr->pm,
				    graphPtr->copyGC,
				    0, win_height-dy < 0 ? 0 : win_height-dy,
				    win_width, win_height);
		} else if (graphPtr->new_y <= graphPtr->y){
		    if (win_height + dy > 0) {
			XCopyArea (Tk_Display(tkwin), graphPtr->pm,
				   graphPtr->pm, graphPtr->copyGC,
				   0, 0,
				   win_width, win_height + dy,
				   0, -dy);
		    }
		    XFillRectangle (display, graphPtr->pm,
				    graphPtr->copyGC,
				    0, 0,
				    win_width, -dy > win_height ? win_height :-dy);
		}
	    } else {
		XFillRectangle (display, graphPtr->pm,
				graphPtr->copyGC,
				0, 0,
				win_width, win_height);
#if 1
		if (graphPtr->new_y > graphPtr->y) {
		    /* Copying bottom chunk to top, clear bottom edge */
		    XCopyArea (Tk_Display(tkwin), graphPtr->pm, graphPtr->pm,
			       graphPtr->copyGC,
			       0, dy,
			       win_width, win_height - dy,
			       0, 0);
		} else if (graphPtr->new_y <= graphPtr->y){
		    /* Copying top chunk to bottom, clear top edge */
		    XCopyArea (Tk_Display(tkwin), graphPtr->pm, graphPtr->pm,
			       graphPtr->copyGC,
			       0, 0,
			       win_width, win_height + dy,
			       0, -dy);
		}
#endif
	    }
#ifdef DEBUG
	    printf("graphscrolly\n");
#endif

	    GraphScrollY(canvas, display, graphPtr, graphPtr->new_y, graphPtr->y);

	}
    } else {
	/* don't do anything */
#ifdef DEBUG
	printf("DON'T DISPLAY x %d y %d graphX %d drawableX %d x0 %d graphWidth %d origin %d w %d h %d graphY %d graphHeight %d win_width %d win_height %d\n", x, y, graphX,
	       drawableX, graphPtr->x0, graphWidth, drawableX - (graphX - graphPtr->x0), width, height, graphY, graphHeight, win_width, win_height);
#endif

	XSetClipOrigin(Tk_Display(tkwin), graphPtr->paintGC,
		       drawableX - (graphX - graphPtr->x0),
		       drawableY - (graphY - graphPtr->y0));

	XSetClipMask(Tk_Display(tkwin), graphPtr->paintGC, graphPtr->pm);

	XFillRectangle (display, drawable, graphPtr->paintGC,
			drawableX, drawableY, graphWidth, graphHeight);

	XSetClipOrigin(Tk_Display(tkwin), graphPtr->paintGC, 0, 0);
	/* must reset clipmask since tk is still using paintgc */
	XSetClipMask(Tk_Display(tkwin), graphPtr->paintGC, None);
	return;
    }
    /* copies entire graph pixmap on to the canvas */

    /* NB must remember that the drawable pixmap given to me by tk is only the
       size of graphX, graphY, graphWidth, graphHeight and it's top left hand
       corner is always drawableX,drawableY. I want to paste from pm the
       correct data so I set the clip origin accordingly.

                 drawableX, drawableY
       ***************************
       *         |----|          *
       *         |    |          *
       *         |    |          *
       *         |----|          *
       ***************************
       <--graphX->

       so origin of pm (****) relative to drawable (----) is
       drawableX - (graphX - graph->x0)
    */
    XSetClipOrigin(Tk_Display(tkwin), graphPtr->paintGC,
		   drawableX - (graphX - graphPtr->x0),
		   drawableY - (graphY - graphPtr->y0));

#ifdef DEBUG
    XWriteBitmapFile(Tk_Display(tkwin), "pm.bitmap", graphPtr->pm,
		     win_width, win_height, -1, -1);
#endif

    XSetClipMask(Tk_Display(tkwin), graphPtr->paintGC, graphPtr->pm);

    XFillRectangle (display, drawable, graphPtr->paintGC,
		    drawableX, drawableY, graphWidth, graphHeight);

    XSetClipOrigin(Tk_Display(tkwin), graphPtr->paintGC, 0, 0);

    /* must reset clipmask since tk is still using paintgc */

    XSetClipMask(Tk_Display(tkwin), graphPtr->paintGC, None);

    graphPtr->x = graphPtr->new_x;
    graphPtr->y = graphPtr->new_y;
}

/*
 *--------------------------------------------------------------
 *
 * GraphToPoint --
 *
 *	Computes the distance from a given point to a given
 *	graph, in canvas units.
 *
 * Results:
 *	The return value is 0 if the point whose x and y coordinates
 *	are coordPtr[0] and coordPtr[1] is inside the graph.  If the
 *	point isn't inside the graph then the return value is the
 *	distance from the point to the graph.
 *
 * Side effects:
 *	None.
 *
 *--------------------------------------------------------------
 */

	/* ARGSUSED */
static double
GraphToPoint(Tk_Canvas canvas,
	     Tk_Item *itemPtr,
	     double *coordPtr)
{
    GraphItem *graphPtr = (GraphItem *) itemPtr;
    g_pt pixel, world;
    double dist1, dist2, dist3, tmp;
    double min, min1, min2, min3;
    int i;
    int index;

    dist1 = dist2 = dist3 = DBL_MAX;
    min = min1 = min2 = min3 = DBL_MAX;

    pixel.x = coordPtr[0];
    pixel.y = coordPtr[1];

    world.x = coordPtr[0];
    world.y = coordPtr[1];

    if (graphPtr->inverty) {
#ifdef DEBUG
	printf("coord %f %f ", world.x, world.y);
#endif
	world.y = graphPtr->max_y - (world.y + graphPtr->y0) + graphPtr->min_y;
	pixel.y = world.y;

#ifdef DEBUG
	printf("max_y %f max_x %f coord %f %f min_y %f min_x %f\n",
	       graphPtr->max_y, graphPtr->max_x,
	       world.x, world.y, graphPtr->min_y, graphPtr->min_x);
#endif

    }

    /* need to convert coordPtr into world coords */
    if (graphPtr->vertical) {
	/* swap x and y over */
	tmp = pixelx_to_world(world.x, graphPtr);
	world.x = pixely_to_world(world.y, graphPtr);
	world.y = tmp;
    } else {
	world.x = pixelx_to_world(world.x, graphPtr);
	world.y = pixely_to_world(world.y, graphPtr);
    }

#ifdef DEBUG
    printf("GraphToPoint %f %f %f %f %f x %d y %d\n", graphPtr->line_width,
	   coordPtr[0], coordPtr[1],
	   world.x, world.y, graphPtr->x0, graphPtr->y0);
#endif

    /* FIXME - may have several arrays */
    for (i = 0; i < graphPtr->graph->n_darrays; i++) {
	if (graphPtr->graph->d_arrays[i].type == G_LINES || 
	    graphPtr->graph->d_arrays[i].type == G_RECTANGLEFILL) {
	    if (graphPtr->graph->d_arrays[i].n_dlines > 0) {
		index = b_search_nearest_line((int)world.x,
					      graphPtr->graph->d_arrays[i].d_array,
					      graphPtr->graph->d_arrays[i].n_dlines);
		find_nearest_line(graphPtr, graphPtr->graph->d_arrays[i].d_array,
				  graphPtr->graph->d_arrays[i].n_dlines, index,
				  pixel, index, 1, &min1);
		if (min1 < dist1)
		    dist1 = min1;
	    }
	}
    }

    for (i = 0; i < graphPtr->graph->n_parrays; i++) {
	if (graphPtr->graph->p_arrays[i].n_pts > 0) {
	    index = b_search_nearest_pt((int)world.x,
					graphPtr->graph->p_arrays[i].p_array,
					graphPtr->graph->p_arrays[i].n_pts);
	    if (graphPtr->graph->p_arrays[0].type == G_LINE) {
		find_nearest_cline(graphPtr, graphPtr->graph->p_arrays[i].p_array,
				   graphPtr->graph->p_arrays[i].n_pts, index,
				   pixel, index, 1, &min2);
		if (min2 < dist2)
		    dist2 = min2;
	    } else if (graphPtr->graph->p_arrays[0].type == G_DOT) {
		find_nearest_dot(graphPtr, graphPtr->graph->p_arrays[i].p_array,
				 graphPtr->graph->p_arrays[i].n_pts, index, pixel, index, 1, &min3);
		if (min3 < dist3)
		    dist3 = min3;
	    }
	}
    }

    min = dist1;

    if (dist2 < min)
	min = dist2;

    if (dist3 < min)
	min = dist3;

#ifdef DEBUG
    printf("graphtopoint nearest %f\n", min);
#endif
    return min;
}

/*
 *--------------------------------------------------------------
 *
 * GraphToArea --
 *
 *	This procedure is called to determine whether an item
 *	lies entirely inside, entirely outside, or overlapping
 *	a given rectangle.
 *
 * Results:
 *	-1 is returned if the item is entirely outside the area
 *	given by rectPtr, 0 if it overlaps, and 1 if it is entirely
 *	inside the given area.
 *
 * Side effects:
 *	None.
 *
 *--------------------------------------------------------------
 */

	/* ARGSUSED */
static int
GraphToArea(Tk_Canvas canvas,
	    Tk_Item *itemPtr,
	    double *rectPtr)
{

    printf("GraphToArea\n");
    return 0;
}

/*
 *--------------------------------------------------------------
 *
 * ScaleGraph --
 *
 *	This procedure is invoked to rescale a graph item in a
 *	canvas.  It is one of the standard item procedures for
 *	graph items, and is invoked by the generic canvas code.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The item referred to by itemPtr is rescaled so that the
 *	following transformation is applied to all point coordinates:
 *		x' = originX + scaleX*(x-originX)
 *		y' = originY + scaleY*(y-originY)
 *
 *--------------------------------------------------------------
 */

static void
ScaleGraph(Tk_Canvas canvas,
	   Tk_Item *itemPtr,
	   double originX,
	   double originY,
	   double scaleX,
	   double scaleY)
{
    GraphItem *graphPtr = (GraphItem *) itemPtr;
    ZoomValues *zoom_values;
    double M1[3][3], T[3][3];

#ifdef DEBUG
    printf("ScaleGraph %f %f %f %f\n", originX, originY, scaleX, scaleY);
    printf("BEFORE MIN3 %f %f MAX %f %f\n", graphPtr->min_x, graphPtr->min_y,
	   graphPtr->max_x, graphPtr->max_y);
#endif
    /* set redraw flag */
    graphPtr->redraw = 1;

    /* pixels */
    graphPtr->min_x = originX + scaleX*(graphPtr->min_x-originX);
    graphPtr->max_x = originX + scaleX*(graphPtr->max_x-originX);

    graphPtr->min_y = originY + scaleY*(graphPtr->min_y-originY);
    graphPtr->max_y = originY + scaleY*(graphPtr->max_y-originY);

    zoom_values = (ZoomValues *)malloc(sizeof(ZoomValues));
    zoom_values->scale_x = scaleX;
    zoom_values->scale_y = scaleY;
    zoom_values->origin_x = originX;
    zoom_values->origin_y = originY;

    make_zoom_matrix(M1, zoom_values);

    /* accumulate changes in graphPtr->M */
    copy_matrix(T, graphPtr->M);
    multiply_matrices(graphPtr->M, M1, T);

#ifdef DEBUG
    printf("MIN3 %f %f MAX %f %f\n", graphPtr->min_x, graphPtr->min_y,
	   graphPtr->max_x, graphPtr->max_y);
#endif
    ComputeGraphBbox(canvas, graphPtr);
    free(zoom_values);
}

/*
 *--------------------------------------------------------------
 *
 * TranslateGraph --
 *
 *	This procedure is called to move an item by a given amount.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The position of the item is offset by (xDelta, yDelta), and
 *	the bounding box is updated in the generic part of the item
 *	structure.
 *
 *--------------------------------------------------------------
 */

static void
TranslateGraph(Tk_Canvas canvas,
	       Tk_Item *itemPtr,
	       double deltaX,
	       double deltaY)
{
    GraphItem *graphPtr = (GraphItem *) itemPtr;
    double M1[3][3], T[3][3];

#ifdef DEBUG
    printf("TranslateGraph %f %f\n", deltaX, deltaY);
#endif

    /* set redraw flag */
    graphPtr->redraw = 1;

    /* pixels */
    /* FIXME - removed these lines because the only time I use move is for the
       stop codon plot and the movement is within the range of min & max so
       I don't need to add on delta. There maybe cases when this is not true
       in which case I need to also store the min and max of the actual data
    */
#if 0
    graphPtr->min_x = graphPtr->min_x + deltaX;
    graphPtr->max_x = graphPtr->max_x + deltaX;

    graphPtr->min_y = graphPtr->min_y + deltaY;
    graphPtr->max_y = graphPtr->max_y + deltaY;
#endif

    make_translation_matrix(M1, deltaX, deltaY);

    /* accumulate changes in graphPtr->M */
    copy_matrix(T, graphPtr->M);
    multiply_matrices(graphPtr->M, M1, T);

#ifdef DEBUG
    printf("MIN4 %f %f MAX %f %f\n", graphPtr->min_x, graphPtr->min_y,
	   graphPtr->max_x, graphPtr->max_y);
#endif
    ComputeGraphBbox(canvas, graphPtr);

}

/*
 *--------------------------------------------------------------
 *
 * GraphToPostscript --
 *
 *	This procedure is called to generate Postscript for
 *	bitmap items.
 *
 * Results:
 *	The return value is a standard Tcl result.  If an error
 *	occurs in generating Postscript then an error message is
 *	left in the interp's result, replacing whatever used to be there.
 *	If no error occurs, then Postscript for the item is appended
 *	to the result.
 *
 * Side effects:
 *	None.
 *
 *--------------------------------------------------------------
 */

static int
GraphToPostscript(Tcl_Interp *interp,
		  Tk_Canvas canvas,
		  Tk_Item *itemPtr,
		  int prepass)
{

    return 0;
}

/* convert pixel into canvas value */
int pixel_to_canvas(Tcl_Interp *interp,
		    Tk_Window tkwin,
		    char direction,
		    int pixel)
{
    char buf[1024];

    if (direction == 'x') {
	sprintf(buf, "%s canvasx %d", Tk_PathName(tkwin), pixel);
    } else {
	sprintf(buf, "%s canvasy %d", Tk_PathName(tkwin), pixel);
    }
    Tcl_Eval(interp, buf);
    return (atoi(Tcl_GetStringResult(interp)));

}

void ResizeGraph(Tk_Canvas canvas,
		 Display *display,
		 Tk_Window tkwin,
		 GraphItem *graphPtr,
		 int old_width,
		 int old_height,
		 int new_width,
		 int new_height)
{

    Pixmap newpm;

    newpm = Tk_GetPixmap (Tk_Display(tkwin),
			  Tk_WindowId(tkwin),
			  new_width, new_height,
			  1);

    XFillRectangle (display, newpm, graphPtr->copyGC,
		    0, 0, new_width, new_height);

    XCopyArea (display, graphPtr->pm, newpm,
	       graphPtr->copyGC, 0, 0, old_width, old_height, 0, 0);
    Tk_FreePixmap (display, graphPtr->pm);
    graphPtr->pm = newpm;

    /*
    GraphPlotFunc(canvas, display, graphPtr, graphPtr->x0, graphPtr->x1);
    */
}

double find_magnitude(g_pt *pt1,
		      g_pt *pt2)
{
    g_pt vector;

    vector.x = pt2->x - pt1->x;
    vector.y = pt2->y - pt1->y;

    return (sqrt((vector.x * vector.x) + (vector.y * vector.y)));
}

/* convert array[i] line into canvas coords */
void canvas_array_line(GraphItem *graphPtr,
		       gd_line array,
		       g_pt *start,
		       g_pt *end)
{
    m_coords new_pt, old_pt;
    double x0, y0, x1, y1;

    if (graphPtr->vertical) {
	y0 = array.x0;
	x0 = array.y0;
	y1 = array.x1;
	x1 = array.y1;
    } else {
	x0 = array.x0;
	y0 = array.y0;
	x1 = array.x1;
	y1 = array.y1;
    }

    make_coord(&old_pt, x0, y0);
    multiply_coords_by_matrix(&new_pt, graphPtr->M, &old_pt);

    (*start).x = new_pt.x;
    (*start).y = new_pt.y;

    make_coord(&old_pt, x1, y1);
    multiply_coords_by_matrix(&new_pt, graphPtr->M, &old_pt);

    (*end).x = new_pt.x;
    (*end).y = new_pt.y;
}

void canvas_array_pt(GraphItem *graphPtr,
		     g_pt array,
		     g_pt *pt)
{
    m_coords new_pt, old_pt;
    double x, y;

    if (graphPtr->vertical) {
	y = array.x;
	x = array.y;
    } else {
	x = array.x;
	y = array.y;
    }

    make_coord(&old_pt, x, y);
    multiply_coords_by_matrix(&new_pt, graphPtr->M, &old_pt);

    (*pt).x = new_pt.x;
    (*pt).y = new_pt.y;

}

void nearest_dist(g_pt pt,
		  g_pt start,
		  g_pt end,
		  double *min)
{
    double dist, dist1, dist2, d;
    g_pt pt3;

    d = (((pt.x - start.x) * (end.x - start.x)) +
	((pt.y - start.y) * (end.y - start.y))) /
	(find_magnitude(&end, &start) * find_magnitude(&end, &start));

#ifdef DEBUG
    printf("u %f\n", u);
#endif

    if (d < 0.0 || d > 1.0) {
#ifdef DEBUG
	printf("line is not perpendicular to line segment\n");
#endif
	dist1 = sqrt(((pt.x - start.x) * (pt.x - start.x)) +
		     ((pt.y - start.y) * (pt.y - start.y)));
	dist2 = sqrt(((pt.x - end.x) * (pt.x - end.x)) +
		     ((pt.y - end.y) * (pt.y - end.y)));
#ifdef DEBUG
	printf("dist1 %f dist2 %f\n", dist1, dist2);
#endif
	if (dist1 < *min)
	    *min = dist1;

	if (dist2 < *min)
	    *min = dist2;
    } else {

	pt3.x = start.x + d * (end.x - start.x);
	pt3.y = start.y + d * (end.y - start.y);

#ifdef DEBUG
	printf("Intersection %f %f\n", pt3.x, pt3.y);
#endif
	dist = find_magnitude(&pt, &pt3);

#ifdef DEBUG
	printf("dist %f\n", dist);
#endif
	if (dist < *min) {
	    *min = dist;
	}
    }
#ifdef DEBUG
    printf("min %f\n", *min);
#endif

}

/* finds the shortest distance from pt to a line in array */
void find_nearest_line(GraphItem *graphPtr,
		       gd_line *array,
		       int num_lines,
		       int index,
		       g_pt pt,
		       int start_x,
		       int counter,
		       double *min)
{
    g_pt start, end;

    canvas_array_line(graphPtr, array[index], &start, &end);

#ifdef DEBUG
    printf("index %d pt.x %f pt.y %f start.x %f start.y %f end.x %f end.y %f\n",
	   index, pt.x, pt.y, start.x, start.y, end.x, end.y);
#endif

    nearest_dist(pt, start, end, min);

#ifdef DEBUG
    printf("start.x %f pt.x %f %f %f\n", start.x, pt.x, pt.x + *min, pt.x - *min);
#endif
    if (start.x > (pt.x + *min)) {
#ifdef DEBUG
	printf("reached right hand edge\n");
#endif
	return;
    }

    if (start.x < (pt.x - *min)) {
#ifdef DEBUG
	printf("reached left hand edge\n");
#endif
	return;
    }

    counter *= -1;
    index = counter + start_x;
    if (counter > 0)
	counter++;

    if (index < 0) {
#ifdef DEBUG
	printf("reached left hand edge\n");
#endif
	return;
    }
    if (index >= num_lines) {
	return;
    }

    find_nearest_line(graphPtr, array, num_lines, index, pt, start_x, counter, min);
    return;
}

/* finds the shortest distance from pt to a continuous line in array */
void find_nearest_cline(GraphItem *graphPtr,
			g_pt *array,
			int n_pts,
			int index,
			g_pt pt,
			int start_x,
			int counter,
			double *min)
{
    g_pt start, end;

    canvas_array_pt(graphPtr, array[index], &start);
    if (index < n_pts-1)
	canvas_array_pt(graphPtr, array[index+1], &end);
    else
	return;

#ifdef DEBUG
    printf("index %d pt.x %f pt.y %f start.x %f start.y %f end.x %f end.y %f\n",
	   index, pt.x, pt.y, start.x, start.y, end.x, end.y);
#endif

    nearest_dist(pt, start, end, min);

#ifdef DEBUG
    printf("start.x %f pt.x %f min %f %f %f\n", start.x, pt.x, *min, pt.x + *min, pt.x - *min);
#endif
    if (start.x > (pt.x + *min)) {
#ifdef DEBUG
	printf("reached right hand edge\n");
#endif
	return;
    }

    if (start.x < (pt.x - *min)) {
#ifdef DEBUG
	printf("reached left hand edge\n");
#endif
	return;
    }

    counter *= -1;
    index = counter + start_x;
    if (counter > 0)
	counter++;

    if (index < 0) {
#ifdef DEBUG
	printf("reached left hand edge\n");
#endif
	return;
    }

    if (index >= n_pts) {
	return;
    }

    find_nearest_cline(graphPtr, array, n_pts, index, pt, start_x, counter, min);
    return;
}

void find_nearest_dot(GraphItem *graphPtr,
		      g_pt *array,
		      int n_pts,
		      int index,
		      g_pt pt,
		      int start_x,
		      int counter,
		      double *min)
{
    g_pt start;
    double x, y, square;

    canvas_array_pt(graphPtr, array[index], &start);

    x = (pt.x - start.x);
    y = (pt.y - start.y);
    square = x * x + y * y;

    if (square < *min) {
	*min = square;
    }

    if (start.x > (pt.x + *min)) {
#ifdef DEBUG
	printf("reached right hand edge\n");
#endif
	return;
    }

    if (start.x < (pt.x - *min)) {
#ifdef DEBUG
	printf("reached left hand edge\n");
#endif
	return;
    }

    counter *= -1;
    index = counter + start_x;
    if (counter > 0)
	counter++;

    if (index < 0) {
#ifdef DEBUG
	printf("reached left hand edge\n");
#endif
	return;
    }

    if (index >= n_pts)
	return;

    find_nearest_dot(graphPtr, array, n_pts, index, pt, start_x, counter, min);

}
