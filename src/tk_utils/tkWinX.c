/*
 * USED ONLY IN THE WINDOWS VERSION
 *
 * A rewrite of various X functions used by our own Tk widgets that are not
 * available in the Tk X emulation library as used by windows/mac.
 * The routines have been rewritten as other primitives that do exist. This
 * slows them down, but as a first attempt it's valid.
 */


#include <tk.h>
#include <X11/Xlib.h>
#include <tkWinInt.h>
#include "tkWinX.h"



/* 14/1/99 johnt - Optimised Drawing:
   added  XDrawSegments
   re-implemented XDrawPoints
*/

typedef BOOL (CALLBACK *PolyPolyFunc)(HDC,CONST POINT *,CONST DWORD *, DWORD);

int tkpWinRopModes[] = {
    R2_BLACK,			/* GXclear */
    R2_MASKPEN,			/* GXand */
    R2_MASKPENNOT,		/* GXandReverse */
    R2_COPYPEN,			/* GXcopy */
    R2_MASKNOTPEN,		/* GXandInverted */
    R2_NOT,			/* GXnoop */
    R2_XORPEN,			/* GXxor */
    R2_MERGEPEN,		/* GXor */
    R2_NOTMERGEPEN,		/* GXnor */
    R2_NOTXORPEN,		/* GXequiv */
    R2_NOT,			/* GXinvert */
    R2_MERGEPENNOT,		/* GXorReverse */
    R2_NOTCOPYPEN,		/* GXcopyInverted */
    R2_MERGENOTPEN,		/* GXorInverted */
    R2_NOTMASKPEN,		/* GXnand */
    R2_WHITE			/* GXset */
};

/*
 * The following two raster ops are used to copy the foreground and background
 * bits of a source pattern as defined by a stipple used as the pattern.
 */

#define COPYFG		0x00CA0749 /* dest = (pat & src) | (!pat & dst) */
#define COPYBG		0x00AC0744 /* dest = (!pat & src) | (pat & dst) */

/*
 * Macros used later in the file.
 */

#define MIN(a,b)	((a>b) ? b : a)
#define MAX(a,b)	((a<b) ? b : a)


#ifdef NOTDEF
static POINT *
ConvertPoints(points, npoints, mode, bbox)
    XPoint *points;
    int npoints;
    int mode;			/* CoordModeOrigin or CoordModePrevious. */
    RECT *bbox;			/* Bounding box of points. */
{
    static POINT *winPoints = NULL; /* Array of points that is reused. */
    static int nWinPoints = -1;	    /* Current size of point array. */
    int i;

    /*
     * To avoid paying the cost of a malloc on every drawing routine,
     * we reuse the last array if it is large enough.
     */

    if (npoints > nWinPoints) {
	if (winPoints != NULL) {
	    ckfree((char *) winPoints);
	}
	winPoints = (POINT *) ckalloc(sizeof(POINT) * npoints);
	if (winPoints == NULL) {
	    nWinPoints = -1;
	    return NULL;
	}
	nWinPoints = npoints;
    }

    bbox->left = bbox->right = points[0].x;
    bbox->top = bbox->bottom = points[0].y;

    if (mode == CoordModeOrigin) {
	for (i = 0; i < npoints; i++) {
	    winPoints[i].x = points[i].x;
	    winPoints[i].y = points[i].y;
	    bbox->left = MIN(bbox->left, winPoints[i].x);
	    bbox->right = MAX(bbox->right, winPoints[i].x);
	    bbox->top = MIN(bbox->top, winPoints[i].y);
	    bbox->bottom = MAX(bbox->bottom, winPoints[i].y);
	}
    } else {
	winPoints[0].x = points[0].x;
	winPoints[0].y = points[0].y;
	for (i = 1; i < npoints; i++) {
	    winPoints[i].x = winPoints[i-1].x + points[i].x;
	    winPoints[i].y = winPoints[i-1].y + points[i].y;
	    bbox->left = MIN(bbox->left, winPoints[i].x);
	    bbox->right = MAX(bbox->right, winPoints[i].x);
	    bbox->top = MIN(bbox->top, winPoints[i].y);
	    bbox->bottom = MAX(bbox->bottom, winPoints[i].y);
	}
    }
    return winPoints;
}
#endif

/* We reuse these within the XDrawPoints and XDrawSegments to save on malloc time and memory usage */
static POINT* segpoints= NULL; /* array of start/end points of each segment */
static DWORD* seglens = NULL;   /* array detailing number of points in each segment (always 2) */
static int nsegpoints = -1;     /* number of segments currently allocated for */
static int nseglens = -1;       /* number of segment lengths currently allocated for */

static void
RenderPolyObject(dc, gc, wPoints, nPoints, mode, pen, func,pRect)
    HDC dc;
    GC gc;
    POINT* wPoints;
    int nPoints;
    int mode;
    HPEN pen;
    PolyPolyFunc func;
    RECT *pRect;
{
    HPEN oldPen;
    HBRUSH oldBrush;

    /* grow the segment lengths array if necessary */
    if (nPoints > (nseglens*2) ) {
	int i;
	if (seglens != NULL) {
	    ckfree((char *) seglens);
	}
	seglens   = (DWORD *) ckalloc(sizeof(DWORD)*nPoints/2);
	if (seglens == NULL) {
	    nseglens = -1;
	    return;
	}
	nseglens = nPoints/2;
	/* fill the segments lengths (always 2) */
	for(i=0;i<nseglens;i++)
	    seglens[i]=2;
    }

    if ((gc->fill_style == FillStippled
	    || gc->fill_style == FillOpaqueStippled)
	    && gc->stipple != None) {

	TkWinDrawable *twdPtr = (TkWinDrawable *)gc->stipple;
	HDC dcMem;
	LONG width, height;
	HBITMAP oldBitmap;
	int i;
	HBRUSH oldMemBrush;

	if (twdPtr->type != TWD_BITMAP) {
	    panic("unexpected drawable type in stipple");
	}

	/*
	 * Grow the bounding box enough to account for wide lines.
	 */

	if (gc->line_width > 1) {
	    pRect->left -= gc->line_width;
	    pRect->top -= gc->line_width;
	    pRect->right += gc->line_width;
	    pRect->bottom += gc->line_width;
	}

	width = pRect->right - pRect->left;
	height = pRect->bottom - pRect->top;

	/*
	 * Select stipple pattern into destination dc.
	 */

	SetBrushOrgEx(dc, gc->ts_x_origin, gc->ts_y_origin, NULL);
	oldBrush = SelectObject(dc, CreatePatternBrush(twdPtr->bitmap.handle));

	/*
	 * Create temporary drawing surface containing a copy of the
	 * destination equal in size to the bounding box of the object.
	 */

	dcMem = CreateCompatibleDC(dc);
	oldBitmap = SelectObject(dcMem, CreateCompatibleBitmap(dc, width,
		height));
	oldPen = SelectObject(dcMem, pen);
	BitBlt(dcMem, 0, 0, width, height, dc, pRect->left, pRect->top, SRCCOPY);

	/*
	 * Translate the object for rendering in the temporary drawing
	 * surface.
	 */

	for (i = 0; i < nPoints; i++) {
	    wPoints[i].x -= pRect->left;
	    wPoints[i].y -= pRect->top;
	}

	/*
	 * Draw the object in the foreground color and copy it to the
	 * destination wherever the pattern is set.
	 */

	SetPolyFillMode(dcMem, (gc->fill_rule == EvenOddRule) ? ALTERNATE
		: WINDING);
	oldMemBrush = SelectObject(dcMem, CreateSolidBrush(gc->foreground));
	(*func)(dcMem, wPoints, seglens, nPoints/2);
	BitBlt(dc, pRect->left, pRect->top, width, height, dcMem, 0, 0, COPYFG);

	/*
	 * If we are rendering an opaque stipple, then draw the polygon in the
	 * background color and copy it to the destination wherever the pattern
	 * is clear.
	 */

	if (gc->fill_style == FillOpaqueStippled) {
	    DeleteObject(SelectObject(dcMem,
		    CreateSolidBrush(gc->background)));
	    (*func)(dcMem, wPoints, seglens, nPoints/2);
	    BitBlt(dc, pRect->left, pRect->top, width, height, dcMem, 0, 0,
		    COPYBG);
	}

	SelectObject(dcMem, oldPen);
	DeleteObject(SelectObject(dcMem, oldMemBrush));
	DeleteObject(SelectObject(dcMem, oldBitmap));
	DeleteDC(dcMem);
    } else {
	oldPen = SelectObject(dc, pen);
	oldBrush = SelectObject(dc, CreateSolidBrush(gc->foreground));
	SetROP2(dc, tkpWinRopModes[gc->function]);

	SetPolyFillMode(dc, (gc->fill_rule == EvenOddRule) ? ALTERNATE
		: WINDING);

	(*func)(dc, wPoints, seglens, nPoints/2);

	SelectObject(dc, oldPen);
    }
    DeleteObject(SelectObject(dc, oldBrush));
}


void
XDrawPoints(display, d, gc, points, npoints, mode)
    Display* display;
    Drawable d;
    GC gc;
    XPoint* points;
    int npoints;
    int mode;
{
    HPEN pen;
    TkWinDCState state;
    HDC dc;
    int i,j;
    RECT bbox;

    if (d == None) {
	return;
    }

    /* allocate polypoints array if existing memory is not big enough */
    /* each X point takes 2 Windows points, as an Xpoint is drawn as a windows line */
    if (npoints > (nsegpoints/2)) {
	if (segpoints != NULL)
	    ckfree((char *) segpoints);
	segpoints = (POINT*) ckalloc(sizeof(POINT)*npoints*2);
	if (segpoints == NULL ) {
	    nsegpoints = -1;
	    return;
	}
	nsegpoints = npoints*2;
    }

    /* Fill in the array */
    bbox.left = bbox.right  = points[0].x;
    bbox.top  = bbox.bottom = points[0].y;
    if (mode == CoordModeOrigin) {
	for (i=0,j=0; i < npoints; i++,j+=2) {
	    segpoints[j].x = points[i].x;
	    segpoints[j].y = points[i].y;
	    segpoints[j+1].x = points[i].x+1;
	    segpoints[j+1].y = points[i].y+1;
	    bbox.left = MIN(bbox.left, points[i].x);
	    bbox.right = MAX(bbox.right, points[i].x+1);
	    bbox.top = MIN(bbox.top, points[i].y);
	    bbox.bottom = MAX(bbox.bottom, points[i].y+1);
	}
    } else {
	segpoints[0].x = points[0].x;
	segpoints[0].y = points[0].y;
	segpoints[1].x = points[0].x+1;
	segpoints[1].y = points[0].y+1;
	for (i=1,j=2; i < npoints; i++,j+=2) {
	    segpoints[j].x = points[i-1].x + points[i].x;
	    segpoints[j].y = points[i-1].y + points[i].y;
	    segpoints[j+1].x = segpoints[j].x+1;
	    segpoints[j+1].y = segpoints[j].y+1;
	    bbox.left = MIN(bbox.left, points[i].x);
	    bbox.right = MAX(bbox.right, points[i].x+1);
	    bbox.top = MIN(bbox.top, points[i].y);
	    bbox.bottom = MAX(bbox.bottom, points[i].y+1);
	}
    }

    dc = TkWinGetDrawableDC(display, d, &state);

    pen = CreatePen(PS_SOLID, 1, gc->foreground);

    RenderPolyObject(dc, gc, segpoints, npoints*2, mode, pen, PolyPolyline,&bbox);
    DeleteObject(pen);

    TkWinReleaseDrawableDC(d, dc, &state);
}

void
XDrawSegments(display, d, gc, segments, nsegs)
    Display* display;
    Drawable d;
    GC gc;
    XSegment* segments;
    int nsegs;
{
    HPEN pen;
    TkWinDCState state;
    HDC dc;
    int i,j;
    RECT bbox;

    if (d == None) {
	return;
    }

    /* allocate polypoints array if existing memory is not big enough */
    /* each X segment takes 2 Windows points */
    if (nsegs > (nsegpoints/2)) {
	if (segpoints != NULL)
	    ckfree((char *) segpoints);
	segpoints = (POINT*) ckalloc(sizeof(POINT)*nsegs*2);
	if (segpoints == NULL ) {
	    nsegpoints = -1;
	    return;
	}
	nsegpoints = nsegs*2;
    }
    /* Fill in the array */
    bbox.left = bbox.right  = segments[0].x1;
    bbox.top  = bbox.bottom = segments[0].y1;
    for (i=0,j=0; i < nsegs; i++,j+=2) {
	segpoints[j].x = segments[i].x1;
	segpoints[j].y = segments[i].y1;
	segpoints[j+1].x = segments[i].x2;
	segpoints[j+1].y = segments[i].y2;
	bbox.left = MIN(bbox.left, segments[i].x1);
	bbox.right = MAX(bbox.right, segments[i].x1);
	bbox.top = MIN(bbox.top, segments[i].y1);
	bbox.bottom = MAX(bbox.bottom, segments[i].y1);
	bbox.left = MIN(bbox.left, segments[i].x2);
	bbox.right = MAX(bbox.right, segments[i].x2);
	bbox.top = MIN(bbox.top, segments[i].y2);
	bbox.bottom = MAX(bbox.bottom, segments[i].y2);
    }

    dc = TkWinGetDrawableDC(display, d, &state);

    if ( gc->line_width > 1) {
	LOGBRUSH lb;
	DWORD style;

	lb.lbStyle = BS_SOLID;
	lb.lbColor = gc->foreground;
	lb.lbHatch = 0;

	style = PS_GEOMETRIC|PS_COSMETIC;
	switch (gc->cap_style) {
	    case CapNotLast:
	    case CapButt:
		style |= PS_ENDCAP_FLAT;
		break;
	    case CapRound:
		style |= PS_ENDCAP_ROUND;
		break;
	    default:
		style |= PS_ENDCAP_SQUARE;
		break;
	}
	pen = ExtCreatePen(style, gc->line_width, &lb, 0, NULL);
    } else {
	pen = CreatePen(PS_SOLID, gc->line_width, gc->foreground);
    }
    RenderPolyObject(dc, gc, segpoints, nsegs*2, CoordModeOrigin, pen, PolyPolyline,&bbox);
    DeleteObject(pen);

    TkWinReleaseDrawableDC(d, dc, &state);
}

/* XDrawPoint - use XDrawLine */
void XDrawPoint(Display *dsp, Drawable d, GC gc, int x, int y) {
    XDrawLine(dsp, d, gc, x, y, x+1, y+1);
}


/* XDrawRectangles - use XDrawRectangle */
void XDrawRectangles(Display *dsp, Drawable d, GC gc,
		     XRectangle *rectangles, int nrectangles) {
    int i;

    for (i=0; i<nrectangles; i++, rectangles++)
	XDrawRectangle(dsp, d, gc,
		       rectangles->x, rectangles->y,
		       rectangles->width, rectangles->height);
}

/* XSetClipRectangles - unimplemented */
int XSetClipRectangles( Display* d, GC gc, int x, int y, XRectangle* r, int w, int z )
{
   return 0;
}


/*
 * MJ: A version of this function is already provided by TK, so we'll defer
 * to its implementation, instead of this one to avoid link errors.
 */
#ifdef USE_CUSTOM_XCHANGEGC
void XChangeGC(Display *dsp, GC gp, unsigned long mask, XGCValues *values) {
    if (mask & GCFunction)          gp->function = values->function;
    if (mask & GCPlaneMask)         gp->plane_mask = values->plane_mask;
    if (mask & GCForeground)        gp->foreground = values->foreground;
    if (mask & GCBackground)        gp->background = values->background;
    if (mask & GCLineWidth)         gp->line_width = values->line_width;
    if (mask & GCLineStyle)         gp->line_style = values->line_style;
    if (mask & GCCapStyle)          gp->cap_style = values->cap_style;
    if (mask & GCJoinStyle)         gp->join_style = values->join_style;
    if (mask & GCFillStyle)         gp->fill_style = values->fill_style;
    if (mask & GCArcMode)           gp->arc_mode = values->arc_mode;
    if (mask & GCTile)              gp->tile = values->tile;
    if (mask & GCStipple)           gp->stipple = values->stipple;
    if (mask & GCTileStipXOrigin)   gp->ts_x_origin = values->ts_x_origin;
    if (mask & GCTileStipYOrigin)   gp->ts_y_origin = values->ts_y_origin;
    if (mask & GCFont)              gp->font = values->font;
    if (mask & GCSubwindowMode)
	gp->subwindow_mode = values->subwindow_mode;
    if (mask & GCGraphicsExposures)
	gp->graphics_exposures = values->graphics_exposures;
    if (mask & GCClipXOrigin)       gp->clip_x_origin = values->clip_x_origin;
    if (mask & GCClipYOrigin)       gp->clip_y_origin = values->clip_y_origin;
    if (mask & GCDashOffset)        gp->dash_offset = values->dash_offset;
    if (mask & GCDashList)          gp->dashes = values->dashes;
}
#endif
