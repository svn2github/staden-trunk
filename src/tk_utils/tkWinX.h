#ifndef __TKWINX_H__
#define __TKWINX_H__


void XDrawPoint( Display* display, Drawable d, GC gc, int x, int y );
int  XSetClipRectangles( Display*, GC, int, int, XRectangle*, int, int );
void XDrawSegments(Display* display, Drawable d, GC gc, XSegment* segments, int nsegs );
void XDrawPoints( Display* display, Drawable d, GC gc, XPoint* points, int npoints, int mode );
void XDrawRectangles( Display* display, Drawable d, GC gc, XRectangle *rectangles, int nrectangles );


#endif
