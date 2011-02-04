/*
 * USED ONLY IN THE MacOS X "Aqua" VERSION
 *
 * A rewrite of various X functions used by our own Tk widgets that are not
 * available in the Tk X emulation library as used by the native aqua mac.
 * The routines have been rewritten as other primitives that do exist. This
 * slows them down, but as a first attempt it's valid.
 */


#include <tk.h>
#include <X11/Xlib.h>

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

