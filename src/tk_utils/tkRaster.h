/* tkRaster.h (Header for module tkRaster)
 * ----------
 *
 */

#ifndef TK_RASTER_H
#define TK_RASTER_H

#include <tk.h>

/*--------------------------------------------------------------------------- 
 *
 *  A Tk_Raster is here defined as void. This was done on purpose -
 *  it is a simple way of saying that Tk_Raster should only be accessed
 *  from inside the tkRaster.c module. All other modules should 
 *  only make use of the procedures defined here.
 */
typedef void Tk_Raster;

/* Flags */
/* Redraw only the borders */
#define RASTER_BORDER  (1<<0)

/* Redraw all of text and trace widgets, but not borders */
#define RASTER_REDRAW  (1<<1)

/* Only redraw the bits changed by scrolling */
#define RASTER_SCROLL  (1<<2)

/* Forces a redraw of everything */
#define RASTER_RESCALE (1<<3)

/* A display request has been registered. */
#define RASTER_WAITING (1<<4)

/* The cursor needs redrawing */
#define RASTER_CURSOR  (1<<5)

#define RASTER_INIT          0
#define RASTER_REPLOT_ALL    1
#define RASTER_REPLOT_SLIVER 2
#define RASTER_REPLOT_ZOOM   3

/* 
 * A Drawing Environment is a set of arguments that affect the
 * drawing of primitives (color, line thickness, fill pattern etc). It
 * is implemented as a  XGCValues structure plus a bitmask. A raster may
 * have many of these loaded at the same time
 */
typedef struct {
   XGCValues gcValues;
   unsigned long valMask;
   XColor * fgColor, * bgColor;
   int index;
   GC drawGC;
} DrawEnvironment;

/*---------------------------------------------------------------------------
 *
 *  Below are defined built-in drawing functions for rasters. They should
 *  be enough for extending the Raster to draw any object (see the
 *  function RasterAddPrimitive). All object coordinates are specified in
 *  "world" units (but see functions RasterToWorld and WorldToRaster). 
 *  Attributes such as color, linestyle etc using for drawing are those
 *  of the current drawing environment of the raster (see also functions 
 *  DrawEnvIndex and SetDrawEnv).
 */
int Raster_Init (Tcl_Interp *interp);
extern int RasterCmd(ClientData clientData, Tcl_Interp *interp, int argc, 
	      char **argv);

extern void RasterDrawPoints _ANSI_ARGS_((Tk_Raster*, double* coord,int npts));
/* 
 *  Draws the points (coord [2i], coord [2i+1]) for 0 <= i < npts
 */
extern void RasterDrawPoint _ANSI_ARGS_((Tk_Raster*, int x, int y));

extern void RasterDrawLines _ANSI_ARGS_((Tk_Raster*,double* coord,int npts));
/* 
 *  Draws the npts-1 line segments connecting the successive
 *  'npts' points (coord [2i], coord [2i+1]), for 0 <= i < npts. 
 */
extern void RasterDrawLine _ANSI_ARGS_((Tk_Raster*, int x1, double y1, int x2, 
					double y2));

extern void RasterDrawRectangles _ANSI_ARGS_ ((Tk_Raster*, double* coord, 
					       int nrects));
/*     
 *  Draws rectangles with sides aligned with the coordinate axes and 
 *  with diagonal corners at (coord [4i], coord [4i+1]) and 
 *  (coord [4i+2], coord [4i+3]) for  0 <= i < nrects, 
 */

extern void RasterFillRectangles _ANSI_ARGS_ ((Tk_Raster*, double* coord, 
					       int nrects));
/*     
 *  Draws filled rectangles with sides aligned with the coordinate axes and 
 *  with diagonal corners at (coord [4i], coord [4i+1]) and 
 *  (coord [4i+2], coord [4i+3]) for  0 <= i < nrects, 
 */


extern void RasterFillPolygon _ANSI_ARGS_ ((Tk_Raster*, double* coord,
					   int npts));
/* 
 *  Fills the polygon bounded by the line segments connecting the successive
 *  'npts' points (coord [2i], coord [2i+1]), for 0 <= i < npts. 
 *  If the last point is not equal to the first point, a line segment 
 *  connecting them is added.
 */

/*--------------------------------------------------------------------------
 *
 *  The user may specify many 'drawing environments' which are essentially 
 *  Graphics Contexts. These are assigned numbers (indices), being 0 the number
 *  of the default. DrawEnvIndex performs a search in the Raster structure
 *  for a given drawing environment index and, if successful,
 *  returns a pointer to the structure that holds that drawing environment's
 *  data (the ClientData argument). 
 *  SetDrawEnv then can be used to load the parameters of the
 *  drawing environment into the Raster's GC.
 *  GetDrawEnv is used to get the Raster's current drawing environment.
 */

extern int DrawEnvIndex _ANSI_ARGS_((Tcl_Interp*, Tk_Raster*,int,ClientData*));
extern int SetDrawEnv _ANSI_ARGS_ ((Tcl_Interp*, Tk_Raster*, ClientData));
extern void GetDrawEnv _ANSI_ARGS_ ((Tk_Raster*, ClientData*));

/*--------------------------------------------------------------------------
 *
 *  SetRasterModifiedArea may be used by primitive implementations
 *  to let the raster know what part of the pixmap was modified. 
 *  If not called, raster will assume that all pixmap was modified.
 *  If implementing a raster command that does not actually modify the
 *  pixmap, this routine should be called with 0 0 0 0 as args.
 */
extern void SetRasterModifiedArea _ANSI_ARGS_((Tk_Raster* raster,
					       int rx0, int ry0, 
					       int rx1, int ry1));

/*---------------------------------------------------------------------------
 *
 *  Utility functions for handling world-to-raster-and-back  transformations
 */

extern void SetRasterCoords _ANSI_ARGS_((Tk_Raster*, 
					 double x0, double y0,
					 double x1, double y1));
/*    
 *  (x0,y0) are the coordinates of the upper-left corner and 
 *  (x1,y1) are the coordinates of the lower-right corner
 */

extern void GetWorldToRasterConversion(Tk_Raster *raster, double *ax, double *ay, 
    	    	    	    	    	double *bx, double *by);

/*
 * get the conversion factors for use outside
 * for speed rather than anything else
 */

extern void WorldToRaster _ANSI_ARGS_ ((Tk_Raster*, double wx, double wy,
					int* rx, int* ry));
/*
 *  (wx,wy) are  world coordinates
 *  (*rx,*ry) are set to the corresponding raster coordinates
 */

extern void RasterToWorld _ANSI_ARGS_ ((Tk_Raster*, int rx, int ry, 
					double* wx, double* wy));
/*
 *  (rx,ry) are raster coordinates
 *  (*wx,*wy) are set to the corresponding world coordinates
 */

/*--------------------------------------------------------------------------
 * 
 *  A new primitive is to be implemented as three procedures:
 *
 *  . A RasterPrimInitProc, which sets up some storage to be maintained
 *     	 	by the raster module.
 *
 *  . A RasterPrimDrawProc, which actually draws the primitive.
 *
 *  . A RasterPrimFreeProc, which frees up the storage allocated by 
 *           	RasterPrimInitProc.
 * 
 *  The 'Init' and 'Free' procs are called just after the raster is created
 *  and just before the raster is destroyed respectively. In order
 *  to implement a new primitive, one must write these three routines (or
 *  just a RasterPrimDrawProc if the primitive does not require
 *  storage to be kept per window) and then put a call to
 *  RasterAddPrimitive in the application initialization
 *  routine (usually called Tcl_AppInit), just after the call to RasterInit.
 */

typedef int RasterPrimDrawProc _ANSI_ARGS_((Tcl_Interp * interp,
					    Tk_Raster * raster,
					    ClientData data,
					    int argc, 
					    char* argv [])); 

typedef int RasterPrimInitProc _ANSI_ARGS_((Tcl_Interp * interp,
					    Tk_Raster* raster,
					    ClientData * dataptr));

typedef void RasterPrimFreeProc _ANSI_ARGS_((Tk_Raster* raster,
					     ClientData data));
					    
extern int RasterAddPrimitive _ANSI_ARGS_((Tcl_Interp * interp,
					   char * primitivename, 
					   RasterPrimDrawProc * drawproc,
					   RasterPrimInitProc * initproc,
					   RasterPrimFreeProc * freeproc));

/*---------------------------------------------------------------------------
 *
 *  RasterInit: initializes the Raster Widget module. This involves basically
 *  setting up the built in raster primitives. This procedure
 *  should be called just once (During the application initialization)
 */
extern int RasterInit _ANSI_ARGS_((Tcl_Interp * interp));

/*---------------------------------------------------------------------------
 *  RasterCmd: implements the raster widget command, i.e., 
 *
 *     raster <pathName> ?option? ... ?option?
 *
 */
extern int RasterCmd _ANSI_ARGS_((ClientData clientData, Tcl_Interp *interp, 
				  int argc,  char **argv));


/*---------------------------------------------------------------------------
 *
 *  The following are utility functions to access elements of a Tk_Raster.
 *  Normally, there's no need to use them - drawing should be done
 *  by calling the RasterDrawXXX functions. 
 *
 *  (The functions below are needed only if you are implementing something
 *   that has to use some weird X drawing function)
 */


extern Drawable GetRasterDrawable _ANSI_ARGS_(( Tk_Raster* ));
extern Display* GetRasterDisplay _ANSI_ARGS_(( Tk_Raster* ));
extern Tk_Window GetRasterTkWin _ANSI_ARGS_(( Tk_Raster* ));
extern GC GetRasterGC _ANSI_ARGS_(( Tk_Raster*));

char *
GetRasterColour (Tcl_Interp *interp,
		 Tk_Raster *RasterPtr,
		 int index);

int GetRasterLineWidth (Tcl_Interp *interp,
			Tk_Raster *RasterPtr,
			int index);
/*
 * C interface to tcl routine envcreate
 */
int
CreateDrawEnviron (Tcl_Interp *interp,
		   Tk_Raster *RasterPtr,
		   int argc, char **argv);

/*
 * C interface to tcl routine envset
 */
int
SetDrawEnviron (Tcl_Interp *interp,
		Tk_Raster *RasterPtr,
		int index);

void
tk_RasterClear(Tk_Raster *RasterPtr);

int
RasterSetWorldScroll(Tk_Raster *rasterptr,
		     double x0, double y0, 
		     double x1, double y1);

void
RasterGetWorldScroll(Tk_Raster *rasterptr,
		     double *x0, double *y0, 
		     double *x1, double *y1);

void GetRasterCoords(Tk_Raster *rasterptr, 
		     double *wx0, double *wy0, double *wx1, double *wy1);

void RasterWinSize(Tk_Raster *rasterptr, int *width, int *height);


void tk_RasterRefresh(Tk_Raster *RasterPtr);
void RasterInitPlotFunc(Tk_Raster *rasterptr,
		   void (*plot_func)(Tk_Raster *raster, char *raster_win, 
				     int redraw_y, int x0, int y0, int x1, 
				     int y1));
void RasterDeletePlotFunc(Tk_Raster *rasterptr);
void RasterResetWorldScroll(Tk_Raster *rasterptr);
double RasterGetXMag(Tcl_Interp *interp, double old_range, double new_range);
double RasterGetYMag(Tcl_Interp *interp, double old_range, double new_range);
double RasterSetXMag(Tcl_Interp *interp, double old_range, int value);
double RasterSetYMag(Tcl_Interp *interp, double old_range, int value);
int GetBgPixel(Tk_Raster* RasterPtr);
int GetFgPixel(Tcl_Interp* interp, Tk_Raster* RasterPtr, int index);
int SetFgPixel(Tcl_Interp* interp, Tk_Raster* RasterPtr, int index, 
	       int fg_pixel);
void RasterCallPlotFunc(Tk_Raster *rasterptr, int job, int x0, int y0, 
			int x1, int y1);

/* removed static kfs 10.05.00 */
int ConfigInfoDrawEnv (Tcl_Interp* interp, Tk_Raster* RasterPtr, 
				   DrawEnvironment* drawEnv,
				   int argc, char **argv);
void RasterDrawSegments (Tk_Raster * raster, double* coord, int nsegs);

#endif /* TK_RASTER_H */



