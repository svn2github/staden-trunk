#include <tcl.h>
#include <tk.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "tkRaster.h"
#include "tkRasterBuiltIn.h"

static int ParseCoords (Tcl_Interp *interp, int argc, char **argv,
			double **coordptr, int *ncoordsptr)
/*         -----------
 * 
 *  Parses the argv strings as doubles. Tries to detect when a string
 *  represents a number and stops when the string is not numeric. Stores
 *  in *coordptr a pointer to an array (dynamically allocatd) with
 *  the parsed values and in *ncoordsptr the number of values found.
 *  Returns a standard Tcl Result.
 */
{
   int i;
   double* coord = (double*) malloc (sizeof (double) * argc);
   for (i = 0; 
	i < argc && (isdigit (argv [i][0]) || 
		     (argv [i][0]=='-' && isdigit (argv [i][1])));
	++i) {
      if (Tcl_GetDouble (interp, argv [i], & coord [i]) != TCL_OK) {
	 free (coord);
	 return TCL_ERROR;
      }
   }
   * coordptr = coord;
   * ncoordsptr = i;
   return TCL_OK;
}

static int DrawLine (Tcl_Interp *interp, Tk_Raster *raster,
		     ClientData clientdata, int argc, char **argv)
/*         --------
 *
 *  Implements raster drawing commands of the form:
 *
 *     <raster> draw_line <x0> <y0> <x1> <y1> ?... <xN> <yN>? ?-style <env>?
 *
 *     Draws a polygonal line linking the successive points with
 *     coordinates (<xi>,<yi>) using the color and linestyle of the raster's
 *     current drawing environment or the one specified in the -style option.
 *
 *  or
 *     <raster> fill_polygon <x0> <y0> <x1> <y1> ?... <xN> <yN>? ?-style <env>?
 * 
 *     Fills the polygon instead of drawing its outline. 
 */
{
   double *coords = NULL;
   int n;
   ClientData prevdrawenv = NULL;
   ClientData drawenv;
   int envno;

   if (ParseCoords (interp, argc-2, argv+2, &coords, &n) != TCL_OK) 
      goto error;

   if (n < 4 || n%2 == 1) {
      Tcl_AppendResult (interp, "insufficient # of coordinates", (char*)NULL);
      goto error;
   }

   if (argc > n+2) {
      if (strncmp (argv [n+2], "-style", strlen (argv [n+2])) == 0) {
	 if (argc != n+4) {
	    Tcl_AppendResult (interp, "expecting style #", (char*) NULL);
	    goto error;
	 }
	 GetDrawEnv (raster, &prevdrawenv);
	 if (Tcl_GetInt (interp, argv [n+3], &envno) != TCL_OK ||
	     DrawEnvIndex (interp, raster, envno, &drawenv) != TCL_OK ||      
	     SetDrawEnv (interp, raster, drawenv) != TCL_OK) {
	    goto error;
	 }
      }
      else {
	 Tcl_AppendResult (interp, "unrecognized option: ", argv [n+2], 
			   (char*) NULL);
	 goto error;
      }
   }
   
   if (argv [1][0] == 'f') {
      RasterFillPolygon (raster, coords, n/2);
   }
   else {
      RasterDrawLines (raster, coords, n/2);
   }
   free (coords);
   if (prevdrawenv != NULL) SetDrawEnv (interp, raster, prevdrawenv);
   return TCL_OK;

error:
   if (coords != NULL) free (coords);
   if (prevdrawenv != NULL) SetDrawEnv (interp, raster, prevdrawenv);
   return TCL_ERROR;   
}

static int DrawPoint (Tcl_Interp *interp, Tk_Raster *raster,
		      ClientData clientdata, int argc, char **argv)
/*         --------
 *
 *  Implements a raster drawing command of the form:
 *
 *     <raster> draw_point <x0> <y0> ? ... <xN> <yN>? ?-style <env>?
 *
 *  Draws the given points using the color of the raster's
 *  current drawing environment or the one specified in the -style option.
 */
{
   double *coords = NULL;
   int n;
   ClientData prevdrawenv = NULL;
   ClientData drawenv;
   int envno;

   if (ParseCoords (interp, argc-2, argv+2, &coords, &n) != TCL_OK) 
      goto error;

   if (n < 2 || n%2 == 1) {
      Tcl_AppendResult (interp, "insufficient # of coordinates", (char*)NULL);
      goto error;
   }

   if (argc > n+2) {
      if (strncmp (argv [n+2], "-style", strlen (argv [n+2])) == 0) {
	 if (argc != n+4) {
	    Tcl_AppendResult (interp, "expecting style #", (char*) NULL);
	    goto error;
	 }
	 GetDrawEnv (raster, &prevdrawenv);
	 if (Tcl_GetInt (interp, argv [n+3], &envno) != TCL_OK ||
	     DrawEnvIndex (interp, raster, envno, &drawenv) != TCL_OK ||      
	     SetDrawEnv (interp, raster, drawenv) != TCL_OK) {
	    goto error;
	 }
      }
      else {
	 Tcl_AppendResult (interp, "unrecognized option: ", argv [n+2], 
			   (char*) NULL);
	 goto error;
      }
   }
   
   RasterDrawPoints (raster, coords, n/2);
   free (coords);
   if (prevdrawenv != NULL) SetDrawEnv (interp, raster, prevdrawenv);
   return TCL_OK;

error:
   if (coords != NULL) free (coords);
   if (prevdrawenv != NULL) SetDrawEnv (interp, raster, prevdrawenv);
   return TCL_ERROR;   
}

static int DrawRect (Tcl_Interp *interp, Tk_Raster *raster,
		     ClientData clientdata, int argc, char **argv)
/*         ---------
 *
 *  Implements a raster drawing command of the form:
 *
 *     <raster> xxx_rectangle <x0> <y0> <x1> <y1> ?... <xN> <yN>?
 *              ?-style <env>?
 *
 *  (where xxx is 'draw' or 'fill)
 *
 *  Draws (or fills) the given rectangles using the settings of the raster's
 *  current drawing environment or the one specified in the -style option.
 */
{
   double *coords = NULL;
   int n;
   ClientData prevdrawenv = NULL;
   ClientData drawenv;
   int envno;

   if (ParseCoords (interp, argc-2, argv+2, &coords, &n) != TCL_OK) 
      goto error;

   if (n < 4 || n%4 != 0) {
      Tcl_AppendResult (interp, "insufficient # of coordinates", (char*)NULL);
      goto error;
   }

   if (argc > n+2) {
      if (strncmp (argv [n+2], "-style", strlen (argv [n+2])) == 0) {
	 if (argc != n+4) {
	    Tcl_AppendResult (interp, "expecting style #", (char*) NULL);
	    goto error;
	 }
	 GetDrawEnv (raster, &prevdrawenv);
	 if (Tcl_GetInt (interp, argv [n+3], &envno) != TCL_OK ||
	     DrawEnvIndex (interp, raster, envno, &drawenv) != TCL_OK ||      
	     SetDrawEnv (interp, raster, drawenv) != TCL_OK) {
	    goto error;
	 }
      }
      else {
	 Tcl_AppendResult (interp, "unrecognized option: ", argv [n+2], 
			   (char*) NULL);
	 goto error;
      }
   }
   
   if (argv [1][0] == 'f') {
      RasterFillRectangles (raster, coords, n/4);
   }
   else {
      RasterDrawRectangles (raster, coords, n/4);
   }
   free (coords);
   if (prevdrawenv != NULL) SetDrawEnv (interp, raster, prevdrawenv);
   return TCL_OK;

error:
   if (coords != NULL) free (coords);
   if (prevdrawenv != NULL) SetDrawEnv (interp, raster, prevdrawenv);
   return TCL_ERROR;   
}

int RasterBuiltInInit (interp)
/*  -----------------
 *
 *  Initializes the built in raster drawing commands
 */
     Tcl_Interp * interp;
{
   return (RasterAddPrimitive (interp, "draw_line", DrawLine,
			   (RasterPrimInitProc*) NULL,
			   (RasterPrimFreeProc*) NULL) == TCL_OK &&
	   RasterAddPrimitive (interp, "fill_polygon", DrawLine,
			   (RasterPrimInitProc*) NULL,
			   (RasterPrimFreeProc*) NULL) == TCL_OK &&
	   RasterAddPrimitive (interp, "draw_point", DrawPoint,
			   (RasterPrimInitProc*) NULL,
			   (RasterPrimFreeProc*) NULL) == TCL_OK &&
	   RasterAddPrimitive (interp, "draw_rectangle", DrawRect,
			   (RasterPrimInitProc*) NULL,
			   (RasterPrimFreeProc*) NULL) == TCL_OK &&
	   RasterAddPrimitive (interp, "fill_rectangle", DrawRect,
			   (RasterPrimInitProc*) NULL,
			   (RasterPrimFreeProc*) NULL) == TCL_OK) 
      ? TCL_OK : TCL_ERROR;
}
