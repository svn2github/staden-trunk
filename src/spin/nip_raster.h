#ifndef _NIP_RASTER_H_
#define _NIP_RASTER_H_

#include <tcl.h>
#include "tkRaster.h"
#include "seq_reg.h"
#include "seq_raster.h"

#define NEW_WINDOW 0
#define SEQ_WINDOW 1
#define MAX_NUM_SEQ 100 /* max num of seq to be displayed in single raster */

void nip_raster_callback(int seq_num, void *obj, seq_reg_data *jdata);
char *nip_get_raster_frame(Tcl_Interp *interp, int seq_id, int type, int frame);
int DeregisterRasterWindow(int seq_id, RasterResult *old_result, 
			   char *raster_win);
int NipReSetRasterWindowSize(Tcl_Interp *interp, char *raster_win);
int nip_raster_select_cursor(int raster_id, Tk_Raster *raster, 
			     char *raster_win, int pos, int max_dist);
int GetRasterWindowSize(Tcl_Interp *interp, char *raster_win,
			double *x0, double *y0, double *x1, double *y1);
int NipAddRasterToWindow(Tcl_Interp *interp, char *raster_new);
int nip_raster_reg(Tcl_Interp *interp, char *raster_win, seq_id_dir *seq_array,
		   int num_seq_id);
int NipDeregisterRasterWindow(int seq_id, RasterResult *old_result,
			      char *raster_win);
#endif
