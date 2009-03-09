#ifndef _TEMPLATE_DISPLAY_H_
#define _TEMPLATE_DISPLAY_H_

#include <tcl.h>
#include <tkRaster.h>
#include "tg_gio.h"

typedef struct {
    GapIO *io;
    contig_t *contig;
    int crec;
    Tk_Raster *raster;
    Tk_Window tkwin;
    Tcl_Interp *interp;
    Tk_OptionTable optionTable;
    int map_col[256];
    int single_col;
    int span_col;
    int inconsistent_col;
    int fwd_col;
    int rev_col;
    int logy;
    int cmode;
    int ymode;
    int accuracy;
    int spread;
    int reads_only;
    double yzoom;
    double xzoom;
    int sep_by_strand;
    int yoffset;
    int ymin, ymax; /* visible extents of data in Y */
    int filter; /* bitmask */
} template_disp_t;

#define FILTER_PAIRED      (1<<0)
#define FILTER_CONSISTENT  (1<<1)

int TDisp_Init(Tcl_Interp *interp);

template_disp_t *template_new(GapIO *io, int cnum,
			      Tcl_Interp *interp,
			      Tk_Raster *raster);
void template_destroy(template_disp_t *t);
int template_replot(template_disp_t *t);

#endif /* _TEMPLATE_DISPLAY_H_ */
