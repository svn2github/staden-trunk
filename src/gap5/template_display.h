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
    int fwd_col3;
    int rev_col3;
    int xhair_col;
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
    int *tdepth; /* paired template depth (consistent only) */
    int *sdepth; /* sequence depth */
    int depth_width;
    int plot_depth;
    double xhair_pos;
    double yhair_pos;
    int min_qual, max_qual; /* Filter parameters */
} template_disp_t;

/* If bit set we filter our this data type */
#define FILTER_PAIRED        (1<<0)
#define FILTER_SINGLE        (1<<1)
#define FILTER_CONSISTENT    (1<<2)
#define FILTER_INCONSISTENT  (1<<3)
#define FILTER_SPANNING      (1<<4)
#define FILTER_NONSPANNING   (1<<5)

int TDisp_Init(Tcl_Interp *interp);

template_disp_t *template_new(GapIO *io, int cnum,
			      Tcl_Interp *interp,
			      Tk_Raster *raster);
void template_destroy(template_disp_t *t);
int template_replot(template_disp_t *t);

#endif /* _TEMPLATE_DISPLAY_H_ */
