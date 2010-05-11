#ifndef _DEPTH_TRACK_H_
#define _DEPTH_TRACK_H_

#include <tcl.h>
#include <tkRaster.h>
#include "tg_gio.h"
#include "gap_range.h"
#include "template_display.h"

typedef struct {
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
    int background;
    int logy;
    int cmode;
    int ymode;
    int accuracy;
    int spread;
    int reads_only;
    double yzoom;
    double xzoom;
    double yz;
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
    int min_sz;  /* For stacking mode */
    image_t *image;
    gap_range_t *gr;
    // remember the old parameters for reuse
    double wx0; // the old x range
    double wx1;
    tline *tl;
    int ntl;
    int old_filter;
    int old_min_qual, old_max_qual;
    int old_cmode;
    int old_accuracy;
    
} depth_track_t;

int DTrack_Init(Tcl_Interp *interp);

depth_track_t *depth_track_new(Tcl_Interp *interp, Tk_Raster *raster);
void depth_track_destroy(depth_track_t *t);
int depth_track_replot(depth_track_t *t);

#endif /* _DEPTH_TRACK_H_ */
