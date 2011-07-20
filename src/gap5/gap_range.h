/*
 *
 * gap_range.h - a tcl command to hold the x range for tracks.  For example
 *           template_display and depth_display
 *
 */
 
#ifndef _GAP_RANGE_H
#define _GAP_RANGE_H
 
#include <tcl.h>
#include <tk.h>
#include "tg_gio.h"

typedef struct {
    double y;          /* y coord */
    int col[3];        /* colours */
    int x[4];          /* coordinates */
    int t_strand;      /* template strand */
    int mq;            /* mapping qual */
    tg_rec rec;	       /* used for yspread */
} tline;


typedef struct {
    int filter;
    int min_qual;
    int max_qual;
    int c_mode;
    int accuracy;
} gap_filter_t;

typedef struct {
    double t;
    double s;
} gap_depth_t;
    

typedef struct {
    rangec_t *r;             // the range structure
    tline    *tl;            // line structure
    int ntl;                 // line size
    int nr;                  // number of entries in range
    GapIO *io;               // gap db
    contig_t *contig;        // contig
    double wx0;              // left most part of range (world coordinates)  
    double wx1;              // right most part of range (world coordinates)
    int width;	    	     // window width
    int max_height;          // depth plot height
    int template_mode;       // the drawing mode (for template_display)
    gap_filter_t old_filter; // old filter settings for comparison
    gap_filter_t new_filter; // new filter settings
    gap_depth_t  *depth;     // store depth
    tg_rec crec;
} gap_range_t;

/* If bit set we filter our this data type */
#define FILTER_PAIRED        (1<<0)
#define FILTER_SINGLE        (1<<1)
#define FILTER_CONSISTENT    (1<<2)
#define FILTER_INCONSISTENT  (1<<3)
#define FILTER_SPANNING      (1<<4)
#define FILTER_NONSPANNING   (1<<5)

/* define range beyond window size used in track calc */
#define GR_WINDOW_RANGE 1000

/* Global gap range_option */
extern Tk_CustomOption range_option;

int  GRange_Init(Tcl_Interp *interp);
void gap_range_destroy(gap_range_t *gr);
int  gap_range_test(gap_range_t *gr);
int  gap_range_recalculate(gap_range_t *gr, int width, double new_wx0, double new_wx1, int new_mode, int force); 
void gap_range_reset(gap_range_t *gr);
void set_filter(gap_range_t *gr, int filter, int min, int max, int mode, int accuracy);
int  gap_range_x(gap_range_t *gr, double ax_conv, double bx_conv, 
    	    	int forward_col, int reverse_col, int single_col, int span_col, int inconsistent_col,
		int force, int reads_only);


#endif
