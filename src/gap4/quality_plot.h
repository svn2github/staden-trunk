#ifndef QUALITY_PLOT_H_
#define QUALITY_PLOT_H_

#include <tcl.h>
#include "ruler_display.h"
#include "io-reg.h"

typedef struct c_qual_t {
    int contig;
    int length;
    char *qual;
    int start;
    int end;
} c_qual;

typedef struct obj_t_qual_t {
    Tcl_Interp *interp;
    float cons_cutoff;
    int qual_cutoff;
    char window[100];
    char frame[100];
    int template_id;
    c_qual *quality;
    int num_contigs;
} obj_t_qual;

typedef struct obj_qual_t {
    Tcl_Interp *interp;
    float cons_cutoff;
    int qual_cutoff;
    char window[100];
    char frame[100];
    int template_id;
    int id;
    c_qual quality;
    ruler_s *ruler;
    cursor_s xhair;
    win **win_list;
    int num_wins;
    WorldPtr *world;
    CanvasPtr *canvas;
    StackPtr *zoom;
    cursor_t *cursor;
    int cursor_visible;
} obj_qual;

/*
 * Registers and initialises a quality buffer for a particular contig.
 */
int template_quality_reg(GapIO *io, 
			 Tcl_Interp *interp,
			 int *contig_array,
			 int num_contigs,
			 float cons_cutoff, 
			 int qual_cutoff, 
			 char *frame,
			 char *win_quality,
			 int template_id);

int quality_reg(GapIO *io, 
		Tcl_Interp *interp,
		int contig,
		int lreg,
		int rreg,
		float cons_cutoff, 
		int qual_cutoff, 
		char *frame,
		char *win_quality,
		ruler_s *ruler,
		cursor_s cursor);

#endif
