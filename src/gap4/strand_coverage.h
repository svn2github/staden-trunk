#ifndef _STRAND_COVERAGE_H
#define _STRAND_COVERAGE_H

#include "canvas_box.h"

typedef struct {
    Tcl_Interp *interp;
    int **forward;
    int **reverse;
    int forward_offset;
    int reverse_offset;
    int strand; /* 1=for 2=rev 3=both */
    int problems; /* 1=coverage 2=problems */
    char frame[100];
    char c_win[100];
    int id;
    int cons_id;
    int linewidth;
    char colour1[30];
    char colour2[30];
} obj_strand_coverage;

int strand_coverage_reg(GapIO *io, Tcl_Interp *interp, char *frame,char	*win, 
			int cons_id, int strand, int problems);

#endif
