#if !defined(TRACE_PRINT_H)
#define TRACE_PRINT_H

#include <stdio.h>
#include "misc.h"
#include "tkTrace.h"
#include "postscript.h"

int ps_configure_trace(DNATrace *t, int argc, char **argv);
int ps_configure_trace_line(DNATrace *t, int argc, char **argv);
int free_ps_trace(PS_TRACE *ps_trace);
int trace_print(DNATrace *t, char *file_name);
int ps_trace_draw_trace(DNATrace *t, FILE *ps_file);
PS_POINT *ps_trace_segment(TRACE *trace, int first, int NPoints,
			   double x_scale, double y_scale, int y_max);
int ps_sequence_segment(DNATrace *t, int first, int NPoints,
			 PS_TEXT **A, PS_TEXT **C, PS_TEXT **G, PS_TEXT **T, PS_TEXT **gap,
			 int *nA, int *nC, int *nG, int *nT, int *nGap);
int ps_numbers_segment(DNATrace *t, int first, int NPoints, PS_TEXT **numbers, int *n_num);

#endif
