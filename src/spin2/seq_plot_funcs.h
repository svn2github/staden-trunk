#ifndef _NIP_PLOT_FUNCS_H_
#define _NIP_PLOT_FUNCS_H_

#include "seq_reg.h"

void graph_plot_func(void *obj, seq_reg_plot *plot);
void emboss_graph_plot_func(void *obj, seq_reg_plot *plot);
void gene_search_plot_func(void *obj, seq_reg_plot *plot);
void stick_plot_func(void *obj, seq_reg_plot *plot);
void stick_pair_plot_func(void *obj, seq_reg_plot *plot);
void dot_plot_middot_func(void *obj, seq_reg_plot *plot);
void dot_plot_dot_func(void *obj, seq_reg_plot *plot);
void dot_plot_line_func(void *obj, seq_reg_plot *plot);
void dot_plot_scoreline_func(void *obj, seq_reg_plot *plot);

#endif
