#ifndef _CONFIDENCE_GRAPH_H
#define _CONFIDENCE_GRAPH_H

typedef struct {
    Tcl_Interp *interp;
    float **qual;
    float *max;
    float *min;
    float t_max;
    float t_min;
    char frame[100];
    char c_win[100];
    int id;
    int cons_id;
    int linewidth;
    char colour[30];
    ruler_s *ruler;
} obj_confidence_graph;

int confidence_graph_reg(GapIO *io, Tcl_Interp *interp, char *frame,char 
			 *win_quality, int cons_id, ruler_s *ruler);

#endif
