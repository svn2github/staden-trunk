#ifndef _CONFIDENCE_GRAPH_H
#define _CONFIDENCE_GRAPH_H

#define CONFIDENCE_GRAPH_QUALITY  0
#define CONFIDENCE_GRAPH_SECOND   1
#define CONFIDENCE_GRAPH_DISCREP  2
#define CONFIDENCE_GRAPH_DISCREP2 3
#define CONFIDENCE_GRAPH_DNP      4

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
    int mode;
    ruler_s *ruler;
} obj_confidence_graph;

int confidence_graph_reg(GapIO *io, Tcl_Interp *interp, char *frame,char 
			 *win_quality, int cons_id, ruler_s *ruler,
			 int mode /* QUALITY or DISCREP */);

#endif
