#ifndef _READPAIR_COVERAGE_H
#define _READPAIR_COVERAGE_H
/*
typedef struct {
    Tcl_Interp *interp;
    int **histogram;
    int *max;
    int *min;
    int t_max;
    int t_min;
    char frame[100];
    char c_win[100];
    int id;
    int cons_id;
    int linewidth;
    char colour[30];
    ruler_s ruler;
} obj_readpair_coverage;
*/
int readpair_coverage_reg(GapIO *io, Tcl_Interp *interp, char *frame,
			  char *rcov_win, int cons_id, ruler_s *ruler);

#endif
