#ifndef _READING_COVERAGE_H
#define _READING_COVERAGE_H

typedef struct {
    Tcl_Interp *interp;
    int **histogram1; /* forward */
    int **histogram2; /* reverse */
    int *max;
    int *min;
    int t_max;
    int t_min;
    int strand;
    char frame[100];
    char c_win[100];
    int id;
    int cons_id;
    int linewidth;
    char colour1[30]; /* forward */
    char colour2[30]; /* reverse */
    ruler_s *ruler;
} obj_reading_coverage;

int reading_coverage_reg(GapIO *io, Tcl_Interp *interp, char *frame,char 
			 *rcov_win, int cons_id, ruler_s *ruler, int strand);

#endif
