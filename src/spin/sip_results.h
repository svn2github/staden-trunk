#ifndef _SIP_RESULT_H_
#define _SIP_RESULT_H_

#include "array.h"
#include "seq_results.h"
#include "seq_reg.h"
#include "seq_raster.h"

#ifdef REMOVE
typedef struct sipresult {
    void (*op_func)(int seq_num,
		     void *obj,
		     seq_reg_data *data);
    void (*pr_func)(void *obj);
    void (*txt_func)(void *obj);
    Array ptarray;
    int n_pts;
    void *input;
    void *output;
    int index;
    int seq[2];
    int replot;
    double sf_m;
    double sf_c;
} SipResult;
#endif

typedef struct matname {
    char *name;
    int **matrix;
} mat_name;

void SipFreeResult(seq_result *result);

#ifdef REMOVE
/*
 * returns the nearest match coords to position start
 */ 
d_point FindNearestMatch(seq_result *result, d_point start, double x_scale);
d_point FindNearestLine(seq_result *result, d_point start, double x_scale);
#endif

/* setting and getting global variables */
void set_replot_temp(int response);
int get_replot_temp(void);
void set_max_matches(int num);
int get_max_matches(void);
void set_def_matches(int num);
int get_def_matches(void);

int get_remove_dup(void);
void set_remove_dup(int yn);

/*
 * score, probabilites, expected and observed scores output in the
 * text output window
 */
int CalcProbs(seq_result *result, int span_length, int min_score);
int CalcIdentityProbs(seq_result *result, int word_len);

int set_matrix_file(char *file, int type);
int ** get_matrix_file(int type);
char *get_matrix_name(int type);
int set_matrix_identity(int type);
int **get_matrix_identity(int type);

void DisplaySequences(Tcl_Interp *interp, char *raster_win, int result_index, 
		      int seq_id_h, int seq_id_v);
void DestroyDisplaySequences(Tcl_Interp *interp, int result_index);

void SipRasterPlotFunc(Tk_Raster *raster, char *raster_win, int job,
		       int x0, int y0, int x1, int y1);

void SipUpdateRasterWindow(Tcl_Interp *interp, char *raster_old,
			   char *raster_new, int new_id, int old_id, int job);
void SipUpdateResultWindow(Tcl_Interp *interp, char *raster_old,
			   char *raster_new, int old_id, int new_id,
			   int result_id, int job);
int SipReSetRasterWindowSize(Tcl_Interp *interp, char *raster_win);
int sip_find_result(char *raster_win, int seq_id_h, int seq_id_v);
#endif
