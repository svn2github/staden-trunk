#ifndef _SEQ_RESULTS_H_
#define _SEQ_RESULTS_H_
#include <tcl.h>

#include "sequence_formats.h"
#include "seq_reg.h"

#define SEQ_TYPE_STRINGSEARCH   (1<<0)
#define SEQ_TYPE_RESTRICTION    (1<<1)
#define SEQ_TYPE_BASECOMP       (1<<2)
#define SEQ_TYPE_CODONPREF      (1<<3)
#define SEQ_TYPE_AUTHOR         (1<<4)
#define SEQ_TYPE_BASEBIAS       (1<<5)
#define SEQ_TYPE_TRNA           (1<<6)
#define SEQ_TYPE_STOPCODON      (1<<7)
#define SEQ_TYPE_STARTCODON     (1<<8)
#define SEQ_TYPE_SPLICE         (1<<9)
#define SEQ_TYPE_WTMATRIXSEARCH (1<<10)
#define SEQ_TYPE_GRAPH_PLOT     (1<<11) /* emboss plots */
#define SEQ_TYPE_DOT_PLOT       (1<<12) /* sip/emboss plots */
#define SEQ_TYPE_RULER          (1<<13)

/* compound definitions */
#define SEQ_TYPE_GENESEARCH (SEQ_TYPE_CODONPREF | SEQ_TYPE_AUTHOR)

/* data type definitions, stored in seq_result->graph */
#define SEQ_GRAPH   0  /* r_graph */
#define SEQ_DOT     1  /* d_plot */
#define SEQ_STICK   2  /* stick */
#define SEQ_GENE    3  /* gene_search */
#define SEQ_RENZ    4
#define SEQ_E_DOT   5  /* emboss e_graph */

#define MAXLINE 1024

/* result information jobs */
#define INPUT  0
#define OUTPUT 1
#define NPTS   2
#define INDEX  3
#define RESULT 4

/* sequence structure */
#define LINEAR 0
#define CIRCULAR 1

/* emboss object types */
#define E_RECTANGLE     0
#define E_RECTANGLEFILL 1
#define E_TEXT          2
#define E_LINE          3

typedef struct seq_info {
    int library;        /* library index */
    int seq_len;        /* length */
    int type;           /* 1 = DNA, 2 = PROTEIN */
    int seq_structure;  /* 0 = linear, 1 = circular */
    int seq_id;         /* unique identifier */
    int count;          /* number of sub-sequences */
    char *sequence;     /* sequence */
    char *name;         /* name */
    char *raster;       /* default raster associated with sequences */
} SeqInfo;

typedef struct subseq_ {
    SeqInfo *seq;
    int start;
    int end;
    int seq_id;
    char *name;
    char *identifier;
    Featcds **key_index;
} SubSeq;

typedef struct range_ {
  int start;
  int end;
} range;

typedef struct d_line_ {
    double x0;
    double y0;
    double x1;
    double y1;
} d_line;

typedef struct p_score_ {
    int pos;
    double score;
} p_score;

typedef struct a_score_ {
    p_score *p_array;
    int n_pts;
    d_line dim;
} a_score;

typedef struct stick_ {
    a_score *ap_array;
    int n_pts;
} stick;

typedef struct r_graph_ {
    p_score *p_array;
    int n_pts;
    d_line dim;
} r_graph;

typedef struct gene_search_ {
    p_score *p_array;
    int n_pts;
    d_line dim;
    char *top;
} gene_search;

typedef struct d_point_ {
    int x;
    double y;
} d_point;

typedef struct pt_score_ {
    int x;
    int y;
    int score;
} pt_score;

typedef struct d_plot_ {
    pt_score *p_array;
    int n_pts;
    d_line dim;
    int win_len;
} d_plot;

/* emboss graph structure */
typedef struct e_obj_ {
    int type;
    d_line pos;
    char colour[20];
} e_obj;

typedef struct e_graph_ {
    p_score *p_array;
    int n_pts;
    d_line dim;
    int n_data_obj;
    e_obj *d_obj;
    int n_graph_obj;
    e_obj *g_obj;
    char *title;
    char *maintitle;
    char *subtitle;
    char *xtitle;
    char *ytitle;
} e_graph;

typedef struct seq_result_ {
    void (*op_func)(int seq_num, void *obj, seq_reg_data *data);
    void (*pr_func)(void *obj, seq_reg_plot *plot);
    void (*txt_func)(void *obj);
    void *data;         /* data to plot result */
    void *input;        /* input parameters */
    void *output;       /* general info about plotting the result */
    int id;             /* unique result identifier */
    int seq_id[2];      /* unique sequence identifier */
    int type;           /* unique name for each result type */
    int frame;          /* reading frame of a result or 0 if not applicable */
    void *text_data;    /* specific data required for textual output */
    int graph;          /* defines the data type */
} seq_result;

/*
 * add a new sequence to sequence list
 * add a new sequence to registration scheme
 * update the sequence list box
 */
int AddSequence(Tcl_Interp *interp, int direction, int library, char *entry,
		char *sequence, int seq_structure, int seq_type, 
		Featcds **key_index, char *identifier);

/*
 * delete a sequence from sequence list
 * delete a sequence from registration scheme
 * update the sequence list box
 */
void DeleteSequence(Tcl_Interp *interp, int seq_num);

/*
 * create or extend array containing the sequence info created each time a
 * new sequence is read in using the 'get sequences', 'get horizontal seq'
 * or 'get vertical seq' options
 * returns the next free index
 */
int SeqCreate(void);

void Delete_Seq(int seq_num);
int GetSeqId(int seq_num);
int NumSequences(void);

/*
 * set of functions that return info on a particular sequence in the global 
 * seqs array
 */

int GetSeqNum(int id);

/* added for parsing features */
char *GetSeqIdentifier(int seq_num);
char *GetSeqCdsExpr(int seq_num, int i);
char *GetSeqKeyIndexCds(int seq_num, int idx);
int GetSeqNumberCds(int seq_num);
char *GetSeqProteinId(int seq_num, int k);
Featcds **GetSeqKeyIndex(int seq_num);
/* added for parsing features */

char *GetSeqName(int seq_num);
char *GetSeqBaseName(int seq_num);
char *GetSeqSequence(int seq_num);
int GetSeqLength(int seq_num);
int GetSeqType(int seq_num);
int GetSeqStructure(int seq_num);
void SetSeqStructure(int seq_num, int seq_structure);
char *GetRaster(int seq_num);
void SetRaster(int seq_num, char *raster);
int GetSeqLibrary(int seq_num);
char *GetSeqLibraryName(int seq_num);
int GetSeqDirection(int seq_num);
int GetParentalSeqId(int seq_num);
char *GetParentalSeqName(int seq_num);

/*
 * set the "active" sequence in a given direction
 */
int Set_Active_Seq(int seq_index, int direction);

int GetActiveSeqNumber(int direction);

int CopyRange(Tcl_Interp *interp, int seq_id, int start, int end);
int ComplementSeq(Tcl_Interp *interp, int index);

int TranslateSeq(Tcl_Interp *interp, int index, int rf, int start, int end);

int ScrambleSeq(Tcl_Interp *interp, int index);
int RotateSeq(Tcl_Interp *interp, int index, int origin);
int RnaSeq(Tcl_Interp *interp, int index);
int GetSubSeqStart(int seq_num);
int GetSubSeqEnd(int seq_num);
int GetSubSeqLength(int seq_num);
int TranslateTogether(Tcl_Interp *interp, int index);

int SetRange(Tcl_Interp *interp, int seq_id, int start, int end);

int GetSeqIdFromName(char *name);
seq_result *seq_id_to_result(int result_id);
void SeqRasterPlotFunc(Tk_Raster *raster, char *raster_win,int job,
		       int x0, int y0, int x1, int y1);
void SeqSuperimposeResult(Tcl_Interp *interp, char *raster_win, int result_id,
			  double o_wx0, double o_wy0, double o_wx1,
			  double o_wy1);
int comparison2(void *v_result, int type);
d_point FindNearestMatch(seq_result *result, d_point start, double x_scale);
d_point FindNearestLine(seq_result *result, d_point start, double x_scale);
int seq_find_result(char *raster_win, int seq_id_h, int seq_id_v);

#endif


