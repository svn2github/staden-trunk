#ifndef _GAP_DISPLAY_H
#define _GAP_DISPLAY_H

#include <tk.h>
#include "canvas_box.h"
#include "io-reg.h"
#include "hash.h"
#include "gap-dbstruct.h"
#include "IO.h"
#include "list.h"
#include "template.h"
#include "ruler_display.h"

#define FALSE 0
#define TRUE 1

#define ONE_READING (1<<0)         /* one reading on the template */
#define DIFF_CONTIGS (1<<1)        /* readings on different contigs */
#define FORW_REV_READINGS (1<<2)   /* both for/rev readings for single template */
#define CONTRADICTORY (1<<3)       /* readings on the template are contradictory */
#define SPAN_CONTIG     (1<<4)     /* template spans two contigs used in plot */
#define SPAN_CONTIG_INCONS (1<<5)  /* templates spans 2 contigs but inconsistent */

#define TASK_TEMPLATE_REDRAW   0

/* configure array for what to display on the template display */
#define NUM_CONFIGS        9
#define TEMPLATES          0
#define READINGS           1
#define MULTI_TEMPLATES    2
#define READ_PAIRS         3   
#define RULER              4
#define TICKS              5
#define SPAN_READ_PAIRS    6
#define CALC_CONTIG_POS    7
#define CONSIST_READ_PAIRS 8

typedef struct {
    list_t *gel_cont;
    int consistency;
    int flags;
} td_template;

/* Template positions */
typedef struct template_ps {
    template_c *t;
    int contig;
    int t_num;
    int diff;
    int length;
    int start;
    int end;
    int consist;
    int num_r;
} template_p;

/* Offsets and distances between contigs? */
typedef struct template_os {
    int gap;
    int cnt;
    float average;
} template_o;

/* Temporary struct for read-pair/template data */
typedef struct template_ds {
    int start;
    int end;
    int diff;
    int consist;
    int readpair;
} template_d;

#define SCRIPT_LEN 1024
typedef struct {
    Tcl_Interp *interp;
    c_offset *contig_offset;
    int *contig;
    int num_contigs;
    char frame[100];
    char window[100];
    char t_win[100];
    int id;
    ruler_s *ruler;
    cursor_s xhair;
    win **win_list;
    int num_wins;
    WorldPtr *world;
    CanvasPtr *canvas;
    StackPtr *zoom;
    PlotRec *readings;
    int num_readings;
    PlotRec *ruler_coord;
    template_c **tarr;
    int configs[NUM_CONFIGS];
    int line_width;
    int line_bold;
    int do_update;
    int buffer_count;
    cursor_t **cursor;
    int *cursor_visible;
} obj_template_disp;


/*****************************************************************************/
/*                               String2Int                                  */
/*****************************************************************************/
/* converts tcl string argument into a C integer */
int
String2Int(Tcl_Interp *interp, char *arg, int *num);

/* Convert 'position' values of contig, templates and readings to pixels.    */
/*****************************************************************************/
/*                                  CalcXFactor                              */
/*****************************************************************************/
float CalcXFactor(int max_len); 

void CalcYDepth(int num, PlotRec *PArray, int max_depth, int *num_levels);

/*****************************************************************************/
/*                             CalcContigs                                   */
/*****************************************************************************/
void CalcContigs(GapIO *io, int contig_num, int max_len, int ruler_offset,
		 char *colour, PlotRec CArray[]);

void CalcContigTicks(int max_len, int x1, int x2, int tick_dist, 
		     float *tick_interval, int *num_ticks);

void CalcTags(GapIO *io, int max_len, int item_num, int tag_db_count,
	      PlotRec *ItemArray[], int *ann_index,
	      GAnnotations AnnotationsArray[], PlotRec TagArray[],
	      int *num_tags);

/*****************************************************************************/
/*                             CalcTotalContigLen                            */
/*****************************************************************************/
/* return the total length of all the contigs in a database */
int64_t
CalcTotalContigLen(GapIO *io);                                         /* in */

/*****************************************************************************/
/*                             CalcTotalReadingLen                           */
/*****************************************************************************/
/* return the total length of all the readings in a database */
int64_t
CalcTotalReadingLen (GapIO *io, int ngels); 

/*****************************************************************************/
/*                             CalcLongContig                                */
/*****************************************************************************/
/* return the contig num of the longest contig in the database */
int64_t
CalcLongContig(GapIO *io);                                             /* in */

/*****************************************************************************/
/*                             CalcContigNumLen                              */
/*****************************************************************************/
/* find the contig length & number of a contig containing reading_num */
void CalcNumLenContig(GapIO *io, int reading_num, int *contig_num,
		      int *contig_len);

void CalcConsTags(GapIO *io, int max_len, int c_num, int *ann_index,
		  GAnnotations AnnotationsArray[], PlotRec ConsArray[],
		  int *num_tags);

/*
 * Registers the template display.
 */
int template_reg(Tcl_Interp *interp, GapIO *io, int *contig_array, 
		 int num_contigs, char *frame, char *t_win, char *r_win, 
		 ruler_s *ruler, cursor_s xhair, int line_width, int line_bold);

int update_template_display(Tcl_Interp *interp, GapIO *io, 
			    obj_template_disp *t, int recalculate);

int display_reading_tags(Tcl_Interp *interp, GapIO *io, obj_template_disp *t);

double TemplateLocalCursor(int id, c_offset *contig_offset, int *contig_array,
		    int num_contigs, double wx);

void template_update_cursors(GapIO *io, obj_template_disp *t, int show_cursor);

int FindTemplatePositions(GapIO *io, c_offset *contig_offset,
			  int *contig_array, int num_contigs,
			  template_c **tarr, template_d **t_changes);

void update_template_contig_order(Tcl_Interp *interp, GapIO *io,
				  int template_id, int cx, int *contig_list,
				  int num_contigs);

void refresh_contig_order(Tcl_Interp *interp, GapIO *io, int template_id);

void plot_lines(Tcl_Interp *interp, PlotRec *array, int num, char *win_name,
		int line_width);

int template_find_left_position(GapIO *io, int *contig_array, int num_contigs,
				c_offset *contig_offset, double wx);

void SetReadingPosLen(int whole_reading, GapIO  *io, int reading_num, 
		      int *r_pos, int *r_len); 

void CalcXCoords(int pos, int len, int *x1, int *x2);

double TemplateLocalCursor(int id, c_offset *contig_offset, int *contig_array,
			   int num_contigs, double wx);

/*
 * Set status (corresponds to colour in template display).
 * The order here counts as it sets the precedence for which
 * colour should be displayed.
 */
int getStatus(template_c *t);

int
inContigList(int *contig_array,
	     int num_contigs,
	     int contig);

#endif /* _GAP_DISPLAY_H */
