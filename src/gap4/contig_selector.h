#ifndef _CONTIG_SELECT_H_
#define _CONTIG_SELECT_H_

#include "cs-object.h"
#include "canvas_box.h"

#define TASK_CS_REDRAW       0

extern HTablePtr csplot_hash[HASHMODULUS];

typedef struct {
    int buffer_count;
    int do_update;
    char hori[100];
    char vert[100];
    tag_s tag;
    tick_s *tick;
    cursor_s cursor;
    int line_width;
    char *line_colour;

    char frame[100];
    char window[100];
    int spare; /* RE-USE this */
    win **win_list;
    int num_wins;
    WorldPtr *world;
    CanvasPtr *canvas;
    StackPtr *zoom;
} obj_cs;


/* determines the position of a base in terms of the entire database */
int
find_position_in_DB(GapIO *io, 
		 int c_num, 
		 int position);

void
PlotRepeats(GapIO *io, 
	    mobj_repeat *repeat);

/* draw the contig lines of the contig selector */
int
Draw_CS_Contigs(Tcl_Interp *interp,                                    /* in */
		GapIO *io,                                             /* in */
		char *win_name,                                        /* in */
		char *colour,                                          /* in */
		int width,                                             /* in */
		int tick_wd,                                           /* in */
		int tick_ht,                                           /* in */
		int offset,                                            /* in */
		char *direction);                                      /* in */

void
update_contig_order(Tcl_Interp *interp,
		    GapIO *io,
		    int cs_id,
		    int *order_array,
		    int num_contigs,
		    int cx);

void
ReOrder(GapIO *io, 
	GCardinal *order,
	int c_from,
	int c_to); 

void
move_matches(int item_num);

void
plot_cs_tags(Tcl_Interp *interp,                                       /* in */
	     GapIO *io,                                                /* in */
	     char *win_name,                                           /* in */
	     int width,                                                /* in */
	     int offset);                                              /* in */
/*
 * variation on vtagget. In this case, the tag_num is returned 
 */
GAnnotations *
get_tag_num(GapIO *io,                                                 /* in */
	    int gel,                                                   /* in */
	    int num_t,                                                 /* in */
	    char **type,                                               /* in */
	    int *tag_num);                                            /* out */
int
contig_selector_reg(Tcl_Interp *interp, 
		    GapIO *io, 
		    char *frame,
		    char *csh_win,
		    tag_s tag,
		    cursor_s cursor,
		    tick_s *tick);
int
contig_comparator_reg(Tcl_Interp *interp, 
		      GapIO *io, 
		      obj_cs *cs,
		      char *csp_win,
		      char *csv_win);

#endif
