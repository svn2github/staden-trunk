#ifndef _RESTRICTION_ENZYME_H_
#define _RESTRICTION_ENZYME_H_

#include <tcl.h>
#include <tk.h>

#include "renz_utils.h"
#include "tg_gio.h"
#include "template_display.h"

#define TASK_RENZ_INFO   0

typedef struct s_match_t {
    tg_rec contig;
    int length;
    R_Match *match;
    int num_match;
} s_match;

typedef struct obj_renz_t { 
    Tcl_Interp *interp;
    int num_enzymes;
    R_Enz *r_enzyme;
    int num_match;
    s_match c_match;  /* matches for contig */
    int start;
    int end;
    int sequence_type;
    int yoffset;
    tick_s *tick;
    cursor_s xhair;
    int text_offset;
    char *text_colour;

    int id;
    char window[100];
    char names_win[100];
    char frame[100];
    ruler_s *ruler;
    win **win_list;
    int num_wins;
    WorldPtr *world;
    CanvasPtr *canvas;
    StackPtr *zoom;
    cursor_t *cursor;
    int cursor_visible;
    int max_overlap;
} obj_renz;

/* restriction enzyme structure */
typedef struct obj_t_renz_t {
    Tcl_Interp *interp;
    int num_enzymes;
    R_Enz *r_enzyme;
    int num_match;
    R_Match *match; 
    s_match *c_match;  /* matches for each contig */
    char frame[100];
    char window[100];
    GapIO *io;
    int sequence_type;
    int yoffset;
    tick_s *tick;
    int text_offset;
    int template_id;
    int num_contigs;
    int start;
    int end;
} obj_t_renz;

int r_enz_1d_reg(GapIO *io, 
		 tg_rec contig, 
		 char *disp_script, 
		 char *quit_script, 
		 char *window); 

int
template_renz_reg(Tcl_Interp *interp,                                  /* in */
		  GapIO *io,                                           /* in */
		  int *contig_array,
		  int num_contigs,
		  char *filename,                                      /* in */
		  char *frame,                                         /* in */
		  char *plot,                                          /* in */
		  int template_id,                                     /* in */
		  char *inlist,                                        /* in */
		  int num_items,                                       /* in */
		  tick_s *tick,                                        /* in */
		  int yoffset);                                        /* in */

int renz_reg(Tcl_Interp *interp, GapIO *io, char *filename, char *frame,
	     char *names_win, char *re_win, char *inlist, int num_items,
	     tg_rec contig_num, int lreg, int rreg, int text_offset,
	     char *text_fill, tick_s *tick, int yoffset, ruler_s *ruler,
	     cursor_s xhair);

int
plot_multi_r_enz(GapIO *io,
		 char *filename,
		 char *win_name,
		 char *plot,
		 char *inlist,
		 int list_size,
		 tg_rec contig_num,
		 int start,
		 int end,
		 int text_offset,
		 int yoffset,
		 int tick_ht,
		 int tick_wd);

int
PrintEnzymeByEnzyme(R_Enz *r_enzyme,                                   /* in */
		    R_Match *match,                                    /* in */
		    int total_matches,                                 /* in */
		    int num_enzymes,                                   /* in */
		    char *sequence,                                    /* in */
		    int sequence_len,                                  /* in */
		    int sequence_type,                                 /* in */
		    int lreg,
		    int do_all);                                         /* in */

int
OrderOnPosition (R_Enz *r_enzyme,
		 R_Match *match,
		 int total_matches,
		 char *sequence,
		 int sequence_len,
		 int sequence_type,
		 int lreg);

int
Create_REnz_Tags(GapIO *io,
		 tg_rec contig_num,
		 obj_renz *r,
		 char *enz_list,
		 char **enz_ids,
		 int num_ids);

int REnz_Init(Tcl_Interp *interp);

#endif
