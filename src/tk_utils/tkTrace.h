#ifndef _TK_TRACE_H
#define _TK_TRACE_H

#include <io_lib/Read.h>

#include "tk.h"
#include "postscript.h"
#include "cli_arg.h"

/* Pseudo functions... */
/* Trace *t, double point */
#define point_to_pixel(t, point) \
    ((int)((point) * (t)->scale_x) - (int)((t)->disp_offset * (t)->scale_x))

#define pixel_to_point(t, pixel) \
    (((int)((t)->disp_offset * (t)->scale_x) + (pixel)) / (t)->scale_x)

/* Positional information */
typedef struct {
    int bd; /* border width */
    int y;  /* Y coordinate */
    int h;  /* height */
} Trace_pos;

/* Postscript information */
typedef struct {
    FILE		*ps_file;
    LINE_OPTIONS	options_A, options_C, options_G, options_T, options_N;
    double		scale_y, scale_x;
    int			*basePos_lookup;
    int			trace_height, seq_y_pos, num_y_pos;
    int			first_base, last_base;
    char		*title;
} PS_TRACE;

/*
 * A data structure of the following type is kept for each
 * frame that currently exists for this process:
 */
typedef struct {
    Tk_Window tkwin;            /* Window that embodies the frame.  NULL
                                 * means that the window has been destroyed
                                 * but the data structures haven't yet been
                                 * cleaned up.*/
    Display *display;           /* Display containing widget.  Used, among
                                 * other things, so that resources can be
                                 * freed even after tkwin has gone away. */
    Tcl_Interp *interp;         /* Interpreter associated with
                                 * widget.  Used to delete widget
                                 * command.  */
    Tk_3DBorder border;
    int borderWidth;
    int relief;

    int width;                  /* Width to request for window.  <= 0 means
                                 * don't request any size. */
    int height;                 /* Height to request for window.  <= 0 means
                                 * don't request any size. */
    int flags;

    int cursor_pos;		/* Cursor position in base coordinates */
    Read *read;
    XColor *Acol, *Ccol, *Gcol, *Tcol;
    XColor *CursorCol;
    XColor *CutoffCol;
    XColor *VectorCol;
    XColor *ConfCol;
    XColor *ConfNegCol;
    GC Agc, Cgc, Ggc, Tgc;
    GC CursorGC;
    GC CutoffGC;
    GC VectorGC;
    GC ConfGC;
    GC ConfNegGC;
    Pixmap stipple;
    GC QualVecGC;
    int disp_offset;		/* Sample coordinates */
    int disp_offset_old;	/* Sample coordinates */
    int disp_width;		/* Sample coordinates */
    double scale_y;		/* Multiplies the y coordinates */
    double scale_x;		/* Number of pixels per sample */
    double scale_conf;		/* Confidence plot scaling */
    char *xScrollCmd;

    uint_2 *tracePos;		/* Index of each sample in the base array */
    uint_2 *tracePosE;		/* Index of each sample in the base array */

    Trace_pos pos[4];
    int show_numbers;
    int show_sequence;
    int show_edits;
    int show_trace;
    int show_conf;

    Tk_Font font;		/* Sequence + base number font */
    Tk_FontMetrics fm;
    int font_width;

    Tk_Font conf_font;		/* Confidence number font */
    int conf_font_width;
    int conf_font_height;

    int Ned;			/* Number of edited bases */
    int MaxNed;			/* Maximum number of edited bases */
    char *edBases;		/* Our new sequence */
    int_2 *edPos;		/* An array of length 'read->NBases' showing
				 * which original base corresponds to which
				 * edited base.
				 */

    Pixmap pn;			/* The 'numbers' pixmap */
    Pixmap ps;			/* The 'sequence' pixmap */
    Pixmap pe;			/* The 'edit' pixmap */
    Pixmap pt;			/* The 'trace' pixmap */

    int comp;

    int leftVector;		/* The location of various vector clips. */
    int rightVector;		/*   Set to -1 when none are found.      */

    int show_ends;		/* Whether to allow scrolling off the ends
				 * of the trace. Used to allow the editor to
				 * force any part of the trace (even an end)
				 * to be centred in the screen.
				 */

    int line_width;             /* Width of traces */
    int cursor_pos_old;		/* Old cursor position (for fast redrawing) */

    int1 *edConf;		/* 1 deep confidence values (for exp files) */

    int last_disp_offset;	/* For optimisation of redraw code */
    int last_disp_width;	/* For optimisation of redraw code */
    int last_disp_height[4];	/* For optimisation of redraw code */
    int visible;		/* Visibility mode, for optimisation */
    int old_scroll_mode;	/* Style of scrollbar callback arguments */

    PS_OPTIONS		ps_options;
    PS_TRACE		ps_trace;

    int trace_scale;		/* Height in trace coordinates. zero => auto */
    int style;			/* 0=chromtogram
				 * 1=filled chromatogram (not very useful)
				 * 2=historgram (pyro)
				 * 3=stick chart (Solexa etc)
				 */
    int yticks;			/* Horizontal lines in plot. zero => none */

} DNATrace;

/* Flags */
/* Redraw only the borders */
#define TRACE_BORDER  (1<<0)

/* Redraw all of text and trace widgets, but not borders */
#define TRACE_REDRAW  (1<<1)

/* Only redraw the bits changed by scrolling */
#define TRACE_SCROLL  (1<<2)

/* Forces a redraw of everything */
#define TRACE_RESCALE (1<<3)

/* A display request has been registered. */
#define TRACE_WAITING (1<<4)

/* The cursor needs redrawing */
#define TRACE_CURSOR  (1<<5)

#define TRACEP_N 0
#define TRACEP_S 1
#define TRACEP_E 2
/* TRACEP_T should be the last defined */
#define TRACEP_T 3

/* Maximum number of edits */
#define TRACE_EDITS 200

/* Styles */
#define STYLE_CHROMA 0
#define STYLE_FILLED 1
#define STYLE_PYRO   2
#define STYLE_STICK  3

/* ---- Prototypes ---- */
extern int pixel_to_base(DNATrace *tracePtr, int pixel, int allow_end);

extern void trace_draw_trace(DNATrace *t, Display *d, Pixmap p,
			     int x0, int x1, int yoff, int height);

extern void trace_draw_sequence(DNATrace *t, Display *d, Pixmap p,
				int x0, int x1, int yoff, int height);

extern void trace_draw_numbers(DNATrace *t, Display *d, Pixmap p,
			       int x0, int x1, int yoff, int height);

extern int trace_get_pos(DNATrace *t, int ind);

extern void trace_draw_edits(DNATrace *t, Display *d, Pixmap p,
			     int x0, int x1, int yoff, int height);

extern void complement_trace(DNATrace *t); /* from tkTraceComp.c */

extern void trace_flash(DNATrace *t);

int    	*trace_index_to_basePos(uint_2 *basePos, int NBases, int NPoints);
int   visible_region(DNATrace *t);

#endif /* _TK_TRACE_H */
