#ifndef _CANVAS_BOX_H
#define _CANVAS_BOX_H

#include <math.h>
#include <tcl.h>
#include <inttypes.h>

#define TASK_CANVAS_SCROLLX         1000
#define TASK_CANVAS_SCROLLY         1001
#define TASK_CANVAS_ZOOMBACK        1002
#define TASK_CANVAS_ZOOM            1003
#define TASK_CANVAS_CURSOR_X        1004
#define TASK_CANVAS_CURSOR_Y        1005
#define TASK_CANVAS_CURSOR_DELETE   1006
#define TASK_CANVAS_RESIZE          1007
#define TASK_CANVAS_REDRAW          1008
#define TASK_CANVAS_WORLD           1009
#define TASK_WINDOW_ADD             1010
#define TASK_WINDOW_DELETE          1011
#define TASK_DISPLAY_TICKS          1012
#define TASK_DISPLAY_RULER          1013
#define TASK_CONS_WORLD             1014
#define TASK_CONS_JOIN              1015
#define TASK_CONS_CURSOR_DELETE     1016
#define TASK_CONS_ID                1017

#define MAX_NUM_WINS 11
/* HACK: should move to os.h - can use rint on unix machines */
#define ROUND(x)	((x) < 0.0 ? ceil((x) - 0.5) : floor((x) + 0.5))

typedef struct line {
    int x1;
    int x2;
    int y1;
    int y2;
} Line;

typedef struct plotrec {
    Line l;		/* line structure */
    int num;		/* index number */
    char *type;	        /* line tags */
    char *colour;	/* display colour */
    char arrow[6];      /* arrow type: none, first, last, both */
} PlotRec;

/* 'Double' version of Line */
typedef struct dline {
    double x1;
    double x2;
    double y1;
    double y2;
} DLine;

/* 'Double' version of PlotRec */
typedef struct dplotrec {
    DLine l;		/* line structure */
    int num;		/* index number */
    char *type;	        /* line tags */
    char *colour;	/* display colour */
    char arrow[6];      /* arrow type: none, first, last, both */
} DPlotRec;

typedef struct coffset {
    int offset;
    int gap;
} c_offset;

typedef struct box_s{
    int x1;
    int y1;
    int x2;
    int y2;
} box;

typedef struct d_box_s{
    double x1;
    double y1;
    double x2;
    double y2;
} d_box;

typedef struct cir_s_ {
    int x1;
    int y1;
    int diameter;
} cir_s;

typedef struct s_zoom_s {
    box *zoom;
    float amount;
    char scroll;
    int r_id; /* result id */
} s_zoom;

typedef struct stack_s StackPtr;

typedef struct stack_s {
    d_box *data;
    StackPtr *next;
} Stack;


typedef struct WorldPtr_s {
    d_box *visible;
    d_box *total;
} WorldPtr;

typedef struct CanvasPtr_s {
    int width;
    int height;
    double ax;
    double ay;
    double bx;
    double by;
    int64_t x;    /* canvas pixel value of left hand side */
    int64_t y;
} CanvasPtr;

typedef struct win_s {
    WorldPtr *world;
    CanvasPtr *canvas;
    StackPtr *zoom;
    char *window;
    char scroll;
    int id;
} win;

typedef struct tag_t {
    int width;
    int offset;
} tag_s;

typedef struct cursor_t {
    int width;
    char *colour;
} cursor_s;

typedef struct tick_t {
    int line_width;
    int ht;
    char *colour;
} tick_s;

typedef struct rtick_t {
    tick_s t;
    int offset;
    int num;
} rtick_s;

typedef struct ruler_t {
    rtick_s tick;
    char *window;
    int offset;
    char *colour;
    int line_width;
    tag_s tag;
    int start;
    int end;
} ruler_s;

ruler_s *ruler_struct(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *win_name, 
		      int tags);

void
printCanvas(CanvasPtr *canvas);

void
initCanvas(Tcl_Interp *interp,
	   CanvasPtr *canvas,        /* in, out */
	   char *window);

void
CanvasToWorld(CanvasPtr *canvas,
	      int64_t cx, int64_t cy,
	      double *wx, double *wy);

void
WorldToCanvas(CanvasPtr *canvas,
	      double wx, double wy, 
	      double *cx, double *cy);

void
SetCanvasCoords(Tcl_Interp *interp,
		double x1, double y1, double x2, double y2,
		CanvasPtr *canvas);  /* out */

void scrollRegion(Tcl_Interp *interp, win **win_list, int num_wins,
	     d_box *total, CanvasPtr *canvas);

void scaleCanvas(Tcl_Interp *interp, win **win_list, int num_wins,
	    char *tags, d_box *current, CanvasPtr *canvas);

void scaleSingleCanvas(Tcl_Interp *interp, WorldPtr *world, CanvasPtr *canvas,
		  char *window, char scroll, char *tag);

void resizeCanvas(Tcl_Interp *interp, char *window, win **win_list,
	     int num_wins, d_box *visible, d_box *total, CanvasPtr *canvas);

void
createZoom(StackPtr **zoom);

d_box *
examineZoom(StackPtr *zoom);

void
popZoom(StackPtr **zoom);

void
pushZoom(StackPtr **zoom, d_box *visible);

void freeZoom(StackPtr **zoom);

void canvasZoom(Tcl_Interp *interp, CanvasPtr *canvas, char *window,
	   WorldPtr *world, win **win_list, int num_wins,
	   StackPtr **zoom_list, box *zoom, char scroll);

void canvasZoomback(Tcl_Interp *interp, CanvasPtr *canvas, char *window,
	       WorldPtr *world, win **win_list, int num_wins,
		    StackPtr **zoom_list);
tag_s
tag_struct(Tcl_Interp *interp,
	   Tcl_Obj *defs_ptr,
	   char *win_name,
	   int width,
	   int offset);

cursor_s
cursor_struct(Tcl_Interp *interp,
	      Tcl_Obj *defs_ptr,
	      char *win_name,
	      int width,
	      char *colour);

tick_s *tick_struct(Tcl_Interp *interp, Tcl_Obj *defs_ptr, char *win_name,
		    int width, int height, char *colour);

int addWindow(win **win_list, int *num_wins, char *window, char scroll,
	      int id);

void deleteWindow(win **win_list, int *num_wins, char *window);

void free_win_list(win **win_list, int num_wins);

void canvasCursorX(Tcl_Interp *interp, CanvasPtr *canvas, char *frame,
	      char *label, char *colour, int line_width, int64_t cx,
	      double wx, win **win_list, int num_wins);

void canvasCursorY(Tcl_Interp *interp, CanvasPtr *canvas, char *frame,
	      char *label, char *colour, int line_width, int64_t cy,
	      double wy, win **win_list, int num_wins);

void
canvasScrollX(Tcl_Interp *interp,
	      char *window,
	      win **win_list,
	      int num_wins,
	      d_box *visible,
	      CanvasPtr *canvas,
	      char *scroll_args);

void
canvasScrollY(Tcl_Interp *interp,
	      char *window,
	      win **win_list,
	      int num_wins,
	      d_box *visible,
	      CanvasPtr *canvas,
	      char *scroll_args);

void
freeZoom(StackPtr **zoom);

int
lengthZoom(StackPtr *zoom);
void copyZoom(StackPtr **zoom2, StackPtr *zoom1);

void draw_single_ruler(Tcl_Interp *interp, ruler_s *ruler, CanvasPtr *canvas,
		       double start, double end, int disp_ticks);

void draw_single_ruler_vertical(Tcl_Interp *interp, ruler_s *ruler,
				CanvasPtr *canvas, double start, double end,
				int disp_ticks);
double
calc_zoom_origin(double min1, double min2, double max1, double max2);
double
calc_zoom_sf(double min1, double min2, double max1, double max2);

int canvasx(Tcl_Interp *interp, char *window, double val);
int canvasy(Tcl_Interp *interp, char *window, double val);
void free_ruler_struct(ruler_s *r);

#endif


