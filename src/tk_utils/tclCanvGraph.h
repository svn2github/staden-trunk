#ifndef _TCLCANVGRAPH_H
#define _TCLCANVGRAPH_H

#include <tcl.h>

#define G_DOT           1  /* single dots */
#define G_LINE          2  /* contiguous line */
#define G_LINES         3  /* discrete lines */
#define G_RECTANGLE     4
#define G_RECTANGLEFILL 5
#define G_TEXT          6

typedef struct gd_line_ {
    double x0;
    double y0;
    double x1;
    double y1;
} gd_line;

typedef struct g_pt_ {
    double x;
    double y;
} g_pt;

typedef struct darray_ {
    gd_line *d_array;
    int n_dlines;
    int type;
} darray;

typedef struct parray_ {
    g_pt *p_array;
    int n_pts;
    int type;
} parray;

typedef struct Graph_s {
    darray *d_arrays; /* discrete line array */
    int n_darrays;
    parray *p_arrays; /* discrete points array */
    int n_parrays;
    gd_line dim;     /* dimensions of data range */
} Graph;

void Tcl_GraphInit(Tcl_Interp *interp);

Tcl_Obj *Tcl_NewGraphObj (Graph *grPtr);
int Tcl_SetGraphObj(Tcl_Obj *objPtr, Graph *grPtr);
Graph *Tcl_GetGraphFromObj(Tcl_Obj *objPtr);
void PrintGraphObj(Tcl_Obj *objPtr);
int IsGraphType(Tcl_Obj *objPtr);
void PrintObj(Tcl_Obj *obj);
void Tcl_InitGraph(Graph **graph);

#endif
