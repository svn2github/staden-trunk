#include <float.h>
#include <string.h>
#include "tclCanvGraph.h"


static void
FreeGraphInternalRep (Tcl_Obj *grPtr);
static void
DupGraphInternalRep (Tcl_Obj *srcPtr, Tcl_Obj *copyPtr);
static int
SetGraphFromAny (Tcl_Interp *interp, Tcl_Obj *objPtr) ;
static void
UpdateStringOfGraph (Tcl_Obj  *grPtr);


/*
 * Type definition.
 */
static Tcl_ObjType graphType = {
    "graph",              /* name */
    FreeGraphInternalRep, /* freeIntRepProc */
    DupGraphInternalRep,  /* dupIntRepProc */
    UpdateStringOfGraph,  /* updateStringProc */
    SetGraphFromAny       /* setFromAnyProc */
};

void FreeGraphData(Graph *graph)
{
    int i;

    if (graph->n_darrays > 0) {
	for (i = 0; i < graph->n_darrays; i++) {
	    if (graph->d_arrays[i].n_dlines > 0) {
		ckfree((VOID*)graph->d_arrays[i].d_array);
	    }
	}
	ckfree((VOID*)graph->d_arrays);
    }
    if (graph->n_parrays > 0) {
	for (i = 0; i < graph->n_parrays; i++) {
	    if (graph->p_arrays[i].n_pts > 0) {
		ckfree((VOID*)graph->p_arrays[i].p_array);
	    }
	}
	ckfree((VOID*)graph->p_arrays);
    }
    ckfree((VOID*)graph);
}

/*-----------------------------------------------------------------------------
 * FreeGraphInternalRep --
 *   Free the internal representation of a graph object.
 *
 * Parameters:
 *   o grPtr - graph object being deleted.
 *-----------------------------------------------------------------------------
 */
static void
FreeGraphInternalRep (Tcl_Obj *grPtr)
{
    printf("FreeGraphInternalRep\n");

    FreeGraphData((Graph *) (grPtr)->internalRep.otherValuePtr);
}

/*-----------------------------------------------------------------------------
 * DupGraphInternalRep --
 *   Duplicate the internal representation of a graph.
 *
 * Parameters:
 *   o srcPtr - graph object to copy.
 *   o copyPtr - Target object to copy internal representation to.
 *-----------------------------------------------------------------------------
 */
static void
DupGraphInternalRep (Tcl_Obj *srcPtr, Tcl_Obj *copyPtr)
{

    printf("DupGraphInternalRep\n");

}

/*-----------------------------------------------------------------------------
 * SetGraphFromAny --
 *   Convert an object to a graph from its string representation.
 *
 * Parameters:
 *   o objPtr - Object to convert to a graph object.
 *-----------------------------------------------------------------------------
 */
static int
SetGraphFromAny (Tcl_Interp *interp, Tcl_Obj *objPtr)
{
#ifdef DEBUG
    printf("SetGraphFromAny\n");
#endif
    return TCL_OK;
}

/*-----------------------------------------------------------------------------
 * UpdateStringOfGraph --
 *    Update the string representation of a graph.
 *
 * Parameters:
 *   o objPtr - Object to convert to a graph.
 *-----------------------------------------------------------------------------
 */
static void
UpdateStringOfGraph (Tcl_Obj  *objPtr)
{
    char buffer[10];
    int len;

#ifdef DEBUG
    printf("UpdateStringOfGraph \n");
#endif
    len = 5;
    sprintf(buffer, "hello");




    objPtr->bytes = ckalloc((unsigned) len + 1);
    strcpy(objPtr->bytes, buffer);
    objPtr->length = len;

}

void Tcl_InitGraph(Graph **graph)
{
    Graph *grPtr = *graph;

    grPtr->n_darrays = 0;
    grPtr->d_arrays = NULL;
    grPtr->n_parrays = 0;
    grPtr->p_arrays = NULL;
}


/*-----------------------------------------------------------------------------
 * Tcl_NewGraphObj --
 *   Create and initialize a new graph object.
 *
 * Returns:
 *    A pointer to the object.
 *-----------------------------------------------------------------------------
 */
Tcl_Obj *Tcl_NewGraphObj (Graph *graph)
{
    Tcl_Obj *objPtr = Tcl_NewObj ();
    Graph *grPtr;
    int i, j;

    /* HACK FIXME - better way of copying over memory? */
#ifdef DEBUG
    printf("Tcl_NewGraphObj\n");
#endif
    if (NULL == (grPtr = (Graph *)ckalloc(sizeof(Graph))))
	return NULL;

    if (graph->n_darrays > 0) {
	if (NULL == (grPtr->d_arrays = (darray *)ckalloc(graph->n_darrays * sizeof(darray))))
	    return NULL;
	for (i = 0; i < graph->n_darrays; i++) {

	    if (NULL == (grPtr->d_arrays[i].d_array = (gd_line *)ckalloc(sizeof(gd_line) *
							  graph->d_arrays[i].n_dlines)))
		return NULL;
	    grPtr->d_arrays[i].n_dlines = graph->d_arrays[i].n_dlines;
	    grPtr->d_arrays[i].type = graph->d_arrays[i].type;
	}
    }

    if (graph->n_parrays > 0) {
	if (NULL == (grPtr->p_arrays = (parray *)ckalloc(graph->n_parrays * sizeof(parray))))
	    return NULL;
	for (i = 0; i < graph->n_parrays; i++) {

	    if (graph->p_arrays[i].n_pts > 0) {
		if (NULL == (grPtr->p_arrays[i].p_array = (g_pt *)ckalloc(sizeof(g_pt) *
									  graph->p_arrays[i].n_pts))) {
		    return NULL;
		}
	    }
	    grPtr->p_arrays[i].n_pts = graph->p_arrays[i].n_pts;
	    grPtr->p_arrays[i].type = graph->p_arrays[i].type;
	}
    }

    grPtr->n_darrays = graph->n_darrays;
    grPtr->n_parrays = graph->n_parrays;

    for (j = 0; j < graph->n_darrays; j++) {
	for (i = 0; i < graph->d_arrays[j].n_dlines; i++) {
	    grPtr->d_arrays[j].d_array[i].x0 = graph->d_arrays[j].d_array[i].x0;
	    grPtr->d_arrays[j].d_array[i].y0 = graph->d_arrays[j].d_array[i].y0;
	    grPtr->d_arrays[j].d_array[i].x1 = graph->d_arrays[j].d_array[i].x1;
	    grPtr->d_arrays[j].d_array[i].y1 = graph->d_arrays[j].d_array[i].y1;
	}
    }

    for (j = 0; j < graph->n_parrays; j++) {
	for (i = 0; i < graph->p_arrays[j].n_pts; i++) {
	    grPtr->p_arrays[j].p_array[i].x = graph->p_arrays[j].p_array[i].x;
	    grPtr->p_arrays[j].p_array[i].y = graph->p_arrays[j].p_array[i].y;
	}
    }

    grPtr->dim.x0 = graph->dim.x0;
    grPtr->dim.y0 = graph->dim.y0;
    grPtr->dim.x1 = graph->dim.x1;
    grPtr->dim.y1 = graph->dim.y1;

    objPtr->bytes = NULL;
    objPtr->internalRep.otherValuePtr = (VOID *) grPtr;
    objPtr->typePtr = &graphType;

    return objPtr;
}

int Tcl_SetGraphObj(Tcl_Obj *objPtr, Graph *grPtr)
{
    Tcl_InvalidateStringRep(objPtr);
    objPtr->internalRep.otherValuePtr = (VOID *)(grPtr);
    return TCL_OK;
}

Graph *Tcl_GetGraphFromObj(Tcl_Obj *objPtr)
{
    if( objPtr->typePtr != &graphType )
	SetGraphFromAny( NULL, objPtr );
    return (Graph *) (objPtr)->internalRep.otherValuePtr;
}

void Tcl_GraphInit(Tcl_Interp *interp)
{

    Tcl_RegisterObjType (&graphType);

}

void PrintGraphObj(Tcl_Obj *objPtr)
{
    Graph *grIntPtr;

    printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$PrintGraphObj\n");
    if (objPtr->typePtr != &graphType) {
	printf("not graphType\n");
    } else {
	printf("graphType\n");
    }

    grIntPtr = (Graph *) (objPtr)->internalRep.otherValuePtr;

    printf("n_darrays %d n_parrys %d dim %f %f %f %f\n",
	   grIntPtr->n_darrays, grIntPtr->n_parrays,
	   grIntPtr->dim.x0,  grIntPtr->dim.x1,
	    grIntPtr->dim.y0,  grIntPtr->dim.y1);
#ifdef DEBUG
    for (i = 0; i < grIntPtr->n_dlines; i++) {
	printf("d_array %d %f\n", grIntPtr->d_array[i].x0,
	       grIntPtr->d_array[i].x1,
	       grIntPtr->d_array[i].y0,
	       grIntPtr->d_array[i].y1);
    }

    for (i = 0; i < grIntPtr->n_pts; i++) {
	printf("p_array %d %f\n", grIntPtr->p_array[i].x,
	       grIntPtr->p_array[i].y);

    }
#endif
}

int IsGraphType(Tcl_Obj *objPtr)
{
    if (objPtr->typePtr == &graphType) {
	return 1;
    }
    return 0;

}

void PrintObj(Tcl_Obj *obj)
{
    /*
    printf("tcl_obj refcount %d bytes %s length %d ",
	   obj->refCount, obj->bytes, obj->length);
    */
    if (obj->typePtr)
	printf("name %s\n", obj->typePtr->name);
}


