#include <stdlib.h>
#include <string.h>
#include <tcl.h>

#include "canvas_box.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"
#include "misc.h"
#include "ruler_tick.h"
#include "element_canvas.h"

/*
 * debugging only
 */
void
printCanvas(CanvasPtr *c)
{

    printf("wd %d ht %d x %"PRId64" y %"PRId64" ax %.20f ay %.20f bx %.20f by %.20f\n",
	   c->width, c->height, c->x, c->y, c->ax, c->ay, c->bx, c->by);

}

void printWorld(WorldPtr *world)
{
    printf("Total x1 %f y1 %f x2 %f y2 %f\n"
	   "Visible x1 %f y1 %f x2 %f y2 %f\n",
	   world->total->x1, world->total->y1, 
	   world->total->x2, world->total->y2,
	   world->visible->x1, world->visible->y1, 
	   world->visible->x2, world->visible->y2);	   

}

void
initCanvas(Tcl_Interp *interp,
	   CanvasPtr *canvas,        /* in, out */
	   char *window)
{
    Tcl_VarEval(interp, "winfo width ", window, NULL);
    canvas->width = atoi(Tcl_GetStringResult(interp)) - 1;

    Tcl_VarEval(interp, "winfo height ", window, NULL);
    canvas->height = atoi(Tcl_GetStringResult(interp)) - 1;
    
    canvas->x = 0;
    canvas->y = 0;
    canvas->ax = 0;
    canvas->ay = 0;
    canvas->bx = 0;
    canvas->by = 0;
#ifdef DEBUG
    printf("--\nINIT\n");
    printCanvas(canvas);
#endif
}

/*
 * convert between canvas (pixel) coords and world (base) coords
 * canvas = world * m + c 
 * world = (canvas - c) / m
 */
void
CanvasToWorld(CanvasPtr *canvas,
	      int64_t cx, int64_t cy,
	      double *wx, double *wy)
{
   *wx = (cx - canvas->bx) / canvas->ax;
   *wy = (cy - canvas->by) / canvas->ay; 

#ifdef DEBUG
   printf("CanvasToWorld wx %f wy %f cx %d cy %d\n", *wx, *wy, cx, cy); 
#endif

}

/*
 * convert between world (base) coords and canvas (pixel) coords
 * canvas = world * m + c 
 */
void
WorldToCanvas(CanvasPtr *canvas,
	      double wx, double wy, 
	      double *cx, double *cy)
{
    /*
    *cx = (int)ROUND((wx * canvas->ax) + canvas->bx);
    *cy = (int)ROUND((wy * canvas->ay) + canvas->by);
    */
    *cx = (wx * canvas->ax) + canvas->bx;
    *cy = (wy * canvas->ay) + canvas->by;
#ifdef DEBUG

    printf("WorldToCanvas %f %f ay %f by %f\n", wy, *cy, canvas->ay, canvas->by);

    printf("WorldToCanvas wx %f cx %f ax %f bx %f %f\n", 
	   wx, *cx, canvas->ax, canvas->bx, (wx * canvas->ax) + canvas->bx);
#endif
}

double
calc_zoom_origin(double min1, double min2, double max1, double max2)
{
#ifdef DEBUG
    printf("or min1 %f min2 %f max1 %f max2 %f\n", min1, min2, max1, max2);
#endif
    if ((max2 + min1 - max1 - min2) == 0.0) {
	return 0.0;
    } else {
	return (((min1 * max2) - (max1 * min2)) / (max2 + min1 - max1 - min2));
    }
}

double
calc_zoom_sf(double min1, double min2, double max1, double max2)
{
#ifdef DEBUG
    printf("sf min1 %f min2 %f max1 %f max2 %f\n", min1, min2, max1, max2);
#endif
    if (max2 - min2 == 0.0) 
	return 1.0;
    else
	return ((max1 - min1) / (max2 - min2));
}

void
scrollRegion(Tcl_Interp *interp,
	     win **win_list,
	     int num_wins,
	     d_box *total,
	     CanvasPtr *canvas)
{
    char cmd[1024];
    int i;
    double x1, y1, x2, y2;
    x1 = x2 = y1 = y2 = 0;

    for (i = 0; i < num_wins; i++) {

	WorldToCanvas(canvas, total->x1, total->y1, &x1, &y1);
	WorldToCanvas(canvas, total->x2, total->y2, &x2, &y2);

	/* scrollregion on right hand side is not inclusive so must add 1 */
	x2++;

	if (win_list[i]->scroll == 'x') {
	    y1 = 0;
	    y2 = 0;
	}
	if (win_list[i]->scroll == 'y') {
	    x1 = 0;
	    x2 = 0;
	}
	if (win_list[i]->scroll == 'n') {
	    x1 = 0;
	    x2 = 0;
	    y1 = 0;
	    y2 = 0;
	}
	sprintf(cmd, "%s configure -scrollregion \"%.20f %.20f %.20f %.20f\"", 
		win_list[i]->window, x1, y1, x2, y2);
#ifdef DEBUG
	printf("--\n");
	printf("scroll region %s\n", cmd);
#endif
	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    verror(ERR_WARN, "scrollRegion", "%s\n",
		   Tcl_GetStringResult(interp));
    }
}

/*
 * scale all canvases by the current scale factor
 */
void
scaleCanvas(Tcl_Interp *interp,
	    win **win_list,
	    int num_wins,
	    char *tags,
	    d_box *current, /* world coords */
	    CanvasPtr *canvas) 
{
    char cmd[1024];
    double x_origin, y_origin;
    double sf_x, sf_y;
    int64_t c_x1, c_y1, c_x2, c_y2;
    int i;

    c_x1 = canvas->x;
    c_y1 = canvas->y;
    c_x2 = canvas->x + canvas->width;
    c_y2 = canvas->y + canvas->height;

    x_origin = calc_zoom_origin(current->x1, (double)c_x1, 
				current->x2, (double)c_x2);

    y_origin = calc_zoom_origin(current->y1, (double)c_y1, 
				current->y2, (double)c_y2);

    sf_x = calc_zoom_sf((double)c_x1, current->x1, 
			(double)c_x2, current->x2);
    sf_y = calc_zoom_sf((double)c_y1, current->y1, 
			(double)c_y2, current->y2);

    /*
     * check if the window the scale was executed in can only scroll in 
     * x or y or neither
     */

    for (i = 0; i < num_wins; i++) {

	if (win_list[i]->scroll == 'x') {
	    if (current->x1 == c_x1 && current->x2 == c_x2) {
		/* do nothing */
		sprintf(cmd, "%s scale %s %.20f %.20f %.20f %.20f", 
			win_list[i]->window, tags, 0.0, 0.0, 1.0, 1.0);

	    } else if (x_origin == 0.0 && sf_x == 1.0) {
		sprintf(cmd, "%s move %s %"PRId64" %d",
			win_list[i]->window, tags, canvas->x, 0);
	    } else {
		sprintf(cmd, "%s scale %s %.20f %.20f %.20f %.20f", 
			win_list[i]->window, tags, x_origin, 0.0, sf_x, 1.0);
	    }
	} else if (win_list[i]->scroll == 'y') {
	    if (current->y1 == c_y1 && current->y2 == c_y2) {
		/* do nothing */
		sprintf(cmd, "%s scale %s %.20f %.20f %.20f %.20f", 
			win_list[i]->window, tags, 0.0, 0.0, 1.0, 1.0);
	    } else if (y_origin == 0.0 && sf_y == 1.0) {
		sprintf(cmd, "%s move %s %d %"PRId64,
			win_list[i]->window, tags, 0, canvas->y);
	    } else {
		sprintf(cmd, "%s scale %s %.20f %.20f %.20f %.20f", 
			win_list[i]->window, tags, 0.0, y_origin, 1.0, sf_y);
	    }
	} else if (win_list[i]->scroll == 'n') {
	    sprintf(cmd, "%s scale %s %.20f %.20f %.20f %.20f", 
		    win_list[i]->window, tags, 0.0, 0.0, 1.0, 1.0);
	} else {
	    if (current->x1 == c_x1 && current->x2 == c_x2 && current->y1 == c_y1 && current->y2 == c_y2) {
		sprintf(cmd, "%s scale %s %.20f %.20f %.20f %.20f", 
			win_list[i]->window, tags, 0.0, 0.0, 1.0, 1.0);
	    } else if (x_origin == 0.0 && sf_x == 1.0 && y_origin == 0.0 && sf_y == 1.0){ 
		sprintf(cmd, "%s move %s %"PRId64" %"PRId64,
			win_list[i]->window, tags, canvas->x, canvas->y);
		if (TCL_ERROR == Tcl_Eval(interp, cmd))
		    verror(ERR_WARN, "moveCanvas", "%s\n",
			   Tcl_GetStringResult(interp));
	    } else {
		sprintf(cmd, "%s scale %s %.20f %.20f %.20f %.20f", 
			win_list[i]->window, tags, x_origin, y_origin, sf_x, sf_y);
	    }
	}
#ifdef DEBUG
	printf("--\nSCALE %s\n", cmd);
	printCanvas(canvas);
#endif
	if (TCL_ERROR == Tcl_Eval(interp, cmd))
	    verror(ERR_WARN, "scaleCanvas", "%s\n",
		   Tcl_GetStringResult(interp));
    }
}

void
scaleSingleCanvas(Tcl_Interp *interp,
		  WorldPtr *world,
		  CanvasPtr *canvas,
		  char *window,
		  char scroll,
		  char *tag) 
{
     win *win_list;

    if (NULL == (win_list = (win *)xmalloc(sizeof(win))))
        return;

     win_list->window = strdup(window);
     win_list->scroll = scroll;

     scaleCanvas(interp, &win_list, 1, tag, world->visible, canvas);
     scrollRegion(interp, &win_list, 1, world->total, canvas);

     free(win_list->window);
     xfree(win_list);
}

void
canvasScrollX(Tcl_Interp *interp,
	      char *window,
	      win **win_list,
	      int num_wins,
	      d_box *visible,
	      CanvasPtr *canvas,
	      char *scroll_args)
{
    char cmd[1024];
    int i;
    double wy;

    if (num_wins == 0) 
	return;

    for (i = 0; i < num_wins; i++) {

	if ((win_list[i]->scroll == 'x') || (win_list[i]->scroll == 'b')) {
	    sprintf(cmd, "eval %s xview %s", win_list[i]->window, scroll_args);
	    Tcl_Eval(interp, cmd);
	}
    }

    /* find new left hand edge of canvas in canvasx coords */
    canvas->x = canvas_x(interp, window, 0.0);

    /* find new left and right edges of canvas in world coords */
    CanvasToWorld(canvas, canvas->x, 0, &visible->x1, &wy);
    CanvasToWorld(canvas, canvas->x + canvas->width, 0, &visible->x2, &wy);
    SetCanvasCoords(interp, visible->x1, visible->y1, 
		    visible->x2, visible->y2, canvas);
}

void
canvasScrollY(Tcl_Interp *interp,
	      char *window,
	      win **win_list,
	      int num_wins,
	      d_box *visible,
	      CanvasPtr *canvas,
	      char *scroll_args)
{
    char cmd[1024];
    int i;
    double wx;

    if (num_wins == 0) 
	return;

    for (i = 0; i < num_wins; i++) {

	if ((win_list[i]->scroll == 'y') || (win_list[i]->scroll == 'b')) {
	    sprintf(cmd, "eval %s yview %s", win_list[i]->window, scroll_args);
	    Tcl_Eval(interp, cmd);
	}
    }

    /* find new top edge of canvas in canvasy coords */
    Tcl_VarEval(interp, window, " canvasy 0", NULL);
    canvas->y = atoi(Tcl_GetStringResult(interp));

    /* find new top and bottom edges of canvas in world coords */
    CanvasToWorld(canvas, 0, canvas->y, &wx, &visible->y1);
    CanvasToWorld(canvas, 0, canvas->y + canvas->height, &wx, &visible->y2);
    SetCanvasCoords(interp, visible->x1, visible->y1, 
		    visible->x2, visible->y2, canvas);
}

void
resizeCanvas(Tcl_Interp *interp,
	     char *window,
	     win **win_list,
	     int num_wins,
	     d_box *visible,
	     d_box *total,
	     CanvasPtr *canvas)
{
    d_box *bbox;
    int height, width;

    if (NULL == (bbox = (d_box *)xmalloc(sizeof(d_box))))
	return;

    bbox->x1 = (double)canvas->x;
    bbox->y1 = (double)canvas->y;
    bbox->x2 = (double)(canvas->x + canvas->width);
    bbox->y2 = (double)(canvas->y + canvas->height);

    Tcl_VarEval(interp, "winfo width ", window, NULL);
    width = atoi(Tcl_GetStringResult(interp)) - 1;

    Tcl_VarEval(interp, "winfo height ", window, NULL);
    height = atoi(Tcl_GetStringResult(interp)) - 1;

#ifdef DEBUG
    printf("resizeCanvas \n");
    printf("window %s width %d %d height %d %d\n", window, width, canvas->width, height, canvas->height);
#endif
    /* if changed height or width */
    if ((height != canvas->height) || (width != canvas->width)) {
	canvas->height = height;
	canvas->width = width;
	SetCanvasCoords(interp, visible->x1, visible->y1, 
			visible->x2, visible->y2, canvas);

	scaleCanvas(interp, win_list, num_wins, "all", bbox, canvas);
	scrollRegion(interp, win_list, num_wins, total, canvas);
    }
    xfree(bbox);
}


void
listZoom(StackPtr *zoom)
{
    int i = 0;
    StackPtr *stack;

    stack = zoom;
    while (stack != NULL) {
	printf("list %d x1 %f y1 %f x2 %f y2 %f\n", i, stack->data->x1, 
	       stack->data->y1, stack->data->x2, stack->data->y2);
	stack = stack->next;
	i++;
    }

}

int
lengthZoom(StackPtr *zoom)
{
    int i = 0;
    StackPtr *stack;

    stack = zoom;
    while (stack != NULL) {
	stack = stack->next;
	i++;
    }
    return i;
}

void
createZoom(StackPtr **zoom)
{
    *zoom = NULL;
#ifdef DEBUG
    listZoom(*zoom);
#endif
}

d_box * examineZoom(StackPtr *zoom)
{

#ifdef DEBUG
    listZoom(zoom);
#endif
    if (zoom == NULL) {
	return (NULL);
    } else {    
	return(zoom->data);
    }
}

/*
 * get previous zoom values
 */
void popZoom(StackPtr **zoom)
{
    StackPtr *stack;

    /* don't remove the last item from the stack */
    if (*zoom == NULL || (*zoom)->next == NULL) {
#ifdef DEBUG
	printf("end of stack \n");
#endif
    } else {
	stack = *zoom;
	*zoom = (*zoom)->next;
	if (stack->data) xfree(stack->data);
	xfree(stack);
    }
}

/*
 * save current zoom values
 */
void pushZoom(StackPtr **zoom, d_box *visible)
{
    StackPtr *stack;

#ifdef DEBUG
    listZoom(*zoom);
#endif
    stack = (StackPtr *)xmalloc(sizeof(StackPtr));
    stack->data = (d_box *)xmalloc(sizeof(d_box));
    
    memcpy(stack->data, visible, sizeof(d_box));
    stack->next = *zoom;
    *zoom = stack;

#ifdef DEBUG
    listZoom(*zoom);
#endif
}

/* copy zoom1 to zoom2 */
void copyZoom(StackPtr **zoom2,
	      StackPtr *zoom1)
{
    StackPtr *stack1p;
    StackPtr *stack2, *stack2p;
    StackPtr *tmp;

    createZoom(zoom2);

    stack1p = zoom1;
    stack2 = *zoom2;
    stack2p = stack2;

    /* copy stack1 to stack2 */
    while (stack1p != NULL) {
        
        /* create a new item for stack2 */
        tmp = (StackPtr *)xmalloc(sizeof(StackPtr));
	tmp->data = (d_box *)xmalloc(sizeof(d_box));
	memcpy(tmp->data, stack1p->data, sizeof(d_box));

	/* first item in stack2 */
	if (stack2 == NULL) {
	    stack2 = tmp;
	    stack2p = tmp;
	} else {
	    /* subsequence items in stack2 */
	    stack2p->next = tmp;
	    stack2p = tmp;
	}
	stack1p = stack1p->next;
    }
    stack2p->next = *zoom2;
    *zoom2 = stack2;
}

void freeZoom(StackPtr **zoom)
{
    StackPtr *stack;

    while (*zoom != NULL) {
	stack = *zoom;
	*zoom = (*zoom)->next;
	if (stack->data) xfree(stack->data);
	xfree(stack);
    }
}

/*
 * canvas = world * m + c 
 */
void
SetCanvasCoords(Tcl_Interp *interp,
		double x1, double y1, double x2, double y2,
		CanvasPtr *canvas)     /* out */
{
    if (x2 - x1 == 0) {
	canvas->ax = 1;
    } else {
	canvas->ax = canvas->width / (x2-x1);
    }
    if (y2 - y1 == 0) {
	canvas->ay = 1;
    } else {
	canvas->ay = canvas->height / (y2-y1);
    }
    canvas->bx = canvas->x - (canvas->ax * x1);
    canvas->by = canvas->y - (canvas->ay * y1);

#ifdef DEBUG
    printf("********************SETCANVAS COORDS \n");
    printf("x1 %f x2 %f y1 %f y2 %f \n", x1, x2, y1, y2);
    printf("x %d y %d ax %f ay %f bx %f by %f width %d\n", canvas->x, canvas->y, canvas->ax, 
	   canvas->ay, canvas->bx, canvas->by, canvas->width);
#endif
}

void
canvasCursorX(Tcl_Interp *interp,
	      CanvasPtr *canvas,
	      char *frame,
	      char *label,
	      char *colour,
	      int line_width,
	      int64_t cx,
	      double wx,
	      win **win_list,
	      int num_wins)
{
    char cmd[1024];
    int i;

    sprintf(cmd, "%s%s configure -text %"PRId64"\n",
	    frame, label, (int64_t)wx);
    Tcl_Eval(interp, cmd);
    
    /* draw the cursor in each window */
    for (i = 0; i < num_wins; i++) {
	/* only draw cursors in x direction */
	if ((win_list[i]->scroll == 'x') || (win_list[i]->scroll == 'b')) {
	    sprintf(cmd, "DrawCanvasCursorX %s %s %"PRId64" %s %d\n", 
		    frame,win_list[i]->window, cx, colour, line_width);
	    if (TCL_ERROR == Tcl_Eval(interp, cmd))
		verror(ERR_WARN, "canvasCursorX", "%s\n",
		       Tcl_GetStringResult(interp));
	}
    }
}

void
canvasCursorY(Tcl_Interp *interp,
	      CanvasPtr *canvas,
	      char *frame,
	      char *label,
	      char *colour,
	      int line_width,
	      int64_t cy,
	      double wy,
	      win **win_list,
	      int num_wins)
{
    char cmd[1024];
    int i;

    sprintf(cmd, "%s%s configure -text %"PRId64"\n",
	    frame, label, (int64_t)wy);
    Tcl_Eval(interp, cmd);
    
    /* draw the cursor in each window */
    for (i = 0; i < num_wins; i++) {
	/* only draw cursors in y direction */
	if ((win_list[i]->scroll == 'y') || (win_list[i]->scroll == 'b')) {
	    sprintf(cmd, "DrawCanvasCursorY %s %s %"PRId64" %s %d\n", 
		    frame,win_list[i]->window, cy, colour, line_width);
	    if (TCL_ERROR == Tcl_Eval(interp, cmd))
		verror(ERR_WARN, "canvasCursorY", "%s\n",
		       Tcl_GetStringResult(interp));
	}
    }
}

void
canvasZoom(Tcl_Interp *interp,
	   CanvasPtr *canvas,
	   char *window,
	   WorldPtr *world,
	   win **win_list,
	   int num_wins,
	   StackPtr **zoom_list,
	   box *zoom,
	   char scroll)
{	
    d_box *bbox;
    double x1, y1, x2, y2;

    if (num_wins <= 0)
	return;

    x1 = world->visible->x1;
    y1 = world->visible->y1;
    x2 = world->visible->x2;
    y2 = world->visible->y2;

    CanvasToWorld(canvas, zoom->x1, zoom->y1, 
		  &world->visible->x1, &world->visible->y1);
    CanvasToWorld(canvas, zoom->x2, zoom->y2, 
		  &world->visible->x2, &world->visible->y2);
    
    if (NULL == (bbox = (d_box *)xmalloc(sizeof(d_box))))
	return;

    bbox->x1 = (double)zoom->x1;
    bbox->y1 = (double)zoom->y1;
    bbox->x2 = (double)zoom->x2;
    bbox->y2 = (double)zoom->y2;

    /*
     * need to take into account how the originating zooming window scrolled
     * so that eg a window that only scrolls in x only zooms in x
     */  
    if (scroll == 'x' || scroll == 'n') {
	world->visible->y1 = y1;
	world->visible->y2 = y2;
	bbox->y1 = 0.0;
	bbox->y2 = 0.0;
    }
    if (scroll == 'y' || scroll == 'n') {
	world->visible->x1 = x1;
	world->visible->x2 = x2;
	bbox->x1 = 0.0;
	bbox->x2 = 0.0;
    }
    SetCanvasCoords(interp, world->visible->x1, world->visible->y1, 
		    world->visible->x2, world->visible->y2, canvas);
    
    scaleCanvas(interp, win_list, num_wins, "all", bbox, canvas);
    scrollRegion(interp, win_list, num_wins, world->total, canvas);
    
    pushZoom(zoom_list, world->visible);
    
    canvas->x = canvas_x(interp, window, 0.0);

    xfree(bbox);
}

void
canvasZoomback(Tcl_Interp *interp,
	       CanvasPtr *canvas,
	       char *window,
	       WorldPtr *world,
	       win **win_list,
	       int num_wins,
	       StackPtr **zoom_list)
{
    d_box *zoom;
    
    if (num_wins <= 0)
	return;

    if (NULL == (zoom = (d_box *)xmalloc(sizeof(d_box))))
	return;
	
    popZoom(zoom_list);
    if (examineZoom(*zoom_list) == NULL) {
	return;
    }

    memcpy(world->visible, examineZoom(*zoom_list), sizeof(d_box));
    WorldToCanvas(canvas, world->visible->x1, world->visible->y1, 
		  &zoom->x1, &zoom->y1);
    WorldToCanvas(canvas, world->visible->x2, world->visible->y2, 
		  &zoom->x2, &zoom->y2);
    
    scaleCanvas(interp, win_list, num_wins, "all", zoom, canvas);
		
    /* set new canvas coords from zoom */
    SetCanvasCoords(interp, world->visible->x1, world->visible->y1, 
		    world->visible->x2, world->visible->y2, canvas);

    scrollRegion(interp, win_list, num_wins, world->total, canvas);
 
    canvas->x = canvas_x(interp, window, 0.0);

    xfree(zoom);
}

tag_s tag_struct(Tcl_Interp *interp,
		 Tcl_Obj *defs_ptr,
		 char *win_name,
		 int width,
		 int offset)
{
    tag_s tag_config;

    if (width == -1) {
	tag_config.width =
	    get_default_int(interp, defs_ptr, 
			    vw("%s.TAG_WIDTH", win_name));
    } else {
	tag_config.width = width;
    }

    if (offset == -1) {
	tag_config.offset =
	    get_default_int(interp, defs_ptr, 
			    vw("%s.TAG_OFFSET", win_name));
    } else {
	tag_config.offset = offset;
    }

    return tag_config;
}

cursor_s cursor_struct(Tcl_Interp *interp,
		       Tcl_Obj *defs_ptr,
		       char *win_name,
		       int width,
		       char *colour)
{
    cursor_s cursor_config;

    if (width == -1) {
	cursor_config.width =
	    get_default_int(interp, defs_ptr, 
			    vw("%s.CURSOR_WIDTH", win_name));
    } else {
	cursor_config.width = width;
    }

    if (strcmp(colour, "") == 0) {
	cursor_config.colour =
	    get_default_astring(interp, defs_ptr, 
				vw("%s.CURSOR_COLOUR", win_name));
    } else {
	cursor_config.colour = strdup(colour);
    }
 
    return cursor_config;
}


tick_s *tick_struct(Tcl_Interp *interp,
		    Tcl_Obj *defs_ptr,
		    char *win_name,
		    int width,
		    int height,
		    char *colour)
{
    tick_s *tick_config;

   if (NULL == (tick_config = (tick_s *)xmalloc(sizeof(tick_s))))
        return NULL;

    if (width == -1) {
	tick_config->line_width =
	    get_default_int(interp, defs_ptr, 
			    vw("%s.TICK_WIDTH", win_name));
    } else {
	tick_config->line_width = width;
    }

    if (height == -1) {
	tick_config->ht =
	    get_default_int(interp, defs_ptr,
			    vw("%s.TICK_HEIGHT", win_name));
    } else {
	tick_config->ht = height;
    }

    if (strcmp(colour, "") == 0) {
	tick_config->colour =
	    get_default_astring(interp, defs_ptr,
				vw("%s.TICK_COLOUR", win_name));
    } else {
	tick_config->colour = strdup(colour);
    }
    return tick_config;
}

ruler_s *ruler_struct(Tcl_Interp *interp,
		     Tcl_Obj *defs_ptr,
		     char *win_name,
		     int tags)
{
    ruler_s *ruler;
 
    if (NULL == (ruler = (ruler_s *)xmalloc(sizeof(ruler_s))))
        return NULL;
 
    if (NULL == (ruler->window = (char *)xmalloc(100 * sizeof(char))))
      return NULL;

    ruler->tick.t.ht =
	get_default_int(interp, defs_ptr,
			vw("%s.RULER.TICK_HEIGHT", win_name));

    ruler->tick.t.line_width =
	get_default_int(interp, defs_ptr,
			vw("%s.RULER.TICK_WIDTH", win_name));

    ruler->tick.t.colour =
	get_default_astring(interp, defs_ptr,
			    vw("%s.RULER.TICK_COLOUR", win_name));

    ruler->tick.offset =
	get_default_int(interp, defs_ptr,
			vw("%s.RULER.TEXT_OFFSET", win_name));

    ruler->offset =
	get_default_int(interp, defs_ptr,
			vw("%s.RULER.OFFSET", win_name));

    ruler->colour =
	get_default_astring(interp, defs_ptr,
			    vw("%s.RULER.COLOUR", win_name));

    ruler->line_width =
	get_default_int(interp, defs_ptr,
			vw("%s.RULER.LINE_WIDTH", win_name));

    /* not alway have these values set */
    if (tags) {
	ruler->tag.offset =
	    get_default_int(interp, defs_ptr,
			    vw("%s.RULER.TAG_OFFSET", win_name));
	
	ruler->tag.width =
	    get_default_int(interp, defs_ptr,
			    vw("%s.RULER.TAG_WIDTH", win_name));
	ruler->tick.num = 0;
    }

    return ruler;
}

void free_ruler_struct(ruler_s *r) {
    if (!r)
	return;

    if (r->window) free(r->window);
    if (r->colour) free(r->colour);
    if (r->tick.t.colour) free(r->tick.t.colour);
    xfree(r);
}

int addWindow(win **win_list,
	      int *num_wins,
	      char *window,
	      char scroll,
	      int id)
{
    int i;

    if (*num_wins == MAX_NUM_WINS) {
	verror(ERR_WARN, "addWindow", "too many windows \n");
	return -1;
    }
    /* check the window doesn't already exist in the win_list */
    for (i = 0; i < *num_wins; i++) {
	if (strcmp(win_list[i]->window, window) == 0) {
	    return 0;
	}
    }

    if (NULL == (win_list[*num_wins] = (win *)xmalloc(sizeof(win))))
        return -1;
    
    win_list[*num_wins]->window = strdup(window);
    win_list[*num_wins]->scroll = scroll;
    win_list[*num_wins]->id = id;
    (*num_wins)++;
    return 0;
}

void deleteWindow(win **win_list,
		  int *num_wins,
		  char *window)
{
    int i, length;

    for (i = 0; i < *num_wins; i++) {
	if (strcmp(win_list[i]->window, window) == 0) {
	    xfree(win_list[i]->window);
	    xfree(win_list[i]);
	    length = *num_wins - i - 1;
	    if (length > 0) 
	      memmove(&(win_list)[i], &(win_list)[i+1], length * sizeof(win*));
	    (*num_wins)--;
	}
    }
}

void free_win_list(win **win_list, int num_wins) {
    int i;

    if (!win_list)
	return;

    for (i = 0; i < num_wins; i++)
      if (win_list[i]->window) {
	  free(win_list[i]->window);
	  xfree(win_list[i]);
      }
    free(win_list);
}

void draw_single_ruler(Tcl_Interp *interp,
		       ruler_s *ruler,
		       CanvasPtr *canvas,
		       double start,
		       double end,
		       int disp_ticks)
{
    char cmd[1024];

    /* remove current ruler before drawing the new ruler */
    Tcl_VarEval(interp, ruler->window, " delete all", NULL);

    sprintf(cmd, "%s create line %.20f %d %.20f %d -fill %s -width %d",
	    ruler->window, start, ruler->offset, end,
	    ruler->offset, ruler->colour, ruler->line_width);
    
    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "draw_single_ruler", "%s\n",
	       Tcl_GetStringResult(interp));
    
    if (disp_ticks)
	display_ruler_ticks(interp, canvas, 0, ruler->offset, ruler, 
			    start, end);
}

void draw_single_ruler_vertical(Tcl_Interp *interp,
				ruler_s *ruler,
				CanvasPtr *canvas,
				double start,
				double end,
				int disp_ticks)
{
    char cmd[1024];

    /* remove current ruler before drawing the new ruler */
    Tcl_VarEval(interp, ruler->window, " delete all", NULL);

    sprintf(cmd, "%s create line %d %.20f %d %.20f -fill %s -width %d",
	    ruler->window, ruler->offset, start, ruler->offset,
	    end, ruler->colour, ruler->line_width);

#ifdef DEBUG
    printf("draw_single_ruler_vertical %s\n", cmd);
#endif

    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "draw_single_ruler_vertical", "%s\n",
	       Tcl_GetStringResult(interp));

    if (disp_ticks)
	display_ruler_ticks_v(interp, canvas, ruler, start, end);
}

int canvasx(Tcl_Interp *interp,
	    char *window,
	    double val)
{
  char cmd[1024];

  sprintf(cmd, "%s canvasx %.20f", window, val);
  Tcl_Eval(interp, cmd);
  return (atoi(Tcl_GetStringResult(interp)));
}

int canvasy(Tcl_Interp *interp,
	    char *window,
	    double val)
{
  char cmd[1024];

  sprintf(cmd, "%s canvasy %.20f", window, val);
  Tcl_Eval(interp, cmd);
  return (atoi(Tcl_GetStringResult(interp)));
}
