#include <math.h>
#include <tk.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#include "canvas_box.h"
#include "misc.h"
#include "ruler_tick.h"
#include "container.h"

#define TICK_FACTOR         10/100
#define UFLOOR(x,u)         (floor((x)/(u))*(u))
#define UCEIL(x,u)          (ceil((x)/(u))*(u))

/* HACK - move to somewhere else */
#define ROUND(x)      ((x) < 0.0 ? ceil((x) - 0.5) : floor((x) + 0.5))

/*
 * ----------------------------------------------------------------------
 *
 * NiceNum --
 *
 * 	Taken from Paul Heckbert's "Nice Numbers for Graph Labels" in
 * 	Graphics Gems (pp 61-63).  Finds a "nice" number approximately
 * 	equal to x.
 *
 * ----------------------------------------------------------------------
 */
static double
NiceNum(double x, 
	int round)
{
    /* If "round" is non-zero, round. Otherwise take ceiling
				 * of value. */    
    double exponX;		/* exponent of x */
    double fractX;		/* fractional part of x */
    double nf;			/* nice, rounded fraction */

    exponX = floor(log10(x));

    fractX = (float)(x / (pow(10.0,(exponX)))); /* between 1 and 10 */
    if (round) {
	if (fractX < 1.5) {
	    nf = 1.0;
	} else if (fractX < 3.0) {
	    nf = 2.0;
	} else if (fractX < 7.0) {
	    nf = 5.0;
	} else {
	    nf = 10.0;
	}
    } else if (fractX <= 1.0) {
	nf = 1.0;
    } else if (fractX <= 2.0) {
	nf = 2.0;
    } else if (fractX <= 5.0) {
	nf = 5.0;
    } else {
	nf = 10.0;
    }
#ifdef DEBUG
    printf("x %.20f exp %.40f fract %.20f nf %f\n", x, exponX, fractX, nf);
    printf("return %.20f %.20f\n", pow(10.0,(exponX)), nf * pow(10.0,(exponX)));
#endif

    return (nf * pow(10.0,(exponX)));

}

void
ruler_ticks(double ruler_min, 
	    double ruler_max,
	    int def_num_ticks,
	    double *firstTick,
	    double *step,
	    int *numTicks)
{
    double range;
    double min, max;
    double majorStep;
    int numMajor;
    double tickMin, tickMax;

    min = ruler_min;
    max = ruler_max;

    range = max - min;

    if (range <= 0) {
	*firstTick = 0;
	*step = 0;
	*numTicks = 0;
	return;
    }

    /*
     * Calculate the major step.  If the user provided a step, use that value,
     * but only if it fits within the current range of the axis.
     */
    range = NiceNum(range, 0);
    majorStep = NiceNum(range / def_num_ticks, 1);

#ifdef DEBUG
    printf("range %.20f step %.20f num_ticks %d \n", range, majorStep, def_num_ticks);
#endif
    /*
     * Find the outer tick values in terms of the major step interval.  Add
     * +0.0 to preclude the possibility of an IEEE -0.0.
     */
/*
    tickMin = UFLOOR(min, majorStep) + 0.0;
    tickMax = UCEIL(max, majorStep) + 0.0; 
*/

    tickMin = UCEIL(min, majorStep) + 0.0;
    tickMax = UFLOOR(max, majorStep) + 0.0;
    range = tickMax - tickMin;
    numMajor = (int)ROUND(range / majorStep) + 1;

    *firstTick = tickMin;
    *step = majorStep;
    *numTicks = numMajor;
#ifdef DEBUG
    printf("firstTick %.20f step %.20f numTicks %d\n", *firstTick, *step, *numTicks);
#endif
}

static int FindNumTicks(CanvasPtr *canvas,
			int direction,
			double w1,
			double w2)
{
    double c1, c2, c;
    double w = 1.0;
    double range;
    int num_ticks;

    if (direction == HORIZONTAL) {
        WorldToCanvas(canvas, w1, w, &c1, &c);  
	WorldToCanvas(canvas, w2, w, &c2, &c);  
    } else {
        WorldToCanvas(canvas, w, w1, &c, &c1);  
	WorldToCanvas(canvas, w, w2, &c, &c2);  
    }

    range = c2 - c1;
    num_ticks = ROUND(range * TICK_FACTOR);
    return num_ticks;
}

static void
PlotTicks(Tcl_Interp *interp, 
	  ruler_s *ruler,
	  int offset,
	  int yoffset,
	  double firstTick,
	  double step,
	  int numTicks)
{
    int i;
    double x;
    char cmd[1024];
    x = firstTick;
    
#ifdef DEBUG
    printf("PlotTicks step %f ticks %d\n", step, numTicks);
#endif
    for (i = 0; i < numTicks; i++) {
	int height = (i % 5 == 4) ? ruler->tick.t.ht : ruler->tick.t.ht / 2;
	sprintf(cmd, "%s create line %.20f %d %.20f %d "
		"-fill %s -width %d -tag tick\n",
		ruler->window,
		x+offset, yoffset, x+offset, yoffset + height,
		ruler->tick.t.colour,
		ruler->tick.t.line_width);

#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	Tcl_Eval(interp, cmd);

	if (i % 5 == 4) {
	    sprintf(cmd, "%s create text %.20f %d -text %g -tag tick\n",
		    ruler->window, x+offset, yoffset+ruler->tick.offset, x);
#ifdef DEBUG
	    printf("%s\n", cmd);
#endif
	    Tcl_Eval(interp, cmd);
	}
	x += step;
    }
}

static void PlotTicks_v(Tcl_Interp *interp, 
			ruler_s *ruler,
			double wy1,
			double wy2,
			double firstTick,
			double step,
			int numTicks)
{
    int i;
    double y, tick_y;
    char cmd[1024];
    
    /* plot from the bottom upwards */
    tick_y = wy2 - firstTick + wy1;
    y = firstTick;

#ifdef DEBUG
    printf("PlotTicks step %d ticks %d wy1 %f wy2 %f\n", step, numTicks, wy1, wy2);
#endif
    for (i = 0; i < numTicks; i++) {
	int height = (i % 5 == 4) ? ruler->tick.t.ht : ruler->tick.t.ht / 2;

	sprintf(cmd, "%s create line %d %.20f %d %.20f "
		"-fill %s -width %d -tag tick\n",
		ruler->window,
		ruler->offset, tick_y, ruler->offset-height, tick_y,
		ruler->tick.t.colour,
		ruler->tick.t.line_width);

#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	Tcl_Eval(interp, cmd);

	if (i % 5 == 4) {
	    sprintf(cmd, "%s create text %d %.20f -text %g -tag tick\n",
		    ruler->window, ruler->tick.offset, tick_y, y);
#ifdef DEBUG
	    printf("%s\n", cmd);
#endif
	    Tcl_Eval(interp, cmd);
	}
	y += step;
	tick_y -= step;
    }
}

static void PlotTicks_c(Tcl_Interp *interp, 
			ruler_s *ruler,
			int wx1,
			int wx2,
			double origin,
			cir_s circle,
			double firstTick,
			double step,
			int numTicks)
{
    int i;
    double x, text;
    char cmd[1024];
    int seq_len = wx2 - wx1 + 1;
    double angle, theta;
    int circle_x, circle_y;
    double radius = circle.diameter/2;
    double x1, y1, x2, y2, x3, y3;
    int height;

    circle_x = circle.x1 + radius;
    circle_y = circle.y1 + radius;

    /* larger tick at the origin */
    theta = M_PI * origin / 180;
    height = 1.5 * ruler->tick.t.ht;
    x1 = circle_x - (radius * cos(theta));
    y1 = circle_y - (radius * sin(theta));
    x2 = circle_x - ((radius+height) * cos(theta));
    y2 = circle_y - ((radius+height) * sin(theta));
    x3 = circle_x - ((radius-height) * cos(theta));
    y3 = circle_y - ((radius-height) * sin(theta));

    sprintf(cmd, "%s create line %.20f %.20f %.20f %.20f "
	    "-fill %s -width %d -tag tick\n",
	    ruler->window, x1, y1, x2, y2, ruler->tick.t.colour,
	    ruler->tick.t.line_width);
    Tcl_Eval(interp, cmd);
    sprintf(cmd, "%s create text %.20f %.20f -text %.3g -tag tick\n",
	    ruler->window, x3, y3, (double)wx1);
    Tcl_Eval(interp, cmd);

    text = firstTick;
    x = firstTick-ruler->start;

#ifdef DEBUG
    printf("PlotTicks_c step %f ticks %d wx1 %d wx2 %d origin %f seq_len %d \n", step, numTicks, wx1, wx2, origin, seq_len);
#endif
    for (i = 0; i < numTicks; i++) {
	height = (i % 5 == 4) ? ruler->tick.t.ht : ruler->tick.t.ht / 2;

	/* angle = ((double)x / seq_len * 360) + (180- origin); */
	angle = ((double)x / seq_len * -360) + origin;
	theta = M_PI * angle / 180; /* in radians! */

	x1 = circle_x - (radius * cos(theta));
	y1 = circle_y - (radius * sin(theta));

	x2 = circle_x - ((radius+height) * cos(theta));
	y2 = circle_y - ((radius+height) * sin(theta));

	x3 = circle_x - ((radius-height) * cos(theta));
	y3 = circle_y - ((radius-height) * sin(theta));
	sprintf(cmd, "%s create line %.20f %.20f %.20f %.20f "
		"-fill %s -width %d -tag tick\n",
		ruler->window,
		x1, y1, x2, y2,
		ruler->tick.t.colour,
		ruler->tick.t.line_width);

#ifdef DEBUG
	printf("%s\n", cmd);
#endif
	Tcl_Eval(interp, cmd);

	if (i % 5 == 4) {
	    sprintf(cmd, "%s create text %.20f %.20f -text %.3g -tag tick\n",
		    ruler->window, x3, y3, text);
#ifdef DEBUG
	    printf("%s\n", cmd);
#endif
	    Tcl_Eval(interp, cmd);
	}
	x += step;
	text += step;
    }
}

void 
display_ruler_ticks(Tcl_Interp *interp, 
		    CanvasPtr *canvas, 
		    int x_offset,
		    int y_offset,
		    ruler_s *ruler,
		    int start,
		    int end)
{
    int def_num_ticks;
    double firstTick;
    double step;
    int numTicks;

    def_num_ticks = FindNumTicks(canvas, HORIZONTAL,
				 (double)start, (double)end);

    /* if no ticks then don't try to draw any! */
    if (def_num_ticks < 1) {
	return;
    }
    ruler_ticks(start, end, def_num_ticks, &firstTick, &step, 
		&numTicks);
   
    PlotTicks(interp, ruler, x_offset, y_offset,
	      (int)firstTick, MAX(1, (int)step), numTicks);
    
}

void display_ruler_ticks_v(Tcl_Interp *interp, 
			   CanvasPtr *canvas, 
			   ruler_s *ruler,
			   double wy1,
			   double wy2)
{
    int def_num_ticks;
    double firstTick;
    double step;
    int numTicks;

    def_num_ticks = FindNumTicks(canvas, VERTICAL, wy1, wy2);

    /* if no ticks then don't try to draw any! */
    if (def_num_ticks < 1) {
	return;
    }
    ruler_ticks(wy1, wy2, def_num_ticks, &firstTick, &step, 
		&numTicks);

    PlotTicks_v(interp, ruler, wy1, wy2, firstTick, step, numTicks);   
}

/* circular ruler ticks */
void display_ruler_ticks_c(Tcl_Interp *interp, 
			   CanvasPtr *canvas, 
			   ruler_s *ruler,
			   int wx1,
			   int wx2,
			   double origin,
			   cir_s circle)
{
    int def_num_ticks;
    double firstTick;
    double step;
    int numTicks;

    /* circumference of circle in canvas coords */
    def_num_ticks =  ROUND(circle.diameter * M_PI * TICK_FACTOR);

    /* if no ticks then don't try to draw any! */
    if (def_num_ticks < 1) {
	return;
    }
    ruler_ticks((double)wx1, (double)wx2, def_num_ticks, &firstTick, &step, 
		&numTicks);

    PlotTicks_c(interp, ruler, wx1, wx2, origin, circle, firstTick, step, numTicks);   
}

