#include <string.h>
#include <stdio.h>

#include <tk.h>
#include <X11/Xlib.h>

#ifdef _WIN32
#include "tkWinport.h"
#endif

#include <io_lib/Read.h>

#include "tkTrace.h"
#include "tkTraceIO.h"
#include "xalloc.h"
#include "misc.h"


int pixel_to_base(DNATrace *tracePtr, int pixel, int allow_end) {
    int cur, nearest, dist, t;
    int pos = pixel_to_point(tracePtr, pixel - tracePtr->borderWidth - 1);

    if (pos < 0)
	pos = 0;

    if (pos >= tracePtr->read->NPoints)
	pos = tracePtr->read->NPoints -1;

    nearest = tracePtr->tracePosE[pos];

    if (!allow_end) {
	while (nearest < tracePtr->Ned- 1) {
	    nearest++;
	    if (tracePtr->edPos[nearest])
		break;
	}
    } else {
	while (nearest < tracePtr->Ned) {
	    nearest++;
	    if (tracePtr->edPos[nearest])
		break;
	}
    }

    dist = trace_get_pos(tracePtr, nearest) - pos;
    if (dist < 0)
	dist = 9999;

    cur = nearest;
    do {
	cur--;

	if (cur >= 0) {
	    t = trace_get_pos(tracePtr, cur) - pos;

	    if (t <= 0)
		t = 9999;

	    if (t < dist) {
		dist = t;
		nearest = cur;
	    }
	}
    } while (cur >= 0 && t != 9999);

    return nearest;
}

void trace_update_extents(DNATrace *t, int *x0p, int *xnp) {
    int ind, pos, x1, fw, cfw, new_x1;
    int min =  999999;
    int max = -999999;
    int x0 = *x0p, xn = *xnp;
    int new_x0, new_xn;
    
    if (t->Ned <= 0)
	return;

    if (x0 < 0) x0 = 0;
    if (x0 >= t->read->NPoints) x0 = t->read->NPoints-1;
    x1 = x0 + xn < t->read->NPoints ? x0+xn : t->read->NPoints - 1;
    x1 = t->tracePos[x1] + 1;
    if (x1 >= t->read->NBases) x1--;
    x1 = t->read->basePos[x1];
    ind = t->tracePosE[x0];       /* Index of first base on screen */

    fw = t->font_width/2 + 1;
    cfw = t->conf_font_width;

    while (ind < t->read->NBases && (pos = trace_get_pos(t, ind)) <= x1) {
	pos = point_to_pixel(t, pos) - fw;
	if (pos < min)
	    min = pos;
	if (pos + cfw > max)
	    max = pos + cfw;
	ind++;
    }

    new_x0 = (int)pixel_to_point(t, min-cfw/2-1);
    new_xn = (int)pixel_to_point(t, max+cfw/2+1) - new_x0;

    x1 = x0 + xn;
    new_x1 = new_x0 + new_xn;
    x0 = MIN(x0, new_x0);
    x1 = MAX(x1, new_x1);
    xn = x1 - x0;

    if (x0 < 0) {
	xn += x0;
	x0 = 0;
    }
    if (x0 + xn > t->read->NPoints)
	xn = t->read->NPoints - x0;

    *x0p = x0;
    *xnp = xn;
}

void trace_draw_confidence(DNATrace *t, Display *d, Pixmap p,
			   int x0, int xn, int height) {
    int ind, pos, x1, fw, cfw, cfh;
    char b[5];

    if (!p || t->Ned <= 0)
	return;

    x1 = x0 + xn < t->read->NPoints ? x0+xn : t->read->NPoints - 1;
    x1 = t->tracePos[x1] + 1;
    if (x1 >= t->read->NBases) x1--;
    x1 = t->read->basePos[x1];
    ind = t->tracePosE[x0];       /* Index of first base on screen */

    fw = t->font_width/2 + 1;
    cfh = t->conf_font_height;
    cfw = t->conf_font_width;

    while (ind < t->read->NBases && (pos = trace_get_pos(t, ind)) <= x1) {
	int conf;

	conf = t->edConf[ind];
	if (conf < 100)
	    sprintf(b, "%02d", conf);
	else
	    strcpy(b, "XX");

	pos = point_to_pixel(t, pos) - fw;
	XFillRectangle(d, p, t->ConfGC, pos, cfh,
		       cfw, (int)(conf * t->scale_conf));

	Tk_DrawChars(d, p, t->ConfGC, t->conf_font, b, 2,
		     pos, cfh);
	ind++;
    }
}

/*
 * Draws the 4 confidence values produced by Solexa in their +/- log-odds
 * encoding.
 *
 * NB: Only works when trace editing hasn't been performed
 */
void trace_draw_confidence4(DNATrace *t, Display *d, Pixmap p,
			    int x0, int xn, int height) {
    int ind, x1, fw, cfh, cfw;
    char b[5];
    double pos;

    if (!p || t->Ned <= 0)
	return;

    x1 = x0 + xn < t->read->NPoints ? x0+xn : t->read->NPoints - 1;
    x1 = t->tracePos[x1] + 1;
    if (x1 >= t->read->NBases) x1--;
    x1 = t->read->basePos[x1];
    ind = t->tracePosE[x0];       /* Index of first base on screen */

    fw = t->font_width/2 + 1;
    cfh = t->conf_font_height;
    cfw = t->conf_font_width;

    while (ind < t->read->NBases && (pos = trace_get_pos(t, ind)) <= x1) {
	int conf;
	int edind;
	double xoff;

	/* Called confidence value */
	conf = (signed char)(t->edConf[ind]);
	if (conf < 100)
	    sprintf(b, "%02d", conf);
	else 
	    strcpy(b, "XX");

	xoff = 0;
	switch (t->read->base[ind]) {
	case 'A': case 'a': xoff = 0.00; break;
	case 'C': case 'c': xoff = 0.15; break;
	case 'G': case 'g': xoff = 0.30; break;
	case 'T': case 't': xoff = 0.45; break;
	}

	Tk_DrawChars(d, p, t->ConfGC, t->conf_font, b, 2,
		     point_to_pixel(t, pos + xoff)-fw, cfh);

	/* Graphics bars for 4 confidence values */
	edind = t->edPos[ind];
	pos = point_to_pixel(t, pos);
	if (edind) {
	    int i;
	    for (i = 0; i < 4; i++, pos += .15 * t->scale_x) {
		int v;
		switch(i) {
		case 0: v = t->read->prob_A[ind]; break;
		case 1: v = t->read->prob_C[ind]; break;
		case 2: v = t->read->prob_G[ind]; break;
		case 3: v = t->read->prob_T[ind]; break;
		}
		if (v >= 0)
		    XFillRectangle(d, p, t->ConfGC, pos-fw, cfh+30-v, cfw, v);
		else
		    XFillRectangle(d, p, t->ConfNegGC, pos, cfh+30, 3, -v);
		/*
		if (v >= 0)
		    XFillRectangle(d, p, t->ConfGC, pos-1, cfh, 5, v);
		else
		    XFillRectangle(d, p, t->ConfNegGC, pos, cfh, 3, -v);
		*/
	    }
	}

	ind++;
    }
}

/*
 * Draws a segment of a single trace. We draw 'samples' points from the 'tr'
 * trace array to pixel coordinates (x,y) with scale xs, ys.
 *
 * The confidence display is quite wide, and draws more than x0-xn.
 * To compensate, we need to increase x0-xn.
 */
static void trace_draw2(DNATrace *t, TRACE *tr, Display *d, Pixmap p, int max,
			GC gc, int x0, int xn, int yoff, int height, double ys,
			int off, int sign) {
    int i, h = height-1, o;
    XPoint *xp, *xp2;

    if (xn <= 0 || NULL == (xp = (XPoint *)xmalloc(xn * sizeof(XPoint)))) {
	return;
    }

    o = t->disp_offset * t->scale_x;

    if (t->read->maxTraceVal)
	h -= h*(double)off/t->read->maxTraceVal;
    for (i = 0; i < xn; i++, tr++) {
	xp[i].x = (int)((x0 + i) * t->scale_x) - o;
	if (sign) {
	    xp[i].y = h - ys * ((int_2)*tr-off) + yoff;
	} else {
	    xp[i].y = h - ys * (*tr-off) + yoff;
	}
    }

    for (i = 0; i < xn; i++)
        if (xp[i].x >= 0)
            break;
    xp2 = &xp[i];
    xn -= i;

    if (xn > 0) {
	XDrawLines(d, p, gc, xp2, xn, CoordModeOrigin);
    }

    xfree(xp);
}

/*
 * Draws a segment of a single trace using histogram bars, suitable for
 * pyrosequencing (eg 454, solexa). We draw 'samples' points from the 'tr'
 * trace array to pixel coordinates (x,y) with scale xs, ys.
 *
 * The confidence display is quite wide, and draws more than x0-xn.
 * To compensate, we need to increase x0-xn.
 */
static void trace_draw_hist(DNATrace *t, TRACE *tr, Display *d, Pixmap p,
			    int max, GC gc, double x0, int xn, int yoff,
			    int height, double ys, int off) {
    int i, h = height-1, o;
    int x, y;

    if (xn <= 0)
	return;

    o = t->disp_offset * t->scale_x;

    if (t->read->maxTraceVal)
	h -= h*(double)off/t->read->maxTraceVal;

    for (i = 0; i < xn; i++, tr++) {
	x = (int)((x0 + i) * t->scale_x) - o;
	y = h - ys * (*tr-off) + yoff;
	if (h >= y)
	    XFillRectangle(d, p, gc, x-1, y, 3, h-y);
	else
	    XFillRectangle(d, p, gc, x-1, h, 3, y-h);
    }

}

static void trace_draw_yticks(DNATrace *t, Display *d, Pixmap p, GC gc,
			      int x0, int xn, int yoff, int height,
			      double ys, int off) {
    int i, h = height-1, o;
    int x1, x2, y;

    if (xn <= 0)
	return;

    o = t->disp_offset * t->scale_x;

    y = h - ys * (0*t->yticks-off) + yoff;
    x1 = (int)(x0 * t->scale_x) - o;
    x2 = (int)((x0+xn) * t->scale_x) - o;
    i = 0;
    do {
	XDrawLine(d, p, gc, x1, y, x2, y);
	y = h - ys * (++i*t->yticks-off) + yoff;
    } while (y >= 0);
}

static void trace_draw_pyro(DNATrace *t, Display *d, Pixmap p,
			    int x0, int xn, int yoff,
			    int height, double ys) {
    int i, h = height-1, o;

    if (xn <= 0)
	return;

    o = t->disp_offset * t->scale_x;
    ys = 30;

    printf("drawing %d to %d (%.*s)\n", x0, x0+xn, xn, &t->read->flow_order[x0]);

    for (i = 0; i < xn; i++) {
	int x = (int)((x0 + i) * t->scale_x) - o;
	int inrange = (x0+i >= 0 && x0+i < t->read->nflows);
	float y =  inrange ? t->read->flow[x0+i] : 0;
	char base = inrange ? t->read->flow_order[x0+i] : 0;
	GC gc;
	    
	printf("%d/%d: %c %f\n", x0+i, x, base, y);
	if (y < .05)
	    y = .05;
	y = h - ys * y + yoff;
	switch (base) {
	case 'A':
	case 'a':
	    gc = t->Agc;
	    break;
	case 'C':
	case 'c':
	    gc = t->Cgc;
	    break;
	case 'G':
	case 'g':
	    gc = t->Ggc;
	    break;
	case 'T':
	case 't':
	    gc = t->Tgc;
	    break;
	default:
	    gc = t->CursorGC;
	}
	
	XFillRectangle(d, p, gc, x, y, 3, h-y);
    }
}

#if 0
#define MAX4(t,i) ((t)[0][i].y > (t)[1][i].y \
    ? ((t)[2][i].y > (t)[3][i].y \
        ? ((t)[0][i].y > (t)[2][i].y \
            ? 0 \
            : 2) \
        : ((t)[0][i].y > (t)[3][i].y \
            ? 0 \
            : 3)) \
    : ((t)[2][i].y > (t)[3][i].y \
        ? ((t)[1][i].y > (t)[2][i].y \
            ? 1 \
            : 2) \
        : ((t)[1][i].y > (t)[3][i].y \
            ? 1 \
            : 3)))

#define MIN4(t,i) ((t)[0][i].y < (t)[1][i].y \
    ? ((t)[2][i].y < (t)[3][i].y \
        ? ((t)[0][i].y < (t)[2][i].y \
            ? 0 \
            : 2) \
        : ((t)[0][i].y < (t)[3][i].y \
            ? 0 \
            : 3)) \
    : ((t)[2][i].y < (t)[3][i].y \
        ? ((t)[1][i].y < (t)[2][i].y \
            ? 1 \
            : 2) \
        : ((t)[1][i].y < (t)[3][i].y \
            ? 1 \
            : 3)))

static void trace_draw2filled(DNATrace *t, Display *d, Pixmap p, int max,
			int x0, int xn, int height, double ys,
			int off) {
    int i, h = height-1, o;
    XPoint *xp[4], *xp2[4];
    GC gc[4];
    TRACE *tr[4];
    int ch;

    gc[0] = t->Agc;
    gc[1] = t->Cgc;
    gc[2] = t->Ggc;
    gc[3] = t->Tgc;

    tr[0] = &t->read->traceA[x0];
    tr[1] = &t->read->traceC[x0];
    tr[2] = &t->read->traceG[x0];
    tr[3] = &t->read->traceT[x0];

    for (ch = 0; ch < 4; ch++) {
	if (NULL == (xp[ch] = (XPoint *)xmalloc((xn+2) * sizeof(XPoint)))) {
	    return;
	}
	xp[ch]++;
    }

    o = t->disp_offset * t->scale_x;

    if (t->read->maxTraceVal)
	h -= h*(double)off/t->read->maxTraceVal;
    for (ch = 0; ch < 4; ch++) {
	TRACE *trp = tr[ch];
	XPoint *xpp = xp[ch];
	for (i = 0; i < xn; i++, trp++) {
	    xpp[i].x = (int)((x0 + i) * t->scale_x) - o;
	    xpp[i].y = h - ys * (*trp-off);
	}
	tr[ch] = trp;
    }

    /* Clip off negative X coordinates as this crashes Linux X servers */
    for (i = 0; i < xn; i++) {
	if (xp[0][i].x >= 0 &&
	    xp[1][i].x >= 0 &&
	    xp[2][i].x >= 0 &&
	    xp[3][i].x >= 0)
	    break;
    }
    xp2[0] = &xp[0][i];
    xp2[1] = &xp[1][i];
    xp2[2] = &xp[2][i];
    xp2[3] = &xp[3][i];
    xn -= i;

#if 1
    /* Draw first peak as a line and all lower ones as solid */
    {
	int cur_max;
	int start = 0;

	i = 1;
	cur_max = MIN4(xp2,0);

	do {
	    while (MIN4(xp2,i) == cur_max && i < xn)
		i++;

	    for (ch = 0; ch < 4; ch++) {
		if (ch == cur_max) {
		    XDrawLines(d, p, gc[ch], &xp2[ch][start], i-start,
			       CoordModeOrigin);
		} else {
		    XPoint tmp;
		    tmp = xp2[ch][i];
		    xp2[ch][start-1].x = xp2[ch][start].x;
		    xp2[ch][start-1].y = h;
		    xp2[ch][i].x = xp2[ch][i-1].x;
		    xp2[ch][i].y = h;
		    XFillPolygon(d, p, gc[ch], &xp2[ch][start-1], i-start+2,
				 Nonconvex, CoordModeOrigin);
		    xp2[ch][i] = tmp;
		}
	    }
	    start = i;
	    cur_max = MIN4(xp2,i);
	} while (i < xn);
    }
#endif

#if 0
    /* Draw all peaks as solid */
    for (ch = 0; ch < 4; ch++) {
	/*
	XDrawLines(d, p, gc[ch], xp2[ch], xn,
		   CoordModeOrigin);
	*/
	xp2[ch][0].y = h;
	xp2[ch][xn-1].y = h;
	XFillPolygon(d, p, gc[ch], xp2[ch], xn,
		     Nonconvex, CoordModeOrigin);
    }
#endif

    for (ch = 0; ch < 4; ch++) {
	xp[ch]--;
	xfree(xp[ch]);
    }
}
#endif

/*
 * Draw a trace from pixel coordinates (x0,0) to (x0+xn, height).
 */
void trace_draw_trace(DNATrace *t, Display *d, Pixmap p,
		      int x0, int xn, int yoff, int height) {
    double yscale;
    int m = t->read->maxTraceVal;

    if (x0 < 0) {
	xn += x0;
	x0 = 0;

	if (xn <= 0)
	    return;
    }

    if (x0 + xn > t->read->NPoints)
	xn = t->read->NPoints - x0;

    trace_update_extents(t, &x0, &xn);
    if (xn <= 0)
	return;

    if (t->show_conf) {
	if (t->style == STYLE_STICK)
	    trace_draw_confidence4(t, d, p, x0, xn, height);
	else
	    trace_draw_confidence(t, d, p, x0, xn, height);
    }

    if (t->read->traceA == 0 || p == 0)
	return;

    if (t->read->maxTraceVal) {
	int c = (t->show_conf && t->style == STYLE_STICK) ? 20 : 0;
	yscale = (double)(t->scale_y * (height-1-c)) /
	    (double)(t->trace_scale ? t->trace_scale : t->read->maxTraceVal);
    } else {
	yscale = 0;
    }

    switch (t->style) {
    case STYLE_CHROMA:
	trace_draw2(t, &t->read->traceA[x0], d, p, m,
		    t->Agc, x0, xn, yoff, height, yscale,
		    t->read->baseline, 1);
	trace_draw2(t, &t->read->traceC[x0], d, p, m,
		    t->Cgc, x0, xn, yoff, height, yscale,
		    t->read->baseline, 1);
	trace_draw2(t, &t->read->traceG[x0], d, p, m,
		    t->Ggc, x0, xn, yoff, height, yscale,
		    t->read->baseline, 1);
	trace_draw2(t, &t->read->traceT[x0], d, p, m,
		    t->Tgc, x0, xn, yoff, height, yscale,
		    t->read->baseline, 1);
	break;
#if 0
    case STYLE_FILLED:
	trace_draw2filled(t, d, p, m, x0, xn, height, yscale,
			  t->read->baseline);
	break;
#endif

    case STYLE_PYRO:
	trace_draw_hist(t, &t->read->traceA[x0], d, p, m,
			t->Agc, x0, xn, yoff, height, yscale,
			t->read->baseline);
	trace_draw_hist(t, &t->read->traceC[x0], d, p, m,
			t->Cgc, x0, xn, yoff, height, yscale,
			t->read->baseline);
	trace_draw_hist(t, &t->read->traceG[x0], d, p, m,
			t->Ggc, x0, xn, yoff, height, yscale,
			t->read->baseline);
	trace_draw_hist(t, &t->read->traceT[x0], d, p, m,
			t->Tgc, x0, xn, yoff, height, yscale,
			t->read->baseline);
	break;

    case STYLE_STICK:
	trace_draw_hist(t, &t->read->traceA[x0], d, p, m,
			t->Agc, x0, xn, yoff, height, yscale,
			t->read->baseline);
	trace_draw_hist(t, &t->read->traceC[x0], d, p, m,
			t->Cgc, x0+.15, xn, yoff, height, yscale,
			t->read->baseline);
	trace_draw_hist(t, &t->read->traceG[x0], d, p, m,
			t->Ggc, x0+.3, xn, yoff, height, yscale,
			t->read->baseline);
	trace_draw_hist(t, &t->read->traceT[x0], d, p, m,
			t->Tgc, x0+.45, xn, yoff, height, yscale,
			t->read->baseline);
	/* trace_draw_pyro(t, d, p, x0, xn, yoff, height, yscale); */
	break;
    }

    if (t->yticks)
	trace_draw_yticks(t, d, p, t->CursorGC, x0, xn, yoff, height,
			  yscale, t->read->baseline);

    if (t->show_edits == 0) {
	int pos;

	pos = point_to_pixel(t, trace_get_pos(t, t->cursor_pos));

	XFillRectangle(d, p, t->CursorGC, pos-1, yoff, 1, height);
    }
}

/*
 * Draw the sequence from pixel coordinates (x,y) to (x+width, y+height).
 */
void trace_draw_sequence(DNATrace *t, Display *d, Pixmap p,
			 int x0, int xn, int yoff, int height) {
    int ind, pos, x1, fw, fh;

    if( !p || !t || !t->read || (t->read->NBases==0) )
	return;

    x1 = x0 + xn < t->read->NPoints ? x0+xn : t->read->NPoints - 1;
    x1 = t->tracePos[x1] + 1;
    if (x1 >= t->read->NBases) x1--;
    x1 = t->read->basePos[x1];
    ind = t->tracePos[x0];       /* Index of first base on screen */

    fw = t->font_width/2 + 1;
    fh = t->fm.ascent + yoff;

    while (ind < t->read->NBases && (pos = t->read->basePos[ind]) <= x1) {
	char base = t->read->base[ind];
	GC gc;
	double xoff = 0;

	switch (base) {
	case 'A':
	case 'a':
	    gc = t->Agc;
	    break;
	case 'C':
	case 'c':
	    gc = t->Cgc;
	    xoff = 0.15;
	    break;
	case 'G':
	case 'g':
	    gc = t->Ggc;
	    xoff = 0.3;
	    break;
	case 'T':
	case 't':
	    gc = t->Tgc;
	    xoff = 0.45;
	    break;
	default:
	    gc = t->CursorGC;
	}

	/* XDrawString(d, p, gc, point_to_pixel(t, pos) - fw, fh, &base, 1); */
	if (t->style != STYLE_STICK)
	    xoff = 0;
	pos = point_to_pixel(t, (double)pos + xoff) - fw;
	Tk_DrawChars(d, p, gc, t->font, &base, 1, pos, fh);

	ind++;
    }
}

/*
 * Draw the numbers from pixel coordinates (x,y) to (x+width, y+height).
 */
void trace_draw_numbers(DNATrace *t, Display *d, Pixmap p,
			int x0, int xn, int yoff, int height) {
    int ind, pos, x1, fh, i;
    float fwid, fw;
    char number[10];
    double xoff;

    if (!p)
	return;

    x1 = x0 + xn < t->read->NPoints ? x0+xn : t->read->NPoints - 1;
    x1 = t->tracePos[x1] + 1;
    if (x1 >= t->read->NBases) x1--;
    x1 = t->read->basePos[x1];

    fw = t->font_width/2 + 1;

    x0 -= 2*fw;
    if (x0 < 0) x0 = 0;
    ind = t->tracePos[x0];       /* Index of first base on screen */
    if (ind < 1) ind = 1;

    fwid = t->font_width;
    fh = t->fm.ascent + yoff;

    while (ind < t->read->NBases && (pos = t->read->basePos[ind-1]) <= x1) {
	int cind = t->comp ? t->read->NBases - ind + 1: ind;
	if (cind % 10 == 0) {
	    /* Base numbers are right aligned */
	    i = ABS(cind);
	    if (i >= 1000) {
		fw = fwid * 3.5;
	    } else if (i >= 100) {
		fw = fwid * 2.5;
	    } else if (i >= 10) {
		fw = fwid * 1.5;
	    } else {
		fw = fwid * .5;
	    }

	    /* Stick plots have A,C,G,T separated out, so move the number */
	    xoff = 0;
	    if (t->style == STYLE_STICK) {
		switch (t->read->base[ind-1]) {
		case 'A': case 'a': xoff = 0.00; break;
		case 'C': case 'c': xoff = 0.15; break;
		case 'G': case 'g': xoff = 0.30; break;
		case 'T': case 't': xoff = 0.45; break;
		}
	    }

	    sprintf(number, "%d", cind);
	    Tk_DrawChars(d, p, t->CursorGC, t->font, number, strlen(number),
			 point_to_pixel(t, pos + xoff) - fw, fh);
	}

	ind++;
    }
}

/*
 * Find the position of an edited base
 */
#define tc(x) (t->comp ? t->read->NBases - (x) - 1 : (x))
int trace_get_pos(DNATrace *t, int ind) {
    int back, forw;
    int backp, forwp;
    double avg;

    if (t->Ned <= 0)
	return 0;

    avg = (t->read->basePos[t->read->NBases-1] - t->read->basePos[0]) /
      (double)t->read->NBases;
    if (ind < 0) {
	/* Interpolate backwards */
	return trace_get_pos(t, 0) + ind*avg;
    }

    if (ind >= t->Ned) {
	/* Interpolate forwards */
	return trace_get_pos(t, t->Ned-1) + (ind - (t->Ned-1))*avg;
    }

    /* The simple case */
    if (t->edPos[ind])
	return t->read->basePos[t->edPos[tc(ind)]-1];

    for (back = ind-1; back >= 0 && t->edPos[back] == 0; back--)
	;
    if (back < 0)
	back = 0;

    for (forw = ind+1; forw < t->Ned && t->edPos[forw] == 0; forw++)
	;

    if (t->edPos[forw])
	forwp = t->read->basePos[t->edPos[tc(forw)]-1];
    else
	forwp = t->read->NPoints;

    if (t->edPos[back])
	backp = t->read->basePos[t->edPos[tc(back)]-1];
    else
	backp = forwp / 4;


/*
    printf("ind = %d, back = %d, forw = %d, backp = %d, forwp = %d res = %d\n",
	   ind, back, forw, backp, forwp,
	   (forwp - backp) * (ind - back) / (forw - back) + backp);
*/
    return (forwp - backp) * (ind - back) / (forw - back) + backp;
}

/*
 * Draw the edited seq. from pixel coordinates (x,y) to (x+width, y+height).
 */
void trace_draw_edits(DNATrace *t, Display *d, Pixmap p,
		      int x0, int xn, int yoff, int height) {
    int ind, pos, x1, fw, fh;

    if( !p || !t || !t->read || (t->read->NBases==0) )
	return;

    x0 -= 4;
    if (x0 < 0) x0 = 0;
    xn += 8;
    x1 = x0 + xn < t->read->NPoints ? x0+xn : t->read->NPoints - 1;
    x1 = t->tracePos[x1] + 1;
    if (x1 >= t->read->NBases) x1--;
    x1 = t->read->basePos[x1];
    ind = t->tracePosE[x0];       /* Index of first base on screen */

    fw = t->font_width/2 + 1;
    fh = t->fm.ascent + yoff;

    while (ind < t->Ned && (pos = trace_get_pos(t, ind)) <= x1) {
	char base;
	GC gc;

	base = t->edBases[ind];

	switch (base) {
	case 'A':
	case 'a':
	    gc = t->Agc;
	    break;
	case 'C':
	case 'c':
	    gc = t->Cgc;
	    break;
	case 'G':
	case 'g':
	    gc = t->Ggc;
	    break;
	case 'T':
	case 't':
	    gc = t->Tgc;
	    break;
	default:
	    gc = t->CursorGC;
	}

	/* XDrawString(d, p, gc, point_to_pixel(t, pos) - fw, fh, &base, 1); */
	Tk_DrawChars(d, p, gc, t->font, &base, 1,
		     point_to_pixel(t, pos) - fw, fh);
	ind++;
    }

    /* Cursor */
    if (t->cursor_pos > 0)
	pos = trace_get_pos(t, t->cursor_pos-1);
    else
	pos = 0;

    XFillRectangle(d, p, t->CursorGC, point_to_pixel(t, pos)+4,
		   t->fm.linespace-3, 8, 3);
}

void trace_flash (DNATrace *t) {
    int x0, i;
    Display *d = t->display;
    Window w = Tk_WindowId(t->tkwin);
    Pixmap tmp;

    if (!Tk_IsMapped(t->tkwin) || Tk_WindowId(t->tkwin) == 0)
	return;

    x0 = point_to_pixel(t, trace_get_pos(t, t->cursor_pos));

    tmp = Tk_GetPixmap(d, w, 24, t->pos[TRACEP_T].h, Tk_Depth(t->tkwin));
    XCopyArea(d, w, tmp, t->CursorGC,
	      x0-12, t->pos[TRACEP_T].y, 24, t->pos[TRACEP_T].h, 0, 0);

    for (i=12; i > 0; i-=3) {
	XCopyArea(d, tmp, w, t->CursorGC,
		  0, 0, 24, t->pos[TRACEP_T].h, x0-12, t->pos[TRACEP_T].y);

	XFillRectangle(d, w,
		       t->CursorGC,
		       x0 - i, t->pos[TRACEP_T].y,
		       i, t->pos[TRACEP_T].h);

	XSync(d, False);

	myusleep(20000);
    }
    XCopyArea(d, tmp, w, t->CursorGC,
	      0, 0, 24, t->pos[TRACEP_T].h, x0-12, t->pos[TRACEP_T].y);


    Tk_FreePixmap(d, tmp);
}
