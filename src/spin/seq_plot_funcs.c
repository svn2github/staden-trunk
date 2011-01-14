#include <tcl.h>
#include <stdlib.h>

#include "xalloc.h"
#include "seq_raster.h"
#include "tkRaster.h"
#include "text_output.h"

/*
static int b_compare(const void *p1,
		     const void *p2)
{
    int *i1 = (int *) p1;
    p_score *i2 = (p_score *) p2;

    return ((*i1) - i2->pos);

}
*/

static int b_search(int key,
		    p_score *array,
		    int num_pts)
{
    int found = 0;
    int middle, low, high;

    low = 0;
    high = num_pts-1;

    while (!found) {

	middle = (low + high) / 2;

	if (key < array[middle].pos) {
	    high = middle - 1;
	} else if (key > array[middle].pos) {
	    low = middle + 1;
	} else {
	    found = 1;
	    return middle;
	}
	if (low > high) {
	    break;
	}
    }
    return middle;
}

void graph_plot_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_raster *output = result->output;
    r_graph *data = result->data;
    int num_pts = data->n_pts;
    double sf_m = output->sf_m;
    double sf_c = output->sf_c;
    int i, j;
    double coords[2]; 
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    p_score prev;
    double x0, y0, x1, y1;
    int first = -1;
    int last = -1;

    if (output->hidden) {
	return;
    }
    Tcl_GetCommandInfo(output->interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;
    SetDrawEnviron(output->interp, raster, output->env_index);
    RasterGetWorldScroll(raster, &x0, &y0, &x1, &y1);

    if (num_pts == 1) {
	coords[0] = data->p_array[0].pos;
	coords[1] = y1 - (sf_m * data->p_array[0].score + sf_c);
	RasterDrawPoints(raster, coords, 1);
    } else { 
	double *coords = xmalloc(sizeof(*coords) *
				 (2 * (num_pts + 4)));

	/*
	 * it isn't necessarily the case that I will have as many points as
	 * sequence bases so need to check if each base is within my range
	 */

	/* need to find array indices for first and last points */
	if (data->p_array[0].pos >= plot->x0) {
	    first = 0;
	} else {
	    first = b_search(plot->x0, data->p_array, num_pts);
	    if (first > 0)
		first--;
	}

	last = b_search(plot->x1, data->p_array, num_pts);

	last+=2;
	if (last > num_pts)
	    last = num_pts;
#ifdef DEBUG
	printf("first %d %d last %d %d\n", first, data->p_array[first].pos,
	       last, data->p_array[last-1].pos);
#endif
	prev = data->p_array[first];

	coords[j=0] = prev.pos;
	coords[j+1] = y1 - (sf_m * prev.score + sf_c) + y0;
	j += 2;

	for (i = first+1; i < last; i++, j+=2) {
#ifdef DEBUG    
	    printf("i %d prev1 %d score %f\n", i, data->p_array[i].pos, 
		   data->p_array[i].score);
#endif
	    coords[j] = data->p_array[i].pos;
	    coords[j+1] = y1 - (sf_m * data->p_array[i].score + sf_c) + y0;
	}
	RasterDrawLines(raster, coords, last-first);
	xfree(coords);
    }
}

void emboss_graph_plot_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_raster *output = result->output;
    e_graph *data = result->data;
    int num_pts = data->n_pts;
    double sf_m = output->sf_m;
    double sf_c = output->sf_c;
    int i, j;
    double *coords; 
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    p_score prev;
    double x0, y0, x1, y1;
    int first = -1;
    int last = -1;
    double coords_4[4];

#ifdef DEBUG
    printf("emboss_graph_plot_func\n");
#endif

    if (output->hidden) {
	return;
    }
    Tcl_GetCommandInfo(output->interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;
    SetDrawEnviron(output->interp, raster, output->env_index);
    RasterGetWorldScroll(raster, &x0, &y0, &x1, &y1);

    if (num_pts == 1) {
	coords = xmalloc(sizeof(*coords) * 2);
	coords[0] = data->p_array[0].pos;
	coords[1] = y1 - (sf_m * data->p_array[0].score + sf_c);
	RasterDrawPoints(raster, coords, 1);
    } else if (num_pts > 1) { 
	coords = xmalloc(sizeof(*coords) *
			 (2 * (num_pts + 4)));

	/*
	 * it isn't necessarily the case that I will have as many points as
	 * sequence bases so need to check if each base is within my range
	 */

	/* need to find array indices for first and last points */
	if (data->p_array[0].pos >= plot->x0) {
	    first = 0;
	} else {
	    first = b_search(plot->x0, data->p_array, num_pts);
	    if (first > 0)
		first--;
	}

	last = b_search(plot->x1, data->p_array, num_pts);
	last+=2;
	if (last > num_pts)
	    last = num_pts;
#ifdef DEBUG
	printf("first %d %d last %d %d\n", first, data->p_array[first].pos,
	       last, data->p_array[last-1].pos);
#endif
	prev = data->p_array[first];

	coords[j=0] = prev.pos;
	coords[j+1] = y1 - (sf_m * prev.score + sf_c) + y0;
	j += 2;

	for (i = first+1; i < last; i++, j+=2) {
#ifdef DEBUG    
	    printf("prev1 %d score %f\n", data->p_array[i].pos, 
		   data->p_array[i].score);
#endif
	    coords[j] = data->p_array[i].pos;
	    coords[j+1] = y1 - (sf_m * data->p_array[i].score + sf_c) + y0;
#ifdef DEBUG
	    printf("coords %f %f\n", coords[j], coords[j+1]);
#endif
	}
	RasterDrawLines(raster, coords, last-first);
	xfree(coords);
    }
    
    /* do data and graph objects */

    for (i = 0; i < data->n_data_obj; i++) {
	if (data->d_obj[i].type == E_LINE) {
	    RasterDrawLine(raster, (int)data->d_obj[i].pos.x0, 
			   y1 - (sf_m * data->d_obj[i].pos.y0 + sf_c) + y0,
			   (int)data->d_obj[i].pos.x1, 
			   y1 - (sf_m * data->d_obj[i].pos.y1 + sf_c) + y0);
	} else if (data->d_obj[i].type == E_RECTANGLE) {
#ifdef DEBUG
	    printf("data rectangle \n");
#endif
	    coords_4[0] = data->d_obj[i].pos.x0;
	    coords_4[1] = y1 - (sf_m * data->d_obj[i].pos.y0 + sf_c) + y0;
	    coords_4[2] = data->d_obj[i].pos.x1;
	    coords_4[3] = y1 - (sf_m * data->d_obj[i].pos.y1 + sf_c) + y0;
	    RasterDrawRectangles(raster, coords_4, 1);
	} else if (data->d_obj[i].type == E_RECTANGLEFILL) {
#ifdef DEBUG
	    printf("data rectangle filled\n");
#endif
	    coords_4[0] = data->d_obj[i].pos.x0;
	    coords_4[1] = y1 - (sf_m * data->d_obj[i].pos.y0 + sf_c) + y0;
	    coords_4[2] = data->d_obj[i].pos.x1;
	    coords_4[3] = y1 - (sf_m * data->d_obj[i].pos.y1 + sf_c) + y0;
	    RasterFillRectangles(raster, coords_4, 1);
	} else if (data->d_obj[i].type == E_TEXT) {
#ifdef DEBUG
	    printf("data text \n");
#endif
	}

    }

    /* 
     * need to do plotting one item at a time because I can't guarentee all
     * objects will be of the same type
     */
    for (i = 0; i < data->n_graph_obj; i++) {
	if (data->g_obj[i].type == E_LINE) {
	    RasterDrawLine(raster, (int)data->g_obj[i].pos.x0, 
			   y1 - (sf_m * data->g_obj[i].pos.y0 + sf_c) + y0,
			   (int)data->g_obj[i].pos.x1, 
			   y1 - (sf_m * data->g_obj[i].pos.y1 + sf_c) + y0);
	} else if (data->g_obj[i].type == E_RECTANGLE) {
#ifdef DEBUG
	    printf("graph rectangle \n");
#endif
	    coords_4[0] = data->g_obj[i].pos.x0;
	    coords_4[1] = y1 - (sf_m * data->g_obj[i].pos.y0 + sf_c) + y0;
	    coords_4[2] = data->g_obj[i].pos.x1;
	    coords_4[3] = y1 - (sf_m * data->g_obj[i].pos.y1 + sf_c) + y0;
	    RasterDrawRectangles(raster, coords_4, 1);
	} else if (data->g_obj[i].type == E_RECTANGLEFILL) {
	    coords_4[0] = data->g_obj[i].pos.x0;
	    coords_4[1] = y1 - (sf_m * data->g_obj[i].pos.y0 + sf_c) + y0;
	    coords_4[2] = data->g_obj[i].pos.x1;
	    coords_4[3] = y1 - (sf_m * data->g_obj[i].pos.y1 + sf_c) + y0;
#ifdef DEBUG
	    printf("graph rectangle filled %f %f %f %f\n", coords_4[j],
		    coords_4[j+1],  coords_4[j+2],  coords_4[j+3]);
#endif
	    RasterFillRectangles(raster, coords_4, 1);
	    
	} else if (data->g_obj[i].type == E_TEXT) {
#ifdef DEBUG
	    printf("graph text \n");
#endif
	}
    }
}

/*
 * A more specific type of graph plot to account for the "top" line
 * required for the gene search functions
 */
void gene_search_plot_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_raster *output = result->output;
    gene_search *data = result->data;
    int num_pts = data->n_pts;
    double sf_m = output->sf_m;
    double sf_c = output->sf_c;
    int type = result->type;
    int i, j;
    double coords[2]; 
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    p_score prev;
    double x0, y0, x1, y1;
    int plot_x0, plot_x1;
    double middle;
    double y_76;

    if (output->hidden) {
	return;
    }

    Tcl_GetCommandInfo(output->interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;
    SetDrawEnviron(output->interp, raster, output->env_index);
    RasterGetWorldScroll(raster, &x0, &y0, &x1, &y1);

    /* 
     * convert x world coords into an index into the array containing the
     * y values ie subtract the x "start" value and divide by 3
     */
    if (plot->x0 <= data->dim.x0) {
	plot_x0 = 0;
    } else {
	plot_x0 = ((plot->x0 - data->dim.x0) / 3) - 1;
    }

    if (plot->x0 > data->dim.x1) {
	plot_x0 = num_pts;
    }

    if (plot_x0 >= num_pts)
	plot_x0 = num_pts-1;

    if (plot->x1 >= data->dim.x1) {
	plot_x1 = num_pts;
    } else {
	plot_x1 = ((plot->x1 - data->dim.x0) / 3) + 2;
    }

    /* at end of sequence, adding 2 can be too much */
    if (plot_x1 > num_pts)
	plot_x1 = num_pts;
    
    if (plot_x1 < 0) {
	plot_x1 = 0;
    }

#ifdef DEBUG    
    printf("plot_x0 %d plot_x1 %d\n", plot_x0, plot_x1);
#endif
    middle = (y0 + y1) / 2;
    if (num_pts == 1) {
	coords[0] = data->p_array[0].pos;
	coords[1] = y1 - (sf_m * data->p_array[0].score + sf_c);
	RasterDrawPoints(raster, coords, 1);
    } else { 
	/*
	 * Buffer lines and points into single X requests.
	 * This uses about 68% of the CPU and 41% of the X traffic compared
	 * to separate requests.
	 */
	double *coords = xmalloc(sizeof(*coords) *
				 (2 * (plot_x1 - plot_x0 + 2)));

	/* 7/1/99 johnt - use a single line rather than many points, to show max peaks, as Windows is VERY slow at drawing XPoints */
#ifdef USEPOINTS
	double *points = xmalloc(sizeof(*points) *
				 (2 * (plot_x1 - plot_x0 + 2)));
#else
	double *points = xmalloc(sizeof(*points) *4);
#endif
	int npoints = 0;

	prev = data->p_array[plot_x0];
	coords[j=0] = prev.pos;
	coords[j+1] = y1 - (sf_m * prev.score + sf_c) + y0;
	j += 2;
	for (i = plot_x0; i < plot_x1; i++, j+=2) {
	    coords[j] = data->p_array[i].pos;
	    coords[j+1] = y1 - (sf_m * data->p_array[i].score + sf_c) + y0;

	    /* 7/1/99 johnt - use a single line rather than many points, to show max peaks, as Windows is VERY slow at drawing XPoints */
	    /* draw dot in middle of graph to denote max peak */
	    if (data->top) {
		if (result->frame == data->top[i])
#ifdef USEPOINTS
		{
		    points[npoints] = data->p_array[i].pos;
		    points[npoints+1] = middle;
		    npoints += 2;
		}
#else
		{
		    switch(npoints){
		    case 0: /* got no points yet */
		    case 2: /* got first point, but not last */
    			points[npoints] = data->p_array[i].pos;
			points[npoints+1] = middle;
			npoints += 2;
			break;
		    case 4: /* got an end point, so replace it */
			points[2] = data->p_array[i].pos;
			break;
		    }
		}
		else{
		    switch(npoints){
		    case 0: /* no lines to draw */
			break;
		    case 2: /* a single point to draw */
			points[2] = points[0];
			points[3] = middle;
			/* fall through */
		    case 4: /* a line to draw */
			RasterDrawLine(raster,(int)points[0],points[1],(int)points[2]+1,points[3]);
			npoints=0;
			break;
		    }
		}
#endif
	    }
	}

	RasterDrawLines(raster, coords, plot_x1 - plot_x0 + 1);

        /* 7/1/99 johnt - use a single line rather than many points, to show max peaks, as Windows is VERY slow at drawing XPoints */
#ifdef USE_POINTS
	RasterDrawPoints(raster, points, npoints/2);
#else
	switch(npoints){
	case 0: /* no lines to draw */
	    break;
	case 2: /* draw to the end of the display */
    	    points[2] = data->p_array[i-1].pos;
	    points[3] = middle;
	    /* fall through */
	case 4: /* a line to draw */
	    RasterDrawLine(raster,(int)points[0],points[1],(int)points[2],points[3]);
	    npoints=0;
	    break;
	}
#endif
	xfree(coords);
	xfree(points);
    }

    /* only want to draw if displaying base bias results */
    /* turn 0.76 upside down */
    if (type & SEQ_TYPE_BASEBIAS) {
      y_76 = y1 - 0.76 + y0;
      RasterDrawLine(raster, (int)data->dim.x0, y_76, 
		     (int)data->dim.x1, y_76); 
    }
}

/* eg splice junction plot where there are 2 stick plots, of 1/2 height */
void stick_pair_plot_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_raster *output = result->output;
    stick *data = result->data;
    config **configure = output->configure;
    int i, j;
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    double sx0, sy0, sx1, sy1;
    double y0, y1;
    double tick_ht;
    int plot_x0, plot_x1;
    int raster_id;
    int x;
    double score, scaled_score;
    RasterResult *raster_result;
    double m = 0, c = 0;

    if (output->hidden) {
	return;
    }
    Tcl_GetCommandInfo(output->interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;

    SetDrawEnviron(output->interp, raster, output->env_index);

    /* world scroll limits */
    RasterGetWorldScroll(raster, &sx0, &sy0, &sx1, &sy1);
    Tcl_VarEval(output->interp, "GetRasterId ", output->raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(output->interp));
    raster_result = raster_id_to_result(raster_id);

    
    if (plot->x0 > data->ap_array[0].dim.x0) {
	plot_x0 = plot->x0;
    } else {
	plot_x0 = data->ap_array[0].dim.x0;
    }

    if (plot->x1 < data->ap_array[0].dim.x1) {
	plot_x1 = plot->x1;
    } else {
	plot_x1 = data->ap_array[0].dim.x1;
    }

    for (j = 0; j < data->n_pts; j++) {

	FindRasterResultY0(raster, raster_id, configure[j], 
			   raster_result->num_results, &y0, &tick_ht);

	if (j > 0) {
	    /* 
	     * scaling factors to ensure all plots are on the same scale as
	     * the first plot
	     */
	    m = (data->ap_array[0].dim.y0 - data->ap_array[0].dim.y1) / 
		(data->ap_array[j].dim.y0 - data->ap_array[j].dim.y1);
	    c = data->ap_array[0].dim.y1 - (data->ap_array[j].dim.y1 * m);
	}
	for (i = 0; i < data->ap_array[j].n_pts; i++) {
	    x = data->ap_array[j].p_array[i].pos;
	    score = (double) data->ap_array[j].p_array[i].score;

	    /* need to scale the score to fit on the donor scale */
	    if (j == 0) {
		scaled_score = score;
	    } else {
		scaled_score = m * (double)score + c;
	    }

	    if (x >= plot_x0 && x <= plot_x1) {
#ifdef DEBUG
		printf("DrawLine %d x %d score %f\n", 
		       i, x, (double)score);
#endif
		if ((configure[j]->zoom == 1 && 
		     raster_result->num_results == 1) ||
		    configure[j]->zoom == 2) {
		    if (configure[j]->y_direction == '+') {
			y1 = scaled_score + sy0;
			/* printf("1 y1 %f\n", y1); */
		    } else {
			y1 = scaled_score + sy0;
			y1 = sy1 - y1 + sy0;
			/* printf("2 y1 %f\n", y1); */
		    }
		} else {
		    if (configure[j]->y_direction == '+') {
			y1 = y0 + (tick_ht * (score-data->ap_array[j].dim.y0)/ 
				   ((data->ap_array[j].dim.y1/2) - 
				    data->ap_array[j].dim.y0));
			/* printf("3 y1 %f\n", y1); */
		    } else {
			y1 = y0 - (tick_ht * (score-data->ap_array[j].dim.y0)/ 
				   ((data->ap_array[j].dim.y1/2) - 
				    data->ap_array[j].dim.y0));
			/* printf("4 y1 %f\n", y1); */
		    }
		}
#ifdef DEBUG
		printf("x %d y0 %f y1 %f\n", x, sy1 - y0 + sy0,sy1 - y1 + sy0);
#endif
		RasterDrawLine(raster, x, sy1 - y0 + sy0, x, sy1 - y1 + sy0);
	    }
	}
    }
}

/* 
 * draws sticks eg weight matrix search
 */
void stick_plot_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_raster *output = result->output;
    stick *data = result->data;
    config **configure = output->configure;
    int i, j;
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    double sx0, sy0, sx1, sy1;
    double y0, y1;
    double tick_ht;
    int plot_x0, plot_x1;
    int raster_id;
    int x;
    double score, scaled_score;
    RasterResult *raster_result;
    double m = 0, c = 0;

    if (output->hidden) {
	return;
    }
    Tcl_GetCommandInfo(output->interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;

    SetDrawEnviron(output->interp, raster, output->env_index);

    /* world scroll limits */
    RasterGetWorldScroll(raster, &sx0, &sy0, &sx1, &sy1);
    Tcl_VarEval(output->interp, "GetRasterId ", output->raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(output->interp));
    raster_result = raster_id_to_result(raster_id);

    
    if (plot->x0 > data->ap_array[0].dim.x0) {
	plot_x0 = plot->x0;
    } else {
	plot_x0 = data->ap_array[0].dim.x0;
    }

    if (plot->x1 < data->ap_array[0].dim.x1) {
	plot_x1 = plot->x1;
    } else {
	plot_x1 = data->ap_array[0].dim.x1;
    }

    for (j = 0; j < data->n_pts; j++) {

	FindRasterResultY0(raster, raster_id, configure[j], 
			   raster_result->num_results, &y0, &tick_ht);

	if (j > 0) {
	    /* 
	     * scaling factors to ensure all plots are on the same scale as
	     * the first plot
	     */
	    m = (data->ap_array[0].dim.y0 - data->ap_array[0].dim.y1) / 
		(data->ap_array[j].dim.y0 - data->ap_array[j].dim.y1);
	    c = data->ap_array[0].dim.y1 - (data->ap_array[j].dim.y1 * m);
	}

	for (i = 0; i < data->ap_array[j].n_pts; i++) {
	    x = data->ap_array[j].p_array[i].pos;
	    score = (double) data->ap_array[j].p_array[i].score;

	    /* need to scale the score to fit on the donor scale */
	    if (j == 0) {
		scaled_score = score;
	    } else {
		scaled_score = m * (double)score + c;
	    }

	    if (x >= plot_x0 && x <= plot_x1) {
#ifdef DEBUG
		printf("DrawLine %d x %d score %f\n", 
		       i, x, (double)score);
#endif
		if ((configure[j]->zoom == 1 && 
		     raster_result->num_results == 1) ||
		    configure[j]->zoom == 2) {
		    if (configure[j]->y_direction == '+') {
			/* y1 = (scaled_score + sy0);  */
			y1 = scaled_score;
		    } else {
			/* y1 = (scaled_score + sy0); */
			y1 = scaled_score;
			y1 = sy1 - y1 + sy0;
		    }
		} else {
		    if (configure[j]->y_direction == '+') {
			y1 = y0 + (tick_ht * (score-data->ap_array[j].dim.y0)/ 
				   (data->ap_array[j].dim.y1 - 
				    data->ap_array[j].dim.y0));
		    } else {
			y1 = y0 - (tick_ht * (score-data->ap_array[j].dim.y0)/ 
				   (data->ap_array[j].dim.y1 - 
				    data->ap_array[j].dim.y0));
		    }
		}
#ifdef DEBUG
		printf("x %d y0 %f y1 %f\n", x, sy1 - y0 + sy0,sy1 - y1 + sy0);
#endif
		RasterDrawLine(raster, x, sy1 - y0 + sy0, x, sy1 - y1 + sy0);
	    }
	}
    }
}

/* 
 * draws a single dot at the midpt of a line 
 * used in similar spans plot
 */
void dot_plot_middot_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_raster *output = result->output;
    d_plot *data = result->data;
    int num_pts = data->n_pts;
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    double x0, y0, x1, y1;
    int mid_pt;
    int i;
    double coords[2]; 
  
    if (output->hidden) {
	return;
    }
    Tcl_GetCommandInfo(output->interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;
    SetDrawEnviron(output->interp, raster, output->env_index);

    RasterGetWorldScroll(raster, &x0, &y0, &x1, &y1);

    /* ensure pts off the edge are plotted at the edge */
    for (i = 0; i < num_pts; i++) {
	mid_pt = (int) (data->win_len / 2);

	if ((data->p_array[i].x + mid_pt) > x1) {
	    coords[0] = x1;
	} else if ((data->p_array[i].x + mid_pt) < 1) {
	    coords[0] = 1;
	} else {
	    coords[0] = data->p_array[i].x + mid_pt;
	}
	if (rasterY(raster, data->p_array[i].y + mid_pt) < 1) {
	    coords[1] = 1;
	} else if (rasterY(raster, data->p_array[i].y + mid_pt) > (y1 - 1)) {
	    coords[1] = y1 - 1;
	} else {
	    coords[1] = rasterY(raster, data->p_array[i].y + mid_pt);
	}
#ifdef DEBUG
	printf("%d %f y1=%f\n", data->p_array[i].x+mid_pt, data->p_array[i].y+mid_pt, y1);
	printf("%f %f \n", coords[0], coords[1]);
#endif

	RasterDrawPoints(raster, coords, 1);

    }
}

/* 
 * draws a single dot at x,y
 */
void dot_plot_dot_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_raster *output = result->output;
    d_plot *data = result->data;
    int num_pts = data->n_pts;
    int i,j;
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    double wx0, wy0, wx1, wy1;
    double *points;
      
    if (output->hidden) {
	return;
    }
    Tcl_GetCommandInfo(output->interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;

    SetDrawEnviron(output->interp, raster, output->env_index);
    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

    if(!(points = malloc(sizeof(double) * 2 *num_pts))) 
	return;
  
    /* draw line from x1, y1 to x1+score-1, y1+score-1 */
    for (i = 0, j = 0; i < num_pts; i++, j+=2) {
	points[j] = data->p_array[i].x;
	points[j+1] = rasterY(raster, data->p_array[i].y);
    }
    RasterDrawPoints(raster, points, num_pts);
    free(points);
    
    tk_RasterRefresh(raster);
}

/*
 * draws lines from x1, y1 to x2, y2
 */
void dot_plot_line_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_raster *output = result->output;
    d_plot *data = result->data;
    int num_pts = data->n_pts;
    int i,j;
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    double wx0, wy0, wx1, wy1;
    double coords[2]; 

    if (output->hidden) {
	return;
    }
    Tcl_GetCommandInfo(output->interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;

    SetDrawEnviron(output->interp, raster, output->env_index);
    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

    if (num_pts == 1) {
	coords[0] = data->p_array[0].x;
	coords[1] = (int)rasterY(raster, data->p_array[0].y);
	RasterDrawPoints(raster, coords, 1);
    } else {
	/* 21/1/99 johnt - modified to use RasterDrawLines for performance on WINNT */
        double* points=malloc(sizeof(double) * num_pts * 2);
	for(i = 0, j = 0; i < num_pts; i++){
	    /* found delimiter between alignments, draw everything upto here */
	    if (data->p_array[i].x == -1 && data->p_array[i].y == -1 && 
		data->p_array[i].score == -1){
		switch(j){
		case 0: /* nothing to draw */
		    break;
		case 2: /* draw last point as a dot */
		    RasterDrawPoints(raster, points, 1);
		    break;
		default: /* draw lines upto last point */
    		    RasterDrawLines(raster,points,j/2);
		}
		j=0; /* start line again */
		continue; /* skip this point */
	    }
	    /* add this point to the list */
	    points[j++]=data->p_array[i].x;
    	    points[j++]=rasterY(raster, data->p_array[i].y);
	}
	switch(j){
	case 0: /* nothing to draw */
	    break;
	case 2: /* draw last point as a dot */
	    RasterDrawPoints(raster, points, 1);
	    break;
	default: /* draw lines upto last point */
    	    RasterDrawLines(raster, points, j/2);
	}   
	free(points);
    }
    tk_RasterRefresh(raster);
}

/*
 * draws lines from x1, y1 to x1+score-1, y1+score-1
 */
void dot_plot_scoreline_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_raster *output = result->output;
    d_plot *data = result->data;
    int num_pts = data->n_pts;
    int i,j;
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    double wx0, wy0, wx1, wy1;
    double *points;

    if (output->hidden) {
	return;
    }
    Tcl_GetCommandInfo(output->interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;

    SetDrawEnviron(output->interp, raster, output->env_index);
    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

/* 29/1/99 johnt - use RasterDrawSegments for performance on WINNT */
#ifdef USE_DRAW_LINE
    /* draw line from x1, y1 to x1+score-1, y1+score-1 */
    for (i = 0; i < num_pts; i++) {
	RasterDrawLine(raster, data->p_array[i].x, 
		       rasterY(raster, data->p_array[i].y),
		       data->p_array[i].x + (data->p_array[i].score - 1), 
		       rasterY(raster, data->p_array[i].y + 
			       (data->p_array[i].score - 1)));
    }
#else
    points = malloc(sizeof(double) * num_pts * 4);
    for(i = 0, j = 0; i < num_pts; i++){
	points[j++] = data->p_array[i].x;
	points[j++] = rasterY(raster, data->p_array[i].y);
	points[j++] = data->p_array[i].x + (data->p_array[i].score - 1);
	points[j++] = rasterY(raster, data->p_array[i].y + 
			    (data->p_array[i].score - 1));
    }
    RasterDrawSegments(raster, points, num_pts);
    free(points);
#endif
    tk_RasterRefresh(raster);
}
