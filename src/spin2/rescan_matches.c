#include <string.h>

#include "xalloc.h"
#include "seq_results.h"
#include "sip_results.h"
#include "readpam.h"
#include "dna_utils.h"
#include "sip_similar_spans.h"

#if 0
void SipRescanMatches(Tcl_Interp *interp,
		      seq_result *s_result,
		      int id,
		      int min_score) 
{
    int i, j;
    double coords[2];
    int index;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int score;
    int start;
    int x, y;
    int index1, index2;
    double wx0, wy0, wx1, wy1;

    if (output->hidden) {
	return;
    }

    index1 = GetSeqNum(s_result->seq_id[HORIZONTAL]);
    index2 = GetSeqNum(s_result->seq_id[VERTICAL]);

    /* if either index is -1; the corresponding seq must have been deleted */
    if ((index1 == -1) || (index2 == -1)) {
	return;
    }
    seq1 = GetSeqSequence(index1);
    seq2 = GetSeqSequence(index2);
    seq1_len = GetSeqLength(index1);
    seq2_len = GetSeqLength(index2);
    

    


    Tcl_GetCommandInfo(interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;
    opts[0] = "-fg";
    opts[1] = "purple";
    opts[2] = NULL;
    index = CreateDrawEnviron(interp, raster, 2, opts);

    SetDrawEnviron(output->interp, raster, index);

    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

    start = (int)(data->win_len / 2);

    /* printing */
    for (i = 0; i < num_pts; i++) {

	x = data->p_array[i].x - start;
	y = data->p_array[i].y - start;

	for (j = 0; j < data->win_len; j++, x++, y++) {
	    
	    if ((x < 1) || (y < 1) || (x > seq1_len) || (y > seq2_len)) {
		continue;
	    }
	    score = score_matrix[char_lookup[seq1[x - 1]]]
	      [char_lookup[seq2[y - 1]]];
    
	    if (score >= min_score) {
		/* printf("j %d x %d y %d score %d \n", j, x, y, score); */
		coords[0] =  x;
		coords[1] =  (int)wy1 - y + wy0;
		RasterDrawPoints(raster, coords, 1);
	    }
	}
    }
    tk_RasterRefresh(raster);
}
#endif

void sip_rescan_matches(Tcl_Interp *interp,
			seq_result *s_result,
			int id,
			int min_score) 
{
    int i, j;
    double coords[2];
    int index;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int score;
    int start;
    int x, y;
    int index1, index2;
    double wx0, wy0, wx1, wy1;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
    element *e = s_result->e;
    int result_id = s_result->id;
    in_sim_spans *input = s_result->input;
    double invert_y;
    char cmd[1024];
    plot_data *result;
    Graph *graph;
    double m, c, pos;
    d_box bbox;

    result = find_plot_data(e, result_id);

    graph = Tcl_GetGraphFromObj(graph_obj);

    if (result->hidden) {
	return;
    }

    index1 = GetSeqNum(s_result->seq_id[HORIZONTAL]);
    index2 = GetSeqNum(s_result->seq_id[VERTICAL]);

    /* if either index is -1; the corresponding seq must have been deleted */
    if ((index1 == -1) || (index2 == -1)) {
	return;
    }
    seq1 = GetSeqSequence(index1);
    seq2 = GetSeqSequence(index2);
    seq1_len = GetSeqLength(index1);
    seq2_len = GetSeqLength(index2);
    
#if 0
    Tcl_GetCommandInfo(interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;
    opts[0] = "-fg";
    opts[1] = "purple";
    opts[2] = NULL;
    index = CreateDrawEnviron(interp, raster, 2, opts);

    SetDrawEnviron(output->interp, raster, index);

    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);
#endif

    start = (int)(input->win_length / 2);

    /* printing */
    for (i = 0; i < graph->p_arrays[0].n_pts; i++) {

	x = graph->p_arrays[0].p_array[i].x - start;
	y = graph->p_arrays[0].p_array[i].y - start;

	for (j = 0; j < input->win_length; j++, x++, y++) {
	    
	    if ((x < 1) || (y < 1) || (x > seq1_len) || (y > seq2_len)) {
		continue;
	    }
	    score = score_matrix[char_lookup[seq1[x - 1]]]
	      [char_lookup[seq2[y - 1]]];

	    if (score >= min_score) {
		/* printf("j %d x %d y %d score %d \n", j, x, y, score); */

		printf("x %d y %d %d %d\n", x, y, score, min_score);

		/*
		coords[0] =  x;
		coords[1] =  (int)wy1 - y + wy0;
		*/
		invert_y = graph->dim.y1 - y + graph->dim.y0;

		sprintf(cmd, "%s create line %d %f %d %f -fill purple -width %d -tag rescan", e->win, x, invert_y,
			x+1, invert_y+1, result->line_width);

		if (TCL_ERROR == Tcl_Eval(interp, cmd))
		    printf("ERROR %s\n", interp->result);

		/* RasterDrawPoints(raster, coords, 1); */
		
	    }
	}
    }
    /*    tk_RasterRefresh(raster); */

    printf("dim %f %f\n", graph->dim.y1, graph->dim.y0);

    sprintf(cmd, "%s create line %d %f %d %f -fill pink -width 10 -tag rescan", e->win, 500, graph->dim.y1 - 500+ graph->dim.y0, 500, graph->dim.y1 - 500+ graph->dim.y0);
    Tcl_Eval(interp, cmd);
 
    bbox = scale_box(e);

    /* use the data min and max so all graphs are to their own scale */
    if (e->orientation & HORIZONTAL) {
	m = (e->world->visible->y2 - e->world->visible->y1) /
	    (e->world->total->y2 - e->world->total->y1);
	
	c = e->world->visible->y2 - (m * e->world->total->y2);
	
	bbox.y1 = m * graph->dim.y0 + c;
	bbox.y2 = m * graph->dim.y1 + c;
    } 
    if (e->orientation & VERTICAL) {
	m = (e->world->visible->x2 - e->world->visible->x1) /
	    (e->world->total->x2 - e->world->total->x1);
	
	c = e->world->visible->x2 - (m * e->world->total->x2);
	
	/* note that graph->dim never changes orientation */
	if (!(e->orientation & HORIZONTAL)) {
	    /* vertical only */
	    bbox.x1 = m * graph->dim.y0 + c;
	    bbox.x2 = m * graph->dim.y1 + c;
	} else {
	    /* dot plot */
	    bbox.x1 = m * graph->dim.x0 + c;
	    bbox.x2 = m * graph->dim.x1 + c;
	}
    }

    {
    double x_origin, y_origin;
    double sf_x, sf_y;
    int p_x1, p_y1, p_x2, p_y2;
    int i;
    double w_x1, w_x2, w_y1, w_y2;

    /* initialise world coords */
    w_x1 = bbox.x1;
    w_x2 = bbox.x2;
    w_y1 = bbox.y1;
    w_y2 = bbox.y2;

    p_x1 = e->pixel->x;
    p_x2 = e->pixel->x + e->pixel->width;
    p_y1 = (double)e->pixel->y;
    p_y2 = (double)e->pixel->y + e->pixel->height;

    if (e->orientation & HORIZONTAL) {
	p_x1 = e->c->column[e->column_index]->pixel->x;
	p_x2 = e->c->column[e->column_index]->pixel->x + e->c->column[e->column_index]->pixel->width;
    }

    if (e->orientation & VERTICAL) {
	p_y1 = e->c->row[e->row_index]->pixel->y;
	p_y2 = e->c->row[e->row_index]->pixel->y + e->c->row[e->row_index]->pixel->height;
    }
     
    x_origin = calc_zoom_origin(w_x1, p_x1, w_x2, p_x2);
    sf_x = calc_zoom_sf((double)p_x1, w_x1, p_x2, w_x2);
    y_origin = calc_zoom_origin(w_y1, p_y1, w_y2, p_y2);
    sf_y = calc_zoom_sf(p_y1, w_y1, p_y2, w_y2);
    sprintf(cmd, "%s scale rescan %f %f %f %f \n", 
	    e->win, x_origin, y_origin, sf_x, sf_y);	

    printf("scale %s\n", cmd);
    Tcl_Eval(interp, cmd);
    }
 
    /*
    e->scrollregion_func(interp, e, e->world->total, 
			 e->c->column[e->column_index]->pixel,
			 e->c->row[e->row_index]->pixel);    
    */
}
