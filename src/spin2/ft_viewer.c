#include <stdio.h>
#include <string.h>
#include <tcl.h>
#include <tclXkeylist.h>

#include "misc.h"
#include "seq_reg.h"
#include "text_output.h"
#include "sequence_formats.h"
#include "xalloc.h"
#include "seq_results.h"
#include "canvas_box.h"
#include "nip_results.h"
#include "nip_globals.h"
#include "tcl_utils.h"
#include "ft_viewer.h"
#include "spin_globals.h"
#include "feature_selector.h"
#include "ruler_tick.h"
#include "tclCanvGraph.h"
#include "seq_element.h"
#include "element_canvas.h"

void
ft_viewer_callback(int seq_num, void *obj, seq_reg_data *jdata);

void ft_viewer_shutdown(Tcl_Interp *interp,
			seq_result *s_result,
			element *e,
			int seq_num)
{
    in_ft_viewer *input = s_result->input;
    seq_reg_key_name info;
    static char buf[80];

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    seq_deregister(seq_num, ft_viewer_callback, (seq_result *)s_result);

    if (e->num_results > 0) {
	e->replot_func(e);
    }
    free(input->params);
    xfree(s_result->input);
    xfree(s_result);

}

void ft_viewer_text_func(void *obj)
{

}

void circular_params(CanvasPtr *canvas,
		     int type,
		     int ft_offset,
		     int ft_height,
		     cir_s *circle)
{
    /* case for ruler */
    circle->diameter = (MIN(canvas->width, canvas->height)/2);
    circle->x1 = canvas->width/4;
    circle->y1 = canvas->height/4;

    if (type == FT_FORWARD) {
	circle->diameter = circle->diameter + (2*ft_offset);
	circle->x1 = circle->x1-ft_offset;
	circle->y1 = circle->y1-ft_offset;
    } else if (type == FT_REVERSE) {
	circle->diameter = circle->diameter - (2*ft_offset);
	circle->x1 = circle->x1+ft_offset;
	circle->y1 = circle->y1+ft_offset;
    }

}

void draw_circular_ruler(Tcl_Interp *interp,
			 ruler_s *ruler,
			 char *tags,
			 CanvasPtr *canvas,
			 double origin,
			 int disp_ticks) 
{
    char cmd[1024];
    cir_s circle;

    circular_params(canvas, FT_RULER, 0, 0, &circle);

    sprintf(cmd, "%s create oval %d %d %d %d -outline %s -width %d -tags {%s}",
	    ruler->window, circle.x1, circle.y1, circle.x1+circle.diameter, 
	    circle.y1+circle.diameter,
	    ruler->colour, ruler->line_width, tags);
    
    if (TCL_ERROR == Tcl_Eval(interp, cmd))
	verror(ERR_WARN, "draw_circular_ruler", "%s\n", interp->result);
    
    if (disp_ticks) {
	display_ruler_ticks_c(interp, canvas, ruler, ruler->start, ruler->end,
			      origin, circle);
    }
    
}

void draw_join_line(Tcl_Interp *interp,
		    feature_db_struct key,
		    char *window,
		    BasePos *loca1,
		    BasePos *loca2,
		    int offset,
		    int height)
{
    char buf[1024];
    float b1, t1, b2;

    b1 = (loca1->end_pos + loca1->start_pos)/2;
    b2 = (loca2->end_pos + loca2->start_pos)/2;
    t1 = (b2 + b1)/2;

    sprintf(buf, "%s create line %f %d %f %d %f %d -fill %s\n",
	    window, b1, offset, t1, height, b2, offset, key.fill);
    if (TCL_ERROR == Tcl_Eval(interp, buf)) {
	verror(ERR_WARN, "ft_viewer", "draw_join_line %s\n", 
	       interp->result);
    }
}

void draw_l_feature(Tcl_Interp *interp,
		    feature_db_struct key,
		    char *window,
		    char *tags,
		    int x1,
		    int y1,
		    int x2,
		    int y2)
{
    char buf[1024];
 
    if (strcmp(key.shape, "oval") == 0) {
	sprintf(buf, "%s create oval %d %d %d %d -fill %s -outline %s -width %d-tag {%s %s}\n", 
		window, x1, y1, x2, y2, key.fill, key.outline, key.width,
		key.type, tags);

    } else if (strcmp(key.shape, "diamond") == 0) {
	sprintf(buf, "%s create polygon %d %d %d %d %d %d %d %d -fill %s -outline %s -width %d -tag {%s %s}\n", 
		window, x1, (y1+y2)/2, (x1+x2)/2, y1, x2,(y1+y2)/2,
		(x1+x2)/2, y2, key.fill, key.outline, key.width, key.type, tags);
    } else {
	/* default to rectangle */
	sprintf(buf, "%s create rectangle %d %d %d %d -fill %s -outline %s -width %d -tag {%s %s}\n", 
		window, x1, y1, x2, y2, key.fill, key.outline, key.width, key.type, tags);
    }

    if (TCL_ERROR == Tcl_Eval(interp, buf)) {
	printf("DRAW %s\n", buf);

	verror(ERR_WARN, "ft_viewer", "draw_l_feature %s\n", 
	       interp->result);
    }
}

void draw_horizontal_features(Tcl_Interp *interp,
			      BasePos *loca,
			      element *e,
			      Graph *graph,
			      plot_data *result,
			      int strand,
			      char *tags,
			      double y_offset,
			      char *type_loca,
			      feature_db_struct key)
{
    static BasePos *loca1 = NULL;
    int y1, y2;

    y1 = e->pixel->height - (result->configure[0]->position + 
			       y_offset);
    y2 = e->pixel->height - (result->configure[0]->position +
			       result->configure[0]->height + y_offset);

    if (strcmp(type_loca, "n") == 0) {
	if (strand == TOP_S) {
	    draw_l_feature(interp, key, e->win, tags,
			   loca->start_pos, 
			   y1,
			   loca->end_pos, 
			   y2);
	}
    } else if (strcmp(type_loca, "c") == 0)  {
	if (strand == BOTTOM_S) {
	    draw_l_feature(interp, key, e->win, tags, loca->start_pos, 
			   y1,
			   loca->end_pos, 
			   y2);
	}
    } else if (strcmp(type_loca, "j") == 0) {
	if (strcmp(loca->type_range, "n") == 0) {
	    if (strand == TOP_S) {
		draw_l_feature(interp, key, e->win, tags,
			       loca->start_pos, y1, loca->end_pos, y2);
	    
		if (loca1) {
		    draw_join_line(interp, key, e->win, loca1, loca, 
				   result->configure[0]->position, 0);
		    loca1 = loca;
		} else {
		    loca1 = loca;
		}
	    }
	} else if (strcmp(type_loca, "c") == 0)  {
	    if (strand == BOTTOM_S) {
		draw_l_feature(interp, key, e->win, tags,
			       loca->start_pos, y1, loca->end_pos, y2);
		if (loca1) {
		    draw_join_line(interp, key, e->win, loca1, loca, 
				   result->configure[0]->position, 0);
		} else {
		    loca1 = loca;
		}
	    }
	} 
    }
}

void draw_vertical_features(Tcl_Interp *interp,
			    BasePos *loca,
			    element *e,
			    Graph *graph,
			    double y_offset_f,
			    double y_offset_r,
			    char *type_loca,
			    feature_db_struct key)
{
#if 0
    static BasePos *loca1 = NULL;
    CanvasPtr *canvas = data->canvas;

    if (strcmp(type_loca, "n") == 0) {
	draw_l_feature(interp, key, data->forward.win, 
		       data->ft_offset, 
		       loca->start_pos, 
		       data->ft_offset+data->ft_height,
		       loca->end_pos);

    } else if (strcmp(type_loca, "c") == 0)  {
	draw_l_feature(interp, key, data->reverse.win, 
		       data->ft_offset, 
		       loca->start_pos, 
		       data->ft_offset+data->ft_height,
		       loca->end_pos);
    } else if (strcmp(type_loca, "j") == 0)  {
	if (strcmp(loca->type_range, "n") == 0) {
	    draw_l_feature(interp, key, data->forward.win, 
			   data->ft_offset+y_offset_f, loca->start_pos, 
			   data->ft_offset+data->ft_height+y_offset_f,
			   loca->end_pos);
	    if (loca1) {
		draw_join_line(interp, key, data->forward.win, loca1, loca, 
			       data->ft_offset, 0);
		loca1 = loca;
	    } else {
		loca1 = loca;
	    }
	} else if (strcmp(type_loca, "c") == 0)  {
	    draw_l_feature(interp, key, data->reverse.win, 
			 data->ft_offset, loca->start_pos, 
			 data->ft_offset+data->ft_height, loca->end_pos);
	    if (loca1) {
		draw_join_line(interp, key, data->reverse.win, loca1, loca, 
			       data->ft_offset, 0);
	    } else {
		loca1 = loca;
	    }
	} 
    }
#endif
}

void draw_c_feature(Tcl_Interp *interp,
		    feature_db_struct key,
		    char *win,
		    char *tags,
		    cir_s circle,
		    double origin,
		    int x1,
		    int x2,
		    CanvasPtr *canvas,
		    int seq_len)
{
    char buf[1024];
    double start;
    double extent;

    /* * -360 for clockwise rotation */
    start = ((double)x1 / seq_len * -360) + origin;
    extent = (double)(x2-x1) / seq_len * -360;

    /* awful HACK to get round bug in arc which won't draw a complete circle */
    if ((x2-x1+1) == seq_len) {
	sprintf(buf, "%s create oval %d %d %d %d -outline %s -width 20 -tag {%s %s}\n",
		win, circle.x1, circle.y1, 
		circle.x1+circle.diameter, 
		circle.y1+circle.diameter, key.outline, key.type, tags);
    } else {
	sprintf(buf, "%s create arc %d %d %d %d -style arc -extent %f -start %f -outline %s -width 20 -tag {%s %s}\n",
		win, circle.x1, circle.y1, 
		circle.x1+circle.diameter, 
		circle.y1+circle.diameter, extent, start, 
		key.outline, key.type, tags);
    }

    if (TCL_ERROR == Tcl_Eval(interp, buf)) {
	verror(ERR_WARN, "ft_viewer", "draw_c_feature %s\n", 
	       interp->result);
    }
}

int feature_pos(BasePos *loca,
		int x0,
		int x1,
		int *start_pos,
		int *end_pos)
{
    *start_pos = loca->start_pos;
    *end_pos = loca->end_pos;

    if (loca->end_pos < x0) {
	return 0;
    }

    if (loca->start_pos > x1) {
	return 0;
    }

    if (loca->start_pos < x0)
	*start_pos = x0;

    if (loca->end_pos > x1)
	*end_pos = x1;

    return 1;
}


void draw_circular_features(Tcl_Interp *interp,
			    BasePos *loca,
			    element *e,
			    int x0,
			    int x1,
			    plot_data *result,
			    int strand,
			    char *tags,
			    double y_offset,
			    cir_s forw,
			    cir_s rev,
			    char *type_loca,
			    feature_db_struct key)
{
    static BasePos *loca1 = NULL;
    int seq_len = x1 - x0 + 1;
    int start_pos, end_pos;
    feats *extra_data = result->extra_data;

    if (!feature_pos(loca, x0, x1, &start_pos, &end_pos))
	return;
    start_pos -= x0 - 1;
    end_pos -= x0 - 1;

    if (strcmp(type_loca, "n") == 0) {
	if (strand & TOP_S) {
	    draw_c_feature(interp, key, e->win, tags, forw, 
			   extra_data->origin, start_pos, 
			   end_pos, e->pixel, seq_len);
	}
    } else if (strcmp(type_loca, "c") == 0)  {
	if (strand &  BOTTOM_S) {
	    draw_c_feature(interp, key, e->win, tags, rev,
			   extra_data->origin, start_pos, end_pos, 
			   e->pixel, seq_len);
	}
    } else if (strcmp(type_loca, "j") == 0)  {
	if (strcmp(loca->type_range, "n") == 0) {
	    if (strand & TOP_S) {
		draw_c_feature(interp, key, e->win, tags, forw,
			       extra_data->origin, start_pos, end_pos, 
			       e->pixel, seq_len);
		if (loca1) {
		    draw_join_line(interp, key, e->win, loca1, loca, 
				   result->configure[0]->position, 0);
		    loca1 = loca;
		} else {
		    loca1 = loca;
		}
	    }
	} else if (strcmp(type_loca, "c") == 0)  {
	    if (strand & BOTTOM_S) {
		draw_c_feature(interp, key, e->win, tags, rev,
			       extra_data->origin, start_pos, end_pos, 
			       e->pixel, seq_len);
		if (loca1) {
		    draw_join_line(interp, key, e->win, loca1, loca, 
			       result->configure[0]->position, 0);
		} else {
		    loca1 = loca;
		}
	    }
	} 
    }
}

/* convert a tcl list of selected feature keys into a feature_db_struct */
feature_db_struct* ft_list_to_key(Tcl_Interp *interp,
				  char *items,
				  int *num_keys)
{
    char **array;
    int num_items;
    feature_db_struct *keys;

    if (Tcl_SplitList(interp, items, &num_items, &array) != TCL_OK) {
	return NULL;
    }
    keys = klist_feature(interp, array, num_items);
    *num_keys = num_items;
    return keys;
}

void ft_viewer_plot_func(void *obj, seq_reg_plot *plot)
{
    seq_result *s_result = (seq_result *) obj;
    int i, j, ii;
    Featcds **key_index = NULL;
    int seq_num;
    BasePos *loca;
    char cmd[1024];
    char buf[1024];
    feature_db_struct *keys;
    int num_keys;
    Tcl_Interp *interp;
    cir_s forw, rev;
    double y_offset;
    element *e = s_result->e;
    Graph *graph;
    d_box bbox;
    plot_data *result;
    feats *extra_data;
    in_ft_viewer *input = s_result->input;
    double pos;
    int win_height;

    interp = e->c->interp;

    graph = Tcl_GetGraphFromObj(s_result->data);
    result = find_plot_data(e, s_result->id);

    extra_data = result->extra_data;
    y_offset = 0.0;

    key_index = (Featcds **)xmalloc(number_keys * sizeof(Featcds*));
    if (NULL == key_index)
        goto fail;

    if (e->orientation == CIRCLE) {
	seq_num = GetSeqNum(s_result->seq_id[HORIZONTAL]);
    } else {
	seq_num = GetSeqNum(s_result->seq_id[e->orientation]);
    }
    key_index = GetSeqKeyIndex(seq_num);

    sprintf(buf, "get_default_feats %s %d",  
	    get_default_string(interp, spin_defs, "FT.NAME"), s_result->id);
    
    if (TCL_ERROR == Tcl_Eval(interp, buf)) {
	verror(ERR_WARN, "ft_viewer", "plot_ft_viewer %s\n", 
	       interp->result);
    }

    if (e->orientation == CIRCLE) {
	if (e->ruler == NULL) {
	    e->ruler = ruler_struct(interp, tk_utils_defs, "CONTAINER", 0);
	    e->ruler->start = input->start;
	    e->ruler->end = input->end;
	    e->ruler->window = strdup(e->win);
	}
	circular_params(e->pixel, FT_FORWARD, 
			result->configure[0]->position, 
			result->configure[0]->height,
			&forw);
	circular_params(e->pixel, FT_REVERSE, 
			result->configure[0]->position, 
			result->configure[0]->height,
			&rev);
    }

    keys = ft_list_to_key(interp, interp->result, &num_keys);

    sprintf(cmd, "%s delete all", e->win);
    Tcl_Eval(interp, cmd);

    for (ii = 0; ii < num_keys; ii++) {
	i = keys[ii].index;

	for (j = 1; j <= key_index[i]->id; j++) {
	    
            /* printf("    %s instance %d %d\n", feat_key[i], j, i); */
	    
            loca = key_index[i][j].loca;
	    if (extra_data->display == FT_ALL) {
		y_offset += (result->configure[0]->position * 
			     result->configure[0]->height);
	    }

            /* Iterate around the location for this feature */
            while (loca) {
		/*
                printf("        loca=%d..%d type='%s'\n",
                       loca->start_pos,
                       loca->end_pos,
                       loca->type_range);
		*/		
		if (e->orientation == CIRCLE) {
		    draw_circular_features(interp, loca, e, input->start,
					   input->end, result,
					   s_result->strand, result->tags, 
					   y_offset, forw, rev, 
					   key_index[i][j].type_loca, 
					   keys[ii]);
	    
		} else if (e->orientation == VERTICAL) {
		    draw_vertical_features(interp, loca, e, graph, y_offset,
					   y_offset,
					   key_index[i][j].type_loca, 
					   keys[ii]);
		} else {
		    draw_horizontal_features(interp, loca, e, graph, result,
					     s_result->strand, result->tags, 
					     y_offset,
					     key_index[i][j].type_loca, 
					     keys[ii]);
		}
		/* raw_feature(interp, loca); */
                loca = loca->next;
            }
            /* printf("        type_loca='%s'\n", key_index[i][j].type_loca); */
	}
    }

    if (e->orientation == CIRCLE) {
	draw_circular_ruler(interp, e->ruler, result->tags, e->pixel, extra_data->origin, 1);
    }
    
    bbox = scale_box(e);

    e->scale_func(interp, e, -1, &bbox, e->pixel); 

    e->scrollregion_func(interp, e, e->world->total, 
			 e->c->column[e->column_index]->pixel,
			 e->c->row[e->row_index]->pixel);
    
    if (e->orientation & HORIZONTAL) {
	win_height = e->element_height(interp, e->win);
    } else {
	win_height = e->element_width(interp, e->win);
    }
     for (i = 0; i < result->n_configure; i++) {
	if (result->configure[i]->position != -1.0) {
	    if (e->orientation & HORIZONTAL) {
		pos = win_height * result->configure[i]->position;
		if (result->configure[i]->y_direction == '+') {
		    pos = pos - win_height;
		}
#ifdef DEBUG
		printf("POS %f %f %f\n", pos, graph->dim.y1, graph->dim.y0);
#endif
		canvas_move(interp, e, s_result->id, 0, pos);
	    }
	    if (e->orientation & VERTICAL) {
		printf("TODO\n");
	    }
	}
    }
   
    

    return;
 fail:
    printf("FAILED\n");
    if (key_index) {
        for (i = 0; i < number_keys; i++) {
            if (key_index[i])
                xfree(key_index[i]);
        }
        xfree(key_index);
    }
}

void
ft_viewer_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    element *e = s_result->e;
    Tcl_Interp *interp;
    int result_id = s_result->id;
    plot_data *result;

    if (e) {
	result = find_plot_data(e, result_id);
	interp = e->c->interp;
    }

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Feature table");
	break;    

    case SEQ_GET_BRIEF:
	/*
	sprintf(jdata->name.line, "ft_viewer: seq=%s", GetSeqName(GetSeqNum(s_result->seq_id[HORIZONTAL])));
	*/
	Tcl_VarEval(interp, "highlight_ft_item ", e->win, NULL);
	sprintf(jdata->name.line, "%s", interp->result);

	break;

    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "ft_viewer #%d", result_id);
	break;

    case SEQ_GET_OPS:
	jdata->get_ops.ops = "Display mode\0Remove\0";
	break;
    case SEQ_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0: /* Display mode */
	    {
		printf("display mode\n");
	    }
	    break;
	case 1: /* Remove */
	    {
		ft_viewer_shutdown(interp, s_result, e, seq_num);
		remove_result_from_element(e, result_id);
		break;
	    }
	}
	break;
    case SEQ_PLOT:
	s_result->pr_func(s_result, (seq_reg_plot *)jdata);
	break;
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case RESULT:
	    jdata->info.result = (void *)s_result;
	    break;
	case WIN_NAME:
	    {
		char *r_win = e->win;
		jdata->info.result = (void *)r_win;
		break;
	    }
	case WIN_SIZE:
	    {
		static d_point pt;
		pt.x = get_default_int(interp, spin_defs, 
					w("FT.PLOT_WIDTH"));
		pt.y = get_default_double(interp, spin_defs,
					   w("FT.SINGLE.PLOT_HEIGHT"));

		jdata->info.result = (void *)&pt;
		break;
	    } /* WIN_SIZE */
	} /* switch SEQ_RESULT_INFO */
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: 
	{
	    ft_viewer_shutdown(interp, s_result, e, seq_num);
	    break;
	} /* SEQ_DELETE */
    case SEQ_GENERIC:
	break; /* SEQ_GENERIC */
    } /* switch */
}


int ft_viewer_reg(Tcl_Interp *interp,
		  int seq_id,
		  char *frame,
		  char *forward,
		  char *reverse,
		  int start,
		  int end,
		  int plot_type,
		  int ft_offset,
		  int ft_height,
		  ruler_s *ruler,
		  cursor_s cursor,
		  double origin)
{
    seq_result *result;
    out_canvas *output;
    char *seq;
    int seq_len;
    int seq_num;
    int sequence_type;
    seq_cursor_notify cn;
    int line_width;
    int id;
    ft_viewer_res *data;

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

   if (NULL == (data = (ft_viewer_res *)xmalloc(sizeof(ft_viewer_res))))
	return -1;

    if (NULL == (output = (out_canvas *)xmalloc(sizeof(out_canvas))))
	return TCL_OK;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    sequence_type = GetSeqStructure(seq_num);

    seq = &seq[start-1];
    seq_len = end - start + 1;

    result->data = data;

    if (plot_type == VERTICAL) {
	result->seq_id[HORIZONTAL] = -1;
	result->seq_id[VERTICAL] = GetSeqId(seq_num);
    } else {
	/* if horizontal or circular */
	result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
	result->seq_id[VERTICAL] = -1;
    }

    id = get_reg_id();

    result->id = id; 

    /* result->input = (void *)input; */
    result->output = (void *)output;

    result->pr_func = ft_viewer_plot_func;
    result->op_func = ft_viewer_callback;
    result->txt_func = ft_viewer_text_func;

    /* data structure */
    data->frame = strdup(frame);

    data->forward.win = strdup(forward);
    data->forward.display = FT_SINGLE;
    data->forward.offset = 1.0;
    data->forward.min = 0;
    data->forward.max = ft_height + ft_offset;

    data->reverse.win = strdup(reverse);
    data->reverse.display = FT_SINGLE;
    data->reverse.offset = 1.0;
    data->reverse.max = ft_height + ft_offset;
    data->reverse.min = 0;

    data->ft_offset = ft_offset;
    data->ft_height = ft_height;
    data->plot_type = plot_type;
    
    data->ruler = ruler;
    data->sequence_len = GetSeqLength(seq_num);
    data->cursor = cursor;
    data->origin = origin;

    /* create list of windows in the feature table display */
    if (NULL == (data->win_list = (win **)xmalloc(MAX_NUM_WINS * sizeof(win*))))
	return -1;

    data->num_wins = 0;
    
    addWindow(data->win_list, &data->num_wins, data->forward.win, 'b', id);
    addWindow(data->win_list, &data->num_wins, data->ruler->window, 'x', id);
    addWindow(data->win_list, &data->num_wins, data->reverse.win, 'b', id);
    

   if (NULL == (data->canvas = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    if (NULL == (data->world= (WorldPtr *)xmalloc(sizeof(WorldPtr))))
	return -1;

    if (NULL == (data->world->visible = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    if (NULL == (data->world->total = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    if (plot_type == CIRCLE) {
	initCanvas(interp, data->canvas, data->ruler->window);
    } else {
	initCanvas(interp, data->canvas, data->forward.win);
    }
    createZoom(&data->zoom);

    data->sequence_type = sequence_type;

    output->interp = interp;
    line_width = get_default_int(interp, nip_defs, w("NIP.CURSOR.LINE_WIDTH"));

    output->cursor = create_cursor(seq_num, 0, NULL, line_width, 1, 
				   HORIZONTAL, 1);
    output->cursor_visible = 0;

    /* move cursor to start position if this is our own cursor */
    if (output->cursor->refs == 1)
	output->cursor->abspos = start;
    
    seq_register(seq_num, ft_viewer_callback, (void *)result, SEQ_PLOT_PERM, id);

    ft_viewer_plot_func(result, NULL);

    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = output->cursor;
    cn.cursor->job = CURSOR_MOVE;
    cn.show_all = 0;
    seq_notify(seq_num, (seq_reg_data *)&cn); 
    return id;
    
}

int store_ft_viewer(Tcl_Interp *interp,
		    int seq_num,
		    in_ft_viewer *input,
		    int start,
		    int end,
		    int strand,
		    int orientation,
		    Tcl_Obj **graph_obj)
{
   seq_result *s_result;
   int id, i;
   Graph *graph;

   if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
       return -1;

    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;

    Tcl_InitGraph(&graph);

    graph->dim.x0 = start;
    graph->dim.x1 = end;

    *graph_obj = Tcl_NewGraphObj(graph);

    id = get_reg_id();

    /* initialise seq_id */
    for (i = 0; i <= CIRCLE; i++) {
	s_result->seq_id[i] = -1;
    }

    if (orientation == CIRCLE) {
	s_result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    } else {
	s_result->seq_id[orientation] = GetSeqId(seq_num);
    }

    s_result->id = id; 
    s_result->input = (void *)input; 
    s_result->type = SEQ_TYPE_FT_VIEWER;
    s_result->frame = 0;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj;
    s_result->e = NULL; /* initialise */
    s_result->strand = strand;
    s_result->gr_type = SEQ_TYPE_FT_VIEWER;

    s_result->pr_func = ft_viewer_plot_func;
    s_result->op_func = ft_viewer_callback;
    s_result->txt_func = ft_viewer_text_func;

    seq_register(seq_num, ft_viewer_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);
    return id;
}

int init_ft_viewer_create(Tcl_Interp *interp,
			  int seq_id,
			  int start,
			  int end,
			  int strand,
			  int orientation,
			  Tcl_Obj **graph_obj,
			  int *id)
{
    char strand_s[7];
    Tcl_DString input_params;
    in_ft_viewer *input;
    char *seq;
    int seq_len;
    int seq_num;

    vfuncheader("feature table viewer");

    if (NULL == (input = (in_ft_viewer *)xmalloc(sizeof(in_ft_viewer))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);

    if (strand & TOP_S) {
        strcpy(strand_s, "top");
    } else if (strand & BOTTOM_S) {
        strcpy(strand_s, "bottom");
    } else {
        strcpy(strand_s, "both");
    }

    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "strand %s\n",
		       GetSeqName(seq_num), start, end, strand_s);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    input->start = start;
    input->end = end;

    if (-1 == (*id = store_ft_viewer(interp, seq_num, input, start, end,
				     strand, orientation, graph_obj))){
	verror(ERR_FATAL,"base composition", "error in saving matches\n");
	return -1;
    }
    return 0;

}

int init_ft_viewer_plot(Tcl_Interp *interp,
			int seq_id,
			int result_id,
			char *e_win,
			char *c_win, 
			Tcl_Obj *results,
			int container_id,
			int element_id,
			char *colour,
			int orientation,
			int tick_ht,
			float offset,
			char *element_type,
			double origin,
			int display_mode)
{
    seq_result *s_result;
    Graph *graph;
    configs *configure;
    plot_data *result;
    feats *extra_data;
    in_ft_viewer *input;

    s_result = seq_id_to_result(result_id);
    
   if (NULL == (result = (plot_data *)xmalloc(sizeof(plot_data))))
	return -1;

    if (NULL == (result->configure = (configs **)xmalloc(sizeof(configs*))))
	return -1;
    
    if (NULL == (configure = (configs *)xmalloc(sizeof(configs))))
	return -1;

    if (NULL == (extra_data = (feats *)xmalloc(sizeof(feats))))
	return -1;
    
    configure->position = offset;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = tick_ht;
    configure->zoom = 2;
    configure->scroll = 1;

    result->configure[0] = configure;
    result->n_configure = 1;
    result->sf_m = 1.0;
    result->sf_c = 0.0;
    result->result_id = result_id;
    if (orientation == CIRCLE) {
	result->scale = SCALE_X|SCALE_Y;
    } else {
	result->scale = SCALE_X;
    }
    result->hidden = 0;
    result->colour = strdup(colour);
    result->len_ruler = 1;
    result->amp_ruler = 0;
    
    extra_data->origin = origin;
    extra_data->display = display_mode;
    
    result->extra_data = (feats *)extra_data;

    input = s_result->input;

    sprintf(result->tags, "id%d", result_id);

    graph = Tcl_GetGraphFromObj(results);

    if (orientation == CIRCLE) {
	graph->dim.x0 = 0;
	graph->dim.x1 = canvas_width(interp, e_win);
	graph->dim.y0 = 0;
	graph->dim.y1 = canvas_height(interp, e_win);
    } else {
	graph->dim.y0 = 1;
	graph->dim.y1 = tick_ht;
    }

    init_seq_element(interp, s_result, seq_id, -1, result, container_id, 
		     element_id, e_win, c_win, orientation, HORIZONTAL, graph, 
		     element_type);

    return 0;

}
