#include <tcl.h>
#include <string.h>
#include <limits.h>

#include "xalloc.h"
#include "seq_reg.h"
#include "seq_results.h"
#include "nip_results.h"
#include "text_output.h"
#include "nip_globals.h"
#include "tcl_utils.h"
#include "dna_utils.h"
#include "nip_base_comp.h"
#include "base_comp.h"
#include "seq_plot_funcs.h"
#include "tclCanvGraph.h"
#include "container.h"
#include "seq_element.h"

void plot_base_comp_callback(int seq_num, void *obj, seq_reg_data *jdata);

void plot_base_comp_shutdown(Tcl_Interp *interp,
			     seq_result *s_result,
			     element *e,
			     int seq_num)
{
    in_plot_base_comp *input = s_result->input;
    seq_reg_key_name info;
    static char buf[80];

#ifdef DEBUG  
    printf("plot_base_comp_shutdown %s %d\n", e->win, e->num_results);
#endif
    
    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    seq_deregister(seq_num, plot_base_comp_callback, (seq_result *)s_result);

    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);

    if (e->num_results > 0) {
	e->replot_func(e);
    }

    free(input->params);
    xfree(s_result->input);
    xfree(s_result);
}

void plot_base_comp_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_plot_base_comp *input = s_result->input;
    int result_id = s_result->id;
    char cmd[1024];
    element *e = s_result->e;
    plot_data *result;
    Tcl_Interp *interp;

    if (e) {
	result = find_plot_data(e, result_id);
	interp = e->c->interp;
    }

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Plot base composition");
	break;    

    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "base comp #%d", result_id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "base comp: seq=%s #%d", 
		GetSeqName(GetSeqNum(s_result->seq_id[HORIZONTAL])), result_id);
	break;
    case SEQ_GET_OPS:
	if (result->hidden) {
	    jdata->get_ops.ops = "Information\0List results\0"
		"PLACEHOLDER\0PLACEHOLDER\0Reveal\0SEPARATOR\0Remove\0";
	} else {
	    jdata->get_ops.ops = "Information\0List results\0Configure\0"
	       "Hide\0PLACEHOLDER\0SEPARATOR\0Remove\0";
	}
	break;
    case SEQ_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0: /* information */
	    vfuncheader("input parameters");
	    vmessage("%s\n", input->params);
	    break;
	case 1: /* results */
	    Tcl_Eval(interp, "SetBusy");
	    vfuncheader("results");
	    s_result->txt_func(s_result);
	    Tcl_Eval(interp, "ClearBusy");
	    break;
	case 2: /* configure */
	    {
		sprintf(cmd, "result_config %d %d %s %d %s", 
			result_id, result->line_width, result->colour, e->id, 
			e->win);
		if (TCL_OK != Tcl_Eval(interp, cmd)){
		    puts(interp->result);
		}
		break;
	    }
	case 3: /* hide all */
	    result->hidden = 1;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 4: /* reveal all */
	    result->hidden = 0;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 5: /* remove */ 
	    {
		plot_base_comp_shutdown(interp, s_result, e, seq_num);
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
		char *r_win;
		if (e) {
		   r_win = e->win; 
		   jdata->info.result = (void *)r_win;
		}
		break;
	    }
	case WIN_SIZE:
	    {
		static d_point pt;
		pt.x = get_default_int(interp, tk_utils_defs, 
					w("ELEMENT.PLOT_WIDTH"));
		pt.y = get_default_double(interp, tk_utils_defs,
					   w("ELEMENT.PLOT_HEIGHT"));

		jdata->info.result = (void *)&pt;
		break;
	    } /* WIN_SIZE */
	} /* switch SEQ_RESULT_INFO */
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: 
	{
	    plot_base_comp_shutdown(interp, s_result, e, seq_num);
	    break;
	} /* SEQ_DELETE, SEQ_QUIT */
    } /* switch */
}


void
plot_base_comp_text_func(void *obj)
{
    seq_result *s_result = (seq_result *) obj;
    int i;
    Graph *graph;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
    graph = Tcl_GetGraphFromObj(graph_obj);
 
    for (i = 0; i < graph->p_arrays[0].n_pts; i++) {
	UpdateTextOutput();
	vmessage("Position %10d score %10d\n", 
		 (int)graph->p_arrays[0].p_array[i].x, (int)graph->p_arrays[0].p_array[i].y);
    }
}

int store_base_comp(Tcl_Interp *interp,
		    int seq_num,
		    int window_length,
		    in_plot_base_comp *input,
		    double *match, 
		    int start,
		    int end,
		    int num_elements,
		    double min_match,
		    double max_match,
		    int strand,
		    Tcl_Obj **graph_obj)
{
    int i, j;
    seq_result *s_result;
    int id;
    Graph *graph;

    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;

    Tcl_InitGraph(&graph);

    if (NULL == (graph->p_arrays = (parray *)xmalloc(sizeof(parray))))
	 return -1;

    if (NULL == (graph->p_arrays[0].p_array = (g_pt *)xmalloc(sizeof(g_pt) * 
							      num_elements)))
	 return -1;

    id = get_reg_id();

    /* data assignment */
    if (strand & TOP_S) {
	for (i = 0; i < num_elements; i++) {
	    graph->p_arrays[0].p_array[i].x = i + start;
	    graph->p_arrays[0].p_array[i].y = match[i];
	}
    } else {
	for (i = 0, j = num_elements-1; i < num_elements; i++, j--) {
	    graph->p_arrays[0].p_array[i].x = i + start;
	    graph->p_arrays[0].p_array[i].y = match[j];
	}
    }
    graph->n_parrays = 1;
    graph->p_arrays[0].n_pts = num_elements;
    graph->p_arrays[0].type = G_LINE;
    graph->dim.x0 = start;
    graph->dim.x1 = end;
    graph->dim.y0 = min_match;
    graph->dim.y1 = max_match;

    *graph_obj = Tcl_NewGraphObj(graph);

    /* these are copied over in Tcl_NewGraphObj so don't need them anymore */
    xfree(graph->p_arrays[0].p_array);
    xfree(graph->p_arrays);
    xfree(graph);

    s_result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    s_result->seq_id[VERTICAL] = -1;

    s_result->id = id; 
    s_result->input = (void *)input; 
    s_result->type = SEQ_TYPE_BASECOMP;
    s_result->gr_type = SEQ_TYPE_GRAPH_PLOT;
    s_result->frame = 0;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */

    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = plot_base_comp_callback;
    s_result->txt_func =plot_base_comp_text_func;

    seq_register(seq_num, plot_base_comp_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);
    return id;
}

int init_nip_base_comp_create(Tcl_Interp *interp, 
			      int seq_id, 
			      int start, 
			      int end, 
			      int win_len, 
			      int a,
			      int c, 
			      int g, 
			      int t, 
			      int strand,
			      Tcl_Obj **graph_obj,
			      int *id)
{
    in_plot_base_comp *input;
    char *seq;
    int seq_len;
    double *match;
    int i;
    double score[5];
    int seq_num;
    int num_char;
    double min_match, max_match;
    Tcl_DString input_params;
    int irs;
    char strand_s[7];

    vfuncheader("plot base composition");
    set_char_set(DNA);

    if (NULL == (input = (in_plot_base_comp *)xmalloc
		 (sizeof(in_plot_base_comp))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }

    num_char = 5;

    for (i = 0; i < num_char; i++) {
	score[i] = 0.0;
    }
    if (a) 
	score[char_lookup['a']] = 1.0;
    if (c) 
	score[char_lookup['c']] = 1.0;
    if (g) 
	score[char_lookup['g']] = 1.0;
    if (t) 
	score[char_lookup['t']] = 1.0;

    if (NULL == (match = (double *)xmalloc((seq_len+2) * sizeof(double))))
	return -1;

    /* Rodger's version */

    /* NOTE that strand should only be TOP_S or BOTTOM_S, not both */
    if (strand == TOP_S) {
	irs = get_base_comp_res(seq, seq_len, win_len, start, end, score, 
				match, &min_match, &max_match);
    } else {
	complement_seq(seq, seq_len);
	irs = get_base_comp_res(seq, seq_len, win_len, start, end, score, 
				match, &min_match, &max_match);
	complement_seq(seq, seq_len);
	
    }
    if (irs == -1 || (min_match == 0 && max_match == 0)) {

	verror(ERR_WARN,"plot base composition", 
	       "error in calculating the base composition \n");
	xfree(input);
	xfree(match);
	return -1;
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
		       "window length %d A %d C %d G %d T %d strand %s\n",
		       GetSeqName(seq_num), start, end,
		       win_len, a, c, g, t, strand_s);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (-1 == (*id = store_base_comp(interp, seq_num, win_len, 
				     input, match, start, end,
				     end - start + 1, min_match, max_match,
				     strand, graph_obj))){
      verror(ERR_FATAL,"base composition", "error in saving matches\n");
      return -1;
    }

    xfree(match);
    return 0;
}

int init_nip_base_comp_plot(Tcl_Interp *interp,
			    int seq_id,
			    int result_id,
			    char *e_win,
			    char *c_win, 
			    Tcl_Obj *results,
			    int container_id,
			    int element_id,
			    char *element_type,
			    int line_width,
			    char *colour,
			    int orientation)
{
    seq_result *s_result;
    Graph *graph;
    configs *configure;
    plot_data *result;
    int seq_id_h, seq_id_v;

    s_result = seq_id_to_result(result_id);
    
   if (NULL == (result = (plot_data *)xmalloc(sizeof(plot_data))))
	return -1;

    if (NULL == (result->configure = (configs **)xmalloc(sizeof(configs*))))
	return -1;
    
    if (NULL == (configure = (configs *)xmalloc(sizeof(configs))))
	return -1;


    configure->position = -1.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = 0.0;
    configure->zoom = 2;
    configure->scroll = 1;

    result->configure[0] = configure;
    result->n_configure = 1;
    result->sf_m = 1.0;
    result->sf_c = 0.0;
    result->result_id = result_id;
    result->scale = SCALE_X | SCALE_Y;
    result->hidden = 0;
    result->line_width = line_width;
    result->colour = strdup(colour);
    result->len_ruler = 1;
    result->amp_ruler = 1;
    sprintf(result->tags, "id%d", result_id);

    if (orientation == HORIZONTAL) {
	seq_id_h = seq_id;
	seq_id_v = -1;
    } else {
	seq_id_h = -1;
	seq_id_v = seq_id;
    }

    graph = Tcl_GetGraphFromObj(results);
    
    init_seq_element(interp, s_result, seq_id_h, seq_id_v, result, 
		     container_id, 
		     element_id, e_win, c_win, orientation, 
		     HORIZONTAL|VERTICAL, graph, element_type);

    return 0;
}
