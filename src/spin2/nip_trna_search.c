#include <tcl.h>

#include "xalloc.h"
#include "seq_reg.h"
#include "seq_results.h"
#include "nip_results.h"
#include "text_output.h"
#include "nip_globals.h"
#include "trna_search.h"
#include "nip_trna_search.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "string.h"
#include "misc.h"
#include "seq_plot_funcs.h"
#include "dna_utils.h"
#include "tclCanvGraph.h"
#include "seq_element.h"

void trna_search_callback(int seq_num, void *obj, seq_reg_data *jdata);

void trna_search_shutdown(Tcl_Interp *interp,
			  seq_result *s_result, 
			  element *e,
			  int seq_num)
{
    in_trna_search *input = s_result->input;
    TrnaRes **text_data = s_result->text_data;
    seq_reg_key_name info;
    static char buf[80];

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    seq_deregister(seq_num, trna_search_callback, (seq_result *)s_result);
    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
    xfree(text_data);
    xfree(s_result->text_data);
    free(input->params);
    xfree(input->t);
    xfree(s_result->input);
    xfree(s_result);
}

void trna_search_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_trna_search *input = s_result->input;
    int result_id = s_result->id;
    char cmd[1024];
    plot_data *result;
    Tcl_Interp *interp;
    element *e = s_result->e;

    if (e) {
	result = find_plot_data(e, result_id);
	interp = e->c->interp;
    }

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "tRNA search");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "trna #%d", result_id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "trna: seq=%s", 
		GetSeqName(GetSeqNum(s_result->seq_id[HORIZONTAL])));
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
	    sprintf(cmd, "result_config %d %d %s %d %s", 
		    result_id, result->line_width, result->colour, e->id, 
		    e->win);
	    if (TCL_OK != Tcl_Eval(interp, cmd)){
		puts(interp->result);
	    }
	    break;
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
		trna_search_shutdown(interp, s_result, e, seq_num);
		remove_result_from_element(e, result_id);
		break;
	    }
	}
	break;
    case SEQ_PLOT:
	s_result->pr_func(s_result, (seq_reg_plot *)jdata);
	break;
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) 
	    {
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
		    pt.x = get_default_int(interp, tk_utils_defs,
					   w("ELEMENT.PLOT_WIDTH"));
		    pt.y = get_default_double(interp, tk_utils_defs,
					      w("ELEMENT.SINGLE.PLOT_HEIGHT"));
		    
		    jdata->info.result = (void *)&pt;
		    break;
		}
	    }
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: 
	{
	    trna_search_shutdown(interp, s_result, e, seq_num);
	    break;
	}
    }
}

void
trna_search_text_func(void *obj)
{
    seq_result *s_result = (seq_result *) obj;
    TrnaRes **text_data = s_result->text_data;
    int i;
    in_trna_search *input = s_result->input;
    TrnaSpec *t = input->t;
    int x, score;
    Graph *graph;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
    graph = Tcl_GetGraphFromObj(graph_obj);
    
    
    for (i = 0; i < graph->d_arrays[0].n_dlines; i++) {
	if ( text_data[i]->total_cb_score >= t->min_total_cb_score ) 
	    draw_trna(text_data[i]);
    }
    
    
    for (i = 0; i < graph->d_arrays[0].n_dlines; i++) {
	UpdateTextOutput();
	x = text_data[i]->aa_left+1;
	score = text_data[i]->total_bp_score;

	vmessage("Position %10d score %10d\n", x, score);
    }
}

int store_trna_search(int seq_num,
		      in_trna_search *input,
		      int start,
		      int end,
		      int strand,
		      TrnaRes **results, 
		      int num_elements,
		      TrnaSpec *t,
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

    if (NULL == (graph->d_arrays = (darray *)xmalloc(sizeof(darray))))
	return -1;

    if (NULL == (graph->d_arrays[0].d_array = (gd_line *)xmalloc(num_elements * sizeof(gd_line))))
	return -1;

    graph->d_arrays[0].n_dlines = num_elements;
    graph->d_arrays[0].type = G_LINES;
    graph->n_darrays = 1;

    graph->dim.x0 = start;
    graph->dim.x1 = end;
    graph->dim.y0 = 0;
    graph->dim.y1 = MAX_TRNA_BP;

    for (i = 0; i < num_elements; i++) {
	graph->d_arrays[0].d_array[i].x0 = results[i]->aa_left+1;
	graph->d_arrays[0].d_array[i].x1 = results[i]->aa_left+1;
	graph->d_arrays[0].d_array[i].y0 = 0;
	graph->d_arrays[0].d_array[i].y1 = results[i]->total_bp_score;
    }

    *graph_obj = Tcl_NewGraphObj(graph);
    /* these are copied over in Tcl_NewGraphObj so don't need them anymore */
    xfree(graph->d_arrays[0].d_array);
    xfree(graph->d_arrays);
    xfree(graph);

    id = get_reg_id();

    s_result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    s_result->seq_id[VERTICAL] = -1;

    s_result->id = id; 
    s_result->input = (void *)input; 
    s_result->text_data = results;
    s_result->type = SEQ_TYPE_TRNA;
    s_result->gr_type = SEQ_TYPE_GRAPH_PLOT;
    s_result->frame = 0;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */
    s_result->strand = strand;

    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = trna_search_callback;
    s_result->txt_func = trna_search_text_func;
   
    seq_register(seq_num, trna_search_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);
    return id;
}

int init_nip_trna_search_create(Tcl_Interp *interp,
				int seq_id,
				int start,
				int end,
				int strand,
				Tcl_Obj **graph_obj,
				int *id)
{
    in_trna_search *input;
    char *seq;
    int seq_len;
    int i;
    int seq_num;
    TrnaRes **results;
    int nmatch;
    TrnaSpec *t;
    int max_score = 0;
    Tcl_DString input_params;
    char strand_s[8];

    vfuncheader("trna search");
    set_char_set(DNA);

    if (NULL == (input = (in_trna_search *)xmalloc
		 (sizeof(in_trna_search))))
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

    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\nstrand %s\n",
		       GetSeqName(seq_num), start, end, strand_s);
    
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (NULL == (results = (TrnaRes **)xmalloc(MAX_TRNA * sizeof(TrnaRes *))))
	return -1;

    if (strand == TOP_S) {
	trna_search(seq, seq_len, start, end, &results, &nmatch, 
		    &max_score, &t);
    } else {
	complement_seq(seq, seq_len);
	trna_search(seq, seq_len, start, end, &results, &nmatch, 
		    &max_score, &t);
	complement_seq(seq, seq_len);
    }
    if (nmatch == 0) {
	verror(ERR_WARN, "trna search", "no matches found");
	
	for (i = 0; i < MAX_TRNA; i++) {
	    xfree(results[i]);
	}
	xfree(results);
	xfree(t);
	xfree(input->params);
	xfree(input);
	return -1;
    }

    input->t = t;

    *id = store_trna_search(seq_num, input, start, end, strand, results, 
			    nmatch, t, graph_obj);

    for (i = 0; i < nmatch; i++) {
	if ( results[i]->total_cb_score >= t->min_total_cb_score ) 
	    draw_trna(results[i]);
    }

    return 0;
}

int init_nip_trna_search_plot(Tcl_Interp *interp,
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
			      float tick_ht,
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

    configure->position = 0.5;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = tick_ht;
    configure->zoom = 1;
    configure->scroll = 0;

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
    result->amp_ruler = 0;
    sprintf(result->tags, "id%d", result_id);

    graph = Tcl_GetGraphFromObj(results);

    if (orientation == HORIZONTAL) {
	seq_id_h = seq_id;
	seq_id_v = -1;
    } else {
	seq_id_h = -1;
	seq_id_v = seq_id;
    }
    init_seq_element(interp, s_result, seq_id_h, seq_id_v, result, 
		     container_id, 
		     element_id, e_win, c_win, orientation, 
		     HORIZONTAL|VERTICAL,
		     graph, element_type);
    return 0;
}
