#include <tcl.h>
#include <string.h>
#include <math.h>

#include "xalloc.h"
#include "seq_reg.h"
#include "seq_results.h"
#include "nip_results.h"
#include "text_output.h"
#include "nip_globals.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "dna_utils.h"
#include "misc.h"
#include "nip_string_search.h"
#include "seq_plot_funcs.h"
#include "tclCanvGraph.h"
#include "seq_element.h"
#include "element_canvas.h"

void nip_string_search_callback(int seq_num, void *obj, seq_reg_data *jdata);

void nip_string_search_shutdown(Tcl_Interp *interp,
				seq_result *s_result,
				element *e,
				int seq_num)
{
    in_string_search *input = s_result->input;
    seq_reg_key_name info;
    static char buf[80];

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    seq_deregister(seq_num, nip_string_search_callback,(seq_result *)s_result);
		
    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);

    if (e->num_results > 0) {
	e->replot_func(e);
    }

    free(input->params);
    free(input->string);
    xfree(s_result->input);
    xfree(s_result);
}

void nip_string_search_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_string_search *input = s_result->input;
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
	sprintf(jdata->name.line, "string search");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "string #%d", result_id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "string: seq=%s", 
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
	  {
	      vfuncheader("input parameters");
	      vmessage("%s\n", input->params);
	      break;
	  }
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
		nip_string_search_shutdown(interp, s_result, e, seq_num);
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
	    nip_string_search_shutdown(interp, s_result, e, seq_num);
	break;
	}
    }
}

void
nip_string_search_text_func(void *obj)
{
    int i;
    int string_length;
    char *seq_match;
    seq_result *s_result = (seq_result *) obj;
    in_string_search *input = s_result->input;
    char *seq_name;
    char *sequence;
    int seq_num;
    int pos;
    double score;
    Graph *graph;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;

    graph = Tcl_GetGraphFromObj(graph_obj);
    
    seq_num = GetSeqNum(s_result->seq_id[HORIZONTAL]);
    seq_name = GetSeqName(seq_num);
    sequence = GetSeqSequence(seq_num);
    string_length = strlen(input->string);
    if (NULL == (seq_match = (char *)xcalloc((string_length + 1),
					     sizeof(char))))
	return;

    for (i = 0; i < graph->d_arrays[0].n_dlines; i++) {
	pos = graph->d_arrays[0].d_array[i].x0;
	score = graph->d_arrays[0].d_array[i].y1;

	vmessage("Position %d score %f", pos, score);

	strncpy(seq_match, &sequence[pos-1], string_length);
	iubc_list_alignment(input->string, seq_match, "string", seq_name, 1, 
			    pos, "");
    }
    xfree(seq_match);
}

int store_string_search(int seq_num,
			in_string_search *input,
			int start,
			int end,
			int *pos,
			int *score,
			int n_matches,
			int string_length,
			Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    int id;
    int i;
    Graph *graph;
    
    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;
    
    Tcl_InitGraph(&graph);

    if (NULL == (graph->d_arrays = (darray *)xmalloc(sizeof(darray))))
	return -1;

    if (NULL == (graph->d_arrays[0].d_array = (gd_line *)xmalloc(n_matches * sizeof(gd_line))))
	return -1;

    graph->d_arrays[0].n_dlines = n_matches;
    graph->d_arrays[0].type = G_LINES;
    graph->n_darrays = 1;

    graph->dim.x0 = start;
    graph->dim.x1 = end;
    graph->dim.y0 = 0;
    graph->dim.y1 = 100;

    for (i = 0; i < n_matches; i++) {
	graph->d_arrays[0].d_array[i].x0 = (double)pos[i] + start - 1;
	graph->d_arrays[0].d_array[i].x1 = (double)pos[i] + start - 1; 
	graph->d_arrays[0].d_array[i].y0 = 0;
	graph->d_arrays[0].d_array[i].y1 = (double)score[i] / string_length * 100;
#ifdef DEBUG
	printf("pos %f score %f\n", (double)pos[i], graph->d_arrays[0].d_array[i].y1);
#endif
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
    s_result->output = NULL;
    s_result->type = SEQ_TYPE_STRINGSEARCH;
    s_result->gr_type = SEQ_TYPE_GRAPH_PLOT;
    s_result->frame = 0;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */

    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = nip_string_search_callback;
    s_result->txt_func = nip_string_search_text_func;
    
    seq_register(seq_num, nip_string_search_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);
    return id;
}

int init_nip_string_search_create(int strand_sym,
				  float match,
				  char *string,
				  int use_iub_code,
				  int start,
				  int end,
				  int seq_id,
				  Tcl_Obj **graph_obj,
				  int *id)
{
    in_string_search *input;
    char *seq;
    int seq_len;
    int seq_num;
    int string_length;
    int max_matches, min_match;
    int n_matches;
    int *pos;
    int *score;
    Tcl_DString input_params;
    char strand[8];

    vfuncheader("string search");

    if (NULL == (input = (in_string_search *)xmalloc
		 (sizeof(in_string_search))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }

    seq_len = end - start + 1;

    max_matches = seq_len;
    string_length = strlen(string);

    if (NULL == (pos = (int *)xmalloc((max_matches + 1) * sizeof(int))))
	return -1;
    
    if (NULL == (score = (int *)xmalloc((max_matches + 1) * sizeof(int))))
	return -1;

    /* convert percentage mis-matches into min matches */
    min_match = (int) (ceil(string_length * match / 100)); 

    /* complement */
    if (strand_sym & BOTTOM_S) {
	complement_seq(string, string_length);
    }
    
#ifdef DEBUG
    printf("min_match %d match %f string %d\n", min_match, match, string_length);
#endif

    n_matches = iubc_inexact_match(&seq[start-1], seq_len, string, 
				   string_length, min_match, use_iub_code, 
				   pos, score, max_matches);

    if (n_matches <= 0) {
	vmessage("String search: no matches found\n");
	xfree(input);
	xfree(pos);
	xfree(score);
	return -1;
    }

    input->string = strdup(string);
    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    if (strand_sym & TOP_S && strand_sym & BOTTOM_S) 
	strcpy(strand, "both");
    else if (strand_sym & TOP_S) {
        strcpy(strand, "forward");
    } else {
        strcpy(strand, "reverse");
    }

    {
      char tmp[10];
      if (use_iub_code)
	strcpy(tmp, "iub");
      else
	strcpy(tmp, "literal");
      vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
			 "strand %s\nuse %s code\nminimum percent match %f\nstring %s\n",
			 GetSeqName(seq_num), start, end,
			 strand, tmp, match, string);
    }
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (-1 == (*id = store_string_search(seq_num, input, start, end, pos, 
					 score, n_matches, string_length,
					 graph_obj))){
	verror(ERR_FATAL,"string search", "error in saving matches\n");
	return -1;
    }

    xfree(pos);
    xfree(score);
    return 0;
}

int init_nip_string_search_plot(Tcl_Interp *interp,
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

    configure->position = 1.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = tick_ht;
    configure->zoom = 1;
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
		     container_id, element_id, e_win, c_win, orientation, 
		     HORIZONTAL|VERTICAL, graph, element_type);
    
    return 0;
}
