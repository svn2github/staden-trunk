#include <tcl.h>
#include <string.h>
#include <float.h>

#include "misc.h"
#include "xalloc.h"
#include "seq_reg.h"
#include "nip_results.h"
#include "text_output.h"
#include "codon_content.h"
#include "nip_globals.h"
#include "tcl_utils.h"
#include "array_arith.h"
#include "seq_results.h"
#include "nip_gene_search.h"
#include "seq_plot_funcs.h"
#include "tclCanvGraph.h"
#include "seq_element.h"
#include "dna_utils.h"

void plot_gene_search_callback(int seq_num, void *obj, seq_reg_data *jdata);

void plot_gene_search_shutdown(Tcl_Interp *interp,
			       seq_result *s_result,
			       element *e,
			       int seq_num)
{
    in_plot_gene_search *input = s_result->input;
    seq_reg_key_name info;
    static char buf[80];

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    /* need to deregister sequence */
    seq_deregister(seq_num, plot_gene_search_callback, (seq_result *)s_result);

    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);

    if (e->num_results > 0) {
	e->replot_func(e);
    }

    free(input->params);
    xfree(s_result->input);
    xfree(s_result);
  
}

void
plot_gene_search_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_plot_gene_search *input = s_result->input;
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
	sprintf(jdata->name.line, "Plot gene search");
	break;
	
    case SEQ_KEY_NAME:
	if (s_result->frame) {
	    sprintf(jdata->name.line, "gene f%d #%d", s_result->frame, s_result->id);
	} else {
	    /* single plot encompassing all frames */
	    sprintf(jdata->name.line, "gene #%d", result_id);
	}
	break;

    case SEQ_GET_BRIEF:
	if (s_result->frame) {
	    sprintf(jdata->name.line, "gene: seq=%s frame=%d", 
		    GetSeqName(GetSeqNum(s_result->seq_id[HORIZONTAL])), s_result->frame);
	} else {
	    sprintf(jdata->name.line, "gene: seq=%s", 
		    GetSeqName(GetSeqNum(s_result->seq_id[HORIZONTAL])));

	}
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
	    vmessage("%s \n", input->method);
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
		plot_gene_search_shutdown(interp, s_result, e, seq_num);
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
	    }
	}
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: {
	plot_gene_search_shutdown(interp, s_result, e, seq_num);
	break;
	}
    }
}

void plot_gene_search_text_func(void *obj)
{
    seq_result *s_result = (seq_result *) obj;
    Tcl_Obj *graph_obj = s_result->data;
    int i;
    Graph *graph;

    graph = Tcl_GetGraphFromObj(graph_obj);

    for (i = 0; i < graph->p_arrays[0].n_pts; i++) {
	UpdateTextOutput();
	vmessage("Position %10f score %.5g \n", graph->p_arrays[0].p_array[i].x,
		 graph->p_arrays[0].p_array[i].y);
    }

}

int DoCodonPref(char *seq,
		int seq_length,
		char *codon_table,
		int window_length,
		int start,
		int end,
		int option,
		CodRes **matches)
{
    CodRes *results;
    int res;
    int num_results;
    double codon_usage_table[4][4][4];

    num_results = end - start + 1;
    num_results = 1 + num_results/3;

    if ((results = init_CodRes(num_results)) == NULL) return -2;

    results->num_results = num_results;
    results->window_length = window_length;
    results->user_start = start;
    results->user_end = end;

    res = init_codon_pref (codon_table, codon_usage_table, option);
    if ( res ) {
	free_CodRes ( results );
	return -1;
    }

    res = do_codon_pref ( seq, seq_length, codon_usage_table, results);

    if ( res ) {
	free_CodRes ( results );
	return -1;
    }

    *matches = results;
    return 0;
}

int DoAuthorTest(char *seq,
		 int seq_length,
		 char *codon_table,
		 double percent_error,
		 int start,
		 int end,
		 CodRes **matches)
{
    CodRes *results;
    int res;
    int num_results;
    double codon_usage_table[4][4][4];
    int window_length;

    num_results = end - start + 1;
    num_results = 1 + num_results/3;

    if ( ( results = init_CodRes ( num_results ) ) == NULL ) return -2;

    results->num_results = num_results;
    results->error = percent_error;
    results->user_start = start;
    results->user_end = end;

    res = init_author_test(codon_table, seq, seq_length, codon_usage_table,
			   percent_error, &window_length);
    if ( res ) {
	free_CodRes ( results );
	return -1;
    }
    results->window_length = window_length;

    res = do_author_test(seq, seq_length, codon_usage_table, results);

    if ( res ) {
	free_CodRes ( results );
	return -1;
    }
   
    *matches = results;
    return 0;
}

int DoPosBaseBias(char *seq,
		  int seq_length,
		  int window_length,
		  int start,
		  int end,
		  CodRes1 **matches)
{
    CodRes1 *results;
    int res;
    int num_results;

    num_results = end - start + 1;
    num_results = 1 + num_results/3;

    if ( ( results = init_CodRes1 ( num_results ) ) == NULL ) return -2;

    results->num_results = num_results;
    results->window_length = window_length;
    results->user_start = start;
    results->user_end = end;

    res = do_pos_base_bias(seq, seq_length, results);

    if ( res ) {
	free_CodRes1 ( results );
	return -1;
    }

    *matches = results;
    return 0;
}

int store_gene_search(int seq_num,
		      int start,
		      int end,
		      int frame,
		      in_plot_gene_search *input,
		      double *score, 
		      char *top,
		      int num_elements,
		      double min,
		      double max,
		      int type,
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

    if (NULL == (graph->p_arrays = (parray *)xmalloc(2 * sizeof(parray))))
	 return -1;

    if (NULL == (graph->p_arrays[0].p_array = (g_pt *)xmalloc(sizeof(g_pt) * 
							      num_elements)))
	 return -1;
    
    if (NULL == (graph->p_arrays[1].p_array = (g_pt *)xmalloc(sizeof(g_pt) * 
							      num_elements)))
	 return -1;
 
    id = get_reg_id();

    graph->p_arrays[1].n_pts = 0;

    /* data assignment */
    if (strand == TOP_S) {
	for (i = 0; i < num_elements; i++) {
	    graph->p_arrays[0].p_array[i].x = (start-1) + (i * 3) + frame;
	    graph->p_arrays[0].p_array[i].y = score[i];
	    
	    if (top && frame == top[i]) {
		graph->p_arrays[1].p_array[graph->p_arrays[1].n_pts].x = 
		    graph->p_arrays[0].p_array[i].x;
		graph->p_arrays[1].p_array[graph->p_arrays[1].n_pts].y = (max+min)/2;
		graph->p_arrays[1].n_pts++;
	    }
	}
    } else {
	for (i = 0, j = num_elements-1; i < num_elements; i++, j--) {
	    graph->p_arrays[0].p_array[i].x = (start-1) + (i * 3) + frame;
	    graph->p_arrays[0].p_array[i].y = score[j];
	    
	    if (top && frame == top[j]) {
		graph->p_arrays[1].p_array[graph->p_arrays[1].n_pts].x = 
		    graph->p_arrays[0].p_array[i].x;
		graph->p_arrays[1].p_array[graph->p_arrays[1].n_pts].y = (max+min)/2;
		graph->p_arrays[1].n_pts++;
	    }
	}
    }

    graph->p_arrays[0].n_pts = num_elements;
    graph->p_arrays[0].type = G_LINE;
    graph->p_arrays[1].type = G_DOT;
    graph->n_parrays = 2;
    graph->dim.x0 = start;
    graph->dim.x1 = end;
    graph->dim.y0 = min;
    graph->dim.y1 = max;

    *graph_obj = Tcl_NewGraphObj(graph);
    /* these are copied over in Tcl_NewGraphObj so don't need them anymore */
    xfree(graph->p_arrays[0].p_array);
    xfree(graph->p_arrays[1].p_array);
    xfree(graph->p_arrays);
    xfree(graph);

    s_result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    s_result->seq_id[VERTICAL] = -1;

    s_result->id = id; 
    s_result->input = (void *)input; 
    s_result->type = type;
    s_result->gr_type = SEQ_TYPE_GRAPH_PLOT;
    s_result->frame = frame;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */

    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = plot_gene_search_callback;
    s_result->txt_func =plot_gene_search_text_func;
    
    seq_register(seq_num, plot_gene_search_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);
    return id;
}

int init_nip_codon_pref_create(Tcl_Interp *interp,
			       int seq_id,
			       int start,
			       int end,
			       char *codon_table,
			       int win_len,
			       int option,
			       int strand,
			       Tcl_Obj **graph_obj,
			       int *id)
{
    char *seq;
    int seq_len, seq_num;
    in_plot_gene_search *input1;
    in_plot_gene_search *input2;
    in_plot_gene_search *input3;
    CodRes *results;
    Tcl_DString input_params;
    char tmp[1024];
    char strand_s[7];

    vfuncheader("plot codon pref");

    if (NULL == (input1 = (in_plot_gene_search *)xmalloc
		 (sizeof(in_plot_gene_search))))
	return -1;

    if (NULL == (input2 = (in_plot_gene_search *)xmalloc
		 (sizeof(in_plot_gene_search))))
	return -1;

    if (NULL == (input3 = (in_plot_gene_search *)xmalloc
		 (sizeof(in_plot_gene_search))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    if (end == -1) {
	end = seq_len;
    }

    if (option == 2) {
	strcpy(tmp, get_default_string(interp, nip_defs, w("NIP.PGS.MODE.BUTTON.3")));
    } else if (option == 4) {
	strcpy(tmp, get_default_string(interp, nip_defs, w("NIP.PGS.MODE.BUTTON.4")));
    } else if (option == 6) {
	sprintf(tmp, "%s\n%s\n", get_default_string(interp, nip_defs, w("NIP.PGS.MODE.BUTTON.3")), get_default_string(interp, nip_defs, w("NIP.PGS.MODE.BUTTON.4")));
    } else {
	strcpy(tmp, "");
    }

    win_len *= 3;

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
		       "window length %d strand %s codon table %s\n%s\n",
		       GetSeqName(seq_num), start, end,
		       win_len, strand_s, codon_table, tmp);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input1->params = strdup(Tcl_DStringValue(&input_params)); 
    input2->params = strdup(Tcl_DStringValue(&input_params)); 
    input3->params = strdup(Tcl_DStringValue(&input_params)); 
    input1->method = "codon preference";
    input2->method = "codon preference";
    input3->method = "codon preference";
    Tcl_DStringFree(&input_params);

    if (strand == TOP_S) {
	if (-1 == DoCodonPref(seq, seq_len, codon_table, win_len, 
			      start, end, option, &results)) {
	    verror(ERR_WARN, "CodonPref", "Failed DoCodonPref\n");
	    goto fail;
	}
    } else {
	complement_seq(seq, seq_len);
	if (-1 == DoCodonPref(seq, seq_len, codon_table, win_len, 
			      start, end, option, &results)) {
	    verror(ERR_WARN, "CodonPref", "Failed DoCodonPref\n");
	    goto fail;
	}
	complement_seq(seq, seq_len);	
    }

    id[0] = store_gene_search(seq_num, start, end, 1, input1, 
			      results->frame1, 
			      results->top, results->num_results, 
			      results->min, results->max,
			      SEQ_TYPE_CODONPREF, strand,
			      &graph_obj[0]);
    id[1] = store_gene_search(seq_num, start, end, 2, input2, 
			      results->frame2, 
			      results->top, results->num_results, 
			      results->min, results->max,
			      SEQ_TYPE_CODONPREF, strand,
			      &graph_obj[1]);
    id[2] = store_gene_search(seq_num, start, end, 3, input3, 
			      results->frame3, 
			      results->top, results->num_results, 
			      results->min, results->max,
			      SEQ_TYPE_CODONPREF, strand,
			      &graph_obj[2]);

     free_CodRes(results);
     return 0;

fail:
    xfree(input1);
    xfree(input2);
    xfree(input3);
    return -1;
}

int init_nip_author_test_create(Tcl_Interp *interp,
				int seq_id,
				int start,
				int end,
				char *codon_table,
				double error,
				int strand,
				Tcl_Obj **graph_obj,
				int *id)
{
    char *seq;
    int seq_len, seq_num;
    in_plot_gene_search *input1;
    in_plot_gene_search *input2;
    in_plot_gene_search *input3;
    CodRes *results;
    Tcl_DString input_params;
    char strand_s[7];

    vfuncheader("plot author_test");

    if (NULL == (input1 = (in_plot_gene_search *)xmalloc
		 (sizeof(in_plot_gene_search))))
	return -1;

    if (NULL == (input2 = (in_plot_gene_search *)xmalloc
		 (sizeof(in_plot_gene_search))))
	return -1;

    if (NULL == (input3 = (in_plot_gene_search *)xmalloc
		 (sizeof(in_plot_gene_search))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

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
		       "percent error %f codon table %s\n",
		       GetSeqName(seq_num), start, end,
		       error, codon_table);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input1->method = "author test";
    input2->method = "author test";
    input3->method = "author test";
    input1->params = strdup(Tcl_DStringValue(&input_params));
    input2->params = strdup(Tcl_DStringValue(&input_params));
    input3->params = strdup(Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);
    
    if (strand == TOP_S) {
	if (-1 == DoAuthorTest(seq, seq_len, codon_table, error, 
			       start, end, &results)) {
	    verror(ERR_WARN, "AuthorTest", "Failed DoAuthorTest\n");
	    goto fail;
	}
    } else {
	complement_seq(seq, seq_len);
	if (-1 == DoAuthorTest(seq, seq_len, codon_table, error, 
			       start, end, &results)) {
	    verror(ERR_WARN, "AuthorTest", "Failed DoAuthorTest\n");
	    goto fail;
	}
	complement_seq(seq, seq_len);
    }
    id[0] = store_gene_search(seq_num, start, end, 1, input1, results->frame1, 
			      results->top, results->num_results, 
			      results->min, results->max,
			      SEQ_TYPE_AUTHOR, strand, &graph_obj[0]);
    id[1] = store_gene_search(seq_num, start, end, 2, input2, results->frame2, 
			      results->top, results->num_results, 
			      results->min, results->max,
			       SEQ_TYPE_AUTHOR, strand, &graph_obj[1]);
    id[2] = store_gene_search(seq_num, start, end, 3, input3, results->frame3, 
			      results->top, results->num_results, 
			      results->min, results->max,
			      SEQ_TYPE_AUTHOR, strand, &graph_obj[2]);
     free_CodRes(results);
     return 0;

fail:
    xfree(input1);
    xfree(input2);
    xfree(input3);
    return -1;
}

int init_nip_base_bias_create(Tcl_Interp *interp,
			      int seq_id,
			      int start,
			      int end,
			      int win_len,
			      int strand,
			      Tcl_Obj **graph_obj,
			      int *id)
{
    char *seq;
    int seq_len, seq_num;
    in_plot_gene_search *input;
    CodRes1 *results;
    Tcl_DString input_params;
 
    vfuncheader("plot base bias");

    if (NULL == (input = (in_plot_gene_search *)xmalloc
		 (sizeof(in_plot_gene_search))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    if (end == -1) {
	end = seq_len;
    }

    win_len *= 3;

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "window length %d\n",
		       GetSeqName(seq_num), start, end, win_len);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->method = "base bias";
    input->params = strdup(Tcl_DStringValue(&input_params));
    Tcl_DStringFree(&input_params);
    
    if (strand == TOP_S) {
	if (-1 == DoPosBaseBias(seq, seq_len, win_len, start, end, &results)) {
	    verror(ERR_WARN, "BaseBias", "Failed DoPosBaseBias\n");
	    goto fail;
	}
    } else {
	complement_seq(seq, seq_len);
	if (-1 == DoPosBaseBias(seq, seq_len, win_len, start, end, &results)) {
	    verror(ERR_WARN, "BaseBias", "Failed DoPosBaseBias\n");
	    goto fail;
	}
	complement_seq(seq, seq_len);
    }

    *id = store_gene_search(seq_num, start, end, 0, input, results->frame1, 
			    NULL, results->num_results, 
			    results->min, results->max,
			    SEQ_TYPE_BASEBIAS, strand, graph_obj);

    free_CodRes1(results);
    return 0;
fail:
    xfree(input);
    return -1;
}

int init_nip_gene_search_plot(Tcl_Interp *interp,
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


