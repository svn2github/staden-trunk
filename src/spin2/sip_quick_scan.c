#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "sip_hash.h"
#include "seq_results.h"
#include "sip_results.h"
#include "compare_spans.h"
#include "sip_find_identity.h"
#include "misc.h"
#include "sequence_formats.h"    /* PROTEIN DNA */
#include "tcl_utils.h"
#include "sip_globals.h"
#include "dna_utils.h"
#include "text_output.h"
#include "readpam.h"
#include "sip_quick_scan.h"
#include "seq_plot_funcs.h"
#include "sequence_pair_display.h"
#include "tclCanvGraph.h"
#include "nip_results.h"
#include "seq_reg.h"
#include "seq_element.h"
#include "element_canvas.h"

#include "restriction_enzyme_map.h"
#include "nip_restriction_enzymes.h"

void quick_scan_callback(int seq_num, void *obj, seq_reg_data *jdata);

void quick_scan_shutdown(Tcl_Interp *interp,
			 seq_result *s_result,
			 element *e)
{
    seq_reg_key_name info;
    static char buf[80];
    in_quick_scan *input = s_result->input;

#ifdef DEBUG
    printf("START quick_scan_delete_result\n");
#endif

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

     /* deregister horizontal sequence */
    seq_deregister(GetSeqNum(s_result->seq_id[HORIZONTAL]), quick_scan_callback, (seq_result *)s_result);

    /* deregister vertical sequence */
    seq_deregister(GetSeqNum(s_result->seq_id[VERTICAL]), quick_scan_callback, (seq_result *)s_result);

    if (e->num_results > 1) {
	Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	e->replot_func(e);
    } else {
	Tcl_VarEval(e->c->interp, "result_list_update ", NULL);
    }
    DestroySequencePairDisplay(interp, s_result->id);
    free(input->params);
}

void quick_scan_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_quick_scan *input = s_result->input;
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
	sprintf(jdata->name.line, "Find best diagonals");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "diagonals #%d", result_id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "diagonals: hori=%s vert=%s", 
		GetSeqBaseName(GetSeqNum(s_result->seq_id[HORIZONTAL])),
		GetSeqBaseName(GetSeqNum(s_result->seq_id[VERTICAL])));
	break;

    case SEQ_GET_OPS:
	if (result->hidden) {
	    jdata->get_ops.ops = "Information\0List results\0PLACEHOLDER\0"
		"PLACEHOLDER\0PLACEHOLDER\0Reveal\0SEPARATOR\0Remove\0";
	} else if (!get_replot_temp()) {
	    jdata->get_ops.ops = "Information\0PLACEHOLDER\0PLACEHOLDER\0"
	       "Display sequences\0PLACEHOLDER\0PLACEHOLDER\0SEPARATOR\0Remove\0";
	} else {
	    jdata->get_ops.ops = "Information\0List results\0Configure\0"
	       "Display sequences\0Hide\0PLACEHOLDER\0SEPARATOR\0Remove\0";
	}
	break;
    case SEQ_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0: /* information */
	    vfuncheader("input parameters");
	    vmessage("%s\n", input->params);
	    break;
	case 1: /* list results */
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
	case 3: /* display sequences */
	    SequencePairDisplay(interp, e->win, result_id, 
				s_result->seq_id[HORIZONTAL], 
				s_result->seq_id[VERTICAL]);
	    break;
	case 4: /* hide all */
	    result->hidden = 1;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 5: /* reveal all */
	    result->hidden = 0;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 6: /* remove */
	    {
		quick_scan_shutdown(interp, s_result, e);
		remove_result_from_element(e, result_id);

	    }
	    break;
	}
	break;
    case SEQ_PLOT:
	{
	    int save_results = 1;

	    /* 
	     * save_results is always true. If need to go back to the 
	     * previous code, must look at the backups (1999.1p1)
	     */
	    if (save_results) {
		s_result->pr_func(s_result, NULL);
	    } else {
	    }
	    break;
	}
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case RESULT:
	    jdata->info.result = (void *)s_result;
	    break;
	case WIN_NAME:
	    {
		char *r_win =  e->win;
		jdata->info.result = (void *)r_win;
		break;
	    }
	case WIN_SIZE:
	    {
		static d_point pt;
		pt.x = get_default_int(interp, tk_utils_defs, 
					w("DOT.PLOT_WIDTH"));
		pt.y = get_default_double(interp, tk_utils_defs,
					   w("DOT.PLOT_HEIGHT"));

		jdata->info.result = (void *)&pt;
		break;
	    } /* WIN_SIZE */
	}
	break;
    case SEQ_QUIT:
    case SEQ_DELETE:
	{    
	    quick_scan_shutdown(interp, s_result, e);
	}
    }
}

void quick_scan_text_func(void *obj)
{
    seq_result *s_result = (seq_result *) obj;
    int i;
    char *seq1;
    char *seq2;
    int seq1_len, seq2_len;
    int max_str;
    int seq_num_h, seq_num_v;
    Graph *graph;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;    

    graph = Tcl_GetGraphFromObj(graph_obj);

    seq_num_h = GetSeqNum(s_result->seq_id[HORIZONTAL]);
    seq_num_v = GetSeqNum(s_result->seq_id[VERTICAL]);

    seq1 = GetSeqSequence(seq_num_h);
    seq1_len = GetSeqLength(seq_num_h);
    seq2 = GetSeqSequence(seq_num_v);
    seq2_len = GetSeqLength(seq_num_v);

    if (seq1_len > seq2_len)
	max_str = seq1_len;
    else
	max_str = seq2_len;

    for (i = 0; i < graph->p_arrays[0].n_pts; i++) {
	UpdateTextOutput();
	vmessage("Positions %10d h %10d v \n", 
	       (int)graph->p_arrays[0].p_array[i].x, 
		 (int)graph->p_arrays[0].p_array[i].y);
    }
}

void quick_scan_recalc_func(void *obj, seq_reg_plot *plot)
{
}

int store_quick_scan(int seq1_num,
		     int seq2_num,
		     int start_h,
		     int end_h,
		     int start_v,
		     int end_v,
		     in_quick_scan *input,
		     int *seq1_match, 
		     int *seq2_match, 
		     int num_elements,
		     int save_results,
		     int strand,
		     Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    int i, id;
    Graph *graph;

    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (save_results) {

	if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	    return -1;
    
	Tcl_InitGraph(&graph);
	
	if (NULL == (graph->p_arrays = (parray *)xmalloc(sizeof(parray))))
	    return -1;

	if (NULL == (graph->p_arrays[0].p_array = (g_pt *)xmalloc(num_elements * sizeof(g_pt))))
	    return -1;

	/* data assignment */
	for (i = 0; i < num_elements; i++) {
	    graph->p_arrays[0].p_array[i].x = seq1_match[i];
	    graph->p_arrays[0].p_array[i].y = seq2_match[i];

	}
	/* array must be in ascending x order */
	qsort((void *) graph->p_arrays[0].p_array, num_elements, sizeof(g_pt), 
	      compare_g_pt);
	
#ifdef DEBUG
	for (i = 0; i < num_elements; i++) {
	    printf("x %f y %f\n",  graph->p_arrays[0].p_array[i].x,  graph->p_arrays[0].p_array[i].y);
	}
#endif
	graph->p_arrays[0].n_pts = num_elements;
	graph->p_arrays[0].type = G_DOT;
	graph->n_parrays = 1;
	graph->dim.x0 = start_h;
	graph->dim.x1 = end_h;
	graph->dim.y0 = start_v;
	graph->dim.y1 = end_v;
	*graph_obj = Tcl_NewGraphObj(graph);

	/* these are copied over in Tcl_NewGraphObj so don't need them anymore */
	xfree(graph->p_arrays[0].p_array);
	xfree(graph->p_arrays);
	xfree(graph);

    } else {
    }
    id = get_reg_id();

    s_result->seq_id[HORIZONTAL] = GetSeqId(seq1_num);
    s_result->seq_id[VERTICAL] = GetSeqId(seq2_num);

    s_result->id = id; 
    s_result->input = (void *)input; 
    s_result->type = SEQ_TYPE_DOT_PLOT;
    s_result->gr_type = SEQ_TYPE_DOT_PLOT;
    s_result->frame = 0;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */

    s_result->op_func = quick_scan_callback;
    s_result->txt_func = quick_scan_text_func;
    
    if (save_results) {
	s_result->pr_func = seq_plot_graph_func;
	seq_register(seq1_num, quick_scan_callback, (void *)s_result, 
		     SEQ_PLOT_PERM, id);
	seq_register(seq2_num, quick_scan_callback, (void *)s_result, 
		     SEQ_PLOT_PERM, id);
    } else {
	s_result->pr_func =  quick_scan_recalc_func;
	seq_register(seq1_num, quick_scan_callback, (void *)s_result, 
		     SEQ_PLOT_TEMP, id);
	seq_register(seq2_num, quick_scan_callback, (void *)s_result, 
		     SEQ_PLOT_TEMP, id);
    }
    return id;
}

int init_sip_best_diagonals_create(Tcl_Interp *interp, 
				   int seq_id_h,
				   int seq_id_v, 
				   int start_h,
				   int end_h, 
				   int start_v,
				   int end_v, 
				   int win_len,
				   int min_match,
				   int word_len,
				   float min_sd,
				   int strand,
				   Tcl_Obj **graph_obj,
				   int *id)
{
    char *seq1;
    char *seq2;
    int seq1_len, seq2_len;
    int same_seq;
    int max_matches = get_max_matches();
    int seq1_num, seq2_num;
    int seq1_type, seq2_type;
    int sub1_len, sub2_len;
    in_quick_scan *input = NULL;
    int n_matches;
    int *seq1_match = NULL;
    int *seq2_match = NULL; 
    Tcl_DString input_params;
    char strand_s[7];

    /* 
     * FIXME: decide whether to implement this - if false, will only plot 
     * results, & not save them ie zooming, scrolling, moving etc will not work
     */
    int save_results = 1;

    vfuncheader("Find best diagonals");

     if (-1 == (seq1_num = GetSeqNum(seq_id_h))) {
	verror(ERR_WARN, "find best diagonals", 
	       "horizontal sequence undefined");
	goto error;
    }
     
    if (-1 == (seq2_num = GetSeqNum(seq_id_v))) {
	verror(ERR_WARN, "find best diagonals", 
	       "vertical sequence undefined");
	goto error;
    }
   
    if (NULL == (input = (in_quick_scan *)xmalloc(sizeof(in_quick_scan))))
	goto error;

     
    seq1 = GetSeqSequence(seq1_num);
    seq1_len = GetSeqLength(seq1_num);
    seq2 = GetSeqSequence(seq2_num);
    seq2_len = GetSeqLength(seq2_num);
    seq1_type = GetSeqType(seq1_num);
    seq2_type = GetSeqType(seq2_num);

     if (end_h == -1)
	end_h = seq1_len;

    if (end_v == -1)
	end_v = seq2_len;

    if (seq1_type != seq2_type) {
	verror(ERR_WARN, "quick scan", "sequences must both be either DNA or protein");
	return TCL_OK;
    } else if (seq1_type == PROTEIN) {
	set_char_set(PROTEIN);
        set_score_matrix(get_matrix_file(PROTEIN));
    } else if (seq1_type == DNA) {
	set_char_set(DNA);
        set_score_matrix(get_matrix_file(DNA));
    }
    /*
     * first check if seq lengths are equal, if not the seqs cannot be the
     * same
     */

    sub1_len = end_h - start_h + 1;
    sub2_len = end_v - start_v + 1;

    if (sub1_len == sub2_len) {
	if (strncmp(seq1 + start_h - 1, seq2 + start_v - 1, sub1_len) == 0) {
	    same_seq = 1;
	} else {
	    same_seq = 0;
	}
    } else {
	same_seq = 0;
    }
    if (!get_remove_dup() && same_seq)
	same_seq = 0;

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    if (strand & TOP_S) {
        strcpy(strand_s, "top");
    } else if (strand & BOTTOM_S) {
        strcpy(strand_s, "bottom");
    } else {
        strcpy(strand_s, "both");
    }

    vTcl_DStringAppend(&input_params, "horizontal %s: %s\nvertical %s: %s\n"
		       "strand %s window length %d minimum score %d word length %d minimum sd %f", 
		       GetSeqLibraryName(seq1_num), 
		       GetSeqName(seq1_num), 
		       GetSeqLibraryName(seq2_num), 
		       GetSeqName(seq2_num), strand_s, win_len, min_match, 
		       word_len, min_sd); 
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);
    
    if (NULL == (seq1_match = (int *)xmalloc(max_matches * sizeof(int))))
      goto error;
    if (NULL == (seq2_match = (int *)xmalloc(max_matches * sizeof(int))))
      goto error;

    if (save_results) {
	set_replot_temp(1);
	
	if (strand == TOP_S) {
	    n_matches = quick_scan(seq1, seq2, start_h, end_h, start_v, end_v, 
				   seq1_type, max_matches, same_seq, win_len, 
				   min_match, word_len, min_sd, save_results,
				   &seq1_match, &seq2_match, NULL, NULL);
	    if (n_matches == -1) {
		goto error;
	    }
	} else {
	    if (seq2_type == DNA)
		complement_seq(seq2, seq2_len);

	    n_matches = quick_scan(seq1, seq2, start_h, end_h, start_v, end_v, 
				   seq1_type, max_matches, same_seq, win_len, 
				   min_match, word_len, min_sd, save_results,
				   &seq1_match, &seq2_match, NULL, NULL);
	    if (seq2_type == DNA)
		complement_seq(seq2, seq2_len);

	    if (n_matches == -1) {
		goto error;
	    }
	}
    }

    if (n_matches > 0) {
 	if (-1 == (*id = store_quick_scan(seq1_num, seq2_num, start_h, 
					  end_h, start_v, end_v, input, 
					  seq1_match, seq2_match, 
					  n_matches, save_results, strand,
					  graph_obj))) {
	    goto error;
	}
    } else {
	verror(ERR_WARN, "Find best diagonals", "no matches found\n"); 
	if (seq1_match)
	    xfree(seq1_match);
	if (seq2_match) 
	    xfree(seq2_match);
	return -1;
    }

    if (seq1_match)
      xfree(seq1_match);
    if (seq2_match) 
      xfree(seq2_match);
    return 0;

 error:
    verror(ERR_WARN, "Find best diagonals", "failure in find best diagonals");
    if (seq1_match)
      xfree(seq1_match);
    if (seq2_match) 
      xfree(seq2_match);
    return -1;
    

}

int init_sip_best_diagonals_plot(Tcl_Interp *interp,
				 int seq_id_h,
				 int seq_id_v,
				 int result_id,
				 char *e_win, 
				 int element_id,
				 char *c_win, 
				 int container_id,
				 Tcl_Obj *results,
				 int line_width,
				 char *colour,
				 char *element_type)

{
    seq_result *s_result;
    Graph *graph;
    plot_data *result;
    configs *configure;
 
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
    configure->height = 1.0;
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
    result->amp_ruler = 0;
    sprintf(result->tags, "id%d", result_id);

    graph = Tcl_GetGraphFromObj(results);

    init_seq_element(interp, s_result, seq_id_h, seq_id_v, result, 
		     container_id, element_id, e_win, c_win, 
		     HORIZONTAL|VERTICAL, HORIZONTAL|VERTICAL, graph, 
		     element_type);
    return 0;
}
