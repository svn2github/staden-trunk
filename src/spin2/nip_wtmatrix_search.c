#include <tcl.h>
#include <string.h>

#include "xalloc.h"
#include "seq_reg.h"
#include "seq_results.h"
#include "nip_results.h"
#include "text_output.h"
#include "nip_globals.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "dna_utils.h"
#include "splice_search.h"
#include "nip_wtmatrix_search.h"
#include "misc.h"
#include "seq_plot_funcs.h"
#include "tclCanvGraph.h"
#include "seq_element.h"

void nip_wtmatrix_search_callback(int seq_num, void *obj, seq_reg_data *jdata);

void nip_wtmatrix_search_shutdown(Tcl_Interp *interp,
				  seq_result *s_result,
				  element *e,
				  int seq_num)
{
    in_wtmatrix_search *input = s_result->input;
    seq_reg_key_name info;
    static char buf[80];

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    seq_deregister(seq_num, nip_wtmatrix_search_callback, 
		   (seq_result *)s_result);
		
    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
   
    if (e->num_results > 1) {
	e->replot_func(e);
    }
 
    free(input->params);
    xfree(s_result->input);
    xfree(s_result->text_data);
    xfree(s_result);
}

void nip_wtmatrix_search_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_wtmatrix_search *input = s_result->input;
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
	sprintf(jdata->name.line, "wtmatrix search");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "wtmatrix #%d", result_id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "wtmatrix: seq=%s", 
		GetSeqName(GetSeqNum(s_result->seq_id[HORIZONTAL])));
	break;
	
    case SEQ_GET_BRIEF_TAG:
	{
	      Graph *graph;
	      Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
	      text_wtmatrix *text_data = s_result->text_data;
	      int i, j;
	      char *sequence;
	      int seq_num, seq_len;

	      graph = Tcl_GetGraphFromObj(graph_obj);
	      seq_num = GetSeqNum(s_result->seq_id[HORIZONTAL]);
	      sequence = GetSeqSequence(seq_num);
	      seq_len = GetSeqLength(seq_num);

	      if (s_result->strand == BOTTOM_S) {
		  complement_seq(sequence, seq_len);
	      }
	      if (strcmp(jdata->brief_tag.array_type, "d") == 0) {
		  i = jdata->brief_tag.arrays;
		  j = jdata->brief_tag.array;
		  sprintf(jdata->brief_tag.line, "pos %d score %f %.*s", 
			  (int)graph->d_arrays[i].d_array[j].x0, 
			  graph->d_arrays[i].d_array[j].y1, text_data->length, 
			  &sequence[(int)graph->d_arrays[i].d_array[j].x0 - 1 - 
			  text_data->mark_pos]);
	      }
	      if (s_result->strand == BOTTOM_S) {
		  complement_seq(sequence, seq_len);
	      }
	    
	break;
	}
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
		nip_wtmatrix_search_shutdown(interp, s_result, e, seq_num);
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
					      w("ELEMENT.SINGLE.PLOT_HEIGHT"));
		    jdata->info.result = (void *)&pt;
		    break;
		}
	    }
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: {
	nip_wtmatrix_search_shutdown(interp, s_result, e, seq_num);
	break;
    }
    }
}

void
nip_wtmatrix_search_text_func(void *obj)
{
    int i;
    seq_result *s_result = (seq_result *) obj;
    text_wtmatrix *text_data = s_result->text_data;
    char *seq_name;
    char *sequence;
    int seq_num, seq_len;
    Graph *graph;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;

    seq_num = GetSeqNum(s_result->seq_id[HORIZONTAL]);
    seq_name = GetSeqName(seq_num);
    sequence = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);


    graph = Tcl_GetGraphFromObj(graph_obj);

    if (s_result->strand == BOTTOM_S) {
	complement_seq(sequence, seq_len);
    }
    for (i = 0; i < graph->d_arrays[0].n_dlines; i++) {
	UpdateTextOutput();

	vmessage("Position %8d %8d score %f %.*s\n", 
		 (int)graph->d_arrays[0].d_array[i].x0 - text_data->mark_pos, 
		 (int)graph->d_arrays[0].d_array[i].x0 + 1, 
		 graph->d_arrays[0].d_array[i].y1, 
		 text_data->length, 
		 &sequence[(int)graph->d_arrays[0].d_array[i].x0 - 1 - 
			  text_data->mark_pos]);
    }
    if (s_result->strand == BOTTOM_S) {
	complement_seq(sequence, seq_len);
    }
}

int store_wtmatrix_search(int seq_num,
			  in_wtmatrix_search *input,
			  int start,
			  int end,
			  int strand,
			  WtmatrixRes *results,
			  Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    int id, i, j;
    text_wtmatrix *text_data;
    Graph *graph;
    int seq_len = end - start + 1;

    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;
    
    Tcl_InitGraph(&graph);
    
    if (NULL == (graph->d_arrays = (darray *)xmalloc(sizeof(darray))))
	 return -1;

     if (NULL == (graph->d_arrays[0].d_array = (gd_line *)xmalloc(results->number_of_res * sizeof(gd_line))))
	 return -1;

    if (NULL == (text_data = (text_wtmatrix *)xmalloc(sizeof(text_wtmatrix))))
	return -1;
    
    graph->d_arrays[0].n_dlines = results->number_of_res;
    graph->d_arrays[0].type = G_LINES;
    graph->n_darrays = 1;
    graph->dim.x0 = start;
    graph->dim.x1 = end;
    graph->dim.y0 = results->min;
    graph->dim.y1 = results->max;
    
    if (strand == TOP_S) {
	for (i = 0; i < results->number_of_res; i++) {
	    graph->d_arrays[0].d_array[i].x0 = results->match[i]->pos + 1;
	    graph->d_arrays[0].d_array[i].x1 = results->match[i]->pos + 1;
	    graph->d_arrays[0].d_array[i].y0 = 0.0;
	    graph->d_arrays[0].d_array[i].y1 = (double)results->match[i]->score;
	}    
    } else {
	for (i = 0, j = results->number_of_res -1 ; 
	     i < results->number_of_res; i++, j--) {
	    graph->d_arrays[0].d_array[i].x0 = seq_len - (results->match[j]->pos + 1);
	    graph->d_arrays[0].d_array[i].x1 = seq_len - (results->match[j]->pos + 1);
	    graph->d_arrays[0].d_array[i].y0 = 0.0;
	    graph->d_arrays[0].d_array[i].y1 = (double)results->match[j]->score;
	    printf("i %d x0 %f %d y1 %f\n", i, graph->d_arrays[0].d_array[i].x0,
		   results->mark_pos, graph->d_arrays[0].d_array[i].y1);

	}    
    }
    
    *graph_obj = Tcl_NewGraphObj(graph);
    /* these are copied over in Tcl_NewGraphObj so don't need them anymore */
    xfree(graph->d_arrays[0].d_array);
    xfree(graph->d_arrays);
    xfree(graph);

    id = get_reg_id();
    
    s_result->text_data = text_data; 
    text_data->length = results->length;
    text_data->mark_pos = results->mark_pos;
    
    s_result->id = id; 
    s_result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    s_result->seq_id[VERTICAL] = -1;

    s_result->input = (void *)input; 
    s_result->type = SEQ_TYPE_WTMATRIXSEARCH;
    s_result->gr_type = SEQ_TYPE_GRAPH_PLOT;
    s_result->frame = 0;
    s_result->graph = CANV_LINE;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */
    s_result->strand = strand;

    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = nip_wtmatrix_search_callback;
    s_result->txt_func = nip_wtmatrix_search_text_func;
    
    seq_register(seq_num, nip_wtmatrix_search_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);

    if (results) free_WtmatrixRes (results);
    return id;
}

int init_nip_wtmatrix_search_create(int start,
				    int end,
				    int seq_id,
				    char *wt_matrix,
				    int strand,
				    Tcl_Obj **graph_obj,
				    int *id)
{
    in_wtmatrix_search *input;
    char *seq;
    int seq_len;
    int seq_num;
    WtmatrixRes *results = NULL;
    Tcl_DString input_params;
    int irs;
    char strand_s[7];

    vfuncheader("weight matrix search");
    set_char_set(DNA);

    if (NULL == (input = (in_wtmatrix_search *)xmalloc
		 (sizeof(in_wtmatrix_search))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }

    if (strand == TOP_S) {
	irs = weight_search(seq, seq_len, start, end, wt_matrix, &results);
	if (irs == -1) {
	    verror(ERR_WARN, "weight matrix search", "error in weight matrix search");
	    return -1;
	}
    } else {
	complement_seq(seq, seq_len);
	irs = weight_search(seq, seq_len, start, end, wt_matrix, &results);
	if (irs == -1) {
	    verror(ERR_WARN, "weight matrix search", "error in weight matrix search");
	    return -1;
	}
	complement_seq(seq, seq_len);
    }

    if (results->number_of_res == 0) {
	verror(ERR_WARN, "weight matrix search", "no matches found");
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
		       "weight matrix %s\n",
		       GetSeqName(seq_num), start, end, wt_matrix);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);
    
    if (-1 == (*id = store_wtmatrix_search(seq_num, input, start, end, strand,
					   results, graph_obj))){
	verror(ERR_FATAL,"weight matrix search", "error in saving matches\n");
	return -1;
    }

  return 0;
}

int init_nip_wtmatrix_search_plot(Tcl_Interp *interp, 
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
    configs *configure;
    seq_result *s_result;
    Graph *graph;
    plot_data *result;
    int seq_id_h, seq_id_v;
    
    /* no results found */
    if (result_id == -1) {
	return 0;
    }
    
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
		     container_id, 
		     element_id, e_win, c_win, orientation, HORIZONTAL|VERTICAL,
		     graph, element_type);

    return 0;
}


