#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "seq_reg.h"
#include "seq_results.h"
#include "nip_results.h"
#include "splice_search.h"
#include "nip_splice_search.h"
#include "nip_globals.h"
#include "text_output.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "xalloc.h"
#include "seq_plot_funcs.h"
#include "dna_utils.h"
#include "nip_wtmatrix_search.h"
#include "tclCanvGraph.h"
#include "element_canvas.h"
#include "seq_element.h"
#include "misc.h"

void splice_search_callback(int seq_num, void *obj, seq_reg_data *jdata);

void splice_search_shutdown(Tcl_Interp *interp,
			    seq_result *s_result,
			    element *e,
			    int seq_num)
{
    seq_reg_key_name info;
    static char buf[80];
    in_splice *input = s_result->input;

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);
	
    /* need to deregister sequence */
    seq_deregister(seq_num, splice_search_callback, (seq_result *)s_result);

    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
   
     if (e->num_results > 1) {
	e->replot_func(e);
    }

    free(input->params);
    xfree(s_result->input);
    if (s_result->text_data) {
	if (((text_wtmatrix **)s_result->text_data)[0])
	    xfree(((text_wtmatrix **)s_result->text_data)[0]);
	xfree(s_result->text_data);
    }

    xfree(s_result);
}

void splice_search_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_splice *input = s_result->input;
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
	sprintf(jdata->name.line, "Splice search");
	break;    

    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "splice f%d #%d", s_result->frame,
		result_id);
	break;
	
    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "splice search: seq=%s frame=%d", 
		GetSeqName(GetSeqNum(s_result->seq_id[HORIZONTAL])), 
		s_result->frame);
	break;

    case SEQ_GET_BRIEF_TAG:
	{
	      Graph *graph;
	      Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
	      text_wtmatrix **text_data = s_result->text_data;
	      int i, j;
	      char type[10];
	      char *sequence;
	      int seq_num;

	      graph = Tcl_GetGraphFromObj(graph_obj);
	      seq_num = GetSeqNum(s_result->seq_id[HORIZONTAL]);
	      sequence = GetSeqSequence(seq_num);

	      if (strcmp(jdata->brief_tag.array_type, "d") == 0) {
		  i = jdata->brief_tag.arrays;
		  j = jdata->brief_tag.array;

		  if (s_result->type == SEQ_TYPE_SPLICE_DONOR) {
		      sprintf(type, "donor");
		  } else {
		      sprintf(type, "acceptor");
		  }

		  sprintf(jdata->brief_tag.line, "%s frame %d pos %d score %f %.*s", 
			  type, s_result->frame, 
			  (int)graph->d_arrays[0].d_array[j].x0, 
			  graph->d_arrays[0].d_array[j].y1, 
			  text_data[0]->length, 
			  &sequence[(int)graph->d_arrays[0].d_array[j].x0 - 1 - 
				   text_data[0]->mark_pos]);
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
	    vfuncheader("splice search results frame %d", s_result->frame);
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
		splice_search_shutdown(interp, s_result, e, seq_num);
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
					   w("ELEMENT.SINGLE.PLOT_HEIGHT"));
		jdata->info.result = (void *)&pt;
		break;
	    }
	}
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: 
	{
	    splice_search_shutdown(interp, s_result, e, seq_num);
	    break;
	}
    }
}

void
splice_search_text_func(void *obj)
{
    seq_result *s_result = (seq_result *) obj;
    text_wtmatrix **text_data = s_result->text_data;
    int i;
    char *sequence;
    int seq_num;
    Graph *graph;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;

    graph = Tcl_GetGraphFromObj(graph_obj);

    seq_num = GetSeqNum(s_result->seq_id[HORIZONTAL]);
    sequence = GetSeqSequence(seq_num);

    if (s_result->type == SEQ_TYPE_SPLICE_DONOR) 
	vmessage("Donor\n");
    else 
	vmessage("Acceptor\n");

    for (i = 0; i < graph->d_arrays[0].n_dlines; i++) {
	UpdateTextOutput();
	vmessage("Position %8d %8d score %f %.*s\n", 
		(int)graph->d_arrays[0].d_array[i].x0 - text_data[0]->mark_pos, 
		 (int)graph->d_arrays[0].d_array[i].x0 + 1, 
		 graph->d_arrays[0].d_array[i].y1, 
		 text_data[0]->length, 
		 &sequence[(int)graph->d_arrays[0].d_array[i].x0 - 1 - 
			  text_data[0]->mark_pos]);
    }
}

int store_splice_search(int seq_num,
			WtmatrixRes *matrix,
			in_splice *input,
			int start,
			int end,
			int frame,
			int strand,
			int type,
			Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    Graph *graph;
    int id, i;
    text_wtmatrix **text_data;

    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;

    Tcl_InitGraph(&graph);

    if (NULL == (text_data = (text_wtmatrix **)xmalloc(sizeof(text_wtmatrix*))))
	return -1;
 
    if (NULL == (text_data[0] = (text_wtmatrix *)xmalloc(sizeof(text_wtmatrix))))
	return -1;
    
     if (NULL == (graph->d_arrays = (darray *)xmalloc(sizeof(darray))))
	 return -1;

     if (NULL == (graph->d_arrays[0].d_array = (gd_line *)xmalloc(matrix->number_of_res * sizeof(gd_line))))
	 return -1;

    graph->d_arrays[0].n_dlines = matrix->number_of_res;
    graph->d_arrays[0].type = G_LINES;
    graph->n_darrays = 1;
    graph->dim.x0 = start;
    graph->dim.x1 = end;

    graph->dim.y0 = matrix->min;
    graph->dim.y1 = 2 * matrix->max;

    for (i = 0; i < matrix->number_of_res; i++) {
	graph->d_arrays[0].d_array[i].x0 = matrix->match[i]->pos + 1;
	graph->d_arrays[0].d_array[i].x1 = matrix->match[i]->pos + 1;
	graph->d_arrays[0].d_array[i].y0 = 0.0;
	graph->d_arrays[0].d_array[i].y1 = (double)matrix->match[i]->score;
    }    
      
    *graph_obj = Tcl_NewGraphObj(graph);
    /* these are copied over in Tcl_NewGraphObj so don't need them anymore */
    xfree(graph->d_arrays[0].d_array);
    xfree(graph->d_arrays);
    xfree(graph);

    id = get_reg_id();
    s_result->text_data = text_data; 
    text_data[0]->length = matrix->length;
    text_data[0]->mark_pos = matrix->mark_pos;

    s_result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    s_result->seq_id[VERTICAL] = -1;

    s_result->id = id; 
    s_result->input = (void *)input; 
    s_result->type = type;
    s_result->gr_type = SEQ_TYPE_GRAPH_PLOT;
    s_result->frame = frame;
    s_result->graph = CANV_LINE;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */
    s_result->strand = strand;

    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = splice_search_callback;
    s_result->txt_func = splice_search_text_func;
    
    seq_register(seq_num, splice_search_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);

    if (matrix) free_WtmatrixRes (matrix);
    return id;
}

int init_splice_search_create(int seq_id,
			      int start,
			      int end,
			      char *donor,
			      char *acceptor,
			      int strand,
			      Tcl_Obj **graph_obj,
			      int *id) /* out */
{
    in_splice *input1, *input2, *input3;
    in_splice *input4, *input5, *input6;
    char *seq;
    int seq_len;
    int seq_num;
    SpliceResults *splice_result;
    Tcl_DString input_params;
    int irs;
    char strand_s[7];

    vfuncheader("splice search");
    set_char_set(DNA);

    if (NULL == (input1 = (in_splice *)xmalloc (sizeof(in_splice))))
	return -1;

    if (NULL == (input2 = (in_splice *)xmalloc (sizeof(in_splice))))
	return -1;

    if (NULL == (input3 = (in_splice *)xmalloc (sizeof(in_splice))))
	return -1;

    if (NULL == (input4 = (in_splice *)xmalloc (sizeof(in_splice))))
	return -1;

    if (NULL == (input5 = (in_splice *)xmalloc (sizeof(in_splice))))
	return -1;

    if (NULL == (input6 = (in_splice *)xmalloc (sizeof(in_splice))))
	return -1;

    if (NULL == (splice_result = (SpliceResults *)xmalloc 
		 (sizeof(SpliceResults))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }

    if (strand == TOP_S) {
	irs = splice_search(seq, seq_len, start, end, donor, acceptor, 
			    splice_result);
    } else {
	complement_seq(seq, seq_len);
	irs = splice_search(seq, seq_len, start, end, donor, acceptor, 
			    splice_result);
	complement_seq(seq, seq_len);	
    }

    if (irs == -1) {
      xfree(splice_result);
      xfree(input1);
      xfree(input2);
      xfree(input3);
      xfree(input4);
      xfree(input5);
      xfree(input6);
      verror(ERR_WARN, "splice search",
	     "error in splice search (maybe none found)");
      return -1;
    }

    if (splice_result->ied_f1->number_of_res == 0 &&
	splice_result->ied_f2->number_of_res == 0 &&
	splice_result->ied_f3->number_of_res == 0 &&
	splice_result->eia_f1->number_of_res == 0 &&
	splice_result->eia_f2->number_of_res == 0 &&
	splice_result->eia_f3->number_of_res == 0) {
	verror(ERR_WARN, "splice search", "no matches found");
	xfree(splice_result);
	xfree(input1);
	xfree(input2);
	xfree(input3);
	xfree(input4);
	xfree(input5);
	xfree(input6);
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
		       "strand %s donor weight matrix %s\nacceptor weight matrix %s\n",
		       GetSeqName(seq_num), start, end, strand_s, donor, 
		       acceptor);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input1->params = strdup(Tcl_DStringValue(&input_params)); 
    input2->params = strdup(Tcl_DStringValue(&input_params)); 
    input3->params = strdup(Tcl_DStringValue(&input_params)); 
    input4->params = strdup(Tcl_DStringValue(&input_params)); 
    input5->params = strdup(Tcl_DStringValue(&input_params)); 
    input6->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (-1 == (id[0] = store_splice_search(seq_num, splice_result->ied_f1,
					   input1, 
					   start, end, 1, strand, SEQ_TYPE_SPLICE_DONOR, &graph_obj[0]))){
	verror(ERR_FATAL,"nip splice search", "error in saving matches\n");
	return -1;
    }
    if (-1 == (id[1] = store_splice_search(seq_num, splice_result->ied_f2, 
					   input2, 
					   start, end, 2, strand, SEQ_TYPE_SPLICE_DONOR, &graph_obj[1]))){
	verror(ERR_FATAL,"nip splice search", "error in saving matches\n");
	return -1;
    }


    if (-1 == (id[2] = store_splice_search(seq_num, splice_result->ied_f3,
					   input3, 
					   start, end, 3, strand, SEQ_TYPE_SPLICE_DONOR, &graph_obj[2]))){
	verror(ERR_FATAL,"nip splice search", "error in saving matches\n");
	return -1;
    }

    if (-1 == (id[3] = store_splice_search(seq_num, 
					   splice_result->eia_f1, input4, 
					   start, end, 1, strand, SEQ_TYPE_SPLICE_ACCEPTOR, &graph_obj[3]))){
	verror(ERR_FATAL,"nip splice search", "error in saving matches\n");
	return -1;
    }

    if (-1 == (id[4] = store_splice_search(seq_num, splice_result->eia_f2, 
					   input5, 
					   start, end, 2, strand, SEQ_TYPE_SPLICE_ACCEPTOR, &graph_obj[4]))){
	verror(ERR_FATAL,"nip splice search", "error in saving matches\n");
	return -1;
    }

    if (-1 == (id[5] = store_splice_search(seq_num, 
					   splice_result->eia_f3, input6, 
					   start, end, 3, strand, SEQ_TYPE_SPLICE_ACCEPTOR, &graph_obj[5]))){
	verror(ERR_FATAL,"nip splice search", "error in saving matches\n");
	return -1;
    }
    xfree(splice_result);
    return 0;
}

int init_splice_search_plot(Tcl_Interp *interp,
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

    if (s_result->type == SEQ_TYPE_SPLICE_DONOR) {
	configure->position = 1.0;
	configure->x_direction = '+';
	configure->y_direction = '+';
	configure->height = tick_ht;
	configure->zoom = 1;
	configure->scroll = 1;
    } else {
	configure->position = 0.0;
	configure->x_direction = '+';
	configure->y_direction = '-';
	configure->height = tick_ht;
	configure->zoom = 1;
	configure->scroll = 1;
    }

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
