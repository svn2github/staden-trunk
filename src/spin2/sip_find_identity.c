#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "seq_reg.h"
#include "sip_hash.h"
#include "probs.h"
#include "seq_results.h"
#include "sip_results.h"
#include "sip_find_identity.h"
#include "misc.h"
#include "sequence_formats.h"
#include "text_output.h"
#include "readpam.h"
#include "tcl_utils.h"
#include "sip_globals.h"
#include "dna_utils.h"
#include "seq_plot_funcs.h"
#include "sequence_pair_display.h"
#include "tclCanvGraph.h"
#include "seq_element.h"

void identities_callback(int seq_num, void *obj, seq_reg_data *jdata);

/*
 * score, probabilites, expected and observed scores output in the
 * text output window
 */
int
CalcIdentityProbs(seq_result *s_result, int word_len)
{
    char *seq1;
    char *seq2;
    int seq1_type, seq2_type;
    int seq_num1, seq_num2;
    int i, j;
    int max = 0;
    int cum = 0;
    Graph *graph;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
    int *score_hist;
    int score;

     graph = Tcl_GetGraphFromObj(graph_obj);

    for (i = 0; i < graph->d_arrays[0].n_dlines; i++) {
	score = graph->d_arrays[0].d_array[i].x1 - graph->d_arrays[0].d_array[i].x0;

#ifdef DEBUG
	printf("i %d score %d \n", i, score);
#endif
	if (score > max) 
	    max = score;
    }
#ifdef DEBUG
    printf("min %d max %d \n", word_len, max);
#endif
    if (NULL == (score_hist = (int *)xcalloc(max - word_len + 1, sizeof(int))))
	return -1;
    
    /* create histogram of scores between word_len and max */
    for (i = 0; i < graph->d_arrays[0].n_dlines; i++) {
	for (j = word_len; j <= max; j++) {
	    score = graph->d_arrays[0].d_array[i].x1 - graph->d_arrays[0].d_array[i].x0;
	    if (score == j) {
		score_hist[j - word_len]++;
		break;
	    }
	}
    }
/*
    for (j = 0; j <= max-word_len ; j++) {
	printf("j %d score %d \n", j, score_hist[j]);
    }
*/

    for (j = max-word_len; j >= 0; j--) {
	cum += score_hist[j];
	score_hist[j] = cum;
    }

    seq_num1 = GetSeqNum(s_result->seq_id[HORIZONTAL]);
    seq_num2 = GetSeqNum(s_result->seq_id[VERTICAL]);

    /* if either seq_num is -1; the corresponding seq must have been deleted */
    if ((seq_num1 == -1) || (seq_num2 == -1)) {
	return 0;
    }
    seq1 = GetSeqSequence(seq_num1);
    seq2 = GetSeqSequence(seq_num2);
    seq1_type = GetSeqType(seq_num1);
    seq2_type = GetSeqType(seq_num2);
    
    if (seq1_type != seq2_type) {
	verror(ERR_FATAL, "calc probs", "sequences must both be either DNA or protein");
	return -1;
    }  else if (seq1_type == PROTEIN) {
	set_char_set(PROTEIN);
	/* set identity matrix */

	if (-1 == (set_matrix_identity(PROTEIN))) {
	    verror(ERR_WARN, "set score matrix", 
		   "unable to set identity score matrix");
	    return TCL_OK;
	}
	set_score_matrix(get_matrix_identity(PROTEIN));
    } else if (seq1_type == DNA) {
	set_char_set(DNA);
	if (-1 == (set_matrix_identity(DNA))) {
	    verror(ERR_WARN, "set score matrix", 
		   "unable to set identity score matrix");
	    return TCL_OK;
	}
	set_score_matrix(get_matrix_identity(DNA));

    }
    ListIdentityProbs(seq1, seq2, graph->dim.x0, graph->dim.x1, graph->dim.y0, 
		      graph->dim.y1, seq1_type, word_len, max, score_hist);

    xfree(score_hist);
    return 0;
}

void identities_shutdown(Tcl_Interp *interp,
			 seq_result *s_result,
			 element *e)
{
    seq_reg_key_name info;
    static char buf[80];

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    /* deregister horizontal sequence */
    seq_deregister(GetSeqNum(s_result->seq_id[HORIZONTAL]), identities_callback,
		   (seq_result *)s_result);

    /* deregister vertical sequence */
    seq_deregister(GetSeqNum(s_result->seq_id[VERTICAL]), identities_callback, 
		   (seq_result *)s_result);

    if (e->num_results > 1) {
	Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	e->replot_func(e);
    } else {
	Tcl_VarEval(e->c->interp, "result_list_update ", NULL);
    }
}

void identities_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_find_identities *input = s_result->input;
    text_find_identities *text = s_result->text_data;
    int result_id = s_result->id;
    char cmd[1024];
    element *e = s_result->e;
    plot_data *result;
    Tcl_Interp *interp;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
    Graph *graph;

    if (e) {
	result = find_plot_data(e, result_id);
	interp = e->c->interp;
    }

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Find matching words");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "matching words #%d", result_id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "matching words: hori=%s vert=%s", 
		GetSeqBaseName(GetSeqNum(s_result->seq_id[HORIZONTAL])),
		GetSeqBaseName(GetSeqNum(s_result->seq_id[VERTICAL])));
	break;
    case SEQ_GET_OPS:
	if (result->hidden) {
	    jdata->get_ops.ops = "Information\0List results\0"
		"Tabulate scores\0PLACEHOLDER\0PLACEHOLDER\0PLACEHOLDER\0"
		"Reveal\0SEPARATOR\0Remove\0";
	} else if ((seq_get_type(result_id) == SEQ_PLOT_TEMP) && !get_replot_temp()) {
	    jdata->get_ops.ops = "Information\0List results\0PLACEHOLDER\0"
		"PLACEHOLDER\0"
	       "Display sequences\0PLACEHOLDER\0PLACEHOLDER\0SEPARATOR\0Remove\0";
	} else {
	    jdata->get_ops.ops = "Information\0List results\0"
		"Tabulate scores\0Configure\0Display sequences\0Hide\0"
		"PLACEHOLDER\0SEPARATOR\0Remove\0";
	}
	break;
    case SEQ_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0: /* information */
	    graph = Tcl_GetGraphFromObj(graph_obj);
	    vfuncheader("input parameters");
	    vmessage("%s Number of matches %d\n", input->params, graph->d_arrays[0].n_dlines);
	    break;
	case 1: /* results */
	    Tcl_Eval(interp, "SetBusy");
	    vfuncheader("results");
	    s_result->txt_func(s_result);
	    Tcl_Eval(interp, "ClearBusy");
	    break;
	case 2: /* scores */
	    {
		Tcl_Eval(interp, "SetBusy");
		vfuncheader("scores");
		CalcIdentityProbs(s_result, text->word_len);
		Tcl_Eval(interp, "ClearBusy");
		break;
	    }
	case 3: /* configure */
	    sprintf(cmd, "result_config %d %d %s %d %s", 
		    result_id, result->line_width, result->colour, e->id, 
		    e->win);
	    
	    if (TCL_OK != Tcl_Eval(interp, cmd)){
		puts(interp->result);
	    }
	    break;
	case 4: /* display sequences */
	    SequencePairDisplay(interp, e->win, result_id, 
				s_result->seq_id[HORIZONTAL], 
				s_result->seq_id[VERTICAL]);
	    break;
	case 5: /* hide all */
	    result->hidden = 1;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 6: /* reveal all */
	    result->hidden = 0;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 7: /* remove */
	    {
		identities_shutdown(interp, s_result, e);
		remove_result_from_element(e, result_id);
		DestroySequencePairDisplay(interp, result_id);
		free(input->params);
		xfree(s_result->text_data);
		break;
	    }
	}
	break;
    case SEQ_PLOT: 
	{
	    s_result->pr_func(s_result, NULL);
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
	    identities_shutdown(interp, s_result, e);
	    DestroySequencePairDisplay(interp, result_id);
	    free(input->params);
	    xfree(s_result->text_data);
	    break;
	}
    }
}


void identities_text_func(void *obj)
{
    seq_result *s_result = (seq_result *) obj;
    int i;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    char *tmp_str;
    int max_str;
    int seq_num_h, seq_num_v;
    Graph *graph;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
    int score;
     
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
	
    if (NULL == (tmp_str = (char *)xmalloc(max_str * sizeof(char))))
	return;

    for (i = 0; i < graph->d_arrays[0].n_dlines; i++) {
	UpdateTextOutput();
	score = graph->d_arrays[0].d_array[i].x1 - graph->d_arrays[0].d_array[i].x0;
	vmessage("Positions %10d h %10d v and length %10d\n", 
		 (int)graph->d_arrays[0].d_array[i].x0, (int)graph->d_arrays[0].d_array[i].y0, 
		 score);
	strncpy(tmp_str, &seq1[(int)graph->d_arrays[0].d_array[i].x0 - 1], 
			      score);
	tmp_str[score] = '\0';
	vmessage("%s\n", tmp_str);
    }
    xfree(tmp_str);
}

void identities_recalc_func(void *obj, seq_reg_plot *plot)
{
}

int store_matching_words(int seq1_num,
			 int seq2_num,
			 int start_h,
			 int end_h,
			 int start_v,
			 int end_v,
			 int word_length,
			 in_find_identities *input,
			 int *seq1_match, 
			 int *seq2_match, 
			 int *match_score,
			 int num_elements,
			 int strand,
			 int too_many_matches,
			 Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    int i, id;
    text_find_identities *text_data;
    Graph *graph;

    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;
    
    Tcl_InitGraph(&graph);

    if (NULL == (text_data = (text_find_identities *)xmalloc(sizeof(text_find_identities))))
	return -1;


    if (!too_many_matches) {
	if (NULL == (graph->d_arrays = (darray *)xmalloc(sizeof(darray))))
	    return -1;
    
	if (NULL == (graph->d_arrays[0].d_array = (gd_line *)xmalloc(num_elements * sizeof(gd_line))))
	    return -1;
    
	/* data assignment */
	for (i = 0; i < num_elements; i++) {
	    graph->d_arrays[0].d_array[i].x0 = seq1_match[i];	
	    graph->d_arrays[0].d_array[i].x1 = seq1_match[i] + match_score[i];

	    if (strand == TOP_S) {
		graph->d_arrays[0].d_array[i].y0 = seq2_match[i];
		graph->d_arrays[0].d_array[i].y1 = seq2_match[i] + match_score[i];
	    } else {
		graph->d_arrays[0].d_array[i].y0 = seq2_match[i] + match_score[i];
		graph->d_arrays[0].d_array[i].y1 = seq2_match[i];
	    }
	}


	/* array must be in ascending x order */
	qsort((void *) graph->d_arrays[0].d_array, num_elements, sizeof(gd_line), 
	      compare_gd_line);
	
    }

    id = get_reg_id();
    
    graph->d_arrays[0].n_dlines = num_elements;
    graph->d_arrays[0].type = G_LINES;
    graph->n_darrays = 1;
    graph->dim.x0 = start_h;
    graph->dim.x1 = end_h;
    graph->dim.y0 = start_v;
    graph->dim.y1 = end_v;
    *graph_obj = Tcl_NewGraphObj(graph);

    /* these are copied over in Tcl_NewGraphObj so don't need them anymore */
    xfree(graph->d_arrays[0].d_array);
    xfree(graph->d_arrays);
    xfree(graph);

    s_result->text_data = text_data;
    text_data->word_len = word_length;

    s_result->seq_id[HORIZONTAL] = GetSeqId( seq1_num);
    s_result->seq_id[VERTICAL] = GetSeqId(seq2_num);

    s_result->input = (void *)input; 
    s_result->id = id; 
    s_result->type = SEQ_TYPE_DOT_PLOT;
    s_result->gr_type = SEQ_TYPE_DOT_PLOT;
    s_result->frame = 0;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */

    s_result->op_func = identities_callback;
    s_result->txt_func = identities_text_func;

    /* too_many_matches is always false at the present time */
    if (too_many_matches) {
	s_result->pr_func = identities_recalc_func;
	seq_register(seq1_num, identities_callback, (void *)s_result, 
		     SEQ_PLOT_TEMP, id);
	seq_register(seq2_num, identities_callback, (void *)s_result, 
		     SEQ_PLOT_TEMP, id);
    } else {
	s_result->pr_func = seq_plot_graph_func;
	seq_register(seq1_num, identities_callback, (void *)s_result, 
		     SEQ_PLOT_PERM, id);
	seq_register(seq2_num, identities_callback, (void *)s_result, 
		     SEQ_PLOT_PERM, id);
    }

    return id;
}

int init_sip_matching_words_create(Tcl_Interp *interp, 
				   int seq_id_h,
				   int seq_id_v, 
				   int start_h,
				   int end_h, 
				   int start_v,
				   int end_v, 
				   int word_len,
				   int strand,
				   Tcl_Obj **graph_obj,
				   int *id)
{
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int max_matches;
    int n_matches;
    int *seq1_match = NULL;
    int *seq2_match = NULL; 
    int *len_match = NULL;
    int max_hash, word_length;
    in_find_identities *input = NULL;
    int same_seq;
    int too_many_matches = 0;
    int seq1_num, seq2_num;
    int seq1_type, seq2_type;
    int sub1_len, sub2_len;
    double prob;
    Tcl_DString input_params;    
    char strand_s[7];

    vfuncheader("find matching words");

    if (-1 == (seq1_num = GetSeqNum(seq_id_h))) {
	verror(ERR_WARN, "find matching words", 
	       "horizontal sequence undefined");
	goto error;
    }
     
    if (-1 == (seq2_num = GetSeqNum(seq_id_v))) {
	verror(ERR_WARN, "find matching words", 
	       "vertical sequence undefined");
	goto error;
    }

    if (NULL == (input = (in_find_identities *)
		          xmalloc(sizeof(in_find_identities))))
	goto error;

    seq1 = GetSeqSequence(seq1_num);
    seq1_len = GetSeqLength(seq1_num);
    seq2 = GetSeqSequence(seq2_num);
    seq2_len = GetSeqLength(seq2_num);

    /* returns 1 for dna, 2 for protein, 0 for anything else */
    seq1_type = GetSeqType(seq1_num);
    seq2_type = GetSeqType(seq2_num);
    
    if (end_h == -1)
	end_h = seq1_len;

    if (end_v == -1)
	end_v = seq2_len;
    
    if (seq1_type != seq2_type) {
	verror(ERR_WARN, "find matching words", "sequences must both be either DNA or protein");
	goto error;
    } else if (seq1_type == PROTEIN) {
	set_char_set(PROTEIN);
	/* set identity matrix */

	if (-1 == (set_matrix_identity(PROTEIN))) {
	    verror(ERR_WARN, "set score matrix", 
		   "unable to set identity score matrix");
	    goto error;
	}
	set_score_matrix(get_matrix_identity(PROTEIN));
	max_hash = 3;

    } else if (seq1_type == DNA) {
	set_char_set(DNA);
	if (-1 == (set_matrix_identity(DNA))) {
	    verror(ERR_WARN, "set score matrix", 
		   "unable to set identity score matrix");
	    goto error;
	}
	set_score_matrix(get_matrix_identity(DNA));

	max_hash = 7;
    }

    word_length = MIN(word_len, max_hash);

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

    /* complement if necessary */
    if (strand == BOTTOM_S) {
	if (seq2_type == DNA)
	    complement_seq(seq2, seq2_len);
    }
   
    /* 
     * decide whether to save matches or not depending on number of expected
     * matches 
     */
    prob = FindExpectedProb(seq1, seq2, start_h, end_h, start_v,
			    end_v, word_len, seq1_type);
    max_matches = get_max_matches();

    /* 
     * next bit has been disabled for consistency with other searches and
     * therefore I didn't bother putting it here. Must look at backups if ever
     * want to reintroduce permanent and temporary matches
     */

    /* 
     * allocate match arrays here if saving the results.
     */
    
    if (prob > max_matches) {
	max_matches = (int) prob;
    } else {
	max_matches = get_max_matches();
    }
    if (NULL == (seq1_match = (int *)xmalloc(max_matches * sizeof(int))))
	goto error;
    if (NULL == (seq2_match = (int *)xmalloc(max_matches * sizeof(int))))
	goto error;
    if (NULL == (len_match = (int *)xmalloc(max_matches * sizeof(int))))
	goto error;
    
    sip_hash(seq1, seq2, start_h, end_h, start_v, end_v, max_matches, 
	     word_len, word_length, seq1_type, same_seq, &seq1_match, 
	     &seq2_match, &len_match, &n_matches, NULL, NULL);

    /* complement back again */
    if (strand == BOTTOM_S) {
	if (seq2_type == DNA)
	    complement_seq(seq2, seq2_len);
    }

    if (n_matches < 0) {
	verror(ERR_WARN, "find matching words", 
	       "failed in find matching words \n");
	goto error;
    } else if (n_matches == 0) {
	verror(ERR_WARN, "find matching words", "no matches found \n");
	if (seq1_match)
	    xfree(seq1_match);
	if (seq2_match) 
	    xfree(seq2_match);
	if (len_match)
	    xfree(len_match);
	if (input)
	    xfree(input);
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
    vTcl_DStringAppend(&input_params, "horizontal %s: %s\nvertical %s: %s\n"
	    "strand %s word length %d \n", 
	    GetSeqLibraryName(seq1_num), 
	    GetSeqName(seq1_num), 
	    GetSeqLibraryName(seq2_num),
	    GetSeqName(seq2_num), strand_s, word_len);
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);
    
    if (-1 == (*id = store_matching_words(seq1_num, seq2_num,start_h,
					  end_h, start_v, end_v,word_length,
					  input, seq1_match, seq2_match, 
					  len_match, n_matches, strand,
					  too_many_matches, graph_obj))) {
	goto error;
    }
    if (seq1_match)
	xfree(seq1_match);
    if (seq2_match) 
	xfree(seq2_match);
    if (len_match)
	xfree(len_match);
    return 0;

 error:
    verror(ERR_WARN, "find matching words", "failure find matching words");
    if (seq1_match)
	xfree(seq1_match);
    if (seq2_match) 
	xfree(seq2_match);
    if (len_match)
	xfree(len_match);
    if (input)
	xfree(input);
    return -1;

}

int init_sip_matching_words_plot(Tcl_Interp *interp,
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
