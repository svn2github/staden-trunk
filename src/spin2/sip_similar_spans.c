#include <string.h>

#include "seq_reg.h"
#include "seq_results.h"
#include "compare_spans.h"
#include "sip_hash.h"
#include "sip_results.h"
#include "seq_reg.h"
#include "readpam.h"
#include "misc.h" /* MIN */
#include "text_output.h"  /* UpdateTextOutput */
#include "dna_utils.h"
#include "tcl_utils.h"
#include "sip_globals.h"
#include "rescan_matches.h"
#include "sip_similar_spans.h"
#include "sequence_formats.h"
#include "seq_plot_funcs.h"
#include "sequence_pair_display.h"
#include "tclCanvGraph.h"
#include "seq_element.h"
#include "probs.h"

void similar_spans_callback(int seq_num, void *obj, seq_reg_data *jdata);

/************************************************************/
int
Compare_Spans(char *seq_1,                                             /* in */
	      char *seq_2,                                             /* in */
	      int seq1_len,                                            /* in */
	      int seq2_len,                                            /* in */
	      int seq1_start,
	      int seq1_end,
	      int seq2_start,
	      int seq2_end,
	      int max_matches,                                         /* in */
	      int same_seq,                                            /* in */
	      int window_length,                                       /* in */
	      int min_match,                                           /* in */
	      int fopt,                                                /* in */
	      int ropt,                                                /* in */
	      int **seq1_match,                                       /* out */
	      int **seq2_match,                                       /* out */
	      int **match_score,                                      /* out */
	      int *n_matches)                                         /* out */
{
    char dir;
    int num_files;

    num_files = 2;

    if ( 0 != num_files ) {
	if ( fopt ) {
	    dir = 'f';
	    *n_matches = cmpspn (&dir,
				 &min_match,
				 seq1_match,
				 seq2_match,
				 match_score,
				 &max_matches,
				 &window_length,
				 seq_1,
				 seq_2,
				 &seq1_len,
				 &seq2_len,
				 seq1_start,
				 seq1_end,
				 seq2_start,
				 seq2_end,
				 same_seq);

/*
	    for (i = 0; i < *n_matches; i++) {
		printf ( "%d %d %d\n",seq1_match[i], seq2_match[i], match_score[i]);
	    }
*/
	}

	if ( ropt ) {

	    dir = 'r';
	    *n_matches = cmpspn (
				 &dir,
				 &min_match,
				 seq1_match,
				 seq2_match,
				 match_score,
				 &max_matches,
				 &window_length,
				 seq_1,
				 seq_2,
				 &seq1_len,
				 &seq2_len,
				 seq1_start,
				 seq1_end,
				 seq2_start,
				 seq2_end,
				 same_seq
				 );

/*
	    for (i = 0; i< n_matches; i++) {
		printf ( "%d %d %d\n",seq1_match[i], seq2_match[i], match_score[i]);
	    }
*/
	}
	return 0;
    }
    return -1;
}
/*
 * score, probabilites, expected and observed scores output in the
 * text output window
 */
int
CalcProbs(seq_result *s_result, 
	  int *match_score,
	  int span_length, 
	  int min_score)
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

    graph = Tcl_GetGraphFromObj(graph_obj);

    for (i = 0; i < graph->p_arrays[0].n_pts; i++) {
#ifdef DEBUG
	printf("i %d score %d \n", i, match_score[i]);
#endif
	if (match_score[i] > max) 
	    max = match_score[i];
    }
#ifdef DEBUG
    printf("min %d max %d \n", min_score, max);
#endif
    if (NULL == (score_hist = (int *)xcalloc(max - min_score + 1, sizeof(int))))
	return -1;
    
    /* create histogram of scores between min_score and max */
    for (i = 0; i < graph->p_arrays[0].n_pts; i++) {
	for (j = min_score; j <= max; j++) {
	    if (match_score[i] == j) {
		score_hist[j - min_score]++;
		break;
	    }
	}
    }
/*
    for (j = 0; j <= max-min_score ; j++) {
	printf("j %d score %d \n", j, score_hist[j]);
    }
*/

    for (j = max-min_score; j >= 0; j--) {
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
	return -1;;
    }  else if (seq1_type == PROTEIN) {
	set_char_set(PROTEIN);
        set_score_matrix(get_matrix_file(PROTEIN));
    } else if (seq1_type == DNA) {
	set_char_set(DNA);
        set_score_matrix(get_matrix_file(DNA));
    }
    ListProbs(seq1, seq2, graph->dim.x0, graph->dim.x1, 
	      graph->dim.y0, graph->dim.y1, span_length, seq1_type, 
	      min_score, max, score_hist);

    xfree(score_hist);
    return 0;
}


void similar_spans_shutdown(Tcl_Interp *interp,
			    seq_result *s_result,
			    element *e)
{
    seq_reg_key_name info;
    static char buf[80];
    in_sim_spans *input = s_result->input;
    text_sim_spans *text_data = s_result->text_data;

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    /* deregister horizontal sequence */
    seq_deregister(GetSeqNum(s_result->seq_id[HORIZONTAL]), 
		   similar_spans_callback, (seq_result *)s_result);
    
    /* deregister vertical sequence */
    seq_deregister(GetSeqNum(s_result->seq_id[VERTICAL]), 
		   similar_spans_callback, (seq_result *)s_result);
    
    if (e->num_results > 1) {
	Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	e->replot_func(e);
    } else {
	Tcl_VarEval(e->c->interp, "result_list_update ", NULL);
    }

    DestroySequencePairDisplay(interp, s_result->id);
    free(input->params);
    xfree(text_data->match_score);
    xfree(s_result->text_data);
}

void similar_spans_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_sim_spans *input = s_result->input;
    text_sim_spans *text_data = s_result->text_data;
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
	sprintf(jdata->name.line, "Find similar spans");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "similar spans #%d", result_id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "similar spans: hori=%s vert=%s", 
		GetSeqBaseName(GetSeqNum(s_result->seq_id[HORIZONTAL])),
		GetSeqBaseName(GetSeqNum(s_result->seq_id[VERTICAL])));
	break;
    case SEQ_GET_OPS:
	if (result->hidden) {
	    jdata->get_ops.ops = "Information\0List results\0Tabulate scores\0PLACEHOLDER\0"
		"PLACEHOLDER\0PLACEHOLDER\0Reveal\0SEPARATOR\0Remove\0";
	} else {
	    jdata->get_ops.ops = "Information\0List results\0Tabulate scores\0Configure\0"
	       "Display sequences\0Hide\0PLACEHOLDER\0SEPARATOR\0Remove\0";
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
	case 2: /* scores */
	    Tcl_Eval(interp, "SetBusy");
	    vfuncheader("scores");
	    CalcProbs(s_result, text_data->match_score, input->win_length, 
		      text_data->min_score);
	    Tcl_Eval(interp, "ClearBusy");
	    break;
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
		similar_spans_shutdown(interp, s_result, e);
		remove_result_from_element(e, result_id);
		break;
	    }
	}
	break;
    case SEQ_PLOT:
	s_result->pr_func(s_result, NULL);
	break;
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
	    similar_spans_shutdown(interp, s_result, e);
	    break;
	}
    }
}

void similar_spans_text_func(void *obj)
{
    seq_result *s_result = (seq_result *) obj;
    in_sim_spans *input = s_result->input;
    int i;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int max_str;
    char *r_seq1, *r_seq2;
    int r;
    int num_spaces;
    int seq_num_h, seq_num_v;
    int seq1_type;
    text_sim_spans *text_data = s_result->text_data;
    Graph *graph;
    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
    double mid_pt, x, y;

    mid_pt = input->win_length/2;

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

    if (seq1_len >= input->win_length) {
	if (NULL == (r_seq1 = (char *)xmalloc(seq1_len + 1)))
	    return;
    } else {
	if (NULL == (r_seq1 = (char *)xmalloc(input->win_length + 1)))
	    return;
    }
    if (seq2_len >= input->win_length) {
	if (NULL == (r_seq2 = (char *)xmalloc(seq2_len + 1)))
	    return;
    } else {
	if (NULL == (r_seq2 = (char *)xmalloc(input->win_length + 1)))
	    return;
    }

    for (i = 0; i < graph->p_arrays[0].n_pts; i++) {
	x = graph->p_arrays[0].p_array[i].x - mid_pt;
	y = graph->p_arrays[0].p_array[i].y - mid_pt;

	UpdateTextOutput();
	vmessage("Positions %10d h %10d v and score %10d\n", 
	       (int)x, (int)y, text_data->match_score[i]);	
	if (x <= 0) {
	    num_spaces = abs(x) + 1;
	    memset(r_seq1, ' ', num_spaces);
	    r_seq1[num_spaces] = '\0';
	    strncat(r_seq1, seq1, input->win_length - num_spaces);
	} 
	if (y <= 0) {
	    num_spaces = abs(y) + 1;
	    memset(r_seq2, ' ', num_spaces);
	    r_seq2[num_spaces] = '\0';
	    strncat(r_seq2, seq2, input->win_length - num_spaces);
	} 
	if (x > 0) {
	    strncpy(r_seq1, &seq1[(int)x-1], input->win_length);
	}
	if (y > 0) {
	    strncpy(r_seq2, &seq2[(int)y-1], input->win_length);
	}
	r_seq1[input->win_length] = '\0';
	r_seq2[input->win_length] = '\0';
	/*yy*/
	seq1_type = GetSeqType(seq_num_h);
	r = spin_list_alignment(r_seq1, r_seq2, "H", "V", x, y, "", seq1_type);
#ifdef DEBUG
	printf("i %d \n%s\n%s\n", i, r_seq1, r_seq2);
#endif
	r_seq1[0] = '\0';
	r_seq2[0] = '\0';
    }
    xfree(r_seq1);
    xfree(r_seq2);
}

int store_sip_similar_spans(int seq1_num,
			    int seq2_num,
			    char *seq1, 
			    char *seq2,
			    int seq1_len,
			    int seq2_len,
			    int window_length,
			    int min_score,
			    int start_h,
			    int end_h,
			    int start_v,
			    int end_v,
			    int *seq1_match, 
			    int *seq2_match, 
			    int *match_score,
			    int num_elements,
			    int strand,
			    int min_char_score,
			    in_sim_spans *input,
			    Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    text_sim_spans *text_data;
    int i, j, id;
    Graph *graph;
    double mid_pt;
    double x, y;
    int score;
    int cnt = 0;

    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;
    
    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;
    
    Tcl_InitGraph(&graph);

    if (NULL == (graph->p_arrays = (parray *)xmalloc(sizeof(parray))))
	return -1;
    
    if (min_char_score == -1) {
	if (NULL == (graph->p_arrays[0].p_array = (g_pt *)xmalloc(num_elements * sizeof(g_pt))))
	    return -1;

    } else {
	if (NULL == (graph->p_arrays[0].p_array = (g_pt *)xmalloc(num_elements * window_length * sizeof(g_pt))))
	    return -1;
    }	

    if (NULL == (text_data = (text_sim_spans *)xmalloc(sizeof(text_sim_spans))))
	return -1;

    id = get_reg_id();
    mid_pt = window_length/2;

     /* data assignment */
    for (i = 0; i < num_elements; i++) {
	x = seq1_match[i];
	y = seq2_match[i];

	if (min_char_score == -1) {
	    graph->p_arrays[0].p_array[i].x = x + mid_pt;
	    graph->p_arrays[0].p_array[i].y = y + mid_pt;
	} else {
	    for (j = 0; j < window_length; j++, x++, y++) {
		
		if ((x < 1) || (y < 1) || (x > start_h + seq1_len) || (y > start_v + seq2_len)) {
		    continue;
		}
		score = score_matrix[char_lookup[seq1[(int)x - 1]]]
		    [char_lookup[seq2[(int)y - 1]]];
		
		if (score >= min_char_score) {
		    graph->p_arrays[0].p_array[cnt].x = x;
		    graph->p_arrays[0].p_array[cnt].y = y;
		    cnt++;
		}
	    }
	}
    }

    if (min_char_score != -1) 
	num_elements = cnt;

    /* array must be in ascending x order */
    qsort((void *) graph->p_arrays[0].p_array, num_elements, sizeof(g_pt), 
	  compare_g_pt);
    
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

    s_result->text_data = text_data;
    text_data->min_score = min_score;
    text_data->match_score = match_score;
    text_data->num = num_elements;

    /* save global sip_seq into sip_result structure */
    s_result->seq_id[HORIZONTAL] = GetSeqId(seq1_num);
    s_result->seq_id[VERTICAL] = GetSeqId(seq2_num);

    s_result->input = (void *)input; 
    s_result->id = id; 
    s_result->type = SEQ_TYPE_DOT_PLOT;
    s_result->gr_type = SEQ_TYPE_DOT_PLOT;
    s_result->frame = 0;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */

    /* s_result->pr_func = dot_plot_middot_func; */
    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = similar_spans_callback;
    s_result->txt_func = similar_spans_text_func;
    
    seq_register(seq1_num, similar_spans_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);
    seq_register(seq2_num, similar_spans_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);

    return id;
}

int init_sip_similar_spans_create(Tcl_Interp *interp, 
				  int seq_id_h,
				  int seq_id_v, 
				  int start_h,
				  int end_h, 
				  int start_v,
				  int end_v, 
				  int win_len,
				  int min_match, 
				  int strand,
				  int score,
				  Tcl_Obj **graph_obj,
				  int *id)
{
    in_sim_spans *input = NULL;
    int *seq1_match = NULL;
    int *seq2_match = NULL;
    int *match_score = NULL;
    int n_matches;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int same_seq;
    int max_matches = get_max_matches();
    int seq1_num, seq2_num;
    int seq1_type, seq2_type;
    int sub1_len, sub2_len;
    Tcl_DString input_params;
    char strand_s[7];

    vfuncheader("find similar spans");
    
    if (NULL == (seq1_match = (int *)xmalloc(max_matches * sizeof(int))))
	goto error;
    if (NULL == (seq2_match = (int *)xmalloc(max_matches * sizeof(int))))
	goto error;
    if (NULL == (match_score = (int *)xmalloc(max_matches * sizeof(int))))
	goto error;
    if (NULL == (input = (in_sim_spans *)xmalloc(sizeof(in_sim_spans))))
	goto error;
    
    /* get first and second sequence saved using extract_sequence */
    seq1_num = GetSeqNum(seq_id_h);
    seq2_num = GetSeqNum(seq_id_v);
    
    if (seq1_num == -1) {
	verror(ERR_WARN, "find similar spans", "horizontal sequence undefined");
	goto error;
    } else if (seq2_num == -1) {
	verror(ERR_WARN, "find similar spans", "vertical sequence undefined");
	goto error;
    }

    seq1 = GetSeqSequence(seq1_num);
    seq2 = GetSeqSequence(seq2_num);
    seq1_len = GetSeqLength(seq1_num);
    seq2_len = GetSeqLength(seq2_num);
    seq1_type = GetSeqType(seq1_num);
    seq2_type = GetSeqType(seq2_num);

    if (end_h == -1)
	end_h = seq1_len;

    if (end_v == -1)
	end_v = seq2_len;

    if (seq1_type != seq2_type) {
	verror(ERR_WARN, "find similar spans", "sequences must both be either DNA or protein");
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

    /*
     * Should check length of sub sequences only. These lengths are not
     * stored, so have to calculate them here. Not storing them in
     * seq1_len and seq2_len as I'm unsure whether subsequent functions
     * expect the length of the whole sequence. Anyway, the compare_spans
     * function recalculates the lengths of the sub sequences before doing
     * the comparison.
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

    if (strand == TOP_S) {
	Compare_Spans(seq1, seq2, seq1_len, seq2_len, start_h, end_h, 
		      start_v, end_v, max_matches, same_seq, 
		      win_len, min_match, 1, 0,
		      &seq1_match, &seq2_match, &match_score, &n_matches);
    } else {
	if (seq2_type == DNA)
	    complement_seq(seq2, seq2_len);

	Compare_Spans(seq1, seq2, seq1_len, seq2_len, start_h, end_h, 
		      start_v, end_v, max_matches, same_seq, 
		      win_len, min_match, 1, 0,
		      &seq1_match, &seq2_match, &match_score, &n_matches);

	if (seq2_type == DNA)
	    complement_seq(seq2, seq2_len);
    }


    /* n_matches == -1 if malloc problem or -2 if too many matches */
    if (n_matches == -2) {
	verror(ERR_WARN, "find similar spans", "too many matches");
	goto error;
    } else if (n_matches == -1) {
	goto error;
    } else if (n_matches == 0) {
	verror(ERR_WARN, "Find similar spans", "no matches found\n"); 
	if (seq1_match)
	    xfree (seq1_match);
	if (seq2_match)
	    xfree (seq2_match);
	if (match_score)
	    xfree(match_score);
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

    vTcl_DStringAppend(&input_params, "horizontal %s: %s \nvertical %s: %s\n"
	    "strand %s window length %d min match %d number of matches %d", 
	    GetSeqLibraryName(seq1_num), 
	    GetSeqName(seq1_num), 
	    GetSeqLibraryName(seq2_num), 
	    GetSeqName(seq2_num), 
	    strand_s, win_len, min_match, n_matches);
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    input->win_length = win_len;
    if (-1 == (*id = store_sip_similar_spans(seq1_num, seq2_num, seq1, seq2,
					     sub1_len, sub2_len, win_len,
					     min_match, start_h, end_h, 
					     start_v, end_v,
					     seq1_match, seq2_match, 
					     match_score, n_matches,
					     strand, score, input, graph_obj))) {
	goto error;
    }
    
    if (seq1_match)
	xfree (seq1_match);
    if (seq2_match)
	xfree (seq2_match);

    return 0;
    
 error:
    verror(ERR_WARN, "find similar spans", "failure in find similar spans");
    if (seq1_match)
	xfree (seq1_match);
    if (seq2_match)
	xfree (seq2_match);
    if (match_score)
	xfree(match_score);
    if (input)
      xfree(input);
    return -1;
}

int init_sip_similar_spans_plot(Tcl_Interp *interp,
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
