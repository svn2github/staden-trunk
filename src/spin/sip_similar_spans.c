#include <string.h>

#include "seq_reg.h"
#include "seq_results.h"
#include "seq_raster.h"
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

void similar_spans_shutdown(Tcl_Interp *interp,
			    seq_result *result,
			    char *raster_win,
			    int seq_num,
			    int id)
{
    char *tmp;
    seq_reg_key_name info;
    static char buf[80];
    int raster_id;
    double wx0, wy0, wx1, wy1;
    Tk_Raster *raster;
    Tcl_CmdInfo info1;
    RasterResult *raster_result;
    in_comp_spans *input = result->input;

    /* determine raster_id and raster_result structure */
    Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));
    raster_result = raster_id_to_result(raster_id);

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(result->id, (seq_reg_data *)&info, 0);

    /* deregister horizontal sequence */
    seq_deregister(GetSeqNum(result->seq_id[HORIZONTAL]), 
		   similar_spans_callback, (seq_result *)result);
    
    /* deregister vertical sequence */
    seq_deregister(GetSeqNum(result->seq_id[VERTICAL]), 
		   similar_spans_callback, (seq_result *)result);
    
    /* 
     * only bother replotting the raster if there are still results in the
     * raster
     */
    if (raster_result && raster_result->num_results > 1) {
	
	tmp = get_default_string(interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));

	ReplotAllCurrentZoom(interp, raster_win);

	Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
	raster_id = atoi(Tcl_GetStringResult(interp));
	if (TCL_OK != Tcl_VarEval(interp, "RemoveRasterResultKey ", raster_win,
				  " {", info.line, "}", NULL))
	    verror(ERR_WARN, "similar spans_shutdown1", "%s \n", 
		   Tcl_GetStringResult(interp));
	
	/* find original y before reset size */
	Tcl_GetCommandInfo(interp, raster_win, &info1);
	raster = (Tk_Raster*)info1.clientData;
	RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

	/* update the scale of the remaining results in the raster window */
	SeqReSetRasterWindowSize(interp, raster_win, result->graph);
	ReSetRasterWindowWorld(interp, raster_win, wy1, result->graph);
	ReplotAllRasterWindow(interp, raster_win);

	if (TCL_OK != Tcl_VarEval(interp, "seq_result_list_update ", tmp, NULL)){
	    verror(ERR_WARN, "similar_spans_shutdown2", "%s\n", Tcl_GetStringResult(interp));
	}
    }
    DestroySequencePairDisplay(interp, id);
    free(input->params);
    xfree(result->text_data);

    SipFreeResult(result);

    if (raster_result) 
	DeleteResultFromRaster(raster_result);
}

void similar_spans_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_comp_spans *input = result->input;
    out_raster *output = result->output;
    d_plot *data = result->data;
    text_sim_spans *text = result->text_data;

    int id = result->id;
    char cmd[1024];

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Find similar spans");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "similar spans #%d", result->id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "similar spans: hori=%s vert=%s", 
		GetSeqBaseName(GetSeqNum(result->seq_id[HORIZONTAL])),
		GetSeqBaseName(GetSeqNum(result->seq_id[VERTICAL])));
	break;
    case SEQ_GET_OPS:
	if (output->hidden) {
	    jdata->get_ops.ops = "Information\0List results\0Tabulate scores\0PLACEHOLDER\0PLACEHOLDER\0"
		"PLACEHOLDER\0PLACEHOLDER\0Reveal\0SEPARATOR\0Remove\0";
	} else {
	    jdata->get_ops.ops = "Information\0List results\0Tabulate scores\0Rescan matches\0Configure\0"
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
	    Tcl_Eval(output->interp, "SetBusy");
	    vfuncheader("results");
	    result->txt_func(result);
	    Tcl_Eval(output->interp, "ClearBusy");
	    break;
	case 2: /* scores */
	    Tcl_Eval(output->interp, "SetBusy");
	    vfuncheader("scores");
	    CalcProbs(result, data->win_len, text->min_score);
	    Tcl_Eval(output->interp, "ClearBusy");
	    break;
	case 3: /* rescan matches */
	    {
		Tcl_Eval(output->interp, "sip_rescan_matches");
		Tcl_Eval(output->interp, "SetBusy");
		SipRescanMatches(output->interp, result, id, atoi(Tcl_GetStringResult(output->interp)));
		Tcl_Eval(output->interp, "ClearBusy");
		break;
	    }
	case 4: /* configure */
	    sprintf(cmd, "RasterConfig %d", id);
	    if (TCL_OK != Tcl_Eval(output->interp, cmd)){
		puts(Tcl_GetStringResult(output->interp));
	    }
	    break;
	case 5: /* display sequences */

	    SequencePairDisplay(output->interp, output->raster_win, id, 
				result->seq_id[HORIZONTAL], 
				result->seq_id[VERTICAL]);
	    break;
	case 6: /* hide all */
	    output->hidden = 1;
	    ReplotAllCurrentZoom(output->interp, output->raster_win);
	    break;
	case 7: /* reveal all */
	    output->hidden = 0;
	    ReplotAllCurrentZoom(output->interp, output->raster_win);
	    break;
	case 8: /* remove */ 
	    {
		Tcl_Interp *interp = output->interp;
		similar_spans_shutdown(interp, result, output->raster_win, 
				       seq_num, id);
		break;
	    }
	}
	break;
    case SEQ_PLOT:
	result->pr_func(result, NULL);
	break;
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case OUTPUT:
	    jdata->info.result = (void *)output;
	    break;
	case INPUT: 
	    jdata->info.result = (void *)input;
	    break;
	case DIMENSIONS: 
	    jdata->info.result = (void *)&data->dim;
	    break;
	case INDEX: 
	    jdata->info.result = (void *)id;
	    break;
	case RESULT:
	    jdata->info.result = (void *)result;
	    break;
	case WIN_NAME:
	    {
		char *r_win = output->raster_win;
		jdata->info.result = (void *)r_win;
		break;
	    }
	case WIN_SIZE:
	    {
		static d_point pt;
		Tcl_Interp *interp = output->interp;
		pt.x = get_default_int(interp, sip_defs, 
					w("RASTER.PLOT_WIDTH"));
		pt.y = get_default_double(interp, sip_defs,
					   w("RASTER.PLOT_HEIGHT"));

		jdata->info.result = (void *)&pt;
		break;
	    } /* WIN_SIZE */
	}
	break;
    case SEQ_HIDE: 
	output->hidden = 1;
	break;
    case SEQ_REVEAL: 
	output->hidden = 0;
	break;
    case SEQ_QUIT: 
    case SEQ_DELETE: 
	{
	    Tcl_Interp *interp = output->interp;
	    similar_spans_shutdown(interp, result, output->raster_win, 
				   seq_num, id);
	    break;
	}
    }
}

void similar_spans_text_func(void *obj)
{
    seq_result *result = (seq_result *) obj;
    d_plot *data = result->data;
    int n_pts = data->n_pts;
    int i;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int max_str;
    char *r_seq1, *r_seq2;
    int r;
    int num_spaces;
    int seq_num_h, seq_num_v;
    int seq1_type;

    seq_num_h = GetSeqNum(result->seq_id[HORIZONTAL]);
    seq_num_v = GetSeqNum(result->seq_id[VERTICAL]);

    seq1 = GetSeqSequence(seq_num_h);
    seq1_len = GetSeqLength(seq_num_h);
    seq2 = GetSeqSequence(seq_num_v);
    seq2_len = GetSeqLength(seq_num_v);

    if (seq1_len > seq2_len)
	max_str = seq1_len;
    else
	max_str = seq2_len;

    if (seq1_len >= data->win_len) {
	if (NULL == (r_seq1 = (char *)xmalloc(seq1_len + 1)))
	    return;
    } else {
	if (NULL == (r_seq1 = (char *)xmalloc(data->win_len + 1)))
	    return;
    }
    if (seq2_len >= data->win_len) {
	if (NULL == (r_seq2 = (char *)xmalloc(seq2_len + 1)))
	    return;
    } else {
	if (NULL == (r_seq2 = (char *)xmalloc(data->win_len + 1)))
	    return;
    }

    for (i = 0; i < n_pts; i++) {
	UpdateTextOutput();
	vmessage("Positions %10d h %10d v and score %10d\n", 
	       data->p_array[i].x, (int)data->p_array[i].y, 
		 data->p_array[i].score);	
	if (data->p_array[i].x <= 0) {
	    num_spaces = abs(data->p_array[i].x) + 1;
	    memset(r_seq1, ' ', num_spaces);
	    r_seq1[num_spaces] = '\0';
	    strncat(r_seq1, seq1, data->win_len - num_spaces);
	} 
	if (data->p_array[i].y <= 0) {
	    num_spaces = abs(data->p_array[i].y) + 1;
	    memset(r_seq2, ' ', num_spaces);
	    r_seq2[num_spaces] = '\0';
	    strncat(r_seq2, seq2, data->win_len - num_spaces);
	} 
	if (data->p_array[i].x > 0) {
	    strncpy(r_seq1, &seq1[data->p_array[i].x-1], data->win_len);
	}
	if (data->p_array[i].y > 0) {
	    strncpy(r_seq2, &seq2[(int)data->p_array[i].y-1], data->win_len);
	}
	r_seq1[data->win_len] = '\0';
	r_seq2[data->win_len] = '\0';
	/*yy*/
	seq1_type = GetSeqType(seq_num_h);
	r = spin_list_alignment(r_seq1, r_seq2, "H", "V", data->p_array[i].x, 
			   data->p_array[i].y, "", seq1_type);
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
			    in_comp_spans *input)
{
    seq_result *sip_result;
    d_plot *data;
    text_sim_spans *text_data;
    int i, id;

    if (NULL == (sip_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (d_plot *)xmalloc(sizeof(d_plot))))
	return -1;
    
    if (NULL == (data->p_array = (pt_score *)xmalloc(sizeof(pt_score) * 
						     num_elements)))
	 return -1;

    if (NULL == (text_data = (text_sim_spans *)xmalloc(sizeof(text_sim_spans))))
	return -1;

    id = get_reg_id();
    sip_result->data = data;

     /* data assignment */
    for (i = 0; i < num_elements; i++) {
	data->p_array[i].x = seq1_match[i];
	data->p_array[i].y = seq2_match[i];
	data->p_array[i].score = match_score[i];
    }

    data->win_len = window_length;
    data->n_pts = num_elements;
    data->dim.x0 = start_h;
    data->dim.x1 = end_h;
    data->dim.y0 = start_v;
    data->dim.y1 = end_v;

    sip_result->text_data = text_data;
    text_data->min_score = min_score;

    /* need to initialise this here */
    sip_result->output = NULL;

    /* save global sip_seq into sip_result structure */
    sip_result->seq_id[HORIZONTAL] = GetSeqId(seq1_num);
    sip_result->seq_id[VERTICAL] = GetSeqId(seq2_num);

    sip_result->input = (void *)input; 
    sip_result->id = id; 
    sip_result->graph = SEQ_DOT;

    sip_result->pr_func = dot_plot_middot_func;
    sip_result->op_func = similar_spans_callback;
    sip_result->txt_func = similar_spans_text_func;
    
    seq_register(seq1_num, similar_spans_callback, (void *)sip_result, 
		 SEQ_PLOT_PERM, id);
    seq_register(seq2_num, similar_spans_callback, (void *)sip_result, 
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
				  int *id)
{
    in_comp_spans *input = NULL;
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
   
    vfuncheader("find similar spans");
    
    if (NULL == (seq1_match = (int *)xmalloc(max_matches * sizeof(int))))
	goto error;
    if (NULL == (seq2_match = (int *)xmalloc(max_matches * sizeof(int))))
	goto error;
    if (NULL == (match_score = (int *)xmalloc(max_matches * sizeof(int))))
	goto error;
    if (NULL == (input = (in_comp_spans *)xmalloc(sizeof(in_comp_spans))))
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

    Compare_Spans(seq1, seq2, seq1_len, seq2_len, start_h, end_h, 
		  start_v, end_v, max_matches, same_seq, 
		  win_len, min_match, 1, 0,
		  &seq1_match, &seq2_match, &match_score, &n_matches);

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
    vTcl_DStringAppend(&input_params, "horizontal %s: %s \nvertical %s: %s\n"
	    "window length %d min match %d number of matches %d", 
	    GetSeqLibraryName(seq1_num), 
	    GetSeqName(seq1_num), 
	    GetSeqLibraryName(seq2_num), 
	    GetSeqName(seq2_num), 
	    win_len, min_match, n_matches);
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (-1 == (*id = store_sip_similar_spans(seq1_num, seq2_num, win_len,
					     min_match, start_h, end_h, 
					     start_v, end_v,
					     seq1_match, seq2_match, 
					     match_score, n_matches,
					     input))) {
	goto error;
    }
    
    if (seq1_match)
	xfree (seq1_match);
    if (seq2_match)
	xfree (seq2_match);
    if (match_score)
	xfree(match_score);
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
				char *raster_win, 
				int raster_id,
				char *colour, 
				int line_width)
{
    char *opts[5];
    seq_result *result;
    d_plot *data;
   
    if (NULL == (opts[1] = (char *)xmalloc((strlen(colour)+1) * 
					   sizeof(char))))
	return -1;
    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return -1;

    opts[0] = "-fg";
    strcpy(opts[1], colour);
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", line_width);
    opts[4] = NULL;

    result = result_data(result_id, GetSeqNum(seq_id_h));
    data = result->data;    

    init_dot_plot(interp, seq_id_h, seq_id_v, result_id, "similar spans",
		  raster_win, raster_id, opts, 4, DOT, data->dim);
    
    xfree(opts[1]);
    xfree(opts[3]);

    return 0;
}
