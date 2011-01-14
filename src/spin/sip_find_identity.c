#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "seq_raster.h"
#include "seq_reg.h"
#include "sip_hash.h"
#include "probs.h"
#include "seq_results.h"
#include "sip_results.h"
#include "sip_find_identity.h"
#include "misc.h"
#include "sequence_formats.h"
#include "tkRaster.h"
#include "text_output.h"
#include "readpam.h"
#include "tcl_utils.h"
#include "sip_globals.h"
#include "dna_utils.h"
#include "seq_plot_funcs.h"
#include "sequence_pair_display.h"

void identities_callback(int seq_num, void *obj, seq_reg_data *jdata);

void identities_shutdown(Tcl_Interp *interp,
			 int seq_num,
			 seq_result *result,
			 char *raster_win,
			 RasterResult *raster_result)
{
    seq_reg_key_name info;
    static char buf[80];
    double wx0, wy0, wx1, wy1;
    Tk_Raster *raster;
    Tcl_CmdInfo info1;

    Tcl_GetCommandInfo(interp, raster_win, &info1);
    raster = (Tk_Raster*)info1.clientData;
    
    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(result->id, (seq_reg_data *)&info, 0);

    /* deregister horizontal sequence */
    seq_deregister(GetSeqNum(result->seq_id[HORIZONTAL]), identities_callback,
		   (seq_result *)result);

    /* deregister vertical sequence */
    seq_deregister(GetSeqNum(result->seq_id[VERTICAL]), identities_callback, 
		   (seq_result *)result);

    /* 
     * only bother replotting the raster if there are still results in the
     * raster
     */
    if (raster_result && raster_result->num_results > 1) {
	if (TCL_OK != Tcl_VarEval(interp, "RemoveRasterResultKey ", raster_win,
				  " {", info.line, "}", NULL))
	    verror(ERR_WARN, "quick_scan_shutdown", "%s \n", 
		   Tcl_GetStringResult(interp));
	
	/* find original y before reset size */
	RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

	/* update the scale of the remaining results in the raster window */
	SeqReSetRasterWindowSize(interp, raster_win, result->graph);
	ReSetRasterWindowWorld(interp, raster_win, wy1, result->graph);

	ReplotAllRasterWindow(interp, raster_win);
    }
}

void identities_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    d_plot *data = result->data;
    in_find_identities *input = result->input;
    out_raster *output = result->output;
    /* int num_pts = data->n_pts;*/
    text_find_identities *text = result->text_data;
    int id = result->id;
    char cmd[1024];
    char *tmp;

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Find matching words");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "matching words #%d", result->id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "matching words: hori=%s vert=%s", 
		GetSeqBaseName(GetSeqNum(result->seq_id[HORIZONTAL])),
		GetSeqBaseName(GetSeqNum(result->seq_id[VERTICAL])));
	break;
    case SEQ_GET_OPS:
	if (output->hidden) {
	    jdata->get_ops.ops = "Information\0List results\0"
		"Tabulate scores\0PLACEHOLDER\0PLACEHOLDER\0PLACEHOLDER\0"
		"Reveal\0SEPARATOR\0Remove\0";
	} else if ((seq_get_type(id) == SEQ_PLOT_TEMP) && !get_replot_temp()) {
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
	    vfuncheader("input parameters");
	    vmessage("%sNumber of matches %d\n", input->params, data->n_pts);
	    break;
	case 1: /* results */
	    Tcl_Eval(output->interp, "SetBusy");
	    vfuncheader("results");
	    result->txt_func(result);
	    Tcl_Eval(output->interp, "ClearBusy");
	    break;
	case 2: /* scores */
	    {
		Tcl_Eval(output->interp, "SetBusy");
		vfuncheader("scores");
		CalcIdentityProbs(result, text->word_len);
		Tcl_Eval(output->interp, "ClearBusy");
		break;
	    }
	case 3: /* configure */
	    sprintf(cmd, "RasterConfig %d", id);
	    if (TCL_OK != Tcl_Eval(output->interp, cmd)){
		puts(Tcl_GetStringResult(output->interp));
	    }
	    break;
	case 4: /* display sequences */
	    SequencePairDisplay(output->interp, output->raster_win, id, 
				result->seq_id[HORIZONTAL], 
				result->seq_id[VERTICAL]);
	    break;
	case 5: /* hide all */
	    output->hidden = 1;
	    ReplotAllCurrentZoom(output->interp, output->raster_win);
	    break;
	case 6: /* reveal all */
	    output->hidden = 0;
	    ReplotAllCurrentZoom(output->interp, output->raster_win);
	    break;
	case 7: /* remove */
	    {
		int raster_id;
		RasterResult *raster_result;

		Tcl_VarEval(output->interp, "GetRasterId ", output->raster_win, NULL);
		raster_id = atoi(Tcl_GetStringResult(output->interp));
		raster_result = raster_id_to_result(raster_id);

		identities_shutdown(output->interp, seq_num, result, 
				    output->raster_win, raster_result);

		if (raster_result && raster_result->num_results > 1) {
		    
		    ReplotAllCurrentZoom(output->interp, output->raster_win);
		    tmp = get_default_string(output->interp, tk_utils_defs, 
					     w("RASTER.RESULTS.WIN"));
		    Tcl_VarEval(output->interp, "seq_result_list_update ", tmp,
				NULL);
		}
		DestroySequencePairDisplay(output->interp, id);
		free(input->params);
		xfree(result->text_data);
		SipFreeResult(result);
		if (raster_result) 
		    DeleteResultFromRaster(raster_result);
		break;
	    }
	}
	break;
    case SEQ_PLOT: 
	{
	    result->pr_func(result, NULL);
	    break;
	}
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case OUTPUT:
	    jdata->info.result = (void *)output;
	    break;
	case INPUT: 
	    jdata->info.result = (void *)input;
	    break;
	case DIMENSIONS: 
	    {
#ifdef REMOVE
		d_line *dim;
		if (NULL == (dim = (d_line *)xmalloc(sizeof(d_line))))
		    return;
		dim = &data->dim;
		jdata->info.result = (void *)dim;
#endif
		jdata->info.result = (void *)&data->dim;
		break;
	    }
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
	{
	    RasterResult *raster_result;
	    int raster_id;

	    Tcl_VarEval(output->interp, "GetRasterId ", output->raster_win, NULL);
	    raster_id = atoi(Tcl_GetStringResult(output->interp));
	    raster_result = raster_id_to_result(raster_id);

	    /* called when hide all in dot plot */
	    if ((seq_get_type(id) == SEQ_PLOT_TEMP) && (!get_replot_temp())) {
		identities_shutdown(output->interp, seq_num, result, 
				    output->raster_win, raster_result);
	    } else {
		output->hidden = 1;
	    }
	    break;
	}
    case SEQ_REVEAL: 
	output->hidden = 0;
	break;
    case SEQ_QUIT:
    case SEQ_DELETE:
	{
	    RasterResult *raster_result;
	    int raster_id;

	    Tcl_VarEval(output->interp, "GetRasterId ", output->raster_win, NULL);
	    raster_id = atoi(Tcl_GetStringResult(output->interp));
	    raster_result = raster_id_to_result(raster_id);

	    identities_shutdown(output->interp, seq_num, result, 
				output->raster_win, raster_result);
	    if (raster_result && raster_result->num_results > 1) {
		ReplotAllCurrentZoom(output->interp, output->raster_win);
		tmp = get_default_string(output->interp, tk_utils_defs, 
					 w("RASTER.RESULTS.WIN"));
		Tcl_VarEval(output->interp, "seq_result_list_update ", tmp, NULL);
	    }
	    DestroySequencePairDisplay(output->interp, id);
	    free(input->params);
	    xfree(result->text_data);
	    SipFreeResult(result);
	    if (raster_result) 
		DeleteResultFromRaster(raster_result);

	    break;
	}
    }
}


void identities_text_func(void *obj)
{
    seq_result *result = (seq_result *) obj;
    d_plot *data = result->data;
    int n_pts = data->n_pts;
    int i;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    char *tmp_str;
    int max_str;
    int seq_num_h, seq_num_v;

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
	
    if (NULL == (tmp_str = (char *)xmalloc(max_str * sizeof(char))))
	return;

    for (i = 0; i < n_pts; i++) {
	UpdateTextOutput();
	vmessage("Positions %10d h %10d v and length %10d\n", 
		 data->p_array[i].x, (int)data->p_array[i].y, 
		 data->p_array[i].score);
	strncpy(tmp_str, &seq1[data->p_array[i].x - 1], 
			      data->p_array[i].score);
	tmp_str[data->p_array[i].score] = '\0';
	vmessage("%s\n", tmp_str);
    }
    xfree(tmp_str);
}

void identities_recalc_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    d_plot *data = result->data;
    out_raster *output = result->output;
    text_find_identities *text = result->text_data;
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    /*int index;*/
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int seq1_type, seq2_type;
    Tk_Raster *printData;
    int max_matches = get_max_matches();
    int max_hash = 7, word_length;
    int *seq1_match, *seq2_match; 
    int *len_match;
    int n_matches;
    int same_seq;
    int index1, index2;

    if (output->hidden) {
	return;
    }
    Tcl_GetCommandInfo(output->interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;

    SetDrawEnviron(output->interp, raster, output->env_index);

    printData = raster;

    index1 = GetSeqNum(result->seq_id[HORIZONTAL]);
    index2 = GetSeqNum(result->seq_id[VERTICAL]);

    /* if either index is -1; the corresponding seq must have been deleted */
    if ((index1 == -1) || (index2 == -1)) {
	return;
    }
    seq1 = GetSeqSequence(index1);
    seq2 = GetSeqSequence(index2);
    seq1_len = GetSeqLength(index1);
    seq2_len = GetSeqLength(index2);
    seq1_type = GetSeqType(index1);
    seq2_type = GetSeqType(index2);

    if (seq1_len == seq2_len) {
	if (strcmp(seq1, seq2) == 0) {
	    same_seq = 1;
	} else {
	    same_seq = 0;
	}
    } else {
	same_seq = 0;
    }    

    if (seq1_type != seq2_type) {
	verror(ERR_WARN, "find matching words", 
	       "sequences must both be either DNA or protein");
	return;
    } else if (seq1_type == PROTEIN) {
	set_char_set(PROTEIN);
	/* set identity matrix */

	if (-1 == (set_matrix_identity(PROTEIN))) {
	    verror(ERR_WARN, "set score matrix", 
		   "unable to set identity score matrix");
	    return;
	}
	set_score_matrix(get_matrix_identity(PROTEIN));

	max_hash = 3;
    } else if (seq1_type == DNA) {
	set_char_set(DNA);
	if (-1 == (set_matrix_identity(DNA))) {
	    verror(ERR_WARN, "set score matrix", 
		   "unable to set identity score matrix");
	    return;
	}
	set_score_matrix(get_matrix_identity(DNA));

	max_hash = 7;
    }
/*
    max_matches = (int)FindExpectedProb(seq1, seq2, seq1_len, seq2_len, 
				   input->word_len, seq1_type);

    if (NULL == (seq1_match = (int *)xmalloc(max_matches * sizeof(int)))) {
	return;
    }
    if (NULL == (seq2_match = (int *)xmalloc(max_matches * sizeof(int)))) {
	return;
    }
    if (NULL == (len_match = (int *)xmalloc(max_matches * sizeof(int)))) {
	return;
    }
*/
    word_length = MIN(text->word_len, max_hash);

    sip_hash(seq1, seq2, GetSubSeqStart(index1), GetSubSeqEnd(index1),
	     GetSubSeqStart(index2), GetSubSeqEnd(index2),
	     max_matches, text->word_len, 
	     word_length, seq1_type, same_seq, &seq1_match, 
	     &seq2_match, &len_match, &n_matches, RasterDrawPoint, 
	     (void *)printData);
    if (n_matches < 0) {
	verror(ERR_WARN, "find matching words", 
	       "failed in find matching words \n");
	return;
    }

    /* save n_matches here incase of too_many_matches */
    data->n_pts = n_matches;
    tk_RasterRefresh(raster);
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
			 int too_many_matches)
{
    seq_result *result;
    d_plot *data;
    int i, id;
    text_find_identities *text_data;
    

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (d_plot *)xmalloc(sizeof(d_plot))))
	return -1;
    
    if (NULL == (text_data = (text_find_identities *)xmalloc(sizeof(text_find_identities))))
	return -1;

    if (!too_many_matches) {
	if (NULL == (data->p_array = (pt_score *)xmalloc(sizeof(pt_score) * 
							 num_elements)))
	    return -1;

	/* data assignment */
	for (i = 0; i < num_elements; i++) {
	    data->p_array[i].x = seq1_match[i];
	    data->p_array[i].y = seq2_match[i];
	    data->p_array[i].score = match_score[i];
	}
    } else {
	data->p_array = NULL;
    }

    id = get_reg_id();
    result->data = data;
    
    data->n_pts = num_elements;
    data->dim.x0 = start_h;
    data->dim.x1 = end_h;
    data->dim.y0 = start_v;
    data->dim.y1 = end_v;

    result->text_data = text_data;
    text_data->word_len = word_length;

    /* need to initialise this here */
    result->output = NULL;

    result->seq_id[HORIZONTAL] = GetSeqId( seq1_num);
    result->seq_id[VERTICAL] = GetSeqId(seq2_num);

    result->input = (void *)input; 
    result->id = id; 
    result->graph = SEQ_DOT;

    result->op_func = identities_callback;
    result->txt_func = identities_text_func;

    /* too_many_matches is always false at the present time */
    if (too_many_matches) {
	result->pr_func = identities_recalc_func;
	seq_register(seq1_num, identities_callback, (void *)result, 
		     SEQ_PLOT_TEMP, id);
	seq_register(seq2_num, identities_callback, (void *)result, 
		     SEQ_PLOT_TEMP, id);
    } else {
	result->pr_func = dot_plot_scoreline_func;
	seq_register(seq1_num, identities_callback, (void *)result, 
		     SEQ_PLOT_PERM, id);
	seq_register(seq2_num, identities_callback, (void *)result, 
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
				   int *id)
{
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int max_matches;
    int n_matches;
    int *seq1_match = NULL;
    int *seq2_match = NULL; 
    int *len_match = NULL;
    int max_hash = 7, word_length;
    in_find_identities *input = NULL;
    int same_seq;
    int too_many_matches = 0;
    int seq1_num, seq2_num;
    int seq1_type, seq2_type;
    int sub1_len, sub2_len;
    double prob;
    Tcl_DString input_params;    
   
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
    vTcl_DStringAppend(&input_params, "horizontal %s: %s\nvertical %s: %s\n"
	    "word length %d \n", 
	    GetSeqLibraryName(seq1_num), 
	    GetSeqName(seq1_num), 
	    GetSeqLibraryName(seq2_num),
	    GetSeqName(seq2_num), word_len);
    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);
    
    if (-1 == (*id = store_matching_words(seq1_num, seq2_num,start_h,
					  end_h, start_v, end_v,word_length,
					  input, seq1_match, seq2_match, 
					  len_match, n_matches,
					  too_many_matches))) {
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
				 char *raster_win, 
				 int raster_id,
				 char *colour, 
				 int line_width)
{
    char *opts[7];
    seq_result *result;
    d_plot *data;

    if (NULL == (opts[1] = (char *)xmalloc((strlen(colour)+1) * 
					   sizeof(char))))
	return -1;
    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return -1;
    if (NULL == (opts[5] = (char *)xmalloc(15 * sizeof(char))))
	return -1;

    opts[0] = "-fg";
    strcpy(opts[1], colour);
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", line_width);
    opts[4] = "-capstyle";
    strcpy(opts[5], "round");
    opts[6] = NULL;

    result = result_data(result_id, GetSeqNum(seq_id_h));
    data = result->data;    
 
    init_dot_plot(interp, seq_id_h, seq_id_v, result_id, "matching words", 
		  raster_win, raster_id, opts, 6, LINE, data->dim);
    
    xfree(opts[1]);
    xfree(opts[3]);
    xfree(opts[5]);

    return 0;
}
