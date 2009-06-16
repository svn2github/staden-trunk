#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "seq_raster.h"
#include "seq_reg.h"
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

void quick_scan_callback(int seq_num, void *obj, seq_reg_data *jdata);

void quick_scan_shutdown(Tcl_Interp *interp,
			 int seq_num,
			 seq_result *result,
			 char *raster_win,
			 RasterResult *raster_result)
{
    int raster_id;
    seq_reg_key_name info;
    static char buf[80];
    double wx0, wy0, wx1, wy1;
    Tk_Raster *raster;
    Tcl_CmdInfo info1;

    Tcl_GetCommandInfo(interp, raster_win, &info1);
    raster = (Tk_Raster*)info1.clientData;

#ifdef DEBUG
    printf("START quick_scan_delete_result\n");
#endif

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(result->id, (seq_reg_data *)&info, 0);

     /* deregister horizontal sequence */
    seq_deregister(GetSeqNum(result->seq_id[HORIZONTAL]), quick_scan_callback, (seq_result *)result);

    /* deregister vertical sequence */
    seq_deregister(GetSeqNum(result->seq_id[VERTICAL]), quick_scan_callback, (seq_result *)result);

    Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));

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

void quick_scan_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_quick_scan *input = result->input;
    out_raster *output = result->output;
    d_plot *data = result->data;
    int id = result->id;
    char cmd[1024];
    char *tmp;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Find best diagonals");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "diagonals #%d", result->id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "diagonals: hori=%s vert=%s", 
		GetSeqBaseName(GetSeqNum(result->seq_id[HORIZONTAL])),
		GetSeqBaseName(GetSeqNum(result->seq_id[VERTICAL])));
	break;

    case SEQ_GET_OPS:
	if (output->hidden) {
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
	    Tcl_Eval(output->interp, "SetBusy");
	    vfuncheader("results");
	    result->txt_func(result);
	    Tcl_Eval(output->interp, "ClearBusy");
	    break;
	case 2: /* configure */
	    sprintf(cmd, "RasterConfig %d", id);
	    if (TCL_OK != Tcl_Eval(output->interp, cmd)){
		puts(Tcl_GetStringResult(output->interp));
	    }
	    break;
	case 3: /* display sequences */
	    SequencePairDisplay(output->interp, output->raster_win, id, 
				result->seq_id[HORIZONTAL], 
				result->seq_id[VERTICAL]);
	    break;
	case 4: /* hide all */
	    output->hidden = 1;
	    ReplotAllCurrentZoom(output->interp, output->raster_win);
	    break;
	case 5: /* reveal all */
	    output->hidden = 0;
	    ReplotAllCurrentZoom(output->interp, output->raster_win);
	    break;
	case 6: /* remove */
	    {
		Tcl_Interp *interp = output->interp;
		int raster_id;
		RasterResult *raster_result;

		Tcl_VarEval(interp, "GetRasterId ", output->raster_win, NULL);
		raster_id = atoi(Tcl_GetStringResult(interp));
		raster_result = raster_id_to_result(raster_id);

		quick_scan_shutdown(output->interp, seq_num, result, 
				    output->raster_win, raster_result);
		if (raster_result && raster_result->num_results > 1) {
		    
		    ReplotAllCurrentZoom(output->interp, output->raster_win);
		    tmp = get_default_string(output->interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
		    Tcl_VarEval(output->interp, "seq_result_list_update ", tmp, NULL);
		}
		DestroySequencePairDisplay(output->interp, id);
		free(input->params);
		SipFreeResult(result);
		if (raster_result) 
		    DeleteResultFromRaster(raster_result);

	    }
	    break;
	}
	break;
    case SEQ_PLOT:
	{
	    RasterResult *raster_result;
	    int raster_id; 
	    int save_results = 1;
	    /* determine raster_id and raster_result structure */
	    Tcl_VarEval(output->interp, "GetRasterId ", output->raster_win, 
			NULL);
	    raster_id = atoi(Tcl_GetStringResult(output->interp));
	    raster_result = raster_id_to_result(raster_id);

	    /* 
	     * save_results is always true. If need to go back to the 
	     * previous code, must look at the backups (1999.1p1)
	     */
	    if (save_results) {
		result->pr_func(result, NULL);
	    } else {
		if (Tcl_GetCommandInfo(output->interp, output->raster_win, 
				       &info) == 0) 
		    return;
		raster = (Tk_Raster*)info.clientData;
		if (!get_replot_temp() && output->replot) {
		    /* plot result only once */
		    
		    /* reset plot func to null */
		    RasterDeletePlotFunc(raster);
		    result->pr_func(result, NULL);
		    output->replot = 0;
		} else if (!get_replot_temp()) {
		    
		    /* reset plot func to null */
		    RasterDeletePlotFunc(raster);
		    quick_scan_shutdown(output->interp, seq_num, result, 
					output->raster_win, raster_result);
		} else {
		    /* necessary to allow scrolling */
		    RasterInitPlotFunc(raster, SeqRasterPlotFunc);
		    result->pr_func(result, NULL);
		    output->replot = 0;
		}
	    }
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
	{
	    RasterResult *raster_result;
	    int raster_id; 

	    /* determine raster_id and raster_result structure */
	    Tcl_VarEval(output->interp, "GetRasterId ", output->raster_win, 
			NULL);
	    raster_id = atoi(Tcl_GetStringResult(output->interp));
	    raster_result = raster_id_to_result(raster_id);

	    if (get_replot_temp()) {
		output->hidden = 1;
	    } else {
		quick_scan_shutdown(output->interp, seq_num, result, 
				    output->raster_win, raster_result);
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

	    /* determine raster_id and raster_result structure */
	    Tcl_VarEval(output->interp, "GetRasterId ", output->raster_win, 
			NULL);
	    raster_id = atoi(Tcl_GetStringResult(output->interp));
	    raster_result = raster_id_to_result(raster_id);

	    quick_scan_shutdown(output->interp, seq_num, result, 
				output->raster_win, raster_result);
	    if (raster_result && raster_result->num_results > 1) {

		ReplotAllCurrentZoom(output->interp, output->raster_win);
		tmp = get_default_string(output->interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
		Tcl_VarEval(output->interp, "seq_result_list_update ", tmp, NULL);
	    }
	    DestroySequencePairDisplay(output->interp, id);
	    free(input->params);
	    SipFreeResult(result);
	    if (raster_result) 
		DeleteResultFromRaster(raster_result);
	    
	}
    }
}

void quick_scan_text_func(void *obj)
{
    seq_result *result = (seq_result *) obj;
    d_plot *data = result->data;
    int n_pts = data->n_pts;
    int i;
    char *seq1;
    char *seq2;
    int seq1_len, seq2_len;
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

    for (i = 0; i < n_pts; i++) {
	UpdateTextOutput();
	vmessage("Positions %10d h %10d v \n", 
	       data->p_array[i].x, (int)data->p_array[i].y);
    }
}

void quick_scan_recalc_func(void *obj, seq_reg_plot *plot)
{
#ifdef RECALC
    seq_result *result = (seq_result *) obj;
    in_quick_scan *input = result->input;
    out_raster *output = result->output;
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int seq1_type, seq2_type;
    Tk_Raster *printData;
    int max_matches = get_max_matches();
    int *seq1_match, *seq2_match; 
    int same_seq;
    int index1, index2;
    int n_matches;

#ifdef START_DEBUG
    printf("START quick_scan_recalc_func hidden %d\n", output->hidden);
#endif
    if (output->hidden) {
	return;
    }
    if (NULL == (seq1_match = (int *)xmalloc(max_matches * sizeof(int))))
	return;
    if (NULL == (seq2_match = (int *)xmalloc(max_matches * sizeof(int))))
	return;

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
	verror(ERR_WARN, "find best diagonals", 
	       "sequences must both be either DNA or protein");
	return;
    } else if (seq1_type == PROTEIN) {
	set_char_set(PROTEIN);
        set_score_matrix(get_matrix_file(PROTEIN));
    } else if (seq1_type == DNA) {
	set_char_set(DNA);
        set_score_matrix(get_matrix_file(DNA));
    }
   
    n_matches = quick_scan(seq1, seq2, GetSubSeqStart(index1),
			   GetSubSeqEnd(index1), GetSubSeqStart(index2),
			   GetSubSeqEnd(index2), seq1_type, 
			   max_matches, same_seq, input->win_length, 
			   input->min_score, input->word_length, input->sd, 
			   input->save_results, &seq1_match, &seq2_match, 
			   RasterDrawPoint, (void *)printData);
#endif
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
		     int save_results)
{
    seq_result *result;
    d_plot *data;
    int i, id;

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (d_plot *)xmalloc(sizeof(d_plot))))
	return -1;

    if (save_results) {
	if (NULL == (data->p_array = (pt_score *)xmalloc(sizeof(pt_score) * 
							 num_elements)))
	    return -1;
	/* data assignment */
	
	for (i = 0; i < num_elements; i++) {
	    data->p_array[i].x = seq1_match[i];
	    data->p_array[i].y = seq2_match[i];
	    data->p_array[i].score = 0;
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

    /* need to initialise this here */
    result->output = NULL;

    result->seq_id[HORIZONTAL] = GetSeqId(seq1_num);
    result->seq_id[VERTICAL] = GetSeqId(seq2_num);

    result->input = (void *)input; 
    result->id = id; 
    result->graph = SEQ_DOT;

    result->op_func = quick_scan_callback;
    result->txt_func = quick_scan_text_func;
    
    if (save_results) {
	result->pr_func = dot_plot_dot_func;
	seq_register(seq1_num, quick_scan_callback, (void *)result, 
		     SEQ_PLOT_PERM, id);
	seq_register(seq2_num, quick_scan_callback, (void *)result, 
		     SEQ_PLOT_PERM, id);
    } else {
	result->pr_func =  quick_scan_recalc_func;
	seq_register(seq1_num, quick_scan_callback, (void *)result, 
		     SEQ_PLOT_TEMP, id);
	seq_register(seq2_num, quick_scan_callback, (void *)result, 
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
    vTcl_DStringAppend(&input_params, "horizontal %s: %s\nvertical %s: %s\n"
		       "window length %d minimum score %d word length %d minimum sd %f", 
		       GetSeqLibraryName(seq1_num), 
		       GetSeqName(seq1_num), 
		       GetSeqLibraryName(seq2_num), 
		       GetSeqName(seq2_num), win_len, min_match, 
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
	n_matches = quick_scan(seq1, seq2, start_h, end_h, start_v, end_v, 
			       seq1_type, max_matches, same_seq, win_len, 
			       min_match, word_len, min_sd, save_results,
			       &seq1_match, &seq2_match, NULL, NULL);
	if (n_matches == -1) {
	  goto error;
	}
    }

    if (n_matches > 0) {
	if (-1 == (*id = store_quick_scan(seq1_num, seq2_num, start_h, 
					  end_h, start_v, end_v, input, 
					  seq1_match, seq2_match, 
					  n_matches, save_results))) {
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

    init_dot_plot(interp, seq_id_h, seq_id_v, result_id, "diagonals", 
		  raster_win, raster_id, opts, 6, DOT, data->dim);
    
    xfree(opts[1]);
    xfree(opts[3]);
    xfree(opts[5]);

    return 0;
}
