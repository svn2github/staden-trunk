#include <tcl.h>
#include <string.h>
#include <float.h>

#include "seq_raster.h"
#include "misc.h"
#include "xalloc.h"
#include "seq_reg.h"
#include "nip_results.h"
#include "tkRaster.h"
#include "text_output.h"
#include "codon_content.h"
#include "renz_utils.h"
#include "nip_globals.h"
#include "tcl_utils.h"
#include "nip_raster.h"
#include "array_arith.h"
#include "seq_results.h"
#include "nip_gene_search.h"
#include "seq_plot_funcs.h"

void plot_gene_search_callback(int seq_num, void *obj, seq_reg_data *jdata);

void plot_gene_search_shutdown(Tcl_Interp *interp,
			       seq_result *result,
			       char *raster_win,
			       int seq_num)
{
    gene_search *data = result->data;
    out_raster *output = result->output;
    char *tmp;
    seq_reg_key_name info;
    static char buf[80];
    int raster_id;
    Tk_Raster *raster;
    Tcl_CmdInfo info1;
    double wx0, wy0, wx1, wy1;
    RasterResult *raster_result;
    in_plot_gene_search *input;

    input = result->input;

    /* determine raster_id and raster_result structure */
    Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));
    raster_result = raster_id_to_result(raster_id);

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(result->id, (seq_reg_data *)&info, 0);

    /* need to deregister sequence */
    seq_deregister(seq_num, plot_gene_search_callback, (seq_result *)result);

    /* 
     * only bother replotting the raster if there are still results in the
     * raster
     */
    if (raster_result && raster_result->num_results > 1) {
   
	ReplotAllCurrentZoom(interp, raster_win);
	
	tmp = get_default_string(interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
	if (TCL_OK != Tcl_VarEval(interp, "seq_result_list_update ", 
				  tmp, NULL)){
	    verror(ERR_WARN, "gene search", "shutdown: %s \n", Tcl_GetStringResult(interp));
	}
	
	Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
	raster_id = atoi(Tcl_GetStringResult(interp));

	if (TCL_OK != Tcl_VarEval(interp, "RemoveRasterResultKey ", raster_win,
				  " {", info.line, "}", NULL))
	    verror(ERR_WARN, "gene search", "shutdown %s \n", Tcl_GetStringResult(interp));
	

	Tcl_GetCommandInfo(interp, raster_win, &info1);
	raster = (Tk_Raster*)info1.clientData;

	/* find original y before reset size */
	RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

	/* update the scale of the remaining results in the raster window */
	SeqReSetRasterWindowSize(interp, raster_win, result->graph);
	ReSetRasterWindowWorld(interp, raster_win, wy1, result->graph);
	ReplotAllRasterWindow(interp, raster_win);
    }

    free(data->top); /* strdup */
    
    xfree(data->p_array);
    xfree(result->data);

    xfree(output->configure[0]);
    if (output->n_configure == 2) 
	xfree(output->configure[1]);
    xfree(output->configure);
    free(input->params);
    xfree(result->input);
    xfree(result->output);
    xfree(result);

    if (raster_result) 
	DeleteResultFromRaster(raster_result);
    
}

void
plot_gene_search_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_plot_gene_search *input = result->input;
    out_raster *output = result->output;
    gene_search *data = result->data;
    int id = result->id;
    char cmd[1024];

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Plot gene search");
	break;
	
    case SEQ_KEY_NAME:
	if (result->frame) {
	    sprintf(jdata->name.line, "gene f%d #%d", result->frame, result->id);
	} else {
	    /* single plot encompassing all frames */
	    sprintf(jdata->name.line, "gene #%d", result->id);
	}
	break;

    case SEQ_GET_BRIEF:
	if (result->frame) {
	    sprintf(jdata->name.line, "gene: seq=%s frame=%d", 
		    GetSeqName(GetSeqNum(result->seq_id[0])), result->frame);
	} else {
	    sprintf(jdata->name.line, "gene: seq=%s", 
		    GetSeqName(GetSeqNum(result->seq_id[0])));

	}
	break;

    case SEQ_GET_OPS:
	if (output->hidden) {
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
	case 3: /* hide all */
	    output->hidden = 1;
	    ReplotAllCurrentZoom(output->interp, output->raster_win);
	    break;
	case 4: /* reveal all */
	    output->hidden = 0;
	    ReplotAllCurrentZoom(output->interp, output->raster_win);
	    break;
	case 5: /* remove */ 
	    {
		Tcl_Interp *interp = output->interp;
		plot_gene_search_shutdown(interp, result, 
					  output->raster_win, seq_num);
		break;
	    }
	}
	break;
    case SEQ_PLOT:
	result->pr_func(result, (seq_reg_plot *)jdata);
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
	    {
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
		pt.x = get_default_int(interp, nip_defs,
					w("RASTER.PLOT_WIDTH"));
		pt.y = get_default_double(interp, nip_defs,
					   w("RASTER.PLOT_HEIGHT"));

		jdata->info.result = (void *)&pt;
		break;
	    }
	}
	break;
    case SEQ_HIDE: {
	Tcl_Interp *interp = output->interp;
	output->hidden = 1;
	ReplotAllCurrentZoom(interp, output->raster_win);
	break;
    }
    case SEQ_REVEAL: 
	output->hidden = 0;
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: {
	Tcl_Interp *interp = output->interp;
	
	plot_gene_search_shutdown(interp, result, output->raster_win,
				  seq_num);
	break;
	}
    }
}

void plot_gene_search_text_func(void *obj)
{
    seq_result *result = (seq_result *) obj;
    gene_search *data = result->data;
    int n_pts = data->n_pts;
    int i;

    for (i = 0; i < n_pts; i++) {
	UpdateTextOutput();
	vmessage("Position %10d score %.5g \n", data->p_array[i].pos,
		 data->p_array[i].score);
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
		      int type)
{
    int i;
    seq_result *result;
    int id;
    gene_search *data;

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (gene_search *)xmalloc(sizeof(gene_search))))
	return -1;

    if (NULL == (data->p_array = (p_score *)xmalloc(sizeof(p_score) * 
						    num_elements * 3)))
	return -1;

    result->data = data;

    id = get_reg_id();

    /* data assignment */
    for (i = 0; i < num_elements; i++) {
	data->p_array[i].pos = (start-1) + (i * 3) + frame;
	data->p_array[i].score = score[i];
    }

    data->n_pts = num_elements;
    data->dim.x0 = start;
    data->dim.x1 = end;
    data->dim.y0 = min; /* min of all three frames */
    data->dim.y1 = max; /* max of all three frames */
    
    result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    result->seq_id[VERTICAL] = -1;

    result->id = id; 
    result->graph = SEQ_GENE;
    result->input = (void *)input; 

    /* need to initialise this here */
    result->output = NULL;

    if (top) {
	data->top = strdup(top);
    } else {
	data->top = NULL;
    }

    result->pr_func = gene_search_plot_func;
    result->op_func = plot_gene_search_callback;
    result->txt_func =plot_gene_search_text_func;
    result->type = type;
    result->frame = frame;
    
    seq_register(seq_num, plot_gene_search_callback, (void *)result, 
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
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "window length %d codon table %s\n%s\n",
		       GetSeqName(seq_num), start, end,
		       win_len, codon_table, tmp);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input1->params = strdup(Tcl_DStringValue(&input_params)); 
    input2->params = strdup(Tcl_DStringValue(&input_params)); 
    input3->params = strdup(Tcl_DStringValue(&input_params)); 
    input1->method = "codon preference";
    input2->method = "codon preference";
    input3->method = "codon preference";
    Tcl_DStringFree(&input_params);

    if (-1 == DoCodonPref(seq, seq_len, codon_table, win_len, 
			  start, end, option, &results)) {
	verror(ERR_WARN, "CodonPref", "Failed DoCodonPref\n");
	goto fail;
    }

    
    id[0] = store_gene_search(seq_num, start, end, 1, input1, 
			      results->frame1, 
			      results->top, results->num_results, 
			      results->min, results->max,
			      SEQ_TYPE_CODONPREF);
    id[1] = store_gene_search(seq_num, start, end, 2, input2, 
			      results->frame2, 
			      results->top, results->num_results, 
			      results->min, results->max,
			      SEQ_TYPE_CODONPREF);
    id[2] = store_gene_search(seq_num, start, end, 3, input3, 
			      results->frame3, 
			      results->top, results->num_results, 
			      results->min, results->max,
			      SEQ_TYPE_CODONPREF);

     free_CodRes(results);
     return 0;

fail:
    xfree(input1);
    xfree(input2);
    xfree(input3);
    return -1;
}

int init_nip_gene_search_plot(Tcl_Interp *interp, 
			      char *rasters, 
			      char *raster_ids,
			      int seq_id, 
			      char *r_id,
			      char *colours,
			      int line_width)
{

    int num_id;
    int retval        = -1;
    char **result_id  = NULL;
    char **raster_win = NULL;
    char **raster_id  = NULL;
    char **colour     = NULL;
    

    if (Tcl_SplitList(interp, rasters, &num_id, &raster_win) != TCL_OK)
	goto cleanup;
    if (Tcl_SplitList(interp, raster_ids, &num_id, &raster_id) != TCL_OK)
	goto cleanup;
    if (Tcl_SplitList(interp, colours, &num_id, &colour) != TCL_OK)
	goto cleanup;
    if (Tcl_SplitList(interp, r_id, &num_id, &result_id) != TCL_OK)
	goto cleanup;


    init_gene_search_raster(interp, num_id, raster_win, raster_id, seq_id, 
			    result_id, colour, line_width);
    retval = 0;


cleanup:
    if(result_id)  Tcl_Free((char *)result_id);
    if(raster_win) Tcl_Free((char *)raster_win);
    if(raster_id)  Tcl_Free((char *)raster_id);
    if(colour)     Tcl_Free((char *)colour);
    return retval;
}

int init_nip_author_test_create(Tcl_Interp *interp,
				int seq_id,
				int start,
				int end,
				char *codon_table,
				double error,
				int *id)
{
    char *seq;
    int seq_len, seq_num;
    in_plot_gene_search *input1;
    in_plot_gene_search *input2;
    in_plot_gene_search *input3;
    CodRes *results;
    Tcl_DString input_params;

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

    if (-1 == DoAuthorTest(seq, seq_len, codon_table, error, 
			   start, end, &results)) {
	verror(ERR_WARN, "AuthorTest", "Failed DoAuthorTest\n");
	goto fail;
    }

    id[0] = store_gene_search(seq_num, start, end, 1, input1, results->frame1, 
			      results->top, results->num_results, 
			      results->min, results->max,
			      SEQ_TYPE_AUTHOR);
    id[1] = store_gene_search(seq_num, start, end, 2, input2, results->frame2, 
			      results->top, results->num_results, 
			      results->min, results->max,
			       SEQ_TYPE_AUTHOR);
    id[2] = store_gene_search(seq_num, start, end, 3, input3, results->frame3, 
			      results->top, results->num_results, 
			      results->min, results->max,
			      SEQ_TYPE_AUTHOR);

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

    if (-1 == DoPosBaseBias(seq, seq_len, win_len, start, end, &results)) {
	verror(ERR_WARN, "BaseBias", "Failed DoPosBaseBias\n");
	goto fail;
    }

    *id = store_gene_search(seq_num, start, end, 0, input, results->frame1, 
			    NULL, results->num_results, 
			    results->min, results->max,
			    SEQ_TYPE_BASEBIAS);

    free_CodRes1(results);
    return 0;
fail:
    xfree(input);
    return -1;
}

int init_nip_base_bias_plot(Tcl_Interp *interp, 
			    char *raster1_win, 
			    char *raster_id1,
			    int seq_id, 
			    char *r_id,
			    char *colour1,
			    int line_width)
{
    init_gene_search_raster(interp, 1, &raster1_win, &raster_id1, seq_id, 
			    &r_id, &colour1, line_width);
    return 0;
}



