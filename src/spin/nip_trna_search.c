#include <tcl.h>

#include "seq_raster.h"
#include "xalloc.h"
#include "seq_reg.h"
#include "seq_results.h"
#include "nip_results.h"
#include "tkRaster.h"
#include "text_output.h"
#include "nip_globals.h"
#include "trna_search.h"
#include "nip_trna_search.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "nip_raster.h"
#include "string.h"
#include "misc.h"
#include "seq_plot_funcs.h"
#include "dna_utils.h"

void trna_search_callback(int seq_num, void *obj, seq_reg_data *jdata);

void trna_search_shutdown(Tcl_Interp *interp,
			  seq_result *result, 
			  char *raster_win,
			  int seq_num)
{
    stick *data = result->data;
    out_raster *output = result->output;
    in_trna_search *input = result->input;
    TrnaRes **results = result->text_data;
    int i;
    char *tmp;
    seq_reg_key_name info;
    static char buf[80];
    int raster_id;
    RasterResult *raster_result;
     
    /* determine raster_id and raster_result structure */
    Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));
    raster_result = raster_id_to_result(raster_id);

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(result->id, (seq_reg_data *)&info, 0);

    seq_deregister(seq_num, trna_search_callback, (seq_result *)result);

    /* 
     * only bother replotting the raster if there are still results in the
     * raster
     */
    if (raster_result && raster_result->num_results > 1) {
	
	ReplotAllCurrentZoom(interp, raster_win);
	/* DestroyDisplaySequences(interp, id); */
	tmp = get_default_string(interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
	if (TCL_OK != Tcl_VarEval(interp, "seq_result_list_update ", 
				  tmp, NULL)){
	    puts(Tcl_GetStringResult(interp));
	}

	if (TCL_OK != Tcl_VarEval(interp, "RemoveRasterResultKey ", raster_win,
				  " {", info.line, "}", NULL))
	    verror(ERR_WARN, "trna search", "shutdown %s \n", Tcl_GetStringResult(interp));
    }
    for (i = 0; i < MAX_TRNA; i++) {
	xfree(results[i]);
    }
    xfree(results);
    xfree(data->ap_array[0].p_array);
    xfree(data->ap_array);
    xfree(result->data);

    free(input->params);
    xfree(input->t);
    xfree(output->configure[0]);
    xfree(output->configure);
    xfree(result->input);
    xfree(result->output);
    xfree(result);

    if (raster_result) 
	DeleteResultFromRaster(raster_result);
    
}

void trna_search_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_trna_search *input = result->input;
    out_raster *output = result->output;
    stick *data = result->data;
    int id = result->id;
    char cmd[1024];

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "tRNA search");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "trna #%d", result->id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "trna: seq=%s", 
		GetSeqName(GetSeqNum(result->seq_id[0])));
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
		trna_search_shutdown(interp, result, output->raster_win, 
				     seq_num);
		break;
	    }
	}
	break;
    case SEQ_PLOT:
	result->pr_func(result, (seq_reg_plot *)jdata);
	break;
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) 
	    {
	    case OUTPUT:
		jdata->info.result = (void *)output;
		break;
	    case INPUT: 
		jdata->info.result = (void *)input;
		break;
	    case DIMENSIONS: 
		jdata->info.result = (void *)&data->ap_array[0].dim;
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
		    pt.x = get_default_int(interp, nip_defs,
					   w("RASTER.PLOT_WIDTH"));
		    pt.y = get_default_double(interp, nip_defs,
					      w("NIP.TRNA.PLOT_HEIGHT"));
		    
		    jdata->info.result = (void *)&pt;
		    break;
		}
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
	    trna_search_shutdown(interp, result, output->raster_win, seq_num);
	    break;
	}
    }
}

void
trna_search_text_func(void *obj)
{
    seq_result *result = (seq_result *) obj;
    stick *data = result->data;
    TrnaRes **results = result->text_data;
    int n_pts = data->ap_array[0].n_pts;
    int i;
    in_trna_search *input = result->input;
    TrnaSpec *t = input->t;
    int x, score;
    
    for (i = 0; i < n_pts; i++) {
	if ( results[i]->total_cb_score >= t->min_total_cb_score ) 
	    draw_trna(results[i]);
    }
    
    
    for (i = 0; i < n_pts; i++) {
	UpdateTextOutput();
	x = results[i]->aa_left+1;
	score = results[i]->total_bp_score;

	vmessage("Position %10d score %10d\n", x, score);
    }
}

int store_trna_search(int seq_num,
		       in_trna_search *input,
		       int start,
		       int end,
		       TrnaRes **results, 
		       int num_elements,
		       TrnaSpec *t)
{
    seq_result *result;
    stick *data;
    int id, i;
  
    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (stick *)xmalloc(sizeof(stick))))
	return -1;

    if (NULL == (data->ap_array = (a_score *)xmalloc(sizeof(a_score))))
	return -1;

    if (NULL == (data->ap_array[0].p_array = (p_score *)xmalloc(num_elements * 
							    sizeof(p_score))))
	return -1;
    
    result->data = data;
    data->n_pts = 1;
    data->ap_array[0].n_pts = num_elements;
    data->ap_array[0].dim.x0 = start;
    data->ap_array[0].dim.x1 = end;
    data->ap_array[0].dim.y0 = 0;
    data->ap_array[0].dim.y1 = MAX_TRNA_BP;

    for (i = 0; i < num_elements; i++) {
	data->ap_array[0].p_array[i].pos = results[i]->aa_left+1;
	data->ap_array[0].p_array[i].score = results[i]->total_bp_score;
    }

    id = get_reg_id();

    result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    result->seq_id[VERTICAL] = -1;

    result->id = id; 
    result->input = (void *)input; 
    result->text_data = results;
    result->type = SEQ_TYPE_TRNA;
    result->frame = 0;
    result->graph = SEQ_STICK;

    result->pr_func = stick_plot_func;
    result->op_func = trna_search_callback;
    result->txt_func = trna_search_text_func;
   
    seq_register(seq_num, trna_search_callback, (void *)result, 
		 SEQ_PLOT_PERM, id);
    return id;
}

int init_nip_trna_search_create(Tcl_Interp *interp,
				int seq_id,
				int start,
				int end,
				int *id)
{
    in_trna_search *input;
    char *seq;
    int seq_len;
    int i;
    int seq_num;
    TrnaRes **results;
    int nmatch;
    TrnaSpec *t;
    int max_score = 0;
    Tcl_DString input_params;
    
    vfuncheader("trna search");
    set_char_set(DNA);

    if (NULL == (input = (in_trna_search *)xmalloc
		 (sizeof(in_trna_search))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n",
		       GetSeqName(seq_num), start, end);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (NULL == (results = (TrnaRes **)xmalloc(MAX_TRNA * sizeof(TrnaRes *))))
	return -1;

    trna_search(seq, seq_len, start, end, &results, &nmatch, 
		&max_score, &t);

    if (nmatch == 0) {
	verror(ERR_WARN, "trna search", "no matches found");
	
	for (i = 0; i < MAX_TRNA; i++) {
	    xfree(results[i]);
	}
	xfree(results);
	xfree(t);
	xfree(input->params);
	xfree(input);
	return -1;
    }

    input->t = t;

    *id = store_trna_search(seq_num, input, start, end, results, nmatch, t);

    for (i = 0; i < nmatch; i++) {
	if ( results[i]->total_cb_score >= t->min_total_cb_score ) 
	    draw_trna(results[i]);
    }

    return 0;
}

int init_nip_trna_search_plot(Tcl_Interp *interp,
			      int seq_id,
			      int result_id,
			      char *raster_win,
			      int raster_id,
			      char *colour,
			      int line_width,
			      int tick_ht)
{
    config *configure;

    if (NULL == (configure = (config *)xmalloc(sizeof(config))))
	return -1;

    configure->position = 0.5;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = tick_ht;
    configure->zoom = 1;
    configure->scroll = 0;

    init_stick_raster(interp, seq_id, result_id, raster_win, raster_id, 
		      configure, colour, line_width, tick_ht);

    return 0;
}
