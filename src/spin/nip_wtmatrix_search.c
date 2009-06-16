#include <tcl.h>
#include <string.h>

#include "seq_raster.h"
#include "xalloc.h"
#include "seq_reg.h"
#include "seq_results.h"
#include "nip_results.h"
#include "tkRaster.h"
#include "text_output.h"
#include "nip_globals.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "nip_raster.h"
#include "dna_utils.h"
#include "splice_search.h"
#include "nip_wtmatrix_search.h"
#include "misc.h"
#include "seq_plot_funcs.h"

void nip_wtmatrix_search_callback(int seq_num, void *obj, seq_reg_data *jdata);

void nip_wtmatrix_search_shutdown(Tcl_Interp *interp,
				  seq_result *result,
				  char *raster_win,
				  int seq_num)
{
    in_wtmatrix_search *input = result->input;
    stick *data = result->data;
    out_raster *output = result->output;
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

    seq_deregister(seq_num, nip_wtmatrix_search_callback, 
		   (seq_result *)result);
		
    /* 
     * only bother replotting the raster if there are still results in the
     * raster
     */
    if (raster_result && raster_result->num_results > 1) {
	ReplotAllCurrentZoom(interp, raster_win);
	tmp = get_default_string(interp, tk_utils_defs, 
				 w("RASTER.RESULTS.WIN"));
	if (TCL_OK != Tcl_VarEval(interp, "seq_result_list_update ", 
				  tmp, NULL)){
	    puts(Tcl_GetStringResult(interp));
	}
	
	if (TCL_OK != Tcl_VarEval(interp, "RemoveRasterResultKey ", raster_win,
				  " {", info.line, "}", NULL))
	    verror(ERR_WARN, "wtmatrix_search", "shutdown %s \n", Tcl_GetStringResult(interp));
    }
    xfree(data->ap_array[0].p_array);
    xfree(data->ap_array);
    xfree(result->data);

    free(input->params);
    xfree(result->input);
    xfree(result->text_data);
    xfree(output->configure[0]);
    xfree(output->configure);
    xfree(result->output);

    xfree(result);

    if (raster_result) 
	DeleteResultFromRaster(raster_result);    
}

void nip_wtmatrix_search_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_wtmatrix_search *input = result->input;
    out_raster *output = result->output;
    stick *data = result->data;
    int id = result->id;
    char cmd[1024];

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "wtmatrix search");
	break;
	
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "wtmatrix #%d", result->id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "wtmatrix: seq=%s", 
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
		nip_wtmatrix_search_shutdown(interp, result, 
					   output->raster_win, seq_num);
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
					      w("NIP.WTMATRIX_SEARCH.PLOT_HEIGHT"));
		    
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
    case SEQ_DELETE: {
	Tcl_Interp *interp = output->interp;
	nip_wtmatrix_search_shutdown(interp, result, output->raster_win, 
				   seq_num);
	
	break;
    }
    }
}

void
nip_wtmatrix_search_text_func(void *obj)
{
    int i;
    seq_result *result = (seq_result *) obj;
    stick *data = result->data;
    text_wtmatrix *text_data = result->text_data;
    char *seq_name;
    char *sequence;
    int seq_num;
    seq_num = GetSeqNum(result->seq_id[0]);
    seq_name = GetSeqName(seq_num);
    sequence = GetSeqSequence(seq_num);

    for (i = 0; i < data->ap_array[0].n_pts; i++) {
	UpdateTextOutput();

	vmessage("Position %8d %8d score %f %.*s\n", 
		 data->ap_array[0].p_array[i].pos - text_data->mark_pos, 
		 data->ap_array[0].p_array[i].pos + 1, 
		 data->ap_array[0].p_array[i].score, 
		 text_data->length, 
		 &sequence[data->ap_array[0].p_array[i].pos - 1 - 
			  text_data->mark_pos]);
	/*
	vmessage("Position %8d %8d score %f %.*s\n", 
		 data->match[i]->pos+1-data->mark_pos, data->match[i]->pos+2, 
		 data->match[i]->score, data->length, 
		 &sequence[data->match[i]->pos-data->mark_pos]);
	*/
    }
}

int store_wtmatrix_search(int seq_num,
			  in_wtmatrix_search *input,
			  int start,
			  int end,
			  WtmatrixRes *results)
{
    seq_result *result;
    int id, i;
    stick *data;
    text_wtmatrix *text_data;

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (stick *)xmalloc(sizeof(stick))))
	return -1;
    
    if (NULL == (data->ap_array = (a_score *)xmalloc(sizeof(a_score))))
	return -1;

    if (NULL == (data->ap_array[0].p_array = (p_score *)xmalloc(results->number_of_res * 
							    sizeof(p_score))))
	return -1;
    
    if (NULL == (text_data = (text_wtmatrix *)xmalloc(sizeof(text_wtmatrix))))
	return -1;
    
    result->data = data;
    data->n_pts = 1;
    data->ap_array[0].n_pts = results->number_of_res;
    data->ap_array[0].dim.x0 = start;
    data->ap_array[0].dim.x1 = end;
    data->ap_array[0].dim.y0 = results->min;
    data->ap_array[0].dim.y1 = results->max;

    for (i = 0; i < results->number_of_res; i++) {
	data->ap_array[0].p_array[i].pos = results->match[i]->pos + 1;
	data->ap_array[0].p_array[i].score = results->match[i]->score;
    }

    result->text_data = text_data; 
    text_data->length = results->length;
    text_data->mark_pos = results->mark_pos;
    
    id = get_reg_id();
    result->id = id; 
    result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    result->seq_id[VERTICAL] = -1;

    result->input = (void *)input; 
    result->output = NULL;
    result->type = SEQ_TYPE_WTMATRIXSEARCH;
    result->frame = 0;
    result->graph = SEQ_STICK;

    result->pr_func = stick_plot_func;
    result->op_func = nip_wtmatrix_search_callback;
    result->txt_func = nip_wtmatrix_search_text_func;
    
    seq_register(seq_num, nip_wtmatrix_search_callback, (void *)result, 
		 SEQ_PLOT_PERM, id);

    if (results) free_WtmatrixRes (results);
    return id;
}

int init_nip_wtmatrix_search_create(int start,
				    int end,
				    int seq_id,
				    char *wt_matrix,
				    int *id)
{
    in_wtmatrix_search *input;
    char *seq;
    int seq_len;
    int seq_num;
    WtmatrixRes *results = NULL;
    Tcl_DString input_params;
    int irs;

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

    irs = weight_search(seq, seq_len, start, end, wt_matrix, &results);
    if (irs == -1) {
	verror(ERR_WARN, "weight matrix search", "error in weight matrix search");
	return -1;
    }
    if (results->number_of_res == 0) {
	verror(ERR_WARN, "weight matrix search", "no matches found");
	return -1;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "weight matrix %s\n",
		       GetSeqName(seq_num), start, end, wt_matrix);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);
    
    if (-1 == (*id = store_wtmatrix_search(seq_num, input, start, end, results))){
	verror(ERR_FATAL,"weight matrix search", "error in saving matches\n");
	return -1;
    }

  return 0;
}

int init_nip_wtmatrix_search_plot(Tcl_Interp *interp, 
				  int seq_id,  
				  int result_id,
				  char *raster_win, 
				  int raster_id,
				  char *colour,
				  int line_width,
				  int tick_ht)
{
    config *configure;
    
    /* no results found */
    if (result_id == -1) {
	return 0;
    }
    
    if (NULL == (configure = (config *)xmalloc(sizeof(config))))
	return -1;

    configure->position = 0.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = tick_ht;
    configure->zoom = 1;
    configure->scroll = 0;

    init_stick_raster(interp, seq_id, result_id, raster_win, raster_id, 
		      configure, colour, line_width, tick_ht);

    return 0;
}


