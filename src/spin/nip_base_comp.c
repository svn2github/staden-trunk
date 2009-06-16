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
#include "tcl_utils.h"
#include "nip_raster.h"
#include "dna_utils.h"
#include "nip_base_comp.h"
#include "base_comp.h"
#include "seq_plot_funcs.h"

void plot_base_comp_callback(int seq_num, void *obj, seq_reg_data *jdata);

void plot_base_comp_shutdown(Tcl_Interp *interp,
			     seq_result *result,
			     char *raster_win,
			     int seq_num)
{
    in_plot_base_comp *input = result->input;
    r_graph *data = result->data;
    out_raster *output = result->output;
    char *tmp;
    seq_reg_key_name info;
    static char buf[80];
    int raster_id;
    double wx0, wy0, wx1, wy1;
    Tk_Raster *raster;
    Tcl_CmdInfo info1;
    RasterResult *raster_result;

    /* determine raster_id and raster_result structure */
    Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));
    raster_result = raster_id_to_result(raster_id);

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(result->id, (seq_reg_data *)&info, 0);

    seq_deregister(seq_num, plot_base_comp_callback, (seq_result *)result);

    /* 
     * only bother replotting the raster if there are still results in the
     * raster
     */
    if (raster_result && raster_result->num_results > 1) {
	ReplotAllCurrentZoom(interp, raster_win);
	tmp = get_default_string(interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
    
	if (TCL_OK != Tcl_VarEval(interp, "seq_result_list_update ", 
				  tmp, NULL)){
	    verror(ERR_WARN, "base composition", "base_comp shutdown %s \n", 
		   Tcl_GetStringResult(interp));
	}

	if (TCL_OK != Tcl_VarEval(interp, "RemoveRasterResultKey ", raster_win,
				  " {", info.line, "}", NULL))
	    verror(ERR_WARN, "base composition", "base_comp remove %s \n", 
		   Tcl_GetStringResult(interp));
	
	Tcl_GetCommandInfo(interp, raster_win, &info1);
	raster = (Tk_Raster*)info1.clientData;

	/* find original y before reset size */
	RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

	/* update the scale of the remaining results in the raster window */
	SeqReSetRasterWindowSize(interp, raster_win, result->graph);
	ReSetRasterWindowWorld(interp, raster_win, wy1, result->graph);
	ReplotAllRasterWindow(interp, raster_win);
    }

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

void plot_base_comp_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_plot_base_comp *input = result->input;
    out_raster *output = result->output;
    r_graph *data = result->data;
    int id = result->id;
    char cmd[1024];

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Plot base composition");
	break;    

    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "base comp #%d", result->id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "base comp: seq=%s", 
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
		plot_base_comp_shutdown(interp, result, output->raster_win, 
					seq_num);
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
	    } /* WIN_SIZE */
	} /* switch SEQ_RESULT_INFO */
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
#ifdef DEBUG
	    printf("BASECOMP delete\n");
#endif
	    plot_base_comp_shutdown(interp, result, output->raster_win,
				    seq_num);
	    break;
	} /* SEQ_DELETE, SEQ_QUIT */
    } /* switch */
}


void
plot_base_comp_text_func(void *obj)
{
    seq_result *result = (seq_result *) obj;
    r_graph *data = result->data;
    int n_pts = data->n_pts;
    int i;

    for (i = 0; i < n_pts; i++) {
	UpdateTextOutput();
	vmessage("Position %10d score %10d\n", 
		 data->p_array[i].pos, (int)data->p_array[i].score);
    }
}

int store_base_comp(int seq_num,
		    int window_length,
		    in_plot_base_comp *input,
		    double *match, 
		    int start,
		    int end,
		    int num_elements,
		    double min_match,
		    double max_match)
{
    int i;
    seq_result *result;
    r_graph *data;
    int id;
    out_raster *output;

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (r_graph *)xmalloc(sizeof(r_graph))))
	return -1;

    if (NULL == (data->p_array = (p_score *)xmalloc(sizeof(p_score) * 
						    num_elements)))
	 return -1;

    result->data = data;

    if (NULL == (output = (out_raster *)xmalloc(sizeof(out_raster))))
	return -1;

    id = get_reg_id();

    /* data assignment */
    for (i = 0; i < num_elements; i++) {
	data->p_array[i].pos = i + start;
	data->p_array[i].score = match[i];
    }
    data->n_pts = num_elements;
    data->dim.x0 = start;
    data->dim.x1 = end;
    data->dim.y0 = min_match;
    data->dim.y1 = max_match;

    result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    result->seq_id[VERTICAL] = -1;

    result->id = id; 
    result->input = (void *)input; 
    result->output = (void *)output;
    result->type = SEQ_TYPE_BASECOMP;
    result->frame = 0;
    result->graph = SEQ_GRAPH;

    result->pr_func = graph_plot_func;
    result->op_func = plot_base_comp_callback;
    result->txt_func =plot_base_comp_text_func;

    seq_register(seq_num, plot_base_comp_callback, (void *)result, 
		 SEQ_PLOT_PERM, id);
    return id;
}

int init_nip_base_comp_create(Tcl_Interp *interp, 
			      int seq_id, 
			      int start, 
			      int end, 
			      int win_len, 
			      int a,
			      int c, 
			      int g, 
			      int t, 
			      int *id)
{
    in_plot_base_comp *input;
    char *seq;
    int seq_len;
    double *match;
    int i;
    double score[5];
    int seq_num;
    int num_char;
    double min_match, max_match;
    Tcl_DString input_params;
    int irs;

    vfuncheader("plot base composition");
    set_char_set(DNA);

    if (NULL == (input = (in_plot_base_comp *)xmalloc
		 (sizeof(in_plot_base_comp))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }

    num_char = 5;

    for (i = 0; i < num_char; i++) {
	score[i] = 0.0;
    }
    if (a) 
	score[char_lookup['a']] = 1.0;
    if (c) 
	score[char_lookup['c']] = 1.0;
    if (g) 
	score[char_lookup['g']] = 1.0;
    if (t) 
	score[char_lookup['t']] = 1.0;

    if (NULL == (match = (double *)xmalloc((seq_len+2) * sizeof(double))))
	return -1;

    /* Rodger's version */
    irs = get_base_comp_res(seq, seq_len, win_len, start, end, score, 
			    match, &min_match, &max_match);

    if (irs == -1 || (min_match == 0 && max_match == 0)) {

	verror(ERR_WARN,"plot base composition", 
	       "error in calculating the base composition \n");
	xfree(input);
	xfree(match);
	return -1;
    }

    /* My version */
    /*
      Plot_Base_Comp(win_len, score, seq, seq_len, start, end,
      match, &max_match);
      */

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "window length %d A %d C %d G %d T %d\n",
		       GetSeqName(seq_num), start, end,
		       win_len, a, c, g, t);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (-1 == (*id = store_base_comp(seq_num, win_len, 
				     input, match, start, end,
				     end - start + 1, min_match, max_match))){
      verror(ERR_FATAL,"base composition", "error in saving matches\n");
      return -1;
    }

    xfree(match);
    return 0;
}

int init_nip_base_comp_plot(Tcl_Interp *interp, 
			    int seq_id, 
			    int result_id,
			    char *raster_win,
			    int raster_id,
			    char *colour, 
			    int line_width) 
{
    config *configure;
    if (NULL == (configure = (config *)xmalloc(sizeof(config))))
	return -1;

    configure->position = 0.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = 1.0;
    configure->zoom = 2;
    configure->scroll = 1;

    init_graph_raster(interp, seq_id, result_id, raster_win, raster_id,
		      configure, colour, line_width);
    return 0;
}



