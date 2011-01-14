#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "seq_raster.h"
#include "seq_reg.h"
#include "seq_results.h"
#include "nip_results.h"
#include "nip_raster.h"
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

void splice_search_callback(int seq_num, void *obj, seq_reg_data *jdata);

void splice_search_shutdown(Tcl_Interp *interp,
			    seq_result *result,
			    char *raster_win,
			    int seq_num)
{
    char *tmp;
    seq_reg_key_name info;
    static char buf[80];
    stick *data = result->data;
    out_raster *output = result->output;
    int raster_id;
    RasterResult *raster_result;
    in_splice *input;
    int i;

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
    seq_deregister(seq_num, splice_search_callback, result);
    
    /* 
     * only bother replotting the raster if there are still results in the
     * raster
     */
    if (raster_result && raster_result->num_results > 1) {
	ReplotAllCurrentZoom(interp, raster_win);
    
	tmp = get_default_string(interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
	if (TCL_OK != Tcl_VarEval(interp, "seq_result_list_update ", 
				  tmp, NULL)){
	    puts(Tcl_GetStringResult(interp));
	}
	
	if (TCL_OK != Tcl_VarEval(interp, "RemoveRasterResultKey ", raster_win,
			      " {", info.line, "}", NULL))
	    verror(ERR_WARN, "splice search", "shutdown %s \n", Tcl_GetStringResult(interp));

    }

    for (i = 0; i < data->n_pts; i++) {
	xfree(data->ap_array[i].p_array);
    }
    xfree(data->ap_array);
    xfree(result->data);

    free(input->params);
    xfree(result->input);
    xfree(output->configure[0]);
    xfree(output->configure[1]);
    xfree(output->configure);
    xfree(result->output);

    if (result->text_data) {
	if (((text_wtmatrix **)result->text_data)[0])
	    xfree(((text_wtmatrix **)result->text_data)[0]);
	if (((text_wtmatrix **)result->text_data)[1])
	    xfree(((text_wtmatrix **)result->text_data)[1]);
	xfree(result->text_data);
    }

    xfree(result);

    if (raster_result) 
	DeleteResultFromRaster(raster_result);
}

void splice_search_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_splice *input = result->input;
    out_raster *output = result->output;
    stick *data = result->data;
    int id = result->id;
    char cmd[1024];

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Splice search");
	break;    

    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "splice f%d #%d", result->frame,
		result->id);
	break;
	
    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "splice search: seq=%s frame=%d", 
		GetSeqName(GetSeqNum(result->seq_id[0])), result->frame);

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
	    vfuncheader("splice search results frame %d", result->frame);
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
		splice_search_shutdown(interp, result, output->raster_win, 
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
					   w("NIP.SPLICE.PLOT_HEIGHT"));

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
	    splice_search_shutdown(interp, result, output->raster_win,
				   seq_num);
	    break;
	}
    }
}

void
splice_search_text_func(void *obj)
{
    seq_result *result = (seq_result *) obj;
    stick *data = result->data;
    text_wtmatrix **text_data = result->text_data;
    int i;
    char *sequence;
    int seq_num;

    seq_num = GetSeqNum(result->seq_id[0]);
    sequence = GetSeqSequence(seq_num);

    vmessage("Donor\n");
    for (i = 0; i < data->ap_array[0].n_pts; i++) {
	UpdateTextOutput();
	vmessage("Position %8d %8d score %f %.*s\n", 
		 data->ap_array[0].p_array[i].pos - text_data[0]->mark_pos, 
		 data->ap_array[0].p_array[i].pos + 1, 
		 data->ap_array[0].p_array[i].score, 
		 text_data[0]->length, 
		 &sequence[data->ap_array[0].p_array[i].pos - 1 - 
			  text_data[0]->mark_pos]);
    }

    vmessage("Acceptor\n");
    for (i = 0; i < data->ap_array[1].n_pts; i++) {
	UpdateTextOutput();
	vmessage("Position %8d %8d score %f %.*s\n", 
		 data->ap_array[1].p_array[i].pos - text_data[1]->mark_pos, 
		 data->ap_array[1].p_array[i].pos + 1, 
		 data->ap_array[1].p_array[i].score, 
		 text_data[1]->length, 
		 &sequence[data->ap_array[1].p_array[i].pos - 1 - 
			  text_data[1]->mark_pos]);
    }
}

int StoreSpliceSearch(int seq_num,
		      WtmatrixRes *ied,
		      WtmatrixRes *eia,
		      in_splice *input,
		      int start,
		      int end,
		      int frame)
{
    seq_result *result;
    stick *data;
    int id, i;
    text_wtmatrix **text_data;

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (stick *)xmalloc(sizeof(stick))))
	return -1;

    if (NULL == (data->ap_array = (a_score *)xmalloc(2 * sizeof(a_score))))
	return -1;

    if (NULL == (data->ap_array[0].p_array = (p_score *)xmalloc(ied->number_of_res * 
						       sizeof(p_score))))
	return -1;
    if (NULL == (data->ap_array[1].p_array = (p_score *)xmalloc(eia->number_of_res * 
						       sizeof(p_score))))
	 return -1;

    if (NULL == (text_data = (text_wtmatrix **)xmalloc(2*sizeof(text_wtmatrix*))))
	return -1;
 
    if (NULL == (text_data[0] = (text_wtmatrix *)xmalloc(sizeof(text_wtmatrix))))
	return -1;
    if (NULL == (text_data[1] = (text_wtmatrix *)xmalloc(sizeof(text_wtmatrix))))
	return -1;
    


    result->data = data;
    data->n_pts = 2;
    data->ap_array[0].n_pts = ied->number_of_res;
    data->ap_array[1].n_pts = eia->number_of_res;
    data->ap_array[0].dim.x0 = start;
    data->ap_array[0].dim.x1 = end;
    data->ap_array[0].dim.y0 = ied->min;
    data->ap_array[0].dim.y1 = 2*ied->max;
    data->ap_array[1].dim.x0 = start;
    data->ap_array[1].dim.x1 = end;
    data->ap_array[1].dim.y0 = eia->min;
    data->ap_array[1].dim.y1 = 2*eia->max;

#ifdef DEBUG
    printf("ied %f %f eia %f %f\n", ied->min, ied->max, eia->min, eia->max);
#endif

    id = get_reg_id();
    for (i = 0; i < ied->number_of_res; i++) {
	data->ap_array[0].p_array[i].pos = ied->match[i]->pos + 1;
	data->ap_array[0].p_array[i].score = (double)ied->match[i]->score;
    }    
      
    for (i = 0; i < eia->number_of_res; i++) {
	data->ap_array[1].p_array[i].pos = eia->match[i]->pos + 1;
	data->ap_array[1].p_array[i].score = (double)eia->match[i]->score;
    }    
      
    result->text_data = text_data; 
    text_data[0]->length = ied->length;
    text_data[0]->mark_pos = ied->mark_pos;
    text_data[1]->length = eia->length;
    text_data[1]->mark_pos = eia->mark_pos;

    result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    result->seq_id[VERTICAL] = -1;

    result->id = id; 
    result->input = (void *)input; 
    result->output = NULL;
    result->type = SEQ_TYPE_SPLICE;
    result->frame = frame;
    result->graph = SEQ_STICK;

    result->pr_func = stick_pair_plot_func;
    result->op_func = splice_search_callback;
    result->txt_func = splice_search_text_func;
    
    seq_register(seq_num, splice_search_callback, (void *)result, 
		 SEQ_PLOT_PERM, id);

    if (ied) free_WtmatrixRes (ied);
    if (eia) free_WtmatrixRes (eia);
    return id;
}

int NipSpliceSearchPlot(Tcl_Interp *interp,
			int result_id,
			int seq_num,
			char *raster_win,
			char *colour,
			int line_width,
			float tick_ht,
			int frame)
{
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    char *opts[5];
    out_raster *output; 
    config *configure_d;
    config *configure_a;
    int raster_id;
    int superimpose = 1;
    RasterResult *raster_result;
    stick *data;
    seq_result *nip_result;
    
    /* no results were found */
    if (result_id == -1)
	return 0;

    nip_result = result_data(result_id, seq_num);
    data = nip_result->data;

    if (NULL == (output = (out_raster *)xmalloc(sizeof(out_raster))))
	return -1;

    if (NULL == (opts[1] = (char *)xmalloc(100 * sizeof(char))))
	return -1;

    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return -1;

    if (NULL == (output->configure = (config **)xmalloc(2*sizeof(config*))))
	return -1;

    if (NULL == (configure_a = (config *)xmalloc(sizeof(config))))
      return -1;
    if (NULL == (configure_d = (config *)xmalloc(sizeof(config))))
      return -1;
    configure_a->position = 0.0;
    configure_a->x_direction = '+';
    configure_a->y_direction = '-';
    configure_a->height = tick_ht;
    configure_a->zoom = 1;
    configure_a->scroll = 0;
    
    configure_d->position = 0.0;
    configure_d->x_direction = '+';
    configure_d->y_direction = '+';
    configure_d->height = tick_ht;
    configure_d->zoom = 1;
    configure_d->scroll = 0;

    if (Tcl_GetCommandInfo(interp, raster_win, &info) == 0) 
	return -1;
    raster = (Tk_Raster*)info.clientData;

    RasterInitPlotFunc(raster, SeqRasterPlotFunc);
 
    strcpy(output->raster_win, raster_win);
    output->interp = interp;
    output->hidden = 0;

    /* need to check if superimposing result on another plot */
    Tcl_VarEval(interp, "GetRasterId ", output->raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));
    raster_result = raster_id_to_result(raster_id);
    if (raster_result->num_results == 0) {
	superimpose = 0;
    }

    if (!superimpose) {
	RasterSetWorldScroll(raster, (double)data->ap_array[0].dim.x0, 
			     data->ap_array[0].dim.y0, 
			     (double)data->ap_array[0].dim.x1, 
			     data->ap_array[0].dim.y1);
	
	SeqAddRasterToWindow(interp, raster_win, nip_result->graph);
    }

    opts[0] = "-fg";
    strcpy(opts[1], colour);
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", line_width);
    opts[4] = NULL;

    output->env_index = CreateDrawEnviron(interp, raster, 4, opts);

    nip_result->output = (void *)output;
    output->scroll = 'x';
    output->configure[0] = configure_d;
    output->configure[1] = configure_a;
    output->sf_m = 1.0;
    output->sf_c = 0.0;

    if (superimpose) {
	SeqSuperimposeResult(interp, output->raster_win, result_id, 
			     data->ap_array[0].dim.x0, 
			     data->ap_array[0].dim.y0, 
			     data->ap_array[0].dim.x1, 
			     data->ap_array[0].dim.y1);
    }
    
    ReplotAllCurrentZoom(interp, raster_win);
    xfree(opts[1]);
    xfree(opts[3]);
    return 0;
}

int init_splice_search_create(int seq_id,
			      int start,
			      int end,
			      char *donor,
			      char *acceptor,
			      int *id) /* out */
{
    in_splice *input1;
    in_splice *input2;
    in_splice *input3;
    char *seq;
    int seq_len;
    int seq_num;
    SpliceResults *splice_result;
    Tcl_DString input_params;
    int irs;

    vfuncheader("splice search");
    set_char_set(DNA);

    if (NULL == (input1 = (in_splice *)xmalloc (sizeof(in_splice))))
	return -1;

    if (NULL == (input2 = (in_splice *)xmalloc (sizeof(in_splice))))
	return -1;

    if (NULL == (input3 = (in_splice *)xmalloc (sizeof(in_splice))))
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

    irs = splice_search(seq, seq_len, start, end, donor, acceptor, 
			splice_result);
    
    if (irs == -1) {
      xfree(splice_result);
      xfree(input1);
      xfree(input2);
      xfree(input3);
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
	return -1;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "donor weight matrix %s\nacceptor weight matrix %s\n",
		       GetSeqName(seq_num), start, end, donor, acceptor);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input1->params = strdup(Tcl_DStringValue(&input_params)); 
    input2->params = strdup(Tcl_DStringValue(&input_params)); 
    input3->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (-1 == (id[0] = StoreSpliceSearch(seq_num, splice_result->ied_f1,
					 splice_result->eia_f1, input1, 
					 start, end, 1))){
	verror(ERR_FATAL,"nip splice search", "error in saving matches\n");
	return -1;
    }

    if (-1 == (id[1] = StoreSpliceSearch(seq_num, splice_result->ied_f2, 
					 splice_result->eia_f2, input2, 
					 start, end, 2))){
	verror(ERR_FATAL,"nip splice search", "error in saving matches\n");
	return -1;
    }

    if (-1 == (id[2] = StoreSpliceSearch(seq_num, splice_result->ied_f3,
					 splice_result->eia_f3, input3, 
					 start, end, 3))){
	verror(ERR_FATAL,"nip splice search", "error in saving matches\n");
	return -1;
    }
    xfree(splice_result);
    return 0;
}

int init_splice_search_plot(Tcl_Interp *interp,
			    char *raster_win, 
			    int raster_id, 
			    char *r_id,
			    int seq_id, 
			    char *colours, 
			    int line_width, 
			    float tick_ht)
{
    int seq_num;
    int num_id;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    RasterResult *raster_result;
    cursor_t *cursor;
    seq_result *nip_result = NULL;
    seq_cursor_notify cn;
    stick *data;
    int i, cnt;
    int retval       = -1;
    char **result_id = NULL;
    char **colour    = NULL;


    seq_num = GetSeqNum(seq_id);


    if (Tcl_SplitList(interp, colours, &num_id, &colour) != TCL_OK)
	goto cleanup;
    if (Tcl_SplitList(interp, r_id, &num_id, &result_id) != TCL_OK)
	goto cleanup;


    if (Tcl_GetCommandInfo(interp, raster_win, &info) == 0) 
	goto cleanup;
    raster = (Tk_Raster*)info.clientData;

    RasterInitPlotFunc(raster, SeqRasterPlotFunc);

    raster_result = raster_id_to_result(raster_id);
    cursor = find_raster_result_cursor(raster_result, seq_id, HORIZONTAL);

    cnt = 0;
    for (i = 0; i < num_id; i++) {
	if (atoi(result_id[i]) > -1) {
	    nip_result = result_data(atoi(result_id[i]), seq_num);
	} else {
	    cnt++;
	}
    }
    /* no results were found */
    if (cnt == num_id) {
	retval = 0;
	goto cleanup;
    }
    data = nip_result->data;

    /* move cursor to start position if no cursor is yet present */
    if (raster_result->cursor_array[cursor->id].prev_pos == -1 
	&& data->ap_array[0].dim.x0 > raster_result->cursor_array[cursor->id].prev_pos) {
	cursor->abspos = data->ap_array[0].dim.x0;
    }

    for (i = 0; i < num_id; i++) {
	if (-1 == NipSpliceSearchPlot(interp, atoi(result_id[i]), 
				      seq_num, raster_win, colour[i], 
				      line_width, tick_ht, i+1)) {
	    verror(ERR_FATAL,"nip splice search", "error in saving matches\n");
	    goto cleanup;
	}
    }

    /* 
     * need to ensure done all the necessary plotting before adding the 
     * cursor
     */
    Tcl_VarEval(interp, "update idletasks ", NULL);
    cn.job = SEQ_CURSOR_NOTIFY;

    cn.cursor = cursor;
    cn.cursor->job = CURSOR_MOVE;
    seq_notify(seq_num, (seq_reg_data *)&cn); 
    raster_result = raster_id_to_result(raster_id);

    /* adding 3 results to a single raster */
    AddResultToRaster(raster_result);
    AddResultToRaster(raster_result);
    AddResultToRaster(raster_result);
    retval = 0;


cleanup:
    if(result_id) Tcl_Free((char *)result_id);
    if(colour)    Tcl_Free((char *)colour);
    return retval;
}
