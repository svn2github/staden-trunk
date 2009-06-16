#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "seq_raster.h"
#include "nip_raster.h"
#include "seq_results.h"
#include "nip_results.h"
#include "seq_reg.h"
#include "nip_stop_codon.h"
#include "nip_globals.h"
#include "text_output.h"
#include "tclXkeylist.h"
#include "tcl_utils.h"
#include "xalloc.h"
#include "tcl_utils.h"
#include "renz_utils.h"
#include "genetic_code.h"
#include "dna_utils.h"
#include "seq_plot_funcs.h"

#define MAXMATCHES 10000

void nip_stop_codons_callback(int seq_num, void *obj, seq_reg_data *jdata);

void nip_stop_codons_shutdown(Tcl_Interp *interp,
			      seq_result *result,
			      char *raster_win,
			      int seq_num)
{
    stick *data = result->data;
    out_raster *output = result->output;
    char *tmp;
    seq_reg_key_name info;
    static char buf[80];
    int raster_id;
    RasterResult *raster_result;
    in_s_codon *input;

    input = result->input;

    /* determine raster_id and raster_result structure */
    Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));
    raster_result = raster_id_to_result(raster_id);

#ifdef DEBUG
    printf("nip_stop_codons_shutdown\n");
#endif

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(result->id, (seq_reg_data *)&info, 0);

    seq_deregister(seq_num, nip_stop_codons_callback, (seq_result *)result);

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
	verror(ERR_WARN, "stop codons", "shutdown %s \n", Tcl_GetStringResult(interp));
    }
    
    xfree(data->ap_array[0].p_array);
    if (data->n_pts == 2) 
	xfree(data->ap_array[1].p_array);
    xfree(data->ap_array);
    xfree(data);

    free(input->params);
    xfree(result->input);

    xfree(output->configure[0]);
    if (output->n_configure == 2) 
	xfree(output->configure[1]);
    xfree(output->configure);

    xfree(result->output);
    xfree(result);

    if (raster_result) 
	DeleteResultFromRaster(raster_result);
}


void nip_stop_codons_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_s_codon *input = result->input;
    out_raster *output = result->output;
    stick *data = result->data;
    int id = result->id;
    char cmd[1024];

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	if (result->type == SEQ_TYPE_STOPCODON) {
	    sprintf(jdata->name.line, "Plot stop codons");
	} else {
	    sprintf(jdata->name.line, "Plot start codons");
	}
	break;    

    case SEQ_KEY_NAME:
	if (result->type == SEQ_TYPE_STOPCODON) {
	    sprintf(jdata->name.line, "stop f%d #%d", result->frame,
		    result->id);
	} else {
	    sprintf(jdata->name.line, "start f%d #%d", result->frame,
		    result->id);
	}
	break;
	
    case SEQ_GET_BRIEF:
	if (result->type == SEQ_TYPE_STOPCODON) {
	    sprintf(jdata->name.line, "stop codons: seq=%s frame=%d", 
		    GetSeqName(GetSeqNum(result->seq_id[0])), result->frame);
	} else {
	    sprintf(jdata->name.line, "start codons: seq=%s frame=%d", 
		    GetSeqName(GetSeqNum(result->seq_id[0])), result->frame);
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
		nip_stop_codons_shutdown(interp, result, output->raster_win, 
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
					   w("RASTER.SINGLE.PLOT_HEIGHT"));

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
	    nip_stop_codons_shutdown(interp, result, output->raster_win, 
				     seq_num);
	    break;
	}
    }
}


/*static int compare_int(const void *p1, const void *p2)
{
    int *i1 = (int *)p1, *i2 = (int *)p2;
    return (*i1 - *i2);

     HACK - maybe needed for very large numbers to avoid overflow  
     if (*i1 < *i2) 
     return (-1);
     else if (*i1 == *i2) 
     return (0);
     else 
     return (1);    
}*/

static int compare_p_score(const void *p1,
			   const void *p2)
{
    p_score *i1 = (p_score *) p1, *i2 = (p_score *) p2;
    return i1->pos - i2->pos;
}


int
NipFindStopCodons(char *strand,                                        /* in */
		  char *sequence,                                      /* in */
		  int sequence_len,                                    /* in */
		  int sequence_type,                                   /* in */
		  int start,                                           /* in */
		  int end,
		  int num_codons,                                      /* in */
		  char **codon,                                        /* in */
		  s_codon_res *stop_codon)                            /* out */
{
#define CODON_LEN 3

    int *seq_hash_values;
    int *matches, max_matches;
    int num_matches;
    int i, k;
    int cnt1 = 0;
    int cnt2 = 0;
    int cnt3 = 0;
    int start_codon = 0;
    int end_codon;
    int frame;
    int last_word[SIZE_HASH];
    int word_count[SIZE_HASH];

    end_codon = num_codons - 1;
    
    if (strcmp(strand, "-") == 0) {
	start_codon = num_codons;
    }
    if ((strcmp(strand, "-") == 0) || (strcmp(strand, "both") == 0)) {
	end_codon = num_codons * 2 - 1;
    }
    
    if (NULL == (stop_codon->match1 = (int *) xmalloc (sizeof(int) * 
						       sequence_len/3 + 1)))
	return -1;
    if (NULL == (stop_codon->match2 = (int *) xmalloc (sizeof(int) * 
						       sequence_len/3 + 1)))
	return -1;
    if (NULL == (stop_codon->match3 = (int *) xmalloc (sizeof(int) * 
						       sequence_len/3 + 1)))
	return -1;

    if ( ! (seq_hash_values = (int *) xmalloc ( sizeof(int)*sequence_len))) {
	return -2;
    }

    max_matches = sequence_len+3;
    if ( ! (matches = (int *) xmalloc ( sizeof(int)*(max_matches) ))) {
	return -1;
    }
    hash_dna(&sequence[start-1], sequence_len, seq_hash_values, last_word,
	     word_count);

    for (i = start_codon; i <= end_codon; i++) {
	/* find the matches between rec seq and dna */
	/*
	error = search_dna (&sequence[start-1], sequence_len, codon[i], 
			    CODON_LEN, sequence_type, 
			    matches, max_matches, 
			    &num_matches, seq_hash_values);
			    */
	dna_search(&sequence[start-1], sequence_len, codon[i], CODON_LEN, 
		   sequence_type, seq_hash_values, last_word, word_count,
		   matches, max_matches, &num_matches);

	/* store the match results in array 'match' */
	for (k = 0; k < num_matches; k++) {

	    /* find reading frame */
	    frame = matches[k]%CODON_LEN;

	    /* set reading frame 0 to be 3 */
	    if (frame == 0) {
		frame = CODON_LEN;
	    }

	    /* 
	     * convert matches position into a absolute position of the
	     * sequence, where the first base is at position 1
	     */
	    matches[k] += (start-1);

#ifdef DEBUG
	    printf("i %d k %d matches %d frame %d \n", i, k, matches[k], frame);
#endif	    
	    if (frame == 1) {
		stop_codon->match1[cnt1] = matches[k];
		cnt1++;
	    } else if (frame == 2) {
		stop_codon->match2[cnt2] = matches[k];
		cnt2++;
	    } else if (frame == 3) {
		stop_codon->match3[cnt3] = matches[k];
		cnt3++;
	    }
	}
    } 
    stop_codon->n_match1 = cnt1;
    stop_codon->n_match2 = cnt2;
    stop_codon->n_match3 = cnt3;

    xfree(seq_hash_values);
    xfree(matches);
    return 1;
}

void nip_stop_codons_text_func(void *obj)
{
    seq_result *result = (seq_result *) obj;
    stick *data = result->data;
    int i;

    qsort((void *) data->ap_array[0].p_array, data->ap_array[0].n_pts,
	  sizeof(p_score), compare_p_score);

    for (i = 0; i < data->ap_array[0].n_pts; i++) {
	UpdateTextOutput();
	vmessage("Position %10d\n", data->ap_array[0].p_array[i].pos);
    }
    if (data->ap_array[1].n_pts > 0) {
	qsort((void *) data->ap_array[1].p_array, data->ap_array[1].n_pts,
	      sizeof(p_score), compare_p_score);
	
	for (i = 0; i < data->ap_array[1].n_pts; i++) {
	    UpdateTextOutput();
	    vmessage("Position %10d\n", data->ap_array[1].p_array[i].pos);
	}
    }
}

int
nip_stop_codons(char *sequence,
		int sequence_type,
		int start, 
		int end,
		char *strand, 
		s_codon_res *stop_codon)
{
    int max_codons = 125;
    char **codon;
    int sequence_len;
    int i, j, k, cnt = 0, num_codons;
    char bases[] = {"tcag-"};
    char (*genetic_code)[5][5] = get_global_genetic_code();

    if (NULL == (codon = (char **)xmalloc((max_codons * 2) * 
					  sizeof(char*)))) {
	return -1;
    }
    
    for (i = 0; i < max_codons; i++) {
	if (NULL == (codon[i] = (char *)xmalloc(3 * sizeof(char*)))) {
	    return -1;
	}
    }

    for (i = 0; i < 5; i++ ) {
	for (j = 0; j < 5; j++ ) {
	    for (k = 0; k < 5; k++ ) {
		if (genetic_code[i][j][k] == '*') {
		    sprintf(codon[cnt++], "%c%c%c", 
			   bases[i], bases[j], bases[k]);
#ifdef DEBUG
		    printf("i %d j %d k %d %c%c%c\n", 
			   i, j, k, bases[i], bases[j], bases[k]); 
#endif
		}
	    }
	}
    }

    num_codons = cnt;
    for (i = 0; i < num_codons; i++) {
	strcpy(codon[cnt], codon[i]);
	complement_seq(codon[cnt++], 3);
    }

#ifdef DEBUG
    printf("num_codons %d\n", num_codons);
    for (i = 0; i < cnt; i++) {
	printf("codon[%d]=%s\n", i, codon[i]);
    }
#endif

    sequence_len = end - start + 1;

    NipFindStopCodons(strand, sequence, sequence_len, sequence_type, 
		      start, end, num_codons, codon, stop_codon);


    for (i = 0; i < max_codons; i++) {
	xfree(codon[i]);
    }
    
    xfree(codon);
    return 0;
}

int nip_start_codons(char *sequence,
		     int sequence_type,
		     int start, 
		     int end,
		     char *strand, 
		     s_codon_res *start_codon)
{
    char **codon;
    int sequence_len;
    int i, j, k, cnt = 0, num_codons;
    char bases[] = {"tcag-"};
    int max_codons = 125;
    char (*genetic_code)[5][5] = get_global_genetic_code();

    if (NULL == (codon = (char **)xmalloc((max_codons * 2) * 
					  sizeof(char*)))) {
	return -1;
    }
    for (i = 0; i < max_codons; i++) {
	if (NULL == (codon[i] = (char *)xmalloc(3 * sizeof(char*)))) {
	    return -1;
	}
    }

    for (i = 0; i < 5; i++ ) {
	for (j = 0; j < 5; j++ ) {
	    for (k = 0; k < 5; k++ ) {
		if (genetic_code[i][j][k] == 'M') {
		    sprintf(codon[cnt++], "%c%c%c", 
			   bases[i], bases[j], bases[k]);
#ifdef DEBUG
		    printf("i %d j %d k %d %c%c%c\n", 
			   i, j, k, bases[i], bases[j], bases[k]); 
#endif
		}
	    }
	}
    }

    num_codons = cnt;
    for (i = 0; i < num_codons; i++) {
	strcpy(codon[cnt], codon[i]);
	complement_seq(codon[cnt++], 3);
    }

#ifdef DEBUG
    printf("num_codons %d\n", num_codons);
    for (i = 0; i < cnt; i++) {
	printf("codon[%d]=%s\n", i, codon[i]);
    }
#endif

    sequence_len = end - start + 1;

    NipFindStopCodons(strand, sequence, sequence_len, sequence_type, 
		      start, end, num_codons, codon, start_codon);
    
    for (i = 0; i < max_codons; i++) {
	xfree(codon[i]);
    }
    xfree(codon);
    return 0;
}

int InitStopCodonPlot(Tcl_Interp *interp,
		      char *raster_win,	
		      d_line dim,
		      float tick_ht,
		      char *colour,
		      int line_width,
		      int superimpose,
		      int raster_type,
		      out_raster **outputs)
{
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    out_raster *output; 
    char *opts[5];

    if (NULL == (opts[1] = (char *)xmalloc(100 * sizeof(char))))
	return -1;

    if (NULL == (opts[3] = (char *)xmalloc(5 * sizeof(char))))
	return -1;

    if (NULL == (output = (out_raster *)xmalloc(sizeof(out_raster))))
	return -1;

    if (Tcl_GetCommandInfo(interp, raster_win, &info) == 0) 
	return -1;
    raster = (Tk_Raster*)info.clientData;

    RasterInitPlotFunc(raster, SeqRasterPlotFunc);
 
    strcpy(output->raster_win, raster_win);
    output->interp = interp;
    output->hidden = 0;

    if (!superimpose) {
	RasterSetWorldScroll(raster, dim.x0, dim.y0, dim.x1, dim.y1);
	Tcl_VarEval(interp, "rasterInitZoom ", raster_win, NULL);
	SeqAddRasterToWindow(interp, raster_win, raster_type);
    } 
    opts[0] = "-fg";
    strcpy(opts[1], colour);
    opts[2] = "-linewidth";
    sprintf(opts[3], "%d", line_width);
    opts[4] = NULL;

    output->env_index = CreateDrawEnviron(interp, raster, 4, opts);

    xfree(opts[1]);
    xfree(opts[3]);
    *outputs = output;
    return 0;
}

int NipStopCodonsPlot(Tcl_Interp *interp,
		      int result_id,
		      int seq_num,
		      char *raster_win,
		      char *colour,
		      int line_width,
		      float tick_ht)
{
    out_raster *output; 
    config **configure;
    int superimpose;
    RasterResult *raster_result;
    seq_result *nip_result;
    stick *data;
    int i;

    nip_result = result_data(result_id, seq_num);
    data = nip_result->data;

    data->ap_array[0].dim.y0 = 1;
    data->ap_array[0].dim.y1 = tick_ht;
    data->ap_array[1].dim.y0 = 1;
    data->ap_array[1].dim.y1 = tick_ht;
    
    for (i = 0; i < data->ap_array[0].n_pts; i++) {
	data->ap_array[0].p_array[i].score = tick_ht;
    }    
      
    superimpose = 1;

    /* need to check if superimposing result on another plot */
    raster_result = raster_name_to_result(interp, raster_win);
    if (raster_result->num_results == 0) {
	superimpose = 0;
    }
    if (NULL == (configure = (config **)xmalloc(sizeof(config *))))
	return -1;

    if (NULL == (configure[0] = (config *)xmalloc(sizeof(config))))
	return -1;
    
    if (-1 == InitStopCodonPlot(interp, raster_win, data->ap_array[0].dim,
				tick_ht, colour, line_width, superimpose, 
				nip_result->graph, &output))
	return -1;
 
    configure[0]->position = 0.5;
    configure[0]->x_direction = '+';
    if (nip_result->type == SEQ_TYPE_STOPCODON) {
	configure[0]->y_direction = '+';
    } else {
	configure[0]->y_direction = '-';
    }
    configure[0]->height = tick_ht;
    configure[0]->zoom = 0;
    configure[0]->scroll = 0;

    nip_result->output = (void *)output;
    output->configure = configure;
    output->n_configure = 1;
    output->scroll = 'x';

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
    return 0;
}

int NipStopCodonsPlotBoth(Tcl_Interp *interp,
			  int result_id,
			  int seq_num,
			  char *raster_win,
			  char *colour,
			  int line_width,
			  float tick_ht)
{
    out_raster *output; 
    config **configure;
    int i;
    int superimpose;
    RasterResult *raster_result;
    seq_result *nip_result;
    stick *data;

    nip_result = result_data(result_id, seq_num);
    data = nip_result->data;

    data->ap_array[0].dim.y0 = 1;
    data->ap_array[0].dim.y1 = tick_ht;
    data->ap_array[1].dim.y0 = 1;
    data->ap_array[1].dim.y1 = tick_ht;
    for (i = 0; i < data->ap_array[0].n_pts; i++) {
	data->ap_array[0].p_array[i].score = tick_ht;
    }    
      
    for (i = 0; i < data->ap_array[1].n_pts; i++) {
	data->ap_array[1].p_array[i].score = tick_ht;
    }    

    superimpose = 1;

    /* need to check if superimposing result on another plot */
    raster_result = raster_name_to_result(interp, raster_win);
    if (raster_result->num_results == 0) {
	superimpose = 0;
    }
    if (NULL == (configure = (config **)xmalloc(2*sizeof(config*))))
	return -1;

    for (i = 0; i < 2; i++) {
	if (NULL == (configure[i] = (config*)xmalloc(sizeof(config))))
	    return -1;
    }

    if (-1 == InitStopCodonPlot(interp, raster_win, data->ap_array[0].dim,
				tick_ht, colour, line_width, superimpose, 
				nip_result->graph, &output))
	return -1;
  
    configure[0]->position = 0.0;
    configure[0]->x_direction = '+';
    configure[0]->y_direction = '-';
    configure[0]->height = tick_ht;
    configure[0]->zoom = 0;
    configure[0]->scroll = 0;

    configure[1]->position = 0.0;
    configure[1]->x_direction = '+';
    configure[1]->y_direction = '+';
    configure[1]->height = tick_ht;
    configure[1]->zoom = 0;
    configure[1]->scroll = 0;

    nip_result->output = (void *)output;
    output->configure = configure;
    output->n_configure = 2;
    output->scroll = 'x';

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
    return 0;
}

int store_stop_codons(int seq_num,
		      in_s_codon *input,
		      int start,
		      int end,
		      int *stop_codon,
		      int num_stop_codons,
		      int *stop_codon_rev,
		      int num_stop_codons_rev,
		      int frame,
		      int type)
{
    seq_result *result;
    stick *data;
    int id, i;

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (stick *)xmalloc(sizeof(stick))))
	return -1;

    if (NULL == (data->ap_array = (a_score *)xmalloc(2 * sizeof(a_score))))
	return -1;

    if (NULL == (data->ap_array[0].p_array = (p_score *)xmalloc(num_stop_codons * 
						       sizeof(p_score))))
	return -1;

    if (num_stop_codons_rev) {
	if (NULL == (data->ap_array[1].p_array = (p_score *)xmalloc(num_stop_codons_rev * 
								    sizeof(p_score))))
	    return -1;
    }

    result->data = data;

    if (num_stop_codons_rev == 0) {
	data->n_pts = 1;
    } else {
	data->n_pts = 2;
    }
    data->ap_array[0].n_pts = num_stop_codons;
    data->ap_array[1].n_pts = num_stop_codons_rev;
    data->ap_array[0].dim.x0 = start;
    data->ap_array[0].dim.x1 = end;
    data->ap_array[1].dim.x0 = start;
    data->ap_array[1].dim.x1 = end;
   

    for (i = 0; i < num_stop_codons; i++) {
	data->ap_array[0].p_array[i].pos = stop_codon[i];
	data->ap_array[0].p_array[i].score = 0;
    }    
      
    for (i = 0; i < num_stop_codons_rev; i++) {
	data->ap_array[1].p_array[i].pos = stop_codon_rev[i];
	data->ap_array[1].p_array[i].score = 0;
    }    

    id = get_reg_id();

    result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    result->seq_id[VERTICAL] = -1;

    result->id = id; 
    result->input = (void *)input; 
    result->output = NULL;
    result->frame = frame;
    result->pr_func = stick_plot_func;
    result->op_func = nip_stop_codons_callback;
    result->txt_func = nip_stop_codons_text_func;
    result->graph = SEQ_STICK;

    if (type == STOP_CODON) {
	result->type = SEQ_TYPE_STOPCODON;
    } else {
	result->type = SEQ_TYPE_STARTCODON;
    }

    seq_register(seq_num, nip_stop_codons_callback, (void *)result, 
		 SEQ_PLOT_PERM, id);

    xfree(stop_codon);
    if (num_stop_codons_rev > 0) 
	xfree(stop_codon_rev);

    return id;
}

int init_nip_stop_codons_create(int seq_id,
				int start,
				int end,
				char *strand_sym,
				int *id)
{
    in_s_codon *input1;
    in_s_codon *input2;
    in_s_codon *input3;
    char *seq;
    int seq_len;
    int sequence_type;
    int seq_num;
    s_codon_res *stop_codon, *stop_codon_rev;
    int type = STOP_CODON;
    Tcl_DString input_params;
    char strand[8];
    /*  seq_reg_delete del;*/

    vfuncheader("plot stop codons");

    if (NULL == (input1 = (in_s_codon *)xmalloc
		 (sizeof(in_s_codon))))
	return -1;
    if (NULL == (input2 = (in_s_codon *)xmalloc
		 (sizeof(in_s_codon))))
	return -1;
    if (NULL == (input3 = (in_s_codon *)xmalloc
		 (sizeof(in_s_codon))))
	return -1;

    if (NULL == (stop_codon = (s_codon_res *)xmalloc(sizeof(s_codon_res))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    sequence_type = GetSeqStructure(seq_num);

    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }

    if (strcmp(strand_sym, "both") == 0) {

	if (NULL == (stop_codon_rev = (s_codon_res *)xmalloc(sizeof(s_codon_res))))
	    return -1;

	nip_stop_codons(seq, sequence_type, start, end, "+",
			stop_codon);
	nip_stop_codons(seq, sequence_type, start, end, "-",
			stop_codon_rev);
	
    } else {
	nip_stop_codons(seq, sequence_type, start, end, strand_sym,
			stop_codon);
	stop_codon_rev = NULL;
    }

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);

    if (strcmp(strand_sym, "+") == 0) {
        strcpy(strand, "forward");
    } else if (strcmp(strand_sym, "-") == 0) {
        strcpy(strand, "reverse");
    } else {
        strcpy(strand, "both");
    }
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "strand %s\n",
		       GetSeqName(seq_num), start, end, strand);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input1->params = strdup(Tcl_DStringValue(&input_params)); 
    input2->params = strdup(Tcl_DStringValue(&input_params)); 
    input3->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    /* remove results if no matches */
    if (!stop_codon_rev) {
	if (stop_codon->n_match1 > 0) {
	    if (-1 == (id[0] = store_stop_codons(seq_num, input1, start, end,
						 stop_codon->match1, 
						 stop_codon->n_match1,
						 NULL, 0, 1, type))){
		verror(ERR_FATAL,"nip stop codons", "error in saving matches\n");
		return -1;
	    }
	} else {
	    id[0] = -1;
	    free(input1->params);
	    xfree(input1);
	    xfree(stop_codon->match1);
	    verror(ERR_WARN, "nip stop codons", "no matches found for frame 1");
	}

	if (stop_codon->n_match2 > 0) {
	    if (-1 == (id[1] = store_stop_codons(seq_num, input2, start, end,
						 stop_codon->match2, 
						 stop_codon->n_match2,
						 NULL, 0, 2, type))){
		verror(ERR_FATAL,"nip stop codons", "error in saving matches\n");
		return -1;
	    }
	} else {
	    id[1] = -1;
	    free(input2->params);
	    xfree(input2);
	    xfree(stop_codon->match2);
	    verror(ERR_WARN, "nip stop codons", "no matches found for frame 2");
	}

	if (stop_codon->n_match3 > 0) {
	    if (-1 == (id[2] = store_stop_codons(seq_num, input3, start, end,
						 stop_codon->match3, 
						 stop_codon->n_match3,
						 NULL, 0, 3, type))){
		verror(ERR_FATAL,"nip stop codons", "error in saving matches\n");
		return -1;
	    }
	} else {
	    id[2] = -1;
	    free(input3->params);
	    xfree(input3);
	    xfree(stop_codon->match3);
	    verror(ERR_WARN, "nip stop codons", "no matches found for frame 3");
	}
    } else {
	/* both strands */
	if (stop_codon->n_match1 > 0 || stop_codon_rev->n_match1 > 0) {
	    if (-1==(id[0]=store_stop_codons(seq_num, input1, start, end,
					     stop_codon->match1, 
					     stop_codon->n_match1,
					     stop_codon_rev->match1,
					     stop_codon_rev->n_match1, 1, type))){
		verror(ERR_FATAL,"nip stop codons", "error in saving matches\n");
		return -1;
	    }
	} else {
	    id[0] = -1;
	    free(input1->params);
	    xfree(input1);
	    xfree(stop_codon->match1);
	    xfree(stop_codon_rev->match1);
	    verror(ERR_WARN, "nip stop codons", "no matches found for frame 1");
	}
	if (stop_codon->n_match2 > 0 || stop_codon_rev->n_match2 > 0) {
	    if (-1==(id[1]=store_stop_codons(seq_num, input2, start, end,
					     stop_codon->match2, 
					     stop_codon->n_match2,
					     stop_codon_rev->match2,
					     stop_codon_rev->n_match2, 2, type))){
		verror(ERR_FATAL,"nip stop codons", "error in saving matches\n");
		return -1;
	    }
	} else {
	    id[1] = -1;
	    free(input2->params);
	    xfree(input2);
	    xfree(stop_codon->match2);
	    xfree(stop_codon_rev->match2);
	    verror(ERR_WARN, "nip stop codons", "no matches found for frame 2");
	}
	if (stop_codon->n_match3 > 0 || stop_codon_rev->n_match3 > 0) {
	    if (-1==(id[2]=store_stop_codons(seq_num, input3, start, end,
					     stop_codon->match3, 
					     stop_codon->n_match3,
					     stop_codon_rev->match3,
					     stop_codon_rev->n_match3, 3, type))){
		verror(ERR_FATAL,"nip stop codons", "error in saving matches\n");
		return -1;
	    }
	} else {
	    id[2] = -1;
	    free(input3->params);
	    xfree(input3);
	    xfree(stop_codon->match3);
	    xfree(stop_codon_rev->match3);
	    verror(ERR_WARN, "nip stop codons", "no matches found for frame 3");
	}
    }

    xfree(stop_codon);
    if (stop_codon_rev)
	xfree(stop_codon_rev);
    return 0;
}

int init_nip_start_codons_create(int seq_id,
				 int start,
				 int end,
				 char *strand_sym,
				 int *id)
{
    in_s_codon *input1;
    in_s_codon *input2;
    in_s_codon *input3;
    char *seq;
    int seq_len;
    int sequence_type;
    int seq_num;
    s_codon_res *start_codon;
    int type = START_CODON;
    Tcl_DString input_params;
    char strand[8];
    /* seq_reg_delete del;*/

    vfuncheader("plot start codons");

    if (NULL == (input1 = (in_s_codon *)xmalloc
		 (sizeof(in_s_codon))))
	return -1;
    if (NULL == (input2 = (in_s_codon *)xmalloc
		 (sizeof(in_s_codon))))
	return -1;
    if (NULL == (input3 = (in_s_codon *)xmalloc
		 (sizeof(in_s_codon))))
	return -1;

    if (NULL == (start_codon = (s_codon_res *)xmalloc(sizeof(s_codon_res))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    sequence_type = GetSeqStructure(seq_num);

    /* if the end has not been defined, set it to be the sequence length */
    if (end == -1) {
	end = seq_len;
    }

    nip_start_codons(seq, sequence_type, start, end, strand_sym, start_codon);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);

    if (strcmp(strand_sym, "+") == 0) {
        strcpy(strand, "forward");
    } else if (strcmp(strand_sym, "-") == 0) {
        strcpy(strand, "reverse");
    } else {
        strcpy(strand, "both");
    }
    vTcl_DStringAppend(&input_params, "sequence %s: from %d to %d\n"
		       "strand %s\n",
		       GetSeqName(seq_num), start, end, strand);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input1->params = strdup(Tcl_DStringValue(&input_params)); 
    input2->params = strdup(Tcl_DStringValue(&input_params)); 
    input3->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (start_codon->n_match1 > 0) {
	if (-1 == (id[0] = store_stop_codons(seq_num, input1, start, end,
					     start_codon->match1, 
					     start_codon->n_match1,
					     NULL, 0, 1, type))){
	    verror(ERR_FATAL,"nip start codons", "error in saving matches\n");
	    return -1;
	}
    } else {
	id[0] = -1;
	free(input1->params);
	xfree(input1);
	xfree(start_codon->match1);
	verror(ERR_WARN, "nip start codons", "no matches found for frame 1");
    }
    if (start_codon->n_match2 > 0) {
	if (-1 == (id[1] = store_stop_codons(seq_num, input2, start, end,
					     start_codon->match2, 
					     start_codon->n_match2,
					     NULL, 0, 2, type))){
	    verror(ERR_FATAL,"nip start codons", "error in saving matches\n");
	    return -1;
	}
    } else {
	id[1] = -1;
	free(input2->params);
	xfree(input2);
	xfree(start_codon->match2);
	verror(ERR_WARN, "nip start codons", "no matches found for frame 2");
    }
    if (start_codon->n_match3 > 0) {
	if (-1 == (id[2] = store_stop_codons(seq_num, input3, start, end,
					     start_codon->match3, 
					     start_codon->n_match3,
					     NULL, 0, 3, type))){
	    verror(ERR_FATAL,"nip start codons", "error in saving matches\n");
	    return -1;
	}
    } else {
	id[2] = -1;
	free(input3->params);
	xfree(input3);
	xfree(start_codon->match3);
	verror(ERR_WARN, "nip start codons", "no matches found for frame 3");
    }
    xfree(start_codon);
    return 0;
}

int init_nip_stop_codons_plot(Tcl_Interp *interp, 
			      char *rasters, 
			      char *raster_ids,
			      int seq_id, 
			      char *r_id,
			      char *colours,
			      int line_width,
			      int tick_ht)
{
    char *seq;
    int seq_len;
    int sequence_type;
    int seq_num;
    seq_cursor_notify cn;
    RasterResult *raster_result;
    cursor_t *cursor;
    seq_result *result;
    int num_id;
    stick *data;
    int i;
    int retval = -1;
    char **result_id = NULL;
    char **raster_win = NULL;
    char **raster_id = NULL;
    char **colour = NULL;


    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    sequence_type = GetSeqStructure(seq_num);

    if (Tcl_SplitList(interp, rasters, &num_id, &raster_win) != TCL_OK)
	goto cleanup;
    if (Tcl_SplitList(interp, raster_ids, &num_id, &raster_id) != TCL_OK)
	goto cleanup;
    if (Tcl_SplitList(interp, colours, &num_id, &colour) != TCL_OK)
	goto cleanup;
    if (Tcl_SplitList(interp, r_id, &num_id, &result_id) != TCL_OK)
	goto cleanup;

    raster_result = raster_id_to_result(atoi(raster_id[0]));
    cursor = find_raster_result_cursor(raster_result, seq_id, HORIZONTAL);
    
    result = result_data(atoi(result_id[0]), seq_num);
    data = result->data;
  
    /* move cursor to start position if no cursor is yet present */
    if (raster_result->cursor_array[cursor->id].prev_pos == -1 
	&& data->ap_array[0].dim.x0 > raster_result->cursor_array[cursor->id].prev_pos) {
	cursor->abspos = data->ap_array[0].dim.x0;
    }

    /* if only single strand */
    if (data->ap_array[1].n_pts == 0) {
	for (i = 0; i < num_id; i++) {
	    if (-1 == NipStopCodonsPlot(interp, atoi(result_id[i]), seq_num,
					raster_win[i], colour[i], 
					line_width, tick_ht)){
		verror(ERR_FATAL,"nip stop codons", "error in saving matches\n");
		goto cleanup;
	    }

	}
    } else {
	for (i = 0; i < num_id; i++) {
	    if (-1 == NipStopCodonsPlotBoth(interp, atoi(result_id[i]), seq_num,
					    raster_win[i], colour[i], 
					    line_width, tick_ht)){
		verror(ERR_FATAL,"nip stop codons", "error in saving matches\n");
		goto cleanup;
	    }
	}
    }

    /* 
     * need to ensure done all the necessary plotting before adding the 
     * cursor
     */
    Tcl_VarEval(interp, "update idletasks ", NULL);
    cn.job = SEQ_CURSOR_NOTIFY;

    for (i = 0; i < num_id; i++) {
	raster_result = raster_id_to_result(atoi(raster_id[i]));
	cursor = find_raster_result_cursor(raster_result, seq_id, HORIZONTAL);
	cn.cursor = cursor;
	cn.cursor->job = CURSOR_MOVE;
	seq_notify(seq_num, (seq_reg_data *)&cn); 
	AddResultToRaster(raster_result);
    }
    retval = 0;


cleanup:
    if(result_id)  Tcl_Free((char *)result_id);
    if(raster_win) Tcl_Free((char *)raster_win);
    if(raster_id)  Tcl_Free((char *)raster_id);
    if(colour)     Tcl_Free((char *)colour);
    return retval;
}

