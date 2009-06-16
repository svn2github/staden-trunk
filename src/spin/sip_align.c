#include <string.h>

#include "seq_results.h"
#include "seq_raster.h"
#include "sip_results.h"
#include "seq_reg.h"
#include "sip_align.h"
#include "misc.h" /* MIN */
#include "text_output.h"  /* UpdateTextOutput */
#include "tcl_utils.h"
#include "sip_globals.h"
#include "sequence_formats.h"
#include "dna_utils.h"
#include "readpam.h"
#include "seq_plot_funcs.h"
#include "sequence_pair_display.h"

static char *aa_order = "ABCDEFGHIKLMNPQRSTVWYZX*-";
static char *nt_order = "ACGT";

void align_callback(int seq_num, void *obj, seq_reg_data *jdata);

void align_shutdown(Tcl_Interp *interp,
		    int seq_num,
		    seq_result *result,
		    char *raster_win,
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
    in_align *input = result->input;

    /* determine raster_id and raster_result structure */
    Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));
    raster_result = raster_id_to_result(raster_id);

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(result->id, (seq_reg_data *)&info, 0);

    /* deregister horizontal sequence */
    seq_deregister(GetSeqNum(result->seq_id[HORIZONTAL]), align_callback, 
		   (seq_result *)result);

    /* deregister vertical sequence */
    seq_deregister(GetSeqNum(result->seq_id[VERTICAL]), align_callback, 
		   (seq_result *)result);
    
    /* 
     * only bother replotting the raster if there are still results in the
     * raster
     */
    if (raster_result && raster_result->num_results > 1) {

	ReplotAllCurrentZoom(interp, raster_win);

	Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
	raster_id = atoi(Tcl_GetStringResult(interp));
    
	if (TCL_OK != Tcl_VarEval(interp, "RemoveRasterResultKey ", raster_win,
				  " {", info.line, "}", NULL))
	    verror(ERR_WARN, "align_shutdown", "%s \n", Tcl_GetStringResult(interp));
	
	Tcl_GetCommandInfo(interp, raster_win, &info1);
	raster = (Tk_Raster*)info1.clientData;

	/* find original y before reset size */
	RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

	/* update the scale of the remaining results in the raster window */
	SeqReSetRasterWindowSize(interp, raster_win, result->graph);
	ReSetRasterWindowWorld(interp, raster_win, wy1, result->graph);

	ReplotAllRasterWindow(interp, raster_win);
	tmp = get_default_string(interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
	if (TCL_OK != Tcl_VarEval(interp, "seq_result_list_update ", tmp, NULL)){
	    verror(ERR_WARN, "align_shutdown", "%s\n", Tcl_GetStringResult(interp));
	}
    }
    DestroySequencePairDisplay(interp, id);
    free(input->params);
    SipFreeResult(result);

    if (raster_result) 
	DeleteResultFromRaster(raster_result);
}

void align_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_align *input = result->input;
    out_raster *output = result->output;
    d_plot *data = result->data;
    int id = result->id;
    char cmd[1024];

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Align sequences");
	break;
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "align #%d", result->id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "align: hori=%s vert=%s", 
		GetSeqBaseName(GetSeqNum(result->seq_id[HORIZONTAL])),
		GetSeqBaseName(GetSeqNum(result->seq_id[VERTICAL])));
	break;
    case SEQ_GET_OPS:
	if (output->hidden) {
	    jdata->get_ops.ops = "Information\0List results\0PLACEHOLDER\0"
		"PLACEHOLDER\0PLACEHOLDER\0Reveal\0SEPARATOR\0Remove\0";
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
	case 1: /* results */
	    Tcl_Eval(output->interp, "SetBusy");
	    vfuncheader("results");
	    /* result->txt_func(result); */
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
	      align_shutdown(interp, seq_num, result, output->raster_win, id);
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
    case SEQ_DELETE: {
	Tcl_Interp *interp = output->interp;
	align_shutdown(interp, seq_num, result, output->raster_win, id);
    }
    }
}

void align_text_func(void *obj)
{
#ifdef REMOVE
    seq_result *result = (seq_result *) obj;
    in_align *input = result->input;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    char *r_seq1, *r_seq2;
    int pos1, pos2;
    int seq_num_h, seq_num_v;

    seq_num_h = GetSeqNum(result->seq_id[HORIZONTAL]);
    seq_num_v = GetSeqNum(result->seq_id[VERTICAL]);

    seq1 = GetSeqSequence(seq_num_h);
    seq1_len = GetSeqLength(seq_num_h);
    seq2 = GetSeqSequence(seq_num_v);
    seq2_len = GetSeqLength(seq_num_v);

    if (NULL == (r_seq1 = (char *)xcalloc(((seq1_len+seq2_len)*2+1),
					  sizeof(char)))) {
	return;
    }
    if (NULL == (r_seq2 = (char *)xcalloc(((seq1_len+seq2_len)*2+1),
					  sizeof(char)))) {
	return;
    }
/*
    {
        int r_len1, r_len2;
        cexpand(seq1, seq2, seq1_len, seq2_len, r_seq1, r_seq2,
	        &r_len1, &r_len2, ALIGN_J_SSH | ALIGN_J_PADS, res);
    }
*/
    pos1 = input->min1;
    pos2 = input->min2;
/*
    spin_list_alignment(r_seq1, r_seq2, name1, name2, pos1, pos2, "");
*/    
    xfree(r_seq1);
    xfree(r_seq2);
#endif
}

int store_align(int seq1_num,
		int seq2_num,
		int start_h,
		int end_h,
		int start_v,
		int end_v,
		int seq1_len,
		int seq2_len,
		in_align *input,
		align_int *res)
{
    seq_result *result;
    d_plot *data;
    int i, j, k, d, id, x, y;
    int cnt = 0;
    int num_elements = seq1_len + seq2_len + 1;
    int n_pts = 0;
    
    /*
    for (i = 0; i < (seq1_len + seq2_len + 1); i++) {
	printf("res[%d]=%d\n", i, res[i]);
    }
   */

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (d_plot *)xmalloc(sizeof(d_plot))))
	return -1;
    
    if (NULL == (data->p_array = (pt_score *)xmalloc(sizeof(pt_score) * 
						     num_elements)))
	return -1;

    /* data assignment */
    i = 0;
    j = 0;
    k = 0;
    d = 0;
    x = y = 0;

    i = start_h;
    j = start_v;

    data->p_array[k].x = i;
    data->p_array[k].y = j;
    while (i <= end_h || j <= end_v) { 
	if (res[cnt] == 0) {
	    i++;
	    j++;
	    if (cnt == 0)
		k++;

	    if (cnt && res[cnt-1])
		k++;
	    data->p_array[k].x = i;
	    data->p_array[k].y = j;
	} else if (res[cnt] < 0) {
	    i -= res[cnt];
	    data->p_array[++k].x = i;
	    data->p_array[k].y = j;
	} else if (res[cnt] > 0) {
	    j += res[cnt];
	    data->p_array[++k].x = i;
	    data->p_array[k].y = j;
	}
	cnt++;
    }
    k++;

    n_pts = k;

#ifdef DEBUG
    for (i = 0; i < n_pts; i++) {
	printf("i %d x %d y %d\n", i, data->p_array[i].x, data->p_array[i].y);
    }
#endif
    id = get_reg_id();
    result->data = data;
    
    data->n_pts = n_pts;
    data->dim.x0 = start_h;
    data->dim.x1 = end_h;
    data->dim.y0 = start_v;
    data->dim.y1 = end_v;

    /* need to initialise this here */
    result->output = NULL;

    result->seq_id[HORIZONTAL] = GetSeqId( seq1_num);
    result->seq_id[VERTICAL] = GetSeqId(seq2_num);

    result->input = (void *)input; 
    result->id = id; 
    result->graph = SEQ_DOT;

    result->pr_func = dot_plot_line_func;
    result->op_func = align_callback;
    result->txt_func = align_text_func;

    seq_register(seq1_num, align_callback, (void *)result, SEQ_PLOT_PERM, id);
    seq_register(seq2_num, align_callback, (void *)result, SEQ_PLOT_PERM, id);

    return id;
}

/*
 * save default name of alignment based on the parent name with _uid appended
 */
int
SipSaveAlignment(Tcl_Interp *interp,
		 int seq_num,
		 char *sequence,
		 char *name)
{
    static int id = 0;
    char *name1;
    char *seq;

    /* 12/1/99 johnt - made SRS support conditional - define USE_SRS for SRS support */
    int seq_len = strlen(sequence);

    if (NULL == (name1 = (char *)xmalloc((strlen(name)+10)*sizeof(char)))) {
	return -1;
    }
    if (NULL == (seq = (char *)xmalloc((seq_len + 1)*sizeof(char)))) {
	return -1;
    }

    memcpy(seq, sequence, seq_len);
    seq[seq_len] = '\0';
    sprintf(name1, "%s_a%d", name, id);
    if (AddSequence(interp, -1, GetSeqLibrary(seq_num), name1, seq, 
		    GetSeqStructure(seq_num), GetSeqType(seq_num), NULL, " ") == -1)
	return -1;

    id++;
    xfree(name1);
    return 0;
}


int init_sip_global_align_create(Tcl_Interp *interp, 
				 int seq_id_h,
				 int seq_id_v, 
				 int start_h,
				 int end_h, 
				 int start_v,
				 int end_v, 
				 int match,
				 int mis_match,
				 int start_gap,
				 int cont_gap,
				 int *id)
{
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    align_int *res;
    int r_len1, r_len2;
    char *r_seq1, *r_seq2;
    char *name1, *name2;
    int i, j;
    int seq1_type, seq2_type;
    int seqtype = 0;
    char *order;
    int seq1_num;
    int seq2_num;
    in_align *input;
    int **matrix = NULL;
    Tcl_DString input_params;

    vfuncheader("align sequences");

    if (-1 == (seq1_num = GetSeqNum(seq_id_h))) {
	verror(ERR_WARN, "align sequences", 
	       "horizontal sequence undefined");
	goto error;
    }
     
    if (-1 == (seq2_num = GetSeqNum(seq_id_v))) {
	verror(ERR_WARN, "align sequences", 
	       "vertical sequence undefined");
	goto error;
    }

     /* only align dna or protein */
    seq1_type = GetSeqType(seq1_num);
    seq2_type = GetSeqType(seq2_num);

    if (seq1_type != seq2_type) {
	verror(ERR_FATAL, "align sequences", "sequences must both be either DNA or protein");
	return TCL_OK;
    } else if (seq1_type == PROTEIN) {
	set_char_set(PROTEIN);
        set_score_matrix(get_matrix_file(PROTEIN));
	order = aa_order;
    } else if (seq1_type == DNA) {

	if (NULL == (matrix = (int **)xmalloc(MAX_SCORE_MATRIX * sizeof(int*))))
	    return TCL_OK;
	for (i = 0; i < MAX_SCORE_MATRIX; i++) {
	    if (NULL == (matrix[i] = (int *)xmalloc(MAX_SCORE_MATRIX * sizeof(int))))
		return TCL_OK;
	}
	set_char_set(DNA);

	for (i = 0; i < 5; i++) {
	    for (j = 0; j < 5; j++) {
		if (i == j && i < 4) {
		    matrix[i][j] = match;
		} else {
		    matrix[i][j] = mis_match;
		}
	    }
	}
        set_score_matrix(matrix);
	order = nt_order;
    }

    if (seq1_type == PROTEIN) {
	seqtype = 1;
    }
    seq1 = GetSeqSequence(seq1_num);
    if ((seq1_len = end_h - start_h + 1) < 1) {
	verror(ERR_WARN, "align sequences", "negative length");
	return TCL_OK;
    }
    seq2 = GetSeqSequence(seq2_num);
    if ((seq2_len = end_v - start_v + 1) < 1) {
	verror(ERR_WARN, "align sequences", "negative length");
	return TCL_OK;
    }

    if (NULL == (res = (align_int *)xcalloc((seq1_len+seq2_len+1),
					    sizeof(align_int)))) {
	return TCL_OK;
    }
     if (NULL == (r_seq1 = (char *)xcalloc(((seq1_len+seq2_len)*2+1),
					  sizeof(char)))) {
	return TCL_OK;
    }
    if (NULL == (r_seq2 = (char *)xcalloc(((seq1_len+seq2_len)*2+1),
					  sizeof(char)))) {
	return TCL_OK;
    }
    if (NULL == (input = (in_align *)xmalloc(sizeof(in_align))))
	return TCL_OK;

    name1 = GetSeqBaseName(seq1_num);
    name2 = GetSeqBaseName(seq2_num);

    init_W128(score_matrix, order, score_matrix[char_lookup['-']][char_lookup['-']]);
   
     calign(&seq1[start_h-1], &seq2[start_v-1], seq1_len, seq2_len, 
	   NULL, NULL, NULL, NULL,
	   0, 0,                            /* low, high band */
	   start_gap, cont_gap,   /* gap_open, gap_extend */
	   3,                               /* job */
	   seqtype,                        /* is protein */
	   res);
    cexpand(&seq1[start_h-1], &seq2[start_v-1], seq1_len, seq2_len, 
	    r_seq1, r_seq2, &r_len1, &r_len2, 
	    ALIGN_J_SSH | ALIGN_J_PADS, res);
  
    Tcl_DStringInit(&input_params);
    vTcl_DStringAppend(&input_params, "horizontal %s: %s from %d to %d\n"
		       "vertical %s: %s from %d to %d\nscore for match %d\n"
		       "score for mis-match %d\npenalty for starting gap %d\n"
		       "penalty for each residue in gap %d\n", 
		       GetSeqLibraryName(seq1_num), GetSeqName(seq1_num), 
		       start_h, end_h, GetSeqLibraryName(seq2_num), 
		       GetSeqName(seq2_num), start_v, end_v, match, mis_match, 
		       start_gap, cont_gap);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (-1 == (*id = store_align(seq1_num, seq2_num, start_h,
				 end_h, start_v, end_v, seq1_len, seq2_len,
				 input, res))) {
	goto error;
    }

    SipSaveAlignment(interp, seq1_num, r_seq1, name1);
    SipSaveAlignment(interp, seq2_num, r_seq2, name2);

    i = spin_list_alignment(r_seq1, r_seq2, name1, name2, start_h, 
		       start_v, "",seq1_type);
    
    xfree(r_seq1);
    xfree(r_seq2);
    xfree(res);

    if (matrix) {
	for (i = 0; i < MAX_SCORE_MATRIX; i++) {
	    xfree(matrix[i]);
	}
	xfree(matrix);
    }

    return 0;

 error:
    return -1;
}

int init_sip_global_align_plot(Tcl_Interp *interp, 
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

    init_dot_plot(interp, seq_id_h, seq_id_v, result_id, "align", raster_win,
		  raster_id, opts, 6, LINE, data->dim);
    
    xfree(opts[1]);
    xfree(opts[3]);
    xfree(opts[5]);

    return 0;
}


