#include   <stdio.h>
#include   <ctype.h>
#include <math.h>
#include <string.h>


#include "sim.h"
#include "seq_results.h"
#include "seq_reg.h"
#include "seq_raster.h"
#include "sip_results.h"
#include "sip_sim.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "sip_globals.h"
#include "align.h"
#include "sequence_formats.h"
#include "readpam.h"
#include "dna_utils.h"
#include "seq_plot_funcs.h"
#include "sequence_pair_display.h"

#define ROUND(x)	((x) < 0.0 ? ceil((x) - 0.5) : floor((x) + 0.5))

#define DEFAULT_M 1.0
#define DEFAULT_I -1.0
#define DEFAULT_V -1.0
#define DEFAULT_O 6.0
#define DEFAULT_E 0.2
#define DEFAULT_PAM_O 12.0
#define DEFAULT_PAM_E 4.0

void sim_callback(int seq_num, void *obj, seq_reg_data *jdata);

void sim_shutdown(Tcl_Interp *interp,
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
    in_sim *input = result->input;

    /* determine raster_id and raster_result structure */
    Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));
    raster_result = raster_id_to_result(raster_id);

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(result->id, (seq_reg_data *)&info, 0);

    /* deregister horizontal sequence */
    seq_deregister(GetSeqNum(result->seq_id[HORIZONTAL]), sim_callback, 
		   (seq_result *)result);

    /* deregister vertical sequence */
    seq_deregister(GetSeqNum(result->seq_id[VERTICAL]), sim_callback, 
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
	    verror(ERR_WARN, "sim_shutdown", "%s \n", Tcl_GetStringResult(interp));
	
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
	    verror(ERR_WARN, "sim_shutdown", "%s\n", Tcl_GetStringResult(interp));
	}
    }
    DestroySequencePairDisplay(interp, id);
    free(input->params);
    SipFreeResult(result);

    if (raster_result) 
	DeleteResultFromRaster(raster_result);
}

void sim_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    in_sim *input = result->input;
    out_raster *output = result->output;
    d_plot *data = result->data;
    int id = result->id;
    char cmd[1024];
    int seq1_num, seq2_num;

    seq1_num = GetSeqNum(result->seq_id[HORIZONTAL]);
    seq2_num = GetSeqNum(result->seq_id[VERTICAL]);

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Local alignment");
	break;
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "local #%d", result->id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "local alignment: hori=%s vert=%s", 
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
	      sim_shutdown(interp, seq_num, result, output->raster_win, id);
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
	sim_shutdown(interp, seq_num, result, output->raster_win, id);
    }
    }

}

void sim_text_func(void *obj)
{

}

void sim_align(char *seq1,
	       char *seq2,
	       int seq1_len,
	       int seq2_len,
	       int seq_type,
	       int *num_alignments, /* in, out */
	       float score_align,
	       float match,
	       float transition,
	       float transversion,
	       float start_gap,
	       float cont_gap,
	       align_int **res,     /* out */
	       long *start1,        /* out */
	       long *start2,        /* out */
	       long *end1,          /* out */
	       long *end2)          /* out */
{ 
    long  M, N, K;		/* Sequence lengths and k     */
    char  *A,  *B;		/* The two sequences      */
    int i, j;
    long V[128][128], Q,R;		/* Converted integer weights  */
    float parm_M, parm_I, parm_V, parm_O, parm_E;
    int same_seq;
    char achars[] = "ARNDCQEGHILKMFPSTWYVBZX";	/* amino acid names */
    int	alpha = 23;			/* alphabet size */

    A = seq1;
    --A;	/* subscripts start with 1 */
    B = seq2;
    --B;	/* subscripts start with 1 */
    M = seq1_len;
    N = seq2_len;
    K = *num_alignments;

    parm_M = DEFAULT_M;
    parm_I = DEFAULT_I;
    parm_V = DEFAULT_V;
    parm_O = DEFAULT_O;
    parm_E = DEFAULT_E;

    parm_M = match;
    parm_I = transition;
    parm_V = transversion;
    parm_O = start_gap;
    parm_E = cont_gap;

#ifdef DEBUG
    printf("#:lav\n\n");
    printf("d {\n  \"SIM output with parameters:\n");
    if (seq_type == DNA) 
	printf("  M = %g, I = %g, V = %g\n", parm_M, parm_I, parm_V);
    printf("  O = %g, E = %g\"\n}\n", parm_O, parm_E);
#endif

    if (seq_type == PROTEIN) {
	set_char_set(PROTEIN);
	set_score_matrix(get_matrix_file(PROTEIN));

        for (i = 0; i < 128; i++) {
	    for (j = 0; j < 128; j++) {
	        V[i][j] = score_matrix[char_lookup['-']][char_lookup['-']];
	    }
	}

	for (i = 0; i < alpha; ++i) {
	    for (j = 0; j < alpha; ++j) {
		V[achars[i]][achars[j]] = 10*score_matrix[char_lookup[achars[i]]][char_lookup[achars[j]]];
	    }
	}
    } else {
        int j;
	parm_V += (parm_V > 0) ? 0.05 : -0.05;
	/* initialise V - HACK - don't know what to! */
        for (i = 0; i < 128; i++) {
	  for (j = 0; j < 128; j++) {
	        V[i][j] = parm_V;
	  }
	}

	parm_M += (parm_M > 0) ? 0.05 : -0.05;
	V['A']['A'] = V['C']['C'] = V['G']['G'] = V['T']['T']
	    = 10*parm_M;
	V['a']['a'] = V['c']['c'] = V['g']['g'] = V['t']['t']
	    = 10*parm_M;
	V['a']['A'] = V['c']['C'] = V['g']['G'] = V['t']['T']
	    = 10*parm_M;
	V['A']['a'] = V['C']['c'] = V['G']['g'] = V['T']['t']
	    = 10*parm_M;

	parm_I += (parm_I > 0) ? 0.05 : -0.05;
	V['A']['G'] = V['G']['A'] = V['C']['T'] = V['T']['C']
	    = 10*parm_I;
	V['a']['g'] = V['g']['a'] = V['c']['t'] = V['t']['c']
	    = 10*parm_I;
	V['a']['G'] = V['g']['A'] = V['c']['T'] = V['t']['C']
	    = 10*parm_I;
	V['A']['g'] = V['G']['a'] = V['C']['t'] = V['T']['c']
	    = 10*parm_I;

	V['A']['C'] = V['A']['T'] = V['C']['A'] = V['C']['G'] =
	    V['G']['C'] = V['G']['T'] = V['T']['A'] = V['T']['G']
	    = 10*parm_V;
	V['a']['c'] = V['a']['t'] = V['c']['a'] = V['c']['g'] =
	    V['g']['c'] = V['g']['t'] = V['t']['a'] = V['t']['g']
	    = 10*parm_V;
	V['a']['C'] = V['a']['T'] = V['c']['A'] = V['c']['G'] =
	    V['g']['C'] = V['g']['T'] = V['t']['A'] = V['t']['G']
	    = 10*parm_V;
	V['A']['c'] = V['A']['t'] = V['C']['a'] = V['C']['g'] =
	    V['G']['c'] = V['G']['t'] = V['T']['a'] = V['T']['g']
	    = 10*parm_V;
    }

    parm_O += (parm_O > 0) ? 0.05 : -0.05;
    Q = 10 * parm_O;
    parm_E += (parm_E > 0) ? 0.05 : -0.05;
    R = 10 * parm_E;
    
    if (seq1_len == seq2_len) {
	if (strcmp(seq1, seq2) == 0) {
	    same_seq = 1;
	} else {
	    same_seq = 0;
	}
    } else {
	same_seq = 0;
    }

    if (!same_seq)
        *num_alignments = SIM(A,B,M,N,K,V,Q,R,2L, score_align, res, start1, start2, end1, end2);
    else {
	/* self-comparison; insert trivial diagonal */
        start1[0] = 1;
	start2[0] = 1;
	end1[0] = M;
	end2[0] = M;
	*res[0] = 0;
	/* 
	 * if there is only 1 alignment to be found, it must be the main
	 * diagonal
	 */
	if (*num_alignments == 1)
	    return;
	*num_alignments = SIM(A,A,M,M,K-1,V,Q,R,1L, score_align, &res[1], 
			      &start1[1], &start2[1], &end1[1], &end2[1]);
	/* add one to num_alignments for the main diagonal */
	(*num_alignments)++;
    }
}

void find_seq_lengths(align_int *S,
		      long M,
		      long N,
		      int *len1,
		      int *len2)
{
    int i, j;
    long op;

    i = j = 0;
    while (i < M && j < N) {
	++i;
	++j;
	op = *S;
	if (op > 0) {
	    (*len2) += op;
	} else if (op < 0) {
	    (*len1) -= op;
	} else {
	    (*len1)++;
	    (*len2)++;
	}
	S++;
    }
}

void store_sim1(char *A,
		char *B,
		long seq1_len,
		long seq2_len,
		long M,
		long N,
		align_int *S,
		long AP,
		long BP,
		pt_score *p_array,
		int *k)
{
    long   i, j, op, start_i, start_j, match, mis_match;

    for (i = j = 0; i < M || j < N; ) {
	start_i = i;
	start_j = j;
	match = mis_match = 0;
	while (i < M && j < N && *S == 0) {
	    ++i;
	    ++j;
	    if (A[i] == B[j])
		++match;
	    else
		++mis_match;
	    S++;
	}

	p_array[*k].x = (int)AP + start_i;
	p_array[(*k)++].y = (double)BP + start_j;
	
	p_array[*k].x = (int)AP + i - 1;
	p_array[(*k)++].y = (double)BP + j - 1;

#ifdef DEBUG
	printf("  l %d %d %d %d %1.1f\n", AP+start_i, BP+start_j, AP+i-1,
	       BP+j-1, (float)(100*match)/(float)(match+mis_match));
#endif
	if (i < M || j < N)
	    if ((op = *S++) > 0)
		j += op;
	    else
		i -= op;
    }

    /* delimiter between sets of alignments */
    p_array[*k].x = -1;
    p_array[*k].y = -1;
    p_array[(*k)++].score = -1;
}

int store_sim2(int seq1_num,
	       int seq2_num,
	       int start_h,
	       int end_h,
	       int start_v,
	       int end_v,
	       in_sim *input,
	       d_plot *data,
	       int n_pts)
{
    seq_result *result;
    int id;

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    /* realloc p_array to be the actual number of lines */
    if (NULL == (data->p_array = (pt_score *)realloc(data->p_array, (n_pts+1)* 
					       sizeof(pt_score))))
	return -1;

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
    result->op_func = sim_callback;
    result->txt_func = sim_text_func;

    seq_register(seq1_num, sim_callback, (void *)result, SEQ_PLOT_PERM, id);
    seq_register(seq2_num, sim_callback, (void *)result, SEQ_PLOT_PERM, id);

    return id;
}

int init_sip_local_align_create(Tcl_Interp *interp, 
				int seq_id_h,
				int seq_id_v, 
				int start_h,
				int end_h, 
				int start_v,
				int end_v, 
				int num_align,
				float score_align,
				float match,
				float transition,
				float transversion,
				float start_gap,
				float cont_gap,
				int *id)
{
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int r_len1, r_len2;
    char *r_seq1, *r_seq2;
    char *name1, *name2;
    int i;
    int seq1_type, seq2_type;
    int seq1_num;
    int seq2_num;
    in_sim *input;
    align_int **res;
    long *start1, *start2, *end1, *end2;
    int cnt = 0;
    int max_align;
    Tcl_DString input_params;
    int num_elements;
    d_plot *data;

#define NUM_SCORES 100
    
    vfuncheader("local alignment");

    if (-1 == (seq1_num = GetSeqNum(seq_id_h))) {
	verror(ERR_WARN, "local alignment", 
	       "horizontal sequence undefined");
	goto error;
    }
     
    if (-1 == (seq2_num = GetSeqNum(seq_id_v))) {
	verror(ERR_WARN, "local alignment", 
	       "vertical sequence undefined");
	goto error;
    }

    /* only align dna or protein */
    seq1_type = GetSeqType(seq1_num);
    seq2_type = GetSeqType(seq2_num);

    if (seq1_type != seq2_type) {
	verror(ERR_FATAL, "sim alignment", "sequences must both be either DNA or protein");
	goto error;
    }

    seq1 = GetSeqSequence(seq1_num);
    if ((seq1_len = end_h - start_h + 1) < 1) {
	verror(ERR_WARN, "align sequences", "negative length");
	goto error;
    }
    seq2 = GetSeqSequence(seq2_num);
    if ((seq2_len = end_v - start_v + 1) < 1) {
	verror(ERR_WARN, "align sequences", "negative length");
	goto error;
    }

    if (NULL == (input = (in_sim *)xmalloc(sizeof(in_sim))))
	goto error;

     Tcl_DStringInit(&input_params);
     vTcl_DStringAppend(&input_params, "horizontal %s: %s from %d to %d\n"
	    "vertical %s: %s from %d to %d\n", 
	    GetSeqLibraryName(seq1_num), 
	    GetSeqName(seq1_num), 
	    start_h, end_h, 
	    GetSeqLibraryName(seq2_num), 
	    GetSeqName(seq2_num), 
	    start_v, end_v);

     if (score_align == -1) {
	 vTcl_DStringAppend(&input_params, "number of alignments %d \n",
			    num_align);
     } else {
	 vTcl_DStringAppend(&input_params, "alignments above score %g\n",
			    score_align);
     }
     
     if (GetSeqType(seq1_num) == DNA) {
	 vTcl_DStringAppend(&input_params, "score for match %g\n"
		  "score for transition %g\n"
		  "score for transversion %g\n",
		  match, transition, transversion);
     }
     vTcl_DStringAppend(&input_params, "penalty for starting gap %g\n"
			"penalty for each residue in gap %g\n",
			start_gap, cont_gap);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);
     
    if (NULL == (r_seq1 = (char *)xcalloc(((seq1_len+seq2_len)*2+1),
					  sizeof(char)))) {
	goto error;
    }
    if (NULL == (r_seq2 = (char *)xcalloc(((seq1_len+seq2_len)*2+1),
					  sizeof(char)))) {
	goto error;
    }
    
    if (score_align != -1) {
	num_align = NUM_SCORES;
    }
    max_align = num_align;

    if (NULL == (start1 = (long *)xmalloc(max_align * sizeof(long)))) {
      goto error;
    }
    if (NULL == (start2 = (long *)xmalloc(max_align * sizeof(long)))) {
	goto error;
    }
    if (NULL == (end1 = (long *)xmalloc(max_align * sizeof(long)))) {
	goto error;
    }
    if (NULL == (end2 = (long *)xmalloc(max_align * sizeof(long)))) {
	goto error;
    }
    if (NULL == (res = (align_int **)xmalloc(max_align *sizeof(align_int*)))) {
	goto error;
    }

    for (i = 0; i < max_align; i++) {
	if (NULL == (res[i] = (align_int *)xcalloc((seq1_len+seq2_len+1),
						   sizeof(align_int)))) {
	    goto error;
	}
    }
    
    /* 
     * if finding all alignments above a certain score, the return value of
     * num_align is the number of alignments found
     */
    sim_align(&seq1[start_h-1], &seq2[start_v-1], seq1_len,
	      seq2_len, seq1_type, &num_align, score_align,
	      match, transition, transversion, start_gap, cont_gap, res, 
	      start1, start2, end1, end2);

    if (num_align <= 0) {
	verror(ERR_WARN, "local alignment", "no matches found\n");
	goto error;
    }
    name1 = GetSeqBaseName(seq1_num);
    name2 = GetSeqBaseName(seq2_num);

    num_elements = (seq1_len + seq2_len + 1) * num_align;
    if (NULL == (data = (d_plot *)xmalloc(sizeof(d_plot))))
	goto error;
    
    if (NULL == (data->p_array = (pt_score *)xmalloc(sizeof(pt_score) * 
						     num_elements)))
	goto error;
    
    for (i = 0; i < num_align; i++) {
	store_sim1(&seq1[start_h+start1[i]-2], 
		  &seq2[start_v+start2[i]-2], seq1_len, 
		  seq2_len, end1[i]-start1[i]+1, end2[i]-start2[i]+1,
		  res[i], start_h+start1[i]-1, start_v+start2[i]-1,
		  data->p_array, &cnt);

	cexpand(&seq1[start_h+start1[i]-2], 
		&seq2[start_v+start2[i]-2], end1[i]-start1[i]+1,
		end2[i]-start2[i]+1, r_seq1, r_seq2, &r_len1, &r_len2, 
		ALIGN_J_SSH | ALIGN_J_PADS, res[i]);
	
	spin_list_alignment(r_seq1, r_seq2, name1, name2, start_h+start1[i]-1, 
		       start_v+start2[i]-1, "", seq1_type);
	
    }
    *id = store_sim2(seq1_num, seq2_num, start_h, end_h, start_v, end_v, 
		     input, data, cnt);

    xfree(r_seq1);
    xfree(r_seq2);
    xfree(start1);
    xfree(start2);
    xfree(end1);
    xfree(end2);
    for (i = 0; i < max_align; i++) {
	xfree(res[i]);
    }
    xfree(res);
    
    return 0;

 error:
    return -1;
}

int init_sip_local_align_plot(Tcl_Interp *interp, 
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

    init_dot_plot(interp, seq_id_h, seq_id_v, result_id, "local", raster_win,
		  raster_id, opts, 6, LINE, data->dim);
    
    xfree(opts[1]);
    xfree(opts[3]);
    xfree(opts[5]);

    return 0;
}
