#include   <stdio.h>
#include   <ctype.h>
#include <math.h>
#include <string.h>


#include "sim.h"
#include "seq_results.h"
#include "seq_reg.h"
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
#include "tclCanvGraph.h"
#include "nip_results.h"
#include "seq_element.h"


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
		  seq_result *s_result,
		  element *e)
{
    seq_reg_key_name info;
    static char buf[80];
    in_sim *input = s_result->input;

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    /* deregister horizontal sequence */
    seq_deregister(GetSeqNum(s_result->seq_id[HORIZONTAL]), sim_callback, 
		   (seq_result *)s_result);

    /* deregister vertical sequence */
    seq_deregister(GetSeqNum(s_result->seq_id[VERTICAL]), sim_callback, 
		   (seq_result *)s_result);
    
    if (e->num_results > 0) {
	Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);     
	e->replot_func(e);
    } else {
	Tcl_VarEval(e->c->interp, "result_list_update ", NULL);
    }
    DestroySequencePairDisplay(interp, s_result->id);
    free(input->params);
}

void sim_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_sim *input = s_result->input;
    int result_id = s_result->id;
    char cmd[1024];
    int seq1_num, seq2_num;
    element *e = s_result->e;
    plot_data *result;
    Tcl_Interp *interp;

    if (e) {
	result = find_plot_data(e, result_id);
	interp = e->c->interp;
    }

    seq1_num = GetSeqNum(s_result->seq_id[HORIZONTAL]);
    seq2_num = GetSeqNum(s_result->seq_id[VERTICAL]);

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Local alignment");
	break;
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "local #%d", result_id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "local alignment: hori=%s vert=%s", 
		GetSeqBaseName(GetSeqNum(s_result->seq_id[HORIZONTAL])),
		GetSeqBaseName(GetSeqNum(s_result->seq_id[VERTICAL])));
	break;
    case SEQ_GET_OPS:
	if (result->hidden) {
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
	    Tcl_Eval(interp, "SetBusy");
	    vfuncheader("results");
	    /* result->txt_func(result); */
	    Tcl_Eval(interp, "ClearBusy");
	    break;
	case 2: /* configure */
	    sprintf(cmd, "result_config %d %d %s %d %s", 
		    result_id, result->line_width, result->colour, e->id, 
		    e->win);
	    
	    if (TCL_OK != Tcl_Eval(interp, cmd)){
		puts(interp->result);
	    }	    
	    break;
	case 3: /* display sequences */
	    SequencePairDisplay(interp, e->win, result_id, 
				s_result->seq_id[HORIZONTAL], 
				s_result->seq_id[VERTICAL]);
	    break;
	case 4: /* hide all */
	    result->hidden = 1;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 5: /* reveal all */
	    result->hidden = 0;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 6: /* remove */ 
	    {
		sim_shutdown(interp, s_result, e);
		remove_result_from_element(e, result_id);
	      break;
	    }
	}
	break;
    case SEQ_PLOT:
	s_result->pr_func(s_result, NULL);
	break;
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case RESULT:
	    jdata->info.result = (void *)s_result;
	    break;
	case WIN_NAME:
	    {
		char *r_win =  e->win;
		jdata->info.result = (void *)r_win;
		break;
	    }
	case WIN_SIZE:
	    {
		static d_point pt;
		pt.x = get_default_int(interp, tk_utils_defs, 
					w("DOT.PLOT_WIDTH"));
		pt.y = get_default_double(interp, tk_utils_defs,
					   w("DOT.PLOT_HEIGHT"));
		jdata->info.result = (void *)&pt;
		
		break;
	    } /* WIN_SIZE */
	}
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: {
	sim_shutdown(interp, s_result, e);
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
		int strand,
		parray *p,
		int *k)
{
    long   i, j, op, start_i, start_j, match, mis_match;
    int num = 0;

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

	if (strand == TOP_S) {
	    p->p_array[num].x = (int)AP + start_i;
	    p->p_array[num].y = (double)BP + start_j;
	    (*k)++;
	    num++;
	    p->p_array[num].x = (int)AP + i - 1;
	    p->p_array[num].y = (double)BP + j - 1;
	} else {
	    p->p_array[num].x = (int)AP + start_i;
	    p->p_array[num].y = (double)BP + j - 1;
	    (*k)++;
	    num++;
	    p->p_array[num].x = (int)AP + i - 1;
	    p->p_array[num].y = (double)BP + start_j;

	}

	(*k)++;
	num++;
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
	
    /* array must be in ascending x order */
    qsort((void *) p->p_array, num, sizeof(g_pt),  compare_g_pt);

    p->n_pts = num;
    p->type = G_LINE;
}

int store_sim2(int seq1_num,
	       int seq2_num,
	       int start_h,
	       int end_h,
	       int start_v,
	       int end_v,
	       in_sim *input,
	       Graph *graph,
	       int n_pts,
	       Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    int id;

    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    id = get_reg_id();

    graph->dim.x0 = start_h;
    graph->dim.x1 = end_h;
    graph->dim.y0 = start_v;
    graph->dim.y1 = end_v;

    *graph_obj = Tcl_NewGraphObj(graph);
    /* these are copied over in Tcl_NewGraphObj so don't need them anymore */
    xfree(graph->p_arrays[0].p_array);
    xfree(graph->p_arrays);
    xfree(graph);

    s_result->seq_id[HORIZONTAL] = GetSeqId( seq1_num);
    s_result->seq_id[VERTICAL] = GetSeqId(seq2_num);

    s_result->input = (void *)input; 
    s_result->id = id; 
    s_result->type = SEQ_TYPE_DOT_PLOT;
    s_result->gr_type = SEQ_TYPE_DOT_PLOT;
    s_result->frame = 0;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */

    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = sim_callback;
    s_result->txt_func = sim_text_func;

    seq_register(seq1_num, sim_callback, (void *)s_result, SEQ_PLOT_PERM, id);
    seq_register(seq2_num, sim_callback, (void *)s_result, SEQ_PLOT_PERM, id);

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
				int strand,
				Tcl_Obj **graph_obj,
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
    Graph *graph;
    char strand_s[7];

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
     if (strand & TOP_S) {
	 strcpy(strand_s, "top");
     } else if (strand & BOTTOM_S) {
	 strcpy(strand_s, "bottom");
     } else {
	 strcpy(strand_s, "both");
     }
     vTcl_DStringAppend(&input_params, "horizontal %s: %s from %d to %d\n"
	    "vertical %s: %s from %d to %d\n strand %s", 
	    GetSeqLibraryName(seq1_num), 
	    GetSeqName(seq1_num), 
	    start_h, end_h, 
	    GetSeqLibraryName(seq2_num), 
	    GetSeqName(seq2_num), 
	    start_v, end_v, strand_s);

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
    if (strand == TOP_S) {
	
	sim_align(&seq1[start_h-1], &seq2[start_v-1], seq1_len,
		  seq2_len, seq1_type, &num_align, score_align,
		  match, transition, transversion, start_gap, cont_gap, res, 
		  start1, start2, end1, end2);
	if (num_align <= 0) {
	    verror(ERR_WARN, "local alignment", "no matches found\n");
	    goto error;
	}
    } else {
	if (seq2_type == DNA)
	    complement_seq(seq2, seq2_len);

	sim_align(&seq1[start_h-1], &seq2[start_v-1], seq1_len,
		  seq2_len, seq1_type, &num_align, score_align,
		  match, transition, transversion, start_gap, cont_gap, res, 
		  start1, start2, end1, end2);
	if (num_align <= 0) {
	    verror(ERR_WARN, "local alignment", "no matches found\n");
	    goto error;
	}
	
	if (seq2_type == DNA)
	    complement_seq(seq2, seq2_len);
    }

    name1 = GetSeqBaseName(seq1_num);
    name2 = GetSeqBaseName(seq2_num);

    num_elements = (seq1_len + seq2_len + 1) * num_align;
    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;
    
    Tcl_InitGraph(&graph);

    if (NULL == (graph->p_arrays = (parray *)xmalloc(num_align * sizeof(parray))))
	 return -1;

    graph->n_parrays = num_align;
    for (i = 0; i < num_align; i++) {
	if (NULL == (graph->p_arrays[i].p_array = (g_pt *)xmalloc(sizeof(g_pt) * 
								  (seq1_len + seq2_len + 1))))
	    return -1;

	store_sim1(&seq1[start_h+start1[i]-2], 
		  &seq2[start_v+start2[i]-2], seq1_len, 
		  seq2_len, end1[i]-start1[i]+1, end2[i]-start2[i]+1,
		  res[i], start_h+start1[i]-1, start_v+start2[i]-1,
		  strand, &graph->p_arrays[i], &cnt);

	cexpand(&seq1[start_h+start1[i]-2], 
		&seq2[start_v+start2[i]-2], end1[i]-start1[i]+1,
		end2[i]-start2[i]+1, r_seq1, r_seq2, &r_len1, &r_len2, 
		ALIGN_J_SSH | ALIGN_J_PADS, res[i]);
	
	spin_list_alignment(r_seq1, r_seq2, name1, name2, start_h+start1[i]-1, 
		       start_v+start2[i]-1, "", seq1_type);
	
    }
    *id = store_sim2(seq1_num, seq2_num, start_h, end_h, start_v, end_v, 
		     input, graph, cnt, graph_obj);

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
			      char *e_win, 
			      int element_id,
			      char *c_win, 
			      int container_id,
			      Tcl_Obj *results,
			      int line_width,
			      char *colour,
			      char *element_type)
{
    seq_result *s_result;
    Graph *graph;
    plot_data *result;
    configs *configure;

    s_result = seq_id_to_result(result_id);

    if (NULL == (result = (plot_data *)xmalloc(sizeof(plot_data))))
	return -1;

    if (NULL == (result->configure = (configs **)xmalloc(sizeof(configs*))))
	return -1;
    
    if (NULL == (configure = (configs *)xmalloc(sizeof(configs))))
	return -1;

    configure->position = -1.0;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = 1.0;
    configure->zoom = 2;
    configure->scroll = 1;

    result->configure[0] = configure;
    result->n_configure = 1;
    result->sf_m = 1.0;
    result->sf_c = 0.0;
    result->result_id = result_id;
    result->scale = SCALE_X | SCALE_Y;
    result->hidden = 0;
    result->line_width = line_width;
    result->colour = strdup(colour);
    result->len_ruler = 1;
    result->amp_ruler = 0;
    sprintf(result->tags, "id%d", result_id);

    graph = Tcl_GetGraphFromObj(results);

    init_seq_element(interp, s_result, seq_id_h, seq_id_v, result, 
		     container_id, element_id, e_win, c_win, 
		     HORIZONTAL|VERTICAL, HORIZONTAL|VERTICAL, graph, 
		     element_type);
    return 0;
   
}
