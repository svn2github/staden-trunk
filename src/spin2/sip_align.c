#include <string.h>

#include "seq_results.h"
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
#include "tclCanvGraph.h"
#include "nip_results.h"
#include "seq_element.h"


static char *aa_order = "ABCDEFGHIKLMNPQRSTVWYZX*-";
static char *nt_order = "ACGT";

void align_callback(int seq_num, void *obj, seq_reg_data *jdata);

void align_shutdown(Tcl_Interp *interp,
		    seq_result *s_result,
		    element *e)
{
    seq_reg_key_name info;
    static char buf[80];
    in_align *input = s_result->input;

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    /* deregister horizontal sequence */
    seq_deregister(GetSeqNum(s_result->seq_id[HORIZONTAL]), align_callback, 
		   (seq_result *)s_result);

    /* deregister vertical sequence */
    seq_deregister(GetSeqNum(s_result->seq_id[VERTICAL]), align_callback, 
		   (seq_result *)s_result);
    
    if (e->num_results > 1) {
	Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	e->replot_func(e);
    } else {
	Tcl_VarEval(e->c->interp, "result_list_update ", NULL);
    }
    DestroySequencePairDisplay(interp, s_result->id);
    free(input->params);
}

void align_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_align *input = s_result->input;
    int result_id = s_result->id;
    char cmd[1024];
    element *e = s_result->e;
    plot_data *result;
    Tcl_Interp *interp;

    if (e) {
	result = find_plot_data(e, result_id);
	interp = e->c->interp;
    }

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Align sequences");
	break;
    case SEQ_KEY_NAME:
	sprintf(jdata->name.line, "align #%d", result_id);
	break;

    case SEQ_GET_BRIEF:
	sprintf(jdata->name.line, "align: hori=%s vert=%s", 
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
	      align_shutdown(interp, s_result, e);
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
	align_shutdown(interp, s_result, e);
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
    int r_len1, r_len2;
    int seq1_type;
    int start_h, start_v;
    char *name1, *name2;

    seq_num_h = GetSeqNum(result->seq_id[HORIZONTAL]);
    seq_num_v = GetSeqNum(result->seq_id[VERTICAL]);

    seq1 = GetSeqSequence(seq_num_h);
    seq1_len = GetSeqLength(seq_num_h);
    seq2 = GetSeqSequence(seq_num_v);
    seq2_len = GetSeqLength(seq_num_v);

    seq1_type = GetSeqType(seq_num_h);

    if (NULL == (r_seq1 = (char *)xcalloc(((seq1_len+seq2_len)*2+1),
					  sizeof(char)))) {
	return;
    }
    if (NULL == (r_seq2 = (char *)xcalloc(((seq1_len+seq2_len)*2+1),
					  sizeof(char)))) {
	return;
    }

    start_h = ;
    start_v = ;


    cexpand(&seq1[start_h-1], &seq2[start_v-1], seq1_len, seq2_len, 
	    r_seq1, r_seq2, &r_len1, &r_len2, 
	    ALIGN_J_SSH | ALIGN_J_PADS, res);

    spin_list_alignment(r_seq1, r_seq2, name1, name2, start_h, start_v, "", 
			seq1_type);
    
    xfree(r_seq1);
    xfree(r_seq2);
#endif

}

/*
 * save default name of alignment based on the parent name with _uid appended
 */
int
sip_save_alignment(Tcl_Interp *interp,
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

int store_align(int seq1_num,
		int seq2_num,
		int start_h,
		int end_h,
		int start_v,
		int end_v,
		int seq1_len,
		int seq2_len,
		in_align *input,
		align_int *res,
		int strand,
		Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    int i, j, k, d, id;
    int cnt = 0;
    int num_elements = seq1_len + seq2_len + 1;
    int n_pts = 0;
    Graph *graph;
    /*
    for (i = 0; i < (seq1_len + seq2_len + 1); i++) {
	printf("res[%d]=%d\n", i, res[i]);
    }
    */
    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;
    
    Tcl_InitGraph(&graph);

    if (NULL == (graph->p_arrays = (parray *)xmalloc(sizeof(parray))))
	 return -1;

    if (NULL == (graph->p_arrays[0].p_array = (g_pt *)xmalloc(num_elements * sizeof(g_pt))))
	return -1;
    

    /* data assignment */
    i = 0;
    j = 0;
    k = 0;
    d = 0;

    if (strand == TOP_S) {
	i = start_h;
	j = start_v;

	graph->p_arrays[0].p_array[k].x = i;
	graph->p_arrays[0].p_array[k].y = j;
	while (i <= end_h || j <= end_v) {
	    if (res[cnt] == 0) {
		i++;
		j++;
		if (cnt == 0)
		    k++;
		
		if (cnt && res[cnt-1])
		    k++;
		graph->p_arrays[0].p_array[k].x = i;
		graph->p_arrays[0].p_array[k].y = j;
	    } else if (res[cnt] < 0) {
		i -= res[cnt];
		graph->p_arrays[0].p_array[++k].x = i;
		graph->p_arrays[0].p_array[k].y = j;
	    } else if (res[cnt] > 0) {
		j += res[cnt];
		graph->p_arrays[0].p_array[++k].x = i;
		graph->p_arrays[0].p_array[k].y = j;
	    }
	    cnt++;
	}
    } else {
	i = start_h;
	j = end_v;
	graph->p_arrays[0].p_array[0].y = j;
	while (i <= end_h || j >= start_v) {
	    if (res[cnt] == 0) {
		i++;
		j--;
		if (cnt == 0)
		    k++;
		
		if (cnt && res[cnt-1])
		    k++;
		graph->p_arrays[0].p_array[k].x = i;
		graph->p_arrays[0].p_array[k].y = j;
	    } else if (res[cnt] < 0) {
		i -= res[cnt];
		graph->p_arrays[0].p_array[++k].x = i;
		graph->p_arrays[0].p_array[k].y = j;

	    } else if (res[cnt] > 0) {
		j -= res[cnt];
		graph->p_arrays[0].p_array[++k].x = i;
		graph->p_arrays[0].p_array[k].y = j;
	    }
	    cnt++;
	}
    }
    k++;
    n_pts = k;

     /* array must be in ascending x order */
    qsort((void *) graph->p_arrays[0].p_array, n_pts, sizeof(g_pt), 
	  compare_g_pt);

    graph->n_parrays = 1;
    graph->p_arrays[0].n_pts = n_pts;
    graph->p_arrays[0].type = G_LINE;
    graph->dim.x0 = start_h;
    graph->dim.x1 = end_h;
    graph->dim.y0 = start_v;
    graph->dim.y1 = end_v;

    *graph_obj = Tcl_NewGraphObj(graph);
    /* these are copied over in Tcl_NewGraphObj so don't need them anymore */
    xfree(graph->p_arrays[0].p_array);
    xfree(graph->p_arrays);
    xfree(graph);

    id = get_reg_id();
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
    s_result->op_func = align_callback;
    s_result->txt_func = align_text_func;

    seq_register(seq1_num, align_callback, (void *)s_result, SEQ_PLOT_PERM, id);
    seq_register(seq2_num, align_callback, (void *)s_result, SEQ_PLOT_PERM, id);

    return id;
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
				 int strand, 
				 Tcl_Obj **graph_obj,
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
    char strand_s[7];

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
   
    if (strand == TOP_S) {
	calign(&seq1[start_h-1], &seq2[start_v-1], seq1_len, seq2_len, 
	       NULL, NULL, NULL, NULL,
	       0, 0,                            /* low, high band */
	       start_gap, cont_gap,             /* gap_open, gap_extend */
	       3,                               /* job */
	       seqtype,                         /* is protein */
	       res);
	cexpand(&seq1[start_h-1], &seq2[start_v-1], seq1_len, seq2_len, 
		r_seq1, r_seq2, &r_len1, &r_len2, 
		ALIGN_J_SSH | ALIGN_J_PADS, res);
    } else {
	if (seq2_type == DNA)
	    complement_seq(seq2, seq2_len);
	calign(&seq1[start_h-1], &seq2[start_v-1], seq1_len, seq2_len, 
	       NULL, NULL, NULL, NULL,
	       0, 0,                            /* low, high band */
	       start_gap, cont_gap,             /* gap_open, gap_extend */
	       3,                               /* job */
	       seqtype,                         /* is protein */
	       res);
	cexpand(&seq1[start_h-1], &seq2[start_v-1], seq1_len, seq2_len, 
		r_seq1, r_seq2, &r_len1, &r_len2, 
		ALIGN_J_SSH | ALIGN_J_PADS, res);
	
	if (seq2_type == DNA)
	    complement_seq(seq2, seq2_len);
    }
    Tcl_DStringInit(&input_params);
    if (strand & TOP_S) {
        strcpy(strand_s, "top");
    } else if (strand & BOTTOM_S) {
        strcpy(strand_s, "bottom");
    } else {
        strcpy(strand_s, "both");
    }
    
    vTcl_DStringAppend(&input_params, "horizontal %s: %s from %d to %d\n"
		       "vertical %s: %s from %d to %d\nstrand %s\nscore for match %d\n"
		       "score for mis-match %d\npenalty for starting gap %d\n"
		       "penalty for each residue in gap %d\n", 
		       GetSeqLibraryName(seq1_num), GetSeqName(seq1_num), 
		       start_h, end_h, GetSeqLibraryName(seq2_num), 
		       GetSeqName(seq2_num), start_v, end_v, strand_s, match, 
		       mis_match, 
		       start_gap, cont_gap);

    vfuncparams("%s", Tcl_DStringValue(&input_params));
    input->params = strdup(Tcl_DStringValue(&input_params)); 
    Tcl_DStringFree(&input_params);

    if (-1 == (*id = store_align(seq1_num, seq2_num, start_h,
				 end_h, start_v, end_v, seq1_len, seq2_len,
				 input, res, strand, graph_obj))) {
	goto error;
    }

    sip_save_alignment(interp, seq1_num, r_seq1, name1);
    sip_save_alignment(interp, seq2_num, r_seq2, name2);

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

