#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

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
#include "tclCanvGraph.h"
#include "seq_element.h"
#include "element_canvas.h"

#define MAXMATCHES 10000

void nip_stop_codons_callback(int seq_num, void *obj, seq_reg_data *jdata);

void nip_stop_codons_shutdown(Tcl_Interp *interp,
			      seq_result *s_result,
			      element *e,
			      int seq_num)
{
    in_s_codon *input = s_result->input;
    seq_reg_key_name info;
    static char buf[80];

#ifdef DEBUG
    printf("nip_stop_codons_shutdown\n");
#endif

    /* find key name BEFORE deregister */
    info.job = SEQ_KEY_NAME;
    info.line = buf;
    seq_result_notify(s_result->id, (seq_reg_data *)&info, 0);

    seq_deregister(seq_num, nip_stop_codons_callback, (seq_result *)s_result);

    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);

    if (e->num_results > 0) {
	e->replot_func(e);
    }

    free(input->params);
    xfree(s_result->input);
    xfree(s_result);
}


void nip_stop_codons_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *s_result = (seq_result *) obj;
    in_s_codon *input = s_result->input;
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
	if (s_result->type == SEQ_TYPE_STOPCODON) {
	    sprintf(jdata->name.line, "Plot stop codons");
	} else {
	    sprintf(jdata->name.line, "Plot start codons");
	}
	break;    

    case SEQ_KEY_NAME:
	if (s_result->type == SEQ_TYPE_STOPCODON) {
	    sprintf(jdata->name.line, "stop f%d #%d", s_result->frame,
		    result_id);
	} else {
	    sprintf(jdata->name.line, "start f%d #%d", s_result->frame,
		    result_id);
	}
	break;
	
    case SEQ_GET_BRIEF:
	if (s_result->type == SEQ_TYPE_STOPCODON) {
	    sprintf(jdata->name.line, "stop codons: seq=%s frame=%d", 
		    GetSeqName(GetSeqNum(s_result->seq_id[HORIZONTAL])), 
		    s_result->frame);
	} else {
	    sprintf(jdata->name.line, "start codons: seq=%s frame=%d", 
		    GetSeqName(GetSeqNum(s_result->seq_id[HORIZONTAL])), 
		    s_result->frame);
	}
	break;

    case SEQ_GET_OPS:
	if (result->hidden) {
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
	    Tcl_Eval(interp, "SetBusy");
	    vfuncheader("results");
	    s_result->txt_func(s_result);
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
	case 3: /* hide all */
	    result->hidden = 1;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 4: /* reveal all */
	    result->hidden = 0;
	    Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);
	    e->replot_func(e);
	    break;
	case 5: /* remove */ 
	    {
		nip_stop_codons_shutdown(interp, s_result, e, seq_num);
		remove_result_from_element(e, result_id);
		break;
	    }
	}
	break;
    case SEQ_PLOT:
	{
	    Tcl_Obj *graph_obj = (Tcl_Obj *) s_result->data;
	    Graph *graph;
	    
	    /* HACK - but the only way I have discovered to making stop 
	       codons work. Need to set the dim.y1 to be the height of the
	       current canvas (for instance moving codons from gene pref plot
	       to own plot
	    */
	    graph = Tcl_GetGraphFromObj(graph_obj);
	    graph->dim.y1 = canvas_height(interp, e->win);
	    s_result->pr_func(s_result, (seq_reg_plot *)jdata);
	    break;
	}
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case RESULT:
	    jdata->info.result = (void *)s_result;
	    break;
	case WIN_NAME:
	    {
		char *r_win = e->win;
		jdata->info.result = (void *)r_win;
		break;
	    }
	case WIN_SIZE:
	    {
		static d_point pt;
		pt.x = get_default_int(interp, tk_utils_defs,
					w("ELEMENT.PLOT_WIDTH"));
		pt.y = get_default_double(interp, tk_utils_defs,
					   w("ELEMENT.SINGLE.PLOT_HEIGHT"));

		jdata->info.result = (void *)&pt;
		break;
	    }
	}
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: 
	{
	    nip_stop_codons_shutdown(interp, s_result, e, seq_num);
	    break;
	}
    }
}

int
NipFindStopCodons(int strand,                                          /* in */
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
    int *matches, max_matches=MAXMATCHES;
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
    
    if (strand & BOTTOM_S) {
	start_codon = num_codons;
    }
    /* if ((strcmp(strand, "-") == 0) || (strcmp(strand, "both") == 0)) { */
    if (strand & BOTTOM_S) {
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
    seq_result *s_result = (seq_result *) obj;
    Tcl_Obj *graph_obj = s_result->data;
    int i;
    Graph *graph;

    graph = Tcl_GetGraphFromObj(graph_obj);

    for (i = 0; i < graph->d_arrays[0].n_dlines; i++) {
	UpdateTextOutput();
	vmessage("Position %10f\n", graph->d_arrays[0].d_array[i].x0);
    }
}

int
nip_stop_codons(char *sequence,
		int sequence_type,
		int start, 
		int end,
		int strand, 
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
		     int strand, 
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
static int compare_darray(const void *p1,
			  const void *p2)
{
    gd_line *i1 = (gd_line *) p1, *i2 = (gd_line *) p2;
    return i1->x0 - i2->x0;
}

int store_stop_codons(int seq_num,
		      in_s_codon *input,
		      int start,
		      int end,
		      int *stop_codon,
		      int num_stop_codons,
		      int frame,
		      int type,
		      Tcl_Obj **graph_obj)
{
    seq_result *s_result;
    Graph *graph;
    int id, i;

    if (NULL == (s_result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (graph = (Graph *)xmalloc(sizeof(Graph))))
	return -1;

    Tcl_InitGraph(&graph);

    if (NULL == (graph->d_arrays = (darray *)xmalloc(sizeof(darray))))
	return -1;
    if (NULL == (graph->d_arrays[0].d_array = (gd_line *)xmalloc(num_stop_codons * sizeof(gd_line))))
	return -1;

    graph->d_arrays[0].n_dlines = num_stop_codons;
    graph->n_darrays = 1;
    graph->d_arrays[0].type = G_LINES;
    graph->dim.x0 = start;
    graph->dim.x1 = end;
    
    for (i = 0; i < num_stop_codons; i++) {
	graph->d_arrays[0].d_array[i].x0 = stop_codon[i];
	graph->d_arrays[0].d_array[i].x1 = stop_codon[i]; 
    }    

    /* must sort the array on x0 */
    qsort((void *) graph->d_arrays[0].d_array, graph->d_arrays[0].n_dlines,
	  sizeof(gd_line), compare_darray);

    *graph_obj = Tcl_NewGraphObj(graph);

    /* these are copied over in Tcl_NewGraphObj so don't need them anymore */
    xfree(graph->d_arrays[0].d_array);
    xfree(graph->d_arrays);
    xfree(graph);

    id = get_reg_id();

    s_result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    s_result->seq_id[VERTICAL] = -1;

    s_result->id = id; 
    s_result->input = (void *)input; 
    s_result->frame = frame;
    s_result->graph = GRAPH;
    s_result->data = *graph_obj; 
    s_result->e = NULL; /* initialise */

    s_result->pr_func = seq_plot_graph_func;
    s_result->op_func = nip_stop_codons_callback;
    s_result->txt_func = nip_stop_codons_text_func;

    if (type == STOP_CODON) {
	s_result->type = SEQ_TYPE_STOPCODON;
    } else {
	s_result->type = SEQ_TYPE_STARTCODON;
    }
    s_result->gr_type = SEQ_TYPE_GRAPH_PLOT;

    seq_register(seq_num, nip_stop_codons_callback, (void *)s_result, 
		 SEQ_PLOT_PERM, id);

    xfree(stop_codon);
    return id;
}

int init_nip_stop_codons_create(int seq_id,
				int start,
				int end,
				int strand_sym,
				Tcl_Obj **graph_obj,
				int *id)
{
    in_s_codon *input1;
    in_s_codon *input2;
    in_s_codon *input3;
    char *seq;
    int seq_len;
    int sequence_type;
    int seq_num;
    s_codon_res *stop_codon;
    int type = STOP_CODON;
    Tcl_DString input_params;
    char strand[8];

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

    nip_stop_codons(seq, sequence_type, start, end, strand_sym,
		    stop_codon);

    /* create inputs parameters */
    Tcl_DStringInit(&input_params);

    if (strand_sym & TOP_S) {
        strcpy(strand, "top");
    } else if (strand_sym & BOTTOM_S) {
        strcpy(strand, "bottom");
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

    if (stop_codon->n_match1 > 0) {
	if (-1 == (id[0] = store_stop_codons(seq_num, input1, start, end,
					     stop_codon->match1, 
					     stop_codon->n_match1,
					     1, type, &graph_obj[0]))){
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
					     stop_codon->n_match2, 2, type,
					     &graph_obj[1]))){
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
					     stop_codon->n_match3, 3, type,
					     &graph_obj[2]))){
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

    xfree(stop_codon);
    return 0;
}

int init_nip_start_codons_create(int seq_id,
				 int start,
				 int end,
				 int strand_sym,
				 Tcl_Obj **graph_obj,
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

    if (strand_sym & TOP_S) {
        strcpy(strand, "forward");
    } else if (strand_sym & BOTTOM_S) {
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
					     1, type, &graph_obj[0]))){
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
					     2, type, &graph_obj[1]))){
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
					     3, type, &graph_obj[2]))){
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
			      int seq_id,
			      int result_id,
			      char *e_win,
			      char *c_win, 
			      Tcl_Obj *results,
			      int container_id,
			      int element_id,
			      char *element_type,
			      int line_width,
			      char *colour,
			      float tick_ht,
			      int orientation)
{
    seq_result *s_result;
    Graph *graph;
    configs *configure;
    plot_data *result;
    int i;
    int height;
    int seq_id_h, seq_id_v;

    s_result = seq_id_to_result(result_id);
    
    if (NULL == (result = (plot_data *)xmalloc(sizeof(plot_data))))
	return -1;
    
    if (NULL == (result->configure = (configs **)xmalloc(sizeof(configs*))))
	return -1;
    
    if (NULL == (configure = (configs *)xmalloc(sizeof(configs))))
	return -1;

    configure->position = 0.5;
    configure->x_direction = '+';
    configure->y_direction = '+';
    configure->height = tick_ht;
    configure->zoom = 0;
    configure->scroll = 0;

    result->configure[0] = configure;
    result->n_configure = 1;
    result->sf_m = 1.0;
    result->sf_c = 0.0;
    result->result_id = result_id;
    result->hidden = 0;
    result->line_width = line_width;
    result->colour = strdup(colour);
    result->len_ruler = 1;
    result->amp_ruler = 0;
    sprintf(result->tags, "id%d", result_id);

    graph = Tcl_GetGraphFromObj(results);
   
    if (orientation == HORIZONTAL) {
	seq_id_h = seq_id;
	seq_id_v = -1;
	height = canvas_height(interp, e_win);
	result->scale = SCALE_X;
    } else {
	seq_id_h = -1;
	seq_id_v = seq_id;
	height = canvas_width(interp, e_win);
	result->scale = SCALE_Y;
    }

    graph->dim.y0 = 1;
    graph->dim.y1 = height;

    /* graph->dim.y1 = tick_ht; */

    for (i = 0; i < graph->d_arrays[0].n_dlines; i++) {
	graph->d_arrays[0].d_array[i].y0 = 1;
	graph->d_arrays[0].d_array[i].y1 = tick_ht;

    }
    
    init_seq_element(interp, s_result, seq_id_h, seq_id_v, result, 
		     container_id, 
		     element_id, e_win, c_win, orientation, orientation,
		     graph, element_type);
    
    return 0;
}

