#include <limits.h>
#include <string.h>

#include "seq_results.h"
#include "nip_results.h"
#include "restriction_enzyme_map.h"
#include "nip_restriction_enzymes.h"
#include "nip_globals.h"
#include "tcl_utils.h"
#include "misc.h"
#include "text_output.h"
#include "renz_utils.h"
#include "nip_canvas_box.h"

void nip_renz_callback(int seq_num, void *obj, seq_reg_data *jdata);

void nip_renz_shutdown(Tcl_Interp *interp,
		       seq_result *result,
		       int seq_num)
{
    renz_res *data = result->data;
    out_canvas *output = result->output;
    char cmd[1024];
    char *tmp; 
    
    /* need to deregister sequence */
    seq_deregister(seq_num, nip_renz_callback, (seq_result *)result);
    delete_cursor(seq_num, output->cursor->id, 0);
 
    tmp = get_default_string(interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
    if (TCL_OK != Tcl_VarEval(interp, 
			      "seq_result_list_update ", 
			      tmp, NULL)){
	verror(ERR_WARN, "restriction enzymes", "shutdown %s \n", 
	       Tcl_GetStringResult(interp));
    }
    sprintf(cmd, "DeleteREnzPlot %s %s\n", data->frame, data->re_win);
    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	verror(ERR_WARN, "restriction enzymes", "shutdown %s\n", 
	       Tcl_GetStringResult(interp));
    }

    renz_shutdown(data->r_enzyme, data->num_enzymes, data->match,
		  data->canvas, data->world, data->zoom);
    free(data->text_colour);
    free(data->tick->colour);
    free(data->cursor.colour);
    free(data->ruler->tick.t.colour);
    free(data->ruler->colour);
    xfree(data->ruler);
    free_win_list(data->win_list, data->num_wins);

    xfree(result->data);
    /* xfree(result->input); */
    xfree(result->output);
    xfree(result);
}

void nip_renz_info(int seq_num,
		   renz_res *data,
		   int start,
		   int print_opt)
{
    int seq_len;
    char *seq;
    int sequence_type;
    
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    sequence_type = GetSeqStructure(seq_num);

    vfuncheader("Restriction enzymes result list");
    
    vmessage("Sequence %s\n", GetSeqName(seq_num));

    vmessage("Number of enzymes = %d\n", data->num_enzymes);
    vmessage("Number of matches = %d\n", data->num_match);

    if (print_opt == 0) {
	start_message();
	PrintEnzymeByEnzyme(data->r_enzyme, data->match, data->num_match, 
			    data->num_enzymes, seq, seq_len, 
			    sequence_type, start, 1);
	end_message(data->re_win);
    } else {
	start_message();
	OrderOnPosition(data->r_enzyme, data->match, data->num_match,
			seq, seq_len, sequence_type, start);
	end_message(data->re_win);
    }
}

void
nip_renz_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    seq_result *result = (seq_result *) obj;
    renz_res *data = result->data;

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Restriction enzyme map");
	break;    

    case SEQ_GET_OPS:
	jdata->get_ops.ops = "Output enzyme by enzyme\0"
	    "Output ordered on position\0"
		"SEPARATOR\0"
	     "Remove\0";
	break;
    case SEQ_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0: /* Information print enzyme by enzyme */
	    {
		int seq_num;

		seq_num = GetSeqNum(result->seq_id[0]);
		nip_renz_info(seq_num, data, data->ruler->start, 0);
		break;
	    }
	case 1: /* Information order */
	    {
		int seq_num;

		seq_num = GetSeqNum(result->seq_id[0]);
		nip_renz_info(seq_num, data, data->ruler->start, 1);
		break;
	    }
	case 2: /* Remove */
	    {
		out_canvas *output = result->output;
		nip_renz_shutdown(output->interp, result, seq_num);
		break;
	    }
	}
	break;
    case SEQ_PLOT:
	result->pr_func(result, (seq_reg_plot *)jdata);
	break;
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
#ifdef DEBUG
	case OUTPUT:
	    jdata->info.result = (void *)output;
	    break;
	case INPUT: 
	    jdata->info.result = (void *)input;
	    break;
	case NPTS: 
	    jdata->info.result = (void *)num_pts;
	    break;
	case INDEX: 
	    jdata->info.result = (void *)id;
	    break;
#endif
	case RESULT:
	    jdata->info.result = (void *)result;
	    break;
	} /* switch SEQ_RESULT_INFO */
	break;
#ifdef DEBUG
    case SEQ_HIDE: 
	output->hidden = 1;
	break;
    case SEQ_REVEAL: 
	output->hidden = 0;
	break;
#endif
    case SEQ_QUIT:
    case SEQ_DELETE: 
	{
	    out_canvas *output = result->output;
	    nip_renz_shutdown(output->interp, result,seq_num);
	    break;
	} /* SEQ_DELETE */
    case SEQ_GENERIC:
	switch (jdata->generic.task) {
	case TASK_NIP_RENZ_INFO: 
	    {
		/* info on all cuts by a specific restriction enzyme */
		renz_res *r;
		int i;
	        R_Match *tmp_match;
		int cnt = 0;
		int *enz_name = (int *)jdata->generic.data;
		int seq_len;
		char *seq;
		int sequence_type;
		int seq_num;

		r = result->data;
#ifdef DEBUG
		printf("TASK_RENZ_INFO name %d\n", *enz_name);
#endif
		seq_num = GetSeqNum(result->seq_id[0]);
		seq = GetSeqSequence(seq_num);
		seq_len = GetSeqLength(seq_num);
		sequence_type = GetSeqStructure(seq_num);

		if (NULL == (tmp_match = (R_Match *)malloc(r->num_match * 
							   sizeof(R_Match))))
		    return;
		/* 
		 * make tmp match array containing only those matches for the 
		 * selected enzyme
		 */
		for (i = 0; i < r->num_match; i++) {
		    if (*enz_name == r->match[i].enz_name) {
			tmp_match[cnt++] = r->match[i];
		    }
		}
		start_message();
		PrintEnzymeByEnzyme(r->r_enzyme, tmp_match, cnt,
				    r->num_enzymes, seq, seq_len, 
				    sequence_type, r->ruler->start, 0);
		end_message(r->re_win);
		xfree(tmp_match);
		break;
	    }
	}
	break; /* SEQ_GENERIC */
    case SEQ_CURSOR_NOTIFY:
	{
	    renz_res *r = result->data;
	    out_canvas *output = result->output;
#ifdef DEBUG	    
	    printf("SEQ_CURSOR_NOTIFY: restriction enzymes\n");
#endif
	    nip_canvas_cursor_refresh(output->interp, GetSeqId(seq_num),
				      jdata->cursor_notify.cursor, 
				      output->cursor, r->canvas, r->win_list, 
				      r->num_wins, result->id, 
				      &output->cursor_visible, r->world, 1);
	    break;
	}
    } /* switch */
}

void
nip_renz_print_func(void *obj, seq_reg_plot *plot)
{
    seq_result *result = (seq_result *) obj;
    out_canvas *output = result->output;
    renz_res *data = result->data;

#ifdef DEBUG_START
    printf("START nip_renz_print_func\n");
#endif

    plot_renz_matches(output->interp, data->re_win, data->names_win,
		      data->text_offset, data->text_colour,
		      data->yoffset, data->num_enzymes, data->r_enzyme,
		      data->ruler, data->sequence_len,
		      data->num_match, data->match,
		      data->tick, data->frame, data->world,
		      data->canvas, data->win_list, data->num_wins,
		      &data->zoom);


}

void
nip_renz_text_func(void *obj)
{

}

int nip_renz_reg(Tcl_Interp *interp,
		 int seq_id,
		 out_canvas *output,
		 char *filename,
		 char *frame,
		 char *names_win,
		 char *re_win,
		 char *inlist,
		 int num_items,
		 int start,
		 int end,
		 int text_offset,
		 char *text_fill,
		 tick_s *tick,
		 int yoffset,
		 ruler_s *ruler,
		 cursor_s cursor)
{
    seq_result *result;
    int num_enzymes;
    R_Enz *r_enzyme;
    renz_res *data;
    int id, i;
    R_Match *match;
    int total_matches;
    char *seq;
    int seq_len;
    int seq_num;
    int sequence_type;
    seq_cursor_notify cn;
    int line_width;
    int max_overlap;
    int new_start, new_end;

    if (NULL == (result = (seq_result *)xmalloc(sizeof(seq_result))))
	return -1;

    if (NULL == (data = (renz_res *)xmalloc(sizeof(renz_res))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);
    sequence_type = GetSeqStructure(seq_num);

    result->data = data;

    result->seq_id[HORIZONTAL] = GetSeqId(seq_num);
    result->seq_id[VERTICAL] = -1;

    id = get_reg_id();

    result->id = id; 

    /* nip_result->input = (void *)input; */
    result->output = (void *)output;

    result->pr_func = nip_renz_print_func;
    result->op_func = nip_renz_callback;
    result->txt_func = nip_renz_text_func;

    /* data structure */
    strcpy(data->re_win, re_win);
    strcpy(data->frame, frame);
    strcpy(data->names_win, names_win);
    
    data->tick = tick;
    data->ruler = ruler;
    data->sequence_len = GetSeqLength(seq_num);
    data->cursor = cursor;

    /* create list of windows in the restriction enzyme display */
    if (NULL == (data->win_list = (win **)xmalloc(MAX_NUM_WINS * sizeof(win*))))
	return -1;

    data->num_wins = 0;
    addWindow(data->win_list, &data->num_wins, data->re_win, 'b', id);
    addWindow(data->win_list, &data->num_wins, data->ruler->window, 'x', id);
    addWindow(data->win_list, &data->num_wins, data->names_win, 'y', id);

    if (NULL == (data->canvas = (CanvasPtr *)xmalloc(sizeof(CanvasPtr))))
	return -1;

    if (NULL == (data->world= (WorldPtr *)xmalloc(sizeof(WorldPtr))))
	return -1;

    if (NULL == (data->world->visible = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    if (NULL == (data->world->total = (d_box *)xmalloc(sizeof(d_box))))
	return -1;

    initCanvas(interp, data->canvas, data->re_win);
    createZoom(&data->zoom);

    /*
     * open the file of enzyme names and parse the selected enzymes into the
     * r_enzyme array of size num_enzymes
     */
    open_renz_file(filename, inlist, num_items, &r_enzyme, &num_enzymes);
    max_overlap = find_max_cut_dist(r_enzyme, num_enzymes);

    new_start = MAX(1, start - max_overlap);
    new_end = MIN(seq_len, end + max_overlap);

    seq = &seq[new_start-1];
    seq_len = new_end - new_start + 1;

    data->r_enzyme = r_enzyme;
    data->num_enzymes = num_enzymes;
    data->tick = tick;
    data->yoffset = yoffset;
    data->sequence_type = sequence_type;
    data->text_offset = text_offset;
    data->text_colour = strdup(text_fill);

    if (NULL == (match = (R_Match*)xcalloc(MAXMATCHES, sizeof(R_Match))))
	return -1;

    FindMatches(r_enzyme, num_enzymes, seq, seq_len,
		sequence_type, &match, &total_matches);

    for (i = 0; i < total_matches; i++) {
	match[i].cut_pos = match[i].cut_pos - (start - new_start);
    }
    data->match = match;
    data->num_match = total_matches;
    
    line_width = get_default_int(interp, nip_defs, w("NIP.CURSOR.LINE_WIDTH"));

    output->cursor = create_cursor(seq_num, 0, NULL, line_width, 1, 
				       HORIZONTAL);
    output->cursor_visible = 0;

    /* move cursor to start position if this is our own cursor */
    if (output->cursor->refs == 1)
	output->cursor->abspos = start;
    
    seq_register(seq_num, nip_renz_callback, (void *)result, SEQ_PLOT_PERM, id);

    plot_renz_matches(output->interp, data->re_win, data->names_win,
		      data->text_offset, data->text_colour, 
		      data->yoffset, data->num_enzymes, data->r_enzyme,
		      data->ruler, data->sequence_len,
		      data->num_match, data->match,
		      data->tick, data->frame, data->world,
		      data->canvas, data->win_list, data->num_wins,
		      &data->zoom);

    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = output->cursor;
    cn.cursor->job = CURSOR_MOVE;
    seq_notify(seq_num, (seq_reg_data *)&cn); 
    
    if (total_matches == 0) {
      nip_renz_shutdown(interp, result, seq_num);
      return -1;
    }
    return id;
}

