#include <float.h>
#include <limits.h>

#include "seq_reg.h"
#include "xalloc.h"
#include "tkSeqed.h"
#include "seq_results.h"
#include "seqed.h"
#include "tkSeqedUtils.h"
#include "tcl_utils.h"
#include "seqed_restriction_enzymes.h"
#include "text_output.h"
#include "spin_globals.h"
#include "seq_element.h"

int seqed_reg(Tcl_Interp *interp, char *seqed_win, int seq_num, tkSeqed *se,
	      int pos, int container_id, char *c_win, int element_id);

void seqed_shutdown(Tcl_Interp *interp,
		    SeqedResult *result,
		    element *e)
{
    Tcl_CmdInfo info;
    tkSeqed *se;

#ifdef DEBUG
    printf("seqed shutdown \n");
#endif

    Tcl_GetCommandInfo(interp, result->seqed_win, &info);
    se = (tkSeqed*)info.clientData;

   if (se->renzDisplayed) {
	free_lines();
	free_r_enzyme(se->r_enzyme, se->num_enzymes);
    }
   Tcl_VarEval(e->c->interp, "result_list_update ", e->c->win, NULL);

    xfree(result);
}

tkSeqed *seqed_id_to_se(Tcl_Interp *interp,
			 int seqed_id)
{
    Tcl_CmdInfo info;
    tkSeqed *se;
    seq_reg_info info1;
    char *seqed_win;

    info1.job = SEQ_RESULT_INFO;
    info1.op = WIN_NAME;
    info1.result = NULL;

    seq_result_notify(seqed_id, (seq_reg_data *)&info1, 0);
    seqed_win = (char *)info1.result;
    Tcl_GetCommandInfo(interp, seqed_win, &info);
    se = (tkSeqed*)info.clientData;
    return(se);
}

/*
 * add sequence to seqed widget
 */
int add_seq_seqed(Tcl_Interp *interp,
		  char *sequence,
		  char *seqed_win,
		  int seq_num,
		  int pos,
		  int container_id,
		  char *c_win,
		  int element_id)
{

    Tcl_CmdInfo info;
    tkSeqed *se;
    char *seq_name;
    int sequence_type;
    int seqed_id;

    Tcl_GetCommandInfo(interp, seqed_win, &info);
    se = (tkSeqed*)info.clientData;

    seq_name = GetSeqName(seq_num);
    sequence_type = GetSeqStructure(seq_num);

    seqed_add_sequence(se, strlen(sequence), sequence, seq_name,
		       sequence_type, GetSeqId(seq_num), 0, 0);
    seqed_id = seqed_reg(interp, seqed_win, seq_num, se, pos, container_id,
			 c_win, element_id);
    return seqed_id;
}

/*
 * update sequence position
 */
void update_seqed(Tcl_Interp *interp,
		  char *seqed_win,
		  int pos)
{

    Tcl_CmdInfo info;
    tkSeqed *se;

    Tcl_GetCommandInfo(interp, seqed_win, &info);
    se = (tkSeqed*)info.clientData;

    /* seqed_redisplay_seq(se, pos); */
    seqed_setCursorPos(se, pos);
}

void seqed_move_cursor(Tcl_Interp *interp,
		       char *seqed_win,
		       int pos)
{
    Tcl_CmdInfo info;
    tkSeqed *se;

    /* printf("seqed move cursor %s %d %d \n", seqed_win, seqed_id, pos); */

    Tcl_GetCommandInfo(interp, seqed_win, &info);
    se = (tkSeqed*)info.clientData;

    se->cursorPos = pos;
    seqed_showCursor(se, se->cursorSeq, se->cursorPos);
}


void seqed_set_cursor_pos(Tcl_Interp *interp,
			  int seqed_id,
			  int pos)
{
    tkSeqed *se;

    se = seqed_id_to_se(interp,seqed_id);
    seqed_showCursor(se, se->cursorSeq, pos);
    seqed_setCursorPos(se, pos);
}


void
seqed_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    SeqedResult *result = (SeqedResult *) obj;
    char *seqed_win = result->seqed_win;
    tkSeqed *se;
    Tcl_CmdInfo info;
    element *e = result->e;

    Tcl_GetCommandInfo(result->interp, seqed_win, &info);
    se = (tkSeqed*)info.clientData;

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "sequence editor");
	    return;
	}
    case SEQ_QUIT:
    case SEQ_DELETE:
	{
	    seq_deregister(seq_num, seqed_callback, (SeqedResult *)result);

	    /* need to set prev_pos */

	    se->prev_pos = se->cursor->abspos;
	    delete_cursor(seq_num, se->cursor->id, 1);
	    seqed_shutdown(result->interp, result, e);
	    return;
	}
    case SEQ_GET_OPS:
	{
	    jdata->get_ops.ops = "Remove\0";
	    break;
	}

    case SEQ_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0:
	    printf("seqed Remove menu\n");

	    /* need to set prev_pos before shutdown which frees se */
	    se->prev_pos = se->cursor->abspos;
	    delete_cursor(seq_num, se->cursor->id, 1);

	    seq_register_dump();
	    seq_deregister(seq_num, seqed_callback, (SeqedResult *)result);
	    seq_register_dump();

	    seqed_shutdown(result->interp, result, e);
	    remove_result_from_element(e, result->index);
	    return;
	}
	break;
    case SEQ_CURSOR_NOTIFY:
	{
	    cursor_t *cursor = (cursor_t *)jdata->cursor_notify.cursor;

#ifdef DEBUG
	    printf("SEQ_CURSOR_NOTIFY seqed pos %d id %d result id %d\n",
		   cursor->abspos, cursor->id, se->cursor->id);
#endif
	    if (se->cursor->id == cursor->id) {
		seqed_move_cursor(result->interp, result->seqed_win,
				  cursor->abspos);
	    }
	    break;
	}
    case SEQ_SETCURSOR:
	{
	    seqed_move_cursor(result->interp, result->seqed_win,
			      jdata->cursor_move.pos);

	    break;
	}
    case SEQ_GENERIC:
	switch(jdata->generic.task) {
	case TASK_SEQED_SETCURSOR:
	    {
		task_seqed_setcursor *tc =
		    (task_seqed_setcursor *)jdata->generic.data;

		seqed_move_cursor(result->interp, result->seqed_win,
				      tc->pos);

		break;
	    }
	case TASK_SEQED_GETCURSOR:
	    {
		Tcl_CmdInfo info;
		tkSeqed *se;
		task_seqed_setcursor *tc;

		if (NULL == (tc = (task_seqed_setcursor *)xmalloc
			     (sizeof(task_seqed_setcursor))))
		    return;

		Tcl_GetCommandInfo(result->interp, result->seqed_win, &info);
		se = (tkSeqed*)info.clientData;

		tc->pos = se->cursorPos;
		tc->col = se->cursorCol;
		jdata->info.result = (void *)tc;
		break;
	    }
	}
	break;

    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case WIN_NAME:
	    jdata->info.result = (void *)seqed_win;
	    break;
	}
	break;
    case SEQ_SEQUENCE_TYPE:
	{
	    int *type = (int *)jdata->sequence_type.data;

	    se->sequence_type = *type;
	    break;
	}
    }
}

int seqed_reg(Tcl_Interp *interp,
	      char *seqed_win,
	      int seq_num,
	      tkSeqed *se,
	      int pos,
	      int container_id,
	      char *c_win,
	      int element_id)
{
    int id;
    SeqedResult *seqed_result;
    char *colour;
    seq_cursor_notify cn;
    int line_width;
    element *e;
    plot_data *result;

    if (NULL == (seqed_result = (SeqedResult *)xmalloc(sizeof(SeqedResult))))
	return -1;

   if (NULL == (result = (plot_data *)xmalloc(sizeof(plot_data))))
	return -1;

    seqed_result->op_func = seqed_callback;
    seqed_result->seq_id = GetSeqId(seq_num);
    strcpy(seqed_result->seqed_win, seqed_win);
    seqed_result->interp = interp;

    id = get_reg_id();
    seqed_result->index = id;

    line_width = get_default_int(interp, spin_defs,
				 w("SEQ.CURSOR.LINE_WIDTH"));

    se->rid = id;
    se->seq_id = GetSeqId(seq_num);
    colour = get_cursor_colour(se->rid);
    strcpy(se->cursorCol, colour);
    se->cursor = create_cursor(seq_num, 1, NULL, line_width, 1, HORIZONTAL, pos);
    se->cursor_visible = 0;

    se->prev_pos = se->cursor->abspos;
    se->cursor->abspos = se->cursorPos;


    if (NULL == (e = (get_element(element_id)))) {
	e = create_element(interp, container_id, element_id, seqed_win, HORIZONTAL, -1);
    }
    seqed_result->e = e;

    result->result_id = id;
    result->amp_ruler = 0;
    result->len_ruler = 0;
    result->n_configure = 0;
    result->sf_m = 1.0;
    result->sf_c = 0.0;
    result->scale = 0;
    result->hidden = 0;
    result->colour = NULL;
    result->cursor = NULL;

    set_seq_element_funcs(e);

    add_seq_id_to_element(e, GetSeqId(seq_num), HORIZONTAL);

    add_result_to_element(e, result, DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX, HORIZONTAL, SEQ_EDITOR);

    add_element_to_container(interp, container_id, c_win, e,
			     INT_MAX, INT_MIN, INT_MAX, INT_MIN);

    seq_register(seq_num, seqed_callback, (void *)seqed_result, SEQ_SEQED,
		 id);

    seq_element_reg(interp, e);

    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = se->cursor;
    cn.cursor->job = CURSOR_MOVE;
    cn.show_all = 0;
    seq_notify(seq_num, (seq_reg_data *)&cn);

    return id;
}


