#include <tcl.h>
#include <stdlib.h>
#include <string.h>

#include "seq_reg.h"
#include "misc.h"
#include "seq_results.h"
#include "tcl_utils.h"

typedef struct { 
    char *rid;
    Tcl_Interp *interp;
    int communicating;
} sender_result;
void sender_callback(int seq_num, void *obj, seq_reg_data *jdata);

void sender_shutdown(int seq_num,
		     sender_result *send)
{
    char cmd[1024];
    char *tmp;

    sprintf(cmd, "upvar #0 commn_[list %s] commn; "
	    "eval $commn(command) EventHandler {{{%s}}} "
	    "STOP_SEQUENCE",
	    send->rid, send->rid);

    send->communicating = 1;
    seq_deregister(seq_num, sender_callback, (sender_result *)send);

    if (TCL_ERROR == Tcl_Eval(send->interp, cmd)) 
	verror(ERR_WARN, "sender_shutdown", "%s\n", Tcl_GetStringResult(send->interp));
    Tcl_VarEval(send->interp, "unset commn", NULL);
    send->communicating = 0;

    /* check if no more sequences - shutdown */
#ifdef REMOVE
    if (NumSequences() == 0) {
	if (TCL_ERROR == (Tcl_VarEval(send->interp, "ExitSip", NULL)))
	    verror(ERR_WARN, "sender_shutdown",  "%s\n", Tcl_GetStringResult(send->interp));
    }
#endif

    tmp = get_default_string(send->interp, tk_utils_defs, 
			     w("RASTER.RESULTS.WIN"));
    if (TCL_OK != Tcl_VarEval(send->interp, "seq_result_list_update ", 
			      tmp, NULL)){
	verror(ERR_WARN, "sender shutdown", "%s \n", Tcl_GetStringResult(send->interp));
    }
}

void sender_callback(int seq_num, void *obj, seq_reg_data *jdata)
{
    sender_result *send = (sender_result *)obj;
    
    if (send->communicating) {
	return;
    }

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	sprintf(jdata->name.line, "Send to Gap4, %s", GetSeqName(seq_num));
	break;    

    case SEQ_GET_OPS:
	jdata->get_ops.ops = "Information\0Remove\0";
	break;
    case SEQ_INVOKE_OP:
	switch (jdata->invoke_op.op) {
	case 0: /* information */
	    {
		char cmd[1024];
		char *command;
		vfuncheader("input parameters");
		sprintf(cmd, "upvar #0 commn_%s commn;", send->rid);
		Tcl_Eval(send->interp, cmd);
		command = Tcl_GetVar2(send->interp, "commn", "command", 
				      TCL_GLOBAL_ONLY);
		vmessage("Sequence %s\nCommand \"%s\"\n", 
			 GetSeqName(seq_num), command);
		break;
	    }
	case 1: /* remove */ 
	    {
#ifdef FIXME
		SipDeleteSequence(send->interp, seq_num);
#endif
		/* delete_cursor -seq_id $commn(seq_id) -id $commn(cursor_id) */
		sender_shutdown(seq_num, send);
		break;
	    }
	}
	break;
    case SEQ_QUIT:
    case SEQ_DELETE: 
	{
	    sender_shutdown(seq_num, send);

	    break;
	} /* SEQ_DELETE, SEQ_QUIT */
    case SEQ_CURSOR_NOTIFY:
	{
	    char cmd[1024];
	    cursor_t *cursor = (cursor_t *)jdata->cursor_notify.cursor;
	    char id[100];
	    char c[1024];
	    char job[1024];
	    int first = 1;
	    int num;
	    char **list;

	    *job = 0;
	    strcat(job, "{");
	    if (jdata->cursor_notify.cursor->job & CURSOR_MOVE) {
		strcat(job, "MOVE");
		first = 0;
	    }
	    if (jdata->cursor_notify.cursor->job & CURSOR_INCREMENT) {
		strcat(job, first ? "INCREMENT" : " INCREMENT");
		first = 0;
	    }
	    if (jdata->cursor_notify.cursor->job & CURSOR_DECREMENT) {
		strcat(job, first ? "DECREMENT" : " DECREMENT");
		first = 0;
	    }
	    if (jdata->cursor_notify.cursor->job & CURSOR_DELETE) {
		strcat(job, first ? "DELETE" : " DELETE");
		first = 0;
	    }
	    strcat(job, "}");

	    sprintf(cmd, "upvar #0 commn_[list %s] commn; "
		    "eval $commn(command) EventHandler {{{%s}}} "
		    "CURSOR_NOTIFY [list {{id %d} {pos %d} {seq 0} "
		    "{abspos %d} {job %s}}]",
		    send->rid, send->rid, cursor->id,
		    cursor->abspos, cursor->abspos, job);

	    send->communicating = 1;
	    Tcl_SetVar2(send->interp, "communicating", send->rid, "1",
			TCL_GLOBAL_ONLY);
	    if (TCL_ERROR == Tcl_Eval(send->interp, cmd)) 
		verror(ERR_WARN, "sender_callback", "%s\n", 
		       Tcl_GetStringResult(send->interp));

#ifdef DEBUG
	    printf("SIP CURSORS %s\n", Tcl_GetStringResult(send->interp));
#endif
	    if (strcmp(Tcl_GetStringResult(send->interp), "") != 0) {
		sprintf(id, "%d", cursor->id);
		if (cursor->direction == HORIZONTAL) {
		    /* sip */
		    sprintf(c, "cursor_h_%s", send->rid); 
		} else if (cursor->direction == VERTICAL)  {
		    sprintf(c, "cursor_v_%s", send->rid);
		} else {
		    /* nip */
		    sprintf(c, "cursor_%s", send->rid); 
		}
		if (Tcl_SplitList(send->interp, Tcl_GetStringResult(send->interp), &num, &list) 
		    != TCL_OK) {
		    return;
		}
#ifdef DEBUG
		printf("cid %d refs %d\n", atoi(list[0]), atoi(list[1]));
#endif
		Tcl_SetVar2(send->interp, c, list[0], id, TCL_GLOBAL_ONLY);
#ifdef DEBUG
		printf("SIP refs %d GAP refs %d\n", cursor->refs, atoi(list[1]));
#endif
		if (atoi(list[1]) > cursor->refs) {
		    cursor->refs = atoi(list[1]);
		}
		Tcl_Free((char*)list);
#ifdef DEBUG
		printf("SENDTO %s %s %s refs %d\n", c, Tcl_GetStringResult(send->interp),
		       id, cursor->refs);
#endif
	    }

	    send->communicating = 0;
	    Tcl_SetVar2(send->interp, "communicating", send->rid, "0",
			TCL_GLOBAL_ONLY);
	    break;
	}
    }

}


int seq_sender(Tcl_Interp *interp,
	       char *rid,
	       int seq_id)
{
    int seq_num;
    int id;
    sender_result *send;

    if (NULL == (send = (sender_result *)xmalloc(sizeof(sender_result))))
	return -1;

    seq_num = GetSeqNum(seq_id);
    id = get_reg_id();
    send->rid = strdup(rid);
    send->interp = interp;
    send->communicating = 0;
    seq_register(seq_num, sender_callback, (void *)send, SEQ_SENDER, id);

    return id;
}
