#include <staden_config.h>

#include <tcl.h>
#include <stdlib.h>
#include <string.h>

#include "seq_results.h"
#include "seq_reg.h"
#include "seq_raster.h"
#include "xalloc.h"
#include "seq_reg_structs.h"
#include "cli_arg.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "sequence_formats.h"
#include "array.h"
#include "seq_sendto.h"
#include "spin_globals.h"
#include "emboss_input_funcs.h"
#include "limits.h"
#include "sequence_pair_display.h"
#include "seqed.h"
#include "sequtils.h"
#ifdef USE_SEQLIB
#include "seq_cmds.h"
#endif
#include "sip_cmds.h"
#include "nip_cmds.h"
#include "raster_cmds.h"
#include "misc.h"
#include "renz_utils.h"

static int seq_file_save_ft(FILE *pw, Featcds **key_index, int start, int end,
			    char *seq,  char *identifier);

/*
 * comparison function for qsort to sort the result names into id order
 */
int
compare_rnames(const void *vdata1, const void *vdata2) {
    seq_reg_name *data1 = (seq_reg_name *)vdata1;
    seq_reg_name *data2 = (seq_reg_name *)vdata2;

    if ((*data1).id < (*data2).id)
	return -1;
    else if ((*data1).id == (*data2).id)
	return 0;
    else
	return 1;

}

/*
 * return a list of available functions to tcl results manager
 * the order of the results in order of creation (on unique id number)
 */
int 
tcl_seq_result_names(ClientData clientData, 
		     Tcl_Interp *interp, 
		     int argc, 
		     char *argv[]) 
{    
    int i;
    char buf[1024];
    int num_elements, num_results;
    seq_reg_name *data;
    seq_result_names_arg args;
    char *win_name;
    seq_reg_info info;
    RasterResult *raster_result;

    cli_args a[] = {
	{"-raster_id", ARG_INT, 1, "-1", offsetof(seq_result_names_arg, raster_id)},
	{"-result_id", ARG_INT, 1, "-1", offsetof(seq_result_names_arg, result_id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (data = seq_result_names(&num_elements))) {
	/* if there are no registered results */
	return TCL_OK;
    }

    num_results = seq_num_results();
    if (args.raster_id == -1 && args.result_id == -1) {
	/* return list of all result names ie for results manager */
	/* order results according to id number */
	qsort((void *) data, num_elements, sizeof(seq_reg_name ), compare_rnames);
	
	Tcl_ResetResult(interp);
	for (i = 0; i < num_elements; i++) {
	    sprintf(buf, "%s : %s (#%d)", data[i].time, data[i].line, data[i].id);
	    Tcl_AppendElement(interp, buf);
	}
    } else if (args.raster_id > -1 && args.result_id == -1) {
	/* return list of all result names in a raster */
	/* find all those results in a particular raster */
	raster_result = raster_id_to_result(args.raster_id);
	
	info.job = SEQ_RESULT_INFO;
	info.op = WIN_NAME;
	info.result = NULL;
	
	Tcl_ResetResult(interp);
	for (i = 0; i < num_elements; i++) {
	    seq_result_notify(data[i].id, (seq_reg_data *)&info, 0);
	    
	    win_name = (char *)info.result;
	    if (strcmp(win_name, raster_result->raster_win) == 0) {
		sprintf(buf, "%s : %s (#%d)", data[i].time, data[i].line, data[i].id);
		Tcl_AppendElement(interp, buf);
	    }
	}
    } else {
	/* return specific result name */
	info.job = SEQ_RESULT_INFO;
	info.op = WIN_NAME;
	info.result = NULL;
	
	for (i = 0; i < num_elements; i++) {
	    if (data[i].id == args.result_id) {
		break;
	    }
	}
	Tcl_ResetResult(interp);
	seq_result_notify(data[i].id, (seq_reg_data *)&info, 0);
	sprintf(buf, "%s : %s (#%d)", data[i].time, data[i].line, data[i].id);
	Tcl_AppendElement(interp, buf);    
    }
    
    for (i = 0; i < num_results; i++) {
	xfree(data[i].line);
	xfree(data[i].time);
    }
    xfree(data);

    return TCL_OK;
}

/*
 * operations specific to a result
 */
int
tcl_seq_get_ops(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{    
    int l;
    seq_reg_get_ops ro;
    char *ops;
    get_ops_arg args;

    cli_args a[] = {
	{"-index", ARG_INT, 1, NULL, offsetof(get_ops_arg, id)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ro.job = SEQ_GET_OPS;
    ro.ops = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&ro, 0);

    if (!ro.ops) 
	return TCL_OK;

    ops = ro.ops;
    Tcl_ResetResult(interp);
    while(l = strlen(ops)) {
	Tcl_AppendElement(interp, ops);
	ops += l+1;
    }
    return TCL_OK;

}

/*
 * invoke operation on result
 */
int
tcl_seq_invoke_op(ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, 
		  char *argv[]) 
{    
    invoke_arg args;
    seq_reg_invoke_op inv;

    cli_args a[] = {
	{"-index", ARG_INT, 1, NULL, offsetof(invoke_arg, id)},
	{"-job",   ARG_INT, 1, NULL, offsetof(invoke_arg, option)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    inv.job = SEQ_INVOKE_OP;
    inv.op = args.option;
    seq_result_notify(args.id, (seq_reg_data *)&inv, 0);

    return TCL_OK;
}

/*
 * get the operations to be listed in the sequence manager list box
 */
int
tcl_seq_get_seq_ops(ClientData clientData, 
		    Tcl_Interp *interp, 
		    int argc, 
		    char *argv[]) 
{    
    seq_get_arg args;
    int l;
    char *ops;
    
    cli_args a[] = {
	{"-seq_num", ARG_INT, 1, NULL, offsetof(seq_get_arg, seq_num)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (GetSeqType(args.seq_num) == DNA) {
	ops = "Horizontal\0Vertical\0Set range\0Copy\0Complement\0Interconvert t and u\0Translate\0Scramble\0Type\0Rotate\0SEPARATOR\0Save\0Delete\0";
    } else {
	ops = "Horizontal\0Vertical\0Set range\0Copy\0PLACEHOLDER\0PLACEHOLDER\0PLACEHOLDER\0Scramble\0PLACEHOLDER\0Rotate\0SEPARATOR\0Save\0Delete\0";
    }
    Tcl_ResetResult(interp);
    while(l = strlen(ops)) {
	Tcl_AppendElement(interp, ops);
	ops += l+1;
    }
    return TCL_OK;

}

/*
 * operations invoked via sequence manager list box
 */
int
tcl_seq_invoke_seq_op(ClientData clientData, 
		      Tcl_Interp *interp, 
		      int argc, 
		      char *argv[]) 
{
    seq_invoke_arg args;
    char cmd[1024];

    cli_args a[] = {
	{"-seq_num", ARG_INT, 1, NULL, offsetof(seq_invoke_arg, seq_num)},
	{"-job",     ARG_INT, 1, NULL, offsetof(seq_invoke_arg, option)},
	{"-data",    ARG_INT, 1, "0",  offsetof(seq_invoke_arg, data)},
	{NULL,       0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    switch(args.option) {
    case 0: /* horizontal */
	Set_Active_Seq(args.seq_num, HORIZONTAL);
	Tcl_VarEval(interp, "sequence_list_update", NULL);
	break;
    case 1: /* vertical */
	Set_Active_Seq(args.seq_num, VERTICAL);
	Tcl_VarEval(interp, "sequence_list_update", NULL);
	break;
    case 2: /* set range */
	sprintf(cmd, "set_range_d %d", GetSeqId(args.seq_num));
	if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	    printf("SeqInvokeSeqOp: %s\n", Tcl_GetStringResult(interp));
	}
	break;
    case 3: /* copy range */
	sprintf(cmd, "copy_range_d %d", GetSeqId(args.seq_num));
	if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	    printf("SeqInvokeSeqOp: %s\n", Tcl_GetStringResult(interp));
	}
	break;
    case 4: /* complement */
	vfuncheader("complement sequence");
	if (GetSeqType(args.seq_num) == PROTEIN) {
	    verror(ERR_WARN, "Complement sequence", "Unable to complement a protein sequence");
	    break;
	}
	if (0 == ComplementSeq(interp, args.seq_num))
	    Tcl_Eval(interp, "sequence_list_update");
	
	break;
    case 5: /* Interconvert t and u */
        vfuncheader("interconvert t and u");
	if (GetSeqType(args.seq_num) == PROTEIN) {
	    verror(ERR_WARN, "Interconvert sequence", "Unable to interconvert t and u for a protein sequence");
	    break;
	}
	if (0 == RnaSeq(interp, args.seq_num))
	    Tcl_Eval(interp, "sequence_list_update");
	break;

    case 6: /* translate */
	{
	    sprintf(cmd, "translate_d %d", GetSeqId(args.seq_num));
	    if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
		verror(ERR_WARN, "SeqInvokeSeqOp", "%s\n", Tcl_GetStringResult(interp));
	    }
	}
	break;
    case 7: /* scramble */
        vfuncheader("scramble sequence");
	if (ScrambleSeq(interp, args.seq_num) == 0)
	    Tcl_Eval(interp, "sequence_list_update");
	break;
    case 8: /* type: change sequence type ie linear or circular */
	{
	    seq_reg_sequence_type info;
	    int *data;
	    data = &args.data;

	    SetSeqStructure(args.seq_num, args.data);
	    info.job = SEQ_SEQUENCE_TYPE;
	    info.data = (void *)data;
	    seq_notify(args.seq_num, (seq_reg_data *)&info);
	    Tcl_Eval(interp, "sequence_list_update");
	    break;
	}
    case 9: /* rotate */
	sprintf(cmd, "rotate_d %d", GetSeqId(args.seq_num));
	if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	    verror(ERR_WARN, "SeqInvokeSeqOp", "%s\n", Tcl_GetStringResult(interp));
	}
	break;
    case 10: {/* save */
	sprintf(cmd, "file_save_d %d", GetSeqId(args.seq_num));
	if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
	    verror(ERR_WARN, "SeqInvokeSeqOp", "%s\n", Tcl_GetStringResult(interp));
	}
	break;
    }
    case 11: /* delete */
	{
	    char cmd[100];
	    sprintf(cmd, "seq_shutdown %d\n", GetSeqId(args.seq_num));
	    DeleteSequence(interp, args.seq_num);
	    Tcl_Eval(interp, cmd);
	    Tcl_VarEval(interp, "sequence_list_update", NULL);
	}
	break;
    }
    return TCL_OK;
}

/*
 * hides or deletes all results on a single raster widget
 */
int
SeqResultUpdate(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{
    update_arg args;
    seq_reg_info info;

    cli_args a[] = {
	{"-index", ARG_INT, 1, "-1", offsetof(update_arg, id)},
	{"-job",   ARG_STR, 1, NULL, offsetof(update_arg, option)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (strcmp(args.option, "HIDE") == 0) {
	info.job = SEQ_HIDE;
    } else if (strcmp(args.option, "REVEAL") == 0) {
	info.job = SEQ_REVEAL;
    } else if (strcmp(args.option, "DELETE") == 0) {
	info.job = SEQ_DELETE;
    } else if (strcmp(args.option, "QUIT") == 0) {
	info.job = SEQ_QUIT;
    } else {
	verror(ERR_FATAL, "seq_result_notify_all", "invalid command");
	return TCL_OK;
    }

    if (args.id == -1) {
	seq_result_notify_all((seq_reg_data *)&info);
    } else {
	seq_result_notify(args.id, (seq_reg_data *)&info, 1);
    }
    return TCL_OK;
}

/*
 * gets the name associated with a result on the raster key
 */
int
SeqResultKeyName(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{
    update_arg args;
    seq_reg_key_name info;
    static char buf[80];

    cli_args a[] = {
	{"-index", ARG_INT, 1, NULL, offsetof(update_arg, id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_KEY_NAME;
    info.line = buf;

    seq_result_notify(args.id, (seq_reg_data *)&info, 0);

    vTcl_SetResult(interp, "%s", info.line);
    return TCL_OK;
}

int SeqGetBrief(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{
    update_arg args;
    seq_reg_brief info;
    static char buf[1024];

    cli_args a[] = {
	{"-index", ARG_INT, 1, NULL, offsetof(update_arg, id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_GET_BRIEF;
    info.line = buf;

    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    vTcl_SetResult(interp, "%s", info.line);
    return TCL_OK;
   
}

int
tcl_seq_result_info(ClientData clientData, 
		    Tcl_Interp *interp, 
		    int argc, 
		    char *argv[]) 
{

    result_info_arg args;
    out_raster *output;
    seq_reg_info info;
    d_point *pt;
    Tk_Raster *raster;
    Tcl_CmdInfo cmd_info;
    seq_result *result;
    int seq_num;

    cli_args a[] = {
	{"-index",     ARG_INT, 1, NULL, offsetof(result_info_arg, id)},
	{"-option",    ARG_STR, 1, NULL, offsetof(result_info_arg, option)},
	{"-direction", ARG_INT, 1, "0", offsetof(result_info_arg, direction)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;

    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    result = (seq_result *)info.result;

    if (result == NULL) {
	printf("result is null\n");
	return TCL_OK;
    }

    seq_num = GetSeqNum(result->seq_id[args.direction]);

    info.op = OUTPUT;
    info.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);

    if (info.result == NULL) {
	verror(ERR_WARN, "SeqResultInfo", "Result %d no longer exists", args.id);
	return TCL_OK;
    } else {
	output = (out_raster *)info.result;
    }

    info.op = WIN_SIZE;
    info.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info, 0);
    if (info.result == NULL) {
	verror(ERR_WARN, "SeqResultInfo", "Result %d no longer exists", args.id);
	return TCL_OK;
    } else {
	pt = (d_point *)info.result;
    }

    if (strcmp(args.option, "length") == 0) {
	vTcl_SetResult(interp, "%d", GetSeqLength(seq_num));
    } else if (strcmp(args.option, "type") == 0) {
	vTcl_SetResult(interp, "%d", seq_get_type(args.id));
    } else if (strcmp(args.option, "name") == 0) {
	vTcl_SetResult(interp, "%s", GetSeqName(seq_num));
    } else if (strcmp(args.option, "basename") == 0) {
	vTcl_SetResult(interp, "%s", GetSeqBaseName(seq_num));
    } else if (strcmp(args.option, "colour") == 0) {
	Tcl_GetCommandInfo(interp, output->raster_win, &cmd_info);
	raster = (Tk_Raster*)cmd_info.clientData;
	vTcl_SetResult(interp, "%s", GetRasterColour(interp, raster, 
						     output->env_index));
    } else if (strcmp(args.option, "raster") == 0) {
	vTcl_SetResult(interp, "%s", output->raster_win);
    } else if (strcmp(args.option, "win_size") == 0) {
	vTcl_SetResult(interp, "%d %f", pt->x, pt->y);
    } else {
	verror(ERR_WARN, "seq_result_info", "unknown option: %s\n", args.option);
    }
    return TCL_OK;
}

int CreateCursor(ClientData clientData, 
		 Tcl_Interp *interp, 
		 int argc, 
		 char *argv[]) 
{
    create_cursor_arg args;
    cursor_t *cp;

    cli_args a[] = {
	{"-seq_num",    ARG_INT, 1, NULL, offsetof(create_cursor_arg, seq_num)},
	{"-line_width", ARG_INT, 1, "2", offsetof(create_cursor_arg, line_width)},
	{"-direction", ARG_INT, 1, "0", offsetof(create_cursor_arg, direction)},
	{"-private", ARG_INT, 1, "0", offsetof(create_cursor_arg, private)},
	{NULL,       0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (cp = create_cursor(args.seq_num, args.private, NULL,
				    args.line_width, 1, args.direction))) {
	Tcl_SetResult(interp, "-1", TCL_STATIC);
	return TCL_OK;
    }
    
    vTcl_SetResult(interp, "%d", cp->id);
    return TCL_OK;
}


int DeleteCursor(ClientData clientData, 
		 Tcl_Interp *interp, 
		 int argc, 
		 char *argv[]) 
{
    delete_cursor_arg args;

    cli_args a[] = {
	{"-seq_num", ARG_INT, 1, NULL, offsetof(delete_cursor_arg, seq_num)},
	{"-id",      ARG_INT, 1, NULL, offsetof(delete_cursor_arg, id)},
	{"-private", ARG_INT, 1, "0",  offsetof(delete_cursor_arg, private)},
	{NULL,       0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    delete_cursor(args.seq_num, args.id, args.private);
    return TCL_OK;
}

int CursorNotify(ClientData clientData, 
		 Tcl_Interp *interp, 
		 int argc, 
		 char *argv[]) 
{
    cursor_notify_arg args;
    seq_cursor_notify cn;

    cli_args a[] = {
	{"-seq_num", ARG_INT, 1, NULL, offsetof(cursor_notify_arg, seq_num)},
	{"-id",      ARG_INT, 1, NULL, offsetof(cursor_notify_arg, id)},
	{"-pos",     ARG_INT, 1, NULL, offsetof(cursor_notify_arg, pos)},
	{"-direction", ARG_INT, 1, "0", offsetof(cursor_notify_arg, direction)},
	{NULL,       0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (cn.cursor = find_cursor(&args.seq_num, args.id, args.direction))) {
#ifdef DEBUG
	printf("CursorNotify: unable to find cursor for seq_num %d cursor_id %d\n",
	       args.seq_num, args.id);
#endif
	return TCL_OK;
    }
    cn.cursor->abspos = args.pos;
    cn.cursor->job = CURSOR_MOVE;
    cn.cursor->sent_by = -1; /* FIXME: Don't know correct value. */
    cn.job = SEQ_CURSOR_NOTIFY;

    seq_notify(args.seq_num, (seq_reg_data *)&cn);

    return TCL_OK;
}

int CursorRef(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    cursor_ref_arg args;
    cursor_t *cursor;
    seq_cursor_notify cn;

    cli_args a[] = {
	{"-seq_num",   ARG_INT, 1, NULL, offsetof(cursor_ref_arg, seq_num)},
	{"-id",        ARG_INT, 1, NULL, offsetof(cursor_ref_arg, id)},
	{"-ref",       ARG_INT, 1, NULL, offsetof(cursor_ref_arg, ref)},
	{"-direction", ARG_INT, 1, "0", offsetof(cursor_ref_arg, direction)},
	{NULL,       0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (NULL == (cursor = find_cursor(&args.seq_num, args.id, args.direction))) {
#ifdef DEBUG
	printf("CursorRef: unable to find cursor for seq_num %d cursor_id %d\n",
	       args.seq_num, args.id);
#endif
	return TCL_OK;
    }

    cursor->refs += args.ref;

    cn.cursor = cursor;
    cn.cursor->job = CURSOR_MOVE;
    cn.job = SEQ_CURSOR_NOTIFY;

    seq_notify(args.seq_num, (seq_reg_data *)&cn);
    return TCL_OK;
}

int QueryCursor(ClientData clientData, Tcl_Interp *interp,
		int argc, char **argv)
{
    qc_arg args;

    cli_args a[] = {
	{"-cursorid", ARG_INT, 1, NULL, offsetof(qc_arg, cursorid)},
	{"-seq_num",  ARG_INT, 1, "0",  offsetof(qc_arg, seq_num)},
        {NULL,        0,       0, NULL, 0}
    };
    cursor_t *gc;

    if (-1 == parse_args(a, &args, argc, argv))
        return TCL_ERROR;

    if (NULL == (gc = find_cursor(&args.seq_num, args.cursorid, -1)))
	return TCL_OK;

    vTcl_SetResult(interp,
		   "{id %d} {refs %d} {private %d} {abspos %d}",
		   gc->id, gc->refs, gc->private, gc->abspos);
    return TCL_OK;
}

/*
 * read a file of dna or protein sequences and extract their identifiers
 */
int
GetArchiveList(ClientData clientData, 
	       Tcl_Interp *interp, 
	       int argc, 
	       char *argv[]) 
{
    archive_list_arg args;
    char **identifier;
    int num_identifiers;
    int i;
    Array id_array;
    int num_elements = 100;
    char word[20];

    cli_args a[] = {
	{"-file",  ARG_STR, 1, NULL, offsetof(archive_list_arg, file)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (strcmp(args.file, "") == 0) 
	return TCL_OK;

    if (NULL == (id_array = (ArrayCreate(sizeof(word), 
					 num_elements))))
	 return TCL_OK;
    
    if (get_identifiers(args.file, &identifier, &num_identifiers) != 0) {
	verror(ERR_WARN, "reading archive list", "error reading %s", 
	       args.file);
	return TCL_OK;
    }

    Tcl_ResetResult(interp);
    for (i = 0; i < num_identifiers; i++) {
	Tcl_AppendElement(interp, identifier[i]);
    }

    for (i = 0; i < num_identifiers; i++) {
	xfree(identifier[i]);
    }
    xfree(identifier);

    return TCL_OK;
}

/*
 * seq_info seq_id type
 * seq_info seq_id length
 * seq_info seq_id number
 * seq_info seq_id name
 * seq_info seq_id sequence
 * seq_info seq_id is_sub_seq
 * seq_info seq_id structure
 */
int
tcl_seq_info(ClientData clientData, 
	     Tcl_Interp *interp, 
	     int argc, 
	     char *argv[]) 
{
    int seq_num; 
    if (argc < 3) {
	goto seq_error;
    }

    seq_num = GetSeqNum(atoi(argv[1]));

    if (seq_num == -1) {
      verror(ERR_WARN, "tcl seq info", "Invalid sequence number");
      vTcl_SetResult(interp, "%d", -1);
      return TCL_OK;
    }

    if (strcmp(argv[2], "type") == 0 ) {
	vTcl_SetResult(interp, "%d", GetSeqType(seq_num));
    } else if (strcmp(argv[2], "structure") == 0 ) {
	vTcl_SetResult(interp, "%d", GetSeqStructure(seq_num));
/* added for reading feature tables */
    } else if (strcmp(argv[2], "key_index_cds") == 0 ) {
      vTcl_SetResult(interp, "%s", GetSeqCdsExpr(seq_num,atoi(argv[3])));
    } else if (strcmp(argv[2], "numbercds") == 0 ) {
	vTcl_SetResult(interp, "%d", GetSeqNumberCds(seq_num));
/* added for reading feature tables */
    } else if (strcmp(argv[2], "start") == 0 ) {
	vTcl_SetResult(interp, "%d", GetSubSeqStart(seq_num));
    } else if (strcmp(argv[2], "end") == 0 ) {
	vTcl_SetResult(interp, "%d", GetSubSeqEnd(seq_num));
    } else if (strcmp(argv[2], "length") == 0 ) {
	vTcl_SetResult(interp, "%d", GetSeqLength(seq_num));
    } else if (strcmp(argv[2], "mass") == 0 ) {
	vTcl_SetResult(interp, "%f", get_seq_mass(seq_num));
    } else if (strcmp(argv[2], "number") == 0 ) {
	vTcl_SetResult(interp, "%d", GetSeqNum(seq_num));
    } else if (strcmp(argv[2], "name") == 0 ) {
	vTcl_SetResult(interp, "%s", GetSeqName(seq_num));
    } else if (strcmp(argv[2], "library") == 0 ) {
	vTcl_SetResult(interp, "%s", GetSeqLibraryName(seq_num));
    } else if (strcmp(argv[2], "sequence") == 0 ) {
	vTcl_SetResult(interp, "%s", GetSeqSequence(seq_num));
    } else if (strcmp(argv[2], "is_sub_seq") == 0 ) {
	if (strcmp(GetParentalSeqName(seq_num),
		   GetSeqName(seq_num)) == 0) {
	    vTcl_SetResult(interp, "%d", 0);
	} else {
	    vTcl_SetResult(interp, "%d", 1);
	}
    } else {
	goto seq_error;
    }
    return TCL_OK;
 
seq_error:
    Tcl_SetResult(interp,
		  "wrong # args: should be \"tcl_seq_info seq_id command\"\n",
		  TCL_STATIC);
    return TCL_ERROR;
}

typedef struct {
    int seq_id;
    int structure;
} seq_structure_arg;

int
tcl_set_seq_structure(ClientData clientData, 
		      Tcl_Interp *interp, 
		      int argc, 
		      char *argv[]) 
{
    seq_structure_arg args;
    char str[20];

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, "-1", offsetof(seq_structure_arg, seq_id)},
	{"-structure", ARG_INT, 1, "-1", offsetof(seq_structure_arg, structure)},
	{NULL,      0,       0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    SetSeqStructure(GetSeqNum(args.seq_id), args.structure);

    if (args.structure == 0) {
	strcpy(str, "linear");
    } else {
	strcpy(str, "circular");
    }

    vfuncheader("Sequence structure");
    vmessage("Sequence %s is %s\n", GetSeqName(GetSeqNum(args.seq_id)), str);
    return TCL_OK;
}

/*
 * return seq_id given a seq_num
 */
int
tcl_GetSeqId(ClientData clientData, 
	     Tcl_Interp *interp, 
	     int argc, 
	     char *argv[]) 
{
    vTcl_SetResult(interp, "%d", GetSeqId(atoi(argv[1])));
    return TCL_OK;
}

/*
 * return seq_num given a seq_id
 */
int
tcl_GetSeqNum(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    vTcl_SetResult(interp, "%d", GetSeqNum(atoi(argv[1])));
    return TCL_OK;
}

/*
 * find the seqid from the name 
 * return -1 if none exist
 */
int
NameToSeqId(ClientData clientData, 
	    Tcl_Interp *interp, 
	    int argc, 
	    char *argv[]) 
{

    vTcl_SetResult(interp, "%d", GetSeqIdFromName(argv[1]));

    return TCL_OK;
}

int
GetActiveSeqId(ClientData clientData, 
	       Tcl_Interp *interp, 
	       int argc, 
	       char *argv[]) 
{
    int direction = 0; /* default direction is horizontal */;

    if (argc > 1) {
	direction = atoi(argv[1]);
    }

    vTcl_SetResult(interp, "%d", GetSeqId(GetActiveSeqNumber(direction)));
    return TCL_OK;
}

int
GetActiveSeqName(ClientData clientData, 
		 Tcl_Interp *interp, 
		 int argc, 
		 char *argv[]) 
{
    int direction = -1;

    if (argc > 1) {
	direction = atoi(argv[1]);
    }
    if (GetActiveSeqNumber(direction) > -1) {
	vTcl_SetResult(interp, "%s", 
		       GetSeqName(GetActiveSeqNumber(direction)));
    } else {
	vTcl_SetResult(interp, "");
    }
    return TCL_OK;
}

int
tcl_NumSequences(ClientData clientData, 
		 Tcl_Interp *interp, 
		 int argc, 
		 char *argv[])
{
    vTcl_SetResult(interp, "%d", NumSequences());	    
    return TCL_OK;
}


int
tcl_s_length(ClientData clientData, 
	     Tcl_Interp *interp, 
	     int argc, 
	     char *argv[]) 
{
    seq_id_arg args;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, "-1", offsetof(seq_id_arg, seq_id)},
	{NULL,      0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (args.seq_id == -1) {
	vTcl_SetResult(interp, "%d", GetSeqLength(GetActiveSeqNumber(0)));
    } else {
	vTcl_SetResult(interp, "%d", GetSeqLength(GetSeqNum(args.seq_id)));
    }
    return TCL_OK;
}

int SeqFileSave(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{
    file_save_arg args;
    FILE *fp;
    char *seq;
    char *identifier;
    int line_len = 60;
    int seqnum;
    int start, end, i;
    Featcds **keys;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(file_save_arg, seq_id)},
	{"-start",  ARG_INT, 1, "0",  offsetof(file_save_arg, start)},
	{"-end",    ARG_INT, 1, "0",  offsetof(file_save_arg, end)},
	{"-format", ARG_INT, 1, "0",  offsetof(file_save_arg, format)},
	{"-file",   ARG_STR, 1, NULL, offsetof(file_save_arg, file)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    vfuncheader("save sequence");

    if (NULL == (fp = fopen(args.file, "w"))) {
	verror(ERR_WARN, "save sequence", "Unable to save sequence");
	return TCL_OK;
    }

    seqnum = GetSeqNum(args.seq_id);
    seq = GetSeqSequence(seqnum);

    /* identifier = GetSeqIdentifier(seqnum); */

    identifier = GetSeqName(seqnum);
    keys = GetSeqKeyIndex( seqnum );

    if ((start = args.start) == 0) {
	start = 1;
    }
    if ((end = args.end) == 0) {
	end = strlen(seq);
    }
 
    if(args.format == 2){

	 /* write out sequence in EMBL format */

	seq_file_save_ft(fp, keys, start, end, seq, identifier);	 
    } else {

	/* write out sequence in FASTA format */
	fprintf(fp,">%s\n", identifier);
	fputc(seq[start-1], fp);
	for (i = start; i < end; i++) {
	    if ((i-start+1) % line_len == 0) {
	    fputc('\n', fp);
	    }
	    fputc(seq[i], fp);
	}
	fprintf(fp,"\n");
    }
    fclose(fp);
    return TCL_OK;
}

static int seq_file_save_ft(FILE *pw, Featcds **key_index, int start, int end, 
			    char *seq, char *identifier){

    int i, j, k, h;
    int l = 0;
    int line_len = 60;
    int line_seg = 10;
    int num_comma = 0;
    int cds_len, qua_len;
    char *cdsexpr, *quaexpr;
    char qua;

    fprintf(pw,"ID   %s\n", identifier);

    if(key_index != NULL) {
    if( strlen(seq) == end - start + 1) {
      for(i = 0; i < number_keys; i++){
	for(k = 1; k <= key_index[i]->id; k++){

	  cds_len = strlen(key_index[i][k].cdsexpr);
	  cdsexpr = key_index[i][k].cdsexpr;

	  if(cdsexpr != NULL && cds_len < 60){
	    fprintf(pw,"FT   %-16s%s", feat_key[i], cdsexpr);
	  } else if(cdsexpr != NULL) {
	    fprintf(pw,"FT   %-16s", feat_key[i]);	       
	    for(j = 0; j < cds_len; j++){
	      fputc(cdsexpr[j], pw);
	      if(cdsexpr[j] == ',') num_comma++;
	      if( j > 1 && cdsexpr[j] == ','&& num_comma == 2 ) {
		fprintf(pw, "\nFT                   ");
		num_comma = 0;
	      }  
	    }
	  }
 
	  for(j = 0; j < number_quas; j++){
	    qua_len = strlen(key_index[i][k].qualifier[j]);
	    quaexpr = key_index[i][k].qualifier[j];
	    l = 0;
	    if (qua_len > 1 ){
	      fprintf(pw, "\nFT                   ");
	      for(h = 0; h < qua_len; h++){
		l++;	  
		if( l > 1 && l % 60 == 0 ||  quaexpr[h] == '?' ){
		  if( l > 0 && quaexpr[h] == '?'){
		    l = 0;
		    h++;
		  }
		  fprintf(pw, "\nFT                   ");
		}
		fputc(quaexpr[h], pw);
		qua = quaexpr[h];
	      }
	    }
	  }
	  fprintf(pw, "\n");
	}
      }
    } 
    }    
    fprintf(pw,"SQ   \n");
    fprintf(pw, "    ");
    j = 0;   
    for (i = start-1; i < end; i++) {
      if (i > start && (i-start+1) % line_len == 0) {
	  j = 0;
	  fprintf(pw, "%10d\n", (i-start+1 ));
	  fprintf(pw, "    ");
	}
        if ((i-start+1) % line_seg == 0) {
	  fputc(' ', pw);
	  j = j+1;
	}
	fputc(seq[i], pw);
	j = j+1;
    }
    for (i = 1; i <= 66-j; i++)
      fprintf(pw, "%c", ' ');
    fprintf(pw, "%10d\n", end - start + 1);
    fprintf(pw, "//\n"); 
    return 1;
}

int SeqFileDelete(ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, 
		  char *argv[]) 
{
    seq_id_arg args;
    int seq_num;
    char cmd[100];

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(seq_id_arg, seq_id)},
	{NULL,     0,       0, NULL, 0}
    };

    vfuncheader("delete sequence");
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    seq_num = GetSeqNum(args.seq_id);    
    DeleteSequence(interp, seq_num);
    sprintf(cmd, "seq_shutdown %d\n", args.seq_id);
    Tcl_Eval(interp, cmd);
    return TCL_OK;
}

int SeqComplement(ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, 
		  char *argv[]) 
{
    seq_id_arg args;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(seq_id_arg, seq_id)},
	{NULL,     0,       0, NULL, 0}
    };

    vfuncheader("complement sequence");
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ComplementSeq(interp, GetSeqNum(args.seq_id));
    return TCL_OK;
}

int SeqSetActiveSeq(ClientData clientData, 
		    Tcl_Interp *interp, 
		    int argc, 
		    char *argv[]) 
{
    active_seq_arg args;

    cli_args a[] = {
	{"-seq_id",    ARG_INT, 1, NULL, offsetof(active_seq_arg, seq_id)},
	{"-direction", ARG_INT, 1, "0", offsetof(active_seq_arg, direction)},
	{NULL,      0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    vfuncheader("set active sequence");
    Set_Active_Seq(GetSeqNum(args.seq_id), args.direction);

    return TCL_OK;
}

int SeqInterconvert(ClientData clientData, 
		    Tcl_Interp *interp, 
		    int argc, 
		    char *argv[]) 
{
    seq_id_arg args;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(seq_id_arg, seq_id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    vfuncheader("interconvert t and u");
    RnaSeq(interp, GetSeqNum(args.seq_id));

    return TCL_OK;
}

int SeqScramble(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{
    seq_id_arg args;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(seq_id_arg, seq_id)},
	{NULL,     0,       0, NULL, 0}
    };

    vfuncheader("scramble sequence");
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    ScrambleSeq(interp, GetSeqNum(args.seq_id));
    	
    return TCL_OK;
}

int SeqRotate(ClientData clientData, 
	      Tcl_Interp *interp, 
	      int argc, 
	      char *argv[]) 
{
    rotate_arg args;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(rotate_arg, seq_id)},
	{"-origin", ARG_INT, 1, NULL, offsetof(rotate_arg, origin)},
	{NULL,     0,       0, NULL, 0}
    };

    vfuncheader("rotate sequence");
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    RotateSeq(interp, GetSeqNum(args.seq_id), args.origin);
    return TCL_OK;
}

int 
SeqSetRange(ClientData clientData, 
	    Tcl_Interp *interp, 
	    int argc, 
	    char *argv[]) 
{
    set_range_arg args;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(set_range_arg, seq_id)},
	{"-start",  ARG_INT, 1, "1",  offsetof(set_range_arg, start)},
	{"-end",    ARG_INT, 1, "-1", offsetof(set_range_arg, end)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    vfuncheader("set range");
    
    /* if the end has not been defined, set it to be the sequence length */
    if (args.end == -1) {
	args.end = GetSeqLength(GetSeqNum(args.seq_id));
    }

    SetRange(interp, args.seq_id, args.start, args.end);
    return TCL_OK;
}

int 
SeqCopyRange(ClientData clientData, 
	     Tcl_Interp *interp, 
	     int argc, 
	     char *argv[]) 
{
    set_range_arg args;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(set_range_arg, seq_id)},
	{"-start",  ARG_INT, 1, "1",  offsetof(set_range_arg, start)},
	{"-end",    ARG_INT, 1, "-1", offsetof(set_range_arg, end)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    vfuncheader("set range");
    
    /* if the end has not been defined, set it to be the sequence length */
    if (args.end == -1) {
	args.end = GetSeqLength(GetSeqNum(args.seq_id));
    }

    CopyRange(interp, args.seq_id, args.start, args.end);
    return TCL_OK;
}

/*
 * translate sequence in sequence manager
 */
int 
SeqTranslateSequence(ClientData clientData, 
		     Tcl_Interp *interp, 
		     int argc, 
		     char *argv[]) 
{
    translate_seq_arg args;
    int num[3], seq_num;
    char buf[100];

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(translate_seq_arg, seq_id)},
	{"-f1",     ARG_INT, 1, "0",  offsetof(translate_seq_arg, f1)},
	{"-f2",     ARG_INT, 1, "0",  offsetof(translate_seq_arg, f2)},
	{"-f3",     ARG_INT, 1, "0",  offsetof(translate_seq_arg, f3)},
	{"-all",    ARG_INT, 1, "0",  offsetof(translate_seq_arg, all)},
	{"-start",  ARG_INT, 1, "0",  offsetof(translate_seq_arg, start)},
	{"-end",    ARG_INT, 1, "0",  offsetof(translate_seq_arg, end)},
	{NULL,      0,       0, NULL, 0}
    };

    vfuncheader("translate sequence");

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (args.start == 0) 
	args.start = 1;

    num[0] = num[1] = num[2] = -1;
    seq_num = GetSeqNum(args.seq_id);

    if (args.end == 0) 
	args.end = GetSeqLength(seq_num);

    if (args.f1) 
	num[0] = TranslateSeq(interp, seq_num, 0, args.start, args.end);
    if (args.f2) 
	num[1] = TranslateSeq(interp, seq_num, 1, args.start, args.end);
    if (args.f3) 
	num[2] = TranslateSeq(interp, seq_num, 2, args.start, args.end);
    if (args.all)
	TranslateTogether(interp, seq_num);

    /* return seq id */
    Tcl_ResetResult(interp);
    if (num[0] > -1) {
	sprintf(buf, "%d", GetSeqId(num[0]));
	Tcl_AppendElement(interp, buf);
    }
    if (num[1] > -1) {
	sprintf(buf, "%d", GetSeqId(num[1]));
	Tcl_AppendElement(interp, buf);
    }
    if (num[2] > -1) {
	sprintf(buf, "%d", GetSeqId(num[2]));
	Tcl_AppendElement(interp, buf);
    }

    return TCL_OK;
}


int
SeqSender(ClientData clientData, 
	  Tcl_Interp *interp, 
	  int argc, 
	  char *argv[]) 
{

    sender_arg args;
    cli_args a[] = {
	{"-rid",    ARG_STR, 1, NULL, offsetof(sender_arg, rid)},
	{"-seq_id", ARG_INT, 1, NULL, offsetof(sender_arg, seq_id)},
	{NULL,      0,  0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    vTcl_SetResult(interp, "%d", seq_sender(interp, args.rid, args.seq_id));
    return TCL_OK;
}

/*
 * extracts sequences from either a personal library or personal file or
 * sequence library.
 * personal library:- check if library == num_libs; set file, entry
 * read in sequence:- check if library == -1, sequence != NULL
 * library:- check if library > -1; set library, entry_mode, entry
 * sets these to be the default horizontal and vertical sequences for future
 * operations within sip.
 */
int
tcl_read_sequence(ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, 
		  char *argv[]) 
{    
#define MAX_SEQ 100000 /* not actually needed, get_seq reallocs if necessary */

#define SET_LIBRARY 0
#define SET_SINGLE_PERSONAL 1
#define SET_MULTI_PERSONAL 2
#define SET_SEQUENCE 3

    /* added for parsing feature tables */
    Featcds **key_index = NULL;
    char *identifier = NULL;
    int err = 0;
    /* added for parsing feature tables */

    set_seq_arg args;
    char *sequence = NULL;
    int seq_len, max_len;
    char *seq = NULL;
    char direction[11];
    static int initialised_seq_reg = 0;
    int seq_num;
    int num_seqs;
    int i, j;
    static int unique_name = 0;
    char *name = NULL;
    int num_entry = 0;
    char **entry;
    char buf[10];
    Tcl_DString id_list;
    int mode;
    int type;
    int err_code;


    cli_args a[] = {
	{"-library",    ARG_INT, 1, "-1", offsetof(set_seq_arg, library)},
	{"-entry_mode", ARG_INT, 1, "0",  offsetof(set_seq_arg, entry_mode)},
	{"-file",       ARG_STR, 1, "",   offsetof(set_seq_arg, file)},
	{"-entry",      ARG_STR, 1, "",   offsetof(set_seq_arg, entry)},
	{"-sequence",   ARG_STR, 1, "",   offsetof(set_seq_arg, sequence)},
	{"-direction",  ARG_INT, 1, "0",  offsetof(set_seq_arg, direction)},
	{"-name",	ARG_STR, 1, "",   offsetof(set_seq_arg, name)},
	{NULL,          0,       0, NULL, 0}
    };

   
  /* added for parsing feature tables */
    /*  if (NULL == (key_index = (Featcds **)xmalloc(number_keys*sizeof(Featcds *))))
      return -1;
  for (i=0; i < number_keys; i++){
   if (NULL == (key_index[i] = (Featcds *)xmalloc(sizeof(Featcds ))))  
     return -1;  
   key_index[i]->id=0;
   } */
  /* added for parsing feature tables */
    
    vfuncheader("read sequence");
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

#ifdef DEBUG
    printf("args %d %d @%s@ @%s@ @%s@ %d\n", args.library, args.entry_mode,
	   args.file, args.entry, args.sequence, args.direction);
#endif

    /* initialise sequence registration if first time thro' */
    
    if (!initialised_seq_reg) {
	seq_register_init(interp);
	initialised_seq_reg = 1;
    }
   
    /* set up a string in case of error messages */
    if (args.direction == HORIZONTAL) {
	strcpy(direction, "horizontal");
    } else if (args.direction == VERTICAL) {
	strcpy(direction, "vertical");
    }

    /* args.entry can contain multiple entries */
    if (Tcl_SplitList(interp, args.entry, &num_entry, &entry) != TCL_OK) {
	return TCL_ERROR;
    }

    if (args.library > -1) {
	vmessage("Reading in library files \n");
	mode = SET_LIBRARY;
#ifdef DEBUG
	printf("SET_LIBRARY\n");
#endif
	if (strcmp(args.entry, "") == 0) {
	    verror(ERR_WARN, "read_sequence", "need to specify -entry entryname(s)\n");
	    vTcl_SetResult(interp, "%d", -1);
	    return TCL_OK;
	}
    } else if (strcmp(args.sequence, "") != 0) {
	vmessage("Reading in sequence\n");
	mode = SET_SEQUENCE;
#ifdef DEBUG
	printf("SET_SEQUENCE\n");
#endif
	if (strcmp(args.entry, "") == 0) {
	    verror(ERR_WARN, "read_sequence", "need to specify -entry entryname(s)\n");
	    vTcl_SetResult(interp, "%d", -1);	    
	    return TCL_OK;
	}
    } else if (strcmp(args.entry, "") == 0) {
	vmessage("Reading in single personal file\n");
	mode = SET_SINGLE_PERSONAL;
#ifdef DEBUG
	printf("SET_SINGLE_PERSONAL\n");
#endif
	num_entry = 1;
	entry[0] = '\0';
	if (strcmp(args.file, "") == 0) {
	    verror(ERR_WARN, "read_sequence", "need to specify -file filename\n");
	    vTcl_SetResult(interp, "%d", -1);	    
	    return TCL_OK;
	}
    } else if (strcmp(args.entry, "") != 0) {
	vmessage("Reading in multiple personal files \n");
	mode = SET_MULTI_PERSONAL;
#ifdef DEBUG
	printf("SET_MULTI_PERSONAL\n");
#endif
	if (strcmp(args.file, "") == 0) {
	    verror(ERR_WARN, "read_sequence", "need to specify -file filename\n");
	    vTcl_SetResult(interp, "%d", -1);	    
	    return TCL_OK;
	}
    }

    /* 
     * loop round for each entryname selected. Maybe not the most efficient
     * method, but simple for now
     */
    Tcl_DStringInit(&id_list);
    for (j = 0; j < num_entry; j++) {

	if (NULL == (key_index = (Featcds **)xmalloc(number_keys*sizeof(Featcds *))))
	    return -1;
	for (i=0; i < number_keys; i++){
	    if (NULL == (key_index[i] = (Featcds *)xmalloc(sizeof(Featcds ))))  
		return -1;  
	    key_index[i]->id=0;
	} 

	if (mode == SET_SEQUENCE) {
	    /* read in personal sequence */
	    sequence = strdup(args.sequence);
	    
	} else if (mode == SET_SINGLE_PERSONAL || mode == 
		   SET_MULTI_PERSONAL) {
	    /* get personal sequence from file */
	    
	    max_len = MAX_SEQ;

  /* modified for parsing feature tables */
	    err_code = get_seq_ft( key_index ,&seq, max_len, &seq_len, args.file, entry[j], &identifier, &err);

	    if( err_code || !seq_len )
	    {
		switch( err_code )
		    {
		    case 1: verror(ERR_WARN, "read_sequence", "failed to open file %s", args.name);  break;
		    case 2: verror(ERR_WARN, "read_sequence", "file %s is not dna or protein", args.name);  break;
		    case 3: verror(ERR_WARN, "read_sequence", "file %s is of unrecognised format", args.name); break;
		    case 4: verror(ERR_WARN, "read_sequence", "failed to seek to start of file %s", args.name); break;
		    case 5: /* Memory allocation error, xmalloc prints message */
		    default:  break;
		}
		Tcl_ResetResult(interp);
		vTcl_SetResult(interp, "%d", -1);
		return TCL_OK;
	    }

	    if (NULL == (sequence = (char *)xmalloc((seq_len+1) * sizeof(char)))) {
		vTcl_SetResult(interp, "%d", -1);	    
		return TCL_OK;
	    }
	    strncpy(sequence, seq, seq_len);
	    sequence[seq_len] = '\0';
	    xfree(seq);
	    seq = NULL;
	    /* printf("%s \n %s \n %d\n", seq, sequence, seq_len);  */

#ifdef USE_SEQLIB
	} else if (mode == SET_LIBRARY) {
	    /* get sequence from library */
	    if (-1 == ExtractSequence(interp, args.library, args.entry_mode, 
				      entry[j], &sequence)) {
		verror(ERR_WARN, "Extract Sequence", "%s sequence %s not found", 
		       direction, args.entry); 
		Tcl_ResetResult(interp);
		vTcl_SetResult(interp, "%d", -1);
		return TCL_OK;
	    }
#endif
	} else {
	    verror(ERR_WARN, "read_sequence", "invalid arguments"); 
	    Tcl_ResetResult(interp);
	    vTcl_SetResult(interp, "%d", -1);	    
	    return TCL_OK;
	}

	if (0 == (type = get_seq_type(sequence, strlen(sequence)))) {
	    verror(ERR_WARN, "read_sequence", "sequence is neither dna nor protein");
	    vTcl_SetResult(interp, "%d", -1);
	    return TCL_OK;
	}
       
	if (err == 1 && type == 1) {
	       verror(ERR_WARN, "Load Sequence",
		      "Problem parsing feature table\n");
	}


	/* 
	 * FIXME - this needs some more thought - perhaps shouldn't merely 
	 * produce a unique name but update the sequence properly
	 */

	/*
	 * need to check if entry_name already exists in the database. If it does
	 * then must create a unique name for that sequence.
	 */
	/* if name hasn't been set yet, ie first time through */
	if (!name) {
	    if (mode == SET_SINGLE_PERSONAL) {
		name = strdup(*args.name ? args.name : args.file);
	    } else {
		name = strdup(*args.name ? args.name : entry[j]);
	    }
	}

	/* set up identifier */
	
	if (!identifier) {
	    identifier = strdup(name);
	}

	num_seqs = NumSequences();
	for (i = 0; i < num_seqs; i++) {
	    if (strcmp(identifier, GetSeqName(i)) == 0) {

		if (NULL == (identifier = (char *)realloc(identifier, 
						    (strlen(identifier)+10)))){ 
		    verror(ERR_WARN, "read_sequence", "Unable to reallocate "
			   "memory\n");
		    vTcl_SetResult(interp, "%d", -1);	    
		    return TCL_OK;
		}
		sprintf(identifier, "%s#%d", identifier, unique_name++);
	    }
	}

	/* saved sequences */
	/* default structure be LINEAR and default type is 0 (unknown) */
	if (AddSequence(interp, args.direction, args.library, identifier, 
			   sequence,  LINEAR, 0, key_index, identifier) == -1) {
	    verror(ERR_WARN, "adding sequence", 
		   "failure to add sequence to registration scheme");
	    vTcl_SetResult(interp, "-1");
	    return TCL_OK;
	}
	
	
	seq_num = GetActiveSeqNumber(args.direction);

	sequence_info(GetSeqName(seq_num), GetSeqSequence(seq_num),
		      1, GetSeqLength(seq_num), 0, GetSeqType(seq_num));

	sprintf(buf, "%d", GetSeqId(seq_num));
	Tcl_DStringAppendElement(&id_list, buf);
	if (seq)
	    xfree(seq);
	if (identifier)
	    xfree(identifier);
	free(name);
	name = NULL;
    }

    vTcl_SetResult(interp, "%s", Tcl_DStringValue(&id_list));
    Tcl_DStringFree(&id_list);
    Tcl_Free((char *)entry);
    return TCL_OK;
}

/*
 * current sequences to display in the sequence manager list box
 */
int
tcl_sequence_names(ClientData clientData, 
		   Tcl_Interp *interp, 
		   int argc, 
		   char *argv[]) 
{
    int num_seqs;
    char *name;
    int i;
    int direction, type, seq_structure;
    char dir, seqtype, seqstruct;
    Tcl_DString ds;

    Tcl_DStringInit( &ds );
    num_seqs = NumSequences();
    Tcl_ResetResult(interp);
    for (i = 0; i < num_seqs; i++) {
	name = GetSeqName(i);
	direction = GetSeqDirection(i);
	if (direction == HORIZONTAL)
	    dir = 'H';
	else if (direction == VERTICAL)
	    dir = 'V';
	else 
	    dir = ' ';
	type = GetSeqType(i);
	if (type == DNA) 
	    seqtype = 'D';
	else if (type == PROTEIN) 
		seqtype = 'P';
	else 
	    seqtype = ' ';

	seq_structure = GetSeqStructure(i);
	if (seq_structure == LINEAR)
	    seqstruct = 'L';
	else if (seq_structure == CIRCULAR)
	    seqstruct = 'C';
	else 
	    seqstruct = ' ';

	Tcl_DStringStartSublist(&ds);
	vTcl_DStringAppendElement( &ds, "%c", dir );
	vTcl_DStringAppendElement( &ds, "%s", name );
	vTcl_DStringAppendElement( &ds, "%d..%d", GetSubSeqStart(i), GetSubSeqEnd(i) );
	vTcl_DStringAppendElement( &ds, "%d", GetSubSeqLength(i) );
	vTcl_DStringAppendElement( &ds, "%c", seqtype );
	vTcl_DStringAppendElement( &ds, "%c", seqstruct );
	Tcl_DStringEndSublist(&ds);
    }
    Tcl_DStringResult(interp, &ds);
    return TCL_OK;
}

/*
 * return all seq_id's registered with a result
 */
int
tcl_get_result_seq_id(ClientData clientData, 
		      Tcl_Interp *interp, 
		      int argc, 
		      char *argv[]) 
{

    result_id_arg args;
    int i;
    char str[1024];
    RasterResult *result;
    seq_result *seq_result;
    int type;

    cli_args a[] = {
	{"-result_id", ARG_INT, 1, NULL, offsetof(result_id_arg, result_id)},
	{NULL,      0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    type = seq_get_type(args.result_id);
    Tcl_ResetResult(interp);
    if (type == SEQ_RASTER) {
	result = raster_id_to_result(args.result_id);
	for (i = 0; i < result->num_seq_id; i++) {
	    sprintf(str, "%d %d", result->seq[i].seq_id, result->seq[i].direction);
	    Tcl_AppendElement(interp, str); 
	}
    } else if (type == SEQ_PLOT_PERM || type == SEQ_PLOT_TEMP) {
	seq_result = seq_id_to_result(args.result_id);
	sprintf(str, "%d %d", seq_result->seq_id[HORIZONTAL], HORIZONTAL);
	Tcl_AppendElement(interp, str); 
	if (seq_result->seq_id[VERTICAL] != -1) {
	    sprintf(str, "%d %d", seq_result->seq_id[VERTICAL], VERTICAL);
	    Tcl_AppendElement(interp, str); 
	}
    }
    return TCL_OK;
}

/*
 * finds the nearest match from x, y coords from a result
 */
int 
tcl_nearest_match(ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, 
		  char *argv[]) 
{    
    nearest_match_arg args;
    d_point start;
    d_point nearest;
    seq_reg_info info1, info2;
    seq_result *result;
    out_raster *output;
    Tk_Raster *raster;
    Tcl_CmdInfo cmd_info;
    double wx0, wy0, wx1, wy1;    
    double sx0, sy0, sx1, sy1;

    cli_args a[] = {
	{"-x",  ARG_INT, 1, NULL, offsetof(nearest_match_arg, pt_x)},
	{"-y",  ARG_INT, 1, NULL, offsetof(nearest_match_arg, pt_y)},
	{"-result_index", ARG_INT, 1, NULL, offsetof(nearest_match_arg, id)},
	{"-match", ARG_INT, 1, NULL, offsetof(nearest_match_arg, match)},
	{NULL,     0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

   /* printf("left1 %d left2 %d \n", args.pt_x, args.pt_y); */
    start.x = args.pt_x;
    start.y = (double)args.pt_y;

    info1.job = SEQ_RESULT_INFO;
    info1.op = RESULT;
    info1.result = NULL;
    
    seq_result_notify(args.id, (seq_reg_data *)&info1, 0);
    result = (void *)info1.result;
	
    info2.job = SEQ_RESULT_INFO;
    info2.op = OUTPUT;
    info2.result = NULL;
    seq_result_notify(args.id, (seq_reg_data *)&info2, 0);
    output = (void *)info2.result;
    Tcl_GetCommandInfo(interp, output->raster_win, &cmd_info);
    raster = (Tk_Raster*)cmd_info.clientData;
    GetRasterCoords(raster, &wx0, &wy0, &wx1, &wy1);
    RasterGetWorldScroll(raster, &sx0, &sy0, &sx1, &sy1);

    if (output->plot_style == LINE) {
      if (!args.match) {
        nearest = FindNearestLine(result, start, (wx1-wx0)/(wy1-wy0));	
      } else {
        nearest = FindNearestLine(result, start, 1);
      }
    } else {
        /* looking at visual dots ie take scale into account */

      /* NB workshop flags start as rui but is due to structure packing */
        if (!args.match) {
	    nearest = FindNearestMatch(result, start, (wx1-wx0)/(wy1-wy0));
        } else {
            nearest = FindNearestMatch(result, start, 1);
	}
    }
    vTcl_SetResult(interp, "%d %d", nearest.x, (int)nearest.y);
    return TCL_OK;
}


int
tcl_seq_find_result(ClientData clientData, 
		    Tcl_Interp *interp, 
		    int argc, 
		    char *argv[]) 
{
    find_result_arg args;
    int result_id;

    cli_args a[] = {
	{"-raster",   ARG_STR, 1, NULL, offsetof(find_result_arg, raster)},
	{"-seq_id_h", ARG_INT, 1, NULL, offsetof(find_result_arg, seq_id_h)},
	{"-seq_id_v", ARG_INT, 1, NULL, offsetof(find_result_arg, seq_id_v)},
	{NULL,      0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    result_id = seq_find_result(args.raster, args.seq_id_h, args.seq_id_v);
    vTcl_SetResult(interp, "%d", result_id);
    return TCL_OK;
}

int seq_list(ClientData clientData, 
	     Tcl_Interp *interp, 
	     int argc, 
	     char *argv[]) 
{
    list_arg args;
    seq_result *result;
    int seq_num;
    int num_id;
    char **result_id;
    int i;

    cli_args a[] = {
	{"-seq_id",    ARG_INT,   1, NULL, offsetof(list_arg, seq_id)},
	{"-result_id", ARG_STR,   1, NULL, offsetof(list_arg, result_id)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1]))
	return -1;

    if (Tcl_SplitList(interp, args.result_id, &num_id, &result_id) != TCL_OK) {
	return -1;
    }

    seq_num = GetSeqNum(args.seq_id);;
    for (i = 0; i < num_id; i++) {
	result = result_data(atoi(result_id[i]), seq_num);
	result->txt_func(result);
    }

    Tcl_Free((char *)result_id);
    return 0;
}

int tcl_seq_quit_displays(ClientData clientData, 
			  Tcl_Interp *interp, 
			  int argc, 
			  char *argv[]) 
{
    int i;
    seq_reg_delete del;

    del.job = SEQ_QUIT;

    for (i = 0; i < NumSequences(); i++) {
	seq_notify(i, (seq_reg_data *)&del);
    }
    return TCL_OK;
}

/*
 * score matrix, max matches etc setup (formally in sip)
 */

/* 
 * sequence pair display routines
 */
int tcl_seq_pair_display(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{
    display_arg args;
    int seq_disp_id;

    cli_args a[] = {
	{"-window",   ARG_STR, 1, NULL, offsetof(display_arg, window)},
	{"-seq_id_h", ARG_INT, 1, NULL, offsetof(display_arg, seq_id_h)},
	{"-seq_id_v", ARG_INT, 1, NULL, offsetof(display_arg, seq_id_v)},
	{"-cursor_id_h", ARG_INT, 1, NULL, offsetof(display_arg, cursor_id_h)},
	{"-cursor_id_v", ARG_INT, 1, NULL, offsetof(display_arg, cursor_id_v)},
	{"-result_id", ARG_INT, 1, NULL, offsetof(display_arg, result_id)},
	{"-x",        ARG_INT, 1, NULL, offsetof(display_arg, wx)},
	{"-y",        ARG_INT, 1, NULL, offsetof(display_arg, wy)},
	{NULL,        0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    seq_disp_id = seq_disp_reg(interp, args.window, args.seq_id_h, 
			       args.seq_id_v, args.cursor_id_h, 
			       args.cursor_id_v, args.result_id, args.wx, 
			       args.wy);

    vTcl_SetResult(interp, "%d", seq_disp_id);
    return TCL_OK;
}

int
tcl_seq_pair_move_cursor(ClientData clientData, 
			 Tcl_Interp *interp, 
			 int argc, 
			 char *argv[]) 
{
    move_cursor_arg args;
    seq_pair_disp_result *result;
    seq_cursor_notify cn;
    seq_reg_info info;

    cli_args a[] = {
	{"-seqdisp_id", ARG_INT, 1, NULL, offsetof(move_cursor_arg, seqdisp_id)},
	{"-direction", ARG_INT, 1, NULL, offsetof(move_cursor_arg, direction)},
	{"-pos",       ARG_INT, 1, NULL, offsetof(move_cursor_arg, pos)},
	{NULL,      0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;    

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(args.seqdisp_id, (seq_reg_data *)&info, 0);
    if (!info.result) {
	return -1;
    }
    result = (seq_pair_disp_result *)info.result;
   
#ifdef REMOVE
    result->cursor[args.direction]->prev_pos = result->cursor[args.direction]->abspos;
#endif
    
    result->prev_pos[args.direction] = result->cursor[args.direction]->abspos;

    result->cursor[args.direction]->abspos = args.pos;
    result->cursor[args.direction]->job = CURSOR_MOVE;

    /* result->cursorPos = args.pos; */
    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = result->cursor[args.direction];
    cn.cursor->job = CURSOR_MOVE;
    seq_notify(GetSeqNum(result->seq_id[args.direction]), (seq_reg_data *)&cn); 

    return TCL_OK;
}


/*
 * updates the sequence displays when the cursor is moved within the raster
 * display
 */
int
tcl_update_seq_pair(ClientData clientData, 
		    Tcl_Interp *interp, 
		    int argc, 
		    char *argv[]) 
{
    update_seqs_arg args;
    char *seq1;
    char *seq2;
    int seq1_len;
    int seq2_len;
    int s1;
    int s2;
    seq_reg_info info;
    seq_result *result;
    int index1, index2;
    int seq1_type;

    cli_args a[] = {
	{"-win_diff",  ARG_STR, 1, NULL, offsetof(update_seqs_arg, win_diff)},
	{"-win_1", ARG_STR, 1, NULL, offsetof(update_seqs_arg, win_1)},
	{"-win_2", ARG_STR, 1, NULL, offsetof(update_seqs_arg, win_2)},
	{"-left1",  ARG_INT, 1, NULL, offsetof(update_seqs_arg, left1)},
	{"-left2",  ARG_INT, 1, NULL, offsetof(update_seqs_arg, left2)},
	{"-win_len", ARG_INT, 1, NULL, offsetof(update_seqs_arg, win_len)},
	{"-result_index", ARG_INT, 1, NULL, offsetof(update_seqs_arg, 
						     result_index)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

#ifdef DEBUG
    printf("SipUpdateSeqs\n");
#endif

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(args.result_index, (seq_reg_data *)&info, 0);
    result = (seq_result *)info.result;
    if (result == NULL) {
	verror(ERR_WARN, "update_seq_pair", "no result data found for %d",
	     args.result_index);
	return TCL_OK;
    }

    index1 = GetSeqNum(result->seq_id[0]);
    index2 = GetSeqNum(result->seq_id[1]);
    seq1_type=GetSeqType(index1);

    /* if either index is -1; the corresponding seq must have been deleted */
    if ((index1 == -1) || (index2 == -1)) {
	return TCL_OK;
    }
    /* get sequences */
    seq1 = GetSeqSequence(index1);
    seq2 = GetSeqSequence(index2);
    seq1_len = GetSeqLength(index1);
    seq2_len = GetSeqLength(index2);

    /* change from BASE coords to ARRAY coords */
    s1 = args.left1 - 1;
    s2 = args.left2 - 1;

    update_seqs(interp, args.win_1, args.win_2, args.win_diff, seq1, seq2, 
		seq1_len, seq2_len, s1, s2, args.win_len, seq1_type);
    return TCL_OK;
}

/*
 * raster graph plot functions (formally nip raster)
 */

int tcl_raster_select_cursor_graph(ClientData clientData, 
				   Tcl_Interp *interp, 
				   int argc, 
				   char *argv[]) 
{
    raster_scursor_arg args;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    int cursor_id;
    int max_dist;

    cli_args a[] = {
	{"-id",     ARG_INT, 1, NULL, offsetof(raster_scursor_arg, raster_id)},
	{"-raster", ARG_STR, 1, NULL, offsetof(raster_scursor_arg, raster)},
	{"-x",      ARG_INT, 1, NULL, offsetof(raster_scursor_arg, rx)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    if (Tcl_GetCommandInfo(interp, args.raster, &info) == 0) 
	return TCL_ERROR;
    raster = (Tk_Raster*)info.clientData;

    max_dist = get_default_int(interp, spin_defs, 
			       w("GRAPH.CURSOR.SENSITIVE"));
    cursor_id = raster_select_cursor_graph(args.raster_id, raster, args.raster,
					 args.rx, max_dist);
    vTcl_SetResult(interp, "%d", cursor_id);
    return TCL_OK;
}

int tcl_get_raster_frame_graph(ClientData clientData, 
			       Tcl_Interp *interp, 
			       int argc, 
			       char *argv[]) 
{    
    raster_win_arg args;
    char *raster_win;
    int type;

    cli_args a[] = {
	{"-seq_id", ARG_INT, 1, NULL, offsetof(raster_win_arg, seq_id)},
	{"-type",   ARG_STR, 1, "",   offsetof(raster_win_arg, type)},
	{"-frame",  ARG_INT, 1, "0",  offsetof(raster_win_arg, frame)},
	{NULL,      0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (strcmp(args.type, "") == 0) {
	type = -1;
    } else {
	if (strcmp(args.type, "GENESEARCH") == 0) {
	    type = SEQ_TYPE_GENESEARCH;
	} else if (strcmp(args.type, "STRINGSEARCH") == 0) {
	    type = SEQ_TYPE_STRINGSEARCH;
	} else if (strcmp(args.type, "RESTRICTION") == 0) {
	    type = SEQ_TYPE_RESTRICTION;
	} else if (strcmp(args.type, "BASECOMP") == 0) {
	    type = SEQ_TYPE_BASECOMP;
	} else if (strcmp(args.type, "CODONPREF") == 0) {
	    type = SEQ_TYPE_CODONPREF;
	} else if (strcmp(args.type, "AUTHOR") == 0) {
	    type = SEQ_TYPE_AUTHOR;
	} else if (strcmp(args.type, "BASEBIAS") == 0) {
	    type = SEQ_TYPE_BASEBIAS;
	} else if (strcmp(args.type, "TRNA") == 0) {
	    type = SEQ_TYPE_TRNA;
	} else if (strcmp(args.type, "STOPCODON") == 0) {
	    type = SEQ_TYPE_STOPCODON;
	} else if (strcmp(args.type, "STARTCODON") == 0) {
	    type = SEQ_TYPE_STARTCODON;
	} else if (strcmp(args.type, "SPLICE") == 0) {
	    type = SEQ_TYPE_SPLICE;
	} else if (strcmp(args.type, "WTMATRIXSEARCH") == 0) {
	    type = SEQ_TYPE_WTMATRIXSEARCH;
	} else {
	    verror(ERR_WARN, "NipGetRasterFrame", "Unrecognised type");
	    return TCL_OK;
	}
    }
    raster_win = get_raster_frame_graph(interp, args.seq_id, type, 
					args.frame); 
    vTcl_SetResult(interp, "%s", raster_win);
    xfree(raster_win);
    return TCL_OK;
}

int tcl_seqed_display(ClientData clientData, 
		      Tcl_Interp *interp, 
		      int argc, 
		      char *argv[])
{
    seq_disp_arg args;
    char *seq;
    char *sequence;
    int seq_len;
    int seq_num;
    int seqed_id;

    cli_args a[] = {
	{"-window", ARG_STR, 1, NULL, offsetof(seq_disp_arg, seqed_win)},
	{"-seq_id", ARG_INT, 1, NULL, offsetof(seq_disp_arg, seq_id)},
	{NULL,      0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    seq_num = GetSeqNum(args.seq_id);
    seq = GetSeqSequence(seq_num);
    seq_len = GetSeqLength(seq_num);

    if (NULL == (sequence = (char *)xmalloc((seq_len+1) * sizeof(char))))
	return TCL_OK;
    strncpy(sequence, seq, seq_len);
    sequence[seq_len] = '\0';
    seqed_id = add_seq_seqed(interp, sequence, args.seqed_win, seq_num); 

    xfree(sequence);
    vTcl_SetResult(interp, "%d", seqed_id);
    return TCL_OK;
}


/* 
 * raster dot plot functions (formally sip raster) 
 */

int tcl_raster_select_cursor_dot(ClientData clientData, 
				 Tcl_Interp *interp, 
				 int argc, 
				 char *argv[]) 
{
    raster_scursor_arg args;
    Tcl_CmdInfo info;
    Tk_Raster *raster;
    int cursor_id_h, cursor_id_v;
    int max_dist;

    cli_args a[] = {
	{"-id",     ARG_INT, 1, NULL, offsetof(raster_scursor_arg, raster_id)},
	{"-raster", ARG_STR, 1, NULL, offsetof(raster_scursor_arg, raster)},
	{"-x",      ARG_INT, 1, NULL, offsetof(raster_scursor_arg, rx)},
	{"-y",      ARG_INT, 1, NULL, offsetof(raster_scursor_arg, ry)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;
    
    if (Tcl_GetCommandInfo(interp, args.raster, &info) == 0) 
	return TCL_ERROR;
    raster = (Tk_Raster*)info.clientData;

    max_dist = get_default_int(interp, spin_defs, 
			       w("DOT.CURSOR.SENSITIVE"));
    if (-1 == (raster_select_cursor_dot(args.raster_id, raster, 
					args.raster,
					args.rx, args.ry, max_dist, 
					&cursor_id_h, &cursor_id_v))) {
	verror(ERR_WARN, "raster_select_cursor_dot", "unable to allocate memory");
	vTcl_SetResult(interp, "-1 -1 -1");
    } else {
	vTcl_SetResult(interp, "%d %d", cursor_id_h, cursor_id_v);
    }
    return TCL_OK;
}

int tcl_get_raster_frame_dot(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc, 
			     char *argv[]) 
{    
    raster_frame_arg args;
    char raster_win[1024];

    cli_args a[] = {
	{"-seq_id_h", ARG_INT, 1, NULL, offsetof(raster_frame_arg, seq_id_h)},
	{"-seq_id_v", ARG_INT, 1, NULL, offsetof(raster_frame_arg, seq_id_v)},
	{NULL,      0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;    

    if (-1 == (get_raster_frame_dot(interp, args.seq_id_h, args.seq_id_v, raster_win))) {
	verror(ERR_WARN, "Failure in get_raster_frame_dot", "Unable to allocate "
	       "a raster frame");
	return TCL_OK;
    }

    vTcl_SetResult(interp, "%s", raster_win);
    return TCL_OK;
}

int tcl_seq_raster_reg(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, 
		       char *argv[]) 
{    
    reg_arg args;
    char **array;
    char **array1;
    int num_seq_id;
    seq_id_dir *seq_array;
    int i, id, num;

    cli_args a[] = {
	{"-window",  ARG_STR, 1, NULL, offsetof(reg_arg, raster)},
	{"-seq_id",  ARG_STR, 1, NULL, offsetof(reg_arg, seq_id)},
	{NULL,         0,       0, NULL, 0}
    };
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if (Tcl_SplitList(interp, args.seq_id, &num_seq_id, &array) != TCL_OK) {
	return TCL_ERROR;
    }
#ifdef DEBUG
    printf("max_num_Seq %d\n", MAX_NUM_SEQ);
#endif
    if (NULL == (seq_array = (seq_id_dir *)xmalloc(MAX_NUM_SEQ * 
						   sizeof(seq_id_dir))))
	return TCL_ERROR;

    for (i = 0; i < num_seq_id; i++) {
	if (Tcl_SplitList(interp, array[i], &num, &array1) != TCL_OK) {
	    return TCL_ERROR;
	}
	if (num != 2) {
	    printf("Error in RasterReg\n");
	    return TCL_OK;
	}

	seq_array[i].seq_id = atoi(array1[0]);
	seq_array[i].direction = atoi(array1[1]);
    }

    id = seq_raster_reg(interp, args.raster, seq_array, num_seq_id);
    vTcl_SetResult(interp, "%d", id);

    Tcl_Free((char *)array);
    Tcl_Free((char *)array1);
    return TCL_OK;
}

int tcl_add_seq_to_raster(ClientData clientData, 
			  Tcl_Interp *interp, 
			  int argc, 
			  char *argv[]) 
{
    seqtoraster_arg args;
    int seq_num;
    RasterResult *result;

    cli_args a[] = {
	{"-raster_id", ARG_INT, 1, NULL, offsetof(seqtoraster_arg, raster_id)},
	{"-seq_id",    ARG_INT, 1, NULL, offsetof(seqtoraster_arg, seq_id)},
	{"-direction", ARG_INT, 1, NULL, offsetof(seqtoraster_arg, direction)},
	{"-line_width",ARG_INT, 1, NULL, offsetof(seqtoraster_arg, line_width)},
	{NULL,     0,  0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    result = raster_id_to_result(args.raster_id);
    seq_num = GetSeqNum(args.seq_id);
    add_seq_to_raster(result, args.seq_id, seq_num, args.direction, 
		      args.line_width, seq_raster_callback);

    return TCL_OK;
}


/*
 * emboss functions
 */

int emboss_create(ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, 
		  char *argv[]) 
{
    emboss_arg args;
    int id;

    cli_args a[] = {
	{"-seq_id_h",  ARG_INT, 1, NULL, offsetof(emboss_arg, seq_id_h)},
	{"-start_h",   ARG_INT, 1, "1",  offsetof(emboss_arg, start_h)},
	{"-end_h",     ARG_INT, 1, "-1", offsetof(emboss_arg, end_h)},
	{"-seq_id_v",  ARG_INT, 1, "-1", offsetof(emboss_arg, seq_id_v)},
	{"-start_v",   ARG_INT, 1, "1",  offsetof(emboss_arg, start_v)},
	{"-end_v",     ARG_INT, 1, "-1", offsetof(emboss_arg, end_v)},
	{"-graph",     ARG_INT, 1, NULL, offsetof(emboss_arg, graph)},
	{"-data",      ARG_STR, 1, NULL, offsetof(emboss_arg, data)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "emboss", "unable to parse input\n");
	return TCL_ERROR;
    }
 
    if (args.graph == SEQ_GRAPH) {
	if (-1 == init_emboss_graph_create(interp, args.seq_id_h, 
					   args.start_h, args.end_h,  
					   args.data, &id)) {
	    vTcl_SetResult(interp, "%d", -1);
	    return TCL_OK;
	}
    } else if (args.graph == SEQ_DOT) {
	if (-1 == init_emboss_dot_create(interp, args.seq_id_h, args.start_h,
					 args.end_h, args.seq_id_v, 
					 args.start_v, args.end_v, args.data, 
					 &id)) {
	    vTcl_SetResult(interp, "%d", -1);
	    return TCL_OK;
	}
    } else if (args.graph == SEQ_STICK) {
	if (-1 == init_emboss_stick_create(interp, args.seq_id_h, args.start_h,
					   args.end_h, args.data, &id)) {
	    vTcl_SetResult(interp, "%d", -1);
	    return TCL_OK;
	}
    }

    vTcl_SetResult(interp, "%d", id);
    return TCL_OK;
}

int emboss_plot(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{
    plot_arg args;

    cli_args a[] = {
	{"-seq_id_h",  ARG_INT, 1, NULL, offsetof(plot_arg, seq_id_h)},
	{"-seq_id_v",  ARG_INT, 1, "-1", offsetof(plot_arg, seq_id_v)},
	{"-result_id", ARG_INT, 1, NULL, offsetof(plot_arg, result_id)},
	{"-graph",     ARG_INT, 1, NULL, offsetof(plot_arg, graph)},
	{"-name",      ARG_STR, 1, NULL, offsetof(plot_arg, name)},
	{"-window",    ARG_STR, 1, NULL, offsetof(plot_arg, raster)},
	{"-window_id", ARG_INT, 1, NULL, offsetof(plot_arg, raster_id)},
	{"-fill",      ARG_STR, 1, NULL, offsetof(plot_arg, colour)},
	{"-width",     ARG_INT, 1, "1",  offsetof(plot_arg, line_width)},
	{NULL,     0,       0, NULL, 0}
    };

    if (-1 == parse_args(a, &args, argc-1, &argv[1])) {
	verror(ERR_WARN, "emboss", "unable to parse input\n");
	return TCL_ERROR;
    }
    if (args.graph == 0) {
	if (-1 == init_emboss_graph_plot(interp, args.seq_id_h, args.result_id,
					 args.name, args.raster, 
					 args.raster_id,
					 args.colour, args.line_width)) {
	    vTcl_SetResult(interp, "%d", -1);
	    return TCL_OK;
	}
    } else if (args.graph == 1) {
	if (-1 == init_emboss_dot_plot(interp, args.seq_id_h, args.seq_id_v,
				       args.result_id, args.name,
				       args.raster, args.raster_id,
				       args.colour, args.line_width)) {
	    vTcl_SetResult(interp, "%d", -1);
	    return TCL_OK;
	}
    }
    return TCL_OK;
}

int tcl_emboss(ClientData clientData, 
	       Tcl_Interp *interp, 
	       int argc, 
	       char *argv[]) 
{
    if (strcmp(argv[1], "create") == 0) {
	emboss_create(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "list") == 0) {
	seq_list(clientData, interp, argc, argv);
    } else if (strcmp(argv[1], "plot") == 0) {
	emboss_plot(clientData, interp, argc, argv);
    }
    return TCL_OK;
}

int tcl_INT_MAX(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{    
    vTcl_SetResult(interp, "%d", INT_MAX);
    return TCL_OK;
}

int tcl_INT_MIN(ClientData clientData, 
		Tcl_Interp *interp, 
		int argc, 
		char *argv[]) 
{    
    vTcl_SetResult(interp, "%d", INT_MIN);
    return TCL_OK;
}

typedef struct {
    char *filename;
} read_enz_arg;

int tcl_read_enz_file(ClientData clientData,
		      Tcl_Interp *interp, 
		      int argc, 
		      char *argv[])
{
    read_enz_arg args;
    int num_enzymes;
    int i;
    char **names;

    cli_args a[] = {
	{"-file",   ARG_STR, 1, NULL, offsetof(read_enz_arg, filename)},
	{NULL,	    0,	     0, NULL, 0}
    };
    
    if (-1 == parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    /* open_r_enz_file(args.filename, &r_enzyme, &num_enzymes); */
    if (0 == r_enz_file_names(args.filename, &names, &num_enzymes))
	return TCL_OK;
    
    for (i =0; i < num_enzymes; i++) {
	Tcl_AppendElement(interp, names[i]);
	xfree(names[i]);
    }

    if (num_enzymes)
	xfree(names);

    return TCL_OK;
}

/* set up tcl commands which call C procedures */
/*****************************************************************************/
/*                        Spin_Init                                          */
/*****************************************************************************/
int
Spin_Init(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "seq_result_names", tcl_seq_result_names, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_get_ops", tcl_seq_get_ops, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_invoke_op", tcl_seq_invoke_op, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_get_seq_ops", tcl_seq_get_seq_ops, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_invoke_seq_op", tcl_seq_invoke_seq_op, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_result_update", SeqResultUpdate, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_result_key_name", SeqResultKeyName, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_get_brief", SeqGetBrief, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_result_info", tcl_seq_result_info, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "create_cursor", CreateCursor,
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "delete_cursor", DeleteCursor,
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "cursor_notify", CursorNotify, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "cursor_ref", CursorRef, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "query_cursor", QueryCursor,
		      (ClientData) NULL,
		      NULL);
    Tcl_CreateCommand(interp, "get_archive_list", GetArchiveList, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_info", tcl_seq_info, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "set_seq_structure", tcl_set_seq_structure, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_seq_id", tcl_GetSeqId, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_seq_num", tcl_GetSeqNum, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "name_to_seq_id", NameToSeqId, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_active_seq_id", GetActiveSeqId, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_active_seq_name", GetActiveSeqName, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "num_sequences", tcl_NumSequences, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "s_length", tcl_s_length, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_file_save", SeqFileSave, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_file_delete", SeqFileDelete, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_complement", SeqComplement, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_set_active_seq", SeqSetActiveSeq, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_interconvert", SeqInterconvert, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_scramble", SeqScramble, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_rotate", SeqRotate, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_translate_seq", SeqTranslateSequence, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_set_range", SeqSetRange, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_copy_range", SeqCopyRange, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_sender", SeqSender, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "read_sequence", tcl_read_sequence, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "sequence_names", tcl_sequence_names, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_result_seq_id", tcl_get_result_seq_id,
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_quit_displays", tcl_seq_quit_displays,
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "nearest_match", tcl_nearest_match, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_find_result", 
		      tcl_seq_find_result, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_pair_display", 
		      tcl_seq_pair_display, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_pair_move_cursor", 
		      tcl_seq_pair_move_cursor, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "update_seq_pair", tcl_update_seq_pair, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "raster_select_cursor_graph", 
		      tcl_raster_select_cursor_graph, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_raster_frame_graph", 
		      tcl_get_raster_frame_graph, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seqed_display", 
		      tcl_seqed_display,
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "raster_select_cursor_dot", 
		      tcl_raster_select_cursor_dot, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_raster_frame_dot", 
		      tcl_get_raster_frame_dot, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "seq_raster_reg", 
		      tcl_seq_raster_reg, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "add_seq_to_raster", 
		      tcl_add_seq_to_raster,
		      (ClientData) NULL, 
	              (Tcl_CmdDeleteProc *) NULL);
    /*yy
       Tcl_CreateCommand(interp, "set_score_matrix", SetScoreMatrix, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);   
    */
    Tcl_CreateCommand(interp, "emboss", tcl_emboss, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "INT_MAX", tcl_INT_MAX, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "INT_MIN", tcl_INT_MIN, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "read_enz_file", tcl_read_enz_file,
		      (ClientData) NULL,
		      NULL);

   {
        char *s, c[10];
	/*
	 * Set packages(name). This is done to prevent subsequent reloading
	 * of this library (for efficiency reasons). The only reason that this
	 * is necessary is that currently gap4 dynamically links with some
	 * libraries at link time. When they're all at run time this won't
	 * be necessary.
	 */
	if (s = Tcl_GetVar2(interp, "packages", "spin", 
			    TCL_GLOBAL_ONLY))
	    sprintf(c, "%d", atoi(s)|2);
	else
	    strcpy(c, "2");
	Tcl_SetVar2(interp, "packages", "spin", c, TCL_GLOBAL_ONLY);
    }

    spin_init_globals(interp);
    Sip_Init(interp);
    Nip_Init(interp);
    return RasterUtils_Init(interp);
}

int
Spin_SafeInit(Tcl_Interp *interp) {
    return Spin_Init(interp);
}

int Spin_Unload(Tcl_Interp *interp, int flags) {
    Tcl_SetResult(interp, "Pkg_Unload() function not implemented",
		  TCL_STATIC);
    return TCL_ERROR;
}

int Spin_SafeUnload(Tcl_Interp *interp, int flags) {
    Tcl_SetResult(interp, "Pkg_SafeUnload() function not implemented",
		  TCL_STATIC);
    return TCL_ERROR;
}
