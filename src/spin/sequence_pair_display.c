#include <string.h>
#include <ctype.h>

#include "seq_reg.h"
#include "sequence_pair_display.h"
#include "seq_raster.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "seq_results.h"
#include "xalloc.h"
#include "spin_globals.h"

/*yy*/
#include "sequence_formats.h"
#include "readpam.h"
#include "sip_results.h"

extern int **score_matrix;    
/*yy*/

void seq_disp_callback(int seq_num, void *obj, seq_reg_data *jdata);

void seq_disp_move_cursor(Tcl_Interp *interp,
			  char *seq_disp_win,
			  int result_id,
			  int pos,
			  int direction)
{
    char cmd[1024];
    sprintf(cmd, "seq_disp_show_cursor %s %d %d %d\n", seq_disp_win, result_id,
	    pos, direction);
    if (TCL_OK != Tcl_Eval(interp, cmd)) {
	printf("seq_disp_move_cursor: %s\n",Tcl_GetStringResult(interp));
    }
}

void seq_disp_shutdown(Tcl_Interp *interp,
		       seq_pair_disp_result *result,
		       int seq_num)
{
    int seq_num_h;
    int seq_num_v;
    char *tmp;
    
    /* need to set prev_pos */
#ifdef REMOVE
    result->cursor[HORIZONTAL]->prev_pos = result->cursor[HORIZONTAL]->abspos;
    result->cursor[VERTICAL]->prev_pos = result->cursor[VERTICAL]->abspos;
#endif
    result->prev_pos[HORIZONTAL] = result->cursor[HORIZONTAL]->abspos;
    result->prev_pos[VERTICAL] = result->cursor[VERTICAL]->abspos;

    
    seq_num_h = GetSeqNum(result->seq_id[HORIZONTAL]);
    seq_num_v = GetSeqNum(result->seq_id[VERTICAL]);
    
    seq_deregister(seq_num_h, seq_disp_callback, (seq_pair_disp_result *)result);
    seq_deregister(seq_num_v, seq_disp_callback, (seq_pair_disp_result *)result);

    delete_cursor(seq_num_h, result->cursor[HORIZONTAL]->id, 1);
    delete_cursor(seq_num_v, result->cursor[VERTICAL]->id, 1);
    
    tmp = get_default_string(result->interp, tk_utils_defs, w("RASTER.RESULTS.WIN"));
    if (TCL_OK != Tcl_VarEval(result->interp, "seq_result_list_update ", 
			      tmp, NULL)){
	verror(ERR_WARN, "seq disp shutdown", "%s \n", Tcl_GetStringResult(result->interp));
    }
  
    xfree(result);
}

void seq_disp_callback(int seq_num, void *obj, seq_reg_data *jdata) 
{
    seq_pair_disp_result *result = (seq_pair_disp_result *) obj;
    char *seq_disp_win = result->seq_disp_win;

    switch(jdata->job) {
    case SEQ_QUERY_NAME:
	{
	    sprintf(jdata->name.line, "sequence display");
	    return;
	}
    case SEQ_QUIT:
    case SEQ_DELETE:
	{
	    seq_disp_shutdown(result->interp, result, seq_num);
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
	    seq_disp_shutdown(result->interp, result, seq_num);
	    return;
	}
	break;
    case SEQ_CURSOR_NOTIFY: 
	{
	    cursor_t *cursor = (cursor_t *)jdata->cursor_notify.cursor;

#ifdef DEBUG
	    printf("SEQ_CURSOR_NOTIFY seqed pos %d id %d result id %d\n", 
		   cursor->abspos, cursor->id, result->cursor[HORIZONTAL]->id);
#endif	    
	    if (result->cursor[HORIZONTAL]->id == cursor->id) {
		seq_disp_move_cursor(result->interp, result->seq_disp_win, 
				     result->result_id, cursor->abspos,
				     HORIZONTAL);
	    }
	    if (result->cursor[VERTICAL]->id == cursor->id) {
		seq_disp_move_cursor(result->interp, result->seq_disp_win, 
				     result->result_id, cursor->abspos,
				     VERTICAL);
	    }
	    break;
	}
    case SEQ_RESULT_INFO:
	switch (jdata->info.op) {
	case WINDOW:
	    jdata->info.result = (void *)seq_disp_win;
	    break;
	case RESULT:
	    jdata->info.result = (void *)result;
	    break;
	}
	break;
    }
}

/*
 * register sequence display widget
 */
int seq_disp_reg(Tcl_Interp *interp,
		 char *seq_disp_win,
		 int seq_id_h,
		 int seq_id_v,
		 int cursor_id_h,
		 int cursor_id_v,
		 int result_id,
		 int wx,
		 int wy)
{
    int id;
    seq_pair_disp_result *seq_disp_result;
    seq_cursor_notify cn;
    cursor_t *cursor;
    int seq_num;
    int line_width;

    if (NULL == (seq_disp_result = (seq_pair_disp_result *)xmalloc(sizeof(seq_pair_disp_result))))
	return -1;

    seq_disp_result->op_func = seq_disp_callback;
    seq_disp_result->seq_id[HORIZONTAL] = seq_id_h;
    seq_disp_result->seq_id[VERTICAL] = seq_id_v;
    seq_disp_result->result_id = result_id;

    strcpy(seq_disp_result->seq_disp_win, seq_disp_win);
    seq_disp_result->interp = interp;

    id = get_reg_id();
    seq_disp_result->index = id;

    line_width = get_default_int(interp, spin_defs, w("SEQ.CURSOR.LINE_WIDTH"));

    if (cursor_id_h < 0) {
	/* no horizontal cursor */
	seq_disp_result->cursor[HORIZONTAL] = create_cursor(GetSeqNum(seq_id_h), 1, NULL, line_width, 1, HORIZONTAL);
	seq_disp_result->cursor_visible[HORIZONTAL] = 0;
    } else {
	seq_num = GetSeqNum(seq_id_h);
	cursor = find_cursor(&seq_num, cursor_id_h, -1);
#ifdef DEBUG
	printf("USING EXISTING CURSOR private %d\n", cursor->private);
#endif
	cursor->private = 1;
	seq_disp_result->cursor[HORIZONTAL] = cursor;
    }

    if (cursor_id_v < 0) {
	/* no vertical cursor */
	seq_disp_result->cursor[VERTICAL] = create_cursor(GetSeqNum(seq_id_v), 1, NULL, line_width, 1, VERTICAL);
	seq_disp_result->cursor_visible[VERTICAL] = 0;
    } else {
	seq_num = GetSeqNum(seq_id_v);
	cursor = find_cursor(&seq_num, cursor_id_v, -1);
#ifdef DEBUG
	printf("USING EXISTING CURSOR private %d\n", cursor->private);
#endif
	cursor->private = 1;
	seq_disp_result->cursor[VERTICAL] = cursor;
    }
 
#ifdef REMOVE
    seq_disp_result->cursor[HORIZONTAL]->prev_pos = seq_disp_result->cursor[HORIZONTAL]->abspos;
    seq_disp_result->cursor[VERTICAL]->prev_pos = seq_disp_result->cursor[VERTICAL]->abspos;
#endif

    seq_disp_result->prev_pos[HORIZONTAL] = seq_disp_result->cursor[HORIZONTAL]->abspos;
    seq_disp_result->prev_pos[VERTICAL] = seq_disp_result->cursor[VERTICAL]->abspos;

    seq_disp_result->cursor[HORIZONTAL]->abspos = wx;
    seq_disp_result->cursor[VERTICAL]->abspos = wy;
    
    seq_register(GetSeqNum(seq_id_h), seq_disp_callback, 
		 (void *)seq_disp_result, SEQ_SEQED, id);
    seq_register(GetSeqNum(seq_id_v), seq_disp_callback, 
		 (void *)seq_disp_result, SEQ_SEQED, id);

    cn.job = SEQ_CURSOR_NOTIFY;
    cn.cursor = seq_disp_result->cursor[HORIZONTAL];
    cn.cursor->job = CURSOR_MOVE;
    seq_notify(GetSeqNum(seq_id_h), (seq_reg_data *)&cn);

    cn.cursor = seq_disp_result->cursor[VERTICAL];
    cn.cursor->job = CURSOR_MOVE;
    seq_notify(GetSeqNum(seq_id_v), (seq_reg_data *)&cn); 

    return id;
}

/*
 * scrolling procedure in sequence display widget
 */
int update_seqs(Tcl_Interp *interp,
		char *window1,
		char *window2,
		char *win_diff,
		char *seq1,
		char *seq2,
		int seq1_len,
		int seq2_len,
		int seq1_left,
		int seq2_left,
		int width,
                int type)
     
{ 
    /* yypp */
    int row;
    int col;
    int row_pam;
    int col_pam; 
    char chap_pam[] = {'A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','Z','X','*','-'};   
    int chap_pam_size = sizeof(chap_pam);

   /* yypp */
   
    char tmp[10];
    int start;
    int end;
    int i;

    static int cur_width = 0;
    static char *tmp2 = NULL;
    static char *tmp_str = NULL;
    int s1_left_tmp, s2_left_tmp;
    int **score_display;

    score_display = get_matrix_file(PROTEIN);
   
    /* 
     * check if current sequence window width is the same as new width ie 
     * if the same amount of sequence is to be displayed
     */
    if (cur_width != width) {
	/* printf("%d %d \n", cur_width, width); */
	if (tmp2) 
	    xfree(tmp2);
	if (tmp_str)
	    xfree(tmp_str);
	if (NULL == (tmp2 = (char *)xmalloc((width+21) * sizeof(char))))
	    return TCL_OK;
	if (NULL == (tmp_str = (char *)xmalloc((width+2) * sizeof(char))))
	    return TCL_OK;
    }
    /* insert padding in line 1 and line 2 */
    if (seq1_left < 0) {
	s1_left_tmp = seq1_left < -width ? width : -seq1_left;
	memset(tmp_str, ' ', s1_left_tmp);
	tmp_str[s1_left_tmp] = '\0';
	Tcl_VarEval(interp, window1, " insert 1.end ", "\"", tmp_str, "\"", NULL);
	Tcl_VarEval(interp, window1, " insert 2.end ", "\"", tmp_str, "\"", NULL);
    }
    if (seq2_left < 0) {
	s2_left_tmp = seq2_left < -width ? width : -seq2_left;
	memset(tmp_str, ' ', s2_left_tmp);
	tmp_str[s2_left_tmp] = '\0';
	Tcl_VarEval(interp, window2, " insert 1.end ", "\"", tmp_str, "\"", NULL);
	Tcl_VarEval(interp, window2, " insert 2.end ", "\"", tmp_str, "\"", NULL);
    }

    /* insert numbers in line 1 of win1 and line 2 of win2 */
    tmp[0] = '\0';
    tmp2[0] = '\0';

    if (seq1_left < 0) 
	start = 1;
    else
	start = seq1_left / 10;

/*
    if ((seq1_left + width) > seq1_len)
	end = seq1_len/ 10; 
    else

*/
    end = (width + seq1_left) / 10;

    for (i = start; i <= end; i++) {
	sprintf(tmp, "%10d", i*10);
	strcat(tmp2, tmp);
    }
    tmp_str[0] = '\0';

    if (seq1_left < 0)
	strncpy(tmp_str, &tmp2[0], width);
     else
	strncpy(tmp_str, &tmp2[seq1_left%10+10], width);
    
    tmp_str[width] = '\0';
    Tcl_VarEval(interp, window1, " insert 1.end ", "\"", tmp_str, "\"", NULL);

    /* sequence 2 numbers */
    tmp2[0] = '\0';
    tmp_str[0] = '\0';
    if (seq2_left <= 0) 
	start = 1;
    else 
	start = seq2_left / 10;
/*
    if ((seq2_left + width) > seq2_len)
	end = seq2_len/ 10; 
    else
*/
    end = (width + seq2_left) / 10;

    for (i = start; i <= end; i++) {
	sprintf(tmp, "%10d", i*10);
	strcat(tmp2, tmp);
    }
    if (seq2_left < 0) 
	strncpy(tmp_str, &tmp2[0], width);
    else
	strncpy(tmp_str, &tmp2[seq2_left%10+10], width);
    tmp_str[width] = '\0';

    Tcl_VarEval(interp, window2, " insert 2.end ", "\"", tmp_str, "\"", NULL);

    /* insert seq in line 2 of win1 and line 1 of win2 */
    tmp_str[0] = '\0';
    memset(tmp_str, ' ', width);

    if (seq1_left < 0) {
	if (seq1_left > -width) {
	    strncpy(tmp_str, &seq1[0], width + seq1_left);
	    tmp_str[width + seq1_left] = '\0';
	} else {
	    tmp_str[width] = '\0';
	}
    } else { 
	if (seq1_left + width <= seq1_len) {
	    strncpy(tmp_str, &seq1[seq1_left], width);
	} else {
	    if (seq1_left < seq1_len) {
		strcpy(tmp_str, &seq1[seq1_left]);
		tmp_str[strlen(tmp_str)] = ' ';
	    }
	}
	tmp_str[width] = '\0';
    }
    Tcl_VarEval(interp, window1, " insert 2.end {", tmp_str, "}", NULL);

    /* printf("%s\n", tmp_str); */
    /* sequence 2 */
    tmp_str[0] = '\0';
    memset(tmp_str, ' ', width);

    if (seq2_left < 0) {
	if (seq2_left > -width) {
	    strncpy(tmp_str, &seq2[0], width + seq2_left);
	    tmp_str[width + seq2_left] = '\0';
	} else {
	    tmp_str[width] = '\0';
	}
    } else {
	if (seq2_left + width <= seq2_len) {
	    strncpy(tmp_str, &seq2[seq2_left], width);
	} else {
	    if (seq2_left < seq2_len) {
		strcpy(tmp_str, &seq2[seq2_left]);
		tmp_str[strlen(tmp_str)] = ' ';
	    }
	}
	tmp_str[width] = '\0';
    }
    Tcl_VarEval(interp, window2, " insert 1.end {", tmp_str, "}", NULL);

    /* diffs line */
    tmp_str[0] = '\0';
    for (i = 0; i < width; i++, seq1_left++, seq2_left++) {
	if ((seq1_left < 0) || (seq2_left < 0) || (seq1_left > (seq1_len-1)) 
	    || (seq2_left > (seq2_len-1))){
	      tmp_str[i] = '!'; 
	} else {

	  /*yypp */

	    if ( type ==PROTEIN ){
		if (toupper(seq1[seq1_left]) == toupper(seq2[seq2_left])) {
		    tmp_str[i] = symbol_align0[0];
		} else { row_pam = -1; col_pam = -1;
		for( row=0; row < chap_pam_size; row++){
		    if(toupper(seq1[seq1_left]) == toupper(chap_pam[row])){
			row_pam = row;
			break;
		    }
		}      
		if(row_pam < 0){
		      verror(ERR_WARN, "Update Seqs", 
                      "Sequence contains a character that is not mentioned in score matrix\n");
           return -1;
		}
		for( col=0; col < chap_pam_size; col++){
		    if(toupper(seq2[seq2_left]) == toupper(chap_pam[col])){
			col_pam = col;
			break;
		    }
		}      
		if(col_pam < 0){
		    verror(ERR_WARN, "Update Seqs", 
                   "Sequence contains a character that is not mentioned in score matrix\n");
		    return -1;
		}
		if(score_display[row_pam][col_pam] > cutoff_align3 && 
                   score_display[row_pam][col_pam] <= cutoff_align2 )
		     tmp_str[i] = symbol_align3[0];
		else if (score_display[row_pam][col_pam] > cutoff_align2 && 
                         score_display[row_pam][col_pam] <= cutoff_align1 )
		    tmp_str[i] = symbol_align2[0];
		else if (score_display[row_pam][col_pam] > cutoff_align1 )
		    tmp_str[i] = symbol_align1[0];
		else
		    tmp_str[i] = ' ';
		}
	    }
	    if(type == DNA){ 
		if (toupper(seq1[seq1_left]) == toupper(seq2[seq2_left])) 
		    tmp_str[i] = '|'; 
		else
		    tmp_str[i] = ' '; 
   
	    }
	}
    }
    tmp_str[width] = '\0';
    Tcl_VarEval(interp, win_diff, " insert end ", "\"", tmp_str, "\"", NULL);
    cur_width = width;
    return TCL_OK;
}

/*
 * create the sequence display widget
 */
void SequencePairDisplay(Tcl_Interp *interp,
			 char *raster_win,
			 int result_index,
			 int seq_id_h,
			 int seq_id_v)
{
    char cmd[1024];
    int raster_id;
    int height;

    Tcl_VarEval(interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(interp));
    
    Tcl_VarEval(interp, "winfo height ", raster_win, NULL);
    height = atoi(Tcl_GetStringResult(interp));
   
    sprintf(cmd, "SequencePairDisplay 1 1 %d %d -1 -1 %d\n", 
	    seq_id_h, seq_id_v, result_index);

       /* dmalloc_verify(NULL); */
 
    if (TCL_OK != Tcl_GlobalEval(interp, cmd)) {
	printf("DisplaySequences: %s\n", Tcl_GetStringResult(interp));
    }
}

/*
 * destroy the sequence display widget
 */
void DestroySequencePairDisplay(Tcl_Interp *interp,
				int result_index)
{
    char cmd[1024];
    char *seq_disp_win;

    seq_disp_win = get_default_string(interp, spin_defs, 
				      "SEQ_DISP.WIN");
    sprintf(cmd, "SeqDispStartShutdown %s%d", seq_disp_win, result_index);
    if (TCL_ERROR == (Tcl_Eval(interp, cmd))) {
	printf("DestroyDisplaySequences %s\n", Tcl_GetStringResult(interp));
    }
}

/********/
int spin_list_alignment ( char *seq1, char *seq2, char *name1, char *name2,
		     int pos1, int pos2, char *title, int type )

{
    char *match_syms ;
    int i,j,k,seq_len,p1,p2,line_length=60;
    int l, spads1,spads2, p11,p22;
    int row,col;
    int row_pam,col_pam;

    char chap_pam[] = {'A','B','C','D','E','F','G','H','I','K','L','M','N',
                       'P','Q','R','S','T','V','W','Y','Z','X','*','-','.', ' '}; 
    int chap_pam_size = sizeof(chap_pam);

     p11=pos1;
     p22=pos2;
     seq_len = strlen(seq1);
     if ( ! (match_syms = (char *) xmalloc ( sizeof(char) * (seq_len) + 1))) {
		return -1;
     }	
     for (i=0; i<seq_len; i++){
	 if (type == PROTEIN){
	     if (toupper(seq1[i]) == toupper(seq2[i])) {
		 match_syms[i] = symbol_align0[0];
	     }else 
		 if(seq1[i]=='.'||seq2[i]=='.'){
		     match_syms[i] = ' ';
		 } else { row_pam = -1; col_pam = -1;
		 for( row=0; row < chap_pam_size; row++){
		     if(toupper(seq1[i]) == toupper(chap_pam[row])){
			 row_pam = row;
			 break;
		     }
		 }
     
		 if(row_pam < 0){
		      verror(ERR_WARN, "Update Seqs", 
                      "Sequence contains a character that is not mentioned in score matrix\n");
		     return -1;
		 }
		 for( col=0; col < chap_pam_size; col++){
		     if(toupper(seq2[i]) == toupper(chap_pam[col])){
			 col_pam = col;
			 break;
		     }
		 }
      
		 if(col_pam < 0){
		      verror(ERR_WARN, "Update Seqs", 
                      "Sequence contains a character that is not mentioned in score matrix\n");
		     return -1;
		 }
		 if(score_matrix[row_pam][col_pam] > cutoff_align3 && 
		    score_matrix[row_pam][col_pam] <= cutoff_align2 )
		     match_syms[i] = symbol_align3[0];
		 else if (score_matrix[row_pam][col_pam] > cutoff_align2 && 
			  score_matrix[row_pam][col_pam] <= cutoff_align1 )
		     match_syms[i] = symbol_align2[0];
		 else if (score_matrix[row_pam][col_pam] > cutoff_align1 )
		     match_syms[i] = symbol_align1[0];
		 else
		     match_syms[i] = ' ';
		 }
	 }else
	     if (type == DNA){
		 if (toupper(seq1[i]) == toupper(seq2[i])) {
		     match_syms[i] = '|';
		 }else 
		     match_syms[i] = ' ';  
	     }
     }
     for ( i=0,p1=pos1,p2=pos2;i<seq_len;i+=line_length) {
	 vmessage("        ");
	 vmessage("%4s","    ");
	 for (j=0;j<6 && p1<pos1+seq_len;j++,p1+=10) {
	     spads1=0;
	     for(l=0; l<10 && j*10+l+i < seq_len; l++){
		 if (seq1[j*10+l+i]=='.')
		     spads1++;
	     }
	     if ( seq1[p1-pos1]=='.'){
		 vmessage("%10c", '-');
		 p11=p11-spads1+10;
	     }
	     else{
		 vmessage("%10d",p11);
		 p11=p11-spads1+10;
	     }
	}
	vmessage("\n%20.16s %.*s\n                 ",
		 name1,
		 i+line_length < seq_len ? 60 : seq_len - i,
		 &seq1[i]);

	vmessage("%4s", "    ");
 
     	for (k=i;k<seq_len && k<i+line_length;k++) {
	    vmessage("%c",match_syms[k]);
	}
	vmessage("\n%20.16s %.*s\n        ",
		 name2,
		 i+line_length < seq_len ? 60 : seq_len - i,
		 &seq2[i]);

	vmessage("%4s","    ");

	for (j=0;j<6 && p2<pos2+seq_len;j++,p2+=10) {
	    spads2=0;

            for(l=0; l<10 && j*10+l+i < seq_len; l++){
		if(seq2[j*10+l+i]=='.')
		    spads2++;
	    } 
	    if (seq2[p2-pos2]=='.'){
	      vmessage("%10c", '-');
	      p22=p22-spads2+10;
	    }
	    else{
		vmessage("%10d",p22);
		p22=p22-spads2+10;
	    }
	}
	vmessage("\n\n");
     }
     free ( match_syms );
     return 0;
}

      















