#include <tcl.h>
#include <tk.h>
#include <itcl.h>
#include <itk.h>

#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dna_utils.h"
#include "misc.h"
#include "editor.h"
#include "tcl_utils.h"
#include "feature_colour.h"
#include "genetic_code.h"
#include "text_editor.h"
#include "editor.h"
#include "end_editor.h"
#include "renzyme_search.h"
#include "read_sequence.h"
#include "feature_table.h"
#include "graphic_editor.h"
#include "editor_reg.h"


static int *pos = NULL;
static int *score = NULL;

int GetSeqInfo (ClientData clientData,
		Tcl_Interp *interp, 
		int argc, 
		char **argv)
{
    int i, num_seq;
    EDITOR_RECORD *er;
    int editor_id, group_id, num_editor;
    
    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "editor_id\"", (char*)NULL);
        return TCL_ERROR;
    } 
    Tcl_ResetResult(interp);

    editor_id = atoi (argv[1]);
    num_editor = GetEditorNum ();
    if (editor_id + 1 > num_editor) {
      Tcl_AppendResult(interp, "wrong #editor number \"", 
		       (char*)NULL);
      return TCL_OK;
    }
    er = GetEditor (editor_id);
    group_id = 0; /*FIXME*/
    num_seq = er->seq_group[group_id]->nmembers;
    for (i = 0; i <= num_seq; i++) {
	char *seq, *name;
	Tcl_DString dstr;
	Tcl_DStringInit(&dstr);   
	seq = GetEditorSeq (er, group_id, i);
	name = GetEditorSeqName (er, group_id, i);
	Tcl_DStringAppendElement(&dstr, name);
	Tcl_DStringAppendElement(&dstr, seq);
	Tcl_AppendElement (interp, Tcl_DStringValue(&dstr));
	Tcl_DStringFree(&dstr);
    }
    return TCL_OK;
}

int GetSequence (ClientData clientData,
		Tcl_Interp *interp, 
		int argc, 
		char **argv)
{
    char *seq;
    int i, num_seq;
    EDITOR_RECORD *er;
    int editor_id, group_id;

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "editor_id\"", (char*)NULL);
        return TCL_ERROR;
    } 
    Tcl_ResetResult(interp);

    editor_id = atoi (argv[1]);
    er = editor_records->editor_record[editor_id];
    group_id = 0; /*FIXME*/
    num_seq = er->seq_group[group_id]->nmembers;
    for (i = 0; i <= num_seq; i++) {
      seq = er->seq_group[group_id]->members[i]->data->sequence->seq;
      Tcl_AppendElement (interp, seq);
    }
    return TCL_OK;
}

int GetNumSeq (ClientData clientData,
		Tcl_Interp *interp, 
		int argc, 
		char **argv)
{
    EDITOR_RECORD *er;
    int num_seq;
    int editor_id, group_id;
    
    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "editor_id\"",(char*)NULL);
        return TCL_ERROR;
    }

    editor_id = atoi (argv[1]);
    er = editor_records->editor_record[editor_id];
    group_id = 0; /*FIXME*/

    Tcl_ResetResult(interp);
    num_seq = er->seq_group[group_id]->nmembers;
    vTcl_SetResult(interp, "%d", num_seq);

    return TCL_OK;
}

int GetBufferInfo (ClientData clientData,
		Tcl_Interp *interp, 
		int argc, 
		char **argv)
{
    Tcl_ResetResult(interp);

    if (editor_records->clipboard == NULL)
	vTcl_SetResult(interp, "%d", 0);
    else 
	vTcl_SetResult(interp, "%d", 1);

    return TCL_OK;
}


int GetBaseNum (ClientData clientData,
		Tcl_Interp *interp, 
		int argc, 
		char **argv)
{
    
    int start, end, i;
    int flag, p0;
    char tmp[10];
    char buf[1024] = { 0 };
    char buf1[1024] = { 0 };
  
    if (argc != 3) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "start_pos \"", "end_pos\"", (char*)NULL);
        return TCL_ERROR;
    }
    
    flag = 0;
    start = atoi (argv[1]);
    end = atoi (argv[2]) + 10;
 
    p0 = start - (start/10)*10;
    if (p0 == 0 && start ) p0 = 10;
    for (i = start; i <= end; i++) {
      if ( (i/10) * 10 == i && i != 0) {
	sprintf(tmp, "%10d", i);
	    strcat(buf, tmp);   
      }
    }
    strcpy (buf1, &buf[p0]);
    Tcl_ResetResult(interp);
    vTcl_SetResult(interp, "%s", buf1);
    return TCL_OK;
}

void seqed_string_search_free(void)
{
    xfree(pos);
    xfree(score);
    pos = NULL;
    score = NULL;
}

int string_search_in_sequence(EDITOR_RECORD *er, int group_id, int seq_id, char *string, int direction, 
			      int strand, int use_iub_code, double per_match)
{
    int min_match, max_matches ;
    char *string_match;
    int string_length;
    char *sequence, *seq;
    int n_matches;
    int seq_len;
    
    if (seq_id > er->seq_group[group_id]->nmembers) return -1;
    string_length = strlen(string);
    min_match = ceil(string_length * per_match / 100);
    max_matches = er->seq_group[group_id]->members[seq_id]->data->sequence->length;
    seq_len = er->seq_group[group_id]->members[seq_id]->data->sequence->length;

    /* if previously allocated pos and score, then need to free them here */
    if (pos != NULL) {
	seqed_string_search_free();
    } 
    if (NULL == (pos = (int *)xmalloc((max_matches + 1) * sizeof(int))))
	goto error;
    if (NULL == (score = (int *)xmalloc((max_matches + 1) * sizeof(int))))
	goto error;
    if (NULL == (string_match = (char *)xmalloc((string_length + 1) * sizeof(char))))
	goto error;
    if (NULL == (sequence = (char *)xmalloc((seq_len + string_length + 1) * sizeof(char))))
	goto error;

    seq = er->seq_group[group_id]->members[seq_id]->data->sequence->seq;
    strcpy (string_match, string);
    string_match[string_length] = 0;

    strcpy (sequence, er->seq_group[group_id]->members[seq_id]->data->sequence->seq);  
    memmove (&sequence[seq_len], &seq[0], string_length - 1);
    sequence[seq_len + string_length - 1] = 0;

    /* reverse & complement to search from 5' to 3' on complementary strand */   
    if (strand == 1) {
	complement_seq(string_match, string_length);
    }   
    n_matches = iubc_inexact_match(sequence, seq_len+string_length-1, string_match, string_length,
				   min_match, use_iub_code, pos, score, max_matches);
    /*printf("n_matches=%d\n", n_matches);*/
                      ;
    xfree (string_match);
    xfree (sequence);
    return n_matches;
 error:
    if (string_match) xfree (string_match);
    if (sequence) xfree (sequence);
    return -1;
}

int ExecStringSearch (ClientData clientData,
		      Tcl_Interp *interp, 
		      int argc, 
		      char **argv)
{
    int direction, strand, sa;
    int num_seq, k;
    float mpm; 
    char *string, *str;
    int i, j, num_match;
    char buf[1024];
    EDITOR_RECORD *er;
    int editor_id, num_group;
  
    if (argc != 7) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "direction\"", "strand\"", "search algorithm\"",
			 "Minimum percent match\"","string\"","editor_id\"",(char*)NULL);
        return TCL_ERROR;
    }

    direction = atoi (argv[1]);
    strand = atoi (argv[2]);
    sa = atoi (argv[3]);
    mpm = atof (argv[4]);
    str = argv[5];
    editor_id = atoi (argv[6]);
    er = editor_records->editor_record[editor_id];
    num_group = er->num_group;

    if ( ( NULL == (string = (char * ) xmalloc ((strlen(str) + 1)*sizeof(char) ))))
    return TCL_OK;

    strcpy (string, str);
    string[strlen(string)] = 0;
   
   if (direction == 0) {
       for (j = 0; j < num_group; j++) {
	   num_seq = er->seq_group[j]->nmembers;
	   for (k = 0; k <= num_seq; k++) {
	       num_match = string_search_in_sequence (er, j, k, string, direction, strand, sa, mpm);
	       if (num_match <= 0) {
		   /*printf("No match has been found\n");*/
		   return TCL_OK;
	       }	
	       for (i = 0; i < num_match; i++) {
		   sprintf(buf, "%d %d %d", k, pos[i], score[i]);
		   Tcl_AppendElement (interp, buf);
	       }
	   }
       }
   }
   	
   if (direction == 1) {
       for (j = 0; j < num_group; j++) {
	   num_seq = er->seq_group[j]->nmembers;
	   for (k = num_seq; k >= 0; k--) {
	       num_match = string_search_in_sequence (er, j, k, string, direction, strand, sa, mpm);
	       if (num_match <= 0) {
		   /*printf("No match has been found\n");*/
		   return TCL_OK;
	       }	 
	       for (i = num_match - 1; i >= 0; i--) {
		   sprintf(buf, "%d %d %d", k, pos[i], score[i]);
		   Tcl_AppendElement (interp, buf);
	       }
	   }
       }
   }

   return TCL_OK;
}

int check_over_lapping (ft_entry *e, int start, int end) {

    ft_range *r;
     
    for (r = e->range; r; r = r->next) {
	if (r->left->min <= end && (r->right && (r->right->min >= start) ) )
	    return 1;
    }
    if ((e->range->right && (e->range->right->min <= start)) && e->range->next) {
	for (r = e->range; r; r = r->next) {
	    if ((r->next == NULL) && r->right->min >= end) return 1;
	}
    }
    return 0;
}

char *get_feat_colour (char *type) {

    int i;

    if (!fcol_db) {
	get_fcol_types ();
    }
    for (i = 0; i < fcol_db_count; i++) {
	if (!strcmp (type, fcol_db[i].type)) {
	    return fcol_db[i].bg_colour;
	}
    }
    return NULL;
}

char *get_gene_name (ft_entry *e) {

    ft_value_element *ele;
    char *gname;
    int glen;

    ele = search_ft_qual_hash(e, "gene");
    if (ele && ele->value) {
	glen = strlen (ele->value);
	if (NULL == (gname = (char *)xmalloc(( glen + 1)*sizeof(char))))
	    return NULL;	
	strncpy (gname, &ele->value[0], glen);
	gname[glen] = 0;
    } else {
	glen = 7;
	if (NULL == (gname = (char *)xmalloc(( glen + 1)*sizeof(char))))
	    return NULL;	
	strcpy (gname, "unknown");
	gname[glen] = 0;
    } 
    return gname;
}

char *reformf(char *tr, int length, int codon_start) {

    char *tr_new;
    int i, len, str_len;
    int j = 0;
    int pos1 = codon_start + 1;

    str_len = strlen (tr);
    if (NULL == (tr_new = (char *)xmalloc ((length+1)*sizeof(char))))
	    return NULL;
    memset(tr_new, ' ', length);

    for (i = 0; i < str_len; i++) {
	if (isgraph (tr[i]) ) {
	    tr_new[pos1 + 3*j] = tr[i];
	    j++;
	} else {
	    len = len - 3;
	}
    }
    tr_new[length] = 0;
    return tr_new;
}

char *TranslateSubseq(char *subseq, int rf, int strand) {
    
    int i, cnt = 0;
    char *prot_seq;

    int length = strlen(subseq);
  
    if (NULL == (prot_seq = (char *)xmalloc(((length/3)+3) * sizeof(char))))
	return NULL;

    for (i = rf; i < length-2; i+=3) {
	if (strand == -1) {
	    prot_seq[cnt++] = codon_to_cacid1(&subseq[i]);
	} else if (strand == 1) {
	    prot_seq[cnt++] = codon_to_acid1(&subseq[i]);
	}
    }
    prot_seq[cnt] = '\0';
    return (prot_seq);
}


char *get_CDS_translation (ft_entry *e, char *seq, int length, int start, int end) {

    ft_value_element *ele;
    ft_range *r;
    char *buf, *prot_seq;
    int left, start_pos, end_pos;
    int transl_table_number, codon_start = 1;
    char *sub_seq = NULL;
   
    /* Look for /condon_start */
    ele = search_ft_qual_hash(e, "codon_start");
    if (ele && ele->value)
	codon_start = atoi(ele->value);
    /* Look for /transl_table */
    transl_table_number = 1; 		
    ele = search_ft_qual_hash(e, "transl_table");
    if (ele && ele->value)
	transl_table_number = atoi(ele->value);
    if (load_genetic_code_number(transl_table_number) == -1) {
	printf(" Failed to load code %d; using standard code",transl_table_number);
	load_genetic_code_number(1);
    }

    left = e->range->left->min;
    for (r = e->range; r; r = r->next) {
      char *range_seq;
      int rlen;
      start_pos = r->left->min;
      if (r->right) {
	end_pos = r->right->min;
      } else end_pos = start_pos;

      rlen = end_pos - start_pos + 1;
      if (NULL == (range_seq = (char*)xmalloc((rlen + 1)*sizeof(char))))
	goto error;
      if (sub_seq == NULL) {
	if (NULL == (sub_seq = (char*)xmalloc((rlen + 1)*sizeof(char))))
	  goto error;
	sub_seq[0] = '\0';
      } else {
	sub_seq = (char *)xrealloc (sub_seq, strlen(sub_seq) + 1 + rlen + 1); 
      }
      /*extract corresponding sub_sequence and to form coding sequence */  
      strncpy (range_seq, &seq[start_pos-1], rlen);
      range_seq[rlen] = 0;
	if (r->complemented)
	    (void) complement_seq (range_seq, rlen);
	strcat (sub_seq, range_seq);
	xfree (range_seq);
    }	
    prot_seq = TranslateSubseq (sub_seq, codon_start - 1, 1);
    buf = reformf (prot_seq, length, codon_start - 1);
    xfree (prot_seq);
    xfree (sub_seq);
    return buf;
 error:
    if (buf) xfree (buf);
    if (sub_seq) xfree (sub_seq);
    if (prot_seq) xfree (prot_seq);
    return NULL;
}
    



int GetFeats (ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {
   
    EDITOR_RECORD *er;
    FEATURE_TABLE *ft;    
    ft_entry *e;
    ft_range *r;   
    int i;
    char buf[1024];
    int length;
    char *seq, *prot;
    int editor_id, group_id, member_id, start, end;
    int overlapping ;
    int left, ll, rr;
    char *gname = NULL;
    char *colour;
    int next_range_start, flag;
        
    if (argc != 5) {
        Tcl_AppendResult(interp, "wrong # args: should be \"", argv[0],"editor_id\"", 
			"member_id\"","start\"","end\"",(char*)NULL);
        return TCL_ERROR;
    }

    editor_id = atoi (argv[1]);    
    member_id = atoi (argv[2]);
    start = atoi (argv[3]) + 1;
    end = atoi (argv[4]);
    group_id = 0; /*FIXME*/
    
    er = GetEditor (editor_id);
    ft = GetEditorSeqFt (er, group_id, member_id);
    length = GetEditorSeqLength (er, group_id, member_id);
    seq = GetEditorSeq (er, group_id, member_id);

    for (i = 0; i < ft->num_entry; i++) {
	e = ft->entry[i];
	overlapping = check_over_lapping (e, start, end);
	ll = 0;
	rr = 0;
        if (overlapping && strcmp (e->type, "source")) {
	  Tcl_DString dstr;
	  Tcl_DStringInit(&dstr);
	    sprintf(buf, "%s", e->type);
	    Tcl_AppendElement (interp, buf);

	    /* getfeature colour */
	    colour = get_feat_colour (e->type);	    
	    if (colour) {
		sprintf(buf, "%s", colour);
		Tcl_AppendElement (interp, buf);
	    } else {
		sprintf(buf, "%s", "red");
		Tcl_AppendElement (interp, buf);
	    }

	    /* Look for /gene qualifier */
	    gname = get_gene_name (e);
	    sprintf(buf, "%s", gname);
	    Tcl_AppendElement (interp, gname);

	    /* get CDS onthefly translation */
	    if (!strcmp (e->type, "CDS")) {	
		prot = get_CDS_translation (e, seq, length, start, end);
	    }
	    left = e->range->left->min;
	    next_range_start = 0;
	    for (r = e->range; r; r = r->next) {
		flag = 0;
		if (r->left->min <= start && r->right->min >= end) {
		    ll = 0;
		    rr = end - start;
		    flag = 1;	    
		}
		if (r->left->min < end && r->left->min >= start && r->right->min > end) {
		    ll = r->left->min - start;
		    rr = end - start;
		    flag = 1;  
		} 
		if (r->left->min >= start && r->right->min < end) {
		    ll = r->left->min - start;
		    rr = r->right->min - start;
		    flag = 1;
		} 
		if (r->left->min <= start && r->right->min < end && r->right->min > start) {
		    ll = 0;
		    rr = r->right->min - start;
		    flag = 1;	      
		}
		if (flag && !strcmp (e->type, "CDS")) {
		  Tcl_DString dstr1;
		  char *tr;
		  int rl = rr - ll + 1;
		  Tcl_DStringInit(&dstr1);
		  if (NULL == (tr = (char *)xmalloc(( rl+ 1)*sizeof(char))))
		    goto error;
		  if (start <= r->left->min) {
		    strncpy (&tr[0], &prot[next_range_start], rl);
		  }
		  else {
		    strncpy (&tr[0], &prot[start-left-next_range_start], rl);
		  }
		  tr[rl] = 0; 
		  vTcl_DStringAppendElement(&dstr1, "%s", tr);
		  sprintf(buf, "%d %d %s", ll, rr, Tcl_DStringValue(&dstr1));
		  vTcl_DStringAppendElement(&dstr, buf);
		  Tcl_DStringFree(&dstr1);
		  xfree (tr);
		 
		} else if (flag) {
		  sprintf(buf, "%d %d", ll, rr);
		  vTcl_DStringAppendElement(&dstr, buf);
		}
		 next_range_start += (r->right->min - r->left->min + 1);
	    }
	    Tcl_AppendElement (interp, Tcl_DStringValue(&dstr));
	    Tcl_DStringFree(&dstr);   
	}
    }
    /*if (prot) xfree (prot);*/
    return TCL_OK;    
    error:
    if (prot) xfree (prot);
    /*if (tr) xfree (tr);*/
    return TCL_OK;
}

int GetTrans (ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {
   
    EDITOR_RECORD *er; 
    char *buf, *prot;
    int i, length;
    char *seq;
    char *buf1, *buf2;
    int editor_id, group_id, member_id, start, end;
    int strand, seq_type;
    char **str = NULL;
    int num_str;
    char *sequence;

            
    if (argc != 6) {
        Tcl_AppendResult(interp, "wrong # args: should be \"", argv[0],"editor_id\"", 
			"member_id\"","start\"","end\"",(char*)NULL);
        return TCL_ERROR;
    }
    editor_id = atoi (argv[1]);    
    member_id = atoi (argv[2]);
    start = atoi (argv[3]) + 1;
    end = atoi (argv[4]);

    er = editor_records->editor_record[editor_id];
    group_id = 0;
    length = er->seq_group[group_id]->members[member_id]->data->sequence->length;
    seq_type = er->seq_group[group_id]->members[member_id]->data->sequence->type;
    /*seq_type = 1;*/
    seq_type = 2;
    if (start >= length) {
	return TCL_OK;
    }
    seq = er->seq_group[group_id]->members[member_id]->data->sequence->seq;
    load_genetic_code_number(1);
    
    if (NULL == (sequence = (char *)xmalloc((length + 4)*sizeof(char))))
	goto err;
    if (NULL == (buf1 = (char *)xmalloc((end - start + 2)*sizeof(char))))
	goto err;
    if (NULL == (buf2 = (char *)xmalloc((end - start + 2)*sizeof(char))))
	goto err;     
    if (Tcl_SplitList(interp, argv[5], &num_str, &str) != TCL_OK)
	     return TCL_ERROR;
   
    strcpy (&sequence[1], seq);
    memmove(&sequence[0], &seq[length - 1], 1);
    sequence[length+1] = 0;
    if (seq_type == 2) {
	memmove(&sequence[length + 1], &seq[0], 2);
	sequence[length+3] = 0;    
    }    
   
    for (i = 0; i < num_str; i++) {
	if (atoi (str[i]) != 0 && i <= 2 ) {
	    strand = 1;	    
	    if (seq_type == 2 && i == 2 ) {
		prot = TranslateSubseq (sequence, 0, strand);
		buf = reformf (prot, length, -1);
	    } else {
		prot = TranslateSubseq (&sequence[1], i, strand);
		buf = reformf (prot, length, i);
	    }
	    buf[strlen(buf)] = 0;
	    strncpy (buf1, &buf[start - 1], end - start + 1);
	    buf1[end-start+1] = 0;
	    sprintf(buf2, "Frame %d+", i+1);
	    Tcl_AppendElement (interp, buf2);  
	    sprintf(buf2, "%s", buf1);
	  
	    Tcl_AppendElement (interp, buf2);
	}
	if (atoi (str[i]) != 0 && i > 2 ) {  
	    strand = -1;
	    if (seq_type == 2 && i == 5 ) {
		prot = TranslateSubseq (sequence, 0, strand);
		buf = reformf (prot, length, -1);
	    } else {
		prot = TranslateSubseq (&sequence[1], i-3, strand);
		buf = reformf (prot, length, i-3);
	    }
	    buf[strlen(buf)] = 0;
	    strncpy (buf1, &buf[start-1], end-start+1);
	    buf1[end-start+1] = 0;
	    sprintf(buf2, "Frame %d-", i-2);
	    Tcl_AppendElement (interp, buf2);  
	    sprintf(buf2, "%s", buf1);
	    Tcl_AppendElement (interp, buf2);
	}
    }
    xfree (buf1);
    xfree (buf2);
    xfree (sequence);
    ckfree((char*)str);
    return TCL_OK;
 err:
    if (buf1) xfree (buf1);
    if (buf2) xfree (buf2);
    if (sequence) xfree (sequence);
    if(str) ckfree((char*)str);
    return TCL_OK;
}
  
int FindFeatKey (ClientData clientData,
		  Tcl_Interp *interp, 
		  int argc, 
		  char **argv)
{
    EDITOR_RECORD *er;
    FEATURE_TABLE *ft;
    ft_entry *e;
    char buf[1024];
    char *key;
    int i, j, k;
    int num_seq, num_entry, num_group;
    int editor_id;

    if (argc != 3) {
        Tcl_AppendResult(interp, "wrong # args: should be \"", argv[0], "ed_id\"","key\"",(char*)NULL);
        return TCL_ERROR;
    }
    editor_id = atoi (argv[1]);
    er = editor_records->editor_record[editor_id];    
    num_group = er->num_group;
    key = argv[2];
    
    for (k = 0; k < num_group; k++) {
	num_seq = er->seq_group[k]->nmembers;
	for (i = 1; i <= num_seq; i++) { /* FIXME: cons? */
	    ft = er->seq_group[k]->members[i]->data->sequence->feature_table;
	    num_entry = ft->num_entry;
	    for (j = 0; j < num_entry; j++) {
		e = ft->entry[j];
		if (!strcmp (e->type, key)) {
		    sprintf(buf, "%d %d", i, e->range->left->min);
		    Tcl_AppendElement (interp, buf);
		}
	    }
	}
    }
    return TCL_OK;
}

int FindFeatQual (ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {
    
    EDITOR_RECORD *er;
    FEATURE_TABLE *ft;
    ft_entry *e;
    char buf[1024];
    char *qual;
    int i, j, k;
    int num_seq, num_group, num_entry;
    int editor_id;

    if (argc != 3) {
        Tcl_AppendResult(interp, "wrong # args: should be \"", argv[0], "ed_id\"","key\"",(char*)NULL);
        return TCL_ERROR;
    }
    editor_id = atoi (argv[1]);
    er = editor_records->editor_record[editor_id];
    num_group = er->num_group;
    qual = argv[2];
    for (k = 0; k < num_group; k++) {
	num_seq = er->seq_group[k]->nmembers;
	for (i = 1; i <= num_seq; i++) { /* FIXME: cons? */
	    ft = er->seq_group[k]->members[i]->data->sequence->feature_table;
	    num_entry = ft->num_entry;
	    for (j = 0; j < num_entry; j++) {
		ft_value_element *ele;
		e = ft->entry[j];	  
		ele = search_ft_qual_hash(e, qual);
		if (ele && ele->value) {
		    sprintf(buf, "%d %d", i, e->range->left->min);
		    Tcl_AppendElement (interp, buf);
		}
	    }
	}
    }
    return TCL_OK;
}

int GetUndo (ClientData clientData,
		Tcl_Interp *interp, 
		int argc, 
		char **argv)
{
  int undo;
  EDITOR_RECORD *er;
  int editor_id;
 
  if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "editor_id\"", (char*)NULL);
        return TCL_ERROR;
    }

  Tcl_ResetResult(interp);

  editor_id = atoi (argv[1]);
  er = GetEditor (editor_id);

  if (er->edits->head) undo = 0;
  else undo = 1; /* no more undo */
    
  vTcl_SetResult(interp, "%d", undo);
  return TCL_OK;
}

int remove_edit(EDITS *edits) {
  
    EDIT *tmp;

    if(edits->head == NULL) return 1; /* special return value?? */
    tmp = edits->head;
    edits->head = edits->head->next;
    if(tmp) xfree(tmp);
    return 0;
}

int undo_edit (EDITOR_RECORD *editor_record, int *seq_id, int *pos, int *len) {

  int err = 0;
  EDIT *edit = NULL;
  int add_to_undo;
  int i, num_seq, group_id;
      
  add_to_undo = 0;
  edit = editor_record->edits->head;
  group_id = 0;
  num_seq = editor_record->seq_group[group_id]->nmembers;

  if ( NULL == edit ) {
    err = 1;
    goto bail_out;
  }

  if (edit->operation == INSERT_FRAGMENT ) {
      if (edit->consensus == 1) {
	 for (i = 0; i <= num_seq; i++) {
            *seq_id = edit->seq_id;
	    *pos = edit->position;
	    *len = -(strlen (edit->sequence->seq));
	    err = delete_fragment_in_editor (editor_record, edit, add_to_undo);
	    if (err == -1) goto bail_out;
	    err = remove_edit(editor_record->edits);
	    if (err) goto bail_out;
	    if (editor_record->edits->head != NULL) edit = editor_record->edits->head;
	 }
      } else {
          *seq_id = edit->seq_id;
	  *pos = edit->position;
	  *len = -(strlen (edit->sequence->seq));
	  err = delete_fragment_in_editor (editor_record, edit, add_to_undo);
	  if (err == -1) goto bail_out;
	  err = remove_edit(editor_record->edits);
	  if (err) goto bail_out;
	  err = update_consensus (editor_record->seq_group[group_id], 1);
      }     
  }
  if ( edit->operation == DELETE_FRAGMENT ) {
      if (edit->consensus == 1) {
	 for (i = 0; i <= num_seq; i++) {
             *seq_id = edit->seq_id;
	     *pos = edit->position;
	     *len = strlen (edit->sequence->seq);
	     err = insert_fragment_in_editor (editor_record, edit, add_to_undo);
	     if (err == -1) goto bail_out;
	     err = remove_edit(editor_record->edits);
	     if (err) goto bail_out;
	     if (editor_record->edits->head != NULL) edit = editor_record->edits->head; 
	 } 
      } else {  
          *seq_id = edit->seq_id;
	  *pos = edit->position;
	  *len = strlen (edit->sequence->seq);
	  err = insert_fragment_in_editor (editor_record, edit, add_to_undo);
	  if (err == -1) goto bail_out;
	  err = update_consensus (editor_record->seq_group[group_id], 1);
	  err = remove_edit(editor_record->edits);
	  if (err) goto bail_out;
      }
  }
    
  if ( edit->operation == CHANGE_FRAGMENT ) {
    if (edit->consensus == 1) {
	 for (i = 0; i <= num_seq; i++) {

	     printf ("undo:member_id=%i\n", i);
	     
	     *seq_id = edit->seq_id;
	     *pos = edit->position;
	     *len = strlen (edit->sequence->seq);
	     /* delete */
	     err = delete_fragment_in_editor (editor_record, edit, add_to_undo);
	     if (err == -1) goto bail_out;
	     err = remove_edit(editor_record->edits);
	     if (err) goto bail_out;
	     if (editor_record->edits->head != NULL) edit = editor_record->edits->head;
	     /*insert*/
	     err = insert_fragment_in_editor (editor_record, edit, add_to_undo);
	     err = remove_edit(editor_record->edits);
	     if (err) goto bail_out;
	     if (editor_record->edits->head != NULL) edit = editor_record->edits->head;
	 } 
    } else {
	*seq_id = edit->seq_id;
	*pos = edit->position;
	*len = strlen (edit->string);
	err = delete_fragment_in_editor (editor_record, edit, add_to_undo);
	err = remove_edit(editor_record->edits);
	if (err) goto bail_out;
	if (editor_record->edits->head != NULL) edit = editor_record->edits->head;
	/*delete*/
	err = insert_fragment_in_editor (editor_record, edit, add_to_undo);
	err = remove_edit(editor_record->edits);
	if (err) goto bail_out;
	/*insert*/
	err = update_consensus (editor_record->seq_group[group_id], 1);
	err = remove_edit(editor_record->edits);
	if (err) goto bail_out;
      }    
  }
  if (edit->operation == INSERT_FRAGMENT_COMPLEMENT) {
      *seq_id = edit->seq_id;
      *pos = edit->position;
      *len = -(strlen (edit->sequence->seq));
      err = delete_fragment_in_editor (editor_record, edit, add_to_undo);
      if (err == -1) goto bail_out;
      err = remove_edit(editor_record->edits);
      if (err) goto bail_out;
      err = update_consensus (editor_record->seq_group[group_id], 1);
      /* close complement fragment window */
      if (editor_record->edits->head != NULL) edit = editor_record->edits->head;
      editor_complement_shutdown (editor_record, edit, add_to_undo);
      err = remove_edit(editor_record->edits);
      if (err) goto bail_out;
  }
  if (edit->operation == TRIM) {
      *seq_id = edit->seq_id;
      *pos = edit->position;
      *len = strlen (edit->sequence->seq);
      err = un_trim_sequence (editor_record, edit, add_to_undo);   
      if (err) goto bail_out;
  }
  if (edit->operation == FILL) {
      *seq_id = edit->seq_id;
      *pos = edit->position;
      *len = strlen (edit->sequence->seq);
      err = un_fill_sequence (editor_record, edit, add_to_undo);   
      if (err) goto bail_out;
  }

  /*if (edit->operation == ADD_SEQUENCE) {
    err = remove_sequence_in_editor (editor_record, edit, add_to_undo);
    if (err) goto bail_out;
  }
  if (edit->operation == REMOVE_SEQUENCE) {
    err = add_sequence_in_editor (editor_record, edit, add_to_undo);
    if (err) goto bail_out;
  }
  if (edit->operation == ADD_FEATURE) {
    err = remove_feature_in_editor (editor_record, edit, add_to_undo);
    if (err) goto bail_out;
  }
  if (edit->operation == REMOVE_FEATURE) {
    err = add_feature_in_editor (editor_record, edit, add_to_undo);
    if (err) goto bail_out;
    }*/
  
  /*err = remove_edit(editor_record->edits);
    return err;*/
 bail_out:
  return err;
}



int ExecUndo (ClientData clientData,
		Tcl_Interp *interp, 
		int argc, 
		char **argv)
{ 
    int seq_id = 0, pos = 0, len = 0;
    int err;
    char buf[1024];
    EDITOR_RECORD *er;
    int editor_id;

    if (argc != 2) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], "editor_id\"", (char*)NULL);
	return TCL_ERROR;
    }
    
    editor_id = atoi (argv[1]);
    er = GetEditor (editor_id);

    err = undo_edit (er, &seq_id, &pos, &len);
    /*printf(" seq_id=%d pos=%d  len=%d\n", seq_id, pos, len);*/
     
    sprintf(buf, "%d %d %d", seq_id, pos, len);
    vTcl_SetResult (interp, "%s", buf);

    return TCL_OK;
}

int exec_undo (int ed_id, int *member_id, int *p) {
 
     EDITOR_RECORD *er;
    int seq_id = 0, pos = 0, len = 0;
    int err;
   
    er = GetEditor (ed_id);
    err = undo_edit (er, &seq_id, &pos, &len);
    *member_id = seq_id;
    *p = pos + 1;

    return 0;
}

int editor_save (int ed_id) {

    EDITOR_RECORD *er;
    SEQ_GROUP *sg;
    int num_group, num_member;
    int i, j;
    
    er = GetEditor (ed_id);
    num_group = er->num_group;
    
    for ( i = 0; i < num_group; i++) {
	sg = er->seq_group[i];
	num_member = sg->nmembers;
	for (j = 1; j <= num_member; j++) {
	    if (-1 == sequence_save (sg->members[j]->data->sequence) )
		return -1;
	}
    }
    return 0;
}

int SequenceRedisplay (ClientData clientData,
		       Tcl_Interp *interp, 
		       int argc, 
		       char **argv)
{
   
    char *operation, *frame_name;
    int ed_id, group_id = 0, member_id = 0, pos;
    int err;
     
    if (argc < 3) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], " operation\"" , "ed_id\"", (char*)NULL);
        return TCL_ERROR;
    } 
    operation = argv[1];    
    ed_id = atoi (argv[2]);

    if (!strcmp (operation, "insert")) {

	SEQUENCE *s;
	char *string;	
	if (argc != 7) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " editor_id\"", "member_id\"", 
			     "pos\"", "string\"",(char*)NULL);
	    return TCL_ERROR;
	}
	member_id = atoi (argv[3]);
	pos = atoi (argv[4]);
	string = argv[5];
	frame_name = strdup (argv[6]);
  
	if (string != NULL) s = init_sequence ();
	s->seq = strdup (string);
	err = insert_update_sequence (ed_id, group_id, member_id, pos, s);
	if (err != 0) goto bail;
	/* new cursor position after editing */
	pos = pos + strlen(string);
    }

    /* The difference between "insert" and "paste" is where is the 
       insertion string come from, paste may include FT*/ 
    if (!strcmp (operation, "paste")) {

	SEQUENCE *s;	
	if (argc != 6) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " editor_id\"", "member_id\"", 
			     "pos\"", (char*)NULL);
	    return TCL_ERROR;
	}
	member_id = atoi (argv[3]);
	pos = atoi (argv[4]);
	frame_name = strdup (argv[5]);
       
	s = copy_sequence (editor_records->clipboard);
	if (s == NULL) goto bail;
	err = insert_update_sequence (ed_id, group_id, member_id, pos, s);
	if (err != 0) goto bail;
	/* new position after editing */
	pos = pos + strlen(s->seq);
    }

    if (!strcmp (operation, "delete")) {
	SEQUENCE *s;
	char *str;
	if (argc != 7) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " editor_id\"", "member_id\"", 
			     "pos\"", "string\"",(char*)NULL);
	    return TCL_ERROR;
	}
	member_id = atoi (argv[3]);
	pos = atoi (argv[4]);
	str = argv[5];
	frame_name = strdup (argv[6]);

	if (str != NULL) s = init_sequence ();
	s->seq = strdup (str);
	s->length = strlen (s->seq);
	if (s->length > 1) {
	    FEATURE_TABLE *ft = NULL;
	    SEQUENCE *ss;
	    ss = editor_records->editor_record[ed_id]->seq_group[group_id]->members[member_id]->data->sequence;
	    ft = get_feature_table (ss, s->length, pos);
	    if (ft) {
		change_feature_table_location ( ft, -(pos) );
		s->feature_table = ft;
	    }
	}
	err = delete_update_sequence (ed_id, group_id, member_id, pos, s);
	if (err != 0) goto bail;
	pos++;
    }

    /* The difference between "delete" and "cut" is cut will save cut sequence 
       to editor buffer(clipboard) */ 
    if (!strcmp (operation, "cut")) {
	SEQUENCE *s;
	char *str;
	if (argc != 7) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " editor_id\"", "member_id\"", 
			     "pos\"", "string\"",(char*)NULL);
	    return TCL_ERROR;
	}
	member_id = atoi (argv[3]);
	pos = atoi (argv[4]);
	str = argv[5];
	frame_name = strdup (argv[6]);
	if (str != NULL) s = init_sequence ();
	s->seq = strdup (str);
	err = copy_update_buffer (ed_id, group_id, member_id, pos, s);
	if (err != 0) goto bail;
	
	if (s->length > 1) {
	    FEATURE_TABLE *ft = NULL;
	    SEQUENCE *ss;
	    ss = editor_records->editor_record[ed_id]->seq_group[group_id]->members[member_id]->data->sequence;
	    ft = get_feature_table (ss, s->length, pos);
	    if (ft) {
		change_feature_table_location ( ft, -(pos) );
		s->feature_table = ft;
	    }
	}
	err = delete_update_sequence (ed_id, group_id, member_id, pos, s);
	if (err != 0) goto bail;
	pos++;
    }
    if (!strcmp (operation, "replace")) {
	char *str_replace, *str_with;
	SEQUENCE *sr, *sw;

	if (argc != 8) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " editor_id\"", "member_id\"", "pos\"", "etc\"",(char*)NULL);
	    return TCL_ERROR;
	}
	member_id = atoi (argv[3]);
	pos = atoi (argv[4]);
	str_replace = argv[5];
	str_with = argv[6];
	frame_name = strdup (argv[7]);
	sr = init_sequence ();
	sr->seq = strdup (str_replace);
	sw = init_sequence ();
	sw->seq = strdup (str_with);
	err = replace_update_sequence (ed_id, group_id, member_id, pos, sr, sw);
	if (err != 0) goto bail;
	pos = pos + strlen(str_with) + 1;
    }

    if (!strcmp (operation, "copy")) {

	SEQUENCE *s;
	char *string;	
	if (argc != 7) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " editor_id\"", "member_id\"", 
			     "pos\"", "etc\"",(char*)NULL);
	    return TCL_ERROR;
	}
	member_id = atoi (argv[3]);
	pos = atoi (argv[4]);
	string = (argv[5]);
	string[strlen(string)] = 0;
	frame_name = strdup (argv[6]);

	if (string != NULL) s = init_sequence ();
	s->seq = strdup (string);
	err = copy_update_buffer (ed_id, group_id, member_id, pos, s);
	if (err != 0) goto bail;
	err = 1; /* copy does not need notify */
    }
    if (!strcmp (operation, "select")) {

	seq_reg_selected info;
        SELECTION *sel;
	int seq_id, seq_num;

	if (argc != 7) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " editor_id\"", "member_id\"", 
			     "pos\"", "etc\"",(char*)NULL);
	    return TCL_ERROR;
	}
	sel = init_selection ();
	sel->ed_id = ed_id;
	sel->member_id = atoi(argv[3]);
	sel->sel_first = atoi(argv[4]) + 1;
	sel->sel_last = atoi(argv[5]) + 1;
	if (editor_records->selection != NULL) 
	    free_selection (editor_records->selection);
	editor_records->selection = sel;

	/* selection notify */	
	info.job = SEQ_SELECTED;
        info.selection = sel;
	info.frame_name = strdup (argv[6]);
	seq_id = GetEditorSeqId (ed_id, group_id, sel->member_id);
     
	if (seq_id != -1) seq_num = GetEdenNum (seq_id);
	else goto bail;
	editor_notify (seq_num, (editor_reg_data *)&info);
	err = 1;
    }
   
    if (!strcmp (operation, "undo")) {

	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], "editor_id\"", "etc\"",(char*)NULL);
	    return TCL_ERROR;
	}
	
	exec_undo (ed_id, &member_id, &pos);
	frame_name = strdup (argv[3]);
	
	err = 0; 
    }
    
    if (!strcmp (operation, "icursor_move")) {
	
	seq_reg_cursor_notify info;
	seq_reg_changed sc;
	SELECTION *sel;
	int seq_id, seq_num;

	if (argc != 6) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], "ed_id\"", "cursorY\"", 
			     "cursorX\"", "win_name\"", (char*)NULL);
	    return TCL_ERROR;
	}

	member_id = atoi (argv[3]);
	if (member_id == 0) { /*cons*/
	    member_id = 1;
	}
	pos = atoi (argv[4]);
	seq_id = GetEditorSeqId (ed_id, group_id, member_id);
	seq_num = GetEdenNum (seq_id);

	info.cursor = get_editor_cursor (ed_id, member_id, argv[5]);
	if (info.cursor == NULL) goto bail;
	
	info.cursor->abspos = pos;
	info.cursor->posy = member_id;
	info.cursor->job = CURSOR_MOVE;
	info.cursor->sent_by = -1;
	sel = init_selection ();
	sel->ed_id = ed_id;
	sel->member_id = member_id;
	sel->sel_first = 0;
	sel->sel_last = 0;

	sc.job = SEQ_CHANGED;
	sc.selection = sel;
	info.job = SEQ_CURSOR_NOTIFY;
	info.selection = sel;
	/*editor_notify (seq_num, (editor_reg_data *)&sc);*/
	editor_notify (seq_num, (editor_reg_data *)&info);
	err = 1;
    }
    
    if (!strcmp (operation, "save")) {
	
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], "editor_id\"", (char*)NULL);
	    return TCL_ERROR;
	}
	
	err = editor_save (ed_id);
	if (err == -1) goto bail;
	err = 1;
    }
    if (!strcmp (operation, "exit")) {

	seq_reg_exit info;
	int seq_id, seq_num;

	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], "editor_id\"", (char*)NULL);
	    return TCL_ERROR;
	}
	frame_name = strdup (argv[3]);
	
	info.job = TEXT_EDITOR_QUIT;
	info.frame_name = strdup (frame_name);
	/* select first member in editor to get register number */
	seq_id = GetEditorSeqId (ed_id, group_id, 1);
	seq_num = GetEdenNum (seq_id);
	editor_notify (seq_num, (editor_reg_data *)&info);
        err = 1;
    }

    if (!err) {
	seq_reg_changed info;
	seq_reg_cursor_notify cn;
	SELECTION *sel;
	int seq_id, seq_num;
	
	sel = init_selection ();
	sel->ed_id = ed_id;
	sel->member_id = member_id;
	info.job = SEQ_CHANGED;
        info.selection = sel;
	seq_id = GetEditorSeqId (ed_id, group_id, sel->member_id);
	if (seq_id != -1) seq_num = GetEdenNum (seq_id);
	else goto bail;

	editor_notify (seq_num, (editor_reg_data *)&info);

	/* cursor notify */
	cn.cursor = get_editor_cursor (ed_id, member_id, frame_name);
	if (cn.cursor == NULL) goto bail;
	 
	cn.cursor->abspos = pos;
	cn.cursor->posy = member_id; 
	cn.cursor->job = CURSOR_MOVE;
	cn.cursor->sent_by = -1;
	cn.job = SEQ_CURSOR_NOTIFY;
	editor_notify (seq_num, (editor_reg_data *)&cn);
    }
    return TCL_OK;

 bail:
    Tcl_AppendResult(interp, "text editor operation error \"", 
		       (char*)NULL);
    return TCL_OK;    
}

int max_renzyme_cut_pos (RENZYMES *sel_rs) {

    int i, num_ren;
    int over_lap = 0;
    
    num_ren = sel_rs->used;
    for ( i = 0; i < num_ren; i++) {
	over_lap = MAX (over_lap, sel_rs->renzyme[i]->cut_pos_1);
    }
    return over_lap;
}

int save_sites_in_sequence (SEQUENCE *s, RENZYMES *sel_rs, R_MATCH *match, int total_matches) {

    int i, num_site = 0;
    SITES *ss;

    ss = s->sites;

    if (ss == NULL) {
	if (NULL == (ss = init_sites ()))
	    goto err;
    }
    num_site = ss->used;
 
    for (i = total_matches-1; i >= 0; i--) {
	SITE  *st;
	st = init_site();
	if (!st) goto err;
	st->pos1 = match[i].cut_pos1 - 1;	
	if (match[i].cut_pos2 != 0) 
	    st->pos2 = match[i].cut_pos2 - 1;

	if (st->pos1 > s->length)
	    st->pos1 = st->pos1 - s->length;
	if (st->pos2 > s->length)
	    st->pos2 = st->pos2 - s->length;
	st->name = strdup (sel_rs->renzyme[match[i].enz_name]->name);
	st->rec_seq = strdup (sel_rs->renzyme[match[i].enz_name]->rec_seq);
         
	if (num_site != 0) {
	    ss = realloc_sites (ss, num_site);
	}
	ss->site[num_site] = st;
	num_site++;
	ss->used = num_site;
    }
    s->sites = ss;
    
    return 0;
 err:
    return -1;
}


int RenzSearch (ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {

    EDITOR_RECORD *er;
    RENZYMES *sel_rs;
    SEQUENCE *seque;
    int num_sel;
    char **sel = NULL;
    int i, j;
    int ed_id, group_id, member_id;
    char *seq;
    int seq_len, subseq_len;
    int start, end;
    int seq_type;
    int total_matches;
    char buf[1024];
    char *subseq;
    int overlap, max_cut_pos, max_rec_length;
    int line_len, start_from, start_to, end_from;
   
    if (argc != 6) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "ed_id\"", "selected items\"", "start\"", "end\"",(char*)NULL);
        return TCL_ERROR;
    }    
    ed_id = atoi (argv[1]);
    /* create selecing Renzyme name array */
    if (Tcl_SplitList(interp, argv[2], &num_sel, &sel) != TCL_OK)
	goto err;
    sel_rs = get_selected_renzyme (num_sel, sel); 
    member_id = atoi (argv[3]);
    start = atoi (argv[4]);
    end = atoi (argv[5]);
    line_len = end - start;
    er = GetEditor (ed_id);
    group_id = 0;
    seque = er->seq_group[group_id]->members[member_id]->data->sequence;
    seq = seque->seq;
    seq_len = seque->length;
    seq_type = seque->type;
    /*seq_type = 1;*/ /* DNA linear */
    seq_type = 2; /* DNA circular */
    
    max_cut_pos = max_renzyme_cut_pos (sel_rs);
    max_rec_length = max_rec_seq_length (sel_rs);
    overlap =  max_cut_pos + max_rec_length;
    
    if (NULL == (subseq = (char *)malloc( (line_len + 2*overlap +1) * sizeof(char))))
	goto err; /* FIXME: do we need set max_line_length */
    
    if (start == 0 && er->seq_group[group_id]->members[member_id]->data->sequence->sites == NULL) {
	char *sequence;
	int i, total_match_all, sequence_len;
	
	if (NULL == (sequence = (char *)malloc( (seq_len + overlap +1) * sizeof(char))))
	    goto err;
	memmove (sequence, seq, seq_len);
	sequence_len = seq_len;
    
	if (seq_type == 2 || seq_type == 4) { 
	    memmove (&sequence[seq_len], &seq[0], overlap - 1);
	    sequence_len = seq_len + overlap-1;
	}
	sequence[sequence_len] = 0;
	
	for (i = 0; i < sel_rs->used; i++) {
	    R_MATCH *match_all;
	    if (NULL == (match_all = (R_MATCH*)xcalloc(MAXMATCHES, sizeof(R_MATCH))))
		goto err;
	    find_matches(sel_rs->renzyme[i], sequence, sequence_len, 
			 seq_type,i, &match_all, &total_match_all);
	    if (total_match_all != 0) {	
		save_sites_in_sequence (seque, sel_rs, match_all, total_match_all);
	    }
	    xfree (match_all);  
	}
	/*save_fragment_in_sequence (seque, sequence);*/
    }
    start_from = start - overlap;
    start_to = 0;
    end_from = start + line_len + overlap; 
    memset(subseq, 'N', line_len + 2*overlap);
    if (start_from < 1) {
	start_to = overlap - (start - 1);
	start_from = 1;
    } 
   
    if (end_from > seq_len) {	
	end_from = seq_len;
    }
    if (start_from >= seq_len)
	subseq = "";	
    else {
	memmove(&subseq[start_to], &seq[start_from-1], end_from - start_from + 1);
	subseq[line_len + 2*overlap] = 0;
    }
    if (seq_type == 2 || seq_type == 4) {
	start_from = start - overlap;
	end_from = start + line_len + overlap;
	if (start_from >= seq_len)
	    subseq = "";
	else {
	    if (start_from < 1) {  
		memmove(&subseq[0], &seq[seq_len - (overlap - start +1)], overlap-start + 1);
	    }
	    if (end_from > seq_len) {
		int end_to = overlap + seq_len - (start - 1);
		memmove(&subseq[end_to], &seq[0], end_from - seq_len - 1);
	    }
	}
    }	
    subseq_len = strlen(subseq); 
    if (subseq != "") {
	for (i = 0; i < sel_rs->used; i++) {
	    R_MATCH *match;
	    if (NULL == (match = (R_MATCH*)xcalloc(MAXMATCHES, sizeof(R_MATCH))))
		goto err;
	    find_matches(sel_rs->renzyme[i], subseq, subseq_len, seq_type, i, &match, &total_matches);
	    if (total_matches != 0) {
		sprintf(buf, "%s", sel_rs->renzyme[i]->name);
		Tcl_AppendElement(interp, buf);	
		for (j = total_matches-1; j >= 0; j--) {
		    match[j].cut_pos1 = match[j].cut_pos1 - 1 -overlap - 1;
		    if (match[j].cut_pos2 != 0)
		    match[j].cut_pos2 = match[j].cut_pos2 - 1 - overlap - 1;
		    if (match[j].cut_pos1 > (seq_len - start))
			match[j].cut_pos1 = match[j].cut_pos1 - seq_len ;
		    if (match[j].cut_pos2 > (seq_len - start))
			match[j].cut_pos2 = match[j].cut_pos2 - seq_len ;
		    if (match[j].cut_pos1 < 0 && match[j].cut_pos1 > start_from && start_from < 0)
			match[j].cut_pos1 = match[j].cut_pos1 + seq_len ;
		    if (match[j].cut_pos2 < 0 && match[j].cut_pos2 > start_from && start_from < 0)
			match[j].cut_pos2 = match[j].cut_pos2 + seq_len ;
		    /*printf ("cutpos=%d  %d\n", match[j].cut_pos1, match[j].cut_pos2);*/
		    sprintf(buf, "%d %d", match[j].cut_pos1, match[j].cut_pos2);
		    Tcl_AppendElement(interp, buf);		    	
		}
     

	    }
	    xfree (match);
	}
    }
    ckfree((char*)sel);
    return TCL_OK;
 err:
    if(sel)   ckfree((char*)sel);
    if (subseq) xfree(subseq);  
    return TCL_OK;
}

int TextEditor_Init (Tcl_Interp *interp) {

    if (Itcl_RegisterC(interp, "get_buffer_info", GetBufferInfo, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "get_seq_info", GetSeqInfo, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "get_sequence", GetSequence, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "get_num_seq", GetNumSeq, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "get_base_num", GetBaseNum, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "exec_string_search", ExecStringSearch, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "get_feature", GetFeats, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "get_translation", GetTrans, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "feat_key_found", FindFeatKey, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "feat_qual_found", FindFeatQual, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "get_fcol", tcl_get_fcol_array, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "get_undo", GetUndo, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "exec_undo", ExecUndo, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "seq_redisplay", SequenceRedisplay, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "renzyme_search", RenzSearch, NULL, NULL) != TCL_OK) {
      return TCL_ERROR;
    }
    return TCL_OK;      
}


