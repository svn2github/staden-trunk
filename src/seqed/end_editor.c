#include <tcl.h>
#include <tk.h>
#include <itcl.h>
#include <itk.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "misc.h"
#include "tcl_utils.h"
#include "renzyme_search.h"
#include "editor.h"
#include "dna_utils.h"
#include "feature_table.h"
#include "editor_reg.h"
#include "end_editor.h"
#include "text_editor.h"

#define END_LEN 6

int get_sequence_end (int seq_id, char **sl, char **sr, 
		      char **sl_c, char **sr_c, int *sl_len, int *sr_len) {
    
    int seq_num;
    SEQUENCE *s;
    char *ssl, *ssr, *ssl_c, *ssr_c;
    int cut_site;

    seq_num = GetEdenNum(seq_id);
    s = GetEdenSequence (seq_num);
    cut_site = s->cut_site_1->pos;

    if (s->cut_site_1 != NULL && s->cut_site_1->len != 0) {
      
	int l = END_LEN + ABS(s->cut_site_1->len);
	*sl_len = s->cut_site_1->len;
	if (s->cut_site_1->len > 0) {
	  
	    if(NULL == (ssl = (char *)xmalloc((l + 1) * sizeof(char))))
		goto err;
	    if(NULL == (ssl_c = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
		goto err;
	    strncpy (ssl, &s->seq[cut_site - l], l);
	    ssl[l] = 0;
	    *sl = ssl;
	    strncpy (ssl_c, &s->seq[cut_site - l], END_LEN);
	    ssl_c[END_LEN] = 0;
	    complement_dna (ssl_c, END_LEN);
	    *sl_c = ssl_c;
	} else {
	    
	    if(NULL == (ssl = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
		goto err;
	    if(NULL == (ssl_c = (char *)xmalloc((l + 1) * sizeof(char))))
		goto err;
	    strncpy (ssl, &s->seq[cut_site - END_LEN ], END_LEN);
	    ssl[END_LEN] = 0;
	    *sl = ssl;
	    strcpy (ssl_c, *sl);
	    strcat (ssl_c, s->cut_site_1->hang);
	    ssl_c[l] = 0;
	    complement_dna (ssl_c, l);
	    *sl_c = ssl_c;
	}
    } else {
	int ll = cut_site - END_LEN;
	*sl_len = 0;
	if(NULL == (ssl = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
	    goto err;
	if(NULL == (ssl_c = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
	    goto err;
	strncpy (ssl, &s->seq[ll], END_LEN);
	ssl[END_LEN] = 0;
	*sl = ssl;
	strncpy (ssl_c, &s->seq[ll], END_LEN);
	ssl_c[END_LEN] = 0;
	complement_dna (ssl_c, END_LEN);
	*sl_c = ssl_c;

    }
    if (s->cut_site_2 != NULL && s->cut_site_2->len != 0) {
	int l = END_LEN + ABS(s->cut_site_2->len);
	*sr_len = s->cut_site_2->len;
	if (s->cut_site_2->len > 0) {
	    
	    if(NULL == (ssr = (char *)xmalloc((l + 1) * sizeof(char))))
		goto err;
	    if(NULL == (ssr_c = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
		goto err;
	    strncpy (ssr, &s->seq[cut_site], l);
	    ssr[l] = 0;
	    *sr = ssr;
	  
	    strncpy (ssr_c, &s->seq[cut_site + s->cut_site_2->len], END_LEN);
	    ssr_c[END_LEN] = 0;
	    complement_dna (ssr_c, END_LEN);
	    *sr_c = ssr_c;
	} else {
	    if(NULL == (ssr = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
		goto err;
	    if(NULL == (ssr_c = (char *)xmalloc((l + 1) * sizeof(char))))
		goto err;
	    strncpy (ssr, &s->seq[cut_site], END_LEN);
	    ssr[END_LEN] = 0;
	    *sr = ssr;
	    strcpy (ssr_c, s->cut_site_2->hang);
	    strcat (ssr_c, ssr);
	    ssr_c[l] = 0;
	    complement_dna (ssr_c, l);
	    *sr_c = ssr_c;
	}
    } else {	
	*sr_len = 0;
	if(NULL == (ssr = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
	    goto err;
	if(NULL == (ssr_c = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
	    goto err;
	strncpy (ssr, &s->seq[cut_site], END_LEN);
	ssr[END_LEN] = 0;
	*sr = ssr;
	strncpy (ssr_c, &s->seq[cut_site], END_LEN);
	ssr_c[END_LEN] = 0;
	complement_dna (ssr_c, END_LEN);
	*sr_c = ssr_c;
    }
    
    return 0;
 err: 
    if (ssl) xfree (ssl);
    if (ssl_c) xfree (ssl_c);
    if (ssr) xfree (ssr);
    if (ssr_c) xfree (ssr_c);
    return -1;
}

int get_fragment_end (char **fl, char **fr, char **fl_c, char **fr_c, int *fl_len, int *fr_len) {

    SEQUENCE *s;
    char *ffl, *ffr, *ffl_c, *ffr_c;

    s = get_editor_buffer ();
  
    if (s->left_end != NULL && s->left_end->len != 0) {
	int l = END_LEN + ABS(s->left_end->len);
	*fl_len = s->left_end->len;
	if (s->left_end->len > 0) {
	    if(NULL == (ffl = (char *)xmalloc((l + 1) * sizeof(char))))
		goto err;
	    if(NULL == (ffl_c = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
		goto err;
	    strncpy (ffl, &s->seq[0], l);
	    ffl[l] = 0;
	    *fl = ffl;
	    strncpy (ffl_c, &s->seq[s->left_end->len], END_LEN);
	    ffl_c[END_LEN] = 0;
	    complement_dna (ffl_c, END_LEN);
	    *fl_c = ffl_c;
	} else {
	    if(NULL == (ffl = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
		goto err;
	    if(NULL == (ffl_c = (char *)xmalloc((l + 1) * sizeof(char))))
		goto err;
	    strncpy (ffl, &s->seq[0], END_LEN);
	    ffl[END_LEN] = 0;
	    *fl = ffl;
	 
	    strcpy (ffl_c, s->left_end->hang);
	    strcat (ffl_c, ffl);
	    ffl_c[l] = 0;
	    complement_dna (ffl_c, l);
	    *fl_c = ffl_c;
	}
    } else {
	*fl_len = 0;
	
	if(NULL == (ffl = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
	    goto err;
	if(NULL == (ffl_c = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
	    goto err;
	strncpy (ffl, &s->seq[0], END_LEN);
	ffl[END_LEN] = 0;
	*fl = ffl;

	strncpy (ffl_c, &s->seq[0], END_LEN);
	ffl_c[END_LEN] = 0;
	complement_dna (ffl_c, END_LEN);
	*fl_c = ffl_c;
    }
    if (s->right_end != NULL && s->right_end->len != 0) {
	int l = END_LEN + ABS(s->right_end->len);
	int ll;
	*fr_len = s->right_end->len;
	if (s->right_end->len > 0) {
	 
	    if(NULL == (ffr = (char *)xmalloc((l + 1) * sizeof(char))))
		goto err;
	    if(NULL == (ffr_c = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
		goto err;
	    ll = s->length - l; 
	    strncpy (ffr, &s->seq[ll], l);
	    ffr[l] = 0;
	    *fr = ffr;
	    strncpy (ffr_c, &s->seq[ll], END_LEN);
	    ffr_c[END_LEN] = 0;
	    complement_dna (ffr_c, END_LEN);
	    *fr_c = ffr_c;
	} else {
	  
	    if(NULL == (ffr = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
		goto err;
	    if(NULL == (ffr_c = (char *)xmalloc((l + 1) * sizeof(char))))
		goto err;
	    ll = s->length - END_LEN;
	    strncpy (ffr, &s->seq[ll], END_LEN);
	    ffr[END_LEN] = 0;
	    *fr = ffr;
	    strcpy (ffr_c, *fr);
	    strcat (ffr_c, s->right_end->hang);
	    ffr_c[l] = 0;
	    complement_dna (ffr_c, l);
	    *fr_c = ffr_c;
	}
    } else {
	int ll;
	*fr_len = 0;
	if(NULL == (ffr = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
	    goto err;
	if(NULL == (ffr_c = (char *)xmalloc((END_LEN + 1) * sizeof(char))))
	    goto err;
	ll = s->length - END_LEN;	
	strncpy (ffr, &s->seq[ll], END_LEN);
	ffr[END_LEN] = 0;
	*fr = ffr;
	strncpy (ffr_c, &s->seq[ll], END_LEN);
	ffr_c[END_LEN] = 0;
	complement_dna (ffr_c, END_LEN);
	*fr_c = ffr_c;
    }
    
    return 0;
 err:
    if (ffl) xfree (ffl);
    if (ffl_c) xfree (ffl_c);
    if (ffr) xfree (ffr);
    if (ffr_c) xfree (ffr_c);
    return -1;
}

/*void convert_to_upper (char *s) {
    
    int i, len;
    len = strlen (s);
    for (i = 0; i < len; i++) {
	s[i] = toupper(s[i]);
    }
    }*/

int get_sequence_site1_end (SEQUENCE *s, SITE_HANG *sh, int len, char **sl, char **sl_c){

  char *ssl, *ssl_c;
  int lll = END_LEN + len;
  int ll = sh->pos - END_LEN;
 
  if(NULL == (ssl = (char *)xmalloc((lll + 1) * sizeof(char))))
    goto err;
  if(NULL == (ssl_c = (char *)xmalloc((lll + 1) * sizeof(char))))
    goto err;
  strncpy (ssl, &s->seq[ll], lll);
  ssl[lll] = 0;
  *sl = ssl;
  strncpy (ssl_c, &s->seq[ll], lll);
  ssl_c[lll] = 0;
  complement_dna (ssl_c, lll);
  *sl_c = ssl_c;
 
  return 0;
 err: 
  if (ssl) xfree(ssl);
  if (ssl_c) xfree(ssl_c);
  return -1;
}

int get_sequence_site2_end (SEQUENCE *s, SITE_HANG *sh, int len, char **sr, char **sr_c){

  char *ssr, *ssr_c;
  int ll = END_LEN + len;
  int cut_site = sh->pos;
 
  if(NULL == (ssr = (char *)xmalloc((ll + 1) * sizeof(char))))
    goto err;
  if(NULL == (ssr_c = (char *)xmalloc((ll + 1) * sizeof(char))))
    goto err;
  strncpy (ssr, &s->seq[cut_site], ll);
  ssr[ll] = 0;
  *sr = ssr;
  strncpy (ssr_c, &s->seq[cut_site], ll);
  ssr_c[ll] = 0;
  complement_dna (ssr_c, ll);
  *sr_c = ssr_c;
 
   return 0;
 err: 
  if (ssr) xfree(ssr);
  if (ssr_c) xfree(ssr_c);
  return -1;
}

EDIT *create_trim_edit (EDITOR_RECORD *er, 
			int member_id, 
			int pos,
			SEQUENCE *fd,
			int flag) {

    EDIT *edit;

    edit = create_edit();
    edit->seq_id    = member_id;
    edit->position  = pos;
    edit->consensus = flag; /*FIXME*/  
    edit->operation = TRIM;
    edit->sequence  = copy_sequence (fd);
    
    return edit;
}

void reset_sequence_hang_trim (EDITOR_RECORD *er, EDIT *edit) {

    SEQUENCE *s, *cp;
    SITE_HANG *sh1, *sh2;
    int group_id = 0;

    if (edit->seq_id != -1) {
	s = er->seq_group[group_id]->members[edit->seq_id]->data->sequence;
    }
    cp = get_editor_buffer ();

    if (edit->consensus == 1) {
	if (edit->seq_id != -1) sh1 = s->cut_site_1;
	if (edit->seq_id == -1) sh1 = cp->right_end;
	sh1->len = strlen (edit->sequence->seq);
	sh1->hang = strdup (edit->sequence->seq);
	sh1->pos += sh1->len;
    }
    if (edit->consensus == -1) {
	if (edit->seq_id != -1) sh1 = s->cut_site_1;
	if (edit->seq_id == -1) sh1 = cp->right_end;
	sh1->len = -strlen (edit->sequence->seq);
	sh1->hang = strdup (edit->sequence->seq);
    }
    if (edit->consensus == 2) {
	if (edit->seq_id != -1) sh2 = s->cut_site_2;
	if (edit->seq_id == -1) sh2 = cp->left_end;
	sh2->len = strlen (edit->sequence->seq);
	sh2->hang = strdup (edit->sequence->seq);
    }
    if (edit->consensus == -2) {
	if (edit->seq_id != -1) sh2 = s->cut_site_2;
	if (edit->seq_id == -1) sh2 = cp->left_end;
	sh2->len = -strlen (edit->sequence->seq);
	sh2->hang = strdup (edit->sequence->seq);
    }
}

int un_trim_sequence (EDITOR_RECORD *er, EDIT *edit, int add_to_undo) {

    int err;
    int group_id = 0; /*FIXME*/

    if (edit->consensus > 0) {
	if (edit->seq_id != -1) {
	    err = insert_fragment_in_editor (er, edit, add_to_undo);
	    if (err) goto bail_out;
	    err = update_consensus (er->seq_group[group_id], 1);
	    if (err) goto bail_out;
	}
	if (edit->seq_id == -1) {
	    SEQUENCE *cp;
	    cp = get_editor_buffer ();
	    err = insert_fragment_in_sequence (cp, edit->sequence, edit->position, 0);
	}
    }
    reset_sequence_hang_trim (er, edit);
    err = remove_edit (er->edits);
    if (err) goto bail_out;
    return err;
 bail_out: 
    return -1;
}
int trim_update_sequence (int ed_id, int member_id, int start, SEQUENCE *fd, int flag) {

    int  err, add_to_undo = 1;
    EDIT *edit;
    EDITOR_RECORD *er;
    
    er = editor_records->editor_record[ed_id];
    edit = create_trim_edit (er, member_id, start, fd, flag);
    if (flag > 0) {
	if (member_id != -1) { /* editor */ 
	    err = delete_fragment_in_editor (er, edit, add_to_undo);
	    return err;
	}
	if (member_id == -1) { /* clipboard */
	    SEQUENCE *cp;
	    FEATURE_TABLE *ft;
	 
	    cp = get_editor_buffer ();
	    ft = delete_fragment_in_sequence (cp, fd, start);
	    if (ft && edit->seq_id != 0) {
		err = add_ft_in_sequence (edit->sequence, ft);
		if (err == -1) return -1;;
	    }
	    if (add_to_undo) {
		err = add_edit (er->edits, edit);
		if (err) return -1;;
	    }
	    return err;
	}
    }
    if (flag < 0) { /* add to undo */
	if (add_to_undo) {
	    err = add_edit (er->edits, edit);
	    if (err) return -1;
	}
    }
    return 0;
}

/* flag: 1:cut_site_1 & right_end; -1:cut_site_1 & right_end(reverse sequence);
         2:cut_site_2 & left_end;  -2:cut_site_2 & left_end(reverse sequence);
*/
int trim_sequence_site_1 (int ed_id, int group_id, int member_id, SITE_HANG *sh) {

    SEQUENCE *fd;
    EDITOR_RECORD *er;
    int seq_id, seq_num;
    int start, end;
    int err;

    if (member_id == -1) {
	seq_num = -1;
    } else {
	er = GetEditor (ed_id);
	seq_id = er->seq_group[group_id]->members[member_id]->data->sequence->id;
	seq_num = GetEdenNum (seq_id);
    }
    if (sh->len > 0) {
	start = sh->pos - sh->len;
	end = sh->pos;
	fd = create_fragment_from_sequence (seq_num, start, end, "", "");
	err = trim_update_sequence (ed_id, member_id, start, fd, 1);
	sh->pos -= sh->len;
    }
    if (sh->len < 0) {
	
	if (sh->hang != NULL) fd = init_sequence (); 
	fd->seq = strdup (sh->hang);	
	err = trim_update_sequence (ed_id, member_id, sh->pos, fd, -1);
    }
    if (err == 0) {
	sh->len = 0;
	sh->hang = NULL;
    }
    return err;
}

void trim_sequence_site_2 (int ed_id, int group_id, int member_id, SITE_HANG *sh) {

    SEQUENCE *fd;
    EDITOR_RECORD *er;
    int seq_id, seq_num;
    int start, end;
    int err;
    
    if (member_id == -1) {
	seq_num = -1;
    }
    else {
	er = GetEditor (ed_id);
	seq_id = er->seq_group[group_id]->members[member_id]->data->sequence->id;
	seq_num = GetEdenNum (seq_id);
    }
    
    if (sh->len > 0) { /* forward */
	start = sh->pos;
	end = sh->pos + sh->len;
	fd = create_fragment_from_sequence (seq_num, start, end, "", "");
	err = trim_update_sequence (ed_id, member_id, start, fd, 2);
    }
    if (sh->len < 0) {
	if (sh->hang != NULL) fd = init_sequence (); 
	fd->seq = strdup (sh->hang);
	err = trim_update_sequence (ed_id, member_id, sh->pos, fd, -2);
    }
    if (err == 0) {
	sh->len = 0;
	sh->hang = NULL;
    }
}

char *extend_string (char *sl, int dir, char *dot) {

  int sl_l = strlen(sl);
  int dot_l = strlen(dot);

  if (dir == 0) { 
    sl = (char *)xrealloc(sl, sl_l + 1 + dot_l + 1);
    memmove (&sl[dot_l], &sl[0], sl_l);
    memmove (&sl[0], &dot[0], dot_l);
    sl[sl_l + dot_l] = 0;
    return sl;
  }
  if (dir == 1) {
    sl = (char *)xrealloc(sl, sl_l + 1 + dot_l + 1);
    strncpy(&sl[sl_l], dot, dot_l);
    sl[sl_l + dot_l] = 0;
    return sl;
  }
  return NULL;
}

int trim_sequence_hang (int ed_id, int group_id, int member_id, char *type, int *pos) {

    SEQUENCE *s, *cp;

    if (member_id != -1 )
	s = editor_records->editor_record[ed_id]->seq_group[group_id]->members[member_id]->data->sequence;
    
    cp = get_editor_buffer ();
    if (!strcmp (type, "sl")) {
	if (s->cut_site_1 != NULL) {
	    if (s->cut_site_1->len > 0) 
		s->cut_site_2->pos -= s->cut_site_1->len;
	    trim_sequence_site_1 (ed_id, group_id, member_id, s->cut_site_1);
	    *pos = s->cut_site_1->pos;  
	} else goto err;
    }
    if (!strcmp (type, "sr")) {
	if (s->cut_site_2 != NULL) {
	    trim_sequence_site_2 (ed_id, group_id, member_id, s->cut_site_2); 
	    *pos = s->cut_site_2->pos;
	} else goto err;	
    }
    if (!strcmp (type, "fl")) {
	if (cp->left_end != NULL) {
	    cp->right_end->pos = cp->length;
	    cp->left_end->pos = 0;
	    trim_sequence_site_2 (ed_id, group_id, -1, cp->left_end);
	    *pos = 0; /* FIXME return -1 means editing clipboard */
	
	} else goto err;	
    }
    if (!strcmp (type, "fr")) {
	if (cp->right_end != NULL) {
	    cp->right_end->pos = cp->length;
	    cp->left_end->pos = 0;
	    trim_sequence_site_1 (ed_id, group_id, -1, cp->right_end); 
	    *pos = 0;
	} else goto err;
    }
    return 0;
 err:
    return -1;
}

void reset_sequence_hang_fill (EDITOR_RECORD *er, EDIT *edit) {

    
    SEQUENCE *s, *cp;
    int group_id = 0; /*FIXME */

    if (edit->seq_id != -1)
	s = er->seq_group[group_id]->members[edit->seq_id]->data->sequence;
    
    cp = get_editor_buffer ();
    if (edit->consensus == 1) {
	if (edit->seq_id != -1) {
	    s->cut_site_1->len = -strlen (edit->sequence->seq);
	    s->cut_site_1->hang = strdup (edit->sequence->seq);
	    s->cut_site_1->pos -= abs (s->cut_site_1->len);
	    s->cut_site_2->pos -= abs (s->cut_site_1->len);
	}
	if (edit->seq_id == -1) {
	    cp->right_end->len = -strlen (edit->sequence->seq);
	    cp->right_end->hang = strdup (edit->sequence->seq);
	    cp->right_end->pos -= abs (cp->right_end->len);
	}
    }
    if (edit->consensus == -2) {
	if (edit->seq_id != -1) {
	    s->cut_site_2->len = strlen (edit->sequence->seq);
	    s->cut_site_2->hang = strdup (edit->sequence->seq);
	}
	if (edit->seq_id == -1) {
	    cp->left_end->len = strlen (edit->sequence->seq);
	    cp->left_end->hang = strdup (edit->sequence->seq);
	}
    }
}

int un_fill_sequence (EDITOR_RECORD *er, EDIT *edit, int add_to_undo) {

    int err;
    int group_id = 0; /* FIXME: use edit to get it*/

    if (edit->consensus > 0) {
	if (edit->seq_id != -1) {
	    err = delete_fragment_in_editor (er, edit, add_to_undo);
	    if (err) goto bail_out;
	    err = update_consensus (er->seq_group[group_id], 1);
	    if (err) goto bail_out;
	}
	if (edit->seq_id == -1) {
	    SEQUENCE *cp;
	    FEATURE_TABLE *ft;
	    
	    cp = get_editor_buffer ();
	    ft = delete_fragment_in_sequence (cp, edit->sequence, edit->position);
	}
    }
    /* for filling */
    reset_sequence_hang_fill (er, edit);
    err = remove_edit (er->edits);
    if (err) goto bail_out;
    return err;
 bail_out: 
    return -1;
}

EDIT *create_fill_edit (EDITOR_RECORD *er, int member_id, int pos, SEQUENCE *fi, int flag) {

    EDIT *edit;

    edit = create_edit();
    edit->seq_id    = member_id;
    edit->position  = pos;
    edit->consensus = flag; /*FIXME*/  
    edit->operation = FILL;
    edit->sequence  = copy_sequence (fi);
    
    return edit;
}

int fill_update_sequence (int ed_id, int member_id, int pos, SEQUENCE *fi, int flag) {

    int  err, add_to_undo = 1;
    EDIT *edit;
    EDITOR_RECORD *er;
  
    er = GetEditor (ed_id);
    edit = create_fill_edit (er, member_id, pos, fi, flag);   
    if (flag > 0) {
	if (member_id != -1) { /* editor */ 
	    
	    err = insert_fragment_in_editor (er, edit, add_to_undo);
	}
	if (member_id == -1) { /* clipboard */
	    SEQUENCE *cp;
	    SITE_HANG *sh;
	 
	    cp = get_editor_buffer ();
	    sh = cp->right_end;
	  
	    err = insert_fragment_in_sequence (cp, fi, sh->pos, 0);
	    if (add_to_undo) {
		err = add_edit (er->edits, edit);
		if (err) return -1;;
	    } 
	} 	
    }
    if (flag < 0) { /* for reverse sequence just modify hang, add to undo */
	if (add_to_undo) {
	    err = add_edit (er->edits, edit);
	    if (err) return -1;;
	}
    }
    return err;
}

/*flag: 1 for site_1 & right_end; -1 for site_1 & right_end (reverse)
        2 for site_2 & left_end;  -2 for site_2 & left_end (reverse)
*/
int fill_sequence_site_1 (int ed_id, int member_id, SITE_HANG *sh, int flag) {

  SEQUENCE *s;
  int err;
   
  if (sh->hang != NULL) s = init_sequence ();
  s->seq = strdup (sh->hang);
  err = fill_update_sequence (ed_id, member_id, sh->pos, s, flag);

  /* reset hang */
  if (err == 0) {
      sh->pos += ABS(sh->len); 
      sh->len = 0;
      sh->hang = NULL;
  }
  return 0;
}

int fill_sequence_site_2 (int ed_id, int member_id, SITE_HANG *sh, int flag) {

    SEQUENCE *si;
    int err;
       
    if (sh->hang != NULL) si = init_sequence ();
    si->seq = strdup (sh->hang);
    err = fill_update_sequence (ed_id, member_id, sh->pos, si, flag);

    /* reset hang */
    if (err == 0) {
	sh->len = 0;
	sh->hang = NULL;
    }
    return 0;
}

int fill_sequence_hang (int ed_id, int group_id, int member_id, char *type, int *pos) {

    SEQUENCE *s, *cp;

    s = editor_records->editor_record[ed_id]->seq_group[group_id]->members[member_id]->data->sequence;
  
    cp = get_editor_buffer ();

    if (!strcmp (type, "sl")) {
	if (s->cut_site_1 != NULL) {
	    if (s->cut_site_1->len < 0) {
		s->cut_site_2->pos += ABS(s->cut_site_1->len);
		fill_sequence_site_1 (ed_id, member_id, s->cut_site_1, 1);
		*pos = s->cut_site_1->pos;
	    } else goto err1;
	} else goto err2;
    }
    if (!strcmp (type, "sr")) {
	if (s->cut_site_2 != NULL) {
	    if (s->cut_site_2->len > 0) {
		fill_sequence_site_2 (ed_id, member_id, s->cut_site_2, -2);
		*pos = s->cut_site_2->pos;
	    } else goto err1;
	} else goto err2;
    }
    if (!strcmp (type, "fl")) {
	if (cp->left_end != NULL) {
	    if (cp->left_end->len > 0) { 
		cp->left_end->pos = 0;	
		fill_sequence_site_2 (ed_id, -1, cp->left_end, -2);
		*pos = 0; /* return -1 means editing clipboard */
	    } else goto err1;
	}else goto err2;
    }
    if (!strcmp (type, "fr")) {
	if (cp->right_end != NULL) {
	    if (cp->right_end->len < 0) {
		cp->right_end->pos = cp->length;
		fill_sequence_site_1 (ed_id, -1, cp->right_end, 1);
		*pos = 0;
	    } else goto err1;
	} else goto err2;
    }
    return 0;
 err1:
    printf ("error: can't fill\n");	
    return -1;
 err2:
    printf ("error: blunt end\n");	
    return -1;
}

int EditSequenceHang (ClientData clientData, 
		      Tcl_Interp *interp, 
		      int argc, 
		      char *argv[]) 
{
    SEQUENCE *s;
    int ed_id, seq_id, seq_num, group_id, member_id;    
    char *operation, *type;
    int cut_pos;
    int err;

    if (argc != 4) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "selected sequence\"",(char*)NULL);
        return TCL_ERROR;
    }

    operation = argv[1];
    seq_id = atoi(argv[2]);
    type = argv[3];
    seq_num = GetEdenNum(seq_id);
    ed_id = GetEdIdFromSeqId (seq_id);
    group_id = GetGroupIdFromSeqId (seq_id);
    member_id = GetMemberIdFromSeqId (seq_id);
    s = GetEdenSequence (seq_num);

    if (!strcmp (operation, "trim")) {
	err = trim_sequence_hang (ed_id, group_id, member_id, type, &cut_pos);
    }
    if (!strcmp (operation, "fill")) {
	err = fill_sequence_hang (ed_id, group_id, member_id, type, &cut_pos);
    }
    /* notify */
    if (err == 0) {
	if (cut_pos != -1) {
	    sequences_notify_graphic (seq_num, cut_pos);
	}
    } 
    /*else {
	printf ("can't fill or trim\n"); 
	}*/
    return TCL_OK; 
}

int GetHangParameter(ClientData clientData, 
		     Tcl_Interp *interp, 
		     int argc, 
		     char *argv[]) 
{
  int seq_id;
  char *sl, *sr, *sl_c, *sr_c;
  int sl_len, sr_len;
  char *fl, *fr, *fl_c, *fr_c;
  int fl_len, fr_len;
  char buf[1024];
  
  if (argc != 2) {
    Tcl_AppendResult(interp, "wrong # args: should be \"",
		     argv[0], "selected sequence\"",(char*)NULL);
    return TCL_ERROR;
  }    
  seq_id = atoi(argv[1]);
  
  get_sequence_end (seq_id, &sl, &sr, &sl_c, &sr_c, &sl_len, &sr_len);
  get_fragment_end (&fl, &fr, &fl_c, &fr_c, &fl_len, &fr_len);	    

  convert_to_upper (fl);
  convert_to_upper (fl_c);
  convert_to_upper (fr);
  convert_to_upper (fr_c);
  sl = extend_string (sl, 0, "...");
  sl_c = extend_string (sl_c, 0, "...");
  sr = extend_string (sr, 1, "...");
  sr_c = extend_string (sr_c, 1, "...");
  fl = extend_string (fl, 1, "......");
  fl_c = extend_string (fl_c, 1, "......");

  sprintf (buf, "%s %s %d %s %s %d %s %s %d %s %s %d", 
      sl, sl_c, sl_len, sr, sr_c, sr_len, fl, fl_c, fl_len, 
      fr, fr_c, fr_len);

  vTcl_SetResult(interp, "%s", buf);
  return TCL_OK; 
}

void end_editor_redisplay (Tcl_Interp *interp, int seq_num, text_editor_result *result) {

    char seqid[10];
    char *win_name = result->win_name;
    int ed_id, seq_id;
    
    seq_id = GetEdenId (seq_num);
    ed_id = GetEdIdFromSeqId (seq_id);
    sprintf (seqid, "%d", seq_id);

    if (result->ed_id == ed_id) {
	if (TCL_OK != Tcl_VarEval(interp, "EndEditorRedisplay ", win_name, " ", seqid, " ", NULL)) {
	    fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
	}
    }
}

void static end_editor_callback(int seq_num, void *fdata, editor_reg_data *jdata) {

    text_editor_result *result = (text_editor_result *) fdata;

    switch(jdata->job) {
    case SEQ_CHANGED:
	{	    
	    end_editor_redisplay (result->interp, seq_num, result);
	    /*printf ("end_editor_callback\n");*/
	    break;
	}
     case END_EDITOR_QUIT:
	{    
	    int ed_id = result->ed_id;
	    char cmd[1024];
	    editor_deregister(ed_id, end_editor_callback, (text_editor_result *)result);
	    
	    sprintf(cmd, "DeleteEndEditor %s", result->frame_name);
	    if (TCL_ERROR == Tcl_Eval(result->interp, cmd)) {
		verror(ERR_WARN, "end text editor ", "shutdown %s\n", 
		       result->interp->result);
	    }
	    /*printf ("end_editor:END_EDITOR_QUIT\n");*/
	    break;
	}
    }	
}

int EndEditorRegister (ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc, char **argv) {

   char *frame_name, *ee_win;
   text_editor_result *result;
   int id;
   int seq_num, seq_id, ed_id;
    
   if (argc != 4) {
       Tcl_AppendResult(interp, "wrong # args: should be \"",
			argv[0], "seq_ids\"", (char*)NULL);
       return TCL_ERROR;
   }

   frame_name = argv[1];
   ee_win = argv[2];
   seq_id = atoi (argv[3]);
   ed_id = GetEdIdFromSeqId (seq_id); 
   seq_num = GetEdenNum(seq_id);

   result = init_text_editor_result ();
   
   strcpy(result->frame_name, frame_name);
   strcpy(result->win_name, ee_win);
   result->seq_id = seq_id;
   result->ed_id = ed_id;/*editor ID*/
   result->op_func = end_editor_callback;
   result->interp = interp;  
   id = get_editor_reg_id();
   result->index = id;/*register ID*/
			
   editor_register (seq_num, end_editor_callback, (void *)result, TEXT, id);    
  
   return TCL_OK;  
}


int Hang_editor_Init(Tcl_Interp *interp) {
    
    Tcl_CreateCommand(interp, "end_editor_register", EndEditorRegister, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    if (Itcl_RegisterC(interp, "get_hang_para", GetHangParameter, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "edit_sequence_hang", EditSequenceHang, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    return TCL_OK;
}


