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
#include <ctype.h>

#include "tcl_utils.h"
#include "misc.h"
#include "cli_arg.h"
#include "dna_utils.h"
#include "sequence_formats.h"
#include "nip_globals.h"
#include "spin_globals.h"
#include "read_sequence.h"
#include "editor.h"
#include "read_sequence.h"
#include "feature_table.h"
#include "renzyme_search.h"
#include "editor_reg.h"
#include "end_editor.h"
#include "seq_results.h"
#include "nip_results.h"
#include "seq_reg_structs.h"
#include "text_editor.h"
#include "feature_colour.h"
#include "feature_editor.h"
#include "graphic_editor.h"
#include "renzyme_map_canvas.h"


EDITOR_RECORDS *editor_records = NULL; /* text & graphical editor */
SEQUENCES *eden = NULL; /* registration for text and graphic editor */

void static editor_callback(int seq_num, void *fdata, editor_reg_data *jdata);

EDITS *create_edits(void) {

    EDITS *edits;

    if(NULL == (edits = (EDITS *) xmalloc(sizeof(EDITS)))) {
	return NULL;
    }
    edits->head = NULL;
    return edits;
}

EDIT *create_edit(void) {

    EDIT *edit;

    if(NULL == (edit = (EDIT *) xmalloc(sizeof(EDIT)))) {
	return NULL;
    }
    edit->seq_id    = -1;
    edit->position  = -1;
    edit->operation = -1;
    edit->consensus = -1;
    edit->string    = NULL;  
    edit->sequence  = NULL;  
    edit->next      = NULL;
    return edit;
}

/* add an edit to the list */

int add_edit(EDITS *edits, EDIT *edit) {
 
    if(edits->head) edit->next = edits->head;
    edits->head = edit;
    return 0;
}

void free_edits (EDITS *edits) {

    EDIT *p; 
    for (p = edits->head; p; p = p->next) {
	if (p->string) xfree (p->string);
	if (p->sequence) free_sequence (p->sequence);
	xfree(p);
    }
    xfree (edits);
}

cursors *init_cursors (void) {

    cursors *cs;

    if(NULL == (cs = (cursors *) xmalloc(sizeof(cursors)))) {
	return NULL;
    }
    cs->cursor = NULL;
    cs->next = NULL;
    return cs;
}

void free_cursors (cursors *cs) {
    
    
}

DATA *init_data_sequence (void) {

    DATA *d;
    
    if ( ( NULL == ( d = (DATA* ) xmalloc ( sizeof (DATA) )))) 
	return NULL;
    d->sequence = NULL;
    d->group = NULL;
    return d;
}

GROUP_MEMBER *init_group_member (void) {

  GROUP_MEMBER *gm;

  if ( ( NULL == ( gm = (GROUP_MEMBER* ) xmalloc ( sizeof (GROUP_MEMBER) )))) 
	return NULL;
  gm->data = NULL;
  gm->type = SEQ;
  return gm;
}

void free_group_member (GROUP_MEMBER *gm) {
    if (gm->type == SEQ) {
	free_sequence (gm->data->sequence);
	xfree (gm->data);
    } else if (gm->type == GROUP) {
	free_seq_group (gm->data->group);
    }
    xfree (gm);
}

void free_seq_group (SEQ_GROUP *sg) {

    int i, num_member;
    
    if (sg) {
	num_member = sg->nmembers;
	for (i = 0; i < num_member; i++) {
	    free_group_member (sg->members[i]);
	}
	if (sg->group_name) xfree (sg->group_name);
    }
}

/*void free_seq_group (SEQ_GROUP *sg) {

    int i, num_member;
    
    if (sg) {
	num_member = sg->nmembers;
	for (i = 0; i < num_member; i++) {
	    if (sg->members[i]->type == SEQ) {
		free_sequence (sg->members[i]->data->sequence);
	    } else if (sg->members[i]->type == GROUP) {
		printf("free group_member\n");
		free_seq_group (sg->members[i]->data->group);
	    }  
	}
	if (sg->group_name) xfree (sg->group_name);
    }
    }*/


SEQ_GROUP *init_seq_group (void) {

    SEQ_GROUP *sg;
    GROUP_MEMBER **gm;
    GROUP_MEMBER *gmm;
    DATA *d;

    if ( NULL == ( sg = (SEQ_GROUP* ) xmalloc ( sizeof (SEQ_GROUP) ))) return NULL; 
    if ( NULL == ( gm = (GROUP_MEMBER** ) xmalloc ( sizeof (GROUP_MEMBER *) ))) return NULL;
    if ( NULL == ( gmm = init_group_member())) return NULL;
    if ( NULL == ( d = init_data_sequence())) return NULL;
    
    sg->members = gm;
    sg->nmembers = 0; 
    sg->group_id = 0;
    sg->group_parent_id = 0;   
    sg->group_name = NULL;
    gmm->type = SEQ;
    gmm->data = d;
    sg->members[0] = gmm; /*first for storing consensus sequence*/
    return sg;
}

int realloc_seq_group (SEQ_GROUP *seq_group, int num_members ) {

    GROUP_MEMBER  **gm;
 
    gm = seq_group->members;
    if ((NULL == ( gm = (GROUP_MEMBER **)xrealloc(gm, sizeof(GROUP_MEMBER *)*(num_members+1))))) 
	return -1;
    seq_group->members = gm;
    seq_group->nmembers = num_members;
    return 0;
}

/* free an editor_record structure */
void free_editor_record ( EDITOR_RECORD *er) {

    int i, num_group;

    num_group = er->num_group;
    for (i = 0; i < num_group; i++) {
	free_seq_group ( er->seq_group[i] );
    }
    free_edits ( er->edits);    
    xfree (er);
}

EDITOR_RECORD *init_editor_record (void) {
    
    EDITOR_RECORD *er;
    SEQ_GROUP  **sg;
    EDITS *edits;
    
    if ( (NULL == (er = (EDITOR_RECORD* ) xmalloc ( sizeof (EDITOR_RECORD) )))) 
	return NULL;
    if ( (NULL == (sg = (SEQ_GROUP** ) xmalloc ( sizeof (SEQ_GROUP *) )))) 
	return NULL;
    if ( NULL == (edits = create_edits () ))
       return NULL;
    er->seq_group = sg;
    er->num_group = 0;
    er->ft_imode = 0;
    er->text = 0;
    er->graphical = 0;
    er->edits = edits;
    er->cursors = NULL;
    return er;

}

int realloc_editor_record (EDITOR_RECORD *er, int num_group) {

    SEQ_GROUP **sg;

    sg = er->seq_group;
    if ((NULL == (sg = (SEQ_GROUP **)xrealloc (sg, sizeof (SEQ_GROUP *) * (num_group + 1) )))) 
	return -1;
    er->seq_group = sg;
    er->num_group = num_group;
    return 0;
}

SELECTION *init_selection (void) {

    SELECTION *s;
    
    if ( NULL == (s = (SELECTION *) xmalloc ( sizeof (SELECTION)))) 
	return NULL;
    s->ed_id = 0;
    s->member_id = 0;
    s->sel_first = 0;
    s->sel_last = 0;
    s->name_first = NULL;
    s->name_last = NULL;
    return s;
}

void free_selection (SELECTION *s) {
    
    if (s->name_first) xfree (s->name_first);
    if (s->name_last) xfree (s->name_last);
    xfree (s);
}

void free_editor_records (EDITOR_RECORDS *ers) {

    int i, num_editor;
    
    num_editor = ers->num_editor;
    for (i = 0; i < num_editor; i++) {
	free_editor_record (ers->editor_record[i]);	
    }
    if (ers->clipboard) free_sequence (ers->clipboard);
    if (ers->selection) free_selection (ers->selection);
    xfree (ers);
}

EDITOR_RECORDS *init_editor_records (void) {
    
    EDITOR_RECORDS *ers;
    EDITOR_RECORD  **er;
    
    if ( (NULL == (ers = (EDITOR_RECORDS* ) xmalloc ( sizeof (EDITOR_RECORDS) )))) 
	return NULL;
    /* need to allow zero capacity! */
    if ( (NULL == (er = (EDITOR_RECORD** ) xmalloc ( sizeof (EDITOR_RECORD *) )))) 
	return NULL;
    ers->editor_record = er;
    ers->clipboard = NULL;
    ers->selection = NULL;
    ers->num_editor = 0;
    return ers;
}

int realloc_editor_records (EDITOR_RECORDS *ers, int num_editor) {

  EDITOR_RECORD  **er;
 
  er = ers->editor_record;
  if (NULL == (er = (EDITOR_RECORD **)xrealloc (er, sizeof(EDITOR_RECORD *)*(num_editor + 1)))) 
      return -1;
  ers->editor_record = er;
  ers->num_editor = num_editor;

  return 0;
}

SITE_HANG *copy_site_hang (SITE_HANG *sh) {

    SITE_HANG *h;
    h = init_site_hang ();
        
    h->pos = sh->pos;
    h->len = sh->len;
    if (sh->hang) h->hang = strdup (sh->hang); 
    return h;
}

/* made a copy for sequence */
SEQUENCE *copy_sequence (SEQUENCE *f) {

  SEQUENCE  *s;

  if (!f) return NULL;
  if (NULL == (s = init_sequence())) return NULL;

  if (f->seq) s->seq = strdup(f->seq);
  if (f->name) s->name = strdup(f->name);
  if (f->source) s->source = strdup(f->source);
  if (f->feature_table) s->feature_table = copy_feature_table(f->feature_table);
  s->length = f->length;
  s->start = f->start;
  s->end = f->end;
  s->type = f->type;
  s->direction = f->direction;
  if (f->left_end) s->left_end = copy_site_hang(f->left_end);
  if (f->right_end) s->right_end = copy_site_hang (f->right_end);
  if (f->cut_site_1) s->cut_site_1 = copy_site_hang (f->cut_site_1);
  if (f->cut_site_2) s->cut_site_2 = copy_site_hang (f->cut_site_2);
  s->id = f->id;
  s->genetic_code = f->genetic_code;
  return s;
}


int extend_unique_name (SEQUENCE *s) {

    char *name;
    int static next_complement_sequence = 0;

    next_complement_sequence++;
    name = s->name;
    if (NULL == (name = (char *)realloc(name, (strlen(name)+10))))
	return -1;
    sprintf(name, "%s_c%d", name, next_complement_sequence);
    s->name = name;
    return 0;
}

/** simple consensus generator
 *  consensus is highest scorer if >74% otherwise N
 *  DNA upper case only (FIXME)
 */

int sequence_adder (SEQ_GROUP *sg) {

    int i, j, column, row, element, length, top, top_index, depth, *freqs;
    int consensus_length;
    int num_member;
    char *seq;
    char chars[]={"ACGTN"};

    num_member = sg->nmembers;
    consensus_length = sg->members[0]->data->sequence->length;

    if (NULL == (freqs = ( int *)xcalloc(char_set_size*consensus_length,sizeof(int)))) return -1;
    
    for (i = 1; i <= num_member; i++) {
	if (sg->members[i]->type == SEQ) {
	    seq = sg->members[i]->data->sequence->seq;
	    column = sg->members[i]->data->sequence->offset;
	    length = sg->members[i]->data->sequence->length;
	    for (j = 0; j < length; j++) {
		row = char_lookup[(int)seq[j]]; /*char_lookup=dna_lookup, dna_lookup['a'] = 0 */
		element = row + char_set_size * column++;
		freqs[element] += 1;
	    }
	}
    }
    length = sg->members[0]->data->sequence->length;
    for (i = 0; i < length; i++) {
      for (row = 0, top = -1, depth = 0; row < char_set_size; row++) {
	  element = row + char_set_size * i;
	  if (freqs[element] > 0) {
	      depth += freqs[element];
	      if (freqs[element] > top) {
		  top = freqs[element];
		  top_index = row;
	      }
	  }
      }
       if ((float)top / depth > .60) {
	  sg->members[0]->data->sequence->seq[i] = chars[top_index];
      }
      else {
	  sg->members[0]->data->sequence->seq[i] = chars[char_set_size-1];
      }
  }
  sg->members[0]->data->sequence->seq[length] = 0;
  xfree(freqs);
  return 0;
}

int update_consensus (SEQ_GROUP *sg, int job ) {

    SEQUENCE *consensus;

    char *seq;
    char *name_0 = "CONSENSUS";
    int num_member;
    int i, consensus_length;

    if ( job == 1 ) {
	/* do it all */
	consensus_length = 0;
	num_member = sg->nmembers;
	for (i = 1;i <= num_member;i++) {
	    if (sg->members[i]->type == SEQ) {
		consensus_length = MAX(consensus_length, 
				       sg->members[i]->data->sequence->offset 
				       + sg->members[i]->data->sequence->length);
	    } else if (sg->members[i]->type == GROUP) {
		update_consensus (sg->members[i]->data->group, job);
	    } 
	}
	if (sg->members[0]->data->sequence) {  
	    if (sg->members[0]->data->sequence->seq ) 
		xfree (sg->members[0]->data->sequence->seq);
	    if ( ( NULL == ( seq = (char * ) xmalloc ( consensus_length + 1) ))) return -1;
	    sg->members[0]->data->sequence->seq = seq;
	}else {
	    /* initialise consensus */
	    if (NULL == (consensus = init_sequence())) return -1;
	    if (NULL == (consensus->seq = (char *)xmalloc((consensus_length + 1)*sizeof(char))))
		return -1;
	    sg->members[0]->data->sequence = consensus;
	    sg->members[0]->data->sequence->name = strdup(name_0);
	}
	sg->members[0]->data->sequence->length = consensus_length;
    }
    i = sequence_adder ( sg );
    return 0;
}

int create_copy_for_editor (Tcl_Interp *interp, int seq_id) {

    int i, num_seq;
    static int initialised_editor_reg = 0;

    if (!eden) {
	if (NULL == (eden = init_sequences( ))) 
	    goto err;
    }
    num_seq = eden->num_seq;

    /* initialise editor registration if first time thro' */    
    if (!initialised_editor_reg) {
	editor_register_init(interp);
	initialised_editor_reg = 1;
    }

    /* check if already in editor_sequence array */
    for (i = 0; i < num_seq; i++) {
	if (!strcmp (GetSequenceName (seq_id), eden->sequence[i]->name)) {
	    return i;
	}
    }
    /* if not, make a copy in editor_sequence array */
    if (num_seq != 0) {
	    eden = realloc_sequences (eden, num_seq);
    }
    eden->sequence[num_seq] = copy_sequence (GetSequencesSequence (seq_id));
 
    eden->num_seq ++;
    
    add_editor_reg(num_seq);

    return num_seq;

 err: 
    return -1;
}

EDIT *copy_edit (EDIT *e) {

    EDIT *e_new;
    e_new = create_edit();
    e_new->seq_id = e->seq_id;
    e_new->position = e->position;
    e_new->operation = e->operation;
    if (e->string) e_new->string = strdup (e->string);
    if (e->sequence) e_new->sequence = copy_sequence (e->sequence);
    e_new->next = NULL;
    return e_new;
}

EDITS *copy_edits (EDITS *es) {

    EDITS *es_new;
    EDIT *e;
    
    es_new = create_edits();
    for (e = es->head; e; e = e->next) {
	EDIT *ee;
	ee = copy_edit (e);
	add_edit (es_new, ee);	
    }
    return es_new;
}

/*EDITS *init_undo_list (EDITOR_RECORD *er) {

    EDITS *ul;
    int i, j, k, num_editor, num_seq, n_seq;
 
    num_editor = editor_records->num_editor;
    num_editor--;
    num_seq = er->seq_group->nmembers; 

    for (i = 1; i <= num_seq; i++) {
	printf ("num_seq=%d\n", num_seq);
	printf ("sss=%d\n", er->seq_group->members[i]->data->sequence->id);
	for (j = num_editor; j >= 0; j--) {
	     n_seq = editor_records->editor_record[j]->seq_group->nmembers;
	     printf ("j=%d n_seq=%d\n", j, n_seq);
	     for (k = 1; k <= n_seq; k++) {
		 if (er->seq_group->members[i]->data->sequence->id == 
		     editor_records->editor_record[j]->seq_group->members[k]->data->sequence->id 
		     && editor_records->editor_record[j]->edits->head != NULL) {
		     ul = copy_edits (editor_records->editor_record[j]->edits);
		     return ul;
		 }
	    }
	}
    }
    return NULL;
    }*/


int check_sequence_in_editor (char **sel_seq, int num_sel_seq) {

    int i, num_editor;
    int j, n_seq;
    int k, h, num_group;;
    	
    num_editor = editor_records->num_editor;
    if (num_editor == 0) return -2;

    for (i = 0; i < num_editor; i++) {
	EDITOR_RECORD *er;
	int m = 0;
	er = editor_records->editor_record[i];
	num_group = er->num_group;
	for (h = 0; h < num_group; h++) {
	    n_seq = er->seq_group[h]->nmembers;
	    for (j = 1; j <= n_seq; j++) {
		for (k = 0; k < num_sel_seq; k++) {
		    if (!strcmp (sel_seq[k], er->seq_group[h]->members[j]->data->sequence->name)) {
		    m++;
		}
	    }
	}
	}
	if (m == n_seq && m == num_sel_seq) return i; /*editor_id */
	if (m != 0) return -1; /* one of the selected sequence already in other editor */ 
    }
    return -2; /* create new editor */
}


int check_if_already_in_editor (int seq_id) {

    EDITOR_RECORD *er;
    int i, num_editor;
    int k, num_group;
    int j, n_seq;

    num_editor = editor_records->num_editor;
    if (num_editor == 0) return -1;

    for (i = 0; i < num_editor; i++) {
	er = editor_records->editor_record[i];
	num_group = er->num_group;
	for (k = 0; k < num_group; k++ ) {
	    n_seq = er->seq_group[k]->nmembers;
	    for (j = 1; j <= n_seq; j++) {
		if (er->seq_group[k]->members[j]->data->sequence->id == seq_id)  
		    return i;
	    }
	}
    }
	return -1;
}

int add_new_editor (Tcl_Interp *interp, int seq_id) {

    EDITOR_RECORD *er;
    int i, ee_id;
    int num_editor, num_group, num_seq, num_member;
    int editor_num;

    if (!editor_records) {
	if (NULL == (editor_records = init_editor_records( ))) 
	    goto err;
    }

    editor_num = check_if_already_in_editor (seq_id);
    
    if (editor_num == -1) {
	SEQ_GROUP *sg;
	num_seq = GetSequenceNums ();/*the number of the sequence in SEQUENCES*/
	if (NULL == (er = init_editor_record( ))) goto err;
	if (NULL == (sg = init_seq_group( ))) goto err;
       
	for (i = 0; i < num_seq; i++) {
	    if ( GetSequenceId (i) == seq_id ) {
		GROUP_MEMBER *gm;
		DATA *d;
		
		if (NULL == ( gm = init_group_member())) 
		    goto err;
		if (NULL == ( d = init_data_sequence())) 
		    goto err;
	       
		ee_id = create_copy_for_editor (interp, seq_id); /* registration */
		d->sequence = GetEdenSequence (ee_id); /*point to registration*/
		gm->type = SEQ;
		gm->data = d;
		sg->nmembers++;
		num_member = sg->nmembers;
		realloc_seq_group (sg, num_member);
		sg->members[num_member] = gm; 	
	    }
	}
	num_group = er->num_group;
	if (num_group != 0) {
	    realloc_editor_record (er, num_group);	
	}
	er->seq_group[num_group] = sg;
	er->num_group++;

	er->graphical = 1;
	/* create consensus in editor */
	if (update_consensus (er->seq_group[num_group], 1) == -1)
	    goto err;
	num_editor = editor_records->num_editor;
	if (num_editor != 0) {
	    realloc_editor_records (editor_records, num_editor);	
	}
	editor_records->editor_record[num_editor] = er;
	editor_records->num_editor++;
	return num_editor;
    } else {
	return editor_num;
    }    
 err:
    if (editor_records) free_editor_records(editor_records);
    return -1;
}

int create_new_editor (Tcl_Interp *interp, char **sel_seq, int num_sel_seq) {

    EDITOR_RECORD *er;
    int i, j, k, ee_id;
    int num_editor, num_group, num_seq, num_member;
    int editor_num;
    
    if (!editor_records) {
	if (NULL == (editor_records = init_editor_records( ))) 
	    goto err;
    }

    editor_num = check_sequence_in_editor (sel_seq, num_sel_seq);
    
    if (editor_num == -2) {
	SEQ_GROUP *sg;
	num_seq = GetSequenceNums ();/*the number of the sequence in SEQUENCES*/
	if (NULL == (er = init_editor_record( ))) goto err;
	if (NULL == (sg = init_seq_group( ))) goto err;
	k = 0;
	for (i = 0; i < num_seq; i++) {
	    for (j = 0; j < num_sel_seq; j++) {
		if (!strcmp (GetSequenceName (i), sel_seq[j])) {
		    GROUP_MEMBER *gm;
		    DATA *d;
		    int seq_id;

		    if (NULL == ( gm = init_group_member())) 
			goto err;
		    if (NULL == ( d = init_data_sequence())) 
			goto err;
		    k++;
		    seq_id = GetSequenceId (i);
		    ee_id = create_copy_for_editor (interp, seq_id); /* registration */
		    d->sequence = GetEdenSequence (ee_id); /*point to registration*/
		    gm->type = SEQ;
		    gm->data = d; 
		    sg->nmembers++;
		    num_member = sg->nmembers;
		    realloc_seq_group (sg, num_member);
		    sg->members[num_member] = gm; 
		}
	    }
	}
	num_group = er->num_group;
	if (num_group != 0) {
	    realloc_editor_record (er, num_group);	
	}
	er->seq_group[num_group] = sg;
	er->num_group++;
	er->text++ ; /*open text_editor*/
	/* create consensus in editor */
	
	if (update_consensus (er->seq_group[num_group], 1) == -1)
	    goto err;
	num_editor = editor_records->num_editor;
	if (num_editor != 0) {
	    realloc_editor_records (editor_records, num_editor);	
	}
	editor_records->editor_record[num_editor] = er;
	editor_records->num_editor++;
	return num_editor;
    } else 
	return editor_num;  
 err:
    if (editor_records) free_editor_records(editor_records);
    return -1;
}

int GetEdenIdByName (char *name) {

    int i, num_seq, name_len;

    num_seq = eden->num_seq;
    for (i = 0; i < num_seq; i++) {
	name_len = strlen (eden->sequence[i]->name);
	if (!strncmp(eden->sequence[i]->name, name, name_len)) {
	    return eden->sequence[i]->id;
	}
    }
    return -1;
}

char *GetEdenSeq (int seq_num) {

    return eden->sequence[seq_num]->seq;
}

int GetEdenLength (int seq_num) {

    return eden->sequence[seq_num]->length;
}

int GetEdenType (int seq_num) {
    
    return eden->sequence[seq_num]->type;
}

char *GetEdenName (int seq_num) {
    
    return eden->sequence[seq_num]->name;
}
int GetEdenNums (void) {
    
    return eden->num_seq;
}

SEQUENCE *GetEdenSequence (int seq_num) {

    return eden->sequence[seq_num];
}

void replace_sequence (int seq_num, SEQUENCE *s) {

    free_sequence (eden->sequence[seq_num]);
    eden->sequence[seq_num] = s;

}
/*
 * return the unique id num of sequence at position seq_num in eden.
 */
int GetEdenId (int seq_num) {

    int num_seq;

    num_seq = eden->num_seq;
    if (seq_num < num_seq && seq_num > -1) {
	return eden->sequence[seq_num]->id;
    } 
    return -1;
}
/*
 * return the position of seq in eden array from a given unique id num
 */
int GetEdenNum(int seq_id) {

    int i, num_seq;

    num_seq = eden->num_seq;
   
    for (i = 0; i < num_seq; i++) {
	if (eden->sequence[i]->id == seq_id) {
	    return i;
	}
    }
    return -1;
}

/*EDITOR_RECORD *er;
    cursor_e *c;

    er = editor_records->editor_record[ed_id];
    for (c = er->cursor; c; c = c->next) { 
	if (c == NULL) {
	    c = cursor;
	}
	}*/

cursor_e *get_editor_cursor (int ed_id, int member_id, char *frame_name) {

    EDITOR_RECORD *er;
    cursor_e *ce = NULL;
    cursors *c;
    
    er = editor_records->editor_record[ed_id];
    if (er== NULL)  return NULL;
    if (er->cursors==NULL) return NULL; 
    c = er->cursors;

    
    while (c != NULL) {
	if (c->cursor->frame_name != NULL) {
	    if (!strcmp (c->cursor->frame_name, frame_name)) {
		ce = c->cursor;
		return ce;
	    }
	} else {   
	    if (c->cursor->posy == member_id) {
		ce = c->cursor;
		return ce;
	    }
	}
	c = c->next;
    }
    return ce;
}

int GetEditorRegNum (ClientData clientData, 
		     Tcl_Interp *interp, 
		     int argc, char **argv) {

    int seq_id;

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "seq_id\"", (char*)NULL);
        return TCL_ERROR;
    }
    seq_id = atoi (argv[1]);
    Tcl_ResetResult(interp);  
  
    vTcl_SetResult(interp, "%d", GetEdenNum(seq_id));
    return TCL_OK;
}

int AddEditor (ClientData clientData, 
	       Tcl_Interp *interp, 
	       int argc, char **argv) 
{
    char **sel_seq = NULL;
    int num_sel_seq;
    int ed_id;

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "seq_ids\"", (char*)NULL);
        return TCL_ERROR;
    }

    /* create selected sequence name array */
    if (Tcl_SplitList(interp, argv[1], &num_sel_seq, &sel_seq) != TCL_OK)
	return TCL_ERROR;

    ed_id = create_new_editor (interp, sel_seq, num_sel_seq);
    ckfree((char*)sel_seq);

    Tcl_ResetResult(interp);    
    vTcl_SetResult(interp, "%d", ed_id);
    return TCL_OK;
}


/* return the editor for giving editor id */
EDITOR_RECORD *GetEditor (int id) {

  return editor_records->editor_record[id];
}

int GetEditorNum (void) {
  
  return editor_records->num_editor;
}

char *GetEditorSeq (EDITOR_RECORD *er, int group_id, int member_id) {

    return er->seq_group[group_id]->members[member_id]->data->sequence->seq;
}
char *GetEditorSeqName (EDITOR_RECORD *er, int group_id, int member_id) {

    return er->seq_group[group_id]->members[member_id]->data->sequence->name;
}

int GetEditorSeqLength (EDITOR_RECORD *er, int group_id, int member_id) {

    return er->seq_group[group_id]->members[member_id]->data->sequence->length;
}

int GetEditorSeqId (int ed_id, int group_id, int member_id) {

    if (ed_id >= 0 && group_id >= 0 && member_id >= 0) {
	return editor_records->editor_record[ed_id]->seq_group[group_id]->members[member_id]->data->sequence->id;
    } else 
	return -1; 
}

int GetEdIdFromSeqId (int seq_id) {

    int i, num_editor;
    int k, num_group;

    num_editor = editor_records->num_editor;
    for (i = 0; i < num_editor; i++) {
	EDITOR_RECORD *er;
	er = editor_records->editor_record[i];
	num_group = er->num_group;
	for (k = 0; k < num_group; k++) {
	    int j, num_seq;
	    num_seq = er->seq_group[k]->nmembers;
	    for (j = 1; j <= num_seq; j++) {
		if (seq_id == er->seq_group[k]->members[j]->data->sequence->id) return i;
	    } 
	}
    }
    return -1;
}

int GetGroupIdFromSeqId (int seq_id) {

    int i, num_editor;

    num_editor = editor_records->num_editor;
    for (i = 0; i < num_editor; i++) {
	EDITOR_RECORD *er;
	int k, num_group;
	er = editor_records->editor_record[i];
	num_group = er->num_group;
	for (k = 0; k < num_group; k++) {
	    int j, num_seq;
	    num_seq = er->seq_group[k]->nmembers;
	    for (j = 1; j <= num_seq; j++) {
		if (seq_id == er->seq_group[k]->members[j]->data->sequence->id) return k;
	    } 
	}
    }
    return -1;
}

int GetMemberIdFromSeqId (int seq_id) {

    int i, num_editor;

    num_editor = editor_records->num_editor;
    for (i = 0; i < num_editor; i++) {
	EDITOR_RECORD *er;
	int k, num_group;
	er = editor_records->editor_record[i];
	num_group = er->num_group;
	for (k = 0; k < num_group; k++) {
	    int j, num_seq;
	    num_seq = er->seq_group[k]->nmembers;
	    for (j = 1; j <= num_seq; j++) {
		if (seq_id == er->seq_group[k]->members[j]->data->sequence->id) return j;
	    }
	}
    }
    return -1;
}

FEATURE_TABLE *GetEditorSeqFt (EDITOR_RECORD *er, int group_id, int member_id) {

    
    return er->seq_group[group_id]->members[member_id]->data->sequence->feature_table; 
}

int GetSeqMaxLen (ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, char **argv) {

    int len, id;
    int i, num_group;

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "editor_id\"", (char*)NULL);
        return TCL_ERROR;
    }
    Tcl_ResetResult(interp);
    id = atoi (argv[1]);
    num_group = editor_records->editor_record[id]->num_group;
    len = 0;
    for (i = 0; i < num_group; i++) {
	len = MAX (len, editor_records->editor_record[id]->seq_group[i]->members[0]->data->sequence->length);
    }    
    vTcl_SetResult(interp, "%d", len);
    return TCL_OK;
}

int GetNumSeqsInEditor (ClientData clientData, 
			Tcl_Interp *interp, 
			int argc, char **argv) 
{
    int id, num_seqs;
    int i, num_group;

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "editor_id\"", (char*)NULL);
        return TCL_ERROR;
    }
    Tcl_ResetResult(interp);
    id = atoi (argv[1]);

    num_group = editor_records->editor_record[id]->num_group;
    num_seqs = 0;
    for (i = 0; i < num_group; i++) {
	num_seqs += editor_records->editor_record[id]->seq_group[i]->nmembers;
    }
    vTcl_SetResult(interp, "%d", num_seqs);
    return TCL_OK;
}

int GetEdIdBySeqNum (int seq_num, int ed_id) {

    int i, j, k;
    int num_editor, num_group, seq_id, num_memb;
    
    seq_id = GetEdenId (seq_num);
    if (editor_records != NULL) {
	num_editor = editor_records->num_editor;
    } else num_editor = -1;
    if (ed_id >= num_editor ) return -1;
    for (i = ed_id; i < num_editor; i++) {
	num_group = editor_records->editor_record[i]->num_group;
	for (k = 0; k < num_group; k++) {
	    num_memb = editor_records->editor_record[i]->seq_group[k]->nmembers;
	    for (j = 1; j <= num_memb; j++) {
		if (seq_id == editor_records->editor_record[i]->seq_group[k]->members[j]->data->sequence->id)  
		    return i;
	    }
	}
    }
    return -1;
}

void editor_update (Tcl_Interp *interp, int ed_id, text_editor_result *result) {

    char cmd[1024];
    char *win_name = result->win_name;
    
    if (result->ed_id == ed_id) {
	sprintf(cmd, "%s update_display", win_name);    
	Tcl_Eval(interp, cmd);
    }
}

void editor_redisplay (Tcl_Interp *interp, int seq_num, text_editor_result *result) {

    int e_id = 0;
    
    while (e_id != -1) {
	e_id = GetEdIdBySeqNum (seq_num, e_id);	
	if (e_id != -1) {
	    editor_update (interp, e_id, result);
	    e_id ++;
	}
    }
}

void editor_icursor_move (Tcl_Interp *interp, text_editor_result *result, int posx, int posy) {

    char cmd[1024], px[20], py[20], win_name[100];

    sprintf (px, "%d", posx);
    sprintf (py, "%d", posy);
    sprintf (win_name, "%s", result->win_name);

    sprintf (cmd, "%s move_icursor %s %s", win_name, px, py);    
    Tcl_Eval(interp, cmd);
}

void text_editor_draw_selection (Tcl_Interp *interp, 
				 text_editor_result *result, 
				 int sel_first, 
				 int sel_last) {

    char cmd[1024], win_name[100];
   
    sprintf (win_name, "%s", result->win_name);
    if (sel_first != 0 && sel_last != 0) {
	sprintf (cmd, "%s text_draw_selection %d %d", win_name, sel_first, sel_last);
    } else {
	sprintf (cmd, "%s select clear", win_name);
	
    }
    Tcl_Eval(interp, cmd);
}

void delete_editor (int ed_id) {

    EDITOR_RECORDS *ers;
    int i, num_editor;

    ers = editor_records;
    num_editor = ers->num_editor;
    
    for (i = 0; i < num_editor; i++) {
	if (i == ed_id){
	    memmove(&ers->editor_record[i], &ers->editor_record[i+1], 
		    (num_editor - i - 1) * sizeof(ers->editor_record[i]));
	    ers->num_editor--;
	    i--;
	    num_editor--;
	}
    }
}

/*remove window cursor from editor cursor list */
int remove_editor_cursor (int ed_id, char *frame) {

    EDITOR_RECORD *er;
    cursors *c;
    
    er = editor_records->editor_record[ed_id];
    c = er->cursors;

    while (c != NULL) {
	if (c->cursor->frame_name != NULL) {
	    if (!strcmp (c->cursor->frame_name, frame)) {
		c = c->next;
		return 0;
	    }
	} else {   
	    /*if (c->cursor->posy == member_id) {
		printf ("get_editor_cursor: from renz \n");
		ce = c->cursor;
		return ce;
		}*/
	}
	c = c->next;
    }
    return 0;
}

void text_editor_shutdown (Tcl_Interp *interp, text_editor_result *result, int ed_id, int seq_num, char *frame) {

   char cmd[1024];
   cursor_e *c;
   int i, j, num_seq, seq_id, num_group, seq_nums;
  
   /* need to deregister sequence */
   num_group = editor_records->editor_record[ed_id]->num_group;
   for (j = 0; j < num_group; j++) {
       num_seq = editor_records->editor_record[ed_id]->seq_group[j]->nmembers;
       for (i = 1; i <= num_seq; i++) {
	   seq_id = GetEditorSeqId (ed_id, j, i);
	   seq_nums = GetEdenNum (seq_id);
	   editor_deregister (seq_nums, editor_callback, (text_editor_result *)result);
	   c = get_editor_cursor (ed_id, i, frame);
	   if (c != NULL) {
	       ed_delete_cursor (seq_nums, c->id, 1);
	   }
       }
   }
   /*remove editor cursor from list */
   /*remove_editor_cursor (ed_id, frame);*/

   editor_records->editor_record[ed_id]->text--;
   if (editor_records->editor_record[ed_id]->text == 0 
       && editor_records->editor_record[ed_id]->graphical == 0) {
       delete_editor (ed_id);
   }
   
   sprintf(cmd, "DeleteTextEditor %s", result->frame_name);
   if (TCL_ERROR == Tcl_Eval(interp, cmd)) {
       verror(ERR_WARN, "text editor ", "shutdown %s\n", 
	      interp->result);
   }
}

void static editor_callback (int seq_num, void *fdata, editor_reg_data *jdata) {

    text_editor_result *result = (text_editor_result *) fdata;

    switch(jdata->job) {
    case SEQ_CHANGED:
	{
	    
	    editor_redisplay (result->interp, seq_num, result);
	    /*printf ("text_editor:SEQ_CHANGED\n");*/
	    break;
	}
    case SEQ_CURSOR_NOTIFY:
	{	 
	    int posx = jdata->cursor_moved.cursor->abspos;
	    int posy = jdata->cursor_moved.cursor->posy;
	    
	    posx--;
	    if (result->cursor_id == jdata->cursor_moved.cursor->id) {
		editor_icursor_move (result->interp, result, posx, posy);
		/*printf ("text_editor:SEQ_CURSOR_NOTIFY\n");*/
	    }
	    break;
	}
    case SEQ_SELECTED:
	{
	    int sel_first, sel_last;

	    sel_first = jdata->selected.selection->sel_first;
	    sel_last = jdata->selected.selection->sel_last;
   
	    if (sel_first != 0 && sel_last != 0) {
		if (!strcmp (result->frame_name, jdata->selected.frame_name)) {
		    text_editor_draw_selection (result->interp, result, sel_first, sel_last);
		}
	    }
	    /*printf ("text_editor:SEQ_SELECTED\n");*/
	    break;
	}
    case TEXT_EDITOR_QUIT:
	{    
	    int ed_id = result->ed_id;
	    char *frame_name;

	    frame_name = jdata->exit.frame_name;
	    
	    if (!strcmp (result->frame_name, frame_name)) {
		text_editor_shutdown (result->interp, result, ed_id, seq_num, frame_name);
		/*printf ("text_editor:TEXT_QUIT\n");*/
	    }
	    break;
	} 
    }
}

out_canvas_e *init_output (void) {

    out_canvas_e *o;
    
    if (NULL == (o = (out_canvas_e *)xmalloc(sizeof(out_canvas_e))))
	return NULL;
    o->interp = NULL;
    o->cursor = NULL;
    o->cursor_visible = 0;

    return o;
}

editor_data *init_editor_reg_data (void) {

    editor_data *e;

    if (NULL == (e = (editor_data *)xmalloc(sizeof(editor_data))))
	return NULL;
    e->win_name = NULL;
    e->ed_id = 0;
    return e;
}

editor_cursor *init_editor_cursor (void) {

    editor_cursor *d;

    if (NULL == (d = (editor_cursor *)xmalloc(sizeof(editor_cursor))))
	return NULL;

    d->rid = 0;
    d->ed_id = 0;
    d->seq_id = 0;
    d->cursorPos = 0;
    d->cursor = NULL;
    d->prev_pos = 0;
    d->cursor_visible = 0;
    return d;
}

void set_editor_cursor (int ed_id, cursor_e *cursor) {
    
    EDITOR_RECORD *er;
    cursors *c;

    er = editor_records->editor_record[ed_id];

    c = init_cursors();
    c->cursor = cursor;
    c->next = er->cursors;
    er->cursors = c;
}

int SeqedRegister (ClientData clientData, 
		   Tcl_Interp *interp, 
		   int argc, 
		   char **argv) {

    char *frame_name, *ed_win;
    int ed_id, cur_pos;
    text_editor_result *result;
    seq_reg_cursor_notify cn;
    int i, id, num_seq, line_width = 1;
    int seq_num, seq_id;
    int j, num_group;
 
    if (argc != 5) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], "seq_ids\"", (char*)NULL);
        return TCL_ERROR;
    }
  
    frame_name = argv[1];
    ed_win = argv[2];
    ed_id = atoi (argv[3]);
    cur_pos = atoi (argv[4]);
  
    result = init_text_editor_result ();   
    strcpy (result->win_name, ed_win);
    strcpy (result->frame_name, frame_name);
    result->ed_id = ed_id; /*editor ID*/
    result->op_func = editor_callback;
    result->interp = interp;
    id = get_editor_reg_id();
    result->index = id; /*register ID*/

    /*for every sequence do registration using same id */
    num_group = editor_records->editor_record[ed_id]->num_group;
    for (j = 0; j < num_group; j++) {
	num_seq = editor_records->editor_record[ed_id]->seq_group[j]->nmembers;
	for (i = 1; i <= num_seq; i++) {
       	    editor_cursor *sd;
	    sd = init_editor_cursor();
	    seq_id = GetEditorSeqId (ed_id, j, i);
	    seq_num = GetEdenNum (seq_id);
	    sd->cursor = editor_create_cursor (seq_num, 1, NULL, line_width, 1, HORIZONTAL);
	    sd->cursor->editor_reg_id = id;
	    result->cursor_id = sd->cursor->id;
	    editor_register (seq_num, editor_callback, (void *)result, TEXT, id);/*FIXME */
	    sd->prev_pos = sd->cursor->abspos;
	    sd->cursor->abspos = cur_pos;

	    sd->cursor->frame_name = strdup (frame_name);
	    set_editor_cursor (ed_id, sd->cursor);
	   
	    cn.job = SEQ_CURSOR_NOTIFY;
	    cn.cursor = sd->cursor;
	    cn.cursor->job = CURSOR_MOVE;
	    editor_notify(seq_num, (editor_reg_data *)&cn);
	}
    }
    /*editor_redisplay (interp, ed_id, result);*/
    vTcl_SetResult(interp, "%d", id);
    return TCL_OK; 
}

char *get_editor_name (int id) {

    EDITOR_RECORD *er;
    int group_id = 0, member_id = 1; /*give first group and first member's name */

    er = editor_records->editor_record[id];     
    return er->seq_group[group_id]->members[member_id]->data->sequence->name;
}

int GetEditorName (ClientData clientData, 
		   Tcl_Interp *interp, 
		   int argc, char **argv) {

    int ed_id;
    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "ed_id\"", (char*)NULL);
        return TCL_ERROR;
    }
    ed_id = atoi (argv[1]);
    Tcl_ResetResult(interp);
    vTcl_SetResult(interp, "%s", get_editor_name(ed_id));
    return TCL_OK;
}

EDIT *create_insert_fragment_edit (EDITOR_RECORD *er, int member_id, SEQUENCE *fi, int pos, int cons) {

    EDIT *edit;

    edit = create_edit();
    edit->seq_id    = member_id;
    edit->position  = pos;
    edit->consensus = cons;
    edit->operation = INSERT_FRAGMENT;
    edit->sequence  = copy_sequence (fi);
    
    return edit;
}

int dis_seq ( SEQUENCE *sequence) {
    
    int i, num_entry;
    
    if(sequence->feature_table) {
	num_entry = sequence->feature_table->num_entry;
	for (i = 0; i < num_entry; i++) {
	    print_entry(sequence->feature_table->entry[i]);	  
	}	
    }
    return 0;
}

/* prepare to insert in string */
void make_hole_in_string (char *string, int position, int string_len) {

  /* paranoia: caller checks!!! make sure we keep 0 on end*/

  memmove((void *)&string[position+string_len],(void *)&string[position], strlen(string) - position + 1);

}

void convert_to_upper (char *string) {
    
    int i, len;

    len = strlen (string);

    for (i = 0; i <len; i++) {
	string[i] = toupper (string[i]);
    }
}

void convert_to_lower (char *string) {

    int i, len;

    len = strlen (string);

    for (i = 0; i <len; i++) {
	string[i] = tolower (string[i]);
    }
}

int copy_string_to_sequence (SEQUENCE *sequence, char *string, int position) {

    int string_len;
    
    string_len = strlen(string);
  
    if (NULL == (sequence->seq = (char *)xrealloc(sequence->seq, sequence->length + 1 + string_len + 1)))
	return -1;
 
    make_hole_in_string(sequence->seq, position, string_len);

    if (isupper (sequence->seq[0])) convert_to_upper (string);
    else convert_to_lower (string);
    strncpy (&(sequence->seq)[position], string, string_len);
    sequence->length += string_len;

    return 0;
}

int insert_string_in_sequence (SEQUENCE *sequence, 
			       char *string, 
			       int position, 
			       int imode,
			       FEATURE_TABLE **ft) {
    FEATURE_TABLE *dft = NULL;
    int err;

    if (sequence->length >= position) {	
	err = copy_string_to_sequence (sequence, string, position);
	if (err == -1) return -1;
	if (imode == 0) {
	   insert_string_extend_ft_range (sequence, string, position);
	}
	if (imode == 1) {
	    insert_string_break_ft_range (sequence, string, position);
	}
	if (imode == 2) {
	    dft = insert_string_delete_ft_range (sequence, position);
	}
    } 
    *ft = dft;
    return 0;
}

int insert_fragment_in_sequence (SEQUENCE *sequence, SEQUENCE *segment, int position, int imode) {

    FEATURE_TABLE *ft, *dft;
    int err;
   
    err = insert_string_in_sequence (sequence, segment->seq, position, imode, &dft);
    if (err == -1) return -1;

    /* add segment feature table to sequence feature table */
    if (segment->feature_table) {
	ft = copy_feature_table (segment->feature_table);
	change_feature_table_location (ft, position - 1);
	err = insert_sub_ft_in_feature_table (sequence, ft);
	if (err) return -1; 
    }
    /*if (dft != NULL) {
      FIXME: how to store this message in UNDO list? 
     }*/
    /*if (dft != NULL) {
      insert_sub_ft_in_feature_table (sequence, dft);	
      }*/
  return 0;
}

int insert_fragment_in_editor ( EDITOR_RECORD *er, EDIT *edit, int add_to_undo ) {

    int group_id = 0; /*fixme*/
    int err;
    
    err = insert_fragment_in_sequence (er->seq_group[group_id]->members[edit->seq_id]->data->sequence,
				     edit->sequence, edit->position, er->ft_imode);
    if (err) goto bail_out;
    if (add_to_undo) {
	err = add_edit(er->edits, edit);
	if (err) goto bail_out;
    }
 bail_out:
    return err;
}

int insert_update_sequence (int ed_id, int group_id, int member_id, int pos, SEQUENCE *sequence) {

    int i, num_seq;
    int cons, err, add_to_undo = 1;
    EDIT *edit;
    EDITOR_RECORD *er;

    er = editor_records->editor_record[ed_id];
    num_seq = er->seq_group[group_id]->nmembers;
    pos -- ;
    if (member_id == 0) {
	cons = 1;
	for (i = 0; i <= num_seq; i++) {         
	    edit = create_insert_fragment_edit (er, i, sequence, pos, cons);
	    err = insert_fragment_in_editor (er, edit, add_to_undo);
	    if (err != 0) return -1;
	}
    } else {
	cons = 0;
        edit = create_insert_fragment_edit (er, member_id, sequence, pos, cons);
	err = insert_fragment_in_editor (er, edit, add_to_undo);
	if (err != 0) return -1;
	err = update_consensus (er->seq_group[group_id], 1);
	if (err != 0) return -1;
    }
    return 0; 
}

EDIT *create_delete_fragment_edit (EDITOR_RECORD *er, 
				   int member_id, 
				   SEQUENCE *fd, 
				   int pos, 
				   int cons) {

    EDIT *edit;

    edit = create_edit();
    edit->seq_id    = member_id;
    edit->position  = pos;
    edit->consensus = cons;
    edit->operation = DELETE_FRAGMENT;
    edit->sequence  = copy_sequence (fd);
    
    return edit;
}

/*delete in string */
void fill_hole_in_string (char *string, int position, int string_len) {

  /* paranoia: caller checks!!! make sure we keep 0 on end*/
   
    memmove((void *)&string[position],(void *)&string[position+string_len], 
	    strlen(string) - position - string_len + 1);

}

FEATURE_TABLE *delete_fragment_in_sequence (SEQUENCE *sequence, SEQUENCE *segment, int position) {

    int string_len;
    FEATURE_TABLE *ft = NULL;
 
    string_len = strlen (segment->seq);
    
    if (sequence->length >= position) {
	ft = get_feature_table (sequence, string_len, position);
	fill_hole_in_string(sequence->seq, position, string_len);
	sequence->length -= string_len;
	sequence->end = sequence->length;
    }
    if (sequence->feature_table != NULL)
      ft = delete_string_modify_feature_table (sequence, string_len, position);

    if (ft == NULL) return NULL;
    if (ft) {
	change_feature_table_location (ft, -(position - 1) );
    } 
    return ft;
}

int delete_fragment_in_editor (EDITOR_RECORD *er, EDIT *edit, int add_to_undo) {

    FEATURE_TABLE *ft;
    int group_id = 0; /*Fixme: get this value from edit */
    int err;

    ft = delete_fragment_in_sequence (er->seq_group[group_id]->members[edit->seq_id]->data->sequence,
				     edit->sequence, edit->position);
 
    if (ft && edit->seq_id != 0) {
	err = add_ft_in_sequence (edit->sequence, ft);
	if (err == -1) goto bail_out;
    }
    if (add_to_undo) {
	err = add_edit(er->edits, edit);
	if (err) goto bail_out;
    }
    return 0;
 bail_out:
    if (ft) free_feature_table (ft);
    return err;
}


int delete_update_sequence (int ed_id, int group_id, int member_id, int pos, SEQUENCE *sequence) {

    int i, num_seq;
    int cons, err, add_to_undo = 1;
    EDIT *edit;
    EDITOR_RECORD *er;

    er = editor_records->editor_record[ed_id];
    num_seq = er->seq_group[group_id]->nmembers;
    
    if (member_id == 0) {
	cons = 1;
	for (i = 0; i <= num_seq; i++) {
	    edit = create_delete_fragment_edit (er, i, sequence, pos, cons);
	    err = delete_fragment_in_editor (er, edit, add_to_undo);
	    if (err != 0) return -1;
	}
    } else {
	cons = 0;
	edit = create_delete_fragment_edit (er, member_id, sequence, pos, cons);
	err = delete_fragment_in_editor (er, edit, add_to_undo);
	if (err != 0) return -1;
	err = update_consensus (er->seq_group[group_id], 1);
	if (err != 0) return -1; 
    }
    
    return 0;
}

EDIT *create_change_fragment_edit (EDITOR_RECORD *er, int member_id, SEQUENCE *s_with, int pos, int cons) {
    EDIT *edit;
    
    edit = create_edit();
    edit->seq_id    = member_id;
    edit->position  = pos;
    edit->consensus = cons;
    edit->operation = CHANGE_FRAGMENT;
    edit->sequence  = copy_sequence (s_with);
    return edit;
}
/* no longer used */
int change_fragment_in_editor (EDITOR_RECORD *er, EDIT *edit, int add_to_undo) {

    int err;
    FEATURE_TABLE *ft;
    
    if (add_to_undo) {
	err = add_edit(er->edits, edit);
	if (err) goto bail_out;
    }
    return 0;
 bail_out:
    if (ft) free_feature_table (ft);
    return err;
}

int replace_update_sequence (int ed_id, 
			     int group_id, 
			     int member_id, 
			     int pos, 
			     SEQUENCE *s_replace, 
			     SEQUENCE *s_with) 
{
    int i, num_seq, cons;
    int err, add_to_undo = 1;
    EDIT *edit;
    EDITOR_RECORD *er;

    er = editor_records->editor_record[ed_id];   
    num_seq = er->seq_group[group_id]->nmembers;
  
    if (member_id == 0) {
	cons = 1;
	for (i = 0; i <= num_seq; i++) {
	    edit = create_change_fragment_edit (er, i, s_replace, pos, cons);    
	    err = delete_fragment_in_editor (er, edit, add_to_undo);
	    if (err != 0) return -1;
	    edit = create_change_fragment_edit (er, i, s_with, pos, cons);
	    err = insert_fragment_in_editor (er, edit, add_to_undo);	    
	    if (err != 0) return -1;
	}
    } else {
	cons = 0;
	edit = create_delete_fragment_edit (er, i, s_replace, pos, cons);    
	err = delete_fragment_in_editor (er, edit, add_to_undo);
	if (err != 0) return -1;
	edit = create_insert_fragment_edit (er, i, s_with, pos, cons);
	err = insert_fragment_in_editor (er, edit, add_to_undo);
	if (err != 0) return -1;
	err = update_consensus (er->seq_group[group_id], 1);
    }
    return 0;
}

int copy_update_buffer (int ed_id, int group_id, int member_id, int pos, SEQUENCE *sequence) {

    SEQUENCE *s;
    FEATURE_TABLE *ft;
    SITE_HANG *sh_left; /* added for merage text editor and graphical editor */
    SITE_HANG *sh_right;  
    int len;

    s = editor_records->editor_record[ed_id]->seq_group[group_id]->members[member_id]->data->sequence;
    len = strlen (sequence->seq);
    sequence->start = 1;
    sequence->end = len;
    sequence->length = len;
    sequence->id = -1; /* means clipboard */
    if (len > 1) {
	ft = get_feature_table (s, len, pos);
	if (ft) {
	    change_feature_table_location (ft, -(pos - 1) );
	} 
	sequence->feature_table = ft;
    }
    sh_left = init_site_hang ();
    sh_right = init_site_hang ();
    sequence->left_end = sh_left;
    sequence->right_end = sh_right;
    if (editor_records->clipboard != NULL) 
	free_sequence (editor_records->clipboard);
    editor_records->clipboard = sequence;
   
    return 0;
}

SEQUENCE *get_editor_buffer (void) {

    return editor_records->clipboard;
}

void set_editor_buffer (SEQUENCE *s) {

    editor_records->clipboard = s;
}
/*
 * deletes all results on a single widget
 */
int TextEditorResultUpdate (ClientData clientData, 
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

    if (strcmp(args.option, "QUIT") == 0) {
	info.job = TEXT_EDITOR_QUIT;
    } else {
	verror(ERR_FATAL, "text_editor_result_notify", "invalid command");
	return TCL_OK;
    }

    if (args.id == -1) {
	verror(ERR_FATAL, "text_editor_result_notify", "invalid command");
	return TCL_OK;
    } else {
	editor_result_notify(args.id, (editor_reg_data *)&info, 1);
    }
    return TCL_OK;
}

int SetFtInsertionMode (ClientData clientData,
			Tcl_Interp *interp, 
			int argc, 
			char **argv)

{
    int ed_id;
    char *selection;

    if (argc != 3) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], "editor_id\"", "selection\"", (char*)NULL);
	return TCL_ERROR;
    }
    
    ed_id = atoi (argv[1]);
    selection = argv[2];

    if (!strcmp (selection, "extend_feature")) {
	editor_records->editor_record[ed_id]->ft_imode = 0;
    }
    if (!strcmp (selection, "break_feature")) {
	editor_records->editor_record[ed_id]->ft_imode = 1;
    }
    if (!strcmp (selection, "delete_feature")) {
	editor_records->editor_record[ed_id]->ft_imode = 2;
    }
   return TCL_OK; 
}

int GetEditorSelection (ClientData clientData,
			Tcl_Interp *interp, 
			int argc, 
			char **argv)
{
    EDITOR_RECORD *er;
    int seq_id;
    int ed_id, group_id, member_id, sel_first, sel_last;

    if (argc != 2) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], "seq_id\"",(char*)NULL);
	return TCL_ERROR;
    }
  
    seq_id = GetSequenceIdByName (argv[1]);   
    if (editor_records == NULL) {
	vTcl_SetResult(interp, "%d", 0);
    } else {
	if (editor_records->selection == NULL) {
	    vTcl_SetResult(interp, "%d", 0);
	}
	if (editor_records->selection != NULL) {
	    ed_id = editor_records->selection->ed_id;
	    member_id = editor_records->selection->member_id;
	    sel_first = editor_records->selection->sel_first;
	    sel_last = editor_records->selection->sel_last;
	    er = editor_records->editor_record[ed_id];
	    group_id = 0; /*FIXME*/
	    if (seq_id == er->seq_group[group_id]->members[member_id]->data->sequence->id) {
	  	vTcl_SetResult(interp, "%d %d", sel_first, sel_last);  
	    } else {
		vTcl_SetResult(interp, "%d", 0);
	    }
	}
    }
    return TCL_OK; 
}

void set_editor_selection (SELECTION *sel) {

    if (editor_records->selection != NULL)
	free_selection (editor_records->selection);
    editor_records->selection = sel;
}

int SpinEditorExit (ClientData clientData,
		    Tcl_Interp *interp, 
		    int argc, 
		    char **argv)
{
    /*free_editor_registration_list*/
    free_sequences (eden);
    /*free_spin_sequence_list*/ 
    free_sequences (sequences);

    return TCL_OK;
}

char *get_qualifier (char *qualifiers, int begin, int *next) {

    int i, len, ll, b;
    char *qua;
    int nonspace = 0;

    len = strlen (qualifiers);
    b = begin;
    for (i = begin; i <= len; i++) {
	if (qualifiers[i] && !isspace(qualifiers[i]))
	    nonspace=1;

	if (nonspace &&
	    (qualifiers[i] == '\n' ||
	     (qualifiers[i] == '\0' && qualifiers[i-1] != '\n')
	     )) {
	    ll = i - b + 1;
	    if (NULL == (qua = (char *)xmalloc((ll + 1)*sizeof(char) )))
		return NULL;
	    strncpy (qua, &qualifiers[b], ll);
	    if (qualifiers[i] == '\n')
		ll--;
	    qua[ll] = 0;
	    *next = i+1;
	    return qua;
	}
    }
    *next = 0;
    return NULL;
}

void editor_file_save (FILE *pw, SEQUENCE *s) {

    FEATURE_TABLE *ft;
    char *seq;
    int num_entry, start, end;
    int line_len = 60;
    int line_seg = 10;
    int i, j;
    
    ft = s->feature_table;
    seq = s->seq;
   
    num_entry = ft->num_entry;
    start = 1;
    end = s->length;
 
    fprintf(pw,"ID   %s\n", s->name);
   
    /*start to write feature table*/
    if (ft != NULL) {
	for (i = 0; i < num_entry; i++) {
	    char *qua, *locc;
	    int begin, next;

	    fprintf (pw,"FT   %-16s", ft->entry[i]->type);
	    /* to write location */
	    begin = 0;
	    locc = get_qualifier (ft->entry[i]->location, begin, &next);
	    fprintf (pw, "%s\n", locc);
	    while (begin < end && locc != NULL) {
		begin = next;
		locc = get_qualifier (ft->entry[i]->location, begin, &next);
		if (locc != NULL) {
		    fprintf (pw, "FT   %s\n", locc);
		}
	    }	 
	    /* to write qualifier */
	    begin = 0;
	    qua = get_qualifier (ft->entry[i]->qualifiers, begin, &next);
	    while (begin < end && qua != NULL) {
		fprintf (pw, "FT                   %s\n", qua);
		begin = next;
		qua = get_qualifier (ft->entry[i]->qualifiers, begin, &next);
	    }	    
	}
    }
    /*start to write sequence*/
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
	fprintf (pw, "%c", ' ');
    fprintf (pw, "%10d\n", end - start + 1);
    fprintf (pw, "//\n"); 
    fclose (pw);
}

int EditorSaveFile (ClientData clientData,
		    Tcl_Interp *interp, 
		    int argc, 
		    char **argv)
{
    char *file_in, *file_out;
    SEQUENCE *s;
    FILE *fp;
    int seq_id;
    
    if (argc != 3) {
	Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], "file_name_in\"", "file_name_out\"", 
			 (char*)NULL);
	return TCL_ERROR;
    }
    file_in = argv[1];
    file_out = argv[2];

    if (NULL == (fp = fopen(file_out, "w"))) {
	Tcl_AppendResult (interp, "Unable to save sequence\"", (char*)NULL);
	return TCL_OK;
    }

    seq_id = GetSequenceIdByName (file_in);
    s = GetSequencesSequence (seq_id);
  
    editor_file_save (fp, s);

    return TCL_OK; 
}

int Editor_init(Tcl_Interp *interp) {
 
    Tcl_CreateCommand(interp, "get_seq_max_len", GetSeqMaxLen, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_num_seqs_in_editor", GetNumSeqsInEditor, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_editor_reg_num", GetEditorRegNum, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "add_editor", AddEditor, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_editor_name", GetEditorName, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);    
    Tcl_CreateCommand(interp, "seqed_register", SeqedRegister, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "text_editor_result_update", TextEditorResultUpdate, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "set_ft_insertion_mode", SetFtInsertionMode, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL); 
    Tcl_CreateCommand(interp, "get_qualifier", tcl_get_qual_array, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_keyword", tcl_get_key_array, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_editor_selection", GetEditorSelection, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);   
    Tcl_CreateCommand(interp, "spin_editor_exit", SpinEditorExit, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "editor_save_file", EditorSaveFile, 
		      (ClientData) NULL, 
		      (Tcl_CmdDeleteProc *) NULL);
    return TCL_OK;
}

int Seqed_Init(Tcl_Interp *interp) {

  Read_sequence_Init (interp);
  REnzyme_box_Init (interp);
  REnzyme_search_Init (interp);
  RenzymeCmds_Init (interp);
  RenzCanvas_Init (interp);
  TextEditor_Init (interp);
  Editor_init (interp);
  Hang_editor_Init (interp);
  FeatureEditor_Init (interp);
  spin_init_globals (interp);
  nip_init_globals (interp);
  
  return TCL_OK;
}
