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

#include "misc.h"
#include "parse_feature.h"
#include "feature_colour.h"
#include "feature_table.h"

FEATURE_TABLE *init_feature_table ( void ) {

  FEATURE_TABLE *ft;
  ft_entry  **e;

  if ( ( NULL == ( ft = (FEATURE_TABLE* ) xmalloc ( sizeof (FEATURE_TABLE) )))) 
    return NULL; 
  if ( ( NULL == ( e = (ft_entry** ) xmalloc ( sizeof (ft_entry *) )))) return NULL;  
  ft->entry = e;
  ft->num_entry = 0;
  return ft;
}

FEATURE_TABLE *realloc_feature_table ( FEATURE_TABLE *ft, int num_entry ) {

  ft_entry  **e;
  e = ft->entry;

  if ((NULL == ( e = 
		(ft_entry **)xrealloc (e, sizeof (ft_entry *) * (num_entry))))) 
      return NULL;
  ft->entry = e;
  
  ft->num_entry = num_entry;
  return ft;
}

void free_feature_table (FEATURE_TABLE *ft) {
    
    int i, num_entry;

    if (!ft) return;
    num_entry = ft->num_entry;

    for (i = 0; i < num_entry; i++) {
	del_ft_entry (ft->entry[i]);     
    }
    free(ft);
}

ft_range * copy_ft_range (ft_range *r_sou) {

    ft_range *r, *r_des = NULL, *r_tmp;
    for (r = r_sou; r; r = r->next) {
	ft_range *tmp;

	if (NULL == (tmp = new_ft_range()))
	    goto error;

	if (NULL == (tmp->left = new_ft_location()))
	    goto error;	   
	tmp->left->min = r->left->min;
	
	tmp->left->max = r->left->max;

	tmp->left->min_lt = r->left->min_lt;
	
	tmp->left->max_lt = r->left->max_lt;

	tmp->left->type = r->left->type;
        
	if (r->right != NULL) {
	    if (NULL == (tmp->right = new_ft_location()))
		goto error;
	    tmp->right->min = r->right->min;
	    tmp->right->max = r->right->max; 
	    tmp->right->min_lt = r->right->min_lt;
	    tmp->right->max_lt = r->right->max_lt;
	    ; 
	    tmp->right->type = r->right->type;		
	}
	tmp->next = NULL;

	if (!r_des) { 
	    r_des = tmp; /* first */   
	} else {
	    r_tmp = r_des;
	    while (r_tmp->next) {
		r_tmp = r_tmp->next;
	    }
	    r_tmp->next = tmp;	  
	}
	r_des->complemented = r_sou->complemented;
    }
    return r_des;
 error:
    return NULL;
}

ft_entry *copy_ft_entry ( ft_entry *entry) {

    ft_entry *e;
    ft_range *r, *r_tmp;
    int loc_l, qua_l;
    
    if (NULL == (e = new_ft_entry()))
	goto error;
    
    strcpy (e->type, entry->type);
    /*printf ("entry->location=%s\n", entry->location);*/
    loc_l = strlen(entry->location);
    if (NULL == (e->location = (char *)xmalloc((loc_l + 1)*sizeof(char))))
	    goto error;	 
    strcpy (e->location, entry->location);

    if (entry->qualifiers != NULL) {
	qua_l = strlen(entry->qualifiers);
	if (NULL == (e->qualifiers = (char *)xmalloc((qua_l + 1)*sizeof(char))))
	    goto error;
	strcpy (e->qualifiers, entry->qualifiers);
    }

    e->qual_hash_init = 0;
    init_ft_qual_hash (e, e->qualifiers);

     for (r = entry->range; r; r = r->next) {
	ft_range *tmp;

	if (NULL == (tmp = new_ft_range()))
	    goto error;
	if (NULL == (tmp->left = new_ft_location()))
	    goto error;	   
	tmp->left->min = r->left->min;
	tmp->left->max = r->left->max;
	tmp->left->min_lt = r->left->min_lt;
	tmp->left->max_lt = r->left->max_lt;
	tmp->left->type = r->left->type;
	if (r->right != NULL) {
	    if (NULL == (tmp->right = new_ft_location()))
		goto error;
	    tmp->right->min = r->right->min;
	    tmp->right->max = r->right->max;
	    tmp->right->min_lt = r->right->min_lt;
	    tmp->right->max_lt = r->right->max_lt;
	    tmp->right->type = r->right->type;		
	}

	tmp->next = NULL;
	if (!e->range) { 
	    e->range = tmp;	   
	} else {
	    r_tmp = e->range;
	    while (r_tmp->next) {
		r_tmp = r_tmp->next;
	    }	    
	    r_tmp->next = tmp;	  
	}
	e->range->complemented = r->complemented;
	}  
    /*e->range = copy_ft_range (entry->range);*/
    /* print_entry (e);*/
    return e;
 error:
    if (e) del_ft_entry (e);
    return NULL;
}  

FEATURE_TABLE *copy_feature_table ( FEATURE_TABLE *ft ) {

    FEATURE_TABLE *ft_new;
    int i, num_entry;
    
    ft_new = init_feature_table();
    num_entry = ft->num_entry;

    for ( i = 0; i < num_entry; i++) {
	if ( ( NULL == ( ft_new->entry = 
			 (ft_entry **)xrealloc (ft_new->entry, sizeof (ft_entry *) * (i + 1))))) 
	    goto error;
	ft_new->entry[i] = copy_ft_entry (ft->entry[i]);
    } 
     ft_new->num_entry = num_entry;
    return ft_new;
 error:
    if (ft_new)
	free_feature_table (ft_new);
    return NULL;
}

int extend_deleted_feature_table (FEATURE_TABLE *ft, ft_entry *entry, int entry_id) {
   
    if (entry_id != 1) {
	if ((NULL == ( ft->entry = 
		   (ft_entry **)xrealloc (ft->entry, sizeof (ft_entry*) * (entry_id)))))
	    return -1;
     }	    
    ft->entry[entry_id - 1] = copy_ft_entry (entry);
    ft->num_entry = entry_id; 
    return 0;
}

/* return  1:  r1 great than r2. */ 
/* return -1:  r1 less than r2.  */
/* return  0:  r1 equal r2.      */
/* return  2:  no use at moment  */

int compare_ft_single_range (ft_range *r1, ft_range *r2) {

    if (r1->left->min > r2->left->min) return 1;
    else if (r1->left->min < r2->left->min) return -1;
    else {
	if (r1->left->type == BASE && r2->left->type == BASE 
	    || r1->left->type == SIT && r2->left->type == SIT)
	    if (r1->left->max != r2->left->max) return 2; 
	if (r1->right && r2->right){
	    if (r1->right->min != r2->right->min) return 2;
	    if (r1->right->type == BASE && r2->right->type == BASE 
		|| r1->right->type == SIT && r2->right->type == SIT)
		if (r1->right->max != r2->right->max) return 2;
	}
    } 
    return 0;
}

void set_single_range_null (ft_range *r) {
    
    r->left->min = -1;
    if (r->left->type == BASE || r->left->type == SIT ) 
	r->left->max = -1;
    if (r->right) {
	r->right->min = -1;
	if (r->right->type == BASE || r->right->type == SIT)
	    r->right->max = -1;
    }
}

void copy_single_ft_range (ft_range *r_des, ft_range *r_sou) {
    
    r_des->left->min = r_sou->left->min;
    if (r_sou->left->type == BASE || r_sou->left->type == SIT ) { 
	r_des->left->type = r_sou->left->type;
	r_des->left->max = r_sou->left->max;
    }
    if (r_sou->right) {
	r_des->right->min = r_sou->right->min;
	if (r_sou->right->type == BASE || r_sou->right->type == SIT) {
	    r_des->right->type = r_sou->right->type;
	    r_des->right->max = r_sou->right->max;
	}
    }
}

void change_single_ft_range (ft_range *r, int string_len) {

     r->left->min -= string_len;
     if(r->left->type == BASE || r->left->type == SIT)
	 r->left->max -= string_len;
     if(r->right != NULL) {
	 r->right->min -= string_len;
	 if(r->right->type == BASE || r->right->type == SIT) 
	     r->right->max -= string_len;
     }
}

void change_single_ft_range_insert (ft_range *r, int string_len) {

     r->left->min += string_len;
     if(r->left->type == BASE || r->left->type == SIT)
	 r->left->max += string_len;
     if(r->right != NULL) {
	 r->right->min += string_len;
	 if(r->right->type == BASE || r->right->type == SIT) 
	     r->right->max += string_len;
     }
}
/* return 1: there is no "-1" or max same as min in the linked range, return 0: others */
int check_ft_range (ft_range *r) {

    ft_range *r_tmp;

    for (r_tmp = r; r_tmp; r_tmp = r_tmp->next) {
	if (r_tmp->left->min == -1) return 0;
	if (r_tmp->left->type == BASE || r_tmp->left->type == SIT) {
	    if (r_tmp->left->max == -1) return 0;
	    if (r_tmp->left->min == r_tmp->left->max) return 0;
	}
	if (r->right) {
	    if (r->right->min == -1) return 0;
	    if (r->right->type == BASE || r->right->type == SIT) {
		if (r->right->max == -1) return 0;
		if (r_tmp->right->min == r_tmp->right->max) return 0;
	    }
	}
    }
    return 1;
}

int add_to_deleted_feature_table (FEATURE_TABLE *ft, ft_entry *e, ft_entry *e_start, int num_del) {
    
    int err;

    err = check_ft_range (e->range);
    if (!err) {
	num_del ++;		
	err = extend_deleted_feature_table (ft, e_start, num_del);
	if (err == -1) return -1; 
    }
    return num_del;
}

void mv_right_to_left (ft_range *r) {

  if (r->right) {
    r->left->min = r->right->min;
    if (r->right->type == BASE || r->right->type == SIT) {
    r->left->max = r->right->max;
    }
    r->left->type = r->right->type; 
    r->right = NULL;
  }
}

/** 
 *  got a new feature_key and add it to feature_table
 *  the sequence will not be changed
 */

int add_feature_to_feature_table (SEQUENCE *sequence, ft_entry *entry) {

    FEATURE_TABLE *ft;
    int num_entry;

    ft = sequence->feature_table;
    if (ft == NULL) ft = init_feature_table ();
    sequence->feature_table = ft;

    num_entry = ft->num_entry;
    num_entry ++;
  
    if (NULL == realloc_feature_table (ft, num_entry ))
	return -1;;
    ft->entry[num_entry - 1] = entry;
    
    return 0;
}

int add_ft_in_sequence (SEQUENCE *s, FEATURE_TABLE *ft) {

    int i, nf;
    
    nf = ft->num_entry;

    for (i = 0; i < nf; i++) {
	if (-1 == add_feature_to_feature_table (s, ft->entry[i]))
	    return -1;;
    }
    return 0;
}

int check_single_ft_range (ft_range *r) {
    
    if (r->left->min != -1) return 0;
    if (r->left->type == BASE || r->left->type == SIT) {
	if (r->left->max != -1) 
	    return 0;
    }
    if (r->right) {
	if (r->right->min != -1) return 0;
	if (r->right->type == BASE || r->right->type == SIT) {
	    if (r->right->max != -1) 
	    return 0;
	}
    }
    return 1;
}

int compare_ft_entry (ft_entry *e1, ft_entry *e2) {

    if ( strcmp (e1->type, e2->type)) return 1;
    if ( strcmp (e1->location, e2->location) )return 1;
    if ( strcmp (e1->qualifiers, e2->qualifiers ) ) return 1;
    return 0;
}

int add_entry_to_feature_table (FEATURE_TABLE *ft, ft_entry *e) {
    
    int num_delete;

    num_delete = ft->num_entry;
    num_delete++;
    if (num_delete != 1)
	ft = realloc_feature_table (ft, num_delete); 
    ft->entry[num_delete - 1] = copy_ft_entry (e);
    ft->num_entry = num_delete;
    
    return 0;
}

int insert_sub_ft_in_feature_table ( SEQUENCE *seq_des, FEATURE_TABLE *ft_sou) {

    FEATURE_TABLE *ft_des;
     ft_entry *e_sou;
    int num_entry_des, num_entry_sou;
    int i;
    
    ft_des = seq_des->feature_table;
    num_entry_des = ft_des->num_entry;
    num_entry_sou = ft_sou->num_entry;
    for (i = 0; i < ft_sou->num_entry; i++) {
	e_sou = ft_sou->entry[i];
	add_entry_to_feature_table (ft_des, e_sou);
    }
    return 0;
}

void change_feature_table_location (FEATURE_TABLE *ft, int position) {

    int i, num_entry;
    ft_entry *e;
    ft_range *r;
    
    num_entry = ft->num_entry;
    for (i = 0; i < num_entry; i++) {
	e = ft->entry[i];
	for (r = e->range; r; r = r->next) {
	    r->left->min += position;
	    if (r->left->type == BASE || r->left->type == SIT)
		r->left->max += position;
	    else r->left->max = r->left->min;
	    if (r->right){
		r->right->min += position;
	        if (r->right->type == BASE || r->right->type == SIT)
		    r->right->max += position;
		else r->right->max = r->right->min;
	    }
	}
    }
}

ft_entry *get_ft_range_start (ft_entry *e, ft_range *r) {

  ft_entry * e_del;

  e_del = copy_ft_entry (e); 
  del_ft_range (e_del->range);   /* delete original one */
  e_del->range = copy_ft_range (r); /* new range start from r */
  return e_del;
}

void get_ft_range_Aa_start (ft_range *r, int position) {
    
    r->left->min = position;
    r->next = NULL;
}

void get_ft_range_Ab_end (ft_range *r, int del_end_pos) {

  r->left->max = del_end_pos;
  r->right = NULL;
  r->next = NULL; 
}

void get_ft_range_Ba_start (ft_range *r, int position) {
    r->left->min = position;
}

void get_ft_range_Ba_end (ft_range *r) {

  r->right = NULL;
  r->next = NULL;

}

void get_ft_range_Bb_end (ft_range *r, int del_end_pos) {

  r->right->max = del_end_pos;
  r->next = NULL;

}

void get_ft_range_Bb_start ( ft_range *r, int position) {

  mv_right_to_left (r);
  r->left->min = position;
}

void get_ft_range_middle (ft_range *r, int position, int del_end_pos) {

  r->left->min = position;
  r->right->min = del_end_pos;
  r->left->type = EXACT;
  r->right->type = EXACT;
  r->next =  NULL;
}





/*****************************************************************/
/*      A .. B           */
/*      a .. a     EXACT */
/*     a.b..a.b    BASE */
/*     a^b..a^b    SIT */ 
     
/******************************************************************/

/* extend the range of the feature when insertion is made inside the range */    
void insert_string_extend_ft_range (SEQUENCE *sequence, char *string, int position ) {

    int i, string_len, num_entry;
    ft_entry *e;
    ft_range *r, *rr;

    string_len = strlen(string);
   
    if(sequence->feature_table) {
	num_entry = sequence->feature_table->num_entry;
	for (i = 0; i < num_entry; i++) {
	    e = sequence->feature_table->entry[i];
	    for (r = e->range; r; r = r->next) {
		if (position < r->left->min) {
		    /*(a.b a.b) (a.b a.b) (a.b a.b)*/
		    /*         ^                   */
		    for (rr = r; rr; rr = rr->next) {
		    change_single_ft_range_insert (rr, string_len);
		}
		break;
		} else if ((r->left->type == BASE || r->left->type == SIT) && position < r->left->max) {
		    /*(a.b a.b)(a.b a.b) (a.b a.b)*/
		    /*           ^                */
		    for (rr = r; rr; rr = rr->next) {
			change_single_ft_range_insert (rr, string_len);
		    } 
		    r->left->min -= string_len;
		    break;
		} else if (r->right && position < r->right->min) {
		    /*(a.b a.b)(a.b a.b) (a.b a.b)*/
		    /*             ^              */		
		    for (rr = r; rr; rr = rr->next) {
			change_single_ft_range_insert (rr, string_len);
		    } 
		    r->left->min -= string_len;
		    if (r->left->type != EXACT) r->left->max -= string_len;
		    break;
		} else if (r->right && (r->right->type == BASE || r->right->type == SIT) 
			   && position < r->right->max) {
		    /*(a.b a.b)(a.b a.b) (a.b a.b)*/
		    /*               ^            */
		    
		    for (rr = r; rr; rr = rr->next) {
			change_single_ft_range_insert (rr, string_len);
		    }
		    r->left->min -= string_len;
		    if (r->left->type != EXACT) r->left->max -= string_len;
		    r->right->min -= string_len;
		    break;
		}
	    }
	}
    }
}

/* to break feature range into two parts when insertion is made in the feature range, */
/* and do not expand the range of the feature range, LINK??? YES!*/

void insert_string_break_ft_range (SEQUENCE *s, char *string, int position) {

    int frag_len;
    int i, num_entry;    
    ft_entry *e;
    ft_range *r, *rr;

    frag_len = strlen (string);
    num_entry = s->feature_table->num_entry;

    for (i = 0; i < num_entry; i++) {
	e = s->feature_table->entry[i];
	for (r = e->range; r; r = r->next) {
	    if (position < r->left->min) {
		/*(a.b a.b) (a.b a.b) (a.b a.b)*/
		/*         ^                   */
		for (rr = r; rr; rr = rr->next) {
		    change_single_ft_range_insert (rr, frag_len);
		}
		break;	 
	    } else if ((r->left->type == BASE || r->left->type == SIT) && position < r->left->max) {
		/*(a.b a.b)(a.b a.b) (a.b a.b)*/
		/*           ^                */
		ft_range *r_new;
		r_new = copy_ft_range (r);
		r_new->left->min = position + 1;
		for (rr = r_new; rr; rr = rr->next) {
		    change_single_ft_range_insert (rr, frag_len);
		} 
		r->left->max = position;
		r->right = NULL;
		r->next = r_new;
		break;
	    } else if (r->right && position < r->right->min) {
		/*(a.b a.b)(a.b a.b) (a.b a.b)*/
		/*             ^              */		
		ft_range *r_new;
		r_new = copy_ft_range (r);    
		r_new->left->min = position + 1;
		for (rr = r_new; rr; rr = rr->next) {
		    change_single_ft_range_insert (rr, frag_len);
		} 
		r->right->min = position;
		r->next = r_new;			
		break;
	    } else if (r->right && (r->right->type == BASE || r->right->type == SIT) 
		       && position < r->right->max) {
		/* (a.b a.b)(a.b a.b) (a.b a.b) */
		/*                ^             */
		ft_range *r_new;
		r_new = copy_ft_range (r);
		r_new->left->min = position + 1;
		r_new->left->max = r_new->right->max;
		r_new->right = NULL;
		for (rr = r_new; rr; rr = rr->next) {
		    change_single_ft_range_insert (rr, frag_len);
		} 
		r->right->max = position;
		r->next = r_new;
		break;
	    }
	    /*break;*/
	}
    }
}

void delete_entry_from_feature_table (FEATURE_TABLE *ft, int entry_id) {

    int i, num_entry;

    num_entry = ft->num_entry;
    for (i = 0; i < num_entry; i++) {
	if ( i == entry_id) {   
	    memmove(&ft->entry[i], &ft->entry[i+1], (num_entry - i - 1) * sizeof(ft->entry[i]));
	    ft->num_entry--;
	}
    }
}

FEATURE_TABLE *insert_string_delete_ft_range (SEQUENCE *sequence, int position) {

    int i, num_entry;    
    ft_entry *e;
    ft_range *r;
    FEATURE_TABLE *ft;
    int left, right;

    ft = init_feature_table ();
    num_entry = sequence->feature_table->num_entry;
    for (i = 0; i < num_entry; i++) {
	e = sequence->feature_table->entry[i];
	left = e->range->left->min;
	for (r = e->range; r; r = r->next) {
	    if (r->next == NULL) {
		if (r->right) {
		    if (r->right->type != EXACT)
			right = r->right->max;
		    else right = r->right->min;
		} else {
		   if (r->left->type != EXACT) 
		       right = r->left->max;
		   else right = r->left->min;
		}	
	    }
	}
	if (position >= left && position <= right) {
	   add_entry_to_feature_table (ft, e);
	   delete_entry_from_feature_table (sequence->feature_table, i);
	   i--;
	   num_entry--;
	}
    }
    return ft;
}

int get_ft_range_check (ft_entry *e, int start_pos, int end_pos) {

    ft_range *r;
    int r_min, r_max;

    r_min = e->range->left->min;
    for (r = e->range; r; r = r->next){
	    if (r->next == NULL){
		if (r->right) {
		    if (r->right->type != EXACT)
			r_max = r->right->max;
		    r_max = r->right->min;
		} else {
		   if (r->left->type != EXACT) 
		       r_max = r->left->max;
		   r_max = r->left->min;
		}
	    }
    }
    if (start_pos <= r_min && end_pos >= r_max) {
	return 1;
    }
    return 0;
}

int add_feature_source (FEATURE_TABLE *ft, int start, int end, char *name){

    ft_entry *e;
    ft_range *rg;
    ft_location *l;
    ft_location *r;
    char clone[100];
    char location[100];   
    int num_entry;

    
    if (NULL == (l = new_ft_location()))
	goto err;
    l->min = start+1;
    l->max = start+1;
    l->type = EXACT;
    if (NULL == (r = new_ft_location()))
	goto err;
    r->min = end+1;
    r->max = end+1;
    r->type = EXACT;

    if (NULL == (rg = new_ft_range()))
	goto err;
    rg->left = l;
    rg->right = r;
    rg->complemented = 0;
    rg->next = NULL;

    if (NULL == (e = new_ft_entry()))
	goto err;

    sprintf (e->type, "source");
    e->range = rg;
    sprintf (clone, "clone=\"%s:%d..%d\"", name, start, end);
    e->qualifiers = strdup (clone);
    sprintf(location, "%d..%d", start, end);
    e->location = strdup (location);
    /*FIXME: store qual_hash?*/
    init_ft_qual_hash (e, e->qualifiers);

    num_entry = ft->num_entry;
    num_entry ++;
    realloc_feature_table (ft, num_entry);
    ft->entry[num_entry - 1] = e;

    return 0;
 err:
    if (l) del_ft_location (l);
    if (r) del_ft_location (r);
    if (rg) del_ft_range (rg);
    if (e) del_ft_entry (e);
    return -1;
}


/* To create fragment from sequence range of (start..end), if the feature key inside
   this range, then fragment will have this feature. BUT will always have source key
   with /clone=name:ID, start..end */
 
FEATURE_TABLE *get_feature_table (SEQUENCE *sequence, int string_len, int position) {

  FEATURE_TABLE *ft;
  ft_entry *e, *e_del;
  ft_entry *e_start = NULL; 
  ft_range *r_start, *r_end, *r_del;
  int i, num_entry;
  int num_del = 0;
  int del_end_pos;
  int err;

  if (sequence->feature_table == NULL) return NULL;
  ft = init_feature_table ();
  del_end_pos = position + string_len; 
  num_entry = sequence->feature_table->num_entry;
 
  for (i = 0; i < num_entry; i++) {
	e = sequence->feature_table->entry[i];	
	if (get_ft_range_check (e, position, del_end_pos) ){
	e_start = copy_ft_entry (e);
	for (r_start = e_start->range; r_start; r_start = r_start->next){
	    if (position < r_start->left->min) {
		e_del = get_ft_range_start (e, r_start);
		for (r_end = r_start; r_end; r_end = r_end->next){
		    if (del_end_pos < r_end->left->min) {
			/*printf("position<A(a) & del_end_pos<A(a)\n");*/
			err = compare_ft_single_range (r_start, r_end);
			if (err != 0) {
			    for (r_del = e_del->range; r_del; r_del = r_del->next) {
				if (r_del->next != NULL) {
				    err = compare_ft_single_range (r_del->next, r_end);
				    if (err == 0) r_del->next = NULL;   
				}
			    }
			} else e_del = NULL;
		    }
		    else if ((r_end->left->type == BASE || r_end->left->type == SIT) 
			      && del_end_pos < r_end->left->max) {
			/*printf("del_start_pos<A(a) and del_end_pos<A(b)\n");*/ 	      
			for (r_del = e_del->range; r_del; r_del = r_del->next) {
			    if (compare_ft_single_range (r_del, r_end) == 0) {
				get_ft_range_Ab_end (r_del, del_end_pos);
			    }  
			}   
		    }
		    else if (r_end->right && del_end_pos < r_end->right->min) {
			/*printf("del_start_pos<A(a) and del_end_pos<B(a)\n");*/ 
			for (r_del = e_del->range; r_del; r_del = r_del->next) {
			    err = compare_ft_single_range (r_del, r_end);
			    if (err == 0) {
				r_del->next = NULL; 
				r_del->right = NULL;
			    }  
			}						
			
		    }
		    else if ( r_end->right && (r_end->right->type == BASE || r_end->right->type == SIT) 
			      && del_end_pos < r_end->right->max) {
			/*printf("del_start_pos<A(a) and del_end_pos<B(b)\n");*/
			for (r_del = e_del->range; r_del; r_del = r_del->next) {
			    if (compare_ft_single_range (r_del, r_end)== 0) {
				get_ft_range_Bb_end (r_del, del_end_pos);	    
			    }
			}				 		    
		
		    }	   
		}  
		if (e_del != NULL) {
		    num_del ++;	
		    /* printf("num_del=%d\n", num_del);*/
		    err = extend_deleted_feature_table (ft, e_del, num_del);
		    if (err == -1) goto error; 
		}   	
		break;
	    }
	    /********/	    
	    else if( (r_start->left->type == BASE || r_start->left->type == SIT) 
		     && position < r_start->left->max) {
		e_del = get_ft_range_start (e, r_start);
		e_del->range->left->min = position + 1;
		for (r_end = r_start; r_end; r_end = r_end->next){
		    if (del_end_pos < r_end->left->min) {
			/*printf("position<A(b) & del_end_pos<A(a)(next range)\n");*/
			err = compare_ft_single_range (r_start, r_end);		  
			if (err == 0) {
			    /*printf ("error\n");*/
			    e_del = NULL;
			}		    
			for (r_del = e_del->range; r_del; r_del = r_del->next) {
			    if (r_del->next != NULL) {
				err = compare_ft_single_range (r_del->next, r_end);
				if (err == 0) r_del->next = NULL;
			    }
			}     					  		    
		    }
		    else if ( (r_end->left->type == BASE || r_end->left->type == SIT) 
			  && del_end_pos < r_end->left->max) {
			/*printf("position<A(b) & del_end_pos<A(b)\n");*/  	
			err = compare_ft_single_range (r_start, r_end);	    
			if (!err) 
			    get_ft_range_middle (e_del->range, position+1, del_end_pos);
			else {
			    for (r_del = e_del->range; r_del; r_del = r_del->next) {
				err = compare_ft_single_range (r_del, r_end);
				if (err == 0) {
				    get_ft_range_Ab_end (r_del, del_end_pos) ;		
				}  
			    }
			}		  			    	  
		    }	
		    else if (r_end->right && del_end_pos < r_end->right->min) {
			/*printf("position<A(b) & del_end_pos<B(a)\n");*/	
			for (r_del = e_del->range; r_del; r_del = r_del->next) {
			    err = compare_ft_single_range (r_del, r_end);
			    if (err == 0) {
				get_ft_range_Ba_end (r_del) ;		
			    }  
			}				    
		    }
		    else if ( r_end->right && (r_end->right->type == BASE || r_end->right->type == SIT) 
			      && del_end_pos < r_end->right->max) {
			/*printf("position<A(b) & del_end_pos<B(b)\n");*/  
			for (r_del = e_del->range; r_del; r_del = r_del->next) {
			    err = compare_ft_single_range (r_del, r_end);
			    if (err == 0) {
				get_ft_range_Bb_end (r_del, del_end_pos) ;		
			    }  
			}
		    }
		}
		/*printf("position<A(b) & del_end_pos>last B(b)\n");*/
		if (e_del != NULL) {
		    num_del ++;
		    /*printf("num_del=%d\n", num_del);*/
		    err = extend_deleted_feature_table (ft, e_del, num_del);
		    if (err == -1) goto error;
		}                  		
		    break;
	    }
/********/
	    else if (r_start->right && position < r_start->right->min) {
		e_del = get_ft_range_start (e, r_start);
		get_ft_range_Ba_start (e_del->range, position+1); 
		for (r_end = r_start; r_end; r_end = r_end->next){
		    if (del_end_pos < r_end->left->min) {
			/*printf("position<B(a) & del_end_pos<A(a)(next range)\n");*/ 
			for (r_del = e_del->range; r_del; r_del = r_del->next) {
			    if (r_del->next != NULL) {
				err = compare_ft_single_range (r_del->next, r_end);
				if (err == 0) r_del->next = NULL;
			    }
			}    		  
			mv_right_to_left (e_del->range);					
		    }			    
		    else if ((r_end->left->type == BASE || r_end->left->type == SIT) &&
			     del_end_pos < r_end->left->max) {
			/*printf("position<B(a) & del_end_pos<A(b)(next range)\n");*/
			for (r_del = e_del->range; r_del; r_del = r_del->next) {
			    if (compare_ft_single_range (r_del, r_end) == 0) {
				get_ft_range_Ab_end (r_del, del_end_pos) ;		
			    }  
			}
			mv_right_to_left (e_del->range);	
		    }
		    else if (r_end->right && del_end_pos < r_end->right->min) {
			/*printf("position<B(a) & del_end_pos<B(a))\n");*/
			err = compare_ft_single_range (r_start, r_end); 		 
			if (!err) {
			    get_ft_range_middle (e_del->range, position+1, del_end_pos);	
			} else {
			    for (r_del = e_del->range; r_del; r_del = r_del->next) {
				if (compare_ft_single_range (r_del, r_end) == 0) {
				    get_ft_range_Ba_end (r_del) ;		
				}  
			    }
			    mv_right_to_left (e_del->range); 
			}	
		    }
		    else if ( r_end->right && (r_end->right->type == BASE || r_end->right->type == SIT) 
			  && del_end_pos < r_end->right->max) {
			/*printf("position<B(a) & del_end_pos<B(b))\n");*/	
			for (r_del = e_del->range; r_del; r_del = r_del->next) {
			    err = compare_ft_single_range (r_del, r_end);
			    if (err == 0) {
				get_ft_range_Bb_end (r_del, del_end_pos) ;		
			    }  
			}
			mv_right_to_left (e_del->range);	
		    }	     					
		}
	      /* del_end_pos > last r_end->right->min(max) */
		if (e_del != NULL) {
		    num_del ++;
		    /*printf("num_del=%d\n", num_del);*/	
		    err = extend_deleted_feature_table (ft, e_del, num_del);
		    if (err == -1) goto error;
		}	 
	      break;
	    }
   /*******/
	    else if (r_start->right && (r_start->right->type == BASE 
					|| r_start->right->type == SIT) 
		     && position < r_start->right->max) {	
		e_del = get_ft_range_start (e, r_start);
		get_ft_range_Bb_start (e_del->range, position+1);	
		for (r_end = r_start; r_end; r_end = r_end->next){	
		    if (del_end_pos < r_end->left->min) {				
			for (r_del = e_del->range; r_del; r_del = r_del->next) {
			    if (r_del->next != NULL) {
				err = compare_ft_single_range (r_del->next, r_end);
				if (err == 0) r_del->next = NULL;
			    }
			} 	       		    
		    }
		    else if ((r_end->left->type == BASE || r_end->left->type == SIT) &&
			     del_end_pos < r_end->left->max) {
			/*printf("position<B(b) & del_end_pos<A(b)(next range)\n");*/
			for (r_del = e_del->range; r_del; r_del = r_del->next) {      
			    err = compare_ft_single_range (r_del, r_end);
			    if(!err) {
				get_ft_range_Ab_end (r_del, del_end_pos);
			    }      
			} 	 		    			    
		    }
		    else if (r_end->right && del_end_pos < r_end->right->min) {
			/*printf("position<B(b) & del_end_pos<B(a))\n");*/
		
			for (r_del = e_del->range; r_del; r_del = r_del->next) {      
			    err = compare_ft_single_range (r_del, r_end);
			    if(!err) {
				get_ft_range_Ba_end (r_del);
			    }      
			}	
		    }
		    else if (r_end->right && (r_end->right->type == BASE || r_end->right->type == SIT) 
			 && del_end_pos < r_end->right->max) {
			/* printf("position<B(b) & del_end_pos<B(b)\n");*/
			err = compare_ft_single_range (r_start, r_end);
			if (err == 0) {
			    get_ft_range_middle (e_del->range, position+1, del_end_pos);
			} else {
			    for (r_del = e_del->range; r_del; r_del = r_del->next) {      
				err = compare_ft_single_range (r_del, r_end);
				if(!err) {
				    get_ft_range_Bb_end (r_del, del_end_pos);
				}      
			    } 
			}		    
		    }
		    if (e_del != NULL) {
			num_del++;	
			err = extend_deleted_feature_table (ft, e_del, num_del);
			if (err == -1) goto error;
		    }
		   
		}
		break;	
	    }	   
	}
	/*break;*/	
	}
  }
  if (e_start) del_ft_entry (e_start);
  add_feature_source (ft, position, del_end_pos, sequence->name);
  return ft;

 error: 
  if (e_start) del_ft_entry (e_start);
  return NULL; 
}


/* The modified feature table to be returned. it is used for UNDO */
FEATURE_TABLE *delete_string_modify_feature_table (SEQUENCE *sequence, int string_len, int position) {

    int i, j, num_entry;
    ft_entry *e, *e_start, *e_tmp;
    ft_range *r, *r_start, *r_end;
    FEATURE_TABLE *ft;
  
    int num_del;
    int del_end_pos;
    int err, err_start;
    int deletion, lost;
    int min = 0, max = 0, min_r;

    if (!sequence->feature_table) goto error;
    ft = init_feature_table ();

    del_end_pos = position + string_len - 1;   
    num_entry = sequence->feature_table->num_entry;
    num_del = 0;
    for (i = 0; i < num_entry; i++) {
	e = sequence->feature_table->entry[i];
	e_start = copy_ft_entry (e);
	e_tmp = copy_ft_entry (e);
	if (!e_start) goto error;
	deletion = 0;
	lost = 0;

	for (r_start = e->range; r_start; r_start = r_start->next){
	    if (position < r_start->left->min) {
		for (r_end = r_start; r_end; r_end = r_end->next){
		    if (del_end_pos < r_end->left->min) {
			/*printf("position<A(a) & del_end_pos<A(a)\n");	*/	
			deletion = 1;		
			err = compare_ft_single_range (r_start, r_end);
			err_start = compare_ft_single_range (e->range, r_start);
			if (err == 0) {			    
			    for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
			    }
			}
			if (err != 0) {
			    lost = 1;
			    if (err_start == 0) {
				e->range = copy_ft_range (r_end);
				for (r = e->range; r; r = r->next) {	
				    change_single_ft_range (r, string_len);
				}
			    } 
			    if (err_start != 0) {
				for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
				}
				for (r = e->range; r; r = r->next) {
				    if (r->next != NULL) {
					err = compare_ft_single_range (r->next, r_start);
					if (!err) r->next = copy_ft_range (r_end);
				    }
				}
			    }
			}
			if (lost) {			
			    num_del ++;		
			    err = extend_deleted_feature_table (ft, e_start, num_del);
			    if (err == -1) goto error;
			} 	
			break;
		    }
		    else if ( (r_end->left->type == BASE || r_end->left->type == SIT) 
			      && del_end_pos < r_end->left->max) {
			/*printf("del_start_pos<A(a) and del_end_pos<A(b)\n");*/ 
			deletion = 1;
			err = compare_ft_single_range (r_start, r_end);
			err_start = compare_ft_single_range (e->range, r_start);
			if (err == 0) {			    
			    for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
			    }
			    r_end->left->min = position ; 
			}
			if (err != 0) {
			    lost = 1;
			    if (err_start == 0) {
				for (r = r_end; r; r = r->next) {	
				    change_single_ft_range (r, string_len);
				}
				r_end->left->min = position ; 
				e->range = copy_ft_range (r_end);
			    } 
			    if (err_start != 0) {
				for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
				}
				r_end->left->min = position;
				for (r = e->range; r; r = r->next) {
				    if (r->next != NULL) {
					err = compare_ft_single_range (r->next, r_start);
					if (!err) r->next = copy_ft_range (r_end);
				    }
				}
			    }
			}
			if (lost) {
			num_del ++;
			err = extend_deleted_feature_table (ft, e_start, num_del);
			if (err == -1) goto error;
			}
			break;
		    }			
		    else if (r_end->right && del_end_pos < r_end->right->min) {
			/*printf("del_start_pos<A(a) and del_end_pos<B(a)\n");*/
			deletion = 1;
			err = compare_ft_single_range (r_start, r_end);
			err_start = compare_ft_single_range (e->range, r_start);
			if (err == 0) {			    
			    for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
			    }
			    mv_right_to_left (r_end);    
			}
			if (err != 0) {
			    lost = 1;
			    if (err_start == 0) {
				for (r = r_end; r; r = r->next) {	
				    change_single_ft_range (r, string_len);
				}
				mv_right_to_left (r_end);
				e->range = copy_ft_range (r_end);
			    } 
			    if (err_start != 0) {
				for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
				}
				mv_right_to_left (r_end);
				for (r = e->range; r; r = r->next) {
				    if (r->next != NULL) {
					err = compare_ft_single_range (r->next, r_start);
					if (!err) r->next = copy_ft_range (r_end);
				    }
				}
			    }
			}
			if (lost) {
			    num_del ++;
			    err = extend_deleted_feature_table (ft, e_start, num_del);
			    if (err == -1) goto error;
			}
			break;
		    }
		    else if ( r_end->right && (r_end->right->type == BASE || r_end->right->type == SIT) 
			      && del_end_pos < r_end->right->max) {
			/*printf("del_start_pos<A(a) and del_end_pos<B(b)\n");*/
			deletion = 1;
			err = compare_ft_single_range (r_start, r_end);
			err_start = compare_ft_single_range (e->range, r_start);
			if (err == 0) {
			    mv_right_to_left (r_end);
			    r_end->left->min = del_end_pos + 1;
			    for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
			    }
			}
			if (err != 0) {
			    lost = 1;
			    if (err_start == 0) {
				mv_right_to_left (r_end);
				r_end->left->min = del_end_pos + 1;
				for (r = r_end; r; r = r->next) {	
				    change_single_ft_range (r, string_len);
				}
				e->range = copy_ft_range (r_end);
			    } 
			    if (err_start != 0) {
				for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
				}
				mv_right_to_left (r_end);
				 r_end->left->min = del_end_pos + 1;
				for (r = e->range; r; r = r->next) {
				    if (r->next != NULL) {
					err = compare_ft_single_range (r->next, r_start);
					if (!err) r->next = copy_ft_range (r_end);
				    }
				}
			    }
			}
			if (lost) {
			    num_del ++;
			    err = extend_deleted_feature_table (ft, e_start, num_del);
			    if (err == -1) goto error;
			}		
			break;		
		    }			   
		}
		/* del_end_pos > last r_end->right->min(max) */
		if (!deletion) {
		    /*printf("del_start_pos<A(a) and del_end_pos>last B(b)\n");*/
		    err_start = compare_ft_single_range (e->range, r_start);
		    if (err_start == 0) {
			del_ft_entry (e);
			for (j = i; j < num_entry - 1; j++) {
			    sequence->feature_table->entry[j] = sequence->feature_table->entry[j+1];
			}
			num_entry --;
			sequence->feature_table->num_entry --;
			i --;
		    }
		    if (err_start != 0) {
			for (r = e->range; r; r = r->next) {
			    if (r->next != NULL) {
				err = compare_ft_single_range (r->next, r_start);
				if (!err) r->next = NULL;
			    }
			}
		    }
		    /*printf("num_del=%d\n", num_del);*/
		    num_del ++;
		    err = extend_deleted_feature_table (ft, e_start, num_del);
		    if (err == -1) goto error;
		}
		break;  	
	    }
	    /*******/	    
	    else if( (r_start->left->type == BASE || r_start->left->type == SIT) 
		     && position < r_start->left->max) {
		for (r_end = r_start; r_end; r_end = r_end->next){
		    if (del_end_pos < r_end->left->min) {
			/*printf("position<A(b) & del_end_pos<A(a)(next range)\n");*/
			deletion = 1;
			err = compare_ft_single_range (r_start, r_end);
			err_start = compare_ft_single_range (e->range, r_start);
			/*if (err == 0) {			    
			    printf ("err\n");
			}*/
			if (err != 0) {
			    lost = 1;
			    for (r = r_end; r; r = r->next) {
				change_single_ft_range (r, string_len);
			    } 
			    for (r = e->range; r; r = r->next) {			
				err = compare_ft_single_range (r, r_start);
				if (!err) r->next = copy_ft_range (r_end); 		
			    }
			    r_start->left->max = position - 1;
			    if (r_start->right) r_start->right = NULL;
			}
			if (lost) {
			    num_del ++;
			    err = extend_deleted_feature_table (ft, e_start, num_del);
			    if (err == -1) goto error;
			}
			break;			    
		    }
		    else if ( (r_end->left->type == BASE || r_end->left->type == SIT) 
			      && del_end_pos < r_end->left->max) {
			/*printf("position<A(b) & del_end_pos<A(b)\n");*/
			deletion = 1;
			err = compare_ft_single_range (r_start, r_end);
			err_start = compare_ft_single_range (e->range, r_start);
			if (err == 0) {
			    min = r_start->left->min;
			    for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
			    }
			    r_start->left->min = min;
			}
			if (err != 0) {
			    lost = 1;
			    for (r = r_end; r; r = r->next) {
				change_single_ft_range (r, string_len);
			    } 
			     r_end->left->min = position;
			    for (r = e->range; r; r = r->next) {			
				err = compare_ft_single_range (r, r_start);
				if (!err) r->next = copy_ft_range (r_end);
			    }
			    r_start->left->max = position - 1;
			    if (r_start->right) r_start->right = NULL;
			}		
			if (lost) {
			    num_del ++;
			    err = extend_deleted_feature_table (ft, e_start, num_del);
			    if (err == -1) goto error;
			}
			break;			    
		    }			
		    else if (r_end->right && del_end_pos < r_end->right->min) {
			/*printf("position<A(b) & del_end_pos<B(a)\n");*/
		   	deletion = 1;	
			err = compare_ft_single_range (r_start, r_end);
			err_start = compare_ft_single_range (e->range, r_start);
			if (err == 0) {
			    min = r_start->left->min;
			    for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
			    }
			    r_start->left->min = min;
			    r_start->left->max = position - 1;
			}
			if (err != 0) {
			    lost = 1;
			    for (r = r_end; r; r = r->next) {
				change_single_ft_range (r, string_len);
			    } 
			    mv_right_to_left (r_end);
			    for (r = e->range; r; r = r->next) {			
				err = compare_ft_single_range (r, r_start);
				if (!err) r->next = copy_ft_range (r_end);				
			    }
			    r_start->left->max = position - 1;
			    if (r_start->right) r_start->right = NULL;
			}		
			if (lost) {
			    num_del ++;
			    err = extend_deleted_feature_table (ft, e_start, num_del);
			    if (err == -1) goto error;
			}	
			break;					    
		    }
		    else if ( r_end->right && (r_end->right->type == BASE || r_end->right->type == SIT) 
			      && del_end_pos < r_end->right->max) {
			/*printf("position<A(b) & del_end_pos<B(b)\n");*/
			deletion = 1;	 
			err = compare_ft_single_range (r_start, r_end);
			err_start = compare_ft_single_range (e->range, r_start);
			if (err == 0) {
			    min = r_start->left->min;
			    for (r = r_end; r; r = r->next) {	
				change_single_ft_range (r, string_len);
			    }
			    r_start->left->min = min;
			    r_start->left->max = position - 1;
			    r_start->right->min = position;
			}
			if (err != 0) {
			    lost = 1;
			    for (r = r_end; r; r = r->next) {
				change_single_ft_range (r, string_len);
			    } 
			    mv_right_to_left (r_end);
			    r_end->left->min = position;
			    for (r = e->range; r; r = r->next) {			
				err = compare_ft_single_range (r, r_start);
				if (!err) r->next = copy_ft_range (r_end);
			    }
			    r_start->left->max = position - 1;
			    if (r_start->right) r_start->right = NULL;
			}		
			if (lost) {
			    num_del ++;
			    err = extend_deleted_feature_table (ft, e_start, num_del);
			    if (err == -1) goto error;
			}	
			break;	
		    }
		}
	
		/* del_end_pos > last r_end->right->min(max) */
		if (!deletion) {
		    /*printf("position<A(b) & del_end_pos>last B(b)\n");*/
		    r_start->left->max = position - 1;
		     if (r_start->right) r_start->right = NULL;
		     r_start->next = NULL;
		}
		num_del ++;
		err = extend_deleted_feature_table (ft, e_start, num_del);
		if (err == -1) goto error;
		break;
	    }
/********/
	    else if (r_start->right && position < r_start->right->min) {
		    for (r_end = r_start; r_end; r_end = r_end->next){
			if (del_end_pos < r_end->left->min) {
			    /*printf("position<B(a) & del_end_pos<A(a)(next range)\n");*/
			    deletion = 1;
			    err = compare_ft_single_range (r_start, r_end);
			    err_start = compare_ft_single_range (e->range, r_start);
			    /*if (err == 0) {			    
				printf ("err\n");
				}*/
			    if (err != 0) {
				lost = 1;
				for (r = r_end; r; r = r->next) {
				    change_single_ft_range (r, string_len);
				} 
				for (r = e->range; r; r = r->next) {			
				    err = compare_ft_single_range (r, r_start);
				    if (!err) r->next = copy_ft_range (r_end);
				}
				if (r_start->right) r_start->right = NULL;
			    }
			    if (lost) {
				num_del ++;
				err = extend_deleted_feature_table (ft, e_start, num_del);
				if (err == -1) goto error;
			    }
			    break;			    
			}
			else if ((r_end->left->type == BASE || r_end->left->type == SIT) &&
				 del_end_pos < r_end->left->max) {
			    /*printf("position<B(a) & del_end_pos<A(b)(next range)\n");*/
			    deletion = 1;
			    err = compare_ft_single_range (r_start, r_end);
			    err_start = compare_ft_single_range (e->range, r_start);
			    /* if (err == 0) {			    
				printf ("err\n");
			    }*/
			    if (err != 0) {
				lost = 1;
				for (r = r_end; r; r = r->next) {
				    change_single_ft_range (r, string_len);
				}
				r_end->left->min = position;
				for (r = e->range; r; r = r->next) {			
				    err = compare_ft_single_range (r, r_start);
				    if (!err) r->next = copy_ft_range (r_end);
				}
				if (r_start->right) r_start->right = NULL;
			    }
			    if (lost) {
				num_del ++;
				err = extend_deleted_feature_table (ft, e_start, num_del);
				if (err == -1) goto error;
			    }
			    break;			    
			}
   /*tttttt*/
			else if (r_end->right && del_end_pos < r_end->right->min) {
			    /* printf("position<B(a) & del_end_pos<B(a))\n");*/
			    deletion = 1;
			    err = compare_ft_single_range (r_start, r_end);
			    err_start = compare_ft_single_range (e->range, r_start);
			    if (err == 0) {			    
				min = r_start->left->min;
				if (r_start->left->type != EXACT) max = r_start->left->max;
				for (r = r_end; r; r = r->next) {	
				    change_single_ft_range (r, string_len);
				}
				r_start->left->min = min;	
				if (r_start->left->type != EXACT) r_start->left->max = max;
			    }
			    if (err != 0) {
				lost = 1;
				for (r = r_end; r; r = r->next) {
				    change_single_ft_range (r, string_len);
				}
				mv_right_to_left (r_end);
				for (r = e->range; r; r = r->next) {			
				    err = compare_ft_single_range (r, r_start);
				    if (!err) r->next = copy_ft_range (r_end);
				}
				if (r_start->right) r_start->right = NULL;
			    }
			    /* printf("num_del=%d\n", num_del);*/
			    
			    if (lost) {
				num_del ++;
				err = extend_deleted_feature_table (ft, e_start, num_del);
				if (err == -1) goto error;
			    }
			    break;	
			}
			else if ( r_end->right && (r_end->right->type == BASE || r_end->right->type == SIT) 
				  && del_end_pos < r_end->right->max) {
			    /*printf("position<B(a) & del_end_pos<B(b))\n");*/
			    deletion = 1;
			    err = compare_ft_single_range (r_start, r_end);
			    err_start = compare_ft_single_range (e->range, r_start);
			    if (err == 0) {			    
				min = r_start->left->min;
				if (r_start->left->type != EXACT) max = r_start->left->max;
				   for (r = r_end; r; r = r->next) {	
				   change_single_ft_range (r, string_len);
				}
				r_start->left->min = min;	
				if (r_start->left->type != EXACT) r_start->left->max = max;	
			    }
			    if (err != 0) {
				lost = 1;
				for (r = r_end; r; r = r->next) {
				    change_single_ft_range (r, string_len);
				}
				mv_right_to_left (r_end);
				r_end->left->min = position;
				for (r = e->range; r; r = r->next) {			
				    err = compare_ft_single_range (r, r_start);
				    if (!err) r->next = copy_ft_range (r_end);	
				}
				if (r_start->right) r_start->right = NULL;
			    }
			    if (lost) {
				num_del ++;
				err = extend_deleted_feature_table (ft, e_start, num_del);
				if (err == -1) goto error;
			    }
			    break;	
			}		     					
		    }
		    /* del_end_pos > last r_end->right->min(max) */
		    if (!deletion) {
			/*printf("position<B(a) & del_end_pos> last B(b))\n");*/
			if (r_start->right) r_start->right->min = position;
			r_start->next = NULL;
			num_del ++;
			err = extend_deleted_feature_table (ft, e_start, num_del);
			if (err == -1) goto error;
		    }
		    break;
	    }
   /*******/
	    else if (r_start->right && (r_start->right->type == BASE || r_start->right->type == SIT) 
		     && position < r_start->right->max) {
		for (r_end = r_start; r_end; r_end = r_end->next){
		    if (del_end_pos < r_end->left->min) {
			/*printf("position<B(b) & del_end_pos<A(a)(next range)\n");*/
			deletion = 1;				
			err = compare_ft_single_range (r_start, r_end);
			    err_start = compare_ft_single_range (e->range, r_start);
			    /*if (err == 0) {			    
				printf ("err\n");
			    }*/
			    if (err != 0) {
				lost = 1;
				for (r = r_end; r; r = r->next) {
				    change_single_ft_range (r, string_len);
				} 
				for (r = e->range; r; r = r->next) {			
				    err = compare_ft_single_range (r, r_start);
				    if (!err) {
					r->next = copy_ft_range (r_end);
					break;
				    }
				}
				r_start->right->max = position - 1;
				if (r_start->right->min > r_start->right->max) r_start->right = NULL;	
			    }
			    if (lost) {
				num_del ++;
				err = extend_deleted_feature_table (ft, e_start, num_del);
				if (err == -1) goto error;
			    }
			break;			    
		    }
		    else if ((r_end->left->type == BASE || r_end->left->type == SIT) &&
			     del_end_pos < r_end->left->max) {
			/*printf("position<B(b) & del_end_pos<A(b)(next range)\n");*/
			deletion = 1;
			err = compare_ft_single_range (r_start, r_end);
			    err_start = compare_ft_single_range (e->range, r_start);
			    /*if (err == 0) {			    
				printf ("err\n");
			    }*/
			    if (err != 0) {
				lost = 1;
				for (r = r_end; r; r = r->next) {
				    change_single_ft_range (r, string_len);
				}
				r_end->left->min = position;
				for (r = e->range; r; r = r->next) {			
				    err = compare_ft_single_range (r, r_start);
				    if (!err) r->next = copy_ft_range (r_end);
				}
				r_start->right->max = position - 1;
				if (r_start->right->min > r_start->right->max) r_start->right = NULL;
			    }
			    if (lost) {
				num_del ++;
				err = extend_deleted_feature_table (ft, e_start, num_del);
				if (err == -1) goto error;
			    }
			break;			    
		    }
		 
		    else if (r_end->right && del_end_pos < r_end->right->min) {
			/*printf("position<B(b) & del_end_pos<B(a))\n");*/
			deletion = 1;
			err = compare_ft_single_range (r_start, r_end);
			    err_start = compare_ft_single_range (e->range, r_start);
			    /* if (err == 0) 
				printf("err\n");*/
			    if (err != 0) {
				lost = 1;
				for (r = r_end; r; r = r->next) {
				    change_single_ft_range (r, string_len);
				}
				mv_right_to_left (r_end);
				for (r = e->range; r; r = r->next) {			
				    err = compare_ft_single_range (r, r_start);
				    if (!err) r->next = copy_ft_range (r_end);
				}
				r_start->right->max = position - 1;
				if (r_start->right->min > r_start->right->max) r_start->right = NULL;
			    }
			    if (lost) {
				num_del ++;
				err = extend_deleted_feature_table (ft, e_start, num_del);
				if (err == -1) goto error;
			    }	
			break;	
		    }
		    else if (r_end->right && (r_end->right->type == BASE || r_end->right->type == SIT) 
			     && del_end_pos < r_end->right->max) {
			/* printf("position<B(b) & del_end_pos<B(b)\n");*/
			deletion = 1;
			err = compare_ft_single_range (r_start, r_end);
			    err_start = compare_ft_single_range (e->range, r_start);
			    if (err == 0) {			    
				min = r_start->left->min;
				if (r_start->left->type != EXACT) max = r_start->left->max;
				min_r = r_start->right->min;
				for (r = r_end; r; r = r->next) {	
				    change_single_ft_range (r, string_len);
				}
				r_start->left->min = min;	
				if (r_start->left->type != EXACT) r_start->left->max = max;
				r_start->right->min = min_r;
			    }
			    if (err != 0) {
				lost = 1;
				for (r = r_end; r; r = r->next) {
				    change_single_ft_range (r, string_len);
				}
				mv_right_to_left (r_end);
				r_end->left->min = position;
				for (r = e->range; r; r = r->next) {			
				    err = compare_ft_single_range (r, r_start);
				    if (!err) r->next = copy_ft_range (r_end);	
				}
				r_start->right->max = position - 1;
				if (r_start->right->min > r_start->right->max) r_start->right = NULL;
			    }
			    if (lost) {
				num_del ++;
				err = extend_deleted_feature_table (ft, e_start, num_del);
				if (err == -1) goto error;
			    }		
			break;	
		    }			    		    
		}
		/* del_end_pos > last r_end->right->min(max) */
		    if (!deletion) {
			/*printf("position<B(b) & del_end_pos> last B(b)\n");*/
			r_start->right->max = position - 1;
			if (r_start->right->min > r_start->right->max) r_start->right = NULL;	
			r_start->next = NULL;
			num_del++;
			err = extend_deleted_feature_table (ft, e_start, num_del);
			if (err == -1) goto error;	
		    }		
		break;	
	    }
/*******/  		   	        	    		   	       
	}
}
    del_ft_entry (e_start);
    return ft;
 error:
    if(ft) free_feature_table (ft);
    if(e_start) del_ft_entry (e_start);
    return NULL;
}




