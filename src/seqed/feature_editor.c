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
#include "sequence_formats.h"
#include "parse_feature.h"
#include "parse_db.h"
#include "feature_editor.h"

#define TAGDB "FEATQUALDB"

feat_qual_db *fqual_db = NULL;
int fqual_db_count = 0;

SEQUENCES static *fsequence = NULL; /* feature_editor */

static void fqualdb_parse(char *fn) {
    pf_spec a[] = {
	{"id",	     PF_STR, offsetof(feat_qual_db, qual)},
	{"value",     PF_STR, offsetof(feat_qual_db, value)},
	{"entry",    PF_STR, offsetof(feat_qual_db, entry)},
	{NULL,	     0,	      0}
    };

    fqual_db = (feat_qual_db *)parse_file(fn, a, fqual_db, &fqual_db_count,
					 sizeof(*fqual_db), NULL);
}

/*static void tidyUpFqualDBFields(int fqual) {

    int len;
    
    if (fqual_db[fqual].qual == NULL)
	fqual_db[fqual].qual = fqual_db[fqual].type;

    len =  strlen(fqual_db[fqual].qual);
    if (len < 4)
	strncpy(fqual_db[fqual].qual,"    ",4);
    else
	len = 4;
    strncpy(fqual_db[fqual].id,fqual_db[fqual].search_id,len);
    
    }*/



void readInFqualDB(void)
{
    char *path, *p;
    char tmp_path[2000];
    int gotfile = 0;

    /* 22/1/99 johnt - fallback to to use STADTABL is GFCOLDB isn't defined */
    if (NULL == (path = (char *)getenv(TAGDB))){
        if(getenv("STADTABL")) {
	   strcpy(tmp_path,getenv("STADTABL"));
           strcat(tmp_path,"/");
	   strcat(tmp_path,TAGDB);
           path = tmp_path;
	} else {
	    path = TAGDB;
	}
    }

    do {
	p = strrchr(path, PATHSEP); /* 22/1/99 johnt - path seperator now macro to allow WINNT support */
	if (p) {
	    *p = 0;
	    p++;
	} else
	    p = path;
	    
	if (file_exists(p)) {
	    fqualdb_parse(p);
            gotfile++;
	}
    } while (p != path);
    /*for (i = 0; i < fcol_db_count; i++) {
	tidyUpFqualDBFields(i);
	}*/
    if( !gotfile )
	verror(ERR_WARN, "Tag DB", "No Files found - please check GTAGDB environment variable.");

}

void get_fqual_types (void) {

    static int done = 0;

    if (!done) {
	readInFqualDB();
	done = 1;
    }
}

int tcl_get_qual_array(ClientData clientData, 
		       Tcl_Interp *interp,
		       int argc, 
		       char **argv) {
    Tcl_DString fqual;
    int i;

    if (!fqual_db) {
	get_fqual_types();
    }
    Tcl_DStringInit(&fqual);
    for (i = 0; i < fqual_db_count; i++) {	
	Tcl_DStringStartSublist(&fqual);
	Tcl_DStringAppendElement(&fqual, fqual_db[i].type);
	Tcl_DStringEndSublist(&fqual);
    }
    Tcl_DStringResult(interp, &fqual);

    return TCL_OK;
}

int GetFseqNum(int seq_id) {

    int i, num_seq;

    num_seq = fsequence->num_seq;
   
    for (i = 0; i < num_seq; i++) {
	if (fsequence->sequence[i]->id == seq_id) {
	    return i;
	}
    }
    return -1;
}


int GetFtList (ClientData clientData,
	       Tcl_Interp *interp, 
	       int argc, 
	       char **argv)
{
    FEATURE_TABLE *ft;
    int n;

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], "seq_num \"", (char*)NULL);
        return TCL_ERROR;
    } 
    n = atoi (argv[1]);
    ft = fsequence->sequence[n]->feature_table;
    if (ft != NULL) {

	ft_entry *e;
	int i, num_entry;

	num_entry = ft->num_entry;

	for (i = 0; i < num_entry; i++) {
	    
	    Tcl_DString dstr, loc, qua;

	    e = ft->entry[i];
	    Tcl_DStringInit(&dstr);
	    Tcl_DStringInit(&loc);
	    Tcl_DStringInit(&qua);

	    vTcl_DStringAppendElement(&dstr, "%s",
				      e->type ? e->type : "");
	    vTcl_DStringAppendElement(&dstr, "%s",
				      e->location ? e->location : "");
	    vTcl_DStringAppendElement(&dstr, "%s",
				      e->qualifiers ? e->qualifiers : "");

	    Tcl_AppendElement(interp, Tcl_DStringValue(&dstr));
	    Tcl_DStringFree(&dstr);
	}	
    }
    return TCL_OK;  
}

int MadeCopy (ClientData clientData,
	      Tcl_Interp *interp, 
	      int argc, 
	      char **argv)
{
   int i, seq_id, num_seq;
 
    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], "seq_id \"", (char*)NULL);
        return TCL_ERROR;
    }
    
    if (!fsequence) {
	if (NULL == (fsequence = init_sequences( ))) 
	    goto err;
    }
    num_seq = fsequence->num_seq;
    seq_id = GetSequenceIdByName (argv[1]);

    /* check if already in editor_sequence array */
    for (i = 0; i < num_seq; i++) {
	if (fsequence->sequence[i]->id == seq_id) {
	    vTcl_SetResult(interp, "%d", i);
	    return TCL_OK;
	}
    }
    /* if not, make a copy in editor_sequence array */
    if (num_seq != 0) {
	    fsequence = realloc_sequences (fsequence, num_seq);
    }
    fsequence->sequence[num_seq] = copy_sequence (GetSequencesSequence (seq_id));
    fsequence->num_seq ++;
    
    vTcl_SetResult(interp, "%d", num_seq);
    return TCL_OK; 
 err:
    return TCL_OK;
}
/* return 0: pass.
   return 1: complement range.
   return 2: there are more than one range.
   return -1:fail, not correct format.
   return -2:fail, out of range.
*/
int check_location (char *s, int length) {

    char *l;
    int complement = 0;
    int join = 0;
    int join_coma = 0;
    int base;

    enum state_t {
	STATE_TEXT,
        STATE_COMP,
	STATE_JOIN,    
	STATE_LEFT_START,
        STATE_LEFT_MIN,
	STATE_LEFT_SINGLE_DOT, 
	STATE_LEFT_CARET,
	STATE_LEFT_MAX,
	STATE_SINGLE_DOT,
	STATE_DOUBLE_DOT,
	STATE_RIGHT_START,
	STATE_RIGHT_MIN,
	STATE_RIGHT_SINGLE_DOT,
        STATE_RIGHT_CARET,
	STATE_RIGHT_MAX,
	STATE_NEXT,
	STATE_EXIT    	    /* 6 Newline or nul char */
    } state;
   
    purify_range (s);

    state = STATE_TEXT;
    while (state != STATE_EXIT) {
	switch (state) {
	case STATE_TEXT:
	    
	    l = NULL;
	    if (*s == '<' || *s == '>' || *s == '(' || *s == 'c' || *s == 'j' || isdigit(*s)) {
		state = STATE_LEFT_START;
	    } else {
		printf ("STATE_TEXT:\n");	
		return -1;
	    }
	case STATE_LEFT_START:
	    
	    if (*s == '<' || *s == '>' || *s == '(') {
	       
		state = STATE_LEFT_MIN;
	    } else if (*s == 'c') {
		state = STATE_COMP;
	    } else if (*s == 'j') {
		state = STATE_JOIN;
	    } else if (isdigit(*s)) {
		l = s;
		state = STATE_LEFT_MIN;	
	    } else {
		printf ("STATE_LEFT_START:\n");
		return -1;
	    }
	    break;
	  case STATE_JOIN:
	      join = 1;
	      if (isdigit(*s)) {
		 l = s;
		 state = STATE_LEFT_MIN; 
	      }
	      break;
	   case STATE_COMP:
	       complement = 1;
	       if (isdigit(*s)) {
		   l = s;
		   state = STATE_LEFT_MIN;	
	       }
	       break;
	  
	case STATE_LEFT_MIN:
	    
	    if (isdigit(*s)) {
		state = STATE_LEFT_MIN;
	    } else if (*s == '.') {
		base = atoi(l);
		if (base < 1 || base > length) return -2; 
		state = STATE_LEFT_SINGLE_DOT;
	    } else if (*s == '^') {
		base = atoi(l);
		if (base < 1 || base > length) return -2; 
		state = STATE_LEFT_CARET;
	    } else if (*s == 'j') { /*fixme, may be no needed*/
		l = NULL;
		state = STATE_JOIN;
	    
	    } else {
		printf ("STATE_LEFT_MIN:\n");
		return -1;
	    }	    
	    break;
	case STATE_LEFT_SINGLE_DOT:
	    
	    if (*s == '.') {
		state = STATE_RIGHT_START;
	    } else if (isdigit(*s)) {
		l = s;
		state = STATE_LEFT_MAX;
	    } else {
		printf ("STATE_LEFT_SINGLE_DOT:\n");	
		return -1;
	    }
	    break;
	case STATE_LEFT_CARET:
	    
	    if (isdigit(*s)) {
		l = s;
		state = STATE_LEFT_MAX;	
	    } else {
		printf (":STATE_LEFT_CARET:\n");
		return -1;
	    }
	    break;
	case STATE_LEFT_MAX:
	    
	    if (isdigit(*s)) {
		state = STATE_LEFT_MAX;
	    } else if (*s == ')' || *s == '\0') {
		base = atoi(l);
		if (base < 1 || base > length) return -2; 
		state = STATE_SINGLE_DOT; 
	    } else {
		printf ("STATE_LEFT_MAX:\n");
		return -1;
	    }	    
	    break;

	case STATE_SINGLE_DOT:
	    
	    if (*s == '.') {
		state = STATE_DOUBLE_DOT;
	    } else {
		printf ("STATE_SINGLE_DOT:\n");
		return -1;
	    }
	    break;

	case STATE_DOUBLE_DOT:
	    
	    if (*s == '.') {
		state = STATE_RIGHT_START;
	    } else {
		printf ("STATE_DOUBLE_DOT\n");
		return -1;
	    }
	    break;

	case STATE_RIGHT_START:
	    
	    if (*s == '<' || *s == '>' || *s == '(') {
		state = STATE_RIGHT_MIN;
	    } else if (isdigit(*s)) {
		l = s;
		state = STATE_RIGHT_MIN;
	    } else {	
		printf ("STATE_RIGHT_START:\n");
		return -1;
	    }
	    break;

	case STATE_RIGHT_MIN:
	    
	    if (isdigit(*s)) {
		state = STATE_RIGHT_MIN;
	    } else if (*s == '.') {
		base = atoi(l);
		if (base < 1 || base > length) return -2; 
		state = STATE_RIGHT_SINGLE_DOT; 
	    } else if (*s == '^') {
		base = atoi(l);
		if (base < 1 || base > length) return -2; 
		state = STATE_RIGHT_CARET;
	    } else if (*s == ')') {
		base = atoi(l);
		if (base < 1 || base > length) return -2; 
		state = STATE_NEXT;
	    } else if (*s == ',') {
		base = atoi(l);
		if (base < 1 || base > length) return -2; 
		join_coma = 1;
		state = STATE_LEFT_START;
	    } else {
		printf ("STATE_RIGHT_MIN:\n");
		return -1;
	    }	    
	    break;

	case STATE_RIGHT_SINGLE_DOT:
	   
	    if (isdigit(*s)) {
		l = s;
		state = STATE_RIGHT_MAX;
	   
	    } else {
		printf ("STATE_RIGHT_SINGLE_DOT:\n");
		return -1;
	    }
	    break;

	case STATE_RIGHT_CARET:
	   
	    if (isdigit(*s)) {
		l = s;
		state = STATE_RIGHT_MAX;	
	    } else {
		printf (":STATE_RIGHT_CARET:\n");
		return -1;
	    }
	    break;

	case STATE_RIGHT_MAX:
	    
	    if (isdigit(*s)) {
		state = STATE_RIGHT_MAX;
	    } else if (*s == ')') {
		base = atoi(l);
		if (base < 1 || base > length) return -2; 
		state = STATE_NEXT; 
	    } else if (*s == ',') {
		base = atoi(l);
		if (base < 1 || base > length) return -2; 
		 join_coma = 1;
		state = STATE_LEFT_START;
	    } else {
		printf ("STATE_RIGHT_MAX:\n");
		return -1;
	    }	    
	    break;
	case STATE_NEXT:
	    
	    if (*s == ',') {
		 join_coma = 1;
		state = STATE_LEFT_START;
	    } else if (*s == ')') {
		state = STATE_EXIT;
	    } else {
		printf ("STATE_NEXT:\n");	
		return -1;
	    }	   
	    
	    break;
	case STATE_EXIT:
	    printf ("STATE_EXIT\n");
	    
	}
	s++;
	if (!*s || *s == '/') {
	    base = atoi(l);
	    if (base < 1 || base > length) return -2; 
	    state = STATE_EXIT;
	}
    }
    if (complement) return 1;
    if (join && !join_coma) return -1;
    if (join && join_coma) return 0;
    if (!join && join_coma) return 2;
    return 0;
     
}

int get_value_format (char *qualifier) {

    int i;
    
    for (i = 0; i < fqual_db_count; i++) {	
	if (!strcmp (fqual_db[i].type, qualifier)) {
	    return atoi(fqual_db[i].value);  
	}
    }
    return -1;
}

int CheckValueFormat(ClientData clientData,
		     Tcl_Interp *interp, 
		     int argc, 
		     char **argv) 
{    
    char *q;
    int v;

    if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], "qualifier\"", (char*)NULL);
	    return TCL_ERROR;
    }
    q = argv[1];
    v = get_value_format (q);
    vTcl_SetResult(interp, "%d", v);
    return TCL_OK; 
}


int check_validation (char *s) {

    int i, len;

    len = strlen (s);
    /* delete extra space */
    while (isspace (s[len - 1])) {	
	s [len - 1] = 0;
	len --;
    }
    /* check validation */
    for (i = 0; i < fqual_db_count; i++) {	
	if (!strcmp (fqual_db[i].type, s)) {
	    return 0;
	}
    }
    return -1;
}

/*
return 0: successful
return 1: invalid qualifier
return 2: expect a equal after the qualifier
return 3: expect a quote after the equal
return 4: expect a quote at the end of the comment
return 5: unexpect a quoat after the equal
return 6: unexpect a quoat at the end of the comment
 */
int check_qualifier (char *str, char **name) {

    char *cp, *qualifier, *cpt;
    int not_done = 1, cv;
    enum state_t {
	STATE_START,		/* start */
	STATE_QUALIFIER,	/* start a qualifier */
	STATE_EQUALS,		/* found an = */
	STATE_QUOTE,		/* text format*/
	STATE_NEXT_QUOTE,
	STATE_NQUOTE,           /* string format */
	STATE_NEXT_NQUOTE
    } state;

    if (!str)
	return -1;
    cp = strdup (str);
    state = STATE_START;
   
    while (not_done) {
	switch(state) {
	case STATE_START:
	    qualifier = NULL;
	    if (*cp == '/') { /* start a new qualifier */
		state = STATE_QUALIFIER; 
		cp++;
		qualifier = cp; 
	    } else if (*cp) { /* next line for same qualifier */
		cp++;
	     
	    } else {
		not_done = 0; /* end for this entry */
	    }
	    break;
	    
	case STATE_QUALIFIER:
	   
	    if (*cp == '=') {/* start adding comment */
		*cp = 0;
		cv = check_validation (qualifier);
		if (cv == 0 ) {	
		    state = STATE_EQUALS; 
		    cp++;		  
		} else {/* invalid qualifier */
		    *name = qualifier;
		    return 1;
		}
	    } else if (*cp == '\n' || *cp == '\0') { /* value format: none */
		*cp = 0;
		cv = check_validation (qualifier);
		if (cv == 0) {
		    if (get_value_format (qualifier) == 2) {
			state = STATE_START; /* start next qualifier */
			cp++;
		     
		    } else { /* expect a equal */
			*name = qualifier;
			return 2; 
		    }
		} else { /* invalid qualifier */ 
		   
		    *name = qualifier;
		    return 1;
		}
	    } else if (*cp == '"' ) {
		*cp = 0;
		cv = check_validation (qualifier);
		if (cv == 0) {/* expect a equal */	  
		  *name = qualifier;
		  return 2;   
		}
	    }
	    else {
	        cp++; /* continue get qualifier ID */	       
	    }
	    break;

	case STATE_EQUALS:
	    if (get_value_format (qualifier) == 0) {
		state = STATE_QUOTE; 
	    } else if (get_value_format (qualifier) == 1) {
		state = STATE_NQUOTE;
	    }
	    break;

	case STATE_QUOTE:
	    if (*cp == '"') {
		state = STATE_NEXT_QUOTE;
		cp++;
	     
	    } else if (*cp) {/*	expect a quote-start */	
		*name = qualifier;
		return 3;
	    }
	    break;

	case STATE_NEXT_QUOTE: 
	   if (*cp == '\0' || *cp == '\n') {	       
	       cpt = cp;
	       cpt++;
	       if (!*cpt) {
		   while (*cp == '\0' || *cp == '\n' || *cp == ' ') {
		       cp--;
		   }
		   if (*cp =='"') {
		     not_done = 0;
		   } else {
		       *name = qualifier;
		       return 4;
		   }  
	       } else if (*cpt == '/') {
		   while (*cp == '\0' || *cp == '\n' || *cp == ' ') {
		       cp--;
		   }
		   if (*cp =='"') {
		     
		     state = STATE_START; 
		     cp++; 

		   } else { /*expect a quote-end */
		       *name = qualifier;
		       return 4;
		   }  
	       } else {
		   state = STATE_NEXT_QUOTE;
		   cp++;
	       }
	   } else {
	       cp++;
	   }
	   break;
	
	case STATE_NQUOTE:
	    if (*cp == '"') { /* unexpect a quote-start */
		*name = qualifier;
		return 5;
	    } else if (*cp == '\0' || *cp == '\n') {
		state = STATE_START; /* start a new qualifier */
		cp++;
	    } else if (*cp) {
		state = STATE_NEXT_NQUOTE;
		cp++;
	    }
	    break;
	case STATE_NEXT_NQUOTE: 
	    if (*cp == '\0' || *cp == '\n') {	       
	       cpt = cp;
	       cpt++;
	       if (!*cpt) {
		   while (*cp == '\0' || *cp == '\n' || *cp == ' ') {
		       cp--;
		   }
		   if (*cp =='"') { /*unexpect a quote-end */  
		       *name = qualifier;
		       return 6;    
		   } else {
		       not_done = 0;  
		   }  
	       } else if (*cpt == '/') {
		   while (*cp == '\0' || *cp == '\n' || *cp == ' ') {
		       cp--;
		   }
		   if (*cp =='"') { /*unexpect a quote-end */ 
		       *name = qualifier;
		       return 6;
		   } else { 
		      state = STATE_START; 
		     cp++;  
		   }  
	       } else {
		   state = STATE_NEXT_NQUOTE;
		   cp++;
	       }
	   } else {
	       cp++;
	   }
	    break;
	}
    }
    *name = "all";
    return 0;
}

int delete_feature_from_sequence (char *seq_name, int entry_id) {

    FEATURE_TABLE *ft;
    int i, seq_id, n, num_entry;

    seq_id = GetSequenceIdByName (seq_name);
    n = GetFseqNum (seq_id);
    ft = fsequence->sequence[n]->feature_table;
    
    if (ft == NULL) return -1;
    num_entry = ft->num_entry;

    for (i = 0; i < num_entry; i++) {
	if ( i == entry_id) {   
	    memmove(&ft->entry[i], &ft->entry[i+1], (num_entry - i - 1) * sizeof(ft->entry[i]));
	    ft->num_entry--;
	    if (entry_id == ft->num_entry) return i-1;
	    return i;
	}
    }
    return -1;
}

int add_feature_to_sequence (char *seq_name, char *key, char *location, char *qualifier) {

    FEATURE_TABLE *ft;
    ft_entry *e;
    int num_entry;
    int seq_id, n;
  
    /* create new entry */
    if (NULL == (e = new_ft_entry())) goto error;
    sprintf (e->type, "%s", key);
    e->location = strdup (location);
    e->qualifiers = strdup (qualifier);
    parse_ft_location (e, location);
    init_ft_qual_hash (e, qualifier);
    
    /* add to sequence */
    seq_id = GetSequenceIdByName (seq_name);
    n = GetFseqNum (seq_id);
    ft = fsequence->sequence[n]->feature_table;
    /* fixme: use new function: add_entry_to_feature_table */ 
    num_entry = ft->num_entry;   
    num_entry ++;
    if (NULL == realloc_feature_table (ft, num_entry)) goto error;
    ft->entry[num_entry - 1] = e;
    ft->num_entry = num_entry;

    return num_entry - 1;
 error:
    if (e) del_ft_entry (e);
    return -1;
}

int change_feature_in_sequence (char *seq_name, int entry_id, char *key, char *loc, char *qua) {

    FEATURE_TABLE *ft;
    ft_entry *e;
    int seq_id, n;
  
    /* create new entry */
    if (NULL == (e = new_ft_entry())) goto error;
    sprintf (e->type, "%s", key);
    e->location = strdup (loc);
    e->qualifiers = strdup (qua);
    parse_ft_location (e, loc);
    init_ft_qual_hash (e, qua);

    /* replace */
    seq_id = GetSequenceIdByName (seq_name);
    n = GetFseqNum (seq_id);

    ft = fsequence->sequence[n]->feature_table;
    del_ft_entry (ft->entry[entry_id]);
    ft->entry[entry_id] = e;
    
     return entry_id;
 error:
    if (e) del_ft_entry (e);
    return -1;
    
}

void remove_sequence (int n) {

    int i, num_seq;
    num_seq = fsequence->num_seq;

    for (i = 0; i < num_seq; i++) {
	if (i == n) {
	    memmove(&fsequence->sequence[n], &fsequence->sequence[n+1], 
		    (num_seq - n - 1) * sizeof(fsequence->sequence[n]));
	    fsequence->num_seq--;
	    i--;
	    num_seq--;
	}    
    }
	 
}

int FeatureEditor (ClientData clientData,
		   Tcl_Interp *interp, 
		   int argc, 
		   char **argv) 
{    
    char *op;
    int err;
    
    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], "operation\"", (char*)NULL);
        return TCL_ERROR;
    } 
    op = argv[1]; 
    if (!strcmp (op, "apply")) {

	ft_entry *entry = NULL;
	char *seq_name, *key, *loc, *qua;
	int mode, entry_id;
	if (argc != 8) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			 argv[0], " operation\"", "apply_mode\"", "entry_id\"",   
			   "keyword\"", "location\"", "comments\"", (char*)NULL);
	    return TCL_ERROR;
	} 
	if (NULL == (entry = new_ft_entry()))
	    goto error;
	mode = atoi (argv[2]);
	seq_name = argv[3];
	entry_id = atoi (argv[4]);
	key = argv[5];
	loc = argv[6];
	qua = argv[7];
 
	purify_range (loc);
	if (mode == 0) {
	    printf ("change\n");
	    err = change_feature_in_sequence (seq_name, entry_id, key, loc, qua);
	} 
	if (mode == 1) {
	    printf ("create\n");
	    err = add_feature_to_sequence (seq_name, key, loc, qua);
	}
	vTcl_SetResult(interp, "%d", err);
    }
    
    if (!strcmp (op, "delete")) {

	int entry_id;
	char *seq_name;

	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " operation\"", "entry_id\"", (char*)NULL);
	    return TCL_ERROR;
	}
	seq_name = argv[2];
	entry_id = atoi (argv[3]);
      
	err = delete_feature_from_sequence (seq_name, entry_id); 
	vTcl_SetResult(interp, "%d", err);

    }
    
    if (!strcmp (op, "save")) {
	
	int seq_id, seq_num;
	SEQUENCE *s;

	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " operation\"", "seq_id\"", (char*)NULL);
	    return TCL_ERROR;
	}
	seq_id = GetSequenceIdByName (argv[2]);
	seq_num = GetFseqNum (seq_id);
	s = fsequence->sequence[seq_num];
	/*err = save_change_to_sequence (seq_id, s);*/
	err = sequence_save (s);
	if (!err) remove_sequence (seq_num);
	vTcl_SetResult(interp, "%d", err);	
    }

    if (!strcmp (op, "location")) {
	
	char *t;
	int seq_id, seq_len;

	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], "operation\"", "location\"", (char*)NULL);
	    return TCL_ERROR;
	} 

	if (NULL == (t = (char *)xmalloc((strlen (argv[3]) + 1)*sizeof(char)))) 
	    goto error;
	strcpy (t, argv[3]);
	seq_id = GetSequenceIdByName (argv[2]);
	seq_len = GetSequenceLength (seq_id);
	err = check_location (t, seq_len);  
	if (err == 1) { /* already complement */
	    t = strchr (t, '(');
	    memmove (&t[0], &t[1], (strlen(t)-1));
	    t[strlen(t) - 2] =0;
	}
	vTcl_SetResult(interp, "%d %s", err, t);
    }

    if (!strcmp (op, "qualifier")) {
	char *qua, *name;
	if (argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], "operation\"", "comments\"",(char*)NULL);
	    return TCL_ERROR;
	} 
	if (NULL == (qua = (char *)xmalloc((strlen (argv[2]) + 1)*sizeof(char)))) 
	    goto error;
	strcpy (qua, argv[2]);
	err = check_qualifier (qua, &name);
	
	vTcl_SetResult(interp, "%d %s", err, name);
    }

   return TCL_OK;
 error:
  
   return TCL_OK;
}

int GetFseqNumTcl (ClientData clientData,
		Tcl_Interp *interp, 
		int argc, 
		char **argv) 
{   
    int seq_id, n;
 
    if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], "seq_name\"",(char*)NULL);
	    return TCL_ERROR;
	} 

    seq_id = GetSequenceIdByName (argv[1]);
    n = GetFseqNum (seq_id);

    vTcl_SetResult(interp, "%d", n);
    return TCL_OK;
}


int FeatureEditor_Init(Tcl_Interp *interp) {

    if (Itcl_RegisterC(interp, "ft_editor", FeatureEditor, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "check_format", CheckValueFormat, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "made_copy", MadeCopy, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "get_ft",  GetFtList, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "get_num",  GetFseqNumTcl, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    
    return TCL_OK;
}
