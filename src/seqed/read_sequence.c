#include <tcl.h>
#include <tk.h>
#include <itcl.h>
#include <itk.h>

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "misc.h"
#include "tcl_utils.h"
#include "sequence_formats.h"
#include "read_sequence.h"
#include "editor_reg.h"
#include "feature_table.h"
#include "renzyme_search.h"


SEQUENCES  *sequences = NULL;

/* initialise a SEQUENCE structure */
 
SEQUENCE *init_sequence (void) {

  SEQUENCE  *s; 
  if ( NULL == ( s = (SEQUENCE* ) xmalloc ( sizeof (SEQUENCE) )))
      return NULL;
  s->seq = NULL;
  s->name = NULL;
  s->source = NULL;
  s->feature_table = NULL;
  s->id = 0;
  s->parent_id = 0;
  s->offset = 0;
  s->length = 0;
  s->start = 0;
  s->end = 0;
  s->type = 0;
  s->direction = 0;
  s->genetic_code = 0;
  s->group_id = 0;
  s->group_parent_id = 0;
  s->left_end = NULL;
  s->right_end = NULL;
  s->cut_site_1 = NULL;
  s->cut_site_2 = NULL;
  s->sites = NULL;
  s->fragment = NULL;
  return s;
}

/**
 * free a SEQUENCE structure
 */

void free_sequence (SEQUENCE *sequence) {
  if (sequence) {
    if (sequence->seq) xfree(sequence->seq);
    if (sequence->name) xfree(sequence->name);
    if (sequence->source) xfree(sequence->source);   
    if (sequence->feature_table) free_feature_table(sequence->feature_table);
    if (sequence->left_end) free_site_hang (sequence->left_end);
    if (sequence->right_end) free_site_hang (sequence->right_end);
    if (sequence->cut_site_1) free_site_hang (sequence->cut_site_1);
    if (sequence->cut_site_2) free_site_hang (sequence->cut_site_2);
    if (sequence->fragment) free_sequence (sequence->fragment);
    if (sequence->sites) free_sites (sequence->sites);
    xfree (sequence);
  }
}

/* initialise a SEQUENCES structure */
SEQUENCES *init_sequences (void) {

    SEQUENCES *ss;
    SEQUENCE **s;
    
    if ( NULL == (ss = (SEQUENCES *) xmalloc ( sizeof (SEQUENCES)) ) ) 
	goto err;
    if ( NULL == (s = (SEQUENCE **) xmalloc (sizeof (SEQUENCE*))))
	goto err;
    ss->sequence = s;
    ss->num_seq = 0;
    return  ss;
 err:
    if (ss) xfree (ss);
    if (s) xfree (s);
    return NULL;
}

/* realloc the sequence array in the sequences structure */
SEQUENCES *realloc_sequences ( SEQUENCES *ss, int num_seq) {

    SEQUENCE  **s;
    s = ss->sequence;
    if ((NULL == ( s = (SEQUENCE **)xrealloc (s, sizeof(SEQUENCE *)*(num_seq+1)))))
	return NULL;
    ss->sequence = s;
    ss->num_seq = num_seq;
    return ss;
}

/* free a SEQUENCES structure */
void free_sequences ( SEQUENCES *ss ) {
    
    int i;
    if (ss) {
	for ( i = 0; i < ss->num_seq; i++ ) {
	    free_sequence (ss->sequence[i]);
	}
    }
    xfree (ss);
}

SITE_HANG *init_site_hang (void) {

    SITE_HANG *sh;
    
    if ( NULL == (sh = (SITE_HANG* ) xmalloc (sizeof (SITE_HANG))))
	return NULL;
    sh->pos = 0;
    sh->len = 0;
    sh->hang = NULL;
    return sh;
}

void free_site_hang (SITE_HANG *sh) {
    
    if (sh->hang) xfree (sh->hang);
    free (sh);
}

char *check_unique_name ( char *identifier) {

    int i, num_seq;
    static int unique_name = 0;

    num_seq = sequences->num_seq;
    for (i = 0; i < num_seq; i++) {
	if (strcmp(identifier, GetSequenceName(i)) == 0) {
	    if (NULL == (identifier = (char *)realloc(identifier, 
						    (strlen(identifier)+10))))
		return NULL;
	    sprintf(identifier, "%s#%d", identifier, unique_name++);
	}
    }
    return identifier;
}

int write_sequence_to_sequences (FEATURE_TABLE *feature_table, 
				 char *seq, 
				 char *identifier, 
				 int seq_len)
{
    SEQUENCE *s;
    int num_seq;
       
    num_seq = sequences->num_seq;
    if (NULL == ( s = init_sequence())) return -1;
    /*if (NULL == (s->seq = (char *)xmalloc((seq_len+1)*sizeof(char))))	
	return -1;
	strcpy (s->seq, seq);*/
    s->seq = seq;
    s->name = strdup (identifier);
    s->length = seq_len;
    s->start = 1;
    s->end = seq_len;
    s->type = get_seq_type (seq, seq_len);
    s->feature_table = feature_table;
    s->id = num_seq;
 
    if (num_seq != 0) {
	sequences = realloc_sequences (sequences, num_seq);
    }
    sequences->sequence[num_seq] = s;
    sequences->num_seq ++;
   
    add_reg_seq(num_seq);

    return 0;
}

int add_sequence_to_sequences (SEQUENCE *s) {

    int num_seq, err;

    num_seq = sequences->num_seq;
    if (num_seq != 0) {
	sequences = realloc_sequences (sequences, num_seq);
    }
    extend_unique_name (s); 
    s->id = num_seq;
    sequences->sequence[num_seq] = s;
 
    sequences->num_seq ++;
    err = add_reg_seq(num_seq);
    if (err == -1) return -1;
    return num_seq;
}

/* realloc the sequence buffer when necessary */
int realloc_sequence_editor(char **array,
		     int *array_size,
		     int increment)
{
    (*array_size) += increment;

    if (NULL == ((*array) = (char *)xrealloc((*array), 
					     *array_size * sizeof(char)))) {
	return -1;
    }

    return 0;
}


void write_sequence_editor(char *line,
		    char **seq,
		    int *seq_len,
		    int *buf_size)
{
    int j;
    int increment = 50000;
    
    for (j = 0; j < MAX_SEQ_LINE && line[j]; j++) {
	if ( isalpha ( (int) line[j]) || (int) line[j] == '-') {
	    if ( *seq_len >= (*buf_size)) {
		realloc_sequence_editor(seq, buf_size, increment);
	    }
	    (*seq)[*seq_len] = line[j];
	    *seq_len += 1;
	    (*seq)[*seq_len] = 0;
	}
    }
} 


/* return 0: OK; return -1: ERROR; return 1: to read next ID */

FEATURE_TABLE *read_embl_format_seq ( FILE *fp, char *line, char **seq, int *seq_len) 
{
    int i;
    int num_entry = 0;
    int sequence_line = 0;
    char *entry_loca_qua;
    char *last_line;
    int buf_size = 0;
    FEATURE_TABLE *feature_table;

    **seq = 0;
    *seq_len = 0; 
    feature_table = init_feature_table();

 next_entry: 
    while ((strncmp("//",line,2)) ) {
	if(!strncmp(&line[0], "FT", 2)) {
	    i = 0;
	    while (isgraph(line[5+i]) && i < 20) {
		i++;
	    }
	    if (i != 0) {		
		entry_loca_qua = strdup(get_ft_entry(fp, line, &last_line));
		num_entry ++;
		if (num_entry != 1)	
		feature_table = realloc_feature_table(feature_table, num_entry);
		
		feature_table->entry[num_entry - 1] = parse_ft_entry(entry_loca_qua);
		feature_table->num_entry = num_entry;

		line = last_line;
		
	       	goto next_entry;
	    }
	    feature_table->num_entry = num_entry; 
	} else if ( 0 == (strncmp("SQ", line, 2))) {
	    sequence_line = 1;
	} else if (sequence_line == 1) {
	    line[strlen(line)] = 0;	    
	    write_sequence_editor(line, seq, seq_len, &buf_size);
	}
	fgets(line, MAX_SEQ_LINE, fp);
    }
   
    return feature_table; 
}

int check_sequences_list (FILE *fp, char *line, char *iden) {

    int i, num_seq;
   
    num_seq = sequences->num_seq;
    for (i = 0; i < num_seq; i++) {
	if (!strcmp(iden, sequences->sequence[i]->name)) {   
	    while (strncmp("//", line, 2) ) {
		fgets(line, MAX_SEQ_LINE, fp);
	    }
	    return 1;
	} 
    }
    return 0;
}

int read_multiple_embl_format_sequence ( FILE *fp ) {

    char line[MAX_SEQ_LINE] = { 0 };
    char *identifier;
    int err;
  
    while ( fgets( line,sizeof(line), fp ) != NULL ) {
	char *identifier;
	char *seq = 0;
	int seq_len = 0;
	FEATURE_TABLE *feature_table;

	if (NULL == (identifier = (char *)xmalloc((MAX_SEQ_LINE+1)*sizeof(char))))
	    goto error;
	
	if ( 0 == (strncmp("ID", line, 2)) ) {
	    sscanf(line, "ID %20s\n", identifier);
	    /* check if already in sequences list */
	    err = check_sequences_list (fp, line, identifier);
	    if (err == 0) { 
		if (NULL == (seq = (char *)xmalloc((MAX_SEQ_LINE+1)*sizeof(char))))	
		    goto error;
		feature_table = read_embl_format_seq (fp, line, &seq, &seq_len);
		write_sequence_to_sequences (feature_table, seq, identifier, seq_len);
	    } else {
		if (strncmp("//", line, 2) ) {
		   fgets(line, MAX_SEQ_LINE, fp); 
		}
	    }
	} else {
	    printf ("read_sequence error: sequences in unknown format.\n");
	    return -2;
	}
    }
    return 0;
 error:
    if (identifier) xfree(identifier);
    return -1;
}

int read_file_to_sequences(char *fname) {
  
    FILE *pr;

    pr = fopen(fname, "r");
    if (pr == NULL){
	printf("%s cannot be opened\n", fname);
	return -1;
    }
    if (sequences == NULL) {
	if (NULL == (sequences = init_sequences( ))) return -1;
    }    
    if ( read_multiple_embl_format_sequence ( pr ) != 0) return -1; 
    
    return 0;
}

int ReadFileInfo (ClientData clientData,
		  Tcl_Interp *interp, 
		  int argc, 
		  char **argv)
{
    int err;
    static int initialised_seq_reg = 0;

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], " filename\"", (char*)NULL);
        return TCL_ERROR;
    }

    /* initialise sequence registration if first time thro' */    
    if (!initialised_seq_reg) {
	seq_register_init(interp);
	initialised_seq_reg = 1;
    }
    
    err = read_file_to_sequences(argv[1]);
    vTcl_SetResult(interp, "%d", err);    
    return TCL_OK;
}

int GetSeqIden (ClientData clientData,
		Tcl_Interp *interp, 
		int argc, 
		char **argv)
{
    char *name;
    int num_seq, i;
   
    if (argc != 1) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0],  (char*)NULL);
        return TCL_ERROR;
    }

    num_seq = sequences->num_seq;
    
    for (i = 0; i < num_seq; i++) {
      name = sequences->sequence[i]->name;
      Tcl_AppendElement (interp, name);
    }
    return TCL_OK;    
}

int GetSequenceIdByName (char *name) {

    int i, num_seq, name_len;

    num_seq = sequences->num_seq;
    for (i = 0; i < num_seq; i++) {
	name_len = strlen (sequences->sequence[i]->name);
	if (!strncmp(sequences->sequence[i]->name, name, name_len)) {
	    return sequences->sequence[i]->id;
	}
    }
    return -1;
}

char *GetSequenceSeq (int seq_id) {

    int i, num_seq;
    num_seq = sequences->num_seq;
    for (i = 0; i < num_seq; i++) {
	if (seq_id == sequences->sequence[i]->id)
	    return sequences->sequence[seq_id]->seq;
    }
    return NULL;
}

int GetSequenceLength (int seq_id) {

    int i, num_seq;
    num_seq = sequences->num_seq;
    for (i = 0; i < num_seq; i++) {
	if (seq_id == sequences->sequence[i]->id)
	  return sequences->sequence[seq_id]->length;  
    }
    return 0;
}

FEATURE_TABLE *GetSequenceFt (int seq_id) {

    int i, num_seq;

    num_seq = sequences->num_seq;
    for (i = 0; i < num_seq; i++) {
	if (seq_id == sequences->sequence[i]->id)
	    return sequences->sequence[seq_id]->feature_table;
    }
    return NULL;
}

int GetSequenceType (int seq_id) {
    
    int i, num_seq;
    num_seq = sequences->num_seq;
    for (i = 0; i < num_seq; i++) {
	if (seq_id == sequences->sequence[i]->id)
	    return sequences->sequence[seq_id]->type;
    }
    return 0;
}

char *GetSequenceName (int seq_id) {
    
    int i, num_seq;
    num_seq = sequences->num_seq;
    for (i = 0; i < num_seq; i++) {
	if (seq_id == sequences->sequence[i]->id)
	    return sequences->sequence[seq_id]->name;  
    }
    return NULL;   
}

int GetSequenceNums (void) {
    
    return sequences->num_seq;
}

SEQUENCE *GetSequencesSequence (int seq_id) {

    int i, num_seq;
    num_seq = sequences->num_seq;
    for (i = 0; i < num_seq; i++) {
	if (seq_id == sequences->sequence[i]->id)
	    return sequences->sequence[i];  
    }
    return NULL;   
}

/* make a copy for editing sequence in sequences array */
int make_copy_for_editor (int seq_id) {

    int num_seq;
    num_seq = sequences->num_seq;
    
    sequences = realloc_sequences (sequences, num_seq);    
    sequences->sequence[num_seq] = copy_sequence (sequences->sequence[seq_id]);
    extend_unique_name (sequences->sequence[num_seq]);
    sequences->num_seq ++;
    return  num_seq;
}

int GetIdFromName(ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, char *argv[]) {

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "sequence_identifier\"", (char*)NULL);
        return TCL_ERROR;
    }    
    vTcl_SetResult(interp, "%d", GetSequenceIdByName(argv[1]));
    return TCL_OK;
}

int GetNameFromId(ClientData clientData, 
		  Tcl_Interp *interp, 
		  int argc, char *argv[]) {
    int seq_id;

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "id\"", (char*)NULL);
        return TCL_ERROR;
    }    
    seq_id = atoi (argv[1]);
    vTcl_SetResult(interp, "%s", sequences->sequence[seq_id]->name);
    return TCL_OK;
}

/*
 * return the position of sequence in sequences array from a given unique id num
 */
int GetSequenceNum(int id) {

    int i, num_seq;
    
    num_seq = sequences->num_seq;
    for (i = 0; i < num_seq; i++) {
	if (sequences->sequence[i]->id == id) {
	    return i;
	}
    }
    return -1;
}

int GetSequenceId(int num) {

    return sequences->sequence[num]->id; 
}

void remove_sequence_from_sequence_list (int seq_id) {

    int i, num_seq;
    num_seq = sequences->num_seq;
    
    for (i = 0; i < num_seq; i++) {
	if (seq_id == sequences->sequence[i]->id) {
	    memmove(&sequences->sequence[i], &sequences->sequence[i+1], 
		    (num_seq - i - 1) * sizeof(sequences->sequence[i]));
	    sequences->num_seq--;
	    i--;
	    num_seq--;
	}
    }
}

int sequence_save (SEQUENCE *s) {

    int i, num_seq;

    num_seq = sequences->num_seq;

    for (i = 0; i < num_seq; i++) {
	if (s->id == sequences->sequence[i]->id) {
	    free_sequence (sequences->sequence[i]);
	    sequences->sequence[i] = copy_sequence (s);
	    if (sequences->sequence[i] == NULL) return -1;
	}
    }
    return 0;
}

int save_change_to_sequence (int seq_id, SEQUENCE *s) {

    int i, num_seq;
    num_seq = sequences->num_seq;

    for (i = 0; i < num_seq; i++) {
	if (seq_id == sequences->sequence[i]->id) {
	    free_sequence (sequences->sequence[i]);
	    sequences->sequence[i] = copy_sequence (s);
	    if (sequences->sequence[i] == NULL) return -1;
	}
    }
    return 0;
}

int Read_sequence_Init(Tcl_Interp *interp) {

    Tcl_CreateCommand(interp, "get_sequences_iden", GetSeqIden,
		      (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);   
    Tcl_CreateCommand(interp, "get_id_from_name", GetIdFromName,
		      (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "get_name_from_id", GetNameFromId,
		      (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
    if (Itcl_RegisterC(interp, "read_file_info", ReadFileInfo, NULL, NULL) != TCL_OK) {
	return TCL_ERROR;
    }
    return TCL_OK;
}
