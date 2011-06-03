#include <string.h>
#include <math.h>
#include <float.h>

/* added for parsing features */
#include <stdio.h>
#include "codon_content.h"
/* added for parsing features */

#include "seq_results.h"
#include "seq_raster.h"
#include "misc.h"
#ifdef USE_SEQLIB
#include "seqlib_globals.h"
#endif
#ifdef USE_SRS /* 11/1/99 johnt - made SRS support conditional - define USE_SRS for SRS support */
#include "seqlib_file_srs.h"
#endif
#ifdef USE_SEQLIB
#include "seqlib_file.h"
#endif
#include "sequence_formats.h"
#include "dna_utils.h"
#include "genetic_code.h"
#include "scramble.h"

static SubSeq *seqs = NULL;
static int num_seqs = 0;
static int horizontal = -1;
static int vertical = -1;
static int active_seq = -1;

int Set_Seqs(int seq_num,
	     int direction,
	     int l_index,
	     char *name,
	     char *sequence, 
	     int seq_structure, 
	     int seq_type,
	     Featcds **key_index,
             char *identifier);

/* check if adding same sequence name, if are, remove previous seq */
void RemoveDuplicateSeq(Tcl_Interp *interp,
			char *name)
{
    int n_seqs;
    int i;
    /*
     * need to check if entry_name already exists in the database. If it does
     * then must delete current sequence (and all its results) and add the
     * new entry_name
     */
    n_seqs = NumSequences();
    for (i = 0; i < n_seqs; i++) {
	if (strcmp(name, GetSeqName(i)) == 0) {
	    verror(ERR_WARN, "RemoveDuplicateSeq", 
		   "%s already exists in. Removing previous sequence and adding new sequence\n", name); 
	    DeleteSequence(interp, i--);
	    n_seqs--;
	}
    }
}

/*
 * see if sequence already exists ie name and sequence are identical
 * return seq_num if does exist else return -1
 */
int CheckSeqExists(char *name,
		   char *seq)
{
    int i;

    for (i = 0; i < num_seqs; i++) {
	if ((strcmp(seqs[i].name, name) == 0) && 
	    (strcmp(seqs[i].seq->sequence, seq) == 0)) {
	    return i;
	}
    } 
    return -1;
}

/*
 * add a new sequence to sequence list
 * add a new sequence to registration scheme
 * update the sequence list box
 * return the sequence number just created
 */
int AddSequence(Tcl_Interp *interp,
		int direction,
		int library,
		char *entry,
		char *sequence, 
		int seq_structure,
		int seq_type,
		Featcds **key_index,
                char *identifier)
{
    int seq_num;

#ifdef FIXME
    /* 
     * don't do this yet - this ends up removing plots when I do translations
     * twice. All needs careful thinking about!
     */
    /* check if adding same sequence name, if are, remove previous seq */
    RemoveDuplicateSeq(interp, entry);
#endif

    seq_num = CheckSeqExists(entry, sequence);
    if (seq_num > -1) {
	xfree(sequence);
	return seq_num;
    }

    /* add seq_num to seqs array */
    seq_num = SeqCreate();
    if (-1 == Set_Seqs(seq_num, direction, library, entry, sequence, 
		       seq_structure, seq_type, key_index, identifier)) {
	Delete_Seq(seq_num);
	return -1;
    }
    /* add cur_seq_num to registration scheme */
    if (-1 == add_reg_seq(seq_num)) {
	Delete_Seq(seq_num);
	return -1;
    }

    vmessage("Added sequence %s\n", entry);
    return seq_num;
}

/*
 * delete a sequence from sequence list
 * delete a sequence from registration scheme
 * update the sequence list box
 */
void DeleteSequence(Tcl_Interp *interp,
		   int seq_num)
{

    /* delete seq_num from registration scheme */
    delete_reg_seq(seq_num);

    /* Do this after delete_reg_seq which need seq_id */
    /* delete seq_num from seqs array */
    Delete_Seq(seq_num);
    /*    
	  {
	  char cmd[100];
	  sprintf(cmd, "SipSeqShutdown %d\n", GetSeqId(seq_num));
	  Tcl_Eval(interp, cmd);
	  Tcl_VarEval(interp, "sip_sequence_list_update", NULL);
	  }
	  */
}

int CreateSeqid(void)
{
    static int id = 0;
    return id++;
}

/*
 * add new sequence range into seqs array
 */
int Set_SubSeqs(int seq_id,
		 int seq_num,
		 int start,
		 int end, 
		 char *name,
		 Featcds **key_index,
	         char *identifier)
{
    int num;
    if (-1 == (num = GetSeqNum(seq_id)))
        return -1;

    seqs[seq_num].seq = seqs[num].seq;
    seqs[seq_num].start = start;
    seqs[seq_num].end = end;
    seqs[seq_num].name = name;
    seqs[seq_num].seq_id = CreateSeqid();
    seqs[seq_num].key_index = key_index;

    if (NULL == (seqs[seq_num].identifier = (char *)xmalloc((strlen(identifier)+1) * 
				 sizeof(char))))
	return -1;

    strcpy(seqs[seq_num].identifier, identifier);

    /* increment counter */
    seqs[seq_num].seq->count++;
    return 0;
}

/*
 * add a new sequence range to sequence list
 * add a new sequence range to registration scheme
 * update the sequence list box
 * return the sequence number just created
 */
static int AddSubSequence(Tcl_Interp *interp,
			  int seq_id,
			  int start,
			  int end,
			  char *name)
{
    int seq_num;
    Featcds **key_index = NULL;
    char *identifier = " ";

    /* add seq_num to seqs array */
    seq_num = SeqCreate();
    if (-1 == (Set_SubSeqs(seq_id, seq_num, start, end, name, key_index, identifier)))
       return -1;

    /* add cur_seq_num to registration scheme */
    if (-1 == add_reg_seq(seq_num)) {
	Delete_Seq(seq_num);
	return -1;
    }
    /* Tcl_Eval(interp, "sip_sequence_list_update"); */
    return seq_num;
}

/*
 * create or extend array containing the sequence info created each time a
 * new sequence is read in using the 'get sequences', 'get horizontal seq'
 * or 'get vertical seq' options
 * returns the next free index
 */
int SeqCreate()
{
    static int total_size = 0;
    static int num_results = 0;
    int increment = 5;

    if (total_size <= num_seqs) { 
	num_results = num_seqs + increment;
	if (NULL == (seqs = (SubSeq *)xrealloc(seqs, 
					   num_results * sizeof(SubSeq)))) {
	    xfree(seqs);
	    return -1;
	}
	total_size = num_results;
    }
    seqs[num_seqs].seq = NULL;
    seqs[num_seqs].identifier = NULL;
    num_seqs++;
    return (num_seqs - 1);
}

/*
 * add new sequence into seqs array
 */
int Set_Seqs(int seq_num,
	     int direction,
	     int l_index,
	     char *name,
	     char *sequence, 
	     int seq_structure, 
	     int seq_type,
	     Featcds **key_index,
             char *identifier)
{

    /* only work out type when it is unknown. */
    if (seq_type == 0) {
	if (0 == (seq_type = get_seq_type(sequence, strlen(sequence))))
	    return -1;
    }

    if (NULL == (seqs[seq_num].seq = (SeqInfo *)xcalloc(1, sizeof(SeqInfo))))
	return -1;

    if (NULL == (seqs[seq_num].seq->name = (char *)xmalloc((strlen(name)+1) * 
				 sizeof(char))))
	return -1;

    if (NULL == (seqs[seq_num].identifier = (char *)xmalloc((strlen(identifier)+1) * 
				 sizeof(char))))
	return -1;

    strcpy(seqs[seq_num].identifier, identifier);
    strcpy(seqs[seq_num].seq->name, name);
    seqs[seq_num].seq->library = l_index;
    seqs[seq_num].seq->sequence = sequence;
    seqs[seq_num].seq->seq_len = strlen(sequence);
    seqs[seq_num].seq->type = seq_type;
    seqs[seq_num].seq->seq_id = CreateSeqid();
    seqs[seq_num].seq->count = 1;
    seqs[seq_num].seq->raster = NULL;
    seqs[seq_num].seq->seq_structure = seq_structure;

    seqs[seq_num].start = 1;
    seqs[seq_num].end = seqs[seq_num].seq->seq_len;
    seqs[seq_num].name = strdup(seqs[seq_num].seq->name);
    seqs[seq_num].seq_id = seqs[seq_num].seq->seq_id;

  /* added for parsing feature tables */
    if(key_index != NULL){
	seqs[seq_num].key_index = key_index;
    } else {
	seqs[seq_num].key_index = NULL;
    }

    /* added for parsing feature tables */
    
    if (direction == HORIZONTAL) {
	/* sip */
	horizontal = seq_num;
    } else if (direction == VERTICAL) {
	vertical = seq_num;
    } else {
	/* nip */
	active_seq = seq_num;
    }
    return 0;
}

/*
 * delete sequence from seqs array
 */
void Delete_Seq(int seq_num)
{
    /* this happens if Set_Seqs fails but still need to reduce num_seqs */
    if (!seqs[seq_num].seq) {
	num_seqs--;
	return;
    }
    seqs[seq_num].seq->count--;
    /* if last substructure referencing seq struct, delete seq struct */
    if (seqs[seq_num].seq->count == 0) {

	if (seqs[seq_num].seq->name)
	    xfree(seqs[seq_num].seq->name);
	if (seqs[seq_num].seq->sequence)
	    xfree(seqs[seq_num].seq->sequence);
	if (seqs[seq_num].seq->raster)
	    xfree(seqs[seq_num].seq->raster);
	if (seqs[seq_num].seq)
	    xfree(seqs[seq_num].seq);
	if (seqs[seq_num].name)
	    xfree(seqs[seq_num].name);
	    
	free_key_index(seqs[seq_num].key_index);

	if (seqs[seq_num].identifier)
	    xfree(seqs[seq_num].identifier);    
    }
    /* only move if not the last sequence in array */
    if (seq_num < num_seqs - 1) {
	memmove(&seqs[seq_num], &seqs[seq_num + 1], 
		(num_seqs - seq_num - 1) * sizeof(SubSeq));
    }

    num_seqs--;
    /* change active seq numbers accordingly */

    /* sip */
    if (horizontal > -1 && seq_num < horizontal) 
	horizontal--;
    else if (seq_num == horizontal)
      /*horizontal = -1;*/
     horizontal = 0;
 
    if (vertical > -1 && seq_num < vertical)
	vertical--;
    else if (seq_num == vertical)
	vertical = -1;

    /* nip */
    if (active_seq > -1 && seq_num < active_seq) 
	active_seq--;
    else if (seq_num == active_seq)
	active_seq = -1;
}

/*
 * return the position of seq in seqs array from a given unique id num
 */
int GetSeqNum(int id)
{
    int i;

    for (i = 0; i < num_seqs; i++) {
	if (seqs[i].seq_id == id) {
	    return i;
	}
    }
    return -1;
}

/*
 * return the unique id num of seq at position seq_num in seqs
 */
int GetSeqId(int seq_num)
{
    if (seq_num < num_seqs && seq_num > -1) {
	return seqs[seq_num].seq_id;
    } 
    return -1;
}

/*
 * return the number of stored sequences
 */
int NumSequences()
{
    return num_seqs;
}

/*
 * set of functions that return info on a particular sequence in the global 
 * seqs array
 */

/* added for parsing features */
char *GetSeqKeyIndexCds(int seq_num, int idx)
{

    BasePos *current;
    char tmp[1024];
    char *key_index_subcds;

    if (NULL == (key_index_subcds = (char *)xmalloc(sizeof(char ))))
	return NULL;
    sprintf(tmp, "CDS %3d %2s ",idx, seqs[seq_num].key_index[0][idx].type_loca);
     if (NULL == (key_index_subcds = (char *)xrealloc(key_index_subcds,
                 (strlen(tmp)+1)*sizeof(char ))))   
	return NULL;    
    strcpy(key_index_subcds, tmp);

    for (current = seqs[seq_num].key_index[0][idx].loca; current!=NULL; 
                                     current=current->next){
	sprintf(tmp, " %2s %d..%d ", current->type_range,
		current->start_pos,current->end_pos );
	if (NULL == (key_index_subcds=(char *)xrealloc(key_index_subcds,
                    (strlen(tmp)+1)*sizeof(char ))))   
     return NULL;       
	strcat(key_index_subcds, tmp);
    }
    return key_index_subcds;
}

int GetSeqNumberCds(int seq_num)
{
    if (seqs[seq_num].key_index != NULL){
	return seqs[seq_num].key_index[0]->id;
    }
    else  return 0;
}

char *GetSeqProteinId(int seq_num, int k)
{ 
    int i;
    for(i=0; i < number_quas; i++) {
	if(seqs[seq_num].key_index[0][k].qualifier[i] != NULL 
                && !strncmp(&seqs[seq_num].key_index[0][k].qualifier[i][0],
                "/protein_id", 11))
	    return seqs[seq_num].key_index[0][k].qualifier[i];  
    }
    return NULL;/* no protein id */ 
}

char *GetSeqIdentifier(int seq_num)
{
    return seqs[seq_num].identifier;
}

char *GetSeqCdsExpr(int seq_num, int i)
{
    return seqs[seq_num].key_index[0][i].cdsexpr;
}

Featcds **GetSeqKeyIndex(int seq_num)
{
    return seqs[seq_num].key_index;
}

/* added for parsing features */

char *GetSeqName(int seq_num)
{
    return seqs[seq_num].name;
}

char *GetSeqBaseName(int seq_num)
{
  char *cp;
  if (cp = strrchr(seqs[seq_num].name, '/'))
    return cp+1;
  else
    return seqs[seq_num].name;
}

char *GetSeqSequence(int seq_num)
{
    return seqs[seq_num].seq->sequence;
}

int GetSeqLength(int seq_num)
{
    return seqs[seq_num].seq->seq_len;
}

int GetSeqType(int seq_num)
{
    return seqs[seq_num].seq->type;
}

int GetSeqLibrary(int seq_num)
{
    return seqs[seq_num].seq->library;
}

int GetSeqStructure(int seq_num)
{
    return seqs[seq_num].seq->seq_structure;
}

void SetSeqStructure(int seq_num,
		     int seq_structure)
{
    seqs[seq_num].seq->seq_structure = seq_structure;
}

char *GetSeqLibraryName(int seq_num)
{
    if (seqs[seq_num].seq->library == -1) {
	return "PERSONAL";
    } else {
#ifdef USE_SRS
	if (get_seqlib_lib() == SRS51_LIB) {
	    return GetSeqLibName_SRS(seqs[seq_num].seq->library);
	}
	else
#endif
#ifdef USE_SEQLIB
	{
	    return GetSeqLibName(seqs[seq_num].seq->library);
	}
#endif
    }
    return NULL;
}

char *GetRaster(int seq_num)
{
    return (seqs[seq_num].seq->raster);
}

void SetRaster(int seq_num,
	      char *raster)
{

    seqs[seq_num].seq->raster = raster;

}

int GetSeqDirection(int seq_num)
{
    if (seq_num == horizontal) {
	return HORIZONTAL;
    } else if (seq_num == vertical) {
	return VERTICAL;
    } else {
	return NEITHER;
    }
}

int GetParentalSeqId(int seq_num)
{
    return (seqs[seq_num].seq->seq_id);
}

char *GetParentalSeqName(int seq_num)
{
    return (seqs[seq_num].seq->name);
}

int GetSubSeqLength(int seq_num)
{
    return (seqs[seq_num].end - seqs[seq_num].start + 1);
}

int GetSubSeqStart(int seq_num)
{
    return (seqs[seq_num].start);
}

int GetSubSeqEnd(int seq_num)
{
    return (seqs[seq_num].end);
}

char *GetSubSeqName(int seq_num)
{
    return (seqs[seq_num].name);
}

int GetSeqIdFromName(char *name) {
    int i;

    for (i = 0; i < num_seqs; i++) {
	if (strcmp(seqs[i].name, name) == 0) {
	    return seqs[i].seq_id;
	}
    }
    return -1;
}

/*
 * set the "active" sequence in a given direction
 */
int
Set_Active_Seq(int seq_num,
	       int direction)
{
    /* nip */
    if (direction == -1)
	active_seq = seq_num;

    if (direction == HORIZONTAL) {
	horizontal = seq_num;
    } else if (direction == VERTICAL) {
	vertical = seq_num;
    } else {
	return -1;
    }
    /* 
     * HACK to deal with h and v set to be the same value, default to the
     * latest change and set the other to be -1
     */
    /*   if (horizontal == vertical) {
	if (direction == HORIZONTAL){ 
	    vertical = -1;
	} else if (direction == VERTICAL) {
	    horizontal = -1;
	}
	}*/

     /*yy:Set horizontal sequence as active sequence  in sequence array
     */
    if (horizontal == vertical){
	if (direction == HORIZONTAL) { 
	    vertical = -1;
	}else if (direction == VERTICAL) { 
	    if( NumSequences() > 1 ){	       
	    horizontal = 0;                   
	    }else if (NumSequences()==1) {
		vertical = -1;
		horizontal = 0;
	    }
	}
    }

    return 0;
}

/*
 * set of functions to return information about the active horizontal or
 * vertical sequences
 */
int GetActiveSeqNumber(int direction)
{
    if (direction == HORIZONTAL && horizontal > -1)
	return horizontal;
    else if (direction == VERTICAL && vertical > -1)
	return vertical;
    
    return -1;
}

/*
 * set active range for sequence
 * 
 */
int SetRange(Tcl_Interp *interp,
	     int seq_id,
	     int start,
	     int end)
{
    char *name;
    int seq_num;
    static int count = 1;

    seq_num = GetSeqNum(seq_id);
    if (NULL == (name = (char *)xmalloc((strlen(GetSeqName(seq_num))+20) 
					 * sizeof(char))))
	return -1;

    sprintf(name, "%s_s%d", GetSeqName(seq_num), count++);

    return (AddSubSequence(interp, seq_id, start, end, name));
}

/*
 * copy range of sequence
 * add new sequence name to sequence manager
 */
int CopyRange(Tcl_Interp *interp,
	      int seq_id,
	      int start,
	      int end)
{
    int seq_num = GetSeqNum(seq_id);
    char *seq1 = GetSeqSequence(seq_num);
    char *seq2;
    char *name;
    int length = end - start + 1;
    int new_seq_num;
    char *parental_name, *child_name;
    static int count = 1;


    if (NULL == (seq2 = (char *)xmalloc((length+2) * sizeof(char))))
	return -1;
     
    strncpy(seq2, &seq1[start-1], end - start + 1);
    seq2[end - start + 1] = '\0';

    parental_name = GetParentalSeqName(seq_num);
    child_name = GetSeqName(seq_num);

    if (NULL == (name = (char *)xmalloc((strlen(parental_name)+20) * 
					sizeof(char))))
	return -1;
    sprintf(name, "%s_n%d", parental_name, count++);
    if (-1 == (new_seq_num = AddSequence(interp, -1, GetSeqLibrary(seq_num), 
					 name, seq2, GetSeqStructure(seq_num), 
					 GetSeqType(seq_num), NULL, " ")))
	return -1;
    xfree(name);

    /* don't think I need to deal with subsequences in copying a sequence. */
#if 0
    if (strcmp(parental_name, child_name) != 0) {
	/* sub-sequence */
	/* 
	 * need to get seq num from seq_id instead of using seq_num incase
	 * AddSequence has deleted duplicate names
	 */
	start = GetSubSeqStart(GetSeqNum(seq_id));
	end = GetSubSeqEnd(GetSeqNum(seq_id));

	if (NULL == (name = (char *)xmalloc((strlen(child_name)+3) * 
					    sizeof(char))))
	    return -1;

	sprintf(name, "%s_n%d", child_name, count++);
	if (-1 == (AddSubSequence(interp, GetSeqId(new_seq_num), 
				  start, end, name)))
	    return -1;
    }
#endif
    return 0;
}

/*
 * complement sequence
 * add new sequence name to sequence manager
 */
int ComplementSeq(Tcl_Interp *interp,
		  int seq_num)
{
    char *seq1 = GetSeqSequence(seq_num);
    char *seq2;
    char *name;
    int length = GetSeqLength(seq_num);
    int seq_id = GetSeqId(seq_num);
    int new_seq_num, start, end;
    char *parental_name, *child_name;

    if (NULL == (seq2 = (char *)xmalloc((length+1) * sizeof(char))))
	return -1;

    memcpy(seq2, seq1, length);
    (void) complement_seq(seq2, length);
    seq2[length] = '\0';

    parental_name = GetParentalSeqName(seq_num);
    child_name = GetSeqName(seq_num);

    if (NULL == (name = (char *)xmalloc((strlen(parental_name)+3) * 
					sizeof(char))))
	return -1;
    sprintf(name, "%s_c", parental_name);
    if (-1 == (new_seq_num = AddSequence(interp, -1, GetSeqLibrary(seq_num), 
					 name, seq2, GetSeqStructure(seq_num), 
					 GetSeqType(seq_num), NULL, " ")))
	return -1;

    xfree(name);

    if (strcmp(parental_name, child_name) != 0) {
	/* sub-sequence */
	/* 
	 * need to get seq num from seq_id instead of using seq_num incase
	 * AddSequence has deleted duplicate names
	 */
	start = GetSubSeqStart(GetSeqNum(seq_id));
	end = GetSubSeqEnd(GetSeqNum(seq_id));

	if (NULL == (name = (char *)xmalloc((strlen(child_name)+3) * 
					    sizeof(char))))
	    return -1;

	sprintf(name, "%s_c", child_name);
	if (-1 == (AddSubSequence(interp, GetSeqId(new_seq_num), 
				  length - end + 1, 
				  length - start + 1, name)))
	    return -1;
    }
    return 0;
}

/*
 * convert t to u or u to t
 * add new sequence name to sequence manager
 */
int RnaSeq(Tcl_Interp *interp,
	   int seq_num)

{ 
    char *seq1 = GetSeqSequence(seq_num);
    int seq_id = GetSeqId(seq_num);
    char *seq2;
    char *name;
    int length = GetSeqLength(seq_num);
    int i;
    char *parental_name, *child_name;
    int new_seq_num, start, end;

    if (NULL == (seq2 = (char *)xmalloc((length+1) * sizeof(char))))
	return -1;

    memcpy(seq2, seq1, length);
    for (i = 0; i < length; i++) {
	if (seq2[i] == 't') {
	    seq2[i] = 'u';
	} else if (seq2[i] == 'T') {
	     seq2[i] = 'U';
	} else if (seq2[i] == 'u') {
	    seq2[i] = 't';
	} else if (seq2[i] == 'U') {
	    seq2[i] = 'T';
	}
    }
    seq2[length] = '\0';

    parental_name = GetParentalSeqName(seq_num);
    child_name = GetSeqName(seq_num);

    if (NULL == (name = (char *)xmalloc((strlen(parental_name)+3) * 
					sizeof(char))))
	return -1;
    sprintf(name, "%s_r", parental_name);
    if (-1 == (new_seq_num = AddSequence(interp, -1, GetSeqLibrary(seq_num), 
					 name, seq2, GetSeqStructure(seq_num),
					 GetSeqType(seq_num), NULL, " "))) 
	return -1;
    xfree(name);

    if (strcmp(parental_name, child_name) != 0) {
	/* sub-sequence */
	/* 
	 * need to get seq num from seq_id instead of using seq_num incase
	 * AddSequence has deleted duplicate names
	 */
	start = GetSubSeqStart(GetSeqNum(seq_id));
	end = GetSubSeqEnd(GetSeqNum(seq_id));

	if (NULL == (name = (char *)xmalloc((strlen(child_name)+3) * 
					    sizeof(char))))
	    return -1;

	sprintf(name, "%s_r", child_name);
	if (-1 == (AddSubSequence(interp, GetSeqId(new_seq_num), start, end, name)))
	    return -1;
    }

    return 0;
}

/*
 * translate seq given by seq_num into reading frames given by rf (0, 1, 2) or
 * 3 for all 3.
 * add new sequence name to sequence manager
 * return the sequence number of new sequence
 */

int TranslateSeq(Tcl_Interp *interp,
		 int seq_num,
		 int rf,
		 int start,
		 int end)
{
    int i;
    char *name;
    char *dna_seq;
    char *prot_seq;
    int cnt = 0;
    int seq_id = GetSeqId(seq_num);
    int new_seq_num;
    char *ptr;
    char *parental_name, *child_name;
    char *new_name;
    int length = end - start + 1;
    static int num = 0;

#ifdef DEBUG
    printf("START translate seq %d to %d\n", start, end);
#endif
    dna_seq = GetSeqSequence(seq_num);
    if (NULL == (prot_seq = (char *)xmalloc(((length/3)+3) * sizeof(char))))
	return -1;
    if (NULL == (new_name = (char *)xmalloc(strlen(GetSeqName(seq_num))
					    * sizeof(char))))
	return -1;

    for (i = rf+start-1; i < end-2; i+=3) {
	prot_seq[cnt++] = codon_to_acid1(&dna_seq[i]);
    }
    prot_seq[cnt] = '\0';
#ifdef DEBUG
    printf("%s\n", prot_seq);
#endif
    /* 
     * special case: remove _rf123 from name before adding _rfx to end 
     */
    parental_name = GetParentalSeqName(seq_num);
    child_name = GetSeqName(seq_num);
    ptr = strstr(parental_name, "_rf123");

    if (NULL == (name = (char *)xmalloc((strlen(parental_name)+28) 
					 * sizeof(char))))
	return -1;
    if (ptr) {
	strncpy(new_name, parental_name, (ptr - parental_name));
	new_name[ptr - parental_name] = '\0';
	strcat(new_name, ptr+6);
	sprintf(name, "%s_rf%d_%d", new_name, rf+1, num);
    } else {
	sprintf(name, "%s_rf%d_%d", parental_name, rf+1, num);
    }

    /* proteins can only be LINEAR ! */
    if (-1 == (new_seq_num = AddSequence(interp, -1, GetSeqLibrary(seq_num), 
					 name, prot_seq, LINEAR, PROTEIN, NULL, " ")))
	return -1;
    xfree(name);
    xfree(new_name);

    if (strcmp(parental_name, child_name) != 0) {
	/* sub-sequence */
	/* 
	 * need to get seq num from seq_id instead of using seq_num incase
	 * AddSequence has deleted duplicate names
	 */
	start = ceil((GetSubSeqStart(GetSeqNum(seq_id))-1)/3.0 + 1);
	end = (GetSubSeqEnd(GetSeqNum(seq_id)) - rf) / 3;
	if (NULL == (name = (char *)xmalloc((strlen(child_name)+15) * 
					    sizeof(char))))
	    return -1;
	if (NULL == (new_name = (char *)xmalloc(strlen(GetSeqName(seq_num))
						* sizeof(char))))
	    return -1;

	ptr = strstr(child_name, "_rf123");
	
	if (ptr) {
	    strncpy(new_name, child_name, (ptr - child_name));
	    new_name[ptr - child_name] = '\0';
	    strcat(new_name, ptr+6);
	    sprintf(name, "%s_rf%d_%d", new_name, rf+1, num);
	} else {
	    sprintf(name, "%s_rf%d_%d", child_name, rf+1, num);
	}

	/* sprintf(name, "%s_rf%d", child_name, rf+1); */
	new_seq_num = AddSubSequence(interp, GetSeqId(new_seq_num), start, end, 
				     name);
	xfree(new_name);
    }
    num++;
    return new_seq_num;
}

/*
 * add new sequence to sequence manager, but keep as dna but with an extension
 * of _rf123 which signifies to the comparison functions that the sequence
 * is to translated into it's 3 reading frames, each of which will be used
 * in the comparison routine
 */
int TranslateTogether(Tcl_Interp *interp,
		      int seq_num)
{
    char *name;
    char *dna_seq;
    char *prot_seq;
    int seq_id = GetSeqId(seq_num);
    int new_seq_num;
    char *parental_name, *child_name;
    int start, end;

#ifdef DEBUG
    printf("START translate together \n");
#endif
    dna_seq = GetSeqSequence(seq_num);
    if (NULL == (prot_seq = strdup(dna_seq)))
	return -1;
    
    parental_name = GetParentalSeqName(seq_num);
    child_name = GetSeqName(seq_num);

    if (NULL == (name = (char *)xmalloc((strlen(parental_name)+7) * 
					sizeof(char))))
	return -1;
    sprintf(name, "%s_rf123", parental_name);
    if (-1 == (new_seq_num = AddSequence(interp, -1, GetSeqLibrary(seq_num), 
					 name, prot_seq, LINEAR, PROTEIN, NULL, " ")))
	return -1;

    xfree(name);

    if (strcmp(parental_name, child_name) != 0) {
	/* sub-sequence */
	/* 
	 * need to get seq num from seq_id instead of using seq_num incase
	 * AddSequence has deleted duplicate names
	 */
	start = GetSubSeqStart(GetSeqNum(seq_id));
	end = GetSubSeqEnd(GetSeqNum(seq_id));

	if (NULL == (name = (char *)xmalloc((strlen(child_name)+7) * 
					    sizeof(char))))
	    return -1;

	sprintf(name, "%s_rf123", child_name);
	new_seq_num = AddSubSequence(interp, GetSeqId(new_seq_num), start, end, 
				     name);
    }
    return new_seq_num;
}

/*
 * create a random sequence from either dna or protein input sequence
 * add sequence to sequence manager
 */
int ScrambleSeq(Tcl_Interp *interp,
		int seq_num)
{
    char *seq1 = GetSeqSequence(seq_num);
    int length = GetSeqLength(seq_num);
    int seq_id = GetSeqId(seq_num);
    char *seq2;
    char *name;
    time_t tim;
    int seed;
    char *parental_name, *child_name;
    int new_seq_num, start, end;
    static int num = 0;


    if (NULL == (seq2 = (char *)xmalloc((length+1) * sizeof(char))))
	return -1;

    memcpy(seq2, seq1, length);

    tim = time(NULL);
    seed = (int) tim;

    scramble_seq(seq2, length, seed);
    seq2[length] = '\0';

    parental_name = GetParentalSeqName(seq_num);
    child_name = GetSeqName(seq_num);

    if (NULL == (name = (char *)xmalloc((strlen(parental_name)+13) * 
					sizeof(char))))
	return -1;
    sprintf(name, "%s_x%d", parental_name, num);
    if (-1 == (new_seq_num = AddSequence(interp, -1, GetSeqLibrary(seq_num), 
					 name, seq2, GetSeqStructure(seq_num), 
					 GetSeqType(seq_num), NULL , " ")))
	return -1;
    xfree(name);

    if (strcmp(parental_name, child_name) != 0) {
	/* sub-sequence */
	/* 
	 * need to get seq num from seq_id instead of using seq_num incase
	 * AddSequence has deleted duplicate names
	 */
	start = GetSubSeqStart(GetSeqNum(seq_id));
	end = GetSubSeqEnd(GetSeqNum(seq_id));

	if (NULL == (name = (char *)xmalloc((strlen(child_name)+13) * 
					    sizeof(char))))
	    return -1;

	sprintf(name, "%s_x%d", child_name, num);
	if (-1 == (AddSubSequence(interp, GetSeqId(new_seq_num), start, end, name)))
	    return -1;
    }
    num++;
    return 0;
}


/*
 * rotate sequence
 * add new sequence name to sequence manager
 */
int RotateSeq(Tcl_Interp *interp,
	      int seq_num,
	      int origin)
{
    char *seq1 = GetSeqSequence(seq_num);
    char *seq2;
    char *name;
    int length = GetSeqLength(seq_num);
    int seq_id = GetSeqId(seq_num);
    int new_seq_num, start, end;
    char *parental_name, *child_name;
    static int num = 0;

    if (NULL == (seq2 = (char *)xmalloc((length+1) * sizeof(char))))
	return -1;

    memcpy(seq2, seq1, length);
    (void) rotate_seq(seq2, length, origin);
    seq2[length] = '\0';

    parental_name = GetParentalSeqName(seq_num);
    child_name = GetSeqName(seq_num);

    if (NULL == (name = (char *)xmalloc((strlen(parental_name)+13) * 
					sizeof(char))))
	return -1;
    sprintf(name, "%s_o%d", parental_name, num);
    if (-1 == (new_seq_num = AddSequence(interp, -1, GetSeqLibrary(seq_num), 
					 name, seq2, GetSeqStructure(seq_num),
					 GetSeqType(seq_num), NULL, " ")))
	return -1;
    xfree(name);

    if (strcmp(parental_name, child_name) != 0) {
	/* sub-sequence */
	/* 
	 * need to get seq num from seq_id instead of using seq_num incase
	 * AddSequence has deleted duplicate names
	 */
	start = GetSubSeqStart(GetSeqNum(seq_id));
	end = GetSubSeqEnd(GetSeqNum(seq_id));

	if (NULL == (name = (char *)xmalloc((strlen(child_name)+13) * 
					    sizeof(char))))
	    return -1;

	sprintf(name, "%s_o", child_name);
	if (-1 == (AddSubSequence(interp, GetSeqId(new_seq_num), 
				  length - end + 1, 
				  length - start + 1, name)))
	    return -1;
    }
    num++;
    return 0;
}

seq_result *seq_id_to_result(int result_id)
{
    seq_reg_info info;

    info.job = SEQ_RESULT_INFO;
    info.op = RESULT;
    info.result = NULL;
    seq_result_notify(result_id, (seq_reg_data *)&info, 0);
    if (!info.result) {
#ifdef DEBUG
	printf("NO result info\n");
#endif
	return NULL;
    }
    return (seq_result *)info.result;
}

int comparison2(void *v_result,
		int type)
{
    if (type == SEQ_PLOT_PERM || type == SEQ_PLOT_TEMP) {
	return 1;
    }
    return 0;
}

/*
 * replots each result in a raster
 */
void SeqReplotResults(Tk_Raster *raster,
		      char *raster_win,
		      int replot_all,
		      int zoom, 
		      int x0, int y0,
		      int x1, int y1)
{
    int num_elements;
    seq_result **data;
    seq_result *result;
    int num_funcs;
    int i;
    seq_reg_plot jdata;
    out_raster *output;
    RasterResult *raster_result;
    int raster_id;

#ifdef START_DEBUG
    printf("START SipReplotResults window %s x0 %d x1 %d y0 %d y1 %d \n",
	   raster_win, x0, x1, y0, y1);
#endif

    jdata.job = SEQ_PLOT;
    jdata.x0 = x0;
    jdata.x1 = x1;
    jdata.y0 = y0;
    jdata.y1 = y1;

    num_elements = seq_num_results();
    if (num_elements == 0)
	return;

    data = (seq_result **)xmalloc(num_elements * sizeof(seq_result *));
 	
    if (-1 == search_reg_data(comparison2, (void **)data, &num_funcs)) {
	xfree(data);
	return;	
    }
    if (num_funcs == 0) {
	xfree(data);
	return;
    }
    result = data[0];
    output = result->output;

    if (zoom) {
	if (TCL_OK != Tcl_VarEval(output->interp, "rasterRescaleZoom ", 
				  raster_win, NULL))
	    verror(ERR_WARN, "SeqReplotResults", "%s\n", 
		   Tcl_GetStringResult(output->interp));
    }

    Tcl_VarEval(output->interp, "GetRasterId ", raster_win, NULL);
    raster_id = atoi(Tcl_GetStringResult(output->interp));

    if (NULL == (raster_result = raster_id_to_result(raster_id))) {
	xfree(data);
	return;
    }
    if (replot_all) {
	remove_all_raster_cursors(output->interp, raster, raster_result);
	tk_RasterClear(raster);
    }

    for (i = 0; i < num_funcs; i++) {
	result = data[i];
	output = result->output;

	if (output != NULL) {
	    if (strcmp(output->raster_win, raster_win) == 0) {
		seq_result_notify(result->id, (seq_reg_data *)&jdata, 0);
	    }
	}
    }
    raster_update_cursors(raster_result, raster);
    xfree(data);
}

void SeqRasterPlotFunc(Tk_Raster *raster,
		       char *raster_win,
		       int job,
		       int x0, int y0,
		       int x1, int y1)
{

    switch (job) {
    case RASTER_INIT: 
	{
	    int num_elements;
	    seq_result **data;
	    seq_result *result;
	    int num_funcs;
	    out_raster *output;
	    RasterResult *raster_result;
	    int raster_id;
   
	    num_elements = seq_num_results();
	    if (num_elements == 0)
		return;
	    
	    data = (seq_result **)xmalloc(num_elements * sizeof(seq_result *));
	    if (-1 == search_reg_data(comparison2, (void **)data, &num_funcs)){
		xfree(data);
		return;	
	    }
	    if (num_funcs == 0) {
		xfree(data);
		return;
	    }
    
	    result = data[0];
	    output = result->output;

	    Tcl_VarEval(output->interp, "GetRasterId ", raster_win, NULL);
	    raster_id = atoi(Tcl_GetStringResult(output->interp));

	    if (NULL == (raster_result = raster_id_to_result(raster_id))) {
		xfree(data);
		return;
	    }
	    
	    remove_all_raster_cursors(output->interp, raster, raster_result);
	    xfree(data);

	    break;
	}
    case RASTER_REPLOT_ALL:
	SeqReplotResults(raster, raster_win, 1, 0, x0, y0, x1, y1);
	break;
    case RASTER_REPLOT_SLIVER:
	SeqReplotResults(raster, raster_win, 0, 0, x0, y0, x1, y1);
	break;
    case RASTER_REPLOT_ZOOM:
	/* NOTE: this is different in nip4 */
	SeqReplotResults(raster, raster_win, 0, 1, x0, y0, x1, y1);
	break;
    }
}

/*
 * Superimpose a result onto an existing window 
 */
void SeqSuperimposeResult(Tcl_Interp *interp,
			  char *raster_win,
			  int result_id,
			  double o_wx0,
			  double o_wy0,
			  double o_wx1,
			  double o_wy1)
{
    seq_result *result;
    out_raster *output;
    Tcl_CmdInfo cmd_info;
    Tk_Raster *raster;
    double wx0, wy0, wx1, wy1;
    double p2, q2;
    double m, c;
    d_line *dim;
    seq_reg_info info;

#ifdef DEBUG
    printf("SeqSuperimposeResult %d\n", result_id);
#endif

    result = seq_id_to_result(result_id);
    output = result->output;

    if (Tcl_GetCommandInfo(interp, raster_win, &cmd_info) == 0) 
	return;
    raster = (Tk_Raster*)cmd_info.clientData;
 
    /* 
     * get the current scroll region for the raster
     */
    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);
    
    /* find dimensions of result */
    info.job = SEQ_RESULT_INFO;
    info.op = DIMENSIONS;
    info.result = NULL;
    seq_result_notify(result_id, (seq_reg_data *)&info, 0);
    if (!info.result) {
	return;
    }
    dim = (d_line *)info.result;

    /* superimpose, update an exising raster */
    p2 = (wy1 - wy0) * (dim->y0 - o_wy0) / (o_wy1 - o_wy0) + wy0;
    q2 = (wy1 - wy0) * (dim->y1 - o_wy0) / (o_wy1 - o_wy0) + wy0;
    
    m = (p2 - q2) / (dim->y0 - dim->y1);
    c = p2 - (m * dim->y0);

    output->sf_c = (m * output->sf_c) + (c);
    output->sf_m *= m;

    /* enlarge x scroll region if necessary */
    RasterSetWorldScroll(raster, o_wx0, wy0, o_wx1, wy1);
}

d_point E_FindNearestLine(seq_result *result,
			  d_point start,
			  double x_scale)
{
    e_graph *e_data;
    int n_pts;
    e_obj *lines;
    int i;
    double min = DBL_MAX;
    d_point nearest_pt;
    double a, b, c;
    double x1, y1, x2, y2;
    double dist, dist1, dist2;
    double start_x;

    e_data = result->data; 
    lines = e_data->d_obj;
    n_pts = e_data->n_data_obj;


    for (i = 0; i < n_pts; i++) {
	x1 = lines[i].pos.x0/x_scale;
	y1 = lines[i].pos.y0;
	x2 = lines[i].pos.x1/x_scale;
	y2 = lines[i].pos.y1;

      start_x = start.x / x_scale;

      a = (y1 - y2) / (x2 - x1);
      
      b = 1.0;
      c = (-1 * a * x2) - y2;
      
#ifdef DEBUG
      printf("a %f b %f c %f x1y1 %f x2y2 %f x_scale %f\n", a, b, c, (a*x1)+(b*y1)+c, (a*x2)+(b*y2)+c, x_scale);
#endif

      /* 
       * check boundary conditions: if x,y of point is perpendicular to the
       * line then find the dist of perpendicular 
       * else find the dist of point to the nearest end of the line
       */
      if (start_x >= x1 && start_x <= x2 && start.y >= y1 && start.y <= y2) {
	dist = fabs(((a * start_x) + (b * start.y) + c) /
	  (sqrt((a * a) + (b * b))));

	if (dist < min) {
	  min = dist;
	  nearest_pt.x = lines[i].pos.x0;
	  nearest_pt.y = lines[i].pos.y0;
	}
      } else {
	dist1 = sqrt(((start_x - x1) * (start_x - x1)) + 
	  ((start.y - y1) * (start.y - y1)));
	dist2 = sqrt(((start_x - x2) * (start_x - x2)) + 
	  ((start.y - y2) * (start.y - y2)));
	
	if (dist1 < min) {
	  min = dist1;
	  nearest_pt.x = lines[i].pos.x0;
	  nearest_pt.y = lines[i].pos.y0;
	}
	if (dist2 < min) {
	  min = dist2;
	  nearest_pt.x = lines[i].pos.x0;
	  nearest_pt.y = lines[i].pos.y0;

	}
      }
    }
#ifdef DEBUG
    printf("START %d %f start_x %f END %d %f\n", start.x, start.y, start_x, nearest_pt.x,nearest_pt.y);
#endif
    return nearest_pt;
}


/*
 * returns the nearest match coords to position start
 */ 
d_point FindNearestLine(seq_result *result,
			d_point start,
			double x_scale)
{
    d_plot *data = result->data;
    int i;
    double min = DBL_MAX;
    d_point nearest_pt;
    int n_pts = data->n_pts;
    double a, b, c;
    double x1, y1, x2, y2;
    double dist, dist1, dist2;
    double start_x;

    /* emboss dot plot are always lines */
    if (result->graph == SEQ_E_DOT) {
	return(E_FindNearestLine(result, start, x_scale));
    }

    nearest_pt.x = 0;
    nearest_pt.y = 0;

    /* equation of line is ax + by + c = 0 */
    /* if (x1,y1) and (x2,y2) define the ends of the line, and set b = 1 then:
     * a = (y1 - y2) / (x2 - x1)
     * b = 1
     * c = (-1 * a * x2) - y2
     */
    for (i = 0; i < n_pts; i++) {
      x1 = data->p_array[i].x/x_scale;
      y1 = data->p_array[i].y;
      x2 = (data->p_array[i].x + (data->p_array[i].score - 1))/x_scale;
      y2 = data->p_array[i].y + (data->p_array[i].score - 1);

      start_x = start.x / x_scale;

      a = (y1 - y2) / (x2 - x1);
      
      b = 1.0;
      c = (-1 * a * x2) - y2;
      
#ifdef DEBUG
      printf("a %f b %f c %f x1y1 %f x2y2 %f x_scale %f\n", a, b, c, (a*x1)+(b*y1)+c, (a*x2)+(b*y2)+c, x_scale);
#endif

      /* 
       * check boundary conditions: if x,y of point is perpendicular to the
       * line then find the dist of perpendicular 
       * else find the dist of point to the nearest end of the line
       */
      if (start_x >= x1 && start_x <= x2 && start.y >= y1 && start.y <= y2) {
	dist = fabs(((a * start_x) + (b * start.y) + c) /
	  (sqrt((a * a) + (b * b))));

	if (dist < min) {
	  min = dist;
	  nearest_pt.x = data->p_array[i].x;
	  nearest_pt.y = data->p_array[i].y;
	}
      } else {
	dist1 = sqrt(((start_x - x1) * (start_x - x1)) + 
	  ((start.y - y1) * (start.y - y1)));
	dist2 = sqrt(((start_x - x2) * (start_x - x2)) + 
	  ((start.y - y2) * (start.y - y2)));
	
	if (dist1 < min) {
	  min = dist1;
	  nearest_pt.x = data->p_array[i].x;
	  nearest_pt.y = data->p_array[i].y;
	}
	if (dist2 < min) {
	  min = dist2;
	  nearest_pt.x = data->p_array[i].x;
	  nearest_pt.y = data->p_array[i].y;

	}
      }
    }
#ifdef DEBUG
    printf("START %d %f start_x %f END %d %f\n", start.x, start.y, start_x, nearest_pt.x,nearest_pt.y);
#endif
    return nearest_pt;
}

d_point FindNearestMatch(seq_result *result,
			 d_point start,
			 double x_scale)
{
    d_plot *data = result->data;
    int i;
    int x, y;
    double min = DBL_MAX;
    double square;
    d_point nearest_pt;
    int n_pts = data->n_pts;

    nearest_pt.x = 0;
    nearest_pt.y = 0;
    
    for (i = 0; i < n_pts; i++) {
	x = (start.x - data->p_array[i].x)/x_scale;
	y = (start.y - data->p_array[i].y);
#ifdef DEBUG
	printf("x %d y %f \n", data->p_array[i].x, data->p_array[i].y);
#endif
	square = (double)x * x + y * y;
#ifdef DEBUG
	printf("square %f min %f \n", square, min);
#endif
	if (square < min) {
	    min = square;
	    nearest_pt.x = data->p_array[i].x;
	    nearest_pt.y = data->p_array[i].y;
	}
    }
#ifdef DEBUG
    printf("START %d %f END %d %f\n", start.x, start.y, nearest_pt.x,nearest_pt.y);
#endif
    return nearest_pt;
}

/*
 * find the first sip result on raster_id and seq_id of direction.
 * used when double click on raster to bring up a sequence display
 */
int seq_find_result(char *raster_win,
		    int seq_id_h,
		    int seq_id_v)
{
    int num_elements;
    seq_result **data;
    seq_result *result;
    int num_funcs;
    out_raster *output;
    int i;

    num_elements = seq_num_results();
    if (num_elements == 0)
	return -1;

    data = (seq_result **)xmalloc(num_elements * sizeof(seq_result *));
 	
    if (-1 == search_reg_data(comparison2, (void **)data, &num_funcs)) {
	xfree(data);
	return -1;	
    }
    if (num_funcs == 0) {
	xfree(data);
	return -1;
    }
    /* 
     * go through all seq results and return the first one with the correct
     * raster and sequence
     */
    for (i = 0; i < num_funcs; i++) {
	result = data[i];
	output = result->output;

	if ((strcmp(output->raster_win, raster_win) == 0) && (seq_id_h == result->seq_id[HORIZONTAL]) && (seq_id_v == result->seq_id[VERTICAL])) {
	    xfree(data);
	    return result->id;
	}
    }
    xfree(data);
    return -1;
}

