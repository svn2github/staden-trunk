#if !defined(READ_FILE_H)
#define READ_FILE_H

#include "parse_feature.h"

#define MAX_SEQ_LINE 101 

typedef struct SITE_ {
    int pos1;
    int pos2;
    char *name;
    char *rec_seq;
    int direction; /* ie forward or reverse of recognition sequence */
} SITE;

typedef struct _SITES {
  SITE **site;
  int capacity;
  int used;
} SITES;

typedef struct FEATURE_TABLE_ {
  ft_entry **entry;
  int num_entry;
} FEATURE_TABLE;

typedef struct SITE_HANG_ {
    int  pos;    /* position to left of an internal cut */
    int  len;    /* +2 = overhang, -2 = underhang */
    char *hang;  /* a string which is overhang or underhang */
} SITE_HANG;


typedef struct _SEQUENCE {
    char *seq;
    char *name;
    char *source;
    int  id;         /* do we need 2 ids: one for SEQUENCES, one for members of SEQUENCES ? */
    int  parent_id;  /* id of sequence of which it is a fragment */
    int  offset;     /* offset in sequence of which it is a fragment */
    int  length;
    int  start;      /* active region */
    int  end;        /* active region */
    int  type;       /* DNA , RNA, PROTEIN */
    int  direction;  /* orientation relative to original */
    int  genetic_code;
    int  group_id;
    int  group_parent_id;
    SITE_HANG *left_end;   /* the hang infor when a fragment has been made */  
    SITE_HANG *right_end;  
    SITE_HANG *cut_site_1; /* the hang infor when doing insertion and replacement */
    SITE_HANG *cut_site_2;    
    FEATURE_TABLE *feature_table;
    SITES *sites;    /* current list of restriction sites */
    struct _SEQUENCE *fragment;
} SEQUENCE;

/* for type 0 = not set, 1 = DNA linear, 2 = DNA circular
   3 = RNA linear, 4 = RNA circular, 5 = protein. 
*/


typedef struct _SEQUENCES { 
    SEQUENCE **sequence;
    int num_seq;
} SEQUENCES;

SEQUENCE *init_sequence (void);
void free_sequence (SEQUENCE *sequence);
SEQUENCES *init_sequences (void);
void free_sequences (SEQUENCES *sequences);
char *GetSequenceSeq (int seq_id);
int GetSequenceLength (int seq_id);
int GetSequenceType (int seq_id);
int GetSequenceIdByName (char *name);
char *GetSequenceName (int seq_id);
int GetSequenceId(int num);
int GetSequenceNum(int seq_id);
int GetSequenceNums (void);
SEQUENCE *GetSequencesSequence (int id);
int make_copy_for_editor (int seq_id);
int Read_sequence_Init(Tcl_Interp *interp);
SITE_HANG *init_site_hang (void);
void free_site_hang (SITE_HANG *sh);
int add_sequence_to_sequences (SEQUENCE *s);
void remove_sequence_from_sequence_list (int seq_id);
FEATURE_TABLE *GetSequenceFt (int seq_id);
int save_change_to_sequence (int seq_id, SEQUENCE *s);
FEATURE_TABLE *read_embl_format_seq ( FILE *fp, char *line, char **seq, int *seq_len);
SEQUENCES *realloc_sequences ( SEQUENCES *ss, int num_seq);
int ReadFileInfo (ClientData cData, Tcl_Interp *interp, int argc,  char **argv);
int sequence_save (SEQUENCE *s);

extern SEQUENCES *sequences;

#endif
 
