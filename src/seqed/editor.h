#if !defined(EDITOR_H)
#define EDITOR_H

#include "read_sequence.h"
#include "seq_reg.h"

#define INSERT 1
#define DELETE 2
#define CHANGE 3
#define INSERT_FRAGMENT 4
#define DELETE_FRAGMENT 5
#define CHANGE_FRAGMENT 6
#define INSERT_FRAGMENT_COMPLEMENT 7
#define FILL 8
#define TRIM 9
#define ADD_FEATURE 10
#define REMOVE_FEATURE 11
#define ADD_SEQUENCE 12
#define REMOVE_SEQUENCE 13

#define TEXT   0 
#define GRAP   1

typedef struct EDIT_ {
    int seq_id;
    int position;
    int operation;
    int consensus;
    char *string;  
    SEQUENCE *sequence;
    struct EDIT_ *next;
}EDIT;

typedef struct EDITS_ {
  EDIT *head;
} EDITS;

enum member_type {GROUP, SEQ};

typedef union {
    struct SEQ_GROUP_ *group;
    SEQUENCE *sequence;
} DATA;

typedef struct GROUP_MEMBER_ {
    enum  member_type type;
    DATA  *data;
} GROUP_MEMBER;

typedef struct SEQ_GROUP_ {
    int   nmembers;
    int   group_id;
    int   group_parent_id;
    char  *group_name;
    GROUP_MEMBER **members;
} SEQ_GROUP;

/* coursor_e come from SPIN cursor_t, it added a new element posy */
typedef struct _cursor_e {
    int  editor_reg_id; /* editor registration ID */ 
    char *frame_name;   /* the current cursor's win name */
    int  id;		/* Cursor identification number */
    int  refs;		/* Number of concurrent uses for this cursor */
    int  private;	/* Whether this is a private cursor */
    int  abspos; 	/* Absolute position in sequence */
    int  posy;          /* vertical position in editor */ 
    int  job;    	/* move, reused, or delete */
    char *colour;       /* colour of cursor */
    int  line_width;    /* line width of cursor */
    int  direction;     /* direction of cursor */
    int  sent_by;
    struct _cursor_e *next;
} cursor_e;

typedef struct _cursors_ {
    cursor_e *cursor;
    struct _cursors_ *next;
} cursors;

typedef struct EDITOR_RECORD_ {
    SEQ_GROUP **seq_group;
    int       num_group;
    int       ft_imode; /* FT insertion mode; 0:extend; 1:break; 2:delete */
    int       text; /* n: numbers of opened */
    int       graphical; /* 1:opened */
    EDITS     *edits;
    cursors   *cursors;
    int       prev_pos;
} EDITOR_RECORD;

typedef struct SELECTION_ {
    int ed_id;
    int member_id;
    int sel_first;
    int sel_last;
    char *name_first;
    char *name_last;
} SELECTION;

typedef struct EDITOR_RECORDS_ { 
    EDITOR_RECORD **editor_record;
    SEQUENCE *clipboard;
    SELECTION *selection;/*FIXME: should in here or EDITOR_RECORD?*/
    int num_editor;
} EDITOR_RECORDS;

typedef struct {
    char *sel_seq;
} editor_seq_arg;

typedef struct {
    char *win_name;
    int  ed_id;
} editor_data;




typedef struct editor_cursor_ {
    int rid;  /*registration id*/
    int ed_id;
    int seq_id; /*fixme*/
    int cursorPos;
    int prev_pos;
    int cursor_visible;
    cursor_e *cursor;
} editor_cursor;


EDITOR_RECORDS *init_editor_records (void);
EDITOR_RECORD *init_editor_record (void);
GROUP_MEMBER *init_group_member (void);
SEQ_GROUP *init_seq_group (void);
EDITS *create_edits(void);
EDIT *create_edit(void);
int add_edit(EDITS *edits, EDIT *edit);
SEQUENCE *copy_sequence (SEQUENCE *f);
DATA *init_data_sequence (void);
int realloc_seq_group (SEQ_GROUP *seq_group, int num_members );
int realloc_editor_records (EDITOR_RECORDS *ers, int num_editor);
int extend_unique_name (SEQUENCE *s);
int create_copy_for_editor (Tcl_Interp *interp, int seq_id);
void free_editor_records (EDITOR_RECORDS *ers);
void free_seq_group (SEQ_GROUP *sg);
SELECTION *init_selection (void);
void free_selection (SELECTION *s);

/* for text and graphic editor */
int Editor_init(Tcl_Interp *interp);
int GetEdenLength (int seq_id);
int GetEdenType (int seq_id);
int GetEdenNums (void);
int GetEdenId(int seq_num);
int GetEdenNum(int id);
char *GetEdenSeq (int seq_id);
char *GetEdenName (int seq_id);
SEQUENCE *GetEdenSequence (int id);
void replace_sequence (int seq_num, SEQUENCE *s);

int delete_update_sequence (int ed_id, int group_id, int seq_id, int pos, SEQUENCE *sequence);
int insert_update_sequence (int ed_id, int group_id, int member_id, int pos, SEQUENCE *sequence);
int replace_update_sequence (int ed_id, int group_id, int member_id, int pos, SEQUENCE *s_replace, SEQUENCE *s_with);
int copy_update_buffer (int ed_id, int group_id, int member_id, int pos, SEQUENCE *s);
SEQUENCE *get_editor_buffer (void);
void set_editor_buffer (SEQUENCE *s);
int check_if_already_in_editor (int seq_id);
int add_new_editor (Tcl_Interp *interp, int seq_id);
int insert_fragment_in_editor ( EDITOR_RECORD *er, EDIT *edit, int add_to_undo );
int delete_fragment_in_editor ( EDITOR_RECORD *er, EDIT *edit, int add_to_undo );
FEATURE_TABLE *delete_fragment_in_sequence ( SEQUENCE *sequence, SEQUENCE *segment, int position);
int insert_fragment_in_sequence (SEQUENCE *sequence, SEQUENCE *segment, int position, int imode);
/* for text editor */
EDITOR_RECORD *GetEditor (int seq_id);
int GetEditorNum (void);
char *GetEditorSeq (EDITOR_RECORD *er, int group_id, int member_id);
char *GetEditorSeqName (EDITOR_RECORD *er, int group_id, int member_id);
int GetEditorSeqLength (EDITOR_RECORD *er, int group_id, int member_id);
FEATURE_TABLE *GetEditorSeqFt (EDITOR_RECORD *er, int group_id, int member_id);
int GetEdIdFromSeqId (int seq_id);

int GetMemberIdFromSeqId (int seq_id);
int GetEdIdBySeqNum (int seq_num, int ed_id);
int GetGroupIdFromSeqId (int seq_id);
void delete_editor (int ed_id);
EDIT *create_insert_fragment_edit (EDITOR_RECORD *er, int member_id, SEQUENCE *fi, int pos, int cons);
int update_consensus (SEQ_GROUP *sg, int job );
int GetEditorSeqId (int ed_id, int group_id, int member_id);
cursor_e *get_editor_cursor (int ed_id, int member_id, char *frame_name);
void set_editor_selection (SELECTION *sel);
void set_editor_cursor (int ed_id, cursor_e *cursor);
void convert_to_upper (char *string);
void convert_to_lower (char *string);

#endif
