#ifndef _EDITOR_REG_H_
#define _EDITOR_REG_H_

#include "seq_reg.h"
#include "editor.h"

#define SEQ_CHANGED          16  /* sequence has been changed via editors*/
#define END_EDITOR_QUIT      17  /* delete an end_editor window */
#define TEXT_EDITOR_QUIT     18  /* delete a text_editor window */
#define GRAPHIC_EDITOR_QUIT  19  /* delete a graphic_editor window */
#define SEQ_SELECTED         20  /* a selection has been made */
#define EDITOR_CURSOR_NOTIFY 21  /* the icursor has been moved */ 


/* added for editor */
typedef struct {
    int job;	/* SEQ_CHANGED, implies sequenced is changed */
    SELECTION *selection;
} seq_reg_changed;

typedef struct {
    int job;
    char *frame_name;
    SELECTION *selection;
} seq_reg_selected;

typedef struct {
    int job;
    cursor_e *cursor;
    SELECTION *selection;
} seq_reg_cursor_notify;

typedef struct {
    int job;
    char *frame_name;
} seq_reg_exit;

/* editor_reg_data from seq_reg_data, so may be change back? */ 
typedef union _editor_reg_data {
    /* MUST be first here and in job data structs */
    int job;
    seq_reg_generic	  generic;
    seq_reg_query_name 	  name;
    seq_reg_get_ops	  get_ops;
    seq_reg_invoke_op	  invoke_op;
    seq_reg_plot	  plot;
    seq_reg_info	  info;
    seq_reg_delete	  delete;
    seq_cursor_notify     cursor_notify;
    seq_reg_cursor_delete cursor_delete;
    seq_reg_cursor_move   cursor_move;
    seq_reg_sequence_type sequence_type;
    /* added for editor */
    seq_reg_changed       changed;
    seq_reg_selected      selected;
    seq_reg_exit          exit; 
    seq_reg_cursor_notify cursor_moved;
} editor_reg_data;

typedef struct {
    void (*func)(int seq_num, void *fdata, editor_reg_data *jdata);
    void *fdata;
    time_t time;
    int type;
    int id;
} edit_reg;

typedef struct text_editor_result_ {
    void (*op_func)(int seq_num,
		    void *obj,
		    editor_reg_data *data);
    Tcl_Interp *interp;
    char frame_name[1024];
    char win_name[1024];
    char *colour;
    int seq_id;
    int ed_id; 
    int index; /*register ID */
    int cursor_id;
} text_editor_result;

typedef struct editor_result_ {
    void (*op_func)(int seq_num, void *obj, editor_reg_data *data);
    void (*pr_func)(void *obj, seq_reg_plot *plot);
    void (*txt_func)(void *obj);
    void *data;         /* data to plot result */
    void *input;        /* input parameters */
    void *output;       /* general info about plotting the result */
    int id;             /* unique result identifier */
    int seq_id[2];      /* unique sequence identifier */
    int type;           /* unique name for each result type */
    int frame;          /* reading frame of a result or 0 if not applicable */
    void *text_data;    /* specific data required for textual output */
    int graph;          /* defines the data type */
} editor_result;

/*
 *----------------------------------------------------------------------------
 * Function prototypes
 *----------------------------------------------------------------------------
 */

/* array of functions associated with seq_num */
#define editor_func_array(seq_num) (arr(Array, editor_reg, (seq_num)))
#define editor_Nfuncs(seq_num)     (ArrayMax(editor_func_array(seq_num)))
#define editor_cursor(seq_num)     (arr(cursor_e *, editor_cursor_reg, (seq_num))) /*FIXME */

void ed_delete_cursor(int seq_num, int id, int private);
cursor_e *editor_create_cursor(int seq_num, int private, char *colour, 
			int line_width, int new_cursor, int direction);
/*
 * Initialise the sequence register lists
 *
 * Returns 0 on success and -1 for error;
 */
int editor_register_init (Tcl_Interp *interp);

/*
 * Add a new sequence to the registration
 */
int add_editor_reg(int index);

/*
 * Delete a sequence from the registration list
 */
void remove_sequence_from_registration (int seq_num);

/*
 * return a unique identifier for each set of results
 */
int get_editor_reg_id(void);

/*
 * Registers func(seq_num, fdata, jdata) with sequence 'seq_num'.
 * Doesn't check if the (func,fdata) pair are already existant.
 *
 * Returns 0 on success, and -1 for error. 
 */
int editor_register(int seq_num,
		 void (*func)(int seq_num, void *fdata, editor_reg_data *jdata),
		 void *fdata, int type, int id);

/*
 * Deregisters func(seq_num, data, jdata) from sequence 'seq_num'.
 *
 * Returns 0 for success, and -1 for error.
 */
int editor_deregister(int seq_num,
		   void (*func)(int seq_num, void *fdata, editor_reg_data *jdata),
		   void *fdata);
/*
 * Uses the register list for a given seq_num to call a particular job.
 */

void  editor_notify(int seq_num, editor_reg_data *jdata);

/*
 * Uses the register list for a given result to call a particular job.
 */
void editor_result_notify(int id, editor_reg_data *jdata, int all);

/*
 * Uses the register list for all results to call a particular job.
 */

cursor_e *ed_find_cursor(int *seq_num, int id, int direction);

text_editor_result *init_text_editor_result (void);
#endif
