#ifndef _SEQ_REG_H_
#define _SEQ_REG_H_

#include <time.h>

#include "tkRaster.h"
/* registration functions */
#define SEQ_QUERY_NAME    0 
#define SEQ_GET_OPS       1 
#define SEQ_INVOKE_OP     2
#define SEQ_PLOT          3
#define SEQ_RESULT_INFO   4
#define SEQ_HIDE          5
#define SEQ_DELETE        6   /* delete a sequence */
#define SEQ_QUIT          7   /* delete a result */
#define SEQ_REVEAL        8
#define SEQ_CURSOR_NOTIFY 9
#define SEQ_CURSOR_DELETE 10
#define SEQ_GENERIC       11
#define SEQ_KEY_NAME      12
#define SEQ_GET_BRIEF     13
#define SEQ_SEQUENCE_TYPE 14
#define SEQ_SETCURSOR     15

/* types */
#define SEQ_PLOT_PERM   0
#define SEQ_PLOT_TEMP   1
#define SEQ_RASTER      2
#define SEQ_SEQED       3
#define SEQ_SENDER      4

/*
 * Tasks defined in the reg_generic request. These apply to one specific type
 * of window.
 */
#define TASK_SEQED_SETCURSOR	0
#define TASK_SEQED_GETCURSOR	1

/* result information jobs */
#define WINDOW 0

/* sequence direction */
#define NEITHER -1
#define HORIZONTAL 0
#define VERTICAL 1

/*
 * The cursor contains information about the current display/editing cursors.
 * Cursors can be shared. A private cursor is private for one controlling
 * display, and can only be shared with other displays that do not need
 * access to their own private cursor.
 */
typedef struct _cursor_t {
    int id;		/* Cursor identification number */
    int refs;		/* Number of concurrent uses for this cursor */
    int private;	/* Whether this is a private cursor */
    int abspos; 	/* Absolute position in sequence */
    int job;    	/* move, reused, or delete */
    char *colour;       /* colour of cursor */
    int line_width;     /* line width of cursor */
    int direction;      /* direction of cursor */
    int sent_by;
    struct _cursor_t *next;
} cursor_t;

#define CURSOR_MOVE      (1<<0)
#define CURSOR_INCREMENT (1<<1)
#define CURSOR_DECREMENT (1<<2)
#define CURSOR_DELETE    (1<<3)
#define NUM_CURSORS    100

typedef struct {
    int job;	/* SEQ_GENERIC */
    int task;	/* Some specific task */
    void *data;	/* And data associated with the task */
    void *result; /* result of operation */
} seq_reg_generic;

typedef struct {
    int id;
    char *line;
    char *time;
} seq_reg_name;

typedef struct {
    int job;	/* REG_QUERY_NAME */
    char *line;	
} seq_reg_query_name, seq_reg_key_name, seq_reg_brief;

typedef struct {
    int job;	/* REG_DELETE, REG_QUIT REG_PLOT */
} seq_reg_delete, seq_reg_quit;

typedef struct {
    int job;
    int x0;
    int x1;
    int y0;
    int y1;
} seq_reg_plot;

typedef struct {/* Lists valid operations for the results window */
    int job;	/* REG_GET_OPS */
    char *ops;	/* Somewhere to place ops in, unalloced to start with */
} seq_reg_get_ops;

typedef struct {
    int job;	      /* REG_INVOKE_OP */
    int op;           /* Operation to perform */
} seq_reg_invoke_op;

typedef struct {
    int job;	   /* SEQ_INFO */
    int direction; /* direction of sequence */
    int op;	   /* Operation to perform */
    void *result;  /* result of operation */
} seq_reg_info;

typedef struct {
    int job;
    void *data;
} seq_reg_sequence_type;

typedef struct {
    int job;	      /* REG_CURSOR_NOTIFY */
    cursor_t *cursor; /* Cursor structure */
} seq_cursor_notify;

typedef struct {
    int job;
    int seqed_id;
    int pos;
    int env;
} seq_reg_cursor_delete;

typedef struct {
    int job;
    int seqed_id;
    int prev_pos;
    int env;
    int pos;
} seq_reg_cursor_move;

typedef union _seq_reg_data {
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
} seq_reg_data;


/*
 *----------------------------------------------------------------------------
 * Function prototypes
 *----------------------------------------------------------------------------
 */

/* array of functions associated with seq_num */
#define seq_func_array(seq_num)    (arr(Array, sequence_reg, (seq_num)))
#define seq_Nfuncs(seq_num)        (ArrayMax(seq_func_array(seq_num)))
#define seq_cursor(seq_num)        (arr(cursor_t *, seq_cursor_reg, (seq_num)))

typedef struct {
    void (*func)(int seq_num, void *fdata, seq_reg_data *jdata);
    void *fdata;
    time_t time;
    int type;
    int id;
} seq_reg;

typedef struct {
    int seqed_id;  /* Which sequence editor */
    int pos;	   /* Position in sequence */
    char *col;     /* colour of cursor */
} task_seqed_setcursor;

int get_num_cursors(void);
cursor_t *create_cursor(int seq_num, int private, char *colour, 
			int line_width, int new_cursor, int direction);
void delete_cursor(int seq_num, int id, int private);
int find_nearest_cursor(Tk_Raster *raster, int seq_num, int pos,
			int max_dist, int direction, int *cursor_pos);

/*
 * Initialise the sequence register lists
 *
 * Returns 0 on success and -1 for error;
 */
int seq_register_init(Tcl_Interp *interp);

/*
 * Add a new sequence to the registration
 */
int add_reg_seq(int index);

/*
 * Delete a sequence from the registration
 */
void delete_reg_seq(int index);

/*
 * return a unique identifier for each set of results
 */
int get_reg_id(void);

/*
 * Registers func(seq_num, fdata, jdata) with sequence 'seq_num'.
 * Doesn't check if the (func,fdata) pair are already existant.
 *
 * Returns 0 on success, and -1 for error. 
 */
int seq_register(int seq_num,
		 void (*func)(int seq_num, void *fdata, seq_reg_data *jdata),
		 void *fdata, int type, int id);

/*
 * Deregisters func(seq_num, data, jdata) from sequence 'seq_num'.
 *
 * Returns 0 for success, and -1 for error.
 */
int seq_deregister(int seq_num,
		   void (*func)(int seq_num, void *fdata, seq_reg_data *jdata),
		   void *fdata);

/*
 * Uses the register list for a given seq_num to call a particular job.
 */
void  seq_notify(int seq_num, 
		 seq_reg_data *jdata);

/*
 * Uses the register list for a given result to call a particular job.
 */
void seq_result_notify(int id, seq_reg_data *jdata, int all);

/*
 * Uses the register list for all results to call a particular job.
 */
void seq_result_notify_all(seq_reg_data *jdata);

void seq_type_notify(int seq_num, int type, seq_reg_data *jdata);

/*
 * return the type of result ie PERMenant or TEMPorary
 */
int seq_get_type(int id);

void *result_data(int id, int seq_num);

seq_reg_name *seq_result_names(int *num_elements);

/*
 * Returns a time (string) that a given id was registered. Assumes all
 * contig registrations for a particular id are registered together.
 * 'contig' isn't really needed here, but it's currently known by tcl
 * and speeds up our search.
 */
char *seq_result_time(int seq_num, int id);

/*
 * Dump lists
 */
void seq_register_dump(void);

/*
 * returns the number of sequences registered
 */
int seq_num_seqs(void);

/*
 * total number of results registered
 */
int seq_num_results(void);

/*
 * creates an array of data that succeed the comparison function
 */
int search_reg_data(int (*comparison)(void *fdata, int type),
		    void **array, int *num_elements);

int *result_to_seq_nums(int id, int *num_seqs);
int *result_to_seq_nums(int id, int *num_seqs);
seq_reg **result_to_regs(int id);
cursor_t *find_cursor(int *seq_num, int id, int direction);

#endif
