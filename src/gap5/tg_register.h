/*
 * File: tg_register.h: derived from gap4/io-reg.h
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *         Hills Road
 *         Cambridge CB2 2QH
 *         United Kingdom
 *
 * Description: Handles registration of gap data with display routines.
 *
 * Created: 9 March 1995
 */

#ifndef _IO_REG_H_
#define _IO_REG_H_

#include <time.h>
#include <tg_gio.h>

/*
 *----------------------------------------------------------------------------
 * Definitions
 *----------------------------------------------------------------------------
 */

/*
 * Known job numbers. Used in flags as well as requests.
 * They're not all requests to perform actions - most are used only in
 * notification of an action already performed.
 */
#define REG_GENERIC		(1<<0)
#define REG_NUMBER_CHANGE	(1<<1)
#define REG_JOIN_TO		(1<<2)
#define REG_ORDER		(1<<3)
#define REG_LENGTH		(1<<4)
#define REG_QUERY_NAME		(1<<5)
/* The contig has been deleted */
#define REG_DELETE		(1<<6)
#define REG_GET_LOCK		(1<<7)
#define REG_SET_LOCK		(1<<8)
#define REG_COMPLEMENT		(1<<9)
#define REG_PARAMS	       (1<<10)
/* Quit this display */
#define REG_QUIT	       (1<<11)
#define REG_CURSOR_NOTIFY      (1<<12)
#define REG_GET_OPS	       (1<<13)
#define REG_INVOKE_OP	       (1<<14)
#define REG_ANNO	       (1<<15)
#define REG_REGISTER	       (1<<16)
#define REG_DEREGISTER	       (1<<17)
#define REG_HIGHLIGHT_READ     (1<<18)
/* Buffering */
#define REG_BUFFER_START       (1<<19)
#define REG_BUFFER_END         (1<<20)
#define REG_NOTE	       (1<<21)
#define REG_CHILD_EDIT	       (1<<22)
#define REG_RENAME	       (1<<23)


/*
 * Useful flags (musn't coincide with the job numbers)
 */
#define REG_FLAG_INVIS	       (1<<29)
#define REG_FLAG_INACTIVE      (1<<30)


/*
 * Some useful compound definitions
 */
#define REG_REQUIRED	(REG_QUERY_NAME | REG_DELETE | REG_QUIT | REG_PARAMS)
#define REG_DATA_CHANGE (REG_JOIN_TO | REG_LENGTH | REG_COMPLEMENT | REG_RENAME)
#define REG_OPS		(REG_GET_OPS | REG_INVOKE_OP)
#define REG_LOCKS	(REG_GET_LOCK | REG_SET_LOCK)
#define REG_REGISTERS	(REG_REGISTER | REG_DEREGISTER)
#define REG_BUFFER	(REG_BUFFER_START | REG_BUFFER_END)
#define REG_ALL 	(REG_REQUIRED | REG_DATA_CHANGE | REG_OPS | REG_LOCKS \
			 | REG_ORDER | REG_CURSOR_NOTIFY | REG_NUMBER_CHANGE \
			 | REG_ANNO | REG_REGISTERS | REG_HIGHLIGHT_READ \
			 | REG_BUFFER | REG_NOTE)


/*
 * Registered types. Each function that maintains registered data has
 * a unique type. This is to facilitate inter-data communication.
 */
#define REG_TYPE_UNKNOWN	    0
#define REG_TYPE_EDITOR		    1
#define REG_TYPE_FIJ		    2
#define REG_TYPE_READPAIR	    3
#define REG_TYPE_REPEAT		    4
#define REG_TYPE_QUALITY	    5
#define REG_TYPE_TEMPLATE	    6 
#define REG_TYPE_RESTRICTION	    7
#define REG_TYPE_STOPCODON	    8
#define REG_TYPE_CONTIGSEL	    9
#define REG_TYPE_CHECKASS          10
#define REG_TYPE_OLIGO             11
#define REG_TYPE_CONSISTENCY_DISP  12
#define REG_TYPE_CONFIDENCE        13
#define REG_TYPE_READING_COVERAGE  14
#define REG_TYPE_READPAIR_COVERAGE 15
#define REG_TYPE_STRAND_COVERAGE   16

/*
 * Tasks defined in the reg_generic request. These apply to one specific type
 * of window.
 */
#define TASK_GENERIC            0
#define TASK_EDITOR_SETCURSOR	1
#define TASK_EDITOR_GETCON	2

/*
 *----------------------------------------------------------------------------
 * Structures
 *----------------------------------------------------------------------------
 */

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
    tg_rec seq;		/* Cursor sequence number (0 for consensus) */
    int pos;		/* Position from start of sequence */
    int abspos; 	/* Absolute position in contig */
    int sent_by;	/* reg id of display that created cursor */
    int job;    	/* move, reused, or delete */
    struct _cursor_t *next;
} cursor_t;

#define CURSOR_MOVE      (1<<0)
#define CURSOR_INCREMENT (1<<1)
#define CURSOR_DECREMENT (1<<2)
#define CURSOR_DELETE    (1<<3)

/*
 * REG_GENERIC is used for requests specific to a single registered item.
 * As such it'll typically be used in a result_notify or type_contig_notify
 * request.
 */
typedef struct {
    int job;	/* REG_GENERIC */
    int task;	/* Some specific task */
    void *data;	/* And data associated with the task */
} reg_generic;

typedef struct {
    int job;	/* REG_NUMBER_CHANGE */
    int number;	/* New number */
} reg_number;
    
typedef struct {
    int job;	   /* REG_JOIN_TO */
    tg_rec contig; /* New contig number */
    int offset;	   /* Offset of old contig into new contig */
} reg_join;

typedef struct {
    int job;	/* REG_ORDER */
    int pos;	/* New order */
} reg_order;
    
typedef struct {
    int job;	/* REG_LENGTH, implies data change too */
    int length;	/* New length */
} reg_length;
    
typedef struct {
    int job;	/* REG_QUERY_NAME */
    char *line;	/* char[80] */
} reg_query_name;

typedef struct {
    int job;	/* REG_RENAME */
    char *name;	/* char[1024] */
} reg_rename;

typedef struct {
    int job;	/* REG_DELETE, REG_COMPLEMENT, etc */
} reg_delete, reg_complement, reg_anno, reg_buffer_start, reg_buffer_end;

#define REG_NOTE_CREATE	0
#define REG_NOTE_DELETE 1
#define REG_NOTE_EDIT   2
typedef struct {
    int job;	/* REG_NOTE */
    int note;	/* Note number */
    int task;	/* REG_NOTE_{CREATE,DELETE,EDIT} */
} reg_note;

#define REG_LOCK_READ	1
#define REG_LOCK_WRITE	2
typedef struct {
    int job;	/* REG_GET_LOCK */
    int lock;	/* Sends lock requirements, returns locks allowed */
} reg_get_lock, reg_set_lock, reg_quit;

typedef struct {/* Lists valid operations for the results window */
    int job;	/* REG_GET_OPS */
    char *ops;	/* Somewhere to place ops in, unalloced to start with */
} reg_get_ops;

typedef struct {
    int job;	/* REG_INVOKE_OP */
    int op;	/* Operation to perform */
} reg_invoke_op;

typedef struct {
    int job;	  /* REG_PARAMS */
    char *string; /* Private string (ie do not free) containing params */
} reg_params;

typedef struct {
    int job;	      /* REG_CURSOR_NOTIFY */
    cursor_t *cursor; /* Cursor structure */
} reg_cursor_notify;

typedef struct {
    int job;	   /* REG_REGISTER, REG_DEREGISTER */
    int id;	   /* Registration id */
    int type;	   /* Registration type */
    tg_rec contig; /* Contig number */
} reg_register, reg_deregister;

typedef struct {
    int job;	   /* REG_HIGHLIGHT_READ */
    tg_rec seq;	   /* Gel reading number (-ve == contig consensus) */
    int val;	   /* 1==highlight, 0==dehighlight */
} reg_highlight_read;

typedef union _reg_data {
    /* MUST be first here and in job data structs */
    int job;
    
    reg_generic		generic;
    reg_number		number;
    reg_join		join;
    reg_order		order;
    reg_length		length;
    reg_query_name	name;
    reg_delete		delete;
    reg_complement	complement;
    reg_get_lock	glock;
    reg_set_lock	slock;
    reg_quit		quit;
    reg_get_ops		get_ops;
    reg_invoke_op	invoke_op;

    reg_params		params;
    reg_cursor_notify	cursor_notify;
    reg_anno		annotations;
    reg_register	c_register;
    reg_deregister	c_deregister;
    reg_highlight_read	highlight;
    reg_buffer_start	buffer_start;
    reg_buffer_end	buffer_end;
    reg_note		note;
    reg_rename		rename;
} reg_data;


typedef struct {
    void (*func)(GapIO *io, tg_rec contig, void *fdata, reg_data *jdata);
    void *fdata;
    int id;
    time_t time;
    int flags;
    int type;
    int uid; /* A _unique_ identifier for this contig_reg_t */
    int ref_count; /* Auto destroy when this hits zero */
    HacheItem *hi[2]; /* Possible entries to hache for this reg */
} contig_reg_t;


/*
 * Task structures to be used with the above TASK_*_* definitions.
 */
typedef struct {
    int position;
    int seq;
} task_editor_setcursor;

typedef struct {
    char *con; /* Alloced by the contig editor */
    int lreg; /* Set lreg and rreg to 0 for all consensus */
    int rreg;
    int con_cut;
    int qual_cut;
} task_editor_getcon;


/* Used by result_names() */
typedef struct {
    char name[80];
    int id;
    tg_rec contig;
    contig_reg_t *r;
} result_name_t;

/*
 *----------------------------------------------------------------------------
 * Function prototypes
 *----------------------------------------------------------------------------
 */

#define io_contig_reg(io)    ((io)->contig_reg)
#define io_cursor_reg(io)    ((io)->contig_cursor)


/*
 * Initialise the contig register lists
 *
 * Returns 0 on success and -1 for error;
 */
int contig_register_init(GapIO *io);


/*
 * Deallocates memory used by the contig registration scheme.
 */
void contig_register_destroy(GapIO *io);


/*
 * Allocates new register id numbers as requires. No protection for wrap-
 * around, but this ought to take ages to ever happen.
 */
int register_id(void);


/*
 * Registers func(io, contig, fdata, jdata) with contig 'contig'.
 * Doesn't check if the (func,fdata) pair are already existant.
 *
 * Returns 0 on success, and -1 for error. 
 */
int contig_register(GapIO *io, tg_rec contig,
		    void (*func)(GapIO *io, tg_rec contig, void *fdata,
				 reg_data *jdata),
		    void *fdata,
		    int id, int flags, int type);


/*
 * Deregisters func(io, contig, data, jdata) from contig 'contig'.
 *
 * Returns 0 for success, and -1 for error.
 */
int contig_deregister(GapIO *io, tg_rec contig,
		      void (*func)(GapIO *io, tg_rec contig, void *fdata,
				   reg_data *jdata),
		      void *fdata);


/*
 * Uses the register list for a given contig to call a particular job.
 */
void contig_notify(GapIO *io, tg_rec contig, reg_data *jdata);


/*
 * Joins two registers lists. This doesn't check for duplicate entries.
 *
 * Returns 0 for success, and -1 for error.
 */
int contig_register_join(GapIO *io, tg_rec cfrom, tg_rec cto);


/*
 * Event handling for destroying a contig.
 * Call this just after we delete it. It sends out notification events and
 * updated the contig_reg hash tables.
 */
void contig_register_delete(GapIO *io, tg_rec contig);

/*
 * Iterates through all contigs_reg_t registered with 'contig' and having
 * a id of 'id'. start_from maybe NULL, but if not then it's used to
 * store the last hache table position to facilitate iteration.
 *
 * An 'id' of zero implies no filtering on ID is needed - return all.
 *
 * Returns contig_reg_t* on success
 *         NULL on failure
 */
contig_reg_t *get_reg_by_contig_id(GapIO *io, tg_rec contig, int id,
				   HacheItem **start_from);
contig_reg_t *get_reg_by_id(GapIO *io, int id,
			    HacheItem **start_from);

/*
 * Converts a registration id to an array of contig_reg_t structures,
 * suitable for further interrogation.
 *
 * Returns NULL for failure or a NULL terminated 'contig_reg_t *' array which
 * is expected to be xfree()d by the calling function.
 */
contig_reg_t **result_to_regs(GapIO *io, int id);


/*
 * Generates description of registered functions.
 */
result_name_t *result_names(GapIO *io, int *nres);

/*
 * Returns a time (string) that a given id was registered. Assumes all
 * contig registrations for a particular id are registered together.
 * 'contig' isn't really needed here, but it's currently known by tcl
 * and speeds up our search.
 */
char *result_time(GapIO *io, int id);

/*
 * Uses the register list for a given result to call a particular job.
 * As per contig_notify except on a result basis rather than contig basis.
 * 'all' declares whether to send message to all registrations with this id
 * or just the first found.
 */
void result_notify(GapIO *io, int id, reg_data *jdata, int all);

/*
 * Returns the data component of a contig_reg_t for a specific id.
 *
 * We return only the first data for id found, or NULL if none found.
 */
void *result_data(GapIO *io, int id);

/*
 * Returns the id for a registered item of a particular type, If contig is
 * non zero, we search only this contig. Otherwise we scan all.
 *
 * We return only the first data for id found, or 0 if none found.
 */
int type_to_result(GapIO *io, int type, tg_rec contig);

/*
 * Notifies all (or the first) registered items of a given type.
 * Works across all contigs.
 */
int type_notify(GapIO *io, int type, reg_data *jdata);

/*
 * Notifies all (or the first) registered items of a given type within
 * a specified contig.
 */
int type_contig_notify(GapIO *io, tg_rec contig, int type,
		       reg_data *jdata, int all);

/*
 * Attempts to lock a contig for exclusive write access.
 * No record of the lock is kept; that's done implicitly by registration.
 *
 * Returns 0 for success and -1 for failure.
 */
int contig_lock_write(GapIO *io, tg_rec contig);

/*
 * Create a cursor for this contig.
 * If private == 1 then always create a new cursor.
 * Otherwise use an existing one if available
 * Returns the cursor pointer, or NULL for failure.
 */
cursor_t *create_contig_cursor(GapIO *io, tg_rec contig, int private, int sent_by);

/*
 * Given a cursor identifier, return the cursor structure or NULL if not
 * found.
 */
cursor_t *find_contig_cursor(GapIO *io, tg_rec contig, int id);

/*
 * Deletes a contig cursor for this option. If the cursor is in use
 * more than once then this simply decrements the reference count.
 * 'private' indicates whether this cursor was "your private one".
 */
void delete_contig_cursor(GapIO *io, tg_rec contig, int id, int private);

#endif
