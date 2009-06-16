/*
 * TO DO
 *
 * - Add io registration callbacks for NOTE_CREATE, NOTE_DELETE, etc.
 *
 * - more checks - note hand holding, note prev type checks.
 */

#include <string.h>
#include <stdlib.h>
#include <tk.h>
#include <time.h>
#include <ctype.h>

#include "IO.h"
#include "cli_arg.h"
#include "tagUtils.h"
#include "io-reg.h"
#include "gap_globals.h"
#include "text_output.h"
#include "gap_cli_arg.h"
#include "notes.h"
#include "io_utils.h"

/* Tcl interface to new_note() */
typedef struct {
    GapIO *io;
    char *type;
    char *to;
    int num;
} new_note_args;

int tcl_new_note(ClientData clientData, Tcl_Interp *interp,
		 int objc, Tcl_Obj *CONST objv[])
{
    new_note_args args;
    GCardinal itype;
    int ito;

    cli_args a[] = {
	{"-io",		ARG_IO,	 1, NULL,   offsetof(new_note_args, io)},
	{"-type",	ARG_STR, 1, "COMM", offsetof(new_note_args, type)},
	{"-to",		ARG_STR, 1, NULL,   offsetof(new_note_args, to)},
	{"-number",	ARG_INT, 1, "1",    offsetof(new_note_args, num)},
	{NULL,		0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    itype = str2type(args.type);
    if (strcmp(args.to, "database") == 0) {
	ito = GT_Database;
    } else if (strcmp(args.to, "reading") == 0) {
	ito = GT_Readings;
    } else if (strcmp(args.to, "contig") == 0) {
	ito = GT_Contigs;
    } else {
	return TCL_ERROR;
    }

    vTcl_SetResult(interp, "%d", new_note(args.io, itype, ito, args.num));

    flush2t(args.io);

    return TCL_OK;
}

/* Tcl interface to delete_note() */
typedef struct {
    GapIO *io;
    int num;
} del_note_args;

int tcl_delete_note(ClientData clientData, Tcl_Interp *interp,
		    int objc, Tcl_Obj *CONST objv[])
{
    del_note_args args;

    cli_args a[] = {
	{"-io",		ARG_IO,	 1, NULL,   offsetof(del_note_args, io)},
	{"-note",	ARG_INT, 1, NULL,   offsetof(del_note_args, num)},
	{NULL,		0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    vTcl_SetResult(interp, "%d", delete_note(args.io, args.num));

    flush2t(args.io);

    return TCL_OK;
}

/* Tcl interface to delete_note() */
typedef struct {
    GapIO *io;
    int num;
    char *type;
    char *comment;
} edit_note_args;

int tcl_edit_note(ClientData clientData, Tcl_Interp *interp,
		  int objc, Tcl_Obj *CONST objv[])
{
    edit_note_args args;

    cli_args a[] = {
	{"-io",		ARG_IO,	 1, NULL,   offsetof(edit_note_args, io)},
	{"-note",	ARG_INT, 1, NULL,   offsetof(edit_note_args, num)},
	{"-type",	ARG_STR, 1, "",     offsetof(edit_note_args, type)},
	{"-comment",	ARG_STR, 1, NULL,   offsetof(edit_note_args, comment)},
	{NULL,		0,	 0, NULL, 0}
    };

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    vTcl_SetResult(interp, "%d", edit_note(args.io, args.num, args.type,
					   args.comment));

    flush2t(args.io);

    return TCL_OK;
}

/*
 * ---------------------------------------------------------------------------
 * Note manipulation functions - add, delete, etc.
 * ---------------------------------------------------------------------------
 */

/*
 * Strips a note off the free_note list or allocates a new note.
 * Returns the note number.
 */
int get_free_note(GapIO *io) {
    GNotes n;
    GCardinal nn;

    /* Get note number 'nn' */
    if (io->db.free_notes) {
	nn = io->db.free_notes;
	note_read(io, nn, n);
	io->db.free_notes = n.next;
	DBDelayWrite(io);
	if (io->db.free_notes) {
	    note_read(io, io->db.free_notes, n);
	    n.prev = 0;
	    n.prev_type = 0;
	    note_write(io, io->db.free_notes, n);
	}
    } else {
	io_init_note(io, Nnotes(io)+1);
	nn = Nnotes(io);
    }

    return nn;
}


/*
 * Adds a new note to a reading, database, or contig list.
 *
 * ntype is the note type as in integer (eg from COMM, EXPT, etc)
 * gtype is the type to add a note to - GT_Database, GT_Contigs, GT_Readings
 * num   is the number of record to add to (contig number or reading number -
 *          ignored for GT_Database type).
 *
 * Returns the new note number or -1 for failure.
 */
int new_note(GapIO *io, int ntype, int gtype, int num) {
    int nn, nlist;
    GContigs c;
    GReadings r;
    GNotes n, nnext;
    time_t tt;
    reg_note rn;

    /* Get a free note */
    nn = get_free_note(io);

    /* Find the current note list */
    switch (gtype) {
    case GT_Database:
	nlist = io->db.notes;
	break;

    case GT_Contigs:
	contig_read(io, num, c);
	nlist = c.notes;
	break;

    case GT_Readings:
	gel_read(io, num, r);
	nlist = r.notes;
	break;

    default:
	return -1;
    }

    /* Initialise new note and link to nlist */
    note_read(io, nn, n);
    time(&tt);
    n.ctime_top = 0;
    n.ctime = (GCardinal)tt;
    n.mtime_top = 0;
    n.mtime = (GCardinal)tt;
    n.type = ntype;
    n.prev = (gtype == GT_Database) ? 0 : num;
    n.next = nlist;
    n.prev_type = gtype;
    n.annotation = 0;
    note_write(io, nn, n);

    /* Link nlist back to new note */
    if (nlist) {
	note_read(io, nlist, nnext);
	nnext.prev = nn;
	nnext.prev_type = GT_Notes;
	note_write(io, nlist, nnext);
    }

    /* Update database/contig/reading note list */
    switch (gtype) {
    case GT_Database:
	io->db.notes = nn;
	DBDelayWrite(io);
	break;

    case GT_Contigs:
	c.notes = nn;
	contig_write(io, num, c);
	break;

    case GT_Readings:
	r.notes = nn;
	gel_write(io, num, r);
	break;
    }

    /* Notify */
    rn.job = REG_NOTE;
    rn.note = nn;
    rn.task = REG_NOTE_CREATE;
    contig_notify(io, 0, (reg_data *)&rn);

    return nn;
}


/*
 * Deletes a note from a note list.
 *
 * nnum  is the note number.
 *
 * Returns 0 for success, -1 for failure.
 */
int delete_note(GapIO *io, int nnum) {
    GNotes n, nnext, nprev;
    GContigs c;
    GReadings r;
    reg_note rn;

    note_read(io, nnum, n);

    /* Link next note to our prev */
    if (n.next) {
	note_read(io, n.next, nnext);
	nnext.prev = n.prev;
	nnext.prev_type = n.prev_type;
	note_write(io, n.next, nnext);
    }

    /* Link prev note to our next */
    switch (n.prev_type) {
    case GT_Notes:
	note_read(io, n.prev, nprev);
	nprev.next = n.next;
	note_write(io, n.prev, nprev);
	break;

    case GT_Database:
	io->db.notes = n.next;
	DBDelayWrite(io);
	break;

    case GT_Contigs:
	contig_read(io, n.prev, c);
	c.notes = n.next;
	contig_write(io, n.prev, c);
	break;

    case GT_Readings:
	gel_read(io, n.prev, r);
	r.notes = n.next;
	gel_write(io, n.prev, r);
	break;
    }

    /* Deallocate note text */
    if (n.annotation) {
	deallocate(io, n.annotation);
	n.annotation = 0;
    }

    /* Link deleted note to free list */
    n.prev = 0;
    n.prev_type = 0; /* on free list */
    n.next = io->db.free_notes;
    note_write(io, nnum, n);
    io->db.free_notes = nnum;
    DBDelayWrite(io);

    /*
     * Link freelist back - not essential, but it helps to spot possible
     * database inconsistencies.
     */
    if (n.next) {
	note_read(io, n.next, nnext);
	nnext.prev = nnum;
	nnext.prev_type = GT_Notes;
	note_write(io, n.next, nnext);
    }

    /* Notify */
    rn.job = REG_NOTE;
    rn.note = nnum;
    rn.task = REG_NOTE_DELETE;
    contig_notify(io, 0, (reg_data *)&rn);

    return 0;
}

/*
 * Edits a note by chaning its type and/or comment. This also updates the
 * modification date, which is a mandatory edit.
 *
 * nnum is the note number.
 * type is the new type, or ""/NULL for no change.
 * comment is the new comment, or "" for blank comm.
 *
 * Returns 0 for success, -1 for failure.
 */
int edit_note(GapIO *io, int nnum, char *type, char *comment) {
    int itype;
    GNotes n;
    time_t tt;
    reg_note rn;

    note_read(io, nnum, n);

    if (type && *type != 0) {
	itype = str2type(type);
	n.type = itype;
    }

    if (comment) {
	if (*comment &&
	    strcmp(comment, " -- No text attached to this note --\n") != 0) {
	    if (!n.annotation) {
		n.annotation = allocate(io, GT_Text);
	    }
	    TextWrite(io, n.annotation, comment, strlen(comment));
	} else {
	    if (n.annotation) {
		deallocate(io, n.annotation);
		n.annotation = 0;
	    }
	}
    }

    time(&tt);
    n.mtime = (GCardinal)tt;
    note_write(io, nnum, n);

    /* Notify */
    rn.job = REG_NOTE;
    rn.note = nnum;
    rn.task = REG_NOTE_EDIT;
    contig_notify(io, 0, (reg_data *)&rn);

    return 0;
}

/*
 * Edits a note by chaning its time stamps. Specify the time stamps to be
 * zero for no change.
 *
 * Returns 0 for success, -1 for failure.
 */
int set_note_time(GapIO *io, int nnum, time_t ctime, time_t mtime) {
    GNotes n;
    reg_note rn;

    /* Update */
    note_read(io, nnum, n);
    if (ctime) {
	n.ctime_top = 0;
	n.ctime = (GCardinal)ctime;
    }
    if (mtime) {
	n.mtime_top = 0;
	n.mtime = (GCardinal)mtime;
    }
    note_write(io, nnum, n);

    /* Notify */
    rn.job = REG_NOTE;
    rn.note = nnum;
    rn.task = REG_NOTE_EDIT;
    contig_notify(io, 0, (reg_data *)&rn);

    return 0;
}

/*
 * Deletes an entire note list - eg when destroying a contig or a reading.
 *
 * Returns 0 for success, -1 for failure.
 */
int delete_note_list(GapIO *io, int nnum) {
    GNotes n;
    GContigs c;
    GReadings r;
    int last;

    if (0 == nnum)
	return 0;
    note_read(io, nnum, n);

    /* Disconnect from the predecessor */
    switch (n.prev_type) {
    case GT_Database:
	io->db.notes = 0;
	/* write done later */
	break;

    case GT_Contigs:
	contig_read(io, n.prev, c);
	c.notes = 0;
	contig_write(io, n.prev, c);
	break;

    case GT_Readings:
	gel_read(io, n.prev, r);
	r.notes = 0;
	gel_write(io, n.prev, r);
	break;

    case GT_Notes:
	break;

    default:
	return -1;
    }

    /* Set nnum's prev/prev_type to be 'free list' */
    n.prev_type = 0; /* Free list */
    n.prev = 0;
    note_write(io, nnum, n);

    /* Find the end of this note list, deallocating note text as we go */
    last = nnum;
    do {
	if (n.annotation) {
	    deallocate(io, n.annotation);
	    n.annotation = 0;
	    note_write(io, last, n);
	}
	if (n.next) {
	    last = n.next;
	    note_read(io, last, n);
	}
    } while (n.next);

    /* Link the end of this to the free notes list & vice versa. */
    n.next = io->db.free_notes;
    note_write(io, last, n);
    io->db.free_notes = nnum;
    DBDelayWrite(io);
    nnum = n.next;
    if (nnum) {
	note_read(io, nnum, n);
	n.prev = last;
	n.prev_type = GT_Notes;
	note_write(io, nnum, n);
    }

    return 0;
}

/*
 * Duplicate an entire note list. This is used when breaking a contig.
 *
 * Note that we assume cto has no existing note list - this is a valid
 * assumption when called from break contig.
 *
 * Returns 0 for success, -1 for failure.
 */
int dup_contig_notes(GapIO *io, int cfrom, int cto) {
    int first, new_first;
    int n_old, n_new;
    int next = 0;
    GContigs c;
    GNotes n1, n2;
    int prev_note, prev_type;

    /* Find first note to copy from */
    contig_read(io, cfrom, c);
    first = c.notes;

    if (!first)
	return 0;

    /* Create a new note list, starting from first (copied as new_first) */
    n_new = new_first = get_free_note(io);
    n_old = first;
    prev_note = cto;
    prev_type = GT_Contigs;

    do {
	/* Copy n_old to n_new */
	note_read(io, n_old, n1);
	memcpy(&n2, &n1, sizeof(n1));
	if (n1.annotation) {
	    char *tmp;
	    n2.annotation = allocate(io, GT_Text);
	    tmp = TextAllocRead(io, n1.annotation);
	    TextWrite(io, n2.annotation, tmp, strlen(tmp));
	    xfree(tmp);
	}
	n2.prev = prev_note;
	n2.prev_type = prev_type;
	if (n1.next) {
	    next = get_free_note(io);
	    n2.next = next;
	}
	note_write(io, n_new, n2);
	prev_note = n_new;
	prev_type = GT_Notes;
	n_new = next;
	n_old = n1.next;
    } while (n_old);

    /* Now update the 'cto' contig to link to our first note */
    contig_read(io, cto, c);
    c.notes = new_first;
    contig_write(io, cto, c);

    return 0;
}

/*
 * Fortran interface to dup_contig_notes().
 */
f_proc_ret dupnot_(f_int *HANDLE, f_int *CONT1, f_int *CONT2)
{
    GapIO *io;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();
    dup_contig_notes(io, *CONT1, *CONT2);

    f_proc_return();
}


/*
 * Merges the notes from two contigs.
 *
 * Returns 0 for success, -1 for failure.
 */
int merge_contig_notes(GapIO *io, int cfrom, int cto) {
    GContigs cf, ct;
    GNotes n;
    int last;

    contig_read(io, cfrom, cf);
    if (cf.notes == 0)
	return 0;
    contig_read(io, cto, ct);

    /* Find end of 'cto' notes */
    if (ct.notes) {
	last = ct.notes;
	do {
	    note_read(io, last, n);
	    if (n.next)
		last = n.next;
	} while (n.next);
    } else {
	last = 0;
    }

    /* Copy from 'cfrom' to 'cto' */
    if (last) {
	n.next = cf.notes;
	note_write(io, last, n);
	note_read(io, cf.notes, n);
	n.prev = last;
	n.prev_type = GT_Notes;
	note_write(io, cf.notes, n);
    } else {
	note_read(io, cf.notes, n);
	n.prev = cto;
	n.prev_type = GT_Contigs;
	note_write(io, cf.notes, n);
	ct.notes = cf.notes;
	contig_write(io, cto, ct);
    }

    /* Unlink note list from 'cfrom' */
    cf.notes = 0;
    contig_write(io, cfrom, cf);

    return 0;
}


/*
 * Fortran interface to merge_contig_notes().
 */
f_proc_ret mrgnot_(f_int *HANDLE, f_int *CONT1, f_int *CONT2)
{
    GapIO *io;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();
    merge_contig_notes(io, *CONT1, *CONT2);

    f_proc_return();
}


/*
 * ---------------------------------------------------------------------------
 * Special note types, and the actions to perform on them.
 * ---------------------------------------------------------------------------
 */

void execute_database_notes(GapIO *io, char *type) {
    GNotes n;
    int note;
    int itype = str2type(type);

    if (!exec_notes)
	return;

    if (!(note = io->db.notes))
	return;

    while(note) {
	note_read(io, note, n);
	if (n.type == itype && n.annotation) {
	    char *txt;

	    if (NULL == (txt = TextAllocRead(io, n.annotation)))
		break;

	    if (Tcl_GlobalEval(GetInterp(), txt) != TCL_OK) {
		verror(ERR_WARN, "execute_database_note",
		       "Note '%s' failed with message \"%s\"",
		       type, GetInterpResult());
	    }

	    xfree(txt);
	}
	note = n.next;
    }
}

void process_rawdata_note(GapIO *io) {
    GNotes n;
    int note;
    int itype = str2type("RAWD");
    static char *orig_path = NULL;
    static int orig = 0;

    if (!rawdata_note)
	return;

    /*
     * Purify flags this xmalloc as a memory leak. However see the linux
     * manpage on putenv() for a good description of the different putenv()
     * standards.
     * The code below is corect and backwards compatible with older putenv()
     * semantics. (It is better to leak a few bytes occasionally instead of
     * corrupting the environment.)
     */
    if (!orig) {
	char *p;
	orig = 1;
	if (p = getenv("RAWDATA")) {
	    if (NULL == (orig_path = xmalloc(strlen(p)+100)))
		return;
	    sprintf(orig_path, "RAWDATA=%s", p);
	}
    }

    if (!(note = io->db.notes)) {
	/* Revert to original RAWDATA setting */
	if (orig_path)
	    putenv(orig_path);
	else
	    putenv("RAWDATA=.");
	return;
    }

    while(note) {
	note_read(io, note, n);
	if (n.type == itype && n.annotation) {
	    char *rawd, *p;
	    char *buf;

	    if (!(rawd = TextAllocRead(io, n.annotation)))
		break;

	    /* Sanitise things, incase of any errors */
	    for (p = rawd; *p; p++) {
		if (*p == '\n' || *p == '\r') {
		    *p = 0;
		    break;
		}
		if (!(isalnum(*p) || ispunct(*p) || isspace(*p))) {
		    verror(ERR_WARN, "rawdata_note", "Malformed RAWD note");
		    xfree(rawd);
		    return;
		}
	    }
	    if (NULL == (buf = xmalloc(strlen(rawd) + 100))) {
		xfree(rawd);
		break;
	    }

	    sprintf(buf, "RAWDATA=%s", rawd);
	    putenv(buf);
	    xfree(rawd);
	}
	note = n.next;
    }
}

/*
 * Fix any notes created during the Gap4 4.4beta days. These had only
 * 32-bit time_t values, which isn't compliant to the newer standards
 * introduced in Solaris 7.
 */
void fix_notes(GapIO *io) {
    GViewInfo vi;
    int view;
    GNotes n;
    int i;

    if (io->db.Nnotes == 0)
	return;

    view = arr(GView, io->views, arr(GCardinal, io->notes, 0));
    if (view == -INT_MAX) {
	puts("View not found");
	return;
    }
    g_view_info(io->client, view, &vi);

    if (vi.used != 32)
	return;

    for (i = 1; i <= Nnotes(io); i++) {
	note_read(io, i, n);

	n.prev_type = n.next;
	n.prev = n.annotation;
	n.next = n.mtime;
	n.annotation = n.mtime_top;
	n.mtime = n.ctime;
	n.mtime_top = 0;
	n.ctime = n.ctime_top;
	n.ctime_top = 0;

	note_write(io, i, n);
    }

    return;
}

/*
 * Invokes the note selector for a specified reading, contig or db.
 * gtype is one of GT_Readings, GT_Contigs or GT_Database.
 * num is the read/contig number (ignored for GT_Database).
 */
void select_note(GapIO *io, int gtype, int num) {
    char *type;
    char ident[100];
    char cmd[1024];

    switch (gtype) {
    case GT_Readings:
	type = "reading";
	sprintf(ident, "#%d", num);
	break;
    case GT_Contigs:
	type = "contig";
	sprintf(ident, "=%d", num);
	break;
    default:
	type="database";
	*ident = 0;
    }

    sprintf(cmd, "NoteSelector %d %s %s", *handle_io(io), type, ident);
    if (Tcl_Eval(GetInterp(), cmd) != TCL_OK)
	verror(ERR_WARN, "select_note", "%s\n", GetInterpResult());
}

/*
 * ---------------------------------------------------------------------------
 * Note and experiment files; parsing, IO, etc.
 * ---------------------------------------------------------------------------
 */

/*
 * Converts a time_t value to a string form.
 */
char *time_t2str(time_t t) {
    char buf2[1024];
    static char buf[1024];

    strftime(buf2, sizeof(buf2)-1, "%c %Z", localtime(&t));
    sprintf(buf, "%s (%ld)", buf2, (long)t);

    return buf;
}


/*
 * Converts the string output of time_t2str() into a time_t value.
 */
time_t str2time_t(char *str) {
    char *p;
    time_t t;
    long l;

    if (p = strchr(str, '(')) {
	/* Parse the "("seconds")" section */
	sscanf(p+1, "%ld", &l);
	t = l;
    } else {
#ifndef NO_STRPTIME
	struct tm tm;
	char *strptime(const char *buf,
		       const char *format,
		       struct tm *tm);
	/* Use strptime instead, if it exists */
	memset(&tm, 0, sizeof(tm));
	strptime(str, "%c %Z", &tm);
	t = mktime(&tm);
#else
	/* otherwise use the current time */
	t = time(NULL);
#endif
    }

    return t;
}

/*
 * Converts a note into a string format.
 * The format is:
 *
 * TYPE ctime=<str_time> (time_t)
 * mtime=<str_time> (time_t)
 * from=[database|reading <name>|contig <name>]
 * comment=<comment>
 *
 * Eg within an exp file:
 * NT   REFT ctime=Mon Mar 25 14:28:33 2002 GMT (1017066513)
 * NT        mtime=Mon Mar 25 14:28:33 2002 GMT (1017066513)
 * NT        from=reading wtr
 * NT        comment=control -ve
 */
char *note2str(GapIO *io, GNotes n, int source_type, int source_num) {
    char *comment, *str, *strp;
    char type[5];
    char ctime[100], mtime[100];


    /* Read note comment */
    if (n.annotation) {
	if (NULL == (comment = TextAllocRead(io, n.annotation)))
	    return NULL;
    } else {
	comment = NULL;
    }

    /* read type */
    type[0] = (n.type >> 030) & 0xff;
    type[1] = (n.type >> 020) & 0xff;
    type[2] = (n.type >> 010) & 0xff;
    type[3] = (n.type >> 000) & 0xff;
    type[4] = '\0';

    /* Lots of room for dates, type, etc */
    strp = str = xmalloc(2*(comment ? strlen(comment): 0) + 1000);
    if (str == NULL)
	return NULL;

    /* Convert c/m times. */
    strcpy(ctime, time_t2str(n.ctime));
    strcpy(mtime, time_t2str(n.mtime));

    strp += sprintf(str, "%s ctime=%s\nmtime=%s",
		    type, ctime, mtime);

    switch(source_type) {
    case GT_Database:
	strp += sprintf(strp, "\nfrom=database");
	break;

    case GT_Readings:
	strp += sprintf(strp, "\nfrom=reading %s",
			get_read_name(io, source_num));
	break;

    case GT_Contigs:
	strp += sprintf(strp, "\nfrom=contig %s",
			get_contig_name(io, source_num));
	break;
    }

    if (comment) {
	char *c2, *c2p, *commentp;

	/* Replace \n in comment with \\\n. */
	if (strchr(comment, '\n')) {
	    if (!(c2 = xmalloc(strlen(comment) * 2)))
		return NULL;
	    commentp = comment;
	    c2p = c2;
	    while(*commentp) {
		if (*commentp == '\n') {
		    *c2p++ = '\\';
		}
		*c2p++ = *commentp++;
	    }
	    *c2p = '\0';
	} else {
	    c2 = comment;
	}

	strp += sprintf(strp, "\ncomment=%s", c2);
	if (c2 != comment)
	    xfree(c2);
	xfree(comment);
    }

    /* Shrink some of the unused space */
    str = xrealloc(str, strlen(str)+1);

    return str;
}

/*
 * Parses a note in format 'str' as defined by the above note2str
 * function back into the constituent parts. All OUTPUT pointers have
 * to be non-NULL. The comment pointer points into the supplied 'str'
 * pointer so it should not be freed.
 *
 * Returns 0 for success,
 *        -1 for failure.
 */
int str2note(/* INPUT */
	     GapIO *io, char *str,
	     /* OUTPUT */
	     int *type,
	     time_t *ctime, time_t *mtime,
	     int *source_type, int *source_number,
	     char **comment) {
    time_t t;
    char *cp, *cp2, tmp;
    char type_s[1024], name_s[1024];

    /* type */
    *type = str2type(str);
    cp = str+5;

    *source_number = 0;
    *source_type = 0;
    *comment = NULL;
    *ctime = 0;
    *mtime = 0;

    while (*cp) {
	if (strncmp("ctime=", cp, 6) == 0) {
	    /* ctime */
	    if (NULL == (cp = strchr(cp, '(')))
		return -1;
	    cp++;
	    sscanf(cp, "%ld", (long*)&t);
	    *ctime = t;
	}

	else if (strncmp("mtime=", cp, 6) == 0) {
	    /* mtime */
	    if (NULL == (cp = strchr(cp, '(')))
		return -1;
	    cp++;
	    sscanf(cp, "%ld", (long*)&t);
	    *mtime = t;
	}

	else if (strncmp("from=", cp, 6) == 0) {
	    if (NULL == (cp = strchr(cp, '\n')))
		return -1;
	    cp++;
	    *name_s = 0;
	    if (sscanf(cp, "from=%s %s\n", type_s, name_s) < 1)
		return -1;
	    if (NULL == (cp2 = strchr(cp, '\n')))
		return -1;
	    tmp = *cp2;
	    *cp2 = 0;
	    *cp2 = tmp;

	    if (strcmp(type_s, "database") == 0) {
		*source_type = GT_Database;
	    } else if (strcmp(type_s, "reading") == 0) {
		*source_type = GT_Readings;
		if (*name_s)
		    *source_number = get_gel_num(io, name_s, GGN_ID);
	    } else if (strcmp(type_s, "contig") == 0) {
		*source_type = GT_Contigs;
		if (*name_s)
		    *source_number = get_gel_num(io, name_s, GGN_ID);
		if (*source_number)
		    *source_number = rnumtocnum(io, *source_number);
	    } else {
		return -1;
	    }
	    *cp2 = tmp;
	    cp = cp2+1;
	}

	else if (strncmp(cp, "comment=", 8) == 0) {
	    /* comment */
	    *comment = cp+8;
	    break; /* Comment is up to the end of the record, not line */
        }

	cp++;
	while (*cp && *cp != '\n')
	    cp++;
	if (*cp == '\n')
	    cp++;
    }

    return 0;
}

/*
 * Adds a note in string form (str) to the specified reading.
 *
 * Returns 0 for success,
 *        -1 for failure.
 */
int create_note_for_gel(GapIO *io, int rnum, char *str) {
    int type, source_type, source_num;
    time_t ctime, mtime;
    char *comment;
    int nnote;

    if (str2note(io, str,
		 &type, &ctime, &mtime, &source_type, &source_num,
		 &comment) == -1) {
	verror(ERR_WARN, "create_note_for_gel",
	       "Malformed note '%s'", str);
	return -1;
    }
    if (!source_type)
	source_type = GT_Readings;
    if (!source_num)
	source_num = rnum;

    if ((nnote = new_note(io, type, source_type, source_num)) == -1)
	return -1;

    if (comment)
	edit_note(io, nnote, NULL, comment);

    set_note_time(io, nnote, ctime, mtime);

    return 0;
}

/*
 * Searches for a note of a given type.
 * Returns note first number found of that type
 *         0 if not found.
 */
int find_note(GapIO *io, int rnum, char *str_type) {
    GReadings r;
    GNotes n;
    int nn;
    int itype = str2type(str_type);

    gel_read(io, rnum, r);
    for (nn = r.notes; nn; nn = n.next) {
	note_read(io, nn, n);
	if (n.type == itype)
	    return nn;
    }

    return 0;
}
