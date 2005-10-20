#ifndef _NOTES_H_
#define _NOTES_H_

#include <tk.h>
int tcl_new_note(ClientData clientData, Tcl_Interp *interp,
		 int objc, Tcl_Obj *CONST objv[]);

int tcl_delete_note(ClientData clientData, Tcl_Interp *interp,
		    int objc, Tcl_Obj *CONST objv[]);

int tcl_edit_note(ClientData clientData, Tcl_Interp *interp,
		  int objc, Tcl_Obj *CONST objv[]);


int new_note(GapIO *io, int ntype, int gtype, int num);

int delete_note(GapIO *io, int nnum);

int edit_note(GapIO *io, int nnum, char *type, char *comment);

int delete_note_list(GapIO *io, int nnum);

void execute_database_notes(GapIO *io, char *type);

void process_rawdata_note(GapIO *io);

void fix_notes(GapIO *io);

void select_note(GapIO *io, int gtype, int num);

/*
 * Merges the notes from two contigs.
 *
 * Returns 0 for success, -1 for failure.
 */
int merge_contig_notes(GapIO *io, int cfrom, int cto);

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
char *note2str(GapIO *io, GNotes n, int source_type, int source_num);

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
	     char **comment);

/*
 * Adds a note in string form (str) to the specified reading.
 *
 * Returns 0 for success,
 *        -1 for failure.
 */
int create_note_for_gel(GapIO *io, int rnum, char *str);

/*
 * Searches for a note of a given type.
 * Returns note first number found of that type
 *         0 if not found.
 */
int find_note(GapIO *io, int rnum, char *str_type);

#endif /* _NOTES_H_ */
