#ifndef _TK_IO_REG_H_
#define _TK_IO_REG_H_

#include <tcl.h>

/*
 * Create a list of results and associated contig numbers. The format of the
 * list is "{contig regnum id string} ?{contig regnum id string}? ..."
 *
 * The list is stored in the variable specified on the command line.
 */
int tk_result_names(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv);

/*
 * Tk interface to register_id() - all this work for such a miniscule 2 line
 * function!!!
 */
int tk_register_id(ClientData clientData, Tcl_Interp *interp,
		   int argc, char **argv);

/*
 * Tk interface to result_time()
 */
int tk_result_time(ClientData clientData, Tcl_Interp *interp,
		   int argc, char **argv);

/*
 * Fetches a list of operations for a particular result (id)
 */
int tk_reg_get_ops(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv);

/*
 * Performs an operation for a particular result (id)
 */
int tk_reg_invoke_op(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv);

/*
 * Notifies a contig of a general update. Contig '0' implies all contigs
 */
int tk_reg_notify_update(ClientData clientData, Tcl_Interp *interp,
			 int argc, char **argv);

/*
 * Sends a highlight reading notification
 */
int tk_reg_notify_highlight(ClientData clientData, Tcl_Interp *interp,
			    int argc, char **argv);

/*
 * Updates the result manager window. For efficiency reasons we only bother
 * updating when things become idle again.
 */
void update_results_(ClientData clientData);
void update_results(GapIO *io);


/*
 * Displays a tk dialogue box explaining that a contig is busy - it's all
 * to easy to not notice that comment in the error window.
 */
void busy_dialog(GapIO *io, int contig);


/*
 * Configures a result.
 */
int tk_matchresult_configure(ClientData clientData, Tcl_Interp *interp,
			     int argc, char **argv);

/*
 * Attempts to shut down all active displays.
 * As used by alter relationships, assembly, etc.
 */
int tcl_quit_displays(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv);

/*
 * Send a delete request to a specific result
 */
int tk_result_delete(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv);
/*
 * Send a quit request to a specific result
 */
int tk_result_quit(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv);

/*
 * Delete all contig comparator displays.
 */
int tk_clear_cp(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv);

/*
 * Delete all template plot displays.
 */
int tk_clear_template(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv);

/*
 * Delete all consistency displays.
 */
int tk_clear_consistency(ClientData clientData, Tcl_Interp *interp,
			 int argc, char **argv);

/*
 * Determines whether a result is a component of the 2D contig comparator
 * display. Returns 1 or 0.
 */
int tk_result_is_2d(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv);

/*
 * Determines whether a result is a component of the consistency
 * display. Returns 1 or 0.
 */
int tk_result_is_consistency(ClientData clientData, Tcl_Interp *interp,
			     int argc, char **argv);

/*
 * Register Tcl procedures with contig/notifications
 */
int tk_contig_register(ClientData clientData, Tcl_Interp *interp,
		       int argc, char **argv);

/*
 * Deregister Tcl procedures.
 */
int tk_contig_deregister(ClientData clientData, Tcl_Interp *interp,
			 int argc, char **argv);

/*
 * An arbitray Tcl interface to contig event notification
 */
int tk_contig_notify_obj(ClientData clientData, Tcl_Interp *interp,
			 int objc, Tcl_Obj **objv);
int tk_contig_notify(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv);

/*
 * An arbitray Tcl interface to contig event notification
 */
int tk_query_cursor(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv);

/*
 * A tcl interface to create_contig_cursor().
 */
int tk_create_cursor(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv);

/*
 * A tcl interface to delete_contig_cursor().
 */
int tk_delete_cursor(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv);

int tk_cursor_ref(ClientData clientData, Tcl_Interp *interp,
		  int argc, char **argv);

/*
 * Interface to contig_lock_write()
 */
int tk_contig_lock_write(ClientData clientData, Tcl_Interp *interp,
			 int argc, char **argv);

#endif /* _TK_IO_REG_H_ */
