/*  Last edited: Jan  7 10:35 2004 (mng) */
#include <tk.h>

#include "io-reg.h"
#include "tcl_utils.h"
#include "edStructs.h"
#include "tkEditor.h"
#include "tkEdNames.h"
#include "edUtils.h"
#include "contigEditor.h"
#include "IO.h"
#include "locks.h"
#include "tkSheet.h"
#include "fort.h"
#include "gap_cli_arg.h"
#include "gap_globals.h"
#include "misc.h"
#include "tclXkeylist.h"
#include "select.h"
#include "tman_interface.h"

/*
 * ============================================================================
 * Tk interface bits and pieces. Parses the command line arguments and
 * calls the appropriate C-only interfaces.
 * ============================================================================
 */

typedef struct {
    GapIO *io;
    char *contig;
    char *reading;
    int pos;
    int reuse;
    int nojoin;
} ec_arg;

int tk_edit_contig(ClientData clientData, Tcl_Interp *interp,
		   int argc, char **argv) {
    ec_arg args;
    int contig, reading = 0;
    cli_args a[] = {
	{"-io",      ARG_IO,  1, NULL, offsetof(ec_arg, io)},
	{"-contig",  ARG_STR, 1, NULL, offsetof(ec_arg, contig)},
	{"-reading", ARG_STR, 1, "",   offsetof(ec_arg, reading)},
	{"-pos",     ARG_INT, 1, "1",  offsetof(ec_arg, pos)},
	{"-reuse",   ARG_INT, 1, "0",  offsetof(ec_arg, reuse)},
        {"-nojoin",  ARG_INT, 1, "0",  offsetof(ec_arg, nojoin)},
	{NULL,      0,       0, NULL, 0}
    };

    vfuncheader("edit contig");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    if ((contig = get_contig_num(args.io, args.contig, GGN_ID)) < 0)
	return TCL_ERROR;

    if (*args.reading)
	reading = get_gel_num(args.io, args.reading, GGN_ID);
    if (reading <= 0)
	reading = io_clnbr(args.io, contig);

    if (args.reuse) {
      int id;
      
      if (-1 != (id = editor_available(contig, args.nojoin))) {
        if (*args.reading)
          move_editor(id, reading, args.pos);
        else
          move_editor(id, 0, args.pos);

        Tcl_SetResult(interp, Tk_PathName(EDTKWIN(editor_id_to_edstruct(id)->ed)), NULL);
        return TCL_OK;
      }
    }
    
    return edit_contig(interp, args.io, contig, reading, args.pos,
		       consensus_cutoff, quality_cutoff, 0);
}

typedef struct {
    GapIO *io;
    char *contig[2];
    char *reading[2];
    int pos[2];
} jc_arg;

int tk_join_contig(ClientData clientData, Tcl_Interp *interp,
		   int argc, char **argv) {
    jc_arg args;
    int contig[2], reading[2], i;
    cli_args a[] = {
	{"-io",       ARG_IO,  1, NULL, offsetof(jc_arg, io)},
	{"-contig1",  ARG_STR, 1, NULL, offsetof(jc_arg, contig[0])},
	{"-reading1", ARG_STR, 1, "",   offsetof(jc_arg, reading[0])},
	{"-pos1",     ARG_INT, 1, "1",  offsetof(jc_arg, pos[0])},
	{"-contig2",  ARG_STR, 1, NULL, offsetof(jc_arg, contig[1])},
	{"-reading2", ARG_STR, 1, "",   offsetof(jc_arg, reading[1])},
	{"-pos2",     ARG_INT, 1, "1",  offsetof(jc_arg, pos[1])},
	{NULL,      0,       0, NULL, 0}
    };

    vfuncheader("join contigs");

    if (-1 == gap_parse_args(a, &args, argc, argv))
	return TCL_ERROR;

    for (i = 0; i < 2; i++) {
	if ((contig[i] = get_contig_num(args.io, args.contig[i], GGN_ID)) < 0)
	    return TCL_ERROR;

	reading[i] = 0;
	if (*args.reading[i])
	    reading[i] = get_gel_num(args.io, args.reading[i], GGN_ID);
	if (reading[i] <= 0)
	    reading[i] = io_clnbr(args.io, contig[i]);
    }

    return join_contig(interp, args.io, contig, reading, args.pos,
		       consensus_cutoff, quality_cutoff);
}

int Ced_Init(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "edit_contig", tk_edit_contig,
		      NULL, NULL);

    Tcl_CreateCommand(interp, "join_contig", tk_join_contig,
		      NULL, NULL);

    Tcl_CreateCommand(interp, "edid_to_editor", tk_edid_to_editor,
		      NULL, NULL);

    /* Initialise semaphoreSystem */
    activeLock = semaphoreCreate(65535);

    return TCL_OK;
}


/*
 * ============================================================================
 * Actual editor stuff (global; for access by edUtils etc)
 * ============================================================================
 */

int inJoinMode(EdStruct *xx) {
    return xx->editorMode == JOINMODE;
}

int editorLocked(EdStruct *xx) {
    if (xx->editorMode == JOINMODE && xx->link) {
	return xx->link->locked;
    } else
	return 0;
}

/*
 * Find out the locked position.
 */
int editorLockedPos(EdStruct *xx[2], int force) {
    if (force) {
        return xx[1]->displayPos - xx[0]->displayPos;
    } else
        return (xx[0]->link->lockOffset);
}


void bell(void) {
    extern Tcl_Interp *GetInterp(void);

    Tcl_Eval(GetInterp(), "bell");
}

static char *next_editor(Tcl_Interp *interp) {
    static int editor_id = 0;
    static char buf[100];
    char *pname, *spec;

    pname = get_default_string(interp, gap_defs, "CONTIG_EDITOR.WIN");
    sprintf(buf, "%s%d", pname, editor_id++);

    /*
     * Initialise auto-display-traces setup. It's done here as this is a
     * convenient location where we know the contents of gap_defs, and it's
     * executed before bringing up the editor. (Although strictly speaking
     * initialising just once would suffice.)
     */
    spec = get_default_string(interp, gap_defs,
			      "CONTIG_EDITOR.AUTO_DISPLAY_TRACES_CONF");
    if (spec)
	tman_init_problem_traces(spec);

    return buf;
}


/*
 *----------------------------------------------------------------------------
 * The DB callback mechanism
 *---------------------------------------------------------------------------
 */
void db_callback_tk(void *xxv, int type, int seq, int pos, void *pointer) {
    EdStruct *xx = xxv;

    if (xx->editorState == StateDown)
	return;

    switch (type) {
    case DBCALL_REDISPLAY:
	tk_redisplaySequences(xx);
	break;

    case DBCALL_INSERT:
	selectInsertBase(xx, seq, pos);
	break;

    case DBCALL_DELETE:
	selectDeleteBase(xx, seq, pos);
	break;

    case DBCALL_CURSOR:
	setCursorPosSeq(xx, pos, seq);
	redisplayWithCursor(xx);
	break;

    case DBCALL_ADJUST_START:
	xx->displayPos += pos;
	break;

    case DBCALL_REINIT: {
	GapIO *io = DBI_io(xx);
	GContigs c;

	invalidate_consensus(xx);

	/*
	 * If we've already got this contig registered then we don't need
	 * to re-register / re-initialise. This occurs when we have multiple
	 * editors for the same contig. A single deregister is performed,
	 * but DBCALL_REINIT is needed for each editor viewing this contig
	 * (to set the extents).
	 *
	 * To prevent multiple registering, before calling REINIT we negate
	 * the saved registration id. This is then renegated only once.
	 * This also means the the registration id for the contig editor does
	 * not change, so any other function that has queried and stored this
	 * value does not need updating. A hacky but trivial solution.
	 */
	if (DBI_registration_id(xx) < 0) {
	    DBI_registration_id(xx) = -DBI_registration_id(xx);
	    contig_read(io, DBI_contigNum(xx), c);
	    initialiseDB(xx, io, DBI_contigNum(xx), io_dbsize(io), c.left);
	}

	getExtents(xx);

	xx->refresh_flags |= ED_DISP_ALL;
	redisplaySequences(xx, 0);
	break;
    }

    case DBCALL_JOIN_SHIFT:
	setDisplayPosP(xx, xx->displayPos + pos);
	setCursorPosSeq(xx,
			positionInContig(xx, xx->cursorSeq, xx->cursorPos) +
			pos, 0);
	if (xx->cursorPos > DB_Length(xx, 0))
	    setCursorPos(xx, DB_Length(xx, 0));

	invalidate_consensus(xx);

	break;

    case DBCALL_QUIT:
	/* db_callback_tk(xx, DBCALL_CURSOR_NOTIFY, -1, 0, NULL); */

	if (xx->link) {
	    xx->link->xx[0]->editorState = StateDown;
	    xx->link->xx[1]->editorState = StateDown;
	    delete_contig_cursor(DBI_io(xx->link->xx[0]),
				 DBI_contigNum(xx->link->xx[0]),
				 xx->link->xx[0]->cursor->id, 1);
	    delete_contig_cursor(DBI_io(xx->link->xx[1]),
				 DBI_contigNum(xx->link->xx[1]),
				 xx->link->xx[1]->cursor->id, 1);
        } else {
	    delete_contig_cursor(DBI_io(xx), DBI_contigNum(xx),
				 xx->cursor->id, 1);
	    xx->editorState = StateDown;
	}
	if (TCL_OK != (Tcl_VarEval(EDINTERP(xx->ed), "editor_quit_internal ",
				   Tk_PathName(EDTKWIN(xx->ed)), NULL)))
	    fprintf(stderr, "%s\n", EDINTERP(xx->ed)->result);

	if (xx->link) {
	    tman_shutdown_traces(xx->link->xx[0], 0);
	    tman_shutdown_traces(xx->link->xx[1], 0);
	} else {
	    tman_shutdown_traces(xx, 0);
	}

	break;

    case DBCALL_CURSOR_NOTIFY:
	{
	    reg_cursor_notify cn;
	    dstring_t *ds;
	    char var[1024];

	    /*
	     * Update cursor structure using a (seq,pos) pair.
	     * (0,0)  expands to current cursor and sequence number.
	     * otherwise it's the coordinates.
	     */
	    if (0 == seq && 0 == pos) {
		seq = xx->cursorSeq;
		pos = xx->cursorPos;
	    }
	    xx->cursor->seq = seq;
	    xx->cursor->pos = pos;
	    xx->cursor->abspos = positionInContig(xx, xx->cursorSeq, pos);
	    xx->cursor->job = CURSOR_MOVE;
	    xx->cursor->sent_by = DBI_registration_id(xx);

	    /* Send the notification */
	    cn.cursor = xx->cursor;
	    cn.job = REG_CURSOR_NOTIFY;
	    contig_notify(DBI_io(xx), DBI_contigNum(xx), (reg_data *)&cn);

	    /* Update the Tcl 'tags under cursor' variable. */
	    ds = listAnnotation(xx);
	    sprintf(var, "%s.Tags", Tk_PathName(EDTKWIN(xx->ed)));
	    Tcl_SetVar(EDINTERP(xx->ed),
		       var,
		       dstring_str(ds),
		       TCL_GLOBAL_ONLY|TCL_LEAVE_ERR_MSG);
	    dstring_destroy(ds);
	}

	break;

    case DBCALL_RELINK:
	/* Unlink xx from the DB and free DB if unused by anything else */
	freeDB(xx, 0);

	/* Link xx to the new DB */
	DBI(xx) = (DBInfo *)pointer;
	DBI_dispFunc(xx)[DBI_nextDisp(xx)] = db_callback_tk;
	DBI_dispData(xx)[DBI_nextDisp(xx)] = xx;
	DBI_nextDisp(xx)++;

	break;

    default:
	verror(ERR_FATAL, "db_callback_tk",
	       "Unknown callback - %d, seq %d, pos %d\n", type, seq, pos);
    }
}

int edit_contig(Tcl_Interp *interp, GapIO *io, int cnum, int llino, int pos,
		float con_cut, int qual_cut, int reveal_cutoffs) {
    EdStruct *xx;
    char ccut[10], qcut[10], rev[10], *edname, dbptr[50];
    int i;
    Tcl_CmdInfo cmdinfo;
    char *ptr;

    sprintf(ccut, "%d", (int)(con_cut * 100 + 0.1));
    sprintf(qcut, "%d", qual_cut);
    sprintf(rev,  "%d", reveal_cutoffs);

    /* Create edstruct */
    if (NULL == (xx = getFreeEdStruct(io, cnum, db_callback_tk)))
	return TCL_ERROR;

    sprintf(dbptr, "%p", (void *)DBI(xx));

    if (TCL_OK != Tcl_VarEval(interp, "create_editor ",
			      edname = next_editor(interp), /* toplevel name */
			      " 0",                /* editor subname */
			      " 0 ",               /* is not a join editor */
			      rev, " ",            /* reveal cutoffs */
			      ccut, " ", qcut, " ",
			      dbptr,
			      NULL)) {
	fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
    }

    ptr = strchr(interp->result, ' ');
    if (ptr)
	*ptr++ = 0;

    if (0 == Tcl_GetCommandInfo(interp, interp->result, &cmdinfo)) {
	verror(ERR_FATAL, "edit_contig", "No Editor structure!");
	return TCL_ERROR;
    }
    xx->ed = (Editor *)cmdinfo.clientData;
    xx->ed->xx = xx;

    if (0 == Tcl_GetCommandInfo(interp, ptr, &cmdinfo)) {
	verror(ERR_FATAL, "edit_contig", "No Names structure!");
	return TCL_ERROR;
    }
    xx->names = (edNames *)cmdinfo.clientData;
    xx->names->xx = xx;

    xx->editorState = StateUp;
    xx->editorMode = EDITMODE;

    xx->cursor = create_contig_cursor(io, cnum, 1, 0);

    /*
     * Set up data structures
     */
    if (DBI_nextDisp(xx) < 2 &&
	initialiseDB(xx, io, cnum, io->db.actual_db_size, io_clnbr(io,cnum))) {
	return TCL_ERROR;
    }

    xx->cursor->sent_by = DBI_registration_id(xx);

    xx->con_cut = con_cut;
    xx->qual_cut = qual_cut;
    for (i=0; i<10; i++)
	xx->qual_bg[i] = xx->ed->qual_bg[i]->pixel;
    for (i=0; i<4; i++)
	xx->edit_bg[i] = xx->ed->edit_bg[i]->pixel;
    for (i=0; i<5; i++)
	xx->tmpl_bg[i] = xx->ed->tmpl_bg[i]->pixel;
    xx->qual_below = xx->ed->qual_below->pixel;
    xx->diff_bg = xx->ed->diff_bg->pixel;

    getExtents(xx);
    /*
     * Yuk - why can't we get the default tag list within C without having
     * to resort to Tcl_VarEval!
     */
    if (TCL_OK !=
	Tcl_VarEval(interp, "eval ",
		    Tk_PathName(EDTKWIN(xx->ed)),
		    " set_displayed_annos [GetDefaultTags CONTIG_EDITOR.TAGS ",
		    Tk_PathName(EDTKWIN(xx->ed)),
		    "]",
		    NULL)) {
	fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
    }

    if (TCL_OK != Tcl_VarEval(interp, "wm title ",
			      " [winfo toplevel ",
			      Tk_PathName(EDTKWIN(xx->ed)),
			      "] {Contig Editor: ",
			      edGetGelName(xx, 1),
			      "}",
			      NULL)) {
	fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
    }

    if (reveal_cutoffs)
	edSetRevealCutoffs(xx, 1);
    else {
	/*
	 * Not ideal. This assumes pos relates to the consensus sequence
	 * (which it may not) and that cutoffs aren't shown.
	 * The ideal fix is to edit createEdDisplay().
	 * However, this fixes just enough to remove bugs triggered using the
	 * template display.
	 */
	if (pos < 1)
	    pos = 1;
	if (pos > io_clength(io, cnum)+1)
	    pos = io_clength(io, cnum)+1;
    }
    createEdDisplay(xx, llino, pos);

    if (TCL_OK != Tcl_VarEval(interp, "init_editor_states ",
			      edname, ".0 ",
			      Tk_PathName(EDTKWIN(xx->ed)),
			      " ", dbptr,
			      NULL)) {
	fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
    }

    {
	char c_io[10];
	sprintf(c_io, "%d", *handle_io(io));
	Tcl_VarEval(interp, "SelectReadingList ", c_io, NULL);
    }

    Tcl_SetResult(interp, Tk_PathName(EDTKWIN(xx->ed)), NULL);

    return TCL_OK;
}


int join_contig(Tcl_Interp *interp, GapIO *io, int cnum[2], int llino[2],
		int pos[2], float con_cut, int qual_cut) {
    EdStruct *xx[2];
    int i, j;
    char edn[100];
    char ccut[10], qcut[10], rev[10], dbptr[50];
    int reveal_cutoffs[2];
    Tcl_CmdInfo cmdinfo;
    char *ptr;
    int seq;

    strcpy(edn, next_editor(interp));
    sprintf(ccut, "%d", (int)(con_cut * 100 + 0.1));
    sprintf(qcut, "%d", qual_cut);

    /* Create edstructs */
    for (i=0; i<2; i++) {
	if (NULL == (xx[i] = getFreeEdStruct(io, cnum[i], db_callback_tk)))
	    return TCL_ERROR;

	xx[i]->editorState = StateUp;
	xx[i]->editorMode = JOINMODE;
	xx[i]->cursor = create_contig_cursor(io, cnum[i], 1, 0);

	/*
	 * Set up data structures
	 */
	if (DBI_nextDisp(xx[i]) < 2 &&
	    initialiseDB(xx[i], io, cnum[i], io->db.actual_db_size,
			 io_clnbr(io, cnum[i]))) {
	    return TCL_ERROR;
	}

	xx[i]->cursor->sent_by = DBI_registration_id(xx[i]);

	xx[i]->con_cut = con_cut;
	xx[i]->qual_cut = qual_cut;

	if (pos[i] <= 0 || pos[i] > io_clength(io, cnum[i]))
	    reveal_cutoffs[i] = 1;
	else
	    reveal_cutoffs[i] = 0;
    }
    /* None or both, not one only */
    if (reveal_cutoffs[0] || reveal_cutoffs[1]) {
	reveal_cutoffs[0] = 1;
	reveal_cutoffs[1] = 1;
    }

    if (NULL == CreateEdLink(xx[0], xx[1]))
	return TCL_ERROR;

    /* Create the tk bits */
    sprintf(rev,  "%d", reveal_cutoffs[0]);
    sprintf(dbptr, "%p", (void *)DBI(xx[0]));
    if (TCL_OK != Tcl_VarEval(interp, "create_editor ",
			      edn,   /* toplevel name */
			      " 0",  /* editor subname */
			      " 1 ", /* is a join editor */
			      rev, " ", /* reveal cutoffs */
			      ccut, " ", qcut, " ",
			      dbptr,
			      NULL)) {
	fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
    }

    ptr = strchr(interp->result, ' ');
    if (ptr)
	*ptr++ = 0;

    if (0 == Tcl_GetCommandInfo(interp, interp->result, &cmdinfo)) {
	verror(ERR_FATAL, "edit_contig", "No Editor structure!");
	return TCL_ERROR;
    }
    xx[0]->ed = (Editor *)cmdinfo.clientData;
    xx[0]->ed->xx = xx[0];

    if (0 == Tcl_GetCommandInfo(interp, ptr, &cmdinfo)) {
	verror(ERR_FATAL, "edit_contig", "No Names structure!");
	return TCL_ERROR;
    }
    xx[0]->names = (edNames *)cmdinfo.clientData;
    xx[0]->names->xx = xx[0];

    if (TCL_OK != Tcl_VarEval(interp, "create_editor_diff ", edn, " d",
			      " 0", NULL))
	puts(interp->result);

    if (0 == Tcl_GetCommandInfo(interp, interp->result, &cmdinfo)) {
	verror(ERR_FATAL, "edit_contig", "No 'diff' sheet structure!");
	return TCL_ERROR;
    }
    xx[0]->link->diffs = (tkSheet *)cmdinfo.clientData;

    sprintf(rev,  "%d", reveal_cutoffs[1]);
    sprintf(dbptr, "%p", (void *)DBI(xx[1]));
    if (TCL_OK != Tcl_VarEval(interp, "create_editor ",
			      edn,   /* toplevel name */
			      " 1",  /* editor subname */
			      " 1 ", /* is a join editor */
			      rev, " ", /* reveal cutoffs */
			      ccut, " ", qcut, " ",
			      dbptr,
			      NULL))
	puts(interp->result);

    ptr = strchr(interp->result, ' ');
    if (ptr)
	*ptr++ = 0;

    if (0 == Tcl_GetCommandInfo(interp, interp->result, &cmdinfo)) {
	verror(ERR_FATAL, "edit_contig", "No Editor structure!");
	return TCL_ERROR;
    }
    xx[1]->ed = (Editor *)cmdinfo.clientData;
    xx[1]->ed->xx = xx[1];

    if (0 == Tcl_GetCommandInfo(interp, ptr, &cmdinfo)) {
	verror(ERR_FATAL, "edit_contig", "No Names structure!");
	return TCL_ERROR;
    }
    xx[1]->names = (edNames *)cmdinfo.clientData;
    xx[1]->names->xx = xx[1];

    /* Needed to ensure that xx->displayWidth isn't fudged too much! */
    Tcl_Eval(interp, "update idletasks");

    /*
     * Due to the update, some awkward so'n'so may have already quitted the
     * first editor!
     */
    if (EDTKWIN(xx[0]->ed) == NULL || EDTKWIN(xx[1]->ed) == NULL)
	return TCL_OK;

    if (TCL_OK != Tcl_VarEval(interp, "wm title ",
			      " [winfo toplevel ",
			      Tk_PathName(EDTKWIN(xx[0]->ed)),
			      "] {Join Editor: ",
			      edGetGelName(xx[0], 1),
			      " ",
			      edGetGelName(xx[1], 1),
			      "}",
			      NULL)) {
	fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
    }

    for (i = 0; i < 2; i++) {
	for (j=0; j<10; j++)
	    xx[i]->qual_bg[j] = xx[i]->ed->qual_bg[j]->pixel;
	for (j=0; j<4; j++)
	    xx[i]->edit_bg[j] = xx[i]->ed->edit_bg[j]->pixel;
	for (j=0; j<4; j++)
	    xx[i]->tmpl_bg[j] = xx[i]->ed->tmpl_bg[j]->pixel;
	xx[i]->qual_below = xx[i]->ed->qual_below->pixel;
	xx[i]->diff_bg = xx[i]->ed->diff_bg->pixel;

	getExtents(xx[i]);
	if (TCL_OK !=
	    Tcl_VarEval(interp, "eval ", Tk_PathName(EDTKWIN(xx[i]->ed)),
		    " set_displayed_annos [GetDefaultTags CONTIG_EDITOR.TAGS ",
			Tk_PathName(EDTKWIN(xx[i]->ed)),
			"]",
			NULL)) {
	    fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
	}

	createEdDisplay(xx[i], llino[i], pos[i]);
	seq = rnum_to_edseq(xx[i], llino[i]);
	if (seq <= 0 || pos[i] > DB_Length(xx[i], seq) || pos[i] < 1)
	    seq = 0;
	if (reveal_cutoffs[i])
	    edSetRevealCutoffs(xx[i], 1);
	setCursorPosSeq(xx[i], pos[i], seq);
	if (seq == 0)
	    setDisplayPos(xx[i], pos[i]);

	sprintf(dbptr, "%p", (void *)DBI(xx[i]));
	if (TCL_OK != Tcl_VarEval(interp, "init_editor_states ",
				  edn, i == 0 ? ".0 " : ".1 ",
				  Tk_PathName(EDTKWIN(xx[i]->ed)),
				  " ", dbptr,
				  NULL)) {
	    fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
	}
    }

    if (xx[0]->link) {
	xx[0]->link->locked = 1;
	edSetJoinLock(xx[0], xx[0]->link->locked);
    } else {
	verror(ERR_FATAL, "join_contig", "link failed");
	return TCL_ERROR;
    }

    return TCL_OK;
}
