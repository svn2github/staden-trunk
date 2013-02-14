/*
 * The EdNames tk widget.
 * As with the Editor widget this is an inheritance of the tkSheet widget.
 *
 * See tkEditor for better commenting of this method.
 */

#include "tkSheet_common.h"
#include "tkEdNames.h"
#include "tk_defs.h"
//#include "edUtils.h"
#include "xalloc.h"
#include "tcl_utils.h"
//#include "notes.h"

#undef TKSHEET
#define TKSHEET(en)   ((tkSheet *)(en))

/* ---- Local function prototypes ---- */
static int NamesCmd(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv);
static int NamesWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			   int argc, char **argv);
static int NamesConfigure(Tcl_Interp *interp, edNames *en,
			  int argc, char **argv, int flags);
static void EdNamesSheetExtension(ClientData clientData, int job, void *data);

/* ---- Local static data ---- */


#define offset(field) Tk_Offset(edNames, field)
static Tk_ConfigSpec configSpecs[] = {
#   include "tkSheet_config.h"
    {TK_CONFIG_STRING,
         "-xscrollcommand", "xScrollCommand", "ScrollCommand",
         "", offset(xScrollCmd),        TK_CONFIG_NULL_OK, NULL},
    {TK_CONFIG_END,
	 (char *)NULL,	(char *)NULL,	(char *)NULL,	(char *) NULL,
         0,	0,	NULL},
};
#undef offset

/* ---- Global functions ---- */

/*
 * Init.
 */
int EdNames_Init(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "ednames", NamesCmd,
		      NULL, (Tcl_CmdDeleteProc *)NULL);

    return TCL_OK;
}



/* ---- Local functions ---- */

/*
 * Our class command procedure. This is different from usual class command
 * procs as we actually use the SheetCmd() procedure and then carefully
 * subvert their actions to our own uses.
 */
static int NamesCmd(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {
    Tk_Window tkwin;
    edNames *en;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " pathName ?options?\"", (char *)NULL);
        return TCL_ERROR;
    }

    /*
     * Allocate
     */
    if (NULL == (en = (edNames *)xmalloc(sizeof(edNames)))) {
        return TCL_ERROR;
    }


    /*
     * Call common sheet initialisation code
     */
    tkwin = SheetCmdCommon(interp, Tk_MainWindow(interp),
			   configSpecs, (tkSheet *)en,
			   argv[1], "EdNames");
    if (NULL == tkwin) {
	xfree(en);
	return TCL_ERROR;
    }

    /*
     * Initialised rest of edNames structure.
     */
    en->xx = NULL;
    en->xScrollCmd = NULL;
    TKSHEET(en)->extensionData = (ClientData)en;
    TKSHEET(en)->extension = EdNamesSheetExtension;


    /*
     * Register our instance of the widget class
     */
    Tcl_CreateCommand(interp, Tk_PathName(tkwin),
                      NamesWidgetCmd, (ClientData)en,
                      (Tcl_CmdDeleteProc *)NULL);


    /*
     * And process our arguments - send them to Configure
     */
    if (NamesConfigure(interp, en, argc-2, argv+2, 0) != TCL_OK) {
        Tk_DestroyWindow(tkwin);
        return TCL_ERROR;
    }


    Tcl_SetResult(interp, Tk_PathName(tkwin), TCL_VOLATILE);
    return TCL_OK;
}





/*
 * Configure.
 */
static int NamesConfigure(Tcl_Interp *interp, edNames *en,
			  int argc, char **argv, int flags) {
    if (TCL_ERROR == SheetConfigureCommon(interp, (tkSheet *)en,
					  argc, argv, flags))
	return TCL_ERROR;

    if (en->xx) {
	ed_set_nslider_pos(en->xx, en->xx->names_xPos);
    }

    return TCL_OK;
}


/*
 * EdNames widget command
 */
static int NamesWidgetCmd(ClientData clientData, Tcl_Interp *interp,	
			  int argc, char **argv) {
    edNames *en = (edNames *)clientData;
    int result = TCL_OK;
    int index;
    Tcl_Obj *obj;
    edview *xx = en->xx;

    static char *optionStrings[] = {
	"configure", "io",   "xview",	"yview",  "get_number",
	"recs_between",
	NULL
    };
    enum options {
	CONFIGURE, IO, XVIEW,	YVIEW,     GET_NUMBER,
	RECS_BETWEEN,
    };

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " option ?arg arg ...?\"",
                         (char *) NULL);
        return TCL_ERROR;
    }

    obj = Tcl_NewStringObj(argv[1], -1);
    Tcl_IncrRefCount(obj);
    if (Tcl_GetIndexFromObj(interp, obj, optionStrings, "option", 0,
	    &index) != TCL_OK) {
	Tcl_DecrRefCount(obj);
	return TCL_ERROR;
    }
    Tcl_DecrRefCount(obj);

    switch ((enum options)index) {
    default:
	return TCL_ERROR;

    case CONFIGURE:
	result = SheetWidgetCmdConfig(interp, TKSHEET(en), argc, argv);
	break;
	
    case IO:
	Tcl_SetResult(interp, io_obj_as_string(xx->io) , TCL_VOLATILE);
	break;

    case XVIEW: {
	double f1;
	int type, count, offset;

	//type = Tk_GetScrollInfo(interp, objc, objv, &f1, &count);
	type = Tk_GetScrollInfo(interp, argc, argv, &f1, &count);
	switch (type) {
	default:
	    return TCL_ERROR;

	case TK_SCROLL_MOVETO:
	    offset = f1 * MAX_NAME_LEN + 0.5;
	    break;

	case TK_SCROLL_PAGES:
	    count *= 4;
	case TK_SCROLL_UNITS:
	    offset = xx->names_xPos + count;
	    break;
	}

	if (offset < 0)
	    offset = 0;
	xx->names_xPos = offset;
	xx->refresh_flags = ED_DISP_NAMES;
	edview_redraw(xx);

	ed_set_nslider_pos(xx, xx->names_xPos);
	break;
    }
	
    case YVIEW: {
	double f1;
	int type, count, offset;
	rangec_t *r;
	int nr;
	contig_t *c = cache_search(xx->io, GT_Contig, xx->cnum);

	/* Cache elsewhere - ie last number of seqs visible? */
	r = contig_seqs_in_range(xx->io, &c, xx->displayPos,
				 xx->displayPos + xx->displayWidth,
				 CSIR_COUNT_ONLY, &nr);
	free(r);

	//type = Tk_GetScrollInfoObj(interp, objc, objv, &f1, &count);
	type = Tk_GetScrollInfo(interp, argc, argv, &f1, &count);
	switch (type) {
	default:
	    return TCL_ERROR;

	case TK_SCROLL_MOVETO:
	    offset = f1 * nr;
	    break;

	case TK_SCROLL_PAGES:
	    offset = xx->displayYPos + count * 0.9 * xx->displayHeight;
	    break;

	case TK_SCROLL_UNITS:
	    offset = xx->displayYPos + count;
	    break;
	}

	if (offset >= 0)
	    xx->displayYPos = offset;
	xx->refresh_flags = ED_DISP_ALL;
	edview_redraw(xx);
	break;
    }

    /* Get reading number under at x,y coord (eg mouse pointer), or
     * the current editior cursor if no x,y specified.
     */
    case GET_NUMBER: {
	char buf[100];
	int x, y, type, pos;
	tg_rec rec;

	if (argc != 2 && argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_number ?xpos ypos?\"",
			     (char *) NULL);
	    return TCL_ERROR;
	}

	if (argc == 4) {
	    sheet_arg_x(TKSHEET(en), argv[2], &x); /* cell coordinates */
	    sheet_arg_y(TKSHEET(en), argv[3], &y); y++;

	    if (-1 != (type = edview_item_at_pos(xx, y, x, 1,
						 xx->ed->stack_mode,
						 1, &rec,&pos))) {
		sprintf(buf, "%d %"PRIrec" %d", type, rec, pos);
		Tcl_AppendResult(interp, buf, NULL);
	    } /* otherwise return a blank */
	} else {
	    sprintf(buf, "%d %"PRIrec" %d",
		    xx->cursor_type,
		    xx->cursor_rec,
		    xx->cursor_pos);
	    Tcl_AppendResult(interp, buf, NULL);
	}
	break;
    }

    case RECS_BETWEEN: {
	tg_rec from_rec, to_rec;
	Array a;
	int i;
	Tcl_Obj *rec_list;

	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " recs_between from_rec to_rec\"",
			     (char *) NULL);
	    return TCL_ERROR;
	}

	from_rec = atol(argv[2]);
	to_rec   = atol(argv[3]);

	a = edview_items_between(xx, from_rec, to_rec);
	if (!a)
	    return TCL_ERROR;

	rec_list = Tcl_NewListObj(0, NULL);
	for (i = 0; i < ArrayMax(a); i++) {
	    Tcl_Obj *o;
	    tg_rec rec = arr(tg_rec, a, i);

	    o = Tcl_NewWideIntObj(rec);
	    Tcl_ListObjAppendElement(interp, rec_list, o);
	}
	Tcl_SetObjResult(interp, rec_list);

	ArrayDestroy(a);
	break;
    }
    }

    return result;
}

#if 0
static int NamesWidgetCmd(ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv) {
    edNames *en = (edNames *)clientData;
    int result = TCL_OK;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " option ?arg arg ...?\"",
                         (char *) NULL);
        return TCL_ERROR;
    }

    Tcl_Preserve((ClientData)TKSHEET(en));

    if ('c' == *argv[1] && strcmp(argv[1], "configure") == 0) {
	result = SheetWidgetCmdConfig(interp, TKSHEET(en), argc, argv);

    } else if (!en->xx) {
	en->xx = 0; /* Do nothing and flow through to 'success' */

    } else if ('u' == *argv[1] && strcmp(argv[1], "update_brief") == 0) {
	int x, y;

	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " update_brief xpos ypos\"",
			     (char *) NULL);
	    goto fail;
	}
	
	sheet_arg_x(TKSHEET(en), argv[2], &x);
	sheet_arg_y(TKSHEET(en), argv[3], &y); y++;

	//edSetBriefNameStatus(en->xx, x, y);

    } else if ('h' == *argv[1] && strcmp(argv[1], "highlight") == 0) {
#if 0
	int y, x, mode = -1;
	if (argc != 4 && argc != 5) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " highlight mode ?xpos? ypos\"",
			     (char *) NULL);
	    goto fail;
	}

	Tcl_GetInt(interp, argv[2], &mode);

	if (argc == 5) {
	    sheet_arg_x(TKSHEET(en), argv[3], &x);
	    y = 4; /* argv number of y, rather than y itself */
	} else {
	    x = 2;
	    y = 3;
	}

	if (*argv[y] == '=') {
	    y = get_gel_num(DBI_io(en->xx), &argv[y][1], GGN_ID);
	} else {
	    sheet_arg_y(TKSHEET(en), argv[y], &y); y++;
	    y = edGetGelNumber(en->xx, 0, y);
	}

	if (-1 != y) {
	    if (x != 0) {
		edSelectRead(en->xx, y, mode);

		vTcl_SetResult(interp, "%d",
			       (DB_Flags(en->xx, y) & DB_FLAG_SELECTED) ?1:0);
	    } else {
		if (y)
		    select_note(DBI_io(en->xx), GT_Seq,
				DB_Number(en->xx,y));
		else
		    select_note(DBI_io(en->xx), GT_Contig,
				DBI_contigNum(en->xx));
	    }
	}

#endif
    } else if ('g' == *argv[1] && strcmp(argv[1], "get_number") == 0) {
	int x, y, num;

	if (argc != 2 && argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_number ?xpos ypos?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 4) {
	    sheet_arg_x(TKSHEET(en), argv[2], &x);
	    sheet_arg_y(TKSHEET(en), argv[3], &y); y++;

	    if (-1 != (num = edGetGelNumber(en->xx, x, y))) {
		vTcl_SetResult(interp, "%d", num);
	    } /* otherwise return a blank */
	} else {
	    vTcl_SetResult(interp, "%d", en->xx->cursorSeq);
	}

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_contig_number") == 0) {
	if (argc != 2) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_contig_number",
			     (char *) NULL);
	    goto fail;
	}

	vTcl_SetResult(interp, "%"PRIrec, en->xx->cnum);

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_name") == 0) {
	char *name;
	int num;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_name ?editor_gel_number?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3)
	    Tcl_GetInt(interp, argv[2], &num);
	else
	    num = en->xx->cursorSeq;

	Tcl_ResetResult(interp);
	if (NULL != (name = seq_name(en->xx->io, num))) {
	    Tcl_AppendResult(interp, name, NULL);
	} /* otherwise return a blank */

    } else if ('g' == *argv[1] && strcmp(argv[1], "get_names_to_right") == 0) {
	int num;
	dstring_t *ds = NULL;

	if (argc != 2 && argc != 3) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " get_name ?gel_number?\"",
			     (char *) NULL);
	    goto fail;
	}

	if (argc == 3)
	    Tcl_GetInt(interp, argv[2], &num);
	else
	    num = en->xx->cursorSeq;

	if (NULL != (ds = edGetGelNamesToRight(en->xx, num))) {
	    Tcl_SetResult(interp, dstring_str(ds), TCL_VOLATILE);
	    dstring_destroy(ds);
	} else {
	    /* otherwise return a blank */
	    Tcl_ResetResult(interp);
	}

    } else if ('x' == *argv[1] && strcmp(argv[1], "xview") == 0) {
	int offset, count, type;
	double fraction;

	if (argc == 2) {
	    /* query */
	    vTcl_SetResult(interp, "%d", en->xx->names_xpos);
	} else {
	    if (argc == 3) {
		/* old syntax */
		if (Tcl_GetInt(interp, argv[2], &offset) != TCL_OK) {
		    goto fail;
		}
		setNamePos(en->xx, offset);
	    } else {
		/* new syntax */
		type = Tk_GetScrollInfo(interp, argc, argv, &fraction, &count);
		switch (type) {
                case TK_SCROLL_ERROR:
                    goto fail;
                case TK_SCROLL_MOVETO:
                    offset = (int) (fraction*DB_NAMELEN + 0.5);
                    break;
                case TK_SCROLL_PAGES:
		    count *= 4;
                case TK_SCROLL_UNITS:
		    offset = en->xx->names_xpos + count;
                    break;
		}
		if (offset < 0)
		    offset = 0;
		if (offset + (en->sw.columns - (DB_GELNOLEN + 1)) > DB_NAMELEN)
		    offset = DB_NAMELEN - (en->sw.columns - (DB_GELNOLEN + 1));
		setNamePos(en->xx, offset);
	    }
	    
	    /* Update X scrollbar */
	    ed_set_nslider_pos(en->xx, en->xx->names_xpos);
	}

    } else if ('c' == *argv[1] && strcmp(argv[1], "coords") == 0) {
	int x, y;
	if (argc != 4) {
	    Tcl_AppendResult(interp, "wrong # args: should be \"",
			     argv[0], " update_brief xpos ypos\"",
			     (char *) NULL);
	    goto fail;
	}
	
	sheet_arg_x(TKSHEET(en), argv[2], &x);
	sheet_arg_y(TKSHEET(en), argv[3], &y); y++;

	vTcl_SetResult(interp, "%d %d", x, y);

    } else {
	Tcl_AppendResult(interp, "bad option \"", argv[1], "\": must be FIXME",
			 NULL);
	goto fail;
    }

    Tcl_Release((ClientData)TKSHEET(en));
    return result;

 fail:
    Tcl_Release((ClientData)TKSHEET(en));
    return TCL_ERROR;
}
#endif

/*
 * Callbacks from the sheet widget for our registered extension function.
 * See tkSheet.h for job numbers available and tkSheet.c for where they're
 * called.
 */
static void EdNamesSheetExtension(ClientData clientData, int job, void *data) {
    edNames *en = (edNames *)clientData;

    if (!en->xx)
	return;

    switch (job) {
    case SHEET_JOB_RESIZE:
	en->xx->refresh_flags |= ED_DISP_NAMES;
	edview_redraw(en->xx);
	ed_set_nslider_pos(en->xx, en->xx->names_xPos);
	break;
    }
}
