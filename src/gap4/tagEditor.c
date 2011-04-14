#include <stdio.h>
#include <tk.h>
#include <string.h>

#include "misc.h"
#include "edUtils.h"
#include "tagUtils.h"
#include "undo.h"
#include "gap_globals.h"
#include "select.h"
#include "tkEditor.h"

typedef struct {
    int status;		/* 0 = editing, 1 = saved, 2 = quit */
    char window[100];	/* Name of the tag editor window */
    char array[100]; 	/* Name of the tag data array in tcl*/
    char command[256]; 	/* Name of the tag command in tcl */

    EdStruct *xx;
    tagStruct *tag;	/* Current tag structure, NULL if none */
    tag_id id;
    int  seq;
    int  len;
    int  pos;
    int  sense;
    char type[5];
    char *anno;
} TagEd;


/*
 * Removes a tag editor and tidies up memory
 */
static void TagEdDestroy(Tcl_Interp *interp, TagEd *te) {
    Tcl_UnsetVar(interp, te->array, TCL_GLOBAL_ONLY);
    Tcl_DeleteCommand(interp, te->command);
    Tcl_VarEval(interp, "destroy ", te->window, NULL);

    xfree(te);
}


/*
 * Updates an existing tag or creates a new one
 */
static void TagEdSave(Tcl_Interp *interp, TagEd *te) {
    EdStruct *xx = te->xx;
    tagStruct *tag;
    char *ncomment;

    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return;
    }

    U_adjust_cursor(xx, 0);

    if (te->tag) {
	/* Update an existing tag */

	/* find tag prior to our edited one */
	tag = findPreviousTag(xx, te->seq, te->tag);

	U_delete_annotation(xx, te->seq, tag);

	ncomment = (char *)TAG_MALLOC(strlen(te->anno)+1);
	strcpy(ncomment, te->anno);
	U_create_annotation(xx, te->seq, te->pos, te->len, te->type,
			    ncomment, tag,
			    normaliseSense(xx, te->seq, te->sense));

    } else {
	/* Create a new tag */
	int pos;

	pos = normalisePos2(xx, te->seq, te->pos, te->len);
	
	/*
	 * find position to insert tag, and insert it there.
	 */
	tag = findTagPos(xx, te->seq, pos);

	ncomment = (char *)TAG_MALLOC(strlen(te->anno)+1);
	strcpy(ncomment, te->anno);
	U_create_annotation(xx, te->seq, pos, te->len, te->type,
			    ncomment, tag,
			    normaliseSense(xx, te->seq, te->sense));

	if (tag == NULL)
	    xx->select_tag = DBgetTags(DBI(xx), te->seq);
	else
	    xx->select_tag = tag->next;
    }
    U_adjust_cursor(xx, 0);

    redisplaySequences(xx, 1);
    DBsetFlags(xx, te->seq, DB_Flags(xx, te->seq) | DB_FLAG_TAG_MODIFIED);
}

/* Saves a point (one base) annotation */
int saveAnnotation(EdStruct *xx, char *type, char *anno, int strand) {
    tagStruct *tag;
    char *ncomment;
    int pos;
    int seq;
    int len = 1;

    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return -1;
    }

    openUndo(DBI(xx));

    if (! getSelection(xx, &seq, &pos, &len, &tag) || len == 0) {
	/* default selection is current cursor position */
	seq = xx->cursorSeq;
	pos = xx->cursorPos + DB_Start(xx, seq);
	len = 1;
	if (pos > DB_Length2(xx, seq)) {
	    bell();
	    return -1;
	}
    }

    /* Complement positions if required */
    pos = normalisePos2(xx, seq, pos, len);
    
    /*
     * find position to insert tag, and insert it there.
     */
    tag = findTagPos(xx, seq, pos);
    
    ncomment = (char *)TAG_MALLOC(strlen(anno)+1);
    strcpy(ncomment, anno);
    U_adjust_cursor(xx, 0);
    U_create_annotation(xx, seq, pos, len, type,
			ncomment, tag,
			normaliseSense(xx, seq, strand));
    U_adjust_cursor(xx, 0);

    if (tag == NULL)
	xx->select_tag = DBgetTags(DBI(xx), seq);
    else
	xx->select_tag = tag->next;

    redisplaySequences(xx, 1);
    DBsetFlags(xx, seq, DB_Flags(xx, seq) | DB_FLAG_TAG_MODIFIED);

    closeUndo(xx, DBI(xx));
    return 0;
}

/*
 * Callbacks from tcl/tk for tag manipulation.
 */
static int TagEdCommand(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    TagEd *te = (TagEd *)clientData;

    if (argc < 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " option\"",
                         (char *) NULL);
        return TCL_ERROR;
    }
    
    if (strcmp(argv[1], "quit") == 0) {
	TagEdDestroy(interp, te);

    } else if (strcmp(argv[1], "save") == 0 ||
	       strcmp(argv[1], "move") == 0 ||
	       strcmp(argv[1], "copy") == 0) {
	char *p;
	EdStruct *xx = te->xx;

	openUndo(DBI(xx));

	/* Actually more like a move to */
	if (strcmp(argv[1], "move") == 0 || strcmp(argv[1], "copy") == 0) {
	    int seq, pos, len;
	    Tcl_CmdInfo info;
	    EdStruct *xx2;

	    if (0 == Tcl_GetCommandInfo(interp, argv[2], &info)) {
		bell();
		return TCL_ERROR;
	    }
	    xx2 = ((Editor *)info.clientData)->xx;

	    /* Find the new selection */
	    if (! getSelection(xx2, &seq, &pos, &len, NULL)) len = 0;
	    if (!len) {
		seq = xx2->cursorSeq;
		pos = xx2->cursorPos + DB_Start(xx2, seq);
		len = 1;
		if (pos > DB_Length2(xx2, seq)) {
		    bell();
		    return TCL_ERROR;
		}
	    }

	    /* Remove the existing tag */
	    if (strcmp(argv[1], "move") == 0 && te->tag) {
		int tmp = xx->cursorSeq;
		U_adjust_cursor(xx, 0);
		xx->cursorSeq = te->seq;
		deleteAnnotation(xx, te->tag);
		te->tag = NULL;
		xx->cursorSeq = tmp;
	    }

	    /* And update the location for the new tag */
	    te->tag = NULL;
	    te->xx = xx2;
	    te->seq = seq;
	    te->pos = pos;
	    te->len = len;
	}

	p = Tcl_GetVar2(interp, te->array, "type", TCL_GLOBAL_ONLY);
	if (p)
	    strncpy(te->type, p, 4);
	else
	    fprintf(stderr, "Error at %s:%d\n", __FILE__, __LINE__);

	p = Tcl_GetVar2(interp, te->array, "anno", TCL_GLOBAL_ONLY);
	if (p)
	    te->anno = p;
	else
	    fprintf(stderr, "Error at %s:%d\n", __FILE__, __LINE__);

	p = Tcl_GetVar2(interp, te->array, "strand", TCL_GLOBAL_ONLY);
	if (p)
	    te->sense = atoi(p);
	else
	    fprintf(stderr, "Error at %s:%d\n", __FILE__, __LINE__);

	if (xx != te->xx) {
	    openUndo(DBI(te->xx));
	}
	TagEdSave(interp, te);
	if (xx != te->xx)
	    closeUndo(te->xx, DBI(te->xx));
	TagEdDestroy(interp, te);

	closeUndo(xx, DBI(xx));
    }

    return TCL_OK;
}


/*
 * Brings up a tag editor.
 * 
 * Returns the pathname of the tag editor window or NULL for failure.
 */
char *invokeTagEditor(EdStruct *xx, tag_id id, int seq, int pos, int length,
		      int sense, char *comment, char *type_id, tagStruct *tag)
{
    Tcl_Interp *interp = EDINTERP(xx->ed);
    char buf[2], *pname;
    TagEd *te;

    /* Setup our TagEd structure */
    if (NULL == (te = xmalloc(sizeof(TagEd))))
	return NULL;
    te->xx = xx;
    pname = Tk_PathName(EDTKWIN(xx->ed));
    if (tag) {
	sprintf(te->window,  "%s.tag%d%p",         pname, id, tag);
	sprintf(te->array,   "%s.tag%d%p.data",    pname, id, tag);
	sprintf(te->command, "%s.tag%d%p.command", pname, id, tag);
    } else {
	sprintf(te->window,  "%s.tag%d",         pname, id);
	sprintf(te->array,   "%s.tag%d.data",    pname, id);
	sprintf(te->command, "%s.tag%d.command", pname, id);
    }
    te->id = id;
    te->status = 0;
    te->pos = pos;
    te->len = length;
    te->sense = sense;
    strncpy(te->type, type_id, 4);
    te->type[4] = 0;
    te->anno = comment;
    te->seq = seq;
    te->tag = tag;

    /* Setup our tcl/tk data structure for this tag */
    Tcl_SetVar2(interp, te->array, "type", te->type, TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, te->array, "anno", te->anno, TCL_GLOBAL_ONLY);
    sprintf(buf, "%d", sense);
    Tcl_SetVar2(interp, te->array, "strand", buf, TCL_GLOBAL_ONLY);
    Tcl_SetVar2(interp, te->array, "default", tag ? "0" : "1",
		TCL_GLOBAL_ONLY);


    /* Create our editor window and add some commands to it */
    if (TCL_OK != Tcl_VarEval(interp, "create_tag_editor ", te->window,
			      " ", te->command, " ", te->array, " ", NULL))
	fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));

    Tcl_CreateCommand(interp, te->command, TagEdCommand, (ClientData)te, NULL);

    return te->window;
}

