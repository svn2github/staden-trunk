/* See end for discussion on this code */

#include <staden_config.h>

#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <tcl.h>
/* #include <tclInt.h> *//* Tcl_GetCommandFromObj */
#include <tcl_utils.h>

#include "array.h"
#include "misc.h"
#include "tg_gio.h"
#include "tg_check.h"
#include "gap_cli_arg.h"
#include "consensus.h"

extern Tcl_Command Tcl_GetCommandFromObj(Tcl_Interp *interp,
					 Tcl_Obj *objPtr);

/* ------------------------------------------------------------------------ */
/* Some standard argument types */
typedef struct {
    GapIO *io;
} io_arg;

typedef struct {
      GapIO *io;
      int contig;
} contig_arg;

typedef struct {
    GapIO *io;
    char *inlist;
} list2_arg;

/* ------------------------------------------------------------------------ */
/* Tcl_Obj "gapio" type implementation */

static int tcl_database_read(GapIO *io, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]);
static int tcl_contig_read(GapIO *io, Tcl_Interp *interp,
			   int objc, Tcl_Obj *CONST objv[]);
static int tcl_contig_order(GapIO *io, Tcl_Interp *interp,
			    int objc, Tcl_Obj *CONST objv[]);
static int tcl_sequence_read(GapIO *io, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]);
static int tcl_anno_ele_read(GapIO *io, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]);
static int tcl_library_read(GapIO *io, Tcl_Interp *interp,
			   int objc, Tcl_Obj *CONST objv[]);

static void io_update_string(Tcl_Obj *obj);
int io_from_any(Tcl_Interp *interp, Tcl_Obj *obj);

static Tcl_ObjType io_obj_type = {
    "gapio",
    (Tcl_FreeInternalRepProc*)NULL,
    (Tcl_DupInternalRepProc*)NULL,
    io_update_string,
    io_from_any
};

char *io_obj_as_string(GapIO *io) {
    static char buf[80];
    sprintf(buf, "io=%p", io);
    return buf;
}

/*
 * No access to the string representation, but for the sake of clarity we
 * produce a dummy string to indicate the type of the item we're printing
 * and to allow comparison of strings.
 */
static void io_update_string(Tcl_Obj *obj) {
    GapIO *io = obj->internalRep.otherValuePtr;
    obj->bytes = ckalloc(30);
    obj->length = sprintf(obj->bytes, "%s", io_obj_as_string(io));
}

/*
 * If we do things like:
 * "set io [g5::open_database -name foo]; puts [$io database]"
 * then $io is now a cmdName objType instead of gapio as that's the last
 * context it was used in.
 *
 * Here we provide the necessary mechanism to convert back from string form
 * to gapio type again.
 */
int io_from_any(Tcl_Interp *interp, Tcl_Obj *obj) {
    char *bytes;
    int length;
    GapIO *io;

    if (NULL == (bytes = Tcl_GetStringFromObj(obj, &length)))
	return TCL_ERROR;

    if (0 != strncmp(bytes, "io=", 3))
	return TCL_ERROR;

    /* Free the old internalRep before setting the new one. */
    if (obj->typePtr && obj->typePtr->freeIntRepProc)
	(*obj->typePtr->freeIntRepProc)(obj);

    /* Convert the hex value to a pointer once more */
    if (1 != sscanf(bytes+3, "%p", &io))
	return TCL_ERROR;

    obj->internalRep.otherValuePtr = io;
    obj->typePtr = &io_obj_type;
    return TCL_OK;
}

/* Returns a GapIO from any GapIO convertable object */
GapIO *io_from_obj(Tcl_Obj *obj) {
    if (obj->typePtr != &io_obj_type)
	if (TCL_ERROR == io_from_any(NULL /* unused */, obj))
	    return NULL;

    return (GapIO *)obj->internalRep.otherValuePtr;
}

static int io_cmd(ClientData clientData, Tcl_Interp *interp,
		  int objc, Tcl_Obj *CONST objv[]) {
    int index;
    GapIO *io = (GapIO *)clientData;

    static char *options[] = {
	"flush",       "close",	       "debug_level",
	"get_contig",  "get_sequence", "get_database", "get_anno_ele",
	"contig_order","num_contigs",  "seq_name2rec", "child",
	"get_library", "db_version",   "name",         "read_only",
	"new_contig",  "new_sequence", "new_anno_ele", "rec_exists",
	"seq_name_iter","seq_name_next","seq_name_end","check",
	"contig_name2rec", "base",
	(char *)NULL,
    };

    enum options {
	IO_FLUSH,     IO_CLOSE,       IO_DEBUG_LEVEL,
	IO_CONTIG,    IO_SEQUENCE,    IO_DATABASE,    IO_ANNO_ELE,
	IO_CORDER,    NUM_CONTIGS,    SEQ_NAME2REC,   IO_CHILD,
	IO_LIBRARY,   IO_DB_VERSION,  IO_NAME,        IO_READ_ONLY,
	NEW_CONTIG,   NEW_SEQUENCE,   NEW_ANNO_ELE,   IO_REC_EXISTS,
	SEQ_NAME_ITER,SEQ_NAME_NEXT,  SEQ_NAME_END,   CHECK,
	CONTIG_NAME2REC, IO_BASE
    };

    if (objc < 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "option arg ?arg ...?");
	return TCL_ERROR;
    }

    if (Tcl_GetIndexFromObj(interp, objv[1], options, "option", 0,
            &index) != TCL_OK) {
        return TCL_ERROR;
    }

    switch ((enum options)index) {
    case IO_CLOSE:
	gio_close(io);
	break;

    case IO_DEBUG_LEVEL: {
	int level;
	if (objc == 3) {
	    Tcl_GetIntFromObj(interp, objv[2], &level);
	    gio_debug_level(io, level);
	}
	break;
    }

    case IO_FLUSH:
	cache_flush(io);
	break;

    case IO_CONTIG:
	return tcl_contig_read(io, interp, objc-1, objv+1);
	break;

    case IO_SEQUENCE:
	return tcl_sequence_read(io, interp, objc-1, objv+1);
	break;

    case IO_DATABASE:
	return tcl_database_read(io, interp, objc-1, objv+1);

    case IO_ANNO_ELE:
	return tcl_anno_ele_read(io, interp, objc-1, objv+1);
	break;

    case IO_LIBRARY:
	return tcl_library_read(io, interp, objc-1, objv+1);
	break;

    case IO_CORDER:
	return tcl_contig_order(io, interp, objc-1, objv+1);
	break;

    case NUM_CONTIGS:
	vTcl_SetResult(interp, "%d", io->db->Ncontigs);
	break;

    case IO_DB_VERSION:
	vTcl_SetResult(interp, "%d", io->db->version);
	break;

    case IO_NAME:
	vTcl_SetResult(interp, "%s", io->name);
	break;

    case IO_READ_ONLY:
	vTcl_SetResult(interp, "%d", io->read_only);
	break;

    case SEQ_NAME2REC: {
	char *txt = Tcl_GetStringFromObj(objv[2], NULL);
	vTcl_SetResult(interp, "%"PRIrec, sequence_index_query(io, txt));
	break;
    }

    case CONTIG_NAME2REC: {
	char *txt = Tcl_GetStringFromObj(objv[2], NULL);
	vTcl_SetResult(interp, "%"PRIrec, contig_index_query(io, txt));
	break;
    }

    case NEW_CONTIG: {
	contig_t *c = contig_new(io, "contig");
	vTcl_SetResult(interp, "%"PRIrec, c->rec);
	break;
    }

    case NEW_SEQUENCE: {
	seq_t s;
	memset(&s, 0, sizeof(s));
	vTcl_SetResult(interp, "%"PRIrec, cache_item_create(io, GT_Seq, &s));
	break;
    }

    case NEW_ANNO_ELE: {
	int obj_type, start, end;
	Tcl_WideInt obj_rec;
	char dir = ANNO_DIR_NUL;

	if (objc != 6 && objc != 7) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s obj_type obj_rec start end ?dir?\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &obj_type);
	Tcl_GetWideIntFromObj(interp, objv[3], &obj_rec);
	Tcl_GetIntFromObj(interp, objv[4], &start);
	Tcl_GetIntFromObj(interp, objv[5], &end);
	if (objc == 7)
	    dir = *Tcl_GetStringFromObj(objv[6], NULL);

	vTcl_SetResult(interp, "%"PRIrec,
		       anno_ele_add(io, obj_type, obj_rec,
				    0 /* anno_rec */,
				    str2type("COMM"), "",
				    start, end, dir));
	break;
    }

    case IO_CHILD: {
	GapIO *child;
	Tcl_Obj *iobj;

	if (!(child = gio_child(io)))
	    return TCL_ERROR;
	if (NULL == (iobj = Tcl_NewObj()))
	    return TCL_ERROR;

	iobj->internalRep.otherValuePtr = child;
	iobj->typePtr = &io_obj_type;
	io_update_string(iobj);
	
	/* Register the string form as a new command */
	if (NULL == Tcl_CreateObjCommand(interp, iobj->bytes, io_cmd,
					 (ClientData)child,
					 (Tcl_CmdDeleteProc *)NULL))
	    return TCL_ERROR;
	
	Tcl_SetObjResult(interp, iobj);
	break;
    }

    case IO_BASE: {
	if (!io->base) {
	    Tcl_SetObjResult(interp, objv[0]);
	} else {
	    char *cmd = io_obj_as_string(io->base);
	    /* This *will* exist already in the Tcl Interpreter */
	    vTcl_SetResult(interp, "%s", cmd);
	}
	break;
    }

    case IO_REC_EXISTS: {
	int obj_type;
	Tcl_WideInt obj_rec;

	if (objc != 4) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s obj_type obj_rec\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &obj_type);
	Tcl_GetWideIntFromObj(interp, objv[3], &obj_rec);

	vTcl_SetResult(interp, "%d", cache_exists(io, obj_type, obj_rec));
	break;
    }
	
    case SEQ_NAME_ITER: {
	char *name = Tcl_GetStringFromObj(objv[2], NULL);
	io->seq_name_iter = sequence_index_iter(io, name);
	break;
    }

    case SEQ_NAME_NEXT: {
	BTRec rec;
	char *key;
	
	if (io->seq_name_iter && (key = btree_next(io->seq_name_iter, &rec))) {
	    Tcl_Obj *objv[2];
	    objv[0] = Tcl_NewStringObj(key, -1);
	    objv[1] = Tcl_NewLongObj(rec);
	    Tcl_SetObjResult(interp, Tcl_NewListObj(2, objv));
	} else {
	    Tcl_ResetResult(interp);
	}

	break;
    }

    case SEQ_NAME_END:
	if (io->seq_name_iter) {
	    btree_iter_del(io->seq_name_iter);
	    io->seq_name_iter = NULL;
	}
	break;

    case CHECK:  {
	int fix = 0, level = 2;
	if (objc >= 3)
	    Tcl_GetIntFromObj(interp, objv[2], &fix);
	if (objc >= 4)
	    Tcl_GetIntFromObj(interp, objv[3], &level);
	

	Tcl_SetObjResult(interp,
			 Tcl_NewIntObj(check_database(io, fix, level)));
	break;
    }
    }

    return TCL_OK;
}

/* ------------------------------------------------------------------------ */
/* Tcl_Obj "contig" type implementation */

static void contig_update_string(Tcl_Obj *obj);
static int contig_from_any(Tcl_Interp *interp, Tcl_Obj *obj);

typedef struct {
    GapIO *io;
    contig_t *contig;
    btree_iter_t *iter;
} tcl_contig;

static Tcl_ObjType contig_obj_type = {
    "contig",
    (Tcl_FreeInternalRepProc*)NULL,
    (Tcl_DupInternalRepProc*)NULL,
    contig_update_string,
    contig_from_any
};

static void contig_update_string(Tcl_Obj *obj) {
    tcl_contig *c = obj->internalRep.otherValuePtr;
    obj->bytes = ckalloc(30);
    obj->length = sprintf(obj->bytes, "contig=%p", c);
}

static int contig_from_any(Tcl_Interp *interp, Tcl_Obj *obj) {
    char *bytes;
    int length;
    tcl_contig *c;

    if (NULL == (bytes = Tcl_GetStringFromObj(obj, &length)))
	return TCL_ERROR;

    if (0 != strncmp(bytes, "contig=", 3))
	return TCL_ERROR;

    /* Free the old internalRep before setting the new one. */
    if (obj->typePtr && obj->typePtr->freeIntRepProc)
	(*obj->typePtr->freeIntRepProc)(obj);

    /* Convert the hex value to a pointer once more */
    if (1 != sscanf(bytes+3, "%p", &c))
	return TCL_ERROR;

    obj->internalRep.otherValuePtr = c;
    obj->typePtr = &contig_obj_type;
    return TCL_OK;
}

static int sort_range(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    return r1->start - r2->start;
}

static int tcl_contig_seqs_range(tcl_contig *tc, Tcl_Interp *interp,
				 int objc, Tcl_Obj *CONST objv[]) {
    GapIO *io = tc->io;
    contig_t *c = tc->contig;
    int start, end;
    Tcl_Obj *items;
    rangec_t *r;
    int nr, i;

    if (objc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s start end\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetIntFromObj(interp, objv[1], &start);
    Tcl_GetIntFromObj(interp, objv[2], &end);
    r = (rangec_t *)contig_seqs_in_range(io, &c, start, end,
					 CSIR_SORT_BY_X | CSIR_PAIR,
					 &nr);
    qsort(r, nr, sizeof(*r), sort_range);

    items = Tcl_NewListObj(0, NULL);
    for (i = 0; i < nr; i++) {
	Tcl_Obj *ele, *e4[15];

	e4[0]  = Tcl_NewIntObj(r[i].start);
	e4[1]  = Tcl_NewIntObj(r[i].end);
	e4[2]  = Tcl_NewWideIntObj(r[i].rec);
	e4[3]  = Tcl_NewIntObj(r[i].mqual);
	e4[4]  = Tcl_NewIntObj((r[i].flags & GRANGE_FLAG_COMP1) ? 1 : 0);
	e4[5]  = Tcl_NewIntObj((r[i].flags & GRANGE_FLAG_END_MASK) ? 1 : 0);
	e4[6]  = Tcl_NewIntObj(r[i].pair_start);
	e4[7]  = Tcl_NewIntObj(r[i].pair_end);
	e4[8]  = Tcl_NewWideIntObj(r[i].pair_rec);
	e4[9]  = Tcl_NewIntObj(r[i].pair_mqual);
	e4[10] = Tcl_NewIntObj((r[i].flags & GRANGE_FLAG_COMP2) ? 1 : 0);
	e4[11] = Tcl_NewIntObj((r[i].flags & GRANGE_FLAG_PEND_MASK) ? 1 : 0);
	e4[12] = Tcl_NewIntObj(r[i].flags & GRANGE_FLAG_TYPE_MASK);
	e4[13] = Tcl_NewIntObj((r[i].flags & GRANGE_FLAG_CONTIG) ? 1 : 0);
	e4[14] = Tcl_NewIntObj(r[i].comp);
	ele = Tcl_NewListObj(15, e4);

	Tcl_ListObjAppendElement(interp, items, ele);
    }

    Tcl_SetObjResult(interp, items);

    free(r);

    return TCL_OK;
}

static int tcl_contig_anno_range(tcl_contig *tc, Tcl_Interp *interp,
				 int objc, Tcl_Obj *CONST objv[]) {
    GapIO *io = tc->io;
    contig_t *c = tc->contig;
    int start, end;
    Tcl_Obj *items;
    rangec_t *r;
    int nr, i;

    if (objc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s start end\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetIntFromObj(interp, objv[1], &start);
    Tcl_GetIntFromObj(interp, objv[2], &end);
    r = (rangec_t *)contig_anno_in_range(io, &c, start, end,
					 CSIR_SORT_BY_X,
					 &nr);

    items = Tcl_NewListObj(0, NULL);
    for (i = 0; i < nr; i++) {
	Tcl_Obj *ele, *e4[15];

	e4[0]  = Tcl_NewIntObj(r[i].start);
	e4[1]  = Tcl_NewIntObj(r[i].end);
	e4[2]  = Tcl_NewWideIntObj(r[i].rec);
	e4[3]  = Tcl_NewIntObj(r[i].mqual);
	e4[4]  = Tcl_NewIntObj((r[i].flags & GRANGE_FLAG_COMP1) ? 1 : 0);
	e4[5]  = Tcl_NewIntObj((r[i].flags & GRANGE_FLAG_END_MASK) ? 1 : 0);
	e4[6]  = Tcl_NewIntObj(r[i].pair_start);
	e4[7]  = Tcl_NewWideIntObj(r[i].pair_end);
	e4[8]  = Tcl_NewWideIntObj(r[i].pair_rec);
	e4[9]  = Tcl_NewIntObj(r[i].pair_mqual);
	e4[10] = Tcl_NewIntObj((r[i].flags & GRANGE_FLAG_COMP2) ? 1 : 0);
	e4[11] = Tcl_NewIntObj((r[i].flags & GRANGE_FLAG_PEND_MASK) ? 1 : 0);
	e4[12] = Tcl_NewIntObj(r[i].flags & GRANGE_FLAG_TYPE_MASK);
	e4[13] = Tcl_NewIntObj((r[i].flags & GRANGE_FLAG_CONTIG) ? 1 : 0);
	e4[14] = Tcl_NewIntObj(r[i].comp);
	ele = Tcl_NewListObj(15, e4);

	Tcl_ListObjAppendElement(interp, items, ele);
    }

    Tcl_SetObjResult(interp, items);

    free(r);

    return TCL_OK;
}

static int tcl_contig_pileup(tcl_contig *tc, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]) {
    GapIO *io = tc->io;
    contig_t *c = tc->contig;
    int pos;
    Tcl_Obj *items;
    rangec_t *r;
    int nr, i;

    if (objc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s pos\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    /* Fetch all sequences covering this point */
    Tcl_GetIntFromObj(interp, objv[1], &pos);
    r = (rangec_t *)contig_seqs_in_range(io, &c, pos, pos,
					 CSIR_SORT_BY_X | CSIR_PAIR,
					 &nr);
    qsort(r, nr, sizeof(*r), sort_range);

    /* Produce a tcl list of elements consisting of seq rec, pos, base, qual */
    items = Tcl_NewListObj(0, NULL);
    for (i = 0; i < nr; i++) {
	Tcl_Obj *ele, *e6[6];
	seq_t *s = cache_search(io, GT_Seq, r[i].rec);
	char base;
	int conf, ret, cut;

	ret = sequence_get_base(io, &s, pos - r[i].start, &base, &conf,
				&cut, 1);
	if (-1 == ret) {
	    base = '?';
	    conf = 1;
	    fprintf(stderr, "ERROR: failed to read base at position %d "
		    "in seq #%"PRIrec"\n", pos, r[i].rec);
	}
	e6[0] = Tcl_NewWideIntObj(r[i].rec);
	e6[1] = Tcl_NewIntObj(pos - r[i].start);
	e6[2] = Tcl_NewStringObj(&base, 1);
	e6[3] = Tcl_NewIntObj(conf);
	e6[4] = Tcl_NewIntObj(cut);
	e6[5] = Tcl_NewIntObj(r[i].start);

	ele = Tcl_NewListObj(6, e6);

	Tcl_ListObjAppendElement(interp, items, ele);
    }

    Tcl_SetObjResult(interp, items);

    free(r);
    return TCL_OK;
}

static int tcl_read_depth(tcl_contig *tc, Tcl_Interp *interp,
			  int objc, Tcl_Obj *CONST objv[]) {
    GapIO *io = tc->io;
    contig_t *c = tc->contig;
    int start, end, i;
    double bpv;
    Tcl_Obj *items;
    track_t *track;

    if (objc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s start end bpv\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetIntFromObj(interp, objv[1], &start);
    Tcl_GetIntFromObj(interp, objv[2], &end);
    Tcl_GetDoubleFromObj(interp, objv[3], &bpv);

    track = contig_get_track(io, &c, start, end, TRACK_READ_DEPTH, bpv);

    items = Tcl_NewListObj(0, NULL);
    for (i = 0; i < track->nitems; i++) {
	int d = arr(int, track->data, i);
	Tcl_ListObjAppendElement(interp, items, Tcl_NewIntObj(d));
    }

    track_free(track);

    Tcl_SetObjResult(interp, items);

    return TCL_OK;
}

static int contig_cmd(ClientData clientData, Tcl_Interp *interp,
		      int objc, Tcl_Obj *CONST objv[]) {
    int index;
    tcl_contig *tc = (tcl_contig *)clientData;

    static char *options[] = {
	"delete",       "io",           "dump_ps",
	"get_start",    "get_end",      "get_len",      "get_length",
	"get_name",     "seqs_in_range","get_rec",      "read_depth",
	"insert_base",  "delete_base",  "remove_sequence","add_sequence",
	"nseqs",	"anno_in_range","get_pileup",   "ref_to_padded",
	"nrefpos",	"nanno",        "shift_base",   "move_anno",
	"check",        "move_seq",
	"get_visible_start", "get_visible_end", "get_visible_length",
	"set_visible_start", "invalidate_consensus",    "set_name",
	(char *)NULL,
    };

    enum options {
	DELETE,         IO,     	DUMP_PS,
	GET_START,      GET_END,        GET_LEN,        GET_LENGTH,
	GET_NAME,       SEQS_IN_RANGE,  GET_REC,        READ_DEPTH,
	INSERT_BASE,    DELETE_BASE,    REMOVE_SEQUENCE,ADD_SEQUENCE,
	NSEQS,          ANNO_IN_RANGE,  GET_PILEUP,     REF_TO_PADDED,
	NREFPOS,        NANNO,	        SHIFT_BASE,     MOVE_ANNO,
	CHECK,          MOVE_SEQ,
	GET_VISIBLE_START, GET_VISIBLE_END, GET_VISIBLE_LENGTH,
	SET_VISIBLE_START, INVALIDATE_CONSENSUS,        SET_NAME
    };

    if (objc < 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "option arg ?arg ...?");
	return TCL_ERROR;
    }

    if (Tcl_GetIndexFromObj(interp, objv[1], options, "option", 0,
            &index) != TCL_OK) {
        return TCL_ERROR;
    }

    switch ((enum options)index) {
    case DELETE:
	Tcl_DeleteCommandFromToken(interp,
				   Tcl_GetCommandFromObj(interp, objv[0]));

	break;

    case IO:
	Tcl_SetResult(interp, io_obj_as_string(tc->io) , TCL_VOLATILE);
	break;

    case DUMP_PS:
	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s dump_ps filename\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}
	contig_dump_ps(tc->io, &tc->contig,
		       Tcl_GetStringFromObj(objv[2], NULL));
	Tcl_ResetResult(interp);
	break;

    case GET_REC:
	Tcl_SetWideIntObj(Tcl_GetObjResult(interp), ci_ptr(tc->contig)->rec);
	break;

    case GET_START:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), tc->contig->start);
	break;

    case GET_END:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), tc->contig->end);
	break;

    case GET_LEN:
    case GET_LENGTH:
	Tcl_SetIntObj(Tcl_GetObjResult(interp),
		      tc->contig->end - tc->contig->start + 1);
	break;

    case GET_VISIBLE_START: {
	int st;
	if (-1 == consensus_valid_range(tc->io, tc->contig->rec, &st, NULL))
	    return TCL_ERROR;
	Tcl_SetIntObj(Tcl_GetObjResult(interp), st);
	break;
    }

    case GET_VISIBLE_END: {
	int en;
	if (-1 == consensus_valid_range(tc->io, tc->contig->rec, NULL, &en))
	    return TCL_ERROR;
	Tcl_SetIntObj(Tcl_GetObjResult(interp), en);
	break;
    }

    case GET_VISIBLE_LENGTH: {
	int st, en;
	if (-1 == consensus_valid_range(tc->io, tc->contig->rec, &st, &en))
	    return TCL_ERROR;
	Tcl_SetIntObj(Tcl_GetObjResult(interp), en-st+1);
	break;
    }

    case SET_VISIBLE_START: {
	int pos, r;
	tg_rec crec;

	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_visible_start pos\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}
	Tcl_GetIntFromObj(interp, objv[2], &pos);
	
	crec = tc->contig->rec;
	cache_decr(tc->io, tc->contig);

	r = contig_set_visible_start(tc->io, crec, pos);
	Tcl_SetIntObj(Tcl_GetObjResult(interp), r);

	tc->contig = cache_search(tc->io, GT_Contig, crec);
	cache_incr(tc->io, tc->contig);

	break;
    }

    case GET_NAME:
	Tcl_SetStringObj(Tcl_GetObjResult(interp), tc->contig->name, -1);
	break;

    case SET_NAME: {
	char *name;
	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_name new_name\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}
	name = Tcl_GetStringFromObj(objv[2], NULL);

	Tcl_SetIntObj(Tcl_GetObjResult(interp),
		      contig_set_name(tc->io, &tc->contig, name));
	break;
    }

    case SEQS_IN_RANGE:
	return tcl_contig_seqs_range(tc, interp, objc-1, objv+1);

    case ANNO_IN_RANGE:
	return tcl_contig_anno_range(tc, interp, objc-1, objv+1);

    case READ_DEPTH:
	return tcl_read_depth(tc, interp, objc-1, objv+1);

    case INSERT_BASE: {
	int pos, qual = -1;
	char base = '*';
	if (objc < 3 || objc > 5) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s insert_base position ?base ?qual??\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &pos);
	if (objc >= 4) {
	    base = *Tcl_GetStringFromObj(objv[3], NULL);
	}
	if (objc == 5) {
	    Tcl_GetIntFromObj(interp, objv[4], &qual);
	}
	contig_insert_base(tc->io, &tc->contig, pos, base, qual);
	break;
    }

    case DELETE_BASE: {
	int pos;
	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s delete_base position\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &pos);
	contig_delete_base(tc->io, &tc->contig, pos);
	break;
    }

    case SHIFT_BASE: {
	int pos, dir;
	if (objc != 4) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s insert_base position dir(+1/-1)\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &pos);
	Tcl_GetIntFromObj(interp, objv[3], &dir);
	vTcl_SetResult(interp, "%d",
		       contig_shift_base(tc->io, &tc->contig, pos, dir));
	break;
    }

    case REMOVE_SEQUENCE: {
	Tcl_WideInt rec;
	seq_t *s;
	bin_index_t *b;
	range_t *r;

	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s remove_sequence rec\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetWideIntFromObj(interp, objv[2], &rec);

	/* Get old range and pair data */
	s = cache_search(tc->io, GT_Seq, rec);
	b = cache_search(tc->io, GT_Bin, s->bin);
	r = arrp(range_t, b->rng, s->bin_index);
	assert(r->rec == s->rec);
	assert(ABS(r->end - r->start) + 1 == ABS(s->len));

	vTcl_SetResult(interp, "%"PRIrec" %d", r->pair_rec, r->flags);

	bin_remove_item(tc->io, &tc->contig, GT_Seq, rec);
	break;
    }

    case ADD_SEQUENCE: {
	Tcl_WideInt rec;
	int pos;
	range_t r, *r_out;
	seq_t *s;
	bin_index_t *bin;
	Tcl_WideInt pair_rec;
	int flags;

	if (objc < 4 || objc > 6) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s add_sequence rec pos ?pair_rec ?flags??\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetWideIntFromObj(interp, objv[2], &rec);
	Tcl_GetIntFromObj(interp, objv[3], &pos);

	s = (seq_t *)cache_search(tc->io, GT_Seq, rec);
	memset(&r, 0, sizeof(r));
	r.start = pos;
	r.end   = pos + (s->len > 0 ? s->len : -s->len) - 1;
	r.rec   = rec;
	r.mqual = s->mapping_qual;

	if (objc >= 5) {
	    Tcl_GetWideIntFromObj(interp, objv[4], &pair_rec);
	    r.pair_rec = pair_rec;
	} else {
	    /* Insufficient in most cases, but we can override it */
	    if (s->parent_type == GT_Seq)
		r.pair_rec = s->parent_rec;
	    else
		r.pair_rec = 0;
	}

	
	/* What about other flags? Can't guess */
	if (objc >= 6) {
	    Tcl_GetIntFromObj(interp, objv[5], &flags);
	    r.flags = flags;
	} else {
	    r.flags = 0;
	    if (s->flags & SEQ_END_REV)
		r.flags |= GRANGE_FLAG_END_REV;
	    if (s->flags & SEQ_END_FWD)
		r.flags |= GRANGE_FLAG_END_FWD;
	    if (s->len < 0)
		r.flags |= GRANGE_FLAG_COMP1;
	}

	bin = bin_add_range(tc->io, &tc->contig, &r, &r_out, NULL, 0);
	if (s->bin != bin->rec) {
	    int old_comp = bin_get_orient(tc->io, s->bin);
	    int new_comp = bin_get_orient(tc->io, bin->rec);

	    //printf("New seq bin (%d)%d->(%d)%d\n",
	    //	   old_comp, s->bin, new_comp, bin->rec);

	    /* Bin number changed - update seq too */
	    s = cache_rw(tc->io, s);
	    s->bin = bin->rec;
	    s->bin_index = r_out - ArrayBase(range_t, bin->rng);

	    /* Check if the new bin has a different complemented status too */
	    if (new_comp != old_comp) {
		s->len *= -1;
		s->flags ^= SEQ_COMPLEMENTED;
		//tmp = s->left;
		//s->left  = ABS(s->len) - (s->right-1);
		//s->right = ABS(s->len) - (tmp-1);
	    }
	}
	break;
    }

    case MOVE_SEQ: {
	/* A combination of remove and add */
	Tcl_WideInt rec;
	seq_t *s;
	bin_index_t *bin;
	range_t r, *r_out;
	int dist, dir;

	if (objc < 4) {
	    vTcl_SetResult(interp, "wrong # args: should be \"%s "
			   "move_seq rec distance\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetWideIntFromObj(interp, objv[2], &rec);
	Tcl_GetIntFromObj(interp, objv[3], &dist);

	/* Get old range coords and convert from relative to absolute */
	s = cache_search(tc->io, GT_Seq, rec);
	cache_incr(tc->io, s);

	bin = cache_search(tc->io, GT_Bin, s->bin);
	r = arr(range_t, bin->rng, s->bin_index);
	assert(r.rec == s->rec);
	assert(ABS(r.end - r.start) + 1 == ABS(s->len));
	sequence_get_position(tc->io, s->rec, NULL, &r.start, &r.end, &dir);

	bin_remove_item(tc->io, &tc->contig, GT_Seq, rec);

	/* Add it back at the new range */
	r.start += dist;
	r.end += dist;
	bin = bin_add_range(tc->io, &tc->contig, &r, &r_out, NULL, 0);

	/* Update seq if parent has changed */
	if (s->bin != bin->rec) {
	    int old_comp = bin_get_orient(tc->io, s->bin);
	    int new_comp = bin_get_orient(tc->io, bin->rec);

	    s = cache_rw(tc->io, s);
	    s->bin = bin->rec;
	    s->bin_index = r_out - ArrayBase(range_t, bin->rng);

	    /* Check if the new bin has a different complemented status too */
	    if (new_comp != old_comp) {
		s->len *= -1;
		s->flags ^= SEQ_COMPLEMENTED;
		//tmp = s->left;
		//s->left  = ABS(s->len) - (s->right-1);
		//s->right = ABS(s->len) - (tmp-1);
	    }
	}

	cache_decr(tc->io, s);

	break;
    }

    case MOVE_ANNO: {
	Tcl_WideInt rec;
	int start, end;
	range_t r, *r_out;
	anno_ele_t *a;
	bin_index_t *bin;
	Tcl_WideInt obj_rec;
	int obj_type;
	tg_rec seq_bin = 0;

	/* Parse args */
	if (objc < 4 || objc > 7) {
	    vTcl_SetResult(interp, "wrong # args: should be \"%s "
			   "move_anno rec start end ?obj_type obj_rec?\" or "
			   "\"%s move_anno rec distance\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL),
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetWideIntFromObj(interp, objv[2], &rec);
	a = (anno_ele_t *)cache_search(tc->io, GT_AnnoEle, rec);

	if (objc == 4) {
	    int dist;

	    Tcl_GetIntFromObj(interp, objv[3], &dist);
	    anno_get_position(tc->io, rec, NULL, &start, &end, NULL);
	    start += dist;
	    end   += dist;
	} else {
	    Tcl_GetIntFromObj(interp, objv[3], &start);
	    Tcl_GetIntFromObj(interp, objv[4], &end);
	}

	if (objc >= 5) {
	    Tcl_GetIntFromObj(interp, objv[5], &obj_type);
	} else {
	    obj_type = a->obj_type;
	}

	if (objc >= 6) {
	    Tcl_GetWideIntFromObj(interp, objv[6], &obj_rec);
	} else {
	    if (obj_type == GT_Contig) {
		obj_rec = tc->contig->rec;
	    } else {
		obj_rec = a->obj_rec;
	    }
	}


	/* Remove from old place */
	bin_remove_item(tc->io, &tc->contig, GT_AnnoEle, rec);

	if (obj_type == GT_Seq) {
	    int st, en;
	    cache_incr(tc->io, a);
	    sequence_get_position2(tc->io, obj_rec, NULL, &st, &en, NULL,
				   &seq_bin, NULL, NULL);
	    cache_decr(tc->io, a);

	    start += st;
	    end += st;
	}

	/* Add back to new location */
	memset(&r, 0, sizeof(r));
	r.start    = start;
	r.end      = end;
	r.rec      = rec;
	r.mqual    = a->tag_type;
	r.pair_rec = obj_rec;
	r.flags    = GRANGE_FLAG_ISANNO;
	if (GT_Seq == obj_type)
	    r.flags |= GRANGE_FLAG_TAG_SEQ;

	if (seq_bin)
	    bin = bin_add_to_range(tc->io, &tc->contig, seq_bin,
				   &r, &r_out, NULL, 0);
	else
	    bin = bin_add_range(tc->io, &tc->contig, &r, &r_out, NULL, 0);


	/* The move may have changed bin, if so update anno pointer too */
	if (a->bin != bin->rec) {
	    a = cache_rw(tc->io, a);
	    a->bin = bin->rec;
	}
	break;
    }

    case NSEQS: {
	int nseqs;
	if (!tc->contig->bin) {
	    nseqs = 0;
	} else {
	    bin_index_t *bin;
	    bin = (bin_index_t *)cache_search(tc->io, GT_Bin, tc->contig->bin);
	    nseqs = bin->nseqs;
	}

	Tcl_SetObjResult(interp, Tcl_NewIntObj(nseqs));
	break;
    } 

    case NREFPOS: {
	int nrefpos;
	if (!tc->contig->bin) {
	    nrefpos = 0;
	} else {
	    bin_index_t *bin;
	    bin = (bin_index_t *)cache_search(tc->io, GT_Bin, tc->contig->bin);
	    nrefpos = bin->nrefpos;
	}

	Tcl_SetObjResult(interp, Tcl_NewIntObj(nrefpos));
	break;
    } 

    case NANNO: {
	int nanno;
	if (!tc->contig->bin) {
	    nanno = 0;
	} else {
	    bin_index_t *bin;
	    bin = (bin_index_t *)cache_search(tc->io, GT_Bin, tc->contig->bin);
	    nanno = bin->nanno;
	}

	Tcl_SetObjResult(interp, Tcl_NewIntObj(nanno));
	break;
    } 

    case GET_PILEUP: 
	return tcl_contig_pileup(tc, interp, objc-1, objv+1);
	break;
	
    case REF_TO_PADDED: {
	int rpos, ppos;

	if (objc < 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s ref_to_padded ref_pos ?padded_start_pos?\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &rpos);
	if (objc >= 4) {
	    /* From known starting point */
	    Tcl_GetIntFromObj(interp, objv[3], &ppos);
	    if (!reference_to_padded_pos2(tc->io, tc->contig->rec, -1,
					  rpos, ppos, &ppos))
		Tcl_SetIntObj(Tcl_GetObjResult(interp), ppos);
	    else
		return TCL_ERROR;
	} else {
	    if (!reference_to_padded_pos(tc->io, tc->contig->rec, -1,
					 rpos, &ppos))
		Tcl_SetIntObj(Tcl_GetObjResult(interp), ppos);
	    else
		return TCL_ERROR;
	}

	break;
    }

    case CHECK: {
	int fix = 0, level = 2, fixed = 0;
	int ret;

	if (objc >= 3)
	    Tcl_GetIntFromObj(interp, objv[2], &fix);
	if (objc >= 4)
	    Tcl_GetIntFromObj(interp, objv[3], &level);
	
	ret = check_contig(tc->io, tc->contig->rec, fix, level, NULL,
			   &fixed, NULL);
	vTcl_SetResult(interp, "%d %d", ret, fixed);
	break;
    }

    case INVALIDATE_CONSENSUS: {
	int start = tc->contig->start;
	int end   = tc->contig->end;
	int ret;

	if (objc >= 3)
	    Tcl_GetIntFromObj(interp, objv[2], &start);
	if (objc >= 4)
	    Tcl_GetIntFromObj(interp, objv[3], &end);

	ret = bin_invalidate_consensus(tc->io, tc->contig->rec, start, end);
	vTcl_SetResult(interp, "%d", ret);
    }
    }

    return TCL_OK;
}

static void _cmd_delete(ClientData clientData) {
    /* Could be any of the types in this file as they're compatible structs */
    tcl_contig *tc = (tcl_contig *)clientData;
    if (tc->contig)
	cache_decr(tc->io, tc->contig);
    if (tc->iter)
	btree_iter_del(tc->iter);
}

static int tcl_contig_read(GapIO *io, Tcl_Interp *interp,
			   int objc, Tcl_Obj *CONST objv[]) {
    contig_t *c;
    Tcl_WideInt cnum;
    tcl_contig *tc;
    Tcl_Obj *res;

    if (objc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s contig_id\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetWideIntFromObj(interp, objv[1], &cnum);
    if (cnum > 0)
	c = (contig_t *)cache_search(io, GT_Contig, cnum);
    else
	gio_read_contig(io, -(int)cnum, &c); /* contig order lookup */

    if (!c) {
	vTcl_SetResult(interp, "Unable to read contig =%"PRIrec, cnum);
	return TCL_ERROR;
    }

    if (NULL == (tc = (tcl_contig *)ckalloc(sizeof(*tc))))
	return TCL_ERROR;
    tc->io = io;
    tc->contig = c;
    tc->iter = NULL;

    if (NULL == (res = Tcl_NewObj()))
	return TCL_ERROR;

    res->internalRep.otherValuePtr = tc;
    res->typePtr = &contig_obj_type;
    contig_update_string(res);

    /* Register the string form as a new command */
    cache_incr(io, c);
    if (NULL == Tcl_CreateObjCommand(interp, res->bytes, contig_cmd,
				     (ClientData)tc,
				     (Tcl_CmdDeleteProc *)_cmd_delete))
	return TCL_ERROR;

    Tcl_SetObjResult(interp, res);

    return TCL_OK;
}

static int tcl_contig_order(GapIO *io, Tcl_Interp *interp,
			    int objc, Tcl_Obj *CONST objv[]) {
    int cind;
    tg_rec crec;

    if (objc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s contig_index\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetIntFromObj(interp, objv[1], &cind);
    crec = arr(tg_rec, io->contig_order, cind);

    vTcl_SetResult(interp, "%"PRIrec, crec);
    return TCL_OK;
}

/* ------------------------------------------------------------------------ */
/* Tcl_Obj "sequence" type implementation */

static void sequence_update_string(Tcl_Obj *obj);
static int sequence_from_any(Tcl_Interp *interp, Tcl_Obj *obj);

typedef struct {
    GapIO *io;
    seq_t *seq;
    btree_iter_t *iter;
} tcl_sequence;

static Tcl_ObjType sequence_obj_type = {
    "sequence",
    (Tcl_FreeInternalRepProc*)NULL,
    (Tcl_DupInternalRepProc*)NULL,
    sequence_update_string,
    sequence_from_any
};

static void sequence_update_string(Tcl_Obj *obj) {
    tcl_sequence *ts = obj->internalRep.otherValuePtr;
    obj->bytes = ckalloc(30);
    obj->length = sprintf(obj->bytes, "sequence=%p", ts);
}

static int sequence_from_any(Tcl_Interp *interp, Tcl_Obj *obj) {
    char *bytes;
    int length;
    tcl_sequence *ts;

    if (NULL == (bytes = Tcl_GetStringFromObj(obj, &length)))
	return TCL_ERROR;

    if (0 != strncmp(bytes, "sequence=", 3))
	return TCL_ERROR;

    /* Free the old internalRep before setting the new one. */
    if (obj->typePtr && obj->typePtr->freeIntRepProc)
	(*obj->typePtr->freeIntRepProc)(obj);

    /* Convert the hex value to a pointer once more */
    if (1 != sscanf(bytes+3, "%p", &ts))
	return TCL_ERROR;

    obj->internalRep.otherValuePtr = ts;
    obj->typePtr = &sequence_obj_type;
    return TCL_OK;
}

static int sequence_cmd(ClientData clientData, Tcl_Interp *interp,
			int objc, Tcl_Obj *CONST objv[]) {
    int index;
    tcl_sequence *ts = (tcl_sequence *)clientData;

    static char *options[] = {
	"delete",       "io",
	"get_rec",      "get_len",      "get_length",   "get_pair",
	"get_left",     "get_right",    "get_name",     "get_seq",
	"get_conf",	"get_conf4",    "get_contig",   "get_position",
	"get_clipped_position",         "get_orient",   "get_mapping_qual",
	"get_base",     "insert_base",  "delete_base",  "replace_base",
	"get_clips",    "set_clips",    "move_annos",   "get_template_orient",
	"set_clips_no_invalidate",
	(char *)NULL,
    };

    enum options {
	DELETE,         IO,
	GET_REC,        GET_LEN,        GET_LENGTH,     GET_PAIR,
	GET_LEFT,	GET_RIGHT,      GET_NAME,       GET_SEQ,
	GET_CONF,       GET_CONF4,      GET_CONTIG,     GET_POSITION,
	GET_CLIPPED_POSITION,           GET_ORIENT,     GET_MAPPING_QUAL,
	GET_BASE,       INSERT_BASE,    DELETE_BASE,    REPLACE_BASE,
	GET_CLIPS,      SET_CLIPS,      MOVE_ANNOS,     GET_TEMPLATE_ORIENT,
	SET_CLIPS_NO_INVALIDATE,
    };

    if (objc < 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "option arg ?arg ...?");
	return TCL_ERROR;
    }

    if (Tcl_GetIndexFromObj(interp, objv[1], options, "option", 0,
            &index) != TCL_OK) {
        return TCL_ERROR;
    }

    switch ((enum options)index) {
    case DELETE:
	Tcl_DeleteCommandFromToken(interp,
				   Tcl_GetCommandFromObj(interp, objv[0]));
	break;

    case IO:
	Tcl_SetResult(interp, io_obj_as_string(ts->io) , TCL_VOLATILE);
	break;

    case GET_REC:
	Tcl_SetWideIntObj(Tcl_GetObjResult(interp), ts->seq->rec);
	break;

    case GET_LEN:
    case GET_LENGTH:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), ts->seq->len);
	break;

    case GET_LEFT:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), ts->seq->left);
	break;

    case GET_RIGHT:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), ts->seq->right);
	break;

    case GET_NAME:
	Tcl_SetStringObj(Tcl_GetObjResult(interp),
			 ts->seq->name, ts->seq->name_len);
	break;

    case GET_SEQ:
	Tcl_SetStringObj(Tcl_GetObjResult(interp),
			 ts->seq->seq, ABS(ts->seq->len));
	break;

    case GET_CONF:
	if (ts->seq->format != SEQ_FORMAT_CNF4) {
	    Tcl_SetStringObj(Tcl_GetObjResult(interp),
			     ts->seq->conf, ABS(ts->seq->len));
	} else {
	    int len = ABS(ts->seq->len);
	    char *buf = malloc(len);
	    int i;
	    for (i = 0; i < len; i++) {
		switch(ts->seq->seq[i]) {
		case 'A': case 'a':
		    buf[i] = ts->seq->conf[i*4+0];
		    break;
		case 'C': case 'c':
		    buf[i] = ts->seq->conf[i*4+1];
		    break;
		case 'G': case 'g':
		    buf[i] = ts->seq->conf[i*4+2];
		    break;
		case 'T': case 't':
		    buf[i] = ts->seq->conf[i*4+3];
		    break;
		default:
		    buf[i] = -5;
		}
	    }
	    Tcl_SetStringObj(Tcl_GetObjResult(interp), buf, len);
	    free(buf);
	}
	break;

    case GET_CONF4:
	if (ts->seq->format == SEQ_FORMAT_CNF4) {
	    Tcl_SetStringObj(Tcl_GetObjResult(interp),
			     ts->seq->conf, ABS(ts->seq->len)*4);
	} else {
	    int len = ABS(ts->seq->len);
	    char *buf = malloc(len*4);
	    int i;
	    for (i = 0; i < len; i++) {
		/* Hack for now */
		switch(ts->seq->seq[i]) {
		case 'A': case 'a':
		    buf[i*4+0] = ts->seq->conf[i];
		    buf[i*4+1] = 0;
		    buf[i*4+2] = 0;
		    buf[i*4+3] = 0;
		    break;
		case 'C': case 'c':
		    buf[i*4+0] = 0;
		    buf[i*4+1] = ts->seq->conf[i];
		    buf[i*4+2] = 0;
		    buf[i*4+3] = 0;
		    break;
		case 'G': case 'g':
		    buf[i*4+0] = 0;
		    buf[i*4+1] = 0;
		    buf[i*4+2] = ts->seq->conf[i];
		    buf[i*4+3] = 0;
		    break;
		case 'T': case 't':
		    buf[i*4+0] = 0;
		    buf[i*4+1] = 0;
		    buf[i*4+2] = 0;
		    buf[i*4+3] = ts->seq->conf[i];
		    break;
		default:
		    buf[i*4+0] = -5;
		    buf[i*4+1] = -5;
		    buf[i*4+2] = -5;
		    buf[i*4+3] = -5;
		}
	    }
	    Tcl_SetStringObj(Tcl_GetObjResult(interp), buf, len*4);
	    free(buf);
	}
	break;

    case GET_CONTIG: {
	tg_rec rec = ts->seq->rec;
	tg_rec cnum = sequence_get_contig(ts->io, rec);
	Tcl_SetWideIntObj(Tcl_GetObjResult(interp), cnum);
	break;
    }

    case GET_POSITION: {
	tg_rec rec = ts->seq->rec, cnum;
	int pos;
	sequence_get_position(ts->io, rec, &cnum, &pos, NULL, NULL);
	Tcl_SetIntObj(Tcl_GetObjResult(interp), pos);
	break;
    }

    case GET_CLIPPED_POSITION: {
	tg_rec rec = ts->seq->rec, cnum;
	int pos;
	sequence_get_clipped_position(ts->io, rec, &cnum,
				      NULL, NULL, &pos, NULL, NULL);
	Tcl_SetIntObj(Tcl_GetObjResult(interp), pos);
	break;
    }

    case GET_ORIENT: {
	tg_rec rec = ts->seq->rec, cnum;
	int pos, dir;
	range_t r;
	seq_t *s;
	sequence_get_position2(ts->io, rec, &cnum, &pos, NULL, &dir, NULL,
			       &r, &s);
	if (!s) {
	    Tcl_SetResult(interp, "Unable to find sequence position.",
			  TCL_STATIC);
	    return TCL_ERROR;
	}

	Tcl_SetIntObj(Tcl_GetObjResult(interp), (s->len < 0) ^ dir);
	cache_decr(ts->io, s);
	break;
    }

    case GET_TEMPLATE_ORIENT: {
	/*
	 * Computes the orientation of the template. Useful for checking
	 * consistency.
	 *
	 * Eg for a standard --> <-- style read-pair:
	 * 0  -fwd-> <-rev-
	 * 1  -rev-> <-fwd-
	 *
	 * Hence -fwd-> -rev-> gives orient 0/1 and a conflict.
	 */
	tg_rec rec = ts->seq->rec, cnum;
	int pos, dir;
	range_t r;
	seq_t *s;
	int tdir;
	int lib_type = LIB_T_INWARD;

	sequence_get_position2(ts->io, rec, &cnum, &pos, NULL, &dir, NULL,
			       &r, &s);
	if (!s) {
	    Tcl_SetResult(interp, "Unable to find sequence position.",
			  TCL_STATIC);
	    return TCL_ERROR;
	}

	if (s->parent_type == GT_Library)
	    get_library_stats(ts->io, s->parent_rec,
			      NULL, NULL, &lib_type, NULL);

	switch (lib_type) {
	case LIB_T_INWARD:
	case LIB_T_OUTWARD:
	    tdir = (s->flags & SEQ_END_MASK) == SEQ_END_REV;
	    tdir ^= (s->len < 0) ^ dir;
	    break;

	case LIB_T_SAME:
	default:
	    tdir = (s->len < 0) ^ dir;
	}

	Tcl_SetIntObj(Tcl_GetObjResult(interp), tdir);
	cache_decr(ts->io, s);
	break;
    }

    case GET_MAPPING_QUAL:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), ts->seq->mapping_qual);
	break;

    case GET_PAIR:
	Tcl_SetWideIntObj(Tcl_GetObjResult(interp),
			  sequence_get_pair(ts->io, ts->seq));
	break;

    case GET_BASE: {
	char base;
	int pos, conf, cutoff;

	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s get_base position\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &pos);
	if (-1 == sequence_get_base(ts->io, &ts->seq, pos,
				    &base, &conf, &cutoff, 1))
	    return TCL_ERROR;

	vTcl_SetResult(interp, "%c %d %d", base, conf, cutoff);
	break;
    }

    case INSERT_BASE: {
	char base;
	int pos, conf;

	if (objc != 5) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s insert_base position call confidence\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &pos);
	base = *Tcl_GetString(objv[3]);
	Tcl_GetIntFromObj(interp, objv[4], &conf);
	sequence_insert_base(ts->io, &ts->seq, pos, base, conf, 1);
	sequence_range_length(ts->io, &ts->seq);
	break;
    }
	
    case DELETE_BASE: {
	int pos;

	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s delete_base position\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &pos);
	sequence_delete_base(ts->io, &ts->seq, pos, 1);
	sequence_range_length(ts->io, &ts->seq);
	break;
    }
	
    case REPLACE_BASE: {
	char base;
	int pos, conf;

	if (objc != 5) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s replace_base position call confidence\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &pos);
	base = *Tcl_GetString(objv[3]);
	Tcl_GetIntFromObj(interp, objv[4], &conf);
	cache_decr(ts->io, ts->seq);
	sequence_replace_base(ts->io, &ts->seq, pos, base, conf, 1);
	cache_incr(ts->io, ts->seq);
	break;
    }

    case GET_CLIPS:
	vTcl_SetResult(interp, "%d %d",
		       sequence_get_left(&ts->seq),
		       sequence_get_right(&ts->seq));
	break;

    case SET_CLIPS: {
	int left, right;

	if (objc != 4) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_clips left right\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &left);
	Tcl_GetIntFromObj(interp, objv[3], &right);
	cache_decr(ts->io, ts->seq);
	sequence_set_left (ts->io, &ts->seq, left);
	sequence_set_right(ts->io, &ts->seq, right);
	cache_incr(ts->io, ts->seq);

	break;
    }

    case SET_CLIPS_NO_INVALIDATE: {
	int left, right;

	if (objc != 4) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_clips left right\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &left);
	Tcl_GetIntFromObj(interp, objv[3], &right);
	cache_decr(ts->io, ts->seq);
	sequence_set_left_no_invalidate (ts->io, &ts->seq, left);
	sequence_set_right_no_invalidate(ts->io, &ts->seq, right);
	cache_incr(ts->io, ts->seq);

	break;
    }

    case MOVE_ANNOS: {
	int dist;

	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s move_annos dist\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &dist);
	cache_decr(ts->io, ts->seq);
	sequence_move_annos(ts->io, &ts->seq, dist);
	cache_incr(ts->io, ts->seq);

	break;
    }
    }

    return TCL_OK;
}

static int tcl_sequence_read(GapIO *io, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]) {
    seq_t *s;
    Tcl_WideInt snum;
    tcl_sequence *ts;
    Tcl_Obj *res;

    if (objc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s sequence_id\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetWideIntFromObj(interp, objv[1], &snum);
    s = (seq_t *)cache_search(io, GT_Seq, snum);

    if (NULL == (ts = (tcl_sequence *)ckalloc(sizeof(*ts))))
	return TCL_ERROR;
    ts->io = io;
    ts->seq = s;
    ts->iter = NULL;

    if (NULL == (res = Tcl_NewObj()))
	return TCL_ERROR;

    res->internalRep.otherValuePtr = ts;
    res->typePtr = &sequence_obj_type;
    sequence_update_string(res);

    /* Register the string form as a new command */
    if (s) cache_incr(io, s);
    if (!s ||
	NULL == Tcl_CreateObjCommand(interp, res->bytes, sequence_cmd,
				     (ClientData)ts,
				     (Tcl_CmdDeleteProc *)_cmd_delete))
	return TCL_ERROR;

    Tcl_SetObjResult(interp, res);

    return TCL_OK;
}

/* ------------------------------------------------------------------------ */
/* Tcl_Obj "anno_ele" type implementation */

static void anno_ele_update_string(Tcl_Obj *obj);
static int anno_ele_from_any(Tcl_Interp *interp, Tcl_Obj *obj);

typedef struct {
    GapIO *io;
    anno_ele_t *anno;
    btree_iter_t *iter;
} tcl_anno_ele;

static Tcl_ObjType anno_ele_obj_type = {
    "anno_ele",
    (Tcl_FreeInternalRepProc*)NULL,
    (Tcl_DupInternalRepProc*)NULL,
    anno_ele_update_string,
    anno_ele_from_any
};

static void anno_ele_update_string(Tcl_Obj *obj) {
    tcl_anno_ele *te = obj->internalRep.otherValuePtr;
    obj->bytes = ckalloc(30);
    obj->length = sprintf(obj->bytes, "anno_ele=%p", te);
}

static int anno_ele_from_any(Tcl_Interp *interp, Tcl_Obj *obj) {
    char *bytes;
    int length;
    tcl_anno_ele *te;

    if (NULL == (bytes = Tcl_GetStringFromObj(obj, &length)))
	return TCL_ERROR;

    if (0 != strncmp(bytes, "anno_ele=", 3))
	return TCL_ERROR;

    /* Free the old internalRep before setting the new one. */
    if (obj->typePtr && obj->typePtr->freeIntRepProc)
	(*obj->typePtr->freeIntRepProc)(obj);

    /* Convert the hex value to a pointer once more */
    if (1 != sscanf(bytes+3, "%p", &te))
	return TCL_ERROR;

    obj->internalRep.otherValuePtr = te;
    obj->typePtr = &anno_ele_obj_type;
    return TCL_OK;
}

static int anno_ele_cmd(ClientData clientData, Tcl_Interp *interp,
			int objc, Tcl_Obj *CONST objv[]) {
    int index;
    tcl_anno_ele *te = (tcl_anno_ele *)clientData;

    static char *options[] = {
	"delete",       "io",           "get_rec",
	"get_contig",   "get_position", "get_comment",
	"get_obj_type", "get_obj_rec",  "get_anno_rec",
	"get_type",     "get_direction","move",
	"set_contig",   "set_position", "set_comment",
	"set_obj_type", "set_obj_rec",  "set_anno_rec",
	"set_type",     "set_direction","remove",
	(char *)NULL,
    };

    enum options {
	DELETE,         IO,             GET_REC,
	GET_CONTIG,     GET_POSITION,   GET_COMMENT,
	GET_OBJ_TYPE,   GET_OBJ_REC,    GET_ANNO_REC,
	GET_TYPE,       GET_DIRECTION,  MOVE,
	SET_CONTIG,     SET_POSITION,   SET_COMMENT,
	SET_OBJ_TYPE,   SET_OBJ_REC,    SET_ANNO_REC,
	SET_TYPE,	SET_DIRECTION,  REMOVE,
    };

    if (objc < 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "option arg ?arg ...?");
	return TCL_ERROR;
    }

    if (Tcl_GetIndexFromObj(interp, objv[1], options, "option", 0,
            &index) != TCL_OK) {
        return TCL_ERROR;
    }

    /* Get read/write if appropriate */
    switch ((enum options)index) {
	anno_ele_t *t;
    case SET_CONTIG:
    case SET_POSITION:
    case SET_OBJ_TYPE:
    case SET_OBJ_REC:
    case SET_COMMENT:
    case SET_ANNO_REC:
    case SET_TYPE:
    case MOVE:
	if (NULL == (t = cache_rw(te->io, te->anno)))
	    return TCL_ERROR;
	te->anno = t;
	break;
    default:
	;
    }

    /* Perform the command proper */
    switch ((enum options)index) {
    case REMOVE: {
	tg_rec contig;

	if (anno_get_range(te->io, te->anno->rec, &contig, 0)) {
	    contig_t *c = cache_search(te->io, GT_Contig, contig);
	    bin_remove_item(te->io, &c, GT_AnnoEle, te->anno->rec);
	}

	anno_ele_destroy(te->io, te->anno);
	/* Deliberate flow through to DELETE */
    }

    case DELETE:
	Tcl_DeleteCommandFromToken(interp,
				   Tcl_GetCommandFromObj(interp, objv[0]));
	break;

    case IO:
	Tcl_SetResult(interp, io_obj_as_string(te->io) , TCL_VOLATILE);
	break;

    case GET_REC:
	Tcl_SetWideIntObj(Tcl_GetObjResult(interp), te->anno->rec);
	break;

    case GET_CONTIG: {
	tg_rec contig;
	if (NULL == anno_get_range(te->io, te->anno->rec, &contig, 0))
	    return TCL_ERROR;

	Tcl_SetWideIntObj(Tcl_GetObjResult(interp), contig);
	
	break;
    }

    case GET_POSITION: {
	tg_rec contig;
	range_t *r = anno_get_range(te->io, te->anno->rec, &contig, 1);
	if (NULL == r)
	    return TCL_ERROR;

	vTcl_SetResult(interp, "%d %d %"PRIrec, 
		       r->start, r->end, contig);
	break;
    }

    case SET_CONTIG:
    case SET_POSITION:
	puts("Unimplemented\n");
	break;

    case MOVE: {
	int otype, start, end;
	tg_rec orec, contig;
	range_t r, *r_out;
	tg_rec seq_bin = 0, crec;
	contig_t *c;
	bin_index_t *bin;

	if (objc != 6) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_position obj_type obj_rec start end\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &otype);
	Tcl_GetWideIntFromObj(interp, objv[3], &orec);
	Tcl_GetIntFromObj(interp, objv[4], &start);
	Tcl_GetIntFromObj(interp, objv[5], &end);
	
	if (otype != GT_Seq && otype != GT_Contig) {
	    vTcl_SetResult(interp, "obj_type should be %d or %d\n",
			   GT_Seq, GT_Contig);
	    return TCL_ERROR;
	}

	if (anno_get_range(te->io, te->anno->rec, &contig, 0)) {
	    contig_t *c = cache_search(te->io, GT_Contig, contig);
	    bin_remove_item(te->io, &c, GT_AnnoEle, te->anno->rec);
	}

	memset(&r, 0, sizeof(r));
	r.flags    = GRANGE_FLAG_ISANNO;
	r.rec      = te->anno->rec;
	r.mqual    = te->anno->tag_type;
	r.pair_rec = orec;

	if (otype == GT_Seq) {
	    int st, en;
	    sequence_get_position2(te->io, orec, &crec, &st, &en, NULL,
				   &seq_bin, NULL, NULL);
	    start += st;
	    end   += st;
	    
	    r.flags |= GRANGE_FLAG_TAG_SEQ;
	} else {
	    crec = orec;
	}

	r.start    = start;
	r.end      = end;

	c = cache_search(te->io, GT_Contig, crec);
	if (!c)
	    return TCL_ERROR;
	cache_incr(te->io, c);

	/* Place tag in new location */
	if (seq_bin)
	    bin = bin_add_to_range(te->io, &c, seq_bin, &r, &r_out, NULL, 0);
	else
	    bin = bin_add_range(te->io, &c, &r, &r_out, NULL, 0);

	te->anno = cache_rw(te->io, te->anno);
	te->anno->obj_type = otype;
	te->anno->obj_rec  = orec;

	/* Update bin backref */	
	if (te->anno->bin != bin->rec) {
	    te->anno->bin = bin->rec;
	}

	cache_decr(te->io, c);

	break;
    }

    case GET_COMMENT:
	Tcl_SetStringObj(Tcl_GetObjResult(interp),
			 te->anno->comment ? te->anno->comment : "",
			 -1);
	break;

    case SET_COMMENT: {
	char *str;

	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_comment string\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	str = Tcl_GetStringFromObj(objv[2], NULL);
	anno_ele_set_comment(te->io, &te->anno, str);
	break;
    }

    case GET_DIRECTION:
	Tcl_SetStringObj(Tcl_GetObjResult(interp), &te->anno->direction ,1);
	break;

    case SET_DIRECTION: {
	char *dir;

	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_direction +|-|.|?\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	dir = Tcl_GetStringFromObj(objv[2], NULL);
	if (-1 == anno_ele_set_direction(te->io, &te->anno, *dir))
	    return TCL_ERROR;

	break;
    }

    case GET_TYPE: {
	char type[5];
	(void)type2str(te->anno->tag_type, type);
	Tcl_SetStringObj(Tcl_GetObjResult(interp), type ,4);
	break;
    }

    case SET_TYPE:
	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_type type_str\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}
	
	if (-1 == anno_ele_set_type(te->io, &te->anno,
				    Tcl_GetStringFromObj(objv[2], NULL)))
	    return TCL_ERROR;

	break;

    case GET_OBJ_TYPE:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), te->anno->obj_type);
	break;

    case SET_OBJ_TYPE: {
	int otype;

	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_obj_type integer_type_code\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetIntFromObj(interp, objv[2], &otype);
	te->anno->obj_type = otype;
	break;
    }

    case GET_OBJ_REC:
	if (te->anno->obj_type == GT_Contig) {
	    tg_rec contig;
	    if (NULL == (anno_get_range(te->io, te->anno->rec, &contig, 0)))
		return TCL_ERROR;
	    
	    Tcl_SetWideIntObj(Tcl_GetObjResult(interp), contig);
	} else {
	    Tcl_SetWideIntObj(Tcl_GetObjResult(interp), te->anno->obj_rec);
	}
	break;

    case SET_OBJ_REC: {
	Tcl_WideInt orec;

	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_obj_rec record_no\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetWideIntFromObj(interp, objv[2], &orec);
	te->anno->obj_rec = (tg_rec)orec;
	break;
    }

    case GET_ANNO_REC:
	Tcl_SetWideIntObj(Tcl_GetObjResult(interp), te->anno->anno_rec);
	break;

    case SET_ANNO_REC: {
	Tcl_WideInt rec;

	if (objc != 3) {
	    vTcl_SetResult(interp, "wrong # args: should be "
			   "\"%s set_anno_rec record_no\"\n",
			   Tcl_GetStringFromObj(objv[0], NULL));
	    return TCL_ERROR;
	}

	Tcl_GetWideIntFromObj(interp, objv[2], &rec);
	te->anno->anno_rec = (tg_rec)rec;
	break;
    }
    }

    return TCL_OK;
}

static int tcl_anno_ele_read(GapIO *io, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]) {
    anno_ele_t *e;
    Tcl_WideInt elenum;
    tcl_anno_ele *te;
    Tcl_Obj *res;

    if (objc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s anno_ele_id\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetWideIntFromObj(interp, objv[1], &elenum);
    e = (anno_ele_t *)cache_search(io, GT_AnnoEle, elenum);

    if (NULL == (te = (tcl_anno_ele *)ckalloc(sizeof(*te))))
	return TCL_ERROR;
    te->io = io;
    te->anno = e;
    te->iter = NULL;

    if (NULL == (res = Tcl_NewObj()))
	return TCL_ERROR;

    res->internalRep.otherValuePtr = te;
    res->typePtr = &anno_ele_obj_type;
    anno_ele_update_string(res);

    /* Register the string form as a new command */
    cache_incr(io, e);
    if (NULL == Tcl_CreateObjCommand(interp, res->bytes, anno_ele_cmd,
				     (ClientData)te,
				     (Tcl_CmdDeleteProc *)_cmd_delete))
	return TCL_ERROR;

    Tcl_SetObjResult(interp, res);

    return TCL_OK;
}

/* ------------------------------------------------------------------------ */
/* Tcl_Obj "library" type implementation */

static void library_update_string(Tcl_Obj *obj);
static int library_from_any(Tcl_Interp *interp, Tcl_Obj *obj);

typedef struct {
    GapIO *io;
    library_t *library;
    btree_iter_t *iter;
} tcl_library;

static Tcl_ObjType library_obj_type = {
    "library",
    (Tcl_FreeInternalRepProc*)NULL,
    (Tcl_DupInternalRepProc*)NULL,
    library_update_string,
    library_from_any
};

static void library_update_string(Tcl_Obj *obj) {
    tcl_library *l = obj->internalRep.otherValuePtr;
    obj->bytes = ckalloc(30);
    obj->length = sprintf(obj->bytes, "library=%p", l);
}

static int library_cmd(ClientData clientData, Tcl_Interp *interp,
		       int objc, Tcl_Obj *CONST objv[]) {
    int index;
    tcl_library *tl = (tcl_library *)clientData;

    static char *options[] = {
	"delete",         "io",           "get_rec",
	"get_orient",     "get_machine",  "get_dist",
	"get_insert_size","get_insert_sd","get_count",
	"get_name",	  "update_stats",
	(char *)NULL,
    };

    enum options {
	DELETE,          IO,             GET_REC,
	GET_ORIENT,      GET_MACHINE,    GET_DIST,
	GET_INSERT_SIZE, GET_INSERT_SD,  GET_COUNT,
	GET_NAME,	 UPDATE_STATS
    };

    if (objc < 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "option arg ?arg ...?");
	return TCL_ERROR;
    }

    if (Tcl_GetIndexFromObj(interp, objv[1], options, "option", 0,
            &index) != TCL_OK) {
        return TCL_ERROR;
    }

    switch ((enum options)index) {
    case DELETE:
	Tcl_DeleteCommandFromToken(interp,
				   Tcl_GetCommandFromObj(interp, objv[0]));
	break;

    case IO:
	Tcl_SetResult(interp, io_obj_as_string(tl->io) , TCL_VOLATILE);
	break;

    case GET_REC:
	Tcl_SetWideIntObj(Tcl_GetObjResult(interp), tl->library->rec);
	break;

    case GET_ORIENT:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), tl->library->lib_type);
	break;

    case GET_MACHINE:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), tl->library->machine);
	break;

    case GET_INSERT_SIZE: {
	Tcl_Obj *obj[3];
	obj[0] = Tcl_NewIntObj(tl->library->insert_size[0]);
	obj[1] = Tcl_NewIntObj(tl->library->insert_size[1]);
	obj[2] = Tcl_NewIntObj(tl->library->insert_size[2]);
	Tcl_SetObjResult(interp, Tcl_NewListObj(3, obj));
	break;
    }

    case GET_INSERT_SD: {
	Tcl_Obj *obj[3];
	obj[0] = Tcl_NewDoubleObj(tl->library->sd[0]);
	obj[1] = Tcl_NewDoubleObj(tl->library->sd[1]);
	obj[2] = Tcl_NewDoubleObj(tl->library->sd[2]);
	Tcl_SetObjResult(interp, Tcl_NewListObj(3, obj));
	break;
    }

    case GET_DIST: {
	Tcl_Obj *lo[3];
	int i, j;
	for (j = 0; j < 3; j++) {
	    lo[j] = Tcl_NewListObj(0, NULL);
	    for (i = 0; i < LIB_BINS; i++) {
		int b;
		double d;
		if (!tl->library->size_hist[j][i]) continue;

		b = ibin2isize(i);
		d = (double)tl->library->size_hist[j][i] / ibin_width(i);
		Tcl_ListObjAppendElement(interp, lo[j], Tcl_NewIntObj(b));
		Tcl_ListObjAppendElement(interp, lo[j], Tcl_NewDoubleObj(d));
	    }
	}
	Tcl_SetObjResult(interp, Tcl_NewListObj(3, lo));
	break;
    }

    case GET_COUNT: {
	Tcl_Obj *obj[3];
	int i, j;

	for (j = 0; j < 3; j++) {
	    int c = 0;
	    for (i = 0; i < LIB_BINS; i++) 
		c += tl->library->size_hist[j][i];
	    obj[j] = Tcl_NewIntObj(c);
	}
        Tcl_SetObjResult(interp, Tcl_NewListObj(3, obj));
	break;
    }

    case GET_NAME:
	if (tl->library->name) {
	    Tcl_SetStringObj(Tcl_GetObjResult(interp), tl->library->name, -1);
	} else {
	    char buf[100];
	    sprintf(buf, "rec#%"PRIrec, tl->library->rec);
	    Tcl_SetStringObj(Tcl_GetObjResult(interp), buf, -1);
	}
	break;

    case UPDATE_STATS:
	update_library_stats(tl->io, tl->library->rec, 100, NULL, NULL, NULL);
	break;
    }

    return TCL_OK;
}

static int library_from_any(Tcl_Interp *interp, Tcl_Obj *obj) {
    char *bytes;
    int length;
    tcl_library *l;

    if (NULL == (bytes = Tcl_GetStringFromObj(obj, &length)))
	return TCL_ERROR;

    if (0 != strncmp(bytes, "library=", 3))
	return TCL_ERROR;

    /* Free the old internalRep before setting the new one. */
    if (obj->typePtr && obj->typePtr->freeIntRepProc)
	(*obj->typePtr->freeIntRepProc)(obj);

    /* Convert the hex value to a pointer once more */
    if (1 != sscanf(bytes+3, "%p", &l))
	return TCL_ERROR;

    obj->internalRep.otherValuePtr = l;
    obj->typePtr = &library_obj_type;
    return TCL_OK;
}

static int tcl_library_read(GapIO *io, Tcl_Interp *interp,
			   int objc, Tcl_Obj *CONST objv[]) {
    Tcl_WideInt rec;
    library_t *l;
    tcl_library *tl;
    Tcl_Obj *res;

    if (objc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s contig_id\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetWideIntFromObj(interp, objv[1], &rec);
    l = (library_t *)cache_search(io, GT_Library, rec);

    if (NULL == (tl = (tcl_library *)ckalloc(sizeof(*tl))))
	return TCL_ERROR;
    tl->io = io;
    tl->library = l;
    tl->iter = NULL;

    if (NULL == (res = Tcl_NewObj()))
	return TCL_ERROR;

    res->internalRep.otherValuePtr = tl;
    res->typePtr = &library_obj_type;
    library_update_string(res);

    /* Register the string form as a new command */
    cache_incr(io, l);
    if (NULL == Tcl_CreateObjCommand(interp, res->bytes, library_cmd,
				     (ClientData)tl,
				     (Tcl_CmdDeleteProc *)_cmd_delete))
	return TCL_ERROR;

    Tcl_SetObjResult(interp, res);

    return TCL_OK;
}

/* ------------------------------------------------------------------------ */
/* Tcl_Obj "database" type implementation */

static void database_update_string(Tcl_Obj *obj);
static int database_from_any(Tcl_Interp *interp, Tcl_Obj *obj);

static Tcl_ObjType database_obj_type = {
    "database",
    (Tcl_FreeInternalRepProc*)NULL,
    (Tcl_DupInternalRepProc*)NULL,
    database_update_string,
    database_from_any
};

static void database_update_string(Tcl_Obj *obj) {
    GapIO *io = obj->internalRep.otherValuePtr;
    obj->bytes = ckalloc(30);
    obj->length = sprintf(obj->bytes, "database=%p", io);
}

static int database_from_any(Tcl_Interp *interp, Tcl_Obj *obj) {
    char *bytes;
    int length;
    GapIO *io;

    if (NULL == (bytes = Tcl_GetStringFromObj(obj, &length)))
	return TCL_ERROR;

    if (0 != strncmp(bytes, "database=", 3))
	return TCL_ERROR;

    /* Free the old internalRep before setting the new one. */
    if (obj->typePtr && obj->typePtr->freeIntRepProc)
	(*obj->typePtr->freeIntRepProc)(obj);

    /* Convert the hex value to a pointer once more */
    if (1 != sscanf(bytes+3, "%p", &io))
	return TCL_ERROR;

    obj->internalRep.otherValuePtr = io;
    obj->typePtr = &database_obj_type;
    return TCL_OK;
}

static int database_cmd(ClientData clientData, Tcl_Interp *interp,
			int objc, Tcl_Obj *CONST objv[]) {
    int index;
    GapIO *io = (GapIO *)clientData;

    static char *options[] = {
	"get_num_contigs",  "flush", "get_num_libraries",
	"get_library_rec",
	(char *)NULL,
    };

    enum options {
	GET_NUM_CONTIGS,     FLUSH,   GET_NUM_LIBRARIES,
	GET_LIBRARY_REC
    };

    if (objc < 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "option arg ?arg ...?");
	return TCL_ERROR;
    }

    if (Tcl_GetIndexFromObj(interp, objv[1], options, "option", 0,
            &index) != TCL_OK) {
        return TCL_ERROR;
    }

    switch ((enum options)index) {
    case GET_NUM_CONTIGS:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), io->db->Ncontigs);
	break;

    case GET_NUM_LIBRARIES:
	Tcl_SetIntObj(Tcl_GetObjResult(interp), io->db->Nlibraries);
	break;

    case GET_LIBRARY_REC: {
	int idx;
	if (objc != 3) {
	    Tcl_WrongNumArgs(interp, 1, objv, "get_library_rec library_index");
	    return TCL_ERROR;
	}
	
	Tcl_GetIntFromObj(interp, objv[2], &idx);
	Tcl_SetWideIntObj(Tcl_GetObjResult(interp),
			  arr(tg_rec, io->library, idx));
	break;
    }

    case FLUSH:
	cache_flush(io);
	break;
    }

    return TCL_OK;
}

static int tcl_database_read(GapIO *io, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]) {
    Tcl_Obj *res;

    if (objc != 1) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    if (NULL == (res = Tcl_NewObj()))
	return TCL_ERROR;

    res->internalRep.otherValuePtr = io;
    res->typePtr = &database_obj_type;
    database_update_string(res);

    /* Register the string form as a new command */
    if (NULL == Tcl_CreateObjCommand(interp, res->bytes, database_cmd,
				     (ClientData)io,
				     (Tcl_CmdDeleteProc *)NULL))
	return TCL_ERROR;

    Tcl_SetObjResult(interp, res);
    return TCL_OK;
}


/* ------------------------------------------------------------------------ */
/*
 * Our low-level Tcl interface functions to the GapIO type.
 */

typedef struct {
    char *db_name;
    char *access_str;
    int create;
} open_db_arg;
static int tcl_open_database(ClientData clientData, Tcl_Interp *interp,
			     int objc, Tcl_Obj *CONST objv[]) {
    GapIO *io;
    Tcl_Obj *iobj;
    int ro = 1;
    

    open_db_arg args;
    cli_args a[] = {
	{"-name",    ARG_STR, 1, NULL, offsetof(open_db_arg, db_name)},
	{"-create",  ARG_INT, 1, "0",  offsetof(open_db_arg, create)},
	{"-access",  ARG_STR, 1, "r",  offsetof(open_db_arg, access_str)},
	{NULL,	     0,	      0, NULL, 0}
    };

    vfuncheader("open database");

    if (-1 == gap_parse_obj_args(a, &args, objc, objv))
	return TCL_ERROR;

    if (strcmp(args.access_str, "rw") == 0) {
	ro = 0;
    }

    if (NULL == (io = gio_open(args.db_name, ro, args.create)))
	return TCL_ERROR;

    if (NULL == (iobj = Tcl_NewObj()))
	return TCL_ERROR;
    
    iobj->internalRep.otherValuePtr = (VOID *)io;
    iobj->typePtr = &io_obj_type;
    io_update_string(iobj);

    /* Register the string form as a new command */
    if (NULL == Tcl_CreateObjCommand(interp, iobj->bytes, io_cmd,
				     (ClientData)io,
				     (Tcl_CmdDeleteProc *)NULL))
	return TCL_ERROR;

    Tcl_SetObjResult(interp, iobj);
    
    return TCL_OK;
}


/* ------------------------------------------------------------------------ */
/*
 * The registration function called by tcl when the library is dynamically
 * loaded.
 */
int G5_Init(Tcl_Interp *interp) {
    Tcl_RegisterObjType(&io_obj_type);
    Tcl_RegisterObjType(&database_obj_type);
    Tcl_RegisterObjType(&contig_obj_type);
    Tcl_RegisterObjType(&sequence_obj_type);
    Tcl_RegisterObjType(&anno_ele_obj_type);
    Tcl_RegisterObjType(&library_obj_type);

    if (NULL == Tcl_CreateObjCommand(interp, "g5::open_database",
				     tcl_open_database,
				     (ClientData)NULL,
				     (Tcl_CmdDeleteProc *)NULL))
	return TCL_ERROR;

    return TCL_OK;
}

int G5_SafeInit(Tcl_Interp *interp) {
    return G5_Init(interp);
}

int G5_Unload(Tcl_Interp *interp, int flags) {
    Tcl_SetResult(interp, "Pkg_Unload() function not implemented",
		  TCL_STATIC);
    return TCL_ERROR;
}

int G5_SafeUnload(Tcl_Interp *interp, int flags) {
    Tcl_SetResult(interp, "Pkg_SafeUnload() function not implemented",
		  TCL_STATIC);
    return TCL_ERROR;
}


/*

Notes on the object system used here.

We initially create an io object via:

    set io [g5::open_database -name foo]

We then run methods on io to obtain new objects:

    set seq [$io get_sequence 10]
    set cnt [$io get_contig 3]

Finally we query sequence and contig objects to obtain data about them:

    puts "sequence length = [$seq get_length]"
    set s_ids [$cnt seqs_in_range 1 1000]

Note that in all these cases the tcl object exists as both a function
name at global scope and a scalar variable containing the name of that
function. Eg:

proc foo {} {
    set seq [$io get_sequence 10]
    puts seq=$seq
}
foo

Running this may print up seq=0x276363, which is the contents of $seq
and the name of the function.

After foo exits seq is no longer in scope and is freed, but that's
just a variable holding the name of the global function, which itself
where the real data is held in memory. If this was interactive we
could even type "seq=0x276363 get_length" now.

To remove the memory leak and function pollution caused by this we
have a remove method. Ie the above function should be:

proc foo {} {
    set seq [$io get_sequence 10]
    puts seq=$seq
    $seq delete
}
foo

(At the time of writing this the delete method exists for sequences and
contigs only.)



Suggested change
----------------

A procedural method instead of object would work. Eg

set seq [g5::io::get_sequence $io 10]
puts "Length = [g5::sequence::get_length $seq]

This is how perl works, but the addition of knowing which package a scalar is
in (via the "bless" command) and the -> operator means that
"$seq->get_length" is translated to "sequence::get_length($seq)" automatically.

Is it possible to implement this in Tcl?



*/
