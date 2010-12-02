/*
 * Interfaces to the raw data structures, plus a few utility routines for
 * manipulating these.
 *
 * The general approach is to use TclX keyed lists for translation from
 * C structures to lists and back again.
 */

#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <tclXkeylist.h>

#include "IO.h"
#include "cli_arg.h"
#include "gap_globals.h"
#include "gap-error.h"
#include "tagUtils.h"
#include "tcl_utils.h"
#include "clones.h"
#include "gap-tcl.h"

static int auto_flush = 0;

/*
 *-----------------------------------------------------------------------------
 * Converting C structures to Keyed Lists
 *-----------------------------------------------------------------------------
 */
#define setl(s,x) \
do { \
    Tcl_Obj *o = Tcl_NewIntObj((s)->x); \
    TclX_KeyedListSet(interp, list, w(#x), o); \
} while (0)


#define setl_tag_type(s,x) \
do { \
    char buf[5]; \
    Tcl_Obj *o = Tcl_NewStringObj(type2str((s)->x, buf), -1); \
    TclX_KeyedListSet(interp, list, w(#x), o); \
} while (0)

Tcl_Obj *GAnnotations_klist(Tcl_Interp *interp, GapIO *io, GAnnotations *a) {
    Tcl_Obj *list = TclX_NewKeyedListObj();

    setl_tag_type(a, type);
    setl(a, position); 
    setl(a, length); 
    setl(a, strand); 
    setl(a, annotation); 
    setl(a, next);

    return list;
}

Tcl_Obj *GNotes_klist(Tcl_Interp *interp, GapIO *io, GNotes *a) {
    Tcl_Obj *list = TclX_NewKeyedListObj();

    setl_tag_type(a, type);
    setl(a, ctime); 
    setl(a, mtime); 
    setl(a, annotation); 
    setl(a, next); 
    setl(a, prev); 
    setl(a, prev_type); 

    return list;
}

Tcl_Obj *GVectors_klist(Tcl_Interp *interp, GapIO *io, GVectors *v) {
    Tcl_Obj *list = TclX_NewKeyedListObj();

    setl(v, name);             /* vector name */
    setl(v, level);            /* 1=clone, 2=subclone, etc */

    return list;
}

Tcl_Obj *GReadings_klist(Tcl_Interp *interp, GapIO *io, GReadings *r) {
    Tcl_Obj *list = TclX_NewKeyedListObj();

    setl(r, name);
    setl(r, trace_name);
    setl(r, trace_type);
    setl(r, left);             /* left neighbour */
    setl(r, right);            /* right neighbour */
    setl(r, position);         /* position in contig */
    setl(r, length);           /* total length of reading */
    setl(r, sense);            /* 0 = original, 1 = reverse */
    setl(r, sequence);
    setl(r, confidence);
    setl(r, orig_positions);
    setl(r, chemistry);        /* bit 0 = primer/terminator, bits 1-4 = type */
    setl(r, annotations);      /* start of annotation list */
    setl(r, sequence_length);  /* clipped length */
    setl(r, start);            /* last base of left cutoff */
    setl(r, end);              /* first base of right cutoff */
    setl(r, template);         /* aka subclone */
    setl(r, strand);           /* 0 = forward, 1 = reverse */
    setl(r, primer);           /* 0 = unknown, 1 = forw, 2 = rev, 3 = custom */
    setl(r, notes);	       /* start of notes list */

    return list;
}

Tcl_Obj *GContigs_klist(Tcl_Interp *interp, GapIO *io, GContigs *c) {
    Tcl_Obj *list = TclX_NewKeyedListObj();

    setl(c, left); 
    setl(c, right); 
    setl(c, length);
    setl(c, annotations);      /* start of annotation list */
    setl(c, notes);	       /* start of notes list */

    return list;
}

Tcl_Obj *GDatabase_klist(Tcl_Interp *interp, GapIO *io, GDatabase *d) {
    Tcl_Obj *list = TclX_NewKeyedListObj();

    setl(d, version);
    setl(d, maximum_db_size);  /* MAXDB */
    setl(d, actual_db_size);   /* */
    setl(d, max_gel_len);      /* 4096 */
    setl(d, data_class);       /* 0 = DNA, 1 = protein */

    setl(d, num_contigs);      /* number of contigs used */
    setl(d, num_readings);     /* number of readings used */

    /* Bitmaps */
    setl(d, Nfreerecs);        /* number of bits */
    setl(d, freerecs);

    /* Arrays */
    setl(d, Ncontigs);         /* elements in array */
    setl(d, contigs);          /* records that are GT_Contigs */

    setl(d, Nreadings);        /* elements in array */
    setl(d, readings);         /* records that are GT_Readings */

    setl(d, Nannotations);     /* elements in array */
    setl(d, annotations);      /* records that are GT_Annotations */
    setl(d, free_annotations); /* head of list of free annotations */

    setl(d, Ntemplates);       /* elements in array */
    setl(d, templates);        /* records that are GT_Templates */

    setl(d, Nclones);          /* elements in array */
    setl(d, clones);           /* records that are GT_Templates */

    setl(d, Nvectors);         /* elements in array */
    setl(d, vectors);          /* records that are  GT_Vectors */

    setl(d, contig_order);     /* Array record number */

    setl(d, Nnotes);	       /* elements in array */
    setl(d, notes_a);          /* records that are GT_Annotations */
    setl(d, notes);            /* start of database notes list */
    setl(d, free_notes);       /* head of list of free notes */

    return list;
}


Tcl_Obj *GTemplates_klist(Tcl_Interp *interp, GapIO *io, GTemplates *t) {
    Tcl_Obj *list = TclX_NewKeyedListObj();

    setl(t, name);
    setl(t, strands);
    setl(t, vector);
    setl(t, clone);
    setl(t, insert_length_min);
    setl(t, insert_length_max);

    return list;
}

Tcl_Obj *GClones_klist(Tcl_Interp *interp, GapIO *io, GClones *c) {
    Tcl_Obj *list = TclX_NewKeyedListObj();

    setl(c, name);
    setl(c, vector);

    return list;
}

#undef setl
#undef setl_tag_type

/*
 *-----------------------------------------------------------------------------
 * Converting Keyed lists to C structures
 *-----------------------------------------------------------------------------
 */
#define getl(s, x) \
do { \
    if (TCL_OK == TclX_KeyedListGet(interp, list, w(#x), &ptr)) { \
        (void)Tcl_GetIntFromObj(interp, ptr, &(s)->x); \
    } \
} while (0);

#define getl_tag_type(s, x) \
do { \
    if (TCL_OK == TclX_KeyedListGet(interp, list, w(#x), &ptr)) { \
	(s)->x = str2type(Tcl_GetStringFromObj(ptr, NULL)); \
    } \
} while (0);

int klist_GAnnotations(Tcl_Interp *interp, GapIO *io, GAnnotations *a,
		       Tcl_Obj *list) {
    Tcl_Obj *ptr;

    getl_tag_type(a, type);
    getl(a, position); 
    getl(a, length); 
    getl(a, strand); 
    getl(a, annotation); 
    getl(a, next);

    return 0;
}

int klist_GNotes(Tcl_Interp *interp, GapIO *io, GNotes *a,
		 Tcl_Obj *list) {
    Tcl_Obj *ptr;

    getl_tag_type(a, type);
    getl(a, ctime); 
    getl(a, mtime); 
    getl(a, annotation); 
    getl(a, next);
    getl(a, prev);
    getl(a, prev_type); 

    return 0;
}

int klist_GVectors(Tcl_Interp *interp, GapIO *io, GVectors *v, Tcl_Obj *list) {
    Tcl_Obj *ptr;

    getl(v, name);             /* vector name */
    getl(v, level);            /* 1=clone, 2=subclone, etc */

    return 0;
}

int klist_GReadings(Tcl_Interp *interp, GapIO *io, GReadings *r,
		    Tcl_Obj *list) {
    Tcl_Obj *ptr;

    getl(r, name);
    getl(r, trace_name);
    getl(r, trace_type);
    getl(r, left);             /* left neighbour */
    getl(r, right);            /* right neighbour */
    getl(r, position);         /* position in contig */
    getl(r, length);           /* total length of reading */
    getl(r, sense);            /* 0 = original, 1 = reverse */
    getl(r, sequence);
    getl(r, confidence);
    getl(r, orig_positions);
    getl(r, chemistry);        /* bit 0 = primer/terminator, bits 1-4 = type */
    getl(r, annotations);      /* start of annotation list */
    getl(r, sequence_length);  /* clipped length */
    getl(r, start);            /* last base of left cutoff */
    getl(r, end);              /* first base of right cutoff */
    getl(r, template);         /* aka subclone */
    getl(r, strand);           /* 0 = forward, 1 = reverse */
    getl(r, primer);           /* 0 = unknown, 1 = forw, 2 = rev, 3 = custom */
    getl(r, notes);	       /* start of notes list */

    return 0;
}

int klist_GContigs(Tcl_Interp *interp, GapIO *io, GContigs *c, Tcl_Obj *list) {
    Tcl_Obj *ptr;

    getl(c, left); 
    getl(c, right); 
    getl(c, length);
    getl(c, annotations);      /* start of annotation list */
    getl(c, notes);	       /* start of notes list */

    return 0;
}

int klist_GDatabase(Tcl_Interp *interp, GapIO *io, GDatabase *d,
		    Tcl_Obj *list) {
    Tcl_Obj *ptr;

    getl(d, version);
    getl(d, maximum_db_size);  /* MAXDB */
    getl(d, actual_db_size);   /* */
    getl(d, max_gel_len);      /* 4096 */
    getl(d, data_class);       /* 0 = DNA, 1 = protein */

    getl(d, num_contigs);      /* number of contigs used */
    getl(d, num_readings);     /* number of readings used */

    /* Bitmaps */
    getl(d, Nfreerecs);        /* number of bits */
    getl(d, freerecs);

    /* Arrays */
    getl(d, Ncontigs);         /* elements in array */
    getl(d, contigs);          /* records that are GT_Contigs */

    getl(d, Nreadings);        /* elements in array */
    getl(d, readings);         /* records that are GT_Readings */

    getl(d, Nannotations);     /* elements in array */
    getl(d, annotations);      /* records that are GT_Annotations */
    getl(d, free_annotations); /* head of list of free annotations */

    getl(d, Ntemplates);       /* elements in array */
    getl(d, templates);        /* records that are GT_Templates */

    getl(d, Nclones);          /* elements in array */
    getl(d, clones);           /* records that are GT_Templates */

    getl(d, Nvectors);         /* elements in array */
    getl(d, vectors);          /* records that are  GT_Vectors */

    getl(d, contig_order);     /* Array record number */

    getl(d, Nnotes);	       /* elements in array */
    getl(d, notes_a);          /* records that are GT_Annotations */
    getl(d, notes);            /* start of database notes list */
    getl(d, free_notes);       /* head of list of free notes */

    return 0;
}


int klist_GTemplates(Tcl_Interp *interp, GapIO *io, GTemplates *t,
		     Tcl_Obj *list) {
    Tcl_Obj *ptr;

    getl(t, name);
    getl(t, strands);
    getl(t, vector);
    getl(t, clone);
    getl(t, insert_length_min);
    getl(t, insert_length_max);

    return 0;
}

int klist_GClones(Tcl_Interp *interp, GapIO *io, GClones *c, Tcl_Obj *list) {
    Tcl_Obj *ptr;

    getl(c, name);
    getl(c, vector);

    return 0;
}

#undef getl
#undef getl_tag_type


/*
 *-----------------------------------------------------------------------------
 * Tcl interfaces to the above commands.
 *-----------------------------------------------------------------------------
 */
int tcl_read_annotation(ClientData clientData, Tcl_Interp *interp,
			int argc, char **argv) {
    GAnnotations a;
    GapIO *io;
    int handle;
    int err;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    err = GT_Read(io, arr(GCardinal, io->annotations, atoi(argv[2])-1),
		  &a, sizeof(a), GT_Annotations);

    if (err) {
	Tcl_ResetResult(interp);
    } else {
	Tcl_SetObjResult(interp, GAnnotations_klist(interp, io, &a));
    }

    return TCL_OK;
}

int tcl_read_note(ClientData clientData, Tcl_Interp *interp,
		  int argc, char **argv) {
    GNotes n;
    GapIO *io;
    int handle;
    int err;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    err = GT_Read(io, arr(GCardinal, io->notes, atoi(argv[2])-1),
		  &n, sizeof(n), GT_Notes);

    if (err) {
	Tcl_ResetResult(interp);
    } else {
	Tcl_SetObjResult(interp, GNotes_klist(interp, io, &n));
    }

    return TCL_OK;
}

int tcl_read_vector(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {
    GVectors v;
    GapIO *io;
    int handle;
    int err;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    err = GT_Read(io, arr(GCardinal, io->vectors, atoi(argv[2])-1),
		  &v, sizeof(v), GT_Vectors);

    if (err) {
	Tcl_ResetResult(interp);
    } else {
	Tcl_SetObjResult(interp, GVectors_klist(interp, io, &v));
    }

    return TCL_OK;
}

int tcl_read_reading(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv) {
    GReadings r;
    GapIO *io;
    int handle;
    int err;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    err = gel_read(io, atoi(argv[2]), r);

    if (err) {
	Tcl_ResetResult(interp);
    } else {
	Tcl_SetObjResult(interp, GReadings_klist(interp, io, &r));
    }

    return TCL_OK;
}

int tcl_read_contig(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {
    GContigs c;
    GapIO *io;
    int handle;
    int err;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    err = GT_Read(io, arr(GCardinal, io->contigs, atoi(argv[2])-1),
		  &c, sizeof(c), GT_Contigs);

    if (err) {
	Tcl_ResetResult(interp);
    } else {
	Tcl_SetObjResult(interp, GContigs_klist(interp, io, &c));
    }

    return TCL_OK;
}

int tcl_read_database(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    GapIO *io;
    int handle;

    if (argc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    Tcl_SetObjResult(interp, GDatabase_klist(interp, io, &io->db));

    return TCL_OK;
}

int tcl_read_template(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    GTemplates t;
    GapIO *io;
    int handle;
    int err;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    err = GT_Read(io, arr(GCardinal, io->templates, atoi(argv[2])-1),
		  &t, sizeof(t), GT_Templates);

    if (err) {
	Tcl_ResetResult(interp);
    } else {
	Tcl_SetObjResult(interp, GTemplates_klist(interp, io, &t));
    }

    return TCL_OK;
}

int tcl_read_clone(ClientData clientData, Tcl_Interp *interp,
		   int argc, char **argv) {
    GClones c;
    GapIO *io;
    int handle;
    int err;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    err = GT_Read(io, arr(GCardinal, io->clones, atoi(argv[2])-1),
		  &c, sizeof(c), GT_Clones);

    if (err) {
	Tcl_ResetResult(interp);
    } else {
	Tcl_SetObjResult(interp, GClones_klist(interp, io, &c));
    }

    return TCL_OK;
}

int tcl_write_annotation(ClientData clientData, Tcl_Interp *interp,
			 int argc, char **argv) {
    GAnnotations a;
    GapIO *io;
    int handle;
    int err;
    Tcl_Obj *o;

    if (argc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number keyedlist\"\n", *argv);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    o = Tcl_NewStringObj(argv[3], -1);
    klist_GAnnotations(interp, io, &a, o);
    err = GT_Write(io, arr(GCardinal, io->annotations, atoi(argv[2])-1),
		   &a, sizeof(a), GT_Annotations);

    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}

int tcl_write_note(ClientData clientData, Tcl_Interp *interp,
		   int argc, char **argv) {
    GNotes n;
    GapIO *io;
    int handle;
    int err;
    Tcl_Obj *o;

    if (argc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number keyedlist\"\n", *argv);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    o = Tcl_NewStringObj(argv[3], -1);
    klist_GNotes(interp, io, &n, o);
    err = GT_Write(io, arr(GCardinal, io->notes, atoi(argv[2])-1),
		   &n, sizeof(n), GT_Notes);

    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}

int tcl_write_vector(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv) {
    GVectors v;
    GapIO *io;
    int handle;
    int err;
    Tcl_Obj *o;

    if (argc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number keyedlist\"\n", *argv);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    o = Tcl_NewStringObj(argv[3], -1);
    klist_GVectors(interp, io, &v, o);
    err = GT_Write(io, arr(GCardinal, io->vectors, atoi(argv[2])-1),
		   &v, sizeof(v), GT_Vectors);

    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}

int tcl_write_reading(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    GReadings r;
    GapIO *io;
    int handle;
    int err, num;
    Tcl_Obj *o;

    if (argc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number keyedlist\"\n", *argv);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    num = atoi(argv[2]);
    o = Tcl_NewStringObj(argv[3], -1);
    klist_GReadings(interp, io, &r, o);
    err = gel_write(io, num, r);

    /* Keep IO structure up to date */
    io_relpos(io, num) = r.position;
    io_length(io, num) = r.sequence_length * (r.sense ? -1 : 1);
    io_lnbr  (io, num) = r.left;
    io_rnbr  (io, num) = r.right;

    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}

int tcl_write_contig(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv) {
    GContigs c;
    GapIO *io;
    int handle;
    int err, num;
    Tcl_Obj *o;

    if (argc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number keyedlist\"\n", *argv);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    num = atoi(argv[2]);
    o = Tcl_NewStringObj(argv[3], -1);
    klist_GContigs(interp, io, &c, o);
    err = GT_Write(io, arr(GCardinal, io->contigs, num-1),
		   &c, sizeof(c), GT_Contigs);

    /* Keep IO structure up to date */
    io_clength(io, num) = c.length;
    io_clnbr  (io, num) = c.left;
    io_crnbr  (io, num) = c.right;

    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}

int tcl_write_database(ClientData clientData, Tcl_Interp *interp,
		       int argc, char **argv) {
    GapIO *io;
    int handle;
    int err;
    Tcl_Obj *o;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io keyedlist\"\n", *argv);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    o = Tcl_NewStringObj(argv[2], -1);
    klist_GDatabase(interp, io, &io->db, o);
    err = DBDelayWrite(io);

    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}

int tcl_write_template(ClientData clientData, Tcl_Interp *interp,
		       int argc, char **argv) {
    GTemplates t;
    GapIO *io;
    int handle;
    int err;
    Tcl_Obj *o;

    if (argc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number keyedlist\"\n", *argv);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    o = Tcl_NewStringObj(argv[3], -1);
    klist_GTemplates(interp, io, &t, o);
    err = GT_Write(io, arr(GCardinal, io->templates, atoi(argv[2])-1),
		   &t, sizeof(t), GT_Templates);

    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}

int tcl_write_clone(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {
    GClones c;
    GapIO *io;
    int handle;
    int err;
    Tcl_Obj *o;

    if (argc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number keyedlist\"\n", *argv);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    o = Tcl_NewStringObj(argv[3], -1);
    klist_GClones(interp, io, &c, o);
    err = GT_Write(io, arr(GCardinal, io->clones, atoi(argv[2])-1),
		   &c, sizeof(c), GT_Clones);

    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}



/*
 *-----------------------------------------------------------------------------
 * Tcl interfaces to reading and writing reading names. This is different
 * from the generic TextRead interface as it makes use of name caching.
 *-----------------------------------------------------------------------------
 */
int tcl_read_reading_name(ClientData clientData, Tcl_Interp *interp,
			     int argc, char **argv) {
    GapIO *io;
    int handle, rnum;
    char *name;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    rnum = atoi(argv[2]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    name = io_rname(io, rnum);
    Tcl_SetResult(interp, name, TCL_STATIC);
    return TCL_OK;
}

int tcl_write_reading_name(ClientData clientData, Tcl_Interp *interp,
			   int argc, char **argv) {
    GapIO *io;
    int handle, rnum;

    if (argc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io number name\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    rnum = atoi(argv[2]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    write_rname(io, rnum, argv[3]);

    if (auto_flush) flush2t(io);

    Tcl_SetResult(interp, argv[3], TCL_VOLATILE);
    return TCL_OK;
}


/*
 *-----------------------------------------------------------------------------
 * Tcl interfaces to the C {Text,Data,Array}{Read,Write} functions
 *-----------------------------------------------------------------------------
 */
int tcl_io_read_text(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv) {
    GapIO *io;
    int handle, record;
    char *buf;
    char *cp;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io record\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);
    record = atoi(argv[2]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

/*
    printf("Length is %d\n",
	   arr(View, io->client->local.server->local.gdb->view,
	       arr(GView, io->views, record)).cache->used);
*/

    /* Buf will be NULL terminated by TextAllocRead */
    buf = TextAllocRead(io, record);
    if (NULL == buf)
	Tcl_ResetResult(interp);
    else {
	cp = &buf[strlen(buf)];

	/* Convert fortran space-padded strings to C ones */
	if (cp != buf /* strlen != 0 */) {
	    for (cp--; cp >= buf && *cp == ' '; cp--)
		;
	    *(cp+1) = 0;
	}
    
	Tcl_SetResult(interp, buf, TCL_VOLATILE);
	free(buf);
    }
    return TCL_OK;
}


int tcl_io_write_text(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    GapIO *io;
    int handle, record, err;
    
    if (argc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io record string\"\n", *argv);
	return TCL_ERROR;
    }

    handle = atoi(argv[1]);
    record = atoi(argv[2]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    err = TextWrite(io, record, argv[3], strlen(argv[3]));
    
    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}

int tcl_io_read_data(ClientData clientData, Tcl_Interp *interp,
		     int objc, Tcl_Obj *CONST objv[]) {
    GapIO *io;
    int handle, record, len, size;
    char *buf;
    Tcl_Obj *obj;
    int err;
    
    if (objc != 5) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io record numbytes datasize\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetIntFromObj(interp, objv[1], &handle);
    Tcl_GetIntFromObj(interp, objv[2], &record);
    Tcl_GetIntFromObj(interp, objv[3], &len);
    Tcl_GetIntFromObj(interp, objv[4], &size);
    
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }
    
    buf = Tcl_Alloc(len+1);
    err = DataRead(io, record, buf, len, size);

    if (err)
	Tcl_ResetResult(interp);
    else {
	obj = Tcl_NewStringObj(buf, len);
	Tcl_SetObjResult(interp, obj);
    }

    Tcl_Free(buf);

    return TCL_OK;
}


int tcl_io_write_data(ClientData clientData, Tcl_Interp *interp,
		      int objc, Tcl_Obj *CONST objv[]) {

    GapIO *io;
    int handle, record, len, size;
    int err;
    
    if (objc != 6) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io record numbytes datasize datastring\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetIntFromObj(interp, objv[1], &handle);
    Tcl_GetIntFromObj(interp, objv[2], &record);
    Tcl_GetIntFromObj(interp, objv[3], &len);
    Tcl_GetIntFromObj(interp, objv[4], &size);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    err = DataWrite(io, record, Tcl_GetStringFromObj(objv[5], NULL),
		    len, size);
    
    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}


int tcl_io_read_array(ClientData clientData, Tcl_Interp *interp,
		      int objc, Tcl_Obj *CONST objv[]) {
    GapIO *io;
    int handle, record, len;
    Array arr;
    Tcl_Obj *obj;
    
    if (objc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io record numelements\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetIntFromObj(interp, objv[1], &handle);
    Tcl_GetIntFromObj(interp, objv[2], &record);
    Tcl_GetIntFromObj(interp, objv[3], &len);
    
    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }
    
    arr = ArrayRead(io, record, len);

    if (!arr)
	Tcl_ResetResult(interp);
    else {
	int i;
	Tcl_Obj **o;

	o = calloc(len, sizeof(*o));
	for (i = 0; i < len; i++) {
	    o[i] = Tcl_NewIntObj(arr(GCardinal, arr, i));
	}
	obj = Tcl_NewListObj(len, o);
	Tcl_SetObjResult(interp, obj);
    }

    return TCL_OK;
}


int tcl_io_write_array(ClientData clientData, Tcl_Interp *interp,
		       int objc, Tcl_Obj *CONST objv[]) {

    GapIO *io;
    int handle, record, len;
    int err, i;
    Array arr;
    
    if (objc != 4) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io record list\"\n",
		       Tcl_GetStringFromObj(objv[0], NULL));
	return TCL_ERROR;
    }

    Tcl_GetIntFromObj(interp, objv[1], &handle);
    Tcl_GetIntFromObj(interp, objv[2], &record);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    if (TCL_OK != Tcl_ListObjLength(interp, objv[3], &len))
	return TCL_ERROR;

    if (NULL == (arr = ArrayCreate(sizeof(GCardinal), len)))
	return TCL_ERROR;

    for (i = 0; i < len; i++) {
	int ele;
	Tcl_Obj *i_obj;
	if (TCL_OK != Tcl_ListObjIndex(interp, objv[3], i, &i_obj))
	    return TCL_ERROR;

	if (TCL_OK != Tcl_GetIntFromObj(interp, i_obj, &ele))
	    return TCL_ERROR;

	arr(GCardinal, arr, i) = ele;
    }
    
    err = ArrayWrite(io, record, len, arr);
    if (auto_flush) flush2t(io);

    /* Also copy over in-memory copy if the record matches standard ones */
    if (record == io->db.contigs)
	memcpy(ArrayBase(GCardinal, io->contigs),
	       ArrayBase(GCardinal, arr), len * sizeof(GCardinal));
    else if (record == io->db.readings)
	memcpy(ArrayBase(GCardinal, io->readings),
	       ArrayBase(GCardinal, arr), len * sizeof(GCardinal));
    else if (record == io->db.annotations)
	memcpy(ArrayBase(GCardinal, io->annotations),
	       ArrayBase(GCardinal, arr), len * sizeof(GCardinal));
    else if (record == io->db.templates)
	memcpy(ArrayBase(GCardinal, io->templates),
	       ArrayBase(GCardinal, arr), len * sizeof(GCardinal));
    else if (record == io->db.clones)
	memcpy(ArrayBase(GCardinal, io->clones),
	       ArrayBase(GCardinal, arr), len * sizeof(GCardinal));
    else if (record == io->db.vectors)
	memcpy(ArrayBase(GCardinal, io->vectors),
	       ArrayBase(GCardinal, arr), len * sizeof(GCardinal));
    else if (record == io->db.notes_a)
	memcpy(ArrayBase(GCardinal, io->notes),
	       ArrayBase(GCardinal, arr), len * sizeof(GCardinal));
    else if (record == io->db.contig_order)
	memcpy(ArrayBase(GCardinal, io->contig_order),
	       ArrayBase(GCardinal, arr), len * sizeof(GCardinal));

    ArrayDestroy(arr);

    vTcl_SetResult(interp, "%d", err != 0);

    return TCL_OK;
}


/*
 *-----------------------------------------------------------------------------
 * Tcl interfaces to the io_init_* (IO.c) and add_* (clone.c) functions.
 * Also an interface to the allocate routine.
 *-----------------------------------------------------------------------------
 */
int tcl_io_add_reading(ClientData clientData, Tcl_Interp *interp,
		       int argc, char **argv) {
    GapIO *io;
    int handle, err, n;
    GReadings r;

    if (argc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    io_init_reading(io, n = NumReadings(io)+1);

    if (err = gel_read(io, n, r)) {
	GAP_ERROR("GT_Read (gel extend)");
	return TCL_ERROR;
    }

    if (r.name == 0)
	r.name = allocate(io, GT_Text);
    TextWrite(io, r.name, "uninitialised", sizeof("uninitialised"));
    gel_write(io, n, r);
    io_wname(io, n, "uninitialised");
    
    io_write_rd(io, n, 
		"uninitialised", sizeof("uninitialised"),
		"uninitialised", sizeof("uninitialised"));
    if (auto_flush) flush2t(io);

    vTcl_SetResult(interp, "%d", n);

    return TCL_OK;
}

int tcl_io_add_contig(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    GapIO *io;
    int handle;

    if (argc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    io_init_contig(io, NumContigs(io)+1);
    if (auto_flush) flush2t(io);
    vTcl_SetResult(interp, "%d", NumContigs(io));

    return TCL_OK;
}

int tcl_io_add_annotation(ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv) {
    GapIO *io;
    int handle;

    if (argc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    io_init_annotations(io, Nannotations(io)+1);
    if (auto_flush) flush2t(io);
    vTcl_SetResult(interp, "%d", Nannotations(io));

    return TCL_OK;
}

int tcl_io_add_note(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {
    GapIO *io;
    int handle;

    if (argc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    io_init_note(io, Nnotes(io)+1);
    if (auto_flush) flush2t(io);
    vTcl_SetResult(interp, "%d", Nnotes(io));

    return TCL_OK;
}

/*
 * Find a vector and use as a default.
 *
 * If we do not supply a valid vector for the add_*() routines then an
 * extra vector will automatically be created. We do not wish this to be
 * so.
 */
static char *default_vector(GapIO *io) {
    static char buf[1024];
    GVectors v;

    if (Nvectors(io) == 0 ||
        GT_Read(io, arr(GCardinal, io->vectors, 0),
		&v, sizeof(v), GT_Vectors) != 0 ||
        TextRead(io, v.name, buf, 1023) != 0)
	*buf = 0;

    return buf;
}

int tcl_io_add_template(ClientData clientData, Tcl_Interp *interp,
			int argc, char **argv) {
    GapIO *io;
    int handle;

    if (argc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    add_template(io, "uninitialised", default_vector(io), "0", "0..0", 0);
    if (auto_flush) flush2t(io);
    vTcl_SetResult(interp, "%d", Ntemplates(io));

    return TCL_OK;
}

int tcl_io_add_clone(ClientData clientData, Tcl_Interp *interp,
		     int argc, char **argv) {
    GapIO *io;
    int handle;

    if (argc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    add_clone(io, "uninitialised", default_vector(io));
    if (auto_flush) flush2t(io);
    vTcl_SetResult(interp, "%d", Nclones(io));

    return TCL_OK;
}

int tcl_io_add_vector(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    GapIO *io;
    int handle;

    if (argc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    add_vector(io, "uninitialised", 0);
    if (auto_flush) flush2t(io);
    vTcl_SetResult(interp, "%d", Nvectors(io));

    return TCL_OK;
}

int tcl_io_allocate(ClientData clientData, Tcl_Interp *interp,
		    int argc, char **argv) {
    GapIO *io;
    int handle, type;
    int rec;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		"\"%s io type\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    if (strcmp(argv[2], "text") == 0) {
	type = GT_Text;
    } else {
	Tcl_SetResult(interp, "Only \"text\" type supported at present\n",
		      TCL_STATIC);
	return TCL_ERROR;
    }

    if (auto_flush) flush2t(io);

    rec = allocate(io, type);
    GT_Write(io, rec, NULL, 0, type);
    vTcl_SetResult(interp, "%d", rec);

    return TCL_OK;
}

int tcl_io_deallocate(ClientData clientData, Tcl_Interp *interp,
		      int argc, char **argv) {
    GapIO *io;
    int handle;
    int rec;

    if (argc != 3) {
	vTcl_SetResult(interp, "wrong # args: should be "
		"\"%s io record\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    rec = atoi(argv[2]);

    vTcl_SetResult(interp, "%d", deallocate(io, rec));

    if (auto_flush) flush2t(io);

    return TCL_OK;
}

/*
 *-----------------------------------------------------------------------------
 * Flushing - not done automatically
 *-----------------------------------------------------------------------------
 */
int tcl_io_flush(ClientData clientData, Tcl_Interp *interp,
		 int argc, char **argv) {
    GapIO *io;
    int handle;

    if (argc != 2) {
	vTcl_SetResult(interp, "wrong # args: should be "
		       "\"%s io\"\n", *argv);
	return TCL_ERROR;
    }
    
    handle = atoi(argv[1]);

    if (NULL == (io = io_handle(&handle))) {
	Tcl_SetResult(interp, "Invalid io handle\n", TCL_STATIC);
	return TCL_ERROR;
    }

    flush2t(io);
    
    return TCL_OK;
}


/*
 *-----------------------------------------------------------------------------
 * Tcl initialisation - register our commands
 *-----------------------------------------------------------------------------
 */
int Db_Init(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "io_read_annotation",
		      tcl_read_annotation, NULL, NULL);
    Tcl_CreateCommand(interp, "io_read_note",
		      tcl_read_note, NULL, NULL);
    Tcl_CreateCommand(interp, "io_read_vector",
		      tcl_read_vector, NULL, NULL);
    Tcl_CreateCommand(interp, "io_read_reading",
		      tcl_read_reading, NULL, NULL);
    Tcl_CreateCommand(interp, "io_read_contig",
		      tcl_read_contig, NULL, NULL);
    Tcl_CreateCommand(interp, "io_read_database",
		      tcl_read_database, NULL, NULL);
    Tcl_CreateCommand(interp, "io_read_template",
		      tcl_read_template, NULL, NULL);
    Tcl_CreateCommand(interp, "io_read_clone",
		      tcl_read_clone, NULL, NULL);

    Tcl_CreateCommand(interp, "io_write_annotation",
		      tcl_write_annotation, NULL, NULL);
    Tcl_CreateCommand(interp, "io_write_note",
		      tcl_write_note, NULL, NULL);
    Tcl_CreateCommand(interp, "io_write_vector",
		      tcl_write_vector, NULL, NULL);
    Tcl_CreateCommand(interp, "io_write_reading",
		      tcl_write_reading, NULL, NULL);
    Tcl_CreateCommand(interp, "io_write_contig",
		      tcl_write_contig, NULL, NULL);
    Tcl_CreateCommand(interp, "io_write_database",
		      tcl_write_database, NULL, NULL);
    Tcl_CreateCommand(interp, "io_write_template",
		      tcl_write_template, NULL, NULL);
    Tcl_CreateCommand(interp, "io_write_clone",
		      tcl_write_clone, NULL, NULL);

    Tcl_CreateCommand(interp, "io_read_reading_name",
		      tcl_read_reading_name, NULL, NULL);
    Tcl_CreateCommand(interp, "io_write_reading_name",
		      tcl_write_reading_name, NULL, NULL);

    Tcl_CreateCommand(interp, "io_read_text", tcl_io_read_text, NULL, NULL);
    Tcl_CreateCommand(interp, "io_write_text", tcl_io_write_text, NULL, NULL);
    Tcl_CreateObjCommand(interp, "io_read_data", tcl_io_read_data, NULL, NULL);
    Tcl_CreateObjCommand(interp, "io_write_data", tcl_io_write_data, NULL, NULL);
    Tcl_CreateObjCommand(interp, "io_read_array", tcl_io_read_array, NULL, NULL);
    Tcl_CreateObjCommand(interp, "io_write_array", tcl_io_write_array, NULL, NULL);

    Tcl_CreateCommand(interp, "io_flush", tcl_io_flush, NULL, NULL);
    Tcl_LinkVar(interp, "gap_auto_flush", (char *)&auto_flush,
		TCL_LINK_BOOLEAN);

    Tcl_CreateCommand(interp, "io_add_reading",
		      tcl_io_add_reading, NULL, NULL);
    Tcl_CreateCommand(interp, "io_add_annotation",
		      tcl_io_add_annotation, NULL, NULL);
    Tcl_CreateCommand(interp, "io_add_note",
		      tcl_io_add_note, NULL, NULL);
    Tcl_CreateCommand(interp, "io_add_contig",
		      tcl_io_add_contig, NULL, NULL);
    Tcl_CreateCommand(interp, "io_add_clone",
		      tcl_io_add_clone, NULL, NULL);
    Tcl_CreateCommand(interp, "io_add_template",
		      tcl_io_add_template, NULL, NULL);
    Tcl_CreateCommand(interp, "io_add_vector",
		      tcl_io_add_vector, NULL, NULL);
    Tcl_CreateCommand(interp, "io_allocate",
		      tcl_io_allocate, NULL, NULL);
    Tcl_CreateCommand(interp, "io_deallocate",
		      tcl_io_deallocate, NULL, NULL);

    return TCL_OK;
}
