#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "os.h"
#include "IO.h"
#include "xalloc.h"
#include "copy_db.h"
#include "text_output.h"
#include "misc.h"

/*
 * A bunch of useful defines to tidy up the following function.
 * Some of these could be functions, but others would require a lot of
 * extra bits and so it's not really worth it.
 */

#define COPY_STR(iof, iot, rec_from, rec_to, buf, len) \
if (rec_from > 0 && rec_from <= iof->Nviews) { \
    if (TextRead(iof, rec_from, buf, len) && errs) \
	goto error; \
    if (-1 == (rec_to = allocate(iot, GT_Text))) \
	goto error; \
    if (TextWrite(iot, rec_to, buf, len) && errs) \
	goto error; \
} else { \
    rec_to = 0; \
}

#if 0
#define COPY_ALLOC_STR(iof, iot, rec_from, rec_to) \
if (rec_from > 0 && rec_from <= iof->Nviews) { \
    char *str = TextAllocRead(iof, rec_from); \
    if (str == NULL || -1 == (rec_to = allocate(iot, GT_Text)) && errs) \
	goto error; \
    if (TextWrite(iot, rec_to, str, strlen(str)) && errs) \
	goto error; \
    xfree(str); \
} else { \
    rec_to = 0; \
}
#endif
#define COPY_ALLOC_STR(iof, iot, rec_from, rec_to) \
if (rec_from && rec_from <= iof->Nviews) { \
    char *str = TextAllocRead(iof, rec_from); \
    if (str == NULL || -1 == (rec_to = allocate(iot, GT_Text)) && errs) { \
	printf("fail on record %d\n", rec_from); rec_to = 0; \
    } else if (TextWrite(iot, rec_to, str, strlen(str)) && errs) \
	goto error; \
    xfree(str); \
} else { \
    rec_to = 0; \
}


#define COPY_DATA(iof, iot, rec_from, rec_to, buf, len, size) \
if (rec_from && rec_from <= iof->Nviews) { \
    if (DataRead(iof, rec_from, buf, len*size, size) && errs) \
	goto error; \
    if (-1 == (rec_to = allocate(iot, GT_Data)) && errs) \
	goto error; \
    if (DataWrite(iot, rec_to, buf, len*size, size) && errs) \
	goto error; \
} else { \
    rec_to = 0; \
}

/* Update the GDatabase structure and one of the io arrays */
#define WRITE_ARRAY(iot, array_name) \
do { \
    if (GT_Write(iot, GR_Database, &iot->db, sizeof(iot->db), GT_Database) && errs) \
	goto error; \
    if (ArrayWrite(iot, iot->db.array_name, N##array_name(iot), iot->array_name) && errs) \
	goto error; \
} while(0)

#define VERB_NUM(num) \
if (verbose) { \
    printf("%.*s",back, "\b\b\b\b\b\b"); \
    back = printf("%d", i); \
    fflush(stdout); \
}

#define NEW_RECORD(iot, type, array, index) \
do { \
    int freerec; \
\
    if (-1 == (freerec = allocate(iot, type)) && errs) \
	goto error; \
    arr(GCardinal, array, index-1) = freerec; \
} while (0);


/*
 * Copies the reading name for sequence rnum from DB iof to DB iot. Normally
 * the name is the same, but it may be "unknown.%d" or "oldname#%d" if the
 * sequence name is duplicated.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int copy_read_name(GapIO *iof, GapIO *iot, int rnum, int start_r,
		   GReadings *rf, GReadings *rt, int *num_unknowns) {
    char buf[DB_NAMELEN+1];

    /* Find the name we are going to use */
    if (rf->name) {
	TextRead(iof, rf->name, buf, DB_NAMELEN);

	/* If it's already in iot, generate a new name */
	if (get_gel_num(iot, buf, GGN_NAME) != -1) {
	    int iter = 1;
	    size_t l = strlen(buf);
	    char new_name[DB_NAMELEN+1];

	    do {
		char num[10];
		sprintf(num, "#%d", iter++);
		sprintf(new_name, "%.*s%s",
			(int)(MIN(l, DB_NAMELEN-strlen(num))),
			buf, num);
	    } while (get_gel_num(iot, new_name, GGN_NAME) != -1);
	    
	    printf("Fixed duplicate reading %s, given new name %s\n",
		    buf, new_name);

	    strcpy(buf, new_name);
	}
    } else {
	/* Name does not exist, so make one up */
	do {
	    sprintf(buf, "unknown.%d", (*num_unknowns)++);
	} while (get_gel_num(iot, buf, GGN_NAME) != -1);

	printf("Fixed unknown reading name for #%d (now %s)\n",
	       rnum, buf);
    }

    /* Allocate and copy it */
    if (-1 == (rt->name = allocate(iot, GT_Text)))
	return -1;

    TextWrite(iot, rt->name, buf, DB_NAMELEN+1);
    io_wname(iot, rnum + start_r, buf);

    return 0;
}

/*
 * Copies database pointed to by 'iof' to the database pointed to by 'iot'.
 *
 * If 'verbose' is set then a running summary of its actions is displayed.
 */
int copy_database(GapIO *iof, GapIO *iot, int verbose, int errs) {
    GReadings r, rt;
    GContigs c;
    GAnnotations a, at;
    GNotes n, nt;
    GTemplates t, tt;
    GClones s, st;
    GVectors v, vt;
    int i, back, old;
    int start_c, start_r, start_v, start_s, start_a, start_t, start_f;
    int start_n;
    int num_unknowns = 0;
    int max_len = find_max_gel_len(iof, 0, 0)+1;
    char *seq = NULL;
    int1 *conf = NULL;
    int2 *opos = NULL;
    int retcode = -1;

    extern int gap_fatal_errors;
    old = gap_fatal_errors;
    gap_fatal_errors = 0;

    start_c = NumContigs(iot);
    start_r = NumReadings(iot);
    start_v = Nvectors(iot);
    start_s = Nclones(iot);
    start_a = Nannotations(iot);
    start_n = Nnotes(iot);
    start_t = Ntemplates(iot);
    start_f = BIT_NBITS(iot->freerecs);


    /*
     * Write the GReadings structures.
     */
    if (verbose) (printf("Copying reading "), back = 0);
    fflush(stdout);
    if (-1 == io_init_reading(iot, NumReadings(iof)+start_r) && errs)
	goto error;

    if (NULL == (seq = (char *)xmalloc(max_len * sizeof(*seq))))
	goto error;
    if (NULL == (conf = (int1 *)xmalloc(max_len * sizeof(*conf))))
	goto error;
    if (NULL == (opos = (int2 *)xmalloc(max_len * sizeof(*opos))))
	goto error;

    for (i=1; i<=NumReadings(iof); i++) {
	char t_type[5];

	VERB_NUM(i);
	if (gel_read(iof, i, r) && errs) goto error;
	memcpy(&rt, &r, sizeof(r));

	rt.sequence_length = rt.end - rt.start - 1;
	if (rt.annotations > iof->db.Nannotations || rt.annotations < 0) {
	    printf("Fixed annotations pointer for reading %d\n", i);
	    rt.annotations = 0;
	}
	if (rt.notes > iof->db.Nnotes || rt.notes < 0) {
	    printf("Fixed notes pointer for reading %d\n", i);
	    rt.notes = 0;
	}
	if (rt.template > iof->db.Ntemplates || rt.template < 0) {
	    printf("Fixed template pointer for reading %d\n", i);
	    rt.template = 0;
	}
	if (rt.left > iof->db.num_readings || rt.left < 0) {
	    printf("Fixed left pointer for reading %d\n", i);
	    rt.left = 0;
	}
	if (rt.right > iof->db.num_readings || rt.right < 0) {
	    printf("Fixed right pointer for reading %d\n", i);
	    rt.right = 0;
	}

	COPY_STR(iof, iot, r.sequence, rt.sequence, seq, r.length);
	COPY_DATA(iof, iot, r.confidence, rt.confidence, conf,
		  r.length, sizeof(int1));
	COPY_DATA(iof, iot, r.orig_positions, rt.orig_positions, opos,
		  r.length, sizeof(int2));
	if (copy_read_name(iof, iot, i, start_r, &r, &rt, &num_unknowns) == -1)
	    goto error;
	COPY_ALLOC_STR(iof, iot, r.trace_name, rt.trace_name);
	COPY_STR(iof, iot, r.trace_type, rt.trace_type, t_type, 5);

	if (rt.strand < 0 || rt.strand > 1) {
	    printf("Fixed strand for reading %d\n", i);
	    rt.strand = 0;
	}
	if (rt.primer < 0 || rt.primer > 4) {
	    printf("Fixed primer for reading %d\n", i);
	    rt.primer = 0;
	}

	if (rt.left)
	    rt.left += start_r;
	if (rt.right)
	    rt.right += start_r;
	if (rt.template)
	    rt.template += start_t;
	if (rt.annotations)
	    rt.annotations += start_a;
	if (rt.notes)
	    rt.notes += start_n;
	
	if (gel_write(iot, i + start_r, rt) && errs) goto error;
    }
    if (verbose) printf("\n");


    /*
     * Write the GContig structures
     */
    if (verbose) (printf("Copying contig "), back = 0);
    fflush(stdout);
    if (io_init_contig(iot, NumContigs(iof)+start_c) && errs) goto error;
    for (i=1; i<=NumContigs(iof); i++) {
	VERB_NUM(i);
	if (contig_read(iof, i, c) && errs) goto error;
	c.left += start_r;
	c.right += start_r;
	if (c.annotations > iof->db.Nannotations || c.annotations < 0) {
	    printf("Fixed annotations pointer for contig %d\n", i);
	    c.annotations = 0;
	}
	if (c.notes > iof->db.Nnotes || c.notes < 0) {
	    printf("Fixed notes pointer for contig %d\n", i);
	    c.notes = 0;
	}
	if (c.annotations)
	    c.annotations += start_a;
	if (c.notes)
	    c.notes += start_n;
	if (contig_write(iot, i + start_c, c) && errs) goto error;
    }
    if (verbose) printf("\n");


    /*
     * Write the GAnnotations structures
     *
     * NOTE: An optimisation would be to only write out used annotations.
     * However this requires renumbering of annotations as those on the free
     * list are not necessary the end numbers. This in turn requires a lot
     * of annotation relinking and work with reading and contig structures.
     */
    if (verbose) (printf("Copying annotation "), back = 0);
    fflush(stdout);
    if (io_init_annotations(iot, Nannotations(iof)+start_a) && errs)
	goto error;
    for (i=1; i<=Nannotations(iof); i++) {
	VERB_NUM(i);
	if (tag_read(iof, i, a) && errs) goto error;
	memcpy(&at, &a, sizeof(a));
	if (errs==0) {
	    if (a.annotation >= 0 && a.annotation < iof->Nviews) {
		COPY_ALLOC_STR(iof, iot, a.annotation, at.annotation);
	    } else {
		printf("Fixed a.annotation (%d) value for anno %d\n",
		       a.annotation, i);
		at.annotation = 0;
	    }
	} else {
	    COPY_ALLOC_STR(iof, iot, a.annotation, at.annotation);
	}
	if (at.next> iof->db.Nannotations) {
	    printf("Fixed next pointer for annotation %d\n", i);
	    at.next = 0;
	}
	if (at.next)
	    at.next += start_a;
	if (tag_write(iot, i + start_a, at) && errs) goto error;
    }
    /* Join free lists */
    if (iot->db.free_annotations) {
	int anno;

	if (anno = iof->db.free_annotations) {
	    for (;;) {
		tag_read(iof, anno, a);
		if (!a.next)
		    break;
		anno = a.next;
	    }
	    a.next = iot->db.free_annotations;
	    tag_write(iot, anno + start_a, a);
	    iot->db.free_annotations = iof->db.free_annotations + start_a;
	}
    } else {
	if (iof->db.free_annotations)
	    iot->db.free_annotations = iof->db.free_annotations + start_a;
    }
    if (verbose) printf("\n");


    /*
     * Write the GNotes structures
     *
     * NOTE: An optimisation would be to only write out used notes.
     * However this requires renumbering of notes as those on the free
     * list are not necessary the end numbers. This in turn requires a lot
     * of note relinking and work with reading, contig and database
     * structures.
     */
    if (verbose) (printf("Copying note "), back = 0);
    fflush(stdout);
    if (io_init_note(iot, Nnotes(iof)+start_n) && errs)
	goto error;
    for (i=1; i<=Nnotes(iof); i++) {
	VERB_NUM(i);
	if (note_read(iof, i, n) && errs) goto error;
	memcpy(&nt, &n, sizeof(n));
	if (errs==0) {
	    if (n.annotation >= 0 && n.annotation < iof->Nviews) {
		COPY_ALLOC_STR(iof, iot, n.annotation, nt.annotation);
	    } else {
		printf("Fixed n.annotation (%d) value for note %d\n",
		       n.annotation, i);
		nt.annotation = 0;
	    }
	} else {
	    COPY_ALLOC_STR(iof, iot, n.annotation, nt.annotation);
	}
	if (nt.next > iof->db.Nnotes || nt.next < 0) {
	    printf("Fixed next pointer for note %d\n", i);
	    nt.next = 0;
	}
	if (nt.next)
	    nt.next += start_n;
	if (nt.prev_type == GT_Notes) {
	    if (nt.prev)
		nt.prev += start_n;
	} else if (nt.prev_type == GT_Readings) {
	    nt.prev += start_r;
	} else if (nt.prev_type == GT_Contigs) {
	    nt.prev += start_c;
	} /* else nt.prev_type == GT_Database or 0 (freelist) */

	if (note_write(iot, i + start_n, nt) && errs) goto error;
    }
    /* Join free lists */
    if (iot->db.free_notes && iof->db.free_notes) {
	int note, next;

	/* Find last note in free list */
	for (note = iot->db.free_notes;
	     note_read(iot, note, n), n.next;
	     note = n.next)
	    ;

	/* And then join forward... */
	next = n.next = iof->db.free_notes + start_n;
	note_write(iot, note, n);

	/* ...and backward */
	note_read(iot, next, n);
	n.prev = note;
	n.prev_type = GT_Notes;
	note_write(iot, next, n);
    } else {
	/* No notes in "io_to", so copy value from "io_from". */
	if (iof->db.free_notes)
	    iot->db.free_notes = iof->db.free_notes + start_n;
    }
    if (verbose) printf("\n");


    /*
     * Write the GVectors structures
     */
    if (verbose) (printf("Copying vector "), back = 0);
    fflush(stdout);
    /* Extend arrays */
    Nvectors(iot) = Nvectors(iof) + start_v;
    ArrayRef(iot->vectors, Nvectors(iot)-1);
    for (i=1; i<=Nvectors(iof); i++) {
	VERB_NUM(i);
	if (vector_read(iof, i, v) && errs) goto error;
	memcpy(&vt, &v, sizeof(v));
	NEW_RECORD(iot, GT_Vectors, iot->vectors, i + start_v);
	COPY_ALLOC_STR(iof, iot, v.name, vt.name);
	if (vector_write(iot, i + start_v, vt) && errs) goto error;
    }
    WRITE_ARRAY(iot, vectors);
    if (verbose) printf("\n");


    /*
     * Write the GTemplates structures
     */
    if (verbose) (printf("Copying template "), back = 0);
    fflush(stdout);
    /* Extend arrays */
    Ntemplates(iot) = Ntemplates(iof) + start_t;
    ArrayRef(iot->templates, Ntemplates(iot)-1);
    for (i=1; i<=Ntemplates(iof); i++) {
	VERB_NUM(i);
	if (template_read(iof, i, t) && errs) goto error;
	memcpy(&tt, &t, sizeof(t));
	NEW_RECORD(iot, GT_Templates, iot->templates, i + start_t);
	COPY_ALLOC_STR(iof, iot, t.name, tt.name);
	if (tt.clone > iof->db.Nclones || tt.clone < 0) {
	    printf("Fixed clone pointer for template %d\n", i);
	    tt.clone = 0;
	}
	if (tt.vector > iof->db.Nvectors || tt.vector < 0) {
	    printf("Fixed vector pointer for template %d\n", i);
	    tt.vector = 0;
	}
	if (tt.insert_length_min > tt.insert_length_max) {
	    int max = tt.insert_length_max;
	    tt.insert_length_max = tt.insert_length_min;
	    tt.insert_length_min = max;
	    printf("Fixed insert size for template %d\n", i);
	}
	if (tt.clone)
	    tt.clone += start_s;
	if (tt.vector)
	    tt.vector += start_v;
	if (template_write(iot, i + start_t, tt) && errs) goto error;
    }
    WRITE_ARRAY(iot, templates);
    if (verbose) printf("\n");


    /*
     * Write the GClones structures
     */
    if (verbose) (printf("Copying clone "), back = 0);
    fflush(stdout);
    /* Extend arrays */
    Nclones(iot) = Nclones(iof) + start_s;
    ArrayRef(iot->clones, Nclones(iot)-1);
    for (i=1; i<=Nclones(iof); i++) {
	VERB_NUM(i);
	if (clone_read(iof, i, s) && errs) goto error;
	memcpy(&st, &s, sizeof(s));
	NEW_RECORD(iot, GT_Clones, iot->clones, i + start_s);
	COPY_ALLOC_STR(iof, iot, s.name, st.name);
	if (st.vector > iof->db.Nvectors || st.vector < 0) {
	    printf("Fixed vector pointer for clone %d\n", i);
	    st.vector = 0;
	}
	if (st.vector)
	    st.vector += start_v;
	if (clone_write(iot, i + start_s, st) && errs) goto error;
    }
    WRITE_ARRAY(iot, clones);
    if (verbose) printf("\n");

    /*
     * Write the freerecs array
     */
    if (verbose) printf("Copying free records bitmap\n");
    if (BitmapWrite(iot, iot->db.freerecs, iot->freerecs) && errs)
	goto error;

    /*
     * Write the contig order
     */
    if (iof->db.contig_order) {
	int i;

	if (verbose) printf("Copying contig order\n");
	ArrayRef(iot->contig_order, start_c + iof->db.Ncontigs);
	for (i = 0; i < iof->db.num_contigs; i++) {
	    arr(GCardinal, iot->contig_order, i+start_c) =
		start_c + arr(GCardinal, iof->contig_order, i);
	}
	ArrayWrite(iot, iot->db.contig_order, iot->db.Ncontigs,
		   iot->contig_order);
    }

    /*
     * Update the database notes
     */
    if (iof->db.notes) {
	int note;

	if (iof->db.notes > iof->db.Nnotes || iof->db.notes < 0) {
	    printf("Fixed notes pointer for database record\n");
	    iof->db.notes = 0;
	}

	for (note = iof->db.notes;;) {
	    note_read(iof, note, n);
	    if (!n.next)
		break;
	    note = n.next;
	}
	note_read(iot, note + start_n, n);
	n.next = iot->db.notes;
	note_write(iot, note + start_n, n);
	if (iot->db.notes) {
	    note_read(iot, iot->db.notes, n);
	    n.prev = note + start_n;
	    n.prev_type = GT_Notes;
	    note_write(iot, iot->db.notes, n);
	}
	iot->db.notes = iof->db.notes + start_n;

	if (GT_Write(iot, GR_Database, &iot->db, sizeof(iot->db), GT_Database))
	    goto error;
    }

    retcode = 0;
 error:
    if (seq)
	xfree(seq);
    if (conf)
	xfree(conf);
    if (opos)
	xfree(opos);
    gap_fatal_errors = old;
    return retcode;
}

/*
 * Interface to copy_database. This is called when we already have a database
 * open and wish to copy it to a new database. Checking for the file existing
 * is already done by the calling routine.
 */
int copy_database_from(GapIO *iof, char *name, char *version) {
    GapIO *iot;
    int status;
    char buf[1024];

    sprintf(buf, "%s.%s", name, version);
    remove(buf);
    sprintf(buf, "%s.%s.aux", name, version);
    remove(buf);
    sprintf(buf, "%s.%s.BUSY", name, version);
    remove(buf);

    if (NULL == (iot = open_db(name, version, &status, 1, 0))) {
	return -1;
    }

    if (-1 == copy_database(iof, iot, 0, 1)) {
	close_db(iot);
	return -1;
    }
    
    close_db(iot);

    /* Close existing log file and reopen original one */
    {
	char log_fn[1024];
	sprintf(log_fn, "%s.log", io_name(iof));

	log_file(log_fn, NULL);
    }

    return 0;
}
