/*
 * Check the database for consistency. Based upon Check Database Consistency
 * code from dbsyscommon.f
 *
 * Created: 31 January 1994 - jkb
 *
 * TODO:
 *     check that the name fields of template, vector, clone can be loaded up
 *
 * NOTES:
 *     Where applicable we use the memory arrays instead of disk. This makes
 *     things easier for hand holding checking etc, but does not check the
 *     integrity of memory vs disk. This is performed as a separate check
 *     for each gel.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "fortran.h"
#include "IO.h"

#include "misc.h"
#include "gap-error.h"
#include "xalloc.h"
#include "dbcheck.h"
#include "gap_globals.h"
#include "text_output.h"

/* A 'level' of vector greater than this is classed as a likely corruption */
#define MAX_LEVEL 10

#define abs(a) ((a) >= 0 ? (a) : -(a))

/*
 * Checks a gel no. is within the correct range.
 * Returns:
 *     0 ok
 *     1 is zero
 *     2 is outside of valid range
 */
static int bad_gel(GapIO *io, int gel) {
    if (gel == 0)
	return 1;
    else if (gel > NumReadings(io) || gel < 0)
	return 2;
    else
	return 0;
}

/*---------------------------------------------------------------------------*/
/* check routines. One per database type (although they are not strictly
 * separate checks and some need to be called in the correct order).
 */

int check_database(GapIO *io, f_int idbsiz, f_int ngels,
		    f_int nconts, int *note_used, int *minor) {
    int status = 0;
    int note;
    GNotes n;

    /* check number of contigs */
    if (NumContigs(io) > Ncontigs(io)) {
	vmessage("Database: more contigs used (%d) than allocated (%d).\n",
	       NumContigs(io), Ncontigs(io));
	status++;
    }
    if (NumContigs(io) != nconts) {
	vmessage("Database: number of contigs used in memory (%d) and disk (%d) "
	       "differ.\n", nconts, NumContigs(io));
	status++;
    }

    /* check number of gels */
    if (NumReadings(io) > Nreadings(io)) {
	vmessage("Database: more readings used (%d) than allocated (%d).\n",
	       NumReadings(io), Nreadings(io));
	status++;
    }
    if (NumReadings(io) != ngels) {
	vmessage("Database: number of readings used in memory (%d) and disk (%d)"
	       " differ.\n", ngels, NumReadings(io));
	status++;
    }

    /* database size */
    if (io->db.actual_db_size != idbsiz) {
	vmessage("Database: database size in memory (%d) and disk (%d) "
	       "differ.\n", idbsiz, io->db.actual_db_size);
	status++;
    }
    if (io->db.actual_db_size > io->db.maximum_db_size) {
	vmessage("Database: actual database size (%d) is greater than the "
	       "maximum (%d).\n",
	       io->db.actual_db_size, io->db.maximum_db_size);
	status++;
    }

    /* other valid fields */
    if (io->db.data_class < 0 || io->db.data_class > 1) {
	vmessage("Database: data_class (%d) is not 0 or 1.\n",
	       io->db.data_class);
	status++;
    }
    if (io->db.free_annotations < 0 ||
	io->db.free_annotations > Nannotations(io)) {
	vmessage("Database: invalid number of free annotations (%d).\n"
	       "          Total number of annotations = %d.\n",
	       io->db.free_annotations, Nannotations(io));
	status++;
    }
    if (io->db.free_notes < 0 ||
	io->db.free_notes > Nnotes(io)) {
	vmessage("Database: invalid number of free notes (%d).\n"
	       "          Total number of notes = %d.\n",
	       io->db.free_notes, Nnotes(io));
	status++;
    }

    /* Fill in note_used[] array */
    note = io->db.notes;
    if (note) {
	note_read(io, note, n);
	/* Check prev/prev_type of 1st item */
	if (n.prev_type != GT_Database ||
	    n.prev != 0) {
	    vmessage("Database note %d links back to prev=%d prev_type=%d\n",
		     note, n.prev, n.prev_type);
	    status++;
	}
    }
    while (note) {
	if (note_used[note]) {
	    vmessage("Database note %d used more than once (loop?).\n", note);
	    status++;
	    break;
	}
	
	note_used[note]++;
	if (note = n.next)
	    note_read(io, note, n);
    }

    return status;
}

int check_contigs(GapIO *io, f_int *relpg, f_int *lngthg, f_int *lnbr,
		  f_int *rnbr, int *gel_used, int *tag_used,
		  int *note_used, int *minor) {
    int contig, gel, len, loop = 0, last;
    GContigs c;
    GAnnotations a;
    int anno;
    int status = 0;
    int last_pos, last_anno;
    int note;
    GNotes n;

    for (contig = 1; contig <= NumContigs(io); contig++) {
	(void)GT_Read(io, arr(GCardinal, io->contigs, contig - 1),
			  &c, sizeof(c), GT_Contigs);

	/* check memory vs disk */
	if (io_clnbr(io, contig) != c.left) {
	    vmessage("Contig %d: Memory left = %d, disk left = %d.\n",
		     contig, io_clnbr(io, contig), c.left);
	    status++;
	}
	if (io_crnbr(io, contig) != c.right) {
	    vmessage("Contig %d: Memory right = %d, disk right = %d.\n",
		     contig, io_crnbr(io, contig), c.right);
	    status++;
	}
	if (io_clength(io, contig) != c.length) {
	    vmessage("Contig %d: Memory length = %d, disk length = %d.\n",
		     contig, io_clength(io, contig), c.length);
	    status++;
	}

	/*
	 * check left && right end
	 */
	if (c.left == 0) {
	    vmessage("Contig %d: no left gel number.\n", contig);
	    status++;
	}
	if (c.right == 0) {
	    vmessage("Contig %d: no right gel number.\n", contig);
	    status++;
	}

	if (c.left && lnbr[c.left] != 0) {
	    vmessage("Contig %d: left gel (%d) has leftward neighbour.\n",
		   contig, c.left);
	    status++;
	}
	if (c.right && rnbr[c.right] != 0) {
	    vmessage("Contig %d: right gel (%d) has rightward neighbour.\n",
		   contig, c.right);
	    status++;
	}

	/*
	 * Chain right checking we end up at the correct place. Also check
	 * for loops and holes.
	 */
	gel = c.left;
	last = 0;
	len = 2;
	while (!bad_gel(io, gel)) {
	    /* check loop */
	    if (gel_used[gel] > 0) {
		vmessage("Contig %d: reading %d used twice (loop) to right.\n",
		       contig, gel);
		loop = 1;
		status++;
		break;
	    } else {
		gel_used[gel]++;
	    }

	    /* check for a 'hole' in contig */
	    if (relpg[gel] >= len) {
		vmessage("Contig %d: not contiguous between gel %d and %d.\n",
		       contig, gel, last);
		if (relpg[gel] == len)
		    (*minor)++;
		else
		    status++;
	    }

	    /* update max length so far */
	    if (relpg[gel] + abs(lngthg[gel]) > len)
		len = relpg[gel] + abs(lngthg[gel]);

	    last = gel;
	    gel = rnbr[gel];
	}

	if (loop) {
	    vmessage("Contig %d: aborting further right-checks.\n", contig);
	} else {
	    /* check length */
	    if (len-1 != c.length) {
		vmessage("Contig %d: has length %d, but chaining right gives "
		       "length %d.\n", contig, c.length, len-1);
		status++;
	    }
	    
	    /* check right end is correct */
	    if (last != c.right) {
		vmessage("Contig %d: right gel (%d) is not found by chaining "
		       "right from left gel.\n", contig, c.right);
		status++;
	    }

	    /* check termination conditions */
	    if (gel != 0) {
		/* invalid gel no. */
		vmessage("Contig %d: invalid gel no %d. (rnbr[%d] = %d).\n",
		       contig, gel, last, rnbr[last]);
		status++;
	    }
	}

	/*
	 * Chain left checking that we end up in the correct place. We do
	 * not bother to perform the extra checks during right chaining as
	 * these are unlikely to change. Chaining left may seem redundant, but
	 * in the case of incorrect hand holding this may yield more
	 * information than right chaining alone.
	 */
	gel = c.right;
	last = gel;
	loop = 0;
	while (!bad_gel(io, gel)) {
	    /* loop check */
	    if (gel_used[gel] > 1) {
		vmessage("Contig %d: reading %d used twice (loop) to left.\n",
		       contig, gel);
		status++;
		loop = 1;
		break;
	    } else {
		gel_used[gel]++;
	    }

	    last = gel;
	    gel = lnbr[gel];
	}

	if (loop) {
	    vmessage("Contig %d: aborting further left-checks.\n", contig);
	} else {
	    /* check left end is correct */
	    if (last != c.left) {
		vmessage("Contig %d: left gel (%d) is not found by chaining "
		       "left from right gel.\n", contig, c.left);
		status++;
	    }

	    /* check termination conditions */
	    if (gel != 0) {
		vmessage("Contig %d: invalid gel no %d. (lnbr[%d] = %d.\n",
		       contig, gel, last, lnbr[last]);
		status++;
	    }
	}

	/* check annotations list */
	anno = c.annotations;
	last_pos = 1;
	last_anno = 0;
	while (anno && GT_Read(io, arr(GCardinal, io->annotations, anno-1),
			       &a, sizeof(a), GT_Annotations) == 0) {
	    if (tag_used[anno]) {
		vmessage("Contig %d: annotation %d used more than once "
			 "(loop?).\n", contig, anno);
		status++;
		break;
	    }

	    tag_used[anno]++;

	    if (a.position < 1 || a.position + a.length > c.length + 1) {
		vmessage("Annotation %d: Pos (%d-%d), outside of contig %d.\n",
		       anno, a.position, a.position + a.length, contig);
		status++;
	    }

	    if (a.position < last_pos) {
		vmessage("Annotation %d: Pos (%d), leftwards of previous tag "
		       "%d (Pos %d).\n",
		       anno, a.position, last_anno, last_pos);
		status++;
	    }

	    last_anno = anno;
	    last_pos = a.position;
	    anno = a.next;
	    if (anno <= 0 || anno > Nannotations(io))
		anno = 0;
	}

	/* Fill in note_used[] array */
	note = c.notes;
	if (note) {
	    note_read(io, note, n);
	    /* Check prev/prev_type of 1st item */
	    if (n.prev_type != GT_Contigs ||
		n.prev != contig) {
		vmessage("Contig %d: note %d links back to prev=%d prev_type=%d\n",
			 contig, note, n.prev, n.prev_type);
		status++;
	    }
	}
	while (note) {
	    if (note_used[note]) {
		vmessage("Contig %d: note %d used more than once (loop?).\n",
			 contig, note);
		status++;
		break;
	    }
	    
	    note_used[note]++;
	    if (note = n.next)
		note_read(io, note, n);
	}

    }

    return status;
}

int check_readings(GapIO *io, f_int *relpg, f_int *lngthg, f_int *lnbr,
		   f_int *rnbr, int *gel_used, int *tag_used,
		   int *note_used, int *minor) {
    int gel, right, left;
    GReadings r;
    GAnnotations a;
    int anno;
    int status = 0;
    int last_pos, last_anno;
    int note;
    GNotes n;
    char *seq;

    for (gel = 1; gel <= NumReadings(io); gel++) {
	gel_read(io, gel, r);
#if GAP_CACHE!=0
	{
	    GReadings r_disk;
	    GT_Read(io, arr(GCardinal, io->readings, gel-1),
		    &r_disk, sizeof(r_disk), GT_Readings);

	    if (memcmp(&r, &r_disk, sizeof(r)) != 0) {
		vmessage("Gel %d: Cached copy is not same as disk copy\n",
			 gel);
		status++;
	    }
	}
#endif
	/* check memory vs disk */
	if (lnbr[gel] != r.left) {
	    vmessage("Gel %d: Memory left = %d, disk left = %d.\n",
		   gel, lnbr[gel], r.left);
	    status++;
	}
	if (rnbr[gel] != r.right) {
	    vmessage("Gel %d: Memory right = %d, disk right = %d.\n",
		   gel, rnbr[gel], r.right);
	    status++;
	}
	if (relpg[gel] != r.position) {
	    vmessage("Gel %d: Memory position = %d, disk position = %d.\n",
		   gel, relpg[gel], r.position);
	    status++;
	}
	if (lngthg[gel] != (r.sense
			    ? -r.sequence_length
			    : +r.sequence_length)) {
	    vmessage("Gel %d: Memory length = %d, disk sense;length = %d;%d.\n",
		   gel, lngthg[gel], r.sense, r.sequence_length);
	    status++;
	}

	/* check validity of left/right */
	if (bad_gel(io, left = lnbr[gel]) == 2) {
	    vmessage("Gel %d: left neighbour (%d) is invalid.\n",
		   gel, left);
	    left = -1;
	    status++;
	}

	if (bad_gel(io, right = rnbr[gel]) == 2) {
	    vmessage("Gel %d: right neighbour (%d) is invalid.\n",
		   gel, right);
	    right = -1;
	    status++;
	}

	/* how many times used */
	switch(gel_used[gel]) {
	case 0:
	    vmessage("Gel %d: never used.\n", gel);
	    (*minor)++;
	    break;
	case 1:
	    vmessage("Gel %d: used only in one direction.\n", gel);
	    status++;
	    break;
	case 2:
	    break;
	default:
	    vmessage("Gel %d: used %d times.\n", gel, gel_used[gel]-1);
	    status++;
	}
	
	/* hand holding - need only check either left or right (not both) */
	if (right > 0 && lnbr[right] != gel) {
	    status++;
	    vmessage("Gel %d: hand holding problem.\n", gel);
	    vmessage("    gel:%04d left:%04d right:%04d\n",
		   gel, left, right);
	    vmessage("    gel:%04d left:%04d right:%04d\n",
		   right, lnbr[right], rnbr[right]);
	}

	/* relative positioning */
	if (left > 0 && relpg[gel] < relpg[left]) {
	    vmessage("Gel %d: positioned leftwards of its left neighbour, %d\n",
		   gel, left);
	    status++;
	}

	/* valid lengths */
	if (lngthg[gel] == 0) {
	    vmessage("Gel %d: has zero length.\n", gel);
	    status++;
	}
	if (r.sequence_length +1 != r.end - r.start) {
	    vmessage("Gel %d: start and end of clips do not correspond with "
		   "used sequence length.\n", gel);
	    status++;
	}
	if (r.sequence_length < 0) {
	    vmessage("Gel %d: sequence_length is less than zero.\n", gel);
	    status++;
	}

	/* valid fields for strand, primer, & sense */
	if (r.strand < 0 || r.strand > 1) {
	    vmessage("Gel %d: invalid value for strand field, %d\n",
		   gel, r.strand);
	    status++;
	}
	if (r.primer < 0 || r.primer > 4) {
	    vmessage("Gel %d: invalid value for primer field, %d\n",
		   gel, r.primer);
	    status++;
	}
	if (r.sense < 0 || r.sense > 1) {
	    vmessage("Gel %d: invalid value for sense field, %d\n",
		   gel, r.sense);
	    status++;
	}

	/* fill in tag_used[] array */
	anno = r.annotations;
	last_pos = 1;
	last_anno = 0;
	while (anno && GT_Read(io, arr(GCardinal, io->annotations, anno-1),
			       &a, sizeof(a), GT_Annotations) == 0) {
	    if (tag_used[anno]) {
		vmessage("Gel %d: annotation %d used more than once (loop?).\n"
			 , gel, anno);
		status++;
		break;
	    }

	    tag_used[anno]++;

	    if (a.position < 1 ||
		a.position + a.length > r.length+1) {
		vmessage("Annotation %d: Pos (%d-%d), outside of gel %d.\n",
		       anno, a.position, a.position + a.length, gel);
		(*minor)++;
	    }

	    if (a.position < last_pos) {
		vmessage("Annotation %d: Pos (%d), leftwards of previous tag "
		       "%d (Pos %d).\n",
		       anno, a.position, last_anno, last_pos);
		(*minor)++;
	    }

	    last_anno = anno;
	    last_pos = a.position;
	    anno = a.next;
	    if (anno <= 0 || anno > Nannotations(io))
		anno = 0;
	}

	/* Fill in note_used[] array */
	note = r.notes;
	if (note) {
	    note_read(io, note, n);
	    /* Check prev/prev_type of 1st item */
	    if (n.prev_type != GT_Readings ||
		n.prev != gel) {
		vmessage("Gel %d: note %d links back to prev=%d prev_type=%d\n"
			 , gel, note, n.prev, n.prev_type);
		status++;
	    }
	}
	while (note) {
	    if (note_used[note]) {
		vmessage("Gel %d: note %d used more than once (loop?).\n",
			 gel, note);
		status++;
		break;
	    }

	    note_used[note]++;
	    if (note = n.next)
		note_read(io, note, n);
	}


	/* Check sequence for valid chars */
	seq = TextAllocRead(io, r.sequence);
	if (!seq) {
	    vmessage("Gel %d: sequence not readable\n", gel);
	    status++;
	} else {
	    int i;
	    for (i = 0; i < r.length; i++) {
		if (!isprint(seq[i])) {
		    vmessage("Gel %d: contains non-printable characters\n",
			     gel);
		    status++;
		    break;
		}
	    }
	    xfree(seq);
	}
    }

    return status;
}

int check_annotations(GapIO *io, int *tag_used, int *minor) {
    GAnnotations a;
    int t;
    int status = 0;
    int *tag_free;

    /* allocate temporary memory */
    if (NULL == (tag_free=(int *)xmalloc((Nannotations(io)+1) * sizeof(int)))){
	vmessage("Out of memory.\n");
	verror(ERR_WARN, "check_database", "Out of memory");
	return 1;
    }

    memset(tag_free, 0, (Nannotations(io)+1) * sizeof(*tag_free));

    /* scan free annotations list */
    t = io->db.free_annotations;
    while (t) {
	if (tag_free[t] != 0) {
	    vmessage("Annotation %d: loop detected in free list.\n", t);
	    status++;
	    break;
	}
	tag_free[t]++;
	if (GT_Read(io, arr(GCardinal, io->annotations, t-1),
		    &a, sizeof(a), GT_Annotations)) {
	    GAP_ERROR("reading annotation");
	    status++;
	    break;
	}
	t = a.next;
    }

    /* check on used annotations ? */
    for (t = 1; t <= Nannotations(io); t++) {
	(void)GT_Read(io, arr(GCardinal, io->annotations, t-1),
		      &a, sizeof(a), GT_Annotations);

	/* how many times used */
	if (tag_used[t] == 0 && tag_free[t] == 0) {
	    vmessage("Annotation %d: Neither used or free.\n", t);
	    (*minor)++;
	} else if (tag_used[t] > 1) {
	    vmessage("Annotation %d: used %d times.\n", t, tag_used[t]);
	    status++;
	}
	if (tag_used[t] != 0 && tag_free[t] != 0) {
	    vmessage("Annotation %d: used %d time%s, yet is on the free list.\n",
		   t, tag_used[t], tag_used[t] == 1 ? "" : "s");
	    status++;
	}

	/* check valid fields */
	if (a.length < 0) {
	    vmessage("Annotation %d: negative length (%d).\n", t, a.length);
	    status++;
	}
	if (a.strand < 0 || a.strand > 2) {
	    vmessage("Annotation %d: invalid value for strand field, %d.\n",
		   t, a.strand);
	    status++;
	}
    }

    xfree(tag_free);

    return status;
}

int check_notes(GapIO *io, int *note_used, int *minor) {
    GNotes n;
    int t;
    int status = 0;
    int *note_free;
    int *next, *prev;

    /* allocate temporary memory */
    if (NULL == (note_free=(int *)xmalloc((Nnotes(io)+1) * sizeof(int))) ||
	NULL == (next = (int *)xmalloc((Nnotes(io)+1) * sizeof(int))) ||
	NULL == (prev = (int *)xmalloc((Nnotes(io)+1) * sizeof(int)))) {
	vmessage("Out of memory.\n");
	verror(ERR_WARN, "check_database", "Out of memory");
	return 1;
    }

    memset(note_free, 0, (Nnotes(io)+1) * sizeof(*note_free));
    memset(next, 0, (Nnotes(io)+1) * sizeof(*next));
    memset(prev, 0, (Nnotes(io)+1) * sizeof(*prev));

    /* scan free notes list */
    t = io->db.free_notes;
    while (t) {
	if (note_free[t] != 0) {
	    vmessage("Note %d: loop detected in free list.\n", t);
	    status++;
	    break;
	}
	note_free[t]++;
	if (note_read(io, t, n)) {
	    GAP_ERROR("reading note");
	    status++;
	    break;
	}
	t = n.next;
    }

    /* check on used notes ? */
    for (t = 1; t <= Nnotes(io); t++) {
	note_read(io, t, n);

	/* initialise hand holding arrays */
	next[t] = n.next;
	prev[t] = n.prev;

	/* how many times used */
	if (note_used[t] == 0 && note_free[t] == 0) {
	    vmessage("Note %d: Neither used or free.\n", t);
	    (*minor)++;
	} else if (note_used[t] > 1) {
	    vmessage("Note %d: used %d times.\n", t, note_used[t]);
	    status++;
	}
	if (note_used[t] != 0 && note_free[t] != 0) {
	    vmessage("Note %d: used %d time%s, yet is on the free list.\n",
		   t, note_used[t], note_used[t] == 1 ? "" : "s");
	    status++;
	}

	/* check valid fields */
    }

    /* Check hand holding between notes */
    for (t = 1; t <= Nnotes(io); t++) {
	if (next[t] && prev[next[t]] != t) {
	    vmessage("Note %d: hand holding problem.\n", t);
	    vmessage("    note %04d left:%04d right:%04d\n",
		     t, next[t], prev[t]);
	    vmessage("    note %04d left:%04d right:%04d\n",
		     next[t], next[next[t]], prev[next[t]]);
	    status++;
	}
    }

    xfree(note_free);
    xfree(next);
    xfree(prev);

    return status;
}

int check_templates(GapIO *io, int *minor) {
    int i;
    GTemplates t;
    int status = 0;

    for (i = 1; i <= Ntemplates(io); i++) {
	(void)GT_Read(io, arr(GCardinal, io->templates, i-1),
		      &t, sizeof(t), GT_Templates);

	/* check insert lengths */
	if (t.insert_length_min > t.insert_length_max) {
	    vmessage("Template %d: minimum insert length (%d) greater than "
		   "the maximum (%d).\n",
		   i, t.insert_length_min, t.insert_length_max);
	    status++;
	}

	/* check valid fields */
	if (t.vector > Nvectors(io) || t.vector < 0) {
	    vmessage("Template %d: invalid vector number %d.\n",
		   i, t.vector);
	    status++;
	}

	if (t.clone > Nclones(io) || t.clone < 1) {
	    vmessage("Template %d: invalid clone number %d.\n",
		   i, t.clone);
	    status++;
	}

	if (t.strands < 1 || t.strands > 2) {
	    vmessage("Template %d: invalid strand %d.\n",
		   i, t.strands);
	    status++;
	}
    }

    return status;
}

int check_vectors(GapIO *io, int *minor) {
    int i;
    GVectors v;
    int status = 0;

    for (i = 1; i <= Nvectors(io); i++) {
	(void)GT_Read(io, arr(GCardinal, io->vectors, i-1),
		      &v, sizeof(v), GT_Vectors);

	/* valid level */
	if (v.level < 0) {
	    vmessage("Vector %d: Invalid level %d.\n", i, v.level);
	    status++;
	}
	if (v.level > MAX_LEVEL) {
	    vmessage("Vector %d: Absurdly large level %d.\n", i, v.level);
	    status++;
	}
    }

    return status;
}

int check_clones(GapIO *io, int *minor) {
    int i;
    GClones c;
    int status = 0;

    for (i = 1; i <= Nclones(io); i++) {
	(void)GT_Read(io, arr(GCardinal, io->clones, i-1),
		      &c, sizeof(c), GT_Clones);

	/* valid vector */
	if (c.vector < 1 || c.vector > Nvectors(io)) {
	    vmessage("Clone %d: invalid vector number %d.\n", i, c.vector);
	    status++;
	}
    }

    return status;
}

int check_order(GapIO *io, int *minor) {
    int i, nc = NumContigs(io);
    int *used;

    if (NULL == (used = (int *)xcalloc(nc+1, sizeof(int)))) {
	vmessage("Out of memory.\n");
	verror(ERR_WARN, "check_database", "Out of memory");
	return 1;
    }

    for (i=0; i<nc; i++) {
	int x = arr(GCardinal, io->contig_order, i);
	if (x >= 0 && x <= nc) {
	    used[x-1]++;
	}
    }

    /* Each contig should be used once in the order array */
    for (i=0; i<nc; i++)
	if (used[i] != 1) {
	    vmessage("Database: Contig order is inconsistent.\n");
	    xfree(used);
	    return 1;
	}

    xfree(used);
    return 0;
}

/*---------------------------------------------------------------------------*/
/* The common entry point. Returns 0 for OK, 1 for minor, 2 for fatal */
int db_check_common(GapIO *io, int dbsize, int ngels, int nconts,
		    int *relpg, int *lngthg, int *lnbr, int *rnbr) {
    int status = 0;
    int minor = 0;
    int *gel_used;
    int *tag_used;
    int *note_used;

    if (Nreadings(io) == 0 && Ncontigs(io) == 0)
	return 0;

    if (NULL == (gel_used = (int *)xmalloc((Nreadings(io)   +1)*sizeof(int)))||
	NULL == (tag_used = (int *)xmalloc((Nannotations(io)+1)*sizeof(int)))||
	NULL == (note_used = (int*)xmalloc((Nnotes(io)      +1)*sizeof(int)))){
	if (gel_used) xfree(gel_used);
	if (tag_used) xfree(tag_used);
	return 2;
    }

    memset(gel_used,  0, (Nreadings(io)   +1) * sizeof(*gel_used));
    memset(tag_used,  0, (Nannotations(io)+1) * sizeof(*tag_used));
    memset(note_used, 0, (Nnotes(io)      +1) * sizeof(*note_used));

    set_gap_fatal_errors(0);
    log_vmessage(1);

    status += check_database(io, dbsize, ngels, nconts, note_used, &minor);
    status += check_order(io, &minor);
    status += check_contigs(io, relpg, lngthg, lnbr, rnbr,
			    gel_used, tag_used, note_used, &minor);
    status += check_readings(io, relpg, lngthg, lnbr, rnbr,
			     gel_used, tag_used, note_used, &minor);
    status += check_annotations(io, tag_used, &minor);
    status += check_templates(io, &minor);
    status += check_vectors(io, &minor);
    status += check_clones(io, &minor);
    status += check_notes(io, note_used, &minor);

    log_vmessage(0);
    set_gap_fatal_errors(1);

    if (status) {
	vmessage("Database is not consistent. %d problems detected.\n",
	       status + minor);
	verror(ERR_WARN, "check_database",
	       "Database is not consistent. %d problems detected.\n",
	       status + minor);
	status = ignore_checkdb ? 1 : 2;
    } else if (minor) {
	vmessage("Database is not consistent. %d minor problems detected.\n",
	       minor);
	verror(ERR_WARN, "check_database",
	       "Database is not consistent. %d minor problems detected.\n",
	       minor);
	status = 1;
    } else {
	vmessage("Database is logically consistent\n");
	status = 0;
    }

    xfree(gel_used);
    xfree(tag_used);
    xfree(note_used);

    return status;
}


/*---------------------------------------------------------------------------*/
/* The fortran entry point. */

f_proc_ret dbchek_(f_int *handle, f_int *relpg, f_int *lngthg, f_int *lnbr,
		   f_int *rnbr, f_int *idm, f_int *idbsiz, f_int *ngels,
		   f_int *nconts, f_int *ierr) {
    GapIO *io;

    if (NULL == (io = io_handle(handle))) {
	verror(ERR_FATAL, "check_database", "invalid io handle");
	printf("Invalid file handle '%d'\n", handle ? *handle : -1);
	f_proc_return();
    }

    *ierr = db_check_common(io, *idbsiz, *ngels, *nconts,
			    relpg-1, lngthg-1, lnbr-1, rnbr-1);

    f_proc_return();
}


/*---------------------------------------------------------------------------*/
/* The C entry point. */

int db_check(GapIO *io) {
    return db_check_common(io, io_dbsize(io), NumReadings(io), NumContigs(io),
			   io->relpos, io->length, io->lnbr, io->rnbr);
}
