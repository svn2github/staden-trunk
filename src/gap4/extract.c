/*
 * File: extract.c:
 * Version: 2.0 (1.0 == FORTRAN version)
 *
 * Author: James Bonfield
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: "Extract gel readings" code.
 *
 * Created: 28 March 1994
 * Updated: Replaced strerror calls with perror as strerror. Despite being a
 *              ANSI defined function, SunOS does not have it.
 * 06/05/94 Added consensus extraction too.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

#include "IO.h"
#include "fort.h"
#include "FtoC.h"
#include "edUtils.h"
#include <io_lib/expFileIO.h>
#include "qual.h"
#include "xalloc.h"
#include "io_utils.h"
#include "misc.h"
#include "extract.h"
#include "text_output.h"
#include "notes.h"

/* 6/1/99 johnt - S_ISDIR not defined with Visual C++ */
#ifndef S_ISDIR
#  define S_ISDIR(m)	(((m)&S_IFMT) == S_IFDIR)
#endif

/*
 * Defines:
 *
 * EXTRACT_ORIG_ORDER
 * 	The experiment files are created with the sequence in the original
 *	sense. That is - we try and recreate the original experiment file
 *	as closely as possible.
 *
 * NB: Without this defined the output_vector function will need some
 * modifications.
 */
#define EXTRACT_ORIG_ORDER

/*
 * LibC
 */
/*
 * Evaluates a pathname.
 *
 * Replaces ~ with users home directory
 *          ~usr with user 'usr's home directory
 *          $var with the value of environment variable 'var'
 */
/*
FILE *open_path() {
}
*/

/*
 * Output TG lines.
 * Common code used by both the gel and consensus experiment file routines
 *
 * Offset is a value to add to all annotations. Used by the consensus output
 * to shift gel annotations onto the consensus coordinates.
 *
 * If orig is set, then we output the tags relative to their original positions
 * in the sequences. Otherwise we swap around so that they are in their
 * positions as seen in the contig editor. g_sense and g_len need only be sent
 * over when orig is 0.
 * 
 * clip_left and clip_right, when non zero values, are used to clip annotations
 * outside of a given range. This currently only works (as it's only needed)
 * when not shuffling the tags to their original senses.
 *
 * If cc_line != NULL, a CC line is outputed before the first annotation.
 *
 * We don't output SVEC and CVEC tags as these are included in the SL, SR
 * and CS lines.
 *
 * Returns 0 for success, -1 for failure.
 */
int output_annotations(GapIO *io, Exp_info *e, int anno, int offset,
		       int orig, int g_sense, int g_len, int consensus,
		       int clip_left, int clip_right, char *cc_line,
		       int *pads, int npads) {
    char *buf;
    int pos, err = 0;
    GAnnotations a;
    char type[5];
    int svec = str2type("SVEC");
    int cvec = str2type("CVEC");
    char *comment;

    for (; anno; anno = a.next) {
	GT_Read(io, arr(GCardinal, io->annotations, anno-1),
		&a, sizeof(a), GT_Annotations);

	if (a.type == svec || a.type == cvec)
	    continue;

	if (orig) {
	    pos = a.position;
	} else {
	    if (g_sense == GAP_SENSE_REVERSE) {
		pos = g_len - (a.position + a.length - 1) + 1;
		if (a.strand != 2)
		    a.strand = !a.strand;
	    } else {
		pos = a.position;
	    }
	}

	if (clip_left || clip_right) {
	    if (pos < clip_left + 1) {
		a.length -= clip_left + 1 - pos;
		pos = clip_left + 1;
	    }
	    
	    if (pos + a.length > clip_right) {
		a.length = clip_right - pos;
	    }

	    if (a.length <= 0) {
		continue;
	    }
	}

	/* read comment */
	if (a.annotation)
	    comment = TextAllocRead(io, a.annotation);
	else
	    comment = NULL;
	
	/* read type */
	type[0] = (a.type >> 030) & 0xff;
	type[1] = (a.type >> 020) & 0xff;
	type[2] = (a.type >> 010) & 0xff;
	type[3] = (a.type >> 000) & 0xff;
	type[4] = '\0';

	/*
	 * If we've provided a pads array, strip out the pads from the
	 * tag.
	 */
	if (pads) {
	    int p, p2;

	    p = pos + offset - 1;
	    p2 = p + a.length - 1;

	    if (p < 0) p = 0;
	    if (p >= npads) p = npads-1;
	    if (p2 >= npads) p2 = npads-1;
	    pos -= pads[p];
	    if (p2 >= 0) {
		if (p > 0) {
		    a.length -= pads[p2] - pads[p-1];
		    if (pads[p-1] != pads[p])
			pos++;
		} else {
		    a.length -= pads[p2];
		    if (0 != pads[p])
			pos++;
		}
	    }
	}

	if (a.length - 1 >= 0 && pos + offset > 0) {
	    /* create tag string */
	    buf = (char *)xmalloc(100 + (comment ? strlen(comment) : 0));
	    if (!buf) {
		if (comment)
		    xfree(comment);
		return -1;
	    }

	    values2tag(buf, type,
		       pos + offset, pos + offset + a.length - 1,
		       a.strand, comment);
	    
	    if (cc_line) {
		err |= exp_put_str(e, EFLT_CC, cc_line, strlen(cc_line));
		cc_line = NULL;
	    }

	    err |= exp_put_str(e, consensus ? EFLT_TC : EFLT_TG,
			       buf, strlen(buf));

	    xfree(buf);
	}

	if (comment) {
	    xfree(comment);
	    comment = NULL;
	}
    }

    return err;
}


/*
 * Output NT lines (notes).
 * Common code used by both the gel and consensus experiment file routines
 *
 * Returns 0 for success, -1 for failure.
 */
int output_notes(GapIO *io, Exp_info *e, int note,
		 int source_type, int source_num) {
    GNotes n;
    char *str;

    for (; note; note = n.next) {
	note_read(io, note, n);
	str = note2str(io, n, source_type, source_num);
	exp_put_str(e, EFLT_NT, str, strlen(str));
    }

    return 0;
}

/*
 * Outputs the SL, SR and CS vector information lines.
 *
 * Note: If there is no SR tag, no SR line is outputted, even if we have
 * outputted SL. This is not the same as vepe.
 */
int output_vector(GapIO *io, Exp_info *e, int gel, int glen) {
    char *types[] = {"SVEC", "CVEC"};
    int svec = str2type("SVEC");
    GAnnotations *a;

    /* SL and SR */
    a = vtagget(io, gel, sizeof(types)/sizeof(*types), types);

    while (a && a != (GAnnotations *)-1) {
	if (a->type == svec) {
	    if (a->position == 1) { /* template start */
		exp_put_int(e, EFLT_SL, &a->length);
	    } else if (a->position + a->length == glen + 1) {
		exp_put_int(e, EFLT_SR, &a->position);
	    }
	} else /* CVEC */ {
	    int st, en;

	    st = a->position;
	    en = a->position + a->length - 1;
	    exp_put_rng(e, EFLT_CS, &st, &en);
	}

	a = vtagget(io, 0, sizeof(types)/sizeof(*types), types);
    }

    return 0;
}

/*
 * Temporary structure for sorting readings into contig and position.
 */
typedef struct {
    GapIO *io;
    char *name;
    int rnum;
    int contig;
    int position;
    int anchor;
} sort_struct;

/*
 * Creates an experiment file from a given reading number
 *
 * Notes:
 * o      We do not know, from within gap, whether the left and right cutoffs
 *        are due to quality or vector. We therefore assume that left cutoff
 *        is from vector, and right is quality (which is normally true).
 *
 * o      Only outputs AV, ON and PC lines when format == 1.
 *
 * Formats:
 * 0 Just the sequence bits without any assembly details
 * 1 The old 'preassembled data' format (AV, ON, PC lines)
 * 2 Disassemble readings format, anchored to next left (AP)
 * 3 Disassemble readings format, anchored to first reading (AP)
 *
 * returns 0 for success, -1 for failure.
 */
static int create_exp_file(GapIO *io, char *name, int format,
			   sort_struct *readinfo) {
    char *buf = NULL;
    GReadings r;
    GTemplates t;
    GVectors v;
    GClones c;
    Exp_info *e;
    int err = 0, tmp;
    mFILE *fp;
    int1 *conf = NULL;
    int2 *opos = NULL;
    char *seq = NULL;
    int retcode = -1;
    int buf_size;
    int gel;
    
    gel = readinfo->rnum;
    gel_read(io, gel, r);

    /* The largest thing in buf is opos array */
    buf_size = r.length * 11 + 1;
    if (NULL == (buf = (char *)xmalloc(buf_size))) 
	goto error;

    if (NULL == (fp = mfopen(name, "w"))) {
	perror(name);
	goto error;
    }

    /* allocate experiment file info */
    e = exp_create_info();
    e->fp = fp;

    /*
     * Fill experiment file with data
     *
     * exp_put_int(e, int id, int *val);              "id   %d"
     * exp_put_rng(e, int id, int *from, int *to);    "id   %d..%d";
     * exp_put_str(e, int id, char *string, int len); "id   %s";
     */
    /* reading */
    if (!r.sequence) {
	/*
	 * Corrupted sequence - skip it.
	 * This check is only here for using extract_readings for recovery
	 * from nasty problems.
	 */
	goto error;
    }

    strcpy(buf, io_rname(io, gel));
    err |= exp_put_str(e, EFLT_ID, buf, strlen(buf));
    err |= exp_put_str(e, EFLT_EN, buf, strlen(buf));

    err |= exp_put_str(e, EFLT_DR,
		       STRAND(r) == GAP_STRAND_FORWARD ? "+" : "-", 1);

    seq = TextAllocRead(io, r.sequence);
    if (r.confidence)
	conf = (int1 *)DataAllocRead(io, r.confidence, sizeof(*conf));

    if (r.orig_positions)
	opos = (int2 *)DataAllocRead(io, r.orig_positions, sizeof(*opos));

#ifdef EXTRACT_ORIG_ORDER
#   define ORDER 1
    /*
     * If complemented, the recomplement sequence back
     */
    if (GAP_SENSE_REVERSE == r.sense) {
	int length = r.length, start = r.start, end = r.end;

	io_complement_seq(&length, &start, &end, seq, conf, opos);
	r.length = length;
	r.start = start;
	r.end = end;
    }
#else
#   define ORDER 0
#endif

    err |= exp_put_int(e, EFLT_QL, &r.start);
    err |= exp_put_int(e, EFLT_QR, &r.end);
    
    if (r.trace_name) {
	TextRead(io, r.trace_name, buf, buf_size);
	err |= exp_put_str(e, EFLT_LN, buf, strlen(buf));
    }

    if (r.trace_type) {
	TextRead(io, r.trace_type, buf, buf_size);
	err |= exp_put_str(e, EFLT_LT, buf, strlen(buf));
    }

    /* reading.annotation+ */
    err |= output_annotations(io, e, r.annotations, 0, ORDER,
			      r.sense, r.length, 0, 0,0, NULL, NULL, 0);
    
    /* reading.notes */
    err |= output_notes(io, e, r.notes, GT_Readings, gel);

    /* reading.template */
    if (r.template != 0) {
	GT_Read(io, arr(GCardinal, io->templates, r.template-1),
		&t, sizeof(t), GT_Templates);
	err |= exp_put_int(e, EFLT_ST, &t.strands);
	err |= exp_put_rng(e, EFLT_SI, &t.insert_length_min,
			   &t.insert_length_max);

	TextRead(io, t.name, buf, buf_size);
	Fstr2Cstr(buf, DB_NAMELEN, buf, buf_size);
	err |= exp_put_str(e, EFLT_TN, buf, strlen(buf));

	/* reading.template.vector - sequence vector */
	if (t.vector) {
	    GT_Read(io, arr(GCardinal, io->vectors, t.vector-1),
		    &v, sizeof(v), GT_Vectors);
	    TextRead(io, v.name, buf, buf_size);
	    Fstr2Cstr(buf, DB_NAMELEN, buf, buf_size);
	    err |= exp_put_str(e, EFLT_SV, buf, strlen(buf));
	}

	/* reading.template.clone */
	if (t.clone) {
	    GT_Read(io, arr(GCardinal, io->clones, t.clone-1),
		    &c, sizeof(c), GT_Clones);
	    TextRead(io, c.name, buf, buf_size);
	    Fstr2Cstr(buf, DB_NAMELEN, buf, buf_size);
	    err |= exp_put_str(e, EFLT_CN, buf, buf_size);

	    /* reading.template.clone.vector - cloning vector */
	    if (c.vector) {
		GT_Read(io, arr(GCardinal, io->vectors, c.vector-1),
			&v, sizeof(v), GT_Vectors);
		TextRead(io, v.name, buf, buf_size);
		Fstr2Cstr(buf, DB_NAMELEN, buf, buf_size);
		err |= exp_put_str(e, EFLT_CV, buf, strlen(buf));
	    }
	}
    }

    /* Chemistry type */
    err |= exp_put_int(e, EFLT_CH, &r.chemistry);

    /* Primer type */
    tmp = PRIMER_TYPE(r);
    err |= exp_put_int(e, EFLT_PR, &tmp);

    /* SL, SR and CS lines */
    err |= output_vector(io, e, gel, r.length);

    if (format == 1 || format == 2 || format == 3) {
	/* Original positions */
	if (r.orig_positions && r.length <= 32767) {
	    opos2str(opos, r.length, buf);
	    err |= exp_put_str(e, EFLT_ON, buf, strlen(buf));
	}

	/* Confidence values */
	if (r.confidence) {
	    conf2str(conf, r.length, buf);
	    err |= exp_put_str(e, EFLT_AV, buf, strlen(buf));
	}

	if (format == 1) {
	    /* PC line for preassemble data */
	    err |= exp_put_int(e, EFLT_PC, &r.position);

	    /* Sense */
	    err |= exp_put_int(e, EFLT_SE, &r.sense);
	} else if (format == 2 || format == 3) {
	    int left = readinfo->anchor;

	    if (left) 
		sprintf(buf, "%s %c %d -1",
			get_read_name(io, left),
			"+-"[r.sense],
			r.position - io_relpos(io, left));
	    else
		sprintf(buf, "*new* %c", "+-"[r.sense]);

	    err |= exp_put_str(e, EFLT_AP, buf, strlen(buf));

	    if (r.right == 0) {
		int offset;
		GContigs c;
		/* Last reading, so output consensus tags too */

		c.annotations = 0;
		contig_read(io, rnumtocnum(io, gel), c);

		if (r.sense) {
		    offset = r.position - (r.length - r.end) - 2;
		} else {
		    offset = -(r.position - r.start - 1);
		}
		output_annotations(io, e, c.annotations,
				   offset,
				   0, /* orig */
				   r.sense, /* sense */
				   r.length, /* g_length */
				   1, /* consensus */
				   0, /* clip_left */
				   0, /* clip_right */
				   NULL /* cc_line */,
				   NULL /* pads */,
				   0 /* npads */);
	    }
	}
    }

    /*
     * Sequence with // terminator always goes last
     */
    err |= exp_put_str(e, EFLT_SQ, seq, r.length);

    /* deallocate experiment file info - will close fp for us */
    exp_destroy_info(e);

    retcode = err ? -1 : 0;

 error:
    if (buf)
	xfree(buf);
    if (seq)
	xfree(seq);
    if (conf)
	xfree(conf);
    if (opos)
	xfree(opos);

    return retcode;
}

static int sort_readings_callback(const void *v1, const void *v2) {
    const sort_struct *s1 = (const sort_struct *)v1;
    const sort_struct *s2 = (const sort_struct *)v2;

    if (s1->contig < s2->contig)
	return -1;
    if (s1->contig > s2->contig)
	return 1;
    if (s1->position < s2->position)
	return -1;
    if (s1->position > s2->position)
	return 1;
    return 0;
}

/*
 * Sorts reading array by contig and position in contig.
 * Returns 0 for success,
 *        -1 for failure.
 */
static sort_struct *sort_readings(GapIO *io,
				  int num_readings,
				  char **reading_array,
				  int leftmost)
{
    sort_struct *sortme;
    int i;
    int anchor;
    int curr_contig = 0;

    sortme = (sort_struct *)xmalloc(num_readings * sizeof(*sortme));
    if (NULL == sortme)
	return NULL;

    for (i = 0; i < num_readings; i++) {
	int gel = get_gel_num(io, reading_array[i], GGN_ID);
	if (gel) {
	    sortme[i].io = io;
	    sortme[i].rnum = gel;
	    sortme[i].position = io_relpos(io, gel);
	    sortme[i].contig = rnumtocnum(io, gel);
	    sortme[i].name = reading_array[i];
	} else {
	    sortme[i].io = io;
	    sortme[i].rnum = 0;
	    sortme[i].position = 0;
	    sortme[i].contig = 0;
	    sortme[i].name = reading_array[i];
	}
    }

    qsort(sortme, num_readings, sizeof(*sortme), sort_readings_callback);

    for (anchor = i = 0; i < num_readings; i++) {
	if (curr_contig != sortme[i].contig) {
	    sortme[i].anchor = 0;
	    curr_contig = sortme[i].contig;
	    anchor = sortme[i].rnum;
	} else {
	    sortme[i].anchor = anchor;
	    if (!leftmost)
		anchor = sortme[i].rnum;
	}
    }

    return sortme;
}

/*
 * C interface
 *
 * Returns:
 * -1 failure
 *  0 success
 */
int extract_readings(GapIO *io,
		     int num_readings,
		     char **reading_array,
		     char *dir,
		     int format)
{
    char path[1024], *tmp;
    struct stat statbuf;
    int i;
    FILE *fp;
    sort_struct *sorted_reads;
    
    /*
     * Sort readings into contig position so that they can be correctly
     * pieced together again by the likes of directed assembly.
     */
    sorted_reads = sort_readings(io, num_readings, reading_array, format==3);

    /* check for directory existance, and create if needed */
    if (-1 == stat(dir, &statbuf)) {
	if (-1 == mkdir(dir, 0777)) {
	    verror(ERR_WARN, "extract_readings",
		   "Could not make directory '%s': %s",
		   dir, strerror(errno));
	    return -1;
	}
    } else {
	if (!S_ISDIR(statbuf.st_mode)) {
	    verror(ERR_WARN, "extract_readings",
		   "%s already exists and is not a directory", dir);
	    return -1;
	}
    }

    sprintf(path, "%s/fofn", dir);
    if (NULL == (fp = fopen(path, "w+"))) {
	verror(ERR_WARN, "extract_readings", "Couldn't create 'fofn' file");
	return -1;
    }

    for (i= 0; i < num_readings; i++) {
	if (NULL == (tmp = get_read_name(io, sorted_reads[i].rnum)))
	    continue;
	sprintf(path, "%s/%s", dir, tmp);
	fprintf(fp, "%s\n", tmp);
	vmessage("Creating experiment file %s\n", path);
	if (-1 == create_exp_file(io, path, format, &sorted_reads[i])){
	    verror(ERR_WARN, "extract_readings",
		   "Writing experiment file %s failed.",
		   tmp);
	}

	UpdateTextOutput();
    }

    xfree(sorted_reads);

    fclose(fp);
    return 0;
}
