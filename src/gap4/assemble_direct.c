/*
 * This file contains the code for the directed assembly mode.
 *
 * More importantly, it contains code for aligning and entering readings.
 * These functions are analogous to the preassemble code and to the stikit_()
 * function, but should be considered as superceding these. We should merge
 * in preassemble and stikit_() at some future stage.
 */

#include <string.h>
#include <tcl.h>
#include <stddef.h>

#include "io_utils.h"
#include "list_proc.h"
#include "assemble_direct.h"
#include "expFileIO.h"
#include "consen.h"
#include "gap_globals.h"
#include "seqInfo.h"
#include "list_proc.h"
#include "text_output.h"
#include "IO.h"
#include "IO2.h"
#include "misc.h"
#include "align.h"
#include "tagUtils.h"
#include "dna_utils.h"
#include "clones.h"
#include "notes.h"

typedef struct {
    int start;
    int end;
    int strand;
    char type[5];
    char *comment;
} anno_info;

typedef struct {
    int pos;
    int size;
} edit_info;

int remove_contig_holes_all(GapIO *io);

/*
 * Produces the consensus sequence for all contigs in the database.
 *
 * Returns consen_info struct for success, NULL for failure.
 */
consen_info *all_consensus(GapIO *io, float percd) {
    Contig_parms *clist;
    int task_mask = NORMALCONSENSUS | ADDTITLE;
    consen_info *ci = NULL;
    int i, *cends = NULL, *cnums = NULL;
    Hidden_params p;

    p.min = p.max = p.verbose = p.use_conf = p.qual_val = p.window_len =0;
    p.test_mode = 0;
    p.start = 0;
    p.lwin1 = 0;
    p.lcnt1 = 0;
    p.rwin1 = 0;
    p.rcnt1 = 0;
    p.do_it = 0;

    /* Alloc memory */
    if (NULL == (ci = (consen_info *)xcalloc(1, sizeof(consen_info))))
	return NULL;

    if (NULL == (ci->con_all = (char *)xmalloc(maxseq)))
	goto error;

    /* silly special case */
    if (NumContigs(io) == 0) {
	ci->con_len = 0;
	return ci;
    }

    ci->num_contigs = NumContigs(io);
    if (NULL == (ci->con_item = (char **)xmalloc(NumContigs(io) *
						 sizeof(char *))))
	goto error;

    if (NULL == (cends = (int *)xmalloc((NumContigs(io)+1) * sizeof(int))))
	goto error;

    if (NULL == (cnums = (int *)xmalloc((NumContigs(io)+1) * sizeof(int))))
	goto error;


    /* Generate the contig list */
    clist = get_contig_list(io_dbsize(io), io,
			    0 /* num_contigs (unused) */,
			    NULL /* contig_array => all */);


    /* Create con_all */
    if (0 != make_consensus(task_mask, io, ci->con_all, NULL,
			    clist, NumContigs(io),
			    &ci->con_len, max_gel_len(io), maxseq,
			    p, percd)) 
	goto error;

    /* Generate con_item */
    find_contig_ends(ci->con_all, ci->con_len, cends, cnums);
    ci->con_item[0] = &ci->con_all[20];
    for (i = 1; i < NumContigs(io); i++)
	ci->con_item[i] = &ci->con_all[cends[i]+20];

    xfree(cends);
    xfree(cnums);
    xfree(clist);

    return ci;

 error:
    if (ci) {
	if (ci->con_all)
	    xfree(ci->con_all);
	if (ci->con_item)
	    xfree(ci->con_item);
	xfree(ci);
    }

    if (cends)
	xfree(cends);
    if (cnums)
	xfree(cnums);

    return NULL;
}

void free_all_consensus(consen_info *ci) {
    if (ci) {
	if (ci->con_all)
	    xfree(ci->con_all);
	if (ci->con_item)
	    xfree(ci->con_item);
	xfree(ci);
    }

    return;
}

int realloc_consensus(consen_info *ci, int new_size) {
    int i;

    /* Turn con_item pointer array into offset array */
    for (i = 0; i < ci->num_contigs; i++)
	ci->con_item[i] = (char *)(ci->con_item[i] - ci->con_all);

    /* Realloc */
    maxseq = new_size * 1.5;
    if (NULL == (ci->con_all = xrealloc(ci->con_all, maxseq)))
	return -1;

    /* Turn con_item array back into pointers */
    for (i = 0; i < ci->num_contigs; i++)
	ci->con_item[i] = (ptrdiff_t)ci->con_item[i] + ci->con_all;

    return 0;
}

int recalc_consensus(GapIO *io, consen_info *ci, int contig, int pos,
		     int read_len, int old_clen, int new_clen) {
    char *seq;
    int i, diff, len;

    if (NumContigs(io) > ci->num_contigs) {
	ci->con_item = (char **)xrealloc(ci->con_item,
					 NumContigs(io) * sizeof(char *));
	if (NULL == ci->con_item)
	    return -1;
	for (i=ci->num_contigs; i<NumContigs(io); i++)
	    ci->con_item[i] = NULL;
	ci->num_contigs = NumContigs(io);
    }

    /* Create a new contig if necessary */
    if (ci->con_item[contig-1] == NULL) {
	char *cp = NULL;

	for (i = contig; i < NumContigs(io); i++) {
	    if (ci->con_item[i]) {
		cp = ci->con_item[i] - 20;
		break;
	    }
	}
	if (i == NumContigs(io))
	    cp = ci->con_all + ci->con_len;

	len = ci->con_all + ci->con_len - cp;
	if (cp + 20 + len - ci->con_all >= maxseq) {
	    ptrdiff_t cp_offset = cp - ci->con_all;

	    verror(ERR_WARN, "directed_assembly",
		   "consensus too large - increasing maxseq");
	    if (-1 == realloc_consensus(ci, cp + 20 + len - ci->con_all)) {
		verror(ERR_WARN, "directed_assembly",
		       "consensus too large");
		return -1;
	    }

	    cp = cp_offset + ci->con_all;
	}
	if (len > 0)
	    memmove(cp + 20, cp, len);
	add_contig_title(cp, "t.b", io_clnbr(io, contig));
	ci->con_item[contig-1] = cp + 20;
	ci->con_len += 20;
	for (i = contig; i < NumContigs(io); i++)
	    if (ci->con_item[i])
		ci->con_item[i] += 20;
    }

    diff = new_clen - old_clen;

    /* Shift data */
    if (pos <= 0)
	pos = 1;
    seq = ci->con_item[contig-1] + pos - 1;
    if (diff) {
        len = ci->con_all + ci->con_len - seq;
	if (seq + diff + len - ci->con_all >= maxseq) {
	    ptrdiff_t seq_offset = seq - ci->con_all;

	    verror(ERR_WARN, "directed_assembly",
		   "consensus too large - increasing maxseq");
	    if (-1 == realloc_consensus(ci, seq + diff + len - ci->con_all)) {
		verror(ERR_WARN, "directed_assembly",
		       "consensus too large");
		return -1;
	    }
	    seq = seq_offset + ci->con_all;
	}
	if (len > 0)
	    memmove(seq + diff, seq, len);
	ci->con_len += diff;

	for (i = contig; i < NumContigs(io); i++) {
	    if (ci->con_item[i])
		ci->con_item[i] += diff;
	}
    }

    /* Calculate new consensus */
    calc_consensus(contig, pos, pos + read_len, CON_SUM,
		   seq, NULL, NULL, NULL,
		   consensus_cutoff, quality_cutoff,
		   database_info, (void *)io);
    return 0;
}


/*
 * The AP line is of format:
 *
 * "reading_name orientation positition tolerance".
 * (eg "fred.seq - 250 30").
 *
 * The reading name is translated into a reading number. If the name is the
 * same as the reading name then reading number 0 is stored instead.
 * The position returned is an offset from the start of the contig, rather
 * than relative to the reading name (as is held in the AP line itself).
 * Contig refers to the contig number the reading is found in (zero if not).
 *
 * Returns 0 for a successful parse
 *        -1 for parse error
 *        -2 for unknown reading
 */
static int AP_parse(GapIO *io, char *readname, char *line, int *reading,
		    int *orient, int *position, int *contig, int *tolerance) {
    char s_orient[100], s_reading[100];
    int scanned;

    /* Get line words */
    scanned = sscanf(line, "%s %s %d %d", s_reading, s_orient,
		     position, tolerance);
    if (!(scanned >= 2 && strcmp(s_reading, "*new*") == 0) &&
	scanned != 4) {
	return -1;
    }

    /* Translate orient */
    *orient = *s_orient == '-' ? GAP_STRAND_REVERSE : GAP_STRAND_FORWARD;

    /* Translate reading name */
    if (strncmp(readname, s_reading, DB_NAMELEN) == 0 ||
	strcmp(s_reading, "*new*") == 0) {
	*reading = 0;
    } else if (-1 == (*reading = get_gel_num(io, s_reading, GGN_NAME))) {
	return -2;
    }

    /* Shift relative position to absolute position */
    if (*reading) {
	*position += io_relpos(io, *reading);
	*contig = rnumtocnum(io, *reading);
	if (-1 == *contig)
	    *contig = 0;
    } else {
	*contig = 0;
    }

    return 0;
}


/*
 * Creates an anno list in memory. In memory copies are faster to deal with
 * when padding - we only write them to disk once they've been fully edited.
 *
 * We also handle creation of additional tags here when exp_type == EFLT_TG
 * for marking the sequencing and cosmid vector.
 */
static anno_info *create_anno_list(SeqInfo *si, int exp_type, int *lenp,
				   int read_length) {
    int i, len, fudge;
    anno_info *ap, *a;

    *lenp = len = exp_Nentries(si->e, exp_type);
    fudge = (exp_type == EFLT_TG) ? 3 : 0;

    if (len+fudge == 0)
	return NULL;

    if (NULL == (ap = (anno_info *)xmalloc((len+fudge) * sizeof(anno_info)))) {
	return NULL;
    }

    for (i = 0, a = ap; i < len; i++) {
	char *comment, *tag;

	tag = arr(char *, si->e->entries[exp_type], i);
	if (NULL == (comment = (char *)xmalloc(strlen(tag))))
	    continue;

	if (-1 == tag2values(tag, a->type, &a->start, &a->end, &a->strand,
			     comment)) {
	    verror(ERR_WARN, "enter_reading", "Failed to parse annotation");
	    (*lenp)--;
	    continue;
	}

	if (*comment) {
	    a->comment = comment;
	} else {
	    a->comment = NULL;
	    xfree(comment);
	}
	
	a++;
    }

    if (exp_type == EFLT_TG) {
	/* Sequencing vector */
	if (exp_Nentries(si->e, EFLT_SL)) {
	    int end = atoi(exp_get_entry(si->e, EFLT_SL));
	    if (end > 0) {
	        strcpy(a->type, "SVEC");
	        a->start = 1;
	        a->end = end;
	        a->strand = 0;
	        a->comment = NULL;
	    
	        a++;
	        (*lenp)++;
	    }
	}

	/*
	 * Note vepe in the past (and maybe still) always created an SR line
	 * when an SL line was required. In these cases SR could be marking
	 * a zero length component of sequence - we check.
	 */
	if (exp_Nentries(si->e, EFLT_SR)) {
	    int start = atoi(exp_get_entry(si->e, EFLT_SR));

	    if (start <= si->length) {
		strcpy(a->type, "SVEC");
		a->start = start;
		a->end = read_length;
		a->strand = 0;
		a->comment = NULL;

		a++;
		(*lenp)++;
	    }
	}
	
	/* Cosmid vector */
	if (exp_Nentries(si->e, EFLT_CS)) {
	    exp_get_rng(si->e, EFLT_CS, &a->start, &a->end);
	    if (a->start < 1)
		a->start = 1;
	    if (a->end > read_length)
		a->end = read_length;
	    strcpy(a->type, "CVEC");
	    a->strand = 0;
	    a->comment = NULL;

	    a++;
	    (*lenp)++;
	}
    }

    if (*lenp == 0 && ap) {
	xfree(ap);
	return NULL;
    }

    return ap;
}


/*
 * Writes an annotation list to disk. If reading is a negative number then
 * it is assumed to be a contig number instead.
 *
 * Comp, start, end and length refer to the parameters for the reading that
 * this tag is from.
 *
 * Offset should contain the offset of this reading in the contig when writing
 * consensus tags.
 */
static void write_anno_list(GapIO *io, anno_info *ap, int len, int reading,
			    int offset, int comp, int start, int end,
			    int length) {
    int i, start2, end2;

    for (i = 0; i < len; i++) {
	/*
	 * Clip to used length of sequence when translating tags to the
	 * consensus.
	 */
	if (reading < 0) {
	    if (comp) {
		start2 = length - ap[i].end + 1;
		end2 = length - ap[i].start + 1;
	    } else {
		start2 = ap[i].start;
		end2 = ap[i].end;
	    }

	    if (start2 + offset < 0)
		start2 = 0 - offset;

	    /* Can't check yet as we don't know the contig length */
#if 0
	    if (end2 >= end)
		end2 = end-1;

	    if (end2 < start2)
		continue;
#endif

	} else {
	    start2 = ap[i].start;
	    end2 = ap[i].end;

	    if (start2 + offset < 1 || end2 + offset > length) {
		verror(ERR_WARN, "enter_reading",
		       "Tag on reading %s overlaps gel reading ends - not "
		       "entered", get_read_name(io, reading));
		continue;
	    } else if (end2 < start2) {
		verror(ERR_WARN, "enter_reading",
		       "Tag on reading %s has negative length - not entered",
		       get_read_name(io, reading));
		continue;
	    }
	}

	insert_NEW_tag(io, (tag_id)reading, start2 + offset, end2 - start2 + 1,
		       ap[i].type, ap[i].comment, ap[i].strand);
    }

    return;
}


/*
 * Frees a temporary anno list.
 */
static void free_anno_list(anno_info *ap, int len) {
    int i;

    for (i = 0; i < len; i++) {
	if (ap[i].comment)
	    free(ap[i].comment);
    }

    if (ap)
	xfree(ap);
}


/*
 * Modifies the tag start and end positions simulating an insertion at a
 * specific point. Tags will either stay the same (insert is rightwards of tag
 * end), grow (insert is within tag) or move right (insert is leftwards of tag
 * start).
 *
 * Uses the same algorithm as tag_shift_for_insert(), but works on memory
 * data rather than requiring database access. This is required as we must
 * perform edits on the new consensus tags before they are added to the
 * existing consensus tags.
 */
static void insert_to_anno_list(anno_info *ap, int len, int pos, int size) {
    int i;

    if (!ap)
	return;

    for (i = 0; i < len; i++) {
	anno_info *a = &ap[i];
	if (a->start >= pos) {
	    a->start += size;
	    a->end += size;
	} else if (a->end >= pos) {
	    a->end += size;
	}
    }
}


/*
 * Using an alignment buffer this performs edits on a sequence (and connected
 * information) and the consensus.
 */
static void edit_sequence(GapIO *io, align_info *ai,
			  int *alength, int *length, int *start, int *end,
			  char **seq_p, int1 **conf_p, int2 **opos_p,
			  int contig, int comp,
			  anno_info *anno_r, int anno_rl,
			  anno_info *anno_c, int anno_cl) {
    int i, j, k, M, N, op = 0, rem, div;
    int i2 = 0, j2 = 0;
    align_int *res = ai->res;
#define BLOCK 20
    char pads[BLOCK+1] = "********************";
    edit_info *ei, *eip;
    char *seq = *seq_p;
    int1 *conf = *conf_p;
    int2 *opos = *opos_p;

    i = ai->start1;
    j = ai->start2;
    M = ai->len1 + ai->start1;
    N = ai->len2 + ai->start2;

    if (NULL == (ei = eip = (edit_info *)xmalloc((M-i+1) * sizeof(edit_info))))
	return;

    while (i < M && j < N) {
	if (op == 0 && *res == 0) {
	    /* match */
	    op = *res++;
	    i++;
	    j++;
	} else {
	    if (op == 0)
		op = *res++;

	    if (op > 0) {
		/* insertion in sequence */
		rem = op % BLOCK;
		div = op / BLOCK;

		/*
		 * Remember edit info for later use with the annotations.
		 * Can't do this here as we need to process these backwards
		 * for complemented readings.
		 */
		eip->pos  = i+1+i2;
		eip->size = op;
		eip++;

		/* Realloc if required */
		if (*length + op >= *alength-1) {
		    *alength = *length + op + 100;
		    *seq_p  = (char *)xrealloc(*seq_p,  *alength * sizeof(**seq_p));
		    *conf_p = (int1 *)xrealloc(*conf_p, *alength * sizeof(**conf_p));
		    *opos_p = (int2 *)xrealloc(*opos_p, *alength * sizeof(**opos_p));
		    seq = *seq_p;
		    conf = *conf_p;
		    opos = *opos_p;
		}

		for (k=0; k < div; k++, i2 += BLOCK)
		    io_insert_seq(length, start, end, seq,
				  conf, opos, i+1+i2, pads, NULL, NULL,
				  BLOCK);
		if (rem) {
		    io_insert_seq(length, start, end, seq,
				  conf, opos, i+1+i2, pads, NULL, NULL,
				  rem);
		    i2 += rem;
		}

		j += op;
		op = 0;
	    } else {
		/* insertion in consensus */
		pad_consensus(io, contig, j+1+j2, -op);
		i -= op;
		j2 -= op;
		op = 0;
	    }
	}
    }

    /*
     * Adjust annotation positions
     */
    if (eip != ei) {
	edit_info *e = ei;

	if (comp) {
	    for (e = eip-1; e >= ei; e--) {
		insert_to_anno_list(anno_r, anno_rl,
				    *length - (e->pos + e->size) + 2, e->size);
		insert_to_anno_list(anno_c, anno_cl,
				    *length - (e->pos + e->size) + 2, e->size);
	    }
	} else {
	    for (e = ei; e < eip; e++) {
		insert_to_anno_list(anno_r, anno_rl, e->pos, e->size);
		insert_to_anno_list(anno_c, anno_cl, e->pos, e->size);
	    }
	}
    }

    xfree(ei); 
}


/*
 * Enters a reading. This does not perform any linking of neighbours and
 * updating of the contig. If an alignment buffer 'res' is non NULL then
 * it is used as a list of edits. See the calign() code for details of
 * the format for this buffer.
 *
 * It does however edit and enter the tags (on both the reading and
 * consensus).
 *
 * Returns reading number for success, -1 for failure.
 */
int enter_reading(GapIO *io, SeqInfo *si, int comp, align_info *ai,
		  int contig, int position) {
    GReadings r;
    int reading;
    char *name;
    int start, end, length, alength;
    char *seq = NULL;
    int1 *conf = NULL;
    int2 *opos = NULL;
    anno_info *anno_r, *anno_c;
    int anno_rl, anno_cl;
    int retcode = -1;
    int i;

    /*
     * Allocate
     */
    io_init_reading(io, NumReadings(io)+1);
    reading = NumReadings(io);


    /*
     * Write the reading name
     */
    if (NULL == (name = read_sequence_name(si)))
	goto end;
    write_rname(io, reading, name);


    /*
     * length, start, end, sense
     */
    length = si->length;
    start = si->start;
    end = si->end;
    alength = length + 100;

    /* Allocate temporary buffers */
    seq = (char *)xmalloc(alength * sizeof(*seq));
    conf = (int1 *)xmalloc(alength * sizeof(*conf));
    opos = (int2 *)xmalloc(alength * sizeof(*opos));
    if (!seq || !conf || !opos)
	goto end;
    
    /*
     * Obtain the sequence, original positions, and confidence values.
     */
    strcpy(seq, exp_get_entry(si->e, EFLT_SQ));
    SeqInfo_opos(si, opos, length);
    SeqInfo_conf(si, conf, length);


    /*
     * Complement if necessary
     */
    if (comp)
	io_complement_seq(&length, &start, &end, seq, conf, opos);


    /*
     * Load the tags into memory
     */
    anno_r = create_anno_list(si, EFLT_TG, &anno_rl, length);
    anno_c = create_anno_list(si, EFLT_TC, &anno_cl, 0);
    
    /*
     * Edit the sequence, consensus and tags
     */
    if (ai)
	edit_sequence(io, ai, &alength, &length, &start, &end,
		      &seq, &conf, &opos, contig,
		      comp, anno_r, anno_rl, anno_c, anno_cl);


    /*
     * write sequence to file
     */
    if (io_write_seq(io, reading, &length, &start, &end, seq, conf, opos)) {
	verror(ERR_WARN, "enter_reading",
	       "Problem writing new sequence to database: %s", name);
	return -1;
    }
    /* Shouldn't this be in io_write_seq? */
    gel_read(io, reading, r);
    r.sequence_length = end - start - 1;
    if (comp) {
	io_length(io, reading) = -r.sequence_length;
	r.sense = GAP_SENSE_REVERSE;
    } else {
	io_length(io, reading) = r.sequence_length;
	r.sense = GAP_SENSE_ORIGINAL;
    }
    gel_write(io, reading, r);


    /*
     * write raw data info
     */
    if (exp_Nentries(si->e, EFLT_LT) &&
	exp_Nentries(si->e, EFLT_LN)) {
	if (io_write_rd(io,
			reading,
			exp_get_entry(si->e, EFLT_LN),
			strlen(exp_get_entry(si->e, EFLT_LN)),
			exp_get_entry(si->e, EFLT_LT),
			strlen(exp_get_entry(si->e, EFLT_LT)))) {
	    verror(ERR_WARN, "enter_reading",
		   "Problem writing raw data information to database: %s",
		   name);
	    return -1;
	}
    }

    /*
     * Write the tag data
     */
    write_anno_list(io, anno_r, anno_rl, reading, 0, comp,
		    start, end, length);
    write_anno_list(io, anno_c, anno_cl, -contig, position-1-r.start,
		    comp, start, end, length);
    free_anno_list(anno_r, anno_rl);
    free_anno_list(anno_c, anno_cl);


    /*
     * Add NoTes
     */
    for(i = 0; i < exp_Nentries(si->e, EFLT_NT); i++) {
	create_note_for_gel(io, reading,
			    arr(char *, si->e->entries[EFLT_NT], i));
    }

    /*
     * write everything else
     */
    retcode = add_seq_details(io, reading, si) ? -1 : reading;

 end:
    if (seq)
	xfree(seq);
    if (conf)
	xfree(conf);
    if (opos)
	xfree(opos);

    return retcode;
}


/*
 * Links a reading into a contig at a specified position, updating the
 * necessary information to maintain consistency.
 *
 * A positon <= 0 implies that this reading should go leftwards of the first
 * reading in this  contig. In this case we shift everything else rightwards.
 *
 * Returns 0 for success, -1 for failure.
 */
int link_reading(GapIO *io, int from_rnum, int rnum,
		 int contig, int position) {
    GContigs c;
    GReadings r, rt;
    int right = 0, left = 0;

    if (gel_read(io, rnum, r) || contig_read(io, contig, c))
	return -1;

    io_relpos(io, rnum) = r.position = position;

    /*
     * Chain left/right from anchor seq looking for a suitable place to insert.
     */
    if (from_rnum == 0)
	from_rnum = io_clnbr(io, contig);
    if (position >= io_relpos(io, from_rnum)) {
	left = io_lnbr(io, from_rnum);
	for (right = from_rnum; right; right = io_rnbr(io, right)) {
	    if (io_relpos(io, right) > position)
		break;
	    left = right;
	}
    } else {
	right = io_rnbr(io, from_rnum);
	for (left = from_rnum; left; left = io_lnbr(io, left)) {
	    if (io_relpos(io, left) <= position)
		break;
	    right = left;
	}
    }

    /*
     * Link the neighbours
     */
    if (left) {
	gel_read(io, left, rt);
	io_rnbr(io, left) = rt.right = rnum;
	gel_write(io, left, rt);
    } else {
	io_clnbr(io, contig) = c.left = rnum;
    }

    io_lnbr(io, rnum) = r.left = left;

    if (right) {
	gel_read(io, right, rt);
	io_lnbr(io, right) = rt.left = rnum;
	gel_write(io, right, rt);
    } else {
	io_crnbr(io, contig) = c.right = rnum;
    }

    io_rnbr(io, rnum) = r.right = right;

    /*
     * If this reading gets added onto the left end of a contig then
     * we need to shift all other readings right.
     */
    if (position < 1) {
	io_relpos(io, rnum) = r.position = 1;

	for (; right; right = io_rnbr(io, right)) {
	    gel_read(io, right, rt);
	    io_relpos(io, right) = rt.position = rt.position - position + 1;
	    gel_write(io, right, rt);
	}

	/* We also need to shift all tags right */
	shift_contig_tags(io, contig, 0, -position + 1);
    }

    /*
     * Update contig length. This is trivial when we've simply added a reading
     * somewhere in a contig other than the start, otherwise we need to
     * rescan.
     */
    if (position < 1) {
	int len = 0, end;

	for (right = io_clnbr(io, contig); right; right = io_rnbr(io, right)) {
	    end = io_relpos(io, right) + ABS(io_length(io, right));
	    if (end > len)
		len = end;
	}

	io_clength(io, contig) = c.length = len - 1;
    } else {
	if (r.sequence_length + r.position - 1 > c.length) {
	    io_clength(io, contig) = r.sequence_length + r.position - 1;
	    c.length = r.sequence_length + r.position - 1;
	}
    }

    /*
     * Write back
     */
    if (gel_write(io, rnum, r) || contig_write(io, contig, c))
	return -1;

    return 0;
}


/*
 * Clips a sequence to remove the end padding, updating positions and lengths
 * accordingly.
 * 'which' specifies which read to clip. 1==1st, 2==2nd, 3==both, 0==none.
 *
 * This also returns the percentage mismatch (a double between 0 and 100).
 */
static double clip_score(char *seq1, char *seq2,
			 int *len1, int *len2,
			 align_int *S,
			 int *pos, int *start1, int *start2,
			 int which) {
    register int i = 0, j = 0, op = 0, mat = 0;
    register int mism = 0, l1 = *len1, l2 = *len2;
    
    /*
     * Clip excess padding from left end.
     */
    if (*S > 0 && (which&1)) {
	*pos     = *S + *start2 + 1;
	*start2 += *S;
	seq2    += *S;
	l2      -= *S;
	memmove(S, S+1, (*len1 + *len2) * sizeof(*S));
    } else if (*S < 0 && (which&2)) {
	*pos     = *S + *start2 + 1;
	*start1 -= *S;
	seq1    -= *S;
	l1      += *S;
	memmove(S, S+1, (*len1 + *len2) * sizeof(*S));
    }

    /*
     * Count mismatches and find first sequence to end
     */
    while (i < l1 || j < l2) {
	if (((which&1) && i >= l1) || ((which&2) && j >= l2)) {
	    break;
	}

	if (op == 0 && *S == 0) {
	    op = *S++;
	    if (seq1[i++] != seq2[j++])
		mism++;
	    else
		mat++;
	} else {
	    if (op == 0)
		op = *S++;

	    if (op > 0)
		j+=op, mism+=op;
	    else
		i-=op, mism-=op;

	    op = 0;
	}
    }

    /*
     * Clip excess padding from right end
     */
    *len1 = i;
    *len2 = j;

    return 100 * ((mism + mat) ? ((double)mism / (mism + mat)) : 1);
}


/*
 * Produces an alignment buffer between sequence 'si' (in 'dir' direction)
 * and the contig at position 'pos', tolerance 'tol' and allowing for maximum
 * padding of 'maxpads'.
 *
 * Returns NULL for failure, otherwise an alloced buffer containing
 * the alignment edits.
 * 'ierr' is set to contain an error code. 0 means OK, 1 means fatal error,
 * 2 means failed mismatch score, 3 means no overlap.
 */
align_info *assemble_align(GapIO *io, SeqInfo *si, consen_info *ci,
			   int contig, int *pos, int dir, int tol,
			   int display, int maxpads, double max_mism,
			   int *ierr) {
    char *seq = NULL;
    int start, end, len, orig_start;
    int s_start, s_end, s_length;
    align_info *ai = NULL;
    double mism;
    int orig_pos = *pos;

    end = 0;
    *ierr = 0;

    if (NULL == (ai = (align_info *)xmalloc(sizeof(align_info)))) {
	*ierr = 1;
	goto error;
    }
    ai->res = NULL;

    s_length = si->length;
    s_start = si->start;
    s_end = si->end;
    seq = (char *)xmalloc(s_length * sizeof(*seq));
    strncpy(seq, exp_get_entry(si->e, EFLT_SQ), s_length);

    /* Complement if requested */
    if (dir == GAP_STRAND_REVERSE)
	io_complement_seq(&s_length, &s_start, &s_end, seq, NULL, NULL);

    orig_start = s_start;

    /* Calculate ranges to align */    
    start = *pos-1 - tol;

    if (start < 0) {
	if (-*pos - tol > 0) {
	    s_start += -*pos - tol;
	}
	start = 0;
	end = -1; /* flag */
    }

    if (s_start < 0)
	s_start = 0;

    if ((len = s_end - s_start - 1) < 1) {
#ifdef DEBUG
	printf("len %d s_end %d s_start %d\n", len, s_end, s_start);
#endif
	*ierr = 3;
	goto error;
    }

    /* KFB 14.09.01 FIXME - check with James
     * this is wrong if start is set to 0
     *    end   = *pos-1 + tol + maxpads + len;
     */
    if (end == 0) {
	end = *pos-1 + tol + maxpads + len;
    } else {
	end = tol + maxpads + len;
    }
#ifdef DEBUG
    printf("len %d end %d pos %d tol %d maxpads %d\n", len, end, *pos, tol, maxpads);
#endif

    if (end >= io_clength(io, contig)) {
	/* KFB 06.09.01 found end was one out */
	/* end = io_clength(io, contig)-1; */
	end = io_clength(io, contig);
    }

    if (end <= start) {
#ifdef DEBUG
	printf("end %d start %d\n", end, start);
#endif
	*ierr = 3;
	goto error;
    }


    /* Perform the alignment */
    if (NULL == (ai->res = (align_int *)xcalloc((len + end - start + 1),
						sizeof(align_int)))) {
	*ierr = 1;
	goto error;
    }

    if (-1 == calign(&seq[s_start], &ci->con_item[contig-1][start], /* seqs */
		     len, end - start,       /* lengths */
		     NULL, NULL, NULL, NULL, /* returned sequences */
		     0, 0,		     /* bands */
		     gopenval, gextendval,   /* penalties */
		     3, 0,		     /* job, protein */
		     ai->res)) {             /* result */
	*ierr = 1;
	goto error;
    }

    ai->start1 = s_start;
    ai->start2 = start;
    ai->len1 = len;
    ai->len2 = end - start;
    
    /*
     * Clip ends and derive score
     */
    mism = clip_score(&seq[ai->start1], &ci->con_item[contig-1][ai->start2],
 		      &ai->len1, &ai->len2, ai->res,
		      pos, &ai->start1, &ai->start2, 3 /* clip both ends */);

    if (display) {
	char *seq1, *seq2;
	int seq1l, seq2l;
	if (NULL == (seq1 = (char *)xmalloc(2 * s_length * sizeof(*seq))))
	    goto error;

	if (NULL == (seq2 = (char *)xmalloc(2 * s_length * sizeof(*seq)))) {
	    xfree(seq1);
	    goto error;
	}

	cexpand(&seq[ai->start1], &ci->con_item[contig-1][ai->start2],
		ai->len1, ai->len2,
		seq1, seq2, &seq1l, &seq2l,
		3 | ALIGN_J_PADS, ai->res);
	
	/* KFB: out by one error
	 * list_alignment(seq1, seq2, "Reading", "Consensus",
	 *       ai->start1 - orig_start + 1, ai->start2, "");
	 */
	list_alignment(seq1, seq2, "Reading", "Consensus",
		       ai->start1 - orig_start + 1, ai->start2+1, "");

	xfree(seq1);
	xfree(seq2);
    }

    if (ABS(orig_pos - (*pos - (s_start - orig_start))) > tol) {
#ifdef DEBUG
	printf("TOL %d %d orig_pos %d pos %d s_start %d orig_start %d\n",
	       ABS(orig_pos - (*pos - (s_start - orig_start))), tol, 
	       orig_pos, *pos, s_start, orig_start);
#endif
	*ierr = 4;
	goto error;
    }

    if (max_mism >= 0 && mism > max_mism) {
	*ierr = 2;
	goto error;
    }

    /* was *pos -= s_start - orig_start; */
    *pos = ai->start2 - (ai->start1 - orig_start) + 1;

    xfree(seq);
    return ai;

 error:
    if (ai->res)
	xfree(ai->res);
    if (ai)
	xfree(ai);
    if (seq)
	xfree(seq);

    return NULL;
}


char *assemble_direct(GapIO *io, int display, double max_mism,
		      char *inlist, int do_alignments, int enter_all) {
    consen_info *ci = NULL;
    SeqInfo *si;
    int ierr;
    char *file;
    int name, dir, pos, tol, maxpads = 25;
    int contig, rnum;
    align_info *ai;
    char *rname;
    int old_clen, count = 1;
    void *dl;
    char *res;
    
    dl = alloc_dlist();

    if (do_alignments)
	/* Calculate consensus */
	if (NULL == (ci = all_consensus(io, consensus_cutoff)))
	    return NULL;

    /* Loop around each input file */
    if (-1 == set_active_list(inlist))
        return NULL;
    while (file = get_active_list_item()) {
	UpdateTextOutput();
	vmessage("Processing number %8d: %s\n", count++, file);
	
	/*
	 * The last argument controls whether to set the QL/QR based on
	 * SL/SR. "1" disables SL/SR use, which is what we need when copying
	 * data in generated by extract_seq. FIXME: this could break things
	 * for other data sources though, unless explicit QL/QR values have
	 * been used.
	 */
	if (NULL == (si = read_sequence_details(file, 1))) {
	    verror(ERR_WARN, "directed_assembly", "couldn't read '%s'", file);
	    add_to_dlist(dl, file);
	    vmessage("  failed\n");
	    continue;
	}
	if (si->start < 0)
	    si->start = 0;
	if (si->end > si->length+1)
	    si->end = si->length+1;

	if (si->end - si->start <= 1 || si->length < 1) {
	    verror(ERR_WARN, "directed assembly", "sequence '%s' too short",
		   file);
	    add_to_dlist(dl, file);
	    freeSeqInfo(si);
	    vmessage("  failed\n");
	    continue;
	}

	if (NULL == (rname = read_sequence_name(si))) {
	    verror(ERR_WARN, "directed_assembly", "no name found for '%s'",
		   file);
	    add_to_dlist(dl, file);
	    freeSeqInfo(si);
	    vmessage("  failed\n");
	    continue;
	}

	/* FIXME: perform caching of names (also template and vector names) */
	if (get_gel_num(io, rname, GGN_NAME) > 0) {
	    verror(ERR_WARN, "directed_assembly",
		   "reading '%s' already exists", rname);
	    add_to_dlist(dl, file);
	    vmessage("  failed\n");
	    freeSeqInfo(si);
	    continue;
	}

	if (exp_Nentries(si->e, EFLT_AP) == 0) {
	    verror(ERR_WARN, "directed_assembly", "no AP line in '%s'", file);
	    add_to_dlist(dl, file);
	    vmessage("  failed\n");
	    freeSeqInfo(si);
	    continue;
	}

	if (0 != AP_parse(io, rname, exp_get_entry(si->e, EFLT_AP),
			  &name, &dir, &pos, &contig, &tol)) {
	    verror(ERR_WARN, "directed_assembly", "invalid AP line in '%s'.",
		   file);
	    name = 0;

	    if (!enter_all) {
		verror(ERR_WARN, "directed_assembly",
		       "invalid AP line in '%s'", file);
		add_to_dlist(dl, file);
		vmessage("  failed\n");
		freeSeqInfo(si);
		continue;
	    }
	}

	ai = NULL;

	if (name != 0) {
	    if (tol >= 0 && do_alignments) {
		/* Align with existing contig */
		ai = assemble_align(io, si, ci, contig, &pos, dir, tol,
				    display, maxpads, max_mism, &ierr);
		if (NULL == ai) {
		    if (ierr == 2) {
			vmessage("  Percentage mismatch is too high\n");
			name = 0;
		    } else if (ierr == 3) {
			vmessage("  Reading does not overlap\n");
			name = 0;
		    } else if (ierr == 4) {
			vmessage("  Reading does not overlap "
				 "within tolerance\n");
			name = 0;
		    } else {
			verror(ERR_WARN, "directed_assembly",
			       "failed to align reading '%s'", file);
			name = 0;
		    }

		    if (!enter_all) {
			add_to_dlist(dl, file);
			freeSeqInfo(si);
			continue;
		    }
		}
	    }
	}

	if (name == 0) {
	    /* new contig */
	    vmessage("  Creating new contig\n");
	    if (-1 == io_init_contig(io, NumContigs(io)+1))
		return NULL;

	    contig=NumContigs(io);
	    io_clnbr(io, contig) = 0;
	    io_crnbr(io, contig) = 0;
	    io_clength(io, contig) = 0;
	    pos = 0;
	}

	old_clen = io_clength(io, contig);
	if (-1 == (rnum = enter_reading(io, si, dir, ai, contig, pos)))
	    return NULL;
	
	if (-1 == link_reading(io, name, rnum, contig, pos))
	    return NULL;

	freeSeqInfo(si);
	if (ai) {
	    xfree(ai->res);
	    xfree(ai);
	}

	if (do_alignments) {
	    if (-1 == recalc_consensus(io, ci, contig, pos,
				       ABS(io_length(io, rnum)),
				       old_clen, io_clength(io, contig))) {
		verror(ERR_WARN, "directed_assembly",
		       "failed to recalculate consensus - quitting");
		return NULL;
	    }
	}
    }

    remove_contig_holes_all(io);

    flush2t(io);
    if (do_alignments)
	free_all_consensus(ci);

    res = strdup(read_dlist(dl));
    free_dlist(dl);

    return res;
}

/*
 * ---------------------------------------------------------------------------
 * Contig tidying functions.
 *
 * If we break the logical consistency of the databse then these will try to
 * fix it. This is useful to do for simplicities sake. For example to enter
 * directed assembly data we can just slurp up all the data and then scan
 * through afterwards to make sure we have no contigs and that the contig
 * starts at base 1.
 */

/*
 * remove_contig_holes - checks for gaps in a contig at the start, end or
 * internal. Internal holes require splitting the contig in two, while start
 * and end gaps simply require adjustment of the length and shuffling sequence
 * positions.
 *
 * This function obtains all information from the disk (or cached) database 
 * structures rather than the internal arrays incase these are not
 * consistent (yet). It does however assume that the reading positions
 * are sorted left to right.
 *
 * Returns 0 for success
 *        -1 for error
 */
int remove_contig_holes(GapIO *io, int cnum) {
    int rnum, prev_rnum, new_contig;
    int furthest;
    int cstart;
    int shift;
    GContigs c;
    GReadings r;
    
    vfuncheader("remove_contig_holes()");

    if (contig_read(io, cnum, c))
	return -1;

    do {
	new_contig = 0;
	
	rnum = c.left;
	prev_rnum = 0;
	cstart = 1;
	furthest = 1;

	while (rnum) {
	    if (gel_read(io, rnum, r))
		return -1;
	    
	    if (cstart) {
		if (r.position != 1)
		    vmessage("Gap at start, shifting by %d bases\n", shift);
		shift = r.position - 1;
	    }
	    
	    r.position -= shift;
	    io_relpos(io, rnum) -= shift;

	    /* If there's a gap - start a new contig */
	    if (!cstart && r.position > furthest) {
		vmessage("Hole from %d to %d, breaking contig\n",
			 furthest, r.position);
		new_contig = 1;
		break;
	    }

	    /* keep track of rightmost sequence end, to spot gaps */
	    if (furthest < r.position + r.sequence_length - 1)
		furthest = r.position + r.sequence_length - 1;

	    if (shift) {
		gel_write(io, rnum, r);
		io_relpos(io, rnum) = r.position;
	    }

	    prev_rnum = rnum;
	    rnum = r.right;
	    cstart = 0;
	}

	/* Update contig size etc */
	c.length = furthest;
	c.right = prev_rnum;
	contig_write(io, cnum, c);
	io_crnbr(io, cnum) = c.right;
	io_clength(io, cnum) = c.length;

	if (new_contig) {
	    if (-1 == io_init_contig(io, cnum = NumContigs(io)+1))
		return -1;

	    contig_read(io, cnum, c);
	    c.left = rnum;
	    io_clnbr(io, cnum) = c.left;

	    /* Terminate read link list for existing contig */
	    gel_read(io, prev_rnum, r);
	    r.right = 0;
	    io_rnbr(io, prev_rnum) = 0;
	    gel_write(io, prev_rnum, r);

	    /* Start read link list for new contig */
	    gel_read(io, rnum, r);
	    r.left = 0;
	    io_lnbr(io, rnum) = 0;
	    gel_write(io, rnum, r);
	}
    } while (new_contig);

    return 0;
}

/*
 * Calls remove_contig_holes for all contigs in the DB.
 * Returns the same codes.
 */
int remove_contig_holes_all(GapIO *io) {
    int i, ret = 0;
    for (i = 1; i <= NumContigs(io); i++) {
	ret |= remove_contig_holes(io, i);
    }

    return ret;
}
