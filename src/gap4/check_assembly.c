#include "IO.h"
#include "io-reg.h"
#include "cs-object.h"
#include "gap_globals.h"
#include "misc.h"
#include "qual.h"
#include "align.h"
#include "contig_selector.h"
#include "contigEditor.h"
#include "dna_utils.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"
#include "consen.h"

/*
 * Checks the left end of an alignment and clips to the point where both
 * sequences start. Ie it relaces:
 *
 * *****AGGATTT
 * GGGGCAGGATTT
 *
 * with:
 *
 * AGGATTT
 * AGGATTT
 *
 * Modifies seq1, seq2, len1, len2 and the contents of S.
 * Returns the distance shifted (0 for none, >0 for above example, <0 for
 * opposite case).
 */
int align_clip_left(char **seq1, char **seq2, int *len1, int *len2,
		    align_int *S) {
    if (*S > 0) {
	*seq2 += *S;
	*len2 -= *S;
	memmove(S, S+1, (*len1 + *len2) * sizeof(*S));
    } else if (*S < 0) {
	*seq1 -= *S;
	*len1 += *S;
	memmove(S, S+1, (*len1 + *len2) * sizeof(*S));
    }

    return *S;
}

/*
 * Checks the right end of an alignment and clips to the point where
 * both sequences end. Ie it replaces:
 *
 * AGGATTT*****
 * AGGATTTCGGGG
 *
 * with:
 *
 * AGGATTT
 * AGGATTT
 *
 * Modifies only len1, len2.
 * Returns the amount removed (0 for none, >0 for above example, <0 for
 * opposite case).
 */
int align_clip_right(char **seq1, char **seq2, int *len1, int *len2,
		     align_int *S) {
    register int i = 0, j = 0, op;
    register int l1 = *len1, l2 = *len2;

    while (i < l1 && j < l2) {
	if ((op = *S++) == 0)
	    i++, j++;
	else if (op > 0)
	    j += op;
	else
	    i -= op;
    }

    *len1 = i;
    *len2 = j;

    return l1 == i ? l2 - j : i - l1;
}


/*
 * Returns the percentage mismatch of the alignment between seq1 and seq2.
 * Pads vs pads are not counted within the mismatch.
 *
 * Also fills in the mism_p and mat_p pointers if non NULL.
 */
double align_score(char *seq1, char *seq2, int len1, int len2,
		  int *mism_p, int *mat_p, align_int *S) {
    register int i = 0, j = 0, op;
    register int l1 = len1, l2 = len2;
    register int mism = 0, len = 0;

    while (i < l1 || j < l2) {
	if ((op = *S++) == 0) {
	    if (seq1[i++] != seq2[j++])
		mism++;
	    len++;
	} else if (op > 0) {
	    len += op;
	    do {
		if (seq2[j++] != '*')
		    mism++;
	    } while (--op);
	} else {
	    len -= op;
	    do {
		if (seq1[i++] != '*')
		    mism++;
	    } while (++op);
	}
    }

    if (mism_p) *mism_p = mism;
    if (mat_p)  *mat_p = len - mism;

    return (double)mism / len;
}

/*
 * Checks a single reading for correct assembly by analysing the used data.
 *
 * Returns -1 for system error, otherwise a score (0-1000000)
 * 'pos_p' and 'len_p' are filled in with the position and length of the match
 * within the consensus.
 */
int check_uassembly_single(GapIO *io, char *con, int contig, int rn,
			   int *pos_p, int *len_p, float maxperc,
			   int win_len) {
    int length, start, end;
    char *seq = NULL, tmp;
    int i, j, mism = 0;
    int worst, worst_pos = -1;

    /* Get sequence */
    if (-1 == io_aread_seq(io, rn, &length, &start, &end, &seq, NULL, NULL, 1))
    {
	if (seq)
	    xfree(seq);
	return -1;
    }

    /* Initialise scoring for divergence checks */
    if (end - start - 1 < win_len) {
	win_len = end - start - 1;
    }
    worst = 0.5 + maxperc * win_len;

    for (i=start, j=io_relpos(io, rn)-1; i < start + win_len; i++, j++) {
	if (!same_char(seq[i], con[j]))
	    mism++;
    }

    /* Loop through sequence looking for problems */
    do {
	if (mism >= worst) {
	    worst_pos = i;
	    worst = mism;
	}

	mism -= !same_char(seq[i++-win_len], con[j++-win_len]);
	if (i < end-2)
	    mism += !same_char(seq[i], con[j]);
    } while (i < end - 1);

    /* Display problem, listing worst score */
    if (worst_pos != -1) {
	*pos_p = io_relpos(io, rn);
	*len_p = end - start + 1;

	vmessage("\nReading %d(%s) has a local percentage mismatch of %2.1f\n",
		 rn, get_read_name(io, rn),
		 100 * (float)worst / win_len);

	/* Display alignment (ensuring it's only the used portion) */
	seq[end-1] = 0;
	tmp = con[io_relpos(io, rn)+end-start-2];

	con[io_relpos(io, rn)+end-start-2] = 0;

	list_alignment(&seq[start], &con[io_relpos(io, rn)-1],
		       "Reading", "Consensus",
		       1, io_relpos(io, rn), "");
	con[io_relpos(io, rn)+end-start-2] = tmp;

	xfree(seq);
	return 10000 * (float)worst / win_len;
    }

    xfree(seq);
    return 0;
}

/*
 * Checks a single reading for correct assembly by analysing the cutoff
 * data.
 *
 * Returns -1 for system error, otherwise a score (0-1000000)
 * 'pos_p' and 'len_p' are filled in with the position and length of the match
 * within the consensus.
 */
int check_cassembly_single(GapIO *io, char *con, int contig, int rn,
			   int *pos_p, int *len_p,
			   int minlen, int winsize, int maxdash,
			   float maxperc) {
    char *cutoff, *cutoffp;
    int cutlen, cutlenc, pos, sense, mism, mat;
    double score;
    int job;
    GReadings r;
    align_int *res = NULL;

    gel_read(io, rn, r);
    cutoffp = cutoff = (char *)xmalloc(r.length + 1);
    if (NULL == cutoff)
	return -1;

    /* Obtain extension */
    if (-1 == io_get_extension(io, rn, cutoff, r.length, &cutlen, &sense)) {
	xfree(cutoff);
	return -1;
    }

    if (cutlen < minlen) {
	xfree(cutoff);
	return 0;
    }

    /* Deduce region of reading to slide window along */
    pos = io_relpos(io, rn) +
	(io_length(io, rn) < 0 ? -cutlen : io_length(io, rn) - 1);

    if (pos + cutlen > io_clength(io, contig)) {
	cutlen = io_clength(io, contig) - pos;
    } else if (pos < 1) {
	cutlen -= 1-pos;
	cutoffp += 1-pos;
	pos = 0;
    }

    if (cutlen < minlen) {
	xfree(cutoff);
	return 0;
    }

    /* Slide window and deduce region of consensus/reading to align */
    if (io_length(io, rn) < 0) {
	char *seq = (char *)xmalloc(cutlen+1);
	int e;

	if (!seq) {
	    xfree(cutoff);
	    return -1;
	}

	strncpy(seq, cutoffp, cutlen);
	seq[cutlen] = 0;
	complement_seq(seq, cutlen);
	e = end_of_good(seq, 1, winsize, maxdash);
	cutoffp += cutlen - e;
	pos += cutlen - e;
	cutlen = e;

	/* Give the consensus a little bit more to align with */
	pos -= winsize / maxdash + 1;
	cutlenc = cutlen + winsize / maxdash + 1;
	if (pos < 0) {
	    cutlenc -= -pos;
	    pos = 0;
	}

	xfree(seq);
    } else {
	cutlen = end_of_good(cutoffp, 1, winsize, maxdash);

	/* Give the consensus a little bit more to align with */
	cutlenc = cutlen + winsize / maxdash + 1;
	if (cutlenc > io_clength(io, contig) - pos)
	    cutlenc = io_clength(io, contig) - pos;
    }

    if (cutlen < minlen) {
	xfree(cutoff);
	return 0;
    }

    /* Align */
    cutoffp[cutlen] = 0;
    con += pos;
    job = 3 | ((io_length(io,rn) < 0)
	       ? (ALIGN_GAP_E1 | ALIGN_GAP_E2)
	       : (ALIGN_GAP_S1 | ALIGN_GAP_S2));
    /* calloc in order to stop workshop producing false warnings */
    if (NULL == (res = (align_int *)xcalloc((cutlen + cutlenc + 1),
					    sizeof(align_int)))) {
	xfree(cutoff);
	return -1;
    }
    calign(cutoffp, con, cutlen, cutlenc,
	   NULL, NULL, NULL, NULL,
	   0, 0, gopenval, gextendval, job, 0,
	   res);

    if (io_length(io, rn) >= 0)
	align_clip_right(&cutoffp, &con, &cutlen, &cutlenc, res);
    else
	align_clip_left(&cutoffp, &con, &cutlen, &cutlenc, res);

    score = align_score(cutoffp, con, cutlen, cutlenc, &mism, &mat, res);

    if (score > maxperc) {
	char *seq1, *seq2;
	int seq1l, seq2l;

	seq1 = (char *)xmalloc(cutlen + cutlenc + 1);
	seq2 = (char *)xmalloc(cutlen + cutlenc + 1);
	if (!seq1 || !seq2) {
	    if (seq1) xfree(seq1);
	    if (seq2) xfree(seq1);
	    xfree(cutoff);
	    xfree(res);
	    return -1;
	}

	vmessage("\nReading %d(%s) has percentage mismatch of %2.1f\n",
		 rn, get_read_name(io, rn), score * 100);

	cexpand(cutoffp, con, cutlen, cutlenc,
		seq1, seq2, &seq1l, &seq2l,
		ALIGN_J_SSH | ALIGN_J_PADS, res);

	list_alignment(seq1, seq2, "Reading", "Consensus",
		       io_length(io, rn) > 0
		       ? io_length(io, rn) : (-seq1l + 1),
		       pos+1, "");

	*len_p = mism + mat;
	*pos_p = pos+1;

	xfree(seq1);
	xfree(seq2);
	xfree(cutoff);
	xfree(res);
	return (int) (score * 10000);
    }

    xfree(cutoff);
    xfree(res);
    return 0;
}


/*
 * Match callback.
 * 'obj' is a match contained within the check assembly list.
 */
void *checkass_obj_func(int job, void *jdata, obj_checkass *obj,
			mobj_checkass *ca) {
    static char buf[80];
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(ca->io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(ca->io, cs_id, 0);

    switch(job) {
    case OBJ_LIST_OPERATIONS:
	return "Information\0Hide\0Invoke contig editor *\0"
	    "SEPARATOR\0Remove\0";

    case OBJ_INVOKE_OPERATION:
	switch(*((int *)jdata)) {
	case 0: /* Information */
	    vfuncgroup(1, "2D plot matches");
	case -1: /* Information from results manager */
	    start_message();
	    if (ca->cutoffs) {
		vmessage("check_assembly match: hidden data\n");
	    } else {
		vmessage("check_assembly match: used data\n");
	    }
	    vmessage("    From contig %s(#%d,%d) at position %d\n",
		     get_contig_name(ca->io, ABS(obj->c1)),
		     io_clnbr(ca->io, ABS(obj->c1)), obj->c1,obj->pos1);
	    vmessage("    From reading %s(#%d) at position %d\n",
		     get_read_name(ca->io, obj->read),
		     obj->read, obj->pos1 - io_relpos(ca->io, obj->read));
	    vmessage("    Length %d, mismatch %2.2f%%\n\n",
		     obj->length, ((float)obj->score)/10000);
	    end_message(cs->window);
	    break;

	case 1: /* Hide */
	    obj_hide(GetInterp(), cs->window, obj,
		     (mobj_repeat *)ca, csplot_hash);
	    break;

	case -2: /* default */
	case 2: /* Invoke contig editor */ {
	    int cnum, llino, pos, id;

	    obj->flags |= OBJ_FLAG_VISITED;
	    ca->current = obj - ca->match;
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(ca), NULL);

	    cnum  = abs(obj->c1);
	    llino = obj->read;
	    pos   = obj->pos1 - io_relpos(ca->io, llino);
	    if (pos < 1)
		pos = 1;
	    if (pos > ABS(io_length(ca->io, llino)))
		pos = ABS(io_length(ca->io, llino));

	    if ((id = editor_available(cnum, 1)) != -1) {
		move_editor(id, llino, pos);
	    } else {
		edit_contig(GetInterp(), ca->io, cnum, llino, pos,
			    consensus_cutoff, quality_cutoff, ca->cutoffs, 
			    NULL);
	    }
	    break;
	}

	case 3: /* Remove */
	    obj_remove(GetInterp(), cs->window, obj,
		       (mobj_repeat *)ca, csplot_hash);
	    break;

        }
        break;

    case OBJ_GET_BRIEF:
	sprintf(buf,
		"check_assembly: %c#%d@%d len %d, mis %2.2f%%",
		io_length(ca->io, obj->read) > 0 ? '+' : '-',
		obj->read, obj->pos1,
		obj->length, ((float)obj->score)/10000);
	return buf;
    }

    return NULL;
}


static int sort_func(const void *p1, const void *p2) {
    obj_checkass *m1 = (obj_match *)p1, *m2 = (obj_checkass *)p2;
    return m2->score - m1->score;
}

void check_assembly_callback(GapIO *io, int contig, void *fdata,
			     reg_data *jdata) {
    mobj_checkass *r = (mobj_checkass *)fdata;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id, 0);

    switch (jdata->job) {

    case REG_QUERY_NAME:

	sprintf(jdata->name.line, "Check Assembly");
	break;


    case REG_JOIN_TO:

	csmatch_join_to(io, contig, &jdata->join, (mobj_repeat *)r,
			csplot_hash, cs->window);
	break;


    case REG_COMPLEMENT:

	csmatch_complement(io, contig, (mobj_repeat *)r, csplot_hash,
			   cs->window);
	break;


    case REG_GET_OPS:

	if (r->all_hidden)
	    jdata->get_ops.ops = "PLACEHOLDER\0PLACEHOLDER\0Information\0"
		"PLACEHOLDER\0Hide all\0Reveal all\0Sort Matches\0"
		    "SEPARATOR\0Remove\0";
	else
	    jdata->get_ops.ops = "Use for 'Next'\0Reset 'Next'\0Information\0"
		"Configure\0Hide all\0Reveal all\0Sort Matches\0"
		    "SEPARATOR\0Remove\0";
	break;


    case REG_INVOKE_OP:

	switch (jdata->invoke_op.op) {
	case 0: /* Next */
	    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(r), NULL);
	    break;
	case 1: /* Reset Next */
	    csmatch_reset_next((mobj_repeat *)r);
	    break;
	case 2: /* Information */
	    csmatch_info((mobj_repeat *)r, "Check Assembly");
	    break;
	case 3: /* Configure */
	    csmatch_configure(io, cs->window, (mobj_repeat *)r);
	    break;
	case 4: /* Hide all */
	    csmatch_hide(GetInterp(), cs->window, (mobj_repeat *)r,
			 csplot_hash);
	    break;
	case 5: /* Reveal all */
	    csmatch_reveal(GetInterp(), cs->window, (mobj_repeat *)r,
			   csplot_hash);
	    break;
	case 6: /* Sort */
	    qsort(r->match, r->num_match, sizeof(obj_checkass), sort_func);
	    csmatch_reset_hash(csplot_hash, (mobj_repeat *)r);
	    r->current = -1;
	    break;
	case 7: /* Remove */
	    csmatch_remove(io, cs->window,
			   (mobj_repeat *)r,
			   csplot_hash);
	    break;
	}
	break;


    case REG_PARAMS:

	jdata->params.string = r->params;
	break;


    case REG_NUMBER_CHANGE:

	csmatch_renumber(io, contig, jdata->number.number,
			 (mobj_repeat *)r, csplot_hash, cs->window);
	break;


    case REG_ORDER:

	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;

    case REG_QUIT:

	csmatch_remove(io, cs->window,
		       (mobj_repeat *)r,
		       csplot_hash);
	break;

    case REG_DELETE:

	csmatch_contig_delete(io, (mobj_repeat *)r, contig,
			      cs->window, csplot_hash);
	break;

    case REG_LENGTH:
	csmatch_replot(io, (mobj_repeat *)r, csplot_hash, cs->window);
	break;
    }
}

int check_assembly_plot(GapIO *io, int *reads, int *conts, int *score,
			int *pos, int *length, int count, int cutoffs) {
    int i, id;
    mobj_checkass *ca;
    obj_checkass *matches;
    char *val;

    if (count == 0)
	return 0;

    if (NULL == (ca = (mobj_checkass *)xmalloc(sizeof(mobj_checkass)))) {
	return -1;
    }

    if (NULL == (matches = (obj_checkass *)xmalloc(count *
						   sizeof(obj_checkass)))) {
	xfree(ca);
	return -1;
    }

    /* Setup the meta-object */
    ca->num_match = count;
    ca->match = (obj_match *)matches;
    ca->io = io;
    ca->cutoffs = cutoffs;
    strcpy(ca->tagname, CPtr2Tcl(ca));

    val = get_default_string(GetInterp(), gap_defs, "CHECK_ASSEMBLY.COLOUR");
    strcpy(ca->colour, val);
    ca->linewidth = get_default_int(GetInterp(), gap_defs,
				    "CHECK_ASSEMBLY.LINEWIDTH");

    ca->params = (char *)xmalloc(100);
    if (ca->params)
	sprintf(ca->params, "Unknown at present");
    ca->all_hidden = 0;
    ca->current = -1;
    ca->current = -1;
    ca->reg_func = check_assembly_callback;
    ca->match_type = REG_TYPE_CHECKASS;

    /* Set up each object */
    for (i=0; i<count; i++) {
	matches[i].func = (void *(*)(int, void *, struct obj_match_t *,
				    struct mobj_repeat_t *))checkass_obj_func;
	matches[i].data = (mobj_repeat *)ca;
	matches[i].c1 = matches[i].c2 = conts[i];
	matches[i].pos1 = matches[i].pos2 = pos[i];
	matches[i].length = length[i];
	matches[i].score = score[i];
	matches[i].flags = 0;
	matches[i].read = reads[i];
    }

    /* Sort matches */
    qsort(ca->match, ca->num_match, sizeof(obj_match), sort_func);

    PlotRepeats(io, (mobj_repeat *)ca);
    Tcl_VarEval(GetInterp(), "CSLastUsed ", CPtr2Tcl(ca), NULL);

    /*
     * Register the repeat search with each of the contigs used.
     * Currently we assume that this is all.
     */
    id = register_id();
    for (i = 1; i <= NumContigs(io); i++) {
	contig_register(io, i, check_assembly_callback, (void *)ca, id,
			REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
			REG_NUMBER_CHANGE | REG_ORDER, REG_TYPE_CHECKASS);
    }

    return 0;
}

/*
 * Scans through one or more contigs checking each reading for correct
 * assembly.
 * This is done by obtaining the cutoff data and aligning this with
 * the consensus sequence.
 *
 * 'cutoff' determines whether we're checking the cutoff data or the used
 * data against the consensus sequence. The minlen, winsize and maxdash
 * parameters are only needed when cutoff == 1.
 *
 * Returns -1 for failure, 0 for success.
 */
int check_assembly(GapIO *io, int num_contigs, int *contigs, int cutoff,
		   int minlen, int winsize, int maxdash, float maxperc) {
    int i, rn, sc, count = 0;
    char *con;
    int *reads = NULL, *conts = NULL, *score = NULL, *length = NULL;
    int *pos = NULL;

    if (NULL == (reads = (int *)xcalloc(NumReadings(io), sizeof(int))))
	goto error;

    if (NULL == (conts = (int *)xcalloc(NumReadings(io), sizeof(int))))
	goto error;

    if (NULL == (score = (int *)xcalloc(NumReadings(io), sizeof(int))))
	goto error;

    if (NULL == (length = (int *)xcalloc(NumReadings(io), sizeof(int))))
	goto error;

    if (NULL == (pos = (int *)xcalloc(NumReadings(io), sizeof(int))))
	goto error;

    for (i = 0; i < num_contigs; i++) {
	int cnum = contigs[i], len, posn;

	if (NULL == (con = (char *)xmalloc(io_clength(io, cnum)+1)))
	    return -1;

	calc_consensus(contigs[i], 1, io_clength(io, cnum), CON_SUM,
		       con, NULL,  NULL, NULL,
		       consensus_cutoff, quality_cutoff,
		       database_info, (void *)io);

	/* nul terminate - only needed to spot WorkShop warnings */
	con[io_clength(io, cnum)] = 0;

	for (rn = io_clnbr(io, cnum); rn; rn = io_rnbr(io, rn)) {
	    UpdateTextOutput();
	    if (cutoff) {
		sc = check_cassembly_single(io, con, cnum, rn, &posn, &len,
					    minlen, winsize, maxdash,
					    maxperc);
	    } else {
		sc = check_uassembly_single(io, con, cnum, rn, &posn, &len,
					    maxperc, winsize);
	    }

	    if (sc > 0) {
		reads[count] = rn;
		score[count] = sc * 100;
		pos[count] = posn;
		length[count] = len;
		conts[count++] = cnum;
	    }
	}

	xfree(con);
    }

    if (-1 == check_assembly_plot(io, reads, conts, score, pos, length, count,
				  cutoff))
	goto error;

    xfree(reads);
    xfree(conts);
    xfree(pos);
    xfree(length);
    xfree(score);

    return 0;

 error:
    if (reads)
	xfree(reads);
    if (conts)
	xfree(conts);
    if (pos)
	xfree(pos);
    if (length)
	xfree(length);
    if (score)
	xfree(score);

    return -1;
}
