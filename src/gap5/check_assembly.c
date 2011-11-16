#include "IO.h"
#include "cs-object.h"
#include "gap_globals.h"
#include "misc.h"
#include "qual.h"
#include "align.h"
#include "contig_selector.h"
#include "editor_view.h"
#include "dna_utils.h"
#include "text_output.h"
#include "tcl_utils.h"
#include "tclXkeylist.h"
#include "consen.h"

/*
 * Checks a single reading for correct assembly by analysing the used data.
 *
 * Returns -1 for system error, otherwise a score (0-1000000)
 * 'pos_p' and 'len_p' are filled in with the position and length of the match
 * within the consensus.
 */
int check_uassembly_single(GapIO *io, char *con, int contig, rangec_t *r,
			   float maxperc, int win_len, int ignore_N) {
    int start, end;
    char *seq = NULL, tmp;
    int i, j, mism = 0;
    int worst, worst_pos = -1;
    seq_t *s, *sorig;
    static int lookup[256];
    static int lookup_done = 0;

    if (!lookup_done) {
	for (i = 0; i < 256; i++)
	    lookup[i] = 0;
	lookup['A'] = lookup['a'] = 1;
	lookup['C'] = lookup['c'] = 2;
	lookup['G'] = lookup['g'] = 3;
	lookup['T'] = lookup['t'] = 4;
	lookup['U'] = lookup['u'] = 4;
	lookup['*'] = lookup[','] = lookup['-'] = 5;
	lookup_done = 1;
    }

    /* Get sequence */
    if (!(sorig = s = cache_search(io, GT_Seq, r->rec)))
	return -1;

    /* Complement data on-the-fly */
    if ((s->len < 0) ^ r->comp) {
	s = dup_seq(s);
	complement_seq_t(s);
    }
    seq = s->seq;

    start = s->left;
    end = s->right;

    /* Initialise scoring for divergence checks */
    if (end - start - 1 < win_len) {
	win_len = end - start - 1;
    }
    worst = 0.5 + maxperc * win_len;

    for (i=start-1, j=r->start + start-1; i < start-1 + win_len; i++, j++) {
	if (ignore_N) {
	    if (lookup[seq[i]] && lookup[seq[i]] != lookup[con[j]])
		mism++;
	} else {
	    if (lookup[seq[i]] != lookup[con[j]])
		mism++;
	}
    }

    /* Loop through sequence looking for problems */
    if (ignore_N) {
	do {
	    if (mism >= worst) {
		worst_pos = i;
		worst = mism;
	    }
	    
	    mism -= lookup[seq[i-win_len]] &&
		lookup[seq[i-win_len]] != lookup[con[j-win_len]];
	    i++; j++;

	    if (i < end-1)
		mism += lookup[seq[i]] && lookup[seq[i]] != lookup[con[j]];
	} while (i < end);
    } else {
	do {
	    if (mism >= worst) {
		worst_pos = i;
		worst = mism;
	    }
	    
	    mism -= lookup[seq[i++-win_len]] != lookup[con[j++-win_len]];
	    if (i < end-1)
		mism += lookup[seq[i]] != lookup[con[j]];
	} while (i < end);
    }

    /* Display problem, listing worst score */
    if (worst_pos != -1) {
	//*pos_p = io_relpos(io, rn);
	//*len_p = end - start + 1;

	vmessage("\nReading #%"PRIrec"(%s) has a local percentage "
		 "mismatch of %2.1f\n",
		 s->rec, s->name, 100 * (float)worst / win_len);
	vmessage("SEQ: %.*s\n", end-start+1, &seq[start-1]);
	vmessage("CON: %.*s\n", end-start+1, &con[r->start + start-1]);

	if (sorig != s)
	    xfree(s);

	return 10000 * (float)worst / win_len;
    }

    if (sorig != s)
	xfree(s);

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
    cs = result_data(ca->io, cs_id);

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
	    vmessage("check_assembly match:\n");
	    vmessage("    From contig %s(#%"PRIrec") at %d\n",
		     get_contig_name(ca->io, ABS(obj->c1)),
		     io_clnbr(ca->io, ABS(obj->c1)), obj->pos1);
	    vmessage("    With contig %s(#%"PRIrec") at %d\n",
		     get_contig_name(ca->io, ABS(obj->c2)),
		     io_clnbr(ca->io, ABS(obj->c2)), obj->pos2);
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
	    tg_rec cnum, llino;
	    int pos, id;

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

//	    if ((id = editor_available(cnum, 1)) != -1) {
//		move_editor(id, llino, pos);
//	    } else {
		edit_contig(ca->io, cnum, llino, pos);
//	    }
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
		"check_assembly: #%"PRIrec"@%d len %d, mis %2.2f%%",
		obj->read, obj->pos1, obj->length, ((float)obj->score)/10000);
	return buf;
    }

    return NULL;
}


static int sort_func(const void *p1, const void *p2) {
    obj_checkass *m1 = (obj_match *)p1, *m2 = (obj_checkass *)p2;
    return m2->score - m1->score;
}

void check_assembly_callback(GapIO *io, tg_rec contig, void *fdata,
			     reg_data *jdata) {
    mobj_checkass *r = (mobj_checkass *)fdata;
    obj_cs *cs;
    int cs_id;

    cs_id = type_to_result(io, REG_TYPE_CONTIGSEL, 0);
    cs = result_data(io, cs_id);

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
	    csmatch_remove(io, cs->window, (mobj_repeat *)r,
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

	csmatch_remove(io, cs->window, (mobj_repeat *)r, csplot_hash);
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

int check_assembly_plot(GapIO *io, tg_rec *reads, tg_rec *conts, int *score,
			int *pos, int *length, int count) {
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
    ca->cutoffs = 0;
    strcpy(ca->tagname, CPtr2Tcl(ca));

    val = get_default_string(GetInterp(), gap5_defs, "CHECK_ASSEMBLY.COLOUR");
    strcpy(ca->colour, val);
    ca->linewidth = get_default_int(GetInterp(), gap5_defs,
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
    contig_register(io, 0, check_assembly_callback, (void *)ca, id,
		    REG_REQUIRED | REG_DATA_CHANGE | REG_OPS |
		    REG_NUMBER_CHANGE | REG_ORDER, REG_TYPE_CHECKASS);
    update_results(io);

    return 0;
}

/*
 * Scans through one or more contigs checking each reading for correct
 * assembly. This is simply a check for misaligned data, not looking into
 * cutoff data. (The gap4 method did this, but it hasn't yet been implemented
 * in gap5).
 *
 * Returns -1 for failure, 0 for success.
 */
int check_assembly(GapIO *io, int num_contigs, contig_list_t *contigs,
		   int winsize, float maxperc, int ignore_N) {
    int i, rn, sc, count = 0, allocated = 0;
    char *con;
    tg_rec *reads = NULL, *conts = NULL;
    int *score = NULL, *length = NULL, *pos = NULL;

    for (i = 0; i < num_contigs; i++) {
	tg_rec crec = contigs[i].contig;
	contig_iterator *ci;
	rangec_t *r;
	int start = contigs[i].start, end = contigs[i].end;

	if (NULL == (con = (char *)xmalloc(end-start+1)))
	    return -1;

	calculate_consensus_simple(io, crec, start, end, con, NULL);

	ci = contig_iter_new(io, crec, 0, CITER_FIRST, start, end);
	while (r = contig_iter_next(io, ci)) {
	    UpdateTextOutput();
	    sc = check_uassembly_single(io, con - start, crec, r,
					maxperc, winsize, ignore_N);
	    if (count >= allocated) {
		allocated = allocated ? allocated * 2 : 256;
		reads  = xrealloc(reads, allocated * sizeof(*reads));
		conts  = xrealloc(conts, allocated * sizeof(*conts));
		score  = xrealloc(score, allocated * sizeof(*score));
		length = xrealloc(length, allocated * sizeof(*length));
		pos    = xrealloc(pos, allocated * sizeof(*pos));
		if (!reads || !conts || !score || !length || !pos)
		    goto error;
	    }

	    if (sc > 0) {
		reads[count]   = r->rec;
		score[count]   = sc * 100;
		pos[count]     = r->start;
		length[count]  = r->end - r->start+1;
		conts[count++] = crec;
	    }
	}

	contig_iter_del(ci);
	xfree(con);
    }

    if (-1 == check_assembly_plot(io, reads, conts, score, pos, length, count))
	goto error;

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
