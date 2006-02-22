#include <stdio.h>
#include <math.h>

#include "IO.h"
#include "gap_globals.h"
#include "finish.h"
#include "finish_utils.h"
#include "finish_long.h"
#include "finish_walk.h"
#include "finish_reverse.h"
#include "finish_filter.h"
#include "find_fragments.h"
#include "xalloc.h"
#include "misc.h"
#include "tagUtils.h"

static void classify_callback(GapIO *io, int contig, int start, int end,
			      seq_frag *frag, int num_frags,
			      void *clientdata);

static int analyse_templates(finish_t *fin);

/*
 * ---------------------------------------------------------------------------
 * The main finishing code
 * ---------------------------------------------------------------------------
 */

typedef struct {
    con_bits_t *con_bits;	/* Consensus classification table */
    int nbits;			/* Size of con_bits */
    int start;			/* Start position in contig */
    unsigned int *bits;		/* Our main return bit-pattern */
    float *conf;		/* Consensus confidence */
    char *con1;			/* Consensus top strand */
    char *filtered;		/* Low-complex filtered cons */
    vcontig_t *vc;		/* Virtual contig */
    int cvec_left;		/* Is their cosmid at the start? */
    int cvec_right;		/* Is their cosmid at the end? */
    int *template_dup;		/* Linked-list (as array) of dup. templates */
    int *virtual;		/* Number of virtual seqs at this pos? */
} depth_t;


/*
 * Returns 1 if template t1 and t2 have been identified as duplicates of
 *     one another, due to probable cross-contamination.
 * Returns 0 otherwise.
 */
static int template_dupped(int *template_dup, int t1, int t2) {
    int t;

    if (!template_dup)
	return 0;
    if (!template_dup[t1] || !template_dup[t2])
	return 0;

    t = template_dup[t1];
    while (t != t2 && t != t1)
	t = template_dup[t];

    return t == t2;
}


/*
 * classify_callback
 *
 * We are given a set of fragments where each fragment covers a range
 * of consensus such that we have the same number of fragments
 * covering the entire range. (Ie no extra sequences start or end within
 * that consnesus range.)
 *
 * The clientdata argument passed here is a depth_t structure. This contains
 * the main input and output to our parent function (actually the function
 * that called find_fragments).
 *
 * We process these according to the classification bit-pattern stored in
 * depth_t->con_bits.
 *
 * Arguments:
 *	io		gap4 IO handle
 *	contig		contig number
 *	start,end	range of bases in contig that all these fragments span
 *	frag		an array of sequence frags (type seq_frag)
 *	num_frags	size of 'frag' array
 *	clientdata	passed from find_fragments; cast into depth_t struct.
 *
 * Returns:
 *	No return value, but the clientdata (depth_t) is modified.
 */
/*ARGSUSED*/
static void classify_callback(GapIO *io, int contig, int start, int end,
			      seq_frag *frag, int num_frags,
			      void *clientdata) {
    int i, j, k, bit;
    depth_t *d = (depth_t *)clientdata;
    static int *templates = NULL;
    static int ntemplates = 0;
    int template_depth;
    int j_end;
    int do_left = 1;
    int do_right = 1;

    /* --- Sequence depth --- */
    /* Trivial. (== num_frags) */


    /* --- Template depth --- */
    /* (more below) */
    if (num_frags > ntemplates) {
	ntemplates = num_frags;
	templates = (int *)xrealloc(templates, ntemplates * sizeof(int));
    }


    /* --- Perform our strand, confidence and chemistry classifications --- */
    /* These are things which vary from sequence to sequence */
    for (i = 0; i < num_frags; i++) {
	vrseq_t *vr;
	GReadings r;
	GNotes n;
	int virtual, fake, nnote;
	int cons_type = str2type("FAKE");

	vr = vrseq_index2ptr(d->vc, frag[i].num);
	
	if (vr->vseq) {
	    r = vr->vseq->r;
	    virtual = 1;
	} else {
	    gel_read(io, frag[i].num, r);
	    virtual = 0;
	}


	/*
	 * Check if it has a FAKE note, so we can avoid considering it to
	 * be verification of top or bottom strand.
	 */
	for (nnote = r.notes; nnote; nnote = n.next) {
	    note_read(io, nnote, n);
	    if (n.type == cons_type) {
		break;
	    }
	}
	fake = nnote ? 1 : 0;


	/* Also store template numbers here for counting below */
	templates[i] = r.template;

	j_end = end - d->start;
	for (j = start - d->start; j <= j_end; j++) {
	    int cj = d->conf[j];
	    int bitval = d->bits[j];
	    con_bits_t *cbits, *cend = &d->con_bits[d->nbits];
	    int contig_len = io_clength(io, contig);

	    if (d->virtual)
		d->virtual[j] += virtual;

	    for (cbits = d->con_bits; cbits != cend; cbits++) {
		switch(cbits->type) {
		case CLASS_STRAND_TOP:
		    bit = !fake && (r.sense == GAP_SENSE_ORIGINAL);
		    break;
		    
		case CLASS_STRAND_BOTTOM:
		    bit = !fake && (r.sense == GAP_SENSE_REVERSE);
		    break;

		case CLASS_SEQ_DEPTH_GT:
		case CLASS_SEQ_DEPTH_GE:
		    bit = 0; /* Do this later */
		    break;

		case CLASS_TEMP_DEPTH_GT:
		case CLASS_TEMP_DEPTH_GE:
		    bit = 0; /* Do this later */
		    break;

		case CLASS_CONFIDENCE_GT:
		    bit = cj > cbits->arg;
		    break;

		case CLASS_CONFIDENCE_GE:
		    bit = cj >= cbits->arg;
		    break;

		case CLASS_CHEMISTRY:
		    bit = (r.chemistry == cbits->arg);
		    break;

		case CLASS_CONTIG_LEFT_END:
		    bit = do_left
			? (d->start == 1 && j == 0 && !d->cvec_left)
			: 0;
		    /*
		     * If this seq extends the contig, clear any previous
		     * setting of CLASS_CONTIG_LEFT and and disable subsequent
		     * setting.
		     */
		    if (d->start == 1 &&
			j == 0 &&
			r.position < 1) {
			bitval &= ~(1 << cbits->bit);
			do_left = 0;
		    }
		    break;

		case CLASS_CONTIG_RIGHT_END:
		    bit = do_right
			? (d->start - 1 + j == contig_len-1 && !d->cvec_right)
			: 0;
		    /*
		     * If this seq extends the contig, clear any previous
		     * setting of CLASS_CONTIG_RIGHT and and disable
		     * subsequent setting.
		     */
		    if (d->start - 1 + j == contig_len-1 && 
			r.position + r.sequence_length-1 > contig_len) {
			bitval &= ~(1 << cbits->bit);
			bit = 0;
			do_right = 0;
		    }
		    break;

		default:
		    bit = 0;
		    break;
		}

		bitval |= bit << cbits->bit;
	    }

	    d->bits[j] = bitval;
	}
    }


    /* --- More template depth */
    /* Check for duplicate templates; num_frags^2 / 2 ops */
    template_depth = 0;
    for (i = 0; i < num_frags; i++) {
	for (j = i+1; j < num_frags; j++) {
	    if (templates[i] == templates[j] ||
		template_dupped(d->template_dup, templates[i], templates[j]))
		break;
	}
	if (j == num_frags)
	    template_depth++;
    }

    /* --- Now perform our sequence and template depth classifications --- */
    /* These are things which vary from consensus base to consensus base */
    for (i = start; i <= end; i++) {
	for (k = 0; k < d->nbits; k++) {
	    int bit;
	    int i2 = i - d->start;

	    switch (d->con_bits[k].type) {
	    case CLASS_SEQ_DEPTH_GT:
		bit = (num_frags > d->con_bits[k].arg);
		break;

	    case CLASS_SEQ_DEPTH_GE:
		bit = (num_frags >= d->con_bits[k].arg);
		break;

	    case CLASS_TEMP_DEPTH_GT:
		bit = (template_depth > d->con_bits[k].arg);
		break;

	    case CLASS_TEMP_DEPTH_GE:
		bit = (template_depth >= d->con_bits[k].arg);
		break;

	    case CLASS_LOW_COMPLEXITY:
		bit = (d->filtered[i2] == FILTER_LOWCOMPLEX);
		break;
		
	    case CLASS_POLY_A:
		bit = (d->filtered[i2] == FILTER_POLYA);
		break;
		
	    case CLASS_POLY_C:
		bit = (d->filtered[i2] == FILTER_POLYC);
		break;
		
	    case CLASS_POLY_G:
		bit = (d->filtered[i2] == FILTER_POLYG);
		break;
		
	    case CLASS_POLY_T:
		bit = (d->filtered[i2] == FILTER_POLYT);
		break;
		
	    case CLASS_POLY_K:
		bit = (d->filtered[i2] == FILTER_POLYK);
		break;
		
	    case CLASS_POLY_M:
		bit = (d->filtered[i2] == FILTER_POLYM);
		break;
		
	    case CLASS_POLY_R:
		bit = (d->filtered[i2] == FILTER_POLYR);
		break;
		
	    case CLASS_POLY_S:
		bit = (d->filtered[i2] == FILTER_POLYS);
		break;
		
	    case CLASS_POLY_W:
		bit = (d->filtered[i2] == FILTER_POLYW);
		break;
		
	    case CLASS_POLY_Y:
		bit = (d->filtered[i2] == FILTER_POLYY);
		break;
		
	    default:
		bit = 0;
	    }

	    d->bits[i2] |= bit << d->con_bits[k].bit;
	}
    }
}

unsigned int *classify_bases(finish_t *fin,
			     int start,
			     int end,
			     int **virtual,
			     int (*info_func)(int        job,
					      void       *mydata,
					      info_arg_t *theirdata),
			     void *info_data) {
    int len = end - start + 1;
    depth_t d;

    if (start < 1)
	start = 1;
    if (end > io_clength(fin->io, fin->contig))
	end = io_clength(fin->io, fin->contig);

    analyse_templates(fin); /* Does nothing after the first time through */

    /* Allocate */
    d.start = start;
    d.bits = NULL;
    d.conf = fin->qual+start-1;
    d.con1 = fin->cons+start-1;
    d.filtered = fin->filtered+start-1;
    d.con_bits = fin->classify;
    d.nbits = fin->nclassify;
    d.vc = fin->vc;
    d.cvec_left = fin->cvec_left;
    d.cvec_right = fin->cvec_right;
    d.template_dup = fin->template_dup;
    if (virtual) {
	*virtual = d.virtual = (int *)xcalloc(len, sizeof(int));
	if (!d.virtual)
	    return NULL;
    } else {
	d.virtual = NULL;
    }

    if (NULL == (d.bits = (unsigned int *)xcalloc(len, sizeof(*d.bits))))
	return NULL;

    /*
     * Fill out the template and sequence depth arrays, using a more
     * parallel stepping through the sequence assembly.
     */
    find_fragments(fin->io, fin->contig, start, end,
		   info_func, info_data,
		   classify_callback, &d);

    return d.bits;
}

experiments_t *generate_experiments(finish_t *fin, int solution_ind,
				    int *nexp_p, int *cover_end,
				    int *extending) {
    int i;
    experiments_t *exp = NULL, *exp_tmp;
    int nexp = 0, nexp_tmp;
    int solution = fin->solution_bits[solution_ind];
    int pos = fin->start + solution_ind;
    int lowest_cend = INT_MAX;
    int first_try;

    *extending = 0;

    /*
     * Our preferred solution is to walk leftwards from the end
     * of a region requiring the bottom strand and walk rightwards
     * from the start of a region requiring the top strand.
     * However sometimes this simply isn't possible (if it implies
     * placing a primer off the end of the contig for example).
     *
     * So if the first try does not find any solutions then we
     * try from the less optimal end to see if this produces any
     * results instead.
     */
    for (first_try = 1; first_try >= 0; first_try--) {
	for (i = 0; i < 16; i++) {
	    exp_tmp = NULL;
	    
	    if (solution & (1<<i)) {
		int sol_start, sol_end;
		int strand = (solution >> 16) & 0xff;
		int chem   = (solution >> 24) & 0xff;

		/*
		 * Identify the extents over which this solution also extends.
		 * This may be useful for careful placement of experiments.
		 *
		 * Special case for contig start and end - make sure we do not
		 * drift too far from the end; for walking experiments.
		 */
		if (fin->prob_bits[solution_ind] & 3) {
		    /* extend left or right end */
		    sol_start = solution_ind;
		    sol_end   = solution_ind;
		    *extending = 1;
		} else {
		    for (sol_start = solution_ind;
			 sol_start > 0 &&
			     fin->solution_bits[sol_start-1] & (1<<i);
			 sol_start--)
			;
		    for (sol_end = solution_ind;
			 sol_end < fin->length-1 &&
			     fin->solution_bits[sol_end+1] & (1<<i);
		     sol_end++)
			;
		}
		
		sol_start += fin->start;
		sol_end += fin->start;
		
		if (first_try) {
		    if (strand == 2)
			pos = sol_end;
		    else if (strand == 1)
			pos = sol_start;
		    /* else pos stay as left base */
		} else {
		    if (strand == 2 && sol_start != sol_end)
			pos = sol_start;
		    else if (strand == 1 && sol_start != sol_end)
			pos = sol_end;
		    else
			continue;
		}

		/*
		 * If we're near the right end and the end needs extending then
		 * we should also set extending to be true to ensure that we
		 * score extending experiments correctly.
		 */
		if ((fin->prob_bits[fin->end - fin->start] & 3) &&
		    io_clength(fin->io, fin->contig) - sol_end < 500 &&
		    *extending != 1) {
		    *extending = 2;
		} else if ((fin->prob_bits[0] & 3) &&
			   sol_start - 1 < 500 &&
			   *extending != 1) {
		    *extending = 2;
		}

		if (fin->opts.debug[FIN_DEBUG_SCORE])
		    printf("Solution %d strand %d chem %d (prob %d, bits %d) "
			   "at pos %d: covers %d..%d, extending=%d\n",
			   i, strand, chem, fin->prob_bits[solution_ind],
			   fin->base_bits[solution_ind], pos, sol_start,
			   sol_end, *extending);
		if (sol_end < lowest_cend)
		    lowest_cend = sol_end;

		switch (i) {
		case 0:
		    /* Do nothing */
		    exp_tmp = NULL;
		    break;
		case 1:
		    /* Resequence */
		    exp_tmp = experiment_reseq(fin, pos, chem, strand,
					       &nexp_tmp, 0 /* normal */);
		    break;
		case 2:
		    /* Primer walk off sequencing vector */
		    exp_tmp = experiment_walk(fin, pos, chem, strand,
					      sol_start, sol_end, &nexp_tmp,
					      EXPERIMENT_VPWALK);
		    break;
		case 3:
		    /* Long-gel */
		    exp_tmp = experiment_reseq(fin, pos, chem, strand,
					       &nexp_tmp, 1 /* long */);
		    break;
		case 4:
		    /* PCR */
		    fprintf(stderr, "PCR - not implemented yet\n");
		    exp_tmp = NULL;
		    break;
		case 5:
		    /* Primer walk off BAC, YAC, PAC or real chromosome */
		    exp_tmp = experiment_walk(fin, pos, chem, strand,
					      sol_start, sol_end, &nexp_tmp,
					      EXPERIMENT_CPWALK);
		    break;
		case 6:
		    /* Reverse sequence (or fwd if reverse present) */
		    exp_tmp = experiment_reverse(fin, pos, chem, strand,
						 sol_start, sol_end,
						 &nexp_tmp);
		    break;
		default:
		    /* Unknown */
		    printf("Unknown experiment bit %d (0x%x)\n", i, 1<<i);
		    exp_tmp = NULL;
		}
	    }
	
	    if (exp_tmp) {
		int j;
		for (j = 0; j < nexp_tmp; j++) {
		    exp = xrealloc(exp, ++nexp * sizeof(*exp));
		    /* WORKSHOP
		     * 4 byte padding between exp_tmp[j].log_func and .data.
		     */
		    exp[nexp-1] = exp_tmp[j];
		}
		xfree(exp_tmp);

		/*
		 * Whether first or second try; we have potential solutions
		 * so stop looking.
		 */
		break;
	    }
	}
    }

    if (cover_end) {
	*cover_end = (lowest_cend != INT_MAX) ? lowest_cend : 0;
    }

    *nexp_p = nexp;
    return exp;
}

int experiment_score_sort(const void *p1, const void *p2) {
    experiments_t *e1 = (experiments_t *)p1;
    experiments_t *e2 = (experiments_t *)p2;

    /* "e2->score - e1->score" will not work as this is floating point */
    return e2->score > e1->score ? 1 : (e2->score == e1->score ? 0 : -1);
}

static void score_experiments(Tcl_Interp *interp,
			      finish_t *fin,
			      GapIO *io,
			      experiments_t *exp,
			      int nexp,
			      vcontig_t *vc, 
			      int extending) {
    int i, j;
    vrseq_t *vrseq;
    unsigned int *bits, *probs1, *probs2;
    double nprobs1, nprobs2;
    double mand_probs1, mand_probs2;
    int pos, len;
    int extended;
    char *old_cons;
    float *old_qual;
    int new_cvec_left = fin->cvec_left, new_cvec_right = fin->cvec_right;
    int tmp_cvec_left, tmp_cvec_right;
    int shift, nearby_probs;
    int *virtual;
    int virtual_count=0;
    int group_size;
    
    for (i = 0; i < nexp; i++)
	exp[i].score = 666;

    for (i = 0; i < nexp; i++) {
	pos = exp[i].r.position;
	len = exp[i].r.sequence_length;

	if (len <= 0) {
	    exp[i].score = 0;
	    continue;
	}

	/*
	 * Identify how many other experiments are part of this same group.
	 */
	group_size = 0;
	for (j = 0; j < nexp; j++) {
	    if (exp[i].group_id == exp[j].group_id)
		group_size++;
	}

	extended = 0;
	nprobs1 = nprobs2 = mand_probs1 = mand_probs2 = 0;

	/*
	 * Special case. If extending is flagged as true then we are
	 * attempting to extend the contig. In this case we score differently
	 * than when not extending. If this sequence extends the contig
	 * then we decrement the extension amount from nprobs2 (we
	 * treat it as a negative problem).
	 * When extending is false we skip this stage.
	 * We also modify the cvec_left and cvec_right flags to fool the
	 * classification into considering that there is no desire to extend
	 * this contig. This is used (see below) to count the number of
	 * remaining mandatory problems left.
	 *
	 * short extensions are given a greatly reduced bonus, so that we don't
	 * consider experiments that only extend by a short amount as worthy
	 * (unless they also happen to solve other problems). Given
	 * everything else as equal, we still want to prefer experiments that
	 * extend (even if by a short amount) as better than ones that don't.
	 */
	if (pos-1 < fin->left_extent) {
	    extended = fin->left_extent - (pos-1);
	    if (extending)
		if (extended >= fin->opts.min_extension)
		    nprobs2 -= extended * fin->pscore[0];
		else
		    nprobs2 -= extended * fin->pscore[0] / 10;
	    new_cvec_left = 1;
	}
	if (pos+len-1 > fin->right_extent) {
	    extended = (pos+len-1) - fin->right_extent;
	    if (extending)
		if (extended >= fin->opts.min_extension)
		    nprobs2 -= extended * fin->pscore[1];
		else
		    nprobs2 -= extended * fin->pscore[1] / 10;
	    new_cvec_right = 1;
	}

	if (pos < fin->start) {
	    len += pos-fin->start;
	    shift = fin->start-pos;
	    pos = 1;
	} else {
	    shift = 0;
	}

	if (len <= 0) {
	    exp[i].score = 0;
	    continue;
	}

	/*
	 * Generally we already know the classifications, except when
	 * dealing with contig extensions.
	 */
	if (pos + len - 1 <= fin->end) {
	    probs1 = (unsigned int *)xmalloc(len * sizeof(int));
	    if (!probs1)
		return;
	    memcpy(probs1, &fin->prob_bits[pos-fin->start], len * sizeof(int));
	} else {
	    bits = classify_bases(fin, pos, pos + len, NULL,
				  virtual_info_func, vc);

	    probs1 = finishing_rules(interp, fin, pos-fin->start,
				     fin->prob_script, bits, len);
	    xfree(bits);
	}

	/* Add new sequence */
	vrseq = new_vrseq(vc);
	vrseq->vseq->r = exp[i].r;
	link_vrseq(vc, vrseq, exp[i].r.position);

	if (pos + len >= fin->alloc_len) {
	    fin->alloc_len = pos + len + 1;
	    fin->cons = (char *)xrealloc(fin->cons, fin->alloc_len);
	    fin->qual = (float *)xrealloc(fin->qual, fin->alloc_len *
					  sizeof(float));
	    fin->vc->cons = fin->cons;
	}
	old_cons = (char *)xmalloc(len+1);
	old_qual = (float *)xmalloc((len+1)*sizeof(float));
	if (!old_cons || !old_qual) {
	    del_vrseq(vc, vrseq);
	    return;
	}
	/* WORKSHOP
	 * When extending a contig with the realloc above the extension
	 * is not initialised, but for simplicity of code we take a temporary
	 * backup of it anyway. Under some cases WORKSHOP will flag this as
	 * an error.
	 */
	memcpy(old_cons, &fin->cons[pos-1], len+1);
	memcpy(old_qual, &fin->qual[pos-1], (len+1)*sizeof(float));
	if (fin->opts.no_consensus) {
	    /* Fake consensus */
	    int i;
	    memset(&fin->cons[pos-1], 'A', len);
	    for (i = pos-1; i < pos+len; i++) {
		fin->qual[i] = 1.0;
	    }
	} else {
	    calc_consensus(fin->contig, pos, pos + len, CON_SUM,
			   &fin->cons[pos-1], NULL,
			   &fin->qual[pos-1], NULL,
			   gap4_global_get_consensus_cutoff(),
			   gap4_global_get_quality_cutoff(),
			   virtual_info_func, (void *)vc);
	}

	/*
	 * Classify, but not counting the end-extension types if we
	 * have succeeding in extending the contig.
	 */
	tmp_cvec_left = fin->cvec_left;
	tmp_cvec_right = fin->cvec_right;
	
	fin->cvec_left = new_cvec_left;
	fin->cvec_right = new_cvec_right;
	bits = classify_bases(fin, pos, pos + len, &virtual,
			      virtual_info_func, vc);
	fin->cvec_left = tmp_cvec_left;
	fin->cvec_right = tmp_cvec_right;

	probs2 = finishing_rules(interp, fin, pos-fin->start,
				 fin->prob_script, bits, len);
	xfree(bits);

	/* Count problems */
	for (j = 0; j < len; j++) {
	    unsigned int bit, biti;
	    /*
	     * FIXME: Add seq/conf to vseq so that we can compute the virtual
	     * confidence, and hence rate the problems by the chance of their
	     * resolution.
	     */
	    double err;

	    if (!probs1[j] && !probs2[j])
		continue; /* a trivial optimisation */

	    if (vrseq->vseq->conf)
		err = pow(10, vrseq->vseq->conf[j+shift] / -10.0);
	    else
		err = 0;
	    /* err = pow(10, fin->qual[pos-1+j] / -10.0); */

	    for (biti = 0, bit = 1; biti < 32; bit *= 2, biti++) {
		/* Count the number of problems before/after experiment */
		if (probs1[j] & bit) {
		    nprobs1 += fin->pscore[biti];
		    /* Fixed - except for a minor error component */
		    if (!(probs2[j] & bit) && fin->pscore[biti] != 0) {
			nprobs2 += err;
		    }
		}
		if (probs2[j] & bit) {
		    /* Not fixed */
		    nprobs2 += fin->pscore[biti];
		}

		/* Count mandatory problems that haven't been solved */
		if (j+pos >= fin->start && j+pos <= fin->end) {
		    if ((probs1[j] & fin->opts.prob_mandatory) & bit) {
			mand_probs1 += fin->mscore[biti];
		    }
		    if ((probs2[j] & fin->opts.prob_mandatory) & bit) {
			mand_probs2 += fin->mscore[biti];
		    }
		}
	    }
	}

	/* Count the number of bases already covered by a virtual sequence */
	for (virtual_count = j = 0; j < len; j++) {
	    if (virtual[j] > 1)
		virtual_count++;
	}

	xfree(virtual);
	xfree(probs1);
	xfree(probs2);

	if (extending == 1 && !extended) {
	    /*
	     * Regardless if it solved other things, if it didn't extend
	     * we'll reject it.
	     */
	    if (fin->opts.debug[FIN_DEBUG_SCORE])
		printf("%c%d: reject - no extension\n", "+-"[exp[i].r.sense],
		       i);
	    exp[i].score = 0;
	} else if (mand_probs2 && mand_probs1 &&
		   mand_probs2/mand_probs1 > fin->opts.mandatory_ratio) {
	    /*
	     * If we have lots of mandatory problems still left then we'll
	     * reject this experiment anyway.
	     * An example is when we need two be covered by two templates. Even
	     * if this experiment happens to solve lots of problems (eg double
	     * stranded, low confidence, etc), if it doesn't solve our key
	     * requirement of double-templates then we'll still need to do a
	     * further experiment later (which will probably solve the other
	     * problems too).
	     */
	    if (fin->opts.debug[FIN_DEBUG_SCORE])
		printf("%c%d: reject %d - too many mandatory problems not solved."
		       " %f/%f = %f\n",
		       "+-"[exp[i].r.sense], i, exp[i].r.template,
		       mand_probs2, mand_probs1,
		       mand_probs2 / mand_probs1);
	    exp[i].score = 0;

	} else {
	    if (fin->opts.debug[FIN_DEBUG_SCORE])
		printf("%c%d (t:%d, pos:%d..%d): probs %f -> %f / %f => %f "
		       "mprobs %f / %f\n",
		       "+-"[exp[i].r.sense], i, exp[i].r.template,
		       exp[i].r.position,
		       exp[i].r.position + exp[i].r.sequence_length-1,
		       nprobs1, nprobs2,
		       exp[i].cost,
		       (nprobs1-nprobs2) / exp[i].cost,
		       mand_probs1, mand_probs2);
	    exp[i].score = (nprobs1-nprobs2) / exp[i].cost;
	}

	/*
	 * Decrease score based on number of solutions in this group. (Eg
	 * the number of templates chosen for a primer-walk experiment.)
	 * At the requested number of solutions no decrease is made.
	 */
	exp[i].score *= pow(MIN(group_size, exp[i].group_num)/
			    (double)exp[i].group_num, 0.3);

	/* Penalise for over use of a template */
	if (fin->template_used[exp[i].r.template] >=
	    fin->opts.pwalk_use_template) {
	    exp[i].score *= fin->opts.pwalk_use_template_score;
	    if (fin->opts.debug[FIN_DEBUG_SCORE])
		printf("Penalty (overuse). score=%f\n", exp[i].score);
	}

	/*
	 * Minor addition to score for determining new bases on this
	 * template (regardless of whether they solve problems.
	 * This is simply to improve the sorting when everything else is
	 * considered as equal.
	 * Only do if score is already > 0 and it's got a template (ie is
	 * not a chromosomal primer-walk).
	 */
	if (exp[i].score > 0 && exp[i].r.template) {
	    exp[i].score +=
		(len - template_covered(io, fin->tarr[exp[i].r.template],
					fin->contig, pos, pos + len-1))
		*.00005;
	    if (fin->opts.debug[FIN_DEBUG_SCORE])
		printf("Bonus (extra bases). score=%f\n", exp[i].score);
	}

	/*
	 * Check for nearby problems. If there are no mandatory
	 * problems within +/- 100 bases either side of this sequence then
	 * give the score an extra bonus. This encourages use of one
	 * experiment instead of 2.
	 */
	nearby_probs = 0;
	for (j = -100; j < 0; j++) {
	    int bit, biti;
	    int p = j+pos;
	    if (p < fin->start || p > fin->end)
		continue;

	    for (biti = 0, bit = 1; biti < 32; bit *= 2, biti++) {
		if ((fin->prob_bits[p-1] & fin->opts.prob_mandatory) & bit)
		    nearby_probs += fin->mscore[biti];
	    }
	}
	for (j = len; j < len+100; j++) {
	    int bit, biti;
	    int p = j+pos;
	    if (p < fin->start || p > fin->end)
		continue;

	    for (biti = 0, bit = 1; biti < 32; bit *= 2, biti++) {
		if ((fin->prob_bits[p-1] & fin->opts.prob_mandatory) & bit)
		    nearby_probs += fin->mscore[biti];
	    }
	}
	if (nearby_probs == 0) {
	    exp[i].score *= 1.5;
	    printf("Bonus due to solving all nearby problems. "
		   "Score=%f\n", exp[i].score);
	}

	/* Penalise for sequences overlapping over virtual sequences */
	if (virtual_count && fin->opts.penalise_overlaps) {
	    exp[i].score *= MIN(1, 3.0*(1.0-(double)virtual_count/len));
	    printf("Penalty from other virtual seqs (%d of %d bases)"
		   " => Score = %f\n",
		   virtual_count, len, exp[i].score);
	}

	del_vrseq(vc, vrseq);

#if 0
	calc_consensus(fin->contig, pos, pos + len, CON_SUM,
		       &fin->cons[pos-1], NULL,
		       &fin->qual[pos-1], NULL,
		       gap4_global_get_consensus_cutoff(),
		       gap4_global_get_quality_cutoff(),
		       virtual_info_func, (void *)vc);
#endif
	memcpy(&fin->cons[pos-1], old_cons, len+1);
	memcpy(&fin->qual[pos-1], old_qual, (len+1)*sizeof(float));
	xfree(old_cons);
	xfree(old_qual);

    }

    qsort(exp, nexp, sizeof(*exp), experiment_score_sort);

    return; /* As we've sorted them, the best is always element 0 */
}

static void add_experiment(Tcl_Interp *interp,
			   finish_t *fin,
			   GapIO *io,
			   int start,
			   int end,
			   unsigned int *base_bits,
			   unsigned int *prob_bits,
			   unsigned int *solution_bits,
			   experiments_t *best_exp,
			   vcontig_t *vc) {
    vrseq_t *vrseq;
    unsigned int *bits, *probs, *soln;
    int i;
    int pos, len;

    for (i = 0; i < 2; i++) {
	FILE *fp = i ? stdout : fin->out_fp;
	experiments_t *e = best_exp;

	if (!fp)
	    continue;

	fprintf(fp, "=== EXPERIMENT %d/%d: pos %d-%d\n",
		e->expt_id, e->group_id,
		e->r.position,
		e->r.position + e->r.sequence_length);
	fprintf(fp, "=== Score: %f\n", e->score);
	fprintf(fp, "=== Chemistry %d\n", e->r.chemistry);
	if (e->r.template)
	    fprintf(fp, "=== Template name %s\n",
		    get_template_name(io, e->r.template));
	fprintf(fp, "=== Template score %f\n", e->t_score);
	fprintf(fp, "=== Template direction %c (%d)\n",
		"?+-"[e->t_dir + 1], e->t_dir);
	if (e->log_func)
	    e->log_func(fp, fin, e, fin->contig, &fin->tag_list);
	fprintf(fp, "\n");

	fflush(fp);
    }

    /* Increment template usage */
    fin->template_used[best_exp->r.template]++;

    /* Add in new sequence */
    vrseq = new_vrseq(vc);
    vrseq->vseq->r = best_exp->r;
    link_vrseq(vc, vrseq, best_exp->r.position);

    /* Compute new consensus and confidence values */
    pos = best_exp->r.position;
    len = best_exp->r.sequence_length;

    if (pos-1 < fin->left_extent)
	fin->left_extent = pos-1;
    if (pos+len-1 > fin->right_extent)
	fin->right_extent = pos+len-1;

    if (pos-1 < 0) {
	len += pos-1;
	pos = 1;
    }

    if (pos + len >= fin->alloc_len) {
	fin->alloc_len = pos + len + 1;
	fin->cons = (char *)xrealloc(fin->cons, fin->alloc_len);
	fin->qual = (float *)xrealloc(fin->qual, fin->alloc_len *
				      sizeof(float));
	fin->vc->cons = fin->cons;
    }

    if (fin->opts.no_consensus) {
	/* Fake consensus */
	int i;
	memset(&fin->cons[pos-1], 'A', len);
	for (i = pos-1; i < pos+len; i++) {
	    fin->qual[i] = 1.0;
	}
    } else {
	calc_consensus(fin->contig, pos, pos + len, CON_SUM,
		       &fin->cons[pos-1], NULL,
		       &fin->qual[pos-1], NULL,
		       gap4_global_get_consensus_cutoff(),
		       gap4_global_get_quality_cutoff(),
		       virtual_info_func, (void *)vc);
    }

    /* Compute new classification with this sequence in place */
    bits = classify_bases(fin, best_exp->r.position,
			  best_exp->r.position + best_exp->r.sequence_length,
			  NULL, virtual_info_func, vc);

    /* Recompute problems */
    probs = finishing_rules(interp, fin, pos-fin->start, fin->prob_script,
			    bits, best_exp->r.sequence_length);

    /* Recompute solutions */
    soln = finishing_solutions(interp, fin->solu_script, bits, probs,
			       best_exp->r.sequence_length);

    for (i = 0; i < len; i++) {
	int p = pos + i;
	if (p >= start && p <= end) {
	    base_bits[p-start] = bits[i];
	    prob_bits[p-start] = probs[i];
	    solution_bits[p-start] = soln[i];
	}
    }

    xfree(bits);
    xfree(probs);
    xfree(soln);
}


static int qsort_template_c(const void *p1, const void *p2) {
    template_c **t1 = (template_c **)p1;
    template_c **t2 = (template_c **)p2;
    return (*t1)->start - (*t2)->start;
}

/*
 * Searches for possible cross-contamination intemplates by checking the
 * start and end positions. If these match close enough (currently +/-1 bp)
 * then we consider the templates to be the same thing and we group them
 * together.
 *
 * It fills out the fin->template_dup array. Normally this contains zeros, but
 * where several templates are considered to be duplicates they contain the
 * template numbers of the other templates in a non-bounded circular linked-
 * list (so that 3 or more templates can be grouped together if required).
 *
 * Returns 0 for success
 *        -1 for failure
 */
static int identify_contaminant_templates(finish_t *fin) {
    template_c **temps;
    int i, j, ntemps;
    int tol = fin->opts.find_dup_templates;
    
    /* Option not enabled? return now if so. */
    if (!tol)
	return 0;

    if (fin->opts.debug[FIN_DEBUG])
	printf("Identify contaminant templates...\n");

    if (fin->template_dup)
	xfree(fin->template_dup);
    
    fin->template_dup = (int *)xcalloc(Ntemplates(fin->io)+1, sizeof(int));
    if (!fin->template_dup)
	return -1;

    temps = (template_c **)xcalloc(Ntemplates(fin->io)+1, sizeof(*temps));
    if (!temps)
	return -1;

    /* Sort an array of template_c pointers into start position order */
    for (i = j = 0; i <= Ntemplates(fin->io); i++) {
	if (fin->tarr[i])
	    temps[j++] = fin->tarr[i];
    }
    ntemps = j;
    qsort(temps, ntemps, sizeof(*temps), qsort_template_c);

    /* Identify template groups */
    for (i = 0; i < ntemps-1; ) {
	int last = i;
	j = i+1;

	while (j < ntemps && temps[j]->start - temps[i]->start < tol) {
	    if (ABS(temps[j]->end - temps[i]->end) < tol) {
		fin->template_dup[temps[last]->num] = temps[j]->num;
		last = j;
	    }
	    j++;
	}
	if (last != i) {
	    fin->template_dup[temps[last]->num] = temps[i]->num;
	}
	i = j;
    }

    for (i = 0; i <= Ntemplates(fin->io); i++) {
	if (fin->template_dup[i]) {
	    if (fin->opts.debug[FIN_DEBUG] > 1)
		printf("  dup[%d]=%d (%s)\n", i, fin->template_dup[i],
		       get_template_name(fin->io, fin->template_dup[i]));
	}
    }

    if (fin->opts.debug[FIN_DEBUG] > 1)
	printf("...Done\n");
    xfree(temps);

    return 0;
}

int analyse_templates(finish_t *fin) {
    int i;

    /* Analyse all existing templates */
    if (!fin->tarr) {
	/*
	 * Pick templates for contigs we are interested in, and then flesh
	 * out this info for all contigs using only 'connected' templates.
	 */
	fin->tarr = init_template_checks(fin->io, 1, &fin->contig, 1);
	if (!fin->tarr)
	    return -1;

	for (i = 0; i <= Ntemplates(fin->io); i++) {
	    if (fin->tarr[i]) {
		fin->tarr[i]->oflags |= TEMP_OFLAG_CVEC;
		if (fin->opts.chk_template_stat == 0)
		    fin->tarr[i]->oflags |= TEMP_OFLAG_IGNORE_PTYPE;
		fin->tarr[i]->min_vector_len = fin->opts.min_vector_len;
		if (!fin->opts.use_avg_insert) {
		    fin->tarr[i]->oflags |= TEMP_OFLAG_MINMAX_SIZE;
		}
	    }
	}

	check_all_templates(fin->io, fin->tarr);


	for (i = 0; i <= Ntemplates(fin->io); i++) {
	    if (fin->tarr[i]) {
		if (fin->tarr[i]->flags & TEMP_FLAG_SPANNING) {
		    get_template_positions(fin->io, fin->tarr[i], fin->contig);
		}
		printf("Template %c%d, span %d, pos=%d-%d, %d-%d %d-%d len %d consist 0x%x flags 0x%x score %f\n",
		       "?+-"[fin->tarr[i]->direction+1],
		       i,
		       (fin->tarr[i]->flags & TEMP_FLAG_SPANNING) ? 1 : 0,
		       fin->tarr[i]->start, fin->tarr[i]->end,
		       fin->tarr[i]->start2, fin->tarr[i]->end2,
		       fin->tarr[i]->min, fin->tarr[i]->max,
		       fin->tarr[i]->computed_length,
		       fin->tarr[i]->consistency,
		       fin->tarr[i]->flags,
		       fin->tarr[i]->score);
	    }
	}

	return identify_contaminant_templates(fin);
    }

    return 0;
}

void implement_solutions(Tcl_Interp *interp, finish_t *fin) {
    int i;
    experiments_t *exp;
    int nexp;
    int found_soln, soln_type;
    int exp_ind;

    printf("In implement_solutions\n");
    fflush(stdout);

    /* Step through contig identifying problems */
    for (i = 0; i < fin->length; i++) {
	int cover_end;
	int extending;
	
	if (!fin->solution_bits[i])
	    continue;

	exp = generate_experiments(fin, i, &nexp, &cover_end, &extending);
	while (Tcl_DoOneEvent(TCL_FILE_EVENTS | TCL_DONT_WAIT));

	if (!exp) {
	    /*
	     * No experiments. Skip to next 'cover_end' point as up to there
	     * we are also unlikely to solve these same problems.
	     */
	    if (cover_end) {
		if (cover_end > i + fin->opts.reseq_length)
		    cover_end = i + fin->opts.reseq_length;

		if (cover_end == fin->length && i < fin->length-1)
		    i = fin->length-2;
		else
		    i = cover_end-1;
	    }
	    continue;
	}

	/*
	 * score_experiments(exp, nexp); # Sorts exp
	 * perform_experiment(exp[0]);
	 *
	 * To score we need to simulate doing the experiment by producing
	 * a new classify_bases() output for the region spanned by that
	 * experiment. Then we need to recompute the finishing_rules over
	 * that region to determine how many problems would be solved.
	 *
	 * Problems: Need to predict new consensus confidence. How?
	 */
	score_experiments(interp, fin, fin->io, exp, nexp, fin->vc,
			  extending);

	/*
	 * We find the best solution, and that will inform us how many
	 * implementations of that solution are required (eg 3 x reseqs).
	 * For each implementation we may have multiple members of a group,
	 * eg for each primer we may chose 4 different templates.
	 *
	 * Hence group_num and nsolutions are independent and not to be
	 * confused.
	 */
	found_soln = 0;
	soln_type = exp[0].type;
	for (exp_ind = 0;
	     exp_ind < nexp && found_soln < exp[0].nsolutions;
	     exp_ind++) {
	    int grp, grpn, e, found;
	    double max_score;

	    /* All nsolutions should be the same type */
	    if (exp[exp_ind].type != soln_type)
		continue;

	    /* And above the min score */
	    if (exp[exp_ind].score < fin->opts.min_score) {
		if (fin->opts.debug[FIN_DEBUG_SCORE])
		    printf("Skipping experiment(s) with score below threshold."
			   " %f < %f\n", exp[exp_ind].score,
			   fin->opts.min_score);
		/* Sorted, so we'll never a better one */
		break;
	    }

	    grpn = 0;
	    grp = exp[exp_ind].group_id;
	    max_score = exp[exp_ind].score;
	    fin->count[exp[exp_ind].type]++;

	    /* Find all the group members for this solution... */
	    found = 0;
	    for (e = 0; e < nexp; e++) {
		if (exp[e].group_id != grp)
		    continue;

		/* ... but only add the best 'group_num' members */
		if (grpn < exp[exp_ind].group_num &&
		    exp[e].score > 0 &&
		    exp[e].score / max_score >= fin->opts.max_score_drop &&
		    exp[e].t_score >= fin->opts.min_template_score) {
		    add_experiment(interp, fin, fin->io,
				   fin->start, fin->end,
				   fin->base_bits, fin->prob_bits,
				   fin->solution_bits,
				   &exp[e], fin->vc);
		    grpn++;
		    found=1;
		}

		/*
		 * Mark this group item as type 0 to prevent picking this
		 * group again when looking for the next solution of the
		 * same type.
		 */
		exp[e].type = 0;
	    }

	    if (found)
		found_soln++;
	}

	if (!found_soln) {
	    if (fin->opts.debug[FIN_DEBUG])
		printf("--- Failed to find solution for position %d..%d\n",
		       i, cover_end ? cover_end : i);

	    if (cover_end) {
		if (cover_end > i + fin->opts.reseq_length)
		    cover_end = i + fin->opts.reseq_length;

		i = cover_end-1;
	    }
	}

	if (exp)
	    xfree(exp);
    }

    Tcl_DStringResult(interp, &fin->tag_list);

    puts("");
    printf("Totals so far:\n");
    printf("    number of long reads:         %d\n",
	   fin->count[EXPERIMENT_LONG]);
    printf("    number of resequences:        %d\n",
	   fin->count[EXPERIMENT_RESEQ]);
    printf("    number of vector walks:       %d\n",
	   fin->count[EXPERIMENT_VPWALK]);
    printf("    number of chromosomal walks:  %d\n",
	   fin->count[EXPERIMENT_CPWALK]);
    printf("Total number of reverse sequences:  %d\n",
	   fin->count[EXPERIMENT_REVERSE]);
    puts("");

    fflush(stdout);

    return;
}

