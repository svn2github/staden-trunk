#include <stdio.h>
#include <math.h>

#include "IO.h"
#include "finish_walk.h"
#include "finish_utils.h"
#include "finish_filter.h"
#include "gap_globals.h"
#include "xalloc.h"
#include "misc.h"
#include "search_utils.h"
#include "dna_utils.h"
#include "primlib.h"

#ifdef USE_OSP
#    define SANGER_HACK
#    define PNAME "OSP"
#else
#    define PNAME "Primer3"
#endif


/* Compile time options to tweak */
#define MAX_PRIMER_SEQ 1024

/*
 * FIXME:
 * Additional cost for short reads is divided into this. Larger numbers
 * implies reduced cost of short reads.
 */
#define LENGTH_SCALE 5.0

/*
 * Returns the probability of a sequence from 'start' to 'end' containing
 * an error somewhere within it.
 */
static double primer_err(int start, int end, float *unpadded_qual,
			   finish_t *fin, int *mapping) {
    int j;
    int last_map_pos = 0;
    double success = 1.0;

    for (j = start; j <= end; j++) {
	/* Confidence for this base */
	success *= 1-pow(10, unpadded_qual[j]/-10.0);

	/* Confidence for any pads between bases */
	if (last_map_pos) {
	    int cur_pos = mapping[j];
	    if (cur_pos > last_map_pos + 1) {
		for (last_map_pos++; last_map_pos < cur_pos; last_map_pos++) {
		    success *= 1-pow(10, fin->orig_qual[last_map_pos]/-10.0);
		}
	    } else if (cur_pos < last_map_pos - 1) {
		for (last_map_pos--; last_map_pos > cur_pos; last_map_pos--) {
		    success *= 1-pow(10, fin->orig_qual[last_map_pos]/-10.0);
		}
	    }
	}

	last_map_pos = mapping[j];
    }


    return 1.0 - success;
}


/*
 * near_low_complexity()
 *
 * Checks whether a given position in a given direction is approaching a
 * low complexity region.
 *
 * Arguments:
 *	fin		Finish_t structure
 *	pos		Position
 *	dir		Direction (1 = fwd, 2 = rev)
 *
 * Returns:
 *	1 if low complexity near
 *	0 otherwise
 */
static int near_low_complexity(finish_t *fin, int pos, int dir) {
    int i;
    int i_end;
    int count_filt = 0, count_total = 0;
    
    if (!fin->filtered)
	return 0;

    /*
     * Count how many filtered and non-filtered symbols there are between
     * where the pos and additional 100 basepairs beyond where we
     * expect the sequence to start. If it's got too many then reject.
     */
    if (dir == 1) {
	i_end = MIN(io_clength(fin->io, fin->contig),
		    pos + fin->opts.pwalk_seq_gap+100);
	for (i = pos; i < i_end; i++) {
	    if (is_filtered(fin->filtered[i]))
		count_filt++;
	}
	count_total = i_end - pos;
    } else {
	i_end = MAX(-1, pos - fin->opts.pwalk_seq_gap-100);
	for (i = pos; i > i_end; i--) {
	    if (is_filtered(fin->filtered[i]))
		count_filt++;
	}
	count_total = pos - i_end;
    }

    if (!count_total)
	return 1;

    /*
     * Reject if over 80% is filtered or there are not at least
     * pwalk_seq_gap unique bases.
     */
    if (count_filt / (double)count_total > 0.8 ||
	count_total - count_filt < fin->opts.pwalk_seq_gap)
	return 1;

    return 0;
}

/*
 * Checks whether there is a run of identical bases.
 * Returns the maximum homopolymer length.
 */
static int base_run(char *primer) {
    int count = 1;
    int type = 0;
    int run = 1;
    while (*primer) {
	if (*primer == type) {
	    if (++count > run)
		run = count;
	} else {
	    count = 1;
	    type = *primer;
	}
	primer++;
    }

    return run;
}

/*
 * find_primers()
 *
 * Identifies potential sequencing primers in a particular region of
 * consenus sequence.
 *
 * Arguments:
 *	fin		Finish_t structure
 *	cons_len	Length of cons and qual arrays
 *	start, end	First and last base in cons/qual to search within
 *	dir		Direction of primer (1 == fwd, 2 == rev)
 *	nprimers	Returned: number of suitable primers found
 *
 * Returns:
 *	Success: a malloced array of experiment_walk_t structures
 *	         containing the primer details. "nprimers" contains the
 *	         dimension of this array.
 *	Failure: NULL
 */
#define MAX_NPRIMERS 30
static experiment_walk_t *find_primers(finish_t *fin,
				       int cons_len,
				       int start,
				       int end,
				       int dir /* 1 or 2 */,
				       int *nprimers) {
    int i, j, np;
    experiment_walk_t *prim;
    char unpadded_cons[MAX_PRIMER_SEQ];
    float unpadded_qual[MAX_PRIMER_SEQ];
    int mapping[MAX_PRIMER_SEQ];
    int min_qual;
    double max_err;
    primlib_state *pstate;

    /* Allocate memory */
    prim = (experiment_walk_t *)xmalloc(MAX_NPRIMERS * sizeof(*prim));
    if (!prim)
	return NULL;

    /* Determine min_qual and max_err; depends on proximity to contig ends */
    if (end <= fin->opts.pwalk_end_dist ||
	cons_len - start <= fin->opts.pwalk_end_dist) {
	/* End */
	min_qual = fin->opts.pwalk_min_qual2;
	max_err  = fin->opts.pwalk_max_err2;
    } else {
	/* Middle */
	min_qual = fin->opts.pwalk_min_qual;
	max_err  = fin->opts.pwalk_max_err;
    }

    /* Produce an unpadded copy of the sequence and quality */
    for (i = start, j = 0; i <= end; i++) {
	if (fin->cons[i] != '*') {
	    if (fin->orig_qual[i] < min_qual ||
		fin->orig_prob_bits[i] & fin->opts.pwalk_prob_mask)
		unpadded_cons[j] = '-';
	    else
		unpadded_cons[j] = is_filtered(fin->filtered[i])
		    ? '-'
		    : fin->cons[i];
	    unpadded_qual[j] = fin->orig_qual[i];
	    mapping[j] = i;
	    j++;
	}
    }
    unpadded_cons[j] = 0;


    /* Complement unpadded_* and mapping arrays, if needed */
    if (dir == 2)
	complement_seq_qual_mapping(j /* length */, unpadded_cons,
				    unpadded_qual, mapping);


    /* Initialise Primer Design package */
    pstate = primlib_create();
    if (NULL == pstate)
	return NULL;
    primlib_set_args(pstate, &fin->opts.primer_args);
    
    /* Invoke it */
    printf("Seq=%s\n", unpadded_cons);
    if (primlib_choose(pstate, unpadded_cons) == -1) {
	primlib_destroy(pstate);
	return NULL;
    }

    if (pstate->nprimers > 0)
	if (fin->opts.debug[EXPERIMENT_VPWALK])
	    printf("Best primer score = %f\n", pstate->primers[0].quality);

    /* Parse OSP results and remap to padded positions */
    for (i = np = 0;
	 i < pstate->nprimers && np < MAX_NPRIMERS;
	 i++)
    {
	int len = pstate->primers[i].length;

	/* Check suitability of OSP score */
#ifdef USE_OSP
#ifdef SANGER_HACK
	if (pstate->primers[i].quality < 16)
	    continue;
#endif

	if (pstate->primers[i].quality > fin->opts.pwalk_osp_score)
	    break;
#endif

	/* Create prim[] structure */
	prim[np].osp_score = pstate->primers[i].quality;
	prim[np].start = dir == 1
	    ? mapping[pstate->primers[i].start]
	    : mapping[pstate->primers[i].start + pstate->primers[i].length -1];
	prim[np].end   = dir == 1
	    ? mapping[pstate->primers[i].start + pstate->primers[i].length -1]
	    : mapping[pstate->primers[i].start];
	strncpy(prim[np].primer, &unpadded_cons[pstate->primers[i].start],
		len);
	prim[np].primer[len] = 0;
	prim[np].dir = dir;
	prim[np].secbind = 0;

	/* Does it contain a run of bases? */
	prim[np].homopolymer = base_run(prim[np].primer);
	if (prim[np].homopolymer >= 6) {
	    if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
		printf("Reject "PNAME" primer due to homopolymer >= 6\n");
	    continue;
	}

	/* Is the oligo too close to a low-complexity region? */
	if (near_low_complexity(fin,
				dir == 1 ? prim[np].end : prim[np].start,
				dir)) {
	    if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
		printf("Reject "PNAME" primer %s due to proximity to "
		       "low-complexity sequence\n",
		       prim[np].primer);
	    continue;
	}
	
	/* Compute chance of a sequencing error */
	prim[np].p_err = primer_err(pstate->primers[i].start,
				    pstate->primers[i].start +
				    pstate->primers[i].length - 1,
				    unpadded_qual,
				    fin,
				    mapping);
	if (prim[np].p_err > max_err) {
	    if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
		printf("Reject "PNAME" primer due to probability of "
		       "sequencing error (%f > %f\n",
		       prim[np].p_err, max_err);
	    continue;
	}
#if 0
	if (prim[np].p_err < 0) {
	    if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
		printf("Reject "PNAME" primer due to being in a "
		       "mandatory-resolve region\n");
	    continue;
	}
#endif


	/* Search vector for matches */
	prim[np].secbind = secondary_primer_match(fin, 0/*contig*/,
						  0, 0, 0, 0, 1/*ext*/,
						  prim[np].primer);
	if (prim[np].secbind >= fin->opts.pwalk_max_match) {
	    if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
		printf("Reject "PNAME" primer due to external match "
		       "elsewhere\n");
	    continue;
	}


	/*
	 * FIXME: need to combine scores to take into account percentage
	 * match elsewhere, consensus quality (chance of oligo being
	 * incorrect), and proximity to desired position.
	 */
	/*prim[np].score = prim[np].osp_score * (1+prim[np].p_err); */
#ifdef USE_OSP
	prim[np].score = prim[np].osp_score + 1000*prim[np].p_err; 
#else
	prim[np].score = prim[np].osp_score + 10*prim[np].p_err; 
#endif
	if (fin->opts.debug[EXPERIMENT_VPWALK])
	    printf(PNAME" score=%g p_err=%g => %g %s\n",
		   prim[np].osp_score, prim[np].p_err, prim[np].score,
		   prim[np].primer);
	np++;
    }
    *nprimers = np;
    if (fin->opts.debug[EXPERIMENT_VPWALK])
	printf("Nprimers=%d of which %d are valid\n", pstate->nprimers, np);


    /* Tidy up */
    if (!*nprimers) {
	xfree(prim);
	prim = NULL;
    }
    primlib_destroy(pstate);

    return prim;
}

/*
 * log_walk
 *
 * Outputs experiment information pertinent only to primer-walks.
 *
 * Arguments:
 *	fp		Where to fprintf the information.
 *	io		An open GapIO pointer.
 *	e		The experiment structure to print.
 *	contig		Current contig number.
 */
static void log_walk(FILE *fp, finish_t *fin, experiments_t *e, int contig,
		     Tcl_DString *tagds) {
    static FILE *tagfp = NULL;
    char buf[1024];
    GapIO *io = fin->io;

    if (!tagfp) {
	tagfp = fopen("tags", "w+");
    }

    if (e->type == EXPERIMENT_VPWALK)
	fprintf(fp, "=== Type Vector Primer walk\n");
    else
	fprintf(fp, "=== Type Chromosomal Primer walk\n");
    fprintf(fp, "=== Primer sequence %s\n", e->data.e_walk.primer);
    fprintf(fp, "=== Primer position %d..%d\n",
	    e->data.e_walk.start+1, e->data.e_walk.end+1);
    fprintf(fp, "=== Primer 2ndary=%f\n", e->data.e_walk.secbind);
    fprintf(fp, "=== Primer score=%f, direction '%c'\n",
	    e->data.e_walk.score, "?+-"[e->data.e_walk.dir]);

    /* Add tag on to Tcl dstring */
    if (tagds) {

	if (!fin->opts.pwalk_tag_cons) {
	    int st = e->data.e_walk.start+1;
	    int en = e->data.e_walk.end+1;
	    int rnum;
	    GReadings r;

	    if ((rnum = tag_template(io, contig, e->r.template, &st, &en))) {
		gel_read(io, rnum, r);
		if (r.sense == 0) {
		    st = st - r.position+1 + r.start;
		    en = en - r.position+1 + r.start;
		} else {
		    int len = en - st;
		    en = r.length - (st - r.position+1 + r.start)+1;
		    st = en - len;
		}
		sprintf(buf,
			"%d\n"
			"%s %c %d..%d\n"
			"Sequence %s\n"
			"Template %s",
			rnum,
			fin->opts.pwalk_tag_type
			  ? fin->opts.pwalk_tag_type
			  : "PRIM",
			"?+-"[e->data.e_walk.dir], st, en,
			e->data.e_walk.primer,
			e->type == EXPERIMENT_VPWALK
			? get_template_name(io, e->r.template)
			: "chromosomal/[BPY]AC");
	    } else {
		*buf = 0;
	    }
	}

	if (fin->opts.pwalk_tag_cons || *buf == 0) {
	    sprintf(buf,
		    "-%d\n"
		    "PRIM %c %d..%d\n"
		    "Sequence %s\n"
		    "Template %s",
		    contig, "?+-"[e->data.e_walk.dir],
		    e->data.e_walk.start+1, e->data.e_walk.end+1,
		    e->data.e_walk.primer,
		    e->type == EXPERIMENT_VPWALK
		    ? get_template_name(io, e->r.template)
		    : "chromosomal/[BPY]AC");
	}
	Tcl_DStringAppendElement(tagds, buf);
    }
}

/*
 * generate_chr_exp
 *
 * Having found some primers, we produce experiments for priming directly
 * of the [BPY]AC or chromosome.
 */
experiments_t *generate_chr_exp(finish_t *fin, experiment_walk_t *prim,
				int nprimers, int primer_dir,
				experiments_t *exp, int *nexp_p, int chem) {
    int p;
    int nexp = *nexp_p;

    for (p = 0; p < nprimers && p < fin->opts.pwalk_noligos; p++) {
	int primer_pos_start, primer_pos_end;
	int group_id;
	int s_start;
	double pscore;

	primer_pos_start = prim[p].start+1;
	primer_pos_end = prim[p].end+1;
	group_id = finish_next_group_id(0);

	pscore = secondary_primer_match(fin, -1 /* all */,
					0, 0, 0, 0, 0, prim[p].primer);
	if (pscore > prim[p].secbind)
	    prim[p].secbind = pscore;
	if (pscore >= fin->opts.pwalk_max_match) {
	    if (fin->opts.debug[EXPERIMENT_CPWALK] > 1)
		printf("Reject "PNAME" primer due to consensus match "
		       "elsewhere\n");
	    continue;
	}

	/*
	 * Compute sequence start and end points.
	 *
	 * For now we ignore the chance that this is right at the end of
	 * the BAC (etc). Ideally we need to adjust score for this.
	 */
	s_start = primer_dir == 1 ?
	    primer_pos_end   + fin->opts.pwalk_seq_gap :
	    primer_pos_start - fin->opts.pwalk_seq_gap -fin->opts.pwalk_length;


	/* Create experiment */
	exp = xrealloc(exp, ++nexp * sizeof(*exp));
	memset(&exp[nexp-1].r, 0, sizeof(GReadings));

	/* Add the read-specific bits. */
	exp[nexp-1].r.position = s_start;
	exp[nexp-1].r.sequence_length = fin->opts.pwalk_length;
	exp[nexp-1].r.start = 1;
	exp[nexp-1].r.end = 1 + exp[nexp-1].r.sequence_length + 1;
	exp[nexp-1].r.strand = 0;
	exp[nexp-1].r.sense = primer_dir == 1 ? 0 : 1;
	exp[nexp-1].r.primer = primer_dir == 1 ? 3 : 4; /* custom for/rev */
	exp[nexp-1].r.template = 0;
	exp[nexp-1].r.chemistry = chem;

	/* Generic component */
	exp[nexp-1].type = EXPERIMENT_CPWALK;
	exp[nexp-1].nsolutions = fin->opts.pwalk_nsolutions;
	exp[nexp-1].score = 0;
	exp[nexp-1].cost = fin->cost[EXPERIMENT_CPWALK];
	exp[nexp-1].expt_id = finish_next_expt_id(0);
	exp[nexp-1].group_id = group_id;
	exp[nexp-1].group_num = fin->opts.pwalk_ntemplates;
	exp[nexp-1].t_score = 1; /* max score */
	exp[nexp-1].t_dir = -1; /* unknown 'template direction' */
	exp[nexp-1].log_func = log_walk;

	/* e_walk component */
	exp[nexp-1].data.e_walk = prim[p];

	if (fin->opts.debug[EXPERIMENT_CPWALK])
	    printf("chromsomal read %d: %d-%d (primer at %d)\n",
		   nexp-1, 
		   exp[nexp-1].r.position,
		   exp[nexp-1].r.position +
		   exp[nexp-1].r.sequence_length-1,
		   primer_pos_start);
    }

    *nexp_p = nexp;
    return exp;
}

/*
 * find_templates
 *
 * Having found some primers, we iterate around all the primers checking for
 * suitable templates.
 */
experiments_t *find_templates(finish_t *fin, experiment_walk_t *prim,
			      int nprimers, int primer_dir,
			      experiments_t *exp, int *nexp_p,
			      int prob_start, int prob_end,
			      int prob_end_orig, int chem) {
    int i, j;
    int nexp = *nexp_p;
    int templates_picked[20];
    double pscore;
    int noligos;

    if (fin->opts.pwalk_ntemplates > 10) {
	fin->opts.pwalk_ntemplates = 10;
	fprintf(stderr, "Restricting pwalk_ntemplates parameter to 10\n");
    }

    /*
     * Keep track of how many primers have valid templates. We only need to
     * process up to fin->opts.pwalk_noligos valid ones.
     */
    noligos = 0;

    /* Find suitable templates covering position primer_pos_* */
    for (j = 0; j < nprimers && noligos < fin->opts.pwalk_noligos; j++) {
	int tnum;
	int group_id;
	int round;
	int primer_pos_start, primer_pos_end;
	int t_start, t_end;
	int t_start1, t_end1;
	int t_start2, t_end2;
	int s_start, s_end;
	item_t *clist;

	primer_pos_start = prim[j].start+1;
	primer_pos_end = prim[j].end+1;
	
	if (fin->opts.debug[EXPERIMENT_VPWALK])
	if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
	    printf("primer %d, score = %f, at pos %d-%d 2ndary %f %s\n",
		   j, prim[j].score, primer_pos_start, primer_pos_end,
		   prim[j].secbind, prim[j].primer);

	group_id = finish_next_group_id(0);

	    /*
	     * Two rounds - first consistent templates and then the
	     * inconsistent ones (penalised by increasing their cost) if we
	     * do not have enough consistent ones.
	     */
	tnum = 0;
	memset(templates_picked, 0, fin->opts.pwalk_ntemplates * sizeof(int));
	for (round = 0; round < 2; round++) {
	    if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
		printf("Round %d: %sconsistent templates\n",
		       round+1, round ? "in" : "");
	    for (i = 1;
		 i <= Ntemplates(fin->io) &&
		     tnum < 2*fin->opts.pwalk_ntemplates;
		 i++) {
		char *tname;
		GTemplates t;
		GReadings r;
		int temp_fwd;
		double cost_mult;

		if (!fin->tarr[i])
		    continue;

		if (fin->template_skip && fin->template_skip[i]) {
		    if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
			printf("Template %d(%s) skipped/!available\n", i,
			       get_template_name(fin->io, i));
		    continue;
		}

		if (0 == round &&
		    (fin->tarr[i]->consistency ||
		     fin->template_used[i] >= fin->opts.pwalk_use_template))
		    continue;

		if (1 == round &&
		    !(fin->tarr[i]->consistency ||
		      fin->template_used[i] >= fin->opts.pwalk_use_template))
		    continue;

		/* Increase cost as template score diminishes */
		cost_mult = 1/fin->tarr[i]->score;

#ifdef USE_OSP
		/* Marginally increase cost for higher scoring OSP results */
		cost_mult *= 1+prim[j].score / 1000.0;
#else
		cost_mult *= 1+prim[j].score / 20;
#endif

		/* Increase cost for long homopolymers */
		if (prim[j].homopolymer > 4) {
		    cost_mult *= prim[j].homopolymer-3;
		}

		/*
		 * Is the template in the correct position?
		 * tarr[i]->start/end assume using the minimum insert size
		 * tarr[i]->start2/end2 assumed the max insert size.
		 */
		if (!fin->tarr[i]->consistency) {
		    t_start1 = MIN(fin->tarr[i]->start, fin->tarr[i]->end);
		    t_end1 = MAX(fin->tarr[i]->start, fin->tarr[i]->end);
		    t_start2 = MIN(fin->tarr[i]->start2, fin->tarr[i]->end2);
		    t_end2 = MAX(fin->tarr[i]->start2, fin->tarr[i]->end2);

		    t_start = MIN(t_start1, t_start2);
		    t_end   = MAX(t_end1,   t_end2);
		} else {
		    /*
		     * Start/end are unreliable for inconsistent templates.
		     * min/max isn't correct, but should catch problems it
		     * causes later with our insistence that primers have to
		     * be covered by a read on that template when using
		     * inconsistent templates.
		     */
		    t_start1 = t_start2 = t_start = fin->tarr[i]->min;
		    t_end1   = t_end2   = t_end   = fin->tarr[i]->max;
		}

		/* Does the template cover the problem region? */
		if (t_start > prob_end || t_end < prob_start)
		    continue;

		/* Does the primer exist on this template (max length) */
		if (t_start > primer_pos_end || t_end < primer_pos_start)
		    continue;

		/* What about in min length? If not, increase cost */
		if (t_start1 > primer_pos_end || t_end1 < primer_pos_start)
		    {
			cost_mult /= 0.000001 +
			    template_exists_chance(primer_pos_start,
						   primer_pos_end,
						   t_start1, t_end1,
						   t_start2, t_end2);
		    }

		/* Is it a contaminant template? */
		if (fin->template_dup && fin->template_dup[i]) {
		    if (template_is_dup(fin, templates_picked, tnum, i)) {
			cost_mult *= fin->opts.pwalk_dup_template_cost;
		    }
		}

		/*
		 * For inconsistent templates we make sure that at least
		 * one sequence covers our primer position, otherwise
		 * we're simply pushing our luck TOO far!
		 */
		if (1 == round &&
		    !template_covered(fin->io, fin->tarr[i], fin->contig,
				      primer_pos_start, primer_pos_end)) {
		    if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
			printf("Skipping template %d as it is not covered "
			       "by a sequence\n", fin->tarr[i]->num);
		    continue;
		}
							 
		temp_fwd = (fin->tarr[i]->start >= fin->tarr[i]->end);
	    
		if (fin->opts.debug[EXPERIMENT_VPWALK] > 1) {
		    template_read(fin->io, i, t);
		    tname = TextAllocRead(fin->io, t.name);
		    printf("Pos %c%d: Template %c %d/%s (%d-%d). Reads:",
			   primer_dir == 1 ? '+' : '-',
			   primer_pos_start, "-+"[temp_fwd],
			   i, tname, t_start, t_end);
		    xfree(tname);
		    
		    for (clist = head(fin->tarr[i]->gel_cont); clist;
			 clist = clist->next) {
			gel_cont_t *gc = (gel_cont_t *)(clist->data);
			printf("%d ", gc->read);
		    }
		    putchar('\n');
		}
	    
		clist = head(fin->tarr[i]->gel_cont);
		gel_read(fin->io, ((gel_cont_t *)(clist->data))->read, r);

		/*
		 * Ok, so the template looks like a valid proposition, but
		 * does it have false matches? For inconsistencies due to
		 * large sizes, this may falsely reject the template. (However
		 * poor data _will_ give poor results.)
		 */
		{
		    int st = fin->tarr[i]->start2;
		    int en = fin->tarr[i]->end2;
		    if (st > en) {
			int tmp;
			tmp = st;
			st = en;
			en = tmp;
		    }
		    pscore = secondary_primer_match(fin, fin->contig, st, en,
						    1,
						    primer_dir == 1 ? 0 : 1, 0,
						    prim[j].primer);
		    pscore = MAX(pscore, prim[j].secbind);
		    if (pscore >= fin->opts.pwalk_max_match) {
			if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
			    printf("Reject "PNAME" primer due to secondary "
				   "match on this template (%f >= %f)\n",
				   pscore, fin->opts.pwalk_max_match);
			continue;
		    }

		    /* Increase cost to account for secondary binding */
		    cost_mult *= 1+pscore / 100;
		}

		/*
		 * Compute how much of this reading is valid in this
		 * template. We have two start-end ranges. One for the
		 * minimum sized template and one for the maximum.
		 * We use the maximum, but score based on the difference
		 * of maximum to minimum so that we'd rather pick readings
		 * on templates entirely fitting in the minimum size.
		 */
		s_start = primer_dir == 1
		    ? primer_pos_end   + fin->opts.pwalk_seq_gap
		    : primer_pos_start - fin->opts.pwalk_seq_gap
		    - fin->opts.pwalk_length;
		s_end = s_start + fin->opts.pwalk_length - 1;

		/*
		 * Check that the experiment overlaps the problem and isn't
		 * simply solving something different.
		 */
		if (s_start > prob_end || s_end < prob_start) {
		    printf("Experiment rejected as it does not overlap "
			   "the problem area\n");
		    continue;
		}

		/* Create experiment */
		exp = xrealloc(exp, ++nexp * sizeof(*exp));
		memset(&exp[nexp-1].r, 0, sizeof(GReadings));

		/*
		 * Compute new expected length.
		 */
		{
		    int alen;
		    int s_start2, s_end2;
		    alen = finish_avg_length(s_start, s_end, primer_dir,
					     t_start1, t_end1,
					     t_start2, t_end2,
					     &s_start2, &s_end2);
		    /*
		     * cost_mult *= ((fin->opts.pwalk_length + 1) /
		     *     (double)(alen + 1)-1)/LENGTH_SCALE+1;
		     */
		    s_start = s_start2;
		    s_end = s_end2;
		}

		/*
		 * Update cost to include distance from the problem.
		 * This should ideally be part of the quality estimation,
		 * so that we position the experiment such that the bases
		 * solving the problem(s) have the highest confidence.
		 */
		{
		    int dist = MIN(s_start - prob_start,
				   prob_end_orig - s_end);
		    if (dist > 0) {
			if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
			    printf("Solution not adjoining problem end => "
				   "Adjust cost by %f\n", 0.2 + 0.001*dist);
			cost_mult += 0.2 + 0.001*dist;
		    }
		}
		
		/* GReading component */
		exp[nexp-1].r.position = s_start;
		exp[nexp-1].r.sequence_length = s_end - s_start + 1;
		exp[nexp-1].r.start = 1;
		exp[nexp-1].r.end = 1 + exp[nexp-1].r.sequence_length + 1;
		exp[nexp-1].r.strand = r.strand; /* FIXME */
		exp[nexp-1].r.sense = primer_dir == 1 ? 0 : 1;
		exp[nexp-1].r.primer = temp_fwd
		    ? primer_dir == 1 ? 3 : 4
		    : primer_dir == 1 ? 4 : 3;
		exp[nexp-1].r.template = r.template;
		exp[nexp-1].r.chemistry = chem;

		/* Generic component */
		exp[nexp-1].type = EXPERIMENT_VPWALK;
		exp[nexp-1].nsolutions = fin->opts.pwalk_nsolutions;
		exp[nexp-1].score = 0;
		exp[nexp-1].cost = fin->cost[EXPERIMENT_VPWALK] * cost_mult;
		exp[nexp-1].expt_id = finish_next_expt_id(0);
		exp[nexp-1].group_id = group_id;
		exp[nexp-1].group_num = fin->opts.pwalk_ntemplates;
		exp[nexp-1].t_score = fin->tarr[r.template]->score;
		exp[nexp-1].t_dir = fin->tarr[r.template]->direction;
		exp[nexp-1].log_func = log_walk;

		/* e_walk component */
		/* WORKSHOP:
		 * prim[j] has 4 bytes padding characters at end.
		 */
		exp[nexp-1].data.e_walk = prim[j];
		exp[nexp-1].data.e_walk.secbind = pscore;

		if (fin->opts.debug[EXPERIMENT_VPWALK])
		    printf("read %d template %d: %d-%d (primer at %d)\n",
			   nexp-1, exp[nexp-1].r.template,
			   exp[nexp-1].r.position,
			   exp[nexp-1].r.position +
			   exp[nexp-1].r.sequence_length-1,
			   primer_pos_start);

		templates_picked[tnum] = i;

		tnum++;
	    }

	    if (tnum >= 2*fin->opts.pwalk_ntemplates)
		round = 2; /* no need for 2nd round */

	    if (fin->opts.pwalk_consistent_only)
		break;
	}

	if (tnum == 0) {
	    if (fin->opts.debug[EXPERIMENT_VPWALK]) {
		printf("No suitable templates found\n");
	    }
	} else {
	    noligos++;
	}
    }

    *nexp_p = nexp;
    return exp;
}

/*
 * experiment_walk
 *
 * The main externable function to produce a walk experiment for the problem
 * starting a position 'pos'. The problem region (for this solution type) may
 * extend either side of pos to prob_start and prob_end. This is used for
 * better positioning of solution experiments.
 *
 * Arguments:
 *	fin		Finish_t structure
 *	pos		Position in contig of problem base (1..LEN)
 *	chem		Suggested chemistry type
 *	dir		Suggested dir (1 == fwd, 2 == rev, 0 == either)
 *	prob_start	Extension of this problem region to the left
 *	prob_end	Extension of this problem region to the right
 *	nexp_p		Returned: number of experiments produced
 *
 * Returns:
 *	A malloced array of experiments if successful. Otherwise NULL.
 *	The size of the allocated array is returned in nexp_p.
 */
/*ARGSUSED*/
experiments_t *experiment_walk(finish_t *fin, int pos, int chem, int dir,
			       int prob_start, int prob_end, int *nexp_p,
			       int etype) {
    experiments_t *exp = NULL;
    int nexp = 0;
    int j;
    int primer_pos_start, primer_pos_end;
    int nprimers;
    experiment_walk_t *prim;
    int dir_ind[2], dir_loop;
    int prob_end_orig;

    printf(">>> PROBLEM AT %d - PRIMER WALK <<<\n", pos);

    /* Loop: dir_loop = 1 (dir==1), 2 (dir==2) or 1 and 2 (dir == 0) */
    dir_ind[0] = dir ? dir : 1;
    dir_ind[1] = dir ? 0 : 2;

    for (dir_loop = 0; dir_loop < 2; dir_loop++) {
	int primer_dir = dir_ind[dir_loop];
	if (!primer_dir)
	    break;

	if (fin->opts.debug[EXPERIMENT_VPWALK] > 1)
	    printf("primer_dir = %d\n", primer_dir);

	/* Find location to search for oligos */
	switch (primer_dir) {
	case 1: /* + strand */
	    primer_pos_start = pos - fin->opts.pwalk_offset1;
	    primer_pos_end = pos - fin->opts.pwalk_offset2;
	    break;

	case 2: /* - strand */
	    prob_end_orig = prob_end;
	    if (prob_end > pos + fin->opts.pwalk_length - fin->opts.pwalk_offset1)
		prob_end = pos + fin->opts.pwalk_length - fin->opts.pwalk_offset1;

	    primer_pos_start = prob_end + fin->opts.pwalk_offset2;
	    primer_pos_end   = prob_end + fin->opts.pwalk_offset1;
	    break;
	
	default: /* error */
	    fprintf(stderr, "Invalid primer direction\n");
	    return NULL;
	}

	/* Change positions from base numbers (1..N) to array
	 * indices (0..N-1)
	 */
	primer_pos_start--;
	primer_pos_end--;

	/* Generate the primers, moving the search window as appropriate */
	prim = NULL;
	for (j = 0; j < 10; j++) {
	    if (primer_pos_start < 0)
		primer_pos_start = 0;
	    if (primer_pos_end < 0)
		primer_pos_end = primer_pos_start + 40;
	    if (primer_pos_start >= io_clength(fin->io, fin->contig))
		primer_pos_start = io_clength(fin->io, fin->contig)-1;
	    if (primer_pos_end >= io_clength(fin->io, fin->contig))
		primer_pos_end = io_clength(fin->io, fin->contig)-1;

	    if (primer_pos_start >= primer_pos_end)
		break;

	    if (fin->opts.debug[EXPERIMENT_VPWALK])
		printf("Searching for primers between %d and %d\n",
		       primer_pos_start, primer_pos_end);

	    prim = find_primers(fin, io_clength(fin->io, fin->contig),
				primer_pos_start, primer_pos_end,
				primer_dir, &nprimers);
	    if (prim) {
		int new_nexp = nexp;
		int k;

		if (etype == EXPERIMENT_VPWALK) {
		    exp = find_templates(fin, prim, nprimers, primer_dir,
					 exp, &new_nexp, prob_start, prob_end,
					 prob_end_orig, chem);
		} else {
		    exp = generate_chr_exp(fin, prim, nprimers, primer_dir,
					   exp, &new_nexp, chem);
		}
		for (k = nexp; k < new_nexp; k++) {
		    /* prefer (just) closer primers */
		    exp[k].cost += ABS(j-1) * 0.01;
		}
		nexp = new_nexp;
		xfree(prim);
	    }
	    
	    if (fin->opts.debug[EXPERIMENT_VPWALK])
		printf("Expanding search range.\n");
	    if (primer_dir == 1) {
		if (primer_pos_start <= 0)
		    break;
		primer_pos_start -= 50;
		primer_pos_end -= 50;
	    } else {
		if (primer_pos_end >= io_clength(fin->io, fin->contig) - 1)
		    break;
		primer_pos_start += 50;
		primer_pos_end += 50;
	    }
	}
    }
    *nexp_p = nexp;
    return exp;
}
