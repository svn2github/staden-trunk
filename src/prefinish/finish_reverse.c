#include <stdio.h>
#include <math.h>

#include "IO.h"
#include "finish_reverse.h"
#include "finish_utils.h"
#include "xalloc.h"
#include "misc.h"

static void log_reverse(FILE *fp, finish_t *fin, experiments_t *e, int contig,
			Tcl_DString *tagds) {
    fprintf(fp, "=== Type Reverse reading\n");
    fprintf(fp, "=== Direction %c\n", "+-"[e->r.sense]);
}

/*
 * Computes the chance of a sequence from position 'seq_min' to 'seq_max'
 * within a template claiming to have insert sizes ranging from is_min
 * to is_max covering a problem located at position p_min to p_max.
 */
static double insert_end_chance(int is_min, int is_max,
				int seq_min, int seq_max,
				int p_min, int p_max,
				int end /* 0 == left, 1 == right */)
{
#define SIZED 10000
    double size_dist[SIZED];
    int i, pos;
    double prob_start;
    int probs_covered;
    double total_prob;
    int seq_len = seq_max - seq_min + 1;
    int spos, epos;

    /* For now, frabricate a distribution */
    for (i = 0; i < SIZED; i++) {
	size_dist[i] = 0;
    }
    {
	int size = is_max - is_min; 
	int mid = is_min + size/2; 
	double tot;
	int i; 
 
	for (tot = 0, i = is_min; i < is_max; i++)
	    tot += sqrt(size*size/4 - (i-mid)*(i-mid));

	for (i = is_min; i < is_max; i++)
	    size_dist[i] = sqrt(size*size/4 - (i-mid)*(i-mid)) / tot;
    }

    /*
     * Iterate around all positions 'pos' within possible template
     * sizes. Foreach pos that covers p_min to p_max add probability of this
     * being a valid starting point (size_dist[pos]) multiplied by the
     * percentage of p_min to p_max covered (in the case when it's longer than
     * one base) and the chance of a high quality.
     */
    /* TODO - need reorganisation of primary code in finish_main.c */
    
    total_prob = 0;
    if (end == 0) {
	spos = MAX(p_min - seq_len, 0);
	epos = MIN(p_max, SIZED);
    } else {
	spos = MAX(p_min, 0);
	epos = MIN(p_max + seq_len, SIZED);
    }
    for (pos = spos; pos < epos; pos++) {
	/* Work out probability of sequence starting at this point */
	prob_start = size_dist[pos];

	if (prob_start != 0) {
	    probs_covered = 0;
	    for (i = pos; i <= pos + seq_len; i++) {
		if (i >= p_min && i <= p_max)
		    probs_covered++;
	    }
	} else {
	    probs_covered = 0;
	}

	total_prob += prob_start * probs_covered / (double)(p_max - p_min + 1);
    }

    return total_prob;
}

experiments_t *experiment_reverse(finish_t *fin, int pos,
				  int chem, int dir,
				  int prob_start, int prob_end,
				  int *nexp_p) {
    experiments_t *exp = NULL;
    int nexp = 0;
    int i;
    int tnum;
    int templates_picked[100];

    printf("reverse: pos=%d, chem=%d, dir=%d\n", pos, chem, dir);

    /*
     * Find templates that have only one end known and the other end likely
     * covering this region.
     */
    tnum = 0;
    memset(templates_picked, 0, 100 * sizeof(int));
    for (i = 1; i <= Ntemplates(fin->io); i++) {
	item_t *clist;
	int have_fwd, have_rev;
	int seq_start, seq_end; /* contig position of candidate sequence */
	int seq_left, seq_right;
	GTemplates t;
	GReadings r;
	char *tname;
	double cost = fin->cost[EXPERIMENT_REVERSE];
	int t_start1, t_start2, t_start;
	int t_end1, t_end2, t_end;
	double cov_chance;

	/* Template has to at least be in this contig in part */
	if (!fin->tarr[i])
	    continue;

	/*
	 * Only use consistent templates.
	 * Allowing inconsistent templates will mean the TEMP_DIRECTION
	 * macros below fail.
	 */
	if (fin->tarr[i]->consistency)
	    continue;

	/* Already have both ends used? */
	have_fwd = have_rev = 0;
	for (clist = head(fin->tarr[i]->gel_cont); clist;
	     clist = clist->next) {
	    gel_cont_t *gc = (gel_cont_t *)(clist->data);

	    gel_read(fin->io, gc->read, r);
	    switch (PRIMER_TYPE_GUESS(r)) {
	    case GAP_PRIMER_FORWARD:
		have_fwd = 1;
		break;
	    case GAP_PRIMER_REVERSE:
		have_rev = 1;
		break;
	    }
	}
	if (have_fwd && have_rev)
	    continue;

	/* Don't overuse the template */
	if (fin->template_used[i] >= fin->opts.pwalk_use_template)
	    continue;

	/*
	 * Is the template in the correct position?
	 * tarr[i]->start/end assume using the minimum insert size
	 * tarr[i]->start2/end2 assumed the max insert size.
	 */
	t_start1 = MIN(fin->tarr[i]->start, fin->tarr[i]->end);
	t_end1 = MAX(fin->tarr[i]->start, fin->tarr[i]->end);
	t_start2 = MIN(fin->tarr[i]->start2, fin->tarr[i]->end2);
	t_end2 = MAX(fin->tarr[i]->start2, fin->tarr[i]->end2);

	t_start = MIN(t_start1, t_start2);
	t_end   = MAX(t_end1,   t_end2);

	/* Does the template cover the problem region? */
	if (t_start > prob_end || t_end < prob_start)
	    continue;

	/*
	 * Find the most probable position for the sequence, which is
	 * computed by averaging the min and max template sizes.
	 *
	 * FIXME: this isn't what we need...
	 */
	if (have_fwd) {
	    if (TEMP_DIRECTION(fin->tarr[i]) == 0 /* +ve */) {
		seq_end = (fin->tarr[i]->end + fin->tarr[i]->end2)/2;
		seq_start = seq_end - fin->opts.reseq_length;
	    } else {
		seq_start = (fin->tarr[i]->end + fin->tarr[i]->end2)/2;
		seq_end = seq_start + fin->opts.reseq_length;
	    }
	} else {
	    if (TEMP_DIRECTION(fin->tarr[i]) == 0 /* +ve */) {
		seq_start = (fin->tarr[i]->start + fin->tarr[i]->start2)/2;
		seq_end = seq_start + fin->opts.reseq_length;
	    } else {
		seq_end = (fin->tarr[i]->start + fin->tarr[i]->start2)/2;
		seq_start = seq_end - fin->opts.reseq_length;
	    }
	}

	/*
	 * When seeing if this solution covers the problem we need to allow
	 * for tolerance. We maybe do not know the template position exactly,
	 * so seq_start/end is the most probable position only. We score
	 * highest for this and lower the score the further away the solution
	 * is (by increasing cost).
	 * Obviously tight size ranges on templates are better than large
	 * ranges, so we take this into account too.
	 *
	 * We have a specific problem we want to solve, so we work out the
	 * probability of this sequence covering this position and score it
	 * base on what other problems it also solves, rather than scoring
	 * base on the most probable position for this sequence or by summing
	 * all the possible positions with their associated probabilities and
	 * the problems solved at these positions. (This latter technique
	 * sounds the most accurate, but it doesn't fit our "have problem,
	 * find solution" mode of working.)
	 */

	/* Pick sequence position centred on problem */
	if ((have_fwd && TEMP_DIRECTION(fin->tarr[i]) == 0 /* +ve */) ||
	    (!have_fwd && TEMP_DIRECTION(fin->tarr[i]) != 0 /* -ve */))
	{
	    seq_right = prob_start + (prob_end - prob_start) / 2
		+ fin->opts.reseq_length / 2;
	    seq_left = seq_right - fin->opts.reseq_length;
	} else {
	    seq_left = prob_start + (prob_end - prob_start) / 2
		- fin->opts.reseq_length / 2;
	    seq_right = seq_left + fin->opts.reseq_length;
	}

	/* And shift to be within the template if necessary */
	if (seq_left < t_start) {
	    int dist = t_start - seq_left;
	    seq_left += dist;
	    seq_right += dist;
	}
	if (seq_right > t_end) {
	    int dist = seq_right - t_end;
	    seq_left -= dist;
	    seq_right -= dist;
	}

	printf("prob at %d-%d Template %d ends %d to %d, exp %d-%d\n",
	       prob_start, prob_end, i, t_start, t_end,
	       seq_left, seq_right);
	/*
	 * Compute score for likelihood of this experiment covering this
	 * position.
	 */
	template_read(fin->io, i, t);
	
	if (have_fwd) {
	    if (fin->tarr[i]->direction == 0) {
		cov_chance = 
		    insert_end_chance(t.insert_length_min,
				      t.insert_length_max,
				      seq_start - t_start,
				      seq_end - t_start,
				      prob_start - t_start,
				      prob_end - t_start,
				      1);
	    } else {
		cov_chance = 
		    insert_end_chance(t.insert_length_min,
				      t.insert_length_max,
				      t_end - seq_end,
				      t_end - seq_start,
				      t_end - prob_start,
				      t_end - prob_end,
				      1);
	    }
	} else {
	    if (fin->tarr[i]->direction == 0) {
		cov_chance = 
		    insert_end_chance(t.insert_length_min,
				      t.insert_length_max,
				      seq_start - t_start,
				      seq_end - t_start,
				      prob_start - t_start,
				      prob_end - t_start,
				      0);
	    } else {
		cov_chance = 
		    insert_end_chance(t.insert_length_min,
				      t.insert_length_max,
				      t_end - seq_end,
				      t_end - seq_start,
				      t_end - prob_start,
				      t_end - prob_end,
				      0);
	    }
	}
	printf("chance of coverage = %f\n", cov_chance);
	if (cov_chance == 0)
	    continue;
	cost /= cov_chance;

	/* Debug info */
	tname = TextAllocRead(fin->io, t.name);
	printf("Template %c %d/%s. Reads:",
	       "-+"[TEMP_DIRECTION(fin->tarr[i])],
	       i, tname);
	xfree(tname);

	printf("i=%d\n", i);
	for (clist = head(fin->tarr[i]->gel_cont); clist;
	     clist = clist->next) {
	    gel_cont_t *gc = (gel_cont_t *)(clist->data);
	    printf("%d ", gc->read);
	}
	putchar('\n');
	    
	printf("reverses: %s seq at pos %d-%d\n",
	       have_fwd ? "rev" : "fwd", seq_start, seq_end);

	/*
	 * From this point on it sounds like a viable experiment, so we
	 * now start working out the score and groupings.
	 */

	/* Increase cost as template score diminishes */
	cost *= 1/fin->tarr[i]->score;

	/* Is it a contaminant template? */
	if (fin->template_dup && fin->template_dup[i]) {
	    if (template_is_dup(fin, templates_picked, tnum, i)) {
		printf("Template dup %d, lowering score\n", i);
		cost *= fin->opts.pwalk_dup_template_cost;
	    }
	}

	printf("read %c%d: %d-%d\n",
	       "+-"[r.sense], nexp, 
	       seq_start, seq_end);

	/* Expected to cover region */
	exp = xrealloc(exp, ++nexp * sizeof(*exp));
	exp[nexp-1].r.name = 0;
	exp[nexp-1].r.position = seq_start;
	exp[nexp-1].r.sequence_length = seq_end - seq_start + 1;
	exp[nexp-1].r.start = 1;
	exp[nexp-1].r.end = 1 + exp[nexp-1].r.sequence_length + 1;
	exp[nexp-1].r.strand = have_fwd ? 0 : 1;
	exp[nexp-1].r.sense = (have_fwd != exp[nexp-1].r.strand)
	    ? GAP_SENSE_ORIGINAL : GAP_SENSE_REVERSE;
	exp[nexp-1].r.primer = have_fwd
	    ? GAP_PRIMER_REVERSE : GAP_PRIMER_FORWARD;
	exp[nexp-1].r.template = r.template;
	exp[nexp-1].r.chemistry = chem;

	exp[nexp-1].type = EXPERIMENT_REVERSE;
	exp[nexp-1].score = 0;
	exp[nexp-1].cost = cost;
	exp[nexp-1].log_func = log_reverse;
	exp[nexp-1].expt_id = finish_next_expt_id(0);
	exp[nexp-1].group_id = finish_next_group_id(0);
	exp[nexp-1].group_num = 1;
	exp[nexp-1].t_score = fin->tarr[r.template]->score;
	exp[nexp-1].t_dir = TEMP_DIRECTION(fin->tarr[i]);

	templates_picked[tnum] = i;
	tnum++;

	if (tnum >= 100)
	    break;
    }

    *nexp_p = nexp;
    return exp;
}
