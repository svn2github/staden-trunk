#include <stdio.h>

#include "IO.h"
#include "finish_long.h"
#include "finish_utils.h"
#include "xalloc.h"
#include "misc.h"

static void log_long(FILE *fp, finish_t *fin, experiments_t *e, int contig,
		     Tcl_DString *tagds) {
    char name[DB_NAMELEN+1];
    fprintf(fp, "=== Type Long gel reading\n");
    fprintf(fp, "=== Direction %c\n", "+-"[e->r.sense]);
    TextRead(fin->io, e->r.name, name, DB_NAMELEN);
    name[DB_NAMELEN] = 0;
    fprintf(fp, "=== Reading name %s\n", name);
}

static void log_reseq(FILE *fp, finish_t *fin, experiments_t *e, int contig,
		      Tcl_DString *tagds) {
    char name[DB_NAMELEN+1];
    fprintf(fp, "=== Type Resequence\n");
    fprintf(fp, "=== Direction %c\n", "+-"[e->r.sense]);
    TextRead(fin->io, e->r.name, name, DB_NAMELEN);
    name[DB_NAMELEN] = 0;
    fprintf(fp, "=== Reading name %s\n", name);
}

experiments_t *experiment_reseq(finish_t *fin, int pos,
				int chem, int dir, int *nexp_p,
				int islong) {
    experiments_t *exp = NULL;
    int nexp = 0;
    int rnum;
    GReadings r;
    int s_start, s_end;
    int s_start2, s_end2;
    int alen;
    int t_end1, t_start1;
    int t_end2, t_start2;
    int t_end, t_start;
    int e_len = islong ? fin->opts.long_length : fin->opts.reseq_length;

    if (fin->opts.debug[EXPERIMENT_RESEQ])
	printf("%s: pos=%d, chem=%d, dir=%d\n",
	       islong ? "long" : "reseq",
	       pos, chem, dir);

    /*
     * Find sequences that look to cover this region when performed as a
     * long gel reading.
     */
    for (rnum = io_clnbr(fin->io, fin->contig);
	 rnum;
	 rnum = io_rnbr(fin->io, rnum)) {
	double cost =
	    islong ? fin->cost[EXPERIMENT_LONG] : fin->cost[EXPERIMENT_RESEQ];

	if (io_relpos(fin->io, rnum) + e_len < pos)
	    continue; /* Not covering it yet */

	if (io_relpos(fin->io, rnum) - e_len > pos)
	    break; /* Beyond the region */

	/* Check sequence directions */
	s_start = io_length(fin->io, rnum) > 0
	    ? io_relpos(fin->io, rnum)
	    : io_relpos(fin->io, rnum) - io_length(fin->io, rnum) - 1;

	if (s_start < pos && io_length(fin->io, rnum) < 0)
	    continue; /* Going in the wrong direction */

	if (s_start > pos && io_length(fin->io, rnum) > 0)
	    continue; /* Going in the wrong direction */

	/* Check conditions */
	gel_read(fin->io, rnum, r);
	if (r.sense == 0 && dir == 2)
	    continue;
	if (r.sense == 1 && dir == 1)
	    continue;

	/* Skip specific templates */
	if (fin->template_skip && fin->template_skip[r.template])
	    continue;

	/* FIXME: Is it correct to only allow universal primer types? */
	if (PRIMER_TYPE_GUESS(r) == GAP_PRIMER_CUSTFOR ||
	    PRIMER_TYPE_GUESS(r) == GAP_PRIMER_CUSTREV)
	    continue;

	/* Increase cost as template score diminishes */
	cost *= 1 / fin->tarr[r.template]->score;

	/* Compute template positions */
	t_start1 = MIN(fin->tarr[r.template]->start,
		       fin->tarr[r.template]->end);
	t_end1 = MAX(fin->tarr[r.template]->start,
		     fin->tarr[r.template]->end);
	t_start2 = MIN(fin->tarr[r.template]->start2,
		       fin->tarr[r.template]->end2);
	t_end2 = MAX(fin->tarr[r.template]->start2,
		     fin->tarr[r.template]->end2);

	t_start = MIN(t_start1, t_start2);
	t_end   = MIN(t_end1, t_end2);

	/* Compute expected new sequence length */
	s_start = r.sense
	    ? r.position + r.sequence_length - 1 - e_len
	    : r.position;
	s_end = s_start + e_len - 1;

	/* Adjust expected sequence length based on template info */
	alen = finish_avg_length(s_start, s_end, r.sense,
				 t_start1, t_end1, t_start2, t_end2,
				 &s_start2, &s_end2);
	s_start = s_start2;
	s_end = s_end2;

	/* Adjust expected sequence length based on vector tags */
	finish_clip_svec(fin->io, &s_start, &s_end, rnum);

	if (fin->opts.debug[EXPERIMENT_RESEQ])
	    printf("read %c%d (%d): %d-%d\n",
		   "+-"[r.sense], nexp, rnum, 
		   s_start, s_end);

	/* Expected to cover region */
	exp = xrealloc(exp, ++nexp * sizeof(*exp));
	exp[nexp-1].r.name = r.name;
	exp[nexp-1].r.position = s_start;
	exp[nexp-1].r.sequence_length = s_end - s_start + 1;
	exp[nexp-1].r.start = 1;
	exp[nexp-1].r.end = 1 + exp[nexp-1].r.sequence_length + 1;
	exp[nexp-1].r.strand = r.strand;
	exp[nexp-1].r.sense = r.sense;
	exp[nexp-1].r.primer = r.primer;
	exp[nexp-1].r.template = r.template;
	exp[nexp-1].r.chemistry = chem;

	exp[nexp-1].type = islong
	    ? EXPERIMENT_LONG
	    : EXPERIMENT_RESEQ;
	exp[nexp-1].nsolutions = islong
	    ? fin->opts.reseq_nsolutions
	    : fin->opts.long_nsolutions;
	exp[nexp-1].score = 0;
	exp[nexp-1].cost = cost;
	exp[nexp-1].log_func = islong ? log_long : log_reseq;
	exp[nexp-1].expt_id = finish_next_expt_id(0);
	exp[nexp-1].group_id = finish_next_group_id(0);
	exp[nexp-1].group_num = 1;
	exp[nexp-1].t_score = fin->tarr[r.template]->score;
	exp[nexp-1].t_dir = fin->tarr[r.template]->direction;
    }

    *nexp_p = nexp;
    return exp;
}
