#include "haplo.h"
#include "IO.h"
#include "vseqs.h"
#include "gap_globals.h"
#include "xalloc.h"

/*
 * This file implements a consensus algorithm that operates only on a
 * specific set of templates within a contig. It is based on the virtual
 * contig code from Gap4 (and shared by prefinish).
 */

static vcontig_t *vcontig_from_templates(GapIO *io, int contig,
					 int *templates, int ntemplates) {
    vcontig_t *vc;
    int *present;
    int i;
    vrseq_t *vr, *next;

    /*
     * Create a vcontig_t, which initially has every reading in the contig
     * listed.
     */
    vc = new_vcontig(io, contig);

    /* Map the templates array to a lookup table indexed by template number */
    present = (int *)xcalloc(Ntemplates(io)+1, sizeof(int));
    for (i = 0; i < ntemplates; i++) {
	present[templates[i]] = 1;
    }

    /*
     * Now remove any vrseq in vc that isn't one of our listed templates.
     */
    for (vr = vc->left; vr; vr = next) {
	GReadings r;
	next = vr->right;
	gel_read(io, vr->rnum, r);
	if (!present[r.template]) {
	    del_vrseq(vc, vr);
	}
    }

    xfree(present);

    return vc;
}

/**
 * Computes the consensus of a contig using only readings from a specific
 * set of templates.
 *
 * start and/or end may be zero in which case they refer to the first and
 * last bases in the contig respectively.
 */
int calc_template_consensus(GapIO *io, int contig,
			    int start, int end,
			    int *templates, int ntemplates,
			    char **cons, float **qual) {
    int clen;
    vcontig_t *vc;

    if (!start) start = 1;
    if (!end) end = io_clength(io, contig);
    clen = end - start + 1;

    *cons = (char *)xmalloc(clen+1);
    if (qual) {
	*qual = (float *)xcalloc(clen+1, sizeof(float));
    }
    if (!*cons || (qual && !*qual))
	return -1;

    vc = vcontig_from_templates(io, contig, templates, ntemplates);

    calc_consensus(contig, start, end, CON_SUM,
		   *cons, NULL, qual ? *qual : NULL, NULL,
		   gap4_global_get_consensus_cutoff(),
		   gap4_global_get_quality_cutoff(),
		   virtual_info_func, (void *)vc);
    (*cons)[clen] = 0;

    del_vcontig(vc);

    return 0;
}

/**
 * Computes the template depth in 'contig' for coordinates start to end
 * inclusive. If start and/or end are zero then these are substituted for the
 * first and last base in the contig.
 *
 * The template depth is returned in tdepth which should be allocated by
 * the caller to be end-start+1 elements long.
 *
 * Returns: max depth on success
 *         -1 on failure
 */
int calc_template_depth(GapIO *io, int contig, int start, int end,
			int *tdepth) {
    int rnum, max = 0;
    int *t_end; /* Next valid position to increment for this template */

    if (start == 0)
	start = 1;
    if (end == 0)
	end = io_clength(io, contig);

    if (NULL == (t_end = (int *)xcalloc(Ntemplates(io)+1, sizeof(int))))
	return -1;

    for (rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum)) {
	GReadings r;
	int tnum, r_start, r_end, i;

	gel_read(io, rnum, r);
	tnum = r.template;

	r_start = r.position;
	r_end = r.position + r.sequence_length-1;

	/* Optimisation for readings outside of start..end range */
	if (r_start > end)
	    break;
	if (r_end < start)
	    continue;

	/*
	 * Compute absolute coordinates for this reading. We clip based on the
	 * input start/end and also on the last position written to for this
	 * template. This means that overlapping sequences on a single template
	 * will not be counted more than once.
	 */
	if (r_start < start)
	    r_start = start;
	if (r_end > end)
	    r_end = end;
	if (r_start < t_end[tnum])
	    r_start = t_end[tnum];

	for (i = r_start; i <= r_end; i++) {
	    tdepth[i - start]++;
	    if (max < tdepth[i - start])
		max = tdepth[i - start];
	}
	
	t_end[tnum] = i;
    }

    xfree(t_end);

    return max;
}
