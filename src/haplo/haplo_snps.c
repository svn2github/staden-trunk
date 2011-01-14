/*
 * This produces a Tcl list of predicted SNP locations for a given contig
 * It does so by using the gap4 consensus algorithm to obtain 4 scores per
 * base call. This is then adjusted for by checking for mononucleotide
 * repeas as these often have either misalignments or fluctuations in
 * copy-number from PCR.
 */

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <unistd.h>

#include "IO.h"
#include "qual.h"
#include "gap_globals.h"
#include "misc.h"
#include "template.h"
#include "dstring.h"
#include "haplo.h"
#include "notes.h"

/* ------------------------------------------------------------------------ */
/**
 * Use calc_discrepancies() to identify candidate SNP locations within
 * a specific contig.
 * 'Cutoff' is a threshold above (or equal) with the combined scores from
 * calc_discrepancies should be. A value of '40' works well.
 * start and end may independently be zero, in which case they default to
 * the contig start (1) and end (length). Otherwise they are the inclusive
 * ranges in which to look for SNPs.
 *
 * Two_alleles is a boolean flag to control whether to the input data
 * should match a 50/50 ratio from two alleles. If true a simple
 * binomial rule is used to modify the 2nd-highest confidence value,
 * otherwise this value is used as-is.
 *
 * Returns a malloced array of snp_t structures holding the SNPs.
 *         NULL on failure.
 */
static snp_t *candidate_snps(GapIO *io, int contig, int start, int end,
			     int *nsnps_p, int cutoff,
			     int min_base_qual, int two_alleles) {
    float *q1, *q2;
    int clen;
    int i;
    snp_t *snps = NULL;
    int nsnps = 0;
    int *costs;

    clen = io_clength(io, contig);
    if (!start) start = 1;
    if (!end) end = clen;

    /* Calculate discrepancies matrix */
    if (NULL == (q1 = (float *)xmalloc((end-start+1) * sizeof(float))))
	return NULL;

    if (NULL == (q2 = (float *)xmalloc((end-start+1) * sizeof(float)))) {
	xfree(q1);
	return NULL;
    }

    calc_discrepancies(contig, start, end, q1, q2,
		       consensus_cutoff, min_base_qual,
		       database_info, (void *)io);


    /* Search for high scoring SNPs and record position/score */
    for (i = 0; i < (end-start+1); i++) {
	float combined;

	/* combined = (q1[i] + log(MAX(.00001, 1-q2[i]/3))*-500)/2; */
	combined = two_alleles ? MIN(1000, q1[i] * q2[i]) : q1[i];
	if (combined >= cutoff) {
	    /* Grow snps array in blocks of 32 at a time */
	    if (nsnps % 32 == 0) {
		snps = (snp_t *)xrealloc(snps, (nsnps + 32) * sizeof(snp_t));
		if (!snps)
		    return NULL;
	    }
	    snps[nsnps].pos   = i+start;
	    snps[nsnps].score = combined;
	    snps[nsnps].seqs  = NULL;
	    snps[nsnps].nseqs = 0;
	    nsnps++;
	}
    }

    xfree(q1);
    xfree(q2);

    /*
     * Reduce score of neighbouring snps. Each SNP is a triangular impulse
     * of 0, 1, 2, 3, 2, 1, 0 centred on the SNP location. We then divide
     * the SNP score by 3/costs[?] so that neighbouring SNPs get reduced.
     */
    if (NULL == (costs = (int *)xcalloc(clen+1, sizeof(int))))
	return NULL;

    for (i = 0; i < nsnps; i++) {
	int j;
	for (j = 0; j < 3; j++) {
	    if (snps[i].pos - j >= 1 && snps[i].pos - j <= clen)
		costs[snps[i].pos - j] += 3-j;
	    if (j != 0 && (snps[i].pos + j >= 1 && snps[i].pos + j <= clen))
		costs[snps[i].pos + j] += 3-j;
	}
    }

    for (i = 0; i < nsnps; i++) {
	snps[i].score /= costs[snps[i].pos] / 3.0;
    }
    xfree(costs);


    *nsnps_p = nsnps;
    return snps;
}


/**
 * Produce a list of all templates overlapping base 'pos' in 'contig'.
 * If a sequence at the specific position has a quality lower than the global
 * quality_cutoff then it is ignored. If a template is represented more than
 * once (ie by more than one sequence) then the sequences are checked for
 * matches and mismatches. Mismatching templates will be rejected if the
 * match is high. Matching sequences on a templates get coalesced into one
 * record.
 *
 * Returns an array of seq_base_t structs on success, of length *nseqs_p.
 *            NULL on failure
 */
static seq_base_t *seqs_at_region(GapIO *io, int contig, template_c **tarr,
				  int pos, int *nseqs_p, int min_qual) {
    int rnum;
    seq_base_t *seqs = NULL;
    int nseqs = 0;
    GReadings r;
    int i, j;

    /*
     * Loop through all sequences in the contig. Could be optimised as
     * this gets called multiple times, but for now take the easy
     * approach.
     */
    for (rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum)) {
	double tscore;
	int1 conf;
	char *seq;

	/* Read pos are sorted - optimise */
	if (io_relpos(io, rnum) > pos)
	    break;

	if (io_relpos(io, rnum) + ABS(io_length(io, rnum)) <= pos)
	    continue;

	/* Skip if it's a reference sequence */
	if (find_note(io, rnum, "REFS"))
	    continue;

	/* Read reading and sequence details */
	gel_read(io, rnum, r);
	if (r.confidence) {
	    int1 *c = (int1 *)xmalloc(r.length);
	    if (c) {
		DataRead(io, r.confidence, c, r.length, 1);
		conf = c[pos - r.position + r.start];
		xfree(c);
	    } else {
		conf = 0;
	    }
	} else {
	    conf = 0;
	}
	seq = SeqReadStatic(io, r.sequence, r.length);

	/* Filter out sequences of poor quality */
	if (conf < min_qual)
	    continue;

	if (!r.template)
	    continue;

	/* Penalise for poor quality templates */
	tscore = 1;
	if (r.template) {
	    int consist;

	    consist = tarr[r.template]->consistency;
	    if (consist & (TEMP_CONSIST_STRAND | TEMP_CONSIST_PRIMER))
		tscore /= 10;
	    else if (consist)
		tscore /= 3;
	}

	/* Grow in blocks of 8 */
	if (nseqs % 8 == 0) {
	    seqs = (seq_base_t *)xrealloc(seqs,
					  (nseqs + 8) * sizeof(seq_base_t));
	    if (!seqs)
		return NULL;
	}

	seqs[nseqs].tmplate = r.template;
	seqs[nseqs].tscore  = tscore;
	seqs[nseqs].conf    = conf;
	seqs[nseqs].base    = toupper(seq[pos - r.position + r.start]);
	nseqs++;
    }

    /*
     * We now have a list of sequences, but we want to remove any duplicate
     * templates. We mark sequences for removal by setting tscore to 0.
     */
    for (i = 0; i < nseqs; i++) {
	int t = seqs[i].tmplate;
	int A = 0, C = 0, G = 0, T = 0, X = 0;

	if (!seqs[i].tscore)
	    continue;

	/*
	 * Compute a simple consensus + score.
	 */
	for (j = i; j < nseqs; j++) {
	    if (seqs[j].tmplate == t) {
		switch (seqs[j].base) {
		case 'A': case 'a':
		    A += seqs[j].conf;
		    break;
		case 'C': case 'c':
		    C += seqs[j].conf;
		    break;
		case 'G': case 'g':
		    G += seqs[j].conf;
		    break;
		case 'T': case 't':
		    T += seqs[j].conf;
		    break;
		case '*':
		    X += seqs[j].conf;
		    break;
		}
		/* Mark all other uses of this template as to be filtered */
		if (j > i)
		    seqs[j].tscore = 0;
	    }
	}

	/* And update first occurance of this template */
	if (A -C-G-T-X > 0) {
	    seqs[i].base = 'A';
	    seqs[i].conf = A-C-G-T-X;
	} else if (C -A-G-T-X > 0) {
	    seqs[i].base = 'C';
	    seqs[i].conf = C-A-G-T-X;
	} else if (G -A-C-T-X > 0) {
	    seqs[i].base = 'G';
	    seqs[i].conf = G-A-C-T-X;
	} else if (T -A-C-G-X > 0) {
	    seqs[i].base = 'T';
	    seqs[i].conf = T-A-C-G-X;
	} else if (X -A-C-G-T > 0) {
	    seqs[i].base = '*';
	    seqs[i].conf = X-A-C-G-T;
	} else {
	    seqs[i].tscore = 0;
	}

	if (seqs[i].conf < min_qual)
	    seqs[i].tscore = 0;
	/*
	printf("Pos %d Template %d A%d C%d G%d T%d X%d =>%c %d\n",
	       pos, t, A, C, G, T, X, seqs[i].base, seqs[i].conf);
	*/
    }

    /* Finally remove the templates with tscore == 0 */
    for (i = j = 0; i < nseqs; i++) {
	if (seqs[i].tscore) {
	    seqs[j++] = seqs[i];
	}
    }

    *nseqs_p = j;
    return seqs;
}


/**
 * Given a set of snps within a contig, this fills out the seqs/nseqs
 * component of the snp_t structure to list the sequences at that
 * base position.
 */
static void add_snp_bases(GapIO *io, int contig, template_c **tarr,
			  snp_t *snp, int nsnps, int min_qual) {
    int i;

    for (i = 0; i < nsnps; i++) {
	int nseqs;
	snp[i].seqs  = seqs_at_region(io, contig, tarr, snp[i].pos, &nseqs,
				      min_qual);
	snp[i].nseqs = snp[i].seqs ? nseqs : 0;
    }
}

/**
 * Rescaore SNPs based on mononucleotide runs.
 * SNPs occuring within these may be false for two reasons.
 * 1) Bad alignments causing pads to be placed in different locations.
 * 2) Variation in the copy-number of a base produced via PCR issues.
 *
 * The snps array is modified inline.
 */
static void mononucleotide_scores(GapIO *io, int contig,
				  snp_t *snps, int nsnps) {
    int i, j;
    char *cons;
    int clen = io_clength(io, contig);

    cons = (char *)xmalloc(clen+1);
    calc_consensus(contig, 1, clen, CON_SUM,
		   cons, NULL, NULL, NULL,
		   consensus_cutoff, quality_cutoff,
		   database_info, (void *)io);

    for (i = 0; i < nsnps; i++) {
	char base;
	int count;

	/*
	 * Count the number of bases matching the SNP consensus either side.
	 * We skip pads too. The reason for this is to work out if the SNP
	 * occurs within a run of the same base type so we can downweight the
	 * impact if it is a base vs pad SNP. The alignments are often poor
	 * in such cases and also there seems to be copy-number variations
	 * introduced by PCR that are not necessarily a real component of the
	 * two alleles.
	 */

	/* Pads are a special case - compute base type as longest from l/r*/
	if ((base = cons[snps[i].pos-1]) == '*') {
	    char basel = '*', baser = '*';
	    int countl = 0, countr = 0;

	    /* Find left and right base type */
	    j = snps[i].pos-2;
	    while (j >= 0 && cons[j] == '*')
		j--;
	    if (j >= 0)
		basel = cons[j];
	    while (j >= 0 && cons[j] == basel)
		j--, countl++;

	    j = snps[i].pos;
	    while (j < clen && cons[j] == '*')
		j++;
	    if (j < clen)
		baser = cons[j];
	    while (j < clen && cons[j] == baser)
		j++, countr++;

	    if (countl > countr && basel != '*')
		base = basel;
	    else if (countr > countl && baser != '*')
		base = baser;
	    else if (basel != '*')
		base = basel;
	    else
		base = baser;
	} else {
	    base = cons[snps[i].pos-1];
	}

	count = 0;
	for (j = snps[i].pos-1; j >= 0; j--) {
	    if (cons[j] == '*')
		continue;
	    
	    if (cons[j] == base)
		count++;
	    else
		break;
	}
	for (j = snps[i].pos; j < clen; j++) {
	    if (cons[j] == '*')
		continue;
	    
	    if (cons[j] == base)
		count++;
	    else
		break;
	}

	snps[i].score /= count;
    }

    xfree(cons);
}

/**
 * Returns a dstring containing "{position score <template>...} ..."
 * where <template> is "{template_num template_score basecall conf}".
 */
static dstring_t *list_snps(GapIO *io, int contig, snp_t *snps, int nsnps,
			     template_c **tarr, int snp_cutoff) {
    int i, j;
    dstring_t *ds;

    ds = dstring_create(NULL);
    for (i = 0; i < nsnps; i++) {
	if (snps[i].score < snp_cutoff)
	    continue;

	dstring_appendf(ds, "{%d %f", snps[i].pos, snps[i].score);

	for (j = 0; j < snps[i].nseqs; j++) {
	    dstring_appendf(ds, " {%d %f %c %d}",
			    snps[i].seqs[j].tmplate,
			    snps[i].seqs[j].tscore,
			    snps[i].seqs[j].base,
			    snps[i].seqs[j].conf);
	}

	dstring_appendf(ds, "} ");
    }

    return ds;
}

/**
 * Frees a snp_t array of size nsnps, also freeing linked seq_base_t structs
 * too.
 */
void free_snps(snp_t *snps, int nsnps) {
    if (snps) {
	int i;

	for (i = 0; i < nsnps; i++) {
	    if (snps[i].seqs)
		xfree(snps[i].seqs);
	}
	xfree(snps);
    }
}


/* ------------------------------------------------------------------------ */
/*
 * MAIN ENTRY POINT.
 *
 * Given a contig containing a mixed assembly from two alleles, this
 * function attempts to split templates into two sets corresponding
 * to each of the two alleles.
 *
 * Returns a dstring containing "{position score <template>...} ..."
 * where <template> is "{template_num template_score basecall conf}".
 * Returns NULL on failure.
 */
dstring_t *haplo_snps(GapIO *io, int contig, int start, int end,
		      int discrep_cutoff, int snp_cutoff,
		      int min_base_qual, int two_alleles)
{
    template_c **tarr = NULL;
    snp_t *snps = NULL;
    int nsnps = 0;
    dstring_t *ds = NULL;

    /* Compute the discrepancies and obtain possible SNP locations & scores */
    snps = candidate_snps(io, contig, start, end, &nsnps, discrep_cutoff,
			  min_base_qual, two_alleles);
    if (!nsnps)
	goto error;

    /* Compute template consistency status; used for scoring SNPs */
    if (NULL == (tarr = init_template_checks(io, 1, &contig, 1)))
	goto error;
    check_all_templates(io, tarr);

    /* Identify template basecalls at each snp site */
    add_snp_bases(io, contig, tarr, snps, nsnps, min_base_qual);

    /* Adjust scores in mononucleotide runs */
    mononucleotide_scores(io, contig, snps, nsnps);

    /* Summary list of snps & templates */
    ds = list_snps(io, contig, snps, nsnps, tarr, snp_cutoff);

    /* Free up memory */
 error:
    if (snps)
	free_snps(snps, nsnps);

    if (tarr)
	uninit_template_checks(io, tarr);

    return ds;
}
