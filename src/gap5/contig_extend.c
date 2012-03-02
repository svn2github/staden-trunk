/*
 * Assemblies produced by NGS tools often collapse data to an internal
 * representation during assembly, and then map reads back to the 
 * consensus sequences produced during assembly.
 *
 * At that stage we see oddities though, such as a consensus ending early
 * and dozens of reads with cutoff that agrees perfectly. For example,
 * with lowercase being cutoff:
 *
 * Consensus: AGCTATGAGGCGATACGATCCGTAA
 * Read1:     AGCTATGAGGCGATACGATCCGTAAgata
 * Read2:     AGCTATGAGGCGATACGATCCGTAAgatat
 * Read3:     AGCTATGTGGCGATACGATCCGTaagatct
 * Read4:     AGCTATGAGGCGATACGATCCGTAAgtatg
 * Read5:     AGCTATGAGGCGATACG-TCCGTAAgatatgcgtt
 * Read6:     AGCTATGAGGCGATACGATCCGTAAgatatgcgttaa
 * Read7:     AGCTATGAGGCGATACGATCCGTAAgatatgcgttaagatcg
 *
 * In this case we'd be happy to extend to gatatgcgtt say (3 fold depth),
 * but leaving the trailing end clipped still as it's not got deep enough
 * coverage.
 *
 * Also we'd extend reads 1,2,3,5,6,7 but not read4 due to an indel.
 * Note read3 has a different clip point than the others, but that is fine
 * as it aligns still.
 */

#include <string.h>
#include <stdio.h>
#include <tg_gio.h>
#include <ctype.h>

#include "dna_utils.h"
#include "contig_extend.h"
#include "consensus.h"
#include "text_output.h"
#include "gap4_compat.h"  /* complement_contig */

#define CSZ 1024 /* consensus size - how far back in contig we'll look */
#define ESZ 1024 /* extension size - how far we can extend */

/*
 * Extends the right hand end of a single contig.
 *
 * Min_depth is the minimum depth for extension. If lower then even if the
 * data matches we'll not extend further.
 *
 * Match_score (+ve) and mismatch_score (-ve) are accumulated during
 * extension to ensure that we don't extend into junk mismatching DNA.
 */
static int contig_extend_single(GapIO *io, tg_rec crec, int dir, int min_depth,
				int match_score, int mismatch_score) {
    int end;
    rangec_t *r;
    int nr, i;
    contig_t *c;
    char cons[CSZ], new_cons[ESZ];
    int freqs[ESZ][4], depth[ESZ];
    double score, best_score;
    int best_pos, nseq;

    vmessage("Processing contig #%"PRIrec", %s end\n",
	     crec, dir ? "left" : "right");

    for (i = 0; i < ESZ; i++) {
	freqs[i][0] = freqs[i][1] = freqs[i][2] = freqs[i][3] = 0;
	depth[i] = 0;
    }

    c = cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

    if (consensus_valid_range(io, crec, NULL, &end) != 0) {
	cache_decr(io, c);
	return -1;
    }

    calculate_consensus_simple(io, crec, end-(CSZ-1), end, cons, NULL);

    /* Start */
    /* Not implemented for now: rev complement and go again! */

    /* End */
    r = contig_seqs_in_range(io, &c, end, end, 0, &nr);
    if (!r) {
	cache_decr(io, c);
	return -1;
    }

    for (i = 0; i < nr; i++) {
	seq_t *s = cache_search(io, GT_Seq, r[i].rec);
	seq_t *sorig = s;
	int cstart, cend;
	int j, k, slen;

	if ((s->len < 0) ^ r[i].comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	cstart = r[i].start + s->left-1;
	cend   = r[i].start + s->right-1;

	/* Does cutoff extend to contig end, if so does it match cons? */
	if (cend < end) {
	    int mis = 0, len = 0;
	    if (end - cend >= CSZ) {
		/*
		fprintf(stderr,"Skipping #%"PRIrec" due to length of cutoff\n",
			r[i].rec);
		*/
		if (sorig != s)
		    free(s);
		r[i].rec = 0; /* Mark for removal */
		continue;
	    }

	    for (k = s->right, j = cend+1; j <= end; j++, k++) {
		//printf("%d: %c %c\n", j, s->seq[k], cons[j-(end-(CSZ-1))]);
		if (s->seq[k] != cons[j-(end-(CSZ-1))])
		    mis++;
	    }
	    len = end - cend;
	    if (100*mis/len > 5) {
		/*
		fprintf(stderr, "Skipping #%"PRIrec" due to high disagreement "
			"with consensus.\n", r[i].rec);
		*/
		if (sorig != s)
		    free(s);
		r[i].rec = 0;
		continue;
	    }
	}

	/* So we got here, let's accumulate extension stats */
	slen = ABS(s->len);
	for (k = 0, j = end+1 - r[i].start; j < slen && k < ESZ; j++, k++) {
	    //printf("%d: %c\n", j + r[i].start, s->seq[j]);
	    if(s->seq[j] == 'N')
		continue;

	    freqs[k][dna_lookup[s->seq[j]]]++;
	    depth[k]++;
	}

	if (sorig != s)
	    free(s);
    }

    score = best_score = 0;
    best_pos = 0;
    
    for (i = 0; i < ESZ; i++) {
	int call, best = 0, j;
	double dd;

	if (depth[i] < min_depth)
	    break;

	for (j = 0; j < 4; j++) {
	    if (best < freqs[i][j]) {
		best = freqs[i][j];
		call = j;
	    }
	}
	new_cons[i] = "ACGT"[call];

	dd = (double)depth[i];
	switch (call) {
	case 0:
	    score +=  freqs[i][0] / dd;
	    score -= (freqs[i][1] + freqs[i][2] + freqs[i][3]) / dd;
	    break;
	case 1:
	    score +=  freqs[i][1] / dd;
	    score -= (freqs[i][0] + freqs[i][2] + freqs[i][3]) / dd;
	    break;
	case 2:
	    score +=  freqs[i][2] / dd;
	    score -= (freqs[i][0] + freqs[i][1] + freqs[i][3]) / dd;
	    break;
	case 3:
	    score +=  freqs[i][3] / dd;
	    score -= (freqs[i][0] + freqs[i][1] + freqs[i][2]) / dd;
	    break;
	}

	if (best_score <= score) {
	    best_score = score;
	    best_pos = i+1;
	}
	/*
	printf("%3d %3d\t%c\t%3d %3d %3d %3d %7.1f\n",
	       i, depth[i], "ACGT"[call],
	       freqs[i][0], freqs[i][1], freqs[i][2], freqs[i][3],
	       score);
	*/
    }
    /* printf("Best score is %f at %d\n", best_score, best_pos); */

    /* Extend */
    nseq = 0;
    if (best_pos > 0) {
	int furthest_left = end;

	for (i = 0; i < nr; i++) {
	    seq_t *s;
	    int right;
	    int r_pos;
	    int score;

	    if (r[i].rec == 0)
		continue;

	    s = cache_search(io, GT_Seq, r[i].rec);
	    s = cache_rw(io, s);

	    if (furthest_left > r[i].start)
		furthest_left = r[i].start;

	    /*
	     * end + best_pos is the furthest right we can go, but this
	     * specific read may not be justified in reaching that far
	     * if it has too many disagreements.
	     */
	    if ((s->len > 0) ^ r[i].comp) {
		int best_r = 0, j, k;
		int len = ABS(s->len);

		//printf(">%s\t", s->name);

		r_pos = 0;
		score = 0;
		//for (k = s->right, j = 0; j < best_pos && k < len; j++, k++) {
		for (k = end - r[i].start + 1, j = 0; j < best_pos && k < len; j++, k++) {
		    if (new_cons[j] == toupper(s->seq[k])) {
			score += match_score;
			if (best_r <= score) {
			    best_r  = score;
			    r_pos = k+1;
			}
		    } else {
			score += mismatch_score;
		    }

		    //putchar(new_cons[j] == toupper(s->seq[k])
		    //	    ? toupper(s->seq[k])
		    //        : tolower(s->seq[k]));
		}
		//putchar('\n');

		if (s->right != r_pos) {
		    s->right  = r_pos;
		    nseq++;
		}
	    } else {
		int best_r = 0, j, k;
		int len = ABS(s->len);

		//printf("<%s\t", s->name);

		r_pos = 0;
		score = 0;
		//for (k = s->left-2, j = 0; j < best_pos && k >= 0; j++, k--) {
		for (k = r[i].end - end - 1, j = 0; j < best_pos && k >= 0; j++, k--) {
		    char b = complement_base(s->seq[k]);
		    if (new_cons[j] == b) {
			score += match_score;
			if (best_r <= score) {
			    best_r  = score;
			    r_pos = k-1;
			}
		    } else {
			score += mismatch_score;
		    }

		    //putchar(new_cons[j] == toupper(b)
		    //	    ? toupper(b)
		    //	    : tolower(b));
		}
		//putchar('\n');

		if (s->left != r_pos+2) {
		    s->left  = r_pos+2;
		    nseq++;
		}
	    }
	}

	vmessage("    Extended by %d, adjusting %d sequence clip%s\n",
		 best_pos, nseq, nseq == 1 ? "" : "s");

	bin_invalidate_consensus(io, crec, furthest_left, end + best_pos);
    } else {
	vmessage("    Unable to extend contig\n");
    }
    free(r);

    cache_decr(io, c);
    cache_flush(io);
    return 0;
}

/*
 * The main contig_extend interface.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_extend(GapIO *io, tg_rec *contig, int ncontigs, int min_depth,
		  int match_score, int mismatch_score) {
    int i, err = 0;

    for (i = 0; i < ncontigs; i++) {
	/* Left end */
	UpdateTextOutput();
	complement_contig(io, contig[i]);
	err |= contig_extend_single(io, contig[i], 1, min_depth,
				    match_score, mismatch_score);

	/* Right end */
	UpdateTextOutput();
	complement_contig(io, contig[i]);
	err |= contig_extend_single(io, contig[i], 0, min_depth,
				    match_score, mismatch_score);
    }

    return err ? -1 : 0;
}

