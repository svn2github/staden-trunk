#include <tcl.h>

#include "IO.h"
#include "gap_cli_arg.h"
#include "primlib.h"
#include "text_output.h"
#include "io_utils.h"
#include "qual.h"
#include "misc.h"
#include "gap_globals.h"
#include "dna_utils.h"
#include "dstring.h"
#include "finish.h"
#include "finish_pcr.h"
#include "finish_filter.h"
#include "finish_utils.h"

/* ------------------------------------------------------------------ */

/*
 * Filters primers based on strong matches to secondary binding sites.
 * The data to screen against for repeats has already been setup in the
 * fin record. Ie fin->all_cons_h and fin->external_seq_h. For this to work
 * we need to know which strand this primer was chosen for so the self-match
 * (at 100% obviously) can be ignored.
 *
 * Returns 1 for false match.
 *         0 for no match.
 */
int filter_primers(finish_t *fin, int strand, char *seq) {
    double second;

    second = secondary_primer_match(fin,
				    -1, 0, 0, /* all of all contigs */
				    1, strand, 1, seq);

    return (second >= fin->opts.pwalk_max_match) ? 1 : 0;
}


/*
 * Main picking function.
 * Picks primers from the right end of contig1 and the left end of contig2
 * suitable for a PCR reaction.
 * After calling the returned value is a an array of structures linking into
 * the primer3 primer_pair structure along with gap4 sanitised copies holding
 * padded position and length in each contig and the depadded sequence.
 * The number of elements in this array can be fetched from pstate->npairs.
 *
 * Returns g4_primer_pair array pointer for success,
 *        NULL for failure
 */
static g4_primer_pair *pick_pcr_primers2(finish_t *fin, primlib_state *pstate,
					 int contig1, int contig2)
{
    char *cons1 = NULL, *cons2 = NULL, *cons_joined = NULL;
    int pos1l, pos1r, pos2l, pos2r;
    int len1, len2;
    int *depad1 = NULL, *depad2 = NULL;
    char *upcons1 = NULL, *upcons2 = NULL;
    g4_primer_pair *pp = NULL;
    int i, j;

    /* Compute contig ranges */
    pos1l = MAX(1, io_clength(fin->io, contig1) - (fin->opts.pcr_offset1-1));
    pos1r = MAX(1, io_clength(fin->io, contig1) - (fin->opts.pcr_offset2-1));
    len1 = pos1r - pos1l + 1;
    if (len1 < 25)
	return NULL;

    pos2l = MIN(io_clength(fin->io, contig2), fin->opts.pcr_offset2);
    pos2r = MIN(io_clength(fin->io, contig2), fin->opts.pcr_offset1);
    len2 = pos2r - pos2l + 1;
    if (len2 < 25)
	return NULL;


    /* Get the depadded consensus */
    cons1 = (char *)xmalloc(len1+1);
    cons2 = (char *)xmalloc(len2+1);
    if (!cons1 || !cons2)
	goto error;

    calc_consensus(contig1, pos1l, pos1r, CON_SUM, cons1, NULL,
		   NULL, NULL, consensus_cutoff, quality_cutoff,
		   database_info, (void *)fin->io);
    calc_consensus(contig2, pos2l, pos2r, CON_SUM, cons2, NULL,
		   NULL, NULL, consensus_cutoff, quality_cutoff,
		   database_info, (void *)fin->io);
    cons1[pos1r-pos1l+1] = 0;
    cons2[pos2r-pos2l+1] = 0;

    upcons1 = strdup(cons1);
    upcons2 = strdup(cons2);

    if (!(depad1 = (int *)xmalloc((len1+1)*sizeof(int))))
	goto error;
    if (!(depad2 = (int *)xmalloc((len2+1)*sizeof(int))))
	goto error;
    depad_seq(cons1, &len1, depad1);
    depad_seq(cons2, &len2, depad2);

    /* Filter low complexity data from the consensus */
    finish_filter(fin, cons1, len1);
    finish_filter(fin, cons2, len2);
    

    /*
     * For primer3 we join our two sequences together thus, with 20 Ns:
     *
     * <CONS1>NNNNNNNNNNNNNNNNNNNN<CONS2>
     *        ^                  ^
     *        x                  y
     *
     * Points x and y define the target ("TARGET=x,y-x" in the normal
     * boulder-io input file).
     * We also need to redefine the product size range to be 20 to
     * 20 + 2*(len2-len1) allowing for full flexibility of primer
     * positioning within the two consensus fragments.
     *
     * PCR primers will then be chosen using one primer within <CONS1>
     * and the other within <CONS2>.
     */
    if (NULL == (cons_joined = (char *)xmalloc(2*(len1 + len2 +2)+20)))
	goto error;
    sprintf(cons_joined, "%sNNNNNNNNNNNNNNNNNNNN%s", cons1, cons2);

    {
	size_t l = strlen(cons_joined);
	for (i = 0; i < l; i++)
	    if (cons_joined[i] != 'A' &&
		cons_joined[i] != 'C' &&
		cons_joined[i] != 'G' &&
		cons_joined[i] != 'T')
		cons_joined[i] = 'N';
    }
    puts(cons_joined);
    printf("target = %"PRId64",%d\n", (uint64_t)strlen(cons1)+1, 20);

    /* Tweak arguments */
    pstate->p3args.primer_task = pick_pcr_primers;
    pstate->p3args.num_return = 20;

    /* Pick the primer pairs */
    if (-1 == primlib_choose_pcr(pstate, cons_joined, strlen(cons1)+1, 20))
	goto error;

    if (!(pp = (g4_primer_pair *)xmalloc(pstate->npairs * sizeof(*pp))))
	goto error;



    /* Store the primer pairs in the return structures */
    for (i = j = 0; i < pstate->npairs; i++) {
	int p1, p2, len;

	/*
	 * Only pick pairs that have not had one or both primers rejected
	 * by the secondary primer-site detection code.
	 */
	/*
	if (pstate->pairs[i].left->excl ||
	    pstate->pairs[i].right->excl) {
	    continue;
	}
	*/

	pp[j].pair = &pstate->pairs[i];

	/* Compute padded start + length for these primers. */
	p1 = depad1[pstate->pairs[i].left->start];
	p2 = depad1[pstate->pairs[i].left->start +
		    pstate->pairs[i].left->length-1];
	pp[j].contig[0] = contig1;
	pp[j].pos[0] = pos1l + p1;
	pp[j].len[0] = p2-p1+1;
	
	p1 = depad2[pstate->pairs[i].right->start-
		    pstate->pairs[i].right->length+1 - len1 -20];
	p2 = depad2[pstate->pairs[i].right->start - len1 -20];
	pp[j].contig[1] = contig2;
	pp[j].pos[1] = pos2l + p1;
	pp[j].len[1] = p2-p1+1;

	/* Copy over depadded primer sequence */
	len = MIN(pstate->pairs[i].left->length, MAX_PRIMER_LEN);
	strncpy(pp[j].seq[0], &cons_joined[pstate->pairs[i].left->start], len);
	pp[j].seq[0][len] = '\0';

	len = MIN(pstate->pairs[i].right->length, MAX_PRIMER_LEN);
	strncpy(pp[j].seq[1],
		&cons_joined[pstate->pairs[i].right->start-
			     pstate->pairs[i].right->length+1],
		len);
	pp[j].seq[1][len] = '\0';
	complement_seq(pp[j].seq[1], len);
	
	/*
	 * Check if left/right primers have secondary binding sites, caching
	 * the result (in primer_rec.excl) to avoid subsequent searches.
	 */
	if (pstate->pairs[i].left->excl == 0) {
	    if (filter_primers(fin, 0, pp[j].seq[0]))
		pstate->pairs[i].left->excl = 1;
	    else
		pstate->pairs[i].left->excl = -1;
	}

	if (pstate->pairs[i].right->excl == 0) {
	    if (filter_primers(fin, 1, pp[j].seq[1]))
		pstate->pairs[i].right->excl = 1;
	    else
		pstate->pairs[i].right->excl = -1;
	}

	/* Use only if both L & R have no 2ndary match */
	if (pstate->pairs[i].left->excl == -1 &&
	    pstate->pairs[i].right->excl == -1)
	    j++;
    }

    pstate->npairs = j;
    if (!pstate->npairs) {
	xfree(pp);
	pp = NULL;
    }


    xfree(cons1);
    xfree(cons2);
    xfree(upcons1);
    xfree(upcons2);
    xfree(cons_joined);
    xfree(depad1);
    xfree(depad2);

    return pp;
    
 error:
    if (cons1)
	xfree(cons1);
    if (cons2)
	xfree(cons2);
    if (upcons1)
	xfree(upcons1);
    if (upcons2)
	xfree(upcons2);
    if (cons_joined)
	xfree(cons_joined);
    if (depad1)
	xfree(depad1);
    if (depad2)
	xfree(depad2);
    if (pp)
	xfree(pp);

    return NULL;
}


/* Debug output */
void pcr_list_primers(primlib_state *pstate, g4_primer_pair *pp) {
    int i;
    for (i = 0; i < pstate->npairs; i++) {
	printf("pair %d: qual %f, cmpl %f, difftm %f, prodtm %f pdtm %f\n",
	       i, 
	       pp[i].pair->pair_quality,
	       pp[i].pair->compl_measure,
	       pp[i].pair->diff_tm,
	       pp[i].pair->product_tm,
	       pp[i].pair->product_tm_oligo_tm_diff);

	printf("pair %d: left pos %3d+%d/%d+%d, tm %f, gc %f %s\n",
	       i,
	       pp[i].pair->left->start,
	       pp[i].pair->left->length,
	       pp[i].pos[0], pp[i].len[0],
	       pp[i].pair->left->temp,
	       pp[i].pair->left->gc_content,
	       pp[i].seq[0]);

	printf("pair %d: right pos %d+%d/%d+%d, tm %f, gc %f %s\n",
	       i,
	       pp[i].pair->right->start,
	       pp[i].pair->right->length,
	       pp[i].pos[1], pp[i].len[1],
	       pp[i].pair->right->temp,
	       pp[i].pair->right->gc_content,
	       pp[i].seq[1]);

	printf("\n");
    }
}


/*
 * Converts a g4_primer_pair array into a dstring_t dynamic string suitable
 * for passing back to Tcl.
 */
dstring_t *g4_pp2dstring(dstring_t *ds, g4_primer_pair *pp, int npairs)
{
    int i;

    if (!ds)
	ds = dstring_create(NULL);

    for (i = 0; i < npairs; i++) {
	dstring_append(ds, "{ ");
	/* pair info */
	dstring_appendf(ds, "{ %d %f %f %f %f %f } ",
			i == 0 ? 1 : 0, /* First is always picked */
			pp[i].pair->pair_quality,
			pp[i].pair->compl_measure,
			pp[i].pair->diff_tm,
			pp[i].pair->product_tm,
			pp[i].pair->product_tm_oligo_tm_diff);
	
	/* left primer */
	dstring_appendf(ds, "{ %s %d %d %d %f %f } ",
			pp[i].seq[0],
			pp[i].contig[0], pp[i].pos[0], pp[i].len[0],
			pp[i].pair->left->temp,
			pp[i].pair->left->gc_content);

	/* right primer */
	dstring_appendf(ds, "{ %s %d %d %d %f %f } ",
			pp[i].seq[1],
			pp[i].contig[1], pp[i].pos[1], pp[i].len[1],
			pp[i].pair->right->temp,
			pp[i].pair->right->gc_content);
	dstring_append(ds, "} ");
    }

    return ds;
}

/* ------------------------------------------------------------------ */
/* Tcl interface */


dstring_t *finish_pcr_primers(finish_t *fin, char *pdefs,
			      contig_list_t *contigs, int ncontigs) {
    primlib_state *pstate;
    primlib_args *pargs;
    dstring_t *ds = NULL;
    g4_primer_pair *pp;
    int i;

    /* Create and configure the primlib_state */
    pstate = primlib_create();
    if (NULL == (pargs = primlib_str2args(pdefs)))
	return NULL;
    primlib_set_args(pstate, pargs);
    free(pargs);

    /* Linear mode, rather than combinatorial */
    for (i = 0; i < ncontigs-1; i++) {
	pp = pick_pcr_primers2(fin, pstate,
			       contigs[i].contig,
			       contigs[i+1].contig);

	if (pp) {
	    ds = g4_pp2dstring(ds, pp, pstate->npairs);
	    pcr_list_primers(pstate, pp);
	    xfree(pp);
	} else {
	    if (!ds)
		ds = dstring_create(NULL);

	    dstring_appendf(ds, "{ { 0 0 0 0 0 0 } "
			    "{ none %d 0 0 0 0 } "
			    "{ none %d 0 0 0 0 } } ",
			    contigs[i].contig,
			    contigs[i+1].contig);
	}
    }

    primlib_destroy(pstate);

    return ds;
}
    
