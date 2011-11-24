/*
 * This file contains functions that interact with the consensus algorithm in
 * a Gap4 specific manner. These used to be in qual.c, but have been moved
 * here to allow qual.c to compile and be used without linking in the rest of
 * Gap4.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "gap_globals.h"
#include "tg_gio.h"
#include "qual.h"
#include "qualP.h"
#include "misc.h"
#include "xalloc.h"
#include "gap4_compat.h"

/*
 * `Info' function. Obtains various information from either the database
 * or from the contig editor structures. For this reason a pointer to the
 * relevant routine is supplied by the calling code.
 * This is the database version.
 */
int database_info(int job, void *mydata, info_arg_t *theirdata) {
    GapIO *io = (GapIO *)mydata;

    if (io == NULL)
	return -1;

    switch (job) {
    case GET_SEQ:
	{
	    gel_seq_t *gel_seq = &theirdata->gel_seq;
	    seq_t *s;
	    int free_it = 0;

	    if (NULL == (s = (seq_t *)cache_search(io, GT_Seq, gel_seq->gel)))
		return -1;

	    if (s->len < 0) {
		s = dup_seq(s);
		complement_seq_t(s);
		free_it = 1;
	    }

	    gel_seq->gel_start  = s->left-1;
	    gel_seq->gel_end    = s->right+1;
	    gel_seq->gel_conf   = (unsigned char *)s->conf;
	    gel_seq->gel_opos   = NULL;
	    gel_seq->gel_length = ABS(s->len);

	    gel_seq->gel_seq    = malloc(gel_seq->gel_length+1);
	    memcpy(gel_seq->gel_seq, s->seq, gel_seq->gel_length);
	    gel_seq->gel_seq[gel_seq->gel_length] = 0;

	    gel_seq->gel_conf   = malloc(gel_seq->gel_length);
	    memcpy(gel_seq->gel_conf, s->conf, gel_seq->gel_length);
	    
	    if (free_it)
		free(s);

	    return 0;
	}
    case DEL_SEQ:
	{
	    gel_seq_t *gel_seq = &theirdata->gel_seq;

	    free(gel_seq->gel_seq);
	    free(gel_seq->gel_conf);

	    return 0;
	}
    case GET_CONTIG_INFO:
	{
	    contig_info_t *contig_info = &theirdata->contig_info;
	    contig_t *c;
	    contig_iterator *ci;
	    rangec_t *r;
	    
	    c = (contig_t *)cache_search(io, GT_Contig, contig_info->contig);
	    ci = contig_iter_new(io, contig_info->contig, 1, CITER_FIRST,
				 contig_info->range_start,
				 contig_info->range_end);

	    contig_info->length   = c->end - c->start + 1;
	    contig_info->iterator = ci;
	    r = contig_iter_next(io, ci);
	    contig_info->gel      = r ? r->rec : 0;
	    
	    return 0;
	}
    case CONTIG_INFO_NEXT:
	{
	    rangec_t *r;
	    contig_info_t *contig_info = &theirdata->contig_info;
	    r = contig_iter_next(io, contig_info->iterator);
	    contig_info->gel = r ? r->rec : 0;

	    return 0;
	}
    case DEL_CONTIG_INFO:
	{
	    contig_info_t *contig_info = &theirdata->contig_info;

	    contig_iter_del(contig_info->iterator);

	    return 0;
	}
    case GET_GEL_INFO:
	{
	    gel_info_t *gel_info = &theirdata->gel_info;
	    seq_t *s;
	    tg_rec cnum;
	    int pos;
	    
	    if (NULL == (s = (seq_t *)cache_search(io, GT_Seq, gel_info->gel)))
		return -1;
	    if (-1 == sequence_get_position(io, gel_info->gel, &cnum, &pos, NULL, NULL))
		verror(ERR_FATAL, "database_info",
		       "Cannot find bin for sequence %"PRIrec, gel_info->gel);

	    gel_info->complemented = s->len < 0 ? 1 : 0;
	    gel_info->position     = pos;
	    gel_info->as_double	   = 0;
	    gel_info->start	   = s->len < 0
		                         ? -s->len - s->right
		                         : s->left - 1;
	    gel_info->length       = s->right - s->left + 1;
	    gel_info->unclipped_len= ABS(s->len);
	    gel_info->template     = 0;

	    if (s->len >= 0)
		gel_info->position += s->left - 1;
	    else
		gel_info->position += -s->len - s->right;

	    return 0;
	}
    case DEL_GEL_INFO:
	{
	    return 0;
	}
    case GET_GEL_LEN:
	{
	    puts("FIXME: GET_GEL_LEN");
	    return 666666;
	}

    default:
	verror(ERR_FATAL, "database_info", "Unknown job number (%d)", job);
	return -1;
    }
}

/*
 * Count the frequence of each confidence value, returning an array of
 * 101 integers for frequences of values 0 to 100 inclusive.
 */
int *count_confidence(GapIO *io, tg_rec contig, int start, int end)
{
    char *con;
    float *qual;
    int i;
    static int freqs[101];

    for (i = 0; i <= 100; i++)
	freqs[i] = 0;

    qual = (float *)xmalloc((end - start + 1) * sizeof(*qual));
    con = (char *)xmalloc((end - start + 1) * sizeof(*con));
    if (!qual || !con)
	return NULL;
    
    /* Get the consensus along with the confidence values */
    calc_consensus(contig, start, end, CON_SUM, con, NULL, qual, NULL,
		   consensus_cutoff, quality_cutoff,
		   database_info, (void *)io);

    for (i = 0; i < end - start + 1; i++) {
	if (qual[i] < 0) qual[i] = 0;
	if (qual[i] >= 100) qual[i] = 99;
	freqs[(int)(qual[i] + 0.499)]++;
    }

    xfree(qual);
    xfree(con);

    return freqs;
}

/*
 * List information about confidence values held in freqs.
 */
int list_confidence(int *freqs, int length)
{
    double cum_errs, err, err_rate;
    int cum_freq, i;
    double total_errs;
    char num[100];

    cum_errs = 0.0;
    cum_freq = 0;
    total_errs = 0.0;
    for (i = 0; i <= 99; i++) {
	err = freqs[i] * pow(10.0, -i/10.0);
	total_errs += err;
    }

    vmessage("Sequence length = %d bases.\n"
	     "Expected errors = %7.2f bases (1/%d error rate).\n",
	     length, total_errs,
	     total_errs ? (int)(length / total_errs) : 0);
    vmessage("Value	Frequencies	Expected  Cumulative	"
	     "Cumulative	Cumulative\n");
    vmessage("			errors    frequencies	errors	"
	     "	error rate\n");
    vmessage("---------------------------------------------------"
	     "-----------------------\n");

    for (i = 0; i <= 99; i++) {
	cum_freq += freqs[i];
	err = freqs[i] * pow(10.0, -i/10.0);
	cum_errs += err;
	if (total_errs - cum_errs > 0)
	    err_rate = (length) / (total_errs - cum_errs);
	else
	    err_rate = 0;
	if (err_rate) {
	    sprintf(num, "%g", err_rate);
	} else {
	    strcpy(num, "-");
	}
	vmessage("%3d\t%6d\t\t%7.2f     %5d\t%7.2f\t\t1/%s\n",
		 i, freqs[i], err, cum_freq, cum_errs, num);
    }
    vmessage("\n");
    
    return 0;
}



/*
 * Adds to match_freqs and mismatch_freqs (indexed by confidence from 0 to 255
 * inclusive) based on the number of bases which agree or disagree with the
 * consensus.
 *
 * The buffers are not initialised to zero; this is expected to be done by the
 * caller (if desired).
 *
 * Returns 0 on success
 *        -1 on failure
 */
int get_base_confidences(GapIO *io, tg_rec contig,
			 int *match_freqs, int *mismatch_freqs) {
    char *con;
    contig_iterator *ci;
    rangec_t *r;
    int clen = io_clength(io, contig);

    /* Get the consensus along with the confidence values */
    con = (char *)xmalloc((clen+1) * sizeof(*con));
    if (!con)
	return -1;

    calc_consensus(contig, 1, clen, CON_SUM,
		   con, NULL, NULL, NULL,
		   consensus_cutoff, quality_cutoff,
		   database_info, (void *)io);

    /* Loop through all sequences in this contig */
    ci = contig_iter_new(io, contig, 1, CITER_FIRST, CITER_CSTART, CITER_CEND);

    while ((r = contig_iter_next(io, ci))) {
	seq_t *s = cache_search(io, GT_Seq, r->rec);
	seq_t *origs = s;
	int i, p;

	if ((s->len < 0) ^ r->comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	for (i = s->left-1, p = r->start + i - 1; i < s->right; i++, p++) {
	    char con_base = p-1 >= 0 && p-1 < clen ? con[p] : 'N';

	    /* Skip pads for now */
	    if (con_base == '*' || s->seq[i] == '*')
		continue;

	    if (tolower(s->seq[i]) == tolower(con_base))
		match_freqs[s->conf[i]]++;
	    else
		mismatch_freqs[s->conf[i]]++;
	}

	if (s != origs)
	    free(s);
    }

    xfree(con);

    return 0;
}


/*
 * Produces a textual report of the match and mismatch frequencies
 * Returns the total "problem" score.
 */
double list_base_confidence(int *matfreqs, int *misfreqs)
{
    int max = 0, i;
    double var_denominator = 0, var_numerator = 0;

    /*
     * Weighted variance, skipping confidences outside of valid ranges.
     * There's very little real probabilistic background to any of this,
     * so I do not want to actually state it as a "variance" figure in the
     * output. It is simply a convenient single figure which has the property
     * of increasing as more discrepant bases appear.
     */
    for (i = 3; i < 100; i++) {
	int t = matfreqs[i] + misfreqs[i];
	double expected = t * pow(10.0, -i/10.0);
	double ratio;

	if (!t)
	    continue;

	/* Plus 1 to correct for small sample sizes giving excessive ratios */
	ratio =  misfreqs[i] > expected
	    ? (misfreqs[i]+1) / (expected+1)
	    : (expected+1) / (misfreqs[i]+1);

	var_denominator += t;
	var_numerator += t * (ratio - 1) * (ratio - 1);
    }
    vmessage("Total bases considered : %d\n", (int)var_denominator);
    vmessage("Problem score          : %f\n",
	     var_numerator / var_denominator);
    vmessage("\n");

    /* Headings */
    vmessage("Conf.        Match        Mismatch           Expected      Over-\n");
    vmessage("value         freq            freq               freq  representation\n");
    vmessage("---------------------------------------------------------------------\n");

    /* Find the highest used confidence value */
    for (i = 0; i < 256; i++) {
	if (matfreqs[i] || misfreqs[i])
	    max = i;
    }

    /* Dump out the frequencies */
    for (i = 0; i <= max; i++) {
	int tot = matfreqs[i] + misfreqs[i];
	double expected = tot * pow(10.0, -i/10.0);

	vmessage("%3d\t%10d\t%10d\t%13.2f\t%7.2f\n",
		 i, matfreqs[i], misfreqs[i],
		 expected, expected ? (double)misfreqs[i] / expected : 0);
    }

    return var_numerator / var_denominator;
}
