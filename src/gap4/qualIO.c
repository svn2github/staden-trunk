/*
 * This file contains functions that interact with the consensus algorithm in
 * a Gap4 specific manner. These used to be in qual.c, but have been moved
 * here to allow qual.c to compile and be used without linking in the rest of
 * Gap4.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gap_globals.h"
#include "IO.h"
#include "IO2.h"
#include "qual.h"
#include "qualP.h"
#include "misc.h"
#include "xalloc.h"

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
	    char *g_seq  = NULL;
	    int1 *g_conf = NULL;
	    int g_start, g_end, g_len;
	    
	    if (0 != io_aread_seq(io, gel_seq->gel, &g_len, &g_start, &g_end,
				  &g_seq, &g_conf, /*g_opos*/ NULL, 0)) {
		if (g_seq)  xfree(g_seq);
		if (g_conf) xfree(g_conf);
		return -1;
	    }

	    gel_seq->gel_start  = g_start;
	    gel_seq->gel_end    = g_end;
	    gel_seq->gel_seq    = g_seq;
	    gel_seq->gel_conf   = g_conf;
	    gel_seq->gel_opos   = /*g_opos*/NULL;
	    gel_seq->gel_length = g_len;
	    
	    return 0;
	}
    case DEL_SEQ:
	{
	    gel_seq_t *gel_seq = &theirdata->gel_seq;

	    if (gel_seq->gel_seq ) xfree(gel_seq->gel_seq);
	    if (gel_seq->gel_conf) xfree(gel_seq->gel_conf);
	    /*if (gel_seq->gel_opos) xfree(gel_seq->gel_opos);*/
	    
	    return 0;
	}
    case GET_CONTIG_INFO:
	{
	    contig_info_t *contig_info = &theirdata->contig_info;
	    GContigs c;

	    GT_Read(io, arr(GCardinal, io->contigs, contig_info->contig-1),
		    &c, sizeof(c), GT_Contigs);

	    contig_info->length  = c.length;
	    contig_info->leftgel = c.left;
	    
	    return 0;
	}
    case DEL_CONTIG_INFO:
	{
	    return 0;
	}
    case GET_GEL_INFO:
	{
	    gel_info_t *gel_info = &theirdata->gel_info;
	    GReadings r;

	    gel_read(io, gel_info->gel, r);

	    gel_info->complemented = r.sense;
	    gel_info->position     = r.position;
	    gel_info->next_right   = r.right;
	    gel_info->as_double	   = r.chemistry & GAP_CHEM_TERMINATOR;
	    gel_info->start	   = r.start;
	    gel_info->length       = r.end - r.start - 1;
	    gel_info->unclipped_len= r.length;
	    /*
	    gel_info->as_double	   =
		((r.chemistry & GAP_CHEM_TERMINATOR) == GAP_CHEM_TERMINATOR) &&
		((r.chemistry & GAP_CHEM_TYPE_MASK) == GAP_CHEM_TYPE_BIGDYE);
	    */
	    
	    return 0;
	}
    case DEL_GEL_INFO:
	{
	    return 0;
	}
    case GET_GEL_LEN:
	{
	    return find_max_gel_len(io, 0, 0);
	}

    case SEQ_INS:
	{
	    /* Horribly inefficient at the moment! */
	    seq_ins_t *seq_ins = &theirdata->seq_ins;
	    int i;
	    int pos = seq_ins->position;
	    char *base = seq_ins->bases;
	    
	    for (i = 0; i < seq_ins->length; i++) {
		io_insert_base(io, seq_ins->gel, pos++, *base++);
	    }
	    return 0;
	}

    case SEQ_DEL:
	{
	    /* Horribly inefficient at the moment! */
	    seq_del_t *seq_del = &theirdata->seq_del;
	    int i;
	    int pos = seq_del->position;
	    
	    for (i = 0; i < seq_del->length; i++) {
		io_delete_base(io, seq_del->gel, pos);
	    }
	    return 0;
	}

    case CONS_INS:
	{
	    cons_ins_t *cons_ins = &theirdata->cons_ins;
	    /* Only allows insertion of pads */
	    printf("PADCON contig %d at %d+%d\n",
		   cons_ins->contig, cons_ins->position,
		   cons_ins->length);
	    pad_consensus(io, cons_ins->contig, cons_ins->position,
			  cons_ins->length);
	    return 0;
	}

    case IF_FLUSH:
	{
	    flush2t(io);
	    return 0;
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
int *count_confidence(GapIO *io, int contig, int start, int end)
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
	if (qual[i] > 100) qual[i] = 100;
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
