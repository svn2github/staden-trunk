/*
 * File: oligo_sel.c
 * Version:
 *
 * Author: James Bonfield
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */

/* #define DEBUG_OLIGO_SEL */
#define DSTR stdout

#define SUBVERSION

#include <stdio.h>
#include <stdlib.h>

#include "list_proc.h"
#include "io_utils.h"
#include "qual.h"
#include "IO.h"
#include "tagUtils.h"
#include "edUtils.h"
#include "fort.h"
#include "xalloc.h"
#include "misc.h"
#include "gap_globals.h"
#include "io-reg.h"
#include "text_output.h"
#include "oligo_sel.h"
#include "dstrand.h"
#include "FtoC.h"
#include "primlib.h"

#define MAXCOMLEN (1024)

/* Static global variables */
static int_f last_gel;	/* remembers where we last were for		*/
    			/* for efficiency (defaults leftmost)		*/
static void *dl;	/* dlist - replaces old "list" mechanism	*/

/* --- hooks --- */

int insert_size(GapIO *io, int N)
{
    int insert_min;
    int insert_max;

    (void) get_read_info(io,
			 N,
			 NULL,0,
			 NULL,0,
			 NULL,0,
			 NULL,0,
			 NULL,
			 &insert_min,
			 &insert_max,
			 NULL,
			 NULL,
			 NULL,
			 NULL,
			 NULL,
			 NULL,
			 NULL);

    return insert_min;
}


int avg_read_len(GapIO *io) {
    static int avg_len = 0;

    if (!avg_len) {
	int i, end = NumReadings(io);
	int len_tot = 0;

	for (i = 1; i <= end; i++) {
	    len_tot += ABS(io_length(io, i));
	}
	
	avg_len = (int)(0.5 + len_tot/end);
    }

    return avg_len;
}

/*
 * Do the actual oligo selection.
 * Returns the new offset to search from. (This takes into account the
 * average gel reading length so that we do not finish we lots of suggested
 * oligos next to each other.)
 */
static int find_oligos(GapIO *io, int offset, int lreg, char *con,
		       int olilen, int olibak, char  *comment, char *oligo,
		       char *consensus, char sense, char *listname,
		       int *temnum, int olinum, int contig,
		       primlib_state *state)
{
    int cur_gel, oligostart, choice, oligolen, oligoend;
    int bestscore = 0, bestgel = 0, template = 0, score, olis;
    char *gelname, tmpbuf[100];
    char tempname[DB_NAMELEN + 1];
    char dbname[DB_FILELEN + 1], *cp;

    strcpy(dbname, io_name(io));
    if (cp = strrchr(dbname, '.'))
	*cp = 0;

    /*
     * Find oligo based on this position.
     * We choose a region starting OLISTART back from here, extending
     * for OLILEN bases.
     */
    olis = offset - (olilen + olibak);
    if (olis < 1)
	olis = 1;
    memcpy(consensus, &con[olis-lreg], olilen);
    consensus[olilen] = '\0';
    
    if (-1 == primlib_choose(state, consensus))
	return -1;

    if (state->nprimers <= 0)
	return avg_read_len(io);

    /*
     * Choose the first in the list (for the time being). This should be
     * the oligo with the best score.
     */

    /* choice = 0; */
    
    oligostart = 0;
    /* Choose the first 'olinum' oligos from the list */
    for (choice = 0; choice < olinum && choice < state->nprimers; choice++) {

	memset(oligo, 0, olilen);
	memcpy(oligo, &consensus[state->primers[choice].start],
	       state->primers[choice].length);

#ifdef DEBUG_OLIGO_SEL
    fprintf(DSTR, "Choosing oligo '%s'\n", oligo);
#endif

	/* find sequences that could yield a valid primer */
	cur_gel = last_gel;
	oligostart = state->primers[choice].start + offset - (olilen+olibak);
	oligolen = state->primers[choice].length;

#ifdef DEBUG_OLIGO_SEL
	fprintf(DSTR, "oligo starts %d, ends %d\n",
		oligostart, oligostart + oligolen);
#endif

	do {
	    /* only look for positive directing sequences */
	    if (io_length(io, cur_gel) >= 0) {
		/* stop when we've gone past our oligo */
		if (io_relpos(io, cur_gel) >= oligostart)
		    break;
		/*
		 * If the sequence ends past end of oligo then we have a
		 * candidate template.
		 */
		if (io_relpos(io, cur_gel) + insert_size(io, cur_gel) -
		    avg_read_len(io) >= oligostart + oligolen) {
		    /* remember possible template */
		    if (io_relpos(io, cur_gel) + io_length(io, cur_gel) >
			oligostart + oligolen)
			template = cur_gel;
#ifdef DEBUG_OLIGO_SEL
		    fprintf(DSTR, "Covered by gel no. %d (ends %d)\n",
			    cur_gel, io_relpos(io, cur_gel) + io_length(io, cur_gel));
#endif

		    /*
		     * Our ideal chosen gel should have greater than the
		     * average reading length left in this insert; should
		     * not be a short read (could signify bad clone); should
		     * have the oligo is a good quality region (ie the further
		     * away the 3' end the better).
		     *
		     * FIXME: currently we cannot implement the insert size
		     * part of this as we do not know the position of the gel
		     * in its insert.
		     */

#ifdef DEBUG_OLIGO_SEL
		    fprintf(DSTR, "insert_size(gel no. %d) = %d\n",
			    cur_gel, insert_size(io, cur_gel));
#endif		    
		    score = io_length(io, cur_gel) + io_relpos(io, cur_gel)
			+ io_length(io, cur_gel) - oligostart;

		    if (score > bestscore)
			bestgel = cur_gel, bestscore = score;
		}
	    }
	    cur_gel = io_rnbr(io, cur_gel);
	} while (cur_gel != 0);
	
	if (bestgel == 0) {
	    return -1;
	}
#ifdef DEBUG_OLIGO_SEL
	fprintf(DSTR, "Choosing gel no. %d for template\n", bestgel);
#endif

	/*
	 * Create tag - far more convoluted than it needs to be. We cannot tag
	 * the consensus so we have to tag a sequence instead. Preferably we
	 * use the template, failing that a sequence in the same direction, or
	 * failing that any sequence that at least covers the region.
	 */
	oligoend = oligostart + oligolen;
	if (bestgel && io_relpos(io, bestgel) + io_length(io, bestgel) >= oligoend)
	    template = bestgel;
	else if (template == 0) {
	    /* no template in this direction */
	    cur_gel = last_gel;
	    while (cur_gel != 0) {
		if ((io_relpos(io, cur_gel) + 
		     (io_length(io, cur_gel)>0?io_length(io, cur_gel):-io_length(io, cur_gel)
		      )) > oligoend) {
		    template = cur_gel;
		    break;
		}
		cur_gel = io_rnbr(io, cur_gel);
	    }
	}
	if (template) {
	    GReadings r;
	    GTemplates t;

	    /*
	     * In the comment we use 'bestgel' as the template even if we
	     * are tagging a different sequence.
	     */
	    gel_read(io, bestgel, r);
	    gelname = io_rname(io, bestgel);
	    if (0 == template_read(io, r.template, t)) {
		TextRead(io, t.name, tempname, DB_NAMELEN);
		Fstr2Cstr(tempname, DB_NAMELEN, tempname, DB_NAMELEN+1);
	    } else {
		strcpy(tempname, gelname);
	    }

#ifdef DEBUG_OLIGO_SEL
	    fprintf(DSTR, "Using gel no. %d for template\n", template);
#endif
	    sprintf(comment,
		    "Template=%s\nReading=%s\nName=%s.%d\nSequence=%s\n",
		    tempname, gelname, dbname, *temnum, oligo);

	    {
		GReadings r;

		gel_read(io, template, r);

		insert_NEW_tag(io, template, oligostart -
			       io_relpos(io, template) + 1 + r.start,
			       oligolen, "OLIG", comment, sense == '-');
	    }

#ifdef DEBUG_OLIGO_SEL
	    puts(comment);
#endif
	    sprintf(tmpbuf, "%s %s %s.%d %s %d %c",
		    tempname, gelname, dbname, *temnum, oligo, oligostart,
		    sense);
	    /* add_to_list(listname, tmpbuf); */
	    add_to_dlist(dl, tmpbuf);

	    vmessage("At %5d - template %s, primer %s, number %d\n",
		    (sense == '+') ? offset
		     : io_clength(io, contig) - offset +1,
		    tempname, oligo, (*temnum)++);
	    fflush(stdout);
	    UpdateTextOutput();
	} else {
	    vmessage("At %d - no suitable oligos found\n",
		    (sense == '+') ? offset
		     : io_clength(io, contig) - offset +1);
	    return -1;
	}
	UpdateTextOutput();
    }

    return oligostart + avg_read_len(io);
}


/*
 * Returns the offset in 'con' for the oligo, or -1 for none. Sets len too.
 */
static int find_oligo_for_gel(GapIO *io, int search_st, int search_en,
			      int gel, char *con, int *len, int l,
			      char *consensus, primlib_state *state) {
    GReadings r;
    GTemplates t;
    int con_st, con_len;
    
    gel_read(io, gel, r);
    template_read(io, r.template, t);
    
    /* Check the strand is correct */
    if (r.sense != GAP_SENSE_ORIGINAL && t.strands != 2)
	return -1;
    
    /* Check length/pos, must overlap with search_st-search_en */
    if (r.position + r.sequence_length <= l - search_en ||
	r.position > l - search_st)
	return -1;
    
    /*
     * Get the overlapping chunk of consensus.
     *
     * From MAX(r.position, l-search_en)
     * To   MIN(r.position+r.sequence_length, l-search_st)
     */
    con_st = (r.position < l - search_en)
	? l - search_en : r.position;
    con_len = ((r.position + r.sequence_length > l - search_st)
	       ? l - search_st : r.position + r.sequence_length) - con_st;
    if (con_len < 10)
	return -1;
    strncpy(consensus, &con[con_st], con_len);
    consensus[con_len] = 0;
    
    /* Analyse it */
    if (-1 == primlib_choose(state, consensus))
	return -1;
    if (state->nprimers <= 0)
	return -1;
    
    *len = state->primers[0].length;
    return state->primers[0].start + con_st;
}

/*
 * Find a primer suitable for walking. This is sometimes found by the single
 * stranded search, but not always.
 */
static void find_oligo_walk(GapIO *io, int search_st, int search_en,
			    int contig, char *con, char sense,
			    int *oli_start, int oli_num, char *listname,
			    primlib_state *state) {
    int g, l;
    char *consensus, *cp;
    char tmpbuf[1024], *gelname, dbname[DB_FILELEN+1];
    char tempname[DB_NAMELEN+1], comment[MAXCOMLEN];
    int max_len;

    if (NULL == (consensus = (char *)xmalloc(search_en - search_st + 1))) {
	return;
    }

    strcpy(dbname, io_name(io));
    if (cp = strrchr(dbname, '.'))
	*cp = 0;

    /*
     * Find the nearest forward facing reading at the 3' end of this contig.
     * Use this to suggest an oligo. We ought to check others too to find
     * the best, but this is a lot more work.
     */

    max_len = find_max_gel_len(io, contig, 0);
    for (g = io_crnbr(io, contig), l = io_clength(io, contig);
	 g && l - io_relpos(io, g) < max_len;
	 g = io_lnbr(io, g)) {
	int offset, len;
	GReadings r;
	GTemplates t;

	if (-1 == (offset = find_oligo_for_gel(io, search_st, search_en,
					       g, con, &len, l, consensus,
					       state)))
	    continue;

	/* Get information */
	gel_read(io, g, r);
	gelname = io_rname(io, g);
	if (0 == template_read(io, r.template, t)) {
	    TextRead(io, t.name, tempname, DB_NAMELEN);
	    Fstr2Cstr(tempname, DB_NAMELEN, tempname, DB_NAMELEN+1);
	} else {
	    strcpy(tempname, gelname);
	}

	/* Create tag */
	sprintf(comment,
		"Template=%s\nReading=%s\nName=%s.%d\nSequence=%.*s\n",
		tempname, gelname, dbname, *oli_start,
		len, &con[offset]);

	insert_NEW_tag(io, g, offset - io_relpos(io, g) + 2 + r.start,
		       len, "OLIG", comment, sense == '-');

	/* Output to primer file */
	sprintf(tmpbuf, "%s %s %s.%d %.*s %d %c",
		tempname, gelname, dbname, *oli_start,
		len, &con[offset],
		offset, sense);
	/* add_to_list(listname, tmpbuf); */
	add_to_dlist(dl, tmpbuf);

	vmessage("At %5d - template %s, primer %.*s, number %d\n",
		 (sense == '+') ? offset : io_clength(io, contig) - offset + 1,
		 tempname, len, &con[offset], (*oli_start)++);
	UpdateTextOutput();

	if (--oli_num == 0)
	    break;
    }

    xfree(consensus);
}


static void primer_top(GapIO *io, int contig, int lreg, int rreg,
		       int search_st, int search_en, int oli_num,
		       int *oli_start, char *out_list,
		       char *con, char *qual, char sense,
		       primlib_state *state) {
    char *oligo, *tmp;
    char comment[MAXCOMLEN];
    register int i, j = 0, l;

    if (NULL == (oligo = (char *)xmalloc(2 * (search_en - search_st + 1))))
	return;
    tmp = &oligo[search_en - search_st + 1];

    last_gel = io_clnbr(io, contig);

    /* Do not look at the very left hand end. */
    if ((i = lreg) < search_en + 1)
	i = search_en + 1;

    /* scan through quality buffer */
    for (; i<rreg; i++) {
	register int ind = i - lreg;

	/* strong negative strand, but no positive strand */
	if (qual[ind] == R_NONE_GOOD || qual[ind] == R_NONE_BAD) {
	    
	    /* find length of single stranded section */
	    j = ind;
	    while((j <= rreg - lreg) &&
		  (qual[j] == R_NONE_GOOD || qual[j] == R_NONE_BAD))
		j++;

#ifdef DEBUG_OLIGO_SEL
	    fprintf(DSTR, "Single strand at %d - %d, len %d\n",
		    i, j+lreg, j-ind);
#endif

	    l = find_oligos(io, i, lreg, con, search_en - search_st, search_st,
			    comment, oligo, tmp, sense, out_list, oli_start,
			    oli_num, contig, state);
#ifdef DEBUG_OLIGO_SEL
	    fprintf(DSTR, "Oligo created should double strand to %d\n", l);
#endif

	    if (l != -1)
		i = l;
	    if (i < j+lreg)
		i = j+lreg;
	}
    }

    /*
     * Now find an oligo for walking.
     *
     * Check that we haven't already dealt with this case.
     */
    if (j+lreg != rreg) {
	find_oligo_walk(io, search_st, search_en, contig, con, sense,
			oli_start, oli_num, out_list, state);
    }

    xfree(oligo);
}

/*
 * Main C interface to the suggest primers function.
 * This looks for single stranded regions and finds suitable primers to be
 * used for walking to double strand this region.
 *
 * Internally the routine works only on the forward strand. We then complement
 * and repeat to check the other strand (and finally complement back again).
 */
void suggest_primers(GapIO *io, int contig, int lreg, int rreg,
		     int search_st, int search_en, int oli_num,
		     int *oli_start, char *out_list,
		     float con_cut, float qual_cut, char *primer_defs) {
    char *con, *qual;
    primlib_state *state;
    primlib_args *args;
    int clen = io_clength(io, contig);

    state = primlib_create();
    args = primlib_str2args(primer_defs);
    if (!args) {
	verror(ERR_WARN, "suggest_primers",
	       "Failed to parse primer arguments");
	return;
    }
    primlib_set_args(state, args);
    free(args);

    /* Find maximum extents */

    if (NULL == (con = (char *)xmalloc(clen+1)))
	return;

    if (NULL == (qual = (char *)xmalloc(clen+1))) {
	xfree(con);
	return;
    }

    /* Do forward strand */
    calc_consensus(contig, 1, clen, CON_SUM, con, NULL, NULL, NULL,
		   con_cut, qual_cut, database_info, (void *)io);
    calc_quality(contig, 1, clen, qual, con_cut, qual_cut,
		 database_info, (void *)io);
    primer_top(io, contig, lreg, rreg, search_st, search_en, oli_num,
	       oli_start, out_list,
	       &con[lreg-1], &qual[lreg-1], '+', state);
    flush2t(io);

    /* Complement */
    dbl_complement(io, &lreg, &rreg, contig); /* dstrand.c */
    
    /* Do reverse strand */
    calc_consensus(contig, 1, clen, CON_SUM, con, NULL, NULL, NULL,
		   con_cut, qual_cut, database_info, io);
    calc_quality(contig, 1, clen, qual, con_cut, qual_cut,
		 database_info, (void *)io);
    primer_top(io, contig, lreg, rreg, search_st, search_en, oli_num,
	       oli_start, out_list,
	       &con[lreg-1], &qual[lreg-1], '-', state);
    
    /* Complement again */
    dbl_complement(io, &lreg, &rreg, contig);
    flush2t(io);

    xfree(con);
    xfree(qual);
    
    primlib_destroy(state);
}

void suggest_primers_single(GapIO *io, int contig, int lreg, int rreg,
			int search_st, int search_en, int oli_num,
			int *oli_start, char *out_list, char *primer_defs) {
    reg_anno ra;

    if (contig_lock_write(io, contig) == -1) {
	verror(ERR_WARN, "suggest_primers", "Contig is busy");
	return;
    }

    if (lreg == 0)
	lreg = 1;
    if (rreg == 0)
	rreg = io_clength(io, contig);

    vmessage("Selecting oligos for contig %s between %d and %d\n",
	     get_contig_name(io, contig), lreg, rreg);
    UpdateTextOutput();
    suggest_primers(io, contig, lreg, rreg, search_st, search_en,
		    oli_num, oli_start, out_list,
		    consensus_cutoff, quality_cutoff, primer_defs);
    vmessage("\n");

    /*
     * Notify of database change.
     */
    ra.job = REG_ANNO;

    contig_notify(io, contig, (reg_data *)&ra);
}


char *suggest_primers_list(GapIO *io, int num_contigs, contig_list_t *contigs,
			   int search_st, int search_en, int oli_num,
			   int oli_start, char *primer_defs) {
    int i;
    char *res;
    reg_buffer_start rs;
    reg_buffer_end re;

    /* Notify of the start of the flurry of updates */
    rs.job = REG_BUFFER_START;
    for (i = 0; i < num_contigs; i++) {
	contig_notify(io, contigs[i].contig, (reg_data *)&rs);
    }

    dl = alloc_dlist();

    for (i = 0; i<num_contigs; i++)
	suggest_primers_single(io, contigs[i].contig,
			       contigs[i].start, contigs[i].end,
			       search_st, search_en,
			       oli_num, &oli_start, "dummy", primer_defs);

    res = strdup(read_dlist(dl));
    free_dlist(dl);

    /* Notify the end of our updates */
    re.job = REG_BUFFER_END;
    for (i = 0; i < num_contigs; i++) {
	contig_notify(io, contigs[i].contig, (reg_data *)&re);
    }

    return res;
}
