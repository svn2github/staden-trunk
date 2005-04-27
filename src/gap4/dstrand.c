/*
 * File: dstrand.c
 * Version:
 *
 * Author: James Bonfield
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Created:
 * Updated:
 *
 */

/*#define DEBUG_DSTRAND 1*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "io_utils.h"
#include "qual.h"
#include "edUtils.h"
#include "align.h"
#include "fort.h"
#include "misc.h"
#include "xalloc.h"
#include "gap_globals.h"
#include "io-reg.h"
#include "text_output.h"
#include "dstrand.h"
#include "FtoC.h"
#include "IO2.h"
#include "complement.h"

/* Static global variables */
static int last_gel;	/* remembers where we last were for		*/
    			/* for efficiency (defaults leftmost)		*/

static int countdb;	/* how many bases have been double stranded	*/
static int consins;	/* how many insertions into consensus		*/
static int counthl;	/* total number of holes filled 		*/

#define MINHOLELEN (5)	/* minimum single stranded size to bother with	*/
			/* not used anymore */

#define ALEXTRA (40)	/* no. extra bases to align with past hole end	*/
#define DSTR stderr


#ifdef notdef
/*----------------------------------------------------------------------*/

#define MAXINS (3)	/* maximum number of neighbouring inserts	*/
#define INSFREQ (15)	/* lowest len/inserts frequency allowed		*/

/*
 * Evaluates how good an alignment is.
 * The integer returned is how much of the alignment to trust (the length
 * of this trusted part).
 */
static int evalal(char *seq1, char *seq2, int *albuf, int_f seqlen) {
    int continserts, freqinserts, scan;
    
    /* display the alignment as it stands */
#if DEBUG_DSTRAND
    cdisplay(seq1, seq2, seqlen, seqlen, albuf, 0, 0);
#endif

    /* end if more than MAXINS inserts or len/inserts < INSFREQ */
    continserts = 0;
    freqinserts = -1;
    for (scan = 0; scan < seqlen; scan++) {
	if (albuf[scan] != 0) {
	    continserts++;
	    freqinserts++;
	    if (continserts >= MAXINS) {
		scan-=continserts;
		break;
	    }
	    if (freqinserts && scan/freqinserts <= INSFREQ) {
		scan--;
		break;
	    }
	} else
	    continserts = 0;
    }
#if DEBUG_DSTRAND
    fprintf(DSTR, "scan:%d, conti:%d, freqi:%d\n",
	    scan, continserts, freqinserts);
#endif

    return scan;
}
#endif

#ifdef notdef
/*----------------------------------------------------------------------*/

/*
 * Evaluates how good an alignment is.
 * The integer returned is how much of the alignment to trust (the length
 * of this trusted part).
 */
static int evalal(char *seq1, char *seq2, int *albuf, int_f seqlen) {
    int op = 0, scan = 0, score = 20;
    float mismatch = 0, cutgaps = 0, congaps = 0;
    int *alptr = albuf;
    char *s1ptr = seq1, *s2ptr = seq2;
    
    /* display the alignment as it stands */
#if DEBUG_DSTRAND
    cdisplay(seq1, seq2, seqlen, seqlen, 0, albuf, 0, 0);
#endif

    /*
     * Need to add max_score - cut the score off when it gets too high so
     * we cannot accumalate good data and hence resistance to bad data.
     */
    while (scan++ < seqlen) {
	mismatch*=.95;
	cutgaps*=.95;
	congaps*=.95;
	if (op == 0 && *alptr == 0) {
	    if (*s1ptr == *s2ptr) {
#if DEBUG_DSTRAND == 2
		fprintf(DSTR, "'%c' = '%c'", *s1ptr, *s2ptr);
#endif
		score++;
	    } else {
#if DEBUG_DSTRAND == 2
		fprintf(DSTR, "'%c' ! '%c'", *s1ptr, *s2ptr);
#endif
		score-=11;
		mismatch++;
	    }
	    alptr++;
	    s1ptr++;
	    s2ptr++;
	} else {
	    if (op == 0) {
		op = *alptr++;
	    }

	    if (op > 0) {
#if DEBUG_DSTRAND == 2
		fprintf(DSTR, "' ' - '%c'", *s2ptr);
#endif
		op--;
		s2ptr++;
		cutgaps++;
		mismatch++;
		score-=(cutgaps*2);
	    } else /* op < 0 */ {
#if DEBUG_DSTRAND == 2
		fprintf(DSTR, "'%c' - ' '", *s1ptr);
#endif
		op++;
		s1ptr++;
		congaps++;
		mismatch++;
		score-=(congaps*2);
	    }
	}
#if DEBUG_DSTRAND == 2
	fprintf(DSTR," score=%4d, mismatch=%2.2f, gaps(cut)=%2.2f, gaps(con)=%2.2f, off=%d\n",
	score, mismatch, cutgaps, congaps, scan); 
#endif
	
	if (mismatch > 3.5 || cutgaps > 3 || congaps > 1)
	    break;
    }
    return scan -1;
}
#endif

/*----------------------------------------------------------------------*/

#define FILL_BONUS 10
/*
 * Evaluates how good an alignment is.
 * The integer returned is how much of the alignment to trust (the length
 * of this trusted part).
 */
static int old_evalal(char *seq1, char *seq2, align_int *albuf, int_f seqlen,
		  int_f maxmis, int_f missc, int_f matchsc, int_f padsc,
		  int_f plen) {
    int score = 0, op = 0, scan = 0, miscount = 0, bestscan = 0, bestscore = 0;
    align_int *alptr = albuf;
    char *s1ptr = seq1, *s2ptr = seq2;

    /* display the alignment as it stands */
#if DEBUG_DSTRAND
    fprintf(DSTR, "Using old evaluation scheme instead.\n");
    fprintf(DSTR, "cutlen=%d, plen=%d\n", seqlen, plen);
    cdisplay(seq1, seq2, seqlen, seqlen, 0, albuf, 0, 0);
#endif

    while (scan < seqlen && miscount <= maxmis) {
	if (s1ptr - seq1 == plen) {
	    score += FILL_BONUS;
	    /* only continue extending if perfect match now */
	    missc = -999;
	    padsc = -999;
	}

	if (score >= bestscore) {
	    bestscore = score;
	    bestscan = scan;
	}
	if (op == 0 && *alptr == 0) {
	    if (*s1ptr == *s2ptr) {
		/* correct match */
		/*
		 * Note *FAILURE* here to do reasonable things if both the
		 * aligned sequences happen to contain a '-'. It will simply
		 * assume that both the sequences agree with each other and
		 * give them a positive score. However for speeds sake it's
		 * one less thing to check and the end result is that very
		 * ocassionally we'll double strand something we should have.
		 * This will show up in the quality checks so it'll still
		 * get dealt with.
		 */
#if DEBUG_DSTRAND
		fprintf(DSTR, "'%c' = '%c'", *s1ptr, *s2ptr);
#endif
		score += matchsc;
	    } else {
		/* incorrect match */
		/*
		 * Should we count less for failing to match our new sequence
		 * against the existing one (consensus in this case) then the
		 * existing one happens to be a '-' (ie unknown). In this case
		 * if the consensus is unknown should we penalise the sequence
		 * regardless? Currently we do....
		 */
#if DEBUG_DSTRAND
		fprintf(DSTR, "'%c' ! '%c'", *s1ptr, *s2ptr);
#endif
		score += missc;
		miscount++;
	    }
	    alptr++;
	    s1ptr++;
	    s2ptr++;
	} else {
	    miscount++;
	    if (op == 0) {
		op = *alptr++;
	    }
	    if (op > 0) {
#if DEBUG_DSTRAND
		fprintf(DSTR, "' ' - '%c'", *s2ptr);
#endif
		/*
		 * Pad in sequence - don't penalise if there's also a pad
		 * in the consensus.
		 */
		if (*s1ptr != '*')
		    score += padsc;
		op--;
		s2ptr++;
	    } else {
#if DEBUG_DSTRAND
		fprintf(DSTR, "'%c' - ' '", *s1ptr);
#endif
		/* pad in consensus */
		score += padsc;
		op++;
		s1ptr++;
	    }
	}
#if DEBUG_DSTRAND
	fprintf(DSTR, " score=%4d, off=%3d, miscount=%d\n",
		score, scan, miscount);
#endif
	scan++;
    }

    return bestscan+1;
}

/*----------------------------------------------------------------------*/

/*
 * Evaluates how good an alignment is.
 * The integer returned is how much of the alignment to trust (the length
 * of this trusted part).
 * FIXME: easy to optimise this
 */
static int evalal(char *seq1, char *seq2, align_int *albuf, int_f seqlen,
		  int_f maxmis, int_f misperc, int_f plen, int *plug,
		  int *pads, int *problems) {
    int scan = 0, miscount = 0, op = 0, con_pad = 0, bestscan = 0, seq_pad = 0;
    char *s1ptr = seq1, *s2ptr = seq2;
    align_int *alptr = albuf;

#if DEBUG_DSTRAND
    fprintf(DSTR, "cutlen=%d, plen=%d\n", seqlen, plen);
    cdisplay(seq1, seq2, seqlen, seqlen, 0, albuf, 0, 0);
#endif

    scan++;
    while (scan <= seqlen && scan <= plen + con_pad) {
	if ((miscount <= maxmis || (100 * (float)miscount / scan) <= misperc))
	    bestscan = scan;
#if DEBUG_DSTRAND
	else
	    fprintf(DSTR, "Beyond criteria @ scan=%d\n", scan);
#endif
	/* save some needless work */
	if (miscount > maxmis)
	    break;

	if (op == 0 && *alptr == 0) {
	    op = *alptr++;
	    if (*s1ptr++ == *s2ptr++) {		/* correct match */
	    } else {				/* incorrect match */
		miscount++;
	    }
#if DEBUG_DSTRAND
	    fprintf(DSTR, "'%c' . '%c'", *(s1ptr-1), *(s2ptr-1));
#endif
	} else {
	    miscount++;
	    if (op == 0) {
		op = *alptr++;
	    }
	    if (op > 0) {
#if DEBUG_DSTRAND
		fprintf(DSTR, "' ' - '%c'", *s2ptr);
#endif
		seq_pad++;
		op--;
		s2ptr++;
	    } else {
#if DEBUG_DSTRAND
		fprintf(DSTR, "'%c' - ' '", *s1ptr);
#endif
		con_pad++;
		op++;
		s1ptr++;
	    }
	}
#if DEBUG_DSTRAND
	fprintf(DSTR, " scan=%3d, miscount=%d, misperc=%5.2f\n", 
		scan, miscount, 100 * (float)miscount / scan);
#endif
	scan++;
    }

/*    bestscan; */
    if (bestscan < plen + con_pad) {
	*problems = 999;
	*pads = 0;
	return old_evalal(seq1, seq2, albuf, seqlen, maxmis, -5, 1, -5, plen);
    } else {
	*plug = 1;
	*problems = miscount + seq_pad + con_pad;
	*pads = con_pad;
	return bestscan;
    }
}

/*----------------------------------------------------------------------*/

/*
 * Generates two buffers from an alignment. Each buffer represents a sequence
 * of edits to perform on the contig in an easily parsable fashion.
 */
static void dstrform(char *seq1, align_int *albuf, int_f seqlen, char *new1,
		     char *new2) {
    int_f scan = 0;
    int op = 0;
    char *ptr1 = new1, *ptr2 = new2, *seq1ptr = seq1;
    align_int *alptr = albuf;

    while (scan++ < seqlen) {
	*ptr1 = '.';
	*ptr2 = '.';
	if (op == 0 && *alptr == 0) {
	    alptr++;
	    *ptr1 = *seq1ptr++;
	} else {
	    if (op == 0)
		op = *alptr++;
	    if (op > 0) {
		*ptr1 = '*';
		op--;
	    } else {
		*ptr2 = '*';
		*ptr1 = *seq1ptr++;
		op++;
	    }
	}
	ptr1++;
	ptr2++;
    }
    *ptr1 = '\0';
    *ptr2 = '\0';
}

/*----------------------------------------------------------------------*/

/*
 * patch() will attempt to 'patch' up a single stranded section into a double
 * stranded section. The integer returned is the actual length we managed to
 * double strand.
 *
 * Yes - I know it's ugly!
 */
static int patch(
     GapIO *io,
     int_f *idev_p,	/* int - passed back to fortran                 */
     int    off,	/* current relative offset in contig		*/
     int    plen,	/* length of hole to patch			*/
     char **cons_pp,	/* the consensus				*/
     char **qual_pp,	/* the consensus quality			*/
     int   *cons_e,	/* number of bases of extension to consensus    */
     int    lreg,	/* Start of region				*/
     int    rreg,	/* End of region				*/
     int    miscnt,	/* int - maximum number of mismatches		*/
     int    misprc,	/* int - maximum percentage of mismatches	*/
     int    sense,	/* int - direction of double stranding		*/
     int    contignum,
     int_f *ngels_p,	/* int - passed back to padcon			*/
     int_f *lincon_p,	/* int - passed back to padcon			*/
     int_f *nconts_p,	/* int - passed back to padcon			*/
     int   *gel_l,	/* int - length of *gel_p (max gel len)		*/
     int   *n_ext	/* int - max extension to a sequence 		*/
     )
{
    /* Local variables */
    int sc;		/* a score for how good the alignment is 	*/
    int cur_gel;	/* record of which gel we're looking at 	*/
    int bestgeln = 0;	/* the best gel num found so far		*/
    int bestlen = 0;	/* longest 'patch' coverage found so far	*/
    int bestuse = 0;	/* longest needed gel extension found so far	*/
    align_int *albuf = (align_int *)xcalloc(*gel_l * 2, sizeof(align_int));
    align_int *salbuf= (align_int *)xcalloc(*gel_l * 2, sizeof(align_int));
			/* buffer for aligned sequences			*/
    int fail;		/* does getext() succeed?			*/
    char *cutbuf = (char *)xmalloc(*gel_l);
			/* buffer for cutoff data			*/
    int cutlen;		/* length of cutoff data			*/
    int tgelend, gelend;/* end of current gel (t = including cutoff)	*/
    int tmp, len, len2, tmp_f;
    int_f posn, one = 1;
    int plug, num_pads = 0, pads, problems, num_problems = 999;
    char *newb1, *newb2;
    char gelname[DB_NAMELEN + 1];
    char *cons_p, *qual_p;
    char *gel_p = (char *)xmalloc(*gel_l + 1);

    if (albuf == NULL || salbuf == NULL || cutbuf == NULL || gel_p == NULL)
	return plen;

    /*
     * We're going to try to 'patch' this single strand section.
     * So we need to find which gels (in the positive direction)
     * exist at the start of the 'hole'.
     */
    
    cur_gel = last_gel;
    plug = 0;
    cons_p = *cons_pp;
    qual_p = *qual_pp;
    do {
	/* positive only */
	if (io_length(io, cur_gel) >= 0) {
	    GReadings r;

	    /* are we two far past? */
	    if (io_relpos(io, cur_gel) >= off)
		break;

	    gelend = io_relpos(io, cur_gel) + io_length(io, cur_gel);

	    gel_read(io, cur_gel, r);
	    if (r.position + r.length/*unclipped*/ - r.start > off) {
		/* get cutoff data so we can compute total length */
		cutlen = *gel_l; /* allocated length of cutbuf */
		fail = (-1 == cgetext(io, cur_gel, cutbuf, &cutlen)) ? 1 : 0;
		tgelend = gelend + cutlen;
	    } else {
		fail = 1;
		tgelend = 0; /* stops compiler warnings */
	    }

	    if (!fail && tgelend >= off) {

#if DEBUG_DSTRAND
		fprintf(DSTR, "covered by gel %d\n", cur_gel);
#endif
		/*
		 * now do a quick(?) alignment check to see if better than
		 * any we've found before.
		 *
		 * we need to only align with as little as possible. That is
		 * do not align further than 10 bases past the end of the
		 * 'hole'. (10 is some arbitrary amount to account for the
		 * insertion of padding characters).
		 */
		if (tgelend >= (tmp_f = off + plen + ALEXTRA)) {
		    cutlen = tmp_f - gelend + 1;
		    tgelend = gelend + cutlen;
		}
		if (tgelend > (tmp_f = off + plen + ALEXTRA)) {
		    cutlen -= tgelend - tmp_f -1;
		    tgelend = gelend + cutlen;		    
		}
		
		/* Protect against over-running consensus buffer */
		if (gelend-lreg-1 + cutlen-1 > rreg) {
		    cutlen = rreg+1 - (gelend-lreg-1);
		}

		/*
		 * align maximum possible and 'evaluate' alignment (walking
		 * left to right along albuf[]).
		 */
		sc = calign(cutbuf, &cons_p[gelend-lreg-1], cutlen, cutlen,
			    NULL, NULL, NULL, NULL, 0, 0,
			    gopenval, gextendval, 0, 0, albuf);
		len = evalal(cutbuf, &cons_p[gelend-lreg-1], albuf, cutlen,
			     miscnt, misprc, plen + off - gelend,
			     &plug, &pads, &problems);
#if DEBUG_DSTRAND
		fprintf(DSTR, "aligned %d bases @ %d = %d, eval = %d\n",
			cutlen, gelend, sc, len);
#endif
		len2 = len + (int)gelend - (int)off - pads;
#if DEBUG_DSTRAND
		fprintf(DSTR, "overlap of %d (/%d)\n", len2, bestlen);
#endif
		/*
		 * Pick best coverage length.
		 * If equal pick gel with fewer problems
		 * If still equal, then pick last gel (likely to also be
		 * the shortest extension).
		 */
		if (len2 > bestlen ||
		    (len2 == bestlen && problems <= num_problems)) {
		    bestlen = len2;
		    bestuse = len;
		    num_pads = pads;
		    num_problems = problems;
		    bestgeln = cur_gel;
		    memcpy(salbuf, albuf, *gel_l * 2 * sizeof(align_int));
		    /* bestgels = ... */
#if DEBUG_DSTRAND
		    fprintf(DSTR, "New best gel! %d(%d)\n", bestgeln, bestlen);
#endif
		}
	    }
	    last_gel = cur_gel;
	}
	/* jump to next element in list */
	cur_gel = io_rnbr(io, cur_gel);
    } while (cur_gel != 0);

    /*
     * When we've got here, bestlen is the best overlap (upto max 10 more
     * than the hole length) and bestuse is the amount of data needed to
     * be extended for the gel (bestgeln).
     */

    bestuse--;

    /*
     * find best gel to use - only bother if the extra data would save an
     * experiment. This is when we either totally double strand the section
     * or if we extend by approx the average gel reading length.
     */
    /*
     * Currently cheating - use if over half the hole or longer than 20 bases.
     */
/*    if (bestgeln && (bestlen >= 20 || plen/bestlen < 2)) {*/
#if DEBUG_DSTRAND
    fprintf(DSTR, "bestgeln=%d, bestlen=%d, plen=%d\n",
	    bestgeln, bestlen, plen);
#endif

    if (bestgeln && (bestlen >= 20 || bestlen >= plen + num_pads)) {

#if DEBUG_DSTRAND
	fprintf(DSTR, "Shall use gel no. %d (len %d x=%d)\n",
		bestgeln, bestlen, bestuse);
#endif

	if (plug)
	    counthl++;

	if (bestgeln != cur_gel) {
	    /* get cutoff data so we can compute total length */
	    cutlen = *gel_l;
	    /* do we need to get 'cutlen' amount? why not 'bestuse'? */
	    fail = (-1 == cgetext(io, bestgeln, cutbuf, &cutlen)) ? 1 : 0;
	}

	xfree(albuf);
	/* switch to best (saved) alignment */
        gelend = io_relpos(io, bestgeln) + io_length(io, bestgeln);
	albuf=salbuf;
	salbuf=NULL;

	/* format data in a fashion that is easy to use */
	if (NULL == (newb1 = (char *)xmalloc(bestuse * 2 +1))) {
	    xfree(albuf);
	    xfree(gel_p);
	    return plen;
	}
	if (NULL == (newb2 = (char *)xmalloc(bestuse * 2 +1))) {
	    xfree(albuf);
	    xfree(newb1);
	    xfree(gel_p);
	    return plen;
	}
	
	dstrform(cutbuf, albuf, bestuse, newb1, newb2);
#if DEBUG_DSTRAND
	fprintf(DSTR, "'%s'\n'%s'\n\n", newb1, newb2);
#endif

	countdb += bestuse;
	readn_(idev_p, &bestgeln, gelname, (int_fl)(DB_NAMELEN));
	Fstr2Cstr(gelname, DB_NAMELEN, gelname, (int_fl)(DB_NAMELEN+1));

	{
	    int_f offset = io_relpos(io, bestgeln) + io_length(io, bestgeln);

	    if (sense == 1) /* negative direction */
		offset = io_clength(io, contignum) - offset + 1;

#if DEBUG_DSTRAND
	   printf("Double stranded %s by %d base%s at offset %d (was %d) %c\n",
		   gelname, bestuse, bestuse==1 ? "" : "s",
		   offset, offset - consins, "ny"[plug]);
#else
	    vmessage("Double stranded %s by %d base%s at offset %d%s\n",
		   gelname, bestuse, bestuse==1 ? "" : "s",
		   offset, plug ? " - Filled" : "");
#endif
	    UpdateTextOutput();
	}
	bestlen = bestuse; /* Yuk! bestlen changes meaning here :-( */
	for (tmp = 0; tmp<bestuse; tmp++) {
	    if (newb2[tmp] == '*') {
		/* pad in consensus */
#if DEBUG_DSTRAND
		fprintf(DSTR, "Inserting pad (c) at offset %d\n", gelend+tmp);
#endif
		posn = gelend+tmp;
		/* If we're past our rreg (or the end of the contig
		 * then stop now
		 */
		if (posn >= rreg - 1)
		    break;

		padcon_(&io_relpos(io, 1), &io_length(io, 1),
			&io_lnbr(io, 1), &io_rnbr(io, 1),
			ngels_p, nconts_p, gel_p, lincon_p, &posn, &one,
			&io_dbsize(io), idev_p, gel_l, *gel_l);
		(*n_ext)++;
		rreg++;
		gel_p = xrealloc(gel_p, *gel_l + 1 + *n_ext);

		/* Insertion to consensus, so realloc & shuffle the arrays */
		*cons_pp = (char *)xrealloc(*cons_pp,
					    io_clength(io, contignum) + 1);
		*qual_pp = (char *)xrealloc(*qual_pp,
					    io_clength(io, contignum) + 1);

		cons_p = *cons_pp;
		qual_p = *qual_pp;
		
		memmove(&cons_p[posn-lreg+1], &cons_p[posn-lreg], rreg-posn);
		memmove(&qual_p[posn-lreg+1], &qual_p[posn-lreg], rreg-posn);


		(*cons_e)++;
	    }
	    /*
	     * At the same time we're computing how much extra to leave on
	     * the cutoff data. (This occurs when we pad out the sequence
	     * and hence have a longer sequence than before.
	     */
	    if (newb1[tmp] == '*')
		bestlen--;
	}

	/* shrink the cutoff data - should ideally check the return. */
	(void)modext(io, (int)bestgeln, bestlen);

	/*
	 * Add our new end of sequence onto the existing one.
	 * This requires reading the sequence ('w'orking version), 
	 * adding onto the end, and writing it back. Similarly for
	 * 'r'elationships.
	 */
	readw_(idev_p, &bestgeln, gel_p, gel_l, *gel_l);
	strncpy(gel_p + io_length(io, bestgeln), newb1, (size_t)bestuse);

	/* jkb 1/3/94. Extend by bestlen not bestuse - bestlen includes the
	 * number of pads added to its length and so is the correct value to
	 * extend r.sequence_length by
	 */
	{
	    GReadings r;
	    
	    gel_read(io, bestgeln, r);
	    r.sequence_length += bestuse;
	    (*n_ext) += bestuse;
	    gel_write(io, bestgeln, r);
	    io_length(io, bestgeln) += bestuse;

	    if (r.position + r.sequence_length-1 > io_clength(io, contignum)) {
		GContigs c;
		contig_read(io, contignum, c);
		c.length = r.position + r.sequence_length-1;
		contig_write(io, contignum, c);
		io_clength(io, contignum) = c.length;
	    }
	}

	/* create the necessary tags for any pads in our extension. */
	for (tmp = 0; tmp<bestuse; tmp++) {
	    if (newb1[tmp] == '*') {
#if DEBUG_DSTRAND
		fprintf(DSTR, "Inserting pad (g) at offset %d)\n",
			gelend+tmp);
#endif
		posn = gelend+tmp - io_relpos(io, bestgeln)+1;
		io_insert_base(io, bestgeln, posn, '*');

		countdb++;
	    }
	}

	xfree(newb1);
	xfree(newb2);
    } else {
#if DEBUG_DSTRAND
	fprintf(DSTR, "No suitable gel.\n");
#endif
    }

    /* tidy up memory */
    xfree(albuf);
    if (salbuf)
	xfree(salbuf);
    xfree(cutbuf);

    if (gel_p)
	xfree(gel_p);

    /*
     * We could either return the length of the original single stranded
     * section, or we could return the length of the amount we managed
     * to double strand. (if any). In this case must make sure we do not
     * return 0 for a failed patch and hence get into infinite loops.
     * Taking the easy solution currently...
     */
    return plen;
}

/*----------------------------------------------------------------------*/

#ifdef notdef
/*
 * Calculate the average length of 'used' data in the gel readings.
 */
int avggellen(int_f *lngthg_p, int_f *ngels_p) {
    int_f i, len = 0;

    for (i=0; i<*ngels_p; i++)
	len += lngthg_p[i];
    
    return (int)(len / *ngels_p);
}
#endif

/*----------------------------------------------------------------------*/

/*
 * Searches for 'holes' on a specific strand.
 * The con and qual arrays hold the consensus and quality information. They
 * are the full arrays extending from 1 to contig length inclusive.
 * gel_p is a pointer to a buffer holding the maximum known length of a
 * sequence. This will be realloced as necessary on a worst-case scenario (all
 * pads get inserted into the longest reading). The size of this buffer is in
 * *gel_l.
 */
void dstrand_top(GapIO *io, int contig, int lreg, int rreg, int miscount,
		 int misperc, char **con, char **qual, int sense, int *gel_l) {
    register int i, j;
    f_int ngels, lincon, nconts;
    int cons_e;
    static int countdbt, consinst, counthlt;
    int n_ext;

    if (sense == 0)
	countdbt =  consinst = counthlt = 0;
    countdb = consins = counthl = 0;


    /* Setup our horrid Fortran parameters (for padcon) */
    ngels = NumReadings(io);
    nconts = NumContigs(io);
    lincon = io_dbsize(io) - contig;

    /* Initialise last remembered gel to the left most one for this contig */
    last_gel = io_clnbr(io, contig);

    for (i = lreg; i <= rreg; i++) {
	register int ind = i - lreg;

	/* Look for strong negative strand, but no positive strand */
	if ((*qual)[ind] == R_NONE_GOOD || (*qual)[ind] == R_NONE_BAD) {

	    /* Find length of single stranded section */
	    j = ind;
	    while ((j <= rreg - lreg) &&
		   ((*qual)[j] == R_NONE_GOOD || (*qual)[j] == R_NONE_BAD))
		j++;

	    if (j > rreg - lreg)
		j = rreg - lreg + 1;

#ifdef ndef
	    /*
	     * We have a minimum length of single strand to patch.
	     * This saves us having to do too many checks on strands
	     * with padding characters in etc.
	     */
	    if ((j-ind) < MINHOLELEN) {
		i = j-1+lreg;
		continue;
	    }
#endif

#if DEBUG_DSTRAND
	    fprintf(DSTR, "Single strand at %d - %d, (was %d - %d) len %d\n",
		    i, j+lreg, i-consins, j+lreg-consins, j+lreg-i);
#endif

	    /*
	     * Perform the 'operation' on the contig.
	     * And leap forward to next potential problem.
	     */
	    cons_e = 0;
	    /*
	     * We pass over 'j-ind+1' as the length due to an (as yet) unfound
	     * 'feature'. For some reason it appears that the actual length
	     * double stranded is not always the same as the amount we asked
	     * for. So just to make sure we don't leave any 1 base gaps we
	     * cheat a bit. Note: this bug appears to be indeterminate in that
	     * running double strand twice on the same data with the same args
	     * doesn't always do the same thing! (jkb 23/12/92)
	     */
	    n_ext = 0;
	    j = patch(io, handle_io(io), i, j-ind+2,
		      con, qual, &cons_e, lreg, rreg, miscount, misperc,
		      sense, contig, &ngels, &lincon, &nconts, gel_l, &n_ext);
	    (*gel_l) += n_ext;
	    /* move back quality buffer to ensure alignment with consensus */
#if DEBUG_DSTRAND
	    fprintf(DSTR, "Inserted %d pads into consensus\n", cons_e);
#endif
	    /*
	     * Take into account number of additional consensus entries.
	     * This involves shifting our consensus buffer along to the
	     * right by a bit, and changing our right margin (rreg) for the
	     * region.
	     */

	    /*
	     * If we've extended the consensus then we need to shift
	     * rreg along by the appropriate amount.
	     */
	    rreg += cons_e;
	    consins += cons_e;

	    /* Skip over the hole we just patched (& take account of cons_e) */
	    i += j + cons_e;
	}
    }

    countdbt += countdb;
    consinst += consins;
    counthlt += counthl;
    vmessage("%s strand :\n"
	   "\tDouble stranded %d base%s with %d insert%s into consensus\n"
	   "\tFilled %d hole%s\n",
	   sense?"Negative":"Positive",
	   countdb, countdb==1 ? "" : "s",
	   consins, consins==1 ? "" : "s",
	   counthl, counthl==1 ? "" : "s");
    if (sense)
	vmessage("Total :\n"
	       "\tDouble stranded %d bases with %d inserts\n"
	       "\tFilled %d holes\n",
	       countdbt, consinst, counthlt);
    UpdateTextOutput();
}

void dbl_complement(GapIO *io, int *lreg, int *rreg, int contig) {
    int tmp;

    /* Complement - need to juggle the active region firstly */
    tmp     = *rreg;
    *rreg   = io_clength(io, contig) - *lreg + 1;
    *lreg   = io_clength(io, contig) - tmp  + 1;

    complement_contig(io, contig);
}

/*
 * Main C interface to the double stranding function.
 * This looks for single stranded regions and attempts to double strand these
 * by comparing suitable cutoff data with the opposite strand. If a close match
 * is found then we extend as far as possible, within the limits of the
 * miscount and misperc parameters.
 *
 * Internally the routine works only on the forward strand. We then complement
 * and repeat to check the other strand (and finally complement back again).
 */

/*
 * FIXME: we broke this so that lreg and rreg no longer work unless they are
 * the entire contig! This entire file really needs a complete rewrite, along
 * with replacing the fortran code (padcon_) with some C.
 */
void double_strand(GapIO *io, int contig, int lreg, int rreg, int miscount,
		   int misperc, float con_cut, int qual_cut) {
    char *con = NULL, *qual = NULL;
    int gel_l = find_max_gel_len(io, 0, 0) + 1000; /* extend later if needed */
    int clen = io_clength(io, contig);
    int right_end, left_end;

    if (NULL == (con = (char *)xmalloc(clen+1)))
	goto end;

    if (NULL == (qual = (char *)xmalloc(clen+1)))
	goto end;

    /* Do forward strand */
    left_end = lreg == 1;
    right_end = rreg == clen;

    lreg--; /* count from 0 instead of 1 */
    rreg--;
    calc_consensus(contig, 1, clen, CON_SUM, con, NULL, NULL, NULL,
		   con_cut, qual_cut, database_info, (void *)io);
    calc_quality(contig, 1, clen, qual, con_cut, qual_cut,
		 database_info, (void *)io);
    dstrand_top(io, contig, lreg, rreg, miscount, misperc,
		&con, &qual, 0, &gel_l);
    lreg++;
    rreg++;

    flush2t(io);

    /* Complement contig and our region pointers within it */
    dbl_complement(io, &lreg, &rreg, contig);

    /* Just incase we've managed to stretch our contig */
    clen = io_clength(io, contig);
    if (left_end) lreg = 1;
    if (right_end) rreg = clen;

    if (NULL == (qual = (char *)xrealloc(qual, clen+1)))
	goto end;

    if (NULL == (con = (char *)xrealloc(con, clen+1)))
	goto end;

    /* Do reverse strand */
    lreg--; /* count from 0 instead of 1 */
    rreg--;
    calc_consensus(contig, 1, clen, CON_SUM, con, NULL, NULL, NULL,
		   con_cut, qual_cut, database_info, (void *)io);
    calc_quality(contig, 1, clen, qual, con_cut, qual_cut,
		 database_info, (void *)io);
    dstrand_top(io, contig, lreg, rreg, miscount, misperc,
		&con, &qual, 1, &gel_l);
    lreg++;
    rreg++;

    /* Complement again */
    dbl_complement(io, &lreg, &rreg, contig);

    flush2t(io);

 end:
    if (con)
	xfree(con);
    if (qual)
	xfree(qual);
}

void double_strand_single(GapIO *io, int contig, int lreg, int rreg,
			  int miscount, int misperc) {
    reg_length rl;

    if (contig_lock_write(io, contig) == -1) {
	verror(ERR_WARN, "double_strand", "Contig is busy");
	return;
    }

    if (lreg == 0)
	lreg = 1;
    if (rreg == 0)
	rreg = io_clength(io, contig);

    vmessage("Double stranding contig %s between %d and %d\n",
	     get_contig_name(io, contig), lreg, rreg);
    double_strand(io, contig, lreg, rreg, miscount, misperc,
		  consensus_cutoff, quality_cutoff);
    vmessage("\n");

    /*
     * Notify of database change.
     */
    rl.job = REG_LENGTH;
    rl.length = io_clength(io, contig);

    contig_notify(io, contig, (reg_data *)&rl);
}

void double_strand_list(GapIO *io, int num_contigs,
			contig_list_t *contigs,
			int miscount, float misperc) {
    int i;
    reg_buffer_start rs;
    reg_buffer_end re;

    /* Notify of the start of the flurry of updates */
    rs.job = REG_BUFFER_START;
    for (i = 0; i < num_contigs; i++) {
	contig_notify(io, contigs[i].contig, (reg_data *)&rs);
    }

    for (i = 0; i<num_contigs; i++)
	double_strand_single(io, contigs[i].contig,
			     contigs[i].start, contigs[i].end,
			     miscount, (int)misperc);

    /* Notify the end of our updates */
    re.job = REG_BUFFER_END;
    for (i = 0; i < num_contigs; i++) {
	contig_notify(io, contigs[i].contig, (reg_data *)&re);
    }
}
