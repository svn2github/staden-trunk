/*
 * File: searchUtils.c
 * Version:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: Search functions for the contig editor
 *
 * This file is split into two parts, probably implying that it should be
 * divided into two files. The first half consists of all the search
 * engines, the second half the user interface code.
 *
 * Observations:
 *   Only one search window can be up at one time.
 *
 * Created:
 * Updated:
 *  7-Nov-1991 SD  Select the tag when searching on tag.
 * 20-Feb-1992 SD  Renamed "OK" button to "Quit"
 * 29-Apr-1992 SD  Changes relevant to general speedup in edUtils.c
 * 11-May-1992 SD  Search by name now dependant on cursor position
 * 15-May-1992 SD  now NOSTRSTR and NOSTRDUP
 * 20-Jul-1993 JKB Added system V regex compatibility
 * 01-Mar-1994 JKB Merge in confidence value code
 * 31-Mar-1999 JohnT Added interp to regexp macros to alow tcl regexps
 */


#define REDISPLAY(X) showCursor((X),(X)->cursorSeq, (X)->cursorPos);\
    		     setDisplayPos((X),positionInContig((X),\
							(X)->cursorSeq,\
							(X)->cursorPos)); \
                     repositionTraces((X));


/*
 * The first half
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "qual.h"
#include "edUtils.h"
#include "tagUtils.h"
#include "select.h"
#include "string.h"
#include "ctype.h"
#include "misc.h"
#include "reg_exp.h"
#include "contigEditor.h"
#include "io_utils.h"
#include "xalloc.h"
#include "active_tags.h"
#include "text_output.h"
#include "tman_interface.h"
#include "search_utils.h"
#include "dna_utils.h"

/* Is C between A and B */
#define in_interval(A,B,C) ( ((A)<(B))?((A)<=(C) && (C)<=(B)):((B)<=(C) && (C)<=(A)) )

#define SEARCH_CHUNKS MAX_DISPLAY_WIDTH

static int findGelByNumber (EdStruct *xx, char *s)
/*
 * Position cursor on left end of gel sequence matched by string s.
 */
{
    int i;
    
    int gel;
    gel = atoi(s);
    for (i=1; i <= DBI_gelCount(xx); i++) {
	if (DB_Number(xx,i) == gel) {
	    setCursorPosSeq(xx, 1, i);
	    edSetBriefSeqBase(xx, -1, -1, 1);
	    REDISPLAY(xx);
	    return 1;
	}
    }
    
    return 0;
}

static int findNextByName(EdStruct *xx, char *s)
/*
 * Search forwards from the cursor position until the sequence specified
 * is found. The cursor is positioned at the left end of the sequence,
 * if found.
 *
 */
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)+1;
    int i;
    int n;
    
    if (!*s) return 0;
    n = strlen(s);
    
    if ((i = posToIndex(xx,spos))==0) return 0;
    
    for (; i <= DBI_gelCount(xx); i++) {
	if (strncmp(DBgetGelName(xx,DBI_order(xx)[i]),s,n)==0) {
	    setCursorPosSeq(xx, 1, DBI_order(xx)[i]);
	    edSetBriefSeqBase(xx, -1, -1, 1);
	    REDISPLAY(xx);
	    return 1;
	}
    }
    
    
    return 0;
}


static int findPrevByName(EdStruct *xx, char *s)
/*
 * Search forwards from the cursor position until the sequence specified
 * is found. The cursor is positioned at the left end of the sequence,
 * if found.
 *
 */
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)-1;
    int i;
    int n;
    
    if (!*s) return 0;
    n = strlen(s);
    
    if ((i = posToIndex(xx,spos))==0) return 0;
    
    if (i == xx->cursorSeq && xx->cursorPos == 1)
	i--;

    for (; i >= 0 ; i--) {
	if (strncmp(DBgetGelName(xx,DBI_order(xx)[i]),s,n)==0) {
	    setCursorPosSeq(xx, 1, DBI_order(xx)[i]);
	    edSetBriefSeqBase(xx, -1, -1, 1);
	    REDISPLAY(xx);
	    return 1;
	}
    }
    
    
    return 0;
}

static int findNextGelByName (EdStruct *xx, char *s)
/*
 * Position cursor on left end of gel sequence matched by string s.
 * If s starts with a hash assume a number is specified.
 * Otherwise assume a gel name is specified.
 */
{
    if (*s) {
	if (*s == '#') {
	    return findGelByNumber(xx,++s);
	} else
	    return findNextByName(xx,s);
    }
    return 0;
}


static int findPrevGelByName (EdStruct *xx, char *s)
/*
 * Position cursor on left end of gel sequence matched by string s.
 * If s starts with a hash assume a number is specified.
 * Otherwise assume a gel name is specified.
 */
{
    if (*s) {
	if (*s == '#') {
	    return findGelByNumber(xx,++s);
	} else
	    return findPrevByName(xx,s);
    }
    return 0;
}





static int findNextTagByType (EdStruct *xx, char *type)
/*
 * Search forwards from the cursor position until a tag of a specified
 * type is encountered. The cursor is positioned at the left end of the
 * tag, if found.
 */
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)+1;
    int epos = DB_Length(xx,0);
    int fpos,fseq,i, seq;
    int fseqpos;
    int *seqList;
    tagStruct *found_tag = NULL;
    
    seqList = sequencesInRegion(xx,spos, epos);
    fseq = -1;
    fseqpos = 0;
    fpos = INT_MAX;
    i=-1;
    seq=0;

    do {
	/* search through tag list for sequence, starting with consensus */
	/* int seq = seqList[i]; */
	tagStruct *t;
	
	/* force sequence in memory */
	(void) DBgetSeq(DBI(xx), seq);

	t = (tagStruct *) DBgetTags(DBI(xx),seq);
	while (t != NULL) {
	    int normpos, tagpos, taglen;
	    normpos = normalisePos2(xx, seq, t->tagrec.position,
				    t->tagrec.length) - DB_Start(xx, seq);
	    if (!xx->reveal_cutoffs && normpos <= 0) {
		taglen = t->tagrec.length - 2 + normpos;
		normpos = 1;
	    } else {
		taglen = t->tagrec.length - 1;
	    }
	    tagpos=positionInContig(xx,seq,normpos);

	    if (in_interval(spos,fpos,tagpos) &&
		strncmp(t->tagrec.type.c,type,4)==0 &&
		(xx->reveal_cutoffs ||
		 (!xx->reveal_cutoffs &&
		  normpos + taglen > 0 &&
		  normpos <= DB_Length(xx, seq)))) {
		fseq = seq;
		fseqpos = normpos;
		fpos = tagpos;
		found_tag = t;
		if (DB_Comp(xx,seq) == COMPLEMENTED) 
		    /* Keep looking */
		    t = t->next;
		else
		    /* Stop now */
		    t = NULL;
	    } else
		t = t->next;
	}

	if (0 == (seq = seqList[++i]))
	    break;

    } while (DB_RelPos(xx,seqList[i]) - DB_Start(xx, seqList[i]) < fpos);
    
    
    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos, fseq);
	_select_tag(xx,fseq,found_tag);
	REDISPLAY(xx);
    }
    
    return fseq != -1;
}




static int findPrevTagByType (EdStruct *xx, char *type)
/*
 * Search backwards from the current cursor position until a tag of the
 * specified type is encountered. The cursor is positioned at the left
 * end of the tag, if found.
 */
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)-1;
    int epos = 1;
    int fpos,fseq,i;
    int fseqpos;
    int *seqList;
    tagStruct *found_tag = NULL;
    int max_len = dbi_max_gel_len(DBI(xx),0);
    
    seqList = sequencesInRegion(xx,epos, spos);
    fseq = -1;
    fseqpos = 0;
    fpos = -INT_MAX;
    
    for (i=0; seqList[i]; i++) ;
    for (; i>=0 && DB_RelPos(xx,seqList[i])+max_len > fpos ; i--) {
	/* search through tag list for sequence */
	int seq = seqList[i];
	tagStruct *t;
	
	/* force sequence in memory */
	(void) DBgetSeq(DBI(xx), seq);

	t = (tagStruct *) DBgetTags(DBI(xx),seq);
	while (t != NULL) {
	    int normpos, tagpos, taglen;
	    normpos = normalisePos2(xx, seq, t->tagrec.position,
				    t->tagrec.length) - DB_Start(xx, seq);
	    if (!xx->reveal_cutoffs && normpos <= 0) {
		taglen = t->tagrec.length - 2 + normpos;
		normpos = 1;
	    } else {
		taglen = t->tagrec.length - 1;
	    }
	    tagpos = positionInContig(xx,seq,normpos);
	    if (in_interval(spos,fpos,tagpos) &&
		strncmp(t->tagrec.type.c,type,4)==0 &&
		(xx->reveal_cutoffs ||
		 (!xx->reveal_cutoffs &&
		  normpos + taglen > 0 &&
		  normpos <= DB_Length(xx, seq)))) {
		fseq = seq;
		fseqpos = normpos;
		fpos = tagpos;
		found_tag = t;
		if (DB_Comp(xx,seq) != COMPLEMENTED) 
		    /* Keep looking */
		    t = t->next;
		else
		    /* Stop now */
		    t = NULL;
	    } else
		t = t->next;
	}
	
    }
    
    
    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos, fseq);
	_select_tag(xx,fseq,found_tag);
	REDISPLAY(xx);
    }
    
    return fseq != -1;
}


static int findNextSequence(EdStruct *xx, char *s, int mismatches, int strand,
			    int where)
/*
 * Search forwards from the cursor position until the sequence specified
 * is found. The cursor is positioned at the left end of the sequence,
 * if found.
 *
 * The search is made in the sequences and/or in the consensus.
 * The search is not case sensitive.
 *
 * Where=1 => readings
 * Where=2 => consensus
 * where=3 => both
 */
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)+1;
    int epos = DB_Length(xx,0);
    int fpos,fseq,i;
    int fseqpos;
    int *seqList;
    char *reading;
    int maxlen;
    char *uppert, *upperb;
    int patlen, len;
    char *cons = NULL;
    int cons_len = 0;
    int wloop;

    /* uppercase search string and store fwd/rev copies */
    patlen = strlen(s);
    depad_seq(s, &patlen, NULL);
    if (NULL == (uppert = (char *)xmalloc(patlen + 1)))
	return 0;
    if (NULL == (upperb = (char *)xmalloc(patlen + 1)))
	return 0;

    uppert[patlen] = upperb[patlen] = 0;
    for (i = patlen-1; i >= 0; i--) {
	upperb[i] = uppert[i] = toupper(s[i]);
    }
    complement_seq(upperb, patlen);
   
    seqList = sequencesInRegion(xx,spos, epos);
    fseq = -1;
    fseqpos = 0;
    fpos = INT_MAX;
    
    /* Identify maximum sequence length - we scan sequence by sequence */
    for (maxlen = 0, i=0; seqList[i] ; i++)
	maxlen = max(maxlen, xx->reveal_cutoffs ? DB_Length2(xx,seqList[i])
		     : DB_Length(xx,seqList[i]));

    if (NULL == (reading = (char *)xmalloc(maxlen+1)))
	return 0;

    /*
     * Note, bailing out of the search early is not possible due to
     * the fact that the sequences are sorted by "used" start position and
     * not cutoff start position.
     * Eg:
     * A ---------------------X-----------......
     * B                ...X....-----------------................
     * C                           .....-------------------...........
     * D           ....X...................----------.....
     *
     * We find A first, then B, but C looks tempting to abort the search on.
     * However then we'd miss D.
     */
    for (i=0; seqList[i]; i++) {
	/* search through tag list for sequence */
	int seq = seqList[i];
	char *str = NULL;
	char *indt, *indb, *ind;
	int offset;
	int diff;
	int cpos = 0;
	
	/* readings vs consensus */
	for (wloop = 0; wloop < 2; wloop++) {
	    char *data = NULL;

	    if (0 == wloop && 0 == (where & 1))
		continue;

	    if (1 == wloop && 0 == (where & 2))
		continue;

	    /* Search readings */
	    if (wloop == 0) {
		if (DB_RelPos(xx, seq) - DB_Start(xx, seq) > fpos)
		    continue;

		str = DBgetSeq(DBI(xx),seq);
		offset = 0;

		/*
		 * Fill reading[] array with sequence, including cutoffs if
		 * shown. Set 'str' to point to used sequence position
		 * in reading[].
		 */
		if (xx->reveal_cutoffs) {
		    char *lower = DB_Seq(xx, seq);
		    len = DB_Length2(xx, seq);
		    reading[len] = '\0';
		    for (len--; len >= 0; len--) {
			reading[len] = toupper(lower[len]);
		    }

		    str = reading + DB_Start(xx, seq);

		    if (in_interval(-DB_Start(xx, seq) + 1,
				    DB_Length2(xx, seq),
				    spos - DB_RelPos(xx, seq) + 1))
			offset = spos - (DB_RelPos(xx, seq)
				      - DB_Start(xx, seq));
		} else {
		    len = DB_Length(xx, seq);
		    reading[len] = '\0';
		    /* strncpy(reading,str,DB_Length(xx,seq)); */
		    for (len--; len >= 0; len--) {
			reading[len] = toupper(str[len]);
		    }

		    str = reading;

		    if (in_interval(1, DB_Length(xx,seq),
				    spos - DB_RelPos(xx,seq)+1))
			offset = spos - DB_RelPos(xx, seq);
		}

		data = reading + offset;
	    }

	    /* Also search consensus from this relpos to next relpos */
	    else if (wloop == 1) {
		int npads = 0;
		int toadd = 10;

		if (DB_RelPos(xx, seq) > fpos)
		    continue;

		do {
		    int p;

		    cpos = MAX(spos, DB_RelPos(xx, seqList[i]));
		    npads = toadd*2;
		    diff = patlen + (seqList[i+1]
				     ?  DB_RelPos(xx, seqList[i+1])
				     : DB_Length(xx, 0))
			          - DB_RelPos(xx, seq) + npads;
		    if (cons_len < diff) {
			cons_len = diff;
			cons = (char *)xrealloc(cons, cons_len+1);
		    }
		    DBcalcConsensus(xx, cpos, diff, cons,
				    NULL,BOTH_STRANDS);
		    cons[diff] = 0;

		    /*
		     * Count pads. Imagine a contig like this:
		     * R1        ------------
		     * R2                  -------------
		     * pattern            xxxx
		     * search    ..............
		     * 
		     * The search range is the distance from R1 to R2 plus
		     * the pattern length, so that patterns overlapping
		     * sequences are found.
		     * But the search length needs to overlap R2 by "patlen"
		     * real bases (excluding pads). We cycle around
		     * hoping to get this many.
		     */
		    for (toadd = 0, p = diff - patlen - npads; cons[p]; p++) {
			if (cons[p] == '*')
			    toadd++;
		    }
		} while (npads < toadd);

		data = cons;
	    }

	    indb = indt = NULL;
	    if (strand == '+' || strand == '=')
		indt = pstrstr_inexact(data, uppert, mismatches, NULL);
	    if (strand == '-' || strand == '=')
		indb = pstrstr_inexact(data, upperb, mismatches, NULL);
	
	    if (indt && indb)
		ind = MIN(indt, indb);
	    else if (indt)
		ind = indt;
	    else if (indb)
		ind = indb;
	    else
		ind = NULL;
	    
	    if (ind != NULL) {
		int pos;
		if (wloop == 0) {
		    pos = positionInContig(xx,seq,(int) (ind - str) + 1);
		    if (in_interval(spos,fpos,pos)) {
			fseqpos = (int) (ind - str) + 1;
			fpos = pos;
			fseq = seq;
		    }
		} else {
		    pos = ind - cons + cpos;
		    if (in_interval(spos,fpos,pos)) {
			fseqpos = pos;
			fpos = pos;
			fseq = 0;
		    }
		}
	    }
	}
    }
    
    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos, fseq);
	REDISPLAY(xx);
    }
    
    xfree(reading);
    xfree(uppert);
    xfree(upperb);
    xfree(cons);

    return fseq != -1;
}

static int findPrevSequence(EdStruct *xx, char *s, int mismatches, int strand,
			    int where)
/*
 * Search backwards from the cursor position until the sequence specified
 * is found. The cursor is positioned at the left end of the sequence,
 * if found.
 *
 * The search is made in the sequences and/or in the consensus.
 * The search is not case sensitive.
 *
 * Where=1 => readings
 * Where=2 => consensus
 * where=3 => both
 */
{
    int spos;
    int epos = 1;
    int fpos,fseq,i;
    int fseqpos;
    int *seqList;
    char *reading;
    char *uppert, *upperb;
    int maxlen;
    int patlen;
    int len;
    char *cons = NULL;
    int cons_len = 0;
    
    /* uppercase search string and store fwd/rev copies */
    patlen = strlen(s);
    depad_seq(s, &patlen, NULL);
    if (NULL == (uppert = (char *)xmalloc(patlen + 1)))
	return 0;
    if (NULL == (upperb = (char *)xmalloc(patlen + 1)))
	return 0;

    uppert[patlen] = upperb[patlen] = 0;
    for (i = patlen-1; i >= 0; i--) {
	upperb[i] = uppert[i] = toupper(s[i]);
    }
    complement_seq(upperb, patlen);
    
    spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)+strlen(s)-2;
    
    seqList = sequencesInRegion(xx,epos, spos);
    fseq = -1;
    fseqpos = 0;
    fpos = -INT_MAX;
     
    /* Identify maximum sequence length - we scan sequence by sequence */
    for (maxlen = 0, i=0; seqList[i] ; i++)
	maxlen = max(maxlen, xx->reveal_cutoffs ? DB_Length2(xx,seqList[i])
		     : DB_Length(xx,seqList[i]));

    if (NULL == (reading = (char *)xmalloc(maxlen+2)))
	return 0;

    for (i=0; seqList[i]; i++) ;
    for (i--; i>=0; i--) {
	/* search through tag list for sequence */
	int seq = seqList[i];
	char *str = NULL;
	char *indt, *indb, *ind;
	int wloop;
	int cpos, diff;
	
	for (wloop = 0; wloop < 2; wloop++) {
	    char *data = NULL;

	    if (0 == wloop && 0 == (where & 1))
		continue;

	    if (1 == wloop && 0 == (where & 2))
		continue;
	    

	    /* Search readings */
	    if (wloop == 0) {
		if (xx->reveal_cutoffs) {
		    if (DB_RelPos(xx, seq) - DB_Start(xx, seq)
			+ DB_Length2(xx, seq) <= fpos)
			continue;
		} else {
		    if (DB_RelPos(xx, seq) + DB_Length(xx, seq) <= fpos)
			continue;
		}

		str = DBgetSeq(DBI(xx),seq);

		/*
		 * Fill reading[] array with sequence, including
		 * cutoffs if shown. We null terminate it if it overlaps
		 * our search right-end position (spos).
		 */
		if (xx->reveal_cutoffs) {
		    char *lower = DB_Seq(xx, seq);
		    len = DB_Length2(xx, seq);
		    reading[len] = '\0';
		    for (len--; len >= 0; len--) {
			reading[len] = toupper(lower[len]);
		    }

		    str = reading + DB_Start(xx, seq);

		    if (in_interval(1, DB_Length2(xx, seq),
				    spos+1 - (DB_RelPos(xx, seq)
					      - DB_Start(xx, seq))))
			reading[spos+1 - (DB_RelPos(xx, seq)
					- DB_Start(xx, seq))] = 0;
		} else {
		    len = DB_Length(xx, seq);
		    reading[len] = '\0';
		    for (len--; len >= 0; len--) {
			reading[len] = toupper(str[len]);
		    }

		    str = reading;

		    if (in_interval(1, DB_Length(xx, seq),
				    spos+1 - DB_RelPos(xx, seq)))
			reading[spos+1 - DB_RelPos(xx, seq)] = 0;
		}

		data = reading;
	    }

	    /* Search consensus */
	    else if (wloop == 1) {
		int npads = 0;
		int toadd = 10;
		
		cpos = seqList[i+1]
		    ? DB_RelPos(xx, seqList[i+1])
		    : DB_Length(xx, 0);
		if (cpos < fpos)
		    continue;

		do {
		    int p;
		    
		    npads = toadd*2;
		    diff = MIN(spos+1, patlen + cpos + npads)
			- DB_RelPos(xx, seqList[i]);
		    if (diff <= 0)
			break;
		    if (cons_len < diff) {
			cons_len = diff;
			cons = (char *)xrealloc(cons, cons_len+1);
		    }
		    DBcalcConsensus(xx, DB_RelPos(xx, seqList[i]), diff, cons,
				    NULL,BOTH_STRANDS);
		    cons[diff] = 0;

		    /*
		     * Count pads. Imagine a contig like this:
		     * R1        ------------
		     * R2                  -------------
		     * pattern            xxxx
		     * search    ..............
		     * 
		     * The search range is the distance from R1 to R2 plus
		     * the pattern length, so that patterns overlapping
		     * sequences are found.
		     * But the search length needs to overlap R2 by "patlen"
		     * real bases (excluding pads). We cycle around
		     * hoping to get this many.
		     */
		    for (toadd = 0, p = diff - patlen - npads; cons[p]; p++) {
			if (cons[p] == '*')
			    toadd++;
		    }
		} while (npads < toadd);
		if (diff <= 0)
		    continue;

		data = cons;
	    }

	    indb = indt = NULL;
	    if (strand == '+' || strand == '=')
		indt = prstrstr_inexact(data, uppert, mismatches, NULL);
	    if (strand == '-' || strand == '=')
		indb = prstrstr_inexact(data, upperb, mismatches, NULL);
	    
	    if (indt && indb)
		ind = MAX(indt, indb);
	    else if (indt)
		ind = indt;
	    else if (indb)
		ind = indb;
	    else
		ind = NULL;

	    if (ind != NULL) {
		int pos;
		if (wloop == 0) {
		    pos = positionInContig(xx,seq,(int) (ind - str) + 1);
		    if (in_interval(spos,fpos,pos)) {
			fseqpos = (int) (ind - str) + 1;
			fpos = pos;
			fseq = seq;
		    }
		} else {
		    pos = ind - cons + DB_RelPos(xx, seq);
		    if (in_interval(spos,fpos,pos)) {
			fseqpos = pos;
			fpos = pos;
			fseq = 0;
		    }
		}
	    }
	}
    }
    
    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos, fseq);
	REDISPLAY(xx);
    }
    
    xfree(reading);
    xfree(uppert);
    xfree(upperb);
    return fseq != -1;
}





static int findNextAnno(EdStruct *xx, char *anno)
/*
 * Search forwards from the cursor position until a tag containing the
 * specified annotation is found. The cursor is positioned at the left
 * end of the tag, if found.
 *
 * Observations:
 *   A regular expression search is found, giving unpredictable results
 *   to people unfamiliar with such searches
 */
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)+1;
    int epos = DB_Length(xx,0);
    int fpos,fseq,i,seq;
    int fseqpos;
    int *seqList;
    char *r_exp;
    char *find_all = "$";
    tagStruct *found_tag = NULL;
    
    if (! *anno) anno = find_all;
    
    if (NULL == (r_exp = REGCMP(EDINTERP(xx->ed),anno))) {
	verror(ERR_WARN, "findNextAnno", "invalid regular expression");
	return 0;
    }

    seqList = sequencesInRegion(xx,spos, epos);
    fseq = -1;
    fseqpos = 0;
    fpos = INT_MAX;
    
    i=-1;
    seq=0;

    do {
	/* search through tag list for sequence, starting with consensus */
	tagStruct *t;
	
	/* force sequence in memory */
	(void) DBgetSeq(DBI(xx), seq);

	t = (tagStruct *) DBgetTags(DBI(xx),seq);
	while (t != NULL) {
	    int normpos, tagpos, taglen;

	    /* Check this tag type is displayed */
	    if (!xx->tag_list[idToIndex(t->tagrec.type.c)]) {
		t = t->next;
		continue;
	    }

	    normpos = normalisePos2(xx, seq, t->tagrec.position,
				    t->tagrec.length) - DB_Start(xx, seq);
	    if (!xx->reveal_cutoffs && normpos <= 0) {
		taglen = t->tagrec.length - 2 + normpos;
		normpos = 1;
	    } else {
		taglen = t->tagrec.length - 1;
	    }
	    tagpos=positionInContig(xx,seq,normpos);

	    /*
	     * Get annotation
	     */
	    force_comment(DBI_io(xx), t);
	    if (in_interval(spos,fpos,tagpos) &&
		t->tagrec.type.c[0] != '*' &&  /* avoid special tags */
		t->tagrec.position != 0 &&
		(xx->reveal_cutoffs ||
		 (!xx->reveal_cutoffs &&
		  normpos + taglen > 0 &&
		  normpos <= DB_Length(xx, seq))) &&
		(REGEX(EDINTERP(xx->ed),t->newcomment, r_exp) == 1)) {
		fseq = seq;
		fseqpos = normpos;
		fpos = tagpos;
		found_tag = t;
		if (DB_Comp(xx,seq) == COMPLEMENTED) 
		    /* Keep looking */
		    t = t->next;
		else
		    /* Stop now */
		    t = NULL;
	    } else
		t = t->next;
	}

	if (0 == (seq = seqList[++i]))
	    break;

    } while (DB_RelPos(xx,seqList[i]) - DB_Start(xx, seqList[i]) < fpos);

        
    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos, fseq);
	_select_tag(xx,fseq,found_tag);
	REDISPLAY(xx);
    }
    
    REGFREE(EDINTERP(xx->ed),r_exp);
    return fseq != -1;
}

static int findPrevAnno(EdStruct *xx, char *anno)
/*
 * Search backwards from the cursor position until a tag containing the
 * specified annotation is found. The cursor is positioned at the left
 * end of the tag, if found.
 *
 * Observations:
 *   A regular expression search is found, giving unpredictable results
 *   to people unfamiliar with such searches
 */
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)-1;
    int epos = 1;
    int fpos,fseq,i;
    int fseqpos;
    int *seqList;
    int max_len = dbi_max_gel_len(DBI(xx),0);
    char *r_exp;
    char *find_all = "$";
    tagStruct *found_tag = NULL;
    
    if (! *anno) anno = find_all;
    
    if (NULL == (r_exp = REGCMP(EDINTERP(xx->ed),anno))) {
	verror(ERR_WARN, "findPrevAnno", "invalid regular expression");
	return 0;
    }

    seqList = sequencesInRegion(xx,epos, spos);
    fseq = -1;
    fseqpos = 0;
    fpos = -INT_MAX;
    
    for (i=0; seqList[i]; i++) ;
    for (; i>=0 && DB_RelPos(xx,seqList[i])+max_len > fpos ; i--) {
	/* search through tag list for sequence */
	int seq = seqList[i];
	tagStruct *t;
	
	/* force sequence in memory */
	(void) DBgetSeq(DBI(xx), seq);

	t = (tagStruct *) DBgetTags(DBI(xx),seq);
	while (t != NULL) {
	    int normpos, tagpos, taglen;

	    /* Check this tag type is displayed */
	    if (!xx->tag_list[idToIndex(t->tagrec.type.c)]) {
		t = t->next;
		continue;
	    }

	    normpos = normalisePos2(xx, seq,t->tagrec.position,
				    t->tagrec.length) - DB_Start(xx, seq);
	    if (!xx->reveal_cutoffs && normpos <= 0) {
		taglen = t->tagrec.length - 2 + normpos;
		normpos = 1;
	    } else {
		taglen = t->tagrec.length - 1;
	    }
	    tagpos=positionInContig(xx,seq,normpos);

	    /*
	     * Get annotation
	     */
	    force_comment(DBI_io(xx), t);
	    if (in_interval(spos,fpos,tagpos) &&
		t->tagrec.type.c[0] != '*' &&  /* avoid special tags */
		t->tagrec.position != 0 &&
		(xx->reveal_cutoffs ||
		 (!xx->reveal_cutoffs &&
		  normpos + taglen > 0 &&
		  normpos <= DB_Length(xx, seq))) &&
		(REGEX(EDINTERP(xx->ed),t->newcomment, r_exp) == 1)) {
		fseq = seq;
		fseqpos = normpos;
		fpos = tagpos;
		found_tag = t;
		if (DB_Comp(xx,seq) != COMPLEMENTED) 
		    /* Keep looking */
		    t = t->next;
		else
		    /* Stop now */
		    t = NULL;
	    } else
		t = t->next;
	}
	
    }
    
    
    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos, fseq);
	_select_tag(xx,fseq,found_tag);
	REDISPLAY(xx);
    }
    
    REGFREE(EDINTERP(xx->ed),r_exp);
    return fseq != -1;
}



static int findPosition(EdStruct *xx, char *text_pos)
/*
 * Position the cursor at the position specified.
 * There are four modes:
 *   1. By position in contig.       eg 30717
 *   2. By position in current gel.  eg @100
 *   3. By position in defined gel.  eg @100/xy66a2.s1
 *   4. By a relative offset.        eg +1000 eg -1000
 *
 * In addition to this, the prefix of a 'u' indicates an unpadded
 * consensus position. eg u100/xy66a2.s1. This doesn't apply for
 * positions relative to the readings.
 *
 * Observations:
 *   The cursor is positioned in the same gel if possible.
 *   If it possible to specify negative or large numbers when 
 *   specifying position in gel.
 */
{
    
    int pos;
    int cseq = xx->cursorSeq;
    int cpos = xx->cursorPos;
    char *cp;
    int unpadded;
    
    for(; *text_pos && isspace((int)*text_pos) ; text_pos++) ;
    
    if (*text_pos == 'u' || *text_pos == 'U') {
	unpadded = 1;
	text_pos++;
    } else {
	unpadded = 0;
    }

    switch (*text_pos) {
    case '\0':
	return 0;
    case '+':
    case '-':
	pos = positionInContig(xx,cseq,cpos) + atoi(text_pos);
	break;
    case '@':
	if (NULL != (cp = strchr(text_pos, '/'))) {
	    int i;
	    int rnum;

	    rnum = get_gel_num(DBI_io(xx), cp+1, GGN_ID);

	    /* Find the sequence, and take pos relative to this one */
	    cseq = 0;
	    for (i=0; i <= DBI_gelCount(xx); i++) {
		if (DB_Number(xx, i) == rnum) {
		    cseq = i;
		    break;
		}
	    }
	    if (!cseq)
 		return 0;
	    xx->cursorSeq = cseq;
	}
	pos = DB_RelPos(xx,cseq) + atoi(++text_pos) - 1;
	break;
    default:
	pos = atoi(text_pos);
	if (DBI(xx)->reference_seq) {
	    char *bases = DBgetSeq(DBI(xx), DBI(xx)->reference_seq);
	    int i, j;

	    pos -= DBI(xx)->reference_offset-1;
	    if (DBI(xx)->reference_len)
		pos = (pos-1) % DBI(xx)->reference_len;
	    while (pos < 0)
		pos += DBI(xx)->reference_len;
	    pos++;

	    for (j = i = 0; j < pos &&
		 i < DB_Length2(xx, DBI(xx)->reference_seq); i++) {
		if (bases[i] != '*')
		    j++;
	    }
	    pos = i;
	}
	break;
    }

    if (unpadded && !DBI(xx)->reference_seq) {
	char *con;
	int i, upos, len;

	len = DB_Length(xx, 0);

	if (con = xmalloc(len+1)) {
	    DBcalcConsensus(xx, 1, len, con, NULL, BOTH_STRANDS);

	    for (i = upos = 0; i < len && upos != pos; i++) {
		if (con[i] != '*')
		    upos++;
	    }
	    
	    xfree(con);
	    pos = i;
	}
    }
    
    if ((xx->reveal_cutoffs &&
	 in_interval(DB_RelPos(xx, cseq) - DB_Start(xx, cseq),
		     DB_RelPos(xx, cseq) + DB_Length2(xx, cseq) -
		     DB_Start(xx, cseq),
		     pos)) ||
	(!xx->reveal_cutoffs &&
	 in_interval(DB_RelPos(xx, cseq),
		     DB_RelPos(xx, cseq) + DB_Length(xx, cseq),
		     pos ))) {
	setCursorPos(xx, pos - DB_RelPos(xx, cseq) + 1);
    } else {
	if (pos > 0 && pos <= DB_Length(xx, 0)) {
	    setCursorPosSeq(xx, pos, 0);
	    edSetBriefSeqBase(xx, -1, -1, 1);
	} else
	    return 0;
    }

    REDISPLAY(xx);
    return 1;
}



static int findNextProblem (EdStruct *xx)
/*
 * Search forward from the cursor position until a the consensus is not
 * A, C, G or T. The cursor is positioned on the problem base, if found.
 */
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)+1;
    int epos = DB_Length(xx,0);
    char buffer[SEARCH_CHUNKS+1];
    float buffer2[SEARCH_CHUNKS+1];
    int i,j,width;
    
    if (spos < 1)
	spos = 1;

    for (i=spos; i<= epos; i+=SEARCH_CHUNKS) {
	width = min(epos-i+1,SEARCH_CHUNKS);
	DBcalcConsensus (xx,i,width, buffer, buffer2, BOTH_STRANDS);

	for (j=0; j<width; j++) {
	    if ((buffer[j] == '*' && (buffer2[j] != 1.0 || xx->qual_cut < 0))
		 || (buffer[j] == '-')
		 || (buffer[j] == 'N')) {
		/* we have a problem! */
		
		/* _cursor_set(DBI(xx), 0, i+j); */
		setCursorPosSeq(xx, i+j, 0);
		edSetBriefSeqBase(xx, -1, -1, 1);
		REDISPLAY(xx);
		
		if (xx->display_traces)
		    tman_problem_traces(xx, i+j);

		return 1;
	    }
	}
    }

    return 0;
}



static int findPrevProblem (EdStruct *xx)
/*
 * Search forward from the cursor position until a the consensus is not
 * A, C, G or T. The cursor is positioned on the problem base, if found.
 */
{
#define SEARCH_CHUNKS MAX_DISPLAY_WIDTH
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)-1;
    int epos = 1;
    char buffer[SEARCH_CHUNKS+1];
    float buffer2[SEARCH_CHUNKS+1];
    int i,width;
    int j;
    
    if (spos > DB_Length(xx, 0))
	spos = DB_Length(xx, 0);

    for (i=spos; i>= epos; i-=SEARCH_CHUNKS) {
	width = min(i-epos+1,SEARCH_CHUNKS);
	DBcalcConsensus (xx,i-width+1,width, buffer, buffer2, BOTH_STRANDS);

	for (j=width-1; j>=0; j--) {
	    if ((buffer[j] == '*' && (buffer2[j] != 1.0 || xx->qual_cut < 0))
		 || (buffer[j] == '-')
		 || (buffer[j] == 'N')) {
		/* we have a problem! */

		setCursorPosSeq(xx, i-width+1+j, 0);
		edSetBriefSeqBase(xx, -1, -1, 1);
		REDISPLAY(xx);

		if (xx->display_traces)
		    tman_problem_traces(xx, i-width+1+j);

		return 1;
	    }
	}
    }
    return 0;
}



static int findNextQualProb (EdStruct *xx)
/*
 * Search forwards from the cursor position until a problem relating to
 * quality is found. The cursor is positioned on the problematic base, if found.
 *
 * Observations:
 *   Large stretches of sequence on one strand only could cause frustration
 */
{
#define SEARCH_CHUNKS MAX_DISPLAY_WIDTH
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)+1;
    int epos = DB_Length(xx,0);
    char buffer[SEARCH_CHUNKS+1];
    int i,width;
    int j;
    
    if (spos < 1)
	spos = 1;

    for (i=spos; i<= epos; i+=SEARCH_CHUNKS) {
	width = min(epos-i+1,SEARCH_CHUNKS);
	calc_quality(0, i, i+width, buffer, xx->con_cut, xx->qual_cut, 
		     contEd_info, (void *)xx);
	for (j=0;j<width;j++) {
	    if ((xx->consensus_mode == CONSENSUS_MODE_CONFIDENCE &&
		 buffer[j] != R_GOOD_GOOD_EQ &&
		 buffer[j] != R_GOOD_BAD &&
		 buffer[j] != R_BAD_GOOD &&
		 buffer[j] != R_BAD_BAD) ||
		(xx->consensus_mode != CONSENSUS_MODE_CONFIDENCE &&
		 (buffer[j] != R_GOOD_GOOD_EQ))) {
		/* we have problem! */

		setCursorPosSeq(xx, i+j, 0);
		edSetBriefSeqBase(xx, -1, -1, 1);
		REDISPLAY(xx);

		if (xx->display_traces)
		    tman_problem_traces(xx, i+j);

		return 1;
	    }
	}
    }
    return 0;
}



static int findPrevQualProb (EdStruct *xx)
/*
 * Search backwards from the cursor position until a problem relating to
 * quality is found. The cursor is positioned on the problematic base, if found.
 *
 * Observations:
 *   Large stretches of sequence on one strand only could cause frustration
 */
{
#define SEARCH_CHUNKS MAX_DISPLAY_WIDTH
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)-1;
    int epos = 1;
    char buffer[SEARCH_CHUNKS+1];
    int i,width;
    int j;
    
    if (spos > DB_Length(xx, 0))
	spos = DB_Length(xx, 0);

    for (i=spos; i>= epos; i-=SEARCH_CHUNKS) {
	width = min(i-epos+1,SEARCH_CHUNKS);
	calc_quality(0, i - width + 1, i + 1, buffer,
		     xx->con_cut, xx->qual_cut, 
		     contEd_info, (void *)xx);
	for (j=width-1; j>=0; j--) {
	    if ((xx->consensus_mode == CONSENSUS_MODE_CONFIDENCE &&
		 buffer[j] != R_GOOD_GOOD_EQ &&
		 buffer[j] != R_GOOD_BAD &&
		 buffer[j] != R_BAD_GOOD &&
		 buffer[j] != R_BAD_BAD) ||
		(xx->consensus_mode != CONSENSUS_MODE_CONFIDENCE &&
		 (buffer[j] != R_GOOD_GOOD_EQ))) {
		/* we have problem! */

		setCursorPosSeq(xx, i-width+1+j, 0);
		edSetBriefSeqBase(xx, -1, -1, 1);
		REDISPLAY(xx);

		if (xx->display_traces)
		    tman_problem_traces(xx, i-width+1+j);

		return 1;
	    }
	}
    }
    return 0;
}


/*
 * Returns 1 or 0 to indicate whether position 'pos' of sequence 'seq' has
 * been edited. We are scanning in direction 'dir' and 'comp' tells us whether
 * the sequence has been complimented or not.
 *
 * Call will 'xx == NULL' to initialise.
 */
static int is_edit(EdStruct *xx, int seq, int pos, int dir, int comp) {
    static int last;
    int c1, o1, o2 = 0, tmp = last;
    char b1;
    
    if (!xx) {
	last = 0;
	return 0;
    }

    /*
     * At present the opos array is 16 bit, so we do not report edits in
     * longer sequences.
     */
    if (ABS(DB_Length2(xx, seq)) >= 32768)
	return 0;

    if (pos >= 0 && pos < DB_Length2(xx, seq)) {
	c1 = DB_Conf(xx, seq)[pos];

	/* remember last non zero original position */
	o1 = DB_Opos(xx, seq)[pos];
	if (o1) last = o1;
    
	b1 = DB_Seq(xx, seq)[pos];

	if (pos + dir >= 0 && pos + dir < DB_Length2(xx, seq)) {
    
	    /* last non zero orig position doesn't tally with next base */
	    o2 = DB_Opos(xx, seq)[pos+dir];
    
	    /* if both bases are non zero - next base must tally with us */
	    if (o1 && o2 && o1 != o2 + comp) {
		vmessage("%d base(s) to the right of the cursor deleted\n",
			 abs(o1 - (o2 + comp)));
		return 1;
	    }
	}
	
	if (!o1) {
	    if (o2 && tmp && o2 + comp != tmp) {
		vmessage("Base type or confidence changed\n");
		return 1;
	    } else if (b1 != '*') {
		vmessage("Base inserted (or changed from pad)\n");
		return 1;
	    } else if (b1 == '*' && (c1 == 0 || c1 == 100)) {
		vmessage("Pad confidence changed\n");
		return 1;
	    }
	}
    }

    /* all clear */
    return 0;
}

/*
 * Search forward from the cursor position until we find an edited base.
 * This is a base with a original position out of step, or with a
 * quality value of 0 or 100.
 */
static int findNextEdit(EdStruct *xx) {
    int spos = positionInContig(xx, xx->cursorSeq, xx->cursorPos)+1;
    int epos = DB_Length(xx, 0);
    int fseq,fpos,fseqpos, i, j;
    int *seqList;
    
    seqList = sequencesInRegion(xx,spos, epos);
    fseq = -1;
    fseqpos = 0;
    fpos = INT_MAX;
    
    for (i=0; seqList[i] && DB_RelPos(xx,seqList[i]) < fpos;
	 i++) {
	int seq = seqList[i], e, dir;

	/* force sequence in memory */
	(void) DBgetSeq(DBI(xx), seq);

	j = DB_RelPos(xx, seq) < spos
	    ? spos - DB_RelPos(xx, seq) + DB_Start(xx, seq)
		: DB_Start(xx, seq);

	e = DB_RelPos(xx, seq) + DB_Length(xx, seq) > fpos
	    ? fpos - DB_RelPos(xx, seq) + DB_Start(xx, seq)
		: DB_Length(xx, seq) + DB_Start(xx, seq);

	dir = DB_Comp(xx, seq) == -1 ? 1 : -1;

	(void)is_edit(0, 0, 0, 0, 0);
	for (; j < e; j++) {
	    if (is_edit(xx, seq, j, 1, dir)) {
		fseqpos = j + 1 - DB_Start(xx, seq);
		fseq = seq;
		fpos = positionInContig(xx, seq, fseqpos);

		break;
	    }
	}
    }

    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos, fseq);
	edSetBriefSeqBase(xx, -1, -1, 1);
	REDISPLAY(xx);
    }

    return fseq != -1;
}

/*
 * Search backwards from the cursor position until we find an edited base.
 * This is a base with a original position out of step, or with a
 * quality value of 0 or 100.
 */
static int findPrevEdit(EdStruct *xx) {
    int spos = positionInContig(xx, xx->cursorSeq, xx->cursorPos)-1;
    int epos = 1;
    int fseq,fpos,fseqpos, i, j;
    int *seqList;
    int max_len = dbi_max_gel_len(DBI(xx),0);
    
    seqList = sequencesInRegion(xx,epos, spos);
    fseq = -1;
    fseqpos = 0;
    fpos = epos;
    
    for (i=0; seqList[i]; i++) ;
    for (i--; i >= 0 && DB_RelPos(xx,seqList[i]) >= epos &&
	     DB_RelPos(xx, seqList[i]) + max_len > fpos; i--) {
	int seq = seqList[i], e, dir;

	/* force sequence in memory */
	(void) DBgetSeq(DBI(xx), seq);

	j = spos - DB_RelPos(xx, seq) + DB_Start(xx, seq);
	if (j > DB_Start(xx, seq) + DB_Length(xx, seq) - 1)
	    j = DB_Start(xx, seq) + DB_Length(xx, seq) - 1;
	if (j < 0)
	    j = 0;

	e = DB_RelPos(xx, seq) < fpos
	    ? fpos - DB_RelPos(xx, seq) + DB_Start(xx, seq)
		: DB_Start(xx, seq);

	dir = DB_Comp(xx, seq) == -1 ? -1 : 1;

	(void)is_edit(0, 0, 0, 0, 0);
	for (; j >= e; j--) {
	    if (is_edit(xx, seq, j, -1, dir)) {
		fseqpos = j + 1 - DB_Start(xx, seq);
		fseq = seq;
		fpos = positionInContig(xx, seq, fseqpos);
		
		break;
	    }
	}
    }

    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos, fseq);
	edSetBriefSeqBase(xx, -1, -1, 1);
	REDISPLAY(xx);
    }

    return fseq != -1;
}

#define CSEARCH_CHUNKS 4096

/*
 * Checks for a problem. Returns the position if there is a problem, or -1
 * if there is no problem.
 * Mode is 0 for checking both strands independantly or 1 for
 * both strands together.
 */
static int cop_check(EdStruct *xx, int spos, char *ok, int from, int to,
	      int mode, int from_p, int to_p, int from_m, int to_m) {
    int j;

    from   -= spos;
    to     -= spos;
    from_p -= spos;
    to_p   -= spos;
    from_m -= spos;
    to_m   -= spos;

    if (mode) {
	for (j = from; j < to; j++) {
	    if (ok[j] == 0) {
		return j + spos;
	    }
	}
    } else {

	for (j = from; j < to; j++) {
	    if ((ok[j] & 3) != 3) {
		char c1, c2;
		float q1, q2;

		/*
		 * Check if this is a single stranded region, in which case
		 * we'll ignore this problem. (Verify is designed to spot
		 * conflicts, not lack of information.)
		 */
		if (ok[j] != 0) {
		    if ((ok[j] & 1) == 0) { /* _p */
			if (j <= from_p || j >= to_p) {
			    continue;
			}
		    }

		    if ((ok[j] & 2) == 0) { /* _m */
			if (j <= from_m || j >= to_m) {
			    continue;
			}
		    }
		}

		/*
		 * Inefficient to use calc_consensus lots of time on single
		 * bases rather than just once.
		 *
		 * However cop is designed to be ran at the end of a project
		 * and so the above conditional ought to mean this code is
		 * executed rarely.
		 */
		calc_consensus(0, j + spos, j + spos, CON_SUM, &c1, &c2,
			       &q1, &q2, xx->con_cut, xx->qual_cut,
			       contEd_info, (void *)xx);

		if (!(ok[j] & 1) && c1 != '*')
		    return j + spos;
		
		if (!(ok[j] & 2) && c2 != '*')
		    return j + spos;
	    }
	}
    }

    return -1;
}

/*
 * Search forwards from the cursor position until we find a consensus base
 * with little justification.
 *
 * ok[] is an array with bit 0 being set if the positive strand is ok and
 * bit 1 being set if the negative strand is ok.
 * 
 * We use the 'mode' argument to determine what we consider a problem.
 * mode == 1 implies check both strands together; mode == 0 implies
 * check each independently.
 */
static int findNextCop(EdStruct *xx, int mode) {
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos) + 1;
    int epos = DB_Length(xx,0);
    char *ok;
    int i,j;
    int *seqList;
    int from, to;
    int from_p, to_p; /* these hold the last covered data for the plus */
    int from_m, to_m; /* and minus strands */
    int fseq = -1, fseqpos = 0;
    int bit_[3], *bit = &bit_[1];

    if (spos > epos)
	return 0;

    bit[+1] = 0x01;
    bit[-1] = 0x02;

    /* FIXME: we should do this buffer at a time rather than all at once */
    if (NULL == (ok = (char *)xcalloc(epos - spos + 1, 1))) {
	return -1;
    }

    seqList = sequencesInRegion(xx,spos, epos);

    from = from_p = from_m = spos;
    to_p = to_m = spos;
    for (i = 0; seqList[i]; i++) {
	int s = seqList[i];
	int start, end, dir, last;

	/*
	 * At present the opos array is 16 bit, so we do not report edits in
	 * longer sequences.
	 */
	if (ABS(DB_Length2(xx, s)) >= 32768)
	    continue;

	to = DB_RelPos(xx, s);

	if (from > from_m)
	    from_m = from;
	if (from > from_p)
	    from_p = from;

	if (to_m < from_m)
	    to_m = from_m;
	if (to_p < from_p)
	    to_p = from_p;


	if (-1 != (j = cop_check(xx, spos, ok, from, to, mode,
				 from_p, to_p, from_m, to_m))) {
	    fseqpos = j;
	    fseq = 0;
	    break;
	}

	/* force opos to be in memory */
	(void)DBgetSeq(DBI(xx), s);

	start = (spos > DB_RelPos(xx, s) ? spos - DB_RelPos(xx, s) : 0)
	    + DB_Start(xx, s);
	end = epos < DB_RelPos(xx, s) + DB_Length(xx, s) ?
	    epos - DB_RelPos(xx, s) + DB_Start(xx, s) : DB_End(xx, s) - 1;

	if (DB_Comp(xx, s) == COMPLEMENTED) {
	    dir = -1;
	    from_m = to;
	    to_m = MAX(to_m, DB_RelPos(xx, s) + DB_Length(xx, s));
	} else {
	    dir = 1;
	    from_p = to;
	    to_p = MAX(to_p, DB_RelPos(xx, s) + DB_Length(xx, s));
	}

	/*
	 * Original trace base for base x (start <= x <= end)
	 *     seq->base[DB_Opos(xx, s)[x] - 1];
	 * provided that DB_Opos(xx, s)[x] != 0.
	 *
	 * This base corresponds to the consensus base at position
	 *    pos = DB_RelPos(xx, s) + x - start;
	 *    cons[pos-spos];
	 *
	 * However we can compute the differences looking only at original
	 * positions. The three cases are:
	 *
	 * Insertion:
	 *   Opos  1 2 0 3 4
	 *   Edit  A C a G T
	 *   Orig  A C   G T
	 *
	 * Deletion:
	 *   Opos  1 2 4
	 *   Edit  A C T
	 *   Orig  A C T
	 *
	 * Change:
	 *   Opos  1 2 0 4
	 *   Edit  A C a T
	 *   Orig  A C G T
	 */
	last = DB_Opos(xx, s)[start] - dir;
	for (j = start; j < end; j++) {
	    int posc = DB_RelPos(xx, s) -DB_Start(xx, s) + j - spos;
	    int post = DB_Opos(xx, s)[j];

	    if (post != 0) {
		if (post == last + dir)
		    ok[posc] |= bit[dir];

		last = post;
	    }
	}
	
	if (to > spos)
	    from = to;
    }

    if (fseq == -1) {
	to = epos;
	
	if (from > from_m)
	    from_m = from;
	if (from > from_p)
	    from_p = from;

	if (to_m < from_m)
	    to_m = from_m;
	if (to_p < from_p)
	    to_p = from_p;

	if (-1 != (j = cop_check(xx, spos, ok, from, to, mode,
				 from_p, to_p, from_m, to_m))) {
	    fseqpos = j;
	    fseq = 0;
	}
    }

    xfree(ok);

    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos > 0 ? fseqpos : 1, fseq);
	edSetBriefSeqBase(xx, -1, -1, 1);
	REDISPLAY(xx);

	if (xx->display_traces)
	    tman_problem_traces(xx, DB_RelPos(xx, fseq)-1 +
				(fseqpos > 0 ? fseqpos : 1));

	return 1;
    }

    return 0;
}

/*
 * Search backwards from the cursor position until we find a consensus base
 * with little justification.
 */
static int findPrevCop(EdStruct *xx, int mode) {
    vmessage("Backwards verify searches are currently unimplemented.\n");

    return 0;
}


/*
 * Search forward from the cursor position until a the consensus quality is
 * below a particular threshold
 * The cursor is positioned on the problem base, if found.
 */
static int findNextConsQual (EdStruct *xx, int value)
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)+1;
    int epos = DB_Length(xx,0);
    char buffer[SEARCH_CHUNKS+1];
    float buffer2[SEARCH_CHUNKS+1];
    int i,j,width;
    float fval = value + 0.5;
    
    if (spos < 1)
	spos = 1;

    for (i=spos; i<= epos; i+=SEARCH_CHUNKS) {
	width = min(epos-i+1,SEARCH_CHUNKS);
	DBcalcConsensus (xx,i,width, buffer, buffer2, BOTH_STRANDS);

	for (j=0; j<width; j++) {
	    if (buffer2[j] < fval) {
		setCursorPosSeq(xx, i+j, 0);
		edSetBriefSeqBase(xx, -1, -1, 1);
		REDISPLAY(xx);

		if (xx->display_traces)
		    tman_problem_traces(xx, i+j);

		return 1;
	    }
	}
    }

    return 0;
}


/*
 * Search backward from the cursor position until a the consensus quality is
 * below a particular threshold
 * The cursor is positioned on the problem base, if found.
 */
static int findPrevConsQual (EdStruct *xx, int value)
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)-1;
    int epos = 1;
    char buffer[SEARCH_CHUNKS+1];
    float buffer2[SEARCH_CHUNKS+1];
    int i,width;
    int j;
    float fval = value + 0.5;
    
    if (spos > DB_Length(xx, 0))
	spos = DB_Length(xx, 0);

    for (i=spos; i>= epos; i-=SEARCH_CHUNKS) {
	width = min(i-epos+1,SEARCH_CHUNKS);
	DBcalcConsensus (xx,i-width+1,width, buffer, buffer2, BOTH_STRANDS);

	for (j=width-1; j>=0; j--) {
	    if (buffer2[j] < fval) {
		setCursorPosSeq(xx, i-width+1+j, 0);
		edSetBriefSeqBase(xx, -1, -1, 1);
		REDISPLAY(xx);

		if (xx->display_traces)
		    tman_problem_traces(xx, i-width+1+j);

		return 1;
	    }
	}
    }
    return 0;
}


/*
 * Search forward from the cursor position until a the consensus discrepancy
 * value is >= a particular threshold
 * The cursor is positioned on the problem base, if found.
 */
static int findNextConsDiscrep (EdStruct *xx, int value)
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)+1;
    int epos = DB_Length(xx,0);
    float buffer2[SEARCH_CHUNKS+1];
    int i,j,width;
    float fval = value + 0.5;
    
    if (spos < 1)
	spos = 1;

    for (i=spos; i<= epos; i+=SEARCH_CHUNKS) {
	width = min(epos-i+1,SEARCH_CHUNKS);
	DBcalcDiscrepancies (xx,i,width, buffer2);

	for (j=0; j<width; j++) {
	    if (buffer2[j] >= fval) {
		setCursorPosSeq(xx, i+j, 0);
		edSetBriefSeqBase(xx, -1, -1, 1);
		REDISPLAY(xx);

		if (xx->display_traces)
		    tman_problem_traces(xx, i+j);

		return 1;
	    }
	}
    }

    return 0;
}


/*
 * Search backward from the cursor position until a the consensus discrepancy
 * value is >= a particular threshold
 * The cursor is positioned on the problem base, if found.
 */
static int findPrevConsDiscrep (EdStruct *xx, int value)
{
    int spos = positionInContig(xx,xx->cursorSeq,xx->cursorPos)-1;
    int epos = 1;
    float buffer2[SEARCH_CHUNKS+1];
    int i,width;
    int j;
    float fval = value + 0.5;
    
    if (spos > DB_Length(xx, 0))
	spos = DB_Length(xx, 0);

    for (i=spos; i>= epos; i-=SEARCH_CHUNKS) {
	width = min(i-epos+1,SEARCH_CHUNKS);
	DBcalcDiscrepancies (xx,i-width+1,width, buffer2);

	for (j=width-1; j>=0; j--) {
	    if (buffer2[j] >= fval) {
		setCursorPosSeq(xx, i-width+1+j, 0);
		edSetBriefSeqBase(xx, -1, -1, 1);
		REDISPLAY(xx);

		if (xx->display_traces)
		    tman_problem_traces(xx, i-width+1+j);

		return 1;
	    }
	}
    }
    return 0;
}


/*
 * Search forward from the cursor position until we find a place with
 * two or more bases with a quality > 'value' that are not the same base type.
 */
static int findNextDiscrepancy (EdStruct *xx, int value)
{
    int spos = positionInContig(xx, xx->cursorSeq, xx->cursorPos)+1;
    int epos = DB_Length(xx, 0);
    int fpos,i, j;
    int fseqpos = 0;
    int fseq = -1;
    int bestpos = epos+1;
    int *seqList;
    signed char *quals;
    int clookup[256];
    int type;
    int gel_len = dbi_max_gel_len(DBI(xx),0);
    int last_pos;
    int st;
    float fval = value + 0.5;

    /* Create a base type to index table in the order of acgt*- */
    for (i=0; i<256; i++)
	clookup[i] = 5;
    for (i=0; i<9; i++)
	clookup[(int)"ACGT*acgt"[i]] = i%5;

    quals = (signed char *)xmalloc(gel_len);
    if (NULL == quals)
	return 0;
    memset(quals, -1, gel_len);
    
    seqList = sequencesInRegion(xx,spos, epos);
    fpos = INT_MAX;
    
    last_pos = spos;
    for (i=0; seqList[i]; i++) {
	int seq = seqList[i], e, diff;

	if (fseq != -1 && DB_RelPos(xx, seq) > bestpos)
	    break;

	/* Slide quals buffer left */
	diff = DB_RelPos(xx, seq)-last_pos;
	if (diff > 0) {
	    memmove(quals, &quals[diff], gel_len - diff);
	    memset(&quals[gel_len - diff], -1, diff);
	}

	/* force sequence in memory */
	(void) DBgetSeq(DBI(xx), seq);

	j = DB_RelPos(xx, seq) < spos
	    ? spos - DB_RelPos(xx, seq) + DB_Start(xx, seq)
		: DB_Start(xx, seq);

	e = DB_RelPos(xx, seq) + DB_Length(xx, seq) > fpos
	    ? fpos - DB_RelPos(xx, seq) + DB_Start(xx, seq)
		: DB_Length(xx, seq) + DB_Start(xx, seq);

	/* Add reading to quals array */
	st = DB_Start(xx, seq);
	for (; j < e; j++) {
	    if (DB_Conf(xx, seq)[j] >= fval) {
		type = clookup[(int)DB_Seq(xx, seq)[j]];
		if (quals[j-st] != -1 && quals[j-st] != type) {
		    if (j-st+DB_RelPos(xx, seq) < bestpos) {
			bestpos = j-st+DB_RelPos(xx, seq);
			fseqpos = j-st+1;
			fseq = seq;
		    }
		} else {
		    quals[j-st] = type;
		}
	    }
	}

	last_pos = DB_RelPos(xx, seq);
    }

    xfree(quals);

    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos, fseq);
	edSetBriefSeqBase(xx, -1, -1, 1);
	REDISPLAY(xx);

	if (xx->display_traces)
	    tman_problem_traces(xx, DB_RelPos(xx, fseq)-1 + fseqpos);

	return 1;
    }

    return 0;
}

/*
 * Search backward from the cursor position until we find a place with
 * two or more bases with a quality > 'value' that are not the same base type.
 */
static int findPrevDiscrepancy (EdStruct *xx, int value)
{
    int spos = positionInContig(xx, xx->cursorSeq, xx->cursorPos)+1;
    int epos = 1;
    int fpos,i, j;
    int fseqpos = 0;
    int fseq = -1;
    int bestpos = epos-1;
    int *seqList;
    signed char *quals;
    int clookup[256];
    int type;
    int gel_len = dbi_max_gel_len(DBI(xx),0);
    int last_pos;
    int st;
    float fval = value + 0.5;

    /* Create a base type to index table in the order of acgt*- */
    for (i=0; i<256; i++)
	clookup[i] = 5;
    for (i=0; i<9; i++)
	clookup[(int)"ACGT*acgt"[i]] = i%5;

    quals = (signed char *)xmalloc(gel_len);
    if (NULL == quals)
	return 0;
    memset(quals, -1, gel_len);
    
    seqList = sequencesInRegion(xx,epos, spos);
    fpos = epos;
    
    last_pos = spos;
    for (i=0; seqList[i]; i++); /* Find end of seqlist */
    for (i--; i >= 0 && DB_RelPos(xx, seqList[i]) >= epos; i--) {
	int seq = seqList[i], e, diff;

	if (fseq != -1 && DB_RelPos(xx, seq) + gel_len < bestpos)
	    break;

	/* Slide quals buffer right */
	diff = last_pos - DB_RelPos(xx, seq);
	if (diff > 0) {
	    memmove(&quals[diff], quals, gel_len - diff);
	    memset(quals, -1, diff);
	}

	/* force sequence in memory */
	(void) DBgetSeq(DBI(xx), seq);

	j = DB_RelPos(xx, seq) + DB_Length(xx, seq) > spos
	    ? spos - DB_RelPos(xx, seq) + DB_Start(xx, seq) - 2
		: DB_RelPos(xx, seq) + DB_Length(xx, seq) + DB_Start(xx, seq);

	if (j >= DB_Start(xx, seq) + DB_Length(xx, seq))
	    j = DB_Start(xx, seq) + DB_Length(xx, seq) - 1;

	e = DB_RelPos(xx, seq) < fpos
	    ? fpos - DB_RelPos(xx, seq) + DB_Start(xx, seq)
		: DB_Start(xx, seq);

	/* Add reading to quals array */
	st = DB_Start(xx, seq);

	for (; j >= e; j--) {
	    if (DB_Conf(xx, seq)[j] >= fval) {
		type = clookup[(int)DB_Seq(xx, seq)[j]];
		if (quals[j-st] != -1 && quals[j-st] != type) {
		    if (j-st+DB_RelPos(xx, seq) > bestpos) {
			bestpos = j-st+DB_RelPos(xx, seq);
			fseqpos = j-st;
			fseq = seq;
		    }
		} else {
		    quals[j-st] = type;
		}
	    }
	}

	last_pos = DB_RelPos(xx, seq);
    }

    xfree(quals);
    if (fseq != -1) {
	setCursorPosSeq(xx, fseqpos+1, fseq);
	edSetBriefSeqBase(xx, -1, -1, 1);
	REDISPLAY(xx);

	if (xx->display_traces)
	    tman_problem_traces(xx, DB_RelPos(xx, fseq)-1 + fseqpos+1);

	return 1;
    }

    return 0;
}

/*
 *----------------------------------------------------------------------------
 * Main search interface
 *----------------------------------------------------------------------------
 */
int edDoSearch(EdStruct *xx, int forwards, int strand, char *type, char *value)
{
    int found = 0;
    editor_sort_t group_mode = xx->group_mode;

    /*
     * The "Group templates" option causes sequencesInRegion to return
     * sequences out of positional order. The searching is optimised to assume
     * positional order, so we temporarily disable this while looking for
     * matches.
     */
    xx->group_mode = POSITION;

    if (forwards) {
	if (strcmp(type, "name") == 0)
	    found = findNextGelByName(xx, value);
	
	else if (strcmp(type, "anno") == 0)
	    found = findNextAnno(xx, value);

	else if (strcmp(type, "sequence") == 0) {
	    char *p = strchr(value, '#');
	    int mismatches = 0;
	    int where = 3;

	    if (p) {
		mismatches = atoi(p+1);
		*p = 0;
	    }
	    p = strchr(p+1, '#');
	    if (p) {
		where = atoi(p+1);
	    }

	    found = findNextSequence(xx, value, mismatches, strand, where);
	}

	else if (strcmp(type, "tag") == 0)
	    found = findNextTagByType(xx, value);

	else if (strcmp(type, "position") == 0)
	    found = findPosition(xx, value);

	else if (strcmp(type, "problem") == 0)
	    found = findNextProblem(xx);

	else if (strcmp(type, "quality") == 0)
	    found = findNextQualProb(xx);

	else if (strcmp(type, "edit") == 0)
	    found = findNextEdit(xx);

	else if (strcmp(type, "verifyand") == 0)
	    found = findNextCop(xx, 1);

	else if (strcmp(type, "verifyor") == 0)
	    found = findNextCop(xx, 0);

	else if (strcmp(type, "consquality") == 0)
	    found = findNextConsQual(xx, atoi(value));

	else if (strcmp(type, "discrepancy") == 0)
	    found = findNextDiscrepancy(xx, atoi(value));

	else if (strcmp(type, "consdiscrepancy") == 0)
	    found = findNextConsDiscrep(xx, atoi(value));
    } else {
	/* backwards */
	if (strcmp(type, "name") == 0)
	    found = findPrevGelByName(xx, value);

	else if (strcmp(type, "anno") == 0)
	    found = findPrevAnno(xx, value);

	else if (strcmp(type, "sequence") == 0) {
	    char *p = strchr(value, '#');
	    int mismatches = 0;
	    int where = 3;

	    if (p) {
		mismatches = atoi(p+1);
		*p = 0;
	    }
	    p = strchr(p+1, '#');
	    if (p) {
		where = atoi(p+1);
	    }

	    found = findPrevSequence(xx, value, mismatches, strand, where);
	}

	else if (strcmp(type, "tag") == 0)
	    found = findPrevTagByType(xx, value);

	else if (strcmp(type, "position") == 0)
	    found = findPosition(xx, value);

	else if (strcmp(type, "problem") == 0)
	    found = findPrevProblem(xx);

	else if (strcmp(type, "quality") == 0)
	    found = findPrevQualProb(xx);

	else if (strcmp(type, "edit") == 0)
	    found = findPrevEdit(xx);

	else if (strcmp(type, "verifyand") == 0)
	    found = findPrevCop(xx, 1);

	else if (strcmp(type, "verifyor") == 0)
	    found = findPrevCop(xx, 0);

	else if (strcmp(type, "consquality") == 0)
	    found = findPrevConsQual(xx, atoi(value));

	else if (strcmp(type, "discrepancy") == 0)
	    found = findPrevDiscrepancy(xx, atoi(value));

	else if (strcmp(type, "consdiscrepancy") == 0)
	    found = findPrevConsDiscrep(xx, atoi(value));
    }

    xx->group_mode = group_mode;

    return found;
}


