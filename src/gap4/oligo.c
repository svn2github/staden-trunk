/*
 * File: oligo.c
 * Version:
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: oligo selection module
 *
 * Created: 1991
 * Updated: 06 Nov 92 - changes for distribution
 *	    24 May 95 - major rewrite to generalise suitable for tcl/tk
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "edStructs.h"
#include "IO.h"
#include "edUtils.h"
#include "tagUtils.h"
#include "oligo.h"
#include "misc.h"
#include "reg_exp.h"
#include "xalloc.h"
#include "primlib.h"
#include "dna_utils.h"

/*
 * Compilation modes:
 *
 * VERBOSENESS - Allow for varying degress of verbose output, including debugging information
 */
/* #define VERBOSENESS */

/*
 * Useful #defines
 */
#define FORWARDS  0
#define BACKWARDS 1
extern void messagef(char *format, ...);

/*
 * Parameters for template selection
 */
#ifdef VERBOSENESS
static int verbose = 2;			/* verbose output is required */
char verbosity[10];                     /* space for string form of verbose */
#define verbose_debug (verbose==2 || verbose==3)
#define verbose_panic (verbose==3)
#endif /*VERBOSENESS*/

/*
 *----------------------------------------------------------------------------
 * Add and remove temporary tags to the consensus display
 *----------------------------------------------------------------------------
 */

/*
 * Flag temporary tag as deleted
 */
static void destroy_temporary_tag(EdStruct *xx)
{
    tagStruct *t;

    if (xx->tmp_tag) {
	t = findPreviousTag(xx, 0, xx->tmp_tag);
	(void) _destroy_annotation(DBI(xx), xx, 0/*consensus*/, t,
				   DB_Flags(xx,0) );
	xx->tmp_tag = NULL;
    }

    xx->refresh_flags |= ED_DISP_CONS;
}


/*
 * Create a temporary tag in the consensus to show position of oligo under
 * consideration
 */
static void create_temporary_tag(EdStruct *xx,int pos, int len,
				 int oligo_sense)
{
    tagStruct *tag;
    char *defComment = strdup("*** Temporary Annotation ***\n");
    int sense = (oligo_sense == BACKWARDS);

    destroy_temporary_tag(xx);
    tag = findTagPos(xx, 0, pos);
    xx->tmp_tag = _create_annotation(xx, 0/*consensus*/, pos, len, "OLIG",
				     defComment, tag, sense, DB_Flags(xx,0));


    /*
     * Jiggle about if tag is off screen
     */
    if (xx->displayPos > pos ||	xx->displayPos + xx->displayWidth < pos + len) 
	xx->displayPos = (pos+pos+len-xx->displayWidth)/2;

    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 1);
}



/*
 *----------------------------------------------------------------------------
 * Template checking code. The main function is check_template_for_oligo()
 * which calls each of the appropriate checks.
 *----------------------------------------------------------------------------
 */

/*
 * templates that are in the wrong sense should fail
 * returns 0 - ok,   1 - wrong sense
 */
static int check_sense(EdStruct *xx, int i/*template number*/, int sense)
{
    return (DB_Comp(xx,i) == COMPLEMENTED && sense == FORWARDS ||
	    DB_Comp(xx,i) == UNCOMPLEMENTED && sense == BACKWARDS);
}


/*
 * templates that have their 5' end after the oligo position are not usable
 */
static int check_5prime(EdStruct *xx, int i/*template number*/,
			int sense, int pos, int len)
{
    int relpos = DB_RelPos(xx,i);
    int length = DB_Length(xx,i);

    if (sense == FORWARDS)
	return (relpos > pos);
    else
	return (relpos+length < pos + len);

}


/*
 * Checks position of template
 */
static int check_template_suitability(EdStruct *xx, int i, int sense,
				      int pos, int len)
{
    int relpos = DB_RelPos(xx,i);
    int length = DB_Length(xx,i);
    int near_dist;


    /*
     * Get size information
     */
    get_read_info(DBI_io(xx),
		  i,
		  NULL, 0,	/* clone */
		  NULL, 0,	/* cvector */
		  NULL, 0,	/* subclone */
		  NULL, 0,	/* scvector */
		  NULL, /* length */
		  &near_dist,	/* insert_min */
		  NULL,	/* insert_max */
		  NULL,	/* direction */
		  NULL,	/* strands */
		  NULL,	/* primer */
		  NULL,	/* clone_id */
		  NULL,	/* subclone_id */
		  NULL,	/* cvector_id */
		  NULL		/* scvector_id */
		  );
    
    /* reject ones that are not near the required interval */
    if (sense == FORWARDS)
	return (pos - relpos > (near_dist-xx->sel_oli->ave_read_len));
    else
	return (relpos+length-pos-len > (near_dist-xx->sel_oli->ave_read_len));

}


/*
 * The template number here is actually the number of an existing gel for
 * that template
 *
 * The template:
 *   1. must be in the appropriate sense.
 *   2. must be found "near" the interval required.
 *   3. need not have a sequenced gel over the interval required,
 *	but must be past the 5' end.
 */
static int check_template_for_oligo(EdStruct *xx,int pos, int len,
				    int sense, int i/*template no*/)
{
#define reject_wrong_sense    1
#define reject_mapped_after   2
#define reject_not_close      3


    /* reject ones that have a sense reverse to the one required */
    if (check_sense(xx, i, sense)) {
#ifdef VERBOSENESS
	if (verbose_panic)
	    messagef(" %s rejected because it is in the wrong sense\n",
		     DBgetName(DBI(xx),i));
#endif /*VERBOSENESS*/
	return reject_wrong_sense;
    }


    /* reject ones with 5' end after our oligo position */
    if (check_5prime(xx, i, sense, pos, len)){
#ifdef VERBOSENESS
	if (verbose_panic)
	    messagef(" %s rejected because template starts after oligo "
		     "primer\n", DBgetName(DBI(xx),i));
#endif /*VERBOSENESS*/
	return reject_mapped_after;
    }


    /*
     * Check suitability of template according to position
     */
    if (check_template_suitability(xx, i, sense, pos, len)){
#ifdef VERBOSENESS
	if (verbose_panic)
	    messagef(" %s rejected because template isn't near oligo "
		     "primer position\n", DBgetName(DBI(xx),i));
#endif /*VERBOSENESS*/
	return reject_not_close;
    }


#ifdef VERBOSENESS
    if (verbose_debug)
	messagef(" %s selected\n", DBgetName(DBI(xx),i));
#endif /*VERBOSENESS*/

    return 0;
}



/*
 *----------------------------------------------------------------------------
 *----------------------------------------------------------------------------
 */


static int strncmp_ignorecase(char *a, char *b, int n)
{
    for ( ; n && tolower(*a) == tolower(*b); a++, b++, n--)
	if (*a == '\0') return 0;
    if (!n)
	return 0;
    else
	return tolower(*a) - tolower(*b);
}


/*
 * Check the template name is valid
 *   + that it exists
 *   + that is in the correct sense etc etc
 *
 * Return the chosen gel number for this template or 0 if none found.
 */
static int check_template_name(EdStruct *xx, char *template_name,
			       int pos, int len, int sense) {
    int i, template_index;
    int found_index;
    int template_len;
    int found;

    /*
     * Check template_name exists
     */
    found = found_index = 0;
    template_len = strlen(template_name);

    for (i=1; i<=DBI_gelCount(xx) ; i++) {
	char *name;

	name = DBgetGelName(xx,i);

	if (strncmp_ignorecase(template_name, DBgetGelName(xx,i),
			       template_len) == 0) {
	    found++;
	    found_index = i;
#ifdef VERBOSENESS
	    if (verbose_debug)
		messagef("%s matches %s\n",
			template_name, DBgetName(DBI(xx), i));
#endif
	}
    }


    template_index = 0;

    if (!found)
	messagef("template %s not found\n", template_name);
    else {
	if (found > 1) {
#ifdef VERBOSENESS
	    if (verbose_debug)
		messagef("template %s found, but is not unique\n",
			template_name);
#endif
	} else {
	    if (check_sense(xx, found_index, sense)) {
		messagef("template %s in the wrong sense\n", template_name);
		return 0;
	    }

	    if (check_5prime(xx,found_index,sense,pos,len)) {
		messagef("template %s starts after oligo position\n",
			 template_name);
		return 0;
	    }

	    template_index = found_index;
	}
    }

    return template_index;
}


/*
 * If gel overlaps region in contig,
 * return true
 */
static int gel_ok(EdStruct *xx, int pos, int len, int seq)
{
    return (pos >= DB_RelPos(xx,seq) &&
	    pos+len <= DB_RelPos(xx,seq)+DB_Length(xx,seq));
}





/*
 * Find a gel number for this oligo. Try
 * (A) gel trySeq first, (if valid)
 * (B) then any in the correct sense,
 * (C) otherwise any at all??
 * The position here is a position in the contig
 */
static int find_gel_for_oligo(EdStruct *xx, int pos, int len, int sense,
			      int trySeq) {
    int i;

#ifdef VERBOSENESS
    if (verbose_debug)
	messagef("Trying to find gel for oligo: pos=%d, len=%d, sense=%d\n",
		 pos, len, sense);
#endif

    /**A**/
    if (trySeq > 0) {
#ifdef VERBOSENESS
	if (verbose_debug)
	    messagef("find_gel_for_oligo: Trying gel %s (%d)...\n",
		    DBgetName(DBI(xx),trySeq),trySeq);
#endif

	if ( gel_ok(xx,pos,len,trySeq) ) {
#ifdef VERBOSENESS
	if (verbose_debug)
	    messagef("Using gel %s (%d)...\n", 
		    DBgetName(DBI(xx),trySeq),trySeq);
#endif
	    return trySeq;
	}
    }

    /**B**/
#ifdef VERBOSENESS
    if (verbose_debug)
	message("find_gel_for_oligo: Trying gels in correct sense\n");
#endif /*VERBOSENESS*/

    for (i=1; i<=DBI_gelCount(xx); i++) {
	if (DB_Comp(xx,i) == COMPLEMENTED && sense == BACKWARDS ||
	    DB_Comp(xx,i) == UNCOMPLEMENTED && sense == FORWARDS) {
	    if (gel_ok(xx,pos,len,i)) {
#ifdef VERBOSENESS
		if (verbose_debug)
		    messagef("Using gel %s (%d)...\n",
			     DBgetName(DBI(xx),trySeq),trySeq);
#endif /*VERBOSENESS*/
		return i;
	    }
	}
    }

    /**C**/
#ifdef VERBOSENESS
    if (verbose_debug) message("find_gel_for_oligo: Trying any gel\n");
#endif /*VERBOSENESS*/
    for (i=1; i<=DBI_gelCount(xx); i++) {
	    if (gel_ok(xx,pos,len,i)) {
#ifdef VERBOSENESS
		if (verbose_debug)
		    messagef("Using gel %s (%d)...\n",
			     DBgetName(DBI(xx),trySeq),trySeq);
#endif /*VERBOSENESS*/
		return i;
	    }
    }

#ifdef VERBOSENESS
    if (verbose_debug) message("find_gel_for_oligo: Failed to find any suitable gel\n");
#endif /*VERBOSENESS*/
    return 0;

}	 





static char *generate_oligo_comment(EdStruct *xx, int oligo, char *tname)
{
    char s[200];
    char seq[100];
    char *c;
    double tm, gc, qual;

    int pos,len;
#ifdef VERBOSENESS
    if (verbose_debug) message("creating comment for oligo:\n");
#endif /*VERBOSENESS*/

    pos = xx->sel_oli->state->primers[oligo].start;
    len = xx->sel_oli->state->primers[oligo].length;
    tm  = xx->sel_oli->state->primers[oligo].temp;
    gc  = xx->sel_oli->state->primers[oligo].gc_content;
    qual= xx->sel_oli->state->primers[oligo].quality;
    strncpy(seq,&xx->sel_oli->consensus[pos],len);
    seq[len]='\0';

    
    sprintf(s,"serial#=\ntemplate=%s\nsequence=%s\nTm=%.2f\n"
	    "GC%%=%.2f\nquality=%.2f\n",
	    tname, seq, tm, gc, qual);

#ifdef VERBOSENESS
    if (verbose_debug) messagef("(%s)\n",s);
#endif /*VERBOSENESS*/

    c = TAG_MALLOC(strlen(s)+1);
    strcpy(c,s);
    return c;
}





/*
 * This routine creates a new oligo tag, prior to leaving the
 * oligo selection window
 */
static int create_new_oligo_tag(EdStruct *xx, int oligo, int pos, int len,
				int sense, char *template_name) {
    tagStruct *tag;
    int seq;
    int npos;

    /* Find a suitable reading number */
    seq = check_template_name(xx, template_name, pos, len, sense);
    seq = find_gel_for_oligo(xx, pos, len, sense, seq);

    if (! seq) {
	messagef("NO SUITABLE GEL FOR THIS OLIGO TAG POSITION %d LENGTH %d\n",
		 pos, len);
	return 1;
    }

    openUndo(DBI(xx));
    pos = pos - DB_RelPos(xx, seq) + 1 + DB_Start(xx, seq);
    npos = normalisePos2(xx, seq, pos, len);
    tag = findTagPos(xx,seq,npos);
    U_create_annotation(xx,
			seq,
			npos,
			len,
			"OLIG",
			generate_oligo_comment(xx, oligo, template_name),
			tag,
			normaliseSense(xx,seq,(sense==BACKWARDS)));
    U_adjust_cursor(xx, 0);
    closeUndo(xx, DBI(xx));

    return 0;
}



/*
 *----------------------------------------------------------------------------
 *----------------------------------------------------------------------------
 */


/*
 * Once an oligo has been found for the consensus at position pos, length len,
 * search for a suitable template within the contig.
 *
 */
static int *find_templates_for_oligo(EdStruct *xx, int pos, int len, int sense)
{
    static int *templateList = NULL;
    int i, count;

    if (NULL == (templateList = (int *)xmalloc((DBI_gelCount(xx)+1)*sizeof(int))))
	return NULL;

#ifdef VERBOSENESS
    if (verbose_debug) {
	message("Finding template for oligo:\n");
	messagef("position = %d, length = %d, forward sense=%d\n",
		 pos,len,sense);
    }
#endif /*VERBOSENESS*/

    count = 0;
    if (sense == BACKWARDS) {
	for(i=1; i<=DBI_gelCount(xx); i++) {
	    char *name;
	    name = DBgetGelName(xx,i);

	    if (!check_template_for_oligo(xx, pos, len, sense, i))
		templateList[count++] = i;
	    
	}
    } else {
	int ind;
	ind=posToIndex(xx,pos); /* optimise a bit */
	if(!ind) ind=DBI_gelCount(xx);
	for(;ind>0;ind--) {
	    char *name;
	    i = DBI_order(xx)[ind];
	    name = DBgetGelName(xx,i);
	    
	    if (!check_template_for_oligo(xx, pos, len, sense, i))
		templateList[count++] = i;
	    
	}
    }
    
    templateList[count] = 0;

    return templateList;
}




static void get_temp_name(char *name, EdStruct *xx, int i)
/*
 * Get the template name for gel number 'i'
 */
{
    GTemplates t;
    GReadings r;

    gel_read(DBI_io(xx), DB_Number(xx,i), r);
    template_read(DBI_io(xx), r.template, t);

    if (!t.name || 0 != TextRead(DBI_io(xx), t.name, name, DB_NAMELEN)) {
	(void) strcpy(name, DBgetGelName(xx,i) );
    }
}




static int score_template(EdStruct *xx, int seq)
/*
 * Score this template
 */
{
    return 1; /* this should force the first one to be chosen */
}



/*
 * Pick a default template from the list of available templates
 */
static char *get_default_template(EdStruct *xx, int *templateList)
{
    int template_index;
    static char template_name[DB_NAMELEN+1];

    *template_name = 0;

    if (templateList[0]) {
	int i;
	int score, high_score;
	high_score = 0;
	/* search */
	for (i=0; templateList[i]; i++) {
	    score = score_template(xx, templateList[i]);
	    if (high_score < score) {
		high_score = score;
		template_index = templateList[i];
	    }
	}

	/* set template name */
	get_temp_name(template_name, xx, template_index);
	template_name[DB_NAMELEN] = 0;
    }

    return template_name;
}




/*
 * Builds the return string for nextOligo. See nextOligo() for the format
 * of this string.
 */
static char *nextOligo_return(EdStruct *xx, int *list)
{
    int i;
    char *ret, *ptr;

    /* Alloc memory needed */
    for (i = 0; list[i]; i++)
	;

    if (NULL == (ret = (char *)xmalloc((DB_NAMELEN+1)*(i+1)+1))) {
	return NULL;
    }
    ptr = ret;

    /* Add the default */
    sprintf(ptr, "%s ", get_default_template(xx, list));
    ptr[DB_NAMELEN] = 0;
    ptr += strlen(ptr);

    /* Now add the list of templates to choose from */
    for (i = 0; i < list[i]; i++) {

#ifdef VERBOSENESS
	if (verbose_debug)
	    messagef("Adding %s to return list\n",DBgetName(DBI(xx),list[i]));
#endif /*VERBOSENESS*/

	/*
	 * Prepare clone name
	 */
	get_temp_name(ptr, xx, list[i]);
	ptr[DB_NAMELEN] = 0;
	ptr += strlen(ptr);
	*ptr++ = ' ';
    }

    *ptr++ = 0;

    return ret;
}






/*
 * Display important oligo details
 *   a. sequence
 *   b. position, length
 *   c. primer score, gc, tm
 */
static void display_oligo_details(EdStruct *xx, int oligo)
{
    char seq[100];
    int pos,len;
    int i;
    char *a;
    int consensus_start, consensus_end;
    char statline[1024];

    pos = xx->sel_oli->state->primers[oligo].start;
    len = xx->sel_oli->state->primers[oligo].length;

    /*
     * Get consensus sequence with context
     *
     * "...conSENSUS_SEQUEnce..."
     */
    a = seq;
    *a++ = '.'; *a++ = '.'; *a++ = '.';

    for (i = (pos>2)?2:pos; i; i--)
	*a++ = tolower(xx->sel_oli->consensus[pos-i]);
    strncpy(a,&xx->sel_oli->consensus[pos],len); a+=len;
    for (i=0;i<2&&xx->sel_oli->consensus[pos+len+i];i++)
	*a++ = tolower(xx->sel_oli->consensus[pos+len+i]);

    *a++ = '.'; *a++ = '.'; *a++ = '.';
    *a++ = '\0';


    /*
     * Positions in contig
     */
    if (xx->sel_oli->sense == BACKWARDS) {
	consensus_start = xx->sel_oli->r -
	    (xx->sel_oli->state->primers[oligo].start + 
	     xx->sel_oli->state->primers[oligo].length-1);
	consensus_end   = consensus_start + len - 1;
    } else {
	consensus_start = xx->sel_oli->l +
	    xx->sel_oli->state->primers[oligo].start;
	consensus_end   = consensus_start + len - 1;
    }

    vfuncgroup(4, "Oligo details");
    vmessage("Oligo: %s\n",seq);
    vmessage("Primer # %2d                                      PRIMER-SELF\n"
	     "5' end   3' end    length  Score G+C(%%)  Tm      3'  Internal\n"
	     "%6d   %6d      %4d    %4.1f  %4.1f  %4.1f    %4.1f    %4.1f\n",
	     oligo+1,
	     consensus_start,
	     consensus_end,
	     len,
	     xx->sel_oli->state->primers[oligo].quality,
	     xx->sel_oli->state->primers[oligo].gc_content,
	     xx->sel_oli->state->primers[oligo].temp,
	     xx->sel_oli->state->primers[oligo].self_any/100.0,
	     xx->sel_oli->state->primers[oligo].self_end/100.0);

    sprintf(statline, "Oligo %s, Len %d, Score %4.1f, Tm %4.1f, "
	    "GC %4.1f%%C, self-3' %4.1f, self-internal %4.1f",
	    seq, len,
	    xx->sel_oli->state->primers[oligo].quality,
	    xx->sel_oli->state->primers[oligo].temp,
	    xx->sel_oli->state->primers[oligo].gc_content,
	    xx->sel_oli->state->primers[oligo].self_any/100.0,
	    xx->sel_oli->state->primers[oligo].self_end/100.0);
    tk_update_brief_line(xx, statline);
}




/*
 *----------------------------------------------------------------------------
 * Publically callable functions below here
 *----------------------------------------------------------------------------
 */

/*
 * Find suitable oligos using Primer3 with current parameter settings.
 * Return number of oligos found (or -1 for error).
 */
int edSelectOligoGenerate(EdStruct *xx, int sense, int bkwd_width,
			  int fwd_width, int avg_read_len, char *primer_defs) {
    int seq = xx->cursorSeq;
    int pos = xx->cursorPos;
    int position; /* in contig */
    int contigLength; /* length of contig */
    int consensusLength;
    int i, j;
    select_oligo_t *so;
    primlib_args *args;
    int ok;

    if (xx->editorState == StateDown)
	return -1;

#ifdef VERBOSENESS
    if (verbose_debug) message("Finding oligos...\n");
#endif /*VERBOSENESS*/

    position =  positionInContig(xx,seq,pos);
    contigLength = DB_Length(xx,0/*consensus*/);

    /*
     * Allocate a select_oligo_t structure to store our semi-persistent
     * results. These stay allocated until select oligos is quitted.
     */
    if (NULL == (so = (select_oligo_t *)xmalloc(sizeof(select_oligo_t)))) {
	bell();
	return -1;
    }

    so->consensus = NULL;
    so->opos = NULL;
    so->state = NULL;
    so->opos_start = NULL;
    so->opos_end = NULL;

    /* Create an initialise primlib state */
    so->state = primlib_create();
    if (NULL == (args = primlib_str2args(primer_defs)))
	return -1;

    primlib_set_args(so->state, args);
    free(args);

    /*
     * Conceptually we select off consensus.
     * Determine consensus around point...
     *    the oligo will be selected from this region
     */
    if (sense == FORWARDS) {
	so->l = (position > bkwd_width)
	    ? position - bkwd_width : 1;
	so->r = (position + fwd_width < contigLength)
	    ? position + fwd_width : contigLength;
    } else {
	so->l = (position > fwd_width)
	    ? position - fwd_width : 1;
	so->r = (position + bkwd_width < contigLength)
	    ? position + bkwd_width : contigLength;
    }

    /*
     * Allocated and calculate the consensus sequence.
     */
    consensusLength = (so->r - so->l)+1;
    /* allocate space for consensus */
    if (NULL == (so->consensus = (char *)xmalloc(consensusLength+1)) ||
	NULL == (so->opos = (int *)xmalloc((consensusLength+1)*sizeof(int)))) {
	bell();
	xfree(so);
	return -1;
    }
    so->consensus[consensusLength] = 0;
    DBcalcConsensus(xx, so->l, consensusLength, so->consensus,
		    NULL, BOTH_STRANDS);
    /*
     * Complement if necessary
     */
    if (sense == BACKWARDS) {
	complement_seq(so->consensus, consensusLength);
    }

    for (j = i = 0; i < consensusLength; i++) {
	so->opos[i] = j;
	if (so->consensus[i] != '*') {
	    so->consensus[j] = so->consensus[i];
	    j++;
	}
    }
    so->consensus[j] = 0;

#ifdef VERBOSENESS
    if (verbose_debug) {
	messagef("Cursor position = %d\n",pos);
	messagef("Seqence = %s (%d)\n",DBgetName(DBI(xx),seq),seq);
	messagef("Consensus for region %d-%d = (%s)\n",
		 so->l, so->r, so->consensus);
	messagef("Sense = %s\n", (sense==FORWARDS)?"forward":"reverse");
    }
#endif /*VERBOSENESS*/

    ok = primlib_choose(so->state, so->consensus);
    if (ok == -1 || so->state->nprimers == 0)
	return 0;

#ifdef VERBOSENESS
    if (verbose) messagef("primlib_choose returned with status %d\n", ok);
#endif /*VERBOSENESS*/


#ifdef VERBOSENESS
    if (verbose) messagef("%d suitable oligos found\n", so->state->nprimers);
#endif /*VERBOSENESS*/

    so->opos_start = (int *)xcalloc(so->state->nprimers, sizeof(int));
    so->opos_end = (int *)xcalloc(so->state->nprimers, sizeof(int));
    for (i = 0; i < so->state->nprimers; i++) {
	int st = so->state->primers[i].start, st2 = st;
	int en = st + so->state->primers[i].length-1, en2 = en;
	for (j = st; j < consensusLength; j++) {
	    if (so->opos[j] == st)
		st2 = j;
	    if (so->opos[j] == en)
		en2 = j;
	}
	so->opos_start[i] = st2;
	so->opos_end[i] = en2;
    }

    so->curr_oligo = -1;
    so->ave_read_len = avg_read_len;
    so->sense = sense;

    if (xx->sel_oli) {
	xfree(xx->sel_oli->consensus);
	xfree(xx->sel_oli->opos);
	xfree(xx->sel_oli);
    }
    xx->sel_oli = so;

    return so->state->nprimers;
}



/*
 * We cycle through the oligo list
 * curr-oligo gives the current oligo entry under consideration.
 *
 * We return a list of strings in the form:
 * default_tname tname ...
 *
 * The string returned is malloced and so is expected to be freed by the
 * calling routine.
 */
char *edSelectOligoSwitch(EdStruct *xx) {
    char *ret;
    int *templateList;
    int sense = xx->sel_oli->sense;
    int i = xx->sel_oli->curr_oligo;

    /*
     * Print out information on oligos
     */
#ifdef VERBOSENESS
    if (verbose_debug) messagef("************** %d ************\n",i);
    if (verbose_debug) messagef("stp %d, len %d, qual %f,  gc %f, tm %f\n",
				xx->sel_oli->state->primers[i].start,
				xx->sel_oli->state->primers[i].length,
				xx->sel_oli->state->primers[i].quality,
				xx->sel_oli->state->primers[i].gc_content,
				xx->sel_oli->state->primers[i].temp);
#endif /*VERBOSENESS*/

    /*
     * Convert position returned from oligo selection to
     * position in contig
     */
    if (sense == BACKWARDS) {
	templateList =
	    find_templates_for_oligo(xx, 
				     xx->sel_oli->r -
				     xx->sel_oli->opos_end[i],
				     xx->sel_oli->opos_end[i] -
				     xx->sel_oli->opos_start[i] + 1,
				     sense);
	create_temporary_tag(xx,
			     xx->sel_oli->r -
			     xx->sel_oli->opos_end[i],
			     xx->sel_oli->opos_end[i] -
			     xx->sel_oli->opos_start[i] + 1,
			     sense);
    } else {
	templateList =
	    find_templates_for_oligo(xx, 
				     xx->sel_oli->l +
				     xx->sel_oli->opos_start[i],
				     xx->sel_oli->opos_end[i] -
				     xx->sel_oli->opos_start[i] + 1,
				     sense);
	create_temporary_tag(xx,
			     xx->sel_oli->l +
			     xx->sel_oli->opos_start[i],
			     xx->sel_oli->opos_end[i] -
			     xx->sel_oli->opos_start[i] + 1,
			     sense);
    }

    display_oligo_details(xx,i);

    if (NULL == templateList) {
	return NULL;
    }

    ret = nextOligo_return(xx, templateList);
    xfree(templateList);

    return ret;
}

char *edSelectOligoNext(EdStruct *xx)
{
    if (xx->editorState == StateDown)
	return NULL;

    if (xx->sel_oli->curr_oligo+1 >= xx->sel_oli->state->nprimers) {
	return NULL;
    } else {
	++xx->sel_oli->curr_oligo;
    }

    return edSelectOligoSwitch(xx);
}

char *edSelectOligoPrev(EdStruct *xx)
{
    if (xx->editorState == StateDown)
	return NULL;

    if (xx->sel_oli->curr_oligo == 0) {
	return NULL;
    } else {
	--xx->sel_oli->curr_oligo;
    }

    return edSelectOligoSwitch(xx);
}



/*
 * Select the current oligo.
 * Creats a tag for it and everything else...
 *
 * Returns a status line containing the template and sequence
 */
char *edSelectOligoAccept(EdStruct *xx, char *template_name) {
    static char status[100];
    int i = xx->sel_oli->curr_oligo, r;
    int st, en;

    if (xx->editorState == StateDown)
	return NULL;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	bell();
	return NULL;
    }
    
    st = xx->sel_oli->opos_start[i];
    en = xx->sel_oli->opos_end[i];

    if (xx->sel_oli->sense == BACKWARDS) {
	r = create_new_oligo_tag(xx, i,
				 xx->sel_oli->r - en,
				 en - st + 1, 
				 xx->sel_oli->sense,
				 template_name);
    } else {
	r = create_new_oligo_tag(xx,i,
				 xx->sel_oli->l + st, 
				 en - st + 1,
				 xx->sel_oli->sense,
				 template_name);
    }

    if (r)
	bell();

    redisplaySequences(xx, 1);

    st = xx->sel_oli->state->primers[i].start;
    en = st + xx->sel_oli->state->primers[i].length - 1;
    sprintf(status, "%s %.*s",
	    *template_name ? template_name : "<None>",
	    en - st + 1, &xx->sel_oli->consensus[st]);

    return status;
}



/*
 * Frees up the temporary structures used by the oligo selection code.
 */
void edSelectOligoQuit(EdStruct *xx) {
    destroy_temporary_tag(xx);

    if (xx->editorState == StateDown)
	return;

    if (xx->sel_oli) {
	xfree(xx->sel_oli->consensus);
	xfree(xx->sel_oli->opos);
	if (xx->sel_oli->opos_start)
	    xfree(xx->sel_oli->opos_start);
	if (xx->sel_oli->opos_end)
	    xfree(xx->sel_oli->opos_end);
	if (xx->sel_oli->state)
	    primlib_destroy(xx->sel_oli->state);
	xfree(xx->sel_oli);
	xx->sel_oli = NULL;
    }

    redisplaySequences(xx, 1);
}



/*************************************************************
 * New features - unused at present
 *************************************************************/
#include "array.h"
#include "gap-dbstruct.h"

#define REL_UNKNOWN 0
#define REL_EQ      1
#define REL_GT      2
#define REL_GE      3
#define REL_LT      4
#define REL_LE      5
#define REL_NE      6
/*#define REL_UNKNOWN	7*/


#if 0
static char *REL_STR[8] = { "?","",">",">=","<","<=","<>","<=>" };
#define rel2str(A) (REL_STR[A])
#endif

typedef struct {
	int id;
	int rel[2];		/* relationship of ends */
	int end[2];		/* mapping of ends */
	int sense;		/* direction of insert is 3'->5' ? */
	int strands;
	int insert_min;
	int insert_max;
	int flag;		/* in use */
} Templates;


#if 0
static void squid(EdStruct *xx, int pos, int len, int sense)
{
    Array templates;
    int Ntemplates;
    int i,j;

    Ntemplates = 0;
    templates = ArrayCreate(sizeof(Templates),Ntemplates);

    for(i=1; i<=DBI_gelCount(xx); i++) {
	int imin, imax, tnum, strand, strands, primer;
	int hit;
	int gstart, gend;
	int gorient;


	/* Gel information */
	gstart = DB_RelPos(xx,i);
	gend = gstart + DB_Length(xx,i) - 1;
	gorient = DB_Comp(xx,i);
	


	/* Determine Template ID */
	get_read_info(DBI_io(xx),
		      DB_Number(xx,i),
		      NULL, 0,	/* clone */
		      NULL, 0,	/* cvector */
		      NULL, 0,	/* subclone */
		      NULL, 0,	/* scvector */
		      NULL,	/* length */
		      &imin,	/* insert_min */
		      &imax,	/* insert_max */
		      &strand,	/* direction */
		      &strands,	/* strands */
		      &primer,	/* primer */
		      NULL,	/* clone_id */
		      &tnum,	/* subclone_id */
		      NULL,	/* cvector_id */
		      NULL	/* scvector_id */
		      );
	
	/* Is template new ? */
	
	hit = -1;
	for(j=0;j<Ntemplates && hit<0 ;j++) {
	    if(arr(Templates,templates,j).id == tnum) hit = j;
	}

	/* NO */
	if (hit<0) {
	    hit = Ntemplates++;
	    (void) ArrayRef(templates,hit);

	    arr(Templates,templates,hit).id = tnum;
	    arr(Templates,templates,hit).strands = strands;
	    arr(Templates,templates,hit).insert_min = imin;
	    arr(Templates,templates,hit).insert_max = imax;
	    arr(Templates,templates,hit).flag = 0;

	    /*
	     * Template.sense
	     * This is a little complicated to explain...
	     * We need to know where the forward and reverse primers are
	     * relative to each other. A value of 0 means forward primer is to
	     * the left of reverse primer.
	     */
	    arr(Templates,templates,hit).sense =
		 (gorient==UNCOMPLEMENTED && strand==GAP_STRAND_REVERSE ||
		  gorient==COMPLEMENTED && strand==GAP_STRAND_FORWARD);

	    /*
	     * Let's pretend we don't know anything about this for the time being
	     */
	    arr(Templates,templates,hit).rel[0] = REL_UNKNOWN;
	    arr(Templates,templates,hit).end[0] = 0;
	    arr(Templates,templates,hit).rel[1] = REL_UNKNOWN;
	    arr(Templates,templates,hit).end[1] = 0;
	}



	/*
	 * Check for consistency between primer and direction(strand) of read
	 */
	if (primer==GAP_PRIMER_UNKNOWN ||
	    primer==GAP_PRIMER_REVERSE && strand==GAP_STRAND_FORWARD ||
	    primer==GAP_PRIMER_FORWARD && strand==GAP_STRAND_REVERSE) {
	    printf(primer == GAP_PRIMER_UNKNOWN
		   ? "Primer direction unknown for gel reading %d\n"
		   : "Inconsistant sense information for primer "
		     "(gel reading %d)\n", DB_Number(xx,i));
	    continue;
	}


	/*
	 * Modify number of strands if we are given contradictory information
	 */

	if (primer==GAP_PRIMER_REVERSE)
	    arr(Templates,templates,hit).strands = 2;



	/*
	 * Let's add information about this reading to the mapping
	 * information.
	 */

	/*
	 * Primer	gorient Comment
	 * fwd or rev	0	L==gstart,R>=gend
	 * fwd or rev	1	L<=gstart,R==gend
	 * 2(cust)	any	L<=gstart,R>=gend
	 */


#define yip(A,B,C) do { printf(">>>>>> %d should be %s %d but it isn't (gel reading %d)\n", (A),rel2str(B),(C),DB_Number(xx,i)); } while(0)
#define yap(A) printf(">>>>>> %d is an invalid relationship for this end\n",(A))

	/*
	 * Set left end limits
	 */
	if ((primer==GAP_PRIMER_FORWARD ||
	     primer==GAP_PRIMER_REVERSE) &&
	    gorient==UNCOMPLEMENTED) {
	    /* Lend == start */
	    int Lrel = arr(Templates,templates,hit).rel[0];
	    int Lend = arr(Templates,templates,hit).end[0];
	    switch(Lrel) {
	    case REL_EQ:
		if (gstart<Lend)
		    arr(Templates,templates,hit).end[0] = gstart;
		break;
	    case REL_LE:
		if (gstart>Lend) yip(gstart,REL_LE,Lend);
		arr(Templates,templates,hit).rel[0] = REL_EQ;
		arr(Templates,templates,hit).end[0] = gstart;
		break;
	    case REL_UNKNOWN:
		arr(Templates,templates,hit).rel[0] = REL_EQ;
		arr(Templates,templates,hit).end[0] = gstart;
		break;
	    default:
		yap(Lrel);
		break;
	    }
	} else {
	    /* Lend <= start */
	    int Lrel = arr(Templates,templates,hit).rel[0];
	    int Lend = arr(Templates,templates,hit).end[0];
	    switch(Lrel) {
	    case REL_EQ:
		if (gstart<Lend) yip(gstart,REL_GE,Lend);
		break;
	    case REL_LE:
		if (gstart<Lend) {
		    arr(Templates,templates,hit).rel[0] = REL_LE;
		    arr(Templates,templates,hit).end[0] = gstart;
		}
		break;
	    case REL_UNKNOWN:
		arr(Templates,templates,hit).rel[0] = REL_LE;
		arr(Templates,templates,hit).end[0] = gstart;
		break;
	    default:	
		yap(Lrel);
		break;
	    }
	}




	/*
	 * Set Right end limits
	 */
	if ((primer==GAP_PRIMER_FORWARD ||
	     primer==GAP_PRIMER_REVERSE) &&
	    gorient==COMPLEMENTED) {
	    /* Rend == end */
	    int Rrel = arr(Templates,templates,hit).rel[1];
	    int Rend = arr(Templates,templates,hit).end[1];
	    switch(Rrel) {
	    case REL_EQ:
		if (gend>Rend)
		    arr(Templates,templates,hit).end[1] = gend;
		break;
	    case REL_GE:
		if (gend<Rend) yip(gend,REL_GE,Rend);
		arr(Templates,templates,hit).rel[1] = REL_EQ;
		arr(Templates,templates,hit).end[1] = gend;
		break;
	    case REL_UNKNOWN:
		arr(Templates,templates,hit).rel[1] = REL_EQ;
		arr(Templates,templates,hit).end[1] = gend;
		break;
	    default:
		yap(Rrel);
		break;
	    }
	} else {
	    /* Rend <= end */
	    int Rrel = arr(Templates,templates,hit).rel[1];
	    int Rend = arr(Templates,templates,hit).end[1];
	    switch(Rrel) {
	    case REL_EQ:
		if (gend<Rend) yip(gend,REL_GE,Rend);
		break;
	    case REL_GE:
		if (gend>Rend) {
		    arr(Templates,templates,hit).rel[1] = REL_GE;
		    arr(Templates,templates,hit).end[1] = gend;
		}
		break;
	    case REL_UNKNOWN:
		arr(Templates,templates,hit).rel[1] = REL_GE;
		arr(Templates,templates,hit).end[1] = gend;
		break;
	    default:
		yap(Rrel);
		break;
	    }
	}
    }



    /*
     * Merge in insert length information
     */
    for(i=0;i<Ntemplates;i++) {
	if (arr(Templates,templates,i).rel[0] == REL_EQ &&
	    arr(Templates,templates,i).rel[1] == REL_EQ) {
	    /*
	     * Confirm our mapping
	     */
	    int len;
	    len = arr(Templates,templates,i).end[1] -
		arr(Templates,templates,i).end[0]+1;
	    if (len < arr(Templates,templates,i).insert_min ||
		len > arr(Templates,templates,i).insert_max) {
		printf("?????? Template %d - length==%d, %d<=insert<=%d\n",
		       i,len,
		       arr(Templates,templates,i).insert_min,
		       arr(Templates,templates,i).insert_max);
	    }

	    /* set */
	    arr(Templates,templates,i).insert_min = len;
	    arr(Templates,templates,i).insert_max = len;

	} else  {
	    int min_len;
	    min_len = arr(Templates,templates,i).end[1] -
		arr(Templates,templates,i).end[0]+1;

	    if (min_len > arr(Templates,templates,i).insert_min)
		arr(Templates,templates,i).insert_min = min_len;

	    if (min_len > arr(Templates,templates,i).insert_max) {
		printf("?????? Template %d - length>=%d, %d<=insert<=%d\n",
		       i,
		       arr(Templates,templates,i).end[1]-
		       arr(Templates,templates,i).end[0]+1,
		       arr(Templates,templates,i).insert_min,
		       arr(Templates,templates,i).insert_max);
		
		arr(Templates,templates,i).insert_max = min_len;
	    }
	    

	    if (arr(Templates,templates,i).rel[0] == REL_EQ) {
		arr(Templates,templates,i).end[1] =
		    arr(Templates,templates,i).end[0] +
			arr(Templates,templates,i).insert_min - 1;
	    } else if (arr(Templates,templates,i).rel[1] == REL_EQ) {
		arr(Templates,templates,i).end[0] =
		    arr(Templates,templates,i).end[1] -
			arr(Templates,templates,i).insert_min + 1;
	    }
	}
    }
    

#ifdef VERBOSENESS
    if (verbose_debug) {
    /* dump out mapping information so far */
    printf("Template Mapping ... %d templates\n",Ntemplates);
    printf("%-20.20s %2s%6s %2s%6s %7s %7s insert %6s %6s\n",
	   "Template Name", " ", "5'", " ", "3'", "(sense)", "strands", "min", "max");

    for(i=0;i<Ntemplates;i++) {
	char name[20];

	(void) get_subclone_info(DBI_io(xx),
				 arr(Templates,templates,i).id,
				 NULL, 0,
				 NULL, 0,
				 name, sizeof(name),
				 NULL, 0,
				 NULL, NULL, NULL, NULL, NULL, NULL);

	
	/* template-name 5'-end 3'-end (sense) strands insert-min insert-max */

	printf("%-20.20s %2s%6d %2s%6d %7d %7d        %6d %6d\n",
	       name,
	       rel2str(arr(Templates,templates,i).rel[0]),
	       arr(Templates,templates,i).end[0],
	       rel2str(arr(Templates,templates,i).rel[1]),
	       arr(Templates,templates,i).end[1],
	       arr(Templates,templates,i).sense,
	       arr(Templates,templates,i).strands,
	       arr(Templates,templates,i).insert_min,
	       arr(Templates,templates,i).insert_max);

    }

    printf("\n\n");
    }
#endif /* VERBOSENESS */

    ArrayDestroy(templates);


}
#endif

/*
 * Creates a temporary annotation view, much like the suggest primers and the
 * pick pcr primer functions.
 * It bypasses the undo mechanism and will not be saved.
 *
 * NB: Currently this ignores many of the arguments.
 */
void createTmpAnnotation(EdStruct *xx, int seq, int pos, int len,
			 char *type, char *text, int strand)
{
    tagStruct *tag;
    tag = findTagPos(xx, seq, pos);
    if (xx->tmp_tag)
	destroy_temporary_tag(xx);
    create_temporary_tag(xx, pos, len, strand);
}


