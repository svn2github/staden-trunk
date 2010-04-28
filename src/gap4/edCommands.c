#include <staden_config.h>

#include <ctype.h>

#include "edUtils.h"
#include "align.h"
#include "xalloc.h"
#include "select.h"
#include "tman_interface.h"
#include "gap_globals.h"

/*
 *----------------------------------------------------------------------------
 * Dump Contig code
 *---------------------------------------------------------------------------
 */

/*
 * get part of a sequence from its `pos' base for `width' bases
 * Bases number from 0?
 */
static void dumpSequence(EdStruct *xx, int seq, int pos, int width, char *str)
{
    char *src;
    int length = DB_Length(xx,seq);
    int i,j;
    
    src = DBgetSeq(DBI(xx),seq);
    
    /* Lefthand cut off */
    if (pos<0) {
	i = (width<-pos)?width:-pos;
	getLCut(xx,seq, -pos, i, str);
	for(j=0;j<i;j++) if (isupper(str[j])) str[j] = tolower(str[j]);
    } else
	i=0;
    
    /*copy sequence*/
    for (; i<width && (pos+i)<length; i++) {
	str[i]=src[pos+i]; 
    }
    
    /* Righthand cut off */
    if (i<width) {
	getRCut(xx,seq, pos+i-length, width-i, &str[i]);
	for(j=i;j<width;j++) if (isupper(str[j])) str[j] = tolower(str[j]);
    }
    
    str[width]='\0';
}


/* If the tcl dialogue changes, this also needs changing */
#define MAX_LINE_LENGTH 1000

/*
 * Print out a section
 */
static void dumpLine(EdStruct *xx, FILE *fp, int pos, int width, int nwidth)
{
    int *seqList;
    int i;
    char spare[MAX_LINE_LENGTH+25];
    char con[MAX_LINE_LENGTH+5], *consensus = &con[2];
    int displayHeight;
    int namelen = nwidth + DB_GELNOLEN + 1;
    int same_cons;
    
    displayHeight = linesOnScreen(xx,pos,width);
    seqList = sequencesOnScreen(xx,pos, width);
    
    
    if (xx->rulerDisplayed) {
	char *k = spare;
	int j,lower,times;

	lower = (pos - pos%10);
	times = width/10 + 2;

	for (j=0; j<times; j++,k+=10,lower+=10)
	    sprintf(k,"%10d",lower);
	fprintf(fp,"%*.*s   %-*.*s\n",
		namelen, namelen, " ",
		width,width,&spare[9+pos%10]);
    }
    DBcalcConsensus(xx,pos-2,width+4,con,NULL,BOTH_STRANDS);
    
    for (i=0 ; i < displayHeight ; i++ ) {
	char * ptr;
	
	if (DB_Flags(xx,seqList[i]) & DB_FLAG_SELECTED)
	    fprintf(fp,"%*.*s * ",namelen,namelen, DBgetName(DBI(xx),
							     seqList[i]));
	else
	    fprintf(fp,"%*.*s   ",namelen,namelen, DBgetName(DBI(xx),
							     seqList[i]));
	
	if (seqList[i]==0){
	    ptr = consensus;
	}else{
	    dumpSequence(xx,seqList[i],pos-DB_RelPos(xx,seqList[i]),width,
			 spare);
	    ptr = spare;
	}
	if (xx->showDifferences) {
	    int j;

	    for (j=0;j<width;j++) if (spare[j]==consensus[j])
		spare[j]='.';
	}
	fprintf(fp,"%*.*s\n",width,width,ptr);
    }

    same_cons = 0;

    /* Status lines */
    tk_redisplaySeqStatusCompute(xx, pos, width);
    for (i = 0; i < xx->status_depth; i++) {
	fprintf(fp, "%.*s   %.*s\n",
		namelen, xx->status_lines[i].name+1,
		width, xx->status_lines[i].line);
    }

    fprintf(fp,"\n");
}


static void dumpRegion(EdStruct *xx, FILE *fp, int start, int end, int width,
		       int nwidth)
{
    for(;start<=end;start+=width)
	dumpLine(xx, fp, start, (end-start+1<width)?end-start+1:width,
		 nwidth);
}


void dumpContig(EdStruct *xx, char *fn, int left, int right, int llength,
		int nwidth)
{
    /* int left,right;
       static int i = 0;
       char fn[1024]; 
       i++;
       sprintf(fn,"dump.%d.%d",getpid(),i);
    */

    FILE *fp;

    if (xx->editorState == StateDown)
	return;

    if (llength > MAX_LINE_LENGTH)
	llength = MAX_LINE_LENGTH;

    if ( (fp = fopen(fn,"w")) != NULL ) {
	/* extents(xx, &left, &right); 
	   bell();
	   */
	dumpRegion(xx,fp,left,right,llength, nwidth);
	fclose(fp);
    }
}


/*
 *----------------------------------------------------------------------------
 * The consensus_trace command
 *---------------------------------------------------------------------------
 */
int save_consensus_trace(EdStruct *xx, char *fn, int left, int right,
			 int strand, int matching) {
    Read *r;
    int ret;

    if (NULL == (r = cons_trace(xx, left, right, strand, matching, 0)))
	return -1;
    
    ret = write_reading(fn, r, TT_SCF);

    read_deallocate(r);

    return ret;
}

/*
 *----------------------------------------------------------------------------
 * The align command
 *---------------------------------------------------------------------------
 */

/*
 * Aligns a section of reading with the consensus.
 * Currently we align the same length of data from each, this means that the
 * endings of alignments are often incorrect. To align 'x' bases of reading
 * with 'x' + 'y' bases of consensus implies knowing about moving readings,
 * analysing pad insert positions, etc.
 */
int align_read(EdStruct *xx) {
    int seq, start, length, score, c_pos, c_length;
    align_int *res;
    char *con, *read;
    char *con_, *read_;
    int max_pads = 0; /* 20; */
    int tmp_l;

    if (xx->editorState == StateDown)
	return 1;

    /* check in valid edit mode */
    if (!(DBI_flags(xx) & DB_ACCESS_UPDATE)) {
	verror(ERR_WARN, "contig_editor", "Editor is in read-only mode");
	return 1;
    }
    
    /* find the selected region */
    if (!getSelection(xx, &seq, &start, &length, NULL) || seq == 0) {
	bell();
	return 1;
    }

    vfuncheader("Align reading (contig editor)");

    /* Determine the strings to align */
    start --;
    c_pos = start-DB_Start(xx, seq) + DB_RelPos(xx, seq) - max_pads;
    c_length = length + 2*max_pads;
    if (c_pos < 1) {
	c_length -= 1 - c_pos;
	c_pos = 1;
    }

    tmp_l = MAX(length, c_length);
    if (NULL == (con_ = (char *)xcalloc(tmp_l+1+10,1)) ||
	NULL == (read_ = (char *)xcalloc(tmp_l+1+10,1)) ||
	NULL == (res = (align_int *)xcalloc((tmp_l*2+1),sizeof(align_int)))) {
	return 1;
    }

    con = &con_[5];
    read = &read_[5];

    DBcalcConsensus(xx, c_pos, c_length, con, NULL,  BOTH_STRANDS);
    strncpy(read, &(DB_Seq(xx, seq)[start]), length);
    read[length] = '\0';
    con[c_length] = '\0';

    /* align the two */
    score = calign(read, con, length, c_length,
		   NULL, NULL, NULL, NULL,
		   0, 0, gopenval, gextendval, 0, 0, res);
    vmessage("alignment returned %d\n",score);
    cdisplay(read, con, length, c_length, 0, res, start, c_pos);
    vmessage("\n\n");

    openUndo(DBI(xx));

    /* Add the pads */
    {
	int i, j, k, l, ic = 0, ir = 0;
	int tmp_conf_n;
	align_int *S;
	/*
	 * Padding characters are inserted into contig in blocks of up to
	 * BLOCK at a time
	 */
#define BLOCK 20
	char pads[BLOCK+1] = "********************";
	
	tmp_conf_n = xx->default_conf_n;
	xx->default_conf_n = -1;

	S = res;
	start -= DB_Start(xx, seq);

	i = j = 0;
	while (i < c_length || j < length) {
	    if (*S == 0) {
		i++;j++;
		S++;
	    } else if (*S < 0) {
		/*
		 * Insert pad to consensus (all reads except the one we're
		 * aligning with)
		 */
		k = -*S;
		do {
		    int pos;

		    l = k > BLOCK ? BLOCK : k;
		    pos = c_pos + i + ic;

		    insertBasesConsensus(xx, pos, l, pads);

		    if ((pos <= DB_RelPos(xx, seq) + DB_Length(xx, seq)) &&
			(pos >= DB_RelPos(xx, seq))) {
			deleteBases(xx, seq, start + j + 1 + ir, l);
		    } else if (pos < DB_RelPos(xx, seq)) {
                        shiftLeft(xx, seq, l);
                    }

		    k -= l;
		    ic += l;
		} while (k);

		j -= *S;
		S++;
	    } else {
		/* insert pad to reading */
		k = *S;
		do {
		    l = k > BLOCK ? BLOCK : k;

		    insertBases(xx, seq, start + j + 1 + ir, l, pads);
                    if (start + j + 1 + ir <= 0)
                        shiftRight(xx, seq, l);
		    else
			ir += l;
		    k -= l;
		} while (k);

		i += *S;
		S++;
	    }
	}

	xx->default_conf_n = tmp_conf_n;
    }

    closeUndo(xx, DBI(xx));

    xx->refresh_flags |= ED_DISP_ALL;
    redisplaySequences(xx, 1);

    xfree(con_);
    xfree(read_);
    xfree(res);

    return 0;
}


/*
 *----------------------------------------------------------------------------
 * The strip pads command
 *---------------------------------------------------------------------------
 */


/*
 * Removes consensus pads where there is a COMPLETE column of pads.
 *
 * Returns the number of pads removed
 *         -1 for error.
 */
static int remove_pad_columns(EdStruct *xx, int sp_consensus_mode,
			      float sp_consensus_cutoff) {
#define CHUNKSIZE 8192
    char buffer[CHUNKSIZE+1];
    int i, j, width, count=0;
    int tmp;
    int npads = 0;

    xx->cursorPos = positionInContig(xx, xx->cursorSeq, xx->cursorPos);
    xx->cursorSeq = 0;

    for (i = DB_Length(xx, 0); i >= 1; i -= CHUNKSIZE) {
	width = min(i, CHUNKSIZE);
	tmp = consensus_mode;
	consensus_mode = sp_consensus_mode;
	calc_consensus(0, i-width+1, i, CON_SUM, buffer, NULL, NULL, NULL,
		       sp_consensus_cutoff, -1, contEd_info, (void *)xx);
	consensus_mode = tmp;

	for (j = width-1; j >= 0; j--) {
	    if (buffer[j] == '*') {
		npads++;
		deleteBasesConsensus(xx, i-width+1 + j, 1);

		if (i-width+1+j < xx->cursorPos)  count--;
		if (i-width+1+j < xx->displayPos) xx->displayPos--;
	    }
	}
    }

    U_adjust_cursor(xx, count);
    return npads;
}

/*
 * The original editor-based code, copied back here from a 20th April
 * 2001 version of edCommands.c. It's largely redundant due to the new
 * re-alignment code but is left here primarily for those that like to
 * see pad shuffling happen inside the editor. In time this needs to
 * be replaced by the same re-alignment code.
 */
int strip_pads(EdStruct *xx, int sp_consensus_mode,
	       float sp_consensus_cutoff)
{
    int npads;
    int old_undo = DBI_store_undo(xx);

    /* Don't attempt to store undo data for large contigs. */
    if (DB_Length2(xx, 0) > 1000000) {
	verror(ERR_WARN, "contig_editor",
	       "Disabling undo data as pad stripping produces too many edits");
        freeAllUndoLists(xx);
	DBI_store_undo(xx) = 0;
    }

    openUndo(DBI(xx));
    npads = remove_pad_columns(xx, sp_consensus_mode, sp_consensus_cutoff);
    if (npads > 0) {
	invalidate_consensus(xx);

	xx->refresh_flags |= ED_DISP_ALL;
	redisplaySequences(xx, 0);
    }
    closeUndo(xx,DBI(xx));

    if (npads <= 0)
	undoLastCommand(xx);

    DBI_store_undo(xx) = old_undo;

    return 0;
}
