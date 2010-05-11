/*
 * File: seqInfo.c
 * Version:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: misc routines in a temp file
 *
 * Created:
 * Updated:
 *
 */

#include <stdio.h>
#include <string.h> /* IMPORT: strchr, strdup */

#include <os.h>
#include <io_lib/expFileIO.h>
#include <io_lib/traceType.h>

#include "seqInfo.h"
#include "array.h"
#include <io_lib/scf_extras.h>

#include "misc.h"   /* IMPORT: strdup */
#include "FtoC.h"
#include "xalloc.h"
#include "IO.h"
#include "parse_ft.h"
#include "tagUtils.h"
#include "tagdb.h"
#include <io_lib/open_trace_file.h>

#ifdef USE_BIOLIMS
#include "spBiolims.h"
#endif

/*************************************************************
 * Miscellaneous
 *************************************************************/
void trim_white_space(char *s)
/*
 * trim trailing white space off the string s
 */
{
    int i;
    for(i = strlen(s)-1; i>=0 && isspace(s[i]); i--) s[i]='\0';
}


/*************************************************************
 * Create and free routines
 *************************************************************/
SeqInfo *allocSeqInfo(void)
/*
 *
 */
{
    SeqInfo *si;
    if ( (si = (SeqInfo *)xmalloc(sizeof(SeqInfo))) != NULL) {
	si->confidence = NULL;
	si->origpos = NULL;
	si->length = -1;
	si->start = -1;
	si->end = -1;
    }

    return si;
}



void freeSeqInfo(SeqInfo *si)
/*
 *
 */
{
#define FREE(F) \
    if (si->F != NULL) {\
        xfree(si->F);\
	si->F = NULL;\
    }

    if (si != NULL) {
	if(si->e != NULL) { exp_destroy_info(si->e); si->e = NULL; }
	FREE(confidence);
	FREE(origpos);
	xfree(si);
    }
}


/*************************************************************
 * Staden format files
 *************************************************************/

Exp_info *exp_read_staden_info(mFILE *fp, char *filename)
/*
 * Read a staden file into an Exp_info data structure
 */
{
    Exp_info *e;
    char *seq;
    char *fn;

    /* find last / in filename */
    for(fn=filename+strlen(filename)-1;fn>filename && *fn!='/';fn--);
    if (*fn=='/') fn++;
    
    
    if ( (e = exp_create_info()) != NULL ) {
	char line[128];
	/* an upper limit for the file size */
	int max_seq_size = file_size(filename);
	int left, right, len, i;
	int lineno;
	int formatIsBap;
	int CS_from, CS_to;
	int SR;

	CS_from = CS_to = SR = 0;

	/*ID*/
	(void)ArrayRef(e->entries[EFLT_ID],e->Nentries[EFLT_ID]++);
	exp_get_entry(e,EFLT_ID) = (char *)strdup(fn);

	/*EN*/
	(void)ArrayRef(e->entries[EFLT_EN],e->Nentries[EFLT_EN]++);
	exp_get_entry(e,EFLT_EN) = (char *)strdup(fn);

	/*CC*/
	(void)ArrayRef(e->entries[EFLT_CC],e->Nentries[EFLT_CC]++);
	exp_get_entry(e,EFLT_CC) = (char *)strdup("Created from a staden format sequence assembly file");

	seq = (char *) xmalloc(max_seq_size+1);
	if (!seq)
	    return NULL;

	left = 0;
	right = 0;
	len = 0;
	    
	lineno = 0;
	formatIsBap = 1;
	while (mfgets(line,sizeof(line),fp)!=NULL) {
	    char *s;
	    lineno++;
	    if (lineno == 1) {
		int pos, ret;
		char *cp;
		/*
		 * This maybe a fasta format file.
		 */
		if (line[0] == '>') {
		    if (cp = strchr(line, ' '))
			*cp = 0;
		    if (cp = strchr(line, '\t'))
			*cp = 0;
		    if (cp = strchr(line, '\n'))
			*cp = 0;
		    exp_set_entry(e, EFLT_ID, strdup(line+1));
		    exp_set_entry(e, EFLT_EN, strdup(line+1));
		    continue;
		}

		/*
		 * This maybe a file created from 'output consensus',
		 * in which case it'll have the Staden format style of:
		 * " <-----T2.00047----->" header on the first line.
		 *
		 * We strip this off. Argh - why don't Suns have the
		 * ANSI function memmove() ?
		 */
		ret = sscanf(line, " <%*18s>%n", &pos);
		if (ret && pos == 21) {
		    int i, e = sizeof(line)-21;

		    for (i=0; i < e; i++)
			line[i] = line[i+21];
		    /* memmove((void *)line,
		     *	(void *)(line + 21),
		     *	sizeof(line)-21); 
		     */
		}
	    }
	    if (line[0] == ';') {
		/********************************************
		 * Title line parsing
		 *******************************************/
		if (lineno==1 &&
		    !(line[1] == ';' || line[1] == '<' || line[1] == '>')
		    ) {
		    /* format of line is:
		     * <ln> ::= ;<i6><i6><i6><c4><c+>
		     * <i6> ::= <s><s><s><s><s><d> | ... <d><d><d><d><d><d>
		     * <c4> ::= <C><C><C><s>
		     * <c+> ::= <a><c+> | <a>
		     * <a> is [a-zA-Z0-9.]
		     * <s> is ' '
		     * <d> is [0-9]
		     * <C> is [A-Z]
		     */
		    int d;
		    if (sscanf(line,";%6d%6d%6d",&d,&d,&d)==3 &&
			strlen(line)>23) {
			/* trim off trailing white space */
			trim_white_space(line+23);
			/*LN*/
			(void)ArrayRef(e->entries[EFLT_LN],e->Nentries[EFLT_LN]++);
			exp_get_entry(e,EFLT_LN) = (char *)strdup(line+23); line[23]='\0';
			/*LT*/
			(void)ArrayRef(e->entries[EFLT_LT],e->Nentries[EFLT_LT]++);
			/* trim off trace type white space */
			trim_white_space(line+19);
			exp_get_entry(e,EFLT_LT) = (char *)strdup(line+19);
		    }
		} else if (formatIsBap) {
		    switch (line[1]) {
		    case '<':
			for(s=&line[2];*s;s++) {
			    if(!isspace(*s) && isprint(*s))
				seq[left++] = *s;
			}
			break;
		    case '>':
			for(s=&line[2];*s;s++) {
			    if(!isspace(*s) && isprint(*s))
				seq[max_seq_size-right++] = *s;
			}
			break;
		    case ';':
			/*TG*/
#if 0
			trim_white_space(line);
			(void)ArrayRef(e->entries[EFLT_TG],e->Nentries[EFLT_TG]++);
			/* convert format from Staden format to
			 * Experiment file format
			 */
			{
			    char *cp;
			    int pos, len;
			    char type[5];

			    cp = (char *)xmalloc(strlen(line)+20);
			    if (cp == NULL)
				break;

			    sscanf(line, ";;%4s %6d %6d",
				   type, &pos, &len);
			    /*
			     * Need to add 'left' to each start position
			     * in tag. ASSUMPTION: 'left' has already been
			     * defined. Ie that the ;< lines are before
			     * any ;; lines.
			     */
			    pos += left;
			    values2tag(cp, type, pos, pos + len - 1,
				       2, &line[20]);

			    exp_get_entry(e,EFLT_TG) = cp;
			}

			if (strncmp(line+2,"IGNC",4)==0 ||
			    strncmp(line+2,"CVEC",4)==0) {
			    CS_from = atoi(&line[2+4+1]);
			} else if (strncmp(line+2,"IGNS",4)== 0 ||
				   strncmp(line+2,"SVEC",4)== 0) {
			    SR = 1;
			}
#endif
			break;

		    default:
			break;
		    }
		}
	    } else {
		/********************************************
		 * The actual sequence bit
		 *******************************************/
		formatIsBap = 0; /* turn off title line parsing stuff */
		for (s=line;*s;s++) {
		    if(!isspace(*s) && isprint(*s))
			seq[left+len++] = *s;
		}
		    
	    }
	}

	/*
	 * The right cutoff has been stashed into the end of the array
	 * Move to correct place
	 */
	for(i=(max_seq_size-(left+len))/2;i>=0;i--) {
	    char temp;
	    /* swap */
	    temp = seq[left+len+i];
	    seq[left+len+i] = seq[max_seq_size-i];
	    seq[max_seq_size-i] = temp;
	}
	/* null terminate */
	seq[left+len+right] = '\0';

	/*SQ*/
	(void)ArrayRef(e->entries[EFLT_SQ],e->Nentries[EFLT_SQ]++);
	exp_get_entry(e,EFLT_SQ) = seq;

	/*SL*/
	sprintf(line,"%d",left);
	(void)ArrayRef(e->entries[EFLT_SL],e->Nentries[EFLT_SL]++);
	exp_get_entry(e,EFLT_SL) = (char *)strdup(line);

	/*SR*/
	if (SR) {
	    sprintf(line,"%d",left+len+1);
	    (void)ArrayRef(e->entries[EFLT_SR],e->Nentries[EFLT_SR]++);
	    exp_get_entry(e,EFLT_SR) = (char *)strdup(line);
	}

	/*CS*/
	if (CS_from) {
	    if (CS_from == 1) {
		CS_to = left;
	    } else {
		CS_from = left + len + 1;
		CS_to = left + len + right;
	    }
	    sprintf(line,"%d..%d",CS_from,CS_to);
	    (void)ArrayRef(e->entries[EFLT_CS],e->Nentries[EFLT_CS]++);
	    exp_get_entry(e,EFLT_CS) = (char *)strdup(line);
	}

	/*QR*/
	if (!SR && !CS_from) {
	    sprintf(line,"%d",left+len+1);
	    (void)ArrayRef(e->entries[EFLT_QR],e->Nentries[EFLT_QR]++);
	    exp_get_entry(e,EFLT_QR) = (char *)strdup(line);
	}

#if 0
	/*TG*/
	{
	    int i;
	    /* need to add LEFT to each start position in tag */
	    for(i=0;i<e->Nentries[EFLT_TG];i++) {
		sprintf(line,"%4.4s %6d%s",
			arr(char *,e->entries[EFLT_TG],i),
			atoi(arr(char *,e->entries[EFLT_TG],i)+5)+left,
			arr(char *,e->entries[EFLT_TG],i)+11);
		xfree(arr(char *,e->entries[EFLT_TG],i));
		arr(char *,e->entries[EFLT_TG],i) = (char *)strdup(line);

	    }
	}
#endif
    }
    return e;
}



/*************************************************************
 * Experiment file files... use expFileIO.c
 *************************************************************/





/*************************************************************
 * Utilities
 *************************************************************/

static void determine_active_region(Exp_info *e, int *left, int *right,
				    int ignore_vec)
{
    /*
     * need to determine active region
     */
    int SQlen;
    int CSfrom, CSto, SL, SR, QL, QR;
    int l,r;
    int noCS; /* if there is no CS line */

    SQlen = strlen(exp_get_entry(e,EFLT_SQ));
    noCS = exp_get_rng(e,EFLT_CS,&CSfrom,&CSto);
    if (exp_get_int(e,EFLT_SL,&SL)) SL = 0;
    if (exp_get_int(e,EFLT_SR,&SR)) SR = SQlen+1;
    if (exp_get_int(e,EFLT_QL,&QL)) QL = 0;
    if (exp_get_int(e,EFLT_QR,&QR)) QR = SQlen+1;

    if (ignore_vec) {
	*left = QL;
	*right = QR;
	return;
    }

    /* find where the good bits are */
    if (SL>QL) l = SL; else l = QL;
    if (SR<QR) r = SR; else r = QR;

    /*
     * NOTE: tmp fix to remove CS clipping.
     */
    noCS = -1;
    if (!noCS) {
	/*
	 * Six posibilities:
	 *           l                r
	 * 1  -----  >                <
	 * 2  -------|------->        <
	 * 3  -------|<---------------|-------
	 * 4         >    <--------   |
	 * 5         >    <-----------|-------
	 * 6         >                <  -----
	 *
	 * key: [-] CS [>] new left c.o [<] new right c.o [|] orig c.o
	 *
	 * Coded for understanding not efficiency!!!
	 */
	if (/*1*/ CSfrom <= l && CSto <= l)
	    l = l; /* null statement */
	else if (/*2*/ CSfrom <= l+1 && CSto < r)
	    l = CSto;
	else if (/*3*/ CSfrom <= l+1 && CSto >= r)
	    r = l+1;
	else if (/*4*/ CSfrom < r && CSto < r)
	    r = CSfrom;
	else if (/*5*/ CSfrom < r && CSto >= r)
	    r = CSfrom;
	else if (/*6*/ CSfrom >= r && CSto >= r)
	    l = l; /* null statement */
    }
    if (l>r) l = r-1; /* just a quick consistency check */

    *left = l;
    *right = r;
}


/* ------------------------------------------------------------------------- */

/*
 * Converts FT lines to TG lines.
 * This is a workaround for the fact that gap4 (currently) doesn't
 * have a feature table structure, so we generate (splitting if
 * required) tags instead.
 */
void parse_features(Exp_info *e) {
    int i;
    ft_entry *entry;
    ft_range *range;
    char *tag;
    int comment_len;
    int tag_len;
    char *comment = NULL;
    int t;
    char tag_type[5];
    int feat_num = 0;
    int ele_num;

    for (i = 0; i < e->Nentries[EFLT_FT]; i++) {
	entry = parse_ft_entry(arr(char *, e->entries[EFLT_FT], i));

	if (!entry)
	    continue;

	/*
	 * Our tag comment will be type, location and full qualifiers
	 */
	comment_len  =
	    15 + 1 + 12 + 1 + /* "#FEATURE 000000 ELEMENT 000" */
	    strlen(entry->type) + 1 + 
	    strlen(entry->location ) + 1 +
	    (entry->qualifiers ? strlen(entry->qualifiers) : 0) + 1
	    + 5 /* 1 for null and 4 for paranoia */;

	/*
	 * The TG string comprises of (lengths in brackets):
	 * tag type (4 + 1)
	 * strand (1 + 1)
	 * location (max 10+2+10 + 1)
	 * newline (1)
	 * comment (comment_len from above)
	 * + 1 for null, and another 10 for paranoia! (newlines hopefully
	 * accounted for, but...).
	 */
	tag_len = 4+1 + 1+1 + 10+2+10+1 + 1 + comment_len + 11;
	comment = (char *)xmalloc(comment_len);
	if (!comment)
	    return;
	sprintf(comment, "#FEATURE 000000 ELEMENT 000\n%s\n%s\n%s",
		entry->type,
		entry->location,
		entry->qualifiers ? entry->qualifiers : "");
	feat_num++;
	ele_num = 0;
	for (range = entry->range; range; range = range->next) {
	    int start = INT_MAX, end = -INT_MAX;

	    if (!range->left) {
		verror(ERR_WARN, "parse_features", "invalid range");
		continue;
	    }

	    if (range->left) {
		start = range->left->min;
		end = range->left->max;
	    }

	    if (range->right) {
		start = MIN(start, range->right->min);
		end = MAX(end, range->right->max);
	    }

	    tag = (char *)xmalloc(tag_len);
	    if (!tag)
		continue; /* verror has already been called */

	    /* Find the 4-letter tag code from the feature key name */
	    strcpy(tag_type, "F---"); /* default */
	    for (t = 0; t < tag_db_count; t++) {
		char cmp[1024];
		if (!tag_db[t].type)
		    continue;

		sprintf(cmp, "FEATURE: %s", entry->type);
		if (strcmp(tag_db[t].type, cmp) == 0) {
		    memcpy(tag_type, tag_db[t].id, 4);
		    break;
		}
	    }

	    /* Set feature number at start of comment */
	    sprintf(comment+9, "%06d", feat_num);
	    comment[15] = ' ';

	    /* Set element number */
	    sprintf(comment+24, "%03d", ele_num);
	    comment[27] = '\n';
	    ele_num++;

	    /*
	     * Each entry may contain multiple locations. We duplicate the tag
	     * for each location, keeping the comment identical.
	     */
	    if (-1 == values2tag(tag, tag_type, start, end,
				 range->complemented, comment)) {
		verror(ERR_WARN, "parse_features",
		       "couldn't create tag from feature table entry");
		continue;
	    }

	    exp_set_entry(e, EFLT_TG, tag);
	    xfree(tag);	
	}

	xfree(comment);
    }
}
/* ------------------------------------------------------------------------- */


/*************************************************************
 *
 *************************************************************/



/*
 * Reads in sequence from file filname.
 * Format of this file must be either:
 * 	staden (80 char lines, title lines prefixed with ";")
 * or   experiment file
 *
 * "ignore_vec" specifies whether to ignore the vector sequence lines (SL,
 * SR, CL, CL, CS) when determining the active region. This is necessary
 * when inputting extracted data to preassembly.
 *
 * returns length read in
 */
SeqInfo *read_sequence_details(char *filename, int ignore_vec)
{
    SeqInfo *si = NULL;
    Exp_info *e;
    mFILE *fp;
    int format;

    /*
     * read sequence details into experiment file format
     */
    if (NULL == (fp = open_exp_mfile(filename, NULL)))
	return NULL;

    format = fdetermine_trace_type(fp);
    mrewind(fp);
    switch(format) {
    case TT_PLN:
	e = exp_read_staden_info(fp, filename);
	mfclose(fp);
	break;
    case TT_EXP:
	e = exp_mfread_info(fp);
	mfclose(fp);
	if (e)
	    exp_close(e);
	break;
#ifdef USE_BIOLIMS
    case TT_BIO:
        e = spBiolims2exp(filename);
	mfclose(fp);
	break;
#endif
    case TT_ERR:
	verror(ERR_WARN, "read_sequence_details",
	       "Failed to read file %s", filename);
	e = NULL;
	mfclose(fp);
	break;
    default:
	verror(ERR_WARN, "read_sequence_details",
	       "File %s is not in plain or Experiment File format", filename);
	e = NULL;
	mfclose(fp);
	break;
    }

    if ( e != NULL ) {
	if (e->Nentries[EFLT_SQ] == 0 || ( si = allocSeqInfo() ) == NULL ) {
	    exp_destroy_info(e);
	    return NULL;
	} else {
	    si->e = e;
	    si->length = strlen(exp_get_entry(e,EFLT_SQ));
	    determine_active_region(e, &si->start, &si->end, ignore_vec);
	}

	/* orig pos and conf. values */
	if (e->Nentries[EFLT_ON]) {
	    int2 *opos;

	    if (NULL != (opos = (int2 *)xmalloc((si->length+1)*
						sizeof(int2)))) {
		if (str2opos(opos, si->length+1, exp_get_entry(e, EFLT_ON))
		    != si->length)
		    verror(ERR_WARN, "read_sequence_details",
			   "Experiment file %s - 'ON' line has wrong number "
			   "of items", filename);
		si->origpos = opos;
	    } else
		si->origpos = NULL;
	}

	if (e->Nentries[EFLT_AV]) {
	    int1 *conf;

	    if (NULL != (conf = (int1 *)xmalloc((si->length+1)*
						sizeof(int1)))) {
		if (str2conf(conf, si->length+1, exp_get_entry(e, EFLT_AV))
		    != si->length)
		    verror(ERR_WARN, "read_sequence_details",
			   "Experiment file %s - 'AV' line has wrong number "
			   "of items", filename);
		si->confidence = conf;
	    } else
		si->confidence = NULL;
	}

	if (e->Nentries[EFLT_FT]) {
	    parse_features(e);
	}
    }

    return si;
}


/*
 * Reads an ID line from an experiment file, which has already been loaded
 * into a SeqInfo structure.
 *
 * Returns name for success
 *         NULL for failure
 *
 * The name returned is valid only until the next call to idline().
 */
char *read_sequence_name(SeqInfo *si) {
    static char name[DB_NAMELEN+1];
    char *namep;
    int i;

    if (exp_Nentries(si->e, EFLT_ID) < 1) {
	verror(ERR_WARN, "read_sequence_name", "No ID line in experiment file");
	if (exp_Nentries(si->e, EFLT_EN) < 1) {
	    verror(ERR_WARN, "read_sequence_name", "Not even an EN line!");
	    return NULL;
	} else {
	    namep = exp_get_entry(si->e, EFLT_EN);
	}
    } else {
	namep = exp_get_entry(si->e, EFLT_ID);
    }

    /* Copy up to the first white space */
    i = 0;
    do {
	name[i++] = *namep++;
    } while (i < DB_NAMELEN &&
	     *namep != ' ' && *namep != '\t' &&
	     *namep != '\n' && *namep != '\r' &&
	     *namep != '\0');
    name[i] = 0;

    return name;
}


/*
 * Fortran interface to read_sequence_name. Only now used in dbas.f
 */
f_proc_ret idline_(char *NAMARC,
		  char *eID,
		  f_implicit NAMARC_l,
		  f_implicit eID_l) {
    char filename[1024], *name;
    SeqInfo *si;
    
    Fstr2Cstr(NAMARC, NAMARC_l, filename, 1023);

    /*
     * read in whole of original sequence
     */
    if ( (si = read_sequence_details(filename, 0)) == NULL ) f_proc_return();

    if (NULL == (name = read_sequence_name(si))) {
	freeSeqInfo(si);
	f_proc_return();
    }

    Cstr2Fstr(name, eID, eID_l);
    freeSeqInfo(si);

    f_proc_return();
}

/*
 * Fills a buffer 'conf' of length 'length' with the confidence values.
 *
 * These are either read from si, read from the trace file, or guessed.
 */
void SeqInfo_conf(SeqInfo *si, int1 *conf, int length) {
    if (si->confidence) {
	memcpy(conf, si->confidence, sizeof(int1) * length);
    } else {
	/* Read from SCF file */
	if (0 != get_read_conf(si->e, length, NULL, conf)) {
	    int i;

	    for (i = 0; i < length; i++)
		conf[i] = 2;
	}
    }
}


/*
 * Fills a buffer 'opos' of length 'length' with the original positions.
 *
 * These are either read from si, or guessed.
 */
void SeqInfo_opos(SeqInfo *si, int2 *opos, int length) {
    if (si->origpos) {
	memcpy(opos, si->origpos, sizeof(int2) * length);
    } else {
	int i, j;
	char *seq;

	seq = exp_get_entry(si->e, EFLT_SQ);

	for (i = j = 0; i < length; i++)
	    if (seq[i] != '*')
		opos[i] = ++j;
	    else
		opos[i] = 0;
    }
}


/*
 * Returns a pointer to the used component of the sequence, along with the
 * length (if a pointer is supplied).
 */
char *SeqInfo_used_seq(SeqInfo *si, int *length) {
    char *seq;

    if (length)
	*length = si->end - si->start - 1;
    
    seq = exp_get_entry(si->e, EFLT_SQ);

    return seq ? &seq[si->start] : NULL;
}
