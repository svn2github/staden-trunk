/*
 * File: IO2.c
 * Version:
 *
 * Author: Simon Dear
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

#include <stdio.h>
#include <ctype.h>

#include "IO.h"
#include "IO2.h"

#include "tagUtils.h" /* IMPORT: create_tag_for_gel */
#include "misc.h"
#include "gap-create.h" /* IMPORT: cpdb */
#include "io-reg.h"
#include "Read.h"
#include "scf_extras.h"
#include "FtoC.h"
#include "contig_selector.h"
#include "clones.h"
#include "notes.h"

/*************************************************************
 * Sequence manipulation routines
 *************************************************************/

int io_complement_seq(int  *length,	/* length of complete string */
		      int  *start,	/* start */
		      int  *end,	/* end */
		      char *seq,	/* complete sequence */
		      int1 *conf,	/* confidence vals */
		      int2 *opos)	/* original pos */
/*
 * We can ignore updating the sense because the calling routine does this
 * explicitly.
 */
{
    f_int len;

    int1 tint1;			/* temporary int1 */
    int2 tint2;			/* temporary int2 */
    int i;
    int orig_start;
    int orig_end;

    /* complement sequence */
    len = *length;
    (void) sqcom_(seq,&len,(f_implicit)1);
    (void) sqrev_(seq,&len,(f_implicit)1);

    /* reverse start and end pointers */
    orig_start = *start;
    orig_end = *end;
    *start = *length - orig_end + 1;
    *end = *length - orig_start + 1;

    if (conf && opos) {
	/* reverse conf and opos*/
	for (i=len/2;i>0;i--) {
	    tint1 = conf[i-1], conf[i-1] = conf[len-i], conf[len-i] = tint1;
	    tint2 = opos[i-1], opos[i-1] = opos[len-i], opos[len-i] = tint2;
	}
    }

    return 0;
}



/*
 * default confidence value stored here
 */
int DEFAULT_CONFIDENCE = 100;

static void assign_conf(char *seq, int1 *conf, int pos, int length) {
    int j, left, right;

    left=right=-1;
    
    /* Find previous and next non pad characters */
    for (j=pos-2; j>=0; j--) {
	if (seq[j] != '*') {
	    left=conf[j];
	    break;
	}
    }
    
    for (j=pos; j<length; j++) {
	if (seq[j] != '*') {
	    right=conf[j];
	    break;
	}
    }
    
    /* Derive conf as average of left and right values */
    if (left == -1)
	if (right == -1)
	    conf[pos-1] = 0;
	else
	    conf[pos-1] = right;
    else
	if (right == -1)
	    conf[pos-1] = left;
	else
	    conf[pos-1] = (left + right)/2;
}

/*
 * Insert bases into a sequence
 *
 * Assumes that the caller knows that they are doing and that the seq, conf
 * and opos arrays have been allocated large enough to accomodate Nbases extra
 * bases.
 */
int io_insert_seq(int  *length,	/* length of complete string */
		  int  *start,	/* start */
		  int  *end,	/* end */
		  char *seq,	/* complete sequence */
		  int1 *conf,	/* confidence vals */
		  int2 *opos,	/* original pos */
		  /**********/
		  int  pos,
		  char *bases,
		  int1 *newconf,
		  int2 *newopos,
		  int  Nbases
		  )
{
    int i;

    /* shift existing bases */
    for (i=*length-1 + Nbases; i>=pos+Nbases-1; i--) {
	seq[i] = seq[i-Nbases];
	conf[i] = conf[i-Nbases];
	opos[i] = opos[i-Nbases];
    }

    /* insert new bases */
    for(i=0; i<Nbases; i++) {
	seq[i+pos-1] = bases[i];
	if (newconf==NULL)
	    conf[i+pos-1] = DEFAULT_CONFIDENCE;
	else {
	    if (newconf[i] != 255)
		conf[i+pos-1] = newconf[i];
	    else
		assign_conf(seq, conf, i+pos, *length);
	}
	if (newopos==NULL)
	    opos[i+pos-1] = 0;
	else
	    opos[i+pos-1] = newopos[i];
	    
    }

    /* adjust lengths and things */
    *length += Nbases;
    if (pos <= *start) *start += Nbases;
    if (pos <= *end) *end += Nbases;
    if (*start > *length+1) *start = *length + 1;
    if (*end > *length+1) *end = *length + 1;

    if (newconf == NULL) {
	/* Update confidence values if these bases are pads or conf == -1 */
	for (i=0; i < Nbases; i++) {
	    if (bases[i] == '*' || (signed char)conf[i] == -1)
		assign_conf(seq, conf, i+pos, *length);
	}
    }

    return 0;
}

int io_delete_seq(int  *length,	/* length of complete string */
		  int  *start,	/* start */
		  int *end,	/* end */
		  char *seq,	/* complete sequence */
		  int1 *conf,	/* confidence vals */
		  int2 *opos,	/* original pos */
		  /**********/
		  int  pos,
		  int  Nbases
		  )
/*
 * Delete bases from a sequence
 */
{
    int i;

    /* shift existing bases */
    for (i=pos+Nbases;i<=*length;i++) {
	seq[i-Nbases-1] = seq[i-1];
	conf[i-Nbases-1] = conf[i-1];
	opos[i-Nbases-1] = opos[i-1];
    }

    /*
     * Six possibilities:
     *
     *  pppppp S              E			S-=Nbases,E-=Nbases
     *  pppppppSpppppppppp    E			S =pos   ,E-=Nbases
     *  pppppppSppppppppppppppEppppppp		S =pos   ,E-=Nbases
     *         S   ppppppp    E			S =S     ,E-=Nbases
     *         S   pppppppppppEppppppp          S =S     ,E =pos
     *         S              E pppppp          S =S     ,E =E
     *
     * (S)tart (E)nd (p)osition/length
     */

    /* adjust length */
    *length -= Nbases;

    /* adjust start */
    if (pos <= *start) {
	if (pos + Nbases <= *start+1)
	    *start -= Nbases;
	else
	    *start = pos;
    }

    /* adjust end */
    if (pos < *end) {
	if (pos + Nbases < *end)
	    *end -= Nbases;
	else
	    *end = pos;
    }


    return 0;
}



int io_replace_seq(int  *length, /* length of complete string */
		   int  *start,	 /* start */
		   int  *end,	 /* end */
		   char *seq,	 /* complete sequence */
		   int1 *conf,	 /* confidence vals */
		   int2 *opos,	 /* original pos */
		   /**********/
		   int  pos,
		   char *bases,
		   int1 *newconf,
		   int2 *newopos,
		   int  Nbases,
		   int  diff_only, /* Only change those that are different */
		   int  conf_only  /* Only modify confidence  */
		   )
/*
 * Insert bases into a sequence
 */
{
    int i, j, left, right;

    /* replace new bases */
    for(i=0;i<Nbases;i++) {
	if (diff_only) {
	    if (toupper(seq[i+pos-1]) == toupper(bases[i]))
		continue;
	}

	if (newconf==NULL)
	    conf[i+pos-1] = DEFAULT_CONFIDENCE;
	else
	    conf[i+pos-1] = newconf[i];

	if (conf_only)
	    continue;

	seq[i+pos-1] = bases[i];

	if (newopos==NULL)
	    opos[i+pos-1] = 0;
	else
	    opos[i+pos-1] = newopos[i];
    }

    if (newconf == NULL) {
	/* Update confidence values if these bases are pads */
	for (i=0; i < Nbases; i++) {
	    if (bases[i] == '*') {
		left=right=-1;
		
		/* Find previous and next non pad characters */
		for (j=i+pos-2; j>=0; j--) {
		    if (seq[j] != '*') {
			left=conf[j];
			break;
		    }
		}
		
		for (j=i+pos; j<*length; j++) {
		    if (seq[j] != '*') {
			right=conf[j];
			break;
		    }
		}
		
		/* Derive conf as average of left and right values */
		if (left == -1)
		    if (right == -1)
			conf[i+pos-1] = 0;
		    else
			conf[i+pos-1] = right;
		else
		    if (right == -1)
			conf[i+pos-1] = left;
		    else
			conf[i+pos-1] = (left + right)/2;
	    }
	}
    }

    return 0;
}

/*
 * Pads in the consensus - rewrite in C. Currently uses PADCON().
 */
int pad_consensus(GapIO *io, int contig, int pos, int npads) {
    char *tmp_gel;
    int max_len;
    f_int cnum = io_dbsize(io) - contig;

    max_len = find_max_gel_len(io, contig, 0);
    max_len += npads + 1;
    if (NULL == (tmp_gel = (char *)xmalloc(max_len)))
	return -1;

    padcon_(&io_relpos(io, 1), &io_length(io, 1), &io_lnbr(io, 1),
	    &io_rnbr(io, 1), &NumReadings(io), &NumContigs(io),
	    tmp_gel, &cnum, &pos, &npads,
	    &io_dbsize(io), handle_io(io), &max_len, max_len);

    xfree(tmp_gel);

    return 0;
}


/*************************************************************
 * Code for writing a new gel to a database
 ************************************************************/

#include "seqInfo.h"
#include "array.h"



#define MAX_ARCHIVE_LENGTH 16
f_proc_ret stikit_(f_int *HANDLE,
		   char *NAMARC,
		   f_int *NGEL,
		   f_int *LENGTH,
		   char *GEL,
		   f_int *MAXGEL/*unused*/,
		   f_int *IOK,
		   f_int *CONTIG,
		   f_int *RELPG,
		   /*fortran implicits*/
		   f_implicit NAMARC_l,
		   f_implicit GEL_l)
/*
 * Write gel reading to database.
 * We assume that writeg_ and writec_ will have been called (if needed) before
 * here. (This is required for adding the contig tags.)
 */
{
    int err;
    int i,j;
    char namarc[F_NAMLEN+1], *name;
    SeqInfo *si;
    GapIO *io;
    int gel_len;

    /*
     * The first of the working versions
     */
    int length;	
    int alength; /* allocated length */
    int start;
    int end;
    char *seq = NULL;
    int1 *conf = NULL;
    int2 *opos = NULL;

    *IOK = 1;

    if ( (io = io_handle(HANDLE)) == NULL) goto end;

    /*
     * convert FORTRAN file name into a clean C string.
     */
    for (i=0;i<NAMARC_l && i<F_NAMLEN && !isspace(NAMARC[i]);i++)
	namarc[i]=NAMARC[i];
    namarc[i]='\0';



    /*
     * read in whole of original sequence
     */
    if ( (si = read_sequence_details(namarc, 0)) == NULL ) goto end;

    /* Get the reading name */
    if (NULL == (name = read_sequence_name(si))) {
	goto end;
    }

    /* initialise records for this gel */
    (void) io_init_reading(io, (int) *NGEL);

    /*
     * blank tag record
     */
    /*blank_tag_rec(io, (tag_id) *NGEL);*/

    /* Write the reading name */
    writen_(HANDLE, NGEL, name, strlen(name)+1);
  

    /*
     * Create default tags
     * Firstly for gel readings.
     */
    for(i=0;i<exp_Nentries(si->e,EFLT_TG);i++) {
	create_tag_for_gel(io, (int)*NGEL, si->length, 
			   arr(char *, si->e->entries[EFLT_TG], i),
			   NULL, 0, NULL);
    }

    /*
     * And secondly consensus tags
     */
    for (i=0; i < exp_Nentries(si->e, EFLT_TC); i++) {
	char *comment;
	char type[5];
	int start, end, strand;
	char *tag;

	tag = arr(char *, si->e->entries[EFLT_TC], i);
	if (NULL == (comment = (char *)xmalloc(strlen(tag))))
	    goto end;

	if (-1 == tag2values(tag, type, &start, &end, &strand, comment))
	    continue;

	if (*LENGTH >= 0) {
	    start += *RELPG - si->start - 1;
	    end   += *RELPG - si->start - 1;
	} else {
	    int len = end - start;

	    start = *RELPG + si->end - end - 1;
	    end   = start + len;
	}

	type[4] = '\0';

	insert_NEW_tag(io, (tag_id)-*CONTIG, start, end-start+1, type, comment,
		       strand);

	xfree(comment);
    }


    /*
     * And lastly, the SVEC/CVEC tags derived from the SL, SR, and CS lines
     */
    if (exp_Nentries(si->e, EFLT_SL)) {
	int start = 1, len;
	
	len = atoi(exp_get_entry(si->e, EFLT_SL));
	if (len > 0) {
	    insert_NEW_tag(io, (tag_id)*NGEL, start, len, "SVEC", NULL, 0);
	}
    }

    if (exp_Nentries(si->e, EFLT_SR)) {
	int start, end;
	
	/*
	 * Vepe is currently bugged and always creates an SR when we have an
	 * SL. Often this SR is simply the far right end of the reading, in
	 * which case we don't really have any right-end sequencing vector
	 */
	start = atoi(exp_get_entry(si->e, EFLT_SR));
	if (start < si->length) {
	    end = si->length;
	    
	    insert_NEW_tag(io, (tag_id)*NGEL, start, end - start + 1, "SVEC",
			   NULL, 0);
	}
    }

    if (exp_Nentries(si->e, EFLT_CS)) {
	int start, end;
	
	exp_get_rng(si->e, EFLT_CS, &start, &end);
	if (start < 1) start = 1;

	insert_NEW_tag(io, (tag_id)*NGEL, start, end - start + 1, "CVEC",
		       NULL, 0);
    }

    /*
     * Create NoTes.
     */
    for(i=0;i<exp_Nentries(si->e,EFLT_NT);i++) {
	create_note_for_gel(io, (int)*NGEL,
			    arr(char *, si->e->entries[EFLT_NT], i));
    }

    /*
     * Transfer to working versions
     */
    length = si->length;
    start = si->start;
    end = si->end;
    alength = length + 100;

    /* Allocate temporary buffers */
    seq = (char *)xmalloc(alength * sizeof(*seq));
    conf = (int1 *)xmalloc(alength * sizeof(*conf));
    opos = (int2 *)xmalloc(alength * sizeof(*opos));
    if (!seq || !conf || !opos)
	goto end;

    SeqInfo_conf(si, conf, length);
    SeqInfo_opos(si, opos, length);

    memcpy(seq, exp_get_entry(si->e,EFLT_SQ), length);

    /*
     * Complement if necessary
     */
    if (*LENGTH<0) io_complement_seq(&length,&start,&end,seq,conf,opos);



    /*
     * Set gel_len
     */
    if ( (gel_len = *LENGTH) < 0) gel_len = -gel_len;



    /*
     * Consistency check
     */
    for(i=0,j=start;i<gel_len;i++) {
	if (GEL[i]!=seq[j]) {
	    if (GEL[i]!='*') {
		freeSeqInfo(si);
		crash("Serious problem encountered entering sequence into database: %s\n",namarc);
	    }
	} else
	    j++;
    }


    /*
     * Adjust tags
     * NOTE: Must always traverse reading in reverse of original sense
     */
    if ( *LENGTH < 0 ) {
	for(i=0,j=start;i<gel_len;i++) {
	    if (GEL[i]!=seq[j]) {
		tag_shift_for_insert(io, (int)*NGEL,length-j+1);
	    } else
		j++;
	    
	}
    } else {
	for(i=gel_len-1,j=end-2;i>=0;i--) {
	    if (GEL[i]!=seq[j]) {
		tag_shift_for_insert(io, (int)*NGEL,j+2);
	    } else
		j--;
	}
    }



    /*
     * pad original sequence by examining pads in GEL.
     * ASSUMPTION: pads are only in the used data (valid assumption at time
     * of writing).
     */
    for(i=0;i<gel_len;i++) {
	if (GEL[i]!=seq[start+i]) {
	    if (length >= alength-1) {
		alength += 100; /* realloc sequence buffers */
		seq = xrealloc(seq, alength * sizeof(*seq));
		conf = xrealloc(conf, alength * sizeof(*conf));
		opos = xrealloc(opos, alength * sizeof(*opos));
		if (!seq || !conf || !opos)
		    goto end;
	    }
	    io_insert_seq(&length, &start, &end, seq,
			  conf, opos, start+i+1, "*", NULL, NULL, 1);
	}
    }


    
    /*
     * io_write_seq now updates io_length(), which breaks things as
     * r.sense hasn't been set yet (this is done by writeg_().
     * So we now have to manually do this first.
     */
    {
	GReadings r;
	gel_read(io, *NGEL, r);
	r.sense = *LENGTH < 0 ? 1 : 0;
	gel_write(io, *NGEL, r);
    }

    /*
     * write sequence to file
     */
    if (io_write_seq(io,(int)*NGEL,&length,&start,&end,seq,conf,opos)) {
	verror(ERR_WARN, "stikit",
	       "Problem writing new sequence to database: %s", namarc);
	goto end;
    }


    /*
     * write raw data info
     */
    if (exp_Nentries(si->e,EFLT_LT) &&
	exp_Nentries(si->e,EFLT_LN)) {
	if (io_write_rd(io,
			(int)*NGEL,
			exp_get_entry(si->e,EFLT_LN),
			strlen(exp_get_entry(si->e,EFLT_LN)),
			exp_get_entry(si->e,EFLT_LT),
			strlen(exp_get_entry(si->e,EFLT_LT)))) {
	    verror(ERR_WARN, "stikit",
		   "Problem writing raw data information to database: %s",
		   namarc);
	    goto end;
	}
    }

    /*
     * write everything else
     */
    err = add_seq_details(io, (int) *NGEL, si);

    /* exit */
    *IOK = 0;

 end:
    if (si) freeSeqInfo(si);
    if (seq) xfree(seq);
    if (conf) xfree(conf);
    if (opos) xfree(opos);
    
    f_proc_return();
}


/*************************************************************
 * Utility routines
 **************************************************************/

/* ARGSUSED */
/*
 * Complement a sequence
 */
f_proc_ret cplseq_(f_int *HANDLE,
		   f_int *N,
		   f_int *MAXGEL/*unused*/)
{
    GapIO *io;

    /*
     * The first of the working versions
     */
    int length;	
    int start;
    int end;
    char *seq = NULL;
    int1 *conf = NULL;
    int2 *opos = NULL;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    if (io_aread_seq(io, (int)*N, &length, &start, &end, &seq, &conf, &opos, 0)
	== 0) {
	io_complement_seq(&length, &start, &end, seq, conf, opos);
	(void)io_write_seq(io, (int)*N, &length, &start, &end, seq, conf, opos);
    }

    if (seq)	xfree(seq);
    if (conf)	xfree(conf);
    if (opos)	xfree(opos);

    f_proc_return();
}



/* ARGSUSED */
f_proc_ret arrfio_(char *NAMARC,
		   char *SEQ,
		   f_int *MAXSEQ,
		   f_int *MODE,	/* =1 clipped sequence =0 whole sequence */
		   f_int *IOK,
		   /* fortran implicits */
		   f_implicit NAMARC_l,
		   f_implicit SEQ_l)
/*
 * Read in a sequence
 *
 * Returns IOK as the error code.
 * 0 == success.
 * 1 == file not found
 * 2 == invalid file
 */
{
    int len;
    SeqInfo *si;
    char filename[F_NAMLEN+1];
    int i;

    for (i=0;i<NAMARC_l && i<F_NAMLEN && !isspace(NAMARC[i]);i++)
	filename[i]=NAMARC[i];
    filename[i]='\0';

    *IOK = 0;

    if ( (si = read_sequence_details(filename, 0)) == NULL) {
	*IOK = 1;
	len = 0;
    } else {

	if (exp_Nentries(si->e, EFLT_ID) < 1 &&
	    exp_Nentries(si->e, EFLT_EN) < 1) {
	    verror(ERR_WARN, "arrfio", "Invalid file format (No ID line)");
	    *IOK = 1;
	    freeSeqInfo(si);
	    f_proc_return();
	}

	if ((int)*MODE==0) {
	    /* include cutoff data */
	    len = si->length - si->start;
	} else {
	    /* ignore cutoff data */
	    len = si->end - si->start - 1;
	}
	if ( len > (int)*MAXSEQ) {
	    verror(ERR_WARN, "arrfio",
		   "Too much data. Maximum possible = %d, input stopped there",
		   (int)*MAXSEQ);
	    len = (int)*MAXSEQ;
	}
	if (len >= 0) {
	    /* copy over string */
	    strncpy(SEQ,exp_get_entry(si->e,EFLT_SQ)+si->start,len);
	}

	freeSeqInfo(si);
    }

    *MAXSEQ = (f_int)len;

    /* return len; */
    f_proc_return();
}


/************************
 * Base editing routines
 ************************/

/*
 * Inserts a single base into a gel reading ensuring that tags are shifted
 * and/or extended appropriately.
 *
 * Returns 0 for success, -1 for error.
 */
int io_insert_base(GapIO *io, int gel, int pos, char base) {
    /*
     * The first of the working versions
     */
    int length;	
    int start;
    int end;
    int npos;
    char *seq = NULL;
    int1 *conf = NULL;
    int2 *opos = NULL;
    int retcode = -1;

    if (io_aread_seq(io, gel, &length, &start, &end, &seq, &conf, &opos, 1)
	== 0) {
	io_insert_seq(&length,&start,&end, seq, conf, opos,
		      start + pos, &base, NULL, NULL, 1);
	(void)io_write_seq(io, gel, &length, &start, &end, seq, conf, opos);

	/* 
	 * Adjust tags
	 */
	if (io_length(io, gel) < 0)
	    npos = length - (pos + start) + 1;
	else
	    npos = pos + start;

	tag_shift_for_insert(io, gel, npos);

	retcode = 0;
    }

    if (seq)	xfree(seq);
    if (conf)	xfree(conf);
    if (opos)	xfree(opos);

    return retcode;
}

f_proc_ret insbas_(f_int *HANDLE,
		   f_int *N,
		   f_int *POS,
		   char *BASE,
		   /* fortran implicits */
		   f_implicit BASE_l)
{
    GapIO *io;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();
    (void)io_insert_base(io, *N, *POS, *BASE);

    f_proc_return();
}


/*
 * Replace a single base into a gel reading in the database.
 *
 * Returns 0 for success, -1 for error.
 */
int io_modify_base(GapIO *io, int gel, int pos, char base) {
    /*
     * The first of the working versions
     */
    int length;	
    int start;
    int end;
    char *seq = NULL;
    int1 *conf = NULL;
    int2 *opos = NULL;

    if (io_aread_seq(io, gel, &length, &start, &end, &seq, &conf, &opos, 0)
	== 0) {
	io_replace_seq(&length,&start,&end,seq,conf,opos,
		       start + pos, &base, NULL, NULL, 1, 0, 0);
	(void)io_write_seq(io, gel, &length, &start, &end, seq, conf, opos);
    }

    if (seq)	xfree(seq);
    if (conf)	xfree(conf);
    if (opos)	xfree(opos);

    return 0;
}


f_proc_ret modbas_(f_int *HANDLE,
		   f_int *N,
		   f_int *POS,
		   char *BASE,
		   /* fortran implicits */
		   f_implicit BASE_l)
{
    GapIO *io;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();
    io_modify_base(io, *N, *POS, *BASE);

    f_proc_return();
}




/*
 * Deletes a single base from a gel reading ensuring that tags are shifted
 * and/or extended appropriately.
 *
 * Returns 0 for success, -1 for error.
 */
int io_delete_base(GapIO *io, int gel, int pos) {
    /*
     * The first of the working versions
     */
    int length;	
    int start;
    int end;
    char *seq = NULL;
    int1 *conf = NULL;
    int2 *opos = NULL;
    int npos;
    int retcode = -1;

    if (io_aread_seq(io, gel, &length, &start, &end, &seq, &conf, &opos, 0)
	== 0) {
	io_delete_seq(&length,&start,&end,seq,conf,opos,
		      start + pos, 1);
	(void)io_write_seq(io, gel, &length, &start, &end, seq, conf, opos);

	/* 
	 * Adjust tags
	 */
	if (io_length(io, gel) < 0)
	    npos = length - (pos + start) + 1;
	else
	    npos = pos + start;

	tag_shift_for_delete(io, gel, npos);

	retcode = 0;
    }

    if (seq)	xfree(seq);
    if (conf)	xfree(conf);
    if (opos)	xfree(opos);

    return retcode;
}


f_proc_ret delbas_(f_int *HANDLE,
		   f_int *N,
		   f_int *POS)
{
    GapIO *io;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();
    (void)io_delete_base(io, *N, *POS);

    f_proc_return();
}



/*************************************************************
 * Get Extension
 *************************************************************/

/*
 * Gets an extension and complements it if needed.
 * Cutlen holds the length of 'cutoff' both before (length allowed to
 * return) and after (length actually returned) the call.
 */
int cgetext(GapIO *io, int gel, char *cutoff, int *cutlen) {
    int len, comp = 0;

    if (-1 == io_get_extension(io, gel, cutoff, *cutlen, &len, &comp))
	return -1;

    *cutlen = len;

    if (comp) {
	f_int flen = len;
	sqcom_(cutoff,&flen,(f_implicit)1);
	sqrev_(cutoff,&flen,(f_implicit)1);
    }
	
    return 0;
}


f_proc_ret getext_(f_int *handle,
		   f_int *gel,
		   char *cutoff,
		   f_int *lcutoff,
		   f_int *ok,
		   /* fortran implicits */
		   f_implicit l_cutoff)
/*
 * This routine is for backwards compatablity only
 */
{
    GapIO *io;

    if ( (io = io_handle(handle)) == NULL) f_proc_return();
    if (0 == cgetext(io, *gel, cutoff, (int *)lcutoff))
	*ok = 1;
    else
	*ok = 0;

    f_proc_return();
}


/*
 * shorted an extension by 'shorten_by' characters
 * Returns:
 *    0 - modification successful
 *  !=0 - an error has occurred
 */
int modext(GapIO *io, int gel, int shorten_by)
{
    return io_mod_extension(io,gel,shorten_by);
}

/*************************************************************
 * Read and Write contig annotation list. Used by fortran from
 * within the join code
 *************************************************************/
f_proc_ret getctg_(f_int *HANDLE, f_int *contig, f_int *anno) {
    GapIO *io;
    GContigs c;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    GT_Read(io, arr(GCardinal, io->contigs, *contig-1),
	    &c, sizeof(c), GT_Contigs);

    *anno = c.annotations;

    f_proc_return();
}

f_proc_ret putctg_(f_int *HANDLE, f_int *contig, f_int *anno) {
    GapIO *io;
    GContigs c;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    GT_Read(io, arr(GCardinal, io->contigs, *contig-1),
	    &c, sizeof(c), GT_Contigs);

    c.annotations = *anno;

    GT_Write(io, arr(GCardinal, io->contigs, *contig-1),
	     &c, sizeof(c), GT_Contigs);

    f_proc_return();
}

/*
 * delete a contig, contig_num from database and reorder order array.
 *
 * Sorry it's so verbose - but this is more confusing than it originally
 * appears!
 */
int io_delete_contig(GapIO *io, int contig_num) {
    int i;
    GContigs c;
    GCardinal *order = ArrayBase(GCardinal, io->contig_order);
    int el;
    reg_number rn;
    reg_delete rd;
    Array ar;


    /* remove annotations for contig contig_num */
    remove_contig_tags(io, contig_num, 0, 0);

    /* remove notes for contig contig_num */
    contig_read(io, contig_num, c);
    delete_note_list(io, c.notes);

    /* move the last contig into the space left by the deleted contig */
    contig_read(io, NumContigs(io), c);
    contig_write(io, contig_num, c);
    
    /* update the 3 fortan contig arrays */
    io_clnbr(io, contig_num) = c.left;
    io_crnbr(io, contig_num) = c.right;
    io_clength(io, contig_num) = c.length;

    /* Update the 1st note's "prev" pointer */
    if (c.notes) {
	GNotes n;
	note_read(io, c.notes, n);
	n.prev = contig_num;
	note_write(io, c.notes, n);
    }
      
    /*
     * Find contig_num in our order array and shuffle the array down
     * to fill this hole.
     *
     * Assuming a consistent ordering, upon termination of the loop 'el'
     * will hold the element in order for contig_num.
     */
    for (el = 0; el < NumContigs(io); el++) {
	if (order[el] == contig_num)
	    break;
    }
    
    /*
     * We've renumber the now last contig to be our old deleted contig.
     * In the case where contig_num == NumContigs(io) then we've simply
     * renumbered contig_num to contig_num - which is to do nothing.
     *
     * To keep the order up to date we now look for the highest contig.
     * and replace it's element in the order array with contig_num
     */
    for (i = 0; i < NumContigs(io); i++) {
	if (order[i] == NumContigs(io)) {
	    order[i] = contig_num;
	    break;
	}
    }

    /*
     * We now have shuffle down to fill in where contig_num used to be in
     * our order array
     */
    memmove(&order[el], &order[el+1], sizeof(int) * (NumContigs(io)-el-1));

    NumContigs(io)--;
    ArrayDelay(io, io->db.contig_order, io->db.Ncontigs, io->contig_order);
    DBDelayWrite(io);
    flush2t(io);
    NumContigs(io)++;

    /* Notify the deletion */
    rd.job = REG_DELETE;
    contig_notify(io, contig_num, (reg_data *)&rd);

    /* Notify our change of contig number */
    rn.job = REG_NUMBER_CHANGE;
    rn.number = contig_num;
    contig_notify(io, NumContigs(io), (reg_data *)&rn);

    /*
     * Decrement number of contigs.
     *
     * We need to do this after the REG_NUMBER_CHANGE request as plots such
     * as find_repeats do a search for the registered result for CONTIG_SEL.
     * If we've decremented NumContigs previously, and we only have two
     * contigs, the last of which is the one being renamed (we joined the first
     * to the second, and are then renaming the second), then it fails to
     * find the registered result as it only looks at contig 1.
     * However it also needs temporarily doing before writing out the
     * new GT_Database data otherwise the disk copy is not updated correctly.
     * You have been warned!
     */
    NumContigs(io)--;

    /*
     * Move our registration arrays down. We shuffle the deleted contig to the
     * hole made and set the number of registered items in it to zero.
     */
    ar = io_reg(io, contig_num);
    memcpy(&io_reg(io, contig_num), &io_reg(io, NumContigs(io)+1),
           sizeof(io_reg(io, contig_num)));
    memcpy(&io_cursor(io, contig_num), &io_cursor(io, NumContigs(io)+1),
           sizeof(io_cursor(io, contig_num)));
    io_reg(io, NumContigs(io)+1) = ar;
    io_Nreg(io, NumContigs(io)+1) = 0;
    io_cursor(io, NumContigs(io)+1) = NULL;

    return 0;
} 

f_proc_ret remcnl_(f_int *relpg, f_int *lngthg, f_int *lnbr, f_int *rnbr,
		      f_int *ngels, f_int *ncontigs, f_int *idbsize, 
		      f_int *contig_num, f_int *handle) 
{
    GapIO *io;

/*
      SUBROUTINE REMCNL(RELPG,LNGTHG,LNBR,RNBR,NGELS,NCONTS,IDBSIZ,
     +REMME,IDEVR)
*/

    if ( (io = io_handle(handle)) == NULL) f_proc_return();

    /* Yuk! Assumed to be done as part of the old WRITRN call by REMCON */
    NumReadings(io) = *ngels;

    io_delete_contig(io, *idbsize - *contig_num);
    *ncontigs = NumContigs(io);

    f_proc_return();
}

/*
 * Moves contig number 'from' to adjacent (right) to contig 'to'.
 */
f_proc_ret movec_(f_int *handle, f_int *from, f_int *to) 
{
    GapIO *io;
    int from_pos=0, to_pos=0, i;
    GCardinal *order;

    if ( (io = io_handle(handle)) == NULL) f_proc_return();

    /* Find positions in order of the two contigs */
    order = ArrayBase(GCardinal, io->contig_order);
    for (i = 0; i < NumContigs(io); i++) {
	if (order[i] == *from)
	    from_pos = i;
	if (order[i] == *to)
	    to_pos = i;
    }

    ReOrder(io, order, from_pos, to_pos+1);
    ArrayDelay(io, io->db.contig_order, io->db.Ncontigs, io->contig_order);
    flush2t(io);

    f_proc_return();
}


/*
 * Deallocate a reading
 */
f_proc_ret delgel_(f_int *HANDLE, f_int *N) {
    GapIO *io;
    GReadings r;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    /* remove notes */
    gel_read(io, *N, r);
    delete_note_list(io, r.notes);

    (void)io_deallocate_reading(io, *N);

    f_proc_return();
}
