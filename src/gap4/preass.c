#include <stdio.h>
#include <stdlib.h>

#include "misc.h"
#include "IO.h"
#include "seqInfo.h"
#include "array.h"
#include "io_utils.h"
#include "tagUtils.h"
#include "scf_extras.h"
#include "preass.h"
#include "text_output.h"
#include "clones.h"
#include "dbcheck.h"

struct reads_t {
    int position;
    int gel;
};

int add_reading(GapIO *io, SeqInfo *si, int contig, int position, int sense);

/*
 * Our local reading sort function. Sorts on position.
 */
int sort_reads(const void *r1, const void *r2) {
    return ((struct reads_t *)r1)->position - ((struct reads_t *)r2)->position;
}

/*
 * Adds a single preassembled contig to the database.
 * Returns 0 for success, -1 for failure.
 */
int load_preassembled(GapIO *io, 
		      int num,
		      char **reading_array) 
{
    int leftpos = 0;
    struct reads_t *reads;
    int i, j, left, right, contig, length;
    GReadings r;
    GContigs c;
    int count = 0, failed = 0;

    invalidate_rnumtocnum(io, 1);

    /* ---- Check if we've enough room in the database ---- */
    if (num + 1 + (NumReadings(io) + NumContigs(io) + 2) >=
	io->db.actual_db_size) {
	verror(ERR_FATAL, "enter_preassembled",
	      "Not enough free database slots - aborting");
	return -1;
    }


    /* ---- Allocate tmp arrays ---- */
    if (NULL == (reads = (struct reads_t *)xcalloc(num,
						   sizeof(struct reads_t))))
	return -1;


    /* ---- Create a contig ---- */
    contig = NumContigs(io) + 1;
    vmessage("Creating contig\n");
    if (-1 == io_init_contig(io, contig)) {
	xfree(reads);
	return -1;
    }
    UpdateTextOutput();


    /* ---- Add reads to the database ---- */
    for (i = 0; i < num; i++) {
	SeqInfo *si;
	char *cp, *n;
	int pos, sense = 0;
	
	n = reading_array[i];
	/* n = read_fofn(fp); */
	vmessage("Adding reading %s\n", n);
	UpdateTextOutput();

	if (NULL == (si = read_sequence_details(n, 1))) {
	    verror(ERR_WARN, "enter_preassembled",
		  "Failed to enter - couldn't process exp. file");
	    failed++;
	    continue;
	}

	if (exp_Nentries(si->e, EFLT_PC) == 0 ||
	    NULL == (cp = exp_get_entry(si->e, EFLT_PC))) {
	    freeSeqInfo(si);
	    verror(ERR_WARN, "enter_preassembled", 
		   "Failed to enter - no gel position information");
	    failed++;
	    continue;
	}

	pos = atoi(cp);
	if (exp_Nentries(si->e, EFLT_SE) &&
	    (cp = exp_get_entry(si->e, EFLT_SE)))
	    sense = atoi(cp);

	j = add_reading(io, si, contig, pos, sense);

	if (j > 0) {
	    reads[i].position = pos;
	    reads[i].gel = j;
	} else {
	    failed++;
	}

	freeSeqInfo(si);
    }

    /* ---- Sort readings in positional order ---- */
    qsort(reads, num, sizeof(struct reads_t), sort_reads);
    

    /* ---- Link readings together to form a contig ---- */
    left = 0;
    length = 0;
    vmessage("Linking readings\n");
    UpdateTextOutput();
    for (i = 0; i < num; left = reads[i].gel, i++) {
	if (!reads[i].gel)
	    continue;

	count++;

	gel_read(io, reads[i].gel, r);

	r.left = left;
	r.right = (i < num - 1) ? reads[i+1].gel : 0;

	/* Shift left if readings mark a region of a contig */
	if (!leftpos)
	    leftpos = reads[i].position;
	r.position -= leftpos - 1;

	if (r.position + r.sequence_length > length)
	    length = r.position + r.sequence_length;

	gel_write(io, reads[i].gel, r);
    }


    /* ---- Add contig details ---- */
    vmessage("Linking contig\n");
    UpdateTextOutput();
    GT_Read(io, arr(GCardinal, io->contigs, contig-1),
	    &c, sizeof(c), GT_Contigs);
    
    /* Possibly there is no contig (eg all the readings were previously
     * entered. We'd best check for that.
     */

    for (left = right = i = 0; i < num; i++) {
	if (reads[i].gel) {
	    right = reads[i].gel;

	    if (!left)
		left = reads[i].gel;
	}
    }

    if (left) {
	c.left = left;
	c.right = right;
	c.length = length - 1;

	GT_Write(io, arr(GCardinal, io->contigs, contig-1),
		 &c, sizeof(c), GT_Contigs);
    } else {
	NumContigs(io)--;
	DBDelayWrite(io);
    }

    /* ---- Tidy up ---- */
    xfree(reads);

    vmessage("\n%4d sequences processed\n", num);
    vmessage("%4d sequences entered into database\n\n", num - failed);
    UpdateTextOutput();

    invalidate_rnumtocnum(io, 0);

    return 0;
}


/*
 * Adds a single reading to the database. No linking of neighbours or update
 * of contig lengths is done, however the contig needs to be known in
 * order to add any consensus tags.
 */
int add_reading(GapIO *io, SeqInfo *si, int contig, int position, int sense) {
    char *t, *seq;
    GReadings r;
    int length, start, end;
    int geln, i;
    int2 *opos = NULL;
    int1 *conf = NULL;


    /* ---- Check for existance already ---- */
    if (get_gel_num(io, read_sequence_name(si), GGN_NAME) > 0) {
	verror(ERR_WARN, "enter_preassembled",
	       "ERROR!!! Reading already exists in database");
	return -1;
    }

    /* ---- Size checks ---- */
    length = si->length;
    start = si->start;
    end = si->end;

    /* ---- Initialise the reading ---- */
    geln = NumReadings(io)+1;
    if (-1 == io_init_reading(io, geln))
	return -1;

    gel_read(io, geln, r);
    seq = exp_get_entry(si->e, EFLT_SQ);


    /* ---- Original positions & Confidence values ---- */
    opos = (int2 *)xmalloc(length * sizeof(int2));
    if (!opos)
	return -1;
    if (si->origpos) {
	memcpy(opos, si->origpos, sizeof(int2) * length);
    } else {
	int i;

	for (i = 0; i < length; i++)
	    opos[i] = i+1;
    }

    conf = (int1 *)xmalloc(length);
    if (!conf) {
	xfree(opos);
	return -1;
    }
    if (si->confidence) {
	memcpy(conf, si->confidence, sizeof(int1) * length);
    } else {
	/* Read from SCF, CTF or ZTR file */
	if (0 != get_read_conf(si->e, length, opos, conf)) {
	    int i;

	    for (i = 0; i < length; i++)
		conf[i] = 99;
	}
    }

    /* ---- Add the reading name ---- */
    if (t = read_sequence_name(si)) {
	if (-1 == (r.name = allocate(io, GT_Text))) {
	    freeSeqInfo(si);
	    xfree(opos);
	    xfree(conf);
	    return -1;
	}

	if (-1 == TextWrite(io, r.name, t, strlen(t)+1)) {
	    freeSeqInfo(si);
	    xfree(opos);
	    xfree(conf);
	    return -1;
	}
	io_wname(io, geln, t);
    }


    /* ---- Gel tags ---- */
    if (gel_write(io, geln, r)) {
	verror(ERR_FATAL, "enter_preassembled",
	       "Problem writing reading to database");
	freeSeqInfo(si);
	xfree(opos);
	xfree(conf);
	return -1;
    }

    for(i = 0; i < exp_Nentries(si->e, EFLT_TG); i++) {
	create_tag_for_gel(io, geln, si->length, 
			   arr(char *, si->e->entries[EFLT_TG], i),
			   NULL, 0, NULL, 0);
    }


    /* ---- Consensus tags ---- */
    for (i=0; i < exp_Nentries(si->e, EFLT_TC); i++) {
	char *comment, *tag;
	char type[5];
	int start, end, strand;

	tag = arr(char *, si->e->entries[EFLT_TC], i);
	if (NULL == (comment = (char *)xmalloc(strlen(tag))))
	    continue;

	if (-1 == tag2values(tag, type, &start, &end, &strand, comment))
	    continue;

	if (sense == 0) {
	    start += position - si->start - 1;
	    end   += position - si->start - 1;
	} else {
	    int len = end - start;

	    start = position + si->end - end - 1;
	    end   = start + len;
	}

	type[4] = '\0';

	insert_NEW_tag(io, (tag_id)-contig, start, end-start+1, type, comment,
		       strand);

	xfree(comment);
    }


    /* ---- SVEC/CVEC tags derived from the SL, SR and CS lines ---- */
    if (exp_Nentries(si->e, EFLT_SL)) {
	int start = 1, len;
	
	len = atoi(exp_get_entry(si->e, EFLT_SL));

	insert_NEW_tag(io, (tag_id)geln, start, len, "SVEC", NULL, 0);
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
	    
	    insert_NEW_tag(io, (tag_id)geln, start, end - start + 1, "SVEC",
			   NULL, 0);
	}
    }

    if (exp_Nentries(si->e, EFLT_CS)) {
	int start, end;
	
	exp_get_rng(si->e, EFLT_CS, &start, &end);

	insert_NEW_tag(io, (tag_id)geln, start, end - start + 1, "CVEC", NULL, 0);
    }

    if (gel_read(io, geln, r)) {
	verror(ERR_FATAL, "enter_preassembled",
	       "Problem loading reading from database");
	freeSeqInfo(si);
	xfree(opos);
	xfree(conf);
	return -1;
    }

    
    /* ---- Complement if needed ---- */
    r.sense = sense;
    if (sense == GAP_SENSE_REVERSE) {
	io_complement_seq(&length, &start, &end, seq, conf, opos);
    }

    
    /* ---- Reading positions and lengths ---- */
    r.position = position;
    r.sequence_length = end - start - 1;

    
    /* ---- Write GReadings structure once more ---- */
    if (gel_write(io, geln, r)) {
	verror(ERR_FATAL, "enter_preassembled",
	       "Problem writing reading to database");
	freeSeqInfo(si);
	xfree(opos);
	xfree(conf);
	return -1;
    }


    /* ---- Add sequence, length, start, end, opos, and conf ---- */
    if (io_write_seq(io, geln, &length, &start, &end, seq, conf, opos)) {
	verror(ERR_FATAL, "enter_preassembled",
	       "Problem writing sequence to database");
	freeSeqInfo(si);
	xfree(opos);
	xfree(conf);
	return -1;
    }


    /* ---- Trace information ---- */
    if (exp_Nentries(si->e,EFLT_LT) && exp_Nentries(si->e,EFLT_LN)) {
	i = io_write_rd(io, geln,
			exp_get_entry(si->e,EFLT_LN),
			strlen(exp_get_entry(si->e,EFLT_LN)),
			exp_get_entry(si->e,EFLT_LT),
			strlen(exp_get_entry(si->e,EFLT_LT)));
    } else {
	verror(ERR_WARN, "enter_preassembled",
	       "no trace filename and filetype information found");
	i = io_write_rd(io, geln,
			"unknown", strlen("unknown"),
			"unknown", strlen("unknown"));
    }

    if (i) {
	verror(ERR_FATAL, "enter_preassembled",
	       "Problem writing raw data information to database");
	freeSeqInfo(si);
	xfree(opos);
	xfree(conf);
	return -1;
    }

    /* ---- Everything else: clone & primer info ---- */
    add_seq_details(io, geln, si);

    xfree(opos);
    xfree(conf);

    return geln;
}

    

/* ------------------------ Fortran interface ------------------------------ */

void update_fortran_arrays(GapIO *io,
			   int *ngels, 
			   int *nconts,
			   int *idbsiz, 
			   int *relpg, 
			   int *lngthg,
			   int *lnbr, 
			   int *rnbr) 
{
    int i;
    GReadings r;
    GContigs c;
    
    *ngels = NumReadings(io);
    *nconts = NumContigs(io);

    for (i = 1; i <= *ngels; i++) {
	gel_read(io, i, r);
	io_relpos(io, i) = r.position;
	io_length(io, i) = (r.sense == GAP_SENSE_REVERSE)
	    ? -r.sequence_length
	    : +r.sequence_length;
	io_lnbr(io, i)   = r.left;
	io_rnbr(io, i)   = r.right;
    }

    for (i = 1; i <= *nconts; i++) {
	contig_read(io, i, c);
	io_clength(io, i) = r.length;
	io_clnbr(io, i)   = r.left;
	io_crnbr(io, i)   = r.right;
    }
}

int
pre_assemble(int handle,
	     int num_readings,
	     char **reading_array)
{
    GapIO *io;
    int ngels;
    int nconts;
    int idbsiz;
    int *relpg;
    int *lngthg;
    int *lnbr;
    int *rnbr;

    if (NULL == (io = io_handle(&handle))) return -1;

    /* initialise fortran arrays */
    idbsiz = io_dbsize(io);
    relpg = &io_relpos(io,1);
    lngthg = &io_length(io,1);
    lnbr  = &io_lnbr(io,1);
    rnbr = &io_rnbr(io,1);
    
    if (-1 == load_preassembled(io, num_readings, reading_array)) {
	verror(ERR_WARN, "enter_preassembled", "failed");
    }
    
    update_fortran_arrays(io, &ngels, &nconts, &idbsiz,
			  relpg, lngthg, lnbr, rnbr);
/*    
    dbchek_(handle, relpg, lngthg, lnbr, rnbr, idm, idbsiz, ngels, nconts,
	    &ierr);
*/
    if (db_check(io) != 0) {
	verror(ERR_FATAL, "enter_preassembled",
	       "The database is now inconsistent.\n"
	       "You may wish to revert to a copy or to disassemble the newly "
	       "assembled contig.");
    }

    flush2t(io);
    return 0;
} /* end pre_assemble */
		  
#if 0
f_proc_ret preass_(f_int *handle, f_int *ngels, f_int *nconts,
		   f_int *idbsiz, f_int *relpg, f_int *lngthg,
		   f_int *lnbr, f_int *rnbr, f_int *idm) {
    char fofn[1024];
    GapIO *io;
    f_int ierr;
    
    if (NULL == (io = io_handle(handle))) f_proc_return();

    if (gtstr("File of filenames", "", fofn, 1024) != 0) {
	f_proc_return();
    }

    if (-1 == load_preassembled(io, fofn)) {
	puts("Failed");
    }

    update_fortran_arrays(io, ngels, nconts, idbsiz,
			  relpg, lngthg, lnbr, rnbr);
    
    dbchek_(handle, relpg, lngthg, lnbr, rnbr, idm, idbsiz, ngels, nconts,
	    &ierr);

    if (ierr) {
	verror(ERR_FATAL, "enter_preassembled",
	       "The database is now inconsistent.\n"
	       "You may wish to revert to a copy or to disassmble the newly "
	       "assembled contig.");
    }

    f_proc_return();
}
#endif
