/*
 * File: IO.c
 * Version:
 *
 * Author: Simon Dear
 *	   MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: C replacements for FORTRAN versions of Staden Bap I/O routines
 *
 * Created: 23 February 1993
 * Updated:
 *
 */


#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>		/* IMPORT: access */
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <sys/types.h>
#include <time.h>
#include <pwd.h>
#include <tcl.h>

#include "misc.h"
#include "array.h"
#include "bitmap.h"
#include "licence.h"

/*
 * Low level gap server stuff
 */
#include "g-defs.h"		/* IMPORT: G_LOCK_EX */
#include "g-error.h"		/* IMPORT: gerr_print */

/*
 * Gap interface level stuff
 */
#include "gap-if.h"
#include "gap-init.h"
#include "gap-dbstruct.h"
#include "gap-create.h"
#include "gap-io.h"

#include "gap-error.h"

#include "IO.h"
#include "actf.h"
#include "edUtils.h"		/* IMPORT: DB_NAMELEN */
#include "fortran.h"
#include "FtoC.h"
#include "xalloc.h"
#include "io-reg.h"
#include "text_output.h"
#include "notes.h"

static void convert_db(GapIO *io, char *project, char *version);

/*************************************************************
 * Low-ish level IO routines
 *************************************************************/

int GT_Read(GapIO *io, int rec, void *buf, int len, GCardinal type_check)
/*
 * Read in a GAP database record
 */
{
    int err;

    err = GAP_READ(io->client,
		   arr(GView,io->views,rec),
		   buf,len,type_check,sizeof(GCardinal));

    if (err) {
	GAP_ERROR_FATAL("reading record %d", rec);
    }

    return err;
}


int GT_Write(GapIO *io, int rec, void *buf, int len, GCardinal type)
/*
 * Write a GAP database record
 */
{
    int err;

    BIT_SET(io->updaterecs, rec);
    err = GAP_WRITE(io->client,
		    arr(GView,io->views,rec),
		    buf,len,type,sizeof(GCardinal));

    if (err) {
	GAP_ERROR_FATAL("writing record %d", rec);
    }

    return err;
}

#if GAP_CACHE!=0
/*
 * Write a GReadings structure to the database and cache
 */
int GT_Write_cached(GapIO *io, int read, GReadings *r) {
    int err;
    int rec = arr(GCardinal, io->readings, read-1);

    BIT_SET(io->updaterecs, rec);
    err = GAP_WRITE(io->client,
		    arr(GView, io->views, rec),
		    (char *)r, sizeof(GReadings),
		    GT_Readings, sizeof(GCardinal));

    if (err) {
	GAP_ERROR_FATAL("writing record %d", rec);
    }

    memcpy(arrp(GReadings, io->reading, read-1), r, sizeof(GReadings));

    return err;
}
#endif


int TextRead(GapIO *io, int rec, char *buf, int len)
/*
 * Read in a text string
 */
{
    int err;
    GViewInfo vi;
    int view;
    int len2;

    /* find out length of text */
    view = arr(GView,io->views,rec);
    if (view == -INT_MAX)
	return GAPERR_NOT_FOUND;

    /* Read only as much as required to avoid the memset in readv_image_() */
    g_view_info(io->client, view, &vi);
    len2 = vi.used - sizeof(GCardinal);

    err = GAP_READ(io->client, view, buf, MIN(len,len2), GT_Text,sizeof(char));

    /* NULL terminate if room */
    if (len > len2)
	buf[len2] = 0;

    if (err) {
	GAP_ERROR_FATAL("reading text %d", rec);
    }

    return err;
}






int TextWrite(GapIO *io, int rec, char *buf, int len)
/*
 * Write a text string
 */
{
    int err;

    BIT_SET(io->updaterecs, rec);
    err = GAP_WRITE(io->client,
		    arr(GView,io->views,rec),
		    buf,strnlen(buf, len), GT_Text, sizeof(char));

    if (err) {
	GAP_ERROR_FATAL("writing text %d", rec);
    }

    return err;
}




int DataRead(GapIO *io, int rec, void *buf, int len, int size)
/*
 * Read in data
 */
{
    int err;

    err = GAP_READ(io->client,
		   arr(GView,io->views,rec),
		   buf,len,GT_Data,size);

    if (err) {
	GAP_ERROR_FATAL("reading data %d", rec);
    }

    return err;
}

void *DataAllocRead(GapIO *io, int rec, int size)
/*
 * Read in data. Size is the size of each element in data (1, 2 or 4).
 */
{
    void *buf;
    int len;
    int err;
    GViewInfo vi;
    int view;

    /* find out length of data */
    view = arr(GView,io->views,rec);
    if (view == -INT_MAX)
	return NULL;

    err = g_view_info(io->client,view,&vi);
    len = vi.used - sizeof(GCardinal) /*type*/;

    /* allocate */
    buf = (void *)xmalloc(len+1);

    if (buf != NULL) {
	err = GAP_READ(io->client,
		       arr(GView,io->views,rec),
		       buf,len,GT_Data,size);
    }

    if (err) {
	GAP_ERROR_FATAL("reading data %d", rec);
    }

    return buf;
}






int DataWrite(GapIO *io, int rec, void *buf, int len, int size)
/*
 * Write data
 */
{
    int err;

    BIT_SET(io->updaterecs, rec);
    err = GAP_WRITE(io->client,
		    arr(GView,io->views,rec),
		    buf,len,GT_Data,size);
    
    if (err) {
	GAP_ERROR_FATAL("writing data %d", rec);
    }

    return err;
}



Array ArrayRead(GapIO *io, int rec, int elements)
/*
 * Allocate and read in an array
 */
{
    Array a;
    int err;

    a = ArrayCreate(sizeof(GCardinal),elements);

    if (a == NULL)
	GAP_ERROR_FATAL("creating array");

    if (ArrayRef(a,elements-1) == NULL)
	GAP_ERROR_FATAL("resizing array");

    err = GAP_READ(io->client,
		   arr(GView,io->views,rec),
		   arrp(GCardinal,a,0),sizeof(GCardinal)*elements,GT_Array,sizeof(GCardinal));
    if (err) {
	GAP_ERROR_FATAL("reading array %d", rec);
    }

    return a;
}



int ArrayWrite(GapIO *io, int rec, int elements, Array a)
/*
 * Write an array
 */
{
    BIT_SET(io->updaterecs, rec);
    return GAP_WRITE(io->client,
		     arr(GView,io->views,rec),
		     arrp(GCardinal,a,0),sizeof(GCardinal)*elements,GT_Array,sizeof(GCardinal));

}

/* ARGSUSED */
int ArrayDelay(GapIO *io, int rec, int elements, Array a)
/*
 * Set an array for writing, but delay it for now. See checks in flush2t().
 */
{
    BIT_SET(io->updaterecs, rec);
    return 0;
}

/* Mark the GR_Database record for writing at the next flush */
int DBDelayWrite(GapIO *io) {
    BIT_SET(io->updaterecs, GR_Database);
    return 0;
}



Bitmap BitmapRead(GapIO *io, int rec, int elements)
/*
 * Allocate and read a bitmap
 */
{
    Bitmap b;
    int err;

    b = BitmapCreate(elements);

    if (b == NULL)
	GAP_ERROR_FATAL("creating bitmap");

    err = GAP_READ(io->client,
		   arr(GView,io->views,rec),
		   b->base,b->Nbitmap*CHR_ELE,GT_Bitmap,CHR_ELE);

    if (err) {
	GAP_ERROR_FATAL("reading bitmap %d", rec);
    }

    return b;
}






int BitmapWrite(GapIO *io, int rec, Bitmap b)
/*
 * write a bitmap
 */
{
    BIT_SET(io->updaterecs, rec);
    return GAP_WRITE(io->client,
		     arr(GView,io->views,rec),
		     b->base,b->Nbitmap*CHR_ELE,GT_Bitmap,CHR_ELE);
}




/*************************************************************
 * Read utilities
 *************************************************************/



char *SeqReadStatic(GapIO *io, GCardinal rec, GCardinal length)
/*
 * Read the whole of a sequence into a static buffer
 */
{
    static int seql = 0;
    static char *seq = NULL;
    int err;

    if (seql < length) {
	seql = length;
	if (seq == NULL)
	    seq = (char *)xmalloc(seql);
	else
	    seq = (char *)xrealloc(seq,seql);
    }
    err = TextRead(io,rec,seq,seql);

    return seq;
}





char *TextAllocRead(GapIO *io, int rec)
/*
 * Allocate a buffer big enough and read in a text string
 */
{
    char *buf;
    int len;
    int err;
    GViewInfo vi;
    int view;

    /* find out length of text */
    view = arr(GView,io->views,rec);
    if (view == -INT_MAX)
	return NULL;

    err = g_view_info(io->client,view,&vi);
    len = vi.used - sizeof(GCardinal) /*type*/;

    /* allocate */
    buf = (char *)xmalloc(len+1);

    if (buf != NULL) {
	GAP_READ(io->client, view, buf, len, GT_Text, sizeof(char));
	buf[len] = '\0';
    }

    return buf;
}






/*************************************************************
 * Routines to initialise new reading and contig records
 *************************************************************/

/*
 * Reallocates io buffers when the 'maxdb' database size needs increasing.
 * db_size is the new size; io holds the existing size. The contents of io
 * will be valid after calling this function, but we make an assumption that
 * other objects have not kept their own internal pointers to io (which
 * _should_ be true).
 *
 * Returns 0 for success, -1 for failure.
 */
static int change_db_size(GapIO *io, int db_size) {
    int old_size = io->db.actual_db_size;
    int *tmp;
    int nc;

    if (db_size < old_size)
	return 0;

    /* Don't allow marginal resizes as this is an expensive function */
    if (db_size < old_size * 2)
	db_size = old_size * 2;

    /* Realloc the 4 'fortran' arrays */
    if (NULL == (tmp = (int *)xrealloc(io->relpos, db_size * sizeof(int))))
	return -1;
    io->relpos = tmp;

    if (NULL == (tmp = (int *)xrealloc(io->length, db_size * sizeof(int))))
	return -1;
    io->length = tmp;

    if (NULL == (tmp = (int *)xrealloc(io->lnbr, db_size * sizeof(int))))
	return -1;
    io->lnbr = tmp;

    if (NULL == (tmp = (int *)xrealloc(io->rnbr, db_size * sizeof(int))))
	return -1;
    io->rnbr = tmp;

    /* Update dbsize and copy the 'contig' parts of the arrays down */
    nc = NumContigs(io);
    memcpy(&io->relpos[db_size-nc], &io->relpos[old_size-nc], nc*sizeof(int));
    memcpy(&io->length[db_size-nc], &io->length[old_size-nc], nc*sizeof(int));
    memcpy(&io->lnbr  [db_size-nc], &io->lnbr  [old_size-nc], nc*sizeof(int));
    memcpy(&io->rnbr  [db_size-nc], &io->rnbr  [old_size-nc], nc*sizeof(int));
    io->db.actual_db_size = io->db.maximum_db_size = db_size;

    maxdb = db_size; /* update global and tcl variable */

    return 0;
}

/*
 * Calls change_db_size if we're looking close to running out of "database"
 * size.
 * Returns 0 (success) or -1 (failure).
 */
static int check_db_size(GapIO *io) {
    if (NumReadings(io) + NumContigs(io) + 1 >= io_dbsize(io))
	return change_db_size(io, io_dbsize(io)*2);
    else
	return 0;
}


int allocate(GapIO *io,GCardinal type)
/*
 * Find a new record
 */
{
    int freerec;

    /* allocate a record */
    freerec = BitmapFree(io->freerecs);

    if (freerec<0)
	GAP_ERROR_FATAL("allocating free record (BitmapFree)");
    if (BitmapExtend(io->updaterecs, freerec+1))
	GAP_ERROR_FATAL("allocating updaterecs record (BitmapExtend)");
    if (BitmapExtend(io->tounlock, freerec+1))
	GAP_ERROR_FATAL("allocating tounlock record (BitmapExtend)");

    /* BitmapExtend clears items by default */
    BIT_SET(io->freerecs,freerec);
    io->freerecs_changed = 1;

    /* write db record */
    io->db.Nfreerecs = io->freerecs->Nbitmap;
    DBDelayWrite(io);

    /* extend views array*/
    if (freerec >= io->Nviews) {
	(void) ArrayRef(io->views,freerec);
	if (freerec > io->Nviews) {
	    /* Shouldn't happen, but just in case */
	    int i;
	    puts("Warning - skipping views");
	    for (i = io->Nviews; i < freerec; i++)
		arr(GView, io->views, i) = -INT_MAX;
	}
	io->Nviews = freerec+1;
    } else {
	if (arr(GView, io->views, freerec) != -INT_MAX)
	    GAP_ERROR_FATAL("locking an inuse record %d, view %d",
			    freerec, arr(GView, io->views, freerec));
    }

    /* lock record */
    arr(GView,io->views,freerec) = g_lock_N(io->client,GAP_DATABASE_FILE,freerec,G_LOCK_EX);

    if (arr(GView,io->views,freerec)==-1) {
	GAP_ERROR_FATAL("could not lock new record %d", freerec);
    }

    /* GT_Write(io,freerec,NULL,0,type); */

    return freerec;
}


int deallocate(GapIO *io, int rec) {
    int err = 0;

    if (!BIT_CHK(io->freerecs, rec)) {
	GAP_ERROR_FATAL("deallocating an already free record %d", rec);
    }

    BIT_SET(io->updaterecs, rec);
    /* remove record */
    err |= g_remove(io->client, arr(GView, io->views, rec));

    /* Mark record for unlocking at the next flush */
    BIT_SET(io->tounlock, rec);
    io->freerecs_changed = 1;

    if (err) {
	GAP_ERROR_FATAL("deallocate() failed");
	return 1;
    }

    return 0;
}


/**
 * Deallocates a reading. Note that this doesn't care which reading is
 * being deallocated - it just frees the space. Other parts of the program
 * require that the reading numbers be consecutive. These other parts deal
 * with making sure that this is so!
 *
 * This neither deletes notes or tags. See delete_note_list() and
 * remove_gel_tags() for this functionality.
 *
 * Returns 0 for success
 *         non-zero for error (1)
 */
int io_deallocate_reading(GapIO *io, int N) {
    GReadings r;
    int err = 0;

    gel_read(io, N, r);
    
    if (r.name) {
	cache_delete_read_name(io, N);
	err += deallocate(io, r.name);
    }

    update_rnumtocnum(io, N, 0);

    if (r.trace_name)
	err += deallocate(io, r.trace_name);

    if (r.trace_type)
	err += deallocate(io, r.trace_type);

    if (r.sequence)
	err += deallocate(io, r.sequence);

    if (r.confidence)
	err += deallocate(io, r.confidence);

    if (r.orig_positions)
	err += deallocate(io, r.orig_positions);

    return err;
}


/*
 * Initialises a new reading.
 */
int io_init_reading(GapIO *io,	/*  */
		   int N	/* record/gel number  */
		   )
{
    int i;
    int err;

    if (check_db_size(io)) {
	verror(ERR_FATAL, "io_init_reading", "Couldn't grow database");
	return -1;
    }

    if (N > NumReadings(io)) {
/*
	if (N != NumReadings(io)+1)
	    verror(ERR_WARN, "io_init_reading", "NumReadings = %d, N = %d",
		  NumReadings(io),N);
*/

#if GAP_CACHE!=0
	/*
	 * Extend reading and names array
	 */
	ArrayRef(io->reading, N);
#endif
#if GAP_CACHE!=0
	ArrayRef(io->read_names, N);
#endif

	/*
	 * reclaim space
	 */
	for (i=NumReadings(io)+1;i<=Nreadings(io) && i<=N ;i++) {
	    GReadings r;
	    /*
	     * Need to recover old record! Read, clear, and write.
	     */

	    gel_read(io, i, r);
	    memset(&r, 0, sizeof(r));
	    gel_write(io, i, r);
	    io_wname(io, i, "");
#if GAP_CACHE!=0
	    memset(arrp(GReadings, io->reading, i-1), 0, sizeof(r));
#endif
	    update_rnumtocnum(io, i, 0); /* initialise to unknown contig */
	}
	NumReadings(io) = N;


	/*
	 * resize arrays, ensure records are allocated and initialised
	 */
	if (N > Nreadings(io)) {
	    /* extend arrays */
	    (void)ArrayRef(io->readings,N-1);
	    for (i=Nreadings(io)+1;i<=N;i++) {
		int freerec;
		
		freerec = allocate(io,GT_Readings);
		arr(GCardinal,io->readings,i-1) = freerec;

		/* initialise record - only need to write the type */
		err = GT_Write(io,freerec,NULL,0,GT_Readings);
#if GAP_CACHE!=0
		memset(arrp(GReadings, io->reading, i-1), 0,
		       sizeof(GReadings));
#endif
		io_wname(io, i, "");
	    }
	    Nreadings(io) = N;
	}

	/* YUK! should write database record back */
	DBDelayWrite(io);
	/* write array */
	err = ArrayDelay(io, io->db.readings, io->db.Nreadings, io->readings);
    }

    return 0;
}


int io_init_contig(GapIO *io,	/*  */
		   int N	/* record/gel number  */
		   )
{
    int i;
    int err;

    if (check_db_size(io)) {
	verror(ERR_FATAL, "io_init_contig", "Couldn't grow database");
	return -1;
    }

/*
    if (isrd_()) {
	puts("Error! Attempted write in read-only mode");
	return -1;
    }
*/
    if (N > NumContigs(io)) {
/*
	if (N != NumContigs(io)+1)
	    verror(ERR_WARN, "io_init_contig", "NumContigs = %d, N = %d",
		  NumContigs(io),N);
*/

	/*
	 * reclaim space
	 */
	for (i=NumContigs(io)+1;i<=Ncontigs(io) && i<=N;i++) {
	    GContigs c;
	    /*
	     * Need to recover old record! Read, clear, and write.
	     */

	    GT_Read(io, arr(GCardinal, io->contigs, i-1),
		    &c, sizeof(c), GT_Contigs);
	    memset(&c, 0, sizeof(c));
	    GT_Write(io, arr(GCardinal, io->contigs, i-1),
		     &c, sizeof(c), GT_Contigs);

	    arr(GCardinal, io->contig_order, i-1) = i;
	    if (io_reg(io, i) == NULL)
		io_reg(io, i) = ArrayCreate(sizeof(contig_reg_t), 0);
	    io_Nreg(io, i) = 0;
	    if (io_cursor(io, i) != NULL)
		xfree(io_cursor(io, i));
	    io_cursor(io, i) = NULL;
	}
	NumContigs(io) = N;


	/*
	 * resize arrays, ensure records are allocated and initialised
	 */
	if (N > Ncontigs(io)) {
	    /* extend arrays */
	    (void)ArrayRef(io->contigs,N-1);
	    (void)ArrayRef(io->contig_order,N-1);
	    (void)ArrayRef(io->contig_reg,N);
	    (void)ArrayRef(io->contig_cursor,N-1);
	    for (i=Ncontigs(io)+1;i<=N;i++) {
		int freerec;
		
		freerec = allocate(io,GT_Contigs);
		arr(GCardinal,io->contigs,i-1) = freerec;
		

		/* initialise record - only need to write the type */
		err = GT_Write(io,freerec,NULL,0,GT_Contigs);

		arr(GCardinal, io->contig_order, i-1) = i;

		io_reg(io, i) = ArrayCreate(sizeof(contig_reg_t), 0);
		io_Nreg(io, i) = 0;
		io_cursor(io, i) = NULL;
	    }
	    Ncontigs(io) = N;
	    
	}
	/*
	 * Write contig order to disk.
	 */
	ArrayDelay(io, io->db.contig_order, N, io->contig_order);

	/* YUK! should write database record back */
	DBDelayWrite(io);
	/* write array */
	err = ArrayDelay(io, io->db.contigs, io->db.Ncontigs, io->contigs);

    }

    return 0;
}




int io_init_annotations(GapIO *io, /*  */
			int N	   /* */
			)
/*
 * YUK! Not nice... should add new records onto free list
 */
{
    int i;
    int err;

/*
    if (isrd_()) {
	puts("Error! Attempted write in read-only mode");
	return -1;
    }
*/
    if (N > Nannotations(io)) {
/*
	if (N != Nannotations(io)+1)
	    verror(ERR_WARN, "io_init_annotations",
		   "Nannotations = %d, N = %d", Nannotations(io),N);
*/

	/* extend arrays */
	(void)ArrayRef(io->annotations,N-1);

	for (i=Nannotations(io)+1;i<=N;i++) {
	    int freerec;

	    freerec = allocate(io,GT_Annotations);
	    arr(GCardinal,io->annotations,i-1) = freerec;

	    /* initialise record - only need to write the type */
	    err = GT_Write(io,freerec,NULL,0,GT_Annotations);
	}

	Nannotations(io) = N;
	DBDelayWrite(io);
	/* write array */
	err = ArrayDelay(io, io->db.annotations, io->db.Nannotations, io->annotations);
    }

    return 0;
}


int io_init_note(GapIO *io, /*  */
		 int N	   /* */
		 )
{
    int i;
    int err;

    if (N > Nnotes(io)) {
	/* extend arrays */
	(void)ArrayRef(io->notes,N-1);

	for (i=Nnotes(io)+1;i<=N;i++) {
	    int freerec;

	    freerec = allocate(io,GT_Notes);
	    arr(GCardinal,io->notes,i-1) = freerec;

	    /* initialise record - only need to write the type */
	    err = GT_Write(io,freerec,NULL,0,GT_Notes);
	}

	Nnotes(io) = N;
	DBDelayWrite(io);
	/* write array */
	err = ArrayDelay(io, io->db.notes_a, io->db.Nnotes, io->notes);

    }

    return 0;
}



/*************************************************************
 * Higher level C interface routines
 *************************************************************/




int io_read_seq(GapIO *io,	/*  */
		int N,		/* record/gel number  */
		int  *length,	/* length of complete string */
		int  *start,	/* start */
		int  *end,	/* end */
		char *seq,	/* complete sequence */
		int1 *conf,	/* confidence vals (NULL if not needed) */
		int2 *opos)	/* original pos (NULL if not needed) */
{
    int err;
    GReadings r;

    if (N>Nreadings(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid reading %d", N);
    }

    /* read record */
    err = gel_read(io, N, r);

    /* update */
    *length = r.length;
    *start = r.start;
    *end = r.end;

    if (r.sequence==0)
	memset(seq,'?',*length);
    else
	err = TextRead(io,r.sequence,seq,*length);

    if (conf) {
	if (r.confidence==0)
	    memset(conf, 0, *length * sizeof(int1));
	else
	    err = DataRead(io,r.confidence,conf,*length * sizeof(int1),
			   sizeof(int1));
    }

    if (opos) {
	if (r.orig_positions==0)
	    memset(opos, 0, *length * sizeof(int2));
	else
	    err = DataRead(io,r.orig_positions,opos,*length * sizeof(int2),
			   sizeof(int2));
    }

    return 0;
}


/*
 * Reads the sequence and some associated details.
 * The various returned parameters may be NULL if desired in which case they
 * will not be filled out (or allocated).
 * Memory is allocated by this function, but is the responsibility of the
 * caller to xfree.
 *
 * Arguments in:
 *	io	GapIO pointer
 * 	N	The sequence reading number
 *
 * Arguments out:
 *	length	Complete length of sequence
 *	start	Last base of 3' clip (counting from 1)
 *	end	First base of 5' clip (counting from 1)
 *	seqp	Pointer to allocated buffer for sequence
 *	confp	Pointer to allocated buffer for confidence values
 *      oposp	Pointer to allocated buffer for original positions
 *
 * Returns:
 * 	0 for success
 *	-1 for failure
 *	May not return (GAP_ERROR_FATAL call).
 */
int io_aread_seq(GapIO *io,	/*  */
		 int    N,	/* record/gel number */
		 int   *length,	/* length of complete string */
		 int   *start,	/* start */
		 int   *end,	/* end */
		 char **seqp,	/* complete sequence */
		 int1 **confp,	/* confidence vals (NULL if not needed) */
		 int2 **oposp,	/* original pos (NULL if not needed) */
		 int    extra)  /* any extra allocation size required */
{
    int err;
    GReadings r;
    char *seq = NULL;
    int1 *conf = NULL;
    int2 *opos = NULL;

    if (N>Nreadings(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid reading %d", N);
    }

    /* read record */
    err = gel_read(io, N, r);

    /* Allocate buffers */
    if (seqp)  seq  = (char *)xmalloc((r.length + extra) * sizeof(*seq));
    if (confp) conf = (int1 *)xmalloc((r.length + extra) * sizeof(*conf));
    if (oposp) opos = (int2 *)xmalloc((r.length + extra) * sizeof(*opos));

    if ((seqp && !seq) ||
	(confp && !conf) ||
	(oposp && !opos)) {
	if (seq)  {xfree(seq);  *seqp  = NULL;}
	if (conf) {xfree(conf); *confp = NULL;}
	if (opos) {xfree(opos); *oposp = NULL;}
	return -1;
    }

    if (seq) {
	if (r.sequence==0)
	    memset(seq, '?', r.length);
	else
	    err = TextRead(io, r.sequence, seq, r.length);
    }

    if (conf) {
	if (r.confidence==0)
	    memset(conf, 0, r.length * sizeof(int1));
	else
	    err = DataRead(io,r.confidence,conf, r.length * sizeof(int1),
			   sizeof(int1));
    }

    if (opos) {
	if (r.orig_positions==0)
	    memset(opos, 0, r.length * sizeof(int2));
	else
	    err = DataRead(io,r.orig_positions,opos, r.length * sizeof(int2),
			   sizeof(int2));
    }

    if (length) *length = r.length;
    if (start)  *start  = r.start;
    if (end)    *end    = r.end;
    if (seqp)   *seqp   = seq;
    if (confp)  *confp  = conf;
    if (oposp)  *oposp  = opos;

    return 0;
}







int io_write_seq(GapIO *io,	/*  */
		 int N,		/* record/gel number  */
		 int  *length,	/* length of complete string */
		 int  *start,	/* start */
		 int  *end,	/* end */
		 char *seq,	/* complete sequence */
		 int1 *conf,	/* confidence vals */
		 int2 *opos)	/* original pos */
{
    int err;
    GReadings r;

    if (N>Nreadings(io)) err = io_init_reading(io,N);
/* GBUG - ought to be  if (N>NumReadings(io)) err = io_init_reading(io,N); */

    /* read record */
    err = gel_read(io, N, r);

    /* update */
    r.length = *length;
    r.start = *start;
    r.end = *end;

    r.sequence_length = r.end - r.start -1;

    if (r.sense == GAP_SENSE_REVERSE) {
	io_length(io, N) = -r.sequence_length;
    } else {
	io_length(io, N) = r.sequence_length;
    }

    if (r.sequence==0) r.sequence = allocate(io,GT_Text);
    err = TextWrite(io,r.sequence,seq,*length);

    if (r.confidence==0) r.confidence = allocate(io,GT_Data);
    err = DataWrite(io,r.confidence,conf,*length * sizeof(int1), sizeof(int1));

    if (r.orig_positions==0) r.orig_positions = allocate(io,GT_Data);
    err = DataWrite(io,r.orig_positions,opos,*length * sizeof(int2), sizeof(int2));

    /* write record back */
    err = gel_write(io, N, r);

    return 0;
}








int io_write_rd(GapIO *io,	/*  */
		int N,		/* record/gel number  */
		char *file,	/* trace file */
		int filelen,
		char *type,	/* trace file type */
		int typelen)
{
    int err;
    GReadings r;

    if (N>Nreadings(io)) err = io_init_reading(io,N);

    /* read record */
    err = gel_read(io, N, r);

    /* update */
    if (r.trace_name==0)
	r.trace_name = allocate(io,GT_Text);

    if (r.trace_type==0)
	r.trace_type = allocate(io,GT_Text);

    if (get_licence_type() == LICENCE_DEMO) {
	if (r.trace_name < r.trace_type) {
	    int tmp = r.trace_type;
	    r.trace_type = r.trace_name;
	    r.trace_name = tmp;
	}
    }

    err = TextWrite(io,r.trace_name,file,filelen);
    err = TextWrite(io,r.trace_type,type,typelen);

    /* write record back */
    err = gel_write(io, N, r);

    return 0;
}

int io_read_rd(GapIO *io,	/*  */
	       int N,		/* record/gel number  */
	       char *file,	/* trace file */
	       int filelen,
	       char *type,	/* trace file type */
	       int typelen)
{
    int err;
    GReadings r;
    int ret = 0;

    if (N>Nreadings(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid reading %d", N);
    }

    /* read record */
    err = gel_read(io, N, r);

    /* update */
    if (r.trace_name==0) {
	memset(file,' ',filelen);
	ret = 1;
    } else
	err = TextRead(io,r.trace_name,file,filelen);

    if (r.trace_type==0) {
	strncpy(type, "ANY", typelen);
    } else
	err = TextRead(io,r.trace_type,type,typelen);

    return ret;
}


int io_read_free_annotation(GapIO *io, int *f)
{
    *f = io->db.free_annotations;
    return 0;
}




int io_write_free_annotation(GapIO *io, int *f)
{
    io->db.free_annotations = *f;
    DBDelayWrite(io);

    return 0;
}




int io_read_annotation(GapIO *io,
		       int N,
		       int *anno)
{
    int err;
    GReadings r;
    GContigs c;

    if (N < 0) {
	N = -N;
	if (N > Ncontigs(io)) {
	    (void)gaperr_set(GAPERR_NOT_FOUND);
	    GAP_ERROR_FATAL("invalid contig %d", N);
	    *anno = 0;
	    return 1;
	}

	/* read record */
	err = GT_Read(io,arr(GCardinal,io->contigs,N-1),
		      &c,sizeof(c),GT_Contigs);

	*anno = c.annotations;
    } else {
	if (N>Nreadings(io)) {
	    (void)gaperr_set(GAPERR_NOT_FOUND);
	    GAP_ERROR_FATAL("invalid reading %d", N);
	    *anno = 0;
	    return 1;
	}

	/* read record */
	err = gel_read(io, N, r);

	*anno = r.annotations;
    }

    return 0;
}


int io_write_annotation(GapIO *io,
			int N,
			int *anno)
{
    int err;
    GReadings r;
    GContigs c;

    if (N < 0) {
	N = -N;
	if (N > Ncontigs(io)) err = io_init_contig(io, N);

	/* read record */
	err = GT_Read(io,arr(GCardinal,io->contigs,N-1),
		      &c,sizeof(c),GT_Contigs);

	/* update */
	c.annotations = *anno;

	/* write record back */
	err = GT_Write(io,arr(GCardinal,io->contigs,N-1),
		       &c,sizeof(c),GT_Contigs);
    } else {
	if (N>Nreadings(io)) err = io_init_reading(io,N);
	
	/* read record */
	err = gel_read(io, N, r);

	/* update */
	r.annotations = *anno;
	
	/* write record back */
	err = gel_write(io, N, r);
    }

    return 0;
}

/*
 * Open all database files
 */
GapIO *open_db(char *project, char *version, int *status, int create,
	       int read_only) {
    GapIO *io;
    int err, i;
    GHeaderInfo h;
    char db_fn[1024];

    if (get_licence_type() == LICENCE_VIEWER) {
	read_only = 1;
	create = 0;
    }

    *status = OK;

    /* Check if database is available for read-write */
    if (!read_only) {
	sprintf(db_fn, "%s.%s", project, version);
	errno = 0;
	if (access(db_fn, R_OK | W_OK) != 0 && errno == EACCES)
	    read_only = 1;
	sprintf(db_fn, "%s.%s.aux", project, version);
	errno = 0;
	if (access(db_fn, R_OK | W_OK) != 0 && errno == EACCES)
	    read_only = 1;

	if (read_only)
	    *status = IO_READ_ONLY;
    }

    if (0 != (err = (actf_lock(read_only, project, version, create)))) {
	if (err != 5 && err != 3) {
	    *status = ERROR;
	    return NULL;
	} else {
	    vmessage("Opening database in read only mode instead.\n");
	    read_only = 1;
	    *status = IO_READ_ONLY;
	}
    }

    /* Create a database if requested */
    if (create) {
	int i;

	if (!gap_server_is_local()) {
	    (void)gaperr_set(GAPERR_TRUSTME);
	    GAP_ERROR_FATAL("cannot create new database on a remote server");
	    *status = ERROR;
	    return NULL;
	}

	/*
	 * Check if database exists, exit with status==2 if so.
	 */
	for (i = 0; i < GAP_FILES; i++) {
	    char fn[1024];

	    (void)gap_construct_file(project, file_list[i], version, fn);
	    if (!(access(fn, F_OK) == -1 && errno == ENOENT)) {
		*status = 2;
		return NULL;
	    }
	    
	    strcat(fn, ".aux");
	    if (!(access(fn, F_OK) == -1 && errno == ENOENT)) {
		*status = 2;
		return NULL;
	    }
	}

	if (gap_new_db(project,version,read_only)) {
	    /* not enough memory for GapIO */
	    (void)gaperr_set(GAPERR_NO_ERROR);
	    GAP_ERROR_FATAL("cannot create database");
	    *status = ERROR;
	    return NULL;
	}
    }
    
    if ( (io = (GapIO *)xmalloc(sizeof(GapIO))) == NULL ) {
	/* not enough memory for GapIO */
	*status = ERROR;
	return NULL;
    }

    /*
     * Start a new server
     */
    if ( (io->server = gap_open_server(project,version, read_only)) == NULL ){
	GAP_ERROR("cannot open database");
	*status = NO_FILE;
	return NULL;
    }

#if !defined(__MINGW32__) && !defined(_MSC_VER)
    /* Start logging */
    {
	char log_buf[256], *user;
	struct passwd *pw;
	pw = getpwuid(getuid());
	user = pw ? pw->pw_name : "unknown";
	sprintf(log_buf, "opening r%c... by %s(%d)",
		read_only ? 'o' : 'w',
		user, getuid());

	(void)gap_construct_file(project, file_list[0], version, db_fn);
	strcat(db_fn, ".log");
	log_file(get_licence_type() == LICENCE_FULL ? db_fn : NULL,
		 log_buf);
    }
#endif

    /*
     * Connect client
     */
    io->client = g_connect_client(io->server,
				  (GLock)(read_only ? G_LOCK_RO : G_LOCK_EX));
    if (io->client == NULL)
	GAP_ERROR_FATAL("cannot connect client");
    if (!read_only && io->client->generic.mode == G_LOCK_RO) {
	vmessage("Couldn't open database in read/write mode, opened in read only instead.\n");
	read_only = 1;
	*status = IO_READ_ONLY;
	actf_unlock(0, project, version);
    }

    /*
     * IMPORTANT NOTE
     *	All the code here assumes that we are the only person accessing the database.
     *	I need to put some thought into how a multi user database would be done.
     */
    err = g_lock_file_N(io->client,GAP_DATABASE_FILE);
    if (err) GAP_ERROR_FATAL("locking database file");



    /*
     * Lock first records in database - these are always present.
     */
    err = g_header_info(io->client, GAP_DATABASE_FILE, &h);	
    io->Nviews = h.num_records;
    io->views = ArrayCreate(sizeof(GView),io->Nviews);
    (void)ArrayRef(io->views,io->Nviews);
    for(i=0;i<=GR_Unknown;i++) {
	arr(GView,io->views,i) = g_lock_N(io->client,GAP_DATABASE_FILE,i,
					  (GLock)(read_only
						  ? G_LOCK_RO
						  : G_LOCK_EX));
	if (arr(GView,io->views,i) < 0)
	    GAP_ERROR_FATAL("cannot lock record %d", i);
    }

    /*
     * Read GDatabase structure
     */
    err = GT_Read(io, GR_Database, &io->db, sizeof(io->db),GT_Database);
    if (err) GAP_ERROR_FATAL("cannot read database record");

    /* dumpGDatabase(&io->db); */

    /*
     * Fix up max_gel_len.
     * This protects against corrupt databases where max_gel_len is still
     * set at 4096, but longer sequences have been added. Shouldn't happen
     * with our code - but CAFtools seems to do this.
     *
     * max_gel_len is now an indicator of the maximum sequence length of a
     * trace file and hence the default allocation size. However sequences
     * may legitimately be longer than this.
     */
    io->db.max_gel_len = GAP_READ_LEN;

    /*
     * Read in freerecs bitmap, intitialise other bitmaps
     */
    io->freerecs = BitmapRead(io,io->db.freerecs,io->db.Nfreerecs * BIT_ELE);
    io->freerecs_changed = 0;
    io->updaterecs = BitmapCreate(io->db.Nfreerecs * BIT_ELE);
    io->tounlock = BitmapCreate(io->db.Nfreerecs * BIT_ELE);

    /*
     * We lock every other used record
     */
    for(i=GR_Unknown+1;i<io->Nviews;i++) {
	if (BIT_CHK(io->freerecs, i)) {
	    arr(GView,io->views,i) = g_lock_N(io->client,GAP_DATABASE_FILE, i,
					      (GLock)(read_only
						      ? G_LOCK_RO
						      : G_LOCK_EX));
	    if (arr(GView,io->views,i) < 0)
		GAP_ERROR_FATAL("cannot lock record %d", i);
	} else {
	    arr(GView, io->views, i) = -INT_MAX;
	}
    }

    if (io->db.version == 0) { /* earlier tag positioning */
	*status += 999;
    }

    /*
     * if (io->db.version != GAP_DB_VERSION && io->db.version != 0)
     *	   GAP_ERROR_FATAL("invalid database","old format not supported");
     */

    /* Set maxdb to be at least 50% larger than the current database */
    if (maxdb < 1.5 * (Nreadings(io) + Ncontigs(io) + 1))
	maxdb = 1.5 * (Nreadings(io) + Ncontigs(io) + 1);

    io->db.actual_db_size = maxdb;
    io->db.maximum_db_size = maxdb;


    /*************************************************************
     * Read in contigs array
     *************************************************************/
    io->contigs = ArrayRead(io,io->db.contigs,io->db.Ncontigs);


    /*************************************************************
     * Read in readings array
     *************************************************************/
    io->readings = ArrayRead(io,io->db.readings,io->db.Nreadings);


    /*************************************************************
     * Read in annotations array
     *************************************************************/
    io->annotations = ArrayRead(io,io->db.annotations,io->db.Nannotations);


    /*************************************************************
     * Read in templates array
     *************************************************************/
    io->templates = ArrayRead(io,io->db.templates,io->db.Ntemplates);

    /*************************************************************
     * Read in clones array
     *************************************************************/
    io->clones = ArrayRead(io,io->db.clones,io->db.Nclones);


    /*************************************************************
     * Read in vectors array
     *************************************************************/
    io->vectors = ArrayRead(io,io->db.vectors,io->db.Nvectors);


    /*************************************************************
     *
     *************************************************************/
    err = g_unlock_file_N(io->client,GAP_DATABASE_FILE);
    if (err) GAP_ERROR_FATAL("unlocking database file");


    /*************************************************************
     * Read in contig order array
     *************************************************************/
    if (io->db.contig_order) {
	io->contig_order = ArrayRead(io,io->db.contig_order,io->db.Ncontigs);
    } else {
	/* Create and extend order array */
	if (read_only) {
	    io->contig_order = ArrayCreate(sizeof(GCardinal), 0);
	    ArrayRef(io->contig_order, io->db.Ncontigs-1);
	    
	    for (i = 0; i < io->db.Ncontigs; i++)
		arr(GCardinal, io->contig_order, i) = i+1;
	    
	} else {
	    if (-1 == (io->db.contig_order = allocate(io, GT_Array)))
		GAP_ERROR_FATAL("Initialising contig order array");
	    io->contig_order = ArrayCreate(sizeof(GCardinal), io->db.Ncontigs);
	    ArrayDelay(io, io->db.contig_order, io->db.Ncontigs,
		       io->contig_order);
	    
	    for (i = 0; i < io->db.Ncontigs; i++)
		arr(GCardinal, io->contig_order, i) = i+1;
	    
	    ArrayDelay(io, io->db.contig_order, io->db.Ncontigs,
		       io->contig_order);
	    DBDelayWrite(io);
	}
    }

    /*************************************************************
     * Read in notes array
     *************************************************************/
    io->notes = NULL; /* 03/02/99 johnt - mark notes as unused so
                         it isn't mistakenly freed later */
    if (io->db.notes_a) {
	io->notes = ArrayRead(io, io->db.notes_a, io->db.Nnotes);
	fix_notes(io);
    } else if (!read_only) {
	/* Create notes array */
	if (-1 == (io->db.notes_a = allocate(io, GT_Notes)))
	    GAP_ERROR_FATAL("Initialising notes array");
	io->db.Nnotes = 0;
	io->notes = ArrayCreate(sizeof(GCardinal), io->db.Nnotes);
	ArrayDelay(io, io->db.notes_a, io->db.Nnotes, io->notes);
	DBDelayWrite(io);
    }

    /* Store database name */
    {
	char *p, buf[DB_FILELEN+1];
	
	if (p = strrchr(project, '/'))
	    p++;
	else
	    p = project;

	sprintf(buf, "%s.%s", p, version);
	strncpy(io_name(io), buf, DB_FILELEN-1);
	io_name(io)[DB_FILELEN-1] = 0;
    }

    /*
     * Allocate the reading/contig arrays.
     */
    io->relpos = (int *)xcalloc(io->db.actual_db_size, sizeof(int));
    io->length = (int *)xcalloc(io->db.actual_db_size, sizeof(int));
    io->lnbr   = (int *)xcalloc(io->db.actual_db_size, sizeof(int));
    io->rnbr   = (int *)xcalloc(io->db.actual_db_size, sizeof(int));

#if GAP_CACHE!=0
    io->reading= ArrayCreate(sizeof(GReadings), NumReadings(io));
#if GAP_CACHE==1
    io->read_names = ArrayCreate(sizeof(name_t), NumReadings(io));
#else
    io->read_names = ArrayCreate(sizeof(Tcl_HashEntry *), NumReadings(io));
    Tcl_InitHashTable(&io->rname_hash, TCL_STRING_KEYS);
    Tcl_InitHashTable(&io->tname_hash, TCL_STRING_KEYS);
#endif

    /* Intialise reading no. to contig no. mapping to zero (unknown) */
    io->rnum2cnum = ArrayCreate(sizeof(int), NumReadings(io));
    io->cached_rnum2cnum = 1;
    for (i = 0; i < NumReadings(io); i++)
	arr(int, io->rnum2cnum, i) = 0;

    if (NULL == io->reading) {
	*status = ERROR;
	return NULL;
    }
#endif

    if (NULL == io->relpos ||
	NULL == io->length ||
	NULL == io->lnbr ||
	NULL == io->rnbr) {
	*status = ERROR;
	return NULL;
    }

    /*
     * Initialise them
     */
    for (i = 1; i <= NumReadings(io); i++) {
#if GAP_CACHE!=0
	GReadings *r = arrp(GReadings, io->reading, i-1);
#else
	GReadings ra, *r = &ra;
#endif


	GT_Read(io, arr(GCardinal, io->readings, i-1),
		r, sizeof(*r), GT_Readings);

	io_wname(io, i, "");
#if GAP_CACHE==2
	/* Force load into cache */
	get_read_name(io, i);
#endif

	io_relpos(io,i) = r->position;
	io_length(io,i) = (r->sense == GAP_SENSE_REVERSE)
	    ? -r->sequence_length : r->sequence_length;
	io_lnbr	 (io,i) = r->left; 
	io_rnbr	 (io,i) = r->right;
    }

    for (i = 1; i <= NumContigs(io); i++) {
	GContigs c;

	GT_Read(io, arr(GCardinal, io->contigs, i-1),
		&c, sizeof(c), GT_Contigs);
	
	io_clength(io,i) = c.length;
	io_clnbr  (io,i) = c.left;
	io_crnbr  (io,i) = c.right;
    }

    contig_register_init(io);

#if GAP_CACHE == 2
    /* Hash template names */
    for (i = 1; i <= Ntemplates(io); i++) {
	char *name = get_template_name(io, i);
	cache_template_name(io, i, name);
    }
#endif

    /*
     * Check we are licenced to use this database.
     * We are allowed to open new databases and we're allowed to open
     * databases consisting of "valid" sequences. These may be edited though
     * so it's not easy to spot them...
     */
    if (get_licence_type() == LICENCE_DEMO) {
	int i;
	for (i = 1; i <= NumReadings(io); i++) {
	    GReadings r;
	    gel_read(io, i, r);
	    /* Obscure? :-) */
	    if (!r.trace_name|(r.trace_name<=r.trace_type)|!r.trace_type) {
		close_db(io);
		viewer_mode();
		return open_db(project, version, status, create, read_only);
	    }
	}
    }

    if (*status >= 999) {
	convert_db(io, project, version);
	*status -= 999;
    }

    flush2t(io); /* Incase we've changed things */
    log_file(NULL, "...opened");
    
    /* Execute any OPEN and RAWD notes */
    execute_database_notes(io, "OPEN");
    process_rawdata_note(io);

    return io;
}

/*
 * Deletes a database.
 * Returns -1 for failure, 0 for success.
 */
int del_db(char *project, char *version) {
    int i;
    char fn[1024];

    for (i = 0; i < GAP_FILES; i++) {

	(void)gap_construct_file(project, file_list[i], version, fn);

	if (-1 == remove(fn) ||
	    (strcat(fn, ".aux"), -1 == remove(fn))) {
	    verror(ERR_FATAL, "del_db", "Failed to remove old database");
	    return -1;
	}
    }

    return 0;
}

/*
 * Closes a database.
 * Returns -1 for failure, 0 for success.
 */
int close_db(GapIO *io) {
    int i, err, ro;

    flush2t(io);

    /* Execute any CLOS notes */
    execute_database_notes(io, "CLOS");

    /* Destroy registration lists */
    contig_register_destroy(io);

    log_file(NULL, "closing...");

    ro = io->client->generic.mode == G_LOCK_RO ? 1 : 0;

    /* Unlock everything */
    err = g_lock_file_N(io->client, GAP_DATABASE_FILE);

    for (i=0; i<io->Nviews; i++) {
	if (BIT_CHK(io->freerecs, i))
	    err |= g_unlock(io->client, arr(GView, io->views, i));
    }

    err |= g_unlock_file_N(io->client, GAP_DATABASE_FILE);

    /* Disconnect client */
    if (g_disconnect_client(io->client)) {
	GAP_ERROR("problem disconnecting");
	return -1;
    }

    /*
     * Shutdown server (YUK! why do we need to do this?)
     */
    gap_shutdown_server(io->server);

    /*
     * Destroy views array
     */
    ArrayDestroy(io->views);

    /*
     * Free up data structures
     */
    ArrayDestroy(io->contigs);
    ArrayDestroy(io->readings);
    ArrayDestroy(io->annotations);
    ArrayDestroy(io->templates);
    ArrayDestroy(io->clones);
    ArrayDestroy(io->vectors);
    ArrayDestroy(io->contig_cursor);
    ArrayDestroy(io->contig_order);
#if GAP_CACHE!=0
    ArrayDestroy(io->reading);
#if GAP_CACHE==2
    Tcl_DeleteHashTable(&io->rname_hash);
    Tcl_DeleteHashTable(&io->tname_hash);
#endif
    ArrayDestroy(io->read_names);
#endif
    ArrayDestroy(io->notes); /* 03/02/99 johnt - free notes if present */
    BitmapDestroy(io->freerecs);
    BitmapDestroy(io->updaterecs);
    BitmapDestroy(io->tounlock);

    {
	char project[DB_FILELEN], *version, *p;

	if (p = strrchr(io_name(io), '.')) {
	    version=p+1;
	    strncpy(project, io_name(io), p-io_name(io));
	    project[p-io_name(io)] = 0;

	    (void)actf_unlock(ro, project, version);
	}
    }

    xfree(io->relpos);
    xfree(io->length);
    xfree(io->lnbr);
    xfree(io->rnbr);
    xfree(io);

    log_file("", "...closed");

    return err ? -1 : 0;
}

static void convert_db(GapIO *io, char *project, char *version) {
    char *vers = "~";
    GReadings r;
    GAnnotations a;
    int gel, anno;

    vmessage("Your database is in an old style.\n");
    vmessage("I will automatically backup the existing database and convert\n");
    vmessage("the current database for you.\n\n");
    
    vmessage("Copying database to \"%s.%s\"...", project, vers);
    UpdateTextOutput();

    if (-1 == cpdb(project, version, vers)) {
	verror(ERR_WARN, "convert_db",
	    "Database copy failed! WARNING: not converting");
	return;
    }

    vmessage("Done\nConverting database...");
    UpdateTextOutput();

    /* convert gel annotations - all others need not be modified */
    for (gel = 1; gel <= NumReadings(io); gel++) {
	gel_read(io, gel, r);

	anno = r.annotations;

	while (anno && GT_Read(io, arr(GCardinal, io->annotations, anno-1),
			       &a, sizeof(a), GT_Annotations) == 0) {

	    a.position += r.sense == GAP_SENSE_ORIGINAL
		? r.start : r.length - r.end + 1;

	    GT_Write(io, arr(GCardinal, io->annotations, anno-1),
		     &a, sizeof(a), GT_Annotations);

	    anno = a.next;
	}
    }

    io->db.version = GAP_DB_VERSION;
    DBDelayWrite(io);
    flush2t(io);

    vmessage("Done\n\n");

    return;
}

void flush2t(GapIO *io) {
/*
 * Flush the database
 * This sets a checkpoint at which we can restart from
 */
    int i;
    int err;

    if (io->freerecs_changed) {
	err = BitmapWrite(io, io->db.freerecs, io->freerecs);
	if (err)
	    GAP_ERROR_FATAL("writing freerecs bitmap (flushing)");
	io->freerecs_changed = 0;
    }

    /*
     * Flush everything
     */
    err = g_lock_file_N(io->client,GAP_DATABASE_FILE);
    if (err)
	GAP_ERROR_FATAL("locking database file (to flush)");
    if (BIT_CHK(io->updaterecs, io->db.contigs)) {
	ArrayWrite(io, io->db.contigs, io->db.Ncontigs,
		   io->contigs);
    }
    if (BIT_CHK(io->updaterecs, io->db.readings)) {
	ArrayWrite(io, io->db.readings, io->db.Nreadings,
		   io->readings);
    }
    if (BIT_CHK(io->updaterecs, io->db.annotations)) {
	ArrayWrite(io, io->db.annotations, io->db.Nannotations,
		   io->annotations);
    }
    if (BIT_CHK(io->updaterecs, io->db.templates)) {
	ArrayWrite(io, io->db.templates, io->db.Ntemplates,
		   io->templates);
    }
    if (BIT_CHK(io->updaterecs, io->db.clones)) {
	ArrayWrite(io, io->db.clones, io->db.Nclones,
		   io->clones);
    }
    if (BIT_CHK(io->updaterecs, io->db.vectors)) {
	ArrayWrite(io, io->db.vectors, io->db.Nvectors,
		   io->vectors);
    }
    if (BIT_CHK(io->updaterecs, io->db.notes_a)) {
	ArrayWrite(io, io->db.notes_a, io->db.Nnotes,
		   io->notes);
    }
    if (BIT_CHK(io->updaterecs, io->db.contig_order)) {
	ArrayWrite(io, io->db.contig_order, io->db.Ncontigs,
		   io->contig_order);
    }
    if (BIT_CHK(io->updaterecs, GR_Database)) {
	GT_Write(io,GR_Database,&io->db,sizeof(io->db),GT_Database);
    }
    for(i=0;i<io->Nviews;i++) {
	if (BIT_CHK(io->updaterecs, i)) {
	    BIT_CLR(io->updaterecs, i);
	    if (BIT_CHK(io->tounlock, i)) {
		err = g_unlock(io->client, arr(GView, io->views, i));
		/* Mark as free in freerecs, and clear the tounlock flag */
		BIT_CLR(io->freerecs, i);
		BIT_CLR(io->tounlock, i);
		/* Mark view as unused - needed for error checking only */
		arr(GView, io->views, i) = -INT_MAX;
	    } else {
		err = g_flush(io->client,arr(GView, io->views,i));
	    }
	    if (err)
		GAP_ERROR_FATAL("flushing database file, rec %d", i);
	}
    }
    err = g_unlock_file_N(io->client,GAP_DATABASE_FILE);
    if (err)
	GAP_ERROR_FATAL("unlocking database file (flushed)");

    return;
}

/*
 * Returns the length of the longest gel in the contig, or in all contigs
 * if contig is specified as zero.
 *
 * If 'clipped' is true then we only look at the used quality-clipped portion,
 * otherwise we consider the hidden data too.
 */
int find_max_gel_len(GapIO *io, int contig, int clipped) { 
    int cnum, cstart, cend;
    int max_length = 0;

    if (!io)
	return -1;

    if (contig) {
	cstart = cend = contig;
    } else {
	cstart = 1;
	cend = NumContigs(io);
    }

    if (clipped) {
	/* Quality clipped portion */
	for (cnum = cstart; cnum <= cend; cnum++) {
	    int rnum, rlen;
	    
	    for (rnum = io_clnbr(io, cnum); rnum; rnum = io_rnbr(io, rnum)) {
		rlen = ABS(io_length(io, rnum));
		if (max_length < rlen)
		    max_length = rlen;
	    }
	}
    } else {
	/* Full length sequence */
	for (cnum = cstart; cnum <= cend; cnum++) {
	    int rnum;
	    
	    for (rnum = io_clnbr(io, cnum); rnum; rnum = io_rnbr(io, rnum)) {
		GReadings r;
		gel_read(io, rnum, r);
		if (max_length < r.length)
		    max_length = r.length;
	    }
	}
    }

    return max_length;
}


