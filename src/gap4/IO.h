/*
 * File: IO.h
 * Version:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
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

#ifndef _IO_H_
#define _IO_H_

#include "g-defs.h"		/* IMPORT: G_LOCK_EX */

/*
 * Defining this as 0 disables caching of reading names and structures.
 * Defining this as 1 enables caching of reading names/structures, but it
 * is not possible to search them to quickly obtain numbers from names.
 * Defining this as 2 enables caches of reading names/structures as well as
 * template names. The names are cached in a hash table and so provides
 * very quick access via the get_gel_num (etc) functions.
 */
#define GAP_CACHE 2

/*
 * Error return values from open_db
 */
#define OK           0	/* Open was successful for the mode requested */
#define ERROR        1	/* Error! */
#define FILE_EXISTS  2	/* Database exists (create failed) */
#define NO_FILE      3	/* Couldn't find it */
#define IO_READ_ONLY 4	/* Open was successful, but at reduced access */

#include "os.h"

#include "array.h"
#include "bitmap.h"
#include "gap-dbstruct.h"
#include "gap-if.h"

#if GAP_CACHE==2
#include <tcl.h>
#endif

#define DB_FILELEN 256 /* size of the filenames on disk (project.version) */
		       /* See also select_contig.tcl CheckDBFilename */
#define FILE_NAME_LENGTH DB_FILELEN /* Length of trace file names */

#define DB_NAMELEN 40  /* size of records in AR file. See also db_namelen
			  in gap.tcl and the *40's in legacy.f */
#define F_NAMLEN  40  /* length for filenames - namelen param in fortran */

/* typedefs */

typedef struct {
    GapServer *server;		/* our server */
    GapClient *client;		/* ourselves */

    int Nviews;			/* number of locked views */
    Array views;		/* all locked views */

    GDatabase db;		/* main database record */
    Bitmap freerecs;		/* bitmap of unused */
    Array contigs;		/* list of contig */
    Array readings;		/* list of reading records */
    Array annotations;		/* list of annotation records */
    Array templates;		/* list of template records */
    Array clones;		/* list of clone records */
    Array vectors;		/* list of vector records */
    Array notes;		/* list of note records */

    int4 *relpos;		/* relpg[] */
    int4 *length;		/* length[] */
    int4 *lnbr;			/* lnbr[] */
    int4 *rnbr;			/* rnbr[] */

    char db_name[DB_FILELEN];	/* database "file.version" */

    Array contig_order;		/* order of contigs */
    Array contig_reg;		/* Registration arrays for each contig */

#if GAP_CACHE==1
    Array reading;		/* Array of GReading _structures_ */
    Array read_names;		/* Array of reading names, type char **/
#elif GAP_CACHE==2
    Array reading;		/* Array of GReading _structures_ */
    Array read_names;		/* Array of reading names, type Tcl_HashEntry*/
    Tcl_HashTable rname_hash;	/* Tcl hash table of reading names */
    Tcl_HashTable tname_hash;	/* Tcl hash table of template names */
#endif
    int freerecs_changed;	/* Whether to flush freerecs bitmap */
    Bitmap updaterecs;		/* bitmap of updated records */
    Bitmap tounlock;		/* bitmap of records to unlock at next flush */

    Array contig_cursor;	/* array (NumContigs) of cursor_t lists */
} GapIO;

/* Required to enable usage with the Array macros */
typedef struct {
    char name[DB_NAMELEN+1];
} name_t;




/*
 * Useful macros
 */
#define io_dbsize(IO)    ( (IO)->db.actual_db_size )
#define io_maxdbsize(IO) ( (IO)->db.maximum_db_size )
#define max_gel_len(IO)  ( (IO)->db.max_gel_len )
#define NumContigs(IO)   ( (IO)->db.num_contigs )
#define NumReadings(IO)  ( (IO)->db.num_readings )
#define Ncontigs(IO)     ( (IO)->db.Ncontigs )
#define Nreadings(IO)    ( (IO)->db.Nreadings )
#define Nannotations(IO) ( (IO)->db.Nannotations )
#define Nnotes(IO)       ( (IO)->db.Nnotes )
#define Ntemplates(IO)   ( (IO)->db.Ntemplates )
#define Nclones(IO)      ( (IO)->db.Nclones )
#define Nvectors(IO)     ( (IO)->db.Nvectors )
#define io_relpos(IO,g)  ( (IO)->relpos[(g)] )
#define io_length(IO,g)	 ( (IO)->length[(g)] )
#define io_lnbr(IO,g)	 ( (IO)->lnbr[(g)] )
#define io_rnbr(IO,g)	 ( (IO)->rnbr[(g)] )
#define io_clength(IO,c) ( (IO)->relpos[io_dbsize(IO) - (c)] )
#define io_clnbr(IO,c)	 ( (IO)->lnbr[io_dbsize(IO) - (c)] )
#define io_crnbr(IO,c)	 ( (IO)->rnbr[io_dbsize(IO) - (c)] )
#define io_name(IO)	 ( (IO)->db_name )
#define io_rdonly(IO)	 ( (IO)->client->generic.mode == G_LOCK_RO )

/*
 * Code for reading and writing reading names directly from the cached
 * and/or hashed copies.
 */
#if GAP_CACHE!=0
    #define io_rname(IO,g)   ( get_read_name((IO), (g)) )
    #define io_wname(IO,g,n) ( cache_read_name((IO),(g),(n)) )
#else
    #define io_rname(IO,g)   ( get_read_name((IO), (g)) )
    #define io_wname(IO,g,n)
#endif

extern int primer_type_arr[5][2];
extern int primer_type_guess_arr[5][2];
extern int strand_arr[5][2];

/*
 * When we have primer and strand info...
 */
/*
#define PRIMER_TYPE(r) (((r).primer == GAP_PRIMER_CUSTFOR && \
			 (r).strand == GAP_STRAND_REVERSE) \
			? GAP_PRIMER_CUSTREV : (r).primer)
*/
#define PRIMER_TYPE(r) (primer_type_arr[(r).primer][(r).strand])

/*
 * But we sometimes (often) have primer type as unknown, so we can guestimate
 * in these cases.
 */
/*
#define PRIMER_TYPE_GUESS(r) ((r).primer == GAP_PRIMER_UNKNOWN ? \
			      ((r).strand + GAP_PRIMER_CUSTFOR) : \
			      (((r).primer == GAP_PRIMER_CUSTFOR && \
				(r).strand == GAP_STRAND_REVERSE) \
			       ? GAP_PRIMER_CUSTREV : (r).primer))
*/
#define PRIMER_TYPE_GUESS(r) (primer_type_guess_arr[(r).primer][(r).strand])

/*
 * Computes strand from the primer information. At present, this may compute
 * the primer from the strand (but won't in the future). The result is
 * taken from primer in preference when primer and strand disagree.
 */
/*
#define STRAND(r) (1-((PRIMER_TYPE_GUESS(r)) & 1))
*/
#define STRAND(r) (strand_arr[(r).primer][(r).strand])


/*************************************************************
 * Low-ish level IO routines
 *************************************************************/

int GT_Read(GapIO *io, int rec, void *buf, int len, GCardinal type_check);
/*
 * Read in a GAP database record
 */

int GT_Write(GapIO *io, int rec, void *buf, int len, GCardinal type);
/*
 * Write a GAP database record
 */

/*
 * Write a GReadings structure to the database and cache
 */
int GT_Write_cached(GapIO *io, int read, GReadings *r);

int TextRead(GapIO *io, int rec, char *buf, int len);
/*
 * Read in a text string
 */

char *TextAllocRead(GapIO *io, int rec);
/*
 * Read in and alloc a text string
 */

int TextWrite(GapIO *io, int rec, char *buf, int len);
/*
 * Write a text string
 */

int DataRead(GapIO *io, int rec, void *buf, int len, int size);
void *DataAllocRead(GapIO *io, int rec, int size);
/*
 * Read in data
 */

int DataWrite(GapIO *io, int rec, void *buf, int len, int size);
/*
 * Write data
 */

Array ArrayRead(GapIO *io, int rec, int elements);
/*
 * Allocate and read in an array
 */

int ArrayWrite(GapIO *io, int rec, int elements, Array a);
/*
 * Write an array
 */

int ArrayDelay(GapIO *io, int rec, int elements, Array a);
/*
 * Set an array for writing, but delay it for now. See checks in flush2t().
 */

int DBDelayWrite(GapIO *io);
/* Mark the GR_Database record for writing at the next flush */

Bitmap BitmapRead(GapIO *io, int rec, int elements);
/*
 * Allocate and read a bitmap
 */

int BitmapWrite(GapIO *io, int rec, Bitmap b);
/*
 * write a bitmap
 */


/*************************************************************
 * Utilities
 *************************************************************/

extern f_proc_ret cloz2t_(f_int *HANDLE);
/*
 * Close all
 */

extern void flush2t(GapIO *io);
extern f_proc_ret flus2t_(f_int *HANDLE);
/*
 * Flush a database file
 */

extern f_proc_ret readr0_(f_int *HANDLE, f_int *IMSIZ, f_int *IASIZ, f_int *MXG, f_int *IDA);
/*
 * Read first record
 */

extern f_proc_ret writr0_(f_int *HANDLE, f_int *IMSIZ, f_int *IASIZ, f_int *MXG, f_int *IDA);
/*
 * Write first record
 */

extern f_proc_ret readrn_(f_int *HANDLE, f_int *NGELS, f_int *NCONTS);
/*
 * Read ngels and ncontigs
 */

extern f_proc_ret writrn_(f_int *HANDLE, f_int *NGELS, f_int *NCONTS);
/*
 * Write ngels and ncontigs
 */

extern f_proc_ret readg_(f_int *HANDLE, f_int *N, f_int *RELPG, f_int *LNGTHG, f_int *LNBR, f_int *RNBR);
/*
 * Read a gel line
 */

extern f_proc_ret writeg_(f_int *HANDLE, f_int *N, f_int *RELPG, f_int *LNGTHG, f_int *LNBR, f_int *RNBR);
/*
 * Write a gel line
 */

extern f_proc_ret readc_(f_int *HANDLE, f_int *N, f_int *LNGTHC, f_int *LGEL, f_int *RGEL);
/*
 * Write a contig line
 */

extern f_proc_ret writec_(f_int *HANDLE, f_int *N, f_int *LNGTHC, f_int *LGEL, f_int *RGEL);
/*
 * Write a contig line
 */

extern f_proc_ret readtg_(f_int *HANDLE, f_int *N, f_int *LPOS, f_int *LLEN, f_int *LCOM, f_int *LTYPE, f_int *NEXT, f_int *SENSE);
/*
 * Read a tag record
 */

extern f_proc_ret writtg_(f_int *HANDLE, f_int *N, f_int *LPOS, f_int *LLEN, f_int *LCOM, f_int *LTYPE, f_int *NEXT, f_int *SENSE);
/*
 * Write a tag record
 */

extern f_proc_ret readcc_(f_int *HANDLE, f_int *N, f_int *ICNT, f_int *NEXT, char *NOTE, f_implicit NOTE_l);
/*
 * Read a comment record
 */

extern f_proc_ret writcc_(f_int *HANDLE, f_int *N, f_int *ICNT, f_int *NEXT, char *NOTE, f_implicit NOTE_l);
/*
 * Write a comment record
 */

extern f_proc_ret readn_(f_int *HANDLE, f_int *N, char *NAME, f_implicit NAME_l);
/*
 * Read a name
 */

extern f_proc_ret writen_(f_int *HANDLE, f_int *N, char *NAME, f_implicit NAME_l);
/*
 * Write a name
 */

extern f_proc_ret readr_(f_int *HANDLE, f_int *N, f_int *A1, f_int *A2, f_int *A3, f_int *A4);
/*
 * No longer used
 */

extern f_proc_ret writer_(f_int *HANDLE, f_int *N, f_int *A1, f_int *A2, f_int *A3, f_int *A4);
/*
 * No longer used
 */

extern f_proc_ret readw_(f_int *HANDLE, f_int *N, char *GEL, f_int *MAXGEL, f_implicit GEL_l);
/*
 * Read a squence
 */

extern f_proc_ret writew_(f_int *HANDLE, f_int *N, char *GEL, f_int *MAXGEL, f_implicit GEL_l);
/*
 * Write a squence
 */

extern int readrd(f_int *HANDLE, int N, char *type, char *trace, int type_l, int trace_l);
/*
 * Rrace trace details
 */

extern f_proc_ret open2t_(
		   /* takes */
		   char *PROJECT, f_int *PLEN, char *VERSION, f_int *VLEN, f_int *JOB,
		   f_int *RELPG,
		   f_int *LNGTHG,
		   f_int *LNBR,
		   f_int *RNBR,
		   /* returns */
		   f_int *HANDLE, f_int *IOK,
		   f_int *READ_ONLY,
		   /* implicit args */
		   f_implicit PROJECT_l, f_implicit VERSION_l);
/*
 * Open all database files
 */

extern GapIO *open_db(char *project, char *version, int *status, int create,
	       int read_only);

extern f_proc_ret del2t_(
		  /* takes */
		  char *PROJECT, f_int *PLEN, char *VERSION, f_int *VLEN,
		  /* returns */
		  f_int *IOK,
		  /* implicit args */
		  f_implicit PROJECT_l, f_implicit VERSION_l);
/*
 * Delete all database files for a project/version
 */


extern int close_db(GapIO *io);


/*
 * Read template for a given gel
 */
extern f_int rdtmpl_(f_int *HANDLE, f_int *N);


extern char *SeqReadStatic(GapIO *io, GCardinal rec, GCardinal length);

int io_get_extension(GapIO *io,
		     int N,	      /* gel number */
		     char *seq,	      /* buffer for sequence */
		     int max_seq,     /* size of buffer */
		     /* returns */
		     int *length,     /* length returned */
		     int *complement);/* reading is complemented */

int io_mod_extension(GapIO *io,
		     int N,	      /* gel number */
		     int shorten_by); /* ammount to shorten by */

/*************************************************************
 * Lower level IO routines
 *************************************************************/

extern GapIO *io_handle(f_int *HANDLE);


extern
    int io_read_seq(GapIO *io,	/*  */
		    int N,	/* record/gel number  */
		    int *length, /* length of complete string */
		    int *start, /* start */
		    int *end,	/* end */
		    char *seq,	/* complete sequence */
		    int1 *conf,	/* confidence vals */
		    int2 *opos); /* original pos */


extern
    int io_write_seq(GapIO *io,	/*  */
		     int N,		/* record/gel number  */
		     int *length,	/* length of complete string */
		     int *start,	/* start */
		     int *end,	/* end */
		     char *seq,	/* complete sequence */
		     int1 *conf,	/* confidence vals */
		     int2 *opos);	/* original pos */


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
		 int2 **oposp,  /* original pos (NULL if not needed) */
		 int    extra);	/* any extra allocation size required */

extern
    int get_read_info(GapIO *io,
		      int N,
		      char *clone, int l_clone,
		      char *cvector, int l_cvector,
		      char *subclone, int l_subclone, /* aka template */
		      char *scvector, int l_scvector,
		      int *length,
		      int *insert_min,
		      int *insert_max,
		      int *direction, /* 0 = forward, 1 = reverse */
		      int *strands,
		      int *primer,
		      int *clone_id,
		      int *subclone_id,
		      int *cvector_id,
		      int *scvector_id);


extern int get_vector_info(GapIO *io, int vector_id,
			   char *vector, int l_vector);

extern int get_clone_info(GapIO *io, int clone_id,
			  char *clone, int l_clone,
			  char *cvector, int l_cvector,
			  int *cvector_id);

extern int get_subclone_info(GapIO *io,
			     int subclone_id,
			     char *clone, int l_clone,
			     char *cvector, int l_cvector,
			     char *subclone, int l_subclone, /* aka template */
			     char *scvector, int l_scvector,
			     int *insert_min,
			     int *insert_max,
			     int *strands,
			     int *clone_id,
			     int *cvector_id,
			     int *scvector_id);

extern int io_init_reading(GapIO *io,
			   int N);

extern int io_init_contig(GapIO *io,
			  int N);

extern int io_init_annotations(GapIO *io,
			       int N);

extern int io_init_note(GapIO *io,
			int N);

int io_read_seq(GapIO *io,	/*  */
		int N,		/* record/gel number  */
		int *length,	/* length of complete string */
		int *start,	/* start */
		int *end,	/* end */
		char *seq,	/* complete sequence */
		int1 *conf,	/* confidence vals (NULL if not needed) */
		int2 *opos);	/* original pos (NULL if not needed) */

extern int io_read_annotation(GapIO *io,
			      int N,
			      int *anno);

int io_write_annotation(GapIO *io,
			int N,
			int *anno);

extern int io_write_rd(GapIO *io,
		       int N,
		       char *file,
		       int filelen,
		       char *type,
		       int typelen);

extern int io_deallocate_reading(GapIO *io, int N);

extern int allocate(GapIO *io, GCardinal type);

extern int deallocate(GapIO *io, int rec);

/*
 * Deletes a database.
 * Returns -1 for failure, 0 for success.
 */
int del_db(char *project, char *version);

/*
 * Write a name. Returns 0 for success, -1 for failure.
 */
int write_rname(GapIO *io, int reading, char *name);

int io_read_rd(GapIO *io,	/*  */
	       int N,		/* record/gel number  */
	       char *file,	/* trace file */
	       int filelen,
	       char *type,	/* trace file type */
	       int typelen);

int io_complement_seq(int *length,	/* length of complete string */
		      int *start,	/* start */
		      int *end,	/* end */
		      char *seq,	/* complete sequence */
		      int1 *conf,	/* confidence vals */
		      int2 *opos);	/* original pos */

/*
 * Gets an extension and complements it if needed.
 * Cutlen holds the length of 'cutoff' both before (length allowed to
 * return) and after (length actually returned) the call.
 */
int cgetext(GapIO *io, int gel, char *cutoff, int *cutlen);

/*
 * shorted an extension by 'shorten_by' characters
 * Returns:
 *    0 - modification successful
 *  !=0 - an error has occurred
 */
int modext(GapIO *io, int gel, int shorten_by);

/*
 * Inserts a single base into a gel reading ensuring that tags are shifted
 * and/or extended appropriately.
 *
 * Returns 0 for success, -1 for error.
 */
int io_insert_base(GapIO *io, int gel, int pos, char base);

int io_read_free_annotation(GapIO *io, int *f);
int io_write_free_annotation(GapIO *io, int *f);

/*
 * Returns the length of the longest gel in the contig, or in all contigs
 * if contig is specified as zero.
 *
 * If 'clipped' is true then we only look at the used quality-clipped portion,
 * otherwise we consider the hidden data too.
 */
int find_max_gel_len(GapIO *io, int contig, int clipped);

/*
 * swap readings N and M
 */
int swap_read(GapIO *io, int M, int N);

#include "io_handle.h"
#include "io_utils.h"

#endif /*_IO_H_*/
