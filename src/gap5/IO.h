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

#include <tg_gio.h>
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
//#define GAP_CACHE 0

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

#define DB_NAMELEN 40  /* size of records in AR file. See also db_namelen
			  in gap.tcl and the *40's in legacy.f */
#define F_NAMLEN  40  /* length for filenames - namelen param in fortran */

/* Required to enable usage with the Array macros */
typedef struct {
    char name[DB_NAMELEN+1];
} name_t;



#if 0
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
#endif

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


/*
 * Open all database files
 */

extern GapIO *open_db(char *project, char *version, int *status, int create,
	       int read_only);

/*
 * Delete all database files for a project/version
 */


extern int close_db(GapIO *io);



/*************************************************************
 * Lower level IO routines
 *************************************************************/

extern GapIO *io_handle(f_int *HANDLE);



#include "io_handle.h"
#include "io_utils.h"

#endif /*_IO_H_*/
