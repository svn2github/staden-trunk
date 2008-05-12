/*
 * File: IO2.h
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
 * Created: 23 February 1993
 * Updated:
 *
 */

#ifndef _IO2_H_
#define _IO2_H_





extern
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
		      );

extern
    int io_delete_seq(int  *length,	/* length of complete string */
		      int  *start,	/* start */
		      int  *end,	/* end */
		      char *seq,	/* complete sequence */
		      int1 *conf,	/* confidence vals */
		      int2 *opos,	/* original pos */
		      /**********/
		      int  pos,
		      int  Nbases
		      );

extern
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
		       int  diff_only, /* Only change those that differ */
		       int  conf_only  /* Only modify confidence  */
		       );




extern f_proc_ret insbas_(f_int *HANDLE,
			  f_int *N,
			  f_int *POS,
			  char *BASE,
			  /* fortran implicits */
			  f_implicit BASE_l);

extern f_proc_ret modbas_(f_int *HANDLE,
			  f_int *N,
			  f_int *POS,
			  char *BASE,
			  /* fortran implicits */
			  f_implicit BASE_l);

extern f_proc_ret delbas_(f_int *HANDLE,
			  f_int *N,
			  f_int *POS);

extern f_proc_ret getext_(f_int *handle,
			  f_int *gel,
			  char *cutoff,
			  f_int *lcutoff,
			  f_int *ok,
			  /* fortran implicits */
			  f_implicit l_cutoff);

extern int modext(GapIO *io, int gel, int shorten_by);

extern f_proc_ret getctg_(f_int *HANDLE, f_int *contig, f_int *anno);
extern f_proc_ret putctg_(f_int *HANDLE, f_int *contig, f_int *anno);

/*
 * Pads in the consensus - rewrite in C. Currently uses PADCON().
 */
int pad_consensus(GapIO *io, int contig, int pos, int npads);

/*
 * delete a contig, contig_num from database and reorder order array.
 */
int io_delete_contig(GapIO *io, int contig_num);

/*
 * Deletes a single base from a gel reading ensuring that tags are shifted
 * and/or extended appropriately.
 *
 * Returns 0 for success, -1 for error.
 */
int io_delete_base(GapIO *io, int gel, int pos);

#endif /*_IO2_H_*/
