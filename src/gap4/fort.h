/*
 * File:
 * Version:
 *
 * Author:
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
#ifndef _FORT_H_
#define _FORT_H_

/*
 * When passing INTEGERs from Fortran to C we need to know the type they'll be
 * at the C end. We could define int_f as 'int *' (or 'long *') but this may
 * cause obscureness and unreadability of C.
 *
 * Also we define the extra argument given to C by fortran for the length of
 * an array passed. This may not necessarily be the same as int_f (either
 * now or the future). So we define it as a separate type (int_fl).
 */

#include "os.h"
#include "io_utils.h"

/*
 * ---------------------------------------------------------------------------
 * Prototypes for Fortran functions. This helps to establish when we pass the
 * wrong arguments to a function.
 *
 * Terminology:
 * name_l : name is the length of the string name (or an element in a string
 *          array?)
 * ---------------------------------------------------------------------------
 */

/* dbas.f */
f_proc_ret dbauto_(f_int *relpg,  f_int *lngthg, f_int *lnbr,   f_int *rnbr,
		   f_int *maxdb,  f_int *idbsiz, f_int *ngels,  f_int *nconts,
		   f_int *clist,  f_int *maxgel, char  *seq1,   char  *seq2,
		   char  *seq3,   char  *seq4,   char  *seq5,   char  *seqc2,
		   char  *seqg2,  char  *seqg3,  char  *seqc3,  char  *rnames,
		   f_int *maxseq, f_int *maxglm, f_int *sav1,   f_int *sav2,
		   f_int *sav3,   f_int *maxsav, f_int *cends,  f_int *nends,
		   f_int *maxcon, f_int *idev1,  char  *namarc, char  *nampro,
		   float *percd,  f_int *iokent, f_int *ishow,  f_int *minmat,
		   f_int *maxpg,  float *permax, f_int *irepsc, f_int *iopt,
		   f_int *nopt,   f_int *ansjok, f_int *ansfe,  f_int *iwing,
		   f_int *nbad,   char  *list,   f_int *iok,    f_int *minovr,
		   f_implicit seq1_l,	/* maxseq */
		   f_implicit seq2_l,	/* maxglm */
		   f_implicit seq3_l,	/* maxglm */
		   f_implicit seq4_l,	/* maxglm */
		   f_implicit seq5_l,	/* maxglm */
		   f_implicit seqc2_l,	/* maxglm,2 */
		   f_implicit seqg2_l,	/* maxglm,2 */
		   f_implicit seqg3_l,	/* maxglm */
		   f_implicit seqc3_l,	/* maxglm */
		   f_implicit rnames_l,	/* DB_NAMELEN */
		   f_implicit namarc_l,	/* F_NAMLEN */
		   f_implicit nampro_l,	/* DB_FILELEN */
		   f_implicit list_l   /* list_l */
		   );

f_proc_ret merge_(f_int *relpg,  f_int *lngthg, f_int *lnbr,   f_int *rnbr,
		  f_int *lincon, f_int *idbsiz);

f_proc_ret padcon_(f_int *relpg,  f_int *lngthg, f_int *lnbr,   f_int *rnbr,
		   f_int *ngels,  f_int *nconts, char  *gel,    f_int *lincon,
		   f_int *posn,   f_int *nc,     f_int *idbsiz, f_int *idevr,
		   f_int *maxgel,
		   f_implicit gel_l	/* maxgel */
		   );
		   

/* dbsr.f */
f_proc_ret baprep_(f_int *relpg,  f_int *lngthg, f_int *lnbr,   f_int *rnbr,
		   f_int *idbsiz, f_int *ngels,  f_int *nconts, f_int *cnum,
		   contig_list_t *clist,
		   f_int *idevn,  char  *nampro, char  *seq1,   f_int *maxseq,
		   f_int *maxgel, float *percd,  f_int *cends,  f_int *nends,
		   f_int *maxcon, f_int *sav1,   f_int *sav2,   f_int *sav3,
		   f_int *sav4,   f_int *sav5,   f_int *maxmat, f_int *idir,
		   f_int *minmat, f_int *nres,   f_int *iok,    f_int *mask,
		   f_implicit nampro_l,	/* DB_FILELEN */
		   f_implicit seq1_l	/* maxseq */
		   );

f_proc_ret inits_(void);

f_proc_ret shiftc_(f_int *relpg,  f_int *lngthg, f_int *lnbr,   f_int *rnbr,
		   f_int *ngels,  f_int *nconts, f_int *idevr,  f_int *idbsiz,
		   f_int *ign,    f_int *ncont,  f_int *dist);

f_proc_ret autoj_(f_int *relpg,  f_int *lngthg, f_int *lnbr,   f_int *rnbr,
		  f_int *idbsiz, f_int *ngels,  f_int *nconts, f_int *cnum,
		  contig_list_t *clist,
		  f_int *maxgel, f_int *iladd,  f_int *iradd,
		  char  *seq1,   f_int *maxseq, char  *seq2,   char  *seq3,
		  char  *seq4,   char  *seq5,   f_int *maxglm, f_int *sav1,
		  f_int *sav2,   f_int *sav3,   f_int *maxsav, f_int *cends,
		  f_int *nends,  f_int *maxcon, f_int *idevr,  char  *nampro,
		  float *percd,  f_int *idm,    char  *seqg3,  char  *seqc3,
		  f_int *iok,    f_int *lincon, f_int *llino,  f_int *mode,
		  f_int *mask,   f_int *lreg,   f_int *window, f_int *iwing,
		  f_int *bad,    f_int *minmat, f_int *maxp,   float *permax,
		  f_implicit seq1_l,	/* maxseq */
		  f_implicit seq2_l,	/* maxglm */
		  f_implicit seq3_l,	/* maxglm */
		  f_implicit seq4_l,	/* maxglm */
		  f_implicit seq5_l,	/* maxglm */
		  f_implicit nampro_l,	/* DB_FILELEN */		  
		  f_implicit seqg3_l,	/* maxglm */
		  f_implicit seqc3_l	/* maxglm */
		  );

f_int clinno_(f_int *lnbr,   f_int *idbzon, f_int *nconts, f_int *iin);

f_int chainl_(f_int *relpg,  f_int *lngthg, f_int *lnbr,   f_int *rnbr,
	      f_int *ngels,  f_int *nconts, f_int *idbsiz, f_int *iin);
		   
f_int nameno_(char  *name,   f_int *ngels,  f_int *idevn,
	      f_implicit name_l		/* DB_NAMELEN */
	      );


/* remgbc.f */
f_proc_ret remgbc_(f_int *relpg,  f_int *lngthg, f_int *lnbr,   f_int *rnbr,
		   f_int *ngels,  f_int *nconts, f_int *idbsiz,	char  *gel,
		   f_int *maxgel, f_int *idevr,  f_int *iok,    f_int *array,
		   f_int *iall,   f_int *iopt,
		   f_implicit gel_l	/* maxgel */
		   );

f_proc_ret breakc_(f_int *relpg,  f_int *lngthg, f_int *lnbr,   f_int *rnbr,
		   f_int *ngels,  f_int *nconts, f_int *idbsiz,	f_int *idevr,
		   f_int *ir,     f_int *iok);


/* initlu.f */
f_proc_ret initlu_(f_int *idm);


/* remcon.f */
f_proc_ret remcon_(f_int *relpg,  f_int *lngthg, f_int *lnbr,   f_int *rnbr,
		   f_int *ngels,  f_int *nconts, f_int *idbsiz,	f_int *icont,
		   char  *gel,    f_int *llino,  f_int *idevr,  f_int *maxgel,
		   f_int *temp,
		   f_implicit gel_l	/* maxgel */
		   );
		   
		   
/* fmtdb.f */
f_proc_ret fmtdb_(char  *seq1,   f_int *idim,    f_int *isw,    f_int *ise,
		  f_int *linlen, f_int *iden,
		  f_implicit seq_l     /* idim */
		  );

/*
 * ---------------------------------------------------------------------------
 * Prototypes for Fortran functions written in C that are (unfortunately)
 * also used by C. We should rewrite these to be true C routines with a stub
 * Fortran interface.
 * ---------------------------------------------------------------------------
 */


#endif /*_FORT_H_*/
