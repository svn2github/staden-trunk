/*
 * File: IO.c
 * Version:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: IO that doesn't belong in the C level gap IO library.
 *
 * Created: 23 February 1993
 * Updated: 12 July 1994. Split from IO.c
 *
 */

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

#include "IO.h"
#include "misc.h"
#include "edUtils.h"
#include "tagUtils.h"
#include "xalloc.h"
#include "gap-error.h"
#include "FtoC.h"


/*************************************************************
 *
 *************************************************************/


int get_vector_info(GapIO *io, int vector_id,
		    char *vector, int l_vector)
/*
 * Read vector information
 * Pointers for requested information can be NULL.
 */
{
    int err;
    GVectors v;

    if ( ! (vector&&l_vector>0) ) return 0;

    err = GT_Read(io,arr(GCardinal,io->vectors,vector_id-1),&v,sizeof(v),GT_Vectors);
    err = TextRead(io,v.name,vector,l_vector);

    return 0;
}




int get_clone_info(GapIO *io, int clone_id,
		   char *clone, int l_clone,
		   char *cvector, int l_cvector,
		   int *cvector_id)
/*
 * Read clone information
 * Pointers for requested information can be NULL.
 */
{
    int err;
    GClones c;

    if ( ! (clone&&l_clone>0 ||
	    cvector&&l_cvector>0 ||
	    cvector_id) ) return 0;

    err = GT_Read(io,arr(GCardinal,io->clones,clone_id-1),&c,sizeof(c),GT_Clones);

    if (clone && l_clone>0)
	err = TextRead(io,c.name,clone,l_clone);
	
    if (cvector_id)
	*cvector_id = c.vector;

    /* vector */
    get_vector_info(io, c.vector, cvector, l_cvector);
    

    return 0;
}





int get_subclone_info(GapIO *io,
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
		      int *scvector_id)
/*
 * Read subclone and clone information
 * Pointers for requested information can be NULL.
 */
{
    int err;
    GTemplates t;

    if (!subclone_id)
	return 0;

    if ( ! (clone&&l_clone>0 ||
	    cvector&&l_cvector>0 ||
	    subclone&&l_subclone>0 ||
	    scvector&&l_scvector>0 ||
	    insert_min ||
	    insert_max ||
	    strands ||
	    clone_id ||
	    cvector_id ||
	    scvector_id) ) return 0;

    /* template */
    err = GT_Read(io,arr(GCardinal,io->templates,subclone_id-1),&t,sizeof(t),GT_Templates);
    if (subclone && l_subclone>0)
	err = TextRead(io,t.name,subclone,l_subclone);

    if (insert_min)
	*insert_min = t.insert_length_min;
    if (insert_max)
	*insert_max = t.insert_length_max;
    if (strands)
	*strands = t.strands;
    if (scvector_id)
	*scvector_id = t.vector;
    if (clone_id)
	*clone_id = t.clone;


    /* vector */
    (void) get_vector_info(io,t.vector, scvector, l_scvector);


    /* clone */
    (void) get_clone_info(io, t.clone, clone, l_clone, cvector, l_cvector, cvector_id);

    return 0;
}






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
		  int *scvector_id)
/*
 * Read clone information
 * Pointers for requested information can be NULL.
 */
{
    int err;
    GReadings r;

    if ( ! (clone&&l_clone>0 ||
	    cvector&&l_cvector>0 ||
	    subclone&&l_subclone>0 ||
	    scvector&&l_scvector>0 ||
	    length ||
	    insert_min ||
	    insert_max ||
	    direction ||
	    strands ||
	    primer ||
	    clone_id ||
	    subclone_id ||
	    cvector_id ||
	    scvector_id) ) return 0;

    if (N>Nreadings(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid reading %d", N);
    }

    /* read record */
    err = gel_read(io, N, r);

    if (direction)
	*direction = STRAND(r);
    if (primer)
	*primer = PRIMER_TYPE(r);
    if (subclone_id)
	*subclone_id = r.template;
    if (length)
	*length = r.length;

    if (!r.template)
	return 0;

    /* template */
    (void) get_subclone_info(io,
			     r.template,
			     clone, l_clone,
			     cvector, l_cvector,
			     subclone, l_subclone, /* aka template */
			     scvector, l_scvector,
			     insert_min,
			     insert_max,
			     strands,
			     clone_id,
			     cvector_id,
			     scvector_id);


    return 0;
}
 

int io_get_extension(GapIO *io,
		     int N,	      /* gel number */
		     char *seq,	      /* buffer for sequence */
		     int max_seq,     /* size of buffer */
		     /* returns */
		     int *length,     /* length returned */
		     int *complement) /* reading is complemented */
/*
 * Get right cutoff for lowly Fortran Users
 * If a IGN[SC] tag exists we use the old method (ignore cutoff)
 * otherwise we look for the [SC]VEC tag and return only as much as is
 * available (ie not tagged).
 * We return 0 for success, or 1 for failure. Failure here is considered
 * finding an IGN tag (ie 'no cutoff available').
 */
{
    tagRecord rec;
    tag_id next;
    int err;
    GReadings r;
    char *wholeSeq;
    int t_st, t_end;

    if (N>Nreadings(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid reading %d", N);
    }

    /* read record */
    err = gel_read(io, N, r);

    wholeSeq = SeqReadStatic(io,r.sequence,r.length);


    /* set t_st/end to be _original_ (ie tag order) position of cutoff */
    t_st  = r.sense ? r.length - r.start : r.end - 1;
    t_end = r.length;
    
    /*
     * Look for ignore tag.
     * NB - we currently look through all the tags for this reading. Should
     * we stop as soon as we find the 3' cutoff tag? (ie can we assume there's
     * only ever one?)
     */
    next = first_tag(io, N);
    while (next) {
	(void) read_tag(io, next, &rec);
	
	if (strncmp(rec.type.c, "IGN", 3) == 0) {
	    *length = 0;
	    return 1; /* NO CUTOFF */
	}

	if (strncmp(&rec.type.c[1], "VEC", 3) == 0
	    && rec.position + rec.length >= t_st
	    && rec.position < t_end) {
	    t_end = rec.position - 1;
	}
	next = rec.next;
    }

    /* Prune down to length */
    *length = t_end - t_st;
    if (*length < 0) {
	*length = 0;
    } else if (*length > max_seq) {
	*length = max_seq;
	t_end = t_st + *length;
    }

    /* Switch around for complemented reads (no longer interested in t_end) */
    if (r.sense) {
	t_st  = r.length - t_end;
    }

    /* copy */
    *complement = r.sense;
    memcpy(seq, &wholeSeq[t_st], *length);
    seq[*length] = 0;

    return 0;
}






int io_mod_extension(GapIO *io,     /*  */
		     int N,	    /* gel number */
		     int shorten_by) /* ammount to shorten by */
/*
 * Modify right cutoff for lowly Fortran Users
 */
{
    int err;
    GReadings r;

    if (N>Nreadings(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid reading %d", N);
    }

    /* read record */
    err = gel_read(io, N, r);

    /*
     * if complemented,
     *    adjust "left" end
     * else
     *    adjust "right" end
     */
    if (r.sense) {
	/* complemented */
	r.start -= shorten_by;
	if (r.start<0) r.start = 0;
    } else {
	/* not complemented */
	r.end += shorten_by;
	if (r.end > r.length) r.end = r.length + 1;
    }

    /* write record back */
    err = gel_write(io, N, r);

    return 0;
}





/*************************************************************
 * The real stuff - the FORTRAN interface!!
 *************************************************************/


f_proc_ret readrn_(f_int *HANDLE, f_int *NGELS, f_int *NCONTS)
/*
 * Read ngels and ncontigs
 */
{
    GapIO *io;
    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    *NGELS  = NumReadings(io);
    *NCONTS = NumContigs(io);
    
    f_proc_return();
}

f_proc_ret writrn_(f_int *HANDLE, f_int *NGELS, f_int *NCONTS)
/*
 * Write ngels and ncontigs
 */
{
    GapIO *io;
    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    NumReadings(io) = *NGELS;
    NumContigs(io) = *NCONTS;

    /* write database record */
    DBDelayWrite(io);

    f_proc_return();
}




f_proc_ret readg_(f_int *HANDLE, f_int *N, f_int *RELPG, f_int *LNGTHG, f_int *LNBR, f_int *RNBR)
/*
 * Read a gel line
 */
{
    GapIO *io;
    int err;
    GReadings r;
    int l;
    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    if (*N>Nreadings(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid reading %d", *N);
    }

    /* read record */
    err = gel_read(io, *N, r);

    /* update */
    *LNBR = r.left;
    *RNBR = r.right;
    *RELPG = r.position;
    l = r.end - r.start - 1;
    *LNGTHG = (r.sense==GAP_SENSE_REVERSE)?-l:l;

    f_proc_return();

}

f_proc_ret writeg_(f_int *HANDLE, f_int *N, f_int *RELPG, f_int *LNGTHG, f_int *LNBR, f_int *RNBR)
/*
 * Write a gel line
 */
{
    GapIO *io;
    int err;
    GReadings r;
    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    if (*N>Nreadings(io)) err = io_init_reading(io,*N);

    /* read record */
    err = gel_read(io, *N, r);

    /* update */
    r.left = *LNBR;
    r.right = *RNBR;
    r.position = *RELPG;
    r.sequence_length = *LNGTHG>0?*LNGTHG:-*LNGTHG;
    r.sense = (*LNGTHG < 0) ? GAP_SENSE_REVERSE : GAP_SENSE_ORIGINAL;

    /* write record back */
    err = gel_write(io, *N, r);

    f_proc_return();

}

f_proc_ret readc_(f_int *HANDLE, f_int *N, f_int *LNGTHC, f_int *LGEL, f_int *RGEL)
/*
 * Read a contig line
 */
{
    GapIO *io;
    int err;
    GContigs r;
    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    if (*N>NumContigs(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid contig %d", *N);
    }

    /* read record */
    err = GT_Read(io,arr(GCardinal,io->contigs,*N-1),&r,sizeof(r),GT_Contigs);

    /* update */
    *LGEL = r.left;
    *RGEL = r.right;
    *LNGTHC = r.length;

    f_proc_return();
}



f_proc_ret writec_(f_int *HANDLE, f_int *N, f_int *LNGTHC, f_int *LGEL, f_int *RGEL)
/*
 * Write a contig line
 */
{
    GapIO *io;
    int err;
    GContigs r;
    int reset_anno = 0;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    if (*N>NumContigs(io)) {
	err = io_init_contig(io,*N);
	reset_anno = 1;
    }

    /* read record */
    err = GT_Read(io,arr(GCardinal,io->contigs,*N-1),&r,sizeof(r),GT_Contigs);

    /* update */
    r.left = *LGEL;
    r.right = *RGEL;
    r.length = *LNGTHC;

    if (reset_anno)
	r.annotations = 0;


    /* write record back */
    err = GT_Write(io,arr(GCardinal,io->contigs,*N-1),&r,sizeof(r),GT_Contigs);

    f_proc_return();
}





f_proc_ret readtg_(f_int *HANDLE, f_int *N, f_int *LPOS, f_int *LLEN, f_int *LCOM, f_int *LTYPE, f_int *NEXT, f_int *SENSE)
/*
 * Read a tag record
 */
{
    GapIO *io;
    int err;
    GAnnotations r;
    char *t;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    if (*N>Nannotations(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid annotation %d", *N);
    }

    /* read record */
    err = GT_Read(io,arr(GCardinal,io->annotations,*N-1),&r,sizeof(r),GT_Annotations);

    t = (char *)&r.type;
    *LPOS  = r.position;
    *LLEN  = r.length;
    *LCOM  = r.annotation;
    *LTYPE = (t[0]<<24) + (t[1]<<16) + (t[2]<<8) + (t[3]<<0);
    *NEXT  = r.next;
    *SENSE = r.strand;

    f_proc_return();
}


f_proc_ret writtg_(f_int *HANDLE, f_int *N, f_int *LPOS, f_int *LLEN, f_int *LCOM, f_int *LTYPE, f_int *NEXT, f_int *SENSE)
/*
 * Write a tag record
 */
{
    GapIO *io;
    int err;
    GAnnotations r;
    char *t = (char *)LTYPE;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    if (*N>Nannotations(io)) err = io_init_annotations(io,*N);

    /* no need to read record */

    /* update */
    r.position = *LPOS;
    r.length = *LLEN;
    r.annotation = *LCOM;
    r.type = (t[0]<<24) + (t[1]<<16) + (t[2]<<8) + (t[3]<<0);
    r.next = *NEXT;
    r.strand = *SENSE;

    /* write record back */
    err = GT_Write(io,arr(GCardinal,io->annotations,*N-1),&r,sizeof(r),GT_Annotations);

    f_proc_return();
}



f_proc_ret readcc_(f_int *HANDLE, f_int *N, f_int *ICNT, f_int *NEXT, char *NOTE, f_implicit NOTE_l)
/*
 * Read a comment record
 */
{
    GapIO *io;
    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();


    /*
     * TODO
     */
    verror(ERR_FATAL, "readcc", "not implemented");
    
    f_proc_return();
}


f_proc_ret writcc_(f_int *HANDLE, f_int *N, f_int *ICNT, f_int *NEXT, char *NOTE, f_implicit NOTE_l)
/*
 * Write a comment record
 */
{
    GapIO *io;
    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    /*
     * TODO
     */
    verror(ERR_FATAL, "writecc", "not implemented");
    
    f_proc_return();
}



f_proc_ret readn_(f_int *HANDLE, f_int *N, char *NAME, f_implicit NAME_l)
/*
 * Read a name
 */
{
    GapIO *io;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    if (*N>Nreadings(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid reading %d", *N);
    }

    Cstr2Fstr(io_rname(io, *N), NAME, NAME_l);

    f_proc_return();
}


/*
 * Write a name. Returns 0 for success, -1 for failure.
 */
int write_rname(GapIO *io, int reading, char *name) {
    GReadings r;
    int err = 0, len;

    if (reading > Nreadings(io))
	io_init_reading(io, reading);

    /* load */
    err |= gel_read(io, reading, r);

    /* allocate if needed */
    if (0 == r.name) {
	r.name = allocate(io, GT_Text);
	err |= gel_write(io, reading, r);
    }

    /* update */
    len = strlen(name);
    if (len > DB_NAMELEN+1)
	len = DB_NAMELEN+1;
    err |= TextWrite(io, r.name, name, len);
    io_wname(io, reading, name);

    return err ? -1 : 0;
}

f_proc_ret writen_(f_int *HANDLE, f_int *N, char *NAME, f_implicit NAME_l)
{
    GapIO *io;
    char c_name[DB_NAMELEN+1];

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    Fstr2Cstr(NAME, NAME_l, c_name, DB_NAMELEN+1);
    write_rname(io, *N, c_name);

    f_proc_return();
}




f_proc_ret readw_(f_int *HANDLE, f_int *N, char *GEL, f_int *MAXGEL, f_implicit GEL_l)
/*
 * Read a squence
 */
{
    GapIO *io;
    int err;
    GReadings r;
    char *seq;
    int l;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    if (*N>Nreadings(io)) {
	(void)gaperr_set(GAPERR_NOT_FOUND);
	GAP_ERROR_FATAL("invalid reading %d", *N);
    }

    /* read record */
    err = gel_read(io, *N, r);

    /* read sequence */
    seq = SeqReadStatic(io,r.sequence,r.length);

    l = r.end - r.start - 1;
    if (l > *MAXGEL) l = *MAXGEL;

    memcpy(GEL,seq+r.start,l);

    /*
     * What was this for? Doesn't appear to be required.
     *
     *  {
     *      int i;
     *      for(i=l;i<*MAXGEL;i++) GEL[i]='?';
     *  }
     */
    f_proc_return();
}





int readrd(f_int *HANDLE, int gel, char *type, char *trace,
	   int type_l, int trace_l) {
    GapIO *io;

    if ( (io = io_handle(HANDLE)) == NULL) return 0;

    return io_read_rd(io, gel, trace, trace_l, type, type_l);
}



/*
 * swap readings N and M
 */
int swap_read(GapIO *io, int M, int N) {
    int err = 0;
    GCardinal temp;
    GReadings r;
    char name[DB_NAMELEN+1];

    if (N>Nreadings(io)) err |= io_init_reading(io,N);
    if (M>Nreadings(io)) err |= io_init_reading(io,M);
    if (err)
	GAP_ERROR_FATAL("io_init_reading (swap %d %d)", N, M);

#if GAP_CACHE!=0
    /* swap name elements */
    strcpy(name, io_rname(io, N));
    io_wname(io, N, io_rname(io, M));
    io_wname(io, M, name);
#endif

    /* swap record array elements */
    temp = arr(GCardinal,io->readings,N-1);
    arr(GCardinal,io->readings,N-1) = arr(GCardinal,io->readings,M-1);
    arr(GCardinal,io->readings,M-1) = temp;

#if GAP_CACHE!=0
    /* swap structure array elements */
    r = arr(GReadings, io->reading, N-1);
    arr(GReadings, io->reading, N-1) = arr(GReadings, io->reading, M-1);
    arr(GReadings, io->reading, M-1) = r;
#endif

    /* write array */
    err = ArrayDelay(io, io->db.readings, io->db.Nreadings, io->readings);
    return err;
}


f_proc_ret swapnm_(f_int *HANDLE, f_int *N, f_int *M)
{
    GapIO *io;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();
    swap_read(io, *M, *N);
    f_proc_return();
}

f_proc_ret movnm_(f_int *HANDLE, f_int *N, f_int *M)
/*
 * Part swap, part move N to M.
 *
 * This is mainly a move from N to M, but it still has to swap the entries
 * in the stored (and memory) reading record number array. A straight copy
 * into this array does not work. Why?
 */
{
    GapIO *io;
    int err = 0;
    GCardinal temp;
    int nnote;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    if (*N>Nreadings(io)) err |= io_init_reading(io,*N);
    if (*M>Nreadings(io)) err |= io_init_reading(io,*M);
    if (err)
	GAP_ERROR_FATAL("io_init_reading (swap %d %d)", *N, *M);

#if GAP_CACHE!=0
    /* Copy name elements */
    io_wname(io, *M, io_rname(io, *N));
#endif

    /* Swap reading record numbers in io->db.readings stored array */
    temp = arr(GCardinal,io->readings,*N-1);
    arr(GCardinal,io->readings,*N-1) = arr(GCardinal,io->readings,*M-1);
    arr(GCardinal,io->readings,*M-1) = temp;


#if GAP_CACHE!=0
    /* Copy cached GReadings structure */
    arr(GReadings, io->reading, *M-1) = arr(GReadings, io->reading, *N-1);
#endif

    /* Notes have two-way linked lists, so we update the link back */
    nnote = arr(GReadings, io->reading, *M-1).notes;
    if (nnote) {
	GNotes n;
	note_read(io, nnote, n);
	n.prev = *M;
	note_write(io, nnote, n);
    }

    /* write array */
    err = ArrayDelay(io, io->db.readings, io->db.Nreadings, io->readings);

    f_proc_return();
}
