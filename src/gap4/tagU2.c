/*
 * File: tagU2.c
 * Version:
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: Tag IO routines required by non-X assembly program
 *
 * Created:
 * Updated:
 *
 *
 * 8-Jul-92
 *      getext_() searchs for an IGN tag, to determine if cutoff should be
 *      ignored
 * 7-Aug-92
 *	now initial tags can be specified in the sequence file
 *      format is ";;%4s %6d %6d %s\n",type,position,length,comment
 *      The comment is optional
 * 27-Aug-92
 *      modext() modifies cutoff data
 * 20-Apr-94
 *      Tag format has now changed. Also masses of other unlisted things!
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "tagUtils.h"
#include "misc.h"
#include "contigEditor.h"
#include "fort.h"
#include "edUtils.h"
#include "IO.h"
#include "io_utils.h"
#include "gap-dbstruct.h"
#include "expFileIO.h"
#include "text_output.h"
#include "dbcheck.h"
#include "active_tags.h"
#include "FtoC.h"
#include "qualIO.h"
#include "gap_globals.h"


/*
 * Constants for reading in text files
 */
#define fn_len 256		/* maximum file name length */
#define l_line 100		/* maximum line length */


#define COMMENT_HEAD_ID (1)


/*************************************************************
 *
 *************************************************************/

int read_tag(GapIO *io, tag_id n, tagRecord *t)
{
    readtg_(handle_io(io), &n, &t->position, &t->length, &t->comment,
	    &t->type.i, &t->next, &t->sense);
    return 0;
}

int write_tag(GapIO *io, tag_id n, tagRecord t)
{
    writtg_(handle_io(io), &n, &t.position, &t.length, &t.comment,
	    &t.type. i, &t.next, &t.sense);
    return 0;
}

/*************************************************************
 *
 *************************************************************/

/* comment interface */
tag_id get_free_tag(GapIO *io)
/*
 *
 */
{
    tag_id head;
    tagRecord freerec;
    tag_id free_id;

    (void) io_read_free_annotation(io,&head);
    if (head != 0) {
	/*
	 * if a free slot somewhere, use it
	 */
	free_id = head;
	(void) read_tag(io, free_id,&freerec);
	head = freerec.next;
	(void) io_write_free_annotation(io,&head);
    } else {
	/*
	 * extend file
	 */
	free_id = Nannotations(io)+1;
	io_init_annotations(io,free_id);
    }

    return free_id;
}



void blank_tag_rec(GapIO *io, tag_id t)
/*
 * Blank out fields in tag record t
 */
{
    tagRecord r;
    
    (void) read_tag(io, t, &r);
    
    r.position = 0;
    r.length = 0;
    r.type.i = 0x20202020;
    r.comment = 0;
    r.next = 0;
    r.sense = 0;
    
    (void) write_tag(io, t,r);
}




void delete_tag_rec(GapIO *io, tag_id t)
/*
 * remove t from file, discarding comment if necessary
 */
{
    tag_id head;
    tagRecord freerec;

    /*
     * reclaim comment
     */
    freerec.comment = 0; /* just in case read_tag fails */
    (void) read_tag(io, t,&freerec);
    if (freerec.comment) {
	deallocate(io, freerec.comment);
	freerec.comment = 0;
    }
    
    (void) io_read_free_annotation(io,&head);
    freerec.next = head;
    (void) write_tag(io, t,freerec);
    head = t;
    (void) io_write_free_annotation(io,&head);
}

/*
 * Remove the annotation pointed to by 'last'.next
 * last_type is 0 for contig, 1 for reading or 2 for annotation.
 * Returns the new next annotation for 'last'.
 */
int delete_tag(GapIO *io, int last, int curr, int last_type) {
    GAnnotations a;

    tag_read(io, curr, a);
    delete_tag_rec(io, curr);

    switch(last_type) {
    case 0: {/* contig */
	GContigs c;
	
	contig_read(io, last, c);
	c.annotations = a.next;
	contig_write(io, last, c);
	break;
    }

    case 1: {/* reading */
	GReadings r;
	
	gel_read(io, last, r);
	r.annotations = a.next;
	gel_write(io, last, r);
	break;
    }

    case 2: {/* cannotation */
	GAnnotations t;
	
	tag_read(io, last, t);
	t.next = a.next;
	tag_write(io, last, t);
	break;
    }
    }

    return a.next;
}

comment_id put_comment(GapIO *io, char *c)
/*
 * Allocate a new record, write comment, return record number.
 */
{
    int r;

    r = allocate(io,GT_Text);

    (void) TextWrite(io,r,c,strlen(c)+1);

    return r;
}





tag_id first_tag(GapIO *io, int N)
{
    int anno;

    io_read_annotation(io,N,&anno);
    return anno;
}


void update_tag(GapIO *io, int N, tag_id anno)
{
    io_write_annotation(io,N,&anno);
}



void insert_new_tag2(GapIO *io, int into,
		     int *cache, int cache_len, int *cache_pos,
		     int pos, int length, char *type, char *comment, int sense)
{
    tag_id prev, next, newt;
    tagRecord prev_tag, next_tag, new_tag;
    int i;

    /* Initialise cache if required */
    if (pos >= *cache_pos) {
	next = cache[*cache_pos];
	if (!next)
	    next = first_tag(io, into);
	while (next) {
	    read_tag(io, next, &next_tag);
	    if (next_tag.position < *cache_pos) {
		next = next_tag.next;
		continue;
	    }
	    if (next_tag.position > pos)
		break;
	    for (i = *cache_pos; i < next_tag.position; i++)
		cache[i] = cache[*cache_pos];
	    for (; i <= next_tag.position; i++) {
		cache[i] = next;
	    }
	    *cache_pos = i-1;
	    next = next_tag.next;
	}
	for (i = *cache_pos+1; i <= pos; i++)
	    cache[i] = cache[*cache_pos];
	*cache_pos = pos;
    }

    /* Find previous and next tags - quick lookup in cache */
    prev = cache[pos];
    if (!prev) {
	next = first_tag(io, into);
    } else {
	read_tag(io, prev, &prev_tag);
	next = prev_tag.next;
    }

    /* Create and initialise new tag */
    newt = get_free_tag(io);
    new_tag.position = pos;
    new_tag.length = length;
    strncpy(new_tag.type.c,type,4);
    if (comment!=NULL)
	new_tag.comment = put_comment(io, comment);
    else
	new_tag.comment = 0;
    new_tag.next = next;
    new_tag.sense = sense;
    write_tag(io, newt, new_tag);

    /* Update cache */
    i = pos;
    while (i <= *cache_pos && cache[i] == prev)
	cache[i++] = newt;
    if (pos > *cache_pos) {
	int j, k = cache[*cache_pos];
	for (j = *cache_pos; j < pos; j++)
	    cache[j] = k;
	cache[pos] = newt;
	*cache_pos = pos;
    }
    
    /* Link previous tag */
    if (prev) {
	prev_tag.next = newt;
	write_tag(io, prev, prev_tag);
    } else {
	update_tag(io, into, newt);
    }
}


void insert_NEW_tag(GapIO *io, int N, int pos, int length, char *type,
		    char *comment, int sense)
/*
 * :-)
 */
{
    tag_id last,next,new;
    tagRecord last_tag,next_tag,new_tag;

    last = 0;			/* indicated start of linked list */
    next = first_tag(io, N);

    /*
     * Find position to insert
     */
    while (next) {
	read_tag(io, next,&next_tag);
	if (next_tag.position > pos) break;
	last_tag = next_tag;
	last = next;
	next = next_tag.next;
    }
    

    /*
     * create and initialise new tag
     */
    new = get_free_tag(io);

    new_tag.position = pos;
    new_tag.length = length;
    strncpy(new_tag.type.c,type,4);
    if (comment!=NULL)
	new_tag.comment = put_comment(io, comment);
    else
	new_tag.comment = 0;
    new_tag.next = next;
    new_tag.sense = sense;

    write_tag(io, new,new_tag);

    /*
     * Update things in the chain that require it
     */
    if (last) {
	last_tag.next = new;
	write_tag(io, last,last_tag);
    } else {
	/* update record */
	update_tag(io, N,new);
    }

}


/*
 * Create a tag for a gel
 */
void create_tag_for_gel(GapIO *io, int gel, int gellen, char *tag,
			int *cache, int cache_len, int *cache_pos,
			int unpadded_tags) {
    char type[5], *comment;
    int start, end, strand, npos;

    /*
     * 'tag' includes a few ascii bits and bobs (start, end, etc) plus
     * the comment so we'll simply alloc comment to be the same length.
     */
    if (NULL == (comment = (char *)xmalloc(strlen(tag)))) {
	return;
    }

    if (-1 == tag2values(tag, type, &start, &end, &strand, comment)) {
	verror(ERR_WARN, "create_tag_for_gel",
	       "Failed to parse tag \"%s\".", tag);
	return;
    }


    /*
     * Tag is in depadded sequence. Convert to a padded coordinate.
     */
    if (unpadded_tags && gel > 0) {
	GReadings r;

	gel_read(io, gel, r);
	if (r.sequence) {
	    char *seq;
	    int new_start = start, new_end = end;
	    int i, pads = 0;

	    int st  = r.sense == 0 ? 1          : r.length;
	    int en  = r.sense == 0 ? r.length+1 : 0;
	    int dir = r.sense == 0 ? 1          : -1;

	    seq = TextAllocRead(io, r.sequence);
	    for (i = st; i != en; i+=dir) {
		int j = r.sense ? r.length+1-i : i;
		if (seq[i-1] == '*') {
		    pads++;
		    continue;
		}
		if (j-pads == start) {
		    new_start = start + pads;
		}
		if (j-pads == end) {
		    new_end = end + pads;
		}
	    }

	    start = new_start;
	    end = new_end;
	    xfree(seq);
	}
    } else if (unpadded_tags) {
	/* Unpadded position, but in the contig */
	char *cons;
	int i, pads;
	int clen = io_clength(io, -gel);
	int new_start = start, new_end = end;

	if (!(cons = (char *)xmalloc(clen+1)))
	    return;
	calc_consensus(-gel, 1, clen, CON_SUM, cons, NULL, NULL, NULL,
		       consensus_cutoff, quality_cutoff,
		       database_info, (void *)io);

	for (pads = 0, i = 1; i <= clen; i++) {
	    if (cons[i-1] == '*') {
		pads++;
		continue;
	    }
	    if (i-pads == start)
		new_start = start + pads;
	    if (i-pads == end)
		new_end = end + pads;
	}
	start = new_start;
	end = new_end;

	xfree(cons);
    }


#if 0
    /*
     * This is broken anyway (should be r.length), and nothing seems to use
     * it, so we now make sure that create_tag_for_gel creates a tag
     * at the specified locations regardless of orientation.
     */
    if (gellen<0) {
        /* gel is not in original sense */
	npos = abs(gellen) - end + 1;
    } else {
	npos = start;
    }
#endif
    npos = start;

    /* sanity checks */
    if (start < 1 || end > abs(gellen))
	verror(ERR_WARN, "create_tag_for_gel",
	      "Tag %s overlaps gel reading (#%d) ends (1..%d) - not entered",
	       tag, gel, abs(gellen));
    else if (start > end)
	verror(ERR_WARN, "create_tag_for_gel",
	       "Tag %s has negative length, for gel %d!", tag, gel);
    else {
	if (cache)
	    insert_new_tag2(io, (tag_id)gel,
			    cache, cache_len, cache_pos,
			    npos, end-start+1, type, comment, strand);
	else
	    insert_NEW_tag(io, (tag_id)gel, npos, end-start+1, type, comment,
			   strand);
    }

    xfree(comment);
}


/*
 * >>>>>>>>>>>>>>>>> This routine is no longer used
 */
#ifdef NOTUSED

void movtag_ (int_f *from, int_f *to, int_f *handle)
/*
 * Move tag information of gel ``from'' to gel ``to'',
 * and perform garbage collection on old gel ``to''
 */

/*
 * YUK! obselete?
 * This is ok but it would be better to have a movgel_() routine
 * written in C
 */
 
{
    tagRecord freerec;
    tag_id this,next;
    GapIO *io;
    
    if ( (io = io_handle(handle)) == NULL) return;

    /* Throw away ``to'' tag records */
    (void) read_tag(io, (tag_id) *to,&freerec);
    
    if (freerec.comment) {
	deallocate(io, freerec.comment);
	freerec.comment = 0;
    }
    
    next = freerec.next;
    while ( next ) {
        this = next;
	(void) read_tag(io,  this , &freerec );
        next = freerec.next;
        delete_tag_rec (io, this);
    }
    
    /* copy ``from'' record to ``to'' record */
    (void) read_tag(io, (tag_id) *from,&freerec);
    (void) write_tag(io, (tag_id) *to,freerec);
    
    /* initialise the hole to blank */
    blank_tag_rec(io, *from);
    
}

#endif


char *get_comment(GapIO *io, comment_id cp)
{
    extern char *TextAllocRead(GapIO *io, int rec);

    if (! cp) return NULL;

    return TextAllocRead(io,cp);

}







/*************************************************************
 *
 *************************************************************/



void tagfil_(f_int *relpg, f_int *lngthg, f_int *lnbr, f_int *rnbr,
	     f_int *ngels, f_int *nconts, f_int *idbsiz,
	     char *namarc, f_int *idev, f_int *verb, f_int *unpadded_tags,
	     f_implicit l_namarc)
/*
 * Open up file namarc
 * For each reading ("ID" lines) create tags ("TG/TC" lines)
 */
{
    GapIO *io;
    char fn[fn_len];
    char line[l_line];
    int gel, clen;
    char *cp;
    FILE *fp;
    char *comment = NULL;
    int comment_alloc = 1024, comment_size = 0;

    if ( (io = io_handle(idev)) == NULL) return;

    if (NULL == (comment = (char *)xmalloc(comment_alloc+1)))
	return;

    Fstr2Cstr(namarc,l_namarc,fn,fn_len);

    if (NULL == (fp = fopen(fn, "r")))
	return;

    gel = 0; /* invalid gel no. */
    if (fgets(line, l_line, fp) == NULL) {
	verror(ERR_WARN, "tagfil",
              "Invalid gel no. found when entering tags.");
	return;
    }

    do {
	cp = NULL;

	if (strncmp(line, "ID", 2) == 0) {
	    int i;

	    for (i = 5; line[i]; i++) {
		if (line[i] == '\n' || line[i] == ' ') {
		    line[i] = '\0';
		    break;
		}
	    }

	    /* new gel number */
	    gel = get_gel_num(io, &line[5], GGN_NAME);
	    if (gel == -1)
		verror(ERR_WARN, "tagfil", "Unknown gel name '%s'.", &line[5]);

	} else if (strncmp(line, "TG", 2) == 0 ||
		   strncmp(line, "TC", 2) == 0) {
	    /* tag */

	    int gel_num; /* or contig number! */
	    int len;	 /* gel or contig length */

	    if (gel > 0) {
		if (line[1] == 'C') {
		    gel_num = chainl_(relpg, lngthg, lnbr, rnbr,
				      ngels, nconts, idbsiz, &gel);
		    gel_num = clinno_(lnbr, idbsiz, nconts, &gel_num);
		    len = relpg[gel_num - 1];
		    gel_num = gel_num - *idbsiz;
		} else {
		    GReadings r;

		    gel_num = gel;
		    gel_read(io, gel, r);

		    len = r.length;
		}
	    } else
		gel_num = 0;

	    clen = strlen(&line[5]);
	    if (clen >= comment_alloc) {
		char *ctmp;

		comment_alloc = clen*2 + 1;
		if (NULL == (ctmp = (char *)xrealloc(comment, comment_alloc))){
		    xfree(comment);
		    return;
		}
		comment = ctmp;
	    }
	    strcpy(comment, &line[5]);
	    comment_size = clen;

	    do {
		cp = fgets(line, l_line, fp);
		if (cp && strncmp(&line[2], "        ", 8) == 0) {
		    clen = strlen(&line[10]);
		    if (clen + comment_size >= comment_alloc) {
			char *ctmp;

			comment_alloc = (comment_size + clen)*2 + 1;
			if (NULL == (ctmp = (char *)xrealloc(comment,
							     comment_alloc))){
			    xfree(comment);
			    return;
			}
			comment = ctmp;
		    }
		    strcat(comment, &line[10]);
		    comment_size += clen;
		} else {
		    break;
		}
	    } while (cp);

	    if (gel_num) {
		create_tag_for_gel(io, gel_num, len, comment, NULL, 0, NULL,
				   *unpadded_tags);
		UpdateTextOutput();
	    }
	}

	if (!cp)
	    cp = fgets(line, l_line, fp);

    } while (cp);

    if (comment)
	xfree(comment);

    fclose(fp);
}





/*************************************************************
 * FORTRAN interface for a routine to search on gel type
 ************************************************************/


/*
 * Get all tags of a given type for a gel.
 *
 * First time, call ctagget() with a non-zero gel number and a type.
 * Returns a pointer to the annotation or NULL for no more.
 * Returns -1 for error.
 *
 * To read subsequent tags, call with a gel number of zero.
 */
GAnnotations *ctagget(GapIO *io, int gel, char *type) {
    static GAnnotations a;
    static int anno;
    GCardinal itype = str2type(type);

    if (!gel)
	anno = a.next;
    else
	if (-1 == io_read_annotation(io, gel, &anno))
	    return (GAnnotations *)-1;

    while (anno) {
	GT_Read(io, arr(GCardinal, io->annotations, anno-1),
		&a, sizeof(a), GT_Annotations);

	if (itype == a.type)
	    return &a;

	anno = a.next;
    }

    return (GAnnotations *)NULL;
}

GAnnotations *vtagget(GapIO *io, int gel, int num_t, char **type) {
    static GAnnotations a;
    static int anno;
    int arg;

    if (!gel)
	anno = a.next;
    else
	if (-1 == io_read_annotation(io, gel, &anno))
	    return (GAnnotations *)-1;

    while (anno) {
	tag_read(io, anno, a);
	for (arg = 0; arg < num_t; arg++) {
	    if (str2type(type[arg]) == a.type) {
		return &a;
	    }
	    /*
	     * This was the code at 21st May 2001, but it fails to work
	     * correctly when there are tags in the database that are not
	     * in our current GTAGDB.
	     *
	     * if (idToIndex(type[arg]) == idToIndex(type2str(a.type,str))) {
	     *   return &a;
	     * }
	     */
	}
	anno = a.next;
    }

    return (GAnnotations *)NULL;
}

f_proc_ret tagget_(/* input */
		   int_f *GEL,	/* or 0 for next in gel */
		   char *TYPE,
		   /* returns */
		   int_f *POS,
		   int_f *LEN,
		   int_f *IDEV,
		   int_f *START,
		   /* implicits */
		   int_fl TYPE_l
		   )
/*
 * Get position and length ionformation of all tags of a given type for a gel.
 *
 * First time, call tagget_() with a non-zero gel number and a type.
 * Returns position and length of first tag.
 * If there is no other tag, returns length as -1.
 * To read subsequent tags, call with a gel number of zero. Type is ignored.
 */
{
    /* NOT REENTRANT */
    static char type[5];
    static tagRecord t;
    static int used_start;

    int looking;
    tag_id next_id;
    GapIO *io;

    if ( (io = io_handle(IDEV)) == NULL) f_proc_return();

    if (*GEL) {
	GReadings r;

	Fstr2Cstr(TYPE,TYPE_l,type,sizeof(type));
	next_id = first_tag(io, *GEL);

	gel_read(io, *GEL, r);
	used_start = r.start;

    } else {
	next_id = t.next;
    }

    looking = 1; *POS = -1; *LEN = -1;
    while(next_id && looking) {
	read_tag(io, next_id,&t);

	/*
	 * check to see if we have the correct type
	 */
	if (strncmp(t.type.c,type,4)==0) {
	    looking = 0;
	    *POS = t.position;
	    *LEN = t.length;
	    *START = used_start;
	}

	next_id = t.next;
    }

    f_proc_return();
}



/*************************************************************
 * Move tags left or right by one
 *************************************************************/

void tag_shift_for_insert(GapIO *io, int seq, int pos)
/*
 * Shifts the tags as if we'd inserted a base at position pos.
 */
{
    tagRecord t;
    tag_id next;

    next = first_tag(io, seq);
    while (next) {
	read_tag(io, next,&t);
	/*
	 * Move tags accordingly
	 */
	if (t.position >= pos) {
	    t.position++;
	    write_tag(io, next,t);
	} else if (t.position + t.length > pos) {
	    t.length++;
	    write_tag(io, next,t);
	}
	
	next = t.next;
    }
    
}



void tag_shift_for_delete(GapIO *io, int seq, int pos)
/*
 * Shifts the tags as if we'd deleted a base at position pos.
 */
{
    tagRecord t;
    tag_id next;

    next = first_tag(io, seq);
    while (next) {
	read_tag(io, next,&t);
	/*
	 * Move tags accordingly
	 */
	if (t.position >= pos) {
	    t.position--;
	    write_tag(io, next,t);
	} else if (t.position + t.length > pos) {
	    t.length--;
	    write_tag(io, next,t);
	}
	
	next = t.next;
    }

}



/*************************************************************
 * C interface for routines to handle tags in readings and consensus.
 *************************************************************/

/*
 * Moves tags within a contig right by 'NC'. These tags have just been placed
 * at position 'POSN' by PADCON() or UNLNKR()
 */
void shift_contig_tags(GapIO *io, int contig, int posn, int dist) {
    GContigs c;
    GAnnotations a;
    int anno;

    contig_read(io, contig, c);
    anno = c.annotations;

    while (anno) {
	tag_read(io, anno, a);
	
	if (a.position >= posn) {
	    /* Annotation is after; adjust position */
	    a.position += dist;
	    tag_write(io, anno, a);
	} else if (a.position + a.length > posn) {
	    /* Annotation spans insert - adjust length */
	    a.length += dist;
	    tag_write(io, anno, a);
	}
	/* Else annotation is before; do nothing */

	anno = a.next;
    }

    return;
}

/*
 * Merge tag lists in contigs. We merge tags from CONT2 into CONT1.
 * OFFSET is the relative positions between CONT1 and CONT2.
 */
void merge_contig_tags(GapIO *io, int contig1, int contig2, int off) {
    GContigs c1, c2;
    GAnnotations a1, a2, p;
    int a1n, a2n, pn = 0;

    /* load up contig records */
    contig_read(io, contig1, c1);
    contig_read(io, contig2, c2);

    /*
     * Merge is a simple merge sort - scan through both lists building up
     * a single new list.
     */

    /* special cases */
    a1n = c1.annotations;
    a2n = c2.annotations;

    /* no list for second contig - do nothing (other than copy 1 to 3) */
    if (!a2n) {
	return;
    }

    /* no list for first contig - set it to the second list & shift */
    if (!a1n) {
	c1.annotations = a2n;
	contig_write(io, contig1, c1);
	c2.annotations = 0;
	contig_write(io, contig2, c2);
	shift_contig_tags(io, contig1, 1, off);
	return;
    }

    /*
     * initialise - read the first two annotations and set GContigs.annotation
     * to start with the lowest positioned.
     */
    tag_read(io, a1n, a1);
    tag_read(io, a2n, a2);

    if (a1.position > a2.position + off) {
	c1.annotations = a2n; contig_write(io, contig1, c1);
    }

    /*
     * Now loop for the rest reading from list a1 or a2. We simply link the
     * last annotation to the lowest of the next two, and then step that onto
     * the next in its list.
     */
    for (;;) {
	if (a1.position <= a2.position + off) {
	    /*
	     * link prev = a1
	     * set  prev = a1
	     * set  a1   = a1.next
	     * check end of a1 list
	     */
	    if (pn) {
		if (p.next != a1n) { /* is p on a2 list? */
		    p.next  = a1n;
		    p.position += off;
		    tag_write(io, pn, p);
		}
	    }
	    
	    pn = a1n; tag_read(io, pn, p);
	    
	    if (a1n = a1.next) {
		tag_read(io, a1n, a1);
	    } else {

		/*
		 * link prev = a2
		 * shift a2 tags
		 * end
		 */
		p.next = a2n;
		tag_write(io, pn, p);
		pn = a2n;
		if (off) {
		    while (pn) {
			tag_read(io, pn, p);
			p.position += off;
			tag_write(io, pn, p);
			pn = p.next;
		    }
		}
		break;
	    }
	} else {
	    /*
	     * link prev = a2
	     * set  prev = a2
	     * set  a2   = a2.next
	     * check end of a2 list
	     */
	    if (pn) {
		if (p.next != a2n) {
		    p.next  = a2n;
		} else { /* is p on a2 list? */
		    p.position += off;
		}
		tag_write(io, pn, p);
	    }
	    
	    pn = a2n; tag_read(io, pn, p);
	    
	    if (a2n = a2.next) {
		tag_read(io, a2n, a2);
	    } else {

		/*
		 * link prev = a1
		 * end
		 */
		p.next = a1n;
		p.position += off;
		tag_write(io, pn, p);
		break;
	    }
	}
    }

    /* clear anno list for c2 - needed? */
    c2.annotations = 0;
    contig_write(io, contig2, c2);

    return;
}


typedef struct {
    int n;
    GAnnotations a;
} comp_anno_t;

static int anno_sort_func(const void *v1, const void *v2) {
    comp_anno_t *a1 = (comp_anno_t *)v1;
    comp_anno_t *a2 = (comp_anno_t *)v2;

    /* Ignore any overflow cases - our tags aren't THAT far out! */
    return a1->a.position - a2->a.position;
}

/*
 * Complement the tags in a contig.
 *
 * This simply involves negating the position of each one and reversing the
 * order of the list.
 */
void complement_contig_tags(GapIO *io, int contig) {
    int anno;
    int e, len = io_clength(io, contig);
    int max_anno = 100;
    int num_anno = 0;
    comp_anno_t *annos = NULL;
    int i, i_end;

    if (e = io_read_annotation(io, -contig, &anno)) {
	verror(ERR_FATAL, "complement_contig_tags", "tag read error %d", e);
    }

    if (!anno)
	return;

    if (NULL == (annos = (comp_anno_t *)xmalloc(max_anno * sizeof(*annos))))
	return;

    /*
     * Read in anno list to memory
     */
    while (anno) {
	annos[num_anno].n = anno;
	tag_read(io, anno, annos[num_anno].a);
	anno = annos[num_anno].a.next;
	if (++num_anno == max_anno) {
	    max_anno *= 2;
	    if (NULL == (annos = xrealloc(annos, max_anno * sizeof(*annos))))
		return;
	}
    }
    
    /* Negate positions, strand */
    for (i = 0; i < num_anno; i++) {
	if (annos[i].a.strand != 2)
	    annos[i].a.strand = 1-annos[i].a.strand;
	annos[i].a.position = len+1 -
	    (annos[i].a.position + annos[i].a.length - 1);
    }

    /* Reverse list as a first guess at sorting */
    i_end = num_anno/2;
    for (i = 0; i < i_end; i++) {
	comp_anno_t a;
	a = annos[i];
	annos[i] = annos[num_anno-1-i];
	annos[num_anno-1-i] = a;
    }

    /* sort annotations - just use qsort. */
    qsort(annos, num_anno, sizeof(*annos), anno_sort_func);

    /* Write back annotation data */
    for (i = 0; i < num_anno; i++) {
	annos[i].a.next = (i+1 == num_anno) ? 0 : annos[i+1].n;
	tag_write(io, annos[i].n, annos[i].a);
    }
    
    /* reset contig annotations pointer */
    io_write_annotation(io, -contig, &annos[0].n);
    
    xfree(annos);

    return;
}


/*
 * Called by "break contig". When we're splitting a contig in half we need
 * to move the annotations too.
 * For annotations that overlap the two contigs we duplicate and adjust the
 * lengths.
 *
 * The new contigs used to overlap between POSL and POSR.
 * Ie POSL is the position of the left end of cont2
 *    POSR is the position of the right end of cont1
 * Hence if cont1 and cont2 do not overlap POSL > POSR (eg from disassembly)
 *
 * The logic here is quite tricky as we need to take into account tags that
 * are entirely in left contig, right contig, both contigs, or neither contig,
 * plus all the partial combinations requiring clipping. The implementation
 * is to firstly duplicate if a tag is shared, and then clip left and/or
 * right ends as required.
 *
 * HINT: If you want to change this, draw out all the cases of tag positions
 * on paper, with base numbers assigned, to work out the exact numbering
 * (removing out-by-one errors) for both +ve and -ve overlap sizes.
 */
void split_contig_tags(GapIO *io, int cont1, int cont2, int posl, int posr) {
    GContigs c1, c2;
    GAnnotations l, r, lastr;
    int ln, rn, lastrn = 0, end = 0;

    contig_read(io, cont1, c1);
    contig_read(io, cont2, c2);

    /* easy case */
    if (c1.annotations == 0)
	return;

    /*
     * scan annotation list in c1.
     * (assume c2 has blank list - it's just been created we hope).
     */
    for (ln = c1.annotations; ln; ln = l.next) {
	/* read tag */
	tag_read(io, ln, l);

	/* is tag in both contigs? If so, duplicate it */
	if (l.position <= posr && l.position + l.length-1 >= posl) {
	    /* Get a new tag */
	    rn = get_free_tag(io);
	    tag_read(io, rn, r);

	    /* Copy the annotation */
	    if (l.annotation) {
		char *tbuf = TextAllocRead(io, l.annotation);
		if (tbuf) {
		    r.annotation = allocate(io, GT_Text);
		    TextWrite(io, r.annotation, tbuf, strlen(tbuf));
		    xfree(tbuf);
		}
	    } else {
		r.annotation = 0;
	    }

	    /* set position and length of new tag */
	    r.type = l.type;
	    r.strand = l.strand;
	    r.position = l.position - posl + 1;
	    if (r.position < 1)
		r.position = 1;
	    r.length = l.position + l.length - posl - (r.position-1);
	    r.next = 0;

	    /* save new tag */
	    tag_write(io, rn, r);
	    if (!lastrn) {
		c2.annotations = rn;
		contig_write(io, cont2, c2);
	    } else {
		lastr.next = rn;
		tag_write(io, lastrn, lastr);
	    }
	    lastrn = rn;
	    memcpy(&lastr, &r, sizeof(lastr));

	    if (l.position + l.length > posr) {
		/* truncate length of existing tag */
		l.length = posr - l.position + 1;
		tag_write(io, ln, l);
	    }

	    end = ln;
	}

	/* Is tag in neither contig */
	else if (l.position > posr && l.position + l.length-1 < posl) {
	    delete_tag_rec(io, ln);
	    continue;
	}

	/* Only contig2 => shift it down; may need to clip */
	else if (l.position > posr && l.position + l.length-1 >= posl) {
	    l.position -= posl-1;
	    if (l.position < 1) {
		l.length -= (1-l.position);
		l.position = 1;
	    }
	    tag_write(io, ln, l);

	    if (!lastrn) {
		c2.annotations = ln;
		contig_write(io, cont2, c2);
	    } else {
		lastr.next = ln;
		tag_write(io, lastrn, lastr);
	    }
	    
	    lastrn = ln;
	    memcpy(&lastr, &l, sizeof(lastr));
	}

	/* Only in contig1; may need to clip */
	else if (l.position + l.length-1 < posl) {
	    if (l.position + l.length-1 > posr) {
		l.length -= (l.position + l.length-1)-posr;
		tag_write(io, ln, l);
	    }
	    end = ln;
	}

	else {
	    printf("Tag %d is WHERE?\n", ln);
	}
    }

    /* snip left hand contigs annotation list */
    if (end) {
	tag_read(io, end, l);
	if (l.next != 0) {
	    l.next = 0;
	    tag_write(io, end, l);
	}
    } else {
	c1.annotations = 0;
	contig_write(io, cont1, c1);
    }

    return;
}

/*
 * Removes annotations from within a contig over a specified region
 * (posl <= i < posr).
 * Passing zero for posl and posr implies the entire reading.
 */
void remove_contig_tags(GapIO *io, int contig, int posl, int posr) {
    GContigs c;

    contig_read(io, contig, c);
    if (posl || posr) {
	c.annotations = rmanno(io, c.annotations, posl, posr);
    } else {
	c.annotations = rmanno(io, c.annotations, 1, c.length+1);
    }
    contig_write(io, contig, c);
}

/*
 * Removes annotations from within a gel reading over a specified region
 * (posl <= i < posr).
 * Passing zero for posl and posr implies the entire reading.
 */
void remove_gel_tags(GapIO *io, int gel, int posl, int posr) {
    GReadings r;

    gel_read(io, gel, r);
    if (posl || posr) {
	r.annotations = rmanno(io, r.annotations, posl, posr);
    } else {
	r.annotations = rmanno(io, r.annotations, 1, r.length+1);
    }
    gel_write(io, gel, r);
}

/*
 * Removes annotations from a annotation list over a specified region.
 *
 * The annotation list starts at annotation 'anno'. The new list head (it will
 * have changed if we're deleting the first annotation) is returned.
 *
 * Annotations must have pos >= posl and < posr.
 *
 * Note that annotations overlapping this region may need to have their
 * position or length modified, or may need splitting in two. Consider:
 *
 * a -----------------------------------------------
 *      b ----------   c ---------   d----------------------
 *            ^lpos                         ^rpos
 *
 * should yield (a0 and a1 are duplicates):
 *
 * a1---------                               ------- a2
 *      b ----                               --------------- d
 */
int rmanno(GapIO *io, int anno, int lpos, int rpos) {
    GAnnotations a, a2;
    int newtag, last = 0, right = 0, next;
    int list = anno;

    /* scan annotation list */
    while (anno) {
	tag_read(io, anno, a);
	next = a.next;

	/* past rpos => nothing left to do */
	if (a.position >= rpos)
	    break;

	/* within the 'deletion zone' */
	if (a.position >= lpos) {

	    /* entirely? We can remove it if so */
	    if (a.position + a.length <= rpos) {
		/* unlink */
		if (last) {
		    tag_read(io, last, a2);
		    a2.next = a.next;
		    tag_write(io, last, a2);
		} else {
		    list = next;
		}

		/* add to free list */
		delete_tag_rec(io, anno);

		/* prevent last from updating */
		anno = last;

	    } else {
		/* overlaps right end of hole - adjust position/length */
		a.length -= rpos - a.position;
		a.position = rpos;
		tag_write(io, anno, a);
	    }

	} else if (a.position + a.length >= lpos) {
	    /*
	     * Otherwise we start to the left of the hole. If we cover the
	     * hole entirely then we need to split. Overwise just trim.
	     */

	    /* overlaps both ends */
	    if (a.position + a.length > rpos) {
		
		/* duplicate annotation */
		newtag = get_free_tag(io);
		a2.annotation = allocate(io, GT_Text);
		if (a.annotation) {
		    char *tbuf = TextAllocRead(io, a.annotation);
		    if (tbuf) {
			TextWrite(io, a2.annotation, tbuf, strlen(tbuf));
			xfree(tbuf);
		    }
		}
		a2.type = a.type;
		a2.strand = a.strand;

		/* adjust pos/length of new tag */
		a2.position = rpos;
		a2.length = a.position + a.length - a2.position;
		tag_write(io, newtag, a2);

		/* add to annotation list in correct place */
		if (right) {
		    tag_read(io, right, a2);
		    a2.next = newtag;
		    tag_write(io, right, a2);
		} else
		    right = newtag;
	    }

	    /* now trim position of existing tag */
	    a.length = lpos - a.position;
	    tag_write(io, anno, a);
	}
	/* else (pos+len < lpos) => do nothing */

	last = anno;
	anno = next;
    }

    if (last && right) {
	/*
	 * exit to the while loop means that we've either run out of tags,
	 * or that we're now at position >= rpos. Either way, last (if
	 * set) will indicate the last tag to the left of rpos. In
	 * addition the tag held in 'a' will be this 'last' tag.
	 */
	a.next = right;
	tag_write(io, last, a);
    }

    return list;
}

/*
 * Convert a tag string to the numerical and string values stored within it.
 * The tag format is:
 *
 * 0    5   10
 * |----.----|
 * TYPE S start..end
 * optional comment lines. We can have one
 * or several lines.
 *
 * Where 'S' is the strand - "+", "-" or "b" (both).
 * Each line is terminated with a newline ('\n') except for the last one.
 * "start..end" represents an inclusive range.
 *
 * Returns 0 for success, -1 for failure.
 */
int tag2values(char *tag, char *type, int *start, int *end, int *strand,
	       char *comment) {
    int pos;
    char strc, *p;

    /*
     * Parse the header info
     */
    if (sscanf(tag, "%4c %c %d..%d%n", type, &strc, start, end, &pos)
	!= 4)
	return -1;

    if ('+' == strc)
	*strand = 0;
    else if ('-' == strc)
	*strand = 1;
    else
	*strand = 2;

    /*
     * And now the optional comment lines
     */
    /* Normally starts after the newline, but could just be next word */
    for (p = &tag[pos]; *p && (*p == ' ' || *p == '\t'); p++);
    if (*p == '\n')
	p++;
    strcpy(comment, p);
    
    return 0;
}

/*
 * Given numerical and string values of a tag we turn this into a single tag
 * string.
 * The tag format is described above in tag2values().
 *
 * Returns 0 for success, -1 for failure.
 */
int values2tag(char *tag, char *type, int start, int end, int strand,
	       char *comment) {
    int pos;
    char *tg, strc;

    /*
     * Produce the header line
     */
    strc = "+-b"[strand];
    sprintf(tag, "%4s %c %d..%d%n\n", type, strc, start, end, &pos);

    tg = &tag[pos];

    /*
     * Write out the optional comment lines
     */
    if (comment) {
	while (*comment) {
	    *tg++ = '\n';
	    
	    while ('\n' != *comment && '\0' != *comment)
		*tg++ = *comment++;
	    
	    if ('\n' == *comment)
		comment++;
	}
    }

    *tg = '\0';

    return 0;
}


/*
 * Generate an Array of all annotations of a particular type.
 */
Array anno_list(GapIO *io, int type) {
    GReadings r;
    GContigs c;
    GAnnotations a;
    int i, count = 0;
    GCardinal anno;
    Array l;
    struct anno_list_t *ap;

    if (NULL == (l = ArrayCreate(sizeof(struct anno_list_t), 100)))
	return NULL;

    for (i = 1; i <= NumContigs(io); i++) {
	contig_read(io, i, c);

	anno = c.annotations;
	while (anno) {
	    tag_read(io, anno, a);
	    if (type == a.type) {
		ap = ARRP(struct anno_list_t, l, count++);
		ap->anno = anno;
		ap->type = a.type;
		ap->position = a.position;
		ap->length = a.length;
		ap->strand = a.strand;
	    }

	    anno = a.next;
	}
    }

    for (i = 1; i <= NumReadings(io); i++) {
	gel_read(io, i, r);

	anno = r.annotations;
	while (anno) {
	    tag_read(io, anno, a);
	    if (type == a.type) {
		ap = ARRP(struct anno_list_t, l, count++);
		ap->anno = anno;
		ap->type = a.type;
		ap->position = a.position;
		ap->length = a.length;
		ap->strand = a.strand;
	    }

	    anno = a.next;
	}
    }

    return l;
}

/*
 * Removes all annotations listed in anno_ac and anno_av.
 *
 * Returns -1 for error, 0 for success.
 */
int rmanno_list(GapIO *io, int anno_ac, int *anno_av) {
    int i, anno, last_type, last;
    int *del_me;
    GContigs c;
    GReadings r;
    GAnnotations a;

    /*
     * For fast lookup, scan through our anno_av array once to produce
     * single array for all annotations marking whether they're
     * used or not.
     */
    if (NULL == (del_me = (int *)xcalloc(sizeof(int), Nannotations(io)+1)))
	return -1;

    for (i=0; i<anno_ac; i++)
	del_me[anno_av[i]] = 1;

    /* Scan through contigs */
    for (i=1; i <= NumContigs(io); i++) {
	contig_read(io, i, c);
	last_type = 0;
	last = i;
	anno = c.annotations;

	while (anno) {
	    tag_read(io, anno, a);
	    if (del_me[anno]) {
		anno = delete_tag(io, last, anno, last_type);
	    } else {
		last = anno;
		last_type = 2;
		anno = a.next;
	    }
	}
    }

    /* Scan through readings */
    for (i=1; i <= NumReadings(io); i++) {
	gel_read(io, i, r);
	last_type = 1;
	last = i;
	anno = r.annotations;

	while (anno) {
	    tag_read(io, anno, a);
	    if (del_me[anno]) {
		anno = delete_tag(io, last, anno, last_type);
	    } else {
		last = anno;
		last_type = 2;
		anno = a.next;
	    }
	}
    }

    xfree(del_me);
    flush2t(io);
    db_check(io);

    return 0;
}

/*************************************************************
 * FORTRAN interface for routines to handle tags in readings and consensus
 *************************************************************/

/*
 * Moves tags within a contig right by 'NC'. These tags have just been placed
 * at position 'POSN' by PADCON() or UNLNKR()
 */
f_proc_ret shiftt_(f_int *HANDLE, f_int *cont, f_int *posn, f_int *nc) {
    GapIO *io;        

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();
    shift_contig_tags(io, *cont, *posn, *nc);

    f_proc_return();
}

/*
 * Merge tag lists in contigs. We merge tags from CONT2 into CONT1.
 * OFFSET is the relative positions between CONT1 and CONT2.
 */
f_proc_ret mrgtag_(f_int *HANDLE, f_int *CONT1, f_int *CONT2, f_int *OFF) {
    GapIO *io;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();

    merge_contig_tags(io, *CONT1, *CONT2, *OFF);
    f_proc_return();
}

/*
 * Complement the tags in a contig.
 *
 * This simply involves negating the position of each one and reversing the
 * order of the list.
 */
f_proc_ret comtag_(f_int *HANDLE, f_int *CONT, f_int *LEN) {
    GapIO *io;

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();
    if (*LEN != io_clength(io, *CONT)) {
	fprintf(stderr, "BUG at %s:%d\n", __FILE__, __LINE__);
    }
    complement_contig_tags(io, *CONT);
    
    f_proc_return();
}

/*
 * Called by "break contig". When we're splitting a contig in half we need
 * to move the annotations too.
 * For annotations that overlap the two contigs we duplicate and adjust the
 * lengths.
 * The new contigs used to overlap between POSL and POSR.
 */
f_proc_ret spltag_(f_int *HANDLE, f_int *CONT1, f_int *CONT2,
		   f_int *POSL, f_int *POSR) {
    GapIO *io;    

    if ( (io = io_handle(HANDLE)) == NULL) f_proc_return();
    split_contig_tags(io, *CONT1, *CONT2, *POSL, *POSR);

    f_proc_return();
}

/*
 * Removes annotations from within a contig over a specified region
 * (posl <= i < posr).
 * Passing zero for posl and posr implies the entire reading.
 */
f_proc_ret rmctag_(f_int *handle, f_int *cont, f_int *lpos, f_int *rpos) {
    GapIO *io;

    if ( (io = io_handle(handle)) == NULL) f_proc_return();
    remove_contig_tags(io, *cont, *lpos, *rpos);

    f_proc_return();
}

/*
 * Removes annotations from within a gel reading over a specified region
 * (posl <= i < posr).
 * Passing zero for posl and posr implies the entire reading.
 */
f_proc_ret rmgtag_(f_int *handle, f_int *gel, f_int *lpos, f_int *rpos) {
    GapIO *io;

    if ( (io = io_handle(handle)) == NULL) f_proc_return();
    remove_gel_tags(io, *gel, *lpos, *rpos);
    
    f_proc_return();
}

void write_tags(GapIO *io, char *fname, int num_tags,
		int *read1, int *pos1, int *read2, int *pos2, int *length) {
    register int n = num_tags, i;
    char buf[100], com[100], n1[DB_NAMELEN+1], n2[DB_NAMELEN+1];
    FILE *fp;
    Exp_info *e;

    if (NULL == (fp = fopen(fname, "w"))) {
	verror(ERR_WARN, "write_tag", "Failed to open file %s\n", fname);
	return;
    }
    e = exp_create_info();
    e->fp = fp;

    for (i=0; i<n; i++) {
	/* take the absolute values because reading can be negative if it
	 * contains an inverted repeat
	 */
	f_int r1 = abs(read1[i]);
	f_int r2 = abs(read2[i]);

	/* get reading name from contigs */
	readn_(handle_io(io), &r1, n1, DB_NAMELEN);
	Fstr2Cstr(n1, DB_NAMELEN, n1, DB_NAMELEN+1);
	readn_(handle_io(io), &r2, n2, DB_NAMELEN);
	Fstr2Cstr(n2, DB_NAMELEN, n2, DB_NAMELEN+1);

	sprintf(buf, "Repeat number %d, end 1", i);
	exp_put_str(e, EFLT_CC, buf, strlen(buf));

	/* write tag 1 */
	exp_put_str(e, EFLT_ID, n1, strlen(n1));
	sprintf(com, "Repeats with contig %s, offset %d",
		n2, pos2[i]);
	values2tag(buf, "REPT", pos1[i],
		   pos1[i] + length[i] - 1, 2, com);
	exp_put_str(e, EFLT_TC, buf, strlen(buf));
	    
	sprintf(buf, "Repeat number %d, end 2", i);
	exp_put_str(e, EFLT_CC, buf, strlen(buf));

	/* write tag 2 */
	exp_put_str(e, EFLT_ID, n2, strlen(n2));
	sprintf(com, "Repeats with contig %s, offset %d",
		n1, pos1[i]);
	values2tag(buf, "REPT", pos2[i],
		   pos2[i] + length[i] - 1, 2, com);
	exp_put_str(e, EFLT_TC, buf, strlen(buf));

	exp_put_str(e, EFLT_CC, "", 0);
    }

    /* deallocate experiment file info - will close fp for us */
    exp_destroy_info(e);

    return;
}

/*
 * ---------------------------------------------------------------------------
 * Operations to find and edit tags.
 */

/*
 * Add to an array of tags of a particular type. Searches within one
 * contig within a given range.
 */
int find_tags_contig(GapIO *io, int contig, int start, int end, Array al,
		     int *itypes, int num_t) {
    int rnum, anno;
    GReadings r;
    GAnnotations a;
    int offset;
    int i, prev;

    for (rnum = io_clnbr(io, contig); rnum; rnum = io_rnbr(io, rnum)) {
	/* Prune search to readings within this range */
	if (io_relpos(io, rnum) + ABS(io_length(io, rnum)) - 1 < start)
	    continue;

	if (io_relpos(io, rnum) > end)
	    break;

	if (-1 == io_read_annotation(io, rnum, &anno))
	    continue;
	
	gel_read(io, rnum, r);
	offset = r.position - r.start - 1;

	prev = 0;
	while (anno) {
	    int tstart, tend;
	    tag_read(io, anno, a);

	    if (r.sense) {
		tstart = r.length-1 - a.position + offset;
		tend = r.length-1 - a.position + a.length-1 + offset;
	    } else {
		tstart = a.position + offset;
		tend = a.position + a.length-1 + offset;
	    }
	    if (tend >= start && tstart <= end) {
		for (i = 0; i < num_t; i++) {
		    if (itypes[i] == a.type) {
			anno_ptr *p = ArrayRef(al, ArrayMax(al));
			p->number = anno;
			p->prev = prev;
			p->next = a.next;
			p->from = rnum;
			p->from_type = GT_Readings;
		    }
		}
	    }
	    prev = anno;
	    anno = a.next;
	}
    }

    return 0;
}

/*
 * Return an array of tags of a particular type. Searchs within a list of
 * contig ranges. The Array is of type "GAnnotation *".
 */
Array find_tags(GapIO *io, contig_list_t *contigs, int num_contigs,
		char **type, int num_t) {
    int i;
    Array al;
    int *itypes;
    
    if (NULL == (al = ArrayCreate(sizeof(anno_ptr), 100)))
	return NULL;

    if (NULL == (itypes = (int *)xmalloc(num_t * sizeof(int)))) {
	ArrayDestroy(al);
	return NULL;
    }

    for (i = 0; i < num_t; i++) {
	itypes[i] = str2type(type[i]);
    }

    for (i = 0; i < num_contigs; i++) {
	find_tags_contig(io, contigs[i].contig, contigs[i].start,
			 contigs[i].end, al, itypes, num_t);
    }
    printf("%d tags\n", ArrayMax(al));
    for (i = 0; i < ArrayMax(al); i++) {
	printf("Tag %d, prev %d, next %d, read %d, type %d\n",
	       arr(anno_ptr, al, i).number,
	       arr(anno_ptr, al, i).prev,
	       arr(anno_ptr, al, i).next,
	       arr(anno_ptr, al, i).from,
	       arr(anno_ptr, al, i).from_type);
    }

    xfree(itypes);

    return al;
}
