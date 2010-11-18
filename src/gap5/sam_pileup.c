#include <staden_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam_pileup.h"

/* --------------------------------------------------------------------------
 * The pileup code itself. 
 *
 * This consists of the external pileup_loop() function, which takes a
 * sam/bam samfile_t pointer and a callback function. The callback function
 * is called once per column of aligned data (so once per base in an
 * insertion).
 */

/*
 * Fast conversion from encoded SAM base nibble to a printable character
 */
static char tab[256][2];
static void init_tab(void) {
    int i, j;
    unsigned char b2;
    static int done = 0;

    if (done)
	return;

    for (i = 0; i < 16; i++) {
	for (j = 0; j < 16; j++) {
	    b2 = (i<<4) | j;
	    tab[b2][0] = "NACMGRSVTWYHKDBN"[i];
	    tab[b2][1] = "NACMGRSVTWYHKDBN"[j];
	}
    }

    done = 1;
}


/*
 * Fetches the next base => the nth base at unpadded position pos. (Nth can
 * be greater than 0 if we have an insertion in this column). Do not call this
 * with pos/nth lower than the previous query, although higher is better.
 * (This allows it to be initialised at base 0.)
 *
 * Stores the result in base and also updates is_insert to indicate that
 * this sequence still has more bases in this position beyond the current
 * nth parameter.
 *
 * Returns 1 if a base was fetched
 *         0 if not (eg ran off the end of sequence)
 */
static int get_next_base(pileup_t *p, int pos, int nth, int *is_insert) {
    bam_seq_t *b = p->b;
    enum cigar_op op = p->cigar_op;

    *is_insert = 0;

    /* Find pos first */
    while (p->pos < pos) {
	p->nth = 0;

	if (p->cigar_len == 0) {
	    if (p->cigar_ind >= bam_cigar_len(b)) {
		p->eof = 1;
		return 0;
	    }

	    op=p->cigar_op  = p->b_cigar[p->cigar_ind] & BAM_CIGAR_MASK;
	    p->cigar_len = p->b_cigar[p->cigar_ind] >> BAM_CIGAR_SHIFT;
	    p->cigar_ind++;
	}
	
	if (op == BAM_CMATCH && p->cigar_len <= pos - p->pos) {
	    p->seq_offset += p->cigar_len;
	    p->pos += p->cigar_len;
	    p->cigar_len = 0;
	} else {
	    switch (op) {
	    case BAM_CMATCH:
	    case BAM_CBASE_MATCH:
	    case BAM_CBASE_MISMATCH:
		p->seq_offset++;
		/* Fall through */
	    case BAM_CDEL:
		p->pos++;
		p->cigar_len--;
		break;

	    case BAM_CINS:
	    case BAM_CSOFT_CLIP:
		p->seq_offset += p->cigar_len;
		/* Fall through */
	    case BAM_CPAD:
	    case BAM_CHARD_CLIP:
		p->cigar_len = 0;
		break;

	    case BAM_CREF_SKIP:
	    default:
		fprintf(stderr, "Unhandled cigar_op %d\n", op);
		return -1;
	    }
	}
    }

    /* Now at pos, find nth base */
    while (p->nth < nth) {
	if (p->cigar_len == 0) {
	    if (p->cigar_ind >= bam_cigar_len(b)) {
		p->eof = 1;
		return 0; /* off end of seq */
	    }

	    op=p->cigar_op  = p->b_cigar[p->cigar_ind] & BAM_CIGAR_MASK;
	    p->cigar_len = p->b_cigar[p->cigar_ind] >> BAM_CIGAR_SHIFT;
	    p->cigar_ind++;
	}

	switch (op) {
	case BAM_CMATCH:
	case BAM_CBASE_MATCH:
	case BAM_CBASE_MISMATCH:
	case BAM_CSOFT_CLIP:
	case BAM_CDEL:
	    goto at_nth; /* sorry, but it's fast! */

	case BAM_CINS:
	    p->seq_offset++;
	    /* Fall through */
	case BAM_CPAD:
	    p->cigar_len--;
	    p->nth++;
	    break;

	case BAM_CHARD_CLIP:
	    p->cigar_len = 0;
	    break;

	case BAM_CREF_SKIP:
	default:
	    fprintf(stderr, "Unhandled cigar_op %d\n", op);
	    return -1;
	}
    }
 at_nth:

    /* Fill out base & qual fields */
    if (p->nth < nth && op != BAM_CINS) {
	//p->base = '-';
	p->base = '*';
	p->qual = p->b_qual[p->seq_offset];
	if (p->seq_offset < b->len)
	    p->qual = (p->qual + p->b_qual[p->seq_offset+1])/2;
    } else {
	switch(op) {
	case BAM_CDEL:
	    p->base = '*';
	    p->qual = p->b_qual[p->seq_offset];
	    if (p->seq_offset < b->len)
		p->qual = (p->qual + p->b_qual[p->seq_offset+1])/2;
	    break;

	case BAM_CPAD:
	    //p->base = '+';
	    p->base = '*';
	    p->qual = p->b_qual[p->seq_offset];
	    if (p->seq_offset < b->len)
		p->qual = (p->qual + p->b_qual[p->seq_offset+1])/2;
	    break;

	default:
	    p->qual = p->b_qual[p->seq_offset];
	    p->base = tab[p->b_seq[p->seq_offset/2]][p->seq_offset&1];
	    break;
	}
    }

    /* Check if next op is an insertion of some sort */
    if (p->cigar_len == 0) {
	if (p->cigar_ind < bam_cigar_len(b)) {
	    op=p->cigar_op  = p->b_cigar[p->cigar_ind] & BAM_CIGAR_MASK;
	    p->cigar_len = p->b_cigar[p->cigar_ind] >> BAM_CIGAR_SHIFT;
	    p->cigar_ind++;
	} else {
	    p->eof = 1;
	}
    }

    switch (op) {
    case BAM_CPAD:
    case BAM_CINS:
	*is_insert = p->cigar_len;
	break;

    case BAM_CSOFT_CLIP:
	p->eof = 1;
	break;

    default:
	break;
    }
    
    return 1;
}


/*
 * Loops through a set of supplied ranges producing columns of data.
 * When found, it calls func with clientdata as a callback. Func should
 * return 0 for success and non-zero for failure. seq_init() is called
 * on each new entry before we start processing it. It should return 0 or 1
 * to indicate reject or accept status (eg to filter unmapped data).
 * If seq_init() returns -1 we abort the pileup_loop with an error.
 * seq_init may be NULL.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int pileup_loop(bam_file_t *fp,
		int (*seq_init)(void *client_data,
				bam_file_t *fp,
				pileup_t *p),
		int (*func)(void *client_data,
			    bam_file_t *fp,
			    pileup_t *p,
			    int depth,
			    int pos,
			    int nth),
		void *client_data) {
    int ret = -1;
    pileup_t *phead = NULL, *p, *pfree = NULL, *last, *next;
    pileup_t *pnew = NULL;
    int is_insert, nth = 0;
    int col = 0, r;
    int last_ref = -1;

    /* FIXME: allow for start/stop boundaries rather than consuming all data */

    init_tab();
    if (NULL == (pnew = calloc(1, sizeof(*p))))
	return -1;
    
    do {
	bam_seq_t *b;
	int pos, last_in_contig;

	r = (bam_next_seq(fp, &pnew->b) > 0) ? 1 : -1;
	b = pnew->b;

	/* Force realloc */
	fp->bs = NULL;
	fp->bs_size = 0;

	//r = samread(fp, pnew->b);
	if (r >= 0) {
	    if (bam_flag(b) & BAM_FUNMAP)
		continue;

	    if (b->ref == last_ref) {
		pos = b->pos+1;
		//printf("New seq at pos %d @ %d %s\n", pos, b->ref,
		//       bam_name(b));
		last_in_contig = 0;
	    } else {
		//printf("New ctg at pos %d @ %d\n",b->pos+1,b->ref);
		pos = (b->pos > col ? b->pos : col)+1;
		last_in_contig = 1;
	    }
	} else {
	    last_in_contig = 1;
	    pos = col+1;
	}

	/* Process data between the last column and our latest addition */
	while (col < pos && phead) {
	    int v, ins, depth = 0;

	    /* Pileup */
	    is_insert = 0;
	    for (p = phead; p; p = p->next) {
		if (!get_next_base(p, col, nth, &ins))
		    p->eof = 2;

		if (is_insert < ins)
		    is_insert = ins;
		
		depth++;
	    }

	    /* Call our function on phead linked list */
	    v = func(client_data, fp, phead, depth, col, nth);

	    /* Remove dead seqs */
	    for (p = phead, last = NULL; p; p = next) {
		next = p->next;
		
		p->start = 0;
		if (p->eof) {
		    if (last)
			last->next = p->next;
		    else
			phead = p->next;

		    p->next = pfree;
		    pfree = p;
		    
		    //printf("Del seq %s at pos %d\n", bam_name(p->b), col);
		} else {
		    last = p;
		}
	    }

	    if (v == 1)
		break; /* early abort */
	    
	    if (v != 0)
		goto error;

	    /* Next column */
	    if (is_insert) {
		nth++;
	    } else {
		nth = 0;
		col++;
	    }

	    /* Special case for the last sequence in a contig */
	    if (last_in_contig && phead)
		pos++;
	}

	/* May happen if we have a hole in the contig */
	col = pos;

	/* New contig */
	if (b && b->ref != last_ref) {
	    last_ref = b->ref;
	    pos = b->pos+1;
	    nth = 0;
	    col = pos;
	}

	/* Add this seq */
	if (r >= 0) {
	    p = pnew;
	    p->next       = phead;
	    p->cd         = NULL;
	    p->start      = 1;
	    p->eof        = 0;
	    p->pos        = pos-1;
	    p->cigar_len  = 0;
	    p->cigar_ind  = 0;
	    p->cigar_op   = 'X';
	    p->seq_offset = -1;
	    p->b_strand   = bam_strand(p->b) ? 1 : 0;
	    p->b_qual     = bam_qual(p->b);
	    p->b_seq      = bam_seq(p->b);
	    p->b_cigar    = bam_cigar(p->b);

	    if (seq_init) {
		int v;
		v = seq_init(client_data, fp, p);
		if (v == -1)
		    return -1;
		
		if (v == 1) {
		    /* Keep this seq */
		    phead = p;
		} else {
		    /* Push back on free list */
		    p->next = pfree;
		    pfree = p;
		}
	    } else {
		phead = p;
	    }

	    /* Allocate the next pileup rec */
	    if (pfree) {
		pnew = pfree;
		pfree = pfree->next;
	    } else {
		if (NULL == (pnew = calloc(1, sizeof(*pnew))))
		    goto error;
	    }
	}
    } while (r >= 0);

    ret = 0;
 error:

    if (pnew) {
	//free(pnew->b.data);
	free(pnew);
    }

    /* Tidy up */
    for (p = pfree; p; p = next) {
	next = p->next;
	//free(p->b.data);
	free(p);
    }

    return ret;
}

#ifdef TEST_MAIN
/* --------------------------------------------------------------------------
 * Example usage of the above pileup code
 */

#include <ctype.h>
#define MIN(a,b) ((a)<(b)?(a):(b))

static char *append_int(char *cp, int i) {
    int j;

    if (i < 0) {
	*cp++ = '-';
	i = -i;
    } else if (i == 0) {
	*cp++ = '0';
	return cp;
    }

    //if (i < 10)         goto b0;
    if (i < 100)        goto b1;
    //if (i < 1000)       goto b2;
    if (i < 10000)      goto b3;
    //if (i < 100000)     goto b4;
    if (i < 1000000)    goto b5;
    //if (i < 10000000)   goto b6;
    if (i < 100000000)  goto b7;

 b9: if ((j = i / 1000000000)) {*cp++ = j + '0'; i -= j*1000000000; goto x8;}
 b8: if ((j = i / 100000000))  {*cp++ = j + '0'; i -= j*100000000;  goto x7;}
 b7: if ((j = i / 10000000))   {*cp++ = j + '0'; i -= j*10000000;   goto x6;}
 b6: if ((j = i / 1000000))    {*cp++ = j + '0', i -= j*1000000;    goto x5;}
 b5: if ((j = i / 100000))     {*cp++ = j + '0', i -= j*100000;     goto x4;}
 b4: if ((j = i / 10000))      {*cp++ = j + '0', i -= j*10000;      goto x3;}
 b3: if ((j = i / 1000))       {*cp++ = j + '0', i -= j*1000;       goto x2;}
 b2: if ((j = i / 100))        {*cp++ = j + '0', i -= j*100;        goto x1;}
 b1: if ((j = i / 10))         {*cp++ = j + '0', i -= j*10;         goto x0;}
 b0: if (i)                     *cp++ = i + '0';
    return cp;

 x8: *cp++ = i / 100000000 + '0', i %= 100000000;
 x7: *cp++ = i / 10000000  + '0', i %= 10000000;
 x6: *cp++ = i / 1000000   + '0', i %= 1000000;
 x5: *cp++ = i / 100000    + '0', i %= 100000;
 x4: *cp++ = i / 10000     + '0', i %= 10000;
 x3: *cp++ = i / 1000      + '0', i %= 1000;
 x2: *cp++ = i / 100       + '0', i %= 100;
 x1: *cp++ = i / 10        + '0', i %= 10;
 x0: *cp++ = i             + '0';

    return cp;
}

char strand_char[2][256];
void strand_init(void) {
    int i;
    for (i = 0; i < 256; i++) {
	strand_char[0][i] = toupper((unsigned char)i);
	strand_char[1][i] = tolower((unsigned char)i);
    }
}

/* MAX_DEPTH only has a consequence on the example / test dump output */
#define MAX_DEPTH 4096
static int sam_pileup(void *cd, bam_file_t *fp, pileup_t *p,
		      int depth, int pos, int nth) {
    char seq[MAX_DEPTH*3], *sp = seq, qual[MAX_DEPTH], *qp = qual;
    char buf[MAX_DEPTH*2+100], *cp = buf;
    int ref;

    if (!p)
	return 0;

    ref = p->b->ref;
    for (; p; p = p->next) {
	if (p->start) {
	    *sp++ = '^';
	    *sp++ = MIN(p->qual,93) + '!';
	}
	*sp++ = strand_char[p->b_strand][p->base];
	//*sp++ = strand_char[bam_strand(p->b)][p->base];
	if (p->eof)
	    *sp++ = '$';
	*qp++ = p->qual + '!';
    }

    /* Equivalent to the printf below, but faster */
    cp = buf;
    strcpy(cp, fp->ref[ref].name); cp += strlen(cp);
    *cp++ = '\t';
    cp = append_int(cp, pos);   *cp++ = '\t';
    *cp++ = 'N';
    *cp++ = '\t';
    cp = append_int(cp, depth); *cp++ = '\t';
    memcpy(cp, seq,  sp-seq);  cp += sp-seq;  *cp++ = '\t';
    memcpy(cp, qual, qp-qual); cp += qp-qual; *cp++ = '\0';
    puts(buf);

    //*sp++ = 0;
    //*qp++ = 0;
    //printf("ref\t%d+%d\tN\t%d\t%s\t%s\n", pos, nth, depth, seq, qual);

    return 0;
}

static int basic_pileup(void *cd, bam_file_t *fp, pileup_t *p,
			int depth, int pos, int nth) {
    char seq[MAX_DEPTH*3], *sp = seq, qual[MAX_DEPTH], *qp = qual;
    char buf[MAX_DEPTH*4+100], *cp = buf, *rp;
    int ref;

    if (!p)
	return 0;

    /* Ref, pos, depth */
    ref = p->b->ref;
    rp = fp->ref[ref].name;
    while (*cp++ = *rp++)
	;
    cp--;
    *cp++ = '\t';
    cp = append_int(cp, pos);   *cp++ = '+';
    cp = append_int(cp, nth);   *cp++ = '\t';
    *cp++ = 'N';
    *cp++ = '\t';
    cp = append_int(cp, depth); *cp++ = '\t';

    /* Seq + qual at predetermined offsets */
    qp = cp + depth + 1;
    for (; p; p = p->next) {
	*cp++ = p->base;
	*qp++ = p->qual + '!';
    }
    *cp++ = '\t';
    *qp++ = '\0';
    puts(buf);

    return 0;
}

int null_pileup(void *cd, bam_file_t *fp, pileup_t *p,
		int depth, int pos, int nth) {
    return 0;
}

int main(int argc, char **argv) {
    bam_file_t *fp;

    if (argc != 3) {
	fprintf(stderr, "sam_pileup filename mode\n");
	return 1;
    }

    strand_init();

    fp = bam_open(argv[1], argv[2]);
    if (!fp) {
	perror(argv[1]);
	return 1;
    }
    
    pileup_loop(fp, NULL, basic_pileup, fp);
    //pileup_loop(fp, NULL, sam_pileup, fp);
    //pileup_loop(fp, NULL, null_pileup, NULL);

    bam_close(fp);
    return 0;
}
#endif
