/* ---- Implements importing from BAF (basic alignment format). ---- */

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "tg_gio.h"
#include "hache_table.h"
#include "baf.h"

/*
 * Deallocates aline_t struct.
 */
void free_line(line_t *l) {
    if (l) {
	if (l->str)
	    free(l->str);
	free(l);
    }
}


/*
 * Produces a printable version of a line-type code.
 * Note this returns a static buffer so it should not be freed and it
 * is valid only up until the next call to linetype2str.
 */
char *linetype2str(int lt) {
    static unsigned char type[3];
    type[0] = lt>>8;
    type[1] = lt&0xff;
    type[2] = 0;
    return (char *)type;
}


/*
 * Gets a line of arbitrary length from the opened file.
 * It stores it in the supplied line_t struct or if NULL allocates it's
 * own version.
 *
 * The line of text returned is of arbitrary length and nul terminated.
 * It will not include the newline character.
 *
 * Returns a line_t pointer on success
 *           NULL on failure (check using feof/ferror for why).
 */
#define BLKSZ 1024
line_t *get_line(FILE *fp, line_t *in) {
    line_t *l = in;

    if (!l) {
	l = (line_t *)malloc(sizeof(*l));
	l->str = NULL;
	l->len = 0;
    }

    do {
	size_t pos = 0, len;
	int no_nl = 1;

	do {
	    if (l->len-pos < BLKSZ) {
		l->len = BLKSZ+pos;
		if (NULL == (l->str = realloc(l->str, l->len))) {
		    if (!in) free_line(l);
		    return NULL;
		}
	    }
	    if (NULL == fgets(&l->str[pos], BLKSZ, fp)) {
		if (!in) free_line(l);
		return NULL;
	    }
	    len = strlen(&l->str[pos]);
	    if (l->str[pos+len-1] == '\n') {
		l->str[pos+len-1] = 0;
		len--;
		no_nl = 0;
	    }
	    pos += len;
	} while (no_nl);

	l->len = pos;
    } while (l->str[0] == '#');

    if (l->len == 0) {
	/* Blank line indicates block termination */
	l->assign = 0;
	l->value = NULL;
	l->type = XX;

	return l;
    }

    if (l->len < 3 || (l->str[2] != '=' && l->str[2] != ':')) {
	fprintf(stderr, "Malformed line '%s'\n", l->str);
	if (!in) free_line(l);
	return NULL;
    }

    /* If we allocated this line ourselves, give back any extra memory */
    if (!in) {
	l->str = realloc(l->str, l->len+1);
    }

    l->type = CC2(l->str[0], l->str[1]);
    l->assign = l->str[2];
    l->value = l->str+3;

    return l;
}

baf_block *baf_next_block(FILE *fp) {
    line_t *l;
    baf_block *b;
    int order = 0;

    if (NULL == (l = get_line(fp, NULL))) {
	/* Empty block */
	return NULL;
    }
    if (NULL == (b = (baf_block *)calloc(1, sizeof(*b))))
	return NULL;

    b->type = l->type;
    b->h = HacheTableCreate(0, HASH_DYNAMIC_SIZE | HASH_ALLOW_DUP_KEYS);

    do {
	HacheData hd;
	
	if (l->type == XX) {
	    free_line(l);
	    break;
	}

	l->order = order++;
	hd.p = l;
	HacheTableAdd(b->h, (char *)&l->type, sizeof(l->type), hd, NULL);
    } while (l = get_line(fp, NULL));

    return b;
}

void baf_block_destroy(baf_block *b) {
    if (!b)
	return;

    if (b->h) {
	HacheIter *iter = HacheTableIterCreate();
	HacheItem *hi;
	while (hi = HacheTableIterNext(b->h, iter)) {
	    line_t *l = hi->data.p;
	    if (l) free_line(l);
	}
	HacheTableIterDestroy(iter);
	HacheTableDestroy(b->h, 0);
    }

    free(b);
}

line_t *baf_block_line(baf_block *b, int type) {
    HacheItem *hi = HacheTableSearch(b->h, (char *)&type, sizeof(type));
    return hi ? (line_t *)hi->data.p : NULL;
}

char *baf_block_value(baf_block *b, int type) {
    line_t *l = baf_block_line(b, type);
    return l ? l->value : NULL;
}


int construct_seq_from_block(seq_t *s, baf_block *b, char **tname) {
    int ap, dir, cleft, cright, i, qb, mq, end;
    size_t len;
    char *cp, *seq, *qual, *name, *trace_name, *alignment;
    
    /* Prevent valgrind from complaining about portions of seq being unset */
    memset(s, 0, sizeof(*s));

    /* Extract various sections from the block */
    name = baf_block_value(b, RD);
    seq  = baf_block_value(b, SQ);
    qual = baf_block_value(b, FQ);
    if (NULL == (trace_name = baf_block_value(b, TR)))
	trace_name = "";
    if (NULL == (alignment = baf_block_value(b, AL)))
	alignment = "";

    if (!name || !seq || !qual)
	return -1;

    len  = strlen(seq);

    if ((cp = baf_block_value(b, AP)))
	ap = atoi(cp);
    else
	return -1;

    if (NULL == (*tname = baf_block_value(b, TN)))
	*tname = name;

    if ((cp = baf_block_value(b, QL)))
	cleft = atoi(cp);
    else
	cleft = 0;

    if ((cp = baf_block_value(b, QR)))
	cright = atoi(cp);
    else
	cright = len;

    if ((cp = baf_block_value(b, DR)))
	dir = atoi(cp);
    else
	dir = 1;

    if ((cp = baf_block_value(b, PR)))
	end = atoi(cp);
    else
	end = 0;

    if ((cp = baf_block_value(b, MQ)))
	mq = atoi(cp);
    else
	mq = 50;

    /* From fastq to binary array */
    for (i = 0; i < len; i++)
	qual[i] -= 33;

    qb = 1;
    s->format = SEQ_FORMAT_CNF1; /* pack bytes */
    for (i = 0; i < len; i++) {
	if (seq[i] == '-')
	    seq[i] = '*';
	if (seq[i] == 'n' || seq[i] == 'N')
	    seq[i] = '-';
    }

#if 0
    /* In MAQ style padding qual 0 implies N */
    for (i = 0; i < len; i++) {
	if (qual[i] == 0)
	    qual[i] = 1;
	if (seq[i] != 'A' && seq[i] != 'C' &&
	    seq[i] != 'G' && seq[i] != 'T')
	    qual[i] = 0;
    }
#endif

    s->pos = ap;
    s->len = dir * len;
    s->flags = s->len < 0 ? SEQ_COMPLEMENTED : 0;
    if (end == 1)
	s->flags |= SEQ_END_REV;
    s->left = cleft;
    s->right = cright;
    s->mapping_qual = mq;
    if (dir == 1)
	s->pos -= s->left-1;
    else
	s->pos -= -s->len - s->right;

    s->name_len = strlen(name);
    s->trace_name_len = strlen(trace_name);
    s->alignment_len = strlen(alignment);

    s->name = s->data = (char *)malloc(s->name_len+1 +
				       s->trace_name_len+1 +
				       s->alignment_len+1 +
				       len+qb*len);
    strcpy(s->name, name);

    s->trace_name = s->name + s->name_len + 1;
    strcpy(s->trace_name, trace_name);

    s->alignment = s->trace_name + s->trace_name_len + 1;
    strcpy(s->alignment, alignment);

    s->seq = s->alignment + s->alignment_len + 1;
    memcpy(s->seq, seq, len);

    s->conf = s->seq + len;
    memcpy(s->conf, qual, (s->format == SEQ_FORMAT_CNF4 ? 4 : 1) * len);

    return 0;
}

typedef struct {
    int rec;
    int bin;
    int idx;
    int crec;
} pair_loc_t;

int parse_baf(GapIO *io, char *fn, int no_tree, int pair_reads,
	      int merge_contigs) {
    int nseqs = 0;
    struct stat sb;
    FILE *fp;
    off_t pos;
    contig_t *c;
    HacheTable *pair = NULL;
    int ncontigs = 0;
    baf_block *b, *co = NULL;
    int last_obj_type = 0;
    int last_obj_pos = 0;
    int last_obj_rec = 0;
    int last_obj_orient = 0;
	
    printf("Loading %s...\n", fn);
    if (-1 == stat(fn, &sb) ||
	NULL == (fp = fopen(fn, "r"))) {
	perror(fn);
	return -1;
    }

    if (pair_reads) {
	pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
    }

    /* Loop:
     * Read 1 block of data.
     * If contig, create contig
     * If read, insert it, insert to index.
     * Anything else - reject for now
     */
    pos = 0;
    while (b = baf_next_block(fp)) {
	int delay_destroy = 0;

	switch (b->type) {
	case CO: {
	    char *contig = baf_block_value(b, CO);
	    if (co)
		baf_block_destroy(co);
	    co = b;
	    delay_destroy = 1;

	    ncontigs++;
	    if (!merge_contigs ||
		NULL == (c = find_contig_by_name(io, contig))) {
		c = contig_new(io, contig);
		contig_index_update(io, contig, strlen(contig), c->rec);
	    }
	    cache_incr(io, c);

	    /* For anno */
	    last_obj_type = GT_Contig;
	    last_obj_rec = c->rec;
	    last_obj_pos = c->start + 1;
	    last_obj_orient = 0;

	    break;
	}

	case RD: {
	    seq_t seq;
	    range_t r, *r_out;
	    int recno;
	    bin_index_t *bin;
	    HacheItem *hi;
	    char *tname;
	    pair_loc_t *pl = NULL;

	    /* Construct seq struct */
	    if (-1 == construct_seq_from_block(&seq, b, &tname)) {
		fprintf(stderr, "Failed to parse read block for seq %d\n",
			nseqs);
		break;
	    }

	    /* Create range */
	    r.start = seq.pos;
	    r.end   = seq.pos + (seq.len > 0 ? seq.len : -seq.len) - 1;
	    r.rec   = 0;
	    r.pair_rec = 0;
	    r.mqual = seq.mapping_qual;
	    r.flags = GRANGE_FLAG_TYPE_SINGLE;
	    if (seq.flags & SEQ_END_REV)
		r.flags |= GRANGE_FLAG_END_REV;
	    else
		r.flags |= GRANGE_FLAG_END_FWD;
	    if (seq.len < 0)
		r.flags |= GRANGE_FLAG_COMP1;

	    /* Add the range to a bin, and see which bin it was */
	    bin = bin_add_range(io, &c, &r, &r_out);

	    /* Save sequence */
	    seq.bin = bin->rec;
	    seq.bin_index = r_out - ArrayBase(range_t, bin->rng);
	    recno = sequence_new_from(io, &seq);

	    /* For anno */
	    last_obj_type = GT_Seq;
	    last_obj_rec = recno;
	    if (seq.len >= 0) {
		last_obj_pos = seq.pos;
		last_obj_orient = 0;
	    } else {
		last_obj_pos = seq.pos - seq.len - 1;
		last_obj_orient = 1;
	    }

	    /* Find pair if appropriate */
	    if (pair) {
		int new = 0;
		HacheData hd;

		/* Add data for this end */
		pl = (pair_loc_t *)malloc(sizeof(*pl));
		pl->rec  = recno;
		pl->bin  = bin->rec;
		pl->crec = c->rec;
		pl->idx  = seq.bin_index;
		hd.p = pl;

		hi = HacheTableAdd(pair, tname, strlen(tname), hd, &new);

		/* Pair existed already */
		if (!new) {
		    pair_loc_t *po = (pair_loc_t *)hi->data.p;
		    bin_index_t *bo;
		    range_t *ro;

		    /* We found one so update r_out now, before flush */
		    r_out->flags &= ~GRANGE_FLAG_TYPE_MASK;
		    r_out->flags |=  GRANGE_FLAG_TYPE_PAIRED;
		    r_out->pair_rec = po->rec;

		    /* Link other end to 'us' too */
		    bo = (bin_index_t *)cache_search(io, GT_Bin, po->bin);
		    bo = cache_rw(io, bo);
		    bo->flags |= BIN_RANGE_UPDATED;
		    ro = arrp(range_t, bo->rng, po->idx);
		    ro->flags &= ~GRANGE_FLAG_TYPE_MASK;
		    ro->flags |=  GRANGE_FLAG_TYPE_PAIRED;
		    ro->pair_rec = pl->rec;

#if 0
		    if (c->rec == po->crec) {
			ro->flags    |= GRANGE_FLAG_CONTIG;
			r_out->flags |= GRANGE_FLAG_CONTIG;
		    }
#endif

		    /* And, making an assumption, remove from hache */
		    HacheTableDel(pair, hi, 1);
		    free(pl);
		}
	    }

	    if (!no_tree)
		sequence_index_update(io, seq.name, seq.name_len, recno);
	    free(seq.data);

	    /* Link bin back to sequence too before it gets flushed */
	    r_out->rec = recno;
	    
	    if ((++nseqs & 0x3fff) == 0) {
		int perc = 0;
		if (nseqs & 0xffff == 0)
		    cache_flush(io);
		pos = ftello(fp);
		perc = 100.0 * pos / sb.st_size;
		printf("\r%d%%", perc);
		fflush(stdout);
	    }
	    
	    break;
	}

	case AN: {
	    range_t r;
	    anno_ele_t *e;
	    char *loc = baf_block_value(b, LO);
	    char *len = baf_block_value(b, LL);
	    char *txt = baf_block_value(b, TX);
	    int pos;
	    bin_index_t *bin;

	    if (!loc) {
		pos = last_obj_pos;
	    } else {
		if (*loc == '@') {
		    pos = atoi(loc+1)-1;
		} else {
		    if (last_obj_orient == 0)
			pos = last_obj_pos + atoi(loc)-1;
		    else
			pos = last_obj_pos - (atoi(loc)-1)
			    - (len ? atoi(len)-1 : 0);
		}
	    }

	    r.start = pos;
	    r.end = pos + (len ? atoi(len)-1 : 0);

	    r.mqual    = last_obj_type;  /* obj_type */
	    r.pair_rec = last_obj_rec;   /* obj_rec */

	    r.flags = GRANGE_FLAG_ISANNO;
	    r.rec = anno_ele_new(io, 0, last_obj_type, last_obj_rec, 0, txt);
	    e = (anno_ele_t *)cache_search(io, GT_AnnoEle, r.rec);
	    e = cache_rw(io, e);
	
	    bin = bin_add_range(io, &c, &r, NULL);
	    e->bin = bin->rec;
	    break;
	}

	default:
	    printf("Unsupported block type %s\n",
		   linetype2str(b->type));
	}

	if (!delay_destroy)
	    baf_block_destroy(b);

    }
    if (co)
	baf_block_destroy(co);

#if 0
	if ((nseqs & 0x3fff) == 0) {
	    static int perc = 0;
	    if (perc < 100.0 * pos / sb.st_size) {
		perc = 100.0 * pos / sb.st_size;
		printf("\r%d%%", perc);
		//HacheTableStats(io->cache, stdout);
		//HacheTableStats(((GDB **)io->dbh)[0]->gfile->idx_hash, stdout);
		{
		    static struct timeval last, curr;
		    static int first = 1;
		    static int last_contig = 0;
		    long delta;

		    gettimeofday(&curr, NULL);
		    if (first) {
			last = curr;
			first = 0;
		    }

		    delta = (curr.tv_sec - last.tv_sec) * 1000000
			+ (curr.tv_usec - last.tv_usec);
		    printf(" - %g sec %d contigs\n", delta/1000000.0,
			   ncontigs - last_contig);
		    last = curr;
		    last_contig = ncontigs;
		}
		fflush(stdout);
	    }
	}
#endif

    cache_flush(io);
    fclose(fp);

    if (pair)
	HacheTableDestroy(pair, 1);

    return 0;
}


#ifdef TEST_MAIN
int baf2aln(char *fn) {
    FILE *fp;
    baf_block *b, *co = NULL;

    if (!(fp = fopen(fn, "r"))) { 
	perror(fn);
	return -1;
    }

    while (b = baf_next_block(fp)) {
	int delay_destroy = 0;

	switch (b->type) {
	case CO:
	    if (co)
		baf_block_destroy(co);
	    co = b;
	    delay_destroy = 1;
	    break;

	case RD: {
	    int pos, len;
	    char *seq, *qual;
	    int i;

	    seq = baf_block_value(b, SQ);
	    qual = baf_block_value(b, FQ);
	    pos = atoi(baf_block_value(b, AP));
	    len = strlen(seq);
			   
	    printf("%s %s %d %d %s %s %s 0 0 0 0 %s",
		   baf_block_value(b, RD),
		   baf_block_value(co, CO),
		   pos, len,
		   baf_block_value(b, DR),
		   baf_block_value(b, QL),
		   baf_block_value(b, QR),
		   seq);
	    for (i = 0; i < len; i++) {
		printf(" %d", qual ? qual[i] - 33 : 0);
	    }
	    printf("\n");
	}
	}

	if (!delay_destroy)
	    baf_block_destroy(b);
    }

    if (co)
	baf_block_destroy(co);

    return 0;
}

int main(int argc, char **argv) {
    return baf2aln(argv[1]);
}
#endif
