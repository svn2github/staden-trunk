/* ---- Implements importing from BAF (basic alignment format). ---- */

#include <staden_config.h>

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "tg_gio.h"
#include "tg_index_common.h"
#include "hache_table.h"
#include "baf.h"
#include "zfio.h"

#define CC2(a,b) ((((unsigned char)a)<<8) | ((unsigned char)b))

enum line_type {
    XX=0, /* Blank line */
    CO=CC2('C','O'), /* Contig */
    LN=CC2('L','N'), /*    Length */
    SO=CC2('S','O'), /* DNA Source */
    ST=CC2('S','T'), /*    Source type */
    SI=CC2('S','I'), /*    Insert size mean */
    SS=CC2('S','S'), /*    Insert size standard deviation */
    SV=CC2('S','V'), /*    vector */
    RD=CC2('R','D'), /* Reading */
    SQ=CC2('S','Q'), /*    Sequence */
    FQ=CC2('F','Q'), /*    Fastq quality */
    AP=CC2('A','P'), /*    Contig position */
    QL=CC2('Q','L'), /*    Left quality clip */
    QR=CC2('Q','R'), /*    Right quality clip */
    TN=CC2('T','N'), /*    Template name */
    DR=CC2('D','R'), /*    Direction, 1=>uncomp, -1=>complemented */
    PR=CC2('P','R'), /*    Primer type */
    TR=CC2('T','R'), /*    Trace name */
    MQ=CC2('M','Q'), /*    Mapping quality */
    AL=CC2('A','L'), /* Alignment */
    AN=CC2('A','N'), /*    Annotation */
    LO=CC2('L','O'), /*    Location */
    LL=CC2('L','L'), /*    Length */
    AT=CC2('A','T'), /*    Anno. Type */
    TX=CC2('T','X'), /*    Anno. Text */

    /* Regexp versions of the above */
    ln=CC2('l','n'),
    st=CC2('s','t'),
    si=CC2('s','i'),
    ss=CC2('s','s'),
    sv=CC2('s','v'),
    pa=CC2('p','a'),
    sq=CC2('s','q'),
    fq=CC2('f','q'),
    ap=CC2('a','p'),
    ql=CC2('q','l'),
    qr=CC2('q','r'),
    tn=CC2('t','n'),
    dr=CC2('d','r'),
    pr=CC2('p','r'),
    tr=CC2('t','r'),
    mq=CC2('m','q'),
    al=CC2('a','l'),
    an=CC2('a','n'),
    lo=CC2('l','o'),
    at=CC2('a','t'),
    tx=CC2('t','x'),
};

typedef struct {
    /* Allocated memory and size */
    char *str;
    size_t len;

    /* Key, value and assignment type, eg CO=contig#1 */
    char *value;
    enum line_type type;
    int assign;  /* = or : */

    /* Order so we can reconstruct the file exactly */
    int order;
} line_t;

typedef struct {
    int type;
    HacheTable *h;
} baf_block;

/* The data attached to the hache */
typedef struct {
    char *value;
    int assign;
} baf_data;

void free_line(line_t *l);
char *linetype2str(int lt);
line_t *get_line(zfp *fp, line_t *in);
baf_block *baf_next_block(zfp *fp);
void baf_block_destroy(baf_block *b);
line_t *baf_line_for_type(baf_block *b, int type);

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
line_t *get_line(zfp *fp, line_t *in) {
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
	    if (NULL == zfgets(&l->str[pos], BLKSZ, fp)) {
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


/*
 * Relaces \n with newline and \\ with \.
 * Modifies the line in-situ as it can never grow.
 */
void unescape_line(char *txt) {
    char *cp;
    for (cp = txt; *txt; txt++) {
	if (*txt != '\\') {
	    *cp++ = *txt;
	} else {
	    if (*++txt == 'n')
		*cp++ = '\n';
	    else
		*cp++ = *txt;
	}
    }
    *cp++ = 0;
}


baf_block *baf_next_block(zfp *fp) {
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
    b->h->name = "baf-block";

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

int parse_baf(GapIO *io, char *fn, tg_args *a) {
    int nseqs = 0, nobj = 0, ntags = 0, ncontigs = 0;
    struct stat sb;
    zfp *fp;
    off_t pos;
    contig_t *c = NULL;
    HacheTable *pair = NULL;
    baf_block *b, *co = NULL;
    int last_obj_type = 0;
    int last_obj_pos = 0;
    int last_obj_rec = 0;
    int last_obj_orient = 0;
	
    printf("Loading %s...\n", fn);
    if (-1 == stat(fn, &sb) ||
	NULL == (fp = zfopen(fn, "r"))) {
	perror(fn);
	return -1;
    }

    if (a->pair_reads) {
	pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
	pair->name = "pair";
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
	    
	    create_new_contig(io, &c, contig, a->merge_contigs);

	    /* For anno */
	    last_obj_type = GT_Contig;
	    last_obj_rec = c->rec;
	    last_obj_pos = c->start + 1;
	    last_obj_orient = 0;

	    break;
	}

	case RD: {
	    seq_t seq;
	    int flags;
	    char *tname;
	    int recno;
	    int is_pair = 0;

	    /* Construct seq struct */
	    if (-1 == construct_seq_from_block(&seq, b, &tname)) {
		fprintf(stderr, "Failed to parse read block for seq %d\n",
			nseqs);
		break;
	    }

	    /* Create range, save sequence */
	    flags = GRANGE_FLAG_TYPE_SINGLE;
	    
	    if (seq.flags & SEQ_END_REV)
		flags |= GRANGE_FLAG_END_REV;
	    else
		flags |= GRANGE_FLAG_END_FWD;
	    if (seq.len < 0)
		flags |= GRANGE_FLAG_COMP1;
		
	    if (pair) is_pair = 1;
		
	    recno = save_range_sequence(io, &seq, seq.mapping_qual, pair, is_pair, tname, c, a, flags, NULL);

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

	    nseqs++;
	    
	    break;
	}

	case AN: {
	    range_t r;
	    anno_ele_t *e;
	    char *typ = baf_block_value(b, AN);
	    char *loc = baf_block_value(b, LO);
	    char *len = baf_block_value(b, LL);
	    char *txt = baf_block_value(b, TX);
	    int pos;
	    bin_index_t *bin;

	    if (txt)
		unescape_line(txt);

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

	    //r.mqual    = last_obj_type;  /* obj_type */
	    r.mqual = str2type(typ);
	    r.pair_rec = last_obj_rec;   /* obj_rec */

	    r.flags = GRANGE_FLAG_ISANNO;
	    if (GT_Seq == last_obj_type)
		r.flags |= GRANGE_FLAG_TAG_SEQ;
	    r.rec = anno_ele_new(io, 0, last_obj_type, last_obj_rec, 0,
				 str2type(typ), txt);
	    e = (anno_ele_t *)cache_search(io, GT_AnnoEle, r.rec);
	    e = cache_rw(io, e);
	
	    bin = bin_add_range(io, &c, &r, NULL, NULL);
	    e->bin = bin->rec;

	    ntags++;
	    break;
	}

	case 0:
	    /* blank line */
	    break;

	default:
	    printf("Unsupported block type '%s'\n",
		   linetype2str(b->type));
	}

	if (!delay_destroy)
	    baf_block_destroy(b);

	if ((++nobj & 0xfff) == 0) {
	    int perc = 0;

	    pos = zftello(fp);
	    perc = 100.0 * pos / sb.st_size;
	    printf("\r%d%c", perc, (nobj & 0x3fff) ? '%' : '*');
	    fflush(stdout);
	    if ((nobj & 0x3fff) == 0)
		cache_flush(io);
	}

#if 1
	if ((nobj & 0x3fff) == 0) {
	    static int perc = 0;
	    if (perc < 100.0 * pos / sb.st_size) {
		perc = 100.0 * pos / sb.st_size;
		printf("\r%d%%", perc);
		//HacheTableStats(io->cache, stdout);
		//HacheTableStats(((GDB **)io->dbh)[0]->gfile->idx_hash, stdout);
		{
		    static struct timeval last, curr;
		    static int first = 1;
		    static int last_obj = 0;
		    static int last_contigs = 0;
		    long delta;

		    gettimeofday(&curr, NULL);
		    if (first) {
			last = curr;
			first = 0;
		    }

		    delta = (curr.tv_sec - last.tv_sec) * 1000000
			+ (curr.tv_usec - last.tv_usec);
		    printf(" - %g sec %d obj (%d contigs)\n", delta/1000000.0,
			   nobj - last_obj, ncontigs - last_contigs);
		    last = curr;
		    last_obj = nobj;
		    last_contigs = ncontigs;
		}
		fflush(stdout);
	    }
	}
#endif
    }
    if (co)
	baf_block_destroy(co);

    cache_flush(io);
    zfclose(fp);

    printf("\nLoaded %12d contigs\n",     ncontigs);
    printf("       %12d sequences\n",   nseqs);
    printf("       %12d annotations\n", ntags);

    if (pair)
	HacheTableDestroy(pair, 1);

    return 0;
}
