#include <staden_config.h>
#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tg_gio.h"
#include "tg_struct.h"
#include "tg_index_common.h"
#include "sam_index.h"

#include <staden_config.h>
#ifdef HAVE_SAMTOOLS

#define _IOLIB 2
#include "bam.h"
#include "sam.h"
#include "faidx.h"
#include "bam_maqcns.h"
#include "depad_seq_tree.h"

/*
 * Uncomment this if you want sam auxillary tags to be added as tag in
 * gap5.
 */
/* #define SAM_AUX_AS_TAG */

typedef struct {
    bam1_t *b;
    char *seq;
    char *conf;
    int seq_len;
    int alloc_len;
    int mqual;
    int pos;
    int left;
    int rec;
} bam_seq_t;

typedef struct rec_list {
    int rec;
    struct rec_list *next;
} rec_list_t;

typedef struct {
    GapIO *io;
    const char *fn;
    bam_seq_t *seqs;
    int nseq;
    int max_seq;
    rec_list_t *rec_head;
    rec_list_t *rec_tail;
    HacheTable *pair;
    HacheTable *libs;
    contig_t *c;
    int n_inserts;
    int npads;
    int count;
    int skip;
    bam_header_t *header;
    tg_args *a;
    struct PAD_COUNT *tree; /* re-padding */
    int last_tid;
} bam_io_t;


int bio_extend_seq(bam_io_t *bio, int snum, char base, int conf);

/*
 * Allocates a new bam_seq_t entry in bam_io_t struct.
 * We use a single indexed array counting from 0 to N-1 representing the
 * N active sequences in a bam_pileup1_t struct. This may mean we have
 * O(D^2) complexity for depth D, so if we find this becomes a bottleneck
 * then we can replace the bam_seq_t array with an ordered linked list
 * instead.
 *
 * Returns the sequence index on success
 *         -1 on failure
 */
int bio_new_seq(bam_io_t *bio, const bam_pileup1_t *p, int pos) {
    int i;
    bam_seq_t *s;
    rec_list_t *tmp;

    //printf("New seq %d at %d %s\n", bio->nseq, pos, bam1_qname(p->b));

    if (bio->nseq == bio->max_seq) {
	bio->max_seq = bio->max_seq ? bio->max_seq * 2 : 256;
	bio->seqs = (bam_seq_t *)realloc(bio->seqs,
					 bio->max_seq * sizeof(bam_seq_t));
	if (!bio->seqs)
	    return -1;
    }
    
    s = &bio->seqs[bio->nseq];
    s->alloc_len = p->b->core.l_qseq+10;
    s->seq = (char *)malloc((int)(s->alloc_len * 1.2));
    s->conf = (char *)malloc((int)(s->alloc_len * 1.2));
    for (i = 0; i < p->qpos; i++) {
	s->seq[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), i)];
	s->conf[i] = bam1_qual(p->b)[i];
    }
    s->seq_len = i;
    s->pos = pos - i;
    s->b = p->b;
    s->left = i;

    s->rec = bio->rec_head->rec;
    tmp = bio->rec_head;
    bio->rec_head = bio->rec_head->next;
    xfree(tmp);

    return bio->nseq++;
}

typedef union {
    char  *s;
    int    i;
    float  f;
    double d;
} bam_aux_t;

/*
 * Searches for 'key' in the bam auxillary tags.
 *
 * If found, key type and value are filled out in the supplied
 * 'type' and 'val' pointers. These may be supplied as NULL if the
 * caller simply wishes to test for the presence of a key instead.
 *
 * Returns 0 if found
 *        -1 if not
 */
int bam_aux_find(bam1_t *b, char *key, char *type, bam_aux_t *val) {
    char *s = (char *)bam1_aux(b);
    int match = 0;

    while ((uint8_t *)s < b->data + b->data_len) {
	if (s[0] == key[0] && s[1] == key[1])
	    match = 1;

	switch (s[2]) {
	case 'A':
	    if (type) *type = 'A';
	    if (val) val->i = *(s+3);
	    s+=4;
	    break;

	case 'C':
	    if (type) *type = 'i';
	    if (val) val->i = *(uint8_t *)(s+3);
	    s+=4;
	    break;

	case 'c':
	    if (type) *type = 'i';
	    if (val) val->i = *(int8_t *)(s+3);
	    s+=4;
	    break;

	case 'S':
	    if (type) *type = 'i';
	    if (val) {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		val->i = *(uint16_t *)tmp;
	    }
	    s+=5;
	    break;

	case 's':
	    if (type) *type = 'i';
	    if (val) {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		val->i = *(int16_t *)tmp;
	    }
	    s+=5;
	    break;

	case 'I':
	    if (type) *type = 'i';
	    if (val) {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		val->i = *(uint32_t *)tmp;
	    }
	    s+=7;
	    break;

	case 'i':
	    if (type) *type = 'i';
	    if (val) {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		val->i = *(int32_t *)tmp;
	    }
	    s+=7;
	    break;

	case 'f':
	    if (type) *type = 'f';
	    if (val) memcpy(&val->f, s+3, 4);
	    s+=7;
	    break;

	case 'd':
	    if (type) *type = 'd';
	    if (val) memcpy(&val->d, s+3, 8);
	    s+=11;
	    break;

	case 'Z': case 'H':
	    if (type) *type = s[2];
	    s+=3;
	    if (val) val->s = s;
	    while (*s++);
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return -1;
	}

	if (match)
	    return 0;

    }

    return -1;
}
#endif /* HAVE_SAMTOOLS */

char *sam_aux_stringify(char *s, int len) {
    static char str[8192];
    char *cp = str, *s_end = s+len;
    int first = 1;

    //write(2, s, (int)(b->data + b->data_len - (uint8_t *)s));

    while (s < s_end) {
	if (first)
	    first = 0;
	else
	    *cp++ = '\t';

	switch (s[2]) {
	case 'A':
	    cp += sprintf(cp, "%c%c:A:%c", s[0], s[1], *(s+3));
	    s+=4;
	    break;

	case 'C':
	    cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint8_t *)(s+3));
	    s+=4;
	    break;

	case 'c':
	    cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int8_t *)(s+3));
	    s+=4;
	    break;

	case 'S':
	    {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 's':
	    {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 'I':
	    {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'i':
	    {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'f':
	    {
		float f;
		memcpy(&f, s+3, 4);
		cp += sprintf(cp, "%c%c:f:%f", s[0], s[1], f);
	    }
	    s+=7;
	    break;

	case 'd':
	    {
		double d;
		memcpy(&d, s+3, 8);
		cp += sprintf(cp, "%c%c:d:%f", s[0], s[1], d);
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    cp += sprintf(cp, "%c%c:%c:%s", s[0], s[1], s[2], s+3);
	    s+=3;
	    while (*s++);
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *cp = 0;

    return str;
}

#ifdef HAVE_SAMTOOLS
char *bam_aux_stringify(bam1_t *b, int no_RG) {
    static char str[8192];
    char *s = (char *)bam1_aux(b), *cp = str;
    int first = 1;
    int keep;

    no_RG = 1;

    //write(2, s, (int)(b->data + b->data_len - (uint8_t *)s));

    while ((uint8_t *)s < b->data + b->data_len) {

	keep = (no_RG && s[0] == 'R' && s[1] == 'G') ? 0 : 1;
	if (keep) {
	    if (first)
		first = 0;
	    else
		*cp++ = '\t';
	}

	switch (s[2]) {
	case 'A':
	    if (keep)
		cp += sprintf(cp, "%c%c:A:%c", s[0], s[1], *(s+3));
	    s+=4;
	    break;

	case 'C':
	    if (keep)
		cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint8_t *)(s+3));
	    s+=4;
	    break;

	case 'c':
	    if (keep)
		cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int8_t *)(s+3));
	    s+=4;
	    break;

	case 'S':
	    if (keep) {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 's':
	    if (keep) {
		char tmp[2]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4];
		cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int16_t *)tmp);
	    }
	    s+=5;
	    break;

	case 'I':
	    if (keep) {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		cp += sprintf(cp, "%c%c:i:%u", s[0], s[1], *(uint32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'i':
	    if (keep) {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		cp += sprintf(cp, "%c%c:i:%d", s[0], s[1], *(int32_t *)tmp);
	    }
	    s+=7;
	    break;

	case 'f':
	    if (keep) {
		float f;
		memcpy(&f, s+3, 4);
		cp += sprintf(cp, "%c%c:f:%f", s[0], s[1], f);
	    }
	    s+=7;
	    break;

	case 'd':
	    if (keep) {
		double d;
		memcpy(&d, s+3, 8);
		cp += sprintf(cp, "%c%c:d:%f", s[0], s[1], d);
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    if (keep)
		cp += sprintf(cp, "%c%c:%c:%s", s[0], s[1], s[2], s+3);
	    s+=3;
	    while (*s++);
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *cp = 0;

    return str;
}

/*
 * Filters out specific types from the sam aux records.
 * Returns a new aux record, also in the compact binary form.
 * 'len' holds the size of the returned string.
 *
 * Value returned is statically allocated. Do not free.
 */
char *bam_aux_filter(bam1_t *b, char **types, int ntypes, int *len) {
    static char str[8192];
    char *s = (char *)bam1_aux(b), *cp = str;
    int keep, i;

    while ((uint8_t *)s < b->data + b->data_len) {
	keep = 1;
	for (i = 0; i < ntypes; i++) {
	    if (s[0] == types[i][0] &&
		s[1] == types[i][1]) {
		keep = 0;
		break;
	    }
	}

	if (keep) {
	    *cp++ = s[0];
	    *cp++ = s[1];
	    *cp++ = s[2];
	}

	switch (s[2]) {
	case 'A':
	case 'C':
	case 'c':
	    if (keep)
		*cp++ = s[3];
	    s+=4;
	    break;

	case 'S':
	case 's':
	    if (keep) {
		*cp++ = s[3];
		*cp++ = s[4];
	    }
	    s+=5;
	    break;

	case 'I':
	case 'i':
	case 'f':
	    if (keep) {
		*cp++ = s[3];
		*cp++ = s[4];
		*cp++ = s[5];
		*cp++ = s[6];
	    }
	    s+=7;
	    break;

	case 'd':
	    if (keep) {
		*cp++ = s[3];
		*cp++ = s[4];
		*cp++ = s[5];
		*cp++ = s[6];
		*cp++ = s[7];
		*cp++ = s[8];
		*cp++ = s[9];
		*cp++ = s[10];
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    s+=3;
	    if (keep)
		while ((*cp++ = *s++));
	    else
		while (*s++);
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *len = cp - str;

    return str;
}

/*
 * Creates a new contig and updates the bam_io_t struct.
 */
void bio_new_contig(bam_io_t *bio, int tid) {
    char *cname = bio->header->target_name[tid];

    /* header->target_name[b.core.tid] */
    printf("\n++Processing contig %d / %s\n", tid, cname);
	
    create_new_contig(bio->io, &(bio->c), cname, bio->a->merge_contigs);
    bio->n_inserts = 0;
    bio->npads = 0;
    bio->skip = 0;

    if (bio->a->repad) {
	bio->tree = depad_consensus(bio->io, bio->c->rec);
	//padtree_dump(bio->tree);
    }
	
    bio->last_tid = tid;
}

/*
 * Samtools pileup won't iterate over unmapped reads. Therefore we have a
 * separate function to add these to the database - this one.
 * Although it shares much of the same code so is a candidate for merging
 * at some stage.
 */
int bio_add_unmapped(bam_io_t *bio, bam1_t *b) {
    char type;
    const char *LB;
    HacheItem *hi;
    HacheData hd;
    seq_t s;
    char tname[1024];
    library_t *lib = NULL;
    bam_aux_t val;
    int new = 0;
    char *name;
    int name_len;
    char *aux;
    int i, flags, recno;
    int paired, is_pair = 0;
    char *filter[] = {"RG"};

    bio->count++;

    /* Check if it's a new contig, create if so */
    if (b->core.tid != bio->last_tid) {
	bio_new_contig(bio, b->core.tid);
    }

    /* Fetch read-group and pretend it's a library for now */
    if (0 == bam_aux_find(b, "RG", &type, &val) && type == 'Z') {
	LB = val.s;
    } else {
	LB = bio->fn;
    }

    hd.p = NULL;
    hi = HacheTableAdd(bio->libs, (char *)LB, strlen(LB), hd, &new);
    if (new) {
	int lrec;
	printf("New library %s\n", LB);

	lrec = library_new(bio->io, (char *)LB);
	lib = get_lib(bio->io, lrec);
	hi->data.p = lib;
	cache_incr(bio->io, lib);
    }
    lib = hi->data.p;

    /* Construct a seq_t struct */
    name = bam1_qname(b);
    name_len = strlen(name);

    /*
     * Uncomment one of the two sections below if you wish to allow storing
     * sam auxillary records in gap5.
     * This is experimental and the format for these hasn't been finalised
     * yet.
     */
    aux = NULL;
    s.aux_len = 0;

    if (bio->a->sam_aux)
	aux = bam_aux_filter(b, filter, 1, &s.aux_len);

    //aux = bam_aux_stringify(b, 1);
    //s.aux_len = strlen(aux);
    
    //aux = bam1_aux(b);
    //s.aux_len = (int)(b->data + b->data_len - (uint8_t *)aux);

    s.pos = bio->npads +
	get_padded_coord(bio->tree, b->core.pos + 1 + bio->n_inserts
			 - bio->npads);
    //s.pos = b->core.pos+1;
    s.len = b->core.l_qseq;
    s.rec = 0;
    s.seq_tech = STECH_SOLEXA;
    s.flags = 0;
    s.left  = 1;
    s.right = s.len;
    s.parent_type = 0;
    s.parent_rec = 0;
    if (bio->a->data_type & DATA_NAME) {
	s.name_len = name_len;
	s.name = s.data = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, name);
    } else {
	char *n = "";
	s.name_len = 0;
	s.name = s.data = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, n);
    }
    s.trace_name = s.name + s.name_len + 1;
    *s.trace_name = 0;
    s.trace_name_len = 0;
    s.alignment = s.trace_name + s.trace_name_len + 1;
    *s.alignment = 0;
    s.alignment_len = 0;
    s.seq = s.alignment + s.alignment_len+1;
    s.conf = s.seq+s.len;
    s.mapping_qual = b->core.qual;
    s.format = SEQ_FORMAT_MAQ; /* pack bytes */
    s.anno = NULL;
    s.sam_aux = s.conf + s.len;

    for (i = 0; i < b->core.l_qseq; i++) {
	s.seq[i] = bio->a->data_type & DATA_SEQ
	    ? bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)]
	    : 'N';
	s.conf[i] = bio->a->data_type & DATA_QUAL
	    ? bam1_qual(b)[i]
	    : 0;
    }

    if (bam1_strand(b)) {
	complement_seq_t(&s);
	s.flags |= SEQ_COMPLEMENTED;
    }

    memcpy(s.sam_aux, aux, s.aux_len);

    /* Create the range, save the sequence */
    flags = GRANGE_FLAG_TYPE_SINGLE;
    paired = (b->core.flag & BAM_FPAIRED) ? 1 : 0;

    if (b->core.flag & BAM_FREAD1)
	s.flags |= SEQ_END_FWD;

    if (b->core.flag & BAM_FREAD2)
	s.flags |= SEQ_END_REV;

    strcpy(tname, name);

    if (name_len >= 2 && name[name_len-2] == '/') {
	tname[name_len-2] = 0;

	/* Check validity of name vs bit-fields */
	if ((name[name_len-1] == '1' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_FWD) ||
	    (name[name_len-1] == '2' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_REV)) {
	    fprintf(stderr, "Inconsistent read name vs flags: %s vs 0x%02x\n",
		    name, b->core.flag);
	}
    }

    if (paired)
	flags |= (s.flags & SEQ_END_MASK) == SEQ_END_FWD
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;
    else
	/* Guess work here. For now all <--- are rev, all ---> are fwd */
	flags |= bam1_strand(b)
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;

    if (bam1_strand(b)) {
	flags |= GRANGE_FLAG_COMP1;
    }

    if (bio->pair) is_pair = 1;

    flags |= GRANGE_FLAG_ISUMSEQ;

    recno = save_range_sequence(bio->io, &s, s.mapping_qual, bio->pair,
    				is_pair, tname, bio->c, bio->a, flags, lib);


#ifdef SAM_AUX_AS_TAG
    /* Make an annotation out of the sam auxillary data */
    aux = bam_aux_stringify(b, 1);
    if (aux && *aux) {
	anno_ele_t *e;
	bin_index_t *bin;
	range_t r;

	r.mqual = str2type("SAMX");
	r.start = s.pos;
	r.end = s.pos;
	r.pair_rec = recno;
	r.flags = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;
	r.rec = anno_ele_new(bio->io, 0, GT_Seq, recno, 0, r.mqual, aux);
	e = (anno_ele_t *)cache_search(bio->io, GT_AnnoEle, r.rec);
	e = cache_rw(bio->io, e);
	
	bin = bin_add_range(bio->io, &bio->c, &r, NULL, NULL);
	e->bin = bin->rec;
    }
#endif

    return 0;
}

/*
 * Removes a sequence from the bam_io_t struct.
 * This actually performs the main work of adding a sequence to the gap5
 * database.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bio_del_seq(bam_io_t *bio, const bam_pileup1_t *p, int snum) {
    bam_seq_t *bs;
    bam1_t *b;
    seq_t s;
    HacheItem *hi;
    int recno, i, paired;
    int is_pair = 0;
    int flags;
    char tname[1024];
    library_t *lib = NULL;
    bam_aux_t val;
    char type;
    const char *LB;
    HacheData hd;
    int new = 0;
    char *name;
    int name_len;
    char *aux;
    char *filter[] = {"RG"};

    if (snum < 0 || snum >= bio->nseq)
	return -1;

    bio->count++;

    bs = &bio->seqs[snum];
    b = bs->b;

    /* Fetch read-group and pretend it's a library for now */
    if (0 == bam_aux_find(b, "RG", &type, &val) && type == 'Z') {
	LB = val.s;
    } else {
	LB = bio->fn;
    }

    hd.p = NULL;
    hi = HacheTableAdd(bio->libs, (char *)LB, strlen(LB), hd, &new);
    if (new) {
	int lrec;
	printf("New library %s\n", LB);

	lrec = library_new(bio->io, (char *)LB);
	lib = get_lib(bio->io, lrec);
	hi->data.p = lib;
	cache_incr(bio->io, lib);
    }
    lib = hi->data.p;

    /*
    printf("\nSeq %d @ %6d: '%.*s' '%.*s' => nseq=%d\n",
	   snum, bs->pos, bs->seq_len, bs->seq, bs->seq_len, bs->conf,
	   bio->nseq-1);
    */

    /* Construct a seq_t struct */
    s.right = bs->seq_len;
    for (i = p->qpos+1; i < b->core.l_qseq; i++) {
	int base = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), i)];
	int qual = bam1_qual(p->b)[i];
	bio_extend_seq(bio, snum, base, qual);
    }
    
    name = bam1_qname(b);
    name_len = strlen(name);

    /*
     * Uncomment one of the two sections below if you wish to allow storing
     * sam auxillary records in gap5.
     * This is experimental and the format for these hasn't been finalised
     * yet.
     */
    aux = NULL;
    s.aux_len = 0;

    if (bio->a->sam_aux)
	aux = bam_aux_filter(b, filter, 1, &s.aux_len);
    
    //aux = bam_aux_stringify(b, 1);
    //s.aux_len = strlen(aux);

    //aux = bam1_aux(b);
    //s.aux_len = (int)(b->data + b->data_len - (uint8_t *)aux);
    
    s.rec = bs->rec;
    s.pos = bs->pos;
    s.len = bs->seq_len;
    s.seq_tech = STECH_SOLEXA;
    s.flags = 0;
    s.left = bs->left+1;
    s.parent_type = 0;
    s.parent_rec = 0;
    if (bio->a->data_type & DATA_NAME) {
	s.name_len = name_len;
	s.name = s.data = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, name);
    } else {
	char *n = "";
	s.name_len = 0;
	s.name = s.data = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, n);
    }
    s.trace_name = s.name + s.name_len + 1;
    *s.trace_name = 0;
    s.trace_name_len = 0;
    s.alignment = s.trace_name + s.trace_name_len + 1;
    *s.alignment = 0;
    s.alignment_len = 0;
    s.seq = s.alignment + s.alignment_len+1;
    s.conf = s.seq+s.len;
    s.mapping_qual = b->core.qual;
    s.format = SEQ_FORMAT_MAQ; /* pack bytes */
    s.anno = NULL;
    s.sam_aux = s.conf + s.len;

    if (bio->a->data_type & DATA_SEQ)
	memcpy(s.seq,  bs->seq,  s.len);
    else
	memset(s.seq, 'N', s.len);

    if (bio->a->data_type & DATA_QUAL)
	memcpy(s.conf, bs->conf, s.len);
    else
	memset(s.conf, 0, s.len);
    
    if (bam1_strand(b)) {
	complement_seq_t(&s);
	s.flags |= SEQ_COMPLEMENTED;
    }

    memcpy(s.sam_aux, aux, s.aux_len);

    /* Create the range, save the sequence */
    flags = GRANGE_FLAG_TYPE_SINGLE;
    paired = (b->core.flag & BAM_FPAIRED) ? 1 : 0;

    if (b->core.flag & BAM_FREAD1)
	s.flags |= SEQ_END_FWD;

    if (b->core.flag & BAM_FREAD2)
	s.flags |= SEQ_END_REV;

    strcpy(tname, name);

    if (name_len >= 2 && name[name_len-2] == '/') {
	tname[name_len-2] = 0;

	/* Check validity of name vs bit-fields */
	if ((name[name_len-1] == '1' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_FWD) ||
	    (name[name_len-1] == '2' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_REV)) {
	    fprintf(stderr, "Inconsistent read name vs flags: %s vs 0x%02x\n",
		    name, b->core.flag);
	}
    }

    if (paired)
	flags |= (s.flags & SEQ_END_MASK) == SEQ_END_FWD
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;
    else
	/* Guess work here. For now all <--- are rev, all ---> are fwd */
	flags |= bam1_strand(b)
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;

    if (bam1_strand(b)) {
	flags |= GRANGE_FLAG_COMP1;
    }

    if (bio->pair) is_pair = 1;

    recno = save_range_sequence(bio->io, &s, s.mapping_qual, bio->pair,
    				is_pair, tname, bio->c, bio->a, flags, lib);


#ifdef SAM_AUX_AS_TAG
    /* Make an annotation out of the sam auxillary data */
    aux = bam_aux_stringify(b, 1);
    if (aux && *aux) {
	anno_ele_t *e;
	bin_index_t *bin;
	range_t r;

	r.mqual = str2type("SAMX");
	r.start = s.pos;
	r.end = s.pos;
	r.pair_rec = recno;
	r.flags = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;
	r.rec = anno_ele_new(bio->io, 0, GT_Seq, recno, 0, r.mqual, aux);
	e = (anno_ele_t *)cache_search(bio->io, GT_AnnoEle, r.rec);
	e = cache_rw(bio->io, e);
	
	bin = bin_add_range(bio->io, &bio->c, &r, NULL, NULL);
	e->bin = bin->rec;
    }
#endif

    /* Tidy up */
    if (bs->seq)  free(bs->seq);
    if (bs->conf) free(bs->conf);
    
    return 0;
}

/*
 * Adds a base to a sequence in the bio_io_t struct.
 * Returns 0 on success
 *        -1 on failure
 */
int bio_extend_seq(bam_io_t *bio, int snum, char base, int conf) {
    bam_seq_t *s;

    if (snum < 0 || snum >= bio->nseq)
	return -1;

    s = &bio->seqs[snum];

    /* Extend if appropriate */
    if (s->seq_len >= s->alloc_len) {
	s->alloc_len = (s->alloc_len + 100)*1.5;
	if (NULL == (s->seq  = (char *)realloc(s->seq,  s->alloc_len)))
	    return -1;
	if (NULL == (s->conf = (char *)realloc(s->conf, s->alloc_len)))
	    return -1;
    }
    
    /* And append */
    s->seq [s->seq_len] = base;
    s->conf[s->seq_len] = conf;
    s->seq_len++;

    return 0;
}

/*
 * Called once for every position (pos) in every contig (tid), with an 
 * array of n sequences (pl).
 *
 * We produce our own stack of reads that we append to, spotting new
 * sequences (p->is_head) and terminating sequences (p->is_tail) flushing
 * out our data to disk only on sequence removal. It's possible we could
 * hijack the pl->b bam1_t struct with our own pointer, putting it back to
 * the the bam1_t only on removal.
 */
int bio_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl,
		void *data) {
    bam_io_t *bio = (bam_io_t *)data;
    int i, j, k, insertions = 0;
    int np;

    //printf("\nCallback at pos=%d n=%d tid=%d\n", pos, n, tid);

    /*
     * tid is reference id - aka contig
     * n is number of sequences at this point (pos).
     * pl is the pileup info for these sequences.
     */
    /* Create new contig if appropriate */
    if (tid != bio->last_tid) {
	bio_new_contig(bio, tid);
    }
    
    np = 0;
    if (bio->a->repad) {
	if ((np=padtree_pad_at(bio->tree, pos+1+bio->n_inserts-bio->npads))) {
	    /* Add pads to match existing consensus gaps */
	    for (j = 0; j < np; j++) {
		if (bio->skip) {
		    //printf("Skipping import of pads\n");
		    bio->skip--;
		    continue;
		}

		//printf("Import pads from existing consensus at %d\n",
		//       pos+1+bio->n_inserts-bio->npads);
		for (i = 0; i < n; i++) {
		    bio_extend_seq(bio, i, '*', 0);
		}
	    }
	}
    }

    for (j = 0; j <= insertions; j++) {
	//printf("%d pos: %6d.%d ", tid, pos+1+bio->n_inserts, j);

//	printf("%5d: ", pos+1+bio->n_inserts);
//	for (i = 0; i < n; i++) {
//	    const bam_pileup1_t *p = &pl[i];
//	    printf("%d%d ", p->indel, j);
//	}
//	printf("\n");

	if (j && bio->a->repad /*&& j > np*/) {
	    /* FIXME: And not already in the originally padded version */
	    np=padtree_pad_at(bio->tree, pos+2);

	    //printf("Insert at %d: j=%d np=%d\n", pos+2+bio->n_inserts, j, np);
	    //	    contig_insert_base(bio->io, &bio->c, pos+2+bio->n_inserts,
	    //			       '*', 0);
	    if (j > np) {
		contig_insert_base(bio->io, &bio->c,
				   get_padded_coord(bio->tree, pos+2)
				   +bio->n_inserts,
				   '*', 0);
		bio->npads++;
	    } else {
		bio->n_inserts--;
		bio->skip++;
	    }
	}

	for (i = 0; i < n; i++) {
	    const bam_pileup1_t *p = &pl[i];
	    
	    if (j == 0 && p->is_head) {
		int i2, ppos;
		/* New sequence */
		//printf("^%c", p->b->core.qual+33);
		ppos = bio->npads
		     + get_padded_coord(bio->tree,
					pos + 1 + bio->n_inserts - bio->npads);
		i2 = bio_new_seq(bio, p, ppos);

		/*
		 * The following fails if we have a previous sequence that
		 * ended on a deletion. Eg CIGAR string 36M1D. This causes
		 * samtools tview to break too.
		 */
		assert(i == i2);
	    }

	    if (p->is_del) {
		/* Undercall, aka deleted base compared to ref */
		int q = (bam1_qual(p->b)[p->qpos] + 
			 (p->qpos < p->b->core.l_qseq-1
			  ? bam1_qual(p->b)[p->qpos+1]
			  : bam1_qual(p->b)[p->qpos])) / 2;
		//printf("Undercall in seq %d\n", i);
		bio_extend_seq(bio, i, '*', q);
	    } else {
		if (p->indel >= 0) {
		    int c, q;

		    /* Overcall - edit all other seqs */
		    if (insertions < p->indel)
			insertions = p->indel;

		    if (j <= p->indel) {
			c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b),
							 p->qpos + j)];
			q = bam1_qual(p->b)[p->qpos+j];
		    } else {
			c = '*';
			q = (bam1_qual(p->b)[p->qpos] + 
			     (p->qpos < p->b->core.l_qseq-1
			      ? bam1_qual(p->b)[p->qpos+1]
			      : bam1_qual(p->b)[p->qpos])) / 2;
			if (p->is_tail)
			    continue;
		    }
		    bio_extend_seq(bio, i, c, q);
		} else if (p->indel < 0 /*&& j == 0*/) {
		    /*
		     * This occurs when the next base is an undercall.
		     * However THIS base should still exist, I think.
		     * It's all a little bit confusing if truth be known.
		     */
		    int c = j == 0
		      ? bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)]
		      : '*';
		    int q = bam1_qual(p->b)[p->qpos];
		    //printf("%s extend seq %d base %c\n",
		    //       bam1_qname(p->b), i, c);
		    bio_extend_seq(bio, i, c, q);
		}
	    }
	}
    }

    for (i = j = k = 0; i < n; i++, j++, k++) {
	const bam_pileup1_t *p = &pl[j];
	
	bio->seqs[i] = bio->seqs[j];
	if (p->is_tail) {
	    bio_del_seq(bio, p, i);
	    i--;
	    n--;
	    bio->nseq--;
	}
    }

    bio->n_inserts += insertions;
    //    if (insertions) {
    //	printf("%d insertions (%d) at pos=%d\n", insertions, bio->n_inserts, pos);
    //    }

    return 0;
}

int parse_sam_or_bam(GapIO *io, const char *fn, tg_args *a, char *mode) {
    bam_io_t *bio = (bam_io_t*)calloc(1, sizeof(*bio));
    bam1_t *b;
    int count = 0;
    bam_plbuf_t *plbuf;
    samfile_t *fp;
    rec_list_t *tmp;

    /* for pair data */
    open_tmp_file();

    /* Setup bam_io_t object and create our pileup interface */
    bio->io = io;
    bio->seqs = NULL;
    bio->nseq = 0;
    bio->max_seq = 0;
    bio->a = a;
    bio->c = NULL;
    bio->count = 0;
    bio->fn = fn;
    bio->libs = HacheTableCreate(256, HASH_DYNAMIC_SIZE);
    bio->libs->name = "libs";
    bio->last_tid = -1;
    bio->tree = NULL;
    bio->rec_head = bio->rec_tail = NULL;

    if (a->pair_reads) {
	bio->pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
	bio->pair->name = "pair";
    } else {
	bio->pair = NULL;
    }

    fp = samopen(fn, mode, NULL);
    assert(fp);
    bio->header = fp->header;
    plbuf = bam_plbuf_init(bio_callback, bio);
    //bam_plbuf_set_mask(plbuf, BAM_DEF_MASK /* or BAM_FUNMAP? */);
    bam_plbuf_set_mask(plbuf, 0 /* or BAM_DEF_MASK, or BAM_FUNMAP? */);

    /*
     * Loop through reads - the bulk of the work here is done in the
     * bio_callback function.
     */
    b = (bam1_t*)calloc(1, sizeof(bam1_t));
    while (samread(fp, b) >= 0) {
	if (a->store_unmapped && b->core.flag & BAM_FUNMAP) {
	    bio_add_unmapped(bio, b);
	    if ((++count & 0xffff) == 0) {
	        putchar('.');
	        fflush(stdout);
	        cache_flush(io);
	    }
	    continue;
	} else if (b->core.flag & BAM_FUNMAP)
	    continue;

	//printf("push %s\n", bam1_qname(b));
	/*
	 * Allocate a record number now, so they're used in the same order
	 * as the input data regardless of whether the read is mapped or
	 * unmapped.
	 */
	if (NULL == (tmp = xmalloc(sizeof(*tmp))))
	    return -1;

	tmp->next = NULL;
	if (bio->rec_head) {
	    bio->rec_tail->next = tmp;
	    bio->rec_tail = tmp;
	} else {
	    bio->rec_head = bio->rec_tail = tmp;
	}
	if (bio->a->data_type == DATA_BLANK) {
	    static int fake_recno = 1;
	    tmp->rec = fake_recno++;
	} else {
	    tmp->rec = sequence_new_from(bio->io, NULL);
	}

	bam_plbuf_push(b, plbuf);
	if ((++count & 0xffff) == 0) {
	    putchar('.');
	    fflush(stdout);
	    cache_flush(io);
	}
    }
    putchar('\n');
    bam_plbuf_push(0, plbuf);
    bam_plbuf_destroy(plbuf);

    cache_flush(io);
    vmessage("Loaded %d of %d sequences\n", bio->count, count);

    if (bio->pair && !a->fast_mode) {    
	sort_pair_file();
	
	complete_pairs(io);
	
	close_tmp_file();
    }
 
    /* Tidy up */
    if (b) {
	if (b->data)
	    free(b->data);
	free(b);
    }

    if (fp)
	samclose(fp);

    if (bio) {
	if (bio->pair)
	    HacheTableDestroy(bio->pair, 1);

	if (bio->libs)
	    HacheTableDestroy(bio->libs, 0);

	if (bio->seqs)
	    free(bio->seqs);

	if (bio->tree)
	    depad_seq_tree_free(bio->tree);

	while (bio->rec_head) {
	    tmp = bio->rec_head->next;
	    free(bio->rec_head);
	    bio->rec_head = tmp;
	}

	free(bio);
    }

    return 0;
}

int parse_bam(GapIO *io, const char *fn, tg_args *a) {
    return parse_sam_or_bam(io, fn, a, "rb");
}

int parse_sam(GapIO *io, const char *fn, tg_args *a) {
    return parse_sam_or_bam(io, fn, a, "r");
}

#endif /* HAVE_SAMTOOLS */

