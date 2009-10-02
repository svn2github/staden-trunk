#include <staden_config.h>
#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tg_gio.h"
#include "tg_struct.h"
#include "tg_index_common.h"
#include <sam.h>

#include <staden_config.h>
#ifdef HAVE_SAMTOOLS

#define _IOLIB 2
#include "bam.h"
#include "faidx.h"
#include "bam_maqcns.h"

typedef struct {
    bam1_t *b;
    char *seq;
    char *conf;
    int seq_len;
    int alloc_len;
    int mqual;
    int pos;
    int left;
} bam_seq_t;

typedef struct {
    GapIO *io;
    const char *fn;
    bam_seq_t *seqs;
    int nseq;
    int max_seq;
    HacheTable *pair;
    HacheTable *libs;
    contig_t *c;
    int n_inserts;
    int count;
    bam_header_t *header;
    tg_args *a;
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
    s->seq = (char *)malloc(s->alloc_len);
    s->conf = (char *)malloc(s->alloc_len);
    for (i = 0; i < p->qpos; i++) {
	s->seq[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), i)];
	s->conf[i] = bam1_qual(p->b)[i];
    }
    s->seq_len = i;
    s->pos = pos - i;
    s->b = p->b;
    s->left = i;

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
    char *s = bam1_aux(b);
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
		val->i = *(uint16_t *)tmp;
	    }
	    s+=7;
	    break;

	case 'i':
	    if (type) *type = 'i';
	    if (val) {
		char tmp[4]; /* word aligned data */
		tmp[0] = s[3]; tmp[1] = s[4]; tmp[2] = s[5]; tmp[3] = s[6];
		val->i = *(int16_t *)tmp;
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

	lrec = library_new(bio->io);
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
    
    s.pos = bs->pos;
    s.len = bs->seq_len;
    s.seq_tech = STECH_SOLEXA;
    s.flags = 0;
    s.left = bs->left+1;
    s.parent_type = 0;
    s.parent_rec = 0;
    s.name_len = strlen(bam1_qname(b));
    s.data = (char *)malloc(s.name_len + 3 + 2*s.len);
    s.name = s.data;
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

    strcpy(s.name, bam1_qname(b));
    memcpy(s.seq,  bs->seq,  s.len);
    memcpy(s.conf, bs->conf, s.len);
    
    if (bam1_strand(b)) {
	complement_seq_t(&s);
	s.flags |= SEQ_COMPLEMENTED;
    }

    /* Create the range, save the sequence */

    flags = GRANGE_FLAG_TYPE_SINGLE;
    paired = (b->core.flag & BAM_FPAIRED) ? 1 : 0;

    if (b->core.flag & BAM_FREAD1)
	s.flags |= SEQ_END_FWD;

    if (b->core.flag & BAM_FREAD2)
	s.flags |= SEQ_END_REV;

    strcpy(tname, s.name);

    if (s.name_len >= 2 && s.name[s.name_len-2] == '/') {
	tname[s.name_len-2] = 0;

	/* Check validity of name vs bit-fields */
	if ((s.name[s.name_len-1] == '1' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_FWD) ||
	    (s.name[s.name_len-1] == '2' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_REV)) {
	    fprintf(stderr, "Inconsistent read name vs flags: %s vs 0x%02x\n",
		    s.name, b->core.flag);
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

    save_range_sequence(bio->io, &s, s.mapping_qual, bio->pair, is_pair, tname, bio->c, bio->a, flags, lib);

    /* Tidy up */
    if (bs->seq)  free(bs->seq);
    if (bs->conf) free(bs->conf);
    
    if (snum+1 < bio->nseq)
	memmove(bs, bs+1, (bio->nseq - (snum+1)) * sizeof(*bs));

    bio->nseq--;
	    
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
	s->alloc_len += 100;
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
    GapIO *io = bio->io;
    int i, j, insertions = 0;
    static int last_tid = -1;

    /*
     * tid is reference id - aka contig
     * n is number of sequences at this point (pos).
     * pl is the pileup info for these sequences.
     */
    /* Create new contig if appropriate */
    if (tid != last_tid) {
	char *cname = bio->header->target_name[tid];

	/* header->target_name[b.core.tid] */
	printf("\n++Processing contig %d / %s\n", tid, cname);
	
	create_new_contig(io, &(bio->c), cname, bio->a->merge_contigs);
	
	last_tid = tid;
    }

    //    printf("Callback at pos=%d n=%d tid=%d\n", pos, n, tid);
    for (j = 0; j <= insertions; j++) {
	//printf("%d pos: %6d.%d ", tid, pos+1+bio->n_inserts, j);
	for (i = 0; i < n; i++) {
	    const bam_pileup1_t *p = &pl[i];
	    
	    if (j == 0 && p->is_head) {
		int i2;
		/* New sequence */
		//printf("^%c", p->b->core.qual+33);
		i2 = bio_new_seq(bio, p, pos+1+bio->n_inserts);
		assert(i == i2);
	    }

	    if (p->is_del) {
		/* Undercall, aka deleted base compared to ref */
		int q = (bam1_qual(p->b)[p->qpos] + 
			 (p->qpos < p->b->core.l_qseq-1
			  ? bam1_qual(p->b)[p->qpos+1]
			  : bam1_qual(p->b)[p->qpos])) / 2;
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
		    }
		    bio_extend_seq(bio, i, c, q);
		} else if (p->indel < 0 && j == 0) {
		    /*
		     * This occurs when the next base is an undercall.
		     * However THIS base should still exist, I think.
		     * It's all a little bit confusing if truth be known.
		     */
		    int c = tolower(bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)]);
		    int q = bam1_qual(p->b)[p->qpos];
		    bio_extend_seq(bio, i, c, q);
		}
	    }
	}
    }

    for (i = n-1; i >= 0; i--) {
	const bam_pileup1_t *p = &pl[i];
	
	if (p->is_tail) {
	    bio_del_seq(bio, p, i);
	}
    }
    
    bio->n_inserts += insertions;
    //    if (insertions) {
    //	printf("%d insertions (%d) at pos=%d\n", insertions, bio->n_inserts, pos);
    //    }

    return 0;
}


int parse_bam(GapIO *io, const char *fn, tg_args *a) {
    bam_io_t *bio = (bam_io_t*)calloc(1, sizeof(*bio));
    bam1_t *b;
    int ret, count = 0;
    bam_plbuf_t *plbuf;
    bamFile fp;

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

    if (a->pair_reads) {
	bio->pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
	bio->pair->name = "pair";
    } else {
	bio->pair = NULL;
    }

    fp = bam_open(fn, "r");
    assert(fp);
    bio->header = bam_header_read(fp);
    plbuf = bam_plbuf_init(bio_callback, bio);
    //bam_plbuf_set_mask(plbuf, BAM_DEF_MASK /* or BAM_FUNMAP? */);
    bam_plbuf_set_mask(plbuf, 0 /* or BAM_DEF_MASK, or BAM_FUNMAP? */);

    /*
     * Loop through reads - the bulk of the work here is done in the
     * bio_callback function.
     */
    b = (bam1_t*)calloc(1, sizeof(bam1_t));
    while ((ret = bam_read1(fp, b)) >= 0) {
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
    printf("Loaded %d sequences\n", bio->count);

    /* Tidy up */
    if (b) {
	if (b->data)
	    free(b->data);
	free(b);
    }

    if (fp)
	bam_close(fp);

    if (bio) {
	if (bio->pair)
	    HacheTableDestroy(bio->pair, 1);

	if (bio->libs)
	    HacheTableDestroy(bio->libs, 0);

	if (bio->header)
	    bam_header_destroy(bio->header);

	if (bio->seqs)
	    free(bio->seqs);
	free(bio);
    }

    return 0;
}

#endif /* HAVE_SAMTOOLS */

