#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tg_gio.h"
#include "tg_struct.h"

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
    bam_seq_t *seqs;
    int nseq;
    int max_seq;
    int no_tree;
    int merge_contigs;
    HacheTable *pair;
    contig_t *c;
    int n_inserts;
    int count;
} bam_io_t;

typedef struct {
    int rec;
    int bin;
    int idx;
    int crec;
} pair_loc_t;

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
	s->seq[i] = tolower(bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), i)]);
	s->conf[i] = 'X';
    }
    s->seq_len = i;
    s->pos = pos - i;
    s->b = p->b;
    s->left = i;

    return bio->nseq++;
}


/*
 * Removes a sequence from the bam_io_t struct.
 * Returns 0 on success
 *        -1 on failure
 */
int bio_del_seq(bam_io_t *bio, const bam_pileup1_t *p, int snum) {
    bam_seq_t *bs;
    bam1_t *b;
    seq_t s;
    bin_index_t *bin;
    range_t r, *r_out;
    HacheItem *hi;
    int recno, i, paired;
    GapIO *io = bio->io;
    char tname[1024];

    if (snum < 0 || snum >= bio->nseq)
	return -1;

    bio->count++;

    bs = &bio->seqs[snum];
    b = bs->b;
    /*
    printf("\nSeq %d @ %6d: '%.*s' '%.*s' => nseq=%d\n",
	   snum, bs->pos, bs->seq_len, bs->seq, bs->seq_len, bs->conf,
	   bio->nseq-1);
    */

    /* Construct a seq_t struct */
    s.right = bs->seq_len;
    for (i = p->qpos+1; i < b->core.l_qseq; i++) {
	int base = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), i)];
	int qual = '0';
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

    /* Create the range */
    r.start = s.pos;
    r.end = s.pos + (s.len > 0 ? s.len : -s.len) - 1;
    r.rec = 0;
    r.pair_rec = 0;
    r.mqual = s.mapping_qual;
    r.flags = GRANGE_FLAG_TYPE_SINGLE;

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
	r.flags |= (s.flags & SEQ_END_MASK) == SEQ_END_FWD
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;
    else
	/* Guess work here. For now all <--- are rev, all ---> are fwd */
	r.flags |= s.len > 0
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;
    if (s.len < 0)
	r.flags |= GRANGE_FLAG_COMP1;

    bin = bin_add_range(io, &bio->c, &r, &r_out);

    /* Add the sequence */
    s.bin = bin->rec;
    s.bin_index = r_out - ArrayBase(range_t, bin->rng);
    recno = sequence_new_from(io, &s);

    /* Find the read-pair if appropriate */
    if (bio->pair) {
	int new = 0;
	HacheData hd;
	pair_loc_t *pl;
	    
	/* Add data for this end */
	pl = (pair_loc_t *)malloc(sizeof(*pl));
	pl->rec  = recno;
	pl->bin  = bin->rec;
	pl->crec = bio->c->rec;
	pl->idx  = s.bin_index;
	hd.p = pl;

	hi = HacheTableAdd(bio->pair, tname, strlen(tname), hd, &new);

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

	    /* And, making an assumption, remove from hache */
	    HacheTableDel(bio->pair, hi, 1);
	    free(pl);
	}
    }

    /* Link bin back to sequence too before it gets flushed */
    r_out->rec = recno;

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
	char cname[100];
	sprintf(cname, "Contig=%d", tid);

	/* header->target_name[b.core.tid] */
	printf("\n++Processing contig %d\n", tid);
	if (bio->c)
	    cache_decr(io, bio->c);
	if (!bio->merge_contigs ||
	    (NULL == (bio->c = find_contig_by_name(io, cname)))) {
	    bio->c = contig_new(io, cname);
	    contig_index_update(io, cname, strlen(cname), bio->c->rec);
	}
	cache_incr(io, bio->c);

	bio->n_inserts = 0;

	last_tid = tid;
    }

    //    printf("Callback at pos=%d n=%d tid=%d\n", pos, n, tid);
    for (j = 0; j <= insertions; j++) {
	//printf("%d pos: %6d.%d ", tid, pos+bio->n_inserts, j);
	for (i = 0; i < n; i++) {
	    const bam_pileup1_t *p = &pl[i];
	    
	    if (j == 0 && p->is_head) {
		int i2;
		/* New sequence */
		//printf("^%c", p->b->core.qual+33);
		i2 = bio_new_seq(bio, p, pos+bio->n_inserts);
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


int parse_bam(GapIO *io, const char *fn,
	      int no_tree, int pair_reads, int merge_contigs)
{
    bam_io_t *bio = (bam_io_t*)calloc(1, sizeof(*bio));
    bam1_t *b;
    int ret, count = 0;
    bam_plbuf_t *plbuf;
    bam_header_t *header;
    bamFile fp;

    bio->io = io;
    bio->seqs = NULL;
    bio->nseq = 0;
    bio->max_seq = 0;
    bio->no_tree = no_tree;
    bio->merge_contigs = merge_contigs;
    bio->c = NULL;
    bio->count = 0;

    if (pair_reads) {
	bio->pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
    } else {
	bio->pair = NULL;
    }

    fp = bam_open(fn, "r");
    assert(fp);
    header = bam_header_read(fp);
    plbuf = bam_plbuf_init(bio_callback, bio);
    bam_plbuf_set_mask(plbuf, BAM_DEF_MASK /* or BAM_FUNMAP? */);

    b = (bam1_t*)calloc(1, sizeof(bam1_t));
    while ((ret = bam_read1(fp, b)) >= 0) {
	bam_plbuf_push(b, plbuf);
	if ((++count & 0x1fff) == 0) {
	    putchar('.');
	    fflush(stdout);
	    cache_flush(io);
	}
    }
    putchar('\n');
    bam_plbuf_push(0, plbuf);
    bam_plbuf_destroy(plbuf);

    printf("Loaded %d sequences\n", bio->count);

    cache_flush(io);
    free(b->data);
    free(b);

    return 0;
}


/* ---------------------------------------------------------------------- */
/*
 * Original implementation left here for a while, so it makes it into cvs
 * at least briefly.
 */

#if 0

#define _IOLIB 2
#include "bam.h"
#include "bam_endian.h"

typedef struct {
    int rec;
    int bin;
    int idx;
    int crec;
} pair_loc_t;


int edit_by_cigar(GapIO *io, bam1_t *b, char **eseq_p, char **equal_p,
		  int *eseq_alloc, int *cleft, int *cright) {
    bam1_core_t *c = &b->core;
    int i;
    int p = 0;
    char *eseq = *eseq_p;
    char *equal = *equal_p;
    int len = c->l_qseq;

    *cleft = 1;
    *cright = len;

    for (i = 0; i < c->n_cigar; i++) {
	int op = bam1_cigar(b)[i] & BAM_CIGAR_MASK; // operation
	int l = bam1_cigar(b)[i] >> BAM_CIGAR_SHIFT; // length

	switch (op) {
	case BAM_CMATCH:
	    p += l;
	    break;

	case BAM_CINS:
	    //printf("CINS unsupported at %d+%d len %d\n", c->pos, p, l);
	    break;

	case BAM_CDEL: {
	    int j;
	    int pqual;

	    pqual = p == 0
		? equal[p+1]
		: p == len-1
		  ? equal[p-1]
		  : equal[p-1] < equal[p+1]
		    ? equal[p-1]
		    : equal[p+1];

	    len += l;
	    *cright += len;
	    while (len >= *eseq_alloc) {
		*eseq_alloc *= 2;
		eseq = *eseq_p  = (char *)realloc(eseq,  *eseq_alloc);
		equal= *equal_p = (char *)realloc(equal, *eseq_alloc);
	    }
	    memmove(&eseq[p+1],  &eseq[p],  len-l-p);
	    memmove(&equal[p+1], &equal[p], len-l-p);
	    for (j = 0; j < l; j++, p++) {
		eseq[p] = '*';
		equal[p] = pqual;
	    }
	    break;
	}

	case BAM_CREF_SKIP:
	    break;

	case BAM_CSOFT_CLIP:
	    if (p == 0) {
		*cleft += l;
		c->pos -= l;
	    }
	    if (p + l == c->l_qseq)
		*cright -= l;
	    p += l;
	    break;

	case BAM_CHARD_CLIP:
	    if (p == 0)
		c->pos -= l;
	    break;

	case BAM_CPAD:
	    break;
	}
    }

    c->l_qseq = len;

    return 0;
}


int parse_bam(GapIO *io, const char *fn,
	      int no_tree, int pair_reads, int merge_contigs)
{
    bamFile fp;
    bam_header_t *header;
    int ret = -1;
    bam1_t b;
    contig_t *c = NULL;
    int ncontigs = 0, nseqs = 0;
    int last_tid = -999;
    HacheTable *pair = NULL;
    char tname[1024];
    int i, j = 0;
    int cleft, cright;
    char *eseq = NULL, *equal = NULL;
    int eseq_alloc = 0;
    int count = 0;

    printf("Loading %s...\n", fn);

    if (NULL == (fp = bam_open(fn, "r"))) {
	perror(fn);
	return -1;
    }
    
    if (NULL == (header = bam_header_read(fp)))
	goto fail;

    if (pair_reads) {
	pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
    }

    memset(&b, 0, sizeof(b));
    while (bam_read1(fp, &b) >= 0) {
	seq_t seq;
	range_t r, *r_out;
	int recno;
	bin_index_t *bin;
	HacheItem *hi;
	int paired;
	uint8_t *s, *q;
	int eseq_len;

	/* Reject unmapped data */
	if (b.core.flag & BAM_FUNMAP)
	    continue;

	count++;

	//	bam_view1(header, &b);

	/* Extract our sequence and quality */
	s = bam1_seq(&b);
	q = bam1_qual(&b);
	if (b.core.l_qseq >= eseq_alloc) {
	    eseq_alloc = eseq_alloc ? 2*eseq_alloc : 1024;
	    eseq = (char *)realloc(eseq, eseq_alloc);
	    equal= (char *)realloc(equal, eseq_alloc);
	}
	for (i = 0; i < b.core.l_qseq; i++) {
	    eseq[i] = bam_nt16_rev_table[bam1_seqi(s,i)];
	    equal[i] = q[i];
	}

	/* Cigar string editing */
	edit_by_cigar(io, &b, &eseq, &equal, &eseq_alloc, &cleft, &cright);

	/* Create seq struct */
	memset(&seq, 0, sizeof(seq));
	seq.pos = b.core.pos + 1;
	seq.len = b.core.l_qseq;
	seq.seq_tech = STECH_SOLEXA;
	seq.flags = 0;
	seq.left = cleft;
	seq.right = cright;
	seq.parent_type = 0;
	seq.parent_rec = 0;
	seq.name_len = strlen(bam1_qname(&b));
	seq.data = (char *)malloc(seq.name_len + 3 + 2*b.core.l_qseq);
	seq.name = seq.data;
	seq.trace_name = seq.name + seq.name_len + 1;
	*seq.trace_name = 0;
	seq.trace_name_len = 0;
	seq.alignment = seq.trace_name + seq.trace_name_len + 1;
	*seq.alignment = 0;
	seq.alignment_len = 0;
	seq.seq = seq.alignment + seq.alignment_len+1;
	seq.conf = seq.seq+b.core.l_qseq;
	seq.mapping_qual = b.core.qual;

	strcpy(seq.name, bam1_qname(&b));
	s = bam1_seq(&b);
	q = bam1_qual(&b);
	memcpy(seq.seq,  eseq,  b.core.l_qseq);
	memcpy(seq.conf, equal, b.core.l_qseq);
	/*
	printf("%36s\t%d\t%+d\t%d\t%d\t%02x\n",
	       seq.name, seq.pos, seq.len,
	       bam1_strand(&b), bam1_mstrand(&b),
	       b.core.flag);
	*/
	if (bam1_strand(&b)) {
	    complement_seq_t(&seq);
	    seq.flags |= SEQ_COMPLEMENTED;
	}

	/* Create new contig if appropriate */
	if (b.core.tid != last_tid) {
	    char cname[100];
	    sprintf(cname, "Contig=%d", b.core.tid);

	    /* header->target_name[b.core.tid] */
	    printf("++Processing contig %d\n", b.core.tid);
	    if (c)
		cache_decr(io, c);
	    if (!merge_contigs ||
		(NULL == (c = find_contig_by_name(io, cname)))) {
		c = contig_new(io, cname);
		contig_index_update(io, cname, strlen(cname), c->rec);
	    }
	    cache_incr(io, c);

	    last_tid = b.core.tid;
	}

	/* Create range */
	r.start = seq.pos;
	r.end = seq.pos + (seq.len > 0 ? seq.len : -seq.len) - 1;
	r.rec = 0;
	r.pair_rec = 0;
	r.mqual = seq.mapping_qual;
	r.flags = GRANGE_FLAG_TYPE_SINGLE;

	/* Get direction from flags */
	paired = (b.core.flag & BAM_FPAIRED) ? 1 : 0;
	if (b.core.flag & BAM_FREAD1)
	    seq.flags |= SEQ_END_FWD;
	if (b.core.flag & BAM_FREAD2)
	    seq.flags |= SEQ_END_REV;
	strcpy(tname, seq.name);
	if (seq.name_len >= 2 && seq.name[seq.name_len-2] == '/') {
	    tname[seq.name_len-2] = 0;
	    /* Check validity of name vs bit-fields */
	    if ((seq.name[seq.name_len-1] == '1' &&
		 (seq.flags & SEQ_END_MASK) != SEQ_END_FWD) ||
		(seq.name[seq.name_len-1] == '2' &&
		 (seq.flags & SEQ_END_MASK) != SEQ_END_REV)) {
		printf(stderr, "Inconsistent read name vs flags: %s vs 0x%02x\n",
		       seq.name, b.core.flag);
	    }
	}

	if (paired)
	    r.flags |= (seq.flags & SEQ_END_MASK) == SEQ_END_FWD
		? GRANGE_FLAG_END_FWD
		: GRANGE_FLAG_END_REV;
	else
	    /* Guess work here. For now all <--- are rev, all ---> are fwd */
	    r.flags |= seq.len > 0
		? GRANGE_FLAG_END_FWD
		: GRANGE_FLAG_END_REV;
	if (seq.len < 0)
	    r.flags |= GRANGE_FLAG_COMP1;

	bin = bin_add_range(io, &c, &r, &r_out);

	/* Save sequence */
	seq.bin = bin->rec;
	seq.bin_index = r_out - ArrayBase(range_t, bin->rng);
	recno = sequence_new_from(io, &seq);

	/* Find pair if appropriate */
	if (pair) {
	    int new = 0;
	    HacheData hd;
	    pair_loc_t *pl;
	    
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

	if (((j+1) & 0x1fff) == 0) {
	    //	    fprintf(stderr, "-- %.2f%%\n", 100.0*(j+1)/mm->n_mapped_reads);

	    cache_flush(io);
	}

	++j;
    }

    ret = 0;
 fail:
    if (header)
	bam_header_destroy(header);
    bam_close(fp);

    printf("Loaded %d sequences\n", count);

    return ret;
}

int parse_sam(GapIO *io, const char *hfn, const char *fn,
	      int no_tree, int pair_reads, int merge_contigs)
{
    tamFile fp;
    bam_header_t *header;
    int ret = -1;
    bam1_t b;

    if (NULL == (header = sam_header_read2(hfn)))
	goto fail;

    fp = sam_open(fn);

    memset(&b, 0, sizeof(b));
    while (sam_read1(fp, header, &b) >= 0) {
	bam_view1(header, &b);
    }

    cache_flush(io);

    ret = 0;
 fail:
    if (header)
	bam_header_destroy(header);
    sam_close(fp);
    return ret;
}

#endif

