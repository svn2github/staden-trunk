#include <assert.h>
#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "array.h"
#include "tg_gio.h"
#include "tg_struct.h"
#include "maqmap.h"
#include "maq.h"

typedef struct {
    int rec;
    int bin;
    int idx;
    int crec;
} pair_loc_t;

/* lh3: an analogy of parse_line() */
static int parse_maqmap_aux(seq_t *s,
			    const int m_sz,
			    const maqmap128_t *m,
			    int k)
{
    int i, sz;

    memset(s, 0, sizeof(*s));

    s->pos = (m->pos>>1);
    s->len = (m->pos&1)? -m->size : m->size;
    s->seq_tech = STECH_SOLEXA;
    s->flags = s->len < 0 ? SEQ_COMPLEMENTED : 0;
    s->left = 1;
    s->right = m->size;
    s->parent_type = 0;
    s->parent_rec = 0;
    s->name_len = strlen(m->name);
    s->data = (char*)malloc(s->name_len + 3+ 2*m->size);
    s->name = s->data;
    s->trace_name = s->name + s->name_len + 1;
    *s->trace_name = 0;
    s->trace_name_len = 0;
    s->alignment = s->trace_name + s->trace_name_len + 1;
    *s->alignment = 0;
    s->alignment_len = 0;
    s->seq = s->alignment + s->alignment_len+1;
    s->conf = s->seq+m->size;
    //s->mapping_qual = m->map_qual;
    s->mapping_qual = m->seq[m_sz-1];

    /* fill seq_t::seq && seq_t::conf */
    strcpy(s->name, m->name); /* lh3: read name */
    sz = m->size;
    for (i = 0; i != sz; ++i) {
	if (m->pos&1) {
	    /* reverse strand */
	    bit8_t c = m->seq[m->size-1-i];
	    s->seq[i] = "TGCA"[c>>6];
	    s->conf[i] = (c&0x3f);
	} else {
	    /* forward strand */
	    bit8_t c = m->seq[i];
	    s->seq[i] = "ACGT"[c>>6];
	    s->conf[i] = (c&0x3f);
	}
    }

    /* Paired end data ends with /1 or /2 */
    s->flags &= ~SEQ_END_MASK;
    if (s->name_len >= 2 && s->name[s->name_len-2] == '/') {
	if (s->name[s->name_len-1] == '1')
	    s->flags |= SEQ_END_FWD;
	else if (s->name[s->name_len-1] == '2')
	    s->flags |= SEQ_END_REV;
	else
	    fprintf(stderr, "Unknown name suffix for read %s\n",
		    s->name);
    } else {
	s->flags |= SEQ_END_FWD;
    }

    return 0;
}

/*
 * lh3: an analogy of parse_file().
 * Here no_tree, if true, requests that we do not attempt to index on sequence
 * name.
 * pair_reads, if true, attempts to identify read pairs (based on duplicate
 * sequence names) and points them at each other.
 */
int parse_maqmap(GapIO *io, int max_size, const char *dat_fn,
		 int no_tree, int pair_reads, int merge_contigs,
		 int long_format)
{
    gzFile dat_fp;
    maqmap_t *mm;
    maqmap64_t m64;
    maqmap128_t m128;
    int k = 0, j = 0;
    int curr_contig = -1;
    contig_t *c = NULL;
    HacheTable *pair = NULL;
    char tname[1024];
    int sz;

    fprintf(stderr, "-- Loading %s...\n", dat_fn);
    assert(dat_fp = gzopen(dat_fn, "r"));
    mm = maqmap_read_header(dat_fp);

    fprintf(stderr, "++ The input contains %d sequences.\n", mm->n_ref);

    if (pair_reads) {
	pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
    }

    if (long_format) {
	sz = 128;
    } else {
	sz = maq_detect_size(dat_fp);
	printf("Auto-detected maq sequence size as %d.\n", sz);
	if (sz == -1)
	    return -1;
    }

    while ((sz == 64 && gzread(dat_fp, &m64, sizeof(maqmap64_t)))
	   || gzread(dat_fp, &m128, sizeof(maqmap128_t))) {
	seq_t seq;
	range_t r, *r_out;
	int recno;
	bin_index_t *bin;
	HacheItem *hi;
	int paired;

	/* Convert m64 to m128 if needed */
	if (sz == 64) {
	    memcpy(m128.seq, m64.seq, 64);
	    m128.seq[127] = m64.seq[63];
	    m128.size     = m64.size;
	    m128.map_qual = m64.map_qual;
	    m128.i1	  = m64.i1;
	    m128.i2	  = m64.i2;
	    m128.c[0]	  = m64.c[0];
	    m128.c[1]	  = m64.c[1];
	    m128.flag	  = m64.flag;
	    m128.alt_qual = m64.alt_qual;
	    m128.seqid	  = m64.seqid;
	    m128.pos	  = m64.pos;
	    m128.dist	  = m64.dist;
	    memcpy(m128.name, m64.name, MAX_NAMELEN);
	}

	parse_maqmap_aux(&seq, sz, &m128, k++);

	/* Read is unmapped, but placed along side the pair in the file */
	if (m128.flag == (PAIRFLAG_SW | PAIRFLAG_NOMATCH)) {
	    free(seq.data);
	    continue;
	}

	/* Create new contig if required */
	if (m128.seqid != curr_contig) {
	    if (c) {
		cache_decr(io, c);
	    }
	    if (!merge_contigs ||
		m128.seqid >= mm->n_ref ||
		!mm->ref_name[m128.seqid] ||
		(NULL == (c = find_contig_by_name(io, mm->ref_name[m128.seqid])))) {
		static int cnum=1;
		char name[1024];
		
		if (m128.seqid < mm->n_ref && mm->ref_name[m128.seqid]) {
		    strcpy(name, mm->ref_name[m128.seqid]);
		} else {
		    sprintf(name, "Contig=%d", cnum++);
		}
		c = contig_new(io, name);
		contig_index_update(io, name, strlen(name), c->rec);
	    }
	    cache_incr(io, c);
	    curr_contig = m128.seqid;
	    fprintf(stderr, "++ Processing contig %d\n", m128.seqid);
	}

	/* Create range */
	r.start = seq.pos;
	r.end = seq.pos + (seq.len > 0 ? seq.len : -seq.len) - 1;
	r.rec = 0;
	r.pair_rec = 0;
	r.mqual = seq.mapping_qual;
	r.flags = GRANGE_FLAG_TYPE_SINGLE;

	/* Get direction from name possibly */
	strcpy(tname, seq.name);
	if (seq.name_len >= 2 && seq.name[seq.name_len-2] == '/') {
	    tname[seq.name_len-2] = 0;
	    paired = 1;
	} else {
	    paired = 0;
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

	/* Find pair if requested and we thing the other end worked */
	if (pair && !(m128.flag & PAIRFLAG_NOMATCH)) {
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
	    //if (new) printf("%s ADD %x\n", tname, m128.flag);

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
		//printf("%s DEL %x\n", tname, m128.flag);
		free(pl);
	    }
	}

	if (!no_tree)
	    sequence_index_update(io, seq.name, seq.name_len, recno);
	free(seq.data);
	
	/* Link bin back to sequence too before it gets flushed */
	r_out->rec = recno;

	if (((j+1) & 0xffff) == 0) {
	    fprintf(stderr, "-- %.2f%%\n", 100.0*(j+1)/mm->n_mapped_reads);

	    cache_flush(io);
	}

	++j;

	/* For benchmarking */
	//if (j == 100000)
	//break;
    }

    cache_flush(io);

    fprintf(stderr, "-- %d reads were added.\n", k);
    gzclose(dat_fp);
    maq_delete_maqmap(mm);

    if (pair) {
	HacheTableStats(pair, stdout);
	HacheTableDestroy(pair, 0);
    }

    return 0;
}
