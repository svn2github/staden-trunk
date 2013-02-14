#include <assert.h>
#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#include "array.h"
#include "tg_gio.h"
#include "tg_struct.h"
#include "tg_index_common.h"
#include "maqmap.h"
#include "maq.h"

/* lh3: an analogy of parse_line() */
static int parse_maqmap_aux(tg_args *a,
			    seq_t *s,
			    const int m_sz,
			    const maqmap128_t *m,
			    int k)
{
    int i, sz;

    memset(s, 0, sizeof(*s));

    s->rec = 0;
    s->pos = (m->pos>>1);
    s->len = (m->pos&1)? -m->size : m->size;
    s->seq_tech = STECH_SOLEXA;
    s->flags = s->len < 0 ? SEQ_COMPLEMENTED : 0;
    s->left = 1;
    s->right = m->size;
    s->parent_type = 0;
    s->parent_rec = 0;
    if (a->data_type & DATA_NAME) {
	s->name_len = strlen(m->name);
	s->name = (char*)malloc(s->name_len + 3+ 2*m->size);
	strcpy(s->name, m->name);
    } else {
	char *n = "";
	s->name_len = strlen(n);
	s->name = (char *)malloc(s->name_len + 3 + 2*m->size);
	strcpy(s->name, n);
    }
    s->template_name_len = 0;
    s->trace_name = s->name + s->name_len + 1;
    *s->trace_name = 0;
    s->trace_name_len = 0;
    s->alignment = s->trace_name + s->trace_name_len + 1;
    *s->alignment = 0;
    s->alignment_len = 0;
    s->seq = s->alignment + s->alignment_len+1;
    s->conf = (int8_t *) s->seq+m->size;
    //s->mapping_qual = m->map_qual;
    s->mapping_qual = m->seq[m_sz-1];

    /* fill seq_t::seq && seq_t::conf */
    sz = m->size;
    for (i = 0; i != sz; ++i) {
	if (m->pos&1) {
	    /* reverse strand */
	    bit8_t c = m->seq[m->size-1-i];
	    s->seq[i]  = (a->data_type & DATA_SEQ)  ? "TGCA"[c>>6] : 'N';
	    s->conf[i] = (a->data_type & DATA_QUAL) ? (c&0x3f)     :  0;
	} else {
	    /* forward strand */
	    bit8_t c = m->seq[i];
	    s->seq[i]  = (a->data_type & DATA_SEQ)  ? "ACGT"[c>>6] : 'N';
	    s->conf[i] = (a->data_type & DATA_QUAL) ? (c&0x3f)     :  0;
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
int parse_maqmap(GapIO *io, char *dat_fn, tg_args *a) {
    gzFile dat_fp;
    maqmap_t *mm;
    maqmap64_t m64;
    maqmap128_t m128;
    int k = 0, j = 0;
    int curr_contig = -1;
    contig_t *c = NULL;
    tg_pair_t *pair = NULL;
    HacheTable *libs = HacheTableCreate(256, HASH_DYNAMIC_SIZE);
    char tname[1024];
    int sz;
    int long_format = (a->fmt == 'M');

    fprintf(stderr, "-- Loading %s...\n", dat_fn);
    if (NULL == (dat_fp = gzopen(dat_fn, "r")))
	return -1;
    
    mm = maqmap_read_header(dat_fp);

    fprintf(stderr, "++ The input contains %d sequences.\n", mm->n_ref);

    libs->name = "libs";
    if (a->pair_reads) {
	pair = create_pair(a->pair_queue);
    }

    if (long_format) {
	sz = 128;
    } else {
	sz = maq_detect_size(dat_fp);
	if (sz == -1) {
	    fprintf(stderr, "Failed it identify whether max seq is 64 or 128\n");
	    return -1;
	}
	printf("Auto-detected maq sequence size as %d.\n", sz);
    }

    while ((sz == 64 && gzread(dat_fp, &m64, sizeof(maqmap64_t)))
	   || gzread(dat_fp, &m128, sizeof(maqmap128_t))) {
	seq_t seq;
	int flags;
	HacheItem *hi;
	int paired;
	int is_pair = 0;
	library_t *lib = NULL;
	char *LB = dat_fn;
	HacheData hd;
	int new = 0;

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

	parse_maqmap_aux(a, &seq, sz, &m128, k++);

	/* Read is unmapped, but placed along side the pair in the file */
	if (m128.flag == (PAIRFLAG_SW | PAIRFLAG_NOMATCH))
	    continue;

	/* Fetch library */
	hd.p = NULL;
	hi = HacheTableAdd(libs, (char *)LB, strlen(LB), hd, &new);
	if (new) {
	    tg_rec lrec;
	    printf("New library %s\n", LB);

	    lrec = library_new(io, LB);
	    lib = get_lib(io, lrec);
	    hi->data.p = lib;
	    cache_incr(io, lib);
	}
	lib = hi->data.p;

	/* Create new contig if required */
	if (m128.seqid != curr_contig) {
	    if (c) {
		cache_decr(io, c);
	    }
	    if (!a->merge_contigs ||
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
	    fprintf(stderr, "++ Processing contig %d\n", (int)m128.seqid);
	}


	/* calc range, save sequence */	
	flags = GRANGE_FLAG_TYPE_SINGLE;

	/* Get direction from name possibly */
	strcpy(tname, seq.name);
	if (seq.name_len >= 2 && seq.name[seq.name_len-2] == '/') {
	    tname[seq.name_len-2] = 0;
	    paired = 1;
	} else {
	    paired = 0;
	}

	if (paired)
	    flags |= (seq.flags & SEQ_END_MASK) == SEQ_END_FWD
		? GRANGE_FLAG_END_FWD
		: GRANGE_FLAG_END_REV;
	else
	    /* Guess work here. For now all <--- are rev, all ---> are fwd */
	    flags |= seq.len > 0
		? GRANGE_FLAG_END_FWD
		: GRANGE_FLAG_END_REV;
	if (seq.len < 0)
	    flags |= GRANGE_FLAG_COMP1;
	    
	if (pair && !(m128.flag & PAIRFLAG_NOMATCH)) is_pair = 1;

	save_range_sequence(io, &seq, seq.mapping_qual, pair, is_pair, tname,
			    c, a, flags, lib, NULL);

	if (((j+1) & 0xffff) == 0) {
	    static struct timeval last, curr;
	    static int first = 0;
	    long delta;
	    
	    gettimeofday(&curr, NULL);
	    if (first) {
		last = curr;
		first = 0;
	    }

	    delta = (curr.tv_sec - last.tv_sec) * 1000000
		+ (curr.tv_usec - last.tv_usec);
	    last = curr;

	    fprintf(stderr, "-- %.2f%% %g sec\n",
		    100.0*(j+1)/mm->n_mapped_reads, delta/1000000.0);

	    cache_flush(io);
	}

	++j;

	/* For benchmarking */
	//if (j == 100000)
	//break;
    }
    if (c)
	cache_decr(io, c);

    cache_flush(io);

    fprintf(stderr, "-- %d reads were added.\n", k);
    
    if (!a->fast_mode) {
    	finish_pairs(io, pair, a->link_pairs);
    }
    
    gzclose(dat_fp);
    maq_delete_maqmap(mm);

    if (pair) {
    	delete_pair(pair);
    }

    /* Decr reference count on libs */
    {
	HacheItem *hi;
	HacheIter *iter;

	if (!(iter = HacheTableIterCreate()))
	    return -1;
	
	while (NULL != (hi = HacheTableIterNext(libs, iter))) {
	    library_t *lib = hi->data.p;
	    cache_decr(io, lib);
	}

	HacheTableIterDestroy(iter);
    }
    HacheTableDestroy(libs, 0);

    return 0;
}
