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

/* lh3: an analogy of parse_line() */
static int parse_maqmap_aux(seq_t *s, const maqmap1_t *m, int k)
{
    int i;

    memset(s, 0, sizeof(*s));

    s->pos = (m->pos>>1);
    s->len = (m->pos&1)? -m->size : m->size;
    s->seq_tech = STECH_SOLEXA;
    s->flags = s->len < 0 ? SEQ_COMPLEMENTED : 0;
    s->left = 1;
    s->right = m->size;
    s->parent_type = 0;
    s->parent_rec = 0;
    s->other_end = 0;
    s->name_len = strlen(m->name);
    s->data = (char*)malloc(s->name_len + 2*m->size);
    s->name = s->data;
    s->seq = s->data + s->name_len;
    s->conf = s->seq+m->size;
    //s->mapping_qual = m->map_qual;
    s->mapping_qual = m->seq[MAX_READLEN-1];

    /* fill seq_t::seq && seq_t::conf */
    strcpy(s->name, m->name); /* lh3: read name */
    for (i = 0; i != m->size; ++i) {
	if (m->pos&1) {
	    /* reverse strand */
	    bit8_t c = m->seq[m->size-1-i];
	    s->seq[i] = "ACGT"[(3 - (c>>6))];
	    s->conf[i] = (c&0x3f);
	} else {
	    /* forward strand */
	    s->seq[i] = "ACGT"[m->seq[i] >> 6];
	    s->conf[i] = m->seq[i] & 0x3f;
	}
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
		 int no_tree, int pair_reads, int merge_contigs)
{
    gzFile dat_fp;
    maqmap_t *mm;
    maqmap1_t m;
    int k = 0, j = 0;
    int curr_contig = -1;
    contig_t *c = NULL;
    HacheTable *pair = NULL;

    fprintf(stderr, "-- Loading %s...\n", dat_fn);
    assert(dat_fp = gzopen(dat_fn, "r"));
    mm = maqmap_read_header(dat_fp);

    fprintf(stderr, "++ The input contains %d sequences.\n", mm->n_ref);

    if (pair_reads) {
	pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
    }

    while (gzread(dat_fp, &m, sizeof(maqmap1_t))) {
	seq_t seq;
	range_t r, *r_out;
	int recno;
	bin_index_t *bin;
	HacheItem *hi, *other_hi;

	parse_maqmap_aux(&seq, &m, k++);

	/* Create new contig if required */
	if (m.seqid != curr_contig) {
	    if (c) {
		cache_decr(io, c);
	    }
	    if (!merge_contigs ||
		NULL == (c = find_contig_by_name(io, mm->ref_name[m.seqid]))) {
		c = contig_new(io, mm->ref_name[m.seqid]);
		contig_index_update(io,
				    mm->ref_name[m.seqid],
				    strlen(mm->ref_name[m.seqid]),
				    c->rec);
	    }
	    cache_incr(io, c);
	    curr_contig = m.seqid;
	    fprintf(stderr, "++ Processing contig %d\n", m.seqid);
	}

	/* Create range */
	r.start = seq.pos;
	r.end = seq.pos + (seq.len > 0 ? seq.len : -seq.len) - 1;
	r.rec = 0;

	bin = bin_add_range(io, &c, &r, &r_out);

	/* Find pair if appropriate */
	if (pair) {
	    int new = 0;
	    HacheData hd;
	    
	    hd.i = 0;
	    hi = HacheTableAdd(pair, seq.name, seq.name_len, hd, &new);

	    /* Pair existed already */
	    if (!new) {
		seq.other_end = hi->data.i;

		/* Link other end to 'us' too */
		other_hi = hi;
		hi = NULL;
	    } else {
		other_hi = NULL;
	    }
	} else {
	    hi = NULL;
	    other_hi = NULL;
	}

	/* Save sequence */
	seq.bin = bin->rec;
	recno = sequence_new_from(io, &seq);
	if (hi)
	    hi->data.i = recno;

	if (!no_tree)
	    sequence_index_update(io, seq.name, seq.name_len, recno);
	free(seq.data);
	
	/* Link bin back to sequence too before it gets flushed */
	r_out->rec = recno;

	/* Link other end of pair back to this recno if appropriate */
	if (other_hi) {
	    seq_t *other = (seq_t *)cache_search(io, GT_Seq, other_hi->data.i);
	    sequence_set_other_end(io, &other, recno);
	    HacheTableDel(pair, other_hi, 0);
	}


	if (((j+1) & 0x3fff) == 0) {
	    fprintf(stderr, "-- %.2f%%\n", 100.0*(j+1)/mm->n_mapped_reads);

	    cache_flush(io);
	}

	++j;
    }

    cache_flush(io);

    fprintf(stderr, "-- %d reads are added.\n", k);
    gzclose(dat_fp);
    maq_delete_maqmap(mm);

    if (pair)
	HacheTableDestroy(pair, 0);

    return 0;
}
