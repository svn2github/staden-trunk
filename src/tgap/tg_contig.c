#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "tg_gio.h"
#include "misc.h"

/*
 * Sets the contig start position
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_set_start(GapIO *io, contig_t **c, int value) {
    contig_t *n;
    if (!(n = cache_rw(io, *c)))
	return -1;

    n->start = value;
    *c = n;

    return 0;
}

/*
 * Sets the contig end position
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_set_end(GapIO *io, contig_t **c, int value) {
    contig_t *n;
    if (!(n = cache_rw(io, *c)))
	return -1;

    n->end = value;
    *c = n;

    return 0;
}

/*
 * Sets the contig bin number
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_set_bin(GapIO *io, contig_t **c, int value) {
    contig_t *n;
    if (!(n = cache_rw(io, *c)))
	return -1;

    n->bin = value;
    *c = n;

    return 0;
}

/*
 * Sets a contig name.
 *
 * FIXME: renaming a contig means we need to update the b+tree index
 * on contig name too. TO DO
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_set_name(GapIO *io, contig_t **c, char *name) {
    contig_t *n;

    if (!(n = cache_rw(io, *c)))
	return -1;

    if (NULL == (n = cache_item_resize(n, sizeof(*n) + strlen(name)+1)))
	return -1;

    n->name   = (char *)(&n->name+1);
    *c = n;
    strcpy(n->name, name);

    return 0;
}

/*
 * Returns the offset to apply to the root bin of the contig.
 * This takes into account whether it have been complemented or not.
 */
int contig_offset(GapIO *io, contig_t **c) {
    bin_index_t *bin = (bin_index_t *)cache_search(io, GT_Bin,
						   contig_get_bin(c));
    int offset;

    //    if (bin->flags & BIN_COMPLEMENTED) {
    //	offset = ((*c)->end + (*c)->start) - (bin->pos + bin->size) + 1;
    //    } else {
	offset = bin->pos;
	//    }

    return offset;
}

/*
 * Inserts 'len' bases at position 'pos' in the contig.
 * This edits all sequences stacked up above pos and also updates the bin
 * structure accordingly.
 *
 * Returns 0 for success
 *        -1 for failure
 */
static int contig_insert_base2(GapIO *io, int bnum,
			       int pos, int offset, char base, int conf) {
    int i;
    bin_index_t *bin;

    bin = get_bin(io, bnum);
    cache_incr(io, bin);

    if (!(bin = cache_rw(io, bin)))
	return -1;

    /* Normalise pos for complemented bins */
    if (bin->flags & BIN_COMPLEMENTED) {
	pos = offset + bin->size - pos;
    } else {
	pos -= offset;
    }

    /* Perform the insert into objects first */
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);

	if (MAX(r->start, r->end) >= pos && MIN(r->start, r->end) < pos) {
	    /*
	     * FIXME: This should be a callback function so we can
	     * perform insertion in any bin system - to sequences in contigs,
	     * contigs in super-contigs, or tags in a sequence/consensus.
	     */

	    /* Insert */
	    seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	    sequence_insert_base(io, &s, pos - MIN(r->start, r->end),
				 base, conf, 0);
	    
	    r->end++;
	    bin->flags |= BIN_RANGE_UPDATED;

	} else if (MIN(r->start, r->end) >= pos) {
	    /* Move */
	    r->start++;
	    r->end++;
	    bin->flags |= BIN_RANGE_UPDATED;
	}
    }

    /* Adjust the bin dimensions */
    bin->size++;
    if (bin->start_used != bin->end_used) {
	if (pos <= bin->start_used)
	    bin->start_used++;
	if (pos <= bin->end_used)
	    bin->end_used++;
    }
    bin->flags |= BIN_BIN_UPDATED;

    /* Recurse */
    if (bin->child[0] || bin->child[1]) {
	bin_index_t *ch = NULL;

	/* Find the correct child node */
       	for (i = 0; i < 2; i++) {
	    if (!bin->child[i])
		continue;

	    ch = get_bin(io, bin->child[i]);

	    if (pos >= MIN(ch->pos, ch->pos + ch->size-1) &&
		pos <= MAX(ch->pos, ch->pos + ch->size-1)) {
		contig_insert_base2(io, bin->child[i], pos,
				    MIN(ch->pos, ch->pos + ch->size-1),
				    base, conf);
	    /* All children to the right of this one need pos updating too */
	    } else if (pos < MIN(ch->pos, ch->pos + ch->size-1)) {
		ch = get_bin(io, bin->child[i]);
		if (!(ch = cache_rw(io, ch)))
		    return -1;

		ch->pos++;
		ch->flags |= BIN_BIN_UPDATED;
	    }
	}
    }

    cache_decr(io, bin);

    return 0;
}

int contig_insert_base(GapIO *io, contig_t **c, int pos, char base, int conf) {
    contig_t *n;
    if (!(n = cache_rw(io, *c)))
	return -1;
    *c = n;

    contig_insert_base2(io, contig_get_bin(c), pos, contig_offset(io, c),
			base, conf);

    contig_set_end(io, c, contig_get_end(c)+1);

    return 0;
}

/*
 * FIXME: This cannot be finished until tg_index is rewritten to use the
 * cache instead of array entries.
 *
 * We should always have child and parent values as record numbers and
 * allocate bin recs upon creation. Currently we cannot save parent info
 * as we allocate the parent only after the children have been written.
 */
static int bin_delete(GapIO *io, bin_index_t *bin) {
    bin_index_t *parent;

    if (!bin->parent)
	return 0;

    parent = get_bin(io, bin->parent);
    cache_incr(io, parent);

    if (parent->child[0] == bin->bin_id)
	parent->child[0] = 0;

    if (parent->child[1] == bin->bin_id)
	parent->child[1] = 0;

    parent->flags |= BIN_BIN_UPDATED;

    cache_decr(io, parent);

    return 0;
}

static int contig_delete_base2(GapIO *io, int bnum,
			       int pos, int offset) {
    int i;
    bin_index_t *bin;

    bin = get_bin(io, bnum);
    cache_incr(io, bin);

    if (!(bin = cache_rw(io, bin)))
	return -1;

    /* Normalise pos for complemented bins */
    if (bin->flags & BIN_COMPLEMENTED) {
	pos = offset + bin->size - pos - 1;
    } else {
	pos -= offset;
    }

    /* Perform the deletion into objects first */
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);

	if (MAX(r->start, r->end) >= pos && MIN(r->start, r->end) <= pos) {
	    /*
	     * FIXME: This should be a callback function so we can
	     * perform deletion in any bin system - to sequences in contigs,
	     * contigs in super-contigs, or tags in a sequence/consensus.
	     */
	    seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);

	    if (MIN(r->start, r->end) == MAX(r->start, r->end)) {
		/* Remove object entirely */
		fprintf(stderr, "Delete sequence %s\n", s->name);
		r->start = r->end = -1;
	    } else {
		/* Delete */
		sequence_delete_base(io, &s, pos - MIN(r->start, r->end), 0);
	    
		r->end--;
		bin->flags |= BIN_RANGE_UPDATED;
	    }

	} else if (MIN(r->start, r->end) >= pos) {
	    /* Move */
	    r->start--;
	    r->end--;
	    bin->flags |= BIN_RANGE_UPDATED;
	}
    }

    /* Adjust the bin dimensions */
    if (--bin->size <= 0) {
	/* Remove object entirely */
	fprintf(stderr, "Delete bin bin-%d\n", bin->bin_id);
	bin->size = 0;
	bin_delete(io, bin);
    } else {
	if (bin->start_used != bin->end_used) {
	    if (pos < bin->start_used)
		bin->start_used--;
	    if (pos <= bin->end_used)
		bin->end_used--;
	}
	bin->flags |= BIN_BIN_UPDATED;
    }

    /* Recurse */
    if (bin->child[0] || bin->child[1]) {
	bin_index_t *ch = NULL;

	/* Find the correct child node */
       	for (i = 0; i < 2; i++) {
	    if (!bin->child[i])
		continue;
	    ch = get_bin(io, bin->child[i]);

	    if (pos >= MIN(ch->pos, ch->pos + ch->size-1) &&
		pos <= MAX(ch->pos, ch->pos + ch->size-1)) {
		contig_delete_base2(io, bin->child[i], pos,
				    MIN(ch->pos, ch->pos + ch->size-1));
		//break;
	    }
	}

	/* All children to the right of this one need pos updating too */
	if (i == 0 && bin->child[1-i]) {
	    ch = get_bin(io, bin->child[1-i]);
	    if (!(ch = cache_rw(io, ch)))
		return -1;

	    if (--ch->pos < 0)
		ch->pos = 0;
	    ch->flags |= BIN_BIN_UPDATED;
	}
    }

    cache_decr(io, bin);

    return 0;
}

int contig_delete_base(GapIO *io, contig_t **c, int pos) {
    contig_t *n;

    if (!(n = cache_rw(io, *c)))
	return -1;
    *c = n;

    contig_delete_base2(io, contig_get_bin(c), pos, contig_offset(io, c));

    contig_set_end(io, c, contig_get_end(c)-1);

    return 0;
}

contig_t *contig_new(GapIO *io, char *name) {
    int rec;
    contig_t *c;

    /* Allocate our contig */
    rec = io->iface->contig.create(io->dbh, NULL);

    /* Initialise it */
    c = (contig_t *)cache_search(io, GT_Contig, rec);
    c = cache_rw(io, c);
    c->bin = bin_new(io, 0, MIN_BIN_SIZE, rec, GT_Contig);
    if (name)
	contig_set_name(io, &c, name);
    else
	c->name = NULL;

    /* Add it to the contig order too */
    io->contig_order = cache_rw(io, io->contig_order);
    io->db = cache_rw(io, io->db);
    ARR(GCardinal, io->contig_order, io->db->Ncontigs++) = rec;

    return c;
}

/*
 * Looks for a contig by name and if found it returns that contig object.
 *
 * Returns contig_t ptr on success
 *         NULL on failure (not found)
 */
contig_t *find_contig_by_name(GapIO *io, char *name) {
    int rec = io->iface->contig.index_query(io->dbh, name);
    return rec > 0 ? (contig_t *)cache_search(io, GT_Contig, rec) : NULL;
}

GRec contig_index_query(GapIO *io, char *name) {
    return io->iface->contig.index_query(io->dbh, name);
}

int contig_index_update(GapIO *io, char *name, int name_len, GRec rec) {
    char n2[1024];
    GRec r;
    sprintf(n2, "%.*s", name_len, name);

    r = io->iface->contig.index_add(io->dbh, n2, rec);
    if (r == -1)
	return -1;

    if (r != io->db->contig_name_index) {
	io->db = cache_rw(io, io->db);
	io->db->contig_name_index = r;
    }

    return 0;
}

/*
 * Compute read-pairs for rangec_t structs, where viable.
 * This means when we have both ends visible within our range, we link them
 * together. When we don't we leave the pair information as unknown.
 */
static void pair_rangec(GapIO *io, rangec_t *r, int count) {
    int i;
    HacheTable *h;
    HacheIter *iter;
    HacheItem *hi;

    /* Build a hash on record number */
    h = HacheTableCreate(count, HASH_DYNAMIC_SIZE);
    for (i = 0; i < count; i++) {
	HacheData hd;
	hd.i = i;
	HacheTableAdd(h, (char *)&r[i].rec, sizeof(r[i].rec), hd, NULL);
    }

    /* Iterate through hash linking the pairs together */
    iter = HacheTableIterCreate();
    while (hi = HacheTableIterNext(h, iter)) {
	HacheItem *pair;
	int i = hi->data.i;
	int p;
	assert(i < count && i >= 0);

	pair = HacheTableSearch(h, (char *)&r[i].pair_rec, sizeof(r[i].rec));
	if (pair) {
	    p = pair->data.i;
	    assert(p < count && p >= 0);

	    r[i].pair_ind = p;
	    r[p].pair_ind = i;
	    r[i].pair_start = r[p].start;
	    r[i].pair_end   = r[p].end;
	    r[i].pair_mqual = r[p].mqual;
	    if (r[p].flags &  GRANGE_FLAG_COMP1)
		r[i].flags |= GRANGE_FLAG_COMP2;
	    r[i].flags |= GRANGE_FLAG_CONTIG;
	    r[p].flags |= GRANGE_FLAG_CONTIG;
	    r[i].flags &= ~GRANGE_FLAG_PEND_MASK;
	    if ((r[p].flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD)
		r[i].flags |= GRANGE_FLAG_PEND_FWD;
	    else
		r[i].flags |= GRANGE_FLAG_PEND_REV;
	} else {
	    r[i].pair_ind = -1;
	}
    }

    /* Tidy up */
    HacheTableIterDestroy(iter);
    HacheTableDestroy(h, 0);
}


/* Sort comparison function for range_t; sort by ascending position */
static int sort_range_by_x(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    return r1->start - r2->start;
}


/* Sort comparison function for range_t; sort by ascending position */
static int sort_range_by_y(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    if (r1->y != r2->y)
        return r1->y - r2->y;
    return r1->start - r2->start;
}


/*
 * Determines the Y position for a sequence to display in the contig editor.
 * Returns 0 on success.
 *        -1 on failure
 */
static int compute_ypos(rangec_t *r, int nr, int job) {
    int i;

    if (job & CSIR_ALLOCATE_Y_SINGLE) {
	for (i = 0; i < nr; i++) {
	    r[i].y = i;
	}
    } else {
	/* CSIR_ALLOCATE_Y_MULTIPLE */
	int *yend, j;
	
	if (NULL == (yend = (int *)malloc((nr+1) * sizeof(int))))
	    return -1;
	for (i = 0; i <= nr; i++) {
	    yend[i] = -INT_MAX;
	}

	/*
	 * O(N^2) in depth. Possibly want to keep an ordered list of the
	 * line positions so we can find the next free line rapidly.
	 */
	for (i = 0; i < nr; i++) {
	    for (j = 0; j < nr; j++) {
		if (yend[j] < r[i].start && yend[j+1] < r[i].start) {
		    r[i].y = j;
		    yend[j] = r[i].end + 2;
		    break;
		}
	    }
	}
	free(yend);
    }

    return 0;
}


#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))
static int contig_seqs_in_range2(GapIO *io, int bin_num,
				 int start, int end, int offset,
				 rangec_t **results, int *alloc, int used,
				 int complement) {
    int count = used;
    int i, n, f_a, f_b;
    bin_index_t *bin = get_bin(io, bin_num);
    range_t *l;

    cache_incr(io, bin);
    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }

    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }
    
    if (!(end < NMIN(bin->start_used, bin->end_used) ||
	  start > NMAX(bin->start_used, bin->end_used))
	&& bin->rng) {
	for (n = 0; n < ArrayMax(bin->rng); n++) {
	    l = arrp(range_t, bin->rng, n);
	    if (l->start == -1 && l->end == -1)
		continue;
	    if (NMAX(l->start, l->end) >= start
		&& NMIN(l->start, l->end) <= end) {
		int st, en;
		if (count >= *alloc) {
		    *alloc = *alloc ? *alloc * 2 : 16;
		    *results = (rangec_t *)realloc(*results,
						   *alloc * sizeof(rangec_t));
		}
		(*results)[count].rec   = l->rec;
		st = NORM(l->start);
		en = NORM(l->end);
		if (st <= en) {
		    (*results)[count].start = st;
		    (*results)[count].end   = en;
		} else {
		    (*results)[count].start = en;
		    (*results)[count].end   = st;
		}
		(*results)[count].comp  = complement;
		(*results)[count].mqual = l->mqual;
		(*results)[count].pair_rec   = l->pair_rec;
		(*results)[count].pair_start = 0;
		(*results)[count].pair_end   = 0;
		(*results)[count].pair_mqual = 0;

		(*results)[count].flags = l->flags;
		(*results)[count].y = 0;
		count++;
	    }
	}
    }

    /* Recurse down bins */
    for (i = 0; i < 2 > 0; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;
	ch = get_bin(io, bin->child[i]);
	if (end   >= NMIN(ch->pos, ch->pos + ch->size-1) &&
	    start <= NMAX(ch->pos, ch->pos + ch->size-1)) {
	    count = contig_seqs_in_range2(io, bin->child[i], start, end,
					  NMIN(ch->pos, ch->pos + ch->size-1),
					  results, alloc, count,
					  complement);
	}
    }

    cache_decr(io, bin);
    return count;
}

rangec_t *contig_seqs_in_range(GapIO *io, contig_t **c, int start, int end,
			       int job, int *count) {
    rangec_t *r = NULL;
    int alloc = 0;

    *count = contig_seqs_in_range2(io, contig_get_bin(c), start, end,
    				   contig_offset(io, c), &r, &alloc, 0, 0);
    //*count = contig_seqs_in_range2(io, contig_get_bin(c), start, end,
    //                               0, &r, &alloc, 0, 0);

    if (job & CSIR_PAIR)
	pair_rangec(io, r, *count);

    if (job & (CSIR_SORT_BY_X | CSIR_SORT_BY_Y))
	qsort(r, *count, sizeof(*r), sort_range_by_x);

    if (job & CSIR_ALLOCATE_Y)
	compute_ypos(r, *count, job & CSIR_ALLOCATE_Y);

    if (job & CSIR_SORT_BY_Y)
	qsort(r, *count, sizeof(*r), sort_range_by_y);

    return r;
}

static int contig_bins_in_range2(GapIO *io, int bin_num,
				 int start, int end, int offset,
				 rangec_t **results, int *alloc, int used,
				 int complement) {
    int count = used;
    bin_index_t *bin = get_bin(io, bin_num);
    int i, f_a, f_b;

    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }

    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    if (!count) {
	/* Add root bin */
	if (count >= *alloc) {
	    *alloc = *alloc ? *alloc * 2 : 16;
	    *results = (rangec_t *)realloc(*results, *alloc *sizeof(rangec_t));
	}


	(*results)[count].rec   = bin_num;
	(*results)[count].start = NMIN(bin->pos, bin->pos+bin->size-1);
	(*results)[count].end   = NMAX(bin->pos, bin->pos+bin->size-1);
	(*results)[count].comp  = complement;
	count++;
    }

    /* Items in this case are child bins rather than sequences */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;

	ch = get_bin(io, bin->child[i]);
	if (count >= *alloc) {
	    *alloc = *alloc ? *alloc * 2 : 16;
	    *results = (rangec_t *)realloc(*results, *alloc *sizeof(rangec_t));
	}

	if (end   >= NMIN(ch->pos, ch->pos + ch->size-1) &&
	    start <= NMAX(ch->pos, ch->pos + ch->size-1)) {
	    (*results)[count].rec   = bin->child[i];
	    (*results)[count].start = NMIN(ch->pos, ch->pos+ch->size-1);
	    (*results)[count].end   = NMAX(ch->pos, ch->pos+ch->size-1);
	    (*results)[count].comp  = complement;
	    count++;

	    /* Recurse too */
	    count = contig_bins_in_range2(io, bin->child[i], start, end,
					  NMIN(ch->pos, ch->pos + ch->size-1),
					  results, alloc, count, complement);
	}
    }

    return count;
}

rangec_t *contig_bins_in_range(GapIO *io, contig_t **c, int start, int end,
			       int *count) {
    rangec_t *r = NULL;
    int alloc = 0;

    //    *count = contig_bins_in_range2(io, contig_get_bin(c), start, end,
    //				   contig_offset(io, c), &r, &alloc, 0, 0);
    *count = contig_bins_in_range2(io, contig_get_bin(c), start, end,
				   0, &r, &alloc, 0, 0);
    
    return r;
}


static void contig_bin_dump2(GapIO *io, int bin_num, int level) {
    int i;
    bin_index_t *bin = get_bin(io, bin_num);

    cache_incr(io, bin);

    printf("%*sBin %d\tpos %d\tsize %d\tleft %d\tright %d\tflags %d\n", 
	   level*4, "", bin->rec, bin->pos, bin->size,
	   bin->child[0], bin->child[1], bin->flags);

    for (i = 0; i < 2; i++) {
	if (!bin->child[i])
	    continue;
	contig_bin_dump2(io, bin->child[i], level+1);
    }

    cache_decr(io, bin);
}

void contig_bin_dump(GapIO *io, int cnum) {
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, cnum);
    contig_bin_dump2(io, contig_get_bin(&c), 0);
}

/* ---------------------------------------------------------------------- */
/* iterators on a contig to allow scanning back and forth */

/*
 * Populates a sequence iterator.
 * Returns 0 on success
 *        -1 on failure (typically fallen off the end of the contig)
 */
static int range_populate(GapIO *io, contig_iterator *ci,
			  int cnum, int start, int end) {
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, cnum);
    
    if (ci->r) {
	free(ci->r);
	ci->r = NULL;
    }

    if (!ci->auto_extend) {
	/* cstart..cend are hard limits */
	if (end < ci->cstart)
	    return -1;
	if (start > ci->cend)
	    return -1;

	if (start < ci->cstart)
	    start = ci->cstart;
	if (end > ci->cend)
	    end = ci->cend;
    } else {
	/* otherwise contig boundaries are hard limits */
	if (end < c->start)
	    return -1;
	if (start > c->end)
	    return -1;

	if (start < c->start)
	    start = c->start;
	if (end > c->end)
	    end = c->end;
    }

    ci->cnum = cnum;
    ci->start = start;
    ci->end = end;
    ci->r = contig_seqs_in_range(io, &c, start, end,
				 CSIR_SORT_BY_X, &ci->nitems);
    ci->index = 0;
    return 0;
}

/*
 * Dealloctes memory used by a contig_iterator including the cached ranges.
 */
void contig_iter_del(contig_iterator *ci) {
    if (!ci)
	return;

    if (ci->r)
	free(ci->r);
    free(ci);
}

/*
 * Allocates and initialises a new contig_iterator struct.
 * 'whence' may be either CITER_FIRST or CITER_LAST and it controls whether
 * we start the iteration point from the beginning or end of the list.
 * The start and end parameters dictate the initial region to query. We
 * may specify them as either coordinates or use CITER_CSTART and CITER_CEND
 * as synonyms for the first and last coordinate in the contig.
 * Finally auto_extend controls whether the start..end range is just a
 * location to start iterating from (auto_extend == 1) or a hard limit
 * with no desire to iterate past that range (auto_extend == 0).
 *
 * Returns contig_iterator pointer on success
 *         NULL on failure
 */
contig_iterator *contig_iter_new(GapIO *io, int cnum, int auto_extend,
				 int whence, int start, int end) {
    contig_iterator *ci = (contig_iterator *)malloc(sizeof(*ci));
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, cnum);
    
    if (!ci)
	return NULL;

    if (!c)
	return NULL;

    ci->r = NULL;
    ci->nitems = 0;
    ci->index = 0;
    ci->auto_extend = auto_extend;

    ci->cstart = start == CITER_CSTART ? c->start : start;
    ci->cend   =   end == CITER_CEND   ? c->end   : end;

    if (whence == CITER_FIRST) {
	start = ci->cstart;
	end   = start + 1000;
    } else {
	end   = ci->cend;
	start = end - 1000;
    }

    if (0 != range_populate(io, ci, cnum, start, end)) {
	contig_iter_del(ci);
	return NULL;
    }

    if (whence == CITER_LAST) {
	ci->index = ci->nitems-1;
    }

    return ci;
}

/*
 * Finds the index of snum in a range array.
 * Returns index on success
 *         -1 on failure
 */
static int find_range_index(contig_iterator *ci, int snum) {
    int i;
    for (i = 0; i < ci->nitems; i++) {
	if (ci->r[i].rec == snum)
	    break;
    }

    return (i == ci->nitems) ? -1 : i;
}

/*
 * Returns seq range_t struct pointer on success
 *        NULL on failure (typically fallen off the end of the contig)
 */
rangec_t *contig_iter_next(GapIO *io, contig_iterator *ci) {
    int curr;

    while (ci->index >= ci->nitems) {
	/* Fallen off the range edge */
	//	if (!ci->auto_extend)
	//	    return NULL;

	curr = ci->r ? ci->r[ci->index-1].rec : -1;
	if (-1 == range_populate(io, ci, ci->cnum,
				 ci->start + 1000, ci->end + 1000))
	    return NULL;

	if (ci->r && curr != -1)
	    ci->index = find_range_index(ci, curr)+1;
	else
	    ci->index = 0;
    }

    return &ci->r[ci->index++];
}

/*
 * Returns seq range_t struct pointer on success
 *        NULL on failure (typically fallen off the end of the contig)
 */
rangec_t *contig_iter_prev(GapIO *io, contig_iterator *ci) {
    int curr;

    while (ci->index < 0 || ci->nitems == 0) {
	/* Fallen off the range edge */
	//	if (!ci->auto_extend)
	//	    return NULL;

	curr = ci->r ? ci->r[0].rec : -1;
	if (-1 == range_populate(io, ci, ci->cnum,
				 ci->start - 1000, ci->end - 1000))
	    return NULL;

	if (ci->r && curr != -1) {
	    ci->index = find_range_index(ci, curr);
	    if (ci->index == -1)
		ci->index = ci->nitems-1;
	    else
		ci->index--;
	} else {
	    ci->index = ci->nitems-1;
	}
    }

    return &ci->r[ci->index--];
}

/* Track values prior to resolution resampling */
typedef struct {
    double pos;
    int val;
} tvalues_t;

/*
 * Walks through the bin structure filling out the tv struct with
 * track information at at least bpv resolution (<= bpv).
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int contig_get_track2(GapIO *io, int bin_num, int start, int end,
			     int type, double bpv, int offset,
			     tvalues_t **tv, int *alloc, int used,
			     int complement) {
    int count = used;
    int i, f_a, f_b;
    bin_index_t *bin = get_bin(io, bin_num);
    int recurse = 1;
    static int depth=0;

    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }

    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    printf("%*scontig_get_track2 want %5.1f got ~%5.1f bin %d (%d+%d) l=%d r=%d\n",
	   depth, "",
	   bpv, bin->size / (double)RD_ELEMENTS, bin_num,  offset, bin->size,
	   bin->child[0], bin->child[1]);

    /*
     * Assuming the bin overlaps the region (we check it) then
     * we know the bpv this bin should have, so don't attempt to
     * regenerate it (which ends up in a circular loop) if this
     * isn't capable of holding the bpv required.
     */
    if (!(end < NMIN(0, bin->size-1) || start > NMAX(0, bin->size-1)) &&
	!(bin->size / RD_ELEMENTS > bpv && bin->size > RD_ELEMENTS)) {
	int i;
	track_t *t;

	printf("*query\n");

	if (!(t = bin_query_track(io, bin, type)))
	    return -1;

	printf("*track => %d items => bpv %f\n", t->nitems,
	       (double)bin->size / t->nitems);

	/* Only include data at higher resolution */
	if ((double)bin->size / t->nitems <= bpv) {
	    recurse = 0;

	    for (i = 0; i < t->nitems; i++) {
		while (count + t->nitems > *alloc) {
		    *alloc = *alloc ? *alloc * 2 : 16;
		    *tv = (tvalues_t *)realloc(*tv, *alloc * sizeof(**tv));
		}
		(*tv)[count].pos = NORM((double)i * bin->size / t->nitems);
		(*tv)[count].val = arr(int, t->data, i);
		count++;
	    }
	}
    }

    if (!recurse)
	return count;

    /* Recurse down bins */
    for (i = 0; i < 2 > 0; i++) {
	bin_index_t *ch;
	if (!bin->child[i]) {
	    /* No data available, so fill with zero values */
	    int len, nitems, j, offset;

	    if (i == 0 && bin->child[1]) {
		ch = get_bin(io, bin->child[1]);
		len = bin->size - ch->size;
		offset = 0;
	    } else if (i == 1 && bin->child[0]) {
		ch = get_bin(io, bin->child[0]);
		len = bin->size - ch->size;
		offset = ch->size;
	    } else {
		len = bin->size;
		offset = 0;
	    }
	    nitems = len / bpv;

	    for (j = 0; j < nitems; j++) {
		while (count+1 > *alloc) {
		    *alloc = *alloc ? *alloc * 2 : 16;
		    *tv = (tvalues_t *)realloc(*tv, *alloc * sizeof(**tv));
		}
		(*tv)[count].pos = NORM((double)j * len / nitems + offset);
		(*tv)[count].val = 0;
		count++;
	    }

	    continue;
	}

	ch = get_bin(io, bin->child[i]);
	if (end   >= NMIN(ch->pos, ch->pos + ch->size-1) &&
	    start <= NMAX(ch->pos, ch->pos + ch->size-1)) {
	    depth += 2;
	    count = contig_get_track2(io, bin->child[i], start, end,
				      type, bpv, 
				      NMIN(ch->pos, ch->pos + ch->size-1),
				      tv, alloc, count, complement);
	    depth -= 2;
	}
    }

    return count;
}

#define WL 3
/*
 * Queries and/or creates a track of a given type from a region of the
 * contig at a require resolution (bpv = bases per value).
 *
 * Returns track_t to free with free_track on success
 *         NULL on failure
 */
track_t *contig_get_track(GapIO *io, contig_t **c, int start, int end,
			  int type, double bpv) {
    int nele, i, j;
    track_t *t;
    tvalues_t *tv = NULL;
    int alloc = 0, count = 0;
    double last_pos;
    int last_val, bin_off;
    int *data, *data3;
    bin_index_t *start_bin;
    int start_bin_rec;

    printf("Query %d..%d bpv %f\n", start, end, bpv);

    nele = ceil((end-start+1)/bpv);
    bpv = (end-start+1)/nele;
    t = track_create_fake(type, nele);
    data = ArrayBase(int, t->data);

    start_bin = bin_for_range(io, c, start, end, 0, &bin_off);
    if (start_bin) {
	start_bin_rec = start_bin->rec;
    } else {
	start_bin_rec = contig_get_bin(c);
	bin_off = contig_offset(io, c);
    }

    /* Get out position/value array */
    count = contig_get_track2(io, start_bin_rec, start-bpv, end-bpv, type,
			      bpv/WL < 1 ? 1 : bpv/WL,
			      bin_off, &tv, &alloc, 0, 0);
    //count = square_wave(start-bpv, end+bpv, bpv/3, &tv);
    printf("generated %d pos/val pairs\n", count);

    if (count == 0) {
	/*
	 * We maybe querying beyond the ends of the contig, so just
	 * manufacture some defaults here if so.
	 */
	for (j = 0; j < nele; j++) {	
	    /* FIXME: or whatever value is appropriate for this track */
	    /* Depth => 0, GC% => 0.5? */
	    data[j] = 0;
	}
	free(tv);
	return t;
    }

    /* And now we need to resample it to fit our required resolution */
    last_pos = tv[0].pos;
    last_val = tv[0].val;
    /*
     * All elements in tvalues are higher resolution than bpv, but in
     * theory no more than double. Hence a simple linear downsample is
     * easiest and sufficient.
     */
    data3 = (int *)malloc(nele*sizeof(*data3)*WL);
    for (j = 0; j < count && tv[j].pos <= start; j++)
	;
    if (j > 0)
	j--;

    for (i = 0; i < nele*WL; i++) {
	double p = i * (end-start+1.0) / (nele*WL) + start;

	//	    while (p < tv[j].pos)
	//		j++;
	while (j < count && p > tv[j].pos)
	    j++;

	/*
	 * Boundary conditions - we aim not to hit them though by querying
	 * a large region than we average over.
	 */
	if (j >= count) {
	    data3[i] = count ? tv[count-1].val : 0;
	    continue;
	}

	if (j <= 0) {
	    data3[i] = p < 0 ? 0 : tv[0].val;
	    continue;
	}
	    
	assert(j >= 1 && j < count);
	assert(p >= tv[j-1].pos && p <= tv[j].pos);

#if 1
	/* Linear interpolation */
	data3[i] = (tv[j].val - tv[j-1].val) * (p - tv[j-1].pos) /
	    (double)(tv[j].pos - tv[j-1].pos) + tv[j-1].val;
#endif

#if 0
	/* Nearest neighbour */
	if (ABS(p - tv[j].pos) < ABS(p - tv[j-1].pos)) {
	    data3[i] = tv[j].val;
	} else {
	    data3[i] = tv[j-1].val;
	}
#endif
    }

    /* Filter down to produce data[] */
    for (i = 0; i < nele; i++) {
	/* For WL == 3 */
	//data[i] = (data3[i*3] + data3[i*3+1] + data3[i*3+2]) / 3;
	if (i*3-2 >= 0) {
	    data[i] = (data3[i*3-2] + data3[i*3-1] +
		       data3[i*3] +
		       data3[i*3+1] + data3[i*3+2]) / 5;
	} else {
	    data[i] = (data3[i*3] + data3[i*3+1] + data3[i*3+2]) / 3;
	}
    }

    free(data3);
    free(tv);

    return t;
}

/*
 * Destroys contig 'rec'.
 * Returns 0 for success
 *        -1 for failure
 */
int contig_destroy(GapIO *io, int rec) {
    int i, j;

    printf("Destroy contig rec %d\n", rec);
    io->contig_order = cache_rw(io, io->contig_order);
    io->db = cache_rw(io, io->db);

    for (i = j = 0; i < io->db->Ncontigs; i++) {
	if (arr(GCardinal, io->contig_order, i) == rec) {
	    continue;
	}
	arr(GCardinal, io->contig_order, j) =
	    arr(GCardinal, io->contig_order, i);
	j++;
    }

    if (i == j) {
	fprintf(stderr, "Attempted to remove unknown contig, rec %d\n",
		rec);
	return -1;
    }

    io->db->Ncontigs--;
    return 0;
}

/*
 * bin_num	The current bin being moved or split.
 * pos		The contig break point.
 * offset	The absolute positional offset of this bin in original contig
 * pleft	The parent bin/-contig in the left new contig
 * pright	The parent bin/-contig in the right new contig
 * child_no     0 or 1 - whether this bin is the left/right child of its parent
 */
#define H 20
#define G 20
static int bin_dump_recurse(GapIO *io, contig_t **c,
			    char *fn, int child,
			    int bin_num, int offset,
			    int level, int complement,
			    double rootx, double rooty) {
    int i, f_a, f_b;
    bin_index_t *bin = get_bin(io, bin_num);
    static FILE *gv = NULL;
    static double scale = 0;
    static double level_end[100];
    static double level_st[100];
    double g;

    cache_incr(io, bin);

    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }

    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    if (!gv) {
	int o = 20;
	gv = fopen(fn, "w+");
	if (!gv)
	    return -1;
	scale = 800.0 / bin->size;
	o -= offset * scale;
	fprintf(gv, "%%!\n/Times-Roman findfont\n6 scalefont\nsetfont\n"
		"90 %d translate\n90 rotate\nnewpath\n", o);

	for (i = 0; i < 100; i++) {
	    level_end[i] = -1e10;
	    level_st[i] = 1e10;
	}

	fprintf(gv, "%f %d moveto (%d) show\n",
		scale * (*c)->start, level*(H+G)+H+20+7, (*c)->start);
	fprintf(gv, "%f %d moveto (%d) show\n",
		scale * (*c)->end, level*(H+G)+H+20+7, (*c)->end);
	fprintf(gv, "2 setlinewidth\n"
		"%f %d moveto 0 3 rmoveto 0 -6 rlineto 0 3 rmoveto\n"
		"%f %d lineto 0 3 rmoveto 0 -6 rlineto 0 3 rmoveto stroke\n"
		"1 setlinewidth\n",
		scale * (*c)->start, level*(H+G)+H+20,
		scale * (*c)->end,    level*(H+G)+H+20);
    }

    for (i = -level; i < 100; i++) {
	double p1 = scale * offset;
	double p2 = scale * (offset + bin->size);
	if (p1 >= level_end[i] || p2 <= level_st[i])
	    break;
    }
    level = -i;

    if (level < -12)
	goto skip;

    /* linking lines */
    if (rootx || rooty) {
	g = child == 0 ? 0.6 : 0.9;
	fprintf(gv, "%f setgray %f %f moveto %f %d lineto stroke\n",
		g, rootx, rooty,
		scale*(offset+bin->size/2),
		level*(H+G)+H);
    }
    rootx = scale*(offset+bin->size/2);
    rooty = level*(H+G);

    if (level >= -8) {
	/* Inner grey box showing the used portion */
	fprintf(gv, "0.8 setgray %f %d moveto %f 0 rlineto 0 %d rlineto -%f 0 rlineto 0 -%d rlineto stroke\n",
		scale * NMIN(bin->start_used, bin->end_used), level*(H+G)+3,
		scale * (NMAX(bin->start_used, bin->end_used) -
			 NMIN(bin->start_used, bin->end_used)),
		H-6,
		scale * (NMAX(bin->start_used, bin->end_used) -
			 NMIN(bin->start_used, bin->end_used)),
		H-6);
	
	/* bin offset */
	fprintf(gv, "%f %d moveto (%d) show\n",
		scale * offset, level*(H+G)-7, offset);
    }

    /* Outer bin bounding box */
    g = complement ? 0.5 : 0;
    fprintf(gv, "%f setgray %f %d moveto %f 0 rlineto 0 %d rlineto -%f 0 rlineto 0 -%d rlineto stroke\n",
	    g, scale * offset, level*(H+G), scale * bin->size,
	    H, scale * bin->size, H);

    if (scale * offset < level_st[-level])
	level_st[-level] = scale * offset;
    if (scale * (offset + bin->size) > level_end[-level])
	level_end[-level] = scale * (offset + bin->size);

    /* bin start_used position */
    if (bin->rng)
	fprintf(gv, "%f %d moveto (%d) show\n",
		scale * NMIN(bin->start_used, bin->end_used),
		level*(H+G)+H+2,
		NMIN(bin->start_used, bin->end_used));
    else
	fprintf(gv, "%f %d moveto (empty) show\n",
		scale * offset,
		level*(H+G)+H+2);

    /* Bin record ID */
    fprintf(gv, "%f %d moveto (#%d) show\n",
	    scale * (offset+bin->size/2)-10, level*(H+G)+H/2-2, bin->rec);

    fflush(gv);

 skip:
    /* Recurse */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;

	ch = get_bin(io, bin->child[i]);
	if (0 != bin_dump_recurse(io, c, fn, i, bin->child[i],
				  NMIN(ch->pos, ch->pos + ch->size-1),
				  level-1, complement, rootx, rooty))
	    return -1;
    }

    cache_decr(io, bin);

    if (level == 0) {
	fprintf(gv, "showpage\n");
	fclose(gv);
	gv = NULL;
    }

    return 0;
}

/*
 * Produces a postscript file containing a plot of the contig bin structure.
 */
void contig_dump_ps(GapIO *io, contig_t **c, char *fn) {
    //HacheTableRefInfo(io->cache, stdout);
    bin_dump_recurse(io, c, fn, 0, contig_get_bin(c),
		     contig_offset(io, c), 0, 0, 0.0, 0.0);
    //HacheTableRefInfo(io->cache, stdout);
}
