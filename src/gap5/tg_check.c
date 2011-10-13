/*
 * Checks the validity of data in the database.
 *
 * These can be split into actual corruptions and credibility warnings.
 * Eg a sequence length not matching the start/end range_t is a fatal error,
 * but a sequence being clipped to zero bases is valid, just rather
 * unbelievable.
 */
#include <staden_config.h>
#include <string.h>
#include <stdlib.h>

#include "tg_gio.h"
#include "tg_check.h"
#include "hache_table.h"

/* Enable debugging, which may help track down the location of some errors */
//#define DEBUG_CHECK

#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

typedef struct {
    /* Contig start/end */
    int cstart;
    int cend;
    
    /* Item counts */
    int nseq;
    int nanno;
    int nref;
} bin_stats;

/*
 * Checks a sequence struct is consistent, both internally and with the
 * range that points to it.
 */
static int check_seq(GapIO *io, int fix, bin_index_t *bin, range_t *r,
		     HacheTable *lib_hash, int is_cons, int *fixed) {
    int err = 0;
    int i, len;
    seq_t *s = cache_search(io, GT_Seq, r->rec);

    static char isbase[256], init_done = 0;

    if (!init_done) {
	memset(isbase, 0, 256);
	isbase['A'] = isbase['a'] = 1;
	isbase['C'] = isbase['c'] = 1;
	isbase['G'] = isbase['g'] = 1;
	isbase['T'] = isbase['t'] = 1;
	isbase['N'] = isbase['n'] = 1;
	isbase['*'] = isbase['-'] = 1;
	isbase['0'] = isbase['1'] = 1;
	isbase['2'] = isbase['3'] = 1;
	init_done = 1;
    }

    if (!s) {
	vmessage("Seq %"PRIrec": failed to load\n", r->rec);
	if (fix && !io->base) {
	    /* Returns bin if io->base is NULL */
	    cache_rw(io, bin);
	    r->rec      = 0;
	    r->flags    = GRANGE_FLAG_UNUSED;
	    bin->flags |= BIN_RANGE_UPDATED;
	    if (fixed) (*fixed)++;
	}
	return 1;
    }

    if (!bin->rng ||
	s->bin_index != r - ArrayBase(range_t, bin->rng)) {
	vmessage("Seq %"PRIrec": bin_index does not match range index\n",
		 s->rec);
	err++;
	if (fix) {
	    s = cache_rw(io, s);
	    s->bin_index = r - ArrayBase(range_t, bin->rng);
	    if (fixed) if (fixed) (*fixed)++;
	}
    }

    if (s->bin != bin->rec) {
	vmessage("Seq %"PRIrec": bin does not match observed bin\n", s->rec);
	err++;
	if (fix) {
	    s = cache_rw(io, s);
	    s->bin = bin->rec;
	    if (fixed) (*fixed)++;
	}
    }
    if (r->end - r->start + 1 != ABS(s->len)) {
	vmessage("Seq %"PRIrec": length does not match bin-range\n", s->rec);
	err++;
	if (fix && !io->base) {
	    cache_rw(io, bin);
	    /* Best guess - it shows correctly now so amend appropriate end */
	    if (bin->flags & BIN_COMPLEMENTED) {
		r->start = r->end - (ABS(s->len)-1);
	    } else {
		r->end = ABS(s->len)-1 + r->start;
	    }
	    bin->flags |= BIN_RANGE_UPDATED;
	    if (fixed) (*fixed)++;
	}
    }

    if (s->left < 1 || s->right > ABS(s->len)) {
	vmessage("Seq %"PRIrec": left/right clips outside of sequence "
		 "bounds.\n", s->rec);
	err++;
	if (fix) {
	    s = cache_rw(io, s);
	    if (s->left < 1)
		s->left = 1;
	    if (s->right > ABS(s->len))
		s->right = s->len;
	    if (fixed) (*fixed)++;
	}
    }

    if (s->right < s->left-1) {
	/* Allow 1 diff as this is a seq with no visible data */
	vmessage("Seq %"PRIrec": right clip starts before left clip\n",
		 s->rec);
	err++;
	if (fix) {
	    s = cache_rw(io, s);
	    s->left = s->right;
	    if (fixed) (*fixed)++;
	}
    }

    if (s->mapping_qual != r->mqual) {
	vmessage("Seq %"PRIrec": mapping_qual disagrees with range\n",
		 s->rec);
	err++;
	if (fix) {
	    s = cache_rw(io, s);
	    s->mapping_qual = r->mqual;
	    if (fixed) (*fixed)++;
	}
    }

    /* TODO: Check r->pair_rec? */

    /* TODO: Check r->flags match s->flags */

    if (s->seq_tech > STECH_454) {
	vmessage("Seq %"PRIrec": Unknown seq_tech value\n", s->rec);
	err++;
    }

    if (s->format > 3) {
	vmessage("Seq %"PRIrec": Unknown format value\n", s->rec);
	err++;
    }

    if (s->format != SEQ_FORMAT_CNF1) {
	vmessage("Seq %"PRIrec": Warning - outdated seq format\n", s->rec);
    }

    /* Parent rec / type */
    if (s->parent_rec && !cache_exists(io, s->parent_type, s->parent_rec)) {
	vmessage("Seq %"PRIrec": parent_rec/type do not match\n", s->rec);
	err++;
    }

    if (lib_hash && s->parent_type == GT_Library) {
	HacheItem *hi = HacheTableQuery(lib_hash, (char *)&s->parent_rec,
					sizeof(s->parent_rec));
	if (!hi) {
	    HacheData hd;
	    vmessage("Seq %"PRIrec": parent_rec %"PRIrec" of type "
		     "GT_Library is not in db->library array\n",
		     s->rec, s->parent_rec);
	    hd.i = 0;
	    HacheTableAdd(lib_hash, (char *)&s->parent_rec,
			  sizeof(s->parent_rec), hd, NULL);
	    err++;
	}
    }
    

    /* Check valid sequence characters */
    len = ABS(s->len);
    for (i = 0; i < len; i++) {
	if (!isbase[(unsigned char)s->seq[i]] &&
	    !(s->seq[i] == ' ' && is_cons)) {
	    vmessage("Seq %"PRIrec": seq has unexpected char '%c' (%d)\n",
		     s->rec, isprint(s->seq[i]) ? s->seq[i] : '?', s->seq[i]);
	    err++;
	    break;
	}
    }

    if ((s->name && s->name_len != strlen(s->name)) ||
	(!s->name && s->name_len != 0)) {
	vmessage("Seq %"PRIrec": name and name_len do not match\n", s->rec);
	err++;
    }

    return err;
}

/*
 * Checks a anno_ele struct is consistent, both internally and with the
 * range that points to it.
 */
static int check_anno(GapIO *io, bin_index_t *bin, range_t *r,
		      HacheTable *rec_hash, int db_vers,
		      int valid_ctg_start, int valid_ctg_end) {
    int err = 0;
    anno_ele_t *a = cache_search(io, GT_AnnoEle, r->rec);

    cache_incr(io, a);

    /* Bin records match */
    if (a->bin != bin->rec) {
	vmessage("Anno %"PRIrec": bin does not match observed bin\n", a->rec);
	err++;
    }

    /* Types match */
    switch (a->obj_type) {
    case GT_Contig:
	if (r->flags & GRANGE_FLAG_TAG_SEQ) {
	    vmessage("Anno %"PRIrec": range flags indicate sequence, but anno "
		     "obj_type claims contig\n", a->rec);
	    err++;
	}
	break;

    case GT_Seq: {
	/* Check if annotation is in the same bin as seq */
	if (rec_hash && db_vers >= 3) {
	    HacheItem *hi = HacheTableQuery(rec_hash, (char *)&a->obj_rec,
					    sizeof(a->obj_rec));
	    if (!hi || hi->data.i != a->bin) {
		vmessage("Anno %"PRIrec": attached to seq %"PRIrec" held "
			 "within a different bin.\n", a->rec, a->obj_rec);
		err++;
	    }
	}

	if (!(r->flags & GRANGE_FLAG_TAG_SEQ)) {
	    vmessage("Anno %"PRIrec": range flags indicate contig, but anno "
		     "obj_type claims sequence\n", a->rec);
	    err++;
	}
	break;
    }

    default:
	vmessage("Anno %"PRIrec": unrecognised value (%d) for obj_type\n",
		 a->rec, a->obj_type);
	err++;
    }

    /* That anno links to the correct data type */
    if (a->obj_rec && !cache_exists(io, a->obj_type, a->obj_rec)) {
	vmessage("Anno %"PRIrec": obj_rec/type do not match\n", a->rec);
	err++;
    }

    /* Consistency of cached portions in range */
    if (a->tag_type != r->mqual) {
	vmessage("Anno %"PRIrec": type does not match copy in "
		"range_t struct.\n", a->rec);
	err++;
    }
    if (a->obj_type != GT_Contig && a->obj_rec != r->pair_rec) {
	vmessage("Anno %"PRIrec": obj_rec does not match copy "
		 "in range_t struct.\n", a->rec);
	err++;
    }

    /*
     * Also check position of item being annotated - it should overlap
     * the tag and be in the same contig. This is a time consuming check
     * though.
     */
    {
	int astart, aend, ostart, oend;
	tg_rec acontig, ocontig;

	if (-1 == bin_get_item_position(io, GT_AnnoEle, a->rec,
					&acontig, &astart, &aend, NULL,
					NULL, NULL, NULL)) {
	    cache_decr(io, a);
	    return err+1;
	}

	if (a->obj_type == GT_Contig) {
	    ocontig = a->obj_rec;
	    if (ocontig && ocontig != acontig) {
		/* Previously an error, but removed assumptions on this */
		vmessage("Anno %"PRIrec": non-zero obj_rec with obj_type "
			 "GT_Contig is deprecated\n", a->rec);
	    }
	    ocontig = acontig;
	    ostart  = valid_ctg_start;
	    oend    = valid_ctg_end;
	} else {
	    if (-1 == bin_get_item_position(io, a->obj_type, a->obj_rec,
					    &ocontig, &ostart, &oend, NULL,
					    NULL, NULL, NULL)) {
		cache_decr(io, a);
		return err+1;
	    }
	}

	if (acontig != ocontig || astart < ostart || aend > oend) {
	    vmessage("Anno %"PRIrec": does not overlap annotated object %d/%"
		     PRIrec":\n", a->rec, a->obj_type, a->obj_rec);
	    vmessage("\tTag Ctg %"PRIrec" at %d..%d\n", acontig, astart, aend);
	    vmessage("\tObj Ctg %"PRIrec" at %d..%d\n", ocontig, ostart, oend);
	    
	    err++;
	}
    }

    cache_decr(io, a);
    return err;
}

/*
 * Checks a range of type GRANGE_FLAG_ISREFPOS
 */
static int check_refpos(GapIO *io, range_t *r) {
    int err = 0;

    if (r->start != r->end) {
	vmessage("RefPos %"PRIrec": start and end positions differ\n", r->rec);
	err++;
    }

    return err;
}

/* Adds 'delta' to the start/end coordinate of all ranges. */
static void bin_shift_range(GapIO *io, bin_index_t *bin, int delta) {
    int i;
    int start = INT_MAX;
    int end   = INT_MIN;
    int valid_range = 0;

    printf("Shift range for bin %"PRIrec"\n", bin->rec);

    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);

	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	r->start += delta;
	r->end   += delta;

	if (start > r->start)
	    start = r->start;
	if (end   < r->end)
	    end   = r->end;

	valid_range = 1;
    }

    if (valid_range) {
	bin->start_used = start;
	bin->end_used = end;
    }

    bin->flags |= BIN_RANGE_UPDATED;
}

/* Shifts child bins by delta */
static void bin_shift_children(GapIO *io, bin_index_t *bin, int delta) {
    int i;
    bin_index_t *ch;

    for (i = 0; i < 2; i++) {
	if (bin->child[i] == 0)
	    continue;

	ch = cache_search(io, GT_Bin, bin->child[i]);
	ch = cache_rw(io, ch);
	ch->flags |= BIN_BIN_UPDATED;
	ch->pos += delta;
    }
}

/*
 * Ensures that the parent bin is large enough to cover this bin. Grow it
 * if necessary.
 */
void grow_bin(GapIO *io, bin_index_t *bin) {
    bin_index_t *parent = bin;

    cache_incr(io, bin);

    while (bin->parent_type == GT_Bin) {
	/* Messy */
	int par_par_comp = 0;
	if (parent->parent_type == GT_Bin) {
	    bin_index_t *ppbin = cache_search(io, GT_Bin, parent->parent);
	    if (ppbin->flags & BIN_COMPLEMENTED)
		par_par_comp = 1;
	}


	parent = cache_search(io, GT_Bin, bin->parent);
	cache_incr(io, parent);

	if (parent->flags & BIN_COMPLEMENTED) {
	    int delta = -bin->pos;

	    if (bin->pos < 0) {
		if (par_par_comp == 0) {
		    parent = cache_rw(io, parent);
		    parent->size += delta;
		    parent->flags |= BIN_BIN_UPDATED;
		} else {
		    /* Get r/w copies so pointers don't change in sub-funcs */
		    bin = cache_rw(io, bin);
		    parent = cache_rw(io, parent);

		    parent->pos -= delta;
		    parent->flags |= BIN_BIN_UPDATED;
		    bin_shift_range(io, parent, delta);
		    bin_shift_children(io, parent, delta);
		}
	    }

	    if (bin->pos + bin->size > parent->size) {
		int delta = bin->pos + bin->size - parent->size;

		if (par_par_comp == 0) {
		    bin = cache_rw(io, bin);
		    parent = cache_rw(io, parent);

		    parent->pos -= delta;
		    parent->flags |= BIN_BIN_UPDATED;
		    bin_shift_range(io, parent, delta);
		    bin_shift_children(io, parent, delta);
		} else {
		    parent = cache_rw(io, parent);
		    parent->size += delta;
		    parent->flags |= BIN_BIN_UPDATED;
		}
	    }
	} else {
	    if (bin->pos < 0) {
		int delta = -bin->pos;

		if (par_par_comp == 0) {
		    /* Get r/w copies so pointers don't change in sub-funcs */
		    bin = cache_rw(io, bin);
		    parent = cache_rw(io, parent);

		    parent->pos -= delta;
		    parent->flags |= BIN_BIN_UPDATED;
		    bin_shift_range(io, parent, delta);
		    bin_shift_children(io, parent, delta);
		} else {
		    parent = cache_rw(io, parent);
		    parent->size += delta;
		    parent->flags |= BIN_BIN_UPDATED;
		}
	    }

	    if (bin->pos + bin->size > parent->size) {
		int delta = bin->pos + bin->size - parent->size;

		if (par_par_comp == 0) {
		    parent = cache_rw(io, parent);
		    parent->size += delta;
		    parent->flags |= BIN_BIN_UPDATED;
		} else {
		    bin = cache_rw(io, bin);
		    parent = cache_rw(io, parent);

		    parent->pos -= delta;
		    parent->flags |= BIN_BIN_UPDATED;
		    bin_shift_range(io, parent, delta);
		    bin_shift_children(io, parent, delta);
		}
	    }
	}

	cache_decr(io, bin);
	bin = parent;
    }

    cache_decr(io, parent);
}

/*
 * Walks a contig bin tree, executing callbacks per bin.
 *
 * Returns 0 on success
 *         number of errors on failure
 */
static int bin_walk(GapIO *io, int fix, tg_rec rec, int offset, int complement,
		    int level, HacheTable *lib_hash,
		    HacheTable *rec_hash, bin_stats *bs,
		    int valid_ctg_start, int valid_ctg_end, int *fixed) {
    bin_index_t *bin;
    int i, f_a, f_b, err = 0;
    bin_stats child_stats;
    int start, end, cstart, cend, nthis_seq = 0;
    int valid_range;
    int db_vers = io->base ? io->base->db->version : io->db->version;

    if (!rec)
	return 1;

    if (NULL == (bin = cache_search(io, GT_Bin, rec)))
	return 1;

    /* Set f_a & f_b for NMIN/NMAX macros & handle complementing */
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

#ifdef DEBUG_CHECK
    printf("Bin %"PRIrec"  %d..%d  used %d..%d nseq=%d\n",
	   bin->rec,
	   NMIN(0, bin->size-1),
	   NMAX(0, bin->size-1),
	   NMIN(bin->start_used, bin->end_used),
	   NMAX(bin->start_used, bin->end_used),
	   bin->nseqs);
#endif

    cache_incr(io, bin);

    /* Add recs to the rec_hash */
    if (rec_hash) {
	for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	    range_t *r = arrp(range_t, bin->rng, i);
	    HacheData hd;
	    int new;

	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;

	    if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS)
		continue;

	    hd.i = bin->rec;
	    HacheTableAdd(rec_hash, (char *)&r->rec, sizeof(r->rec), hd, &new);
	    if (!new) {
		vmessage("Rec %"PRIrec" occurs more than once in a bin tree\n",
			 r->rec);
		err++;
	    }
	}
    }

    /* Recurse */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;
	ch = get_bin(io, bin->child[i]);
	if (!ch) {
	    vmessage("bin %"PRIrec" failed to load.\n",
		     bin->child[i]);
	    err++;
	    if (fix) {
		bin = cache_rw(io, bin);
		bin->flags |= BIN_BIN_UPDATED;
		bin->child[i] = 0;
		if (fixed) (*fixed)++;
	    }
	    continue;
	}

	child_stats.cstart = INT_MAX;
	child_stats.cend   = INT_MIN;
	child_stats.nseq   = 0;
	child_stats.nanno  = 0;
	child_stats.nref   = 0;

	if (ch->parent != bin->rec || ch->parent_type != GT_Bin) {
	    vmessage("bin %"PRIrec" parent/type is incorrect\n", ch->rec);
	    err++;
	    if (fix) {
		ch = cache_rw(io, ch);
		ch->parent = bin->rec;
		ch->parent_type = GT_Bin;
		if (fixed) (*fixed)++;
	    }
	}

	if (ch->size <= 0) {
	    vmessage("bin %"PRIrec" has size <= 0 (%d)\n",
		     ch->rec, ch->size);
	    err++;
	    if (fix) {
		/* Just "lose" it */
		bin = cache_rw(io, bin);
		bin->flags |= BIN_BIN_UPDATED;
		bin->child[i] = 0;
		continue;
		if (fixed) (*fixed)++;
	    }

	} else  {
	    if (NMIN(ch->pos, ch->pos + ch->size-1) < NMIN(0, bin->size-1)) {
		vmessage("bin %"PRIrec" has child %d rec %"PRIrec" with left "
			 "edge less than the parent. %d < %d\n",
			 bin->rec, i, bin->child[i],
			 NMIN(ch->pos, ch->pos + ch->size-1),
			 NMIN(0, bin->size-1));
		err++;
		if (fix) {
		    grow_bin(io, bin);
		    if (fixed) (*fixed)++;
		}
	    }

	    if (NMAX(ch->pos, ch->pos + ch->size-1) > NMAX(0, bin->size-1)) {
		vmessage("bin %"PRIrec" has child %d rec %"PRIrec" with right "
			 "edge greater than the parent. %d > %d\n",
			 bin->rec, i, bin->child[i],
			 NMAX(ch->pos, ch->pos + ch->size-1),
			 NMAX(0, bin->size-1));
		err++;
	    }
	}

	err += bin_walk(io, fix, bin->child[i],
			NMIN(ch->pos, ch->pos + ch->size-1) /* offset */,
			complement, level, lib_hash,
			rec_hash, &child_stats,
			valid_ctg_start, valid_ctg_end, fixed);

	bs->nseq  += child_stats.nseq;
	bs->nanno += child_stats.nanno;
	bs->nref  += child_stats.nref;

	if (child_stats.nseq) {
	    if (bs->cstart > child_stats.cstart)
		bs->cstart = child_stats.cstart;
	    if (bs->cend   < child_stats.cend)
		bs->cend   = child_stats.cend;
	}
    }


    /* Check range items */
    start  = INT_MAX;
    end    = INT_MIN;
    cstart = INT_MAX;
    cend   = INT_MIN;
    valid_range = 0;
    if (bin->rng) {
	char *loop;
	int last_range = -1;

	/* Check rng_free list */
	if (bin->rng_free < -1 ||
	    (bin->rng_free != -1 && bin->rng_free >= ArrayMax(bin->rng))) {
	    vmessage("bin %"PRIrec": rng_free has invalid value %d\n",
		     bin->rec, bin->rng_free);
	    err++;
	    if (fix) {
		bin = cache_rw(io, bin);
		bin->rng_free = -1;
		bin->flags |= BIN_BIN_UPDATED;
		if (fixed) (*fixed)++;
	    } else {
		goto skip_checks;
	    }
	}

	loop = (char *)calloc(ArrayMax(bin->rng), sizeof(char));
	i = bin->rng_free;
	while (i != -1) { 
	    range_t *r = arrp(range_t, bin->rng, i);

	    if (!(r->flags & GRANGE_FLAG_UNUSED)) {
		vmessage("bin %"PRIrec" free range item %d: "
			 "flagged as in use\n", bin->rec, i);
		err++;
		if (fix) {
		    bin = cache_rw(io, bin);
		    /* Terminate rng_free list and just accept the
		     * chance of leaking.
		     */
		    if (last_range != -1) {
			r = arrp(range_t, bin->rng, last_range);
			r->rec = -1;
			bin->flags |= BIN_RANGE_UPDATED;
		    } else {
			bin->rng_free = -1;
			bin->flags |= BIN_BIN_UPDATED;
		    }
		    if (fixed) (*fixed)++;
		}
		break;
	    }

	    if ((int)r->rec < -1 || 
		((int)r->rec != -1 && (int)r->rec >= ArrayMax(bin->rng))) {
		vmessage("bin %"PRIrec" free range item %d: "
			 "next rec (r->rec) is invalid\n", bin->rec, i);
		err++;
		if (fix) {
		    bin = cache_rw(io, bin);
		    r->rec = -1;
		    bin->flags |= BIN_RANGE_UPDATED;
		    if (fixed) (*fixed)++;
		}
		break;
	    }

	    /* Loop detection */
	    if (loop[i]) {
		vmessage("bin %"PRIrec": loop detected in free list\n",
			 bin->rec);
		err++;
		if (fix) {
		    bin = cache_rw(io, bin);
		    if (last_range != -1) {
			r = arrp(range_t, bin->rng, last_range);
			r->rec = -1;
			bin->flags |= BIN_RANGE_UPDATED;
		    } else {
			bin->rng_free = -1;
			bin->flags |= BIN_BIN_UPDATED;
		    }
		    if (fixed) (*fixed)++;
		}
		break;
	    }
	    loop[i] = 1;

	    last_range = i;
	    i = (int)r->rec;
	}
	free(loop);

    skip_checks:

	/* Iterate through USED bin items and check them */
	for (i = 0; i < ArrayMax(bin->rng); i++) {
	    range_t *r = arrp(range_t, bin->rng, i);

	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;

	    valid_range = 1;

#ifdef DEBUG_CHECK
	    printf("#%"PRIrec": Range item %d (%"PRIrec" flag %d): %d..%d "
		   "(abs %d..%d)\n",
		   bin->rec, i, r->rec, r->flags, r->start, r->end,
		   NMIN(r->start, r->end), NMAX(r->start, r->end));
#endif
	    if (start > r->start)
		start = r->start;
	    if (end   < r->end)
		end   = r->end;

	    switch (r->flags & GRANGE_FLAG_ISMASK) {
	    case GRANGE_FLAG_ISSEQ: {
		bs->nseq++;
		nthis_seq++;
		if (level > 1)
		    err += check_seq(io, fix, bin, r, lib_hash, 0, fixed);

		if (cstart > r->start)
		    cstart = r->start;
		if (cend   < r->end)
		    cend   = r->end;
		break;
	    }

	    case GRANGE_FLAG_ISANNO:
		bs->nanno++;
		if (level > 1)
		    err += check_anno(io, bin, r, rec_hash, db_vers,
				      valid_ctg_start, valid_ctg_end);
		break;

	    case GRANGE_FLAG_ISREFPOS:
		bs->nref++;
		if (level > 1)
		    err += check_refpos(io, r);
		break;

	    case GRANGE_FLAG_ISCONS:
		if (level > 1)
		    err += check_seq(io, fix, bin, r, NULL, 1, fixed);
		break;

	    case GRANGE_FLAG_ISUMSEQ:
	    case GRANGE_FLAG_ISREF:
		/* FIXME: check */
		break;

	    default:
		/* Unknown ISMASK flag */
		vmessage("bin %"PRIrec" range item %d: "
			 "Unknown GRANGE_FLAG_IS? flag: %d\n",
			 bin->rec, i, r->flags & GRANGE_FLAG_ISMASK);
		err++;
	    }
	}
    }


    /* Check count validity to ensure this + children are correct */
    if (bin->nseqs != bs->nseq) {
	vmessage("bin %"PRIrec": nseqs does not match observed counts\n",
		 bin->rec);
	if (fix) {
	    bin = cache_rw(io, bin);
	    bin->flags |= BIN_BIN_UPDATED;
	    bin->nseqs = bs->nseq;
	    if (fixed) (*fixed)++;
	}
	err++;
    }
    if (db_vers > 1 && bin->nanno != bs->nanno) {
	vmessage("bin %"PRIrec": nanno does not match observed counts\n",
		 bin->rec);
	if (fix) {
	    bin = cache_rw(io, bin);
	    bin->flags |= BIN_BIN_UPDATED;
	    bin->nanno = bs->nanno;
	    if (fixed) (*fixed)++;
	}
	err++;
    } else if (fix && db_vers == 1 && bin->nanno != bs->nanno) {
	vmessage("bin %"PRIrec": fixing nanno\n", bin->rec);
	bin = cache_rw(io, bin);
	bin->flags |= BIN_BIN_UPDATED;
	bin->nanno = bs->nanno;
	if (fixed) (*fixed)++;
    }
    if (db_vers > 1 && bin->nrefpos != bs->nref) {
	vmessage("bin %"PRIrec": nrefpos does not match observed counts\n",
		 bin->rec);
	if (fix) {
	    bin = cache_rw(io, bin);
	    bin->flags |= BIN_BIN_UPDATED;
	    bin->nrefpos = bs->nref;
	    if (fixed) (*fixed)++;
	}
	err++;
    } else if (fix && db_vers == 1 && bin->nrefpos != bs->nref) {
	vmessage("bin %"PRIrec": fixing nrefpos\n", bin->rec);
	bin = cache_rw(io, bin);
	bin->flags |= BIN_BIN_UPDATED;
	bin->nrefpos = bs->nref;
	if (fixed) (*fixed)++;
    }


    /*
     * Check used start/end range, and accumulate absolute positions so we
     * can check the contig later.
     */
    if (valid_range) {
	if (start != bin->start_used ||
	    end   != bin->end_used) {
	    vmessage("bin %"PRIrec": used start/end range are incorrect\n",
		     bin->rec);
	    err++;
	    if (fix) {
		bin = cache_rw(io, bin);
		bin->flags |= BIN_BIN_UPDATED;
		bin->start_used = start;
		bin->end_used = end;
		if (fixed) (*fixed)++;
	    }
	}

	if (nthis_seq) {
	    if (bs->cstart > NMIN(cstart, cend))
		bs->cstart = NMIN(cstart, cend);
	    if (bs->cend   < NMAX(cstart, cend))
		bs->cend   = NMAX(cstart, cend);
	}

	if (bin->start_used < 0 || bin->end_used >= bin->size) {
	    vmessage("bin %"PRIrec": used start/end range beyond the bin "
		     "boundaries (size %d vs start=%d,end=%d).\n",
		     bin->rec, bin->size, bin->start_used, bin->end_used);
	    err++;
	    if (fix) {
		bin = cache_rw(io, bin);
		bin->flags |= BIN_BIN_UPDATED | BIN_RANGE_UPDATED;

		if (bin->start_used < 0) {
		    bin->pos  += bin->start_used;
		    bin_shift_range(io, bin, -bin->start_used);
		}

		if (bin->end_used >= bin->size) {
		    bin->size = bin->end_used+1;
		}

		/* Now check it hasn't outgrown the parent too */
		grow_bin(io, bin);
		if (fixed) (*fixed)++;
	    }
	}

    } else {
	if (bin->start_used != 0 || bin->end_used != 0) {
	    vmessage("bin %"PRIrec": used start/end are non-zero "
		     "in an empty bin\n", bin->rec);
	    err++;
	    if (fix) {
		bin = cache_rw(io, bin);
		bin->flags |= BIN_BIN_UPDATED;
		bin->start_used = 0;
		bin->end_used = 0;
		if (fixed) (*fixed)++;
	    }
	}
    }

    /* Remove items from the rec hash */
    if (rec_hash) {
	for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	    range_t *r = arrp(range_t, bin->rng, i);

	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;

	    if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS)
		continue;

	    if (HacheTableRemove(rec_hash, (char *)&r->rec,
				 sizeof(r->rec), 0)) {
		vmessage("Failed to remove rec %"PRIrec" from rec_hash\n",
			 r->rec);
		err++;
	    }
	}
    }

    cache_decr(io, bin);

    return err;
}


/*
 * Checks a single contig. lib_hash may be NULL, but if not then it is a
 * HacheTable keyed on library record numbers. We then validate all sequences
 * against this to check that no unknown libraries are present.
 *
 * Returns the number of errors found
 *         0 on success.
 */
int check_contig(GapIO *io, tg_rec crec, int fix, int level,
		 HacheTable *lib_hash, int *fixed) {
    contig_t *c;
    bin_stats bs;
    int err = 0;
    bin_index_t *bin;
    HacheTable *rec_hash = NULL;
    int valid_ctg_start, valid_ctg_end;

    consensus_valid_range(io, crec, &valid_ctg_start, &valid_ctg_end);

    if (level > 1)
	rec_hash = HacheTableCreate(65536,
				    HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);

    if (!cache_exists(io, GT_Contig, crec)) {
	vmessage("Record %"PRIrec" is not a contig, but in contig order\n",
		 crec);
	return 1;
    }

    c = cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

    bs.cstart = INT_MAX;
    bs.cend   = INT_MIN;
    bs.nseq   = 0;
    bs.nanno  = 0;
    bs.nref   = 0;

    if (c->bin) {
	bin = cache_search(io, GT_Bin, c->bin);
	if (bin->parent != crec || bin->parent_type != GT_Contig) {
	    vmessage("root bin %"PRIrec" parent/type is incorrect\n", c->bin);
	    err++;
	    if (fix) {
		bin = cache_rw(io, bin);
		bin->parent = crec;
		bin->parent_type = GT_Contig;
		if (fixed) (*fixed)++;
	    }
	}
    }

    err += bin_walk(io, fix, c->bin, contig_offset(io, &c), 0, level,
		    lib_hash, rec_hash, &bs,
		    valid_ctg_start, valid_ctg_end, fixed);

    if (bs.cstart != c->start ||
	bs.cend   != c->end) {
	vmessage("Contig %"PRIrec": used start/end range are incorrect\n",
		 crec);
	err++;
	if (fix) {
	    c = cache_rw(io, c);
	    c->start = bs.cstart;
	    c->end   = bs.cend;
	    if (fixed) (*fixed)++;
	}
    }
    
    cache_decr(io, c);
    if (rec_hash)
	HacheTableDestroy(rec_hash, 0);

    return err;
}

int check_cache(GapIO *io) {
    GapIO *ior = gio_open(io->name, 1, 0);
    HacheTable *h = io->cache;
    int i, j, err = 0;

    if (!ior)
	return 1;

    for (i = 0; i < h->nbuckets; i++) {
	HacheItem *hi, *next;
	void *v;

	for (hi = h->bucket[i]; hi; hi = next) {
	    cached_item *ci = hi->data.p, *ci2;
	    next = hi->next;
	    int mis = 0;

	    v = cache_search(ior, ci->type, ci->rec);
	    if (!v) {
		vmessage("Failed to find rec %"PRIrec" of type %d in disk "
			 "copy, but it is present in memory cache.\n",
			 ci->rec, ci->type);
		err++;
		continue;
	    }
	    ci2 = ci_ptr(v);

	    switch(ci->type) {
	    case GT_RecArray: {
		Array a1 = (Array)&ci->data;
		Array a2 = (Array)&ci2->data;

		if (a1->size != a2->size)
		    mis++;
		if (a1->max != a2->max) {
		    mis++;
		} else {
		    if (memcmp(a1->base, a2->base, a1->size * a1->max))
			mis++;
		}
		break;
	    }

	    case GT_Bin: {
		bin_index_t *b1 = (bin_index_t *)&ci->data;
		bin_index_t *b2 = (bin_index_t *)&ci2->data;

		if (b1->rec         != b2->rec ||
		    b1->pos         != b2->pos ||
		    b1->size        != b2->size ||
		    b1->start_used  != b2->start_used ||
		    b1->end_used    != b2->end_used ||
		    b1->parent_type != b2->parent_type ||
		    b1->parent      != b2->parent ||
		    b1->child[0]    != b2->child[0] ||
		    b1->child[1]    != b2->child[1] ||
		    b1->rng_rec     != b2->rng_rec ||
		    b1->flags       != b2->flags ||
		    b1->track       != b2->track ||
		    b1->track_rec   != b2->track_rec ||
		    b1->nseqs       != b2->nseqs ||
		    b1->rng_free    != b2->rng_free ||
		    b1->nrefpos     != b2->nrefpos ||
		    b1->nanno       != b2->nanno) {
		    mis++;
		} else if (b1->rng && b2->rng) {
		    if (ArrayMax(b1->rng) != ArrayMax(b2->rng)) {
			mis++;
		    } else {
			for (j = 0; j < ArrayMax(b1->rng); j++) {
			    range_t *r1 = arrp(range_t, b1->rng, j);
			    range_t *r2 = arrp(range_t, b2->rng, j);
			    if ((r1->flags & GRANGE_FLAG_UNUSED) !=
				(r2->flags & GRANGE_FLAG_UNUSED)) {
				mis++;
			    } else if ((r1->flags & GRANGE_FLAG_UNUSED) == 0) {
				if (r1->start    != r2->start ||
				    r1->end      != r2->end ||
				    r1->mqual    != r2->mqual ||
				    r1->rec      != r2->rec ||
				    r1->pair_rec != r2->pair_rec ||
				    r1->flags    != r2->flags)
				    mis++;
			    }
			}
		    }
		} else if ((b1->rng && ArrayMax(b1->rng) && !b2->rng) ||
			   (b2->rng && ArrayMax(b2->rng) && !b1->rng)) {
		    mis++;
		}
		break;
	    }

	    case GT_BTree: {
		break;
	    }

	    case GT_Database: {
		database_t *d1 = (database_t *)&ci->data;
		database_t *d2 = (database_t *)&ci2->data;
		if (d1->version != d2->version ||
		    d1->Ncontigs != d2->Ncontigs ||
		    d1->contig_order != d2->contig_order ||
		    d1->Nlibraries != d2->Nlibraries ||
		    d1->library != d2->library ||
		    d1->seq_name_index != d2->seq_name_index ||
		    d1->contig_name_index != d2->contig_name_index)
		    mis++;
		break;
	    }

	    case GT_Library: {
		break;
	    }

	    case GT_Contig: {
		contig_t *c1 = (contig_t *)&ci->data;
		contig_t *c2 = (contig_t *)&ci2->data;
		if (c1->rec      != c2->rec ||
		    c1->start    != c2->start ||
		    c1->end      != c2->end) {
		    mis++;
		} else if (c1->name && c2->name && strcmp(c1->name,c2->name)) {
		    mis++;
		}
		break;
	    }

	    case GT_SeqBlock: {
		seq_block_t *b1 = (seq_block_t *)&ci->data;
		seq_block_t *b2 = (seq_block_t *)&ci2->data;
		seq_t *s1, *s2;

		for (j = 0; j < SEQ_BLOCK_SZ; j++) {
		    if ((b1->seq[j] == NULL) != (b2->seq[j] == NULL)) {
			mis++;
			continue;
		    }

		    if (!b1->seq[j])
			continue;

		    s1 = b1->seq[j];
		    s2 = b2->seq[j];

		    if (s1->len            != s2->len ||
			s1->bin            != s2->bin ||
			s1->bin_index      != s2->bin_index ||
			s1->left           != s2->left ||
			s1->right          != s2->right ||
			s1->parent_rec     != s2->parent_rec ||
			s1->parent_type    != s2->parent_type ||
			s1->rec            != s2->rec ||
			s1->seq_tech       != s2->seq_tech ||
			s1->flags          != s2->flags ||
			s1->format         != s2->format ||
			s1->mapping_qual   != s2->mapping_qual ||
			s1->name_len       != s2->name_len ||
			s1->trace_name_len != s2->trace_name_len ||
			s1->alignment_len  != s2->alignment_len ||
			s1->aux_len        != s2->aux_len ||
			s1->idx            != s2->idx) {
			mis++;
		    } else {
			if (s1->name && s2->name &&
			    memcmp(s1->name, s2->name, s1->name_len))
			    mis++;
			if (s1->trace_name && s2->trace_name &&
			    memcmp(s1->trace_name, s2->trace_name,
				   s1->trace_name_len))
			    mis++;
			if (s1->alignment && s2->alignment &&
			    memcmp(s1->alignment, s2->alignment,
				   s1->alignment_len))
			    mis++;
			if (s1->seq && s2->seq &&
			    memcmp(s1->seq, s2->seq, ABS(s1->len)))
			    mis++;
			if (s1->conf && s2->conf &&
			    memcmp(s1->conf, s2->conf, ABS(s1->len)))
			    mis++;
			if (s1->sam_aux && s2->sam_aux &&
			    memcmp(s1->sam_aux, s2->sam_aux, s1->aux_len))
			    mis++;
		    }
		}
		break;
	    }

	    case GT_AnnoEleBlock: {
		anno_ele_block_t *b1 = (anno_ele_block_t *)&ci->data;
		anno_ele_block_t *b2 = (anno_ele_block_t *)&ci2->data;
		anno_ele_t *a1, *a2;

		for (j = 0; j < SEQ_BLOCK_SZ; j++) {
		    if ((b1->ae[j] == NULL) != (b2->ae[j] == NULL)) {
			mis++;
			continue;
		    }

		    if (!b1->ae[j])
			continue;

		    a1 = b1->ae[j];
		    a2 = b2->ae[j];

		    if (a1->tag_type != a2->tag_type ||
			a1->rec      != a2->rec ||
			a1->bin      != a2->bin ||
			a1->obj_type != a2->obj_type ||
			a1->obj_rec  != a2->obj_rec ||
			a1->anno_rec != a2->anno_rec ||
			a1->idx      != a2->idx) {
			mis++;
		    } else {
			if (a1->comment && a2->comment) {
			    if (strcmp(a1->comment, a2->comment))
				mis++;
			} else if (a1->comment || a2->comment) {
			    mis++;
			}
		    }
		}
		break;
	    }

	    default: {
		vmessage("Rec %"PRIrec" of type %d mismatches\n",
			 ci->rec, ci->type);
	    }
	    }

	    if (mis) {
		vmessage("Rec %"PRIrec" type %d differs in-mem and on-disk\n",
			 ci->rec, ci->type);
		err++;
	    }
	}
    }

    close_db(ior);
    return err;
}


/*
 * Performs a thorough internal consistency check of all on disk data
 * structures. It's therefore quite slow, but can highlight algorithm
 * problems during development.
 *
 * Returns the number of errors found
 *         or 0 on success.
 */
int check_database(GapIO *io, int fix, int level) {
    database_t *db;
    ArrayStruct *contig_order, *library;
    int i;
    int err = 0, fixed = 0;
    HacheTable *hash = NULL;

    vfuncheader("Check Database");
    vmessage("--DB version: %d\n", io->db->version);

    /* Check cache matches disk */
    if (level > 1) {
	vmessage("--Checking in-memory cache against disk\n");
	err += check_cache(io);
    }

    /* Load low level db structs; already in GapIO, but do full reload */
    db = cache_search(io, GT_Database, 0);
    if (!db) {
	vmessage("Failed to read GT_Database record 0\n");
	return ++err;
    }
    cache_incr(io, db);


    /* Check contig_order: 1 per contig and no more than 1 */
    contig_order = cache_search(io, GT_RecArray, db->contig_order);
    if (!contig_order) {
	vmessage("Failed to read contig order array\n");
	cache_decr(io, db);
	return ++err;
    }
    cache_incr(io, contig_order);

    hash = HacheTableCreate(256, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);
    if (db->Ncontigs != ArrayMax(contig_order)) {
	vmessage("Contig order array is not the same size as db->Ncontigs\n");
	err++;
	if (fix) {
	    cache_rw(io, io->contig_order);
	    ArrayMax(io->contig_order) = io->db->Ncontigs;
	    ArrayMax(contig_order) = io->db->Ncontigs;
	    fixed++;
	}
    }

    for (i = 0; i < ArrayMax(contig_order); i++) {
	tg_rec crec = arr(tg_rec, contig_order, i);
	HacheData hd;
	int new;

	hd.i = 0;
	HacheTableAdd(hash, (char *)&crec, sizeof(crec), hd, &new);
	if (!new) {
	    vmessage("Contig %"PRIrec" occurs more than once in the "
		     "contig_order array\n", crec);
	    err++;
	}
    }
    HacheTableDestroy(hash, 0);
    cache_decr(io, contig_order);


    /* Also check contig btree - every contig should have name in index */
    /* TODO */


    /* Check libraries, keep hash for later checking too */
    library = cache_search(io, GT_RecArray, db->library);
    if (!library) {
	vmessage("Failed to read library array\n");
	cache_decr(io, db);
	return ++err;
    }
    cache_incr(io, library);


    hash = HacheTableCreate(256, HASH_POOL_ITEMS | HASH_DYNAMIC_SIZE);
    if (db->Nlibraries != ArrayMax(library)) {
	vmessage("library array is not the same size as db->Nlibraries\n");
	err++;
    }
    for (i = 0; i < ArrayMax(library); i++) {
	tg_rec lrec = arr(tg_rec, library, i);
	HacheData hd;
	int new;

	hd.i = 0;
	HacheTableAdd(hash, (char *)&lrec, sizeof(lrec), hd, &new);
	if (!new) {
	    vmessage("Library %"PRIrec" occurs more than once in the "
		     "library array\n", lrec);
	    err++;
	}
    }
    cache_decr(io, db);
    cache_decr(io, library);

    /* If we're fixing the DB, also migrate from v1 to v2 */
    if (fix && io->db->version == 1) {
	io->db = cache_rw(io, io->db);
	io->iface->vers(io->dbh, 2);
	fixed++;
    }

    /* Check each contig in turn */
    for (i = 0; i < ArrayMax(contig_order); i++) {
	tg_rec crec = arr(tg_rec, contig_order, i);
	vmessage("--Checking contig #%"PRIrec" (%d of %d)\n",
		 crec, i+1, (int)ArrayMax(contig_order));
	UpdateTextOutput();
	err += check_contig(io, crec, fix, level, hash, &fixed);
    }

    if (fix && io->db->version == 1)
	io->db->version = 2;

    HacheTableDestroy(hash, 0);

    vmessage("\n*** Total number of errors: %d ***\n", err);
    if (fix)
	vmessage("*** Attempted to fix:       %d ***\n", fixed);

    return err;
}
