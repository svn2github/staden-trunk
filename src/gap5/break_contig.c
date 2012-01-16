#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "break_contig.h"
#include "misc.h"
#include "dis_readings.h" /* bin_destroy_recurse() */

#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

/*
 * Recursive part of remove_empty_bins.
 * Takes bin record.
 * Removes the bin if it is empty and has no children.
 * Modifies *first to contain the first bin record with data.
 *
 * Returns 1 if removed
 *         0 if not.
 */
static int remove_empty_bins_r(GapIO *io, tg_rec brec, tg_rec *first) {
    bin_index_t *bin = cache_search(io, GT_Bin, brec);
    int i, empty[2]; /* Emptied or non-existant */
    int this_is_empty;
    tg_rec child[2], f[2];

    /* Check if this bin is empty */
    this_is_empty = bin_empty(bin);

    /* Temporary copies to avoid needing cache_incr */
    child[0] = bin->child[0];
    child[1] = bin->child[1];

    f[0] = f[1] = 0;
    empty[0] = child[0] ? remove_empty_bins_r(io, child[0], &f[0]) : 1;
    empty[1] = child[1] ? remove_empty_bins_r(io, child[1], &f[1]) : 1;

    /* Remove this bin if empty and children are too */
    if (empty[0] && empty[1] && this_is_empty) {
	gio_debug(io, 1, "Bin %"PRIrec": this & children are empty "
		  "/ non-existant\n", brec);
	cache_rec_deallocate(io, GT_Bin, brec);
	return 1;
    }


    /* If we removed a child bin but are keeping this, then fix links */
    if ((empty[0] && child[0]) || (empty[1] && child[1])) {
	bin = cache_search(io, GT_Bin, brec);
	bin = cache_rw(io, bin);
	if (empty[0]) {
	    bin->flags |= BIN_BIN_UPDATED;
	    bin->child[0] = 0;
	}
	if (empty[1]) {
	    bin->flags |= BIN_BIN_UPDATED;
	    bin->child[1] = 0;
	}
    }


    /* Track first useful bin */
    if (first && !*first) {
	if ((f[0] && f[1]) || !this_is_empty) {
	    *first = brec;
	} else if (f[0]) {
	    *first = f[0];
	} else if (f[1]) {
	    *first = f[1];
	}
    }

    return 0;
}

/*
 * Tidies up after break contig or disassemble readings, looking for now
 * redundant bins.
 *
 * This has the following functions (not all implemented yet!)
 * 1) If a contig is totally empty, remove the contig.
 * 2) If a bin is empty and all below it, remove the bin.
 * 3) If a bin is empty and all above it, remove parent bins and link
 *    contig to new root.
 */
static void remove_empty_bins(GapIO *io, tg_rec contig) {
    contig_t *c = cache_search(io, GT_Contig, contig);
    tg_rec first = 0;

    cache_incr(io, c);

    if (c->bin) {
	if (remove_empty_bins_r(io, c->bin, &first)) {
	    cache_decr(io, c);
	    contig_destroy(io, contig);
	    return;
	}

	if (first != c->bin) {
	    bin_index_t *bin;
	    tg_rec bp, br, cdummy;
	    int offset, comp;

	    /* Cut out the offending waste */
	    bin = cache_search(io, GT_Bin, first);
	    bin = cache_rw(io, bin);
	    bp = bin->parent;

	    // Find new bin offset
	    if (bin_get_orient(io, bin->rec))
		bin->flags |= BIN_COMPLEMENTED;
	    else
		bin->flags &= ~BIN_COMPLEMENTED;

	    bin_get_position(io, bin, &cdummy, &offset, &comp);
	    assert(cdummy == contig);
	    bin->pos = offset;
	    bin->parent = contig;
	    bin->parent_type = GT_Contig;
	    bin->flags |= BIN_BIN_UPDATED;

	    c = cache_rw(io, c);
	    br = c->bin;
	    c->bin = first;

	    bin = cache_search(io, GT_Bin, bp);
	    bin = cache_rw(io, bin);
	    if (bin->child[0] == first) bin->child[0] = 0;
	    if (bin->child[1] == first) bin->child[1] = 0;

	    /* Recursively remove the bin tree from old root, br */
	    bin_destroy_recurse(io, br);
	}
    }

    cache_decr(io, c);
}

/*
 * Compute the visible statr position of a contig. This isn't just the extents
 * of start_used / end_used in the bins as this can included invisible
 * data such as cached consensus sequences.
 *
 * 'from' is the starting point to search from in a contig iterator.
 * Specify CITER_CSTART if you don't know what this is, otherwise if the
 * contig has been edited (eg shrunk) then you may want to specify an older
 * start coord in order to correctly trim all data.
 */
int contig_visible_start(GapIO *io, tg_rec crec, int from) {
    rangec_t *r;
    contig_iterator *ci;
    int seq_start = 0, seq_clipped_start;
    contig_t *c = cache_search(io, GT_Contig, crec);

    cache_incr(io, c);
    consensus_valid_range(io, crec, &seq_clipped_start, NULL);

    ci = contig_iter_new_by_type(io, crec, 1, CITER_FIRST | CITER_ISTART,
				 CITER_CSTART, CITER_CEND,
				 GRANGE_FLAG_ISANY);
    if (!ci) {
	cache_decr(io, c);
	return c->start;
    }
    
    while (r = contig_iter_next(io, ci)) {
	/* Seq only */
	if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
	    continue;

	seq_start = r->start;
	break;
    }
    contig_iter_del(ci);

    /* Trim annotations to visible portion */
    ci = contig_iter_new_by_type(io, crec, 1, CITER_FIRST | CITER_ISTART,
				 from, CITER_CEND,
				 GRANGE_FLAG_ISANY);
    while (ci && (r = contig_iter_next(io, ci))) {
	if (r->start >= seq_clipped_start)
	    break;

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS) {
	    if (r->end < seq_clipped_start) {
		bin_remove_refpos(io, crec, r->end);
	    }
	    continue;
	} else if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO) {
	    continue;
	}

	/*
	 * Why did we also clip sequence tags (to seq_start) before?
	 * The only time this would have an effect is surely when the tags
	 * started off outside of the sequence to start with, ie the
	 * database was already inconsistent.
	 */
	if (r->flags & GRANGE_FLAG_TAG_SEQ) continue;

	if (r->end < seq_clipped_start) {
	    bin_remove_item(io, &c, GT_AnnoEle, r->rec);
	} else {
	    bin_index_t *bin;
	    anno_ele_t *a;
	    range_t R, *R_out;

	    bin_remove_item(io, &c, GT_AnnoEle, r->rec);
	
	    R.start    = seq_clipped_start;
	    R.end      = MIN(r->end, c->end);
	    R.rec      = r->rec;
	    R.mqual    = r->mqual;
	    R.pair_rec = r->pair_rec;
	    R.flags    = r->flags;

	    bin = bin_add_range(io, &c, &R, &R_out, NULL, 0);

	    /* With luck the bin & index into bin haven't changed */
	    cache_incr(io, bin);
	    a = cache_search(io, GT_AnnoEle, r->rec);
	    if (a->bin != bin->rec /*||
		a->idx != R_out - ArrayBase(range_t, bin->rng)*/) {
		a = cache_rw(io, a);
		a->bin = bin->rec;
		//a->bin_idx = R_out - ArrayBase(range_t, bin->rng);
	    }

	    cache_decr(io, bin);
	}
    }
    contig_iter_del(ci);

    cache_decr(io, c);
    return seq_start;
}

/*
 * Compute the visible end position of a contig. This isn't just the extents
 * of start_used / end_used in the bins as this can included invisible
 * data such as cached consensus sequences.
 */
int contig_visible_end(GapIO *io, tg_rec crec, int to) {
    rangec_t *r;
    contig_iterator *ci;
    int seq_end = 0, seq_clipped_end;
    contig_t *c = cache_search(io, GT_Contig, crec);

    cache_incr(io, c);
    consensus_valid_range(io, crec, NULL, &seq_clipped_end);

    ci = contig_iter_new_by_type(io, crec, 1, CITER_LAST | CITER_IEND,
				 CITER_CSTART, CITER_CEND,
				 GRANGE_FLAG_ISANY);
    if (!ci) {
	cache_decr(io, c);
	return c->end;
    }
    
    while (r = contig_iter_prev(io, ci)) {
	/* Seq only */
	if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
	    continue;

	seq_end = r->end;
	break;
    }
    contig_iter_del(ci);

    /* Trim annotations to visible/clipped portion */
    ci = contig_iter_new_by_type(io, crec, 1, CITER_LAST | CITER_IEND,
				 CITER_CSTART, to,
				 GRANGE_FLAG_ISANY);
    while (ci && (r = contig_iter_prev(io, ci))) {
	if (r->end <= seq_clipped_end)
	    break;

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS) {
	    if (r->start > seq_clipped_end) {
		bin_remove_refpos(io, crec, r->start);
	    }
	    continue;
	} else if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO) {
	    continue;
	}

	/*
	 * Why did we also clip sequence tags (to seq_end) before?
	 * The only time this would have an effect is surely when the tags
	 * started off outside of the sequence to start with, ie the
	 * database was already inconsistent.
	 */
	if (r->flags & GRANGE_FLAG_TAG_SEQ) continue;

	if (r->start > seq_clipped_end) {
	    bin_remove_item(io, &c, GT_AnnoEle, r->rec);
	    continue;
	} else {
	    bin_index_t *bin;
	    anno_ele_t *a;
	    range_t R, *R_out;

	    bin_remove_item(io, &c, GT_AnnoEle, r->rec);
	
	    R.start    = MAX(r->start, c->start);
	    R.end      = seq_clipped_end;
	    R.rec      = r->rec;
	    R.mqual    = r->mqual;
	    R.pair_rec = r->pair_rec;
	    R.flags    = r->flags;

	    bin = bin_add_range(io, &c, &R, &R_out, NULL, 0);

	    /* With luck the bin & index into bin haven't changed */
	    cache_incr(io, bin);
	    a = cache_search(io, GT_AnnoEle, r->rec);
	    if (a->bin != bin->rec/* ||
                a->bin_idx != R_out - ArrayBase(range_t, bin->rng)*/) {
		a = cache_rw(io, a);
		a->bin = bin->rec;
		//a->bin_idx = R_out - ArrayBase(range_t, bin->rng);
	    }

	    cache_decr(io, bin);
	}
    }
    contig_iter_del(ci);

    cache_decr(io, c);
    return seq_end;
}

static int break_contig_move_bin(GapIO *io, bin_index_t *bin, int bpos,
				 contig_t *cfrom, tg_rec pfrom,
				 contig_t *cto,   tg_rec pto,
				 int child_no) {
    /* Add to */
    if (pto == cto->rec) {
	/* Parent is a contig */
	if (bin->rec != cto->bin) {
	    cache_rec_deallocate(io, GT_Bin, cto->rec);
	}
	cto->bin = bin->rec;
	/* Maximum extents - tightened up later */
	cto->start = bpos;
	cto->end = bpos + bin->size;

	bin->parent = cto->rec;
	bin->parent_type = GT_Contig;
	bin->flags |= BIN_BIN_UPDATED;

    } else {
	/* Parent is a bin */
	bin_index_t *pb;

	if (!(pb = get_bin(io, pto)))
	    return -1;
	if (!(pb = cache_rw(io, pb)))
	    return -1;

	pb->child[child_no] = bin->rec;
	pb->flags |= BIN_BIN_UPDATED;

	bin->parent = pto;
	bin->parent_type = GT_Bin;
	bin->flags |= BIN_BIN_UPDATED;
    }

    /* Remove from: NB it may not exist? */
    if (pfrom == cfrom->rec) {
	/* Parent is a contig */
	if (cfrom->bin != bin->rec) {
	    fprintf(stderr, "pfrom incorrect\n");
	    return -1;
	}

	cfrom->bin = 0;
    } else if (pfrom > 0) {
	/* Parent is a bin */
	bin_index_t *pb;

	if (!(pb = get_bin(io, pfrom)))
	    return -1;
	if (!(pb = cache_rw(io, pb)))
	    return -1;

	if (pb->child[0] != bin->rec && pb->child[1] != bin->rec) {
	    fprintf(stderr, "pfrom incorrect\n");
	    return -1;
	}

	if (!(pb = cache_rw(io, pb)))
	    return -1;
	
	if (pb->child[0] == bin->rec)
	    pb->child[0] = 0;
	else
	    pb->child[1] = 0;
	pb->flags |= BIN_BIN_UPDATED;
    }

    return 0;
}

/*
 * Given ranges contained within a bin this makes sure that all sequences
 * referred to in these ranges have their parent listed as the new bin.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int break_contig_reparent_seqs(GapIO *io, bin_index_t *bin) {
    int i, nr = bin->rng ? ArrayMax(bin->rng) : 0;

    for (i = 0; i < nr; i++) {
	range_t *r = arrp(range_t, bin->rng, i);
	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;
	
	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS) {
	    /* No object associated with this range */
	    continue;
	} else if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
	    anno_ele_t *a = (anno_ele_t *)cache_search(io, GT_AnnoEle, r->rec);
	    if (a->bin != bin->rec) {
		a = cache_rw(io, a);
		a->bin = bin->rec;
	    }
	} else {
	    seq_t *seq = (seq_t *)cache_search(io, GT_Seq, r->rec);
	    if (seq->bin != bin->rec) {
		seq = cache_rw(io, seq);
		seq->bin = bin->rec;
		seq->bin_index = i;
	    }
	}
    }

    return 0;
}

/*
 * A recursive break contig function.
 * bin_num	The current bin being moved or split.
 * pos		The contig break point.
 * pos2/iL	The maximum right extent observed for the left hand contig
 * pos3/iR	the minimum left extent observed for the right hand contig
 * offset	The absolute positional offset of this bin in original contig
 * pleft	The parent bin/contig record num in the left new contig
 * pright	The parent bin/contig record num in the right new contig
 * child_no     0 or 1 - whether this bin is the left/right child of its parent
 *
 *                         :
 *                         :        vL     iL(2)
 * ........--------------..:.....   |      |
 *            ......-------+---------.......
 *                     ....:...-----------------...
 *                     |   :   |  ......--------------.....
 *                     |   :   |
 *                  (3)iR  P   vR
 *
 * P     = break point
 * iR/iL = extents of right/left contig using invisible (ie all) data.
 * vR/vL = extents of right/left contig using visible clipped data only.
 */
static int break_contig_recurse(GapIO *io, HacheTable *h,
				contig_t *cl, contig_t *cr, tg_rec bin_num,
				int pos, int pos2, int pos3, int offset,
				int level, tg_rec pleft, tg_rec pright,
				int child_no, int complement) {
    int i, j, k, l, f_a, f_b;
    tg_rec rbin;
    bin_index_t *bin = get_bin(io, bin_num), *bin_dup ;
    //int bin_min, bin_max;
    int nseqs, nrefpos, nanno;
    tg_rec opright; /* old pright, needed if we revert back */

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

    gio_debug(io, 1, "%*sBreak offset %d pos %d/%d/%d => test bin %"PRIrec
	      ": %d..%d\n",
	      level*4, "",
	      offset, pos3, pos, pos2, bin->rec,
	      NMIN(bin->start_used, bin->end_used),
	      NMAX(bin->start_used, bin->end_used));

    bin = cache_rw(io, bin);
    nseqs = bin->nseqs;
    bin->nseqs = 0;
    nrefpos = bin->nrefpos;
    bin->nrefpos = 0;
    nanno = bin->nanno;
    bin->nanno = 0;

    /* Invalidate any cached data */
    bin_invalidate_track(io, bin, TRACK_ALL);
    if (bin->flags & BIN_CONS_VALID) {
	bin->flags |= BIN_BIN_UPDATED;
	bin->flags &= ~BIN_CONS_VALID;
    }

    //bin_min = !bin_empty(bin) ? NMIN(bin->start_used,bin->end_used) : offset;
    //bin_max = !bin_empty(bin) ? NMAX(bin->start_used,bin->end_used) : offset;

    /*
     * Add to right parent if this bin is to the right of pos,
     * or if the used portion is to the right and we have no left child.
     *
     * FIXME: Not a valid assumption!
     * The used portion of a bin is not a placeholder for the used portion
     * of all the the children beneath it. Therefore if the used portion of
     * this bin is > pos (and we have no left child) it still doesn't mean
     * that the absolute positions of the used portion of the right child
     * won't be < pos.
     */
    if (offset > pos2 /*|| (bin_min >= pos && !bin->child[0])*/) {
	gio_debug(io, 1, "%*sADD_TO_RIGHT pl=%"PRIrec" pr=%"PRIrec"\n",
		  level*4, "", pleft, pright);

	if (0 != break_contig_move_bin(io, bin, offset,
				       cl, pleft, cr, pright, 
				       child_no)) {
	    cache_decr(io, bin);
	    return -1;
	}

	bin_incr_nseq(io, bin, nseqs);
	bin_incr_nrefpos(io, bin, nrefpos);
	bin_incr_nanno(io, bin, nanno);
	cache_decr(io, bin);

	return 0;
    }

    /*
     * Add to left parent if this bin is entirely to the left of pos,
     * or if the used portion is to the left and we have no right child.
     */
    if (offset + bin->size < pos3 /*|| (bin_max < pos && !bin->child[1])*/) {
	gio_debug(io, 1, "%*sADD_TO_LEFT\n", level*4, "");

	//if (0 != break_contig_move_bin(io, bin, cr, pright, cl, pleft, child_no))
	//return -1;

	bin_incr_nseq(io, bin, nseqs);
	bin_incr_nrefpos(io, bin, nrefpos);
	bin_incr_nanno(io, bin, nanno);
	cache_decr(io, bin);
	
	return 0;
    }

    /*
     * Nominally the bin overlaps both left and right and so needs duplicating.
     * There are cases though at the roots of our trees where duplicating is
     * unnecessary as it leads to empty bins at the root. In this case
     * we skip creating a duplicate for the right, or alternatively steal
     * the left root bin and use that instead.
     *
     * Similarly the range_t array will either be left where it is, moved to
     * the right contig, or split in half (creating a new one for the right).
     *
     * FIXED: always need this. Eg:
     *
     * |-------------empty--------------|
     * |----------------|---------------|
     * |--------|-------|--------|------|
     *             ^
     *             |
     *             break here
     *
     * In this case we need to duplicate the parent as it overlaps the left
     * bin, which may (or may not) have data that needs to end up in the right
     * hand contig. Just duplicate for now and free later on if needed.
     */
    if (1 /* always! */ || pright != cr->rec ||
	(!bin_empty(bin) && NMAX(bin->start_used, bin->end_used) >= pos)) {
	gio_debug(io, 2, "NMAX=%d >= %d\n",
		  NMAX(bin->start_used, bin->end_used), pos);

	rbin = 0;

	/* Possibly steal left contig's bin */
#if 0
	if (pleft == cl->rec && NMIN(bin->start_used, bin->end_used) >= pos2) {
	    /* Currently this doesn't always work */
	    if (bin->child[1]) {
		bin_index_t *ch = get_bin(io, bin->child[1]);
		if (NMIN(ch->pos, ch->pos + ch->size-1) >= pos) {
		    rbin = cl->bin;
		    cl->bin = bin->child[0];
		}
	    }
	    pleft = bin->rec;
	} else {
	    pleft = bin->rec;
	}
#else
	pleft = bin->rec;
#endif

	/* Create new bin, or use root of contig if it's unused so far */
	if (!rbin && pright == cr->rec) {
	    rbin = cr->bin;
	}

	/* Otherwise we genuingly need a duplicate */
	if (!rbin)
	    rbin = bin_new(io, 0, 0, 0, GT_Bin);

	/* Initialise with duplicate values from left bin */
	bin_dup = get_bin(io, rbin);
	bin_dup = cache_rw(io, bin_dup);
	bin_dup->size = bin->size;
	bin_dup->pos = bin->pos;
	bin_dup->parent = pright;
	bin_dup->parent_type = (pright == cr->rec ? GT_Contig : GT_Bin);
	bin_dup->flags = bin->flags | BIN_BIN_UPDATED;
	bin_dup->start_used = bin->start_used;
	bin_dup->end_used = bin->end_used;

	/*
	 * Shift bin to offset if it's the contig root.
	 * It'll be shifted back by the correct amount later.
	 */
	if (pright == cr->rec) {
	    gio_debug(io, 1, "moving root bin to offset=%d comp=%d\n",
		      offset, complement);
	    bin_dup->pos = offset;
	}

	gio_debug(io, 1, "%*sCreated dup for right, rec %"PRIrec"\n",
		  level*4,"", bin_dup->rec);
	break_contig_move_bin(io, bin_dup, offset,
			      cl, 0, cr, pright, child_no);
	opright = pright;
	pright = bin_dup->rec;
    } else {
	bin_dup = NULL;
	pleft = bin->rec;
    }

    if (bin_empty(bin)) {
	/* Empty bin */
	gio_debug(io, 1, "%*sEMPTY range\n", level*4, "");
	bin->start_used = bin->end_used = 0;
	bin->flags |= BIN_BIN_UPDATED;
	if (bin_dup) {
	    bin_dup->start_used = bin_dup->end_used = 0;
	    bin_dup->flags |= BIN_BIN_UPDATED;
	}
	    
    } else if (NMIN(bin->start_used, bin->end_used) > pos2) {
	/* Move range to right contig */
	gio_debug(io, 1, "%*sDUP %"PRIrec", MOVE Array to right\n",
		  level*4, "", bin_dup->rec);

	bin_dup->rng = bin->rng;
	bin_dup->rng_rec = bin->rng_rec;
	bin_dup->rng_free = bin->rng_free;
	bin_dup->flags |= BIN_BIN_UPDATED;
	if (bin_dup->rng_rec)
	    bin_dup->flags |= BIN_RANGE_UPDATED;

	if (bin->rec != bin_dup->rec) {
	    bin->rng = NULL;
	    bin->rng_rec = 0;
	    bin->rng_free = -1;
	    bin->flags |= BIN_BIN_UPDATED;
	}

	bin->start_used = bin->end_used = 0;
	break_contig_reparent_seqs(io, bin_dup);

	if (bin_dup->rng) {
	    int n = ArrayMax(bin_dup->rng);
	    for (i = j = k = l = 0; i < n; i++) {
		range_t *r = arrp(range_t, bin_dup->rng, i);
		if (r->flags & GRANGE_FLAG_UNUSED)
		    continue;

		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
		    HacheData hd; hd.i = 1;
		    HacheTableAdd(h, (char *)&r->rec, sizeof(r->rec), hd,NULL);
		    j++;
		}
		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS)
		    k++;
		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
		    l++;
	    }
	    bin_incr_nseq(io, bin_dup, j);
	    bin_incr_nrefpos(io, bin_dup, k);
	    bin_incr_nanno(io, bin_dup, l);
	}
    } else if (NMAX(bin->start_used, bin->end_used) < pos3) {
	/* Range array already in left contig, so do nothing */
	gio_debug(io, 1, "%*sMOVE Array to left\n", level*4, "");

	if (bin_dup)
	    bin_dup->start_used = bin_dup->end_used = 0;

	if (bin->rng) {
	    int n = ArrayMax(bin->rng);
	    for (i = j = k = l = 0; i < n; i++) {
		range_t *r = arrp(range_t, bin->rng, i);
		if (r->flags & GRANGE_FLAG_UNUSED)
		    continue;

		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
		    HacheData hd; hd.i = 0;
		    HacheTableAdd(h, (char *)&r->rec, sizeof(r->rec), hd,NULL);
		    j++;
		}
		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS)
		    k++;
		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
		    l++;
	    }
	    bin_incr_nseq(io, bin, j);
	    bin_incr_nrefpos(io, bin, k);
	    bin_incr_nanno(io, bin, l);
	}
    } else {
	/*
	 * Range array covers pos3 to pos2, so split in two.
	 * Optimisation (to do): maybe it doesn't overlap pos though, in
	 * which case we can move the bin, but sift out affected tags.
	 */
	int n;
	int nsl = 0, nsr = 0; /* no. seqs */
	int nrl = 0, nrr = 0; /* no. refpos */
	int nal = 0, nar = 0; /* no. anno */
	int lmin = bin->size, lmax = 0, rmin = bin->size, rmax = 0;

	gio_debug(io, 1, "%*sDUP %"PRIrec", SPLIT array\n",
		  level*4, "", bin_dup->rec);

	bin->flags |= BIN_RANGE_UPDATED;
	bin_dup->flags |= BIN_RANGE_UPDATED;

	bin_dup->rng = ArrayCreate(sizeof(range_t), 0);
	bin_dup->rng_free = -1;

	/* Pass 1 - hash sequences */
	n = ArrayMax(bin->rng);
	for (i = 0; i < n; i++) {
	    range_t *r = arrp(range_t, bin->rng, i);
	    int cstart; /* clipped sequence positions */
	    seq_t *s;

	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;

	    if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
		continue;

	    if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS)
		continue;

	    s = (seq_t *)cache_search(io, GT_Seq, r->rec);
	    if (!s) {
		verror(ERR_WARN, "break_contig", "failed to load seq #%"PRIrec,
		       r->rec);
		continue;
	    }
	    if ((s->len < 0) ^ complement) {
		cstart = NMAX(r->start, r->end) - (s->right-1);
	    } else {
		cstart = NMIN(r->start, r->end) + s->left-1;
	    }
	    
	    if (cstart >= pos)  {
		HacheData hd; hd.i = 1;
		HacheTableAdd(h, (char *)&r->rec, sizeof(r->rec), hd, NULL);
		//printf("Add seq #%"PRIrec" to hash value 1\n", r->rec);
	    } else {
		int end;

		HacheData hd; hd.i = 0;
		HacheTableAdd(h, (char *)&r->rec, sizeof(r->rec), hd, NULL);
		//printf("Add seq #%"PRIrec" to hash value 0\n", r->rec);

		end = NMAX(r->start, r->end);
	    }
	}
	
	/* Pass 2 - do the moving of anno/seqs */
	n = ArrayMax(bin->rng);
	for (i = j = 0; i < n; i++) {
	    range_t *r = arrp(range_t, bin->rng, i), *r2;
	    int cstart, cend; /* clipped sequence positions */

	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;

	    if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
		anno_ele_t *a;

		/* Contig tags moved entirely based on their position */
		a = cache_search(io, GT_AnnoEle, r->rec);
		if (a->obj_type == GT_Contig) {
		    cstart = NMIN(r->start, r->end);
		    cend   = NMAX(r->start, r->end);

		    /*
		     * Ideally we should track the visible potion of
		     * sequences and move consensus tags to wherever they're
		     * in visible sequence. This may be < or > pos depending
		     * on the layout. However for now we take a simpler
		     * approach.
		     */
		    if (cend <= pos2 || cstart < pos3) {
			cstart = pos-1;
		    } else {
			cstart = pos+1;
		    }

		/* Seq tags get moved only if the sequence itself moves */
		} else if (a->obj_type == GT_Seq) {
		    HacheItem *hi;
		    hi = HacheTableSearch(h, (char *)&r->pair_rec,
					  sizeof(r->pair_rec));

		    if (hi) {
			/*
			 * Fake cstart to force move or no-move based on
			 * whether the sequence we belong to is moving
			 * or not.
			 */
			//printf("Found %"PRIrec" val %d\n", r->pair_rec, hi->data.i);
			if (hi->data.i == 0)
			    cstart = pos-1;
			else
			    cstart = pos+1;
		    } else {
			/*
			 * We apparently have a tag in a bin higher up
			 * than the sequence it belongs to. In theory tags
			 * are smaller than their sequences so this could
			 * never happen. However joining contigs can create
			 * overlapping bins and in turn this means new tags
			 * can be placed down a different bin hierarchy
			 * than the seqs. 
			 * In this case query that bin and populate our
			 * hash from that bin too.
			 */
			seq_t *s;
			tg_rec s_ctg;
			int s_start, s_end, s_orient, c_start;
			HacheData hd;

			sequence_get_position(io, r->pair_rec, &s_ctg,
					      &s_start, &s_end, &s_orient);
			s = cache_search(io, GT_Seq, r->pair_rec);

			if ((s->len < 0) ^ s_orient) {
			    c_start = s_end - (s->right-1);
			} else {
			    c_start = s_start + s->left-1;
			}

			if (c_start >= pos) {
			    cstart = pos+1;
			    hd.i = 1;
			} else {
			    cstart = pos-1;
			    hd.i = 0;
			}

			//printf("Add %"PRIrec" val %d\n", s->rec, c_start >= pos);

			/*
			 * Add to hache too so multiple tags on one seq
			 * only need one sequence_get_position call.
			 */
			HacheTableAdd(h, (char *)&s->rec, sizeof(s->rec),
				      hd, NULL);
		    }
		} else {
		    verror(ERR_WARN, "break_contig",
			   "Unknown obj_type for anno_ele #%"PRIrec"\n",
			   r->rec);
		    /* Best guess */
		    cstart = NMIN(r->start, r->end);
		}

	    } else if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS) {
		cstart = NORM(r->start);
	    } else {
		seq_t *s = (seq_t *)cache_search(io, GT_Seq, r->rec);
		if ((s->len < 0) ^ complement) {
		    cstart = NMAX(r->start, r->end) - (s->right-1);
		} else {
		    cstart = NMIN(r->start, r->end) + s->left-1;
		}
	    }

	    if (cstart >= pos) {
		r2 = (range_t *)ArrayRef(bin_dup->rng, ArrayMax(bin_dup->rng));
		*r2 = *r;
		if (rmin > r->start) rmin = r->start;
		if (rmin > r->end)   rmin = r->end;
		if (rmax < r->start) rmax = r->start;
		if (rmax < r->end)   rmax = r->end;
		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ)
		    nsr++;
		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS)
		    nrr++;
		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
		    nar++;

		//printf("Move rec %"PRIrec" bin %"PRIrec" idx %d\n", r->rec, bin->rec, i);

		/* Mark as unused in old bin */
		r->flags = GRANGE_FLAG_UNUSED;
		r->rec = bin->rng_free;
		bin->rng_free = i;
		bin->flags |= BIN_BIN_UPDATED | BIN_RANGE_UPDATED;
	    } else {
		if (lmin > r->start) lmin = r->start;
		if (lmin > r->end)   lmin = r->end;
		if (lmax < r->start) lmax = r->start;
		if (lmax < r->end)   lmax = r->end;
		/*
		if (j != i) {
		    r2 = arrp(range_t, bin->rng, j);
		    *r2 = *r;
		}
		j++;
		*/
		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ)
		    nsl++;
		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS)
		    nrl++;
		if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
		    nal++;
	    }
	}
	bin_incr_nseq(io, bin, nsl);
	bin_incr_nseq(io, bin_dup, nsr);
	bin_incr_nrefpos(io, bin, nrl);
	bin_incr_nrefpos(io, bin_dup, nrr);
	bin_incr_nanno(io, bin, nal);
	bin_incr_nanno(io, bin_dup, nar);
	
	//ArrayMax(bin->rng) = j;
	bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;

#if 0
	/*
	 * Right now this causes problems, but I'm not sure why. Try again
	 * after we've fixed the bin->nseqs issues and other deallocation
	 * woes.
	 */

	if (ArrayMax(bin_dup->rng) == 0 && bin_dup->parent_type == GT_Bin) {
	    /* We didn't need it afterall! Odd. */
	    bin_index_t *pb;

	    gio_debug(io, 1, "Purging bin %d that we didn't need afterall\n",
		      bin_dup->rec);
	    cache_rec_deallocate(io, GT_Bin, bin_dup->rec);
	    pb = cache_search(io, GT_Bin, bin_dup->parent);
	    if (pb->child[0] == bin_dup->rec)
		pb->child[0] = 0;
	    if (pb->child[1] == bin_dup->rec)
		pb->child[1] = 0;
	    bin_dup = NULL;
	    pright = opright;
	}
#endif

	if (bin_dup)
	    break_contig_reparent_seqs(io, bin_dup);

	if (lmax < lmin) {
	    /* No data left in bin */
	    bin->start_used = 0;
	    bin->end_used = 0;
	} else {
	    bin->start_used     = lmin;
	    bin->end_used       = lmax;
	}

	gio_debug(io, 1, "%*sLeft=>%d..%d right=>%d..%d\n", level*4, "",
		  lmin, lmax, rmin, rmax);

	if (bin_dup) {
	    if (rmin <= rmax) {
		bin_dup->start_used = rmin;
		bin_dup->end_used   = rmax;
	    } else {
		/* No data moved in bin */
		bin_dup->start_used = 0;
		bin_dup->end_used   = 0;
	    }
	}
    }

    /* Recurse */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;

	ch = get_bin(io, bin->child[i]);
	if (0 != break_contig_recurse(io, h, cl, cr, bin->child[i],
				      pos, pos2, pos3,
				      NMIN(ch->pos, ch->pos + ch->size-1),
				      level+1, pleft, pright,
				      i, complement)) {
	    cache_decr(io, bin);
	    return -1; 
	}
    }

    cache_decr(io, bin);
    //    if (bin_dup)
    //	cache_decr(io, bin_dup);

    return 0;
}

/*
 * Looks for redundant bins at the root containing no data and just a single
 * child.
 *
 * FIXME: We need to compensate for bin position here. Hence this function
 * is not called for now.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int remove_redundant_bins(GapIO *io, contig_t *c) {
    tg_rec bnum;

    if (!(c = cache_rw(io, c)))
	return -1;

    for (bnum = c->bin; bnum;) {
	bin_index_t *bin = get_bin(io, bnum);
	if (!bin_empty(bin) || (bin->child[0] && bin->child[1]))
	    break;

	/* Empty */
	c->bin = bin->child[0] ? bin->child[0] : bin->child[1];
	gio_debug(io, 1, "Remove bin %"PRIrec"\n", bin->rec);
	bnum = c->bin;
    }

    return 0;
}

/*
 * Given a left and right contig (cl,cr) with an overlap such that
 * cl ends at left_end and cr starts at right_start, we need to duplicate
 * any ISREFPOS type sequences from cr into cl.
 *
 * ie:            |              |left_end
 *                |              v
 * cl ------------|---------------
 *            +   |  x     x    x           x   x      x
 * cr             |      ---------------------------------
 *                |      ^
 *       break pos|      |right start
 *
 * The ISREFPOS markers (x) will have been moved from cl to cr.
 * (Actually right start will be the left-most 'x' and not where the sequence
 * data starts, but this is something we need to fix too.)
 */
int copy_isrefpos_markers(GapIO *io, contig_t *cl, contig_t *cr,
			  int right_start, int left_end) {
    contig_iterator *ci;
    rangec_t *rc;
    int first_seq = left_end;

    gio_debug(io, 1, "Moving ISREFPOS markers from contig %"PRIrec
	      " (%d..%d) to contig %"PRIrec".\n",
	      cl->rec, right_start, left_end, cr->rec);

    ci = contig_iter_new_by_type(io, cr->rec, 0, CITER_FIRST,
				 right_start, left_end,
				 GRANGE_FLAG_ISANY);
    if (!ci)
	return right_start;

    while (rc = contig_iter_next(io, ci)) {
	range_t r;

	if ((rc->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	    if (first_seq > rc->start)
		first_seq = rc->start;
	}

	if ((rc->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISREFPOS)
	    continue;

	if (rc->start < first_seq) {
	    bin_index_t *bin;
	    range_t *r2;

	    gio_debug(io, 1, "** Deleting from cr, bin %"PRIrec" **\n",
		      rc->orig_rec);
	    
	    bin = cache_search(io, GT_Bin, rc->orig_rec);
	    bin = cache_rw(io, bin);

	    r2 = arrp(range_t, bin->rng, rc->orig_ind);
	    assert(r2->mqual == rc->mqual);
	    assert(r2->flags == rc->flags);

	    gio_debug(io, 1, "Mark %d for removal\n", rc->orig_ind);
	    r2->flags = GRANGE_FLAG_UNUSED;
	    r2->rec = bin->rng_free;
	    bin->rng_free = rc->orig_ind;

	    bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;
	    bin_incr_nrefpos(io, bin, -1);
	}

	r.start    = rc->start;
	r.end      = rc->end;
	r.rec      = rc->rec;
	r.pair_rec = rc->pair_rec;
	r.mqual    = rc->mqual;
	r.flags    = rc->flags;
	
	bin_add_range(io, &cl, &r, NULL, NULL, 1);
    }
    bin_add_range(io, NULL, NULL, NULL, NULL, -1);

    gio_debug(io, 1, "First real seq in cr = %d\n", first_seq);

    contig_iter_del(ci);

    return first_seq;
}

/*
 * Checks if the break can do anything. Eg breaking at the very start or
 * very end could yield a zero-read contig, which causes inconsistencies.
 * We just abort if so.
 *
 * Returns 0 for OK (and modifies cpos to pos of next seq),
 *        -1 to abort/error
 */
int break_check_counts(GapIO *io, tg_rec crec, int *cpos_p) {
    contig_iterator *ci;
    rangec_t *r;
    int left_ok = 0, right_ok = 0;
    int cpos = *cpos_p, new_cpos;

    /* Check left end */
    ci = contig_iter_new(io, crec, 1, CITER_LAST | CITER_ISTART,
			 CITER_CSTART, cpos-1);
    if (!ci)
	return -1;

    while (r = contig_iter_prev(io, ci)) {
	seq_t *s = cache_search(io, GT_Seq, r->rec);
	int cstart;
	if (!s)
	    return -1;

	if ((s->len < 0) ^ r->comp) {
	    cstart = r->start + ABS(s->len) - (s->right-1) - 1;
	} else {
	    cstart = r->start + s->left-1;
	}

	if (cstart < cpos) {
	    left_ok = 1;
	    break;
	}
    }

    if (!left_ok) {
	contig_iter_del(ci);
	return -1;
    }


    /* Check right end */
    ci = contig_iter_new(io, crec, 1, CITER_FIRST | CITER_ISTART,
			 cpos-1, CITER_CEND);
    if (!ci)
	return -1;

    new_cpos = INT_MAX;
    while (r = contig_iter_next(io, ci)) {
	seq_t *s;
	int cstart;

	if (new_cpos != INT_MAX && r->start >= new_cpos)
	    break;

	s = cache_search(io, GT_Seq, r->rec);
	if (!s)
	    return -1;

	if ((s->len < 0) ^ r->comp) {
	    cstart = r->start + ABS(s->len) - (s->right-1) - 1;
	} else {
	    cstart = r->start + s->left-1;
	}

	if (cstart >= cpos && new_cpos > cstart)
	    new_cpos = cstart;

	if (cstart >= cpos) {
	    right_ok = 1;
	    //break;
	}
    }
    
    *cpos_p = new_cpos;

    if (!right_ok) {
	contig_iter_del(ci);
	return -1;
    }

    return 0;
}

/*
 * Computes the rightmost extent of sequences with visible start < cpos.
 * (That's total extent, not just clipped data.)
 *
 * We can't rely on doing this as we go along as we may have overlapping bins
 * due to contig joining.
 */
int compute_pos2(GapIO *io, tg_rec crec, int cpos) {
    contig_iterator *ci;
    rangec_t *r;
    int pos2 = cpos;

    ci = contig_iter_new_by_type(io, crec, 0, CITER_FIRST | CITER_ISTART,
				 cpos, CITER_CEND, GRANGE_FLAG_ISSEQ);
    if (!ci) {
	verror(ERR_WARN, "break_contig", "Failed to create contig iterator");
	return cpos;
    }

    while (r = contig_iter_next(io, ci)) {
	if (r->start >= cpos)
	    break;

	if (pos2 < r->end) {
	    seq_t *s = cache_search(io, GT_Seq, r->rec);
	    int start;

	    if ((s->len < 0) ^ r->comp)
		start = r->start + ABS(s->len) - (s->right-1) - 1;
	    else
		start = r->start + s->left - 1;
	
	    if (start < cpos)
		pos2 = r->end;
	}
    }    
    contig_iter_del(ci);

    return pos2;
}

/*
 * As above but computes the first base of seqs with visible start >= cpos.
 */
int compute_pos3(GapIO *io, tg_rec crec, int cpos) {
    contig_iterator *ci;
    rangec_t *r;
    int pos3 = cpos;

    ci = contig_iter_new_by_type(io, crec, 0, CITER_LAST | CITER_IEND,
				 CITER_CSTART, cpos, GRANGE_FLAG_ISSEQ);
    if (!ci) {
	verror(ERR_WARN, "break_contig", "Failed to create contig iterator");
	return cpos;
    }

    while (r = contig_iter_prev(io, ci)) {
	if (r->end < cpos)
	    break;

	if (pos3 > r->start) {
	    seq_t *s = cache_search(io, GT_Seq, r->rec);
	    int start;

	    if ((s->len < 0) ^ r->comp)
		start = r->start + ABS(s->len) - (s->right-1) - 1;
	    else
		start = r->start + s->left - 1;
	
	    if (start >= cpos)
		pos3 = r->start;
	}
    }    
    contig_iter_del(ci);

    return pos3;
}

/*
 * Breaks a contig in two such that snum is the right-most reading of
 * a new contig.
 */
tg_rec break_contig(GapIO *io, tg_rec crec, int cpos, int break_holes) {
    contig_t *cl;
    contig_t *cr;
    int cid;
    char cname[1024], *cname_end;
    int left_end, right_start;
    bin_index_t *bin;
    int do_comp = 0;
    HacheTable *h;
    tg_rec ret = -1;

    /*
     * Check we have at least 1 seq in each half. Also move cpos to be the
     * first base of the right hand contig (and on the right hand edge of
     * any hole if cpos is within one).
     */
    gio_debug(io, 1, "Moved break point from %d ", cpos);
    if (break_check_counts(io, crec, &cpos) == -1) {
	verror(ERR_WARN, "break_contig",
	       "Breaking at %d would create a contig with no sequences. Abort",
	       cpos);
	return -1;
    }
    gio_debug(io, 1, "to %d\n", cpos);

    cl = (contig_t *)cache_search(io, GT_Contig, crec);
    cache_incr(io, cl);

    //contig_dump_ps(io, &cl, "/tmp/tree.ps");

    /*
     * Our hash table is keyed on sequence record numbers for all sequences
     * in all bins spanning the break point. The value is either 0 or 1
     * for left/right contig.
     * 
     * The purpose of this hash is to allow us to work out whether a tag
     * belongs in the left or right contig, as a tag could start beyond the
     * break point but be attached to a sequence before the break point.
     *
     * Further complicating this is that a tag could be in a smaller bin
     * than the sequence as it may not be as long. However we know
     * we'll recurse down these in a logical order so we can be sure
     * we've already "seen" the sequence that the tag has been
     * attached to.
     */
    h = HacheTableCreate(1024, HASH_DYNAMIC_SIZE);

    strncpy(cname, contig_get_name(&cl), 1000);
    cname_end = cname + strlen(cname);
    cid = 1;
    do {
	sprintf(cname_end, "#%d", cid++);
    } while (contig_index_query(io, cname) > 0);

    if (!(cr = contig_new(io, cname))) {
	cache_decr(io, cl);
	verror(ERR_WARN, "break_contig",
	       "Failed to create a new contig with name %s", cname);
	return -1;
    }
    cl = cache_rw(io, cl);
    cr = cache_rw(io, cr);
    gio_debug(io, 1, "Break in contig %"PRIrec", pos %d\n", crec, cpos);

    gio_debug(io, 1, "Existing left bin = %"PRIrec", right bin = %"PRIrec"\n",
	      cl->bin, cr->bin);

    cache_incr(io, cr);

    bin = get_bin(io, cl->bin);
    do_comp = bin->flags & BIN_COMPLEMENTED;

    break_contig_recurse(io, h, cl, cr,
			 contig_get_bin(&cl), cpos,
			 compute_pos2(io, cl->rec, cpos),
			 compute_pos3(io, cl->rec, cpos),
			 contig_offset(io, &cl), 0, cl->rec, cr->rec, 0, 0);

    /* Recompute end positions */
    left_end    = contig_visible_end(io, cl->rec, CITER_CEND);
    right_start = contig_visible_start(io, cr->rec, CITER_CSTART);

    /* Also do other ends, simply to tidy up end tags */
    contig_visible_start(io, cl->rec, CITER_CSTART);
    contig_visible_end(io, cr->rec, CITER_CEND);

    //if (cl->bin) contig_dump_ps(io, &cl, "/tmp/tree_l.ps");
    //if (cr->bin) contig_dump_ps(io, &cr, "/tmp/tree_r.ps");

    /* Duplicate overlapping ISREFPOS markers between right_start & left_end */
    right_start = copy_isrefpos_markers(io, cl, cr, right_start, left_end);

    /* Ensure start/end positions of contigs work out */
    bin = cache_rw(io, get_bin(io, cr->bin));

    //#define KEEP_POSITIONS 1
#ifndef KEEP_POSITIONS
    cr->start = 1;
    cr->end = cl->end - right_start + 1;
    bin->pos -= right_start-1;

    /*
     * END may not be cl->end - right_start + 1 if the previous right hand
     * end was a long read that started before the break point.
     * Brute force and find it. The current cr->end is just a maximum possible
     * extent to make consensus_unclipped_range() efficient.
     */
    consensus_unclipped_range(io, cr->rec, NULL, &cr->end);
#else
    cr->start = right_start;
    cr->end = cl->end;
#endif

    consensus_unclipped_range(io, cl->rec, &cl->start, NULL);

    if ((do_comp && !(bin->flags & BIN_COMPLEMENTED)) ||
	(!do_comp && (bin->flags & BIN_COMPLEMENTED))) {
	bin->flags ^= BIN_COMPLEMENTED;
    }

    cl->end = left_end;

    //    remove_redundant_bins(io, cl);
    //    remove_redundant_bins(io, cr);

    gio_debug(io, 1, "Final left bin = %"PRIrec", right bin = %"PRIrec"\n",
	      cl->bin, cr->bin);

    HacheTableDestroy(h, 0);

    //if (cl->bin) contig_dump_ps(io, &cl, "/tmp/tree_l.ps");
    //if (cr->bin) contig_dump_ps(io, &cr, "/tmp/tree_r.ps");

    cache_flush(io); /* Needed before dellocation? */
    
    remove_empty_bins(io, cl->rec);
    remove_empty_bins(io, cr->rec);

    /* Empty contig? If so remove it completely */
    if (cl->bin == 0) {
	gio_debug(io, 1, "Removing empty contig %"PRIrec"\n", cl->rec);
	contig_destroy(io, cl->rec);
    }
    if (cr->bin == 0) {
	gio_debug(io, 1, "Removing empty contig %"PRIrec"\n", cr->rec);
	contig_destroy(io, cr->rec);
    }

    cache_flush(io);

    /* Check for holes */
    if (break_holes) {
	int st, en;
#ifndef KEEP_POSITIONS
	st = 1;
	en = left_end-right_start+1;
#else
	st = right_start;
	en = left_end;
#endif
	if (0 != remove_contig_holes(io, cr->rec, st, en, 0)) {
	    cache_decr(io, cl);
	    cache_decr(io, cr);

	    verror(ERR_WARN, "break_contig",
		   "Failue in remove_contig_holes(io, cr->rec, %d, %d, 0)",
		   st, en);
	    return -1;
	}
    }

    ret = cr->rec;

    cache_decr(io, cl);
    cache_decr(io, cr);

    return ret;
}
