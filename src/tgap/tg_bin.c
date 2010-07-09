#include <assert.h>
#include <string.h>
#include <math.h>

#include "xalloc.h"
#include "tg_gio.h"
#include "tg_tracks.h"

#define get_bin(io, bnum) ((bin_index_t *)cache_search((io), GT_Bin, (bnum)))

/*
 * Allocates a new bin record.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
#if 0
int bin_new(GapIO *io, int pos, int sz, int parent, int parent_type) {
    int rec;
    bin_index_t bin;

    /* Initialise disk struct */
    bin.pos         = pos;
    bin.size        = sz;
    bin.start_used  = 0;
    bin.end_used    = 0;
    bin.parent      = parent;
    bin.parent_type = parent_type;
    bin.child[0]    = 0;
    bin.child[1]    = 0;
    bin.bin_id      = 0;
    bin.rng	    = NULL;
    bin.rng_rec     = 0;
    bin.flags       = BIN_BIN_UPDATED;
    bin.track       = NULL;
    bin.track_rec   = 0;
    bin.nseqs       = 0;
    bin.rng_free    = -1;

    if (-1 == (rec = io->iface->bin.create(io->dbh, &bin)))
	return -1;

    return rec;
}
#endif

int bin_new(GapIO *io, int pos, int sz, int parent, int parent_type) {
    int rec;
    bin_index_t *bin;

    if (-1 == (rec = io->iface->bin.create(io->dbh, NULL)))
	return -1;

    /* Initialise disk struct */
    bin = get_bin(io, rec);
    bin = cache_rw(io, bin);
    bin->pos         = pos;
    bin->size        = sz;
    bin->start_used  = 0;
    bin->end_used    = 0;
    bin->parent      = parent;
    bin->parent_type = parent_type;
    bin->child[0]    = 0;
    bin->child[1]    = 0;
    bin->bin_id      = 0;
    bin->rng	    = NULL;
    bin->rng_rec     = 0;
    bin->flags       = BIN_BIN_UPDATED;
    bin->track       = NULL;
    bin->track_rec   = 0;
    bin->nseqs       = 0;
    bin->rng_free    = -1;

    return rec;
}


/*
 * Doubles up the number of bins by adding a new root node and duplicating
 * the tree.
 *
 * It takes the old root_id as an argument and returns the new one.
 * Returns -1 on failure.
 */
static bin_index_t *contig_extend_bins_right(GapIO *io, contig_t **c) {
    int old_root_id = contig_get_bin(c);
    bin_index_t *oroot = get_bin(io, old_root_id), *nroot;
    int root_id;
    int sz = oroot->size;

    cache_incr(io, oroot);
    if (!(oroot = cache_rw(io, oroot))) {
	cache_decr(io, oroot);
	return NULL;
    }

    /* Create a new root */
    root_id = bin_new(io, oroot->pos, sz*2, oroot->parent, oroot->parent_type);
    nroot = get_bin(io, root_id);
    cache_incr(io, nroot);
    if (!(nroot = cache_rw(io, nroot))) {
	cache_decr(io, oroot);
	cache_decr(io, nroot);
	return NULL;
    }

    nroot->child[0] = old_root_id;
    nroot->child[1] = 0;
    nroot->nseqs = oroot->nseqs;

    nroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, nroot);

    /* Move old left bin */
    oroot->parent = root_id;
    oroot->parent_type = GT_Bin;
    oroot->pos = 0;

    oroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, oroot);

    contig_set_bin(io, c, root_id);

    return nroot;
}

static bin_index_t *contig_extend_bins_left(GapIO *io, contig_t **c) {
    int old_root_id = contig_get_bin(c);
    bin_index_t *oroot = get_bin(io, old_root_id), *nroot;
    int root_id;
    int sz = oroot->size;

    cache_incr(io, oroot);
    if (!(oroot = cache_rw(io, oroot))) {
	cache_decr(io, oroot);
	return NULL;
    }

    /* Create a new root */
    root_id = bin_new(io, oroot->pos-sz, sz*2, oroot->parent, oroot->parent_type);
    nroot = get_bin(io, root_id);
    cache_incr(io, nroot);
    if (!(nroot = cache_rw(io, nroot))) {
	cache_decr(io, oroot);
	cache_decr(io, nroot);
	return NULL;
    }

    nroot->child[0] = 0;
    nroot->child[1] = old_root_id;
    nroot->nseqs = oroot->nseqs;

    nroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, nroot);

    /* Move old right bin */
    oroot->parent = root_id;
    oroot->parent_type = GT_Bin;
    oroot->pos = sz;

    oroot->flags |= BIN_BIN_UPDATED;
    cache_decr(io, oroot);

    contig_set_bin(io, c, root_id);

    return nroot;
}

/*
 * Allocates and returns next range from a range array
 * Returns -1 for error.
 */
static int next_range(bin_index_t *bin) {
    Array ra = bin->rng;

    if (bin->rng_free == -1) {
	if (NULL == ArrayRef(ra, ArrayMax(ra)))
	    return -1;

	return ArrayMax(ra)-1;

    } else {
	int tmp;
	range_t *r = arrp(range_t, bin->rng, bin->rng_free);
	
	tmp = bin->rng_free;
	bin->rng_free = r->rec;
	return tmp;
    }
}

/*
 * This finds a bin suitable for a given range. If such a bin doesn't
 * exist it can optionally create one if the 'extend' flag is true.
 *
 * When a bin is found the absolute offset of that bin is returned
 * in 'offset_r' (may be NULL).
 *
 * Returns the bin pointer on success
 *         NULL on failure
 */
#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

#define CACHE_LAST_BIN
bin_index_t *bin_for_range(GapIO *io, contig_t **c,
			   int start, int end, int extend,
			   int *offset_r,  int *comp_r) {
    int offset;
    bin_index_t *bin = get_bin(io, contig_get_bin(c));
    int complement = 0;
    int f_a, f_b;

#ifdef CACHE_LAST_BIN
    static int last_c = 0;
    static GapIO *last_io = NULL;
    bin_index_t *last_bin = NULL;
    static int last_bin_rec = 0;
    static int last_start = 0, last_end = 0;
    static int last_offset, last_complement;
#endif

    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }

    offset = bin->pos;
    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    //cache_incr(io, bin);

    /*
     * If we're trying to insert a new node beyond the total bounds of the
     * root node then we need to extend the bin structure first to either
     * the left and/or the right.
     */
    while (end >= bin->pos + bin->size) {
	//cache_decr(io, bin);
	if (extend)
	    bin = contig_extend_bins_right(io, c);
	else
	    return NULL;
	//cache_incr(io, bin);

	complement = 0;
	/*
	 * But root has switched from complemented to uncomplemented, so the
	 * range may need fixing too.
	 * See tg_tcl.c for an example way of detecting and fixing this.
	 */
    }

    while (start < bin->pos) {
	//cache_decr(io, bin);
	if (extend)
	    bin = contig_extend_bins_left(io, c);
	else
	    return NULL;
	//cache_incr(io, bin);

	complement = 0;
    }

    /*
     * In theory we can jump straight to a candidate starting bin, possibly
     * even returning it right here if it's the min bin size, saving about
     * 10% of our CPU time in this function
     */
#ifdef CACHE_LAST_BIN
    if (last_bin_rec && last_c == (*c)->rec && last_io == io) {
	last_bin = cache_search(io, GT_Bin, last_bin_rec);
	if (start >= last_start && end <= last_end) {
	    if (last_bin && last_bin->size <= io->min_bin_size &&
		!last_bin->child[0] && !last_bin->child[1]) {
		/* leaf node, so we can return right now */
		if (offset_r)
		    *offset_r = last_offset;
		if (comp_r)
		    *comp_r = complement;
		return last_bin;
	    }

	    /* Maybe a smaller bin, but start the search from here on */
	    bin = last_bin;
	    offset = last_offset;
	    complement = last_complement;
	    //cache_incr(io, bin);

	    if (complement) {
		f_a = -1;
		f_b = offset + bin->size-1;
	    } else {
		f_a = +1;
		f_b = offset;
	    }

	    goto jump;
	}
    }
#endif

    /* Now recurse down the bin hierachy searching for the smallest bin */
    offset = bin->pos;
 jump:
    for (;;) {
	int i;
	bin_index_t *ch;

	if (complement) {
	    f_a = -1;
	    f_b = offset + bin->size-1;
	} else {
	    f_a = +1;
	    f_b = offset;
	}

	/* Find which child bin is most suitable */
	cache_incr(io, bin);
	for (i = 0; i < 2;) {
	    if (bin->child[i] <= 0) {
		i++;
		continue;
	    }

	    ch = get_bin(io, bin->child[i]);

	    //	    if (start >= offset + ch->pos &&
	    //		end   <= offset + ch->pos + ch->size-1) {
	    if (start >= NMIN(ch->pos, ch->pos + ch->size-1) &&
		end   <= NMAX(ch->pos, ch->pos + ch->size-1)) {
		cache_decr(io, bin);
		bin = ch;
		cache_incr(io, ch);
		offset = NMIN(ch->pos, ch->pos + ch->size-1);

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

		i = 0; /* restart loop */
	    } else {
		i++;
	    }
	}

	if (!extend) {
	    if (offset_r)
		*offset_r = offset;
	    if (comp_r)
		*comp_r = complement;
	    return bin;
	}

	/*
	 * We now have the smallest bin available holding this range, but 
	 * as we delay creation of the sub-bins until required then we
	 * should perhaps create a child bin.
	 */
	if (bin->size > io->min_bin_size &&
	    (!bin->child[0] || !bin->child[1])) {
	    /* Construct left child if needed */
	    if (!bin->child[0]) {
		int pos;
		int sz;

		if (complement) {
		    pos = bin->size/2;
		    
		    if (bin->child[1]) {
			ch = get_bin(io, bin->child[1]);
			pos = ch->size;
		    }
		    sz = bin->size - pos;
		} else {
		    pos = 0;
		    sz = bin->size/2;
		    if (bin->child[1]) {
			ch = get_bin(io, bin->child[1]);
			sz = ch->pos;
		    }
		}

		//		if (start >= offset + pos &&
		//		    end   <= offset + pos + sz-1) {
		if (start >= NMIN(pos, pos+sz-1) &&
		    end   <= NMAX(pos, pos+sz-1)) {
		    /* It would fit - create it and continue recursing */

		    if (!(bin = cache_rw(io, bin)))
			return NULL;

		    bin->child[0] = bin_new(io, pos, sz, bin->rec, GT_Bin);
		    bin->flags |= BIN_BIN_UPDATED;
		    cache_decr(io, bin);

		    bin = get_bin(io, bin->child[0]);
		    offset += bin->pos;
		    continue;
		}
	    }

	    /* Construct right child if needed */
	    if (!bin->child[1]) {
		int pos;
		int sz;

		if (complement) {
		    pos = 0;
		    sz = bin->size/2;
		    if (bin->child[0]) {
			ch = get_bin(io, bin->child[0]);
			sz = ch->pos;
		    }
		} else {
		    pos = bin->size/2;
		    if (bin->child[0]) {
			ch = get_bin(io, bin->child[0]);
			pos = ch->size;
		    }
		    sz = bin->size - pos;
		}

		//		if (start >= offset + pos &&
		//		    end   <= offset + pos + sz-1) {
		if (start >= NMIN(pos, pos+sz-1) &&
		    end   <= NMAX(pos, pos+sz-1)) {
		    /* It would fit - create it and continue recursing */

		    if (!(bin = cache_rw(io, bin)))
			return NULL;

		    bin->child[1] = bin_new(io, pos, sz, bin->rec, GT_Bin);
		    bin->flags |= BIN_BIN_UPDATED;
		    cache_decr(io, bin);

		    bin = get_bin(io, bin->child[1]);
		    offset += bin->pos;
		    continue;
		}
	    }
	}
	
	/* At the smallest bin already */
	break;
    }

    if (offset_r)
	*offset_r = offset;
    if (comp_r)
	*comp_r = complement;

#ifdef CACHE_LAST_BIN
    last_io = io;
    last_c = (*c)->rec;
    last_bin_rec = bin->rec;
    last_start = offset;
    last_end = offset + bin->size-1;
    last_offset = offset;
    last_complement = complement;
#endif

    cache_decr(io, bin);
    return bin;
}

/*
 * Adds 'n' to the nseq counter for a bin and all parent bins chaining up
 * to the root node. 'n' may be negative too.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_incr_nseq(GapIO *io, bin_index_t *bin, int n) {
    while (bin) {
	if (!(bin = cache_rw(io, bin)))
	    return -1;
	bin->nseqs += n;
	bin->flags |= BIN_BIN_UPDATED;

	if (bin->parent_type != GT_Bin)
	    break;

	bin = get_bin(io, bin->parent);
    }

    return 0;
}

/*
 * Adds a range to the contig.
 *
 * Returns the bin we added the range to on success
 *         NULL on failure
 */
bin_index_t *bin_add_range(GapIO *io, contig_t **c, range_t *r,
			   range_t **r_out, int *complemented) {

    bin_index_t *bin;
    range_t *r2;
    int nr, offset, comp;

    if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	if (contig_get_start(c) == contig_get_end(c)) {
	    contig_set_start(io, c, r->start);
	    contig_set_end(io, c, r->end);
	}
	if (contig_get_start(c) > r->start)
	    contig_set_start(io, c, r->start);

	if (contig_get_end(c) < r->end)
	    contig_set_end(io, c, r->end);
    }

    if (!(bin = bin_for_range(io, c, r->start, r->end, 1, &offset,
			      &comp))) //complemented)))
	return NULL;

    if (complemented)
	*complemented = comp;

    if (!(bin = cache_rw(io, bin)))
	return NULL;

    /* Adjust start/end used in bin */
    if (bin->rng) {
	if (bin->start_used > r->start - offset)
	    bin->start_used = r->start - offset;
	if (bin->end_used < r->end - offset)
	    bin->end_used = r->end - offset;
    } else {
	/* Initial case */
	bin->start_used = r->start - offset;
	bin->end_used = r->end - offset;
    }

    /* Update Range array */
    bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;
    bin->flags &= ~BIN_CONS_VALID;
    if (!bin->rng)
	bin->rng = ArrayCreate(sizeof(range_t), 0);

    nr = next_range(bin);
    r2 = arrp(range_t, bin->rng, nr);

    *r2 = *r; /* struct copy */
    r2->start -= offset;
    r2->end -= offset;

    if (comp) {
	int tmp = r2->start;
	r2->start = bin->size-1 - r2->end;
	r2->end   = bin->size-1 - tmp;
    }

    if (r_out)
	*r_out = r2;

    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ)
	bin_incr_nseq(io, bin, 1);

    return bin;
}

/*
 * Finds the contig number and position of a record number.
 *
 * If non-NULL r_out is filled with the associated range_t struct.
 *
 * If non-NULL i_out is filled with a pointer to the object referred to
 * by type/rec. (This is just for minor optimisations.) If returned it
 * will have had cache_incr() run on it, so the caller should
 * use cache_decr() to permit deallocation.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_get_item_position(GapIO *io, int type, GRec rec,
			  int *contig, int *start, int *end, int *orient,
			  int *brec, range_t *r_out, void **i_out) {
    bin_index_t *bin;
    int bnum, i;
    int offset1 = 0, offset2 = 0, found = 0;
    int comp = 0;

    if (type == GT_AnnoEle) {
	anno_ele_t *a = cache_search(io, GT_AnnoEle, rec);
	if (i_out) {
	    cache_incr(io, a);
	    *i_out = a;
	}
	bnum = a->bin;
    } else if (type == GT_Seq) {
	seq_t *s = (seq_t *)cache_search(io, GT_Seq, rec);

	if (i_out) {
	    cache_incr(io, s);
	    *i_out = s;
	}
	bnum = s->bin;
    } else {
	fprintf(stderr, "Unsupported record type %d in bin_get_item_position\n",
		type);
    }
    
    if (brec)
	*brec = bnum;

    /* Find the position of this anno within the bin */
    bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
        range_t *r = arrp(range_t, bin->rng, i);
	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if (r->rec == rec) {
	    found = 1;
	    offset1 = r->start;
	    offset2 = r->end;

	    if (r_out)
		*r_out = *r;
	    break;
	}
    }

    if (!found) {
	if (i_out)
	    *i_out = NULL;
	return -1;
    }

    /* Find the position of this bin relative to the contig itself */
    for (;;) {
	if (bin->flags & BIN_COMPLEMENTED) {
	    offset1 = bin->size-1 - offset1;
	    offset2 = bin->size-1 - offset2;
	    comp ^= 1;
	}
	offset1 += bin->pos;
	offset2 += bin->pos;

	if (bin->parent_type != GT_Bin)
	    break;

	bnum = bin->parent;
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
    }

    assert(bin->parent_type == GT_Contig);

    if (contig)
	*contig = bin->parent;
    if (start)
	*start = offset1 < offset2 ? offset1 : offset2;
    if (end)
	*end = offset1 > offset2 ? offset1 : offset2;
    if (orient)
	*orient = comp;

    return 0;
}

/*
 * Removes a record referenced from a known bin
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_item_from_bin(GapIO *io, contig_t **c, bin_index_t **binp,
			     int type, int rec) {
    bin_index_t *bin;
    int i;

    if (!(bin = cache_rw(io, *binp)))
	return -1;
    *binp = bin;

    /* FIXME: we should check and update bin->start_used and bin->end_used */

    /* FIXME: use seq->bin_index or anno_ele->idx here as a short-cut? */
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);
	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if (r->rec != rec)
	    continue;

	r->flags |= GRANGE_FLAG_UNUSED;
	r->rec = bin->rng_free;
	bin->rng_free = i;
	bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ)
	    bin_incr_nseq(io, bin, -1);

	break;
    }

    return 0;
}

/*
 * Removes a record referenced from a bin. As above but we only know the
 * record and not which bin it's part of
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_remove_item(GapIO *io, contig_t **c, int type, int rec) {
    bin_index_t *bin;
    int start, end;
    int cnum, bnum;

    if (-1 == bin_get_item_position(io, type, rec, &cnum, &start, &end,
				    NULL, &bnum, NULL, NULL))
	return -1;


    //if (!(bin = bin_for_range(io, c, start, end, 0, NULL, NULL)))
    //    return -1;
    bin = cache_search(io, GT_Bin, bnum);

    return bin_remove_item_from_bin(io, c, &bin, type, rec);
}


/*
 * Finds the contig number and position of the start of a bin.
 * The right position is obviously this + bin->size (regardless of
 * whether it has been complemented).
 *
 * Returns 0 on success (plus *contig & *pos)
 *        -1 on failure
 */
int bin_get_position(GapIO *io, bin_index_t *bin, int *contig, int *pos) {
    int bnum;
    int offset = 0;

    /* FIXME: Complemented coordinates need more fixes here */
    if (bin->flags & BIN_COMPLEMENTED) {
	offset = bin->size-1 - offset;
    }

    offset += bin->pos;

    /* Find the position of this bin relative to the contig itself */
    while (bin->parent_type == GT_Bin) {
	bnum = bin->parent;
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
	offset += bin->pos;
    }
    
    assert(bin->parent_type == GT_Contig);
    *contig = bin->parent;

    *pos = offset;

    return 0;
}

/*
 * Some tracks are "official" cached objects and are deallocated as part of
 * the tg_cache system. Others are temporary on-the-fly structs generated
 * for the purpose of temporary display. These we need to free.
 *
 * This function knows which is which and does the appropriate thing.
 */
void track_free(track_t *t) {
    if (t->flag & TRACK_FLAG_FREEME) {
	if (t->data)
	    ArrayDestroy(t->data);
	free(t);
    }
}

/*
 * Creates a fake track struct, to be freed with track_free
 *
 * Returns track_t on success
 *         NULL on failure
 */
track_t *track_create_fake(int type, int size) {
    track_t *t = (track_t *)calloc(1, sizeof(*t));
    if (!t)
	return 0;
    t->type = type;
    t->nitems = size;
    t->item_size = sizeof(int);
    t->data = ArrayCreate(sizeof(int), size);
    t->flag |= TRACK_FLAG_FREEME;

    return t;
}

/*
 * Creates a track of a given type for this bin.
 * Note, this does not actually add it to the bin (but probably should
 * otherwise it's nothing more than the non-bin track_create).
 *
 * Returns track_t pointer on success
 *         NULL on failure
 */
track_t *bin_create_track(GapIO *io, bin_index_t *bin, int type) {
    int rec;
    track_t *t;

    if (-1 == (rec = io->iface->track.create(io->dbh, NULL)))
	return NULL;

    t = (track_t *)cache_search(io, GT_Track, rec);
    track_set_type(io, &t, type);

    return t;
}

/*
 * Adds a track to this bin.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bin_add_track(GapIO *io, bin_index_t **bin, track_t *track) {
    bin_index_t *n;
    int i;
    bin_track_t *bt;

    if (!(n = cache_rw(io, *bin)))
	return -1;
    *bin = n;

    /* Create new bin-track, or error if already found */
    if (!n->track) {
	n->track = ArrayCreate(sizeof(bin_track_t), 0);
	n->flags |= BIN_TRACK_UPDATED;
    }

    for (i = 0; i < ArrayMax(n->track); i++) {
	bt = arrp(bin_track_t, n->track, i);
	if (bt->type == track->type)
	    return -1;
    }

    /* Add the track pointer */
    bt = (bin_track_t *)ArrayRef(n->track, ArrayMax(n->track));
    bt->type  = track->type;
    bt->flags = 1;
    bt->rec   = track->rec;
    bt->track = track;

    return 0;
}

/*
 * Finds the track of a given type for this bin.
 *
 * Returns track_t pointer on success (do not free)
 *         NULL on failure (eg bin too small)
 */
track_t *bin_get_track(GapIO *io, bin_index_t *bin, int type) {
    int i;

    /* If it exists and is up to date, return it */
    if (bin->track) {
	for (i = 0; i < ArrayMax(bin->track); i++) {
	    bin_track_t *bt = arrp(bin_track_t, bin->track, i);
	    if (bt->type == type) {
		if (bt->track)
		    return bt->track;
		
		return (track_t *)cache_search(io, GT_Track, bt->rec);
	    }
	}
    }

    return NULL;
}

/*
 * A bit like bin_get_track, but this is designed to auto-generate and
 * update the track as desired. The expectation is that this will always
 * succeed and anything else is a fatal error.
 */
track_t *bin_query_track(GapIO *io, bin_index_t *bin, int type) {
    track_t *bin_recalculate_track(GapIO *io, bin_index_t *bin, int type);
    int i;

    /* If it exists and is up to date, return it */
    if (bin->track) {
	for (i = 0; i < ArrayMax(bin->track); i++) {
	    bin_track_t *bt = arrp(bin_track_t, bin->track, i);
	    if (bt->type == type && (bt->flags & TRACK_FLAG_VALID))
		return (track_t *)cache_search(io, GT_Track, bt->rec);
	}
    }

    /* Otherwise generate and maybe cache */
    return bin_recalculate_track(io, bin, type);
}

/*
 * Invalidates a track.
 * Returns 0 on success
 *        -1 on failure.
 */
int bin_invalidate_track(GapIO *io, bin_index_t *bin, int type) {
    int i;

    if (!bin->track)
	return 0;

    for (i = 0; i < ArrayMax(bin->track); i++) {
	bin_track_t *bt = arrp(bin_track_t, bin->track, i);
	if (bt->type == type || type == TRACK_ALL) {
	    if (NULL == (bin = cache_rw(io, bin)))
		return -1;
	    printf("Update track for rec %d\n", bin->rec);
	    bin->flags |= BIN_TRACK_UPDATED;
	    bt = arrp(bin_track_t, bin->track, i);
	    bt->flags &= ~TRACK_FLAG_VALID;
	}
    }

    return 0;
}

/*
 * ---------------------------------------------------------------------------
 * Track handling
 */

track_t *bin_recalculate_track(GapIO *io, bin_index_t *bin, int type) {
    int pos, cnum;
    track_t *track, *child;
    int nele;
    double bpv;
    contig_t *c;

    /*
     * So we have a bin of a given size in which we wish to have at least
     * RD_ELEMENTS of track samples, but it's out of date. We query the
     * contig for track data at a resolution of double what we need and
     * then average/downsample it to generate our new bin stats.
     */
    bpv = (double)bin->size / RD_ELEMENTS;
    if (bpv < 1) bpv = 1;
    nele = bin->size / bpv;
    if (nele & 1) nele++;
    bpv = (double)bin->size / nele;

    /*
     * Bottom layer, so no point querying a child - we just create it
     * ourselves now.
     */
    if (bpv <= 2) { /* FIXME: was 1, need 1 to stop aliasing? */
	track_t *fake;
	int rec, *depth;
	fake = track_create_fake(type, bin->size);
	fake->flag = TRACK_FLAG_FREEME;

	/* FIXME: type specific code here - or the size at least (int) */
	switch (type) {
	case TRACK_READ_DEPTH:
	    depth = track_read_depth_r1(io, bin);
	    break;
	default:
	    fprintf(stderr, "Unknown track type %d\n", type);
	    return NULL;
	}
	memcpy(ArrayBase(int, fake->data), depth, bin->size * sizeof(int));
	free(depth);

	rec = io->iface->track.create(io->dbh, fake);
	track = (track_t *)cache_search(io, GT_Track, rec);

	printf("Initialising track %d flag %d in bin %d at bpv of 1\n",
	       rec, track->flag, bin->rec);

	bin_add_track(io, &bin, track);
	track_free(fake);

	return track;
    }

    /* Else use child bin tracks, in ever-decreasing circles */
    if (-1 == bin_get_position(io, bin, &cnum, &pos))
	return NULL;
    c = (contig_t *)cache_search(io, GT_Contig, cnum);
    child = contig_get_track(io, &c, pos, pos + bin->size-1, type, bpv);
    if (NULL == child)
	return NULL;

    track = bin_get_track(io, bin, type);
    if (!track) {
	track = bin_create_track(io, bin, type);
	bin_add_track(io, &bin, track);
    }

    /* Copy child 'fake track' to our real track */
    track_set_data(io, &track, ArrayCreate(sizeof(int), nele));
    track_set_nitems(io, &track, nele);
    track_set_item_size(io, &track, sizeof(int));
    memcpy(ArrayBase(int, track->data), ArrayBase(int, child->data),
	   nele * sizeof(int));

    track_free(child);
    cache_decr(io, track);

    return track;
}
