#include <staden_config.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "tg_gio.h"
#include "misc.h"
#include "tree.h"
#include "dna_utils.h"
#include "consensus.h"
#include "break_contig.h"
#include "tg_check.h"
#include "list_proc.h"

#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

typedef struct bin_list {
  tg_rec rec;
  int parent_comp;
  int parent_offset;
  int parent_size;
  int parent_level;
  int child_index;
  struct bin_list *next;
} bin_list;

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

    //n->timestamp = io_timestamp_incr(io);

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

    //n->timestamp = io_timestamp_incr(io);

    return 0;
}

/*
 * Sets the contig bin number
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_set_bin(GapIO *io, contig_t **c, tg_rec value) {
    contig_t *n;
    if (!(n = cache_rw(io, *c)))
	return -1;

    n->bin = value;
    *c = n;

    // Up to caller to do if required.
    //n->timestamp = io_timestamp_incr(io);

    return 0;
}

/*
 * Sets a contig name.
 * If we try to do this on a child I/O, eg the contig editor, then we
 * only change the name local to this structure.
 * At the time of committing the child I/O to the parent I/O the global
 * name index will then be updated too as part of the flushing process.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_set_name(GapIO *io, contig_t **c, char *name) {
    contig_t *n;
    GapIO *iob = gio_base(io);

    if (!(n = cache_rw(io, *c)))
	return -1;

    /* Delete old name */
    if (n->name && !io->base) {
	tg_rec r = iob->iface->contig.index_del(iob->dbh, n->name, n->rec);
	if (r != -1 && r != io->db->contig_name_index) {
	    io->db = cache_rw(io, io->db);
	    io->db->contig_name_index = r;
	}
    }

    if (NULL == (n = cache_item_resize(n, sizeof(*n) + strlen(name)+1)))
	return -1;

    *c = n;

    /* Add new name */
    n->name   = (char *)(&n->data);
    strcpy(n->name, name);

    if (*name && !io->base) {
	tg_rec r = iob->iface->contig.index_add(iob->dbh, name, n->rec);
	if (r != -1 && r != io->db->contig_name_index) {
	    io->db = cache_rw(io, io->db);
	    io->db->contig_name_index = r;
	}
    }

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
 * Note this code uses some bizarre hybrid of coordinate systems. The original
 * code manipulates offset, comp and pos such that pos changes to be the
 * position relative to that bin. This works fine for working out what should
 * be moved or inserted to. This is passed around in pos / offset / comp.
 *
 * However the hash table for monitoring sequence ranges and whether
 * annotations need moving has to span bins, so it needs an absolute
 * position and offset computed against that absolute pos. These get passed
 * around in apos / aoffset.
 *
 * Returns >=0 for success (0 => no insertion made, 1 => insert done)
 *          -1 for failure
 */
static int contig_insert_base2(GapIO *io, tg_rec crec, tg_rec bnum,
			       int pos, int apos, int start_of_contig,
			       int offset, int aoffset, char base, int conf,
			       int nbases, int comp, HacheTable *hash) {
    int i, ins = 0;
    bin_index_t *bin;
    HacheData hd;
    int f_a, f_b;

    bin = get_bin(io, bnum);
    cache_incr(io, bin);

    if (!(bin = cache_rw(io, bin)))
	return -1;

    /* Normalise pos for complemented bins */
    if (bin->flags & BIN_COMPLEMENTED) {
	comp ^= 1;
	pos = offset + bin->size - pos;
    } else {
	pos -= offset;
    }

    if (comp) {
	f_a = -1;
	f_b = aoffset + bin->size-1;
    } else {
	f_a = +1;
	f_b = aoffset;
    }

    /* FIXME: add end_used (or start_used if complemented?) check here,
     * so we can shortcut looping through the bin if we're inserting beyond
     * any of the bin contents.
     */

    /* Perform the insert into objects first */
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);

	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if ((start_of_contig &&
	     pos <= MAX(r->start, r->end)+1 &&
	     pos >= MIN(r->start, r->end))  ||
	    (!start_of_contig &&
	     pos <= MAX(r->start, r->end)   &&
	     pos >  MIN(r->start, r->end))) {
	    //printf("pos overlap obj #%"PRIrec" %d in %d..%d\n",
	    //	   r->rec, pos,
	    //	   MIN(r->start, r->end),
	    //	   MAX(r->start, r->end));
	    /* Insert */
	    if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISREFPOS &&
		(r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO) {
		/* ISCONS? skip perhaps as we invalidate them anyway */
		seq_t *s = cache_search(io, GT_Seq, r->rec);
		int no_ins;

		if (!s) {
		    verror(ERR_WARN, "contig_insert_base2",
			   "failed to load seq #%"PRIrec, r->rec);
		    continue;
		}

		if (ABS(r->end - r->start) + 1 != ABS(s->len)) {
		    verror(ERR_WARN, "contig_insert_base2", 
			   "Range start/end are inconsistent with seq len. ");
		}

		no_ins = 0;
		if (base) {
		    int r_start, r_end;
		    if (/*comp ^ */(s->len < 0)) {
			r_start = MIN(r->end - (s->right-1),
				      r->end - (s->left-1));
			r_end = MAX(r->end - (s->right-1),
				    r->end - (s->left-1));
		    } else {
			r_start = MIN(r->start + s->left-1,
				      r->start + s->right-1);
			r_end = MAX(r->start + s->right-1,
				    r->start + s->left-1);
		    }
		    if (pos <= r_start || pos > r_end)
			no_ins = 1;

		    // printf("Rec %"PRIrec" visible %d..%d abs %d..%d\n",
		    //	   r->rec, r_start, r_end,
		    //	   NORM(r_start), NORM(r_end));

		    if (!no_ins) {
			//printf("INS %"PRIrec" at %d\n", r->rec,
			//       pos - MIN(r->start, r->end));
			//int i;
			//for (i = 0; i < nbases; i++)
			//    sequence_insert_base(io, &s,
			//			 pos - MIN(r->start, r->end),
			//			 base, conf, 0);
			sequence_insert_bases(io, &s,
					      pos - MIN(r->start, r->end),
					      base, conf, nbases, 0);
			if (hash) {
			    //printf("1 Mov %"PRIrec"\n", r->rec);
			    hd.i = MAX(NMIN(r_start, r_end), apos);
			    HacheTableAdd(hash, (char *)&r->rec,
					  sizeof(r->rec), hd, NULL);
			}
			//printf("1 %"PRIrec"->end++\n", r->rec);
			r->end+=nbases;
			ins = 1;
		    }
		}

		if (!base || no_ins) {
		    /* Shift instead, but only if in left-cutoff data */ 
		    if (/*comp ^ */ (s->len < 0)) {
			if (MIN(r->end - (s->right-1), r->end - (s->left-1))
			    >= pos) {
			    if (hash) {
				//printf("2 Mov %"PRIrec"\n", r->rec);
				hd.i = comp ? INT_MAX : INT_MIN;
				HacheTableAdd(hash, (char *)&r->rec,
					      sizeof(r->rec), hd, NULL);
			    }
			    //printf("2 %"PRIrec"->start/end++\n", r->rec);
			    r->start+=nbases;
			    r->end+=nbases;
			    ins = 1;
			}
		    } else {
			if (MIN(r->start + s->left-1, r->start + s->right-1)
			    >= pos) {
			    if (hash) {
				//printf("3 Mov %"PRIrec"\n", r->rec);
				hd.i = comp ? INT_MAX : INT_MIN;
				HacheTableAdd(hash, (char *)&r->rec,
					      sizeof(r->rec), hd, NULL);
			    }
			    //printf("3 %"PRIrec"->start/end++\n", r->rec);
			    r->start+=nbases;
			    r->end+=nbases;
			    ins = 1;
			}
		    }
		}
	    } else if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
		if (base && !(r->flags & GRANGE_FLAG_TAG_SEQ)) {
		    //printf("grow anno %"PRIrec"\n", r->rec);
		    r->end+=nbases;
		    ins = 1;
		}
	    }
	    
	    bin->flags |= BIN_RANGE_UPDATED;

	} else if (MIN(r->start, r->end) >= pos) {
	    //printf("pos to left of obj #%"PRIrec" %d %d..%d\n", 
	    //	   r->rec, pos,
	    //	   MIN(r->start, r->end),
	    //	   MAX(r->start, r->end));
	    if ( (r->flags&GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO) {
		/* Move */
		//printf("4 %"PRIrec"->start/end++\n", r->rec);
		r->start+=nbases;
		r->end+=nbases;
		ins = 1;
		bin->flags |= BIN_RANGE_UPDATED;

		if (hash) {
		    //printf("4 Mov %"PRIrec"\n", r->rec);
		    hd.i = comp ? INT_MAX : INT_MIN;
		    HacheTableAdd(hash, (char *)&r->rec,
				  sizeof(r->rec), hd, NULL);
		}
	    } else if (!(r->flags & GRANGE_FLAG_TAG_SEQ)) {
		/* Consensus tag */
		r->start+=nbases;
		r->end+=nbases;
		ins = 1;
		bin->flags |= BIN_RANGE_UPDATED;
	    }
	} /* else pos to right of object */
    }

    /* Adjust the bin dimensions */
    {
	bin->size+=nbases;
	ins = 1;
	if (bin->rng && ArrayMax(bin->rng)) {
	    int start = INT_MAX;
	    int end   = INT_MIN;
	    for (i = 0; i < ArrayMax(bin->rng); i++) {
		range_t *r = arrp(range_t, bin->rng, i);
	    
		if (r->flags & GRANGE_FLAG_UNUSED)
		    continue;

		if (start > r->start)
		    start = r->start;
		if (end   < r->end)
		    end   = r->end;
	    }
	    if (start != INT_MAX) {
		bin->start_used = start;
		bin->end_used = end;
	    } else {
		bin->start_used = 0;
		bin->end_used = 0;
	    }
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
		ins |= contig_insert_base2(io, crec, bin->child[i], pos, apos,
					   start_of_contig,
					   MIN(ch->pos, ch->pos + ch->size-1),
					   NMIN(ch->pos, ch->pos + ch->size-1),
					   base, conf, nbases, comp, hash);
		/* Children to the right of this one need pos updating too */
	    } else if (pos < MIN(ch->pos, ch->pos + ch->size-1)) {
		ch = get_bin(io, bin->child[i]);
		if (!(ch = cache_rw(io, ch))) {
		    cache_decr(io, bin);
		    return -1;
		}

		//printf("Mov bin %"PRIrec"\n", ch->rec);

		ch->pos+=nbases;
		if (ch->nseqs)
		    ins=1;
		ch->flags |= BIN_BIN_UPDATED;
	    }
	}
    }

    cache_decr(io, bin);

    return ins;
}

/*
 * Second pass on contig_insert_base(). 
 *
 * As tags can be in bins equal to, above or below the seqs we have to defer
 * tag shifting until we've fully walked the bin hierarchy, and then adjust
 * the tags in a second pass. The first pass filled out the hash table
 * indicating for each sequence at which point tags should be moved.
 */
static int contig_insert_tag2(GapIO *io, tg_rec crec, tg_rec bnum,
			      int pos, int apos, int start_of_contig,
			      int offset, int aoffset, int nbases,
			      int comp, HacheTable *hash) {
    int i;
    bin_index_t *bin;
    int f_a, f_b;

    bin = get_bin(io, bnum);
    cache_incr(io, bin);

    if (!(bin = cache_rw(io, bin)))
	return -1;

    /* Normalise pos for complemented bins */
    if (bin->flags & BIN_COMPLEMENTED) {
	comp ^= 1;
	pos = offset + bin->size - pos;
    } else {
	pos -= offset;
    }

    if (comp) {
	f_a = -1;
	f_b = aoffset + bin->size-1;
    } else {
	f_a = +1;
	f_b = aoffset;
    }

    //printf("Bin %"PRIrec" => pos=%d, offset=%d, comp=%d\n",
    //	   bin->rec, pos, offset, comp);

    /*
     * Loop again through contents if hash is non-null and move tags. Seq
     * tags move if their sequence moved. Consensus tags always move if
     * they're to the right.
     *
     * Assumption: a tag is in the same bin or a lower bin than the
     * sequence it belongs to. This is valid given tags are never
     * larger than the sequence they belong to, but does this hold true
     * after all possible types of editing?  Eg bin A children B,C.
     * Seq+tag put in bin A as they straddle B/C. Seq (& tag) then manually
     * reduced in size by deleting data. They're still in A of course.
     * Then we move seq A, which has the effect of taking it out of the
     * contig and putting it back in at a new location; now B or C.
     * We're safe provided the same logic is also applied to the tag, which
     * I believe it should be. Need to test.
     */
    //printf("Pass 2\n");
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);
	HacheItem *hi;

	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)
	    continue;

	hi = HacheTableSearch(hash, (char *)&r->pair_rec, sizeof(tg_rec));
	if (hi) {
	    //printf("Tag #%"PRIrec" in hash.  %d..%d (%d..%d) vs pos %d\n",
	    //	   r->rec, r->start, r->end,
	    //	   NMIN(r->start, r->end),
	    //	   NMAX(r->start, r->end),
	    //	   hi->data.i);

	    if (comp) {
		if (NMAX(r->start, r->end) <= (int64_t)hi->data.i) {
		    //puts("mov2");
		    r->start+=nbases;
		    r->end+=nbases;
		} else if (NMIN(r->start, r->end) <= (int64_t)hi->data.i) {
		    //puts("grow2");
		    r->end+=nbases;
		}
	    } else {
		if (NMIN(r->start, r->end) >= (int64_t)hi->data.i) {
		    //puts("mov1");
		    r->start+=nbases;
		    r->end+=nbases;
		} else if (NMAX(r->start, r->end) >= (int64_t)hi->data.i) {
		    //puts("grow1");
		    r->end+=nbases;
		}
	    }

	    bin->flags |= BIN_RANGE_UPDATED;
	}
    }

    /* Correct the bin dimensions */
    if (1) {
	if (bin->rng && ArrayMax(bin->rng)) {
	    int start = INT_MAX;
	    int end   = INT_MIN;
	    for (i = 0; i < ArrayMax(bin->rng); i++) {
		range_t *r = arrp(range_t, bin->rng, i);
	    
		if (r->flags & GRANGE_FLAG_UNUSED)
		    continue;

		if (start > r->start)
		    start = r->start;
		if (end   < r->end)
		    end   = r->end;
	    }

	    if (start != INT_MAX) {
		bin->start_used = start;
		bin->end_used = end;
	    } else {
		bin->start_used = 0;
		bin->end_used = 0;
	    }
	    bin->flags |= BIN_BIN_UPDATED;
	}
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
		contig_insert_tag2(io, crec, bin->child[i], pos, apos,
				   start_of_contig,
				   MIN(ch->pos, ch->pos + ch->size-1),
				   NMIN(ch->pos, ch->pos + ch->size-1),
				   nbases, comp, hash);
	    }
	}
    }

    cache_decr(io, bin);

    return 0;
}

int contig_insert_base_common(GapIO *io, contig_t **c,
			      int pos, char base, int conf, int nbases) {
    contig_t *n;
    int rpos, add_indel = 1;
    bin_index_t *bin;
    int bin_idx;
    tg_rec bin_rec;
    rangec_t rc;
    tg_rec ref_id = 0;
    int dir;
    HacheTable *hash = NULL;
    int ret;
    int cstart, cend;
    int nbases2 = nbases;

    if (pos < (*c)->start || pos > (*c)->end)
	return 0;
    
    consensus_valid_range(io, (*c)->rec, &cstart, &cend);

    if (pos-1 < cstart || pos > cend) {
	return 0;
    }
    
    if (!(n = cache_rw(io, *c)))
	return -1;
    *c = n;

    if (base == 0) {
	/*
	 * Shift rather than insert+shift. The difference between consensus
	 * inserts adding to all seqs (even in cutoff - perhaps considered
	 * a bug?) vs shift taking into account cutoff data means we need
	 * to know which seqs moved and which didn't so we can correctly
	 * move annotations.
	 */
    }
    hash = HacheTableCreate(4096, HASH_NONVOLATILE_KEYS | HASH_POOL_ITEMS
			    | HASH_ALLOW_DUP_KEYS | HASH_DYNAMIC_SIZE);

    ret = contig_insert_base2(io, n->rec, contig_get_bin(c), pos, pos,
			      pos == n->start,
			      contig_offset(io, c), contig_offset(io, c),
			      base, conf, nbases, 0, hash);
    ret |= contig_insert_tag2(io, n->rec, contig_get_bin(c), pos, pos,
			      pos == n->start,
			      contig_offset(io, c), contig_offset(io, c),
			      nbases, 0, hash);

    contig_visible_start(io, (*c)->rec, CITER_CSTART);
    contig_visible_end(io, (*c)->rec, CITER_CEND);

    if (1 != ret)
	return 0;

    rpos = padded_to_reference_pos(io, (*c)->rec, pos, &dir, NULL);
    rpos -= dir; /* Compensate for complemented contigs */

    /* Add a pad marker too for reference adjustments */
    if (find_refpos_marker(io, (*c)->rec, pos + 1, // or pos + nbases?
			   &bin_rec, &bin_idx, &rc) == 0) {
	assert((rc.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS);
	if ((rc.flags & GRANGE_FLAG_REFPOS_INDEL) == GRANGE_FLAG_REFPOS_DEL) {
	    /* Existing deletion, so decrement and maybe remove */
	    range_t *r;

	    bin = cache_search(io, GT_Bin, bin_rec);
	    bin = cache_rw(io, bin);
	    r = arrp(range_t, bin->rng, bin_idx);

	    if (rc.pair_rec <= nbases) {
		//printf("Remove DEL of size %d\n", rc.pair_rec);
		r->flags |= GRANGE_FLAG_UNUSED;
		r->rec = bin->rng_free;
		bin->rng_free = bin_idx;
		bin_incr_nrefpos(io, bin, -1);
		if (bin->start_used == r->start || bin->end_used == r->end)
		    bin_set_used_range(io, bin);
		nbases2 -= rc.pair_rec;
	    } else {
		//printf("Decrement DEL of size %d\n", rc.pair_rec);
		r->pair_rec -= nbases;
		nbases2 = 0;
	    }

	    bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;
	    add_indel =  nbases2 ? 1 : 0;
	}
	ref_id = rc.rec;
    }

    if (add_indel &&
	find_refpos_marker(io, (*c)->rec, pos - 1,
			   &bin_rec, &bin_idx, &rc) == 0) {
	assert((rc.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS);
	if ((rc.flags & GRANGE_FLAG_REFPOS_INDEL) == GRANGE_FLAG_REFPOS_DEL) {
	    /* Existing deletion, so decrement and maybe remove */
	    range_t *r;

	    bin = cache_search(io, GT_Bin, bin_rec);
	    bin = cache_rw(io, bin);
	    r = arrp(range_t, bin->rng, bin_idx);

	    if (rc.pair_rec <= nbases2) {
		//printf("Remove DEL of size %d\n", rc.pair_rec);
		r->flags |= GRANGE_FLAG_UNUSED;
		r->rec = bin->rng_free;
		bin->rng_free = bin_idx;
		bin_incr_nrefpos(io, bin, -1);
		if (bin->start_used == r->start || bin->end_used == r->end)
		    bin_set_used_range(io, bin);
		nbases2 -= rc.pair_rec;
	    } else {
		//printf("Decrement DEL of size %d\n", rc.pair_rec);
		r->pair_rec -= nbases2;
		nbases2 = 0;
	    }

	    bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;
	    add_indel =  nbases2 ? 1 : 0;
	}
	ref_id = rc.rec;
    }

    if (dir != -1 && add_indel) {
	range_t r;
	int i;

	for (i = 0; i < nbases2; i++) {
	    //printf("Adding insertion REFPOS at %d/%d\n", pos, rpos);
	    r.start    = pos;
	    r.end      = pos;
	    r.rec      = ref_id;
	    r.pair_rec = 0;
	    r.mqual    = rpos;

	    /* Dir needs checking */
	    r.flags    = GRANGE_FLAG_ISREFPOS
		| GRANGE_FLAG_REFPOS_INS
		| GRANGE_FLAG_REFPOS_FWD;
	    bin_add_range(io, c, &r, NULL, NULL, 0);
	}
    }

    /*
     * Inserting a base *may* grow the contig, but only if the
     * right-most sequence has been inserted into (the right-most base
     * may be a sequence that has cutoff data as 'pos').
     *
     * Similarly inserting a base *may* shrink the contig if the
     * left-most read is in left-cutoff data at pos.
     *
     * For safety, we simply have to recompute start/end.
     */

    /* Best guess */
    contig_set_end(io, c, contig_get_end(c)+nbases);

    /* Verify */
    consensus_unclipped_range(io, (*c)->rec, &cstart, &cend);
    if (contig_get_start(c) != cstart)
	contig_set_start(io, c, cstart);
    if (contig_get_end(c) != cend)
	contig_set_end(io, c, cend);

    (*c)->timestamp = io_timestamp_incr(io);

    if (hash)
	HacheTableDestroy(hash, 0);

    return ret;
}

int contig_insert_base(GapIO *io, contig_t **c, int pos, char base, int conf) {
    return contig_insert_base_common(io, c, pos, base, conf, 1) >= 0 ? 0 : -1;
}

int contig_insert_bases(GapIO *io, contig_t **c, int pos, char base, int conf,
			int nbases) {
    return contig_insert_base_common(io, c, pos, base, conf, nbases)
	>= 0 ? 0 : -1;
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

    if (parent->child[0] == bin->rec)
	parent->child[0] = 0;

    if (parent->child[1] == bin->rec)
	parent->child[1] = 0;

    parent->flags |= BIN_BIN_UPDATED;

    cache_decr(io, parent);

    return 0;
}

static int contig_delete_base2(GapIO *io, tg_rec crec, tg_rec bnum,
			       int pos, int apos, int start_of_contig,
			       int offset, int aoffset,
			       int base, int comp, HacheTable *hash,
			       int bcall) {
    int i, ins = 0;
    bin_index_t *bin;
    HacheData hd;
    int f_a, f_b;
    int min_r = INT_MAX, max_r = INT_MIN;

    bin = get_bin(io, bnum);
    cache_incr(io, bin);

    if (!(bin = cache_rw(io, bin)))
	return -1;

    /* Normalise pos for complemented bins */
    if (bin->flags & BIN_COMPLEMENTED) {
	comp ^= 1;
	pos = offset + bin->size - pos - 1;
    } else {
	pos -= offset;
    }

    if (comp) {
	f_a = -1;
	f_b = aoffset + bin->size-1;
    } else {
	f_a = +1;
	f_b = aoffset;
    }

    /* FIXME: add end_used (or start_used if complemented?) check here,
     * so we can shortcut looping through the bin if we're deleteing beyond
     * any of the bin contents.
     */

    /* Perform the delete into objects first */
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);

	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if (MAX(r->start, r->end) >= pos && MIN(r->start, r->end) <= pos) {
	    //printf("pos overlap obj %d/#%"PRIrec" %d in %d..%d\n",
	    //	   r->flags & GRANGE_FLAG_ISMASK, r->rec, pos,
	    //	   MIN(r->start, r->end),
	    //	   MAX(r->start, r->end));
	    /* Delete */
	    if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISREFPOS &&
		(r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO) {
		/* ISCONS? skip perhaps as we invalidate them anyway */
		seq_t *s = cache_search(io, GT_Seq, r->rec);
		int no_ins;

		if (ABS(r->end - r->start) + 1 != ABS(s->len)) {
		    verror(ERR_WARN, "contig_delete_base2", 
			   "Range start/end are inconsistent with seq len. ");
		}

		no_ins = 0;
		if (base) {
		    int r_start, r_end;
		    if (/*comp ^ */(s->len < 0)) {
			r_start = MIN(r->end - (s->right-1),
				      r->end - (s->left-1));
			r_end = MAX(r->end - (s->right-1),
				    r->end - (s->left-1));
		    } else {
			r_start = MIN(r->start + s->left-1,
				      r->start + s->right-1);
			r_end = MAX(r->start + s->right-1,
				    r->start + s->left-1);
		    }
		    if (pos < r_start || pos > r_end)
			no_ins = 1;

		    // printf("Rec %"PRIrec" visible %d..%d abs %d..%d\n",
		    //	   r->rec, r_start, r_end,
		    //	   NORM(r_start), NORM(r_end));

		    if (!no_ins) {
			if (r->start == r->end) {
			    /* Remove completely */
			    if (hash) {
				//printf("1 Mov %"PRIrec"\n", r->rec);
				hd.i = MAX(NMIN(r_start, r_end), apos);
				HacheTableAdd(hash, (char *)&r->rec,
					      sizeof(r->rec), hd, NULL);
			    }

			    r->flags |= GRANGE_FLAG_UNUSED;
			    r->rec = bin->rng_free;

			    bin->rng_free = i;
			    bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;
			    bin_incr_nseq(io, bin, -1);

			    if (bin->start_used == r->start ||
				bin->end_used == r->end)
				bin_set_used_range(io, bin);
			} else {
			    int bb;
			    //printf("DEL %"PRIrec" at %d\n", r->rec,
			    //       pos - MIN(r->start, r->end));

			    bb = ((r->flags & GRANGE_FLAG_ISMASK) ==
				  GRANGE_FLAG_ISCONS)
				? 0 : bcall;
			    sequence_delete_base2(io, &s,
						  pos - MIN(r->start, r->end),
						  0, bb);
			    if (hash) {
				//printf("1 Mov %"PRIrec"\n", r->rec);
				hd.i = MAX(NMIN(r_start, r_end), apos);
				HacheTableAdd(hash, (char *)&r->rec,
					      sizeof(r->rec), hd, NULL);
			    }
			    //printf("1 %"PRIrec"->end--\n", r->rec);
			    r->end--;
			    ins = 1;
			}
		    }
		}

		if (!base || no_ins) {
		    /* Shift instead, but only if in left-cutoff data */ 
		    if (/*comp ^ */ (s->len < 0)) {
			if (MIN(r->end - (s->right-1), r->end - (s->left-1))
			    >= pos) {
			    if (hash) {
				//printf("2 Mov %"PRIrec"\n", r->rec);
				hd.i = comp ? INT_MAX : INT_MIN;
				HacheTableAdd(hash, (char *)&r->rec,
					      sizeof(r->rec), hd, NULL);
			    }
			    //printf("2 %"PRIrec"->start/end--\n", r->rec);
			    r->start--;
			    r->end--;
			    ins = 1;
			}
		    } else {
			if (MIN(r->start + s->left-1, r->start + s->right-1)
			    >= pos) {
			    if (hash) {
				//printf("3 Mov %"PRIrec"\n", r->rec);
				hd.i = comp ? INT_MAX : INT_MIN;
				HacheTableAdd(hash, (char *)&r->rec,
					      sizeof(r->rec), hd, NULL);
			    }
			    //printf("3 %"PRIrec"->start/end--\n", r->rec);
			    r->start--;
			    r->end--;
			    ins = 1;
			}
		    }
		}
	    } else if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO) {
		if (base && !(r->flags & GRANGE_FLAG_TAG_SEQ)) {
		    /* Shrink consensus tags */
		    if (r->start == r->end) {
			r->flags |= GRANGE_FLAG_UNUSED;
			r->rec = bin->rng_free;

			bin->rng_free = i;
			bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;
			bin_incr_nanno(io, bin, -1);
			
			if (bin->start_used == r->start ||
			    bin->end_used == r->end)
			    bin_set_used_range(io, bin);
		    } else {
			r->end--;
		    }
		    ins = 1;
		}
	    }
	    
	    bin->flags |= BIN_RANGE_UPDATED;

	} else if (MIN(r->start, r->end) >= pos) {
	    //printf("pos to left of obj #%"PRIrec" %d %d..%d\n", 
	    //	   r->rec, pos,
	    //	   MIN(r->start, r->end),
	    //	   MAX(r->start, r->end));
	    if ( (r->flags&GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO) {
		/* Move */
		//printf("4 %"PRIrec"->start/end--\n", r->rec);
		r->start--;
		r->end--;
		ins = 1;
		bin->flags |= BIN_RANGE_UPDATED;

		if (hash) {
		    //printf("4 Mov %"PRIrec"\n", r->rec);
		    hd.i = comp ? INT_MAX : INT_MIN;
		    HacheTableAdd(hash, (char *)&r->rec,
				  sizeof(r->rec), hd, NULL);
		}
	    } else if (!(r->flags & GRANGE_FLAG_TAG_SEQ)) {
		/* Move consensus tags */
		r->start--;
		r->end--;
		ins = 1;
		bin->flags |= BIN_RANGE_UPDATED;
	    }
	} /* else pos to right of object */
    }

    /* Move the annotations too */
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);
	HacheItem *hi;

	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO) {
	    continue;
	}

	hi = HacheTableSearch(hash, (char *)&r->pair_rec, sizeof(tg_rec));
	if (!hi)
	    continue;

	//printf("Rec %"PRIrec" hd.i=%ld comp=%d NORM=%d,%d\n",
	//       r->rec, (long)hd.i, comp,
	//       NORM(r->start), NORM(r->end));

	if (r->start == r->end && NORM(r->start) == (int64_t)hi->data.i) {
	    //puts("del1/2");
	    r->flags |= GRANGE_FLAG_UNUSED;
	    r->rec = bin->rng_free;

	    bin->rng_free = i;
	    bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;
	    bin_incr_nanno(io, bin, -1);

	    if (bin->start_used == r->start || bin->end_used == r->end)
		bin_set_used_range(io, bin);
	} else {
	    if (comp) {
		if (NMAX(r->start, r->end) < (int64_t)hi->data.i) {
		    //puts("mov2");
		    r->start--;
		    r->end--;
		} else if (NMIN(r->start, r->end) <= (int64_t)hi->data.i) {
		    //puts("shrink2");
		    r->end--;
		}
	    } else {
		if (NMIN(r->start, r->end) > (int64_t)hi->data.i) {
		    //puts("mov1");
		    r->start--;
		    r->end--;
		} else if (NMAX(r->start, r->end) >= (int64_t)hi->data.i) {
		    //puts("shrink1");
		    r->end--;
		}
	    }
	}

	bin->flags |= BIN_RANGE_UPDATED;
    }

    /* Check new extents */
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);

	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if (min_r > r->start && r->start >= 0)
	    min_r = r->start;
	if (max_r < r->end)
	    max_r = r->end;
    }
    
    /* Adjust the bin dimensions */
    {
	if (bin->size != max_r+1) {
	    if (--bin->size <= 0) {
		fprintf(stderr, "Delete bin bin-%"PRIrec"\n", bin->rec);
		bin->size = 0;
		bin_delete(io, bin);
	    } else {
		ins = 1;
		if (bin->rng && ArrayMax(bin->rng)) {
		    if (min_r != INT_MAX) {
			bin->start_used = min_r;
			bin->end_used = max_r;
		    } else {
			bin->start_used = 0;
			bin->end_used = 0;
		    }
		}
		bin->flags |= BIN_BIN_UPDATED;
	    }
	} else {
	    /*
	     * Bin didn't shrink, which means we have something stretching
	     * to the right hand edge that wasn't deleted from nor was moved.
	     * (Ie it overlapped pos, but in right hand cutoff data).
	     *
	     * However not shrinking the bin can cause this bin to become
	     * larger (or at least beyond the edges of) its parent bin.
	     * Fix this now.
	     */
	    grow_bin(io, bin);
	}
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
		ins |= contig_delete_base2(io, crec, bin->child[i], pos, apos,
					   start_of_contig,
					   MIN(ch->pos, ch->pos + ch->size-1),
					   NMIN(ch->pos, ch->pos + ch->size-1),
					   base, comp, hash, bcall);
		/* Children to the right of this one need pos updating too */
	    } else if (pos < MIN(ch->pos, ch->pos + ch->size-1)) {
		ch = get_bin(io, bin->child[i]);
		if (!(ch = cache_rw(io, ch))) {
		    cache_decr(io, bin);
		    return -1;
		}

		//printf("Mov bin %"PRIrec"\n", ch->rec);

		ch->pos--;
		if (ch->nseqs)
		    ins=1;
		ch->flags |= BIN_BIN_UPDATED;
	    }
	}
    }

    cache_decr(io, bin);

    return ins;
}


static int contig_delete_base_fix(GapIO *io, tg_rec crec, tg_rec bnum,
				  int pos, int offset, int comp) {
    int i, r = 0;
    bin_index_t *bin;
    int f_a, f_b;
    contig_t *c = cache_search(io, GT_Contig, crec);

    bin = get_bin(io, bnum);
    cache_incr(io, bin);
    cache_incr(io, c);

    if (!(bin = cache_rw(io, bin)))
	return -1;

    /* Normalise pos for complemented bins */
    if (bin->flags & BIN_COMPLEMENTED)
	comp ^= 1;

    if (comp) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    /* Find objects now moved to -ve coords in bins, and re-home them */
    for (i = 0; bin->rng && i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i), r2, *r_out;

	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if (r->start < 0) {
	    bin_index_t *new_bin;

	    //printf("Bin %"PRIrec" obj %"PRIrec" with start %d. Loc %d..%d\n",
	    //	   bin->rec, r->rec, r->start,
	    //	   NORM(r->start), NORM(r->end));

	    /* Convert range coords to absolute locations */
	    r2 = *r;
	    r2.start = NORM(r2.start);
	    r2.end = NORM(r2.end);
	    if (r2.start > r2.end) {
		int t = r2.start;
		r2.start = r2.end;
		r2.end = t;
		r2.flags ^= GRANGE_FLAG_COMP1;
	    }

	    /* Remove from old bin */
	    bin_remove_item_from_bin(io, &c, &bin, 0, r->rec);

	    /* Add back in new location */
	    new_bin = bin_add_range(io, &c, &r2, &r_out, NULL, 0);
	    cache_incr(io, new_bin);

	    /* Fix any back references from the object */
	    switch (r2.flags & GRANGE_FLAG_ISMASK) {
	    case GRANGE_FLAG_ISANNO: {
		anno_ele_t *a = cache_search(io, GT_AnnoEle, r2.rec);
		a = cache_rw(io, a);
		a->bin = new_bin->rec;
		break;
	    }

	    case GRANGE_FLAG_ISCONS:
	    case GRANGE_FLAG_ISSEQ: {
		seq_t *s = cache_search(io, GT_Seq, r2.rec);
		int old_comp = bin_get_orient(io, s->bin);
		int new_comp = bin_get_orient(io, new_bin->rec);
		s = cache_rw(io, s);

		s->bin = new_bin->rec;
		s->bin_index = r_out - ArrayBase(range_t, new_bin->rng);

		//printf("Old bin comp=%d new bin comp=%d\n",
		//       old_comp, new_comp);

		if (new_comp != old_comp) {
		    s->len *= -1;
		    s->flags ^= SEQ_COMPLEMENTED;
		}

		/* Also move annotations */
		sequence_fix_anno_bins(io, &s);
		break;
	    }
	    }

	    cache_decr(io, new_bin);
	}
    }

    if (bin->child[0] || bin->child[1]) {
	bin_index_t *ch = NULL;

	/* Find the correct child node */
       	for (i = 0; i < 2; i++) {
	    if (!bin->child[i])
		continue;

	    ch = get_bin(io, bin->child[i]);

	    if (pos >= NMIN(ch->pos, ch->pos + ch->size-1) &&
		pos <= NMAX(ch->pos, ch->pos + ch->size-1)+1) {
		r |= contig_delete_base_fix(io, crec, bin->child[i], pos,
					    NMIN(ch->pos,
						 ch->pos + ch->size-1),
					    comp);
	    }
	}
    }

    cache_decr(io, bin);
    cache_decr(io, c);

    return r;
}

int contig_delete_base_common(GapIO *io, contig_t **c, int pos, int shift,
			      int base) {
    contig_t *n;
    int bin_idx;
    tg_rec bin_rec;
    rangec_t rc;
    int cur_del = 0;
    int done = 0, reduced;
    HacheTable *hash = NULL;
    int cstart, cend;

    consensus_valid_range(io, (*c)->rec, &cstart, &cend);

    if (pos-1 < cstart || pos > cend) {
	puts("Do nothing");
	return 0;
    }

    if (!(n = cache_rw(io, *c)))
	return -1;
    *c = n;

    /*
     * Look for a pad marker too. If we don't have one, create one to mark
     * a deletion from reference.
     * If we do have one neighbouring, either increment it (if deletion) or
     * remove it (if insertion).
     */
    if (find_refpos_marker(io, (*c)->rec, pos, &bin_rec, &bin_idx, &rc) == 0) {
	range_t *r;
	bin_index_t *bin;

	assert((rc.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS);

	bin = cache_search(io, GT_Bin, bin_rec);
	bin = cache_rw(io, bin);
	r = arrp(range_t, bin->rng, bin_idx);

	if ((rc.flags & GRANGE_FLAG_REFPOS_INDEL) == GRANGE_FLAG_REFPOS_INS) {
	    /* Existing insertion, so remove it */
	    //printf("Remove insertion marker\n");
	    r->flags |= GRANGE_FLAG_UNUSED;
	    r->rec = bin->rng_free;
	    bin->rng_free = bin_idx;
	    bin_incr_nrefpos(io, bin, -1);

	    if (bin->start_used == r->start || bin->end_used == r->end)
		bin_set_used_range(io, bin);
	    done = 1;
	} else {
	    /* Existing deletion, add it to neighbouring position */
	    //printf("DEL %d marker to apply to right\n", r->pair_rec);
	    cur_del = r->pair_rec;
	}

	bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;
    }

    if (!done && find_refpos_marker(io, (*c)->rec, pos+1, &bin_rec,
				    &bin_idx, &rc) == 0) {
	range_t *r;
	bin_index_t *bin;
	int ins_type;

	assert((rc.flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS);

	bin = cache_search(io, GT_Bin, bin_rec);
	bin = cache_rw(io, bin);
	r = arrp(range_t, bin->rng, bin_idx);

	ins_type = (rc.flags & GRANGE_FLAG_REFPOS_INDEL);
	cur_del++;

	if (ins_type == GRANGE_FLAG_REFPOS_INS && cur_del == 1) {
	    //printf("Remove INS marker\n");
	    /* 1 del + 1 ins => remove marker */
	    r->flags |= GRANGE_FLAG_UNUSED;
	    r->rec = bin->rng_free;
	    bin->rng_free = bin_idx;
	    bin_incr_nrefpos(io, bin, -1);

	    if (bin->start_used == r->start || bin->end_used == r->end)
		bin_set_used_range(io, bin);
	} else {
	    /* Otherwise N del + 1 ins => N-1 del */
	    if (ins_type == GRANGE_FLAG_REFPOS_INS) {
		cur_del--;
		r->flags &= ~GRANGE_FLAG_REFPOS_INDEL;
		r->flags |= GRANGE_FLAG_REFPOS_DEL;
		r->pair_rec = cur_del;
		//printf("Replace INS with DEL %d\n", r->pair_rec);
	    } else {
		/* N del + M ndel => N+M del */
		//printf("Inc DEL %d to DEL %d\n",
		//       r->pair_rec, r->pair_rec + cur_del);
		r->pair_rec += cur_del;
	    }
	}
	
	bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;

    } else if (!done) {
	/* No existing marker at pos or pos+1 => a new +1 DEL marker */
	range_t r;
	int dir;
	int rpos = padded_to_reference_pos(io, (*c)->rec, pos+1, &dir, NULL);

	if (dir != -1) {
	    rpos += dir; /* Compensate for complemented contigs */

	    /* Create a deletion marker */
	    //printf("Adding deletion REFPOS at %d/%d\n", pos, rpos);
	    r.start    = pos+1;
	    r.end      = pos+1;
	    r.rec      = 0;
	    r.pair_rec = 1+cur_del;
	    r.mqual    = rpos;

	    /* Dir needs checking */
	    r.flags    = GRANGE_FLAG_ISREFPOS
		| GRANGE_FLAG_REFPOS_DEL
		| GRANGE_FLAG_REFPOS_FWD;
	    bin_add_range(io, c, &r, NULL, NULL, 0);
	} /* else no refpos data => don't keep tracking it */
    }

    /*
     * When shifting only we need to track which sequences moved so we can
     * move annotations too. See comments in insertion code for why.
     */
    hash = HacheTableCreate(4096, /*HASH_NONVOLATILE_KEYS | */HASH_POOL_ITEMS
			    | HASH_ALLOW_DUP_KEYS | HASH_DYNAMIC_SIZE);

    reduced = contig_delete_base2(io, n->rec, contig_get_bin(c), pos, pos,
				  pos == n->start,
				  contig_offset(io, c), contig_offset(io, c),
				  !shift, 0, hash, base);

    /*
     * Deletion can move objects left if the deleted coord is in the left
     * hand cutoff of a sequence. This in turn can push a sequence to
     * negative coords within the bin.
     *
     * We can't fix this in real time very easily as we may move the sequence
     * to a bin we haven't yet visited, causing it to then be found later on
     * in the recursive algorithm and moved yet again. Instead we have a
     * second pass fixing up these issues when they occur.
     */
    contig_delete_base_fix(io, n->rec, contig_get_bin(c),
    			   pos, contig_offset(io, c), 0);

    /*
     * Deleting a base can change contig dimensions in unexpected ways, as
     * sometimes it moves sequences and other times it deletes from sequences,
     * depending on whether we're in cutoff data or not.
     *
     * The easiest solution is just recalculate both ends.
     */
    contig_visible_start(io, (*c)->rec, CITER_CSTART);
    contig_visible_end(io, (*c)->rec, CITER_CEND);
    consensus_unclipped_range(io, (*c)->rec, &cstart, &cend);
    if (contig_get_start(c) != cstart)
	contig_set_start(io, c, cstart);
    if (contig_get_end(c) != cend)
	contig_set_end(io, c, cend);

    (*c)->timestamp = io_timestamp_incr(io);

    if (hash)
	HacheTableDestroy(hash, 0);

    return reduced;
}

int contig_delete_base(GapIO *io, contig_t **c, int pos) {
    return contig_delete_base_common(io, c, pos, 0, 0) >= 0 ? 0 : -1;
}

int contig_delete_pad(GapIO *io, contig_t **c, int pos) {
    return contig_delete_base_common(io, c, pos, 0, '*') >= 0 ? 0 : -1;
}

/*
 * Like an insertion or deletion, but only perform the reading shift stage
 * and not the actual sequence insertion/deletion.
 *
 * The code shares a lot with ins/del, so we just call those functions with
 * specific arguments. Eg insertion base==0 implies shift only.
 *
 * dir +1 => shift right
 *     -1 => shift left
 *
 * Currently 1bp only, but in future we may support +n/-n
 */
int contig_shift_base(GapIO *io, contig_t **c, int pos, int dir) {
    if (dir > 0)
	return contig_insert_base_common(io, c, pos, 0, 0, 1);
    else
	return contig_delete_base_common(io, c, pos+1, 1, 0);
}

contig_t *contig_new(GapIO *io, char *name) {
    tg_rec rec;
    contig_t *c, init_c;

    memset(&init_c, 0, sizeof(contig_t));
    init_c.name = name;

    /* Allocate our contig */
    //rec = io->iface->contig.create(io->dbh, NULL);
    rec = cache_item_create(io, GT_Contig, &init_c);

    /* Initialise it */
    c = (contig_t *)cache_search(io, GT_Contig, rec);
    c = cache_rw(io, c);
    c->bin = bin_new(io, 0, io->min_bin_size, rec, GT_Contig);
    if (name)
        contig_set_name(io, &c, name);
    else
        c->name = NULL;

    /* Add it to the contig order too */
    io->contig_order = cache_rw(io, io->contig_order);
    io->db = cache_rw(io, io->db);
    ARR(tg_rec, io->contig_order, io->db->Ncontigs++) = rec;

    c->timestamp = io_timestamp_incr(io);

    /* Add to the new contigs list */
    if (name)
	add_to_list("new_contigs", name);

    return c;
}

/*
 * Looks for a contig by name and if found it returns that contig object.
 *
 * Returns contig_t ptr on success
 *         NULL on failure (not found)
 */
contig_t *find_contig_by_name(GapIO *io, char *name) {
    tg_rec rec = io->iface->contig.index_query(io->dbh, name, 0);
    return rec > 0 ? (contig_t *)cache_search(io, GT_Contig, rec) : NULL;
}

tg_rec contig_index_query(GapIO *io, char *name) {
    return io->iface->contig.index_query(io->dbh, name, 0);
}

tg_rec contig_index_query_prefix(GapIO *io, char *prefix) {
    return io->iface->contig.index_query(io->dbh, prefix, 0);
}

int contig_index_update(GapIO *io, char *name, int name_len, tg_rec rec) {
    char n2[1024];
    tg_rec r;
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
static void pair_rangec(GapIO *io, tg_rec crec, rangec_t *r, int count) {
    int chunk = 5000; /* number of reads to deal with, 5k works well, not sure why */
    int count_start = 0;
    int count_end;
    HacheTable *unpaired;
    HacheIter *it;
    HacheItem *rec;
    int unpaired_count = 0;
    
    unpaired = HacheTableCreate(count, HASH_NONVOLATILE_KEYS | HASH_POOL_ITEMS);
    unpaired->name = "pair_rangec unpaired";   
    
    if (count < chunk) {
    	count_end = count;
    } else {
    	count_end = chunk;
    }
    
    while (count_start < count) {
	int i;
	HacheTable *h;
	HacheIter *iter;
	HacheItem *hi;

	/* Build a hash on record number */
	h = HacheTableCreate((count_end - count_start), HASH_NONVOLATILE_KEYS | HASH_POOL_ITEMS);
	h->name = "pair_rangec()";

	for (i = count_start; i < count_end; i++) {
	    HacheData hd;
	    hd.i = i;
	    HacheTableAdd(h, (char *)&r[i].rec, sizeof(r[i].rec), hd, NULL);
	}

	/* Iterate through hash linking the pairs together */
	iter = HacheTableIterCreate();

	while (NULL != (hi = HacheTableIterNext(h, iter))) {
	    HacheItem *pair;
	    int p;
	    i = hi->data.i;
	    assert(i < count && i >= 0);

	    pair = HacheTableSearch(h, (char *)&r[i].pair_rec, sizeof(r[i].rec));

	    if (pair) {
		p = pair->data.i;
		assert(p < count && p >= 0);

		if (r[p].rec > r[i].rec)
		    continue;

		r[i].pair_ind    = p;
		r[p].pair_ind    = i;
		r[i].pair_start  = r[p].start;
		r[i].pair_end    = r[p].end;
		r[i].pair_mqual  = r[p].mqual;
		r[i].pair_contig = crec;
		r[p].pair_start  = r[i].start;
		r[p].pair_end    = r[i].end;
		r[p].pair_mqual  = r[i].mqual;
		r[p].pair_contig = crec;

		r[i].pair_timestamp = gio_base(io)->db->timestamp;
		r[p].pair_timestamp = gio_base(io)->db->timestamp;

		/*
		 * Convert relative (to bin) orientation to absolute (to
		 * contig) orientation via r->comp field.
		 *
		 * COMP1 => this read
		 * COMP2 => paired read
		 */
		if (((r[i].flags &  GRANGE_FLAG_COMP1) != 0) ^ r[i].comp) {
		    r[i].flags |=  GRANGE_FLAG_COMP1;
		    r[p].flags |=  GRANGE_FLAG_COMP2;
		} else {
		    r[i].flags &= ~GRANGE_FLAG_COMP1;
		    r[p].flags &= ~GRANGE_FLAG_COMP2;
		}
	    
		if (((r[p].flags &  GRANGE_FLAG_COMP1) != 0) ^ r[p].comp) {
		    r[i].flags |=  GRANGE_FLAG_COMP2;
		    r[p].flags |=  GRANGE_FLAG_COMP1;
		} else {
		    r[i].flags &= ~GRANGE_FLAG_COMP2;
		    r[p].flags &= ~GRANGE_FLAG_COMP1;
		}

		r[i].flags |= GRANGE_FLAG_CONTIG;
		r[p].flags |= GRANGE_FLAG_CONTIG;

		r[i].flags &= ~GRANGE_FLAG_PEND_MASK;
		if ((r[p].flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD)
		    r[i].flags |= GRANGE_FLAG_PEND_FWD;
		else
		    r[i].flags |= GRANGE_FLAG_PEND_REV;

		r[p].flags &= ~GRANGE_FLAG_PEND_MASK;
		if ((r[i].flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD)
		    r[p].flags |= GRANGE_FLAG_PEND_FWD;
		else
		    r[p].flags |= GRANGE_FLAG_PEND_REV;
	    } else {
	    	HacheData hd;
		hd.i = i;
		HacheTableAdd(unpaired,(char *)&r[i].rec, sizeof(r[i].rec), hd, NULL); 
		unpaired_count++;
	    }
	}

	/* Tidy up */
	HacheTableIterDestroy(iter);
	HacheTableDestroy(h, 0);

	count_start = count_end;

	if ((count_end + chunk) > count) {
	    count_end = count;
	} else {
	    count_end += chunk;
	}
    }
    

    /* do the remaining unpaired reads */
    it = HacheTableIterCreate();

    while (NULL != (rec = HacheTableIterNext(unpaired, it))) {
	HacheItem *pair;
	int i = rec->data.i;
	int p;
	assert(i < count && i >= 0);

	pair = HacheTableSearch(unpaired, (char *)&r[i].pair_rec, sizeof(r[i].rec));

	if (pair) {
	    p = pair->data.i;
	    // assert(p < count && p >= 0);

	    if (r[p].rec > r[i].rec)
		continue;

	    r[i].pair_ind    = p;
	    r[p].pair_ind    = i;
	    r[i].pair_start  = r[p].start;
	    r[i].pair_end    = r[p].end;
	    r[i].pair_mqual  = r[p].mqual;
	    r[i].pair_contig = crec;
	    r[p].pair_start  = r[i].start;
	    r[p].pair_end    = r[i].end;
	    r[p].pair_mqual  = r[i].mqual;
	    r[p].pair_contig = crec;

	    if (((r[i].flags & GRANGE_FLAG_COMP1) != 0) ^ r[i].comp) {
		r[i].flags |=  GRANGE_FLAG_COMP1;
		r[p].flags |=  GRANGE_FLAG_COMP2;
	    } else {
		r[i].flags &= ~GRANGE_FLAG_COMP1;
		r[p].flags &= ~GRANGE_FLAG_COMP2;
	    }
	    
	    if (((r[p].flags & GRANGE_FLAG_COMP1) != 0) ^ r[p].comp) {
		r[i].flags |=  GRANGE_FLAG_COMP2;
		r[p].flags |=  GRANGE_FLAG_COMP1;
	    } else {
		r[i].flags &= ~GRANGE_FLAG_COMP2;
		r[p].flags &= ~GRANGE_FLAG_COMP1;
	    }

	    r[i].flags |= GRANGE_FLAG_CONTIG;
	    r[p].flags |= GRANGE_FLAG_CONTIG;


	    r[i].flags &= ~GRANGE_FLAG_PEND_MASK;
	    if ((r[p].flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD)
		r[i].flags |= GRANGE_FLAG_PEND_FWD;
	    else
		r[i].flags |= GRANGE_FLAG_PEND_REV;

	    r[p].flags &= ~GRANGE_FLAG_PEND_MASK;
	    if ((r[i].flags & GRANGE_FLAG_END_MASK) == GRANGE_FLAG_END_FWD)
		r[p].flags |= GRANGE_FLAG_PEND_FWD;
	    else
		r[p].flags |= GRANGE_FLAG_PEND_REV;
	} else {
	    r[i].pair_ind = -1;
	}
    }

    /* Tidy up */
    HacheTableIterDestroy(it);
    HacheTableDestroy(unpaired, 0);
}


/*
 * Sort comparison function for range_t; sort by ascending position of
 * object start position.
 */
static int sort_range_by_x(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    int d;

    /* By X primarily */
    if ((d = r1->start - r2->start))
	return d;
    
#if 0
    /* And then by unmapped first, mapped second */
    {
	int m1, m2;
	m1 = (r1->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISUMSEQ;
	m2 = (r1->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISUMSEQ;
	if (m1 != m2)
	    return m2 - m1;
    }
#endif

    /* And finally by recno, allowing for 64-bit quantities. */
    if (r1->rec > r2->rec)
	return 1;
    else if (r1->rec < r2->rec)
	return -1;
    else
	return 0;
}

/*
 * Sort comparison function for range_t; sort by ascending position of
 * object end position.
 */
static int sort_range_by_x_end(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    int d;

    /* By X primarily */
    if ((d = r1->end - r2->end))
	return d;
    
    /* And finally by recno, allowing for 64-bit quantities. */
    if (r1->rec > r2->rec)
	return 1;
    else if (r1->rec < r2->rec)
	return -1;
    else
	return 0;
}

static GapIO *sort_io = NULL; /* Hack to pass into qsort() */

/* Sort comparison function for range_t; sort by ascending position */
static int sort_range_by_y(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    if (r1->y != r2->y)
        return r1->y - r2->y;
    if (r1->start != r2->start)
    	return r1->start - r2->start;
    return (r1->flags & GRANGE_FLAG_ISMASK) - (r2->flags & GRANGE_FLAG_ISMASK);
}

/**********************************************************************
* TRY A NEW SORTING SYSTEM
***********************************************************************/

/* Variables for holding sort related values,
   not the prettiest of solutions having it here.
*/

static int (*primary)(const void *, const void *);
static int (*secondary)(const void *, const void *);
static int base_sort_pos = 1;


// empty sort
static int no_sort(const void *v1, const void *v2) {
    return 0;
}


// the sort types
static int simple_sort_range_by_x(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;

    return r1->start - r2->start;
}


static int simple_sort_range_by_x_end(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;

    return r1->end - r2->end;
}


static int simple_sort_range_by_x_clipped(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    int r1_start, r2_start;
    
    /* Use clipped coordinate in seqs */
    if ((r1->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	seq_t *s = cache_search(sort_io, GT_Seq, r1->rec);
	if ((s->len < 0) ^ r1->comp) {
	    r1_start = r1->start + ABS(s->len) - (s->right-1) - 1;
	} else {
	    r1_start = r1->start + s->left-1;
	}
    } else {
	r1_start = r1->start;
    }

    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	seq_t *s = cache_search(sort_io, GT_Seq, r2->rec);
	if ((s->len < 0) ^ r2->comp)
	    r2_start = r2->start + ABS(s->len) - (s->right-1) - 1;
	else
	    r2_start = r2->start + s->left-1;
    } else {
	r2_start = r2->start;
    }

    return r1_start - r2_start;
}


static int simple_sort_range_by_x_clipped_end(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    int r1_end, r2_end;

    /* Use clipped coordinate in seqs */
    if ((r1->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	seq_t *s = cache_search(sort_io, GT_Seq, r1->rec);
	if ((s->len < 0) ^ r1->comp)
	    r1_end = r1->start + ABS(s->len) - (s->left-1) - 1;
	else
	    r1_end = r1->start + s->right-1;
    } else {
	r1_end = r1->end;
    }

    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
	seq_t *s = cache_search(sort_io, GT_Seq, r2->rec);
	if ((s->len < 0) ^ r2->comp)
	    r2_end = r2->start + ABS(s->len) - (s->left-1) - 1;
	else
	    r2_end = r2->start + s->right-1;
    } else {
	r2_end = r2->end;
    }

    return r1_end - r2_end;
}


static int simple_sort_range_by_template(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    
    char template1[2048];
    char template2[2048];
    
    memset(template1, '\0', 2048);
    memset(template2, '\0', 2048);
    
    /* use template name if exists, otherwise by name */
    if ((r1->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
    	seq_t *s1 = cache_search(sort_io, GT_Seq, r1->rec);
	
	if (s1->template_name_len > 0) {
	    strncpy(template1, s1->name, s1->template_name_len);
	    template1[s1->template_name_len] = '\0';
	} else if (s1->name_len > 0) {
	    strcpy(template1, s1->name);
	} else {
	    return 0;
	}
    }

    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
    	seq_t *s2 = cache_search(sort_io, GT_Seq, r2->rec);
	
	if (s2->template_name_len > 0) {
	    strncpy(template2, s2->name, s2->template_name_len);
	    template2[s2->template_name_len] = '\0';
	} else if (s2->name_len > 0) {
    	    strcpy(template2, s2->name);	    
	} else {
	    return 0;
	}
    }
    
    return strcmp(template1, template2);
}


/* As we already have rangec_t calculated we can use this instead of the
 * slower sequence_get_base
 *
 * Returns 0 on success and fills *base, *cutoff
 *        -1 on failure
 */
static int get_base(seq_t *s, const rangec_t *r, int pos, char *base, int *cutoff) {
    
    if (pos < 0 || pos >= ABS(s->len)) return -1;

    if ((s->len < 0) ^ r->comp) {             // is complemented
    	pos = ABS(s->len) - 1 - pos;
	*base = complement_base(s->seq[pos]);
    } else {
    	*base = s->seq[pos];
    }
    
    if (pos < s->left - 1 || pos >= s->right) {
    	*cutoff = 1;
    } else {
    	*cutoff = 0;
    }
    
    return 0;
}

/*
 * Collate bases by alphabetical (A, C, G, T), N(123, first value are z),
 * "*"(124), and then any arbitrary unknown/failure(125).
 */
static int simple_sort_by_base(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    char b1, b2;
    
    if ((r1->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
    	seq_t *s = cache_search(sort_io, GT_Seq, r1->rec);
	int cutoff = 0;
	
	if (get_base(s, r1, base_sort_pos - r1->start, &b1, &cutoff) || cutoff) {
	    b1 = 125;
	}

    	if (b1 == 'N') {
	    b1 = 123;
	} else if (b1 == '*'){
	    b1 = 124;
	}
    } else {
	b1 = 125; /* Unknowns are sorted to end */
    }
    
    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
    	seq_t *s = cache_search(sort_io, GT_Seq, r2->rec);
	int cutoff = 0;

	if (get_base(s, r2, base_sort_pos - r2->start, &b2, &cutoff) || cutoff) {
	    b2 = 125;
	}
	
    	if (b2 == 'N') {
	    b2 = 123;
	} else if (b2 == '*') {
	    b2 = 124;
	}
    } else {
	b2 = 125;
    }
    
    return b1 - b2;
}

static int simple_sort_by_strand(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    
    int strand1 = 0;
    int strand2 = 0;
    
    if ((r1->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
    	strand1 = (r1->flags & GRANGE_FLAG_COMP1) ^ r1->comp;
    }
    
    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
    	strand2 = (r2->flags & GRANGE_FLAG_COMP1) ^ r2->comp;
    }
    
    return strand1 - strand2;
}


static int simple_sort_by_template_status(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    int status1 = 0, status2 = 0;
    seq_t *s1, *s2;
    
    if ((r1->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ)
	if ((s1 = cache_search(sort_io, GT_Seq, r1->rec)))
	    status1 = sequence_get_template_info(sort_io, s1, NULL, NULL);

    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ)
	if ((s2 = cache_search(sort_io, GT_Seq, r2->rec)))
	    status2 = sequence_get_template_info(sort_io, s2, NULL, NULL);
    
    return status1 - status2;
}

/*
 * Only on DB format >= 5
 */
static int simple_sort_by_library(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    
    int library1 = 0;
    int library2 = 0;
    
    if ((r1->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
    	library1 = r1->library_rec;
    }
    
    if ((r2->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
    	library2 = r2->library_rec;
    }
    
    return library1 - library2;
}


static int simple_sort_range_by_tech_x(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;

    return r1->seq_tech - r2->seq_tech;
}

static int simple_sort_by_sequence(const void *v1, const void *v2) {
    const rangec_t *r1 = (const rangec_t *)v1;
    const rangec_t *r2 = (const rangec_t *)v2;
    
    if (r2->seq_match != r1->seq_match) {
    	return r2->seq_match - r1->seq_match;
    } else {
    	return r2->seq_hash - r1->seq_hash;
    }
}


/* runs the simple sort choices */
static int chosen_sort(const void *v1, const void *v2) {
    
    int ret;
    
    if ((ret = (*primary)(v1, v2))) {
    	return ret;
    } else if ((ret = (*secondary)(v1, v2))) {
    	return ret;
    } else {
    	return sort_range_by_x(v1, v2);
    }
}


void contig_set_base_sort_point(int pos) {
    base_sort_pos = pos;
}


void contig_set_default_sort(seq_sort_t *set, int pri, int sec) {

    switch (pri) {
    case 1:
	set->p_sort = CSIR_SORT_BY_X | CSIR_SORT_BY_SEQ_TECH;
	break;
	
    case 2:
	set->p_sort = CSIR_SORT_BY_X | CSIR_SORT_BY_CLIPPED;
	break;
	
    case 3:
	set->p_sort = CSIR_SORT_BY_X;
	break;
	
    case 4:
	set->p_sort = CSIR_SORT_BY_TEMPLATE;
	break;
	
    case 5:
	set->p_sort = CSIR_SORT_BY_STRAND;
	break;

    case 6:
	set->p_sort = CSIR_SORT_BY_BASE;
	break;
	
    case 7:
	set->p_sort = CSIR_SORT_BY_SEQUENCE;
	break;

    case 8:
	set->p_sort = CSIR_SORT_BY_TEMPLATE_STATUS;
	break;

    case 9:
	set->p_sort = CSIR_SORT_BY_LIBRARY;
	break;

    default:
	set->p_sort = CSIR_SORT_BY_Y;
    }
    
    switch (sec) {
    case 1:
	set->s_sort = CSIR_SORT_BY_X | CSIR_SORT_BY_SEQ_TECH;
	break;
	
    case 2:
	set->s_sort = CSIR_SORT_BY_X | CSIR_SORT_BY_CLIPPED;
	break;
	
    case 3:
	set->s_sort = CSIR_SORT_BY_X;
	break;
	
    case 4:
	set->s_sort = CSIR_SORT_BY_TEMPLATE;
	break;
	
    case 5:
	set->s_sort = CSIR_SORT_BY_STRAND;
	break;

    case 6:
	set->s_sort = CSIR_SORT_BY_BASE;
	break;

    case 7:
	set->s_sort = CSIR_SORT_BY_SEQUENCE;
	break;

    case 8:
	set->s_sort = CSIR_SORT_BY_TEMPLATE_STATUS;
	break;

    case 9:
	set->s_sort = CSIR_SORT_BY_LIBRARY;
	break;

    default:
	set->s_sort = CSIR_SORT_BY_Y;
    }
}


/* return a pointer to a function */
static int (*set_sort(int job))(const void *, const void *) {

    if (job & CSIR_SORT_BY_TEMPLATE) {
	return simple_sort_range_by_template;
    } else if (job & CSIR_SORT_BY_SEQUENCE) {
    	return simple_sort_by_sequence;
    } else if (job & CSIR_SORT_BY_STRAND) {
    	return simple_sort_by_strand;
    } else if (job & CSIR_SORT_BY_BASE) {
    	return simple_sort_by_base;
    } else if (job & CSIR_SORT_BY_TEMPLATE_STATUS) {
	return simple_sort_by_template_status;
    } else if (job & CSIR_SORT_BY_LIBRARY) {
	return simple_sort_by_library;
    } else if (job & CSIR_SORT_BY_XEND) {
	if (job & CSIR_SORT_BY_CLIPPED) {
	    return simple_sort_range_by_x_clipped_end;
	} else {
	    return simple_sort_range_by_x_end;
	}
    } else if (job & (CSIR_SORT_BY_X | CSIR_SORT_BY_Y)) {
    	if (job & CSIR_SORT_BY_SEQ_TECH) {
	    return simple_sort_range_by_tech_x;
	} else if (job & CSIR_SORT_BY_CLIPPED) {
	    return simple_sort_range_by_x_clipped;
	} else if (job & CSIR_SORT_BY_X || !(job & CSIR_DEFAULT)) {
	    return simple_sort_range_by_x;
	} 
    }
    
    return no_sort;
}


static int find_string_match_values(GapIO *io, rangec_t *r, seq_sort_t *setting, int count) {
    char *seq = NULL;
    int seq_size = (setting->end - setting->start) + 1;
    int i, j, k;
    int c_start = setting->start;
    int c_end   = setting->end;
    
    if (NULL == (seq = (char *)malloc(seq_size * sizeof(char)))) {
    	fprintf(stderr, "Error: memory allocation failure in find_string_match_values\n");
	return 1;
    }
    
    seq[seq_size - 1] = 0;
    
    // get the selection, either from the consensus or a read
    
    if (setting->rec == 0) return 1;

    if (setting->type == GT_Contig) {
    	if (-1 == (calculate_consensus_simple(io, setting->rec,
	    	     setting->start, setting->end, seq, NULL))) {
	    fprintf(stderr, "Error: unable to retrieve consensus data\n");
	    free(seq);
	    return 1;
	}
    } else {
    	seq_t *s = cache_search(io, GT_Seq, setting->rec);
	
	// look for the right record in range for the pos info
	// there is probably a better way but do this for now
	
	for (i = 0; i < count; i++) {
	    if (r[i].rec == setting->rec) break;
	}
	
	if (i == count) { // rec not there, probably out of range
	    free(seq);
	    return 1;
	}
	
	
	k = 0;
	
	// map the select start to the consensus start
	c_start += r[i].start;
	c_end   += r[i].start;;
	
    	for (j = setting->start; j <=  setting->end; j++) {
	    char base;
	    int cutoff;
	    
	    if (get_base(s, &r[i], j, &base, &cutoff) || cutoff) {
	    	seq[k] = '_';
	    } else {
	    	seq[k] = base;
	    }
	    
	    k++;
	}
    }
    
    // go through the reads and assign matching scores
    
    for (i = 0; i < count; i++) {
    	rangec_t *rn = &r[i];
	
    	if ((rn->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
    	    seq_t *s = cache_search(io, GT_Seq, rn->rec);

	    k = 0;
	    rn->seq_match = 0;
	    
	    // create a hash to group like sequences together
	    rn->seq_hash = 1315423911;

	    for (j = (c_start - rn->start); j <= (c_end - rn->start); j++) {
		char base;
		int cutoff;

		if (get_base(s, &r[i], j, &base, &cutoff) || cutoff) {
		    base = '_';
		}
		
		// naive match scoring, +3 match, -1 mismatch, 0 for off the read end
		
		if (seq[k] == base) {
	    	    rn->seq_match += 3;
		} else if (base != '_') {
	    	    rn->seq_match --;
		}

		k++;
		
		rn->seq_hash ^= ((rn->seq_hash << 5) + base + (rn->seq_hash >> 2));
	    }
	    
	    rn->seq_hash &= 0x0FFF;

	} else {
	   rn->seq_match = 0;
	   rn->seq_hash  = 0; 
	} 
    }
    
    free(seq);

    return 0;
}




/* Pick Y coordinates for ranges */
int x_cmp(struct xy_pair *y1, struct xy_pair *y2) {
    int d = (y1->x) - (y2->x);
    return d ? d : y1->y - y2->y;
}

int y_cmp(struct xy_pair *y1, struct xy_pair *y2) {
    return y1->y - y2->y;
}

SPLAY_HEAD(xTREE, xy_pair);
SPLAY_HEAD(yTREE, xy_pair);
SPLAY_PROTOTYPE(xTREE, xy_pair, x_link, x_cmp);
SPLAY_PROTOTYPE(yTREE, xy_pair, y_link, y_cmp);
SPLAY_GENERATE(xTREE, xy_pair, x_link, x_cmp);
SPLAY_GENERATE(yTREE, xy_pair, y_link, y_cmp);

#define KEEP_Y

#if 0
static void xtree_dump(struct xTREE *tree, xy_pair *node) {
    static int level = 0;
    char *space = "                                                  ";

    if (!node) {
	printf("\nxTREE: %p\n", tree);
	node = SPLAY_ROOT(tree);
	printf("Root: %p\n", node);
	level = 0;

	if (!node)
	    return;
    }
    
    printf("%.*sleft=%p\n", level, space, SPLAY_LEFT(node, x_link));
    level++;
    if (SPLAY_LEFT(node, x_link))
	xtree_dump(tree, SPLAY_LEFT(node, x_link));
    level--;

    printf("%.*sright=%p\n", level, space, SPLAY_RIGHT(node, x_link));
    level++;
    if (SPLAY_RIGHT(node, x_link))
	xtree_dump(tree, SPLAY_RIGHT(node, x_link));
    level--;
}

static void ytree_dump(struct yTREE *tree, xy_pair *node) {
    xy_pair *n;
    static int level = 0;
    char *space = "                                                  ";

    if (!node) {
	printf("\nyTREE: %p\n", tree);
	node = SPLAY_ROOT(tree);
	printf("Root: %p\n", node);
	level = 0;

	if (!node)
	    return;
    }
    
    printf("%.*sleft=%p\n", level, space, n=SPLAY_LEFT(node, y_link));
    level++;
    if (n) ytree_dump(tree, n);
    level--;

    printf("%.*sright=%p\n", level, space, n=SPLAY_RIGHT(node, y_link));
    level++;
    if (n) ytree_dump(tree, n);
    level--;
}
#endif

/*
 * This algorithm allocates Y coordinates to the horizontal lines listed
 * in tl[0..ntl-1].
 *
 * To keep track of this we expect sorted data, by X. We then keep track
 * of every display line as either a line currently in use (with each line
 * indexed in two trees - sorted on X and sorted on Y) or a line no longer
 * in use (indexed in a 3rd tree sorted by Y).
 * (NB, change macros SPLAY_ to RB_ for a red-black tree
 * implementation instead.)
 *
 * axtree - active rows sorted by X
 * aytree - active rows sorted by Y
 * iytree - inactive rows sorted by Y
 */
static int compute_ypos(rangec_t *r, int nr, int job) {
    int i, j;
    struct xy_pair *node, *curr, *next;
    int yn = -1;
    
    /* Simple case */
    if (job & CSIR_ALLOCATE_Y_SINGLE) {
	for (i = j = 0; i < nr; i++) {
#ifdef CACHED_CONS_VISIBLE
	    if ((r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISREFPOS &&
		(r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)
#else
		if ((r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISREFPOS&&
		    (r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO &&
		    (r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISCONS)
#endif
		    r[i].y = j++;
		else
		    r[i].y = 0;
	}
	return 0;
    }

    /* Otherwise CSIR_ALLOCATE_Y_MULTIPLE */
#define xgap 3

    /* Create and initialise X and Y trees */
    struct xTREE axtree = SPLAY_INITIALIZER(&axtree);
    struct yTREE aytree = SPLAY_INITIALIZER(&aytree);
    struct yTREE iytree = SPLAY_INITIALIZER(&iytree);

    /* Compute Y coords */
    for (i = 0; i < nr; i++) {
#ifdef CACHED_CONS_VISIBLE
	if ((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS ||
	    (r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
#else
	    if ((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISREFPOS ||
		(r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO ||
		(r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISCONS)
#endif
		{
		    r[i].y = 0;
		    continue;
		}

	//xtree_dump(&axtree, NULL);
	//ytree_dump(&aytree, NULL);

#ifdef KEEP_Y
	/* If the node has a cached y coord already, use it if valid */
	if (r[i].y) {
	    struct xy_pair yline;

	    yline.x = 0;
	    yline.y = r[i].y;

	    /* Ensure we have sufficient yTREE entries */
	    if (r[i].y > yn) {
		int j;
		for (j = yn+1; j <= r[i].y; j++) {
		    node = (struct xy_pair *)malloc(sizeof(*node));
		    node->y = ++yn;
		    node->x = 0;
		    SPLAY_INSERT(yTREE, &iytree, node);
		}
	    }

	    /* See if this Y coord is available */
	    node = SPLAY_FIND(yTREE, &iytree, &yline);
	    if (node) {
		SPLAY_REMOVE(yTREE, &iytree, node);
		node->x = r[i].end + xgap;
		SPLAY_INSERT(xTREE, &axtree, node);
		SPLAY_INSERT(yTREE, &aytree, node);
		//printf("Used cached Y=%d\n", r[i].y);
		continue;
	    }

	    //printf("Looked for x=%d, y=%d, but not in yTREE (0..%d)\n",
	    //r[i].start, r[i].y, yn);

	    /*
	     * Faster equivalent to the above by using aytree
	     */
	    node = SPLAY_FIND(yTREE, &aytree, &yline);
	    if (node) {
		assert(node->y == r[i].y);
		if (r[i].start >= node->x) {
		    SPLAY_REMOVE(xTREE, &axtree, node);
		    node->x = r[i].end + xgap;
		    SPLAY_INSERT(xTREE, &axtree, node);
		    continue;
		}
	    }

	}
#endif

	if ((node = SPLAY_MIN(xTREE, &axtree)) != NULL && r[i].start >= node->x) {
	    //int try_cull = 0;

	    /* We found a node, but is there a smaller Y in the yTREE? */
	    curr = SPLAY_MIN(yTREE, &iytree);
	    if (curr && node->y > curr->y) {
		node = curr;
		r[i].y = node->y;
		SPLAY_REMOVE(yTREE, &iytree, node);
		node->x = r[i].end + xgap;
		SPLAY_INSERT(xTREE, &axtree, node);
		SPLAY_INSERT(yTREE, &aytree, node);
		//printf("1: alloc node=%p, x=%d, y=%d\n", node, node->x, node->y);
	    } else {
		/* Apparently not, what about smaller (in y) in xTREE? */
		curr = SPLAY_NEXT(xTREE, &axtree, node);
		while (curr && r[i].start >= curr->x) {
#ifdef KEEP_Y
		    if (curr->y == r[i].y) {
			/* previous Y, so keep it */
			node = curr;
			break;
		    }
#endif

		    next = SPLAY_NEXT(xTREE, &axtree, curr);
		    if (node->y > curr->y) {
			SPLAY_REMOVE(xTREE, &axtree, node);
			SPLAY_REMOVE(yTREE, &aytree, node);
			SPLAY_INSERT(yTREE, &iytree, node);
			//try_cull = 1;
			node = curr;
		    } else {
			SPLAY_REMOVE(xTREE, &axtree, curr);
			SPLAY_REMOVE(yTREE, &aytree, curr);
			SPLAY_INSERT(yTREE, &iytree, curr);
			//try_cull = 1;
		    }
		    curr = next;
		}

		r[i].y = node->y;
		SPLAY_REMOVE(xTREE, &axtree, node);
		node->x = r[i].end + xgap;
		SPLAY_INSERT(xTREE, &axtree, node);

		//printf("2: alloc node=%p, x=%d, y=%d\n", node, node->x, node->y);
	    }

#if 0
	    /* Cull Y tree if appropriate to remove excess rows */
	    if (try_cull) {
		for (curr = SPLAY_MAX(yTREE, &iytree); curr; curr = next) {
		    next = SPLAY_PREV(yTREE, &iytree, curr);
		    if (curr->y == yn) {
			SPLAY_REMOVE(yTREE, &iytree, curr);
			free(curr);
			yn--;
		    } else {
			break;
		    }
		}
	    }
#endif

	} else {
	    /* Check if we have a free y on ytree */
	    if ((node = SPLAY_MIN(yTREE, &iytree)) != NULL) {
		SPLAY_REMOVE(yTREE, &iytree, node);
	    } else {
		node = (struct xy_pair *)malloc(sizeof(*node));
		node->y = ++yn;
	    }
	    r[i].y = node->y;
	    node->x = r[i].end + xgap;
	    SPLAY_INSERT(xTREE, &axtree, node);
	    SPLAY_INSERT(yTREE, &aytree, node);

	    //printf("0: alloc node=%p, x=%d, y=%d\n", node, node->x, node->y);
	}
    }

    /* Delete trees */
    for (node = SPLAY_MIN(xTREE, &axtree); node; node = next) {
	next = SPLAY_NEXT(xTREE, &axtree, node);
	SPLAY_REMOVE(xTREE, &axtree, node);
	free(node); 
    }

    for (node = SPLAY_MIN(yTREE, &iytree); node; node = next) {
	next = SPLAY_NEXT(yTREE, &iytree, node);
	SPLAY_REMOVE(yTREE, &iytree, node);
	free(node); 
    }

    return 0;
}

/*
 * Given sequences with Y positions, this allocates tags to reside
 * in the same Y coordinate as the sequences they annotate.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int compute_ypos_tags(rangec_t *r, int nr) {
    int i;
    tg_rec key;
    HacheTable *h;
    HacheData hd;

    /* Check for at least one visible tag */
    for (i = 0; i < nr; i++) {
	if ((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISANNO)
	    break;
    }

    if (i == nr)
	return 0;

    /*
     * FIXME:
     *
     * A more optimal approach for this, given large numbers of seqs and
     * small numbers of tags, is to hash the tag ids by record (contents
     * being the r[index] index value) and then iterate through seqs
     * (assuming >= 1 tag found) matching them with the tag record id.
     * That way our hash table is of the size Ntags and not Nseqs.
     */


    /* Build an hash of sequence record numbers to Y coordinates. */
    h = HacheTableCreate(nr/3+1, HASH_DYNAMIC_SIZE);
    if (!h)
	return -1;

    h->name = "compute_ypos_tags()";

    for (i = 0; i < nr; i++) {
	if ((r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
	    continue;

	key = r[i].rec;
	hd.i = r[i].y;
	HacheTableAdd(h, (char *)&key, sizeof(key), hd, NULL);
    }

    /* Iterate through annotations copying over their Y coordinate */
    for (i = 0; i < nr; i++) {
	HacheItem *hi;

	if ((r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISANNO)
	    continue;

	key = r[i].pair_rec;
	if (!(hi = HacheTableSearch(h, (char *)&key, sizeof(key))))
	    /* consensus tag */
	    r[i].y = -1;
	else
	    r[i].y = hi->data.i;
    }

    HacheTableDestroy(h, 0);

    return 0;
}

static int contig_seqs_in_range2(GapIO *io, tg_rec bin_num,
				 int start, int end, int offset,
				 rangec_t **results, int *alloc, int used,
				 int complement, int mask, int val) {
    int count = used;
    int i, n, f_a, f_b;
    bin_index_t *bin = get_bin(io, bin_num);
    range_t *l;
    
    if (NULL == bin) return -1;

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

	    if (l->flags & GRANGE_FLAG_UNUSED)
		continue;

	    /* Select type: seqs, annotations, or maybe both (mask=val=0) */
	    if ((l->flags & mask) != val)
		continue;

	    if (NMAX(l->start, l->end) >= start
		&& NMIN(l->start, l->end) <= end) {
		int st, en;

		if (count >= *alloc) {
		    rangec_t *new_range;
		
		    *alloc = *alloc ? *alloc * 2 : 16;
		    new_range = (rangec_t *)realloc(*results,
						    *alloc * sizeof(rangec_t));
						   
		    if (new_range == NULL) {
		    	if (*results) free(*results);
			*results = NULL;
			*alloc = 0;
			cache_decr(io, bin);
			return -1;
		    } else {
		    	*results = new_range;
		    }
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
		(*results)[count].comp           = complement;
		(*results)[count].mqual          = l->mqual;
		(*results)[count].pair_rec       = l->pair_rec;
		(*results)[count].pair_start     = l->pair_start;
		(*results)[count].pair_end       = l->pair_end;
		(*results)[count].pair_mqual     = l->pair_mqual;
		(*results)[count].pair_contig    = l->pair_contig;
		(*results)[count].pair_timestamp = l->pair_timestamp;
		(*results)[count].library_rec    = l->library_rec;
		
		(*results)[count].flags = l->flags;
		(*results)[count].y = l->y;
		(*results)[count].orig_rec = bin->rec;
		(*results)[count].orig_ind = n;
		(*results)[count].seq_match = 0;
		(*results)[count].seq_hash  = 0;
		count++;
	    }
	}
    }

    /* Recurse down bins */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;
	ch = get_bin(io, bin->child[i]);
	if (!ch) {
	    verror(ERR_WARN, "contig_seqs_in_range2", "Bin #%"PRIrec
		   " failed to load", bin->child[i]);
	    continue;
	}

	if (end   >= NMIN(ch->pos, ch->pos + ch->size-1) &&
	    start <= NMAX(ch->pos, ch->pos + ch->size-1)) {
	    count = contig_seqs_in_range2(io, bin->child[i], start, end,
					  NMIN(ch->pos, ch->pos + ch->size-1),
					  results, alloc, count,
					  complement, mask, val);
	    if (count == -1) {
		cache_decr(io, bin);
		return -1;
	    }
	}
    }

    cache_decr(io, bin);
    return count;
}


/*
 * Writes back the computed Y position into the cached range arrays
 * in the tg_cache.c data structs. (If not still in cache, do nothing.)
 */
void update_range_y(GapIO *io, rangec_t *r, int count) {
    int i;
    tg_rec last_bin = -1;

    for (i = 0; i < count; i++) {
	bin_index_t *bin = NULL;
	range_t *rng;

	if (last_bin != r[i].orig_rec) {
	    last_bin  = r[i].orig_rec;
	    bin = cache_search_no_load(io, GT_Bin, r[i].orig_rec);
	}

	if (!bin)
	    continue;
	
	rng = arrp(range_t, bin->rng, r[i].orig_ind);
	assert(r[i].rec == rng->rec);

	rng->y = r[i].y;
    }
}


static rangec_t *contig_objects_in_range(GapIO *io, contig_t **c, seq_sort_t *sort_set,
					 int start, int end,
					 int first, int second,
					 int *count, int mask, int val) {

    rangec_t *r = malloc(16 * sizeof(rangec_t));
    int alloc = 16;
    
    if (NULL == r) {
	*count = -1;
	return NULL;
    }

    cache_incr(io, *c);
    *count = contig_seqs_in_range2(io, contig_get_bin(c), start, end,
				   contig_offset(io, c), &r, &alloc, 0, 0,
				   mask, val);
    if (*count <= 0) {
	cache_decr(io, *c);
	if (*count == 0) return r;
	if (NULL != r) free(r);
	return NULL;
    }

    if (r) {
    	int job;
    
    	if (first & CSIR_DEFAULT) {
	    first |= sort_set->p_sort;
	}
	
	if (second & CSIR_DEFAULT ) {
	    second |= sort_set->s_sort;
	}
	
	job = first | second;
    
	if (job & CSIR_PAIR) {
	    pair_rangec(io, (*c)->rec, r, *count);
	}
	
	if (job & (CSIR_DEFAULT | CSIR_SORT_BY_XEND | CSIR_SORT_BY_CLIPPED |
	    CSIR_SORT_BY_X | CSIR_SORT_BY_Y | CSIR_ALLOCATE_Y |
	    CSIR_SORT_BY_SEQ_TECH)) {
	    
	    if (job & CSIR_SORT_BY_SEQ_TECH) {
		int i;

		/* FIXME: Move seq_tech into bin range flags? */
		for (i = 0; i < *count; i++) {
		    if ((r[i].flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ) {
			seq_t *s = cache_search(io, GT_Seq, r[i].rec);
			r[i].seq_tech = s->seq_tech;
		    } else {
			r[i].seq_tech = -1;
		    }
		}
	    }
	    
	    if (job & CSIR_SORT_BY_SEQUENCE) {
	    	if (find_string_match_values(io, r, sort_set, *count)) {
		    // something wrong, possibly selection not set. Default to tech sort
		    first ^= CSIR_SORT_BY_SEQUENCE;
		    first |= CSIR_SORT_BY_X | CSIR_SORT_BY_SEQ_TECH;
		}
	    }
	    
	    sort_io = io;
	    
	    if (sort_set) {
	    	contig_set_base_sort_point(sort_set->base_pos);
	    }
	    
	    primary   = set_sort(first);
	    secondary = set_sort(second);
	    
	    qsort(r, *count, sizeof(*r), chosen_sort);

	    if (job & CSIR_ALLOCATE_Y) {
		if (job & CSIR_SORT_BY_SEQ_TECH) {
		    int i = 0, j, n = *count;
		    int last, start, max_y, shift = 0;

		    do {
			start = i;

			/* Find first real sequence (non-tag) seq_tech */
			for (last = r[i].seq_tech; last == -1 && i < n; i++)
			    last = r[i].seq_tech;

			/* Break into continuous technology fragment */
			while (i < n && (r[i].seq_tech == last ||
					 r[i].seq_tech == -1))
			       i++;

			/* Allocate Y and shift down */
			compute_ypos(&r[start], i - start, job & CSIR_ALLOCATE_Y);

			for (max_y = 0, j = start; j < i; j++) {
			    r[j].y += shift;
			    if (max_y < r[j].y)
				max_y = r[j].y;
			}
			shift = max_y+1;
		    } while (i < n);
		} else {
		    compute_ypos(r, *count, job & CSIR_ALLOCATE_Y);
		}

		compute_ypos_tags(r, *count);
	    }
    	
	    if (job & CSIR_SORT_BY_Y)
 	    	qsort(r, *count, sizeof(*r), sort_range_by_y); 
	}
    } else {
	// count -1 on memory error
	if (*count == -1) {
	    verror(ERR_WARN, "tg_contig", "Out of memory - unable to get objects\n");
	}
    }
    
    cache_decr(io, *c);
    return r;
}


rangec_t *contig_items_in_range(GapIO *io, contig_t **c, seq_sort_t *settings, int start, int end,
				int first_order, int second_order, int *count) {
    
    return contig_objects_in_range(io, c, settings, start, end, first_order, second_order, count, 0, 0);
}


rangec_t *contig_seqs_in_range(GapIO *io, contig_t **c, int start, int end,
			       int job, int *count) {
    
    return contig_objects_in_range(io, c, NULL, start, end, job, 0, count,
    	    GRANGE_FLAG_ISMASK, GRANGE_FLAG_ISSEQ);
}


rangec_t *contig_anno_in_range(GapIO *io, contig_t **c, int start, int end,
			       int job, int *count) {

    return contig_objects_in_range(io, c, NULL, start, end, job, 0, count,
				   GRANGE_FLAG_ISMASK, GRANGE_FLAG_ISANNO);
}

rangec_t *contig_refpos_in_range(GapIO *io, contig_t **c, int start, int end,
				 int job, int *count) {

    return contig_objects_in_range(io, c, NULL, start, end, job, 0, count,
 				   GRANGE_FLAG_ISMASK, GRANGE_FLAG_ISREFPOS);
}


/*
 * As per contig_seqs_in_range, but consensus sequences only.
 */
static int contig_cons_in_range2(GapIO *io, tg_rec bin_num,
				 int start, int end, int offset,
				 rangec_t **results, int *alloc, int used,
				 int complement) {
    int count = used;
    int i, n, f_a, f_b;
    bin_index_t *bin = get_bin(io, bin_num);
    range_t *l;
    int cst = INT_MIN, cend = INT_MIN;

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

	    if (l->flags & GRANGE_FLAG_UNUSED)
		continue;

	    if ((l->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISCONS)
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
		(*results)[count].pair_contig= 0;

		(*results)[count].flags = l->flags;
		(*results)[count].y = l->y;
		(*results)[count].orig_rec = bin->rec;
		(*results)[count].orig_ind = n;

		if ((*results)[count].start > cend) {
		    cst  = (*results)[count].start;
		    cend = (*results)[count].end;
		} else if ((*results)[count].end > cend) {
		    cend = (*results)[count].end;
		}

		count++;
	    }
	}
    }

    /* Recurse down bins */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;
	ch = get_bin(io, bin->child[i]);
	if (end   >= NMIN(ch->pos, ch->pos + ch->size-1) &&
	    start <= NMAX(ch->pos, ch->pos + ch->size-1)) {

	    /*
	     * We don't care about overlapping consensus pieces, so skip
	     * recursion if this portion is already accounted for.
	     */
	    if (NMIN(ch->pos, ch->pos + ch->size-1) >= cst &&
		NMAX(ch->pos, ch->pos + ch->size-1) <= cend)
		continue;

	    count = contig_cons_in_range2(io, bin->child[i], start, end,
					  NMIN(ch->pos, ch->pos + ch->size-1),
					  results, alloc, count,
					  complement);
	}
    }

    cache_decr(io, bin);
    return count;
}

/*
 * As per contig_seqs_in_range, but consensus sequences only.
 */
rangec_t *contig_cons_in_range(GapIO *io, contig_t **c, int start, int end,
			       int job, int *count) {
    rangec_t *r = NULL;
    int alloc = 0;

    cache_incr(io, *c);
    *count = contig_cons_in_range2(io, contig_get_bin(c), start, end,
				   contig_offset(io, c), &r, &alloc, 0, 0);
    cache_decr(io, *c);

    if (job & CSIR_SORT_BY_XEND)
	qsort(r, *count, sizeof(*r), sort_range_by_x_end);
    else if (job & (CSIR_SORT_BY_X | CSIR_SORT_BY_Y))
	qsort(r, *count, sizeof(*r), sort_range_by_x);

    return r;
}


static int contig_bins_in_range2(GapIO *io, tg_rec bin_num,
				 int start, int end, int offset,
				 rangec_t **results, int *alloc, int used,
				 int complement, int min_size,
				 int leaves_only) {
    int count = used;
    bin_index_t *bin = get_bin(io, bin_num);
    int i, f_a, f_b, leaf;

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

    /* Add bin */
    if (count >= *alloc) {
	*alloc = *alloc ? *alloc * 2 : 16;
	*results = (rangec_t *)realloc(*results, *alloc *sizeof(rangec_t));
    }

    (*results)[count].rec   = bin_num;
    (*results)[count].start = offset;
    (*results)[count].end   = offset + bin->size;
    (*results)[count].comp  = complement;

    /* To allow easy extraction of data from an individual bin */
    (*results)[count].pair_start = f_a;
    (*results)[count].pair_end   = f_b;

    if (!leaves_only)
	count++;


    /* Items in this case are child bins rather than sequences */
    leaf = 0;
    for (i = 0; bin->size >= min_size && i < 2; i++) {
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

	    /* Recurse too */
	    leaf++;
	    count = contig_bins_in_range2(io, bin->child[i], start, end,
					  NMIN(ch->pos, ch->pos + ch->size-1),
					  results, alloc, count, complement,
					  min_size, leaves_only);
	}
    }

    cache_decr(io, bin);

    if (leaves_only && !leaf)
	return ++count;

    return count;
}

rangec_t *contig_bins_in_range(GapIO *io, contig_t **c, int start, int end,
			       int job, int min_size, int *count) {
    rangec_t *r = NULL;
    int alloc = 0;

    cache_incr(io, *c);
    *count = contig_bins_in_range2(io, contig_get_bin(c), start, end,
				   contig_offset(io, c),
				   &r, &alloc, 0, 0, min_size,
				   job & CSIR_LEAVES_ONLY);
    cache_decr(io, *c);
    
    if (job & CSIR_SORT_BY_XEND)
	qsort(r, *count, sizeof(*r), sort_range_by_x_end);
    else if (job & CSIR_SORT_BY_X)
	qsort(r, *count, sizeof(*r), sort_range_by_x);

    return r;
}


static void contig_bin_dump2(GapIO *io, tg_rec bin_num, int level) {
    int i;
    bin_index_t *bin = get_bin(io, bin_num);

    cache_incr(io, bin);

    printf("%*sBin %"PRIrec"\tpos %d\tsize %d\tleft %"PRIrec
	   "\tright %"PRIrec"\tflags %d\n", 
	   level*4, "", bin->rec, bin->pos, bin->size,
	   bin->child[0], bin->child[1], bin->flags);

    for (i = 0; i < 2; i++) {
	if (!bin->child[i])
	    continue;
	contig_bin_dump2(io, bin->child[i], level+1);
    }

    cache_decr(io, bin);
}

void contig_bin_dump(GapIO *io, tg_rec cnum) {
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, cnum);
    contig_bin_dump2(io, contig_get_bin(&c), 0);
}

/* ---------------------------------------------------------------------- */
/* iterators on a contig to allow scanning back and forth */

/* Block size for iterator fetches, first set optionally small */
#define CITER_BS 10000
#define CITER_bs 100

/*
 * Similar to contig_bins_in_range(), but checking for the first bin to
 * contain data of a specific type rather than to overlap a specific
 * range.
 *
 * Given start we try to find the first bin with coordinate >= start that
 * contains somewhere within it an object of 'type'. This allows us to
 * rapidly search for needles in haystacks (eg refpos markers in an
 * illumina assembly) by skipping down the bin tree rather than linearly
 * scanning through all objects.
 *
 * It does this using the nseqs, nrefpos or nannos fields in bin_index_t,
 * which act as accumulative figures for the number of these object types
 * in this bin plus children.
 *
 * Returns the next starting position, or best_start if none found.
 */
static int range_next_by_type2(GapIO *io, contig_t *c, int bin_num,
			       int offset, int comp, int start, int type,
			       int best_start) {
    bin_index_t *bin = get_bin(io, bin_num);
    int i, f_a, f_b, n_in_this;
    int child_pos[2], order[2], child_total;

    switch (type) {
    case GRANGE_FLAG_ISSEQ:
	if (bin->nseqs <= 0)
	    return best_start;
	n_in_this = bin->nseqs;
	break;

    case GRANGE_FLAG_ISANNO:
	if (bin->nanno <= 0)
	    return best_start;
	n_in_this = bin->nanno;
	break;

    case GRANGE_FLAG_ISREFPOS:
	if (bin->nrefpos <= 0)
	    return best_start;
	n_in_this = bin->nrefpos;
	break;

    case GRANGE_FLAG_ISANY:
	if (bin->nseqs + bin->nanno + bin->nrefpos <= 0)
	    return best_start;
	n_in_this = bin->nseqs + bin->nanno + bin->nrefpos;
	break;

    default:
	fprintf(stderr, "Unknown type given to range_next_by_type2\n");
	return best_start;
    }

    cache_incr(io, bin);
    if (bin->flags & BIN_COMPLEMENTED) {
	comp ^= 1;
    }

    if (comp) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    /* Determine left to right order of children (varies if complemented) */
    /* Also compute sum of data in children while we're at it */
    child_total = 0;
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;

	if (!bin->child[i]) {
	    child_pos[i] = INT_MAX;
	    continue;
	}
	ch = get_bin(io, bin->child[i]);
	if (!ch) {
	    verror(ERR_FATAL, "range_next_by_type2",
		   "Failed to load bin #%"PRIrec", child of bin #%"PRIrec,
		   bin->child[i], bin->rec);
	    continue;
	}

	child_pos[i] = NMIN(ch->pos,  ch->pos + ch->size-1);

	switch (type) {
	case GRANGE_FLAG_ISSEQ:
	    child_total += ch->nseqs;
	    break;

	case GRANGE_FLAG_ISANNO:
	    child_total += ch->nanno;
	    break;

	case GRANGE_FLAG_ISREFPOS:
	    child_total += ch->nrefpos;
	    break;

	case GRANGE_FLAG_ISANY:
	    child_total += ch->nseqs + ch->nanno + ch->nrefpos;
	    break;
	}
    }
    /* Need to generalise if we have > 2 children ever */
    if (child_pos[0] < child_pos[1]) {
	order[0] = 0;
	order[1] = 1;
    } else {
	order[0] = 1;
	order[1] = 0;
    }

    if (!bin_empty(bin) && 
	best_start > NMIN(bin->start_used, bin->end_used) &&
	start     <= NMAX(bin->start_used, bin->end_used)) {
	/* Bin overlaps, but possibly the items are in children */
	if (n_in_this > child_total)
	    best_start = NMIN(bin->start_used, bin->end_used);
    }

    if (best_start <= start) {
	cache_decr(io, bin);
	return start;
    }

    /* Recurse down bins */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[order[i]])
	    continue;

	ch = get_bin(io, bin->child[order[i]]);
	if (!ch) {
	    verror(ERR_FATAL, "range_next_by_type2",
		   "Failed to load bin #%"PRIrec", child of bin #%"PRIrec,
		   bin->child[order[i]], bin->rec);
	    continue;
	}

	switch (type) {
	case GRANGE_FLAG_ISSEQ:
	    if (ch->nseqs <= 0)
		continue;
	    break;

	case GRANGE_FLAG_ISANNO:
	    if (ch->nanno <= 0)
		continue;
	    break;

	case GRANGE_FLAG_ISREFPOS:
	    if (ch->nrefpos <= 0)
		continue;
	    break;

	case GRANGE_FLAG_ISANY:
	    if (ch->nseqs + ch->nanno + ch->nrefpos <= 0)
		continue;
	    break;
	}

	if (NMAX(ch->pos, ch->pos + ch->size-1) < start)
	    continue;

	if (NMIN(ch->pos, ch->pos + ch->size-1) < best_start) {
	    int b;

	    b = range_next_by_type2(io, c, bin->child[order[i]],
				    NMIN(ch->pos, ch->pos + ch->size-1),
				    comp, start, type, best_start);

	    if (best_start > b)
		best_start = b;
	}
    }

    cache_decr(io, bin);

    return best_start;
}

/*
 * See comments above range_next_by_type2() for the algorithm details.
 *
 * We return a new coordinate >= start such that there are no objects of
 * a given type between this coordinate and start - essentially skipping
 * through the contig faster than a linear scan could achieve.
 */
static int range_next_by_type(GapIO *io, tg_rec cnum, int start, int type) {
    contig_t *c = cache_search(io, GT_Contig, cnum);
    int ret;

    cache_incr(io, c);
    ret = range_next_by_type2(io, c,
			      contig_get_bin(&c), contig_offset(io, &c),
			      0 /* comp */, start, type, INT_MAX);
    cache_decr(io, c);

    return ret;
}

static int range_prev_by_type2(GapIO *io, contig_t *c, int bin_num,
			       int offset, int comp, int start, int type,
			       int best_start) {
    bin_index_t *bin = get_bin(io, bin_num);
    int i, f_a, f_b, n_in_this;
    int child_pos[2], order[2], child_total;

    switch (type) {
    case GRANGE_FLAG_ISSEQ:
	if (bin->nseqs <= 0)
	    return best_start;
	n_in_this = bin->nseqs;
	break;

    case GRANGE_FLAG_ISANNO:
	if (bin->nanno <= 0)
	    return best_start;
	n_in_this = bin->nanno;
	break;

    case GRANGE_FLAG_ISREFPOS:
	if (bin->nrefpos <= 0)
	    return best_start;
	n_in_this = bin->nrefpos;
	break;

    case GRANGE_FLAG_ISANY:
	if (bin->nseqs + bin->nanno + bin->nrefpos <= 0)
	    return best_start;
	n_in_this = bin->nseqs + bin->nanno + bin->nrefpos;
	break;

    default:
	fprintf(stderr, "Unknown type given to range_prev_by_type2\n");
	return best_start;
    }

    cache_incr(io, bin);
    if (bin->flags & BIN_COMPLEMENTED) {
	comp ^= 1;
    }

    if (comp) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

    /* Determine left to right order of children (varies if complemented) */
    /* Also compute sum of data in children while we're at it */
    child_total = 0;
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;

	if (!bin->child[i]) {
	    child_pos[i] = INT_MAX;
	    continue;
	}
	ch = get_bin(io, bin->child[i]);
	child_pos[i] = NMIN(ch->pos,  ch->pos + ch->size-1);

	switch (type) {
	case GRANGE_FLAG_ISSEQ:
	    child_total += ch->nseqs;
	    break;

	case GRANGE_FLAG_ISANNO:
	    child_total += ch->nanno;
	    break;

	case GRANGE_FLAG_ISREFPOS:
	    child_total += ch->nrefpos;
	    break;

	case GRANGE_FLAG_ISANY:
	    child_total += ch->nseqs + ch->nanno + ch->nrefpos;
	    break;
	}
    }
    /* Need to generalise if we have > 2 children ever */
    if (child_pos[0] < child_pos[1]) {
	order[0] = 1;
	order[1] = 0;
    } else {
	order[0] = 0;
	order[1] = 1;
    }

    if (!bin_empty(bin) && 
	best_start < NMAX(bin->start_used, bin->end_used) &&
	start     >= NMIN(bin->start_used, bin->end_used)) {
	/* Bin overlaps, but possibly the items are in children */
	if (n_in_this > child_total)
	    best_start = NMAX(bin->start_used, bin->end_used);
    }

    if (best_start >= start) {
	cache_decr(io, bin);
	return start;
    }

    /* Recurse down bins */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[order[i]])
	    continue;

	ch = get_bin(io, bin->child[order[i]]);

	switch (type) {
	case GRANGE_FLAG_ISSEQ:
	    if (ch->nseqs <= 0)
		continue;
	    break;

	case GRANGE_FLAG_ISANNO:
	    if (ch->nanno <= 0)
		continue;
	    break;

	case GRANGE_FLAG_ISREFPOS:
	    if (ch->nrefpos <= 0)
		continue;
	    break;

	case GRANGE_FLAG_ISANY:
	    if (ch->nseqs + ch->nanno + ch->nrefpos <= 0)
		continue;
	    break;
	}

	if (NMIN(ch->pos, ch->pos + ch->size-1) > start)
	    continue;

	if (NMAX(ch->pos, ch->pos + ch->size-1) > best_start) {
	    int b;

	    b = range_prev_by_type2(io, c, bin->child[order[i]],
				    NMIN(ch->pos, ch->pos + ch->size-1),
				    comp, start, type, best_start);

	    if (best_start < b)
		best_start = b;
	}
    }

    cache_decr(io, bin);

    return best_start;
}

/*
 * See comments above range_next_by_type2() for the algorithm details.
 *
 * We return a new coordinate <= start such that there are no objects of
 * a given type between this coordinate and start - essentially skipping
 * through the contig faster than a linear scan could achieve.
 */
static int range_prev_by_type(GapIO *io, tg_rec cnum, int start, int type) {
    contig_t *c = cache_search(io, GT_Contig, cnum);
    int ret;

    cache_incr(io, c);
    ret = range_prev_by_type2(io, c,
			      contig_get_bin(&c), contig_offset(io, &c),
			      0 /* comp */, start, type, INT_MIN);
    cache_decr(io, c);

    return ret;
}

/*
 * Populates a sequence iterator.
 * Returns 0 on success
 *        -1 on failure (typically fallen off the end of the contig)
 */
static int range_populate(GapIO *io, contig_iterator *ci,
			  tg_rec cnum, int start, int end) {
    contig_t *c = (contig_t *)cache_search(io, GT_Contig, cnum);
    
    if (NULL == c) return -1;

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

	//if (start < c->start)
	//    start = c->start;
	//if (end > c->end)
	//    end = c->end;
    }

    ci->cnum = cnum;
    ci->start = start;
    ci->end = end;
    if (ci->type == GRANGE_FLAG_ISSEQ)
	ci->r = contig_seqs_in_range(io, &c, start, end,
				     ci->sort_mode, &ci->nitems);
    else if (ci->type == GRANGE_FLAG_ISANNO)
	ci->r = contig_anno_in_range(io, &c, start, end,
				     ci->sort_mode, &ci->nitems);
    else if (ci->type == GRANGE_FLAG_ISREFPOS)
	ci->r = contig_refpos_in_range(io, &c, start, end,
				       ci->sort_mode, &ci->nitems);
    else
	ci->r = contig_items_in_range(io, &c, NULL, start, end,
				      ci->sort_mode, 0, &ci->nitems);
    if (NULL == ci->r) return -1;

    /*
     * If we've requested a CLIPPED_START or CLIPPED_END then filter out
     * any sequences that when clipped do not overlap out start..end range.
     *
     * The reason is that these will be brought in during a subsequent
     * 
     */
    if (ci->sort_mode == (CSIR_SORT_BY_X | CSIR_SORT_BY_CLIPPED)) {
	int i, j;
	for (i = j = 0; i < ci->nitems; i++) {
	    int c_start;
	    seq_t *s;

	    if ((ci->r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ) {
		ci->r[j++] = ci->r[i];
		continue;
	    }

	    s = cache_search(io, GT_Seq, ci->r[i].rec);
	    if (NULL == s) return -1;
	    if ((s->len < 0) ^ ci->r[i].comp) {
		c_start = ci->r[i].start + ABS(s->len) - (s->right-1) - 1;
	    } else {
		c_start = ci->r[i].start + s->left-1;
	    }
	    if (c_start > end)
		continue;
	    ci->r[j++] = ci->r[i];
	}
	ci->nitems = j;
    }

    if (ci->sort_mode == (CSIR_SORT_BY_XEND | CSIR_SORT_BY_CLIPPED)) {
	int i, j;
	for (i = j = 0; i < ci->nitems; i++) {
	    int c_end;
	    seq_t *s;

	    if ((ci->r[i].flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ) {
		ci->r[j++] = ci->r[i];
		continue;
	    }

	    s = cache_search(io, GT_Seq, ci->r[i].rec);
	    if (NULL == s) return -1;
	    if ((s->len < 0) ^ ci->r[i].comp) {
		c_end = ci->r[i].start + ABS(s->len) - (s->left-1) - 1;
	    } else {
		c_end = ci->r[i].start + s->right-1;
	    }
	    if (c_end < start)
		continue;
	    ci->r[j++] = ci->r[i];
	}
	ci->nitems = j;
    }

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
 *
 * 'whence' may be either CITER_FIRST or CITER_LAST and it controls whether
 * we start the iteration point from the beginning or end of the list.
 * In addition to this we can add CITER_ISTART or CITER_IEND to whence to
 * control whether we wish to sort our data by item start or end coordinate.
 * CITER_ISTART is the default mode.  (This is particularly useful if we wish
 * to step through sequences in reverse order based on when they become
 * visible.) Note this is the actual start/end and not the clipped good
 * quality portion only. For that use CITER_ICLIPPEDSTART and
 * CITER_ICLIPPEDEND instead.
 *
 * The start and end parameters dictate the initial region to query. We
 * may specify them as either coordinates or use CITER_CSTART and CITER_CEND
 * as synonyms for the first and last coordinate in the contig.
 *
 * 'type' can be either GRANGE_FLAG_ISSEQ or GRANGE_FLAG_ISANNO to only
 * iterate around data of that specific type, or GRANGE_FLAG_ISANY to
 * iterate around all data.
 *
 * Finally auto_extend controls whether the start..end range is just a
 * location to start iterating from (auto_extend == 1) or a hard limit
 * with no desire to iterate past that range (auto_extend == 0).
 *
 * Returns contig_iterator pointer on success
 *         NULL on failure
 */
contig_iterator *contig_iter_new_by_type(GapIO *io, tg_rec cnum,
					 int auto_extend, int whence,
					 int start, int end, int type) {
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
    ci->first_r = 1;
    ci->type = type;
    switch (whence & CITER_SE_MASK) {
    case CITER_ISTART:
	ci->sort_mode = CSIR_SORT_BY_X;
	break;
    case CITER_IEND:
	ci->sort_mode = CSIR_SORT_BY_XEND;
	break;
    case CITER_ICLIPPEDSTART:
	ci->sort_mode = CSIR_SORT_BY_X | CSIR_SORT_BY_CLIPPED;
	break;
    case CITER_ICLIPPEDEND:
	ci->sort_mode = CSIR_SORT_BY_XEND | CSIR_SORT_BY_CLIPPED;
	break;
    }

    /*
     * Generic bug fix - sorry! Expand c->start/end by +/- a small amount
     * so that when we've been inserting or deleting bases we correctly
     * include this data. The problem here is we often use iterators to find
     * the contig extents after munging the contig in some way so we can
     * correct c->start/c->end, but this is a catch-22 situation.
     *
     * Ideally we should pick an outer limit to iterate over by looking at
     * bin dimensions, and then work in from there.
     */
    ci->cstart = start == CITER_CSTART ? c->start - 50 : start;
    ci->cend   =   end == CITER_CEND   ? c->end   + 50 : end;

    if ((whence & CITER_FL_MASK) == CITER_FIRST) {
	start = ci->cstart;
	end   = start + ((whence & CITER_SMALL_BS) ? CITER_bs-1 : CITER_BS-1);
    } else {
	end   = ci->cend;
	start = end - ((whence & CITER_SMALL_BS) ? CITER_bs-1 : CITER_BS-1);
    }

    if (0 != range_populate(io, ci, cnum, start, end)) {
	contig_iter_del(ci);
	return NULL;
    }

    if ((whence & CITER_FL_MASK) == CITER_LAST) {
	ci->index = ci->nitems-1;
    }

    return ci;
}


/*
 * The old interface, which is hardcoded to only iterate through sequences.
 */
contig_iterator *contig_iter_new(GapIO *io, tg_rec cnum, int auto_extend,
				 int whence, int start, int end) {
    return contig_iter_new_by_type(io, cnum, auto_extend, whence,
				   start, end, GRANGE_FLAG_ISSEQ);
}

/*
 * Returns seq range_t struct pointer on success
 *        NULL on failure (typically fallen off the end of the contig)
 */
rangec_t *contig_iter_next(GapIO *io, contig_iterator *ci) {
    rangec_t *r = NULL;

    for (;;) {
	while (ci->index >= ci->nitems) {
	    int candidate;

	    /* Fallen off the range edge */
	    //	if (!ci->auto_extend)
	    //	    return NULL;

	    if (ci->end >= ci->cend)
		return NULL;

	    candidate = range_next_by_type(io, ci->cnum, ci->end+1, ci->type);

	    if (-1 == range_populate(io, ci, ci->cnum,
	    			     candidate, candidate + CITER_BS-1))
				     //ci->start + CITER_BS, ci->end + CITER_BS))
		return NULL;

	    ci->index = 0;
	    ci->first_r = 0;
	}

	if (ci->nitems == 0)
	    return NULL;

	/*
	 * The first range query we start from the first item, even if it's
	 * partially overlapping. Subsequent queries we start from the first
	 * item completely in our next range.
	 */
	while (ci->index < ci->nitems && (r = &ci->r[ci->index++])) {
	    if (r->start >= ci->start)
		return r;
	    else if (ci->first_r && r->end >= ci->start)
		return r;
	}
    }
}

/*
 * Returns seq range_t struct pointer on success
 *        NULL on failure (typically fallen off the end of the contig)
 */
rangec_t *contig_iter_prev(GapIO *io, contig_iterator *ci) {
    rangec_t *r;

    for (;;) {
	while (ci->index < 0 || ci->nitems == 0) {
	    int candidate;

	    /* Fallen off the range edge */
	    //	if (!ci->auto_extend)
	    //	    return NULL;

	    if (ci->start <= ci->cstart)
		return NULL;
	    
	    candidate = range_prev_by_type(io, ci->cnum, ci->start, ci->type);
	    
	    if (-1 == range_populate(io, ci, ci->cnum,
				     candidate - CITER_BS+1, candidate))
	    //if (-1 == range_populate(io, ci, ci->cnum,
	    //ci->start - CITER_BS, ci->end - CITER_BS))
		return NULL;
	    

	    ci->index = ci->nitems-1;
	    ci->first_r = 0;
	}

	while (ci->index >= 0 && (r = &ci->r[ci->index--])) {
	    if (r->end <= ci->end)
		return r;
	    else if (ci->first_r && r->start <= ci->end)
		return r;
	}
    }
}

/*
 * Given an iterator from start..end we'll find sequences that may cover
 * e_start..e_end where e_start and e_end may possibly be larger than start
 * to end. Pictorially:
 *
 *                start      end
 *A-------------- |          |
 *B      -----    |          |
 *C   ---------------------  |
 *D   |         -----------------------------
 *E   |                        --------     |
 *F   |                            ---------------------
 *    |                                     |
 *    e_start                               e_end
 *
 * The original query can fetch back seqs C & D, but annotations entirely
 * outside this could be missed if we're doing GRANGE_FLAG_ISANY queries.
 *
 * So we expand the range to e_start to e_end to pick up annotations.
 * NOTE: the caller then need to manually filter the start/end range
 * itself to avoid then picking up seqs B and E which aren't in the
 * original range.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int iterator_expand_range(GapIO *io, tg_rec crec, int start, int end,
			  int *e_start, int *e_end) {
    int nr;
    rangec_t *r;
    contig_t *c = cache_search(io, GT_Contig, crec);

    if (!c)
	return -1;

    cache_incr(io, c);

    if (e_start) {
	int i;
	r = contig_seqs_in_range(io, &c, start, start, 0, &nr);
	if (NULL == r) {
	    cache_decr(io, c);
	    return -1;
	}
	*e_start = start;

	for (i = 0; i < nr; i++) {
	    if (*e_start > r[i].start)
		*e_start = r[i].start;
	}
	free(r);
    }

    if (e_end) {
	int i;
	r = contig_seqs_in_range(io, &c, end, end, 0, &nr);
	if (NULL == r) {
	    cache_decr(io, c);
	    return -1;
	}
	*e_end = end;

	for (i = 0; i < nr; i++) {
	    if (*e_end < r[i].end)
		*e_end = r[i].end;
	}
	free(r);
    }

    cache_decr(io, c);
    return 0;
}

/* ---------------------------------------------------------------------------
 * Track implementation
 */

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
static int contig_get_track2(GapIO *io, tg_rec bin_num, int start, int end,
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

    printf("%*scontig_get_track2 want %5.1f got ~%5.1f bin %"PRIrec
	   " (%d+%d) l=%"PRIrec" r=%"PRIrec"\n",
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
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i]) {
	    /* No data available, so fill with zero values */
	    int len, nitems, j, offset2;

	    if (i == 0 && bin->child[1]) {
		ch = get_bin(io, bin->child[1]);
		len = bin->size - ch->size;
		offset2 = 0;
	    } else if (i == 1 && bin->child[0]) {
		ch = get_bin(io, bin->child[0]);
		len = bin->size - ch->size;
		offset2 = ch->size;
	    } else {
		len = bin->size;
		offset2 = 0;
	    }
	    nitems = len / bpv;

	    for (j = 0; j < nitems; j++) {
		while (count+1 > *alloc) {
		    *alloc = *alloc ? *alloc * 2 : 16;
		    *tv = (tvalues_t *)realloc(*tv, *alloc * sizeof(**tv));
		}
		(*tv)[count].pos = NORM((double)j * len / nitems + offset2);
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
    int bin_off;
    int *data, *data3;
    bin_index_t *start_bin;
    tg_rec start_bin_rec;

    printf("Query %d..%d bpv %f\n", start, end, bpv);

    nele = ceil((end-start+1)/bpv);
    bpv = (end-start+1)/nele;
    t = track_create_fake(type, nele);
    data = ArrayBase(int, t->data);

    start_bin = bin_for_range(io, c, start, end, 0, &bin_off, NULL);
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
 * Does NOT destroy any data within it, including the root bin.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int contig_destroy(GapIO *io, tg_rec rec) {
    int i, j;
    contig_t *c;

    /* Delete from the btree index */
    if (!(c = cache_search(io, GT_Contig, rec)))
	return -1;
    if (NULL == (c = cache_rw(io, c))) return -1;

    if (c->name) {
	tg_rec r = io->iface->contig.index_del(io->dbh, c->name, rec);
	if (r != -1 && r != io->db->contig_name_index) {
	    io->db = cache_rw(io, io->db);
	    io->db->contig_name_index = r;
	}
    }

    /* Remove from contig order */
    io->contig_order = cache_rw(io, io->contig_order);
    io->db = cache_rw(io, io->db);

    for (i = j = 0; i < io->db->Ncontigs; i++) {
	if (arr(tg_rec, io->contig_order, i) == rec) {
	    continue;
	}
	arr(tg_rec, io->contig_order, j) = arr(tg_rec, io->contig_order, i);
	j++;
    }

    if (i == j) {
	fprintf(stderr, "Attempted to remove unknown contig, rec %"PRIrec"\n",
		rec);
	return -1;
    }

    io->db->Ncontigs--;
    ArrayMax(io->contig_order)--;

    /* Remove from registration system */
    contig_register_delete(io, rec);

    /* And finally, deallocate on disk to allow record to be reused */
    cache_deallocate(io, c);

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

#define DRAW_LINKS 1
#define DRAW_BINS  2
#define DRAW_SEQS  4

static int bin_dump_ps_recurse(FILE *gv, GapIO *io, tg_rec rec,
			       int (*levels)[2], int nlevels,
			       int depth_first, int to_draw) {
    int empty;
    int i;
    int complement;
    int offset;
    int end;
    int level;
    bin_list *head = NULL;
    bin_list *tail = NULL;
    bin_list *item;
    bin_list *kids[2];
    int       nkids;
    int       kid_offsets[2];
    bin_index_t *bin = NULL;
    bin_index_t *kid = NULL;
    
    for (i = 0; i < nlevels; i++) {
	levels[i][0] =  INT_MAX;
	levels[i][1] = -INT_MAX;
    }
    
    head = malloc(sizeof(bin_list));
    if (!head) return -1;
    
    head->rec = rec;
    head->parent_comp   = 0;
    head->parent_offset = 0;
    head->parent_size   = 0;
    head->parent_level  = -1;
    head->child_index   = 0;
    head->next = NULL;
    if (!depth_first) tail = head;
    
    while (head) {
	complement = head->parent_comp;
	
	bin = get_bin(io, head->rec);
	if (bin) {
	    if (bin->flags & BIN_COMPLEMENTED) {
		complement ^= 1;
	    }
	    offset = (head->parent_comp
		      ? (head->parent_offset + head->parent_size
			 - bin->pos - bin->size)
		      : head->parent_offset + bin->pos);
	    end = offset + bin->size;
	    
	    for (level = head->parent_level + 1; level < nlevels; level++) {
		if (offset >= levels[level][1] || end <= levels[level][0]) {
		    break;
		}
	    }
	    assert(level < nlevels);
	    
	    if (levels[level][0] > offset) levels[level][0] = offset;
	    if (levels[level][1] < end)    levels[level][1] = end;
	    
	    empty = bin_empty(bin);
	    
	    if ((to_draw & DRAW_LINKS) && head->parent_level >= 0) {
		/* Linking lines */
		fprintf(gv, "%d %d %d %d %d l\n",
			head->parent_offset + head->parent_size / 2,
			head->parent_level,
			offset + bin->size / 2, level, head->child_index);
	    }
	    
	    /* The bin itself */
	    if (to_draw & DRAW_BINS) {
		fprintf(gv, "%d %d (%"PRIrec") %d %s %d %d %s %d %d sb ",
			bin->nrefpos, bin->nseqs, bin->rec, bin->size,
			complement ? "t" : "f",
			empty ? 0 : bin->start_used,
			empty ? 0 : bin->end_used,
			empty ? "f" : "t",
			level, offset);
		
		if ((to_draw & DRAW_SEQS) && bin->rng && ArrayMax(bin->rng)) {
		    /* Bin contents */
		    int n = ArrayMax(bin->rng);
		    int last_type  = -1;
		    int draw_count = 0;
		    int dup_count = 0;
		    int last_start =  0;
		    int last_end   = -1;
		    int skipped = 0;
		    /* fprintf(gv, "%d %d sqb\n", level, offset); */
		    fprintf(gv, "%d sqb\n", n);
		    for (i = 0; i < n; i++) {
			range_t *r = arrp(range_t, bin->rng, i);
			
			if (dup_count && ((r->flags & GRANGE_FLAG_UNUSED)
					  || (r->flags & GRANGE_FLAG_ISMASK) != last_type
					  || r->start != last_start
					  || r->end != last_end
					  || draw_count > 200)) {
			    if (dup_count > 1) {
				fprintf(gv, "%d { ", dup_count);
			    }
			    fprintf(gv, "%d %d sq%s", last_start, last_end,
				    dup_count > 1 ? " } repeat\n" : "\n");
			    dup_count = 0;
			    if (draw_count > 200) {
				fprintf(gv, "stroke\n");
				draw_count = 0;
			    }
			}
			
			if (r->flags & GRANGE_FLAG_UNUSED) {
			    skipped++;
			    continue;
			}
			if (skipped) {
			    fprintf(gv, "%d add\n", skipped);
			    skipped = 0;
			}
			if ((r->flags & GRANGE_FLAG_ISMASK) != last_type) {
			    last_type = r->flags & GRANGE_FLAG_ISMASK;
			    fprintf(gv, "%d sqg ",
				    last_type == GRANGE_FLAG_ISSEQ ? 1 : 0);
			    draw_count = 0;
			}
			last_start = r->start;
			last_end   = r->end;
			dup_count++;
			draw_count++;
		    }
		    if (dup_count) {
			if (dup_count > 1) { fprintf(gv, "%d { ", dup_count); }
			fprintf(gv, "%d %d sq%s", last_start, last_end,
				dup_count > 1 ? " } repeat\n" : "\n");
		    }
		    fprintf(gv, "sqe\n");
		}
		
		fprintf(gv, "eb\n");
	    }
	    
	    for (i = 0, nkids = 0; i < 2; i++) {
		if (!bin->child[i]) continue;
		kid = get_bin(io, bin->child[i]);
		if (!kid) {
		    fprintf(stderr,
			    "Couldn't get bin %"PRIrec"\n", bin->child[i]);
		    continue;
		}
		item = malloc(sizeof(bin_list));
		if (!item) goto unwind;
		item->rec = bin->child[i];
		item->parent_comp   = complement;
		item->parent_offset = offset;
		item->parent_size   = bin->size;
		item->parent_level  = level;
		item->child_index   = i;
		kid_offsets[nkids] = (complement
				      ? (offset+bin->size-kid->pos-kid->size)
				      : offset + kid->pos);
		kids[nkids++] = item;
	    }
	    if (nkids == 2
		&& (depth_first
		    ? kid_offsets[0] > kid_offsets[1]
		    : kid_offsets[0] < kid_offsets[1])) {
		/* Do rightmost child first */
		item = kids[0]; kids[0] = kids[1]; kids[1] = item;
	    } 
	    for (i = 0; i < nkids; i++) {
		if (depth_first) {
		    /* Put item one down in the list,
		       as head will be removed later */
		    kids[i]->next = head->next;
		    head->next = kids[i];	    
		} else {
		    kids[i]->next = NULL;
		    tail->next = kids[i];
		    tail = kids[i];
		}
	    }
	} else {
	    fprintf(stderr, "Couldn't get bin %"PRIrec"\n", head->rec);
	}
	item = head;
	head = head->next;
	free(item);
    }
    
    return 0;
 unwind:
    while (head) {
	item = head;
	head = head->next;
	free(item);
    }
    return -1;
}

static int count_levels_recurse(GapIO *io, tg_rec rec,
				int (**levels)[2], int *nlevels,
				int depth_first) {
    int maxlevel = 0;
    int i;
    int complement;
    int offset;
    int end;
    int level;
    bin_list *head = NULL;
    bin_list *tail = NULL;
    bin_list *item;
    bin_list *kids[2];
    int       nkids;
    int       kid_offsets[2];
    bin_index_t *bin = NULL;
    bin_index_t *kid = NULL;

    head = malloc(sizeof(bin_list));
    if (!head) return -1;

    head->rec = rec;
    head->parent_comp   = 0;
    head->parent_offset = 0;
    head->parent_size   = 0;
    head->parent_level  = -1;
    head->child_index   = 0;
    head->next = NULL;
    if (!depth_first) tail = head;

    while (head) {
	complement = head->parent_comp;

	bin = get_bin(io, head->rec);
	if (bin) {
	    if (bin->flags & BIN_COMPLEMENTED) {
		complement ^= 1;
	    }
	    offset = (head->parent_comp
		      ? head->parent_offset + head->parent_size - bin->pos - bin->size
		      : head->parent_offset + bin->pos);
	    end = offset + bin->size;

	    for (level = head->parent_level + 1; ; level++) {
		if (level == *nlevels) {
		    int new_n = *nlevels * 2;
		    int (*new_l)[2] = realloc(*levels, sizeof(**levels) * new_n);
	  
		    if (!new_l) {
			perror("count_levels_recurse");
			goto unwind;
		    }
		    for (i = *nlevels; i < new_n; i++) {
			new_l[i][0] =  INT_MAX;
			new_l[i][1] = -INT_MAX;
		    }
		    *levels = new_l;
		    *nlevels = new_n;
		}
		if (offset >= (*levels)[level][1] || end <= (*levels)[level][0]) break;
	    }

	    if ((*levels)[level][0] > offset) (*levels)[level][0] = offset;
	    if ((*levels)[level][1] < end)    (*levels)[level][1] = end;
	    if (level > maxlevel) maxlevel = level;
      
	    for (i = 0, nkids = 0; i < 2; i++) {
		if (!bin->child[i]) continue;
		kid = get_bin(io, bin->child[i]);
		if (!kid) {
		    fprintf(stderr,
			    "Couldn't get bin %"PRIrec"\n", bin->child[i]);
		    continue;
		}
		item = malloc(sizeof(bin_list));
		if (!item) goto unwind;
		item->rec = bin->child[i];
		item->parent_comp   = complement;
		item->parent_offset = offset;
		item->parent_size   = bin->size;
		item->parent_level  = level;
		item->child_index   = i;
		kid_offsets[nkids] = (complement
				      ? (offset+bin->size-kid->pos-kid->size)
				      : offset + kid->pos);
		kids[nkids++] = item;
	    }
	    if (nkids == 2
		&& (depth_first
		    ? kid_offsets[0] > kid_offsets[1]
		    : kid_offsets[0] < kid_offsets[1])) {
		/* Do rightmost child first */
		item = kids[0]; kids[0] = kids[1]; kids[1] = item;
	    }
	    for (i = 0; i < nkids; i++) {
		if (depth_first) {
		    /* Put item one down in the list,
		       as head will be removed later */
		    kids[i]->next = head->next;
		    head->next = kids[i];
		} else {
		    kids[i]->next = NULL;
		    tail->next = kids[i];
		    tail = kids[i];
		}
	    }
	} else {
	    fprintf(stderr, "Couldn't get bin %"PRIrec"\n", head->rec);
	}
	item = head;
	head = head->next;
	free(item);
    }

  return maxlevel;

 unwind:
  while (head) {
    item = head;
    head = head->next;
    free(item);
  }
  return -1;
}

/*
 * Produces a postscript file containing a plot of the contig bin structure.
 */
int contig_dump_ps(GapIO *io, contig_t **c, char *fn,
		   int depth_first, int draw_seqs) {
  FILE *gv = NULL;
  const char *preamble =
    "%!PS-Adobe-3.0\n"
    "%%Pages: 1\n"
    "%%BoundingBox: (atend)\n"
    "%%EndComments\n"
    "%%BeginProlog\n"
    "/f false def\n"
    "/t true def\n"
    "/fnt {\n"
    "    /Times-Roman findfont\n"
    "    exch 2 mul scalefont\n"
    "    setfont\n"
    "} bind def\n"
    "/itoa { 20 string cvs } bind def\n"
    "/sw { stringwidth pop } bind def\n"
    "/ruler {\n"
    "    % <start> <end> ruler\n"
    "    gsave\n"
    "    2 setlinewidth 1 setlinecap\n"
    "    exch dup dup sc mul 27 H add moveto\n"
    "    itoa dup sw neg 0 rmoveto show   % left number\n"
    "    sc mul 20 H add moveto\n"
    "    0 3 rmoveto 0 -6 rlineto 0 3 rmoveto          % left tick\n"
    "    dup sc mul 20 H add lineto                    % line\n"
    "    0 3 rmoveto 0 -6 rlineto 0 3 rmoveto          % right tick\n"
    "    dup sc mul 27 H add moveto itoa show % right number\n"
    "    stroke\n"
    "    grestore\n"
    "} bind def\n"
    "/strscale {\n"
    "    % <string> strscale\n"
    "    dup sw (999999) sw 2 copy gt { pop } { exch pop } ifelse\n"
    "    G 0.707 mul exch div dup scale\n"
    "} bind def\n"
    "/norm {\n"
    "    % <bin_size> <complemented> <pos>\n"
    "    exch\n"
    "    { neg add 1 sub }   % bin_size - 1 - pos\n"
    "    { exch pop }        % pos\n"
    "    ifelse\n"
    "} bind def\n"
    "/nmin {\n"
    "    % <bin_size> <complemented> <pos1> <pos2>\n"
    "    4 2 roll 2 copy 6 -1 roll norm  % norm_pos1\n"
    "    4 1 roll 3 -1 roll norm         % norm_pos1 norm_pos2\n"
    "    2 copy lt { pop } { exch pop } ifelse  % => min\n"
    "} bind def\n"
    "/nmax {\n"
    "    % <bin_size> <complemented> <pos1> <pos2>\n"
    "    4 2 roll 2 copy 6 -1 roll norm  % norm_pos1\n"
    "    4 1 roll 3 -1 roll norm         % norm_pos1 norm_pos2\n"
    "    2 copy gt { pop } { exch pop } ifelse  % => max\n"
    "} bind def\n"
    "/nminmax {\n"
    "    % <bin_size> <complemented> <pos1> <pos2>\n"
    "    4 2 roll 2 copy 6 -1 roll norm  % norm_pos1\n"
    "    4 1 roll 3 -1 roll norm         % norm_pos1 norm_pos2\n"
    "    2 copy gt { exch } if           % => min max\n"
    "} bind def\n"
    "/usedbox {\n"
    "    % <bin_size> <complemented> <used_start> <used_end> usedbox\n"
    "    gsave\n"
    "    0.8 setgray\n"
    "    nminmax 1 index sub sc mul H 6 sub\n"
    "    3 -1 roll sc mul 3 4 2 roll rectstroke\n"
    "    grestore\n"
    "} bind def\n"
    "/usedpos {\n"
    "    %<bin_size> <complemented> <used_start> <used_end> <offset> usedpos\n"
    "    gsave\n"
    "    4 index 6 2 roll nmin dup sc mul H 2 add moveto\n"
    "    dup 4 -1 roll add itoa dup sw 4 2 roll sub sc mul lt\n"
    "    { 0.8 0.8 scale show } { 45 rotate strscale show } ifelse\n"
    "    stroke\n"
    "    grestore\n"
    "} bind def\n"
    "/emptymark {\n"
    "    % <bin_size> emptymark\n"
    "    gsave\n"
    "    0 H 2 add moveto sc mul (empty) dup sw 3 -1 roll lt\n"
    "    { 0.8 0.8 scale show } { 45 rotate strscale show } ifelse\n"
    "    stroke\n"
    "    grestore\n"
    "} bind def\n"
    "/outerbox {\n"
    "    % <bin_size> <complemented> outerbox\n"
    "    gsave\n"
    "    exch sc mul exch\n"
    "    { 0.5 } { 0.0 } ifelse setgray\n"
    "    0 0 3 2 roll H rectstroke\n"
    "    grestore\n"
    "} bind def\n"
    "/binoffset {\n"
    "    % <offset> <bin_size> binoffset\n"
    "    gsave\n"
    "    sc mul exch itoa dup sw 3 -1 roll lt\n"
    "    { 0 -7 moveto 0.8 0.8 scale show }\n"
    "    { 0 -2 moveto 315 rotate strscale show } ifelse\n"
    "    stroke grestore\n"
    "} bind def\n"
    "/strcat {\n"
    "    % [(a) (b) (c)] strcat\n"
    "    0 1 index { length add } forall string\n"
    "    0 3 2 roll\n"
    "    { 3 copy putinterval length add } forall pop\n"
    "} bind def\n"
    "/recid {\n"
    "    % <nrefpos> <nseqs> <recno_str> <bin_size> recid\n"
    "    gsave\n"
    "    dup sc mul 5 1 roll\n"
    "    2 div sc mul H 2 div 2 sub moveto\n"
    "    (#) exch (/) 5 3 roll itoa (/) 3 2 roll itoa\n"
    "    6 array astore strcat\n"
    "    dup sw dup dup\n"
    "    4 index dup H 2 sub gt {\n"
    "        gt { 4 -1 roll exch div dup scale } { 4 -1 roll pop pop } ifelse\n"
    "    } {\n"
    "        6 -1 roll pop pop H 2 sub\n"
    "        gt { H 2 sub exch div dup scale } { pop } ifelse\n"
    "        90 rotate 2 -2.5 rmoveto\n"
    "    } ifelse\n"
    "    2 div neg 0 rmoveto show stroke\n"
    "    grestore\n"
    "} bind def\n"
    "/sqb {\n"
    "    % <bin_size> <complemented> <nseqs> sqb\n"
    "    gsave 1 setlinecap\n"
    "    H 8 sub exch 1 add div /vsc exch def\n"
    "    dup { 1 index 1 sub sc mul 0 translate -1 1 scale } if 1\n"
    "    % stack: <bin_size> <complemented> <vpos = 1>\n"
    "} bind def\n"
    "/sqg {\n"
    "    % <is_seq> sqg\n"
    "    stroke dup 1 exch sub exch 2 div dup setrgbcolor\n"
    "} bind def\n"
    "/sqe {\n"
    "    pop stroke grestore\n"
    "} bind def\n"
    "/sq {\n"
    "    % <vpos> <start> <end> sq\n"
    "    2 index vsc mul 4 add dup 4 1 roll\n"
    "    exch sc mul exch moveto\n"
    "    sc mul exch lineto\n"
    "    1 add\n"
    "} bind def\n"
    "/l {\n"
    "    % <root_pos> <root_level> <bin_pos> <bin_level> <child_index> l\n"
    "    gsave\n"
    "    0 eq { 1.0 0 0 setrgbcolor } { 0 0 1.0 setrgbcolor} ifelse\n"
    "    dup 3 index sub dup 1 gt {\n"
    "      2 div 5 1 roll\n"
    "	   exch sc mul exch\n"
    "      H G add mul H sub neg 2 copy moveto 4 2 roll\n"
    "      exch sc mul exch H G add mul neg 4 2 roll\n"
    "      4 index G mul add 4 2 roll 2 copy 6 index G mul sub\n"
    "      1 index 6 index lt\n"
    "      { 7 -1 roll G mul } { 7 -1 roll G neg mul } ifelse\n"
    "      dup 4 -1 roll add 3 1 roll\n"
    "      7 -1 roll add 6 1 roll 4 2 roll curveto\n"
    "    } {\n"
    "      pop\n"
    "      exch sc mul exch\n"
    "      H G add mul H sub neg moveto\n"
    "      exch sc mul exch\n"
    "      H G add mul neg lineto\n"
    "    } ifelse\n"
    "    stroke\n"
    "    grestore\n"
    "} bind def\n"
    "/sb {\n"
    "    % <nrefpos> <nseqs> <recno_str> <bin_size> <complemented>\n"
    "    %    <used_start> <used_end> <filled> <level> <offset>\n"
    "    gsave\n"
    "    dup 3 1 roll\n"
    "    sc mul exch\n"
    "    H G add mul neg translate\n"
    "    dup 6 index binoffset exch\n"
    "    { 5 copy pop usedbox 5 copy usedpos pop }\n"
    "    { pop 3 index emptymark } ifelse\n"
    "    pop pop 2 copy outerbox\n"
    "    % stack: <nrefpos> <nseqs> <recno_str> <bin_size> <complemented>\n"
    "} bind def\n"
    "/eb {\n"
    "    % <nrefpos> <nseqs> <recno_str> <bin_size> <complemented> eb\n"
    "    pop recid\n"
    "    grestore\n"
    "} bind def\n"
    "%%EndProlog\n"
    "%%Page: 1 1\n";
  
  const int h = 30;
  const int g = 20;
  int (*levels)[2];
  int nlevels = 16;
  int landscape;
  double pwidth;
  double pheight;
  tg_rec rec;
  bin_index_t *bin = NULL;
  int i;
  int maxlevel;
  double scy;

  levels = malloc(nlevels * sizeof(levels[0]));
  if (!levels) {
    perror("bin_dump_ps");
    return -1;
  }

  for (i = 0; i < nlevels; i++) {
    levels[i][0] =  INT_MAX;
    levels[i][1] = -INT_MAX;
  }

  gv = fopen(fn, "w+");
  if (!gv) { perror(fn); return -1; }
  fprintf(gv, "%s", preamble);
  cache_incr(io, *c);
  rec = contig_get_bin(c);
  bin = get_bin(io, rec);
  if (!bin) {
    fprintf(stderr, "Couldn't get bin %"PRIrec"\n", rec);
    cache_decr(io, *c);
    return -1;
  }
  cache_incr(io, bin);

  maxlevel = count_levels_recurse(io, rec, &levels, &nlevels, depth_first);

  if (maxlevel < 0) {
    cache_decr(io, bin);
    fclose(gv);
    return -1;
  }

  landscape = (h + g) * (maxlevel + 1) < 500;
  pwidth  = landscape ? 800.0 : 550.0;
  pheight = landscape ? 550.0 : 800.0;
  scy = ((h + g) * (maxlevel + 1) > pheight
	 ? pheight / ((h + g) * (maxlevel + 1)) : 1.0);
  
  fprintf(gv, 
	  "%%%%PageOrientation: %s\n"
	  "%%%%BeginPageSetup\n"
	  "5 fnt\n"
	  "/H %d def\n"
	  "/G %d def\n"
	  "/sc %f def\n"
	  "%d rotate\n"
	  "%d %d translate\n"
	  "%f %f scale\n"
	  "%d sc mul 0 translate\n"
	  "%%%%EndPageSetup\n"
	  "newpath\n"
	  "%d %d ruler\n",
	  landscape ? "Landscape" : "Portrait",
	  h, g, pwidth / bin->size / scy, 
	  landscape ? 90 : 0,
	  20, landscape ? -h - g - 30 : (int) pheight + 20,
	  scy, scy,
	  -contig_offset(io, c), (*c)->start, (*c)->end);

  bin_dump_ps_recurse(gv, io, rec, levels, nlevels,
		      depth_first, DRAW_LINKS);
  bin_dump_ps_recurse(gv, io, rec, levels, nlevels, depth_first,
		      DRAW_BINS | (draw_seqs ? DRAW_SEQS : 0));
  cache_decr(io, bin);
  cache_decr(io, *c);
  if (landscape) {
    fprintf(gv, "showpage\n"
	    "%%%%Trailer\n"
	    "%%%%BoundingBox: 0 0 %.0f %.0f\n"
	    "%%%%EOF\n",
	    pheight + 40, pwidth + 40);
  } else {
    fprintf(gv, "showpage\n"
	    "%%%%Trailer\n"
	    "%%%%BoundingBox: 0 0 %.0f %.0f\n"
	    "%%%%EOF\n",
	    pwidth + 40, pheight + 40);
  }
  if (fclose(gv)) { perror(fn); return -1; }
  return 0;
}

static int bin_dump_graphviz_recurse(FILE *out, GapIO *io,
				     tg_rec bin_num,  bin_index_t *bin,
				     int offset, int complement) {
  int res = 0;
  int i;

  if (bin->flags & BIN_COMPLEMENTED) {
    complement ^= 1;
  }
  if (bin_empty(bin)) {
    fprintf(out,
	    "%"PRIrec" [label=\"\\N\\l"
	    "offset=%d\\lsize=%d\\ldirn=%s\\lused=empty\\l"
	    "nrefpos=%d\\lnseqs=%d\\l\" shape=hexagon fontsize=10 "
	    "fontname=\"Courier\" style=filled fillcolor=\"%s\"]\n",
	    bin->rec, offset, bin->size, complement ? "rev" : "fwd",
	    bin->nrefpos, bin->nseqs, complement ? "#98fb98" : "#87ceeb");
  } else {
    fprintf(out,
	    "%"PRIrec" [label=\"\\N\\l"
	    "offset=%d\\lsize=%d\\ldirn=%s\\lused=%d..%d\\l"
	    "nrefpos=%d\\lnseqs=%d\\l\" shape=box fontsize=10 "
	    "fontname=\"Courier\" style=filled fillcolor=\"%s\"];\n",
	    bin->rec, offset, bin->size, complement ? "rev" : "fwd",
	    bin->start_used, bin->end_used,
	    bin->nrefpos, bin->nseqs, complement ? "#98fb98" : "#87ceeb");
  }
  for (i = 0; i < 2; i++) {
    bin_index_t *ch;
    int child_offset;
    
    if (!bin->child[i]) continue;
    ch = get_bin(io, bin->child[i]);
    if (!ch) {
      fprintf(stderr, "Couldn't get bin %"PRIrec"\n", bin->child[i]);
      continue;
    }
    
    cache_incr(io, ch);
    child_offset = (complement ? offset + bin->size - ch->pos - ch->size
		    : offset + ch->pos);
    res |= bin_dump_graphviz_recurse(out, io, bin->child[i], ch, child_offset,
				     complement);
    fprintf(out, "%"PRIrec" -> %"PRIrec" [dir=forward arrowhead=normal];\n",
	    bin->rec, bin->child[i]);
    cache_decr(io, ch);
  }

  return res;
}

/*
 * Produces a dump of the bin structure in graphviz format, can be used
 * with dot, dotty etc.
 */

int contig_dump_graph(GapIO *io, contig_t **c, char *fn) {
  FILE *out = NULL;
  bin_index_t *bin = NULL;
  tg_rec rec;

  out = fopen(fn, "w+");
  if (!out) { perror(fn); return -1; }
  cache_incr(io, *c);
  rec = contig_get_bin(c);
  bin = get_bin(io, rec);
  if (!bin) {
    fprintf(stderr, "Couldn't get bin %"PRIrec"\n", rec);
    cache_decr(io, *c);
    fclose(out);
    return -1;
  }
  fprintf(out,
	  "digraph Contig%"PRIrec" {\n"
	  "root=\"%"PRIrec"\";\n",
	  (*c)->rec, rec);
  cache_incr(io, bin);
  bin_dump_graphviz_recurse(out, io, rec, bin, contig_offset(io, c), 0);
  cache_decr(io, bin);
  fprintf(out, "}\n");
  cache_decr(io, *c);
  if (fclose(out)) {
    perror(fn);
    return -1;
  }

  return 0;
}


/* -------------------------------------------------------------------------
 * Padded / reference coordinate mappings.
 */

/*
 * Convert a padded coordinate to a reference coordinate.
 * *dir_p is filled out, if non-null, to hold the direction we count in
 * when going right. -1 => no data (no ref), 0 => upwards, 1 => downwards.
 * ref_id, if non-null, will be returned containing the reference ID used
 * for computing the reference position; -1 if not known
 */
int padded_to_reference_pos(GapIO *io, tg_rec cnum, int ppos, int *dir_p,
			    int *ref_id) {
    contig_iterator *ci;
    rangec_t *r;
    int rpos;
    int dir;

    ci = contig_iter_new_by_type(io, cnum, 1, CITER_FIRST | CITER_ISTART, 
				 ppos, CITER_CEND, GRANGE_FLAG_ISREFPOS);

    if (!ci) {
	if (ref_id) *ref_id = -1;
	if (dir_p)  *dir_p  = -1;
	return ppos;
    }

    r = contig_iter_next(io, ci);

    /*
     * Maybe we're at the end of the contig. If so we search backwards instead
     * of forwards to find the previous REFPOS marker.
     */
    if (!r) {
	contig_iter_del(ci);
	ci = contig_iter_new_by_type(io, cnum, 1, CITER_LAST | CITER_ISTART,
				     CITER_CSTART, ppos, GRANGE_FLAG_ISREFPOS);
	if (!ci) {
	    if (ref_id) *ref_id = -1;
	    if (dir_p)  *dir_p  = -1;
	    return ppos;
	}

	r = contig_iter_prev(io, ci);
	if (!r) {
	    contig_iter_del(ci);
	    if (dir_p)  *dir_p = -1;
	    if (ref_id) *ref_id = -1;
	    return ppos;
	}

	dir = 0 ^ r->comp;
    } else {
	dir = 1 ^ r->comp;
    }

    if (((r->flags & GRANGE_FLAG_REFPOS_DIR) == GRANGE_FLAG_REFPOS_FWD) ^ r->comp)
	rpos = r->mqual + (ppos - r->start + dir);
    else
	rpos = r->mqual - (ppos - r->start - dir);

    if ((r->flags & GRANGE_FLAG_REFPOS_INDEL) == GRANGE_FLAG_REFPOS_DEL) {
	if (dir == 1) /* fwd */
	    rpos -= (ppos < r->start) * r->pair_rec + 1;
    }
    
    if (dir_p)
	*dir_p = r->comp;
    if (ref_id)
	*ref_id = r->rec;

    contig_iter_del(ci);

    return rpos;
}

/*
 * Looks for a refpos marker at position ppos. If it finds one it returns
 * 0 and sets the bin and bin_idx fields.
 *
 * If it doesn't, it returns -1;
 */
int find_refpos_marker(GapIO *io, tg_rec cnum, int ppos,
		       tg_rec *bin_r, int *bin_idx_r, rangec_t *rp) {
    contig_iterator *ci;
    rangec_t *r;

    if (bin_r)
	*bin_r = 0;
    if (bin_idx_r)
	*bin_idx_r = 0;

    ci = contig_iter_new_by_type(io, cnum, 0, CITER_FIRST | CITER_ISTART, 
				 ppos, ppos, GRANGE_FLAG_ISREFPOS);
    if (!ci)
	return -1;

    r = contig_iter_next(io, ci);

    if (!r) {
	contig_iter_del(ci);
	return -1;
    }

    if (r->start == ppos && r->end == ppos) {
	if (bin_r)
	    *bin_r = r->orig_rec;
	if (bin_idx_r)
	    *bin_idx_r = r->orig_ind;
	*rp = *r;
	contig_iter_del(ci);
	return 0;
    }

    contig_iter_del(ci);
    return -1;
}


/*
 * Converts a range of padded coordinates to reference coordinates.
 * Optionally also fills out a reference seq ID array too, if non-NULL.
 *
 * ref_pos and ref_id should be allocated by the caller to be of
 * appropriate size (paddeed_end - padded_start + 1).
 *
 * Returns 0 on success
 *        -1 on failure
 */
int padded_to_reference_array(GapIO *io, tg_rec cnum,
			      int padded_start, int padded_end,
			      int *ref_pos, int *ref_id) {
    int dir, rpos, i, rid;
    contig_iterator *ci;
    int len = padded_end - padded_start + 1;
    rangec_t *r;

    /* Starting point */
    rpos = padded_to_reference_pos(io, cnum, padded_start, &dir, &rid);
    switch (dir) {
    case -1: dir = +1; break; /* guess */
    case  0: dir = +1; break;
    case  1: dir = -1; break;
    }

    /* Now iterate from here up to padded_end */
    ci = contig_iter_new_by_type(io, cnum, 0, CITER_FIRST | CITER_ISTART, 
				 padded_start, padded_end,
				 GRANGE_FLAG_ISREFPOS);
    if (!ci) {
	for (i = 0; i < len; i++) {
	    ref_pos[i] = rpos;
	    rpos += dir;
	    if (ref_id)
		ref_id[i] = rid;
	}
	return 0;
    }

    i = 0;
    while ((r = contig_iter_next(io, ci))) {
	int fwd;

	while (i < len && padded_start < r->start) {
	    ref_pos[i] = rpos;
	    rpos += dir;
	    rid = r->rec;
	    if (ref_id)
		ref_id[i] = rid;
	    i++;
	    padded_start++;
	}

	dir = 1 ^ r->comp;

	fwd = ((r->flags & GRANGE_FLAG_REFPOS_DIR) == GRANGE_FLAG_REFPOS_FWD);
	if (fwd ^ r->comp)
	    rpos = r->mqual + (padded_start - r->start + dir);
	else
	    rpos = r->mqual - (padded_start - r->start - dir);

	if ((r->flags & GRANGE_FLAG_REFPOS_INDEL) == GRANGE_FLAG_REFPOS_DEL) {
	    if (dir == 1) /* fwd */
		rpos -= (padded_start < r->start) * r->pair_rec + 1;
	} else {
	    ref_pos[i] = rpos;
	    if (ref_id)
		ref_id[i] = -1; /* indel */
	    i++;
	    padded_start++;
	}
    }

    while (i < len) {
	ref_pos[i] = rpos;
	rpos += dir;
	if (ref_id)
	    ref_id[i] = rid;
	i++;
    }

    return 0;
}


/*
 * Given a contig record and a reference position, attempt to return
 * the padded coordinate. Note this may not exist, it may in extreme cases
 * exist multiple times (after breaking and rejoining), or it may exist
 * only once but be too hard to discover. If we don't care which specific
 * reference ID to search, pass in ref_id == -1.
 *
 * We only attempt to tackle easy cases of a single reference in this contig
 * with monotonically increasing or decreasing values. (We need more data
 * stored to allow arbitrary queries to be fast.)
 *
 * Returns 0 on success, position stored in *padded_pos.
 *        -1 on failure, *padded_pos undefined.
 */
int reference_to_padded_pos(GapIO *io, tg_rec cnum, int ref_id, int ref_pos,
			    int *padded_pos) {
    int dir_s, dir_e, dir, rid;
    contig_t *c = cache_search(io, GT_Contig, cnum);
    int cstart = c->start;
    int cend = c->end;
    int rstart, rend;
    int cguess, rguess;

    /* FIXME: ignoring ref_id for now */

    rstart = padded_to_reference_pos(io, cnum, cstart, &dir_s, &rid);
    if (ref_id != -1 && rid != ref_id)
	return -1;

    rend   = padded_to_reference_pos(io, cnum, cend,   &dir_e, &rid);
    if (ref_id != -1 && rid != ref_id)
	return -1;

    if (dir_s != dir_e)
	return -1;

    if (ref_pos == rstart) {
	*padded_pos = cstart;
	return 0;
    } else if (ref_pos == rend) {
	*padded_pos = cend;
	return 0;
    }


    while ((dir_s == 0 && ref_pos >= rstart && ref_pos <= rend) ||
	   (dir_s == 1 && ref_pos <= rstart && ref_pos >= rend)) {

	/* Make a guess */
	cguess = (ref_pos - rstart) / (rend - rstart + .0) * (cend - cstart)
	    + cstart;

	/*
	 * Same as previous guess, but still not found it. Maybe it doesn't
	 * exist (eg deletion), so return the guess anyway
	 */
	if (cguess == cstart || cguess == cend) {
	    *padded_pos = cguess;
	    return 0;
	}

	rguess = padded_to_reference_pos(io, cnum, cguess, &dir, &rid);
	if (ref_id != -1 && rid != ref_id)
	    return -1;

	//printf("(%d,%d)..(%d,%d) => guess (%d,%d)\n",
	//       cstart, rstart,  cend, rend,  cguess, rguess);

	if (rguess == ref_pos) {
	    *padded_pos = cguess;
	    return 0;
	}

	if (rguess < ref_pos) {
	    cstart = cguess;
	    rstart = rguess;
	} else {
	    cend = cguess;
	    rend = rguess;
	}
    }

    return -1;
}

/*
 * As above, but starting from a single known point on that reference.
 * This allows for reference positions to occur more than once.
 */
int reference_to_padded_pos2(GapIO *io, tg_rec cnum, int ref_id, int ref_pos,
			     int cur_padded_pos, int *padded_pos) {
    int dir, rid;
    int cguess, rguess, cpos, rpos, try;
    int last1 = INT_MAX, last2 = INT_MAX;

    cpos = cur_padded_pos;
    rpos = padded_to_reference_pos(io, cnum, cpos, &dir, &rid);

    printf("\nLooking for %d\n", ref_pos);
    printf("Starting at %d,%d\n", cpos, rpos);

    if (ref_id != -1 && rid != ref_id)
	return -1;

    for (try = 0; try < 100; try++) {
	if (dir == 0 /* fwd */ || dir == -1) {
	    cguess = cpos + ref_pos - rpos;
	} else {
	    cguess = cpos + rpos - ref_pos;
	}
	
	rguess = padded_to_reference_pos(io, cnum, cguess, &dir, &rid);

	if (ref_id != -1 && rid != ref_id)
	    return -1;

	if (rguess == ref_pos) {
	    *padded_pos = cguess;
	    return 0;
	}

	if (cguess == last2) {
	    printf("Loop detected - guessing\n");
	    /* 2-step loop, likely due to deletion. Just use our best guess */
	    *padded_pos = (last1 + last2) / 2;
	    return 0;
	}

	cpos = cguess;
	rpos = rguess;

	last2 = last1;
	last1 = cpos;
    }

    return -1;
}


/* 
 * Moves an entire contig by a relative amount left (-ve) or right (+ve).
 * Returns 0 on success
 *        -1 on failure
 */
int move_contig(GapIO *io, tg_rec crec, int distance) {
    contig_t *c;
    bin_index_t *bin;
 
    if (!(c = cache_search(io, GT_Contig, crec)))
	return -1;
    if (!(c = cache_rw(io, c)))
	return -1;

    if (!(bin = cache_search(io, GT_Bin, contig_get_bin(&c))))
	return -1;
    if (!(bin = cache_rw(io, bin)))
	return -1;

    bin->pos += distance;
    bin->flags |= BIN_BIN_UPDATED;
    c->start += distance;
    c->end   += distance;

    c->timestamp = io_timestamp_incr(io);

    return 0;
}

/*
 * Sets the visible start of a contig, both in the contig structure and
 * also the root bin to keep everything internally consistent.
 *
 * Ie it's a more end-user version of contig_set_start.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_set_visible_start(GapIO *io, tg_rec contig, int pos) {
    int cstart;

    if (-1 == consensus_valid_range(io, contig, &cstart, NULL))
	return -1;

    return move_contig(io, contig, pos - cstart);
}


/*
 * Copies the nseq, nanno and nrefpos from the root bin to the contig
 * struct.
 *
 * This can be necessary during algorithms that move data around, for example
 * break_contig().
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_fix_nseq(GapIO *io, contig_t *c) {
    if (!c)
	return -1;

    if (c->bin) {
	bin_index_t *bin = cache_search(io, GT_Bin, c->bin);
	if (!bin)
	    return -1;

	c->nseqs = bin->nseqs;
	c->nanno = bin->nanno;
	c->nrefpos = bin->nrefpos;
    } else {
	c->nseqs = 0;
	c->nanno = 0;
	c->nrefpos = 0;
    }
    
    return 0;
}

/*
 * Adds a bi-directional contig link.
 * The data to link from/to is in the passed in abs_link. Coordinates here
 * are all contig absolute coords.
 *
 * Internally this link gets created at both ends (rec1, rec2) and the
 * positions are converted into relative coordinates. (The end1 and end2
 * fields in abs_link are ignored and computed as required.)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int contig_add_link(GapIO *io, contig_link_t *abs_link) {
    contig_t *c1, *c2;
    contig_link_t *l1, *l2;
    int idx1, idx2;
    int pos1 = abs_link->pos1;
    int pos2 = abs_link->pos2;

    /* Create contig link structs */
    if (!(c1 = cache_search(io, GT_Contig, abs_link->rec1)))
	return -1;
    if (!(c1 = cache_rw(io, c1)))
	return -1;

    if (!c1->link)
	c1->link = ArrayCreate(sizeof(contig_link_t), 0);

    ArrayRef(c1->link, idx1 = ArrayMax(c1->link));

    if (!(c2 = cache_search(io, GT_Contig, abs_link->rec2)))
	return -1;
    if (!(c2 = cache_rw(io, c2)))
	return -1;

    if (!c2->link)
	c2->link = ArrayCreate(sizeof(contig_link_t), 0);

    ArrayRef(c2->link, idx2 = ArrayMax(c2->link));


    /* Fill out the links */
    l1 = arrp(contig_link_t, c1->link, idx1);
    l2 = arrp(contig_link_t, c2->link, idx2);
    *l1 = *abs_link;
    *l2 = *abs_link;
    l1->rec2 = abs_link->rec2;
    l2->rec2 = abs_link->rec1;

    printf("Linking %"PRIrec" to %"PRIrec"\n", c1->rec, c2->rec);
    printf("&c1->link = %p, &c2->link = %p\n", &c1->link, &c2->link);
    printf("c1->link = %p, c2->link = %p\n", c1->link, c2->link);

    if (pos1 - c1->start < c1->end - pos1) {
	l1->end1 = 0;
	l1->pos1 = pos1 - c1->start;
	l2->end2 = 0;
	l2->pos2 = pos1 - c1->start;
    } else {
	l1->end1 = 1;
	l1->pos1 = c1->end - pos1;
	l2->end2 = 1;
	l2->pos2 = c1->end - pos1;
    }

    if (pos2 - c2->start < c2->end - pos2) {
	l2->end1 = 0;
	l2->pos1 = pos2 - c2->start;
	l1->end2 = 0;
	l1->pos2 = pos2 - c2->start;
    } else {
	l2->end1 = 1;
	l2->pos1 = c2->end - pos2;
	l1->end2 = 1;
	l1->pos2 = c2->end - pos2;
    }

    return 0;
}

/*
 * Converts a specific link number from relative coordinates to absolute
 * contig coordinates.
 *
 * Input is rel_link, output is abs_link.
 * Returns 0 on success
 *        -1 on failure
 */
int contig_get_link_positions(GapIO *io,
			      contig_link_t *rel_link,
			      contig_link_t *abs_link) {
    contig_t *c1, *c2;

    c1 = cache_search(io, GT_Contig, rel_link->rec1);
    c2 = cache_search(io, GT_Contig, rel_link->rec2);
    
    if (!c1 || !c2)
	return -1;

    *abs_link = *rel_link;

    abs_link->pos1 = abs_link->end1 == 0
	? c1->start + abs_link->pos1
	: c1->end   - abs_link->pos1;

    abs_link->pos2 = abs_link->end2 == 0
	? c2->start + abs_link->pos2
	: c2->end   - abs_link->pos2;

    return 0;
}
