/*
 * TODO:
 * - proper deallocation system (contigs, bins, seqs, annos)
 * - tags on consensus need to move if we're breaking a contig, dup otherwise?
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "dis_readings.h"
#include "misc.h"
#include "break_contig.h"

/*-----------------------------------------------------------------------------
 * Databse checking code.
 */

/* Recursively check linkage on bins and nseqs relationships.
 * brec  = bin record
 * prec  = parent record, to validate against.
 * ptype = type of parent record.
 *
 * Returns nseqs on success
 *         -1 on failure.
 */
//static FILE *errfp = stderr;
static FILE *errfp = NULL;;

static int check_contig_bins_r(GapIO *io, tg_rec brec, int ptype, tg_rec prec){
    bin_index_t *b;
    tg_rec copy_c0, copy_c1;
    int copy_nseq, i, ns;

    /* Check prec/ptype */
    b = cache_search(io, GT_Bin, brec);

    if (b->parent != prec || b->parent_type != ptype) {
	fprintf(errfp, "ERROR: bin parent record/type mismatch for bin %"
		PRIrec" : parent = %"PRIrec"/%"PRIrec" type = %d/%d\n",
		brec, b->parent, prec, b->parent_type, ptype);
	abort();
	return -1;
    }

    /* Count number of sequences in this contig */
    for (ns = i = 0; b->rng && i < ArrayMax(b->rng); i++) {
	range_t *r = arrp(range_t, b->rng, i);
	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if ((r->flags & GRANGE_FLAG_ISMASK) == GRANGE_FLAG_ISSEQ)
	    ns++;
    }

    /* Take copies of the data from b so we don't have to lock it and
     * can allow it to get purged from the hache.
     */
    copy_c0   = b->child[0];
    copy_c1   = b->child[1];
    copy_nseq = b->nseqs;

    /* Recurse and check nseqs */
    if (copy_c0) {
	if ((i = check_contig_bins_r(io, copy_c0, GT_Bin, brec)) == -1)
	    return -1;
	ns += i;
    }
    if (copy_c1) {
	if ((i = check_contig_bins_r(io, copy_c1, GT_Bin, brec)) == -1)
	    return -1;
	ns += i;
    }

    if (ns != copy_nseq) {
	fprintf(errfp, "ERROR: nseq mismatch for bin %"PRIrec" : %d/%d\n",
		brec, ns, copy_nseq);
	abort();
	return -1;
    }

    return ns;
}

int check_contig_bins(GapIO *io) {
    int i, ret = 0;

    errfp = stdout;

    printf("check_contig_bins start, ncontigs=%d\n", io->db->Ncontigs);
    if (io->db->Ncontigs <= 340)
	return 0;

    for (i = 0; i < io->db->Ncontigs; i++) {
	tg_rec crec = arr(tg_rec, io->contig_order, i);
	contig_t *c = cache_search(io, GT_Contig, crec);
	//printf("Check contig %d root %d\n", crec, c->bin);
	if (c->bin) {
	    if (check_contig_bins_r(io, c->bin, GT_Contig, crec) == -1)
		ret = -1;
	}
    }

    printf("check_contig_bins end => %d\n", ret);

    return ret;
}

/* Singular of above */
int check_contig_bin(GapIO *io, tg_rec crec) {
    contig_t *c = cache_search(io, GT_Contig, crec);
    errfp = stdout;

    printf("Check contig %"PRIrec" root %"PRIrec"\n", crec, c->bin);
    if (c->bin) {
	if (check_contig_bins_r(io, c->bin, GT_Contig, crec) == -1)
	    return -1;
    }

    return 0;
}



/*-----------------------------------------------------------------------------
 * Internal functions & data types
 */

typedef struct {
    tg_rec contig;
    int start;
    int end;
    tg_rec rec;       /* sequence record */
    range_t rng;      /* rng in original bin */
    rangec_t *anno;   /* Annotation ranges */
    int n_anno;
} r_pos_t;


/*
 * Part 1/2 of the disassemble_readings implementation (see below).
 *
 * 1. Produce a table of which contig number and region each reading
 *    record is within. Also track the anno_eles for this seq too.
 *
 * 2. Remove these readings from their contigs.
 *
 * 3. If "remove" is true, ensure if it's a read-pair that the other
 *    end gets the pairing information and flags cleared.
 *
 * Writes back to pos.
 *
 * The 'remove' argument is a boolean indicating whether the sequence
 * is about to be removed from the database. (Ie move==0 in
 * disassemble_reads.) We use this to determine whether to update the
 * read-pairing statuses.
 *
 * FIXME: we need to do something with annotations too. Unlink them, but also
 * keep track of their original locations relative to the sequence location
 * so we can add them back if required (move==1 or move==2).
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int unlink_read(GapIO *io, tg_rec rec, r_pos_t *pos, int remove) {
    contig_t *c;
    bin_index_t *bin;
    seq_t *seq;
    tg_rec brec;
    int i, j;

    //printf("%soving record #%d\n", remove ? "Rem" : "M", rec);

    /* Get location */
    if (bin_get_item_position(io, GT_Seq, rec,
			      &pos->contig,
			      &pos->start,
			      &pos->end,
			      NULL, /* orient */
			      &brec,
			      &pos->rng,
			      (void **)&seq)) {
	return -1;
    }
    /* seq already has cache_incr on it */

    /* Find annotations in this region */
    c = cache_search(io, GT_Contig, pos->contig);
    cache_incr(io, c);
    pos->anno = contig_anno_in_range(io, &c, pos->start, pos->end, 0,
				     &pos->n_anno);

    /* Filter for only those from this read */
    for (i = j = 0; i < pos->n_anno; i++) {
	if (pos->anno[i].pair_rec == rec)
	    pos->anno[j++] = pos->anno[i];
    }
    pos->n_anno = j;

    /*
     * Remove from bin range array
     *
     * FIXME: could be done more optimal if we knew bin->rng[]
     * index. Maybe update bin_get_item_position to return this
     * value so we can directly manipulate r->flags.
     *
     * FIXME2: if we're removing lots of items from the same bin,
     * then we could aggregate the nseq adjustment and call
     * bin_incr_nseq with something other than -1 (ie -nseq_in_bin).
     */
    bin = cache_search(io, GT_Bin, brec);
    if (!c || !bin) {
	cache_decr(io, seq);
	cache_decr(io, c);
	return -1;
    }

    if (bin_remove_item_from_bin(io, &c, &bin, GT_Seq, rec)) {
	cache_decr(io, seq);
	cache_decr(io, c);
	return -1;
    }
    cache_decr(io, seq);
    cache_decr(io, c);

    /* For read-pairs, unlink with rest of template */
    if (remove && pos->rng.pair_rec) {
	range_t *r;

	seq = cache_search(io, GT_Seq, pos->rng.pair_rec);
	if (seq->parent_rec == pos->rng.pair_rec) {
	    seq = cache_rw(io, seq);
	    seq->parent_type = 0;
	    seq->parent_rec = 0;
	}

	/* Pair is held in bin range too */
	if (bin_get_item_position(io, GT_Seq, pos->rng.pair_rec,
				  NULL, NULL, NULL, NULL,
				  &brec, NULL, NULL))
	    return -1;
	bin = cache_search(io, GT_Bin, brec);
	r = arrp(range_t, bin->rng, seq->bin_index);
	assert(r->rec == seq->rec);
	bin = cache_rw(io, bin);

	/* Fix other end's range_t */
	r->pair_rec = 0;
	r->flags &= ~GRANGE_FLAG_TYPE_MASK;
	r->flags |=  GRANGE_FLAG_TYPE_SINGLE;
    }

    return 0;
}

void bin_destroy_recurse(GapIO *io, tg_rec rec) {
    bin_index_t *bin = cache_search(io, GT_Bin,rec);

    cache_incr(io, bin);
    if (bin->child[0]) bin_destroy_recurse(io, bin->child[0]);
    if (bin->child[1]) bin_destroy_recurse(io, bin->child[1]);
    cache_decr(io, bin);

    cache_rec_deallocate(io, GT_Bin, rec);
}

/*
 * Checks if a bin is truely empty. It's not just sufficient to check
 * nseqs==0 as we could have tags or refpos markers too.
 * This differs to bin_empty() in that we are working out whether there is
 * no data in this bin and the child bins. Whereas the former is simply
 * interested in the bin->rng contents, so data within that specific bin.
 *
 * Returns 1 if empty.
 *         0 if not.
 */
static int bin_plus_children_empty(bin_index_t *bin) {
    int i;

    if (bin->nseqs || bin->nrefpos || bin->nanno)
	return 0;

    if (!bin->rng)
	return 1;

    for (i = 0; i < ArrayMax(bin->rng); i++) {
	range_t *r = arrp(range_t, bin->rng, i);
	if (!(r->flags & GRANGE_FLAG_UNUSED))
	    return 0;
    }

    return 1;
}

/*
 * Looks for contig gaps between start..end in contig and if it finds them,
 * breaking the contig in two.
 */
int remove_contig_holes(GapIO *io, tg_rec contig, int start, int end,
			int empty_contigs_only) {
    contig_t *c;
    bin_index_t *bin;
    contig_iterator *iter;
    rangec_t *r;
    int last;
    int contig_start, contig_end;

    /* Destroy contigs if they're now entirely empty */
    c = cache_search(io, GT_Contig, contig);
    cache_incr(io, c);

    bin = cache_search(io, GT_Bin, c->bin);
    if (bin_plus_children_empty(bin)) {
	puts("Removing empty contig");

	if (c->bin)
	    bin_destroy_recurse(io, c->bin);
	cache_decr(io, c);
	contig_destroy(io, contig);
	return 0;
    }

    /* Invalidate any cached consensus copies */
    if (bin_invalidate_consensus(io, contig, start, end) != 0) {
	cache_decr(io, c);
	return -1;
    }

    if (empty_contigs_only) {
	cache_decr(io, c);
	return 0;
    }


    /* Hole at left end */
    if (c->start == start) {
	iter = contig_iter_new(io, contig, 1, CITER_FIRST, start, end);
	if (iter) {
	    r = contig_iter_next(io, iter);
	    if (r) {
		c = cache_rw(io, c);
		start = c->start = r->start;
	    }
	    contig_iter_del(iter);
	}
    }

    /* Hole at right end */
    if (c->end == end) {
	iter = contig_iter_new(io, contig, 1, CITER_LAST | CITER_IEND,
			       start, end);
	if (iter) {
	    r = contig_iter_prev(io, iter);
	    if (r) {
		c = cache_rw(io, c);
		end = c->end = r->end;
	    }
	    contig_iter_del(iter);
	}
    }

    /* Make sure start/end are within clipped contig coordinates */
    consensus_valid_range(io, contig, &contig_start, &contig_end);
    if (start < contig_start)
	start = contig_start;
    if (end > contig_end)
	end = contig_end;


    /*
     * Look for holes in the middle, using ICLIPPEDEND mode so the data
     * is sorted by clipped sequence coords rather than just the r->end
     * rangec_t elements we use with most iterators.
     */
    iter = contig_iter_new(io, contig, 0, CITER_LAST | CITER_ICLIPPEDEND,
			   start, end);
    if (!iter) {
	cache_decr(io, c);
	return 0;
    }

    last = end;
    while (r = contig_iter_prev(io, iter)) {
	seq_t *s = cache_search(io, GT_Seq, r->rec);
	int cstart, cend;

	if (!s) {
	    cache_decr(io, c);
	    return -1;
	}

	if ((s->len < 0) ^ r->comp) {
	    cstart = r->start + ABS(s->len) - (s->right-1) - 1;
	    cend   = r->start + ABS(s->len) - (s->left-1) - 1;
	} else {
	    cstart = r->start + s->left-1;
	    cend   = r->start + s->right-1;
	}

	//printf("Seq %d, from %d..%d clipped %d..%d\n",
	//       r->rec, r->start, r->end, cstart, cend);

	if (cend < last) {
	    vmessage("GAP from %d..%d; breaking.\n", cend, last);
	    break_contig(io, contig, last, 0);

	    /* Who knows what impact break_contig has - restart to be safe */
	    contig_iter_del(iter);
	    iter = contig_iter_new(io, contig, 0,
				   CITER_LAST | CITER_ICLIPPEDEND,
	    			   start, last);
	}
	if (last > cstart)
	    last = cstart;
    }
    contig_iter_del(iter);

    cache_decr(io, c);

    return 0;
}

static GapIO *xio = NULL;


/* qsort callback */
static int pos_sort(const void *vp1, const void *vp2) {
    const r_pos_t *p1 = (const r_pos_t *)vp1;
    const r_pos_t *p2 = (const r_pos_t *)vp2;

#if 0
    /* For stable sorting regardless of rec deallocation, use this: */
    if (p1->contig != p2->contig) {
	contig_t *c1 = cache_search(xio, GT_Contig, p1->contig);
	contig_t *c2 = cache_search(xio, GT_Contig, p2->contig);
	return strcmp(c1->name, c2->name);
    }
#endif

    if (p1->contig != p2->contig)
	return p1->contig - p2->contig;

    return p1->start - p2->start;
}

static int pos_sort_end(const void *vp1, const void *vp2) {
    const r_pos_t *p1 = (const r_pos_t *)vp1;
    const r_pos_t *p2 = (const r_pos_t *)vp2;

#if 0
    if (p1->contig != p2->contig) {
	contig_t *c1 = cache_search(xio, GT_Contig, p1->contig);
	contig_t *c2 = cache_search(xio, GT_Contig, p2->contig);
	return strcmp(c1->name, c2->name);
    }
#endif

    if (p1->contig != p2->contig)
	return p1->contig - p2->contig;

    return p1->end - p2->end;
}


/*
 * Part 4 of the disassemble_readings implementation.
 * 
 * 4. Check for holes within the "source" contigs. We use the
 *    initial table for this, along with regions to look around.
 *
 * If empty_only is true we only check for entirely empty contigs (and
 * deallocate them), otherwise we also break contig where holes exist.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int fix_holes(GapIO *io, r_pos_t *pos, int npos,
		     int remove_holes) {
    int i, start, end;
    tg_rec contig;

    if (npos == 0)
	return 0;

    /*
     * pos[] is sorted by start coordinate. We need to loop backwards
     * through it though as break_contig will be producing a new right
     * hand contig, so moving backwards means we're always marching
     * down the un-changing (in contig number terms) left-hand contig.
     *
     * Due to sorted start coord but ragged end coord, we sort by
     * end instead so the overlapping segment detector works.
     */
    qsort(pos, npos, sizeof(*pos), pos_sort_end);

    /*
     * Step through pos finding overlapping reads so we can get entire
     * spans where we've removed data. We use this to limit our hole
     * fixing to just that region, reducing the search time on huge
     * contigs.
     */
    contig = pos[npos-1].contig;
    start  = pos[npos-1].start;
    end    = pos[npos-1].end;
    for (i = npos-1; i >= 0; i--) {
	if (pos[i].contig != contig ||
	    pos[i].end < start) {
	    remove_contig_holes(io, contig, start, end, !remove_holes);
	    contig = pos[i].contig;
	    start  = pos[i].start;
	    end    = pos[i].end;
	} else {
	    if (start > pos[i].start)
		start = pos[i].start;
	}
    }
    remove_contig_holes(io, contig, start, end, !remove_holes);

    return 0;
}


/*
 * Removes a sequence from the database, removing the bin range entry
 * and taking out of the BTree name index too if appropriate.
 */
static int seq_deallocate(GapIO *io, r_pos_t *pos) {
    int i;
    seq_t *s;
    bin_index_t *b;
    range_t *r;
    contig_t *c;

    /* Remove from sequence name btree index */
    if (!(s = cache_search(io, GT_Seq, pos->rec)))
	return -1;
    cache_incr(io, s);
    
    /* Remove from relevant seq_block array */
    cache_item_remove(io, GT_Seq, pos->rec);

    io->iface->seq.index_del(io->dbh, s->name);

    /* Deallocate seq struct itself */
    if (!(b = cache_search(io, GT_Bin, s->bin))) {
	cache_decr(io, s);
	return -1;
    }

    /* Remove from range array */
    r = arrp(range_t, b->rng, s->bin_index);
    if (!(r->flags & GRANGE_FLAG_UNUSED)) {
	assert(r->rec == s->rec);

	b = cache_rw(io, b);
	b->flags |= BIN_RANGE_UPDATED;
	r->flags |= GRANGE_FLAG_UNUSED;
    }

    cache_decr(io, s);

    //Already achieved via cache_item_remove
    //cache_rec_deallocate(io, GT_Seq, pos->rec);


    /* TODO: remove from btree name index */

    /* Remove annotations too */
    c = cache_search(io, GT_Contig, pos->contig);
    cache_incr(io, c);
    for (i = 0; i < pos->n_anno; i++) {
	bin_remove_item(io, &c, GT_AnnoEle, pos->anno[i].rec);
	cache_item_remove(io, GT_AnnoEle, pos->anno[i].rec);
	//cache_rec_deallocate(io, GT_AnnoEle, pos->anno[i].rec);
    }
    cache_decr(io, c);

    return 0;
}

/*
 * Form a new contig from the reads in pos[].
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int create_contig_from(GapIO *io, r_pos_t *pos, int npos) {
    int i;
    contig_t *c_old, *c_new;
    char name[8192];
    static int last_count=1;
    static tg_rec last_contig=-1;
    int offset;
    bin_index_t *bin;
    int old_comp;

    if (npos <= 0)
	return -1;

    vmessage("\n=== new contig ===");
    for (i = 0; i < npos; i++) {
	vmessage("%d\tCtg %"PRIrec"\t%d..%d\tseq %"PRIrec"\n",
		 i, pos[i].contig, pos[i].start, pos[i].end, pos[i].rec);
    }


    /* Pick a new contig name */
    c_old = cache_search(io, GT_Contig, pos[0].contig);
    cache_incr(io, c_old);
    if (last_contig != pos[0].contig) {
	last_contig = pos[0].contig;
	last_count = 1;
    }
    do {
	sprintf(name, "%s%%%d", contig_get_name(&c_old), last_count++);
    } while (contig_index_query(io, name) != -1);

    bin = cache_search(io, GT_Bin, c_old->bin);
    old_comp = (bin->flags & BIN_COMPLEMENTED) ? 1 : 0;


    /* Create the new contig */
    c_new = contig_new(io, name);
    contig_index_update(io, name, strlen(name), c_new->rec);
    cache_incr(io, c_new);



    /* Add in the sequences */
    offset = pos[0].start - 1;
    for (i = 0; i < npos; i++) {
	range_t r, *r_out;
	seq_t *s;
	int j;

	r = pos[i].rng;
	r.start = pos[i].start - offset;
	r.end   = pos[i].end - offset;
	r.y     = 0;
	if (old_comp)
	    r.flags ^= GRANGE_FLAG_COMP1;

	bin = bin_add_range(io, &c_new, &r, &r_out, NULL, 0);

	s = cache_search(io, GT_Seq, pos[i].rec);
	s = cache_rw(io, s);
	s->bin = bin->rec;
	s->bin_index = r_out - ArrayBase(range_t, bin->rng);
	if (old_comp) {
	    s->len = -s->len;
	    s->flags ^= SEQ_COMPLEMENTED;
	}

	/* Similarly move the annotations too; much the same as seqs */
	for (j = 0; j < pos[i].n_anno; j++) {
	    anno_ele_t *a;

	    //printf("Seq %d; pos %d,   anno %d; pos %d-%d\n",
	    //	   pos[i].rec,
	    //	   pos[i].start,
	    //	   pos[i].anno[j].rec,
	    //	   pos[i].anno[j].start,
	    //	   pos[i].anno[j].end);

	    bin_remove_item(io, &c_old, GT_AnnoEle, pos[i].anno[j].rec);
	    r.start    = pos[i].anno[j].start - offset;
	    r.end      = pos[i].anno[j].end - offset;
	    r.rec      = pos[i].anno[j].rec;
	    r.mqual    = pos[i].anno[j].mqual;
	    r.pair_rec = pos[i].anno[j].pair_rec;
	    r.flags    = pos[i].anno[j].flags;
	    
	    bin = bin_add_range(io, &c_new, &r, &r_out, NULL, 0);
	    a = cache_search(io, GT_AnnoEle, pos[i].anno[j].rec);
	    a = cache_rw(io, a);
	    a->bin = bin->rec;
	    //a->bin_idx = r_out - ArrayBase(range_t, bin->rng);
	}
    }

    cache_decr(io, c_old);
    cache_decr(io, c_new);

    return 0;
}

/*
 * Moves a bunch of reads to a new contig. We determine the set of reads
 * based on overlaps. Calls create_contig_from to do the grunt work.
 */
static int move_reads(GapIO *io, r_pos_t *pos, int npos) {
    int i, start, end;
    tg_rec contig;
    int i_start, err = 0;

    if (!npos)
	return 0;

    i_start = 0;
    contig  = pos[0].contig;
    start   = pos[0].start;
    end     = pos[0].end;
    for (i = 1; i < npos; i++) {
	if (pos[i].contig != contig ||
	    pos[i].start > end) {
	    if (create_contig_from(io, &pos[i_start], i - i_start))
		err = 1;
	    i_start = i;
	    contig  = pos[i].contig;
	    start   = pos[i].start;
	    end     = pos[i].end;
	} else {
	    if (end < pos[i].end)
		end = pos[i].end;
	}
    }
    if (create_contig_from(io, &pos[i_start], i - i_start))
	err = 1;

    return err;
}

/*-----------------------------------------------------------------------------
 * External interfaces
 */


/**
 * Removes a set of readings from either the contig or the database.
 *
 * When removing from the database we need to delete everything related to
 * that reading (annotations, template if not used elsewhere, etc). When
 * moving to a new contig, we prefer to keep all readings clustered together.
 *
 * Ie if we remove A & B (overlapping) from one contig and C from another
 * then we create two new contigs containing A & B in one and C in the other.
 *
 * When creating new contigs, we have the option of copying over any
 * overlapping consensus tags to the new contigs. This choice only refers to
 * consensus tags; reading tags are always copied.
 *
 * move == 0   => remove
 * move == 1   => split to new single-read contigs
 * move == 2   => move to new still-joined contigs
 *
 * Returns 0 on success
 *        -1 on failure
 */
int disassemble_readings(GapIO *io, tg_rec *rnums, int nreads, int move,
                         int remove_holes, int duplicate_tags)
{
    int i,err = 0;
    r_pos_t *pos;
    HacheTable *dup_hash;

    vfuncheader("Disassemble_readings");
    //    check_contig_bins(io);

    if (nreads <= 0)
	return 0;

    /*
     * The plan:
     *
     * 1. Produce a table of which contig number and region each reading
     *    record is within.
     *
     * 2. Remove these readings from their contigs.
     * 2b. Ensure if it's a read-pair that the other end gets the
     *     pairing information and flags cleared.
     *
     * 3a. move==0: deallocate the reading record.
     * 3b. move==1: produce new contigs/bins for each read.
     * 3c. move==2: produce 1 new contig for each overlapping set of
     *              disassembled reads.
     *
     * 4. Check for holes within the "source" contigs. We use the
     *    initial table for this, along with regions to look around.
     *
     */

    if (NULL == (pos = calloc(nreads, sizeof(*pos))))
	return -1;

    /*
     * Remove from contig bin.
     * Also handles accumulation of contig/positions and removal of
     * duplicate reads.
     */
    dup_hash = HacheTableCreate(8192, HASH_DYNAMIC_SIZE |
				      HASH_NONVOLATILE_KEYS);
    for (i = 0; i < nreads; i++) {
	HacheData hd;
	int n;

	/* Part 1, 2 and 2b: remove reads */
	hd.i = 0;
	HacheTableAdd(dup_hash, (char *)&rnums[i], sizeof(rnums[i]), hd, &n);
	if (!n) {
	    gio_debug(io, 1, "Skipping duplicate entry %"PRIrec"\n", rnums[i]);
	    pos[i].contig = 0;
	    continue;
	}
	if (rnums[i] == -1) {
	    continue;
	}
	if (unlink_read(io, rnums[i], &pos[i], move == 0))
	    return -1;

	pos[i].rec = rnums[i];
    }

    HacheTableDestroy(dup_hash, 0);


    /* Sort position table and drop the duplicate entries */
    xio = io;
    qsort(pos, nreads, sizeof(*pos), pos_sort);
    for (i = 0; i < nreads; i++)
	if (pos[i].contig)
	    break;
    if (i != 0) {
	memmove(&pos[0], &pos[i], (nreads-i) * sizeof(pos[0]));
	nreads -= i;
    }

    /* Part 3: */
    switch (move) {
    case 0:
	/* 3a. Deallocate the record */
	for (i = 0; i < nreads; i++) {
	    if (seq_deallocate(io, &pos[i]))
		err = 1;
	}
	break;

    case 1:
	/* 3b. produce new contigs/bins for each read. */
	for (i = 0; i < nreads; i++) {
	    if (create_contig_from(io, &pos[i], 1))
		err = 1;
	}
	break;

    case 2:
	/* 3c. produce 1 new contig for each overlapping set of
	       disassembled reads. */
	if (move_reads(io, pos, nreads))
	    err = 1;
	break;

    default:
	fprintf(stderr, "Unexpected 'move' value %d\n", move);
	err = 1;
	return -1;
    }

    cache_flush(io);

    /* Part 4: fix holes in source contigs */
    if (fix_holes(io, pos, nreads, remove_holes))
	/* Too late to undo, so keep going! */
	err = 1;

    cache_flush(io);

    // check_contig_bins(io);

    for (i = 0; i < nreads; i++)
	if (pos[i].anno)
	    free(pos[i].anno);
    free(pos);

    return err ? -1 : 0;
}
