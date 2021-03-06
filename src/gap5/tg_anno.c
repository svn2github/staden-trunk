#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "tg_gio.h"
#include "active_tags.h"
#include "io_utils.h"
#include "io_lib/hash_table.h"
#include "gap4_compat.h" /* NumContigs() */
#include "text_output.h" /* UpdateTextOutput() */

/*
 * Allocates a new annotation element.
 * Returns rec for success
 *        -1 for failure.
 */
tg_rec anno_ele_new(GapIO *io, tg_rec bin,
		    int obj_type, tg_rec obj_rec, tg_rec anno_rec,
		    int type, char dir, char *comment) {
    tg_rec rec;
    anno_ele_t e;

    e.bin       = bin;
    e.obj_type  = obj_type;
    e.obj_rec   = obj_type == GT_Contig ? 0 : obj_rec;
    e.anno_rec  = anno_rec;
    e.tag_type  = type;
    e.direction = dir;
    e.comment   = comment;
    
    rec = cache_item_create(io, GT_AnnoEle, &e);

    return rec;
}

/*
 * Creates an anno_ele as per anno_ele_new, but also adds it to an object
 * and creates the bin Range entry too.
 */
tg_rec anno_ele_add(GapIO *io, int obj_type, tg_rec obj_rec, tg_rec anno_rec,
		    int type, char *comment, int start, int end, char dir) {
    range_t r;
    anno_ele_t *e;
    contig_t *c;
    tg_rec crec;
    bin_index_t *bin;
    tg_rec seq_bin = 0;

    /* Find contig for obj_rec/obj_type */
    if (obj_type == GT_Contig) {
	crec = obj_rec;
    } else {
	int st, en;
	sequence_get_position2(io, obj_rec, &crec, &st, &en, NULL,
			       &seq_bin, NULL, NULL);

	start += st;
	end += st;
    }

    c = (contig_t *)cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

    r.start    = start;
    r.end      = end;
    r.flags    = GRANGE_FLAG_ISANNO;
    r.mqual    = type;
    r.pair_rec = obj_rec;

    if (GT_Seq == obj_type)
	r.flags |= GRANGE_FLAG_TAG_SEQ;

    r.rec = anno_ele_new(io, 0, obj_type, obj_rec, 0, type, dir, comment);
    e = (anno_ele_t *)cache_search(io, GT_AnnoEle, r.rec);
    e = cache_rw(io, e);

    if (seq_bin)
	bin = bin_add_to_range(io, &c, seq_bin, &r, NULL, NULL, 0);
    else
	bin = bin_add_range(io, &c, &r, NULL, NULL, 0);

    if (!bin) 
	verror(ERR_FATAL, "anno_ele_add", "bin_add_to_range returned NULL");

    e->bin = bin ? bin->rec : 0;

    cache_decr(io, c);
    return r.rec;
}

#if 0 /* Obsoleted by bin_remove_item */
/*
 * Removes an anno_ele from the gap database.
 * FIXME: need to deallocate storage too. (See docs/TODO)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int anno_ele_destroy(GapIO *io, anno_ele_t *e) {
    bin_index_t *bin;
    range_t *r;
    int i;

    /* Find the bin range pointing to this object */
    bin = (bin_index_t *)cache_search(io, GT_Bin, e->bin);
    if (!bin || !bin->rng || ArrayMax(bin->rng) == 0)
	return -1;
    if (!(bin = cache_rw(io, bin)))
	return -1;


    for (i = 0; i < ArrayMax(bin->rng); i++) {
	r = arrp(range_t, bin->rng, i);
	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if (r->rec == e->rec)
	    break;
    }
    if (i == ArrayMax(bin->rng))
	return -1;

    /* Mark this bin range as unused */
    r->rec = bin->rng_free;
    r->flags |= GRANGE_FLAG_UNUSED;

    bin->rng_free = i;
    bin->flags |= BIN_RANGE_UPDATED | BIN_BIN_UPDATED;

    bin_incr_nanno(io, bin, -1);

    if (bin->start_used == r->start || bin->end_used == r->end)
	bin_set_used_range(io, bin);

    return 0;
}
#endif

/*
 * Sets the comment for an annotation element.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int anno_ele_set_comment(GapIO *io, anno_ele_t **e, char *comment) {
    anno_ele_t *ae;
    size_t clen;

    if (!(ae = cache_rw(io, *e)))
	return -1;

    clen = comment ? strlen(comment) : 0;
    if (clen > (ae->comment ? strlen(ae->comment) : 0)) {
	ae = cache_item_resize(ae, sizeof(*ae) + clen+1);
	ae->comment = (char *)&ae->data;
    }
    if (clen)
	strcpy(ae->comment, comment);

    *e = ae;

    return 0;
}

/*
 * Sets the annotation direction, one of ANNO_DIR_* macros (+,-,.,?)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int anno_ele_set_direction(GapIO *io, anno_ele_t **e, char dir) {
    anno_ele_t *ae;

    if (!(ae = cache_rw(io, *e)))
	return -1;

    *e = ae;
    ae->direction = dir;

    return 0;
}

/*
 * Sets the annotation type, passed in as a string but held in a 4-byte int.
 * This also attempts to set the cached copy of the type held within the
 * bin range array.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int anno_ele_set_type(GapIO *io, anno_ele_t **e, char *str) {
    int type;
    char stype[5];
    anno_ele_t *ae;

    if (!(ae = cache_rw(io, *e)))
	return -1;

    /* Get integer type */
    memset(stype, 0, 5);
    strncpy(stype, str, 4);
    type = str2type(stype);

    /* Update annotation */
    ae->tag_type = type;

    /* Also update range_t cached copy of type */
    if (ae->bin) {
	bin_index_t *bin = (bin_index_t *)cache_search(io, GT_Bin, ae->bin);
	range_t *r = NULL;
	int i, nranges;

	if (!bin)
	    return -1;
	if (!(bin = cache_rw(io, bin)))
	    return -1;

	/*
	 * Find the index into the bin range.
	 * FIXME: we should add a bin_index element, as seen in seq_t,
	 * to avoid the brute force loop. This doesn't have to be
	 * permanently stored - a cached copy would suffice.
	 */
	nranges = bin->rng ? ArrayMax(bin->rng) : 0;
	for (i = 0; i < nranges; i++) {
	    r = arrp(range_t, bin->rng, i);
	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;

	    if (r->rec == ae->rec)
		break;
	}
	if (i == nranges)
	    return -1;

	bin->flags |= BIN_RANGE_UPDATED;
	r->mqual = type;
    }

    *e = ae;

    return 0;
}

/*
 * Returns the range_t element from the bin holding this annotation.
 * The start and end have been modified to be the absolute position
 * within the contig, unless 'rel' is true.
 *
 * Returns a static range_t pointer on success (valid until next call)
 *         NULL on failure.
 */
range_t *anno_get_range(GapIO *io, tg_rec anno_ele, tg_rec *contig, int rel) {
    anno_ele_t *e = (anno_ele_t *)cache_search(io, GT_AnnoEle, anno_ele);
    bin_index_t *bin;
    tg_rec bnum;
    int offset1, offset2, comp = 0, i;
    range_t *r = NULL;
    static range_t r2;

    if (!e)
	return NULL;

    /* Find and load the associated bin */
    bnum = e->bin;
    bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
    if (!bin->rng || ArrayMax(bin->rng) == 0)
	return NULL;

    /* Find the appropriate range_t element */
    for (i = 0; i < ArrayMax(bin->rng); i++) {
	r = arrp(range_t, bin->rng, i);
	if (r->flags & GRANGE_FLAG_UNUSED)
	    continue;

	if (r->rec == anno_ele)
	    break;
    }
    if (r->rec != anno_ele)
	return NULL;

    /* Find the position and orientation of this bin relative to the contig */
    offset1 = r->start;
    offset2 = r->end;
    r2 = *r; /* Copy now as loop below could push *r out of our cache */
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

    r2.start = offset1;
    r2.end   = offset2;

    if (contig)
	*contig = bin->parent;

    if (rel && e->obj_type == GT_Seq) {
	int st, en, orient;
	sequence_get_position(io, e->obj_rec, NULL, &st, &en, &orient);
	r2.start -= st;
	r2.end   -= st;
    }

    if (r2.start > r2.end) {
	int tmp = r2.start;
	r2.start = r2.end;
	r2.end = tmp;
    }

    return &r2;
}

/*
 * Finds the contig number and position of an anno_ele record number.
 *
 * If non-NULL r_out is filled with the associated range_t struct.
 *
 * If non-NULL a_out is filled with a pointer to the anno_ele_t struct.
 * This will have had cache_incr() run on it, so the caller should
 * use cache_decr() to permit deallocation.
 */
int anno_get_position2(GapIO *io, tg_rec anum, tg_rec *contig,
		       int *start, int *end, int *orient,
		       range_t *r_out, seq_t **a_out) {
    return bin_get_item_position(io, GT_AnnoEle, anum,
				 contig, start, end, orient, NULL,
				 r_out, (void **)a_out);
}

int anno_get_position(GapIO *io, tg_rec anum, tg_rec *contig,
		      int *start, int *end, int *orient) {
    return bin_get_item_position(io, GT_AnnoEle, anum,
				 contig, start, end, orient, NULL,
				 NULL, NULL);
}


#if 0
/*
 * Returns the relative orientation of this annotation vs the contig.
 */
int anno_get_orient(GapIO *io, tg_rec anum) {
    bin_index_t *bin = NULL;
    tg_rec bnum;
    anno_ele_t *a = cache_search(io, GT_AnnoEle, anum);
    int comp = a->orient;

    /* Bubble up bins until we hit the root */
    for (bnum = a->bin; bnum; bnum = bin->parent) {
	bin = (bin_index_t *)cache_search(io, GT_Bin, bnum);
	if (bin->flags & BIN_COMPLEMENTED)
	    comp ^= 1;
	if (bin->parent_type != GT_Bin)
	    break;
    }

    assert(bin && bin->parent_type == GT_Contig);
    return comp;
}
#endif


/*
 * Removes all tags of specific types (hashed in h, or all if h == NULL)
 * from a specified contig.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int delete_tag_single_contig(GapIO *io, tg_rec crec,
				    HashTable *h, int verbose) {
    contig_iterator *ci;
    rangec_t *r;
    contig_t *c;
    int ret = -1;

    ci = contig_iter_new_by_type(io, crec, 1, CITER_FIRST,
				 CITER_CSTART, CITER_CEND,
				 GRANGE_FLAG_ISANNO);
    if (!ci)
	return -1;
    
    if (!(c = cache_search(io, GT_Contig, crec))) {
	contig_iter_del(ci);
	return -1;
    }
    cache_incr(io, c);

    while (NULL != (r = contig_iter_next(io, ci))) {
	char t[5];
	(void)type2str(r->mqual, t);
	if (!h || HashTableSearch(h, t, 4)) {
	    anno_ele_t *e;

	    if (verbose)
		vmessage("Removing anno %s #%"PRIrec"\tContig %s\t%d..%d\n",
			 t, r->rec, c->name, r->start, r->end);
	    if (bin_remove_item(io, &c, GT_AnnoEle, r->rec)) goto fail;
	    /* FIXME: Need to reclaim the GT_AnnoEle record itself */
	}
    }

    ret = 0;
 fail:
    contig_iter_del(ci);
    cache_decr(io, c);

    return ret;
}

/*
 * Removes some or all tags from some or all contigs.
 * If the contig list or tag list is blank it implies all contigs or all tags.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int delete_tags(GapIO *io, int ncontigs, contig_list_t *contigs,
		char *tag_list, int verbose) {
    HashTable *h = NULL;
    int ret = 0;

    /* Hash tag types */
    if (tag_list && *tag_list) {
	int i;
	if (SetActiveTags(tag_list) == -1) {
	    return -1;
	}
	h = HashTableCreate(32, 0);
	for (i = 0; i < number_of_active_tags; i++) {
	    HashData hd;
	    hd.i = 0;
	    HashTableAdd(h, active_tag_types[i], 4, hd, NULL);
	}
    }

    /* Iterate over contig list or all contigs */
    if (verbose)
	vfuncheader("Delete Tags");

    if (ncontigs) {
	int i;

	for (i = 0; i < ncontigs; i++) {
	    contig_t *c = cache_search(io, GT_Contig, contigs[i].contig);
	    vmessage("Scanning contig %d of %d (%s)\n",
		     i+1, ncontigs, c->name);
	    ret |= delete_tag_single_contig(io, contigs[i].contig, h, verbose);
	    UpdateTextOutput();
	    cache_flush(io);
	}

    } else {
	int i;
	tg_rec *order = ArrayBase(tg_rec, io->contig_order);

	for (i = 0; i < NumContigs(io); i++) {
	    contig_t *c = cache_search(io, GT_Contig, order[i]);
	    vmessage("Scanning contig %d of %d (%s)\n",
		     i+1, NumContigs(io), c->name);
	    ret |= delete_tag_single_contig(io, order[i], h, verbose);
	    UpdateTextOutput();
	    cache_flush(io);
	}
    }

    SetActiveTags("");
    if (h)
	HashTableDestroy(h, 0);

    return ret;
}
