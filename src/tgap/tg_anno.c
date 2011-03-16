#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "tg_gio.h"

/*
 * Allocates a new annotation element.
 * Returns rec for success
 *        -1 for failure.
 */
tg_rec anno_ele_new(GapIO *io, tg_rec bin,
		    int obj_type, tg_rec obj_rec, tg_rec anno_rec,
		    int type, char *comment) {
    tg_rec rec;
    anno_ele_t e;

    e.bin      = bin;
    e.obj_type = obj_type;
    e.obj_rec  = obj_rec;
    e.anno_rec = anno_rec;
    e.tag_type = type;
    e.comment  = comment;
    
    rec = cache_item_create(io, GT_AnnoEle, &e);

    return rec;
}

/*
 * Creates an anno_ele as per anno_ele_new, but also adds it to an object
 * and creates the bin Range entry too.
 */
tg_rec anno_ele_add(GapIO *io, int obj_type, tg_rec obj_rec, tg_rec anno_rec,
		    int type, char *comment, int start, int end) {
    range_t r;
    anno_ele_t *e;
    contig_t *c;
    tg_rec crec;
    bin_index_t *bin;

    /* Find contig for obj_rec/obj_type */
    if (obj_type == GT_Contig) {
	crec = obj_rec;
    } else {
	int st, en;
	sequence_get_position(io, obj_rec, &crec, &st, &en, NULL);

	start += st;
	end += st;
    }

    c = (contig_t *)cache_search(io, GT_Contig, crec);

    r.start    = start;
    r.end      = end;
    r.flags    = GRANGE_FLAG_ISANNO;
    r.mqual    = type;
    r.pair_rec = obj_rec;

    if (GT_Seq == obj_type)
	r.flags |= GRANGE_FLAG_TAG_SEQ;

    r.rec = anno_ele_new(io, 0, obj_type, obj_rec, 0, type, comment);
    e = (anno_ele_t *)cache_search(io, GT_AnnoEle, r.rec);
    e = cache_rw(io, e);

    bin = bin_add_range(io, &c, &r, NULL, NULL, 0);
    e->bin = bin->rec;

    return r.rec;
}

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
    if (!bin || !bin->rng)
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

    return 0;
}


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
    strcpy(ae->comment, comment);

    *e = ae;

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
    if (!bin->rng)
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

    r2 = *r;
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

