#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "tg_gio.h"

/*
 * Allocates a new annotation element.
 * Returns 0 for success
 *        -1 for failure.
 */
int anno_ele_new(GapIO *io, int bin,
		 int obj_type, int obj_rec, int anno_rec,
		 int type, char *comment) {
    int rec;
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
 * Returns the range_t element from the bin holding this annotation.
 * The start and end have been modified to be the absolute position
 * within the contig.
 *
 * Returns a static range_t pointer on success (valid until next call)
 *         NULL on failure.
 */
range_t *anno_get_range(GapIO *io, int anno_ele, int *contig) {
    anno_ele_t *e = (anno_ele_t *)cache_search(io, GT_AnnoEle, anno_ele);
    bin_index_t *bin;
    int bnum, offset1, offset2, comp = 0, i;
    range_t *r;
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

    return &r2;
}
