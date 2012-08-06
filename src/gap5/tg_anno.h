#ifndef _TG_ANNO_H_
#define _TG_ANNO_H_

#include "io_utils.h"

/*
 * Allocates a new annotation element.
 * Returns 0 for success
 *        -1 for failure.
 */
tg_rec anno_ele_new(GapIO *io, tg_rec bin,
		    int obj_type, tg_rec obj_rec, tg_rec anno_rec,
		    int type, char dir, char *comment);

/*
 * Removes an anno_ele from the gap database.
 * FIXME: need to deallocate storage too. (See docs/TODO)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int anno_ele_destroy(GapIO *io, anno_ele_t *e);

/*
 * Creates an anno_ele as per anno_ele_new, but also adds it to an object
 * and creates the bin Range entry too.
 */
tg_rec anno_ele_add(GapIO *io, int obj_type, tg_rec obj_rec, tg_rec anno_rec,
		    int type, char *comment, int start, int end, char dir);

/*
 * Returns the range_t element from the bin holding this annotation.
 * The start and end have been modified to be the absolute position
 * within the contig, unless 'rel' is true.
 *
 * Returns a static range_t pointer on success (valid until next call)
 *         NULL on failure.
 */
range_t *anno_get_range(GapIO *io, tg_rec anno_ele, tg_rec *contig, int rel);

/*
 * Sets the comment for an annotation element.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int anno_ele_set_comment(GapIO *io, anno_ele_t **e, char *comment);

/*
 * Sets the annotation type, passed in as a string but held in a 4-byte int.
 * This also attempts to set the cached copy of the type held within the
 * bin range array.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int anno_ele_set_type(GapIO *io, anno_ele_t **e, char *str);

/*
 * Sets the annotation direction, one of ANNO_DIR_* macros (+,-,.,?)
 *
 * Returns 0 on success
 *        -1 on failure
 */
int anno_ele_set_direction(GapIO *io, anno_ele_t **e, char dir);

/*
 * Finds the contig number and position of an anno_ele record number.
 *
 * If non-NULL r_out is filled with the associated range_t struct.
 *
 * If non-NULL a_out is filled with a pointer to the anno_ele_t struct.
 * This will have had cache_incr() run on it, so the caller should
 * use cache_decr() to permit deallocation.
 *
 * FIXME: this is almost identical to seq_get_position2, so maybe they
 * should both call a common function in tg_bin.c instead.
 */
int anno_get_position2(GapIO *io, tg_rec anum, tg_rec *contig,
		       int *start, int *end, int *orient,
		       range_t *r_out, seq_t **a_out);

int anno_get_position(GapIO *io, tg_rec anum, tg_rec *contig,
		      int *start, int *end, int *orient);


/*
 * Removes some or all tags from some or all contigs.
 * If the contig list or tag list is blank it implies all contigs or all tags.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int delete_tags(GapIO *io, int ncontigs, contig_list_t *contigs,
		char *tag_list, int verbose);

#define str2type(s) ((s)[3] + ((s)[2]<<8) + ((s)[1]<<16) + ((s)[0]<<24))
#define type2str(t,s) \
    ( \
     ((s)[0] = ((t) >> 24) & 0xff), \
     ((s)[1] = ((t) >> 16) & 0xff), \
     ((s)[2] = ((t) >>  8) & 0xff), \
     ((s)[3] = ((t) >>  0) & 0xff), \
     ((s)[4] = '\0'), \
     (s))

#endif /* _TG_ANNO_H_ */
