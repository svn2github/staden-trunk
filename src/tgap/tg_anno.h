#ifndef _TG_ANNO_H_
#define _TG_ANNO_H_

/*
 * Allocates a new annotation element.
 * Returns 0 for success
 *        -1 for failure.
 */
int anno_ele_new(GapIO *io, int bin,
		 int obj_type, int obj_rec, int anno_rec,
		 int type, char *comment);

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
int anno_ele_add(GapIO *io, int obj_type, int obj_rec, int anno_rec,
		 int type, char *comment, int start, int end);

/*
 * Returns the range_t element from the bin holding this annotation.
 * The start and end have been modified to be the absolute position
 * within the contig.
 *
 * Returns a static range_t pointer on success (valid until next call)
 *         NULL on failure.
 */
range_t *anno_get_range(GapIO *io, int anno_ele, int *contig);

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
