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
 * Returns the range_t element from the bin holding this annotation.
 * The start and end have been modified to be the absolute position
 * within the contig.
 *
 * Returns a static range_t pointer on success (valid until next call)
 *         NULL on failure.
 */
range_t *anno_get_range(GapIO *io, int anno_ele, int *contig);

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
