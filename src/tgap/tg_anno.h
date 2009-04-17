#ifndef _TG_ANNO_H_
#define _TG_ANNO_H_

/*
 * Allocates a new annotation element.
 * Returns 0 for success
 *        -1 for failure.
 */
int anno_ele_new(GapIO *io, int bin,
		 int obj_type, int obj_rec, int anno_rec,
		 char *comment);

#endif /* _TG_ANNO_H_ */
