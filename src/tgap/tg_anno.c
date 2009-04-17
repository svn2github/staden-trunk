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
		 char *comment) {
    int rec;
    anno_ele_t e;

    e.bin      = bin;
    e.obj_type = obj_type;
    e.obj_rec  = obj_rec;
    e.anno_rec = anno_rec;
    e.comment  = comment;
    
    if (-1 == (rec = io->iface->anno_ele.create(io->dbh, &e)))
	return -1;

    return rec;
}
