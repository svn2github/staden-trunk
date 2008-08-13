#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "tg_gio.h"

/*
 * Sets the track type
 *
 * Returns 0 on success
 *        -1 on failure
 */
int track_set_type(GapIO *io, track_t **t, int value) {
    track_t *n;
    if (!(n = cache_rw(io, *t)))
	return -1;

    n->type = value;
    *t = n;

    return 0;
}

int track_set_flag(GapIO *io, track_t **t, int value) {
    track_t *n;
    if (!(n = cache_rw(io, *t)))
	return -1;

    n->flag = value;
    *t = n;

    return 0;
}

int track_set_bin_size(GapIO *io, track_t **t, int value) {
    track_t *n;
    if (!(n = cache_rw(io, *t)))
	return -1;

    n->bin_size = value;
    *t = n;

    return 0;
}

int track_set_item_size(GapIO *io, track_t **t, int value) {
    track_t *n;
    if (!(n = cache_rw(io, *t)))
	return -1;

    n->item_size = value;
    *t = n;

    return 0;
}

int track_set_nitems(GapIO *io, track_t **t, int value) {
    track_t *n;
    if (!(n = cache_rw(io, *t)))
	return -1;

    n->nitems = value;
    *t = n;

    return 0;
}

int track_set_data(GapIO *io, track_t **t, Array value) {
    track_t *n;
    if (!(n = cache_rw(io, *t)))
	return -1;

    if (n->data)
	ArrayDestroy(n->data);

    n->data = value;
    *t = n;

    return 0;
}
