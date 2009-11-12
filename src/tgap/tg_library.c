#include <assert.h>
#include <string.h>

#include "tg_library.h"
#include "tg_gio.h"

/*
 * integer component of log2(l). NB in this log2(0) = 0;
 *
 * This *modifies* l rather than simply returning it. This is hacky, but
 * it avoids needing temporary variable storage.
 */
#define l2(l) (l="\x0\x1\x1c\x2\x1d\xe\x18\x3\x1e\x16\x14\xf\x19\x11\x4\x8\x1f\x1b\xd\x17\x15\x13\x10\x7\x1a\xc\x12\x6\xb\x5\xa\x9"[((l|=l>>1,l|=l>>2,l|=l>>4,l|=l>>8,l|=l>>16,l/2+1)*0x77CB531U)>>27])

/*
 * Given a signed insert size this returns a bin number in the range from
 * 0 to 3584 (with isize==0 mapping to bin 1792).
 * The scale is such that small inserts map to their own unique bin while
 * larger inserts are binned at lower resolution.
 */
int isize2ibin(int isize) {
    int x, y;
    static int ymap[] = {
	0, 0, 0, 0, 0, 0, 0, 0,
	1, 2, 3, 4, 5, 6, 7, 8,
	9,10,11,12,13,14,15,16
    };

    /* Cap size between 0 and 1<<20 */
    x = MAX(0, isize);
    x = MIN(1<<20, x);

    /* y = (int)log2(x),  but approx 10x faster */
    y = x; l2(y);

    /* Map log2(x) so that small values have the same resolution */
    y = ymap[y];

    /* Finally the clever bit - pick the bin given isize and bin scale */
    y = (y*128+(x>>y));

    return y;
}

/*
 * Converts an insert size bin to the smallest absolute integer size that
 * fits within it (ie closest to zero).
 */
int ibin2isize(int ibin) {
    int y;

    /* Reverse of the map in ibin2isize, y is log2(bin_width) */
    y=ibin/128-1; y+=(y<0);

    /* And now convert from the log size to an insert size itself */
    return (ibin - 128*y) << y;
}

/*
 * Returns the width of an insert-size bin. Ie the range of insert size
 * values that all map to this one bin.
 */
int ibin_width(int ibin) {
    int y;

    /* Reverse of the map in ibin2isize, y is log2(bin_width) */
    y=ABS(ibin)/128-1; y+=(y<0);

    return 1<<y;
}

/*
 * Allocates a new library record
 *
 * Returns record number on success
 *         -1 on failure
 */
int library_new(GapIO *io, char *name) {
    int rec;
    library_t *lib;
    int i;

    /* Allocate a record */
    if (-1 == (rec = io->iface->library.create(io->dbh, NULL)))
	return -1;

    /* Initialise the values */
    lib = get_lib(io, rec);
    lib = cache_rw(io, lib);

    lib->rec = rec;
    lib->machine = 0;
    lib->lib_type = 0;

    if (name && *name) {
	lib = cache_item_resize(lib, sizeof(*lib) + strlen(name) + 1);
	lib->name = (char *)&lib->data;
	strcpy(lib->name, name);
    } else {
	lib->name = NULL;
    }

    for (i = 0; i < 3; i++) {
	lib->insert_size[i] = 0;
	lib->sd[i] = 0.0;
	memset(lib->size_hist[i], 0, LIB_BINS * sizeof(lib->size_hist[i][0]));
    }

    /* Add it to the global library array too */
    io->library = cache_rw(io, io->library);
    io->db = cache_rw(io, io->db);
    ARR(GCardinal, io->library, io->db->Nlibraries++) = rec;

    return rec;
}

int accumulate_library_rec(GapIO *io, int rec, int type, int size) {
    library_t *lib = get_lib(io, rec);

    assert(type >= 0 && type <= 2);

    if (NULL == (lib = cache_rw(io, lib)))
	return -1;

    lib->size_hist[type][isize2ibin(size)]++;
    
    return 0;
}

void accumulate_library(GapIO *io, library_t *lib, int type, int size) {
    assert(type >= 0 && type <= 2);

    lib->size_hist[type][isize2ibin(size)]++;
}

