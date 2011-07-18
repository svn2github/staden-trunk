#include <assert.h>
#include <string.h>
#include <math.h>

#include "tg_library.h"
#include "tg_gio.h"

/*
 * integer component of log2(l). NB in this log2(0) = 0;
 *
 * This *modifies* l rather than simply returning it. This is hacky, but
 * it avoids needing temporary variable storage.
 */
#define l2(l) \
    l|=l>>1;\
    l|=l>>2;\
    l|=l>>4;\
    l|=l>>8;\
    l|=l>>16;\
    l=((l/2+1)*0x77CB531U)>>27;\
    l="\x0\x1\x1c\x2\x1d\xe\x18\x3\x1e\x16\x14\xf\x19\x11\x4\x8\x1f\x1b\xd\x17\x15\x13\x10\x7\x1a\xc\x12\x6\xb\x5\xa\x9"[l]

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
tg_rec library_new(GapIO *io, char *name) {
    tg_rec rec;
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
    lib->flags = 0;

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
	lib->counts[0] = 0;
	memset(lib->size_hist[i], 0, LIB_BINS * sizeof(lib->size_hist[i][0]));
    }

    /* Add it to the global library array too */
    io->library = cache_rw(io, io->library);
    io->db = cache_rw(io, io->db);
    ARR(tg_rec, io->library, io->db->Nlibraries++) = rec;

    return rec;
}

int accumulate_library_rec(GapIO *io, tg_rec rec, int type, int size) {
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

/*
 * Computes the mean and standard deviation of a library along with which
 * type. See LIB_T_* macros in tg_struct.h.
 * 'min_count' specifies a minimum number of paired reads to have before
 * writing back the data to the library struct. Specify this as zero
 * to do it regardless.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int update_library_stats(GapIO *io, tg_rec rec, int min_count,
			 double *mean, double *sd, int *type) {
    library_t *lib = cache_search(io, GT_Library, rec);
    int i, j, k;
    /*
    double sum[3]    = {0, 0, 0};
    double sum_sq[3] = {0, 0, 0};
    double m, s;
    */
    double N[3];
    double isize[3], sd_[3], last_sd[3];

    if (!lib)
	return -1;

    /*
    for (i = 0; i < LIB_BINS; i++) {
	int bsize = ibin2isize(i+1);
	for (j = 0; j < 3; j++) {
	    double d = (double)lib->size_hist[j][i];
	    sum[j]    += bsize * d;
	    sum_sq[j] += bsize * bsize * d;
	    N[j] += d;
	}
    }

    if (sum[0] > sum[1]) {
	j = sum[0] > sum[2] ? 0 : 2;
    } else {
	j = sum[1] > sum[2] ? 1 : 2;
    }

    m = sum[j] / N[j];
    s = sqrt((sum_sq[j]/N[j] - m*m));

    printf("OLD MODE m=%f,%f,%f s=%f,%f,%f\n", 
	   sum[0] / N[0], sum[1] / N[1], sum[2] / N[2],
	   sqrt((sum_sq[0]/N[0] - sum[0]*sum[0]/(N[0]*N[0]))),
	   sqrt((sum_sq[1]/N[1] - sum[1]*sum[1]/(N[1]*N[1]))),
	   sqrt((sum_sq[2]/N[2] - sum[2]*sum[2]/(N[2]*N[2]))));
    */

    /*
     * Try 2: compute mean and s.d. via the interquartile range instead.
     * This seems more stable when given very non gaussian looking data or
     * distributions with very long tails.
     *
     * Skip tiniest bins as these normally in error atm.
     */
#define SMALLEST_BIN 5
    for (k = 0; k < 10; k++) {
	for (j = 0; j < 3; j++) {
	    double count = 0, c2 = 0;
	    double q1, q2, q3;
	    int bin1, bin2;

	    if (k == 0) {
		bin1 = 0;
		bin2 = LIB_BINS;
	    } else {
		if (isize[j] - 3*sd_[j] > 0)
		    bin1 = isize2ibin(isize[j] - 3*sd_[j]);
		else
		    bin1 = 0;
	    
		bin2 = isize2ibin(isize[j] + 3*sd_[j]);
	    }

	    for (i = bin1; i < bin2; i++) {
		count += (double)lib->size_hist[j][i];
	    }
	    if (k == 0)
		N[j] = count;

	    q1 = q2 = q3 = 0;

	    for (i = bin1; i < bin2; i++) {
		q1 = ibin2isize(i+1);
		c2 += (double)lib->size_hist[j][i];
		if (c2 >= .25 * count)
		    break;
	    }
	    for (; i < bin2; i++) {
		q2 = ibin2isize(i+1);
		c2 += (double)lib->size_hist[j][i];
		if (c2 >= .5 * count)
		    break;
	    }
	    for (; i < bin2; i++) {
		q3 = ibin2isize(i+1);
		c2 += (double)lib->size_hist[j][i];
		if (c2 >= .75 * count)
		    break;
	    }

	    isize[j] = q2;
	    sd_[j] = (q3-q1)/1.349;
	}

	/*
	printf("median=%d,%d,%d, sd=%d,%d,%d\n",
	       (int)isize[0], (int)isize[1], (int)isize[2],
	       (int)sd_[0], (int)sd_[1], (int)sd_[2]);
	*/

	if (k && (int)last_sd[0] == (int)sd_[0]
	      && (int)last_sd[1] == (int)sd_[1]
	      && (int)last_sd[2] == (int)sd_[2])
	    break;

	last_sd[0] = sd_[0];
	last_sd[1] = sd_[1];
	last_sd[2] = sd_[2];
    }
    //puts("");

    if (N[0] > N[1]) {
	j = N[0] > N[2] ? 0 : 2;
    } else {
	j = N[1] > N[2] ? 1 : 2;
    }

    if (type) *type = j;
    if (mean) *mean = isize[j];
    if (sd)   *sd   = sd_[j];

    /*
     * We update this data anyway so read only versions get to see
     * on-the-fly computations, but if they've changed and we can lock it rw
     * then we also permit archival of these figures.
     */
    if (N[0] + N[1] + N[2] >= min_count) {
	int edited = 0;

	if (lib->lib_type != j)
	    edited = 1;

	for (i = 0; i < 3; i++) {
	    if (lib->insert_size[i] != isize[i])
		edited = 1;
	    if (fabs(sd_[i] - lib->sd[i]) > 0.01)
		edited = 1;
	}

	if (edited) {
	    library_t *elib;
	    if ((elib = cache_rw(io, lib)) != NULL)
		lib = elib;
	}

	lib->lib_type = j;

	for (i = 0; i < 3; i++) {
	    lib->insert_size[i] = isize[i];
	    lib->sd[i] = sd_[i];
	    lib->counts[i] = N[i];
	}

	lib->flags = 1;
    }

    return 0;
}

int get_library_stats(GapIO *io, tg_rec rec,
		      double *mean, double *sd, int *type, int *count) {
    library_t *lib = cache_search(io, GT_Library, rec);
    int i, j;
    double N[3];

    if (!lib)
	return -1;

    for (j = 0; j < 3; j++) {
	N[j] = 0;
	for (i = 0; i < LIB_BINS; i++) {
	    N[j] += lib->size_hist[j][i];
	}
    }
    
    if (N[0] > N[1]) {
	j = N[0] > N[2] ? 0 : 2;
    } else {
	j = N[1] > N[2] ? 1 : 2;
    }

    if (mean)
	*mean = lib->insert_size[j];
    if (sd)
	*sd = lib->sd[j];
    if (type)
	*type = j;
    if (count)
	*count = N[j];

    return 0;
}
