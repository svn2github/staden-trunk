#ifndef _TG_LIBRARY_
#define _TG_LIBRARY_

#include "misc.h" /* ABS, MIN, MAX */
#include "tg_gio.h"

#define get_lib(io, rec) ((library_t *)cache_search((io), GT_Library, (rec)))

/*
 * Allocates a new library record
 *
 * Returns record number on success
 *         -1 on failure
 */
int library_new(GapIO *io, char *name);

int accumulate_library_rec(GapIO *io, int rec, int type, int size);
void accumulate_library(GapIO *io, library_t *lib, int type, int size);

/*-----------------------------------------------------------------------------
 * Conversions to an from insert-size values to insert-size bins
 */

/*
 * Given a +ve insert size this returns a bin number in the range from
 * 0 to 1792.
 * The scale is such that small inserts map to their own unique bin while
 * larger inserts are binned at lower resolution.
 */
int isize2ibin(int isize);

/*
 * Converts an insert size bin to the smallest absolute integer size that
 * fits within it (ie closest to zero).
 */
int ibin2isize(int ibin);

/*
 * Returns the width of an insert-size bin. Ie the range of insert size
 * values that all map to this one bin.
 */
int ibin_width(int ibin);

#if 0
/*
 * The following are macro implementations of the above C functions, for speed.
 * They're nasty though as they require the names of some temporary variables
 * (working space) to be passed into them.
 *
 * Disabled until we determine if calling the functions is a speed issue.
 */
#  define l2(l) (l="\x0\x1\x1c\x2\x1d\xe\x18\x3\x1e\x16\x14\xf\x19\x11\x4\x8\x1f\x1b\xd\x17\x15\x13\x10\x7\x1a\xc\x12\x6\xb\x5\xa\x9"[((l|=l>>1,l|=l>>2,l|=l>>4,l|=l>>8,l|=l>>16,l/2+1)*0x77CB531U)>>27])
#  define isize2ibin(i,x,y) (x=ABS(i),x=MIN(1<<20,x),y=x,y="00000000123456789:;<=>"[l2(y)]-'0',y=(y*128+(x>>y)),y=i>0?y+1792:1792-y)
#  define ibin2isize(i,x,y) (x=(i)-1792,y=ABS(x)/128-1,y+=(y<0),x=(x+(x>=0?-1:1)*128*y)<<y)
#  define BIN_SZ(i,y) (y=ABS((i)-1792)/128-1,y+=(y<0),y=1<<y)
#endif

/*
 * Computes the mean and standard deviation of a library along with which
 * type. See LIB_T_* macros in tg_struct.h.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int library_stats(GapIO *io, int rec, double *mean, double *sd, int *type);

#endif /* _TG_LIBRARY_ */
