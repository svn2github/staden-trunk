/*
 * File: g-os.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: definitions for basic gap types
 *
 * Created: 23 October 1992
 * Updated:
 *
 */

#ifndef _G_OS_H_
#define _G_OS_H_

#include "os.h"

typedef int4  GImage; /* index to an image record on disk */
typedef uint4 GTimeStamp;  /* time stamp of an image */
typedef int4  GCardinal; /* large value for indexing arrays, files */
typedef int1  GLock;     /* record lock */
typedef int2  GHFlags;	/* header flags */
typedef int1  GFlags;	/* view flags */
typedef int1  GToggle;
typedef int2  GClient;
typedef int2  GFileN;
typedef int4  GView;




/* set maximum values for each type */
#define MAX_GCardinal (MAXINT4)



/* swapping routines */
#define swap_GImage(S,D)     swap_int4(S,D)
#define swap_GTimeStamp(S,D) swap_int4(S,D)
#define swap_GCardinal(S,D)  swap_int4(S,D)
#define swap_GLock(S,D)      /*swap_int2(S,D)*/
#define swap_GFlags(S,D)     /*swap_int2(S,D)*/
#define swap_GHFlags(S,D)    swap_int2(S,D)
#define swap_GToggle(S,D)    /*swap_int2(S,D)*/



/*
 * Variables of the type GToggle are indexes into the array Index,
 * indicating the currently active image-timestamp combination.
 * Their value if one of (G_NO_TOGGLE,0,1).
 * G_NO_TOGGLE is a special case, only ever used in gfile startup only.
 * It indicates that both timestamps for a record are invalid, and will
 * need fixing during g_garbage_collect().
 */

#define G_NO_TOGGLE ((GToggle) -1)
#define G_NO_IMAGE ((GImage) -1)
#define G_YEAR_DOT ((GTimeStamp) 0)


#endif /*_G_OS_H_*/
