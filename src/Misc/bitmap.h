/*
 * File: bitmap.h
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: header for bitmap management routines
 *
 * Created: 2 October 1992
 * Updated:
 *
 */
#ifndef _BITMAP_H_
#define _BITMAP_H_

#include <stdio.h> /* import: FILE */
#include "os.h"
#include "xerror.h"

/* Bitmap code */
#define BITMAP_ERR_START		100
#define BITMAP_NO_ERROR			0
#define BITMAP_FULL			(BITMAP_ERR_START + 0)
#define BITMAP_INVALID_ARGUMENTS	(BITMAP_ERR_START + 1)
#define BITMAP_OUT_OF_MEMORY		(BITMAP_ERR_START + 2)

#define BASE_TYPE uint4

typedef struct {
    BASE_TYPE *base;
    int Nbitmap; /* Number of elements **not bits** allocated */
    int Nbits;   /* Number of bits used */
    int first_free; /* for optimization */
} BitmapStruct, *Bitmap;


/*
 * Bitmap flag modes
 */
#define BITMAP_NONE       0L
#define BITMAP_EXTENDABLE 1L /* ignored */

/*
 * Bitmap constants
 */
#define BIT_CHR 8 /* number of bits per character */
#define CHR_ELE ((int)(sizeof(BASE_TYPE))) /* characters per bitmap element */
#define BIT_ELE ((int)(BIT_CHR * CHR_ELE)) /* bits per bitmap element */
#define ELE_NBITS(B) ((int)(((B)+BIT_ELE-1)/BIT_ELE)) /* number of elements to store B bits */
#define MIN_ELE 16 /* minimum number of elements: MUST BE GREATER THAN 0 */
#define INC_ELE 16 /* when extending bitmap create this number of spare elemanents */

/*
 * Bitmap macros
 */
#define BIT_IDX(I) ((I) / BIT_ELE)
#define BIT_MSK(I) (1 << ((I) % BIT_ELE))
#define BIT_CLR(B,I)  ((B)->first_free=MIN((B)->first_free,I), \
		       (B)->base[BIT_IDX(I)] &= ~BIT_MSK(I))
#define BIT_SET(B,I)  (B)->base[BIT_IDX(I)] |= BIT_MSK(I)
#define BIT_CHK(B,I) ((B)->base[BIT_IDX(I)] & BIT_MSK(I))
#define BIT_FREE(B) (BitmapFree(B))
#define BIT_NBITS(B) (B)->Nbits


#define berr_set(e) ((e)?xerr_set((e), BitmapErrorString((e))):0)

/*
 * Function prototypes
 */

/*
 * Create a bit map for at least Nbits bits
 * Returns: pointer to a Bitmap if one is created, else NULL
 */
extern Bitmap BitmapCreate(int Nbits);

/*
 * Destroy a bitmap
 */
extern int BitmapDestroy(Bitmap bitmap);

/*
 * Extend a bitmap so as to include newNbits.
 * If the bitmap needs to be extended to accomodate newNbits, but
 * cannot because of memory problems, an error is returned and bitmap is
 * unaffected
 */
extern int BitmapExtend(Bitmap bitmap, int newNbits);

/*
 * Return index of an unset bit in bitmap.
 * Assumption: that the bit returned will be used as set (although this
 * routine does not set it).
 */
extern int BitmapFree(Bitmap bitmap);

/*
 * pretty print a bitmap
 */
extern int BitmapPrint(FILE *fp, Bitmap bitmap);

/*
 * Fetch a bitmap error string
 */
extern char *BitmapErrorString(int err);

/*
 * perform bitwise AND, inclusive OR or exclusive OR  between 
 * bitmap1 and bitmap2. Return combined bitmap
 */
extern Bitmap BitmapNOT(Bitmap bitmap1);
extern Bitmap BitmapAND(Bitmap bitmap1, Bitmap bitmap2);
extern Bitmap BitmapOR (Bitmap bitmap1, Bitmap bitmap2);
extern Bitmap BitmapXOR(Bitmap bitmap1, Bitmap bitmap2);

/*
 * find the nth bit that has been set
 */
extern int FindNBitSet(Bitmap bitmap, int n);

#endif /*_BITMAP_H_*/
