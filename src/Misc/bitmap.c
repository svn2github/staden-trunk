/*
 * File: bitmap.c
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: bitmap management routines
 *
 * Created: 2 October 1992
 * Updated:
 *
 */

#include <stdio.h>

#include "bitmap.h"
#include "xalloc.h"

/* number of bits set to represent 0 to 255 */
int nbits[256] = {
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

#define nb(x) (nbits[(x)[0]]+nbits[(x)[1]]+nbits[(x)[2]]+nbits[(x)[3]])

char *BitmapErrorString(int err)
{
    switch(err) {
    case BITMAP_NO_ERROR:          return "No error";
    case BITMAP_FULL:     	   return "Bitmap full";
    case BITMAP_INVALID_ARGUMENTS: return "Invalid arguments";
    case BITMAP_OUT_OF_MEMORY:     return "Out of memory";
    default:			   return "Unknown error";
    }
}


#define default_elements(B) ( ELE_NBITS(B) < MIN_ELE ? MIN_ELE : ELE_NBITS(B) )

Bitmap BitmapCreate(int Nbits)
/*
 * Create a bit map for at least Nbits bits
 * Returns: pointer to a Bitmap if one is created, else NULL
 */
{
    Bitmap bitmap;
    int i;

    /* check for silly things */
    if ( Nbits < 0 ) {
	(void)berr_set(BITMAP_INVALID_ARGUMENTS);
	return NULL;
    }

    /* allocate bitmap header */
    if ( (bitmap = (Bitmap)xmalloc(sizeof(BitmapStruct))) == NULL ) {
	(void)berr_set(BITMAP_OUT_OF_MEMORY);
	return NULL;
    }


    /* allocate bit array */
    bitmap->Nbits = Nbits;
    bitmap->first_free = (BASE_TYPE)0;
    bitmap->Nbitmap = default_elements(Nbits); /* set default size */
    if ( (bitmap->base = (BASE_TYPE*)xmalloc(bitmap->Nbitmap * CHR_ELE))
	== NULL) {
	xfree(bitmap);
	(void)berr_set(BITMAP_OUT_OF_MEMORY);
	return NULL;
    }
    
    /* initialise bit array */
    for(i=0;i<bitmap->Nbitmap;i++) bitmap->base[i] = (BASE_TYPE)0;

    return bitmap;
}


int BitmapDestroy(Bitmap bitmap)
/*
 * Destroy a bitmap
 */
{
    if (bitmap == NULL) {
	return berr_set(BITMAP_INVALID_ARGUMENTS);
    }

    if (bitmap->base != NULL) {
	xfree(bitmap->base);
	bitmap->base = NULL;
    }
    xfree(bitmap);

    return 0;
}


int BitmapExtend(Bitmap bitmap, int new_Nbits)
/*
 * Extend a bitmap so as to include new_Nbits.
 * If the bitmap needs to be extended to accomodate new_Nbits, but
 * cannot because of memory problems, an error is returned and bitmap is unaffected
 */
{
    BASE_TYPE *new_base;
    int new_Nbitmap;
    int i;

    /* check for silly arguments */
    if (bitmap==NULL) {
	return berr_set(BITMAP_INVALID_ARGUMENTS);
    }

    /* do we have room enough? */

    /* YES */

    /* yes and we don't need to modify anything */
    if (new_Nbits < bitmap->Nbits) {
	return 0;
    }

    /* yes but need to adjust */
    if (ELE_NBITS(new_Nbits) <= bitmap->Nbitmap) {
	bitmap->Nbits = new_Nbits;
	return 0;
    }

    /* NO */

    /* allocate just enough and the just a bit more */
    new_Nbitmap = ELE_NBITS(new_Nbits) + INC_ELE;
    new_base = (BASE_TYPE *)xrealloc(bitmap->base,new_Nbitmap * CHR_ELE);
    if (new_base==NULL) {
	return berr_set(BITMAP_OUT_OF_MEMORY);
    }

    /* reset new entries */
    for(i=bitmap->Nbitmap;i<new_Nbitmap;i++) new_base[i] = (BASE_TYPE)0;
    bitmap->base = new_base;
    bitmap->Nbitmap = new_Nbitmap;
    bitmap->Nbits = new_Nbits;

    return 0;
    
}





int BitmapFree(Bitmap bitmap)
/*
 *  Return index of an unset bit in bitmap
 *  returns <0 is there isn't one
 */
{
    int index;   /* index of bitmap where there is a free bit */
    int elements; /* number of whole elements used */
    int bit, element;
    BASE_TYPE mask;

    if (bitmap == NULL) {
	return berr_set(BITMAP_INVALID_ARGUMENTS);
    }

    /*
     * find index in bitmap where there is a bit unset
     *
     * IMPORTANT NOTE:
     *     The last word in the bitmap must be handled as a special case.
     *     There may not be a multiple of BIT_ELE in the bitmap.
     */

    /*
     * We keep track of the first free bit number. This is set to zero
     * on creation, and updated when we clear and set bits. This is a
     * point for the first _possible_ location of a free bit, it may
     * have been used.
     */

    /*
     * First try - if first_free is still free, use it.
     */
    if (bitmap->first_free < bitmap->Nbits &&
	!BIT_CHK(bitmap, bitmap->first_free)) {
	return bitmap->first_free++;
    }

    /*
     * Second try - the first free bit is at the end of the bitmap.
     */
    if (bitmap->first_free >= bitmap->Nbits) {
	if (BitmapExtend(bitmap, bitmap->first_free+1))
	    return -1;
	return bitmap->first_free++;
    }

    /*
     * Third try - We'll search from first free onwards. There'll be
     * no free bit earlier than this one.
     */
    for (elements = ELE_NBITS(bitmap->Nbits)-1,
	 index=BIT_IDX(bitmap->first_free);
	 index<elements && bitmap->base[index]==~(BASE_TYPE)0;
	 index++);

    /*
     * If this is the last element, then check whether or not the element
     * can be found in the allocated portion.
     */
    if (index == elements) {
	mask = ((BASE_TYPE)1<<(bitmap->Nbits%BIT_ELE))-1;
	if (mask == 0) mask = (BASE_TYPE)-1;
	if ((bitmap->base[index] & mask) == mask) {
	    /*
	     * No free bits found in start of element, extend
	     */
	    if (BitmapExtend(bitmap, bitmap->first_free = bitmap->Nbits+1))
		return -1;
	    return bitmap->first_free-1;
	}
    }

    /*
     * We now know that the free bit is in this element. Find and return it.
     */
    element = bitmap->base[index];
    for (bit = 0; (element&1) && element; element>>=1, bit++);
    bitmap->first_free = index * BIT_ELE + bit;
    return bitmap->first_free++;
}



/*
 * pretty print a bitmap
 */
int BitmapPrint(FILE *fp, Bitmap bitmap)
{
    int i;
    int j;

    if (bitmap==NULL) {
	return berr_set(BITMAP_INVALID_ARGUMENTS);
    }

#define BITMAP_CHUNKS 16
    for (i=0;i<bitmap->Nbits;) {
	fprintf(fp,"%6d ",i);
	for (j=0;j<BITMAP_CHUNKS && i<bitmap->Nbits;j++,i++)
	    fprintf(fp,"%c",BIT_CHK(bitmap,i)?'1':'0');
	fprintf(fp,"\n");
    }

    return 0;

}


Bitmap
InitBooleanOp(Bitmap bitmap1,
	      Bitmap bitmap2)
{
    int nrecords;
    Bitmap bitmap;

    if (bitmap1 == NULL) {
	(void)berr_set(BITMAP_INVALID_ARGUMENTS);
	return NULL;
    }
    if (bitmap2 == NULL) {
	(void)berr_set(BITMAP_INVALID_ARGUMENTS);
	return NULL;
    }
    nrecords = BIT_NBITS(bitmap1);
    bitmap = BitmapCreate(nrecords);
 
    return bitmap;
}

Bitmap 
BitmapNOT(Bitmap bitmap1)
{
    Bitmap bitmap;
    int i;
    int nrecords;

    if (bitmap1 == NULL) {
	(void)berr_set(BITMAP_INVALID_ARGUMENTS);
	return NULL;
    }
    nrecords = BIT_NBITS(bitmap1);
    bitmap = BitmapCreate(nrecords);

    for(i = 0; i < bitmap->Nbitmap; i++) 
	bitmap->base[i] = ~bitmap1->base[i];

    return bitmap;
}

Bitmap 
BitmapAND(Bitmap bitmap1,
	  Bitmap bitmap2)
{
    Bitmap bitmap;
    int i;

    if (NULL == (bitmap = InitBooleanOp(bitmap1, bitmap2))){
	return NULL;
    }

    for(i = 0; i < bitmap->Nbitmap; i++) 
	bitmap->base[i] = bitmap1->base[i] & bitmap2->base[i];

    return bitmap;
}

Bitmap 
BitmapOR(Bitmap bitmap1,
	 Bitmap bitmap2)
{
    Bitmap bitmap;
    int i;

    if (NULL == (bitmap = InitBooleanOp(bitmap1, bitmap2))){
	return NULL;
    }
    for(i = 0; i < bitmap->Nbitmap; i++) 
	bitmap->base[i] = bitmap1->base[i] | bitmap2->base[i];

    return bitmap;
}

Bitmap 
BitmapXOR(Bitmap bitmap1,
	  Bitmap bitmap2)
{
    Bitmap bitmap;
    int i;

    if (NULL == (bitmap = InitBooleanOp(bitmap1, bitmap2))){
	return NULL;
    }

    /* initialise bit array */
    for(i = 0; i < bitmap->Nbitmap; i++) 
	bitmap->base[i] = bitmap1->base[i] ^ bitmap2->base[i];

    return bitmap;
}

/* 
 * find the nth bit that has been set
 */
int 
FindNBitSet(Bitmap bitmap,
	    int n) 
{
    int total, i, t;
    uint4 k;
    signed int l;
    unsigned char *cp = (unsigned char *)bitmap->base;

    for (total=0; total + (t=nb(cp)) < n; cp+=4, total+=t)
	;

    t = bitmap->base[i = ((int)(cp-(unsigned char *)bitmap->base))/CHR_ELE];
    for (l=-1, k=1; total < n; k<<=1, l++) {
	if (t & k)
	    total++;
    }
    return i*32+l;
}
