/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

/*
 * File: array.h
 * Version:
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */

#ifndef _ARRAY_H_
#define _ARRAY_H_


#include <stddef.h>
#include "xerror.h"

#define ARRAY_ERR_START 		200
#define ARRAY_NO_ERROR			0
#define ARRAY_FULL			(ARRAY_ERR_START + 0)
#define ARRAY_INVALID_ARGUMENTS		(ARRAY_ERR_START + 1)
#define ARRAY_OUT_OF_MEMORY		(ARRAY_ERR_START + 2)

#define aerr_set(e) ((e)?xerr_set((e), ArrayErrorString((e))):0)

typedef struct {
    int size;			/* element size */
    int dim;			/* allocated number of elements */
    int max;			/* elements accessed */
    void *base;			/* base address of array */
} ArrayStruct, *Array;



extern Array ArrayCreate(size_t size, int dim);

extern int ArrayExtend(Array a, int dim);

extern void *ArrayRef(Array a, int i);

extern int ArrayDestroy(Array a);

#define ArrayMax(a) ( (a)->max )

#define ArrayBase(t,a) ( (t *)((a)->base) )

/*
#define arr(t,a,n) \
    (*(t*)((a)->base + (a)->size*(n)))

#define arrp(t,a,n) \
    ((t*)((a)->base + (a)->size*(n)))
*/




#define arr(t,a,n) \
    ((t*)((a)->base))[n]

#define ARR(t,a,n) \
    (*((t*)ArrayRef((a),(n))))

#define arrp(t,a,n) \
    &((t*)((a)->base))[n]
#define ARRP(t,a,n) \
    ((t*)ArrayRef(a,n))

extern char *ArrayErrorString(int error);


#endif /*_ARRAY_H_*/
