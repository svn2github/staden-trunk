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
 * File: array.c
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

#include <stdio.h>
#include <sys/types.h>

#include "array.h"
#include "xalloc.h"

/*#include <c_asm.h>*/

char *ArrayErrorString(int err)
{
    switch(err) {
    case ARRAY_NO_ERROR:          return "No error";
    case ARRAY_FULL:     	  return "Array full";
    case ARRAY_INVALID_ARGUMENTS: return "Invalid arguments";
    case ARRAY_OUT_OF_MEMORY:     return "Out of memory";
    default:			  return "Unknown error";
    }
}

Array ArrayCreate(size_t size, size_t dim)
/*
 * create a new array
 */
{
    Array a;

/*    printf("CR: 0x%lx:",(long)asm("bis %ra,%ra,%v0")); */

    if ( (a = (Array) xmalloc(sizeof(ArrayStruct)) ) == NULL ) {
	(void)aerr_set(ARRAY_OUT_OF_MEMORY);
    } else {
	a->size = size;
	a->dim = dim?dim:1;
	a->max = 0;
	if ( (a->base = (void *)xmalloc(a->size * a->dim)) == NULL ) {
	    (void)aerr_set(ARRAY_OUT_OF_MEMORY);
	    xfree(a);
	    a = NULL;
	}
    }
    
/*    printf("0x%lx - size %d\n", (long)a, a->size * a->dim); */
    return a;
    
}



int ArrayExtend(Array a, size_t dim)
/*
 * extend array
 */
{
    void *newbase;
    size_t old_dim;

    if (a == NULL) {
	return aerr_set(ARRAY_INVALID_ARGUMENTS);
    }

    if (dim < a->dim) {
	return 0;
    }

    old_dim = a->dim;
    while (dim >= a->dim) {
	a->dim = a->dim * 1.2 + 1;
    }

    if ( (newbase = (void *)xrealloc(a->base, a->size * a->dim)) == NULL ) {
	a->dim = old_dim;
	return aerr_set(ARRAY_OUT_OF_MEMORY);
    } else {
	a->base = newbase;
    }

    return 0;
}


void *ArrayRef(Array a, size_t i)
{
    if (a==NULL) {
	(void)aerr_set(ARRAY_INVALID_ARGUMENTS);
	return NULL;
    }

    if (i >= a->max) {
	if (i >= a->dim) {
	    if (ArrayExtend(a,i+1)) {
		/* ArrayExtend sets aerrnum */
		return NULL;
	    }
	}
	a->max = i+1;
    }

    return (void *) arrp(char,a,i*a->size);
}

int ArrayDestroy(Array a)
/*
 * destroy array
 */
{
/*    printf("DR: 0x%lx:0x%lx\n",(long)asm("bis %ra,%ra,%v0"), (long)a); */

    if (a==NULL) {
	return aerr_set(ARRAY_INVALID_ARGUMENTS);
    }
    
    if (a->base != NULL) xfree(a->base);
    a->base= NULL;
    xfree(a);

    return 0;
}


