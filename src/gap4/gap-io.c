/*
 * File:
 *
 * Author: Simon Dear
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
#include <string.h>

#include "g-os.h"
#include "g-error.h"
#include "g-struct.h"

#include "gap-error.h"
#include "gap-if.h"
#include "xalloc.h"


/*************************************************************
 * Low level IO routines
 *************************************************************/

static void swap2(int2 *src, int2 *dst, int N)
{
    for(;N;N--,src++,dst++)
	swap_int2(*src,*dst);
}

static void swap4(int4 *src, int4 *dst, int N)
{
    for(;N;N--,src++,dst++)
	swap_int4(*src,*dst);
}

static void swap2_inline(int2 *src, int N)
{
    for(;N;N--,src++)
	swap_int2(*src,*src);
}

static void swap4_inline(int4 *src, int N) 
{
    for(;N;N--,src++)
	swap_int4(*src,*src);
}


static void *BUFFER = NULL;
static int BUFFERL = 0;


static void set_buffer(int len)
{
    if (len > BUFFERL) {
	BUFFER = (void *)xrealloc(BUFFER,len);
	BUFFERL = len;
	if (BUFFER == 0) {
	    fprintf(stderr, "FATAL: Couldn't alloc %d bytes\n", len);
	}
    }
}


int (* GAP_READ) (GapClient *s, GView v, void *buf, int len, GCardinal type_check, int size) = NULL;
int (* GAP_WRITE) (GapClient *s, GView v, void *buf, int len, GCardinal type_check, int size) = NULL;



static int GAP_READ_NO_SWAP(GapClient *s, GView v, void *buf, int len, GCardinal type_check, int size)
/*
 * Read in a GAP database record
 */
{
    GCardinal type;
    GCardinal Nvec;
    GIOVec vec[2];
    int err;

    vec[0].buf = &type;
    vec[0].len = sizeof(type);
    Nvec = 1;
    if (len) {
	vec[1].buf = buf;
	vec[1].len = len;
	Nvec++;
    }

    err = g_readv(s, v, vec, Nvec);
    if (err)
	return err;

    /* check we have the correct type */
    if (type != type_check)
	return gaperr_set(GAPERR_INVALID_TYPE);

    return 0;
}


static int GAP_WRITE_NO_SWAP(GapClient *s, GView v, void *buf, int len, GCardinal type, int size)
/*
 * Write a GAP database record
 */
{
    GCardinal Nvec;
    GIOVec vec[2];

    vec[0].buf = &type;
    vec[0].len = sizeof(type);
    Nvec = 1;
    if (len) {
	vec[1].buf = buf;
	vec[1].len = len;
	Nvec++;
    }

    return g_writev(s, v, vec, Nvec);
}




static int GAP_READ_SWAP(GapClient *s, GView v, void *buf, int len, GCardinal type_check, int size)
/*
 * Read in a GAP database record
 */
{
    GCardinal type;
    GCardinal Nvec;
    GIOVec vec[2];
    int err;

    vec[0].buf = &type;
    vec[0].len = sizeof(type);
    Nvec = 1;
    if (len) {
	vec[1].buf = buf;
	vec[1].len = len;
	Nvec++;
    }

    err = g_readv(s, v, vec, Nvec);
    if (err)
	return err;

    /*SWAP*/
    swap_int4(type, type);
    /*SWAP*/
    if (len) {
	switch(size) {
	case 4:
	    swap4_inline(buf,len/size);
	    break;
	case 2:
	    swap2_inline(buf,len/size);
	    break;
	default:
	    break;
	}
    }

    /* check we have the correct type */
    if (type != type_check)
	return gaperr_set(GAPERR_INVALID_TYPE);

    return 0;
}


static int GAP_WRITE_SWAP(GapClient *s, GView v, void *buf, int len, GCardinal type_check, int size)
/*
 * Write a GAP database record
 */
{
    GCardinal Nvec;
    GIOVec vec[2];

    /*SWAP*/
    swap_int4(type_check, type_check);

    vec[0].buf = &type_check;
    vec[0].len = sizeof(type_check);
    Nvec = 1;
    if (len) {
	/*SWAP*/
	set_buffer(len);
	switch(size) {
	case 4:
	    swap4(buf,BUFFER,len/size);
	    break;
	case 2:
	    swap2(buf,BUFFER,len/size);
	    break;
	default:
	    memcpy(BUFFER,buf,len);
	    break;
	}

	vec[1].buf = BUFFER;
	vec[1].len = len;
	Nvec++;
    }

    return g_writev(s, v, vec, Nvec);
}






void gap_io_init(void)
{
    int i=1;


    if ( *(char *)&i ) {
	GAP_READ = GAP_READ_SWAP;
	GAP_WRITE = GAP_WRITE_SWAP;
	BUFFERL = 512;
	BUFFER = (void *)xmalloc(BUFFERL);
    } else {
	GAP_READ = GAP_READ_NO_SWAP;
	GAP_WRITE = GAP_WRITE_NO_SWAP;
    }
}
