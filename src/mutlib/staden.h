/*
 * Copyright (c) Medical Research Council 2001. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written as part of the Staden Package at the MRC Laboratory
 * of Molecular Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 *
 */


#ifndef _MUTLIB_STADEN_H_
#define _MUTLIB_STADEN_H_


#ifdef __cplusplus
extern "C" {
#endif


#include <Misc/misc.h>
#include <seq_utils/dna_utils.h>

/*
 * This is ugly, but io_lib's compress.h defines a pipe2()
 * function. It's only internally used and doesn't really need to be
 * there. However on more recent linuxes (>= 2.6.27) it conflicts with
 * a new pipe2 function added to unistd.h. In theory that's protected
 * behind __USE_GNU, but no matter what I try I cannot undef that.
 *
 * So instead we temporarily rename io_lib's pipe2 to something else.
 * Apologies for the shenanigans.
 */
#define pipe2 pipe_into
#include <io_lib/compress.h>
#undef pipe2

#include <io_lib/Read.h>


#ifdef __cplusplus
}
#endif



#endif

