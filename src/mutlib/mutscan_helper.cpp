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


#include <cassert>
#include <mutlib.h>



void MutScanDestroyResults( mutscan_t* ms )
{
   assert(ms != NULL);


   // Delete result string
   delete [] ms->ResultString;
   ms->ResultString = 0;
   ms->ResultCode   = MUTLIB_RESULT_SUCCESS;


   // Delete tags
   for( int n=0; n<ms->TagCount; n++ )
       delete [] ms->Tag[n].Comment;
   delete [] ms->Tag;
   ms->Tag      = 0;
   ms->TagCount = 0;
}

