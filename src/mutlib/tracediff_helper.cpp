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
#include <trace.hpp>
#include <tracealign.hpp>



void TraceDiffDestroyResults( tracediff_t* td )
{
   // Reset result code/string
   td->ResultCode      = MUTLIB_RESULT_SUCCESS;
   td->ResultString[0] = 0;


   // Delete tags
   for( int n=0; n<td->TagCount; n++ )
       delete [] td->Tag[n].Comment;
   delete [] td->Tag;
   td->Tag      = 0;
   td->TagCount = 0;


   // Delete read structure
   if( td->Difference )
   {
       Trace Diff;
       Diff.Wrap( td->Difference, true );
       td->Difference = 0;
   }
}
