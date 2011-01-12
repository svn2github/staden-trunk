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
#include <cstring>                // For std::memset()
#include <mutlib.h>
#include <trace.hpp>
#include <tracealign_cache.hpp>



void TraceAlignDestroyResults( tracealign_t* ta )
{
   assert(ta != NULL);


   // Reset result code/string
   ta->ResultCode      = MUTLIB_RESULT_SUCCESS;
   ta->ResultString[0] = 0;


   // Delete alignments, trace object destructor disposes of the Read* pointers
   if( ta->Alignment[MUTLIB_INPUT].Trace )
   {
      Trace t1;
      t1.Wrap( ta->Alignment[MUTLIB_INPUT].Trace, true );
      std::memset( &ta->Alignment[MUTLIB_INPUT], 0, sizeof(mutlib_trace_t) );
   }
   if( ta->Alignment[MUTLIB_INPUT_REFERENCE].Trace )
   {
      Trace t2;
      t2.Wrap( ta->Alignment[MUTLIB_INPUT_REFERENCE].Trace, true );
      std::memset( &ta->Alignment[MUTLIB_INPUT_REFERENCE], 0, sizeof(mutlib_trace_t) );
   }
}



void TraceAlignDestroyCache( tracealign_t* ta )
{
   assert(ta != NULL);


   // Delete reference cache data
   TraceAlignCache* pCache = static_cast<TraceAlignCache*>( ta->Cache );
   if( pCache )
      delete pCache;
   ta->Cache = 0;
}

