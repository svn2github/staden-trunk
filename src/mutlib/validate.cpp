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



#include <cstdio>         // For sprintf()
#include <mutlib.h>



/*
   Validates the trace described by the trace descriptor 'td'. If an error occurs
   the message is written out to the result string buffer 'rs'. The string 's' is
   used to customise the error message with trace type. eg: "reference" or "input".
*/
mutlib_result_t MutlibValidateTrace( mutlib_trace_t& td, char* rs, const char* s )
{
   const char* strand = (td.Strand==MUTLIB_STRAND_FORWARD) ? "forward" : "reverse";



   // Check for existence
   if( !td.Trace )
   {
      std::sprintf( rs, "Missing %s %s trace.\n", strand, s );
      return MUTLIB_RESULT_INVALID_INPUT;
   }



   // Check trace length
   if( td.Trace->NBases <= 0 )
   {
      std::sprintf( rs, "Zero length %s %s trace %s.\n", strand, s,
                    td.Trace->trace_name );
      return MUTLIB_RESULT_INVALID_INPUT;
   }
   return MUTLIB_RESULT_SUCCESS;
}


/*
   Validate a trace's clip points.
   Bases are numbered from 1..N. A left clip point at 1 removes the leftmost
   base, a right clip point at N removes the Nth base. Clip points of less than
   zero imply no clipping in that direction.
*/
mutlib_result_t MutlibValidateTraceClipPoints( mutlib_trace_t& td, char* rs, const char* s )
{
   const char* strand = (td.Strand==MUTLIB_STRAND_FORWARD) ? "Forward" : "Reverse";



   // Convert clip points for input trace, limit right clip point if required
   if( td.ClipL < 0 )
      td.ClipL = 0;
   if( td.ClipR < 0 )
      td.ClipR = td.Trace->NBases + 1;
   if( td.ClipR > (td.Trace->NBases+1) )
      td.ClipR = td.Trace->NBases + 1;



   // Check input trace clipping range. Here we check to ensure that the
   // clip range is sufficiently large to be useful - at least 10 bases.
   // Otherwise we will waste lots of time going through the algorithms
   // only to find out at the end that it was a pointless exercise.
   if( (td.ClipR - td.ClipL - 1) < 10 )
   {
      std::sprintf( rs, "%s %s trace clip range of (%d,%d) is too small in %s.\n",
                    strand, s, td.ClipL, td.ClipR, td.Trace->trace_name );
      return MUTLIB_RESULT_INVALID_INPUT;
   }
   return MUTLIB_RESULT_SUCCESS;
}

