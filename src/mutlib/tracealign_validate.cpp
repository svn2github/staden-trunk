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



#include <cstdio>           // For sprintf()
#include <validate.hpp>



/*
   Checks all the input data to ensure it's valid.
*/
mutlib_result_t TraceAlignValidateInput( tracealign_t* ta )
{
   const mutlib_strand_t f = MUTLIB_STRAND_FORWARD;
   const mutlib_strand_t r = MUTLIB_STRAND_REVERSE;



   // Check initialisation flag
   ta->ResultCode = MUTLIB_RESULT_INVALID_INPUT;
   if( !ta->Initialised )
   {
      std::sprintf( ta->ResultString, "Uninitialised input structure.\n" );
      return ta->ResultCode;
   }



   // Check input trace
   if( MutlibValidateTrace(ta->Input,ta->ResultString,"input") )
      return ta->ResultCode;



   // Validate input trace clip points
   if( MutlibValidateTraceClipPoints(ta->Input,ta->ResultString,"input") )
      return ta->ResultCode;



   // Check forward reference trace
   if( ta->Input.Strand == MUTLIB_STRAND_FORWARD )
   {

      // Check trace
      if( MutlibValidateTrace(ta->Reference[f],ta->ResultString,"reference") )
         return ta->ResultCode;


      // Validate trace clip points
      if( MutlibValidateTraceClipPoints(ta->Reference[f],ta->ResultString,"reference") )
         return ta->ResultCode;
   }



   // Check reverse reference trace
   if( ta->Input.Strand == MUTLIB_STRAND_REVERSE )
   {
      // Check trace
      if( MutlibValidateTrace(ta->Reference[r],ta->ResultString,"reference") )
         return ta->ResultCode;


      // Validate trace clip points
      if( MutlibValidateTraceClipPoints(ta->Reference[r],ta->ResultString, "reference") )
         return ta->ResultCode;
   }



   // Success!
   ta->ResultCode = MUTLIB_RESULT_SUCCESS;
   return MUTLIB_RESULT_SUCCESS;
}

