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



#include <cstdio>			// For sprintf()
#include <validate.hpp>
#include <mutscan_parameters.hpp>



/*
   Checks all the input data to ensure it's valid.
*/
mutlib_result_t MutScanValidateInput( mutscan_t* ms, MutScanParameters& p )
{
   const mutlib_strand_t f = MUTLIB_STRAND_FORWARD;
   const mutlib_strand_t r = MUTLIB_STRAND_REVERSE;



   // Check initialisation flag
   ms->ResultCode = MUTLIB_RESULT_INVALID_INPUT;
   if( !ms->Initialised )
   {
      std::strcpy( ms->ResultString, "Uninitialised input structure.\n" );
      return ms->ResultCode;
   }



   // Check parameters
   for( int n=0; n<MUTSCAN_PARAMETERS; n++ )
   {
        if( p[n].IsValid() == false )
        {
            std::sprintf( ms->ResultString, "Invalid %s parameter %.2f. "
                          "Must be in the range %.2f-%.2f.\n",
                          p[n].Name(), p[n].Value(), p[n].Minimum(), p[n].Maximum() );
            return ms->ResultCode;
        }
   }



   // Check input trace
   if( MutlibValidateTrace(ms->InputTrace,ms->ResultString,"input") )
      return ms->ResultCode;



   // Validate input trace clip points
   if( MutlibValidateTraceClipPoints(ms->InputTrace,ms->ResultString,"input") )
      return ms->ResultCode;



   // Check forward reference trace
   if( ms->InputTrace.Strand == MUTLIB_STRAND_FORWARD )
   {

      // Check trace
      if( MutlibValidateTrace(ms->ReferenceTrace[f],ms->ResultString,"reference") )
         return ms->ResultCode;


      // Validate trace clip points
      if( MutlibValidateTraceClipPoints(ms->ReferenceTrace[f],ms->ResultString,"reference") )
         return ms->ResultCode;
   }



   // Check reverse reference trace
   if( ms->InputTrace.Strand == MUTLIB_STRAND_REVERSE )
   {
      // Check trace
      if( MutlibValidateTrace(ms->ReferenceTrace[r],ms->ResultString,"reference") )
         return ms->ResultCode;


      // Validate trace clip points
      if( MutlibValidateTraceClipPoints(ms->ReferenceTrace[r],ms->ResultString, "reference") )
         return ms->ResultCode;
   }



   // Success!
   ms->ResultCode = MUTLIB_RESULT_SUCCESS;
   return MUTLIB_RESULT_SUCCESS;
}

