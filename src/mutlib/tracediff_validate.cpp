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



#include <cstdio>                       // For sprintf()
#include <validate.hpp>
#include <tracediff_parameters.hpp>



/*
   Checks the input parameters to ensure they are valid.
*/
mutlib_result_t TraceDiffValidateParameters( tracediff_t* td, TraceDiffParameters& p )
{
   // Check parameters
   td->ResultCode = MUTLIB_RESULT_SUCCESS;
   for( int n=0; n<TRACEDIFF_PARAMETERS; n++ )
   {
        if( p[n].IsValid() == false )
        {
            std::sprintf( td->ResultString, "Invalid %s parameter %.2f. "
                          "Must be in the range %.2f-%.2f.\n",
                          p[n].Name(), p[n].Value(), p[n].Minimum(), p[n].Maximum() );
            td->ResultCode = MUTLIB_RESULT_INVALID_INPUT;
            return MUTLIB_RESULT_INVALID_INPUT;
        }
   }
   return MUTLIB_RESULT_SUCCESS;
}
