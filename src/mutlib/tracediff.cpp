/*
 * Copyright (c) Medical Research Council 2002. All rights reserved.
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
#include <cstring>                    // For strcpy(), memset()
#include <cstdio>                     // For std::printf(), std::sprintf()
#include <new>                        // For std::bad_alloc
#include <algorithm>                  // For std::max()
#include <list.hpp>                   // For list template
#include <trace.hpp>                  // For Trace object
#include <muttag.hpp>                 // For MutTag object
#include <tagarray.hpp>               // For TagArray object
#include <tracealign.hpp>             // For TraceAlignValidateInput() prototype
#include <tracediff.hpp>              // For helpers
#include <tracediff_parameters.hpp>   // For TraceDiffParameter object



/* #define VERBOSE_DEBUG */



/**
   Initialises an empty tracediff_t structure.
*/
void TraceDiffInit( tracediff_t* td )
{
   assert(td != NULL);
   TraceDiffParameters Parameter;
   std::memset( td, 0, sizeof(tracediff_t) );
   td->ResultString    = new char[ 512 ];
   td->ResultString[0] = 0;
   for( int n=0; n<TRACEDIFF_PARAMETERS; n++ )
       td->Parameter[n] = Parameter[n].Default();
   TraceAlignInit( &td->Alignment );
   td->Initialised = 1;
}



/**
   Frees up all memory allocated by the trace difference algorithm
   in the tracediff_t structure.
*/
void TraceDiffDestroy( tracediff_t* td )
{
   assert(td != NULL);
   assert(td->Initialised);
   try
   {
      // Delete all data
      TraceAlignDestroy( &td->Alignment );
      TraceDiffDestroyResults( td );
      delete [] td->ResultString;
   }
   catch(...)
   {
      // Shouldn't happen, but we musn't throw exceptions outside dll boundary
      assert(0);
   }
}



/**
   Gets the value of the tracediff parameter 'p'.
*/
double TraceDiffGetParameter( tracediff_t* td, tracediff_parameter_t p )
{
   assert(td != NULL);
   assert(td->Initialised);
   assert(p<TRACEDIFF_PARAMETERS);
   return (p<TRACEDIFF_PARAMETERS) ? td->Parameter[p] : 0.0;
}



/**
   Sets the value of the tracediff parameter 'p'.
*/
void TraceDiffSetParameter( tracediff_t* td, tracediff_parameter_t p, double v )
{
   assert(td != NULL);
   assert(td->Initialised);
   assert(p<TRACEDIFF_PARAMETERS);
   if( p < TRACEDIFF_PARAMETERS )
       td->Parameter[p] = v;
}



/**
   Sets the specified reference trace. The clip points are in base units.
   The caller retains ownership of the Read structure.
*/
void TraceDiffSetReference( tracediff_t* td, Read* w, mutlib_strand_t d, int ql, int qr )
{
   assert(td != NULL);
   assert(td->Initialised);
   TraceAlignSetReference( &td->Alignment, d, w, ql, qr );
}



/**
   Sets the current input trace. The clip points are in base units.
   The caller retains ownership of the Read structure.
*/
void TraceDiffSetInput( tracediff_t* td, Read* i, mutlib_strand_t d, int ql, int qr )
{
   assert(td != NULL);
   assert(td->Initialised);
   TraceAlignSetInput( &td->Alignment, d, i, ql, qr );
}



/**
   Returns the result code generated by the Execute() function. This can be
   used to query the result later without having to save the value returned
   by Execute().
*/
mutlib_result_t TraceDiffGetResultCode( tracediff_t* td )
{
   assert(td != NULL);
   assert(td->Initialised);
   return td->ResultCode;
}



/**
   Returns a read-only user friendly error message corresponding to the
   last result code. It may contain more useful information such as a
   filename or the particular intput that was invalid. This can be displayed
   to the user if required.
*/
const char* TraceDiffGetResultString( tracediff_t* td )
{
   assert(td != NULL);
   assert(td->Initialised);
   return td->ResultString;
}



/**
   Returns the difference trace as a Read structure. Ownership of the Read
   structure is retained by the trace difference algorithm. If 'dl' or 'dr'
   is not null, the left and right clip points in base numbers are returned.
*/
Read* TraceDiffGetDifference( tracediff_t* td, int* dl, int* dr )
{
   assert(td != NULL);
   assert(td->Initialised);
   if( dl != 0 )
      *dl = td->DifferenceLeft;
   if( dr != 0 )
      *dr = td->DifferenceRight;
   return td->Difference;
}



/**
   Returns the number of tags generated by the mutation scanning algorithm
   which are available to be read.
*/
int TraceDiffGetTagCount( tracediff_t* td )
{
   assert(td != NULL);
   assert(td->Initialised);
   return td->TagCount;
}



/**
   Returns a pointer to the 'nth' tag in the list. If the tag item is empty,
   a null pointer is returned.
*/
mutlib_tag_t* TraceDiffGetTag( tracediff_t* td, int n )
{
   assert(td != NULL);
   assert(td->Initialised);
   assert(n<td->TagCount);
   if( n < td->TagCount )
   {
      assert(td->Tag != NULL);
      return &td->Tag[n];
   }
   return 0;
}



/**
   Executes the trace difference algorithm optionally followed by the
   mutation detection algorithm.
*/
mutlib_result_t TraceDiffExecute( tracediff_t* td, tracediff_algorithm_t a )
{
    enum { STATE_INITIALISE, STATE_VALIDATE_INPUT, STATE_TRACE_ALIGN,
           STATE_TRACE_DIFFERENCE, STATE_MUTATION_ANALYSIS, STATE_EXIT };
    int                 n;
    mutlib_strand_t     Strand = 0;
    TraceDiffParameters Parameter;
    Trace               AlignedTrace[2];
    int                 AlignedTraceClipL[2];
    int                 AlignedTraceClipR[2];
    List<MutTag>        DiffTagList;
    Trace*              DiffTrace = 0;
    bool                Done      = false;
    int                 State     = STATE_INITIALISE;
    assert(td != NULL);
    try
    {
        while(!Done)
        {
            switch(State)
            {
                case STATE_INITIALISE:
                    // Destroy old results
                    TraceDiffDestroyResults( td );
                    Strand = td->Alignment.Input.Strand;
                    State  = STATE_VALIDATE_INPUT;
                    break;



                case STATE_VALIDATE_INPUT:
                    // Check input values
                    State = STATE_EXIT;
                    for( n=0; n<TRACEDIFF_PARAMETERS; n++ )
                        Parameter[n].Value( td->Parameter[n] );
                    if( TraceDiffValidateParameters(td,Parameter) != MUTLIB_RESULT_SUCCESS )
                        break;
                    if( TraceAlignValidateInput(&td->Alignment) != MUTLIB_RESULT_SUCCESS )
                    {
                        td->ResultCode = td->Alignment.ResultCode;
                        std::strcpy( td->ResultString, td->Alignment.ResultString );
                        break;
                    }
                    State = STATE_TRACE_ALIGN;
                    break;



                case STATE_TRACE_ALIGN:
                    // Align the reference and input traces
                    if( TraceAlignExecute(&td->Alignment) != MUTLIB_RESULT_SUCCESS )
                    {
                        td->ResultCode = TraceAlignGetResultCode( &td->Alignment );
                        std::strcpy( td->ResultString, TraceAlignGetResultString(&td->Alignment) );
                        State = STATE_EXIT;
                        break;
                    }
                    for( int n=0; n<2; n++ )
                        AlignedTrace[n].Wrap(TraceAlignGetAlignment(&td->Alignment,static_cast<mutlib_input_t>(n),&AlignedTraceClipL[n],&AlignedTraceClipR[n]),false);
                    State = STATE_TRACE_DIFFERENCE;
                    break;



                case STATE_TRACE_DIFFERENCE: {
                    // Scale & Subtract the two traces to obtain the difference trace
                    State = (a&TRACEDIFF_ALGORITHM_DEFAULT_DIFFERENCE_ONLY) ? STATE_EXIT : STATE_MUTATION_ANALYSIS;
                    if( Parameter[TRACEDIFF_PARAMETER_YSCALE].Value() > 0.0 )
                        AlignedTrace[MUTLIB_INPUT].ScaleTo( AlignedTrace[MUTLIB_INPUT_REFERENCE] );
                    DiffTrace = AlignedTrace[MUTLIB_INPUT].Subtract( AlignedTrace[MUTLIB_INPUT_REFERENCE] );
                    if( !DiffTrace )
                        throw std::bad_alloc();
                    DiffTrace->AutoDestroy( false );
                    td->Difference      = DiffTrace->Raw();
                    td->DifferenceLeft  = AlignedTraceClipL[MUTLIB_INPUT];
                    td->DifferenceRight = AlignedTraceClipR[MUTLIB_INPUT];
                    break; }



                case STATE_MUTATION_ANALYSIS:
                    // Analyse the difference trace for possible mutations
                    TraceDiffScanForMutations( *DiffTrace, Strand, DiffTrace->IntervalMode(),
                        AlignedTraceClipL[MUTLIB_INPUT], Parameter, DiffTagList );
                    if( DiffTagList.Count() > 0 )
                    {
                        // Construct output mutation tag list
                        TagArray OutTags;
                        bool bComplementBases = Parameter[TRACEDIFF_PARAMETER_COMPLEMENT_TAGS].Value()>0.0 ? true : false;
                        OutTags.Create( DiffTagList.Count() );
                        OutTags.ReadTags( DiffTagList, 1, bComplementBases );
                        OutTags.AutoDestroy( false );
                        td->Tag        = OutTags.Raw();
                        td->TagCount   = DiffTagList.Count();
                    }
                    State = STATE_EXIT;
                    break;



                case STATE_EXIT:
                    // Quit
                    Done = true;
                    break;
            }
        }
    }
    catch( std::bad_alloc& )
    {
        // Memory allocation error
        td->ResultCode = MUTLIB_RESULT_OUT_OF_MEMORY;
        std::strcpy( td->ResultString, "Not enough memory available to complete the operation.\n" );
    }
    catch(...)
    {
        // Unknown exceptions
        td->ResultCode = MUTLIB_RESULT_UNEXPECTED_EXCEPTION;
        std::strcpy( td->ResultString, "An unexpected fatal exception has occurred, please "
                     "report the details to staden-package@mrc-lmb.cam.ac.uk.\n" );
    }



    // Cleanup & exit
    delete DiffTrace;
    return td->ResultCode;
}

