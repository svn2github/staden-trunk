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


#ifndef _TRACEDIFF_HPP_
#define _TRACEDIFF_HPP_


#include <mutlib.h>
#include <array.hpp>                    // For Array<T> objects
#include <trace.hpp>                    // For Trace objects
#include <align.hpp>                    // For Alignment object
#include <muttag.hpp>                   // For MutTag objects
#include <tracediff_parameters.hpp>     // For TraceDiffParameters object


void            TraceDiffDestroyCache( tracediff_t* td );
void            TraceDiffDestroyResults( tracediff_t* td );
mutlib_result_t TraceDiffValidateParameters( tracediff_t* td, TraceDiffParameters& p );
int             TraceDiffFindSeqCutSite( const char* w, const char* i, int nLength, int nPos, int nWindow );
void            TraceDiffInterpolate( char cPad, SimpleArray<char>& Envelope, Trace& Tin, int nClipL, Trace& Tout );
void            TraceDiffInsertBases( char cPad, SimpleArray<char>& Envelope, Trace& Tin, Trace& Tout, int nOverlap[2] );
void            TraceDiffQuantiseEnvelope( NumericArray<int>& e, SimpleArray<char>& qe, int nLevels, int nLower, int nUpper );
void            TraceDiffScanForMutations( Trace& DiffTrace, mutlib_strand_t nStrand, int nBaseInterval, int nFirstBase, TraceDiffParameters& Parameter, List<MutTag>& Mutation );
void            TraceDiffResampleTraces( Trace& win, Trace& iin, Trace& wout, Trace& iout, NumericArray<int>& wpks, NumericArray<int>& ipks, NumericArray<int>& widx, NumericArray<int>& iidx );
int             TraceDiffFindChunkJoinPoint( SimpleArray<char>& WildChunkThis,  SimpleArray<char>& WildChunkPrev, SimpleArray<char>& InputChunkThis, SimpleArray<char>& InputChunkPrev, int nChunkOverlap, int nMatchCount, int& nJoinThis, int& nJoinPrev );


#endif
