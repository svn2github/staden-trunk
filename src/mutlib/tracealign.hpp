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


#ifndef _TRACEALIGN_HPP_
#define _TRACEALIGN_HPP_


#include <mutlib.h>
#include <array.hpp>                    // For Array<T> objects
#include <trace.hpp>                    // For Trace objects
#include <align.hpp>                    // For Alignment object
#include <muttag.hpp>                   // For MutTag objects


void            TraceAlignDestroyCache( tracealign_t* ta );
void            TraceAlignDestroyResults( tracealign_t* ta );
mutlib_result_t TraceAlignValidateInput( tracealign_t* ta );
void            TraceAlignInterpolate( char cPad, SimpleArray<char>& Envelope, Trace& Tin, int nClipL, Trace& Tout );
void            TraceAlignInsertBases( char cPad, SimpleArray<char>& Envelope, Trace& Tin, Trace& Tout, int nOverlap[2] );
void            TraceAlignQuantiseEnvelope( NumericArray<int>& e, SimpleArray<char>& qe, int nLevels, int nLower, int nUpper );


#endif
