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


#ifndef _MUTLIB_CALLER_BASE_HPP_
#define _MUTLIB_CALLER_BASE_HPP_


#include <cassert>
#include <matrix.hpp>
#include <trace.hpp>
#include <caller.hpp>


/**
   Given a peaks matrix containing peak positions and amplitudes, the
   peaks within +/- nAmbiguityWindow columns around the central position
   are ranked in order and a basecall is made.
*/
class BaseCaller : public Caller
{
 public:
   // Constructors
   BaseCaller() { Init(); }
   BaseCaller( Trace& Tr, SimpleMatrix<int>& Peak, int nPos, int nAmbiguityWindow ) { MakeCall(Tr,Peak,nPos,nAmbiguityWindow); }



 public:
   // Services
   bool   IsValid()                            { return m_nCall[0]!='-'; }
   void   Invalidate()                         { m_nCall[0]='-'; }
   char   GetBase( int n=0 ) const             { assert(n>=0); assert(n<3); return m_nCall[n]; }
   int    GetPosition( int n=0 ) const         { assert((n==0)||(n==1)); return m_nPosition[n]; }
   int    GetAmplitude( int n=0 ) const        { assert((n==0)||(n==1)); return m_nAmplitude[n]; }
   double GetPeakRatio() const                 { return m_nPeakRatio; }
   double GetConfidence() const                { return m_nConfidence; }
   void   MakeCall( Trace& Tr, SimpleMatrix<int>& Peak, int nPos, int nAmbiguityWindow );



 private:
    // Methods
    void Init();



 private:
   // Data
   char   m_nCall[3];
   double m_nPeakRatio;
   double m_nConfidence;
   int    m_nPosition[2];
   int    m_nAmplitude[2];
};



#endif

