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


#ifndef _TRACEALIGN_PREPROCESS_HPP_
#define _TRACEALIGN_PREPROCESS_HPP_



#include <trace.hpp>         // For Trace object
#include <array.hpp>         // For SimpleArray object



class TraceAlignPreprocessor
{
 public:
   // Constructor/Destructor
   TraceAlignPreprocessor() { }
  ~TraceAlignPreprocessor() { }



 public:
   // Services
   void               Flush();
   NumericArray<int>& Envelope()                    { return m_oEnvelope; }
   double             BaseIntervalMean() const      { return m_nIntervalMean; }
   int                BaseIntervalMin() const       { return m_nIntervalMin; }
   int                BaseIntervalMax() const       { return m_nIntervalMax; }
   int                BaseIntervalMode() const      { return m_nIntervalMode; }
   double             BaseIntervalStdDev() const    { return m_nIntervalStdDev; }
   void               PreprocessTrace( Trace& t, bool bCacheBaseIntervalStatistics=false );



 private:
   // Data
   int               m_nIntervalMin;
   int               m_nIntervalMax;
   int               m_nIntervalMode;
   double            m_nIntervalMean;
   double            m_nIntervalStdDev;
   NumericArray<int> m_oEnvelope;         // Length = no of samples
};



#endif
