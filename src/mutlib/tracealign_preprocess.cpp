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
#include <new>				// For bad_alloc()
#include <array.hpp>
#include <basecall.hpp>
#include <tracealign_preprocess.hpp>


//-------
// Flush
//-------

void TraceAlignPreprocessor::Flush()
{
    m_oEnvelope.Empty();
}



//------------------
// Preprocess Trace
//------------------

void TraceAlignPreprocessor::PreprocessTrace( Trace& t, bool bCacheBaseIntervalStatistics )
{
   // Sort the basecalls
   t.Sort();



   // Initialisation
   m_nIntervalMin    = 0;
   m_nIntervalMax    = 0;
   m_nIntervalMode   = 0;
   m_nIntervalMean   = 0.0;
   m_nIntervalStdDev = 0.0;



   // Compute base interval statistics
   if( bCacheBaseIntervalStatistics )
   {
       m_nIntervalMin    = t.IntervalMin();
       m_nIntervalMax    = t.IntervalMax();
       m_nIntervalMode   = t.IntervalMode();
       m_nIntervalMean   = t.IntervalMean();
       m_nIntervalStdDev = t.IntervalStdDev();
   }



   // Create a trace envelope in the A channel
   Trace& Envelope = *t.CreateEnvelope();
   if( !(&Envelope) )
      throw std::bad_alloc();

   

   // Make copy of envelope
   m_oEnvelope.Empty();
   m_oEnvelope.Create( t.Samples() );
   for( int n=0; n<t.Samples(); n++ )
      m_oEnvelope[n] = Envelope[0][n];



   // Cleanup
   delete &Envelope;
}



void TraceAlignQuantiseEnvelope( NumericArray<int>& e, SimpleArray<char>& qe, int nLevels, int nLower, int nUpper )
{
/*
   Quantises the envelope 'e' into 'nLevels' between 'nLower' and 'nUpper'.
*/
   assert(nLevels>0);
   assert(nLower<nUpper);



   // Allocate storage for output
   qe.Empty();
   qe.Create( e.Range() );


   // Determine quanta, round up
   int nQuanta = nUpper / nLevels + 1;



   // Do quantisation
   int n, k;
   for( k=0, n=e.RangeLowerLimit(); n<=e.RangeUpperLimit(); k++, n++ )
      qe[k] = char( (e[n] / nQuanta) + nLower );
}
