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


#ifndef _MUTSCAN_PREPROCESS_HPP_
#define _MUTSCAN_PREPROCESS_HPP_



#include <mutlib.h>           // For mutlib types
#include <array.hpp>          // For SimpleArray object
#include <matrix.hpp>         // For SimpleMatrix object
#include <trace.hpp>          // For Trace object



/**
   Preprocesses the trace ready for further analysis by mutscan. One preprocessor
   object is used for each trace.
*/
class MutScanPreprocessor
{
 public:
   // Constructor/Destructor
   MutScanPreprocessor() { }
  ~MutScanPreprocessor() { }



 public:
   // Data
   SimpleMatrix<int> Peak;           // 5*SAMPLES peak matrix, value=amplitude, last row is logical OR of all peaks for every base
   NumericArray<int> Noise;          // 1*SAMPLES noise vector
   NumericArray<int> PeakCount;      // 4*1 peak count for each row of Peak
   int               PeakCountMax;   // Maximum number of peaks for all 4-bases
   double            SearchWindow;   // Search window in samples, derrived from trace stats
   int               PeakInterval;   // Peak interval mode



 public:
   // Services
   mutlib_result_t Execute( mutscan_t* ms, Trace& t, int n, int left, int right );



 private:
   // Helpers
   void Init();
   void CountPeaks();
   void PeakClip();
   void PeakSpacing();
   void EstimateNoiseFloor( Trace& t, int n );
   void PeakFind( Trace& t, int left, int right );



 private:
   // Data
   double m_nNoiseThreshold[2];
};



#endif
