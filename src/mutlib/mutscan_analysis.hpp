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


#ifndef _MUTSCAN_ANALYSIS_HPP_
#define _MUTSCAN_ANALYSIS_HPP_



#include <mutlib.h>                // For mutlib types
#include <array.hpp>               // For SimpleArray object
#include <matrix.hpp>              // For SimpleMatrix object
#include <trace.hpp>               // For Trace object
#include <list.hpp>                // For List object
#include <mutationtag.hpp>         // For MutationTag object
#include <mutscan_preprocess.hpp>  // For MutScanPreprocessor object



/**
   Analyses the preprocessed input data and traces for mutations.
*/
class MutScanAnalyser
{
 public:
   // Constructor/Destructor
   MutScanAnalyser() { }
  ~MutScanAnalyser() { }



 public:
   // Data
   SimpleMatrix<int>    Map;                  // 8*max(PEAKS) peak alignment maps, r[0,1]=A, r[2,3]=C, r[4,5]=G, r[6,7]=T
   int                  MapCount[4];          // 4*1 pair counts for each pair of Map rows
   SimpleMatrix<double> ScaleFactor;          // 4*max(PEAKS) peak pair scale factors, row 0=A, 1=C, 2=G, 3=T
   double               ScaleFactorMean[4];   // 4*1 vector containing mean value of scale factors
   double               ScaleFactorStdDev[4]; // 4*1 vector containing the std deviation of the scale factors
   List<MutationTag>    MutationTagList;



 public:
   // Services
   mutlib_result_t Execute( mutscan_t* ms, MutScanPreprocessor Data[], Trace Tr[], Trace* DifferenceTrace );



 private:
   // Helpers
   bool HasReferencePeak( int base, int pos );
   void AnalysePotentialMutations( Trace Tr[] );
   void AlignPeaks( MutScanPreprocessor Data[] );
   void AllocatePeakMap( MutScanPreprocessor Data[] );
   void ComputeScaleFactors( MutScanPreprocessor Data[] );
   void ValidateMutationsAgainstDifference( Trace& DiffTrace );
   void ComputeScaleFactorLimits( int r, double nosd, double data[2] );
   void DoLevelCall( int pos, Trace& Tr, MutationTag& Tag, bool AllowAmbiguity );
   void ScanForPotentialMutations( MutScanPreprocessor Data[], mutlib_strand_t Strand, Trace Tr[] );



 private:
   // Data
   int    m_nSearchWindow;
   double m_nHetSNRThreshold;
   double m_nUpperPeakDropThreshold;
   double m_nLowerPeakDropThreshold;
};



#endif
