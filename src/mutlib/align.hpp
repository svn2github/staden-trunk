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



#ifndef _MUTLIB_ALIGN_HPP_
#define _MUTLIB_ALIGN_HPP_



#include <cassert>
#include <sp_alignment_structs.h>    // For alignment structures
#include <sp_alignment.h>            // For alignment constants
#include <matrix.hpp>                // For SimpleMatrix<T> object



class Alignment
{
 public:
   // Algorithms
   typedef enum
   {
      ALGORITHM_NORMAL      = SP_ALIGNMENT_ALGORITHM_A,
      ALGORITHM_BLOCK_ALIGN = SP_ALIGNMENT_ALGORITHM_B,
      ALGORITHM_POISSON     = SP_ALIGNMENT_ALGORITHM_C

   }algorithm_t;
   // Edge scoring
   typedef enum
   {
      EDGE_SCORE_LEFT  = SP_ALIGNMENT_LEFT_EDGE_GAPS_COUNT,
      EDGE_SCORE_RIGHT = SP_ALIGNMENT_BEST_RIGHT_EDGE

   }edge_score_t;
   // Support for multiple alignment
   enum
   {
      MAX_INPUT_SEQUENCES = 2
   };



 public:
   // Constructors/destructor
   Alignment();
  ~Alignment();



 public:
   // Services
   void   BandSize( int n )                               { assert(n>=0); m_nBand=n; }
   void   PadSymbol( int p )                              { assert(p);    m_nPadSymbol=p; }
   int    PadSymbol() const                               { return m_nPadSymbol; }
   void   GapPenalty( int nBegin, int nExtend )           { m_nGapPenaltyBegin=nBegin; m_nGapPenaltyExtend=nExtend; }
   void   EdgeScore( edge_score_t es );
   void   InputSequence( int n, const char* s, int l=-1 );
   void   Matrix( int** m, int n, bool AutoDestroy=true );
   int    Execute( algorithm_t a );
   double OutputScore() const;
   char*  OutputSequence( int n ) const;
   int    OutputSequenceLength( int n ) const;
   int    OutputSequenceLeftOverlap( int n ) const;
   int    OutputSequenceRightOverlap( int n ) const;
   void   DumpToFile( const char* s, bool AsInteger=false ) const;



 private:
   // Helpers
   void CreateDefaultMatrix();



 private:
   // Data
   int               m_nBand;
   sp::ALIGN_PARAMS* m_pParams;
   sp::OVERLAP*      m_pOverlap;
   SimpleMatrix<int> m_oMatrix;
   int               m_nPadSymbol;
   int               m_nEdgeScore;
   int               m_nGapPenaltyBegin;
   int               m_nGapPenaltyExtend;
   static bool       m_bDNALookupInitialised;
   const char*       m_pInputSequence[MAX_INPUT_SEQUENCES];
   int               m_nInputSequenceLength[MAX_INPUT_SEQUENCES];
};



#endif

