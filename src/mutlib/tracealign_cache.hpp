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


#ifndef _TRACEALIGN_CACHE_HPP_
#define _TRACEALIGN_CACHE_HPP_



#include <matrix.hpp>
#include <tracealign_preprocess.hpp>



class TraceAlignCache
{
 public:
   // Constructor/Destructor
   TraceAlignCache()  { }
  ~TraceAlignCache()  { }



 public:
   // Services
   void Flush();
   void CreateAlignmentMatrix( int nMatrixSize, int nLevels, int nOffset );



 public:
   // Cached data
   TraceAlignPreprocessor RefData[2];
   SimpleMatrix<int>      AlignmentMatrix;
};



#endif

