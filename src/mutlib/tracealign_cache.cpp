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
#include <tracealign_cache.hpp>



//-------
// Flush
//-------

void TraceAlignCache::Flush()
{
   RefData[0].Flush();
   RefData[1].Flush();
   AlignmentMatrix.Empty();
}



//-------------------------
// Create Alignment Matrix
//-------------------------

void TraceAlignCache::CreateAlignmentMatrix( int nMatrixSize, int nLevels, int nOffset )
{
   int r, c;
   assert(nLevels>0);
   assert(nOffset>=0);



   // Create matrix
   if( AlignmentMatrix.Raw() )
      AlignmentMatrix.Empty();
   AlignmentMatrix.Create( nMatrixSize, nMatrixSize );



   // Fill the score matrix in toeplitz form according to number of
   // characters and offset. For example, a 4-level matrix would look
   // like this (the zero character is not used because it looks like
   // a null terminator to the alignment routines):
   //
   //  0 0 0 0 0
   //  0 4 3 2 1
   //  0 3 4 3 2
   //  0 2 3 4 3
   //  0 1 2 3 4
   //
   // Anything along the diagonal, indicating a perfect match, gets the
   // highest score. The score falls off linearly as the amplitude 
   // difference increases.
   for( r=nOffset; r<nMatrixSize-1; r++ )
   {
      for( c=nOffset; c<nMatrixSize-1; c++ )
      {
         if( (r-c) <= 0 )
         {
            AlignmentMatrix[r][c] = (r-c) + nLevels;
         }
         else
         {
            AlignmentMatrix[r][c] = (c-r) + nLevels;
         }
      }
   }
}
