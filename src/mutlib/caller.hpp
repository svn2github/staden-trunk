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


#ifndef _MUTLIB_CALLER_HPP_
#define _MUTLIB_CALLER_HPP_



#include <matrix.hpp>
#include <trace.hpp>



/**
   Base class for confidence and base callers.
*/
class Caller
{
 public:
    // Destructor
    virtual ~Caller() { }



 protected:
   // Data
   typedef struct
   {
      int Index;
      int Position;
      int Amplitude;

   }call_t;



 protected:
   // Services
   void SortAscending( call_t data[4] );
   void Swap( call_t& a, call_t& b ) const  { call_t tmp=a; a=b; b=tmp; }
   int  LoadPeaks( SimpleMatrix<int>& Peak, int nPos, int nAmbiguityWindow, call_t data[4] );
};



#endif

