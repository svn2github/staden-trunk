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



#ifndef _MUTLIB_BASECALL_HPP_
#define _MUTLIB_BASECALL_HPP_


#include <cstddef>   // For size_t

// Useful structure
typedef struct
{
   int a;
   int b;

}mutlib_pair_t;



class BaseCall
{
/*
   Used to record and manipulate all base call information.
*/
 public:
   // Original data
   int Call;                     // 0=A, 1=C, 2=G, 3=T
   int Position;                 // Base position in samples
   int Amplitude;                // Size of peak
   int Confidence;               // In phred units

   // New Data
   int PeakCall;                 // 0=A, 1=C, 2=G, 3=T
   int PeakPosition[4];          // Centre points for each peak in samples
   int PeakAmplitude[4];         // 0=A, 1=C, 2=G, 3=T
   int PeakWidth[4];             // Useful measure of peak quality



 public:
   // Constructors
   BaseCall();
   BaseCall( int a, int c, int g, int t );
   BaseCall( int call, int pos, int conf=-1 );



 public:
   // Services
   char AsCharacter() const;
   int  PeakCount() const;
   bool PeakPresent() const;
   int  PeakOfRank( int n );
   int  PeakOfRankAsIndex( int n );



 private:
   // Helpers
   void Init();
   int  Rank( std::size_t n, bool bIndex ) const;
   void Swap( mutlib_pair_t& a, mutlib_pair_t& b ) const;
};



#endif

