/*
 * Copyright (c) Medical Research Council 2000. All rights reserved.
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
#include <basecall.hpp>



//--------------
// Constructors
//--------------

void BaseCall::Init()
{
   // Original data
   Call       = -1;
   Position   = -1;
   Amplitude  = -1;
   Confidence = -1;

   // New data
   PeakCall = -1;
   for( int n=0; n<4; n++ )
   {
      PeakPosition[n]  = -1;
      PeakAmplitude[n] = -1;
      PeakWidth[n]     = -1;
   }
}


BaseCall::BaseCall()
{
   Init();
}


BaseCall::BaseCall( int call, int pos, int conf )
{
    Init();
    Call       = call;
    Position   = pos;
    Confidence = conf;
}


BaseCall::BaseCall( int a, int c, int g, int t )
{
   Init();
   PeakAmplitude[0] = a;
   PeakAmplitude[1] = c;
   PeakAmplitude[2] = g;
   PeakAmplitude[3] = t;
}



//--------------
// As Character
//--------------

char BaseCall::AsCharacter() const
{
   const char Table[] = { '-', 'A', 'C', 'G', 'T' };
   assert(Call>-2);
   assert(Call<4);
   return Table[Call+1];
}



//------------
// Peak Stuff
//------------

int BaseCall::PeakCount() const
{
   int pc = 0;
   if( PeakAmplitude[0] >= 0 ) pc++;
   if( PeakAmplitude[1] >= 0 ) pc++;
   if( PeakAmplitude[2] >= 0 ) pc++;
   if( PeakAmplitude[3] >= 0 ) pc++;
   return pc;
}


bool BaseCall::PeakPresent() const
{
   if( (PeakAmplitude[0]>=0) || (PeakAmplitude[1]>=0) ||
       (PeakAmplitude[2]>=0) || (PeakAmplitude[3]>=0) )
      return true;
   return false;
}


int BaseCall::Rank( std::size_t n, bool bIndex ) const
{
   assert(n<4);


   // Setup our call array
   mutlib_pair_t c[4];
   c[0].a = PeakAmplitude[0];    // Amplitude in a, index in b
   c[0].b = 0;
   c[1].a = PeakAmplitude[1];
   c[1].b = 1;
   c[2].a = PeakAmplitude[2];
   c[2].b = 2;
   c[3].a = PeakAmplitude[3];
   c[3].b = 3;



   // Do a quicksort, highest value last
   if( c[0].a > c[1].a )          // 1st pair
      Swap( c[0], c[1] );
   if( c[2].a > c[3].a )          // 2nd pair
      Swap( c[2], c[3] );
   if( c[0].a > c[2].a )          // Tops
      Swap( c[0], c[2] );
   if( c[1].a > c[3].a )          // Bottoms
      Swap( c[1], c[3] );
   if( c[1].a > c[2].a )          // Middle
      Swap( c[1], c[2] );



   // Return nth highest value
   return bIndex ? c[n].b : c[n].a;
}


int BaseCall::PeakOfRank( int n )
{
   return Rank(n,false);
}


int BaseCall::PeakOfRankAsIndex( int n )
{
   return Rank(n,true);
}



//------
// Swap
//------

void BaseCall::Swap( mutlib_pair_t& a, mutlib_pair_t& b ) const
{
   mutlib_pair_t tmp = a;
   a = b;
   b = tmp;
}
