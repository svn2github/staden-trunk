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


#ifndef _MUTLIB_CALLER_LEVEL_HPP_
#define _MUTLIB_CALLER_LEVEL_HPP_



#include <cassert>
#include <trace.hpp>
#include <caller.hpp>



/**
   Given a trace, the amplitude levels are ranked in ascending order.
   GetBase(0) returns the lowest amplitude base and similarly, GetBase(3)
   returns the highest level base.
*/
class LevelCaller : public Caller
{
 public:
   // Constructors
   LevelCaller();
   LevelCaller( Trace& Tr, int nPos )                   { MakeCall(Tr,nPos); }



 public:
   // Services
   char GetBase( int n ) const                          { assert((n>=0)&&(n<4)); return m_cBase[n]; }
   int  GetPosition() const                             { return m_nPosition; }
   int  GetAmplitude( int n ) const                     { assert((n>=0)&&(n<4)); return m_nAmplitude[n]; }
   void MakeCall( Trace& Tr, int nPos );


 private:
    // Data
   char m_cBase[4];
   int  m_nPosition;
   int  m_nAmplitude[4];
};



#endif

