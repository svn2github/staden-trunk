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


#ifndef _MUTLIB_CALLER_SNR_HPP_
#define _MUTLIB_CALLER_SNR_HPP_



#include <trace.hpp>
#include <caller.hpp>



/**
   Given a trace, the SNR at the given position is calculated.
*/
class SNRCaller : public Caller
{
 public:
   // Constructors
   SNRCaller() : m_nPosition(-1), m_nSNR(0.0), m_nRatio(0.0)    { }
   SNRCaller( Trace& Tr, int nPos )                             { MakeCall(Tr,nPos); }



 public:
   // Services
   double GetSNRLog() const                                     { return m_nSNR; }
   double GetSNRLinear() const                                  { return m_nRatio; }
   int    GetPosition() const                                   { return m_nPosition; }
   void   MakeCall( Trace& Tr, int nPos );



 private:
    // Data
   double m_nSNR;
   double m_nRatio;
   int    m_nPosition;
};



#endif

