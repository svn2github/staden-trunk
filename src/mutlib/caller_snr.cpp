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



#include <cmath>                       // For std::log10()
#include <caller_snr.hpp>



void SNRCaller::MakeCall( Trace& Tr, int nPos )
{
   call_t Signal[4];


   // Load in data
   m_nPosition         = nPos;
   Signal[0].Index     = 0;
   Signal[0].Amplitude = Tr[0][nPos];
   Signal[1].Index     = 1;
   Signal[1].Amplitude = Tr[1][nPos];
   Signal[2].Index     = 2;
   Signal[2].Amplitude = Tr[2][nPos];
   Signal[3].Index     = 3;
   Signal[3].Amplitude = Tr[3][nPos];



   // Sort into ascending order of amplitude
   SortAscending( Signal );



   // Compute SNR=confidence, based on 1st and 2nd highest trace levels
   // S/N (db) = 20*log(s/n).
   //
   double S = static_cast<double>( Signal[3].Amplitude );
   double N = static_cast<double>( Signal[2].Amplitude );
   if( N <= 0.0 )
       N = 1.0;
   m_nRatio = S / N;
   m_nSNR   = 20.0 * std::log10( m_nRatio );
}

