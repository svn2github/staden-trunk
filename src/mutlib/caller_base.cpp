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


#include <cassert>
#include <cmath>             // For std::log10()
#include <cstring>           // For std::memset()
#include <algorithm>         // For std::max()
#include <dnatable.hpp>
#include <caller_base.hpp>



void BaseCaller::Init()
{
    m_nPeakRatio  = 0.0;
    m_nConfidence = 0.0;
    for( int n=0; n<3; n++ )
        m_nCall[n] = '-';
    for( int n=0; n<2; n++ )
    {
        m_nPosition[n]  = -1;
        m_nAmplitude[n] =  0;
    }
}



void BaseCaller::MakeCall( Trace& Tr, SimpleMatrix<int>& Peak, int nPos, int nAmbiguityWindow )
{
    assert(nPos>=0);
    assert(nAmbiguityWindow>0);
    call_t Signal[4];


    // Initialisation
    Init();
    // m_nPosition[2] = nPos; jkb 25/06/2003. What should this be?


    // Search for peaks and load them in
    int peaks = LoadPeaks( Peak, nPos, nAmbiguityWindow, Signal );


    // Find biggest peaks position
    if( peaks > 0 )
    {
        int max_sig = 0;
        int max_amp = -1;
        for( int n=3; n>=0; n-- )
        {
            if( Signal[n].Position >= 0 )
            {
                if( Signal[n].Amplitude > max_amp )
                {
                    max_sig = n;
                    max_amp = Signal[n].Amplitude;
                }
            }
        }
        nPos = Signal[max_sig].Position;
    }


    // Load trace amplitudes for peakless bases
    for( int n=0; n<4; n++ )
    {
        if( Signal[n].Position < 0 )
            Signal[n].Amplitude = Tr[n][nPos];
    }


    // Sort the entire lot by amplitude
    SortAscending( Signal );


    // Basecall single peak
    DNATable Table;
    if( peaks == 1 )
    {

        for( int n=3; n>=0; n-- )
        {
            if( Signal[n].Position >= 0 )
            {
                m_nCall[0]      = Table.LookupBase( Signal[n].Index );
                m_nCall[1]      = m_nCall[0];
                m_nPosition[0]  = Signal[n].Position;
                m_nAmplitude[0] = Signal[n].Amplitude;
            }
        }
    }


    // Basecall multiple peaks
    else if( peaks >= 2 )
    {
        call_t highest_signal;
        int    highest_signal_n;
        highest_signal.Index = -1;
        for( int n=3; n>=0; n-- )
        {
            if( Signal[n].Position >= 0 )
            {
                if( highest_signal.Index < 0 )
                {
                    highest_signal   = Signal[n];
                    highest_signal_n = n;
                }
                else
                {
                    m_nCall[0]      = Table.LookupBase( highest_signal.Index, Signal[n].Index );
                    m_nCall[1]      = Table.LookupBase( highest_signal.Index );
                    m_nCall[2]      = Table.LookupBase( Signal[n].Index );
                    m_nPosition[0]  = highest_signal.Position;
                    m_nAmplitude[0] = highest_signal.Amplitude;
                    m_nPosition[1]  = Signal[n].Position;
                    m_nAmplitude[1] = Signal[n].Amplitude;
                }
            }
        }
    }


    // Compute confidence value, just SNR(db) = 20*log(S/N)
    double S = Signal[3].Amplitude;
    double N = Signal[2].Amplitude;
    if( N <= 0.0 )
        N = 1.0;
    m_nPeakRatio  = S / N;
    m_nConfidence = m_nPeakRatio
      ? 20.0 * std::log10( m_nPeakRatio )
      : 0;
}

