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


#include <dnatable.hpp>
#include <caller_level.hpp>


LevelCaller::LevelCaller()
{
    m_nPosition = -1;
    for( int n=0; n<4; n++ )
        m_cBase[n] = '-';
}



void LevelCaller::MakeCall( Trace& Tr, int nPos )
{
    DNATable Table;
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


    // Lookup base characters and store results
    for( int n=0; n<4; n++ )
    {
        m_cBase[n]      = Table.LookupBase( Signal[n].Index );
        m_nAmplitude[n] = Signal[n].Amplitude;
    }
}
