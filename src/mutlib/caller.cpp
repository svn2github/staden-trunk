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


#include <caller.hpp>



int Caller::LoadPeaks( SimpleMatrix<int>& Peak, int nPos, int nAmbiguityWindow, call_t data[4] )
{
    assert(data != NULL);
    assert(nPos>=0);
    assert(nAmbiguityWindow>0);
    int peaks = 0;


    // Load in default data
    for( int n=0; n<4; n++ )
    {
        data[n].Index     =  n;
        data[n].Position  = -1;
        data[n].Amplitude = Peak[n][nPos];
    }


    // For each base that has no peak at nPos, search within ambiguity window
    // for another very close peak. Work outwards from middle position.
    const int rows = 4;
    const int cols = Peak.Cols();
    for( int r=0; r<rows; r++ )
    {
        // If base has a peak at nPos, record the position
        if( data[r].Amplitude )
        {
            data[r].Position = nPos;
            peaks++;
            continue;
        }


        // No peak, search for one nearby, within window
        int j   = nPos;
        int off = 1;
        int win = nAmbiguityWindow;
        while( (win>0) && ((j-off)>=0) && ((j+off)<cols) )
        {
            if( Peak[r][j-off] > 0 )
            {
                data[r].Position  = j-off;
                data[r].Amplitude = Peak[r][j-off];
                peaks++;
                break;
            }
            if( Peak[r][j+off] > 0 )
            {
                data[r].Position  = j+off;
                data[r].Amplitude = Peak[r][j+off];
                peaks++;
                break;
            }
            off++;
            win--;
        }
    }
    return peaks;
}


void Caller::SortAscending( call_t data[4] )
{
    // Sort into amplitude order, lowest to highest

    if( data[0].Amplitude > data[1].Amplitude )          // 1st pair
        Swap( data[0], data[1] );
    if( data[2].Amplitude > data[3].Amplitude )          // 2nd pair
        Swap( data[2], data[3] );
    if( data[0].Amplitude > data[2].Amplitude )          // Tops
        Swap( data[0], data[2] );
    if( data[1].Amplitude > data[3].Amplitude )          // Bottoms
        Swap( data[1], data[3] );
    if( data[1].Amplitude > data[2].Amplitude )          // Middle
        Swap( data[1], data[2] );
}

