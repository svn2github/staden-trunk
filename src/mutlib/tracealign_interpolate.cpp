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
#include <array.hpp>
#include <trace.hpp>
#include <align.hpp>



void TraceAlignInterpolate( char cPad, SimpleArray<char>& Envelope, Trace& Tin, int nClipL, Trace& Tout )
{
/*
    This routine copies the input trace amplitudes to the output trace, interpolating
    where indicated by pads in the envelope alignment sequence. Only the region of the
    trace Tin from n1 to the length of the aligned sequence is copied to Tout.
*/
    assert(Envelope.Length()==Tout.Samples());
    int d, s, k, p, x;
    double m[4];
    int    c[4];



    // Any initial padding is copied over with a trace amplitude of zero
    for( d=0; (d<Envelope.Length())&&(Envelope[d]==cPad); d++ )
    {
        Tout[0][d] = 0;
        Tout[1][d] = 0;
        Tout[2][d] = 0;
        Tout[3][d] = 0;
    }



    // Copy and interpolate
    for( d=d, s=nClipL; d<Envelope.Length(); d++ )
    {
        // If current character is a pad
        if( Envelope[d] == cPad )
        {
            // Determine pad run length
            for( k=d; (k<Envelope.Length())&&(Envelope[k]==cPad); k++ );
            p = k - d;


            // Determine equation of line y=mx+c parameters
            // Slope = (LastNonPad - FirstNonPad) / (Pads+1)
            c[0] = int(Tin[0][s-1]);
            c[1] = int(Tin[1][s-1]);
            c[2] = int(Tin[2][s-1]);
            c[3] = int(Tin[3][s-1]);
            m[0] = double(int(Tin[0][s]) - int(Tin[0][s-1])) / double(p+1);
            m[1] = double(int(Tin[1][s]) - int(Tin[1][s-1])) / double(p+1);
            m[2] = double(int(Tin[2][s]) - int(Tin[2][s-1])) / double(p+1);
            m[3] = double(int(Tin[3][s]) - int(Tin[3][s-1])) / double(p+1);


            // Do linear interpolation over padded region
            for( x=1; x<=p; x++ )
            {
                Tout[0][d+(x-1)] = TRACE( m[0]*x + c[0] );
                Tout[1][d+(x-1)] = TRACE( m[1]*x + c[1] );
                Tout[2][d+(x-1)] = TRACE( m[2]*x + c[2] );
                Tout[3][d+(x-1)] = TRACE( m[3]*x + c[3] );
            }
            d += p - 1;
        }
        else
        {
            // Just copy over values
            Tout[0][d] = Tin[0][s];
            Tout[1][d] = Tin[1][s];
            Tout[2][d] = Tin[2][s];
            Tout[3][d] = Tin[3][s];
            s++;
        }
    }
    Tout.Raw()->maxTraceVal = Tin.Raw()->maxTraceVal;
}



void TraceAlignInsertBases( char cPad, SimpleArray<char>& Envelope, Trace& Tin, Trace& Tout, int nOverlap[2] )
{
/*
    Adds the base calls and base positions from 'Tin' into 'Tout' taking into
    account the alignment.
*/

    // Initialisation
    int     d;
    int     s;
    int     nOrigSamples;
    int     nOverlapL   = nOverlap[0];
    int     nOverlapR   = nOverlap[1];
    int     nBases      = Tin.Raw()->NBases;
    char*   pSrcBase    = Tin.Raw()->base;
    char*   pDstBase    = Tout.Raw()->base;
    uint_2* pSrcBasePos = Tin.Raw()->basePos;
    uint_2* pDstBasePos = Tout.Raw()->basePos;
    uint_2  nBasePos    = 0;


    // Skip over any initial padding
    while( Envelope[nBasePos] == cPad )
        nBasePos++;


    // Insert bases
    for( s=nOverlapL, d=0; s<=nOverlapR && s+1<nBases; s++, d++ )
    {
        pDstBase[d]    = pSrcBase[s];
        pDstBasePos[d] = nBasePos;
        nOrigSamples   = pSrcBasePos[s+1] - pSrcBasePos[s];
        if( s < nOverlapR )
        {
            // Verify bases are not out of order, they can overlap!
            assert(nOrigSamples>=0);


            // Increment sample count, taking into account padding samples
            while( nOrigSamples-- )
            {
                while( Envelope[nBasePos] == cPad )
                    nBasePos++;
                nBasePos++;
            }
        }
    }
}

