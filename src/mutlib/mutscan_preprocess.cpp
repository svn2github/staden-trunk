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
#include <cstring>                  // For std::strcpy(), std::strcat()
#include <cstdio>                   // For std::sprintf()
#include <mutscan_preprocess.hpp>


// #define VERBOSE_DEBUG



/**
   Master routine that invokes each preprocessing stage in turn.
*/
mutlib_result_t MutScanPreprocessor::Execute( mutscan_t* ms, Trace& t, int n, int left, int right )
{
    assert(ms != NULL);
    assert(ms->Initialised);


    // Get algorithm parameters, we set the adaptive noise threshold for
    // the reference to be twice that of the input. This is so we have a
    // rock solid reference signal without missing low level peaks on the
    // input trace.
    PeakInterval         = t.IntervalMode();
    SearchWindow         = ms->Parameter[MUTSCAN_PARAMETER_SEARCH_WINDOW] * PeakInterval;
    m_nNoiseThreshold[0] = ms->Parameter[MUTSCAN_PARAMETER_NOISE_THRESHOLD] * 2.0;
    m_nNoiseThreshold[1] = ms->Parameter[MUTSCAN_PARAMETER_NOISE_THRESHOLD];
    assert(SearchWindow>1.0);
    assert(m_nNoiseThreshold[0]!=0.0);
    assert(m_nNoiseThreshold[1]!=0.0);



    // Start preprocessing
    PeakFind( t, left, right );
    EstimateNoiseFloor( t, n );
    PeakClip();
    PeakSpacing();
    CountPeaks();


    if( PeakCount.Max() < 3 )
    {
        ms->ResultCode = MUTLIB_RESULT_INSUFFICIENT_DATA;
        std::sprintf( ms->ResultString, "Insufficent data to process trace %s.\n", t.Name() );
        return MUTLIB_RESULT_INSUFFICIENT_DATA;
    }



    return MUTLIB_RESULT_SUCCESS;
}



/**
    Scans the trace for peaks and records their amplitude and position
    in the 4 * SAMPLES Peak matrix.

    To avoid discontinuities at the edges, we only search for peaks
    within the two margins, left/right. This avoids false positives.
*/
void MutScanPreprocessor::PeakFind( Trace& Tr, int left, int right )
{
    int   pos;
    int   resume;
    const int rows = 4;
    const int cols = Tr.Samples();



    // Allocate and initialise the peak matrix
    Peak.Create( rows+1, cols );
    Peak.Fill( 0 );



    // Create a peaks matrix from the trace
    for( int r=0; r<rows; r++ )
    {
        resume = left;
        while(1)
        {
            pos = Tr.PosPeakFind( r, resume, right, resume, 1 );
            if( pos > 0 )
                Peak[r][pos] = Tr[r][pos];
            else
                break;
        }
    }
}



/**
    Computes an estimate of the trace noise floor from the envelope
    and records it in the 1*SAMPLES Noise vector.
*/
void MutScanPreprocessor::EstimateNoiseFloor( Trace& Tr, int n )
{
    int       a[2];
    int       pos;
    int       resume;
    const int rows = 4;
    const int cols = Peak.Cols();



    // Allocate and initialise noise floor vector
    Noise.Create( cols );
    Noise.Fill( 0 );



    // Create the trace envelope
    Trace* pEnvelope = Tr.CreateEnvelope();
    Trace& Envelope  = *pEnvelope;



    // Find all peaks in the envelope
    resume = 0;
    while(1)
    {
        pos = Envelope.PosPeakFind( 0, resume, cols-1, resume, 1 );
        if( pos >= 0 )
        {
            // Noise estimate at current position is a fixed percentage of
            // the envelope peak height
            a[0] = Envelope[0][pos];
            a[1] = static_cast<int>( m_nNoiseThreshold[n] * a[0] );
            Noise[pos] = a[1];
        }
        else
            break;
    }



    // Interpolate through non-zero noise points to complete our noise estimate
    int x1 = 0;
    for( int c=1, end=cols-1; c<cols; c++ )
    {
        if( (Noise[c]>0) || (c==end) )
        {
            Noise.Interpolate( x1, c );
            x1 = c;
        }
    }



    #ifdef VERBOSE_DEBUG
    for( int c=0; c<cols; c++ )
    {
        if( (Peak[0][c]>Noise[c]) || (Peak[1][c]>Noise[c]) || (Peak[2][c]>Noise[c]) || (Peak[3][c]>Noise[c]) )
           Envelope[1][c] = Envelope[0][c];
        Envelope[2][c] = Noise[c];
    }
    char name[64];
    std::sprintf( name, "peaks_and_noise%d.ztr", n+1 );
    Envelope.SaveAs( name );
    #endif



    // Cleanup
    delete pEnvelope;
}



/**
    Removes all peaks below noise floor from the 4*SAMPLES Peak matrix
*/
void MutScanPreprocessor::PeakClip()
{
    const int rows = 4;
    const int cols = Peak.Cols();
    for( int r=0; r<rows; r++ )
    {
        for( int c=0; c<cols; c++ )
        {
            // Remove peaks below the noise floor
            if( (Peak[r][c]>0) && (Peak[r][c]<Noise[c]) )
               Peak[r][c] = 0;
        }
    }
}



/**
    Goes through the 1st four rows of the Peak matrix, and produces a 5th
    row containing a logical OR of them to give a picket fence of peaks
    over all bases.

    This data may be used for shoulder and blob detection which looks for
    peak spacing irregularities.
*/
void MutScanPreprocessor::PeakSpacing()
{
    const int rows = 4;
    const int cols = Peak.Cols();
    for( int c=0; c<cols; c++ )
    {
        for( int r=0; r<rows; r++ )
        {
            if( Peak[r][c] > 0 )
            {
                Peak[rows][c] = 1;
                break;
            }
        }
    }
}


/**
    Most blob problems are avoided because they also occur on the reference.

    Elimate peaks caused by gel blobs. A simple peak width measurement
    is not sufficient as this will also catch two identical and adjacent
    bases if one peak happens to be a shoulder.

    A peak width measurement plus a total peak count over that width is
    also inconclusive.

    Blobs often occur very close to other bases, so this is the measure
    we use - to look for unusually close peak spacings and then we examine
    these regions in more depth.

    The total logical OR of peaks are computed in the last row of the peaks
    matrix.

*/



/**
    Counts the peaks for each base, ie in each row of Peaks[][].
*/
void MutScanPreprocessor::CountPeaks()
{
    const int rows = 4;
    const int cols = Peak.Cols();


    // Allocate and initialise peak count vector
    PeakCount.Create( rows );
    PeakCount.Fill( 0 );


    // Count all peaks
    for( int r=0; r<rows; r++ )
    {
        int cnt = 0;
        for( int c=0; c<cols; c++ )
        {
            if( Peak[r][c] > 0 )
               cnt++;
        }
        PeakCount[r] = cnt;
    }


    // Compute maximum number of peaks overall
    PeakCountMax = PeakCount.Max();
}
