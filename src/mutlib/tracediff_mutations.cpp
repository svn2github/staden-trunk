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
#include <cstdio>       // For std::printf(), debugging
#include <cstdlib>      // For std::abs()
#include <mutlib.h>
#include <list.hpp>
#include <array.hpp>
#include <trace.hpp>
#include <muttag.hpp>
#include <peakcall.hpp>
#include <tracediff_parameters.hpp>


//#define VERBOSE_DEBUG


void TraceDiffFindPotentialMutations( Trace& DiffTrace, mutlib_strand_t nStrand,
   int nBaseInterval, int nPosition, int nNoiseThreshold, int nPeakAlignmentThreshold,
   int nPeakTooWideThreshold, double nGlobalMean, List<MutTag>& Mutation )
{
/*
    Forms the core of the double-peak mutation scanning algorithm which is called
    once for each base. We look for double peaks above and below the mean, then
    apply a series of filters to remove the obviously invalid candidates. Each
    candidate is added as a mutation tag to the tag list. Subsequent phases
    whittle these down some.
*/
    int          nP1;
    int          nP2;
    int          x[2];
    int          a[2];
    int          l[2];
    int          r[2];
    int          p[2];
    MutTag*      pTag;
    int          nPeaks;
    int          nMeasurementThreshold;
    PeakCall     PosPeakCall;
    PeakCall     NegPeakCall;
    peak_call_t& PosPeak = PosPeakCall.Data;
    peak_call_t& NegPeak = NegPeakCall.Data;
    MutTag       TmpTag( "MUTA", MUTLIB_MUTATION_NONE, nPosition, nStrand );



    // Set search window limits
    DiffTrace.WindowCentredAt( nPosition, static_cast<int>(nBaseInterval*1.4), l[0], r[0] );



    // Search each trace for +ve and -ve peaks
    for( int j=0; j<4; j++ )
    {
        // No peaks found yet
        PosPeak.Position[j] = -1;
        NegPeak.Position[j] = -1;


        // Scan for nicely formed positive and negative peaks.
        x[0] = DiffTrace.PosPeakFindLargest( j, l[0], r[0], nPeaks, 2 );
        x[1] = DiffTrace.NegPeakFindLargest( j, l[0], r[0], nPeaks, 2 );


        // Save peak positions and amplitudes. The peak amplitude is stored
        // as a signed value centred around the baseline.
        if( x[0] >= 0 )
        {
            PosPeak.Position[j]  = x[0];
            PosPeak.Amplitude[j] = static_cast<int>(double(DiffTrace[j][PosPeak.Position[j]]) - nGlobalMean);
        }
        if( x[1] >= 0 )
        {
            NegPeak.Position[j]  = x[1];
            NegPeak.Amplitude[j] = static_cast<int>(double(DiffTrace[j][NegPeak.Position[j]]) - nGlobalMean);
        }
    }



    // Any valid peaks found at all?
    if( !PosPeakCall.IsValid() || !NegPeakCall.IsValid() )
        return;


    // Get the biggest pos/neg peaks and their corresponding positions
    nP1  = PosPeakCall.MaxAmplitudeAsIndex();
    nP2  = NegPeakCall.MinAmplitudeAsIndex();
    x[0] = PosPeak.Position[nP1];
    x[1] = NegPeak.Position[nP2];
    a[0] = PosPeak.Amplitude[nP1];
    a[1] = NegPeak.Amplitude[nP2];



    // Filter out noise and any single peaks
    if( (nP1==nP2) || (x[0]<0) || (x[1]<0) )
        return;



    // Filter out +ve peaks below the mean and -ve peaks above the mean
    if( (a[0]<=0) || (a[1]>=0) )
        return;



    // Filter out any doublets with either peak below the noise threshold
    a[0] = std::abs( a[0] );
    a[1] = std::abs( a[1] );
    if( (a[0]<nNoiseThreshold) || (a[1]<nNoiseThreshold) )
        return;



    // Measure peak widths at 1/3 of peak amplitude
    const double t = 0.33;
    nMeasurementThreshold = static_cast<int>( nGlobalMean + (a[0] * t) );
    x[0] = DiffTrace.PosPeakWidth( nP1, PosPeak.Position[nP1], l[0], r[0], nMeasurementThreshold );
    nMeasurementThreshold = static_cast<int>( nGlobalMean - (a[1] * t) );
    x[1] = DiffTrace.NegPeakWidth( nP2, NegPeak.Position[nP2], l[1], r[1], nMeasurementThreshold );



    // Save doublet width
    assert(nBaseInterval>0);
    double w = static_cast<double>( x[0]>x[1] ? x[0] : x[1] );
    TmpTag.Width( w/nBaseInterval );



    // Estimate peak centre positions, this is necessary because peaks are
    // often distorted during dynamic programming or because of trace noise
    p[0] = l[0] + ((r[0] - l[0]) / 2);
    p[1] = l[1] + ((r[1] - l[1]) / 2);



    // Filter out any doublets with poor peak centre alignment, save alignment
    int alignment = std::abs( p[0]-p[1] );
    if( alignment > nPeakAlignmentThreshold )
        return;
    TmpTag.Alignment( static_cast<double>(alignment)/static_cast<double>(nBaseInterval) );



    // Filter out any doublets that are too wide
    if( (x[0]>nPeakTooWideThreshold) || (x[1]>nPeakTooWideThreshold) )
        return;



    // Add mutation tag to list, sample position is midway between the two peaks
    pTag = new MutTag( TmpTag );
    pTag->Type( nP1, nP2 );
    pTag->Amplitude( 0, a[0] );
    pTag->Amplitude( 1, a[1] );
    x[0] = PosPeak.Position[nP1];
    x[1] = NegPeak.Position[nP2];
    pTag->Position( 0, (x[0]>x[1]) ? x[1]+((x[0]-x[1])/2) : x[0]+((x[1]-x[0])/2) );
    Mutation.Append( pTag );
}



void TraceDiffComputeLocalEnvelopeStatistics( Trace& DiffTrace, int nPosition, int nNoiseWindow, NumericArray<int>& rEnvelope, double& nMean, double& nSD )
{
/*
   Attempts to compute local envelope statistics over a window. It isn't always
   very successful in the presence of mutations which can cause undue bias.
*/
    int c;
    int n;
    int k;
    int nMin;
    int nMax;
    int nEnd;
    int nBegin;



    // Work out window limits in samples
    DiffTrace.WindowToLeftOf( nPosition, nNoiseWindow, nBegin, nEnd );



    // Ensure we have enough array space
    c = nEnd - nBegin + 1;
    if( c > rEnvelope.Capacity() )
    {
        rEnvelope.Empty();
        rEnvelope.Create( c );
    }
    else
    {
        rEnvelope.Length( c );
    }



    // Extract the envelope
    for( k=0, n=nBegin; n<=nEnd; n++, k++ )
    {
        DiffTrace.MaxAt( n, c, nMax );
        DiffTrace.MinAt( n, c, nMin );
        rEnvelope[k] = nMax - nMin;
    }



    // Compute envelope statistics over window
    nMean = rEnvelope.Mean();
    nSD   = rEnvelope.StandardDeviation( &nMean );
}



void TraceDiffMarkMutationsAboveThreshold( Trace& DiffTrace, double nSensitivity,
   int nNoiseWindow, MutTag& Tag, NumericArray<int>& DiffEnvelope, int& nLastMutation,
   double& nLocalMean, double& nLocalSD )
{
    // Sometimes there is too much variance at the start of a trace difference
    // which can cause us to miss obvious mutations. This might be because the
    // quality clipping is not severe enough, or it could be that the presence
    // of a mutation biases the statistics unfairly. Sometimes the dynamic
    // programming algorithm has problems at the endpoints. So for the first
    // region of the trace, we increase the window over which we compute the
    // statistics in order to get a better estimate.


    // Case 1: First window in trace
    if( Tag.Position() < nNoiseWindow )
    {
        TraceDiffComputeLocalEnvelopeStatistics( DiffTrace, Tag.Position(),
            3*nNoiseWindow, DiffEnvelope, nLocalMean, nLocalSD );
    }



    // Case 2: Normal case, no nearby mutations
    if( Tag.Position()-nLastMutation > nNoiseWindow )
    {
        TraceDiffComputeLocalEnvelopeStatistics( DiffTrace, Tag.Position(),
            nNoiseWindow, DiffEnvelope, nLocalMean, nLocalSD );
    }



    // Case 3: Nearby mutation, use previous window's statistics
    // Values already passed in from last time.



    // Filter out doublets below our threshold
    int nAmplitude        = Tag.Amplitude(0) + Tag.Amplitude(1);
    int nDoubletThreshold = static_cast<int>( nLocalMean + nSensitivity*nLocalSD );
    if( nAmplitude < nDoubletThreshold )
        return;



    // Update mutation tag
    Tag.Confidence( 100 );
    Tag.Sensitivity( (static_cast<double>(nAmplitude)-nLocalMean)/nLocalSD );
    nLastMutation = Tag.Position();
}



void TraceDiffMarkMutationsNearby( Trace& DiffTrace, int nNoiseWindow, MutTag& Tag, MutTag* pTagLast )
{
    // If there was no previous mutation, exit
    if( !pTagLast )
        return;


    // If we are already certain of the current mutation, exit
    if( Tag.Confidence() > 0 )
        return;


    // If distance between previous mutation and current candidate is too great, exit
    int nDistance = Tag.Position() - pTagLast->Position();
    if( nDistance > nNoiseWindow )
        return;


    // Mark as a mutation with reduced confidence
    Tag.Confidence( 50 );
}



void TraceDiffScanForMutations( Trace& DiffTrace, mutlib_strand_t nStrand, int nBaseInterval,
    int nFirstBase, TraceDiffParameters& Parameter, List<MutTag>& Mutation )
{
/*
    The overall automated mutation analysis algorithm.
*/
    assert(nFirstBase>=0);
    assert(nBaseInterval>0);


    // Unpack parameters
    double nSensitivity     = Parameter[TRACEDIFF_PARAMETER_SENSITIVITY].Value();
    double nNoiseThreshold_ = Parameter[TRACEDIFF_PARAMETER_NOISE_THRESHOLD].Value();
    int    nNoiseWindow     = static_cast<int>( Parameter[TRACEDIFF_PARAMETER_NOISE_WINDOW_LENGTH].Value() );
    double nPeakAlignment   = Parameter[TRACEDIFF_PARAMETER_PEAK_ALIGNMENT].Value();
    double nPeakWidthMax    = Parameter[TRACEDIFF_PARAMETER_PEAK_WIDTH_MAXIMUM].Value();



    // Collect some useful statistics
    double nLocalSD    = 0.0;
    double nLocalMean  = 0.0;
    double nGlobalMax  = DiffTrace.Max();
    double nGlobalMean = DiffTrace.Baseline();



    // Create some absolute thresholds/values from the fractional parameters
    nNoiseWindow                = nNoiseWindow * nBaseInterval;
    int nNoiseThreshold         = static_cast<int>( nGlobalMax * nNoiseThreshold_ / 2.0 );
    int nPeakAlignmentThreshold = static_cast<int>( nPeakAlignment * nBaseInterval );
    int nPeakTooWideThreshold   = static_cast<int>( nPeakWidthMax * nBaseInterval );



    // Variable initialisation
    int               n;
    NumericArray<int> DiffEnvelope;
    MutTag*           pTag          = 0;
    MutTag*           pTagLast      = 0;
    int               nSamples      = DiffTrace.Samples();
    int               nLastMutation = -nNoiseWindow;



    // PHASE 1:
    // Scan difference trace for double peaks which fit a sensible mask using a
    // window of size 'nBaseInterval' with 50% overlap. We cannot rely on base
    // call positions to be centred on peaks.
    for( int s=0; s<nSamples; s += nBaseInterval/2 )
    {
        TraceDiffFindPotentialMutations( DiffTrace, nStrand, nBaseInterval, s,
           nNoiseThreshold, nPeakAlignmentThreshold, nPeakTooWideThreshold,
           nGlobalMean, Mutation );
    }



    // PHASE 2:
    // Map mutation doublet sample positions to nearest input-trace base number.
    pTag = Mutation.First();
    while( pTag )
    {
        n = DiffTrace.BaseNumberFromSamplePosition( pTag->Position() );
        pTag->Position( 1, n+nFirstBase+1 );
        pTag = Mutation.Next();
    }



    // PHASE 3:
    // Eliminate duplicate tags that can occur as a result of the window overlap
    // Note, we use the base numbers as a basis for comparison, not sample numbers
    // otherwise we will run into problems.
    pTag = Mutation.First();
    while( pTag )
    {
        if( pTagLast && (pTag->Position(1)==pTagLast->Position(1)) )
        {
            // Keep mutation with largest variance
            n = Mutation.Index();
            if( pTag->Sensitivity() >= pTagLast->Sensitivity() )
               n--;
            delete Mutation.Remove( n );
            pTag = Mutation.Current();
        }
        pTagLast = pTag;
        pTag     = Mutation.Next();
    }



    // PHASE 3A:
    // Print out some useful debugging information
    #ifdef VERBOSE_DEBUG
    pTag = Mutation.First();
    while( pTag )
    {
        std::printf( "%4d (%3d): MUTA Width=%.2f, Alignment=%.2f, Amplitude=%d\n",
                     pTag->Position(), pTag->Position(1), pTag->Width(),
                     pTag->Alignment(), static_cast<int>(pTag->Amplitude(0)+pTag->Amplitude(1)) );
        pTag = Mutation.Next();
    }
    #endif



    // PHASE 4:
    // Compute the background noise in the local region specified by the 'nNoiseWindow'
    // parameter and use this to mark mutations that are above the noise threshold.
    pTag = Mutation.First();
    while( pTag )
    {
        TraceDiffMarkMutationsAboveThreshold( DiffTrace, nSensitivity, nNoiseWindow,
           *pTag, DiffEnvelope, nLastMutation, nLocalMean, nLocalSD );
        pTag = Mutation.Next();
    }



    #ifdef NOT_USED_ANYMORE
    //
    // After a lot of testing, phase 5 was found to have very little
    // beneficial effect - mainly contributing to the false +ve count.
    // It only seemed to be useful in christines LEA9xx dataset where
    // mutations were clustered together. In normal clinical data this
    // is rare.
    //

    // PHASE 5:
    // Reexamine potential mutations just to the right of the ones found in the
    // last phase. These can sometimes be missed due to bias in the statistics
    // caused by mutations.
    pTagLast = 0;
    pTag     = Mutation.First();
    while( pTag )
    {
        if( pTag->Confidence() == 100 )
            pTagLast = pTag;
        TraceDiffMarkMutationsNearby( DiffTrace, nNoiseWindow, *pTag, pTagLast );
        pTag = Mutation.Next();
    }
    #endif



    // PHASE 6:
    // Remove mutation tags from the list with zero confidence.
    pTag = Mutation.First();
    while( pTag )
    {
        if( pTag->Confidence() <= 0 )
        {
            delete Mutation.Remove( Mutation.Index() );
            pTag = Mutation.Current();
        }
        else
        {
            pTag = Mutation.Next();
        }
    }
}
