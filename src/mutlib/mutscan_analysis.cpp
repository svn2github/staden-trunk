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
#include <cmath>                    // For std::sqrt()
#include <limits>                   // For std::numeric_limits<>::max()
#include <cstring>                  // For std::strcpy(), std::strcat(), std::strcmp()
#include <cstdio>                   // For std::fopen(), ... debugging
#include <cstdlib>                  // For std::abs(), std::qsort()
#include <algorithm>                // For std::max(), std::min()
#include <mutlib.h>                 // For mutlib_input_t
#include <dnatable.hpp>             // For DNATable object
#include <caller_base.hpp>          // For BaseCaller object
#include <caller_level.hpp>         // For LevelCaller object
#include <mutscan_analysis.hpp>


// #define VERBOSE_DEBUG


/*
    Module implementation constants.
*/

static const double SDT = 1.1;                          // >= 1.1 misses mutations
static const int    REF = MUTLIB_INPUT_REFERENCE;
static const int    INP = MUTLIB_INPUT;



/**
    Master routine that performs mutation analysis on the preprocessed input
    data and traces.
*/
mutlib_result_t MutScanAnalyser::Execute( mutscan_t* ms, MutScanPreprocessor Data[], Trace Tr[], Trace* DifferenceTrace )
{
    assert(ms != NULL);
    assert(ms->Initialised);


    // Get algorithm parameters
    m_nSearchWindow           = static_cast<int>( Data[REF].SearchWindow / 2.0 + 0.5 );
    m_nHetSNRThreshold        = ms->Parameter[MUTSCAN_PARAMETER_HETSNR_THRESHOLD];
    m_nLowerPeakDropThreshold = ms->Parameter[MUTSCAN_PARAMETER_PEAKDROP_LOWER];
    m_nUpperPeakDropThreshold = ms->Parameter[MUTSCAN_PARAMETER_PEAKDROP_UPPER];
    mutlib_strand_t Strand    = ms->InputTrace.Strand;
    assert(m_nSearchWindow>=1);
    assert(m_nLowerPeakDropThreshold!=0.0);
    assert(m_nUpperPeakDropThreshold!=0.0);



    // Form matrices of aligned peak pairs
    AllocatePeakMap( Data );
    AlignPeaks( Data );


    // Compute the scale factors required to bring the input trace up
    // to the same level as the reference trace.
    ComputeScaleFactors( Data );


    // Do mutation detection
    ScanForPotentialMutations( Data, Strand, Tr );
    AnalysePotentialMutations( Tr );


    // If difference trace supplied, do further analysis and verification
    if( DifferenceTrace )
    {
        // Check that each mutation has a corresponding double peak
        ValidateMutationsAgainstDifference( *DifferenceTrace );
    }

    return MUTLIB_RESULT_SUCCESS;
}



/**
    Creates the two peak map matrices of the appropriate dimensions and sets
    their contents to zero.
*/
void MutScanAnalyser::AllocatePeakMap( MutScanPreprocessor Data[] )
{
    int max_peaks = std::max( Data[REF].PeakCountMax, Data[INP].PeakCountMax );
    Map.Create( 2*4, max_peaks );
    Map.Fill( 0 );
    for( int k=0; k<4; k++ )
        MapCount[k] = 0;
}



/**
    Steps along the reference peaks array peak-by-peak, one base at a time
    and searches for corresponding peaks in the input trace over a small
    window size. As the aligned peaks are found, their positions are added
    into the peak mapping matrix.
*/
void MutScanAnalyser::AlignPeaks( MutScanPreprocessor Data[] )
{
    // Construct a peak alignment map...
    const int cols = Data[REF].Peak.Cols();
    for( int r=0, k=0; r<4; r++, k=0 )
    {
        for( int c=0; c<cols; c++ )
        {
            // If peak in reference?
            if( Data[REF].Peak[r][c] > 0 )
            {
                if( Data[INP].Peak[r][c] > 0 )
                {
                    // Peaks were aligned exactly
                    Map[2*r+1][k] = c;
                }
                else
                {
                    // Look within search window for the largest peak, beginning
                    // in the middle and working our way towards the outsides.
                    int max_c = -1;
                    int max_a =  0;
                    int j     =  c;
                    int off   =  1;
                    int win   =  m_nSearchWindow;
                    while( (win>0) && ((j-off)>=0) && ((j+off)<cols) )
                    {
                        if( Data[INP].Peak[r][j-off] > max_a )
                        {
                            max_a = Data[INP].Peak[r][j-off];
                            max_c = j-off;
                        }
                        if( Data[INP].Peak[r][j+off] > max_a )
                        {
                            max_a = Data[INP].Peak[r][j+off];
                            max_c = j+off;
                        }
                        off++;
                        win--;
                    }
                    Map[2*r+1][k] = (max_a>0) ? max_c : 0;
                }
                Map[2*r][k] = c;
                k++;
            }
        }
        MapCount[r] = k;
    }
    #ifdef VERBOSE_DEBUG
    Map.SaveAs( "peak_alignment_matrix.txt" );
    #endif
}



/**
    Computes the upper and lower scale factor limits based on the mean and
    standard deviation statistics for the current base (r).  2*SD=95% and
    3*SD=99%.
*/
void MutScanAnalyser::ComputeScaleFactorLimits( int r, double nosd, double lim[2] )
{
    lim[0]  = ScaleFactorMean[r] - nosd*ScaleFactorStdDev[r];
    lim[1]  = ScaleFactorMean[r] + nosd*ScaleFactorStdDev[r];

    // Occasionally the SD can be very large, so we clip it to zero
    if( lim[0] < 0.0 )
        lim[0] = 0.0;
}



/**
    Steps along the peak alignment map pairs. For each pair we calculate
    the input peak scale factor. This is the scale factor required to normalise
    the input peak amplitude to that of the reference. This is done for each
    base individually giving 4 rows of scale factors, since they are not usually
    the same due to differences in chemistry.
*/
void MutScanAnalyser::ComputeScaleFactors( MutScanPreprocessor Data[] )
{
    // Initialisation
    const int rows = 4;
    const int cols = Map.Cols();


    // Allocate and initialise scale factor matrices
    ScaleFactor.Create( rows, cols );
    ScaleFactor.Fill( 0.0 );
    for( int n=0; n<rows; n++ )
    {
        ScaleFactorMean[n]   = 1.0;
        ScaleFactorStdDev[n] = 0.0;
    }


    // Create some temporary value storage
    int    idx[2];
    double amp[2];
    NumericArray<double> data;
    data.Create( cols );


    // Compute the scale factors
    for( int r=0; r<rows; r++ )
    {
        int n = 0;
        for( int c=0; c<MapCount[r]; c++ )
        {
            // If an aligned peak pair is found
            idx[REF] = Map[2*r][c];
            idx[INP] = Map[2*r+1][c];
            if( (idx[REF]>0) && (idx[INP]>0) )
            {
                // Compute the scale factor
                amp[REF] = Data[REF].Peak[r][ idx[REF] ];
                amp[INP] = Data[INP].Peak[r][ idx[INP] ];
                ScaleFactor[r][c] = amp[REF] / amp[INP];
                data[n] = ScaleFactor[r][c];
                n++;
            }
        }
        data.Length( n );
        ScaleFactorMean[r]   = data.Mean();
        ScaleFactorStdDev[r] = (n>1) ? std::sqrt(data.Variance(&ScaleFactorMean[r])) : 0.0;
    }


    #ifdef VERBOSE_DEBUG
    ScaleFactor.SaveAs( "scale_factor_matrix.txt", 10 );
    std::printf( "Scale factor means: A=%0.2f, C=%0.2f, G=%0.2f, T=%0.2f\n",
                  ScaleFactorMean[0], ScaleFactorMean[1], ScaleFactorMean[2], ScaleFactorMean[3] );
    std::printf( "Scale factor SD's : A=%0.2f, C=%0.2f, G=%0.2f, T=%0.2f\n",
                  ScaleFactorStdDev[0], ScaleFactorStdDev[1], ScaleFactorStdDev[2], ScaleFactorStdDev[3] );
    #endif
}



/**
    Does a level call on the trace and copies the results into the tag.
*/
void MutScanAnalyser::DoLevelCall( int pos, Trace& Tr, MutationTag& Tag, bool AllowAmbiguity )
{
    DNATable    Table;
    char        base[3];
    LevelCaller LevCaller( Tr, pos );

    if( AllowAmbiguity )
    {
        base[1] = LevCaller.GetBase(3);
        base[2] = LevCaller.GetBase(2);
        base[0] = Table.LookupBase( base[1], base[2] );
        Tag.BaseInput( 0, base[0] );
        Tag.BaseInput( 1, base[1] );
        Tag.BaseInput( 2, base[2] );
    }
    else
    {
        Tag.BaseInput( 0, LevCaller.GetBase(3) );
        Tag.BaseInput( 1, LevCaller.GetBase(3) );
        Tag.BaseInput( 2, '-' );
    }
}


/**
    Searches the peak map for an input peak on 'base' and looks to
    see if it's aligned with a reference peak.
*/
bool MutScanAnalyser::HasReferencePeak( int base, int pos )
{
    assert(base>=0);
    assert(base<=3);
    for( int c=0; c<MapCount[base]; c++ )
    {
        if( Map[2*base+1][c] == pos )
        {
            if( Map[2*base]!=0 )
                return true;
            break;
        }
    }
    return false;
}



/**
    Scans the scale factor matrix looking for potential mutation sites
    and creates a mutation tag for each site.
*/
void MutScanAnalyser::ScanForPotentialMutations( MutScanPreprocessor Data[], mutlib_strand_t Strand, Trace Tr[] )
{
    int       n;
    DNATable  Table;
    double    lim[2];
    const int rows = 4;


    for( int r=0; r<rows; r++ )
    {
        ComputeScaleFactorLimits( r, SDT, lim );
        for( int c=0; c<MapCount[r]; c++ )
        {
            // If we've found a potential mutation site...
            double sf = ScaleFactor[r][c];
            if( (sf<=lim[0]) || (sf>lim[1]) )
            {
                // Gather some information, watch out for single peaks
                bool NoPeak = false;
                int  RefPos = Map[2*r][c];
                int  InpPos = Map[2*r+1][c];
                if( InpPos <= 0 )
                {
                    InpPos = RefPos;
                    NoPeak = true;
                }


                // Do reference basecall. Search window is set to a low value
                // since these peaks should exist exactly at RefPos and we don't
                // want to inadvertantly pick up peaks further away. Heterozygous
                // reference bases are presently ignored.
                BaseCaller RefCaller( Tr[REF], Data[REF].Peak, RefPos, 1 );
                if( Table.IsBaseAmbiguous(RefCaller.GetBase()) )
                    continue;


                // Do input basecall, decide if it's heterozygous based on SNR & noise floor
                BaseCaller InpCaller( Tr[INP], Data[INP].Peak, InpPos, m_nSearchWindow );
                bool Heterozygous = false;
                if( InpCaller.GetConfidence() < m_nHetSNRThreshold )
                {
                    if( static_cast<int>(Tr[INP][r][InpPos]) > Data[INP].Noise[InpPos] )
                        Heterozygous = true;
                }


                // Create a mutation tag, set reference basecall
                MutationTag* pTag = new MutationTag( Heterozygous ? "HETE" : "MUTA" );
                pTag->BaseReference( RefCaller.GetBase() );


                // Tag basecall assignment state machine
                enum { STATE_HETEROZYGOUS, STATE_NO_PEAK, STATE_BASE_CHANGE, STATE_NO_MUTATION };
                bool done  = false;
                int  state = STATE_HETEROZYGOUS;
                while( !done )
                {
                    switch( state )
                    {
                        case STATE_HETEROZYGOUS:
                            // If heterozygous, do a level call to determine tag bases
                            state = STATE_NO_PEAK;
                            if( Heterozygous )
                            {
                                DoLevelCall( InpPos, Tr[INP], *pTag, true );
                                done = true;
                            }
                            // If double peak found, but SNR isn't low, its not heterozygous
                            else if( Table.IsBaseAmbiguous(InpCaller.GetBase()) )
                                InpCaller.Invalidate();
                            break;


                        case STATE_NO_PEAK:
                            // If no peak found, do a level call without ambiguity to determine tag base
                            state = STATE_BASE_CHANGE;
                            if( !InpCaller.IsValid() )
                            {
                                // Discard if base is same as reference
                                DoLevelCall( InpPos, Tr[INP], *pTag, false );
                                if( RefCaller.GetBase() == pTag->BaseInput() )
                                    state = STATE_NO_MUTATION;
                                else
                                    done = true;
                            }
                            break;


                        case STATE_BASE_CHANGE:
                            // Normal case, input call is valid and unambiguous, assign base to tag
                            pTag->BaseInput( 0, InpCaller.GetBase(0) );
                            pTag->BaseInput( 1, InpCaller.GetBase(1) );
                            pTag->BaseInput( 2, InpCaller.GetBase(2) );
                            // Verify that we haven't inadvertently latched onto a peak
                            // next door that is already matched up with the reference.
                            n = Table.LookupIndex( InpCaller.GetBase() );
                            if( HasReferencePeak(n,InpCaller.GetPosition()) )
                            {
                                InpCaller.Invalidate();
                                state = STATE_NO_PEAK;
                            }
                            else
                                done = true;
                            break;


                        case STATE_NO_MUTATION:
                            // Candidate was rejected, not a mutation
                            delete pTag;
                            done = true;
                            break;
                    }
                }


                // Nothing found...
                if( state == STATE_NO_MUTATION )
                    continue;


                // Fill in remaining information fields...
                pTag->Row( r );
                pTag->Col( c );
                pTag->Strand( Strand );
                pTag->Position( 1, InpPos );
                pTag->SNR( InpCaller.GetConfidence() );
                pTag->Amplitude( 0, Data[REF].Peak[r][RefPos] );
                pTag->Amplitude( 1, NoPeak ? Tr[INP][r][InpPos] : Data[INP].Peak[r][InpPos] );
                MutationTagList.Append( pTag );
            }
        }
    }
}



/**
    Does further analysis on the potential mutation tag to determine:

    a) If it's indeed a valid mutation
    b) To determine if it's a base-change or heterozygous mutation
    c) To compute the peak drop wrt to the reference

*/
void MutScanAnalyser::AnalysePotentialMutations( Trace Tr[] )
{
    double lim[2];
    #ifdef VERBOSE_DEBUG
    std::printf("Analysing %d potential mutations.\n", MutationTagList.Count() );
    #endif
    MutationTag* pTag = MutationTagList.First();
    while( pTag )
    {
        // Compute normalised peak ratio wrt to reference
        double sf = ScaleFactorMean[ pTag->Row() ];
        if( pTag->Col() > 0 )
        {
            sf = ScaleFactor[ pTag->Row() ][ pTag->Col()-1 ];
            ComputeScaleFactorLimits( pTag->Row(), SDT, lim );
            if( (sf<lim[0]) || (sf>lim[1]) )
                sf = ScaleFactorMean[ pTag->Row() ];
        }
        assert(pTag->Amplitude(0)!=0.0);
        pTag->Amplitude( 2, pTag->Amplitude(1) * sf / pTag->Amplitude(0) );


        // Heterozyous mutation?
        if( std::strcmp(pTag->Name(),"HETE") == 0 )
        {
            // If peak drop is above/below our thresholds, mark for deletion
            double pkd = pTag->Amplitude(2);
            if( (pkd<m_nLowerPeakDropThreshold)||(pkd>m_nUpperPeakDropThreshold) )
                pTag->Mark( true );
        }
        pTag = MutationTagList.Next();
    }
}



/**
    Use the difference trace to verify that there's a double peak at site of
    each tagged mutation. If there isn't, we mark the mutation for deletion.
*/
void MutScanAnalyser::ValidateMutationsAgainstDifference( Trace& DiffTrace )
{
    // Verify the tag against the difference trace
    const int baseline = DiffTrace.Baseline();
    MutationTag* pTag  = MutationTagList.First();
    while( pTag )
    {
        if( !pTag->Marked() )
        {
            // Check for the presence of a double peak
            int count = 0;
            for( int r=0; r<4; r++ )
            {
                if( DiffTrace[r][pTag->Position(1)] == baseline )
                    continue;
                count++;
            }


            // If no peaks found, mark for deletion
            if( count == 0 )
                pTag->Mark( true );
        }
        pTag = MutationTagList.Next();
    }
}

