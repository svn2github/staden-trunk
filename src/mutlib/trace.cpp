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


#include <new>                  // For bad_alloc()
#include <climits>              // For INT_MIN, INT_MAX
#include <cstdlib>              // For std::qsort(), std::abs()
#include <cstring>              // For strlen(), strcpy(), ...
#include <staden.h>             // For xalloc()
#include <array.hpp>            // For SimpleArray<T>
#include <algorithm>            // For std::min(), std::max()
#include <trace.hpp>



//------
// Init
//------

void Trace::Init()
{
   m_pRead             = 0;
   m_bAutoDestroy      = true;
   m_bExternalTrace    = false;
   m_nLowerLimit       = 0;
   m_nUpperLimit       = 0;
   m_bStatisticsCached = false;
   ZeroTraces();
}



//-------------
// Zero Traces
//-------------

void Trace::ZeroTraces()
{
   m_pTrace[0] = 0;
   m_pTrace[1] = 0;
   m_pTrace[2] = 0;
   m_pTrace[3] = 0;
}



//-------------
// Init Traces
//-------------

void Trace::InitTraces()
{
   // This horrible hack is required due to poor Read structure design
   if( !m_pRead )
      ZeroTraces();
   else
   {
      m_pTrace[0] = m_pRead->traceA;
      m_pTrace[1] = m_pRead->traceC;
      m_pTrace[2] = m_pRead->traceG;
      m_pTrace[3] = m_pRead->traceT;
   }
}



//-----------------
// Base Confidence
//-----------------

int Trace::BaseConfidence( int n ) const
{
   assert(n>=0);
   assert(m_pRead!=0);
   assert(n<m_pRead->NBases);
   switch( m_pRead->base[n] )
   {
      case 'A':
      case 'a':  return m_pRead->prob_A ? m_pRead->prob_A[n] : 0;
      case 'C':
      case 'c':  return m_pRead->prob_C ? m_pRead->prob_C[n] : 0;
      case 'G':
      case 'g':  return m_pRead->prob_G ? m_pRead->prob_G[n] : 0;
      case 'T':
      case 't':  return m_pRead->prob_T ? m_pRead->prob_T[n] : 0;
   }
   return 0;
}



//------
// Wrap
//------

void Trace::Wrap( Read* r, bool AutoDestroy )
{
   assert(r!=0);
   m_pRead          = r;
   m_bAutoDestroy   = AutoDestroy;
   m_bExternalTrace = true;
   InitTraces();
   Range( 0, r->NBases ? r->NBases - 1 : 0 );
}



//------
// Open
//------

bool Trace::Open( const char* pFileName )
{
   if( !m_bExternalTrace )
   {
      m_pRead = read_reading( (char*)pFileName, TT_ANY );
      if( m_pRead )
      {
         InitTraces();
         Range( 0, m_pRead->NBases ? m_pRead->NBases - 1 : 0 );
         return true;
      }
   }
   return false;
}



//-------
// Close
//-------

void Trace::Close()
{
   if( (m_pRead!=0) && m_bAutoDestroy )
      read_deallocate(m_pRead);
   Init();
}



//--------
// Create
//--------

bool Trace::Create( int nSamples, int nBases, const char* pName )
{
   assert(nBases>=0);
   assert(nSamples>=0);

   if( !m_bExternalTrace )
   {
      m_pRead = read_allocate( nSamples, nBases );
      if( m_pRead )
      {
         if( pName )
         {
             m_pRead->trace_name = static_cast<char*>( xmalloc((std::strlen(pName)+1)*sizeof(char)) );
             std::strcpy( m_pRead->trace_name, pName );
         }
         InitTraces();
         Range( 0, nBases ? nBases-1 : 0 );
         return true;
      }
   }
   return false;
}



//---------
// Save As
//---------

bool Trace::SaveAs( const char* pFileName, int fmt )
{
   assert(m_pRead!=0);
   int ret = write_reading( (char*) pFileName, m_pRead, fmt );
   return ret ? false : true;
}



//-------
// MaxAt
//-------

void Trace::MaxAt( int n, int& c, int& m ) const
{
/*
    Computes the maximum value at sample 'n' and returns it in 'm' along with
    it's corresponding trace channel in 'c'.
*/
    c = -1;
    m = INT_MIN;
    if( int(m_pTrace[0][n]) > m )
    {
        m = int( m_pTrace[0][n] );
        c = 0;
    }
    if( int(m_pTrace[1][n]) > m )
    {
        m = int( m_pTrace[1][n] );
        c = 1;
    }
    if( int(m_pTrace[2][n]) > m )
    {
        m = int( m_pTrace[2][n] );
        c = 2;
    }
    if( int(m_pTrace[3][n]) > m )
    {
        m = int( m_pTrace[3][n] );
        c = 3;
    }
}



//-------
// MinAt
//-------

void Trace::MinAt( int n, int& c, int& m ) const
{
/*
    Computes the minimum value at sample 'n' and returns it in 'm' along with
    it's corresponding trace channel in 'c'.
*/
    c = -1;
    m = INT_MAX;
    if( int(m_pTrace[0][n]) < m )
    {
        m = int( m_pTrace[0][n] );
        c = 0;
    }
    if( int(m_pTrace[1][n]) < m )
    {
        m = int( m_pTrace[1][n] );
        c = 1;
    }
    if( int(m_pTrace[2][n]) < m )
    {
        m = int( m_pTrace[2][n] );
        c = 2;
    }
    if( int(m_pTrace[3][n]) < m )
    {
        m = int( m_pTrace[3][n] );
        c = 3;
    }
}




//---------------
// Window Limits
//---------------

void Trace::WindowCentredAt( int nPosition, int nSize, int& l, int& r ) const
{
/*
    Computes the limits, 'l', 'r' in samples for a window of maximum size
    'nSize' positioned at 'nCentre'. The window size may be reduced at each
    end of the trace.
*/
    assert(m_pRead!=0);
    assert(nSize>0);
    assert(nPosition>=0);

    l = nPosition - nSize/2;
    r = nPosition + nSize/2;
    if( l < 0 )
        l = 0;
    if( r >= Samples() )
        r = Samples()-1;
}


void Trace::WindowToLeftOf( int nPosition, int nSize, int& l, int& r ) const
{
/*
    Computes the left and right limits 'l', 'r' in samples for a window of size
    'nSize' samples to the left of 'nPosition'. The window size may be reduced
    near the start of the trace.
*/
    assert(nSize>0);
    assert(nPosition<Samples());



    // Input checking
    l = 0;
    r = 0;
    if( (nPosition>=Samples()) || (nSize<=0) )
        return;



    // If region would be less than window size
    if( nPosition < nSize )
    {
        // Force to be size of window
        l = 0;
        r = nSize - 1;
        if( r >= Samples() )
            r = Samples() - 1;
    }
    else
    {
        // Compute naive limits, left clip
        r = nPosition - 1;
        l = r - nSize + 1;
        if( l < 0 )
            l = 0;
    }
}



//---------------------
// Positive Peak Width
//---------------------

int Trace::PosPeakWidth( int n, int nPos, int& nLeft, int& nRight, int nMeasurementThreshold ) const
{
/*
    Measures the peak width of trace 'n' at 'nPos' stopping when the
   'nMeasurementThreshold' is reached. The peak width is returned and
    the left/right positions are returned in 'left' and 'right'.
*/
    assert(n<4);
    int e, k;
    TRACE* t = m_pTrace[n];



    // Search left
    for( k=nPos; k>0; k-- )
    {
        // If reached start or hit the threshold
        if( (k==1) || (t[k]<=nMeasurementThreshold) )
        {
            nLeft = k;
            break;
        }
    }



    // Search right
    e = Samples() - 1;
    for( k=nPos; k<e; k++ )
    {
        // If reached end or hit the threshold
        if( (k==(e-1)) || (t[k]<=nMeasurementThreshold) )
        {
            nRight = k;
            break;
        }
    }
    return nRight - nLeft;
}



int Trace::PosPeakWidth( int n, int nPos, int& nLeft, int& nRight ) const
{
/*
    Measures the peak width of trace 'n' at 'nPos' stopping when at
    a zero slope or a slope sign change is encountered. The peak width
    is returned and the left/right positions are also returned in 'left'
    and 'right'.
*/
    assert(n<4);
    int k;
    int end;
    int slope;
    TRACE* t = m_pTrace[n];



    // Search left until we hit a slope sign change, zero or the end
    for( k=nPos; k>0; k-- )
    {
        slope = t[k] - t[k-1];
        if( (slope<=0) || (k==1) )
        {
            nLeft = k;
            break;
        }
    }



    // Search right until we hit a slope change, zero or the end
    end = Samples() - 2;
    for( k=nPos; k<end; k++ )
    {
        slope = t[k] - t[k+1];
        if( (slope<=0) || (k==end) )
        {
            nRight = k;
            break;
        }
    }
    return nRight - nLeft;
}



//---------------------
// Negative Peak Width
//---------------------

int Trace::NegPeakWidth( int n, int nPos, int& nLeft, int& nRight, int nMeasurementThreshold ) const
{
/*
    Same as PosPeakWidth, except it looks for inverted peaks.
*/
    assert(n<4);
    int e, k;
    TRACE* t = m_pTrace[n];



    // Search left
    for( k=nPos; k>0; k-- )
    {
        // If reached start or hit the threshold
        if( (k==1) || (t[k]>=nMeasurementThreshold) )
        {
            nLeft = k;
            break;
        }
    }



    // Search right
    e = Samples() - 1;
    for( k=nPos; k<e; k++ )
    {
        // If reached end or hit the threshold
        if( (k==(e-1)) || (t[k]>=nMeasurementThreshold) )
        {
            nRight = k;
            break;
        }
    }
    return nRight - nLeft;
}



//--------------------
// Positive Peak Find
//--------------------

int Trace::PosPeakFind( int n, int nBegin, int nEnd, int& nResume, int nSlopeCount ) const
{
/*
    Scans trace channel 'n' from 'nBegin' to 'nEnd' looking for a well defined
    positive peak. A peak is defined as SLOPE_COUNT consecutive +ve slopes,
    followed by zero or more occurances of a zero slope, followed by SLOPE_COUNT
    consecutive -ve slopes. This ensures that only well defined peaks are
    found. Returns the peak position on success, otherwise -1.
*/
    assert(n<4);
    enum { STATE_START, STATE_CLIMB, STATE_FLAT, STATE_DESCEND, STATE_PEAK, STATE_NOPEAK };
    int    nCtr;
    int    nSlope;
    int    nZeros = 0;
    int    nPeak = 0;
    int    k      = nBegin;
    TRACE* t      = m_pTrace[n];
    int    nState = STATE_START;



    // Search loop
    while(1)
    {
        switch(nState)
        {
            case STATE_START:
                // Search for an initial upwards slope
                nState = STATE_NOPEAK;
                for( k=k; k<nEnd; k++ )
                {
                    nSlope = t[k+1] - t[k];
                    if( nSlope > 0 )
                    {
                        nState = STATE_CLIMB;
                        break;
                    }
                }
                break;



            case STATE_CLIMB:
                // Climb up a rising slope
                nCtr   = 0;
                nZeros = 0;
                nState = STATE_NOPEAK;
                for( k=k; k<nEnd; k++ )
                {
                    nSlope = t[k+1] - t[k];
                    if( nSlope > 0 )
                        nCtr++;
                    else if( (nSlope==0) && (nCtr>=nSlopeCount) )
                    {
                        nState = STATE_FLAT;
                        break;
                    }
                    else if( (nSlope<0) && (nCtr>=nSlopeCount) )
                    {
                        nPeak  = k;
                        nState = STATE_PEAK;
                        break;
                    }
                    else
                    {
                        nState = STATE_START;
                        break;
                    }
                }
                break;



            case STATE_FLAT:
                // Traverse clipped peak. Allow for a bit
                // of noise on this peak.
                nState = STATE_NOPEAK;
                for( k=k; k<nEnd; k++ )
                {
                    nSlope = t[k+1] - t[k];
                    if( (nSlope>-3) && (nSlope<3) )
                        nZeros++;
                    else if( nSlope > 0 )
                    {
                        nState = STATE_CLIMB;
                        break;
                    }
                    else
                    {
                        nPeak  = k - (nZeros/2);
                        nState = STATE_DESCEND;
                        break;
                    }
                }
                break;



            case STATE_DESCEND:
                // Descend down a falling slope
                nCtr   = 0;
                nState = STATE_NOPEAK;
                for( k=k; k<nEnd; k++ )
                {
                    nSlope = t[k+1] - t[k];
                    if( nSlope < 0 )
                    {
                        if( ++nCtr >= nSlopeCount )
                        {
                            nState = STATE_PEAK;
                            break;
                        }
                    }
                    else
                    {
                        nState = STATE_START;
                        break;
                    }
                }
                break;



            case STATE_PEAK:
                // Peak found
                nResume = k + 1;
                return nPeak;



            case STATE_NOPEAK:
                // Peak not found
                nResume = nEnd + 1;
                return -1;

        }
    }
}



int Trace::PosPeakFindLargest( int n, int nBegin, int nEnd, int& nCount, int nSlopeCount ) const
{
/*
    Searches for a positive peak on trace 'n' beginning at 'nStart' and
    continuing until a peak is found or we reach 'nEnd'. Returns the position
    of the largest peak if one is found, otherwise -1. The number of peaks
    between the two limits is also returned in 'nCount'.
*/
    int nCurPos;
    int nCurPeak;
    int nResume;
    int nMaxPos  = -1;
    int nMaxPeak = INT_MIN;
    nCount = 0;
    while(1)
    {
        nCurPos = PosPeakFind( n, nBegin, nEnd, nResume, nSlopeCount );
        if( nCurPos < 0 )
            break;
        nCurPeak = m_pTrace[n][nCurPos];
        if( nCurPeak > nMaxPeak )
        {
            nMaxPos  = nCurPos;
            nMaxPeak = nCurPeak;
        }
        nBegin = nResume;
        nCount++;
    }
    return nMaxPos;
}




//--------------------
// Negative Peak Find
//--------------------

int Trace::NegPeakFind( int n, int nBegin, int nEnd, int& nResume, int nSlopeCount ) const
{
/*
    Scans trace channel 'n' from 'nBegin' to 'nEnd' looking for a well defined
    negative peak. A peak is defined as nSlopeCount consecutive -ve slopes,
    followed by zero or more occurances of a zero slope, followed by nSlopeCount
    consecutive +ve slopes. This ensures that only well defined peaks are
    found. Returns the peak position on success, otherwise -1.
*/
    assert(n<4);
    enum { STATE_START, STATE_CLIMB, STATE_FLAT, STATE_DESCEND, STATE_PEAK, STATE_NOPEAK };
    int    nCtr;
    int    nSlope;
    int    nZeros = 0;
    int    nPeak = 0;
    int    k      = nBegin;
    TRACE* t      = m_pTrace[n];
    int    nState = STATE_START;



    // Search loop
    while(1)
    {
        switch(nState)
        {
            case STATE_START:
                // Search for an initial downwards slope
                nState = STATE_NOPEAK;
                for( k=k; k<nEnd; k++ )
                {
                    nSlope = t[k+1] - t[k];
                    if( nSlope < 0 )
                    {
                        nState = STATE_DESCEND;
                        break;
                    }
                }
                break;



            case STATE_DESCEND:
                // Descend down a falling slope
                nCtr   = 0;
                nZeros = 0;
                nState = STATE_NOPEAK;
                for( k=k; k<nEnd; k++ )
                {
                    nSlope = t[k+1] - t[k];
                    if( nSlope < 0 )
                        nCtr++;
                    else if( (nSlope==0) && (nCtr>=nSlopeCount) )
                    {
                        nState = STATE_FLAT;
                        break;
                    }
                    else if( (nSlope>0) && (nCtr>=nSlopeCount) )
                    {
                        nPeak  = k;
                        nState = STATE_PEAK;
                        break;
                    }
                    else
                    {
                        nState = STATE_START;
                        break;
                    }
                }
                break;



            case STATE_FLAT:
                // Traverse clipped peak. Allow for a bit
                // of noise on this peak.
                nState = STATE_NOPEAK;
                for( k=k; k<nEnd; k++ )
                {
                    nSlope = t[k+1] - t[k];
                    if( (nSlope>-3) && (nSlope<3) )
                        nZeros++;
                    else if( nSlope < 0 )
                    {
                        nState = STATE_DESCEND;
                        break;
                    }
                    else
                    {
                        nPeak  = k - (nZeros/2);
                        nState = STATE_CLIMB;
                        break;
                    }
                }
                break;



            case STATE_CLIMB:
                // Climb up a rising slope
                nCtr   = 0;
                nState = STATE_NOPEAK;
                for( k=k; k<nEnd; k++ )
                {
                    nSlope = t[k+1] - t[k];
                    if( nSlope > 0 )
                    {
                        if( ++nCtr >= nSlopeCount )
                        {
                            nState = STATE_PEAK;
                            break;
                        }
                    }
                    else
                    {
                        nState = STATE_START;
                        break;
                    }
                }
                break;



            case STATE_PEAK:
                // Peak found
                nResume = k + 1;
                return nPeak;



            case STATE_NOPEAK:
                // Peak not found
                nResume = nEnd + 1;
                return -1;

        }
    }
}



int Trace::NegPeakFindLargest( int n, int nBegin, int nEnd, int& nCount, int nSlopeCount ) const
{
/*
    Searches for a negative peak on trace 'n' beginning at 'nStart' and
    continuing until a peak is found or we reach 'nEnd'. Returns the position
    of the largest peak if one is found, otherwise -1. The number of peaks
    between the two limits is also returned in 'nCount'.
*/
    int nCurPos;
    int nCurPeak;
    int nResume;
    int nMaxPos  = -1;
    int nMaxPeak = INT_MAX;
    nCount = 0;
    while(1)
    {
        nCurPos = NegPeakFind( n, nBegin, nEnd, nResume, nSlopeCount );
        if( nCurPos < 0 )
            break;
        nCurPeak = m_pTrace[n][nCurPos];
        if( nCurPeak < nMaxPeak )
        {
            nMaxPos  = nCurPos;
            nMaxPeak = nCurPeak;
        }
        nBegin = nResume;
        nCount++;
    }
    return nMaxPos;
}



//-------
// Clone
//-------

Trace* Trace::Clone( const char* pNewFileName ) const
{
/*
   Creates a copy of this trace object with optionally a new filename. The caller
   is responsible for deleting the returned trace object.
*/
   Read*  r = 0;
   Trace* p = 0;
   try
   {
        // Duplicate the trace and wrap it up
        r = read_dup( m_pRead, pNewFileName );
        if( !r )
            throw std::bad_alloc();
        p = new Trace;
        p->Wrap( r, true );

   }
   catch(...)
   {
      // Tidy up before rethrowing
      if(r) read_deallocate(r);
      if(p) delete p;
      throw;
   }
   return p;
}




//----------------
// CreateEnvelope
//----------------

Trace* Trace::CreateEnvelope() const
{
/*
   Creates a new trace containing the envelope and bases of this trace and
   returns it to the caller. The envelope is stored in the 'A' trace channel.
   The caller is responsible for deleting the envelope.
*/
   TRACE Max1;
   TRACE Max2;
   Trace& Envelope = *Clone();
   if( &Envelope )
   {
      for( int n=0; n<Envelope.Samples(); n++ )
      {
         Max1 = std::max( Envelope[0][n], Envelope[1][n] );
         Max2 = std::max( Envelope[2][n], Envelope[3][n] );
         Envelope[0][n] = std::max( Max1, Max2 );
         Envelope[1][n] = 0;
         Envelope[2][n] = 0;
         Envelope[3][n] = 0;
      }
   }
   return &Envelope;
}



//-------------------
// Update Statistics
//-------------------

extern "C" {

int TraceCompareIntegers( const void* a, const void* b )
{
    // For qsort
    return *((const int*)a) - *((const int*)b);
}

}


void Trace::UpdateStatistics()
{
/*
    Computes the interval statistics of the trace.
*/
    // Already done?
    if( m_bStatisticsCached )
        return;



    // Create an array of base positions
    int n;
    int nBases = Range();
    NumericArray<int> Interval;
    Interval.Create( nBases );
    for( n=0; n<nBases; n++ )
        Interval[n] = int( m_pRead->basePos[n+m_nLowerLimit] );



    // Sort the array, phred can call bases out of order
    std::qsort( Interval.Raw(), nBases, sizeof(int), TraceCompareIntegers );



    // Replace the contents with intervals
    nBases--;
    for( n=0; n<nBases; n++ )
        Interval[n] = Interval[n+1] - Interval[n];
    Interval.Length( nBases );



    // Sort the array again for the mode computation
    std::qsort( Interval.Raw(), nBases, sizeof(int), TraceCompareIntegers );



    // Compute the interval min/max/mean/sd
    m_nIntervalMin    = Interval.Min();
    m_nIntervalMax    = Interval.Max();
    m_nIntervalMean   = Interval.Mean();
    m_nIntervalStdDev = Interval.StandardDeviation( &m_nIntervalMean );



    // Compute the mode by counting the maximum run length
    int nRun        =  0;
    int nRunMax     =  0;
    int nLast       = -1;
    m_nIntervalMode =  0;
    for( n=0; n<nBases; n++ )
    {
        if( Interval[n] == nLast )
            nRun++;
        else
        {
            if( nRun > nRunMax )
            {
                m_nIntervalMode = nLast;
                nRunMax = nRun;
            }
            nRun  = 1;
            nLast = Interval[n];
        }
    }
    m_bStatisticsCached = true;
}



//-------
// Range
//-------

void Trace::Range( int n1, int n2 )
{
    // Range is inclusive of n1 and n2
    assert(n1>=0);
    assert(n1<=n2);
    m_nLowerLimit = n1;
    m_nUpperLimit = n2;
}



//----------
// Subtract
//----------

Trace* Trace::Subtract( Trace& t )
{
/*
    Subtracts trace 't' from this trace object and returns a new trace
    object containing the result. The difference is level shifted for
    optimum display in trev and gap4. The two traces must be the same
    length in samples, otherwise a NULL pointer is returned.
*/
    assert(m_pRead!=0);



    // Check input
    assert(Samples()==t.Samples());
    if( Samples() != t.Samples() )
        return 0;



    // Clone ourself
    Trace* pDifference = Clone( "difference" );
    if( !pDifference )
        return 0;



    // Determine scaling by looking at maxval's. It may be necessary
    // to scale down the trace because some sequencing machines
    // (eg Beckman) generate traces with the full 16-bit resolution.
    double nScale = 1.0;
    int    nMax   = t.Max() >= Max() ? t.Max() : Max();
    if( nMax >= 16384 )
    {
        nScale = 0.5;
        nMax   = nMax / 2;
    }



    // Generate the difference trace centred around nMax
    for( int k=0, d=0; k<Samples(); k++ )
    {
        for( int j=0; j<4; j++ )
        {
            d = static_cast<int>(m_pTrace[j][k]) - static_cast<int>(t[j][k]);
            d = static_cast<int>( d * nScale ) + nMax;
            (*pDifference)[j][k] = static_cast<TRACE>( d );
        }
    }



    // Update metadata
    pDifference->Raw()->baseline    = nMax;
    pDifference->Raw()->maxTraceVal = 2 * nMax;
    pDifference->Raw()->leftCutoff  = 0;
    pDifference->Raw()->rightCutoff = 0;

    return pDifference;
}


void Trace::ScaleTo( Trace& t )
{
/*
    Scales this trace in a pointwise fashion based on the energy at each
    sample position to that of the given trace 't'. The two traces must
    be the same length in samples, otherwise no scaling is done.
*/
    assert(m_pRead!=0);
    assert(Samples()==t.Samples());


    // Check input
    int nLen = Samples();
    if( t.Samples() != nLen )
        return;


    // Create scale factor array
    NumericArray<double> Scale;
    Scale.Create( nLen );


    // Compute pointwise scale factors using signal energy measure
    double e1, e2;
    double last_good_sf = 1.0;
    for( int n=0; n<nLen; n++ )
    {
        e1  = m_pTrace[0][n];
        e1 += m_pTrace[1][n];
        e1 += m_pTrace[2][n];
        e1 += m_pTrace[3][n];
        e2  = t[0][n];
        e2 += t[1][n];
        e2 += t[2][n];
        e2 += t[3][n];
        if( e1 == 0.0 )
            Scale[n] = last_good_sf;
        else
        {
            Scale[n]     = e2 / e1;
            last_good_sf = Scale[n];
        }
    }

#if 1
    // Scale this trace
    for( int n=0; n<nLen; n++ )
    {
        m_pTrace[0][n] = static_cast<TRACE>( m_pTrace[0][n] * Scale[n] );
        m_pTrace[1][n] = static_cast<TRACE>( m_pTrace[1][n] * Scale[n] );
        m_pTrace[2][n] = static_cast<TRACE>( m_pTrace[2][n] * Scale[n] );
        m_pTrace[3][n] = static_cast<TRACE>( m_pTrace[3][n] * Scale[n] );
    }
#else
    // Double this trace - useful for debugging tracediff on heterozygous
    // indels.
    for( int n=0; n<nLen; n++ )
    {
        m_pTrace[0][n] = static_cast<TRACE>( m_pTrace[0][n] * 2 );
        m_pTrace[1][n] = static_cast<TRACE>( m_pTrace[1][n] * 2 );
        m_pTrace[2][n] = static_cast<TRACE>( m_pTrace[2][n] * 2 );
        m_pTrace[3][n] = static_cast<TRACE>( m_pTrace[3][n] * 2 );
    }
#endif
}



//------
// Mean
//------

double Trace::Mean( int n ) const
{
/*
    Computes the mean value of trace 'n' over its entire length. If 'n'
    is less than zero, the mean of all traces combined is computed.
*/
    assert(m_pRead!=0);
    int    k;
    double nMeanVal = 0.0;
    if( n < 0 )
    {
        // All traces
        for( k=0; k<Samples(); k++ )
        {
            nMeanVal += m_pTrace[0][k] + m_pTrace[1][k] +
                        m_pTrace[2][k] + m_pTrace[3][k];
        }
        nMeanVal = nMeanVal / double( 4*Samples() );
    }
    else
    {
        // Single trace
        for( k=0; k<Samples(); k++ )
            nMeanVal += m_pTrace[n][k];
        nMeanVal = nMeanVal / double( Samples() );
    }
    return nMeanVal;
}



//------
// Sort
//------

typedef struct
{
   char cBase;
   int  nBasePos;
   char cProbability[4];

}BASECALL;


extern "C"
int TraceCompareBasePositions( const void *elem1, const void *elem2 )
{
   const BASECALL* p1 = static_cast<const BASECALL*>( elem1 );
   const BASECALL* p2 = static_cast<const BASECALL*>( elem2 );
   return p1->nBasePos - p2->nBasePos;
}


void Trace::Sort()
{
/*
    Phred can swap the order of basecalls which means that the base positions
    are not always monotonic. Here we sort all bases in order of base position.
    This may result in an incorrect sequence (according to Phred), but there's
    not much else we can do without resorting to cutting out bits of trace
    and shifting them around.
*/

   // Initialisation
   assert(m_pRead!=0);
   int n;
   int nBases = Bases();



   // Copy basecall data into an array
   SimpleArray<BASECALL> Data;
   Data.Create( nBases );
   for( n=0; n<nBases; n++ )
   {
      Data[n].cBase    = m_pRead->base[n];
      Data[n].nBasePos = static_cast<int>( m_pRead->basePos[n] );
      if( m_pRead->prob_A !=0 )
      {
         Data[n].cProbability[0] = m_pRead->prob_A[n];
         Data[n].cProbability[1] = m_pRead->prob_C[n];
         Data[n].cProbability[2] = m_pRead->prob_G[n];
         Data[n].cProbability[3] = m_pRead->prob_T[n];
      }
   }



   // Sort the array by base position
   std::qsort( Data.Raw(), nBases, sizeof(BASECALL), TraceCompareBasePositions );



   // Copy basecall data back into the read structure
   for( n=0; n<nBases; n++ )
   {
      m_pRead->base[n]    = Data[n].cBase;
      m_pRead->basePos[n] = static_cast<uint_2>( Data[n].nBasePos );
      if( m_pRead->prob_A != 0 )
      {
         m_pRead->prob_A[n]  = Data[n].cProbability[0];
         m_pRead->prob_C[n]  = Data[n].cProbability[1];
         m_pRead->prob_G[n]  = Data[n].cProbability[2];
         m_pRead->prob_T[n]  = Data[n].cProbability[3];
      }
   }
}



int Trace::BaseNumberFromSamplePosition( int nPosition ) const
{
/*
    Determines the base number from sample position and rounds to the
    nearest base. The baseno is zero based, so you may need to add 1
    to obtain the "end users" base number.
*/
    assert(m_pRead!=0);
    assert(nPosition>=0);
    assert(nPosition<Samples());


    // Check input
    if( (nPosition<0) || (nPosition>=Samples()) )
        return 0;


    // Find approx base number for mutation, assumes bases are in order!
    int b;
    int nBases = Bases() - 1;
    for( b=0; b<nBases; b++ )
    {
       if( BasePosition(b) >= nPosition )
            break;
    }


    // Get the neighbouring base positions
    int bp[2] = { 0, 0 };
    if( b > 0 )
       bp[0] = BasePosition( b-1 );
    bp[1] = BasePosition( b );



    // Return the nearest one
    int bd[2];
    bd[0] = std::abs( nPosition - bp[0] );
    bd[1] = std::abs( bp[1] - nPosition );
    bp[0] = (bd[0]<bd[1]) ? b-1 : b;
    return (bp[0]<0) ? 0 : bp[0];
}



void Trace::SetBase( int n, char b, int pos, int conf )
{
   assert(m_pRead!=0);
   assert(n>=0);
   assert(n<m_pRead->NBases);
   char confidence     = static_cast<char>(conf);
   m_pRead->base[n]    = b;
   m_pRead->basePos[n] = static_cast<uint_2>(pos);
   m_pRead->prob_A[n]  = 0;
   m_pRead->prob_C[n]  = 0;
   m_pRead->prob_G[n]  = 0;
   m_pRead->prob_T[n]  = 0;
   switch( b )
   {
       case 'A':
       case 'a':  m_pRead->prob_A[n] = confidence;  break;

       case 'C':
       case 'c':  m_pRead->prob_C[n] = confidence;  break;

       case 'G':
       case 'g':  m_pRead->prob_G[n] = confidence;  break;

       case 'T':
       case 't':  m_pRead->prob_T[n] = confidence;  break;

       default:   m_pRead->prob_A[n] = confidence;
                  m_pRead->prob_C[n] = confidence;
                  m_pRead->prob_G[n] = confidence;
                  m_pRead->prob_T[n] = confidence;  break;
   }
}


void Trace::Floor( int threshold )
{
    int signal;
    const int cols     = Samples();
    const int baseline = Baseline();
    puts("floor");
    for( int r=0; r<4; r++ )
    {
        for( int c=0; c<cols; c++ )
        {
            signal = std::abs( m_pTrace[r][c] - baseline );
            if( signal < threshold )
                m_pTrace[r][c] = baseline;
        }
    }
}


void Trace::FloorHalfwaves()
{
    const int cols     = Samples();
    const int baseline = Baseline();
    puts("floorhalfwaves");
    for( int c=0; c<cols; c++ )
    {
        // Keep record of number of positive/negative signal values
        int pos = 0;
        int neg = 0;
        for( int r=0; r<4; r++ )
        {
            if( m_pTrace[r][c] != baseline )
            {
                if( (m_pTrace[r][c]-baseline) < 0 )
                    neg++;
                else
                    pos++;
            }
        }
        // If unidirectional, floor it
        if( (neg==0) || (pos==0) )
        {
            m_pTrace[0][c] = baseline;
            m_pTrace[1][c] = baseline;
            m_pTrace[2][c] = baseline;
            m_pTrace[3][c] = baseline;
        }
    }
}


void Trace::FloorNarrowPeaks( int width_threshold )
{
    int pos;
    int left;
    int right;
    int width;
    const int cols     = Samples();
    const int baseline = Baseline();


    puts("floornarrowpeaks");
    for( int r=0; r<4; r++ )
    {
        int resume = 0;
        while(1)
        {
            // Find a positive peak
            pos = PosPeakFind( r, resume, cols-1, resume, 1 );
            if( pos < 0 )
                break;


            // Measure it's width, if less than the threshold, floor that peak
            width = PosPeakWidth( r, pos, left, right, baseline );
            if( width < width_threshold )
            {
                for( int n=left; n<=right; n++ )
                    m_pTrace[r][n] = baseline;
            }
        }
    }


    // Reapply halfwave flooring to remove any negative peaks
    FloorHalfwaves();
}



void Trace::Smooth()
{
    // Apply 3-term moving average
    const int end = Samples() - 1;
    for( int r=0; r<4; r++ )
    {
        if( end >= 2 )
        {
            for( int c=1; c<end; c++ )
                m_pTrace[r][c] = (m_pTrace[r][c-1]+m_pTrace[r][c]+m_pTrace[r][c+1]) / 3;
        }
    }
}


void Trace::FillGaps()
{
    // Fill 1-sample baseline gaps with 3-term average
    const int baseline = Baseline();
    const int end      = Samples() - 1;
    for( int r=0; r<4; r++ )
    {
        if( end >= 2 )
        {
            for( int c=1; c<end; c++ )
            {
                if( (m_pTrace[r][c]==baseline) && (m_pTrace[r][c-1]!=baseline) && (m_pTrace[r][c+1]!=baseline) )
                    m_pTrace[r][c] = (m_pTrace[r][c-1]+m_pTrace[r][c]+m_pTrace[r][c+1]) / 3;
            }
        }
    }
#if 0
    {
      double tot = 0, avg;
      for (int i = 0; i < Samples(); i++) {
	for (int j = 0; j < 4; j++) {
	  tot += ABS(m_pTrace[j][i]-baseline);
	}
      }
      avg = tot / Samples();
      printf("Avg = %f\n", avg);

#define WL 20
      for (int i = Samples()-11; i >= WL; i--) {
	double t = 0;
	for (int x = -WL; x <= WL; x++) {
	  tot = 0;
	  for (int j = 0; j < 4; j++) {
	    tot += ABS(m_pTrace[j][i+x]-baseline);
	  }
	  t += tot;
	}
	t /= 2*WL+1;
	printf("%d %d %f %f\n", i, BaseNumberFromSamplePosition(i), t, t/avg);
      }
    }
#endif

}

void Trace::AvgFilt( double scale ) {
    const int cols     = Samples();
    const int baseline = Baseline();
    double totp = 0, totm = 0;

    for( int c=0; c<cols; c++ ) {
	totp *= 0.98;
	totm *= 0.98;
	double ratio;

	for( int r=0; r<4; r++ ) {
	    if (m_pTrace[r][c] > baseline)
		totp += m_pTrace[r][c] - baseline;
	    else
		totm += baseline - m_pTrace[r][c];
	}

	ratio = (totp+1) / (totm+1);
	if (ratio < 1)
	    ratio = 1/ratio;

	printf("%d %f %f %f %d\n", c, totp, totm, ratio, baseline/2);

	if (ratio > 20 || (totp > baseline*2 && totm > baseline*2)) {
	    for( int r=0; r<4; r++ )
		m_pTrace[r][c] = baseline;
	}
    }    
}
