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



#ifndef _MUTLIB_TRACE_HPP_
#define _MUTLIB_TRACE_HPP_



#include <cassert>
#include <staden.h>     // For Read* structure



class Trace
{
 public:
   // Constructors/destructor
   Trace()                                      { Init(); }
  ~Trace()                                      { Close(); }
   Trace( Read* r, bool AutoDestroy=true )      { Wrap(r,AutoDestroy); }



 public:
    // Services
    void        Sort();
    void        Close();
    void        Smooth();
    void        FillGaps();
    double      Mean( int n=-1 ) const;
    Trace*      CreateEnvelope() const;
    Trace*      Subtract( Trace& t );
    void        ScaleTo( Trace& t );
    void        FloorHalfwaves();
    void        Floor( int threshold );
    void        AvgFilt( double scale );
    void        FloorNarrowPeaks( int width_threshold );
    void        Range( int n1, int n2 );
    bool        Open( const char* pFileName );
    void        MaxAt( int n, int& c, int& m ) const;
    void        MinAt( int n, int& c, int& m ) const;
    bool        SaveAs( const char* pFileName, int fmt=TT_ZTR );
    bool        Create( int nSamples, int nBases, const char* pName=0 );
    void        Wrap( Read* r, bool AutoDestroy=true );
    Trace*      Clone( const char* pNewFileName=0 ) const;
    int         BaseNumberFromSamplePosition( int nPosition ) const;
    int         PosPeakFind( int n, int nBegin, int nEnd, int& nResume, int nSlopeCount=1 ) const;
    int         NegPeakFind( int n, int nBegin, int nEnd, int& nResume, int nSlopeCount=1 ) const;
    int         PosPeakFindLargest( int n, int nBegin, int nEnd, int& nCount, int nSlopeCount=1 ) const;
    int         NegPeakFindLargest( int n, int nBegin, int nEnd, int& nCount, int nSlopeCount=1 ) const;
    int         PosPeakWidth( int n, int nPos, int& nLeft, int& nRight ) const;
    int         PosPeakWidth( int n, int nPos, int& nLeft, int& nRight, int nMeasurementThreshold ) const;
    int         NegPeakWidth( int n, int nPos, int& nLeft, int& nRight, int nMeasurementThreshold ) const;
    void        WindowToLeftOf( int nPosition, int nSize, int& l, int& r ) const;
    void        WindowCentredAt( int nPosition, int nSize, int& l, int& r ) const;
    int         IntervalMin()                       { assert(m_pRead!=0); UpdateStatistics(); return m_nIntervalMin; }
    int         IntervalMax()                       { assert(m_pRead!=0); UpdateStatistics(); return m_nIntervalMax; }
    int         IntervalMode()                      { assert(m_pRead!=0); UpdateStatistics(); return m_nIntervalMode; }
    double      IntervalMean()                      { assert(m_pRead!=0); UpdateStatistics(); return m_nIntervalMean; }
    double      IntervalStdDev()                    { assert(m_pRead!=0); UpdateStatistics(); return m_nIntervalStdDev; }
    void        InvalidateStatistics()              { m_bStatisticsCached=false; }
    bool        AutoDestroy() const                 { return m_bAutoDestroy; }
    void        AutoDestroy( bool s )               { m_bAutoDestroy=s; }
    int         Range() const                       { return m_nUpperLimit-m_nLowerLimit+1; }
    Read*       Raw() const                         { assert(m_pRead!=0); return m_pRead; }
    int         Max() const                         { assert(m_pRead!=0); return m_pRead->maxTraceVal; }
    void        Max( int n )                        { assert(m_pRead!=0); m_pRead->maxTraceVal = n; }
    const char* Name() const                        { assert(m_pRead!=0); return m_pRead->trace_name ? m_pRead->trace_name : ""; }
    int         Bases() const                       { assert(m_pRead!=0); return m_pRead->NBases; }
    void        Bases( int n )                      { assert(m_pRead!=0); m_pRead->NBases=n; }
    int         Samples() const                     { assert(m_pRead!=0); return m_pRead->NPoints; }
    void        Samples( int n )                    { assert(m_pRead!=0); m_pRead->NPoints=n; }
    double      SamplesPerBase() const              { assert(m_pRead!=0); assert(m_pRead->NBases>0); return double(m_pRead->NPoints) / double(m_pRead->NBases); }
    int         Baseline() const                    { assert(m_pRead!=0); return m_pRead->baseline; }
    char        BaseChar( int n ) const             { assert(m_pRead!=0); assert(n>=0); assert(n<m_pRead->NBases); return m_pRead->base[n]; }
    int         BasePosition( int n ) const         { assert(m_pRead!=0); assert(n>=0); assert(n<m_pRead->NBases); return m_pRead->basePos[n]; }
    int         BaseConfidence( int n ) const;
    void        SetBase( int n, char b, int pos, int conf );



 public:
    // Operators
    TRACE* operator[]( int n )                      { assert(n<4); assert(n>=0); return m_pTrace[n]; }



 private:
    // Helpers
    void Init();
    void InitTraces();
    void ZeroTraces();
    void UpdateStatistics();



 private:
    // Data
    Read*  m_pRead;
    TRACE* m_pTrace[4];                         // Required due to poor Read structure design
    int    m_nLowerLimit;
    int    m_nUpperLimit;
    bool   m_bAutoDestroy;
    bool   m_bExternalTrace;
    int    m_nIntervalMin;
    int    m_nIntervalMax;
    int    m_nIntervalMode;
    double m_nIntervalMean;
    double m_nIntervalStdDev;
    bool   m_bStatisticsCached;
};




#endif

