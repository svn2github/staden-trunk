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


#include <climits>             // For INT_MIN, INT_MAX
#include <peakcall.hpp>



//------
// Init
//------

void PeakCall::Init()
{
    for( int n=0; n<4; n++ )
    {
        Data.Width[n]     = -1;
        Data.Position[n]  = -1;
        Data.Amplitude[n] = -1;
    }
    Data.BaseNumber   = -1;
    Data.BasePosition = -1;
}



//-------------
// Constructor
//-------------

PeakCall::PeakCall( int a, int c, int g, int t, int bnum, int bpos )
{
    Data.Amplitude[0] = a;
    Data.Amplitude[1] = c;
    Data.Amplitude[2] = g;
    Data.Amplitude[3] = t;
    Data.BaseNumber   = bnum;
    Data.BasePosition = bpos;
}



//---------
// IsValid
//---------

bool PeakCall::IsValid() const
{
   const int* Amp = Data.Amplitude;
   if( (Amp[0]!=-1) || (Amp[1]!=-1) || (Amp[2]!=-1) || (Amp[3]!=-1) )
       return true;
   return false;
}



//----------
// MaxWidth
//----------

int PeakCall::MaxWidthAsIndex() const
{
    int n    = -1;
    int nMax = INT_MIN;
    if( (Data.Position[0]!=-1) && (Data.Width[0]>nMax) )
    {
        nMax = Data.Width[0];
        n    = 0;
    }
    if( (Data.Position[1]!=-1) && (Data.Width[1]>nMax) )
    {
        nMax = Data.Width[1];
        n    = 1;
    }
    if( (Data.Position[2]!=-1) && (Data.Width[2]>nMax) )
    {
        nMax = Data.Width[2];
        n    = 2;
    }
    if( (Data.Position[3]!=-1) && (Data.Width[3]>nMax) )
        n = 3;
    return n;
}



//--------------
// MaxAmplitude
//--------------

int PeakCall::MaxAmplitudeAsIndex() const
{
    int n    = -1;
    int nMax = INT_MIN;
    if( (Data.Position[0]!=-1) && (Data.Amplitude[0]>nMax) )
    {
        nMax = Data.Amplitude[0];
        n    = 0;
    }
    if( (Data.Position[1]!=-1) && (Data.Amplitude[1]>nMax) )
    {
        nMax = Data.Amplitude[1];
        n    = 1;
    }
    if( (Data.Position[2]!=-1) && (Data.Amplitude[2]>nMax) )
    {
        nMax = Data.Amplitude[2];
        n    = 2;
    }
    if( (Data.Position[3]!=-1) && (Data.Amplitude[3]>nMax) )
        n = 3;
    return n;
}



//--------------
// MinAmplitude
//--------------

int PeakCall::MinAmplitudeAsIndex() const
{
    int n    = -1;
    int nMin = INT_MAX;
    if( (Data.Position[0]!=-1) && (Data.Amplitude[0]<nMin) )
    {
        nMin = Data.Amplitude[0];
        n    = 0;
    }
    if( (Data.Position[1]!=-1) && (Data.Amplitude[1]<nMin) )
    {
        nMin = Data.Amplitude[1];
        n    = 1;
    }
    if( (Data.Position[2]!=-1) && (Data.Amplitude[2]<nMin) )
    {
        nMin = Data.Amplitude[2];
        n    = 2;
    }
    if( (Data.Position[3]!=-1) && (Data.Amplitude[3]<nMin) )
        n = 3;
    return n;
}

