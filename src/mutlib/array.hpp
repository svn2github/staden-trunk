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


#ifndef _MUTLIB_ARRAY_HPP_
#define _MUTLIB_ARRAY_HPP_


#include <cmath>             // For sqrt()
#include <limits>            // For std::numeric_limits<T>
#include <cstring>           // For memcpy()
#include <cassert>
#include <fstream>           // For ofstream
#include <iomanip>           // For setw()


template <typename T>
class SimpleArray
{
 public:
   // Constructor/Destructor
   SimpleArray()                       { Init(); }
  ~SimpleArray()                       { Empty(); }


 public:
   // Services
   T*   Raw()                          { return m_pArray; }
   void Empty();
   void Fill( T v );
   int  Length() const                 { return m_nLength; }
   void Length( int n );
   int  Capacity() const               { return m_nCapacity; }
   void Range( int n1, int n2 );
   int  Range() const                  { return m_nUpperLimit-m_nLowerLimit+1; }
   int  RangeLowerLimit() const        { return m_nLowerLimit; }
   int  RangeUpperLimit() const        { return m_nUpperLimit; }
   void RangeReset()                   { m_nLowerLimit=0; m_nUpperLimit=m_nLength-1; }
   void AutoDestroy( bool s )          { m_bAutoDestroy=s; }
   void Create( int nCapacity );
   void Create( T* p, int nLength );
   void Wrap( T* p, int nCapacity, bool AutoDestroy=false );
   void SaveAs( const char* pFileName, bool bAppend=false, bool bAsInteger=false );
   void Append( T* p, int n );



 public:
   // Operators
   const T& operator[]( int n ) const  { assert(n<m_nCapacity); return m_pArray[n]; }
   T&       operator[]( int n )        { assert(n<m_nCapacity); return m_pArray[n]; }



 protected:
   T*   m_pArray;
   int  m_nLength;
   int  m_nCapacity;
   bool m_bAutoDestroy;
   int  m_nLowerLimit;
   int  m_nUpperLimit;



 private:
   // Helpers
   void Init();
};



//------
// Init
//------

template <typename T>
void SimpleArray<T>::Init()
{
   m_pArray       = 0;
   m_nLength      = 0;
   m_nCapacity    = 0;
   m_nLowerLimit  = 0;
   m_nUpperLimit  = 0;
   m_bAutoDestroy = true;
}



//------
// Fill
//------

template <typename T>
void SimpleArray<T>::Fill( T v )
{
   for( int n=m_nLowerLimit; n<=m_nUpperLimit; n++ )
      m_pArray[n] = v;
}



//--------
// Create
//--------

template <typename T>
void SimpleArray<T>::Create( int nCapacity )
{
   assert(nCapacity>0);
   if( m_pArray )
      Empty();
   m_pArray       = new T[ nCapacity ];
   m_nCapacity    = nCapacity;
   m_nLength      = nCapacity;
   m_nLowerLimit  = 0;
   m_nUpperLimit  = nCapacity - 1;
   m_bAutoDestroy = true;
}


template <typename T>
void SimpleArray<T>::Create( T* p, int nLength )
{
   assert(nLength>0);
   if( m_pArray )
      Empty();
   m_pArray       = new T[ nLength ];
   std::memcpy( m_pArray, p, sizeof(T)*nLength );
   m_nCapacity    = nLength;
   m_nLength      = nLength;
   m_nLowerLimit  = 0;
   m_nUpperLimit  = nLength - 1;
   m_bAutoDestroy = true;
}




//--------
// Length
//--------

template <typename T>
void SimpleArray<T>::Length( int n )
{
    assert(n<=m_nCapacity);
    m_nLength     = n;
    m_nLowerLimit = 0;
    m_nUpperLimit = n-1;
}




//------
// Wrap
//------

template <typename T>
void SimpleArray<T>::Wrap( T* p, int nCapacity, bool AutoDestroy )
{
   assert(p != NULL);
   assert(nCapacity>0);
   if( m_pArray )
      Empty();
   m_pArray       = p;
   m_nLength      = nCapacity;
   m_nCapacity    = nCapacity;
   m_nLowerLimit  = 0;
   m_nUpperLimit  = nCapacity - 1;
   m_bAutoDestroy = AutoDestroy;
}




//-------
// Empty
//-------

template <typename T>
void SimpleArray<T>::Empty()
{
   if( m_bAutoDestroy == true )
      delete [] m_pArray;
   Init();
}



//-------
// Range
//-------

template <typename T>
void SimpleArray<T>::Range( int n1, int n2 )
{
    // Range is inclusive of n1 and n2
    assert(n1>=0);
    assert(n1<=n2);
    assert(n1<m_nLength);
    assert(n2<m_nLength);
    m_nLowerLimit = n1;
    m_nUpperLimit = n2;
}



//--------
// Append
//--------

template <typename T>
void SimpleArray<T>::Append( T* p, int n )
{
    assert(p != NULL);
    assert(n>0);


    // Check input
    if( !p || (n<=0) )
        return;


    // Allocate a new buffer if necessary
    int nSpaceLeft = m_nCapacity - m_nLength;
    if( nSpaceLeft < n )
    {
        // Reallocate buffer, copy old elements
        const int CHUNK_SIZE = 256;
        int nElementsToAllocate = m_nLength;
        nElementsToAllocate += (n<CHUNK_SIZE) ? CHUNK_SIZE : n;
        T* pNewBuffer = new T[ nElementsToAllocate ];
        for( int k=0; k<m_nLength; k++ )
            pNewBuffer[k] = m_pArray[k];
        delete [] m_pArray;
        m_pArray    = pNewBuffer;
        m_nCapacity = nElementsToAllocate;
    }


    // Copy in new elements
    for( int k=0, d=m_nLength; k<n; k++, d++ )
        m_pArray[d] = p[k];
    m_nLength += n;
}




//---------
// Save As
//---------

template <typename T>
void SimpleArray<T>::SaveAs( const char* pFileName, bool bAppend, bool bAsInteger )
{
#ifdef _WIN32
    using namespace std;
    ofstream ofs;
    if( bAppend )
        ofs.open( pFileName, ios_base::out | ios_base::app );
    else
        ofs.open( pFileName, ios_base::out | ios_base::trunc );
    if( bAsInteger )
    {
        for( int n=m_nLowerLimit; n<=m_nUpperLimit; n++ )
            ofs << setw(4) << static_cast<int>(m_pArray[n]);
    }
    else
    {
        for( int n=m_nLowerLimit; n<=m_nUpperLimit; n++ )
            ofs << m_pArray[n] << ' ';
    }
    ofs << endl;
#endif
}




//--------------------------------- Numeric Array -----------------------------



template <typename T>
class NumericArray : public SimpleArray<T>
{
 public:
   // Services
   T      Min() const;
   T      Max() const;
   double Mean() const;
   void   Interpolate( int x1, int x2 );
   double Variance( double* nMean=0 ) const;
   double StandardDeviation( double* nMean=0 ) const    { return std::sqrt( Variance(nMean) ); }
};



//------
// Mean
//------

template <typename T>
double NumericArray<T>::Mean() const
{
   assert(this->m_pArray != NULL);
   double acc = 0.0;
   for( int n=this->m_nLowerLimit; n<=this->m_nUpperLimit; n++ )
      acc += this->m_pArray[n];
   return (this->m_nUpperLimit - this->m_nLowerLimit+1)
     ? acc / (this->m_nUpperLimit - this->m_nLowerLimit+1)
     : 0;
}



//----------
// Variance
//----------

template <typename T>
double NumericArray<T>::Variance( double* nMean ) const
{
   assert(this->m_pArray != NULL);
   double tmp;
   double acc  = 0.0;
   double mean = nMean ? *nMean : Mean();
   for( int n=this->m_nLowerLimit; n<=this->m_nUpperLimit; n++ )
   {
      tmp  = double(this->m_pArray[n]) - mean;
      tmp *= tmp;
      acc += tmp;
   }
   assert(this->m_nUpperLimit - this->m_nLowerLimit!=0);
   return acc / static_cast<double>(this->m_nUpperLimit - this->m_nLowerLimit);
}



//-----
// Min
//-----

template <typename T>
T NumericArray<T>::Min() const
{
   T minval = std::numeric_limits<T>::max();
   for( int n=this->m_nLowerLimit; n<=this->m_nUpperLimit; n++ )
   {
      if( this->m_pArray[n] < minval )
         minval = this->m_pArray[n];
   }
   return minval;
}



//-----
// Max
//-----

template <typename T>
T NumericArray<T>::Max() const
{
   T maxval = std::numeric_limits<T>::min();
   for( int n=this->m_nLowerLimit; n<=this->m_nUpperLimit; n++ )
   {
      if( this->m_pArray[n] > maxval )
         maxval = this->m_pArray[n];
   }
   return maxval;
}



template <typename T>
void NumericArray<T>::Interpolate( int x1, int x2 )
{
   // Linear interpolation between two points using equation y = mx + c
   assert(x1<x2);
   assert(x1>=0);
   assert(x2<this->m_nLength);
   if( x2 > x1 )
   {
      int    points = x2 - x1;
      T      c      = this->m_pArray[x1];
      double m      = static_cast<double>(this->m_pArray[x2] - this->m_pArray[x1]) / static_cast<double>(points);
      for( int x=0; x<points; x++, x1++ )
	  this->m_pArray[x1] = static_cast<T>( m*x + c );
  }
}




//---------------------------------- DNA Array --------------------------------




template <typename CharT>
class DNAArray : public SimpleArray<CharT>
{
 public:
    // Services
    bool IsACGT( int n ) const;
    int  CountPads( char cPad='*' ) const;
    int  GetAlignedPosition( int n, bool bLeft=true, char cPad='*' ) const;
    int  GetOriginalPosition( int i, bool bLeft=true, char cPad='*' ) const;

    // Now obsolete, should be removed eventually!
    int  GotoBase( int nBasePos, char pad='*' ) const;
};




//--------
// IsACGT
//--------

template <typename CharT>
bool DNAArray<CharT>::IsACGT( int nBasePos ) const
{
    assert(nBasePos>=0);
    assert(nBasePos<this->m_nLength);
    switch( this->m_pArray[nBasePos] )
    {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'a':
        case 'c':
        case 'g':
        case 't':
        return true;
    }
    return false;
}



//-----------
// Goto Base
//-----------

template <typename CharT>
int DNAArray<CharT>::GotoBase( int nBasePos, char cPad ) const
{
/*
   Returns the index of base number 'nBasePos' in this sequence obtained by
   counting bases while skipping over any pad characters as specified by 'cPad'.
   Base numbers are assumed to be zero-based.
*/
   assert(nBasePos>=0);
   assert(nBasePos<this->m_nLength);
   int bp = 0;
   for( int n=0; n<this->m_nLength; n++ )
   {
      if( this->m_pArray[n] == cPad )
         continue;
      if( bp == nBasePos )
         return n;
      bp++;
   }
   return -1;
}



//------------
// Count Pads
//------------

template <typename CharT>
int DNAArray<CharT>::CountPads( char cPad ) const
{
/*
    Counts the number of pad characters contained in the current range.
*/
    int nPads = 0;
    for( int n=this->m_nLowerLimit; n<=this->m_nUpperLimit; n++ )
    {
        if( this->m_pArray[n] == cPad )
            nPads++;
    }
    return nPads;
}




//------------------
// Position Mapping - Zero based
//------------------

template <typename CharT>
int DNAArray<CharT>::GetAlignedPosition( int n, bool bLeft, char cPad ) const
{
    int i;


    // If position is relative to the left
    if( bLeft )
    {
        // Scan from left
        i = 0;
        while( (i<this->m_nLength) && (n>0) )
        {
            if( this->m_pArray[i] != cPad )
                n--;
            i++;
        }
    }
    else
    {
        // Scan from right
        i = this->m_nLength - 1;
        while( (i>=0) && (n>0) )
        {
            if( this->m_pArray[i] != cPad )
                n--;
            i--;
        }
    }
    return (n>0) ? -1 : i;
}



template <typename CharT>
int DNAArray<CharT>::GetOriginalPosition( int i, bool bLeft, char cPad ) const
{
    assert(i>=0);
    assert(i<this->m_nLength);


    // Check input
    int n = -1;
    if( i >= this->m_nLength )
        return n;


    // If position is relative to the left
    if( bLeft )
    {
        // Scan to left
        while( i >= 0 )
        {
            if( this->m_pArray[i] != cPad )
                n++;
            i--;
        }
    }
    else
    {
        // Scan to right
        while( i < this->m_nLength )
        {
            if( this->m_pArray[i] != cPad )
                n++;
            i++;
        }
    }
    return n;
}



#endif
