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


#ifndef _MUTLIB_MATRIX_HPP_
#define _MUTLIB_MATRIX_HPP_



#include <cassert>
#include <fstream>           // For std::ofstream
#include <iomanip>           // For std::setw()



template <typename T>
class SimpleMatrix
{
 public:
   // Constructor/Destructor
   SimpleMatrix()                       { Init(); }
  ~SimpleMatrix()                       { Empty(); }



 public:
   // Services
   T**  Raw()                          { return m_pMatrix; }
   void Empty();
   void Fill( T v );
   int  Rows() const                   { return m_nRows; }
   void Rows( int n )                  { m_nRows=n; }
   int  Cols() const                   { return m_nCols; }
   void Cols( int n )                  { m_nCols=n; }
   int  RowCapacity() const            { return m_nRowCapacity; }
   int  ColCapacity() const            { return m_nColCapacity; }
   void Create( int nRows, int nCols );
   void Wrap( T** p, int nRows, int nCols, bool AutoDestroy=false );
   void SaveAs( const char* fname, int field_width=6 );



 public:
   // Operators
   T*& operator[]( int n )             { assert(n<m_nRowCapacity); return m_pMatrix[n]; }



 protected:
   T**  m_pMatrix;
   int  m_nRows;
   int  m_nCols;
   int  m_nRowCapacity;
   int  m_nColCapacity;
   bool m_bAutoDestroy;



 private:
   // Helpers
   void Init();
};



//------
// Init
//------

template <typename T>
void SimpleMatrix<T>::Init()
{
   m_pMatrix      = 0;
   m_nRows        = 0;
   m_nCols        = 0;
   m_nRowCapacity = 0;
   m_nColCapacity = 0;
   m_bAutoDestroy = true;
}



//--------
// Create
//--------

template <typename T>
void SimpleMatrix<T>::Create( int nRows, int nCols )
{ 
   int r;
   assert(nRows>0);
   assert(nCols>0);


   // Discard old matrix
   if( m_pMatrix )
      Empty();



   // Allocate base storage
   m_pMatrix = new T*[ nRows ];
   for( r=0; r<nRows; r++ )
      m_pMatrix[r] = 0;
   m_nRows        = nRows;
   m_nRowCapacity = nRows;



   // Allocate row vectors
   for( r=0; r<nRows; r++ )
      m_pMatrix[r] = new T[ nCols ];
   m_nCols        = nCols;
   m_nColCapacity = nCols;
   m_bAutoDestroy = true;
}



//------
// Wrap
//------

template <typename T>
void SimpleMatrix<T>::Wrap( T** p, int nRows, int nCols, bool AutoDestroy )
{
   assert(p != NULL);
   assert(nRows>0);
   assert(nCols>0);
   if( m_pMatrix )
      Empty();
   m_pMatrix      = p;
   m_nRows        = nRows;
   m_nCols        = nCols;
   m_nRowCapacity = nRows;
   m_nColCapacity = nCols;
   m_bAutoDestroy = AutoDestroy;
}




//-------
// Empty
//-------

template <typename T>
void SimpleMatrix<T>::Empty()
{ 
   if( m_bAutoDestroy == true )
   {
      for( int r=0; r<m_nRows; r++ )
         delete [] m_pMatrix[r];
      delete [] m_pMatrix;
   }
   Init();
}



//------
// Fill
//------

template <typename T>
void SimpleMatrix<T>::Fill( T v )
{
   int r, c;

   // Fill matrix
   for( r=0; r<m_nRows; r++ )
      for( c=0; c<m_nCols; c++ )
         m_pMatrix[r][c] = v;
}



//---------
// Save As
//---------

template <typename T>
void SimpleMatrix<T>::SaveAs( const char* fname, int field_width )
{
#ifdef _WIN32
    using namespace std;
    ofstream ofs;
    ofs.open( fname, ios_base::out | ios_base::trunc );
    for( int r=0; r<m_nRows; r++ )
    {
        for( int c=0; c<m_nCols; c++ )
            ofs << setw(field_width) << m_pMatrix[r][c];
        ofs << endl;
    }
    ofs << endl;
#endif
}


#endif
