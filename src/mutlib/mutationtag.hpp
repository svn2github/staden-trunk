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


#ifndef _MUTLIB_MUTATIONTAG_HPP_
#define _MUTLIB_MUTATIONTAG_HPP_


#include <cassert>
#include <mutlib.h>         // For mutlib_strand_t
#include <listitem.hpp>     // For ListItem template


/*
    This tag class stores the information required to mark a mutation.
*/
class MutationTag : public ListItem<MutationTag>
{
 public:
   // Constructors
   MutationTag( const MutationTag& rhs )         { Clone(rhs); }
   MutationTag( const char* name );



 public:
   // Services
   const char*     Name()                        { return m_pName; }
   void            Name( const char* newname );
   const char*     Comment();
   char            BaseInput( int n=0 ) const    { assert((n>=0)&&(n<3)); return m_cBaseInp[n]; }
   void            BaseInput( int n, char c )    { assert((n>=0)&&(n<3)); m_cBaseInp[n]=c; }
   char            BaseReference() const         { return m_cBaseRef; }
   void            BaseReference( char c )       { m_cBaseRef=c; }
   mutlib_strand_t Strand() const                { return m_nStrand; }
   void            Strand( mutlib_strand_t t )   { m_nStrand=t; }
   int             Position( int n=0 ) const     { assert(n<4); return m_nPosition[n]; }
   void            Position( int n, int v )      { assert(n<4); m_nPosition[n]=v; }
   double          SNR() const                   { return m_nSNR; }
   void            SNR( double n )               { m_nSNR=n; }
   double          Amplitude( int n ) const      { assert(n<3); return m_nAmplitude[n]; }
   void            Amplitude( int n, double a )  { assert(n<3); m_nAmplitude[n]=a; }
   int             Row() const                   { return m_nRow; }
   int             Col() const                   { return m_nCol; }
   void            Row( int n )                  { m_nRow=n; }
   void            Col( int n )                  { m_nCol=n; }
   void            Mark( bool s )                { m_bMarked=s; }
   bool            Marked() const                { return m_bMarked; }



 public:
    // Operators
    void operator=( const MutationTag& rhs )     { Clone(rhs); }



 private:
     // Helpers
     void Clone( const MutationTag& );



 private:
   // Data
   enum { MAX_STRING=80 };
   char            m_cBaseRef;                  // Reference sequence basecall
   char            m_cBaseInp[3];               // Input sequence basecalls
   mutlib_strand_t m_nStrand;                   // Forward or reverse strand
   double          m_nSNR;                      // Signal to noise ratio at site (dB)
   char            m_pName[8];                  // Tag name, HETE, MUTA
   char            m_pComment[MAX_STRING];      // Tag comment field
   int             m_nPosition[4];              // Peak/Base position values
   double          m_nAmplitude[3];             // Storage for peak amplitude info
   bool            m_bMarked;                   // Used to mark tags for processing
   int             m_nRow;                      // Row+Col peak matrix reference
   int             m_nCol;                      //
};



#endif
