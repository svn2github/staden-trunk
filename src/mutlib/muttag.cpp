/*
 * Copyright (c) Medical Research Council 2000. All rights reserved.
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


#include <locale>            // For toupper()
#include <cstring>           // For strcpy(), strlen(), ...
#include <cstdio>            // For sprintf()
#include <muttag.hpp>


//--------
// Tables
//--------

// Must have same layout as enum values in muttag.hpp
static const char* pCommentTable[] = {
"",
"A->A",
"A->C",
"A->G",
"A->T",
"A->*",
"A->?",
"C->A",
"C->C",
"C->G",
"C->T",
"C->*",
"C->?",
"G->A",
"G->C",
"G->G",
"G->T",
"G->*",
"G->?",
"T->A",
"T->C",
"T->G",
"T->T",
"T->*",
"T->?",
"*->A",
"*->C",
"*->G",
"*->T",
"*->*",
"*->?",
"?->A",
"?->C",
"?->G",
"?->T",
"?->*",
"?->?"
};



//-------------
// Constructor
//-------------

MutTag::MutTag( const char* Name, mutlib_mutation_t mt, int p, mutlib_strand_t d )
{
   assert(Name != NULL);
   assert(std::strlen(Name)<5);
   m_nType                   = mt;
   m_nStrand                 = d;
   m_nConfidence             = 0;
   std::strncpy( m_pName, Name, 4 );
   int n;
   for( n=0; n<4; n++ )
       m_pName[n] = std::toupper( m_pName[n] );
   m_pName[n]                = 0;
   m_pComment[0]             = 0;
   m_nPosition[0]            = p;
   m_nPosition[1]            = 0;
   m_nPosition[2]            = 0;
   m_nAmplitude[0]           = 0;
   m_nAmplitude[1]           = 0;
   m_nAmplitudeRatio         = 0.0;
   m_nAmplitudePercentage[0] = 0.0;
   m_nAmplitudePercentage[1] = 0.0;
   m_nWidth                  = 0.0;
   m_nAlignment              = 0.0;
   m_nSensitivity            = 0.0;
}




//-------
// Clone
//-------

void MutTag::Clone( const MutTag& rhs )
{
   // Ensures proper copy of list item pointers and contained strings
   Copy( rhs );
   m_nType                   = rhs.m_nType;
   m_nStrand                 = rhs.m_nStrand;
   m_nConfidence             = rhs.m_nConfidence;
   std::strcpy( m_pName, rhs.m_pName );
   std::strcpy( m_pComment, rhs.m_pComment );
   m_nPosition[0]            = rhs.m_nPosition[0];
   m_nPosition[1]            = rhs.m_nPosition[1];
   m_nPosition[2]            = rhs.m_nPosition[2];
   m_nAmplitude[0]           = rhs.m_nAmplitude[0];
   m_nAmplitude[1]           = rhs.m_nAmplitude[1];
   m_nAmplitudeRatio         = rhs.m_nAmplitudeRatio;
   m_nAmplitudePercentage[0] = rhs.m_nAmplitudePercentage[0];
   m_nAmplitudePercentage[1] = rhs.m_nAmplitudePercentage[1];
   m_nWidth                  = rhs.m_nWidth;
   m_nAlignment              = rhs.m_nAlignment;
   m_nSensitivity            = rhs.m_nSensitivity;
}



//------------
// Complement
//------------

void MutTag::Complement( char* s )
{
/*
   Complements any base calls in the string. ie.

   A -> T
   C -> G
   G -> C
   T -> A
*/

   if(!s) return;
   int slen = std::strlen(s);
   for( int n=0; n<slen; n++ )
   {
       switch( s[n] )
       {
          case 'A':
          case 'a': s[n] = 'T'; break;
          case 'C':
          case 'c': s[n] = 'G'; break;
          case 'G':
          case 'g': s[n] = 'C'; break;
          case 'T':
          case 't': s[n] = 'A'; break;
       }
    }
}



//---------
// Comment
//---------

const char* MutTag::Comment( bool bComplement )
{
    // Copy comment from table
    std::strcpy( m_pComment, pCommentTable[m_nType] );
    if( std::strcmp(m_pName,"HETE") == 0 )
    {
        // Remove arrow "->" for heterozygotes
        m_pComment[1] = m_pComment[3];
        m_pComment[2] = 0;
    }
    if( (m_nStrand==MUTLIB_STRAND_REVERSE) && bComplement )
        Complement( m_pComment );
    int slen = std::strlen( m_pComment );


    // Append mutation parameter measurements
    if( std::strcmp(m_pName,"MUTA") == 0 )
    {
        std::sprintf( &m_pComment[slen], " Sensitivity=%5.2f, Alignment=%4.2f, Width=%4.2f, Amplitude=%d",
                      m_nSensitivity, m_nAlignment, m_nWidth, 
                      static_cast<int>(m_nAmplitude[0]+m_nAmplitude[1]) );
    }


    // Append heterozygote parameter measurements
    else if( std::strcmp(m_pName,"HETE") == 0 )
    {
        std::sprintf( &m_pComment[slen], " Ratio=%4.2f, Alignment=%4.2f, Amplitude1=%4.2f, Amplitude2=%4.2f",
                      m_nAmplitudeRatio, m_nAlignment,
                      m_nAmplitudePercentage[0], m_nAmplitudePercentage[1] );
    }
    assert(std::strlen(m_pComment)<MAX_COMMENT);
    return m_pComment;
}



//------
// Type
//------

void MutTag::Type( int ni, int nw )
{
    m_nType = mutlib_mutation_t( MUTLIB_MUTATION_AA + (nw*6) + ni );
}



//-------------
// BaseToIndex
//-------------

int MutTag::BaseToIndex( int b )
{
    switch( b )
    {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        case '*':           return 4;
        default:            return 5;
    }
}

