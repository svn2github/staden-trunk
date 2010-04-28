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

#include <staden_config.h>

#include <new>                // For std::bad_alloc()
#include <cstdio>             // For fprintf(), fputc()
#include <cctype>             // For tolower()
#include <cstring>            // For strlen()
#include <align.hpp>
#include <locale>


//------------
// Local Data
//------------

bool Alignment::m_bDNALookupInitialised = false;
static const char* CharSet              = "ACGTURYMWSKDHVBN-*";
static const int   ScoreMatrix[18][18]  = {
{ 4,  -8,  -8,  -8,  -8,   2,  -3,   2,   2,  -3,  -3,   1,   1,   1,  -4,   1,   1,   0 },
{-8,   4,  -8,  -8,  -8,  -3,   2,   2,  -3,   2,  -3,  -4,   1,   1,   1,   1,   1,   0 },
{-8,  -8,   4,  -8,  -8,   2,  -3,  -3,  -3,   2,   2,   1,  -4,   1,   1,   1,   1,   0 },
{-8,  -8,  -8,   4,   4,  -3,   2,  -3,   2,  -3,   2,   1,   1,  -4,   1,   1,   1,   0 },
{-8,  -8,  -8,   4,   4,  -3,   2,  -3,   2,  -3,   2,   1,   1,  -4,   1,   1,   1,   0 },
{ 2,  -3,   2,  -3,  -3,   2,  -4,   0,   1,   1,   1,   1,   0,   1,   0,   0,   0,   0 },
{-3,   2,  -3,   2,   2,  -4,   2,   0,   1,   1,   1,   0,   1,   0,   1,   0,   0,   0 },
{ 2,   2,  -3,  -3,  -3,   0,   0,   2,   1,   1,  -3,   0,   1,   1,   0,   0,   0,   0 },
{ 2,  -3,  -3,   2,   2,   1,   1,   1,   2,  -3,   1,   1,   1,   0,   0,   0,   0,   0 },
{-3,   2,   2,  -3,  -3,   1,   1,   1,  -3,   2,   1,   0,   0,   1,   1,   0,   0,   0 },
{-3,  -3,   2,   2,   2,   1,   1,  -3,   1,   1,   2,   1,   0,   0,   1,   0,   0,   0 },
{ 1,  -4,   1,   1,   1,   1,   0,   0,   1,   0,   1,   1,   0,   0,   0,   0,   0,   0 },
{ 1,   1,  -4,   1,   1,   0,   1,   1,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0 },
{ 1,   1,   1,  -4,  -4,   1,   0,   1,   0,   1,   0,   0,   0,   1,   0,   0,   0,   0 },
{-4,   1,   1,   1,   1,   0,   1,   0,   0,   1,   1,   0,   0,   0,   1,   0,   0,   0 },
{ 1,   1,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1 },
{ 1,   1,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1 },
{ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   4 }};



//-------------
// Constructor
//-------------

Alignment::Alignment()
{
   // NOTE. The edge score left/right parameters are important for
   //       the envelope alignments to be performed correctly at the
   //       chunk ends. Do not alter them!
   m_nBand              = 0;
   m_pParams            = 0;
   m_pOverlap           = 0;
   m_nPadSymbol         = '*';
   m_nEdgeScore         = EDGE_SCORE_RIGHT;
   m_nGapPenaltyBegin   = 12;
   m_nGapPenaltyExtend  = 4;
   for( int n=0; n<MAX_INPUT_SEQUENCES; n++ )
   {
      m_pInputSequence[n]       = 0;
      m_nInputSequenceLength[n] = 0;
   }
}



//------------
// Destructor
//------------

Alignment::~Alignment()
{
   // Cleanup
   if( m_pOverlap )
   {
      sp::destroy_overlap(m_pOverlap);
      m_pOverlap = 0;
   }
   if( m_pParams )
   {
      sp::destroy_align_params(m_pParams);
      m_pParams = 0;
   }
}



//------------
// Edge Score
//------------

void Alignment::EdgeScore( edge_score_t es )
{
   m_nEdgeScore = es;
}



//----------------
// Input Sequence
//----------------

void Alignment::InputSequence( int n, const char* s, int l )
{
   assert(n>=0);
   assert(n<MAX_INPUT_SEQUENCES);
   assert(s != NULL);
   assert(*s);
   m_pInputSequence[n]       = s;
   m_nInputSequenceLength[n] = (l<0) ? std::strlen(s) : l;
}



//--------
// Matrix
//--------

void Alignment::Matrix( int** m, int n, bool AutoDestroy )
{
   assert(m != NULL);
   assert(m[0] != NULL);
   assert(n>0);
   m_oMatrix.Wrap( m, n, n, AutoDestroy );
}



//---------
// Execute
//---------

int Alignment::Execute( algorithm_t a )
{
   int n;



   // Check there are at least 2 input sequences
   for( n=0; n<2; n++ )
   {
      if( !m_pInputSequence[n] || !(*(m_pInputSequence[n])) )
         return -1;
   }



   // Initialisation
   if( m_bDNALookupInitialised == false )
   {
      sp::init_DNA_lookup();
      m_bDNALookupInitialised = true;
   }



   // Create a weight matrix if none has been supplied
   if( m_oMatrix.Rows() <= 0 )
      CreateDefaultMatrix();



   // Do memory allocation
   if( !m_pParams )
   {
      m_pParams = sp::create_align_params();
      if( !m_pParams )
         throw std::bad_alloc();
   }
   if( m_pOverlap )
   {
      sp::destroy_overlap(m_pOverlap);
      m_pOverlap = 0;
   }
   m_pOverlap = sp::create_overlap();
   if( !m_pOverlap )
      throw std::bad_alloc();



   // Adjust alignment parameters
   sp::set_align_params( m_pParams, m_nBand, m_nGapPenaltyBegin, m_nGapPenaltyExtend,
                        SP_ALIGNMENT_RETURN_SEQ, 0, 0, m_nPadSymbol, m_nPadSymbol,
                        0, 0, a, 8, 0, m_nEdgeScore, 0.0, m_oMatrix.Raw() );



   // Adjust overlap parameters
   sp::init_overlap( m_pOverlap, (char*)m_pInputSequence[0], (char*)m_pInputSequence[1],
                    m_nInputSequenceLength[0], m_nInputSequenceLength[1] );



   // Do alignment
   return sp::aligner( m_pParams, m_pOverlap );
}



//-----------------
// Output Sequence
//-----------------

char* Alignment::OutputSequence( int n ) const
{
   // This horrible hack is required due to poor Overlap structure design.
   assert(n>=0);
   assert(n<MAX_INPUT_SEQUENCES);
   assert(m_pOverlap != NULL);
   switch(n)
   {
      case 0: return m_pOverlap->seq1_out;
      case 1: return m_pOverlap->seq2_out;
   }
   assert(0);
   return 0;
}



//------------------------
// Output Sequence Length
//------------------------

int Alignment::OutputSequenceLength( int n ) const
{
    assert(n>=0);
    assert(n<MAX_INPUT_SEQUENCES);
    assert(m_pOverlap != NULL);
    return m_pOverlap->seq_out_len;
}



//------------------------------
// Output Sequence Left Overlap
//------------------------------

int Alignment::OutputSequenceLeftOverlap( int n ) const
{
    assert(n>=0);
    assert(n<MAX_INPUT_SEQUENCES);
    assert(m_pOverlap != NULL);
    return m_pOverlap->left;
}



//-------------------------------
// Output Sequence Right Overlap
//-------------------------------

int Alignment::OutputSequenceRightOverlap( int n ) const
{
    assert(n>=0);
    assert(n<MAX_INPUT_SEQUENCES);
    assert(m_pOverlap != NULL);
    return m_pOverlap->right;
}



//--------------
// Output Score
//--------------

double Alignment::OutputScore() const
{
    assert(m_pOverlap != NULL);
    if( m_pOverlap->seq_out_len > 0 )
        return m_pOverlap->score / m_pOverlap->seq_out_len;
    else
        return m_pOverlap->score;
}



//-----------
// Debugging
//-----------

void Alignment::DumpToFile( const char* s, bool AsInteger ) const
{
    int         n;
    int         k;
    int         nLen;
    const char* pSeq;



    // Create file to dump the output to
    std::FILE* pOut = std::fopen( s, "wb" );
    if( !pOut )
        return;


    // Dump the input sequences
    for( k=0; k<2; k++ )
    {
        pSeq = m_pInputSequence[k];
        nLen = m_nInputSequenceLength[k];
        if( AsInteger )
        {
            for( n=0; n<nLen; n++ )
                std::fprintf( pOut, "%3d ", int(pSeq[n]) );
        }
        else
        {
            for( n=0; n<nLen; n++ )
                std::fputc( pSeq[n], pOut );
        }
        std::fprintf( pOut, "\r\n" );
    }


    // Dump the output sequences
    for( k=0; k<2; k++ )
    {
        pSeq = OutputSequence(k);
        nLen = OutputSequenceLength(k);
        if( AsInteger )
        {
            for( n=0; n<nLen; n++ )
                std::fprintf( pOut, "%3d ", int(pSeq[n]) );
        }
        else
        {
            for( n=0; n<nLen; n++ )
                std::fputc( pSeq[n], pOut );
        }
        std::fprintf( pOut, "\r\n" );
    }
    std::fclose( pOut );
}




//---------
// Helpers
//---------

void Alignment::CreateDefaultMatrix()
{
   int  n;
   int  nMax;
   char cr, cc;
   int  r,  c;



   // Compute optimal matrix size for given character set
   n    = 0;
   nMax = 0;
   while( CharSet[n] )
   {
      c = std::tolower( CharSet[n] );
      if( c > nMax )
         nMax = c;
      n++;
   }
   nMax++;



   // Create weight matrix, background score is slightly negative.
   m_oMatrix.Create( nMax, nMax );
   m_oMatrix.Fill( -1 );



   // Initialise weight matrix from score matrix
   int e = std::strlen(CharSet);
   for( r=0; r<e; r++ )
   {
      cr = CharSet[r];
      for( c=0; c<e; c++ )
      {
          cc = CharSet[c];
          m_oMatrix[             cr ][             cc ] = ScoreMatrix[r][c];
          m_oMatrix[std::tolower(cr)][             cc ] = ScoreMatrix[r][c];
          m_oMatrix[             cr ][std::tolower(cc)] = ScoreMatrix[r][c];
          m_oMatrix[std::tolower(cr)][std::tolower(cc)] = ScoreMatrix[r][c];
      }
   }
}

