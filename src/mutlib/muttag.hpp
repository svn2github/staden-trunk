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


#ifndef _MUTLIB_MUTTAG_HPP_
#define _MUTLIB_MUTTAG_HPP_


#include <cassert>
#include <mutlib.h>         // For mutlib_strand_t
#include <listitem.hpp>     // For ListItem template



//-------------
// Definitions
//-------------

typedef enum
{
   // Define some mutation types. MUTLIB_MUTATION_CT means a C to T mutation.
   // P=Padding character, N=unknown/uncertain
   MUTLIB_MUTATION_NONE,
   MUTLIB_MUTATION_AA,
   MUTLIB_MUTATION_AC,
   MUTLIB_MUTATION_AG,
   MUTLIB_MUTATION_AT,
   MUTLIB_MUTATION_AP,
   MUTLIB_MUTATION_AN,
   MUTLIB_MUTATION_CA,
   MUTLIB_MUTATION_CC,
   MUTLIB_MUTATION_CG,
   MUTLIB_MUTATION_CT,
   MUTLIB_MUTATION_CP,
   MUTLIB_MUTATION_CN,
   MUTLIB_MUTATION_GA,
   MUTLIB_MUTATION_GC,
   MUTLIB_MUTATION_GG,
   MUTLIB_MUTATION_GT,
   MUTLIB_MUTATION_GP,
   MUTLIB_MUTATION_GN,
   MUTLIB_MUTATION_TA,
   MUTLIB_MUTATION_TC,
   MUTLIB_MUTATION_TG,
   MUTLIB_MUTATION_TT,
   MUTLIB_MUTATION_TP,
   MUTLIB_MUTATION_TN,
   MUTLIB_MUTATION_PA,
   MUTLIB_MUTATION_PC,
   MUTLIB_MUTATION_PG,
   MUTLIB_MUTATION_PT,
   MUTLIB_MUTATION_PP,
   MUTLIB_MUTATION_PN,
   MUTLIB_MUTATION_NA,
   MUTLIB_MUTATION_NC,
   MUTLIB_MUTATION_NG,
   MUTLIB_MUTATION_NT,
   MUTLIB_MUTATION_NP,
   MUTLIB_MUTATION_NN

}mutlib_mutation_t;



//--------------
// Mutation Tag
//--------------


class MutTag : public ListItem<MutTag>
{
 public:
   // Constructors
   MutTag( const MutTag& rhs )                             { Clone(rhs); }
   MutTag( const char* Name, mutlib_mutation_t mt, int p, mutlib_strand_t d );



 public:
   // Services
   mutlib_mutation_t Type() const                           { return m_nType; }
   void              Type( int ni, int nw );
   int               BaseToIndex( int b );
   const char*       Name() const                           { return m_pName; }
   const char*       Comment( bool bComplement=false );
   mutlib_strand_t   Strand() const                         { return m_nStrand; }
   void              Strand( mutlib_strand_t t )            { m_nStrand=t; }
   int               Position( int n=0 ) const              { assert(n<3); return m_nPosition[n]; }
   void              Position( int n, int v )               { assert(n<3); m_nPosition[n]=v; }
   int               Amplitude( int n )                     { assert(n<2); return m_nAmplitude[n]; }
   void              Amplitude( int n, int a )              { assert(n<2); m_nAmplitude[n]=a; }
   int               Confidence() const                     { return m_nConfidence; }
   void              Confidence( int n )                    { m_nConfidence=n; }
   double            Sensitivity() const                    { return m_nSensitivity; }
   void              Sensitivity( double s )                { m_nSensitivity=s; }
   double            Width() const                          { return m_nWidth; }
   void              Width( double w )                      { m_nWidth=w; }
   double            Alignment() const                      { return m_nAlignment; }
   void              Alignment( double a )                  { m_nAlignment=a; }
   double            AmplitudeRatio() const                 { return m_nAmplitudeRatio; }
   void              AmplitudeRatio( double v )             { m_nAmplitudeRatio=v; }
   void              AmplitudePercentage( int n, double v ) { assert(n<2); m_nAmplitudePercentage[n]=v; }
   double            AmplitudePercentage( int n ) const     { assert(n<2); return m_nAmplitudePercentage[n]; }



 public:
    // Operators
    void operator=( const MutTag& rhs )                      { Clone(rhs); }



 private:
     // Helpers
     void Complement( char* s );
     void Clone( const MutTag& );


 private:
   // Data
   enum { MAX_COMMENT=80 };
   mutlib_mutation_t m_nType;                     // Type of mutation
   mutlib_strand_t   m_nStrand;                   // Forward or reverse strand
   int               m_nConfidence;               // Level of confidence it's a mutation 0-100
   char              m_pName[8];                  // Tag name
   char              m_pComment[MAX_COMMENT];     // Tag comment field
   int               m_nPosition[3];              // Peak position values
   int               m_nAmplitude[2];             // Peak amplitude values
   double            m_nAmplitudeRatio;	          // Double peak amplitude ratio
   double            m_nAmplitudePercentage[2];   // Amplitude of peaks as percentage of max
   double            m_nWidth;                    // Largest peak width in bases
   double            m_nAlignment;                // Alignment of peaks in samples
   double            m_nSensitivity;
};



#endif
