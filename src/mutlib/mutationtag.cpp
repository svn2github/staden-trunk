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


#include <cassert>
#include <cstring>           // For strcpy(), strlen(), ...
#include <cstdio>            // For sprintf()
#include <mutationtag.hpp>



//-------------
// Constructor
//-------------

MutationTag::MutationTag( const char* name )
{
   assert(name != NULL);
   assert(std::strlen(name)<5);
   m_cBaseRef      = '-';
   m_cBaseInp[0]   = '-';
   m_cBaseInp[1]   = '-';
   m_cBaseInp[2]   = '-';
   m_nStrand       = MUTLIB_STRAND_FORWARD;
   m_nSNR          = 0.0;
   m_nPosition[0]  = 0;
   m_nPosition[1]  = 0;
   m_nPosition[2]  = 0;
   m_nPosition[3]  = 0;
   m_nAmplitude[0] = 0.0;
   m_nAmplitude[1] = 0.0;
   m_nAmplitude[2] = 0.0;
   m_pComment[0]   = 0;
   m_bMarked       = false;
   m_nRow          = 0;
   m_nCol          = 0;
   Name( name );
}



//-------
// Clone
//-------

void MutationTag::Clone( const MutationTag& rhs )
{
   // Ensures proper copy of list item pointers and contained strings
   Copy( rhs );
   m_cBaseRef      = rhs.m_cBaseRef;
   m_cBaseInp[0]   = rhs.m_cBaseInp[0];
   m_cBaseInp[1]   = rhs.m_cBaseInp[1];
   m_cBaseInp[2]   = rhs.m_cBaseInp[2];
   m_nStrand       = rhs.m_nStrand;
   m_nSNR          = rhs.m_nSNR;
   m_nPosition[0]  = rhs.m_nPosition[0];
   m_nPosition[1]  = rhs.m_nPosition[1];
   m_nPosition[2]  = rhs.m_nPosition[2];
   m_nPosition[3]  = rhs.m_nPosition[3];
   m_nAmplitude[0] = rhs.m_nAmplitude[0];
   m_nAmplitude[1] = rhs.m_nAmplitude[1];
   m_nAmplitude[2] = rhs.m_nAmplitude[2];
   m_bMarked       = rhs.m_bMarked;
   m_nRow          = rhs.m_nRow;
   m_nCol          = rhs.m_nCol;
   std::strcpy( m_pName, rhs.m_pName );
   std::strcpy( m_pComment, rhs.m_pComment );
}



//---------
// Comment
//---------

const char* MutationTag::Comment()
{
    if( std::strcmp(m_pName,"HETE") == 0 )
    {
        std::sprintf( m_pComment, "%c->%c, SNR=%0.2fdB, PKD=%0.2f",
                      m_cBaseRef, m_cBaseInp[0], m_nSNR, m_nAmplitude[2] );
    }


    if( std::strcmp(m_pName,"MUTA") == 0 )
    {
        std::sprintf( m_pComment, "%c->%c, SNR=%0.2fdB",
                      m_cBaseRef, m_cBaseInp[0], m_nSNR );
    }


    // Check for overflow!
    assert(std::strlen(m_pComment)<MAX_STRING);
    return m_pComment;
}



//------
// Name
//------

void MutationTag::Name( const char* newname )
{
    assert(newname != NULL);
    assert(std::strlen(newname)==4);
    std::strcpy( m_pName, newname );
}
