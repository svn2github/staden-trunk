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
#include <cstring>           // For memset(), strcpy(), strlen()...
#include <tagarray.hpp>


//------
// Init
//------

void TagArray::Init()
{
   m_pArray       = 0;
   m_nLength      = 0;
   m_bAutoDestroy = true;
}



//-------
// Empty
//-------

void TagArray::Empty()
{
   if( m_pArray && m_bAutoDestroy )
   {
      // Delete all comments, then array
      for( int n=0; n<m_nLength; n++ )
         delete [] m_pArray[n].Comment;
      delete [] m_pArray;
   }
   Init();
}



//--------
// Create
//--------

void TagArray::Create( int nLength )
{
   assert(nLength>0);
   if( m_pArray )
      Empty();
   m_pArray  = new mutlib_tag_t[ nLength ];
   m_nLength = nLength;
   std::memset( m_pArray, 0, nLength*sizeof(mutlib_tag_t) );
}



//------
// Wrap
//------

void TagArray::Wrap( mutlib_tag_t* a, int nLength )
{
   assert(nLength>0);
   if( m_pArray )
      Empty();
   m_pArray  = a;
   m_nLength = nLength;
}



//-----------
// Read Tags
//-----------

void TagArray::ReadTags( List<MutTag>& rList, int nPositionIndex, bool bComplement )
{
   const char* s;
   MutTag*     pTag;


   // Copies the content of each tag into the tag array
   pTag = rList.First();
   for( int n=0; pTag && (n<m_nLength); n++ )
   {
      std::strcpy( m_pArray[n].Type, pTag->Name() );
      assert(std::strlen(m_pArray[n].Type)<=4);
      m_pArray[n].Strand   = pTag->Strand();
      m_pArray[n].Position[0] = pTag->Position( nPositionIndex );
      m_pArray[n].Position[1] = 0;
      s = pTag->Comment( bComplement );
      m_pArray[n].Comment = new char[ std::strlen(s)+1 ];
      std::strcpy( m_pArray[n].Comment, s );
      pTag = rList.Next();
   }
}
