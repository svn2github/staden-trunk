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



#ifndef _MUTLIB_TAGARRAY_HPP_
#define _MUTLIB_TAGARRAY_HPP_



#include <mutlib.h>
#include <list.hpp>
#include <muttag.hpp>



class TagArray
{
 public:
   // Constructor/Destructor
   TagArray()                             { Init();  }
  ~TagArray()                             { Empty(); }


 public:
   // Services
   mutlib_tag_t* Raw()                    { return m_pArray; }
   void          Empty();
   int           Length() const           { return m_nLength; }
   void          AutoDestroy( bool s )    { m_bAutoDestroy=s; }
   void          Create( int nLength );
   void          ReadTags( List<MutTag>& rList, int nPositionIndex, bool bComplement );
   void          Wrap( mutlib_tag_t* a, int nLength );



 public:
   // Operators
   mutlib_tag_t& operator[]( int n )      { assert(n<m_nLength); return m_pArray[n]; }


 private:
   // Services
   void Init();


 private:
   // Data
   mutlib_tag_t* m_pArray;
   int           m_nLength;
   bool          m_bAutoDestroy;
};



#endif

