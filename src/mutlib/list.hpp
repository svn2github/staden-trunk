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



#ifndef _MUTLIB_LIST_HPP_
#define _MUTLIB_LIST_HPP_


#include <cassert>



template <typename T>
class List
{
 public:
   // Constructor/Destructor
   List();
  ~List();



 public:
   // Services
   void Empty();
   int  Count() const               { return m_nCount; }
   int  Index() const               { return m_nIndex; }
   T*   Next();
   T*   Prev();
   T*   First();
   T*   Current()                   { return m_pCurrent; }
   void Append( T* p );
   T*   Goto( int n );
   void Pop( int n=0 );
   void Push( int n=0 );
   T*   Remove( int n );
   void Insert( T* p, int n );



 private:
   // Data
   T*  m_pList;
   int m_nIndex;
   int m_nCount;
   T*  m_pCurrent;
   //---
   int m_nStackedIndex[4];
};



//-------------
// Constructor
//-------------

template <typename T>
List<T>::List()
{
   m_pList    = 0;
   m_nIndex   = 0;
   m_nCount   = 0;
   m_pCurrent = 0;
}



//-------
// Empty
//-------

template <typename T>
void List<T>::Empty()
{
   // Delete all tags
   T* pTmp;
   T* p = First();
   while( p )
   {
      pTmp = p->Next();
      delete p;
      p = pTmp;
   }
   m_pList    = 0;
   m_nCount   = 0;
   m_nIndex   = 0;
   m_pCurrent = 0;
}




//------------
// Destructor
//------------

template <typename T>
List<T>::~List()
{
   Empty();
}



//------
// Push
//------

template <typename T>
void List<T>::Push( int n )
{
    assert(n>=0);
    assert(n<4);
    m_nStackedIndex[n] = m_nIndex;
}




//-----
// Pop
//-----

template <typename T>
void List<T>::Pop( int n )
{
    assert(n>=0);
    assert(n<4);
    Goto( m_nStackedIndex[n] );
}




//-------
// First
//-------

template <typename T>
T* List<T>::First()
{
   m_nIndex   = 0;
   m_pCurrent = m_pList;
   return m_pList;
}



//------
// Next
//------

template <typename T>
T* List<T>::Next()
{
   if( m_nCount && m_pCurrent->Next() )
   {
      m_nIndex++;
      m_pCurrent = m_pCurrent->Next();
      return m_pCurrent;
   }
   return 0;
}



//------
// Prev
//------

template <typename T>
T* List<T>::Prev()
{
   if( m_nCount && m_pCurrent->Prev() )
   {
      m_nIndex--;
      m_pCurrent = m_pCurrent->Prev();
      return m_pCurrent;
   }
   return 0;
}



//--------
// Append
//--------

template <typename T>
void List<T>::Append( T* p )
{
   assert(p != NULL);
   if( m_nCount )
   {
       while( Next() );
       m_pCurrent->Next( p );
       p->Next( 0 );
       p->Prev( m_pCurrent );
       m_nIndex++;
   }
   else
   {
       m_pList = p;
   }
   m_nCount++;
   m_pCurrent = p;
}



//------
// Goto
//------

template <typename T>
T* List<T>::Goto( int n )
{
   assert(n>=0);


   // Check argument
   if( (n<0) || (n>=m_nCount) )
       return 0;


   // Traverse list
   int nDistance = m_nIndex - n;
   if( nDistance < 0 )
   {
      // Traverse forwards
      while( nDistance++ )
         Next();
   }
   else if( nDistance > 0 )
   {
      // Traverse backwards
      while( nDistance-- )
         Prev();
   }
   return m_pCurrent;
}



//--------
// Remove
//--------

template <typename T>
T* List<T>::Remove( int n )
{
   assert(n>=0);
   assert(n<m_nCount);


   T* pTag = 0;
   if( (n>=0) && m_nCount && (n<m_nCount) )
   {
      pTag = Goto( n );
      if( pTag->Prev() == 0 )
      {
         // First
         m_pList    = pTag->Next();
         m_pCurrent = pTag->Next();
         if( m_pCurrent )
            m_pCurrent->Prev( 0 );
      }
      else if( pTag->Next() == 0 )
      {
         // Last
         m_pCurrent = pTag->Prev();
         m_pCurrent->Next( 0 );
         m_nIndex--;
      }
      else
      {
         // Middle
         m_pCurrent = pTag->Next();
         m_pCurrent->Prev( pTag->Prev() );
         m_pCurrent->Prev()->Next( m_pCurrent );
      }
      m_nCount--;
      pTag->Next(0);
      pTag->Prev(0);
   }
   return pTag;
}



//--------
// Insert
//--------

template <typename T>
void List<T>::Insert( T* p, int n )
{
   assert(p != NULL);
   assert(n>=0);


   // Check input
   if( !p || (n<0) )
      return;


   // Initialise item
   p->Next(0);
   p->Prev(0);


   // Do insertion
   if( n >= m_nCount )
      Append(p);
   else
   {
      T* pTag = Goto( n );
      if( pTag->Prev() == 0 )
      {
         // First
         m_pList    = p;
         m_pCurrent = p;
         p->Next( pTag );
         pTag->Prev( p );
      }
      else
      {
         // Middle
         m_pCurrent = p;
         p->Next( pTag );
         p->Prev( pTag->Prev() );
         pTag->Prev()->Next( p );
         pTag->Prev( p );
      }
      m_nCount++;
   }
}



#endif

