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


#ifndef _MUTLIB_LISTITEM_HPP_
#define _MUTLIB_LISTITEM_HPP_



template <typename ItemT>
class ListItem
{
 public:
   // Constructor
   ListItem();



 public:
   // Services
   void   Copy( const ItemT& rhs );
   ItemT* Next() const              { return m_pNext; }
   ItemT* Prev() const              { return m_pPrev; }
   void   Next( ItemT* p )          { m_pNext=p; }
   void   Prev( ItemT* p )          { m_pPrev=p; }
   bool   Marked() const            { return m_bMarked; }
   void   Marked( bool f )          { m_bMarked=f; }



 private:
   // Data
   ItemT* m_pNext;
   ItemT* m_pPrev;
   bool   m_bMarked;
};



//-------------
// Constructor
//-------------

template<typename ItemT>
ListItem<ItemT>::ListItem()
{
   m_pNext   = 0;
   m_pPrev   = 0;
   m_bMarked = false;
}



//------
// Copy
//------

template <typename ItemT>
void ListItem<ItemT>::Copy( const ItemT& rhs )
{
   // Ensures proper copy, don't copy prev/next pointers
   m_pNext   = 0;
   m_pPrev   = 0;
   m_bMarked = rhs.Marked();
}



#endif
