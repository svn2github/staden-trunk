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


#include <cassert>
#include <cstring>             // For strcpy(), strlen()
#include <stringlist.hpp>


//------------------
// Node Constructor
//------------------

StringListNode::StringListNode( const char* s )
{
    assert(s != NULL);
    if( !s )
        s = "";
    m_pString = new char[ std::strlen(s)+1 ];
    std::strcpy( m_pString, s );
    m_pNext = 0;
}



//-------------
// Constructor
//-------------

StringList::StringList()
{
    m_pList    = 0;
    m_nLength  = 0;
    m_pCurrent = 0;
}



//------------
// Destructor
//------------

StringList::~StringList()
{
    StringListNode* pNext;
    StringListNode* pNode = m_pList;
    while( pNode )
    {
        pNext = pNode->Next();
        delete pNode;
        pNode = pNext;
    }
}



//-------
// First
//-------

char* StringList::First()
{
    if( m_pList )
    {
        m_pCurrent = m_pList;
        return m_pCurrent->Data();
    }
    return 0;
}



//------
// Next
//------

char* StringList::Next()
{
    if( m_pCurrent )
    {
        if( m_pCurrent->Next() )
        {
            m_pCurrent = m_pCurrent->Next();
            return m_pCurrent->Data();
        }
    }
    return 0;
}



//--------
// Append
//--------

void StringList::Append( const char* s )
{
    StringListNode* pNode = new StringListNode(s);
    if( !m_pList )
        m_pList = pNode;
    else
    {
        StringListNode* pNext = m_pCurrent->Next();
        while( pNext )
        {
            m_pCurrent = pNext;
            pNext = m_pCurrent->Next();
        }
        m_pCurrent->Next(pNode);
    }
    m_nLength++;
    m_pCurrent = pNode;
}



