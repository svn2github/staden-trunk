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



#ifndef _HETSCAN_STRING_LIST_HPP_
#define _HETSCAN_STRING_LIST_HPP_



class StringListNode
{
  public:
    // Constructor    
    StringListNode( const char* s );
   ~StringListNode()                            { delete [] m_pString; }



  public:
    // Services
    char*           Data()                      { return m_pString; }
    StringListNode* Next()                      { return m_pNext; }
    void            Next( StringListNode* p )   { m_pNext=p; }



  private:
    // Data
    char*           m_pString;
    StringListNode* m_pNext;
};



class StringList
{
  public:
    // Constructor/Destructor
    StringList();
   ~StringList();


  public:
    // Services
    int   Length() const                    { return m_nLength; }
    char* First();
    char* Next();
    void  Append( const char* s );


  private:
    // Data
    StringListNode* m_pList;
    int             m_nLength;
    StringListNode* m_pCurrent;
};



#endif


