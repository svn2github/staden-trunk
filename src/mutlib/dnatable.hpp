/*
 * Copyright (c) Medical Research Council 2002. All rights reserved.
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


#ifndef _MUTLIB_DNATABLE_HPP_
#define _MUTLIB_DNATABLE_HPP_


/**
   Lookup table for base letters, including ambiguity codes.
*/
class DNATable
{
 public:
   // Services
   int  LookupIndex( char c ) const;
   char LookupBase( int index ) const;
   char LookupBase( int index1, int index2 ) const;
   char LookupBase( char char1, char char2 ) const;
   bool IsBaseAmbiguous( char c ) const;


 private:
   // Data
   static const char m_IndexTable[4][4];
};



#endif

