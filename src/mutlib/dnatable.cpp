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


#include <cassert>
#include <dnatable.hpp>



const char DNATable::m_IndexTable[4][4] = {
{ 'A', 'M', 'R', 'W'},
{ 'M', 'C', 'S', 'Y'},
{ 'R', 'S', 'G', 'K'},
{ 'W', 'Y', 'K', 'T'}};



int DNATable::LookupIndex( char c ) const
{
    // Convert characters to indices
    int index;
    switch( c )
    {
        case 'A':
        case 'a': index = 0;  break;
        case 'C':
        case 'c': index = 1;  break;
        case 'G':
        case 'g': index = 2;  break;
        case 'T':
        case 't': index = 3;  break;
        case 'K':
        case 'k': index = 4;  break;
        case 'M':
        case 'm': index = 5;  break;
        case 'R':
        case 'r': index = 6;  break;
        case 'S':
        case 's': index = 7;  break;
        case 'W':
        case 'w': index = 8;  break;
        case 'Y':
        case 'y': index = 9;  break;
        default:  index =-1;  break;
    }
    return index;
}


char DNATable::LookupBase( int index ) const
{
    if( (index<0) || (index>3) )
        return '-';
    return m_IndexTable[index][index];
}


char DNATable::LookupBase( int index1, int index2 ) const
{
    if( (index1<0) || (index1>3) )
        return '-';
    if( (index2<0) || (index2>3) )
        return '-';
    return m_IndexTable[index1][index2];
}


char DNATable::LookupBase( char char1, char char2 ) const
{
    // Convert characters to indices
    int index[2] = { char1, char2 };
    for( int n=0; n<2; n++ )
    {
        switch( index[n] )
        {
            case 'A':
            case 'a': index[n] = 0; break;
            case 'C':
            case 'c': index[n] = 1; break;
            case 'G':
            case 'g': index[n] = 2; break;
            case 'T':
            case 't': index[n] = 3; break;
            default:  index[n] =-1;
        }
    }
    return LookupBase( index[0], index[1] );
}


bool DNATable::IsBaseAmbiguous( char c ) const
{
    switch( c )
    {
        case 'K':
        case 'k':
        case 'M':
        case 'm':
        case 'R':
        case 'r':
        case 'S':
        case 's':
        case 'W':
        case 'w':
        case 'Y':
        case 'y':
        return true;

    }
    return false;
}
