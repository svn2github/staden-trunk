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



#include <cstring>                      // For std::strcpy(), std::strlen()
#include <cstdlib>                      // For std::qsort()
#include <mutationtag_utils.hpp>


/*
   Complements any base calls in the tag. ie.

   A -> T
   C -> G
   G -> C
   T -> A
*/
void CompTags( SimpleArray<mutlib_tag_t>& a )
{
    for( int n=0; n<a.Length(); n++ )
    {
        // Find the -> arrow
        char* pos = std::strstr( a[n].Comment, "->" );
        if( pos )
        {
            pos--;
            int p = 2;
            while( p )
            {
                switch( *pos )
                {
                    // Complement the bases
                    case 'A': *pos = 'T'; break;
                    case 'C': *pos = 'G'; break;
                    case 'G': *pos = 'C'; break;
                    case 'T': *pos = 'A'; break;
                    case 'K': *pos = 'M'; break;
                    case 'M': *pos = 'K'; break;
                    case 'R': *pos = 'Y'; break;
                    case 'Y': *pos = 'R'; break;
                }
                pos++;
                pos++;
                pos++;
                p--;
            }
        }
    }
}




/**
    Sorts the tags in the array into ascending order of base position.
*/
extern "C" {
static int TagComparator( const void *t1, const void *t2 )
{
   const mutlib_tag_t* pt1 = static_cast<const mutlib_tag_t*>( t1 );
   const mutlib_tag_t* pt2 = static_cast<const mutlib_tag_t*>( t2 );
   return pt1->Position[0] - pt2->Position[0];
}
}


void SortTags( SimpleArray<mutlib_tag_t>& a )
{
    std::qsort( a.Raw(), a.Length(), sizeof(mutlib_tag_t), TagComparator );
}



/**
    Copies the tags from the taglist into an array ready for output.
*/
void CopyTags( SimpleArray<mutlib_tag_t>& a, List<MutationTag>& l )
{
    assert(a.Length()==l.Count());
    int          slen;
    int          n    = 0;
    MutationTag* pTag = l.First();
    while( pTag )
    {
        std::strcpy( a[n].Type, pTag->Name() );
        a[n].Strand      = pTag->Strand();
        a[n].Position[0] = pTag->Position();
        a[n].Position[1] = (std::strcmp(pTag->Name(),"MCOV")==0) ? pTag->Position(1) : pTag->Position();
        a[n].Marked      = pTag->Marked();
        slen = std::strlen( pTag->Comment() );
        a[n].Comment     = new char[ slen+1 ];
        a[n].Comment[0]  = 0;
        if( slen > 0 )
            std::strcpy( a[n].Comment, pTag->Comment() );
        pTag = l.Next();
        n++;
    }
}



/**
    Prunes the tag array so that there is only one tag at each position
    which has not been marked for deletion. If there's a base-change MUTA
    tag and a heterozygote HETE tag on the same base, the HETE tag is
    given priority. Assumes the tags are already sorted in order of position.
*/
void PruneTags( SimpleArray<mutlib_tag_t>& a )
{
    // Mark duplicate tags for deletion
    const int len = a.Length();
    for( int n=0; n<(len-1); n++ )
    {
        // If not a coverage tag
        if( std::strcmp(a[n].Type,"MCOV") )
        {
            // If positions are identical
            if( a[n].Position[0] == a[n+1].Position[0] )
            {
                if( std::strcmp(a[n].Type,"MUTA") == 0 )
                    a[n].Marked = 1;
                else if( std::strcmp(a[n+1].Type,"MUTA") == 0 )
                    a[n+1].Marked = 1;
                n++;
            }
        }
    }


    // Remove unwanted tags
    int dst = 0;
    for( int n=0; n<len; n++ )
    {
        if( a[n].Marked )
            continue;
        a[dst++] = a[n];
    }
    a.Length( dst );
}
