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



#include <cstring>         // For strlen(), strcpy(), strcat(), strrchr()



void MakeFullPath( const char* pExpFile, char* pTraceFile )
{
/*
    This routine assumes that the trace file and the experiment file are
    in the same directory. We assemble a full trace file name by using the
    path of the experiment file and appending the trace file filename 
    component. Any path on the tracefile is ignored since we cannot tell
    if it's relative or absolute. The routine should work in windows and
    unix since both forward and reverse slashes are handled.
*/
    char pBuffer[512];



    // Point to trace filename component (without path)
    char* p1 = std::strrchr(pTraceFile,'/');
    if(!p1)
    {
        p1 = std::strrchr(pTraceFile,'\\');
        if(!p1)
            p1 = pTraceFile - 1;
    }
    p1++;



    // Point to experiment filename component (without path)
    std::strcpy( pBuffer, pExpFile );
    char* p2 = std::strrchr(pBuffer,'/');
    if(!p2)
    {
        p2 = std::strrchr(pBuffer,'\\');
        if(!p2)
            p2 = pBuffer - 1;
    }
    p2++;



    // Make replacement
    std::strcpy( p2, p1 );
    std::strcpy( pTraceFile, pBuffer );
}




void ReplaceExtension( char* pFileName, const char* pNewExt )
{
/*
    This routine replaces the file extension contained in 'pFileName' with
    'pNewExt'. If there is no existing extension, we append 'pNewExt'.
*/

    char* pExt = std::strrchr( pFileName, '.' );
    if( pExt )
        std::strcpy( pExt, pNewExt );
    else
        std::strcat( pFileName, pNewExt );
}


