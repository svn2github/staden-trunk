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



#ifndef _MUTLIB_PATHUTIL_H_
#define _MUTLIB_PATHUTIL_H_



void MakeFullPath( const char* pExpFile, char* pTraceFile );
void ReplaceExtension( char* pFileName, const char* pNewExt );



#endif


