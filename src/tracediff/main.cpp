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
#include <new>              // For bad_alloc
#include <cstdio>           // For fprintf(), fflush(), etc
#include <cctype>           // For isspace(), toupper()
#include <cstring>          // For strcpy(), strcat(), etc
#include <cstdlib>          // For atof()
#include <staden.h>         // For io_lib and licence stuff
#include <mutlib.h>         // For mutation library
#include <pathutil.h>       // For MakeFullPath()
#include <stringlist.hpp>   // For StringList




//-------------
// Module Data
//-------------

const int BUFSIZE = 512;




//------------------
// Helper Functions
//------------------

void PrintBanner()
{
    std::fprintf( stdout, "Tracediff v2.00.\n"
                  "MRC Laboratory of Molecular Biology, Copyright (C) 2002. All Rights Reserved.\n" );
    std::fflush( stdout );
}



void PrintUsage( tracediff_t& td )
{
    std::fprintf( stdout, "Usage  : tracediff [options] <experiment-files>\n"
    "Options:\n"
    "[-a<peak-alignment>]    = Maximum peak alignment deviation (default=%0.2f)\n"
    "[-c]                    = Complement reverse strand tags   (default=off)\n"
    "[-d]                    = Output difference traces         (default=off)\n"
    "[-f<file-of-filenames>] = File of experiment filenames\n"
    "[-n<noise-window>]      = Background noise window length   (default=%0.0f)\n"
    "[-q]                    = Quiet mode                       (default=off)\n"
    "[-s<sensitivity>]       = Mutation detection sensitivity   (default=%0.2f)\n"
    "[-t<noise-threshold]    = Noise level threshold percentage (default=%0.2f)\n"
    "[-w<peak-width-max>]    = Maximum peak width in bases      (default=%0.2f)\n"
    "[-y]                    = Y-scale traces before difference (default=off)\n"
    "[-z]                    = Enter debugging loop and wait... (default=off)\n",
     td.Parameter[TRACEDIFF_PARAMETER_PEAK_ALIGNMENT],
     td.Parameter[TRACEDIFF_PARAMETER_NOISE_WINDOW_LENGTH],
     td.Parameter[TRACEDIFF_PARAMETER_SENSITIVITY],
     td.Parameter[TRACEDIFF_PARAMETER_NOISE_THRESHOLD],
     td.Parameter[TRACEDIFF_PARAMETER_PEAK_WIDTH_MAXIMUM]
    );
    std::fflush( stdout );
}



//-----------
// Tracediff
//-----------

int main( int argc, char* argv[] )
{
    int             n;
    char*           p;
    tracediff_t     td;
    mutlib_strand_t nStrand;
    StringList      FileList;
    int             nInputClipL;
    int             nInputClipR;
    char            pBuffer[BUFSIZE];
    char            pFileOfFiles[BUFSIZE]  = { 0 };
    Exp_info*       pExpFile               = 0;
    Read*           pInputTrace            = 0;
    Read*           pRefTrace[2]           = { 0,   0 };
    int             nRefClipL[2]           = { -1, -1 };
    int             nRefClipR[2]           = { -1, -1 };
    bool            bDebug                 = false;
    bool            bQuiet                 = false;
    bool            bYScaleTraces          = false;
    bool            bComplementTags        = false;
    bool            bOutputDifferenceTrace = false;
    double          nSensitivity           = -1.0;
    double          nPeakAlignment         = -1.0;
    double          nPeakWidthMaximum      = -1.0;
    double          nNoiseWindowLength     = -1.0;
    double          nNoiseThreshold        = -1.0;
    TraceDiffInit( &td );



    try
    {
        // Request traces to be used in preference to experiment files
        read_experiment_redirect(0);

        // Quick check of input
        if( argc < 2 )
        {
            PrintBanner();
            PrintUsage( td );
            return -1;
        }



        // Parse the command line options
        for( n=1; n<argc; n++ )
        {
            // All options done?
            if( argv[n][0] != '-' )
                break;


            // Copy switch contents to buffer
            char c;
            int  k = 2;
            while( c=argv[n][k] )
            {
                pBuffer[k-2] = c;
                k++;
            }
            pBuffer[k-2] = 0;



            // Process each option
            switch( argv[n][1] )
            {
                case 'a':
                    // Peak alignment threshold
                    nPeakAlignment = std::atof( pBuffer );
                    break;


                case 'c':
                    // Complement reverse strand base tags
                    bComplementTags = true;
                    break;


                case 'd':
                    // Difference traces
                    bOutputDifferenceTrace = true;
                    break;


                case 'f':
                    // File of filenames
                    std::strcpy( pFileOfFiles, pBuffer );
                    break;


                case 'n':
                    // Noise window
                    nNoiseWindowLength = std::atof( pBuffer );
                    break;


                case 'q':
                    // Quiet mode on
                    bQuiet = true;
                    break;


                case 's':
                    // Sensitivity
                    nSensitivity = std::atof( pBuffer );
                    break;


                case 't':
                    // Noise threshold
                    nNoiseThreshold = std::atof( pBuffer );
                    break;


                case 'w':
                    // Peak width maximum
                    nPeakWidthMaximum = std::atof( pBuffer );
                    break;


                case 'y':
                    // Y-scale traces
                    bYScaleTraces = true;
                    break;


                case 'z':
                    // Debug mode
                    bDebug = true;
                    break;


                default:
                    // Unknown switch
                    std::fprintf( stderr, "Unknown option %s, ignoring.", argv[n] );
                    std::fflush( stderr );
                    break;
            }
        }



        // Show banner
        if( !bQuiet )
            PrintBanner();



        // Loop for MSW debugger attachment purposes
        while( bDebug );



        // Read contents of file-of-filenames if supplied
        if( pFileOfFiles[0] )
        {
            std::FILE* pFOFN = std::fopen( pFileOfFiles, "rt" );
            if( pFOFN )
            {
                while( std::fgets(pBuffer,BUFSIZE,pFOFN) )
                {
                    p = std::strrchr(pBuffer,'\n');
                    if(p) *p = 0;
                    p = std::strrchr(pBuffer,'\r');
                    if(p) *p = 0;
                    FileList.Append( pBuffer );
                }
                std::fclose( pFOFN );
            }
            else
            {
                std::fprintf( stderr, "Unable to open file of filenames %s, skipping.\n", pFileOfFiles );
                std::fflush( stderr );
            }
        }



        // Add any experiment files on the command line
        for( n=n; n<argc; n++ )
            FileList.Append( argv[n] );
        if( FileList.Length() < 1 )
        {
            PrintUsage( td );
            return -1;
        }



        // Set parameters if specified, otherwise accept defaults
        if( nSensitivity > 0.0 )
            TraceDiffSetParameter( &td, TRACEDIFF_PARAMETER_SENSITIVITY, nSensitivity );
        if( nNoiseThreshold > 0.0 )
            TraceDiffSetParameter( &td, TRACEDIFF_PARAMETER_NOISE_THRESHOLD, nNoiseThreshold );
        if( nNoiseWindowLength > 0.0 )
            TraceDiffSetParameter( &td, TRACEDIFF_PARAMETER_NOISE_WINDOW_LENGTH, nNoiseWindowLength );
        if( nPeakAlignment > 0.0 )
            TraceDiffSetParameter( &td, TRACEDIFF_PARAMETER_PEAK_ALIGNMENT, nPeakAlignment );
        if( nPeakWidthMaximum > 0.0 )
            TraceDiffSetParameter( &td, TRACEDIFF_PARAMETER_PEAK_WIDTH_MAXIMUM, nPeakWidthMaximum );
        if( bYScaleTraces == true )
            TraceDiffSetParameter( &td, TRACEDIFF_PARAMETER_YSCALE, 1.0 );
        if( bComplementTags == true )
            TraceDiffSetParameter( &td, TRACEDIFF_PARAMETER_COMPLEMENT_TAGS, 1.0 );



        // File processing loop
        for( n=0; n<FileList.Length(); n++ )
        {
            // Cleanup from previous iteration
            if( pExpFile )
            {
                exp_destroy_info(pExpFile);
                pExpFile = 0;
            }
            if( pInputTrace )
            {
                read_deallocate(pInputTrace);
                pInputTrace = 0;
            }
            p = (n==0) ? FileList.First() : FileList.Next();



            // Tell user what we're doing
            if( !bQuiet )
            {
                std::fprintf( stdout, "Processing: %s\n", p );
                std::fflush( stdout );
            }



            // Open experiment file
            pExpFile = exp_read_info( p );
            if( !pExpFile )
            {
                std::fprintf( stderr, "Unable to open experiment file %s, skipping.\n", p );
                std::fflush( stderr );
                continue;
            }



            // Extract PR record, gives us strand direction
            int PR = 0;
            if( exp_get_int(pExpFile,EFLT_PR,&PR) )
            {
                std::fprintf( stderr, "Unable to read PR record from experiment file %s, assuming forward strand.\n", p );
                std::fflush( stderr );
            }
            nStrand = ((PR==0) || (PR==1) || (PR==3)) ? MUTLIB_STRAND_FORWARD : MUTLIB_STRAND_REVERSE;



            // If we don't have a reference trace for this strand
            if( !pRefTrace[nStrand] )
            {
                // Extract reference trace name WT
                if( exp_get_str(pExpFile,EFLT_WT,pBuffer,BUFSIZE) )
                {
                    std::fprintf( stderr, "Unable to read WT record from experiment file %s, skipping.\n", p );
                    std::fflush( stderr );
                    continue;
                }


                // Open reference trace file
                MakeFullPath( p, pBuffer );
                pRefTrace[nStrand] = read_reading( pBuffer, TT_ANY );
                if( !pRefTrace[nStrand] )
                {
                    std::fprintf( stderr, "Unable to open reference trace file %s, skipping.\n", pBuffer );
                    std::fflush( stderr );
                    continue;
                }



                // Extract reference clip points WL/WR
                exp_get_int( pExpFile, EFLT_WL, &nRefClipL[nStrand] );
                exp_get_int( pExpFile, EFLT_WR, &nRefClipR[nStrand] );



                // If default clipping values used for reference
                if( (nRefClipL[nStrand]<0) || (nRefClipR[nStrand]<0) )
                {
                    // Construct the reference experiment file's name and attempt to open it
                    ReplaceExtension( pBuffer, ".exp" );
                    Exp_info* pRefExpFile = exp_read_info( pBuffer );
                    if( pRefExpFile )
                    {
                        // Extract the QL/QR records if they exist
                        if( nRefClipL[nStrand] < 0 )
                            exp_get_int( pRefExpFile, EFLT_QL, &nRefClipL[nStrand] );
                        if( nRefClipR[nStrand] < 0 )
                            exp_get_int( pRefExpFile, EFLT_QR, &nRefClipR[nStrand] );
                        exp_destroy_info(pRefExpFile);
                    }
                    // Otherwise use the entire reference trace
                }


                // Initialise the algorithm
                TraceDiffSetReference( &td, pRefTrace[nStrand], nStrand, nRefClipL[nStrand], nRefClipR[nStrand] );
            }



            // Extract LN record, gives us input trace file name
            if( exp_get_str(pExpFile,EFLT_LN,pBuffer,BUFSIZE) )
            {
                std::fprintf( stderr, "Unable to read LN record from experiment file %s, skipping.\n", p );
                std::fflush( stderr );
                continue;
            }


            // Open input trace
            MakeFullPath( p, pBuffer );
            pInputTrace = read_reading( pBuffer, TT_ANY );
            if( !pInputTrace )
            {
                std::fprintf( stderr, "Unable to open input trace file %s, skipping.\n", pBuffer );
                std::fflush( stderr );
                continue;
            }



            // Extract input clip points QL/QR
            nInputClipL = -1;
            nInputClipR = -1;
            exp_get_int( pExpFile, EFLT_QL, &nInputClipL );
            exp_get_int( pExpFile, EFLT_QR, &nInputClipR );



            // Initialise the algorithm
            TraceDiffSetInput( &td, pInputTrace, nStrand, nInputClipL, nInputClipR );



            // Execute the algorithm
            TraceDiffExecute( &td, TRACEDIFF_ALGORITHM_DEFAULT );
            if( TraceDiffGetResultCode(&td) )
            {
                std::fprintf( stderr, "%s", TraceDiffGetResultString(&td) );
                std::fflush( stderr );
                continue;
            }



            // Output Difference Trace
            if( bOutputDifferenceTrace )
            {
                Read* pDiff = TraceDiffGetDifference( &td, 0, 0 );
                std::strcpy( pBuffer, p );
                ReplaceExtension( pBuffer, "_diff.ztr" );
                int retval = write_reading( pBuffer, pDiff, TT_ZTR );
                if( retval < 0 )
                {
                    std::fprintf( stderr, "Unable to write out difference trace %s.\n", pBuffer );
                    std::fflush( stderr );
                }
            }



            // Output results
            int nTags = TraceDiffGetTagCount( &td );
            for( int i=0; i<nTags; i++ )
            {
                // Write mutation tags to experiment file
                mutlib_tag_t* pTag = TraceDiffGetTag( &td, i );
                assert(pTag != NULL);
                char c = (pTag->Strand==MUTLIB_STRAND_FORWARD) ? '+' : '-';
                std::sprintf( pBuffer, "%s %c %d..%d\n%s", pTag->Type, c, *pTag->Position,
                              *pTag->Position, pTag->Comment );
                exp_put_str(pExpFile, EFLT_TG, pBuffer, std::strlen(pBuffer) );
                if( !bQuiet )
                {
                    std::fprintf( stdout, "%s %5d %s\n", pTag->Type, *pTag->Position, pTag->Comment );
                    std::fflush( stdout );
                }
            }
        }
    }
    catch( std::bad_alloc& )
    {
        // Out of memory
        std::fprintf( stderr, "%s", "Not enough memory to complete the operation.\n" );
        std::fflush( stderr );
    }
    catch(...)
    {
        // Print error message
        std::fprintf( stderr, "%s", "An unexpected fatal exception has occurred, please "
                              "report the details to staden-package@mrc-lmb.cam.ac.uk.\n" );
        std::fflush( stderr );
    }




    // Cleanup & exit
    if(pExpFile)     exp_destroy_info(pExpFile);
    if(pInputTrace)  read_deallocate(pInputTrace);
    if(pRefTrace[0]) read_deallocate(pRefTrace[0]);
    if(pRefTrace[1]) read_deallocate(pRefTrace[1]);
    TraceDiffDestroy( &td );
    return 0;
}
