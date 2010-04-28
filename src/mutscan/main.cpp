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

#include <staden_config.h>

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
#include <string.h>	    // Needed for Solaris C++ compiler?



//-------------
// Module Data
//-------------

const int BUFSIZE = 512;




//------------------
// Helper Functions
//------------------

void PrintBanner()
{
    std::fprintf( stdout, "Mutscan v1.00, by Mark Jordan.\n"
                  "MRC Laboratory of Molecular Biology, Copyright (C) 2002. All Rights Reserved.\n" );
    std::fflush( stdout );
}



void PrintUsage( mutscan_t& ms )
{
    std::fprintf( stdout, "Usage  : mutscan [options] <experiment-files>\n"
    "Options:\n"
    "[-a<mutation-count>]     = Alignment failure report threshold   (default=%0.0f)\n"
    "[-c]                     = Complement reverse strand tags       (default=off)\n"
    "[-f<file-of-filenames>]  = File of experiment filenames\n"
    "[-h<heterzygote-SNR>]    = Heterzygote SNR threshold in dB      (default=%0.2f)\n"
    "[-l<peakdrop-threshold>] = Lower peak drop threshold percentage (default=%0.2f)\n"
    "[-n<noise-threshold]     = Noise level threshold percentage     (default=%0.2f)\n"
    "[-q]                     = Quiet mode                           (default=off)\n"
    "[-u<peakdrop-threshold>] = Upper peak drop threshold percentage (default=%0.2f)\n"
    "[-w<search-window-size>] = Peak search window size in bases     (default=%0.2f)\n"
    "[-p]                     = Proximity filter threshold...        (default=7)\n"
    "[-z]                     = Enter debugging loop and wait...     (default=off)\n",
     ms.Parameter[MUTSCAN_PARAMETER_ALIGNFAIL_THRESHOLD],
     ms.Parameter[MUTSCAN_PARAMETER_HETSNR_THRESHOLD],
     ms.Parameter[MUTSCAN_PARAMETER_PEAKDROP_LOWER],
     ms.Parameter[MUTSCAN_PARAMETER_NOISE_THRESHOLD],
     ms.Parameter[MUTSCAN_PARAMETER_PEAKDROP_UPPER],
     ms.Parameter[MUTSCAN_PARAMETER_SEARCH_WINDOW]
    );
    std::fflush( stdout );
}

// Filter clusters of tags too close to the end of MCOV
void filter_tags ( mutscan_t& ms, int nTags, int threshold ) {
    int *scores_a, *scores;
    int cov_start = 0, cov_end = 0;
    const int win_len = 5;

    // Get coverage range
    for( int i=0; i<nTags; i++ ) {
        mutlib_tag_t* pTag = MutScanGetTag( &ms, i );
        assert(pTag != NULL);

	if (std::strcmp(pTag->Type,"MCOV") == 0) {
	    cov_start = pTag->Position[0];
	    cov_end   = pTag->Position[1];
	    break;
	}
    }
    printf("Cov=%d..%d\n", cov_start, cov_end);


    /*
     * Allocate an array indexed from -win_len to cov_end+win_len.
     * Each tag gets a peak of width 2*win_len+1 centred on it and is
     * accumulated in this array.
     */
    scores_a = new int[win_len*2+cov_end+1];
    scores = &scores_a[win_len];
    memset(scores_a, 0, (win_len*2+cov_end) *sizeof(*scores_a));

    // Convolve tag peaks
    for (int i=0; i < nTags; i++) {
	mutlib_tag_t* pTag = MutScanGetTag( &ms, i );
	int pos = pTag->Position[0];
	if (pos < 0) pos = 0;
	if (pos > cov_end) pos = cov_end;

	scores[pos] += win_len+1;
	for (int j=1; j<=win_len; j++) {
	    scores[pos+j] += win_len+1 - j;
	    scores[pos-j] += win_len+1 - j;
	}
    }

    // Filter tags
    for (int i=0; i < nTags; i++) {
	mutlib_tag_t* pTag = MutScanGetTag( &ms, i );
	int pos = pTag->Position[0];

	if (std::strcmp(pTag->Type,"MCOV") == 0)
	    continue;

	if (scores[pos] >= threshold) {
	    pTag->Type[3] = '?';
	}
    }
    
    delete [] scores_a;

#if 0
    int last = 0;
    int cost = 0, highest_cost = 0;
    int range_start = 0, range_end = 0;
    for ( int i=0; i<nTags; i++) {
        mutlib_tag_t* pTag = MutScanGetTag( &ms, i );

	cost += 4;
	if (i != range_start)
	    cost -= pTag->Position[0] - last;
	printf("Tag %d at %d, cost = %d\n", i, pTag->Position[0], cost);
	if (highest_cost <= cost) {
	    highest_cost = cost;
	    range_end = i;
	}
	
	if ((cost < highest_cost-10 || cost < 0) && highest_cost > 10) {
	    printf("    Tag %d..%d highest_cost=%d\n",
		   range_start, range_end, highest_cost);
	    i = range_end;
	    range_start = i+1;
	    highest_cost = cost = 0;
	}
	last = pTag->Position[0];
    }
    if (highest_cost > 10) {
	printf("    Tag %d..%d highest_cost=%d\n",
	       range_start, range_end, highest_cost);
    }
#endif

#if 0
    // Adjust left end
    int cost = -5, last = cov_start;
    puts("==Left==");
    for( int i=0; i<nTags; i++ ) {
        mutlib_tag_t* pTag = MutScanGetTag( &ms, i );
        assert(pTag != NULL);
	if (std::strcmp(pTag->Type,"MCOV") == 0)
	    continue;

	cost += pTag->Position[0] - last;
	printf("Tag at %d cost=%d\n", pTag->Position[0], cost);
	if (cost <= 0) {
	    printf("Reject tag at %d\n", pTag->Position[0]);
	    *pTag->Type = 0;
	    cost -= 1;
	}
	last = pTag->Position[0];
    }

    // Adjust right end
    cost = -5;
    last = cov_end;
    puts("==Right==");
    for( int i=nTags-1; i>=0; i-- ) {
        mutlib_tag_t* pTag = MutScanGetTag( &ms, i );
        assert(pTag);
	if (std::strcmp(pTag->Type,"MCOV") == 0)
	    continue;

	cost += last - pTag->Position[0];
	printf("Tag at %d cost=%d\n", pTag->Position[0], cost);
	if (cost <= 0) {
	    printf("Reject tag at %d\n", pTag->Position[0]);
	    *pTag->Type = 0;
	    cost -= 2;
	}
	last = pTag->Position[0];
    }
#endif
}

//-----------
// Tracediff
//-----------

int main( int argc, char* argv[] )
{
    int             n;
    char*           p;
    mutscan_t       ms;
    mutlib_strand_t nStrand;
    StringList      FileList;
    int             nInputClipL;
    int             nInputClipR;
    char            pBuffer[BUFSIZE];
    char            pFileOfFiles[BUFSIZE]     = { 0 };
    Exp_info*       pExpFile                  = 0;
    Read*           pInputTrace               = 0;
    int             nRefClipL[2]              = { -1, -1 };
    int             nRefClipR[2]              = { -1, -1 };
    Read*           pRefTrace[2]              = { 0,   0 };
    char*           pRefSequence[2]           = { 0,   0 };
    bool            bDebug                    = false;
    bool            bQuiet                    = false;
    bool            bComplementTags           = false;
    double          nNoiseThreshold           = -1.0;
    double          nAlignFailThreshold       = -1.0;
    double          nSearchWindowSize         = -1.0;
    double          nPeakDropThresholdUpper   = -1.0;
    double          nPeakDropThresholdLower   = -1.0;
    double          nHeterozygoteSNRThreshold = -1.0;
    int 	    proximityThreshold = 7;
    MutScanInit( &ms );



    try
    {
        // Request traces to be used in preference to experiment files
        read_experiment_redirect(0);

        // Quick check of input
        if( argc < 2 )
        {
            PrintBanner();
            PrintUsage( ms );
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
                    // Alignment failure threshold
                    nAlignFailThreshold = std::atof( pBuffer );
                    break;


                case 'c':
                    // Complement reverse strand base tags
                    bComplementTags = true;
                    break;


                case 'f':
                    // File of filenames
                    std::strcpy( pFileOfFiles, pBuffer );
                    break;


                case 'h':
                    // Heterozygote SNR ratio threshold
                    nHeterozygoteSNRThreshold = std::atof( pBuffer );
                    break;


                case 'l':
                    // Lower peakdrop threshold
                    nPeakDropThresholdLower = std::atof( pBuffer );
                    break;


                case 'n':
                    // Noise window
                    nNoiseThreshold = std::atof( pBuffer );
                    break;


                case 'q':
                    // Quiet mode on
                    bQuiet = true;
                    break;


                case 'u':
                    // Upper peakdrop threshold
                    nPeakDropThresholdUpper = std::atof( pBuffer );
                    break;


                case 'w':
                    // Search window size
                    nSearchWindowSize = std::atof( pBuffer );
                    break;


                case 't':
                    // proximity threshold
                    proximityThreshold = std::atoi( pBuffer );
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
            PrintUsage( ms );
            return -1;
        }



        // Set parameters if specified, otherwise accept defaults
        if( nAlignFailThreshold > 0.0 )
            MutScanSetParameter( &ms, MUTSCAN_PARAMETER_ALIGNFAIL_THRESHOLD, nAlignFailThreshold );
        if( nHeterozygoteSNRThreshold > 0.0 )
            MutScanSetParameter( &ms, MUTSCAN_PARAMETER_HETSNR_THRESHOLD, nHeterozygoteSNRThreshold );
        if( nNoiseThreshold > 0.0 )
            MutScanSetParameter( &ms, MUTSCAN_PARAMETER_NOISE_THRESHOLD, nNoiseThreshold );
        if( nSearchWindowSize > 0.0 )
            MutScanSetParameter( &ms, MUTSCAN_PARAMETER_SEARCH_WINDOW, nSearchWindowSize );
        if( nPeakDropThresholdLower > 0.0 )
            MutScanSetParameter( &ms, MUTSCAN_PARAMETER_PEAKDROP_LOWER, nPeakDropThresholdLower );
        if( nPeakDropThresholdUpper > 0.0 )
            MutScanSetParameter( &ms, MUTSCAN_PARAMETER_PEAKDROP_UPPER, nPeakDropThresholdUpper );
        if( bComplementTags == true )
            MutScanSetParameter( &ms, MUTSCAN_PARAMETER_COMPLEMENT_TAGS, 1.0 );



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
                std::fprintf( stdout, "Scanning: %s\n", p );
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



            // If we don't have a reference trace/seq for this strand
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
                MutScanSetReference( &ms, nStrand, pRefTrace[nStrand], nRefClipL[nStrand], nRefClipR[nStrand] );
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
	    if ( !pInputTrace->NPoints ) {
		std::fprintf(stderr, "No trace data found in trace file %s, skipping.\n", pBuffer);
		continue;
	    }



            // Extract input clip points QL/QR
            nInputClipL = -1;
            nInputClipR = -1;
            exp_get_int( pExpFile, EFLT_QL, &nInputClipL );
            exp_get_int( pExpFile, EFLT_QR, &nInputClipR );



            // Initialise the algorithm
            MutScanSetInput( &ms, nStrand, pInputTrace, nInputClipL, nInputClipR );



            // Execute the algorithm
            MutScanExecute( &ms );
            if( MutScanGetResultCode(&ms) )
            {
                // Error
                std::fprintf( stderr, "%s", MutScanGetResultString(&ms) );
                std::fflush( stderr );
                continue;
            }
            int nTags = MutScanGetTagCount( &ms );



            // Find the rightmost tag position for this trace
            int nRightmostTag = -1;
            for( int i=0; i<nTags; i++ )
            {
                mutlib_tag_t* pTag = MutScanGetTag( &ms, i );
                if( pTag->Position[0] > nRightmostTag )
                    nRightmostTag = pTag->Position[0];
            }



            // If mutation tag is beyond the right clip point (as in the case of
            // insertions/deletions), we adjust QR so that it's visible in gap4.
            if( (nTags>0) && (nInputClipR<=nRightmostTag) )
            {
                std::sprintf( pBuffer, "%d", nRightmostTag+1 );
                exp_put_str( pExpFile, EFLT_QR, pBuffer, std::strlen(pBuffer) );
            }

	    // Filter clusters of tags too close to the end of MCOV
	    filter_tags( ms, nTags, proximityThreshold );

            // Output results
            for( int i=0; i<nTags; i++ )
            {
                // Write mutation tags to experiment file & stdout
                mutlib_tag_t* pTag = MutScanGetTag( &ms, i );
                assert(pTag != NULL);

		if (!*pTag->Type)
		  continue;

                bool bCoverageTag = std::strcmp(pTag->Type,"MCOV") == 0;
                char sc = (pTag->Strand==MUTLIB_STRAND_FORWARD) ? '+' : '-';
                if( bCoverageTag )
                {
                    std::sprintf( pBuffer, "%s %c %d..%d", pTag->Type, sc, pTag->Position[0], pTag->Position[1] );
                    if( !bQuiet )
                    {
                        std::fprintf( stdout, "%s\n", pBuffer );
                        std::fflush( stdout );
                    }
                }
                else
                {
                    std::sprintf( pBuffer, "%s %c %d..%d\n%s", pTag->Type, sc, pTag->Position[0], pTag->Position[1], pTag->Comment );
                    if( !bQuiet )
                    {
                        std::fprintf( stdout, "%s %5d %s\n", pTag->Type, pTag->Position[0], pTag->Comment );
                        std::fflush( stdout );
                    }
                }
                exp_put_str(pExpFile, EFLT_TG, pBuffer, std::strlen(pBuffer) );
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
    MutScanDestroy( &ms );
    return 0;
}
