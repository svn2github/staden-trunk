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


#include <mutscan_parameters.hpp>


//-------------
// Constructor
//-------------

MutScanParameters::MutScanParameters()
{
   for( int n=0; n<MUTSCAN_PARAMETERS; n++ )
       m_pParameter[n] = 0;
   m_pParameter[MUTSCAN_PARAMETER_ALIGNFAIL_THRESHOLD] = new NumericParameter<double>( 15.0, 2.0, 100.0, "alignment failure threshold" );
   m_pParameter[MUTSCAN_PARAMETER_COMPLEMENT_TAGS]     = new NumericParameter<double>( -1.0, -2.0, 2.0, "complement strand reverse tags" );
   m_pParameter[MUTSCAN_PARAMETER_HETSNR_THRESHOLD]    = new NumericParameter<double>( 7.6, 1.0, 60.0, "heterozygote SNR threshold" );
   m_pParameter[MUTSCAN_PARAMETER_PEAKDROP_LOWER]      = new NumericParameter<double>( 0.2, 0.01, 1.0, "lower peak drop threshold" );
   m_pParameter[MUTSCAN_PARAMETER_NOISE_THRESHOLD]     = new NumericParameter<double>( 0.25, 0.01, 1.0, "noise threshold" );
   m_pParameter[MUTSCAN_PARAMETER_PEAKDROP_UPPER]      = new NumericParameter<double>( 0.7, 0.01, 1.0, "upper peak drop threshold" );
   m_pParameter[MUTSCAN_PARAMETER_SEARCH_WINDOW]       = new NumericParameter<double>( 0.9, 0.01, 1.5, "search window" );
}



//------------
// Destructor
//------------

MutScanParameters::~MutScanParameters()
{
   for( int n=0; n<MUTSCAN_PARAMETERS; n++ )
       delete m_pParameter[n];
}
