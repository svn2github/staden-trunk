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


#include <tracediff_parameters.hpp>



//-------------
// Constructor
//-------------

TraceDiffParameters::TraceDiffParameters()
{
   for( int n=0; n<TRACEDIFF_PARAMETERS; n++ )
       m_pParameter[n] = 0;
   m_pParameter[TRACEDIFF_PARAMETER_SENSITIVITY]         = new NumericParameter<double>( 5.0, 1.0, 50.0, "sensitivity" );
   m_pParameter[TRACEDIFF_PARAMETER_NOISE_THRESHOLD]     = new NumericParameter<double>( 0.09, 0.01, 0.5, "noise threshold" );
   m_pParameter[TRACEDIFF_PARAMETER_NOISE_WINDOW_LENGTH] = new NumericParameter<double>( 12.0, 2.0, 50.0, "noise window length" );
   m_pParameter[TRACEDIFF_PARAMETER_PEAK_ALIGNMENT]      = new NumericParameter<double>( 0.4, 0.01, 2.0, "peak alignment" );
   m_pParameter[TRACEDIFF_PARAMETER_PEAK_WIDTH_MAXIMUM]  = new NumericParameter<double>( 2.1, 1.0, 4.0, "maximum peak width" );
   m_pParameter[TRACEDIFF_PARAMETER_COMPLEMENT_TAGS]     = new NumericParameter<double>( -1.0, -2.0, 2.0, "complement reverse tags" );
   m_pParameter[TRACEDIFF_PARAMETER_YSCALE]              = new NumericParameter<double>( -1.0, -2.0, 2.0, "y-scale traces" );
}



//------------
// Destructor
//------------

TraceDiffParameters::~TraceDiffParameters()
{
   for( int n=0; n<TRACEDIFF_PARAMETERS; n++ )
       delete m_pParameter[n];
}
