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


#ifndef _TRACEDIFF_PARAMETERS_HPP_
#define _TRACEDIFF_PARAMETERS_HPP_


#include <cassert>
#include <mutlib.h>
#include <parameter.hpp>



class TraceDiffParameters
{
 public:
    // Constructor/Destructor
    TraceDiffParameters();
   ~TraceDiffParameters();



 public:
     // Operators
     NumericParameter<double>& operator[]( int n ) 
        { assert(n<TRACEDIFF_PARAMETERS); return *(m_pParameter[n]); }


 private:
    // Data
    NumericParameter<double>* m_pParameter[TRACEDIFF_PARAMETERS];
};



#endif
