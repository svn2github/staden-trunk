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


#ifndef _MUTSCAN_HPP_
#define _MUTSCAN_HPP_


#include <mutlib.h>
#include <mutscan_parameters.hpp>       // For MutScanParameters object


void            MutScanDestroyCache( mutscan_t* ms );
void            MutScanDestroyResults( mutscan_t* ms );
mutlib_result_t MutScanValidateInput( mutscan_t* ms, MutScanParameters& p );


#endif
