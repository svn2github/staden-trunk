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



#ifndef _MUTLIB_VALIDATE_HPP_
#define _MUTLIB_VALIDATE_HPP_


#include <mutlib.h>


mutlib_result_t MutlibValidateTrace( mutlib_trace_t& td, char* rs, const char* s );
mutlib_result_t MutlibValidateTraceClipPoints( mutlib_trace_t& td, char* rs, const char* s );


#endif

