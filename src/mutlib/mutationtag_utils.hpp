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


#ifndef _MUTLIB_MUTATIONTAG_UTILS_HPP_
#define _MUTLIB_MUTATIONTAG_UTILS_HPP_


#include <mutlib.h>
#include <list.hpp>
#include <array.hpp>
#include <mutationtag.hpp>


void CompTags( SimpleArray<mutlib_tag_t>& a );
void SortTags( SimpleArray<mutlib_tag_t>& a );
void PruneTags( SimpleArray<mutlib_tag_t>& a );
void CopyTags( SimpleArray<mutlib_tag_t>& a, List<MutationTag>& l );



#endif
