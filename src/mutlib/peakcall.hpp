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



#ifndef _MUTLIB_PEAKCALL_HPP_
#define _MUTLIB_PEAKCALL_HPP_



#include <listitem.hpp>     // For ListItem template



// Raw peak-call data
typedef struct
{
    int Width[4];
    int Position[4];
    int Amplitude[4];
    int BaseNumber;
    int BasePosition;

}peak_call_t;




class PeakCall : public ListItem<PeakCall>
{
/*
   Used as a general purpose container for peak-call information.
*/
 public:
    // Data
    peak_call_t Data;



 public:
    // Constructors
    PeakCall() { Init(); }
    PeakCall( int a, int c, int g, int t, int bnum=-1, int bpos=-1 );



 public:
    // Services
    bool IsValid() const;
    int  MaxWidthAsIndex() const;
    int  MaxAmplitudeAsIndex() const;
    int  MinAmplitudeAsIndex() const;



 private:
   // Helpers
   void Init();
};



#endif
