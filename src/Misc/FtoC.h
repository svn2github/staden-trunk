/*
    Title: 	 FtoC

    File: 	 FtoC.h
    Purpose:	 FORTRAN-C string conversion routines
    Last update: Mon Jun 18 1990
*/

#ifndef _FTOC_H_
#define _FTOC_H_

#include "os.h"


extern void Cstr2Fstr(char *Cstr,
		      char *Fstr, int_fl Flen);
/*
    Copy a '\0' terminated C string to a Fortran string, blank padding
    if needed and ignoring excess C characters if needed.

    This function works if the strings are distinct or coincident, but
    not if they overlap in any other way.
*/


extern void Fstr2Cstr(char *Fstr, int_fl Flen,
		      char *Cstr, int_fl Clen);
/*
    Copy the significant characters of a blank padded Fortran string
    to a '\0' terminated C string, ignoring excess characters.

    This function works if the strings are distinct or coincident, but
    not if they overlap in any other way.
*/

extern int_f swapbo_(int_f *i4);
/*
 * Returns the big-endian form of a four byte integer
 */

#endif /*_FTOC_H_*/
