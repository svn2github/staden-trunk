#include "os.h"

void Cstr2Fstr(char *Cstr,
	       char *Fstr, int_fl Flen)
/*
    Copy a '\0' terminated C string to a Fortran string, blank padding
    if needed and ignoring excess C characters if needed.

    This function works if the strings are distinct or coincident, but
    not if they overlap in any other way.
*/
{   int_fl i;

    for (i=0; (i<Flen) && (Cstr[i] != '\0'); i++)
    {   Fstr[i] = Cstr[i];
    }
    for (; i<Flen; i++)
    {   Fstr[i] = ' ';
    }
}




void Fstr2Cstr(char *Fstr, int_fl Flen,
	       char *Cstr, int_fl Clen)
/*
    Copy the significant characters of a blank padded Fortran string
    to a '\0' terminated C string, ignoring excess characters.

    This function works if the strings are distinct or coincident, but
    not if they overlap in any other way.
*/
{
    int_fl FsigLen, i, j;

    /*
     * We need to scan from the start to the end as we cannot trust Flen to
     * really be correct.
     * Eg consider the case of reading an old fortran string of length 16 from
     * the database, when the current code has an array of size 32.
     * For this reason, we work forwards noting the last space character.
     * If we reach Flen or NULL and the last character(s) were spaces, then we
     * trim back.
     */

    /* Find significant length of fortran string */
    for (j = i = 0; i < Flen && Fstr[i]; i++) {
	if (Fstr[i] == ' ')
	    j++;
	else
	    j = 0;
    }
    FsigLen = i - j;

    /*
     * Copy up to (Clen) significant characters, NULL terminating if
     * there is room.
     */
    i=0;
    while ((i < FsigLen) && (i < Clen)) {
	Cstr[i] = Fstr[i];
        i++;
    }

    if (i < Clen)
	Cstr[i] = '\0';
}



int_f swapbo_(int_f *i4)
/*
 * Returns the big-endian form of a four byte integer
 */
{
    int i=1;

    if (*(char*)&i) {

	int_f swapped;

	swap_int4(*i4,swapped);
	return swapped;
    } else
	return *i4;

}
