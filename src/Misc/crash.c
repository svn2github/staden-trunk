#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>  /* varargs needed for v*printf() prototypes */

void crash (char* format,...)
{
    va_list args ;

    va_start (args,format) ;
    vfprintf (stderr,format,args) ;
    va_end (args) ;

    exit (1) ;
}
