#include "misc.h"
#include <string.h>

char *fn_tail(char *fn)
/*
** Return file part (:t) of
** directory path
*/
{
    int len;
    char *s;

    len = strlen(fn);
    for(s=fn+len-1;len && *s != '/'; len--, s--) ;
    s++;

    return s;
}


void fn_toupper (char *s)
/*
** Convert file to upper case
** ignoring directory path head
*/
{
    str_toupper(fn_tail(s));
}



void fn_tolower (char *s)
/*
** Convert file to lower case
** ignoring directory path head
*/
{
    str_tolower(fn_tail(s));
}
