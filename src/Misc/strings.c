#include "misc.h"
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <sys/types.h>

int fstrlen(char *f, int max_f)
{
    for (; max_f > 0 && (isspace(f[max_f-1]) || f[max_f-1]=='\0'); max_f--);
    return max_f;
}





void f2cstr(char *f, int max_f, char *c, int max_c)
{
    int i;

    i = min(fstrlen(f,max_f),max_c);
    strncpy(c,f,i);
    c[i]='\0';
}


void c2fstr(char *c, int max_c, char *f, int max_f)
{
    int i;
    i = min((int)strlen(c),max_f);
    strncpy(f,c,i);
    for( ; i<max_f; i++) f[i]=' ';

}




char *mystrtok(char *s, char *ct)
/*
** When strtok isn't good enough
*/
{
    char *this;
    static char *look;
    static int last;

    if (s == NULL) {
	if (last) return NULL;
    } else {
	look = s;
	last = 0;
    }
    this = look;

    for ( ; *look && strchr(ct,*look)==NULL; look++ ) ;
    last = (! *look);
    *look++ = '\0';
    
    return this;
}


void str_tolower (char *s)
/*
** Convert string to lower case
*/
{
    if (!s) return;
    for ( ; *s ; s++ )
	if (isupper(*s))
	    *s = tolower(*s);
}

void str_toupper (char *s)
/*
** Convert string to upper case
*/
{
    if (!s) return;
    for ( ; *s ; s++ )
	if (islower(*s))
	    *s = toupper(*s);
}

#ifdef NOSTRSTR
/*
** My routines for nice sun ones.
*/
char *strstr(char *cs, char *ct)
/*
** ANSI C has the function strstr().
**
**     strstr() returns a pointer to the first  occurrence  of  the
**     pattern  string  s2  in  s1.   For example, if s1 is "string
**     thing" and s2 is "ing", strstr() returns "ing thing".  If s2
**     does not occur in s1, strstr() returns NULL.
**
** It's not always implemented. Here's my cludge:
*/
{
    int i;
    int len_ct;
    int end;
    len_ct = strlen(ct);
    end = strlen(cs) - len_ct;
    for (i=0;i<=end;i++)
      if (strncmp(&cs[i],ct,len_ct)==0)
	return &cs[i];

    return NULL;
}
#endif

#ifdef NOSTRDUP
char *strdup(const char *str)
/*
** SunOS has a nice strdup() function.
**
**     strdup() returns a pointer to a new string which is a dupli-
**     cate  of the string pointed to by s1.  The space for the new
**     string is obtained using malloc(3V).  If the new string can-
**     not be created, a NULL pointer is returned.
**
** Other ANSI C libraries don't have this. Here is my kludge:
*/
{
    char *newstr;
    int i = strlen(str);

    if ((newstr = (char *)malloc((unsigned int)(i+1))) == NULL)
        return NULL;

    for (; i>=0; i--)
        newstr[i] = str[i];

    return newstr;
}
#endif

#ifdef NOSTRCASECMP
int strcasecmp(const char *s1, const char *s2) {
    while (tolower(*s1) == tolower(*s2)) {
        /* If at the end of the string, then they're equal */
        if (0 == *s1)
	    return 0;
        s1++;
	s2++;
    }
  
    /* One ended before the other, so return 1 or -1 */
    return (*(unsigned char *)s1) < (*(unsigned char *)s2) ? -1 : 1;
}
#endif

#ifdef NOMEMMOVE
void *memmove(void *s1, const void *s2, size_t n) {
  bcopy(s2, s1, n);
  return s1;
}
#endif

/*
 * Strnlen: like strlen(), but with a maximum size so we can do strlen on
 * potentially non-null terminated arrays.
 */
size_t strnlen(const char *buf, size_t n) { 
    size_t i = 0; 
    while (i < n && buf[i]) 
        i++; 
    return i; 
} 

/*
 * Allocates and returns an escaped version of str. This relaces quotes,
 * newlines, and other non-printable characters with backslashed versions of
 * them in a C string style formatting.
 *
 * Returns malloced string on success
 *         NULL on failure.
 */
char *escape_C_string(char *str) {
    size_t l = strlen(str);
    size_t new_l = l*1.1+10;
    char *new = malloc(new_l);
    size_t oi, ni;
    static char type[256];
    static int type_init = 0;

    /* A once-only lookup table to speed up the loop below */
    if (!type_init) {
	int i;
	
	for (i = 0; i < 256; i++) {
	    if (isprint(i) && i != '"' && i != '\\') {
		/* directly printable */
		type[i] = 0;
	    } else {
		switch(i) {
		    /* backslash single-char */
		case '"':
		case '\\':
		    type[i] = i;
		    break;
		case '\n':
		    type[i] = 'n';
		    break;
		case '\r':
		    type[i] = 'r';
		    break;
		case '\t':
		    type[i] = 't';
		    break;
		case '\a':
		    type[i] = 'a';
		    break;

		default:
		    type[i] = 1; /* octal escape */
		}
	    }
	}
	type_init = 1;
    }


    if (!new)
	return NULL;

    for (oi = ni = 0; oi < l; oi++) {
	char c = str[oi];

	/* Make enough room */
	if (ni + 5 >= new_l) {
	    new_l = new_l * 1.2 + 10;
	    if (NULL == (new = realloc(new, new_l)))
		return NULL;
	}

	switch(type[(unsigned char)c]) {
	case 0:
	    new[ni++] = c;
	    break;
	    
	case 1:
	    sprintf(&new[ni], "\\%03o", c);
	    ni+=4;
	    break;

	default:
	    new[ni++] = '\\';
	    new[ni++] = type[(unsigned char)c];
	}
    }
    new[ni++] = 0;

    return new;
}

/*
 * As per escape_C_string but \n and \\ only.
 *
 * Returns malloced string on success
 *         NULL on failure.
 */
char *escape_C_nl(char *str) {
    size_t l = strlen(str);
    size_t new_l = l*1.1+10;
    char *new = malloc(new_l);
    size_t oi, ni;
    static char type[256];
    static int type_init = 0;

    /* A once-only lookup table to speed up the loop below */
    if (!type_init) {
	int i;
	
	for (i = 0; i < 256; i++) {
	    switch(i) {
		/* backslash single-char */
	    case '\\':
		type[i] = '\\';
		break;
	    case '\n':
		type[i] = 'n';
		break;
	    default:
		type[i] = 0;
	    }
	}
	type_init = 1;
    }


    if (!new)
	return NULL;

    for (oi = ni = 0; oi < l; oi++) {
	char c = str[oi];

	/* Make enough room */
	if (ni + 5 >= new_l) {
	    new_l = new_l * 1.2 + 10;
	    if (NULL == (new = realloc(new, new_l)))
		return NULL;
	}

	switch(type[(unsigned char)c]) {
	case 0:
	    new[ni++] = c;
	    break;
	    
	default:
	    new[ni++] = '\\';
	    new[ni++] = type[(unsigned char)c];
	}
    }
    new[ni++] = 0;

    return new;
}

/*
 * Allocates and returns an escaped version of str. This relaces quotes,
 * newlines, and other non-printable characters with %02X hex encoded
 * versions as required by html, gff, etc.
 *
 * 'escape' is a string of additional characters that must be escaped
 * for this string, in addition to obvious unprintables and percent.
 * It may be specified as NULL.
 *
 * Returns malloced string on success
 *         NULL on failure.
 */
char *escape_hex_string(char *str, char *escape) {
    size_t l = strlen(str);
    size_t new_l = l*1.1+10;
    char *new = malloc(new_l);
    size_t oi, ni;
    static int type[256];
    static int type_init = 0;
    int i;

    /* A once-only lookup table to speed up the loop below */
    if (!type_init) {
	for (i = 0; i < 256; i++) {
	    if (isprint(i) && i != '%') {
		/* directly printable */
		type[i] = 0;
	    } else {
		/* hex escape */
		type[i] = 1;
	    }
	}
	type_init = 1;
    }


    /* Per call modifications to the basic escape rules */
    for (i = 0; i < 256; i++) {
	type[i] &= 1;
    }

    if (escape) {
	while (*escape) {
	    type[(unsigned char)*escape] |= 2;
	    escape++;
	}
    }


    if (!new)
	return NULL;

    for (oi = ni = 0; oi < l; oi++) {
	char c = str[oi];

	/* Make enough room */
	if (ni + 4 >= new_l) {
	    new_l = new_l * 1.2 + 10;
	    if (NULL == (new = realloc(new, new_l)))
		return NULL;
	}

	if (type[(unsigned char)c]) {
	    sprintf(&new[ni], "%%%02X", c);
	    ni+=3;
	} else {
	    new[ni++] = c;
	}
    }
    new[ni++] = 0;

    return new;
}

/*
 * Reversal of the escape_hex_string above.
 *
 * Returns a copy of the escaped string on success
 *         NULL on failure.
 *
 * The pointer returned is owned by this function and is valid until the
 * next call (so it is not reentrant). DO NOT FREE the result.
 */
char *unescape_hex_string(char *str) {
    static char *ret = NULL;
    static size_t ret_sz = 0;
    static int hex[256];
    static int hex_init = 0;
    size_t l;
    char *out;


    if (!str)
	return NULL;
    

    /* Initialise lookup tables */
    if (!hex_init) {
	int i;
	memset(hex, 0, 256*sizeof(*hex));
	for (i = 0; i <= 9; i++) {
	    hex['0'+i] = i;
	}
	for (i = 0; i <= 5; i++) {
	    hex['a'+i] = 10+i;
	    hex['A'+i] = 10+i;
	}

	hex_init = 1;
    }


    /* Alloc memory */
    l = strlen(str);
    if (l >= ret_sz) {
	ret_sz = l+1;
	ret = realloc(ret, ret_sz);
	if (!ret) {
	    return NULL;
	    ret_sz = 0;
	}
    }


    /* Decode */
    out = ret;
    while (*str) {
	if (*str == '%') {
	    if (!str[1]) {
		fprintf(stderr,"Truncated %% code in unescape_hex_string()\n");
		return NULL;
	    }
	    *out++ = (hex[str[1]]<<4) | hex[str[2]];
	    str += 3;
	} else {
	    *out++ = *str++;
	}
    }
    *out++ = 0;

    return ret;
}
