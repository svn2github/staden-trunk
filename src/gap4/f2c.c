#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include "os.h"
#include "misc.h"
#include "f2c.h"

#define ADD_BUF(n) \
do { \
  if (width && precision) { \
     buf += sprintf(buf,format,conv_len1,conv_len2,n); \
  } else if (width) { \
     buf += sprintf(buf,format,conv_len1,n); \
  } else if (precision) { \
     buf += sprintf(buf,format,conv_len2,n); \
  } else { \
     buf += sprintf(buf,format,n); \
  } \
} while (0)

#define MIN(A,B) ( ( (A) < (B) ) ? (A) : (B) )
#define MAX(A,B) ( ( (A) > (B) ) ? (A) : (B) )
#define SGN(A) ( (A) ? ( ( (A) < 0 ) ? -1 : 1 ) : 0 )
#define ABS(A) ( (A) < 0 ? -(A) : (A) )

/*
 * This is a fortran callable equivalent to sprintf. This avoids using the
 * fortran WRITE and FORMAT commands, hence freeing up the dependency on
 * needing that part of libI77. (Useful as a first step of removing the
 * fortran code and allowing f2c to work better.)
 *
 * It varies in a few ways to the real printf (IMPORTANT, especially %!):
 *
 * 1. Strings in fortran are not null terminated. The length is added
 * as an implicit argument to the end of the argument list. This makes things
 * VERY hard (if not impossible) to deal with when using varargs. Therefore
 * the end of the format string MUST be specified using "%!".
 *
 * 2. Following on from point 1, strings specified using %s MUST have their
 * size specified either via a literal %.10s or a variable %.*s. The
 * FORTRAN LEN() command is useful here. Eg:
 * CALL WRITEF('msg %.*s%!', LEN(S), S)
 *
 * For now, no mallocing is done so the buffer is assumed to be big enough.
 * Returns 0 for success
 *        -1 for length;
 */
int swritfv(char *sbuf, char *fmt, va_list ap)
{
    char *cp, c;
    long l;
    int i;
    float f;
    char *buf;
    int abort = 0;
    char format[100];

#if defined(NEED_VA_COPY)
    va_list ap_local;
    va_copy(ap_local, ap);
    #define ap ap_local
#endif

    buf = sbuf;

    for(cp = fmt; !abort && *cp; cp++) {
	int width, precision;

	switch(*cp) {
	/* A format specifier */
	case '%': {
	    char *endp, *formatp = format;
	    long conv_len1=0, conv_len2=0, conv_len=0;
	    signed int arg_size;

	    formatp = cp;
	    width = 0;
	    precision = 0;

	    /* Firstly, strip the modifier flags (+-#0 and [space]) */
	    for(; c=*++cp;) {
		if ('!' == c) {
		    abort = 1;
		    break;
		}
		if (!('#' == c || '-' == c || '+' == c || ' ' == c))
		    break;
	    }
	    if (abort)
		break;

	    /* Width specifier */
	    l = strtol(cp, &endp, 10);
	    if (endp != cp) {
		cp = endp;
		conv_len = conv_len1 = l;
	    } else if (*cp == '*') {
		width = 1;
		conv_len = conv_len1 = *((int *)va_arg(ap, int *));
		cp++;
	    }

	    /* Precision specifier */
	    if ('.' == *cp) {
		cp++;
		conv_len2 = strtol(cp, &endp, 10);
		if (endp != cp) {
		    cp = endp;
		} else if (*cp == '*') {
		    precision = 1;
		    conv_len2 = *((int *)va_arg(ap, int *));
		    cp++;
		}
		conv_len = MAX(conv_len1, conv_len2);
	    }

	    /* Short/long identifier */
	    if ('h' == *cp) {
		arg_size = -1; /* short */
		cp++;
	    } else if ('l' == *cp) {
		arg_size = 1; /* long */
		cp++;
	    } else {
		arg_size = 0; /* int */
	    }

	    memcpy(format, formatp, cp-formatp+1);
	    format[cp-formatp+1] = 0;

	    /* The actual type */
	    switch (*cp) {
	    case '%':
		ADD_BUF('%');
		break;

	    case 'd':
	    case 'i':
	    case 'u':
	    case 'a':
	    case 'x':
	    case 'X':
		/* Remember: char and short are sent as int on the stack */
		if (arg_size == -1)
		    l = *((long *)va_arg(ap, int *));
		else if (arg_size == 1)
		    l = *(va_arg(ap, long *)); 
		else 
		    l = *((int *)va_arg(ap, int *));

		ADD_BUF(l);
		break;

	    case 'c':
		i = *(va_arg(ap, int *));
		ADD_BUF(i);
		break;

	    case 'f':
		f = *(va_arg(ap, float *));
		ADD_BUF(f);
		break;

	    case 'e':
	    case 'E':
	    case 'g':
	    case 'G':
		f = *(va_arg(ap, float *));
		ADD_BUF(f);
		break;

	    case 'p':
		l = (long)va_arg(ap, void *);
		ADD_BUF((void *)l);
		break;

	    case 'n':
		/* produces no output */
		break;

	    case 's': {
		char *s = (char *)va_arg(ap, char *);
		ADD_BUF(s);
		break;
	    }

#ifndef NDEBUG
	    default:
		/* wchar_t types of 'C' and 'S' aren't supported */
		printf("Unknown arg is %c\n", *cp);
#endif
	    }
	}

	case '\0':
	    break;

	default:
	    *buf++ = *cp;
	}
    }

    *buf = 0;

    va_end(ap);
    return 0;
}

/*
 * A series of identical interfaces to swritfv above. This is what the fortran
 * code calls. Care has been taken to never call the same function with either
 * a different number of arguments or different types of arguments as this is
 * not supported in fortran (despite the fact that we know it works just
 * fine). This isn't strictly necessary, but it removes a lot of compiler
 * warnings.
 */
int swrt0_(char *sbuf, char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    return swritfv(sbuf, fmt, args);
}

int swrt1_(char *sbuf, char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    return swritfv(sbuf, fmt, args);
}

int swrt2_(char *sbuf, char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    return swritfv(sbuf, fmt, args);
}

int swrt2b_(char *sbuf, char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    return swritfv(sbuf, fmt, args);
}

int swrt3_(char *sbuf, char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    return swritfv(sbuf, fmt, args);
}

int swrt3b_(char *sbuf, char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    return swritfv(sbuf, fmt, args);
}

int swrt4_(char *sbuf, char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    return swritfv(sbuf, fmt, args);
}

int swrt5_(char *sbuf, char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    return swritfv(sbuf, fmt, args);
}

/*
 * Fortran version from the old JFROMC function.
 */
int jfromc_(char *data, integer *length_p, ftnlen l) {
    char tmp[1024];
    memcpy(tmp, data, *length_p);
    tmp[*length_p] = 0;
    return atoi(tmp);
}

/*
 * -----------------------------------------------------------------------------
 * Functions taken from F2C library distribution and tidied up a bit.
 */

/*
 * Compares two fortran strings to see if they are equal.
 */
integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
{
    unsigned char *a, *aend, *b, *bend;
    a = (unsigned char *)a0;
    b = (unsigned char *)b0;
    aend = a + la;
    bend = b + lb;

    if (la <= lb) {
	while(a < aend)
	    if(*a != *b)
		return( *a - *b );
	    else
		++a, ++b;
	
	while(b < bend)
	    if(*b != ' ')
		return( ' ' - *b );
	    else
		++b;

    } else {

	while(b < bend)
	    if(*a == *b)
		++a, ++b;
	    else
		return( *a - *b );

	while(a < aend)
	    if(*a != ' ')
		return(*a - ' ');
	    else
		++a;
    }

    return(0);
}


/*
 * Copies string 'b' to string 'a' with the Fortran 77 limitation
 * that the two strings cannot overlap.
 */
int s_copy(char *a, char *b, ftnlen la, ftnlen lb)
{
    char *aend, *bend;

    aend = a + la;

    if(la <= lb) {
	while(a < aend)
	    *a++ = *b++;
    } else {
	bend = b + lb;
	while(b < bend)
	    *a++ = *b++;

	while(a < aend)
	    *a++ = ' ';
    }

    return 0;
}

/*
 * The length of a string. Just returns the implicit length passed in. How pointless!
 */
integer i_len(char *s, ftnlen n) { return(n); }

/*
 * The fortran SIGN intrinsic. Returns (+/-)ABS(A) depending on if B>=0
 */
integer i_sign(integer *a, integer *b)
{
    integer x;
    x = (*a >= 0 ? *a : - *a);
    return( *b >= 0 ? x : -x);
}
