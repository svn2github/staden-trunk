#ifndef _misc_h
#define _misc_h

#include "os.h"

#include <stdio.h>
#include <stdarg.h>  /* varargs needed for v*printf() prototypes */
#include <sys/types.h>

#include "xalloc.h"

#ifdef _MSC_VER
extern int getopt( int argc, char* const argv[], const char* optstring );
#endif

/*
 * This informs gcc that crash() doesn't return, so it doesn't need to
 * concern itself that code paths going via crash could mean some variables
 * being undefined and then issuing uninitialised variable warnings.
 * This particularly affected convert.
 */
#ifdef __GNUC__
#    define __NORETURN__ __attribute__ ((__noreturn__))
#else
#    define __NORETURN__
#endif

/*
 * Used for printf style argument checking. We can request a function such
 * as vTcl_SetResult does argument checking, avoiding bugs with using
 * %d and passing in a 64-bit record.
 */
#ifdef __GNUC__
#    define __PRINTF_FORMAT__(a,b) __attribute__ ((format (printf, a, b)))
#else
#    define __PRINTF_FORMAT__(a,b)
#endif

extern int is_directory(char * fn);
extern int is_file(char * fn);
extern int file_exists(char * fn);
extern int file_size(char * fn);
extern FILE *open_fofn(char *files);
extern char *read_fofn(FILE *fp);
extern void close_fofn(FILE *fp);
extern int fstrlen(char *f, int max_f);
extern void f2cstr(char *f, int max_f, char *c, int max_c);
extern void c2fstr(char *c, int max_c, char *f, int max_f);
extern char *mystrtok(char *s, char *ct);
extern char *myfind(char *file, char* searchpath, int (*found) (char *) );
extern void crash (char* format,...) __NORETURN__ __PRINTF_FORMAT__(1,2);
extern void str_tolower (char *s);
extern void str_toupper (char *s);
extern char *fn_tail (char *s);
extern void fn_tolower (char *s);
extern void fn_toupper (char *s);
extern void shell_call(char *command, char *output, int len);
extern char *date_str(void);
#ifdef NOSTRDUP
extern char *strdup(const char *s);
#endif
#ifdef NOSTRSTR
extern char *strstr(char *cs, char *ct);
#endif

#ifdef NOMEMMOVE
extern void *memmove(void *s1, const void *s2, size_t n);
#endif
extern int myusleep(unsigned int useconds);

extern void errout(char *fmt, ...) __PRINTF_FORMAT__(1,2);
extern void messout(char *fmt, ...) __PRINTF_FORMAT__(1,2);

/*
 * Useful macros
 */
#define findfile(F,S) myfind((F),(S),file_exists)
/*is_file fails for symbolic links*/
/*#define findfile(F,S) myfind((F),(S),is_file)*/

#if defined(min)
#undef min
#undef max
#endif

#if !defined(__cplusplus)
#define min(A,B) ( ( (A) < (B) ) ? (A) : (B) )
#define max(A,B) ( ( (A) > (B) ) ? (A) : (B) )
#define sgn(A) ( (A) ? ( ( (A) < 0 ) ? -1 : 1 ) : 0 )
#endif

#ifdef MIN
#undef MIN
#endif
#define MIN(A,B) ( ( (A) < (B) ) ? (A) : (B) )
#ifdef MAX
#undef MAX
#endif
#define MAX(A,B) ( ( (A) > (B) ) ? (A) : (B) )
#define SGN(A) ( (A) ? ( ( (A) < 0 ) ? -1 : 1 ) : 0 )
#define ABS(A) ( (A) < 0 ? -(A) : (A) )

/* Number of elements in array */
#define Number(A) ( sizeof(A) / sizeof((A)[0]) )

/*
 * Things taken from the new gap text_output.h. They'll be used globally
 * across all the programs in the end.
 */

/*
 * Usage: verror(priority, format, args...);
 * NB: don't pass more than 8K per call
 */
#define ERR_WARN 0
#define ERR_FATAL 1
void verror(int priority, const char *name, const char *fmt, ...) __PRINTF_FORMAT__(3,4);

/*
 * Usage: vmessage(format, args...);
 * NB: don't pass more than 8K per call
 */
void vmessage(const char *fmt, ...) __PRINTF_FORMAT__(1,2);

/*
 * Adds a new header to the text output window.
 */
void vfuncheader(const char *fmt, ...) __PRINTF_FORMAT__(1,2);

/*
 * As vfuncheader, but only outputting when necessary.
 */
void vfuncgroup(int group, const char *fmt, ...) __PRINTF_FORMAT__(2,3);


#ifdef NOSTRCASECMP

/* Avoid definition in TK/TCL */
#ifdef strcasecmp
#undef strcasecmp
#endif

int strcasecmp(const char *s3, const char *s2);
#endif


/*
 * Strnlen: like strlen(), but with a maximum size so we can do strlen on
 * potentially non-null terminated arrays.
 */
#if !defined(__USE_GNU)
size_t strnlen(const char *buf, size_t n);
#endif

#endif /*_misc_h*/
