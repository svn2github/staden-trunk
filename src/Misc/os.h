/*
 * File: os.h
 *
 * Author: 
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: operating system specific type definitions
 *
 */

#ifndef _OS_H_
#define _OS_H_

#include <limits.h>

/*
 *-----------------------------------------------------------------------------
 * Machine specifics
 *-----------------------------------------------------------------------------
 */
/*
 * SunOS 4.x
 * Even though we use the ANSI gcc, we make use the the standard SunOS 4.x
 * libraries and include files, which are non-ansi
 */
#if defined(__sun__) && !defined(__svr4__)
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#define FOPEN_MAX 64
#define FILENAME_MAX 1024
#define SP_BIG_ENDIAN

/* Missing functions, defined in strings.c */
#define NOMEMMOVE
#define NOSTRERROR
#define BUGGY_SSCANF
#endif

/* 
 * SunOS 5.x - gcc or Sun's cc 
 */ 
#if (defined(__sun__) || defined(__sun)) && (defined(__svr4__) || defined(__SVR4))
#  if defined(__sparc__) || defined(__sparc) 
#    define SP_BIG_ENDIAN 
#  else 
#    define SP_LITTLE_ENDIAN 
#  endif 
#  define IMAGEDISPLAY 
#  define NOSTRDUP 
#endif 


/*
 * Microsoft Visual C++
 * Windows
 */
#if defined(_MSC_VER)
#define popen _popen
#define pclose _pclose
typedef int mode_t;
#define ftruncate(fd,len) _chsize(fd,len)
#define sysconf(x) 512
#define SP_LITTLE_ENDIAN
#define NOPIPE
#define NOLOCKF
#define NOSTRCASECMP
#define NO_STRPTIME
#endif


/*
 * Linux on Intel platforms
 */
#if defined(__linux__)
#  if defined(SP_BIG_ENDIAN)
#    undef SP_BIG_ENDIAN
#  endif
#  define SP_LITTLE_ENDIAN
#endif

/*
 * Linux on AMD64 also needs to use va_copy()
 */
#if defined(__linux__) && defined(__amd64__)
#  define NEED_VA_COPY
#endif

/*
 * DEC Alpha's running Digital UNIX
 */
#if defined(__alpha)
#  if defined(SP_BIG_ENDIAN)
#    undef SP_BIG_ENDIAN
#  endif
#  if !defined(SP_LITTLE_ENDIAN)
#    define SP_LITTLE_ENDIAN
#  endif
#endif

/*
 * Silicon Graphics - Irix
 */
#if defined(__sgi)
#define SP_BIG_ENDIAN
#define NOSTRDUP
#define NO_STRPTIME
#endif

/*
 * Macs (<= OS 9) - yuk!
 */
#if defined(MAC)
#define SP_BIG_ENDIAN
#define NOSTRDUP
#endif

#if defined(__APPLE__)
#define SP_BIG_ENDIAN
#define NO_STRPTIME
#define NOLOCKF
#endif

/*
 *-----------------------------------------------------------------------------
 * Typedefs for data sizes. Note there's umpteen versions of typedefs here
 * due to old code being supported. The ones that should be used everywhere
 * are {u,}int[124].
 *-----------------------------------------------------------------------------
 */

/*
 * One byte integers
 */ 
/*
typedef unsigned char	uint1;
typedef signed char	uint1;
*/
typedef unsigned char	int1;

/*
 * Two byte integers
 */
typedef signed short	int2;
typedef unsigned short	uint2;

/*
 * Four byte integers
 */
typedef signed int	int4;
typedef unsigned int	uint4;


/*
 * Backwards compatibility
 */
typedef signed char	int_1;
typedef unsigned char	uint_1;
typedef signed short	int_2;
typedef unsigned short	uint_2;
typedef signed int	int_4;
typedef unsigned int	uint_4;


/*
 *-----------------------------------------------------------------------------
 * The FORTRAN interface.
 *-----------------------------------------------------------------------------
 */

typedef int4 f_int;
typedef int4 f_implicit;
typedef void f_proc_ret;	/* procedure return value */

/* James Bonfield compatability mode */
typedef int4 int_f;		/* f_int */
typedef int4 int_fl;		/* f_implicit */

#define f_proc_return() return /* (f_proc_ret) 0 */

/*
 * Use when calling/defining a Fortran function from C.
 */
#ifdef VMS
#    define FORT(symbol) (symbol)
#else
#    define FORT(symbol) (_symbol)
#endif


/*
 *-----------------------------------------------------------------------------
 * Some handy definitions.
 *-----------------------------------------------------------------------------
 */

#define MAXINT4 (INT_MAX)
#define MAXINT2 (SHRT_MAX)


/*
 *=============================================================================
 * Anything below here should not be changed.
 *=============================================================================
 */

#define False 0
#define True 1


#ifdef OLD_SWAP
/* See below for reasons not to revert back. We ought to remove these. */

/* copy INT4 from src to dst byteswapping on the way */
#define swap_int4(src, dst) \
    do {\
	((char *)&(dst))[0] = ((char *) &(src))[3];\
	((char *)&(dst))[1] = ((char *) &(src))[2];\
        ((char *)&(dst))[2] = ((char *) &(src))[1];\
        ((char *)&(dst))[3] = ((char *) &(src))[0];\
    } while (0)

/* copy INT2 from src to dst byteswapping on the way */
#define swap_int2(src, dst) \
    do {\
        ((char *) &(dst))[0] = ((char *) &(src))[1];\
        ((char *) &(dst))[1] = ((char *) &(src))[0];\
    } while (0)

#else

/*
 * Our new swap runs at the same speed on Ultrix, but substantially faster
 * (300% for swap_int4, ~50% for swap_int2) on an Alpha (due to the lack of
 * decent 'char' support).
 *
 * They also have the ability to swap in situ (src == dst). Newer code now
 * relies on this so don't change back!
 */
#define swap_int4(src, dst) \
    dst = ((src & 0x000000ff) << 24) + \
          ((src & 0x0000ff00) <<  8) + \
          ((src & 0x00ff0000) >>  8) + \
          ((src & 0xff000000) >> 24)

#define swap_int2(src, dst) \
    dst = ((src & 0x00ff) << 8) + \
          ((src & 0xff00) >> 8)
#endif


/*
 * Slightly updated swap_int? routines that return results rather than
 * swapping from source to destination.
 */
#define iswap_int4(x) \
    (((x & 0x000000ff) << 24) + \
     ((x & 0x0000ff00) <<  8) + \
     ((x & 0x00ff0000) >>  8) + \
     ((x & 0xff000000) >> 24))

#define iswap_int2(x) \
    (((x & 0x00ff) << 8) + \
     ((x & 0xff00) >> 8))

/*
 * Macros to specify that data read in is of a particular endianness.
 * The macros here swap to the appropriate order for the particular machine
 * running the macro and return the new answer. These may also be used when
 * writing to a file to specify that we wish to write in (eg) big endian
 * format.
 *
 * This leads to efficient code as most of the time these macros are
 * trivial.
 */
#ifdef SP_BIG_ENDIAN
#define be_int4(x) (x)
#define be_int2(x) (x)
#define be_int1(x) (x)

#define le_int4(x) iswap_int4((x))
#define le_int2(x) iswap_int2((x))
#define le_int1(x) (x)
#endif

#ifdef SP_LITTLE_ENDIAN
#define be_int4(x) iswap_int4((x))
#define be_int2(x) iswap_int2((x))
#define be_int1(x) (x)

#define le_int4(x) (x)
#define le_int2(x) (x)
#define le_int1(x) (x)
#endif

#ifndef SP_BIG_ENDIAN
#ifndef SP_LITTLE_ENDIAN
#error Must define SP_BIG_ENDIAN or SP_LITTLE_ENDIAN in Makefile
#endif
#endif

#ifdef SP_BIG_ENDIAN
#ifdef SP_LITTLE_ENDIAN
#error Must only define one of SP_BIG_ENDIAN and SP_LITTLE_ENDIAN in Makefile
#endif
#endif

#endif /*_OS_H_*/
