/*
 * File: g-misc.h
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: miscellaneous macro definitions
 *
 * Created: prior to 18 September 1992
 * Updated:
 *
 */

#ifndef _G_MISC_H_
#define _G_MISC_H_

#define G_Number(A) ( sizeof(A) / sizeof(A[0]) )

/* for debugging */
#define GROK fprintf(stderr,"%s: %d\n",__FILE__,__LINE__)
#define YUK do {GROK;exit(1);} while(0)


#define G_Assert(B) \
    do {if (!(B)) fprintf(stderr,"Assertion failed: file %s line %d\n",__FILE__,__LINE__);} while(0)


#endif /*_G_MISC_H_*/
