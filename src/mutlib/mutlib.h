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


#ifndef _MUTLIB_H_
#define _MUTLIB_H_


#include <staden.h>        /* For Read* structure */


#ifdef __cplusplus
extern "C" {
#endif


/*---------------------*/
/* Generic Definitions */
/*---------------------*/

/** Strand type */
typedef enum
{
   MUTLIB_STRAND_FORWARD,
   MUTLIB_STRAND_REVERSE

}mutlib_strand_t;


/** Input type */
typedef enum
{
    MUTLIB_INPUT_REFERENCE,
    MUTLIB_INPUT

}mutlib_input_t;


/** Trace type */
typedef struct
{
   int             New;
   int             ClipL;
   int             ClipR;
   mutlib_strand_t Strand;
   Read*           Trace;

}mutlib_trace_t;


/** Sequence type */
typedef struct
{
   int             New;
   int             ClipL;
   int             ClipR;
   mutlib_strand_t Strand;
   char*           Sequence;

}mutlib_sequence_t;


/** TC/TG experiment file tag */
typedef struct
{
   char            Type[5];
   mutlib_strand_t Strand;
   int             Position[2];
   char*           Comment;
   int             Marked;

}mutlib_tag_t;


/** Algorithm result codes */
typedef enum
{
   MUTLIB_RESULT_SUCCESS,
   MUTLIB_RESULT_INVALID_INPUT,
   MUTLIB_RESULT_INSUFFICIENT_OVERLAP,
   MUTLIB_RESULT_INSUFFICIENT_DATA,
   MUTLIB_RESULT_ALIGNMENT_FAILURE,
   MUTLIB_RESULT_OUT_OF_MEMORY,
   MUTLIB_RESULT_UNEXPECTED_EXCEPTION

}mutlib_result_t;



/*------------------------*/
/* TraceAlign Definitions */
/*------------------------*/

/** Trace align object */
typedef struct
{
   /* Input objects owned by caller */
   mutlib_trace_t  Input;
   mutlib_trace_t  Reference[2];

   /* Output objects owned by tracealign */
   mutlib_trace_t  Alignment[2];
   mutlib_result_t ResultCode;
   char*           ResultString;

   /* Internal objects owned by tracealign */
   void*           Cache;
   int             Initialised;

}tracealign_t;



/*----------------*/
/* TraceAlign API */
/*----------------*/

void            TraceAlignInit( tracealign_t* ta );
void            TraceAlignSetInput( tracealign_t* ta, mutlib_strand_t d, Read* r, int ql, int qr );
void            TraceAlignSetReference( tracealign_t* ta, mutlib_strand_t d, Read* r, int ql, int qr );
mutlib_result_t TraceAlignExecute( tracealign_t* ta );
mutlib_result_t TraceAlignGetResultCode( tracealign_t* ta );
const char*     TraceAlignGetResultString( tracealign_t* ta );
Read*           TraceAlignGetAlignment( tracealign_t* ta, mutlib_input_t i, int* l, int* r );
void            TraceAlignDestroy( tracealign_t* ta );



/*-----------------------*/
/* TraceDiff Definitions */
/*-----------------------*/

/** Maximum number of parameters */
#define TRACEDIFF_PARAMETERS 7


/* Trace difference parameters */
typedef enum
{
   TRACEDIFF_PARAMETER_SENSITIVITY,          /* Number of SD's from mean over noise window */
   TRACEDIFF_PARAMETER_NOISE_THRESHOLD,      /* Threshold below which is noise as % of max */
   TRACEDIFF_PARAMETER_NOISE_WINDOW_LENGTH,  /* Window in bases over which background noise level is computed */
   TRACEDIFF_PARAMETER_PEAK_ALIGNMENT,       /* Upper/lower peak centre alignment threshold in bases */
   TRACEDIFF_PARAMETER_PEAK_WIDTH_MAXIMUM,   /* Threshold in bases above which peak is considered to be too wide */
   TRACEDIFF_PARAMETER_COMPLEMENT_TAGS,      /* Complement reverse strand tags? */
   TRACEDIFF_PARAMETER_YSCALE                /* Scale traces in Y before differencing? */

}tracediff_parameter_t;


/** Trace difference algorithms */
typedef enum
{
   /* New algorithms should be powers of 2 */
   TRACEDIFF_ALGORITHM_DEFAULT                 = 1,
   TRACEDIFF_ALGORITHM_DEFAULT_DIFFERENCE_ONLY = 2

}tracediff_algorithm_t;


/** Trace difference object */
typedef struct
{
   /* Input objects owned by caller */
   tracealign_t    Alignment;
   double          Parameter[TRACEDIFF_PARAMETERS];

   /* Output objects owned by tracediff */
   mutlib_tag_t*   Tag;
   int             TagCount;
   Read*           Difference;
   int             DifferenceLeft;
   int             DifferenceRight;
   mutlib_result_t ResultCode;
   char*           ResultString;

   /* Internal objects owned by tracediff */
   int             Initialised;

}tracediff_t;



/*---------------*/
/* TraceDiff API */
/*---------------*/

void            TraceDiffInit( tracediff_t* td );
void            TraceDiffSetParameter( tracediff_t* td, tracediff_parameter_t p, double v );
double          TraceDiffGetParameter( tracediff_t* td, tracediff_parameter_t p );
void            TraceDiffSetReference( tracediff_t* td, Read* w, mutlib_strand_t d, int ql, int qr );
void            TraceDiffSetInput( tracediff_t* td, Read* i, mutlib_strand_t d, int ql, int qr );
mutlib_result_t TraceDiffExecute( tracediff_t* td, tracediff_algorithm_t a );
void            TraceDiffDestroy( tracediff_t* td );
mutlib_result_t TraceDiffGetResultCode( tracediff_t* td );
const char*     TraceDiffGetResultString( tracediff_t* td );
Read*           TraceDiffGetDifference( tracediff_t* td, int* dl, int* dr );
int             TraceDiffGetTagCount( tracediff_t* td );
mutlib_tag_t*   TraceDiffGetTag( tracediff_t* td, int n );



/*---------------------*/
/* MutScan Definitions */
/*---------------------*/

/** Maximum number of parameters */
#define MUTSCAN_PARAMETERS  7


/** Mutation scanning parameters */
typedef enum
{
   MUTSCAN_PARAMETER_ALIGNFAIL_THRESHOLD,   /* Reports alignment failure after n mutations */
   MUTSCAN_PARAMETER_COMPLEMENT_TAGS,       /* Complement reverse strand tags? */
   MUTSCAN_PARAMETER_HETSNR_THRESHOLD,      /* SNR below which mutations are considered to be heterozygotes */
   MUTSCAN_PARAMETER_PEAKDROP_LOWER,        /* Lower peak height drop threshold in percent */
   MUTSCAN_PARAMETER_NOISE_THRESHOLD,       /* Noise floor threshold as percent of maximum peak height */
   MUTSCAN_PARAMETER_PEAKDROP_UPPER,        /* Upper peak height drop threshold in percent */
   MUTSCAN_PARAMETER_SEARCH_WINDOW          /* Search window to find corresponding peaks in bases */

}mutscan_parameter_t;



/** Mutscan object */
typedef struct
{
   /* Input objects owned by caller */
   mutlib_trace_t    InputTrace;
   mutlib_trace_t    ReferenceTrace[2];
   double            Parameter[MUTSCAN_PARAMETERS];

   /* Output objects owned by mutscan */
   mutlib_tag_t*     Tag;
   int               TagCount;
   mutlib_result_t   ResultCode;
   char*             ResultString;

   /* Internal objects owned by mutscan */
   int               Initialised;

}mutscan_t;



/*-------------*/
/* MutScan API */
/*-------------*/

void            MutScanInit( mutscan_t* ms );
void            MutScanSetParameter( mutscan_t* ms, mutscan_parameter_t p, double v );
double          MutScanGetParameter( mutscan_t* ms, mutscan_parameter_t p );
void            MutScanSetReference( mutscan_t* ms, mutlib_strand_t d, Read* r, int ql, int qr );
void            MutScanSetInput( mutscan_t* ms, mutlib_strand_t d, Read* i, int ql, int qr );
mutlib_result_t MutScanExecute( mutscan_t* ms );
void            MutScanDestroy( mutscan_t* ms );
mutlib_result_t MutScanGetResultCode( mutscan_t* ms );
const char*     MutScanGetResultString( mutscan_t* ms );
int             MutScanGetTagCount( mutscan_t* ms );
mutlib_tag_t*   MutScanGetTag( mutscan_t* ms, int n );




#ifdef __cplusplus
}
#endif


#endif

