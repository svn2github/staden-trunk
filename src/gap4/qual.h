#ifndef _QUAL_H
#define _QUAL_H

#include "os.h"

/*
 * ----------------------------------------------------------------------------
 * Macros
 * ----------------------------------------------------------------------------
 */

/* Job names and structures for the `info' function */
#define GET_SEQ		0
#define DEL_SEQ		1
#define GET_CONTIG_INFO	2
#define DEL_CONTIG_INFO	3 /* not implemented */
#define GET_GEL_INFO	4
#define DEL_GEL_INFO	5 /* not implemented */
#define GET_GEL_LEN	6
#define SEQ_INS		7
#define SEQ_DEL		8
#define CONS_INS	9
#define CONS_DEL	10
#define IF_FLUSH	11

/* quality codes */
#define R_GOOD_GOOD_EQ	'a'
#define R_GOOD_BAD	'b'
#define R_BAD_GOOD	'c'
#define R_GOOD_NONE	'd'
#define R_NONE_GOOD	'e'
#define R_BAD_BAD 	'f'
#define R_BAD_NONE	'g'
#define R_NONE_BAD	'h'
#define R_GOOD_GOOD_NE	'i'
#define R_NONE_NONE	'j'

/* consensus calculation codes */
#define CON_SUM		0
#define CON_WDET	1

#define CONSENSUS_MODE_FREQ	  0
#define CONSENSUS_MODE_WEIGHTED	  1 /* as FREQ, determined by qual_cut >= 0 */
#define CONSENSUS_MODE_CONFIDENCE 2

/* quality defaults - an invalid quality value so we never get confused */
#define QUAL_DEFAULT   -111

/*
 * ----------------------------------------------------------------------------
 * Structures and typedefs
 * ----------------------------------------------------------------------------
 */

typedef struct _gel_seq_t {
    int   gel;
    int  gel_length;
    int  gel_start;
    int  gel_end;
    char *gel_seq;
    int1 *gel_conf;
    int2 *gel_opos;
} gel_seq_t;

typedef struct _gel_info_t {
    int gel;
    int length;
    int complemented;
    int position;
    int as_double;
    int next_right;
    int start;
    int unclipped_len;
    int template;
} gel_info_t;

typedef struct _contig_info_t {
    int contig;
    int length;
    int leftgel;
} contig_info_t;

typedef struct _seq_ins_t {
    int gel;
    int position;
    int length;
    char *bases;
} seq_ins_t;

typedef struct _seq_del_t {
    int gel;
    int position;
    int length;
} seq_del_t;

typedef struct _cons_ins_t {
    int contig;
    int position;
    int length;
    char *bases;
} cons_ins_t;

typedef struct _cons_del_t {
    int contig;
    int position;
    int length;
} cons_del_t;

typedef union _info_arg_t {
    gel_seq_t	  gel_seq;
    gel_info_t	  gel_info;
    contig_info_t contig_info;
    seq_ins_t	  seq_ins;
    seq_del_t	  seq_del;
    cons_ins_t	  cons_ins;
    cons_del_t	  cons_del;
} info_arg_t;


/*
 * ----------------------------------------------------------------------------
 * Function prototypes
 * ----------------------------------------------------------------------------
 */

int calc_consensus(int   contig,
		   int   start,
		   int   end,
		   int   mode,
		   char *con,
		   char *con2,
		   float *qual,
		   float *qual2,
		   float cons_cutoff,
		   int   qual_cutoff,
		   int (*info_func)(int        job,
				    void       *mydata,
				    info_arg_t *theirdata),
		   void *info_data);

int calc_quality(int   contig,
		 int   start,
		 int   end,
		 char *qual,
		 float cons_cutoff,
		 int   qual_cutoff,
		 int (*info_func)(int        job,
				  void       *mydata,
				  info_arg_t *theirdata),
		 void *info_data);

int next_hole(int contig,
	      int position,
	      int rreg,
	      float cons_cutoff,
	      int   qual_cutoff,
	      char **reason,
	      int *len,
	      int (*info_func)(int        job,
			       void       *mydata,
			       info_arg_t *theirdata),
	      void *info_data);

int database_info(int job, void *mydata, info_arg_t *theirdata);

int set_qual_cutoff(int new);

int query_qual_cutoff(void);

int calc_discrepancies(int   contig,
		       int   start,
		       int   end,
		       float *qual1,
		       float *qual2,
		       float cons_cutoff,
		       int   qual_cutoff,
		       int (*info_func)(int        job,
					void       *mydata,
					info_arg_t *theirdata),
		       void *info_data);

#endif /* _QUAL_H */
