#ifndef _SPLICE_SEARCH_H_
#define _SPLICE_SEARCH_H_

typedef struct _Wtmatch {	/* a match has a position, score and seq */
    int    pos;
    double score;
    char   *seq;
} Wtmatch;

typedef struct _WtmatrixRes {
    Wtmatch **match;           /* the matches: position, score, seq */
    int    number_of_res;      /* number of matches */
    int    length;	       /* length of the motif */
    int    mark_pos;	       /* offset to mark position. ie mark position+mark_pos */
    double min;                /* minimum cutoff score */
    double max;		       /* maximum possible score */
} WtmatrixRes;

typedef struct _SpliceResults {
    WtmatrixRes *ied_f1;	/* intron/exon (donor) frame 1 */
    WtmatrixRes *ied_f2;	/* intron/exon (donor) frame 2 */
    WtmatrixRes *ied_f3;	/* intron/exon (donor) frame 3 */

    WtmatrixRes *eia_f1;	/* exon/intron (acceptor) frame 1 */
    WtmatrixRes *eia_f2;	/* exon/intron (acceptor) frame 2 */
    WtmatrixRes *eia_f3;	/* exon/intron (acceptor) frame 3 */
} SpliceResults;

typedef struct splice_res_ {
    WtmatrixRes *ied;
    WtmatrixRes *eia;
    int frame;
} splice_res;

void free_WtmatrixRes ( WtmatrixRes *r );
int splice_search (char seq[], int seq_length, int user_start, int user_end,
		   char *filename_ied, char *filename_eia,
		   SpliceResults *splice_result);
int weight_search (char seq[], int seq_length, int user_start, int user_end,
		   char *filname, WtmatrixRes **results);
#endif
