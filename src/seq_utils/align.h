#ifndef _ALIGN_H
#define _ALIGN_H

#include "os.h"
#include "misc.h"

/*
#ifdef __alpha
typedef long FastInt;
#else
typedef int FastInt;
#endif
*/
typedef int FastInt;
typedef FastInt align_int;

#define MAX_SCORE_MATRIX 30

/* Basic job number */
#define ALIGN_J_MASK    0x0F
#define ALIGN_J_SS	0
#define ALIGN_J_SV	1
#define ALIGN_J_BSV	2
#define ALIGN_J_SSH	3

/* How to expand the sequence up */
#define ALIGN_J_PADS	0x10

/* Which ends to penalise gaps for. Only applies to ALIGN_J_SSH */
#define ALIGN_GAP_S1	0x020
#define ALIGN_GAP_E1	0x040
#define ALIGN_GAP_S2	0x080
#define ALIGN_GAP_E2	0x100

/*
 * Initialise 128x128 weight matrix from our score matrix
 */
extern void init_W128(int **matrix,
		      char *order, int unknown);

extern int calign(void *seq1, void *seq2, int len1, int len2,
		  void *rseq1, void *rseq2, int *rlen1, int *rlen2,
		  int low_band, int high_band, int gap_open, int gap_extend,
		  int job, int is_protein, align_int *res);

extern void cdisplay(void *seq1, void *seq2, int len1, int len2, int job,
		     align_int *results, int pos1, int pos2);

extern void cexpand(void *seq1, void *seq2, int len1, int len2,
		    void *rseq1, void *rseq2, int *rlen1, int *rlen2,
		    int job, align_int *results);

extern int_f align_(void *seq1, int_f *len1, void *seq2, int_f *len2,
		    void *rseq1, int_f *rlen1, void *rseq2, int_f *rlen2,
		    int_f *low_band, int_f *high_band,
		    int_f *gap_open, int_f *gap_extend,
		    int_f *job, int_f *is_protein, int_fl seq1_l,
		    int_fl seq2_l, int_fl rseq1_l, int_fl rseq2_l);

extern FastInt align_ss(char *A, char *B, FastInt M, FastInt N,
			FastInt low, FastInt up, FastInt W[][128],
			FastInt G, FastInt H, FastInt *S,
			FastInt s1, FastInt s2, FastInt e1, FastInt e2);

extern FastInt align_ss2(char *A, char *B, FastInt M, FastInt N,
			FastInt low, FastInt up, FastInt W[][128],
			FastInt G, FastInt H, FastInt *S,
			FastInt s1, FastInt s2, FastInt e1, FastInt e2);

extern FastInt align_sv(char *A, FastInt (*B)[6], FastInt M,
			FastInt N, FastInt low, FastInt up,
			FastInt W[][128], FastInt G, FastInt H, FastInt *S,
			FastInt s1, FastInt s2, FastInt e1, FastInt e2);

extern FastInt balign_sv(char *A, FastInt (*B)[6], FastInt M,
			FastInt N, FastInt low, FastInt up,
			FastInt W[][128], FastInt G, FastInt H, FastInt *S,
			FastInt s1, FastInt s2, FastInt e1, FastInt e2);

extern void display_ss(char *A, char *B, FastInt M, FastInt N,
		       FastInt *S, FastInt AP, FastInt BP);

extern void display_ss2(char *A, char *B, FastInt M, FastInt N,
		       FastInt *S, FastInt AP, FastInt BP);

extern void display_sv(char *A, FastInt (*B)[6], FastInt M, FastInt N,
		       FastInt *S, FastInt AP, FastInt BP);

#include "align_lib_old.h"

#endif /* _ALIGN_H */
