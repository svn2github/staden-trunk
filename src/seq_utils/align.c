#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "os.h"
#include "align.h"
#include "align_lib_old.h"

FastInt W128[128][128] = {{0},{0}};
char base_val[128] = {0};

static FastInt (*align_arr[])() = {align_ss,  
				   align_sv,
				   /*balign_sv*/NULL,
				   align_ss2};
static void (*display_arr[])()  = {display_ss,
				   display_sv,
				   display_sv,
				   display_ss2};
static void (*expand_arr[])()   = {expand,
				   expand_6,
				   expand_6,
				   expand};

/*
 * Initialise 128x128 weight matrix from our score matrix
 */
void init_align_mat(int **matrix,
		    char *order, int unknown,
		    FastInt W128[][128]) {
    int i, j, i_end;
    unsigned char ci, cj;

    for (i = 0; i < 128; i++) {
        for (j = 0; j < 128; j++) {
	    W128[i][j] = unknown;
        }
    }

    i_end = strlen(order);
    for (i = 0; i < i_end; i++) {
	ci = order[i];
	for (j = 0; j < i_end; j++) {
	    cj = order[j];
	    W128[        ci ][        cj ] = matrix[i][j];
	    W128[tolower(ci)][        cj ] = matrix[i][j];
	    W128[        ci ][tolower(cj)] = matrix[i][j];
	    W128[tolower(ci)][tolower(cj)] = matrix[i][j];
	}
    }

    for (i=0; i<128; i++)
	base_val[i] = 5;
    
    base_val['A'] = 0;
    base_val['a'] = 0;
    base_val['C'] = 1;
    base_val['c'] = 1;
    base_val['G'] = 2;
    base_val['g'] = 2;
    base_val['T'] = 3;
    base_val['t'] = 3;
    base_val['U'] = 3;
    base_val['u'] = 3;
    base_val['*'] = 4;
    /* base_val['-'] = 5; the default */
}

void init_W128(int **matrix,
	       char *order, int unknown) {
    init_align_mat(matrix, order, unknown, W128);
}

/*
 * Generic alignment interface:
 *
 * Job 0 = Sequence against sequence global aligment
 * Job 1 = Sequence against vector   global aligment
 * Job 2 = vector   against vector   global aligment (unimplemented)
 * Job 3 = Job 0 except with end gaps no penalised (to replace 0 eventually)
 *
 * Passing low_band == 0 && high_band == 0 implies using the full, rather than
 * band, algorithm.
 * Currently the band algorithm does not work so well always use the full.
 *
 * Passing rseq1 and rseq2 as non NULL implies storing the new aligned
 * sequences in rseq1 and rseq2.
 *
 * Protein alignment only works for sequence against sequence.
 *
 * Return values:
 *   -1 for error
 *   >0 for success (value returned is score)
 */
int calignm(void *seq1, void *seq2, int len1, int len2,
	    void *rseq1, void *rseq2, int *rlen1, int *rlen2,
	    int low_band, int high_band, int gap_open, int gap_extend,
	    int job, int is_protein, align_int *res, FastInt W[][128]) {

    /* int full = low_band == 0 && high_band == 0; */
    int retval = -1;
    FastInt *results;
    int ajob = job & ALIGN_J_MASK;

    if (ajob < 0 || ajob > 3) {
	verror(ERR_FATAL, "align", "unknown job %d", ajob);
	return -1;
    }

    if (res == NULL) {
	if (NULL == (results = (FastInt *)malloc((len1+len2) *
						 sizeof(FastInt)))) {
	    verror(ERR_FATAL, "align", "not enough memory");
	    return -1;
	}
    } else {
	results = (FastInt *)res;
    }

    retval = align_arr[ajob](seq1, seq2, (FastInt)len1, (FastInt)len2,
			     (FastInt)low_band, (FastInt)high_band, W,
			     (FastInt)gap_open, (FastInt)gap_extend,
			     (FastInt *)results,
			     /* End gaps */
			     (job & ALIGN_GAP_S1) ? 1 : 0,
			     (job & ALIGN_GAP_S2) ? 1 : 0,
			     (job & ALIGN_GAP_E1) ? 1 : 0,
			     (job & ALIGN_GAP_E2) ? 1 : 0
			     );

    /* non vectors */
    if (rseq1 && rseq2) {
	if (retval != -1) {
	    expand_arr[ajob](seq1, seq2, len1, len2,
			     rseq1, rseq2, rlen1, rlen2,
			     results, job & ALIGN_J_PADS);
	}
    }

    if (res == NULL)
	free (results);

    return retval;
}

/*
 * As above, but hardcoded to use W128
 */
int calign(void *seq1, void *seq2, int len1, int len2,
	   void *rseq1, void *rseq2, int *rlen1, int *rlen2,
	   int low_band, int high_band, int gap_open, int gap_extend,
	   int job, int is_protein, align_int *res) {
    return calignm(seq1, seq2, len1, len2, rseq1, rseq2, rlen1, rlen2,
		   low_band, high_band, gap_open, gap_extend,
		   job, is_protein, res, W128);
}

/*
 * Generic alignment display interface
 *
 * Job 0 = Sequence against sequence display
 * Job 1 = Sequence against vector   display
 * Job 2 = vector   against vector   display (unimplemented)
 * Job 3 = Job 0 except with end gaps no penalised (to replace 0 eventually)
 *
 */
void cdisplay(void *seq1, void *seq2, int len1, int len2, int job,
	      align_int *results, int pos1, int pos2) {
    display_arr[job & ALIGN_J_MASK](seq1, seq2, len1, len2,
				    results, pos1, pos2);
}

/*
 * Generic alignment expansion interface
 *
 * Job 0 = Sequence against sequence expand
 * Job 1 = Sequence against vector   expand
 * Job 2 = vector   against vector   expand (unimplemented)
 *
 * seq1 and seq2 are the original sequences. Expanded sequences will be
 * returned in rseq1 and rseq2 (it is expected that the caller allocate
 * memory for these buffers) and the lengths in rlen1 and rlen2.
 *
 */
void cexpand(void *seq1, void *seq2, int len1, int len2,
	     void *rseq1, void *rseq2, int *rlen1, int *rlen2,
	     int job, align_int *results) {
    expand_arr[job & ALIGN_J_MASK](seq1, seq2, len1, len2,
		    rseq1, rseq2, rlen1, rlen2,
		    results, job & ALIGN_J_PADS);
}


/*
 * FORTRAN interfaces
 */
int_f align_(void *seq1, int_f *len1, void *seq2, int_f *len2,
	    void *rseq1, int_f *rlen1, void *rseq2, int_f *rlen2,
	    int_f *low_band, int_f *high_band,
	    int_f *gap_open, int_f *gap_extend,
	    int_f *job, int_f *is_protein,
	    int_fl seq1_l, int_fl seq2_l, int_fl rseq1_l, int_fl rseq2_l) {
    return calign(seq1, seq2, *len1, *len2, rseq1, rseq2, rlen1, rlen2,
		  *low_band, *high_band, *gap_open, *gap_extend,
		  *job, *is_protein, NULL);
}
