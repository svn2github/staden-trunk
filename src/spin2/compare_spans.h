#ifndef _COMPARE_SPANS_H_
#define _COMPARE_SPANS_H_
int cmpspn (
	    char *sense,	/* the orientation for seq2 */
	    int *min_mat,	/* the minimum match length */
	    int **seq1_match,	/* positions of matches in seq1 */
	    int **seq2_match,	/* positions of matches in seq2 */
	    int **match_score,	/* match scores */
	    int *max_match,	/* maximum number of matches */
	    int *window_l,	/* length of window */
	    char *seq1,		/* seq1 */
	    char *seq2,		/* seq2 */
	    int *seq1_l, 	/* size of seq1 */
	    int *seq2_l,	/* size of seq2 */
	    int seq1_lreg,      /* left cutoff of seq1 */
	    int seq1_rreg,      /* right cutoff of seq1 */
	    int seq2_lreg,      /* left cutoff of seq2 */
	    int seq2_rreg,      /* right cutoff of seq2 */
	    int same_seq        /* whether seq1 & seq2 are identical */
	    );
#endif
