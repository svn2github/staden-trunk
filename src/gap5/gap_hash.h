#ifndef _GAP_HASH_H_
#define _GAP_HASH_H_


void set_hash8_lookup(void);

int repeat_search (int mode, int min_match, int **seq1_match, int **seq2_match,
		   int **len_match, int max_mat, char *seq1, int seq1_len, 
		   int *num_f_matches, int *num_r_matches);

int repeat_search_depadded(int mode,         /* 1=f, 2=r, 3=b */
			   int min_match,    /* the minimum match length */
			   int **seq1_match, /* positions of matches in seq1 */
			   int **seq2_match, /* positions of matches in seq2 */
			   int **len_match,  /* length of matches */
			   int max_mat,	     /* maximum number of matches */
			   char *seq1,	     /* seq1 */
			   int seq1_len,     /* size of seq1 */
			   int *num_f_matches,
			   int *num_r_matches);

int cmpseq_ (
		   int   *job,		/* the task to perform */
	           char *sense,     	/* the orientation of seq2 */
		   int *min_match,	/* the minimum match length */
		   int *seq1_match,	/* positions of matches in seq1 */
		   int *seq2_match,	/* positions of matches in seq2 */
		   int *len_match,	/* length of matches */
		   int *max_matches,	/* maximum number of matches */
		   char *seq1,		/* seq1 */
		   char *seq2,		/* seq2 */
		   int *seq1_len, 	/* size of seq1 and its hash array */
		   int *seq2_len,	/* size of seq2 */
	           int_fl seq1_l,       /* fortran string length for seq1 */
	           int_fl seq2_l);      /* fortran string length for seq2 */


#endif /* _GAP_HASH_H_ */
