#ifndef _SIP_HASH_H_
#define _SIP_HASH_H_

int sip_hash (char *seq1,
	      char *seq2,
	      int seq1_lreg, 
	      int seq1_rreg, 
	      int seq2_lreg, 
	      int seq2_rreg, 
	      int max_matches,
	      int min_match,
	      int word_len,
	      int sequence_type,
	      int same,
	      int **seq1_match,
	      int **seq2_match,
	      int **len_match,
	      int *n_matches,
	      void (*pr_func)(void *data, int pw, int y),
	      void *data);

    
/*	routine to remove duplicates from a list of repeats * * It
	also removes the self match.
	Input: a list of *n_match match positions and match
	lengths.  * Output: a list in which all duplicates and the selfmatch
	are removed.  * *n_match is set to the new number of matches, or -1
	for error.  */
void sip_remdup ( int **seq1_match, int **seq2_match, int **len_match, int
	     *n_matches );

int quick_scan(
	       char *seq1,		/* seq1 */
	       char *seq2,		/* seq2 */
	       int seq1_lreg,           /* start of seq1 */
	       int seq1_rreg,           /* end of seq1 */
	       int seq2_lreg,           /* start of seq2 */
	       int seq2_rreg,           /* end of seq2 */
	       int sequence_type,
	       int max_matches,	/* maximum number of matches */
	       int same_seq,		/* if set then its an internal comparison */
	       int win_length,
	       int min_match,	        /* the minimum match length */
	       int word_len,
	       double min_sd,              /* minimum std dev */
	       int save_results,
	       int **seq1_match,	/* positions of matches in seq1 */
	       int **seq2_match,	/* positions of matches in seq2 */
	       void (*pr_func)(void *data, int pw, int y),     /* printing function */
	       void *data             /* printing data */
	       );

int
sip_realloc_matches (int **seq1_match, 
		     int **seq2_match,
		     int **len_match,
		     int *max_matches);

#endif
