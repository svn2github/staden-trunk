#ifndef _PROBS_H_
#define _PROBS_H_

int FindProbs(char *seq1, char *seq2,int seq1_lreg, int seq1_rreg,
	      int seq2_lreg, int seq2_rreg, int span_length, int seqtype,
	      int use_av_comp);

double FindExpectedProb(char *seq1, char *seq2, int seq1_lreg, int seq1_rreg, 
			int seq2_lreg, int seq2_rreg, int span_length,
			int seqtype);

/*
 * returns the score for the probability at min_prob for span_length
 */
int FindScore(int seq1_len, int seq2_len, int span_length, int num_matches);

/*
 * writes a list of scores and associated probabilities up to a minimum
 * probability to the output window
 */
void ListProbs(char *seq1, char *seq2, int seq1_lreg, int seq1_rreg,
	       int seq2_lreg, int seq2_rreg, int span_length, int seqtype,
	       int min_score, int max_score, int *score_hist);

/*
 * for the find matching words algorithm:
 * writes a list of scores and associated probabilities over a range between
 * min_score and max_score to the output window
 */
void ListIdentityProbs(char *seq1, char *seq2, int seq1_lreg, int seq1_rreg,
		       int seq2_lreg, int seq2_rreg, int seqtype,
		       int min_score, int max_score, int *score_hist);
#endif

