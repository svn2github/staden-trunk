#ifndef FILTER_WORDS_H
#define FILTER_WORDS_H

#include <stddef.h>

/*
 * Search for occurences of rep in seq (of length len).
 * minsize is the minimum length of consecutive 'rep' matches to use for
 * reporting.
 */
int filter_words(char *seq, char *filt, size_t len, char *rep,
		 int minsize, int maxdrop, char filter_char);

/*
 * Search for occurences of rep in seq (of length len).
 * minsize is the minimum length of region containing 'rep' and minscore
 * is effectively the portion of minsize that is 'rep'.
 * So minsize==minscore implies 100% match.
 *    minsize==minscore+2 implies 2 bases mismatch allowed.
 */
int filter_words_local(char *seq, char *filt, size_t len, char *rep,
		       int minsize, int minscore, char filter_char);


#endif /* FILTER_WORDS_H */
