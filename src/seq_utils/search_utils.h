#ifndef _SEARCH_UTILS_H_
#define _SEARCH_UTILS_H_

/*
 * An implementation of strstr that skips pads ('*'). Pads are only ignored;
 * they will not be introduced to get a better match (as this is not an
 * alignment algorithm). Unlike strstr, our pstrnstr also has a maximum length
 * on both text and query sequences.
 *
 * Arguments:
 *	text		The text to search through
 *	text_len	Length of text
 *	query		The query sequence to find within 'text'
 *	query_len	Length of query
 *
 * Returns:
 *	A pointer to the first match in text when found.
 *	NULL when a match is not found.
 */
char *pstrnstr(char *text, size_t text_len, char *query, size_t query_len);

/*
 * An inexact implementation of pstrnstr.
 * It will find matches allowing at most 'n' mismatches.
 *
 * Arguments:
 *	text		The text to search through
 *	text_len	Length of text
 *	query		The query sequence to find within 'text'
 *	query_len	Length of query
 *	mismatches	Number of allowable mismatches
 *
 * Returns:
 *	A pointer to the first match in text when found.
 *	NULL when a match is not found.
 */
char *pstrnstr_inexact(char *text, size_t text_len,
		       char *query, size_t query_len,
		       int mismatches, int *n_mis);

/*
 * An inexact implementation of prstrnstr.
 * It will find the last match allowing at most 'n' mismatches.
 *
 * Arguments:
 *	text		The text to search through
 *	text_len	Length of text
 *	query		The query sequence to find within 'text'
 *	query_len	Length of query
 *	mismatches	Number of allowable mismatches
 *
 * Returns:
 *	A pointer to the first match in text when found.
 *	NULL when a match is not found.
 */
char *prstrnstr_inexact(char *text, size_t text_len,
			char *query, size_t query_len,
			int mismatches, int *n_mis);

/*
 * An implementation of strstr that skips pads.
 *
 * As strstr, returns a pointer into the buffer when it finds a match, or
 * NULL when it does not.
 */
char *pstrstr(char *text, char *pattern);

/*
 * A version of strstr with inexact string matching.
 *
 * Finds matches between pattern 'p' and text 't' with at most 'm' mismatches.
 * Padding characters are ignored, but this isn't an alignment algorithm - it
 * will not introduce gaps to get a better match.
 */
char *pstrstr_inexact(char *text, char *pattern, int mismatches, int *n_mis);

/* As pstrstr_inexact, but finding the last occurance of pattern */
char *prstrstr_inexact(char *text, char *pattern, int mismatches, int *n_mis);

#endif /* _SEARCH_UTILS_H_ */
