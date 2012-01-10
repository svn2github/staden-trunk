#include <stdlib.h>
#include "search_utils.h"

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
char *pstrnstr(char *text, size_t text_len, char *query, size_t query_len) {
    unsigned int t_ind = 0, q_ind;

    /*
     * Simplest implementation.
     * Not an efficient search, but easy to understand.
     */
    do {
	unsigned int t_ind2;
	for (t_ind2 = t_ind, q_ind = 0;
	     q_ind < query_len && t_ind2 < text_len;) {
	    if (text[t_ind2] == '*')
		t_ind2++;
	    else {
		if (text[t_ind2] != query[q_ind])
		    break;
		t_ind2++;
		q_ind++;
	    }
	}
    
	if (q_ind == query_len)
	    return &text[t_ind];

	t_ind++;
    } while (t_ind < text_len);

    return NULL;
}


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
		       int mismatches, int *n_mis) {
    unsigned int t_ind = 0, q_ind;
    int m;

    if (n_mis)
	*n_mis = 0;

    /*
     * Simplest implementation.
     * Not an efficient search, but easy to understand.
     */
    do {
	unsigned int t_ind2;
	m = 0;
	for (t_ind2 = t_ind, q_ind = 0;
	     q_ind < query_len && t_ind2 < text_len;) {
	    if (text[t_ind2] == '*')
		t_ind2++;
	    else {
		if (text[t_ind2] != query[q_ind]) {
		    if (m++ >= mismatches)
			break;
		}
		t_ind2++;
		q_ind++;
	    }
	}
    
	if (q_ind == query_len) {
	    if (n_mis)
		*n_mis = m;
	    return &text[t_ind];
	}

	t_ind++;
    } while (t_ind < text_len);

    return NULL;
}

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
			int mismatches, int *n_mis) {
    unsigned int t_ind = 0, q_ind;
    int m;
    char *match = NULL;
    char match_mis = 0;

    if (n_mis)
	*n_mis = 0;

    /*
     * Simplest implementation.
     * Not an efficient search, but easy to understand.
     */
    do {
	unsigned int t_ind2;
	m = 0;
	for (t_ind2 = t_ind, q_ind = 0;
	     q_ind < query_len && t_ind2 < text_len;) {
	    if (text[t_ind2] == '*')
		t_ind2++;
	    else {
		if (text[t_ind2] != query[q_ind]) {
		    if (m++ >= mismatches)
			break;
		}
		t_ind2++;
		q_ind++;
	    }
	}
    
	if (q_ind == query_len) {
	    match_mis = m;
	    if (n_mis)
		*n_mis = m;
	    match = &text[t_ind];
	}

	t_ind++;
    } while (t_ind < text_len);

    if (n_mis)
	*n_mis = match_mis;

    return match;
}


/*
 * An implementation of strstr that skips pads.
 *
 * As strstr, returns a pointer into the buffer when it finds a match, or
 * NULL when it does not.
 */
char *pstrstr(char *text, char *pattern) {
    char *text_p, *patt_p;

    /*
     * Simplest implementation.
     * Not an efficient search, but easy to understand.
     */
    do {
	text_p = text;
	patt_p = pattern;
	while (*patt_p && *text_p) {
	    if (*text_p == '*')
		text_p++;
	    else {
		if (*text_p != *patt_p)
		    break;
		text_p++;
		patt_p++;
	    }
	}
	if (*patt_p == 0)
	    return text;

    } while (*text && *++text);

    return NULL;
}


/*
 * A version of strstr with inexact string matching.
 *
 * Finds matches between pattern 'p' and text 't' with at most 'm' mismatches.
 * Padding characters are ignored, but this isn't an alignment algorithm - it
 * will not introduce gaps to get a better match.
 * KFB: 22/11/01 added n_mis to return the number of mismatches found
 */
char *pstrstr_inexact(char *text, char *pattern, int mismatches, int *n_mis) {
    char *text_p, *patt_p;
    int m;

    if (n_mis)
	*n_mis = 0;

    /*
     * Simplest implementation.
     * Not an efficient search, but easy to understand.
     */
    do {
	m = 0;

	/* remove leading pads */
	while (*text && *text == '*')
	    text++;

	text_p = text;
	patt_p = pattern;
	while (*patt_p && *text_p) {
	    if (*text_p == '*')
		text_p++;
	    else {
		if (*text_p != *patt_p) {
		    if (m++ == mismatches)
			break;
		}
		text_p++;
		patt_p++;
	    }
	}
	if (*patt_p == 0) {
	    if (n_mis)
		*n_mis = m;
	    return text;
	}

    } while (*text && *++text);

    return NULL;
}

/*
 * A version of strnstr with inexact string matching, but finding the LAST
 * copy of pattern in text.
 * (Done correctly this would do a reverse search, but this is coded for
 * simplicitly at the moment.)
 *
 * Finds matches between pattern 'p' and text 't' with at most 'm' mismatches.
 * Padding characters are ignored, but this isn't an alignment algorithm - it
 * will not introduce gaps to get a better match.
 *
 * This searches from the end of the string rather than the start, so it
 * will return the last copy of pattern in text.
 */
char *prstrstr_inexact(char *text, char *pattern, int mismatches, int *n_mis) {
    char *text_p, *patt_p;
    int m;
    char *match = NULL;
    char match_mis = 0;

    if (n_mis)
	*n_mis = 0;

    /*
     * Simplest implementation.
     * Not an efficient search, but easy to understand.
     */
    do {
	m = 0;

	/* remove leading pads */
	while (*text && *text == '*')
	    text++;

	text_p = text;
	patt_p = pattern;
	while (*patt_p && *text_p) {
	    if (*text_p == '*')
		text_p++;
	    else {
		if (*text_p != *patt_p) {
		    if (m++ == mismatches)
			break;
		}
		text_p++;
		patt_p++;
	    }
	}
	if (*patt_p == 0) {
	    match_mis = m;
	    match = text;
	}

    } while (*text && *++text);

    if (n_mis)
	*n_mis = match_mis;

    return match;
}

