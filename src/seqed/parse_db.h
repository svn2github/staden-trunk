#ifndef _PARSE_DB_H_
#define _PARSE_DB_H_

#include <stddef.h> /* offsetof() */

#define PF_INT    1
#define PF_STR    2
#define PF_IO     3
#define PF_ARR    4
#define PF_FLOAT  5
#define PF_DOUBLE 6

typedef struct {
    char *name;
    int type;
    int offset;
} pf_spec;

/* 22/1/99 johnt - use ; path seperators for WINNT */
#ifdef _WIN32
#  define PATHSEP ';'
#else
#  define PATHSEP ':'
#endif

/*
 * Parses a file (in the format used by TAGDB and NOTEDB).
 *
 * fn		File to parse
 *
 * spec		The specification of the format for the database. This
 *		consists of a series of identifier,type,offset tuples, ending
 *		with a NULL identifier.
 *
 * store	The database store, with offsets specified in spec.
 *		Specify as NULL for the first file, and the current database
 *		pointer for subsequent files.
 *
 * nitems	The number of items in the database so far.
 *
 * store_size	The size of each structure in store.
 *
 * default	A default structure for the store, or NULL if none.
 *
 * Returns a new copy of the database ('store'), or NULL for failure.
 */
void *parse_file(char *fn, pf_spec *spec, void *store, int *nitems,
		 int store_size, void *default_store);

#endif /* _PARSE_DB_H_ */
