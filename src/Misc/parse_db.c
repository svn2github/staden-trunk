/*
 * File:
 * Version:
 *
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */

#include <staden_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stddef.h>

#include "xalloc.h"
#include "parse_db.h"
#include "misc.h"

/* TOKENS */
#define TOK_NULL	0
#define TOK_END		1
#define TOK_ID		2
#define TOK_SEPARATOR	3
#define TOK_NL		4
#define TOK_EQ		5
    
#define ASCII_TO_HEX(A) (isdigit(A)?(A-'0'):(tolower(A)-'a'+10))

static char *infile;	/* File being parsed */

#define MAXWORD 8192
static char word[MAXWORD];	/* Array containing the TOK_ID word */
static int lineno;	/* Current line number being processed */

/* ---------------- private routines ----------------- */

/*
 * Parser error handler.
 */
static void parse_error(char *s)
{
    fprintf(stderr, "File %s line %d: %s\n", infile, lineno, s);
}

/*
 * Parse out normal escape characters.
 *
 * s = source string
 * t = dest.  string
 *
 * Source and destination may be the same.
 */
static void spring_clean_text (char *s, char *t) {
    while (*s != '\0') {
	if (*s == '\\') {
	    switch (*++s) {
	    case 'a':  *t++ = '\a'; s++; break;
	    case 'b':  *t++ = '\b'; s++; break;
	    case 'f':  *t++ = '\f'; s++; break;
	    case 'n':  *t++ = '\n'; s++; break;
	    case 'r':  *t++ = '\r'; s++; break;
	    case 't':  *t++ = '\t'; s++; break;
	    case 'v':  *t++ = '\v'; s++; break;

	    case '\n':
		while (isspace(*++s));
		break;

	    case 'x': 
		if (isxdigit(s[1]) && isxdigit(s[2])) {
		    *t++ = ASCII_TO_HEX(s[1])*16 + ASCII_TO_HEX(s[2]);
		    s++;s++;s++;
		}
		break;

	    case '\\':
	    case '\?':
	    case '\'':
	    case '\"':
	    case '\0':
		break;
	    }
	} else {
	    *t++ = *s++;
	}
    }
    *t = '\0';
}


/*
 * lexical analyser.
 * Get's the next word from the input stream 'fp' and stores it in the
 * global array 'word'.
 *
 * Returns the token type. (TOK_?)
 */
static int next_word(FILE *fp)
{
    int a;
    int token;
    char *s;
    
    token = TOK_NULL;
    s = word;
    
    while (token == TOK_NULL) {
	switch (a = getc(fp)) {
	case EOF:
	    token = TOK_END;
	    break;

	/* 18/1/99 johnt - handle special case of \r\n on WINNT */
	case '\r':
	    a = getc(fp);
	    if( a != '\n' && a != EOF)
		ungetc(a,fp);
	    /* fall through */
	case '\n':
	    lineno++;
	    token = TOK_NL;
	    break;

	case '=':
	    token = TOK_EQ;
	    break;

	case ':':
	    token = TOK_SEPARATOR;
	    break;

	case '#':
	    /* comment: skip to end of line */
	    for(a=getc(fp); a!=EOF && a!='\n'; a=getc(fp));
	    if (a=='\n') lineno++;
	    if (a==EOF) token = TOK_END;
	    break;

	case '\\':
	    /* back quoted newlines are skipped */
	    /* back quoted "anything else" is "anything else" */
	    a = getc(fp);
	    if( a == '\r' )
		a = getc(fp); /* 18/1/99 johnt - handle \r\n on WINNT */
	    if (a != EOF && a != '\n' && a != '\r') ungetc(a,fp);
	    break;

	case '"':
	    /* quoted string */
	    for(a=getc(fp);a!=EOF && a!='"';a=getc(fp)) {
		if (a=='\n') lineno++;
		if (a=='\\') {
		    *s++ = a;
		    if (EOF == (a=getc(fp)))
			break;
		}
		if( a != '\r' ) /* 18/1/99 johnt - handle \r\n on WINNT */
		  *s++ = a;
	    }
	    token = TOK_ID;
	    break;

	default: /* An identifier */
	    if (isalnum(a)) {
		*s++ = a;
		for(a = getc(fp); a != EOF && isalnum(a); a = getc(fp))
		    *s++ = a;
		if (a != EOF)
		    ungetc(a, fp);

		token = TOK_ID;
	    }
	    break;
	}
    }
    
    *s = '\0';
    spring_clean_text(word,word);
    return token;
}

static int parse_entry(FILE *fp, pf_spec *spec, char *store) {
    int eoe = 0;
    int eof = 0;
    int i;
    char *valp = NULL;
    int type = 0;
    int side = 0; /* 0 == left, 1 == right */
    int done_rhs = 0;

    while (!eoe) {
	switch(next_word(fp)) {
	case TOK_END:
	    eof = 1;
	case TOK_NL:
	    eoe = 1;
	    break;

	case TOK_EQ:
	    side = 1;
	    break;

	case TOK_SEPARATOR:
	    side = 0;
	    valp = NULL;
	    done_rhs = 0;
	    break;

	case TOK_ID:
	    /* printf("\tident%d=\"%s\"\n", side, word); */

	    if (side == 0) {
		for (i = 0; spec[i].name; i++) {
		    if (word[0] == spec[i].name[0] &&
			strcmp(word, spec[i].name) == 0) {
			valp = &store[spec[i].offset];
			type = spec[i].type;

			break;
		    }
		}
		if (!spec[i].name) {
		    char buf[MAXWORD];
		    sprintf(buf, "Warning - unknown identifier \"%s\"\n",word);
		    parse_error(buf);
		}

	    } else { /* side == right */

		if (done_rhs) {
		    parse_error("Syntax error");
		    break;
		}

		done_rhs = 1;

		if (valp) {
		    switch(type) {
		    case PF_STR:
			*((char **)valp) = strdup(word);
			break;
		    case PF_INT:
			*((int *)valp) = atoi(word);
			break;
		    default:
			{
			    char buf[MAXWORD];
			    sprintf(buf, "Unknown type %d\n", type);
			    parse_error(buf);
			}
			break;
		    }
		}
	    }
	    break;

	}
    }

    return eof;
}

static void *find_store(void *store, int store_size, int nitems, char *word) {
    int i;
    void *storep = store;

    for (i = 0; i < nitems; i++) {
	if (**((char **)storep) == word[0] &&
	    0 == strcmp(*((char **)storep), word)) {
	    return storep;
	}
	storep = (char *)storep + store_size;
    }

    return NULL;
}


/* ---------------- external ----------------- */

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
		 int store_size, void *default_store) {
    int eof = 0;
    int count = *nitems;
    char *storep;
    FILE *fp;

    infile = fn;
    lineno = 0;
    fp = fopen(fn, "rb");

    if (!fp) {
	parse_error("Could not open");
	return NULL;
    }

    while (!eof) {
	switch (next_word(fp)) {
	case TOK_END:
	    eof = 1;
	    break;

	case TOK_NL:
	    break;

	case TOK_ID:
	    if (!(storep = find_store(store, store_size, count, word))) {
		count++;
		store = xrealloc(store, count * store_size);
		storep = (char *)store + (count-1) * store_size;
		if (default_store)
		    memcpy(storep, default_store, store_size);
		else
		    memset(storep, 0, store_size);
		*((char **)storep) = strdup(word);
	    }
	    eof = parse_entry(fp, spec, storep);
	    break;

	default:
	    parse_error("Syntax error - stopped parsing");
	    fclose(fp);
	    return NULL;
	}
    }

    fclose(fp);
    
    *nitems = count;
    return store;
}
