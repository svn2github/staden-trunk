#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "xalloc.h"
#include "parse_ft.h"

/*
 * Parses the feature table (held in the FT lines of an EMBL or Genbank file).
 */

ft_location *new_ft_location(void) {
    ft_location *l;

    l = (ft_location *)xmalloc(sizeof(*l));
    if (!l)
	return NULL;

    l->min = 0;
    l->max = 0;
    l->min_lt = 0;
    l->max_lt = 0;
    l->type = EXACT;

    return l;
}

void del_ft_location(ft_location *l) {
    if (l)
	xfree(l);
}

ft_range *new_ft_range(void) {
    ft_range *r;

    r = (ft_range *)xmalloc(sizeof(*r));
    if (!r)
	return NULL;

    r->left = NULL;
    r->right = NULL;
    r->complemented = 0;
    r->next = NULL;

    return r;
}

void del_ft_range(ft_range *r) {
    if (!r)
	return;

    /* Recursively delete the entire list */
    del_ft_range(r->next);

    if (r->left)
	del_ft_location(r->left);
    if (r->right)
	del_ft_location(r->right);
    xfree(r);
}

ft_entry *new_ft_entry(void) {
    ft_entry *e;
    e = (ft_entry *)xmalloc(sizeof(*e));
    if (!e)
	return NULL;

    e->range = NULL;
    e->qualifiers = NULL;
    e->qual_hash_init = 0;

    return e;
}

void del_ft_entry(ft_entry *e) {
    if (!e)
	return;

    if (e->range)
	del_ft_range(e->range);

    if (e->qualifiers)
	xfree(e->qualifiers);

    del_ft_qual_hash(e);
    
    xfree(e);
}

static void print_loc(ft_location *l) {
    if (l->type != EXACT)
	printf("{");

    if (l->min_lt)
	putchar("< >"[l->min_lt+1]);

    printf("%d", l->min);

    if (l->type == EXACT)
	return;

    printf("%c", (l->type == BASE) ? '.' : '^');

    if (l->max_lt)
	putchar("< >"[l->max_lt+1]);

    printf("%d", l->max);

    printf("}");
}

static void print_range(ft_range *r) {
    printf("RANGE='");
    if (r->complemented)
	printf("complement(");
    if (r->left)
	print_loc(r->left);
    if (r->right) {
	printf("..");
	print_loc(r->right);
    }
    if (r->complemented)
	printf(")");
    printf("'\n");
}

void print_entry(ft_entry *e) {
    ft_range *r;

    puts("\n>>>>>");
    if (!e) {
	printf("NULL entry");
	puts("<<<<<\n");
	return;
    }

    printf("Type='%s'\n", e->type);
    printf("Location='%s'\n", e->location);

    for (r = e->range; r; r = r->next) {
	print_range(r);
    }

    printf("Qualifiers='%s'\n", e->qualifiers);

    puts("<<<<<\n");
}


/*
 * Parses just the location portion of a feature table.
 * Returns 0 for success,
 *        -1 for failure.
 */
static int parse_ft_location(ft_entry *entry, char *s) {
    int start, end;
    ft_location *loc;
    ft_range *range = NULL;
    int complemented = 0;

    enum state_t {
	STATE_TEXT,         /* 0 default state when between ranges */
	STATE_START,        /* 1 The 10 in "10.20" */
	STATE_SINGLE_DOT,   /* 2 The dot in "10.20" */
	STATE_END,          /* 3 The 20 in "10.20" */
	STATE_DOUBLE_DOT,   /* 4 The .. in "10..20"; separates start and end */
	STATE_CARET,	    /* 5 The ^ in "10^20" */
	STATE_EXIT    	    /* 6 Newline or nul char */
    } state;

    start = end = 0;
    state = STATE_TEXT;
    loc = NULL;

    while (state != STATE_EXIT) {
	/* printf("STATE=%d, *s = '%c'\n", state, *s); */

	switch (state) {
	case STATE_TEXT:
	    if (*s == '<' || *s == '>' || isdigit(*s)) {
		ft_range *tmp;

		if (NULL == (tmp = new_ft_range()))
		    goto error;

		if (!entry->range) {
		    entry->range = range = tmp;
		} else {
		    range->next = tmp;
		    range = tmp;
		}

		if (NULL == (loc = new_ft_location()))
		    goto error;
		    
		range->left = loc;
		range->complemented = complemented;
		complemented = 0;
		state = STATE_START;
		continue; /* consume no characters */
	    } else if (*s == 'c') {
		/* only valid 'c' in the location string is in "complement" */
		complemented = 1;
	    }
	    break;

	case STATE_START:
	    if (*s == '<')
		loc->min_lt = -1;
	    else if (*s == '>')
		loc->min_lt = 1;
	    else if (isdigit(*s)) {
		loc->max = loc->min = loc->min * 10 + *s - '0';
	    } else if (*s == ')') {
		/* skip brackets */
		break;
	    } else if (*s == '.') {
		state = STATE_SINGLE_DOT;
	    } else if (*s == '^') {
		state = STATE_CARET;
	    } else {
		state = STATE_TEXT;
	    }
	    break;

	case STATE_SINGLE_DOT:
	    if (*s == '.') {
		state = STATE_DOUBLE_DOT;
	    } else if (*s == '<' || *s == '>' || isdigit(*s)) {
		loc->type = BASE;
		loc->max = 0;
		state = STATE_END;
		continue; /* consume no characters */
	    } else if (*s == '(') {
		/* skip brackets */
		break;
	    } else {
		goto error;
	    }
	    break;

	case STATE_END:
	    if (*s == '<')
		loc->max_lt = -1;
	    else if (*s == '>')
		loc->max_lt = 1;
	    else if (isdigit(*s)) {
		loc->max = loc->max * 10 + *s - '0';
	    } else if (*s == ')') {
		/* skip brackets */
		break;
	    } else if (*s == '.') {
		state = STATE_SINGLE_DOT;
	    } else {
		state = STATE_TEXT;
	    }
	    break;

	case STATE_CARET:
	    loc->type = SITE;
	    if (*s == '<' || *s == '>' || isdigit(*s)) {
		loc->max = 0;
		state = STATE_END;
		continue; /* consume no characters */
	    } else {
		goto error;
	    }
	    break;

	case STATE_DOUBLE_DOT:
	    if (*s == '<' || *s == '>' || isdigit(*s)) {
		if (NULL == (loc = new_ft_location()))
		    goto error;

		range->right = loc;
		state = STATE_START;
		continue; /* consume no characters */
	    } else if (*s == '(') {
		/* skip brackets */
		break;
	    } else {
		state = STATE_TEXT;
	    }
	    break;

	case STATE_EXIT:
	    /* Do nothing, but it removes the warning from gcc -Wswitch */
	    break;
	}

	if (!*s || *s == '/') {
	    state = STATE_EXIT;
	}

	++s;
    }

    return 0;

 error:
    return -1;
}


#if 0
int range_sort_func(const void *v1, const void *v2) {
    ft_range **r1 = (ft_range **)v1;
    ft_range **r2 = (ft_range **)v2;

    /* Sort by left min position */
    if ((*r1)->left && (*r2)->left)
	return (*r1)->left->min - (*r2)->left->min;
    else
	return 0;
}

/*
 * Sorts the location ranges in an FT entry into left to right order if
 * not complemented, and right to left if complemented.
 *
 * Returns 0 for success
 *        -1 for failure
 */
static int sort_ft_location(ft_entry *e) {
    int i, n;
    ft_range **ranges, *rptr, **rptrptr;
    
    /* Count number of ranges */
    for (n = 0, rptr = e->range; rptr; rptr = rptr->next)
	n++;

    /* Create an array of range pointers */
    if (NULL == (ranges = (ft_range **)xmalloc(n * sizeof(ft_range *))))
	return -1;

    for (n = 0, rptr = e->range; rptr; rptr = rptr->next) {
	ranges[n++] = rptr;
    }

    /* Sort the array */
    qsort(ranges, n, sizeof(ft_range *), range_sort_func);

    /* Relink the ranges list */
    if (e->strand == 0) {
	rptrptr = &e->range;
	for (i = 0; i < n; i++) {
	    *rptrptr = ranges[i];
	    rptrptr = &(ranges[i]->next);
	}
	*rptrptr = NULL;
    } else {
	rptrptr = &e->range;
	for (n--; n >= 0; n--) {
	    *rptrptr = ranges[n];
	    rptrptr = &(ranges[n]->next);
	}
	*rptrptr = NULL;
    }

    xfree(ranges);

    return 0;
}
#endif

ft_entry *parse_ft_entry(char *str) {
    ft_entry *entry = NULL;
    int i, j;
    int loc_start;

    if (NULL == (entry = new_ft_entry()))
	goto error;

    /* Parse 1st line to get feature type */
    for (i = 0; str[i] && !isspace(str[i]) && i < 19; i++) {
	entry->type[i] = str[i];
    }
    entry->type[i] = 0;

    /* Extract text copy of location */
    while (isspace(str[i]))
	i++;
    loc_start = i;
    while (str[i] && str[i] != '/')
	i++;

    entry->location = (char *)xmalloc(i - loc_start + 1);
    if (!entry->location)
	goto error;
    strncpy(entry->location, &str[loc_start], i - loc_start);
    entry->location[i-loc_start] = 0;

    /* Process qualifiers */
    if (str[i]) {
	size_t len = strlen(str);

	entry->qualifiers = (char *)xmalloc(len-i+1);
	if (!entry->qualifiers)
	    goto error;

	while (i < len && isspace(str[i]))
	    i++;
	for (j = 0; i < len && str[i];) {
	    entry->qualifiers[j++] = str[i++];
	    if (str[i-1] == '\n') {
		while (i < len && isspace(str[i]))
		    i++;
	    }
	}
	entry->qualifiers[j] = 0;
    }

    /* Hash the qualifiers too */
    init_ft_qual_hash(entry, entry->qualifiers);

    /* Parse location */
    if (-1 == parse_ft_location(entry, str))
	goto error;

#if 0
    /* Sort locations in left to right order (right to left if complemented) */
    sort_ft_location(entry);
#endif

    return entry;

 error:
    if (entry)
	del_ft_entry(entry);

    return NULL;
}

ft_value_element *new_ft_value_element(void) {

    ft_value_element* e;
    e = (ft_value_element*) xmalloc( sizeof(*e) );
    if( e ) {
        e->value = NULL;
        e->next  = NULL;
    }
    return e;
}

void del_ft_value_element_list(ft_value_element *e) {
    if (!e)
	return;

    del_ft_value_element_list(e->next);
    if (e->value)
	xfree(e->value);
    xfree(e);
}

/*
 * Frees memory allocated to store the qualifier hash table.
 */
void del_ft_qual_hash(ft_entry *e) {
    Tcl_HashEntry *hash;
    Tcl_HashSearch s;
	
    if (!e->qual_hash_init)
	return;

    /* Delete the malloced value strings contained in the table */
    for (hash = Tcl_FirstHashEntry(&e->qual_hash, &s);
	 hash;
	 hash = Tcl_NextHashEntry(&s)) {
	del_ft_value_element_list((ft_value_element *)Tcl_GetHashValue(hash));
    }

    /* Delete the table itself */
    Tcl_DeleteHashTable(&e->qual_hash);
    e->qual_hash_init = 0;
}

/*
 * This parses the qualifier information held in 'str' and stores it in a Tcl
 * hash table, which may be accessed via search_ft_qual_hash.
 *
 * Memory allocated in this function is owned by these routines. It needs to be
 * freed using del_ft_qual_hash.
 *
 * Qualifiers value paris are of the formats:
 * /qualifier
 * /qualifier=value
 * /qualifier="value"
 *
 * Quotes may be escaped using "". Eg:
 * /qualifier="a quote ("") character".
 *
 * Quoted values may span multiple lines, but non quoted ones do not.
 *
 * http://www.ebi.ac.uk/embl/Documentation/FT_definitions/feature_table.html
 */
void init_ft_qual_hash(ft_entry *e, char *str) {
    char *cp, *cp_alloc;
    char *qualifier = NULL;
    char *value = NULL;
    int not_done = 1;
    enum state_t {
	STATE_START,		/* start */
	STATE_QUALIFIER,	/* inside a qualifier */
	STATE_EQUALS,		/* found an = */
	STATE_VALUE,		/* inside a value */
	STATE_QVALUE,		/* inside a quoted value */
	STATE_QUOTE,		/* " in value, "" escapes " */
	STATE_STORE		/* at end of 'value'; add to hash table */
    } state;

    if (!str)
	return;

    /* Create hash table */
    Tcl_InitHashTable(&e->qual_hash, TCL_STRING_KEYS);
    e->qual_hash_init = 1;

    cp = cp_alloc = strdup(str); /* We want to write to our copy */
    state = STATE_START;
    while (not_done) {
	/*
	fprintf(stderr, "state=%d *cp=%c not_done=%d\n", state, *cp, not_done);
	*/
	switch(state) {
	case STATE_START:
	    qualifier = NULL;
	    value = NULL;
	    if (*cp == '/') {
		state = STATE_QUALIFIER;
		cp++;
		qualifier = cp;
	    } else if (*cp) {
		cp++;
	    } else {
		not_done = 0;
	    }
	    break;

	case STATE_QUALIFIER:
	    if (*cp == '=') {
		*cp = 0;
		state = STATE_EQUALS;
		cp++;
	    } else if (*cp == '\n' || *cp == '\0') {
		state = STATE_STORE;
	    } else {
		cp++;
	    }
	    break;

	case STATE_EQUALS:
	    if (*cp == '"') {
		state = STATE_QVALUE;
		cp++;
		value = cp;
	    } else if (*cp == '\n' || *cp == '\0') {
		state = STATE_STORE;
	    } else {
		value = cp;
		state = STATE_VALUE;
	    }
	    break;

	case STATE_QVALUE:
	    /* Value extends up to trailing quote' */
	    if (*cp == '\0') {
		state = STATE_STORE;
	    } else if (*cp == '"') {
		state = STATE_QUOTE;
		cp++;
	    } else {
		cp++;
	    }
	    break;

	case STATE_QUOTE:
	    if (*cp == '"') {
		state = STATE_QVALUE;
		/* Not particularly efficient, but it works */
		memmove(cp, cp+1, strlen(cp+1));
	    } else {
		cp--; /* back up one byte to real end of string */
		state = STATE_STORE;
	    }
	    break;

	case STATE_VALUE:
	    /* Value extends up to end of line, quotes are literal */
	    if (*cp == '\n' || *cp == '\0') {
		state = STATE_STORE;
	    } else {
		cp++;
	    }
	    break;

	case STATE_STORE:
	    /* Store the (qualifier,value) pair in the hash table */
	    {
		int new;
		Tcl_HashEntry *hash;
		ft_value_element *ele = new_ft_value_element();
		char tmp;

		tmp = *cp;
		if (*cp) {
		    *cp = 0;
		} else {
		    not_done = 0;
		}

		if (value && *value)
		    ele->value = strdup(value);
		else
		    ele->value = NULL;

		hash = Tcl_CreateHashEntry(&e->qual_hash, qualifier, &new);
		*cp = tmp;

		if (!new) {
		    /* Note, this reverses the order of the linked list */
		    ft_value_element* curr = (ft_value_element*) Tcl_GetHashValue(hash);
		    curr->next = ele;
		    ele = curr;
		}

		Tcl_SetHashValue(hash, (ClientData)ele);
		state = STATE_START;
	    }
	    break;
	}
    }

    xfree(cp_alloc);
}

/*
 * Searches for a qualifier in a hashed qualifier table.
 * Returns the pointer to the value, or NULL if none found. The memory returned
 * is owned by this function and should not be freed by the caller.
 *
 * Values are returned as a linked list, so that if multiple /qualifier copies
 * exist then the values can be found.
 */
ft_value_element *search_ft_qual_hash(ft_entry *e, char *qual) {
    Tcl_HashEntry *hash;

    if (!e->qual_hash_init)
	return NULL;

    if (!(hash = Tcl_FindHashEntry(&e->qual_hash, qual)))
	return NULL;

    return (ft_value_element *)Tcl_GetHashValue(hash);
}
