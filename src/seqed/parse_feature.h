#ifndef _PARSE_FEAT_H
#define _PARSE_FEAT_H

#include <tcl.h>

/*
 * A single location type:
 * exact:			"10"
 * single base in a range:	"10.20"
 * site within a range:		"10^20"
 */
typedef enum {
    EXACT,
    BASE,
    SIT
} ft_location_type;

/*
 * A location (half of a range).
 * Locations are either exact or ranges. They may also have qualifiers to
 * indicate whether they are before or after a certain point (eg "<10").
 */
typedef struct ft_location_t {
    int min, max; 		/* minimum and maximum values in "10.20" or
				   "10^20" type */
    int min_lt;			/* -1 = <, 1 = >, 0 = none */
    int max_lt;			/* -1 = <, 1 = >, 0 = none */
    ft_location_type type;	/* Distinguishes "10", "10.20" or "10^20" */
} ft_location;

/*
 * A range: consists of two locations separated by "..".
 * Eg "10..20" or "(<10.20)..(30.40)".
 */
typedef struct ft_range_t {
    /*
     * TODO: Add order(), join() and complement() information. This
     * may require changing this to a hierarchial structure.
     */
    ft_location *left;		/* Two locations for loc..loc. eg
				   "10..(30.40)" => 10 && (30.40) */
    ft_location *right;		/* Right is 2nd location if present or 
				   NULL if not. */
    int complemented;		/* 0 == no, 1 == yes */
    struct ft_range_t *next;	/* Linked list; to next range */
} ft_range;

/*
 * A single feature entry: consists of a feature type, a linked list of
 * ranges (feature locations), and zero or more qualifiers.
 */
typedef struct ft_t {
    char type[20];		/* feature type in string form */
    char *location;		/* textual version of the location string */
    ft_range *range;		/* A linked list of ranges */
    char *qualifiers;		/* Qualifiers as a single block of
				   text. TODO: split this into a
				   suitable data structure. */
    int qual_hash_init;		/* 1 if qual_hash is in use */
    Tcl_HashTable qual_hash;	/* Qualifiers as a Tcl hash table */
} ft_entry;


/*
 * A value item in the qualifier hash table. These are a linked list as the
 * qualifiers may be duplicated.
 */
struct ft_value_element_t;
typedef struct ft_value_element_t {
    char *value;
    struct ft_value_element_t *next;
} ft_value_element;

char *get_ft_entry ( FILE *fp, char *line, char **last_line );
ft_entry *parse_ft_entry(char *str);
void print_entry(ft_entry *e);
void print_range(ft_range *r);
ft_entry *new_ft_entry(void);
ft_range *new_ft_range(void);
ft_location *new_ft_location(void);
void init_ft_qual_hash(ft_entry *e, char *str);
void del_ft_entry(ft_entry *e);
void del_ft_range(ft_range *r);
void del_ft_location(ft_location *l);
ft_value_element *search_ft_qual_hash(ft_entry *e, char *qual);
int parse_ft_location(ft_entry *entry, char *s);



#endif 
