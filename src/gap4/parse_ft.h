#ifndef _PARSE_FT_H
#define _PARSE_FT_H

#include <tcl.h>

/*
 * Parses the EMBL FT lines.
 * A feature table consists of a type, a location line, and qualifiers.
 * For example:
 *
 * FT   CDS             join(complement(115903..116071),complement(95491..95621),
 * FT                   complement(93634..93714),complement(92556..92682),
 * FT                   complement(87329..87423),complement(82749..82887),
 * FT                   complement(79690..79802))
 * FT                   /db_xref="SPTREMBL:O43728"
 * FT                   /evidence=NOT_EXPERIMENTAL
 * FT                   /note="match: proteins: Tr:O43728 Tr:O88872 Sw:P52840
 * FT                   Tr:O75897 Sw:P50224 Sw:P50225 Sw:P17988 Sw:P50226 Sw:P52846
 * FT                   Tr:O43704 Tr:O46503 Tr:O00338 Tr:O70262 Tr:O35403 Tr:O46640
 * FT                   Sw:P19217 Sw:P50237"
 * FT                   /gene="dJ388M5.3"
 * FT                   /product="dJ388M5.3 (novel Sulfotransferase (sulfokinase,
 * FT                   EC 2.8.2.1) like protein)"
 * FT                   /protein_id="CAB09788.1"
 * FT                   /translation="MAESEAETPSTPGEFESKYFEFHGVRLPPFCRGKMEEIANFPVRP
 * FT                   SDVWIVTYPKSGTSLLQEVVYLVSQGADPDEIGLMNIDEQLPVLEYPQPGLDIIKELTS
 * FT                   PRLIKSHLPYRFLPSDLHNGDSKVIYMARNPKDLVVSYYQFHRSLRTMSYRGTFQEFCR
 * FT                   RFMNDKLGYGSWFEHVQEFWEHRMDSNVLFLKYEDMHRDLVTMVEQLARFLGVSCDKAQ
 * FT                   LEALTEHCHQLVDQCCNAEALPVGRGRVGLWKDIFTVSMNEKFDLVYKQKMGKCDLTFD
 * FT                   FYL"
 *
 * Ranges consist of one or more range joined together.
 * For example join(55798..55943,56245..57257). Each range has a left and
 * a right location. Each location is typically expressed exactly, but may
 * itself be a range if the exact location is unknown. Eg: <10, 10^20, 10.20.
 * So a single range may consist of "<100..(150.160)".
 *
 * Qualifiers are stored as /qual, /qual=value or /qual="value".
 * /qual may occur more than once in a single FT entry, as it typical with
 * /dbxref qualifiers.
 *
 * The main data structure for holding this is ft_entry.
 * Within this the type holds "CDS". The location and qualifier fields hold
 * textual information while range and qual_hash contain processed data
 * structures.
 */

/* ---------------------------------------------------------------------- */
/* Data types */


/*
 * A single location type:
 * exact:			"10"
 * single base in a range:	"10.20"
 * site within a range:		"10^20"
 */
typedef enum {
    EXACT,
    BASE,
    SITE
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

/* ---------------------------------------------------------------------- */
/* Function prototypes */

ft_location *new_ft_location(void);
void del_ft_location(ft_location *l);
ft_range *new_ft_range(void);
void del_ft_range(ft_range *r);
ft_entry *new_ft_entry(void);
void del_ft_entry(ft_entry *e);
ft_entry *parse_ft_entry(char *str);
char *get_ft_qualifier(char *str, char *qualifier, int *value_len);

void init_ft_qual_hash(ft_entry *e, char *str);
void del_ft_qual_hash(ft_entry *e);
ft_value_element *search_ft_qual_hash(ft_entry *e, char *qual);

void print_entry(ft_entry *e);

#endif /* _PARSE_FT_H */
