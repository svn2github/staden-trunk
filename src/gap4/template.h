#ifndef _TEMPLATE_H
#define _TEMPLATE_H

#include "IO.h"
#include "list.h"

/* gel_cont_t structure - part of the template_c structure */
typedef struct {
    int read;
    int contig;
} gel_cont_t;

/* Template check structure */
typedef struct {
    /*
     * Template reliability score - computed from consist, flags and the size
     * parameters.
     * 1 = perfect, 0 = utterly wrong.
     */
    double score;

    /*
     * Template direction. Somewhat arbitrary, but it implies that if a
     * reading sequenced using the universal forward primer is present in
     * this template and it has not been complemented then the template itself
     * is not complemented.
     * Note that this can sometimes be incorrect if the template is not
     * consistent (see the "consistency" field and TEMPLATE_DIR macro).
     *
     * 0 == +
     * 1 == -
     * -1 == uninitialised (ie don't know)
     */
    int direction;

    /*
     * These two are initialised when the template_c structure is created.
     */
    list_t *gel_cont;   /* list of gel+contig pairs */
    int num;		/* template number */

    /*
     * See the TEMP_CONSIST_* #defines lower down.
     */
    int consistency;	/* consistency status */

    /*
     * The following fields are specific to a particular contig. When
     * the consistency is checked for the entire template, these fields
     * happen to contain information on the last contig spanned by this
     * template.
     * If TEMP_OFLAG_MINMAX_SIZE is used then we set both start/end and
     * start2/end2. Otherwise we just set start/end to be the average
     * insert size.
     */
    int start;          /* start and end of template, as measured using the */
    int end;		/*   primer information */
    int start2;		/* start and end of template, as measured using the */
    int end2;		/*   primer information */
    int min;		/* start and end of positions of template as */
    int max;		/*   measured using reading positions */
    int flags;		/* flags - see below */

    /*
     * Option flags for controlling the template code. This value defaults
     * to zero for all flags.
     */
    int oflags;		/* option flags - see TEMP_OFLAG_* defines */

    /* Computed length based on templates spanning contigs. */
    int computed_length;
} template_c;


/*
 * Define the consistency values for template_c.consistency.
 * To test, if "t->consistency != 0" then it's invalid. To see why check
 * (eg) "t->consistency & TEMP_CONST_DIST", etc.
 *
 * DIST   == too short or too long.
 * PRIMER == forward/reverse primers not where expected
 * STRAND == inconsistent strand+sense information
 * UNKNOWN== missing primer information (not _necessarily_ inconsistent)
 */
#define TEMP_CONSIST_OK		0
#define TEMP_CONSIST_DIST	(1<<0)
#define TEMP_CONSIST_PRIMER	(1<<1)
#define TEMP_CONSIST_STRAND	(1<<2)
#define TEMP_CONSIST_UNKNOWN	(1<<3)


/*
 * Flags for the template_c.flags field.
 */
#define TEMP_FLAG_NONE		  0
#define TEMP_FLAG_CHECKED_DIST    (1<<0)
#define TEMP_FLAG_CHECKED_PRIMER  (1<<1)
#define TEMP_FLAG_CHECKED_STRAND  (1<<2)
#define TEMP_FLAG_GUESSED_START   (1<<3)
#define TEMP_FLAG_GUESSED_END     (1<<4)
#define TEMP_FLAG_DONE_POSITIONS  (1<<5)
#define TEMP_FLAG_SPANNING	  (1<<6)
#define TEMP_FLAG_EXPECTED	  (1<<7)


/*
 * Flags for the template_c.oflags field.
 */
#define TEMP_OFLAG_NONE		  0
#define TEMP_OFLAG_MINMAX_SIZE	  1 /* min/max insert size instead of avg */
#define TEMP_OFLAG_CVEC		  2 /* also clip cosmid vector */

#define UNKNOWN_POS -999999

/*
 * Length of a template. 't' is a 'template_c *'
 */
#define TEMP_LENGTH(t) (ABS((t)->end - (t)->start))

/*
 * The direction of a template. 0 for +ve, 1 for -ve, -1 for unknown
 */
#define TEMP_DIRECTION(t) (((t)->consistency && TEMP_CONSIST_STRAND) ? -1 \
 : (t)->direction)

/*
 * ---------------------------------------------------------------------------
 * Function Prototypes
 * ---------------------------------------------------------------------------
 */

gel_cont_t *new_gel_cont(void);
void free_gel_cont(gel_cont_t *gc);
template_c **init_template_checks(GapIO *io, int num_contigs,
				  int *contig_array, int connected);
void uninit_template_checks(GapIO *io, template_c **template_check);
int *sort_templates(GapIO *io, template_c **template_check);
void get_template_positions(GapIO *io, template_c *t, int contig);
int last_template_base(GapIO *io, template_c *t, int g);
void check_all_templates(GapIO *io, template_c **template_check);
void remove_single_templates(GapIO *io, template_c **template_check);
void contig_spanning_templates(GapIO *io, template_c **template_check);
void check_template_length(GapIO *io, template_c *t, int overlap);
void check_template_length_overlap(GapIO *io, template_c *t,
				   int c1, int c2, int overlap);

/*
 * Checks whether any of the sequences on template 't' for contig 'c' are 
 * in the range of range_start to range_end inclusive.
 *
 * Returns 0 when not covered.
 *         1 when partially covered.
 *         2 when totally covered.
 */
int template_covered(GapIO *io, template_c *t, int c,
		     int range_start, int range_end);

#endif /* _TEMPLATE_H */
